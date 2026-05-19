#!/usr/bin/env python3
"""
Export LowMachReact-Hex VTU/PVTU output to CSV.

This utility is intended for post-processing parallel PVTU output without
opening ParaView. It can export:

  1. all cell-centered data from one or more snapshots
  2. x-centerline data through a y/z location
  3. compact per-snapshot summary statistics

Examples
--------

Single centerline:

  python tools/export_pvtu_csv.py line \
    --vtu cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu \
    --out cases/counterflow_nonreacting_gmsh/validation/centerline_005000.csv

All matching snapshots, one centerline CSV per snapshot:

  python tools/export_pvtu_csv.py line \
    --vtu-glob 'cases/counterflow_nonreacting_gmsh/output/flow_*.pvtu' \
    --outdir cases/counterflow_nonreacting_gmsh/validation/centerlines

All cells for one snapshot, selected fields only:

  python tools/export_pvtu_csv.py cells \
    --vtu cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu \
    --fields velocity,temperature,Y_CH4,Y_O2,Y_N2,rho,rho_thermo,enthalpy,species_enthalpy_diffusion \
    --out cases/counterflow_nonreacting_gmsh/validation/cells_005000.csv

Summary over all output snapshots:

  python tools/export_pvtu_csv.py summary \
    --vtu-glob 'cases/counterflow_nonreacting_gmsh/output/flow_*.pvtu' \
    --fields temperature,Y_CH4,Y_O2,Y_N2,rho_thermo,enthalpy,species_enthalpy_diffusion \
    --out cases/counterflow_nonreacting_gmsh/validation/output_summary.csv

Requirements
------------
  pip install meshio numpy pandas
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable

import meshio
import numpy as np
import pandas as pd


def parse_step(path: Path) -> int:
    m = re.search(r"flow_(\d+)", path.name)
    return int(m.group(1)) if m else -1


def resolve_vtus(args: argparse.Namespace) -> list[Path]:
    paths: list[Path] = []
    if getattr(args, "vtu", None):
        paths.append(Path(args.vtu))
    if getattr(args, "vtu_glob", None):
        paths.extend(Path().glob(args.vtu_glob))
    paths = sorted(set(paths), key=lambda p: (parse_step(p), str(p)))
    if not paths:
        raise SystemExit("No VTU/PVTU files found.")
    return paths


def resolve_piece_paths(vtu_path: str | Path) -> list[Path]:
    path = Path(vtu_path)
    if path.suffix.lower() == ".pvtu":
        stem = path.stem
        pieces = sorted(path.parent.glob(f"{stem}_P*.vtu"))
        if not pieces:
            raise RuntimeError(
                f"No VTU pieces matched {path.parent}/{stem}_P*.vtu for {path}."
            )
        return pieces
    return [path]


def hex_block_index(mesh: meshio.Mesh) -> int:
    for i, block in enumerate(mesh.cells):
        if block.type in ("hexahedron", "hex8"):
            return i
    return 0


def cell_centers_and_data(mesh: meshio.Mesh) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    bi = hex_block_index(mesh)
    cells = mesh.cells[bi].data
    points = mesh.points
    centers = points[cells].mean(axis=1)

    cell_data: dict[str, np.ndarray] = {}
    for name, blocks in mesh.cell_data.items():
        if len(blocks) > bi:
            cell_data[name] = np.asarray(blocks[bi])
    return centers, cell_data


def read_snapshot(vtu_path: str | Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    centers_parts: list[np.ndarray] = []
    data_parts: dict[str, list[np.ndarray]] = {}

    for p in resolve_piece_paths(vtu_path):
        mesh = meshio.read(str(p))
        centers, cell_data = cell_centers_and_data(mesh)
        centers_parts.append(centers)
        for name, arr in cell_data.items():
            data_parts.setdefault(name, []).append(np.asarray(arr))

    centers_cat = np.vstack(centers_parts)
    merged: dict[str, np.ndarray] = {}
    for name, parts in data_parts.items():
        try:
            merged[name] = np.concatenate(parts, axis=0)
        except Exception as exc:
            shapes = [p.shape for p in parts]
            raise RuntimeError(f"Could not concatenate field {name}; shapes={shapes}") from exc
    return centers_cat, merged


def field_names(data: dict[str, np.ndarray], fields_arg: str | None) -> list[str]:
    if not fields_arg or fields_arg.lower() == "all":
        return sorted(data.keys())
    requested = [s.strip() for s in fields_arg.split(",") if s.strip()]
    missing = [f for f in requested if f not in data]
    if missing:
        print(f"WARNING: missing requested fields: {', '.join(missing)}")
    return [f for f in requested if f in data]


def add_field_columns(df: pd.DataFrame, data: dict[str, np.ndarray], fields: Iterable[str], mask: np.ndarray | None = None) -> None:
    for name in fields:
        arr = np.asarray(data[name])
        if mask is not None:
            arr = arr[mask]

        if arr.ndim == 1:
            df[name] = arr
        elif arr.ndim == 2 and arr.shape[1] == 3:
            df[f"{name}_x"] = arr[:, 0]
            df[f"{name}_y"] = arr[:, 1]
            df[f"{name}_z"] = arr[:, 2]
        elif arr.ndim == 2:
            for j in range(arr.shape[1]):
                df[f"{name}_{j}"] = arr[:, j]
        else:
            flat = arr.reshape(arr.shape[0], -1)
            for j in range(flat.shape[1]):
                df[f"{name}_{j}"] = flat[:, j]


def estimate_centerline_tolerance(centers: np.ndarray) -> float:
    spacings: list[float] = []
    for axis in range(3):
        vals = np.unique(np.round(centers[:, axis], 12))
        vals = np.sort(vals)
        if vals.size < 2:
            continue
        diffs = np.diff(vals)
        span = float(vals[-1] - vals[0])
        diffs = diffs[diffs > max(1.0e-10, 1.0e-6 * span)]
        if diffs.size:
            spacings.append(float(np.median(diffs)))
    return 0.60 * max(spacings) if spacings else 0.05


def out_path_for_snapshot(base_outdir: Path, vtu: Path, suffix: str) -> Path:
    step = parse_step(vtu)
    if step >= 0:
        return base_outdir / f"{suffix}_{step:06d}.csv"
    return base_outdir / f"{vtu.stem}_{suffix}.csv"


def export_cells(args: argparse.Namespace) -> None:
    vtus = resolve_vtus(args)

    if len(vtus) > 1 and not args.outdir:
        raise SystemExit("--outdir is required when exporting multiple snapshots.")

    outdir = Path(args.outdir) if args.outdir else None
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)

    for vtu in vtus:
        centers, data = read_snapshot(vtu)
        fields = field_names(data, args.fields)

        df = pd.DataFrame({
            "step": parse_step(vtu),
            "x": centers[:, 0],
            "y": centers[:, 1],
            "z": centers[:, 2],
        })
        add_field_columns(df, data, fields)

        out = Path(args.out) if args.out and len(vtus) == 1 else out_path_for_snapshot(outdir, vtu, "cells")
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        print(f"Wrote: {out}  rows={len(df)} cols={len(df.columns)}")


def export_line(args: argparse.Namespace) -> None:
    vtus = resolve_vtus(args)

    if len(vtus) > 1 and not args.outdir:
        raise SystemExit("--outdir is required when exporting multiple snapshots.")

    outdir = Path(args.outdir) if args.outdir else None
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)

    for vtu in vtus:
        centers, data = read_snapshot(vtu)
        fields = field_names(data, args.fields)

        mins = centers.min(axis=0)
        maxs = centers.max(axis=0)
        mids = 0.5 * (mins + maxs)

        ymid = float(args.ymid) if args.ymid is not None else float(mids[1])
        zmid = float(args.zmid) if args.zmid is not None else float(mids[2])
        tol = float(args.tol) if args.tol is not None else estimate_centerline_tolerance(centers)

        mask = (np.abs(centers[:, 1] - ymid) <= tol) & (np.abs(centers[:, 2] - zmid) <= tol)
        if not np.any(mask):
            raise RuntimeError(f"No cells within tol={tol} of y={ymid}, z={zmid} for {vtu}")

        c = centers[mask]
        df = pd.DataFrame({
            "step": parse_step(vtu),
            "x": c[:, 0],
            "y": c[:, 1],
            "z": c[:, 2],
        })
        add_field_columns(df, data, fields, mask=mask)

        # Group duplicate x locations across y/z selection.
        df["x_round"] = df["x"].round(args.x_round_digits)
        agg = {}
        for col in df.columns:
            if col == "x_round":
                continue
            if col == "step":
                agg[col] = "first"
            else:
                agg[col] = "mean"
        df = df.groupby("x_round", as_index=False).agg(agg).sort_values("x")
        df = df.drop(columns=["x_round"])

        out = Path(args.out) if args.out and len(vtus) == 1 else out_path_for_snapshot(outdir, vtu, "centerline")
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out, index=False)
        print(f"Wrote: {out}  rows={len(df)} cols={len(df.columns)}")


def export_summary(args: argparse.Namespace) -> None:
    vtus = resolve_vtus(args)
    rows: list[dict[str, float | int | str]] = []

    for vtu in vtus:
        centers, data = read_snapshot(vtu)
        fields = field_names(data, args.fields)
        step = parse_step(vtu)

        for name in fields:
            arr = np.asarray(data[name])
            if arr.ndim == 2 and arr.shape[1] == 3:
                comps = {"x": arr[:, 0], "y": arr[:, 1], "z": arr[:, 2], "mag": np.linalg.norm(arr, axis=1)}
            elif arr.ndim == 1:
                comps = {"": arr}
            else:
                flat = arr.reshape(arr.shape[0], -1)
                comps = {str(j): flat[:, j] for j in range(flat.shape[1])}

            for comp, values in comps.items():
                a = np.asarray(values, dtype=float)
                finite = a[np.isfinite(a)]
                row = {
                    "vtu": str(vtu),
                    "step": step,
                    "field": name if not comp else f"{name}_{comp}",
                    "n": int(a.size),
                    "n_finite": int(finite.size),
                    "min": np.nan,
                    "max": np.nan,
                    "mean": np.nan,
                    "std": np.nan,
                    "l2": np.nan,
                }
                if finite.size:
                    row.update({
                        "min": float(np.min(finite)),
                        "max": float(np.max(finite)),
                        "mean": float(np.mean(finite)),
                        "std": float(np.std(finite)),
                        "l2": float(np.sqrt(np.mean(finite**2))),
                    })
                rows.append(row)

    df = pd.DataFrame(rows)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Wrote: {out}  rows={len(df)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Export LowMachReact-Hex VTU/PVTU output to CSV.")
    sub = parser.add_subparsers(dest="mode", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--vtu", default=None, help="Single .pvtu or .vtu file")
    common.add_argument("--vtu-glob", default=None, help="Glob for .pvtu/.vtu files, quote this argument")
    common.add_argument("--fields", default="all", help="Comma-separated field list, or all")

    p_cells = sub.add_parser("cells", parents=[common], help="Export all cell-centered data.")
    p_cells.add_argument("--out", default=None)
    p_cells.add_argument("--outdir", default=None)
    p_cells.set_defaults(func=export_cells)

    p_line = sub.add_parser("line", parents=[common], help="Export x-centerline data.")
    p_line.add_argument("--ymid", type=float, default=None)
    p_line.add_argument("--zmid", type=float, default=None)
    p_line.add_argument("--tol", type=float, default=None)
    p_line.add_argument("--x-round-digits", type=int, default=9)
    p_line.add_argument("--out", default=None)
    p_line.add_argument("--outdir", default=None)
    p_line.set_defaults(func=export_line)

    p_summary = sub.add_parser("summary", parents=[common], help="Export per-snapshot field statistics.")
    p_summary.add_argument("--out", required=True)
    p_summary.set_defaults(func=export_summary)

    args = parser.parse_args()
    if not args.vtu and not args.vtu_glob:
        raise SystemExit("Provide --vtu or --vtu-glob.")
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
