#!/usr/bin/env python3
"""
Generic PVTU/VTU to CSV exporter.

This tool reads serial `.vtu` files or parallel `.pvtu` files and exports
VTK data arrays to CSV without opening ParaView.

It is intentionally solver-agnostic:
  - works with any `.pvtu` whose pieces are readable by meshio
  - exports cell data, point data, or both
  - can list available arrays
  - can export all fields or selected fields
  - can summarize many snapshots from a glob

Requirements:
  pip install meshio numpy pandas

Common examples
---------------

List available fields:

  python tools/pvtu_export.py list \
    --input cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu

Export all cell data for one snapshot:

  python tools/pvtu_export.py export \
    --input cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu \
    --data cell \
    --out cases/counterflow_nonreacting_gmsh/validation/flow_005000_cells.csv

Export selected cell fields:

  python tools/pvtu_export.py export \
    --input cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu \
    --data cell \
    --fields velocity,temperature,Y_CH4,Y_O2,Y_N2,rho,rho_thermo,nu,enthalpy \
    --out cases/counterflow_nonreacting_gmsh/validation/flow_005000_selected.csv

Export all points/point-data:

  python tools/pvtu_export.py export \
    --input cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu \
    --data point \
    --out cases/counterflow_nonreacting_gmsh/validation/flow_005000_points.csv

Export centerline cells through the domain midpoint in y/z:

  python tools/pvtu_export.py line \
    --input cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu \
    --axis x \
    --out cases/counterflow_nonreacting_gmsh/validation/flow_005000_centerline.csv

Summarize all snapshots:

  python tools/pvtu_export.py summary \
    --glob 'cases/counterflow_nonreacting_gmsh/output/flow_*.pvtu' \
    --data cell \
    --fields temperature,Y_CH4,Y_O2,Y_N2,rho_thermo,enthalpy,species_enthalpy_diffusion \
    --out cases/counterflow_nonreacting_gmsh/validation/summary.csv
"""

from __future__ import annotations

import argparse
import glob as globlib
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable

import meshio
import numpy as np
import pandas as pd


AXIS_TO_INDEX = {"x": 0, "y": 1, "z": 2}
VECTOR_COMPONENT_NAMES = ("x", "y", "z")


def parse_step(path: Path) -> int:
    match = re.search(r"(\d+)(?=\.(?:p)?vtu$)", path.name, flags=re.IGNORECASE)
    if match:
        return int(match.group(1))
    match = re.search(r"flow_(\d+)", path.name, flags=re.IGNORECASE)
    return int(match.group(1)) if match else -1


def collect_inputs(input_path: str | None, glob_pattern: str | None) -> list[Path]:
    paths: list[Path] = []

    if input_path:
        paths.append(Path(input_path))

    if glob_pattern:
        paths.extend(Path(p) for p in globlib.glob(glob_pattern))

    paths = sorted(set(paths), key=lambda p: (parse_step(p), str(p)))
    if not paths:
        raise SystemExit("No input files found. Provide --input and/or --glob.")

    for path in paths:
        if not path.exists():
            raise SystemExit(f"Input file does not exist: {path}")

    return paths


def strip_xml_namespace(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def piece_paths_from_pvtu(pvtu_path: Path) -> list[Path]:
    """
    Read piece paths from a .pvtu file. Fall back to stem_P*.vtu if the XML does
    not contain Piece Source entries.
    """
    try:
        root = ET.parse(pvtu_path).getroot()
        pieces: list[Path] = []
        for elem in root.iter():
            if strip_xml_namespace(elem.tag) == "Piece":
                source = elem.attrib.get("Source")
                if source:
                    pieces.append((pvtu_path.parent / source).resolve())
        if pieces:
            missing = [p for p in pieces if not p.exists()]
            if missing:
                raise RuntimeError(
                    "PVTU references missing piece files:\n"
                    + "\n".join(str(p) for p in missing[:20])
                )
            return pieces
    except ET.ParseError as exc:
        raise RuntimeError(f"Could not parse PVTU XML: {pvtu_path}") from exc

    fallback = sorted(pvtu_path.parent.glob(f"{pvtu_path.stem}_P*.vtu"))
    if fallback:
        return [p.resolve() for p in fallback]

    raise RuntimeError(f"No piece files found for {pvtu_path}")


def resolve_piece_paths(path: Path) -> list[Path]:
    suffix = path.suffix.lower()
    if suffix == ".pvtu":
        return piece_paths_from_pvtu(path)
    if suffix == ".vtu":
        return [path.resolve()]
    raise RuntimeError(f"Expected .vtu or .pvtu file, got: {path}")


def choose_cell_block(mesh: meshio.Mesh, preferred: str | None = None) -> int:
    if preferred:
        for i, block in enumerate(mesh.cells):
            if block.type == preferred:
                return i
        raise RuntimeError(
            f"Requested cell block type {preferred!r} not found. "
            f"Available: {[b.type for b in mesh.cells]}"
        )

    for target in ("hexahedron", "hex8", "quad", "triangle", "tetra", "wedge", "pyramid"):
        for i, block in enumerate(mesh.cells):
            if block.type == target:
                return i

    if not mesh.cells:
        raise RuntimeError("Mesh has no cells.")
    return 0


def cell_centers(mesh: meshio.Mesh, block_index: int) -> np.ndarray:
    cells = mesh.cells[block_index].data
    return mesh.points[cells].mean(axis=1)


def selected_cell_data(mesh: meshio.Mesh, block_index: int) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    for name, blocks in mesh.cell_data.items():
        if len(blocks) > block_index:
            out[name] = np.asarray(blocks[block_index])
    return out


def selected_point_data(mesh: meshio.Mesh) -> dict[str, np.ndarray]:
    return {name: np.asarray(arr) for name, arr in mesh.point_data.items()}


def merge_snapshot(path: Path, data_kind: str, cell_type: str | None = None) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """
    Return coordinates and merged data arrays for a snapshot.

    For cell data, coordinates are cell centers.
    For point data, coordinates are mesh points.
    """
    coords_parts: list[np.ndarray] = []
    data_parts: dict[str, list[np.ndarray]] = {}

    for piece in resolve_piece_paths(path):
        try:
            mesh = meshio.read(str(piece))
        except Exception as exc:
            raise RuntimeError(f"meshio failed while reading {piece}: {exc}") from exc

        if data_kind == "cell":
            bi = choose_cell_block(mesh, preferred=cell_type)
            coords = cell_centers(mesh, bi)
            data = selected_cell_data(mesh, bi)
        elif data_kind == "point":
            coords = np.asarray(mesh.points)
            data = selected_point_data(mesh)
        else:
            raise ValueError(f"unknown data kind: {data_kind}")

        coords_parts.append(coords)
        for name, arr in data.items():
            data_parts.setdefault(name, []).append(np.asarray(arr))

    coords_merged = np.vstack(coords_parts) if coords_parts else np.empty((0, 3))
    data_merged: dict[str, np.ndarray] = {}

    for name, parts in data_parts.items():
        try:
            data_merged[name] = np.concatenate(parts, axis=0)
        except Exception as exc:
            shapes = [p.shape for p in parts]
            raise RuntimeError(f"Could not concatenate field {name!r}; shapes={shapes}") from exc

    return coords_merged, data_merged


def parse_fields(fields_arg: str | None, available: Iterable[str]) -> list[str]:
    available_list = sorted(available)

    if not fields_arg or fields_arg.lower() == "all":
        return available_list

    requested = [field.strip() for field in fields_arg.split(",") if field.strip()]
    missing = [field for field in requested if field not in available_list]
    if missing:
        print(f"WARNING: missing fields skipped: {', '.join(missing)}", file=sys.stderr)

    return [field for field in requested if field in available_list]


def add_array_to_dataframe(df: pd.DataFrame, name: str, arr: np.ndarray) -> None:
    arr = np.asarray(arr)

    if arr.ndim == 1:
        df[name] = arr
        return

    if arr.ndim == 2 and arr.shape[1] == 3:
        for j, comp in enumerate(VECTOR_COMPONENT_NAMES):
            df[f"{name}_{comp}"] = arr[:, j]
        return

    flat = arr.reshape(arr.shape[0], -1)
    for j in range(flat.shape[1]):
        df[f"{name}_{j}"] = flat[:, j]


def write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    print(f"Wrote: {path}  rows={len(df)} cols={len(df.columns)}")


def output_path_for_snapshot(outdir: Path, input_path: Path, data_kind: str, suffix: str) -> Path:
    step = parse_step(input_path)
    if step >= 0:
        return outdir / f"{input_path.stem}_{data_kind}_{suffix}_{step:06d}.csv"
    return outdir / f"{input_path.stem}_{data_kind}_{suffix}.csv"


def command_list(args: argparse.Namespace) -> None:
    paths = collect_inputs(args.input, args.glob)
    for path in paths:
        print(f"\n{path}")
        print("-" * len(str(path)))

        for kind in ("cell", "point"):
            try:
                coords, data = merge_snapshot(path, kind, cell_type=args.cell_type)
            except Exception as exc:
                print(f"{kind} data: ERROR: {exc}")
                continue

            print(f"{kind} coordinates: {coords.shape}")
            if not data:
                print(f"{kind} fields: <none>")
                continue

            print(f"{kind} fields:")
            for name in sorted(data):
                print(f"  {name}: shape={data[name].shape} dtype={data[name].dtype}")


def command_export(args: argparse.Namespace) -> None:
    paths = collect_inputs(args.input, args.glob)

    if len(paths) > 1 and not args.outdir:
        raise SystemExit("--outdir is required when exporting multiple input snapshots.")

    for path in paths:
        coords, data = merge_snapshot(path, args.data, cell_type=args.cell_type)
        fields = parse_fields(args.fields, data.keys())

        df = pd.DataFrame(
            {
                "source": path.name,
                "step": parse_step(path),
                "id": np.arange(coords.shape[0], dtype=int),
                "x": coords[:, 0],
                "y": coords[:, 1],
                "z": coords[:, 2],
            }
        )

        for field in fields:
            add_array_to_dataframe(df, field, data[field])

        if args.data == "point" and args.drop_duplicate_points:
            before = len(df)
            df = df.drop_duplicates(subset=["x", "y", "z"]).reset_index(drop=True)
            print(f"Dropped duplicate points: {before} -> {len(df)}")

        if args.out and len(paths) == 1:
            out = Path(args.out)
        else:
            outdir = Path(args.outdir)
            out = output_path_for_snapshot(outdir, path, args.data, "export")

        write_csv(df, out)


def estimate_axis_tolerance(coords: np.ndarray, axis: int) -> float:
    vals = np.unique(np.round(coords[:, axis], 12))
    vals = np.sort(vals)
    if vals.size < 2:
        return 0.05
    span = float(vals[-1] - vals[0])
    diffs = np.diff(vals)
    diffs = diffs[diffs > max(1.0e-12, 1.0e-8 * max(span, 1.0))]
    if diffs.size == 0:
        return 0.05
    return 0.60 * float(np.median(diffs))


def command_line(args: argparse.Namespace) -> None:
    paths = collect_inputs(args.input, args.glob)

    if len(paths) > 1 and not args.outdir:
        raise SystemExit("--outdir is required when exporting multiple input snapshots.")

    axis = AXIS_TO_INDEX[args.axis]
    fixed_axes = [i for i in range(3) if i != axis]

    for path in paths:
        coords, data = merge_snapshot(path, "cell", cell_type=args.cell_type)
        fields = parse_fields(args.fields, data.keys())

        mins = coords.min(axis=0)
        maxs = coords.max(axis=0)
        mids = 0.5 * (mins + maxs)

        fixed_values = {
            0: args.x if args.x is not None else mids[0],
            1: args.y if args.y is not None else mids[1],
            2: args.z if args.z is not None else mids[2],
        }

        tol = args.tol
        if tol is None:
            tol = max(estimate_axis_tolerance(coords, i) for i in fixed_axes)

        mask = np.ones(coords.shape[0], dtype=bool)
        for fixed_axis in fixed_axes:
            mask &= np.abs(coords[:, fixed_axis] - fixed_values[fixed_axis]) <= tol

        if not np.any(mask):
            fixed_desc = ", ".join(f"{'xyz'[i]}={fixed_values[i]}" for i in fixed_axes)
            raise RuntimeError(f"No cells found near {fixed_desc} with tol={tol} in {path}")

        c = coords[mask]
        df = pd.DataFrame(
            {
                "source": path.name,
                "step": parse_step(path),
                "id": np.arange(c.shape[0], dtype=int),
                "x": c[:, 0],
                "y": c[:, 1],
                "z": c[:, 2],
            }
        )

        for field in fields:
            add_array_to_dataframe(df, field, np.asarray(data[field])[mask])

        coord_col = args.axis
        round_col = f"{coord_col}_round"
        df[round_col] = df[coord_col].round(args.round_digits)

        agg: dict[str, str] = {}
        for col in df.columns:
            if col == round_col:
                continue
            if col in ("source",):
                agg[col] = "first"
            elif col in ("step", "id"):
                agg[col] = "first"
            else:
                agg[col] = "mean"

        df = df.groupby(round_col, as_index=False).agg(agg).sort_values(coord_col)
        df = df.drop(columns=[round_col])

        if args.out and len(paths) == 1:
            out = Path(args.out)
        else:
            outdir = Path(args.outdir)
            out = output_path_for_snapshot(outdir, path, "cell", f"line_{args.axis}")

        write_csv(df, out)


def summarize_array(values: np.ndarray) -> dict[str, float | int]:
    arr = np.asarray(values, dtype=float).reshape(-1)
    finite = arr[np.isfinite(arr)]
    out: dict[str, float | int] = {
        "n": int(arr.size),
        "n_finite": int(finite.size),
        "min": np.nan,
        "max": np.nan,
        "mean": np.nan,
        "std": np.nan,
        "l2": np.nan,
    }
    if finite.size:
        out.update(
            {
                "min": float(np.min(finite)),
                "max": float(np.max(finite)),
                "mean": float(np.mean(finite)),
                "std": float(np.std(finite)),
                "l2": float(np.sqrt(np.mean(finite**2))),
            }
        )
    return out


def command_summary(args: argparse.Namespace) -> None:
    paths = collect_inputs(args.input, args.glob)
    rows: list[dict[str, object]] = []

    for path in paths:
        coords, data = merge_snapshot(path, args.data, cell_type=args.cell_type)
        fields = parse_fields(args.fields, data.keys())

        for field in fields:
            arr = np.asarray(data[field])
            if arr.ndim == 1:
                components = {field: arr}
            elif arr.ndim == 2 and arr.shape[1] == 3:
                components = {
                    f"{field}_x": arr[:, 0],
                    f"{field}_y": arr[:, 1],
                    f"{field}_z": arr[:, 2],
                    f"{field}_mag": np.linalg.norm(arr, axis=1),
                }
            else:
                flat = arr.reshape(arr.shape[0], -1)
                components = {f"{field}_{j}": flat[:, j] for j in range(flat.shape[1])}

            for comp_name, values in components.items():
                row: dict[str, object] = {
                    "source": path.name,
                    "step": parse_step(path),
                    "data": args.data,
                    "field": comp_name,
                }
                row.update(summarize_array(values))
                rows.append(row)

    out = Path(args.out)
    write_csv(pd.DataFrame(rows), out)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generic PVTU/VTU to CSV exporter.")
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--input", default=None, help="Input .pvtu or .vtu file")
    common.add_argument("--glob", default=None, help="Glob pattern for input files. Quote this argument.")
    common.add_argument("--cell-type", default=None, help="Optional meshio cell block type, e.g. hexahedron")

    p_list = sub.add_parser("list", parents=[common], help="List available arrays.")
    p_list.set_defaults(func=command_list)

    p_export = sub.add_parser("export", parents=[common], help="Export cell or point data to CSV.")
    p_export.add_argument("--data", choices=["cell", "point"], default="cell")
    p_export.add_argument("--fields", default="all", help="Comma-separated fields, or all")
    p_export.add_argument("--out", default=None, help="Output CSV for a single input")
    p_export.add_argument("--outdir", default=None, help="Output directory for multiple inputs")
    p_export.add_argument("--drop-duplicate-points", action="store_true")
    p_export.set_defaults(func=command_export)

    p_line = sub.add_parser("line", parents=[common], help="Export a cell-center line cut to CSV.")
    p_line.add_argument("--axis", choices=["x", "y", "z"], default="x", help="Varying coordinate")
    p_line.add_argument("--x", type=float, default=None, help="Fixed x, unless axis=x")
    p_line.add_argument("--y", type=float, default=None, help="Fixed y, unless axis=y")
    p_line.add_argument("--z", type=float, default=None, help="Fixed z, unless axis=z")
    p_line.add_argument("--tol", type=float, default=None, help="Tolerance for fixed coordinates")
    p_line.add_argument("--round-digits", type=int, default=9)
    p_line.add_argument("--fields", default="all", help="Comma-separated fields, or all")
    p_line.add_argument("--out", default=None, help="Output CSV for a single input")
    p_line.add_argument("--outdir", default=None, help="Output directory for multiple inputs")
    p_line.set_defaults(func=command_line)

    p_summary = sub.add_parser("summary", parents=[common], help="Export summary statistics to CSV.")
    p_summary.add_argument("--data", choices=["cell", "point"], default="cell")
    p_summary.add_argument("--fields", default="all", help="Comma-separated fields, or all")
    p_summary.add_argument("--out", required=True)
    p_summary.set_defaults(func=command_summary)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not args.input and not args.glob:
        raise SystemExit("Provide --input and/or --glob.")

    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
