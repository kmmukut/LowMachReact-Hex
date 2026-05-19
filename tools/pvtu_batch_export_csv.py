#!/usr/bin/env python3
"""
Export every PVTU/VTU snapshot in a directory/glob to structured CSV files.

This is a solver-agnostic batch exporter for ParaView XML output. It scans all
matching `.pvtu` / `.vtu` files, reads the parallel pieces, and writes analysis-
ready CSV tables.

Default output is a WIDE table:

  cell_data_wide.csv

with one row per cell per snapshot:

  step, source, piece, local_id, global_id, x, y, z, <field columns...>

Vector/tensor fields are expanded into component columns, for example:

  velocity_x, velocity_y, velocity_z

It can also write a LONG table:

  cell_data_long.csv

with one row per scalar value:

  step, source, piece, local_id, global_id, x, y, z, field, component, value

The LONG table is easier for groupby/pivot workflows, but it is much larger.

It always writes a compact per-snapshot summary table:

  cell_data_summary.csv

with min/max/mean/std/L2 for every exported scalar component.

Requirements
------------
  pip install meshio numpy pandas

Examples
--------

Export all cell data from all PVTU snapshots in an output directory:

  python tools/pvtu_batch_export_csv.py \
    --input-glob 'cases/counterflow_nonreacting_gmsh/output/flow_*.pvtu' \
    --outdir cases/counterflow_nonreacting_gmsh/analysis_csv

Export selected fields only:

  python tools/pvtu_batch_export_csv.py \
    --input-glob 'cases/counterflow_nonreacting_gmsh/output/flow_*.pvtu' \
    --fields velocity,temperature,Y_CH4,Y_O2,Y_N2,rho,rho_thermo,nu,enthalpy,species_enthalpy_diffusion \
    --outdir cases/counterflow_nonreacting_gmsh/analysis_csv

Export both wide and long forms:

  python tools/pvtu_batch_export_csv.py \
    --input-glob 'cases/counterflow_nonreacting_gmsh/output/flow_*.pvtu' \
    --mode both \
    --outdir cases/counterflow_nonreacting_gmsh/analysis_csv

Export point data instead of cell data:

  python tools/pvtu_batch_export_csv.py \
    --input-glob 'cases/foo/output/flow_*.pvtu' \
    --data point \
    --outdir cases/foo/analysis_csv

Inspect fields without exporting:

  python tools/pvtu_batch_export_csv.py \
    --input-glob 'cases/counterflow_nonreacting_gmsh/output/flow_*.pvtu' \
    --list-fields
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


VECTOR_COMPONENT_NAMES = ("x", "y", "z")


def strip_namespace(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def parse_step(path: Path) -> int:
    """Parse timestep from names such as flow_005000.pvtu."""
    m = re.search(r"flow_(\d+)", path.name, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    m = re.search(r"(\d+)(?=\.(?:p)?vtu$)", path.name, flags=re.IGNORECASE)
    return int(m.group(1)) if m else -1


def collect_inputs(pattern: str) -> list[Path]:
    paths = [Path(p) for p in globlib.glob(pattern)]
    paths = sorted(set(paths), key=lambda p: (parse_step(p), str(p)))
    if not paths:
        raise SystemExit(f"No files matched --input-glob: {pattern}")
    return paths


def piece_paths_from_pvtu(path: Path) -> list[Path]:
    """Return piece files referenced by a PVTU file."""
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise RuntimeError(f"Could not parse PVTU XML: {path}") from exc

    pieces: list[Path] = []
    for elem in root.iter():
        if strip_namespace(elem.tag) == "Piece":
            source = elem.attrib.get("Source")
            if source:
                pieces.append((path.parent / source).resolve())

    if pieces:
        missing = [p for p in pieces if not p.exists()]
        if missing:
            raise RuntimeError(
                f"{path} references missing piece files:\n"
                + "\n".join(str(p) for p in missing[:20])
            )
        return pieces

    # Fallback for common LowMachReact-Hex naming.
    fallback = sorted(path.parent.glob(f"{path.stem}_P*.vtu"))
    if fallback:
        return [p.resolve() for p in fallback]

    raise RuntimeError(f"No pieces found for {path}")


def resolve_piece_paths(path: Path) -> list[Path]:
    suffix = path.suffix.lower()
    if suffix == ".pvtu":
        return piece_paths_from_pvtu(path)
    if suffix == ".vtu":
        return [path.resolve()]
    raise RuntimeError(f"Expected .pvtu or .vtu input, got {path}")


def choose_cell_block(mesh: meshio.Mesh, preferred: str | None = None) -> int:
    if preferred:
        for i, block in enumerate(mesh.cells):
            if block.type == preferred:
                return i
        raise RuntimeError(
            f"Requested cell block {preferred!r} not found. "
            f"Available blocks: {[b.type for b in mesh.cells]}"
        )

    preferred_order = ("hexahedron", "hex8", "quad", "triangle", "tetra", "wedge", "pyramid", "line")
    for cell_type in preferred_order:
        for i, block in enumerate(mesh.cells):
            if block.type == cell_type:
                return i

    if not mesh.cells:
        raise RuntimeError("mesh has no cells")
    return 0


def coordinates_and_data(
    mesh: meshio.Mesh,
    data_kind: str,
    cell_type: str | None,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if data_kind == "cell":
        bi = choose_cell_block(mesh, preferred=cell_type)
        cells = mesh.cells[bi].data
        coords = mesh.points[cells].mean(axis=1)
        data = {name: np.asarray(blocks[bi]) for name, blocks in mesh.cell_data.items() if len(blocks) > bi}
        return coords, data

    if data_kind == "point":
        coords = np.asarray(mesh.points)
        data = {name: np.asarray(arr) for name, arr in mesh.point_data.items()}
        return coords, data

    raise ValueError(f"unknown data kind: {data_kind}")


def parse_fields(fields: str | None, available: Iterable[str]) -> list[str]:
    available_set = set(available)
    available_sorted = sorted(available_set)

    if fields is None or fields.strip().lower() == "all":
        return available_sorted

    requested = [x.strip() for x in fields.split(",") if x.strip()]
    missing = [x for x in requested if x not in available_set]
    if missing:
        print(f"WARNING: missing fields skipped: {', '.join(missing)}", file=sys.stderr)

    return [x for x in requested if x in available_set]


def expanded_columns(name: str, arr: np.ndarray) -> dict[str, np.ndarray]:
    arr = np.asarray(arr)

    if arr.ndim == 1:
        return {name: arr}

    if arr.ndim == 2 and arr.shape[1] == 3:
        return {f"{name}_{comp}": arr[:, j] for j, comp in enumerate(VECTOR_COMPONENT_NAMES)}

    flat = arr.reshape(arr.shape[0], -1)
    return {f"{name}_{j}": flat[:, j] for j in range(flat.shape[1])}


def component_arrays(name: str, arr: np.ndarray) -> dict[str, np.ndarray]:
    """Return component name -> values for summary/long format."""
    arr = np.asarray(arr)

    if arr.ndim == 1:
        return {"": arr}

    if arr.ndim == 2 and arr.shape[1] == 3:
        out = {comp: arr[:, j] for j, comp in enumerate(VECTOR_COMPONENT_NAMES)}
        out["mag"] = np.linalg.norm(arr, axis=1)
        return out

    flat = arr.reshape(arr.shape[0], -1)
    return {str(j): flat[:, j] for j in range(flat.shape[1])}


def summarize(values: np.ndarray) -> dict[str, float | int]:
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


def read_piece(piece: Path, data_kind: str, cell_type: str | None) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    try:
        mesh = meshio.read(str(piece))
    except Exception as exc:
        raise RuntimeError(f"meshio failed reading {piece}: {exc}") from exc
    return coordinates_and_data(mesh, data_kind=data_kind, cell_type=cell_type)


def list_fields(paths: list[Path], data_kind: str, cell_type: str | None) -> None:
    first = paths[0]
    print(f"First snapshot: {first}")
    print(f"Snapshots matched: {len(paths)}")

    all_fields: dict[str, tuple[tuple[int, ...], str]] = {}
    for piece in resolve_piece_paths(first):
        coords, data = read_piece(piece, data_kind=data_kind, cell_type=cell_type)
        for name, arr in data.items():
            all_fields[name] = (tuple(arr.shape), str(arr.dtype))

    print(f"\n{data_kind} fields:")
    for name in sorted(all_fields):
        shape, dtype = all_fields[name]
        print(f"  {name}: shape={shape} dtype={dtype}")


def write_wide_rows(
    out_path: Path,
    snapshot: Path,
    piece: Path,
    piece_index: int,
    global_offset: int,
    coords: np.ndarray,
    data: dict[str, np.ndarray],
    fields: list[str],
    append: bool,
) -> int:
    n = coords.shape[0]
    df_dict: dict[str, object] = {
        "step": np.full(n, parse_step(snapshot), dtype=int),
        "source": np.full(n, snapshot.name, dtype=object),
        "piece": np.full(n, piece.name, dtype=object),
        "piece_index": np.full(n, piece_index, dtype=int),
        "local_id": np.arange(n, dtype=int),
        "global_id": np.arange(global_offset, global_offset + n, dtype=int),
        "x": coords[:, 0],
        "y": coords[:, 1],
        "z": coords[:, 2],
    }

    for field in fields:
        for col, values in expanded_columns(field, data[field]).items():
            df_dict[col] = values

    df = pd.DataFrame(df_dict)
    df.to_csv(out_path, mode="a" if append else "w", header=not append, index=False)
    return n


def write_long_rows(
    out_path: Path,
    snapshot: Path,
    piece: Path,
    piece_index: int,
    global_offset: int,
    coords: np.ndarray,
    data: dict[str, np.ndarray],
    fields: list[str],
    append: bool,
) -> int:
    rows_written = 0
    first_chunk = not append

    base = pd.DataFrame(
        {
            "step": np.full(coords.shape[0], parse_step(snapshot), dtype=int),
            "source": np.full(coords.shape[0], snapshot.name, dtype=object),
            "piece": np.full(coords.shape[0], piece.name, dtype=object),
            "piece_index": np.full(coords.shape[0], piece_index, dtype=int),
            "local_id": np.arange(coords.shape[0], dtype=int),
            "global_id": np.arange(global_offset, global_offset + coords.shape[0], dtype=int),
            "x": coords[:, 0],
            "y": coords[:, 1],
            "z": coords[:, 2],
        }
    )

    for field in fields:
        for component, values in component_arrays(field, data[field]).items():
            chunk = base.copy()
            chunk["field"] = field
            chunk["component"] = component
            chunk["value"] = np.asarray(values).reshape(-1)
            chunk.to_csv(out_path, mode="a", header=first_chunk, index=False)
            first_chunk = False
            rows_written += len(chunk)

    return rows_written


def append_summary_rows(
    rows: list[dict[str, object]],
    snapshot: Path,
    piece: Path,
    piece_index: int,
    data_kind: str,
    data: dict[str, np.ndarray],
    fields: list[str],
) -> None:
    step = parse_step(snapshot)

    for field in fields:
        for component, values in component_arrays(field, data[field]).items():
            row: dict[str, object] = {
                "step": step,
                "source": snapshot.name,
                "piece": piece.name,
                "piece_index": piece_index,
                "data": data_kind,
                "field": field,
                "component": component,
            }
            row.update(summarize(values))
            rows.append(row)


def export_batch(args: argparse.Namespace) -> None:
    paths = collect_inputs(args.input_glob)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.list_fields:
        list_fields(paths, args.data, args.cell_type)
        return

    wide_path = outdir / f"{args.prefix}_wide.csv"
    long_path = outdir / f"{args.prefix}_long.csv"
    summary_path = outdir / f"{args.prefix}_summary.csv"
    manifest_path = outdir / f"{args.prefix}_manifest.csv"

    for path in (wide_path, long_path, summary_path, manifest_path):
        if path.exists() and not args.append:
            path.unlink()

    write_wide = args.mode in ("wide", "both")
    write_long = args.mode in ("long", "both")

    wide_append = wide_path.exists()
    long_append = long_path.exists()

    summary_rows: list[dict[str, object]] = []
    manifest_rows: list[dict[str, object]] = []

    total_wide_rows = 0
    total_long_rows = 0

    for snapshot_index, snapshot in enumerate(paths):
        pieces = resolve_piece_paths(snapshot)
        print(f"[{snapshot_index + 1}/{len(paths)}] {snapshot} pieces={len(pieces)}")

        global_offset = 0
        snapshot_rows = 0
        available_fields_for_snapshot: set[str] = set()

        for piece_index, piece in enumerate(pieces):
            coords, data = read_piece(piece, data_kind=args.data, cell_type=args.cell_type)
            fields = parse_fields(args.fields, data.keys())
            available_fields_for_snapshot.update(data.keys())

            if write_wide:
                n = write_wide_rows(
                    out_path=wide_path,
                    snapshot=snapshot,
                    piece=piece,
                    piece_index=piece_index,
                    global_offset=global_offset,
                    coords=coords,
                    data=data,
                    fields=fields,
                    append=wide_append,
                )
                wide_append = True
                total_wide_rows += n

            if write_long:
                n_long = write_long_rows(
                    out_path=long_path,
                    snapshot=snapshot,
                    piece=piece,
                    piece_index=piece_index,
                    global_offset=global_offset,
                    coords=coords,
                    data=data,
                    fields=fields,
                    append=long_append,
                )
                long_append = True
                total_long_rows += n_long

            append_summary_rows(
                rows=summary_rows,
                snapshot=snapshot,
                piece=piece,
                piece_index=piece_index,
                data_kind=args.data,
                data=data,
                fields=fields,
            )

            n_cells_or_points = int(coords.shape[0])
            snapshot_rows += n_cells_or_points
            global_offset += n_cells_or_points

        manifest_rows.append(
            {
                "snapshot_index": snapshot_index,
                "step": parse_step(snapshot),
                "source": snapshot.name,
                "path": str(snapshot),
                "n_pieces": len(pieces),
                "n_rows": snapshot_rows,
                "data": args.data,
                "fields_available": ",".join(sorted(available_fields_for_snapshot)),
                "fields_exported": args.fields,
            }
        )

        # Flush summary incrementally so long runs still leave usable output.
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
        pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)

    print()
    if write_wide:
        print(f"WIDE CSV:    {wide_path} rows={total_wide_rows}")
    if write_long:
        print(f"LONG CSV:    {long_path} rows={total_long_rows}")
    print(f"SUMMARY CSV: {summary_path} rows={len(summary_rows)}")
    print(f"MANIFEST:    {manifest_path} rows={len(manifest_rows)}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch export PVTU/VTU snapshots to structured CSV.")
    p.add_argument(
        "--input-glob",
        required=True,
        help="Quoted glob for .pvtu/.vtu files, e.g. 'cases/foo/output/flow_*.pvtu'",
    )
    p.add_argument("--outdir", default="pvtu_csv", help="Output directory for CSV files.")
    p.add_argument(
        "--prefix",
        default="cell_data",
        help="Output file prefix. Defaults to cell_data, producing cell_data_wide.csv etc.",
    )
    p.add_argument("--data", choices=["cell", "point"], default="cell", help="Export cell or point data.")
    p.add_argument(
        "--fields",
        default="all",
        help="Comma-separated field names, or all. Use --list-fields to inspect names.",
    )
    p.add_argument("--cell-type", default=None, help="Optional meshio cell block type, e.g. hexahedron.")
    p.add_argument(
        "--mode",
        choices=["wide", "long", "both"],
        default="wide",
        help="CSV shape. wide is compact and analysis-friendly; long is tidy but much larger.",
    )
    p.add_argument("--append", action="store_true", help="Append to existing CSV files instead of replacing them.")
    p.add_argument("--list-fields", action="store_true", help="List fields from the first snapshot and exit.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    export_batch(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

