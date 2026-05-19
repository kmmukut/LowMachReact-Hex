#!/usr/bin/env python3
"""Convert a Gmsh/meshio hexahedral mesh to the native v1 text format.

The native format is intentionally boring so the Fortran reader stays small:

  points.dat   : npoints, then "id x y z"
  cells.dat    : ncells, then "id n1..n8 cx cy cz volume"
  faces.dat    : nfaces, then "id owner neighbor patch nx ny nz area cx cy cz"
  patches.dat  : npatches, then "id name nfaces" and one face-id line
  periodic.dat : optional "nlinks", then "face_id pair_face_id neighbor_cell"

Only linear, axis-aligned hexahedral cuboids are accepted.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import meshio
import numpy as np


TOL = 1.0e-10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("msh", type=Path, help="Input .msh file")
    parser.add_argument("outdir", type=Path, help="Output native mesh directory")
    parser.add_argument(
        "--periodic",
        action="append",
        default=[],
        metavar="PATCH_A:PATCH_B",
        help="Pair two boundary patches periodically. May be used multiple times.",
    )
    return parser.parse_args()


def field_name_by_tag(mesh: meshio.Mesh, dim: int) -> Dict[int, str]:
    out: Dict[int, str] = {}
    for name, data in mesh.field_data.items():
        tag = int(data[0])
        entity_dim = int(data[1])
        if entity_dim == dim:
            out[tag] = name
    return out


def physical_tags_for_block(mesh: meshio.Mesh, block_index: int) -> np.ndarray | None:
    if "gmsh:physical" not in mesh.cell_data:
        return None
    return np.asarray(mesh.cell_data["gmsh:physical"][block_index], dtype=int)


def collect_cells(mesh: meshio.Mesh, cell_type: str) -> Tuple[np.ndarray, np.ndarray | None]:
    blocks: List[np.ndarray] = []
    tags: List[np.ndarray] = []

    for ib, block in enumerate(mesh.cells):
        if block.type == cell_type:
            data = np.asarray(block.data, dtype=int)
            blocks.append(data)
            phys = physical_tags_for_block(mesh, ib)
            if phys is not None:
                tags.append(phys)

    if not blocks:
        return np.empty((0, 0), dtype=int), None

    all_cells = np.vstack(blocks)
    all_tags = np.concatenate(tags) if tags else None
    return all_cells, all_tags


def validate_cuboid(points: np.ndarray, node_ids: Iterable[int]) -> Tuple[np.ndarray, float]:
    xyz = points[list(node_ids), :3]
    mins = xyz.min(axis=0)
    maxs = xyz.max(axis=0)
    lengths = maxs - mins
    tol = max(TOL, 1.0e-9 * float(np.max(lengths)))

    if np.any(lengths <= TOL):
        raise ValueError("Degenerate hexahedron with non-positive extent.")

    corner_bits = set()
    for row in xyz:
        bits = []
        for axis in range(3):
            if abs(row[axis] - mins[axis]) <= tol:
                bits.append(0)
            elif abs(row[axis] - maxs[axis]) <= tol:
                bits.append(1)
            else:
                raise ValueError("Hexahedron node is not on a cuboid min/max plane.")
        corner_bits.add(tuple(bits))

    if len(corner_bits) != 8:
        raise ValueError("Non-axis-aligned or non-rectangular hexahedron found.")

    return 0.5 * (mins + maxs), float(np.prod(lengths))


def face_specs(points: np.ndarray, node_ids: np.ndarray) -> List[Tuple[str, Tuple[int, ...], np.ndarray, float, np.ndarray]]:
    xyz = points[node_ids, :3]
    mins = xyz.min(axis=0)
    maxs = xyz.max(axis=0)
    center = 0.5 * (mins + maxs)
    dx, dy, dz = maxs - mins
    tol = max(TOL, 1.0e-9 * float(max(dx, dy, dz)))

    specs = []
    for axis, side, normal, area in [
        (0, mins[0], np.array([-1.0, 0.0, 0.0]), dy * dz),
        (0, maxs[0], np.array([+1.0, 0.0, 0.0]), dy * dz),
        (1, mins[1], np.array([0.0, -1.0, 0.0]), dx * dz),
        (1, maxs[1], np.array([0.0, +1.0, 0.0]), dx * dz),
        (2, mins[2], np.array([0.0, 0.0, -1.0]), dx * dy),
        (2, maxs[2], np.array([0.0, 0.0, +1.0]), dx * dy),
    ]:
        mask = np.abs(xyz[:, axis] - side) <= tol
        ids = tuple(sorted(int(n) for n in node_ids[mask]))
        fcenter = center.copy()
        fcenter[axis] = side
        specs.append((f"axis{axis}", ids, normal, float(area), fcenter))

    return specs


def collect_boundary_patch_names(mesh: meshio.Mesh) -> Dict[frozenset[int], str]:
    quad_cells, quad_tags = collect_cells(mesh, "quad")
    names_by_tag = field_name_by_tag(mesh, 2)
    patch_by_nodes: Dict[frozenset[int], str] = {}

    if quad_cells.size == 0:
        return patch_by_nodes

    if quad_tags is None:
        raise ValueError("Boundary quad cells do not have gmsh physical tags.")

    for conn, tag in zip(quad_cells, quad_tags):
        name = names_by_tag.get(int(tag), f"patch_{int(tag)}")
        patch_by_nodes[frozenset(int(i) for i in conn)] = name

    return patch_by_nodes


def face_match_key(center: np.ndarray, normal: np.ndarray) -> Tuple[int, int, Tuple[int, int]]:
    axis = int(np.argmax(np.abs(normal)))
    other = [i for i in range(3) if i != axis]
    scaled = tuple(int(round(center[i] / TOL)) for i in other)
    return axis, int(np.sign(normal[axis])), scaled


def pair_periodic_faces(
    faces: List[dict],
    patch_name_to_id: Dict[str, int],
    pairs: List[str],
) -> List[Tuple[int, int, int]]:
    links: List[Tuple[int, int, int]] = []

    by_patch: Dict[int, List[dict]] = defaultdict(list)
    for face in faces:
        if face["neighbor"] == 0:
            by_patch[face["patch"]].append(face)

    for spec in pairs:
        if ":" not in spec:
            raise ValueError(f"Invalid --periodic spec {spec!r}; expected A:B")
        a_name, b_name = [part.strip() for part in spec.split(":", 1)]
        if a_name not in patch_name_to_id or b_name not in patch_name_to_id:
            raise ValueError(f"Periodic patches {a_name!r}:{b_name!r} not found in mesh patches.")

        a_faces = by_patch[patch_name_to_id[a_name]]
        b_faces = by_patch[patch_name_to_id[b_name]]

        if len(a_faces) != len(b_faces):
            raise ValueError(f"Periodic patches {a_name}:{b_name} have different face counts.")

        b_by_key = {}
        for face in b_faces:
            axis = int(np.argmax(np.abs(face["normal"])))
            other = [i for i in range(3) if i != axis]
            key = tuple(int(round(face["center"][i] / TOL)) for i in other)
            b_by_key[key] = face

        for face_a in a_faces:
            axis = int(np.argmax(np.abs(face_a["normal"])))
            other = [i for i in range(3) if i != axis]
            key = tuple(int(round(face_a["center"][i] / TOL)) for i in other)
            face_b = b_by_key.get(key)
            if face_b is None:
                raise ValueError(f"Could not find periodic mate for face {face_a['id']}.")
            links.append((face_a["id"], face_b["id"], face_b["owner"]))
            links.append((face_b["id"], face_a["id"], face_a["owner"]))

    return links


def main() -> None:
    args = parse_args()
    mesh = meshio.read(args.msh)
    points = np.asarray(mesh.points[:, :3], dtype=float)

    hex_cells, _ = collect_cells(mesh, "hexahedron")
    if hex_cells.size == 0:
        raise SystemExit("No linear hexahedron cells found.")

    patch_by_nodes = collect_boundary_patch_names(mesh)
    patch_name_to_id: Dict[str, int] = {}
    patch_faces: Dict[int, List[int]] = defaultdict(list)

    cells = []
    face_by_nodes: Dict[frozenset[int], dict] = {}
    faces: List[dict] = []

    for icell, conn in enumerate(hex_cells, start=1):
        center, volume = validate_cuboid(points, conn)
        cells.append(
            {
                "id": icell,
                "nodes": [int(i) + 1 for i in conn],
                "center": center,
                "volume": volume,
            }
        )

        for _, face_nodes0, normal, area, fcenter in face_specs(points, conn):
            key = frozenset(face_nodes0)
            if key in face_by_nodes:
                face = face_by_nodes[key]
                if face["neighbor"] != 0:
                    raise ValueError("Face shared by more than two cells.")
                face["neighbor"] = icell
            else:
                patch_name = patch_by_nodes.get(key, "interior")
                patch_id = 0
                if patch_name != "interior":
                    patch_id = patch_name_to_id.setdefault(patch_name, len(patch_name_to_id) + 1)

                face = {
                    "id": len(faces) + 1,
                    "owner": icell,
                    "neighbor": 0,
                    "patch": patch_id,
                    "normal": normal,
                    "area": area,
                    "center": fcenter,
                }
                faces.append(face)
                face_by_nodes[key] = face

    for face in faces:
        if face["neighbor"] == 0:
            if face["patch"] == 0:
                raise ValueError(f"Boundary face {face['id']} has no physical patch.")
            patch_faces[face["patch"]].append(face["id"])
        else:
            face["patch"] = 0

    periodic_links = pair_periodic_faces(faces, patch_name_to_id, args.periodic)

    args.outdir.mkdir(parents=True, exist_ok=True)

    with (args.outdir / "points.dat").open("w", encoding="utf-8") as f:
        f.write(f"{len(points)}\n")
        for i, xyz in enumerate(points, start=1):
            f.write(f"{i} {xyz[0]:.17e} {xyz[1]:.17e} {xyz[2]:.17e}\n")

    with (args.outdir / "cells.dat").open("w", encoding="utf-8") as f:
        f.write(f"{len(cells)}\n")
        for cell in cells:
            nodes = " ".join(str(n) for n in cell["nodes"])
            c = cell["center"]
            f.write(f"{cell['id']} {nodes} {c[0]:.17e} {c[1]:.17e} {c[2]:.17e} {cell['volume']:.17e}\n")

    with (args.outdir / "faces.dat").open("w", encoding="utf-8") as f:
        f.write(f"{len(faces)}\n")
        for face in faces:
            n = face["normal"]
            c = face["center"]
            f.write(
                f"{face['id']} {face['owner']} {face['neighbor']} {face['patch']} "
                f"{n[0]:.17e} {n[1]:.17e} {n[2]:.17e} {face['area']:.17e} "
                f"{c[0]:.17e} {c[1]:.17e} {c[2]:.17e}\n"
            )

    id_to_patch = {pid: name for name, pid in patch_name_to_id.items()}
    with (args.outdir / "patches.dat").open("w", encoding="utf-8") as f:
        f.write(f"{len(id_to_patch)}\n")
        for pid in sorted(id_to_patch):
            face_ids = patch_faces.get(pid, [])
            f.write(f"{pid} {id_to_patch[pid]} {len(face_ids)}\n")
            f.write(" ".join(str(fid) for fid in face_ids) + "\n")

    if periodic_links:
        with (args.outdir / "periodic.dat").open("w", encoding="utf-8") as f:
            f.write(f"{len(periodic_links)}\n")
            for face_id, pair_face_id, neighbor_cell in periodic_links:
                f.write(f"{face_id} {pair_face_id} {neighbor_cell}\n")

    print(f"Wrote native mesh to {args.outdir}")
    print(f"  points:  {len(points)}")
    print(f"  cells:   {len(cells)}")
    print(f"  faces:   {len(faces)}")
    print(f"  patches: {len(id_to_patch)}")
    print(f"  periodic links: {len(periodic_links)}")


if __name__ == "__main__":
    main()
