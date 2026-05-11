#!/usr/bin/env python3
"""Generate native axis-aligned cuboid block meshes for the low-Mach FV solver.

This is a Level-1 mesh generator: it creates conformal structured/block-structured
cuboid meshes directly in the solver's native format:

  points.dat
  cells.dat
  faces.dat
  patches.dat
  periodic.dat  optional

It does not require Gmsh. The user provides a small Python config file with:
  X, Y, Z               coordinate cut planes
  NX, NY, NZ            number of cells in each interval
  SOLID_BLOCKS          omitted block indices, e.g. [(2, 2, 0)]
  PATCHES               side/solid-to-patch-name map
  PERIODIC              optional patch pairs, e.g. [("zmin", "zmax")]

Only axis-aligned cuboid cells are generated.
"""

from __future__ import annotations

import argparse
import math
import runpy
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

TOL = 1.0e-10


Vec3 = Tuple[float, float, float]
NodeTuple8 = Tuple[int, int, int, int, int, int, int, int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="Python mesh config file")
    parser.add_argument("outdir", type=Path, help="Output native mesh directory")
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove output directory before writing",
    )
    return parser.parse_args()


def require(cfg: Dict[str, object], name: str):
    if name not in cfg:
        raise ValueError(f"Config is missing required variable {name}")
    return cfg[name]


def as_float_list(values, name: str) -> List[float]:
    try:
        out = [float(v) for v in values]
    except Exception as exc:
        raise ValueError(f"{name} must be a list of numbers") from exc
    if len(out) < 2:
        raise ValueError(f"{name} must contain at least two coordinates")
    for a, b in zip(out[:-1], out[1:]):
        if not b > a:
            raise ValueError(f"{name} must be strictly increasing: found {a} then {b}")
    return out


def as_int_list(values, name: str) -> List[int]:
    try:
        out = [int(v) for v in values]
    except Exception as exc:
        raise ValueError(f"{name} must be a list of positive integers") from exc
    if any(v <= 0 for v in out):
        raise ValueError(f"{name} must contain only positive integers")
    return out


def expand_axis(coords: List[float], counts: List[int]) -> Tuple[List[float], List[int]]:
    """Return full node coordinates and block index for each cell."""
    if len(coords) != len(counts) + 1:
        raise ValueError(
            "Axis mismatch: len(coords) must equal len(cells) + 1. "
            f"Got {len(coords)} coords and {len(counts)} cell counts."
        )

    full = [coords[0]]
    block_of_cell: List[int] = []

    for ib, n in enumerate(counts):
        a = coords[ib]
        b = coords[ib + 1]
        for m in range(n):
            block_of_cell.append(ib)
            full.append(a + (b - a) * float(m + 1) / float(n))

    return full, block_of_cell


def close(a: float, b: float) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= TOL * scale


def round_key(value: float) -> int:
    return int(round(value / TOL))


class MeshBuilder:
    def __init__(self, cfg: Dict[str, object]):
        self.Xc = as_float_list(require(cfg, "X"), "X")
        self.Yc = as_float_list(require(cfg, "Y"), "Y")
        self.Zc = as_float_list(require(cfg, "Z"), "Z")

        self.NX = as_int_list(require(cfg, "NX"), "NX")
        self.NY = as_int_list(require(cfg, "NY"), "NY")
        self.NZ = as_int_list(require(cfg, "NZ"), "NZ")

        self.x, self.x_block = expand_axis(self.Xc, self.NX)
        self.y, self.y_block = expand_axis(self.Yc, self.NY)
        self.z, self.z_block = expand_axis(self.Zc, self.NZ)

        self.nx = len(self.x) - 1
        self.ny = len(self.y) - 1
        self.nz = len(self.z) - 1

        solids = cfg.get("SOLID_BLOCKS", [])
        self.solid_blocks = {tuple(map(int, item)) for item in solids}

        self.patches = dict(cfg.get("PATCHES", {}))
        self.periodic = [tuple(pair) for pair in cfg.get("PERIODIC", [])]

        self.default_patches = {
            "xmin": "xmin",
            "xmax": "xmax",
            "ymin": "ymin",
            "ymax": "ymax",
            "zmin": "zmin",
            "zmax": "zmax",
            "solid": "wall",
        }
        for key, value in self.default_patches.items():
            self.patches.setdefault(key, value)

        self.points: List[Vec3] = []
        self.cells: List[dict] = []
        self.faces: List[dict] = []
        self.face_by_nodes: Dict[Tuple[int, ...], dict] = {}
        self.patch_name_to_id: Dict[str, int] = {}
        self.patch_faces: Dict[int, List[int]] = defaultdict(list)
        self.cell_id_by_ijk: Dict[Tuple[int, int, int], int] = {}

    def node_id(self, i: int, j: int, k: int) -> int:
        """1-based native node id."""
        return 1 + i + (self.nx + 1) * (j + (self.ny + 1) * k)

    def cell_nodes(self, i: int, j: int, k: int) -> NodeTuple8:
        n000 = self.node_id(i,     j,     k)
        n100 = self.node_id(i + 1, j,     k)
        n110 = self.node_id(i + 1, j + 1, k)
        n010 = self.node_id(i,     j + 1, k)
        n001 = self.node_id(i,     j,     k + 1)
        n101 = self.node_id(i + 1, j,     k + 1)
        n111 = self.node_id(i + 1, j + 1, k + 1)
        n011 = self.node_id(i,     j + 1, k + 1)
        return (n000, n100, n110, n010, n001, n101, n111, n011)

    def cell_block(self, i: int, j: int, k: int) -> Tuple[int, int, int]:
        return (self.x_block[i], self.y_block[j], self.z_block[k])

    def is_solid_cell(self, i: int, j: int, k: int) -> bool:
        return self.cell_block(i, j, k) in self.solid_blocks

    def ensure_patch_id(self, name: str) -> int:
        if name not in self.patch_name_to_id:
            self.patch_name_to_id[name] = len(self.patch_name_to_id) + 1
        return self.patch_name_to_id[name]

    def build_points(self) -> None:
        for k in range(self.nz + 1):
            for j in range(self.ny + 1):
                for i in range(self.nx + 1):
                    self.points.append((self.x[i], self.y[j], self.z[k]))

    def add_cell(self, i: int, j: int, k: int) -> None:
        x0, x1 = self.x[i], self.x[i + 1]
        y0, y1 = self.y[j], self.y[j + 1]
        z0, z1 = self.z[k], self.z[k + 1]

        dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
        if dx <= 0.0 or dy <= 0.0 or dz <= 0.0:
            raise ValueError("Degenerate cell detected")

        cell_id = len(self.cells) + 1
        nodes = self.cell_nodes(i, j, k)
        center = (0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * (z0 + z1))
        volume = dx * dy * dz

        self.cells.append(
            {
                "id": cell_id,
                "ijk": (i, j, k),
                "block": self.cell_block(i, j, k),
                "nodes": nodes,
                "center": center,
                "volume": volume,
            }
        )
        self.cell_id_by_ijk[(i, j, k)] = cell_id

    def face_definitions_for_cell(self, i: int, j: int, k: int) -> List[dict]:
        x0, x1 = self.x[i], self.x[i + 1]
        y0, y1 = self.y[j], self.y[j + 1]
        z0, z1 = self.z[k], self.z[k + 1]

        dx, dy, dz = x1 - x0, y1 - y0, z1 - z0

        n000, n100, n110, n010, n001, n101, n111, n011 = self.cell_nodes(i, j, k)

        # Node order below is only for defining unique face sets. The native file
        # stores normals/areas/centers explicitly, and the solver flips normals
        # using owner/neighbor when needed.
        return [
            {
                "side": "xmin",
                "nodes": (n000, n010, n011, n001),
                "normal": (-1.0, 0.0, 0.0),
                "area": dy * dz,
                "center": (x0, 0.5 * (y0 + y1), 0.5 * (z0 + z1)),
            },
            {
                "side": "xmax",
                "nodes": (n100, n101, n111, n110),
                "normal": (1.0, 0.0, 0.0),
                "area": dy * dz,
                "center": (x1, 0.5 * (y0 + y1), 0.5 * (z0 + z1)),
            },
            {
                "side": "ymin",
                "nodes": (n000, n001, n101, n100),
                "normal": (0.0, -1.0, 0.0),
                "area": dx * dz,
                "center": (0.5 * (x0 + x1), y0, 0.5 * (z0 + z1)),
            },
            {
                "side": "ymax",
                "nodes": (n010, n110, n111, n011),
                "normal": (0.0, 1.0, 0.0),
                "area": dx * dz,
                "center": (0.5 * (x0 + x1), y1, 0.5 * (z0 + z1)),
            },
            {
                "side": "zmin",
                "nodes": (n000, n100, n110, n010),
                "normal": (0.0, 0.0, -1.0),
                "area": dx * dy,
                "center": (0.5 * (x0 + x1), 0.5 * (y0 + y1), z0),
            },
            {
                "side": "zmax",
                "nodes": (n001, n011, n111, n101),
                "normal": (0.0, 0.0, 1.0),
                "area": dx * dy,
                "center": (0.5 * (x0 + x1), 0.5 * (y0 + y1), z1),
            },
        ]

    def add_faces_for_cell(self, cell: dict) -> None:
        i, j, k = cell["ijk"]
        for spec in self.face_definitions_for_cell(i, j, k):
            key = tuple(sorted(spec["nodes"]))
            if key in self.face_by_nodes:
                face = self.face_by_nodes[key]
                if face["neighbor"] != 0:
                    raise ValueError("Face shared by more than two cells")
                face["neighbor"] = cell["id"]
            else:
                face = {
                    "id": len(self.faces) + 1,
                    "owner": cell["id"],
                    "neighbor": 0,
                    "patch": 0,
                    "normal": spec["normal"],
                    "area": spec["area"],
                    "center": spec["center"],
                    "side": spec["side"],
                    "nodes_key": key,
                    "periodic_face": 0,
                    "periodic_neighbor": 0,
                }
                self.faces.append(face)
                self.face_by_nodes[key] = face

    def patch_name_for_boundary_face(self, face: dict) -> str:
        x, y, z = face["center"]

        if close(x, self.x[0]) and face["normal"][0] < 0.0:
            return self.patches["xmin"]
        if close(x, self.x[-1]) and face["normal"][0] > 0.0:
            return self.patches["xmax"]
        if close(y, self.y[0]) and face["normal"][1] < 0.0:
            return self.patches["ymin"]
        if close(y, self.y[-1]) and face["normal"][1] > 0.0:
            return self.patches["ymax"]
        if close(z, self.z[0]) and face["normal"][2] < 0.0:
            return self.patches["zmin"]
        if close(z, self.z[-1]) and face["normal"][2] > 0.0:
            return self.patches["zmax"]

        # If it is not on the outer box and has no neighbor, it is a solid-block
        # boundary face.
        return self.patches["solid"]

    def assign_boundary_patches(self) -> None:
        for face in self.faces:
            if face["neighbor"] != 0:
                face["patch"] = 0
                continue

            patch_name = self.patch_name_for_boundary_face(face)
            patch_id = self.ensure_patch_id(patch_name)
            face["patch"] = patch_id
            self.patch_faces[patch_id].append(face["id"])

    def periodic_key(self, face: dict) -> Tuple[int, int, int]:
        normal = face["normal"]
        axis = max(range(3), key=lambda a: abs(normal[a]))
        other = [a for a in range(3) if a != axis]
        center = face["center"]
        # Include axis so x-periodic faces cannot accidentally pair with z faces.
        return (axis, round_key(center[other[0]]), round_key(center[other[1]]))

    def pair_periodic_faces(self) -> List[Tuple[int, int, int]]:
        links: List[Tuple[int, int, int]] = []
        id_to_patch = {pid: name for name, pid in self.patch_name_to_id.items()}
        patch_to_faces: Dict[str, List[dict]] = defaultdict(list)

        for face in self.faces:
            if face["neighbor"] == 0:
                patch_to_faces[id_to_patch[face["patch"]]].append(face)

        for a_name, b_name in self.periodic:
            if a_name not in self.patch_name_to_id:
                raise ValueError(f"Periodic patch {a_name!r} was not created")
            if b_name not in self.patch_name_to_id:
                raise ValueError(f"Periodic patch {b_name!r} was not created")

            a_faces = patch_to_faces[a_name]
            b_faces = patch_to_faces[b_name]

            if len(a_faces) != len(b_faces):
                raise ValueError(
                    f"Periodic patches {a_name}:{b_name} have different face counts: "
                    f"{len(a_faces)} vs {len(b_faces)}"
                )

            b_by_key = {}
            for fb in b_faces:
                key = self.periodic_key(fb)
                if key in b_by_key:
                    raise ValueError(f"Duplicate periodic matching key on patch {b_name}")
                b_by_key[key] = fb

            for fa in a_faces:
                key = self.periodic_key(fa)
                fb = b_by_key.get(key)
                if fb is None:
                    raise ValueError(f"No periodic mate for face {fa['id']} in pair {a_name}:{b_name}")

                # Sanity: periodic faces should have opposite normals.
                dot = sum(fa["normal"][m] * fb["normal"][m] for m in range(3))
                if dot > -0.5:
                    raise ValueError(
                        f"Periodic pair {fa['id']}:{fb['id']} does not have opposite normals"
                    )

                fa["periodic_face"] = fb["id"]
                fa["periodic_neighbor"] = fb["owner"]
                fb["periodic_face"] = fa["id"]
                fb["periodic_neighbor"] = fa["owner"]

                links.append((fa["id"], fb["id"], fb["owner"]))
                links.append((fb["id"], fa["id"], fa["owner"]))

        return links

    def check_normals(self) -> None:
        cell_centers = {cell["id"]: cell["center"] for cell in self.cells}

        for face in self.faces:
            owner_center = cell_centers[face["owner"]]
            n = face["normal"]

            if face["neighbor"] > 0:
                nb_center = cell_centers[face["neighbor"]]
                delta = tuple(nb_center[a] - owner_center[a] for a in range(3))
                dot = sum(delta[a] * n[a] for a in range(3))
                if dot <= 0.0:
                    raise ValueError(
                        f"Face {face['id']} normal does not point from owner "
                        f"{face['owner']} to neighbor {face['neighbor']}"
                    )
            else:
                delta = tuple(face["center"][a] - owner_center[a] for a in range(3))
                dot = sum(delta[a] * n[a] for a in range(3))
                if dot <= 0.0:
                    raise ValueError(
                        f"Boundary face {face['id']} normal does not point outward from owner"
                    )

    def build(self) -> List[Tuple[int, int, int]]:
        self.build_points()

        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    if not self.is_solid_cell(i, j, k):
                        self.add_cell(i, j, k)

        if not self.cells:
            raise ValueError("No fluid cells were generated")

        for cell in self.cells:
            self.add_faces_for_cell(cell)

        self.assign_boundary_patches()
        links = self.pair_periodic_faces()
        self.check_normals()
        return links

    def write(self, outdir: Path, links: List[Tuple[int, int, int]], clean: bool = False) -> None:
        if clean and outdir.exists():
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        with (outdir / "points.dat").open("w", encoding="utf-8") as f:
            f.write(f"{len(self.points)}\n")
            for i, xyz in enumerate(self.points, start=1):
                f.write(f"{i} {xyz[0]:.17e} {xyz[1]:.17e} {xyz[2]:.17e}\n")

        with (outdir / "cells.dat").open("w", encoding="utf-8") as f:
            f.write(f"{len(self.cells)}\n")
            for cell in self.cells:
                nodes = " ".join(str(n) for n in cell["nodes"])
                c = cell["center"]
                f.write(
                    f"{cell['id']} {nodes} "
                    f"{c[0]:.17e} {c[1]:.17e} {c[2]:.17e} "
                    f"{cell['volume']:.17e}\n"
                )

        with (outdir / "faces.dat").open("w", encoding="utf-8") as f:
            f.write(f"{len(self.faces)}\n")
            for face in self.faces:
                n = face["normal"]
                c = face["center"]
                f.write(
                    f"{face['id']} {face['owner']} {face['neighbor']} {face['patch']} "
                    f"{n[0]:.17e} {n[1]:.17e} {n[2]:.17e} "
                    f"{face['area']:.17e} "
                    f"{c[0]:.17e} {c[1]:.17e} {c[2]:.17e}\n"
                )

        id_to_patch = {pid: name for name, pid in self.patch_name_to_id.items()}
        with (outdir / "patches.dat").open("w", encoding="utf-8") as f:
            f.write(f"{len(id_to_patch)}\n")
            for pid in sorted(id_to_patch):
                face_ids = self.patch_faces.get(pid, [])
                f.write(f"{pid} {id_to_patch[pid]} {len(face_ids)}\n")
                f.write(" ".join(str(fid) for fid in face_ids) + "\n")

        if links:
            with (outdir / "periodic.dat").open("w", encoding="utf-8") as f:
                f.write(f"{len(links)}\n")
                for face_id, pair_face_id, neighbor_cell in links:
                    f.write(f"{face_id} {pair_face_id} {neighbor_cell}\n")

        self.print_summary(outdir, links)

    def print_summary(self, outdir: Path, links: List[Tuple[int, int, int]]) -> None:
        id_to_patch = {pid: name for name, pid in self.patch_name_to_id.items()}

        print(f"Wrote native mesh to {outdir}")
        print(f"  points:         {len(self.points)}")
        print(f"  cells:          {len(self.cells)}")
        print(f"  faces:          {len(self.faces)}")
        print(f"  patches:        {len(id_to_patch)}")
        print(f"  periodic links: {len(links)}")
        print("  patch face counts:")
        for pid in sorted(id_to_patch):
            print(f"    {id_to_patch[pid]}: {len(self.patch_faces.get(pid, []))}")


def main() -> None:
    args = parse_args()
    cfg = runpy.run_path(str(args.config))

    builder = MeshBuilder(cfg)
    links = builder.build()
    builder.write(args.outdir, links, clean=args.clean)


if __name__ == "__main__":
    main()
