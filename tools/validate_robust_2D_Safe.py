#!/usr/bin/env python3
"""
Validation utilities for hex_lowmach_fv VTU output.

Supports:
  1. Periodic body-force channel validation against analytic Poiseuille flow.
  2. Lid-driven cavity centerline extraction and comparison to Ghia Re=100 data.

This version is safe for:
  - 3D hexahedral VTU output, including thin Nx x Ny x 1 meshes.
  - True 2D quad/triangle VTU output.
  - Cell-centered velocity with 2 or 3 components.
  - Point-centered velocity with 2 or 3 components, averaged to cells.

Requirements:
  pip install numpy pandas matplotlib meshio
"""

from __future__ import annotations

import argparse
from pathlib import Path

import meshio
import numpy as np
import pandas as pd


GHIA_RE100_U_VERTICAL = np.array(
    [
        [1.0000, 1.00000],
        [0.9766, 0.84123],
        [0.9688, 0.78871],
        [0.9609, 0.73722],
        [0.9531, 0.68717],
        [0.8516, 0.23151],
        [0.7344, 0.00332],
        [0.6172, -0.13641],
        [0.5000, -0.20581],
        [0.4531, -0.21090],
        [0.2813, -0.15662],
        [0.1719, -0.10150],
        [0.1016, -0.06434],
        [0.0703, -0.04775],
        [0.0625, -0.04192],
        [0.0547, -0.03717],
        [0.0000, 0.00000],
    ],
    dtype=float,
)


GHIA_RE100_V_HORIZONTAL = np.array(
    [
        [1.0000, 0.00000],
        [0.9688, -0.05906],
        [0.9609, -0.07391],
        [0.9531, -0.08864],
        [0.9453, -0.10313],
        [0.9063, -0.16914],
        [0.8594, -0.22445],
        [0.8047, -0.24533],
        [0.5000, 0.05454],
        [0.2344, 0.17527],
        [0.2266, 0.17507],
        [0.1563, 0.16077],
        [0.0938, 0.12317],
        [0.0781, 0.10890],
        [0.0703, 0.10091],
        [0.0625, 0.09233],
        [0.0000, 0.00000],
    ],
    dtype=float,
)


def default_validation_dir(vtu_path: str | Path) -> Path:
    vtu_path = Path(vtu_path)
    if vtu_path.parent.name == "output":
        return vtu_path.parent.parent / "validation"
    return vtu_path.parent / "validation"


def points_to_3d(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)

    if points.ndim != 2:
        raise RuntimeError(f"Expected points to be 2D array, got shape {points.shape}")

    if points.shape[1] == 3:
        return points

    if points.shape[1] == 2:
        out = np.zeros((points.shape[0], 3), dtype=points.dtype)
        out[:, 0:2] = points
        return out

    raise RuntimeError(f"Expected points with 2 or 3 coordinates, got shape {points.shape}")


def vector_to_3d(vector: np.ndarray, name: str = "vector") -> np.ndarray:
    """
    Convert [u, v] to [u, v, 0] for true 2D output.
    """
    vector = np.asarray(vector, dtype=float)

    if vector.ndim == 1:
        raise RuntimeError(f"Expected {name} to have 2 or 3 components, got shape {vector.shape}")

    if vector.ndim != 2:
        raise RuntimeError(f"Expected {name} to be 2D array, got shape {vector.shape}")

    if vector.shape[1] == 3:
        return vector

    if vector.shape[1] == 2:
        out = np.zeros((vector.shape[0], 3), dtype=vector.dtype)
        out[:, 0:2] = vector
        return out

    raise RuntimeError(f"Expected {name} with 2 or 3 components, got shape {vector.shape}")


def find_case_insensitive_key(data: dict, target: str) -> str | None:
    target_lower = target.lower()
    for key in data.keys():
        if str(key).lower() == target_lower:
            return key
    return None


def choose_cell_block(mesh: meshio.Mesh) -> int:
    if not mesh.cells:
        raise RuntimeError("No cells found in VTU file")

    preferred_types = (
        "hexahedron",
        "hex8",
        "quad",
        "quad4",
        "triangle",
        "triangle3",
        "tetra",
        "wedge",
        "pyramid",
    )

    for cell_type in preferred_types:
        for i, block in enumerate(mesh.cells):
            if block.type == cell_type:
                return i

    return 0


def read_vtu_cell_data(vtu_path: str | Path):
    mesh = meshio.read(vtu_path)

    block_index = choose_cell_block(mesh)
    block = mesh.cells[block_index]
    cells = np.asarray(block.data)

    points = points_to_3d(mesh.points)
    centers = points[cells].mean(axis=1)

    point_mins = points.min(axis=0)
    point_maxs = points.max(axis=0)

    cell_data = {}
    for name, blocks in mesh.cell_data.items():
        if len(blocks) > block_index:
            cell_data[name] = np.asarray(blocks[block_index])

    velocity_key = find_case_insensitive_key(cell_data, "velocity")

    if velocity_key is not None:
        velocity = vector_to_3d(cell_data[velocity_key], name=velocity_key)
        return centers, velocity, cell_data, point_mins, point_maxs, block.type

    point_velocity_key = find_case_insensitive_key(mesh.point_data, "velocity")

    if point_velocity_key is not None:
        point_velocity = vector_to_3d(mesh.point_data[point_velocity_key], name=point_velocity_key)
        velocity = point_velocity[cells].mean(axis=1)
        return centers, velocity, cell_data, point_mins, point_maxs, block.type

    raise RuntimeError(
        "VTU has no CellData or PointData array named 'velocity'. "
        f"Available cell arrays: {list(cell_data.keys())}. "
        f"Available point arrays: {list(mesh.point_data.keys())}."
    )


def estimate_centerline_tolerance(centers: np.ndarray) -> float:
    spacings = []

    for axis in range(3):
        vals = np.unique(np.round(centers[:, axis], 12))
        vals = np.sort(vals)
        diffs = np.diff(vals)
        diffs = diffs[diffs > 1.0e-14]

        if diffs.size > 0:
            spacings.append(float(np.median(diffs)))

    if not spacings:
        return 0.05

    return 0.60 * max(spacings)


def domain_midpoints_from_bounds(point_mins: np.ndarray, point_maxs: np.ndarray):
    mids = 0.5 * (point_mins + point_maxs)
    return point_mins, point_maxs, mids


def has_multiple_layers(centers: np.ndarray, axis: int, eps: float = 1.0e-12) -> bool:
    vals = np.unique(np.round(centers[:, axis], 12))
    return vals.size > 1 and float(vals.max() - vals.min()) > eps


def bin_average(coord: np.ndarray, value: np.ndarray, nbins: int, cmin: float, cmax: float) -> pd.DataFrame:
    edges = np.linspace(cmin, cmax, nbins + 1)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])

    rows = []

    for i in range(nbins):
        lo = edges[i]
        hi = edges[i + 1]

        if i == nbins - 1:
            mask = (coord >= lo) & (coord <= hi)
        else:
            mask = (coord >= lo) & (coord < hi)

        if np.any(mask):
            rows.append(
                {
                    "bin": i,
                    "coord": bin_centers[i],
                    "count": int(mask.sum()),
                    "mean": float(np.mean(value[mask])),
                    "std": float(np.std(value[mask])),
                }
            )
        else:
            rows.append(
                {
                    "bin": i,
                    "coord": bin_centers[i],
                    "count": 0,
                    "mean": np.nan,
                    "std": np.nan,
                }
            )

    return pd.DataFrame(rows)


def validate_channel(args: argparse.Namespace) -> None:
    centers, velocity, _, point_mins, point_maxs, cell_type = read_vtu_cell_data(args.vtu)

    outdir = Path(args.outdir) if args.outdir else default_validation_dir(args.vtu)
    outdir.mkdir(parents=True, exist_ok=True)

    y = centers[:, 1]
    u = velocity[:, 0]
    v = velocity[:, 1]
    w = velocity[:, 2]

    profile = bin_average(y, u, args.nbins, args.ymin, args.ymax)

    H = args.ymax - args.ymin
    yrel = profile["coord"].to_numpy() - args.ymin

    profile["u_numeric"] = profile["mean"]
    profile["u_exact"] = args.body_force_x / (2.0 * args.nu) * yrel * (H - yrel)
    profile["error"] = profile["u_numeric"] - profile["u_exact"]
    profile["abs_error"] = np.abs(profile["error"])

    valid = np.isfinite(profile["u_numeric"].to_numpy())
    err = profile.loc[valid, "error"].to_numpy()
    exact = profile.loc[valid, "u_exact"].to_numpy()

    l2 = float(np.sqrt(np.mean(err**2)))
    linf = float(np.max(np.abs(err)))
    exact_l2 = float(np.sqrt(np.mean(exact**2)))
    rel_l2 = l2 / max(exact_l2, 1.0e-300)

    u_max_exact = args.body_force_x * H**2 / (8.0 * args.nu)
    u_mean_exact = args.body_force_x * H**2 / (12.0 * args.nu)

    summary = {
        "vtu": str(args.vtu),
        "cell_type_used": cell_type,
        "domain_x_min": float(point_mins[0]),
        "domain_x_max": float(point_maxs[0]),
        "domain_y_min": float(point_mins[1]),
        "domain_y_max": float(point_maxs[1]),
        "domain_z_min": float(point_mins[2]),
        "domain_z_max": float(point_maxs[2]),
        "nu": args.nu,
        "body_force_x": args.body_force_x,
        "H": H,
        "u_max_numeric": float(np.nanmax(profile["u_numeric"])),
        "u_max_exact": float(u_max_exact),
        "u_mean_numeric": float(np.mean(u)),
        "u_mean_exact": float(u_mean_exact),
        "L2_error": l2,
        "Linf_error": linf,
        "relative_L2_error": float(rel_l2),
        "max_abs_v": float(np.max(np.abs(v))),
        "max_abs_w": float(np.max(np.abs(w))),
    }

    profile_path = outdir / "channel_profile.csv"
    summary_path = outdir / "channel_summary.txt"
    plot_path = outdir / "channel_profile.png"

    profile.to_csv(profile_path, index=False)

    with summary_path.open("w") as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")

    print("\nChannel validation summary")
    print("--------------------------")
    for key, value in summary.items():
        print(f"{key}: {value}")

    print(f"\nWrote: {profile_path}")
    print(f"Wrote: {summary_path}")

    if args.plot:
        plot_channel(profile, plot_path)


def plot_channel(profile: pd.DataFrame, plot_path: Path) -> None:
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(profile["u_numeric"], profile["coord"], "o", label="numerical")
    plt.plot(profile["u_exact"], profile["coord"], "-", label="analytic")
    plt.xlabel("u")
    plt.ylabel("y")
    plt.title("Channel flow: numerical vs analytic")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200)
    plt.close()

    print(f"Wrote: {plot_path}")


def nearest_line_profile(
    centers: np.ndarray,
    values: np.ndarray,
    fixed: dict[int, float],
    varying_axis: int,
    tol: float,
) -> pd.DataFrame:
    mask = np.ones(centers.shape[0], dtype=bool)

    for axis, target in fixed.items():
        mask &= np.abs(centers[:, axis] - target) <= tol

    if not np.any(mask):
        mask = np.ones(centers.shape[0], dtype=bool)
        nearest_used = {}

        for axis, target in fixed.items():
            vals = np.unique(np.round(centers[:, axis], 12))
            nearest = float(vals[np.argmin(np.abs(vals - target))])
            nearest_used[axis] = nearest
            mask &= np.abs(centers[:, axis] - nearest) <= 1.0e-10

        if not np.any(mask):
            raise RuntimeError(
                f"No cells found within tol={tol} of fixed coordinates {fixed}. "
                f"Nearest layers tried: {nearest_used}. "
                "Increase --tol or check domain coordinates."
            )

        print(
            "Warning: no cells found with requested centerline tolerance. "
            f"Using nearest available cell-center layers instead: {nearest_used}"
        )

    coord = centers[mask, varying_axis]
    val = values[mask]

    order = np.argsort(coord)
    df = pd.DataFrame({"coord": coord[order], "value": val[order]})

    df["coord_round"] = df["coord"].round(12)

    df = (
        df.groupby("coord_round", as_index=False)
        .agg(coord=("coord", "mean"), value=("value", "mean"), count=("value", "count"))
        .sort_values("coord")
    )

    return df[["coord", "value", "count"]]


def add_boundary_point(
    profile: pd.DataFrame,
    coord_name: str,
    value_name: str,
    coord: float,
    value: float,
) -> pd.DataFrame:
    row = pd.DataFrame([{coord_name: coord, value_name: value, "count": 0}])
    out = pd.concat([profile, row], ignore_index=True)
    out = out.sort_values(coord_name)
    out = out.drop_duplicates(subset=[coord_name], keep="last")
    return out.reset_index(drop=True)


def add_cavity_boundary_values(
    u_vertical: pd.DataFrame,
    v_horizontal: pd.DataFrame,
    point_mins: np.ndarray,
    point_maxs: np.ndarray,
    lid_u: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Add known wall values to extracted cavity profiles.

    Cell-centered FV output usually has no samples exactly at y=0, y=1,
    x=0, or x=1, while Ghia includes wall values.
    """
    u_vertical = add_boundary_point(u_vertical, "y", "u", float(point_mins[1]), 0.0)
    u_vertical = add_boundary_point(u_vertical, "y", "u", float(point_maxs[1]), float(lid_u))

    v_horizontal = add_boundary_point(v_horizontal, "x", "v", float(point_mins[0]), 0.0)
    v_horizontal = add_boundary_point(v_horizontal, "x", "v", float(point_maxs[0]), 0.0)

    return u_vertical, v_horizontal


def compare_to_benchmark(
    profile: pd.DataFrame,
    benchmark: np.ndarray,
    coord_col: str,
    value_col: str,
) -> pd.DataFrame:
    prof = profile.sort_values(coord_col)
    coord = prof[coord_col].to_numpy()
    value = prof[value_col].to_numpy()

    bcoord = benchmark[:, 0]
    bval = benchmark[:, 1]

    numerical_at_benchmark = np.interp(bcoord, coord, value)

    return pd.DataFrame(
        {
            coord_col: bcoord,
            f"{value_col}_benchmark": bval,
            f"{value_col}_numeric_interp": numerical_at_benchmark,
            "error": numerical_at_benchmark - bval,
            "abs_error": np.abs(numerical_at_benchmark - bval),
        }
    )


def summarize_error(df: pd.DataFrame) -> dict[str, float]:
    err = df["error"].to_numpy()
    ref_cols = [c for c in df.columns if c.endswith("_benchmark")]
    ref = df[ref_cols[0]].to_numpy()

    l2 = float(np.sqrt(np.mean(err**2)))
    linf = float(np.max(np.abs(err)))
    ref_l2 = float(np.sqrt(np.mean(ref**2)))
    rel_l2 = l2 / max(ref_l2, 1.0e-300)

    return {"L2": l2, "Linf": linf, "relative_L2": rel_l2}


def validate_cavity(args: argparse.Namespace) -> None:
    centers, velocity, _, point_mins, point_maxs, cell_type = read_vtu_cell_data(args.vtu)

    outdir = Path(args.outdir) if args.outdir else default_validation_dir(args.vtu)
    outdir.mkdir(parents=True, exist_ok=True)

    mins, maxs, mids = domain_midpoints_from_bounds(point_mins, point_maxs)

    xmid = args.xmid if args.xmid is not None else float(mids[0])
    ymid = args.ymid if args.ymid is not None else float(mids[1])
    zmid = args.zmid if args.zmid is not None else float(mids[2])

    tol = args.tol if args.tol is not None else estimate_centerline_tolerance(centers)

    u = velocity[:, 0]
    v = velocity[:, 1]
    w = velocity[:, 2]

    use_z_filter = has_multiple_layers(centers, axis=2) and not args.no_z_filter

    fixed_u = {0: xmid}
    fixed_v = {1: ymid}

    if use_z_filter:
        fixed_u[2] = zmid
        fixed_v[2] = zmid

    u_vertical = nearest_line_profile(
        centers=centers,
        values=u,
        fixed=fixed_u,
        varying_axis=1,
        tol=tol,
    ).rename(columns={"coord": "y", "value": "u"})

    v_horizontal = nearest_line_profile(
        centers=centers,
        values=v,
        fixed=fixed_v,
        varying_axis=0,
        tol=tol,
    ).rename(columns={"coord": "x", "value": "v"})

    if not args.no_wall_points:
        u_vertical, v_horizontal = add_cavity_boundary_values(
            u_vertical=u_vertical,
            v_horizontal=v_horizontal,
            point_mins=point_mins,
            point_maxs=point_maxs,
            lid_u=args.lid_u,
        )

    u_path = outdir / "cavity_u_vertical_centerline.csv"
    v_path = outdir / "cavity_v_horizontal_centerline.csv"

    u_vertical.to_csv(u_path, index=False)
    v_horizontal.to_csv(v_path, index=False)

    summary = {
        "vtu": str(args.vtu),
        "cell_type_used": cell_type,
        "domain_x_min": float(mins[0]),
        "domain_x_max": float(maxs[0]),
        "domain_y_min": float(mins[1]),
        "domain_y_max": float(maxs[1]),
        "domain_z_min": float(mins[2]),
        "domain_z_max": float(maxs[2]),
        "xmid_used": xmid,
        "ymid_used": ymid,
        "zmid_used": zmid,
        "tol_used": tol,
        "z_filter_used": bool(use_z_filter),
        "wall_points_added": bool(not args.no_wall_points),
        "u_vertical_points": len(u_vertical),
        "v_horizontal_points": len(v_horizontal),
        "max_abs_w": float(np.max(np.abs(w))),
        "max_speed": float(np.max(np.linalg.norm(velocity, axis=1))),
    }

    if args.re == 100:
        u_compare = compare_to_benchmark(u_vertical, GHIA_RE100_U_VERTICAL, "y", "u")
        v_compare = compare_to_benchmark(v_horizontal, GHIA_RE100_V_HORIZONTAL, "x", "v")

        u_compare_path = outdir / "cavity_Re100_u_vs_ghia.csv"
        v_compare_path = outdir / "cavity_Re100_v_vs_ghia.csv"

        u_compare.to_csv(u_compare_path, index=False)
        v_compare.to_csv(v_compare_path, index=False)

        u_error = summarize_error(u_compare)
        v_error = summarize_error(v_compare)

        summary.update(
            {
                "benchmark": "Ghia et al. Re=100 2D cavity centerline data",
                "u_centerline_L2": u_error["L2"],
                "u_centerline_Linf": u_error["Linf"],
                "u_centerline_relative_L2": u_error["relative_L2"],
                "v_centerline_L2": v_error["L2"],
                "v_centerline_Linf": v_error["Linf"],
                "v_centerline_relative_L2": v_error["relative_L2"],
            }
        )

        print(f"Wrote: {u_compare_path}")
        print(f"Wrote: {v_compare_path}")
    else:
        print(
            "\nNo built-in cavity benchmark for this Re yet. "
            "This script currently includes Ghia Re=100 only."
        )

    summary_path = outdir / "cavity_summary.txt"

    with summary_path.open("w") as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")

    print("\nCavity validation summary")
    print("-------------------------")
    for key, value in summary.items():
        print(f"{key}: {value}")

    print(f"\nWrote: {u_path}")
    print(f"Wrote: {v_path}")
    print(f"Wrote: {summary_path}")

    if args.plot:
        plot_cavity(outdir, u_vertical, v_horizontal, args.re)


def plot_cavity(
    outdir: Path,
    u_vertical: pd.DataFrame,
    v_horizontal: pd.DataFrame,
    re: int,
) -> None:
    import matplotlib.pyplot as plt

    p1 = outdir / "cavity_u_vertical_centerline.png"
    p2 = outdir / "cavity_v_horizontal_centerline.png"

    plt.figure()
    plt.plot(u_vertical["u"], u_vertical["y"], "o-", label="numerical")

    if re == 100:
        plt.plot(GHIA_RE100_U_VERTICAL[:, 1], GHIA_RE100_U_VERTICAL[:, 0], "s", label="Ghia Re=100")

    plt.xlabel("u")
    plt.ylabel("y")
    plt.title("Lid-driven cavity: u along vertical centerline")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p1, dpi=200)
    plt.close()

    plt.figure()
    plt.plot(v_horizontal["x"], v_horizontal["v"], "o-", label="numerical")

    if re == 100:
        plt.plot(GHIA_RE100_V_HORIZONTAL[:, 0], GHIA_RE100_V_HORIZONTAL[:, 1], "s", label="Ghia Re=100")

    plt.xlabel("x")
    plt.ylabel("v")
    plt.title("Lid-driven cavity: v along horizontal centerline")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(p2, dpi=200)
    plt.close()

    print(f"Wrote: {p1}")
    print(f"Wrote: {p2}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate hex_lowmach_fv VTU output.")
    sub = parser.add_subparsers(dest="mode", required=True)

    channel = sub.add_parser("channel", help="Validate periodic body-force channel flow.")
    channel.add_argument("--vtu", required=True)
    channel.add_argument("--nu", type=float, required=True)
    channel.add_argument("--body-force-x", type=float, required=True)
    channel.add_argument("--ymin", type=float, required=True)
    channel.add_argument("--ymax", type=float, required=True)
    channel.add_argument("--nbins", type=int, default=64)
    channel.add_argument("--outdir", default=None)
    channel.add_argument("--plot", action="store_true")
    channel.set_defaults(func=validate_channel)

    cavity = sub.add_parser("cavity", help="Validate lid-driven cavity centerline profiles.")
    cavity.add_argument("--vtu", required=True)
    cavity.add_argument("--re", type=int, default=100)
    cavity.add_argument("--xmid", type=float, default=None)
    cavity.add_argument("--ymid", type=float, default=None)
    cavity.add_argument("--zmid", type=float, default=None)
    cavity.add_argument("--tol", type=float, default=None)
    cavity.add_argument("--outdir", default=None)
    cavity.add_argument("--plot", action="store_true")

    cavity.add_argument(
        "--lid-u",
        type=float,
        default=1.0,
        help="Lid velocity used when adding wall endpoint values for the Ghia comparison.",
    )

    cavity.add_argument(
        "--no-wall-points",
        action="store_true",
        help="Do not append known wall endpoint values to cavity centerline profiles.",
    )

    cavity.add_argument(
        "--no-z-filter",
        action="store_true",
        help="Average over all z layers instead of extracting the z-midplane.",
    )

    cavity.set_defaults(func=validate_cavity)

    args = parser.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())