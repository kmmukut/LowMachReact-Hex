#!/usr/bin/env python3
"""
Validation utilities for lowmach_react_hex VTU output.

Supports:
  1. Periodic body-force channel validation against analytic Poiseuille flow.
  2. Lid-driven cavity centerline extraction and comparison to Ghia Re=100 data.
  3. Non-reacting counterflow: merge parallel VTU pieces, extract x-centerline profiles;
     compare to a Cantera 1D CSV from ``tools/cantera_counterflow_reference.py``.

Requirements:
  pip install numpy pandas matplotlib meshio

Examples:

Channel:
  python tools/validate_flow.py channel \
    --vtu cases/channel_flow/output/flow_001000.vtu \
    --nu 1.0e-2 \
    --body-force-x 1.0e-3 \
    --ymin 0.0 \
    --ymax 1.0 \
    --nbins 64 \
    --plot

Cavity:
  python tools/validate_flow.py cavity \
    --vtu cases/lid_driven_cavity/output/flow_010000.vtu \
    --re 100 \
    --plot

Counterflow (pass the .pvtu; all flow_*_P*.vtu pieces in the same directory are merged).
  Writes CSV reports and images to the case-level ``validation/`` directory by default
  (use ``--no-plot`` to skip images):
  python tools/validate_flow.py counterflow \
    --vtu cases/counterflow_nonreacting_gmsh/output/VTK/flow_050000.pvtu

Counterflow with a Cantera reference CSV:
  python tools/validate_flow.py counterflow \
    --vtu cases/counterflow_nonreacting_gmsh/output/VTK/flow_005000.pvtu \
    --reference-csv cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import meshio
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Ghia et al. 2D lid-driven cavity benchmark data, Re = 100.
#
# u: vertical centerline, x = 0.5, values listed versus y.
# v: horizontal centerline, y = 0.5, values listed versus x.
# ---------------------------------------------------------------------------

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
    """
    Return the case-level validation directory for known solver output layouts.

    Supported layouts:

      cases/foo/output/flow_001000.vtu
        -> cases/foo/validation

      cases/foo/output/VTK/flow_001000.pvtu
        -> cases/foo/validation

    Otherwise, fall back to a validation directory next to the input VTU/PVTU.
    """
    vtu_path = Path(vtu_path)
    parent = vtu_path.parent

    # New layout:
    #   cases/<case>/output/VTK/flow_*.pvtu
    # should write validation products to:
    #   cases/<case>/validation/
    if parent.name == "VTK" and parent.parent.name == "output":
        return parent.parent.parent / "validation"

    # Old layout:
    #   cases/<case>/output/flow_*.vtu
    # should write validation products to:
    #   cases/<case>/validation/
    if parent.name == "output":
        return parent.parent / "validation"

    # Fallback for arbitrary VTU/PVTU locations.
    return parent / "validation"


def _hex_block_index(mesh: meshio.Mesh) -> int:
    for i, block in enumerate(mesh.cells):
        if block.type in ("hexahedron", "hex8"):
            return i
    return 0


def _cell_centers_and_data(mesh: meshio.Mesh) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    bi = _hex_block_index(mesh)
    cells = mesh.cells[bi].data
    points = mesh.points
    centers = points[cells].mean(axis=1)

    cell_data: dict[str, np.ndarray] = {}
    for name, blocks in mesh.cell_data.items():
        if len(blocks) > bi:
            cell_data[name] = np.asarray(blocks[bi])
    return centers, cell_data


def resolve_counterflow_piece_paths(vtu_path: str | Path) -> list[Path]:
    """
    Given a ``.pvtu`` path, return all ``flow_*_P*.vtu`` piece files in the same directory.
    Given a single ``.vtu``, return a one-element list.
    """
    path = Path(vtu_path)
    if path.suffix == ".pvtu":
        stem = path.stem
        pieces = sorted(path.parent.glob(f"{stem}_P*.vtu"))
        if not pieces:
            raise RuntimeError(
                f"No VTU pieces matched {path.parent}/{stem}_P*.vtu for {path}. "
                "Pass a .pvtu next to the per-rank piece files."
            )
        return pieces
    return [path]


def read_vtu_snapshot_merged(vtu_path: str | Path) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """
    Load one or more VTU pieces (merged if ``vtu_path`` is a ``.pvtu``) and return
    cell centers, velocity (nc, 3), and all cell-centered arrays.
    """
    paths = resolve_counterflow_piece_paths(vtu_path)
    centers_parts: list[np.ndarray] = []
    vel_parts: list[np.ndarray] = []
    data_parts: dict[str, list[np.ndarray]] = {}

    for p in paths:
        mesh = meshio.read(str(p))
        if not mesh.cells:
            raise RuntimeError(f"No cells found in {p}")
        centers, cell_data = _cell_centers_and_data(mesh)
        if "velocity" not in cell_data:
            raise RuntimeError(
                f"{p} has no CellData 'velocity'. Available: {list(cell_data.keys())}"
            )
        vel = np.asarray(cell_data["velocity"], dtype=float)
        if vel.ndim != 2 or vel.shape[1] != 3:
            raise RuntimeError(f"{p}: expected velocity shape (n,3), got {vel.shape}")
        centers_parts.append(centers)
        vel_parts.append(vel)
        for name, arr in cell_data.items():
            data_parts.setdefault(name, []).append(np.asarray(arr, dtype=float))

    centers_cat = np.vstack(centers_parts)
    velocity_cat = np.vstack(vel_parts)
    merged: dict[str, np.ndarray] = {}
    for name, parts in data_parts.items():
        merged[name] = np.concatenate(parts, axis=0)
    if merged["velocity"].shape[0] != centers_cat.shape[0]:
        raise RuntimeError("Internal merge error: velocity row count != centers row count")
    return centers_cat, velocity_cat, merged


def read_vtu_cell_data(vtu_path: str | Path):
    mesh = meshio.read(vtu_path)

    if not mesh.cells:
        raise RuntimeError(f"No cells found in {vtu_path}")

    centers, cell_data = _cell_centers_and_data(mesh)

    if "velocity" not in cell_data:
        raise RuntimeError(
            "VTU has no CellData array named 'velocity'. "
            f"Available arrays: {list(cell_data.keys())}"
        )

    velocity = np.asarray(cell_data["velocity"])
    return centers, velocity, cell_data


def estimate_centerline_tolerance(centers: np.ndarray) -> float:
    """
    Estimate a reasonable selection tolerance from cell-center spacing.
    """
    spacings = []

    for axis in range(3):
        vals = np.unique(np.round(centers[:, axis], 12))
        vals = np.sort(vals)
        diffs = np.diff(vals)
        span = float(vals[-1] - vals[0]) if vals.size > 1 else 0.0
        diffs = diffs[diffs > max(1.0e-10, 1.0e-6 * span)]
        if diffs.size > 0:
            spacings.append(float(np.median(diffs)))

    if not spacings:
        return 0.05

    return 0.60 * max(spacings)


def domain_midpoints(centers: np.ndarray):
    mins = centers.min(axis=0)
    maxs = centers.max(axis=0)
    mids = 0.5 * (mins + maxs)
    return mins, maxs, mids


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
    centers, velocity, _ = read_vtu_cell_data(args.vtu)

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
        raise RuntimeError(
            f"No cells found within tol={tol} of fixed coordinates {fixed}. "
            "Increase --tol or check domain coordinates."
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


def compare_to_benchmark(profile: pd.DataFrame, benchmark: np.ndarray, coord_col: str, value_col: str) -> pd.DataFrame:
    """
    Interpolate numerical profile to benchmark coordinates.
    """
    prof = profile.sort_values(coord_col)
    coord = prof[coord_col].to_numpy()
    value = prof[value_col].to_numpy()

    bcoord = benchmark[:, 0]
    bval = benchmark[:, 1]

    numerical_at_benchmark = np.interp(bcoord, coord, value)

    out = pd.DataFrame(
        {
            coord_col: bcoord,
            f"{value_col}_benchmark": bval,
            f"{value_col}_numeric_interp": numerical_at_benchmark,
            "error": numerical_at_benchmark - bval,
            "abs_error": np.abs(numerical_at_benchmark - bval),
        }
    )

    return out


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
    centers, velocity, _ = read_vtu_cell_data(args.vtu)

    outdir = Path(args.outdir) if args.outdir else default_validation_dir(args.vtu)
    outdir.mkdir(parents=True, exist_ok=True)

    mins, maxs, mids = domain_midpoints(centers)

    xmid = args.xmid if args.xmid is not None else float(mids[0])
    ymid = args.ymid if args.ymid is not None else float(mids[1])
    zmid = args.zmid if args.zmid is not None else float(mids[2])

    tol = args.tol if args.tol is not None else estimate_centerline_tolerance(centers)

    u = velocity[:, 0]
    v = velocity[:, 1]
    w = velocity[:, 2]

    u_vertical = nearest_line_profile(
        centers=centers,
        values=u,
        fixed={0: xmid, 2: zmid},
        varying_axis=1,
        tol=tol,
    ).rename(columns={"coord": "y", "value": "u"})

    v_horizontal = nearest_line_profile(
        centers=centers,
        values=v,
        fixed={1: ymid, 2: zmid},
        varying_axis=0,
        tol=tol,
    ).rename(columns={"coord": "x", "value": "v"})

    u_path = outdir / "cavity_u_vertical_centerline.csv"
    v_path = outdir / "cavity_v_horizontal_centerline.csv"

    u_vertical.to_csv(u_path, index=False)
    v_horizontal.to_csv(v_path, index=False)

    summary = {
        "vtu": str(args.vtu),
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


def plot_cavity(outdir: Path, u_vertical: pd.DataFrame, v_horizontal: pd.DataFrame, re: int) -> None:
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


def validate_counterflow(args: argparse.Namespace) -> None:
    """
    Extract near-1D profiles along x through (ymid, zmid) for non-reacting counterflow checks.
    Optional CSV reference (e.g. from Cantera 1D) for qualitative overlay metrics.
    """
    centers, velocity, cell_data = read_vtu_snapshot_merged(args.vtu)

    outdir = Path(args.outdir) if args.outdir else default_validation_dir(args.vtu)
    outdir.mkdir(parents=True, exist_ok=True)

    mins, maxs, mids = domain_midpoints(centers)
    ymid = float(args.ymid) if args.ymid is not None else float(mids[1])
    zmid = float(args.zmid) if args.zmid is not None else float(mids[2])
    tol = float(args.tol) if args.tol is not None else estimate_centerline_tolerance(centers)

    u = velocity[:, 0]
    v = velocity[:, 1]
    w = velocity[:, 2]

    mask = (np.abs(centers[:, 1] - ymid) <= tol) & (np.abs(centers[:, 2] - zmid) <= tol)
    if not np.any(mask):
        raise RuntimeError(
            f"No cells within tol={tol} of line y={ymid}, z={zmid}. "
            "Try larger --tol or set --ymid/--zmid to domain interior."
        )

    x = centers[mask, 0]
    df = pd.DataFrame(
        {
            "x": x,
            "u": u[mask],
            "v": v[mask],
            "w": w[mask],
        }
    )

    if "temperature" in cell_data:
        df["temperature"] = np.asarray(cell_data["temperature"]).reshape(-1)[mask]
    if "pressure" in cell_data:
        df["pressure"] = np.asarray(cell_data["pressure"]).reshape(-1)[mask]
    if "sum_Y" in cell_data:
        df["sum_Y"] = np.asarray(cell_data["sum_Y"]).reshape(-1)[mask]

    # Optional scalar fields useful for non-isothermal Cantera/FV comparisons.
    # Only fields present in the VTU are copied. The Cantera reference generator
    # writes columns such as rho, rho_thermo, cp, and lambda, so keeping the same
    # names here allows the existing overlap-based comparison logic to include them.
    optional_scalar_fields = [
        "rho",
        "nu",
        "rho_thermo",
        "cp",
        "lambda",
        "thermal_conductivity",
        "h",
        "enthalpy",
        "species_enthalpy_diffusion",
        "qrad",
        "divergence",
    ]
    for key in optional_scalar_fields:
        if key in cell_data and key not in df.columns:
            arr = np.asarray(cell_data[key]).reshape(-1)
            if arr.shape[0] == centers.shape[0]:
                df[key] = arr[mask]

    # Alias thermal_conductivity to lambda when the solver writes only the
    # descriptive field name but the Cantera reference uses lambda.
    if "lambda" not in df.columns and "thermal_conductivity" in df.columns:
        df["lambda"] = df["thermal_conductivity"]

    for key in sorted(cell_data.keys()):
        if key.startswith("Y_") and key not in df.columns:
            df[key] = np.asarray(cell_data[key]).reshape(-1)[mask]

    df["x_round"] = df["x"].round(9)
    agg_spec: dict[str, str] = {"x": "mean", "u": "mean", "v": "mean", "w": "mean"}
    for c in df.columns:
        if c not in ("x", "x_round", "u", "v", "w"):
            agg_spec[c] = "mean"
    df = df.groupby("x_round", as_index=False).agg(agg_spec).sort_values("x")
    if "x_round" in df.columns:
        df = df.drop(columns=["x_round"])

    yf = float(args.fuel_y_inlet)
    yo = float(args.fuel_y_oxidizer)
    den = yf - yo
    fuel_col = args.fuel_species_column
    if fuel_col in df.columns and abs(den) > 1.0e-30:
        df["Z_mixture"] = (df[fuel_col] - yo) / den
    else:
        df["Z_mixture"] = np.nan

    profile_path = outdir / "counterflow_centerline.csv"
    df.to_csv(profile_path, index=False)

    summary: dict[str, float | str] = {
        "vtu": str(args.vtu),
        "domain_x_min": float(mins[0]),
        "domain_x_max": float(maxs[0]),
        "ymid_used": ymid,
        "zmid_used": zmid,
        "tol_used": tol,
        "n_points_centerline": int(len(df)),
        "max_abs_v_on_line": float(np.max(np.abs(df["v"]))),
        "max_abs_w_on_line": float(np.max(np.abs(df["w"]))),
        "u_min": float(np.min(df["u"])),
        "u_max": float(np.max(df["u"])),
        "u_left": float(df["u"].iloc[0]),
        "u_right": float(df["u"].iloc[-1]),
        "u_has_sign_change": str(bool(np.any(df["u"].to_numpy()[:-1] * df["u"].to_numpy()[1:] <= 0.0))),
    }
    if "divergence" in cell_data:
        div_all = np.asarray(cell_data["divergence"], dtype=float).reshape(-1)
        summary["max_abs_divergence_all_cells"] = float(np.max(np.abs(div_all)))
        summary["rms_divergence_all_cells"] = float(np.sqrt(np.mean(div_all**2)))
    if "temperature" in df.columns:
        summary["T_min"] = float(np.min(df["temperature"]))
        summary["T_max"] = float(np.max(df["temperature"]))
    if "sum_Y" in df.columns:
        summary["sum_Y_min"] = float(np.min(df["sum_Y"]))
        summary["sum_Y_max"] = float(np.max(df["sum_Y"]))

    ref_path = Path(args.reference_csv) if args.reference_csv else None
    ref_df: pd.DataFrame | None = None
    compare_cols: list[str] = []
    if ref_path is not None:
        if not ref_path.is_file():
            outp = ref_path.resolve()
            raise SystemExit(
                f"reference CSV not found: {ref_path}\n"
                f"(resolved: {outp})\n\n"
                "Create it first with Cantera (from the repository root, with react_env active), e.g.:\n"
                "  python tools/cantera_counterflow_reference.py \\\n"
                f"    --output {ref_path}\n\n"
                "(Default ``--mech gri30.yaml`` uses GRI-30 from Cantera’s data; only pass ``--mech`` if you use a custom YAML/CTI path.)\n\n"
                "If you already ran that, use an absolute --reference-csv path, or run validate_flow "
                "from the repo root so relative paths resolve."
            )
        ref = pd.read_csv(ref_path)
        if "x" not in ref.columns:
            raise RuntimeError("reference CSV must include an 'x' column (same axis as this case).")
        ref_df = ref
        compare_cols = [c for c in ref.columns if c != "x" and c in df.columns]
        if not compare_cols:
            raise RuntimeError(
                f"No overlapping columns between reference and profile besides 'x'. "
                f"ref columns: {list(ref.columns)}, profile: {list(df.columns)}"
            )
        x_num = df["x"].to_numpy()
        rows = []
        xref = ref["x"].to_numpy()
        for col in compare_cols:
            y_num = df[col].to_numpy()
            y_ref = ref[col].to_numpy()
            y_interp = np.interp(xref, x_num, y_num)
            err = y_interp - y_ref
            rows.append(
                {
                    "column": col,
                    "L2": float(np.sqrt(np.mean(err**2))),
                    "Linf": float(np.max(np.abs(err))),
                    "rel_L2": float(np.sqrt(np.mean(err**2)) / max(np.sqrt(np.mean(y_ref**2)), 1.0e-300)),
                }
            )
        cmp_df = pd.DataFrame(rows)
        cmp_path = outdir / "counterflow_vs_reference_errors.csv"
        cmp_df.to_csv(cmp_path, index=False)
        summary["reference_csv"] = str(ref_path)
        summary["compared_columns"] = ",".join(compare_cols)
        print(f"\nWrote: {cmp_path}")
        print(cmp_df.to_string(index=False))
        if summary["u_has_sign_change"] == "False":
            print(
                "\nWARNING: centerline u does not change sign. A Cantera CounterflowDiffusionFlame "
                "reference has opposed-flow/stagnation behavior, so this FV run is not a quantitative "
                "counterflow validation."
            )
        if (
            float(summary.get("max_abs_divergence_all_cells", 0.0)) > 1.0e-3
            or float(summary.get("rms_divergence_all_cells", 0.0)) > 1.0e-5
        ):
            print(
                "WARNING: large cell divergence detected. Check diagnostics.csv and the mass balance "
                "before interpreting Cantera comparison errors."
            )

    summary_path = outdir / "counterflow_summary.txt"
    with summary_path.open("w") as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")

    print("\nCounterflow centerline extraction")
    print("-----------------------------------")
    for key, value in summary.items():
        print(f"{key}: {value}")
    print(f"\nWrote: {profile_path}")
    print(f"Wrote: {summary_path}")

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
        ax_u, ax_T = axes[0, 0], axes[0, 1]
        ax_Z, ax_Y = axes[1, 0], axes[1, 1]

        ax_u.plot(df["x"], df["u"], "-", color="C0", label="FV u")
        if ref_df is not None and "u" in ref_df.columns:
            ax_u.plot(ref_df["x"], ref_df["u"], "--", color="C0", linewidth=1.0, alpha=0.85, label="Cantera u")
        ax_u.set_ylabel("u")
        ax_u.grid(True)
        ax_u.legend(fontsize=8)

        if "temperature" in df.columns:
            ax_T.plot(df["x"], df["temperature"], "-", color="C1", label="FV T")
            if ref_df is not None and "temperature" in ref_df.columns:
                ax_T.plot(
                    ref_df["x"],
                    ref_df["temperature"],
                    "--",
                    color="C1",
                    linewidth=1.0,
                    alpha=0.85,
                    label="Cantera T",
                )
            ax_T.set_ylabel("T")
            ax_T.legend(fontsize=8)
        ax_T.grid(True)

        if np.all(np.isfinite(df["Z_mixture"].to_numpy())):
            ax_Z.plot(df["x"], df["Z_mixture"], "-", color="C2", label="FV Z")
            if ref_df is not None and fuel_col in ref_df.columns and abs(den) > 1.0e-30:
                zr = (ref_df[fuel_col].to_numpy(dtype=float) - yo) / den
                ax_Z.plot(ref_df["x"], zr, "--", color="C2", linewidth=1, alpha=0.85, label="Cantera Z")
            ax_Z.set_ylabel("Z (fuel-based)")
            ax_Z.legend(fontsize=8)
        ax_Z.set_xlabel("x")
        ax_Z.grid(True)

        ycols = [c for c in df.columns if c.startswith("Y_")]
        for c in ycols[:3]:
            ax_Y.plot(df["x"], df[c], "-", label=f"FV {c}")
            if ref_df is not None and c in ref_df.columns:
                ax_Y.plot(ref_df["x"], ref_df[c], "--", linewidth=1, alpha=0.85, label=f"Cantera {c}")
        ax_Y.set_xlabel("x")
        ax_Y.set_ylabel("Y_k")
        ax_Y.legend(fontsize=6, ncol=2)
        ax_Y.grid(True)

        title = "Counterflow: FV centerline vs Cantera 1D" if ref_df is not None else "Counterflow: centerline profiles"
        plt.suptitle(title)
        plt.tight_layout()
        plot_path = outdir / "counterflow_centerline.png"
        plt.savefig(plot_path, dpi=200)
        plt.close()
        print(f"Wrote: {plot_path}")

        if ref_df is not None and compare_cols:
            fig2, axes2 = plt.subplots(2, 2, figsize=(9, 7), sharex=True)
            flat2 = np.ravel(axes2)
            for i, col in enumerate(compare_cols[:4]):
                ax = flat2[i]
                xn = df["x"].to_numpy()
                yn = df[col].to_numpy()
                xr = ref_df["x"].to_numpy()
                yr = ref_df[col].to_numpy()
                y_interp = np.interp(xr, xn, yn)
                ax.plot(xr, y_interp - yr, "-", color="C3")
                ax.axhline(0.0, color="0.5", linewidth=0.8, linestyle=":")
                ax.set_ylabel(f"Δ{col} (FV−Cantera)")
                ax.set_xlabel("x")
                ax.grid(True)
            nplot = min(len(compare_cols), 4)
            for j in range(nplot, 4):
                flat2[j].set_visible(False)
            plt.suptitle("Counterflow: pointwise residual on reference x (FV interpolated − Cantera)")
            plt.tight_layout()
            res_path = outdir / "counterflow_vs_reference_residuals.png"
            plt.savefig(res_path, dpi=200)
            plt.close()
            print(f"Wrote: {res_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate lowmach_react_hex VTU output.")
    sub = parser.add_subparsers(dest="mode", required=True)

    channel = sub.add_parser("channel", help="Validate periodic body-force channel flow.")
    channel.add_argument("--vtu", required=True)
    channel.add_argument("--nu", type=float, required=True)
    channel.add_argument("--body-force-x", type=float, required=True)
    channel.add_argument("--ymin", type=float, required=True)
    channel.add_argument("--ymax", type=float, required=True)
    channel.add_argument("--nbins", type=int, default=64)
    channel.add_argument("--outdir", default=None, help="Validation output directory. Default: case-level validation/ folder.")
    channel.add_argument("--plot", action="store_true")
    channel.set_defaults(func=validate_channel)

    cavity = sub.add_parser("cavity", help="Validate lid-driven cavity centerline profiles.")
    cavity.add_argument("--vtu", required=True)
    cavity.add_argument("--re", type=int, default=100)
    cavity.add_argument("--xmid", type=float, default=None)
    cavity.add_argument("--ymid", type=float, default=None)
    cavity.add_argument("--zmid", type=float, default=None)
    cavity.add_argument("--tol", type=float, default=None)
    cavity.add_argument("--outdir", default=None, help="Validation output directory. Default: case-level validation/ folder.")
    cavity.add_argument("--plot", action="store_true")
    cavity.set_defaults(func=validate_cavity)

    cf = sub.add_parser(
        "counterflow",
        help="Extract x-centerline profiles from counterflow VTU (merge parallel .pvtu pieces).",
    )
    cf.add_argument("--vtu", required=True, help="Path to flow_*.pvtu or a single flow_*_P*.vtu piece.")
    cf.add_argument("--ymid", type=float, default=None)
    cf.add_argument("--zmid", type=float, default=None)
    cf.add_argument("--tol", type=float, default=None)
    cf.add_argument(
        "--reference-csv",
        default=None,
        help="Optional CSV with column 'x' plus fields to compare (e.g. u, temperature, Y_CH4).",
    )
    cf.add_argument(
        "--fuel-y-inlet",
        type=float,
        default=0.2,
        dest="fuel_y_inlet",
        help="Fuel-stream inlet mass fraction of the fuel species (for Z_mixture).",
    )
    cf.add_argument(
        "--fuel-y-oxidizer",
        type=float,
        default=0.0,
        dest="fuel_y_oxidizer",
        help="Oxidizer-stream inlet mass fraction of the fuel species (for Z_mixture).",
    )
    cf.add_argument(
        "--fuel-species-column",
        default="Y_CH4",
        help="Column name in VTU/CSV for fuel tracer (default Y_CH4).",
    )
    cf.add_argument("--outdir", default=None, help="Validation output directory. Default: case-level validation/ folder.")
    cf.add_argument(
        "--plot",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write counterflow_centerline.png (default: true). Use --no-plot to skip.",
    )
    cf.set_defaults(func=validate_counterflow)

    args = parser.parse_args()
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
