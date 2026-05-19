#!/usr/bin/env python3
"""
Plot optional counterflow centerline quantities extracted by tools/validate_flow.py.

Typical use from repo root:

  python tools/validate_flow.py counterflow \
    --vtu cases/counterflow_nonreacting_gmsh/output/flow_005000.pvtu \
    --reference-csv cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal_velocity_matched.csv

  python tools/plot_counterflow_extras.py \
    --profile cases/counterflow_nonreacting_gmsh/validation/counterflow_centerline.csv \
    --reference-csv cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal_velocity_matched.csv

The script only plots columns that exist.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PLOT_GROUPS = [
    ("thermo_density", ["rho", "rho_thermo"]),
    ("viscosity", ["nu"]),
    ("heat_capacity", ["cp"]),
    ("thermal_conductivity", ["lambda", "thermal_conductivity"]),
    ("enthalpy", ["enthalpy", "h"]),
    ("species_enthalpy_diffusion", ["species_enthalpy_diffusion"]),
    ("divergence", ["divergence"]),
    ("sum_Y", ["sum_Y"]),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--profile", required=True, help="counterflow_centerline.csv from validate_flow.py")
    p.add_argument("--reference-csv", default=None, help="Optional Cantera reference CSV")
    p.add_argument("--outdir", default=None, help="Output directory. Defaults to profile parent.")
    p.add_argument("--x-column", default="x")
    return p.parse_args()


def finite_range(values: np.ndarray) -> tuple[float, float] | None:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return None
    return float(np.min(a)), float(np.max(a))


def plot_group(df: pd.DataFrame, ref: pd.DataFrame | None, xcol: str, outdir: Path, name: str, cols: list[str]) -> bool:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    available = [c for c in cols if c in df.columns]
    ref_available = [c for c in cols if ref is not None and c in ref.columns]

    if not available and not ref_available:
        return False

    fig, ax = plt.subplots(figsize=(7, 4.5))

    for c in available:
        ax.plot(df[xcol], df[c], "-", label=f"FV {c}")

    if ref is not None and xcol in ref.columns:
        for c in ref_available:
            ax.plot(ref[xcol], ref[c], "--", linewidth=1.2, label=f"Cantera {c}")

    ax.set_xlabel(xcol)
    ax.set_ylabel(name)
    ax.set_title(f"Counterflow centerline: {name}")
    ax.grid(True)
    ax.legend(fontsize=8)
    fig.tight_layout()

    path = outdir / f"counterflow_extra_{name}.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"Wrote: {path}")
    return True


def write_summary(df: pd.DataFrame, ref: pd.DataFrame | None, xcol: str, outdir: Path) -> None:
    rows = []
    for c in df.columns:
        if c == xcol:
            continue
        try:
            arr = pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        except Exception:
            continue
        rng = finite_range(arr)
        if rng is None:
            continue
        rows.append({"source": "FV", "column": c, "min": rng[0], "max": rng[1]})

    if ref is not None:
        for c in ref.columns:
            if c == xcol:
                continue
            try:
                arr = pd.to_numeric(ref[c], errors="coerce").to_numpy(dtype=float)
            except Exception:
                continue
            rng = finite_range(arr)
            if rng is None:
                continue
            rows.append({"source": "Cantera", "column": c, "min": rng[0], "max": rng[1]})

    if rows:
        path = outdir / "counterflow_extra_summary.csv"
        pd.DataFrame(rows).to_csv(path, index=False)
        print(f"Wrote: {path}")


def main() -> int:
    args = parse_args()
    profile_path = Path(args.profile)
    outdir = Path(args.outdir) if args.outdir else profile_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(profile_path)
    if args.x_column not in df.columns:
        raise SystemExit(f"{profile_path} has no x-column {args.x_column!r}. Columns: {list(df.columns)}")

    ref = None
    if args.reference_csv:
        ref = pd.read_csv(args.reference_csv)

    print("FV columns:")
    print(", ".join(df.columns))
    if ref is not None:
        print("\nCantera columns:")
        print(", ".join(ref.columns))
    print()

    plotted = 0
    for name, cols in PLOT_GROUPS:
        if plot_group(df, ref, args.x_column, outdir, name, cols):
            plotted += 1

    write_summary(df, ref, args.x_column, outdir)

    if plotted == 0:
        print("No optional extra columns were found to plot.")
    else:
        print(f"Plotted {plotted} optional group(s).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
