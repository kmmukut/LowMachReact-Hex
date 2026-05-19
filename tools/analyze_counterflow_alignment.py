#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd


FV_CSV = Path("cases/counterflow_nonreacting_gmsh/validation/counterflow_centerline.csv")
REF_CSV = Path("cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal_velocity_matched.csv")


def find_column(df, candidates):
    for name in candidates:
        if name in df.columns:
            return name
    raise KeyError(f"None of these columns found: {candidates}\nAvailable: {list(df.columns)}")


def x_at_value(x, y, target):
    """
    Return x where y crosses target, using linear interpolation.
    Works best when y is monotonic near the crossing.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    f = y - target
    idx = np.where(f[:-1] * f[1:] <= 0.0)[0]

    if len(idx) == 0:
        raise ValueError(f"No crossing found for target={target}")

    # Use the crossing nearest the domain center.
    xmid = 0.5 * (x.min() + x.max())
    i = idx[np.argmin(np.abs(0.5 * (x[idx] + x[idx + 1]) - xmid))]

    x0, x1 = x[i], x[i + 1]
    y0, y1 = y[i], y[i + 1]

    if abs(y1 - y0) < 1.0e-300:
        return 0.5 * (x0 + x1)

    return x0 + (target - y0) * (x1 - x0) / (y1 - y0)


def fuel_based_mixture_fraction(df):
    """
    Fuel-based mixture fraction using CH4.

    Your Cantera reference script uses:
      fuel side: Y_CH4 = 0.20
      oxidizer side: Y_CH4 = 0.0
    """
    ych4_col = find_column(df, ["Y_CH4", "CH4"])
    y = np.asarray(df[ych4_col], dtype=float)

    y_fuel = y.max()
    y_ox = y.min()

    if abs(y_fuel - y_ox) < 1.0e-300:
        raise ValueError("Cannot compute mixture fraction: CH4 field is constant.")

    z = (y - y_ox) / (y_fuel - y_ox)
    return np.clip(z, 0.0, 1.0)


def mixing_thickness(x, z):
    """
    Thickness between Z=0.9 and Z=0.1.

    For this counterflow, Z decreases from fuel side to oxidizer side,
    so thickness = x(Z=0.1) - x(Z=0.9).
    """
    x_z09 = x_at_value(x, z, 0.9)
    x_z05 = x_at_value(x, z, 0.5)
    x_z01 = x_at_value(x, z, 0.1)
    return {
        "x_Z09": x_z09,
        "x_Z05": x_z05,
        "x_Z01": x_z01,
        "thickness_Z09_Z01": abs(x_z01 - x_z09),
    }


def interp_on_ref(x_src, y_src, x_ref):
    return np.interp(x_ref, x_src, y_src, left=np.nan, right=np.nan)


def error_stats(err):
    err = np.asarray(err, dtype=float)
    err = err[np.isfinite(err)]
    return {
        "L2": float(np.sqrt(np.mean(err**2))),
        "Linf": float(np.max(np.abs(err))),
    }


def main():
    fv = pd.read_csv(FV_CSV)
    ref = pd.read_csv(REF_CSV)

    x_fv_col = find_column(fv, ["x", "X"])
    x_ref_col = find_column(ref, ["x", "X"])

    u_fv_col = find_column(fv, ["u", "U", "velocity_x"])
    u_ref_col = find_column(ref, ["u", "U", "velocity_x"])

    x_fv = np.asarray(fv[x_fv_col], dtype=float)
    x_ref = np.asarray(ref[x_ref_col], dtype=float)

    u_fv = np.asarray(fv[u_fv_col], dtype=float)
    u_ref = np.asarray(ref[u_ref_col], dtype=float)

    # Sort just in case.
    fv_order = np.argsort(x_fv)
    ref_order = np.argsort(x_ref)

    x_fv = x_fv[fv_order]
    u_fv = u_fv[fv_order]
    fv = fv.iloc[fv_order].reset_index(drop=True)

    x_ref = x_ref[ref_order]
    u_ref = u_ref[ref_order]
    ref = ref.iloc[ref_order].reset_index(drop=True)

    # 1. Stagnation locations.
    xstag_fv = x_at_value(x_fv, u_fv, 0.0)
    xstag_ref = x_at_value(x_ref, u_ref, 0.0)
    dx_stag = xstag_fv - xstag_ref

    print("\nStagnation location")
    print("-------------------")
    print(f"x_stag_FV      = {xstag_fv:.12g}")
    print(f"x_stag_Cantera = {xstag_ref:.12g}")
    print(f"dx_stag        = {dx_stag:.12g}  # FV - Cantera")

    # 2. Mixture fraction thickness.
    z_fv = fuel_based_mixture_fraction(fv)
    z_ref = fuel_based_mixture_fraction(ref)

    thick_fv = mixing_thickness(x_fv, z_fv)
    thick_ref = mixing_thickness(x_ref, z_ref)

    print("\nFuel-based mixture-fraction layer")
    print("---------------------------------")
    print("FV:")
    for k, v in thick_fv.items():
        print(f"  {k:20s} = {v:.12g}")

    print("Cantera:")
    for k, v in thick_ref.items():
        print(f"  {k:20s} = {v:.12g}")

    print("\nThickness comparison")
    print("--------------------")
    print(f"FV thickness      = {thick_fv['thickness_Z09_Z01']:.12g}")
    print(f"Cantera thickness = {thick_ref['thickness_Z09_Z01']:.12g}")
    print(
        f"ratio FV/Cantera  = "
        f"{thick_fv['thickness_Z09_Z01'] / thick_ref['thickness_Z09_Z01']:.12g}"
    )

    # 3. Compare raw vs stagnation-aligned errors.
    print("\nRaw vs stagnation-aligned profile errors")
    print("----------------------------------------")

    # To align FV to Cantera, evaluate FV at x_ref + dx_stag.
    # If FV stagnation is to the right of Cantera, this shifts FV left.
    x_fv_sample_for_aligned = x_ref + dx_stag

    for col in ["u", "temperature", "Y_CH4", "Y_O2", "Y_N2"]:
        if col not in fv.columns or col not in ref.columns:
            continue

        y_fv = np.asarray(fv[col], dtype=float)
        y_ref = np.asarray(ref[col], dtype=float)

        fv_raw_on_ref = interp_on_ref(x_fv, y_fv, x_ref)
        fv_aligned_on_ref = interp_on_ref(x_fv, y_fv, x_fv_sample_for_aligned)

        raw = error_stats(fv_raw_on_ref - y_ref)
        aligned = error_stats(fv_aligned_on_ref - y_ref)

        print(
            f"{col:12s} "
            f"raw L2={raw['L2']:.6e}, raw Linf={raw['Linf']:.6e} | "
            f"aligned L2={aligned['L2']:.6e}, aligned Linf={aligned['Linf']:.6e}"
        )


if __name__ == "__main__":
    main()