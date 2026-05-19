#!/usr/bin/env python3
"""
Generate a frozen, isothermal 1D opposed-flow reference using Cantera.

This script is intended for validating the FV counterflow_nonreacting_gmsh case.
It differs from the standard CounterflowDiffusionFlame example in three important
ways:

1. All reaction multipliers are set to zero.
2. The energy equation is disabled on the 1D flow domain.
3. The fixed temperature profile used by Cantera when energy is disabled is
   explicitly set to a flat inlet-temperature profile.

For equal inlet temperatures, the output temperature should therefore be flat.
If it is not flat, the script raises an error instead of silently writing a bad
reference CSV.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def mass_flux_from_velocity_u(rho: float, u: float) -> float:
    """Return mdot = rho * |u| [kg/m^2/s] for a stream normal to the inlet."""
    return float(abs(rho) * abs(u))


def disable_all_reactions(gas) -> None:
    """Set every reaction-rate multiplier to zero on a Cantera Solution."""
    for i in range(gas.n_reactions):
        gas.set_multiplier(0.0, i)


def load_gas(mech_arg: str, t_ref: float, pressure: float):
    """
    Load a Cantera Solution from a filesystem path or from a mechanism name that
    Cantera resolves from its data path, for example 'gri30.yaml'.
    """
    import cantera as ct

    path = Path(mech_arg)
    if path.is_file():
        gas = ct.Solution(str(path.resolve()))
    else:
        try:
            gas = ct.Solution(mech_arg)
        except Exception as exc:
            data = ""
            try:
                data = ct.get_data_directory()
            except Exception:
                pass
            raise SystemExit(
                f"Could not load mechanism {mech_arg!r} ({exc}).\n"
                f"If this is a bare filename, Cantera looks under its data directory "
                f"(often {data!r}).\n"
                "Pass an absolute path to your YAML/CTI, or set CANTERA_DATA."
            ) from exc

    gas.TP = t_ref, pressure
    return gas


def mass_fraction_vector(gas, y_map: dict[str, float]) -> np.ndarray:
    """Build and normalize a Cantera mass-fraction vector from species names."""
    y = np.zeros(gas.n_species)
    for name, val in y_map.items():
        try:
            k = gas.species_index(name)
        except ValueError as exc:
            raise SystemExit(f"Species {name!r} is not present in this mechanism.") from exc
        y[k] = float(val)

    s = y.sum()
    if s <= 0.0:
        raise ValueError("inlet mass fractions must sum to a positive value")
    y /= s
    return y


def set_energy_enabled_robust(flame, enabled: bool) -> None:
    """
    Disable/enable the energy equation across Cantera versions.

    Current Cantera versions generally expose this on flame.flame.energy_enabled,
    while some examples also use flame.energy_enabled. Try both and verify when
    possible.
    """
    changed = False

    if hasattr(flame, "energy_enabled"):
        try:
            flame.energy_enabled = enabled
            changed = True
        except Exception:
            pass

    flow = getattr(flame, "flame", None)
    if flow is not None and hasattr(flow, "energy_enabled"):
        try:
            flow.energy_enabled = enabled
            changed = True
        except Exception:
            pass

    if not changed:
        raise RuntimeError(
            "Could not set energy_enabled on this Cantera flame object. "
            "Inspect your Cantera version's 1D flame API."
        )

    states = []
    if hasattr(flame, "energy_enabled"):
        try:
            states.append(bool(flame.energy_enabled))
        except Exception:
            pass
    if flow is not None and hasattr(flow, "energy_enabled"):
        try:
            states.append(bool(flow.energy_enabled))
        except Exception:
            pass

    if states and any(state != enabled for state in states):
        raise RuntimeError(f"Failed to set all readable energy_enabled flags to {enabled}.")


def set_flat_temperature_profile(flame, t_left: float, t_right: float, n_profile: int = 5) -> None:
    """
    Set both the solver's current T profile and the fixed-T profile used when the
    energy equation is disabled.

    Cantera profile positions are normalized from 0 to 1 across the domain.
    For equal inlet temperatures this creates a flat profile. For unequal inlet
    temperatures this creates a linear imposed profile.
    """
    z = np.linspace(0.0, 1.0, int(n_profile))
    t = np.linspace(float(t_left), float(t_right), int(n_profile))

    try:
        flame.set_profile("T", z, t)
    except Exception:
        pass

    flow = getattr(flame, "flame", None)
    fixed_profile_set = False

    if flow is not None and hasattr(flow, "set_fixed_temp_profile"):
        try:
            flow.set_fixed_temp_profile(z, t)
            fixed_profile_set = True
        except Exception:
            pass

    if hasattr(flame, "set_fixed_temp_profile"):
        try:
            flame.set_fixed_temp_profile(z, t)
            fixed_profile_set = True
        except Exception:
            pass

    if not fixed_profile_set:
        raise RuntimeError(
            "Could not call set_fixed_temp_profile. Your Cantera version may use "
            "a different API for imposed-temperature 1D flames."
        )


def array_or_nan(values, n: int) -> np.ndarray:
    """Return a 1D array, or NaNs if the requested Cantera field is unavailable."""
    try:
        a = np.asarray(values, dtype=float)
        if a.ndim == 0:
            return np.full(n, float(a))
        if a.size == n:
            return a.reshape(n)
    except Exception:
        pass
    return np.full(n, np.nan)


def max_abs_or_nan(values) -> float:
    try:
        a = np.asarray(values, dtype=float)
        if a.size == 0:
            return float("nan")
        return float(np.nanmax(np.abs(a)))
    except Exception:
        return float("nan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--mech",
        default="gri30.yaml",
        help="Mechanism path/name resolved by Cantera (default: gri30.yaml).",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=2.0,
        help="Distance between fuel and oxidizer inlets [m]. Match counterflow.geo Lx.",
    )
    parser.add_argument("--pressure", type=float, default=101325.0, help="Pressure [Pa].")
    parser.add_argument("--t-fuel", type=float, default=300.0, dest="t_fuel")
    parser.add_argument("--t-ox", type=float, default=300.0, dest="t_ox")
    parser.add_argument("--u-fuel", type=float, default=0.5, help="Fuel inlet speed magnitude [m/s].")
    parser.add_argument("--u-ox", type=float, default=0.5, help="Oxidizer inlet speed magnitude [m/s].")
    parser.add_argument("--mdot-fuel", type=float, default=None, help="Override fuel mdot [kg/m^2/s].")
    parser.add_argument("--mdot-ox", type=float, default=None, help="Override oxidizer mdot [kg/m^2/s].")
    parser.add_argument(
        "--transport-model",
        default=None,
        choices=[None, "mixture-averaged", "multicomponent"],
        help="Optional Cantera transport model override.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cases/counterflow_nonreacting_gmsh/validation/cantera_reference.csv"),
        help="Output CSV path.",
    )
    parser.add_argument("--loglevel", type=int, default=1, help="Cantera solve log level, 0-5.")
    parser.add_argument(
        "--no-refine",
        action="store_true",
        help="Disable grid refinement. Useful for debugging fixed profiles.",
    )
    parser.add_argument(
        "--tolerance-temperature",
        type=float,
        default=1.0e-8,
        help="Allowed max deviation from imposed temperature profile [K].",
    )
    parser.add_argument(
        "--write-diagnostics",
        action="store_true",
        help="Include heat_release_rate and net production diagnostics in the CSV when available.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    import cantera as ct

    gas = load_gas(args.mech, args.t_fuel, args.pressure)
    if args.transport_model is not None:
        try:
            gas.transport_model = args.transport_model
        except Exception as exc:
            raise SystemExit(f"Could not set transport_model={args.transport_model!r}: {exc}") from exc

    disable_all_reactions(gas)

    y_fuel = mass_fraction_vector(gas, {"CH4": 0.20, "O2": 0.05, "N2": 0.75})
    y_ox = mass_fraction_vector(gas, {"CH4": 0.0, "O2": 0.23, "N2": 0.77})

    gas.TPY = args.t_fuel, args.pressure, y_fuel
    rho_f = gas.density
    gas.TPY = args.t_ox, args.pressure, y_ox
    rho_o = gas.density

    mdot_f = args.mdot_fuel if args.mdot_fuel is not None else mass_flux_from_velocity_u(rho_f, args.u_fuel)
    mdot_o = args.mdot_ox if args.mdot_ox is not None else mass_flux_from_velocity_u(rho_o, args.u_ox)

    flame = ct.CounterflowDiffusionFlame(gas, width=args.width)

    flame.fuel_inlet.mdot = mdot_f
    flame.fuel_inlet.T = args.t_fuel
    flame.fuel_inlet.Y = y_fuel

    flame.oxidizer_inlet.mdot = mdot_o
    flame.oxidizer_inlet.T = args.t_ox
    flame.oxidizer_inlet.Y = y_ox

    flame.boundary_emissivities = 0.0, 0.0
    flame.radiation_enabled = False

    set_energy_enabled_robust(flame, False)
    set_flat_temperature_profile(flame, args.t_fuel, args.t_ox)

    flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

    try:
        flame.solve(args.loglevel, refine_grid=not args.no_refine, auto=False)
    except Exception as exc:
        raise SystemExit(
            "Cantera frozen/isothermal counterflow solve failed. Try --no-refine, "
            "smaller --width, lower inlet speeds, or explicit --mdot-fuel/--mdot-ox.\n"
            f"Original error: {exc!r}"
        ) from exc

    set_energy_enabled_robust(flame, False)

    x = np.asarray(flame.grid, dtype=float)
    T = np.asarray(flame.T, dtype=float)
    u = np.asarray(flame.velocity, dtype=float)
    n = len(x)

    T_expected = np.interp(x, [0.0, args.width], [args.t_fuel, args.t_ox])
    max_T_error = float(np.max(np.abs(T - T_expected)))
    if max_T_error > args.tolerance_temperature:
        i = int(np.argmax(np.abs(T - T_expected)))
        raise SystemExit(
            "Cantera reference is not following the imposed fixed-temperature profile.\n"
            f"  max |T - T_expected| = {max_T_error:.12e} K at x = {x[i]:.12g}\n"
            f"  T = {T[i]:.12g}, T_expected = {T_expected[i]:.12g}\n"
            "This usually means the energy equation was not actually disabled, the fixed-T "
            "profile was not applied, or an old CSV is being plotted."
        )

    Y = np.asarray(flame.Y, dtype=float)
    if Y.ndim != 2 or Y.shape[0] != gas.n_species:
        raise SystemExit(f"unexpected flame.Y shape {Y.shape}; expected ({gas.n_species}, n_points)")

    rows: dict[str, np.ndarray] = {"x": x, "temperature": T, "u": u}
    for name in ("CH4", "O2", "N2"):
        k = gas.species_index(name)
        rows[f"Y_{name}"] = Y[k, :].copy()

    if args.write_diagnostics:
        rows["heat_release_rate"] = array_or_nan(getattr(flame, "heat_release_rate", np.nan), n)
        try:
            wdot = np.asarray(flame.net_production_rates, dtype=float)
            if wdot.ndim == 2 and wdot.shape[1] == n:
                for name in ("CH4", "O2", "N2"):
                    rows[f"wdot_{name}"] = wdot[gas.species_index(name), :].copy()
        except Exception:
            pass

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    hrr_max = max_abs_or_nan(getattr(flame, "heat_release_rate", np.nan))
    wdot_max = max_abs_or_nan(getattr(flame, "net_production_rates", np.nan))

    print(f"Wrote {args.output.resolve()} with {len(out)} points.")
    print(f"  mdot_fuel={mdot_f:.12g} mdot_ox={mdot_o:.12g} kg/m^2/s  width={args.width:g} m")
    print(f"  rho_fuel={rho_f:.12g} rho_ox={rho_o:.12g} kg/m^3")
    print("  Reactions disabled: all reaction multipliers set to 0.")
    print("  Energy equation disabled; fixed temperature profile imposed.")
    print(f"  T_min={T.min():.12g} T_max={T.max():.12g} max_T_error={max_T_error:.3e} K")
    print(f"  max_abs_heat_release_rate={hrr_max:.3e} W/m^3")
    print(f"  max_abs_net_production_rates={wdot_max:.3e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())