#!/usr/bin/env python3
"""
Generate a frozen non-reacting 1D opposed-flow reference using Cantera.

This script is intended for validating the LowMachReact-Hex
counterflow_nonreacting_gmsh case.

It supports two reference modes:

1. Imposed-temperature / frozen reference, the historical default:
   - all reaction multipliers are set to zero
   - the 1D energy equation is disabled
   - a fixed temperature profile is imposed
   - for equal inlet temperatures, the result should be exactly isothermal

2. Non-isothermal solved-energy / frozen reference:
   - all reaction multipliers are set to zero
   - the 1D energy equation is enabled
   - inlet temperatures are imposed at the two inlets
   - Cantera solves the frozen energy equation with transport terms


Use mode 2 when comparing against a LowMachReact-Hex case with
`enable_energy = .true.` and unequal inlet temperatures.

Examples
--------

Isothermal / imposed-temperature historical reference:

  python tools/cantera_counterflow_reference.py \
    --t-fuel 300 \
    --t-ox 300 \
    --output cases/counterflow_nonreacting_gmsh/validation/cantera_reference_isothermal.csv

Non-isothermal solved-energy reference:

python tools/cantera_counterflow_reference.py \
  --solve-energy \
  --no-auto \
  --t-fuel 400.0 \
  --t-ox 300.0 \
  --u-fuel 0.5 \
  --u-ox 0.5 \
  --width 2.0 \
  --pressure 101325.0 \
  --constant-mdot-rho 1.171533 \
  --transport-model mixture-averaged \
  --output cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal.csv


python tools/cantera_counterflow_reference.py \
  --solve-energy \
  --no-auto \
  --t-fuel 400.0 \
  --t-ox 300.0 \
  --u-fuel 0.5 \
  --u-ox 0.5 \
  --width 2.0 \
  --pressure 101325.0 \
  --constant-mdot-rho 1.171533 \
  --transport-model mixture-averaged \
  --output cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal_massflux_matched.csv

python tools/cantera_counterflow_reference.py \
  --solve-energy \
  --no-auto \
  --t-fuel 400.0 \
  --t-ox 300.0 \
  --u-fuel 0.5 \
  --u-ox 0.5 \
  --width 2.0 \
  --pressure 101325.0 \
  --transport-model mixture-averaged \
  --output cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal_velocity_matched.csv


python tools/cantera_counterflow_reference.py \
  --solve-energy \
  --no-auto \
  --t-fuel 300.0 \
  --t-ox 300.0 \
  --u-fuel 0.5 \
  --u-ox 0.5 \
  --width 2.0 \
  --pressure 101325.0 \
  --transport-model mixture-averaged \
  --output cases/counterflow_nonreacting_gmsh/validation/cantera_reference_nonisothermal_velocity_matched.csv

Notes
-----

`--constant-mdot-rho` is useful when the FV case imposes equal inlet velocities
with a constant projection density. Cantera 1D inlets use mass fluxes, not
velocity boundary values, so this option sets

    mdot = constant_mdot_rho * |u|

for both inlets unless `--mdot-fuel` or `--mdot-ox` explicitly overrides a side.
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


def set_temperature_profile(flame, t_left: float, t_right: float, n_profile: int = 5) -> None:
    """
    Set the current temperature profile as an initial guess.

    Cantera profile positions are normalized from 0 to 1 across the domain.
    """
    z = np.linspace(0.0, 1.0, int(n_profile))
    t = np.linspace(float(t_left), float(t_right), int(n_profile))
    try:
        flame.set_profile("T", z, t)
    except Exception:
        pass


def set_fixed_temperature_profile(flame, t_left: float, t_right: float, n_profile: int = 5) -> None:
    """
    Set both the current T profile and the fixed-T profile used when the energy
    equation is disabled.
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


def set_soret_if_available(flame, enabled: bool) -> None:
    """Set Soret diffusion if the Cantera object exposes the option."""
    flow = getattr(flame, "flame", None)
    for obj in (flame, flow):
        if obj is None:
            continue
        if hasattr(obj, "soret_enabled"):
            try:
                obj.soret_enabled = bool(enabled)
                return
            except Exception:
                pass


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


def get_field_array(flame, names: tuple[str, ...], n: int) -> np.ndarray:
    """Return the first available flame field from a list of attribute names."""
    for name in names:
        if hasattr(flame, name):
            try:
                return array_or_nan(getattr(flame, name), n)
            except Exception:
                pass
    return np.full(n, np.nan)


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
        "--constant-mdot-rho",
        type=float,
        default=None,
        help=(
            "Use this constant density to convert inlet speeds to Cantera mdot. "
            "Useful for matching constant-density FV cases. Explicit --mdot-* overrides this."
        ),
    )
    parser.add_argument(
        "--solve-energy",
        action="store_true",
        help="Enable the Cantera 1D energy equation instead of imposing a fixed T profile.",
    )
    parser.add_argument(
        "--auto",
        action=argparse.BooleanOptionalAction,
        default=None,
        help=(
            "Use Cantera's automatic staged solve. Default is false for both "
            "solved-energy and imposed-temperature frozen-reference modes."
        ),
    )
    parser.add_argument(
        "--soret",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable Soret diffusion if supported by this Cantera version (default: false).",
    )
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
        help="Allowed max deviation from imposed temperature profile [K]. Used only without --solve-energy.",
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

    if args.constant_mdot_rho is not None:
        default_mdot_f = mass_flux_from_velocity_u(args.constant_mdot_rho, args.u_fuel)
        default_mdot_o = mass_flux_from_velocity_u(args.constant_mdot_rho, args.u_ox)
    else:
        default_mdot_f = mass_flux_from_velocity_u(rho_f, args.u_fuel)
        default_mdot_o = mass_flux_from_velocity_u(rho_o, args.u_ox)

    mdot_f = args.mdot_fuel if args.mdot_fuel is not None else default_mdot_f
    mdot_o = args.mdot_ox if args.mdot_ox is not None else default_mdot_o

    flame = ct.CounterflowDiffusionFlame(gas, width=args.width)

    flame.fuel_inlet.mdot = mdot_f
    flame.fuel_inlet.T = args.t_fuel
    flame.fuel_inlet.Y = y_fuel

    flame.oxidizer_inlet.mdot = mdot_o
    flame.oxidizer_inlet.T = args.t_ox
    flame.oxidizer_inlet.Y = y_ox

    try:
        flame.boundary_emissivities = 0.0, 0.0
    except Exception:
        pass
    try:
        flame.radiation_enabled = False
    except Exception:
        pass

    set_soret_if_available(flame, bool(args.soret))
    flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)

    auto_solve = False if args.auto is None else bool(args.auto)

    if args.solve_energy:
        set_energy_enabled_robust(flame, True)
        set_temperature_profile(flame, args.t_fuel, args.t_ox)
        mode = "frozen_solved_energy"
    else:
        set_energy_enabled_robust(flame, False)
        set_fixed_temperature_profile(flame, args.t_fuel, args.t_ox)
        mode = "frozen_imposed_temperature"

    try:
        flame.solve(args.loglevel, refine_grid=not args.no_refine, auto=auto_solve)
    except Exception as exc:
        mode_hint = (
            "For solved-energy mode, try --no-auto or --no-refine if the staged solve fails. "
            if args.solve_energy
            else "For imposed-temperature mode, try --no-refine, smaller --width, lower inlet speeds, or explicit --mdot-fuel/--mdot-ox. "
        )
        raise SystemExit(
            f"Cantera {mode} counterflow solve failed. {mode_hint}\n"
            f"Original error: {exc!r}"
        ) from exc

    set_energy_enabled_robust(flame, bool(args.solve_energy))

    x = np.asarray(flame.grid, dtype=float)
    T = np.asarray(flame.T, dtype=float)
    u = np.asarray(flame.velocity, dtype=float)
    n = len(x)

    T_expected = np.interp(x, [0.0, args.width], [args.t_fuel, args.t_ox])
    max_T_error = float(np.max(np.abs(T - T_expected)))

    if not args.solve_energy and max_T_error > args.tolerance_temperature:
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

    density = get_field_array(flame, ("density",), n)
    cp_mass = get_field_array(flame, ("cp_mass", "cp"), n)
    thermal_conductivity = get_field_array(flame, ("thermal_conductivity",), n)

    rows: dict[str, np.ndarray] = {
        "x": x,
        "temperature": T,
        "u": u,
        "rho": density,
        "rho_thermo": density,
        "cp": cp_mass,
        "lambda": thermal_conductivity,
    }
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
    print(f"  mode={mode} solve_energy={bool(args.solve_energy)} auto={auto_solve} refine={not args.no_refine}")
    print(f"  mdot_fuel={mdot_f:.12g} mdot_ox={mdot_o:.12g} kg/m^2/s  width={args.width:g} m")
    if args.constant_mdot_rho is not None:
        print(f"  mdot computed from constant_mdot_rho={args.constant_mdot_rho:.12g} kg/m^3")
    print(f"  rho_fuel_cantera={rho_f:.12g} rho_ox_cantera={rho_o:.12g} kg/m^3")
    print("  Reactions disabled: all reaction multipliers set to 0.")
    if args.solve_energy:
        print("  Energy equation enabled; inlet temperatures imposed, internal T solved.")
        print(f"  imposed linear-profile deviation diagnostic: max |T - T_linear|={max_T_error:.3e} K")
    else:
        print("  Energy equation disabled; fixed temperature profile imposed.")
        print(f"  fixed-profile check: max |T - T_expected|={max_T_error:.3e} K")
    print(f"  T_min={T.min():.12g} T_max={T.max():.12g}")
    print(f"  max_abs_heat_release_rate={hrr_max:.3e} W/m^3")
    print(f"  max_abs_net_production_rates={wdot_max:.3e}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
