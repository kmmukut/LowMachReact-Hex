#!/usr/bin/env python3
"""
Generate a steady 1D opposed-flow (counterflow diffusion) reference using Cantera,
with **all reaction rates set to zero** so transport is non-reacting.

The default inlet mass fractions and domain width match
``cases/counterflow_nonreacting_gmsh/case.nml`` and ``counterflow.geo`` (Lx = 2 m).

Requirements: ``cantera``, ``numpy``, ``pandas`` (``react_env``).

Example (from repo root with ``react_env`` active). GRI-30 is shipped with Cantera;
``--mech gri30.yaml`` is resolved from Cantera’s data path (no local copy required):

  python tools/cantera_counterflow_reference.py \\
    --output cases/counterflow_nonreacting_gmsh/validation/cantera_reference.csv

Then compare to the FV centerline:

  python tools/validate_flow.py counterflow \\
    --vtu cases/counterflow_nonreacting_gmsh/output/flow_050000.pvtu \\
    --reference-csv cases/counterflow_nonreacting_gmsh/validation/cantera_reference.csv

**Caveat:** The FV case uses constant ``rho`` and ``nu``; Cantera uses variable
thermodynamic density and mixture transport. Treat agreement as **qualitative**
unless you align models (see project validation docs).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def mass_flux_from_velocity_u(rho: float, u: float) -> float:
    """``mdot = rho * |u|`` [kg/m^2/s] for a stream normal to the inlet."""
    return float(abs(rho) * abs(u))


def disable_all_reactions(gas) -> None:
    for i in range(gas.n_reactions):
        gas.set_multiplier(0.0, i)


def load_gas(mech_arg: str, t_ref: float, pressure: float):
    """
    Load a ``Solution`` from a filesystem path, or from a mechanism name / file
    that Cantera resolves via its built-in data (e.g. ``gri30.yaml``).
    """
    import cantera as ct

    path = Path(mech_arg)
    if path.is_file():
        return ct.Solution(str(path.resolve()))

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
            f"If this is a bare filename, Cantera looks under its data directory (often {data!r}).\n"
            "Pass an absolute path to your YAML/CTI, or set CANTERA_DATA."
        ) from exc

    gas.TP = t_ref, pressure
    return gas


def mass_fraction_vector(gas, y_map: dict[str, float]) -> np.ndarray:
    y = np.zeros(gas.n_species)
    for name, val in y_map.items():
        y[gas.species_index(name)] = float(val)
    s = y.sum()
    if s <= 0.0:
        raise ValueError("inlet mass fractions must sum to a positive value")
    y /= s
    return y


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--mech",
        default="gri30.yaml",
        help="Mechanism: path to YAML/CTI, or a name Cantera resolves (default gri30.yaml = GRI-Mech 3.0 in Cantera data).",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=2.0,
        help="Distance between fuel and oxidizer inlets [m], match counterflow.geo Lx.",
    )
    parser.add_argument("--pressure", type=float, default=101325.0, help="Thermodynamic pressure [Pa].")
    parser.add_argument("--t-fuel", type=float, default=300.0, dest="t_fuel")
    parser.add_argument("--t-ox", type=float, default=300.0, dest="t_ox")
    parser.add_argument("--u-fuel", type=float, default=0.5, help="Used with inlet rho to set mdot if mdot not set.")
    parser.add_argument("--u-ox", type=float, default=0.5, help="Used with inlet rho to set mdot if mdot not set.")
    parser.add_argument("--mdot-fuel", type=float, default=None, help="Override fuel inlet mdot [kg/m^2/s].")
    parser.add_argument("--mdot-ox", type=float, default=None, help="Override oxidizer inlet mdot [kg/m^2/s].")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cases/counterflow_nonreacting_gmsh/validation/cantera_reference.csv"),
        help="Output CSV path.",
    )
    parser.add_argument("--loglevel", type=int, default=1, help="Cantera solve log level (0–5).")
    args = parser.parse_args()

    import cantera as ct

    gas = load_gas(args.mech, args.t_fuel, args.pressure)
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

    flame.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2, prune=0.03)
    try:
        flame.solve(args.loglevel, auto=True)
    except Exception as exc:
        raise SystemExit(
            "Cantera counterflow solve failed. Try smaller --width, lower mdot via --u-fuel/--u-ox, "
            "or explicit --mdot-fuel/--mdot-ox.\n"
            f"Original error: {exc!r}"
        ) from exc

    x = np.asarray(flame.grid, dtype=float)
    T = np.asarray(flame.T, dtype=float)
    u = np.asarray(flame.velocity, dtype=float)

    Y = np.asarray(flame.Y, dtype=float)
    if Y.ndim != 2 or Y.shape[0] != gas.n_species:
        raise SystemExit(f"unexpected flame.Y shape {Y.shape}; expected ({gas.n_species}, n_points)")

    rows: dict[str, np.ndarray] = {"x": x, "temperature": T, "u": u}
    for name in ("CH4", "O2", "N2"):
        k = gas.species_index(name)
        rows[f"Y_{name}"] = Y[k, :].copy()

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    print(f"Wrote {args.output.resolve()} with {len(out)} points.")
    print(f"  mdot_fuel={mdot_f:.6g} mdot_ox={mdot_o:.6g} kg/m^2/s  width={args.width} m")
    print("  Reactions disabled (all multipliers 0).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
