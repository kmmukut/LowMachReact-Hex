#!/usr/bin/env python3
"""Audit coupled variable-density species/mass/enthalpy diagnostics.

Patch 070 adds this offline audit tool. It reads existing solver diagnostic CSVs
and checks whether aggregate mass, transported species mass, per-species
balances, and reconciled energy metrics tell a consistent story.

Typical usage from repository root:

    python tools/diagnostics/audit_coupled_transport_conservation.py \
        --output-dir cases/rectangle_2D/output

To write a CSV audit time history:

    python tools/diagnostics/audit_coupled_transport_conservation.py \
        --output-dir cases/rectangle_2D/output \
        --write-csv

Exit status is 0 for PASS and PASS_WITH_WARNINGS. Use --strict-warnings to
return nonzero when WARN metrics are present.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


SPECIES_ENERGY_FILE = "species_energy_conservation.csv"
SPECIES_INTEGRALS_FILE = "species_integrals.csv"
ENTHALPY_FILE = "enthalpy_energy_budget.csv"
CONTINUITY_FILE = "variable_density_continuity_residual.csv"
AUDIT_FILE = "coupled_transport_audit.csv"

TINY = 1.0e-300


@dataclass
class MetricResult:
    group: str
    label: str
    column: str
    value: float | None
    tolerance_text: str
    status: str
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit coupled variable-density species/mass/enthalpy diagnostics."
    )
    location = parser.add_mutually_exclusive_group(required=False)
    location.add_argument(
        "--output-dir",
        help="Simulation output directory containing diagnostics/.",
    )
    location.add_argument(
        "--diagnostics-dir",
        help="Directory containing diagnostic CSV files directly.",
    )

    parser.add_argument(
        "--species-mass-rel-pass",
        type=float,
        default=1.0e-6,
        help="PASS threshold for |sum_k rhoY_k - rho|/rho. Default: 1e-6.",
    )
    parser.add_argument(
        "--species-mass-rel-fail",
        type=float,
        default=1.0e-4,
        help="FAIL threshold for |sum_k rhoY_k - rho|/rho. Default: 1e-4.",
    )
    parser.add_argument(
        "--boundary-flux-rel-pass",
        type=float,
        default=1.0e-6,
        help="PASS threshold for species-flux-sum versus mass-flux relative difference. Default: 1e-6.",
    )
    parser.add_argument(
        "--boundary-flux-rel-fail",
        type=float,
        default=1.0e-4,
        help="FAIL threshold for species-flux-sum versus mass-flux relative difference. Default: 1e-4.",
    )
    parser.add_argument(
        "--species-balance-rel-pass",
        type=float,
        default=1.0e-5,
        help="PASS threshold for summed species-balance defect versus aggregate mass defect. Default: 1e-5.",
    )
    parser.add_argument(
        "--species-balance-rel-fail",
        type=float,
        default=1.0e-3,
        help="FAIL threshold for summed species-balance defect versus aggregate mass defect. Default: 1e-3.",
    )
    parser.add_argument(
        "--species-sum-pass",
        type=float,
        default=1.0e-8,
        help="PASS threshold for species-sum minus one metrics. Default: 1e-8.",
    )
    parser.add_argument(
        "--species-sum-fail",
        type=float,
        default=1.0e-6,
        help="FAIL threshold for species-sum minus one metrics. Default: 1e-6.",
    )
    parser.add_argument(
        "--energy-tol",
        type=float,
        default=1.0e-8,
        help="Strict tolerance for reconciled energy closure metrics. Default: 1e-8.",
    )
    parser.add_argument(
        "--continuity-integral-abs-tol",
        type=float,
        default=1.0e-6,
        help="Optional absolute tolerance for continuity integral context. Default: 1e-6.",
    )
    parser.add_argument(
        "--continuity-relative-l2-warn",
        type=float,
        default=1.0e-3,
        help="Optional PASS threshold for continuity local relative L2. Default: 1e-3.",
    )
    parser.add_argument(
        "--continuity-relative-l2-fail",
        type=float,
        default=1.0e-2,
        help="Optional FAIL threshold for continuity local relative L2. Default: 1e-2.",
    )
    parser.add_argument(
        "--allow-missing-energy",
        action="store_true",
        help="Treat missing enthalpy_energy_budget.csv as SKIP instead of FAIL.",
    )
    parser.add_argument(
        "--allow-missing-species-integrals",
        action="store_true",
        help="Treat missing species_integrals.csv as SKIP instead of WARN.",
    )
    parser.add_argument(
        "--write-csv",
        action="store_true",
        help="Write diagnostics/coupled_transport_audit.csv time history.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit machine-readable JSON instead of text.",
    )
    parser.add_argument(
        "--strict-warnings",
        action="store_true",
        help="Return nonzero if any WARN metrics are present.",
    )
    args = parser.parse_args()

    for name in (
        "species_mass_rel_pass",
        "species_mass_rel_fail",
        "boundary_flux_rel_pass",
        "boundary_flux_rel_fail",
        "species_balance_rel_pass",
        "species_balance_rel_fail",
        "species_sum_pass",
        "species_sum_fail",
        "energy_tol",
        "continuity_integral_abs_tol",
        "continuity_relative_l2_warn",
        "continuity_relative_l2_fail",
    ):
        value = getattr(args, name)
        if not math.isfinite(value) or value <= 0.0:
            parser.error(f"{name.replace('_', '-')} must be a positive finite number")

    if args.species_mass_rel_fail < args.species_mass_rel_pass:
        parser.error("species-mass-rel-fail must be >= species-mass-rel-pass")
    if args.boundary_flux_rel_fail < args.boundary_flux_rel_pass:
        parser.error("boundary-flux-rel-fail must be >= boundary-flux-rel-pass")
    if args.species_balance_rel_fail < args.species_balance_rel_pass:
        parser.error("species-balance-rel-fail must be >= species-balance-rel-pass")
    if args.species_sum_fail < args.species_sum_pass:
        parser.error("species-sum-fail must be >= species-sum-pass")
    if args.continuity_relative_l2_fail < args.continuity_relative_l2_warn:
        parser.error("continuity-relative-l2-fail must be >= continuity-relative-l2-warn")

    return args


def resolve_diagnostics_dir(args: argparse.Namespace) -> Path:
    if args.diagnostics_dir:
        return Path(args.diagnostics_dir).resolve()
    if args.output_dir:
        return (Path(args.output_dir).resolve() / "diagnostics")

    cwd = Path.cwd().resolve()
    if cwd.name == "diagnostics":
        return cwd
    if (cwd / "diagnostics").is_dir():
        return cwd / "diagnostics"
    return cwd


def parse_float(text: Any) -> float:
    return float(str(text).strip().replace("D", "E").replace("d", "e"))


def parse_step(text: Any) -> int | None:
    try:
        return int(round(parse_float(text)))
    except Exception:
        return None


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        rows = [row for row in reader if row and any((value or "").strip() for value in row.values())]
    return list(reader.fieldnames), rows


def existing_csv(diagnostics_dir: Path, filename: str) -> Path | None:
    candidates = [diagnostics_dir / filename]
    parent = diagnostics_dir.parent
    if parent != diagnostics_dir:
        candidates.append(parent / filename)
    for path in candidates:
        if path.is_file():
            return path
    return None


def index_by_step(rows: Iterable[dict[str, str]]) -> dict[int, dict[str, str]]:
    out: dict[int, dict[str, str]] = {}
    for row in rows:
        step = parse_step(row.get("step", ""))
        if step is not None:
            out[step] = row
    return out


def find_key(row: dict[str, str], aliases: tuple[str, ...]) -> str | None:
    exact = {key: key for key in row.keys()}
    lower = {key.strip().lower(): key for key in row.keys()}
    for alias in aliases:
        if alias in exact:
            return exact[alias]
        key = alias.strip().lower()
        if key in lower:
            return lower[key]
    return None


def get_value(row: dict[str, str] | None, aliases: tuple[str, ...]) -> tuple[str | None, float | None]:
    if row is None:
        return None, None
    key = find_key(row, aliases)
    if key is None:
        return None, None
    raw = row.get(key, "")
    try:
        return key, parse_float(raw)
    except Exception:
        return key, None


def rel_diff(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), TINY)


def grade(value: float | None, pass_tol: float, fail_tol: float, label: str, column: str, group: str) -> MetricResult:
    tolerance_text = f"PASS <= {pass_tol:.1e}; WARN <= {fail_tol:.1e}; FAIL > {fail_tol:.1e}"
    if value is None or not math.isfinite(value):
        return MetricResult(group, label, column, value, tolerance_text, "FAIL", "missing or non-finite value")
    if abs(value) <= pass_tol:
        return MetricResult(group, label, column, value, tolerance_text, "PASS", "within pass threshold")
    if abs(value) <= fail_tol:
        return MetricResult(group, label, column, value, tolerance_text, "WARN", "above pass threshold but below fail threshold")
    return MetricResult(group, label, column, value, tolerance_text, "FAIL", "exceeds fail threshold")


def strict(value: float | None, tol: float, label: str, column: str, group: str, missing_status: str = "FAIL") -> MetricResult:
    tolerance_text = f"<= {tol:.1e}"
    if value is None or not math.isfinite(value):
        return MetricResult(group, label, column, value, tolerance_text, missing_status, "missing or non-finite value")
    if abs(value) <= tol:
        return MetricResult(group, label, column, value, tolerance_text, "PASS", "within tolerance")
    return MetricResult(group, label, column, value, tolerance_text, "FAIL", "exceeds tolerance")


def sum_species_balance_defects(row: dict[str, str] | None) -> float | None:
    if row is None:
        return None
    total = 0.0
    found = False
    pattern = re.compile(r"^Y\d+_balance_defect_since_previous$", re.IGNORECASE)
    for key, raw in row.items():
        if not pattern.match(key.strip()):
            continue
        try:
            total += parse_float(raw)
            found = True
        except Exception:
            return None
    return total if found else None


def extract_audit_row(
    aggregate_row: dict[str, str],
    species_integral_row: dict[str, str] | None,
    enthalpy_row: dict[str, str] | None,
    continuity_row: dict[str, str] | None,
) -> dict[str, float | int | str | None]:
    step = parse_step(aggregate_row.get("step", ""))
    _time_key, time = get_value(aggregate_row, ("time",))

    _k, total_mass = get_value(aggregate_row, ("total_mass",))
    _k, species_mass_sum = get_value(aggregate_row, ("transported_species_mass_sum",))
    _k, net_mass_flux = get_value(aggregate_row, ("net_boundary_mass_flux",))
    _k, species_flux_sum = get_value(aggregate_row, ("net_boundary_species_mass_flux_sum",))
    _k, mass_balance_defect = get_value(aggregate_row, ("mass_balance_defect_since_previous",))
    _k, species_sum_abs_max = get_value(aggregate_row, ("species_sum_minus_one_abs_max",))
    _k, species_sum_l2 = get_value(aggregate_row, ("species_sum_minus_one_l2",))
    _k, rho_h_integral = get_value(aggregate_row, ("rho_h_integral",))

    sum_species_defect = sum_species_balance_defects(species_integral_row)

    _k, direct_energy = get_value(enthalpy_row, ("relative_last_energy_update_balance_defect",))
    _k, rel_output_recon = get_value(
        enthalpy_row,
        (
            "rel_output_recon_defect",
            "relative_output_state_budget_defect_after_density_reconciliation",
        ),
    )
    _k, rel_operator_recon = get_value(
        enthalpy_row,
        (
            "rel_operator_recon_defect",
            "relative_operator_consistent_budget_defect_after_density_reconciliation",
        ),
    )
    _k, enthalpy_rho_h_integral = get_value(enthalpy_row, ("rho_h_integral",))
    _k, operator_rho_h_integral = get_value(enthalpy_row, ("operator_consistent_rho_h_integral",))

    _k, continuity_integral = get_value(continuity_row, ("integral_drho_dt_plus_div_mass_flux_dV",))
    _k, continuity_rel_l2 = get_value(continuity_row, ("relative_conservative_residual_l2",))

    species_mass_minus_total = None
    rel_species_mass_minus_total = None
    if total_mass is not None and species_mass_sum is not None:
        species_mass_minus_total = species_mass_sum - total_mass
        rel_species_mass_minus_total = rel_diff(species_mass_sum, total_mass)

    species_flux_minus_mass_flux = None
    rel_species_flux_minus_mass_flux = None
    if net_mass_flux is not None and species_flux_sum is not None:
        species_flux_minus_mass_flux = species_flux_sum - net_mass_flux
        rel_species_flux_minus_mass_flux = rel_diff(species_flux_sum, net_mass_flux)

    species_balance_minus_mass_balance = None
    rel_species_balance_minus_mass_balance = None
    if mass_balance_defect is not None and sum_species_defect is not None:
        species_balance_minus_mass_balance = sum_species_defect - mass_balance_defect
        rel_species_balance_minus_mass_balance = rel_diff(sum_species_defect, mass_balance_defect)

    rho_h_integral_difference = None
    rel_rho_h_integral_difference = None
    if rho_h_integral is not None and enthalpy_rho_h_integral is not None:
        rho_h_integral_difference = rho_h_integral - enthalpy_rho_h_integral
        rel_rho_h_integral_difference = rel_diff(rho_h_integral, enthalpy_rho_h_integral)

    return {
        "step": step,
        "time": time,
        "total_mass": total_mass,
        "transported_species_mass_sum": species_mass_sum,
        "species_mass_sum_minus_total_mass": species_mass_minus_total,
        "rel_species_mass_sum_minus_total_mass": rel_species_mass_minus_total,
        "net_boundary_mass_flux": net_mass_flux,
        "net_boundary_species_mass_flux_sum": species_flux_sum,
        "boundary_species_flux_sum_minus_mass_flux": species_flux_minus_mass_flux,
        "rel_boundary_species_flux_sum_minus_mass_flux": rel_species_flux_minus_mass_flux,
        "mass_balance_defect_since_previous": mass_balance_defect,
        "sum_species_balance_defect_since_previous": sum_species_defect,
        "species_balance_sum_minus_mass_balance_defect": species_balance_minus_mass_balance,
        "rel_species_balance_sum_minus_mass_balance_defect": rel_species_balance_minus_mass_balance,
        "species_sum_minus_one_abs_max": species_sum_abs_max,
        "species_sum_minus_one_l2": species_sum_l2,
        "rho_h_integral_species_energy_csv": rho_h_integral,
        "rho_h_integral_enthalpy_budget_csv": enthalpy_rho_h_integral,
        "rho_h_integral_difference_between_csvs": rho_h_integral_difference,
        "rel_rho_h_integral_difference_between_csvs": rel_rho_h_integral_difference,
        "operator_consistent_rho_h_integral": operator_rho_h_integral,
        "relative_last_energy_update_balance_defect": direct_energy,
        "rel_output_recon_defect": rel_output_recon,
        "rel_operator_recon_defect": rel_operator_recon,
        "integral_drho_dt_plus_div_mass_flux_dV": continuity_integral,
        "relative_conservative_residual_l2": continuity_rel_l2,
    }


def evaluate_audit(row: dict[str, Any], args: argparse.Namespace, has_species_integrals: bool, has_energy: bool) -> list[MetricResult]:
    results: list[MetricResult] = []

    results.append(
        grade(
            row.get("rel_species_mass_sum_minus_total_mass"),
            args.species_mass_rel_pass,
            args.species_mass_rel_fail,
            "transported species mass sum versus total mass",
            "rel_species_mass_sum_minus_total_mass",
            "coupled mass/species",
        )
    )
    results.append(
        grade(
            row.get("rel_boundary_species_flux_sum_minus_mass_flux"),
            args.boundary_flux_rel_pass,
            args.boundary_flux_rel_fail,
            "species boundary flux sum versus mass flux",
            "rel_boundary_species_flux_sum_minus_mass_flux",
            "coupled mass/species",
        )
    )

    if has_species_integrals:
        results.append(
            grade(
                row.get("rel_species_balance_sum_minus_mass_balance_defect"),
                args.species_balance_rel_pass,
                args.species_balance_rel_fail,
                "summed species balance defects versus aggregate mass defect",
                "rel_species_balance_sum_minus_mass_balance_defect",
                "coupled mass/species",
            )
        )
    else:
        status = "SKIP" if args.allow_missing_species_integrals else "WARN"
        results.append(
            MetricResult(
                "coupled mass/species",
                "summed species balance defects versus aggregate mass defect",
                "species_integrals.csv",
                None,
                "n/a",
                status,
                "species_integrals.csv missing",
            )
        )

    results.append(
        grade(
            row.get("species_sum_minus_one_abs_max"),
            args.species_sum_pass,
            args.species_sum_fail,
            "species sum minus one absolute max",
            "species_sum_minus_one_abs_max",
            "species sum",
        )
    )
    results.append(
        grade(
            row.get("species_sum_minus_one_l2"),
            args.species_sum_pass,
            args.species_sum_fail,
            "species sum minus one L2",
            "species_sum_minus_one_l2",
            "species sum",
        )
    )

    if has_energy:
        results.append(
            strict(
                row.get("relative_last_energy_update_balance_defect"),
                args.energy_tol,
                "direct energy update",
                "relative_last_energy_update_balance_defect",
                "energy",
            )
        )
        results.append(
            strict(
                row.get("rel_output_recon_defect"),
                args.energy_tol,
                "reconciled output-state rho*h budget",
                "rel_output_recon_defect",
                "energy",
            )
        )
        results.append(
            strict(
                row.get("rel_operator_recon_defect"),
                args.energy_tol,
                "reconciled operator-consistent rho*h budget",
                "rel_operator_recon_defect",
                "energy",
            )
        )
        results.append(
            grade(
                row.get("rel_rho_h_integral_difference_between_csvs"),
                1.0e-12,
                1.0e-9,
                "rho*h integral agreement between aggregate and enthalpy CSVs",
                "rel_rho_h_integral_difference_between_csvs",
                "energy",
            )
        )
    else:
        status = "SKIP" if args.allow_missing_energy else "FAIL"
        results.append(
            MetricResult("energy", "reconciled energy metrics", ENTHALPY_FILE, None, "n/a", status, "enthalpy_energy_budget.csv missing")
        )

    if row.get("integral_drho_dt_plus_div_mass_flux_dV") is not None:
        results.append(
            strict(
                row.get("integral_drho_dt_plus_div_mass_flux_dV"),
                args.continuity_integral_abs_tol,
                "conservative continuity global integral",
                "integral_drho_dt_plus_div_mass_flux_dV",
                "continuity context",
            )
        )
    if row.get("relative_conservative_residual_l2") is not None:
        results.append(
            grade(
                row.get("relative_conservative_residual_l2"),
                args.continuity_relative_l2_warn,
                args.continuity_relative_l2_fail,
                "conservative continuity relative L2 residual",
                "relative_conservative_residual_l2",
                "continuity context",
            )
        )

    return results


def overall_status(results: list[MetricResult]) -> str:
    if any(result.status == "FAIL" for result in results):
        return "FAIL"
    if any(result.status == "WARN" for result in results):
        return "PASS_WITH_WARNINGS"
    return "PASS"


def load_all(diagnostics_dir: Path) -> dict[str, tuple[list[str], list[dict[str, str]]]]:
    data: dict[str, tuple[list[str], list[dict[str, str]]]] = {}
    required = existing_csv(diagnostics_dir, SPECIES_ENERGY_FILE)
    if required is None:
        raise FileNotFoundError(f"missing required {SPECIES_ENERGY_FILE} in {diagnostics_dir}")
    data[SPECIES_ENERGY_FILE] = read_csv_rows(required)

    for filename in (SPECIES_INTEGRALS_FILE, ENTHALPY_FILE, CONTINUITY_FILE):
        path = existing_csv(diagnostics_dir, filename)
        if path is not None:
            data[filename] = read_csv_rows(path)
    return data


def build_audit_rows(data: dict[str, tuple[list[str], list[dict[str, str]]]]) -> list[dict[str, Any]]:
    aggregate_rows = data[SPECIES_ENERGY_FILE][1]
    species_by_step = index_by_step(data.get(SPECIES_INTEGRALS_FILE, ([], []))[1])
    enthalpy_by_step = index_by_step(data.get(ENTHALPY_FILE, ([], []))[1])
    continuity_by_step = index_by_step(data.get(CONTINUITY_FILE, ([], []))[1])

    audit_rows: list[dict[str, Any]] = []
    for aggregate_row in aggregate_rows:
        step = parse_step(aggregate_row.get("step", ""))
        if step is None:
            continue
        audit_rows.append(
            extract_audit_row(
                aggregate_row,
                species_by_step.get(step),
                enthalpy_by_step.get(step),
                continuity_by_step.get(step),
            )
        )
    return audit_rows


def format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.16e}"
    return str(value)


def write_audit_csv(diagnostics_dir: Path, rows: list[dict[str, Any]], args: argparse.Namespace, data: dict[str, tuple[list[str], list[dict[str, str]]]]) -> Path:
    path = diagnostics_dir / AUDIT_FILE
    if not rows:
        raise ValueError("no audit rows to write")

    has_species_integrals = SPECIES_INTEGRALS_FILE in data
    has_energy = ENTHALPY_FILE in data
    fieldnames = list(rows[0].keys()) + [
        "coupled_audit_overall_status",
        "coupled_audit_fail_count",
        "coupled_audit_warn_count",
    ]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            results = evaluate_audit(row, args, has_species_integrals, has_energy)
            out = dict(row)
            out["coupled_audit_overall_status"] = overall_status(results)
            out["coupled_audit_fail_count"] = sum(1 for result in results if result.status == "FAIL")
            out["coupled_audit_warn_count"] = sum(1 for result in results if result.status == "WARN")
            writer.writerow({key: format_value(out.get(key)) for key in fieldnames})
    return path


def result_to_dict(result: MetricResult) -> dict[str, Any]:
    return {
        "group": result.group,
        "label": result.label,
        "column": result.column,
        "value": result.value,
        "tolerance": result.tolerance_text,
        "status": result.status,
        "message": result.message,
    }


def print_report(
    diagnostics_dir: Path,
    latest_row: dict[str, Any],
    results: list[MetricResult],
    wrote_csv: Path | None,
) -> None:
    overall = overall_status(results)
    print("Coupled transport conservation audit")
    print("====================================")
    print(f"diagnostics_dir: {diagnostics_dir}")
    print(f"step: {latest_row.get('step')}  time: {latest_row.get('time')}")
    print(f"overall: {overall}")
    if wrote_csv is not None:
        print(f"audit_csv: {wrote_csv}")
    print()

    print("Key coupled quantities from latest row:")
    for key in (
        "total_mass",
        "transported_species_mass_sum",
        "rel_species_mass_sum_minus_total_mass",
        "net_boundary_mass_flux",
        "net_boundary_species_mass_flux_sum",
        "rel_boundary_species_flux_sum_minus_mass_flux",
        "mass_balance_defect_since_previous",
        "sum_species_balance_defect_since_previous",
        "rel_species_balance_sum_minus_mass_balance_defect",
        "species_sum_minus_one_abs_max",
        "species_sum_minus_one_l2",
        "relative_last_energy_update_balance_defect",
        "rel_output_recon_defect",
        "rel_operator_recon_defect",
        "integral_drho_dt_plus_div_mass_flux_dV",
        "relative_conservative_residual_l2",
    ):
        value = latest_row.get(key)
        if value is None:
            continue
        if isinstance(value, float):
            print(f"  {key}: {value:.8e}")
        else:
            print(f"  {key}: {value}")
    print()

    current_group: str | None = None
    for result in results:
        if result.group != current_group:
            current_group = result.group
            print(f"{current_group.capitalize()}:")
        value_text = "n/a" if result.value is None else f"{result.value:.8e}"
        print(
            f"  [{result.status:4}] {result.label}\n"
            f"         column={result.column} value={value_text} tolerance={result.tolerance_text}\n"
            f"         {result.message}"
        )
    print()
    print("Interpretation:")
    print("  This audit cross-checks whether the aggregate mass budget, transported species")
    print("  budgets, species-sum constraint, and reconciled rho*h energy budget are")
    print("  mutually consistent. It is an offline validation tool; it does not change")
    print("  solver numerics or diagnostic CSV schemas.")


def main() -> int:
    args = parse_args()
    diagnostics_dir = resolve_diagnostics_dir(args)

    try:
        data = load_all(diagnostics_dir)
        audit_rows = build_audit_rows(data)
    except Exception as exc:  # noqa: BLE001
        if args.json:
            print(json.dumps({"diagnostics_dir": str(diagnostics_dir), "overall": "FAIL", "error": str(exc)}, indent=2))
        else:
            print(f"Coupled transport conservation audit FAILED: {exc}", file=sys.stderr)
        return 1

    if not audit_rows:
        if args.json:
            print(json.dumps({"diagnostics_dir": str(diagnostics_dir), "overall": "FAIL", "error": "no audit rows"}, indent=2))
        else:
            print("Coupled transport conservation audit FAILED: no audit rows", file=sys.stderr)
        return 1

    wrote_csv: Path | None = None
    if args.write_csv:
        try:
            wrote_csv = write_audit_csv(diagnostics_dir, audit_rows, args, data)
        except Exception as exc:  # noqa: BLE001
            if args.json:
                print(json.dumps({"diagnostics_dir": str(diagnostics_dir), "overall": "FAIL", "error": f"could not write audit CSV: {exc}"}, indent=2))
            else:
                print(f"Coupled transport conservation audit FAILED: could not write audit CSV: {exc}", file=sys.stderr)
            return 1

    latest_row = audit_rows[-1]
    has_species_integrals = SPECIES_INTEGRALS_FILE in data
    has_energy = ENTHALPY_FILE in data
    results = evaluate_audit(latest_row, args, has_species_integrals, has_energy)
    status = overall_status(results)

    if args.json:
        payload = {
            "diagnostics_dir": str(diagnostics_dir),
            "overall": status,
            "strict_warnings": bool(args.strict_warnings),
            "audit_csv": str(wrote_csv) if wrote_csv else None,
            "latest": latest_row,
            "results": [result_to_dict(result) for result in results],
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print_report(diagnostics_dir, latest_row, results, wrote_csv)

    if status == "FAIL":
        return 1
    if status == "PASS_WITH_WARNINGS" and args.strict_warnings:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
