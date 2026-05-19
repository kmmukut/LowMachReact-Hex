#!/usr/bin/env python3
"""Validate variable-density low-Mach diagnostic CSV files.

This checker encodes the primary validation metrics documented in
``doc_src/validation_metrics.md``. It intentionally checks the reconciled and
operator-consistent energy metrics, not older unreconciled explanatory columns.

Patch 069b uses category-specific tolerances:

* energy direct/reconciled closure: strict, roundoff-scale relative tolerance
* projection residuals: projection-scale relative tolerance
* continuity: global integral closure is pass/fail; local relative L2 residual
  is a quality metric with PASS/WARN/FAIL thresholds

Typical usage from repository root:

    python tools/diagnostics/check_variable_density_validation.py \
        --output-dir cases/rectangle_2D/output

Useful variants:

    python tools/diagnostics/check_variable_density_validation.py \
        --diagnostics-dir cases/rectangle_2D/output/diagnostics

    python tools/diagnostics/check_variable_density_validation.py \
        --output-dir cases/rectangle_2D/output \
        --energy-only

Exit status is 0 for PASS and PASS_WITH_WARNINGS. Use --strict-warnings to
return nonzero when WARN metrics are present.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class MetricResult:
    group: str
    label: str
    filename: str
    column: str | None
    value: float | None
    tolerance_text: str
    status: str
    message: str


OPTIONAL_CONTEXT_COLUMNS: dict[str, tuple[str, ...]] = {
    "enthalpy_energy_budget.csv": (
        "step",
        "time",
        "relative_last_energy_update_balance_defect",
        "rel_output_recon_defect",
        "rel_operator_recon_defect",
        "cumulative_operator_balance_defect",
        "output_state_density_reconciliation_rate",
        "operator_consistent_density_reconciliation_rate",
        "output_state_budget_defect_after_density_reconciliation",
        "operator_consistent_budget_defect_after_density_reconciliation",
    ),
    "variable_density_compatibility.csv": (
        "step",
        "time",
        "relative_divu_minus_S_projection_l2",
        "relative_divu_minus_S_projection_max",
        "divu_minus_S_projection_l2",
        "divu_minus_S_current_l2",
        "S_current_minus_S_projection_l2",
    ),
    "variable_density_continuity_residual.csv": (
        "step",
        "time",
        "integral_drho_dt_plus_div_mass_flux_dV",
        "conservative_residual_l2",
        "relative_conservative_residual_l2",
    ),
}


ENERGY_FILE = "enthalpy_energy_budget.csv"
PROJECTION_FILE = "variable_density_compatibility.csv"
CONTINUITY_FILE = "variable_density_continuity_residual.csv"


ENERGY_METRICS: tuple[tuple[str, tuple[str, ...]], ...] = (
    (
        "direct energy update",
        ("relative_last_energy_update_balance_defect",),
    ),
    (
        "reconciled output-state rho*h budget",
        (
            "rel_output_recon_defect",
            "relative_output_state_budget_defect_after_density_reconciliation",
        ),
    ),
    (
        "reconciled operator-consistent rho*h budget",
        (
            "rel_operator_recon_defect",
            "relative_operator_consistent_budget_defect_after_density_reconciliation",
        ),
    ),
)


PROJECTION_METRICS: tuple[tuple[str, tuple[str, ...], str], ...] = (
    (
        "projection-source relative L2 residual",
        ("relative_divu_minus_S_projection_l2",),
        "projection_l2_tol",
    ),
    (
        "projection-source relative max residual",
        ("relative_divu_minus_S_projection_max",),
        "projection_max_tol",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check accepted variable-density validation metrics from CSV diagnostics."
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
        "--tol",
        type=float,
        default=None,
        help=(
            "Legacy global tolerance override for strict metrics. When supplied, "
            "it becomes the default for energy, projection, and continuity-integral "
            "checks unless a category-specific option is also supplied."
        ),
    )
    parser.add_argument(
        "--energy-tol",
        type=float,
        default=None,
        help="Relative tolerance for strict energy metrics. Default: 1e-8.",
    )
    parser.add_argument(
        "--projection-tol",
        type=float,
        default=None,
        help="Relative tolerance for both projection metrics. Default: 1e-5.",
    )
    parser.add_argument(
        "--projection-l2-tol",
        type=float,
        default=None,
        help="Relative tolerance for projection-source L2 residual. Default: --projection-tol or 1e-5.",
    )
    parser.add_argument(
        "--projection-max-tol",
        type=float,
        default=None,
        help="Relative tolerance for projection-source max residual. Default: --projection-tol or 1e-5.",
    )
    parser.add_argument(
        "--continuity-integral-abs-tol",
        type=float,
        default=None,
        help="Absolute tolerance for integral d(rho)/dt + div(rho u). Default: 1e-6.",
    )
    parser.add_argument(
        "--continuity-relative-l2-warn",
        type=float,
        default=1.0e-3,
        help="PASS threshold for local continuity relative L2. Default: 1e-3.",
    )
    parser.add_argument(
        "--continuity-relative-l2-fail",
        type=float,
        default=1.0e-2,
        help="FAIL threshold for local continuity relative L2. Values between warn and fail are WARN. Default: 1e-2.",
    )
    parser.add_argument(
        "--continuity-tol",
        type=float,
        default=None,
        help=(
            "Backward-compatible alias for --continuity-relative-l2-fail. "
            "Prefer --continuity-relative-l2-warn/--continuity-relative-l2-fail."
        ),
    )

    parser.add_argument(
        "--energy-only",
        action="store_true",
        help="Only check enthalpy_energy_budget.csv energy metrics.",
    )
    parser.add_argument(
        "--skip-projection",
        action="store_true",
        help="Do not check variable_density_compatibility.csv metrics.",
    )
    parser.add_argument(
        "--skip-continuity",
        action="store_true",
        help="Do not check variable_density_continuity_residual.csv metrics.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit a machine-readable JSON summary instead of text.",
    )
    parser.add_argument(
        "--strict-warnings",
        action="store_true",
        help="Return nonzero if any WARN metrics are present.",
    )
    parser.add_argument(
        "--allow-missing-optional",
        action="store_true",
        help=(
            "Treat missing projection/continuity files as SKIP instead of FAIL. "
            "Energy metrics remain required."
        ),
    )
    args = parser.parse_args()

    global_tol = args.tol
    args.energy_tol = choose_tol(args.energy_tol, global_tol, 1.0e-8)

    projection_default = choose_tol(args.projection_tol, global_tol, 1.0e-5)
    args.projection_l2_tol = choose_tol(args.projection_l2_tol, None, projection_default)
    args.projection_max_tol = choose_tol(args.projection_max_tol, None, projection_default)

    args.continuity_integral_abs_tol = choose_tol(
        args.continuity_integral_abs_tol,
        global_tol,
        1.0e-6,
    )
    if args.continuity_tol is not None:
        args.continuity_relative_l2_fail = args.continuity_tol

    for name in (
        "energy_tol",
        "projection_l2_tol",
        "projection_max_tol",
        "continuity_integral_abs_tol",
        "continuity_relative_l2_warn",
        "continuity_relative_l2_fail",
    ):
        value = getattr(args, name)
        if not math.isfinite(value) or value <= 0.0:
            parser.error(f"{name.replace('_', '-')} must be a positive finite number")

    if args.continuity_relative_l2_fail < args.continuity_relative_l2_warn:
        parser.error("continuity-relative-l2-fail must be >= continuity-relative-l2-warn")

    return args


def choose_tol(primary: float | None, secondary: float | None, default: float) -> float:
    if primary is not None:
        return primary
    if secondary is not None:
        return secondary
    return default


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


def candidate_paths(diagnostics_dir: Path, filename: str) -> list[Path]:
    paths = [diagnostics_dir / filename]
    parent = diagnostics_dir.parent
    if parent != diagnostics_dir:
        paths.append(parent / filename)
    unique: list[Path] = []
    seen: set[Path] = set()
    for path in paths:
        if path not in seen:
            unique.append(path)
            seen.add(path)
    return unique


def find_existing_csv(diagnostics_dir: Path, filename: str) -> Path | None:
    for path in candidate_paths(diagnostics_dir, filename):
        if path.is_file():
            return path
    return None


def parse_float(text: str) -> float:
    return float(text.strip().replace("D", "E").replace("d", "e"))


def read_last_row(path: Path) -> tuple[list[str], dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("CSV file has no header")
        last: dict[str, str] | None = None
        for row in reader:
            if row and any((value or "").strip() for value in row.values()):
                last = row
        if last is None:
            raise ValueError("CSV file has no data rows")
        return list(reader.fieldnames), last


def find_column(fieldnames: list[str], aliases: tuple[str, ...]) -> str | None:
    exact = {name: name for name in fieldnames}
    lower = {name.strip().lower(): name for name in fieldnames}
    for alias in aliases:
        if alias in exact:
            return exact[alias]
        key = alias.strip().lower()
        if key in lower:
            return lower[key]
    return None


def get_csv(
    diagnostics_dir: Path,
    filename: str,
    cache: dict[str, tuple[Path, list[str], dict[str, str]] | None],
) -> tuple[Path, list[str], dict[str, str]] | None:
    if filename not in cache:
        path = find_existing_csv(diagnostics_dir, filename)
        if path is None:
            cache[filename] = None
        else:
            cache[filename] = (path, *read_last_row(path))
    return cache[filename]


def missing_result(group: str, label: str, filename: str, status: str, message: str) -> MetricResult:
    return MetricResult(
        group=group,
        label=label,
        filename=filename,
        column=None,
        value=None,
        tolerance_text="n/a",
        status=status,
        message=message,
    )


def strict_metric(
    group: str,
    label: str,
    filename: str,
    aliases: tuple[str, ...],
    tolerance: float,
    diagnostics_dir: Path,
    cache: dict[str, tuple[Path, list[str], dict[str, str]] | None],
    optional: bool = False,
) -> MetricResult:
    try:
        cached = get_csv(diagnostics_dir, filename, cache)
    except Exception as exc:  # noqa: BLE001
        return missing_result(group, label, filename, "FAIL", f"could not read CSV: {exc}")

    if cached is None:
        status = "SKIP" if optional else "FAIL"
        message = "missing optional CSV file" if optional else "missing CSV file"
        return missing_result(group, label, filename, status, message)

    _path, fieldnames, row = cached
    column = find_column(fieldnames, aliases)
    if column is None:
        return missing_result(group, label, filename, "FAIL", "missing column: " + " or ".join(aliases))

    raw = row.get(column, "")
    try:
        value = parse_float(raw)
    except Exception as exc:  # noqa: BLE001
        return MetricResult(group, label, filename, column, None, f"<= {tolerance:.1e}", "FAIL", f"could not parse value {raw!r}: {exc}")

    if not math.isfinite(value):
        status = "FAIL"
        message = "non-finite value"
    elif abs(value) <= tolerance:
        status = "PASS"
        message = "within tolerance"
    else:
        status = "FAIL"
        message = "exceeds tolerance"

    return MetricResult(
        group=group,
        label=label,
        filename=filename,
        column=column,
        value=value,
        tolerance_text=f"<= {tolerance:.1e}",
        status=status,
        message=message,
    )


def graded_metric(
    group: str,
    label: str,
    filename: str,
    aliases: tuple[str, ...],
    pass_tolerance: float,
    fail_tolerance: float,
    diagnostics_dir: Path,
    cache: dict[str, tuple[Path, list[str], dict[str, str]] | None],
    optional: bool = False,
) -> MetricResult:
    try:
        cached = get_csv(diagnostics_dir, filename, cache)
    except Exception as exc:  # noqa: BLE001
        return missing_result(group, label, filename, "FAIL", f"could not read CSV: {exc}")

    if cached is None:
        status = "SKIP" if optional else "FAIL"
        message = "missing optional CSV file" if optional else "missing CSV file"
        return missing_result(group, label, filename, status, message)

    _path, fieldnames, row = cached
    column = find_column(fieldnames, aliases)
    if column is None:
        return missing_result(group, label, filename, "FAIL", "missing column: " + " or ".join(aliases))

    raw = row.get(column, "")
    try:
        value = parse_float(raw)
    except Exception as exc:  # noqa: BLE001
        return MetricResult(group, label, filename, column, None, f"PASS <= {pass_tolerance:.1e}; FAIL > {fail_tolerance:.1e}", "FAIL", f"could not parse value {raw!r}: {exc}")

    abs_value = abs(value)
    if not math.isfinite(value):
        status = "FAIL"
        message = "non-finite value"
    elif abs_value <= pass_tolerance:
        status = "PASS"
        message = "within pass threshold"
    elif abs_value <= fail_tolerance:
        status = "WARN"
        message = "above pass threshold but below fail threshold"
    else:
        status = "FAIL"
        message = "exceeds fail threshold"

    return MetricResult(
        group=group,
        label=label,
        filename=filename,
        column=column,
        value=value,
        tolerance_text=f"PASS <= {pass_tolerance:.1e}; WARN <= {fail_tolerance:.1e}; FAIL > {fail_tolerance:.1e}",
        status=status,
        message=message,
    )


def evaluate_results(
    diagnostics_dir: Path,
    args: argparse.Namespace,
    cache: dict[str, tuple[Path, list[str], dict[str, str]] | None],
) -> list[MetricResult]:
    results: list[MetricResult] = []

    for label, aliases in ENERGY_METRICS:
        results.append(
            strict_metric(
                "energy",
                label,
                ENERGY_FILE,
                aliases,
                args.energy_tol,
                diagnostics_dir,
                cache,
            )
        )

    if args.energy_only:
        return results

    optional_missing = bool(args.allow_missing_optional)

    if not args.skip_projection:
        for label, aliases, tolerance_attr in PROJECTION_METRICS:
            results.append(
                strict_metric(
                    "projection",
                    label,
                    PROJECTION_FILE,
                    aliases,
                    float(getattr(args, tolerance_attr)),
                    diagnostics_dir,
                    cache,
                    optional=optional_missing,
                )
            )

    if not args.skip_continuity:
        results.append(
            strict_metric(
                "continuity",
                "conservative continuity global integral",
                CONTINUITY_FILE,
                ("integral_drho_dt_plus_div_mass_flux_dV",),
                args.continuity_integral_abs_tol,
                diagnostics_dir,
                cache,
                optional=optional_missing,
            )
        )
        results.append(
            graded_metric(
                "continuity",
                "conservative continuity relative L2 residual (local quality)",
                CONTINUITY_FILE,
                ("relative_conservative_residual_l2",),
                args.continuity_relative_l2_warn,
                args.continuity_relative_l2_fail,
                diagnostics_dir,
                cache,
                optional=optional_missing,
            )
        )

    return results


def collect_context(
    diagnostics_dir: Path,
    cache: dict[str, tuple[Path, list[str], dict[str, str]] | None],
) -> dict[str, dict[str, Any]]:
    context: dict[str, dict[str, Any]] = {}
    for filename, columns in OPTIONAL_CONTEXT_COLUMNS.items():
        try:
            cached = get_csv(diagnostics_dir, filename, cache)
        except Exception:
            cached = None
        if cached is None:
            continue
        path, fieldnames, row = cached
        found: dict[str, Any] = {"path": str(path)}
        for column in columns:
            actual = find_column(fieldnames, (column,))
            if actual is None:
                continue
            raw = row.get(actual, "")
            try:
                found[actual] = parse_float(raw)
            except Exception:
                found[actual] = raw
        context[filename] = found
    return context


def result_to_dict(result: MetricResult) -> dict[str, Any]:
    return {
        "group": result.group,
        "label": result.label,
        "filename": result.filename,
        "column": result.column,
        "value": result.value,
        "tolerance": result.tolerance_text,
        "status": result.status,
        "message": result.message,
    }


def overall_status(results: list[MetricResult]) -> str:
    if any(result.status == "FAIL" for result in results):
        return "FAIL"
    if any(result.status == "WARN" for result in results):
        return "PASS_WITH_WARNINGS"
    return "PASS"


def print_text_report(
    diagnostics_dir: Path,
    results: list[MetricResult],
    context: dict[str, dict[str, Any]],
    args: argparse.Namespace,
) -> None:
    overall = overall_status(results)
    print("Variable-density validation summary")
    print("===================================")
    print(f"diagnostics_dir: {diagnostics_dir}")
    print(f"overall: {overall}")
    print()
    print("Tolerance policy:")
    print(f"  energy relative defects:                  <= {args.energy_tol:.1e}")
    if not args.energy_only and not args.skip_projection:
        print(f"  projection-source relative L2 residual:   <= {args.projection_l2_tol:.1e}")
        print(f"  projection-source relative max residual:  <= {args.projection_max_tol:.1e}")
    if not args.energy_only and not args.skip_continuity:
        print(f"  continuity global integral absolute:      <= {args.continuity_integral_abs_tol:.1e}")
        print(
            "  continuity local relative L2:             "
            f"PASS <= {args.continuity_relative_l2_warn:.1e}; "
            f"WARN <= {args.continuity_relative_l2_fail:.1e}; "
            f"FAIL > {args.continuity_relative_l2_fail:.1e}"
        )
    print()

    current_group: str | None = None
    for result in results:
        if result.group != current_group:
            current_group = result.group
            print(f"{current_group.capitalize()}:")
        if result.value is None:
            value_text = "n/a"
        else:
            value_text = f"{result.value:.8e}"
        column_text = result.column or "n/a"
        print(
            f"  [{result.status:4}] {result.label}\n"
            f"         column={column_text} value={value_text} tolerance={result.tolerance_text}\n"
            f"         {result.message}"
        )
    print()

    if context:
        print("Context from last diagnostic rows:")
        for filename, values in context.items():
            print(f"  {filename}:")
            for key, value in values.items():
                if key == "path":
                    continue
                if isinstance(value, float):
                    print(f"    {key}: {value:.8e}")
                else:
                    print(f"    {key}: {value}")
        print()

    print("Metric interpretation:")
    print("  Energy metrics are strict pass/fail checks of direct and reconciled rho*h closure.")
    print("  Projection metrics use projection-scale tolerances, not roundoff energy tolerances.")
    print("  Continuity global integral closure is pass/fail; local relative L2 is a quality")
    print("  metric that may WARN before it becomes a failure.")
    print("  Large unreconciled rho*h defects are not failures when rel_output_recon_defect")
    print("  and rel_operator_recon_defect pass.")


def main() -> int:
    args = parse_args()
    diagnostics_dir = resolve_diagnostics_dir(args)
    cache: dict[str, tuple[Path, list[str], dict[str, str]] | None] = {}

    results = evaluate_results(diagnostics_dir, args, cache)
    context = collect_context(diagnostics_dir, cache)
    status = overall_status(results)

    if args.json:
        payload = {
            "diagnostics_dir": str(diagnostics_dir),
            "overall": status,
            "strict_warnings": bool(args.strict_warnings),
            "results": [result_to_dict(result) for result in results],
            "context": context,
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print_text_report(diagnostics_dir, results, context, args)

    if status == "FAIL":
        return 1
    if status == "PASS_WITH_WARNINGS" and args.strict_warnings:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
