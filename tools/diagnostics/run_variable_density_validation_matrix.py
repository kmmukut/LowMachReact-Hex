#!/usr/bin/env python3
"""
Run accepted variable-density validation tools over a case matrix.

This wrapper delegates to:
  tools/diagnostics/check_variable_density_validation.py
  tools/diagnostics/audit_coupled_transport_conservation.py

It can audit existing outputs with --case/--output-dir, or read a CSV matrix.
CSV columns:
  name, output_dir, run_command, cwd, enabled, energy_only

run_command entries are skipped unless --run-commands is passed.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence


@dataclass
class Case:
    name: str
    output_dir: Path
    run_command: str = ""
    cwd: Optional[Path] = None
    enabled: bool = True
    energy_only: bool = False


@dataclass
class Result:
    name: str
    output_dir: Path
    overall: str
    run_status: str
    validation_status: str
    audit_status: str
    run_returncode: int
    validation_returncode: int
    audit_returncode: int
    notes: str = ""


def truthy(x: object) -> bool:
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y", "on"}


def falsey(x: object) -> bool:
    return str(x).strip().lower() in {"0", "false", "f", "no", "n", "off", ""}


def repo_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_path(text: str, base: Path) -> Path:
    p = Path(text).expanduser()
    if not p.is_absolute():
        p = (base / p).resolve()
    return p


def row_get(row: dict, *names: str, default: str = "") -> str:
    lowered = {str(k).strip().lower(): v for k, v in row.items()}
    for name in names:
        if name.lower() in lowered and lowered[name.lower()] is not None:
            return str(lowered[name.lower()]).strip()
    return default


def enabled_from(value: object) -> bool:
    if value is None:
        return True
    if falsey(value):
        return False
    return True


def load_matrix(path: Path) -> List[Case]:
    if not path.exists():
        raise FileNotFoundError(f"matrix file not found: {path}")

    rows: List[dict]
    if path.suffix.lower() == ".json":
        data = json.loads(path.read_text(encoding="utf-8"))
        rows = data.get("cases", data) if isinstance(data, dict) else data
        if not isinstance(rows, list):
            raise ValueError("JSON matrix must be a list or an object with a cases list")
    else:
        with path.open("r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))

    base = path.parent.resolve()
    out: List[Case] = []
    for i, row in enumerate(rows, 1):
        enabled = enabled_from(row_get(row, "enabled", default="true"))
        output_dir = row_get(row, "output_dir", "output", "case_output_dir", default="")
        if not output_dir:
            if enabled:
                raise ValueError(f"matrix row {i} is missing output_dir")
            output_dir = "."
        name = row_get(row, "name", "case", default=f"case_{i:03d}")
        cwd_text = row_get(row, "cwd", "workdir", "working_dir", default="")
        out.append(
            Case(
                name=name,
                output_dir=resolve_path(output_dir, base),
                run_command=row_get(row, "run_command", "command", default=""),
                cwd=resolve_path(cwd_text, base) if cwd_text else None,
                enabled=enabled,
                energy_only=truthy(row_get(row, "energy_only", default="false")),
            )
        )
    return out


def parse_case(text: str, base: Path) -> Case:
    if "=" in text:
        name, output = text.split("=", 1)
        name = name.strip()
        output = output.strip()
    else:
        output = text.strip()
        p = Path(output)
        name = p.parent.name if p.name == "output" and p.parent.name else (p.name or "case")
    return Case(name=name or "case", output_dir=resolve_path(output, base))


def parse_status(stdout: str, returncode: int) -> str:
    for line in stdout.splitlines():
        m = re.match(r"^\s*overall\s*:\s*([A-Za-z]+)\s*$", line)
        if m:
            status = m.group(1).upper()
            if status in {"PASS", "WARN", "FAIL", "SKIP"}:
                return status
    return "PASS" if returncode == 0 else "FAIL"


def run_cmd(cmd, cwd=None, shell=False, echo=False):
    if echo:
        if shell:
            print(f"$ {cmd}")
        else:
            print("$ " + " ".join(shlex.quote(str(x)) for x in cmd))
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        shell=shell,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def combine(*statuses: str) -> str:
    sts = [s for s in statuses if s and s != "SKIP"]
    if any(s == "FAIL" for s in sts):
        return "FAIL"
    if any(s == "WARN" for s in sts):
        return "WARN"
    return "PASS" if sts else "SKIP"


def validate_case(case: Case, args, repo_root: Path, checker: Path, audit_tool: Path) -> Result:
    if not case.enabled:
        return Result(case.name, case.output_dir, "SKIP", "SKIP", "SKIP", "SKIP", 0, 0, 0, "disabled")

    notes = []
    run_status = "PASS"
    run_rc = 0

    if case.run_command:
        if args.run_commands:
            r = run_cmd(case.run_command, cwd=case.cwd or repo_root, shell=True, echo=args.echo_commands)
            run_rc = r.returncode
            run_status = "PASS" if r.returncode == 0 else "FAIL"
            if r.returncode != 0:
                notes.append("run_command failed")
        else:
            run_status = "SKIP"
            notes.append("run_command skipped")

    validation_status = "FAIL"
    validation_rc = 1
    if checker.exists():
        cmd = [args.python, str(checker), "--output-dir", str(case.output_dir)]
        for extra in args.checker_arg:
            cmd.append(extra)
        if case.energy_only:
            cmd.append("--energy-only")
        r = run_cmd(cmd, cwd=repo_root, echo=args.echo_commands)
        validation_rc = r.returncode
        validation_status = parse_status(r.stdout, r.returncode)
        if r.returncode != 0 and validation_status == "PASS":
            validation_status = "FAIL"
    else:
        notes.append(f"missing checker: {checker}")

    audit_status = "SKIP"
    audit_rc = 0
    if not args.skip_audit:
        if audit_tool.exists():
            cmd = [args.python, str(audit_tool), "--output-dir", str(case.output_dir)]
            for extra in args.audit_arg:
                cmd.append(extra)
            if args.write_audit_csv:
                cmd.append("--write-csv")
            r = run_cmd(cmd, cwd=repo_root, echo=args.echo_commands)
            audit_rc = r.returncode
            audit_status = parse_status(r.stdout, r.returncode)
            if r.returncode != 0 and audit_status == "PASS":
                audit_status = "FAIL"
        else:
            audit_status = "FAIL"
            audit_rc = 1
            notes.append(f"missing audit tool: {audit_tool}")

    overall = combine(run_status, validation_status, audit_status)
    return Result(
        case.name,
        case.output_dir,
        overall,
        run_status,
        validation_status,
        audit_status,
        run_rc,
        validation_rc,
        audit_rc,
        "; ".join(notes),
    )


def write_summary(path: Path, results: Sequence[Result]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "name",
        "output_dir",
        "overall_status",
        "run_status",
        "validation_status",
        "audit_status",
        "run_returncode",
        "validation_returncode",
        "audit_returncode",
        "notes",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow(
                {
                    "name": r.name,
                    "output_dir": str(r.output_dir),
                    "overall_status": r.overall,
                    "run_status": r.run_status,
                    "validation_status": r.validation_status,
                    "audit_status": r.audit_status,
                    "run_returncode": r.run_returncode,
                    "validation_returncode": r.validation_returncode,
                    "audit_returncode": r.audit_returncode,
                    "notes": r.notes,
                }
            )


def print_summary(results: Sequence[Result], summary: Path) -> None:
    print("Variable-density validation matrix")
    print("==================================")
    print(f"cases: {len(results)}")
    print(f"summary_csv: {summary}")
    print()
    print("case  overall  run  validation  audit  notes")
    print("----  -------  ---  ----------  -----  -----")
    for r in results:
        print(f"{r.name}  {r.overall}  {r.run_status}  {r.validation_status}  {r.audit_status}  {r.notes}")
    print()
    print("Status policy:")
    print("  FAIL if any required validation/audit fails.")
    print("  WARN if no failures occur but any validation/audit reports WARN.")
    print("  PASS only when all enabled checks pass.")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Run variable-density validation across a case matrix.")
    p.add_argument("--matrix", type=Path, help="CSV or JSON matrix.")
    p.add_argument("--case", action="append", default=[], help="NAME=OUTPUT_DIR or OUTPUT_DIR. Repeatable.")
    p.add_argument("--output-dir", action="append", default=[], help="Existing output directory. Repeatable.")
    p.add_argument("--repo-root", type=Path, default=None)
    p.add_argument("--summary-csv", type=Path, default=Path("validation_matrix_summary.csv"))
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--checker", type=Path, default=None)
    p.add_argument("--audit-tool", type=Path, default=None)
    p.add_argument("--run-commands", action="store_true")
    p.add_argument("--skip-audit", action="store_true")
    p.add_argument("--write-audit-csv", action="store_true")
    p.add_argument("--checker-arg", action="append", default=[])
    p.add_argument("--audit-arg", action="append", default=[])
    p.add_argument("--fail-on-warn", action="store_true")
    p.add_argument("--echo-commands", action="store_true")
    args = p.parse_args(argv)

    repo_root = (args.repo_root or repo_root_from_script()).resolve()
    checker = (args.checker or (repo_root / "tools" / "diagnostics" / "check_variable_density_validation.py")).resolve()
    audit_tool = (args.audit_tool or (repo_root / "tools" / "diagnostics" / "audit_coupled_transport_conservation.py")).resolve()

    cases: List[Case] = []
    if args.matrix:
        cases.extend(load_matrix(args.matrix.resolve()))
    for text in args.case:
        cases.append(parse_case(text, Path.cwd()))
    for text in args.output_dir:
        cases.append(parse_case(text, Path.cwd()))
    if not cases:
        p.error("provide --matrix, --case, or --output-dir")

    results = [validate_case(c, args, repo_root, checker, audit_tool) for c in cases]
    summary = args.summary_csv.expanduser()
    if not summary.is_absolute():
        summary = (Path.cwd() / summary).resolve()
    write_summary(summary, results)
    print_summary(results, summary)

    has_fail = any(r.overall == "FAIL" for r in results)
    has_warn = any(r.overall == "WARN" for r in results)
    return 1 if has_fail or (args.fail_on_warn and has_warn) else 0


if __name__ == "__main__":
    raise SystemExit(main())
