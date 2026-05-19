#!/usr/bin/env python3
"""Run MPI scaling sweeps and extract the solver profiler summary."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROFILE_RE = re.compile(
    r"^\s*(?P<name>[A-Za-z0-9_/]+(?:_[A-Za-z0-9_/]+)*)\s+"
    r"(?P<calls>\d+)\s+"
    r"(?P<min>[0-9.Ee+-]+)\s+"
    r"(?P<max>[0-9.Ee+-]+)\s+"
    r"(?P<avg>[0-9.Ee+-]+)"
    r"(?:\s+(?P<avg_pct>[0-9.Ee+-]+))?\s*$"
)
STEP_RE = re.compile(
    r"step=\s*(?P<step>\d+).*?"
    r"max_div=\s*(?P<max_div>[0-9.Ee+-]+).*?"
    r"piter=\s*(?P<piter>\d+).*?"
    r"\|U\|max=\s*(?P<umax>[0-9.Ee+-]+)"
)


def run_command(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def rewrite_output_dir(case_file: Path, output_dir: Path, dest: Path) -> None:
    text = case_file.read_text()
    replacement = f'  output_dir = "{output_dir.as_posix()}"'
    text = re.sub(r'^\s*output_dir\s*=.*$', replacement, text, flags=re.MULTILINE)
    dest.write_text(text)


def parse_run(stdout: str) -> dict[str, str]:
    row: dict[str, str] = {}
    for line in stdout.splitlines():
        step_match = STEP_RE.search(line)
        if step_match:
            row.update(
                final_step=step_match.group("step"),
                max_divergence=step_match.group("max_div"),
                pressure_iterations=step_match.group("piter"),
                max_velocity=step_match.group("umax"),
            )
            continue

        profile_match = PROFILE_RE.match(line)
        if not profile_match:
            continue

        key = profile_match.group("name").lower()
        row[f"{key}_calls"] = profile_match.group("calls")
        row[f"{key}_max_s"] = profile_match.group("max")
        row[f"{key}_avg_s"] = profile_match.group("avg")
        if profile_match.group("avg_pct") is not None:
            row[f"{key}_avg_pct"] = profile_match.group("avg_pct")

    return row


def case_path(case_name: str) -> Path:
    path = Path(case_name)
    if not path.is_absolute():
        path = ROOT / "cases" / case_name
    if path.is_dir():
        path = path / "case.nml"
    if not path.exists():
        raise FileNotFoundError(f"case file not found: {path}")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", default=["lid_driven_cavity"])
    parser.add_argument("--ranks", nargs="+", type=int, default=[1, 2, 4, 8])
    parser.add_argument("--build-target", default="debug")
    parser.add_argument("--output", default="benchmark_summary.csv")
    parser.add_argument("--keep-output", action="store_true")
    args = parser.parse_args()

    build = run_command(["make", args.build_target], ROOT)
    sys.stdout.write(build.stdout)
    if build.returncode != 0:
        return build.returncode

    rows: list[dict[str, str]] = []
    work_parent = ROOT / "tmp_benchmarks" if args.keep_output else None
    if work_parent is not None:
        work_parent.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory(dir=work_parent) as tmp:
        tmp_root = Path(tmp)
        for case_arg in args.cases:
            source_case = case_path(case_arg)
            case_label = source_case.parent.name

            for ranks in args.ranks:
                run_dir = tmp_root / f"{case_label}_np{ranks}"
                output_dir = run_dir / "output"
                run_dir.mkdir(parents=True, exist_ok=True)
                output_dir.mkdir(parents=True, exist_ok=True)
                tmp_case = run_dir / "case.nml"
                rewrite_output_dir(source_case, output_dir, tmp_case)

                cmd = ["mpirun", "-np", str(ranks), "./lowmach_react_hex", str(tmp_case)]
                print(f"running {case_label} np={ranks}", flush=True)
                result = run_command(cmd, ROOT)
                sys.stdout.write(result.stdout)

                row = {
                    "case": case_label,
                    "np": str(ranks),
                    "returncode": str(result.returncode),
                }
                row.update(parse_run(result.stdout))
                rows.append(row)

                if result.returncode != 0:
                    break

        if args.keep_output:
            print(f"kept benchmark outputs under {tmp_root}")
        elif work_parent is not None and work_parent.exists():
            shutil.rmtree(work_parent, ignore_errors=True)

    fieldnames = sorted({key for row in rows for key in row})
    output_csv = Path(args.output)
    if not output_csv.is_absolute():
        output_csv = ROOT / output_csv
    with output_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {output_csv}")
    return 0 if all(row["returncode"] == "0" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
