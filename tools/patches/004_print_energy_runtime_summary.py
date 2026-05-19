#!/usr/bin/env python3
"""
Patch 004: print compact energy diagnostics to the runtime terminal.

This patch does not change the CSV diagnostics format. It only adds one
rank-0 terminal line at the same cadence as write_energy_diagnostics_row(...):

  energy step=1000  Tmin=...  Tmax=...  Tmean=...  dTmax=...  rel_h=...

It reuses the global quantities already computed inside mod_energy, avoiding
any extra MPI reductions or new data paths.

Run from the repository root:

  python tools/patches/004_print_energy_runtime_summary.py --dry-run
  python tools/patches/004_print_energy_runtime_summary.py --apply

Backups of modified files are written to:

  .backups/YYYYMMDD_HHMMSS/
"""

from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PatchResult:
    path: Path
    changed: bool
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Patch 004: print runtime energy summary.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true", help="write changes to disk")
    mode.add_argument("--dry-run", action="store_true", help="show planned changes without writing")
    parser.add_argument("--repo", default=".", help="repository root; default: current directory")
    return parser.parse_args()


def timestamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] required file not found: {path}")


def make_backup(path: Path, repo: Path, backup_root: Path, apply: bool) -> None:
    if not apply:
        return
    rel = path.relative_to(repo)
    dst = backup_root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, dst)


def write_if_changed(path: Path, text: str, repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    old = read_text(path)
    if old == text:
        return PatchResult(path, False, "unchanged")

    if apply:
        make_backup(path, repo, backup_root, apply=True)
        path.write_text(text, encoding="utf-8")

    return PatchResult(path, True, "would update" if not apply else "updated")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if count > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {count} matches")
    return text.replace(old, new, 1)


def patch_mod_energy(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_energy.f90"
    text = read_text(path)

    if "energy step=" in text and "output_unit" in text:
        return PatchResult(path, False, "already patched")

    text = replace_once(
        text,
        "   use mod_kinds, only : rk, zero, tiny_safe, fatal_error\n",
        "   use mod_kinds, only : rk, zero, tiny_safe, fatal_error, output_unit\n",
        "mod_energy output_unit import",
    )

    old_block = """      write(unit_id,'(i0,12(\",\",es16.8))') step, time, global_min_T, global_max_T, mean_T, &
         global_min_h, global_max_h, mean_h, global_min_qrad, global_max_qrad, &
         global_integral_qrad, global_max_delta_T, rel_h_residual

      close(unit_id)
"""

    new_block = """      write(unit_id,'(i0,12(\",\",es16.8))') step, time, global_min_T, global_max_T, mean_T, &
         global_min_h, global_max_h, mean_h, global_min_qrad, global_max_qrad, &
         global_integral_qrad, global_max_delta_T, rel_h_residual

      close(unit_id)

      write(output_unit,'(a,i0,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5)') &
         'energy step=', step, 'Tmin=', global_min_T, 'Tmax=', global_max_T, &
         'Tmean=', mean_T, 'dTmax=', global_max_delta_T, 'rel_h=', rel_h_residual
"""

    text = replace_once(text, old_block, new_block, "energy CSV row terminal summary insertion")

    return write_if_changed(path, text, repo, backup_root, apply)


def sanity_check_repo(repo: Path) -> None:
    required = [repo / "src" / "mod_energy.f90"]
    missing = [p for p in required if not p.exists()]
    if missing:
        for p in missing:
            print(f"[ERROR] missing required file: {p}", file=sys.stderr)
        raise SystemExit(2)


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    sanity_check_repo(repo)

    backup_root = repo / ".backups" / timestamp()
    results = [patch_mod_energy(repo, backup_root, apply)]

    print("\nPatch 004 summary")
    print("=================")
    print(f"repo:   {repo}")
    print(f"mode:   {'apply' if apply else 'dry-run'}")
    if apply:
        print(f"backup: {backup_root}")
    print()

    changed_any = False
    for result in results:
        changed_any = changed_any or result.changed
        try:
            rel = result.path.relative_to(repo)
        except ValueError:
            rel = result.path
        marker = "CHANGE" if result.changed else "SKIP"
        print(f"[{marker}] {rel}: {result.message}")

    if not apply:
        print("\nDry run only. Re-run with --apply to write changes.")
    elif changed_any:
        print("\nApplied patch. Next steps:")
        print("  1. inspect git diff")
        print("  2. rebuild")
        print("  3. run rectangle_2D and confirm energy lines print after flow lines")
    else:
        print("\nNo changes were needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
