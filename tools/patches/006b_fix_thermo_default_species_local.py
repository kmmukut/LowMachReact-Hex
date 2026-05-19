#!/usr/bin/env python3
"""
Patch 006b: fix missing local thermo_default_species declaration.

Problem:
  Patch 006 added thermo_default_species to case_params_t and to the
  /energy_input/ namelist, but in some code states it skipped the local
  declaration inside read_energy_input because the type-level declaration
  already contained the same substring.

Compiler error:
  Symbol 'thermo_default_species' has no IMPLICIT type

Fix:
  Insert:
      character(len=name_len) :: thermo_default_species
  inside subroutine read_energy_input, before the namelist statement.

Run from repository root:

    python tools/patches/006b_fix_thermo_default_species_local.py --dry-run
    python tools/patches/006b_fix_thermo_default_species_local.py --apply
"""

from __future__ import annotations

import argparse
import datetime as dt
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PatchResult:
    path: Path
    changed: bool
    message: str


def parse_args():
    p = argparse.ArgumentParser(description="Patch 006b: fix thermo_default_species local declaration.")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--apply", action="store_true")
    g.add_argument("--dry-run", action="store_true")
    p.add_argument("--repo", default=".")
    return p.parse_args()


def timestamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def backup(path: Path, repo: Path, backup_root: Path) -> Path:
    dst = backup_root / path.relative_to(repo)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, dst)
    return dst


def find_subroutine_bounds(text: str, name: str) -> tuple[int, int]:
    start = text.find(f"subroutine {name}")
    if start < 0:
        raise SystemExit(f"[ERROR] subroutine {name} not found")

    end_marker = f"end subroutine {name}"
    end = text.find(end_marker, start)
    if end < 0:
        raise SystemExit(f"[ERROR] {end_marker} not found")

    line_end = text.find("\n", end)
    if line_end < 0:
        line_end = len(text)
    else:
        line_end += 1

    return start, line_end


def patch_mod_input(repo: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_input.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    sub_start, sub_end = find_subroutine_bounds(text, "read_energy_input")
    body = text[sub_start:sub_end]

    if "thermo_default_species" not in body:
        raise SystemExit("[ERROR] read_energy_input does not contain thermo_default_species in the namelist. Apply Patch 006 first.")

    if "character(len=name_len) :: thermo_default_species" in body:
        return PatchResult(path, False, "local declaration already present")

    anchor_options = [
        "      logical :: enable_energy\n      logical :: enable_cantera_thermo\n",
        "      logical :: enable_cantera_thermo\n",
    ]

    insertion_done = False
    for anchor in anchor_options:
        if anchor in body:
            body = body.replace(anchor, anchor + "      character(len=name_len) :: thermo_default_species\n", 1)
            insertion_done = True
            break

    if not insertion_done:
        # Fallback: insert immediately before the namelist statement.
        namelist_anchor = "      namelist /energy_input/"
        idx = body.find(namelist_anchor)
        if idx < 0:
            raise SystemExit("[ERROR] read_energy_input namelist statement not found")
        body = body[:idx] + "      character(len=name_len) :: thermo_default_species\n" + body[idx:]

    text = text[:sub_start] + body + text[sub_end:]

    if text == original:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    return PatchResult(path, True, "would insert local declaration in read_energy_input")


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_mod_input(repo, apply)

    print("Patch 006b summary")
    print("==================")
    print(f"repo: {repo}")
    print(f"mode: {'apply' if apply else 'dry-run'}")

    marker = "CHANGE" if result.changed else "SKIP"
    try:
        rel = result.path.relative_to(repo)
    except ValueError:
        rel = result.path
    print(f"[{marker}] {rel}: {result.message}")

    if not apply and result.changed:
        print("\nDry run only. Re-run with --apply to write changes.")
    elif apply and result.changed:
        print("\nNext steps:")
        print("  1. make clean")
        print("  2. make release")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
