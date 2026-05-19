#!/usr/bin/env python3
"""
Patch 007d: repair misplaced thermo-default-species runtime print.

Problem:
  Patch 007c inserted the block

      if (params%enable_cantera_thermo ... ) then
         write(output_unit,'(a,a)') ...
      end if

  inside the multi-line energy write(output_unit,...) statement, causing
  Fortran syntax errors.

Fix:
  - Remove the misplaced block wherever it appears in write_energy_diagnostics_row.
  - Reinsert it immediately before the full energy runtime write statement.
  - Preserve the saved one-time-print flag declaration.

Run from repository root:

    python tools/patches/007d_repair_thermo_species_print.py --dry-run
    python tools/patches/007d_repair_thermo_species_print.py --apply
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PatchResult:
    path: Path
    changed: bool
    message: str


PRINT_BLOCK = (
    "      if (params%enable_cantera_thermo .and. .not. params%enable_species .and. &\n"
    "          step > 0 .and. .not. printed_thermo_default_species) then\n"
    "         write(output_unit,'(a,a)') 'Cantera thermo default species: ', trim(params%thermo_default_species)\n"
    "         printed_thermo_default_species = .true.\n"
    "      end if\n"
    "\n"
)


def parse_args():
    p = argparse.ArgumentParser(description="Patch 007d: repair runtime thermo species print placement.")
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


def ensure_output_unit_import(text: str) -> str:
    use_prefix = "   use mod_kinds, only : "
    idx = text.find(use_prefix)
    if idx < 0:
        return text

    line_end = text.find("\n", idx)
    if line_end < 0:
        return text

    line = text[idx:line_end]
    if "output_unit" in line:
        return text

    return text[:idx] + line + ", output_unit" + text[line_end:]


def ensure_saved_flag_declaration(body: str) -> str:
    if "logical, save :: printed_thermo_default_species" in body:
        return body

    anchors = [
        "      real(rk) :: mean_T, mean_h, rel_h_residual\n",
        "      real(rk) :: mean_T, mean_h\n",
    ]
    for anchor in anchors:
        if anchor in body:
            return body.replace(anchor, anchor + "      logical, save :: printed_thermo_default_species = .false.\n", 1)

    raise SystemExit("[ERROR] could not find declaration anchor in write_energy_diagnostics_row")


def remove_existing_print_blocks(body: str) -> str:
    # Remove the exact block if present.
    body = body.replace(PRINT_BLOCK, "")

    # Remove any slightly malformed/spacing-altered version of the same block.
    pattern = re.compile(
        r"\n?      if \(params%enable_cantera_thermo \.and\. \.not\. params%enable_species \.and\. &\n"
        r"          step > 0 \.and\. \.not\. printed_thermo_default_species\) then\n"
        r"         write\(output_unit,'\(a,a\)'\) 'Cantera thermo default species: ', trim\(params%thermo_default_species\)\n"
        r"         printed_thermo_default_species = \.true\.\n"
        r"      end if\n\n?",
        re.MULTILINE,
    )
    body = pattern.sub("\n", body)
    return body


def find_energy_write_statement_start(body: str) -> int:
    energy_pos = body.find("'energy step='")
    if energy_pos < 0:
        raise SystemExit("[ERROR] energy runtime print line containing 'energy step=' not found")

    # Walk backward line by line until the first line of the write statement.
    pos = energy_pos
    while True:
        line_start = body.rfind("\n", 0, pos)
        if line_start < 0:
            line_start = 0
        else:
            line_start += 1

        line = body[line_start:body.find("\n", line_start) if body.find("\n", line_start) >= 0 else len(body)]

        if "write(output_unit" in line:
            return line_start

        if line_start == 0:
            break

        pos = line_start - 1

    raise SystemExit("[ERROR] could not find beginning of energy write(output_unit,...) statement")


def patch_mod_energy(repo: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_energy.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    text = ensure_output_unit_import(text)

    sub_start, sub_end = find_subroutine_bounds(text, "write_energy_diagnostics_row")
    body = text[sub_start:sub_end]

    body = ensure_saved_flag_declaration(body)
    body = remove_existing_print_blocks(body)

    write_start = find_energy_write_statement_start(body)
    body = body[:write_start] + PRINT_BLOCK + body[write_start:]

    text = text[:sub_start] + body + text[sub_end:]

    if text == original:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    return PatchResult(path, True, "would move thermo default species print before energy write statement")


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_mod_energy(repo, apply)

    print("Patch 007d summary")
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
