#!/usr/bin/env python3
"""
Patch 007c: fix thermo-property PVTU metadata and print thermo default species once.

Fixes:
  1. ParaView may not show cp / thermal_conductivity / rho_thermo if Patch 007b
     inserted the VTU piece arrays but skipped the PVTU PDataArray declarations.
     This patch specifically inserts the PVTU declarations after qrad.

  2. Adds a one-time runtime terminal message, printed immediately before the
     first nonzero energy summary line:

       Cantera thermo default species: N2

     It is printed only when:
       enable_energy = .true.
       enable_cantera_thermo = .true.
       enable_species = .false.
       step > 0

Files modified:
  - src/mod_output.f90
  - src/mod_energy.f90

Run from repository root:

    python tools/patches/007c_fix_pvtu_and_print_thermo_species.py --dry-run
    python tools/patches/007c_fix_pvtu_and_print_thermo_species.py --apply
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


PVTU_THERMO_LINES = (
    "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"cp\" format=\"ascii\"/>'\n"
    "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"thermal_conductivity\" format=\"ascii\"/>'\n"
    "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"rho_thermo\" format=\"ascii\"/>'\n"
)


def parse_args():
    p = argparse.ArgumentParser(description="Patch 007c: fix PVTU thermo metadata and print thermo species.")
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


def patch_mod_output(repo: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_output.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    qrad_anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"qrad\" format=\"ascii\"/>'\n"
    qrad_idx = text.find(qrad_anchor)
    if qrad_idx < 0:
        raise SystemExit("[ERROR] PVTU qrad PDataArray anchor not found in src/mod_output.f90")

    # Only consider the local PVTU metadata region after qrad. Patch 007b's bug
    # was using a global string check, which could be satisfied by the VTU piece
    # block while PVTU declarations were still missing.
    local_region = text[qrad_idx:qrad_idx + 1200]
    if 'Name="thermal_conductivity"' in local_region and 'Name="rho_thermo"' in local_region and 'Name="cp"' in local_region:
        return PatchResult(path, False, "PVTU thermo metadata already present")

    text = text[:qrad_idx + len(qrad_anchor)] + PVTU_THERMO_LINES + text[qrad_idx + len(qrad_anchor):]

    if text == original:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    return PatchResult(path, True, "would insert PVTU metadata for cp, thermal_conductivity, rho_thermo")


def ensure_output_unit_import(text: str) -> tuple[str, bool]:
    # The energy runtime print from previous patches usually already imports
    # output_unit. This keeps the new print safe if the current file does not.
    use_prefix = "   use mod_kinds, only : "
    idx = text.find(use_prefix)
    if idx < 0:
        return text, False

    line_end = text.find("\n", idx)
    if line_end < 0:
        return text, False

    line = text[idx:line_end]
    if "output_unit" in line:
        return text, False

    new_line = line + ", output_unit"
    return text[:idx] + new_line + text[line_end:], True


def patch_mod_energy(repo: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_energy.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    text, _ = ensure_output_unit_import(text)

    sub_start, sub_end = find_subroutine_bounds(text, "write_energy_diagnostics_row")
    body = text[sub_start:sub_end]

    if "Cantera thermo default species:" in body:
        if text == original:
            return PatchResult(path, False, "runtime thermo species print already present")
        if apply:
            backup_root = repo / ".backups" / timestamp()
            dst = backup(path, repo, backup_root)
            path.write_text(text, encoding="utf-8")
            return PatchResult(path, True, f"updated import only; backup: {dst}")
        return PatchResult(path, True, "would update import only")

    # Add a saved one-time-print flag in the declaration section.
    decl_inserted = False
    decl_anchors = [
        "      real(rk) :: mean_T, mean_h, rel_h_residual\n",
        "      real(rk) :: mean_T, mean_h\n",
    ]

    for anchor in decl_anchors:
        if anchor in body:
            body = body.replace(anchor, anchor + "      logical, save :: printed_thermo_default_species = .false.\n", 1)
            decl_inserted = True
            break

    if not decl_inserted:
        raise SystemExit("[ERROR] could not find declaration anchor in write_energy_diagnostics_row")

    # Insert immediately before the existing energy runtime summary line.
    runtime_line_idx = body.find("write(output_unit")
    energy_line_idx = body.find("'energy step='")
    if energy_line_idx < 0:
        raise SystemExit("[ERROR] energy runtime print line not found in write_energy_diagnostics_row")

    # Find beginning of the write statement containing 'energy step='.
    line_start = body.rfind("\n", 0, energy_line_idx)
    if line_start < 0:
        line_start = 0
    else:
        line_start += 1

    print_block = (
        "      if (params%enable_cantera_thermo .and. .not. params%enable_species .and. &\n"
        "          step > 0 .and. .not. printed_thermo_default_species) then\n"
        "         write(output_unit,'(a,a)') 'Cantera thermo default species: ', trim(params%thermo_default_species)\n"
        "         printed_thermo_default_species = .true.\n"
        "      end if\n"
        "\n"
    )

    body = body[:line_start] + print_block + body[line_start:]
    text = text[:sub_start] + body + text[sub_end:]

    if text == original:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    return PatchResult(path, True, "would add one-time Cantera thermo default species runtime print")


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()

    results = [
        patch_mod_output(repo, apply),
        patch_mod_energy(repo, apply),
    ]

    print("Patch 007c summary")
    print("==================")
    print(f"repo: {repo}")
    print(f"mode: {'apply' if apply else 'dry-run'}")
    print()

    changed_any = False
    for result in results:
        changed_any = changed_any or result.changed
        marker = "CHANGE" if result.changed else "SKIP"
        try:
            rel = result.path.relative_to(repo)
        except ValueError:
            rel = result.path
        print(f"[{marker}] {rel}: {result.message}")

    if not apply and changed_any:
        print("\nDry run only. Re-run with --apply to write changes.")
    elif apply and changed_any:
        print("\nNext steps:")
        print("  1. make clean")
        print("  2. make release")
        print("  3. delete old VTU/PVTU output or use a fresh output directory")
        print("  4. rerun and open the newly generated .pvd/.pvtu in ParaView")
    elif not changed_any:
        print("\nNo changes were needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
