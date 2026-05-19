#!/usr/bin/env python3
"""
Patch 007b: export Stage-2A thermodynamic properties to VTU/PVTU.

This is a corrected version of Patch 007.

Problem in Patch 007:
  It checked for energy%cp / energy%lambda / energy%rho_thermo inside
  src/mod_output.f90 before inserting the output code. That check is wrong
  because those references do not exist in mod_output.f90 until this patch
  adds them.

Fix:
  Check src/mod_energy.f90 for the Stage-2A fields, then patch
  src/mod_output.f90.

Arrays exported:
  - cp
  - thermal_conductivity
  - rho_thermo

Run from repository root:

    python tools/patches/007b_export_thermo_properties_vtu.py --dry-run
    python tools/patches/007b_export_thermo_properties_vtu.py --apply
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


THERMO_VTU_BLOCK = r"""
         if (allocated(energy%cp)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="cp" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(es24.16)') energy%cp(c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

         if (allocated(energy%lambda)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="thermal_conductivity" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(es24.16)') energy%lambda(c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

         if (allocated(energy%rho_thermo)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="rho_thermo" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(es24.16)') energy%rho_thermo(c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

"""

THERMO_PVTU_BLOCK = r"""         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="cp" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="thermal_conductivity" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_thermo" format="ascii"/>'
"""


def parse_args():
    p = argparse.ArgumentParser(description="Patch 007b: export thermo properties to VTU.")
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


def check_energy_fields(repo: Path) -> None:
    path = repo / "src" / "mod_energy.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")

    missing = []
    for field in ("cp", "lambda", "rho_thermo"):
        if f"{field}(:)" not in text and f"energy%{field}" not in text:
            missing.append(field)

    if missing:
        raise SystemExit(
            "[ERROR] src/mod_energy.f90 does not contain Stage-2A thermo fields: "
            + ", ".join(missing)
            + ". Apply/fix the Stage-2A thermo patch first."
        )


def patch_vtu_piece(text: str) -> tuple[str, bool]:
    sub_start, sub_end = find_subroutine_bounds(text, "write_vtu_unstructured")
    body = text[sub_start:sub_end]

    if 'Name="thermal_conductivity"' in body and 'Name="rho_thermo"' in body and 'Name="cp"' in body:
        return text, False

    if "type(energy_fields_t), intent(in) :: energy" not in body:
        raise SystemExit("[ERROR] write_vtu_unstructured does not take energy argument. Apply energy output patches first.")

    energy_start = body.find("      if (params%enable_energy) then")
    if energy_start < 0:
        raise SystemExit("[ERROR] energy VTU block not found in write_vtu_unstructured")

    qrad_pos = body.find('Name="qrad"', energy_start)
    if qrad_pos < 0:
        raise SystemExit("[ERROR] qrad DataArray not found inside energy VTU block")

    close_marker = "         write(unit_id,'(a)') '        </DataArray>'\n"
    close_pos = body.find(close_marker, qrad_pos)
    if close_pos < 0:
        raise SystemExit("[ERROR] qrad closing DataArray write not found")

    insert_pos = close_pos + len(close_marker)
    new_body = body[:insert_pos] + THERMO_VTU_BLOCK + body[insert_pos:]
    new_text = text[:sub_start] + new_body + text[sub_end:]

    return new_text, True


def patch_pvtu(text: str) -> tuple[str, bool]:
    if 'Name="thermal_conductivity"' in text and 'Name="rho_thermo"' in text and 'Name="cp"' in text:
        return text, False

    anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"qrad\" format=\"ascii\"/>'\n"
    if anchor not in text:
        raise SystemExit("[ERROR] PVTU qrad PDataArray anchor not found")

    return text.replace(anchor, anchor + THERMO_PVTU_BLOCK, 1), True


def patch_mod_output(repo: Path, apply: bool) -> PatchResult:
    check_energy_fields(repo)

    path = repo / "src" / "mod_output.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    text, changed_vtu = patch_vtu_piece(text)
    text, changed_pvtu = patch_pvtu(text)

    if text == original:
        return PatchResult(path, False, "already patched")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    changes = []
    if changed_vtu:
        changes.append("VTU piece arrays")
    if changed_pvtu:
        changes.append("PVTU declarations")

    return PatchResult(path, True, "would update " + ", ".join(changes))


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_mod_output(repo, apply)

    print("Patch 007b summary")
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
        print("  3. remove old output or use a fresh output directory")
        print("  4. run with enable_energy = .true. and enable_cantera_thermo = .true.")
        print("  5. inspect VTU arrays: cp, thermal_conductivity, rho_thermo")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
