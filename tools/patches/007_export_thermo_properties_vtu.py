#!/usr/bin/env python3
"""
Patch 007: export Stage-2A thermodynamic properties to VTU/PVTU.

Adds optional cell-data arrays to the VTU piece files and PVTU collection files:

  - cp
  - thermal_conductivity
  - rho_thermo

These arrays are written when energy is enabled and the corresponding arrays
exist in energy_fields_t. In Stage 2A they are populated from Cantera when
enable_cantera_thermo = .true.; in fallback constant-cp mode they contain the
configured constant values.

Files modified:
  - src/mod_output.f90

Run from repository root:

    python tools/patches/007_export_thermo_properties_vtu.py --dry-run
    python tools/patches/007_export_thermo_properties_vtu.py --apply
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
    p = argparse.ArgumentParser(description="Patch 007: export thermo properties to VTU.")
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


def find_energy_vtu_block_bounds(body: str) -> tuple[int, int]:
    start = body.find("      if (params%enable_energy) then\n         if (.not. allocated(energy%T)")
    if start < 0:
        raise SystemExit("[ERROR] energy VTU block not found in write_vtu_unstructured")

    qrad_pos = body.find('Name="qrad"', start)
    if qrad_pos < 0:
        raise SystemExit("[ERROR] qrad DataArray not found inside energy VTU block")

    end_marker = "      end if\n"
    end = body.find(end_marker, qrad_pos)
    if end < 0:
        raise SystemExit("[ERROR] could not find end of energy VTU block")
    end += len(end_marker)

    if body[end:end+1] == "\n":
        end += 1

    return start, end


def patch_vtu_piece(text: str) -> tuple[str, bool]:
    sub_start, sub_end = find_subroutine_bounds(text, "write_vtu_unstructured")
    body = text[sub_start:sub_end]

    if 'Name="thermal_conductivity"' in body and 'Name="rho_thermo"' in body and 'Name="cp"' in body:
        return text, False

    if "energy%cp" not in text or "energy%lambda" not in text or "energy%rho_thermo" not in text:
        raise SystemExit("[ERROR] energy_fields_t does not appear to contain cp/lambda/rho_thermo. Apply Stage-2A patch first.")

    block_start, block_end = find_energy_vtu_block_bounds(body)
    block = body[block_start:block_end]

    qrad_name_pos = block.find('Name="qrad"')
    if qrad_name_pos < 0:
        raise SystemExit("[ERROR] qrad DataArray not found while patching VTU energy block")

    close_marker = "         write(unit_id,'(a)') '        </DataArray>'\n"
    close_pos = block.find(close_marker, qrad_name_pos)
    if close_pos < 0:
        raise SystemExit("[ERROR] qrad closing DataArray write not found")
    insert_pos = close_pos + len(close_marker)

    new_block = block[:insert_pos] + THERMO_VTU_BLOCK + block[insert_pos:]
    new_body = body[:block_start] + new_block + body[block_end:]
    new_text = text[:sub_start] + new_body + text[sub_end:]
    return new_text, True


def patch_pvtu(text: str) -> tuple[str, bool]:
    if 'Name="thermal_conductivity"' in text and 'Name="rho_thermo"' in text and 'Name="cp"' in text:
        return text, False

    anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"qrad\" format=\"ascii\"/>'\n"
    if anchor not in text:
        raise SystemExit("[ERROR] PVTU qrad PDataArray anchor not found")

    text = text.replace(anchor, anchor + THERMO_PVTU_BLOCK, 1)
    return text, True


def patch_mod_output(repo: Path, apply: bool) -> PatchResult:
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

    print("Patch 007 summary")
    print("=================")
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
