#!/usr/bin/env python3
'''
Patch 002b: fix missing energy arrays in VTU piece files.

Use this if ParaView warns that piece CellData arrays are missing for:
  - temperature
  - enthalpy
  - qrad

Symptom:
  .pvtu files contain PDataArray entries, but .vtu piece files do not contain
  matching DataArray entries.

This patch modifies only src/mod_output.f90 and inserts the energy DataArray
block inside write_vtu_unstructured, immediately before </CellData>.

Run from repository root:

  python tools/patches/002b_fix_vtu_energy_piece_arrays.py --dry-run
  python tools/patches/002b_fix_vtu_energy_piece_arrays.py --apply
'''

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


ENERGY_VTU_BLOCK = r'''
      if (params%enable_energy) then
         if (.not. allocated(energy%T) .or. .not. allocated(energy%h) .or. &
             .not. allocated(energy%qrad)) then
            call fatal_error('output', 'energy output requested but energy arrays are not allocated')
         end if

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="temperature" format="ascii">'
         do n = 1, n_owned
            c = owned_indices(n)
            write(unit_id,'(es24.16)') energy%T(c)
         end do
         write(unit_id,'(a)') '        </DataArray>'

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="enthalpy" format="ascii">'
         do n = 1, n_owned
            c = owned_indices(n)
            write(unit_id,'(es24.16)') energy%h(c)
         end do
         write(unit_id,'(a)') '        </DataArray>'

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="qrad" format="ascii">'
         do n = 1, n_owned
            c = owned_indices(n)
            write(unit_id,'(es24.16)') energy%qrad(c)
         end do
         write(unit_id,'(a)') '        </DataArray>'
      end if

'''


def parse_args():
    p = argparse.ArgumentParser(description="Patch 002b: fix VTU energy piece arrays.")
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

    start, end = find_subroutine_bounds(text, "write_vtu_unstructured")
    body = text[start:end]

    if "type(energy_fields_t), intent(in) :: energy" not in body:
        raise SystemExit("[ERROR] write_vtu_unstructured does not have energy argument. Apply Patch 002 first.")

    if 'Name="temperature"' in body and 'Name="enthalpy"' in body and 'Name="qrad"' in body:
        return PatchResult(path, False, "write_vtu_unstructured already contains energy DataArray entries")

    lines = body.splitlines(keepends=True)
    close_idx = None
    for i, line in enumerate(lines):
        if "</CellData>" in line:
            close_idx = i
            break

    if close_idx is None:
        raise SystemExit("[ERROR] could not find </CellData> write line inside write_vtu_unstructured")

    before_close = "".join(lines[:close_idx])
    if "<CellData" not in before_close:
        raise SystemExit("[ERROR] found </CellData> but no preceding <CellData> inside write_vtu_unstructured")

    new_body = "".join(lines[:close_idx]) + ENERGY_VTU_BLOCK + "".join(lines[close_idx:])
    text = text[:start] + new_body + text[end:]

    if text == original:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    return PatchResult(path, True, "would update")


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_mod_output(repo, apply)

    print("Patch 002b summary")
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
        print("  1. make clean && make -j")
        print("  2. remove old output directory or run into a fresh output directory")
        print("  3. run with enable_energy = .true.")
        print("  4. grep .vtu pieces for Name=\"temperature\", Name=\"enthalpy\", Name=\"qrad\"")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
