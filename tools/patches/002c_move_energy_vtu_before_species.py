#!/usr/bin/env python3
'''
Patch 002c: move energy VTU DataArrays outside the species-output branch.

Problem fixed:
  .pvtu declares temperature/enthalpy/qrad whenever enable_energy = true,
  but .vtu pieces only write those arrays when the energy block is accidentally
  nested inside:

      if (params%enable_species .and. params%nspecies > 0) then

This causes ParaView warnings when energy is enabled but species output is off.

This patch modifies only src/mod_output.f90:
  - finds the energy DataArray block in write_vtu_unstructured
  - removes it from inside the species branch if necessary
  - reinserts it before the species branch, still inside <CellData>

Run from repository root:

  python tools/patches/002c_move_energy_vtu_before_species.py --dry-run
  python tools/patches/002c_move_energy_vtu_before_species.py --apply
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


def parse_args():
    p = argparse.ArgumentParser(description="Patch 002c: move energy VTU arrays before species branch.")
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


def extract_energy_block(body: str) -> tuple[str, int, int]:
    start_marker = "      if (params%enable_energy) then\n         if (.not. allocated(energy%T)"
    start = body.find(start_marker)
    if start < 0:
        raise SystemExit("[ERROR] energy VTU DataArray block not found in write_vtu_unstructured")

    # The block ends at the matching first `      end if` after the qrad DataArray.
    qrad_pos = body.find('Name="qrad"', start)
    if qrad_pos < 0:
        raise SystemExit("[ERROR] qrad DataArray not found inside energy block")

    end_marker = "      end if\n\n"
    end = body.find(end_marker, qrad_pos)
    if end < 0:
        raise SystemExit("[ERROR] could not find end of energy VTU block")
    end += len(end_marker)

    return body[start:end], start, end


def patch_mod_output(repo: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_output.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    sub_start, sub_end = find_subroutine_bounds(text, "write_vtu_unstructured")
    body = text[sub_start:sub_end]

    species_if = "      if (params%enable_species .and. params%nspecies > 0) then\n"
    species_pos = body.find(species_if)
    if species_pos < 0:
        raise SystemExit("[ERROR] species output branch not found")

    energy_block, energy_start, energy_end = extract_energy_block(body)

    if energy_start < species_pos:
        return PatchResult(path, False, "energy VTU block is already before species branch")

    # Remove the block from its current position.
    body_without = body[:energy_start] + body[energy_end:]

    # Recompute species branch position after removal.
    species_pos2 = body_without.find(species_if)
    if species_pos2 < 0:
        raise SystemExit("[ERROR] species branch lost after block removal")

    # Insert energy block immediately before species branch.
    new_body = body_without[:species_pos2] + energy_block + body_without[species_pos2:]
    new_text = text[:sub_start] + new_body + text[sub_end:]

    if new_text == original:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(new_text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    return PatchResult(path, True, "would move energy VTU block before species branch")


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_mod_output(repo, apply)

    print("Patch 002c summary")
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
        print("  2. rm -rf cases/rectangle_2D/output")
        print("  3. rerun with enable_energy = .true.")
        print("  4. grep .vtu pieces for Name=\"temperature\", Name=\"enthalpy\", Name=\"qrad\"")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
