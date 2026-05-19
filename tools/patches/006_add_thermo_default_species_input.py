#!/usr/bin/env python3
"""
Patch 006: add thermo_default_species to &energy_input.

Purpose
-------
Stage 2A currently defaults the no-transported-species Cantera thermo mixture
to N2 when enable_cantera_thermo = .true. and enable_species = .false.

This patch adds an input flag:

    thermo_default_species = "N2"

so the user can switch the inert single-species thermo state to "H2O",
"CO2", "O2", etc., as long as that species exists in the Cantera mechanism.

Files modified:
  - src/mod_input.f90
  - src/mod_transport_properties.f90

Run from repository root:

    python tools/patches/006_add_thermo_default_species_input.py --dry-run
    python tools/patches/006_add_thermo_default_species_input.py --apply
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
    p = argparse.ArgumentParser(description="Patch 006: add thermo_default_species input.")
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


def read_text(path: Path) -> str:
    if not path.exists():
        raise SystemExit(f"[ERROR] required file not found: {path}")
    return path.read_text(encoding="utf-8")


def write_if_changed(path: Path, text: str, repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    old = read_text(path)
    if old == text:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")

    return PatchResult(path, True, "would update" if not apply else "updated")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if n > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {n} matches")
    return text.replace(old, new, 1)


def insert_after_once(text: str, anchor: str, insertion: str, label: str) -> str:
    n = text.count(anchor)
    if n == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if n > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {n} matches")
    return text.replace(anchor, anchor + insertion, 1)


def patch_mod_input(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_input.f90"
    text = read_text(path)
    changed = False

    # 1. Add field to case_params_t near energy controls.
    if "thermo_default_species" not in text:
        anchor = (
            "      logical :: enable_cantera_thermo = .false.                 !< Use Cantera for h(T,Y,p) and T(h,Y,p) in later patches.\n"
        )
        insertion = (
            "      character(len=name_len) :: thermo_default_species = 'N2'   !< Single-species Cantera thermo fallback when species transport is off.\n"
        )
        text = insert_after_once(text, anchor, insertion, "case_params_t thermo_default_species field")
        changed = True

    # 2. Add local variable in read_energy_input.
    if "character(len=name_len) :: thermo_default_species" not in text:
        anchor = (
            "      logical :: enable_energy\n"
            "      logical :: enable_cantera_thermo\n"
        )
        insertion = "      character(len=name_len) :: thermo_default_species\n"
        text = insert_after_once(text, anchor, insertion, "read_energy_input local thermo_default_species")
        changed = True

    # 3. Add namelist item.
    if "enable_energy, enable_cantera_thermo, thermo_default_species" not in text:
        old = (
            "      namelist /energy_input/ enable_energy, enable_cantera_thermo, &\n"
            "                               initial_T, energy_reference_T, energy_reference_h, &\n"
            "                               energy_cp, energy_lambda\n"
        )
        new = (
            "      namelist /energy_input/ enable_energy, enable_cantera_thermo, thermo_default_species, &\n"
            "                               initial_T, energy_reference_T, energy_reference_h, &\n"
            "                               energy_cp, energy_lambda\n"
        )
        text = replace_once(text, old, new, "read_energy_input namelist")
        changed = True

    # 4. Initialize local variable from params.
    if "thermo_default_species = params%thermo_default_species" not in text:
        anchor = (
            "      enable_energy = params%enable_energy\n"
            "      enable_cantera_thermo = params%enable_cantera_thermo\n"
        )
        insertion = "      thermo_default_species = params%thermo_default_species\n"
        text = insert_after_once(text, anchor, insertion, "read_energy_input local default assignment")
        changed = True

    # 5. Copy back into params after successful read.
    if "params%thermo_default_species = trim(thermo_default_species)" not in text:
        anchor = (
            "         params%enable_energy = enable_energy\n"
            "         params%enable_cantera_thermo = enable_cantera_thermo\n"
        )
        insertion = "         params%thermo_default_species = trim(thermo_default_species)\n"
        text = insert_after_once(text, anchor, insertion, "read_energy_input params assignment")
        changed = True

    # 6. Validation.
    if "thermo_default_species must not be empty" not in text:
        anchor = (
            "      if (params%energy_cp <= zero) call fatal_error('input', 'energy_cp must be positive')\n"
            "      if (params%energy_lambda < zero) call fatal_error('input', 'energy_lambda must be non-negative')\n"
        )
        insertion = (
            "      if (params%enable_cantera_thermo .and. len_trim(params%thermo_default_species) == 0) &\n"
            "         call fatal_error('input', 'thermo_default_species must not be empty when enable_cantera_thermo is true')\n"
        )
        text = insert_after_once(text, anchor, insertion, "thermo_default_species validation")
        changed = True

    if not changed:
        return PatchResult(path, False, "already patched")

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_mod_transport_properties(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_transport_properties.f90"
    text = read_text(path)
    changed = False

    # Replace the hard-coded Stage-2A N2 fallback with the input-selected species.
    old = (
        "      else if (params%enable_cantera_thermo .and. params%nspecies == 0) then\n"
        "         ! Stage-2A no-species thermo mode: use N2 as inert bath gas.\n"
        "         params%nspecies = 1\n"
        "         params%species_name(1) = 'N2'\n"
    )
    new = (
        "      else if (params%enable_cantera_thermo .and. params%nspecies == 0) then\n"
        "         ! Stage-2A no-species thermo mode: use user-selected inert/default species.\n"
        "         params%nspecies = 1\n"
        "         params%species_name(1) = trim(params%thermo_default_species)\n"
    )

    if old in text:
        text = text.replace(old, new, 1)
        changed = True
    elif "params%species_name(1) = trim(params%thermo_default_species)" in text:
        pass
    else:
        # Fallback for slightly different comments around the same block.
        block_start = "      else if (params%enable_cantera_thermo .and. params%nspecies == 0) then\n"
        idx = text.find(block_start)
        if idx < 0:
            raise SystemExit("[ERROR] enable_cantera_thermo no-species fallback block not found in mod_transport_properties.f90")

        next_block = text.find("      else if", idx + len(block_start))
        if next_block < 0:
            next_block = text.find("      end if", idx)
        if next_block < 0:
            raise SystemExit("[ERROR] could not find end of enable_cantera_thermo fallback block")

        block = text[idx:next_block]
        if "params%species_name(1)" not in block:
            raise SystemExit("[ERROR] thermo fallback block found but species assignment was not found")

        lines = block.splitlines(True)
        out = []
        replaced = False
        for line in lines:
            if "params%species_name(1)" in line and "=" in line:
                out.append("         params%species_name(1) = trim(params%thermo_default_species)\n")
                replaced = True
            else:
                out.append(line)
        if not replaced:
            raise SystemExit("[ERROR] failed to replace thermo fallback species assignment")

        text = text[:idx] + "".join(out) + text[next_block:]
        changed = True

    if not changed:
        return PatchResult(path, False, "already patched")

    return write_if_changed(path, text, repo, backup_root, apply)


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    backup_root = repo / ".backups" / timestamp()

    results = [
        patch_mod_input(repo, backup_root, apply),
        patch_mod_transport_properties(repo, backup_root, apply),
    ]

    print("Patch 006 summary")
    print("=================")
    print(f"repo: {repo}")
    print(f"mode: {'apply' if apply else 'dry-run'}")
    if apply:
        print(f"backup: {backup_root}")
    print()

    changed_any = False
    for r in results:
        changed_any = changed_any or r.changed
        try:
            rel = r.path.relative_to(repo)
        except ValueError:
            rel = r.path
        marker = "CHANGE" if r.changed else "SKIP"
        print(f"[{marker}] {rel}: {r.message}")

    if not apply:
        print("\nDry run only. Re-run with --apply to write changes.")
    elif changed_any:
        print("\nNext steps:")
        print("  1. make clean")
        print("  2. make release")
        print("  3. add thermo_default_species to &energy_input, for example thermo_default_species = \"CO2\"")
    else:
        print("\nNo changes were needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
