#!/usr/bin/env python3
"""
Patch 005d: fix Cantera HPY inversion for installed Cantera C++ API.

Patch 005c used:

    gas->setState_HPY(h, p, Y)

but the Cantera version in the user's environment exposes ThermoPhase::setState_HP
and composition setters, not setState_HPY. This patch changes the recovery path to:

    gas->setMassFractions(Y)
    gas->setState_HP(h, p)

Run from repository root after Patch 005c:

    python tools/patches/005d_fix_cantera_hpy_api.py --dry-run
    python tools/patches/005d_fix_cantera_hpy_api.py --apply
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
    p = argparse.ArgumentParser(description="Patch 005d: fix Cantera HPY inversion API.")
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


def patch_cantera_interface(repo: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "cantera_interface.cpp"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    old = (
        "                build_cantera_mass_fractions(*gas, nspecies, Y_in, c, sp_map, Y_cantera);\n"
        "                gas->setState_HPY(h_in[c], P[c], Y_cantera.data());\n"
        "                T_out[c] = gas->temperature();\n"
    )
    new = (
        "                build_cantera_mass_fractions(*gas, nspecies, Y_in, c, sp_map, Y_cantera);\n"
        "                // Cantera ThermoPhase in this environment does not provide setState_HPY.\n"
        "                // Set composition first, then invert h,p at fixed composition.\n"
        "                gas->setMassFractions(Y_cantera.data());\n"
        "                gas->setState_HP(h_in[c], P[c]);\n"
        "                T_out[c] = gas->temperature();\n"
    )

    if old in text:
        text = text.replace(old, new, 1)
    elif "gas->setState_HPY" in text:
        raise SystemExit("[ERROR] setState_HPY found, but surrounding anchor did not match expected Patch 005c form")
    elif "gas->setMassFractions(Y_cantera.data());" in text and "gas->setState_HP(h_in[c], P[c]);" in text:
        return PatchResult(path, False, "already patched")
    else:
        raise SystemExit("[ERROR] Cantera temperature-recovery block not found; apply Patch 005c first")

    if text == original:
        return PatchResult(path, False, "unchanged")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    return PatchResult(path, True, "would replace setState_HPY with setMassFractions + setState_HP")


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_cantera_interface(repo, apply)

    print("Patch 005d summary")
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
        print("  3. run the Stage-2A smoke test with enable_cantera_thermo = .true.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
