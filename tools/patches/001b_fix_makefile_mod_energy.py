#!/usr/bin/env python3
"""
Patch 001b: fix Makefile for mod_energy.

Use this if Patch 001 added src/mod_energy.f90 and main.f90 now uses
mod_energy, but the build fails with:

  Fatal Error: Cannot open module file 'mod_energy.mod'

Run from repo root:

  python tools/patches/001b_fix_makefile_mod_energy.py --dry-run
  python tools/patches/001b_fix_makefile_mod_energy.py --apply

This patch only modifies the root Makefile and creates a timestamped backup
under .backups/YYYYMMDD_HHMMSS/Makefile before writing.
"""

from __future__ import annotations

import argparse
import datetime as dt
import shutil
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    g = p.add_mutually_exclusive_group()
    g.add_argument("--apply", action="store_true")
    g.add_argument("--dry-run", action="store_true")
    p.add_argument("--repo", default=".")
    return p.parse_args()


def backup(path: Path, repo: Path) -> Path:
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    dst = repo / ".backups" / stamp / path.relative_to(repo)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, dst)
    return dst


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n == 0:
        raise SystemExit(f"[ERROR] anchor not found: {label}")
    if n > 1:
        raise SystemExit(f"[ERROR] ambiguous anchor for {label}: {n} matches")
    return text.replace(old, new, 1)


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    makefile = repo / "Makefile"
    if not makefile.exists():
        raise SystemExit(f"[ERROR] not found: {makefile}")

    text = makefile.read_text(encoding="utf-8")
    original = text
    changes = []

    # 1. Add mod_energy.f90 to ordered Fortran source list.
    src_line = "  $(SRC_DIR)/mod_energy.f90 \\\n"
    if src_line not in text:
        anchor = "  $(SRC_DIR)/mod_fields.f90 \\\n"
        text = replace_once(text, anchor, anchor + src_line, "F_SRCS mod_fields.f90")
        changes.append("add src/mod_energy.f90 to F_SRCS after mod_fields.f90")
    else:
        changes.append("skip F_SRCS; mod_energy.f90 already present")

    # 2. Add explicit object dependency rule for mod_energy.
    dep_line = "$(BUILD_DIR)/mod_energy.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o\n"
    if dep_line not in text:
        anchor = "$(BUILD_DIR)/mod_fields.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_bc.o\n"
        text = replace_once(text, anchor, anchor + dep_line, "mod_fields.o dependency rule")
        changes.append("add mod_energy.o dependency rule")
    else:
        changes.append("skip dependency rule; mod_energy.o rule already present")

    # 3. Add mod_energy.o as an explicit dependency of main.o.
    main_old = "$(BUILD_DIR)/main.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mesh_io.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_mpi_radiation.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_transport_properties.o $(BUILD_DIR)/mod_species.o $(BUILD_DIR)/mod_output.o $(BUILD_DIR)/mod_profiler.o\n"
    main_new = "$(BUILD_DIR)/main.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mesh_io.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_mpi_radiation.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_energy.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_transport_properties.o $(BUILD_DIR)/mod_species.o $(BUILD_DIR)/mod_output.o $(BUILD_DIR)/mod_profiler.o\n"

    main_lines = [line for line in text.splitlines(True) if line.startswith("$(BUILD_DIR)/main.o:")]
    if not main_lines:
        raise SystemExit("[ERROR] main.o dependency line not found")
    if "$(BUILD_DIR)/mod_energy.o" not in main_lines[0]:
        text = replace_once(text, main_old, main_new, "main.o dependency rule")
        changes.append("add mod_energy.o to main.o dependencies")
    else:
        changes.append("skip main.o dependency; mod_energy.o already present")

    print("Patch 001b Makefile summary")
    print("===========================")
    print(f"repo: {repo}")
    print(f"mode: {'apply' if apply else 'dry-run'}")
    for c in changes:
        print(f"- {c}")

    if text == original:
        print("\nNo Makefile changes needed.")
        return 0

    if apply:
        dst = backup(makefile, repo)
        makefile.write_text(text, encoding="utf-8")
        print(f"\nUpdated Makefile")
        print(f"Backup: {dst}")
    else:
        print("\nDry run only. Re-run with --apply to write changes.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
