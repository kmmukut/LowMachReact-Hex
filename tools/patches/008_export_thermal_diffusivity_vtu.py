#!/usr/bin/env python3
"""
Patch 008: export active thermal diffusivity to VTU/PVTU.

Adds a derived VTU/PVTU cell array:

    thermal_diffusivity = thermal_conductivity / (rho_flow * cp)

where:
  - thermal_conductivity = energy%lambda
  - cp = energy%cp
  - rho_flow = params%rho

This is the thermal diffusivity actually consistent with the current
constant-density Stage-2A energy equation. rho_thermo is diagnostic only and is
not used in the projection/energy mass coefficient yet.

Files modified:
  - src/mod_output.f90

Run from repository root:

    python tools/patches/008_export_thermal_diffusivity_vtu.py --dry-run
    python tools/patches/008_export_thermal_diffusivity_vtu.py --apply
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


THERMAL_DIFFUSIVITY_VTU_BLOCK = r"""
         if (allocated(energy%lambda) .and. allocated(energy%cp)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="thermal_diffusivity" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(es24.16)') energy%lambda(c) / max(params%rho * energy%cp(c), tiny(1.0_rk))
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

"""

THERMAL_DIFFUSIVITY_PVTU_LINE = (
    "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"thermal_diffusivity\" format=\"ascii\"/>'\n"
)


def parse_args():
    p = argparse.ArgumentParser(description="Patch 008: export thermal diffusivity to VTU.")
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


def patch_vtu_piece(text: str) -> tuple[str, bool]:
    sub_start, sub_end = find_subroutine_bounds(text, "write_vtu_unstructured")
    body = text[sub_start:sub_end]

    if 'Name="thermal_diffusivity"' in body:
        return text, False

    if "energy%lambda" not in body or "energy%cp" not in body:
        raise SystemExit(
            "[ERROR] VTU output does not yet reference energy%lambda and energy%cp. "
            "Apply/fix Patch 007b/007c first."
        )

    # Insert after the thermal_conductivity DataArray block.
    cond_pos = body.find('Name="thermal_conductivity"')
    if cond_pos < 0:
        raise SystemExit("[ERROR] thermal_conductivity DataArray not found in write_vtu_unstructured")

    close_marker = "            write(unit_id,'(a)') '        </DataArray>'\n"
    close_pos = body.find(close_marker, cond_pos)
    if close_pos < 0:
        raise SystemExit("[ERROR] thermal_conductivity closing DataArray write not found")

    insert_pos = close_pos + len(close_marker)

    # Include the enclosing end-if block if it immediately follows, so the
    # derived field is inserted after the conductivity block, not inside it.
    following_end_if = "         end if\n"
    if body[insert_pos:insert_pos + len(following_end_if)] == following_end_if:
        insert_pos += len(following_end_if)
        if body[insert_pos:insert_pos + 1] == "\n":
            insert_pos += 1

    new_body = body[:insert_pos] + THERMAL_DIFFUSIVITY_VTU_BLOCK + body[insert_pos:]
    return text[:sub_start] + new_body + text[sub_end:], True


def patch_pvtu(text: str) -> tuple[str, bool]:
    if 'Name="thermal_diffusivity"' in text:
        return text, False

    # Prefer inserting after thermal_conductivity metadata if present.
    anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"thermal_conductivity\" format=\"ascii\"/>'\n"
    if anchor not in text:
        # Fallback to qrad if the thermo metadata patch has not placed conductivity
        # in PVTU yet.
        anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"qrad\" format=\"ascii\"/>'\n"
        if anchor not in text:
            raise SystemExit("[ERROR] PVTU insertion anchor not found")

    return text.replace(anchor, anchor + THERMAL_DIFFUSIVITY_PVTU_LINE, 1), True


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

    parts = []
    if changed_vtu:
        parts.append("VTU piece array")
    if changed_pvtu:
        parts.append("PVTU metadata")
    return PatchResult(path, True, "would add " + ", ".join(parts))


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_mod_output(repo, apply)

    print("Patch 008 summary")
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
        print("  3. delete old output or use a fresh output directory before rerunning")
        print("  4. inspect VTU/PVTU array: thermal_diffusivity")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
