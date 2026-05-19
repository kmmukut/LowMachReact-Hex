#!/usr/bin/env python3
"""
Patch 008b: repair PVTU metadata for thermal_diffusivity.

Problem:
  Patch 008 may add the thermal_diffusivity array to VTU piece files, but skip
  the PVTU metadata because it used a global string check. If the string appears
  in the VTU writer, the patch may incorrectly assume the PVTU declaration also
  exists. ParaView normally opens the .pvd/.pvtu layer, so missing PDataArray
  metadata can hide the array.

Fix:
  Insert the PVTU metadata line locally in the PVTU CellData section, after
  thermal_conductivity when available, otherwise after qrad.

  Also verifies the VTU piece writer has the thermal_diffusivity DataArray.

Files modified:
  - src/mod_output.f90

Run from repository root:

    python tools/patches/008b_fix_thermal_diffusivity_pvtu.py --dry-run
    python tools/patches/008b_fix_thermal_diffusivity_pvtu.py --apply
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


PVTU_LINE = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"thermal_diffusivity\" format=\"ascii\"/>'\n"

VTU_BLOCK = r"""
         if (allocated(energy%lambda) .and. allocated(energy%cp)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="thermal_diffusivity" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(es24.16)') energy%lambda(c) / max(params%rho * energy%cp(c), tiny(1.0_rk))
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

"""


def parse_args():
    p = argparse.ArgumentParser(description="Patch 008b: fix thermal_diffusivity PVTU metadata.")
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


def patch_vtu_piece_if_missing(text: str) -> tuple[str, bool]:
    sub_start, sub_end = find_subroutine_bounds(text, "write_vtu_unstructured")
    body = text[sub_start:sub_end]

    if 'Name="thermal_diffusivity"' in body:
        return text, False

    if "energy%lambda" not in body or "energy%cp" not in body:
        raise SystemExit(
            "[ERROR] write_vtu_unstructured does not contain energy%lambda/energy%cp. "
            "Apply Patch 007b/007c first."
        )

    anchor = 'Name="thermal_conductivity"'
    cond_pos = body.find(anchor)
    if cond_pos >= 0:
        close_marker = "            write(unit_id,'(a)') '        </DataArray>'\n"
        close_pos = body.find(close_marker, cond_pos)
        if close_pos < 0:
            raise SystemExit("[ERROR] thermal_conductivity closing DataArray write not found")
        insert_pos = close_pos + len(close_marker)

        end_if = "         end if\n"
        if body[insert_pos:insert_pos + len(end_if)] == end_if:
            insert_pos += len(end_if)
            if body[insert_pos:insert_pos + 1] == "\n":
                insert_pos += 1
    else:
        # Fallback: insert after qrad DataArray.
        qrad_pos = body.find('Name="qrad"')
        if qrad_pos < 0:
            raise SystemExit("[ERROR] neither thermal_conductivity nor qrad DataArray found in VTU writer")
        close_marker = "         write(unit_id,'(a)') '        </DataArray>'\n"
        close_pos = body.find(close_marker, qrad_pos)
        if close_pos < 0:
            raise SystemExit("[ERROR] qrad closing DataArray write not found")
        insert_pos = close_pos + len(close_marker)

    new_body = body[:insert_pos] + VTU_BLOCK + body[insert_pos:]
    return text[:sub_start] + new_body + text[sub_end:], True


def patch_pvtu_metadata(text: str) -> tuple[str, bool]:
    changed = False

    # Patch each PVTU qrad metadata region independently. Do not use a global
    # string check, because the same name can exist in VTU writer code.
    search_pos = 0
    while True:
        qrad_anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"qrad\" format=\"ascii\"/>'\n"
        qrad_idx = text.find(qrad_anchor, search_pos)
        if qrad_idx < 0:
            break

        # Look only in the local PVTU metadata neighborhood after qrad.
        local = text[qrad_idx:qrad_idx + 1500]
        if 'PDataArray type="Float64" Name="thermal_diffusivity"' in local:
            search_pos = qrad_idx + len(qrad_anchor)
            continue

        # Prefer inserting after thermal_conductivity if it exists locally;
        # otherwise insert after qrad.
        cond_anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"thermal_conductivity\" format=\"ascii\"/>'\n"
        cond_idx = text.find(cond_anchor, qrad_idx, qrad_idx + 1500)

        if cond_idx >= 0:
            insert_at = cond_idx + len(cond_anchor)
        else:
            insert_at = qrad_idx + len(qrad_anchor)

        text = text[:insert_at] + PVTU_LINE + text[insert_at:]
        changed = True
        search_pos = insert_at + len(PVTU_LINE)

    if not changed and "PDataArray" not in text:
        raise SystemExit("[ERROR] no PVTU metadata section found in mod_output.f90")

    if not changed:
        # qrad may not exist in some customized output file. Try after thermal_conductivity.
        cond_anchor = "         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"thermal_conductivity\" format=\"ascii\"/>'\n"
        cond_idx = text.find(cond_anchor)
        if cond_idx >= 0:
            local = text[cond_idx:cond_idx + 1000]
            if 'PDataArray type="Float64" Name="thermal_diffusivity"' not in local:
                insert_at = cond_idx + len(cond_anchor)
                text = text[:insert_at] + PVTU_LINE + text[insert_at:]
                changed = True

    return text, changed


def patch_mod_output(repo: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_output.f90"
    if not path.exists():
        raise SystemExit(f"[ERROR] not found: {path}")

    text = path.read_text(encoding="utf-8")
    original = text

    text, changed_vtu = patch_vtu_piece_if_missing(text)
    text, changed_pvtu = patch_pvtu_metadata(text)

    if text == original:
        return PatchResult(path, False, "already patched")

    if apply:
        backup_root = repo / ".backups" / timestamp()
        dst = backup(path, repo, backup_root)
        path.write_text(text, encoding="utf-8")
        return PatchResult(path, True, f"updated; backup: {dst}")

    parts = []
    if changed_vtu:
        parts.append("VTU piece DataArray")
    if changed_pvtu:
        parts.append("PVTU PDataArray metadata")
    return PatchResult(path, True, "would update " + ", ".join(parts))


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    result = patch_mod_output(repo, apply)

    print("Patch 008b summary")
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
        print("  3. delete old output or use a fresh output directory")
        print("  4. rerun and open the new .pvd/.pvtu in ParaView")
        print("  5. optionally grep output: grep -R \"thermal_diffusivity\" cases/rectangle_2D/output")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
