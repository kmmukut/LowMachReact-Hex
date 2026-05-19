#!/usr/bin/env python3
"""
Patch 002 v2: add energy VTU output and energy diagnostics.

This patch assumes Patch 001 has already added:
  - src/mod_energy.f90
  - energy_fields_t
  - enable_energy in case_params_t
  - main.f90 initialization/finalization of energy

Patch 002 intentionally does NOT add:
  - energy transport
  - temperature boundary conditions
  - Cantera h/T inversion
  - radiation physics

It adds:
  - energy_diagnostics.csv with min/max/mean T, h, and qrad integral
  - VTU/PVTU cell data arrays for temperature, enthalpy, and qrad
  - main.f90 calls to write energy diagnostics
  - Makefile dependency updates for mod_energy/mod_output

Run from repository root:

  python tools/patches/002_add_energy_output_diagnostics.py --dry-run
  python tools/patches/002_add_energy_output_diagnostics.py --apply
"""

from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PatchResult:
    path: Path
    changed: bool
    message: str


ENERGY_DIAGNOSTIC_SUBROUTINES = r"""
   !> Writes the CSV header for energy diagnostics.
   subroutine write_energy_diagnostics_header(params, flow)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow

      integer :: unit_id
      character(len=256 + 32) :: filename

      if (flow%rank /= 0 .or. .not. params%write_diagnostics) return
      if (.not. params%enable_energy) return

      filename = trim(params%output_dir)//'/energy_diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') 'step,time,min_T,max_T,mean_T,min_h,max_h,mean_h,min_qrad,max_qrad,integral_qrad'

      close(unit_id)
   end subroutine write_energy_diagnostics_header


   !> Appends one row of global energy diagnostics.
   !!
   !! The temperature and enthalpy means are volume weighted over owned cells.
   !! qrad integral is the domain integral of qrad dV [W].
   subroutine write_energy_diagnostics_row(mesh, flow, params, energy, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(in) :: energy
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, unit_id, ierr
      character(len=256 + 32) :: filename
      real(rk) :: vol
      real(rk) :: local_min_T, local_max_T, local_sum_T
      real(rk) :: local_min_h, local_max_h, local_sum_h
      real(rk) :: local_min_qrad, local_max_qrad, local_integral_qrad
      real(rk) :: local_volume
      real(rk) :: global_min_T, global_max_T, global_sum_T
      real(rk) :: global_min_h, global_max_h, global_sum_h
      real(rk) :: global_min_qrad, global_max_qrad, global_integral_qrad
      real(rk) :: global_volume
      real(rk) :: mean_T, mean_h

      if (.not. params%write_diagnostics) return
      if (.not. params%enable_energy) return

      if (.not. allocated(energy%T) .or. .not. allocated(energy%h) .or. &
          .not. allocated(energy%qrad)) then
         call fatal_error('energy', 'energy diagnostics requested but energy arrays are not allocated')
      end if

      local_min_T = huge(0.0_rk)
      local_max_T = -huge(0.0_rk)
      local_sum_T = zero
      local_min_h = huge(0.0_rk)
      local_max_h = -huge(0.0_rk)
      local_sum_h = zero
      local_min_qrad = huge(0.0_rk)
      local_max_qrad = -huge(0.0_rk)
      local_integral_qrad = zero
      local_volume = zero

      do c = 1, mesh%ncells
         if (.not. flow%owned(c)) cycle

         vol = mesh%cells(c)%volume
         local_volume = local_volume + vol

         local_min_T = min(local_min_T, energy%T(c))
         local_max_T = max(local_max_T, energy%T(c))
         local_sum_T = local_sum_T + energy%T(c) * vol

         local_min_h = min(local_min_h, energy%h(c))
         local_max_h = max(local_max_h, energy%h(c))
         local_sum_h = local_sum_h + energy%h(c) * vol

         local_min_qrad = min(local_min_qrad, energy%qrad(c))
         local_max_qrad = max(local_max_qrad, energy%qrad(c))
         local_integral_qrad = local_integral_qrad + energy%qrad(c) * vol
      end do

      call MPI_Allreduce(local_min_T, global_min_T, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing min_T')
      call MPI_Allreduce(local_max_T, global_max_T, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing max_T')
      call MPI_Allreduce(local_sum_T, global_sum_T, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing sum_T')

      call MPI_Allreduce(local_min_h, global_min_h, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing min_h')
      call MPI_Allreduce(local_max_h, global_max_h, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing max_h')
      call MPI_Allreduce(local_sum_h, global_sum_h, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing sum_h')

      call MPI_Allreduce(local_min_qrad, global_min_qrad, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing min_qrad')
      call MPI_Allreduce(local_max_qrad, global_max_qrad, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing max_qrad')
      call MPI_Allreduce(local_integral_qrad, global_integral_qrad, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing integral_qrad')

      call MPI_Allreduce(local_volume, global_volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing energy volume')

      if (global_volume > tiny_safe) then
         mean_T = global_sum_T / global_volume
         mean_h = global_sum_h / global_volume
      else
         mean_T = zero
         mean_h = zero
      end if

      if (flow%rank /= 0) return

      filename = trim(params%output_dir)//'/energy_diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='old', position='append', action='write')

      write(unit_id,'(i0,10(",",es16.8))') step, time, global_min_T, global_max_T, mean_T, &
         global_min_h, global_max_h, mean_h, global_min_qrad, global_max_qrad, global_integral_qrad

      close(unit_id)
   end subroutine write_energy_diagnostics_row

"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Patch 002: add energy output/diagnostics.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true", help="write changes to disk")
    mode.add_argument("--dry-run", action="store_true", help="show planned changes without writing")
    parser.add_argument("--repo", default=".", help="repository root; default: current directory")
    parser.add_argument("--skip-build-file", action="store_true", help="do not modify Makefile")
    return parser.parse_args()


def timestamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] required file not found: {path}")


def make_backup(path: Path, repo: Path, backup_root: Path, apply: bool) -> None:
    if not apply:
        return
    rel = path.relative_to(repo)
    dst = backup_root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, dst)


def write_if_changed(path: Path, text: str, repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    old = read_text(path)
    if old == text:
        return PatchResult(path, False, "unchanged")

    if apply:
        make_backup(path, repo, backup_root, apply=True)
        path.write_text(text, encoding="utf-8")

    return PatchResult(path, True, "would update" if not apply else "updated")


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if count > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {count} matches")
    return text.replace(old, new, 1)


def insert_after_once(text: str, anchor: str, insertion: str, label: str) -> str:
    count = text.count(anchor)
    if count == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if count > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {count} matches")
    return text.replace(anchor, anchor + insertion, 1)


def patch_mod_energy(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_energy.f90"
    text = read_text(path)

    if "write_energy_diagnostics_header" in text and "write_energy_diagnostics_row" in text:
        return PatchResult(path, False, "already patched")

    if "   use mpi_f08\n" not in text:
        text = replace_once(
            text,
            "module mod_energy\n   use mod_kinds",
            "module mod_energy\n   use mpi_f08\n   use mod_kinds",
            "add mpi_f08 use to mod_energy",
        )

    if "use mod_mpi_flow, only : flow_mpi_t" not in text:
        text = insert_after_once(
            text,
            "   use mod_input, only : case_params_t\n",
            "   use mod_mpi_flow, only : flow_mpi_t\n",
            "mod_energy use mod_mpi_flow",
        )

    if "public :: write_energy_diagnostics_header, write_energy_diagnostics_row" not in text:
        text = insert_after_once(
            text,
            "   public :: zero_radiation_source\n",
            "   public :: write_energy_diagnostics_header, write_energy_diagnostics_row\n",
            "mod_energy public diagnostics",
        )

    text = replace_once(
        text,
        "\nend module mod_energy\n",
        ENERGY_DIAGNOSTIC_SUBROUTINES + "end module mod_energy\n",
        "insert energy diagnostic subroutines",
    )

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_mod_output(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_output.f90"
    text = read_text(path)

    if "type(energy_fields_t), intent(in) :: energy" in text and 'Name="temperature"' in text:
        return PatchResult(path, False, "already patched")

    if "use mod_energy, only : energy_fields_t" not in text:
        text = insert_after_once(
            text,
            "   use mod_species, only : species_fields_t\n",
            "   use mod_energy, only : energy_fields_t\n",
            "mod_output use mod_energy",
        )

    text = replace_once(
        text,
        "   subroutine write_vtu_unstructured(params, flow, mesh, fields, species, step)\n",
        "   subroutine write_vtu_unstructured(params, flow, mesh, fields, species, energy, step)\n",
        "write_vtu_unstructured signature",
    )

    text = insert_after_once(
        text,
        "      type(species_fields_t), intent(in) :: species\n",
        "      type(energy_fields_t), intent(in) :: energy\n",
        "write_vtu_unstructured energy argument",
    )

    energy_arrays = """\n      if (params%enable_energy) then\n         if (.not. allocated(energy%T) .or. .not. allocated(energy%h) .or. &\n             .not. allocated(energy%qrad)) then\n            call fatal_error('output', 'energy output requested but energy arrays are not allocated')\n         end if\n\n         write(unit_id,'(a)') '        <DataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\">'\n         do n = 1, n_owned\n            c = owned_indices(n)\n            write(unit_id,'(es24.16)') energy%T(c)\n         end do\n         write(unit_id,'(a)') '        </DataArray>'\n\n         write(unit_id,'(a)') '        <DataArray type=\"Float64\" Name=\"enthalpy\" format=\"ascii\">'\n         do n = 1, n_owned\n            c = owned_indices(n)\n            write(unit_id,'(es24.16)') energy%h(c)\n         end do\n         write(unit_id,'(a)') '        </DataArray>'\n\n         write(unit_id,'(a)') '        <DataArray type=\"Float64\" Name=\"qrad\" format=\"ascii\">'\n         do n = 1, n_owned\n            c = owned_indices(n)\n            write(unit_id,'(es24.16)') energy%qrad(c)\n         end do\n         write(unit_id,'(a)') '        </DataArray>'\n      end if\n\n"""
    text = insert_after_once(
        text,
        "      write(unit_id,'(a)') '        </DataArray>'\n\n      if (params%enable_species .and. params%nspecies > 0) then\n",
        energy_arrays,
        "insert energy VTU arrays",
    )

    penergy_arrays = """      if (params%enable_energy) then\n         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"temperature\" format=\"ascii\"/>'\n         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"enthalpy\" format=\"ascii\"/>'\n         write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"qrad\" format=\"ascii\"/>'\n      end if\n"""
    text = insert_after_once(
        text,
        "      write(unit_id,'(a)') '      <PDataArray type=\"Float64\" Name=\"divergence\" format=\"ascii\"/>'\n",
        penergy_arrays,
        "insert energy PVTU arrays",
    )

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_main(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "main.f90"
    text = read_text(path)

    if "write_energy_diagnostics_header" in text and "write_vtu_unstructured(params, flow_mpi, mesh, fields, species, energy" in text:
        return PatchResult(path, False, "already patched")

    text = replace_once(
        text,
        "   use mod_energy, only : energy_fields_t, initialize_energy, finalize_energy\n",
        "   use mod_energy, only : energy_fields_t, initialize_energy, finalize_energy, &\n                           write_energy_diagnostics_header, write_energy_diagnostics_row\n",
        "main use mod_energy diagnostics",
    )

    if "call write_energy_diagnostics_header(params, flow_mpi)" not in text:
        text = insert_after_once(
            text,
            "   call write_diagnostics_header(params, flow_mpi)\n",
            "   call write_energy_diagnostics_header(params, flow_mpi)\n",
            "main energy diagnostics header",
        )

    if "call write_energy_diagnostics_row(mesh, flow_mpi, params, energy, 0, time)" not in text:
        text = insert_after_once(
            text,
            "   call write_diagnostics_row(params, flow_mpi, 0, time, stats)\n",
            "   call write_energy_diagnostics_row(mesh, flow_mpi, params, energy, 0, time)\n",
            "main initial energy diagnostics row",
        )

    text = text.replace(
        "call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, 0)",
        "call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, energy, 0)",
    )
    text = text.replace(
        "call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, step)",
        "call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, energy, step)",
    )

    if "call write_energy_diagnostics_row(mesh, flow_mpi, params, energy, step, time)" not in text:
        text = insert_after_once(
            text,
            "         call write_diagnostics_row(params, flow_mpi, step, time, stats)\n",
            "         call write_energy_diagnostics_row(mesh, flow_mpi, params, energy, step, time)\n",
            "main step energy diagnostics row",
        )

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_makefile(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "Makefile"
    text = read_text(path)
    original = text

    if "$(BUILD_DIR)/mod_energy.o:" in text:
        old = "$(BUILD_DIR)/mod_energy.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o\n"
        new = "$(BUILD_DIR)/mod_energy.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mpi_flow.o\n"
        if old in text:
            text = text.replace(old, new, 1)
        elif new in text:
            pass
        else:
            raise SystemExit("[ERROR] found mod_energy.o rule but it is not in the expected form")
    else:
        anchor = "$(BUILD_DIR)/mod_fields.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_bc.o\n"
        dep_line = "$(BUILD_DIR)/mod_energy.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mpi_flow.o\n"
        text = insert_after_once(text, anchor, dep_line, "add mod_energy.o dependency rule")

    old = "$(BUILD_DIR)/mod_output.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_species.o\n"
    new = "$(BUILD_DIR)/mod_output.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_energy.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_species.o\n"
    if old in text:
        text = text.replace(old, new, 1)
    elif new in text:
        pass
    else:
        raise SystemExit("[ERROR] mod_output.o dependency rule not found in expected form")

    if text == original:
        return PatchResult(path, False, "already patched or unchanged")

    return write_if_changed(path, text, repo, backup_root, apply)


def sanity_check_repo(repo: Path) -> None:
    required = [
        repo / "src" / "main.f90",
        repo / "src" / "mod_input.f90",
        repo / "src" / "mod_energy.f90",
        repo / "src" / "mod_output.f90",
    ]
    missing = [p for p in required if not p.exists()]
    if missing:
        for p in missing:
            print(f"[ERROR] missing required file: {p}", file=sys.stderr)
        raise SystemExit(2)


def main() -> int:
    args = parse_args()
    apply = bool(args.apply)
    if not args.apply and not args.dry_run:
        print("[INFO] neither --apply nor --dry-run supplied; defaulting to --dry-run")
        apply = False

    repo = Path(args.repo).resolve()
    sanity_check_repo(repo)

    backup_root = repo / ".backups" / timestamp()
    results: list[PatchResult] = []

    results.append(patch_mod_energy(repo, backup_root, apply))
    results.append(patch_mod_output(repo, backup_root, apply))
    results.append(patch_main(repo, backup_root, apply))
    if args.skip_build_file:
        results.append(PatchResult(repo / "Makefile", False, "skipped by --skip-build-file"))
    else:
        results.append(patch_makefile(repo, backup_root, apply))

    print("\nPatch 002 summary")
    print("=================")
    print(f"repo:   {repo}")
    print(f"mode:   {'apply' if apply else 'dry-run'}")
    if apply:
        print(f"backup: {backup_root}")
    print()

    changed_any = False
    for r in results:
        changed_any = changed_any or r.changed
        rel = r.path
        try:
            rel = r.path.relative_to(repo)
        except ValueError:
            pass
        marker = "CHANGE" if r.changed else "SKIP"
        print(f"[{marker}] {rel}: {r.message}")

    if not apply:
        print("\nDry run only. Re-run with --apply to write changes.")
    elif changed_any:
        print("\nApplied patch. Next steps:")
        print("  1. inspect git diff")
        print("  2. make clean && make -j")
        print("  3. run with enable_energy = .false. to confirm old behavior")
        print("  4. run with enable_energy = .true. and inspect energy_diagnostics.csv and VTU arrays")
    else:
        print("\nNo changes were needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
