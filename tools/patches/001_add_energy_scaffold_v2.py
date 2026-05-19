#!/usr/bin/env python3
"""
Patch 001: add enthalpy energy scaffolding.

This patch intentionally does NOT add energy transport, Cantera enthalpy
inversion, radiation physics, or output changes. It only adds:

  - src/mod_energy.f90
  - basic &energy_input namelist support in src/mod_input.f90
  - energy field initialize/finalize wiring in src/main.f90
  - Makefile insertion for this repository's explicit F_SRCS/dependency layout

Run from the repository root:

  python tools/patches/001_add_energy_scaffold.py --dry-run
  python tools/patches/001_add_energy_scaffold.py --apply

Backups of modified files are written to:

  .backups/YYYYMMDD_HHMMSS/
"""

from __future__ import annotations

import argparse
import datetime as _dt
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


MOD_ENERGY = """!> Energy field storage and thermodynamic helper routines.
!!
!! This module introduces the enthalpy/temperature/radiation-source storage
!! needed for the staged energy-equation implementation.
!!
!! Patch 001 intentionally provides scaffolding only:
!! - No energy transport is advanced here.
!! - No Cantera enthalpy inversion is called here.
!! - No radiation physics is evaluated here.
!!
!! The first implementation convention is:
!!   h    : mixture sensible enthalpy [J/kg]
!!   T    : temperature [K]
!!   qrad : volumetric radiation source [W/m^3]
!!
!! For now, h and T are related through a constant-cp reference relation:
!!   h = h_ref + cp * (T - T_ref)
!!
!! Later patches should replace the constant-cp conversion with Cantera
!! thermodynamic calls while preserving this field ownership.
module mod_energy
   use mod_kinds, only : rk, zero, tiny_safe, fatal_error
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   implicit none

   private

   !> Cell-centered energy variables.
   type, public :: energy_fields_t
      real(rk), allocatable :: T(:)     !< Temperature [K].
      real(rk), allocatable :: h(:)     !< Mixture sensible enthalpy [J/kg].
      real(rk), allocatable :: h_old(:) !< Previous enthalpy state [J/kg].
      real(rk), allocatable :: qrad(:)  !< Volumetric radiation source [W/m^3].
      logical :: initialized = .false.  !< True after allocation/initialization.
   end type energy_fields_t

   public :: allocate_energy, initialize_energy, finalize_energy
   public :: update_enthalpy_from_temperature_constant_cp
   public :: recover_temperature_constant_cp
   public :: zero_radiation_source

contains

   !> Allocate all energy arrays for the mesh.
   subroutine allocate_energy(mesh, energy)
      type(mesh_t), intent(in) :: mesh
      type(energy_fields_t), intent(inout) :: energy

      call finalize_energy(energy)

      allocate(energy%T(mesh%ncells))
      allocate(energy%h(mesh%ncells))
      allocate(energy%h_old(mesh%ncells))
      allocate(energy%qrad(mesh%ncells))

      energy%T = zero
      energy%h = zero
      energy%h_old = zero
      energy%qrad = zero
      energy%initialized = .true.
   end subroutine allocate_energy


   !> Initialize energy fields from case parameters.
   subroutine initialize_energy(mesh, params, energy)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy

      call allocate_energy(mesh, energy)

      energy%T = params%initial_T
      call update_enthalpy_from_temperature_constant_cp(params, energy)
      energy%h_old = energy%h
      call zero_radiation_source(energy)
   end subroutine initialize_energy


   !> Deallocate all energy arrays.
   subroutine finalize_energy(energy)
      type(energy_fields_t), intent(inout) :: energy

      if (allocated(energy%T)) deallocate(energy%T)
      if (allocated(energy%h)) deallocate(energy%h)
      if (allocated(energy%h_old)) deallocate(energy%h_old)
      if (allocated(energy%qrad)) deallocate(energy%qrad)

      energy%initialized = .false.
   end subroutine finalize_energy


   !> Update h from T using the temporary constant-cp thermodynamic model.
   subroutine update_enthalpy_from_temperature_constant_cp(params, energy)
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy

      if (.not. allocated(energy%T) .or. .not. allocated(energy%h)) then
         call fatal_error('energy', 'energy arrays are not allocated')
      end if

      if (params%energy_cp <= tiny_safe) then
         call fatal_error('energy', 'energy_cp must be positive')
      end if

      energy%h = params%energy_reference_h + &
                 params%energy_cp * (energy%T - params%energy_reference_T)
   end subroutine update_enthalpy_from_temperature_constant_cp


   !> Recover T from h using the temporary constant-cp thermodynamic model.
   subroutine recover_temperature_constant_cp(params, energy)
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy

      if (.not. allocated(energy%T) .or. .not. allocated(energy%h)) then
         call fatal_error('energy', 'energy arrays are not allocated')
      end if

      if (params%energy_cp <= tiny_safe) then
         call fatal_error('energy', 'energy_cp must be positive')
      end if

      energy%T = params%energy_reference_T + &
                 (energy%h - params%energy_reference_h) / params%energy_cp
   end subroutine recover_temperature_constant_cp


   !> Reset the radiation source to zero.
   subroutine zero_radiation_source(energy)
      type(energy_fields_t), intent(inout) :: energy

      if (allocated(energy%qrad)) energy%qrad = zero
   end subroutine zero_radiation_source

end module mod_energy
"""


@dataclass
class PatchResult:
    path: Path
    changed: bool
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Patch 001: add energy scaffolding.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true", help="write changes to disk")
    mode.add_argument("--dry-run", action="store_true", help="show planned changes without writing")
    parser.add_argument("--repo", default=".", help="repository root; default: current directory")
    parser.add_argument(
        "--skip-build-file",
        action="store_true",
        help="do not modify Makefile/build files; useful if you want to edit build ordering manually",
    )
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


def write_if_changed(path: Path, text: str, repo: Path, backup_root: Path, apply: bool, *,
                     must_exist: bool = True) -> PatchResult:
    old = None
    if path.exists():
        old = path.read_text(encoding="utf-8")
        if old == text:
            return PatchResult(path, False, "unchanged")
    elif must_exist:
        raise SystemExit(f"[ERROR] required file not found: {path}")

    if apply:
        if path.exists():
            make_backup(path, repo, backup_root, apply=True)
        path.parent.mkdir(parents=True, exist_ok=True)
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


def patch_mod_input(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_input.f90"
    text = read_text(path)

    if "logical :: enable_energy" in text and "subroutine read_energy_input" in text:
        return PatchResult(path, False, "already patched")

    field_anchor = (
        "      real(rk) :: background_temp = 300.0_rk                      !< Fixed temperature for property evaluation [K].\n"
        "      real(rk) :: background_press = 101325.0_rk                  !< Fixed pressure for property evaluation [Pa].\n"
    )
    field_insert = (
        "\n"
        "      !> @name Enthalpy Energy Equation Controls\n"
        "      logical :: enable_energy = .false.                         !< Enable enthalpy/temperature field storage.\n"
        "      logical :: enable_cantera_thermo = .false.                 !< Use Cantera for h(T,Y,p) and T(h,Y,p) in later patches.\n"
        "      real(rk) :: initial_T = 300.0_rk                           !< Initial gas temperature [K].\n"
        "      real(rk) :: energy_reference_T = 298.15_rk                 !< Reference temperature for temporary constant-cp h model [K].\n"
        "      real(rk) :: energy_reference_h = zero                      !< Reference sensible enthalpy for temporary h model [J/kg].\n"
        "      real(rk) :: energy_cp = 1005.0_rk                          !< Temporary constant heat capacity [J/kg/K].\n"
        "      real(rk) :: energy_lambda = 2.6e-2_rk                      !< Temporary constant thermal conductivity [W/m/K].\n"
    )
    if "logical :: enable_energy" not in text:
        text = insert_after_once(text, field_anchor, field_insert, "case_params_t energy fields")

    call_anchor = "      call read_species_input(filename, params)\n"
    call_insert = "      call read_energy_input(filename, params)\n"
    if "call read_energy_input(filename, params)" not in text:
        text = insert_after_once(text, call_anchor, call_insert, "read_case_params energy call")

    sub_anchor = "   !> Reads the `&output_input` namelist block.\n"
    read_energy_sub = """   !> Reads the `&energy_input` namelist block.
   !!
   !! Patch 001 only enables storage/initialization of energy fields.
   !! Later patches will add energy transport and Cantera thermodynamic inversion.
   subroutine read_energy_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      logical :: enable_energy
      logical :: enable_cantera_thermo
      real(rk) :: initial_T
      real(rk) :: energy_reference_T
      real(rk) :: energy_reference_h
      real(rk) :: energy_cp
      real(rk) :: energy_lambda
      integer :: unit_id, ios

      namelist /energy_input/ enable_energy, enable_cantera_thermo, &
                               initial_T, energy_reference_T, energy_reference_h, &
                               energy_cp, energy_lambda

      enable_energy = params%enable_energy
      enable_cantera_thermo = params%enable_cantera_thermo
      initial_T = params%initial_T
      energy_reference_T = params%energy_reference_T
      energy_reference_h = params%energy_reference_h
      energy_cp = params%energy_cp
      energy_lambda = params%energy_lambda

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=energy_input, iostat=ios)
         close(unit_id)
      end if

      if (ios /= 0 .and. ios /= -1) then
         call fatal_error('input', 'failed reading &energy_input. Check for unknown variables or typos.')
      end if

      if (ios == 0) then
         params%enable_energy = enable_energy
         params%enable_cantera_thermo = enable_cantera_thermo
         params%initial_T = initial_T
         params%energy_reference_T = energy_reference_T
         params%energy_reference_h = energy_reference_h
         params%energy_cp = energy_cp
         params%energy_lambda = energy_lambda
      end if
   end subroutine read_energy_input


"""
    if "subroutine read_energy_input" not in text:
        text = insert_after_once(text, sub_anchor, read_energy_sub, "read_energy_input insertion")

    validate_anchor = (
        "      if (params%pressure_max_iter <= 0) call fatal_error('input', 'pressure_max_iter must be positive')\n"
        "      if (params%pressure_tol <= zero) call fatal_error('input', 'pressure_tol must be positive')\n"
    )
    validate_insert = (
        "\n"
        "      if (params%initial_T <= zero) call fatal_error('input', 'initial_T must be positive')\n"
        "      if (params%energy_reference_T <= zero) call fatal_error('input', 'energy_reference_T must be positive')\n"
        "      if (params%energy_cp <= zero) call fatal_error('input', 'energy_cp must be positive')\n"
        "      if (params%energy_lambda < zero) call fatal_error('input', 'energy_lambda must be non-negative')\n"
    )
    if "initial_T must be positive" not in text:
        text = insert_after_once(text, validate_anchor, validate_insert, "energy validation")

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_main(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "main.f90"
    text = read_text(path)

    if "use mod_energy, only : energy_fields_t" in text:
        return PatchResult(path, False, "already patched")

    use_anchor = (
        "   use mod_species, only : species_fields_t, initialize_species, finalize_species, &\n"
        "                           advance_species_transport\n"
    )
    use_insert = "   use mod_energy, only : energy_fields_t, initialize_energy, finalize_energy\n"
    text = insert_after_once(text, use_anchor, use_insert, "main mod_energy use")

    decl_anchor = "   type(species_fields_t) :: species !< Species mass fractions.\n"
    decl_insert = "   type(energy_fields_t) :: energy   !< Enthalpy/temperature/radiation-source fields.\n"
    text = insert_after_once(text, decl_anchor, decl_insert, "main energy declaration")

    init_anchor = (
        "   if (params%enable_species) then\n"
        "      call initialize_species(mesh, params, species)\n"
        "   end if\n"
    )
    init_insert = (
        "   if (params%enable_energy) then\n"
        "      call initialize_energy(mesh, params, energy)\n"
        "   end if\n"
    )
    text = insert_after_once(text, init_anchor, init_insert, "main energy initialization")

    final_anchor = (
        "   if (params%enable_species) then\n"
        "      call finalize_species(species)\n"
        "   end if\n"
    )
    final_insert = (
        "   if (params%enable_energy) then\n"
        "      call finalize_energy(energy)\n"
        "   end if\n"
    )
    text = insert_after_once(text, final_anchor, final_insert, "main energy finalization")

    return write_if_changed(path, text, repo, backup_root, apply)


def add_mod_energy(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_energy.f90"
    if path.exists():
        current = path.read_text(encoding="utf-8")
        if "module mod_energy" in current and "energy_fields_t" in current:
            return PatchResult(path, False, "already exists")
    return write_if_changed(path, MOD_ENERGY, repo, backup_root, apply, must_exist=False)


def patch_build_lists(repo: Path, backup_root: Path, apply: bool) -> list[PatchResult]:
    """Patch this repository's explicit Makefile only.

    The project Makefile has an ordered F_SRCS list plus explicit Fortran
    module dependencies. Do not use generic build-file probing here.
    """
    path = repo / "Makefile"
    if not path.exists():
        return [PatchResult(path, False, "Makefile not found; build list not changed")]

    text = path.read_text(encoding="utf-8")
    original = text

    src_line = "  $(SRC_DIR)/mod_energy.f90 \\\n"
    src_anchor = "  $(SRC_DIR)/mod_fields.f90 \\\n"
    if src_line not in text:
        if src_anchor not in text:
            raise SystemExit("[ERROR] Makefile source anchor not found: mod_fields.f90")
        text = text.replace(src_anchor, src_anchor + src_line, 1)

    dep_line = "$(BUILD_DIR)/mod_energy.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o\n"
    dep_anchor = "$(BUILD_DIR)/mod_fields.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_bc.o\n"
    if dep_line not in text:
        if dep_anchor not in text:
            raise SystemExit("[ERROR] Makefile dependency anchor not found: mod_fields.o rule")
        text = text.replace(dep_anchor, dep_anchor + dep_line, 1)

    main_old = "$(BUILD_DIR)/main.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mesh_io.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_mpi_radiation.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_transport_properties.o $(BUILD_DIR)/mod_species.o $(BUILD_DIR)/mod_output.o $(BUILD_DIR)/mod_profiler.o\n"
    main_new = "$(BUILD_DIR)/main.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_mesh_io.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_mpi_radiation.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o $(BUILD_DIR)/mod_energy.o $(BUILD_DIR)/mod_flow_projection.o $(BUILD_DIR)/mod_transport_properties.o $(BUILD_DIR)/mod_species.o $(BUILD_DIR)/mod_output.o $(BUILD_DIR)/mod_profiler.o\n"
    if "$(BUILD_DIR)/mod_energy.o" not in next((line for line in text.splitlines(True) if line.startswith("$(BUILD_DIR)/main.o:")), ""):
        if main_old not in text:
            raise SystemExit("[ERROR] Makefile main.o dependency anchor not found or already customized")
        text = text.replace(main_old, main_new, 1)

    if text == original:
        return [PatchResult(path, False, "Makefile already patched")]

    return [write_if_changed(path, text, repo, backup_root, apply)]


def sanity_check_repo(repo: Path) -> None:
    required = [
        repo / "src" / "main.f90",
        repo / "src" / "mod_input.f90",
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

    results.append(add_mod_energy(repo, backup_root, apply))
    results.append(patch_mod_input(repo, backup_root, apply))
    results.append(patch_main(repo, backup_root, apply))
    if args.skip_build_file:
        results.append(PatchResult(repo / "Makefile", False, "skipped by --skip-build-file"))
    else:
        results.extend(patch_build_lists(repo, backup_root, apply))

    print("\nPatch 001 summary")
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
        print("  2. build/compile")
        print("  3. run a baseline case with enable_energy = .false.")
        print("  4. optionally add &energy_input with enable_energy = .true. and run initialization-only check")
    else:
        print("\nNo changes were needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
