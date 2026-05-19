#!/usr/bin/env python3
"""
Patch 003: add Stage-1 constant-property enthalpy transport.

This patch advances the existing energy scaffold from field/output support to
actual constant-property passive enthalpy transport:

  - Adds patch_temperature_type(:) and patch_T(:) boundary inputs.
  - Adds temperature boundary-condition storage/evaluation in mod_bc.
  - Replaces mod_energy with a Stage-1 implementation that advances h using:

        V dh/dt = - sum_f F_f h_f
                 + (1/rho) sum_f lambda A_f (T_nb - T_c) / d_f
                 + (qrad/rho) V

    where h is upwinded for advection and diffusion uses grad(T).
  - Recovers T from h using the temporary constant-cp model.
  - Adds max_delta_T and rel_h_residual to energy diagnostics.
  - Wires advance_energy_transport into main after species transport.

This patch intentionally does NOT add Cantera thermodynamic inversion,
reactions, variable-density projection, or radiation coupling beyond the
already-existing qrad(:) source array.

Run from the repository root:

  python tools/patches/003_add_constant_property_energy_transport.py --dry-run
  python tools/patches/003_add_constant_property_energy_transport.py --apply

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


MOD_ENERGY = r'''!> Energy field storage, diagnostics, and Stage-1 constant-property transport.
!!
!! This module owns the enthalpy/temperature/radiation-source fields used by
!! the staged energy-equation implementation.
!!
!! Patch 003 implements Stage 1 only:
!! - h is transported as a passive constant-property scalar.
!! - T is recovered from h using a constant-cp relation.
!! - thermal diffusion uses grad(T), not grad(h).
!! - qrad is included as an explicit volumetric source, defaulting to zero.
!!
!! It intentionally does NOT add Cantera h/T inversion, chemical heat release,
!! species enthalpy diffusion, or variable-density low-Mach coupling.
module mod_energy
   use mpi_f08
   use mod_kinds, only : rk, zero, tiny_safe, fatal_error
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_mpi_flow, only : flow_mpi_t, flow_exchange_cell_scalar
   use mod_bc, only : bc_set_t, bc_periodic, face_effective_neighbor, &
                      patch_type_for_face, boundary_temperature
   use mod_fields, only : flow_fields_t
   implicit none

   private

   !> Cell-centered energy variables.
   type, public :: energy_fields_t
      real(rk), allocatable :: T(:)     !< Temperature [K].
      real(rk), allocatable :: h(:)     !< Mixture sensible enthalpy [J/kg].
      real(rk), allocatable :: h_old(:) !< Previous enthalpy state [J/kg].
      real(rk), allocatable :: qrad(:)  !< Volumetric radiation/source term [W/m^3].
      logical :: initialized = .false.  !< True after allocation/initialization.
   end type energy_fields_t

   public :: allocate_energy, initialize_energy, finalize_energy
   public :: update_enthalpy_from_temperature_constant_cp
   public :: recover_temperature_constant_cp
   public :: zero_radiation_source
   public :: advance_energy_transport
   public :: write_energy_diagnostics_header, write_energy_diagnostics_row

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


   !> Constant-cp helper for a single boundary/face temperature value.
   pure function enthalpy_from_temperature_value(params, temperature) result(h_value)
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: temperature
      real(rk) :: h_value

      h_value = params%energy_reference_h + &
                params%energy_cp * (temperature - params%energy_reference_T)
   end function enthalpy_from_temperature_value


   !> Reset the radiation/source term to zero.
   subroutine zero_radiation_source(energy)
      type(energy_fields_t), intent(inout) :: energy

      if (allocated(energy%qrad)) energy%qrad = zero
   end subroutine zero_radiation_source


   !> Advance passive sensible enthalpy with constant rho, cp, and lambda.
   !!
   !! The explicit finite-volume update is:
   !!
   !!   V dh/dt = - sum_f F_f h_f
   !!            + (1/rho) sum_f lambda A_f (T_nb - T_c)/d_f
   !!            + (qrad/rho) V
   !!
   !! fields%face_flux is volumetric flux. It is re-oriented outward from the
   !! currently updated cell before applying upwind h advection.
   subroutine advance_energy_transport(mesh, flow, bc, params, fields, energy)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(inout) :: energy

      integer :: c, lf, fid, neighbor
      real(rk) :: flux_out
      real(rk) :: face_area, dist
      real(rk) :: h_cell, h_other, h_face
      real(rk) :: T_cell, T_other
      real(rk) :: rhs_h, diff_term
      logical :: is_dirichlet, do_diffusion

      if (.not. params%enable_energy) return

      if (.not. energy%initialized) then
         call fatal_error('energy', 'advance requested before energy initialization')
      end if

      if (params%enable_cantera_thermo) then
         call fatal_error('energy', 'Patch 003 only supports constant-cp energy transport; disable enable_cantera_thermo')
      end if

      if (params%rho <= tiny_safe) call fatal_error('energy', 'rho must be positive')
      if (params%energy_cp <= tiny_safe) call fatal_error('energy', 'energy_cp must be positive')
      if (params%energy_lambda < zero) call fatal_error('energy', 'energy_lambda must be non-negative')

      ! Ensure off-rank neighbor values are current before taking the old state.
      call flow_exchange_cell_scalar(flow, energy%h)
      call flow_exchange_cell_scalar(flow, energy%T)
      energy%h_old = energy%h

      do c = flow%first_cell, flow%last_cell
         rhs_h = zero
         h_cell = energy%h_old(c)
         T_cell = energy%T(c)

         do lf = 1, mesh%ncell_faces(c)
            fid = mesh%cell_faces(lf, c)

            ! fields%face_flux is oriented owner -> neighbor. Reorient outward
            ! from the current cell.
            if (mesh%faces(fid)%owner == c) then
               flux_out = fields%face_flux(fid)
            else
               flux_out = -fields%face_flux(fid)
            end if

            neighbor = face_effective_neighbor(mesh, bc, fid, c)
            face_area = mesh%faces(fid)%area

            if (neighbor > 0) then
               h_other = energy%h_old(neighbor)
               T_other = energy%T(neighbor)
               do_diffusion = .true.
            else
               call boundary_temperature(mesh, bc, fid, T_cell, T_other, is_dirichlet)
               if (is_dirichlet) then
                  h_other = enthalpy_from_temperature_value(params, T_other)
                  do_diffusion = .true.
               else
                  h_other = h_cell
                  T_other = T_cell
                  do_diffusion = .false.
               end if
            end if

            ! Upwind advection of h using outward volumetric flux.
            if (flux_out >= zero) then
               h_face = h_cell
            else
               h_face = h_other
            end if
            rhs_h = rhs_h - flux_out * h_face

            ! Fourier heat conduction uses grad(T), not grad(h).
            if (do_diffusion .and. params%energy_lambda > zero) then
               dist = energy_face_normal_distance(mesh, bc, fid, c, neighbor)
               diff_term = (params%energy_lambda / params%rho) * &
                           (T_other - T_cell) / dist * face_area
               rhs_h = rhs_h + diff_term
            end if
         end do

         energy%h(c) = energy%h_old(c) + params%dt * &
                       (rhs_h / mesh%cells(c)%volume + energy%qrad(c) / params%rho)
      end do

      call recover_temperature_constant_cp(params, energy)

      ! Synchronize updated owned cells for output and the next step.
      call flow_exchange_cell_scalar(flow, energy%h)
      call flow_exchange_cell_scalar(flow, energy%T)
   end subroutine advance_energy_transport


   !> Outward unit normal from cell_id for face_id.
   function energy_outward_normal(mesh, face_id, cell_id) result(nvec)
      type(mesh_t), intent(in) :: mesh
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      real(rk) :: nvec(3)

      if (mesh%faces(face_id)%owner == cell_id) then
         nvec = mesh%faces(face_id)%normal
      else
         nvec = -mesh%faces(face_id)%normal
      end if
   end function energy_outward_normal


   !> Normal distance used for cell-cell and cell-boundary temperature gradients.
   function energy_face_normal_distance(mesh, bc, face_id, cell_id, nb) result(dist)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      integer, intent(in) :: nb
      real(rk) :: dist

      integer :: pair_face
      integer :: btype
      real(rk) :: nvec(3)

      nvec = energy_outward_normal(mesh, face_id, cell_id)

      if (nb > 0) then
         if (mesh%faces(face_id)%neighbor == 0) then
            btype = patch_type_for_face(mesh, bc, face_id)

            if (btype == bc_periodic) then
               pair_face = mesh%faces(face_id)%periodic_face

               if (pair_face <= 0) then
                  call fatal_error('energy', 'periodic face has no paired face')
               end if

               dist = abs(dot_product(mesh%faces(face_id)%center - &
                                      mesh%cells(cell_id)%center, nvec)) + &
                      abs(dot_product(mesh%cells(nb)%center - &
                                      mesh%faces(pair_face)%center, nvec))

               dist = max(dist, tiny_safe)
               return
            end if
         end if

         dist = abs(dot_product(mesh%cells(nb)%center - &
                                mesh%cells(cell_id)%center, nvec))
      else
         dist = abs(dot_product(mesh%faces(face_id)%center - &
                                mesh%cells(cell_id)%center, nvec))
      end if

      dist = max(dist, tiny_safe)
   end function energy_face_normal_distance


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

      write(unit_id,'(a)') 'step,time,min_T,max_T,mean_T,min_h,max_h,mean_h,min_qrad,max_qrad,integral_qrad,max_delta_T,rel_h_residual'

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
      real(rk) :: vol, old_T
      real(rk) :: local_min_T, local_max_T, local_sum_T
      real(rk) :: local_min_h, local_max_h, local_sum_h
      real(rk) :: local_min_qrad, local_max_qrad, local_integral_qrad
      real(rk) :: local_volume
      real(rk) :: local_max_delta_T, local_delta_h2, local_h_old2
      real(rk) :: global_min_T, global_max_T, global_sum_T
      real(rk) :: global_min_h, global_max_h, global_sum_h
      real(rk) :: global_min_qrad, global_max_qrad, global_integral_qrad
      real(rk) :: global_volume
      real(rk) :: global_max_delta_T, global_delta_h2, global_h_old2
      real(rk) :: mean_T, mean_h, rel_h_residual

      if (.not. params%write_diagnostics) return
      if (.not. params%enable_energy) return

      if (.not. allocated(energy%T) .or. .not. allocated(energy%h) .or. &
          .not. allocated(energy%h_old) .or. .not. allocated(energy%qrad)) then
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
      local_max_delta_T = zero
      local_delta_h2 = zero
      local_h_old2 = zero

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

         old_T = params%energy_reference_T + &
                 (energy%h_old(c) - params%energy_reference_h) / params%energy_cp
         local_max_delta_T = max(local_max_delta_T, abs(energy%T(c) - old_T))
         local_delta_h2 = local_delta_h2 + (energy%h(c) - energy%h_old(c))**2 * vol
         local_h_old2 = local_h_old2 + energy%h_old(c)**2 * vol
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
      call MPI_Allreduce(local_max_delta_T, global_max_delta_T, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing max_delta_T')
      call MPI_Allreduce(local_delta_h2, global_delta_h2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing delta_h2')
      call MPI_Allreduce(local_h_old2, global_h_old2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing h_old2')

      if (global_volume > tiny_safe) then
         mean_T = global_sum_T / global_volume
         mean_h = global_sum_h / global_volume
      else
         mean_T = zero
         mean_h = zero
      end if

      rel_h_residual = sqrt(global_delta_h2) / max(sqrt(global_h_old2), tiny_safe)

      if (flow%rank /= 0) return

      filename = trim(params%output_dir)//'/energy_diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='old', position='append', action='write')

      write(unit_id,'(i0,12(",",es16.8))') step, time, global_min_T, global_max_T, mean_T, &
         global_min_h, global_max_h, mean_h, global_min_qrad, global_max_qrad, &
         global_integral_qrad, global_max_delta_T, rel_h_residual

      close(unit_id)
   end subroutine write_energy_diagnostics_row

end module mod_energy
'''


@dataclass
class PatchResult:
    path: Path
    changed: bool
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Patch 003: add Stage-1 constant-property energy transport.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true", help="write changes to disk")
    mode.add_argument("--dry-run", action="store_true", help="show planned changes without writing")
    parser.add_argument("--repo", default=".", help="repository root; default: current directory")
    parser.add_argument("--skip-build-file", action="store_true", help="do not modify Makefile dependencies")
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


def insert_after_once(text: str, anchor: str, insertion: str, label: str) -> str:
    count = text.count(anchor)
    if count == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if count > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {count} matches")
    return text.replace(anchor, anchor + insertion, 1)


def replace_once(text: str, old: str, new: str, label: str) -> str:
    count = text.count(old)
    if count == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if count > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {count} matches")
    return text.replace(old, new, 1)


def patch_mod_energy(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_energy.f90"
    return write_if_changed(path, MOD_ENERGY, repo, backup_root, apply)


def patch_mod_input(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_input.f90"
    text = read_text(path)

    if "patch_temperature_type" not in text:
        text = insert_after_once(
            text,
            '      character(len=name_len) :: patch_pressure_type(max_patches) = "" !< Override BC type for pressure.\n',
            '      character(len=name_len) :: patch_temperature_type(max_patches) = "" !< Override BC type for temperature/enthalpy.\n',
            "case_params_t patch_temperature_type field",
        )

    if "patch_T(max_patches)" not in text:
        text = insert_after_once(
            text,
            '      real(rk) :: patch_dpdn(max_patches) = zero !< Specified pressure gradient on patch [Pa/m].\n',
            '      real(rk) :: patch_T(max_patches) = 300.0_rk !< Specified temperature on patch [K].\n',
            "case_params_t patch_T field",
        )

    if '      character(len=name_len) :: patch_temperature_type(max_patches)\n' not in text:
        text = insert_after_once(
            text,
            '      character(len=name_len) :: patch_pressure_type(max_patches)\n',
            '      character(len=name_len) :: patch_temperature_type(max_patches)\n',
            "read_boundary_input local patch_temperature_type",
        )

    if '      real(rk) :: patch_T(max_patches)\n' not in text:
        text = insert_after_once(
            text,
            '      real(rk) :: patch_dpdn(max_patches)\n',
            '      real(rk) :: patch_T(max_patches)\n',
            "read_boundary_input local patch_T",
        )

    if "patch_temperature_type, patch_species_type" not in text:
        text = replace_once(
            text,
            '                                 patch_velocity_type, patch_pressure_type, &\n'
            '                                 patch_species_type, &\n',
            '                                 patch_velocity_type, patch_pressure_type, &\n'
            '                                 patch_temperature_type, patch_species_type, &\n',
            "boundary namelist temperature type",
        )

    if "patch_dpdn, patch_T, patch_Y" not in text:
        text = replace_once(
            text,
            '                                 patch_u, patch_v, patch_w, patch_p, patch_dpdn, patch_Y\n',
            '                                 patch_u, patch_v, patch_w, patch_p, patch_dpdn, patch_T, patch_Y\n',
            "boundary namelist patch_T",
        )

    if "patch_temperature_type = params%patch_temperature_type" not in text:
        text = insert_after_once(
            text,
            '      patch_pressure_type = params%patch_pressure_type\n',
            '      patch_temperature_type = params%patch_temperature_type\n',
            "read_boundary_input initialize patch_temperature_type",
        )

    if "patch_T = params%patch_T" not in text:
        text = insert_after_once(
            text,
            '      patch_dpdn = params%patch_dpdn\n',
            '      patch_T = params%patch_T\n',
            "read_boundary_input initialize patch_T",
        )

    if "params%patch_temperature_type = patch_temperature_type" not in text:
        text = insert_after_once(
            text,
            '      params%patch_pressure_type = patch_pressure_type\n',
            '      params%patch_temperature_type = patch_temperature_type\n',
            "read_boundary_input assign patch_temperature_type",
        )

    if "params%patch_temperature_type(i) = trim(lowercase(params%patch_temperature_type(i)))" not in text:
        text = insert_after_once(
            text,
            '         params%patch_pressure_type(i) = trim(lowercase(params%patch_pressure_type(i)))\n',
            '         params%patch_temperature_type(i) = trim(lowercase(params%patch_temperature_type(i)))\n',
            "read_boundary_input lowercase patch_temperature_type",
        )

    if "params%patch_T = patch_T" not in text:
        text = insert_after_once(
            text,
            '      params%patch_dpdn = patch_dpdn\n',
            '      params%patch_T = patch_T\n',
            "read_boundary_input assign patch_T",
        )

    if "patch_T entry must be positive" not in text:
        text = insert_after_once(
            text,
            '         if (len_trim(params%patch_type(i)) == 0) then\n'
            "            call fatal_error('input', 'patch_type entry cannot be empty')\n"
            '         end if\n',
            '\n'
            '         if (params%enable_energy .and. params%patch_T(i) <= zero) then\n'
            "            call fatal_error('input', 'patch_T entry must be positive when energy is enabled')\n"
            '         end if\n',
            "validate_boundary_arrays patch_T",
        )

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_mod_bc(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_bc.f90"
    text = read_text(path)

    if "temperature_type_id" not in text:
        text = insert_after_once(
            text,
            '      integer :: pressure_type_id = bc_unknown   !< BC type for pressure (if different from master).\n',
            '\n'
            '      !> Temperature/enthalpy boundary settings.\n'
            '      integer :: temperature_type_id = bc_unknown !< BC type applied to temperature/enthalpy transport.\n'
            '      real(rk) :: temperature = 300.0_rk          !< Specified boundary temperature [K].\n',
            "bc_patch_t temperature fields",
        )

    if "public :: boundary_temperature" not in text:
        text = insert_after_once(
            text,
            '   public :: boundary_species\n',
            '   public :: boundary_temperature\n',
            "mod_bc public boundary_temperature",
        )

    if "bc%patches(p)%temperature = params%patch_T(q)" not in text:
        text = insert_after_once(
            text,
            '               bc%patches(p)%dpdn = params%patch_dpdn(q)\n',
            '               bc%patches(p)%temperature = params%patch_T(q)\n',
            "build_bc_set patch temperature value",
        )

    if "params%patch_temperature_type(q)" not in text:
        text = insert_after_once(
            text,
            '               if (trim(params%patch_species_type(q)) /= "") then\n'
            '                  bc%patches(p)%species_type_id = parse_bc_type(params%patch_species_type(q))\n'
            '               else\n'
            '                  bc%patches(p)%species_type_id = parse_bc_type(params%patch_type(q))\n'
            '               end if\n'
            '               bc%patches(p)%species_Y(:) = params%patch_Y(:, q)\n',
            '\n'
            '               if (trim(params%patch_temperature_type(q)) /= "") then\n'
            '                  bc%patches(p)%temperature_type_id = parse_bc_type(params%patch_temperature_type(q))\n'
            '               else\n'
            '                  bc%patches(p)%temperature_type_id = parse_bc_type(params%patch_type(q))\n'
            '               end if\n',
            "build_bc_set temperature type",
        )

    if "subroutine boundary_temperature" not in text:
        insertion = r'''
   !> Evaluates temperature at a boundary face.
   !!
   !! Dirichlet/fixed_value boundaries return the configured patch_T.
   !! All other boundary types currently behave as zero-gradient boundaries.
   subroutine boundary_temperature(mesh, bc, face_id, interior_T, ext_T, is_dirichlet)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      real(rk), intent(in) :: interior_T
      real(rk), intent(out) :: ext_T
      logical, intent(out) :: is_dirichlet

      integer :: patch_id

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         ext_T = interior_T
         is_dirichlet = .false.
         return
      end if

      select case (bc%patches(patch_id)%temperature_type_id)
      case (bc_dirichlet)
         ext_T = bc%patches(patch_id)%temperature
         is_dirichlet = .true.
      case default
         ext_T = interior_T
         is_dirichlet = .false.
      end select
   end subroutine boundary_temperature

'''
        text = insert_after_once(
            text,
            '   !> Returns the pressure BC type for a given face.\n',
            insertion,
            "boundary_temperature insertion",
        )

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_main(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "main.f90"
    text = read_text(path)

    if "advance_energy_transport" not in text:
        text = replace_once(
            text,
            '   use mod_energy, only : energy_fields_t, initialize_energy, finalize_energy, &\n'
            '                           write_energy_diagnostics_header, write_energy_diagnostics_row\n',
            '   use mod_energy, only : energy_fields_t, initialize_energy, finalize_energy, &\n'
            '                           advance_energy_transport, &\n'
            '                           write_energy_diagnostics_header, write_energy_diagnostics_row\n',
            "main mod_energy use advance_energy_transport",
        )

    if "call advance_energy_transport" not in text:
        text = insert_after_once(
            text,
            '      if (params%enable_species) then\n'
            "         call profiler_start('Flow_Transport')\n"
            '         call advance_species_transport(mesh, flow_mpi, bc, params, fields, species, transport)\n'
            "         call profiler_stop('Flow_Transport')\n"
            '      end if\n',
            '      \n'
            '      ! E. Advance constant-property enthalpy transport.\n'
            '      if (params%enable_energy) then\n'
            "         call profiler_start('Energy_Transport')\n"
            '         call advance_energy_transport(mesh, flow_mpi, bc, params, fields, energy)\n'
            "         call profiler_stop('Energy_Transport')\n"
            '      end if\n',
            "main energy transport call",
        )

    return write_if_changed(path, text, repo, backup_root, apply)


def patch_makefile(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "Makefile"
    if not path.exists():
        return PatchResult(path, False, "Makefile not found; skipped")

    text = read_text(path)
    old = '$(BUILD_DIR)/mod_energy.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mpi_flow.o\n'
    new = '$(BUILD_DIR)/mod_energy.o: $(BUILD_DIR)/mod_kinds.o $(BUILD_DIR)/mod_mesh_types.o $(BUILD_DIR)/mod_input.o $(BUILD_DIR)/mod_mpi_flow.o $(BUILD_DIR)/mod_bc.o $(BUILD_DIR)/mod_fields.o\n'

    if new in text:
        return PatchResult(path, False, "Makefile dependency already patched")
    if old in text:
        text = text.replace(old, new, 1)
        return write_if_changed(path, text, repo, backup_root, apply)

    return PatchResult(path, False, "mod_energy dependency anchor not found; Makefile not changed")


def sanity_check_repo(repo: Path) -> None:
    required = [
        repo / "src" / "main.f90",
        repo / "src" / "mod_input.f90",
        repo / "src" / "mod_bc.f90",
        repo / "src" / "mod_energy.f90",
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

    results.append(patch_mod_input(repo, backup_root, apply))
    results.append(patch_mod_bc(repo, backup_root, apply))
    results.append(patch_mod_energy(repo, backup_root, apply))
    results.append(patch_main(repo, backup_root, apply))
    if not args.skip_build_file:
        results.append(patch_makefile(repo, backup_root, apply))

    print("\nPatch 003 summary")
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
        print("  2. add patch_temperature_type and patch_T to the case boundary_input")
        print("  3. build/compile")
        print("  4. run a diffusion-only sanity case, then the rectangle_2D advection case")
    else:
        print("\nNo changes were needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
