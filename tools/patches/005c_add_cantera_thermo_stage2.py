#!/usr/bin/env python3
"""
Patch 005c: add Stage-2A Cantera thermodynamics without reactions.

This corrected patch is based on the current Stage-1 source layout. It fixes
Patch 005b's fragile lambda-diffusion anchor by matching the actual current
form:

    diff_term = (params%energy_lambda / params%rho) * &

Stage 2A scope:
  - Cantera initializes when enable_cantera_thermo = .true.
  - no-species smoke tests use N2 as the inert thermodynamic mixture.
  - energy h(T,p,Y), T(h,p,Y), cp, lambda, and rho_thermo come from Cantera.
  - the flow projection remains constant-density and unchanged.
  - reactions and chemical heat release remain disabled.

Run from repository root:

  python tools/patches/005c_add_cantera_thermo_stage2.py --dry-run
  python tools/patches/005c_add_cantera_thermo_stage2.py --apply
"""

from __future__ import annotations

import argparse
import datetime as _dt
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class PatchResult:
    path: Path
    changed: bool
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Patch 005c: add Stage-2A Cantera thermodynamics.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true", help="write changes to disk")
    mode.add_argument("--dry-run", action="store_true", help="show planned changes without writing")
    parser.add_argument("--repo", default=".", help="repository root; default: current directory")
    return parser.parse_args()


def timestamp() -> str:
    return _dt.datetime.now().strftime("%Y%m%d_%H%M%S")


def read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise SystemExit(f"[ERROR] required file not found: {path}")


def make_backup(path: Path, repo: Path, backup_root: Path) -> None:
    rel = path.relative_to(repo)
    dst = backup_root / rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, dst)


def write_if_changed(path: Path, text: str, repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    old = read_text(path) if path.exists() else ""
    if old == text:
        return PatchResult(path, False, "unchanged")
    if apply:
        make_backup(path, repo, backup_root)
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


def insert_before_once(text: str, anchor: str, insertion: str, label: str) -> str:
    n = text.count(anchor)
    if n == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if n > 1:
        raise SystemExit(f"[ERROR] anchor is ambiguous for {label}: found {n} matches")
    return text.replace(anchor, insertion + anchor, 1)


def replace_first(text: str, old: str, new: str, label: str) -> tuple[str, bool]:
    idx = text.find(old)
    if idx < 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    if text[idx:idx + len(new)] == new:
        return text, False
    return text[:idx] + new + text[idx + len(old):], True


def patch_cantera_interface(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "cantera_interface.cpp"
    text = read_text(path)
    changed = False

    if "build_cantera_mass_fractions" not in text:
        helper = r'''
// Build a full Cantera mass-fraction vector from the solver species list.
// If the solver list is empty, or if the sum is below one, the remaining
// mass fraction is assigned to N2 when present. This is the Stage-2A fallback
// for non-reacting no-species thermo smoke tests.
static void build_cantera_mass_fractions(
    Cantera::ThermoPhase& gas_phase,
    int nspecies,
    const double* Y_in,
    int cell_index,
    const std::vector<int>& sp_map,
    std::vector<double>& Y_cantera)
{
    std::fill(Y_cantera.begin(), Y_cantera.end(), 0.0);

    int bath_gas_index = (int)gas_phase.speciesIndex("N2");
    double sum_Y = 0.0;

    for (int k = 0; k < nspecies; ++k) {
        if (sp_map[k] >= 0) {
            double y_val = Y_in[cell_index * nspecies + k];
            if (y_val < 0.0) {
                y_val = 0.0;
            }
            Y_cantera[sp_map[k]] = y_val;
            sum_Y += y_val;
        }
    }

    if (bath_gas_index >= 0) {
        Y_cantera[bath_gas_index] += std::max(0.0, 1.0 - sum_Y);
    }

    double total = 0.0;
    for (double y : Y_cantera) {
        total += y;
    }

    if (total <= 0.0) {
        if (bath_gas_index >= 0) {
            Y_cantera[bath_gas_index] = 1.0;
        } else if (!Y_cantera.empty()) {
            Y_cantera[0] = 1.0;
        }
    }
}

// Map the solver species names to Cantera species indices.
static std::vector<int> build_species_map(
    Cantera::ThermoPhase& gas_phase,
    int nspecies,
    const char* species_names_flat,
    int name_len)
{
    std::vector<int> sp_map(nspecies, -1);

    for (int k = 0; k < nspecies; ++k) {
        std::string sp_name(species_names_flat + k * name_len, name_len);
        size_t last = sp_name.find_last_not_of(" ");
        if (last == std::string::npos) {
            sp_name = "";
        } else {
            sp_name.erase(last + 1);
        }

        size_t c_idx = gas_phase.speciesIndex(sp_name);
        if (c_idx != Cantera::npos) {
            sp_map[k] = (int)c_idx;
        } else {
            std::cerr << "Cantera Bridge: Species " << sp_name
                      << " not found in mechanism!" << std::endl;
        }
    }

    return sp_map;
}

'''
        text = insert_after_once(text, "static std::shared_ptr<Cantera::Transport> trans;\n", helper,
                                 "Cantera helper functions")
        changed = True

    if "cantera_update_thermo_c" not in text:
        funcs = r'''
    // Compute thermodynamic properties from T, P, and composition.
    void cantera_update_thermo_c(int ncells, double* T, double* P, int nspecies, double* Y_in,
                                 double* h_out, double* cp_out, double* lambda_out,
                                 double* rho_thermo_out,
                                 const char* species_names_flat, int name_len) {
        if (!gas || !trans) {
            std::cerr << "Cantera Bridge: Gas not initialized before update_thermo!" << std::endl;
            exit(1);
        }

        try {
            int cantera_nsp = (int)gas->nSpecies();
            std::vector<double> Y_cantera(cantera_nsp, 0.0);
            std::vector<int> sp_map = build_species_map(*gas, nspecies, species_names_flat, name_len);

            for (int c = 0; c < ncells; ++c) {
                build_cantera_mass_fractions(*gas, nspecies, Y_in, c, sp_map, Y_cantera);
                gas->setState_TPY(T[c], P[c], Y_cantera.data());

                h_out[c] = gas->enthalpy_mass();
                cp_out[c] = gas->cp_mass();
                rho_thermo_out[c] = gas->density();
                lambda_out[c] = trans->thermalConductivity();
            }

        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        } catch (std::exception& e) {
            std::cerr << "Standard Exception: " << e.what() << std::endl;
            exit(1);
        }
    }


    // Recover temperature from mixture enthalpy, pressure, and composition.
    void cantera_recover_temperature_from_h_c(int ncells, double* h_in, double* P, int nspecies, double* Y_in,
                                              double* T_out,
                                              const char* species_names_flat, int name_len) {
        if (!gas) {
            std::cerr << "Cantera Bridge: Gas not initialized before recover_temperature!" << std::endl;
            exit(1);
        }

        try {
            int cantera_nsp = (int)gas->nSpecies();
            std::vector<double> Y_cantera(cantera_nsp, 0.0);
            std::vector<int> sp_map = build_species_map(*gas, nspecies, species_names_flat, name_len);

            for (int c = 0; c < ncells; ++c) {
                build_cantera_mass_fractions(*gas, nspecies, Y_in, c, sp_map, Y_cantera);
                gas->setState_HPY(h_in[c], P[c], Y_cantera.data());
                T_out[c] = gas->temperature();
            }

        } catch (Cantera::CanteraError& err) {
            std::cerr << "Cantera Error: " << err.what() << std::endl;
            exit(1);
        } catch (std::exception& e) {
            std::cerr << "Standard Exception: " << e.what() << std::endl;
            exit(1);
        }
    }

'''
        text = insert_before_once(text, "} // extern \"C\"\n", funcs, "Cantera thermo C interface functions")
        changed = True

    if not changed:
        return PatchResult(path, False, "already patched")
    return write_if_changed(path, text, repo, backup_root, apply)


def patch_mod_transport_properties(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_transport_properties.f90"
    text = read_text(path)
    changed = False

    old_cond = "      if (params%enable_cantera_fluid .or. params%enable_cantera_species) then\n"
    new_cond = "      if (params%enable_cantera_fluid .or. params%enable_cantera_species .or. params%enable_cantera_thermo) then\n"
    if new_cond not in text:
        text, did = replace_first(text, old_cond, new_cond, "initialize_transport Cantera thermo condition")
        changed = changed or did

    if "params%enable_cantera_thermo .and. params%nspecies == 0" not in text:
        old_block = (
            "      else if (params%enable_cantera_fluid .and. params%nspecies == 0) then\n"
            "         ! Single Species Fluid Mode: Default to the first species (e.g., N2)\n"
            "         params%nspecies = 1\n"
            "         n_len = len(params%species_name(1))\n"
            "         call cantera_get_species_name_c(0, params%species_name(1), n_len)\n"
            "      end if\n"
        )
        new_block = (
            "      else if (params%enable_cantera_thermo .and. params%nspecies == 0) then\n"
            "         ! Stage-2A no-species thermo mode: use N2 as inert bath gas.\n"
            "         params%nspecies = 1\n"
            "         params%species_name(1) = 'N2'\n"
            "      else if (params%enable_cantera_fluid .and. params%nspecies == 0) then\n"
            "         ! Single Species Fluid Mode: default to N2 when no species list is provided.\n"
            "         params%nspecies = 1\n"
            "         params%species_name(1) = 'N2'\n"
            "      end if\n"
        )
        if old_block in text:
            text = text.replace(old_block, new_block, 1)
            changed = True
        else:
            anchor = "      else if (params%enable_cantera_fluid .and. params%nspecies == 0) then\n"
            insertion = (
                "      else if (params%enable_cantera_thermo .and. params%nspecies == 0) then\n"
                "         ! Stage-2A no-species thermo mode: use N2 as inert bath gas.\n"
                "         params%nspecies = 1\n"
                "         params%species_name(1) = 'N2'\n"
            )
            idx = text.find(anchor)
            if idx < 0:
                raise SystemExit("[ERROR] anchor not found for initialize_cantera_wrapper N2 default")
            text = text[:idx] + insertion + text[idx:]
            changed = True

    if not changed:
        return PatchResult(path, False, "already patched")
    return write_if_changed(path, text, repo, backup_root, apply)


THERMO_INTERFACE = r'''   !> C-Binding interface for Stage-2A Cantera thermodynamics.
   interface
      subroutine cantera_update_thermo_c(ncells, T, P, nspecies, Y_in, h_out, cp_out, lambda_out, &
                                         rho_thermo_out, species_names_flat, name_len) bind(c, name="cantera_update_thermo_c")
         import :: c_char, c_double, c_int
         integer(c_int), value :: ncells
         real(c_double), intent(in) :: T(ncells), P(ncells)
         integer(c_int), value :: nspecies
         real(c_double), intent(in) :: Y_in(*)
         real(c_double), intent(out) :: h_out(ncells), cp_out(ncells), lambda_out(ncells), rho_thermo_out(ncells)
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_update_thermo_c

      subroutine cantera_recover_temperature_from_h_c(ncells, h_in, P, nspecies, Y_in, T_out, &
                                                      species_names_flat, name_len) bind(c, name="cantera_recover_temperature_from_h_c")
         import :: c_char, c_double, c_int
         integer(c_int), value :: ncells
         real(c_double), intent(in) :: h_in(ncells), P(ncells)
         integer(c_int), value :: nspecies
         real(c_double), intent(in) :: Y_in(*)
         real(c_double), intent(out) :: T_out(ncells)
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_recover_temperature_from_h_c
   end interface

'''

THERMO_HELPERS = r'''

   !> Build a flattened C-compatible species-name buffer.
   subroutine build_c_species_names(params, c_names_flat, n_len)
      type(case_params_t), intent(in) :: params
      character(kind=c_char), allocatable, intent(out) :: c_names_flat(:)
      integer, intent(out) :: n_len

      integer :: k, j, nsp

      nsp = max(1, params%nspecies)
      n_len = len(params%species_name(1))
      allocate(c_names_flat(n_len * nsp))
      c_names_flat = ' '

      do k = 1, nsp
         do j = 1, n_len
            c_names_flat((k-1)*n_len + j) = params%species_name(k)(j:j)
         end do
      end do
   end subroutine build_c_species_names


   !> Fill a local inert-mixture mass-fraction array for Stage-2A smoke tests.
   !!
   !! Stage 2A deliberately supports only the no-transported-species Cantera
   !! thermodynamics path. initialize_cantera_wrapper sets species_name(1)=N2
   !! for this path, so Y_local(1,:)=1 represents pure nitrogen.
   subroutine fill_stage2a_Y(params, ncells, Y_local)
      type(case_params_t), intent(in) :: params
      integer, intent(in) :: ncells
      real(rk), allocatable, intent(out) :: Y_local(:,:)

      allocate(Y_local(max(1, params%nspecies), ncells))
      Y_local = zero
      Y_local(1, :) = 1.0_rk
   end subroutine fill_stage2a_Y


   !> Update h, cp, lambda, and diagnostic thermo density from current T.
   subroutine update_thermo_from_temperature_cantera(mesh, params, energy)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy

      integer :: n_len
      real(rk), allocatable :: P_arr(:), Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (.not. params%enable_cantera_thermo) return

      if (params%nspecies <= 0) then
         call fatal_error('energy', 'enable_cantera_thermo requires at least one inert thermo species')
      end if

      allocate(P_arr(mesh%ncells))
      P_arr = params%background_press
      call fill_stage2a_Y(params, mesh%ncells, Y_local)
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_update_thermo_c(mesh%ncells, energy%T, P_arr, params%nspecies, Y_local, &
                                   energy%h, energy%cp, energy%lambda, energy%rho_thermo, &
                                   c_names_flat, n_len)

      deallocate(P_arr)
      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine update_thermo_from_temperature_cantera


   !> Recover T from h using Cantera HPY inversion.
   subroutine recover_temperature_cantera(mesh, params, energy)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy

      integer :: n_len
      real(rk), allocatable :: P_arr(:), Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (.not. params%enable_cantera_thermo) return

      if (params%nspecies <= 0) then
         call fatal_error('energy', 'enable_cantera_thermo requires at least one inert thermo species')
      end if

      allocate(P_arr(mesh%ncells))
      P_arr = params%background_press
      call fill_stage2a_Y(params, mesh%ncells, Y_local)
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_recover_temperature_from_h_c(mesh%ncells, energy%h, P_arr, params%nspecies, Y_local, &
                                                energy%T, c_names_flat, n_len)

      deallocate(P_arr)
      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine recover_temperature_cantera


   !> Convert a boundary temperature to enthalpy using the active thermo model.
   function boundary_enthalpy_from_temperature(params, temperature) result(h_value)
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: temperature
      real(rk) :: h_value

      integer :: n_len
      real(rk) :: T_one(1), P_one(1), h_one(1), cp_one(1), lambda_one(1), rho_one(1)
      real(rk), allocatable :: Y_one(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (.not. params%enable_cantera_thermo) then
         h_value = enthalpy_from_temperature_value(params, temperature)
         return
      end if

      T_one(1) = temperature
      P_one(1) = params%background_press
      call fill_stage2a_Y(params, 1, Y_one)
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_update_thermo_c(1, T_one, P_one, params%nspecies, Y_one, &
                                   h_one, cp_one, lambda_one, rho_one, c_names_flat, n_len)
      h_value = h_one(1)

      deallocate(Y_one)
      deallocate(c_names_flat)
   end function boundary_enthalpy_from_temperature
'''


def replace_regex_once(text: str, pattern: str, repl: str, label: str, flags: int = 0) -> str:
    new_text, n = re.subn(pattern, repl, text, count=1, flags=flags)
    if n == 0:
        raise SystemExit(f"[ERROR] anchor not found for {label}")
    return new_text


def patch_mod_energy(repo: Path, backup_root: Path, apply: bool) -> PatchResult:
    path = repo / "src" / "mod_energy.f90"
    text = read_text(path)
    changed = False

    if "use iso_c_binding, only : c_char, c_double, c_int" not in text:
        text = insert_after_once(text, "   use mpi_f08\n", "   use iso_c_binding, only : c_char, c_double, c_int\n",
                                 "mod_energy iso_c_binding import")
        changed = True

    if "T_old(:)" not in text:
        old = (
            "      real(rk), allocatable :: T(:)     !< Temperature [K].\n"
            "      real(rk), allocatable :: h(:)     !< Mixture sensible enthalpy [J/kg].\n"
            "      real(rk), allocatable :: h_old(:) !< Previous enthalpy state [J/kg].\n"
            "      real(rk), allocatable :: qrad(:)  !< Volumetric radiation/source term [W/m^3].\n"
        )
        new = (
            "      real(rk), allocatable :: T(:)          !< Temperature [K].\n"
            "      real(rk), allocatable :: T_old(:)      !< Previous temperature state [K].\n"
            "      real(rk), allocatable :: h(:)          !< Mixture sensible enthalpy [J/kg].\n"
            "      real(rk), allocatable :: h_old(:)      !< Previous enthalpy state [J/kg].\n"
            "      real(rk), allocatable :: qrad(:)       !< Volumetric radiation/source term [W/m^3].\n"
            "      real(rk), allocatable :: cp(:)         !< Mixture heat capacity at constant pressure [J/kg/K].\n"
            "      real(rk), allocatable :: lambda(:)     !< Thermal conductivity [W/m/K].\n"
            "      real(rk), allocatable :: rho_thermo(:) !< Cantera thermodynamic density [kg/m^3], diagnostic/future-use.\n"
        )
        text = replace_once(text, old, new, "energy_fields_t thermo arrays")
        changed = True

    if "cantera_update_thermo_c" not in text:
        text = insert_before_once(text, "contains\n", THERMO_INTERFACE, "mod_energy Cantera thermo interface")
        changed = True

    if "allocate(energy%T_old" not in text:
        text = insert_after_once(text, "      allocate(energy%T(mesh%ncells))\n", "      allocate(energy%T_old(mesh%ncells))\n",
                                 "allocate T_old")
        text = insert_after_once(text, "      allocate(energy%qrad(mesh%ncells))\n",
                                 "      allocate(energy%cp(mesh%ncells))\n"
                                 "      allocate(energy%lambda(mesh%ncells))\n"
                                 "      allocate(energy%rho_thermo(mesh%ncells))\n", "allocate thermo arrays")
        text = insert_after_once(text, "      energy%T = zero\n", "      energy%T_old = zero\n", "initialize T_old")
        text = insert_after_once(text, "      energy%qrad = zero\n",
                                 "      energy%cp = zero\n"
                                 "      energy%lambda = zero\n"
                                 "      energy%rho_thermo = zero\n", "initialize thermo arrays")
        changed = True

    if "deallocate(energy%T_old)" not in text:
        text = insert_after_once(text, "      if (allocated(energy%T)) deallocate(energy%T)\n",
                                 "      if (allocated(energy%T_old)) deallocate(energy%T_old)\n", "finalize T_old")
        text = insert_after_once(text, "      if (allocated(energy%qrad)) deallocate(energy%qrad)\n",
                                 "      if (allocated(energy%cp)) deallocate(energy%cp)\n"
                                 "      if (allocated(energy%lambda)) deallocate(energy%lambda)\n"
                                 "      if (allocated(energy%rho_thermo)) deallocate(energy%rho_thermo)\n", "finalize thermo arrays")
        changed = True

    if "call update_thermo_from_temperature_cantera(mesh, params, energy)" not in text:
        old_init = (
            "      energy%T = params%initial_T\n"
            "      call update_enthalpy_from_temperature_constant_cp(params, energy)\n"
            "      energy%h_old = energy%h\n"
            "      call zero_radiation_source(energy)\n"
        )
        new_init = (
            "      energy%T = params%initial_T\n"
            "      energy%T_old = energy%T\n"
            "\n"
            "      if (params%enable_cantera_thermo) then\n"
            "         call update_thermo_from_temperature_cantera(mesh, params, energy)\n"
            "      else\n"
            "         call update_enthalpy_from_temperature_constant_cp(params, energy)\n"
            "         energy%cp = params%energy_cp\n"
            "         energy%lambda = params%energy_lambda\n"
            "         energy%rho_thermo = params%rho\n"
            "      end if\n"
            "\n"
            "      energy%h_old = energy%h\n"
            "      call zero_radiation_source(energy)\n"
        )
        text = replace_once(text, old_init, new_init, "initialize_energy Cantera thermo path")
        changed = True

    if "if (allocated(energy%cp)) energy%cp = params%energy_cp" not in text:
        text = insert_after_once(text,
            "      energy%h = params%energy_reference_h + &\n"
            "                 params%energy_cp * (energy%T - params%energy_reference_T)\n",
            "\n"
            "      if (allocated(energy%cp)) energy%cp = params%energy_cp\n"
            "      if (allocated(energy%lambda)) energy%lambda = params%energy_lambda\n"
            "      if (allocated(energy%rho_thermo)) energy%rho_thermo = params%rho\n",
            "constant cp thermo array update")
        changed = True

    if "subroutine update_thermo_from_temperature_cantera" not in text:
        text = insert_after_once(text,
            "   end function enthalpy_from_temperature_value\n",
            THERMO_HELPERS,
            "Cantera thermo helper routines")
        changed = True

    guard = (
        "      if (params%enable_cantera_thermo) then\n"
        "         call fatal_error('energy', 'Patch 003 only supports constant-cp energy transport; disable enable_cantera_thermo')\n"
        "      end if\n"
        "\n"
    )
    if guard in text:
        text = text.replace(guard, "", 1)
        changed = True

    if "Patch 005c Stage-2A supports Cantera thermo only for no-species smoke tests" not in text:
        old_start = (
            "      ! Ensure off-rank neighbor values are current before taking the old state.\n"
            "      call flow_exchange_cell_scalar(flow, energy%h)\n"
            "      call flow_exchange_cell_scalar(flow, energy%T)\n"
            "      energy%h_old = energy%h\n"
        )
        new_start = (
            "      if (params%enable_cantera_thermo) then\n"
            "         if (params%enable_species) then\n"
            "            call fatal_error('energy', 'Patch 005c Stage-2A supports Cantera thermo only for no-species smoke tests')\n"
            "         end if\n"
            "         call update_thermo_from_temperature_cantera(mesh, params, energy)\n"
            "      end if\n"
            "\n"
            "      ! Ensure off-rank neighbor values are current before taking the old state.\n"
            "      call flow_exchange_cell_scalar(flow, energy%h)\n"
            "      call flow_exchange_cell_scalar(flow, energy%T)\n"
            "      if (allocated(energy%lambda)) call flow_exchange_cell_scalar(flow, energy%lambda)\n"
            "      energy%h_old = energy%h\n"
            "      if (allocated(energy%T_old)) energy%T_old = energy%T\n"
        )
        text = replace_once(text, old_start, new_start, "advance energy thermo refresh and T_old")
        changed = True

    if "h_other = boundary_enthalpy_from_temperature(params, T_other)" not in text:
        text = text.replace("                  h_other = enthalpy_from_temperature_value(params, T_other)\n",
                            "                  h_other = boundary_enthalpy_from_temperature(params, T_other)\n", 1)
        changed = True

    if "lambda_face" not in text:
        text = replace_once(text, "      real(rk) :: rhs_h, diff_term\n",
                            "      real(rk) :: rhs_h, diff_term, lambda_face\n", "lambda_face declaration")
        pattern = (
            r"            ! Fourier heat conduction uses grad\(T\), not grad\(h\)\.\n"
            r"            if \(do_diffusion \.and\. params%energy_lambda > zero\) then\n"
            r"               dist = energy_face_normal_distance\(mesh, bc, fid, c, neighbor\)\n"
            r"               diff_term = \(params%energy_lambda / params%rho\) \* &\n"
            r"                           \(T_other - T_cell\) / dist \* face_area\n"
            r"               rhs_h = rhs_h \+ diff_term\n"
            r"            end if\n"
        )
        repl = (
            "            ! Fourier heat conduction uses grad(T), not grad(h).\n"
            "            if (do_diffusion .and. (params%enable_cantera_thermo .or. params%energy_lambda > zero)) then\n"
            "               dist = energy_face_normal_distance(mesh, bc, fid, c, neighbor)\n"
            "\n"
            "               if (allocated(energy%lambda)) then\n"
            "                  if (neighbor > 0) then\n"
            "                     lambda_face = 0.5_rk * (energy%lambda(c) + energy%lambda(neighbor))\n"
            "                  else\n"
            "                     lambda_face = energy%lambda(c)\n"
            "                  end if\n"
            "               else\n"
            "                  lambda_face = params%energy_lambda\n"
            "               end if\n"
            "\n"
            "               if (lambda_face > zero) then\n"
            "                  diff_term = (lambda_face / params%rho) * &\n"
            "                              (T_other - T_cell) / dist * face_area\n"
            "                  rhs_h = rhs_h + diff_term\n"
            "               end if\n"
            "            end if\n"
        )
        text = replace_regex_once(text, pattern, repl, "Cantera lambda diffusion")
        changed = True

    if "call recover_temperature_cantera(mesh, params, energy)" not in text:
        old_recover = (
            "      call recover_temperature_constant_cp(params, energy)\n"
            "\n"
            "      ! Synchronize updated owned cells for output and the next step.\n"
            "      call flow_exchange_cell_scalar(flow, energy%h)\n"
            "      call flow_exchange_cell_scalar(flow, energy%T)\n"
        )
        new_recover = (
            "      if (params%enable_cantera_thermo) then\n"
            "         call recover_temperature_cantera(mesh, params, energy)\n"
            "         call update_thermo_from_temperature_cantera(mesh, params, energy)\n"
            "      else\n"
            "         call recover_temperature_constant_cp(params, energy)\n"
            "      end if\n"
            "\n"
            "      ! Synchronize updated owned cells for output and the next step.\n"
            "      call flow_exchange_cell_scalar(flow, energy%h)\n"
            "      call flow_exchange_cell_scalar(flow, energy%T)\n"
            "      if (allocated(energy%lambda)) call flow_exchange_cell_scalar(flow, energy%lambda)\n"
        )
        text = replace_once(text, old_recover, new_recover, "Cantera temperature recovery")
        changed = True

    if "energy%T_old(c)" not in text:
        old_diag = (
            "         old_T = params%energy_reference_T + &\n"
            "                 (energy%h_old(c) - params%energy_reference_h) / params%energy_cp\n"
        )
        new_diag = (
            "         if (allocated(energy%T_old)) then\n"
            "            old_T = energy%T_old(c)\n"
            "         else\n"
            "            old_T = params%energy_reference_T + &\n"
            "                    (energy%h_old(c) - params%energy_reference_h) / params%energy_cp\n"
            "         end if\n"
        )
        text = replace_once(text, old_diag, new_diag, "energy diagnostics T_old")
        changed = True

    if not changed:
        return PatchResult(path, False, "already patched")
    return write_if_changed(path, text, repo, backup_root, apply)


def sanity_check_repo(repo: Path) -> None:
    required = [
        repo / "src" / "cantera_interface.cpp",
        repo / "src" / "mod_transport_properties.f90",
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
    results = [
        patch_cantera_interface(repo, backup_root, apply),
        patch_mod_transport_properties(repo, backup_root, apply),
        patch_mod_energy(repo, backup_root, apply),
    ]

    print("\nPatch 005c summary")
    print("==================")
    print(f"repo:   {repo}")
    print(f"mode:   {'apply' if apply else 'dry-run'}")
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
        print("\nApplied patch. Next steps:")
        print("  1. inspect git diff")
        print("  2. set enable_cantera_thermo = .true. in &energy_input")
        print("  3. keep enable_species = .false. and enable_reactions = .false. for this Stage-2A smoke test")
        print("  4. rebuild and run rectangle_2D")
    else:
        print("\nNo changes were needed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
