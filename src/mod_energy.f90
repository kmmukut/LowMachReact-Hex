!> Energy field storage, diagnostics, and sensible-enthalpy transport.
!!
!! This module owns the enthalpy/temperature/radiation-source fields used by
!! the current constant-density energy path. The transported thermodynamic
!! state is `h`; temperature is recovered from `h`, composition `Y`, and the
!! background thermodynamic pressure `params%background_press`.
!!
!! In Cantera mode the stored sensible enthalpy is
!! `h_abs(T,Y,p0) - h_abs(T_ref,Y,p0)`, where `p0` is the background pressure
!! and `T_ref` is `params%energy_reference_T`. After species transport changes
!! `Y`, Option A is enforced: preserve transported `h` and recover
!! `T(h,Y,p0)`. Do not rebuild `h` from the old temperature and new
!! composition.
!!
!! Density ownership depends on mode.  In constant-density mode
!! `transport%rho = params%rho` and `energy%rho_thermo` is diagnostic.  In
!! guarded variable-density mode with `density_eos="cantera"`,
!! `transport%rho <- energy%rho_thermo` after Cantera thermo sync, so density
!! comes from the selected Cantera phase and its YAML EOS.  Heat conduction uses
!! `grad(T)`, not `grad(h)`. `qrad` is a volumetric source hook with `qrad > 0`
!! adding energy to the gas; external radiation coupling is future work.
!!
!! Reactions and reaction heat release are not implemented in the supported
!! path. The optional species-enthalpy diffusion correction
!! `-div(sum_k h_k J_k)` is controlled by
!! `enable_species_enthalpy_diffusion`.
module mod_energy
   use mpi_f08
   use iso_c_binding, only : c_char, c_double, c_int
   use mod_kinds, only : rk, zero, tiny_safe, fatal_error, output_unit
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_mpi_flow, only : flow_mpi_t, flow_exchange_cell_scalar
   use mod_transport_properties, only : transport_properties_t
   use mod_bc, only : bc_set_t, bc_periodic, face_effective_neighbor, &
                      patch_type_for_face, boundary_temperature, boundary_species
   use mod_fields, only : flow_fields_t
   use mod_profiler, only : profiler_start, profiler_stop
   implicit none

   private

   !> Cell-centered energy variables.
   type, public :: energy_fields_t
      real(rk), allocatable :: T(:)          !< Temperature [K].
      real(rk), allocatable :: T_old(:)      !< Previous temperature state [K].
      real(rk), allocatable :: h(:)          !< Transported mixture sensible enthalpy [J/kg].
      real(rk), allocatable :: h_old(:)      !< Previous enthalpy state [J/kg].
      real(rk), allocatable :: qrad(:)       !< Volumetric source term [W/m^3]; positive values heat the gas.
      real(rk), allocatable :: cp(:)         !< Mixture heat capacity at constant pressure [J/kg/K].
      real(rk), allocatable :: lambda(:)     !< Thermal conductivity [W/m/K].
      real(rk), allocatable :: rho_thermo(:) !< Cantera thermodynamic density [kg/m^3]; diagnostic unless variable-density mode syncs it into `transport%rho`.
      real(rk), allocatable :: operator_consistent_rho_h(:) !< Cellwise update-consistent rho*h [J/m^3].
      integer :: operator_consistent_rho_h_available = 0 !< 1 when cellwise update-consistent rho*h is valid.
      real(rk), allocatable :: species_enthalpy_diffusion(:) !< Specific enthalpy source from species diffusion [J/kg/s].
            real(rk) :: last_conductive_boundary_flux_out = zero !< Last local outward conductive boundary flux [W].
      integer :: last_conductive_boundary_flux_available = 0 !< 1 after conductive-boundary diagnostic is available.
      real(rk) :: cumulative_boundary_rho_h_advective_flux_out = zero !< Time-integrated outward advective rho*h boundary flux [J].
      real(rk) :: cumulative_boundary_rho_h_conductive_flux_out = zero !< Time-integrated outward conductive boundary flux [J].
      real(rk) :: cumulative_qrad_integral = zero !< Time-integrated qrad volume source [J].
      real(rk) :: cumulative_rho_species_hdiff_integral = zero !< Time-integrated rho*species enthalpy diffusion source [J].
      integer :: cumulative_energy_budget_available = 0 !< 1 after cumulative energy budget terms are available.
      real(rk) :: last_energy_update_delta_rate_integral = zero !< Last global d/dt integral(rho*h) from update [W].
      real(rk) :: last_energy_update_rhs_integral = zero !< Last global energy RHS integral from update [W].
      real(rk) :: last_energy_update_balance_defect = zero !< Last direct energy update closure defect [W].
      real(rk) :: relative_last_energy_update_balance_defect = zero !< Relative direct update closure defect.
      real(rk) :: last_operator_consistent_rho_h_integral = zero !< Last global rho*h integral using update density/time level.
      real(rk) :: cumulative_energy_update_delta_integral = zero !< Time-integrated conservative energy-update delta [J].
      real(rk) :: cumulative_energy_update_rhs_integral = zero !< Time-integrated conservative energy-update RHS [J].
logical :: initialized = .false.  !< True after allocation/initialization.
   end type energy_fields_t

   public :: allocate_energy, initialize_energy, finalize_energy
   public :: update_enthalpy_from_temperature_constant_cp
   public :: recover_temperature_constant_cp
   public :: recover_temperature_cantera
   public :: zero_radiation_source
   public :: advance_energy_transport
   public :: write_energy_diagnostics_header, write_energy_diagnostics_row
   public :: write_cantera_cache_stats

   !> C-Binding interface for Cantera thermodynamics.
   !!
   !! The C++ bridge expects temperatures \(T\), thermodynamic pressures \(p_0\),
   !! and normalized solver species mass fractions \(Y_k\).  It returns sensible
   !! enthalpy using
   !! $$
   !! h_{sens}(T,Y,p_0) =
   !! h_{abs}(T,Y,p_0) - h_{abs}(T_{ref},Y,p_0)
   !! $$
   !! and recovers temperature by adding the same-composition reference
   !! enthalpy before Cantera HP inversion.  Projection pressure is never passed
   !! to these calls.
   interface

      subroutine cantera_get_cache_stats_c(stats_out, nstats) bind(c, name="cantera_get_cache_stats_c")
         import :: c_double, c_int
         real(c_double), intent(out) :: stats_out(*)
         integer(c_int), value :: nstats
      end subroutine cantera_get_cache_stats_c

      subroutine cantera_update_thermo_c(ncells, T, P, nspecies, Y_in, h_out, cp_out, lambda_out, &
                                         rho_thermo_out, T_ref, species_names_flat, name_len) bind(c, name="cantera_update_thermo_c")
         import :: c_char, c_double, c_int
         integer(c_int), value :: ncells
         real(c_double), intent(in) :: T(ncells), P(ncells)
         integer(c_int), value :: nspecies
         real(c_double), intent(in) :: Y_in(*)
         real(c_double), intent(out) :: h_out(ncells), cp_out(ncells), lambda_out(ncells), rho_thermo_out(ncells)
         real(c_double), value :: T_ref
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_update_thermo_c

      subroutine cantera_recover_temperature_from_h_c(ncells, h_in, P, nspecies, Y_in, T_out, &
                                                      T_ref, species_names_flat, name_len) bind(c, name="cantera_recover_temperature_from_h_c")
         import :: c_char, c_double, c_int
         integer(c_int), value :: ncells
         real(c_double), intent(in) :: h_in(ncells), P(ncells)
         integer(c_int), value :: nspecies
         real(c_double), intent(in) :: Y_in(*)
         real(c_double), intent(out) :: T_out(ncells)
         real(c_double), value :: T_ref
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_recover_temperature_from_h_c


      subroutine cantera_species_sensible_enthalpies_c(npoints, T, P, nspecies, Y_in, hk_out, &
                                                           T_ref, species_names_flat, name_len) &
         bind(c, name="cantera_species_sensible_enthalpies_c")
         import :: c_char, c_double, c_int
         integer(c_int), value :: npoints
         real(c_double), intent(in) :: T(npoints), P(npoints)
         integer(c_int), value :: nspecies
         real(c_double), intent(in) :: Y_in(*)
         real(c_double), intent(out) :: hk_out(*)
         real(c_double), value :: T_ref
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_species_sensible_enthalpies_c

      subroutine cantera_recover_temperature_and_update_thermo_c(ncells, h_in, P, nspecies, Y_in, T_out, &
                                                                 cp_out, lambda_out, rho_thermo_out, &
                                                                 T_ref, species_names_flat, name_len) &
         bind(c, name="cantera_recover_temperature_and_update_thermo_c")
         import :: c_char, c_double, c_int
         integer(c_int), value :: ncells
         real(c_double), intent(in) :: h_in(ncells), P(ncells)
         integer(c_int), value :: nspecies
         real(c_double), intent(in) :: Y_in(*)
         real(c_double), intent(out) :: T_out(ncells), cp_out(ncells), lambda_out(ncells), rho_thermo_out(ncells)
         real(c_double), value :: T_ref
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_recover_temperature_and_update_thermo_c

   end interface

contains

   !> Allocate all energy arrays for the mesh.
   subroutine allocate_energy(mesh, energy)
      type(mesh_t), intent(in) :: mesh
      type(energy_fields_t), intent(inout) :: energy

      call finalize_energy(energy)

      allocate(energy%T(mesh%ncells))
      allocate(energy%T_old(mesh%ncells))
      allocate(energy%h(mesh%ncells))
      allocate(energy%h_old(mesh%ncells))
      allocate(energy%qrad(mesh%ncells))
      allocate(energy%cp(mesh%ncells))
      allocate(energy%lambda(mesh%ncells))
      allocate(energy%rho_thermo(mesh%ncells))
      allocate(energy%operator_consistent_rho_h(mesh%ncells))
      allocate(energy%species_enthalpy_diffusion(mesh%ncells))

      energy%T = zero
      energy%T_old = zero
      energy%h = zero
      energy%h_old = zero
      energy%qrad = zero
      energy%cp = zero
      energy%lambda = zero
      energy%rho_thermo = zero
      energy%operator_consistent_rho_h = zero
      energy%operator_consistent_rho_h_available = 0
      energy%species_enthalpy_diffusion = zero
      energy%last_conductive_boundary_flux_out = zero
      energy%last_conductive_boundary_flux_available = 0
      energy%cumulative_boundary_rho_h_advective_flux_out = zero
      energy%cumulative_boundary_rho_h_conductive_flux_out = zero
      energy%cumulative_qrad_integral = zero
      energy%cumulative_rho_species_hdiff_integral = zero
      energy%cumulative_energy_budget_available = 0
      energy%last_energy_update_delta_rate_integral = zero
      energy%last_energy_update_rhs_integral = zero
      energy%last_energy_update_balance_defect = zero
      energy%relative_last_energy_update_balance_defect = zero
      energy%last_operator_consistent_rho_h_integral = zero
      energy%cumulative_energy_update_delta_integral = zero
      energy%cumulative_energy_update_rhs_integral = zero
      energy%initialized = .true.
   end subroutine allocate_energy


   !> Initialize energy fields from case parameters.
   subroutine initialize_energy(mesh, params, energy, species_Y)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy
      real(rk), intent(in), optional :: species_Y(:,:)

      call allocate_energy(mesh, energy)

      energy%T = params%initial_T
      energy%T_old = energy%T

      if (params%enable_cantera_thermo) then
         call update_thermo_from_temperature_cantera(mesh, params, energy, species_Y)
      else
         call update_enthalpy_from_temperature_constant_cp(params, energy)
         energy%cp = params%energy_cp
         energy%lambda = params%energy_lambda
         energy%rho_thermo = params%rho
      end if

      energy%h_old = energy%h
      call zero_radiation_source(energy)
   end subroutine initialize_energy


   !> Deallocate all energy arrays.
   subroutine finalize_energy(energy)
      type(energy_fields_t), intent(inout) :: energy

      if (allocated(energy%T)) deallocate(energy%T)
      if (allocated(energy%T_old)) deallocate(energy%T_old)
      if (allocated(energy%h)) deallocate(energy%h)
      if (allocated(energy%h_old)) deallocate(energy%h_old)
      if (allocated(energy%qrad)) deallocate(energy%qrad)
      if (allocated(energy%cp)) deallocate(energy%cp)
      if (allocated(energy%lambda)) deallocate(energy%lambda)
      if (allocated(energy%rho_thermo)) deallocate(energy%rho_thermo)
      if (allocated(energy%operator_consistent_rho_h)) deallocate(energy%operator_consistent_rho_h)
      if (allocated(energy%species_enthalpy_diffusion)) deallocate(energy%species_enthalpy_diffusion)

      energy%operator_consistent_rho_h_available = 0
      energy%initialized = .false.
   end subroutine finalize_energy


   !> Update h from T using the constant-cp thermodynamic model.
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

      if (allocated(energy%cp)) energy%cp = params%energy_cp
      if (allocated(energy%lambda)) energy%lambda = params%energy_lambda
      if (allocated(energy%rho_thermo)) energy%rho_thermo = params%rho
   end subroutine update_enthalpy_from_temperature_constant_cp


   !> Recover T from h using the constant-cp thermodynamic model.
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
            c_names_flat((k-1)*n_len + j) = char(iachar(params%species_name(k)(j:j)), kind=c_char)
         end do
      end do
   end subroutine build_c_species_names


   !> Build a cellwise thermodynamic composition array for Cantera calls.
   !!
   !! If transported species are available, they are used directly and clipped
   !! to non-negative values. If their sum is positive, the local composition is
   !! normalized before passing to Cantera. If no transported species are
   !! available, use the single-species/default mixture prepared during Cantera
   !! initialization, usually thermo_default_species.
   subroutine build_thermo_Y(params, ncells, Y_local, species_Y)
      type(case_params_t), intent(in) :: params
      integer, intent(in) :: ncells
      real(rk), allocatable, intent(out) :: Y_local(:,:)
      real(rk), intent(in), optional :: species_Y(:,:)

      integer :: c, k, nsp
      real(rk) :: sum_Y

      nsp = max(1, params%nspecies)
      allocate(Y_local(nsp, ncells))
      Y_local = zero

      if (present(species_Y) .and. params%enable_species .and. params%nspecies > 0) then
         if (size(species_Y, 1) < params%nspecies .or. size(species_Y, 2) < ncells) then
            call fatal_error('energy', 'species_Y has incompatible shape for Cantera thermo update')
         end if

         do c = 1, ncells
            sum_Y = zero
            do k = 1, params%nspecies
               Y_local(k, c) = max(zero, species_Y(k, c))
               sum_Y = sum_Y + Y_local(k, c)
            end do

            if (sum_Y > tiny_safe) then
               Y_local(1:params%nspecies, c) = Y_local(1:params%nspecies, c) / sum_Y
            else
               Y_local(1, c) = 1.0_rk
            end if
         end do
      else
         Y_local(1, :) = 1.0_rk
      end if
   end subroutine build_thermo_Y


   !> Build a boundary thermodynamic composition vector for fixed-T boundaries.
   !!
   !! Fixed-composition species boundaries use `patch_Y` through
   !! `boundary_species`; non-Dirichlet species boundaries use the adjacent
   !! interior composition.  Fixed-temperature enthalpy is then evaluated as
   !! \(h_b = h(T_b,Y_b,p_0)\).
   subroutine build_boundary_thermo_Y(mesh, bc, params, face_id, cell_id, species_Y, Y_point)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      real(rk), intent(in) :: species_Y(:,:)
      real(rk), intent(out) :: Y_point(:)

      integer :: k
      real(rk) :: ext_Y, sum_Y
      logical :: is_dirichlet

      if (size(Y_point) < params%nspecies) then
         call fatal_error('energy', 'Y_point too small in build_boundary_thermo_Y')
      end if
      if (size(species_Y, 1) < params%nspecies .or. size(species_Y, 2) < cell_id) then
         call fatal_error('energy', 'species_Y has incompatible shape in boundary thermo composition')
      end if

      sum_Y = zero
      do k = 1, params%nspecies
         call boundary_species(mesh, bc, face_id, k, species_Y(k, cell_id), ext_Y, is_dirichlet)
         Y_point(k) = max(zero, ext_Y)
         sum_Y = sum_Y + Y_point(k)
      end do

      if (sum_Y > tiny_safe) then
         Y_point(1:params%nspecies) = Y_point(1:params%nspecies) / sum_Y
      else
         Y_point = zero
         Y_point(1) = 1.0_rk
      end if
   end subroutine build_boundary_thermo_Y



   !> Update h, cp, lambda, and diagnostic thermo density from current T.
   subroutine update_thermo_from_temperature_cantera(mesh, params, energy, species_Y)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy
      real(rk), intent(in), optional :: species_Y(:,:)

      integer :: n_len
      real(rk), allocatable :: P_arr(:), Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (.not. params%enable_cantera_thermo) return

      if (params%nspecies <= 0) then
         call fatal_error('energy', 'enable_cantera_thermo requires at least one inert thermo species')
      end if

      allocate(P_arr(mesh%ncells))
      P_arr = params%background_press
      call build_thermo_Y(params, mesh%ncells, Y_local, species_Y)
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_update_thermo_c(mesh%ncells, energy%T, P_arr, params%nspecies, Y_local, &
                                   energy%h, energy%cp, energy%lambda, energy%rho_thermo, &
                                   params%energy_reference_T, c_names_flat, n_len)

      deallocate(P_arr)
      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine update_thermo_from_temperature_cantera


   !> Recover T from h using Cantera HPY inversion.
   subroutine recover_temperature_cantera(mesh, params, energy, species_Y)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy
      real(rk), intent(in), optional :: species_Y(:,:)

      integer :: n_len
      real(rk), allocatable :: P_arr(:), Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (.not. params%enable_cantera_thermo) return

      if (params%nspecies <= 0) then
         call fatal_error('energy', 'enable_cantera_thermo requires at least one inert thermo species')
      end if

      allocate(P_arr(mesh%ncells))
      P_arr = params%background_press
      call build_thermo_Y(params, mesh%ncells, Y_local, species_Y)
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_recover_temperature_from_h_c(mesh%ncells, &
                                   energy%h, &
                                   P_arr, &
                                   params%nspecies, &
                                   Y_local, &
                                   energy%T, &
                                   params%energy_reference_T, &
                                   c_names_flat, &
                                   n_len)

      deallocate(P_arr)
      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine recover_temperature_cantera


   !> Recover T from h and refresh cp/lambda/rho_thermo in one Cantera sync.
   !!
   !! This is the preferred energy-step thermo sync routine:
   !! `(T, cp, lambda, rho_thermo) = sync(h, Y, p0)`. It preserves the
   !! transported sensible enthalpy field and updates only derived thermo
   !! state: T, cp, lambda, and rho_thermo.  The thermo-sync cache key in the
   !! C++ bridge is `h`, `p0`, and `Y`; `rho_thermo` becomes active density only
   !! after `sync_active_density_from_thermo` in variable-density mode.
   subroutine recover_temperature_and_update_thermo_cantera(mesh, params, energy, species_Y)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(energy_fields_t), intent(inout) :: energy
      real(rk), intent(in), optional :: species_Y(:,:)

      integer :: n_len
      real(rk), allocatable :: P_arr(:), Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (.not. params%enable_cantera_thermo) return

      if (params%nspecies <= 0) then
         call fatal_error('energy', 'enable_cantera_thermo requires at least one inert thermo species')
      end if

      allocate(P_arr(mesh%ncells))
      P_arr = params%background_press
      call build_thermo_Y(params, mesh%ncells, Y_local, species_Y)
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_recover_temperature_and_update_thermo_c(mesh%ncells, &
                                   energy%h, &
                                   P_arr, &
                                   params%nspecies, &
                                   Y_local, &
                                   energy%T, &
                                   energy%cp, &
                                   energy%lambda, &
                                   energy%rho_thermo, &
                                   params%energy_reference_T, &
                                   c_names_flat, &
                                   n_len)

      deallocate(P_arr)
      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine recover_temperature_and_update_thermo_cantera


   !> Compute species sensible enthalpies h_k(T)-h_k(T_ref) [J/kg_k] for all cells.
   subroutine compute_species_sensible_enthalpies_cantera(mesh, params, T_state, species_Y, hk_species)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: T_state(:)
      real(rk), intent(in) :: species_Y(:,:)
      real(rk), intent(out) :: hk_species(:,:)

      integer :: n_len
      real(rk), allocatable :: P_arr(:), Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (.not. params%enable_cantera_thermo) then
         call fatal_error('energy', 'species sensible enthalpies require Cantera thermo')
      end if
      if (params%nspecies <= 0) then
         call fatal_error('energy', 'species sensible enthalpies require at least one species')
      end if
      if (size(T_state) < mesh%ncells) then
         call fatal_error('energy', 'T_state has incompatible size for species sensible enthalpies')
      end if
      if (size(hk_species, 1) < params%nspecies .or. size(hk_species, 2) < mesh%ncells) then
         call fatal_error('energy', 'hk_species has incompatible shape')
      end if

      allocate(P_arr(mesh%ncells))
      P_arr = params%background_press
      call build_thermo_Y(params, mesh%ncells, Y_local, species_Y)
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_species_sensible_enthalpies_c(mesh%ncells, T_state, P_arr, params%nspecies, &
                                                 Y_local, hk_species, params%energy_reference_T, &
                                                 c_names_flat, n_len)

      deallocate(P_arr)
      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine compute_species_sensible_enthalpies_cantera


   !> Compute species sensible enthalpies for one boundary/face state.
   subroutine species_sensible_enthalpies_from_temperature_cantera_point(params, temperature, hk_value, Y_point)
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: temperature
      real(rk), intent(out) :: hk_value(:)
      real(rk), intent(in) :: Y_point(:)

      integer :: n_len, nsp
      real(rk) :: T_arr(1), P_arr(1)
      real(rk), allocatable :: Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (params%nspecies <= 0) then
         call fatal_error('energy', 'Cantera boundary species enthalpies require at least one species')
      end if
      if (size(hk_value) < params%nspecies) then
         call fatal_error('energy', 'species enthalpy point hk array has incompatible size')
      end if
      if (size(Y_point) < params%nspecies) then
         call fatal_error('energy', 'species enthalpy point Y array has incompatible size')
      end if

      nsp = max(1, params%nspecies)
      allocate(Y_local(nsp, 1))
      Y_local = zero
      Y_local(1:params%nspecies, 1) = max(zero, Y_point(1:params%nspecies))
      if (sum(Y_local(1:params%nspecies, 1)) > tiny_safe) then
         Y_local(1:params%nspecies, 1) = Y_local(1:params%nspecies, 1) / sum(Y_local(1:params%nspecies, 1))
      else
         Y_local(1, 1) = 1.0_rk
      end if

      T_arr(1) = temperature
      P_arr(1) = params%background_press
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_species_sensible_enthalpies_c(1, T_arr, P_arr, params%nspecies, Y_local, &
                                                 hk_value, params%energy_reference_T, &
                                                 c_names_flat, n_len)

      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine species_sensible_enthalpies_from_temperature_cantera_point


   !> Compute h(T,Y,p0) for one boundary state using Cantera.
   subroutine enthalpy_from_temperature_cantera_point(params, temperature, h_value, Y_point)
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: temperature
      real(rk), intent(out) :: h_value
      real(rk), intent(in), optional :: Y_point(:)

      integer :: n_len, nsp
      real(rk) :: T_arr(1), P_arr(1)
      real(rk) :: h_arr(1), cp_arr(1), lambda_arr(1), rho_arr(1)
      real(rk), allocatable :: Y_local(:,:)
      character(kind=c_char), allocatable :: c_names_flat(:)

      if (params%nspecies <= 0) then
         call fatal_error('energy', 'Cantera boundary enthalpy requires at least one species')
      end if

      nsp = max(1, params%nspecies)
      allocate(Y_local(nsp, 1))
      Y_local = zero

      if (present(Y_point)) then
         if (size(Y_point) < params%nspecies) then
            call fatal_error('energy', 'Y_point has incompatible size for Cantera boundary enthalpy')
         end if
         Y_local(1:params%nspecies, 1) = max(zero, Y_point(1:params%nspecies))
         if (sum(Y_local(1:params%nspecies, 1)) > tiny_safe) then
            Y_local(1:params%nspecies, 1) = Y_local(1:params%nspecies, 1) / &
                                            sum(Y_local(1:params%nspecies, 1))
         else
            Y_local(1, 1) = 1.0_rk
         end if
      else
         Y_local(1, 1) = 1.0_rk
      end if

      T_arr(1) = temperature
      P_arr(1) = params%background_press
      call build_c_species_names(params, c_names_flat, n_len)

      call cantera_update_thermo_c(1, T_arr, P_arr, params%nspecies, Y_local, &
                                   h_arr, cp_arr, lambda_arr, rho_arr, &
                                   params%energy_reference_T, c_names_flat, n_len)

      h_value = h_arr(1)

      deallocate(Y_local)
      deallocate(c_names_flat)
   end subroutine enthalpy_from_temperature_cantera_point



   !> Convert a boundary temperature to enthalpy using the active thermo model.
   
   !> Boundary enthalpy from a boundary temperature.
   !!
   !! When a transported boundary composition is available, the fixed-temperature
   !! boundary enthalpy is evaluated as h(T_b,Y_b,p0). This keeps hot/cold
   !! species inlets thermodynamically consistent.
   function boundary_enthalpy_from_temperature(mesh, params, temperature, Y_point) result(h_value)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: temperature
      real(rk), intent(in), optional :: Y_point(:)
      real(rk) :: h_value

      if (.not. params%enable_cantera_thermo) then
         h_value = enthalpy_from_temperature_value(params, temperature)
         return
      end if

      call enthalpy_from_temperature_cantera_point(params, temperature, h_value, Y_point)

      associate(dummy => mesh%ncells)
      end associate
   end function boundary_enthalpy_from_temperature



   !> Reset the volumetric energy source to zero.
   subroutine zero_radiation_source(energy)
      type(energy_fields_t), intent(inout) :: energy

      if (allocated(energy%qrad)) energy%qrad = zero
   end subroutine zero_radiation_source


   !> Advance transported sensible enthalpy one explicit finite-volume step.
   !!
   !! The explicit finite-volume update is:
   !!
   !!   V dh/dt = - sum_f F_f h_f
   !!            + (1/rho) sum_f lambda A_f (T_nb - T_c)/d_f
   !!            + (qrad/rho) V
   !!
   !! `fields%face_flux` is volumetric flux in \( \mathrm{m^3/s} \).
   !! `fields%mass_flux` is density-weighted flux in \( \mathrm{kg/s} \).
   !! Both are stored owner-to-neighbor and re-oriented outward from the
   !! currently updated cell before upwind advection.
   !!
   !! Constant-density mode updates `h` with `params%rho`.  Variable-density
   !! mode uses a conservative `rho*h` update with `fields%mass_flux`,
   !! `transport%rho_old`, and current `transport%rho`, then stores the
   !! transported specific enthalpy back in `energy%h`.
   subroutine advance_energy_transport(mesh, flow, bc, params, fields, energy, species_Y, transport)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(inout) :: energy
      real(rk), intent(in), optional :: species_Y(:,:)
      type(transport_properties_t), intent(in), optional :: transport

      integer :: c, lf, fid, neighbor, k, ierr
      real(rk) :: flux_out
      real(rk) :: face_area, dist, vol
      real(rk) :: h_cell, h_other, h_face
      real(rk) :: T_cell, T_other
      real(rk) :: rhs_h, diff_term, lambda_face
      real(rk) :: conductive_boundary_flux_out_local
      real(rk) :: advective_boundary_rho_h_flux_out_local
      real(rk) :: qrad_integral_step_local
      real(rk) :: rho_species_hdiff_integral_step_local
      real(rk) :: rho_diag
      real(rk) :: energy_update_delta_rate_local, energy_update_rhs_integral_local
      real(rk) :: energy_update_delta_rate_global, energy_update_rhs_integral_global
      real(rk) :: energy_update_balance_global, energy_update_denom
      real(rk) :: update_send(3), update_recv(3)
      real(rk) :: operator_consistent_rho_h_integral_local
      real(rk) :: operator_consistent_rho_h_integral_global
      real(rk) :: rho_cell, rho_old_cell, rho_face
      real(rk) :: sum_diff_flux, diff_face, species_hdiff_rhs
      real(rk) :: Y_cell, Y_other
      logical :: is_dirichlet, do_diffusion, use_species_hdiff, has_species_diffusion_face
      logical :: use_variable_density_energy
      real(rk), allocatable :: Y_point(:)
      ! Explicitly null-initialized pointer scratch avoids false gfortran
      ! release-build descriptor warnings on the conditional h_k diffusion path.
      real(rk), pointer :: hk_species(:,:) => null(), hk_boundary(:) => null()
      real(rk), allocatable :: sp_diff_flux(:), sp_Y_face_lin(:)

      if (.not. params%enable_energy) return

      if (.not. energy%initialized) then
         call fatal_error('energy', 'advance requested before energy initialization')
      end if

      if (params%rho <= tiny_safe) call fatal_error('energy', 'rho must be positive')
      if (params%energy_cp <= tiny_safe) call fatal_error('energy', 'energy_cp must be positive')
      if (params%energy_lambda < zero) call fatal_error('energy', 'energy_lambda must be non-negative')

      use_variable_density_energy = params%enable_variable_density
      if (use_variable_density_energy) then
         if (.not. present(transport)) then
            call fatal_error('energy', 'variable-density energy transport requires transport properties')
         end if
         if (.not. allocated(transport%rho)) then
            call fatal_error('energy', 'variable-density energy transport requires transport%rho')
         end if
         if (.not. allocated(transport%rho_old)) then
            call fatal_error('energy', 'variable-density energy transport requires transport%rho_old')
         end if
         if (.not. allocated(fields%mass_flux)) then
            call fatal_error('energy', 'variable-density energy transport requires fields%mass_flux')
         end if
      end if

      use_species_hdiff = params%enable_species_enthalpy_diffusion
      if (use_species_hdiff) then
         if (.not. present(species_Y)) then
            call fatal_error('energy', 'species enthalpy diffusion requires species_Y')
         end if
         if (.not. present(transport)) then
            call fatal_error('energy', 'species enthalpy diffusion requires transport properties')
         end if
         if (.not. allocated(transport%diffusivity)) then
            call fatal_error('energy', 'species enthalpy diffusion requires transport diffusivity')
         end if
         if (.not. params%enable_cantera_thermo) then
            call fatal_error('energy', 'species enthalpy diffusion requires Cantera thermo')
         end if
         if (params%nspecies <= 0) then
            call fatal_error('energy', 'species enthalpy diffusion requires nspecies > 0')
         end if
         if (size(species_Y, 1) < params%nspecies .or. size(species_Y, 2) < mesh%ncells) then
            call fatal_error('energy', 'species enthalpy diffusion species_Y has incompatible shape')
         end if
         if (size(transport%diffusivity, 1) < params%nspecies .or. &
             size(transport%diffusivity, 2) < mesh%ncells) then
            call fatal_error('energy', 'species enthalpy diffusion diffusivity has incompatible shape')
         end if
      end if

      ! Option A: preserve transported enthalpy across species updates.
      ! Species transport may have changed Y before this energy step.  Do not
      ! recompute h from the old temperature and the new composition.  Instead,
      ! keep h as the transported state, recover T(h,Y,p0), then refresh the
      ! thermo/transport properties used by the enthalpy equation.
      call profiler_start('Energy_Exchange_H')
      call flow_exchange_cell_scalar(flow, energy%h)
      call profiler_stop('Energy_Exchange_H')

      energy%h_old = energy%h

      if (params%enable_cantera_thermo) then
         ! A pre-flux thermo sync is required only when species transport may
         ! have changed the composition since the previous post-flux sync.
         ! Without species, current T/cp/lambda/rho_thermo are already valid
         ! from initialization or the previous energy step's post-sync.
         if (params%enable_species .and. present(species_Y) .and. params%nspecies > 0) then
            call write_variable_density_energy_debug(mesh, flow, params, fields, energy, &
                                                      'pre_cantera_presync', species_Y, transport)
            call profiler_start('Energy_Cantera_PreSync')
            call recover_temperature_and_update_thermo_cantera(mesh, params, energy, species_Y)
            call profiler_stop('Energy_Cantera_PreSync')

            energy%h = energy%h_old
         end if
      end if

      ! Ensure off-rank temperature/property values are current before fluxes.
      call profiler_start('Energy_PreFlux_Exchange')
      call flow_exchange_cell_scalar(flow, energy%T)
      if (allocated(energy%lambda)) call flow_exchange_cell_scalar(flow, energy%lambda)
      call profiler_stop('Energy_PreFlux_Exchange')

      if (params%enable_cantera_thermo .and. params%enable_species .and. &
          present(species_Y) .and. params%nspecies > 0) then
         allocate(Y_point(params%nspecies))
      end if
      if (allocated(energy%T_old)) energy%T_old = energy%T

      if (use_species_hdiff) then
         if (.not. allocated(Y_point)) allocate(Y_point(params%nspecies))
         allocate(hk_species(params%nspecies, mesh%ncells))
         allocate(hk_boundary(params%nspecies))
         allocate(sp_diff_flux(params%nspecies))
         allocate(sp_Y_face_lin(params%nspecies))
         hk_species = zero
         hk_boundary = zero
         sp_diff_flux = zero
         sp_Y_face_lin = zero

         call profiler_start('Energy_Cantera_SpeciesH')
         call compute_species_sensible_enthalpies_cantera(mesh, params, energy%T, species_Y, hk_species)
         call profiler_stop('Energy_Cantera_SpeciesH')
      end if

      if (allocated(energy%species_enthalpy_diffusion)) energy%species_enthalpy_diffusion = zero
      conductive_boundary_flux_out_local = zero
      energy_update_delta_rate_local = zero
      energy_update_rhs_integral_local = zero
      operator_consistent_rho_h_integral_local = zero
      if (allocated(energy%operator_consistent_rho_h)) energy%operator_consistent_rho_h = zero
      energy%operator_consistent_rho_h_available = 0
      advective_boundary_rho_h_flux_out_local = zero
      qrad_integral_step_local = zero
      rho_species_hdiff_integral_step_local = zero

      call profiler_start('Energy_Flux_Update')
      do c = flow%first_cell, flow%last_cell
         rhs_h = zero
         h_cell = energy%h_old(c)
         T_cell = energy%T(c)

         do lf = 1, mesh%ncell_faces(c)
            fid = mesh%cell_faces(lf, c)

            ! `fields%face_flux` and `fields%mass_flux` are oriented
            ! owner -> neighbor. Reorient outward from the current cell.
            ! Variable-density energy uses mass flux for conservative rho*h
            ! advection; constant-density energy uses volumetric flux.
            if (mesh%faces(fid)%owner == c) then
               if (use_variable_density_energy) then
                  flux_out = fields%mass_flux(fid)
               else
                  flux_out = fields%face_flux(fid)
               end if
            else
               if (use_variable_density_energy) then
                  flux_out = -fields%mass_flux(fid)
               else
                  flux_out = -fields%face_flux(fid)
               end if
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
                  if (allocated(Y_point)) then
                     call build_boundary_thermo_Y(mesh, bc, params, fid, c, species_Y, Y_point)
                     h_other = boundary_enthalpy_from_temperature(mesh, params, T_other, Y_point)
                  else
                     h_other = boundary_enthalpy_from_temperature(mesh, params, T_other)
                  end if
                  do_diffusion = .true.
               else
                  h_other = h_cell
                  T_other = T_cell
                  do_diffusion = .false.
               end if
            end if

            ! Upwind advection of h using outward volumetric flux in the
            ! constant-density path and outward mass flux in variable-density
            ! mode.
            if (flux_out >= zero) then
               h_face = h_cell
            else
               h_face = h_other
            end if
            rhs_h = rhs_h - flux_out * h_face

            if (neighbor <= 0) then
               if (use_variable_density_energy) then
                  advective_boundary_rho_h_flux_out_local = advective_boundary_rho_h_flux_out_local + flux_out * h_face
               else
                  advective_boundary_rho_h_flux_out_local = advective_boundary_rho_h_flux_out_local + params%rho * flux_out * h_face
               end if
            end if

            ! Fourier heat conduction uses grad(T), not grad(h).
            if (do_diffusion .and. (params%enable_cantera_thermo .or. params%energy_lambda > zero)) then
               dist = energy_face_normal_distance(mesh, bc, fid, c, neighbor)

               if (allocated(energy%lambda)) then
                  if (neighbor > 0) then
                     lambda_face = 0.5_rk * (energy%lambda(c) + energy%lambda(neighbor))
                  else
                     lambda_face = energy%lambda(c)
                  end if
               else
                  lambda_face = params%energy_lambda
               end if

               if (lambda_face > zero) then
                  if (use_variable_density_energy) then
                     diff_term = lambda_face * (T_other - T_cell) / dist * face_area
                  else
                     diff_term = (lambda_face / params%rho) * &
                                 (T_other - T_cell) / dist * face_area
                  end if
                  rhs_h = rhs_h + diff_term
                  if (neighbor <= 0) then
                     conductive_boundary_flux_out_local = conductive_boundary_flux_out_local - diff_term
                  end if
               end if
            end if


            ! Species-enthalpy diffusion correction:
            ! corrected_diff_flux(k) matches mod_species and equals -J_k,out*A/rho.
            ! Therefore sum_k h_k,face * corrected_diff_flux(k) is added to the
            ! sensible-enthalpy RHS for this cell.
            if (use_species_hdiff) then
               dist = energy_face_normal_distance(mesh, bc, fid, c, neighbor)
               sum_diff_flux = zero
               has_species_diffusion_face = .false.

               do k = 1, params%nspecies
                  Y_cell = species_Y(k, c)

                  if (neighbor > 0) then
                     Y_other = species_Y(k, neighbor)
                     is_dirichlet = .true.
                  else
                     call boundary_species(mesh, bc, fid, k, Y_cell, Y_other, is_dirichlet)
                  end if

                  sp_diff_flux(k) = zero
                  sp_Y_face_lin(k) = 0.5_rk * (Y_cell + Y_other)

                  if (neighbor /= 0 .or. is_dirichlet) then
                     has_species_diffusion_face = .true.
                     if (neighbor == 0) then
                        diff_face = transport%diffusivity(k, c)
                     else
                        diff_face = 0.5_rk * (transport%diffusivity(k, c) + transport%diffusivity(k, neighbor))
                     end if
                     if (use_variable_density_energy) then
                        if (neighbor == 0) then
                           rho_face = transport%rho(c)
                        else
                           rho_face = 0.5_rk * (transport%rho(c) + transport%rho(neighbor))
                        end if
                        sp_diff_flux(k) = rho_face * diff_face * (Y_other - Y_cell) / dist * face_area
                     else
                        sp_diff_flux(k) = diff_face * (Y_other - Y_cell) / dist * face_area
                     end if
                  end if

                  sum_diff_flux = sum_diff_flux + sp_diff_flux(k)
               end do

               species_hdiff_rhs = zero
               if (has_species_diffusion_face) then
                  if (neighbor > 0) then
                     do k = 1, params%nspecies
                        species_hdiff_rhs = species_hdiff_rhs + &
                           0.5_rk * (hk_species(k, c) + hk_species(k, neighbor)) * &
                           (sp_diff_flux(k) - sp_Y_face_lin(k) * sum_diff_flux)
                     end do
                  else
                     call build_boundary_thermo_Y(mesh, bc, params, fid, c, species_Y, Y_point)
                     call species_sensible_enthalpies_from_temperature_cantera_point(params, T_other, hk_boundary, Y_point)

                     do k = 1, params%nspecies
                        species_hdiff_rhs = species_hdiff_rhs + &
                           0.5_rk * (hk_species(k, c) + hk_boundary(k)) * &
                           (sp_diff_flux(k) - sp_Y_face_lin(k) * sum_diff_flux)
                     end do
                  end if
               end if

               rhs_h = rhs_h + species_hdiff_rhs
               if (allocated(energy%species_enthalpy_diffusion)) then
                  if (use_variable_density_energy) then
                     rho_cell = max(transport%rho(c), tiny_safe)
                     energy%species_enthalpy_diffusion(c) = energy%species_enthalpy_diffusion(c) + &
                                                            species_hdiff_rhs / &
                                                            (rho_cell * mesh%cells(c)%volume)
                  else
                     energy%species_enthalpy_diffusion(c) = energy%species_enthalpy_diffusion(c) + &
                                                            species_hdiff_rhs / mesh%cells(c)%volume
                  end if
               end if
            end if
         end do

         if (use_variable_density_energy) then
            rho_cell = max(transport%rho(c), tiny_safe)
            rho_old_cell = max(transport%rho_old(c), tiny_safe)
            energy%h(c) = (rho_old_cell * energy%h_old(c) + params%dt * &
                          (rhs_h / mesh%cells(c)%volume + energy%qrad(c))) / rho_cell
            energy_update_delta_rate_local = energy_update_delta_rate_local + &
               ((rho_cell * energy%h(c) - rho_old_cell * energy%h_old(c)) * mesh%cells(c)%volume) / params%dt
            energy_update_rhs_integral_local = energy_update_rhs_integral_local + &
               rhs_h + energy%qrad(c) * mesh%cells(c)%volume
            operator_consistent_rho_h_integral_local = operator_consistent_rho_h_integral_local + &
               rho_cell * energy%h(c) * mesh%cells(c)%volume
            if (allocated(energy%operator_consistent_rho_h)) then
               energy%operator_consistent_rho_h(c) = rho_cell * energy%h(c)
            end if
         else
            energy%h(c) = energy%h_old(c) + params%dt * &
                          (rhs_h / mesh%cells(c)%volume + energy%qrad(c) / params%rho)
            energy_update_delta_rate_local = energy_update_delta_rate_local + &
               params%rho * (energy%h(c) - energy%h_old(c)) * mesh%cells(c)%volume / params%dt
            energy_update_rhs_integral_local = energy_update_rhs_integral_local + &
               params%rho * rhs_h + energy%qrad(c) * mesh%cells(c)%volume
            operator_consistent_rho_h_integral_local = operator_consistent_rho_h_integral_local + &
               params%rho * energy%h(c) * mesh%cells(c)%volume
            if (allocated(energy%operator_consistent_rho_h)) then
               energy%operator_consistent_rho_h(c) = params%rho * energy%h(c)
            end if
         end if
      end do
      call profiler_stop('Energy_Flux_Update')

      update_send(1) = energy_update_delta_rate_local
      update_send(2) = energy_update_rhs_integral_local
      update_send(3) = operator_consistent_rho_h_integral_local
      call MPI_Allreduce(update_send, update_recv, 3, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing direct energy update closure')
      energy_update_delta_rate_global = update_recv(1)
      energy_update_rhs_integral_global = update_recv(2)
      operator_consistent_rho_h_integral_global = update_recv(3)
      energy_update_balance_global = energy_update_delta_rate_global - energy_update_rhs_integral_global
      energy_update_denom = max(abs(energy_update_delta_rate_global), abs(energy_update_rhs_integral_global), tiny_safe)
      energy%last_energy_update_delta_rate_integral = energy_update_delta_rate_global
      energy%last_energy_update_rhs_integral = energy_update_rhs_integral_global
      energy%last_energy_update_balance_defect = energy_update_balance_global
      energy%relative_last_energy_update_balance_defect = abs(energy_update_balance_global) / energy_update_denom
      energy%cumulative_energy_update_delta_integral = &
         energy%cumulative_energy_update_delta_integral + params%dt * energy_update_delta_rate_global
      energy%cumulative_energy_update_rhs_integral = &
         energy%cumulative_energy_update_rhs_integral + params%dt * energy_update_rhs_integral_global
      energy%last_operator_consistent_rho_h_integral = operator_consistent_rho_h_integral_global
      energy%operator_consistent_rho_h_available = 1

      energy%last_conductive_boundary_flux_out = conductive_boundary_flux_out_local
      energy%last_conductive_boundary_flux_available = 1

      qrad_integral_step_local = zero
      rho_species_hdiff_integral_step_local = zero
      do c = flow%first_cell, flow%last_cell
         vol = mesh%cells(c)%volume
         if (allocated(energy%qrad)) qrad_integral_step_local = qrad_integral_step_local + energy%qrad(c) * vol
         if (allocated(energy%species_enthalpy_diffusion)) then
            if (use_variable_density_energy) then
               rho_diag = max(transport%rho(c), tiny_safe)
            else
               rho_diag = params%rho
            end if
            rho_species_hdiff_integral_step_local = rho_species_hdiff_integral_step_local + &
                                                    rho_diag * energy%species_enthalpy_diffusion(c) * vol
         end if
      end do

      energy%cumulative_boundary_rho_h_advective_flux_out = &
         energy%cumulative_boundary_rho_h_advective_flux_out + params%dt * advective_boundary_rho_h_flux_out_local
      energy%cumulative_boundary_rho_h_conductive_flux_out = &
         energy%cumulative_boundary_rho_h_conductive_flux_out + params%dt * conductive_boundary_flux_out_local
      energy%cumulative_qrad_integral = energy%cumulative_qrad_integral + params%dt * qrad_integral_step_local
      energy%cumulative_rho_species_hdiff_integral = &
         energy%cumulative_rho_species_hdiff_integral + params%dt * rho_species_hdiff_integral_step_local
      energy%cumulative_energy_budget_available = 1


      call write_variable_density_energy_debug(mesh, flow, params, fields, energy, &
                                                'post_flux_pre_cantera_postsync', species_Y, transport)

      if (associated(hk_species)) deallocate(hk_species)
      if (associated(hk_boundary)) deallocate(hk_boundary)
      if (allocated(sp_diff_flux)) deallocate(sp_diff_flux)
      if (allocated(sp_Y_face_lin)) deallocate(sp_Y_face_lin)
      if (allocated(Y_point)) deallocate(Y_point)

      if (params%enable_cantera_thermo) then
         ! Protect transported h from property-refresh roundoff: recover T from
         ! the newly transported h, refresh cp/lambda/rho_thermo, then restore h.
         energy%h_old = energy%h

         call profiler_start('Energy_Cantera_PostSync')
         call recover_temperature_and_update_thermo_cantera(mesh, params, energy, species_Y)
         call profiler_stop('Energy_Cantera_PostSync')

         energy%h = energy%h_old
      else
         call profiler_start('Energy_ConstantCp_RecoverT')
         call recover_temperature_constant_cp(params, energy)
         call profiler_stop('Energy_ConstantCp_RecoverT')
      end if

      ! Synchronize updated owned cells for output and the next step.
      call profiler_start('Energy_Final_Exchange')
      call flow_exchange_cell_scalar(flow, energy%h)
      call flow_exchange_cell_scalar(flow, energy%T)
      if (allocated(energy%lambda)) call flow_exchange_cell_scalar(flow, energy%lambda)
      if (allocated(energy%species_enthalpy_diffusion)) call flow_exchange_cell_scalar(flow, energy%species_enthalpy_diffusion)
      call profiler_stop('Energy_Final_Exchange')
   end subroutine advance_energy_transport



   !> Targeted diagnostics for experimental variable-density energy/Cantera HP failures.
   subroutine write_variable_density_energy_debug(mesh, flow, params, fields, energy, label, species_Y, transport)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(in) :: energy
      character(len=*), intent(in) :: label
      real(rk), intent(in), optional :: species_Y(:,:)
      type(transport_properties_t), intent(in), optional :: transport

      integer :: c, k, ierr
      integer :: local_bad_cell
      real(rk) :: local_min_h, local_max_h, local_min_T, local_max_T
      real(rk) :: local_min_rho, local_max_rho, local_min_rho_old, local_max_rho_old
      real(rk) :: local_min_S, local_max_S
      real(rk) :: global_min_h, global_max_h, global_min_T, global_max_T
      real(rk) :: global_min_rho, global_max_rho, global_min_rho_old, global_max_rho_old
      real(rk) :: global_min_S, global_max_S
      logical :: have_transport, suspicious

      if (.not. params%enable_variable_density) return
      if (.not. params%variable_density_debug) return
      if (.not. allocated(energy%h) .or. .not. allocated(energy%T)) return

      associate(dummy_mesh_ncells => mesh%ncells)
      end associate

      have_transport = present(transport)

      local_min_h = huge(zero)
      local_max_h = -huge(zero)
      local_min_T = huge(zero)
      local_max_T = -huge(zero)
      local_min_rho = huge(zero)
      local_max_rho = -huge(zero)
      local_min_rho_old = huge(zero)
      local_max_rho_old = -huge(zero)
      local_min_S = huge(zero)
      local_max_S = -huge(zero)
      local_bad_cell = -1

      do c = flow%first_cell, flow%last_cell
         local_min_h = min(local_min_h, energy%h(c))
         local_max_h = max(local_max_h, energy%h(c))
         local_min_T = min(local_min_T, energy%T(c))
         local_max_T = max(local_max_T, energy%T(c))

         if (have_transport) then
            if (allocated(transport%rho)) then
               local_min_rho = min(local_min_rho, transport%rho(c))
               local_max_rho = max(local_max_rho, transport%rho(c))
            end if
            if (allocated(transport%rho_old)) then
               local_min_rho_old = min(local_min_rho_old, transport%rho_old(c))
               local_max_rho_old = max(local_max_rho_old, transport%rho_old(c))
            end if
         end if

         if (allocated(fields%divergence_source)) then
            local_min_S = min(local_min_S, fields%divergence_source(c))
            local_max_S = max(local_max_S, fields%divergence_source(c))
         end if

         if (local_bad_cell < 0) then
            if (energy%h(c) < -1.0e5_rk .or. energy%h(c) > 1.0e7_rk) then
               local_bad_cell = c
            end if
         end if
      end do

      if (.not. have_transport) then
         local_min_rho = zero
         local_max_rho = zero
         local_min_rho_old = zero
         local_max_rho_old = zero
      end if
      if (.not. allocated(fields%divergence_source)) then
         local_min_S = zero
         local_max_S = zero
      end if

      call MPI_Allreduce(local_min_h, global_min_h, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug min_h')
      call MPI_Allreduce(local_max_h, global_max_h, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug max_h')
      call MPI_Allreduce(local_min_T, global_min_T, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug min_T')
      call MPI_Allreduce(local_max_T, global_max_T, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug max_T')
      call MPI_Allreduce(local_min_rho, global_min_rho, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug min_rho')
      call MPI_Allreduce(local_max_rho, global_max_rho, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug max_rho')
      call MPI_Allreduce(local_min_rho_old, global_min_rho_old, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug min_rho_old')
      call MPI_Allreduce(local_max_rho_old, global_max_rho_old, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug max_rho_old')
      call MPI_Allreduce(local_min_S, global_min_S, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug min_S')
      call MPI_Allreduce(local_max_S, global_max_S, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI failure reducing vd debug max_S')

      if (flow%rank == 0) then
         write(output_unit,'(a,a,a,2(1x,es13.5),a,2(1x,es13.5),a,2(1x,es13.5),a,2(1x,es13.5),a,2(1x,es13.5))') &
            'vd-energy-debug[', trim(label), '] h=', global_min_h, global_max_h, &
            ' T=', global_min_T, global_max_T, &
            ' rho=', global_min_rho, global_max_rho, &
            ' rho_old=', global_min_rho_old, global_max_rho_old, &
            ' S=', global_min_S, global_max_S
      end if

      suspicious = local_bad_cell > 0
      if (suspicious) then
         c = local_bad_cell
         write(output_unit,'(a,i0,a,i0,a,a,a,es16.8,a,es16.8)') &
            'vd-energy-bad-cell rank=', flow%rank, ' cell=', c, ' label=', trim(label), &
            ' h=', energy%h(c), ' T=', energy%T(c)

         if (have_transport) then
            if (allocated(transport%rho) .and. allocated(transport%rho_old)) then
               write(output_unit,'(a,es16.8,a,es16.8)') &
                  '  rho=', transport%rho(c), ' rho_old=', transport%rho_old(c)
            end if
         end if

         if (allocated(fields%divergence_source)) then
            write(output_unit,'(a,es16.8)') '  divergence_source=', fields%divergence_source(c)
         end if

         if (present(species_Y)) then
            do k = 1, min(params%nspecies, size(species_Y, 1))
               write(output_unit,'(a,i0,a,es16.8)') '  Y(', k, ')=', species_Y(k, c)
            end do
         end if
      end if

   end subroutine write_variable_density_energy_debug


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

      filename = trim(params%output_dir)//'/diagnostics/energy_diagnostics.csv'
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
      logical, save :: printed_thermo_mode = .false.
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

         if (allocated(energy%T_old)) then
            old_T = energy%T_old(c)
         else
            old_T = params%energy_reference_T + &
                    (energy%h_old(c) - params%energy_reference_h) / params%energy_cp
         end if
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

      filename = trim(params%output_dir)//'/diagnostics/energy_diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='old', position='append', action='write')

      write(unit_id,'(i0,12(",",es16.8))') step, time, global_min_T, global_max_T, mean_T, &
         global_min_h, global_max_h, mean_h, global_min_qrad, global_max_qrad, &
         global_integral_qrad, global_max_delta_T, rel_h_residual

      close(unit_id)
      if (params%enable_cantera_thermo .and. step > 0 .and. .not. printed_thermo_mode) then
         if (params%enable_species) then
            write(output_unit,'(a,i0)') 'Cantera thermo mode: transported species composition, nspecies = ', params%nspecies
         else
            write(output_unit,'(a,a)') 'Cantera thermo mode: default single species = ', trim(params%thermo_default_species)
         end if
         printed_thermo_mode = .true.
      end if

      write(output_unit,'(a,i0,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5)') &
         'energy step=', step, 'Tmin=', global_min_T, 'Tmax=', global_max_T, &
         'Tmean=', mean_T, 'dTmax=', global_max_delta_T, 'rel_h=', rel_h_residual
   end subroutine write_energy_diagnostics_row


   !> Print bridge-level Cantera cache statistics, summed over flow ranks.
   !!
   !! The counters are diagnostic only. They help determine whether Cantera
   !! transport, thermo-sync, and species-enthalpy caches are doing useful work
   !! before attempting further performance changes.
   subroutine write_cantera_cache_stats(flow)
      type(flow_mpi_t), intent(in) :: flow

      integer, parameter :: nstats = 16
      real(c_double) :: local_stats(nstats), global_stats(nstats)
      integer :: ierr

      local_stats = 0.0_c_double
      global_stats = 0.0_c_double

      call cantera_get_cache_stats_c(local_stats, nstats)
      call MPI_Reduce(local_stats, global_stats, nstats, MPI_DOUBLE_PRECISION, MPI_SUM, 0, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('energy', 'MPI_Reduce failed for Cantera cache statistics')

      if (flow%rank /= 0) return
      if (sum(global_stats) <= 0.0_c_double) return

      write(output_unit,'(a)') ' ======================================================================'
      write(output_unit,'(a)') '  CANTERA CACHE STATISTICS'
      write(output_unit,'(a)') '  Counts are summed over flow ranks; hits/misses are cell/point states.'
      write(output_unit,'(a)') ' ======================================================================'
      write(output_unit,'(a)') '                         Cache          Calls         Points           Hits         Misses     Hit%'
      write(output_unit,'(a)') ' ----------------------------------------------------------------------'
      call print_cache_row('Transport_mu_Dk', 1)
      call print_cache_row('Energy_ThermoSync', 5)
      call print_cache_row('SpeciesH_Bulk', 9)
      call print_cache_row('SpeciesH_Point', 13)
      write(output_unit,'(a)') ' ======================================================================'

   contains

      subroutine print_cache_row(name, offset)
         character(len=*), intent(in) :: name
         integer, intent(in) :: offset

         real(c_double) :: calls, points, hits, misses, hit_rate

         calls = global_stats(offset)
         points = global_stats(offset + 1)
         hits = global_stats(offset + 2)
         misses = global_stats(offset + 3)

         if (calls <= 0.0_c_double .and. points <= 0.0_c_double) return

         if (hits + misses > 0.0_c_double) then
            hit_rate = 100.0_c_double * hits / (hits + misses)
         else
            hit_rate = 0.0_c_double
         end if

         write(output_unit,'(2x,a24,1x,f12.0,1x,f14.0,1x,f14.0,1x,f14.0,1x,f8.2)') &
            trim(name), calls, points, hits, misses, hit_rate
      end subroutine print_cache_row

   end subroutine write_cantera_cache_stats

end module mod_energy
