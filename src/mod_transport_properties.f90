!> Fluid transport properties and Cantera C++ bridge abstraction.
!!
!! This module manages cell-centered physical properties such as flow density,
!! viscosity, and species diffusivity. It acts as an isolation layer
!! between the high-level flow/species solvers and the external Cantera 
!! C++ library.
!!
!! The module supports two primary modes:
!! 1. **Constant Fallback**: Properties are pulled directly from `case.nml` 
!!    and remain uniform and static.
!! 2. **Cantera Dynamic**: Mixture-averaged transport properties \(\mu\), \(D_k\)
!!    are evaluated at every update interval using the current local 
!!    composition, temperature, and pressure via the `cantera_interface`.
!!    The flow/projection density remains `params%rho`; Cantera thermodynamic
!!    density is stored separately as `energy%rho_thermo` for diagnostics and
!!    future use.
module mod_transport_properties
   use mod_kinds, only : rk, zero, one, fatal_error
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_mpi_flow, only : flow_mpi_t, flow_exchange_cell_scalar, flow_exchange_cell_matrix

   implicit none
   private

   public :: transport_properties_t
   public :: initialize_transport
   public :: finalize_transport
   public :: update_transport_properties

   !> Container for cell-centered transport properties.
   !!
   !! These fields are used by the Momentum and Species transport solvers 
   !! to evaluate diffusion terms and calculate the local Reynolds/Peclet 
   !! numbers.
   type transport_properties_t
      real(rk), allocatable :: rho(:)          !< Constant flow/projection density [kg/m^3]; not Cantera rho_thermo.
      real(rk), allocatable :: mu(:)           !< Dynamic viscosity \(\mu\) [Pa*s].
      real(rk), allocatable :: nu(:)           !< Kinematic viscosity \(\nu = \mu/\rho\) [m^2/s].
      real(rk), allocatable :: lambda(:)       !< Thermal conductivity \(\lambda\) [W/m*K].
      real(rk), allocatable :: diffusivity(:,:) !< Species-specific diffusivity \(D_{k,i}\) [m^2/s]. Indexed as (species, cell).
   end type transport_properties_t

   !> C-Binding interface for the Cantera C++ bridge.
   !!
   !! These routines map directly to the symbols defined in `cantera_interface.cpp`.
   interface
      !> Retrieves the total number of species in the loaded Cantera mechanism.
      function cantera_get_species_count_c() bind(c, name="cantera_get_species_count_c")
         use iso_c_binding, only : c_int
         integer(c_int) :: cantera_get_species_count_c
      end function cantera_get_species_count_c

      !> Retrieves the name of a specific species from Cantera.
      subroutine cantera_get_species_name_c(k, name_out, name_len) bind(c, name="cantera_get_species_name_c")
         use iso_c_binding, only : c_int, c_char
         integer(c_int), value :: k !< 0-indexed species index.
         character(kind=c_char), intent(out) :: name_out(*) !< Output buffer for species name.
         integer(c_int), value :: name_len !< Length of the output buffer.
      end subroutine cantera_get_species_name_c

      !> Initializes the Cantera mixture and loads the mechanism file.
      subroutine cantera_init_c(mech_file, nspecies, species_names_flat, name_len) bind(c, name="cantera_init_c")
         use iso_c_binding, only : c_char, c_int
         character(kind=c_char), intent(in) :: mech_file(*) !< Null-terminated mechanism file path.
         integer(c_int), value :: nspecies !< Number of species to track.
         character(kind=c_char), intent(in) :: species_names_flat(*) !< Flattened array of species names.
         integer(c_int), value :: name_len !< Width of each name field in the flattened array.
      end subroutine cantera_init_c

      !> Performs a bulk update of transport properties for all cells.
      subroutine cantera_update_transport_c(ncells, T, P, nspecies, Y_in, mu_out, diff_out, species_names_flat, name_len) bind(c, name="cantera_update_transport_c")
         use iso_c_binding, only : c_int, c_double, c_char
         integer(c_int), value :: ncells !< Total number of computational cells.
         real(c_double), intent(in) :: T(ncells), P(ncells) !< Local Temperature [K] and Pressure [Pa].
         integer(c_int), value :: nspecies !< Number of species.
         real(c_double), intent(in) :: Y_in(*) !< Input mass fractions (nspecies, ncells).
         real(c_double), intent(out) :: mu_out(ncells), diff_out(*) !< Output viscosity and diffusivity.
         character(kind=c_char), intent(in) :: species_names_flat(*) !< Reference names for validation.
         integer(c_int), value :: name_len !< Width of name fields.
      end subroutine cantera_update_transport_c
   end interface

contains

   !> Allocates property arrays and initializes the Cantera bridge if required.
   !!
   !! @param mesh The computational mesh.
   !! @param params Simulation parameters (contains mechanism and species info).
   !! @param transport The property container to initialize.
   subroutine initialize_transport(mesh, params, transport)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(inout) :: params
      type(transport_properties_t), intent(inout) :: transport

      if (params%enable_cantera_fluid .or. params%enable_cantera_species .or. params%enable_cantera_thermo) then
         call initialize_cantera_wrapper(params)
      end if

      allocate(transport%rho(mesh%ncells))
      allocate(transport%mu(mesh%ncells))
      allocate(transport%nu(mesh%ncells))
      allocate(transport%lambda(mesh%ncells))
      
      if (params%nspecies > 0) then
         allocate(transport%diffusivity(params%nspecies, mesh%ncells))
         transport%diffusivity = zero
      end if

      transport%rho = params%rho
      transport%mu = zero
      transport%nu = zero
      transport%lambda = zero
   end subroutine initialize_transport


   !> Higher-level wrapper to coordinate Cantera initialization.
   !!
   !! Handles string conversion between Fortran and C, and manages the 
   !! automatic discovery of species counts and names from the mechanism file.
   subroutine initialize_cantera_wrapper(params)
      use iso_c_binding, only : c_char, c_null_char
      type(case_params_t), intent(inout) :: params
      character(kind=c_char, len=len(params%cantera_mech_file)+1) :: c_mech_file
      character(kind=c_char), allocatable :: c_names_flat(:)
      integer :: k, n_len, c_nsp

      c_mech_file = trim(params%cantera_mech_file) // c_null_char
      
      ! Initial dummy call to setup mechanism and query species
      call cantera_init_c(c_mech_file, 0, "", 0)

      if (params%enable_reactions) then
         ! Automatic Discovery Mode: Inherit all species from the mechanism
         c_nsp = cantera_get_species_count_c()
         params%nspecies = c_nsp
         n_len = len(params%species_name(1))
         do k = 1, params%nspecies
            call cantera_get_species_name_c(k-1, params%species_name(k), n_len)
         end do
      else if (params%enable_cantera_thermo .and. params%nspecies == 0) then
         ! No-species thermo mode: use user-selected inert/default species.
         params%nspecies = 1
         params%species_name(1) = trim(params%thermo_default_species)
      else if (params%enable_cantera_fluid .and. params%nspecies == 0) then
         ! Single Species Fluid Mode: default to N2 when no species list is provided.
         params%nspecies = 1
         params%species_name(1) = 'N2'
      end if

      ! Re-initialize with the finalized species list for future transport updates
      n_len = len(params%species_name(1))
      allocate(c_names_flat(n_len * params%nspecies))
      c_names_flat = " "
      do k = 1, params%nspecies
         do c_nsp = 1, n_len
            c_names_flat((k-1)*n_len + c_nsp) = char(ichar(params%species_name(k)(c_nsp:c_nsp)))
         end do
      end do

      call cantera_init_c(c_mech_file, params%nspecies, c_names_flat, n_len)
      deallocate(c_names_flat)
   end subroutine initialize_cantera_wrapper


   !> Safely deallocates transport property fields.
   subroutine finalize_transport(transport)
      type(transport_properties_t), intent(inout) :: transport

      if (allocated(transport%rho)) deallocate(transport%rho)
      if (allocated(transport%mu)) deallocate(transport%mu)
      if (allocated(transport%nu)) deallocate(transport%nu)
      if (allocated(transport%lambda)) deallocate(transport%lambda)
      if (allocated(transport%diffusivity)) deallocate(transport%diffusivity)
   end subroutine finalize_transport


   !> Evaluates physical properties for the entire domain.
   !!
   !! If Cantera transport is enabled, it calls the C++ bridge to calculate
   !! mixture-averaged viscosity and diffusivity. Otherwise, it applies the
   !! uniform constant values defined in `params`. In both cases,
   !! `transport%rho` is the constant flow/projection density `params%rho`.
   !!
   !! @param mesh The computational mesh.
   !! @param flow MPI decomposition and halo metadata.
   !! @param params Simulation parameters and background conditions.
   !! @param Y Optional input mass fractions (ncells, nspecies).
   !! @param transport The property container to update.
   subroutine update_transport_properties(mesh, flow, params, Y, transport, T_state)
      use iso_c_binding, only : c_char
      use mod_profiler, only : profiler_start, profiler_stop
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(case_params_t), intent(in) :: params
      real(rk), intent(in), optional :: Y(:,:)
      type(transport_properties_t), intent(inout) :: transport
      real(rk), intent(in), optional :: T_state(:)
      integer :: c, k, i, n_len, nloc
      real(rk), allocatable :: T_arr(:), P_arr(:), Y_local(:,:), Y_pure(:,:)
      real(rk), allocatable :: mu_local(:), diff_local(:,:)
      character(kind=c_char, len=len(params%species_name(1))*params%nspecies) :: c_names_flat

      if (mesh%ncells /= size(transport%rho)) then
         call fatal_error('transport', 'transport field size does not match mesh')
      end if

      ! Constant density is assumed for the standard incompressible solver mode.
      ! This is not a Cantera density update; energy%rho_thermo is diagnostic only.
      transport%rho = params%rho
      
      if (params%enable_cantera_fluid .or. params%enable_cantera_species) then
         call profiler_start('Transport_Setup')

         ! Build owned-cell Temperature and thermodynamic-pressure arrays for the bridge call.
         nloc = flow%nlocal
         allocate(T_arr(nloc))
         allocate(P_arr(nloc))
         allocate(mu_local(nloc))
         allocate(diff_local(max(1, params%nspecies), nloc))
         if (present(T_state)) then
            if (size(T_state) < mesh%ncells) then
               call fatal_error('transport', 'T_state has incompatible size for Cantera transport update')
            end if
            do i = 1, nloc
               c = flow%first_cell + i - 1
               T_arr(i) = T_state(c)
            end do
         else
            T_arr = params%background_temp
         end if
         P_arr = params%background_press
         mu_local = zero
         diff_local = zero

         c_names_flat = ""
         n_len = len(params%species_name(1))
         do k = 1, params%nspecies
            c_names_flat((k-1)*n_len+1:k*n_len) = params%species_name(k)
         end do

         if (present(Y)) then
            allocate(Y_local(max(1, params%nspecies), nloc))
            do i = 1, nloc
               c = flow%first_cell + i - 1
               do k = 1, params%nspecies
                  Y_local(k, i) = Y(k, c)
               end do
            end do

            call profiler_stop('Transport_Setup')
            call profiler_start('Transport_Cantera_Call')
            call cantera_update_transport_c(nloc, T_arr, P_arr, params%nspecies, &
                                            Y_local, mu_local, diff_local, &
                                            c_names_flat, n_len)
            call profiler_stop('Transport_Cantera_Call')

            call profiler_start('Transport_Cleanup')
            deallocate(Y_local)
            call profiler_stop('Transport_Cleanup')
         else
            ! Case where species are not yet initialized: pass pure background gas (Y=1.0)
            allocate(Y_pure(max(1, params%nspecies), nloc))
            Y_pure = one
            call profiler_stop('Transport_Setup')
            call profiler_start('Transport_Cantera_Call')
            call cantera_update_transport_c(nloc, T_arr, P_arr, params%nspecies, &
                                            Y_pure, mu_local, diff_local, &
                                            c_names_flat, n_len)
            call profiler_stop('Transport_Cantera_Call')

            call profiler_start('Transport_Cleanup')
            deallocate(Y_pure)
            call profiler_stop('Transport_Cleanup')
         end if

         call profiler_start('Transport_Unpack')
         do i = 1, nloc
            c = flow%first_cell + i - 1
            transport%mu(c) = mu_local(i)
            do k = 1, params%nspecies
               transport%diffusivity(k, c) = diff_local(k, i)
            end do
         end do
         call profiler_stop('Transport_Unpack')

         call profiler_start('Transport_Exchange')

         ! Apply partial overrides (e.g., if viscosity is constant but diffusivity is dynamic)
         if (.not. params%enable_cantera_fluid) then
            transport%nu = params%nu
            transport%mu = params%rho * params%nu
         else
            do c = flow%first_cell, flow%last_cell
               transport%nu(c) = transport%mu(c) / transport%rho(c)
            end do
            call profiler_start('Transport_Exchange_Mu')
            call flow_exchange_cell_scalar(flow, transport%mu)
            call profiler_stop('Transport_Exchange_Mu')

            call profiler_start('Transport_Exchange_Nu')
            call flow_exchange_cell_scalar(flow, transport%nu)
            call profiler_stop('Transport_Exchange_Nu')
         end if

         if (.not. params%enable_cantera_species .and. params%nspecies > 0) then
            do k = 1, params%nspecies
               transport%diffusivity(k, :) = params%species_diffusivity(k)
            end do
         else if (params%nspecies > 0) then
            call profiler_start('Transport_Exchange_Diff')
            call flow_exchange_cell_matrix(flow, transport%diffusivity)
            call profiler_stop('Transport_Exchange_Diff')
         end if
         call profiler_stop('Transport_Exchange')
         
         call profiler_start('Transport_Cleanup')
         deallocate(T_arr)
         deallocate(P_arr)
         deallocate(mu_local)
         deallocate(diff_local)
         call profiler_stop('Transport_Cleanup')
      else
         ! Constant Fallback Mode
         transport%nu = params%nu
         transport%mu = params%rho * params%nu
         
         if (params%nspecies > 0) then
            do k = 1, params%nspecies
               transport%diffusivity(k, :) = params%species_diffusivity(k)
            end do
         end if
      end if
      
   
      ! Validation-control option: keep flow viscosity/nu constant unless
      ! variable viscosity is explicitly enabled. This still allows Cantera
      ! species diffusivity D_k to be used when species Cantera transport is on.
      if (.not. params%enable_variable_nu) then
         transport%mu = params%rho * params%nu
      end if

end subroutine update_transport_properties

end module mod_transport_properties
