!> Fluid transport properties and Cantera C++ bridge abstraction.
!!
!! This module manages cell-centered physical properties such as flow density,
!! viscosity, and species diffusivity. It acts as an isolation layer
!! between the high-level flow/species solvers and the external Cantera 
!! C++ library.
!!
!! The active density contract is deliberately explicit:
!! - Constant-density mode keeps `transport%rho = params%rho`; any
!!   `energy%rho_thermo` from Cantera is diagnostic.
!! - Guarded variable-density mode with `density_eos="cantera"` copies
!!   `energy%rho_thermo` into `transport%rho`, so density comes from the
!!   selected Cantera phase and its YAML-defined EOS.
!!
!! Cantera transport updates mixture viscosity \(\mu\) and species
!! diffusivities \(D_k\) from \(T\), \(p_0\), and \(Y\).  The transport cache
!! key is \(T,p_0,Y\); the projection pressure field is not used as Cantera
!! thermodynamic pressure.
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
   public :: update_transport_properties, sync_active_density_from_thermo

   !> Container for cell-centered transport properties.
   !!
   !! These fields are used by the Momentum and Species transport solvers 
   !! to evaluate diffusion terms and calculate the local Reynolds/Peclet 
   !! numbers.
   type transport_properties_t
      real(rk), allocatable :: rho(:)          !< Active flow/projection density [kg/m^3]: `params%rho` or synced Cantera `rho_thermo`.
      real(rk), allocatable :: rho_old(:)  !< Previous active flow density [kg/m^3] for low-Mach source and conservative updates.
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
      subroutine cantera_init_c(mech_file, phase_name, nspecies, species_names_flat, name_len) bind(c, name="cantera_init_c")
         use iso_c_binding, only : c_char, c_int
         character(kind=c_char), intent(in) :: mech_file(*) !< Null-terminated mechanism file path.
         character(kind=c_char), intent(in) :: phase_name(*) !< Null-terminated optional Cantera phase name.
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

      allocate(transport%rho_old(mesh%ncells))
      allocate(transport%mu(mesh%ncells))
      allocate(transport%nu(mesh%ncells))
      allocate(transport%lambda(mesh%ncells))
      
      if (params%nspecies > 0) then
         allocate(transport%diffusivity(params%nspecies, mesh%ncells))
         transport%diffusivity = zero
      end if

      transport%rho = params%rho

      if (allocated(transport%rho_old)) transport%rho_old = transport%rho
      ! Initialize the flow transport fields to the constant fallback so
      ! step-0 VTU diagnostics are meaningful before the first dynamic
      ! Cantera transport refresh. Dynamic modes overwrite these values
      ! during update_transport_properties.
      transport%nu = params%nu
      transport%mu = params%rho * params%nu
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
      character(kind=c_char, len=len(params%cantera_phase_name)+1) :: c_phase_name
      character(kind=c_char), allocatable :: c_names_flat(:)
      integer :: k, n_len, c_nsp

      c_mech_file = trim(params%cantera_mech_file) // c_null_char
      c_phase_name = trim(params%cantera_phase_name) // c_null_char
      
      ! Initial dummy call to setup mechanism and query species
      call cantera_init_c(c_mech_file, c_phase_name, 0, "", 0)

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

      call cantera_init_c(c_mech_file, c_phase_name, params%nspecies, c_names_flat, n_len)
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


   !> Evaluates viscosity and species diffusivity for the active configuration.
   !!
   !! If Cantera transport is enabled, it calls the C++ bridge to calculate
   !! mixture-averaged viscosity and diffusivity using \(T\), \(Y\), and
   !! `params%background_press` as \(p_0\). Otherwise, it applies the uniform
   !! constants from `params`. Constant-density mode resets `transport%rho` to
   !! `params%rho`; variable-density mode preserves the active density synced
   !! from `energy%rho_thermo`.
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

      ! Constant-density mode owns active density through params%rho.  In
      ! variable-density mode, active rho is owned by the thermo-density sync and
      ! must not be overwritten during mu/D_k refresh.
      if (.not. params%enable_variable_density) transport%rho = params%rho
      
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

   !> Synchronize active flow density from Cantera thermodynamic density.
   !!
   !! In the default constant-density mode this returns without modifying
   !! `transport`; `transport%rho` remains `params%rho` and `rho_thermo` is
   !! diagnostic.  In guarded variable-density mode,
   !! `density_eos = "cantera"` means the Cantera thermo density
   !! `rho_thermo` becomes the active flow density for projection, face mass
   !! flux, species transport, and conservative `rho*h` energy transport.
   !!
   !! The routine only copies owned cells, then exchanges the replicated scalar
   !! field so all ranks have a globally valid `transport%rho` array.
   subroutine sync_active_density_from_thermo(mesh, flow, params, transport, rho_thermo)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(case_params_t), intent(in) :: params
      type(transport_properties_t), intent(inout) :: transport
      real(rk), intent(in) :: rho_thermo(:)

      integer :: c

      if (.not. params%enable_variable_density) return

      if (trim(params%density_eos) /= 'cantera') then
         call fatal_error('transport', 'variable-density sync currently supports density_eos="cantera" only')
      end if

      if (.not. allocated(transport%rho)) then
         call fatal_error('transport', 'transport rho must be allocated before density sync')
      end if

      ! rho_old stores the active density before the Cantera thermo-density overwrite.
      ! It is used by the staged low-Mach divergence-source update.
      if (.not. allocated(transport%rho_old)) then
         allocate(transport%rho_old(mesh%ncells))
         transport%rho_old = transport%rho
      end if

      if (size(rho_thermo) < mesh%ncells) then
         call fatal_error('transport', 'rho_thermo has incompatible size for active-density sync')
      end if

      do c = 1, mesh%ncells
         if (.not. flow%owned(c)) cycle

         if (rho_thermo(c) <= 0.0_rk) then
            call fatal_error('transport', 'non-positive rho_thermo in active-density sync')
         end if

         transport%rho_old(c) = transport%rho(c)
         transport%rho(c) = rho_thermo(c)

         if (allocated(transport%nu)) then
            if (params%enable_variable_nu .and. allocated(transport%mu)) then
               transport%nu(c) = transport%mu(c) / transport%rho(c)
            else
               transport%nu(c) = params%nu
            end if
         end if
      end do

      call flow_exchange_cell_scalar(flow, transport%rho)
      if (allocated(transport%rho_old)) call flow_exchange_cell_scalar(flow, transport%rho_old)
      if (allocated(transport%nu)) call flow_exchange_cell_scalar(flow, transport%nu)

   end subroutine sync_active_density_from_thermo


end module mod_transport_properties
