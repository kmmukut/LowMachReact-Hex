!> Fluid transport properties and Cantera C++ bridge abstraction.
!!
!! This module manages cell-centered physical properties such as density, 
!! viscosity, and species diffusivity. It acts as an isolation layer 
!! between the high-level flow/species solvers and the external Cantera 
!! C++ library.
!!
!! The module supports two primary modes:
!! 1. **Constant Fallback**: Properties are pulled directly from `case.nml` 
!!    and remain uniform and static.
!! 2. **Cantera Dynamic**: Mixture-averaged properties ($\mu$, $\rho$, $D_k$) 
!!    are evaluated at every update interval using the current local 
!!    composition, temperature, and pressure via the `cantera_interface`.
module mod_transport_properties
   use mod_kinds, only : rk, zero, one, fatal_error
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t

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
      real(rk), allocatable :: rho(:)          !< Mixture density $\rho$ [kg/m^3].
      real(rk), allocatable :: mu(:)           !< Dynamic viscosity $\mu$ [Pa*s].
      real(rk), allocatable :: nu(:)           !< Kinematic viscosity $\nu = \mu/\rho$ [m^2/s].
      real(rk), allocatable :: lambda(:)       !< Thermal conductivity $\lambda$ [W/m*K].
      real(rk), allocatable :: diffusivity(:,:) !< Species-specific diffusivity $D_{k,i}$ [m^2/s]. Indexed as (species, cell).
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

      if (params%enable_cantera_fluid .or. params%enable_cantera_species) then
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

      transport%rho = zero
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
      else if (params%enable_cantera_fluid .and. params%nspecies == 0) then
         ! Single Species Fluid Mode: Default to the first species (e.g., N2)
         params%nspecies = 1
         n_len = len(params%species_name(1))
         call cantera_get_species_name_c(0, params%species_name(1), n_len)
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
   !! If Cantera is enabled, it calls the C++ bridge to calculate 
   !! mixture-averaged viscosity and diffusivity. Otherwise, it applies 
   !! the uniform constant values defined in `params`.
   !!
   !! @param mesh The computational mesh.
   !! @param params Simulation parameters and background conditions.
   !! @param Y Optional input mass fractions (ncells, nspecies).
   !! @param transport The property container to update.
   subroutine update_transport_properties(mesh, params, Y, transport)
      use iso_c_binding, only : c_char
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      real(rk), intent(in), optional :: Y(:,:)
      type(transport_properties_t), intent(inout) :: transport
      integer :: k, n_len
      real(rk), allocatable :: T_arr(:), P_arr(:), Y_pure(:,:)
      character(kind=c_char, len=len(params%species_name(1))*params%nspecies) :: c_names_flat

      ! Constant density is assumed for the standard incompressible solver mode.
      transport%rho = params%rho
      
      if (params%enable_cantera_fluid .or. params%enable_cantera_species) then
         ! Build local Temperature and Pressure arrays for the bridge call
         allocate(T_arr(mesh%ncells))
         allocate(P_arr(mesh%ncells))
         T_arr = params%background_temp
         P_arr = params%background_press

         c_names_flat = ""
         n_len = len(params%species_name(1))
         do k = 1, params%nspecies
            c_names_flat((k-1)*n_len+1:k*n_len) = params%species_name(k)
         end do

         if (present(Y)) then
            call cantera_update_transport_c(mesh%ncells, T_arr, P_arr, params%nspecies, &
                                            Y, transport%mu, transport%diffusivity, &
                                            c_names_flat, n_len)
         else
            ! Case where species are not yet initialized: pass pure background gas (Y=1.0)
            allocate(Y_pure(1, mesh%ncells))
            Y_pure = one
            call cantera_update_transport_c(mesh%ncells, T_arr, P_arr, params%nspecies, &
                                            Y_pure, transport%mu, transport%diffusivity, &
                                            c_names_flat, n_len)
            deallocate(Y_pure)
         end if

         ! Apply partial overrides (e.g., if viscosity is constant but diffusivity is dynamic)
         if (.not. params%enable_cantera_fluid) then
            transport%nu = params%nu
            transport%mu = params%rho * params%nu
         else
            transport%nu = transport%mu / transport%rho
         end if

         if (.not. params%enable_cantera_species .and. params%nspecies > 0) then
            do k = 1, params%nspecies
               transport%diffusivity(k, :) = params%species_diffusivity(k)
            end do
         end if
         
         deallocate(T_arr)
         deallocate(P_arr)
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
      
   end subroutine update_transport_properties

end module mod_transport_properties
