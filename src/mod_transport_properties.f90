!> Fluid transport properties and Cantera bridge abstraction.
!!
!! This module holds cell-centered physical properties such as density, 
!! viscosity, and species diffusivity. It acts as an isolation layer 
!! between the flow/species solvers and the future Cantera C++ bridge.
module mod_transport_properties
   use mod_kinds, only : rk, zero, fatal_error
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t

   implicit none
   private

   public :: transport_properties_t
   public :: initialize_transport
   public :: finalize_transport
   public :: update_transport_properties

   !> Transport property fields available for all cells.
   type transport_properties_t
      real(rk), allocatable :: rho(:)          !< Cell density [kg/m^3]
      real(rk), allocatable :: mu(:)           !< Cell dynamic viscosity [Pa*s]
      real(rk), allocatable :: nu(:)           !< Cell kinematic viscosity [m^2/s]
      real(rk), allocatable :: lambda(:)       !< Cell thermal conductivity [W/m*K]
      real(rk), allocatable :: diffusivity(:,:) !< Cell species diffusivity (nspecies, ncells) [m^2/s]
   end type transport_properties_t

   !> C-Binding stub for the future Cantera bridge.
   !! The implementation of this will eventually reside in a C++ library.
   interface
      function cantera_get_species_count_c() bind(c, name="cantera_get_species_count_c")
         use iso_c_binding, only : c_int
         integer(c_int) :: cantera_get_species_count_c
      end function cantera_get_species_count_c

      subroutine cantera_get_species_name_c(k, name_out, name_len) bind(c, name="cantera_get_species_name_c")
         use iso_c_binding, only : c_int, c_char
         integer(c_int), value :: k
         character(kind=c_char), intent(out) :: name_out(*)
         integer(c_int), value :: name_len
      end subroutine cantera_get_species_name_c

      subroutine cantera_init_c(mech_file, nspecies, species_names_flat, name_len) bind(c, name="cantera_init_c")
         use iso_c_binding, only : c_char, c_int
         character(kind=c_char), intent(in) :: mech_file(*)
         integer(c_int), value :: nspecies
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_init_c

      subroutine cantera_update_transport_c(ncells, T, P, nspecies, Y_in, mu_out, diff_out, species_names_flat, name_len) bind(c, name="cantera_update_transport_c")
         use iso_c_binding, only : c_int, c_double, c_char
         integer(c_int), value :: ncells
         real(c_double), intent(in) :: T(ncells), P(ncells)
         integer(c_int), value :: nspecies
         real(c_double), intent(in) :: Y_in(*)
         real(c_double), intent(out) :: mu_out(ncells), diff_out(*)
         character(kind=c_char), intent(in) :: species_names_flat(*)
         integer(c_int), value :: name_len
      end subroutine cantera_update_transport_c
   end interface

contains

   !> Allocate transport property arrays.
   !!
   !! @param mesh Full replicated mesh.
   !! @param params Case configuration parameters.
   !! @param transport Transport property data structure to initialize.
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

   !> Wrapper to pack strings and call cantera_init_c.
   !!
   !! @param params Case configuration parameters containing mechanism file and species names.
   subroutine initialize_cantera_wrapper(params)
      use iso_c_binding, only : c_char, c_null_char
      type(case_params_t), intent(inout) :: params
      character(kind=c_char, len=len(params%cantera_mech_file)+1) :: c_mech_file
      character(kind=c_char), allocatable :: c_names_flat(:)
      integer :: k, n_len, c_nsp

      c_mech_file = trim(params%cantera_mech_file) // c_null_char
      
      ! Initial dummy call to setup mechanism
      call cantera_init_c(c_mech_file, 0, "", 0)

      if (params%enable_reactions) then
         ! Automatic Mode: Inherit all species from Cantera
         c_nsp = cantera_get_species_count_c()
         params%nspecies = c_nsp
         n_len = len(params%species_name(1))
         do k = 1, params%nspecies
            call cantera_get_species_name_c(k-1, params%species_name(k), n_len)
         end do
      end if

      ! Re-initialize with the finalized species list
      n_len = len(params%species_name(1))
      allocate(c_names_flat(n_len * params%nspecies))
      c_names_flat = " "
      do k = 1, params%nspecies
         do c_nsp = 1, n_len
            c_names_flat((k-1)*n_len + c_nsp) = params%species_name(k)(c_nsp:c_nsp)
         end do
      end do

      call cantera_init_c(c_mech_file, params%nspecies, c_names_flat, n_len)
      deallocate(c_names_flat)
   end subroutine initialize_cantera_wrapper

   !> Deallocate transport property arrays.
   !!
   !! @param transport Transport property data structure to finalize.
   subroutine finalize_transport(transport)
      type(transport_properties_t), intent(inout) :: transport

      if (allocated(transport%rho)) deallocate(transport%rho)
      if (allocated(transport%mu)) deallocate(transport%mu)
      if (allocated(transport%nu)) deallocate(transport%nu)
      if (allocated(transport%lambda)) deallocate(transport%lambda)
      if (allocated(transport%diffusivity)) deallocate(transport%diffusivity)
   end subroutine finalize_transport

   !> Update transport fields (rho, mu, nu, diffusivity).
   !!
   !! If enable_cantera is true, calls the C wrapper to evaluate properties dynamically.
   !! Otherwise, acts as a constant fallback mechanism relying on 
   !! parameters loaded directly from the case config.
   !!
   !! @param mesh Full replicated mesh.
   !! @param params Case configuration parameters.
   !! @param Y Optional species mass fractions (required if Cantera is enabled).
   !! @param transport Transport property data structure to update.
   subroutine update_transport_properties(mesh, params, Y, transport)
      use iso_c_binding, only : c_char
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      real(rk), intent(in), optional :: Y(:,:)
      type(transport_properties_t), intent(inout) :: transport
      integer :: k, n_len
      real(rk), allocatable :: T_arr(:), P_arr(:)
      character(kind=c_char, len=len(params%species_name(1))*params%nspecies) :: c_names_flat

      ! Constant density is assumed for the incompressible formulation.
      transport%rho = params%rho
      
      if (params%enable_cantera_fluid .or. params%enable_cantera_species) then
         ! Build temperature and pressure arrays
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
            call fatal_error('cantera', 'Y array not provided to update_transport_properties but cantera is enabled')
         end if

         ! Apply overrides if specific Cantera features are disabled
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
         ! Constant uniform properties
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
