!> Core flow solver field containers.
!!
!! This module defines the data structure for cell-centered and face-centered
!! flow variables (velocity, pressure, fluxes, etc.).
module mod_fields
   use mod_kinds, only : rk, zero
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_bc, only : bc_set_t
   implicit none

   private

   !> Data structure for primary flow variables.
   type, public :: flow_fields_t
      real(rk), allocatable :: u(:,:)        !< Cell-centered velocity (3, ncells)
      real(rk), allocatable :: u_old(:,:)    !< Previous-step velocity (3, ncells)
      real(rk), allocatable :: u_star(:,:)   !< Intermediate (predicted) velocity (3, ncells)

      real(rk), allocatable :: p(:)          !< Pressure field (ncells)
      real(rk), allocatable :: phi(:)        !< Pressure correction/potential (ncells)
      real(rk), allocatable :: div(:)        !< Divergence of velocity field (ncells)

      !> Conservative face flux, oriented with mesh%faces(f)%normal.
      !! For internal faces this is owner -> neighbor.
      real(rk), allocatable :: face_flux(:)  !< Face-centered mass/volume flux (nfaces)

      !> Previous explicit momentum RHS for AB2 time integration.
      real(rk), allocatable :: rhs_old(:,:)  !< Previous momentum source term (3, ncells)
      logical :: rhs_old_valid = .false.     !< Flag for AB2 initialization
   end type flow_fields_t

   public :: allocate_fields, finalize_fields, initialize_fields

contains

   !> Allocate all flow field arrays.
   !!
   !! @param mesh Mesh data structure for sizing.
   !! @param fields Flow fields structure to allocate.
   subroutine allocate_fields(mesh, fields)
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(inout) :: fields

      call finalize_fields(fields)

      allocate(fields%u(3, mesh%ncells))
      allocate(fields%u_old(3, mesh%ncells))
      allocate(fields%u_star(3, mesh%ncells))

      allocate(fields%p(mesh%ncells))
      allocate(fields%phi(mesh%ncells))
      allocate(fields%div(mesh%ncells))

      allocate(fields%face_flux(mesh%nfaces))
      allocate(fields%rhs_old(3, mesh%ncells))

      fields%u = zero
      fields%u_old = zero
      fields%u_star = zero

      fields%p = zero
      fields%phi = zero
      fields%div = zero

      fields%face_flux = zero
      fields%rhs_old = zero
      fields%rhs_old_valid = .false.
   end subroutine allocate_fields


   !> Initialize flow fields and set initial conditions.
   !!
   !! @param mesh Mesh data structure.
   !! @param params Case parameters.
   !! @param bc Boundary condition set.
   !! @param fields Flow fields structure to initialize.
   subroutine initialize_fields(mesh, params, bc, fields)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(bc_set_t), intent(in) :: bc
      type(flow_fields_t), intent(inout) :: fields

      ! v1 starts from rest. Boundary motion and body force enter through
      ! the FV operators.
      call allocate_fields(mesh, fields)

      associate(dummy => params%nsteps + bc%npatches)
      end associate
   end subroutine initialize_fields


   subroutine finalize_fields(fields)
      type(flow_fields_t), intent(inout) :: fields

      if (allocated(fields%u)) deallocate(fields%u)
      if (allocated(fields%u_old)) deallocate(fields%u_old)
      if (allocated(fields%u_star)) deallocate(fields%u_star)

      if (allocated(fields%p)) deallocate(fields%p)
      if (allocated(fields%phi)) deallocate(fields%phi)
      if (allocated(fields%div)) deallocate(fields%div)

      if (allocated(fields%face_flux)) deallocate(fields%face_flux)
      if (allocated(fields%rhs_old)) deallocate(fields%rhs_old)

      fields%rhs_old_valid = .false.
   end subroutine finalize_fields

end module mod_fields