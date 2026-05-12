!> Allocation and management of primary flow variables (U, P, Fluxes).
!!
!! This module acts as the central repository for all physical fields associated 
!! with the fluid flow. It handles the allocation, initialization, and 
!! deallocation of cell-centered and face-centered variables.
!!
!! The fields are designed to support the **Fractional-Step Projection Method**:
!! 1. **`u`**: The soluable, divergence-free velocity field at the new time level.
!! 2. **`u_star`**: The intermediate predicted velocity (contains advection and diffusion).
!! 3. **`phi`**: The scalar potential used to project `u_star` onto a divergence-free space.
!! 4. **`face_flux`**: The mass/volume flux at cell faces, used for conservative transport.
module mod_fields
   use mod_kinds, only : rk, zero
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_bc, only : bc_set_t
   implicit none

   private

   !> Container for all primary hydrodynamic fields.
   type, public :: flow_fields_t
      real(rk), allocatable :: u(:,:)        !< Current cell-centered velocity vector $(u, v, w)$ [m/s].
      real(rk), allocatable :: u_old(:,:)    !< Velocity vector from the previous time step $n$ [m/s].
      real(rk), allocatable :: u_star(:,:)   !< Intermediate predicted velocity field [m/s].

      real(rk), allocatable :: p(:)          !< Static pressure field $P$ [Pa].
      real(rk), allocatable :: phi(:)        !< Pressure correction potential $\phi$. Used in the Poisson solver.
      real(rk), allocatable :: div(:)        !< Local velocity divergence $\nabla \cdot \mathbf{u}$. Should be $\approx 0$ after projection.

      !> Conservative face-centered flux $U_f = \mathbf{u}_f \cdot \mathbf{n}_f$.
      !! Oriented according to the face normal (owner $\rightarrow$ neighbor).
      real(rk), allocatable :: face_flux(:)  !< Volumetric flux across faces [m^3/s].

      !> Storage for the previous explicit RHS (Advection + Diffusion).
      !! Required for 2nd-order Adams-Bashforth (AB2) time marching.
      real(rk), allocatable :: rhs_old(:,:)  !< Momentum RHS from step $n-1$.
      logical :: rhs_old_valid = .false.     !< False on the first step (triggers Euler fallback).
   end type flow_fields_t

   public :: allocate_fields, finalize_fields, initialize_fields

contains

   !> Dynamically allocates all arrays within the flow fields container.
   !!
   !! Sizing is determined by the number of cells and faces in the provided mesh.
   !! All fields are initialized to zero upon allocation.
   !!
   !! @param mesh The mesh structure defining the domain size.
   !! @param fields The container to be allocated.
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


   !> Initializes flow fields and sets simulation initial conditions.
   !!
   !! Currently implements a "Start from Rest" condition where all 
   !! velocities and pressures are zero. IC overrides for restarts or 
   !! analytical solutions will be implemented here.
   !!
   !! @param mesh The computational mesh.
   !! @param fields The fields container to initialize.
   subroutine initialize_fields(mesh, fields)
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(inout) :: fields

      ! Initial condition: Fluids starts from rest. 
      ! Boundary motion and body forces drive the flow from t=0.
      call allocate_fields(mesh, fields)

   end subroutine initialize_fields


   !> Deallocates all arrays and resets validity flags.
   !!
   !! @param fields The container to be cleared.
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