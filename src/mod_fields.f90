!> Allocation and management of primary flow variables and projection fluxes.
!!
!! This module is the central owner of hydrodynamic state used by the
!! fractional-step solver.  Arrays are globally sized because every MPI rank
!! stores the replicated mesh, but numerical kernels update only owned cells or
!! owner-owned faces before halo exchange.
!!
!! Important flux convention: `fields%face_flux(f)` is the volumetric face flux
!! \(u_f \cdot n_f A_f\) in \( \mathrm{m^3/s} \), oriented from
!! `mesh%faces(f)%owner` to `mesh%faces(f)%neighbor`.  Cell-local loops
!! reorient it to outward-from-cell by changing the sign when the current cell
!! is the neighbor.  `fields%mass_flux(f)` has the same owner-to-neighbor
!! orientation and stores \( \rho_f u_f \cdot n_f A_f \) in
!! \( \mathrm{kg/s} \).
!!
!! Constant-density mode targets \( \nabla \cdot u = 0 \).  Guarded
!! variable-density low-Mach mode targets \( \nabla \cdot u = S \), where
!! `fields%projection_divergence_source` records the source actually consumed
!! by the most recent projection and `fields%divergence_source` may later be
!! advanced by energy/thermo density updates for the next projection.
module mod_fields
   use mod_kinds, only : rk, zero
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_bc, only : bc_set_t
   implicit none

   private

   !> Container for all primary hydrodynamic fields.
   type, public :: flow_fields_t
      real(rk), allocatable :: u(:,:)        !< Current cell-centered velocity vector \((u, v, w)\) [m/s].
      real(rk), allocatable :: u_old(:,:)    !< Velocity vector from the previous time step \(n\) [m/s].
      real(rk), allocatable :: u_star(:,:)   !< Intermediate predicted velocity field [m/s].

      real(rk), allocatable :: p(:)          !< Projection/hydrodynamic pressure-like field [Pa], not Cantera thermodynamic pressure.
      real(rk), allocatable :: phi(:)        !< Projection correction potential \(\phi\) used by the Poisson solve.
      real(rk), allocatable :: div(:)        !< Recomputed \( \nabla \cdot u \) [1/s]; compare with `projection_divergence_source` in variable-density mode.

      !> Conservative face-centered volumetric flux \(U_f = \mathbf{u}_f \cdot \mathbf{n}_f A_f\).
      !! Oriented according to the face normal (owner \(\rightarrow\) neighbor).
      real(rk), allocatable :: face_flux(:)  !< Volumetric flux across faces [m^3/s].

      !> Conservative face-centered mass flux \(\dot{m}_f = \rho_f U_f\).
      !! Oriented according to the face normal (owner \(\rightarrow\) neighbor).
      real(rk), allocatable :: mass_flux(:)  !< Mass flux across faces [kg/s].

      !> Current low-Mach divergence source \(S\) [1/s] prepared for the next projection.
      !!
      !! Constant-density mode keeps this zero.  Variable-density mode updates
      !! it after the energy/thermo density sync from a conservative density
      !! source; it can therefore differ from the source used by the latest
      !! projection.
      real(rk), allocatable :: divergence_source(:)  !< Target divergence source [1/s].
      real(rk), allocatable :: projection_rho(:) !< Active density snapshot at projection start [kg/m^3].

      !> Copy of the low-Mach divergence source actually consumed by the latest projection.
      !! This distinguishes the projection residual `div(u)-S_projection` from
      !! post-energy source-evolution diagnostics using `divergence_source`.
      real(rk), allocatable :: projection_divergence_source(:)  !< Projection-time target divergence source [1/s].

      !> Storage for the previous explicit RHS (Advection + Diffusion).
      !! Required for 2nd-order Adams-Bashforth (AB2) time marching.
      real(rk), allocatable :: rhs_old(:,:)  !< Momentum RHS from step \(n-1\).
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
      allocate(fields%divergence_source(mesh%ncells))
      allocate(fields%projection_rho(mesh%ncells))
      allocate(fields%projection_divergence_source(mesh%ncells))
      allocate(fields%mass_flux(mesh%nfaces))
      allocate(fields%rhs_old(3, mesh%ncells))

      fields%u = zero
      fields%u_old = zero
      fields%u_star = zero

      fields%p = zero
      fields%phi = zero
      fields%div = zero

      fields%face_flux = zero
      fields%divergence_source = zero
      fields%projection_rho = zero
      fields%projection_divergence_source = zero
      fields%mass_flux = zero
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
      if (allocated(fields%divergence_source)) deallocate(fields%divergence_source)
      if (allocated(fields%projection_rho)) deallocate(fields%projection_rho)
      if (allocated(fields%projection_divergence_source)) deallocate(fields%projection_divergence_source)
      if (allocated(fields%mass_flux)) deallocate(fields%mass_flux)
      if (allocated(fields%rhs_old)) deallocate(fields%rhs_old)

      fields%rhs_old_valid = .false.
   end subroutine finalize_fields

end module mod_fields
