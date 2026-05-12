!> Incompressible Navier-Stokes solver using the Fractional-Step Projection method.
!!
!! This module implements the core hydrodynamic solver for the 
!! LowMachReact-Hex framework. It solves the incompressible Navier-Stokes 
!! equations using a semi-implicit fractional-step approach:
!!
!! 1. **Predictor Step**: An intermediate velocity $\mathbf{u}^*$ is 
!!    calculated by advancing the momentum equation explicitly, excluding 
!!     the new pressure gradient.
!!    $$\frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} = -(\mathbf{u} \cdot \nabla) \mathbf{u} + \nu \nabla^2 \mathbf{u} - \frac{1}{\rho} \nabla p^n + \mathbf{f}$$
!!
!! 2. **Poisson Solve**: A pressure correction potential $\phi = p^{n+1} - p^n$ 
!!    is found by solving the Poisson equation derived from the 
!!    continuity constraint ($\nabla \cdot \mathbf{u}^{n+1} = 0$).
!!    $$\nabla^2 \phi = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*$$
!!
!! 3. **Corrector Step**: The final velocity $\mathbf{u}^{n+1}$ and 
!!    pressure $p^{n+1}$ are updated using the potential gradient.
!!    $$\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla \phi$$
!!    $$p^{n+1} = p^n + \phi$$
!!
!! The linear system for the Poisson equation is solved using a 
!! **Preconditioned Conjugate Gradient (PCG)** method with a 
!! diagonal (Jacobi) preconditioner.
module mod_flow_projection
   use mpi_f08
   use mod_kinds, only : rk, zero, one, half, tiny_safe, fatal_error
   use mod_input, only : case_params_t
   use mod_mesh_types, only : mesh_t
   use mod_mpi_flow, only : flow_mpi_t, flow_exchange_cell_scalar, &
                            flow_exchange_cell_matrix, flow_exchange_face_scalar, &
                            flow_global_dot_owned, flow_global_two_dots_owned, &
                            flow_global_max_owned
   use mod_bc, only : bc_set_t, bc_periodic, bc_neumann, bc_dirichlet, &
                      patch_type_for_face, boundary_velocity, &
                      face_effective_neighbor, boundary_pressure_type, &
                      boundary_pressure
   use mod_profiler, only : profiler_start, profiler_stop
   use mod_fields, only : flow_fields_t
   use mod_transport_properties, only : transport_properties_t
   implicit none

   private

   !> Solver diagnostics and global physics statistics.
   !!
   !! Populated at the end of every timestep to monitor convergence 
   !! and physical integrity (e.g., mass balance, energy decay).
   type, public :: solver_stats_t
      integer :: pressure_iterations = 0 !< Number of iterations taken by the most recent PCG solve.
      integer :: pressure_iterations_total = 0 !< Cumulative PCG iterations over all pressure solves.
      integer :: pressure_iterations_max = 0 !< Maximum PCG iterations used by any pressure solve.
      integer :: pressure_solve_count = 0 !< Number of pressure solves completed.
      real(rk) :: pressure_iterations_avg = zero !< Average PCG iterations per pressure solve.
      real(rk) :: pressure_residual = zero !< Final $L_2$ residual norm of the Poisson system.
      real(rk) :: max_divergence = zero    !< Maximum local velocity divergence $\nabla \cdot \mathbf{u}$ [1/s].
      real(rk) :: rms_divergence = zero    !< Root-mean-square divergence across the domain.
      real(rk) :: net_boundary_flux = zero !< Global mass flux imbalance across all boundaries [kg/s].
      real(rk) :: kinetic_energy = zero    !< Total domain kinetic energy $\int \frac{1}{2} \rho |\mathbf{u}|^2 dV$ [J].
      real(rk) :: cfl = zero               !< Current Courant-Friedrichs-Lewy number.
      real(rk) :: wall_time = zero         !< Cumulative solver wall-clock time [s].
      real(rk) :: max_velocity = zero      !< Maximum velocity magnitude $|\mathbf{u}|_{max}$ [m/s].
      real(rk) :: total_mass = zero        !< Total integrated mass in the domain [kg].
      real(rk) :: min_species_y = zero     !< Global minimum species mass fraction (diagnostic for boundedness).
   end type solver_stats_t

   !> Cached sparse matrix coefficients for the Pressure Poisson operator.
   !!
   !! Since the mesh geometry and BC types are static, the Laplacian 
   !! coefficients are computed once at startup to accelerate the PCG matvec.
   type :: pressure_operator_cache_t
      logical :: initialized = .false.      !< Initialization toggle.
      logical :: has_dirichlet_pressure = .false. !< True if at least one patch has a fixed pressure.
      logical :: has_neumann_outlet = .false. !< True if at least one boundary patch can absorb flux imbalance.
      integer :: ncells = 0                 !< Number of cells.
      integer :: max_faces = 0              !< Max faces per cell.
      integer, allocatable :: nb(:,:)       !< Neighbor indices for each cell face.
      real(rk), allocatable :: coeff(:,:)   !< Off-diagonal Laplacian coefficients $A_{ij} = \frac{Area}{dist}$.
      real(rk), allocatable :: diag(:)      !< Diagonal Laplacian coefficients $A_{ii} = \sum A_{ij}$.
   end type pressure_operator_cache_t

   !> Temporary workspace for projection step calculations.
   !!
   !! Holds intermediate vectors (RHS, residuals, search directions) 
   !! to avoid repeated allocation during the simulation loop.
   type :: projection_workspace_t
      logical :: initialized = .false.
      integer :: ncells = 0
      integer :: nfaces = 0

      real(rk), allocatable :: local_vec(:,:)      !< Local momentum RHS.
      real(rk), allocatable :: local_vec_star(:,:) !< Local intermediate velocity.
      real(rk), allocatable :: local_scalar(:)     !< Local divergence scalar.
      real(rk), allocatable :: rhs_poisson(:)      !< Source term for Poisson solve.
      real(rk), allocatable :: local_face_flux(:)  !< Face fluxes on owned faces.
      real(rk), allocatable :: predicted_face_flux(:) !< Global intermediate face fluxes.

      real(rk), allocatable :: r(:)                !< PCG residual vector.
      real(rk), allocatable :: z(:)                !< PCG preconditioned residual.
      real(rk), allocatable :: pvec(:)             !< PCG search direction.
      real(rk), allocatable :: ap(:)               !< PCG operator-applied vector.
      real(rk), allocatable :: local_ap(:)         !< PCG local operator-applied vector.
   end type projection_workspace_t

   type(pressure_operator_cache_t), save :: pressure_cache
   type(projection_workspace_t), save :: projection_work

   public :: advance_projection_step, compute_flow_diagnostics
   public :: compute_and_update_cfl
   public :: finalize_flow_projection_workspace, face_normal_distance

contains

   !> Orchestrates one full fractional-step iteration.
   !!
   !! This routine implements the high-level logic of the projection method. 
   !! It coordinates the explicit momentum update, the elliptic pressure 
   !! solve, and the final divergence-free correction.
   !!
   !! @param mesh The computational mesh.
   !! @param flow MPI decomposition context.
   !! @param bc Boundary condition set.
   !! @param params Case parameters (dt, rho, etc).
   !! @param transport Physical property fields.
   !! @param fields Flow field containers to be updated.
   !! @param stats Diagnostic stats to be populated.
   subroutine advance_projection_step(mesh, flow, bc, params, transport, fields, stats)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(transport_properties_t), intent(in) :: transport
      type(flow_fields_t), intent(inout) :: fields
      type(solver_stats_t), intent(inout) :: stats

      integer :: c
      integer, save :: projection_step_counter = 0
      logical :: do_projection_diagnostics

      projection_step_counter = projection_step_counter + 1
      do_projection_diagnostics = (mod(projection_step_counter, params%output_interval) == 0) .or. &
                                  (projection_step_counter == params%nsteps)

      call ensure_projection_workspace(mesh)
      call ensure_pressure_operator_cache(mesh, flow, bc)

      associate( &
         local_vec => projection_work%local_vec, &
         local_scalar => projection_work%local_scalar, &
         rhs_poisson => projection_work%rhs_poisson, &
         local_face_flux => projection_work%local_face_flux, &
         local_vec_star => projection_work%local_vec_star, &
         predicted_face_flux => projection_work%predicted_face_flux)

         ! 1. Predictor: Compute Explicit momentum RHS and Advance intermediate velocity u*
         call profiler_start('Projection_Momentum_RHS')
         call compute_momentum_rhs(mesh, flow, bc, params, transport, fields%u, fields%p, local_vec)
         call profiler_stop('Projection_Momentum_RHS')
         
         ! 2. Advance intermediate velocity u* locally
         ! This uses the PREVIOUS fields%rhs_old and the CURRENT local_vec
         call profiler_start('Projection_AB2')
         call advance_ab2(mesh, flow, params, fields, local_vec, local_vec_star)
         call profiler_stop('Projection_AB2')
         
         ! Store current RHS for next step (AB2)
         ! rhs_old is effectively partitioned; each rank only needs its owned cell values.
         fields%rhs_old = local_vec
         fields%rhs_old_valid = .true.

         fields%u_star = local_vec_star

         call profiler_start('Projection_Predict_Flux')

         call profiler_start('PredictFlux_Exchange_UStar')
         call flow_exchange_cell_matrix(flow, fields%u_star)
         call profiler_stop('PredictFlux_Exchange_UStar')

         ! 3. Interpolate predicted cell velocity to intermediate face fluxes
         call profiler_start('PredictFlux_Compute')
         call compute_predicted_face_flux(mesh, flow, bc, fields%u_star, predicted_face_flux)
         call profiler_stop('PredictFlux_Compute')

         ! 3b. Balance mass at open boundaries to prevent drift in closed-loop systems
         call profiler_start('PredictFlux_Balance')
         call balance_neumann_outlet_flux(mesh, flow, bc, predicted_face_flux)
         call profiler_stop('PredictFlux_Balance')

         call profiler_start('PredictFlux_Exchange_Face')
         call flow_exchange_face_scalar(flow, predicted_face_flux)
         call profiler_stop('PredictFlux_Exchange_Face')

         call profiler_stop('Projection_Predict_Flux')

         ! 4. Construct the Poisson RHS: b = -rho/dt * div(u*)
         call profiler_start('Projection_Poisson_RHS')
         call compute_flux_divergence(mesh, flow, predicted_face_flux, local_scalar)
         fields%div = local_scalar

         rhs_poisson = zero
         do c = flow%first_cell, flow%last_cell
            rhs_poisson(c) = -params%rho / params%dt * &
                             fields%div(c) * mesh%cells(c)%volume
         end do

         ! Floating pressure handle: pin first cell if no Dirichlet BC exists
         if (.not. pressure_cache%has_dirichlet_pressure) then
            rhs_poisson(1) = zero
         end if
         fields%phi = zero
         call profiler_stop('Projection_Poisson_RHS')

         ! Solve Laplacian system for potential phi
         call profiler_start('Projection_PCG')
         call solve_pressure_correction(mesh, flow, bc, params, rhs_poisson, fields%phi, stats)
         call profiler_stop('Projection_PCG')

         ! 5. Update Pressure: p(n+1) = p(n) + phi
         call profiler_start('Projection_Pressure_Update')
         do c = flow%first_cell, flow%last_cell
            fields%p(c) = fields%p(c) + fields%phi(c)
         end do
         if (.not. pressure_cache%has_dirichlet_pressure) then
            fields%p(1) = zero
         end if
         call flow_exchange_cell_scalar(flow, fields%p)
         call profiler_stop('Projection_Pressure_Update')

         ! 6. Correct face fluxes and cell-centered velocity locally
         call profiler_start('Projection_Correction')
         call correct_face_flux(mesh, flow, bc, params, predicted_face_flux, fields%phi, local_face_flux)
         fields%face_flux = local_face_flux
         call flow_exchange_face_scalar(flow, fields%face_flux)
         call correct_cell_velocity(mesh, flow, bc, params, fields%u_star, fields%phi, local_vec)

         ! 7. Keep velocity ghosts current for the next predictor/species step.
         fields%u = local_vec
         call flow_exchange_cell_matrix(flow, fields%u)

         ! Divergence is only needed for diagnostics/output, so avoid computing
         ! and exchanging it on non-output steps.
         if (do_projection_diagnostics) then
            call compute_flux_divergence(mesh, flow, fields%face_flux, local_scalar)
            fields%div = local_scalar
            call flow_exchange_cell_scalar(flow, fields%div)
         end if
         call profiler_stop('Projection_Correction')

         if (do_projection_diagnostics) then
            call profiler_start('Projection_Diagnostics')
            call compute_flow_diagnostics(mesh, flow, bc, params, fields, stats)
            call profiler_stop('Projection_Diagnostics')
         end if

      end associate
   end subroutine advance_projection_step


   !> Advances velocity using the 2nd-order Adams-Bashforth scheme.
   !!
   !! Falls back to 1st-order Forward Euler on the very first timestep 
   !! when previous RHS data is unavailable.
   subroutine advance_ab2(mesh, flow, params, fields, rhs, local_ustar)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      real(rk), intent(in) :: rhs(:,:)
      real(rk), intent(out) :: local_ustar(:,:)

      integer :: c

      local_ustar = zero

      do c = flow%first_cell, flow%last_cell
         if (fields%rhs_old_valid) then
            local_ustar(:, c) = fields%u(:, c) + params%dt * &
               (1.5_rk * rhs(:, c) - 0.5_rk * fields%rhs_old(:, c))
         else
            local_ustar(:, c) = fields%u(:, c) + params%dt * rhs(:, c)
         end if
      end do

      associate(dummy => mesh%ncells)
      end associate
   end subroutine advance_ab2


   !> Evaluates the advective, diffusive, and pressure terms of the momentum equation.
   !!
   !! Supports both **Upwind** (stable) and **Central** (high-accuracy) 
   !! convection schemes. Advection is computed using a flux-form 
   !! discretization.
   subroutine compute_momentum_rhs(mesh, flow, bc, params, transport, u, p, local_rhs)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(transport_properties_t), intent(in) :: transport
      real(rk), intent(in) :: u(:,:)
      real(rk), intent(in) :: p(:)
      real(rk), intent(out) :: local_rhs(:,:)

      integer :: c, lf, f, nb
      real(rk) :: nvec(3), nhat(3)
      real(rk) :: uf(3), ub(3), advected(3)
      real(rk) :: un_area
      real(rk) :: conv(3), diff(3), gradp(3)
      real(rk) :: dist, diff_face
      logical :: use_central

      use_central = .false.

      select case (trim(params%convection_scheme))
      case ('central', 'central_difference', 'central-difference')
         use_central = .true.
      case ('upwind')
         use_central = .false.
      case default
         call fatal_error('flow', 'unknown convection_scheme '//trim(params%convection_scheme))
      end select

      local_rhs = zero

      do c = flow%first_cell, flow%last_cell
         conv = zero
         diff = zero

         ! Compute local pressure gradient at cell center
         call pressure_gradient_cell(mesh, bc, p, c, gradp)

         do lf = 1, mesh%ncell_faces(c)
            f = mesh%cell_faces(lf, c)

            nvec = outward_normal(mesh, f, c)
            nhat = nvec

            nb = face_effective_neighbor(mesh, bc, f, c)

            if (nb > 0) then
               uf = face_linear_vector(mesh, bc, f, c, nb, u(:, c), u(:, nb))
               ub = u(:, nb)
               dist = face_normal_distance(mesh, bc, f, c, nb)
            else
               call boundary_velocity(mesh, bc, f, u(:, c), ub)
               uf = ub
               dist = face_normal_distance(mesh, bc, f, c, 0)
            end if

            un_area = dot_product(uf, nhat) * mesh%faces(f)%area

            ! Scheme selection for advection
            if (use_central) then
               advected = uf
            else
               if (un_area >= zero) then
                  advected = u(:, c)
               else
                  advected = ub
               end if
            end if

            conv = conv + un_area * advected

            ! Viscous diffusion Term
            if (nb > 0) then
               diff_face = half * (transport%nu(c) + transport%nu(nb))
            else
               diff_face = transport%nu(c)
            end if

            diff = diff + diff_face * mesh%faces(f)%area * (ub - u(:, c)) / dist
         end do

         ! Assemble total RHS
         local_rhs(:, c) = -conv / mesh%cells(c)%volume + &
                            diff / mesh%cells(c)%volume + &
                            params%body_force - gradp / params%rho
      end do
   end subroutine compute_momentum_rhs


   !> Calculates the pressure gradient at a cell center using Gauss's Theorem.
   !!
   !! $$\nabla p = \frac{1}{V_c} \int_S p \mathbf{n} dS$$
   subroutine pressure_gradient_cell(mesh, bc, p, cell_id, gradp)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(in) :: p(:)
      integer, intent(in) :: cell_id
      real(rk), intent(out) :: gradp(3)

      integer :: lf, f, nb
      real(rk) :: nvec(3)
      real(rk) :: pf
      logical :: is_dirichlet

      gradp = zero

      do lf = 1, mesh%ncell_faces(cell_id)
         f = mesh%cell_faces(lf, cell_id)
         nvec = outward_normal(mesh, f, cell_id)

         nb = face_effective_neighbor(mesh, bc, f, cell_id)

         if (nb > 0) then
            pf = face_linear_scalar(mesh, bc, f, cell_id, nb, p(cell_id), p(nb))
         else
            if (boundary_pressure_type(mesh, bc, f) == bc_dirichlet) then
               call boundary_pressure(mesh, bc, f, p(cell_id), pf, is_dirichlet)
            else
               pf = p(cell_id)
            end if
         end if

         gradp = gradp + pf * nvec * mesh%faces(f)%area
      end do

      gradp = gradp / mesh%cells(cell_id)%volume
   end subroutine pressure_gradient_cell


   !> Linearly interpolates cell-centered intermediate velocity to mesh faces.
   subroutine compute_predicted_face_flux(mesh, flow, bc, ustar, face_flux)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(in) :: ustar(:,:)
      real(rk), intent(out) :: face_flux(:)

      integer :: i, f, owner, nb
      real(rk) :: uf(3), ub(3)

      face_flux = zero

      do i = 1, size(flow%owned_faces)
         f = flow%owned_faces(i)
         owner = mesh%faces(f)%owner

         nb = face_effective_neighbor(mesh, bc, f, owner)

         if (nb > 0) then
            uf = face_linear_vector(mesh, bc, f, owner, nb, ustar(:, owner), ustar(:, nb))
         else
            call boundary_velocity(mesh, bc, f, ustar(:, owner), ub)
            uf = ub
         end if

         face_flux(f) = dot_product(uf, mesh%faces(f)%normal) * mesh%faces(f)%area
      end do
   end subroutine compute_predicted_face_flux


   !> Corrects face fluxes using the pressure potential gradient.
   subroutine correct_face_flux(mesh, flow, bc, params, predicted_flux, phi, face_flux)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: predicted_flux(:)
      real(rk), intent(in) :: phi(:)
      real(rk), intent(out) :: face_flux(:)

      integer :: i, f, owner, nb
      real(rk) :: dist

      face_flux = zero

      do i = 1, size(flow%owned_faces)
         f = flow%owned_faces(i)
         owner = mesh%faces(f)%owner

         nb = face_effective_neighbor(mesh, bc, f, owner)

         face_flux(f) = predicted_flux(f)

         if (nb > 0) then
            dist = face_normal_distance(mesh, bc, f, owner, nb)

            face_flux(f) = predicted_flux(f) - params%dt / params%rho * &
                            mesh%faces(f)%area * (phi(nb) - phi(owner)) / dist
         else
            if (boundary_pressure_type(mesh, bc, f) == bc_dirichlet) then
               dist = face_normal_distance(mesh, bc, f, owner, 0)
               face_flux(f) = predicted_flux(f) - params%dt / params%rho * &
                               mesh%faces(f)%area * (zero - phi(owner)) / dist
            end if
         end if
      end do
   end subroutine correct_face_flux


   !> Computes the discrete divergence of face fluxes for each cell.
   subroutine compute_flux_divergence(mesh, flow, face_flux, local_div)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: face_flux(:)
      real(rk), intent(out) :: local_div(:)

      integer :: c, lf, f
      real(rk) :: sgn

      local_div = zero

      do c = flow%first_cell, flow%last_cell
         do lf = 1, mesh%ncell_faces(c)
            f = mesh%cell_faces(lf, c)

            if (mesh%faces(f)%owner == c) then
               sgn = one
            else
               sgn = -one
            end if

            local_div(c) = local_div(c) + sgn * face_flux(f) / mesh%cells(c)%volume
         end do
      end do
   end subroutine compute_flux_divergence


   !> Iteratively solves the Pressure Poisson system using PCG.
   !!
   !! The solver uses a Diagonal Preconditioner to improve convergence 
   !! speed on hexahedral meshes. For purely Neumann problems (closed boxes), 
   !! the first cell is pinned to zero to remove the singularity.
   subroutine solve_pressure_correction(mesh, flow, bc, params, rhs, phi, stats)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: rhs(:)
      real(rk), intent(inout) :: phi(:)
      type(solver_stats_t), intent(inout) :: stats

      real(rk) :: rr, rr_new
      real(rk) :: rz, rz_new
      real(rk) :: pap
      real(rk) :: alpha, beta
      integer :: iter
      integer :: c

      call ensure_projection_workspace(mesh)
      call ensure_pressure_operator_cache(mesh, flow, bc)

      associate( &
         r => projection_work%r, &
         z => projection_work%z, &
         pvec => projection_work%pvec, &
         ap => projection_work%ap, &
         local_ap => projection_work%local_ap)

         r = zero
         z = zero
         pvec = zero
         local_ap = zero

         ! The projection driver resets phi to zero before every pressure solve.
         ! Therefore A*phi is exactly zero for the initial PCG residual, so avoid
         ! one unnecessary halo exchange and pressure matvec per timestep.
         local_ap = zero

         do c = flow%first_cell, flow%last_cell
            r(c) = rhs(c)
         end do

         if (.not. pressure_cache%has_dirichlet_pressure) then
            r(1) = zero
            phi(1) = zero
         end if

         do c = flow%first_cell, flow%last_cell
            if (c == 1 .and. .not. pressure_cache%has_dirichlet_pressure) then
               z(c) = zero
            else
               z(c) = r(c) / max(pressure_cache%diag(c), tiny_safe)
            end if
         end do
          
         block
            real(rk) :: res(2)
            call flow_global_two_dots_owned(flow, r, r, r, z, res)
            rr = res(1); rz = res(2)
         end block

         pvec = z
         if (.not. pressure_cache%has_dirichlet_pressure) then
            pvec(1) = zero
         end if

         stats%pressure_iterations = 0
         stats%pressure_residual = sqrt(rr / max(real(mesh%ncells, rk), one))

         do iter = 1, params%pressure_max_iter
            if (stats%pressure_residual <= params%pressure_tol) exit

            ! a. matvec: q = A * p
            call flow_exchange_cell_scalar(flow, pvec)
            call pressure_matvec(mesh, flow, bc, pvec, local_ap)

            ! b. step length: alpha = (r, z) / (p, A*p)
            pap = flow_global_dot_owned(flow, pvec, local_ap)

            if (abs(pap) <= tiny_safe) exit

            alpha = rz / pap

            ! c. update solution and residual: phi = phi + alpha*p, r = r - alpha*q
            do c = flow%first_cell, flow%last_cell
               phi(c) = phi(c) + alpha * pvec(c)
               r(c) = r(c) - alpha * local_ap(c)
            end do

            if (.not. pressure_cache%has_dirichlet_pressure) then
               phi(1) = zero
               r(1) = zero
            end if

            do c = flow%first_cell, flow%last_cell
               if (c == 1 .and. .not. pressure_cache%has_dirichlet_pressure) then
                  z(c) = zero
               else
                  z(c) = r(c) / max(pressure_cache%diag(c), tiny_safe)
               end if
            end do

            block
               real(rk) :: res(2)
               call flow_global_two_dots_owned(flow, r, r, r, z, res)
               rr_new = res(1); rz_new = res(2)
            end block

            stats%pressure_residual = sqrt(rr_new / max(real(mesh%ncells, rk), one))
            stats%pressure_iterations = iter

            if (stats%pressure_residual <= params%pressure_tol) exit

            beta = rz_new / max(rz, tiny_safe)

            do c = flow%first_cell, flow%last_cell
               pvec(c) = z(c) + beta * pvec(c)
            end do
            if (.not. pressure_cache%has_dirichlet_pressure) then
               pvec(1) = zero
            end if

            rr = rr_new
            rz = rz_new
         end do

         if (.not. pressure_cache%has_dirichlet_pressure) then
            phi(1) = zero
         end if

         ! Accumulate pressure-solver iteration diagnostics.
         ! pressure_iterations remains the most recent solve count.
         stats%pressure_solve_count = stats%pressure_solve_count + 1
         stats%pressure_iterations_total = stats%pressure_iterations_total + stats%pressure_iterations
         stats%pressure_iterations_max = max(stats%pressure_iterations_max, stats%pressure_iterations)
         stats%pressure_iterations_avg = real(stats%pressure_iterations_total, rk) / &
                                         max(real(stats%pressure_solve_count, rk), one)

         call flow_exchange_cell_scalar(flow, phi)

      end associate
   end subroutine solve_pressure_correction


   !> Sparse Matrix-Vector multiplication for the Laplacian operator.
   !!
   !! This routine uses the `pressure_cache` to perform efficient 
   !! off-diagonal products on the MPI-owned partition.
   subroutine pressure_matvec(mesh, flow, bc, x, local_ax)
      use mod_profiler, only : profiler_start, profiler_stop
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(in) :: x(:)
      real(rk), intent(out) :: local_ax(:)

      integer :: c, lf, nb

      call profiler_start('Pressure_Matvec')
      call ensure_pressure_operator_cache(mesh, flow, bc)

      local_ax = zero

      do c = flow%first_cell, flow%last_cell
         if (c == 1 .and. .not. pressure_cache%has_dirichlet_pressure) then
            local_ax(c) = x(c)
            cycle
         end if

         local_ax(c) = pressure_cache%diag(c) * x(c)

         do lf = 1, mesh%ncell_faces(c)
            nb = pressure_cache%nb(lf, c)
            if (nb <= 0) cycle

            local_ax(c) = local_ax(c) - pressure_cache%coeff(lf, c) * x(nb)
         end do
      end do
      call profiler_stop('Pressure_Matvec')
   end subroutine pressure_matvec


   !> Updates cell-centered velocity using the pressure potential gradient.
   subroutine correct_cell_velocity(mesh, flow, bc, params, ustar, phi, local_u)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: ustar(:,:)
      real(rk), intent(in) :: phi(:)
      real(rk), intent(out) :: local_u(:,:)

      integer :: c, lf, f, nb
      real(rk) :: nvec(3), grad(3), phi_face

      local_u = zero

      do c = flow%first_cell, flow%last_cell
         grad = zero

         do lf = 1, mesh%ncell_faces(c)
            f = mesh%cell_faces(lf, c)
            nvec = outward_normal(mesh, f, c)

            nb = face_effective_neighbor(mesh, bc, f, c)

            if (nb > 0) then
               phi_face = face_linear_scalar(mesh, bc, f, c, nb, phi(c), phi(nb))
            else
               if (boundary_pressure_type(mesh, bc, f) == bc_dirichlet) then
                  phi_face = zero
               else
                  phi_face = phi(c)
               end if
            end if

            grad = grad + phi_face * nvec * mesh%faces(f)%area
         end do

         grad = grad / mesh%cells(c)%volume

         local_u(:, c) = ustar(:, c) - params%dt / params%rho * grad
      end do
   end subroutine correct_cell_velocity


   !> Aggregates global residuals, kinetic energy, and mass balance data.
   subroutine compute_flow_diagnostics(mesh, flow, bc, params, fields, stats)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(solver_stats_t), intent(inout) :: stats

      real(rk) :: local_ke, global_ke
      real(rk) :: local_rms, global_rms
      real(rk) :: local_flux, global_flux
      integer :: c, ierr

      stats%max_divergence = flow_global_max_owned(flow, fields%div)

      local_rms = zero
      local_ke = zero

      do c = flow%first_cell, flow%last_cell
         local_rms = local_rms + fields%div(c) * fields%div(c)

         local_ke = local_ke + half * params%rho * mesh%cells(c)%volume * &
                    dot_product(fields%u(:, c), fields%u(:, c))
      end do

      call MPI_Allreduce(local_rms, global_rms, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'diagnostic rms')

      call MPI_Allreduce(local_ke, global_ke, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'diagnostic kinetic energy')

      stats%rms_divergence = sqrt(global_rms / max(real(mesh%ncells, rk), one))
      stats%kinetic_energy = global_ke

      call compute_boundary_flux(mesh, flow, bc, fields%face_flux, local_flux)

      call MPI_Allreduce(local_flux, global_flux, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'diagnostic boundary flux')

      stats%net_boundary_flux = global_flux
   end subroutine compute_flow_diagnostics


   !> Adjusts flux at Neumann outlets to ensure strict global mass balance.
   !!
   !! This is critical for solvers without a pressure-pinned cell, 
   !! as numerical drift in the Poisson RHS can lead to a singular system.
   subroutine balance_neumann_outlet_flux(mesh, flow, bc, face_flux)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(inout) :: face_flux(:)

      integer :: f, owner, btype, ierr
      real(rk) :: local_net_flux
      real(rk) :: global_net_flux
      real(rk) :: local_outlet_area
      real(rk) :: global_outlet_area
      real(rk) :: correction_per_area

      if (pressure_cache%has_dirichlet_pressure) return

      ! If the case has no Neumann outlet faces, there is nothing to balance.
      ! This avoids two unnecessary MPI_Allreduce calls per timestep for closed
      ! domains such as lid-driven cavity.
      if (.not. pressure_cache%has_neumann_outlet) return

      local_net_flux = zero
      local_outlet_area = zero

      do f = 1, mesh%nfaces
         if (mesh%faces(f)%neighbor /= 0) cycle

         btype = patch_type_for_face(mesh, bc, f)

         if (btype == bc_periodic) cycle

         owner = mesh%faces(f)%owner
         if (.not. flow%owned(owner)) cycle

         local_net_flux = local_net_flux + face_flux(f)

         if (btype == bc_neumann) then
            local_outlet_area = local_outlet_area + mesh%faces(f)%area
         end if
      end do

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_net_flux, global_net_flux, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'outlet balance net flux')

      call MPI_Allreduce(local_outlet_area, global_outlet_area, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'outlet balance outlet area')
      call profiler_stop('MPI_Communication')

      if (global_outlet_area <= tiny_safe) then
         if (abs(global_net_flux) > tiny_safe) then
            call fatal_error('flow', 'nonzero boundary flux but no neumann outlet faces found')
         end if
         return
      end if

      correction_per_area = -global_net_flux / global_outlet_area

      do f = 1, mesh%nfaces
         if (mesh%faces(f)%neighbor /= 0) cycle

         btype = patch_type_for_face(mesh, bc, f)
         if (btype /= bc_neumann) cycle
         owner = mesh%faces(f)%owner
         if (.not. flow%owned(owner)) cycle

         face_flux(f) = face_flux(f) + correction_per_area * mesh%faces(f)%area
      end do
   end subroutine balance_neumann_outlet_flux


   !> Helper to integrated boundary flux on MPI-owned partitions.
   subroutine compute_boundary_flux(mesh, flow, bc, face_flux, local_flux)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(in) :: face_flux(:)
      real(rk), intent(out) :: local_flux

      integer :: f, owner, btype

      local_flux = zero

      do f = 1, mesh%nfaces
         if (mesh%faces(f)%neighbor /= 0) cycle

         btype = patch_type_for_face(mesh, bc, f)

         if (btype == bc_periodic) cycle

         owner = mesh%faces(f)%owner

         if (.not. flow%owned(owner)) cycle

         local_flux = local_flux + face_flux(f)
      end do
   end subroutine compute_boundary_flux


   !> Allocates and resets temporary solver vectors.
   subroutine ensure_projection_workspace(mesh)
      type(mesh_t), intent(in) :: mesh

      if (projection_work%initialized) then
         if (projection_work%ncells == mesh%ncells .and. &
             projection_work%nfaces == mesh%nfaces) return

         call finalize_flow_projection_workspace()
      end if

      projection_work%ncells = mesh%ncells
      projection_work%nfaces = mesh%nfaces

      allocate(projection_work%local_vec(3, mesh%ncells))
      allocate(projection_work%local_vec_star(3, mesh%ncells))
      allocate(projection_work%local_scalar(mesh%ncells))
      allocate(projection_work%rhs_poisson(mesh%ncells))
      allocate(projection_work%local_face_flux(mesh%nfaces))
      allocate(projection_work%predicted_face_flux(mesh%nfaces))

      allocate(projection_work%r(mesh%ncells))
      allocate(projection_work%z(mesh%ncells))
      allocate(projection_work%pvec(mesh%ncells))
      allocate(projection_work%ap(mesh%ncells))
      allocate(projection_work%local_ap(mesh%ncells))

      projection_work%local_vec = zero
      projection_work%local_vec_star = zero
      projection_work%local_scalar = zero
      projection_work%rhs_poisson = zero
      projection_work%local_face_flux = zero
      projection_work%predicted_face_flux = zero

      projection_work%r = zero
      projection_work%z = zero
      projection_work%pvec = zero
      projection_work%ap = zero
      projection_work%local_ap = zero

      projection_work%initialized = .true.
   end subroutine ensure_projection_workspace


   !> Deallocate the persistent flow projection workspace.
   subroutine finalize_flow_projection_workspace()
      if (allocated(projection_work%local_vec)) deallocate(projection_work%local_vec)
      if (allocated(projection_work%local_vec_star)) deallocate(projection_work%local_vec_star)
      if (allocated(projection_work%local_scalar)) deallocate(projection_work%local_scalar)
      if (allocated(projection_work%rhs_poisson)) deallocate(projection_work%rhs_poisson)
      if (allocated(projection_work%local_face_flux)) deallocate(projection_work%local_face_flux)
      if (allocated(projection_work%predicted_face_flux)) deallocate(projection_work%predicted_face_flux)

      if (allocated(projection_work%r)) deallocate(projection_work%r)
      if (allocated(projection_work%z)) deallocate(projection_work%z)
      if (allocated(projection_work%pvec)) deallocate(projection_work%pvec)
      if (allocated(projection_work%ap)) deallocate(projection_work%ap)
      if (allocated(projection_work%local_ap)) deallocate(projection_work%local_ap)

      projection_work%initialized = .false.
      projection_work%ncells = 0
      projection_work%nfaces = 0

      if (allocated(pressure_cache%nb)) deallocate(pressure_cache%nb)
      if (allocated(pressure_cache%coeff)) deallocate(pressure_cache%coeff)
      if (allocated(pressure_cache%diag)) deallocate(pressure_cache%diag)

      pressure_cache%initialized = .false.
      pressure_cache%has_dirichlet_pressure = .false.
      pressure_cache%has_neumann_outlet = .false.
      pressure_cache%ncells = 0
      pressure_cache%max_faces = 0

   end subroutine finalize_flow_projection_workspace


   !> Pre-computes Laplacian coefficients for the Poisson operator.
   subroutine ensure_pressure_operator_cache(mesh, flow, bc)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      integer :: ierr, c, lf, f, nb
      real(rk) :: dist, coeff
      integer :: local_dirichlet_flag, global_dirichlet_flag

      if (pressure_cache%initialized) then
         if (pressure_cache%ncells == mesh%ncells) return
         call finalize_flow_projection_workspace()
      end if

      pressure_cache%ncells = mesh%ncells
      pressure_cache%max_faces = size(mesh%cell_faces, 1)

      allocate(pressure_cache%nb(pressure_cache%max_faces, mesh%ncells))
      allocate(pressure_cache%coeff(pressure_cache%max_faces, mesh%ncells))
      allocate(pressure_cache%diag(mesh%ncells))

      pressure_cache%nb = 0
      pressure_cache%coeff = zero
      pressure_cache%diag = zero
      pressure_cache%has_neumann_outlet = .false.

      local_dirichlet_flag = 0

      do c = 1, mesh%ncells
         do lf = 1, mesh%ncell_faces(c)
            f = mesh%cell_faces(lf, c)
            nb = face_effective_neighbor(mesh, bc, f, c)

            if (nb <= 0) then
               if (patch_type_for_face(mesh, bc, f) == bc_neumann) then
                  pressure_cache%has_neumann_outlet = .true.
               end if

               if (boundary_pressure_type(mesh, bc, f) == bc_dirichlet) then
                  dist = face_normal_distance(mesh, bc, f, c, 0)
                  coeff = mesh%faces(f)%area / dist
                  pressure_cache%diag(c) = pressure_cache%diag(c) + coeff
                  local_dirichlet_flag = 1
               end if
               cycle
            end if

            dist = face_normal_distance(mesh, bc, f, c, nb)
            coeff = mesh%faces(f)%area / dist

            pressure_cache%nb(lf, c) = nb
            pressure_cache%coeff(lf, c) = coeff
            pressure_cache%diag(c) = pressure_cache%diag(c) + coeff
         end do

         if (pressure_cache%diag(c) <= tiny_safe) then
            call fatal_error('flow', 'non-positive cached pressure diagonal')
         end if
      end do

      call MPI_Allreduce(local_dirichlet_flag, global_dirichlet_flag, 1, MPI_INTEGER, MPI_MAX, flow%comm, ierr)
      pressure_cache%has_dirichlet_pressure = (global_dirichlet_flag > 0)

      pressure_cache%initialized = .true.
   end subroutine ensure_pressure_operator_cache


   !> Calculates the normal distance between cell centers or cell-to-face.
   !!
   !! Handles periodic boundary logic by accounting for the face pair offset.
   !!
   !! @param mesh The mesh.
   !! @param bc Boundary conditions.
   !! @param face_id Face index.
   !! @param cell_id Source cell index.
   !! @param nb Neighbor cell index (or 0 for boundaries).
   function face_normal_distance(mesh, bc, face_id, cell_id, nb) result(dist)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      integer, intent(in) :: nb
      real(rk) :: dist

      integer :: pair_face
      integer :: btype
      real(rk) :: nvec(3)

      nvec = outward_normal(mesh, face_id, cell_id)

      if (nb > 0) then
         if (mesh%faces(face_id)%neighbor == 0) then
            btype = patch_type_for_face(mesh, bc, face_id)

            if (btype == bc_periodic) then
               pair_face = mesh%faces(face_id)%periodic_face

               if (pair_face <= 0) then
                  call fatal_error('flow', 'periodic face has no paired face')
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
   end function face_normal_distance


   !> Computes the linear interpolation weight for a neighbor cell.
   function face_neighbor_weight(mesh, bc, face_id, cell_id, nb) result(w_nb)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      integer, intent(in) :: nb
      real(rk) :: w_nb

      integer :: pair_face
      integer :: btype
      real(rk) :: nvec(3)
      real(rk) :: d_owner
      real(rk) :: d_nb
      real(rk) :: d_total

      if (nb <= 0) then
         w_nb = zero
         return
      end if

      nvec = outward_normal(mesh, face_id, cell_id)

      d_owner = abs(dot_product(mesh%faces(face_id)%center - &
                                mesh%cells(cell_id)%center, nvec))

      if (mesh%faces(face_id)%neighbor == 0) then
         btype = patch_type_for_face(mesh, bc, face_id)

         if (btype == bc_periodic) then
            pair_face = mesh%faces(face_id)%periodic_face

            if (pair_face <= 0) then
               call fatal_error('flow', 'periodic face has no paired face')
            end if

            d_nb = abs(dot_product(mesh%cells(nb)%center - &
                                   mesh%faces(pair_face)%center, nvec))
         else
            d_nb = abs(dot_product(mesh%cells(nb)%center - &
                                   mesh%faces(face_id)%center, nvec))
         end if
      else
         d_nb = abs(dot_product(mesh%cells(nb)%center - &
                                mesh%faces(face_id)%center, nvec))
      end if

      d_total = max(d_owner + d_nb, tiny_safe)

      w_nb = d_owner / d_total
   end function face_neighbor_weight


   !> Linearly interpolates a scalar field to a face.
   function face_linear_scalar(mesh, bc, face_id, cell_id, nb, owner_value, neighbor_value) result(face_value)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      integer, intent(in) :: nb
      real(rk), intent(in) :: owner_value, neighbor_value
      real(rk) :: face_value

      real(rk) :: w_nb

      w_nb = face_neighbor_weight(mesh, bc, face_id, cell_id, nb)

      face_value = (one - w_nb) * owner_value + w_nb * neighbor_value
   end function face_linear_scalar


   !> Linearly interpolates a vector field to a face.
   function face_linear_vector(mesh, bc, face_id, cell_id, nb, owner_value, neighbor_value) result(face_value)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      integer, intent(in) :: nb
      real(rk), intent(in) :: owner_value(3), neighbor_value(3)
      real(rk) :: face_value(3)

      real(rk) :: w_nb

      w_nb = face_neighbor_weight(mesh, bc, face_id, cell_id, nb)

      face_value = (one - w_nb) * owner_value + w_nb * neighbor_value
   end function face_linear_vector


   !> Determines the outward unit normal relative to a specific cell.
   pure function outward_normal(mesh, face_id, cell_id) result(nvec)
      type(mesh_t), intent(in) :: mesh
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      real(rk) :: nvec(3)

      if (mesh%faces(face_id)%owner == cell_id) then
         nvec = mesh%faces(face_id)%normal
      else
         nvec = -mesh%faces(face_id)%normal
      end if
   end function outward_normal


   !> Internal helper for MPI error checking.
   subroutine check_mpi(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where

      if (ierr /= MPI_SUCCESS) call fatal_error('flow', trim(where)//' MPI failure')
   end subroutine check_mpi


   !> Calculates domain-wide CFL and optionally scales the timestep size.
   !!
   !! For stability in reacting flows, the timestep growth is capped 
   !! to 2% per step when `use_dynamic_dt` is active.
   subroutine compute_and_update_cfl(mesh, flow, params, fields, stats)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(inout) :: params
      type(flow_fields_t), intent(in) :: fields
      type(solver_stats_t), intent(inout) :: stats

      integer :: c, lf, f
      real(rk) :: local_cfl_rate, max_cfl_rate
      real(rk) :: cell_outward_flux
      integer :: ierr

      local_cfl_rate = zero

      do c = flow%first_cell, flow%last_cell
         cell_outward_flux = zero
         do lf = 1, mesh%ncell_faces(c)
            f = mesh%cell_faces(lf, c)
            cell_outward_flux = cell_outward_flux + abs(fields%face_flux(f))
         end do
         cell_outward_flux = half * cell_outward_flux / mesh%cells(c)%volume
         if (cell_outward_flux > local_cfl_rate) then
            local_cfl_rate = cell_outward_flux
         end if
      end do

      call MPI_Allreduce(local_cfl_rate, max_cfl_rate, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      call check_mpi(ierr, 'cfl max rate')

      if (params%use_dynamic_dt) then
         if (max_cfl_rate > tiny_safe) then
            ! Tighten dt growth cap to 1.02x per step for stability starting from rest
            params%dt = min(params%max_cfl / max_cfl_rate, params%dt * 1.02_rk)
         end if
      end if

      stats%cfl = max_cfl_rate * params%dt
   end subroutine compute_and_update_cfl

end module mod_flow_projection
