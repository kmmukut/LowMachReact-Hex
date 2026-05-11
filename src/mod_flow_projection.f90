!> Incompressible Navier-Stokes solver using the projection method.
!!
!! This module implements the fractional-step projection method for 
!! incompressible flows. It handles the momentum predictor step, 
!! pressure Poisson solve using Conjugate Gradient, and velocity correction.
module mod_flow_projection
   use mpi_f08
   use mod_kinds, only : rk, zero, one, half, tiny_safe, fatal_error
   use mod_input, only : case_params_t
   use mod_mesh_types, only : mesh_t
   use mod_mpi_flow, only : flow_mpi_t, flow_allreduce_global_vector, &
                            flow_allreduce_global_scalar, flow_allgather_owned_scalar, &
                            flow_global_dot_owned, flow_global_max_owned
   use mod_bc, only : bc_set_t, bc_periodic, bc_neumann, patch_type_for_face, &
                      boundary_velocity, face_effective_neighbor
   use mod_fields, only : flow_fields_t
   implicit none

   private

   !> Solver diagnostics and performance statistics.
   type, public :: solver_stats_t
      integer :: pressure_iterations = 0 !< Number of CG iterations performed
      real(rk) :: pressure_residual = zero !< Final pressure residual norm
      real(rk) :: max_divergence = zero    !< Maximum velocity divergence [1/s]
      real(rk) :: rms_divergence = zero    !< RMS velocity divergence [1/s]
      real(rk) :: net_boundary_flux = zero !< Net mass flux across all boundaries
      real(rk) :: kinetic_energy = zero    !< Total kinetic energy in domain [J]
   end type solver_stats_t

   type :: pressure_operator_cache_t
      logical :: initialized = .false.
      integer :: ncells = 0
      integer :: max_faces = 0
      integer, allocatable :: nb(:,:)       ! (max_faces,ncells)
      real(rk), allocatable :: coeff(:,:)   ! (max_faces,ncells)
      real(rk), allocatable :: diag(:)      ! (ncells)
   end type pressure_operator_cache_t

   type :: projection_workspace_t
      logical :: initialized = .false.
      integer :: ncells = 0
      integer :: nfaces = 0

      real(rk), allocatable :: local_vec(:,:)
      real(rk), allocatable :: momentum_rhs(:,:)
      real(rk), allocatable :: local_scalar(:)
      real(rk), allocatable :: rhs_poisson(:)
      real(rk), allocatable :: local_face_flux(:)
      real(rk), allocatable :: predicted_face_flux(:)

      real(rk), allocatable :: r(:)
      real(rk), allocatable :: z(:)
      real(rk), allocatable :: pvec(:)
      real(rk), allocatable :: ap(:)
      real(rk), allocatable :: local_ap(:)
   end type projection_workspace_t

   type(pressure_operator_cache_t), save :: pressure_cache
   type(projection_workspace_t), save :: projection_work

   public :: advance_projection_step, compute_flow_diagnostics
   public :: finalize_flow_projection_workspace
   !> Compute distance from cell center to face center along normal.
   !!
   !! @param mesh Full mesh.
   !! @param bc Boundary conditions.
   !! @param face_id Face ID.
   !! @param cell_id ID of the cell from which distance is measured.
   !! @param neighbor_id ID of the neighbor cell (used for periodic handling).
   function face_normal_distance(mesh, bc, face_id, cell_id, neighbor_id) result(dist)

contains

   !> Advance the flow solution by one projection step.
   !!
   !! 1. Compute explicit momentum RHS.
   !! 2. Predict intermediate velocity (u*).
   !! 3. Solve pressure Poisson equation for potential (phi).
   !! 4. Correct velocity and update pressure.
   !!
   !! @param mesh Full replicated mesh.
   !! @param flow MPI flow decomposition data.
   !! @param bc Boundary condition set.
   !! @param params Case configuration parameters.
   !! @param fields Flow fields to update.
   !! @param stats Solver statistics to populate.
   subroutine advance_projection_step(mesh, flow, bc, params, fields, stats)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(inout) :: fields
      type(solver_stats_t), intent(out) :: stats

      integer :: c

      call ensure_projection_workspace(mesh)
      call ensure_pressure_operator_cache(mesh, flow, bc)

      associate( &
         local_vec => projection_work%local_vec, &
         momentum_rhs => projection_work%momentum_rhs, &
         local_scalar => projection_work%local_scalar, &
         rhs_poisson => projection_work%rhs_poisson, &
         local_face_flux => projection_work%local_face_flux, &
         predicted_face_flux => projection_work%predicted_face_flux)

         fields%u_old = fields%u

         ! ------------------------------------------------------------
         ! 1. Explicit momentum RHS using old pressure gradient.
         ! ------------------------------------------------------------
         call compute_momentum_rhs(mesh, flow, bc, params, fields%u, fields%p, local_vec)
         call flow_allreduce_global_vector(flow, local_vec, momentum_rhs)

         ! ------------------------------------------------------------
         ! 2. AB2 update for predicted cell velocity.
         !    First step falls back to forward Euler.
         ! ------------------------------------------------------------
         call advance_ab2(mesh, flow, params, fields, momentum_rhs, local_vec)
         call flow_allreduce_global_vector(flow, local_vec, fields%u_star)

         fields%rhs_old = momentum_rhs
         fields%rhs_old_valid = .true.

         ! ------------------------------------------------------------
         ! 3. Predicted conservative face flux.
         ! ------------------------------------------------------------
         call compute_predicted_face_flux(mesh, flow, bc, fields%u_star, local_face_flux)
         call flow_allreduce_global_scalar(flow, local_face_flux, predicted_face_flux)

         ! ------------------------------------------------------------
         ! 3b. Open-boundary mass balance.
         ! ------------------------------------------------------------
         call balance_neumann_outlet_flux(mesh, flow, bc, predicted_face_flux)

         ! ------------------------------------------------------------
         ! 4. Conservative pressure-correction RHS.
         !
         ! fields%div is flux_sum / cell_volume. For the symmetric
         ! unnormalized FV Poisson operator, use flux_sum.
         ! ------------------------------------------------------------
         call compute_flux_divergence(mesh, flow, predicted_face_flux, local_scalar)
         call flow_allreduce_global_scalar(flow, local_scalar, fields%div)

         rhs_poisson = zero
         do c = 1, mesh%ncells
            rhs_poisson(c) = -params%rho / params%dt * &
                             fields%div(c) * mesh%cells(c)%volume
         end do

         rhs_poisson(1) = zero
         fields%phi = zero

         call solve_pressure_correction(mesh, flow, bc, params, rhs_poisson, fields%phi, stats)

         ! ------------------------------------------------------------
         ! 5. Incremental pressure update.
         ! ------------------------------------------------------------
         fields%p = fields%p + fields%phi
         fields%p(1) = zero

         ! ------------------------------------------------------------
         ! 6. Correct conservative face flux.
         ! ------------------------------------------------------------
         call correct_face_flux(mesh, flow, bc, params, predicted_face_flux, fields%phi, local_face_flux)
         call flow_allreduce_global_scalar(flow, local_face_flux, fields%face_flux)

         ! ------------------------------------------------------------
         ! 7. Correct cell-centered velocity.
         ! ------------------------------------------------------------
         call correct_cell_velocity(mesh, flow, bc, params, fields%u_star, fields%phi, local_vec)
         call flow_allreduce_global_vector(flow, local_vec, fields%u)

         ! ------------------------------------------------------------
         ! 8. Final divergence diagnostic from corrected flux.
         ! ------------------------------------------------------------
         call compute_flux_divergence(mesh, flow, fields%face_flux, local_scalar)
         call flow_allreduce_global_scalar(flow, local_scalar, fields%div)

         call compute_flow_diagnostics(mesh, flow, bc, params, fields, stats)

      end associate
   end subroutine advance_projection_step


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


   subroutine compute_momentum_rhs(mesh, flow, bc, params, u, p, local_rhs)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: u(:,:)
      real(rk), intent(in) :: p(:)
      real(rk), intent(out) :: local_rhs(:,:)

      integer :: c, lf, f, nb
      real(rk) :: nvec(3), nhat(3)
      real(rk) :: uf(3), ub(3), advected(3)
      real(rk) :: un_area
      real(rk) :: conv(3), diff(3), gradp(3)
      real(rk) :: dist
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

            diff = diff + params%nu * mesh%faces(f)%area * (ub - u(:, c)) / dist
         end do

         local_rhs(:, c) = -conv / mesh%cells(c)%volume + &
                            diff / mesh%cells(c)%volume + &
                            params%body_force - gradp / params%rho
      end do
   end subroutine compute_momentum_rhs


   subroutine pressure_gradient_cell(mesh, bc, p, cell_id, gradp)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(in) :: p(:)
      integer, intent(in) :: cell_id
      real(rk), intent(out) :: gradp(3)

      integer :: lf, f, nb
      real(rk) :: nvec(3)
      real(rk) :: pf

      gradp = zero

      do lf = 1, mesh%ncell_faces(cell_id)
         f = mesh%cell_faces(lf, cell_id)
         nvec = outward_normal(mesh, f, cell_id)

         nb = face_effective_neighbor(mesh, bc, f, cell_id)

         if (nb > 0) then
            pf = face_linear_scalar(mesh, bc, f, cell_id, nb, p(cell_id), p(nb))
         else
            pf = p(cell_id)
         end if

         gradp = gradp + pf * nvec * mesh%faces(f)%area
      end do

      gradp = gradp / mesh%cells(cell_id)%volume
   end subroutine pressure_gradient_cell


   subroutine compute_predicted_face_flux(mesh, flow, bc, ustar, local_flux)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(in) :: ustar(:,:)
      real(rk), intent(out) :: local_flux(:)

      integer :: f, owner, nb
      real(rk) :: uf(3), ub(3)

      local_flux = zero

      do f = 1, mesh%nfaces
         owner = mesh%faces(f)%owner

         if (.not. flow%owned(owner)) cycle

         nb = face_effective_neighbor(mesh, bc, f, owner)

         if (nb > 0) then
            uf = face_linear_vector(mesh, bc, f, owner, nb, ustar(:, owner), ustar(:, nb))
         else
            call boundary_velocity(mesh, bc, f, ustar(:, owner), ub)
            uf = ub
         end if

         local_flux(f) = dot_product(uf, mesh%faces(f)%normal) * mesh%faces(f)%area
      end do
   end subroutine compute_predicted_face_flux


   subroutine correct_face_flux(mesh, flow, bc, params, predicted_flux, phi, local_flux)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      real(rk), intent(in) :: predicted_flux(:)
      real(rk), intent(in) :: phi(:)
      real(rk), intent(out) :: local_flux(:)

      integer :: f, owner, nb
      real(rk) :: dist

      local_flux = zero

      do f = 1, mesh%nfaces
         owner = mesh%faces(f)%owner

         if (.not. flow%owned(owner)) cycle

         nb = face_effective_neighbor(mesh, bc, f, owner)

         local_flux(f) = predicted_flux(f)

         if (nb > 0) then
            dist = face_normal_distance(mesh, bc, f, owner, nb)

            local_flux(f) = predicted_flux(f) - params%dt / params%rho * &
                            mesh%faces(f)%area * (phi(nb) - phi(owner)) / dist
         end if
      end do
   end subroutine correct_face_flux


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

         call pressure_matvec(mesh, flow, bc, phi, local_ap)
         call flow_allgather_owned_scalar(flow, local_ap, ap)

         r = rhs - ap

         r(1) = zero
         phi(1) = zero

         do c = 1, mesh%ncells
            if (c == 1) then
               z(c) = zero
            else
               z(c) = r(c) / max(pressure_cache%diag(c), tiny_safe)
            end if
         end do

         pvec = z
         pvec(1) = zero

         rr = flow_global_dot_owned(flow, r, r)
         rz = flow_global_dot_owned(flow, r, z)

         stats%pressure_iterations = 0
         stats%pressure_residual = sqrt(rr / max(real(mesh%ncells, rk), one))

         do iter = 1, params%pressure_max_iter
            if (stats%pressure_residual <= params%pressure_tol) exit

            call pressure_matvec(mesh, flow, bc, pvec, local_ap)
            call flow_allgather_owned_scalar(flow, local_ap, ap)

            pap = flow_global_dot_owned(flow, pvec, ap)

            if (abs(pap) <= tiny_safe) exit

            alpha = rz / pap

            phi = phi + alpha * pvec
            r = r - alpha * ap

            phi(1) = zero
            r(1) = zero

            rr_new = flow_global_dot_owned(flow, r, r)

            stats%pressure_residual = sqrt(rr_new / max(real(mesh%ncells, rk), one))
            stats%pressure_iterations = iter

            if (stats%pressure_residual <= params%pressure_tol) exit

            do c = 1, mesh%ncells
               if (c == 1) then
                  z(c) = zero
               else
                  z(c) = r(c) / max(pressure_cache%diag(c), tiny_safe)
               end if
            end do

            rz_new = flow_global_dot_owned(flow, r, z)

            beta = rz_new / max(rz, tiny_safe)

            pvec = z + beta * pvec
            pvec(1) = zero

            rr = rr_new
            rz = rz_new
         end do

         phi(1) = zero

      end associate
   end subroutine solve_pressure_correction


   subroutine pressure_matvec(mesh, flow, bc, x, local_ax)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc
      real(rk), intent(in) :: x(:)
      real(rk), intent(out) :: local_ax(:)

      integer :: c, lf, nb

      call ensure_pressure_operator_cache(mesh, flow, bc)

      local_ax = zero

      do c = flow%first_cell, flow%last_cell
         if (c == 1) then
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
   end subroutine pressure_matvec


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
               phi_face = phi(c)
            end if

            grad = grad + phi_face * nvec * mesh%faces(f)%area
         end do

         grad = grad / mesh%cells(c)%volume

         local_u(:, c) = ustar(:, c) - params%dt / params%rho * grad
      end do
   end subroutine correct_cell_velocity


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

      call MPI_Allreduce(local_net_flux, global_net_flux, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'outlet balance net flux')

      call MPI_Allreduce(local_outlet_area, global_outlet_area, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'outlet balance outlet area')

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
      allocate(projection_work%momentum_rhs(3, mesh%ncells))
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
      projection_work%momentum_rhs = zero
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
      if (allocated(projection_work%momentum_rhs)) deallocate(projection_work%momentum_rhs)
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
      pressure_cache%ncells = 0
      pressure_cache%max_faces = 0
   end subroutine finalize_flow_projection_workspace


   subroutine ensure_pressure_operator_cache(mesh, flow, bc)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(bc_set_t), intent(in) :: bc

      integer :: c, lf, f, nb
      real(rk) :: dist
      real(rk) :: coeff

      associate(dummy => flow%nlocal)
      end associate

      if (pressure_cache%initialized) then
         if (pressure_cache%ncells == mesh%ncells) return

         if (allocated(pressure_cache%nb)) deallocate(pressure_cache%nb)
         if (allocated(pressure_cache%coeff)) deallocate(pressure_cache%coeff)
         if (allocated(pressure_cache%diag)) deallocate(pressure_cache%diag)

         pressure_cache%initialized = .false.
      end if

      pressure_cache%ncells = mesh%ncells
      pressure_cache%max_faces = size(mesh%cell_faces, 1)

      allocate(pressure_cache%nb(pressure_cache%max_faces, mesh%ncells))
      allocate(pressure_cache%coeff(pressure_cache%max_faces, mesh%ncells))
      allocate(pressure_cache%diag(mesh%ncells))

      pressure_cache%nb = 0
      pressure_cache%coeff = zero
      pressure_cache%diag = zero

      do c = 1, mesh%ncells
         if (c == 1) then
            pressure_cache%diag(c) = one
            cycle
         end if

         do lf = 1, mesh%ncell_faces(c)
            f = mesh%cell_faces(lf, c)
            nb = face_effective_neighbor(mesh, bc, f, c)

            if (nb <= 0) cycle

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

      pressure_cache%initialized = .true.
   end subroutine ensure_pressure_operator_cache


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


   function face_linear_scalar(mesh, bc, face_id, cell_id, nb, owner_value, neighbor_value) result(face_value)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      integer, intent(in) :: nb
      real(rk), intent(in) :: owner_value
      real(rk), intent(in) :: neighbor_value
      real(rk) :: face_value

      real(rk) :: w_nb

      w_nb = face_neighbor_weight(mesh, bc, face_id, cell_id, nb)

      face_value = (one - w_nb) * owner_value + w_nb * neighbor_value
   end function face_linear_scalar


   function face_linear_vector(mesh, bc, face_id, cell_id, nb, owner_value, neighbor_value) result(face_value)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id
      integer, intent(in) :: nb
      real(rk), intent(in) :: owner_value(3)
      real(rk), intent(in) :: neighbor_value(3)
      real(rk) :: face_value(3)

      real(rk) :: w_nb

      w_nb = face_neighbor_weight(mesh, bc, face_id, cell_id, nb)

      face_value = (one - w_nb) * owner_value + w_nb * neighbor_value
   end function face_linear_vector


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


   subroutine check_mpi(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where

      if (ierr /= MPI_SUCCESS) call fatal_error('flow', trim(where)//' MPI failure')
   end subroutine check_mpi

end module mod_flow_projection