!> MPI parallelization and domain decomposition for the flow solver.
!!
!! This module manages the parallel execution of the hydrodynamic solver 
!! using a **Replicated Mesh** strategy. In this approach, every MPI rank 
!! holds the full topological structure of the mesh, but performs physical 
!! updates (fluxes, source terms) only on a specific subset of "owned" cells.
!!
!! Collective communication (Allreduce, Allgather) is then used to synchronize 
!! the full field across all processors. This simplifies the implementation of 
!! complex unstructured operators while allowing multi-CPU acceleration.
module mod_mpi_flow
   use mpi_f08
   use mod_kinds, only : rk, zero, fatal_error
   use mod_mesh_types, only : mesh_t
   implicit none

   private

   public :: flow_mpi_t
   public :: mpi_flow_startup, mpi_flow_shutdown
   public :: flow_mpi_initialize, flow_mpi_finalize
   public :: flow_allreduce_global_vector, flow_allreduce_global_scalar
   public :: flow_global_dot_owned, flow_global_dots_owned, flow_global_sum_owned, flow_global_max_owned
   public :: flow_global_two_dots_owned
   public :: flow_allgather_owned_scalar, flow_allgather_owned_vector
   public :: flow_allgather_owned_matrix, flow_allgather_owned_matrix_inplace
   public :: flow_allgather_owned_v4
   public :: flow_exchange_cell_scalar, flow_exchange_cell_matrix
   public :: flow_exchange_face_scalar
   public :: flow_gather_owned_scalar_root, flow_gather_owned_matrix_root

   !> MPI context for hydrodynamic operations.
   !!
   !! Stores rank information and pre-calculated cell ranges to avoid 
   !! repeated division logic during the simulation loop.
   type :: flow_mpi_t
      type(MPI_Comm) :: comm = MPI_COMM_NULL !< MPI Communicator for flow.
      integer :: rank = -1                   !< Local rank ID (0 to nprocs-1).
      integer :: nprocs = 0                  !< Total number of flow processors.
      integer :: first_cell = 0              !< First cell index owned by this rank.
      integer :: last_cell = -1              !< Last cell index owned by this rank.
      integer :: nlocal = 0                  !< Total number of cells owned locally.
      logical, allocatable :: owned(:)       !< Bitmask for all cells (True if owned by this rank).
      integer, allocatable :: cell_owner(:)  !< Zero-based MPI owner rank for each cell.
      integer, allocatable :: ghost_cells(:) !< Off-rank neighbor cells needed by owned cells.
      integer, allocatable :: owned_faces(:) !< Faces whose mesh owner cell is owned locally.

      ! Cached metadata/buffers for gathering owned cell ranges.
      integer, allocatable :: gather_counts(:)
      integer, allocatable :: gather_displs(:)
      integer, allocatable :: gather_firsts(:)
      real(rk), allocatable :: gather_sendbuf(:)
      real(rk), allocatable :: gather_recvbuf(:)
      integer :: gather_max_components = 0
      integer, allocatable :: gather_matrix_counts(:)
      integer, allocatable :: gather_matrix_displs(:)
      real(rk), allocatable :: gather_matrix_sendbuf(:)
      real(rk), allocatable :: gather_matrix_recvbuf(:)

      ! Cached metadata/buffers for owned-cell halo exchange.
      integer :: cell_halo_max_components = 4
      integer :: ncell_send_ranks = 0
      integer :: ncell_recv_ranks = 0
      integer, allocatable :: cell_send_ranks(:)
      integer, allocatable :: cell_recv_ranks(:)
      integer, allocatable :: cell_send_counts(:)
      integer, allocatable :: cell_recv_counts(:)
      integer, allocatable :: cell_send_displs(:)
      integer, allocatable :: cell_recv_displs(:)
      integer, allocatable :: cell_send_cells(:)
      integer, allocatable :: cell_recv_cells(:)
      real(rk), allocatable :: cell_sendbuf(:)
      real(rk), allocatable :: cell_recvbuf(:)
      type(MPI_Request), allocatable :: cell_requests(:)

      ! Cached metadata/buffers for face-flux exchange.
      integer :: nface_send_ranks = 0
      integer :: nface_recv_ranks = 0
      integer, allocatable :: face_send_ranks(:)
      integer, allocatable :: face_recv_ranks(:)
      integer, allocatable :: face_send_counts(:)
      integer, allocatable :: face_recv_counts(:)
      integer, allocatable :: face_send_displs(:)
      integer, allocatable :: face_recv_displs(:)
      integer, allocatable :: face_send_faces(:)
      integer, allocatable :: face_recv_faces(:)
      real(rk), allocatable :: face_sendbuf(:)
      real(rk), allocatable :: face_recvbuf(:)
      type(MPI_Request), allocatable :: face_requests(:)
   end type flow_mpi_t

   logical :: mpi_started_here = .false.

contains

   !> Initializes the MPI environment if not already active.
   subroutine mpi_flow_startup()
      logical :: initialized
      integer :: ierr

      call MPI_Initialized(initialized, ierr)
      call check_mpi(ierr, 'MPI_Initialized')

      if (.not. initialized) then
         call MPI_Init(ierr)
         call check_mpi(ierr, 'MPI_Init')
         mpi_started_here = .true.
      end if
   end subroutine mpi_flow_startup


   !> Shuts down the MPI environment if it was started by this module.
   subroutine mpi_flow_shutdown()
      logical :: finalized
      integer :: ierr

      call MPI_Finalized(finalized, ierr)
      if (ierr /= MPI_SUCCESS) return
      if (finalized) return

      if (mpi_started_here) then
         call MPI_Finalize(ierr)
      end if
   end subroutine mpi_flow_shutdown


   !> Sets up domain decomposition for a given mesh.
   !!
   !! Splits the total number of cells among available processors 
   !! using a contiguous block decomposition.
   !!
   !! @param mesh The mesh to decompose.
   !! @param flow The MPI context to populate.
   !! @param comm_parent The parent communicator (usually MPI_COMM_WORLD).
   subroutine flow_mpi_initialize(mesh, flow, comm_parent, max_gather_components)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(MPI_Comm), intent(in) :: comm_parent
      integer, intent(in), optional :: max_gather_components

      integer :: ierr, base, rem, c
      integer :: max_components

      call flow_mpi_finalize(flow)

      call MPI_Comm_dup(comm_parent, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Comm_dup flow')

      call MPI_Comm_rank(flow%comm, flow%rank, ierr)
      call check_mpi(ierr, 'MPI_Comm_rank flow')

      call MPI_Comm_size(flow%comm, flow%nprocs, ierr)
      call check_mpi(ierr, 'MPI_Comm_size flow')

      ! Contiguous cell range calculation:
      ! The total ncells are split into roughly equal blocks. If ncells is not 
      ! divisible by nprocs, the first 'rem' ranks get an extra cell to ensure 
      ! all cells are covered.
      base = mesh%ncells / flow%nprocs
      rem = mod(mesh%ncells, flow%nprocs)

      if (flow%rank < rem) then
         flow%nlocal = base + 1
         flow%first_cell = flow%rank * (base + 1) + 1
      else
         flow%nlocal = base
         flow%first_cell = rem * (base + 1) + (flow%rank - rem) * base + 1
      end if

      ! The global cell indices owned by this rank are [first_cell, last_cell].
      flow%last_cell = flow%first_cell + flow%nlocal - 1

      allocate(flow%owned(mesh%ncells))
      flow%owned = .false.

      do c = flow%first_cell, flow%last_cell
         if (c >= 1 .and. c <= mesh%ncells) flow%owned(c) = .true.
      end do

      max_components = 4
      if (present(max_gather_components)) max_components = max(max_components, max_gather_components)
      flow%cell_halo_max_components = max(1, max_components)

      call setup_owned_gather(mesh, flow, max_gather_components)
      call setup_cell_owners(mesh, flow)
      call setup_owned_faces(mesh, flow)
      call setup_cell_halo(mesh, flow)
      call setup_face_halo(mesh, flow)
   end subroutine flow_mpi_initialize


   !> Releases all MPI resources and buffers.
   subroutine flow_mpi_finalize(flow)
      type(flow_mpi_t), intent(inout) :: flow
      integer :: ierr

      if (allocated(flow%owned)) deallocate(flow%owned)
      if (allocated(flow%gather_counts)) deallocate(flow%gather_counts)
      if (allocated(flow%gather_displs)) deallocate(flow%gather_displs)
      if (allocated(flow%gather_firsts)) deallocate(flow%gather_firsts)
      if (allocated(flow%gather_sendbuf)) deallocate(flow%gather_sendbuf)
      if (allocated(flow%gather_recvbuf)) deallocate(flow%gather_recvbuf)
      if (allocated(flow%gather_matrix_counts)) deallocate(flow%gather_matrix_counts)
      if (allocated(flow%gather_matrix_displs)) deallocate(flow%gather_matrix_displs)
      if (allocated(flow%gather_matrix_sendbuf)) deallocate(flow%gather_matrix_sendbuf)
      if (allocated(flow%gather_matrix_recvbuf)) deallocate(flow%gather_matrix_recvbuf)
      if (allocated(flow%cell_owner)) deallocate(flow%cell_owner)
      if (allocated(flow%ghost_cells)) deallocate(flow%ghost_cells)
      if (allocated(flow%owned_faces)) deallocate(flow%owned_faces)
      if (allocated(flow%cell_send_ranks)) deallocate(flow%cell_send_ranks)
      if (allocated(flow%cell_recv_ranks)) deallocate(flow%cell_recv_ranks)
      if (allocated(flow%cell_send_counts)) deallocate(flow%cell_send_counts)
      if (allocated(flow%cell_recv_counts)) deallocate(flow%cell_recv_counts)
      if (allocated(flow%cell_send_displs)) deallocate(flow%cell_send_displs)
      if (allocated(flow%cell_recv_displs)) deallocate(flow%cell_recv_displs)
      if (allocated(flow%cell_send_cells)) deallocate(flow%cell_send_cells)
      if (allocated(flow%cell_recv_cells)) deallocate(flow%cell_recv_cells)
      if (allocated(flow%cell_sendbuf)) deallocate(flow%cell_sendbuf)
      if (allocated(flow%cell_recvbuf)) deallocate(flow%cell_recvbuf)
      if (allocated(flow%cell_requests)) deallocate(flow%cell_requests)
      if (allocated(flow%face_send_ranks)) deallocate(flow%face_send_ranks)
      if (allocated(flow%face_recv_ranks)) deallocate(flow%face_recv_ranks)
      if (allocated(flow%face_send_counts)) deallocate(flow%face_send_counts)
      if (allocated(flow%face_recv_counts)) deallocate(flow%face_recv_counts)
      if (allocated(flow%face_send_displs)) deallocate(flow%face_send_displs)
      if (allocated(flow%face_recv_displs)) deallocate(flow%face_recv_displs)
      if (allocated(flow%face_send_faces)) deallocate(flow%face_send_faces)
      if (allocated(flow%face_recv_faces)) deallocate(flow%face_recv_faces)
      if (allocated(flow%face_sendbuf)) deallocate(flow%face_sendbuf)
      if (allocated(flow%face_recvbuf)) deallocate(flow%face_recvbuf)
      if (allocated(flow%face_requests)) deallocate(flow%face_requests)

      if (flow%comm /= MPI_COMM_NULL) then
         call MPI_Comm_free(flow%comm, ierr)
         call check_mpi(ierr, 'MPI_Comm_free flow')
      end if

      flow%comm = MPI_COMM_NULL
      flow%rank = -1
      flow%nprocs = 0
      flow%first_cell = 0
      flow%last_cell = -1
      flow%nlocal = 0
      flow%gather_max_components = 0
      flow%ncell_send_ranks = 0
      flow%ncell_recv_ranks = 0
      flow%nface_send_ranks = 0
      flow%nface_recv_ranks = 0
   end subroutine flow_mpi_finalize


   !> Pre-calculates MPI gather offsets and counts for allgather operations.
   subroutine setup_owned_gather(mesh, flow, max_gather_components)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      integer, intent(in), optional :: max_gather_components
      integer :: ierr, r, total_count
      integer :: max_components

      allocate(flow%gather_counts(flow%nprocs))
      allocate(flow%gather_displs(flow%nprocs))
      allocate(flow%gather_firsts(flow%nprocs))

      call MPI_Allgather(flow%nlocal, 1, MPI_INTEGER, &
                         flow%gather_counts, 1, MPI_INTEGER, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allgather gather counts')

      call MPI_Allgather(flow%first_cell, 1, MPI_INTEGER, &
                         flow%gather_firsts, 1, MPI_INTEGER, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allgather gather firsts')

      flow%gather_displs(1) = 0
      do r = 2, flow%nprocs
         flow%gather_displs(r) = flow%gather_displs(r - 1) + flow%gather_counts(r - 1)
      end do

      total_count = sum(flow%gather_counts)

      if (total_count /= mesh%ncells) then
         call fatal_error('mpi_flow', 'owned gather counts do not sum to ncells')
      end if

      allocate(flow%gather_sendbuf(flow%nlocal))
      allocate(flow%gather_recvbuf(total_count))

      max_components = 4
      if (present(max_gather_components)) max_components = max(max_components, max_gather_components)
      max_components = max(1, max_components)

      flow%gather_max_components = max_components
      allocate(flow%gather_matrix_counts(flow%nprocs))
      allocate(flow%gather_matrix_displs(flow%nprocs))
      allocate(flow%gather_matrix_sendbuf(flow%nlocal * max_components))
      allocate(flow%gather_matrix_recvbuf(total_count * max_components))
   end subroutine setup_owned_gather


   !> Sum-Allreduce for a 3D global vector field.
   subroutine flow_allreduce_global_vector(flow, local_values, global_values)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: local_values(:,:)
      real(rk), intent(out) :: global_values(:,:)
      integer :: ierr

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_values, global_values, size(local_values), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow vector')
      call profiler_stop('MPI_Communication')
   end subroutine flow_allreduce_global_vector


   !> Sum-Allreduce for a global scalar field.
   subroutine flow_allreduce_global_scalar(flow, local_values, global_values)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: local_values(:)
      real(rk), intent(out) :: global_values(:)
      integer :: ierr

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_values, global_values, size(local_values), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow scalar')
      call profiler_stop('MPI_Communication')
   end subroutine flow_allreduce_global_scalar


   !> Gathers locally-updated cell values and broadcasts to the global mesh.
   !!
   !! This routine uses `MPI_Allgatherv` to synchronize the "owned" 
   !! partition results across all ranks.
   subroutine flow_allgather_owned_scalar(flow, local_global, global)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(in) :: local_global(:)
      real(rk), intent(out) :: global(:)
      integer :: ierr, r, first

      if (.not. allocated(flow%gather_sendbuf)) then
         call fatal_error('mpi_flow', 'owned gather buffers are not initialized')
      end if

      flow%gather_sendbuf = local_global(flow%first_cell:flow%last_cell)

      call profiler_start('MPI_Communication')
      call MPI_Allgatherv(flow%gather_sendbuf, flow%nlocal, MPI_DOUBLE_PRECISION, &
                          flow%gather_recvbuf, flow%gather_counts, flow%gather_displs, &
                          MPI_DOUBLE_PRECISION, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allgatherv owned scalar')
      call profiler_stop('MPI_Communication')

      global = zero

      do r = 1, flow%nprocs
         first = flow%gather_firsts(r)

         global(first:first + flow%gather_counts(r) - 1) = &
            flow%gather_recvbuf(flow%gather_displs(r) + 1: &
                                flow%gather_displs(r) + flow%gather_counts(r))
      end do
   end subroutine flow_allgather_owned_scalar


   !> Gathers locally-updated 3D vector cell values and broadcasts to the global mesh.
   !!
   !! This optimized version performs a single MPI_Allgatherv for all 3 components
   !! to minimize communication latency.
   subroutine flow_allgather_owned_vector(flow, local_global, global)
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(in) :: local_global(:,:)
      real(rk), intent(out) :: global(:,:)

      call flow_allgather_owned_matrix(flow, local_global, global)
   end subroutine flow_allgather_owned_vector


   !> Gathers locally-updated matrix cell values and broadcasts to the global mesh.
   subroutine flow_allgather_owned_matrix(flow, local_global, global)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(in) :: local_global(:,:)
      real(rk), intent(out) :: global(:,:)
      integer :: ierr, ncomp, ncells, nlocal_comp
      integer :: c, k, r, first, pos, recv_pos

      ncomp = size(global, 1)
      ncells = size(global, 2)

      if (size(local_global, 1) /= ncomp .or. size(local_global, 2) /= ncells) then
         call fatal_error('mpi_flow', 'owned matrix gather shape mismatch')
      end if

      call prepare_matrix_gather(flow, ncomp, ncells, nlocal_comp)

      pos = 0
      do c = flow%first_cell, flow%last_cell
         do k = 1, ncomp
            pos = pos + 1
            flow%gather_matrix_sendbuf(pos) = local_global(k, c)
         end do
      end do

      call profiler_start('MPI_Communication')
      call MPI_Allgatherv(flow%gather_matrix_sendbuf, nlocal_comp, MPI_DOUBLE_PRECISION, &
                          flow%gather_matrix_recvbuf, flow%gather_matrix_counts, &
                          flow%gather_matrix_displs, MPI_DOUBLE_PRECISION, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allgatherv owned matrix')
      call profiler_stop('MPI_Communication')

      global = zero

      do r = 1, flow%nprocs
         first = flow%gather_firsts(r)
         recv_pos = flow%gather_matrix_displs(r)

         do c = first, first + flow%gather_counts(r) - 1
            do k = 1, ncomp
               recv_pos = recv_pos + 1
               global(k, c) = flow%gather_matrix_recvbuf(recv_pos)
            end do
         end do
      end do
   end subroutine flow_allgather_owned_matrix


   !> In-place variant for fields that already hold owned-cell updates.
   subroutine flow_allgather_owned_matrix_inplace(flow, field)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(inout) :: field(:,:)
      integer :: ierr, ncomp, ncells, nlocal_comp
      integer :: c, k, r, first, pos, recv_pos

      ncomp = size(field, 1)
      ncells = size(field, 2)

      call prepare_matrix_gather(flow, ncomp, ncells, nlocal_comp)

      pos = 0
      do c = flow%first_cell, flow%last_cell
         do k = 1, ncomp
            pos = pos + 1
            flow%gather_matrix_sendbuf(pos) = field(k, c)
         end do
      end do

      call profiler_start('MPI_Communication')
      call MPI_Allgatherv(flow%gather_matrix_sendbuf, nlocal_comp, MPI_DOUBLE_PRECISION, &
                          flow%gather_matrix_recvbuf, flow%gather_matrix_counts, &
                          flow%gather_matrix_displs, MPI_DOUBLE_PRECISION, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allgatherv owned matrix inplace')
      call profiler_stop('MPI_Communication')

      field = zero

      do r = 1, flow%nprocs
         first = flow%gather_firsts(r)
         recv_pos = flow%gather_matrix_displs(r)

         do c = first, first + flow%gather_counts(r) - 1
            do k = 1, ncomp
               recv_pos = recv_pos + 1
               field(k, c) = flow%gather_matrix_recvbuf(recv_pos)
            end do
         end do
      end do
   end subroutine flow_allgather_owned_matrix_inplace


   !> Gathers 4-component cell values (e.g., Velocity + Scalar) in one call.
   subroutine flow_allgather_owned_v4(flow, local_v, local_s, global_v, global_s)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(in) :: local_v(:,:), local_s(:)
      real(rk), intent(out) :: global_v(:,:), global_s(:)
      integer :: ierr, nlocal4, ncells
      integer :: c, i, r, first, recv_pos

      ncells = size(global_s)
      if (size(local_s) /= ncells .or. size(global_v, 2) /= ncells .or. &
          size(local_v, 2) /= ncells .or. size(local_v, 1) /= 3 .or. &
          size(global_v, 1) /= 3) then
         call fatal_error('mpi_flow', 'owned v4 gather shape mismatch')
      end if

      call prepare_matrix_gather(flow, 4, ncells, nlocal4)

      ! Pack: (U, V, W, S) for owned cells
      i = 0
      do c = flow%first_cell, flow%last_cell
         flow%gather_matrix_sendbuf(i + 1:i + 3) = local_v(:, c)
         flow%gather_matrix_sendbuf(i + 4) = local_s(c)
         i = i + 4
      end do

      call profiler_start('MPI_Communication')
      call MPI_Allgatherv(flow%gather_matrix_sendbuf, nlocal4, MPI_DOUBLE_PRECISION, &
                          flow%gather_matrix_recvbuf, flow%gather_matrix_counts, &
                          flow%gather_matrix_displs, MPI_DOUBLE_PRECISION, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allgatherv owned v4')
      call profiler_stop('MPI_Communication')

      global_v = zero
      global_s = zero

      do r = 1, flow%nprocs
         first = flow%gather_firsts(r)
         recv_pos = flow%gather_matrix_displs(r)

         do c = first, first + flow%gather_counts(r) - 1
            global_v(:, c) = flow%gather_matrix_recvbuf(recv_pos + 1:recv_pos + 3)
            global_s(c) = flow%gather_matrix_recvbuf(recv_pos + 4)
            recv_pos = recv_pos + 4
         end do
      end do
   end subroutine flow_allgather_owned_v4


   !> Updates cached count/displacement arrays for a packed component count.
   subroutine prepare_matrix_gather(flow, ncomp, ncells, nlocal_comp)
      type(flow_mpi_t), intent(inout) :: flow
      integer, intent(in) :: ncomp, ncells
      integer, intent(out) :: nlocal_comp

      if (.not. allocated(flow%gather_matrix_sendbuf)) then
         call fatal_error('mpi_flow', 'owned matrix gather buffers are not initialized')
      end if

      if (ncomp > flow%gather_max_components) then
         call fatal_error('mpi_flow', 'owned matrix gather component count exceeds cached buffer size')
      end if

      if (sum(flow%gather_counts) /= ncells) then
         call fatal_error('mpi_flow', 'owned matrix gather cell count mismatch')
      end if

      nlocal_comp = flow%nlocal * ncomp
      flow%gather_matrix_counts = flow%gather_counts * ncomp
      flow%gather_matrix_displs = flow%gather_displs * ncomp
   end subroutine prepare_matrix_gather


   !> Computes multiple global dot products in a single MPI_Allreduce.
   !!
   !! This batched version reduces MPI latency by combining n_dots synchronizations.
   subroutine flow_global_dots_owned(flow, n_dots, a, b, results)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(in) :: flow
      integer, intent(in) :: n_dots
      real(rk), intent(in) :: a(:,:)  ! (ncells, n_dots)
      real(rk), intent(in) :: b(:,:)  ! (ncells, n_dots)
      real(rk), intent(out) :: results(:) ! (n_dots)
      real(rk) :: local_dots(n_dots)
      integer :: c, i, ierr

      local_dots = zero
      do i = 1, n_dots
         do c = flow%first_cell, flow%last_cell
            local_dots(i) = local_dots(i) + a(c, i) * b(c, i)
         end do
      end do

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_dots, results, n_dots, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce batched dots')
      call profiler_stop('MPI_Communication')
   end subroutine flow_global_dots_owned


   !> Computes two global dot products without constructing temporary full-size batches.
   subroutine flow_global_two_dots_owned(flow, a1, b1, a2, b2, results)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: a1(:), b1(:), a2(:), b2(:)
      real(rk), intent(out) :: results(2)
      real(rk) :: local_dots(2)
      integer :: c, ierr

      local_dots = zero
      do c = flow%first_cell, flow%last_cell
         local_dots(1) = local_dots(1) + a1(c) * b1(c)
         local_dots(2) = local_dots(2) + a2(c) * b2(c)
      end do

      if (flow%nprocs == 1) then
         results = local_dots
         return
      end if

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_dots, results, 2, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce two dots')
      call profiler_stop('MPI_Communication')
   end subroutine flow_global_two_dots_owned


   !> Computes the global dot product of two vectors over owned cells.
   function flow_global_dot_owned(flow, a, b) result(dot)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: a(:), b(:)
      real(rk) :: dot, local_dot
      integer :: c, ierr

      local_dot = 0.0_rk
      do c = flow%first_cell, flow%last_cell
         local_dot = local_dot + a(c) * b(c)
      end do

      if (flow%nprocs == 1) then
         dot = local_dot
         return
      end if

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_dot, dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow dot')
      call profiler_stop('MPI_Communication')
   end function flow_global_dot_owned


   !> Computes the global sum of a field over owned cells.
   function flow_global_sum_owned(flow, a) result(total)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: a(:)
      real(rk) :: total, local_total
      integer :: c, ierr

      local_total = 0.0_rk
      do c = flow%first_cell, flow%last_cell
         local_total = local_total + a(c)
      end do

      if (flow%nprocs == 1) then
         total = local_total
         return
      end if

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_total, total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow sum')
      call profiler_stop('MPI_Communication')
   end function flow_global_sum_owned


   !> Computes the global maximum magnitude of a field over owned cells.
   function flow_global_max_owned(flow, a) result(global_max)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: a(:)
      real(rk) :: global_max, local_max
      integer :: c, ierr

      local_max = 0.0_rk
      do c = flow%first_cell, flow%last_cell
         local_max = max(local_max, abs(a(c)))
      end do

      if (flow%nprocs == 1) then
         global_max = local_max
         return
      end if

      call profiler_start('MPI_Communication')
      call MPI_Allreduce(local_max, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow max')
      call profiler_stop('MPI_Communication')
   end function flow_global_max_owned


   !> Exchanges owned cell scalar values to ranks that keep them as ghosts.
   subroutine flow_exchange_cell_scalar(flow, field)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(inout) :: field(:)
      integer, parameter :: cell_halo_tag = 9281
      integer :: i, j, nreq, ierr, offset, count

      if (.not. allocated(flow%cell_sendbuf)) return
      if (flow%ncell_recv_ranks + flow%ncell_send_ranks == 0) return

      do i = 1, flow%ncell_send_ranks
         offset = flow%cell_send_displs(i)
         count = flow%cell_send_counts(i)
         do j = 1, count
            flow%cell_sendbuf(offset + j) = field(flow%cell_send_cells(offset + j))
         end do
      end do

      call profiler_start('MPI_Communication')
      nreq = 0
      do i = 1, flow%ncell_recv_ranks
         offset = flow%cell_recv_displs(i)
         count = flow%cell_recv_counts(i)
         nreq = nreq + 1
         call MPI_Irecv(flow%cell_recvbuf(offset + 1), count, MPI_DOUBLE_PRECISION, &
                        flow%cell_recv_ranks(i), cell_halo_tag, flow%comm, flow%cell_requests(nreq), ierr)
         call check_mpi(ierr, 'cell scalar halo irecv')
      end do
      do i = 1, flow%ncell_send_ranks
         offset = flow%cell_send_displs(i)
         count = flow%cell_send_counts(i)
         nreq = nreq + 1
         call MPI_Isend(flow%cell_sendbuf(offset + 1), count, MPI_DOUBLE_PRECISION, &
                        flow%cell_send_ranks(i), cell_halo_tag, flow%comm, flow%cell_requests(nreq), ierr)
         call check_mpi(ierr, 'cell scalar halo isend')
      end do
      call MPI_Waitall(nreq, flow%cell_requests(1:nreq), MPI_STATUSES_IGNORE, ierr)
      call check_mpi(ierr, 'cell scalar halo waitall')
      call profiler_stop('MPI_Communication')

      do i = 1, size(flow%cell_recv_cells)
         field(flow%cell_recv_cells(i)) = flow%cell_recvbuf(i)
      end do
   end subroutine flow_exchange_cell_scalar


   !> Exchanges owned cell matrix values to ranks that keep them as ghosts.
   subroutine flow_exchange_cell_matrix(flow, field)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(inout) :: field(:,:)
      integer, parameter :: cell_matrix_halo_tag = 9282
      integer :: ncomp, i, j, k, nreq, ierr, offset, count, pos

      if (.not. allocated(flow%cell_sendbuf)) return
      if (flow%ncell_recv_ranks + flow%ncell_send_ranks == 0) return

      ncomp = size(field, 1)
      if (ncomp > flow%cell_halo_max_components) then
         call fatal_error('mpi_flow', 'cell matrix halo component count exceeds cached buffer size')
      end if

      do i = 1, flow%ncell_send_ranks
         offset = flow%cell_send_displs(i)
         count = flow%cell_send_counts(i)
         do j = 1, count
            pos = (offset + j - 1) * ncomp
            do k = 1, ncomp
               flow%cell_sendbuf(pos + k) = field(k, flow%cell_send_cells(offset + j))
            end do
         end do
      end do

      call profiler_start('MPI_Communication')
      nreq = 0
      do i = 1, flow%ncell_recv_ranks
         offset = flow%cell_recv_displs(i) * ncomp
         count = flow%cell_recv_counts(i) * ncomp
         nreq = nreq + 1
         call MPI_Irecv(flow%cell_recvbuf(offset + 1), count, MPI_DOUBLE_PRECISION, &
                        flow%cell_recv_ranks(i), cell_matrix_halo_tag, flow%comm, flow%cell_requests(nreq), ierr)
         call check_mpi(ierr, 'cell matrix halo irecv')
      end do
      do i = 1, flow%ncell_send_ranks
         offset = flow%cell_send_displs(i) * ncomp
         count = flow%cell_send_counts(i) * ncomp
         nreq = nreq + 1
         call MPI_Isend(flow%cell_sendbuf(offset + 1), count, MPI_DOUBLE_PRECISION, &
                        flow%cell_send_ranks(i), cell_matrix_halo_tag, flow%comm, flow%cell_requests(nreq), ierr)
         call check_mpi(ierr, 'cell matrix halo isend')
      end do
      call MPI_Waitall(nreq, flow%cell_requests(1:nreq), MPI_STATUSES_IGNORE, ierr)
      call check_mpi(ierr, 'cell matrix halo waitall')
      call profiler_stop('MPI_Communication')

      do i = 1, size(flow%cell_recv_cells)
         pos = (i - 1) * ncomp
         do k = 1, ncomp
            field(k, flow%cell_recv_cells(i)) = flow%cell_recvbuf(pos + k)
         end do
      end do
   end subroutine flow_exchange_cell_matrix


   !> Exchanges owner-computed face scalar values to ranks owning the neighbor cell.
   subroutine flow_exchange_face_scalar(flow, face_field)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(inout) :: face_field(:)
      integer, parameter :: face_halo_tag = 9283
      integer :: i, j, nreq, ierr, offset, count

      if (.not. allocated(flow%face_sendbuf)) return
      if (flow%nface_recv_ranks + flow%nface_send_ranks == 0) return

      do i = 1, flow%nface_send_ranks
         offset = flow%face_send_displs(i)
         count = flow%face_send_counts(i)
         do j = 1, count
            flow%face_sendbuf(offset + j) = face_field(flow%face_send_faces(offset + j))
         end do
      end do

      call profiler_start('MPI_Communication')
      nreq = 0
      do i = 1, flow%nface_recv_ranks
         offset = flow%face_recv_displs(i)
         count = flow%face_recv_counts(i)
         nreq = nreq + 1
         call MPI_Irecv(flow%face_recvbuf(offset + 1), count, MPI_DOUBLE_PRECISION, &
                        flow%face_recv_ranks(i), face_halo_tag, flow%comm, flow%face_requests(nreq), ierr)
         call check_mpi(ierr, 'face halo irecv')
      end do
      do i = 1, flow%nface_send_ranks
         offset = flow%face_send_displs(i)
         count = flow%face_send_counts(i)
         nreq = nreq + 1
         call MPI_Isend(flow%face_sendbuf(offset + 1), count, MPI_DOUBLE_PRECISION, &
                        flow%face_send_ranks(i), face_halo_tag, flow%comm, flow%face_requests(nreq), ierr)
         call check_mpi(ierr, 'face halo isend')
      end do
      call MPI_Waitall(nreq, flow%face_requests(1:nreq), MPI_STATUSES_IGNORE, ierr)
      call check_mpi(ierr, 'face halo waitall')
      call profiler_stop('MPI_Communication')

      do i = 1, size(flow%face_recv_faces)
         face_field(flow%face_recv_faces(i)) = flow%face_recvbuf(i)
      end do
   end subroutine flow_exchange_face_scalar


   !> Gathers owned scalar cell values to rank 0 only.
   subroutine flow_gather_owned_scalar_root(flow, field, root_field)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(in) :: field(:)
      real(rk), intent(inout) :: root_field(:)
      integer :: ierr, r, first

      flow%gather_sendbuf = field(flow%first_cell:flow%last_cell)

      call profiler_start('MPI_Communication')
      call MPI_Gatherv(flow%gather_sendbuf, flow%nlocal, MPI_DOUBLE_PRECISION, &
                       flow%gather_recvbuf, flow%gather_counts, flow%gather_displs, &
                       MPI_DOUBLE_PRECISION, 0, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Gatherv owned scalar root')
      call profiler_stop('MPI_Communication')

      if (flow%rank == 0) then
         root_field = zero
         do r = 1, flow%nprocs
            first = flow%gather_firsts(r)
            root_field(first:first + flow%gather_counts(r) - 1) = &
               flow%gather_recvbuf(flow%gather_displs(r) + 1: &
                                   flow%gather_displs(r) + flow%gather_counts(r))
         end do
      end if
   end subroutine flow_gather_owned_scalar_root


   !> Gathers owned matrix cell values to rank 0 only.
   subroutine flow_gather_owned_matrix_root(flow, field, root_field)
      use mod_profiler, only : profiler_start, profiler_stop
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(in) :: field(:,:)
      real(rk), intent(inout) :: root_field(:,:)
      integer :: ierr, ncomp, ncells, nlocal_comp
      integer :: c, k, r, first, pos, recv_pos

      ncomp = size(field, 1)
      ncells = size(field, 2)
      call prepare_matrix_gather(flow, ncomp, ncells, nlocal_comp)

      pos = 0
      do c = flow%first_cell, flow%last_cell
         do k = 1, ncomp
            pos = pos + 1
            flow%gather_matrix_sendbuf(pos) = field(k, c)
         end do
      end do

      call profiler_start('MPI_Communication')
      call MPI_Gatherv(flow%gather_matrix_sendbuf, nlocal_comp, MPI_DOUBLE_PRECISION, &
                       flow%gather_matrix_recvbuf, flow%gather_matrix_counts, &
                       flow%gather_matrix_displs, MPI_DOUBLE_PRECISION, 0, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Gatherv owned matrix root')
      call profiler_stop('MPI_Communication')

      if (flow%rank == 0) then
         root_field = zero
         do r = 1, flow%nprocs
            first = flow%gather_firsts(r)
            recv_pos = flow%gather_matrix_displs(r)
            do c = first, first + flow%gather_counts(r) - 1
               do k = 1, ncomp
                  recv_pos = recv_pos + 1
                  root_field(k, c) = flow%gather_matrix_recvbuf(recv_pos)
               end do
            end do
         end do
      end if
   end subroutine flow_gather_owned_matrix_root


   !> Initializes contiguous-decomposition owner lookup.
   subroutine setup_cell_owners(mesh, flow)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      integer :: r, first, last

      allocate(flow%cell_owner(mesh%ncells))
      flow%cell_owner = -1
      do r = 1, flow%nprocs
         first = flow%gather_firsts(r)
         last = first + flow%gather_counts(r) - 1
         if (last >= first) flow%cell_owner(first:last) = r - 1
      end do
   end subroutine setup_cell_owners


   !> Caches faces whose owner cell belongs to this rank.
   subroutine setup_owned_faces(mesh, flow)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      integer :: f, n

      n = 0
      do f = 1, mesh%nfaces
         if (flow%owned(mesh%faces(f)%owner)) n = n + 1
      end do

      allocate(flow%owned_faces(n))
      n = 0
      do f = 1, mesh%nfaces
         if (.not. flow%owned(mesh%faces(f)%owner)) cycle
         n = n + 1
         flow%owned_faces(n) = f
      end do
   end subroutine setup_owned_faces


   !> Builds cell ghost send/receive metadata for one-ring neighbor stencils.
   subroutine setup_cell_halo(mesh, flow)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      logical, allocatable :: need(:,:)
      integer, allocatable :: recv_counts_all(:), recv_displs_all(:), recv_next(:)
      integer, allocatable :: send_counts_all(:), send_displs_all(:), send_next(:)
      integer :: q, c, lf, f, nb, owner, my_index, first, last
      integer :: total_recv, total_send, total_requests

      allocate(need(mesh%ncells, flow%nprocs))
      allocate(recv_counts_all(flow%nprocs), recv_displs_all(flow%nprocs), recv_next(flow%nprocs))
      allocate(send_counts_all(flow%nprocs), send_displs_all(flow%nprocs), send_next(flow%nprocs))

      need = .false.
      recv_counts_all = 0
      send_counts_all = 0

      do q = 1, flow%nprocs
         first = flow%gather_firsts(q)
         last = first + flow%gather_counts(q) - 1
         do c = first, last
            do lf = 1, mesh%ncell_faces(c)
               f = mesh%cell_faces(lf, c)
               nb = mesh_neighbor_for_cell(mesh, f, c)
               if (nb <= 0) cycle
               owner = flow%cell_owner(nb)
               if (owner /= q - 1) need(nb, q) = .true.
            end do
         end do
      end do

      my_index = flow%rank + 1
      do c = 1, mesh%ncells
         if (.not. need(c, my_index)) cycle
         owner = flow%cell_owner(c)
         recv_counts_all(owner + 1) = recv_counts_all(owner + 1) + 1
      end do

      do q = 1, flow%nprocs
         if (q == my_index) cycle
         do c = 1, mesh%ncells
            if (need(c, q) .and. flow%cell_owner(c) == flow%rank) then
               send_counts_all(q) = send_counts_all(q) + 1
            end if
         end do
      end do

      call prefix_counts(recv_counts_all, recv_displs_all)
      call prefix_counts(send_counts_all, send_displs_all)

      total_recv = sum(recv_counts_all)
      total_send = sum(send_counts_all)

      call pack_rank_metadata(recv_counts_all, recv_displs_all, &
                              flow%ncell_recv_ranks, flow%cell_recv_ranks, &
                              flow%cell_recv_counts, flow%cell_recv_displs)
      call pack_rank_metadata(send_counts_all, send_displs_all, &
                              flow%ncell_send_ranks, flow%cell_send_ranks, &
                              flow%cell_send_counts, flow%cell_send_displs)

      allocate(flow%cell_recv_cells(total_recv))
      allocate(flow%cell_send_cells(total_send))
      allocate(flow%ghost_cells(total_recv))
      allocate(flow%cell_recvbuf(max(1, total_recv * flow%cell_halo_max_components)))
      allocate(flow%cell_sendbuf(max(1, total_send * flow%cell_halo_max_components)))

      recv_next = recv_displs_all
      do c = 1, mesh%ncells
         if (.not. need(c, my_index)) cycle
         owner = flow%cell_owner(c) + 1
         recv_next(owner) = recv_next(owner) + 1
         flow%cell_recv_cells(recv_next(owner)) = c
      end do
      flow%ghost_cells = flow%cell_recv_cells

      send_next = send_displs_all
      do q = 1, flow%nprocs
         if (q == my_index) cycle
         do c = 1, mesh%ncells
            if (.not. need(c, q)) cycle
            if (flow%cell_owner(c) /= flow%rank) cycle
            send_next(q) = send_next(q) + 1
            flow%cell_send_cells(send_next(q)) = c
         end do
      end do

      total_requests = flow%ncell_recv_ranks + flow%ncell_send_ranks
      allocate(flow%cell_requests(max(1, total_requests)))
      flow%cell_recvbuf = zero
      flow%cell_sendbuf = zero

      deallocate(need)
      deallocate(recv_counts_all, recv_displs_all, recv_next)
      deallocate(send_counts_all, send_displs_all, send_next)
   end subroutine setup_cell_halo


   !> Builds face-flux halo metadata for shared internal faces.
   subroutine setup_face_halo(mesh, flow)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      logical, allocatable :: need(:,:)
      integer, allocatable :: recv_counts_all(:), recv_displs_all(:), recv_next(:)
      integer, allocatable :: send_counts_all(:), send_displs_all(:), send_next(:)
      integer :: q, c, lf, f, face_owner_rank, my_index, first, last
      integer :: total_recv, total_send, total_requests

      allocate(need(mesh%nfaces, flow%nprocs))
      allocate(recv_counts_all(flow%nprocs), recv_displs_all(flow%nprocs), recv_next(flow%nprocs))
      allocate(send_counts_all(flow%nprocs), send_displs_all(flow%nprocs), send_next(flow%nprocs))

      need = .false.
      recv_counts_all = 0
      send_counts_all = 0

      do q = 1, flow%nprocs
         first = flow%gather_firsts(q)
         last = first + flow%gather_counts(q) - 1
         do c = first, last
            do lf = 1, mesh%ncell_faces(c)
               f = mesh%cell_faces(lf, c)
               face_owner_rank = flow%cell_owner(mesh%faces(f)%owner)
               if (face_owner_rank /= q - 1) need(f, q) = .true.
            end do
         end do
      end do

      my_index = flow%rank + 1
      do f = 1, mesh%nfaces
         if (.not. need(f, my_index)) cycle
         face_owner_rank = flow%cell_owner(mesh%faces(f)%owner)
         recv_counts_all(face_owner_rank + 1) = recv_counts_all(face_owner_rank + 1) + 1
      end do

      do q = 1, flow%nprocs
         if (q == my_index) cycle
         do f = 1, mesh%nfaces
            if (need(f, q) .and. flow%cell_owner(mesh%faces(f)%owner) == flow%rank) then
               send_counts_all(q) = send_counts_all(q) + 1
            end if
         end do
      end do

      call prefix_counts(recv_counts_all, recv_displs_all)
      call prefix_counts(send_counts_all, send_displs_all)

      total_recv = sum(recv_counts_all)
      total_send = sum(send_counts_all)

      call pack_rank_metadata(recv_counts_all, recv_displs_all, &
                              flow%nface_recv_ranks, flow%face_recv_ranks, &
                              flow%face_recv_counts, flow%face_recv_displs)
      call pack_rank_metadata(send_counts_all, send_displs_all, &
                              flow%nface_send_ranks, flow%face_send_ranks, &
                              flow%face_send_counts, flow%face_send_displs)

      allocate(flow%face_recv_faces(total_recv))
      allocate(flow%face_send_faces(total_send))
      allocate(flow%face_recvbuf(max(1, total_recv)))
      allocate(flow%face_sendbuf(max(1, total_send)))

      recv_next = recv_displs_all
      do f = 1, mesh%nfaces
         if (.not. need(f, my_index)) cycle
         face_owner_rank = flow%cell_owner(mesh%faces(f)%owner) + 1
         recv_next(face_owner_rank) = recv_next(face_owner_rank) + 1
         flow%face_recv_faces(recv_next(face_owner_rank)) = f
      end do

      send_next = send_displs_all
      do q = 1, flow%nprocs
         if (q == my_index) cycle
         do f = 1, mesh%nfaces
            if (.not. need(f, q)) cycle
            if (flow%cell_owner(mesh%faces(f)%owner) /= flow%rank) cycle
            send_next(q) = send_next(q) + 1
            flow%face_send_faces(send_next(q)) = f
         end do
      end do

      total_requests = flow%nface_recv_ranks + flow%nface_send_ranks
      allocate(flow%face_requests(max(1, total_requests)))
      flow%face_recvbuf = zero
      flow%face_sendbuf = zero

      deallocate(need)
      deallocate(recv_counts_all, recv_displs_all, recv_next)
      deallocate(send_counts_all, send_displs_all, send_next)
   end subroutine setup_face_halo


   !> Returns the mesh neighbor, using stored periodic links when present.
   integer function mesh_neighbor_for_cell(mesh, face_id, cell_id) result(nb)
      type(mesh_t), intent(in) :: mesh
      integer, intent(in) :: face_id, cell_id

      if (mesh%faces(face_id)%owner == cell_id) then
         nb = mesh%faces(face_id)%neighbor
      else
         nb = mesh%faces(face_id)%owner
      end if

      if (nb == 0 .and. mesh%faces(face_id)%periodic_neighbor > 0) then
         nb = mesh%faces(face_id)%periodic_neighbor
      end if
   end function mesh_neighbor_for_cell


   !> Converts per-rank counts to zero-based displacements.
   subroutine prefix_counts(counts, displs)
      integer, intent(in) :: counts(:)
      integer, intent(out) :: displs(:)
      integer :: r

      if (size(counts) <= 0) return

      displs(1) = 0
      do r = 2, size(counts)
         displs(r) = displs(r - 1) + counts(r - 1)
      end do
   end subroutine prefix_counts


   !> Packs full per-rank metadata down to active communication partners.
   subroutine pack_rank_metadata(counts_all, displs_all, nactive, ranks, counts, displs)
      integer, intent(in) :: counts_all(:), displs_all(:)
      integer, intent(out) :: nactive
      integer, allocatable, intent(out) :: ranks(:), counts(:), displs(:)
      integer :: r, i

      nactive = count(counts_all > 0)
      allocate(ranks(nactive))
      allocate(counts(nactive))
      allocate(displs(nactive))

      i = 0
      do r = 1, size(counts_all)
         if (counts_all(r) <= 0) cycle
         i = i + 1
         ranks(i) = r - 1
         counts(i) = counts_all(r)
         displs(i) = displs_all(r)
      end do
   end subroutine pack_rank_metadata


   !> Internal helper for MPI error checking.
   subroutine check_mpi(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where
      if (ierr /= MPI_SUCCESS) call fatal_error('mpi_flow', trim(where)//' failed')
   end subroutine check_mpi

end module mod_mpi_flow
