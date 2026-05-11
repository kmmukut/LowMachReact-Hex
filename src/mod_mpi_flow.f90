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
   public :: flow_global_dot_owned, flow_global_sum_owned, flow_global_max_owned
   public :: flow_allgather_owned_scalar

   type :: flow_mpi_t
      type(MPI_Comm) :: comm = MPI_COMM_NULL
      integer :: rank = -1
      integer :: nprocs = 0
      integer :: first_cell = 0
      integer :: last_cell = -1
      integer :: nlocal = 0
      logical, allocatable :: owned(:)

      ! Cached metadata/buffers for gathering owned cell ranges.
      integer, allocatable :: gather_counts(:)
      integer, allocatable :: gather_displs(:)
      integer, allocatable :: gather_firsts(:)
      real(rk), allocatable :: gather_sendbuf(:)
      real(rk), allocatable :: gather_recvbuf(:)
   end type flow_mpi_t

   logical :: mpi_started_here = .false.

contains

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


   subroutine flow_mpi_initialize(mesh, flow, comm_parent)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(MPI_Comm), intent(in) :: comm_parent

      integer :: ierr
      integer :: base, rem
      integer :: c

      call flow_mpi_finalize(flow)

      call MPI_Comm_dup(comm_parent, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Comm_dup flow')

      call MPI_Comm_rank(flow%comm, flow%rank, ierr)
      call check_mpi(ierr, 'MPI_Comm_rank flow')

      call MPI_Comm_size(flow%comm, flow%nprocs, ierr)
      call check_mpi(ierr, 'MPI_Comm_size flow')

      base = mesh%ncells / flow%nprocs
      rem = mod(mesh%ncells, flow%nprocs)

      if (flow%rank < rem) then
         flow%nlocal = base + 1
         flow%first_cell = flow%rank * (base + 1) + 1
      else
         flow%nlocal = base
         flow%first_cell = rem * (base + 1) + (flow%rank - rem) * base + 1
      end if

      flow%last_cell = flow%first_cell + flow%nlocal - 1

      allocate(flow%owned(mesh%ncells))
      flow%owned = .false.

      do c = flow%first_cell, flow%last_cell
         if (c >= 1 .and. c <= mesh%ncells) flow%owned(c) = .true.
      end do

      call setup_owned_gather(mesh, flow)
   end subroutine flow_mpi_initialize


   subroutine flow_mpi_finalize(flow)
      type(flow_mpi_t), intent(inout) :: flow

      integer :: ierr

      if (allocated(flow%owned)) deallocate(flow%owned)

      if (allocated(flow%gather_counts)) deallocate(flow%gather_counts)
      if (allocated(flow%gather_displs)) deallocate(flow%gather_displs)
      if (allocated(flow%gather_firsts)) deallocate(flow%gather_firsts)
      if (allocated(flow%gather_sendbuf)) deallocate(flow%gather_sendbuf)
      if (allocated(flow%gather_recvbuf)) deallocate(flow%gather_recvbuf)

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
   end subroutine flow_mpi_finalize


   subroutine setup_owned_gather(mesh, flow)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow

      integer :: ierr
      integer :: r
      integer :: total_count

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
   end subroutine setup_owned_gather


   subroutine flow_allreduce_global_vector(flow, local_values, global_values)
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: local_values(:,:)
      real(rk), intent(out) :: global_values(:,:)

      integer :: ierr

      call MPI_Allreduce(local_values, global_values, size(local_values), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow vector')
   end subroutine flow_allreduce_global_vector


   subroutine flow_allreduce_global_scalar(flow, local_values, global_values)
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: local_values(:)
      real(rk), intent(out) :: global_values(:)

      integer :: ierr

      call MPI_Allreduce(local_values, global_values, size(local_values), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow scalar')
   end subroutine flow_allreduce_global_scalar


   subroutine flow_allgather_owned_scalar(flow, local_global, global)
      type(flow_mpi_t), intent(inout) :: flow
      real(rk), intent(in) :: local_global(:)
      real(rk), intent(out) :: global(:)

      integer :: ierr
      integer :: r
      integer :: first

      if (.not. allocated(flow%gather_sendbuf)) then
         call fatal_error('mpi_flow', 'owned gather buffers are not initialized')
      end if

      flow%gather_sendbuf = local_global(flow%first_cell:flow%last_cell)

      call MPI_Allgatherv(flow%gather_sendbuf, flow%nlocal, MPI_DOUBLE_PRECISION, &
                          flow%gather_recvbuf, flow%gather_counts, flow%gather_displs, &
                          MPI_DOUBLE_PRECISION, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allgatherv owned scalar')

      global = zero

      do r = 1, flow%nprocs
         first = flow%gather_firsts(r)

         global(first:first + flow%gather_counts(r) - 1) = &
            flow%gather_recvbuf(flow%gather_displs(r) + 1: &
                                flow%gather_displs(r) + flow%gather_counts(r))
      end do
   end subroutine flow_allgather_owned_scalar


   function flow_global_dot_owned(flow, a, b) result(dot)
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: a(:)
      real(rk), intent(in) :: b(:)
      real(rk) :: dot

      real(rk) :: local_dot
      integer :: c, ierr

      local_dot = 0.0_rk
      do c = flow%first_cell, flow%last_cell
         local_dot = local_dot + a(c) * b(c)
      end do

      call MPI_Allreduce(local_dot, dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow dot')
   end function flow_global_dot_owned


   function flow_global_sum_owned(flow, a) result(total)
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: a(:)
      real(rk) :: total

      real(rk) :: local_total
      integer :: c, ierr

      local_total = 0.0_rk
      do c = flow%first_cell, flow%last_cell
         local_total = local_total + a(c)
      end do

      call MPI_Allreduce(local_total, total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow sum')
   end function flow_global_sum_owned


   function flow_global_max_owned(flow, a) result(global_max)
      type(flow_mpi_t), intent(in) :: flow
      real(rk), intent(in) :: a(:)
      real(rk) :: global_max

      real(rk) :: local_max
      integer :: c, ierr

      local_max = 0.0_rk
      do c = flow%first_cell, flow%last_cell
         local_max = max(local_max, abs(a(c)))
      end do

      call MPI_Allreduce(local_max, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      call check_mpi(ierr, 'MPI_Allreduce flow max')
   end function flow_global_max_owned


   subroutine check_mpi(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where

      if (ierr /= MPI_SUCCESS) call fatal_error('mpi_flow', trim(where)//' failed')
   end subroutine check_mpi

end module mod_mpi_flow