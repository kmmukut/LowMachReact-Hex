module mod_mpi_radiation
   use mpi_f08
   use mod_kinds, only : fatal_error
   implicit none

   private

   public :: radiation_mpi_t
   public :: radiation_mpi_initialize, radiation_mpi_finalize
   public :: radiation_task_bounds

   type :: radiation_mpi_t
      type(MPI_Comm) :: comm = MPI_COMM_NULL
      integer :: rank = -1
      integer :: nprocs = 0
      integer :: first_task = 0
      integer :: last_task = -1
      integer :: nlocal_tasks = 0
   end type radiation_mpi_t

contains

   subroutine radiation_mpi_initialize(rad, comm_parent, n_tasks)
      type(radiation_mpi_t), intent(inout) :: rad
      type(MPI_Comm), intent(in) :: comm_parent
      integer, intent(in), optional :: n_tasks

      integer :: ierr

      call radiation_mpi_finalize(rad)

      call MPI_Comm_dup(comm_parent, rad%comm, ierr)
      call check_mpi(ierr, 'MPI_Comm_dup radiation')

      call MPI_Comm_rank(rad%comm, rad%rank, ierr)
      call check_mpi(ierr, 'MPI_Comm_rank radiation')

      call MPI_Comm_size(rad%comm, rad%nprocs, ierr)
      call check_mpi(ierr, 'MPI_Comm_size radiation')

      if (present(n_tasks)) then
         call radiation_task_bounds(rad, n_tasks)
      end if
   end subroutine radiation_mpi_initialize


   subroutine radiation_mpi_finalize(rad)
      type(radiation_mpi_t), intent(inout) :: rad

      integer :: ierr

      if (rad%comm /= MPI_COMM_NULL) then
         call MPI_Comm_free(rad%comm, ierr)
         call check_mpi(ierr, 'MPI_Comm_free radiation')
      end if

      rad%comm = MPI_COMM_NULL
      rad%rank = -1
      rad%nprocs = 0
      rad%first_task = 0
      rad%last_task = -1
      rad%nlocal_tasks = 0
   end subroutine radiation_mpi_finalize


   subroutine radiation_task_bounds(rad, n_tasks)
      type(radiation_mpi_t), intent(inout) :: rad
      integer, intent(in) :: n_tasks

      integer :: base, rem

      if (n_tasks < 0) call fatal_error('mpi_radiation', 'n_tasks cannot be negative')
      if (rad%nprocs <= 0) call fatal_error('mpi_radiation', 'radiation communicator is not initialized')

      base = n_tasks / rad%nprocs
      rem = mod(n_tasks, rad%nprocs)

      if (rad%rank < rem) then
         rad%nlocal_tasks = base + 1
         rad%first_task = rad%rank * (base + 1) + 1
      else
         rad%nlocal_tasks = base
         rad%first_task = rem * (base + 1) + (rad%rank - rem) * base + 1
      end if
      rad%last_task = rad%first_task + rad%nlocal_tasks - 1
   end subroutine radiation_task_bounds


   subroutine check_mpi(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where

      if (ierr /= MPI_SUCCESS) call fatal_error('mpi_radiation', trim(where)//' failed')
   end subroutine check_mpi

end module mod_mpi_radiation

