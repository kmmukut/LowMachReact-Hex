!> MPI infrastructure for the radiation solver.
!!
!! This module manages the parallel execution of the radiation transport 
!! solver. Unlike the flow solver which uses a replicated mesh, the 
!! radiation solver can be configured with its own communicator and 
!! task decomposition logic (e.g., distributing rays or discrete ordinates).
module mod_mpi_radiation
   use mpi_f08
   use mod_kinds, only : fatal_error
   implicit none

   private

   public :: radiation_mpi_t
   public :: radiation_mpi_initialize, radiation_mpi_finalize
   public :: radiation_task_bounds

   !> MPI context for radiation operations.
   type :: radiation_mpi_t
      type(MPI_Comm) :: comm = MPI_COMM_NULL !< MPI Communicator.
      integer :: rank = -1                   !< Local rank ID.
      integer :: nprocs = 0                  !< Total number of radiation processors.
      integer :: first_task = 0              !< First task index owned by this rank.
      integer :: last_task = -1              !< Last task index owned by this rank.
      integer :: nlocal_tasks = 0            !< Total number of tasks owned locally.
   end type radiation_mpi_t

contains

   !> Initializes the radiation MPI context by duplicating a parent communicator.
   !!
   !! @param rad The radiation context to initialize.
   !! @param comm_parent The parent communicator (usually MPI_COMM_WORLD).
   !! @param n_tasks Optional total number of tasks to decompose immediately.
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


   !> Releases radiation MPI resources.
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


   !> Performs contiguous task decomposition for the radiation solver.
   !!
   !! @param rad The radiation context to update.
   !! @param n_tasks Total number of discrete tasks (e.g., rays).
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


   !> Internal helper for MPI error checking.
   subroutine check_mpi(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where

      if (ierr /= MPI_SUCCESS) call fatal_error('mpi_radiation', trim(where)//' failed')
   end subroutine check_mpi

end module mod_mpi_radiation

