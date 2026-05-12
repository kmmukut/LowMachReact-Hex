!> Performance profiling and execution timing utilities.
!!
!! This module provides a simple hierarchical profiling system to measure 
!! wall-clock time spent in different parts of the solver (kernels). 
!! It uses `MPI_Wtime` for high-resolution timing and provides collective 
!! reporting across all MPI ranks.
module mod_profiler
   use mod_kinds, only: rk
   implicit none

   private

   public :: profiler_start, profiler_stop, profiler_report

   !> Maximum number of unique timers that can be tracked.
   integer, parameter :: MAX_TIMERS = 100
   
   !> Maximum length for a timer name.
   integer, parameter :: NAME_LEN = 32
   
   !> Internal structure to store timing data for a specific code block.
   type :: profiler_timer_t
      character(len=NAME_LEN) :: name = '' !< Unique identifier for the kernel.
      integer :: calls = 0                 !< Total number of times this timer was triggered.
      real(rk) :: total_time = 0.0_rk      !< Accumulated wall time [s].
      real(rk) :: start_time = 0.0_rk      !< Timestamp of the most recent start call.
      logical :: active = .false.          !< Status flag to prevent mismatched start/stop.
   end type profiler_timer_t

   !> Global registry of active timers.
   type(profiler_timer_t), save :: timers(MAX_TIMERS)
   
   !> Current count of registered timers.
   integer, save :: ntimers = 0

contains

   !> Starts a timer for a named kernel.
   !!
   !! If the name is seen for the first time, a new entry is created in the 
   !! global registry.
   !!
   !! @param name Human-readable identifier for the code block.
   subroutine profiler_start(name)
      character(len=*), intent(in) :: name
      integer :: i, idx
      real(8) :: mpi_wtime

      idx = 0
      do i = 1, ntimers
         if (trim(timers(i)%name) == trim(name)) then
            idx = i
            exit
         end if
      end do

      if (idx == 0) then
         if (ntimers < MAX_TIMERS) then
            ntimers = ntimers + 1
            idx = ntimers
            timers(idx)%name = trim(name)
         else
            return
         end if
      end if

      timers(idx)%start_time = real(mpi_wtime(), rk)
      timers(idx)%active = .true.
   end subroutine profiler_start


   !> Stops a timer and accumulates the elapsed time.
   !!
   !! @param name Name of the kernel to stop. Must match an earlier start call.
   subroutine profiler_stop(name)
      character(len=*), intent(in) :: name
      integer :: i, idx
      real(rk) :: elapsed
      real(8) :: mpi_wtime

      idx = 0
      do i = 1, ntimers
         if (trim(timers(i)%name) == trim(name)) then
            idx = i
            exit
         end if
      end do

      if (idx > 0 .and. timers(idx)%active) then
         elapsed = real(mpi_wtime(), rk) - timers(idx)%start_time
         timers(idx)%total_time = timers(idx)%total_time + elapsed
         timers(idx)%calls = timers(idx)%calls + 1
         timers(idx)%active = .false.
      end if
   end subroutine profiler_stop


   !> Generates a collective performance report across all MPI ranks.
   !!
   !! Performs global reductions to find the minimum, maximum, and average 
   !! time spent in each kernel across the entire communicator.
   !!
   !! @param comm MPI communicator.
   !! @param rank Local rank ID.
   !! @param nprocs Total number of processors.
   subroutine profiler_report(comm, rank, nprocs)
      use mpi_f08
      type(MPI_Comm), intent(in) :: comm
      integer, intent(in) :: rank
      integer, intent(in) :: nprocs
      
      integer :: i
      integer :: ierr
      real(8) :: global_min, global_max, global_sum
      real(8) :: local_time
      real(8) :: avg_time

      if (rank == 0) then
         write(*,*) '======================================================================'
         write(*,*) ' PERFORMANCE PROFILING REPORT (Wall Time in Seconds)'
         write(*,*) '======================================================================'
         write(*,'(a25,a10,a15,a15,a15)') 'Kernel Name', 'Calls', 'Min', 'Max', 'Avg'
         write(*,*) '----------------------------------------------------------------------'
      end if

      do i = 1, ntimers
         local_time = real(timers(i)%total_time, 8)
         
         call MPI_Reduce(local_time, global_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, comm, ierr)
         call MPI_Reduce(local_time, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
         call MPI_Reduce(local_time, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
         
         if (rank == 0) then
            avg_time = global_sum / real(nprocs, 8)
            write(*,'(a25,i10,f15.6,f15.6,f15.6)') trim(timers(i)%name), timers(i)%calls, global_min, global_max, avg_time
         end if
      end do

      if (rank == 0) then
         write(*,*) '======================================================================'
      end if

   end subroutine profiler_report

end module mod_profiler
