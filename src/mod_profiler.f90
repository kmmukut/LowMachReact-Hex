!> Performance profiling and execution timing utilities.
!!
!! Supports:
!! - Flat timing report for all named regions.
!! - Optional nested report based on profiler_start/profiler_stop call stack.
!! - Runtime enable/disable via profiler_configure().
!!
!! Timings are inclusive. Flat profiler rows are not additive when nested
!! timers are enabled. Current top-level timer names include
!! `Transport_Update`, `Projection_Step`, `Species_Transport`,
!! `Energy_Transport`, `Diagnostics_Write_Flow`, `Diagnostics_Write_Energy`,
!! and `Output_Write_VTU`; energy Cantera sync timers are
!! `Energy_Cantera_PreSync` and `Energy_Cantera_PostSync`.
module mod_profiler
   use mpi_f08
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   use mod_kinds, only : rk
   implicit none

   private

   public :: profiler_start
   public :: profiler_stop
   public :: profiler_report
   public :: profiler_configure
   public :: profiler_reset

   integer, parameter :: MAX_TIMERS = 512
   integer, parameter :: MAX_EDGES = 4096
   integer, parameter :: MAX_STACK = 128
   integer, parameter :: NAME_LEN = 96

   type :: profiler_timer_t
      character(len=NAME_LEN) :: name = ''
      integer :: calls = 0
      real(rk) :: total_time = 0.0_rk
   end type profiler_timer_t

   type :: profiler_edge_t
      integer :: parent = 0
      integer :: child = 0
      integer :: calls = 0
      real(rk) :: total_time = 0.0_rk
   end type profiler_edge_t

   type(profiler_timer_t), save :: timers(MAX_TIMERS)
   type(profiler_edge_t), save :: edges(MAX_EDGES)

   integer, save :: ntimers = 0
   integer, save :: nedges = 0

   integer, save :: stack_ids(MAX_STACK) = 0
   real(rk), save :: stack_start(MAX_STACK) = 0.0_rk
   integer, save :: stack_depth = 0

   logical, save :: profiling_enabled = .true.
   logical, save :: nested_enabled = .true.

contains

   !> Configure profiling behavior at runtime.
   subroutine profiler_configure(enabled, nested)
      logical, intent(in) :: enabled
      logical, intent(in) :: nested

      profiling_enabled = enabled
      nested_enabled = nested

      if (.not. profiling_enabled) then
         call profiler_reset()
      end if
   end subroutine profiler_configure


   !> Reset all profiler state.
   subroutine profiler_reset()
      integer :: i

      ntimers = 0
      nedges = 0
      stack_depth = 0
      stack_ids = 0
      stack_start = 0.0_rk

      do i = 1, MAX_TIMERS
         timers(i)%name = ''
         timers(i)%calls = 0
         timers(i)%total_time = 0.0_rk
      end do

      do i = 1, MAX_EDGES
         edges(i)%parent = 0
         edges(i)%child = 0
         edges(i)%calls = 0
         edges(i)%total_time = 0.0_rk
      end do
   end subroutine profiler_reset


   !> Starts a timer for a named kernel.
   subroutine profiler_start(name)
      character(len=*), intent(in) :: name
      integer :: idx

      if (.not. profiling_enabled) return

      idx = find_or_create_timer(name)

      if (stack_depth >= MAX_STACK) then
         write(error_unit,'(a)') 'profiler: nesting stack overflow'
         error stop 1
      end if

      stack_depth = stack_depth + 1
      stack_ids(stack_depth) = idx
      stack_start(stack_depth) = real(MPI_Wtime(), rk)
   end subroutine profiler_start


   !> Stops a timer and accumulates the elapsed time.
   subroutine profiler_stop(name)
      character(len=*), intent(in) :: name
      integer :: idx, parent
      real(rk) :: elapsed

      if (.not. profiling_enabled) return

      if (stack_depth <= 0) then
         write(error_unit,'(a,a)') 'profiler: stop with empty stack: ', trim(name)
         error stop 1
      end if

      idx = find_or_create_timer(name)

      if (stack_ids(stack_depth) /= idx) then
         write(error_unit,'(a)') 'profiler: mismatched profiler_stop'
         write(error_unit,'(a,a)') '  expected: ', trim(timers(stack_ids(stack_depth))%name)
         write(error_unit,'(a,a)') '  got:      ', trim(name)
         error stop 1
      end if

      elapsed = real(MPI_Wtime(), rk) - stack_start(stack_depth)

      timers(idx)%total_time = timers(idx)%total_time + elapsed
      timers(idx)%calls = timers(idx)%calls + 1

      if (nested_enabled) then
         if (stack_depth > 1) then
            parent = stack_ids(stack_depth - 1)
         else
            parent = 0
         end if
         call record_edge(parent, idx, elapsed)
      end if

      stack_ids(stack_depth) = 0
      stack_start(stack_depth) = 0.0_rk
      stack_depth = stack_depth - 1
   end subroutine profiler_stop


   !> Generates a collective performance report across all MPI ranks.
   subroutine profiler_report(comm, rank, nprocs)
      type(MPI_Comm), intent(in) :: comm
      integer, intent(in) :: rank
      integer, intent(in) :: nprocs

      integer :: i, ierr
      integer :: global_ntimers, global_nedges
      integer :: local_calls, global_calls
      real(rk) :: local_time
      real(rk) :: global_min, global_max, global_sum
      real(rk) :: avg_percent
      real(rk) :: total_avg_time
      real(rk), parameter :: tiny_time = 1.0e-300_rk

      real(rk) :: timer_min(MAX_TIMERS)
      real(rk) :: timer_max(MAX_TIMERS)
      real(rk) :: timer_avg(MAX_TIMERS)
      integer :: timer_calls(MAX_TIMERS)

      real(rk) :: edge_avg(MAX_EDGES)
      integer :: edge_calls(MAX_EDGES)

      if (.not. profiling_enabled) return

      if (stack_depth /= 0) then
         if (rank == 0) then
            write(error_unit,'(a,i0)') 'profiler: warning, nonzero stack depth at report: ', stack_depth
         end if
      end if

      call MPI_Allreduce(ntimers, global_ntimers, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
      call check_mpi(ierr, 'profiler ntimers')

      if (global_ntimers > MAX_TIMERS) then
         if (rank == 0) write(error_unit,'(a)') 'profiler: too many timers in report'
         error stop 1
      end if

      if (nested_enabled) then
         call MPI_Allreduce(nedges, global_nedges, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
         call check_mpi(ierr, 'profiler nedges')
      else
         global_nedges = 0
      end if

      timer_min = 0.0_rk
      timer_max = 0.0_rk
      timer_avg = 0.0_rk
      timer_calls = 0
      total_avg_time = tiny_time

      do i = 1, global_ntimers
         if (i <= ntimers) then
            local_time = timers(i)%total_time
            local_calls = timers(i)%calls
         else
            local_time = 0.0_rk
            local_calls = 0
         end if

         call MPI_Allreduce(local_time, global_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
         call check_mpi(ierr, 'profiler min')

         call MPI_Allreduce(local_time, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
         call check_mpi(ierr, 'profiler max')

         call MPI_Allreduce(local_time, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
         call check_mpi(ierr, 'profiler sum')

         call MPI_Allreduce(local_calls, global_calls, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
         call check_mpi(ierr, 'profiler calls')

         timer_min(i) = global_min
         timer_max(i) = global_max
         timer_avg(i) = global_sum / max(real(nprocs, rk), tiny_time)
         timer_calls(i) = global_calls

         if (i <= ntimers) then
            if (trim(timers(i)%name) == 'Total_Simulation') then
               total_avg_time = max(timer_avg(i), tiny_time)
            end if
         end if
      end do

      edge_avg = 0.0_rk
      edge_calls = 0

      if (nested_enabled) then
         do i = 1, global_nedges
            if (i <= nedges) then
               local_time = edges(i)%total_time
               local_calls = edges(i)%calls
            else
               local_time = 0.0_rk
               local_calls = 0
            end if

            call MPI_Allreduce(local_time, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
            call check_mpi(ierr, 'profiler edge sum')

            call MPI_Allreduce(local_calls, global_calls, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
            call check_mpi(ierr, 'profiler edge calls')

            edge_avg(i) = global_sum / max(real(nprocs, rk), tiny_time)
            edge_calls(i) = global_calls
         end do
      end if

      if (rank /= 0) return

      write(output_unit,'(a)') ' ======================================================================'
      write(output_unit,'(a)') '  PERFORMANCE PROFILING REPORT'
      write(output_unit,'(a)') '  Inclusive wall time in seconds; Avg% is relative to Total_Simulation.'
      write(output_unit,'(a)') '  Flat rows are not additive when nested timers are enabled.'
      write(output_unit,'(a)') ' ======================================================================'
      write(output_unit,'(a)') '                      Kernel Name     Calls            Min            Max            Avg        Avg%'
      write(output_unit,'(a)') ' ----------------------------------------------------------------------'

      do i = 1, global_ntimers
         avg_percent = 100.0_rk * timer_avg(i) / max(total_avg_time, tiny_time)

         if (i <= ntimers) then
            write(output_unit,'(1x,a32,1x,i9,3(1x,f14.6),1x,f10.2)') &
               trim(timers(i)%name), timer_calls(i), timer_min(i), timer_max(i), timer_avg(i), avg_percent
         else
            write(output_unit,'(1x,a32,1x,i9,3(1x,f14.6),1x,f10.2)') &
               'UNKNOWN_TIMER', timer_calls(i), timer_min(i), timer_max(i), timer_avg(i), avg_percent
         end if
      end do

      if (nested_enabled) then
         write(output_unit,'(a)') ' ======================================================================'
         write(output_unit,'(a)') '  NESTED PROFILING REPORT (Inclusive Child Time, Avg Across Ranks)'
         write(output_unit,'(a)') ' ======================================================================'
         write(output_unit,'(a)') '  Region Tree                                              Calls          Avg        Avg%'
         write(output_unit,'(a)') ' ----------------------------------------------------------------------'
         call print_children(0, 0, global_nedges, edge_avg, edge_calls, total_avg_time)
      end if

      write(output_unit,'(a)') ' ======================================================================'

   contains

      recursive subroutine print_children(parent, depth, edge_count, edge_avg_in, edge_calls_in, total_time)
         integer, intent(in) :: parent
         integer, intent(in) :: depth
         integer, intent(in) :: edge_count
         real(rk), intent(in) :: edge_avg_in(:)
         integer, intent(in) :: edge_calls_in(:)
         real(rk), intent(in) :: total_time

         integer :: e, child
         real(rk) :: pct
         character(len=256) :: label

         if (depth > 32) return

         do e = 1, edge_count
            if (e > nedges) cycle
            if (edges(e)%parent /= parent) cycle

            child = edges(e)%child
            if (child <= 0 .or. child > ntimers) cycle
            if (child == parent) cycle

            label = repeat('  ', depth)//'- '//trim(timers(child)%name)
            pct = 100.0_rk * edge_avg_in(e) / max(total_time, tiny_time)

            write(output_unit,'(2x,a52,1x,i9,1x,f12.6,1x,f10.2)') &
               label, edge_calls_in(e), edge_avg_in(e), pct

            call print_children(child, depth + 1, edge_count, edge_avg_in, edge_calls_in, total_time)
         end do
      end subroutine print_children

   end subroutine profiler_report


   integer function find_or_create_timer(name) result(idx)
      character(len=*), intent(in) :: name
      integer :: i

      do i = 1, ntimers
         if (trim(timers(i)%name) == trim(name)) then
            idx = i
            return
         end if
      end do

      if (ntimers >= MAX_TIMERS) then
         write(error_unit,'(a)') 'profiler: MAX_TIMERS exceeded'
         error stop 1
      end if

      ntimers = ntimers + 1
      timers(ntimers)%name = trim(name)
      timers(ntimers)%calls = 0
      timers(ntimers)%total_time = 0.0_rk

      idx = ntimers
   end function find_or_create_timer


   !> Updates the call tree edge statistics between two timers.
   !!
   !! @param parent Index of the calling timer.
   !! @param child Index of the called timer.
   !! @param elapsed Wall time spent in this call.
   subroutine record_edge(parent, child, elapsed)
      integer, intent(in) :: parent
      integer, intent(in) :: child
      real(rk), intent(in) :: elapsed

      integer :: e

      do e = 1, nedges
         if (edges(e)%parent == parent .and. edges(e)%child == child) then
            edges(e)%calls = edges(e)%calls + 1
            edges(e)%total_time = edges(e)%total_time + elapsed
            return
         end if
      end do

      if (nedges >= MAX_EDGES) then
         write(error_unit,'(a)') 'profiler: MAX_EDGES exceeded'
         error stop 1
      end if

      nedges = nedges + 1
      edges(nedges)%parent = parent
      edges(nedges)%child = child
      edges(nedges)%calls = 1
      edges(nedges)%total_time = elapsed
   end subroutine record_edge


   !> Internal utility to check MPI return codes and abort on failure.
   !!
   !! @param ierr The MPI return code.
   !! @param where Descriptive string of where the failure occurred.
   subroutine check_mpi(ierr, where)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: where

      if (ierr /= MPI_SUCCESS) then
         write(error_unit,'(a,a)') 'profiler MPI failure: ', trim(where)
         error stop 1
      end if
   end subroutine check_mpi

end module mod_profiler
