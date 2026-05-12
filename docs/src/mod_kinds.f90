!> Core kind parameters, mathematical constants, and global utility routines.
!!
!! This module serves as the foundational precision manager for the solver. 
!! It centralizes the selection of floating-point precision (real64) and 
!! provides a set of universally available mathematical constants to 
!! avoid hard-coded magic numbers throughout the codebase. 
!!
!! Additionally, it provides critical infrastructure for:
!! 1. **Error Handling**: A standardized `fatal_error` mechanism to ensure 
!!    consistent error reporting and clean termination.
!! 2. **String Manipulation**: Utilities like `lowercase` for case-insensitive 
!!    namelist and file processing.
module mod_kinds
   use, intrinsic :: iso_fortran_env, only : real64, error_unit, output_unit
   implicit none

   private

   public :: rk, zero, one, two, half
   public :: tiny_safe, name_len, path_len
   public :: error_unit, output_unit
   public :: fatal_error, lowercase

   !> Working real precision parameter (Double Precision).
   !! All physical fields and numerical coefficients must use this kind.
   integer, parameter :: rk = real64

   real(rk), parameter :: zero = 0.0_rk !< Double precision constant 0.0
   real(rk), parameter :: one  = 1.0_rk !< Double precision constant 1.0
   real(rk), parameter :: two  = 2.0_rk !< Double precision constant 2.0
   real(rk), parameter :: half = 0.5_rk !< Double precision constant 0.5
   
   !> A very small positive value used for safe divisions and tolerance checks.
   !! Designed to be well above the underflow limit of double precision.
   real(rk), parameter :: tiny_safe = 1.0e-300_rk

   integer, parameter :: name_len = 64  !< Default length for labels, species names, etc.
   integer, parameter :: path_len = 256 !< Default length for directory paths and filenames.

contains

   !> Aborts the simulation with a formatted error message.
   !!
   !! This routine should be used for unrecoverable errors (e.g., missing 
   !! input files, converged solver failures). It writes to standard error 
   !! and calls `stop 1`.
   !!
   !! @param scope The name of the module or feature (e.g., "mod_input").
   !! @param message Descriptive message explaining why the simulation failed.
   subroutine fatal_error(scope, message)
      character(len=*), intent(in) :: scope
      character(len=*), intent(in) :: message

      write(error_unit,'(a,": ",a)') trim(scope), trim(message)
      stop 1
   end subroutine fatal_error


   !> Converts an input string to all lowercase characters.
   !!
   !! Useful for normalizing namelist inputs and performing case-insensitive 
   !! comparison of species names or boundary conditions.
   !!
   !! @param text The source string to convert.
   !! @return The converted lowercase string.
   pure function lowercase(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out

      integer :: i
      integer :: code

      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            out(i:i) = achar(code + iachar('a') - iachar('A'))
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lowercase

end module mod_kinds

