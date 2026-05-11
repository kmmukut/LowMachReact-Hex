!> Core kind parameters and utility routines.
!!
!! This module defines precision parameters, mathematical constants,
!! and error handling utilities.
module mod_kinds
   use, intrinsic :: iso_fortran_env, only : real64, error_unit, output_unit
   implicit none

   private

   public :: rk, zero, one, two, half
   public :: tiny_safe, name_len, path_len
   public :: error_unit, output_unit
   public :: fatal_error, lowercase

   integer, parameter :: rk = real64 !< Working real precision (double precision)

   real(rk), parameter :: zero = 0.0_rk
   real(rk), parameter :: one  = 1.0_rk
   real(rk), parameter :: two  = 2.0_rk
   real(rk), parameter :: half = 0.5_rk
   real(rk), parameter :: tiny_safe = 1.0e-300_rk

   integer, parameter :: name_len = 64
   integer, parameter :: path_len = 256

contains

   !> Print an error message and stop execution.
   !!
   !! @param scope The name of the module or feature where the error occurred.
   !! @param message Descriptive error message.
   subroutine fatal_error(scope, message)
      character(len=*), intent(in) :: scope
      character(len=*), intent(in) :: message

      write(error_unit,'(a,": ",a)') trim(scope), trim(message)
      stop 1
   end subroutine fatal_error


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

