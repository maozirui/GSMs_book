   subroutine time_elapsed(s)

!===============================================================================
!   The standard Fortran 90 routine RTC is used to calculate the elapsed CPU
!===============================================================================

!---    use dfport
    use config_parameter
    implicit none
    real(8) :: s

!---   s = rtc()
    call CPU_TIME(s)

   end subroutine

!   subroutine time_print
!===============================================================================
!   TIME_PRINT                     Print out the current date and time.
!
! . Notes:
!
!   The standard Fortran 90 routine DATE_AND_TIME is used to get the current
!   date and time strings.
!
!   No routines are included for the evaluation of the CPU time as these were
!   not part of the Fortran 90 standard. There is a CPU_TIME function, however,
!   in the Fortran 95 revision and this will be used in future versions of the
!   program.
!===============================================================================

!    implicit none
!    integer(4), parameter :: output = 6

!  Local scalars.
!   character ( len =  8 ) :: datstr
!   character ( len = 10 ) :: timstr

!  Get the current date and time.
!   call date_and_time ( datstr, timstr )

!  Write out the date and time.
!   write ( output, "(/A)"  ) "                  Date = " // datstr(7:8) // "/" // &
!                                          datstr(5:6) // "/" // &
!                                          datstr(1:4)
!   write ( output, "(A)"   ) "                  Time = " // timstr(1:2) // ":" // &
!                                          timstr(3:4) // ":" // &
!                                          timstr(5:10)
!   write ( output, *)

!   end subroutine time_print

