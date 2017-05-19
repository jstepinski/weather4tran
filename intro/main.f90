PROGRAM MAIN
!  Purpose: an example of a how to link subroutines and modules
!    Exclamation points indicate comments

  USE modf !MAIN now has access to all variables declared in modf.f90
  IMPLICIT NONE !do not use implicit data types

!********Declare Local Variables*******************
! You must declare any variables used in the subroutine 
! before you can do calculations
! (We are not using any variables except those from modf, 
!    so no need to define any others)

!********Begin Main Program*******************
  CALL cels2kelv !Converts data stored in modf from C to K

  write(*,*) 'Now find avg of ', num_temps, ' numbers' 
  write(*,*) "======================================"

  CALL my_avg !Call my_avg subroutine

  STOP
END
