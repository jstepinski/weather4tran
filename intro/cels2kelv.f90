SUBROUTINE cels2kelv
! Purpose: 
!    Converts an array of temperatures from Celsius to Kelvin
!    Requires the use of variables defined in modf.f90
!    Changes the values of k_temps

  use modf  !uses the variables defined in modf.f90
  IMPLICIT NONE

!********Declare Local Variables*******************
  INTEGER :: I
!********End Declare Variables*******************

  do I=1,num_temps
    k_temps(I) = c_temps(I) + c2k ! degK =degC + 273.15
    write(*,*) c_temps(I), 'deg C = ', k_temps(I), 'deg K' 
      !example of printing to standard output
  end do

  RETURN
END
