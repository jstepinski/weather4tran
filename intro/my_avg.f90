SUBROUTINE my_avg
! Purpose: 
!    Finds the average of an array
!    Requires the use of variables defined in modf.f90
  USE modf !uses variables declared in modf.f90
  IMPLICIT NONE !do not use implicit data types

! You must declare any variables used in the subroutine 
! before you can do calculations
! These variables are NOT already defined in modf.f90
!********Declare Local Variables*******************
  REAL :: average=0. !average is the output
  REAL :: my_sum=0.
  INTEGER :: I

!********Begin Subroutine*******************
  do I=1,num_temps  !we can use num_temps from modf.f90
    my_sum = my_sum + k_temps(I) !sum the array k_temps defined in modf.f90
  end do
  average = my_sum/real(num_temps) !calculate average

  write(*,300) average !Print using format ID 300 
  write(*,301) average !Print using format ID 301 
  300 format(1X, 'Average = ', ES9.2, ' deg Kelvin') !uses sci. notation
  301 format(1X, 'Average = ', F9.2, ' deg Kelvin') !uses decimal notation


  RETURN
END
