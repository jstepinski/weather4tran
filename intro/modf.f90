MODULE modf
!  Purpose: an example of a module used in subroutines to 
!    share common data. This is the Fortran90 equivalent of 
!    an F77 common block

  IMPLICIT NONE !variable types must be declared in advance, i.e. prevent vars starting with A-H from being reals automatically
  SAVE !Make sure this is always in your module file

!********Declare Variables*******************
! "parameter" means values cannot change (constant)
  real, parameter :: c2k = 273.15 !degK =degC + 273.15
  integer, parameter :: num_temps  =3
! "dimension(n)" makes an array of length n, where n must be defined as
!  a parameter (constant), like num_temps above
! Array of temps in Celsius
  real, dimension(num_temps) :: c_temps = (/22.0, 24.0, 26.3/) 
! Array of temps in Kelvin, initialized to 0., modified in cels2kelv
  real, dimension(num_temps) :: k_temps=0. 
  
END MODULE modf
