SUBROUTINE interp1d(x,y,xi,yi,N)

! Linear interpolation
! Input vectors x and y
! Input a point xi within the domain of x
! Interpolate over x to find yi corresponding to xi

! N is the size of x and y
  
  IMPLICIT NONE
  
  INTEGER :: N
  REAL, DIMENSION(N) :: x
  REAL, DIMENSION(N) :: y
  REAL :: xi
  REAL :: yi

  REAL :: x1 = 0
  REAL :: x2 = 0
  INTEGER :: x1ind = 1
  INTEGER :: x2ind = 1
  REAL :: y1 = 0
  REAL :: y2 = 0
  REAL :: m = 0

  ! Find x1 and x2 within x such that x1 < x < x2
  ! Find the indices of x1 and x2 within x  
  CALL find2els(x,xi,N,x1,x2,x1ind,x2ind)
  
  ! Extract y1 and y2 corresponding to x1 and x2, respectively
  y1 = y(x1ind);
  y2 = y(x2ind);
  
  m = (y2 - y1)/(x2 - x1)    ! Slope
  yi = (xi - x2)*m + y2      ! Result
  
  RETURN
END

