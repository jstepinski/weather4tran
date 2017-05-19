SUBROUTINE find2els(v,x,N,a,b,aindx,bindx)

! Subroutine useful for linear interpolation
! Given a vector v, find the elements a and b such that x is between a and b
! Return a,b, and their indices aindx and bindx in v

  IMPLICIT NONE

  INTEGER :: N
  REAL :: x
  REAL, DIMENSION(N) :: v
  REAL :: a
  REAL :: b
  INTEGER :: aindx
  INTEGER :: bindx
  
  REAL, DIMENSION(N) :: res
  REAL :: mn = 0
  INTEGER :: mn_ind = 1
  INTEGER :: mx_ind = 1
  
  INTEGER :: i
  
  ! Compute absolute difference between x and every element of v
  DO i = 1,N
    res(i) = abs(x-v(i))
  ENDDO
  
  ! Find the minimum difference between x and v. Store the minimum in mn, and 
  ! store its index in mn_ind
  mn = res(1)
  DO i = 2,N
    IF (res(i) .LT. mn) THEN
      mn = res(i)
      mn_ind = i
    ENDIF
  ENDDO
  
  ! Find mx_ind, the index of b
  ! Two Cases:
  IF (x - mn .GE. 0) THEN     ! If x is above mn, then a=mn and b=v(mn_ind+1)
    IF (mn_ind .EQ. N) THEN   ! But if x is not within v at all, we will use
      mn_ind = N - 1          ! the last 2 points of v for our linear interpolation
      mx_ind = N
    ELSE
      mx_ind = mn_ind + 1
    ENDIF
  ELSE                        ! If x is below mn, then a=v(mn_ind-1) and b=mn
    IF (mn_ind .EQ. 1) THEN   ! But if x is not within v at all, we will use
      mx_ind = 2              ! the first 2 points of v for our linear interpolation
    ELSE
      mx_ind = mn_ind
      mn_ind = mn_ind - 1
    ENDIF
  ENDIF
  
  a = v(mn_ind)
  b = v(mx_ind)
  aindx = mn_ind
  bindx = mx_ind
  
  RETURN
END
