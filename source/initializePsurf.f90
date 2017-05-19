SUBROUTINE initializePsurf

! Initialize the surface pressure to be a Gaussian distribution
! Intialize column pressures based on surface pressures and model top pressure

  USE uvars
  USE pvars

  IMPLICIT none

!******************************************************************************
! Local Variables

  INTEGER :: i,j
  REAL :: lon_c, lat_c   ! Grid centers

!******************************************************************************

  ! Compute grid center as average of grid edges. The center does not actually have
  ! to be a point on the grid; we need it to compute the pressure distribution
  lon_c = 0.5*(lon_mdpt(1) + lon_mdpt(Nlon))
  lat_c = 0.5*(lat_mdpt(1) + lat_mdpt(Nlat))

  ! Loop over the grid, and in each cell evaluate surface pressure p_a_surf (hPa) and
  ! column pressure pi_a (hPa). 
  ! Coordinates are stored in radians at this point, so they can be used as-is in COS,
  ! but they must be converted back to degrees for use in the Gaussian function
  DO i = 1,Nlat
    DO j = 1,Nlon
      p_a_surf(i,j) = p_a_base(i,j) + dp_a_peak* &
                      exp( -0.5*( Re1*cos(0.5*(lat_mdpt(i)+lat_c))*(lon_mdpt(j)-lon_c)*180.0/pi )**2 - &
                            0.5*(Re1*(lat_mdpt(i)-lat_c)*180.0/pi)**2)
      pi_a_init(i,j) = p_a_surf(i,j) - p_a_top 
    ENDDO
  ENDDO

  ! Print the surface pressures to a file as a matrix that can 
  ! easily be visualized to verify the computation
  OPEN( UNIT = 229, FILE = 'initP.txt', STATUS = 'UNKNOWN' )
  DO i = 1,Nlat
    DO j = 1,Nlon
      WRITE(229, '(f10.4,2X)', advance='no') p_a_surf(i,j)
    ENDDO
    WRITE(229, *) ''
  ENDDO
  CLOSE(229)

  RETURN
END

