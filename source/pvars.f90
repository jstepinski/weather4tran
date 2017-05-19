MODULE pvars

  ! PROGRAM VARIABLES
  ! Either computed based on user inputs, or set independently of the user

  USE uvars

  IMPLICIT NONE
  SAVE
  
  ! The table of empirical data contains information about 78 altitudes
  INTEGER, PARAMETER :: Ntable = 78
  
  REAL, DIMENSION(Ntable) :: zdat = 0.   ! The altitudes (km)
  REAL, DIMENSION(Ntable) :: gdat = 0.   ! gravitational acceleration (m/s2)
  REAL, DIMENSION(Ntable) :: pdat = 0.   ! pressures (hPa)
  REAL, DIMENSION(Ntable) :: Tdat = 0.   ! temperature (K)
  REAL, DIMENSION(Ntable) :: rhodat = 0. ! density (kg/m3) 
  
  ! We will need to interpolate the empirical data for our model. The following
  ! variables store the results; their names correspond to the "dat" variables
  REAL :: g_int
  REAL :: rho_int

  ! Altitudes and Pressures 
  REAL, DIMENSION(Nvert+1) :: z = 0.            ! (km)
  REAL, DIMENSION(Nvert+1) :: p_a_test = 0.     ! (hPa)
  
  REAL :: z_top_test    ! Store result of equation (7.2) in km

  ! Sigma levels
  REAL, DIMENSION(Nvert+1) :: sigma = 0.
  
  ! Constants
  REAL :: Re1   = 6.371           ! Radius of Earth (millions of meters)
  REAL :: Re    = 6371*1000.0      ! Radius of Earth (meters)
  REAL :: pi    = 3.14159          ! pi
  REAL :: cpd   = 1004.67
  REAL :: kappa = 0.286
  REAL :: omega = 7.292*(10**(-5))

  ! Grid
  REAL, DIMENSION(Nlon) :: lon_mdpt   ! longitudes (rad)
  REAL, DIMENSION(Nlat) :: lat_mdpt   ! latitudes (rad)
  REAL, DIMENSION(Nlat+1) :: lat_edge

  ! Initial surface pressure p_a and corresponding column pressure pi_a, both in hPa
  
  REAL, DIMENSION(Nlat,Nlon) :: p_a_base = 1000.0       ! Background pressure (hPa)
  REAL :: dp_a_peak = 3.0                               ! Peak increment (hPa)
  
  REAL, DIMENSION(Nlat,Nlon) :: p_a_surf
  REAL, DIMENSION(Nlat,Nlon) :: pi_a_init

  ! Surface geopotential, used to read in the DEM
  REAL, DIMENSION(Nlon,Nlat) :: phi_surf_temp
  REAL, DIMENSION(Nlat,Nlon) :: phi_surf
  
  ! Dummy variable used only once as a temperature storage for temperature
  REAL :: Tinit
  
  ! Instantaneous properties of the system: flux and exener pressure
  ! These are entirely reevaluated with each Matsuno step, and they are
  ! used by all the transport equations
  
  REAL, DIMENSION(Nlat,Nlon+1,Nvert) :: fluxF
  REAL, DIMENSION(Nlat+1,Nlon,Nvert) :: fluxG

  REAL, DIMENSION(Nlat,Nlon,Nvert) :: PP_mdpt
  REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: PP_edge
  
END MODULE pvars
