MODULE uvars

  ! USER VARIABLES
  ! The user of this program must set these variables

  IMPLICIT NONE
  SAVE
  
  ! Gride size
  INTEGER, PARAMETER :: Nvert = 15   ! Vertical 
  INTEGER, PARAMETER :: Nlon = 40    ! Longitude
  INTEGER, PARAMETER :: Nlat = 40    ! Latitude

  ! Uniform pressure at model top (hPa)
  REAL :: p_a_top = 250.0

  ! Variables for computing sigma
  REAL :: z_surf_test = 0.0     ! Average altitude at surface (km)
  
  REAL :: z_below = 10.0        ! Altitude (km) of empirical data available at next pressure above p_a_top
  REAL :: p_a_below = 265.0     ! Air pressure (hPa) corresponding to that altitude
  REAL :: g_below = 9.7764      ! As well as gravitational acceleration (m/s2),
  REAL :: rho_a_below = 0.414   ! and density (kg/m3)
  
  ! Grid paramters
  REAL :: dlon = 0.025   ! Longitude increment (deg)
  REAL :: dlat = 0.025   ! Latitude increment (deg)

  ! Initialize grid by setting coordinates of the SW corner (deg)
  REAL :: lat0 = 37.0    
  REAL :: lon0 = -119.0
  
  INTEGER :: topo = 0
  !CHARACTER(LEN=100) :: topofile = 'yosemite_dem_40.txt'
  CHARACTER(LEN=100) :: topofile = ''
  
  !INTEGER :: Nt = 7200
  !REAL :: h = 0.5
  !INTEGER :: Nt = 3600
  !REAL :: h = 1
  INTEGER :: Nt = 1800
  REAL :: h = 2

  ! Print Parameters
  INTEGER, PARAMETER :: npT = 8     ! Number of instances in time to be printed
  INTEGER, PARAMETER :: npL = 3     ! Number of sigma levels to be printed
  INTEGER :: hindex = 3             ! Version of h (time step)
  INTEGER, DIMENSION(npT) :: printTimes = (/ 2, 15, 30, 150, 300, 450, 900, 1800 /) ! Print at these times
  INTEGER, DIMENSION(npL) :: printLevels = (/ 1, 7, 13 /)                           ! Print these sigma levels


END MODULE uvars
