PROGRAM MAIN

  USE uvars     ! import user variables, i.e. variables the user sets
  USE pvars     ! import program variables
  USE funcs
  
  IMPLICIT none
  
!*****************************************************************************! 
! Local Variables
 
  ! Naming convention: base property + suffix
  ! The suffixes describe either the grid over which the property is specified
  ! or the instance in time at which it is recorded.
 
  ! Pressures; not used in the model proper; only to initialize temperature
  REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: p_edge
  REAL, DIMENSION(Nlat,Nlon,Nvert) ::   p_mdpt
  REAL, DIMENSION(Nlat,Nlon) :: Pcurr
  REAL, DIMENSION(Nlat,Nlon) :: Pprev
  REAL, DIMENSION(Nlat,Nlon) :: PP

  ! Column pressures
  REAL, DIMENSION(Nlat,Nlon) :: pi_a_prev
  REAL, DIMENSION(Nlat,Nlon) :: pi_a_prev2
  REAL, DIMENSION(Nlat,Nlon) :: pi_a_est
  REAL, DIMENSION(Nlat,Nlon) :: pi_a_curr

  ! Potential virtual temperatures
  REAL, DIMENSION(Nlat,Nlon,Nvert) :: theta_prev
  REAL, DIMENSION(Nlat,Nlon,Nvert) :: theta_est
  REAL, DIMENSION(Nlat,Nlon,Nvert) :: theta_curr
  REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: theta_interp

  ! Wind velocities
  REAL, DIMENSION(Nlat,Nlon+1,Nvert) :: uwind_prev = 0.
  REAL, DIMENSION(Nlat,Nlon+1,Nvert) :: uwind_est
  REAL, DIMENSION(Nlat,Nlon+1,Nvert) :: uwind_curr

  REAL, DIMENSION(Nlat+1,Nlon,Nvert) :: vwind_prev = 0.
  REAL, DIMENSION(Nlat+1,Nlon,Nvert) :: vwind_est
  REAL, DIMENSION(Nlat+1,Nlon,Nvert) :: vwind_curr

  REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: sigmaDot = 0.

  ! Geopotential
  REAL, DIMENSION(Nlat,Nlon,Nvert) :: phi_mdpt
  REAL, DIMENSION(Nlat,Nlon,Nvert) :: phi_prev
 
  INTEGER :: i,j,k,t
  INTEGER :: printCount = 1

!*****************************************************************************!  
! Read empirical pressure, density, and gravity data from a file

  CALL readTableB1

!*****************************************************************************!
! Compute sigma
  
  ! Estimate the model top altitude with equation (7.2)
  ! Note the unit conversion. p_a is given in hPa, so we convert it to Pa
  !      the result is then in m, but z is in km, so we convert to km
  z_top_test = z_below + (100.0/1000.0)*(p_a_below - p_a_top)/(rho_a_below*g_below)

  ! Set initial values for z and p_a
  z(1) = z_top_test
  p_a_test(1) = p_a_top
  
  ! Iteratively compute each z_k with equation (7.3)
  ! Note that because k and Nvert are integers, we must cast them to reals to maintain
  ! numerical accuracy. 
  ! Compute the pressure at each lower altitude (k+1) by rearranging equation (2.41)
  ! Obtain rho and g at altitude (k+1) by intepolating the empirical data
  DO k = 1,Nvert
    z(k+1)  = z_surf_test + (z_top_test - z_surf_test)*(1.0 - real(k)/real(Nvert))
    CALL interp1d(zdat,rhodat,z(k+1),rho_int,Ntable)
    CALL interp1d(zdat,gdat,z(k+1),g_int,Ntable)
    p_a_test(k+1) = p_a_test(k) + rho_int*g_int*(z(k) - z(k+1))
  ENDDO

  ! Check that the lowest layer is at least 150m thick
  IF (z(Nvert) - z(Nvert+1) .LT. 0.150) THEN
    WRITE(*,*) 'Error: bottom layer is less than 150m thick'
    STOP 1
  ENDIF

  ! Compute the sigma values with equation (7.4)
  DO k = 2,Nvert+1
    sigma(k) = (p_a_test(k) - p_a_top)/(p_a_test(Nvert+1) - p_a_top)
  ENDDO
  
  ! Print sigma values (uncomment the following block to activate)
  !WRITE(*,*) 'sigma'
  !DO k = 1,Nvert+1
  !  WRITE(*,*) sigma(k)
  !ENDDO

!*****************************************************************************!
! Setting up the GRID

  ! Convert all coordinates to radians for the remainder of the program
  dlat = dlat*pi/180.0
  dlon = dlon*pi/180.0

  lat0 = lat0*pi/180.0
  lon0 = lon0*pi/180.0

  ! Compute latitudes at cell edges and midpoints
  lat_mdpt(1) = lat0*pi/180
  lat_edge(1) = (lat0 - dlat/2.0)*pi/180
  DO i = 2,Nlat
    lat_mdpt(i) = lat_mdpt(i-1) + dlat
    lat_edge(i) = lat_edge(i-1) + dlat
  ENDDO
  lat_edge(Nlat+1) = lat_edge(Nlat) + dlat

  ! Compute longitudes only at cell midpoints; they are never used at edges
  lon_mdpt(1) = lon0
  DO i = 2,Nlon
    lon_mdpt(i) = lon_mdpt(i-1) + dlon
  ENDDO
  
!*****************************************************************************!
! TOPOGRAPHY
! Load the input topography here, as it is necessary to adjust the initial pressure  

  ! If the user does not active the "topo" switch, the terrain is flat
  ! Otherwise a DEM is loaded. Fortran reads matrices by column, so the DEM
  ! is initially read into an array of flipped dimensions, i.e. Nlon by Nlat
  ! This temporary array is then transposed into the actual array of elevations
  IF (topo .eq. 0) THEN
    phi_surf = 0.
  ELSE
    OPEN(UNIT=62, FILE=topofile, STATUS='old')
    READ(62,*) phi_surf_temp
    CLOSE(62)
    
    phi_surf = transpose(phi_surf_temp)
    
    ! phi is geopotential, but at this point it contains the elevations
    ! Use the elevations to adjust pressure at ground level. Assume that 
    ! background pressure is 1000 hPa at sea level and falls by 100 hPa
    ! over every km of elevation
    DO i = 1,Nlat
      DO j = 1,Nlon
        p_a_base(i,j) = p_a_base(i,j) - (phi_surf(i,j)/1000)*100
      ENDDO
    ENDDO
    
    ! Convert elevation to geopotential
    phi_surf = phi_surf*9.8
    
  ENDIF
  
!*****************************************************************************!
! Initialize surface pressure

  CALL initializePsurf  

!*****************************************************************************!  
! Initialize temperature throughout the entire grid
! General algorithm:
!   1) compute absolute pressure p at vertical edges 
!      from the definition of column pressure and the sigma values1
!   2) compute exener pressure at vertical edges with Eq 7.10
!   3) interpolate between exener pressures at two edges with Eq 7.11
!      to get pressure at midpoint
!   4) invert Eq 7.10 to recover pressure at midpoint
!   5) interpolate empirical weather data (table B1) to get 
!      an initial temperature estimate at every cell center
!   6) convert temperature to potential virtual temperature (Eq 2.99)
 
  DO k = 1,Nvert+1
    p_edge(:,:,k) = p_a_top + sigma(k)*pi_a_init
  ENDDO

  PP_edge(:,:,1) = (p_a_top/1000.0)**kappa

  DO k = 1,Nvert
    PP_edge(:,:,k+1) = (p_edge(:,:,k+1)/1000.0)**kappa
    Pcurr = PP_edge(:,:,k+1)
    Pprev = PP_edge(:,:,k)
    PP_mdpt(:,:,k) = (1.0/(1.0+kappa))*( (Pcurr*p_edge(:,:,k+1) - Pprev*p_edge(:,:,k)) & 
                                          /(p_edge(:,:,k+1) - p_edge(:,:,k)) )
    PP = PP_mdpt(:,:,k)
    p_mdpt(:,:,k) = 1000.0*PP**(1.0/kappa)
  ENDDO

  DO k = 1,Nvert
    CALL interp1d(pdat,Tdat,p_mdpt(1,1,k),Tinit,Ntable)
    theta_prev(:,:,k) = Tinit/PP_mdpt(:,:,k)
  ENDDO
  
!*****************************************************************************! 
 ! Final initializations
 
 ! Set 2 previous column pressures to initial value; convert all units to Pa
  pi_a_prev  = pi_a_init*100.0
  pi_a_prev2 = pi_a_init*100.0
  p_a_top = p_a_top*100.0
 
!*****************************************************************************!  
! Matsuno Scheme
  
  DO t = 1,Nt
  
    ! Evaluate fluxes using most recent values
    CALL fluxF_calc(pi_a_prev, uwind_prev)
    CALL fluxG_calc(pi_a_prev, vwind_prev)
    
    ! Estimate column pressure from previous value
    pi_a_est = pi_a_calc(pi_a_prev)
    
    ! Estimate vertical velocity from previous values
    sigmaDot = sigmaDot_calc(pi_a_prev, pi_a_est, Nvert+1)
   
    ! Evaluate Exener pressures at the newly estimate column pressure
    CALL PP_calc(pi_a_prev)
   
    ! Interpolate previous temperature onto vertical edges
    theta_interp = theta_interp_calc(theta_prev)
    
    ! Estimate temperature using previous pi in the leading term and estimated in f()
    theta_est = theta_calc(pi_a_prev, pi_a_est, theta_prev, theta_prev, theta_interp, sigmaDot)
    
    ! Compute new geopotential. If the iterations have just begun, initialize the "previous" 
    ! geopotential, i.e. the one at (t-2h) needed for stable boundary conditions in the momentum eq,
    ! to the same value
    phi_mdpt = phi_calc(theta_prev)
    IF (t .eq. 1) THEN
      phi_prev = phi_mdpt
    ENDIF
    
    ! Estimate wind velocities
    uwind_est = uwind_calc(uwind_prev, uwind_prev, vwind_prev, pi_a_prev, pi_a_prev2, pi_a_est, pi_a_est, & 
                           sigmaDot, phi_mdpt, phi_prev, theta_prev)
                        
    vwind_est = vwind_calc(vwind_prev, vwind_prev, uwind_prev, pi_a_prev, pi_a_prev2, pi_a_est, pi_a_est, & 
                           sigmaDot, phi_mdpt, phi_prev, theta_prev)
                           
    !************************
    ! CORRECTION STEP 
    
    ! Update fluxes
    CALL fluxF_calc(pi_a_est, uwind_est)
    CALL fluxG_calc(pi_a_est, vwind_est)
    
    ! Update pressure using "prev," not "est," because "est" is already
    ! incorporated in the fluxes
    pi_a_curr = pi_a_calc(pi_a_prev)
    
    ! Same idea as before. Use original pi from the start of the 
    ! time step in the leading term, but use the newly computed, "current" pi
    ! in the finite difference equation, called f() in the general Matsuno scheme
    sigmaDot = sigmaDot_calc(pi_a_prev, pi_a_curr, Nvert+1)
 
    ! Update exener pressures
    CALL PP_calc(pi_a_est)
    
    ! Interpolate temperature as before
    theta_interp = theta_interp_calc(theta_est)
    
    ! Again, the leading terms are "prev," but those that go in f() are the newest ones
    theta_curr = theta_calc(pi_a_prev, pi_a_curr, theta_prev, theta_est, theta_interp, sigmaDot)
    
    ! Update geopotential
    phi_prev = phi_mdpt
    phi_mdpt = phi_calc(theta_est)
    
    ! Compute wind velocities
    uwind_curr = uwind_calc(uwind_prev, uwind_est, vwind_est, pi_a_prev, pi_a_prev, pi_a_curr, pi_a_est, & 
                            sigmaDot, phi_mdpt, phi_prev, theta_est)
                                                
    vwind_curr = vwind_calc(vwind_prev, vwind_est, uwind_est, pi_a_prev, pi_a_prev, pi_a_curr, pi_a_est, &
                            sigmaDot, phi_mdpt, phi_prev, theta_est)
    
    ! Prints time step to terminal so that you can track progress
    print*, t

    ! Print outputs to file
    IF (t .eq. printTimes(printCount)) THEN
      CALL printToFile(uwind_curr,vwind_curr,theta_curr,pi_a_curr,sigmaDot,printCount)
      printCount = printCount + 1
    ENDIF
    
    ! Update the arrays. Current values become previous ones, and previous ones become prev2's
    pi_a_prev2 = pi_a_prev
    phi_prev = phi_mdpt
    
    uwind_prev = uwind_curr
    vwind_prev = vwind_curr
    theta_prev = theta_curr
    pi_a_prev = pi_a_curr
    
  ENDDO

END






