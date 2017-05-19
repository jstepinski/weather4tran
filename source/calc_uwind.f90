      FUNCTION uwind_calc(uwind_prev, uwind_curr, vwind, pi_a_prev, pi_a_prev2, pi_a_curr, pi_a_est, & 
                           sigmaDot, phi, phi_prev, theta)
      
          USE pvars
          USE xfun

          IMPLICIT none
 
          REAL, DIMENSION(:,:,:), INTENT(IN) :: uwind_prev
          REAL, DIMENSION(:,:,:), INTENT(IN) :: uwind_curr
          REAL, DIMENSION(:,:,:), INTENT(IN) :: vwind
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_prev
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_prev2
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_curr
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_est
          REAL, DIMENSION(:,:,:), INTENT(IN) :: sigmaDot
          REAL, DIMENSION(:,:,:), INTENT(IN) :: phi
          REAL, DIMENSION(:,:,:), INTENT(IN) :: phi_prev
          REAL, DIMENSION(:,:,:), INTENT(IN) :: theta
         
          REAL, DIMENSION(Nlat,Nlon+1,Nvert) :: uwind_calc
          
          REAL, DIMENSION(Nlat,Nlon+1,Nvert) :: uwind
          REAL, DIMENSION(Nlat+2, Nlon+3) :: fluxF_ext
          REAL, DIMENSION(Nlat+3, Nlon+2) :: fluxG_ext
          REAL, DIMENSION(Nlat+2, Nlon+3) :: uwind_curr_ext
          
          REAL :: pp1, pp2, pp3, pp4, pp5, pp6
          REAL :: pc1, pc2, pc3, pc4, pc5, pc6
          REAL :: sdm1, sdm2, sdm3, sdm4, sdm5, sdm6
          REAL :: sdp1, sdp2, sdp3, sdp4, sdp5, sdp6
          REAL :: Re2_dlat_dlon
          REAL :: pi_a_dA_prev, pi_a_dA_curr
          REAL :: pi_a_dA_sigmaDot_minus
          REAL :: pi_a_dA_sigmaDot_plus
          REAL :: B1, B2, C1, C2, D1, D2, E1, E2
          REAL :: u1, u2
          REAL :: fcorio
          
          INTEGER :: i,j,k
          INTEGER :: ii, jj
          
          DO k = 1,Nvert
            DO i = 1,Nlon+1
              DO j = 1,Nlat
              
                !************************************************
                if (i == 1) then
                    if (j == Nlat) then
                        pp1 = pi_a_prev(j,i)*cos(lat_mdpt(j))
                        pc1 = pi_a_curr(j,i)*cos(lat_mdpt(j))
                        sdm1 = sigmaDot(j,i,k)
                        sdp1 = sigmaDot(j,i,k+1)
                    else
                        pp1 = pi_a_prev(j+1,i)*cos(lat_mdpt(j+1))
                        pc1 = pi_a_curr(j+1,i)*cos(lat_mdpt(j+1))
                        sdm1 = sigmaDot(j+1,i,k)
                        sdp1 = sigmaDot(j+1,i,k+1)
                    endif
                else
                    if (j == Nlat) then
                        pp1 = pi_a_prev(j,i-1)*cos(lat_mdpt(j))
                        pc1 = pi_a_curr(j,i-1)*cos(lat_mdpt(j))
                        sdm1 = sigmaDot(j,i-1,k)
                        sdp1 = sigmaDot(j,i-1,k+1)
                    else
                        pp1 = pi_a_prev(j+1,i-1)*cos(lat_mdpt(j+1))
                        pc1 = pi_a_curr(j+1,i-1)*cos(lat_mdpt(j+1))
                        sdm1 = sigmaDot(j+1,i-1,k)
                        sdp1 = sigmaDot(j+1,i-1,k+1)
                    endif
                endif
                !************************************************
                if (i == Nlon+1) then
                    if (j == Nlat) then
                        pp2 = pi_a_prev(j,i-1)*cos(lat_mdpt(j))
                        pc2 = pi_a_curr(j,i-1)*cos(lat_mdpt(j))
                        sdm2 = sigmaDot(j,i-1,k)
                        sdp2 = sigmaDot(j,i-1,k+1)
                    else
                        pp2 = pi_a_prev(j+1,i-1)*cos(lat_mdpt(j+1))
                        pc2 = pi_a_curr(j+1,i-1)*cos(lat_mdpt(j+1))
                        sdm2 = sigmaDot(j+1,i-1,k)
                        sdp2 = sigmaDot(j+1,i-1,k+1)
                    endif
                else
                    if (j == Nlat) then
                        pp2 = pi_a_prev(j,i)*cos(lat_mdpt(j))
                        pc2 = pi_a_curr(j,i)*cos(lat_mdpt(j))
                        sdm2 = sigmaDot(j,i,k)
                        sdp2 = sigmaDot(j,i,k+1)
                    else
                        pp2 = pi_a_prev(j+1,i)*cos(lat_mdpt(j+1))
                        pc2 = pi_a_curr(j+1,i)*cos(lat_mdpt(j+1))
                        sdm2 = sigmaDot(j+1,i,k)
                        sdp2 = sigmaDot(j+1,i,k+1)
                    endif
                endif
                !************************************************
                if (i == 1) then
                    pp3 = pi_a_prev(j,i)*cos(lat_mdpt(j))
                    pc3 = pi_a_curr(j,i)*cos(lat_mdpt(j))
                    sdm3 = sigmaDot(j,i,k)
                    sdp3 = sigmaDot(j,i,k+1)
                else
                    pp3 = pi_a_prev(j,i-1)*cos(lat_mdpt(j))
                    pc3 = pi_a_curr(j,i-1)*cos(lat_mdpt(j))
                    sdm3 = sigmaDot(j,i-1,k)
                    sdp3 = sigmaDot(j,i-1,k+1)
                endif
                !************************************************
                if (i == Nlon+1) then
                    pp4 = pi_a_prev(j,i-1)*cos(lat_mdpt(j))
                    pc4 = pi_a_curr(j,i-1)*cos(lat_mdpt(j))
                    sdm4 = sigmaDot(j,i-1,k)
                    sdp4 = sigmaDot(j,i-1,k+1)
                else
                    pp4 = pi_a_prev(j,i)*cos(lat_mdpt(j))
                    pc4 = pi_a_curr(j,i)*cos(lat_mdpt(j))
                    sdm4 = sigmaDot(j,i,k)
                    sdp4 = sigmaDot(j,i,k+1)
                endif
                !************************************************
                if (i == 1) then
                    if (j == 1) then
                        pp5 = pi_a_prev(j,i)*cos(lat_mdpt(j))
                        pc5 = pi_a_curr(j,i)*cos(lat_mdpt(j))
                        sdm5 = sigmaDot(j,i,k)
                        sdp5 = sigmaDot(j,i,k+1)
                    else
                        pp5 = pi_a_prev(j-1,i)*cos(lat_mdpt(j-1))
                        pc5 = pi_a_curr(j-1,i)*cos(lat_mdpt(j-1))
                        sdm5 = sigmaDot(j-1,i,k)
                        sdp5 = sigmaDot(j-1,i,k+1)
                    endif
                else
                    if (j == 1) then
                        pp5 = pi_a_prev(j,i-1)*cos(lat_mdpt(j))
                        pc5 = pi_a_curr(j,i-1)*cos(lat_mdpt(j))
                        sdm5 = sigmaDot(j,i-1,k)
                        sdp5 = sigmaDot(j,i-1,k+1)
                    else
                        pp5 = pi_a_prev(j-1,i-1)*cos(lat_mdpt(j-1))
                        pc5 = pi_a_curr(j-1,i-1)*cos(lat_mdpt(j-1))
                        sdm5 = sigmaDot(j-1,i-1,k)
                        sdp5 = sigmaDot(j-1,i-1,k+1)
                    endif
                endif
                !************************************************
                if (i == Nlon+1) then
                    if (j == 1) then
                        pp6 = pi_a_prev(j,i-1)*cos(lat_mdpt(j))
                        pc6 = pi_a_curr(j,i-1)*cos(lat_mdpt(j))
                        sdm6 = sigmaDot(j,i-1,k)
                        sdp6 = sigmaDot(j,i-1,k+1)
                    else
                        pp6 = pi_a_prev(j-1,i-1)*cos(lat_mdpt(j-1))
                        pc6 = pi_a_curr(j-1,i-1)*cos(lat_mdpt(j-1))
                        sdm6 = sigmaDot(j-1,i-1,k)
                        sdp6 = sigmaDot(j-1,i-1,k+1)
                    endif
                else
                    if (j == 1) then
                        pp6 = pi_a_prev(j,i)*cos(lat_mdpt(j))
                        pc6 = pi_a_curr(j,i)*cos(lat_mdpt(j))
                        sdm6 = sigmaDot(j,i,k)
                        sdp6 = sigmaDot(j,i,k+1)
                    else
                        pp6 = pi_a_prev(j-1,i)*cos(lat_mdpt(j-1))
                        pc6 = pi_a_curr(j-1,i)*cos(lat_mdpt(j-1))
                        sdm6 = sigmaDot(j-1,i,k)
                        sdp6 = sigmaDot(j-1,i,k+1)
                    endif
                endif
                !************************************************
              
                Re2_dlat_dlon = (Re**2.0)*dlat*dlon
            
                pi_a_dA_prev = (   pp1*Re2_dlat_dlon &
                                 + pp2*Re2_dlat_dlon &
                                 + pp3*Re2_dlat_dlon*2.0 &
                                 + pp4*Re2_dlat_dlon*2.0 &
                                 + pp5*Re2_dlat_dlon &
                                 + pp6*Re2_dlat_dlon  )/8.0
                                 
                pi_a_dA_curr = (   pc1*Re2_dlat_dlon &
                                 + pc2*Re2_dlat_dlon &
                                 + pc3*Re2_dlat_dlon*2.0 &
                                 + pc4*Re2_dlat_dlon*2.0 &
                                 + pc5*Re2_dlat_dlon &
                                 + pc6*Re2_dlat_dlon  )/8.0
                             
                pi_a_dA_sigmaDot_minus = (   pc1*Re2_dlat_dlon*sdm1 &
                                           + pc2*Re2_dlat_dlon*sdm2 &
                                           + pc3*Re2_dlat_dlon*sdm3*2.0 &
                                           + pc4*Re2_dlat_dlon*sdm4*2.0 &
                                           + pc5*Re2_dlat_dlon*sdm5 &
                                           + pc6*Re2_dlat_dlon*sdm6 )/8.0
                
                pi_a_dA_sigmaDot_plus = (    pc1*Re2_dlat_dlon*sdp1 &
                                           + pc2*Re2_dlat_dlon*sdp2 &
                                           + pc3*Re2_dlat_dlon*sdp3*2.0 &
                                           + pc4*Re2_dlat_dlon*sdp4*2.0 &
                                           + pc5*Re2_dlat_dlon*sdp5 &
                                           + pc6*Re2_dlat_dlon*sdp6 )/8.0
                                           
                fluxF_ext = extendMat(fluxF(:,:,k),Nlat,Nlon+1)
                fluxG_ext = extendMat(fluxG(:,:,k),Nlat+1,Nlon)

                jj = j+1
                ii = i+1
                
                B1 = ( fluxF_ext(jj-1,ii-1) + fluxF_ext(jj-1,ii) + 2.0*fluxF_ext(jj,ii-1) + 2.0*fluxF_ext(jj,ii) + fluxF_ext(jj+1,ii-1) + fluxF_ext(jj+1,ii) )/12.0
                B2 = ( fluxF_ext(jj-1,ii) + fluxF_ext(jj-1,ii+1) + 2.0*fluxF_ext(jj,ii) + 2.0*fluxF_ext(jj,ii+1) + fluxF_ext(jj+1,ii) + fluxF_ext(jj+1,ii+1) )/12.0
                
                C1 = ( fluxG_ext(jj-1,ii-1) + fluxG_ext(jj-1,ii) + 2.0*fluxG_ext(jj,ii-1) + 2.0*fluxG_ext(jj,ii) + fluxG_ext(jj+1,ii-1) + fluxG_ext(jj+1,ii) )/12.0
                C2 = ( fluxG_ext(jj,ii-1) + fluxG_ext(jj,ii) + 2.0*fluxG_ext(jj+1,ii-1) + 2.0*fluxG_ext(jj,ii) + fluxG_ext(jj+2,ii-1) + fluxG_ext(jj+2,ii) )/12.0
                
                D1 = ( fluxG_ext(jj-1,ii-1) + 2.0*fluxG_ext(jj,ii-1) + fluxG_ext(jj+1,ii-1) + fluxF_ext(jj-1,ii-1) + fluxF_ext(jj,ii-1) + fluxF_ext(jj-1,ii) + fluxF_ext(jj,ii) )/24.0
                D2 = ( fluxG_ext(jj,ii) + 2.0*fluxG_ext(jj+1,ii) + fluxG_ext(jj+2,ii) + fluxF_ext(jj,ii) + fluxF_ext(jj+1,ii) + fluxF_ext(jj,ii+1) + fluxF_ext(jj+1,ii+1) )/24.0
                
                E1 = ( fluxG_ext(jj-1,ii) + 2.0*fluxG_ext(jj,ii) + fluxG_ext(jj+1,ii) - fluxF_ext(jj-1,ii) - fluxF_ext(jj,ii) - fluxF_ext(jj-1,ii+1) - fluxF_ext(jj,ii+1) )/24.0
                E2 = ( fluxG_ext(jj,ii-1) + 2.0*fluxG_ext(jj+1,ii-1) + fluxG_ext(jj+2,ii-1) - fluxF_ext(jj,ii-1) - fluxF_ext(jj+1,ii-1) - fluxF_ext(jj,ii) - fluxF_ext(jj+1,ii) )/24.0
            
                if (k == Nvert) then
                    u1 = ( (abs(sigma(k+1)-sigma(k)))*uwind_curr(j,i,k-1) + (abs(sigma(k)-sigma(k-1)))*uwind_curr(j,i,k) ) &
                         / ( (abs(sigma(k)-sigma(k-1))) + (abs(sigma(k+1)-sigma(k))) )

                    u2 = 0.0
                elseif (k == 1) then
                    u1 = 0.0

                    u2 = ( (abs(sigma(k+2)-sigma(k+1)))*uwind_curr(j,i,k) + (abs(sigma(k+1)-sigma(k)))*uwind_curr(j,i,k+1) ) &
                         / ( (abs(sigma(k+1)-sigma(k))) + (abs(sigma(k+2)-sigma(k+1))) )
                else
                    u1 = ( (abs(sigma(k+1)-sigma(k)))*uwind_curr(j,i,k-1) + (abs(sigma(k)-sigma(k-1)))*uwind_curr(j,i,k) ) &
                         / ( (abs(sigma(k)-sigma(k-1))) + (abs(sigma(k+1)-sigma(k))) )

                    u2 = ( (abs(sigma(k+2)-sigma(k+1)))*uwind_curr(j,i,k) + (abs(sigma(k+1)-sigma(k)))*uwind_curr(j,i,k+1) ) &
                         / ( (abs(sigma(k+1)-sigma(k))) + (abs(sigma(k+2)-sigma(k+1))) )
                endif
                
                fcorio = 2.0*omega*sin(lat_mdpt(j))
            
            
                uwind_curr_ext = extendMat(uwind_curr(:,:,k),Nlat,Nlon+1)
                
                if (i == 1) then
                    uwind(j,i,k) =  pi_a_dA_prev*uwind_prev(j,i,k)/pi_a_dA_curr &
                                 +  (h/pi_a_dA_curr) * &
                              (     (   B1*(uwind_curr_ext(jj,ii-1) + uwind_curr_ext(jj,ii))/2.0 - B2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj,ii+1))/2.0 &
                                      + C1*(uwind_curr_ext(jj-1,ii) + uwind_curr_ext(jj,ii))/2.0 - C2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii))/2.0 &
                                      + D1*(uwind_curr_ext(jj-1,ii-1) + uwind_curr_ext(jj,ii))/2.0 - D2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii+1))/2.0 &
                                      + E1*(uwind_curr_ext(jj-1,ii+1) + uwind_curr_ext(jj,ii))/2.0 - E2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii-1))/2.0 &
                                    ) &
                                 +  (1.0/abs(sigma(k+1)-sigma(k)))*( pi_a_dA_sigmaDot_minus*u1 - pi_a_dA_sigmaDot_plus*u2 ) &
                                 +  (Re*dlat*dlon/2.0) &
                                 *  (   pi_a_est(j,i)*((vwind(j,i,k) + vwind(j+1,i,k) )/2.0)*(fcorio*Re*cos(lat_mdpt(j)) + (uwind_curr(j,i,k) + uwind_curr(j,i,k))*sin(lat_mdpt(j))/2.0) &
                                      + pi_a_est(j,i)*((vwind(j,i,k) + vwind(j+1,i,k) )/2.0)*(fcorio*Re*cos(lat_mdpt(j)) + (uwind_curr(j,i,k) + uwind_curr(j,i+1,k))*sin(lat_mdpt(j))/2.0) &
                                    ) &
                                 -  Re*dlat &
                                 *  (   (phi_prev(j,i,k) - phi(j,i,k))*pi_a_est(j,i) &
                                      + (pi_a_prev2(j,i) - pi_a_est(j,i))*(cpd/2.0) &
                                      * ( theta(j,i,k)*( sigma(k+1)*(PP_edge(j,i,k+1) - PP_mdpt(j,i,k)) + sigma(k)*(PP_mdpt(j,i,k) - PP_edge(j,i,k)) ) &
                                                         /abs(sigma(k+1)-sigma(k)) &
                                        ) &
                                     ) &
                              )  
                elseif (i == Nlon+1) then
                    uwind(j,i,k) =  pi_a_dA_prev*uwind_prev(j,i,k)/pi_a_dA_curr &
                                 +  (h/pi_a_dA_curr) * &
                              (     (   B1*(uwind_curr_ext(jj,ii-1) + uwind_curr_ext(jj,ii))/2.0 - B2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj,ii+1))/2.0 &
                                      + C1*(uwind_curr_ext(jj-1,ii) + uwind_curr_ext(jj,ii))/2.0 - C2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii))/2.0 &
                                      + D1*(uwind_curr_ext(jj-1,ii-1) + uwind_curr_ext(jj,ii))/2.0 - D2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii+1))/2.0 &
                                      + E1*(uwind_curr_ext(jj-1,ii+1) + uwind_curr_ext(jj,ii))/2.0 - E2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii-1))/2.0 &
                                    ) &
                                 +  (1.0/abs(sigma(k+1)-sigma(k)))*( pi_a_dA_sigmaDot_minus*u1 - pi_a_dA_sigmaDot_plus*u2 ) &
                                 +  (Re*dlat*dlon/2.0) &
                                 *  (   pi_a_est(j,i-1)*((vwind(j,i-1,k) + vwind(j+1,i-1,k) )/2.0)*(fcorio*Re*cos(lat_mdpt(j)) + (uwind_curr(j,i-1,k) + uwind_curr(j,i,k))*sin(lat_mdpt(j))/2.0) &
                                      + pi_a_est(j,i-1)*((vwind(j,i-1,k) + vwind(j+1,i-1,k) )/2.0)*(fcorio*Re*cos(lat_mdpt(j)) + (uwind_curr(j,i,k) + uwind_curr(j,i,k))*sin(lat_mdpt(j))/2.0) &
                                    ) &
                                 -  Re*dlat &
                                 *  (   (phi_prev(j,i-1,k) - phi(j,i-1,k))*pi_a_est(j,i-1) &
                                      + (pi_a_prev2(j,i-1) - pi_a_est(j,i-1))*(cpd/2.0) &
                                      * ( theta(j,i-1,k)*( sigma(k+1)*(PP_edge(j,i-1,k+1) - PP_mdpt(j,i-1,k)) + sigma(k)*(PP_mdpt(j,i-1,k) - PP_edge(j,i-1,k)) ) &
                                                         /abs(sigma(k+1)-sigma(k)) &
                                        ) &
                                     ) &
                              )
                else            
                    uwind(j,i,k) =  pi_a_dA_prev*uwind_prev(j,i,k)/pi_a_dA_curr &
                                 +  (h/pi_a_dA_curr) * &
                              (     (   B1*(uwind_curr_ext(jj,ii-1) + uwind_curr_ext(jj,ii))/2.0 - B2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj,ii+1))/2.0 &
                                      + C1*(uwind_curr_ext(jj-1,ii) + uwind_curr_ext(jj,ii))/2.0 - C2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii))/2.0 &
                                      + D1*(uwind_curr_ext(jj-1,ii-1) + uwind_curr_ext(jj,ii))/2.0 - D2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii+1))/2.0 &
                                      + E1*(uwind_curr_ext(jj-1,ii+1) + uwind_curr_ext(jj,ii))/2.0 - E2*(uwind_curr_ext(jj,ii) + uwind_curr_ext(jj+1,ii-1))/2.0 &
                                    ) &
                                 +  (1.0/abs(sigma(k+1)-sigma(k)))*( pi_a_dA_sigmaDot_minus*u1 - pi_a_dA_sigmaDot_plus*u2 ) &
                                 +  (Re*dlat*dlon/2.0) &
                                 *  (   pi_a_est(j,i-1)*((vwind(j,i-1,k) + vwind(j+1,i-1,k) )/2.0)*(fcorio*Re*cos(lat_mdpt(j)) + (uwind_curr(j,i-1,k) + uwind_curr(j,i,k))*sin(lat_mdpt(j))/2.0) &
                                      + pi_a_est(j,i)*((vwind(j,i,k) + vwind(j+1,i,k) )/2.0)*(fcorio*Re*cos(lat_mdpt(j)) + (uwind_curr(j,i,k) + uwind_curr(j,i+1,k))*sin(lat_mdpt(j))/2.0) &
                                    ) &
                                 -  Re*dlat &
                                 *  (   (phi(j,i,k) - phi(j,i-1,k))*(pi_a_est(j,i-1) + pi_a_est(j,i))/2.0 &
                                      + (pi_a_est(j,i) - pi_a_est(j,i-1))*(cpd/2.0) &
                                      * ( theta(j,i-1,k)*( sigma(k+1)*(PP_edge(j,i-1,k+1) - PP_mdpt(j,i-1,k)) + sigma(k)*(PP_mdpt(j,i-1,k) - PP_edge(j,i-1,k)) ) &
                                                         /abs(sigma(k+1)-sigma(k)) &
                                        + theta(j,i,k)*( sigma(k+1)*(PP_edge(j,i,k+1) - PP_mdpt(j,i,k)) + sigma(k)*(PP_mdpt(j,i,k) - PP_edge(j,i,k)) ) &
                                                         /abs(sigma(k+1)-sigma(k)) &
                                        ) &
                                     ) &
                              )
                endif
                
              ENDDO
            ENDDO
          ENDDO
          
        uwind_calc = uwind  
          
      END    
