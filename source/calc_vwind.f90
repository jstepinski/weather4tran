      FUNCTION vwind_calc(vwind_prev, vwind_curr, uwind, pi_a_prev, pi_a_prev2, pi_a_curr, pi_a_est, & 
                           sigmaDot, phi, phi_prev, theta)
               
          USE pvars
          USE xfun

          IMPLICIT none
          
          REAL, DIMENSION(:,:,:), INTENT(IN) :: vwind_prev
          REAL, DIMENSION(:,:,:), INTENT(IN) :: vwind_curr
          REAL, DIMENSION(:,:,:), INTENT(IN) :: uwind
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_prev
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_prev2
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_curr
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_est
          REAL, DIMENSION(:,:,:), INTENT(IN) :: sigmaDot
          REAL, DIMENSION(:,:,:), INTENT(IN) :: phi
          REAL, DIMENSION(:,:,:), INTENT(IN) :: phi_prev
          REAL, DIMENSION(:,:,:), INTENT(IN) :: theta

          REAL, DIMENSION(Nlat+1,Nlon,Nvert) :: vwind_calc
          
          REAL, DIMENSION(Nlat+1,Nlon,Nvert) :: vwind
          REAL, DIMENSION(Nlat+2, Nlon+3) :: fluxF_ext
          REAL, DIMENSION(Nlat+3, Nlon+2) :: fluxG_ext
          REAL, DIMENSION(Nlat+3, Nlon+2) :: vwind_curr_ext
          
          REAL :: pp1, pp2, pp3, pp4, pp5, pp6
          REAL :: pc1, pc2, pc3, pc4, pc5, pc6
          REAL :: sdm1, sdm2, sdm3, sdm4, sdm5, sdm6
          REAL :: sdp1, sdp2, sdp3, sdp4, sdp5, sdp6
          REAL :: Re2_dlat_dlon
          REAL :: pi_a_dA_prev, pi_a_dA_curr
          REAL :: pi_a_dA_sigmaDot_minus
          REAL :: pi_a_dA_sigmaDot_plus
          REAL :: Q1, Q2, R1, R2, S1, S2, T1, T2
          REAL :: v1, v2
          REAL :: fcorio1, fcorio2
          REAL :: lat1, lat2
          
          INTEGER :: i,j,k
          INTEGER :: ii, jj
          
          DO k = 1,Nvert
            DO i = 1,Nlon
              DO j = 1,Nlat+1
              
                !************************************************
                if (i == Nlon) then
                    if (j == 1) then
                        pp1 = pi_a_prev(j,i)*cos(lat_edge(j))
                        pc1 = pi_a_curr(j,i)*cos(lat_edge(j))
                        sdm1 = sigmaDot(j,i,k)
                        sdp1 = sigmaDot(j,i,k+1)
                    else
                        pp1 = pi_a_prev(j-1,i)*cos(lat_edge(j-1))
                        pc1 = pi_a_curr(j-1,i)*cos(lat_edge(j-1))
                        sdm1 = sigmaDot(j-1,i,k)
                        sdp1 = sigmaDot(j-1,i,k+1)
                    endif
                else
                    if (j == 1) then
                        pp1 = pi_a_prev(j,i+1)*cos(lat_edge(j))
                        pc1 = pi_a_curr(j,i+1)*cos(lat_edge(j))
                        sdm1 = sigmaDot(j,i+1,k)
                        sdp1 = sigmaDot(j,i+1,k+1)
                    else
                        pp1 = pi_a_prev(j-1,i+1)*cos(lat_edge(j-1))
                        pc1 = pi_a_curr(j-1,i+1)*cos(lat_edge(j-1))
                        sdm1 = sigmaDot(j-1,i+1,k)
                        sdp1 = sigmaDot(j-1,i+1,k+1)
                    endif
                endif
                !************************************************
                if (i == Nlon) then
                    if (j == Nlat+1) then
                        pp2 = pi_a_prev(j-1,i)*cos(lat_edge(j-1))
                        pc2 = pi_a_curr(j-1,i)*cos(lat_edge(j-1))
                        sdm2 = sigmaDot(j-1,i,k)
                        sdp2 = sigmaDot(j-1,i,k+1)
                    else
                        pp2 = pi_a_prev(j,i)*cos(lat_edge(j))
                        pc2 = pi_a_curr(j,i)*cos(lat_edge(j))
                        sdm2 = sigmaDot(j,i,k)
                        sdp2 = sigmaDot(j,i,k+1)
                    endif
                else
                    if (j == Nlat+1) then
                        pp2 = pi_a_prev(j-1,i+1)*cos(lat_edge(j-1))
                        pc2 = pi_a_curr(j-1,i+1)*cos(lat_edge(j-1))
                        sdm2 = sigmaDot(j-1,i+1,k)
                        sdp2 = sigmaDot(j-1,i+1,k+1)
                    else
                        pp2 = pi_a_prev(j,i+1)*cos(lat_edge(j))
                        pc2 = pi_a_curr(j,i+1)*cos(lat_edge(j))
                        sdm2 = sigmaDot(j,i+1,k)
                        sdp2 = sigmaDot(j,i+1,k+1)
                    endif
                endif
                !************************************************
                if (j == 1) then
                    pp3 = pi_a_prev(j,i)*cos(lat_edge(j))
                    pc3 = pi_a_curr(j,i)*cos(lat_edge(j))
                    sdm3 = sigmaDot(j,i,k)
                    sdp3 = sigmaDot(j,i,k+1)
                else
                    pp3 = pi_a_prev(j-1,i)*cos(lat_edge(j-1))
                    pc3 = pi_a_curr(j-1,i)*cos(lat_edge(j-1))
                    sdm3 = sigmaDot(j-1,i,k)
                    sdp3 = sigmaDot(j-1,i,k+1)
                endif
                !************************************************
                if (j == Nlat+1) then
                    pp4 = pi_a_prev(j-1,i)*cos(lat_edge(j-1))
                    pc4 = pi_a_curr(j-1,i)*cos(lat_edge(j-1))
                    sdm4 = sigmaDot(j-1,i,k)
                    sdp4 = sigmaDot(j-1,i,k+1)
                else
                    pp4 = pi_a_prev(j,i)*cos(lat_edge(j))
                    pc4 = pi_a_curr(j,i)*cos(lat_edge(j))
                    sdm4 = sigmaDot(j,i,k)
                    sdp4 = sigmaDot(j,i,k+1)
                endif
                !************************************************
                if (i == 1) then
                    if (j == 1) then
                        pp5 = pi_a_prev(j,i)*cos(lat_edge(j))
                        pc5 = pi_a_curr(j,i)*cos(lat_edge(j))
                        sdm5 = sigmaDot(j,i,k)
                        sdp5 = sigmaDot(j,i,k+1)
                    else
                        pp5 = pi_a_prev(j-1,i)*cos(lat_edge(j-1))
                        pc5 = pi_a_curr(j-1,i)*cos(lat_edge(j-1))
                        sdm5 = sigmaDot(j-1,i,k)
                        sdp5 = sigmaDot(j-1,i,k+1)
                    endif
                else
                    if (j == 1) then
                        pp5 = pi_a_prev(j,i-1)*cos(lat_edge(j))
                        pc5 = pi_a_curr(j,i-1)*cos(lat_edge(j))
                        sdm5 = sigmaDot(j,i-1,k)
                        sdp5 = sigmaDot(j,i-1,k+1)
                    else
                        pp5 = pi_a_prev(j-1,i-1)*cos(lat_edge(j-1))
                        pc5 = pi_a_curr(j-1,i-1)*cos(lat_edge(j-1))
                        sdm5 = sigmaDot(j-1,i-1,k)
                        sdp5 = sigmaDot(j-1,i-1,k+1)
                    endif
                endif
                !************************************************
                if (i == 1) then
                    if (j == Nlat+1) then
                        pp6 = pi_a_prev(j-1,i)*cos(lat_edge(j-1))
                        pc6 = pi_a_curr(j-1,i)*cos(lat_edge(j-1))
                        sdm6 = sigmaDot(j-1,i,k)
                        sdp6 = sigmaDot(j-1,i,k+1)
                    else
                        pp6 = pi_a_prev(j,i)*cos(lat_edge(j))
                        pc6 = pi_a_curr(j,i)*cos(lat_edge(j))
                        sdm6 = sigmaDot(j,i,k)
                        sdp6 = sigmaDot(j,i,k+1)
                    endif
                else
                    if (j == Nlat+1) then
                        pp6 = pi_a_prev(j-1,i-1)*cos(lat_edge(j-1))
                        pc6 = pi_a_curr(j-1,i-1)*cos(lat_edge(j-1))
                        sdm6 = sigmaDot(j-1,i-1,k)
                        sdp6 = sigmaDot(j-1,i-1,k+1)
                    else
                        pp6 = pi_a_prev(j,i-1)*cos(lat_edge(j))
                        pc6 = pi_a_curr(j,i-1)*cos(lat_edge(j))
                        sdm6 = sigmaDot(j,i-1,k)
                        sdp6 = sigmaDot(j,i-1,k+1)
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
                
                Q1 = ( fluxF_ext(jj-1,ii-1) + fluxF_ext(jj,ii-1) + 2.0*fluxF_ext(jj-1,ii) + 2.0*fluxF_ext(jj,ii) + fluxF_ext(jj-1,ii+1) + fluxF_ext(jj,ii+1) )/12.0
                Q2 = ( fluxF_ext(jj-1,ii) + fluxF_ext(jj,ii) + 2.0*fluxF_ext(jj-1,ii+1) + 2.0*fluxF_ext(jj,ii+1) + fluxF_ext(jj-1,ii+2) + fluxF_ext(jj,ii+2) )/12.0
                
                R1 = ( fluxG_ext(jj-1,ii-1) + fluxG_ext(jj,ii-1) + 2.0*fluxG_ext(jj-1,ii) + 2.0*fluxG_ext(jj,ii) + fluxG_ext(jj-1,ii+1) + fluxG_ext(jj,ii+1) )/12.0
                R2 = ( fluxG_ext(jj,ii-1) + fluxG_ext(jj+1,ii-1) + 2.0*fluxG_ext(jj,ii) + 2.0*fluxG_ext(jj+1,ii) + fluxG_ext(jj,ii+1) + fluxG_ext(jj+1,ii+1) )/12.0
                
                S1 = ( fluxG_ext(jj-1,ii-1) + fluxG_ext(jj,ii-1) + fluxG_ext(jj-1,ii) + fluxG_ext(jj,ii) + fluxF_ext(jj-1,ii-1) + 2.0*fluxF_ext(jj-1,ii) + fluxF_ext(jj-1,ii+1) )/24.0
                S2 = ( fluxG_ext(jj,ii) + fluxG_ext(jj+1,ii) + fluxG_ext(jj,ii+1) + fluxG_ext(jj+1,ii+1) + fluxF_ext(jj,ii) + 2.0*fluxF_ext(jj,ii+1) + fluxF_ext(jj,ii+2) )/24.0
                
                T1 = ( fluxG_ext(jj-1,ii) + fluxG_ext(jj,ii) + fluxG_ext(jj-1,ii+1) + fluxG_ext(jj,ii+1) - fluxF_ext(jj-1,ii) - 2.0*fluxF_ext(jj-1,ii+1) - fluxF_ext(jj-1,ii+2) )/24.0
                T2 = ( fluxG_ext(jj,ii-1) + fluxG_ext(jj+1,ii-1) + fluxG_ext(jj,ii) + fluxG_ext(jj+1,ii) - fluxF_ext(jj,ii-1) - 2.0*fluxF_ext(jj,ii) - fluxF_ext(jj,ii+1) )/24.0
			
                if (k == Nvert) then
                    v1 = ( (abs(sigma(k+1)-sigma(k)))*vwind_curr(j,i,k-1) + (abs(sigma(k)-sigma(k-1)))*vwind_curr(j,i,k) ) &
                         / ( (abs(sigma(k)-sigma(k-1))) + (abs(sigma(k+1)-sigma(k))) )

                    v2 = 0.0
                elseif (k == 1) then
                    v1 = 0.0

                    v2 = ( (abs(sigma(k+2)-sigma(k+1)))*vwind_curr(j,i,k) + (abs(sigma(k+1)-sigma(k)))*vwind_curr(j,i,k+1) ) &
                         / ( (abs(sigma(k+1)-sigma(k))) + (abs(sigma(k+2)-sigma(k+1))) )
                else
                    v1 = ( (abs(sigma(k+1)-sigma(k)))*vwind_curr(j,i,k-1) + (abs(sigma(k)-sigma(k-1)))*vwind_curr(j,i,k) ) &
                         / ( (abs(sigma(k)-sigma(k-1))) + (abs(sigma(k+1)-sigma(k))) )
             
                    v2 = ( (abs(sigma(k+2)-sigma(k+1)))*vwind_curr(j,i,k) + (abs(sigma(k+1)-sigma(k)))*vwind_curr(j,i,k+1) ) &
                         / ( (abs(sigma(k+1)-sigma(k))) + (abs(sigma(k+2)-sigma(k+1))) )
                endif
                
                vwind_curr_ext = extendMat(vwind_curr(:,:,k),Nlat+1,Nlon)
                
                if (j == 1) then
                
                    lat1 = (lat_edge(j) + lat_edge(j+1))/2.0
                    lat2 = (lat_edge(j) + lat_edge(j+1))/2.0
                    
                    fcorio1 = 2.0*omega*sin(lat1); 
                    fcorio2 = 2.0*omega*sin(lat2);
                
                    vwind(j,i,k) =  pi_a_dA_prev*vwind_prev(j,i,k)/pi_a_dA_curr &
                                 +  (h/pi_a_dA_curr) * &
                              (     (   Q1*(vwind_curr_ext(jj,ii-1) + vwind_curr_ext(jj,ii))/2.0 - Q2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj,ii+1))/2.0 &
                                      + R1*(vwind_curr_ext(jj-1,ii) + vwind_curr_ext(jj,ii))/2.0 - R2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii))/2.0 &
                                      + S1*(vwind_curr_ext(jj-1,ii-1) + vwind_curr_ext(jj,ii))/2.0 - S2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii+1))/2.0 &
                                      + T1*(vwind_curr_ext(jj-1,ii+1) + vwind_curr_ext(jj,ii))/2.0 - T2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii-1))/2.0 &
                                    ) &
                                 +  (1.0/abs(sigma(k+1)-sigma(k)))*( pi_a_dA_sigmaDot_minus*v1 - pi_a_dA_sigmaDot_plus*v2 ) &
                                 -  (Re*dlat*dlon/2.0) &
                                 *  (   pi_a_est(j,i)*((uwind(j,i,k) + uwind(j,i+1,k) )/2.0)*(fcorio1*Re*cos(lat1) + (uwind(j,i,k) + uwind(j,i+1,k))*sin(lat1)/2.0) &
                                      + pi_a_est(j,i)*((uwind(j,i,k) + uwind(j,i+1,k) )/2.0)*(fcorio2*Re*cos(lat2) + (uwind(j,i,k) + uwind(j,i+1,k))*sin(lat2)/2.0) &
                                    ) &
                                 -  Re*dlon*cos(lat_edge(j)) &
                                 *  (   (phi_prev(j,i,k) - phi(j,i,k))*pi_a_est(j,i) &
                                      + (pi_a_prev2(j,i) - pi_a_est(j,i))*(cpd/2.0) &
                                      * ( theta(j,i,k)*( sigma(k+1)*(PP_edge(j,i,k+1) - PP_mdpt(j,i,k)) + sigma(k)*(PP_mdpt(j,i,k) - PP_edge(j,i,k)) ) &
                                                         /abs(sigma(k+1)-sigma(k)) &
                                        ) &
                                     ) &
                              )  
                elseif (j == Nlat+1) then
                
                    lat1 = (lat_edge(j-1) + lat_edge(j))/2.0
                    lat2 = (lat_edge(j-1) + lat_edge(j))/2.0
                    
                    fcorio1 = 2.0*omega*sin(lat1)
                    fcorio2 = 2.0*omega*sin(lat2)
                
                    vwind(j,i,k) =  pi_a_dA_prev*vwind_prev(j,i,k)/pi_a_dA_curr &
                                 +  (h/pi_a_dA_curr) * &
                              (     (   Q1*(vwind_curr_ext(jj,ii-1) + vwind_curr_ext(jj,ii))/2.0 - Q2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj,ii+1))/2.0 &
                                      + R1*(vwind_curr_ext(jj-1,ii) + vwind_curr_ext(jj,ii))/2.0 - R2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii))/2.0 &
                                      + S1*(vwind_curr_ext(jj-1,ii-1) + vwind_curr_ext(jj,ii))/2.0 - S2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii+1))/2.0 &
                                      + T1*(vwind_curr_ext(jj-1,ii+1) + vwind_curr_ext(jj,ii))/2.0 - T2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii-1))/2.0 &
                                    ) &
                                 +  (1.0/abs(sigma(k+1)-sigma(k)))*( pi_a_dA_sigmaDot_minus*v1 - pi_a_dA_sigmaDot_plus*v2 ) &
                                 -  (Re*dlat*dlon/2.0) &
                                 *  (   pi_a_est(j-1,i)*((uwind(j-1,i,k) + uwind(j-1,i+1,k) )/2.0)*(fcorio1*Re*cos(lat1) + (uwind(j-1,i,k) + uwind(j-1,i+1,k))*sin(lat1)/2.0) &
                                      + pi_a_est(j-1,i)*((uwind(j-1,i,k) + uwind(j-1,i+1,k) )/2.0)*(fcorio2*Re*cos(lat2) + (uwind(j-1,i,k) + uwind(j-1,i+1,k))*sin(lat2)/2.0) &
                                    ) &
                                 -  Re*dlon*cos(lat_edge(j)) &
                                 *  (   (phi_prev(j-1,i,k) - phi(j-1,i,k))*pi_a_est(j-1,i) &
                                      + (pi_a_prev2(j-1,i) - pi_a_est(j-1,i))*(cpd/2.0) &
                                      * ( theta(j-1,i,k)*( sigma(k+1)*(PP_edge(j-1,i,k+1) - PP_mdpt(j-1,i,k)) + sigma(k)*(PP_mdpt(j-1,i,k) - PP_edge(j-1,i,k)) ) &
                                                         /abs(sigma(k+1)-sigma(k)) &
                                        ) &
                                     ) &
                              )
                else    

                    lat1 = (lat_edge(j-1) + lat_edge(j))/2.0
                    lat2 = (lat_edge(j) + lat_edge(j+1))/2.0
                    
                    fcorio1 = 2.0*omega*sin(lat1)
                    fcorio2 = 2.0*omega*sin(lat2)
                
                    vwind(j,i,k) =  pi_a_dA_prev*vwind_prev(j,i,k)/pi_a_dA_curr &
                                 +  (h/pi_a_dA_curr) *  &
                              (     (   Q1*(vwind_curr_ext(jj,ii-1) + vwind_curr_ext(jj,ii))/2.0 - Q2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj,ii+1))/2.0 &
                                      + R1*(vwind_curr_ext(jj-1,ii) + vwind_curr_ext(jj,ii))/2.0 - R2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii))/2.0 &
                                      + S1*(vwind_curr_ext(jj-1,ii-1) + vwind_curr_ext(jj,ii))/2.0 - S2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii+1))/2.0 &
                                      + T1*(vwind_curr_ext(jj-1,ii+1) + vwind_curr_ext(jj,ii))/2.0 - T2*(vwind_curr_ext(jj,ii) + vwind_curr_ext(jj+1,ii-1))/2.0 &
                                    ) &
                                 +  (1.0/abs(sigma(k+1)-sigma(k)))*( pi_a_dA_sigmaDot_minus*v1 - pi_a_dA_sigmaDot_plus*v2 ) &
                                 -  (Re*dlat*dlon/2.0) &
                                 *  (   pi_a_est(j-1,i)*((uwind(j-1,i,k) + uwind(j-1,i+1,k) )/2.0)*(fcorio1*Re*cos(lat1) + (uwind(j-1,i,k) + uwind(j-1,i+1,k))*sin(lat1)/2.0) &
                                      + pi_a_est(j,i)*((uwind(j,i,k) + uwind(j,i+1,k) )/2.0)*(fcorio2*Re*cos(lat2) + (uwind(j,i,k) + uwind(j,i+1,k))*sin(lat2)/2.0) &
                                    ) &
                                 -  Re*dlon*cos(lat_edge(j)) &
                                 *  (   (phi(j,i,k) - phi(j-1,i,k))*(pi_a_est(j-1,i) + pi_a_est(j,i))/2.0 &
                                      + (pi_a_est(j,i) - pi_a_est(j-1,i))*(cpd/2.0) &
                                      * ( theta(j-1,i,k)*( sigma(k+1)*(PP_edge(j-1,i,k+1) - PP_mdpt(j-1,i,k)) + sigma(k)*(PP_mdpt(j-1,i,k) - PP_edge(j-1,i,k)) ) &
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
          
          vwind_calc = vwind  
          
      END    
