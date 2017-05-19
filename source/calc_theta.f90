      FUNCTION theta_interp_calc(theta)
      
          USE pvars

          IMPLICIT none
          
          !REAL, DIMENSION(Nlat,Nlon,Nvert), INTENT(IN) :: theta
          REAL, DIMENSION(:,:,:), INTENT(IN) :: theta
          
          REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: theta_interp_calc
          
          REAL, DIMENSION(Nlat,Nlon) :: Pk
          REAL, DIMENSION(Nlat,Nlon) :: Pkp1
          REAL, DIMENSION(Nlat,Nlon) :: Pkp12
          INTEGER :: k
          
          theta_interp_calc(:,:,1) = 0.
          DO k = 2,Nvert
            Pkp12 = PP_edge(:,:,k)
            Pk    = PP_mdpt(:,:,k-1)
            Pkp1  = PP_mdpt(:,:,k)
            
            theta_interp_calc(:,:,k) = ( (Pkp12-Pk)*theta(:,:,k-1) + (Pkp1-Pkp12)*theta(:,:,k) )/(Pkp1 - Pk)
          ENDDO
          theta_interp_calc(:,:,Nvert+1) = 0.
          
      END FUNCTION    
      
!**************************************************************************
!**************************************************************************
!**************************************************************************      

      FUNCTION theta_calc(pi_a_prev, pi_a_curr, theta_prev, theta_curr, theta_interp, sigmaDot)
      
          USE pvars

          IMPLICIT none
          
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_prev
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_curr
          REAL, DIMENSION(:,:,:), INTENT(IN) :: theta_prev
          REAL, DIMENSION(:,:,:), INTENT(IN) :: theta_curr
          REAL, DIMENSION(:,:,:), INTENT(IN) :: theta_interp
          REAL, DIMENSION(:,:,:), INTENT(IN) :: sigmaDot
          
          REAL, DIMENSION(Nlat,Nlon,Nvert) :: theta_calc

!*****************************************************************!
! Local Variables
          
          REAL :: theta_j_im1, theta_j_ip1, theta_jm1_i, theta_jp1_i
          INTEGER :: i,j,k

!*****************************************************************!
! Calculate temperature with eq. 7.27
          
          DO k = 1,Nvert
            DO i = 1,Nlon
              DO j = 1,Nlat
            
                ! Temperatures multipled by fluxes F,G in the equation may fall
                ! outside the grid. In that case, set a value at the virtual node
                ! equal to the value at the nearest interior node
            
                if (i == 1) then
                    theta_j_im1 = theta_curr(j,i,k)
                else
                    theta_j_im1 = theta_curr(j,i-1,k)
                endif
                
                if (i == Nlon) then
                    theta_j_ip1 = theta_curr(j,i,k)
                else
                    theta_j_ip1 = theta_curr(j,i+1,k)
                endif
                
                if (j == 1) then
                    theta_jm1_i = theta_curr(j,i,k)
                else
                    theta_jm1_i = theta_curr(j-1,i,k)
                endif
                
                if (j == Nlat) then
                    theta_jp1_i = theta_curr(j,i,k)
                else
                    theta_jp1_i = theta_curr(j+1,i,k)
                endif
            
                theta_calc(j,i,k) = pi_a_prev(j,i)*theta_prev(j,i,k)/pi_a_curr(j,i) &
                                  + (h/(pi_a_curr(j,i)*(Re**2)*cos(lat_mdpt(j))*dlat*dlon)) &
                                  * ( (   fluxF(j,i,k)*(theta_j_im1+theta_curr(j,i,k))/2 - fluxF(j,i+1,k)*(theta_j_ip1+theta_curr(j,i,k))/2 &
                                        + fluxG(j,i,k)*(theta_jm1_i+theta_curr(j,i,k))/2 - fluxG(j+1,i,k)*(theta_jp1_i+theta_curr(j,i,k))/2 &
                                      ) &
                                      + pi_a_curr(j,i)*(Re**2)*cos(lat_mdpt(j))*dlat*dlon &
                                      *( sigmaDot(j,i,k)*theta_interp(j,i,k) - sigmaDot(j,i,k+1)*theta_interp(j,i,k+1) ) &
                                      /( abs(sigma(k+1) - sigma(k)) ) &
                                    )
            
              ENDDO
            ENDDO
          ENDDO
          
      END FUNCTION
