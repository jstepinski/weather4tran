      FUNCTION pi_a_calc(pi_a_prev)
      
          USE pvars

          IMPLICIT none
          
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_prev
          
          REAL, DIMENSION(Nlat,Nlon) :: pi_a_calc
          
!*****************************************************************!
! Local Variables
          REAL :: S
          INTEGER :: i,j,k
!*****************************************************************!
! Compute column pressure with equation 7.14
          
          DO i = 1,Nlon
            DO j = 1,Nlat
                
                S = 0.
                DO k = 1,Nvert
                    S = S + (abs(sigma(k+1)-sigma(k)))*( fluxF(j,i+1,k) - fluxF(j,i,k) + fluxG(j+1,i,k) - fluxG(j,i,k))
                ENDDO
                
                pi_a_calc(j,i) = pi_a_prev(j,i) - h*S/((Re**2)*cos(lat_mdpt(j))*dlat*dlon)
                
            ENDDO
          ENDDO
          
      END    
