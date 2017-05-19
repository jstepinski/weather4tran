      FUNCTION sigmaDot_calc(pi_a_prev, pi_a_curr, Nvert_edges)
      
          USE pvars

          IMPLICIT none
          
          INTEGER, INTENT(IN) :: Nvert_edges

          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_prev
          REAL, DIMENSION(:,:), INTENT(IN) :: pi_a_curr

          REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: sigmaDot_calc
          
!*****************************************************************!
! Local Variables

          REAL :: S
          INTEGER :: i,j,k,L

!*****************************************************************!
! Calculate vertical velocity with Eq 7.21
          
          sigmaDot_calc(:,:,1) = 0.
          sigmaDot_calc(:,:,Nvert+1) = 0.
      
          DO k = 2,Nvert
            DO i = 1,Nlon
              DO j = 1,Nlat
                    
                S = 0.
                DO L = 1,k-1
                    S = S + (abs(sigma(L+1)-sigma(L)))*( fluxF(j,i+1,L) - fluxF(j,i,L) + fluxG(j+1,i,L) - fluxG(j,i,L))
                ENDDO
                    
                sigmaDot_calc(j,i,k) = -S/(pi_a_curr(j,i)*(Re**2)*cos(lat_mdpt(j))*dlat*dlon) &
                                     - sigma(k)*( (pi_a_curr(j,i) - pi_a_prev(j,i))/(h*pi_a_curr(j,i)) )
                    
              ENDDO
            ENDDO
          ENDDO
      
      
      END
