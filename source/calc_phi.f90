      FUNCTION phi_calc(theta)
      
          USE pvars

          IMPLICIT none
          
          REAL, DIMENSION(:,:,:), INTENT(IN) :: theta          

          REAL, DIMENSION(Nlat,Nlon,Nvert) :: phi_calc
          
!*****************************************************************!
! Local Variables
          REAL, DIMENSION(Nlat,Nlon,Nvert) :: phi_mdpt = 0.
          REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: phi_edge = 0.
          
          INTEGER :: k

!*****************************************************************!
! Calculate geopotential with Equations 7.61-7.63
          
          phi_edge(:,:,Nvert+1) = phi_surf
          phi_mdpt(:,:,Nvert) = phi_surf - cpd*(theta(:,:,Nvert)*(PP_mdpt(:,:,Nvert) - PP_edge(:,:,Nvert+1) ))
         
          phi_edge(:,:,Nvert) = phi_mdpt(:,:,Nvert) - cpd*theta(:,:,Nvert)*(PP_edge(:,:,Nvert) - PP_mdpt(:,:,Nvert))
          
          DO k = Nvert-1,1,-1
            
            phi_mdpt(:,:,k) = phi_edge(:,:,k+1) - cpd*theta(:,:,k)*(PP_mdpt(:,:,k) - PP_edge(:,:,k+1))
            phi_edge(:,:,k) = phi_mdpt(:,:,k)   - cpd*theta(:,:,k)*(PP_edge(:,:,k) - PP_mdpt(:,:,k))
            
          ENDDO
          
          phi_calc = phi_mdpt
          
      END    
