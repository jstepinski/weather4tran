      SUBROUTINE PP_calc(pi_a)
  
          USE pvars

          IMPLICIT none

          REAL, DIMENSION(Nlat,Nlon), INTENT(IN) :: pi_a

      !******************************************************************************
      ! Local Variables
          REAL, DIMENSION(Nlat,Nlon,Nvert+1) :: p_edge
          REAL, DIMENSION(Nlat,Nlon) :: Pprev
          REAL, DIMENSION(Nlat,Nlon) :: Pcurr
          INTEGER :: k
      
      !******************************************************************************
      
          PP_mdpt = 0.
          PP_edge = 0.
          
          DO k = 1,Nvert+1
            p_edge(:,:,k) = p_a_top + sigma(k)*pi_a
          ENDDO
      
          PP_edge(:,:,1) = ((p_a_top/100.0)/1000.0)**kappa

          DO k = 1,Nvert
            PP_edge(:,:,k+1) = ((p_edge(:,:,k+1)/100.0)/1000.0)**kappa
            Pcurr = PP_edge(:,:,k+1);
            Pprev = PP_edge(:,:,k);
            PP_mdpt(:,:,k) = (1.0/(1.0+kappa))*( (Pcurr*p_edge(:,:,k+1) - Pprev*p_edge(:,:,k)) &
                                            /(p_edge(:,:,k+1) - p_edge(:,:,k)) )
          ENDDO

          RETURN
      END
