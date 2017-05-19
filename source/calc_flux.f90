      SUBROUTINE fluxF_calc(pi_a, uwind)
  
          USE pvars

          IMPLICIT none
          
          REAL, DIMENSION(Nlat,Nlon), INTENT(IN) :: pi_a
          REAL, DIMENSION(Nlat,Nlon+1,Nvert), INTENT(IN) :: uwind

      !******************************************************************************
      ! Local Variables

          INTEGER :: i,j, k
      
      !******************************************************************************
      ! Calculate E-W flux with Eq. 7.15 and 7.17
      
          DO k = 1,Nvert
            DO j = 1,Nlat
              DO i = 1,Nlon+1
                    
                IF (i == 1) THEN
                    fluxF(j,i,k) = pi_a(j,i)*Re*dlat*uwind(j,i,k)
                ELSEIF (i == Nlon+1) THEN
                    fluxF(j,i,k) = pi_a(j,i-1)*Re*dlat*uwind(j,i,k)
                ELSE
                    fluxF(j,i,k) = 0.5*(pi_a(j,i-1) + pi_a(j,i))*Re*dlat*uwind(j,i,k)
                ENDIF
                    
              ENDDO
            ENDDO
          ENDDO

          RETURN
      END
      
!******************************************************************************
!******************************************************************************
!******************************************************************************

        SUBROUTINE fluxG_calc(pi_a, vwind)
  
          USE pvars

          IMPLICIT none
          
          REAL, DIMENSION(Nlat,Nlon), INTENT(IN) :: pi_a
          REAL, DIMENSION(Nlat+1,Nlon,Nvert), INTENT(IN) :: vwind

        !*********************************************************************
        ! Local Variables

          INTEGER :: i,j, k
      
        !*********************************************************************
        ! Calculate S-N flux with Eq. 7.16 and 7.18
        
          DO k = 1,Nvert
            DO i = 1,Nlon
              DO j = 1,Nlat+1
                    
                IF (j == 1) THEN
                    fluxG(j,i,k) = pi_a(j,i)*Re*cos(lat_edge(j))*dlon*vwind(j,i,k)
                ELSEIF (j == Nlat+1) THEN
                    fluxG(j,i,k) = pi_a(j-1,i)*Re*cos(lat_edge(j))*dlon*vwind(j,i,k)
                ELSE
                    fluxG(j,i,k) = 0.5*(pi_a(j-1,i) + pi_a(j,i))*Re*cos(lat_edge(j))*dlon*vwind(j,i,k)
                ENDIF
                    
              ENDDO
            ENDDO
          ENDDO

          RETURN
      END
      
