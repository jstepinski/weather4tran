SUBROUTINE printToFile(uwind,vwind,theta,pi_a,sigmaDot,t)
        
        ! This particular version of the file printer was used for my 
        ! error analysis. It would take the 3D vectors at specific times
        ! and print specific levels as 2D arrays. The names were constructed
        ! to help me distinguish tests. "L" is the sigma level. "t" is the 
        ! number being printed. "h" helps me keep track of what time step 
        ! I set.
        
        USE uvars

        IMPLICIT NONE

        REAL, DIMENSION(Nlat,Nlon), INTENT(IN) :: pi_a
        REAL, DIMENSION(Nlat,Nlon+1,Nvert), INTENT(IN) :: uwind
        REAL, DIMENSION(Nlat+1,Nlon,Nvert), INTENT(IN) :: vwind
        REAL, DIMENSION(Nlat,Nlon,Nvert), INTENT(IN) :: theta
        REAL, DIMENSION(Nlat,Nlon,Nvert+1), INTENT(IN) :: sigmaDot

        INTEGER, INTENT(IN) :: t

        CHARACTER(LEN=100) :: fname

        INTEGER i,j,k

        DO k = 1,npL
                WRITE(fname, '( "pi_a_L", I2.2, "_t", I1.1, "_h", I1.1, ".txt" )' ) printLevels(k), t, hindex
                OPEN( UNIT = 220, FILE = fname, STATUS = 'UNKNOWN' )
                DO i = 1,Nlat
                        DO j = 1,Nlon
                        WRITE(220, '(f20.10,2X)', advance='no') pi_a(i,j)
                        ENDDO
                        WRITE(220, *) ''
                ENDDO
                CLOSE(220)


                WRITE(fname, '( "uwind_L", I2.2, "_t", I1.1, "_h", I1.1, ".txt" )' ) printLevels(k), t, hindex
                OPEN( UNIT = 220, FILE = fname, STATUS = 'UNKNOWN' )
                DO i = 1,Nlat
                        DO j = 1,Nlon+1
                        WRITE(220, '(f16.10,2X)', advance='no') uwind(i,j,printLevels(k))
                        ENDDO
                        WRITE(220, *) ''
                ENDDO
                CLOSE(220)


                WRITE(fname, '( "vwind_L", I2.2, "_t", I1.1, "_h", I1.1, ".txt" )' ) printLevels(k), t, hindex
                OPEN( UNIT = 220, FILE = fname, STATUS = 'UNKNOWN' )
                DO i = 1,Nlat+1
                        DO j = 1,Nlon
                        WRITE(220, '(f20.10,2X)', advance='no') vwind(i,j,printLevels(k))
                        ENDDO
                        WRITE(220, *) ''
                ENDDO
                CLOSE(220)


                WRITE(fname, '( "theta_L", I2.2, "_t", I1.1, "_h", I1.1, ".txt" )' ) printLevels(k), t, hindex
                OPEN( UNIT = 220, FILE = fname, STATUS = 'UNKNOWN' )
                DO i = 1,Nlat
                        DO j = 1,Nlon
                        WRITE(220, '(f20.10,2X)', advance='no') theta(i,j,printLevels(k))
                        ENDDO
                        WRITE(220, *) ''
                ENDDO
                CLOSE(220)

                ! Note that because sigmaDot is evaluated on vertical edges, it does not exactly coincide with
                ! the other variables. One can choose to print it at the lower or upper edge of a cell. In this case,
                ! we use the upper edge, because one of the designated print levels is 1.
                WRITE(fname, '( "sigmaDot_L", I2.2, "_t", I1.1, "_h", I1.1, ".txt" )' ) printLevels(k)+1, t, hindex
                OPEN( UNIT = 220, FILE = fname, STATUS = 'UNKNOWN' )
                DO i = 1,Nlat
                        DO j = 1,Nlon
                        WRITE(220, '(f20.10,2X)', advance='no') sigmaDot(i,j,printLevels(k))
                        ENDDO
                        WRITE(220, *) ''
                ENDDO
                CLOSE(220)

        ENDDO


END
