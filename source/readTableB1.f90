SUBROUTINE readTableB1

! Read Table B.1 from "WeatherAtmoData.txt"
 
  USE pvars
  IMPLICIT NONE
  
  INTEGER :: k
  
  OPEN(UNIT=21,FILE='WeatherAtmoData.txt',STATUS='old')
  
  DO k = 1,Ntable
    READ(21,*) zdat(k), gdat(k), pdat(k), Tdat(k), rhodat(k)
  ENDDO

  CLOSE(21)
  
  RETURN

END
