      FUNCTION extendMat(M,r,c)
      
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: r,c
        REAL, DIMENSION(:,:), INTENT(IN) :: M
        
        REAL, DIMENSION(r+2,c+2) :: extendMat
        
        extendMat(2:r+1,2:c+1) = M
        
        extendMat(:,1)           = extendMat(:,2)
        extendMat(:,c+2)         = extendMat(:,c+1)
        extendMat(1,:)           = extendMat(2,:)
        extendMat(r+2,:)         = extendMat(r+1,:)
        extendMat(1,1)           = extendMat(2,2)
        extendMat(1,c+2)         = extendMat(2,c+1)
        extendMat(r+2,1)         = extendMat(r+1,2)
        extendMat(r+2,c+2)       = extendMat(r+1,c+1)
      
      END
