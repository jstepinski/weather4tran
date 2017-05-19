MODULE xfun

! Interface for the extendMat function, which is used in the uwind and vwind 
! calculators to handle boundary conditions

  INTERFACE
            FUNCTION extendMat(M,r,c)
              REAL, INTENT(IN) :: M(:,:)
              INTEGER, INTENT(IN) :: r,c
              REAL extendMat(r+2,c+2)
            END FUNCTION
  END INTERFACE

END MODULE
