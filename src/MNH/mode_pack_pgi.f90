!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
#ifdef MNH_PGI
!###################
MODULE MODE_PACK_PGI
!###################
!
!
!!    AUTHOR
!!    ------
!!      J.ESCOBAR       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!   Juan 24/09/2012: for BUG Pgi rewrite PACK function on mode_pack_pgi

IMPLICIT NONE

  INTERFACE PACK
     MODULE PROCEDURE PACK_I,PACK_L,PACK_X1,PACK_X2,PACK_X2S,PACK_X3
  END INTERFACE PACK

CONTAINS 

  FUNCTION PACK_I(KTAB, MASK)
    IMPLICIT NONE
    
    INTEGER , DIMENSION (:), INTENT(IN) :: KTAB
    LOGICAL , DIMENSION (:), INTENT(IN) :: MASK
    INTEGER , DIMENSION (SIZE(MASK))    :: PACK_I
    
    INTEGER :: JI,JL
    
    JL = 0
    DO JI = 1 , SIZE(MASK)
       IF ( MASK(JI) ) THEN
          JL = JL + 1   
          PACK_I(JL) =  KTAB(JI)
       END IF
    END DO
    
  END FUNCTION PACK_I

 FUNCTION PACK_L(GTAB, MASK)
    IMPLICIT NONE
    
    LOGICAL , DIMENSION (:), INTENT(IN) :: GTAB
    LOGICAL , DIMENSION (:), INTENT(IN) :: MASK
    LOGICAL , DIMENSION (SIZE(MASK))   :: PACK_L
    
    INTEGER :: JI,JL
    
    JL = 0
    DO JI = 1 ,  SIZE(MASK)
       IF ( MASK(JI) ) THEN
          JL = JL + 1   
          PACK_L(JL) =  GTAB(JI)
       END IF
    END DO
    
  END FUNCTION PACK_L

  FUNCTION PACK_X1(PTAB, MASK)
    IMPLICIT NONE
    
    REAL    , DIMENSION (:), INTENT(IN)  :: PTAB
    LOGICAL , DIMENSION (:), INTENT(IN)  :: MASK
    REAL    , DIMENSION (:), ALLOCATABLE :: PACK_X1
  
    REAL    , DIMENSION (SIZE(MASK))    :: PACK_X1TEMP
  
    INTEGER :: JI,JL
    
    JL = 0
    DO JI = 1 , SIZE(MASK)
       IF ( MASK(JI) ) THEN
          JL = JL + 1   
          PACK_X1TEMP(JL) =  PTAB(JI)
       END IF
    END DO
    
    ALLOCATE(PACK_X1(JL))
    PACK_X1(1:JL) = PACK_X1TEMP(1:JL)

  END FUNCTION PACK_X1
  
  FUNCTION PACK_X2(PTAB, MASK)
    IMPLICIT NONE
    
    REAL    , DIMENSION (:,:), INTENT(IN) :: PTAB
    LOGICAL , DIMENSION (:,:), INTENT(IN) :: MASK
    REAL    , DIMENSION (:)  ,ALLOCATABLE :: PACK_X2

    REAL    , DIMENSION (SIZE(MASK))      :: PACK_X2TEMP
    
    INTEGER :: JI,JJ,JL
    
    JL = 0
    DO JJ = 1 , SIZE(MASK,2)
       DO JI = 1 , SIZE(MASK,1)
          IF ( MASK(JI,JJ) ) THEN
             JL = JL + 1   
             PACK_X2TEMP(JL) =  PTAB(JI,JJ)
          END IF
       END DO
    END DO

    ALLOCATE(PACK_X2(JL))
    PACK_X2(1:JL) = PACK_X2TEMP (1:JL)
    
  END FUNCTION PACK_X2
  
  FUNCTION PACK_X2S(PTAB, MASK)
    IMPLICIT NONE
    
    REAL    , DIMENSION (:,:), INTENT(IN) :: PTAB
    LOGICAL                               :: MASK
    REAL    , DIMENSION (:)  ,ALLOCATABLE :: PACK_X2S

    REAL    , DIMENSION (SIZE(PTAB))      :: PACK_X2STEMP
    
    INTEGER :: JI,JJ,JL
    
    JL = 0
    DO JJ = 1 , SIZE(PTAB,2)
       DO JI = 1 , SIZE(PTAB,1)
          IF ( MASK ) THEN
             JL = JL + 1   
             PACK_X2STEMP(JL) =  PTAB(JI,JJ)
          END IF
       END DO
    END DO
    
    ALLOCATE(PACK_X2S(JL))
    PACK_X2S(1:JL) = PACK_X2STEMP (1:JL)

  END FUNCTION PACK_X2S

  FUNCTION PACK_X3(PTAB,MASK)
  IMPLICIT NONE
  
  REAL    , DIMENSION (:,:,:) , INTENT(IN)  :: PTAB
  LOGICAL , DIMENSION (:,:,:) , INTENT(IN)  :: MASK
  REAL    , DIMENSION (:)     ,ALLOCATABLE  :: PACK_X3

  REAL    , DIMENSION (SIZE(MASK))          :: PACK_X3TEMP
  
  INTEGER                                   :: II,IJ,IK, JL
  
  JL = 0
  DO IK=1,SIZE(MASK,3)
     DO IJ=1,SIZE(MASK,2)
        DO II=1,SIZE(MASK,1)
           IF (MASK(II,IJ,IK) ) THEN
              JL = JL+1
              PACK_X3TEMP(JL) = PTAB(II,IJ,IK) 
           END IF
        END DO
     END DO
  END DO

  ALLOCATE(PACK_X3(JL))
  PACK_X3(1:JL) = PACK_X3TEMP (1:JL)
  
END FUNCTION PACK_X3

END MODULE MODE_PACK_PGI
#endif
