!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 interpol 2006/05/18 13:07:25
!-----------------------------------------------------------------
!###################
MODULE MODI_BIKHARDT
!###################
!
INTERFACE BIKHARDT
!
      SUBROUTINE BIKHARDT4D (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                             PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                             KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,   &
                             HLBCX,HLBCY,PFIELD1,PFIELD2)
!
                                    ! interpolation coefficients 
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model domain, relative to the outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
INTEGER,   INTENT(IN)  :: KGRID      ! code of grid point
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:,:,:),         INTENT(IN) :: PFIELD1 ! field of outer model
REAL, DIMENSION(:,:,:,:),         INTENT(OUT):: PFIELD2 ! field of inner model
!
END SUBROUTINE BIKHARDT4D
!
      SUBROUTINE BIKHARDT3D (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                             PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                             KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,   &
                             HLBCX,HLBCY,PFIELD1,PFIELD2)
!
                                    ! interpolation coefficients  
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
INTEGER,   INTENT(IN)  :: KGRID      ! code of grid point
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:,:),           INTENT(IN) :: PFIELD1 ! field of outer model
REAL, DIMENSION(:,:,:),           INTENT(OUT):: PFIELD2 ! field of inner model
!
END SUBROUTINE BIKHARDT3D
!
      SUBROUTINE BIKHARDT2D (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                             PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                             KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,   &
                             HLBCX,HLBCY,PFIELD1,PFIELD2)
!
                                    ! interpolation coefficients  
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model  and outer model
INTEGER,   INTENT(IN)  :: KGRID      ! code of grid point
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:),             INTENT(IN) :: PFIELD1 ! field of outer model
REAL, DIMENSION(:,:),             INTENT(OUT):: PFIELD2 ! field of inner model
!
END SUBROUTINE BIKHARDT2D
!
END INTERFACE
!
END MODULE MODI_BIKHARDT
!
!#####################
MODULE MODI_BIKHARDT4D
!#####################
!
INTERFACE
!
      SUBROUTINE BIKHARDT4D (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                             PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                             KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,   &
                             HLBCX,HLBCY,PFIELD1,PFIELD2)
!
                                    ! interpolation coefficients  
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model  and outer model
INTEGER,   INTENT(IN)  :: KGRID      ! code of grid point
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:,:,:),         INTENT(IN) :: PFIELD1 ! field of outer model
REAL, DIMENSION(:,:,:,:),         INTENT(OUT):: PFIELD2 ! field of inner model
!
END SUBROUTINE BIKHARDT4D
!
END INTERFACE
!
END MODULE MODI_BIKHARDT4D
!
!
!     #########################################################################
      SUBROUTINE BIKHARDT4D (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                             PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                             KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,   &
                             HLBCX,HLBCY,PFIELD1,PFIELD2)
!     #########################################################################
!
!!****  *BIKHARDT4D * - interpolates a 4D field with Bikhardt method
!!
!!    PURPOSE
!!    -------
!!
!     This routine interpolates a field from outer model (PFIELD1) to 
!     inner model (PFIELD2), using Bikhardt interpolation.
!!
!!**  METHOD
!!    ------
!!
!!    The outer model field is extrapolated in each horizontal dimension
!!    in order to allow all cases of KXOR,KYOR,KXEND and KYEND
!!
!!    EXTERNAL
!!    --------
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      Module MODD_PARAMETERS : contains parameters 
!!
!!    REFERENCE
!!    ---------
!!
!!       Book1 of the documentation
!!       Routine BIKHARDT3D (Book2 of the documentation)
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     10/06/96
!!      J.P. Lafore  22/10/96  interpolation coefficients added to the arguments
!!                             list to avoid duplication. 
!!     V. Masson and F. Gheusi (10/10/97) bug in cyclic case
!!     J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!     J.Escobar : 18/12/2015 : set valide default values in corner in // for NHALO <>1 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_PARAMETERS       ! Declarative modules
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
                                    ! interpolation coefficients  
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model  domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
INTEGER,   INTENT(IN)  :: KGRID      ! code of grid point
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:,:,:),         INTENT(IN) :: PFIELD1 ! field of outer model
REAL, DIMENSION(:,:,:,:),         INTENT(OUT):: PFIELD2 ! field of inner model
!
!*       0.2    Declarations of local variables for print on FM file
!
REAL, DIMENSION(0:SIZE(PFIELD1,1)+2, &
                0:SIZE(PFIELD1,2)+2, &
                  SIZE(PFIELD1,3)  , &
                  SIZE(PFIELD1,4)   ) :: ZFIELD1   ! field of outer model 
!
                                    ! interpolation coefficients 
REAL :: ZBMX1,ZBMX2,ZBMX3,ZBMX4     ! at Mass points in X-direc.
REAL :: ZBMY1,ZBMY2,ZBMY3,ZBMY4     ! at Mass points in Y-direc.
REAL :: ZBFX1,ZBFX2,ZBFX3,ZBFX4     ! at Flux points in X-direc.
REAL :: ZBFY1,ZBFY2,ZBFY3,ZBFY4     ! at Flux points in Y-direc.
! 
INTEGER             :: IIU       ! Upper dimension in x direction (inner model)
INTEGER             :: IJU       ! Upper dimension in y direction (inner model)
INTEGER             :: IIB, IIE
INTEGER             :: IJB, IJE
INTEGER             :: IIU1      ! Upper dimension in x direction (outer model)
INTEGER             :: IJU1      ! Upper dimension in y direction (outer model)
INTEGER             :: IIS,IJS   ! indices I and J in x and y dir. for scalars
INTEGER             :: IIF,IJF   ! indices I and J in x and y dir. for flux points
INTEGER             :: JI, JEPSX ! Loop index in x direction
INTEGER             :: JJ, JEPSY ! Loop index in y direction    
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE:
!              ---------
! 
!*       1.1   computes dimensions of arrays and other indices
!
IIU = SIZE(PFIELD2,1)
IJU = SIZE(PFIELD2,2)
IIU1= SIZE(PFIELD1,1)
IJU1= SIZE(PFIELD1,2)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
!*       1.2   extrapolates field of outer model
!
ZFIELD1(:,:,:,:) = 0.
ZFIELD1(1:IIU1,1:IJU1,:,:)=PFIELD1(:,:,:,:)
!
ZFIELD1(   0  ,   :  ,:,:) = 2.*ZFIELD1(  1 ,  : ,:,:) -    ZFIELD1(   2  ,   :  ,:,:)
ZFIELD1(IIU1+1,   :  ,:,:) = 2.*ZFIELD1(IIU1,  : ,:,:) -    ZFIELD1(IIU1-1,   :  ,:,:)
ZFIELD1(IIU1+2,   :  ,:,:) = 3.*ZFIELD1(IIU1,  : ,:,:) - 2.*ZFIELD1(IIU1-1,   :  ,:,:)
!
ZFIELD1(   :  ,   0  ,:,:) = 2.*ZFIELD1(  : ,  1 ,:,:) -    ZFIELD1(   :  ,   2  ,:,:)
ZFIELD1(   :  ,IJU1+1,:,:) = 2.*ZFIELD1(  : ,IJU1,:,:) -    ZFIELD1(   :  ,IJU1-1,:,:)
ZFIELD1(   :  ,IJU1+2,:,:) = 3.*ZFIELD1(  : ,IJU1,:,:) - 2.*ZFIELD1(   :  ,IJU1-1,:,:)
!
IF ( HLBCX(1) == 'CYCL' ) THEN
  ZFIELD1(   0  ,:,:,:) = ZFIELD1(IIU1-2*JPHEXT  ,:,:,:)
  ZFIELD1(IIU1+1,:,:,:) = ZFIELD1(  1 +2*JPHEXT  ,:,:,:)
  IF (SIZE(PFIELD1,1) == 3 ) THEN
    ZFIELD1(IIU1+2,:,:,:) = ZFIELD1(  1            ,:,:,:)
  ELSE
    ZFIELD1(IIU1+2,:,:,:) = ZFIELD1(  1 +2*JPHEXT+1,:,:,:)
  END IF
END IF
!
IF ( HLBCY(1) == 'CYCL' ) THEN
  ZFIELD1(:,   0  ,:,:) = ZFIELD1(:,IJU1-2*JPHEXT  ,:,:)
  ZFIELD1(:,IJU1+1,:,:) = ZFIELD1(:,  1 +2*JPHEXT  ,:,:)
  IF (SIZE(PFIELD1,2) == 3 ) THEN
    ZFIELD1(:,IJU1+2,:,:) = ZFIELD1(:,  1            ,:,:)
  ELSE
    ZFIELD1(:,IJU1+2,:,:) = ZFIELD1(:,  1 +2*JPHEXT+1,:,:)
  END IF
END IF
!-------------------------------------------------------------------------------
!
PFIELD2 = ZFIELD1(1,1,1,1) ! some valide values for missing ones
!
SELECT CASE (KGRID)
!
!*      2.1    Mass points
!
       CASE (1,4)
!
  DO JEPSX = 1,KDXRATIO
    ZBMX1 = PBMX1(JEPSX)
    ZBMX2 = PBMX2(JEPSX)
    ZBMX3 = PBMX3(JEPSX)
    ZBMX4 = PBMX4(JEPSX)
    DO JEPSY = 1,KDYRATIO
      ZBMY1 = PBMY1(JEPSY)
      ZBMY2 = PBMY2(JEPSY)
      ZBMY3 = PBMY3(JEPSY)
      ZBMY4 = PBMY4(JEPSY)
      DO JI = KXOR,KXEND
        IIS = IIB+JEPSX-1+KDXRATIO/2+(JI-KXOR-JPHEXT)*KDXRATIO
        DO JJ = KYOR,KYEND
          IJS = IJB+JEPSY-1+KDYRATIO/2+(JJ-KYOR-JPHEXT)*KDYRATIO
!
          IF (1 <= IIS .AND. IIS <= IIU .AND. 1 <= IJS .AND. IJS <= IJU) THEN
!
            PFIELD2  (IIS,IJS,:,:) = ZBMY1*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ-1,:,:)+ZBMX2*ZFIELD1(JI  ,JJ-1,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ-1,:,:)+ZBMX4*ZFIELD1(JI+2,JJ-1,:,:)) &
                                    +ZBMY2*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ  ,:,:)+ZBMX2*ZFIELD1(JI  ,JJ  ,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ  ,:,:)+ZBMX4*ZFIELD1(JI+2,JJ  ,:,:)) &
                                    +ZBMY3*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ+1,:,:)+ZBMX2*ZFIELD1(JI  ,JJ+1,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ+1,:,:)+ZBMX4*ZFIELD1(JI+2,JJ+1,:,:)) &
                                    +ZBMY4*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ+2,:,:)+ZBMX2*ZFIELD1(JI  ,JJ+2,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ+2,:,:)+ZBMX4*ZFIELD1(JI+2,JJ+2,:,:))
           END IF
         END DO
      END DO
    END DO
  END DO
!
!*      2.2    U points
!
       CASE (2,6)
!
  DO JEPSX = 1,KDXRATIO
    ZBFX1 = PBFX1(JEPSX)
    ZBFX2 = PBFX2(JEPSX)
    ZBFX3 = PBFX3(JEPSX)
    ZBFX4 = PBFX4(JEPSX)
    DO JEPSY = 1,KDYRATIO
      ZBMY1 = PBMY1(JEPSY)
      ZBMY2 = PBMY2(JEPSY)
      ZBMY3 = PBMY3(JEPSY)
      ZBMY4 = PBMY4(JEPSY)
      DO JI = KXOR,KXEND
        IIF = IIB+JEPSX-1           +(JI-KXOR-JPHEXT)*KDXRATIO
        DO JJ = KYOR,KYEND
          IJS = IJB+JEPSY-1+KDYRATIO/2+(JJ-KYOR-JPHEXT)*KDYRATIO

          IF (1 <= IIF .AND. IIF <= IIU .AND. 1 <= IJS .AND. IJS <= IJU) THEN 
!
            PFIELD2  (IIF,IJS,:,:) = ZBMY1*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ-1,:,:)+ZBFX2*ZFIELD1(JI  ,JJ-1,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ-1,:,:)+ZBFX4*ZFIELD1(JI+2,JJ-1,:,:)) &
                                    +ZBMY2*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ  ,:,:)+ZBFX2*ZFIELD1(JI  ,JJ  ,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ  ,:,:)+ZBFX4*ZFIELD1(JI+2,JJ  ,:,:)) &
                                    +ZBMY3*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ+1,:,:)+ZBFX2*ZFIELD1(JI  ,JJ+1,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ+1,:,:)+ZBFX4*ZFIELD1(JI+2,JJ+1,:,:)) &
                                    +ZBMY4*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ+2,:,:)+ZBFX2*ZFIELD1(JI  ,JJ+2,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ+2,:,:)+ZBFX4*ZFIELD1(JI+2,JJ+2,:,:))
          END IF
        END DO
      END DO
    END DO
  END DO
!
!*      2.3    V points
!
       CASE (3,7)
!
  DO JEPSX = 1,KDXRATIO
    ZBMX1 = PBMX1(JEPSX)
    ZBMX2 = PBMX2(JEPSX)
    ZBMX3 = PBMX3(JEPSX)
    ZBMX4 = PBMX4(JEPSX)
    DO JEPSY = 1,KDYRATIO
      ZBFY1 = PBFY1(JEPSY)
      ZBFY2 = PBFY2(JEPSY)
      ZBFY3 = PBFY3(JEPSY)
      ZBFY4 = PBFY4(JEPSY)
      DO JI = KXOR,KXEND
        IIS = IIB+JEPSX-1+KDXRATIO/2+(JI-KXOR-JPHEXT)*KDXRATIO
        DO JJ = KYOR,KYEND
          IJF = IJB+JEPSY-1           +(JJ-KYOR-JPHEXT)*KDYRATIO

          IF (1 <= IIS .AND. IIS <= IIU .AND. 1 <= IJF .AND. IJF <= IJU) THEN 
!
            PFIELD2  (IIS,IJF,:,:) = ZBFY1*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ-1,:,:)+ZBMX2*ZFIELD1(JI  ,JJ-1,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ-1,:,:)+ZBMX4*ZFIELD1(JI+2,JJ-1,:,:)) &
                                    +ZBFY2*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ  ,:,:)+ZBMX2*ZFIELD1(JI  ,JJ  ,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ  ,:,:)+ZBMX4*ZFIELD1(JI+2,JJ  ,:,:)) &
                                    +ZBFY3*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ+1,:,:)+ZBMX2*ZFIELD1(JI  ,JJ+1,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ+1,:,:)+ZBMX4*ZFIELD1(JI+2,JJ+1,:,:)) &
                                    +ZBFY4*                               &
             ( ZBMX1*ZFIELD1(JI-1,JJ+2,:,:)+ZBMX2*ZFIELD1(JI  ,JJ+2,:,:)  &
              +ZBMX3*ZFIELD1(JI+1,JJ+2,:,:)+ZBMX4*ZFIELD1(JI+2,JJ+2,:,:))
          END IF
        END DO
      END DO
    END DO
  END DO
!
!
!*      2.4    vertical vorticity points
!
       CASE (5,8)
!
  DO JEPSX = 1,KDXRATIO
    ZBFX1 = PBFX1(JEPSX)
    ZBFX2 = PBFX2(JEPSX)
    ZBFX3 = PBFX3(JEPSX)
    ZBFX4 = PBFX4(JEPSX)
    DO JEPSY = 1,KDYRATIO
      ZBFY1 = PBFY1(JEPSY)
      ZBFY2 = PBFY2(JEPSY)
      ZBFY3 = PBFY3(JEPSY)
      ZBFY4 = PBFY4(JEPSY)
      DO JI = KXOR,KXEND
        IIF = IIB+JEPSX-1           +(JI-KXOR-JPHEXT)*KDXRATIO
        DO JJ = KYOR,KYEND
          IJF = IJB+JEPSY-1           +(JJ-KYOR-JPHEXT)*KDYRATIO

          IF (1 <= IIF .AND. IIF <= IIU .AND. 1 <= IJF .AND. IJF <= IJU) THEN 
!
            PFIELD2  (IIF,IJF,:,:) = ZBFY1*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ-1,:,:)+ZBFX2*ZFIELD1(JI  ,JJ-1,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ-1,:,:)+ZBFX4*ZFIELD1(JI+2,JJ-1,:,:)) &
                                    +ZBFY2*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ  ,:,:)+ZBFX2*ZFIELD1(JI  ,JJ  ,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ  ,:,:)+ZBFX4*ZFIELD1(JI+2,JJ  ,:,:)) &
                                    +ZBFY3*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ+1,:,:)+ZBFX2*ZFIELD1(JI  ,JJ+1,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ+1,:,:)+ZBFX4*ZFIELD1(JI+2,JJ+1,:,:)) &
                                    +ZBFY4*                               &
             ( ZBFX1*ZFIELD1(JI-1,JJ+2,:,:)+ZBFX2*ZFIELD1(JI  ,JJ+2,:,:)  &
              +ZBFX3*ZFIELD1(JI+1,JJ+2,:,:)+ZBFX4*ZFIELD1(JI+2,JJ+2,:,:))
          END IF
        END DO
      END DO
    END DO
  END DO
!
END SELECT
!-------------------------------------------------------------------------------
!
!
!
END SUBROUTINE BIKHARDT4D
!  
!     #########################################################################
      SUBROUTINE BIKHARDT3D (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                             PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                             KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,   &
                             HLBCX,HLBCY,PFIELD1,PFIELD2)
!     #########################################################################
!
!!****  *BIKHARDT * - interpolates with Bikhardt method
!!
!!    PURPOSE
!!    -------
!!
!     This routine interpolates a field from outer model (PFIELD1) to 
!     inner model (PFIELD2), using Bikhardt interpolation.
!!
!!**  METHOD
!!    ------
!!
!!    The outer model field is extrapolated in each horizontal dimension
!!    in order to allow all cases of KXOR,KYOR,KXEND and KYEND
!!
!!    EXTERNAL
!!    --------
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!! 
!!
!!    REFERENCE
!!    ---------
!!
!!       Book1 of the documentation
!!       Routine BIKHARDT (Book2 of the documentation)
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     10/06/96 
!!      J.P. Lafore  22/10/96  interpolation coefficients added to the arguments
!!                             list to avoid duplication. 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_BIKHARDT4D      
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
                                    ! interpolation coefficients  
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
INTEGER,   INTENT(IN)  :: KGRID      ! code of grid point
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:,:),           INTENT(IN) :: PFIELD1 ! field of outer model
REAL, DIMENSION(:,:,:),           INTENT(OUT):: PFIELD2 ! field of inner model
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PFIELD1,1),SIZE(PFIELD1,2),SIZE(PFIELD1,3),1) :: ZFIELD1 
REAL, DIMENSION(SIZE(PFIELD2,1),SIZE(PFIELD2,2),SIZE(PFIELD2,3),1) :: ZFIELD2 
!
!-------------------------------------------------------------------------------
ZFIELD1(:,:,:,1)=PFIELD1(:,:,:)
CALL BIKHARDT4D(PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
      KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,HLBCX,HLBCY,ZFIELD1,ZFIELD2)
PFIELD2(:,:,:)  =ZFIELD2(:,:,:,1)
!-------------------------------------------------------------------------------
!
END SUBROUTINE BIKHARDT3D
!
!     #########################################################################
      SUBROUTINE BIKHARDT2D (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                             PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                             KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,   &
                             HLBCX,HLBCY,PFIELD1,PFIELD2)
!     #########################################################################
!
!!****  *BIKHARDT * - interpolates with Bikhardt method
!!
!!    PURPOSE
!!    -------
!!
!     This routine interpolates a field from outer model (PFIELD1) to 
!     inner model (PFIELD2), using Bikhardt interpolation.
!!
!!**  METHOD
!!    ------
!!
!!    The outer model field is extrapolated in each horizontal dimension
!!    in order to allow all cases of KXOR,KYOR,KXEND and KYEND
!!
!!    EXTERNAL
!!    --------
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!    REFERENCE
!!    ---------
!!
!!       Book1 of the documentation
!!       Routine BIKHARDT (Book2 of the documentation)
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     10/06/96 
!!      J.P. Lafore  22/10/96  interpolation coefficients added to the arguments
!!                             list to avoid duplication. 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_BIKHARDT4D      
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
                                    ! interpolation coefficients  
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
INTEGER,   INTENT(IN)  :: KGRID      ! code of grid point
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:),             INTENT(IN) :: PFIELD1 ! field of outer model
REAL, DIMENSION(:,:),             INTENT(OUT):: PFIELD2 ! field of inner model
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PFIELD1,1),SIZE(PFIELD1,2),1,1) :: ZFIELD1 
REAL, DIMENSION(SIZE(PFIELD2,1),SIZE(PFIELD2,2),1,1) :: ZFIELD2 
!
!-------------------------------------------------------------------------------
ZFIELD1(:,:,1,1)=PFIELD1(:,:)
CALL BIKHARDT4D(PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
      KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,KGRID,HLBCX,HLBCY,ZFIELD1,ZFIELD2)
PFIELD2(:,:)  =ZFIELD2(:,:,1,1)
!-------------------------------------------------------------------------------
!
END SUBROUTINE BIKHARDT2D
