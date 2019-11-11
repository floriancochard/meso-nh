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
!############################
MODULE MODI_BILINEAR
!############################
!
INTERFACE BILINEAR
!
     SUBROUTINE BILINEAR_4D (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                             PFIELD1,PFIELD2                          )
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PFIELD1 ! FIELD on model 1
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PFIELD2 ! FIELD on model 2
!
END SUBROUTINE BILINEAR_4D
!
     SUBROUTINE BILINEAR_3D (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                             PFIELD1,PFIELD2                          )
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PFIELD1 ! FIELD on model 1
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELD2 ! FIELD on model 2
!
END SUBROUTINE BILINEAR_3D
!
     SUBROUTINE BILINEAR_2D (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                             PFIELD1,PFIELD2                          )
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PFIELD1 ! FIELD on model 1
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD2 ! FIELD on model 2
!
END SUBROUTINE BILINEAR_2D
!
END INTERFACE
!
END MODULE MODI_BILINEAR
!
!######################
MODULE MODI_BILINEAR_4D
!######################
!
INTERFACE
!
     SUBROUTINE BILINEAR_4D (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                             PFIELD1,PFIELD2                          )
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PFIELD1 ! FIELD on model 1
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PFIELD2 ! FIELD on model 2
!
END SUBROUTINE BILINEAR_4D
!
END INTERFACE
!
END MODULE MODI_BILINEAR_4D
!
!     #########################################################################
     SUBROUTINE BILINEAR_4D (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                             PFIELD1,PFIELD2                          )
!     #########################################################################
!
!!****  *BILINEAR * - subroutine to interpolate surface FIELD
!!
!!    PURPOSE
!!    -------
!!
!!      This routine interpolates the FIELD mantel from model 1 grid mesh
!!    to model 2 grid mesh.
!!
!!
!!**  METHOD
!!    ------
!!
!!     Interpolation is bilinear, and uses 9 grid points, located in the
!!   center of model 1 grid mesh, and at the boundaries of this grid
!!   mesh (2 X limits, 2 Y limits and 4 corners).
!!     This implies that the grid mesh values located around the model 1
!!   grid mesh are not used directly. The values at the boundaries of the
!!   grid mesh are defined by the average between the middle point
!!   (this grid mesh value), and the one in the considered direction.
!!   So the eight grid meshes around the considered grid mesh are used
!!   equally.
!!     This is important to note that these average values are erased
!!   and replaced by zero if they are at the limit of any grid
!!   mesh with the zero value. This allows to insure zero value in model 2
!!   grid meshes where there was not the considered class in corresponding
!!   model 1 grid mesh, and to insure continuity of the FIELD type
!!   at such boundaries.
!!
!!
!!    The arrays and array index are defined on the following (model1) grid:
!!
!!
!!        XFIELD                    XFIELD                    XFIELD
!!          *                         *                         *
!!       i-1,j+1                    i,j+1                    i+1,j+1
!!
!!
!!
!!                   ZFIELD_XY     ZFIELD_Y    ZFIELD_XY
!!                       *            *            *
!!                     i,j+1        i,j+1       i+1,j+1
!!
!!
!!
!!        XFIELD      ZFIELD_X      XFIELD      ZFIELD_X      XFIELD
!!          *            *            *            *            *
!!        i-1,j         i,j          i,j         i+1,j        i+1,j
!!
!!
!!
!!                    ZFIELD_XY     ZFIELD_Y    ZFIELD_XY
!!                       *            *            *
!!                      i,j          i,j         i+1,j
!!
!!
!!
!!        XFIELD                    XFIELD                    XFIELD
!!          *                         *                         *
!!       i-1,j-1                    i,j-1                    i+1,j-1
!!
!!
!!
!!
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
!!      Original     12/12/02 
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
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PFIELD1 ! FIELD on model 1
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PFIELD2 ! FIELD on model 2
!
!
!
!*       0.2    Declarations of local variables for print on FM file
!
! 
REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ZW       ! weight. 1 if defined, 
                                                   !         0 if FIELD=XUNDEF
REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ZW_FIELD ! weight. * FIELD
REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ZFIELD_X ! FIELD at mesh interface
REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ZFIELD_Y ! FIELD at mesh interface
REAL, DIMENSION (:,:,:,:), ALLOCATABLE :: ZFIELD_XY! FIELD at mesh corner
!
REAL, DIMENSION(KDXRATIO)            :: ZCOEF_X  ! coefficient array in X direction
REAL, DIMENSION(KDYRATIO)            :: ZCOEF_Y  ! coefficient array in Y direction
!
REAL                                 :: ZC1_X    ! coefficient for left   points
REAL                                 :: ZC2_X    ! coefficient for middle points
REAL                                 :: ZC3_X    ! coefficient for right  points
REAL                                 :: ZC1_Y    ! coefficient for bottom points
REAL                                 :: ZC2_Y    ! coefficient for middle points
REAL                                 :: ZC3_Y    ! coefficient for top    points
!
INTEGER                              :: IIU       ! model 1 X size
INTEGER                              :: IJU       ! model 1 Y size
!
INTEGER                              :: JI,JEPSX  ! Loop index in x direction
INTEGER                              :: JJ,JEPSY  ! Loop index in y direction
!
INTEGER                              :: IIS       ! localization of model 2
INTEGER                              :: IJS       ! grid point
INTEGER                              :: JK, JL    ! loop counters
!
REAL                                 :: ZEPS=1.E-3
!-------------------------------------------------------------------------------
!
IIU=SIZE(PFIELD1,1)
IJU=SIZE(PFIELD1,2)
!
!* weighting factor
!
ALLOCATE(ZW(IIU,IJU,SIZE(PFIELD1,3),SIZE(PFIELD1,4)))
WHERE (PFIELD1/=XUNDEF)
  ZW=1.
ELSEWHERE
  ZW=0.
END WHERE
!
!* weighted FIELD cover
!
ALLOCATE(ZW_FIELD(IIU,IJU,SIZE(PFIELD1,3),SIZE(PFIELD1,4)))
ZW_FIELD=ZW*PFIELD1
!
!-------------------------------------------------------------------------------
!
!*       1.     FIELD type at grid mesh interfaces (in X directions)
!               ----------------------------------
!
ALLOCATE(ZFIELD_X(IIU+1,IJU,SIZE(PFIELD1,3),SIZE(PFIELD1,4)))
!
!*       1.1    Standard case
!               -------------
!
ZFIELD_X(2:IIU  ,:,:,:) =  (ZW_FIELD(1:IIU-1,:,:,:)+ZW_FIELD(2:IIU,:,:,:))  &
                      / MAX(ZW      (1:IIU-1,:,:,:)+ZW      (2:IIU,:,:,:),ZEPS)
ZFIELD_X(1      ,:,:,:) =         ZW_FIELD(1  ,:,:,:)
ZFIELD_X(  IIU+1,:,:,:) =         ZW_FIELD(IIU,:,:,:)
!
!*       1.2    FIELD type value is 0 in the grid mesh
!               --------------------------------------
!
WHERE ( PFIELD1(1:IIU,:,:,:)==0.)
  ZFIELD_X(1:IIU,:,:,:)  = 0.
END WHERE
!
WHERE ( PFIELD1(1:IIU,:,:,:)==0.)
  ZFIELD_X(2:IIU+1,:,:,:)  = 0.
END WHERE
!
!-------------------------------------------------------------------------------
!
!*       2.     FIELD type at grid mesh interfaces (in X directions)
!               ----------------------------------
!
ALLOCATE(ZFIELD_Y(IIU,IJU+1,SIZE(PFIELD1,3),SIZE(PFIELD1,4)))
!
!*       2.1    Standard case
!               -------------
!
ZFIELD_Y(:,2:IJU  ,:,:) =  (ZW_FIELD(:,1:IJU-1,:,:)+ZW_FIELD(:,2:IJU,:,:)) &
                      / MAX(ZW      (:,1:IJU-1,:,:)+ZW      (:,2:IJU,:,:),ZEPS)
ZFIELD_Y(:,1      ,:,:) =       ZW_FIELD(:,1      ,:,:)
ZFIELD_Y(:,  IJU+1,:,:) =       ZW_FIELD(:,  IJU  ,:,:)
!
!*       2.3    FIELD type value is 0 in the grid mesh
!               --------------------------------------
!
WHERE ( PFIELD1(:,1:IJU,:,:)==0.)
  ZFIELD_Y(:,1:IJU,:,:)  = 0.
END WHERE
!
WHERE ( PFIELD1(:,1:IJU,:,:)==0.)
  ZFIELD_Y(:,2:IJU+1,:,:)  = 0.
END WHERE
!
!-------------------------------------------------------------------------------
!
!*       3.     FIELD type at grid mesh corners
!               -------------------------------
!
ALLOCATE(ZFIELD_XY(IIU+1,IJU+1,SIZE(PFIELD1,3),SIZE(PFIELD1,4)))
!
!*       3.1    Standard case
!               -------------
!
ZFIELD_XY(2:IIU  ,2:IJU  ,:,:) =  (  ZW_FIELD(1:IIU-1,1:IJU-1,:,:)   &
                                   + ZW_FIELD(1:IIU-1,2:IJU  ,:,:)   &
                                   + ZW_FIELD(2:IIU  ,1:IJU-1,:,:)   &
                                   + ZW_FIELD(2:IIU  ,2:IJU  ,:,:) ) &
                             / MAX(  ZW      (1:IIU-1,1:IJU-1,:,:)   &
                                   + ZW      (1:IIU-1,2:IJU  ,:,:)   &
                                   + ZW      (2:IIU  ,1:IJU-1,:,:)   &
                                   + ZW      (2:IIU  ,2:IJU  ,:,:) , ZEPS)
!
ZFIELD_XY(1      ,2:IJU  ,:,:) = (  ZW_FIELD(1       ,1:IJU-1,:,:)   &
                                  + ZW_FIELD(1       ,2:IJU  ,:,:) ) &
                            / MAX(  ZW      (1       ,1:IJU-1,:,:)   &
                                  + ZW      (1       ,2:IJU  ,:,:) , ZEPS)
!
ZFIELD_XY(  IIU+1,2:IJU  ,:,:) = (  ZW_FIELD(IIU     ,1:IJU-1,:,:)   &
                                  + ZW_FIELD(IIU     ,2:IJU  ,:,:) ) &
                            / MAX(  ZW      (IIU     ,1:IJU-1,:,:)   &
                                  + ZW      (IIU     ,2:IJU  ,:,:) , ZEPS)
!
ZFIELD_XY(2:IIU  ,1      ,:,:) = (   ZW_FIELD(1:IIU-1,1      ,:,:)   &
                                   + ZW_FIELD(2:IIU  ,1      ,:,:) ) &
                            / MAX(   ZW      (1:IIU-1,1      ,:,:)   &
                                   + ZW      (2:IIU  ,1      ,:,:) , ZEPS)
!
ZFIELD_XY(2:IIU  ,IJU+1  ,:,:) = (   ZW_FIELD(1:IIU-1,IJU    ,:,:)   &
                                   + ZW_FIELD(2:IIU  ,IJU    ,:,:) ) &
                            / MAX(   ZW      (1:IIU-1,IJU    ,:,:)   &
                                   + ZW      (2:IIU  ,IJU    ,:,:) , ZEPS)
!
ZFIELD_XY(1      ,1      ,:,:) =       (   PFIELD1(1      ,1      ,:,:) )
ZFIELD_XY(IIU+1  ,1      ,:,:) =       (   PFIELD1(IIU    ,1      ,:,:) )
ZFIELD_XY(1      ,IJU+1  ,:,:) =       (   PFIELD1(1      ,IJU    ,:,:) )
ZFIELD_XY(IIU+1  ,IJU+1  ,:,:) =       (   PFIELD1(IIU    ,IJU    ,:,:) )
!
!*       3.2    FIELD type value is 0 in one grid mesh
!               --------------------------------------
!
WHERE ( PFIELD1(:    ,1:IJU,:,:)==0.)
  ZFIELD_XY(:    ,1:IJU,:,:)  = 0.
END WHERE
!
WHERE ( PFIELD1(:    ,1:IJU,:,:)==0.)
  ZFIELD_XY(:    ,2:IJU+1,:,:)  = 0.
END WHERE
!
WHERE ( PFIELD1(1:IIU,:    ,:,:)==0.)
  ZFIELD_XY(1:IIU,:    ,:,:)  = 0.
END WHERE
!
WHERE ( PFIELD1(1:IIU,:    ,:,:)==0.)
  ZFIELD_XY(2:IIU+1,:    ,:,:)  = 0.
END WHERE
!
!-------------------------------------------------------------------------------
!
!*       4.     Interpolation coefficients (for left point)
!               --------------------------
!
!*       4.1    X direction
!               -----------
!
DO JEPSX=1,KDXRATIO
  ZCOEF_X(JEPSX) = 1.-2./KDXRATIO*(JEPSX-0.5)
END DO
!
ZCOEF_X(:) = MAX(ZCOEF_X(:),0.)
!
!*       4.2    Y direction
!               -----------
!
DO JEPSY=1,KDYRATIO
  ZCOEF_Y(JEPSY) = 1.-2./KDYRATIO*(JEPSY-0.5)
END DO
!
ZCOEF_Y(:) = MAX(ZCOEF_Y(:),0.)
!
!-------------------------------------------------------------------------------
!
!*       5.     Interpolation
!               -------------
!
PFIELD2(:,:,:,:) = XUNDEF
!
DO JEPSX=1,KDXRATIO
  ZC1_X = ZCOEF_X(           JEPSX)
  ZC3_X = ZCOEF_X(KDXRATIO+1-JEPSX)
  ZC2_X = 1. - ZC1_X - ZC3_X

  DO JEPSY=1,KDYRATIO
    ZC1_Y = ZCOEF_Y(           JEPSY)
    ZC3_Y = ZCOEF_Y(KDYRATIO+1-JEPSY)
    ZC2_Y = 1. - ZC1_Y - ZC3_Y

    DO JI = KXOR,KXEND
      IIS = 1+JPHEXT +JEPSX-1 +(JI-KXOR-1)*KDXRATIO

      DO JJ = KYOR,KYEND
        IJS = 1+JPHEXT +JEPSY-1 +(JJ-KYOR-1)*KDYRATIO
!
        IF (       IIS>=1 .AND. IIS<=SIZE(PFIELD2,1)   &
             .AND. IJS>=1 .AND. IJS<=SIZE(PFIELD2,2) ) THEN
          DO JK=1,SIZE(PFIELD2,3)
            DO JL=1,SIZE(PFIELD2,4)
              IF(PFIELD1(JI,JJ,JK,JL) /= XUNDEF)                               &
              PFIELD2(IIS,IJS,JK,JL) =                                         &
                              ZC1_Y * (   ZC1_X * ZFIELD_XY(JI  ,JJ  ,JK,JL)   &
                                        + ZC2_X * ZFIELD_Y (JI  ,JJ  ,JK,JL)   &
                                        + ZC3_X * ZFIELD_XY(JI+1,JJ  ,JK,JL) ) &
                            + ZC2_Y * (   ZC1_X * ZFIELD_X (JI  ,JJ  ,JK,JL)   &
                                        + ZC2_X * PFIELD1  (JI  ,JJ  ,JK,JL)   &
                                        + ZC3_X * ZFIELD_X (JI+1,JJ  ,JK,JL) ) &
                            + ZC3_Y * (   ZC1_X * ZFIELD_XY(JI  ,JJ+1,JK,JL)   &
                                        + ZC2_X * ZFIELD_Y (JI  ,JJ+1,JK,JL)   &
                                        + ZC3_X * ZFIELD_XY(JI+1,JJ+1,JK,JL) )
            END DO
          END DO
        END IF
!
      END DO
    END DO
  END DO
END DO
!
DEALLOCATE(ZW)
DEALLOCATE(ZW_FIELD)
DEALLOCATE(ZFIELD_X )
DEALLOCATE(ZFIELD_Y )
DEALLOCATE(ZFIELD_XY)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE BILINEAR_4D
!  
!-------------------------------------------------------------------------------
!     #########################################################################
     SUBROUTINE BILINEAR_3D (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                             PFIELD1,PFIELD2                          )
!     #########################################################################
!
!!****  *BILINEAR * - subroutine to interpolate surface FIELD
!!
!!    PURPOSE
!!    -------
!!
!!      This routine interpolates the FIELD mantel from model 1 grid mesh
!!    to model 2 grid mesh.
!!
!!
!!**  METHOD
!!    ------
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
!!      Original     12/12/02 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_BILINEAR_4D
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PFIELD1 ! FIELD on model 1
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PFIELD2 ! FIELD on model 2
!
!
!
!*       0.2    Declarations of local variables
!
! 
REAL, DIMENSION (SIZE(PFIELD1,1),SIZE(PFIELD1,2),SIZE(PFIELD1,3),1) :: ZFIELD1
REAL, DIMENSION (SIZE(PFIELD2,1),SIZE(PFIELD2,2),SIZE(PFIELD2,3),1) :: ZFIELD2
!-------------------------------------------------------------------------------
!
ZFIELD1(:,:,:,1) = PFIELD1(:,:,:)
!
CALL BILINEAR_4D(KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                 ZFIELD1,ZFIELD2                          )
!
PFIELD2(:,:,:) = ZFIELD2(:,:,:,1)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE BILINEAR_3D
!  
!-------------------------------------------------------------------------------
!     #########################################################################
     SUBROUTINE BILINEAR_2D (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                             PFIELD1,PFIELD2                          )
!     #########################################################################
!
!!****  *BILINEAR * - subroutine to interpolate surface FIELD
!!
!!    PURPOSE
!!    -------
!!
!!      This routine interpolates the FIELD mantel from model 1 grid mesh
!!    to model 2 grid mesh.
!!
!!
!!**  METHOD
!!    ------
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
!!      Original     12/12/02 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODI_BILINEAR_4D
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END  
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PFIELD1 ! FIELD on model 1
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD2 ! FIELD on model 2
!
!
!
!*       0.2    Declarations of local variables
!
! 
REAL, DIMENSION (SIZE(PFIELD1,1),SIZE(PFIELD1,2),1,1) :: ZFIELD1
REAL, DIMENSION (SIZE(PFIELD2,1),SIZE(PFIELD2,2),1,1) :: ZFIELD2
!-------------------------------------------------------------------------------
!
ZFIELD1(:,:,1,1) = PFIELD1(:,:)
!
CALL BILINEAR_4D(KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO, &
                 ZFIELD1,ZFIELD2                          )
!
PFIELD2(:,:) = ZFIELD2(:,:,1,1)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE BILINEAR_2D
!  
