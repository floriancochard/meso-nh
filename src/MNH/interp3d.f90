!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
MODULE MODI_INTERP3D
!#################################
!
INTERFACE
      SUBROUTINE INTERP3D(PFIELD,KGRID,PSVAL,PPLEV,PFIELDAP)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD      ! values of the field
INTEGER,                INTENT(IN) :: KGRID       ! Mesonh grid indicator
REAL,                   INTENT(IN) :: PSVAL       ! value for missing data
REAL, DIMENSION(:),     INTENT(IN) :: PPLEV       ! list of vertical levels
REAL, DIMENSION(:,:,:), INTENT(OUT):: PFIELDAP    ! values of the field on the pressure levels
END SUBROUTINE  INTERP3D
END INTERFACE
END MODULE MODI_INTERP3D
!     ######spl
      SUBROUTINE INTERP3D(PFIELD,KGRID,PSVAL,PPLEV,PFIELDAP)
!     #####################
!
!!****  *INTERP3D* - interpole 3D fields on pressure levels
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      Functions MXF, MYF, MZF
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_FIELD1    : contains prognostics  variables 
!!         XPASBM
!!      Module MODD_GRID1
!!         XZZ
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Ducrocq  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/03/97
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
! 
USE MODD_PARAMETERS
USE MODD_DIM_n
USE MODD_FIELD_n
USE MODD_GRID_n
!
USE MODI_SHUMAN ! interface modules 
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD      ! values of the field
INTEGER,                INTENT(IN) :: KGRID       ! Mesonh grid indicator
REAL,                   INTENT(IN) :: PSVAL       ! value for missing data
REAL, DIMENSION(:),     INTENT(IN) :: PPLEV       ! list of vertical levels
REAL, DIMENSION(:,:,:), INTENT(OUT):: PFIELDAP    ! values of the field on the pressure levels
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: JKP,JKLOOP,JJLOOP,JILOOP,IJ,II ! loop indices
INTEGER      :: IIE,IJE,IPU      ! End of usefull area 
INTEGER      :: IIB,IJB,IKB      ! Begining of usefull area 
REAL, DIMENSION(SIZE(XPABST,1),SIZE(XPABST,2),SIZE(XPABST,3)) :: ZPTH  ! pressure for grid points corresponding to KGRID type 
REAL  :: ZREF,ZXP,ZXM,ZDIXEPS ! pressure values and epsilon value
INTEGER :: IKU
!-------------------------------------------------------------------------------
!
!*         1.    
!                 ------------
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IPU=SIZE(PFIELDAP,3)
IKB=1 +JPVEXT
IKU=SIZE(XZHAT)
ZDIXEPS=10.*EPSILON(1.)
!
SELECT CASE (KGRID)
  CASE(1)
    ZPTH=XPABST
  CASE(2)
    ZPTH(:,:,:)=MXM(XPABST(:,:,:))
    ZPTH(1,:,:)=2.*ZPTH(2,:,:) - ZPTH(3,:,:)
  CASE(3)
    ZPTH(:,:,:)=MYM(XPABST(:,:,:))
    ZPTH(:,1,:)=2.*ZPTH(:,2,:) - ZPTH(:,3,:)
  CASE(4)
    ZPTH(:,:,:)=MZM(1,IKU,1,XPABST(:,:,:))
    ZPTH(:,:,1)=2.*ZPTH(:,:,2) - ZPTH(:,:,3)
END SELECT
!
DO JKP= 1, IPU
   ZREF=ALOG10(PPLEV(JKP))
   DO JILOOP = IIB,IIE
      DO JJLOOP = IJB,IJE
         IJ=JJLOOP-IJB+1
         II=JILOOP-IIB+1
         PFIELDAP(II,IJ,JKP)=PSVAL
         DO JKLOOP = 1,NKMAX+2*JPVEXT
            ZXM=ALOG10(ZPTH(JILOOP,JJLOOP,JKLOOP))
            ZXP=ALOG10(ZPTH(JILOOP,JJLOOP,MIN(NKMAX+2*JPVEXT,JKLOOP+1)))
            IF ((ZXP-ZREF)*(ZREF-ZXM) .GE.0.) THEN
               IF (JKLOOP+1 == IKB) THEN
                  CYCLE
               ELSE
                  GO TO 4
               ENDIF
            ELSE IF (ZXP.GE.ZXM-ZDIXEPS.AND.ZXP.LE.ZXM+ZDIXEPS.AND.  &
      ZREF.GE.ZXM-ZDIXEPS.AND.ZREF.LE.ZXM+ZDIXEPS) THEN
               IF(JKLOOP+1 == IKB)THEN
                  CYCLE
               ELSE
                  GO TO 4
               ENDIF
            END IF
         END DO
         GO TO 3
4     CONTINUE
!
!  We interpolate 
         PFIELDAP(II,IJ,JKP)= (PFIELD(II,IJ,JKLOOP)* (ZXP-ZREF)+ &
              PFIELD(II,IJ,MIN(NKMAX+2*JPVEXT,JKLOOP+1))* (ZREF-ZXM)) &
              / MIN(-1.E-08,(ZXP-ZXM))
         GO TO 3
3     CONTINUE
      END DO
   END DO
END DO
!
END SUBROUTINE INTERP3D
