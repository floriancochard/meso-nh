!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      SUBROUTINE SETLB_LG
!     ###################
!
!!****  *SETLB_LG* - routine to set lagragian variable values
!!                  at the LB for outer model (model1)
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set lagragian variable values
!   at the LB for outer model (model1) in the case of constant advective
!   velocity
!
!!**  METHOD
!!    ------
!!    All required values are provided by modeln1 modules.
!     Floating indices are used to set only LBSV variables corresponding to
!     lagrangian variables
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_NSV : NSV_LGBED, NSV_LGEND
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_LG)
!!
!!
!!    AUTHOR
!!    ------
!!  P. Jabouille / J Stein      * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        22/06/01
!!     28/03/2018 P. Wautelet: replace TEMPORAL_DIST by DATETIME_DISTANCE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_TIME
USE MODD_TIME_n
USE MODD_LBC_n
USE MODD_LSFIELD_n
USE MODD_GRID_n
USE MODD_LG
USE MODD_NSV,         ONLY : NSV_LGBEG,NSV_LGEND
!
USE MODE_DATETIME
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.2   declarations of local variables
!
INTEGER      :: IIU,IJU,IKU ! Upper bounds
INTEGER      :: JI,JJ,JK    ! loop index
REAL         :: ZTEMP_DIST  ! time from the begining of simulation
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE
!              --------
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=SIZE(XZZ,3)
!
CALL DATETIME_DISTANCE(TDTEXP,TDTCUR,ZTEMP_DIST)
!
IF ( CLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
  DO JK=1,IKU
    DO JJ=1,IJU
      XLBXSVM(1,JJ,JK,NSV_LGBEG)=0.5*(XXHAT(1)+XXHAT(2))-ZTEMP_DIST*XLGSPEED
    END DO
  END DO
!
  DO JK=1,IKU
    DO JJ=1,IJU-1
      XLBXSVM(1,JJ,JK,NSV_LGBEG+1)=0.5*(XYHAT(JJ)+XYHAT(JJ+1))
    END DO
    XLBXSVM(1,IJU,JK,NSV_LGBEG+1)=1.5*XYHAT(IJU)-0.5*XYHAT(IJU-1)
  END DO
  !
  DO JJ=1,IJU
    DO JK=1,IKU-1
      XLBXSVM(1,JJ,JK,NSV_LGEND)=0.5*(XZZ(1,JJ,JK)+XZZ(1,JJ,JK+1))
    END DO
    XLBXSVM(1,JJ,IKU,NSV_LGEND)=1.5*XZZ(1,JJ,IKU)  -0.5*XZZ(1,JJ,IKU-1)
  END DO
END IF
  !
IF ( CLBCX(1) /= "CYCL" .AND. LEAST_ll()) THEN
  DO JK=1,IKU
    DO JJ=1,IJU
      XLBXSVM(SIZE(XLBXSVM,1),JJ,JK,NSV_LGBEG)=0.5*(XXHAT(IIU-1)+XXHAT(IIU))+ZTEMP_DIST*XLGSPEED
    END DO
  END DO
!
  DO JK=1,IKU
    DO JJ=1,IJU-1
      XLBXSVM(SIZE(XLBXSVM,1),JJ,JK,NSV_LGBEG+1)=0.5*(XYHAT(JJ)+XYHAT(JJ+1))
    END DO
    XLBXSVM(SIZE(XLBXSVM,1),IJU,JK,NSV_LGBEG+1)=1.5*XYHAT(IJU)-0.5*XYHAT(IJU-1)
  END DO
!
  DO JJ=1,IJU
    DO JK=1,IKU-1
      XLBXSVM(SIZE(XLBXSVM,1),JJ,JK,NSV_LGEND)=0.5*(XZZ(IIU,JJ,JK)+XZZ(IIU,JJ,JK+1))
    END DO
    XLBXSVM(SIZE(XLBXSVM,1),JJ,IKU,NSV_LGEND)=1.5*XZZ(IIU,JJ,IKU)-0.5*XZZ(IIU,JJ,IKU-1)
  END DO
  !
ENDIF
  !
IF (SIZE(XLBYSVM,1) .NE. 0 .AND. CLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
!
  DO JK=1,IKU
    DO JI=1,IIU-1
      XLBYSVM(JI,1,JK,NSV_LGBEG)=0.5*(XXHAT(JI)+XXHAT(JI+1))
    END DO
    XLBYSVM(IIU,1,JK,NSV_LGBEG)=1.5*XXHAT(IIU)-0.5*XXHAT(IIU-1)
  END DO
      !
  DO JK=1,IKU
    DO JI=1,IIU
      XLBYSVM(JI,1,JK,NSV_LGBEG+1)=0.5*(XYHAT(1)+XYHAT(2))-ZTEMP_DIST*XLGSPEED
    END DO
  END DO
      !
  DO JI=1,IIU
    DO JK=1,IKU-1
      XLBYSVM(JI,1,JK,NSV_LGEND)=0.5*(XZZ(JI,1,JK)+XZZ(JI,1,JK+1))
    END DO
    XLBYSVM(JI,1,IKU,NSV_LGEND)=1.5*XZZ(JI,1,IKU)  -0.5*XZZ(JI,1,IKU-1)
  END DO

ENDIF
    !
IF (SIZE(XLBYSVM,1) .NE. 0 .AND. CLBCY(1) /= "CYCL" .AND. LNORTH_ll()) THEN
!
  DO JK=1,IKU
    DO JI=1,IIU-1
      XLBYSVM(JI,SIZE(XLBYSVM,2),JK,NSV_LGBEG)=0.5*(XXHAT(JI)+XXHAT(JI+1))
    END DO
    XLBYSVM(IIU,SIZE(XLBYSVM,2),JK,NSV_LGBEG)=1.5*XXHAT(IIU)-0.5*XXHAT(IIU-1)
  END DO
!
  DO JK=1,IKU
    DO JI=1,IIU
      XLBYSVM(JI,SIZE(XLBYSVM,2),JK,NSV_LGBEG+1)=0.5*(XYHAT(IJU-1)+XYHAT(IJU))+ZTEMP_DIST*XLGSPEED
    END DO
  END DO
!
  DO JI=1,IIU
    DO JK=1,IKU-1
      XLBYSVM(JI,SIZE(XLBYSVM,2),JK,NSV_LGEND)=0.5*(XZZ(JI,IJU,JK)+XZZ(JI,IJU,JK+1))
    END DO
    XLBYSVM(JI,SIZE(XLBYSVM,2),IKU,NSV_LGEND)=1.5*XZZ(JI,IJU,IKU)-0.5*XZZ(JI,IJU,IKU-1)
  END DO
ENDIF
!
END SUBROUTINE SETLB_LG
