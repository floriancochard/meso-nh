!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE PUT_SFX_WAVE(S, U, &
                              KLUOUT,PWAVE_CHA,PWAVE_UCU,PWAVE_VCU,PWAVE_HS,PWAVE_TP)
!     ####################################################
!
!!****  *PUT_SFX_WAVE* - routine to put some variables from
!!                       a wave model
!!    PURPOSE
!!    -------
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
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	J. Pianezze      *LPO*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2014
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,   ONLY : NUNDEF, XUNDEF
USE MODD_SFX_OASIS
!
USE MODI_PACK_SAME_RANK
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
INTEGER,           INTENT(IN)  :: KLUOUT
!
REAL, DIMENSION(:), INTENT(IN) :: PWAVE_CHA
REAL, DIMENSION(:), INTENT(IN) :: PWAVE_UCU
REAL, DIMENSION(:), INTENT(IN) :: PWAVE_VCU
REAL, DIMENSION(:), INTENT(IN) :: PWAVE_HS
REAL, DIMENSION(:), INTENT(IN) :: PWAVE_TP
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
CHARACTER(LEN=50)     :: YCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_WAVE',0,ZHOOK_HANDLE)
!
!*       1.0   Initialization
!              --------------
!
!
!*       2.0   Get variable over wave
!              ---------------------
!
IF(U%NSIZE_SEA>0)THEN
! 
  CALL TREAT_WAVE(U%NSIZE_SEA)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_WAVE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE TREAT_WAVE(KLU)
!
USE MODI_PACK_SAME_RANK
!
IMPLICIT NONE
!
INTEGER,     INTENT(IN) :: KLU
!
CHARACTER(LEN=50)       :: YCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_WAVE:TREAT_WAVE',0,ZHOOK_HANDLE)
!
IF(NWAVE_CHA_ID/=NUNDEF)THEN
  YCOMMENT='Charnock coefficient'
  CALL PACK_SAME_RANK(U%NR_SEA(:),PWAVE_CHA(:),S%XCHARN(:))
  CALL CHECK_WAVE(YCOMMENT,S%XCHARN(:))
ENDIF
!
IF(NWAVE_UCU_ID/=NUNDEF)THEN
  YCOMMENT='u-current velocity'
  CALL PACK_SAME_RANK(U%NR_SEA(:),PWAVE_UCU(:),S%XUMER(:))
  CALL CHECK_WAVE(YCOMMENT,S%XUMER(:))
ENDIF
!
IF(NWAVE_VCU_ID/=NUNDEF)THEN
  YCOMMENT='v-current velocity'
  CALL PACK_SAME_RANK(U%NR_SEA(:),PWAVE_VCU(:),S%XVMER(:))
  CALL CHECK_WAVE(YCOMMENT,S%XVMER(:))
ENDIF
!
IF(NWAVE_HS_ID/=NUNDEF)THEN
  YCOMMENT='Significant wave height'
  CALL PACK_SAME_RANK(U%NR_SEA(:),PWAVE_HS(:),S%XHS(:))
  CALL CHECK_WAVE(YCOMMENT,S%XHS(:))
ENDIF
!
IF(NWAVE_TP_ID/=NUNDEF)THEN
  YCOMMENT='Peak period'
  CALL PACK_SAME_RANK(U%NR_SEA(:),PWAVE_TP(:),S%XTP(:))
  CALL CHECK_WAVE(YCOMMENT,S%XTP(:))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_WAVE:TREAT_WAVE',1,ZHOOK_HANDLE)
!
END SUBROUTINE TREAT_WAVE
!
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_WAVE(HCOMMENT,PFIELD)
!
IMPLICIT NONE
!
CHARACTER(LEN=*),   INTENT(IN) :: HCOMMENT
REAL, DIMENSION(:), INTENT(IN) :: PFIELD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_WAVE:CHECK_WAVE',0,ZHOOK_HANDLE)
!
IF(ANY(PFIELD(:)>=XUNDEF))THEN
  WRITE(KLUOUT,*)'PUT_SFX_WAVE: problem after get '//TRIM(HCOMMENT)//' from OASIS'
  WRITE(KLUOUT,*)'PUT_SFX_WAVE: some points not defined = ',COUNT(PFIELD(:)>=XUNDEF)
  CALL ABOR1_SFX('PUT_SFX_WAVE: problem after get '//TRIM(HCOMMENT)//' from OASIS')
ENDIF
!
IF (LHOOK) CALL DR_HOOK('PUT_SFX_WAVE:CHECK_WAVE',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_WAVE
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PUT_SFX_WAVE
