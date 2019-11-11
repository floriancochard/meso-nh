!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE GET_SFX_WAVE(U, DGS, &
                              PWAVE_U10,PWAVE_V10)  
!     ############################################################################
!
!!****  *GET_SFX_WAVE* - routine to get some variables from surfex to
!                        a wave model
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
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DIAG_n, ONLY : DIAG_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DIAG_t), INTENT(INOUT) :: DGS
!
REAL, DIMENSION(:), INTENT(OUT) :: PWAVE_U10  ! 10 meter u-wind (m/s)
REAL, DIMENSION(:), INTENT(OUT) :: PWAVE_V10  ! 10 meter v-wind (m/s)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_WAVE',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.0   Initialization
!              --------------
!
PWAVE_U10(:) = XUNDEF
PWAVE_V10(:) = XUNDEF
!
!*       2.0   Get variable over wave
!              ----------------------
!
IF(U%NSIZE_SEA>0)THEN
!
  CALL UNPACK_SAME_RANK(U%NR_SEA(:),DGS%XZON10M(:),PWAVE_U10(:),XUNDEF)
  CALL UNPACK_SAME_RANK(U%NR_SEA(:),DGS%XMER10M(:),PWAVE_V10(:),XUNDEF)
!
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_WAVE',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_SFX_WAVE
