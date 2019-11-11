!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE INIT_SLT (SLT, &
                     HPROGRAM  &! Program calling unit
       )  
!     ######################################################################
!
!!****  *INIT_SLT* - 
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!    !!
!!    EXTERNAL
!!    --------
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
!!    ?
!!
!!    MODIFICATIONS
!!    -------------
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!! 
!------------------------------------------------------------------------------
!
!
USE MODD_SLT_n, ONLY : SLT_t
!
USE MODD_SLT_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!PASSED VARIABLES
!
TYPE(SLT_t), INTENT(INOUT) :: SLT
!
 CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM              !Passing unit
!
!LOCAL VARIABLES
 CHARACTER(LEN=4)    :: CRGUNITS              ! type of log-normal geometric mean radius
INTEGER             :: JMODE                 ! Counter for sea salt modes
INTEGER             :: JMODE_IDX             ! Index for sea salt modes
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!get output listing unit
IF (LHOOK) CALL DR_HOOK('INIT_SLT',0,ZHOOK_HANDLE)
!
!Allocate memory for the real values which will be used by the model
ALLOCATE(SLT%XEMISRADIUS_SLT(NSLTMDE))
ALLOCATE(SLT%XEMISSIG_SLT   (NSLTMDE))
!
!Get initial size distributions. This is cut and pasted
!from dead routine dstpsd.F90
!Check for different source parameterizations

IF (CEMISPARAM_SLT.eq."Ova14") THEN
  NSLTMDE = 5
!  JORDER_SLT = (/3,2,1,4,5/) !Salt modes in order of importance
  CRGUNITS   = 'NUMB'
  XEMISRADIUS_INI_SLT = (/0.009, 0.021, 0.045, 0.115, 0.415/)
  XEMISSIG_INI_SLT = (/1.37, 1.5, 1.42, 1.53, 1.85/)

ELSE IF(CEMISPARAM_SLT.eq."Vig01") THEN
   NSLTMDE = 5
!  JORDER_SLT = (/3,2,1,4,5/) !Salt modes in order of importance, only three modes
   CRGUNITS   = 'NUMB'
   XEMISRADIUS_INI_SLT =  (/ 0.2,  2.0, 12.,0.,0. /)         ! [um]  Number median radius Viganati et al., 2001
   XEMISSIG_INI_SLT = (/ 1.9, 2.0, 3.00,0.,0.  /)  ! [frc] Geometric standard deviation Viganati et al., 2001

ELSE IF(CEMISPARAM_SLT.eq."Sch04") THEN ! use default of Schultz et al, 2004
   NSLTMDE = 5
!  JORDER_SLT = (/3,2,1,4,5/), only three modes
   CRGUNITS   = 'MASS'
   XEMISRADIUS_INI_SLT = 0.5*(/0.28, 2.25, 15.32, 0., 0./)! [um] Mass median radius
   XEMISSIG_INI_SLT = (/1.59, 2.00, 2.00, 0., 0./) ! [frc] Geometric standard deviation

ENDIF
!
DO JMODE=1,NSLTMDE
  !
  JMODE_IDX = JORDER_SLT(JMODE)
  !
  SLT%XEMISSIG_SLT   (JMODE) = XEMISSIG_INI_SLT   (JMODE_IDX)
  SLT%XEMISRADIUS_SLT(JMODE) = XEMISRADIUS_INI_SLT(JMODE_IDX)
  !
  IF (CRGUNITS=="MASS") &
    SLT%XEMISRADIUS_SLT(JMODE) = SLT%XEMISRADIUS_SLT(JMODE) * EXP(-3.d0 * (LOG(SLT%XEMISSIG_SLT(JMODE)))**2)    
  !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('INIT_SLT',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_SLT
