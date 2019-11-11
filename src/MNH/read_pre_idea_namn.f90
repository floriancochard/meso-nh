!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/06/06 11:10:18
!-----------------------------------------------------------------
!####################################
MODULE MODI_READ_PRE_IDEA_NAM_n
!####################################
INTERFACE 
  SUBROUTINE READ_PRE_IDEA_NAM_n(KLUPRE,KLUOUT)
    INTEGER,INTENT(IN) :: KLUPRE,KLUOUT ! Logical unit numbers for EXPRE file
                                       ! and for output_listing file
  END SUBROUTINE READ_PRE_IDEA_NAM_n
END INTERFACE
!
END MODULE MODI_READ_PRE_IDEA_NAM_n
!
!
!     #############################################
      SUBROUTINE READ_PRE_IDEA_NAM_n(KLUPRE,KLUOUT)
!     #############################################

USE MODD_DIM_n, ONLY: NIMAX_n=>NIMAX, NJMAX_n=>NJMAX
USE MODD_DYN_n, ONLY: CPRESOPT_n=>CPRESOPT,NITR_n=>NITR,XRELAX_n=>XRELAX, &
                      LRES_n=>LRES,XRES_n=>XRES
USE MODD_LBC_n, ONLY: CLBCX_n=>CLBCX,CLBCY_n=>CLBCY
USE MODD_CH_MNHC_n, ONLY:  LCH_INIT_FIELD_n=>LCH_INIT_FIELD, &
                           CCHEM_INPUT_FILE_n=>CCHEM_INPUT_FILE
USE MODD_PARAM_n, ONLY: CSURF_n=>CSURF  
USE MODD_CH_AEROSOL, ONLY : LORILAM
USE MODN_LUNIT_n           ! Namelist modules
USE MODN_CONF_n
USE MODE_POS
!
IMPLICIT NONE
INTEGER,INTENT(IN) :: KLUPRE,KLUOUT ! Logical unit numbers for EXPRE file
                                       ! and for output_listing file
!
LOGICAL :: GFOUND                  ! Return code when searching namelist 
! Namelist variables from $n
INTEGER                :: NIMAX,NJMAX
CHARACTER(LEN=5)       :: CPRESOPT
INTEGER                :: NITR
LOGICAL                :: LRES
REAL                   :: XRES
REAL                   :: XRELAX      
CHARACTER(LEN=4), DIMENSION(2) :: CLBCX,CLBCY
CHARACTER(LEN=5)       :: CSURF
LOGICAL                :: LCH_INIT_FIELD
CHARACTER(LEN=80)      :: CCHEM_INPUT_FILE
!
!*       0.2  Namelist declarations
!
NAMELIST/NAM_DIMn_PRE/ NIMAX,NJMAX          ! Declaration in MODD_DIM1
NAMELIST/NAM_DYNn_PRE/ CPRESOPT,NITR,XRELAX,LRES,XRES ! Declaration in MODD_DYN1
NAMELIST/NAM_LBCn_PRE/ CLBCX,CLBCY          ! Declaration in MODD_LBC1
NAMELIST/NAM_GRn_PRE/  CSURF                ! Declaration in MODD_PARAM1
NAMELIST/NAM_CH_MNHCn_PRE/ LCH_INIT_FIELD, CCHEM_INPUT_FILE, LORILAM
!
!------------------------------------------------------------------------------
!
CALL INIT_NMLVAR
CALL POSNAM(KLUPRE,'NAM_DIMN_PRE',GFOUND,KLUOUT)
IF (GFOUND) READ(UNIT=KLUPRE,NML=NAM_DIMn_PRE)
CALL POSNAM(KLUPRE,'NAM_DYNN_PRE',GFOUND,KLUOUT)
IF (GFOUND) READ(UNIT=KLUPRE,NML=NAM_DYNn_PRE)
CALL POSNAM(KLUPRE,'NAM_LBCN_PRE',GFOUND,KLUOUT)
IF (GFOUND) READ(UNIT=KLUPRE,NML=NAM_LBCn_PRE)
CALL POSNAM(KLUPRE,'NAM_GRN_PRE',GFOUND,KLUOUT)
IF (GFOUND) READ(UNIT=KLUPRE,NML=NAM_GRn_PRE)
CALL POSNAM(KLUPRE,'NAM_CH_MNHCN_PRE',GFOUND,KLUOUT)
IF (GFOUND) READ(UNIT=KLUPRE,NML=NAM_CH_MNHCn_PRE)
CALL UPDATE_MODD_FROM_NMLVAR
!
CALL POSNAM(KLUPRE,'NAM_CONFN',GFOUND,KLUOUT)
IF (GFOUND) THEN 
  CALL INIT_NAM_CONFn
  READ(UNIT=KLUPRE,NML=NAM_CONFn)
  CALL UPDATE_NAM_CONFn
END IF
CALL POSNAM(KLUPRE,'NAM_LUNITN',GFOUND,KLUOUT)
IF (GFOUND) THEN 
  CALL INIT_NAM_LUNITn 
  READ(UNIT=KLUPRE,NML=NAM_LUNITn)
  CALL UPDATE_NAM_LUNITn
END IF

!------------------------------------------------------------------------------
CONTAINS

SUBROUTINE INIT_NMLVAR

NIMAX=NIMAX_n
NJMAX=NJMAX_n
CPRESOPT=CPRESOPT_n
LRES=LRES_n
XRES=XRES_n
NITR=NITR_n
XRELAX=XRELAX_n
CLBCX=CLBCX_n
CLBCY=CLBCY_n
CSURF=CSURF_n
LCH_INIT_FIELD=LCH_INIT_FIELD_n
CCHEM_INPUT_FILE=CCHEM_INPUT_FILE_n
END SUBROUTINE INIT_NMLVAR

SUBROUTINE UPDATE_MODD_FROM_NMLVAR

NIMAX_n=NIMAX
NJMAX_n=NJMAX
CPRESOPT_n=CPRESOPT
LRES_n=LRES
XRES_n=XRES
NITR_n=NITR
XRELAX_n=XRELAX
CLBCX_n=CLBCX
CLBCY_n=CLBCY
CSURF_n=CSURF
LCH_INIT_FIELD_n=LCH_INIT_FIELD
CCHEM_INPUT_FILE_n=CCHEM_INPUT_FILE
END SUBROUTINE UPDATE_MODD_FROM_NMLVAR

END SUBROUTINE READ_PRE_IDEA_NAM_n
