!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ###########################
      MODULE MODI_CH_INIT_MODEL0D
!!    ###########################
!!
INTERFACE
SUBROUTINE CH_INIT_MODEL0D(HNAMELISTFILE)
IMPLICIT NONE
CHARACTER*(*), INTENT(IN) :: HNAMELISTFILE
END SUBROUTINE CH_INIT_MODEL0D
END INTERFACE
END MODULE MODI_CH_INIT_MODEL0D
!!
!!    #########################################
      SUBROUTINE CH_INIT_MODEL0D(HNAMELISTFILE)
!!    #########################################
!!
!!***  *CH_INIT_MODEL0D*
!!
!!    PURPOSE
!!    -------
!     initialisation of the different modules from the namelist file
!     and of the time control variables
!!
!!**  METHOD
!!    ------
!!    simple
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/03/95
!!    27/07/96 (K. Suhre) restructured
!!    19/03/96 (K. Suhre) insert namelists directly for compiler reasons
!!    31/08/96 (K. SUHRE) read MODD_BLANK namelist
!!    01/08/01 (C. Mari)  change NAM_CH_SOLVER to $n 
!!    Philippe Wautelet: 10/01/2019: use newunit argument to open files
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D    ! use whole module for namelist input
!USE MODD_CH_MNHC_n, ONLY: CCHEM_INPUT_FILE
USE MODD_CH_AEROSOL
USE MODD_CH_AEROSOL0D, ONLY : XEMISSIGI, XEMISSIGJ, XEMISRADIUSI, &
                              XEMISRADIUSJ
                       ! for using the dry deposition scheme
                       ! in the box model we need to use some MesoNH
                       ! variables
USE MODD_CH_JVALUES_n
USE MODN_CH_SOLVER_n
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE

CHARACTER*(*), INTENT(IN) :: HNAMELISTFILE ! name of namelist input file

INTEGER :: ILU ! unit number for IO

NAMELIST /NAM_CH_MODEL0D/ XTBEGIN, XTEND, XDTACT,          &
                       XDTOUT, XDTDIAG,                 &
                       CRUNID,                          &
                       CINITFILE, COUTFILE, CMETEOFILE, &
                       CRESULTFILE, CRESULTFORMAT,      &
                       CDIAGFILE, CDIAGFORMAT,          &
                       NVERB,                           &
                       LCH_TUV_ONLINE,                  &
                       CCH_TUV_LOOKUP,                  &
                       CCH_TUV_CLOUDS,                  &
                       XCH_TUV_ALBNEW,                  &
                       XCH_TUV_DOBNEW,                  &
                       XCH_TUV_TUPDATE,                 &
                       LCH_SURFACE0D, LORILAM,          &
                       XN0IMIN, XN0JMIN, LSEDIMAERO,    &
                       LHETEROSO4, CNUCLEATION, XINISIGI, XINISIGJ,  & 
                       XINIRADIUSI, XINIRADIUSJ,  LVARSIGI,&
                       LVARSIGJ, XEMISSIGI, XEMISSIGJ,  &
                       XEMISRADIUSI, XEMISRADIUSJ, CMINERAL, CORGANIC,&
                       XSIGIMIN, XSIGIMAX,XSIGJMIN, XSIGJMAX,  & 
                       XCOEFRADIMAX, XCOEFRADIMIN, XCOEFRADJMAX, XCOEFRADJMIN,&
                       CRGUNIT,  CCHEM_INPUT_FILE
!
!------------------------------------------------------------------------------
!
!*       1.   READ NAMELIST FILE
!        -----------------------
!

CALL INIT_NAM_CH_SOLVERn


PRINT *, 'CH_INIT_MODEL0D: opening namelist file: ', HNAMELISTFILE
OPEN(NEWUNIT=ILU,FILE=HNAMELISTFILE,STATUS="OLD",FORM="FORMATTED")
READ(UNIT=ILU,NML=NAM_CH_MODEL0D)
READ(UNIT=ILU,NML=NAM_CH_SOLVERn)
CLOSE(UNIT=ILU)

! quick fix for use of the surface scheme in the box model
CCHEM_INPUT_FILE = 'CHCONTROL1.nam'

IF (NVERB >= 5) THEN
  PRINT *, '==================================================================='
  PRINT *, 'CH_INIT_MODEL0D: namelist NAM_CH_MODEL0D set to:'
  WRITE(*,NAM_CH_MODEL0D)
  PRINT *, '==================================================================='
  PRINT *, 'CH_INIT_MODEL0D: namelist NAM_CH_SOLVERn set to:'
    WRITE(*,NAM_CH_SOLVERn)
  PRINT *, '==================================================================='
!  PRINT *, 'CH_INIT_MODEL0D: namelist NAM_BLANK set to:'
!    WRITE(*,NAM_BLANK)
  PRINT *, '==================================================================='
  PRINT *, 'RUNID:', CRUNID
  PRINT *, '==================================================================='
END IF
!
!*       2.   SET TIME CONTROL VARIABLES
!        -------------------------------
!
! set XTSIMUL to start time
XTSIMUL = XTBEGIN

! set time for next write to disk
XTNEXTOUT = XTBEGIN + XDTOUT

! set time for next diagnostics
XTNEXTDIAG = XTBEGIN + XDTDIAG
!
!*       3.   PRINT SOME INFORMATION
!        ---------------------------
!
IF (NVERB >= 5) THEN
  PRINT 9, "CH_INIT_MODEL0D: start of simulation at ", XTBEGIN
  PRINT 9, "CH_INIT_MODEL0D: end of simulation at   ", XTEND
  PRINT 9, "CH_INIT_MODEL0D: the time-step is       ", XDTACT
  PRINT 9, "CH_INIT_MODEL0D: number of steps is     ", (XTEND-XTBEGIN)/XDTACT
9 FORMAT(A,F15.0)
END IF
!

CALL UPDATE_NAM_CH_SOLVERn
!
END SUBROUTINE CH_INIT_MODEL0D
