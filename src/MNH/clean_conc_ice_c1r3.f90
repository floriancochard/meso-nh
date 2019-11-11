!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###############################
       MODULE MODI_CLEAN_CONC_ICE_C1R3
!      ###############################
!
INTERFACE
!
      SUBROUTINE CLEAN_CONC_ICE_C1R3 (PR,PSV)
!
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PR     ! microphysical mixing 
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT) :: PSV    ! microphys. concentrations
!
END SUBROUTINE CLEAN_CONC_ICE_C1R3
!
END INTERFACE
!
END MODULE MODI_CLEAN_CONC_ICE_C1R3
!
!     ##############################################
      SUBROUTINE CLEAN_CONC_ICE_C1R3 (PR,PSV)
!     ##############################################
!
!!****  *CLEAN_CONC_ICE_C1R3 * - reinitialize the ice crystal
!!                             concentration before STARTing the C1R3 scheme
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to reinitialize cloud ativated CCN, the
!!    droplet and the rain drop concentrations when the cloud droplet and rain
!!    drop mixing ratios have very low values and when these concentrations
!!    have negative values (after a spawning or prep_real_case operation).
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_ICE_C1R3_DESCR, ONLY : XRTMIN, XCTMIN
!!      Module MODD_CONF,           ONLY : NVERB
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine CLEAN_CONC_ICE_C1R3 )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      P. Jabouille     * CNRM/GMME *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/04/01
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,           ONLY : NVERB
USE MODD_ICE_C1R3_DESCR, ONLY : XRTMIN, XCTMIN
USE MODD_LUNIT_n,        ONLY : TLUOUT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PR     ! microphysical mixing 
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT) :: PSV    ! microphys. concentrations
!
!
!*       0.2   Declarations of local variables :
!
INTEGER    :: ILUOUT  ! Logical unit number of output-listing
!  
!-------------------------------------------------------------------------------
!  
!*       1.    RETRIEVE LOGICAL UNIT NUMBER
!              ----------------------------
!
ILUOUT = TLUOUT%NLU
!
!*       2.    INITIALIZATION
!              --------------
!
PSV(:,:,:,4:5) = MAX( PSV(:,:,:,4:5),0.0 )
WHERE (PR(:,:,:,4) <= XRTMIN(4) .OR. PSV(:,:,:,4) <= XCTMIN(4))
  PSV(:,:,:,5) = 0.0
  PSV(:,:,:,4) = 0.0
  PR (:,:,:,4) = 0.0
ENDWHERE
!
IF( NVERB >= 5 ) THEN
  WRITE (UNIT=ILUOUT,FMT=*) "!INI_MODEL$n: The ice crystal concentration is"
  WRITE (UNIT=ILUOUT,FMT=*) "cleaned to avoid very low and negative values"
END IF
!
END SUBROUTINE CLEAN_CONC_ICE_C1R3
