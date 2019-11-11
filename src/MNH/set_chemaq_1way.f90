!MNH_LIC Copyright 2009-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###########################
       MODULE MODI_SET_CHEMAQ_1WAY
!      ###########################
!
INTERFACE
!
      SUBROUTINE SET_CHEMAQ_1WAY (PRHODREF,PSVTOLD,PSVTNEW)
!
REAL, DIMENSION(:,:,:),    INTENT(IN) :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:,:),  INTENT(IN) :: PSVTOLD    ! chemical concentrations dad model
!
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT):: PSVTNEW  ! chemical concentrations whole array
!
!
END SUBROUTINE SET_CHEMAQ_1WAY
!
END INTERFACE
!
END MODULE MODI_SET_CHEMAQ_1WAY
!
!     ###########################################################
      SUBROUTINE SET_CHEMAQ_1WAY (PRHODREF,PSVTOLD,PSVTNEW)
!     ###########################################################
!
!!****  *SET_CHEMAQ_1WAY * - transfer chemical concentrations in gas phase from
!!                 the dad model in an array with the same dimension than the
!!                 child one, which includes aqueous phase concentrations
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to transfer chemical concentrations in 
!!    gas phase from the dad model in an array with the same dimension than the
!!    child one, which includes aqueous phase concentrations. Aqueous phase
!!    concentrations and gas-phase concentrations for species not considered
!!    in the chemical mechanism used by the dad model are set to zero.
!!      This routine is used to compute the LB tendencies in ONE_WAY$n in case
!!    of grid-nesting when the dad model uses only gas-phase chemistry and the
!!    child uses gas phase and aqueous phase chemistry.
!!
!!    METHOD
!!    ------
!!      The dimension of the old and new xsvt arrays are compared and extra
!!    elements in the new xsvt array are set to zero.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_RAIN_C2R2_DESCR, ONLY : XRTMIN, XCTMIN
!!      Module MODD_RAIN_C2R2_PARAM, ONLY : XCONCC_INI, XCONCR_PARAM_INI
!!      Module MODD_CONF,            ONLY : NVERB
!!
!!    REFERENCE
!!    ---------
!!      Not available
!!
!!    AUTHOR
!!    ------
!!      M. Leriche      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/11/09
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_LUNIT_n, ONLY: TLUOUT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:),    INTENT(IN) :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:,:),  INTENT(IN) :: PSVTOLD    ! chemical concentrations dad model
!
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT):: PSVTNEW  ! chemical concentrations whole array
!
!
!*       0.2   Declarations of local variables :
!
INTEGER    :: IRESP   ! Return code of FM routines
INTEGER    :: J       ! Loop counter
INTEGER    :: ILUOUT  ! Logical unit number of output-listing
!
!-------------------------------------------------------------------------------
!*       1.    RETRIEVE LOGICAL UNIT NUMBER
!              ----------------------------
!
ILUOUT = TLUOUT%NLU
!
!*       2.    TRANSFER DATA FROM OLD TO NEW ARRAY
!              -----------------------------------
!
PSVTNEW(:,:,:,:) = 0.
!
DO J = 1, SIZE(PSVTOLD,4)
  PSVTNEW(:,:,:,J) = PSVTOLD(:,:,:,J)
ENDDO

DO J = SIZE(PSVTOLD,4)+1,SIZE(PSVTNEW,4)
  PSVTNEW(:,:,:,J) = 1.E-36
ENDDO

IF (SIZE(PSVTOLD,4).GE.SIZE(PSVTNEW,4)) THEN
  WRITE (UNIT=ILUOUT,FMT=*) "set_chemaq_1way problem the size of chemical"
  WRITE (UNIT=ILUOUT,FMT=*) "concentrations array from dad model is superior"
  WRITE (UNIT=ILUOUT,FMT=*) "than the size of array from the child model"
ENDIF
!
END SUBROUTINE SET_CHEMAQ_1WAY
