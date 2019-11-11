!MNH_LIC Copyright 2000-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #######################################
       MODULE MODI_SET_CONC_RAIN_C2R2
!      #######################################
!
INTERFACE
!
      SUBROUTINE SET_CONC_RAIN_C2R2 (HGETCLOUD,PRHODREF,PRT,PSVT)
!
CHARACTER (LEN=4),         INTENT(IN) :: HGETCLOUD  ! Get indicator
REAL, DIMENSION(:,:,:),    INTENT(IN) :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PRT     ! microphysical mixing ratios
!
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT):: PSVT     ! microphys. concentrations
!
!
END SUBROUTINE SET_CONC_RAIN_C2R2
!
END INTERFACE
!
END MODULE MODI_SET_CONC_RAIN_C2R2
!
!     ###########################################################
      SUBROUTINE SET_CONC_RAIN_C2R2 (HGETCLOUD,PRHODREF,PRT,PSVT)
!     ###########################################################
!
!!****  *SET_CONC_RAIN_C2R2 * - initialize droplet and raindrop
!!                   concentration for a RESTArt simulation of the C2R2 scheme
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize cloud droplet and rain drop
!!    concentrations when the cloud droplet and rain drop mixing ratios are
!!    only available (generally from a previous run using the Kessler scheme).
!!      This routine is used to initialize the droplet/drop concentrations
!!    using the r_c and r_r of a previous REVE or KESS run but also to compute
!!    the LB tendencies in ONE_WAY$n in case of grid-nesting when the optional
!!    argument PTIME is set (a C2R2 run embedded in a KESS or REVE run).
!!
!!**  METHOD
!!    ------
!!      The method assumes a Csk law for the activation of aerososl with "s"
!!    the supersaturation (here 0.05 % is chosen). C and k are the Twomey
!!    coefficients initialized in INI_RAIN_C2R2. A Marshall-Palmer law with
!!    N_o=10**(-7) m**(-4) is assumed for the rain drop concentration.
!!      The initialization of the PSVT is straightforward for the cloud droplets
!!    while N_r=N_0/Lambda_r with Rho*r_r=Pi*Rho_w*N_0/(Lambda_r**4) is used for
!!    the rain drops. The HGETCLOUD test is used to discriminate between the
!!    'REVE' and 'KESS' options for CCLOUD in the previous run (from which
!!     PRT was calculated).
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_RAIN_C2R2_DESCR, ONLY : XRTMIN, XCTMIN
!!      Module MODD_RAIN_C2R2_KHKO_PARAM, ONLY : XCONCC_INI, XCONCR_PARAM_INI
!!      Module MODD_CONF,            ONLY : NVERB
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine SET_CONC_RAIN_C2R2 )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      P. Jabouille     * CNRM/GMME *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/11/00
!!                        2014 G.Delautier : remplace MODD_RAIN_C2R2_PARAM par MODD_RAIN_C2R2_KHKO_PARAM
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,                 ONLY : NVERB
USE MODD_LUNIT_n,              ONLY : TLUOUT
USE MODD_RAIN_C2R2_DESCR,      ONLY : XRTMIN, XCTMIN
USE MODD_RAIN_C2R2_KHKO_PARAM, ONLY : XCONCC_INI, XCONCR_PARAM_INI
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4),         INTENT(IN) :: HGETCLOUD  ! Get indicator
REAL, DIMENSION(:,:,:),    INTENT(IN) :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PRT     ! microphysical mixing ratios
!
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT):: PSVT     ! microphys. concentrations
!
!
!*       0.2   Declarations of local variables :
!
INTEGER    :: IRESP   ! Return code of FM routines
INTEGER    :: ILUOUT  ! Logical unit number of output-listing
REAL       :: ZRCMIN  ! Minimal value of PRT(,,,2)
REAL       :: ZRRMIN  ! Minimal value of PRT(,,,3)
REAL       :: ZCCMIN  ! Minimal value of PSVT(,,,2)
REAL       :: ZCRMIN  ! Minimal value of PSVT(,,,3)
REAL       :: ZCONCC, ZCONCR
!
!-------------------------------------------------------------------------------
!*       1.    RETRIEVE LOGICAL UNIT NUMBER
!              ----------------------------
!
ILUOUT = TLUOUT%NLU
!
!*       2.    INITIALIZATION
!              --------------
!
ZRCMIN = XRTMIN(2)
ZCCMIN = XCTMIN(2)
ZCONCC = XCONCC_INI
!
WHERE ( PRT(:,:,:,2) > ZRCMIN )
  PSVT(:,:,:,1) = ZCONCC
  PSVT(:,:,:,2) = ZCONCC
END WHERE
WHERE ( PRT(:,:,:,2) <= ZRCMIN )
  PRT(:,:,:,2)  = 0.0
  PSVT(:,:,:,1) = 0.0
  PSVT(:,:,:,2) = 0.0
END WHERE
IF( NVERB >= 5 ) THEN
  WRITE (UNIT=ILUOUT,FMT=*) "!INI_MODEL$n: The droplet concentration has "
  WRITE (UNIT=ILUOUT,FMT=*) "been roughly initialised with a Csk law"
END IF
!
IF (HGETCLOUD == 'INI1') THEN
  PSVT(:,:,:,3) = 0.0
ELSEIF (HGETCLOUD == 'INI2') THEN
  ZRRMIN = XRTMIN(3)
  ZCRMIN = XCTMIN(3)
  ZCONCR = XCONCR_PARAM_INI
!
  WHERE ( PRT(:,:,:,3) > ZRRMIN )
    PSVT(:,:,:,3) = MAX( SQRT(SQRT(PRHODREF(:,:,:)*PRT(:,:,:,3) &
                                     *XCONCR_PARAM_INI)),ZCRMIN )
  END WHERE
  WHERE ( PRT(:,:,:,3) <= ZRRMIN )
    PRT(:,:,:,3)  = 0.0
    PSVT(:,:,:,3) = 0.0
  END WHERE
  IF( NVERB >= 5 ) THEN
    WRITE (UNIT=ILUOUT,FMT=*) "!INI_MODEL$n: The raindrop concentration has "
    WRITE (UNIT=ILUOUT,FMT=*) "been roughly initialised with a MP law"
  END IF
END IF
!
END SUBROUTINE SET_CONC_RAIN_C2R2
