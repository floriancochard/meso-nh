!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 rad 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SUNPOS_n
!     ####################
!
INTERFACE
!
      SUBROUTINE SUNPOS_n (PZENITH, PCOSZEN, PSINZEN, PAZIMSOL)
!
REAL, DIMENSION(:,:),           INTENT(OUT)  :: PZENITH    ! Solar zenithal angle
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)  :: PCOSZEN    ! Cosine of the solar zenithal angle
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)  :: PSINZEN    ! Sinus  of the solar zenithal angle
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)  :: PAZIMSOL   ! Solar azimuthal angle
!
END SUBROUTINE SUNPOS_n
!
END INTERFACE
!
END MODULE MODI_SUNPOS_n
!
!
!     #########################################################
      SUBROUTINE SUNPOS_n (PZENITH, PCOSZEN, PSINZEN, PAZIMSOL)
!     #########################################################
!
!!****  *SUNPOS_n * - routine to compute the position of the sun
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cosine and sinus of the 
!!    solar zenithal angle (angle defined by the local vertical at the position
!!    XLAT, XLON and the direction of the sun) and the azimuthal solar
!!    angle (angle between an horizontal direction (south or north according
!!    to the terrestrial hemisphere) and the horizontal projection of the
!!    direction of the sun.
!!
!!**  METHOD
!!    ------
!!      The cosine and sinus of the zenithal solar angle  and the azimuthal 
!!    solar angle are computed from the true universal time, valid for the (XLAT,
!!    XLON) location, and from the solar declination angle of the day. There
!!    is a special convention to define the azimuthal solar angle.
!!     
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      "Radiative Processes in Meteorology and Climatology"  
!!                          (1976)   Paltridge and Platt 
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             16/10/94 
!!      Revised              12/09/95
!!      (J.Stein)            01:04/96  bug correction for ZZEANG     
!!      (K. Suhre)           14/02/97  bug correction for ZLON0     
!!      (V. Masson)          01/03/03  add zenithal angle output
!!      (V. Masson)          04/01/12  standard definition of Azimuthal angle
!!                                     (from North, clockwise)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF,         ONLY : LCARTESIAN
USE MODD_CST,          ONLY : XPI
USE MODD_GRID,         ONLY : XLAT0, XLON0
USE MODD_GRID_n,       ONLY : XLAT,  XLON
USE MODD_PARAM_RAD_n,  ONLY : XDTRAD
USE MODD_RADIATIONS_n, ONLY : XSINDEL, XCOSDEL, XTSIDER
USE MODD_TIME_n,       ONLY : TDTRAD_FULL
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL, DIMENSION(:,:),           INTENT(OUT)  :: PZENITH    ! Solar zenithal angle
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)  :: PCOSZEN    ! Cosine of the solar zenithal angle
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)  :: PSINZEN    ! Sinus  of the solar zenithal angle
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)  :: PAZIMSOL   ! Solar azimuthal angle
!                                                          ! (radian, from North, clockwise)
!
!*       0.2   declarations of local variables
!
!
REAL :: ZTIME   ! Centered current time for radiation calculations
REAL :: ZUT     ! Universal Time
!
REAL, DIMENSION(SIZE(PZENITH,1),SIZE(PZENITH,2)) :: ZTUT    ,&! True (absolute)
                                                              ! Universal Time
                           ZSOLANG ,&! Hourly solar angle
                           ZSINAZI ,&! Sine of the solar azimuthal angle
                           ZCOSAZI ,&! Cosine of the solar azimuthal angle
                           ZLAT,    &
                           ZLON,    &! Array of latitudes and longitudes
                           ZSINZEN, &!Sine of zenithal angle
                           ZCOSZEN, &!Cosine of zenithal angle
                           ZAZIMSOL  !azimuthal angle
!
!-------------------------------------------------------------------------------
!
!*       1.    LOADS THE ZLAT, ZLON ARRAYS
!              ---------------------------
!
IF(LCARTESIAN) THEN
  ZLAT(:,:) = XLAT0*(XPI/180.)
  ZLON(:,:) = XLON0*(XPI/180.)
  ELSE
  ZLAT(:,:) = XLAT(:,:)*(XPI/180.)
  ZLON(:,:) = XLON(:,:)*(XPI/180.)
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTES THE TRUE SOLAR TIME
!              ----------------------------
!
ZTIME     = TDTRAD_FULL%TIME + 0.5*XDTRAD


ZUT       = MOD( 24.0+MOD(ZTIME/3600.,24.0),24.0 )
!
ZTUT(:,:) = ZUT - XTSIDER + ZLON(:,:)*((180./XPI)/15.0)
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTES THE COSINE AND SINUS OF THE ZENITHAL SOLAR ANGLE
!              ---------------------------------------------------------
!
ZSOLANG(:,:) = (ZTUT(:,:)-12.0)*15.0*(XPI/180.)          ! hour angle in radians
!
ZCOSZEN(:,:) = SIN(ZLAT(:,:))*XSINDEL +                 &! Cosine of the zenithal
               COS(ZLAT(:,:))*XCOSDEL*COS(ZSOLANG(:,:))  !       solar angle
!
ZSINZEN(:,:)  = SQRT( 1. - ZCOSZEN(:,:)*ZCOSZEN(:,:) )
!
!-------------------------------------------------------------------------------
!
!*       4.    ZENITHAL SOLAR ANGLE
!              --------------------
!
PZENITH(:,:) = ACOS(ZCOSZEN(:,:))
!
!-------------------------------------------------------------------------------
!
!*       5.    COMPUTE THE AZIMUTHAL SOLAR ANGLE (PAZIMSOL)
!              --------------------------------------------
!
WHERE (ZSINZEN(:,:)/=0.)
  ZSINAZI(:,:)  = - XCOSDEL * SIN(ZSOLANG(:,:)) / ZSINZEN(:,:)
  ZCOSAZI(:,:)  = (-SIN(ZLAT(:,:))*XCOSDEL*COS(ZSOLANG(:,:))     & 
                   +COS(ZLAT(:,:))*XSINDEL                       &
                  ) / ZSINZEN(:,:)
  ZAZIMSOL(:,:) = ATAN2(ZSINAZI(:,:),ZCOSAZI(:,:))
ELSEWHERE
  ZAZIMSOL(:,:) = XPI
END WHERE
!
!-------------------------------------------------------------------------------
IF (PRESENT(PCOSZEN )) PCOSZEN =ZCOSZEN
IF (PRESENT(PSINZEN )) PSINZEN =ZSINZEN
IF (PRESENT(PAZIMSOL)) PAZIMSOL=ZAZIMSOL
!-------------------------------------------------------------------------------
!
END SUBROUTINE SUNPOS_n

