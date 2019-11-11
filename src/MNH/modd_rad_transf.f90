!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################
      MODULE MODD_RAD_TRANSF
!     ######################
!
!!****  *MODD_RAD_TRANSF* - declaration of variables related to the radiatif 
!!                          transfer diagnostic code
!!
!!    PURPOSE
!!    -------
!
!!
!!**  METHOD
!!    ------
!!  
!!      from main COMMON block inclusions
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      J.-P. Chaboureau       *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       29/03/00     
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: JPGEOST= 5 !number of different geostationnary satellites 
INTEGER, PARAMETER :: JPSAT=7 !number of geostationnary satellites for Meteosat
!
INTEGER, PARAMETER :: JPNINT=225, JPV2=2 , JPV3=3 , JPV10=10 ! 'param.h'
INTEGER, PARAMETER :: JPUABS=14                              !
INTEGER, PARAMETER :: JPWVINT=20  ! number of spectral interval for WV band
INTEGER, PARAMETER :: JPCAN=2     ! number of channels
!
INTEGER, SAVE    :: NMP,NIMP,NULOUT,NULINA,NULNAM           ! COMMON YOMIO
INTEGER, SAVE    :: NABS,NATM,NATMS,NSPWV,NTMP1,NTEMP,    & ! COMMON YOMRADI
              NCH2O,NCCO2,NCO3,NCHAN,NH2O,NCO2,NO3, &
              NCNT,NN2O,NCH4,NCO,NC11,NC12,NCFC,NO2
!
REAL, SAVE  :: X1CO2,XN2O,XCO,XCH4,XF11,XF12,XO2            ! COMMON YOMRADR
REAL,DIMENSION(JPV3), SAVE :: XTEMP,XH2O,X3CO2,XO3          !
REAL,DIMENSION(JPV10), SAVE:: XXLIM                         !
REAL,DIMENSION(JPV2), SAVE :: XCLIM                         !
!
INTEGER, SAVE    :: N_INT,NA                                ! COMMON YOMSPEI
REAL, SAVE       :: XTREF,XTPOLY                            ! 
!
REAL,DIMENSION(JPV2), SAVE   :: XWNU                        ! COMMON YOMSPER
REAL,DIMENSION(6,2), SAVE    :: XPOLPLCK                    !
REAL,DIMENSION(6,8), SAVE    :: XRODWAL                     !
REAL,DIMENSION(JPNINT), SAVE  :: XWVNA                      !
REAL,DIMENSION(JPNINT), SAVE  :: XWVNB                      !
REAL, SAVE :: XALPHA                                        !
REAL, SAVE :: XAIRM,XH2OM,XCO2M,XO3M,XN2OM,XCOM,XCH4M,XO2M  !
REAL, SAVE :: XF11M,XF12M                                   !
!
REAL,DIMENSION(2,20), SAVE   :: XWNUTOT                     ! COMMON YOMSPET
REAL,DIMENSION(6,2,20), SAVE :: XPOLPLCKTOT                 !
REAL,DIMENSION(6,8,20), SAVE :: XRODWALTOT                  !
!
REAL,DIMENSION(2), SAVE      :: XRT1,XWG1                   ! COMMON YOMGOS
!
END MODULE MODD_RAD_TRANSF
