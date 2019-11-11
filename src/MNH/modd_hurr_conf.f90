!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     #####################
      MODULE MODD_HURR_CONF
!     #####################
!
!!****  *MODD_HURR_PARAM* - declaration of the differents parameters
!                           used in the filtering of the large-scale fields
!                           and in the bogussing (inclusion of an
!                           analytical vortex).
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     the parameters for the tropical cyclone simulations (filtering +
!                      bogussing)      
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	O. Nuissier   * L.A *
!!      R. Rogers     * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/12/01                      
!!      Modification   09/08/05 (D.Barbary) Add CDADBOGFILE and CDADATMFILE
!!                     01/02/08 (D.Barbary) Add Holland's parameters
!!                                          XC, XMAX, XRHO_Z, XRHO_ZZ, XB_0, XBETA_Z, XBETA_ZZ
!!                                          and convergence angle
!!                                          XANGCONV0, XANGCONV1000, XANGCONV2000
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
LOGICAL,SAVE       ::   LFILTERING  ! Logical to switch on the filtering of
                                    ! the large-scale vortex
CHARACTER(LEN=5),SAVE  ::   CFILTERING  ! List of fields to be filtered
INTEGER,SAVE       ::   NK          ! Parameter for Barnes filter
REAL   ,SAVE       ::   XLAMBDA     !    "
REAL   ,SAVE       ::   XLATGUESS   ! Initial guess position of the cyclone
REAL   ,SAVE       ::   XLONGUESS   ! for the center determination (via the
                                    ! "simplex-like" or dynamical wind method)
REAL   ,SAVE       ::   XBOXWIND    ! Half-length around the fix position
                                    ! for the dynamical center
REAL   ,SAVE       ::   XRADGUESS   ! First guess of the radius vortex
                                    !of the filter domain (km)
INTEGER,SAVE       ::   NLEVELR0    ! sigma level used to compute R0
INTEGER,SAVE       ::   NPHIL       ! Number of azimuthal directions
INTEGER,SAVE       ::   NDIAG_FILT  ! to write the different components 
                                    !computated in the LFI-file
                                    !
LOGICAL,SAVE       ::   LBOGUSSING  ! Logical to switch on the addition of
                                    ! the bogus vortex in the filtered fields
REAL   ,SAVE       ::   XLATBOG     ! Position of the bogus vortex
REAL   ,SAVE       ::   XLONBOG     !
REAL   ,SAVE       ::   XVTMAXSURF  ! Maximum tangential wind near the
                                    ! surface or about 500 m altitude
REAL   ,SAVE       ::   XRADWINDSURF! Radius of maximum wind near the
                                    ! surface or about 500 m altitude
REAL   ,SAVE       ::   XC          ! standard coefficient for VTmax(z)
REAL   ,SAVE       ::   XRHO_Z, XRHO_ZZ ! standard coefficient for RMW(z)
REAL   ,SAVE       ::   XB_0,XBETA_Z, XBETA_ZZ ! standard coefficient for B(z)
REAL   ,SAVE       ::   XMAX        ! altitude where the tangentiel wind vanishes
REAL   ,SAVE       ::   XANGCONV0   ! Convergence angle near the surface
REAL   ,SAVE       ::   XANGCONV1000! Convergence angle at 1000m altitude
REAL   ,SAVE       ::   XANGCONV2000! Convergence angle at 2000m altitude
CHARACTER(LEN=28),SAVE :: CDADATMFILE  ! Name of the dad of HATMFILE 
CHARACTER(LEN=28),SAVE :: CDADBOGFILE  ! Name of the dad of CINIFILE
 END MODULE MODD_HURR_CONF
