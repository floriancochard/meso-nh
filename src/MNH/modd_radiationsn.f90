!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ########################
       MODULE MODD_RADIATIONS_n
!      ########################
!
!!****  *MODD_RADIATIONS$n* - declaration of parameters and output fields of
!!                            the radiation scheme
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to define astronomical data
!!    and environmental variables for the radiation scheme. The module contains
!!    also the computed radiative fields (tendency and downward surface fluxes).
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_RADIATIONS$n)
!!
!!    AUTHOR
!!    ------
!!     J.-P. Pinty   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      26/02/95
!!      Revised       11/09/95
!!                   Mar.  15,2002 (F.Solmon) add ozone and aerosol fields
!!                    03/03/03 (V. Masson) surface radiative schemes and
!!                                         multiple wavelengths for surface SW
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE RADIATIONS_t
!
!
  INTEGER :: NDLON   ! number of columns where the radiation
		                 ! calculations are performed
  INTEGER :: NFLEV   ! number of vertical levels where the radiation
		                 ! calculations are performed
  INTEGER :: NFLUX   ! number of top and ground fluxes in the output
  INTEGER :: NRAD    ! number of satellite radiances to synthesize
  INTEGER :: NAER    ! number od AERosol classes
  INTEGER :: NSWB_OLD    ! number of SW bands in ECMWF original code (usually 6)
  INTEGER :: NSWB_MNH! number of SW bands practically used (14 if ECRAD, NSWB if original code) 
  INTEGER :: NLWB_MNH! number of LW bands practically used (16 if RRTM) 
  INTEGER :: NSTATM  ! index od the STAndard ATMosphere level just above
		                 ! the model top
!  INTEGER, DIMENSION(:,:), POINTER  :: NCLEARCOL_TM1=>NULL() ! trace of cloud/clear columns at
                                                             ! the previous radiation time step
!  (to be replaced by a logical array when fmread-writ could treat this kind of data)
  REAL    :: XCCO2   ! CO2 mixing ratio
  REAL    :: XSINDEL ! sine   of the solar declination angle
  REAL    :: XCOSDEL ! cosine of the solar declination angle
  REAL    :: XTSIDER ! sideral decimal time correction
  REAL    :: XCORSOL ! daily solar constant correction
  REAL    :: XTSRIS  ! time for sunrise
  REAL    :: XTSSET  ! time for sunset
!
  REAL, DIMENSION(:,:),   POINTER :: XSLOPANG=>NULL() ! slope angle
  REAL, DIMENSION(:,:),   POINTER :: XSLOPAZI=>NULL() ! azimuthal slope angle
  REAL, DIMENSION(:,:),   POINTER :: XSTATM=>NULL()   ! working standard atmosphere
!
  REAL, DIMENSION(:,:,:),   POINTER     :: XOZON=>NULL()  ! ozone mixing ratio 
  REAL, DIMENSION(:,:,:,:), POINTER     :: XAER=>NULL()   ! aerosol optical thickness
  REAL, DIMENSION(:,:,:,:), POINTER     :: XDST_WL=>NULL()   ! Extinction by wavelength
  REAL, DIMENSION(:,:,:,:), POINTER     :: XAER_CLIM=>NULL()   ! aerosol optical thickness
!
!  REAL, DIMENSION(:,:,:), POINTER :: XDTHRAD=>NULL()    ! radiative heating/cooling rate
  REAL, DIMENSION(:,:),   POINTER :: XSRFLWD=>NULL()    ! downward surface LW radiations
  REAL, DIMENSION(:,:),   POINTER :: XSCASRFSWD=>NULL() ! downward surface scaterred SW radiations
!
!  REAL, DIMENSION(:,:),   POINTER :: XFLALWD=>NULL()    ! downward LW radiations on flat surface
  REAL, DIMENSION(:,:),   POINTER :: XFLASWD=>NULL()    ! downward SW radiations on flat surface

  REAL, DIMENSION(:),     POINTER :: XSW_BANDS=>NULL()  ! value of the wavelentgh
!                                                       ! in the middle of each SW spectral band

  REAL, DIMENSION(:),     POINTER :: XLW_BANDS=>NULL()  ! value of the wavelentgh
!                                                       ! in the middle of each LW spectral band
!
!  REAL, DIMENSION(:,:,:), POINTER :: XDIRFLASWD=>NULL() ! downward surface direct    SW radiations
!                                                       ! BEFORE slope and shadows effects
!                                                       ! for each spectral band
!  REAL, DIMENSION(:,:,:), POINTER :: XSCAFLASWD=>NULL() ! downward surface diffuse   SW radiations
!                                                       ! BEFORE slope and shadows effects
!                                                       ! for each spectral band
!  REAL, DIMENSION(:,:,:), POINTER :: XDIRSRFSWD=>NULL() ! downward surface direct    SW radiations
!                                                       ! for each spectral band
!
!  REAL, DIMENSION(:,:,:),       POINTER :: XDIR_ALB=>NULL()   ! direct albedo for each spectral band
!  REAL, DIMENSION(:,:,:),       POINTER :: XSCA_ALB=>NULL()   ! scattered albedo for each spectral band
!  REAL, DIMENSION(:,:,:),       POINTER :: XEMIS=>NULL()      ! emissivity
!  REAL, DIMENSION(:,:),         POINTER :: XTSRAD=>NULL()      ! surface temperature
  REAL, DIMENSION(:,:),         POINTER :: XSEA=>NULL()        ! sea fraction
!
!  REAL, DIMENSION(:,:),         POINTER :: XZENITH=>NULL()    ! zenithal angle  (radian from the vertical)
!  REAL, DIMENSION(:,:),         POINTER :: XAZIM=>NULL()      ! azimuthal angle (radian from N, clockwise)
  REAL, DIMENSION(:,:),         POINTER :: XALBUV=>NULL()      ! UV albedo
  REAL, DIMENSION(:,:,:),       POINTER :: XSWU => NULL()      ! SW_UP
  REAL, DIMENSION(:,:,:),       POINTER :: XSWD => NULL()      ! SW_DOWN
  REAL, DIMENSION(:,:,:),       POINTER :: XLWU => NULL()      ! LW_UP
  REAL, DIMENSION(:,:,:),       POINTER :: XLWD => NULL()      ! LW_DOWN
  REAL, DIMENSION(:,:,:),       POINTER :: XDTHRADSW => NULL() ! DTHRAD SW
  REAL, DIMENSION(:,:,:),       POINTER :: XDTHRADLW => NULL() ! DTHRAD LW
  REAL, DIMENSION(:,:,:),       POINTER :: XRADEFF   => NULL() ! effective radius
!
END TYPE RADIATIONS_t

TYPE(RADIATIONS_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: RADIATIONS_MODEL

INTEGER, POINTER :: NDLON=>NULL()
INTEGER, POINTER :: NFLEV=>NULL()
INTEGER, POINTER :: NFLUX=>NULL()
INTEGER, POINTER :: NRAD=>NULL()
INTEGER, POINTER :: NAER=>NULL()
INTEGER, POINTER :: NSWB_OLD=>NULL()
INTEGER, POINTER :: NSWB_MNH=>NULL()
INTEGER, POINTER :: NLWB_MNH=>NULL()
INTEGER, POINTER :: NSTATM=>NULL()
INTEGER, DIMENSION(:,:), POINTER  :: NCLEARCOL_TM1=>NULL()
REAL, POINTER :: XCCO2=>NULL()
REAL, POINTER :: XSINDEL=>NULL()
REAL, POINTER :: XCOSDEL=>NULL()
REAL, POINTER :: XTSIDER=>NULL()
REAL, POINTER :: XCORSOL=>NULL()
REAL, POINTER :: XTSRIS=>NULL()
REAL, POINTER :: XTSSET=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XSLOPANG=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XSLOPAZI=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XSTATM=>NULL()
REAL, DIMENSION(:,:,:),   POINTER     :: XOZON=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER     :: XAER=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER     :: XDST_WL=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER     :: XAER_CLIM=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDTHRAD=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XSRFLWD=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XSCASRFSWD=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XFLALWD=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XFLASWD=>NULL()
REAL, DIMENSION(:),     POINTER :: XSW_BANDS=>NULL()
REAL, DIMENSION(:),     POINTER :: XLW_BANDS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDIRFLASWD=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSCAFLASWD=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDIRSRFSWD=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDIR_ALB=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSCA_ALB=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XEMIS=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XTSRAD=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XSEA=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XZENITH=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XAZIM=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XALBUV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSWU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSWD=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLWU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLWD=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDTHRADSW=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XDTHRADLW=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRADEFF=>NULL()

CONTAINS

SUBROUTINE RADIATIONS_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!RADIATIONS_MODEL(KFROM)%NCLEARCOL_TM1=>NCLEARCOL_TM1 !Done in FIELDLIST_GOTO_MODEL
RADIATIONS_MODEL(KFROM)%XSLOPANG=>XSLOPANG
RADIATIONS_MODEL(KFROM)%XSLOPAZI=>XSLOPAZI
RADIATIONS_MODEL(KFROM)%XSTATM=>XSTATM
RADIATIONS_MODEL(KFROM)%XOZON=>XOZON
RADIATIONS_MODEL(KFROM)%XAER=>XAER
RADIATIONS_MODEL(KFROM)%XDST_WL=>XDST_WL
RADIATIONS_MODEL(KFROM)%XAER_CLIM=>XAER_CLIM
!RADIATIONS_MODEL(KFROM)%XDTHRAD=>XDTHRAD !Done in FIELDLIST_GOTO_MODEL
RADIATIONS_MODEL(KFROM)%XSRFLWD=>XSRFLWD
RADIATIONS_MODEL(KFROM)%XSCASRFSWD=>XSCASRFSWD
!RADIATIONS_MODEL(KFROM)%XFLALWD=>XFLALWD !Done in FIELDLIST_GOTO_MODEL
RADIATIONS_MODEL(KFROM)%XFLASWD=>XFLASWD
RADIATIONS_MODEL(KFROM)%XSW_BANDS=>XSW_BANDS
RADIATIONS_MODEL(KFROM)%XLW_BANDS=>XLW_BANDS
!RADIATIONS_MODEL(KFROM)%XDIRFLASWD=>XDIRFLASWD !Done in FIELDLIST_GOTO_MODEL
!RADIATIONS_MODEL(KFROM)%XSCAFLASWD=>XSCAFLASWD !Done in FIELDLIST_GOTO_MODEL
!RADIATIONS_MODEL(KFROM)%XDIRSRFSWD=>XDIRSRFSWD !Done in FIELDLIST_GOTO_MODEL
!RADIATIONS_MODEL(KFROM)%XDIR_ALB=>XDIR_ALB !Done in FIELDLIST_GOTO_MODEL
!RADIATIONS_MODEL(KFROM)%XSCA_ALB=>XSCA_ALB !Done in FIELDLIST_GOTO_MODEL
!RADIATIONS_MODEL(KFROM)%XEMIS=>XEMIS !Done in FIELDLIST_GOTO_MODEL
!RADIATIONS_MODEL(KFROM)%XTSRAD=>XTSRAD !Done in FIELDLIST_GOTO_MODEL
RADIATIONS_MODEL(KFROM)%XSEA=>XSEA
!RADIATIONS_MODEL(KFROM)%XZENITH=>XZENITH !Done in FIELDLIST_GOTO_MODEL
!RADIATIONS_MODEL(KFROM)%XAZIM=>XAZIM !Done in FIELDLIST_GOTO_MODEL
RADIATIONS_MODEL(KFROM)%XALBUV=>XALBUV
RADIATIONS_MODEL(KFROM)%XSWU=>XSWU
RADIATIONS_MODEL(KFROM)%XSWD=>XSWD
RADIATIONS_MODEL(KFROM)%XLWU=>XLWU
RADIATIONS_MODEL(KFROM)%XLWD=>XLWD
RADIATIONS_MODEL(KFROM)%XDTHRADSW=>XDTHRADSW
RADIATIONS_MODEL(KFROM)%XDTHRADLW=>XDTHRADLW
RADIATIONS_MODEL(KFROM)%XRADEFF=>XRADEFF
!
! Current model is set to model KTO
NDLON=>RADIATIONS_MODEL(KTO)%NDLON
NFLEV=>RADIATIONS_MODEL(KTO)%NFLEV
NFLUX=>RADIATIONS_MODEL(KTO)%NFLUX
NRAD=>RADIATIONS_MODEL(KTO)%NRAD
NAER=>RADIATIONS_MODEL(KTO)%NAER
NSWB_OLD=>RADIATIONS_MODEL(KTO)%NSWB_OLD
NSWB_MNH=>RADIATIONS_MODEL(KTO)%NSWB_MNH
NLWB_MNH=>RADIATIONS_MODEL(KTO)%NLWB_MNH
NSTATM=>RADIATIONS_MODEL(KTO)%NSTATM
!NCLEARCOL_TM1=>RADIATIONS_MODEL(KTO)%NCLEARCOL_TM1 !Done in FIELDLIST_GOTO_MODEL
XCCO2=>RADIATIONS_MODEL(KTO)%XCCO2
XSINDEL=>RADIATIONS_MODEL(KTO)%XSINDEL
XCOSDEL=>RADIATIONS_MODEL(KTO)%XCOSDEL
XTSIDER=>RADIATIONS_MODEL(KTO)%XTSIDER
XCORSOL=>RADIATIONS_MODEL(KTO)%XCORSOL
XTSRIS=>RADIATIONS_MODEL(KTO)%XTSRIS
XTSSET=>RADIATIONS_MODEL(KTO)%XTSSET
XSLOPANG=>RADIATIONS_MODEL(KTO)%XSLOPANG
XSLOPAZI=>RADIATIONS_MODEL(KTO)%XSLOPAZI
XSTATM=>RADIATIONS_MODEL(KTO)%XSTATM
XOZON=>RADIATIONS_MODEL(KTO)%XOZON
XAER=>RADIATIONS_MODEL(KTO)%XAER
XDST_WL=>RADIATIONS_MODEL(KTO)%XDST_WL
XAER_CLIM=>RADIATIONS_MODEL(KTO)%XAER_CLIM
!XDTHRAD=>RADIATIONS_MODEL(KTO)%XDTHRAD !Done in FIELDLIST_GOTO_MODEL
XSRFLWD=>RADIATIONS_MODEL(KTO)%XSRFLWD
XSCASRFSWD=>RADIATIONS_MODEL(KTO)%XSCASRFSWD
!XFLALWD=>RADIATIONS_MODEL(KTO)%XFLALWD !Done in FIELDLIST_GOTO_MODEL
XFLASWD=>RADIATIONS_MODEL(KTO)%XFLASWD
XSW_BANDS=>RADIATIONS_MODEL(KTO)%XSW_BANDS
XLW_BANDS=>RADIATIONS_MODEL(KTO)%XLW_BANDS
!XDIRFLASWD=>RADIATIONS_MODEL(KTO)%XDIRFLASWD !Done in FIELDLIST_GOTO_MODEL
!XSCAFLASWD=>RADIATIONS_MODEL(KTO)%XSCAFLASWD !Done in FIELDLIST_GOTO_MODEL
!XDIRSRFSWD=>RADIATIONS_MODEL(KTO)%XDIRSRFSWD !Done in FIELDLIST_GOTO_MODEL
!XDIR_ALB=>RADIATIONS_MODEL(KTO)%XDIR_ALB !Done in FIELDLIST_GOTO_MODEL
!XSCA_ALB=>RADIATIONS_MODEL(KTO)%XSCA_ALB !Done in FIELDLIST_GOTO_MODEL
!XEMIS=>RADIATIONS_MODEL(KTO)%XEMIS !Done in FIELDLIST_GOTO_MODEL
!XTSRAD=>RADIATIONS_MODEL(KTO)%XTSRAD !Done in FIELDLIST_GOTO_MODEL
XSEA=>RADIATIONS_MODEL(KTO)%XSEA
!XZENITH=>RADIATIONS_MODEL(KTO)%XZENITH !Done in FIELDLIST_GOTO_MODEL
!XAZIM=>RADIATIONS_MODEL(KTO)%XAZIM !Done in FIELDLIST_GOTO_MODEL
XALBUV=>RADIATIONS_MODEL(KTO)%XALBUV
XSWU=>RADIATIONS_MODEL(KTO)%XSWU
XSWD=>RADIATIONS_MODEL(KTO)%XSWD
XLWU=>RADIATIONS_MODEL(KTO)%XLWU
XLWD=>RADIATIONS_MODEL(KTO)%XLWD
XDTHRADSW=>RADIATIONS_MODEL(KTO)%XDTHRADSW
XDTHRADLW=>RADIATIONS_MODEL(KTO)%XDTHRADLW
XRADEFF=>RADIATIONS_MODEL(KTO)%XRADEFF

END SUBROUTINE RADIATIONS_GOTO_MODEL

END MODULE MODD_RADIATIONS_n
