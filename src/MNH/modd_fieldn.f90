!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_FIELD_n
!     ###################
!
!!****  *MODD_FIELD$n* - declaration of prognostic variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     prognostic variables. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_FIELDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       05/05/94                      
!!      Modifications  03/01/95  (Lafore)  To add the dry mass variables Md  
!!                     09/03/95  (Stein)   eliminate R from the progn. var                    
!!                     15/03/95  (Stein)   add EPS variable
!!      Modifications  21/03/95  (Carriere) To add the subgrid condensation 
!!                                           related parameters
!!                     01/03/96  (J. Stein) add the cloud fraction
!!                     10/10/96  (J. Stein) add XSRCM and XSRCT
!!                     11/04/96  (J.-P. Pinty) add the ice concentration
!!                     25/07/97  (J. Stein) Change the variable pressure
!!                     20/05/06  Remove EPS
!!                        11/11  (C.Lac) FIT version : Remove t-Dt fields except for 
!!                                       radiative cooling (microphysics) +
!!                               add pressure contribution to the tendencies for
!!                               momentum (noted _PRES) + microphysics contrib
!!                               for Theta and r (noted _CLD)
!!                     04/16    (M.Mazoyer) New supersaturation fields
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                   02/2019 C.Lac add rain fraction as an output field
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 06/03/2019: correct XZWS entry
!  P. Wautelet 14/03/2019: add XZWS_DEFAULT parameter
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

REAL, PARAMETER :: XZWS_DEFAULT = 2. ! Default value for XZWS: 2 meters

TYPE FIELD_t
  REAL, DIMENSION(:,:),   POINTER :: XZWS=>NULL()    ! significant sea wave
!  REAL, DIMENSION(:,:,:), POINTER :: XUT=>NULL(),XVT=>NULL(),XWT=>NULL()
                                      ! U,V,W  at time t
  REAL, DIMENSION(:,:,:), POINTER :: XRUS=>NULL(),XRVS=>NULL(),XRWS=>NULL()
                                      ! Source of (rho U), (rho V), (rho w) 
!  REAL, DIMENSION(:,:,:), POINTER :: XRUS_PRES=>NULL(),XRVS_PRES=>NULL(),XRWS_PRES=>NULL()
!  REAL, DIMENSION(:,:,:), POINTER :: XTHT=>NULL()     ! (rho theta) at time t
  REAL, DIMENSION(:,:,:), POINTER :: XRTHS=>NULL()    ! Source of (rho theta)
!  REAL, DIMENSION(:,:,:), POINTER :: XRTHS_CLD=>NULL()    ! Source of (rho theta) from resolved_cloud
!  REAL, DIMENSION(:,:,:), POINTER :: XSUPSAT=>NULL()    ! Sursat
!                                                     ! at time t
!  REAL, DIMENSION(:,:,:), POINTER :: XNACT=>NULL()    ! Sursat
!                                                     ! at time t
!  REAL, DIMENSION(:,:,:), POINTER :: XNPRO=>NULL()    ! Sursat
!                                                     ! at time t
!  REAL, DIMENSION(:,:,:), POINTER :: XSSPRO=>NULL()    ! Sursat
!                                                     ! at time t
!  REAL, DIMENSION(:,:,:), POINTER :: XTKET=>NULL()    ! Kinetic energy
!                                                     ! at time t
  REAL, DIMENSION(:,:,:), POINTER :: XRTKES=>NULL()   ! Source of kinetic energy
                                                     ! (rho e)
!  REAL, DIMENSION(:,:,:), POINTER :: XPABST=>NULL()   ! absolute pressure at
!                                                     ! time t
!  REAL, DIMENSION(:,:,:,:), POINTER :: XRT=>NULL()    ! Moist variables (rho Rn) 
!                                                     ! at time t
  REAL, DIMENSION(:,:,:,:), POINTER :: XRRS=>NULL()   ! Source of Moist variables
                                                     ! (rho Rn) 
!  REAL, DIMENSION(:,:,:,:), POINTER :: XRRS_CLD=>NULL()   ! Source of Moist variables
  REAL, DIMENSION(:,:,:,:), POINTER :: XSVT=>NULL()   ! Additionnal scalar
                                                     ! variables at time t  
  REAL, DIMENSION(:,:,:,:), POINTER :: XRSVS=>NULL()  ! Source of addi. scalar
                                                     !  variables (rho Sn.) 
  REAL, DIMENSION(:,:,:,:), POINTER :: XRSVS_CLD=>NULL() ! Source of (rho Sn) from resolved_cloud
!  REAL                          ::   XDRYMASST    ! Mass of dry air Md
  REAL                          ::   XDRYMASSS    ! LS sources of Md
  REAL, DIMENSION(:,:,:), POINTER :: XSRC=>NULL()     ! turbulent flux <s'Rc'>
!  REAL, DIMENSION(:,:,:), POINTER :: XSRCT=>NULL()    ! turbulent flux <s'Rc'>
!                                                     ! at t
!  REAL, DIMENSION(:,:,:), POINTER :: XSIGS=>NULL()    ! =sqrt(<s's'>) for the
!                                                     ! Subgrid Condensation
!  REAL, DIMENSION(:,:,:), POINTER :: XCLDFR=>NULL()   ! cloud fraction
!  REAL, DIMENSION(:,:,:), POINTER :: XCIT=>NULL()     ! Pristine ice concentration
  REAL, DIMENSION(:,:,:), POINTER :: XTHM=>NULL()      ! Theta at Previous time step 
  REAL, DIMENSION(:,:,:), POINTER :: XRCM=>NULL()      ! Cloud mixing ratio at Previous time step
  REAL, DIMENSION(:,:,:), POINTER :: XPABSM=>NULL()      ! Theta at Previous time step 
!
END TYPE FIELD_t

TYPE(FIELD_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: FIELD_MODEL

REAL, DIMENSION(:,:),   POINTER :: XZWS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XUT=>NULL(),XVT=>NULL(),XWT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRUS=>NULL(),XRVS=>NULL(),XRWS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRUS_PRES=>NULL(),XRVS_PRES=>NULL(),XRWS_PRES=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTHT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRTHS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRTHS_CLD=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSUPSAT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XNACT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XNPRO=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSSPRO=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTKET=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRTKES=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XPABST=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XRT=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XRRS=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XRRS_CLD=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XSVT=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XRSVS=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XRSVS_CLD=>NULL()
REAL, POINTER :: XDRYMASST=>NULL()
REAL, POINTER :: XDRYMASSS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSRC=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSRCT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XSIGS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCLDFR=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRAINFR=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCIT=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTHM=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XPABSM=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRCM=>NULL()

CONTAINS

SUBROUTINE FIELD_GOTO_MODEL(KFROM, KTO)
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
INTEGER :: IID,IRESP
!
! Save current state for allocated arrays
!FIELD_MODEL(KFROM)%XZWS=>XZWS !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XUT=>XUT !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XVT=>XVT !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XWT=>XWT !Done in FIELDLIST_GOTO_MODEL
FIELD_MODEL(KFROM)%XRUS=>XRUS
FIELD_MODEL(KFROM)%XRVS=>XRVS
FIELD_MODEL(KFROM)%XRWS=>XRWS
!FIELD_MODEL(KFROM)%XRUS_PRES=>XRUS_PRES !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XRVS_PRES=>XRVS_PRES !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XRWS_PRES=>XRWS_PRES !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XTHT=>XTHT !Done in FIELDLIST_GOTO_MODEL
FIELD_MODEL(KFROM)%XRTHS=>XRTHS
!FIELD_MODEL(KFROM)%XRTHS_CLD=>XRTHS_CLD !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XSUPSAT=>XSUPSAT !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XNACT=>XNACT !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XNPRO=>XNPRO !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XSSPRO=>XSSPRO !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XTKET=>XTKET !Done in FIELDLIST_GOTO_MODEL
FIELD_MODEL(KFROM)%XRTKES=>XRTKES
!FIELD_MODEL(KFROM)%XPABST=>XPABST !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XRT=>XRT !Done in FIELDLIST_GOTO_MODEL
FIELD_MODEL(KFROM)%XRRS=>XRRS
!FIELD_MODEL(KFROM)%XRRS_CLD=>XRRS_CLD !Done in FIELDLIST_GOTO_MODEL
FIELD_MODEL(KFROM)%XSVT=>XSVT
FIELD_MODEL(KFROM)%XRSVS=>XRSVS
FIELD_MODEL(KFROM)%XRSVS_CLD=>XRSVS_CLD
FIELD_MODEL(KFROM)%XSRC=>XSRC
!FIELD_MODEL(KFROM)%XSRCT=>XSRCT !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XSIGS=>XSIGS !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XCLDFR=>XCLDFR !Done in FIELDLIST_GOTO_MODEL
!FIELD_MODEL(KFROM)%XCIT=>XCIT !Done in FIELDLIST_GOTO_MODEL
FIELD_MODEL(KFROM)%XTHM=>XTHM
FIELD_MODEL(KFROM)%XPABSM=>XPABSM
FIELD_MODEL(KFROM)%XRCM=>XRCM
!
! Current model is set to model KTO
!XZWS=>FIELD_MODEL(KTO)%XZWS !Done in FIELDLIST_GOTO_MODEL
!XUT=>FIELD_MODEL(KTO)%XUT !Done in FIELDLIST_GOTO_MODEL
!XVT=>FIELD_MODEL(KTO)%XVT !Done in FIELDLIST_GOTO_MODEL
!XWT=>FIELD_MODEL(KTO)%XWT !Done in FIELDLIST_GOTO_MODEL
XRUS=>FIELD_MODEL(KTO)%XRUS
XRVS=>FIELD_MODEL(KTO)%XRVS
XRWS=>FIELD_MODEL(KTO)%XRWS
!XRUS_PRES=>FIELD_MODEL(KTO)%XRUS_PRES !Done in FIELDLIST_GOTO_MODEL
!XRVS_PRES=>FIELD_MODEL(KTO)%XRVS_PRES !Done in FIELDLIST_GOTO_MODEL
!XRWS_PRES=>FIELD_MODEL(KTO)%XRWS_PRES !Done in FIELDLIST_GOTO_MODEL
!XTHT=>FIELD_MODEL(KTO)%XTHT !Done in FIELDLIST_GOTO_MODEL
XRTHS=>FIELD_MODEL(KTO)%XRTHS
!XRTHS_CLD=>FIELD_MODEL(KTO)%XRTHS_CLD !Done in FIELDLIST_GOTO_MODEL
!XSUPSAT=>FIELD_MODEL(KTO)%XSUPSAT !Done in FIELDLIST_GOTO_MODEL
!XNACT=>FIELD_MODEL(KTO)%XNACT !Done in FIELDLIST_GOTO_MODEL
!XSSPRO=>FIELD_MODEL(KTO)%XSSPRO !Done in FIELDLIST_GOTO_MODEL
!XNPRO=>FIELD_MODEL(KTO)%XNPRO !Done in FIELDLIST_GOTO_MODEL
!XTKET=>FIELD_MODEL(KTO)%XTKET !Done in FIELDLIST_GOTO_MODEL
XRTKES=>FIELD_MODEL(KTO)%XRTKES
!XPABST=>FIELD_MODEL(KTO)%XPABST !Done in FIELDLIST_GOTO_MODEL
!XRT=>FIELD_MODEL(KTO)%XRT !Done in FIELDLIST_GOTO_MODEL
XRRS=>FIELD_MODEL(KTO)%XRRS
!XRRS_CLD=>FIELD_MODEL(KTO)%XRRS_CLD !Done in FIELDLIST_GOTO_MODEL
XSVT=>FIELD_MODEL(KTO)%XSVT
XRSVS=>FIELD_MODEL(KTO)%XRSVS
XRSVS_CLD=>FIELD_MODEL(KTO)%XRSVS_CLD
!XDRYMASST=>FIELD_MODEL(KTO)%XDRYMASST !Done in FIELDLIST_GOTO_MODEL
XDRYMASSS=>FIELD_MODEL(KTO)%XDRYMASSS
XSRC=>FIELD_MODEL(KTO)%XSRC
!XSRCT=>FIELD_MODEL(KTO)%XSRCT !Done in FIELDLIST_GOTO_MODEL
!XSIGS=>FIELD_MODEL(KTO)%XSIGS !Done in FIELDLIST_GOTO_MODEL
!XCLDFR=>FIELD_MODEL(KTO)%XCLDFR !Done in FIELDLIST_GOTO_MODEL
!XCIT=>FIELD_MODEL(KTO)%XCIT !Done in FIELDLIST_GOTO_MODEL
XTHM=>FIELD_MODEL(KTO)%XTHM
XPABSM=>FIELD_MODEL(KTO)%XPABSM
XRCM=>FIELD_MODEL(KTO)%XRCM

END SUBROUTINE FIELD_GOTO_MODEL

END MODULE MODD_FIELD_n
