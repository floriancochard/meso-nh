!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! NEC0 masdev4_7 2007/06/16 01:41:59
!-----------------------------------------------------------------
!     #################
      MODULE MODD_GET_n
!     #################
!
!!****  *MODD_GET$n* - declaration of variables about getting of variables in
!!                    initialization
!!                     
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the variables 
!     which indicate how and what variables to get in initialization.
!        
!          
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_GETn)     
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        05/05/94  
!!      modification    22/11/94 (J.Stein) add the get indicators for PHIT,PHIM
!!      modification    07/12/94 (J.Stein) add LS fields get indicators  
!!      modification    15/03/95 (J.Stein) remove R from the historical var. 
!!      modification    15/06/95 (J.Stein) add EPS related variables
!!      modification    15/04/96 (J.Stein) add indicator for the cloud fraction
!!      modification    15/04/96 (J.Stein) add indicator for SRCM and T
!!      modification    11/04/96 (J.-P. Pinty) add the CLDFR and ice conc.
!!      modification    25/07/97 (J.Stein) change the pressure variable
!!      modification    25/07/99 (J.Stein) add get indicators for soil, rad and
!!                                         conv schemes
!!      J.-P. Pinty 25/10/00  add get indicator for cloud scheme
!!      V. Masson   01/2004   surface externalization (rm CGETSURF)
!!                  05/2006   Remove EPS and LGETALL
!!      M. Leriche  04/2010   add get indicators for pH in cloud and rain
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX, JPSVMAX
IMPLICIT NONE

TYPE GET_t
!
  CHARACTER (LEN=4) :: CGETUT, CGETVT, CGETWT !  Get indicator for
                                                 ! U,V,W at time t 
  CHARACTER (LEN=4) :: CGETTHT                !  Get indicator for theta
                                                 ! at time t-dt and at time t
  CHARACTER (LEN=4) :: CGETPABST              !  Get indicator for
                                                 ! the absolute pressure at
                                                 ! time t-dt and t
  CHARACTER (LEN=4)  :: CGETTKET              !  Get indicator for TKE
                                                 ! at time t-dt and at time t
  CHARACTER (LEN=4)  :: CGETRVT,CGETRCT,CGETRRT !  Get indicator for Rv
  CHARACTER (LEN=4)  :: CGETRIT,CGETRST,CGETRGT ! Rc,Rr,Ri,Rs,Rg,Rh
  CHARACTER (LEN=4)  :: CGETRHT                 ! at time t 
  CHARACTER (LEN=4)  :: CGETINPRC,CGETINPRR,CGETINPRS !  Get indicator for
                                                 ! 2D precip fields
  CHARACTER (LEN=4)  :: CGETINPRG,CGETINPRH
  CHARACTER (LEN=4)  :: CGETSVCONV ! Get indicator for
                                   ! SV deep convection (LCHTRANS)
  CHARACTER (LEN=4), DIMENSION(:), POINTER :: CGETSVT=>NULL() !  Get indicator 
                                ! for the Scalar Var. at time t
  CHARACTER (LEN=4) :: CGETLSUM, CGETLSVM, CGETLSWM   !  Get indicator for 
                                ! U,V,W for Larger Scales at time t-dt
  CHARACTER (LEN=4) :: CGETLSTHM, CGETLSRVM     !  Get indicator for 
                                ! Theta , Rv for Larger Scales at time t-dt
  CHARACTER (LEN=4)  :: CGETSIGS,CGETSRC      !  Get indicator for SIGS
                                ! and SRC related to the subgrid condensation
  CHARACTER (LEN=4)  :: CGETCLDFR             !  Get indicator for the
                                ! CLouD FRaction
  CHARACTER (LEN=4)  :: CGETSRCT              !  Get indicator for SRCM
                                ! and SRCT related to the subgrid condensation
  CHARACTER (LEN=4)  :: CGETCIT               !  Get indicator for the
                                                 ! primary ice concentration
  CHARACTER (LEN=4)  :: CGETCONV              ! Get indicator for the
                                ! use of a convection scheme
  CHARACTER (LEN=4)  :: CGETRAD               ! Get indicator for the
                                ! use of a radiation scheme
  CHARACTER (LEN=4)  :: CGETCLOUD             ! Get indicator for the
                                ! use of a cloud scheme
  CHARACTER (LEN=4)  :: CGETBL_DEPTH          ! Get indicator for the BL depth
  CHARACTER (LEN=4)  :: CGETSBL_DEPTH         ! Get indicator for the SBL depth
  CHARACTER (LEN=4)  :: CGETPHC,CGETPHR       ! Get indicator for the pH values
  CHARACTER (LEN=4)  :: CGETZWS
                                ! in cloud water and rainwater
!
END TYPE GET_t

TYPE(GET_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: GET_MODEL
LOGICAL    , DIMENSION(JPMODELMAX),         SAVE :: GET_FIRST_CALL = .TRUE.

CHARACTER (LEN=4), POINTER :: CGETUT=>NULL(), CGETVT=>NULL(), CGETWT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETTHT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETPABST=>NULL()
CHARACTER (LEN=4), POINTER :: CGETTKET=>NULL()
CHARACTER (LEN=4), POINTER :: CGETRVT=>NULL(),CGETRCT=>NULL(),CGETRRT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETRIT=>NULL(),CGETRST=>NULL(),CGETRGT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETRHT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETINPRC=>NULL(), CGETINPRR=>NULL(), CGETINPRS=>NULL()
CHARACTER (LEN=4), POINTER :: CGETINPRG=>NULL(), CGETINPRH=>NULL()
CHARACTER (LEN=4), POINTER :: CGETSVCONV=>NULL()
CHARACTER (LEN=4), DIMENSION(:), POINTER :: CGETSVT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETLSUM=>NULL(), CGETLSVM=>NULL(), CGETLSWM=>NULL()
CHARACTER (LEN=4), POINTER :: CGETLSTHM=>NULL(), CGETLSRVM=>NULL()
CHARACTER (LEN=4), POINTER :: CGETSIGS=>NULL(),CGETSRC=>NULL()
CHARACTER (LEN=4), POINTER :: CGETCLDFR=>NULL()
CHARACTER (LEN=4), POINTER :: CGETSRCT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETCIT=>NULL()
CHARACTER (LEN=4), POINTER :: CGETCONV=>NULL()
CHARACTER (LEN=4), POINTER :: CGETRAD=>NULL()
CHARACTER (LEN=4), POINTER :: CGETCLOUD=>NULL()
CHARACTER (LEN=4), POINTER :: CGETBL_DEPTH=>NULL()
CHARACTER (LEN=4), POINTER :: CGETSBL_DEPTH=>NULL()
CHARACTER (LEN=4), POINTER :: CGETPHC=>NULL()
CHARACTER (LEN=4), POINTER :: CGETPHR=>NULL()
CHARACTER (LEN=4), POINTER :: CGETZWS=>NULL()

CONTAINS

SUBROUTINE GET_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
!JUAN
IF (GET_FIRST_CALL(KTO)) THEN
ALLOCATE (GET_MODEL(KTO)%CGETSVT(JPSVMAX))
GET_FIRST_CALL(KTO) = .FALSE.
ENDIF
!JUAN
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
CGETUT=>GET_MODEL(KTO)%CGETUT
CGETVT=>GET_MODEL(KTO)%CGETVT
CGETWT=>GET_MODEL(KTO)%CGETWT
CGETTHT=>GET_MODEL(KTO)%CGETTHT
CGETPABST=>GET_MODEL(KTO)%CGETPABST
CGETTKET=>GET_MODEL(KTO)%CGETTKET
CGETRVT=>GET_MODEL(KTO)%CGETRVT
CGETRCT=>GET_MODEL(KTO)%CGETRCT
CGETRRT=>GET_MODEL(KTO)%CGETRRT
CGETRIT=>GET_MODEL(KTO)%CGETRIT
CGETRST=>GET_MODEL(KTO)%CGETRST
CGETRGT=>GET_MODEL(KTO)%CGETRGT
CGETRHT=>GET_MODEL(KTO)%CGETRHT
CGETINPRC=>GET_MODEL(KTO)%CGETINPRC
CGETINPRR=>GET_MODEL(KTO)%CGETINPRR
CGETINPRS=>GET_MODEL(KTO)%CGETINPRS
CGETINPRG=>GET_MODEL(KTO)%CGETINPRG
CGETINPRH=>GET_MODEL(KTO)%CGETINPRH
CGETSVCONV=>GET_MODEL(KTO)%CGETSVCONV
CGETSVT=>GET_MODEL(KTO)%CGETSVT
CGETLSUM=>GET_MODEL(KTO)%CGETLSUM
CGETLSVM=>GET_MODEL(KTO)%CGETLSVM
CGETLSWM=>GET_MODEL(KTO)%CGETLSWM
CGETLSTHM=>GET_MODEL(KTO)%CGETLSTHM
CGETLSRVM=>GET_MODEL(KTO)%CGETLSRVM
CGETSIGS=>GET_MODEL(KTO)%CGETSIGS
CGETSRC=>GET_MODEL(KTO)%CGETSRC
CGETCLDFR=>GET_MODEL(KTO)%CGETCLDFR
CGETSRCT=>GET_MODEL(KTO)%CGETSRCT
CGETCIT=>GET_MODEL(KTO)%CGETCIT
CGETZWS=>GET_MODEL(KTO)%CGETZWS
CGETCONV=>GET_MODEL(KTO)%CGETCONV
CGETRAD=>GET_MODEL(KTO)%CGETRAD
CGETCLOUD=>GET_MODEL(KTO)%CGETCLOUD
CGETBL_DEPTH=>GET_MODEL(KTO)%CGETBL_DEPTH
CGETSBL_DEPTH=>GET_MODEL(KTO)%CGETSBL_DEPTH
CGETPHC=>GET_MODEL(KTO)%CGETPHC
CGETPHR=>GET_MODEL(KTO)%CGETPHR

END SUBROUTINE GET_GOTO_MODEL

END MODULE MODD_GET_n
