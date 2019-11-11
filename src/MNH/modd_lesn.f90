!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ############### 
      MODULE MODD_LES_n
!     ###############
!
!!****  *MODD_LES* - declaration of prognostic variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     resolved fluxes and the spectra computed in LES mode
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------ 
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_LES)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!
!!    AUTHOR
!!    ------
!!	   J. Cuxart   *INM and Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    March 10, 1995
!!
!!       (J.Stein)  Sept. 25, 1995  add the model number in LES mode
!!       J. Cuxart  Oct.   4, 1996  New time series
!!       V. Masson  Jan.  20, 2000  New LES routines variables & //
!!       V. Masson  Nov.   6, 2002  LES budgets
!!       O.Thouron  June,     2008  New radiation diagnostics
!!                    10/2016 (C.Lac) Add droplet deposition
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!!                    02/2019 (C. Lac) Add rain fraction as a LES diagnostic
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE LES_t
!
!-------------------------------------------------------------------------------
!
!* monitoring variables
!
  INTEGER :: NLES_TIMES         ! number of LES computations in time
  INTEGER :: NLES_DTCOUNT       ! number of time steps between two LES comp.
  INTEGER :: NLES_TCOUNT        ! current time counter for LES comp.
  INTEGER :: NSPECTRA_NI        ! number of wave lengths in I direction
  INTEGER :: NSPECTRA_NJ        ! number of wave lengths in J direction
!
  REAL, DIMENSION(:,:), POINTER :: XLES_DATIME=>NULL() !  date array for diachro
  REAL, DIMENSION(:,:), POINTER :: XLES_TRAJT=>NULL()  !  sampling times array for diachro
!
  REAL, DIMENSION(:),   POINTER :: XLES_Z=>NULL()      ! altitudes
  REAL                              :: XLES_ZS     ! mean orography
!-------------------------------------------------------------------------------
!
  REAL,    DIMENSION(:,:,:), POINTER :: XCOEFLIN_LES=>NULL()
! coefficients for vertical interpolation
!
  INTEGER, DIMENSION(:,:,:), POINTER :: NKLIN_LES=>NULL()
! levels for vertical interpolation
!
  REAL,    DIMENSION(:,:,:), POINTER :: XCOEFLIN_SPEC=>NULL()
! coefficients for vertical interpolation
!
  INTEGER, DIMENSION(:,:,:), POINTER :: NKLIN_SPEC=>NULL()
! levels for vertical interpolation
!
!-------------------------------------------------------------------------------
!
!* averaging control quantities
!
  INTEGER, DIMENSION(:,:,:), POINTER :: NLES_AVG_PTS_ll=>NULL() ! number of points
!                                                         ! used for averaging
  INTEGER, DIMENSION(:,:,:), POINTER :: NLES_UND_PTS_ll=>NULL() ! number of points
!                                                         ! not used for averaging
!                                                         ! because they are
!                                                         ! intersecting
!                                                         ! the orography
!
!-------------------------------------------------------------------------------
!
!* mean quantities
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_U=>NULL()   ! <u>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_V=>NULL()   ! <v>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_W=>NULL()   ! <w>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_P=>NULL()   ! <p>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_RHO=>NULL() ! <density>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Th=>NULL()  ! <theta>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Thv=>NULL() ! <thetav>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Thl=>NULL() ! <thetal>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rt=>NULL()  ! <Rt>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rv=>NULL()  ! <Rv>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rehu=>NULL()  ! <Rh>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Qs=>NULL()  ! saturated spec h 
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rc=>NULL()  ! <Rc>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Cf=>NULL()  ! <CLDFR>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_INDCf=>NULL()  ! <Cf> tq rc>0 (0 OU 1)
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_INDCf2=>NULL()  ! <Cf> tq rc>1E-5 (0 OU 1)
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_RF=>NULL()  ! <RAINFR>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Mf=>NULL()  ! <Mf> 
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_KHt=>NULL()! <Kh for thet
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_KHr=>NULL()! <Kh for qr>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rr=>NULL()  ! <Rr>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Ri=>NULL()  ! <Ri>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rs=>NULL()  ! <Rs>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rg=>NULL()  ! <Rg>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rh=>NULL()  ! <Rh>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_MEAN_Sv=>NULL()! <Sv>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_WIND=>NULL()! <wind modulus>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dUdz=>NULL()   ! <dudz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dVdz=>NULL()   ! <dvdz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dWdz=>NULL()   ! <dwdz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dThldz=>NULL() ! <dThldz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dRtdz=>NULL()  ! <dRtdz>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_MEAN_dSvdz=>NULL()! <dSvdz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_DP=>NULL()   ! <dp>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_TP=>NULL()   ! <tp>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_TR=>NULL()   ! <tr>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_DISS=>NULL()   ! <diss>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_LM=>NULL()   ! <lm>
!
!-------------------------------------------------------------------------------
!
!* resolved quantities
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_U2=>NULL()    ! <u'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_V2=>NULL()    ! <v'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2=>NULL()    ! <w'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_P2=>NULL()    ! <p'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Th2=>NULL()   ! <Th'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Thl2=>NULL()  ! <Thl'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThThv=>NULL() ! <Th'Thv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlThv=>NULL()! <Thl'Thv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Ke=>NULL()    ! 0.5 <u'2+v'2+w'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Rv2=>NULL()   ! <Rv'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Rc2=>NULL()   ! <Rc'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Ri2=>NULL()   ! <Ri'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Rt2=>NULL()   ! <Rt'2>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_Sv2=>NULL() ! <Sv'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThRv=>NULL()  ! <Th'Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThRc=>NULL()  ! <Th'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThRi=>NULL()  ! <Th'Ri'>
  REAL, DIMENSION(:,:,:,:), POINTER ::XLES_RESOLVED_ThSv=>NULL() ! <ThSv>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlRv=>NULL() ! <Thl'Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlRc=>NULL() ! <Thl'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlRi=>NULL() ! <Thl'Ri'>
  REAL, DIMENSION(:,:,:,:), POINTER ::XLES_RESOLVED_ThlSv=>NULL()! <ThlSv>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThvRv=>NULL() ! <Thv'Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThvRc=>NULL() ! <Thv'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThvRi=>NULL() ! <Thv'Ri'>
  REAL, DIMENSION(:,:,:,:), POINTER ::XLES_RESOLVED_ThvSv=>NULL()! <ThvSv>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UV=>NULL()    ! <u'v'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WU=>NULL()    ! <w'u'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WV=>NULL()    ! <w'v'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UP=>NULL()    ! <u'p'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VP=>NULL()    ! <v'p'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WP=>NULL()    ! <w'p'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UTh=>NULL()   ! <u'Th'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VTh=>NULL()   ! <v'Th'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WTh=>NULL()   ! <w'Th'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UThl=>NULL()  ! <u'Thl'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VThl=>NULL()  ! <v'Thl'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThl=>NULL()  ! <w'Thl'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UThv=>NULL()  ! <u'Thv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VThv=>NULL()  ! <v'Thv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThv=>NULL()  ! <w'Thv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_URv=>NULL()   ! <u'Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VRv=>NULL()   ! <v'Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRv=>NULL()   ! <w'Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_URc=>NULL()   ! <u'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VRc=>NULL()   ! <v'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRc=>NULL()   ! <w'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_URi=>NULL()   ! <u'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VRi=>NULL()   ! <v'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRi=>NULL()   ! <w'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRt=>NULL()   ! <w'Rt'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRr=>NULL()   ! <w'Rr'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_INPRR3D=>NULL() ! flux precip
  REAL, DIMENSION(:,:,:), POINTER :: XLES_MAX_INPRR3D=>NULL() ! mass flux precip
  REAL, DIMENSION(:,:,:), POINTER :: XLES_EVAP3D=>NULL() ! evaporation
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_USv=>NULL() ! <u'Sv'>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_VSv=>NULL() ! <v'Sv'>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WSv=>NULL() ! <w'Sv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_U3=>NULL()    ! <u'3>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_V3=>NULL()    ! <v'3>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W3=>NULL()    ! <w'3>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_U4=>NULL()    ! <u'4>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_V4=>NULL()    ! <v'4>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W4=>NULL()    ! <w'4>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Ua_ddxa_P=>NULL() ! <ua'dp'/dxa>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlPz=>NULL() ! <Thl'dp'/dz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThl2=>NULL() ! <w'Thl'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Thl=>NULL() ! <w'2Thl'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRv2=>NULL()  ! <w'Rv'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Rv=>NULL()  ! <w'2Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRt2=>NULL() ! <w'Rt'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Rt=>NULL() ! <w'2Rt'>      
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RvPz=>NULL()  ! <Rv'dp'/dz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RtPz=>NULL()  ! <Rt'dp'/dz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRv=>NULL()! <w'Thl'Rv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRt=>NULL()! <w'Thl'Rt'>      
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRc2=>NULL()  ! <w'Rc'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Rc=>NULL()  ! <w'2Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RCPz=>NULL()  ! <Rc'dp'/dz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRc=>NULL()! <w'Thl'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRvRc=>NULL() ! <w'Rv'Rc'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRi2=>NULL()  ! <w'Ri'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Ri=>NULL()  ! <w'2Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RIPz=>NULL()  ! <Ri'dp'/dz>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRi=>NULL()! <w'Thl'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRvRi=>NULL() ! <w'Rv'Ri'>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WSv2=>NULL()! <w'Sv'2>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_W2Sv=>NULL()! <w'2Sv'>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_SvPz=>NULL()! <Sv'dp'/dz>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WThlSv=>NULL()! <w'Thl'Sv'>
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WRvSv=>NULL() ! <w'Rv'Sv'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_MASSFX=>NULL()! <upward mass flux>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UKe=>NULL()   ! <u'E>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VKe=>NULL()   ! <v'E>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WKe=>NULL()   ! <w'E>
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_KE=>NULL()    ! budget terms
  REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_WThl=>NULL()  ! (including resolved
  REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_Thl2=>NULL()  !  turbulent fluxes)
  REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_WRt=>NULL()   ! for the resolved
  REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_Rt2=>NULL()   ! fluxes, variances
  REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_ThlRt=>NULL() ! and covariances
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_BU_RES_Sv2=>NULL() ! 
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_BU_RES_WSv=>NULL() !
!
!-------------------------------------------------------------------------------
!
!* mixing of resolved and subgrid quantities
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_U_SBG_Tke=>NULL()  ! <u'Tke>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_V_SBG_Tke=>NULL()  ! <v'Tke>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_Tke=>NULL()  ! <w'Tke>
!                                                                ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_WThl=>NULL() !  <w'w'Thl'>
!                                                                _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_WRt=>NULL()  !  <w'w'Rt'>
!                                                                _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_Thl2=>NULL() !  <w'Thl'2>
!                                                                ____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_Rt2=>NULL()  !  <w'Rt'2>
!                                                                _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_ThlRt=>NULL()!  <w'Thl'Rt'>
!                                                                ____
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_RES_W_SBG_WSv=>NULL()  !  <w'w'Sv'>
!                                                                ____
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_RES_W_SBG_Sv2=>NULL()  !  <w'Sv'2>
!                                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_U_SBG_UaU=>NULL()   !  <du'/dxa ua'u'>
!                                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_V_SBG_UaV=>NULL()   !  <dv'/dxa ua'v'>
!                                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaW=>NULL()   !  <dw'/dxa ua'w'>
!                                                                            _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaThl=>NULL() !  <dw'/dxa ua'Thl'>
!                                                                              _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Thl_SBG_UaW=>NULL() !  <dThl'/dxa ua'w'>
!                                                                             ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddz_Thl_SBG_W2=>NULL()   !  <dThl'/dz w'2>
!                                                                            ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaRt=>NULL()  !  <dw'/dxa ua'Rt'>
!                                                                             _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Rt_SBG_UaW=>NULL()  !  <dRt'/dxa ua'w'>
!                                                                            ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddz_Rt_SBG_W2=>NULL()    !  <dRt'/dz w'2>
!                                                                              ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Thl_SBG_UaRt=>NULL()!  <dThl'/dxa ua'Rt'>
!                                                                             _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Rt_SBG_UaThl=>NULL()!  <dRt'/dxa ua'Thl'>
!                                                                              _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Thl_SBG_UaThl=>NULL()! <dThl'/dxa ua'Thl'>
!                                                                             ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Rt_SBG_UaRt=>NULL() !  <dRt'/dxa ua'Rt'>
!                                                                              ______
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaSv=>NULL()  !  <dw'/dxa ua'Sv'>
!                                                                               _____
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddxa_Sv_SBG_UaW=>NULL()  !  <dSv'/dxa ua'w'>
!                                                                              ___
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddz_Sv_SBG_W2=>NULL()    !  <dSv'/dz w'2>
!                                                                               ______
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddxa_Sv_SBG_UaSv=>NULL() !  <dSv'/dxa ua'Sv'>
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_SBG_TKE=>NULL() ! budget terms of Tke
!-------------------------------------------------------------------------------
!
!* subgrid quantities
!                                                            ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_U2=>NULL()    ! <u'2>
!                                                            ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_V2=>NULL()    ! <v'2>
!                                                            ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_W2=>NULL()    ! <w'2>
!                                                            _
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Tke=>NULL()   ! <e>
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Thl2=>NULL()  ! <Thl'2>
!                                                            ____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Rt2=>NULL()   ! <Rt'2>
!                                                            ____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Rc2=>NULL()   ! <Rc'2>
!                                                            ____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Ri2=>NULL()   ! <Ri'2>
!                                                            _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ThlRt=>NULL() ! <Thl'Rt'>
!                                                            ____
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_Sv2=>NULL()   ! <Sv'2>
!                                                            ____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_UV=>NULL()    ! <u'v'>
!                                                            ____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WU=>NULL()    ! <w'u'>
!                                                            ____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WV=>NULL()    ! <w'v'>
!                                                            ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_UThl=>NULL()  ! <u'Thl'>
!                                                            ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VThl=>NULL()  ! <v'Thl'>
!                                                            ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThl=>NULL()  ! <w'Thl'>
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ENTLES=>NULL()  ! Up entr flux
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DETLES=>NULL()  ! Up detr flux      
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_URt=>NULL()   ! <u'Rt'>
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VRt=>NULL()   ! <v'Rt'>
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRt=>NULL()   ! <w'Rt'>
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_URc=>NULL()   ! <u'Rc'>
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VRc=>NULL()   ! <v'Rc'>
!                                                            _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRc=>NULL()   ! <w'Rc'>
!                                                            _____
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_SUBGRID_USv=>NULL() ! <u'Sv'>
!                                                            _____
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_SUBGRID_VSv=>NULL() ! <v'Sv'>
!                                                            _____
  REAL, DIMENSION(:,:,:,:), POINTER :: XLES_SUBGRID_WSv=>NULL() ! <w'Sv'>
!                                                            ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_UTke=>NULL()  ! <u'e>
!                                                            ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VTke=>NULL()  ! <v'e>
!                                                            ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WTke=>NULL()  ! <w'e>
!                                                                 ___
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ddz_WTke=>NULL()  ! <dw'e/dz>
!                                                             ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThv=>NULL()   ! <w'Thv'>
!                                                             ________
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ThlThv=>NULL() ! <Thl'Thv'>
!                                                             _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RtThv=>NULL()  ! <Rt'Thv'>
!                                                             _______
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_SvThv=>NULL()  ! <Sv'Thv'>
!                                                             ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_W2Thl=>NULL()  ! <w'2Thl>
!                                                             _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_W2Rt=>NULL()   ! <w'2Rt>
!                                                             _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThlRt=>NULL() ! <w'ThlRt>
!                                                             _____
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_W2Sv=>NULL()   ! <w'2Sv>
!                                                             ______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThl2=>NULL()  ! <w'Thl2>
!                                                             _____
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRt2=>NULL()   ! <w'Rt2>
!                                                             _____
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_WSv2=>NULL()   ! <w'Sv2>
!                                                               _______
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_Tke=>NULL() ! <epsilon>
!                                                                ____________
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_Thl2=>NULL() ! <epsilon_Thl2>
!                                                                ___________
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_Rt2=>NULL()  ! <epsilon_Rt2>
!                                                                _____________
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_ThlRt=>NULL()! <epsilon_ThlRt>
!                                                                ___________
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_DISS_Sv2=>NULL()  ! <epsilon_Sv2>
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WP=>NULL()     ! <w'p'>
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ThlPz=>NULL()  ! <Thl'dp'/dz>
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RtPz=>NULL()   ! <Rt'dp'/dz>
!
  REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_SvPz=>NULL()   ! <Sv'dp'/dz>
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_PHI3=>NULL()   ! phi3
!      
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_PSI3=>NULL()   ! psi3
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_LMix=>NULL()   ! mixing length
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_LDiss=>NULL()  ! dissipative length
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Km=>NULL()     ! eddy diffusivity for momentum
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Kh=>NULL()     ! eddy diffusivity for heat
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_THLUP_MF=>NULL()  ! Thl of the Updraft
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RTUP_MF=>NULL()   ! Rt of the Updraft
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RVUP_MF=>NULL()   ! Rv of the Updraft
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RCUP_MF=>NULL()   ! Rc of the Updraft
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RIUP_MF=>NULL()   ! Ri of the Updraft
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WUP_MF=>NULL()    ! Thl of the Updraft
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_MASSFLUX=>NULL()  ! Mass Flux
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DETR=>NULL()      ! Detrainment
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ENTR=>NULL()      ! Entrainment
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_FRACUP=>NULL()    ! Updraft Fraction 
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_THVUP_MF=>NULL()  ! Thv of the Updraft
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WTHLMF=>NULL() ! Flux of thl   
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRTMF=>NULL()  ! Flux of rt
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WTHVMF=>NULL() ! Flux of thv 
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WUMF=>NULL()   ! Flux of u
!
  REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WVMF=>NULL()   ! Flux of v

!-------------------------------------------------------------------------------
!
!* updraft quantities
!
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT=>NULL()       ! updraft fraction
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_W=>NULL()     ! <w>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Th=>NULL()    ! <theta>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Thv=>NULL()   ! <thetav>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Thl=>NULL()   ! <thetal>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rv=>NULL()    ! <Rv>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rc=>NULL()    ! <Rc>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rr=>NULL()    ! <Rr>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Ri=>NULL()    ! <Ri>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rs=>NULL()    ! <Rs>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rg=>NULL()    ! <Rg>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rh=>NULL()    ! <Rh>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_Sv=>NULL()  ! <Sv>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Tke=>NULL()   ! <e>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Ke=>NULL()    ! <E>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WTh=>NULL()   ! <w'theta'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WThl=>NULL()  ! <w'thetal'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WThv=>NULL()  ! <w'thetav'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WRv=>NULL()   ! <w'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WRc=>NULL()   ! <w'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WRi=>NULL()   ! <w'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_WSv=>NULL() ! <w'Sv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Th2=>NULL()   ! <th'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Thl2=>NULL()  ! <thl'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThThv=>NULL() ! <th'thv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlThv=>NULL()! <thl'thv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rv2=>NULL()   ! <Rv'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rc2=>NULL()   ! <Rc'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Ri2=>NULL()   ! <Ri'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_Sv2=>NULL() ! <Sv'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThRv=>NULL()  ! <Th'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThRc=>NULL()  ! <Th'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThRi=>NULL()  ! <Th'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_ThSv=>NULL()  ! <Th'Sv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThvRv=>NULL() ! <Thv'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThvRc=>NULL() ! <Thv'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThvRi=>NULL() ! <Thv'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_ThvSv=>NULL() ! <Thv'Sv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlRv=>NULL() ! <Thl'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlRc=>NULL() ! <Thl'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlRi=>NULL() ! <Thl'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_ThlSv=>NULL() ! <Thl'Sv'>
!
!-------------------------------------------------------------------------------
!
!* downdraft quantities
!

  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT=>NULL()       ! updraft fraction
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_W=>NULL()     ! <w>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Th=>NULL()    ! <theta>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Thv=>NULL()   ! <thetav>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Thl=>NULL()   ! <thetal>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rv=>NULL()    ! <Rv>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rc=>NULL()    ! <Rc>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rr=>NULL()    ! <Rr>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Ri=>NULL()    ! <Ri>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rs=>NULL()    ! <Rs>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rg=>NULL()    ! <Rg>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rh=>NULL()    ! <Rh>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_Sv=>NULL()  ! <Sv>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Tke=>NULL()   ! <e>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Ke=>NULL()    ! <E>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WTh=>NULL()   ! <w'theta'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WThl=>NULL()  ! <w'thetal'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WThv=>NULL()  ! <w'thetav'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WRv=>NULL()   ! <w'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WRc=>NULL()   ! <w'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WRi=>NULL()   ! <w'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_WSv=>NULL() ! <w'Sv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Th2=>NULL()   ! <th'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Thl2=>NULL()  ! <thl'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThThv=>NULL() ! <th'thv>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlThv=>NULL()! <thl'thv>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rv2=>NULL()   ! <Rv'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rc2=>NULL()   ! <Rc'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Ri2=>NULL()   ! <Ri'2>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_Sv2=>NULL() ! <Sv'2>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThRv=>NULL()  ! <Th'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThRc=>NULL()  ! <Th'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThRi=>NULL()  ! <Th'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_ThSv=>NULL()  ! <Th'Sv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThvRv=>NULL() ! <Thv'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThvRc=>NULL() ! <Thv'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThvRi=>NULL() ! <Thv'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_ThvSv=>NULL() ! <Thv'Sv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlRv=>NULL() ! <Thl'Rv'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlRc=>NULL() ! <Thl'Rc'>
  REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlRi=>NULL() ! <Thl'Ri'>
  REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_ThlSv=>NULL() ! <Thl'Sv'>
!
!-------------------------------------------------------------------------------
!
!* surface variables
!
  REAL, DIMENSION(:), POINTER :: XLES_UW0=>NULL()       ! uw temporal series
  REAL, DIMENSION(:), POINTER :: XLES_VW0=>NULL()       ! vw temporal series
  REAL, DIMENSION(:), POINTER :: XLES_USTAR=>NULL()     ! u* temporal series
  REAL, DIMENSION(:), POINTER :: XLES_WSTAR=>NULL()     ! w* temporal series
  REAL, DIMENSION(:), POINTER :: XLES_Q0=>NULL()        ! Qo temporal series
  REAL, DIMENSION(:), POINTER :: XLES_E0=>NULL()        ! Eo temporal series
  REAL, DIMENSION(:,:), POINTER :: XLES_SV0=>NULL()     ! scalar surface fluxes
  REAL, DIMENSION(:), POINTER :: XLES_BL_HEIGHT=>NULL() ! BL height
  REAL, DIMENSION(:), POINTER :: XLES_MO_LENGTH=>NULL() ! Monin-Obukhov Length
  REAL, DIMENSION(:), POINTER :: XLES_ZCB=>NULL()       ! cloud base from rc
  REAL, DIMENSION(:), POINTER :: XLES_CFtot=>NULL()     ! total cf 
  REAL, DIMENSION(:), POINTER :: XLES_CF2tot=>NULL()    ! total cf2
  REAL, DIMENSION(:), POINTER :: XLES_LWP=>NULL()       ! lwpath
  REAL, DIMENSION(:), POINTER :: XLES_LWPVAR=>NULL()    ! lwpath variance
  REAL, DIMENSION(:), POINTER :: XLES_RWP=>NULL()       ! rain w path
  REAL, DIMENSION(:), POINTER :: XLES_IWP=>NULL()       ! ice  w path
  REAL, DIMENSION(:), POINTER :: XLES_SWP=>NULL()       ! snow w path
  REAL, DIMENSION(:), POINTER :: XLES_GWP=>NULL()       ! graupel w path
  REAL, DIMENSION(:), POINTER :: XLES_HWP=>NULL()       ! hail w path
  REAL, DIMENSION(:), POINTER :: XLES_INT_TKE=>NULL()   ! vert. integratedtke
  REAL, DIMENSION(:), POINTER :: XLES_INPRR=>NULL()     ! inst precip rate
  REAL, DIMENSION(:), POINTER :: XLES_INPRC=>NULL()     ! inst cloud precip rate
  REAL, DIMENSION(:), POINTER :: XLES_INDEP=>NULL()     ! inst cloud deposition rate
  REAL, DIMENSION(:), POINTER :: XLES_RAIN_INPRR=>NULL()!flux prec rainy cell
  REAL, DIMENSION(:), POINTER :: XLES_ACPRR=>NULL()     ! acc precip rate
  REAL, DIMENSION(:), POINTER :: XLES_ZMAXCF=>NULL()    ! z de max cf
  REAL, DIMENSION(:), POINTER :: XLES_ZMAXCF2=>NULL()   ! z de max de cf2
  REAL, DIMENSION(:), POINTER :: XLES_PRECFR=>NULL()    !% col with sedflux>cr      
!
!-------------------------------------------------------------------------------
!
!* 2-points correlations in I direction
!
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_UU=>NULL()    ! between u and u
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_VV=>NULL()    ! between v and v
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_UV=>NULL()    ! between u and v
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WU=>NULL()    ! between w and u
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WV=>NULL()    ! between w and v
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WW=>NULL()    ! between w and w
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WTh=>NULL()   ! between w and theta
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WThl=>NULL()  ! between w and thetal
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WRv=>NULL()   ! between w and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WRc=>NULL()   ! between w and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WRi=>NULL()   ! between w and Ri
  REAL, DIMENSION(:,:,:,:), POINTER :: XCORRi_WSv=>NULL() ! between w and Sv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThTh=>NULL()  ! between theta and theta
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThRv=>NULL()  ! between theta and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThRc=>NULL()  ! between theta and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThRi=>NULL()  ! between theta and Ri
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlThl=>NULL()! between thetal and thetal
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlRv=>NULL() ! between thetal and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlRc=>NULL() ! between thetal and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlRi=>NULL() ! between thetal and Ri
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_RvRv=>NULL()  ! between Rv and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_RcRc=>NULL()  ! between Rc and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRi_RiRi=>NULL()  ! between Ri and Ri
  REAL, DIMENSION(:,:,:,:), POINTER :: XCORRi_SvSv=>NULL()! between Sv and Sv
!
!-------------------------------------------------------------------------------
!
!* 2-points correlations in J direction
!
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_UU=>NULL()    ! between u and u
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_VV=>NULL()    ! between v and v
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_UV=>NULL()    ! between u and v
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WU=>NULL()    ! between w and u
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WV=>NULL()    ! between w and v
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WW=>NULL()    ! between w and w
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WTh=>NULL()   ! between w and theta
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WThl=>NULL()  ! between w and thetal
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WRv=>NULL()   ! between w and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WRc=>NULL()   ! between w and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WRi=>NULL()   ! between w and Ri
  REAL, DIMENSION(:,:,:,:), POINTER :: XCORRj_WSv=>NULL() ! between w and Sv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThTh=>NULL()  ! between theta and theta
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThRv=>NULL()  ! between theta and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThRc=>NULL()  ! between theta and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThRi=>NULL()  ! between theta and Ri
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlThl=>NULL()! between thetal and thetal
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlRv=>NULL() ! between thetal and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlRc=>NULL() ! between thetal and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlRi=>NULL() ! between thetal and Ri
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_RvRv=>NULL()  ! between Rv and Rv
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_RcRc=>NULL()  ! between Rc and Rc
  REAL, DIMENSION(:,:,:), POINTER :: XCORRj_RiRi=>NULL()  ! between Ri and Ri
  REAL, DIMENSION(:,:,:,:), POINTER :: XCORRj_SvSv=>NULL()! between Sv and Sv
!
!-------------------------------------------------------------------------------
!
!* 2-points correlations in K direction
!
  REAL, DIMENSION(:,:,:), POINTER :: XCORRk_WW=>NULL()   ! between w and w
!
!lw and sw fluxes up and down
  REAL, DIMENSION(:,:), POINTER :: XLES_SWU => NULL()   !mean on the domain of the sw_up flux
                                                 !temporal serie
  REAL, DIMENSION(:,:), POINTER :: XLES_SWD => NULL()   !mean on the domain of the sw_down flux
                                                 !temporal serie
  REAL, DIMENSION(:,:), POINTER :: XLES_LWU => NULL()   !mean on the domain of the lw_up flux
                                                 !temporal serie
  REAL, DIMENSION(:,:), POINTER :: XLES_LWD => NULL()   !mean on the domain of the lw_down flux
                                                 !temporal serie
  REAL, DIMENSION(:,:), POINTER :: XLES_DTHRADSW => NULL()   !mean on the domain of dthrad_sw flux
                                                 !temporal serie
  REAL, DIMENSION(:,:), POINTER :: XLES_DTHRADLW => NULL()   !mean on the domain of dthrad_lw flux
                                                 !temporal serie
  REAL, DIMENSION(:,:), POINTER :: XLES_RADEFF => NULL()   ! effective radius
!-------------------------------------------------------------------------------
!
END TYPE LES_t

TYPE(LES_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: LES_MODEL

INTEGER, POINTER :: NLES_TIMES=>NULL()
INTEGER, POINTER :: NLES_DTCOUNT=>NULL()
INTEGER, POINTER :: NLES_TCOUNT=>NULL()
INTEGER, POINTER :: NSPECTRA_NI=>NULL()
INTEGER, POINTER :: NSPECTRA_NJ=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DATIME=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_TRAJT=>NULL()
REAL, DIMENSION(:),   POINTER :: XLES_Z=>NULL()
REAL, POINTER :: XLES_ZS=>NULL()
REAL,    DIMENSION(:,:,:), POINTER :: XCOEFLIN_LES=>NULL()
INTEGER, DIMENSION(:,:,:), POINTER :: NKLIN_LES=>NULL()
REAL,    DIMENSION(:,:,:), POINTER :: XCOEFLIN_SPEC=>NULL()
INTEGER, DIMENSION(:,:,:), POINTER :: NKLIN_SPEC=>NULL()
INTEGER, DIMENSION(:,:,:), POINTER :: NLES_AVG_PTS_ll=>NULL()
INTEGER, DIMENSION(:,:,:), POINTER :: NLES_UND_PTS_ll=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_U=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_V=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_W=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_P=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_RHO=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Th=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Thv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Thl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rehu=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Qs=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Cf=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_INDCf=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_INDCf2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_RF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Mf=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_KHt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_KHr=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rr=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Ri=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rs=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rg=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_Rh=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_MEAN_Sv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_WIND=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dUdz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dVdz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dWdz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dThldz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_dRtdz=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_MEAN_dSvdz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_DP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_TP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_TR=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_DISS=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MEAN_LM=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_U2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_V2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_P2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Th2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Thl2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Ke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Rv2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Rc2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Ri2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Rt2=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_Sv2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_ThSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_ThlSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThvRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThvRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThvRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_ThvSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UTh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VTh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WTh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_URv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_URc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_URi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRr=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_INPRR3D=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_MAX_INPRR3D=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_EVAP3D=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_USv=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_VSv=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_U3=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_V3=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W3=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_U4=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_V4=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W4=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_Ua_ddxa_P=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_ThlPz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThl2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Thl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRv2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Rv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRt2=>NULL() 
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Rt=>NULL()  
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RvPz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RtPz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRc2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Rc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RCPz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRvRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRi2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_W2Ri=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_RIPz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WThlRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WRvRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WSv2=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_W2Sv=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_SvPz=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WThlSv=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RESOLVED_WRvSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_MASSFX=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_UKe=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_VKe=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RESOLVED_WKe=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_KE=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_WThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_Thl2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_WRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_Rt2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_RES_ThlRt=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_BU_RES_Sv2=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_BU_RES_WSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_U_SBG_Tke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_V_SBG_Tke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_Tke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_WThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_WRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_Thl2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_Rt2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_W_SBG_ThlRt=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_RES_W_SBG_WSv=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_RES_W_SBG_Sv2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_U_SBG_UaU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_V_SBG_UaV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaW=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Thl_SBG_UaW=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddz_Thl_SBG_W2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Rt_SBG_UaW=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddz_Rt_SBG_W2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Thl_SBG_UaRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Rt_SBG_UaThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Thl_SBG_UaThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_RES_ddxa_Rt_SBG_UaRt=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddxa_W_SBG_UaSv=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddxa_Sv_SBG_UaW=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddz_Sv_SBG_W2=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_RES_ddxa_Sv_SBG_UaSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_BU_SBG_TKE=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_U2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_V2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_W2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Tke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Thl2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Rt2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Rc2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_Ri2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ThlRt=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_Sv2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_UV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_UThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ENTLES=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DETLES=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_URt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_URc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRc=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_SUBGRID_USv=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_SUBGRID_VSv=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XLES_SUBGRID_WSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_UTke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_VTke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WTke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ddz_WTke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ThlThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RtThv=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_SvThv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_W2Thl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_W2Rt=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThlRt=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_W2Sv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WThl2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRt2=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_WSv2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_Tke=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_Thl2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_Rt2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DISS_ThlRt=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_DISS_Sv2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WP=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ThlPz=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RtPz=>NULL()
REAL, DIMENSION(:,:,:,:),POINTER:: XLES_SUBGRID_SvPz=>NULL()
REAL, DIMENSION(:,:,:),  POINTER:: XLES_SUBGRID_PHI3=>NULL()
REAL, DIMENSION(:,:,:),  POINTER:: XLES_SUBGRID_PSI3=>NULL()
REAL, DIMENSION(:,:,:),  POINTER:: XLES_SUBGRID_LMix=>NULL()
REAL, DIMENSION(:,:,:),  POINTER:: XLES_SUBGRID_LDiss=>NULL()
REAL, DIMENSION(:,:,:),  POINTER:: XLES_SUBGRID_Km=>NULL()
REAL, DIMENSION(:,:,:),  POINTER:: XLES_SUBGRID_Kh=>NULL()

REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_THLUP_MF=>NULL()  ! Thl of the Updraft
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RTUP_MF=>NULL()   ! Rt of the Updraft
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RVUP_MF=>NULL()   ! Rv of the Updraft
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RCUP_MF=>NULL()   ! Rc of the Updraft
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_RIUP_MF=>NULL()   ! Ri of the Updraft
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WUP_MF=>NULL()    ! Thl of the Updraft
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_MASSFLUX=>NULL()  ! Mass Flux
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_DETR=>NULL()      ! Detrainment
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_ENTR=>NULL()      ! Entrainment
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_FRACUP=>NULL()    ! Updraft Fraction 
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_THVUP_MF=>NULL()  ! Thv of the Updraft
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WTHLMF=>NULL() ! Flux of thl   
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WRTMF=>NULL()  ! Flux of rt
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WTHVMF=>NULL() ! Flux of thv 
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WUMF=>NULL()   ! Flux of u
REAL, DIMENSION(:,:,:), POINTER :: XLES_SUBGRID_WVMF=>NULL()   ! Flux of v

REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_W=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Th=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Thv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Thl=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rr=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Ri=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rs=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rg=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_Sv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Tke=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Ke=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WTh=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WThl=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WThv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_WRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_WSv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Th2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Thl2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThThv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlThv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rv2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Rc2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_Ri2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_Sv2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_ThSv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThvRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThvRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThvRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_ThvSv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_UPDRAFT_ThlRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_UPDRAFT_ThlSv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_W=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Th=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Thv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Thl=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rr=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Ri=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rs=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rg=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_Sv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Tke=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Ke=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WTh=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WThl=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WThv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_WRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_WSv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Th2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Thl2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThThv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlThv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rv2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Rc2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_Ri2=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_Sv2=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_ThSv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThvRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThvRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThvRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_ThvSv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlRv=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlRc=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DOWNDRAFT_ThlRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XLES_DOWNDRAFT_ThlSv=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_UW0=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_VW0=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_USTAR=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_WSTAR=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_Q0=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_E0=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_SV0=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_BL_HEIGHT=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_MO_LENGTH=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_ZCB=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_CFtot=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_CF2tot=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_LWP=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_LWPVAR=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_RWP=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_IWP=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_SWP=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_GWP=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_HWP=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_INT_TKE=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_INPRR=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_INPRC=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_INDEP=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_RAIN_INPRR=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_ACPRR=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_ZMAXCF=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_ZMAXCF2=>NULL()
REAL, DIMENSION(:), POINTER :: XLES_PRECFR=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_UU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_VV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_UV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WW=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WTh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_WRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XCORRi_WSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThTh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_ThlRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_RvRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_RcRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRi_RiRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XCORRi_SvSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_UU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_VV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_UV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WU=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WV=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WW=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WTh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_WRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XCORRj_WSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThTh=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlThl=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_ThlRi=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_RvRv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_RcRc=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRj_RiRi=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XCORRj_SvSv=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCORRk_WW=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_SWU=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_SWD=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_LWU=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_LWD=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DTHRADSW=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_DTHRADLW=>NULL()
REAL, DIMENSION(:,:), POINTER :: XLES_RADEFF=>NULL()

CONTAINS

SUBROUTINE LES_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
LES_MODEL(KFROM)%XLES_DATIME=>XLES_DATIME
LES_MODEL(KFROM)%XLES_TRAJT=>XLES_TRAJT
LES_MODEL(KFROM)%XLES_Z=>XLES_Z
LES_MODEL(KFROM)%XCOEFLIN_LES=>XCOEFLIN_LES
LES_MODEL(KFROM)%NKLIN_LES=>NKLIN_LES
LES_MODEL(KFROM)%XCOEFLIN_SPEC=>XCOEFLIN_SPEC
LES_MODEL(KFROM)%NKLIN_SPEC=>NKLIN_SPEC
LES_MODEL(KFROM)%NLES_AVG_PTS_ll=>NLES_AVG_PTS_ll
LES_MODEL(KFROM)%NLES_UND_PTS_ll=>NLES_UND_PTS_ll
LES_MODEL(KFROM)%XLES_MEAN_U=>XLES_MEAN_U
LES_MODEL(KFROM)%XLES_MEAN_V=>XLES_MEAN_V
LES_MODEL(KFROM)%XLES_MEAN_W=>XLES_MEAN_W
LES_MODEL(KFROM)%XLES_MEAN_P=>XLES_MEAN_P
LES_MODEL(KFROM)%XLES_MEAN_RHO=>XLES_MEAN_RHO
LES_MODEL(KFROM)%XLES_MEAN_Th=>XLES_MEAN_Th
LES_MODEL(KFROM)%XLES_MEAN_Thv=>XLES_MEAN_Thv
LES_MODEL(KFROM)%XLES_MEAN_Thl=>XLES_MEAN_Thl
LES_MODEL(KFROM)%XLES_MEAN_Rt=>XLES_MEAN_Rt
LES_MODEL(KFROM)%XLES_MEAN_Rv=>XLES_MEAN_Rv
LES_MODEL(KFROM)%XLES_MEAN_Rehu=>XLES_MEAN_Rehu
LES_MODEL(KFROM)%XLES_MEAN_Qs=>XLES_MEAN_Qs
LES_MODEL(KFROM)%XLES_MEAN_Rc=>XLES_MEAN_Rc
LES_MODEL(KFROM)%XLES_MEAN_Cf=>XLES_MEAN_Cf
LES_MODEL(KFROM)%XLES_MEAN_INDCf=>XLES_MEAN_INDCf
LES_MODEL(KFROM)%XLES_MEAN_INDCf2=>XLES_MEAN_INDCf2
LES_MODEL(KFROM)%XLES_MEAN_RF=>XLES_MEAN_RF
LES_MODEL(KFROM)%XLES_MEAN_Mf=>XLES_MEAN_Mf
LES_MODEL(KFROM)%XLES_MEAN_KHt=>XLES_MEAN_KHt
LES_MODEL(KFROM)%XLES_MEAN_KHr=>XLES_MEAN_KHr
LES_MODEL(KFROM)%XLES_MEAN_Rr=>XLES_MEAN_Rr
LES_MODEL(KFROM)%XLES_MEAN_Ri=>XLES_MEAN_Ri
LES_MODEL(KFROM)%XLES_MEAN_Rs=>XLES_MEAN_Rs
LES_MODEL(KFROM)%XLES_MEAN_Rg=>XLES_MEAN_Rg
LES_MODEL(KFROM)%XLES_MEAN_Rh=>XLES_MEAN_Rh
LES_MODEL(KFROM)%XLES_MEAN_Sv=>XLES_MEAN_Sv
LES_MODEL(KFROM)%XLES_MEAN_WIND=>XLES_MEAN_WIND
LES_MODEL(KFROM)%XLES_MEAN_dUdz=>XLES_MEAN_dUdz
LES_MODEL(KFROM)%XLES_MEAN_dVdz=>XLES_MEAN_dVdz
LES_MODEL(KFROM)%XLES_MEAN_dWdz=>XLES_MEAN_dWdz
LES_MODEL(KFROM)%XLES_MEAN_dThldz=>XLES_MEAN_dThldz
LES_MODEL(KFROM)%XLES_MEAN_dRtdz=>XLES_MEAN_dRtdz
LES_MODEL(KFROM)%XLES_MEAN_dSvdz=>XLES_MEAN_dSvdz
LES_MODEL(KFROM)%XLES_MEAN_DP=>XLES_MEAN_DP
LES_MODEL(KFROM)%XLES_MEAN_TP=>XLES_MEAN_TP
LES_MODEL(KFROM)%XLES_MEAN_TR=>XLES_MEAN_TR
LES_MODEL(KFROM)%XLES_MEAN_DISS=>XLES_MEAN_DISS
LES_MODEL(KFROM)%XLES_MEAN_LM=>XLES_MEAN_LM
LES_MODEL(KFROM)%XLES_RESOLVED_U2=>XLES_RESOLVED_U2
LES_MODEL(KFROM)%XLES_RESOLVED_V2=>XLES_RESOLVED_V2
LES_MODEL(KFROM)%XLES_RESOLVED_W2=>XLES_RESOLVED_W2
LES_MODEL(KFROM)%XLES_RESOLVED_P2=>XLES_RESOLVED_P2
LES_MODEL(KFROM)%XLES_RESOLVED_Th2=>XLES_RESOLVED_Th2
LES_MODEL(KFROM)%XLES_RESOLVED_Thl2=>XLES_RESOLVED_Thl2
LES_MODEL(KFROM)%XLES_RESOLVED_ThThv=>XLES_RESOLVED_ThThv
LES_MODEL(KFROM)%XLES_RESOLVED_ThlThv=>XLES_RESOLVED_ThlThv
LES_MODEL(KFROM)%XLES_RESOLVED_Ke=>XLES_RESOLVED_Ke
LES_MODEL(KFROM)%XLES_RESOLVED_Rv2=>XLES_RESOLVED_Rv2
LES_MODEL(KFROM)%XLES_RESOLVED_Rc2=>XLES_RESOLVED_Rc2
LES_MODEL(KFROM)%XLES_RESOLVED_Ri2=>XLES_RESOLVED_Ri2
LES_MODEL(KFROM)%XLES_RESOLVED_Rt2=>XLES_RESOLVED_Rt2
LES_MODEL(KFROM)%XLES_RESOLVED_Sv2=>XLES_RESOLVED_Sv2
LES_MODEL(KFROM)%XLES_RESOLVED_ThRv=>XLES_RESOLVED_ThRv
LES_MODEL(KFROM)%XLES_RESOLVED_ThRc=>XLES_RESOLVED_ThRc
LES_MODEL(KFROM)%XLES_RESOLVED_ThRi=>XLES_RESOLVED_ThRi
LES_MODEL(KFROM)%XLES_RESOLVED_ThSv=>XLES_RESOLVED_ThSv
LES_MODEL(KFROM)%XLES_RESOLVED_ThlRv=>XLES_RESOLVED_ThlRv
LES_MODEL(KFROM)%XLES_RESOLVED_ThlRc=>XLES_RESOLVED_ThlRc
LES_MODEL(KFROM)%XLES_RESOLVED_ThlRi=>XLES_RESOLVED_ThlRi
LES_MODEL(KFROM)%XLES_RESOLVED_ThlSv=>XLES_RESOLVED_ThlSv
LES_MODEL(KFROM)%XLES_RESOLVED_ThvRv=>XLES_RESOLVED_ThvRv
LES_MODEL(KFROM)%XLES_RESOLVED_ThvRc=>XLES_RESOLVED_ThvRc
LES_MODEL(KFROM)%XLES_RESOLVED_ThvRi=>XLES_RESOLVED_ThvRi
LES_MODEL(KFROM)%XLES_RESOLVED_ThvSv=>XLES_RESOLVED_ThvSv
LES_MODEL(KFROM)%XLES_RESOLVED_UV=>XLES_RESOLVED_UV
LES_MODEL(KFROM)%XLES_RESOLVED_WU=>XLES_RESOLVED_WU
LES_MODEL(KFROM)%XLES_RESOLVED_WV=>XLES_RESOLVED_WV
LES_MODEL(KFROM)%XLES_RESOLVED_UP=>XLES_RESOLVED_UP
LES_MODEL(KFROM)%XLES_RESOLVED_VP=>XLES_RESOLVED_VP
LES_MODEL(KFROM)%XLES_RESOLVED_WP=>XLES_RESOLVED_WP
LES_MODEL(KFROM)%XLES_RESOLVED_UTh=>XLES_RESOLVED_UTh
LES_MODEL(KFROM)%XLES_RESOLVED_VTh=>XLES_RESOLVED_VTh
LES_MODEL(KFROM)%XLES_RESOLVED_WTh=>XLES_RESOLVED_WTh
LES_MODEL(KFROM)%XLES_RESOLVED_UThl=>XLES_RESOLVED_UThl
LES_MODEL(KFROM)%XLES_RESOLVED_VThl=>XLES_RESOLVED_VThl
LES_MODEL(KFROM)%XLES_RESOLVED_WThl=>XLES_RESOLVED_WThl
LES_MODEL(KFROM)%XLES_RESOLVED_UThv=>XLES_RESOLVED_UThv
LES_MODEL(KFROM)%XLES_RESOLVED_VThv=>XLES_RESOLVED_VThv
LES_MODEL(KFROM)%XLES_RESOLVED_WThv=>XLES_RESOLVED_WThv
LES_MODEL(KFROM)%XLES_RESOLVED_URv=>XLES_RESOLVED_URv
LES_MODEL(KFROM)%XLES_RESOLVED_VRv=>XLES_RESOLVED_VRv
LES_MODEL(KFROM)%XLES_RESOLVED_WRv=>XLES_RESOLVED_WRv
LES_MODEL(KFROM)%XLES_RESOLVED_URc=>XLES_RESOLVED_URc
LES_MODEL(KFROM)%XLES_RESOLVED_VRc=>XLES_RESOLVED_VRc
LES_MODEL(KFROM)%XLES_RESOLVED_WRc=>XLES_RESOLVED_WRc
LES_MODEL(KFROM)%XLES_RESOLVED_URi=>XLES_RESOLVED_URi
LES_MODEL(KFROM)%XLES_RESOLVED_VRi=>XLES_RESOLVED_VRi
LES_MODEL(KFROM)%XLES_RESOLVED_WRi=>XLES_RESOLVED_WRi
LES_MODEL(KFROM)%XLES_RESOLVED_WRt=>XLES_RESOLVED_WRt
LES_MODEL(KFROM)%XLES_RESOLVED_WRr=>XLES_RESOLVED_WRr
LES_MODEL(KFROM)%XLES_INPRR3D=>XLES_INPRR3D     
LES_MODEL(KFROM)%XLES_MAX_INPRR3D=>XLES_MAX_INPRR3D     
LES_MODEL(KFROM)%XLES_EVAP3D=>XLES_EVAP3D     
LES_MODEL(KFROM)%XLES_RESOLVED_USv=>XLES_RESOLVED_USv
LES_MODEL(KFROM)%XLES_RESOLVED_VSv=>XLES_RESOLVED_VSv
LES_MODEL(KFROM)%XLES_RESOLVED_WSv=>XLES_RESOLVED_WSv
LES_MODEL(KFROM)%XLES_RESOLVED_U3=>XLES_RESOLVED_U3
LES_MODEL(KFROM)%XLES_RESOLVED_V3=>XLES_RESOLVED_V3
LES_MODEL(KFROM)%XLES_RESOLVED_W3=>XLES_RESOLVED_W3
LES_MODEL(KFROM)%XLES_RESOLVED_U4=>XLES_RESOLVED_U4
LES_MODEL(KFROM)%XLES_RESOLVED_V4=>XLES_RESOLVED_V4
LES_MODEL(KFROM)%XLES_RESOLVED_W4=>XLES_RESOLVED_W4
LES_MODEL(KFROM)%XLES_RESOLVED_Ua_ddxa_P=>XLES_RESOLVED_Ua_ddxa_P
LES_MODEL(KFROM)%XLES_RESOLVED_ThlPz=>XLES_RESOLVED_ThlPz
LES_MODEL(KFROM)%XLES_RESOLVED_WThl2=>XLES_RESOLVED_WThl2
LES_MODEL(KFROM)%XLES_RESOLVED_W2Thl=>XLES_RESOLVED_W2Thl
LES_MODEL(KFROM)%XLES_RESOLVED_WRv2=>XLES_RESOLVED_WRv2
LES_MODEL(KFROM)%XLES_RESOLVED_W2Rv=>XLES_RESOLVED_W2Rv
LES_MODEL(KFROM)%XLES_RESOLVED_WRt2=>XLES_RESOLVED_WRt2
LES_MODEL(KFROM)%XLES_RESOLVED_W2Rt=>XLES_RESOLVED_W2Rt
LES_MODEL(KFROM)%XLES_RESOLVED_RvPz=>XLES_RESOLVED_RvPz
LES_MODEL(KFROM)%XLES_RESOLVED_RtPz=>XLES_RESOLVED_RtPz
LES_MODEL(KFROM)%XLES_RESOLVED_WThlRv=>XLES_RESOLVED_WThlRv
LES_MODEL(KFROM)%XLES_RESOLVED_WThlRt=>XLES_RESOLVED_WThlRt
LES_MODEL(KFROM)%XLES_RESOLVED_WRc2=>XLES_RESOLVED_WRc2
LES_MODEL(KFROM)%XLES_RESOLVED_W2Rc=>XLES_RESOLVED_W2Rc
LES_MODEL(KFROM)%XLES_RESOLVED_RCPz=>XLES_RESOLVED_RCPz
LES_MODEL(KFROM)%XLES_RESOLVED_WThlRc=>XLES_RESOLVED_WThlRc
LES_MODEL(KFROM)%XLES_RESOLVED_WRvRc=>XLES_RESOLVED_WRvRc
LES_MODEL(KFROM)%XLES_RESOLVED_WRi2=>XLES_RESOLVED_WRi2
LES_MODEL(KFROM)%XLES_RESOLVED_W2Ri=>XLES_RESOLVED_W2Ri
LES_MODEL(KFROM)%XLES_RESOLVED_RIPz=>XLES_RESOLVED_RIPz
LES_MODEL(KFROM)%XLES_RESOLVED_WThlRi=>XLES_RESOLVED_WThlRi
LES_MODEL(KFROM)%XLES_RESOLVED_WRvRi=>XLES_RESOLVED_WRvRi
LES_MODEL(KFROM)%XLES_RESOLVED_WSv2=>XLES_RESOLVED_WSv2
LES_MODEL(KFROM)%XLES_RESOLVED_W2Sv=>XLES_RESOLVED_W2Sv
LES_MODEL(KFROM)%XLES_RESOLVED_SvPz=>XLES_RESOLVED_SvPz
LES_MODEL(KFROM)%XLES_RESOLVED_WThlSv=>XLES_RESOLVED_WThlSv
LES_MODEL(KFROM)%XLES_RESOLVED_WRvSv=>XLES_RESOLVED_WRvSv
LES_MODEL(KFROM)%XLES_RESOLVED_MASSFX=>XLES_RESOLVED_MASSFX
LES_MODEL(KFROM)%XLES_RESOLVED_UKe=>XLES_RESOLVED_UKe
LES_MODEL(KFROM)%XLES_RESOLVED_VKe=>XLES_RESOLVED_VKe
LES_MODEL(KFROM)%XLES_RESOLVED_WKe=>XLES_RESOLVED_WKe
LES_MODEL(KFROM)%XLES_BU_RES_KE=>XLES_BU_RES_KE
LES_MODEL(KFROM)%XLES_BU_RES_WThl=>XLES_BU_RES_WThl
LES_MODEL(KFROM)%XLES_BU_RES_Thl2=>XLES_BU_RES_Thl2
LES_MODEL(KFROM)%XLES_BU_RES_WRt=>XLES_BU_RES_WRt
LES_MODEL(KFROM)%XLES_BU_RES_Rt2=>XLES_BU_RES_Rt2
LES_MODEL(KFROM)%XLES_BU_RES_ThlRt=>XLES_BU_RES_ThlRt
LES_MODEL(KFROM)%XLES_BU_RES_Sv2=>XLES_BU_RES_Sv2
LES_MODEL(KFROM)%XLES_BU_RES_WSv=>XLES_BU_RES_WSv
LES_MODEL(KFROM)%XLES_RES_U_SBG_Tke=>XLES_RES_U_SBG_Tke
LES_MODEL(KFROM)%XLES_RES_V_SBG_Tke=>XLES_RES_V_SBG_Tke
LES_MODEL(KFROM)%XLES_RES_W_SBG_Tke=>XLES_RES_W_SBG_Tke
LES_MODEL(KFROM)%XLES_RES_W_SBG_WThl=>XLES_RES_W_SBG_WThl
LES_MODEL(KFROM)%XLES_RES_W_SBG_WRt=>XLES_RES_W_SBG_WRt
LES_MODEL(KFROM)%XLES_RES_W_SBG_Thl2=>XLES_RES_W_SBG_Thl2
LES_MODEL(KFROM)%XLES_RES_W_SBG_Rt2=>XLES_RES_W_SBG_Rt2
LES_MODEL(KFROM)%XLES_RES_W_SBG_ThlRt=>XLES_RES_W_SBG_ThlRt
LES_MODEL(KFROM)%XLES_RES_W_SBG_WSv=>XLES_RES_W_SBG_WSv
LES_MODEL(KFROM)%XLES_RES_W_SBG_Sv2=>XLES_RES_W_SBG_Sv2
LES_MODEL(KFROM)%XLES_RES_ddxa_U_SBG_UaU=>XLES_RES_ddxa_U_SBG_UaU
LES_MODEL(KFROM)%XLES_RES_ddxa_V_SBG_UaV=>XLES_RES_ddxa_V_SBG_UaV
LES_MODEL(KFROM)%XLES_RES_ddxa_W_SBG_UaW=>XLES_RES_ddxa_W_SBG_UaW
LES_MODEL(KFROM)%XLES_RES_ddxa_W_SBG_UaThl=>XLES_RES_ddxa_W_SBG_UaThl
LES_MODEL(KFROM)%XLES_RES_ddxa_Thl_SBG_UaW=>XLES_RES_ddxa_Thl_SBG_UaW
LES_MODEL(KFROM)%XLES_RES_ddz_Thl_SBG_W2=>XLES_RES_ddz_Thl_SBG_W2
LES_MODEL(KFROM)%XLES_RES_ddxa_W_SBG_UaRt=>XLES_RES_ddxa_W_SBG_UaRt
LES_MODEL(KFROM)%XLES_RES_ddxa_Rt_SBG_UaW=>XLES_RES_ddxa_Rt_SBG_UaW
LES_MODEL(KFROM)%XLES_RES_ddz_Rt_SBG_W2=>XLES_RES_ddz_Rt_SBG_W2
LES_MODEL(KFROM)%XLES_RES_ddxa_Thl_SBG_UaRt=>XLES_RES_ddxa_Thl_SBG_UaRt
LES_MODEL(KFROM)%XLES_RES_ddxa_Rt_SBG_UaThl=>XLES_RES_ddxa_Rt_SBG_UaThl
LES_MODEL(KFROM)%XLES_RES_ddxa_Thl_SBG_UaThl=>XLES_RES_ddxa_Thl_SBG_UaThl
LES_MODEL(KFROM)%XLES_RES_ddxa_Rt_SBG_UaRt=>XLES_RES_ddxa_Rt_SBG_UaRt
LES_MODEL(KFROM)%XLES_RES_ddxa_W_SBG_UaSv=>XLES_RES_ddxa_W_SBG_UaSv
LES_MODEL(KFROM)%XLES_RES_ddxa_Sv_SBG_UaW=>XLES_RES_ddxa_Sv_SBG_UaW
LES_MODEL(KFROM)%XLES_RES_ddz_Sv_SBG_W2=>XLES_RES_ddz_Sv_SBG_W2
LES_MODEL(KFROM)%XLES_RES_ddxa_Sv_SBG_UaSv=>XLES_RES_ddxa_Sv_SBG_UaSv
LES_MODEL(KFROM)%XLES_BU_SBG_TKE=>XLES_BU_SBG_TKE
LES_MODEL(KFROM)%XLES_SUBGRID_U2=>XLES_SUBGRID_U2
LES_MODEL(KFROM)%XLES_SUBGRID_V2=>XLES_SUBGRID_V2
LES_MODEL(KFROM)%XLES_SUBGRID_W2=>XLES_SUBGRID_W2
LES_MODEL(KFROM)%XLES_SUBGRID_Tke=>XLES_SUBGRID_Tke
LES_MODEL(KFROM)%XLES_SUBGRID_Thl2=>XLES_SUBGRID_Thl2
LES_MODEL(KFROM)%XLES_SUBGRID_Rt2=>XLES_SUBGRID_Rt2
LES_MODEL(KFROM)%XLES_SUBGRID_Rc2=>XLES_SUBGRID_Rc2
LES_MODEL(KFROM)%XLES_SUBGRID_Ri2=>XLES_SUBGRID_Ri2
LES_MODEL(KFROM)%XLES_SUBGRID_ThlRt=>XLES_SUBGRID_ThlRt
LES_MODEL(KFROM)%XLES_SUBGRID_Sv2=>XLES_SUBGRID_Sv2
LES_MODEL(KFROM)%XLES_SUBGRID_UV=>XLES_SUBGRID_UV
LES_MODEL(KFROM)%XLES_SUBGRID_WU=>XLES_SUBGRID_WU
LES_MODEL(KFROM)%XLES_SUBGRID_WV=>XLES_SUBGRID_WV
LES_MODEL(KFROM)%XLES_SUBGRID_UThl=>XLES_SUBGRID_UThl
LES_MODEL(KFROM)%XLES_SUBGRID_VThl=>XLES_SUBGRID_VThl
LES_MODEL(KFROM)%XLES_SUBGRID_WThl=>XLES_SUBGRID_WThl
LES_MODEL(KFROM)%XLES_SUBGRID_ENTLES=>XLES_SUBGRID_ENTLES
LES_MODEL(KFROM)%XLES_SUBGRID_DETLES=>XLES_SUBGRID_DETLES
LES_MODEL(KFROM)%XLES_SUBGRID_URt=>XLES_SUBGRID_URt
LES_MODEL(KFROM)%XLES_SUBGRID_VRt=>XLES_SUBGRID_VRt
LES_MODEL(KFROM)%XLES_SUBGRID_WRt=>XLES_SUBGRID_WRt
LES_MODEL(KFROM)%XLES_SUBGRID_URc=>XLES_SUBGRID_URc
LES_MODEL(KFROM)%XLES_SUBGRID_VRc=>XLES_SUBGRID_VRc
LES_MODEL(KFROM)%XLES_SUBGRID_WRc=>XLES_SUBGRID_WRc
LES_MODEL(KFROM)%XLES_SUBGRID_USv=>XLES_SUBGRID_USv
LES_MODEL(KFROM)%XLES_SUBGRID_VSv=>XLES_SUBGRID_VSv
LES_MODEL(KFROM)%XLES_SUBGRID_WSv=>XLES_SUBGRID_WSv
LES_MODEL(KFROM)%XLES_SUBGRID_UTke=>XLES_SUBGRID_UTke
LES_MODEL(KFROM)%XLES_SUBGRID_VTke=>XLES_SUBGRID_VTke
LES_MODEL(KFROM)%XLES_SUBGRID_WTke=>XLES_SUBGRID_WTke
LES_MODEL(KFROM)%XLES_SUBGRID_ddz_WTke=>XLES_SUBGRID_ddz_WTke
LES_MODEL(KFROM)%XLES_SUBGRID_WThv=>XLES_SUBGRID_WThv
LES_MODEL(KFROM)%XLES_SUBGRID_ThlThv=>XLES_SUBGRID_ThlThv
LES_MODEL(KFROM)%XLES_SUBGRID_RtThv=>XLES_SUBGRID_RtThv
LES_MODEL(KFROM)%XLES_SUBGRID_SvThv=>XLES_SUBGRID_SvThv
LES_MODEL(KFROM)%XLES_SUBGRID_W2Thl=>XLES_SUBGRID_W2Thl
LES_MODEL(KFROM)%XLES_SUBGRID_W2Rt=>XLES_SUBGRID_W2Rt
LES_MODEL(KFROM)%XLES_SUBGRID_WThlRt=>XLES_SUBGRID_WThlRt
LES_MODEL(KFROM)%XLES_SUBGRID_W2Sv=>XLES_SUBGRID_W2Sv
LES_MODEL(KFROM)%XLES_SUBGRID_WThl2=>XLES_SUBGRID_WThl2
LES_MODEL(KFROM)%XLES_SUBGRID_WRt2=>XLES_SUBGRID_WRt2
LES_MODEL(KFROM)%XLES_SUBGRID_WSv2=>XLES_SUBGRID_WSv2
LES_MODEL(KFROM)%XLES_SUBGRID_DISS_Tke=>XLES_SUBGRID_DISS_Tke
LES_MODEL(KFROM)%XLES_SUBGRID_DISS_Thl2=>XLES_SUBGRID_DISS_Thl2
LES_MODEL(KFROM)%XLES_SUBGRID_DISS_Rt2=>XLES_SUBGRID_DISS_Rt2
LES_MODEL(KFROM)%XLES_SUBGRID_DISS_ThlRt=>XLES_SUBGRID_DISS_ThlRt
LES_MODEL(KFROM)%XLES_SUBGRID_DISS_Sv2=>XLES_SUBGRID_DISS_Sv2
LES_MODEL(KFROM)%XLES_SUBGRID_WP=>XLES_SUBGRID_WP
LES_MODEL(KFROM)%XLES_SUBGRID_ThlPz=>XLES_SUBGRID_ThlPz
LES_MODEL(KFROM)%XLES_SUBGRID_RtPz=>XLES_SUBGRID_RtPz
LES_MODEL(KFROM)%XLES_SUBGRID_SvPz=>XLES_SUBGRID_SvPz
LES_MODEL(KFROM)%XLES_SUBGRID_PHI3=>XLES_SUBGRID_PHI3
LES_MODEL(KFROM)%XLES_SUBGRID_PSI3=>XLES_SUBGRID_PSI3
LES_MODEL(KFROM)%XLES_SUBGRID_LMix=>XLES_SUBGRID_LMix
LES_MODEL(KFROM)%XLES_SUBGRID_LDiss=>XLES_SUBGRID_LDiss
LES_MODEL(KFROM)%XLES_SUBGRID_Km=>XLES_SUBGRID_Km
LES_MODEL(KFROM)%XLES_SUBGRID_Kh=>XLES_SUBGRID_Kh

LES_MODEL(KFROM)%XLES_SUBGRID_THLUP_MF=>XLES_SUBGRID_THLUP_MF
LES_MODEL(KFROM)%XLES_SUBGRID_RTUP_MF=>XLES_SUBGRID_RTUP_MF
LES_MODEL(KFROM)%XLES_SUBGRID_RVUP_MF=>XLES_SUBGRID_RVUP_MF
LES_MODEL(KFROM)%XLES_SUBGRID_RCUP_MF=>XLES_SUBGRID_RCUP_MF
LES_MODEL(KFROM)%XLES_SUBGRID_RIUP_MF=>XLES_SUBGRID_RIUP_MF
LES_MODEL(KFROM)%XLES_SUBGRID_WUP_MF=>XLES_SUBGRID_WUP_MF
LES_MODEL(KFROM)%XLES_SUBGRID_MASSFLUX=>XLES_SUBGRID_MASSFLUX
LES_MODEL(KFROM)%XLES_SUBGRID_DETR=>XLES_SUBGRID_DETR
LES_MODEL(KFROM)%XLES_SUBGRID_ENTR=>XLES_SUBGRID_ENTR
LES_MODEL(KFROM)%XLES_SUBGRID_FRACUP=>XLES_SUBGRID_FRACUP
LES_MODEL(KFROM)%XLES_SUBGRID_THVUP_MF=>XLES_SUBGRID_THVUP_MF
LES_MODEL(KFROM)%XLES_SUBGRID_WTHLMF=>XLES_SUBGRID_WTHLMF
LES_MODEL(KFROM)%XLES_SUBGRID_WRTMF=>XLES_SUBGRID_WRTMF
LES_MODEL(KFROM)%XLES_SUBGRID_WTHVMF=>XLES_SUBGRID_WTHVMF
LES_MODEL(KFROM)%XLES_SUBGRID_WUMF=>XLES_SUBGRID_WUMF
LES_MODEL(KFROM)%XLES_SUBGRID_WVMF=>XLES_SUBGRID_WVMF

LES_MODEL(KFROM)%XLES_UPDRAFT=>XLES_UPDRAFT
LES_MODEL(KFROM)%XLES_UPDRAFT_W=>XLES_UPDRAFT_W
LES_MODEL(KFROM)%XLES_UPDRAFT_Th=>XLES_UPDRAFT_Th
LES_MODEL(KFROM)%XLES_UPDRAFT_Thv=>XLES_UPDRAFT_Thv
LES_MODEL(KFROM)%XLES_UPDRAFT_Thl=>XLES_UPDRAFT_Thl
LES_MODEL(KFROM)%XLES_UPDRAFT_Rv=>XLES_UPDRAFT_Rv
LES_MODEL(KFROM)%XLES_UPDRAFT_Rc=>XLES_UPDRAFT_Rc
LES_MODEL(KFROM)%XLES_UPDRAFT_Rr=>XLES_UPDRAFT_Rr
LES_MODEL(KFROM)%XLES_UPDRAFT_Ri=>XLES_UPDRAFT_Ri
LES_MODEL(KFROM)%XLES_UPDRAFT_Rs=>XLES_UPDRAFT_Rs
LES_MODEL(KFROM)%XLES_UPDRAFT_Rg=>XLES_UPDRAFT_Rg
LES_MODEL(KFROM)%XLES_UPDRAFT_Rh=>XLES_UPDRAFT_Rh
LES_MODEL(KFROM)%XLES_UPDRAFT_Sv=>XLES_UPDRAFT_Sv
LES_MODEL(KFROM)%XLES_UPDRAFT_Tke=>XLES_UPDRAFT_Tke
LES_MODEL(KFROM)%XLES_UPDRAFT_Ke=>XLES_UPDRAFT_Ke
LES_MODEL(KFROM)%XLES_UPDRAFT_WTh=>XLES_UPDRAFT_WTh
LES_MODEL(KFROM)%XLES_UPDRAFT_WThl=>XLES_UPDRAFT_WThl
LES_MODEL(KFROM)%XLES_UPDRAFT_WThv=>XLES_UPDRAFT_WThv
LES_MODEL(KFROM)%XLES_UPDRAFT_WRv=>XLES_UPDRAFT_WRv
LES_MODEL(KFROM)%XLES_UPDRAFT_WRc=>XLES_UPDRAFT_WRc
LES_MODEL(KFROM)%XLES_UPDRAFT_WRi=>XLES_UPDRAFT_WRi
LES_MODEL(KFROM)%XLES_UPDRAFT_WSv=>XLES_UPDRAFT_WSv
LES_MODEL(KFROM)%XLES_UPDRAFT_Th2=>XLES_UPDRAFT_Th2
LES_MODEL(KFROM)%XLES_UPDRAFT_Thl2=>XLES_UPDRAFT_Thl2
LES_MODEL(KFROM)%XLES_UPDRAFT_ThThv=>XLES_UPDRAFT_ThThv
LES_MODEL(KFROM)%XLES_UPDRAFT_ThlThv=>XLES_UPDRAFT_ThlThv
LES_MODEL(KFROM)%XLES_UPDRAFT_Rv2=>XLES_UPDRAFT_Rv2
LES_MODEL(KFROM)%XLES_UPDRAFT_Rc2=>XLES_UPDRAFT_Rc2
LES_MODEL(KFROM)%XLES_UPDRAFT_Ri2=>XLES_UPDRAFT_Ri2
LES_MODEL(KFROM)%XLES_UPDRAFT_Sv2=>XLES_UPDRAFT_Sv2
LES_MODEL(KFROM)%XLES_UPDRAFT_ThRv=>XLES_UPDRAFT_ThRv
LES_MODEL(KFROM)%XLES_UPDRAFT_ThRc=>XLES_UPDRAFT_ThRc
LES_MODEL(KFROM)%XLES_UPDRAFT_ThRi=>XLES_UPDRAFT_ThRi
LES_MODEL(KFROM)%XLES_UPDRAFT_ThSv=>XLES_UPDRAFT_ThSv
LES_MODEL(KFROM)%XLES_UPDRAFT_ThvRv=>XLES_UPDRAFT_ThvRv
LES_MODEL(KFROM)%XLES_UPDRAFT_ThvRc=>XLES_UPDRAFT_ThvRc
LES_MODEL(KFROM)%XLES_UPDRAFT_ThvRi=>XLES_UPDRAFT_ThvRi
LES_MODEL(KFROM)%XLES_UPDRAFT_ThvSv=>XLES_UPDRAFT_ThvSv
LES_MODEL(KFROM)%XLES_UPDRAFT_ThlRv=>XLES_UPDRAFT_ThlRv
LES_MODEL(KFROM)%XLES_UPDRAFT_ThlRc=>XLES_UPDRAFT_ThlRc
LES_MODEL(KFROM)%XLES_UPDRAFT_ThlRi=>XLES_UPDRAFT_ThlRi
LES_MODEL(KFROM)%XLES_UPDRAFT_ThlSv=>XLES_UPDRAFT_ThlSv
LES_MODEL(KFROM)%XLES_DOWNDRAFT=>XLES_DOWNDRAFT
LES_MODEL(KFROM)%XLES_DOWNDRAFT_W=>XLES_DOWNDRAFT_W
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Th=>XLES_DOWNDRAFT_Th
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Thv=>XLES_DOWNDRAFT_Thv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Thl=>XLES_DOWNDRAFT_Thl
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rv=>XLES_DOWNDRAFT_Rv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rc=>XLES_DOWNDRAFT_Rc
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rr=>XLES_DOWNDRAFT_Rr
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Ri=>XLES_DOWNDRAFT_Ri
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rs=>XLES_DOWNDRAFT_Rs
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rg=>XLES_DOWNDRAFT_Rg
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rh=>XLES_DOWNDRAFT_Rh
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Sv=>XLES_DOWNDRAFT_Sv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Tke=>XLES_DOWNDRAFT_Tke
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Ke=>XLES_DOWNDRAFT_Ke
LES_MODEL(KFROM)%XLES_DOWNDRAFT_WTh=>XLES_DOWNDRAFT_WTh
LES_MODEL(KFROM)%XLES_DOWNDRAFT_WThl=>XLES_DOWNDRAFT_WThl
LES_MODEL(KFROM)%XLES_DOWNDRAFT_WThv=>XLES_DOWNDRAFT_WThv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_WRv=>XLES_DOWNDRAFT_WRv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_WRc=>XLES_DOWNDRAFT_WRc
LES_MODEL(KFROM)%XLES_DOWNDRAFT_WRi=>XLES_DOWNDRAFT_WRi
LES_MODEL(KFROM)%XLES_DOWNDRAFT_WSv=>XLES_DOWNDRAFT_WSv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Th2=>XLES_DOWNDRAFT_Th2
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Thl2=>XLES_DOWNDRAFT_Thl2
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThThv=>XLES_DOWNDRAFT_ThThv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThlThv=>XLES_DOWNDRAFT_ThlThv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rv2=>XLES_DOWNDRAFT_Rv2
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Rc2=>XLES_DOWNDRAFT_Rc2
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Ri2=>XLES_DOWNDRAFT_Ri2
LES_MODEL(KFROM)%XLES_DOWNDRAFT_Sv2=>XLES_DOWNDRAFT_Sv2
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThRv=>XLES_DOWNDRAFT_ThRv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThRc=>XLES_DOWNDRAFT_ThRc
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThRi=>XLES_DOWNDRAFT_ThRi
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThSv=>XLES_DOWNDRAFT_ThSv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThvRv=>XLES_DOWNDRAFT_ThvRv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThvRc=>XLES_DOWNDRAFT_ThvRc
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThvRi=>XLES_DOWNDRAFT_ThvRi
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThvSv=>XLES_DOWNDRAFT_ThvSv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThlRv=>XLES_DOWNDRAFT_ThlRv
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThlRc=>XLES_DOWNDRAFT_ThlRc
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThlRi=>XLES_DOWNDRAFT_ThlRi
LES_MODEL(KFROM)%XLES_DOWNDRAFT_ThlSv=>XLES_DOWNDRAFT_ThlSv
LES_MODEL(KFROM)%XLES_UW0=>XLES_UW0
LES_MODEL(KFROM)%XLES_VW0=>XLES_VW0
LES_MODEL(KFROM)%XLES_USTAR=>XLES_USTAR
LES_MODEL(KFROM)%XLES_WSTAR=>XLES_WSTAR
LES_MODEL(KFROM)%XLES_Q0=>XLES_Q0
LES_MODEL(KFROM)%XLES_E0=>XLES_E0
LES_MODEL(KFROM)%XLES_SV0=>XLES_SV0
LES_MODEL(KFROM)%XLES_BL_HEIGHT=>XLES_BL_HEIGHT
LES_MODEL(KFROM)%XLES_MO_LENGTH=>XLES_MO_LENGTH
LES_MODEL(KFROM)%XLES_ZCB=>XLES_ZCB
LES_MODEL(KFROM)%XLES_CFtot=>XLES_CFtot
LES_MODEL(KFROM)%XLES_CF2tot=>XLES_CF2tot
LES_MODEL(KFROM)%XLES_LWP=>XLES_LWP
LES_MODEL(KFROM)%XLES_LWPVAR=>XLES_LWPVAR
LES_MODEL(KFROM)%XLES_RWP=>XLES_RWP
LES_MODEL(KFROM)%XLES_IWP=>XLES_IWP
LES_MODEL(KFROM)%XLES_SWP=>XLES_SWP
LES_MODEL(KFROM)%XLES_GWP=>XLES_GWP
LES_MODEL(KFROM)%XLES_HWP=>XLES_HWP
LES_MODEL(KFROM)%XLES_INT_TKE=>XLES_INT_TKE
LES_MODEL(KFROM)%XLES_INPRR=>XLES_INPRR
LES_MODEL(KFROM)%XLES_INPRC=>XLES_INPRC
LES_MODEL(KFROM)%XLES_INDEP=>XLES_INDEP
LES_MODEL(KFROM)%XLES_RAIN_INPRR=>XLES_RAIN_INPRR
LES_MODEL(KFROM)%XLES_ACPRR=>XLES_ACPRR
LES_MODEL(KFROM)%XLES_ZMAXCF=>XLES_ZMAXCF
LES_MODEL(KFROM)%XLES_ZMAXCF2=>XLES_ZMAXCF2
LES_MODEL(KFROM)%XLES_PRECFR=>XLES_PRECFR
LES_MODEL(KFROM)%XCORRi_UU=>XCORRi_UU
LES_MODEL(KFROM)%XCORRi_VV=>XCORRi_VV
LES_MODEL(KFROM)%XCORRi_UV=>XCORRi_UV
LES_MODEL(KFROM)%XCORRi_WU=>XCORRi_WU
LES_MODEL(KFROM)%XCORRi_WV=>XCORRi_WV
LES_MODEL(KFROM)%XCORRi_WW=>XCORRi_WW
LES_MODEL(KFROM)%XCORRi_WTh=>XCORRi_WTh
LES_MODEL(KFROM)%XCORRi_WThl=>XCORRi_WThl
LES_MODEL(KFROM)%XCORRi_WRv=>XCORRi_WRv
LES_MODEL(KFROM)%XCORRi_WRc=>XCORRi_WRc
LES_MODEL(KFROM)%XCORRi_WRi=>XCORRi_WRi
LES_MODEL(KFROM)%XCORRi_WSv=>XCORRi_WSv
LES_MODEL(KFROM)%XCORRi_ThTh=>XCORRi_ThTh
LES_MODEL(KFROM)%XCORRi_ThRv=>XCORRi_ThRv
LES_MODEL(KFROM)%XCORRi_ThRc=>XCORRi_ThRc
LES_MODEL(KFROM)%XCORRi_ThRi=>XCORRi_ThRi
LES_MODEL(KFROM)%XCORRi_ThlThl=>XCORRi_ThlThl
LES_MODEL(KFROM)%XCORRi_ThlRv=>XCORRi_ThlRv
LES_MODEL(KFROM)%XCORRi_ThlRc=>XCORRi_ThlRc
LES_MODEL(KFROM)%XCORRi_ThlRi=>XCORRi_ThlRi
LES_MODEL(KFROM)%XCORRi_RvRv=>XCORRi_RvRv
LES_MODEL(KFROM)%XCORRi_RcRc=>XCORRi_RcRc
LES_MODEL(KFROM)%XCORRi_RiRi=>XCORRi_RiRi
LES_MODEL(KFROM)%XCORRi_SvSv=>XCORRi_SvSv
LES_MODEL(KFROM)%XCORRj_UU=>XCORRj_UU
LES_MODEL(KFROM)%XCORRj_VV=>XCORRj_VV
LES_MODEL(KFROM)%XCORRj_UV=>XCORRj_UV
LES_MODEL(KFROM)%XCORRj_WU=>XCORRj_WU
LES_MODEL(KFROM)%XCORRj_WV=>XCORRj_WV
LES_MODEL(KFROM)%XCORRj_WW=>XCORRj_WW
LES_MODEL(KFROM)%XCORRj_WTh=>XCORRj_WTh
LES_MODEL(KFROM)%XCORRj_WThl=>XCORRj_WThl
LES_MODEL(KFROM)%XCORRj_WRv=>XCORRj_WRv
LES_MODEL(KFROM)%XCORRj_WRc=>XCORRj_WRc
LES_MODEL(KFROM)%XCORRj_WRi=>XCORRj_WRi
LES_MODEL(KFROM)%XCORRj_WSv=>XCORRj_WSv
LES_MODEL(KFROM)%XCORRj_ThTh=>XCORRj_ThTh
LES_MODEL(KFROM)%XCORRj_ThRv=>XCORRj_ThRv
LES_MODEL(KFROM)%XCORRj_ThRc=>XCORRj_ThRc
LES_MODEL(KFROM)%XCORRj_ThRi=>XCORRj_ThRi
LES_MODEL(KFROM)%XCORRj_ThlThl=>XCORRj_ThlThl
LES_MODEL(KFROM)%XCORRj_ThlRv=>XCORRj_ThlRv
LES_MODEL(KFROM)%XCORRj_ThlRc=>XCORRj_ThlRc
LES_MODEL(KFROM)%XCORRj_ThlRi=>XCORRj_ThlRi
LES_MODEL(KFROM)%XCORRj_RvRv=>XCORRj_RvRv
LES_MODEL(KFROM)%XCORRj_RcRc=>XCORRj_RcRc
LES_MODEL(KFROM)%XCORRj_RiRi=>XCORRj_RiRi
LES_MODEL(KFROM)%XCORRj_SvSv=>XCORRj_SvSv
LES_MODEL(KFROM)%XCORRk_WW=>XCORRk_WW
LES_MODEL(KFROM)%XLES_SWU=>XLES_SWU
LES_MODEL(KFROM)%XLES_SWD=>XLES_SWD
LES_MODEL(KFROM)%XLES_LWU=>XLES_LWU
LES_MODEL(KFROM)%XLES_LWD=>XLES_LWD
LES_MODEL(KFROM)%XLES_DTHRADSW=>XLES_DTHRADSW
LES_MODEL(KFROM)%XLES_DTHRADLW=>XLES_DTHRADLW
LES_MODEL(KFROM)%XLES_RADEFF=>XLES_RADEFF
!
! Current model is set to model KTO
NLES_TIMES=>LES_MODEL(KTO)%NLES_TIMES
NLES_DTCOUNT=>LES_MODEL(KTO)%NLES_DTCOUNT
NLES_TCOUNT=>LES_MODEL(KTO)%NLES_TCOUNT
NSPECTRA_NI=>LES_MODEL(KTO)%NSPECTRA_NI
NSPECTRA_NJ=>LES_MODEL(KTO)%NSPECTRA_NJ
XLES_DATIME=>LES_MODEL(KTO)%XLES_DATIME
XLES_TRAJT=>LES_MODEL(KTO)%XLES_TRAJT
XLES_Z=>LES_MODEL(KTO)%XLES_Z
XLES_ZS=>LES_MODEL(KTO)%XLES_ZS
XCOEFLIN_LES=>LES_MODEL(KTO)%XCOEFLIN_LES
NKLIN_LES=>LES_MODEL(KTO)%NKLIN_LES
XCOEFLIN_SPEC=>LES_MODEL(KTO)%XCOEFLIN_SPEC
NKLIN_SPEC=>LES_MODEL(KTO)%NKLIN_SPEC
NLES_AVG_PTS_ll=>LES_MODEL(KTO)%NLES_AVG_PTS_ll
NLES_UND_PTS_ll=>LES_MODEL(KTO)%NLES_UND_PTS_ll
XLES_MEAN_U=>LES_MODEL(KTO)%XLES_MEAN_U
XLES_MEAN_V=>LES_MODEL(KTO)%XLES_MEAN_V
XLES_MEAN_W=>LES_MODEL(KTO)%XLES_MEAN_W
XLES_MEAN_P=>LES_MODEL(KTO)%XLES_MEAN_P
XLES_MEAN_RHO=>LES_MODEL(KTO)%XLES_MEAN_RHO
XLES_MEAN_Th=>LES_MODEL(KTO)%XLES_MEAN_Th
XLES_MEAN_Thv=>LES_MODEL(KTO)%XLES_MEAN_Thv
XLES_MEAN_Thl=>LES_MODEL(KTO)%XLES_MEAN_Thl
XLES_MEAN_Rt=>LES_MODEL(KTO)%XLES_MEAN_Rt
XLES_MEAN_Rv=>LES_MODEL(KTO)%XLES_MEAN_Rv
XLES_MEAN_Rehu=>LES_MODEL(KTO)%XLES_MEAN_Rehu
XLES_MEAN_Qs=>LES_MODEL(KTO)%XLES_MEAN_Qs
XLES_MEAN_Rc=>LES_MODEL(KTO)%XLES_MEAN_Rc
XLES_MEAN_Cf=>LES_MODEL(KTO)%XLES_MEAN_Cf
XLES_MEAN_INDCf=>LES_MODEL(KTO)%XLES_MEAN_INDCf
XLES_MEAN_INDCf2=>LES_MODEL(KTO)%XLES_MEAN_INDCf2
XLES_MEAN_RF=>LES_MODEL(KTO)%XLES_MEAN_RF
XLES_MEAN_Mf=>LES_MODEL(KTO)%XLES_MEAN_Mf
XLES_MEAN_KHt=>LES_MODEL(KTO)%XLES_MEAN_KHt
XLES_MEAN_KHr=>LES_MODEL(KTO)%XLES_MEAN_KHr
XLES_MEAN_Rr=>LES_MODEL(KTO)%XLES_MEAN_Rr
XLES_MEAN_Ri=>LES_MODEL(KTO)%XLES_MEAN_Ri
XLES_MEAN_Rs=>LES_MODEL(KTO)%XLES_MEAN_Rs
XLES_MEAN_Rg=>LES_MODEL(KTO)%XLES_MEAN_Rg
XLES_MEAN_Rh=>LES_MODEL(KTO)%XLES_MEAN_Rh
XLES_MEAN_Sv=>LES_MODEL(KTO)%XLES_MEAN_Sv
XLES_MEAN_WIND=>LES_MODEL(KTO)%XLES_MEAN_WIND
XLES_MEAN_dUdz=>LES_MODEL(KTO)%XLES_MEAN_dUdz
XLES_MEAN_dVdz=>LES_MODEL(KTO)%XLES_MEAN_dVdz
XLES_MEAN_dWdz=>LES_MODEL(KTO)%XLES_MEAN_dWdz
XLES_MEAN_dThldz=>LES_MODEL(KTO)%XLES_MEAN_dThldz
XLES_MEAN_dRtdz=>LES_MODEL(KTO)%XLES_MEAN_dRtdz
XLES_MEAN_dSvdz=>LES_MODEL(KTO)%XLES_MEAN_dSvdz
XLES_MEAN_DP=>LES_MODEL(KTO)%XLES_MEAN_DP
XLES_MEAN_TP=>LES_MODEL(KTO)%XLES_MEAN_TP
XLES_MEAN_TR=>LES_MODEL(KTO)%XLES_MEAN_TR
XLES_MEAN_DISS=>LES_MODEL(KTO)%XLES_MEAN_DISS
XLES_MEAN_LM=>LES_MODEL(KTO)%XLES_MEAN_LM
XLES_RESOLVED_U2=>LES_MODEL(KTO)%XLES_RESOLVED_U2
XLES_RESOLVED_V2=>LES_MODEL(KTO)%XLES_RESOLVED_V2
XLES_RESOLVED_W2=>LES_MODEL(KTO)%XLES_RESOLVED_W2
XLES_RESOLVED_P2=>LES_MODEL(KTO)%XLES_RESOLVED_P2
XLES_RESOLVED_Th2=>LES_MODEL(KTO)%XLES_RESOLVED_Th2
XLES_RESOLVED_Thl2=>LES_MODEL(KTO)%XLES_RESOLVED_Thl2
XLES_RESOLVED_ThThv=>LES_MODEL(KTO)%XLES_RESOLVED_ThThv
XLES_RESOLVED_ThlThv=>LES_MODEL(KTO)%XLES_RESOLVED_ThlThv
XLES_RESOLVED_Ke=>LES_MODEL(KTO)%XLES_RESOLVED_Ke
XLES_RESOLVED_Rv2=>LES_MODEL(KTO)%XLES_RESOLVED_Rv2
XLES_RESOLVED_Rc2=>LES_MODEL(KTO)%XLES_RESOLVED_Rc2
XLES_RESOLVED_Ri2=>LES_MODEL(KTO)%XLES_RESOLVED_Ri2
XLES_RESOLVED_Rt2=>LES_MODEL(KTO)%XLES_RESOLVED_Rt2
XLES_RESOLVED_Sv2=>LES_MODEL(KTO)%XLES_RESOLVED_Sv2
XLES_RESOLVED_ThRv=>LES_MODEL(KTO)%XLES_RESOLVED_ThRv
XLES_RESOLVED_ThRc=>LES_MODEL(KTO)%XLES_RESOLVED_ThRc
XLES_RESOLVED_ThRi=>LES_MODEL(KTO)%XLES_RESOLVED_ThRi
XLES_RESOLVED_ThSv=>LES_MODEL(KTO)%XLES_RESOLVED_ThSv
XLES_RESOLVED_ThlRv=>LES_MODEL(KTO)%XLES_RESOLVED_ThlRv
XLES_RESOLVED_ThlRc=>LES_MODEL(KTO)%XLES_RESOLVED_ThlRc
XLES_RESOLVED_ThlRi=>LES_MODEL(KTO)%XLES_RESOLVED_ThlRi
XLES_RESOLVED_ThlSv=>LES_MODEL(KTO)%XLES_RESOLVED_ThlSv
XLES_RESOLVED_ThvRv=>LES_MODEL(KTO)%XLES_RESOLVED_ThvRv
XLES_RESOLVED_ThvRc=>LES_MODEL(KTO)%XLES_RESOLVED_ThvRc
XLES_RESOLVED_ThvRi=>LES_MODEL(KTO)%XLES_RESOLVED_ThvRi
XLES_RESOLVED_ThvSv=>LES_MODEL(KTO)%XLES_RESOLVED_ThvSv
XLES_RESOLVED_UV=>LES_MODEL(KTO)%XLES_RESOLVED_UV
XLES_RESOLVED_WU=>LES_MODEL(KTO)%XLES_RESOLVED_WU
XLES_RESOLVED_WV=>LES_MODEL(KTO)%XLES_RESOLVED_WV
XLES_RESOLVED_UP=>LES_MODEL(KTO)%XLES_RESOLVED_UP
XLES_RESOLVED_VP=>LES_MODEL(KTO)%XLES_RESOLVED_VP
XLES_RESOLVED_WP=>LES_MODEL(KTO)%XLES_RESOLVED_WP
XLES_RESOLVED_UTh=>LES_MODEL(KTO)%XLES_RESOLVED_UTh
XLES_RESOLVED_VTh=>LES_MODEL(KTO)%XLES_RESOLVED_VTh
XLES_RESOLVED_WTh=>LES_MODEL(KTO)%XLES_RESOLVED_WTh
XLES_RESOLVED_UThl=>LES_MODEL(KTO)%XLES_RESOLVED_UThl
XLES_RESOLVED_VThl=>LES_MODEL(KTO)%XLES_RESOLVED_VThl
XLES_RESOLVED_WThl=>LES_MODEL(KTO)%XLES_RESOLVED_WThl
XLES_RESOLVED_UThv=>LES_MODEL(KTO)%XLES_RESOLVED_UThv
XLES_RESOLVED_VThv=>LES_MODEL(KTO)%XLES_RESOLVED_VThv
XLES_RESOLVED_WThv=>LES_MODEL(KTO)%XLES_RESOLVED_WThv
XLES_RESOLVED_URv=>LES_MODEL(KTO)%XLES_RESOLVED_URv
XLES_RESOLVED_VRv=>LES_MODEL(KTO)%XLES_RESOLVED_VRv
XLES_RESOLVED_WRv=>LES_MODEL(KTO)%XLES_RESOLVED_WRv
XLES_RESOLVED_URc=>LES_MODEL(KTO)%XLES_RESOLVED_URc
XLES_RESOLVED_VRc=>LES_MODEL(KTO)%XLES_RESOLVED_VRc
XLES_RESOLVED_WRc=>LES_MODEL(KTO)%XLES_RESOLVED_WRc
XLES_RESOLVED_URi=>LES_MODEL(KTO)%XLES_RESOLVED_URi
XLES_RESOLVED_VRi=>LES_MODEL(KTO)%XLES_RESOLVED_VRi
XLES_RESOLVED_WRi=>LES_MODEL(KTO)%XLES_RESOLVED_WRi
XLES_RESOLVED_WRt=>LES_MODEL(KTO)%XLES_RESOLVED_WRt
XLES_RESOLVED_WRr=>LES_MODEL(KTO)%XLES_RESOLVED_WRr
XLES_INPRR3D=>LES_MODEL(KTO)%XLES_INPRR3D     
XLES_MAX_INPRR3D=>LES_MODEL(KTO)%XLES_MAX_INPRR3D     
XLES_EVAP3D=>LES_MODEL(KTO)%XLES_EVAP3D     
XLES_RESOLVED_USv=>LES_MODEL(KTO)%XLES_RESOLVED_USv
XLES_RESOLVED_VSv=>LES_MODEL(KTO)%XLES_RESOLVED_VSv
XLES_RESOLVED_WSv=>LES_MODEL(KTO)%XLES_RESOLVED_WSv
XLES_RESOLVED_U3=>LES_MODEL(KTO)%XLES_RESOLVED_U3
XLES_RESOLVED_V3=>LES_MODEL(KTO)%XLES_RESOLVED_V3
XLES_RESOLVED_W3=>LES_MODEL(KTO)%XLES_RESOLVED_W3
XLES_RESOLVED_U4=>LES_MODEL(KTO)%XLES_RESOLVED_U4
XLES_RESOLVED_V4=>LES_MODEL(KTO)%XLES_RESOLVED_V4
XLES_RESOLVED_W4=>LES_MODEL(KTO)%XLES_RESOLVED_W4
XLES_RESOLVED_Ua_ddxa_P=>LES_MODEL(KTO)%XLES_RESOLVED_Ua_ddxa_P
XLES_RESOLVED_ThlPz=>LES_MODEL(KTO)%XLES_RESOLVED_ThlPz
XLES_RESOLVED_WThl2=>LES_MODEL(KTO)%XLES_RESOLVED_WThl2
XLES_RESOLVED_W2Thl=>LES_MODEL(KTO)%XLES_RESOLVED_W2Thl
XLES_RESOLVED_WRv2=>LES_MODEL(KTO)%XLES_RESOLVED_WRv2
XLES_RESOLVED_W2Rv=>LES_MODEL(KTO)%XLES_RESOLVED_W2Rv
XLES_RESOLVED_WRt2=>LES_MODEL(KTO)%XLES_RESOLVED_WRt2
XLES_RESOLVED_W2Rt=>LES_MODEL(KTO)%XLES_RESOLVED_W2Rt
XLES_RESOLVED_RvPz=>LES_MODEL(KTO)%XLES_RESOLVED_RvPz
XLES_RESOLVED_RtPz=>LES_MODEL(KTO)%XLES_RESOLVED_RtPz
XLES_RESOLVED_WThlRv=>LES_MODEL(KTO)%XLES_RESOLVED_WThlRv
XLES_RESOLVED_WThlRt=>LES_MODEL(KTO)%XLES_RESOLVED_WThlRt
XLES_RESOLVED_WRc2=>LES_MODEL(KTO)%XLES_RESOLVED_WRc2
XLES_RESOLVED_W2Rc=>LES_MODEL(KTO)%XLES_RESOLVED_W2Rc
XLES_RESOLVED_RCPz=>LES_MODEL(KTO)%XLES_RESOLVED_RCPz
XLES_RESOLVED_WThlRc=>LES_MODEL(KTO)%XLES_RESOLVED_WThlRc
XLES_RESOLVED_WRvRc=>LES_MODEL(KTO)%XLES_RESOLVED_WRvRc
XLES_RESOLVED_WRi2=>LES_MODEL(KTO)%XLES_RESOLVED_WRi2
XLES_RESOLVED_W2Ri=>LES_MODEL(KTO)%XLES_RESOLVED_W2Ri
XLES_RESOLVED_RIPz=>LES_MODEL(KTO)%XLES_RESOLVED_RIPz
XLES_RESOLVED_WThlRi=>LES_MODEL(KTO)%XLES_RESOLVED_WThlRi
XLES_RESOLVED_WRvRi=>LES_MODEL(KTO)%XLES_RESOLVED_WRvRi
XLES_RESOLVED_WSv2=>LES_MODEL(KTO)%XLES_RESOLVED_WSv2
XLES_RESOLVED_W2Sv=>LES_MODEL(KTO)%XLES_RESOLVED_W2Sv
XLES_RESOLVED_SvPz=>LES_MODEL(KTO)%XLES_RESOLVED_SvPz
XLES_RESOLVED_WThlSv=>LES_MODEL(KTO)%XLES_RESOLVED_WThlSv
XLES_RESOLVED_WRvSv=>LES_MODEL(KTO)%XLES_RESOLVED_WRvSv
XLES_RESOLVED_MASSFX=>LES_MODEL(KTO)%XLES_RESOLVED_MASSFX
XLES_RESOLVED_UKe=>LES_MODEL(KTO)%XLES_RESOLVED_UKe
XLES_RESOLVED_VKe=>LES_MODEL(KTO)%XLES_RESOLVED_VKe
XLES_RESOLVED_WKe=>LES_MODEL(KTO)%XLES_RESOLVED_WKe
XLES_BU_RES_KE=>LES_MODEL(KTO)%XLES_BU_RES_KE
XLES_BU_RES_WThl=>LES_MODEL(KTO)%XLES_BU_RES_WThl
XLES_BU_RES_Thl2=>LES_MODEL(KTO)%XLES_BU_RES_Thl2
XLES_BU_RES_WRt=>LES_MODEL(KTO)%XLES_BU_RES_WRt
XLES_BU_RES_Rt2=>LES_MODEL(KTO)%XLES_BU_RES_Rt2
XLES_BU_RES_ThlRt=>LES_MODEL(KTO)%XLES_BU_RES_ThlRt
XLES_BU_RES_Sv2=>LES_MODEL(KTO)%XLES_BU_RES_Sv2
XLES_BU_RES_WSv=>LES_MODEL(KTO)%XLES_BU_RES_WSv
XLES_RES_U_SBG_Tke=>LES_MODEL(KTO)%XLES_RES_U_SBG_Tke
XLES_RES_V_SBG_Tke=>LES_MODEL(KTO)%XLES_RES_V_SBG_Tke
XLES_RES_W_SBG_Tke=>LES_MODEL(KTO)%XLES_RES_W_SBG_Tke
XLES_RES_W_SBG_WThl=>LES_MODEL(KTO)%XLES_RES_W_SBG_WThl
XLES_RES_W_SBG_WRt=>LES_MODEL(KTO)%XLES_RES_W_SBG_WRt
XLES_RES_W_SBG_Thl2=>LES_MODEL(KTO)%XLES_RES_W_SBG_Thl2
XLES_RES_W_SBG_Rt2=>LES_MODEL(KTO)%XLES_RES_W_SBG_Rt2
XLES_RES_W_SBG_ThlRt=>LES_MODEL(KTO)%XLES_RES_W_SBG_ThlRt
XLES_RES_W_SBG_WSv=>LES_MODEL(KTO)%XLES_RES_W_SBG_WSv
XLES_RES_W_SBG_Sv2=>LES_MODEL(KTO)%XLES_RES_W_SBG_Sv2
XLES_RES_ddxa_U_SBG_UaU=>LES_MODEL(KTO)%XLES_RES_ddxa_U_SBG_UaU
XLES_RES_ddxa_V_SBG_UaV=>LES_MODEL(KTO)%XLES_RES_ddxa_V_SBG_UaV
XLES_RES_ddxa_W_SBG_UaW=>LES_MODEL(KTO)%XLES_RES_ddxa_W_SBG_UaW
XLES_RES_ddxa_W_SBG_UaThl=>LES_MODEL(KTO)%XLES_RES_ddxa_W_SBG_UaThl
XLES_RES_ddxa_Thl_SBG_UaW=>LES_MODEL(KTO)%XLES_RES_ddxa_Thl_SBG_UaW
XLES_RES_ddz_Thl_SBG_W2=>LES_MODEL(KTO)%XLES_RES_ddz_Thl_SBG_W2
XLES_RES_ddxa_W_SBG_UaRt=>LES_MODEL(KTO)%XLES_RES_ddxa_W_SBG_UaRt
XLES_RES_ddxa_Rt_SBG_UaW=>LES_MODEL(KTO)%XLES_RES_ddxa_Rt_SBG_UaW
XLES_RES_ddz_Rt_SBG_W2=>LES_MODEL(KTO)%XLES_RES_ddz_Rt_SBG_W2
XLES_RES_ddxa_Thl_SBG_UaRt=>LES_MODEL(KTO)%XLES_RES_ddxa_Thl_SBG_UaRt
XLES_RES_ddxa_Rt_SBG_UaThl=>LES_MODEL(KTO)%XLES_RES_ddxa_Rt_SBG_UaThl
XLES_RES_ddxa_Thl_SBG_UaThl=>LES_MODEL(KTO)%XLES_RES_ddxa_Thl_SBG_UaThl
XLES_RES_ddxa_Rt_SBG_UaRt=>LES_MODEL(KTO)%XLES_RES_ddxa_Rt_SBG_UaRt
XLES_RES_ddxa_W_SBG_UaSv=>LES_MODEL(KTO)%XLES_RES_ddxa_W_SBG_UaSv
XLES_RES_ddxa_Sv_SBG_UaW=>LES_MODEL(KTO)%XLES_RES_ddxa_Sv_SBG_UaW
XLES_RES_ddz_Sv_SBG_W2=>LES_MODEL(KTO)%XLES_RES_ddz_Sv_SBG_W2
XLES_RES_ddxa_Sv_SBG_UaSv=>LES_MODEL(KTO)%XLES_RES_ddxa_Sv_SBG_UaSv
XLES_BU_SBG_TKE=>LES_MODEL(KTO)%XLES_BU_SBG_TKE
XLES_SUBGRID_U2=>LES_MODEL(KTO)%XLES_SUBGRID_U2
XLES_SUBGRID_V2=>LES_MODEL(KTO)%XLES_SUBGRID_V2
XLES_SUBGRID_W2=>LES_MODEL(KTO)%XLES_SUBGRID_W2
XLES_SUBGRID_Tke=>LES_MODEL(KTO)%XLES_SUBGRID_Tke
XLES_SUBGRID_Thl2=>LES_MODEL(KTO)%XLES_SUBGRID_Thl2
XLES_SUBGRID_Rt2=>LES_MODEL(KTO)%XLES_SUBGRID_Rt2
XLES_SUBGRID_Rc2=>LES_MODEL(KTO)%XLES_SUBGRID_Rc2
XLES_SUBGRID_Ri2=>LES_MODEL(KTO)%XLES_SUBGRID_Ri2
XLES_SUBGRID_ThlRt=>LES_MODEL(KTO)%XLES_SUBGRID_ThlRt
XLES_SUBGRID_Sv2=>LES_MODEL(KTO)%XLES_SUBGRID_Sv2
XLES_SUBGRID_UV=>LES_MODEL(KTO)%XLES_SUBGRID_UV
XLES_SUBGRID_WU=>LES_MODEL(KTO)%XLES_SUBGRID_WU
XLES_SUBGRID_WV=>LES_MODEL(KTO)%XLES_SUBGRID_WV
XLES_SUBGRID_UThl=>LES_MODEL(KTO)%XLES_SUBGRID_UThl
XLES_SUBGRID_VThl=>LES_MODEL(KTO)%XLES_SUBGRID_VThl
XLES_SUBGRID_WThl=>LES_MODEL(KTO)%XLES_SUBGRID_WThl
XLES_SUBGRID_ENTLES=>LES_MODEL(KTO)%XLES_SUBGRID_ENTLES
XLES_SUBGRID_DETLES=>LES_MODEL(KTO)%XLES_SUBGRID_DETLES
XLES_SUBGRID_URt=>LES_MODEL(KTO)%XLES_SUBGRID_URt
XLES_SUBGRID_VRt=>LES_MODEL(KTO)%XLES_SUBGRID_VRt
XLES_SUBGRID_WRt=>LES_MODEL(KTO)%XLES_SUBGRID_WRt
XLES_SUBGRID_URc=>LES_MODEL(KTO)%XLES_SUBGRID_URc
XLES_SUBGRID_VRc=>LES_MODEL(KTO)%XLES_SUBGRID_VRc
XLES_SUBGRID_WRc=>LES_MODEL(KTO)%XLES_SUBGRID_WRc
XLES_SUBGRID_USv=>LES_MODEL(KTO)%XLES_SUBGRID_USv
XLES_SUBGRID_VSv=>LES_MODEL(KTO)%XLES_SUBGRID_VSv
XLES_SUBGRID_WSv=>LES_MODEL(KTO)%XLES_SUBGRID_WSv
XLES_SUBGRID_UTke=>LES_MODEL(KTO)%XLES_SUBGRID_UTke
XLES_SUBGRID_VTke=>LES_MODEL(KTO)%XLES_SUBGRID_VTke
XLES_SUBGRID_WTke=>LES_MODEL(KTO)%XLES_SUBGRID_WTke
XLES_SUBGRID_ddz_WTke=>LES_MODEL(KTO)%XLES_SUBGRID_ddz_WTke
XLES_SUBGRID_WThv=>LES_MODEL(KTO)%XLES_SUBGRID_WThv
XLES_SUBGRID_ThlThv=>LES_MODEL(KTO)%XLES_SUBGRID_ThlThv
XLES_SUBGRID_RtThv=>LES_MODEL(KTO)%XLES_SUBGRID_RtThv
XLES_SUBGRID_SvThv=>LES_MODEL(KTO)%XLES_SUBGRID_SvThv
XLES_SUBGRID_W2Thl=>LES_MODEL(KTO)%XLES_SUBGRID_W2Thl
XLES_SUBGRID_W2Rt=>LES_MODEL(KTO)%XLES_SUBGRID_W2Rt
XLES_SUBGRID_WThlRt=>LES_MODEL(KTO)%XLES_SUBGRID_WThlRt
XLES_SUBGRID_W2Sv=>LES_MODEL(KTO)%XLES_SUBGRID_W2Sv
XLES_SUBGRID_WThl2=>LES_MODEL(KTO)%XLES_SUBGRID_WThl2
XLES_SUBGRID_WRt2=>LES_MODEL(KTO)%XLES_SUBGRID_WRt2
XLES_SUBGRID_WSv2=>LES_MODEL(KTO)%XLES_SUBGRID_WSv2
XLES_SUBGRID_DISS_Tke=>LES_MODEL(KTO)%XLES_SUBGRID_DISS_Tke
XLES_SUBGRID_DISS_Thl2=>LES_MODEL(KTO)%XLES_SUBGRID_DISS_Thl2
XLES_SUBGRID_DISS_Rt2=>LES_MODEL(KTO)%XLES_SUBGRID_DISS_Rt2
XLES_SUBGRID_DISS_ThlRt=>LES_MODEL(KTO)%XLES_SUBGRID_DISS_ThlRt
XLES_SUBGRID_DISS_Sv2=>LES_MODEL(KTO)%XLES_SUBGRID_DISS_Sv2
XLES_SUBGRID_WP=>LES_MODEL(KTO)%XLES_SUBGRID_WP
XLES_SUBGRID_ThlPz=>LES_MODEL(KTO)%XLES_SUBGRID_ThlPz
XLES_SUBGRID_RtPz=>LES_MODEL(KTO)%XLES_SUBGRID_RtPz
XLES_SUBGRID_SvPz=>LES_MODEL(KTO)%XLES_SUBGRID_SvPz
XLES_SUBGRID_PHI3=>LES_MODEL(KTO)%XLES_SUBGRID_PHI3
XLES_SUBGRID_PSI3=>LES_MODEL(KTO)%XLES_SUBGRID_PSI3
XLES_SUBGRID_LMix=>LES_MODEL(KTO)%XLES_SUBGRID_LMix
XLES_SUBGRID_LDiss=>LES_MODEL(KTO)%XLES_SUBGRID_LDiss
XLES_SUBGRID_Km=>LES_MODEL(KTO)%XLES_SUBGRID_Km
XLES_SUBGRID_Kh=>LES_MODEL(KTO)%XLES_SUBGRID_Kh

XLES_SUBGRID_THLUP_MF=>LES_MODEL(KTO)%XLES_SUBGRID_THLUP_MF
XLES_SUBGRID_RTUP_MF=>LES_MODEL(KTO)%XLES_SUBGRID_RTUP_MF
XLES_SUBGRID_RVUP_MF=>LES_MODEL(KTO)%XLES_SUBGRID_RVUP_MF
XLES_SUBGRID_RCUP_MF=>LES_MODEL(KTO)%XLES_SUBGRID_RCUP_MF
XLES_SUBGRID_RIUP_MF=>LES_MODEL(KTO)%XLES_SUBGRID_RIUP_MF
XLES_SUBGRID_WUP_MF=>LES_MODEL(KTO)%XLES_SUBGRID_WUP_MF
XLES_SUBGRID_MASSFLUX=>LES_MODEL(KTO)%XLES_SUBGRID_MASSFLUX
XLES_SUBGRID_DETR=>LES_MODEL(KTO)%XLES_SUBGRID_DETR
XLES_SUBGRID_ENTR=>LES_MODEL(KTO)%XLES_SUBGRID_ENTR
XLES_SUBGRID_FRACUP=>LES_MODEL(KTO)%XLES_SUBGRID_FRACUP
XLES_SUBGRID_THVUP_MF=>LES_MODEL(KTO)%XLES_SUBGRID_THVUP_MF
XLES_SUBGRID_WTHLMF=>LES_MODEL(KTO)%XLES_SUBGRID_WTHLMF
XLES_SUBGRID_WRTMF=>LES_MODEL(KTO)%XLES_SUBGRID_WRTMF
XLES_SUBGRID_WTHVMF=>LES_MODEL(KTO)%XLES_SUBGRID_WTHVMF
XLES_SUBGRID_WUMF=>LES_MODEL(KTO)%XLES_SUBGRID_WUMF
XLES_SUBGRID_WVMF=>LES_MODEL(KTO)%XLES_SUBGRID_WVMF

XLES_UPDRAFT=>LES_MODEL(KTO)%XLES_UPDRAFT
XLES_UPDRAFT_W=>LES_MODEL(KTO)%XLES_UPDRAFT_W
XLES_UPDRAFT_Th=>LES_MODEL(KTO)%XLES_UPDRAFT_Th
XLES_UPDRAFT_Thv=>LES_MODEL(KTO)%XLES_UPDRAFT_Thv
XLES_UPDRAFT_Thl=>LES_MODEL(KTO)%XLES_UPDRAFT_Thl
XLES_UPDRAFT_Rv=>LES_MODEL(KTO)%XLES_UPDRAFT_Rv
XLES_UPDRAFT_Rc=>LES_MODEL(KTO)%XLES_UPDRAFT_Rc
XLES_UPDRAFT_Rr=>LES_MODEL(KTO)%XLES_UPDRAFT_Rr
XLES_UPDRAFT_Ri=>LES_MODEL(KTO)%XLES_UPDRAFT_Ri
XLES_UPDRAFT_Rs=>LES_MODEL(KTO)%XLES_UPDRAFT_Rs
XLES_UPDRAFT_Rg=>LES_MODEL(KTO)%XLES_UPDRAFT_Rg
XLES_UPDRAFT_Rh=>LES_MODEL(KTO)%XLES_UPDRAFT_Rh
XLES_UPDRAFT_Sv=>LES_MODEL(KTO)%XLES_UPDRAFT_Sv
XLES_UPDRAFT_Tke=>LES_MODEL(KTO)%XLES_UPDRAFT_Tke
XLES_UPDRAFT_Ke=>LES_MODEL(KTO)%XLES_UPDRAFT_Ke
XLES_UPDRAFT_WTh=>LES_MODEL(KTO)%XLES_UPDRAFT_WTh
XLES_UPDRAFT_WThl=>LES_MODEL(KTO)%XLES_UPDRAFT_WThl
XLES_UPDRAFT_WThv=>LES_MODEL(KTO)%XLES_UPDRAFT_WThv
XLES_UPDRAFT_WRv=>LES_MODEL(KTO)%XLES_UPDRAFT_WRv
XLES_UPDRAFT_WRc=>LES_MODEL(KTO)%XLES_UPDRAFT_WRc
XLES_UPDRAFT_WRi=>LES_MODEL(KTO)%XLES_UPDRAFT_WRi
XLES_UPDRAFT_WSv=>LES_MODEL(KTO)%XLES_UPDRAFT_WSv
XLES_UPDRAFT_Th2=>LES_MODEL(KTO)%XLES_UPDRAFT_Th2
XLES_UPDRAFT_Thl2=>LES_MODEL(KTO)%XLES_UPDRAFT_Thl2
XLES_UPDRAFT_ThThv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThThv
XLES_UPDRAFT_ThlThv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThlThv
XLES_UPDRAFT_Rv2=>LES_MODEL(KTO)%XLES_UPDRAFT_Rv2
XLES_UPDRAFT_Rc2=>LES_MODEL(KTO)%XLES_UPDRAFT_Rc2
XLES_UPDRAFT_Ri2=>LES_MODEL(KTO)%XLES_UPDRAFT_Ri2
XLES_UPDRAFT_Sv2=>LES_MODEL(KTO)%XLES_UPDRAFT_Sv2
XLES_UPDRAFT_ThRv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThRv
XLES_UPDRAFT_ThRc=>LES_MODEL(KTO)%XLES_UPDRAFT_ThRc
XLES_UPDRAFT_ThRi=>LES_MODEL(KTO)%XLES_UPDRAFT_ThRi
XLES_UPDRAFT_ThSv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThSv
XLES_UPDRAFT_ThvRv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThvRv
XLES_UPDRAFT_ThvRc=>LES_MODEL(KTO)%XLES_UPDRAFT_ThvRc
XLES_UPDRAFT_ThvRi=>LES_MODEL(KTO)%XLES_UPDRAFT_ThvRi
XLES_UPDRAFT_ThvSv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThvSv
XLES_UPDRAFT_ThlRv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThlRv
XLES_UPDRAFT_ThlRc=>LES_MODEL(KTO)%XLES_UPDRAFT_ThlRc
XLES_UPDRAFT_ThlRi=>LES_MODEL(KTO)%XLES_UPDRAFT_ThlRi
XLES_UPDRAFT_ThlSv=>LES_MODEL(KTO)%XLES_UPDRAFT_ThlSv
XLES_DOWNDRAFT=>LES_MODEL(KTO)%XLES_DOWNDRAFT
XLES_DOWNDRAFT_W=>LES_MODEL(KTO)%XLES_DOWNDRAFT_W
XLES_DOWNDRAFT_Th=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Th
XLES_DOWNDRAFT_Thv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Thv
XLES_DOWNDRAFT_Thl=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Thl
XLES_DOWNDRAFT_Rv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rv
XLES_DOWNDRAFT_Rc=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rc
XLES_DOWNDRAFT_Rr=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rr
XLES_DOWNDRAFT_Ri=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Ri
XLES_DOWNDRAFT_Rs=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rs
XLES_DOWNDRAFT_Rg=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rg
XLES_DOWNDRAFT_Rh=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rh
XLES_DOWNDRAFT_Sv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Sv
XLES_DOWNDRAFT_Tke=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Tke
XLES_DOWNDRAFT_Ke=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Ke
XLES_DOWNDRAFT_WTh=>LES_MODEL(KTO)%XLES_DOWNDRAFT_WTh
XLES_DOWNDRAFT_WThl=>LES_MODEL(KTO)%XLES_DOWNDRAFT_WThl
XLES_DOWNDRAFT_WThv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_WThv
XLES_DOWNDRAFT_WRv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_WRv
XLES_DOWNDRAFT_WRc=>LES_MODEL(KTO)%XLES_DOWNDRAFT_WRc
XLES_DOWNDRAFT_WRi=>LES_MODEL(KTO)%XLES_DOWNDRAFT_WRi
XLES_DOWNDRAFT_WSv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_WSv
XLES_DOWNDRAFT_Th2=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Th2
XLES_DOWNDRAFT_Thl2=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Thl2
XLES_DOWNDRAFT_ThThv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThThv
XLES_DOWNDRAFT_ThlThv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThlThv
XLES_DOWNDRAFT_Rv2=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rv2
XLES_DOWNDRAFT_Rc2=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Rc2
XLES_DOWNDRAFT_Ri2=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Ri2
XLES_DOWNDRAFT_Sv2=>LES_MODEL(KTO)%XLES_DOWNDRAFT_Sv2
XLES_DOWNDRAFT_ThRv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThRv
XLES_DOWNDRAFT_ThRc=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThRc
XLES_DOWNDRAFT_ThRi=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThRi
XLES_DOWNDRAFT_ThSv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThSv
XLES_DOWNDRAFT_ThvRv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThvRv
XLES_DOWNDRAFT_ThvRc=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThvRc
XLES_DOWNDRAFT_ThvRi=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThvRi
XLES_DOWNDRAFT_ThvSv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThvSv
XLES_DOWNDRAFT_ThlRv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThlRv
XLES_DOWNDRAFT_ThlRc=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThlRc
XLES_DOWNDRAFT_ThlRi=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThlRi
XLES_DOWNDRAFT_ThlSv=>LES_MODEL(KTO)%XLES_DOWNDRAFT_ThlSv
XLES_UW0=>LES_MODEL(KTO)%XLES_UW0
XLES_VW0=>LES_MODEL(KTO)%XLES_VW0
XLES_USTAR=>LES_MODEL(KTO)%XLES_USTAR
XLES_WSTAR=>LES_MODEL(KTO)%XLES_WSTAR
XLES_Q0=>LES_MODEL(KTO)%XLES_Q0
XLES_E0=>LES_MODEL(KTO)%XLES_E0
XLES_SV0=>LES_MODEL(KTO)%XLES_SV0
XLES_BL_HEIGHT=>LES_MODEL(KTO)%XLES_BL_HEIGHT
XLES_MO_LENGTH=>LES_MODEL(KTO)%XLES_MO_LENGTH
XLES_ZCB=>LES_MODEL(KTO)%XLES_ZCB
XLES_CFtot=>LES_MODEL(KTO)%XLES_CFtot
XLES_CF2tot=>LES_MODEL(KTO)%XLES_CF2tot
XLES_LWP=>LES_MODEL(KTO)%XLES_LWP
XLES_LWPVAR=>LES_MODEL(KTO)%XLES_LWPVAR
XLES_RWP=>LES_MODEL(KTO)%XLES_RWP
XLES_IWP=>LES_MODEL(KTO)%XLES_IWP
XLES_SWP=>LES_MODEL(KTO)%XLES_SWP
XLES_GWP=>LES_MODEL(KTO)%XLES_GWP
XLES_HWP=>LES_MODEL(KTO)%XLES_HWP
XLES_INT_TKE=>LES_MODEL(KTO)%XLES_INT_TKE
XLES_INPRR=>LES_MODEL(KTO)%XLES_INPRR
XLES_INPRC=>LES_MODEL(KTO)%XLES_INPRC
XLES_INDEP=>LES_MODEL(KTO)%XLES_INDEP
XLES_RAIN_INPRR=>LES_MODEL(KTO)%XLES_RAIN_INPRR
XLES_ACPRR=>LES_MODEL(KTO)%XLES_ACPRR
XLES_ZMAXCF=>LES_MODEL(KTO)%XLES_ZMAXCF
XLES_ZMAXCF2=>LES_MODEL(KTO)%XLES_ZMAXCF2
XLES_PRECFR=>LES_MODEL(KTO)%XLES_PRECFR
XCORRi_UU=>LES_MODEL(KTO)%XCORRi_UU
XCORRi_VV=>LES_MODEL(KTO)%XCORRi_VV
XCORRi_UV=>LES_MODEL(KTO)%XCORRi_UV
XCORRi_WU=>LES_MODEL(KTO)%XCORRi_WU
XCORRi_WV=>LES_MODEL(KTO)%XCORRi_WV
XCORRi_WW=>LES_MODEL(KTO)%XCORRi_WW
XCORRi_WTh=>LES_MODEL(KTO)%XCORRi_WTh
XCORRi_WThl=>LES_MODEL(KTO)%XCORRi_WThl
XCORRi_WRv=>LES_MODEL(KTO)%XCORRi_WRv
XCORRi_WRc=>LES_MODEL(KTO)%XCORRi_WRc
XCORRi_WRi=>LES_MODEL(KTO)%XCORRi_WRi
XCORRi_WSv=>LES_MODEL(KTO)%XCORRi_WSv
XCORRi_ThTh=>LES_MODEL(KTO)%XCORRi_ThTh
XCORRi_ThRv=>LES_MODEL(KTO)%XCORRi_ThRv
XCORRi_ThRc=>LES_MODEL(KTO)%XCORRi_ThRc
XCORRi_ThRi=>LES_MODEL(KTO)%XCORRi_ThRi
XCORRi_ThlThl=>LES_MODEL(KTO)%XCORRi_ThlThl
XCORRi_ThlRv=>LES_MODEL(KTO)%XCORRi_ThlRv
XCORRi_ThlRc=>LES_MODEL(KTO)%XCORRi_ThlRc
XCORRi_ThlRi=>LES_MODEL(KTO)%XCORRi_ThlRi
XCORRi_RvRv=>LES_MODEL(KTO)%XCORRi_RvRv
XCORRi_RcRc=>LES_MODEL(KTO)%XCORRi_RcRc
XCORRi_RiRi=>LES_MODEL(KTO)%XCORRi_RiRi
XCORRi_SvSv=>LES_MODEL(KTO)%XCORRi_SvSv
XCORRj_UU=>LES_MODEL(KTO)%XCORRj_UU
XCORRj_VV=>LES_MODEL(KTO)%XCORRj_VV
XCORRj_UV=>LES_MODEL(KTO)%XCORRj_UV
XCORRj_WU=>LES_MODEL(KTO)%XCORRj_WU
XCORRj_WV=>LES_MODEL(KTO)%XCORRj_WV
XCORRj_WW=>LES_MODEL(KTO)%XCORRj_WW
XCORRj_WTh=>LES_MODEL(KTO)%XCORRj_WTh
XCORRj_WThl=>LES_MODEL(KTO)%XCORRj_WThl
XCORRj_WRv=>LES_MODEL(KTO)%XCORRj_WRv
XCORRj_WRc=>LES_MODEL(KTO)%XCORRj_WRc
XCORRj_WRi=>LES_MODEL(KTO)%XCORRj_WRi
XCORRj_WSv=>LES_MODEL(KTO)%XCORRj_WSv
XCORRj_ThTh=>LES_MODEL(KTO)%XCORRj_ThTh
XCORRj_ThRv=>LES_MODEL(KTO)%XCORRj_ThRv
XCORRj_ThRc=>LES_MODEL(KTO)%XCORRj_ThRc
XCORRj_ThRi=>LES_MODEL(KTO)%XCORRj_ThRi
XCORRj_ThlThl=>LES_MODEL(KTO)%XCORRj_ThlThl
XCORRj_ThlRv=>LES_MODEL(KTO)%XCORRj_ThlRv
XCORRj_ThlRc=>LES_MODEL(KTO)%XCORRj_ThlRc
XCORRj_ThlRi=>LES_MODEL(KTO)%XCORRj_ThlRi
XCORRj_RvRv=>LES_MODEL(KTO)%XCORRj_RvRv
XCORRj_RcRc=>LES_MODEL(KTO)%XCORRj_RcRc
XCORRj_RiRi=>LES_MODEL(KTO)%XCORRj_RiRi
XCORRj_SvSv=>LES_MODEL(KTO)%XCORRj_SvSv
XCORRk_WW=>LES_MODEL(KTO)%XCORRk_WW
XLES_SWU=>LES_MODEL(KTO)%XLES_SWU
XLES_SWD=>LES_MODEL(KTO)%XLES_SWD
XLES_LWU=>LES_MODEL(KTO)%XLES_LWU
XLES_LWD=>LES_MODEL(KTO)%XLES_LWD
XLES_DTHRADSW=>LES_MODEL(KTO)%XLES_DTHRADSW
XLES_DTHRADLW=>LES_MODEL(KTO)%XLES_DTHRADLW
XLES_RADEFF=>LES_MODEL(KTO)%XLES_RADEFF

END SUBROUTINE LES_GOTO_MODEL

END MODULE MODD_LES_n
