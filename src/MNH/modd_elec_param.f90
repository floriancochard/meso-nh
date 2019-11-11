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
!      #######################
       MODULE  MODD_ELEC_PARAM
!      #######################
!
!!****  *MODD_ELEC_PARAM* - declaration of some electrical factors
!!			    extensively used in the electrical scheme.
!!
!!      PURPOSE
!!      -------
!	  The purpose of this declarative module is to declare some precomputed
!	electrical parameters directly used in routine RAIN_ICE_ELEC.
!!
!!**    IMPLICIT ARGUMENTS
!!      ------------------
!!	  None
!!
!!      REFERENCE
!!      ---------
!!
!!      AUTHOR
!!      ------
!!        Gilles Molinie    * Laboratoire d'Aerologie *
!!        
!!
!!      MODIFICATIONS
!!      -------------
!!        Original      14/11/02
!!
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
REAL, SAVE :: XCOEF_RQ_V, XCOEF_RQ_C, &   ! Constants for proportionality
              XCOEF_RQ_R, XCOEF_RQ_I, &   ! between mass transfer and
              XCOEF_RQ_S, XCOEF_RQ_G, &   ! charge transfer
              XCOEF_RQ_H
!
REAL :: XEGMIN, XEGMAX, XESMIN, XESMAX, &  ! Max and min values for
        XEIMIN, XEIMAX, XECMIN, XECMAX, &  ! e_x in q=e_x D^f_x
        XERMIN, XERMAX, XEHMIN, XEHMAX
!
REAL, SAVE :: XQHON    ! Constant for spontaneous freezing of droplets if T<-35°
!
REAL, SAVE :: XFQSEDR, XEXQSEDR,    &    ! Constant for sedimentation of rain
              XFQSEDI, XEXQSEDI,    &    !  ice
              XFQSEDS, XEXQSEDS,    &    !  snow
              XFQSEDG, XEXQSEDG,    &    !  graupel
              XFQSEDH, XEXQSEDH          !  hail
REAL, SAVE :: XFCI    ! Constant for sedimentation of the mixing ratio of ice
                      ! which the computation is modified in regard of rain_ice.f90
!
REAL, SAVE :: XQSRIMCG, XEXQSRIMCG       ! Constant for riming of cloud droplets
                                         ! on snow 
REAL, DIMENSION(:), SAVE, ALLOCATABLE :: XGAMINC_RIM3  
!
REAL, SAVE :: XQRCFRIG, XEXQRCFRIG       ! Constant for contact freezing between
                                         ! raindrops and pristine ice
REAL, SAVE :: XFQRACCS                   ! Constant in RACCS
!
REAL, SAVE :: XFQIAGGSBH,                &         ! Constant for IAGGS charging
              XFQIAGGSBG, XEXFQIAGGSBG,  &         ! process for HELFA, GARDI,
              XFQIAGGSBS,                &         ! SAUND and TAKAH
              XFQIAGGSBT1, XFQIAGGSBT2, XFQIAGGSBT3 
!
REAL, SAVE :: XLBQRACCS1, XLBQRACCS2, XLBQRACCS3    ! Integral of normalization 
REAL, SAVE :: XLBQSACCRG1, XLBQSACCRG2, XLBQSACCRG3 ! in accretion of raindrops
                                                    ! on snow process 
!                         
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE  &               ! Normalized kernel for
           :: XKER_Q_RACCS, XKER_Q_RACCSS, XKER_Q_SACCRG ! RACCS, RACCSS, SACCRG
!
REAL, SAVE :: XFQSDRYG, XFQSDRYGB, XFQRDRYG        ! Constant in SDRYG and RDRYG
!
! charge separation
!
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE :: XKER_Q_LIMSG
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE  &                    
           :: XKER_Q_SDRYGB, XKER_Q_SDRYGB1, XKER_Q_SDRYGB2 
!
! Helsdon-Farley
!
REAL       :: XHIDRYG     ! Constant charge separated 
REAL       :: XHSDRYG     ! Constant charge separated 
REAL, SAVE :: XLBQSDRYGB4H, XLBQSDRYGB5H, XLBQSDRYGB6H    ! Constants in QIDRYGB
REAL, SAVE :: XFQSDRYGBH                                  !
!
! Gardiner
!
REAL, SAVE :: XLWCC  ! LWC critic in Gardiner NI charging
REAL, SAVE :: XFQIDRYGBG, XLBQIDRYGBG                     ! Constants in QIDRYGB
REAL, SAVE :: XFQSDRYGBG                                  ! Constants in QSDRYGB
REAL, SAVE :: XLBQSDRYGB4G, XLBQSDRYGB5G, XLBQSDRYGB6G    ! 
!
! Saunders
!
REAL, SAVE :: XIMP, XINP, XIKP, &   ! Parameters m, n and k
              XIMN, XINN, XIKN, &   ! for the NI processes
              XSMP, XSNP, XSKP, &   ! following
              XSMN, XSNN, XSKN      ! Saunders et al. (1991)
REAL, SAVE :: XFQIAGGSP, XFQIAGGSN,         & ! Auxiliary parameters
              XFQIDRYGBSP, XFQIDRYGBSN,     & ! containing MOMG function
              XLBQSDRYGB1SP, XLBQSDRYGB1SN, &
              XLBQSDRYGB2SP, XLBQSDRYGB2SN, &
              XLBQSDRYGB3SP, XLBQSDRYGB3SN, &
              XAIGAMMABI
REAL, SAVE :: XIKP_TAK, XIKN_TAK, XSKP_TAK, XSKN_TAK ! Using Takahashi charge
REAL, SAVE :: XFQIAGGSP_TAK, XFQIAGGSN_TAK, XFQIDRYGBSP_TAK, XFQIDRYGBSN_TAK
REAL, SAVE :: XVSCOEF, XVGCOEF
REAL, SAVE :: XFQIDRYGBS,   XLBQIDRYGBS    ! Constants in QIDRYGB
REAL, SAVE :: XFQSDRYGBS                   ! Constants in QSDRYGB
REAL, SAVE :: XLBQSDRYGB1S, XLBQSDRYGB2S   !
!
! Takahashi
!
INTEGER, SAVE :: NIND_TEMP ! number of indexes for temperature
INTEGER, SAVE :: NIND_LWC  ! number of indexes for liquid water content
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE :: XMANSELL ! F(LWC, T) for Takahashi(1978) /Mansell
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE :: XSAUNDER ! F(LWC, T) for SAUN1/SAUN2, BSMP1/BSMP2
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE :: XTAKA_TM ! F(LWC, T) for Takahashi/Tsenova and Mitzeva
!
REAL, SAVE :: XFQIDRYGBT1,  XFQIDRYGBT2,  XFQIDRYGBT3, &  ! IDRYGB
              XFQSDRYGBT1,  XFQSDRYGBT2,  XFQSDRYGBT3, &  ! SDRYGB
              XFQSDRYGBT4,  XFQSDRYGBT5,  XFQSDRYGBT6, &  ! SDRYGB
              XFQSDRYGBT7,  XFQSDRYGBT8,  XFQSDRYGBT9, &  ! SDRYGB
              XFQSDRYGBT10, XFQSDRYGBT11, XFQSDRYGBT12    ! SDRYGB
!
REAL, SAVE :: XLBQRDRYG1, XLBQRDRYG2, XLBQRDRYG3  ! Integral of normalization in
REAL, SAVE :: XLBQSDRYG1, XLBQSDRYG2, XLBQSDRYG3  ! the accretion of graupel on
                                                  ! raindrop and snow process
!
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE &  
           :: XKER_Q_SDRYG, XKER_Q_RDRYG ! Normalized kernel for SDRYG and RDRYG
!
REAL, SAVE :: XFQUPDC, XFQUPDR, XFQUPDI,&         ! Update Q=f(D)
              XEXFQUPDI, XFQUPDS, XFQUPDG, XFQUPDH
!
REAL, SAVE :: XQREVAV1, XQREVAV2                  ! Raindrops evaporation
!
! Add variables to limit the exchanged charge
!
REAL, SAVE :: XAUX_LIM
REAL, SAVE :: XAUX_LIM1, XAUX_LIM2, XAUX_LIM3
!
!
! Inductive charging process
!
REAL, SAVE :: XCOLCG_IND     ! collision effiency
REAL, SAVE :: XEBOUND    ! rebound efficiency
REAL, SAVE :: XALPHA_IND     ! fraction of droplets with grazing trajectories
REAL, SAVE :: XCOS_THETA ! average cosine of the angle of rebounding collision
REAL, SAVE :: XIND1, XIND2, XIND3
!
! lightning
!
REAL :: XFQLIGHTC, XFQLIGHTR, XFQLIGHTI, &
        XFQLIGHTS, XFQLIGHTG, XFQLIGHTH     ! Constant for charge redistribution
REAL :: XEXQLIGHTR, XEXQLIGHTI, &
        XEXQLIGHTS, XEXQLIGHTG, XEXQLIGHTH  ! Exponent for charge redistribution
!
END MODULE  MODD_ELEC_PARAM     
