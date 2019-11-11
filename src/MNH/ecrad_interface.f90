!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/radiations.f90,v $ $Revision: 1.3.2.3.2.2.2.4 $
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!    ########################
     MODULE MODI_ECRAD_INTERFACE  
!    ########################
!
CONTAINS
!
!##############################################################
!OPTION! -Ni
SUBROUTINE ECRAD_INTERFACE ( KLON,KLEV,KRAD_DIAG, KAER, &
     PDZ,HEFRADL, HEFRADI, HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,  &
     PRII0, PAER , PALBD , PALBP, PAPH , PAP,                 &
     PCCO2, PCLFR , PDP  , PEMIS, PEMIW , PLSM , PMU0, POZON, &
     PQ   , PQIWC, PIWC, PQLWC, PLWC,PQS  , PQRWC, PRWC,              &
     PTH  , PT    , PTS, PCCT_C2R2, PCRT_C2R2, PCIT_C1R3,     &
     PFLT  , PFLS , PFCT , PFCS  ,                            &  
     PDTLW, PDTSW ,PFLUX_TOP_GND_IRVISNIR,                    & 
     PSFSWDIR, PSFSWDIF,                                      &
     PFSDWN, PFSUP, PFLUX_LW ,                                &
     PDTLW_CS, PDTSW_CS ,PFLUX_TOP_GND_IRVISNIR_CS,           &
     PFCDWN, PFCUP, PFLUX_CLW,                                & 
     PPLAN_ALB_VIS, PPLAN_ALB_NIR, PPLAN_TRA_VIS, PPLAN_TRA_NIR,&
     PPLAN_ABS_VIS, PPLAN_ABS_NIR, PEFCL_LWD, PEFCL_LWU,      &
     PFLWP,PFIWP,PRADLP,PRADIP, PEFCL_RRTM, PCLSW_TOTAL,      & 
     PTAU_TOTAL, POMEGA_TOTAL, PCG_TOTAL,                     &
     ODUST,PPIZA_DST,PCGA_DST,PTAUREL_DST,PLAT,PLON           )
!##############################################################
!  
!**** *ECRAD* - INTERFACE TO ECRAD RADIATION SCHEMES USING THE SAME INPUTS AS ECMWF_RADIATION_VERS2
!
!     PURPOSE.
!     --------
!     Interface routine to call ECRAD SW and LW radiation schemes in the MesoNH 
!     framework. It is an alternative to ECMWF_RADIATION_VERS2.
!
!     METHOD.
!     -------
!     This routine calls the RADIATION_SCHEME routine from ECRAD package.
!     It converts the standard inputs common with ECMWF_VERSION_2 into those
!     required by RADIATION_SCHEME and then calls the latter routine.

!     It takes as inputs fields from routine RADIATIONS where index 0 is TOA
!     RADIATION_SCHEME takes such arrays and returns similar order. The order is reversed at the end of
!     the routine to return to RADIATIONS arrays from surface to top 
!              
!
!     EXTERNALS.
!     ----------
!     Note that external specific ECMWF module (named "yo==="). They don't follow the MNH MODD norm
!     and are initialised in the ini_radiation step. 
!
!
!     REFERENCE.
!     ----------
!     ECMWF and MESONH radiation part documentation     
!
!     AUTHORS.
!     --------
!        Q. LIBOIS (June 2017)  

!     MODIFICATIONS.
!     --------------     
!       J. Escobar 18/06/2018: bug compile R*4 => REAL(KIND=JPRB) PLAT/LON
!-----------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!ECMWF radiation scheme specific modules 
!
! USE YOMCST   , ONLY : RG ,RD ,RTT ,RPI
! USE YOERAD   , ONLY : NMODE, NSW ,LRRTM ,LINHOM ,LRADIP, LRADLP 
! USE YOELW    , ONLY : NSIL     ,NTRA     ,NUA      ,TSTAND   ,XP
! USE YOESW    , ONLY : RYFWCA   ,RYFWCB   ,RYFWCC   ,RYFWCD   ,&
!             &RYFWCE   ,RYFWCF   ,REBCUA   ,REBCUB   ,REBCUC   ,&
!             &REBCUD   ,REBCUE   ,REBCUF   ,REBCUI   ,REBCUJ   ,&
!             &REBCUG   ,REBCUH   ,RHSAVI   ,RFULIO   ,RFLAA0   ,&
!             &RFLAA1   ,RFLBB0   ,RFLBB1   ,RFLBB2   ,RFLBB3   ,&
!             &RFLCC0   ,RFLCC1   ,RFLCC2   ,RFLCC3   ,RFLDD0   ,&
!             &RFLDD1   ,RFLDD2   ,RFLDD3   ,RFUAA0   ,RFUAA1   ,&
!             &RFUBB0   ,RFUBB1   ,RFUBB2   ,RFUBB3   ,RFUCC0   ,&
!             &RFUCC1   ,RFUCC2   ,RFUCC3   ,RFUETA   ,RASWCA   ,&
!             &RASWCB   ,RASWCC   ,RASWCD   ,RASWCE   ,RASWCF   ,&
!             &RLINLI
! USE YOERDU   , ONLY : RCDAY,  NUAER    ,NTRAER   ,REPLOG   ,REPSC    ,DIFF
! USE YOERDI   , ONLY : REPCLC
! USE YOETHF   , ONLY : RTICE
! USE YOERRTWN , ONLY : NG ,NSPA ,NSPB ,WAVENUM1  ,&
!                       WAVENUM2, DELWAVE, TOTPLNK, TOTPLK16
!
!MESO-NH modules
! USE MODD_PARAMETERS
! USE MODD_RAIN_C2R2_KHKO_PARAM, ONLY : UCREC=>XCREC, UCRER=>XCRER, UFREFFR=>XFREFFR
! USE MODD_RAIN_C2R2_DESCR, ONLY : UAC=>XAC, UAR=>XAR,       &
!                                  ULBEXC=>XLBEXC, ULBEXR=>XLBEXR, &
!                                  URTMIN=>XRTMIN, UCTMIN=>XCTMIN
! USE MODD_ICE_C1R3_PARAM,  ONLY : YFREFFI=>XFREFFI
! USE MODD_ICE_C1R3_DESCR,  ONLY : YLBEXI=>XLBEXI,                      &
!                                  YRTMIN=>XRTMIN, YCTMIN=>XCTMIN
! USE MODD_CST
! USE MODD_CRAD
! 
! USE MODD_PARAM_C2R2,  ONLY : UALPHAC=>XALPHAC,UNUC=>XNUC, &
!                              UALPHAR=>XALPHAR,UNUR=>XNUR !
! USE MODD_PARAM_RAD_n,  ONLY : CAOP                               
! !
! USE MODI_LW
! USE MODI_RRTM_RRTM_140GP
! USE MODI_SW
! !
! ! LIMA
! USE MODD_PARAM_n, ONLY : CCLOUD
! USE MODD_PARAM_LIMA, ONLY : ZRTMIN=>XRTMIN, ZCTMIN=>XCTMIN, &
!                             ZALPHAC=>XALPHAC, ZNUC=>XNUC,   &
!                             ZALPHAR=>XALPHAR, ZNUR=>XNUR
! USE MODD_PARAM_LIMA_WARM, ONLY : ZCREC=>XCREC, ZCRER=>XCRER, ZFREFFR=>XFREFFR, &
!                                  ZAC=>XAC, ZAR=>XAR, ZLBEXC=>XLBEXC, ZLBEXR=>XLBEXR
! USE MODD_PARAM_LIMA_COLD, ONLY : ZFREFFI=>XFREFFI, ZLBEXI=>XLBEXI
#ifdef MNH_ECRAD
USE MODI_RADIATION_SCHEME
USE YOERDU   , ONLY : RCDAY, REPSC
USE YOMCST, ONLY: RG
USE MODD_PARAM_ECRAD_n , ONLY : XCCH4, XCN2O, XCNO2, XCCFC11, XCCFC12, XCCFC22, XCCCL4, & 
                               & XCCNSEA, XCCNLND, LRRTM
#endif                       
USE MODD_RADIATIONS_n , ONLY : NSWB_OLD,NSWB_MNH

USE MODD_GRID_n , ONLY : XLAT, XLON
USE MODD_GRID , ONLY : XLAT0, XLON0
USE MODD_CONF,         ONLY : LCARTESIAN
USE MODD_PARAMETERS , ONLY : XUNDEF
USE MODE_ll , ONLY : GET_INDICE_ll
USE MODD_CST , ONLY : XCPD, XCL, XCI, XCPV, XCPD,XPI

USE PARKIND1 , ONLY : JPIM, JPRB

IMPLICIT NONE
!

!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!  
INTEGER, INTENT(IN) :: KAER !number of aerosol class     
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN) ::PDZ !thickness of the mesh (m)
INTEGER, INTENT(IN) :: KLEV ! number of vertical level for radiation calulation 
INTEGER, INTENT(IN) :: KLON ! number of columns            " 
INTEGER, INTENT(IN) :: KRAD_DIAG ! index for the number of diagnostic fields 
!                                        choice in
CHARACTER (LEN=*), INTENT (IN) :: HEFRADL !cloud water effective radius calculation  
CHARACTER (LEN=*), INTENT (IN) :: HEFRADI !ice water effective radius calculation
CHARACTER (LEN=*), INTENT (IN) :: HOPWSW !cloud water SW optical properties  
CHARACTER (LEN=*), INTENT (IN) :: HOPISW !ice water SW optical properties 
CHARACTER (LEN=*), INTENT (IN) :: HOPWLW !cloud water LW optical properties
CHARACTER (LEN=*), INTENT (IN) :: HOPILW !ice water  LW optical properties
REAL, INTENT(IN)               :: PFUDG  !subgrid cloud inhomogeneity factor
!
REAL(KIND=JPRB), INTENT(INOUT)                     :: PRII0 ! corrected solar constant
REAL, INTENT(IN)                     :: PCCO2 ! CO2 content (Pa/Pa)
REAL(KIND=JPRB), DIMENSION (:,:,:), INTENT (IN) :: PAER  ! aerosol optical thickness
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PALBD ! surface diffuse spectral albedo
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PALBP ! surface direct spectral albedo
REAL(KIND=JPRB), DIMENSION (:), INTENT (IN)     :: PEMIS ! surface emissivity
REAL(KIND=JPRB), DIMENSION (:), INTENT (IN)     :: PEMIW ! surface emissivity in LW window
REAL(KIND=JPRB), DIMENSION (:), INTENT (IN)     :: PLSM  ! land sea mask
REAL(KIND=JPRB), DIMENSION (:), INTENT (IN)     :: PMU0  ! cosine of solar angle 
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: POZON ! ozone content (Pa/Pa)
REAL(KIND=JPRB), DIMENSION (:), INTENT (IN)     :: PTS   ! surfaec temperature
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PT    ! mean layer  temperature (mass point) 
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PAP   ! mean layer  pressure (mass point)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PTH   ! half-level temperature
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PAPH  ! half-level pressure 
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PDP   ! layer pressure thickness
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PQ    ! mean layer specific humidity  (Pa/pa) 
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PQS   ! mean layer saturation spec. humid.
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PQIWC ! mean-layer ice specific water content (kg/kg)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PIWC ! mean-layer ice water content (kg/m3)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PQLWC ! mean-layer liquid specific water content(kg/Kg)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PLWC ! mean-layer liquid water content(kg/m3)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PQRWC ! mean-layer rain specific water content(kg/kg)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PRWC ! mean-layer rain water content(kg/m3)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PCLFR  ! mean-layer cloud fraction
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PCCT_C2R2 ! cloud water concentration (C2R2)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PCRT_C2R2 ! rain water concentration (C2R2)
REAL(KIND=JPRB), DIMENSION (:,:), INTENT (IN)   :: PCIT_C1R3 ! ice crystal concentration (C1R3)
REAL(KIND=JPRB), DIMENSION(:,:,:),INTENT(IN)    :: PPIZA_DST   !Single scattering albedo of dust (wvl dependent)
REAL(KIND=JPRB), DIMENSION(:,:,:),INTENT(IN)    :: PCGA_DST    !Assymetry factor for dust (wvl dependent)
REAL(KIND=JPRB), DIMENSION(:,:,:),INTENT(IN)    :: PTAUREL_DST !Optical depth of dust relative to the one at 550nm
LOGICAL, INTENT (IN)                 :: ODUST  ! flag for dust
!
REAL(KIND=JPRB), DIMENSION (:), INTENT (IN)     :: PLAT ! latitude
REAL(KIND=JPRB), DIMENSION (:), INTENT (IN)     :: PLON  ! longitude
!
! OUTPUTS 
!
REAL, DIMENSION (:,:), INTENT (OUT) :: PDTLW   ! LW temperature tendency
REAL, DIMENSION (:,:), INTENT (OUT) :: PDTSW   ! SW temperature tendency
REAL, DIMENSION (:,:), INTENT (OUT) :: PFLUX_TOP_GND_IRVISNIR ! Top and Ground rad. FLUX.
REAL, DIMENSION (:,:), INTENT (OUT) :: PSFSWDIR ! surface SW direct flux
REAL, DIMENSION (:,:), INTENT (OUT) :: PSFSWDIF ! surface SW diffuse flux
! 
!KRAD_DIAG >=1 --> optional: flux profiles
!
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFLT ! Total LW net flux
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFLS ! Total SW net flux
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFSDWN! Downward SW flux
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFSUP ! Upward SW flux
REAL, DIMENSION (:,:,:), INTENT (OUT):: PFLUX_LW ! LW flux (upward and downward)
!
!KRAD_DIAG >=2 --> optional: clear-sky outputs
!
REAL, DIMENSION (:,:), INTENT (OUT) :: PDTLW_CS   ! LW clear sky temperature tendancy
REAL, DIMENSION (:,:), INTENT (OUT) :: PDTSW_CS   ! SW  clear sky temperature tendancy
REAL, DIMENSION (:,:), INTENT (OUT) :: PFLUX_TOP_GND_IRVISNIR_CS !  Top and
                                                                 !  Ground radiative Clear-sky FLUXes
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFCT ! Clear-sky LW net flux 
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFCS ! Clear-sky SW net flux
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFCDWN  ! Downward SW Clear sky flux 
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFCUP   ! Upward SW Clear sky flux 
REAL, DIMENSION (:,:,:), INTENT (OUT) :: PFLUX_CLW !Clear sky  LW flux (upward and downward)
!
!KRAD_DIAG >=3 --> optional: other macroscpic radiative parameteres
!
REAL, DIMENSION (:), INTENT (OUT) :: PPLAN_ALB_VIS !PLANetary ALBedo in VISible 
REAL, DIMENSION (:), INTENT (OUT) :: PPLAN_ALB_NIR !     "          Near-InfraRed
REAL, DIMENSION (:), INTENT (OUT) :: PPLAN_TRA_VIS !PLANetary TRANsmission in VISible
REAL, DIMENSION (:), INTENT (OUT) :: PPLAN_TRA_NIR !     "          Near-InfraRed
REAL, DIMENSION (:), INTENT (OUT) :: PPLAN_ABS_VIS !PLANetary ABSorption in VISible
REAL, DIMENSION (:), INTENT (OUT) :: PPLAN_ABS_NIR !     "          Near-InfraRed 

!
!KRAD_DIAG >=4 --> optional: more cloud effect radiative parameters 
!
REAL, DIMENSION (:,:), INTENT (OUT) :: PFLWP       ! Liquid water path
REAL, DIMENSION (:,:), INTENT (OUT) :: PFIWP       ! Ice water path
REAL, DIMENSION (:,:), INTENT (OUT) :: PRADLP      ! Cloud water effective radius  
REAL, DIMENSION (:,:), INTENT (OUT) :: PRADIP      ! Cloud ice effective radius
REAL, DIMENSION (:,:), INTENT (OUT) :: PEFCL_LWD   ! effective downward LW nebulosity 
REAL, DIMENSION (:,:), INTENT (OUT) :: PEFCL_LWU   ! effective upward LW nebulosity 
                                                   ! Note: not meaningfull when using RRTM 
REAL, DIMENSION (:,:), INTENT (OUT) :: PEFCL_RRTM ! Effective LW nebulosity ( RRTM case)
REAL, DIMENSION (:,:), INTENT (OUT) :: PCLSW_TOTAL ! Effective SW cloud fraction(mixed phase)
REAL, DIMENSION (:,:,:), INTENT (OUT) :: PTAU_TOTAL !Effective cloud optical thickness
REAL, DIMENSION (:,:,:), INTENT (OUT) :: POMEGA_TOTAL! "   single scattering albedo
REAL, DIMENSION (:,:,:), INTENT (OUT) :: PCG_TOTAL   ! "   asymetry factor 
REAL :: ZLOGREFF  !LOG10of the effective radius 
!  
!
!*       0.2  DECLARATION OF LOCAL VARIABLES
!              -------------------------
!    
! INTEGER       :: IKIDIA, IKFDIA   ! vector indexes in ecmwf code 
! INTEGER       :: IKL, JK, JKL, JKLP1, JL, JNU, JRTM, JSW !loop indices
! INTEGER       :: INDLAY
! 
! REAL          :: ZASYMX, ZDIFFD, ZGI, ZGL, ZGR, ZIWGKG, ZLWGKG,&
!      ZMULTI, ZMSAID, ZMSAIU, ZMSALD, ZMSALU,ZMSARD, ZMSARU, &
!      ZLWFUDG, ZSWFUDG, ZMULTL, ZMULTR, ZOI, ZOL, ZOMGMX, ZOR, &
!      ZRMUZ, ZRWGKG, ZTAUD, ZTAUMX, ZTEMPC, &
!      ZTOI, ZTOL, ZTOR, ZZFIWP, ZZFLWP, ZZFRWP, ZDPOG, ZPODT
     
! REAL          :: ZALND, ZASEA, ZD, ZDEN, ZNTOT, ZNUM, ZRATIO, Z1RADI,&
!      ZBETAI, ZOMGI, ZOMGP, ZFDEL,  ZTCELS, ZFSR,  ZAIWC, ZBIWC,  &
!      ZTBLAY, ZADDPLK,  ZPLANCK, Z1RADL, Z1RADR
! 
! REAL, DIMENSION(KLON) ::  ZTCLEAR, ZDT0, ZEMIS, ZEMIW, &
!      ZFIWP , ZFLWP, ZFRWP, ZIWC,      &
!      ZLWC, ZMU0, ZPSOL, ZVIEW,        &
!      ZBICFU,  ZKICFU1, ZKICFU2,       &
!      ZFSUPN, ZFSUPV,ZFCUPN, ZFCUPV,   &
!      ZFSDNN, ZFSDNV, ZFCDNN, ZFCDNV,  &  
!      ZALFICE, ZGAMICE, ZBICE,         &
!      ZRADIP, ZDESR,ZRES,              &
!      ZRADLP,ZRADRP,ZRADLS,ZRADRS,     &
!      ZTICE , ZEMIT, ZSUDU, ZRHO, ZMR, &
!      ZUVDF, ZPARF
! !cc           , ZRADRD
! !
! REAL, DIMENSION(KLON,NSW)    :: ZALBD , ZALBP , ZDIRFS, ZDIFFS
! REAL, DIMENSION(KLON,KLEV)   :: ZCLFR,  ZCLDLD  , ZCLDLU, ZCLDSW,        &  
!      ZOZON, ZOZ   , ZOZN,    ZTAVE ,  ZDPGCP, &
!      ZCOOLR  , ZCOOLC, ZHEATR  , ZHEATC,      & 
!      ZDFLWT  , ZDFLWC, ZDFSWT  , ZDFSWC
! !
! REAL, DIMENSION(KLON,KLEV+1) :: ZPMB  , ZTL, & 
!      ZFCDWN, ZFCUP, ZFSDWN, ZFSUP, &
!      ZFLT, ZFCT,ZFCS, ZFLS
! !
! REAL, DIMENSION(KLON,NSW,KLEV) :: ZCG ,ZOMEGA, ZTAU
! !
! REAL, DIMENSION(KLON,2,KLEV+1) :: ZFLUX_LW, ZFLUX_CLW      
! !
! REAL, DIMENSION(KLON,KLEV,16)  :: ZTAUCLD
! !
! REAL, DIMENSION(KLON,KAER,KLEV) :: ZAER_SW,ZAER_LW ! Optical aerosol properties
! LOGICAL :: GPROP_OP             !drapeau sur les condition a remplir pour que le
!                                 !calcul des propriétés optiques soit effectué
! !
! REAL, ALLOCATABLE, DIMENSION(:) :: XRTMIN, XCTMIN
! REAL :: XALPHAC,XNUC,XALPHAR,XNUR,XCREC,XCRER,XFREFFR,XAC,XAR,XLBEXC,XLBEXR,XFREFFI,XLBEXI
! !--------------------------------------------------------------

INTEGER       :: IKIDIA, IKFDIA   ! vector indexes in ecmwf code 
INTEGER :: IIB, IJB, IIE, IJE, IIJ, II, JJ, IJ, JI, JL, JK, JKL, IKL,JKLP1 !loop indices

! Old parameters still needed
REAL :: ZIWGKG, ZLWGKG, ZRWGKG 
REAL, DIMENSION(KLON) :: ZFLWP, ZFIWP, ZFRWP
REAL, DIMENSION(KLON,KLEV) :: ZDPGCP, ZDFSWC, ZDFSWT, ZDFLWC, ZDFLWT, ZCLFR
REAL, DIMENSION(KLON,KLEV+1) :: ZFLT, ZFCT, ZFCS, ZFLS

! 
!New local parameters specific to ecrad - Q.L.
REAL(KIND=JPRB), DIMENSION(KLON) :: ZGELAM , ZGEMU
REAL(KIND=JPRB), DIMENSION(KLON) :: ZCCNLND , ZCCNSEA
REAL(KIND=JPRB), DIMENSION(KLON,KLEV) :: ZCCO2, ZCCH4, ZCN2O, ZCNO2, ZCCFC11, ZCCFC12, ZCCFC22, ZCCCL4
INTEGER :: IIDIA , IFDIA
REAL(KIND=JPRB), DIMENSION(KLON,6,KLEV) :: ZAEROSOL_OLD ! initialized to 0 for the time being
INTEGER , PARAMETER :: IAER = 12 ! number of aerosol types in MACC climatology, not used if XAERMACC = 0
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,IAER) :: ZAEROSOL
REAL(KIND=JPRB), DIMENSION(KLON,KLEV) :: ZOZON
REAL(KIND=JPRB), DIMENSION(KLON,KLEV) :: ZQSWC
!---------------------

REAL, DIMENSION(KLON,KLEV) :: ZAP,ZT,ZOZ
REAL, DIMENSION(KLON,KLEV+1) :: ZAPH,ZTH


        
! Outputs of RADIATION_SCHEME
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_SW 
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_LW 
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_SW_CLEAR 
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_LW_CLEAR 
REAL(KIND=JPRB),DIMENSION(KLON,NSWB_MNH) :: ZFLUX_SW_SURF 
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_LW_SURF 
REAL(KIND=JPRB),DIMENSION(KLON,NSWB_MNH) :: ZFLUX_SW_SURF_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_LW_SURF_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON,NSWB_MNH) :: ZFLUX_DIR_SURF
REAL(KIND=JPRB),DIMENSION(KLON,NSWB_MNH) :: ZFLUX_DIR_SURF_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON,NSWB_MNH) :: ZFLUX_DIR_SURF_INTO_SUN
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_UV
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_PAR
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_PAR_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_SW_DN_TOA
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_SW_UP_TOA
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_LW_UP_TOA
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_SW_UP_TOA_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON) :: ZFLUX_LW_UP_TOA_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_SW_DN 
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_LW_DN
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_SW_UP 
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_LW_UP
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_SW_DN_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_LW_DN_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_SW_UP_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZFLUX_LW_UP_CLEAR
REAL(KIND=JPRB),DIMENSION(KLON,KLEV) :: ZRE_LIQ_UM
REAL(KIND=JPRB),DIMENSION(KLON,KLEV) :: ZRE_ICE_UM
REAL(KIND=JPRB),DIMENSION(KLON) :: ZEMIS_OUT   
REAL(KIND=JPRB),DIMENSION(KLON,KLEV+1) :: ZLWDERIVATIVE
REAL(KIND=JPRB),DIMENSION(KLON,NSWB_MNH) :: ZSWDIFFUSEBAND
REAL(KIND=JPRB),DIMENSION(KLON,NSWB_MNH) :: ZSWDIRECTBAND    


! New parameters needed for RADIATION_SCHEME wrt ECMWF_RADIATION_VERS2
! KIDIA, KFDIA (indices of first and last column)
! PCCN_LAND ; PCCN_SEA (set in ini_radiations_ecrad)
! PGELAM ; PGEMU (defined from current latitude and longitude)
! trace gases mixing ratios (climatology in ini_radiations_ecrad)
! PAEROSOL (aerosol mixing ratios from MACC)

! Parameters not entered as inputs to ECMWF_RADIATION_VERS2
! Either defined here, or defined in ini_radiations_ecrad

!

#ifdef MNH_ECRAD
IKIDIA = 1
IKFDIA = KLON
ZGELAM = PLON! longitude (in radians)
ZGEMU = sin(PLAT) ! sin (lat) 
ZCCNLND(:) = XCCNLND
ZCCNSEA(:) = XCCNSEA 
ZCCO2(:,:) = PCCO2
ZCCH4(:,:) = XCCH4
ZCN2O(:,:) = XCN2O
ZCNO2(:,:) = XCNO2
ZCCFC11(:,:) = XCCFC11
ZCCFC12(:,:) = XCCFC12
ZCCFC22(:,:) = XCCFC22
ZCCCL4 = XCCCL4
! ZOZON(:,:) = POZON (:,:) * PDP (:,:) * 46.6968 / RG ! OZONE field = (kg kg) * DP IN RADIATION_SCHEME
ZOZON(:,:) = POZON (:,:) * PDP (:,:) ! OZONE field = (kg kg) * DP IN RADIATION_SCHEME

! ZAEROSOL_OLD(KLON,6,KLEV)
! PAER(KLON,KLEV,6)
! Keep the original formulation for aerosols
DO JJ=1,SIZE(PAER,3)
    DO JI=1,KLEV
        ZAEROSOL_OLD(:,JJ,JI) =  PAER(:,JI,JJ)
    END DO
END DO

! ! Invert all arrays that enters as inputs from MNH_LIC! PAP, PT, PAPH, PTH
! 
! DO JK = 1 , KLEV
!   JKL = KLEV+ 1 - JK
!   DO JL = IKIDIA,IKFDIA
!     ZAPH(JL,JK)=PAPH(JL,JKL+1)
!     ZAP(JL,JK)=PAP(JL,JKL)
!     ZOZ(JL,JK)  = ZOZON(JL,JKL) * 46.6968 / RG
!     ZTH(JL,JK)= PTH(JL,JKL+1)
!     ZT(JL,JK)=PT(JL,JKL)
!   ENDDO
! ENDDO
! !
! DO JL=IKIDIA,IKFDIA
!   ZTH(JL,KLEV+1)= PTH(JL,1)
!   ZAPH(JL,KLEV+1) = PAPH(JL,1) ! Conversion needed? No because all pressures in Pa
! ENDDO

! CAMS climatology


! aerosol mixing ratio from MACC reanalyses if XAERMACC > 0
ZAEROSOL(:,:,:) = 0

ZQSWC(:,:) = 0 ! snow mixing ratio set to 0 for the time being

CALL RADIATION_SCHEME  (IKIDIA, IKFDIA,KLON, KLEV, IAER, &
        PRII0,PMU0, PTS, PALBD, PALBP, &
        PEMIS, PEMIW,ZCCNLND, ZCCNSEA,ZGELAM, ZGEMU, PLSM, &
        PAP, PT, PAPH, PTH, &
        PQ, ZCCO2, ZCCH4, ZCN2O, ZCNO2, ZCCFC11, ZCCFC12, ZCCFC22, ZCCCL4, ZOZON, &
        PCLFR, PQLWC, PQIWC, PQRWC, ZQSWC, &
        ZAEROSOL_OLD, ZAEROSOL, &
        ZFLUX_SW, ZFLUX_LW, ZFLUX_SW_CLEAR, ZFLUX_LW_CLEAR, &
        ZFLUX_SW_SURF, ZFLUX_LW_SURF, ZFLUX_SW_SURF_CLEAR, ZFLUX_LW_SURF_CLEAR, &
        ZFLUX_DIR_SURF, ZFLUX_DIR_SURF_CLEAR, ZFLUX_DIR_SURF_INTO_SUN, &
        ZFLUX_UV, ZFLUX_PAR, ZFLUX_PAR_CLEAR, &
        ZFLUX_SW_DN_TOA, ZFLUX_SW_UP_TOA, ZFLUX_LW_UP_TOA, &
        ZFLUX_SW_UP_TOA_CLEAR, ZFLUX_LW_UP_TOA_CLEAR, &
        ZFLUX_SW_DN, ZFLUX_LW_DN, ZFLUX_SW_UP, ZFLUX_LW_UP, &
        ZFLUX_SW_DN_CLEAR, ZFLUX_LW_DN_CLEAR, ZFLUX_SW_UP_CLEAR, ZFLUX_LW_UP_CLEAR, &
        ZRE_LIQ_UM, ZRE_ICE_UM, &
        ZEMIS_OUT, ZLWDERIVATIVE, &
        ZSWDIFFUSEBAND, ZSWDIRECTBAND                          )
     

! PRINT*,"ZFLUX_DIR_SURF(1,:)",SUM(ZFLUX_DIR_SURF(1,:))
! PRINT*,"ZFLUX_LW_SURF(1)",ZFLUX_LW_SURF(1)
! PAUSE
        
! Once The flux profiles have been computed, compute diagnostics consistent with ECMWF_RADIATION_VERS2
! Requires conversion from RADIATION_SCHEME outputs into usual outputs consistent with the rest of MESONH

!*         5.     FLUXES, RADIATIVE TENDENCIES, DIAGNOSTICS
!               ------------------------------------------------
! note : important
! After the RADIATION_SCHEME routine, radiative fluxes array are consistent with
! ECRAD vertical grid : i.e FLUX (:,1) is the TOA flux. It must be reversed for RADIATIONS
!
!
!          5.0 net fluxes
!
! Reversing fluxes
DO JK = 1 , KLEV+1
  JKL = KLEV+1 -JK +1  
  DO JL = IKIDIA,IKFDIA
    ZFLS(JL,JK) = ZFLUX_SW(JL,JKL)
    ZFLT(JL,JK) = ZFLUX_LW(JL,JKL)
    ZFCS(JL,JK) = ZFLUX_SW_CLEAR(JL,JKL)
    ZFCT(JL,JK) = ZFLUX_LW_CLEAR(JL,JKL)
  ENDDO
ENDDO

!          5.1 Radiative tendencies in (K/DAY)
!
! calculation of DP term (grid orientation consistent with net fluxes,
! ( see note above) 
!computation of CPH (RCDAY is computed with CPD)
DO JKL = 1 , KLEV
  DO JL = IKIDIA,IKFDIA
    JK = KLEV+1  - JKL
    ZDPGCP(JL,JK) = RCDAY/ PDP(JL, JKL)*XCPD/(XCPD+XCPV*PQ(JL, JKL)/(1.-PQ(JL, JKL)) &
                          +XCL*(PQLWC(JL, JKL)+PQRWC(JL, JKL))+XCI*PQIWC(JL, JKL))
  END DO
END DO

! Heating rates prop. to. derivative of net fluxes
DO JK=1,KLEV
  DO JL=IKIDIA,IKFDIA       
    ZDFLWT(JL,JK) = ZFLT(JL,JK+1)-ZFLT(JL,JK)
    ZDFLWC(JL,JK) = ZFCT(JL,JK+1)-ZFCT(JL,JK)
    ZDFSWT(JL,JK) = ZFLS(JL,JK+1)-ZFLS(JL,JK)
    ZDFSWC(JL,JK) = ZFCS(JL,JK+1)-ZFCS(JL,JK)
!
    PDTLW(JL,JK) = ZDFLWT(JL,JK) * ZDPGCP(JL,JK) 
    PDTSW(JL,JK) = ZDFSWT(JL,JK) * ZDPGCP(JL,JK)
    !clear-sky
    PDTLW_CS(JL,JK) = ZDFLWC(JL,JK) * ZDPGCP(JL,JK)
    PDTSW_CS(JL,JK) = ZDFSWC(JL,JK) * ZDPGCP(JL,JK)
  END DO
END DO
!
!         5.2  Top and Ground fluxes (IR-VIS-NIR) - NB no distinction between VIS and NIR so far - Check signs
!

PFLUX_TOP_GND_IRVISNIR(:, 1) =   ZFLUX_LW_UP_TOA (:)   ! TOA LW UP
PFLUX_TOP_GND_IRVISNIR(:, 2) =   ZFLUX_SW_UP_TOA(:)    ! TOA VIS
PFLUX_TOP_GND_IRVISNIR(:, 3) =   ZFLUX_SW_UP_TOA(:)    ! TOA NIR
PFLUX_TOP_GND_IRVISNIR(:, 4) =   ZFLUX_LW_SURF (:)     ! Surface LW down
PFLUX_TOP_GND_IRVISNIR(:, 5) =   ZFLUX_SW_SURF(:,1)        ! Surface VIS down
PFLUX_TOP_GND_IRVISNIR(:, 6) =   ZFLUX_SW_SURF(:,1)        ! Surface NIR down

!Clear-sky
!
PFLUX_TOP_GND_IRVISNIR_CS(:, 1) =   ZFLUX_LW_UP_TOA_CLEAR (:)
PFLUX_TOP_GND_IRVISNIR_CS(:, 2) =   ZFLUX_SW_UP_TOA_CLEAR(:)
PFLUX_TOP_GND_IRVISNIR_CS(:, 3) =   ZFLUX_SW_UP_TOA_CLEAR(:)
PFLUX_TOP_GND_IRVISNIR_CS(:, 4) =   ZFLUX_LW_SURF_CLEAR (:) 
PFLUX_TOP_GND_IRVISNIR_CS(:, 5) =   ZFLUX_SW_SURF_CLEAR(:,1) 
PFLUX_TOP_GND_IRVISNIR_CS(:, 6) =   ZFLUX_SW_SURF_CLEAR(:,1)
!

PSFSWDIR = ZFLUX_DIR_SURF
PSFSWDIF = ZFLUX_SW_SURF - ZFLUX_DIR_SURF

!             5.2  Radiation Diagnostics 
! Total fluxes
DO JK=1,KLEV+1
  JKL = KLEV+1 -JK +1  
  DO JL=IKIDIA,IKFDIA   
    PFSDWN(JL,JK) = ZFLUX_SW_DN(JL,JKL)
    PFSUP (JL,JK) = ZFLUX_SW_UP (JL,JKL)
    PFLUX_LW(JL,2,JK) = -ZFLUX_LW_DN(JL,JKL) ! to be consistent with ECMWF_RADIATION_VERS2 where fluxes are always upwards
    PFLUX_LW(JL,1,JK) = ZFLUX_LW_UP(JL,JKL)
  END DO
END DO   

IF (KRAD_DIAG >= 1) THEN ! net fluxes
    PFLT(:,:) = ZFLT(:,:)
    PFLS(:,:) = ZFLS(:,:)
END IF
!
IF (KRAD_DIAG >=2) THEN ! clear-sky fluxes
    DO JK=1,KLEV+1
        JKL = KLEV+1 -JK +1 
        DO JL=IKIDIA,IKFDIA  
            PFCDWN(JL,JK) = ZFLUX_SW_DN_CLEAR(JL,JKL) 
            PFCUP (JL,JK) = ZFLUX_SW_UP_CLEAR (JL,JKL)
            PFLUX_CLW(JL,2,JK) = ZFLUX_LW_DN_CLEAR(JL,JKL)
            PFLUX_CLW(JL,1,JK) = ZFLUX_LW_UP_CLEAR(JL,JKL)
        END DO
    END DO  
END IF

IF (KRAD_DIAG >= 2) THEN ! net fluxes
    PFCS(:,:) = ZFCS(:,:)
    PFCT(:,:) = ZFCT(:,:)
END IF


! Unchanged
IF (KRAD_DIAG >= 3) THEN 
  ! provisoire 
  PPLAN_ALB_VIS(:) = XUNDEF
  PPLAN_ALB_NIR(:) = XUNDEF
  PPLAN_TRA_VIS(:) = XUNDEF 
  PPLAN_TRA_NIR(:) = XUNDEF
  PPLAN_ABS_VIS(:) = XUNDEF 
  PPLAN_ABS_NIR(:) = XUNDEF

END IF

! QL
IF (KRAD_DIAG >=4) THEN 
  IF (.NOT.LRRTM) THEN
    PEFCL_LWD  (:,:) = XUNDEF ! effective cloudiness LW up and down
    PEFCL_LWU    (:,:) = XUNDEF
    PEFCL_RRTM (:,:) = XUNDEF
  ELSE
    PEFCL_LWD  (:,:) = XUNDEF
    PEFCL_LWU    (:,:) = XUNDEF
    PEFCL_RRTM (:,:)  = XUNDEF
  END IF
   
  PCLSW_TOTAL (:,:)   = XUNDEF
  PTAU_TOTAL  (:,:,:) = XUNDEF
  POMEGA_TOTAL(:,:,:) = XUNDEF
  PCG_TOTAL   (:,:,:) = XUNDEF  ! to be updfated if necessay
END IF
!

!    2.2   cloud ice and liquid content and path - QL
!
DO JL = IKIDIA,IKFDIA
    ZFLWP(JL) = 0
    ZFIWP(JL) = 0
    ZFRWP(JL) = 0
    DO JK = 1 , KLEV
        IKL = KLEV + 1 - JK
        ZCLFR(JL,IKL)=MAX(REPSC,MIN(PCLFR(JL,IKL),1-REPSC))

! --- Liquid Water Content (g.m-3) and Liquid Water Path (g.m-2)
!
        ZLWGKG = MAX(PLWC(JL,IKL)*1000.,0.)
        ZIWGKG = MAX(PIWC(JL,IKL)*1000.,0.)
        ZRWGKG = 0.0   

    IF (ZCLFR(JL,IKL) > 1E-2 ) THEN
        ZLWGKG = ZLWGKG/ZCLFR(JL,IKL) ! contribution of level to water path
        ZIWGKG = ZIWGKG/ZCLFR(JL,IKL)
      
    ELSE
        ZLWGKG=0.
        ZIWGKG=0.
    ENDIF

    ZFLWP(JL) = ZFLWP(JL) + ZLWGKG*PDZ(JL,IKL)
    ZFIWP(JL) = ZFIWP(JL) + ZIWGKG*PDZ(JL,IKL)
    ZFRWP(JL) = ZFRWP(JL) + ZRWGKG*PDZ(JL,IKL)
    
    END DO   
END DO

!   Diagnostic for water path  
IF (KRAD_DIAG >= 4 ) THEN
    DO JL = IKIDIA,IKFDIA
        PFLWP(JL,:) = ZFLWP(JL)
        PFIWP(JL,:) = ZFIWP(JL)
    END DO
    PRADIP(:,:) = ZRE_ICE_UM(:,:)
END IF

PRADLP(:,:) = ZRE_LIQ_UM(:,:)


#else
PRINT *, "ECRAD LIBRARY NOT AVAILABLE = ###ECRAD_INTERFACE###"
#endif
    
END SUBROUTINE ECRAD_INTERFACE    

END MODULE MODI_ECRAD_INTERFACE
