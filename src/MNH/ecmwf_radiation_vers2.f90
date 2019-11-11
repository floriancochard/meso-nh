!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!##############################################################
!OPTION! -Ni
SUBROUTINE ECMWF_RADIATION_VERS2 ( KLON,KLEV,KRAD_DIAG, KAER, &
     PDZ,HEFRADL, HEFRADI, HOPWSW, HOPISW, HOPWLW, HOPILW,PFUDG,  &
     PRII0, PAER , PALBD , PALBP, PAPH , PAP,                 &
     PCCO2, PCLFR , PDP  , PEMIS, PEMIW , PLSM , PMU0, POZON, &
     PQ   , PQIWC, PIWC, PQLWC, PLWC,PQS  , PQRWC, PRWC,              &
     PTH  , PT    , PTS, PCCT_C2R2, PCRT_C2R2, PCIT_C1R3,     &
     PFCT  , PFLT , PFCS , PFLS  ,                            &  
     PDTLW, PDTSW ,PFLUX_TOP_GND_IRVISNIR,                    & 
     PSFSWDIR, PSFSWDIF,                                      &
     PFSDWN, PFSUP, PFLUX_LW ,                                &
     PDTLW_CS, PDTSW_CS ,PFLUX_TOP_GND_IRVISNIR_CS,           &
     PFCDWN, PFCUP, PFLUX_CLW,                                & 
     PPLAN_ALB_VIS, PPLAN_ALB_NIR, PPLAN_TRA_VIS, PPLAN_TRA_NIR,&
     PPLAN_ABS_VIS, PPLAN_ABS_NIR, PEFCL_LWD, PEFCL_LWU,      &
     PFLWP,PFIWP,PRADLP,PRADIP, PEFCL_RRTM, PCLSW_TOTAL,      & 
     PTAU_TOTAL, POMEGA_TOTAL, PCG_TOTAL,                     &
     ODUST,PPIZA_DST,PCGA_DST,PTAUREL_DST                     )
!##############################################################
!  
!**** *ECMWF_RADIATION_VERS2* - INTERFACE TO ECMWF LW AND SW RADIATION SCHEMES
!
!     PURPOSE.
!     --------
!     Interface routine to call ECMWF SW and LW radiation schemes in the MesoNH 
!     framework. It replaces the ECMWF_RADIATION_VERS1 (JP. Pinty) of precedent cycles.
!
!     METHOD.
!     -------
!     This routine is derived from the inital RADLSW driver coming from ECMWF model
!     (Morcrette et al. ). The first part of the routine deals with the set up of 
!      quantities used in SW and LW code and the calculation of cloud optical properties.
!     Then SW and LW scheme are called in function of configuration options (eg RRTM or LW). 
!     At last, pertinent variable (fluxes and radiative tendencies) together with a number 
!     of diagnostics are calculated and passed up to radiations.      
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
!        J.-J. MORCRETTE  *ECMWF* for initial RADLSW routine 
!        F. SOLMON (Mar 2002) for adaptation to MesoNH
!
!     MODIFICATIONS.
!     --------------        
!         P. Jabouille 05/05/03 corrective factor CPD/CPH
!         Y. Seity (Aou 2003) SW spectral 6 bands + dir./sca. separation 
!                             for surface downward fluxes used by the surface
!                             scheme.
!         J. Escobar, JP Pinty 09/11/04 corrective ZRADLP<0. for HEFRADL="PRES"
!         I. Sandu, P Tulet 09/04/05 climatologic aerosol for SSA
!         A. Grini          09/04/05 dust direct effect
!         V.Puygrenier 07/2009 Correction on ice effective radius
!         B. Aouizerats 09/2010 Explicit aerosol optical properties computation
!         G.Delautier 9/2014: remplace MODD_RAIN_C2R2_PARAM par MODD_RAIN_C2R2_KHKO_PARAM
!         M.Mazoyer 2016 :  limit of 100 microns for effective radius 
!         B.VIE 2016 : LIMA
!         J.Escobar 30/03/2017  : Management of compilation of ECMWF_RAD in REAL*8 with MNH_REAL=R4
!!        Q.Libois  02/2018 : ECRAD
!-----------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!ECMWF radiation scheme specific modules 
!
USE PARKIND1 , ONLY : JPRB
USE OYOMCST   , ONLY : RG ,RD ,RTT ,RPI
USE OYOERAD   , ONLY : NMODE, NSW ,LRRTM ,LINHOM ,LRADIP, LRADLP 
USE YOELW    , ONLY : NSIL     ,NTRA     ,NUA      ,TSTAND   ,XP
USE OYOESW    , ONLY : RYFWCA   ,RYFWCB   ,RYFWCC   ,RYFWCD   ,&
            &RYFWCE   ,RYFWCF   ,REBCUA   ,REBCUB   ,REBCUC   ,&
            &REBCUD   ,REBCUE   ,REBCUF   ,REBCUI   ,REBCUJ   ,&
            &REBCUG   ,REBCUH   ,RHSAVI   ,RFULIO   ,RFLAA0   ,&
            &RFLAA1   ,RFLBB0   ,RFLBB1   ,RFLBB2   ,RFLBB3   ,&
            &RFLCC0   ,RFLCC1   ,RFLCC2   ,RFLCC3   ,RFLDD0   ,&
            &RFLDD1   ,RFLDD2   ,RFLDD3   ,RFUAA0   ,RFUAA1   ,&
            &RFUBB0   ,RFUBB1   ,RFUBB2   ,RFUBB3   ,RFUCC0   ,&
            &RFUCC1   ,RFUCC2   ,RFUCC3   ,RFUETA   ,RASWCA   ,&
            &RASWCB   ,RASWCC   ,RASWCD   ,RASWCE   ,RASWCF   ,&
            &RLINLI
USE YOERDU   , ONLY : RCDAY,  NUAER    ,NTRAER   ,REPLOG   ,REPSC    ,DIFF
USE OYOERDI   , ONLY : REPCLC
USE YOETHF   , ONLY : RTICE
USE OYOERRTWN , ONLY : NG ,NSPA ,NSPB ,WAVENUM1  ,&
                      WAVENUM2, DELWAVE, TOTPLNK, TOTPLK16
!
!MESO-NH modules
USE MODD_PARAMETERS
USE MODD_RAIN_C2R2_KHKO_PARAM, ONLY : UCREC=>XCREC, UCRER=>XCRER, UFREFFR=>XFREFFR
USE MODD_RAIN_C2R2_DESCR, ONLY : UAC=>XAC, UAR=>XAR,       &
                                 ULBEXC=>XLBEXC, ULBEXR=>XLBEXR, &
                                 URTMIN=>XRTMIN, UCTMIN=>XCTMIN
USE MODD_ICE_C1R3_PARAM,  ONLY : YFREFFI=>XFREFFI
USE MODD_ICE_C1R3_DESCR,  ONLY : YLBEXI=>XLBEXI,                      &
                                 YRTMIN=>XRTMIN, YCTMIN=>XCTMIN
USE MODD_CST
USE MODD_CRAD

USE MODD_PARAM_C2R2,  ONLY : UALPHAC=>XALPHAC,UNUC=>XNUC, &
                             UALPHAR=>XALPHAR,UNUR=>XNUR !
USE MODD_PARAM_RAD_n,  ONLY : CAOP                               
!
USE MODI_LW
!USE MODI_RRTM_RRTM_140GP
USE MODI_SW
!
! LIMA
USE MODD_PARAM_n, ONLY : CCLOUD
USE MODD_PARAM_LIMA, ONLY : ZRTMIN=>XRTMIN, ZCTMIN=>XCTMIN, &
                            ZALPHAC=>XALPHAC, ZNUC=>XNUC,   &
                            ZALPHAR=>XALPHAR, ZNUR=>XNUR
USE MODD_PARAM_LIMA_WARM, ONLY : ZCREC=>XCREC, ZCRER=>XCRER, ZFREFFR=>XFREFFR, &
                                 ZAC=>XAC, ZAR=>XAR, ZLBEXC=>XLBEXC, ZLBEXR=>XLBEXR
USE MODD_PARAM_LIMA_COLD, ONLY : ZFREFFI=>XFREFFI, ZLBEXI=>XLBEXI
!
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
REAL(KIND=JPRB), INTENT(INOUT)       :: PRII0 ! corrected solar constant
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
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFCT ! Total LW net flux
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFLT ! Total SW net flux
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
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFCS ! Clear-sky LW net flux 
REAL, DIMENSION (:,:), INTENT (OUT)  :: PFLS ! Clear-sky SW net flux
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
REAL::ZLOGREFF  !LOG10of the effective radius 
!  
!
!*       0.2  DECLARATION OF LOCAL VARIABLES
!              -------------------------
!    
INTEGER       :: IKIDIA, IKFDIA   ! vector indexes in ecmwf code 
INTEGER       :: IKL, JK, JKL, JKLP1, JL, JNU, JRTM, JSW !loop indices
INTEGER       :: INDLAY

REAL          :: ZASYMX, ZDIFFD, ZGI, ZGL, ZGR, ZIWGKG, ZLWGKG,&
     ZMULTI, ZMSAID, ZMSAIU, ZMSALD, ZMSALU,ZMSARD, ZMSARU, &
     ZLWFUDG, ZSWFUDG, ZMULTL, ZMULTR, ZOI, ZOL, ZOMGMX, ZOR, &
     ZRMUZ, ZRWGKG, ZTAUD, ZTAUMX, ZTEMPC, &
     ZTOI, ZTOL, ZTOR, ZZFIWP, ZZFLWP, ZZFRWP, ZDPOG, ZPODT
REAL          :: ZALND, ZASEA, ZD, ZDEN, ZNTOT, ZNUM, ZRATIO, Z1RADI,&
     ZBETAI, ZOMGI, ZOMGP, ZFDEL,  ZTCELS, ZFSR,  ZAIWC, ZBIWC,  &
     ZTBLAY, ZADDPLK,  ZPLANCK, Z1RADL, Z1RADR

REAL(KIND=JPRB), DIMENSION(KLON) ::  ZTCLEAR, ZDT0, ZEMIS, ZEMIW, &
     ZFIWP , ZFLWP, ZFRWP, ZIWC,      &
     ZLWC, ZMU0, ZPSOL, ZVIEW,        &
     ZBICFU,  ZKICFU1, ZKICFU2,       &
     ZFSUPN, ZFSUPV,ZFCUPN, ZFCUPV,   &
     ZFSDNN, ZFSDNV, ZFCDNN, ZFCDNV,  &  
     ZALFICE, ZGAMICE, ZBICE,         &
     ZRADIP, ZDESR,ZRES,              &
     ZRADLP,ZRADRP,ZRADLS,ZRADRS,     &
     ZTICE , ZEMIT, ZSUDU, ZRHO, ZMR, &
     ZUVDF, ZPARF
!cc           , ZRADRD
!
REAL(KIND=JPRB), DIMENSION(KLON,NSW)    :: ZALBD , ZALBP , ZDIRFS, ZDIFFS
REAL(KIND=JPRB), DIMENSION(KLON,KLEV)   :: ZCLFR,  ZCLDLD  , ZCLDLU, ZCLDSW,        &  
     ZOZON, ZOZ   , ZOZN,    ZTAVE ,  ZDPGCP, &
     ZCOOLR  , ZCOOLC, ZHEATR  , ZHEATC,      & 
     ZDFLWT  , ZDFLWC, ZDFSWT  , ZDFSWC
!
REAL(KIND=JPRB), DIMENSION(KLON,KLEV+1) :: ZPMB  , ZTL, & 
     ZFCDWN, ZFCUP, ZFSDWN, ZFSUP, &
     ZFLT, ZFCT,ZFCS, ZFLS
!
REAL(KIND=JPRB), DIMENSION(KLON,NSW,KLEV) :: ZCG ,ZOMEGA, ZTAU
!
REAL(KIND=JPRB), DIMENSION(KLON,2,KLEV+1) :: ZFLUX_LW, ZFLUX_CLW      
!
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,16)  :: ZTAUCLD
!
REAL(KIND=JPRB), DIMENSION(KLON,KAER,KLEV) :: ZAER_SW,ZAER_LW ! Optical aerosol properties
LOGICAL :: GPROP_OP             !drapeau sur les condition a remplir pour que le
                                !calcul des propriétés optiques soit effectué
!
REAL, ALLOCATABLE, DIMENSION(:) :: XRTMIN, XCTMIN
REAL :: XALPHAC,XNUC,XALPHAR,XNUR,XCREC,XCRER,XFREFFR,XAC,XAR,XLBEXC,XLBEXR,XFREFFI,XLBEXI
!
REAL(KIND=JPRB) :: ZCCO2_RAD
!--------------------------------------------------------------
ZCCO2_RAD = PCCO2
!
!          0. LIMA
IF ( CCLOUD == "LIMA" ) THEN
   ALLOCATE(XRTMIN(SIZE(ZRTMIN,1)))
   ALLOCATE(XCTMIN(SIZE(ZCTMIN,1)))
   XRTMIN(:)=ZRTMIN(:)
   XCTMIN(:)=ZCTMIN(:)

   XALPHAC = ZALPHAC
   XNUC    = ZNUC
   XALPHAR = ZALPHAR
   XNUR    = ZNUR
   XCREC   = ZCREC
   XCRER   = ZCRER
   XFREFFR = ZFREFFR
   XAC     = ZAC
   XAR     = ZAR
   XLBEXC  = ZLBEXC
   XLBEXR  = ZLBEXR
   XFREFFI = ZFREFFI
   XLBEXI  = ZLBEXI

ELSE IF (CCLOUD == "C3R5") THEN
   ALLOCATE(XRTMIN(SIZE(YRTMIN,1)))
   ALLOCATE(XCTMIN(SIZE(YCTMIN,1)))
   XRTMIN(:)=YRTMIN(:)
   XCTMIN(:)=YCTMIN(:)

   XALPHAC = UALPHAC
   XNUC    = UNUC
   XALPHAR = UALPHAR
   XNUR    = UNUR
   XCREC   = UCREC
   XCRER   = UCRER
   XFREFFR = UFREFFR
   XAC     = UAC
   XAR     = UAR
   XLBEXC  = ULBEXC
   XLBEXR  = ULBEXR
   XFREFFI = YFREFFI
   XLBEXI  = YLBEXI

ELSE IF (CCLOUD == "C2R2" .OR. CCLOUD == "KHKO") THEN
   ALLOCATE(XRTMIN(SIZE(URTMIN,1)))
   ALLOCATE(XCTMIN(SIZE(UCTMIN,1)))
   XRTMIN(:)=URTMIN(:)
   XCTMIN(:)=UCTMIN(:)

   XALPHAC = UALPHAC
   XNUC    = UNUC
   XALPHAR = UALPHAR
   XNUR    = UNUR
   XCREC   = UCREC
   XCRER   = UCRER
   XFREFFR = UFREFFR
   XAC     = UAC
   XAR     = UAR
   XLBEXC  = ULBEXC
   XLBEXR  = ULBEXR
   XFREFFI = YFREFFI
   XLBEXI  = YLBEXI

END IF
!
!*         1.     SET-UP INPUT QUANTITIES FOR RADIATION
!                 -------------------------------------
!  
! Vector characteristics in ecmwf code 
!
IKIDIA=1
IKFDIA= KLON
!
DO JL = IKIDIA,IKFDIA
  ZFCUP(JL,KLEV+1) = 0.
  ZFCDWN(JL,KLEV+1) = REPLOG
  ZFSUP(JL,KLEV+1) = 0.
  ZFSDWN(JL,KLEV+1) = REPLOG
  ZFLUX_LW(JL,1,KLEV+1) = 0.
  ZFLUX_LW(JL,2,KLEV+1) = 0.
  ZFLUX_CLW(JL,1,KLEV+1) = 0.
  ZFLUX_CLW(JL,2,KLEV+1) = 0.
  ZFSDNN(JL) = 0.
  ZFSDNV(JL) = 0.
  ZFCDNN(JL) = 0.
  ZFCDNV(JL) = 0.
  ZFSUPN(JL) = 0.
  ZFSUPV(JL) = 0.
  ZFCUPN(JL) = 0.
  ZFCUPV(JL) = 0.
  ZSUDU (JL) = 0.
  ZRADLP(JL)=0
  ZRADRP(JL)=0
  ZRADIP(JL)=0
  ZPSOL(JL) = PAPH(JL,KLEV+1)
  ZPMB(JL,1) = ZPSOL(JL) / 100.
  ZDT0(JL) = PTS(JL) - PTH(JL,KLEV+1) 
  ZDIRFS(JL,:)=0.
  ZDIFFS(JL,:)=0.
ENDDO
!
!*         1.1    INITIALIZE VARIOUS FIELDS
!                 -------------------------
!
DO JSW=1,NSW
  DO JL = IKIDIA,IKFDIA
    ZALBD(JL,JSW)=PALBD(JL,JSW)
    ZALBP(JL,JSW)=PALBP(JL,JSW)
  ENDDO
ENDDO
DO JL = IKIDIA,IKFDIA
  ZEMIS(JL)  =PEMIS(JL)
  ZEMIW(JL)  =PEMIW(JL)
ENDDO
!
WHERE (PMU0(:) >= 0.) 
  ZMU0(:) = PMU0(:)
ELSEWHERE
  ZMU0(:) = 0.
END WHERE
!
!Care : OZONE field = (concentration) * DP, in initial ECMWF driver routines
!
ZOZON(:,:) = RG / 46.6968 * POZON (:,:) * PDP (:,:) ! to be consistent with the actual input of RADLSW
!
! Aerosol field
!

DO JL = IKIDIA,IKFDIA
 ZAER_SW(JL,:,:) = TRANSPOSE (PAER(JL,:,:))
 ZAER_LW(JL,:,:) = TRANSPOSE (PAER(JL,:,:))
END DO
!
IF (CAOP=='EXPL') THEN
 ZAER_LW(:,3,:) = 0.0
ENDIF

DO JK = 1 , KLEV
  JKL = KLEV+ 1 - JK
  JKLP1 = JKL + 1
  DO JL = IKIDIA,IKFDIA
    ZPMB(JL,JK+1)=PAPH(JL,JKL)/100.
    ZOZ(JL,JK)   = ZOZON(JL,JKL) * 46.6968 / RG
    ZFCUP(JL,JK) = 0.
    ZFCDWN(JL,JK) = 0.
    ZFSUP(JL,JK) = 0.
    ZFSDWN(JL,JK) = 0.
    ZFLUX_LW(JL,1,JK) = 0.
    ZFLUX_LW(JL,2,JK) = 0.
    ZFLUX_CLW(JL,1,JK) = 0.
    ZFLUX_CLW(JL,2,JK) = 0.
  ENDDO
ENDDO
!
DO JK=1,KLEV
  JKL=KLEV+1-JK
  JKLP1=JKL+1
  DO JL=IKIDIA,IKFDIA
    ZTL(JL,JK)=PTH(JL,JKLP1)
    ZTAVE(JL,JK)=PT(JL,JKL)
  ENDDO
ENDDO
DO JL=IKIDIA,IKFDIA
  ZTL(JL,KLEV+1)= PTH(JL,1)
  ZPMB(JL,KLEV+1) = PAPH(JL,1)/100.
ENDDO
!
!     ------------------------------------------------------------------
!
!*         2. CLOUD AND AEROSOL PARAMETERS
!             ----------------------------
!
!          2.0  cloud subgrid inhomogenity factor (Tiedke 1995)
!
! Care : cloud inhomogeneity factor fixed to 1 for MNH 
!        "small" grid-cell conf.
!
ZSWFUDG= PFUDG
ZLWFUDG= PFUDG
!
!          2.1  initialize optical properties to clear sky values
!
DO JK = 1 , KLEV
  IKL = KLEV + 1 - JK
  DO JSW = 1,NSW
    DO JL = IKIDIA,IKFDIA
      ZTAU(JL,JSW,JK)  = 0.
      ZOMEGA(JL,JSW,JK)= 1.
      ZCG(JL,JSW,JK)   = 0.
    ENDDO
  ENDDO
  DO JL = IKIDIA,IKFDIA
    ZCLDSW(JL,JK)  = 0.
    ZCLDLD(JL,JK)  = 0.
    ZCLDLU(JL,JK)  = 0.
    ZFLWP (JL) = 0.
    ZFRWP (JL) = 0.
    ZFIWP (JL) = 0.
  ENDDO
!
!          2.2   cloud ice and liquid content and path
!
  DO JL = IKIDIA,IKFDIA
    ZCLFR(JL,IKL)=MAX(REPSC,MIN(PCLFR(JL,IKL),1-REPSC))


!
! --- Liquid Water Content (g.m-3) and Liquid Water Path (g.m-2)
!
    ZLWGKG=MAX(PLWC(JL,IKL)*1000.,0.)
    ZRWGKG=MAX(PRWC(JL,IKL)*1000.,0.)
    ZIWGKG=MAX(PIWC(JL,IKL)*1000.,0.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Parametrisation des propriete optiques de la pluie non implementee
! Modification to take into account rain below cloud for C2R2 scheme
!    IF (HEFRADL=='C2R2') THEN
!      ZRWGKG=MAX(PRWC(JL,IKL)*1000.,0.)
!      IF (ZRWGKG.ne.0.0) ZCLFR(JL,IKL)=1.0
!    ELSE
!      ZRWGKG=0.0   
!    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (ZCLFR(JL,IKL) > 1E-2 ) THEN
      ZLWGKG=ZLWGKG/ZCLFR(JL,IKL)
      ZRWGKG=ZRWGKG/ZCLFR(JL,IKL)      
      ZIWGKG=ZIWGKG/ZCLFR(JL,IKL)
    ELSE
      ZLWGKG=0.
      ZRWGKG=0.
      ZIWGKG=0.
    ENDIF

!
    ZFLWP(JL)= ZLWGKG*PDZ(JL,IKL)
    ZFIWP(JL)= ZIWGKG*PDZ(JL,IKL)
    ZFRWP(JL)= ZRWGKG*PDZ(JL,IKL)
    ZLWC(JL)=ZLWGKG+ZRWGKG
    ZIWC(JL)=ZIWGKG
!
!   Diagnostic for water path     
!
    IF (KRAD_DIAG >= 4 ) THEN
      PFLWP(JL,JK) = ZFLWP(JL)
      PFIWP(JL,JK) = ZFIWP(JL)
    END IF
!
! --- effective radius for water, rain and small ice particles
!
! WATER
!
    SELECT CASE (HEFRADL)
      CASE ('PRES')
!
! very old parametrization as f(pressure)
!
        ZRADLP(JL)=MAX(4., 10. + (100000.-PAP(JL,IKL))*3.5E-04)
      CASE ('OCLN')
!
! simple distinction between land (10) and ocean (13)
        IF (PLSM(JL) < 0.5) THEN
          ZRADLP(JL)=13.
        ELSE
          ZRADLP(JL)=10.
        END IF
      CASE ('MART')
!
! could be based on Martin et al., 1994, JAS
!
        IF (PLSM(JL) < 0.5) THEN
          ZASEA=150.
          ZD=0.33
          ZNTOT=-1.15E-03*ZASEA*ZASEA+0.963*ZASEA+5.30
        ELSE
          ZALND=900.
          ZD=0.43
          ZNTOT=-2.10E-04*ZALND*ZALND+0.568*ZALND-27.9
        ENDIF
        ZNUM=3.*ZLWC(JL)*(1.+3.*ZD*ZD)**2
        ZDEN=4.*RPI*ZNTOT*(1.+ZD*ZD)**3
        ZRADLP(JL)=100.*(ZNUM/ZDEN)**0.333
        ZRADLP(JL)=MAX(ZRADLP(JL), 4.)
        ZRADLP(JL)=MIN(ZRADLP(JL),16.)
        
!
      CASE ('C2R2')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CALCUL DU RAYON EFFECTIF RESULTANT COMME:
!calcul ds rayons surfacique et effectif pour chaque classe
!ZRADLS,ZRADRS,ZRADLP,ZRADRP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ZRHO(JL) =  (1./XRD)*(PAP(JL,IKL)/ PT(JL,IKL))
          ZRADLP(JL) = 0.
          ZRADLS(JL) = 0.
          ZRADRP(JL) = 0.
          ZRADRS(JL) = 0.
          IF (PCCT_C2R2(JL, IKL)>XCTMIN(2) .AND. PLWC(JL,IKL)>XRTMIN(2).AND.&
                ZFLWP(JL)/=0.0) THEN
            ZNUM =  PLWC(JL,IKL)/ XAC
            ZDEN = XCREC * (PCCT_C2R2(JL,IKL)*(PLWC(JL,IKL))**2)**(XLBEXC)
            ZRADLS(JL) = MIN( 0.5E6 * sqrt(ZDEN/PCCT_C2R2(JL,IKL)) , 100.)
            ZRADLP(JL) = MIN( 0.5E6 * (ZNUM/ZDEN)                  , 100.)
          ENDIF
           IF (PCRT_C2R2(JL, IKL)>XCTMIN(3) .AND. PRWC(JL,IKL)>XRTMIN(3).AND.&
                ZFRWP(JL)/=0.0 ) THEN
            ZNUM =  PRWC(JL,IKL)/ XAR
            ZDEN = XCRER * (PCRT_C2R2(JL,IKL)*(PRWC(JL,IKL))**2)**(XLBEXR)
            ZRADRS(JL) = 0.5E6 * sqrt(ZDEN/PCRT_C2R2(JL,IKL))
            ZRADRP(JL) = 0.5E6 * (ZNUM/ZDEN)
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   END SELECT
   IF (ZFLWP(JL).eq.0.0) ZRADLP(JL)=0.0
!
! --------------
!
! ICE
!
    IF (PT(JL,IKL) < RTICE) THEN
      ZTEMPC=PT(JL,IKL)-RTT
    ELSE
      ZTEMPC=RTICE-RTT
    ENDIF
!
    if (ZFIWP(JL).gt.0.0) then
    SELECT CASE (HEFRADI)
      CASE ('FX40')
!
!  fixed 40 micron effective radius
!
        ZRADIP(JL)= 40. 
        ZDESR(JL)=2.* ZRADIP(JL)            
      CASE ('LIOU')
!
!  ice particle effective radius =f(T) from Liou and Ou (1994)
!
        ZRADIP(JL)=326.3+ZTEMPC*(12.42+ZTEMPC*(0.197+ZTEMPC*0.0012))
        ZRADIP(JL)=MIN(ZRADIP(JL),60.)
        ZRADIP(JL)=MAX(ZRADIP(JL),30.)
        ZDESR(JL)=2.* ZRADIP(JL)   
      CASE ('SURI')
!
!  ice particle effective radius =f(T,IWC) from Sun and Rikus (1999)
!
! Sun, 2001
!
        IF (ZIWC(JL) > 0. ) THEN
          ZTEMPC = PT(JL,IKL)-83.15
          ZTCELS = PT(JL,IKL)-RTT
          ZFSR = 1.2351 +0.0105 * ZTCELS
!
! Sun, 2001 (corrected from Sun & Rikus, 1999)
!
          ZAIWC = 45.8966 * ZIWC(JL)**0.2214
          ZBIWC = 0.7957 * ZIWC(JL)**0.2535
          ZDESR(JL) = ZFSR * (ZAIWC + ZBIWC*ZTEMPC)
          ZDESR(JL) = MIN ( MAX( ZDESR(JL), 45.), 350.)
          ZRADIP(JL)= 0.5 * ZDESR(JL)
        ELSE
          ZDESR(JL) = 80.
          ZRADIP(JL)= 0.5 * ZDESR(JL)
        END IF  
!
      CASE ('C3R5')
!
! based on the prediction of the number concentrations
!
          ZRHO(JL) =  (1./XRD)*(PAP(JL,IKL)/ PT(JL,IKL)) 
!BVIE
!          ZRADIP(JL) = 0.0
          ZRADIP(JL) = 40.0
          IF (PCIT_C1R3(JL, IKL)>XCTMIN(4) .AND. PIWC(JL,IKL)>XRTMIN(4)) THEN
            ZMR(JL) = PIWC(JL,IKL)
            ZRADIP(JL) = 1.0E6 * XFREFFI*(ZMR(JL)/PCIT_C1R3(JL,IKL))**(XLBEXI)
            ZDESR(JL)  = 2.0 * ZRADIP(JL)
          END IF
    END SELECT
    ELSE
     ZRADIP(JL) = 0.0
     ZDESR(JL)  =0.0
    ENDIF  
!
!   Diagnostic for effective radius      
!
    PRADLP(JL,JK) = ZRADLP(JL)
    IF (KRAD_DIAG >= 4 ) THEN
      PRADIP(JL,JK) = ZRADIP(JL)
    END IF
  ENDDO
!
!          2.3    CLOUD SHORTWAVE OPTICAL PROPERTIES
!                 ----------------------------------
!   -------------------------
! --+ SW OPTICAL PARAMETERS +  Water clouds after Fouquart (1987)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)
!
! initialisation of absorbtion coefficient
  ZMSAID = 0.
  ZMSAIU = 0.
  DO JSW=1,NSW
    DO JL = IKIDIA,IKFDIA
      ZTOL = 0.
      ZGL  = 0.
      ZOL  = 0.
      ZTOR = 0.
      ZGR  = 0.
      ZOR  = 0.
      ZTOI = 0.
      ZGI  = 0.
      ZOI  = 0.
      IF ((HEFRADL=='C2R2').and.(HOPWSW=='FOUQ')) THEN
        GPROP_OP=((ZFIWP(JL)/= 0.).OR.(ZRADLS(JL)>1))
      ELSE
        GPROP_OP=((ZFIWP(JL)/= 0.).OR.(ZRADLP(JL)>1))
      ENDIF
      IF (GPROP_OP) THEN
          SELECT CASE ( HOPWSW)
            CASE('SLIN')
!
!-- SW: Slingo, 1989
!
               IF (ZRADLP(JL)>1.0) then
                ZTOL = ZSWFUDG * (ZFLWP(JL)+ZFRWP(JL))*(RASWCA(JSW)+RASWCB(JSW)/ZRADLP(JL))
                ZGL  = RASWCE(JSW)+RASWCF(JSW)*ZRADLP(JL)
                ZOL  = 1. - RASWCC(JSW)-RASWCD(JSW)*ZRADLP(JL)
               ENDIF
            CASE('FOUQ')         
!
!-- SW: Fouquart, 1991
!
              IF (HEFRADL == "C2R2") THEN
                IF (ZRADLS(JL)>1) then
                 ZTOL = ZSWFUDG * (ZFLWP(JL)+ZFRWP(JL))*(RYFWCA(JSW)+RYFWCB(JSW)/ZRADLS(JL))
                 ZGL  = RYFWCF(JSW)
                 ZOL  = RYFWCC(JSW)-RYFWCD(JSW)*EXP(-RYFWCE(JSW)*ZTOL)
                ENDIF
              ELSE IF (ZRADLP(JL)>1.) THEN
                 ZTOL = ZSWFUDG * (ZFLWP(JL)+ZFRWP(JL))*(RYFWCA(JSW)+RYFWCB(JSW)/ZRADLP(JL))
                 ZGL  = RYFWCF(JSW)
                 ZOL  = RYFWCC(JSW)-RYFWCD(JSW)*EXP(-RYFWCE(JSW)*ZTOL)
              END IF
           CASE('MALA')         
!
!-- SW:Optical thickness of Savijarvi (1997), asymetry factor of Fouquart
!(1991) and single scaterring albedo of Slingo 1989. Good compromise for the
!C2R2 and KHKO scheme in regard on the size distribution hypothesis (F. Malavelle
!intership M1, 2007) 
!
              IF (HEFRADL == "C2R2") THEN
                IF (ZRADLP(JL)>1) then
                 ZTOL =ZFLWP(JL)*(XSWSAVIA(JSW)+(XSWSAVIB(JSW)/ZRADLP(JL)))/ZRADLP(JL)
                 ZGL  = RYFWCF(JSW)
                 ZOL  = 1. - RASWCC(JSW)-RASWCD(JSW)*ZRADLP(JL)
! Test for Sc and fog but not to generalize :
! M.Mazoyer, O.Thouron effective radius does not exceed 100 microns
!                ZOL  = 1. - RASWCC(JSW)-RASWCD(JSW)*MIN(ZRADLP(JL),100.0)                 
                ENDIF
              ELSE IF (ZRADLP(JL)>1.) THEN
                 write(*,*)'PROGRAM ERROR: STOP'
                 write(*,*)'YOU USE A PARAMATERESISATION OF THE SW OPTICAL PROPERTIES'
                 write(*,*)'INADAPTED FOR THE 1 MOMENT SCHEME: SEE THE CEFRADL VARIABLE'
                 write(*,*)'IN YOUR NAMELIST'
!callabortstop
CALL ABORT
                 STOP
              END IF
          END SELECT
!
        IF (ZFIWP(JL) /= 0.) THEN
          SELECT CASE (HOPISW)
            CASE ('FULI')

!-- SW: Fu-Liou
!
              Z1RADI = 1./ ZRADIP(JL)
              ZBETAI = RFLAA0(JSW)+Z1RADI* RFLAA1(JSW)
              ZTOI = ZSWFUDG * ZFIWP(JL) * ZBETAI
              ZOMGI= RFLBB0(JSW)+ZRADIP(JL)*(RFLBB1(JSW) + ZRADIP(JL) &
                   &   *(RFLBB2(JSW)+ZRADIP(JL)* RFLBB3(JSW) ))            
              ZOI  = 1 - ZOMGI
              ZOMGP= RFLCC0(JSW)+ZRADIP(JL)*(RFLCC1(JSW) + ZRADIP(JL) &
                   &   *(RFLCC2(JSW)+ZRADIP(JL)* RFLCC3(JSW) )) 
              ZFDEL= RFLDD0(JSW)+ZRADIP(JL)*(RFLDD1(JSW) + ZRADIP(JL) &
                   &   *(RFLDD2(JSW)+ZRADIP(JL)* RFLDD3(JSW) )) 
              ZGI  = ((1.-ZFDEL)*ZOMGP + ZFDEL*3.) / 3.            
            CASE ('EBCU') 
!
!-- SW: Ebert-Curry          
!
              ZTOI = ZSWFUDG * ZFIWP(JL)*(REBCUA(JSW)+REBCUB(JSW)/ZRADIP(JL))
              ZGI  = REBCUE(JSW)+REBCUF(JSW)*ZRADIP(JL)
              ZOI  = 1 - REBCUC(JSW)-REBCUD(JSW)*ZRADIP(JL)

            CASE ('FU96')   
!
!-- SW: Fu 1996 
!
              Z1RADI = 1. / ZDESR(JL)
              ZBETAI = RFUAA0(JSW)+Z1RADI* RFUAA1(JSW)
              ZTOI =  ZSWFUDG * ZFIWP(JL) * ZBETAI
              ZOMGI= RFUBB0(JSW)+ZDESR(JL)*(RFUBB1(JSW) + ZDESR(JL) &
               &   *(RFUBB2(JSW)+ZDESR(JL)* RFUBB3(JSW) ))            
              ZOI  = 1.- ZOMGI
              ZGI  = RFUCC0(JSW)+ZDESR(JL)*(RFUCC1(JSW) + ZDESR(JL) &
               &   *(RFUCC2(JSW)+ZDESR(JL)* RFUCC3(JSW) )) 
          END SELECT
        ENDIF
!
!
!  - MIX of WATER and ICE CLOUDS (notice that R properties are either included
!                                          in L properties     or     disabled)
!

        ZTAUMX= ZTOL         + ZTOI         + ZTOR
        ZOMGMX= ZTOL*ZOL     + ZTOI*ZOI     + ZTOR*ZOR
        ZASYMX= ZTOL*ZOL*ZGL + ZTOI*ZOI*ZGI + ZTOR*ZOR*ZGR
!
        ZASYMX= ZASYMX/ZOMGMX
        ZOMGMX= ZOMGMX/ZTAUMX
!
! --- SW FINAL CLOUD OPTICAL PARAMETERS
!
        ZCLDSW(JL,JK)     = ZCLFR(JL,IKL)
        ZTAU(JL,JSW,JK)   = ZTAUMX
        ZOMEGA(JL,JSW,JK )= ZOMGMX
        ZCG(JL,JSW,JK)    = ZASYMX
      ENDIF
    ENDDO
  ENDDO
!
!          2.4    CLOUD LONGWAVE OPTICAL PROPERTIES FOR EC-OPE
!                 --------------------------------------------
!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Smith and Shi (1992)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)
!
  IF (.NOT.LRRTM) THEN
!
    DO JL = IKIDIA,IKFDIA
      ZALFICE(JL)=0.
      ZGAMICE(JL)=0.
      ZBICE(JL)=0.
      ZTICE(JL)=(PT(JL,IKL)-TSTAND)/TSTAND
      ZBICFU(JL)=0.
      ZKICFU1(JL)=0.
      ZKICFU2(JL)=0.
    ENDDO
    DO JNU= 1,NSIL
      DO JL = IKIDIA,IKFDIA
        ZRES(JL)  = XP(1,JNU)+ZTICE(JL)*(XP(2,JNU)+ZTICE(JL)*(XP(3,&
             &JNU)&
             &+ZTICE(JL)*(XP(4,JNU)+ZTICE(JL)*(XP(5,JNU)+ZTICE(JL)*(XP(6,&
             &JNU)&
             &)))))
        ZBICE(JL) = ZBICE(JL) + ZRES(JL)
        ZGAMICE(JL) = ZGAMICE(JL) + REBCUI(JNU)*ZRES(JL)
        ZALFICE(JL) = ZALFICE(JL) + REBCUJ(JNU)*ZRES(JL)
      ENDDO
    ENDDO
!
!-- Fu et al. (1998) with M'91 LW scheme    
!
    DO JRTM=1,16
      DO JL=IKIDIA,IKFDIA
        IF (PT(JL,IKL) < 339. .AND. PT(JL,IKL) >= 160.) THEN
          INDLAY=PT(JL,IKL)-159.
          ZTBLAY =PT(JL,IKL)-INT(PT(JL,IKL))
        ELSE IF (PT(JL,IKL) >= 339. ) THEN
          INDLAY=180
          ZTBLAY =PT(JL,IKL)-339.
        ELSE IF (PT(JL,IKL) < 160.) THEN
          INDLAY=1
          ZTBLAY =PT(JL,IKL)-160.
        END IF      
        ZADDPLK = TOTPLNK(INDLAY+1,JRTM)-TOTPLNK(INDLAY,JRTM)
        ZPLANCK = DELWAVE(JRTM) * (TOTPLNK(INDLAY,JRTM) + ZTBLAY*ZADDPLK)
        ZBICFU(JL) = ZBICFU(JL) + ZPLANCK
        
        ZMSAID=0.
        ZMSAIU=0.
        IF (ZIWC(JL) > 0.  ) THEN
!
! ice cloud spectral emissivity from Fu & Liou (1993)
!
          ZRATIO= 0.5 / ZRADIP(JL)
          ZMSAID = RFULIO(JRTM,1) + ZRATIO&
             &*(RFULIO(JRTM,2) + ZRATIO*RFULIO(JRTM,3))
          ZKICFU1(JL) = ZKICFU1(JL)+ ZMSAID*ZPLANCK
!
! ice cloud spectral emissivity from Fu et al (1998)
!
          Z1RADI = 1. / ZDESR(JL)
          ZMSAID = RFUETA(JRTM,1) + Z1RADI&
             &*(RFUETA(JRTM,2) + Z1RADI*RFUETA(JRTM,3))
          ZKICFU2(JL) = ZKICFU2(JL)+ ZMSAID*ZPLANCK
        END IF  
      END DO
    END DO
!
    DO JL = IKIDIA,IKFDIA
      ZGAMICE(JL) = ZGAMICE(JL) / ZBICE(JL)
      ZALFICE(JL) = ZALFICE(JL) / ZBICE(JL)
      ZKICFU1(JL) = ZKICFU1(JL) / ZBICFU(JL)
      ZKICFU2(JL) = ZKICFU2(JL) / ZBICFU(JL)
      ZMSALD = 0.
      ZMSARD = 0.
      ZMSALU = 0.
      ZMSARU = 0.

!
      GPROP_OP=((ZFIWP(JL)/= 0.).OR.(ZRADLP(JL)>1))
      IF (GPROP_OP) THEN
        IF(ZRADLP(JL)>1) THEN 
          SELECT CASE (HOPWLW)
          CASE('SAVI')  
!
! water cloud emissivity from Savijarvi, 1997
!
            IF (ZRADLP(JL)>1) THEN
             ZMSALU= 0.2441-0.0105*ZRADLP(JL)
             ZMSALD= 1.2154*ZMSALU
            ENDIF
            
          CASE('SMSH') 
!
! water cloud emissivity from Smith-Shi, 1992
!
            IF (ZRADLP(JL)>1) THEN
             ZMULTL=1.2-0.006*ZRADLP(JL)
             ZMSALD= 0.158*ZMULTL
             ZMSALU= 0.130*ZMULTL
            ENDIF
                   
          END SELECT
        END IF
!
        IF(ZFIWP(JL) /= 0.) THEN
          SELECT CASE (HOPILW) 
          CASE ('SMSH')
!
! ice cloud emissivity from Smith-Shi (1992)
!
            ZMULTI=1.2-0.006*ZRADIP(JL)
            ZMSAID= 0.113*ZMULTI
            ZMSAIU= 0.093*ZMULTI
          CASE ('EBCU')
!
! ice cloud emissivity from Ebert-Curry (1992)
!
            ZMSAID= 1.66*(ZALFICE(JL)+ZGAMICE(JL)/ZRADIP(JL))
            ZMSAIU= ZMSAID
          CASE ('FULI')
!
! ice cloud emissivity from Fu & Liou (1993)
!
            ZMSAID= 1.66 * ZKICFU1(JL)
            ZMSAIU= ZMSAID
          CASE ('FU98')
!
! ice cloud emissivity from Fu et al. (1998)
!
            ZMSAID= 1.66 * ZKICFU2(JL)
            ZMSAIU= ZMSAID
          END SELECT
        END IF
!
! introduce inhomogeneity factor also in LW         
!
        IF (HEFRADL == "C2R2") THEN
         ZZFLWP= ZFLWP(JL) * ZLWFUDG
         ZZFRWP= ZFRWP(JL) * ZLWFUDG
        ELSE
         ZZFLWP= (ZFLWP(JL)+ZFRWP(JL)) * ZLWFUDG
        ENDIF
        ZZFIWP= ZFIWP(JL) * ZLWFUDG
!
! effective cloudiness accounting for condensed water
!
        ZCLDLD(JL,JK) = ZCLFR(JL,IKL)*(1-EXP(-ZMSALD*ZZFLWP-ZMSARD*ZZFRWP-ZMSAID* &
             &ZZFIWP))
        ZCLDLU(JL,JK) = ZCLFR(JL,IKL)*(1-EXP(-ZMSALU*ZZFLWP-ZMSARU*ZZFRWP-ZMSAIU* &
             &ZZFIWP))
      ENDIF
    ENDDO
  ELSE
!
!          2.5    CLOUD LONGWAVE OPTICAL PROPERTIES FOR RRTM
!                 ------------------------------------------
!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Savijarvi (1998)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)
!
! No need for a fixed diffusivity factor, accounted for spectrally below
! The detailed spectral structure does not require defining upward and
! downward effective optical properties
!
    DO JRTM=1,16
      DO JL = IKIDIA,IKFDIA
        ZTAUCLD(JL,JK,JRTM) = 0.
        ZMSALD = 0.
        ZMSARD = 0.
        ZMSAID = 0.
!
      IF (HEFRADL=='C2R2') THEN
          GPROP_OP=((ZFIWP(JL)/= 0.).OR.(ZRADLP(JL)>1).OR.(ZRADRP(JL)>1))
      ELSE
          GPROP_OP=((ZFIWP(JL)/= 0.).OR.(ZRADLP(JL)>1))
      ENDIF    
      IF (GPROP_OP) THEN
          IF(ZRADLP(JL)>1) THEN 
            SELECT CASE (HOPWLW)
            CASE ('LILI') 
!
! water cloud spectral emissivity from Lindner and Li (2000)
!
               IF (ZRADLP(JL)>1) THEN
                Z1RADL = 1. / ZRADLP(JL)
                ZMSALD = XLWLILI(JRTM,1)+ZRADLP(JL)*XLWLILI(JRTM,2)+ Z1RADL*&
                &       (XLWLILI(JRTM,3) + Z1RADL*(XLWLILI(JRTM,4) + Z1RADL*&
                &        XLWLILI(JRTM,5) ))
                ZMSALD=ZMSALD/1.66
               ENDIF
                      
            CASE ('SAVI')  
!
! water cloud spectral emissivity from Savijarvi (1997)
!
            IF (ZRADLP(JL)>1) THEN
              if (JRTM.LE.10) THEN
                  ZMSALD= XLWSAVI(JRTM,1)-XLWSAVI(JRTM,2)*ZRADLP(JL)
                  ZMSALD=ZMSALD  / 1.66
              elseif (JRTM.GE.11) THEN
                  ZMSALD=XLWSAVI(JRTM,1) *exp(-XLWSAVI(JRTM,2)*ZRADLP(JL))
                  ZMSALD=ZMSALD  / 1.66               
              endif
            ENDIF
            CASE ('SMSH') 
!
! water cloud total emissivity from Smith and Shi (1992)
!
              IF (ZRADLP(JL)>1) THEN
               ZMULTL=1.2-0.006*ZRADLP(JL)
               ZMSALD= 0.158*ZMULTL / 1.66
              ENDIF
!Parameterisation based on the size distribution hypothesis 
!use in the KHKO scheme. Also adapted to the C2R2 scheme
            CASE ('MALA') 
              IF (ZRADLP(JL)>1) THEN
               ZLOGREFF=LOG10(ZRADLP(JL))
               ZMSALD=XLWC2R2(JRTM,1)+ZLOGREFF*(XLWC2R2(JRTM,2)+&
                   ZLOGREFF*(XLWC2R2(JRTM,3)+ZLOGREFF*XLWC2R2(JRTM,4)))
               ZMSALD=10**ZMSALD
              ENDIF
            END SELECT
          END IF
!
!Parameterisation of F. Malavelle, intership of M2, 2008
!Adapted to the size distribution of boundary layer cloud
!Adapted to the C2R2 and KHKO scheme
         IF(ZRADRP(JL)>1) THEN 
!INTRODUIRE BOULOT DE FLORENT
         ENDIF
          IF(ZFIWP(JL) /= 0.) THEN
            SELECT CASE (HOPILW)  
            CASE ('FULI')     
!
! ice cloud spectral emissivity from Fu & Liou (1993)
!
              ZRATIO= 0.5 / ZRADIP(JL)
              ZMSAID = RFULIO(JRTM,1) + ZRATIO&
                   &*(RFULIO(JRTM,2) + ZRATIO*RFULIO(JRTM,3))
            CASE('EBCU')     
!
! ice cloud spectral emissivity from Ebert-Curry (1992)
!
              ZMSAID= REBCUH(JRTM)+REBCUG(JRTM)/ZRADIP(JL)
            CASE('FU98')
!
! ice cloud spectral emissivity from Fu et al (1998)
!
              Z1RADI = 1. / ZDESR(JL)
              ZMSAID = RFUETA(JRTM,1) + Z1RADI&
              &*(RFUETA(JRTM,2) + Z1RADI*RFUETA(JRTM,3))
            END SELECT
          END IF
!
          IF (HEFRADL == "C2R2") THEN
            ZTAUD = ZLWFUDG * ZMSALD * ZFLWP(JL)+&
                    ZLWFUDG * ZMSARD * ZFRWP(JL)+&
                    ZLWFUDG * ZMSAID * ZFIWP(JL)
          ELSE
                    ZTAUD = ZLWFUDG * ZMSALD * (ZFLWP(JL)+ZFRWP(JL))+&
                            ZLWFUDG * ZMSAID * ZFIWP(JL)
         ENDIF
!
! Diffusivity correction within clouds from Savijarvi
!          ZDIFFD=MIN(MAX(1.517-0.156*LOG(ZTAUD) , 1) , 2)
!
         ZDIFFD=1.66
         ZTAUCLD(JL,JK,JRTM) = ZTAUD*ZDIFFD
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!
ENDDO
!
NUAER = NUA
NTRAER = NTRA

!
!     ------------------------------------------------------------------
!*         2.6    DIFFUSIVITY FACTOR OR SATELLITE VIEWING ANGLE
!                 ---------------------------------------------
DO JL = IKIDIA,IKFDIA
  ZVIEW(JL) = DIFF
ENDDO
!     ------------------------------------------------------------------
!*         3.     CALL LONGWAVE RADIATION CODE
!                 ----------------------------
!
!*         3.1    FULL LONGWAVE RADIATION COMPUTATIONS
!                 ------------------------------------
IF ( .NOT. LRRTM) THEN
  CALL LW ( IKIDIA , IKFDIA , KLON  , KLEV , NMODE, &
       ZCCO2_RAD , ZCLDLD, ZCLDLU,                     &
       PDP   , ZDT0  , ZEMIS , ZEMIW,              &
       ZPMB  , ZOZON , ZTL,                        &
       ZAER_LW  , ZTAVE , ZVIEW , PQ,                 &
       ZCOOLR, ZCOOLC, ZEMIT , ZFLUX_LW, ZFLUX_CLW )
ELSE
!*         3.2    FULL LONGWAVE RADIATION COMPUTATIONS - RRTM
!                 ------------------------------------   ----
!  i)  pass POZN (ozone mmr concentration) to RRTM; remove pressure
!      weighting applied to ZOZON in driverMC (below)
!  ii) pass ZEMIS and ZEMIW to RRTM; return ZEMIT from RRTM
!  iii)pass ZTAUCLD, cloud optical depths (water+ice) to RRTM, 
!      computed from equations above
!  iv) pass ECRT arrays to RRTM arrays in interface routine ECRTATM
!      in module rrtm_ecrt.f
!
  DO JK = 1, KLEV
    DO JL = IKIDIA,IKFDIA
!       ZOZN(JL,JK) = ZOZON(JL,JK)/PDP(JL,JK) 
      ZOZN(JL,JK) = POZON(JL,JK)/PDP(JL,JK) ! Quentin because climatology in volume mixing ratio, not mmr
    ENDDO
  ENDDO
!
  CALL RRTM_RRTM_140GP(IKIDIA,IKFDIA,KLON,KLEV,   &
       ZAER_LW,PAPH,PAP,PTS,PTH,PT,ZEMIS,ZEMIW,      &
       PQ    , ZCCO2_RAD , ZOZN  , ZCLDSW  , ZTAUCLD, &
       ZEMIT , ZFLUX_LW , ZFLUX_CLW , ZTCLEAR )
ENDIF

!     ------------------------------------------------------------------
!*         4.     CALL SHORTWAVE RADIATION CODE
!                 -----------------------------
ZRMUZ=0.
DO JL = IKIDIA,IKFDIA
  ZRMUZ = MAX (ZRMUZ, ZMU0(JL))
ENDDO
!
IF (ZRMUZ > 0.) THEN
      CALL SW ( IKIDIA , IKFDIA , KLON  , KLEV  , KAER,     &
           PRII0 , ZCCO2_RAD , ZPSOL , ZALBD , ZALBP , PQ   , PQS,  &
           ZMU0  , ZCG   , ZCLDSW, PDP   , ZOMEGA, ZOZ  , ZPMB, &
           ZTAU  , ZTAVE , ZAER_SW,                                &
           ZHEATR, ZFSDWN, ZFSUP , ZHEATC, ZFCDWN, ZFCUP,       &
           ZFSDNN, ZFSDNV, ZFSUPN, ZFSUPV,                      &
           ZFCDNN, ZFCDNV, ZFCUPN, ZFCUPV,                      &
           ZSUDU, ZUVDF, ZPARF,ZDIFFS,ZDIRFS,                   &
           ODUST,PPIZA_DST,PCGA_DST,PTAUREL_DST)
ENDIF

!     ------------------------------------------------------------------
!
!*         5.     FLUXES, RADIATIVE TENDENCIES, DIAGNOSTICS
!               ------------------------------------------------
! note : important
! After the SW and LW routines, radiative fluxes array are consistent with
! the orientation of MNH vertical grid : i.e FLUX (:,1) is the flux at the surface
! and FLUX (:,KLEV+1) is the TOA flux. Remember that it is not the case for input  
! arrays coming from radiation.f90 (e.g PDP(:,1) correponds to TOA) 
!
!          5.0 net fluxes
!
DO JK = 1 , KLEV+1
  DO JL = IKIDIA,IKFDIA
    ZFLS(JL,JK) = ZFSDWN(JL,JK) - ZFSUP(JL,JK)
    ZFLT(JL,JK) = - ZFLUX_LW(JL,1,JK) - ZFLUX_LW(JL,2,JK)
    ZFCS(JL,JK) = ZFCDWN(JL,JK) - ZFCUP(JL,JK)
    ZFCT(JL,JK) = - ZFLUX_CLW(JL,1,JK) - ZFLUX_CLW(JL,2,JK)
  ENDDO
ENDDO
!
!          5.1 Radiative tendencies in (K/DAY)
!
! calculation of DP term (grid orientation consistent with net fluxes,
! ( see note above) 
!computation of CPH (RCDAY is computed with CPD)
DO JKL = 1 , KLEV
  DO JL = IKIDIA,IKFDIA
    JK = KLEV+1  - JKL
    ZDPGCP(JL,JKL) = RCDAY/ PDP(JL, JK)*XCPD/(XCPD+XCPV*PQ(JL, JK)/(1.-PQ(JL, JK)) &
                          +XCL*(PQLWC(JL, JK)+PQRWC(JL, JK))+XCI*PQIWC(JL, JK))
  END DO
END DO
!
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
!         5.2  Top and Ground fluxes (IR-VIS-NIR)
!
PFLUX_TOP_GND_IRVISNIR(:, 1) = - ZFLUX_LW (:,1,KLEV+1)
PFLUX_TOP_GND_IRVISNIR(:, 2) = - ZFSUPV(:)
PFLUX_TOP_GND_IRVISNIR(:, 3) = - ZFSUPN(:)
PFLUX_TOP_GND_IRVISNIR(:, 4) = - ZFLUX_LW (:,2,1)
PFLUX_TOP_GND_IRVISNIR(:, 5) =   ZFSDNV(:)
PFLUX_TOP_GND_IRVISNIR(:, 6) =   ZFSDNN(:)
!
!Clear-sky
!
PFLUX_TOP_GND_IRVISNIR_CS(:, 1) = - ZFLUX_CLW (:,1,KLEV+1)
PFLUX_TOP_GND_IRVISNIR_CS(:, 2) = - ZFCUPV(:)
PFLUX_TOP_GND_IRVISNIR_CS(:, 3) = - ZFCUPN(:)
PFLUX_TOP_GND_IRVISNIR_CS(:, 4) = - ZFLUX_CLW (:,2,1)
PFLUX_TOP_GND_IRVISNIR_CS(:, 5) =   ZFCDNV(:)
PFLUX_TOP_GND_IRVISNIR_CS(:, 6) =   ZFCDNN(:)
!
IF (SIZE(PSFSWDIR,2)==2) THEN
  PSFSWDIR (:,1) = ZFSDNV(:)
  PSFSWDIR (:,2) = ZFSDNN(:)
ELSEIF(SIZE(PSFSWDIR,2)==1) THEN
  PSFSWDIR (:,1) = ZFSDNV(:) + ZFSDNN(:)
  PSFSWDIF (:,:) = 0.
ELSEIF(SIZE(PSFSWDIR,2)==6) THEN
  PSFSWDIR=ZDIRFS
  PSFSWDIF=ZDIFFS
END IF
!
!             5.2  Radiation Diagnostics 
!
PFSDWN(:,:) = ZFSDWN(:,:)
PFSUP (:,:) = ZFSUP (:,:)
PFLUX_LW(:,:,:) = ZFLUX_LW(:,:,:)
IF (KRAD_DIAG >= 1) THEN
  PFLT(:,:) = ZFLT(:,:)
  PFLS(:,:) = ZFLS(:,:)
END IF
!
IF (KRAD_DIAG >=2) THEN
  PFCDWN(:,:) = ZFCDWN(:,:) 
  PFCUP (:,:) = ZFCUP (:,:)
  PFLUX_CLW(:,:,:) = ZFLUX_CLW(:,:,:)
  PFCS(:,:) = ZFCS(:,:)
  PFCT(:,:) = ZFCT(:,:)
END IF
!
IF (KRAD_DIAG >= 3) THEN 
  ! provisoire 
  PPLAN_ALB_VIS(:) = XUNDEF
  PPLAN_ALB_NIR(:) = XUNDEF
  PPLAN_TRA_VIS(:) = XUNDEF 
  PPLAN_TRA_NIR(:) = XUNDEF
  PPLAN_ABS_VIS(:) = XUNDEF 
  PPLAN_ABS_NIR(:) = XUNDEF

END IF
!
IF (KRAD_DIAG >=4) THEN 
  IF (.NOT.LRRTM) THEN
    PEFCL_LWD  (:,:) = ZCLDLD (:,:)
    PEFCL_LWU    (:,:) = ZCLDLU (:,:)
    PEFCL_RRTM (:,:) = XUNDEF
  ELSE
    PEFCL_LWD  (:,:) = XUNDEF
    PEFCL_LWU    (:,:) = XUNDEF
    PEFCL_RRTM (:,:)  = SUM ( ZTAUCLD(:,:,:), DIM=3 ) / 16 ! average on 16 RRTM bands
  END IF
   
  PCLSW_TOTAL (:,:)   = ZCLDSW(:,:) 
  PTAU_TOTAL  (:,:,:) = ZTAU(:,:,:)
  POMEGA_TOTAL(:,:,:) = ZOMEGA(:,:,:)
  PCG_TOTAL   (:,:,:) = ZCG(:,:,:)
END IF
!

END SUBROUTINE ECMWF_RADIATION_VERS2
