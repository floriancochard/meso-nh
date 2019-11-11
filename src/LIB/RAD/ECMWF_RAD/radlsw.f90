!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
SUBROUTINE RADLSW &
 &( KIDIA, KFDIA , KLON , KTDIA, KLEV  , KMODE, KAER, KBOX, NBOX &
 &, NDUMP, KLWRAD &
 &, PRII0 &
 &, PAER , PALBD , PALBP, PAPH , PAP &
 &, PCCO2, PFRCL , PDP  , PEMIS, PEMIW , PLSM , PMU0, POZON &
 &, PQ   , PQIWP , PQLWP, PSQIW, PSQLW , PQS  , PQRAIN, PRAINT &
 &, PRLVRI,PRLVRL, PTH  , PT   , PTS   , PNBAS, PNTOP &
 &, PEMIT, PFCT  , PFLT , PFCS , PFLS  , PFRSOD, PSUDU, PUVDF, PPARF &
 &, PFDCT, PFUCT , PFDLT, PFULT, PFDCS , PFUCS , PFDLS, PFULS &
 &, ASWBOX, OLRBOX, SLWBOX, SSWBOX, TAUBOX, PCLBX &
 &)

!**** *RADLSW* - INTERFACE TO ECMWF LW AND SW RADIATION SCHEMES

!     PURPOSE.
!     --------
!           CONTROLS RADIATION COMPUTATIONS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
! PAER   : (KLON,6,KLEV)     ; OPTICAL THICKNESS OF THE AEROSOLS
! PALBD  : (KLON,NSW)        ; SURF. SW ALBEDO FOR DIFFUSE RADIATION
! PALBP  : (KLON,NSW)        ; SURF. SW ALBEDO FOR PARALLEL RADIATION
! PAPH   : (KLON,KLEV+1)     ; HALF LEVEL PRESSURE
! PAP    : (KLON,KLEV)       ; FULL LEVEL PRESSURE
! PCCO2  :                   ; CONCENTRATION IN CO2 (PA/PA)
! PFRCL  : (KLON,KLEV)       ; CLOUD FRACTIONAL COVER
! PDP    : (KLON,KLEV)       ; LAYER PRESSURE THICKNESS
! PEMIS  : (KLON)            ; SURFACE LW EMISSIVITY
! PEMIW  : (KLON)            ; SURFACE LW WINDOW EMISSIVITY
! PLSM   : (KLON)            ; LAND-SEA MASK
! PMU0   : (KLON)            ; SOLAR ANGLE
! PNBAS  : (KLON)            ; INDEX OF BASE OF CONVECTIVE LAYER
! PNTOP  : (KLON)            ; INDEX OF TOP OF CONVECTIVE LAYER
! POZON  : (KLON,KLEV)       ; CONCENTRATION IN OZONE (PA/PA)
! PQ     : (KLON,KLEV)       ; SPECIFIC HUMIDITY PA/PA
! PQIWP  : (KLON,KLEV)       ; SOLID  WATER KG/KG
! PQLWP  : (KLON,KLEV)       ; LIQUID WATER KG/KG
! PQS    : (KLON,KLEV)       ; SATURATION WATER VAPOR  KG/KG
! PQRAIN : (KLON,KLEV)       ; RAIN WATER KG/KG
! PRAINT : (KLON,KLEV)       ; RAIN RATE (m/s)
! PRLVRI : (KLON,KLEV)       ; RELATIVE VARIANCE OF ICE WATER
! PRLVRL : (KLON,KLEV)       ; RELATIVE VARIANCE OF LIQUID WATER
! PTH    : (KLON,KLEV+1)     ; HALF LEVEL TEMPERATURE
! PT     : (KLON,KLEV)       ; FULL LEVEL TEMPERATURE
! PTS    : (KLON)            ; SURFACE TEMPERATURE
!     ==== OUTPUTS ===
! PFCT   : (KLON,KLEV+1)     ; CLEAR-SKY LW NET FLUXES
! PFLT   : (KLON,KLEV+1)     ; TOTAL LW NET FLUXES
! PFCS   : (KLON,KLEV+1)     ; CLEAR-SKY SW NET FLUXES
! PFLS   : (KLON,KLEV+1)     ; TOTAL SW NET FLUXES
! PFRSOD : (KLON)            ; TOTAL-SKY SURFACE SW DOWNWARD FLUX
! PEMIT  : (KLON)            ; SURFACE TOTAL LONGWAVE EMISSIVITY
! PSUDU  : (KLON)            ; SOLAR RADIANCE IN SUN'S DIRECTION
! PUVDF  : (KLON)            ; SURFACE DOWNWARD U.V. RADIATION
! PPARF  : (KLON)            ; PHOTOSYNTHETICALLY ACTIVE RADIATION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHORS.
!     --------
!        J.-J. MORCRETTE         *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 88-02-04
!        J.-J. MORCRETTE 94-11-15 DIRECT/DIFFUSE SURFACE ALBEDO
!        08/96: J.-J. Morcrette/Ph. Dandin: tests of eff. radius param.
!        9909 : JJMorcrette effect.radius + inhomogeneity factors
!        JJMorcrette 990128 : sunshine duration
!        JJMorcrette : 990831 RRTM-140gp
!-----------------------------------------------------------------------

#include "tsmbkind.h"

!USE YOMCT3   , ONLY : NSTEP
USE OYOMCST   , ONLY : RG       ,RD       ,RTT      ,RPI
USE OYOERAD   , ONLY : NSW      ,LRRTM    ,LINHOM, &
            &LOIFUEC, LTEMPDS, LOWASYF, LOWHSSS, NRADIP, NRADLP, &
            &NICEOPT, NLIQOPT, NOVLP  , NHOWINH, RMINICE
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
USE YOERDU   , ONLY : NUAER    ,NTRAER   ,REPLOG   ,REPSC    ,DIFF
USE OYOERDI   , ONLY : REPCLC
USE YOETHF   , ONLY : RTICE
USE YOEPHLI  , ONLY : LPHYLIN
USE OYOERRTWN , ONLY : NG        ,NSPA      ,NSPB      ,WAVENUM1  ,&
           &WAVENUM2  ,DELWAVE   ,TOTPLNK   ,TOTPLK16
USE YOEDBUG  , ONLY : LDEBUG


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KAER
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KMODE
INTEGER_M :: KTDIA
INTEGER_M :: KBOX
INTEGER_M :: NBOX
INTEGER_M :: NDUMP, KLWRAD

!     DUMMY REAL SCALARS
REAL_B :: PRII0



!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------
REAL_B :: PALBD(KLON,NSW) , PALBP(KLON,NSW)
REAL_B :: PEMIS(KLON)     , PEMIW(KLON)
REAL_B :: PLSM(KLON)      , PMU0(KLON)
REAL_B :: PCCO2           , POZON(KLON,KLEV)
REAL_B :: PTS(KLON)       , PNBAS(KLON)     , PNTOP(KLON)
REAL_B :: PT (KLON,KLEV)  , PAP (KLON,KLEV)
REAL_B :: PTH(KLON,KLEV+1), PAPH(KLON,KLEV+1)
REAL_B :: PDP(KLON,KLEV)
REAL_B :: PQ (KLON,KLEV)  , PQS(KLON,KLEV)
REAL_B :: PQIWP(KLON,KLEV), PQLWP(KLON,KLEV), PQRAIN(KLON,KLEV)
REAL_B :: PRAINT(KLON,KLEV)
REAL_B :: PRLVRI(KLON,KLEV),PRLVRL(KLON,KLEV)
REAL_B :: PSQIW(KLON,KLEV), PSQLW(KLON,KLEV)
REAL_B :: PFRCL(KLON,KLEV), PCLFR(KLON,KLEV), PCLBX(KLON,100,KLEV)
REAL_B :: PAER (KLON,6,KLEV)

!     ==== COMPUTED IN RADLSW ===
REAL_B :: PFCS(KLON,KLEV+1), PFCT(KLON,KLEV+1)
REAL_B :: PFLS(KLON,KLEV+1), PFLT(KLON,KLEV+1)
REAL_B :: PFRSOD(KLON)     , PEMIT(KLON)
REAL_B :: PSUDU(KLON)      , PUVDF(KLON)        , PPARF(KLON)
REAL_B :: PFDCT(KLON,KLEV+1), PFUCT(KLON,KLEV+1)
REAL_B :: PFDLT(KLON,KLEV+1), PFULT(KLON,KLEV+1)
REAL_B :: PFDCS(KLON,KLEV+1), PFUCS(KLON,KLEV+1)
REAL_B :: PFDLS(KLON,KLEV+1), PFULS(KLON,KLEV+1)

REAL_B :: ASWBOX(KLON, 100), OLRBOX(KLON, 100)
REAL_B :: SLWBOX(KLON, 100), SSWBOX(KLON, 100), TAUBOX(KLON, 100)

!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------
!     -----------------------------------------------------------------

!-- ARRAYS FOR LOCAL VARIABLES -----------------------------------------

INTEGER_M :: IBAS(KLON)     , ITOP(KLON)

REAL_B ::&
    &ZALBD(KLON,NSW)    , ZALBP(KLON,NSW)&
  &, ZCG(KLON,NSW,KLEV) , ZOMEGA(KLON,NSW,KLEV)&
  &, ZTAU (KLON,NSW,KLEV) &
  &, ZTAUCLD(KLON,KLEV,16), ZTCLEAR(KLON)
REAL_B ::&
    &ZCLDLD(KLON,KLEV)  , ZCLDLU(KLON,KLEV)&
  &, ZCLDSW(KLON,KLEV)  , ZCLD0(KLON,KLEV)&
  &, ZDT0(KLON)        &
  &, ZEMIS(KLON)        , ZEMIW(KLON)&
  &, ZFLUX (KLON,2,KLEV+1)                 , ZFLUC(KLON,2,KLEV+1)&
  &, ZFIWP(KLON)        , ZFLWP(KLON)      , ZFRWP(KLON)&
  &, ZIWC(KLON)         , ZLWC(KLON)&
  &, ZBICFU(KLON)       , ZKICFU1(KLON)    , ZKICFU2(KLON)&
!cc            , ZRWC(KLON)
  &, ZMU0(KLON)         , ZOZ(KLON,KLEV)   , ZOZN(KLON,KLEV)&
  &, ZOZON(KLON,KLEV)   , ZPMB(KLON,KLEV+1), ZPSOL(KLON)&
  &, ZTAVE (KLON,KLEV)  , ZTL(KLON,KLEV+1)&
  &, ZVIEW(KLON)
REAL_B ::&
    &ZFCDWN(KLON,KLEV+1), ZFCUP(KLON,KLEV+1)&
  &, ZFSDWN(KLON,KLEV+1), ZFSUP(KLON,KLEV+1)&
  &, ZFSUPN(KLON)       , ZFSUPV(KLON)&
  &, ZFCUPN(KLON)       , ZFCUPV(KLON)&
  &, ZFSDNN(KLON)       , ZFSDNV(KLON)&
  &, ZFCDNN(KLON)       , ZFCDNV(KLON)&
  &, ZCOOLR(KLON,KLEV)  , ZCOOLC(KLON,KLEV)&
  &, ZHEATR(KLON,KLEV)  , ZHEATC(KLON,KLEV)
REAL_B ::&
    &ZALFICE(KLON)      , ZGAMICE(KLON)     , ZBICE(KLON),  ZDESR(KLON) &
  &, ZRADIP(KLON)       , ZRADLP(KLON)      , ZCFUDG(KLON)&
!cc           , ZRADRD(KLON)
  &, ZRAINT(KLON)       , ZRES(KLON)&
  &, ZTICE(KLON)        , ZEMIT(KLON)       , ZTAUINT(KLON)
REAL_B :: ZSUDU(KLON)   , ZUVDF(KLON)       , ZPARF(KLON),  ZCOL(KLON) &
  &, ZTCC(KLON)         , ZTCA(KLON)

!-- box-type arrays

REAL_B :: CPFCS(KLON,KLEV+1) , CPFCT(KLON,KLEV+1)
REAL_B :: CPFLS(KLON,KLEV+1) , CPFLT(KLON,KLEV+1)
REAL_B :: CPFRSOD(KLON)      , CPEMIT(KLON)
REAL_B :: CPSUDU(KLON)       , CPUVDF(KLON)       , CPPARF(KLON)
REAL_B :: CPFDCT(KLON,KLEV+1), CPFUCT(KLON,KLEV+1)
REAL_B :: CPFDLT(KLON,KLEV+1), CPFULT(KLON,KLEV+1)
REAL_B :: CPFDCS(KLON,KLEV+1), CPFUCS(KLON,KLEV+1)
REAL_B :: CPFDLS(KLON,KLEV+1), CPFULS(KLON,KLEV+1)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IKL, JAE, JK, JKL, JKLP1, JKP1, JL, JNU, JRTM, JSW &
  &, NBOXL, ICBOX, IMOV, INDLAY

!     LOCAL LOGICAL SCALARS
LOGICAL :: LLINTRP

!     LOCAL REAL SCALARS
REAL_B :: ZASYMX, ZDIFFD, ZGI, ZGL, ZGR, ZIWGKG, ZLWGKG,&
          &ZMSAID, ZMSAIU, ZMSALD, ZMSALU, ZMTCONV, &
          &ZMTFUDG, ZLWFUDG, ZSWFUDG, ZMULTL, ZOI, ZOL, ZOMGMX, ZOR, &
          &ZRMUZ, ZRWGKG, ZTAUD, ZTAUMX, ZTEMPC, &
          &ZTOI, ZTOL, ZTOR, ZZFIWP, ZZFLWP, ZDPOG, ZPODT
REAL_B :: ZALND, ZASEA, ZD, ZDEN, ZNTOT, ZNUM, ZRATIO, ZCOEFF, Z1RADI,&
          &Z1RADL, ZBETAI, ZOMGI, ZOMGP, ZFDEL, ZWGHT, ZVI, ZVL, ZVR
REAL_B :: ZTOI, ZASW, ZOLR, ZSLW, ZSSW, ZMULTI, ZAIWC, ZBIWC,&
          &ZDICE, ZFSR, ZLGIWC, ZTCELS, ZTBLAY, ZADDPLK, ZPLANCK
REAL_B :: ZTOL1, ZTOI1, ZTOR1


!     -----------------------------------------------------------------

!if (NDUMP.LE.3) then
!  JL=KIDIA
!  DO jk=1,klev
!    print 9104,jk,PAPH(JL,JK),PTH(JL,JK),PAP(JL,JK),PT(JL,JK)&
!    &            ,PDP(JL,JK)& 
!    &            ,PQ(JL,JK),PFRCL(JL,JK),PQIWP(JL,JK),PQLWP(JL,JK)&
!    &            ,POZON(JL,JK),PQS(JL,JK)
9104 format(1x,i3,f9.1,f8.2,f9.1,f8.2,f9.1,e10.3,f7.4,4e10.3)
!  ENDDO 
!  jk=klev+1
!  print 9104,jk,PAPH(JL,JK),PTH(JL,JK)
!  print 9105,PTS(JL),(PALBD(JL,JSW),PALBP(JL,JSW),JSW=1,NSW)
9105 FORMAT(13X,f8.2,12f8.4)
!end if

!print *,'NICEOPT, NLIQOPT, NRADIP, NRADLP',NICEOPT,NLIQOPT,NRADIP,NRADLP

!-- compute total cloud cover
DO JL=KIDIA,KFDIA
  ZTCC(JL)=1.-PFRCL(JL,1)
  ZTCA(JL)=0.
END DO
DO JK=2,KLEV
  DO JL=KIDIA,KFDIA  
    ZTCC(JL)=ZTCC(JL)*(1.-MAX(PFRCL(JL,JK),PFRCL(JL,JK-1))) &
    & /(1.-MIN(PFRCL(JL,JK-1),1.-REPCLC))
  END DO
END DO
DO JL=KIDIA,KFDIA
  ZTCC(JL)=1.-ZTCC(JL)
END DO

!JL=KIDIA
!print 9106,ZTCC(JL)
9106 format(1x,'TCC :',F7.4)
!print 9107,LINHOM,NHOWINH
9107 format(1x,'LINHOM=',L8,' NHOWINH=',I2)
  





!*         1.     SET-UP INPUT QUANTITIES FOR RADIATION
!                 -------------------------------------

IF (.NOT.LINHOM) THEN
  ZMTFUDG=1.0_JPRB
  ZMTCONV=1.0_JPRB
  ZSWFUDG=1.0_JPRB
  ZLWFUDG=1.0_JPRB
ELSE IF (LINHOM) THEN
  IF (NHOWINH.EQ.1) THEN  
    ZMTFUDG=0.7_JPRB
    ZMTCONV=0.7_JPRB
    ZSWFUDG=0.7_JPRB
    ZLWFUDG=0.7_JPRB    
  ELSE
    ZMTFUDG=1.0_JPRB
    ZMTCONV=1.0_JPRB
    ZSWFUDG=1.0_JPRB
    ZLWFUDG=1.0_JPRB
  ENDIF    
ENDIF    
!print 9108,LINHOM,NHOWINH,ZSWFUDG
9108 format(1x,'LINHOM=',L8,' NHOWINH=',I2,' FUDG=',f4.2)

DO JL = KIDIA,KFDIA
  ZFCUP(JL,KLEV+1) = _ZERO_
  ZFCDWN(JL,KLEV+1) = REPLOG
  ZFSUP(JL,KLEV+1) = _ZERO_
  ZFSDWN(JL,KLEV+1) = REPLOG
  ZFLUX(JL,1,KLEV+1) = _ZERO_
  ZFLUX(JL,2,KLEV+1) = _ZERO_
  ZFLUC(JL,1,KLEV+1) = _ZERO_
  ZFLUC(JL,2,KLEV+1) = _ZERO_
  ZFSDNN(JL) = _ZERO_
  ZFSDNV(JL) = _ZERO_
  ZFCDNN(JL) = _ZERO_
  ZFCDNV(JL) = _ZERO_
  ZFSUPN(JL) = _ZERO_
  ZFSUPV(JL) = _ZERO_
  ZFCUPN(JL) = _ZERO_
  ZFCUPV(JL) = _ZERO_
  ZPSOL(JL) = PAPH(JL,KLEV+1)
  ZPMB(JL,1) = ZPSOL(JL) / 100._JPRB
  ZDT0(JL) = PTS(JL) - PTH(JL,KLEV+1)
  PSUDU(JL) = _ZERO_
  PUVDF(JL) = _ZERO_
  PPARF(JL) = _ZERO_
  ZSUDU(JL) = _ZERO_
  IBAS(JL) = INT ( 0.01_JPRB + PNBAS(JL) )
  ITOP(JL) = INT ( 0.01_JPRB + PNTOP(JL) )
ENDDO

DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    CPFLS(JL,JK)  = _ZERO_
    CPFLT(JL,JK)  = _ZERO_
    CPFCS(JL,JK)  = _ZERO_
    CPFCT(JL,JK)  = _ZERO_
    CPFDCT(JL,JK) = _ZERO_
    CPFUCT(JL,JK) = _ZERO_
    CPFDLT(JL,JK) = _ZERO_
    CPFULT(JL,JK) = _ZERO_
    CPFDCS(JL,JK) = _ZERO_
    CPFUCS(JL,JK) = _ZERO_
    CPFDLS(JL,JK) = _ZERO_
    CPFULS(JL,JK) = _ZERO_
  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  CPFRSOD(JL) = _ZERO_
  CPEMIT (JL) = _ZERO_ 
  CPSUDU (JL) = _ZERO_
  CPUVDF (JL) = _ZERO_
  CPPARF (JL) = _ZERO_
END DO   


!*         1.1    INITIALIZE VARIOUS FIELDS
!                 -------------------------


DO JSW=1,NSW
  DO JL = KIDIA,KFDIA
    ZALBD(JL,JSW)=PALBD(JL,JSW)
    ZALBP(JL,JSW)=PALBP(JL,JSW)
  ENDDO
ENDDO
DO JL = KIDIA,KFDIA
  ZEMIS(JL)  =PEMIS(JL)
  ZEMIW(JL)  =PEMIW(JL)
  ZMU0(JL)   =PMU0(JL)
  ZUVDF(JL)  = _ZERO_
  ZSUDU(JL)  = _ZERO_
  ZPARF(JL)  = _ZERO_
ENDDO

DO JK = 1 , KLEV
  JKP1 = JK + 1
  JKL = KLEV+ 1 - JK
  JKLP1 = JKL + 1
  DO JL = KIDIA,KFDIA
    ZPMB(JL,JK+1)=PAPH(JL,JKL)/100._JPRB
    ZOZ(JL,JK)   = POZON(JL,JKL) * 46.6968_JPRB / RG
    ZOZON(JL,JK) = POZON(JL,JKL)
    ZCLD0(JL,JK) = _ZERO_
    ZFCUP(JL,JK) = _ZERO_
    ZFCDWN(JL,JK) = _ZERO_
    ZFSUP(JL,JK) = _ZERO_
    ZFSDWN(JL,JK) = _ZERO_
    ZFLUX(JL,1,JK) = _ZERO_
    ZFLUX(JL,2,JK) = _ZERO_
    ZFLUC(JL,1,JK) = _ZERO_
    ZFLUC(JL,2,JK) = _ZERO_
  ENDDO
ENDDO


!** INPUTS ARE FULL LEVEL TEMPERATURES + SURFACE TEMPERATURE
!        INTERPOLATION TO GET HALF-LEVEL TEMPERATURES FOLLOWS
!        WHAT IS DONE IN *RADINT* AND *RADHEAT*

!* LLINTRP=.T.  Half-level temperatures on the coarse grid are
!               vertically interpolated linearly with horizontal
!               sampled pressure from the full-level temperatures
!               of the sampled grid.

!* LLINTRP=.F.  Half-level temperatures are those horizontally
!               sampled on the coarse grid

LLINTRP=.FALSE.
IF (LLINTRP) THEN
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      PTH(JL,JK)=(PT  (JL,JK-1)*PAP  (JL,JK-1)&
       &*(PAP  (JL,JK)-PAPH  (JL,JK))&
       &+PT  (JL,JK)*PAP  (JL,JK)*(PAPH  (JL,JK)-PAP  (JL,JK-1)))&
       &*(_ONE_/(PAPH  (JL,JK)*(PAP  (JL,JK)-PAP  (JL,JK-1))))
    ENDDO
  ENDDO
  IF (LTEMPDS) THEN
    DO JL=KIDIA,KFDIA
      PTH(JL,1)= PT  (JL,1)-PAP  (JL,1)*(PT  (JL,1)-PTH(JL,2))&
        &/(PAP  (JL,1)-PAPH  (JL,2)) 
      PTH(JL,KLEV+1)=PT(JL,KLEV)&
        &            +(PAPH(JL,KLEV+1)-PAP(JL,KLEV))&
        &            *(PT(JL,KLEV)-PTH(JL,KLEV))&
        &            /(PAP(JL,KLEV)-PAPH(JL,KLEV))
    ENDDO
  ELSE      
    DO JL=KIDIA,KFDIA
      PTH(JL,1)= PT  (JL,1)-PAP  (JL,1)*(PT  (JL,1)-PTH(JL,2))&
        &/(PAP  (JL,1)-PAPH  (JL,2)) 
      PTH(JL,KLEV+1)= PTS(JL)
    ENDDO
  ENDIF    
ENDIF

DO JK=1,KLEV
  JKL=KLEV+1-JK
  JKLP1=JKL+1
  DO JL=KIDIA,KFDIA
    ZTL(JL,JK)=PTH(JL,JKLP1)
    ZTAVE(JL,JK)=PT(JL,JKL)
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  ZTL(JL,KLEV+1)= PTH(JL,1)
  ZPMB(JL,KLEV+1) = PAPH(JL,1)/100._JPRB
ENDDO
!***

!     ------------------------------------------------------------------

!*         2.     CLOUD AND AEROSOL PARAMETERS
!                 ----------------------------

NBOXL=1
IF (KBOX.EQ.1) THEN
  CALL COL2BOX &
   & ( KIDIA, KFDIA, KLON, KLEV, NBOX, NOVLP &
   & , PFRCL, PCLBX &
   & )
  NBOXL=NBOX
END IF
ZWGHT=1./FLOAT(NBOXL) 
       
!-- initialise box-type outputs OLR, ASW, SDLW, SDSW, TAU
DO ICBOX=1,NBOXL
  DO JL=KIDIA,KFDIA
    OLRBOX(JL,ICBOX)=_ZERO_             
    ASWBOX(JL,ICBOX)=_ZERO_             
    SLWBOX(JL,ICBOX)=_ZERO_             
    SSWBOX(JL,ICBOX)=_ZERO_ 
    TAUBOX(JL,ICBOX)=_ZERO_ 
  END DO  
END DO             

DO ICBOX=1,NBOXL
  IF (KBOX.EQ.1) THEN
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PCLFR(JL,JK)=PCLBX(JL,ICBOX,JK)
      END DO
    END DO
    
  ELSE       
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PCLFR(JL,JK)=PFRCL(JL,JK)
      END DO
    END DO
  END IF  
  DO JL=KIDIA,KFDIA
    PSUDU(JL) = _ZERO_
    ZTAUINT(JL) = _ZERO_
  END DO  
  
!-- compute total cloud cover for that particular calculation
  DO JL=KIDIA,KFDIA
    ZCOL(JL)=1.-PCLFR(JL,1)
  END DO
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA  
      ZCOL(JL)=ZCOL(JL)*(1.-MAX(PCLFR(JL,JK),PCLFR(JL,JK-1))) &
       & /(1.-MIN(PCLFR(JL,JK-1),1.-REPCLC))
    END DO
  END DO
  DO JL=KIDIA,KFDIA
    ZCOL(JL)=1.-ZCOL(JL)
  END DO
  





DO JK = 1 , KLEV
  IKL = KLEV + 1 - JK

!          2.1    INITIALIZE OPTICAL PROPERTIES TO CLEAR SKY VALUES
!                 -------------------------------------------------

  DO JSW = 1,NSW
    DO JL = KIDIA,KFDIA
      ZTAU(JL,JSW,JK)  = _ZERO_
      ZOMEGA(JL,JSW,JK)= _ONE_
      ZCG(JL,JSW,JK)   = _ZERO_
    ENDDO
  ENDDO
  DO JL = KIDIA,KFDIA
    ZCLDSW(JL,JK)  = _ZERO_
    ZCLDLD(JL,JK)  = _ZERO_
    ZCLDLU(JL,JK)  = _ZERO_
  ENDDO


!          2.2    CLOUD ICE AND LIQUID CONTENT AND PATH
!                 -------------------------------------

  DO JL = KIDIA,KFDIA
    PCLFR(JL,IKL)=MAX( _ZERO_ ,MIN( PCLFR(JL,IKL), _ONE_ ))

! --- LIQUID WATER CONTENT (g.m-3) AND LIQUID WATER PATH (g.m-2)
    ZLWGKG=MAX(PQLWP(JL,IKL)*1000._JPRB,_ZERO_)
    ZIWGKG=MAX(PQIWP(JL,IKL)*1000._JPRB,_ZERO_)
    IF (PCLFR(JL,IKL) > REPCLC) THEN
      ZLWGKG=ZLWGKG/PFRCL(JL,IKL)
      ZIWGKG=ZIWGKG/PFRCL(JL,IKL)
    ELSE
      ZLWGKG=_ZERO_
      ZIWGKG=_ZERO_
    ENDIF

! --- RAIN LIQUID WATER CONTENT (g.m-3) AND LIQUID WATER PATH (g.m-2)
!    IF (PRAINT(JL,IKL).GT.(2.*REPCLC)) THEN
!      ZRWGKG=MAX(PQRAIN(JL,IKL)*1000., 0.0)
!      ZRAINT(JL)=PRAINT(JL,IKL)*3600.*1000.
!- no radiative effect of rain (for the moment)
!      ZRWGKG=0.
!      ZRAINT(JL)=0.
! ===========================================================

    ZRWGKG=_ZERO_
    ZRAINT(JL)=_ZERO_

    IF (IBAS(JL) /= 1.AND. ITOP(JL) /= 1 ) THEN
      ZCFUDG(JL)=ZMTCONV
    ELSE
      ZCFUDG(JL)=ZMTFUDG
    ENDIF

    ZDPOG=PDP(JL,IKL)/RG
    ZFLWP(JL)= ZLWGKG*ZDPOG
    ZFIWP(JL)= ZIWGKG*ZDPOG
    ZFRWP(JL)= ZRWGKG*ZDPOG
    ZPODT=PAP(JL,IKL)/(RD*PT(JL,IKL))
    ZLWC(JL)=ZLWGKG*ZPODT
    ZIWC(JL)=ZIWGKG*ZPODT
!    ZRWC(JL)=ZRWGKG*ZPODT

! --- EFFECTIVE RADIUS FOR WATER, ICE AND RAIN PARTICLES

    IF (NRADLP.EQ.0) THEN
! very old parametrization as f(pressure)
      ZRADLP(JL)=10._JPRB + (100000._JPRB-PAP(JL,IKL))*3.5E-04_JPRB

    ELSE IF (NRADLP.EQ.1) THEN
! simple distinction between land (10) and ocean (13)
      IF (PLSM(JL) < _HALF_) THEN
        ZRADLP(JL)=13._JPRB
      ELSE
        ZRADLP(JL)=10._JPRB
      ENDIF

    ELSE IF (NRADLP.EQ.2) THEN
!--  based on Martin et al., 1994, JAS
      IF (PLSM(JL) < _HALF_) THEN
        ZASEA=150._JPRB
        ZD=0.33_JPRB
        ZNTOT=-1.15E-03_JPRB*ZASEA*ZASEA+0.963_JPRB*ZASEA+5.30_JPRB
      ELSE
        ZALND=900._JPRB
        ZD=0.43_JPRB
        ZNTOT=-2.10E-04_JPRB*ZALND*ZALND+0.568_JPRB*ZALND-27.9_JPRB
      ENDIF
      
      ZNUM=3._JPRB*ZLWC(JL)*(1._JPRB+3._JPRB*ZD*ZD)**2
      ZDEN=4._JPRB*RPI*ZNTOT*(1._JPRB+ZD*ZD)**3
      ZRADLP(JL)=100.*(ZNUM/ZDEN)**0.333_JPRB
      
      ZRADLP(JL)=MAX(ZRADLP(JL), 4._JPRB)
      ZRADLP(JL)=MIN(ZRADLP(JL),16._JPRB)
    END IF  

! ===========================================================
! ___________________________________________________________

! rain drop from          : unused as ZRAINT is 0.
!    ZRADRD(JL)=500._JPRB*ZRAINT(JL)**0.22_JPRB
!    IF (ZFLWP(JL).GT.0.) THEN
!      ZRADRD(JL)=ZRADLP(JL)+ZRADRD(JL)
!    END IF   

  ENDDO
  DO JL = KIDIA,KFDIA

! diagnosing the ice particle effective radius/diameter

!- ice particle effective radius =f(T) from Liou and Ou (1994)
 
    IF (PT(JL,IKL) < RTICE) THEN
      ZTEMPC=PT(JL,IKL)-RTT
    ELSE
      ZTEMPC=RTICE-RTT
    ENDIF
    
    
    IF (NRADIP.EQ. 0) THEN
!-- fixed 40 micron effective radius
      ZRADIP(JL)= 40._JPRB
      ZDESR(JL)=2._JPRB*ZRADIP(JL)
      
    ELSE IF (NRADIP.EQ. 1) THEN 
      ZRADIP(JL)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
     &0.0012_JPRB))
      ZDESR(JL)=2._JPRB*ZRADIP(JL)
!-- old formulation based on Liou & Ou (1994) temperature (40-130microns)    
      ZRADIP(JL)=MAX(ZRADIP(JL),40._JPRB)
      ZDESR(JL)=2._JPRB*ZRADIP(JL)
      
    ELSE IF (NRADIP.EQ. 2) THEN  
      ZRADIP(JL)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
       &0.0012_JPRB))
      ZDESR(JL)=2._JPRB*ZRADIP(JL)
!-- formulation following Jakob, Klein modifications to ice content    
      ZRADIP(JL)=MAX(ZRADIP(JL),30._JPRB)
      ZRADIP(JL)=MIN(ZRADIP(JL),60._JPRB)
      ZDESR(JL)=2._JPRB*ZRADIP(JL)
 
 
    ELSE IF (NRADIP.EQ. 3  ) THEN
      ZRADIP(JL)=326.3_JPRB+ZTEMPC*(12.42_JPRB + ZTEMPC*(0.197_JPRB + ZTEMPC*&
       &0.0012_JPRB))
      ZDESR(JL)=2._JPRB*ZRADIP(JL)
!- ice particle effective radius =f(T,IWC) from Sun and Rikus (1999)
! revised by Sun (2001)
      IF (ZIWC(JL) > _ZERO_ ) THEN
        ZTEMPC = PT(JL,IKL)-83.15_JPRB
        ZTCELS = PT(JL,IKL)-RTT
        ZFSR = 1.2351_JPRB +0.0105_JPRB * ZTCELS
! Sun, 2001 (corrected from Sun & Rikus, 1999)
        ZAIWC = 45.8966_JPRB * ZIWC(JL)**0.2214_JPRB
        ZBIWC = 0.7957_JPRB * ZIWC(JL)**0.2535_JPRB
        ZDESR(JL) = ZFSR * (ZAIWC + ZBIWC*ZTEMPC)
        ZDESR(JL) = MIN ( MAX( ZDESR(JL), 45._JPRB), 350._JPRB)
        ZRADIP(JL)= 0.5 * ZDESR(JL)
      ELSE
        ZDESR(JL) = 80._JPRB
        ZRADIP(JL)= 0.5 * ZDESR(JL)
      END IF  
    END IF  
    
!-- ERA-15 definition of effective radii    
    IF (KLWRAD.EQ.2 .AND. NSW.EQ.2) THEN
      ZRADIP(JL)=40._JPRB
      ZRADLP(JL)=10._JPRB + (100000._JPRB-PAP(JL,IKL))*3.5_JPRB
      LOWASYF=.FALSE.      
      LOIFUEC=.FALSE.
      LRRTM=.FALSE.
      ZDESR(JL)=2._JPRB*ZRADIP(JL)
    END IF  
    
  ENDDO



!          2.3    CLOUD SHORTWAVE OPTICAL PROPERTIES
!                 ----------------------------------

!   -------------------------
! --+ SW OPTICAL PARAMETERS +  Water clouds after Fouquart (1987)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

  DO JSW=1,NSW
    DO JL = KIDIA,KFDIA
      ZTOL=_ZERO_
      ZGL =_ZERO_
      ZOL =_ZERO_
      ZTOI=_ZERO_
      ZGI =_ZERO_
      ZOI =_ZERO_
      ZTOR=_ZERO_
      ZGR =_ZERO_
      ZOR =_ZERO_
      IF (ZFLWP(JL)+ZFIWP(JL)+ZFRWP(JL) > REPSC ) THEN
        IF (ZFLWP(JL) > REPSC ) THEN
          IF (NLIQOPT.NE.0 ) THEN
!-- SW: Slingo, 1989
            ZTOL = ZFLWP(JL)*(RASWCA(JSW)+RASWCB(JSW)/ZRADLP(JL))
            ZGL  = RASWCE(JSW)+RASWCF(JSW)*ZRADLP(JL)
            ZOL  = 1. - RASWCC(JSW)-RASWCD(JSW)*ZRADLP(JL)
          ELSE          
!-- SW: Fouquart, 1991
            ZTOL = ZFLWP(JL)*(RYFWCA(JSW)+RYFWCB(JSW)/ZRADLP(JL))
            ZGL  = RYFWCF(JSW)
            ZOL  = RYFWCC(JSW)-RYFWCD(JSW)*EXP(-RYFWCE(JSW)*ZTOL)
          ENDIF 
        ENDIF

        IF (ZFIWP(JL) > REPSC ) THEN
          IF (NICEOPT.LE.1) THEN
!-- SW: Ebert-Curry          
            ZTOI = ZFIWP(JL)*(REBCUA(JSW)+REBCUB(JSW)/ZRADIP(JL))
            ZGI  = REBCUE(JSW)+REBCUF(JSW)*ZRADIP(JL)
            ZOI  = _ONE_ - REBCUC(JSW)-REBCUD(JSW)*ZRADIP(JL)
            
          ELSE IF (NICEOPT.EQ.2) THEN
!-- SW: Fu-Liou, 1993
            Z1RADI = 0.5 / ZRADIP(JL)
            ZBETAI = RFLAA0(JSW)+Z1RADI* RFLAA1(JSW)
            ZTOI = ZFIWP(JL) * ZBETAI
            ZOMGI= RFLBB0(JSW)+ZRADIP(JL)*(RFLBB1(JSW) + ZRADIP(JL) &
             &   *(RFLBB2(JSW)+ZRADIP(JL)* RFLBB3(JSW) ))            
            ZOI  = _ONE_ - ZOMGI
            ZOMGP= RFLCC0(JSW)+ZRADIP(JL)*(RFLCC1(JSW) + ZRADIP(JL) &
             &   *(RFLCC2(JSW)+ZRADIP(JL)* RFLCC3(JSW) )) 
            ZFDEL= RFLDD0(JSW)+ZRADIP(JL)*(RFLDD1(JSW) + ZRADIP(JL) &
             &   *(RFLDD2(JSW)+ZRADIP(JL)* RFLDD3(JSW) )) 
            ZGI  = ((1.-ZFDEL)*ZOMGP + ZFDEL*3.) / 3.    
                    
          ELSE IF (NICEOPT.EQ.3) THEN
!-- SW: Fu 1996
            Z1RADI = _ONE_ / ZDESR(JL)
            ZBETAI = RFUAA0(JSW)+Z1RADI* RFUAA1(JSW)
            ZTOI = ZFIWP(JL) * ZBETAI
            ZOMGI= RFUBB0(JSW)+ZDESR(JL)*(RFUBB1(JSW) + ZDESR(JL) &
             &   *(RFUBB2(JSW)+ZDESR(JL)* RFUBB3(JSW) ))            
            ZOI  = _ONE_ - ZOMGI
            ZGI  = RFUCC0(JSW)+ZDESR(JL)*(RFUCC1(JSW) + ZDESR(JL) &
             &   *(RFUCC2(JSW)+ZDESR(JL)* RFUCC3(JSW) )) 
             
          ENDIF
        ENDIF

!        IF (ZFRWP(JL) .NE. 0.) THEN
!          ZTOR= ZFRWP(JL)*0.003_JPRB*_JPRBZRAINT(JL)**(-0.22_JPRB)         
!          ZOR = 1._JPRB - RROMA(JSW)*ZRAINT(JL)**RROMB(JSW)
!          ZGR = RRASY(JSW)
!        END IF   

!  - MIX of WATER and ICE CLOUDS
!        ZTAUMX= ZTOL + ZTOI + ZTOR
!        ZOMGMX= ZTOL*ZOL + ZTOI*ZOI + ZTOR*ZOR
!        ZASYMX= ZTOL*ZOL*ZGL + ZTOI*ZOI*ZGI + ZTOR*ZOR*ZGR
!
!        ZASYMX= ZASYMX/ZOMGMX
!        ZOMGMX= ZOMGMX/ZTAUMX

        IF (.NOT.LINHOM .OR. (LINHOM .AND. NHOWINH.EQ.1) ) THEN
          ZVL=ZSWFUDG
          ZVI=ZSWFUDG
          ZVR=0.
          ZTAUMX= ZTOL*ZVL + ZTOI*ZVI + ZTOR*ZVR
          ZOMGMX= ZTOL*ZVL*ZOL + ZTOI*ZVI*ZOI + ZTOR*ZVR*ZOR
          ZASYMX= ZTOL*ZVL*ZOL*ZGL + ZTOI*ZVI*ZOI*ZGI + ZTOR*ZVR*ZOR*ZGR
          ZASYMX= ZASYMX/ZOMGMX
          ZOMGMX= ZOMGMX/ZTAUMX
        ELSE IF (LINHOM .AND. NHOWINH.EQ.2) THEN
          ZVL=PSQLW(JL,IKL)
          ZVI=PSQIW(JL,IKL)
          ZVR=0.
          ZTAUMX= ZTOL*ZVL + ZTOI*ZVI + ZTOR*ZVR
          ZOMGMX= ZTOL*ZVL*ZOL + ZTOI*ZVI*ZOI + ZTOR*ZVR*ZOR
          ZASYMX= ZTOL*ZVL*ZOL*ZGL + ZTOI*ZVI*ZOI*ZGI + ZTOR*ZVR*ZOR*ZGR
          ZASYMX= ZASYMX/ZOMGMX
          ZOMGMX= ZOMGMX/ZTAUMX
        ELSE IF (LINHOM .AND. NHOWINH.EQ.3) THEN
          ZVL=PRLVRL(JL,IKL)
          ZVI=PRLVRI(JL,IKL)
          ZVR=0.
          ZTOL1 = ZTOL/(1.+ZVL)
          ZTOI1 = ZTOI/(1.+ZVI)
          ZTOR1 = ZTOR/(1.+ZVR) 
          ZTAUMX= ZTOL1 + ZTOI1 + ZTOR1
          ZOI=ZOI/(1.+ZVI*(1.-ZOI))
          ZGI=ZGI*(1.+ZVI*(1.-ZOI))/(1.+ZVI*(1.-ZOI*ZGI))
          ZOL=ZOL/(1.+ZVL*(1.-ZOL))
          ZGL=ZGL*(1.+ZVL*(1.-ZOL))/(1.+ZVL*(1.-ZOL*ZGL))
          
          ZOMGMX= ZTOL1*ZOL + ZTOI1*ZOI + ZTOR1*ZOR
          ZASYMX= ZTOL1*ZOL*ZGL + ZTOI1*ZOI*ZGI + ZTOR1*ZOR*ZGR
          ZASYMX= ZASYMX/ZOMGMX
          ZOMGMX= ZOMGMX/ZTAUMX
        END IF  
9009    format(1x,3I3,14E13.6)         
        
! --- SW FINAL CLOUD OPTICAL PARAMETERS

        ZCLDSW(JL,JK)  = PCLFR(JL,IKL)
        ZTAU(JL,JSW,JK)  = ZTAUMX
        ZOMEGA(JL,JSW,JK)= ZOMGMX
        ZCG(JL,JSW,JK)   = ZASYMX
      ENDIF
    ENDDO
  ENDDO
  
  DO JL=KIDIA,KFDIA
    ZTAUINT(JL)=ZTAUINT(JL)+ZTAU(JL,1,JK)
  END DO  


!          2.4    CLOUD LONGWAVE OPTICAL PROPERTIES FOR EC-OPE
!                 --------------------------------------------

!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Smith and Shi (1992)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

  IF (.NOT.LRRTM) THEN

    DO JL = KIDIA,KFDIA
      ZALFICE(JL)=_ZERO_
      ZGAMICE(JL)=_ZERO_
      ZBICE(JL)=_ZERO_
      ZTICE(JL)=(PT(JL,IKL)-TSTAND)/TSTAND
      ZBICFU(JL)=_ZERO_
      ZKICFU1(JL)=_ZERO_
      ZKICFU2(JL)=_ZERO_
    ENDDO
    
    DO JNU= 1,NSIL
      DO JL = KIDIA,KFDIA
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
        
!-- Fu et al. (1998) with M'91 LW scheme    
    DO JRTM=1,16
      DO JL=KIDIA,KFDIA
        IF (PT(JL,IKL) < 339._JPRB .AND. PT(JL,IKL) >= 160._JPRB) THEN
          INDLAY=PT(JL,IKL)-159._JPRB
          ZTBLAY =PT(JL,IKL)-INT(PT(JL,IKL))
        ELSE IF (PT(JL,IKL) >= 339._JPRB ) THEN
          INDLAY=180
          ZTBLAY =PT(JL,IKL)-339._JPRB
        ELSE IF (PT(JL,IKL) < 160._JPRB) THEN
          INDLAY=1
          ZTBLAY =PT(JL,IKL)-160._JPRB
        END IF      
        ZADDPLK = TOTPLNK(INDLAY+1,JRTM)-TOTPLNK(INDLAY,JRTM)
        ZPLANCK = DELWAVE(JRTM) * (TOTPLNK(INDLAY,JRTM) + ZTBLAY*ZADDPLK)
        ZBICFU(JL) = ZBICFU(JL) + ZPLANCK
        
        IF (ZIWC(JL) > _ZERO_ ) THEN
! ice cloud spectral emissivity a la Fu & Liou (1993)
          ZRATIO= 0.5 / ZRADIP(JL)
          ZMSAID = RFULIO(JRTM,1) + ZRATIO&
             &*(RFULIO(JRTM,2) + ZRATIO*RFULIO(JRTM,3))
          ZKICFU1(JL) = ZKICFU1(JL)+ ZMSAID*ZPLANCK
          
! ice cloud spectral emissivity a la Fu et al (1998)
          Z1RADI = _ONE_ / ZDESR(JL)
          ZMSAID = RFUETA(JRTM,1) + Z1RADI&
             &*(RFUETA(JRTM,2) + Z1RADI*RFUETA(JRTM,3))
          ZKICFU2(JL) = ZKICFU2(JL)+ ZMSAID*ZPLANCK
        END IF  
      END DO
    END DO
            
    DO JL = KIDIA,KFDIA
      ZGAMICE(JL) = ZGAMICE(JL) / ZBICE(JL)
      ZALFICE(JL) = ZALFICE(JL) / ZBICE(JL)
      ZKICFU1(JL) = ZKICFU1(JL) / ZBICFU(JL)
      ZKICFU2(JL) = ZKICFU2(JL) / ZBICFU(JL)
      
      IF (ZFLWP(JL)+ZFIWP(JL) /= _ZERO_) THEN

        IF (KLWRAD.EQ.2) THEN        
! ice cloud emissivity a la Smith-Shi
          ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIP(JL)
          ZMSAID= 0.113_JPRB*ZMULTI
          ZMSAIU= 0.093_JPRB*ZMULTI
          ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
          ZMSALD= 0.158_JPRB*ZMULTL
          ZMSALU= 0.130_JPRB*ZMULTL
          ZZFLWP= ZFLWP(JL)
          ZZFIWP= ZFIWP(JL)
          
        ELSE IF (KLWRAD.EQ.0) THEN  
          
          IF (NLIQOPT.EQ.0) THEN
! water cloud emissivity a la Smith & Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
            ZMSALD= 0.158_JPRB*ZMULTL
            ZMSALU= 0.130_JPRB*ZMULTL
          
          ELSE
! water cloud emissivity a la Savijarvi (1997)
            ZMSALU= 0.2441_JPRB-0.0105_JPRB*ZRADLP(JL)
            ZMSALD= 1.2154_JPRB*ZMSALU
          
          END IF  
          
          IF (NICEOPT.EQ.0) THEN          
! ice cloud emissivity a la Smith & Shi (1992)
            ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIP(JL)
            ZMSAID= 0.113_JPRB*ZMULTI
            ZMSAIU= 0.093_JPRB*ZMULTI

          ELSE IF (NICEOPT.EQ.1) THEN
! ice cloud emissivity a la Ebert & Curry (1992)
            ZMSAID= 1.66_JPRB*(ZALFICE(JL)+ZGAMICE(JL)/ZRADIP(JL))
            ZMSAIU= ZMSAID
         
          ELSE IF (NICEOPT.EQ.2) THEN  
! ice cloud emissivity a la Fu & Liou (1993)
            ZMSAID= 1.66_JPRB*ZKICFU1(JL)
            ZMSAIU= ZMSAID
          
          ELSE IF (NICEOPT.EQ.3) THEN  
! ice cloud emissivity a la Fu et al. (1998)
            ZMSAID= 1.66_JPRB*ZKICFU2(JL)
            ZMSAIU= ZMSAID
          END IF  
          
! introduce inhomogeneity factor also in LW         
          ZZFLWP= ZFLWP(JL) * ZLWFUDG
          ZZFIWP= ZFIWP(JL) * ZLWFUDG
        END IF
          
! effective cloudiness accounting for condensed water
        ZCLDLD(JL,JK) = PCLFR(JL,IKL)*(_ONE_-EXP(-ZMSALD*ZZFLWP-ZMSAID* &
          &ZZFIWP))
        ZCLDLU(JL,JK) = PCLFR(JL,IKL)*(_ONE_-EXP(-ZMSALU*ZZFLWP-ZMSAIU* &
          &ZZFIWP))
          
      END IF    
    ENDDO
    
  ELSE

!          2.5    CLOUD LONGWAVE OPTICAL PROPERTIES FOR RRTM
!                 ------------------------------------------

!   -------------------------
! --+ LW OPTICAL PARAMETERS +  Water (and Ice) from Savijarvi (1998)
!   -------------------------  Ice clouds (Ebert, Curry, 1992)

! No need for a fixed diffusivity factor, accounted for spectrally below
! The detailed spectral structure does not require defining upward and
! downward effective optical properties

    DO JRTM=1,16
      DO JL = KIDIA,KFDIA
        ZTAUCLD(JL,JK,JRTM) = _ZERO_
        ZMSALD = _ZERO_
        ZMSAID = _ZERO_
        
        IF (ZFLWP(JL)+ZFIWP(JL) /= _ZERO_) THEN
    
          IF (NLIQOPT.EQ.0) THEN
! water cloud total emissivity a la Smith and Shi (1992)
            ZMULTL=1.2_JPRB-0.006_JPRB*ZRADLP(JL)
            ZMSALD= 0.144_JPRB*ZMULTL / 1.66_JPRB
            
          ELSE IF (NLIQOPT.EQ.1) THEN
! water cloud spectral emissivity a la Savijarvi (1997)
            ZMSALD= RHSAVI(JRTM,1) + ZRADLP(JL)&
             &*(RHSAVI(JRTM,2) + ZRADLP(JL)*RHSAVI(JRTM,3))
             
          ELSE IF (NLIQOPT.EQ.2) THEN
! water cloud spectral emissivity a la Lindner and Li (2000)
            Z1RADL = _ONE_ / ZRADLP(JL)
            ZMSALD = RLINLI(JRTM,1)+ZRADLP(JL)*RLINLI(JRTM,2)+ Z1RADL*&
            &       (RLINLI(JRTM,3) + Z1RADL*(RLINLI(JRTM,4) + Z1RADL*&
            &        RLINLI(JRTM,5) ))
          
          END IF  

          IF (NICEOPT.EQ.0) THEN
! ice cloud emissivity a la Smith & Shi (1992)
            ZMULTI=1.2_JPRB-0.006_JPRB*ZRADIP(JL)
            ZMSAID= 0.108_JPRB*ZMULTI / 1.66_JPRB
                   
          ELSE IF (NICEOPT.EQ.1) THEN
! ice cloud spectral emissivity a la Ebert-Curry (1992)
            ZMSAID= REBCUH(JRTM)+REBCUG(JRTM)/ZRADIP(JL)
            
          ELSE IF (NICEOPT.EQ.2) THEN
! ice cloud spectral emissivity a la Fu & Liou (1993)
            ZRATIO= 0.5 / ZRADIP(JL)
            ZMSAID = RFULIO(JRTM,1) + ZRATIO&
             &*(RFULIO(JRTM,2) + ZRATIO*RFULIO(JRTM,3))
             
          ELSE IF (NICEOPT.EQ.3) THEN
! ice cloud spectral emissivity a la Fu et al (1998)
            Z1RADI = _ONE_ / ZDESR(JL)
            ZMSAID = RFUETA(JRTM,1) + Z1RADI&
             &*(RFUETA(JRTM,2) + Z1RADI*RFUETA(JRTM,3))
             
          END IF    

          IF (.NOT.LINHOM .OR. (LINHOM .AND. NHOWINH.EQ.1) ) THEN
            ZVL=ZLWFUDG
            ZVI=ZLWFUDG
          ELSE IF (LINHOM .AND. NHOWINH.EQ.2) THEN
            ZVL=PSQLW(JL,IKL)
            ZVI=PSQIW(JL,IKL)
          ELSE IF (LINHOM .AND. NHOWINH.EQ.3) THEN
            ZVL=_ONE_/(_ONE_+PRLVRL(JL,IKL))
            ZVI=_ONE_/(_ONE_+PRLVRI(JL,IKL))
          END IF  
          
          ZTAUD = ZVL*ZMSALD*ZFLWP(JL)+ZVI*ZMSAID*ZFIWP(JL)

! Diffusivity correction within clouds a la Savijarvi
!          ZDIFFD=MIN(MAX(1.517_JPRB-0.156_JPRB*LOG(ZTAUD) , _ONE_) , _TWO_)

          ZDIFFD=1.66_JPRB
          ZTAUCLD(JL,JK,JRTM) = ZTAUD*ZDIFFD
        ENDIF
        
      ENDDO
    ENDDO
!  print *,'Radlsw after LW1 cloud optical properties for level JK=',JK

  ENDIF

ENDDO

NUAER = NUA
NTRAER = NTRA

!     ------------------------------------------------------------------

!*         2.6    DIFFUSIVITY FACTOR OR SATELLITE VIEWING ANGLE
!                 ---------------------------------------------


DO JL = KIDIA,KFDIA
  ZVIEW(JL) = DIFF
  ZEMIT(JL) = _ZERO_
ENDDO

!     ------------------------------------------------------------------

!*         3.     CALL LONGWAVE RADIATION CODE
!                 ----------------------------


!*         3.1    FULL LONGWAVE RADIATION COMPUTATIONS
!                 ------------------------------------

IF (.NOT.LPHYLIN) THEN
  IF ( .NOT. LRRTM) THEN

      
    IF (KLWRAD .EQ. 2) THEN
      CALL OLW &
       & ( KIDIA, KFDIA , KLON  , KLEV &
       & , PCCO2, ZCLDLD, ZCLDLU &
       & , PDP  , ZDT0  , ZEMIS  &
       & , PAPH , POZON , PTH &
       & , PAER , PT    , ZVIEW , PQ &
       & , ZCOOLR,ZCOOLC, ZFLUX, ZFLUC &
       & )
       
    ELSE IF (KLWRAD .EQ. 0) THEN  
     
      CALL LW &
       &( KIDIA , KFDIA , KLON  , KLEV , KMODE &
       &, PCCO2 , ZCLDLD, ZCLDLU &
       &, PDP   , ZDT0  , ZEMIS , ZEMIW &
       &, ZPMB  , POZON , ZTL &
       &, PAER  , ZTAVE , ZVIEW , PQ &
       &, ZCOOLR, ZCOOLC, ZEMIT , ZFLUX, ZFLUC &
       &)
       
     END IF  

  ELSE


!*         3.2    FULL LONGWAVE RADIATION COMPUTATIONS - RRTM
!                 ------------------------------------   ----

!  i)  pass POZN (ozone mmr concentration) to RRTM; remove pressure
!      weighting applied to POZON in driverMC (below)
!  ii) pass ZEMIS and ZEMIW to RRTM; return ZEMIT from RRTM
!  iii)pass ZTAUCLD, cloud optical depths (water+ice) to RRTM, 
!      computed from equations above
!  iv) pass ECRT arrays to RRTM arrays in interface routine ECRTATM
!      in module rrtm_ecrt.f

    DO JL = KIDIA,KFDIA
      DO JK = 1, KLEV
        ZOZN(JL,JK) = POZON(JL,JK)/PDP(JL,JK)
      ENDDO
    ENDDO

!    print *,'Just before calling RRTM'

    CALL RRTM_RRTM_140GP &
     &( KIDIA , KFDIA , KLON  , KLEV &
     &, PAER  , PAPH  , PAP   &
     &, PTS   , PTH   , PT     &
     &, ZEMIS , ZEMIW &
     &, PQ    , PCCO2 , ZOZN  , ZCLDSW  , ZTAUCLD &
     &, ZEMIT , ZFLUX , ZFLUC , ZTCLEAR &
     &)
     
!     print *,'just after RRTM'    

  ENDIF
ELSE
  ZCOOLR(:,:) = _ZERO_
  ZCOOLC(:,:) = _ZERO_
  ZEMIT (:)   = _ZERO_
  ZFLUX(:,:,:)= _ZERO_
  ZFLUC(:,:,:)= _ZERO_
ENDIF

!     ------------------------------------------------------------------

!*         4.     CALL SHORTWAVE RADIATION CODE
!                 -----------------------------


ZRMUZ=_ZERO_
DO JL = KIDIA,KFDIA
  ZRMUZ = MAX (ZRMUZ, ZMU0(JL))
ENDDO

IF (ZRMUZ > _ZERO_) THEN

  CALL SW &
   &( KIDIA , KFDIA , KLON  , KLEV  , KAER &
   &, PRII0 , PCCO2 , ZPSOL , ZALBD , ZALBP , PQ   , PQS &
   &, ZMU0  , ZCG   , ZCLDSW, PDP   , ZOMEGA, ZOZ  , ZPMB &
   &, ZTAU  , ZTAVE , PAER &
   &, ZHEATR, ZFSDWN, ZFSUP , ZHEATC, ZFCDWN, ZFCUP &
   &, ZFSDNN, ZFSDNV, ZFSUPN, ZFSUPV &
   &, ZFCDNN, ZFCDNV, ZFCUPN, ZFCUPV &
   &, ZSUDU , ZUVDF , ZPARF &
   &)
   
ENDIF

!     ------------------------------------------------------------------

!*         5.     FILL UP THE MODEL NET LW AND SW RADIATIVE FLUXES
!                 ------------------------------------------------


DO JKL = 1 , KLEV+1
  JK = KLEV+1 + 1 - JKL
  DO JL = KIDIA,KFDIA
  
    CPFLS(JL,JKL) =CPFLS(JL,JKL) +ZWGHT*(ZFSDWN(JL,JK) - ZFSUP(JL,JK))
    CPFLT(JL,JKL) =CPFLT(JL,JKL) +ZWGHT*(- ZFLUX(JL,1,JK) - ZFLUX(JL,2,JK))
    CPFCS(JL,JKL) =CPFCS(JL,JKL) +ZWGHT*(ZFCDWN(JL,JK)  - ZFCUP(JL,JK))
    CPFCT(JL,JKL) =CPFCT(JL,JKL) +ZWGHT*(- ZFLUC(JL,1,JK) - ZFLUC(JL,2,JK))
    CPFDCT(JL,JKL)=CPFDCT(JL,JKL)+ZWGHT*ZFLUC(JL,2,JK)
    CPFUCT(JL,JKL)=CPFUCT(JL,JKL)+ZWGHT*ZFLUC(JL,1,JK)
    CPFDLT(JL,JKL)=CPFDLT(JL,JKL)+ZWGHT*ZFLUX(JL,2,JK)
    CPFULT(JL,JKL)=CPFULT(JL,JKL)+ZWGHT*ZFLUX(JL,1,JK)
    CPFDCS(JL,JKL)=CPFDCS(JL,JKL)+ZWGHT*ZFCDWN(JL,JK)
    CPFUCS(JL,JKL)=CPFUCS(JL,JKL)+ZWGHT*ZFCUP(JL,JK)
    CPFDLS(JL,JKL)=CPFDLS(JL,JKL)+ZWGHT*ZFSDWN(JL,JK)
    CPFULS(JL,JKL)=CPFULS(JL,JKL)+ZWGHT*ZFSUP(JL,JK)
  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  print 9507,ZFSDWN(JL,1),ZSUDU(JL),ZUVDF(JL),ZPARF(JL)
9507 format(1x,'SW Global Normal UV & PAR:' 5f10.3)
  
  CPFRSOD(JL) = CPFRSOD(JL) + ZWGHT*ZFSDWN(JL,1)
  CPEMIT (JL) = CPEMIT (JL) + ZWGHT*ZEMIT (JL)
  CPSUDU (JL) = CPSUDU (JL) + ZWGHT*ZSUDU (JL)
  CPUVDF (JL) = CPUVDF (JL) + ZWGHT*ZUVDF (JL)
  CPPARF (JL) = CPPARF (JL) + ZWGHT*ZPARF (JL)
  
  ASWBOX(JL,ICBOX) = -ZFSDWN(JL,KLEV+1) + ZFSUP(JL,KLEV+1)
  OLRBOX(JL,ICBOX) = -ZFLUX(JL,1,KLEV+1)
  SLWBOX(JL,ICBOX) = -ZFLUX(JL,2,1)
  SSWBOX(JL,ICBOX) = -ZFSDWN(JL,1)
  TAUBOX(JL,ICBOX) = ZTAUINT(JL)
  ZTCA(JL) = ZTCA(JL) + ZWGHT*ZCOL(JL)
9508 format(1x,'radlsw',I3,5F10.3,1x,3F7.4)  
ENDDO


ENDDO
!
!-- end of box-type calculations
!      
 
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFLS(JL,JK)  = CPFLS(JL,JK) 
    PFLT(JL,JK)  = CPFLT(JL,JK) 
    PFCS(JL,JK)  = CPFCS(JL,JK) 
    PFCT(JL,JK)  = CPFCT(JL,JK) 
    PFDCT(JL,JK) = CPFDCT(JL,JK)
    PFUCT(JL,JK) = CPFUCT(JL,JK)
    PFDLT(JL,JK) = CPFDLT(JL,JK)
    PFULT(JL,JK) = CPFULT(JL,JK)
    PFDCS(JL,JK) = CPFDCS(JL,JK)
    PFUCS(JL,JK) = CPFUCS(JL,JK)
    PFDLS(JL,JK) = CPFDLS(JL,JK)
    PFULS(JL,JK) = CPFULS(JL,JK)
  ENDDO
ENDDO

DO JL = KIDIA,KFDIA
  PFRSOD(JL) = CPFRSOD(JL)
  PEMIT (JL) = CPEMIT (JL)
  PSUDU (JL) = CPSUDU (JL)
  PUVDF (JL) = CPUVDF (JL)
  PPARF (JL) = CPPARF (JL)
ENDDO

!-- re-organize the box-tyoe output arrays in decreasing order of TAU
DO JL=KIDIA,KFDIA
  DO ICBOX=2,NBOX
    ZTOI=TAUBOX(JL,ICBOX)
    DO IMOV=ICBOX-1,1,-1
      IF(TAUBOX(JL,IMOV).LE.ZTOI) GO TO 8001
        TAUBOX(JL,IMOV+1)=TAUBOX(JL,IMOV)
    END DO
    IMOV=0
8001 CONTINUE
    TAUBOX(JL,IMOV+1)=ZTOI
  END DO  
END DO

!-- re-organize the box-type output arrays in decreasing order of ASW
DO JL=KIDIA,KFDIA
  DO ICBOX=2,NBOX
    ZASW=ASWBOX(JL,ICBOX)
    DO IMOV=ICBOX-1,1,-1
      IF(ASWBOX(JL,IMOV).LE.ZASW) GO TO 8002
        ASWBOX(JL,IMOV+1)=ASWBOX(JL,IMOV)
    END DO
    IMOV=0
8002 CONTINUE
    ASWBOX(JL,IMOV+1)=ZASW
  END DO  
END DO

!-- re-organize the box-tyoe output arrays in decreasing order of -OLR
DO JL=KIDIA,KFDIA
  DO ICBOX=2,NBOX
    ZOLR=OLRBOX(JL,ICBOX)
    DO IMOV=ICBOX-1,1,-1
      IF(OLRBOX(JL,IMOV).LE.ZOLR) GO TO 8003
        OLRBOX(JL,IMOV+1)=OLRBOX(JL,IMOV)
    END DO
    IMOV=0
8003 CONTINUE
    OLRBOX(JL,IMOV+1)=ZOLR
  END DO  
END DO

!-- re-organize the box-tyoe output arrays in decreasing order of SLW
DO JL=KIDIA,KFDIA
  DO ICBOX=2,NBOX
    ZSLW=SLWBOX(JL,ICBOX)
    DO IMOV=ICBOX-1,1,-1
      IF(SLWBOX(JL,IMOV).LE.ZSLW) GO TO 8004
        SLWBOX(JL,IMOV+1)=SLWBOX(JL,IMOV)
    END DO
    IMOV=0
8004 CONTINUE
    SLWBOX(JL,IMOV+1)=ZSLW
  END DO  
END DO

!-- re-organize the box-type output arrays in decreasing order of -SSW
DO JL=KIDIA,KFDIA
  DO ICBOX=2,NBOX
    ZSSW=SSWBOX(JL,ICBOX)
    DO IMOV=ICBOX-1,1,-1
      IF(SSWBOX(JL,IMOV).LE.ZSSW) GO TO 8005
        SSWBOX(JL,IMOV+1)=SSWBOX(JL,IMOV)
    END DO
    IMOV=0
8005 CONTINUE
    SSWBOX(JL,IMOV+1)=ZSSW
  END DO  
END DO

!-- put all arrays as positive numbers for plotting
DO JL=KIDIA,KFDIA
  DO ICBOX=1,NBOX
    ASWBOX(JL,ICBOX)=-ASWBOX(JL,ICBOX)
    OLRBOX(JL,ICBOX)=-OLRBOX(JL,ICBOX)
    SSWBOX(JL,ICBOX)=-SSWBOX(JL,ICBOX)
  END DO  
END DO

!     --------------------------------------------------------------

RETURN
END SUBROUTINE RADLSW
