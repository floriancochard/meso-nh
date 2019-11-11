!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
PROGRAM RAD1DRIV

#include "tsmbkind.h"

USE PARRTM1D , ONLY : JP_LON   ,JP_IDIA  ,JP_FDIA  ,JP_TDIA  ,&
 &            JP_LEV ,JP_LW    ,JP_SW    ,JP_NUA   ,JP_MODE  ,&
 &            JP_AER ,JP_LEVP1

!     USE YOMLUN   , ONLY : NULNAM

USE OYOMCST   , ONLY : RD       ,RG       ,RTT      ,RSIGMA   ,&
 &            RCPD   ,RPI      ,RDAY     ,RCPD     ,REA      ,&
 &            RI0    ,RSIGMA   ,REPSM
USE YOEAERD  , ONLY : CVDAES   ,CVDAEL   ,CVDAEU   ,CVDAED   ,&
 &            RCAEOPS  ,RCAEOPL  ,RCAEOPU  ,RCAEOPD  ,RCTRBGA  ,&
 &            RCVOBGA  ,RCSTBGA  ,RCTRPT   ,RCAEADM  ,RCAEROS  ,&
 &            RCAEADK
USE YOELW    , ONLY : NSIL     ,NIPD     ,NTRA     ,NUA      ,&
 &            NG1      ,NG1P1    ,WG1
USE YOEOVLP  , ONLY : RA1OVLP
USE YOEPHLI  , ONLY : LPHYLIN
USE OYOERAD   , ONLY : NAER     ,NMODE    ,NOZOCL   ,&
 &            NRADFR   ,NRADPFR  ,NRADPLA  ,NRINT    ,NHOWINH  ,&
 &            NOVLP    ,NRADF2C  ,NRADC2F  ,NLW      ,NSW      ,&
 &            NTSW     ,LERAD6H  ,LERADHS  ,LHVOLCA  ,LNEWAER  ,&
 &            LONEWSW  ,LOWASYF  ,LOWHSSS  ,LOIFUEC  ,LRRTM    ,&
 &            LRADLP   ,LINHOM   ,RAOVLP   ,RBOVLP   ,&
 &            NICEOPT  ,NLIQOPT  ,NRADIP   ,NRADLP   ,RMINICE
USE OYOERDI   , ONLY : RCARDI   ,RCH4     ,RN2O     ,RO3      ,&
 &            RCFC11   ,RCFC12   ,REPCLC  
USE YOERDU   , ONLY : NUAER    ,NTRAER   ,NIMP     ,NOUT     ,&
 &            RCDAY    ,R10E     ,REPLOG   ,REPSC    ,REPSCO   ,&
 &            REPSCQ   ,REPSCT   ,REPSCW   ,DIFF
USE OYOESW    , ONLY : LO3ONLY 
USE YOETHF   , ONLY : R2ES     ,R3LES    ,R3IES    ,R4LES    ,&
 &            R4IES    ,R5LES    ,R5IES    ,R5ALVCP  ,R5ALSCP  ,&
 &            RALVDCP  ,RALSDCP  ,RTWAT    ,RTICE    ,RTICECU
USE YOEDBUG  , ONLY : LDEBUG 

     
!IMPLICIT NONE

!     ------------------------------------------------------------------
!
CHARACTER*132 CBLABLA
CHARACTER*80 CPROF, CSURF
CHARACTER*75 CFILOUT0,CFILOUT1,CFILOUT2,CFILOUT3,CFILOUT4,CFILOUT5,CONFIGUR
CHARACTER*75 CFILOUT6
CHARACTER*40 CTITLE(JP_LON)
CHARACTER*20 CFILE
CHARACTER*7  COP
CHARACTER*3  CBOX, CSTAT, CID(42), CRMINICE, CIRCC
     
!- reference arrays
             
REAL_B :: PAPH0(JP_LON,JP_LEVP1), PAP0(JP_LON,JP_LEV)
REAL_B :: PAPH5(JP_LON,JP_LEVP1), PAP5(JP_LON,JP_LEV)
REAL_B :: PTH0(JP_LON,JP_LEVP1) , PT0(JP_LON,JP_LEV)
REAL_B :: PTH5(JP_LON,JP_LEVP1) , PT5(JP_LON,JP_LEV) , PTS5(JP_LON) 
REAL_B :: PQ0(JP_LON,JP_LEV)    , PRH0(JP_LON,JP_LEV), PQS0(JP_LON,JP_LEV)
REAL_B :: PQ5(JP_LON,JP_LEV)    , PRH5(JP_LON,JP_LEV), PQS5(JP_LON,JP_LEV)

REAL_B :: PAER5(JP_LON,JP_AER,JP_LEV)
REAL_B :: PALBD5(JP_LON,JP_SW)  , PALBP5(JP_LON,JP_SW)
REAL_B :: PCLFR5(JP_LON,JP_LEV) , PDP5(JP_LON,JP_LEV)
REAL_B :: PEMIS5(JP_LON), PEMIW5(JP_LON), PLSM5(JP_LON), PMU05(JP_LON)
REAL_B :: PGELAM5(JP_LON),PGEMU5(JP_LON), PSLON5(JP_LON),PCLON5(JP_LON)
REAL_B :: POZON5(JP_LON,JP_LEV)
REAL_B :: PQIWP5(JP_LON,JP_LEV) , PQLWP5(JP_LON,JP_LEV)
REAL_B :: PSQIW5(JP_LON,JP_LEV) , PSQLW5(JP_LON,JP_LEV)
REAL_B :: PRLVRI5(JP_LON,JP_LEV), PRLVRL5(JP_LON,JP_LEV)
REAL_B :: PQRAIN5(JP_LON,JP_LEV)
REAL_B :: PRAINT5(JP_LON,JP_LEV)
     
REAL_B :: PEMIT5(JP_LON)
REAL_B :: PFCT5(JP_LON,JP_LEVP1), PFLT5(JP_LON,JP_LEVP1)
REAL_B :: PFCS5(JP_LON,JP_LEVP1), PFLS5(JP_LON,JP_LEVP1)
REAL_B :: PFDCT5(JP_LON,JP_LEVP1), PFDLT5(JP_LON,JP_LEVP1)
REAL_B :: PFDCS5(JP_LON,JP_LEVP1), PFDLS5(JP_LON,JP_LEVP1)
REAL_B :: PFUCT5(JP_LON,JP_LEVP1), PFULT5(JP_LON,JP_LEVP1)
REAL_B :: PFUCS5(JP_LON,JP_LEVP1), PFULS5(JP_LON,JP_LEVP1)
REAL_B :: PFRSOD5(JP_LON), PSUDU5(JP_LON), PUVDF5(JP_LON), PPARF5(JP_LON)
REAL_B :: PNBAS5(JP_LON)        , PNTOP5(JP_LON)
REAL_B :: PTCC5(JP_LON)
      
!- extra arrays             
   
REAL_B :: PLAT5(JP_LON), PLON5(JP_LON), ZOZON5(JP_LON,JP_LEV)
REAL_B :: ZETA(JP_LEV),ZETAH(JP_LEVP1), ZAER5(JP_LON,JP_AER,JP_LEV)
REAL_B :: CLOUD(JP_LEV)
REAL_B :: CLOUD1(JP_LEV),CLOUD2(JP_LEV),CLOUD3(JP_LEV),CLOUD4(JP_LEV), &
  & CLOUD5(JP_LEV)

REAL_B :: ZAPH(JP_LEVP1), ZAP(JP_LEV), ZATH(JP_LEVP1), ZAT(JP_LEV)
REAL_B :: ZAOH(JP_LEVP1), ZAO(JP_LEV), ZAQ(JP_LEV)
REAL_B :: ZACFI(JP_LEV) , ZACFL(JP_LEV), ZACF(JP_LEV)  
REAL_B :: ZACQI(JP_LEV) , ZACQL(JP_LEV), ZACQ(JP_LEV)
REAL_B :: ZSCQL(JP_LEV) , ZSCQI(JP_LEV), ZSCQ(JP_LEV)
REAL_B :: AAM(JP_LEVP1), BBM(JP_LEVP1)
REAL_B :: ZAZH(JP_LEVP1), ZAZ(JP_LEV)

!- box-type results

REAL_B :: ASWBOX(JP_LON,100), OLRBOX(JP_LON,100)
REAL_B :: SLWBOX(JP_LON,100), SSWBOX(JP_LON,100)
REAL_B :: TAUBOX(JP_LON,100), CLDBOX(JP_LON,100,JP_LEV)

!- dummy integer variables
!
!     DUMMY INTEGER SCALARS
INTEGER_M :: KAER, KIDIA, KFDIA, KTDIA, KLEV, KLON, KLW, KSW, KMODE

INTEGER_M :: ICLDMIX, IFLAG, IGEO, ILAY, KBOX, KULOUT, NDUMP, NFLUX, &
  & NRAD, NPRINT, NINDAT, NSSSSS,&
  & MONTH, IDAY, IYR, ILWRAD, IFIN, NMM, NDD, NAA, IMINUT, INHOMF, &
  & ILCLOUD, IOP, IOPTPROP, IO3ONLY, ICEWAT, NDAYS, NTOT, ICALC, &
  & IRADIP, IRADLP, NBOX, INOVLP, NWSET
  
INTEGER_M :: JL, JK, JNU, JAER, JA, JSW, ITIMES, ITIMSW, IDFIL, ICBOX, &
  & IS, ISTART, ISEND, IV, IFTIM, IRTIM, NINDATSF, ITROSTRA, ICORSFPR, &
  & JK0
  
INTEGER_M :: IWRT(90:96)
  
!- dummy logicals
LOGICAL :: LCLOUD, LFORCLD, LTIMTST, LMCCLAT
LOGICAL :: LEVOIGT, LLGEOSE, LLGEOSW, LLGMS, LLINDSA, LLMTO

!- dummy real variables

!REAL_B :: ZALB, ZALBD, ZALBP, ZEMIS, ZEMIW, ZLSM,&
!  & RTIMTR,&
!  & RPI, RTWAT, RDAY, RMD, RKBOL, RNAVO, R,&
!  & RCPD, RGAMMAS, RI0, RSIGMA
REAL_B :: ZALB, ZALBD, ZALBP, ZEMIS, ZEMIW, ZLSM,&
  & RTIMTR, RGAMMAS, ZLONP, ZLATP, ZTETA
REAL_B :: PRII05, PCCO25, ZEPAER,  ZTHETOZ, ZANGOZC, ZDEGRAD,&
  & ZTT5, ZCONDW, ZFRACT, ZZMU0, ZPP5, ZRR5, ZDDQQ, ZEE, ZQQ5
REAL_B :: ZCOOLC, ZCOOLT, ZHEATC, ZHEATT  

REAL_B :: RTETA,REL,REM,RRS,RLLS,RLLLS,RDS,RET
REAL_B :: PTIME,PTETA, ZPHOUR
REAL_B :: ZMFDLT,ZMFULT,ZMFDST,ZMFUST,ZMFDLC,ZMFDSC
REAL_B :: ZMFDSTD,ZMFUSTD,ZMFDSCD,ZMFDSI, ZSUDUR
REAL_B :: VAR2D(16)

!-- Interpolation
REAL_B :: ZDDP

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!- load external functions      
!
#include "fctast.h"
#include "fcttim.h"
#include "fcttre.h"
!     ------------------------------------------------------------------
!- define external parameters 
!  CLOUD cloud fraction      
!--- clear-sky
!DATA CLOUD / 60*0. /
DATA CLOUD / 46*0., 0.90, 0.90, 0.90, 11*0. /
!--- low-level overcast water clouds  (MLS between 730 and 870
DATA CLOUD1 / 46*0., 0.90, 0.90, 0.90, 11*0. /
DATA CLOUD2 / 46*0., 0.90, 0.90, 0.90, 11*0. /
DATA CLOUD3 / 46*0., 0.90, 0.90, 0.90, 11*0. /
DATA CLOUD4 / 46*0., 0.90, 0.90, 0.90, 11*0. /
DATA CLOUD5 / 46*0., 0.90, 0.90, 0.90, 11*0. /
!--- high-level overcast ice clouds (TRO 27-28-29    SAW 35-36-37)
!DATA CLOUD1 / 26*0., 0.90, 0.90, 0.90, 31*0. /
!DATA CLOUD2 / 28*0., 0.90, 0.90, 0.90, 29*0. /
!DATA CLOUD3 / 32*0., 0.90, 0.90, 0.90, 25*0. /
!DATA CLOUD4 / 30*0., 0.90, 0.90, 0.90, 27*0. /
!DATA CLOUD5 / 34*0., 0.90, 0.90, 0.90, 23*0. /
!--- test cloud more or less convective tower
!DATA  CLOUD / 42*0., 0.95, 0.90, 0.85, 0.80, 0.6, 0.4, 9*0.15, 7*0. /

!     ------------------------------------------------------------------
!
      DATA AAM / &
     &     0.000000,    20.000000,    38.425343, &
     &    63.647804,    95.636963,   134.483307, &
     &   180.584351,   234.779053,   298.495789, &
     &   373.971924,   464.618134,   575.651001, &
     &   713.218079,   883.660522,  1094.834717, &
     &  1356.474609,  1680.640259,  2082.273926, &
     &  2579.888672,  3196.421631,  3960.291504, &
     &  4906.708496,  6018.019531,  7306.631348, &
     &  8765.053711, 10376.126953, 12077.446289, &
     & 13775.325195, 15379.805664, 16819.474609, &
     & 18045.183594, 19027.695313, 19755.109375, &
     & 20222.205078, 20429.863281, 20384.480469, &
     & 20097.402344, 19584.330078, 18864.750000, &
     & 17961.357422, 16899.468750, 15706.447266, &
     & 14411.124023, 13043.218750, 11632.758789, &
     & 10209.500977,  8802.356445,  7438.803223, &
     &  6144.314941,  4941.778320,  3850.913330, &
     &  2887.696533,  2063.779785,  1385.912598, &
     &   855.361755,   467.333588,   210.393890, &
     &    65.889244,     7.367743,     0.000000, &
     &     0.000000 &
     &         /
!     
      DATA BBM / &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000000000, 0.0000000000, 0.0000000000, &
     & 0.0000758235, 0.0004613950, 0.0018151561, &
     & 0.0050811190, 0.0111429105, 0.0206778757, &
     & 0.0341211632, 0.0516904071, 0.0735338330, &
     & 0.0996746942, 0.1300225109, 0.1643843204, &
     & 0.2024759352, 0.2439331412, 0.2883229554, &
     & 0.3351548910, 0.3838921487, 0.4339629412, &
     & 0.4847715795, 0.5357099175, 0.5861684084, &
     & 0.6355474591, 0.6832686067, 0.7287858129, &
     & 0.7715966105, 0.8112534285, 0.8473749161, &
     & 0.8796569109, 0.9078838825, 0.9319403172, &
     & 0.9518215060, 0.9676452279, 0.9796627164, &
     & 0.9882701039, 0.9940194488, 0.9976301193, &
     & 1.0000000000 &
     &         / 
!
!     ------------------------------------------------------------------
!
REAL_B :: SLAT(42), SLON(42),HEIGHT(42),DPR(42), SHIFT(42)
INTEGER_M :: IDLAT(42), IDLON(42)
CHARACTER*15 CSTA(42)
CHARACTER*3  CIDS(42)

REAL_B :: SLATS(42), SLONS(42),HEIGHTS(42),DPRS(42), SHIFTS(42)
INTEGER_M :: IDLATS(42), IDLONS(42)
CHARACTER*15 CSTAS(42)
CHARACTER*3  CIDSS(42)

REAL_B :: SLATV(42), SLONV(42),HEIGHTV(42),DPRV(42), SHIFTV(42)
INTEGER_M :: IDLATV(42), IDLONV(42)
CHARACTER*15 CSTAV(42)
CHARACTER*3  CIDSV(42)

      DATA (CSTAS(IS),SLATS(IS),SLONS(IS),IDLATS(IS),IDLONS(IS),HEIGHTS(IS) &
& , CIDSS(IS),DPRS(IS) ,IS=1,42 ) / &
&   'NyAlesund'     , 78.560,   11.560,   22,   20,   137., 'NyA', 17. &
& , 'Barrow'        , 71.267, -156.833,   34,  361,     9., 'Bar',  0. &
& , 'Lerwick'       , 60.000,   -1.000,    0,    0,     0., 'Ler',  0. &
& , 'Toravere'      , 58.270,   26.470,    0,    0,     0., 'Tor',  0. &
& , 'Lindenberg'    , 52.210,   14.120,    0,    0,     0., 'Lin',  0. &
& , 'Chibolton'     , 51.400,   -1.400,    0,    0,    20., 'Chi',  0. & 
& , 'Regina'        , 50.120, -104.430,   71,  455,   594., 'Reg',  0. &
& , 'Camborne'      , 50.000,   -5.000,    0,    0,     0., 'Cam',  0. &
& , 'FortPeck'      , 48.310, -105.100,   74,  455,   634., 'FPk', 15. &
& , 'Xilinhot'      , 47.540,  109.590,    0,    0,     0., 'Xil',  0. &
& , 'Budapest'      , 47.500,   19.050,    0,    0,     0., 'Bud',  9. &
& , 'Payerne'       , 46.817,    6.950,   76,   12,   550., 'Pay',  9. &
& , 'Carpentras'    , 44.050,    5.030,   83,    9,   467., 'Car', 25. &
& , 'PennStateU'    , 40.720,  -77.930,    0,    0,   376., 'PSU', 21. &
& , 'TableMountain' , 40.125, -105.227,   88,  454,  1689., 'TMn', 92. &  
& , 'Bondville'     , 40.050,  -88.370,   89,  483,   213., 'Bon',  0. &
& , 'Boulder'       , 40.033, -105.267,   89,  454,  2293., 'Bou', 92. &
& , 'Chesapeake'    , 36.900,  -75.706,    0,    0,     0., 'Che',  0. &
& , 'SGP_Billings'  , 36.370,  -97.300,   96,  468,   332., 'SGP',-11. &
& , 'Tateno'        , 36.030,  140.080,   97,  251,    32., 'Tat', 23. &
& , 'DesertRockNV'  , 36.630, -116.020,    0,    0,  1007., 'DRk', 46. &
& , 'GoodwinCreek'  , 34.250,  -89.870,   99,  481,    98., 'GwC', -4. &
& , 'Bermuda'       , 32.300,  -64.750,  103,  525,    -4., 'Ber',  0. &
& , 'SedeBoger'     , 30.867,   34.766,    0,    0,     0., 'SdB',  0. &
& , 'Riyadh'        , 24.650,   48.767,    0,    0,     0., 'Riy',  0. &
& , 'Kwajalein'     ,  8.720,  167.731,  145,  299,   -10., 'Kwa',  0. &
& , 'Ilorin'        ,  8.533,    4.567,  146,    8,   383., 'Ilo',  0. & 
& , 'Maldives'      ,  5.000,   73.000,    0,    0,     0., 'Mal',  0. &
& , 'ARM_Nauru'     , -0.521,  166.916,  164,  262,     4., 'Nau',  0. &
& , 'ARM_Manus'     , -2.070,  147.430,  164,  262,     4., 'Man',  0. &
& , 'Balbina'       , -3.100,  -60.000,    0,    0,     0., 'Bal',  0. &
& , 'AltaFloresta'  , -9.917,  -56.017,    0,    0,     0., 'AFl',  0. &
& , 'AbracosHill'   ,-10.750,  -62.350,    0,    0,     0., 'AbH',  0. &
& , 'Ndola'         ,-12.983,   28.650,    0,    0,     0., 'Ndo',  0. &
& , 'Zambezi'       ,-13.517,   23.100,    0,    0,     0., 'Zam',  0. &
& , 'Mongu'         ,-15.250,   23.150,    0,    0,     0., 'Mon',  0. &
& , 'AliceSprings'  ,-27.280,  133.520,    0,    0,     0., 'AlS',-13. &
& , 'Florianopolis' ,-27.583,  -48.517,  209,  555,    31., 'Flo', 48. &
& , 'DeAar'         ,-30.400,   23.590,    0,    0,     0., 'DeA',  0. &
& , 'Syowa'         ,-69.000,  -39.350,  283,  571,     5., 'Syo',  0. &
& , 'GeorgVNeumayer',-70.390,   -8.150,  286,  626,    26., 'GvN', -4. &
& , 'SouthPole'     ,-89.980,  -24.480,  321,  597,  2878., 'SPo',-16. &  
& /
       

      data (CSTAV(IS),SLATV(IS),SLONV(IS),IDLATV(IS),IDLONV(IS),HEIGHTV(IS) &
     & ,SHIFTV(IS),CIDSV(IS),DPRV(IS) ,IS=1,22 ) / &
     &   'C1' ,  36.605, -97.485,  0,  0, 318.,  3240., 'C01' ,  0. &
     & , 'E1' ,  38.202, -99.316,  0,  0, 632.,  9060., 'E01' ,  0. &
     & , 'E2' ,  38.306, -97.301,  0,  0, 450.,  9060., 'E02' ,  0. &
     & , 'E3' ,  38.201, -95.597,  0,  0, 338.,  9120., 'E03' ,  0. &
     & , 'E4' ,  37.953, -98.329,  0,  0, 513.,  9180., 'E04' ,  0. &
     & , 'E5' ,  38.114, -97.513,  0,  0, 440.,  9120., 'E05' ,  0. &
     & , 'E6' ,  37.842, -97.020,  0,  0, 409.,  9180., 'E06' ,  0. &
     & , 'E7' ,  37.383, -96.180,  0,  0, 283., 27060., 'E07' ,  0. &
     & , 'E8' ,  37.333, -99.309,  0,  0, 664., 12660., 'E08' ,  0. &
     & , 'E9' ,  37.133, -97.266,  0,  0, 386., 27060., 'E09' ,  0. &
     & , 'E10',  37.068, -95.788,  0,  0, 248., 27060., 'E10' ,  0. &
     & , 'E11',  36.881, -98.285,  0,  0, 360., 26940., 'E11' ,  0. &
     & , 'E12',  36.841, -96.427,  0,  0, 331.,     0., 'E12' ,  0. &
     & , 'E13',  36.605, -97.485,  0,  0, 318.,  3360., 'E13' ,  0. &
     & , 'E15',  36.431, -98.284,  0,  0, 418., 26160., 'E15' ,  0. &
     & , 'E16',  36.061, -99.134,  0,  0, 602.,     0., 'E16' ,  0. &
     & , 'E18',  35.687, -95.856,  0,  0, 217., 26100., 'E18' ,  0. &
     & , 'E19',  35.549, -98.020,  0,  0, 421.,  8160., 'E19' ,  0. &
     & , 'E20',  35.564, -96.988,  0,  0, 309.,     0., 'E20' ,  0. &
     & , 'E22',  35.354, -98.977,  0,  0, 465., 26820., 'E22' ,  0. &
     & , 'E24',  34.883, -98.205,  0,  0, 409.,     0., 'E24' ,  0. &
     & , 'E25',  35.246, -96.736,  0,  0, 277., 26160., 'E25' ,  0. &
     & /
               
!     ------------------------------------------------------------------
!- ZALB SW surface albedo
!  ZEMIS LW emissivity outside the LW window region
!  ZEMIW LW emissivity within the 8-12.5 micron window region
!  ZLSM  land/sea mask   (1.= land   0.=ocean)
ZALB = 0.200
ZEMIS= 1.000
ZEMIW= 1.000
ZLSM = 1.0
!     ------------------------------------------------------------------
!!- basic constants 
           
RPI=3.1415926535898
RTT=273.16
REA=149597870000.
RTWAT=RTT
RTICE=RTT-23.
RDAY=86400.
RG=9.80665
RMD=28.9644
RKBOL=1.380658E-23
RNAVO=6.0221367E+23
R=RNAVO*RKBOL
RD=1000.*R/RMD
RCPD=3.5*RD
RI0=1370.
RSIGMA=5.67E-08
!
RGAMMAS=0.02
ICLDMIX=3
LCLOUD=.FALSE.
LFORCLD=.FALSE.
LTIMTST=.FALSE.
LHVOLCA=.FALSE.
LNEWAER=.FALSE.
IRADIP = 2
LRADLP =.FALSE.
LO3ONLY=.FALSE.
LRRTM  =.FALSE.
LINHOM =.FALSE.
LPHYLIN=.FALSE.
LMCCLAT=.FALSE.
ITROSTRA=-1
ICORSFPR=0
NHOWINH=0

ZMFDLT=0.
ZMFULT=0.
ZMFDST=0.
ZMFUST=0.
ZMFDLC=0.
ZMFDSC=0.
ZMFDSI=0.
ITIMES=0
!
ZMFDSTD=0.
ZMFUSTD=0.
ZMFDSCD=0.
ITIMSW=0

ZMFDLTC=0.
ITIMLWC=0

ZMFDSTC=0.
ITIMSWC=0
!
!     ------------------------------------------------------------------
!
!          1.  DEFINE THE CONFIGURATION
!
!          
KULOUT = 6
    
print *,' Hello we start ' 
!WRITE(0,*)'Enter NDUMP for prints 2=only end, 1=a bit more, 0=all'
!print *,  'Enter NDUMP for prints 2=only end, 1=a bit more, 0=all'
!read (*,9011) NDUMP
!print *,'ndump= ',ndump
NDUMP=2
9011 format(I2)     
if (NDUMP.LT.2) LDEBUG=.TRUE.

!WRITE(0,*)'      NWSET=-1 ICRCCM3    '
!WRITE(0,*)'Enter NWSET= 0 McClatchey '
!WRITE(0,*)'      NWSET= 1 ARM_BSRN_SURFRAD'
!WRITE(0,*)'      NWSET= 2 ARM_SGP_C1_En'
!read (*,9011) NWSET
NWSET=0
print *,'nwset= ',nwset

KLEV=60
if (NWSET.eq.-1) then
  NSTAT=1
  ISTAT=1
  LMCCLAT=.FALSE.
  KLON=1
  KIDIA=1
  KFDIA=1
  WRITE(0,*) 'Be sure than for 1.  ATEX        KLEV=83'
  WRITE(0,*) '             for 2.  BOMEX       KLEV=84'
  WRITE(0,*) '             for 3.  GATE_A      KLEV=46'
  WRITE(0,*) '             for 4.  GATE_B      KLEV=46'
  WRITE(0,*) '             for 5.  GATE_C      KLEV=46'
  WRITE(0,*) '             for 6.  OPEN_CELLS  KLEV=63'
  WRITE(0,*) 'Which ICRCCM3 atmosphere 1 to 6'
  READ (*,9011) ICRCCM3
  if (ICRCCM3.eq.1) then
    KLEV=83
    CPROF='ATEX'
    CIRCC='IC1'
  else if (ICRCCM3.eq.2) then
    KLEV=84
    CPROF='BOMEX'
    CIRCC='IC2'
  else if (ICRCCM3.eq.3) then
    KLEV=46
    CPROF='GATE_A'
    CIRCC='IC3'
  else if (ICRCCM3.eq.4) then
    KLEV=46
    CPROF='GATE_B'
    CIRCC='IC4'
  else if (ICRCCM3.eq.5) then
    KLEV=46
    CPROF='GATE_C'
    CIRCC='IC5'
  else if (ICRCCM3.eq.6) then
    KLEV=63
    CPROF='OPEN_CELLS'
    CIRCC='IC6'
  end if  
  print *,'------------------------'
  print *,CIRCC,CPROF,KLEV
  print *,'------------------------'
  
else if (NWSET.eq.0) then
  KLEV=60
  NSTAT=1
  ISTAT=0 
  LMCCLAT=.TRUE.
  KLON=1
  KIDIA=1
  KFDIA=1

else if (NWSET.eq.1) then
  KLEV=60
  NSTAT=42
  KLON=1
  KIDIA=1
  KFDIA=1
  DO IS=1,NSTAT
    CSTA(IS)=CSTAS(IS)
    SLAT(IS)=SLATS(IS)
    SLON(IS)=SLONS(IS)
    IDLAT(IS)=IDLATS(IS)
    IDLON(IS)=IDLONS(IS)
    HEIGHT(IS)=HEIGHTS(IS)
    CIDS(IS)=CIDSS(IS)
    DPR(IS)=DPRS(IS)
  END DO  
       
WRITE(0,*)'Station names and indices'
WRITE(0,*)  'NyAlesund      NyA =  1'
WRITE(0,*)  'Barrow         Bar =  2'
WRITE(0,*)  'Lerwick        Ler =  3'
WRITE(0,*)  'Toravere       Tor =  4'
WRITE(0,*)  'Lindenberg     Lin =  5'
WRITE(0,*)  'Chibolton      Chi =  6'
WRITE(0,*)  'Regina         Reg =  7'
WRITE(0,*)  'Camborne       Cam =  8'
WRITE(0,*)  'FortPeck       FPk =  9'
WRITE(0,*)  'Xilinhot       Xil = 10'
WRITE(0,*)  'Budapest       Bud = 11'
WRITE(0,*)  'Payerne        Pay = 12'
WRITE(0,*)  'Carpentras     Car = 13'
WRITE(0,*)  'PennStateU     PSU = 14'
WRITE(0,*)  'TableMountain  TMn = 15'  
WRITE(0,*)  'Bondville      Bon = 16'
WRITE(0,*)  'Boulder        Bou = 17'
WRITE(0,*)  'Chesapeake     Che = 18'
WRITE(0,*)  'SGP_Billings   SGP = 19'
WRITE(0,*)  'Tateno         Tat = 20'
WRITE(0,*)  'DesertRockNV   DRk = 21'
WRITE(0,*)  'GoodwinCreek   GwC = 22'
WRITE(0,*)  'Bermuda        Ber = 23'
WRITE(0,*)  'SedeBoger      SdB = 24'
WRITE(0,*)  'Riyadh         Riy = 25'
WRITE(0,*)  'Kwajalein      Kwa = 26'
WRITE(0,*)  'Ilorin         Ilo = 27' 
WRITE(0,*)  'Maldives       Mal = 28'
WRITE(0,*)  'ARM_Nauru      Nau = 29'
WRITE(0,*)  'ARM_Manus      Man = 30'
WRITE(0,*)  'Balbina        Bal = 31'
WRITE(0,*)  'AltaFloresta   AFl = 32'
WRITE(0,*)  'AbracosHill    AbH = 33'
WRITE(0,*)  'Ndola          Ndo = 34'
WRITE(0,*)  'Zambezi        Zam = 35'
WRITE(0,*)  'Mongu          Mon = 36' 
WRITE(0,*)  'AliceSprings   AlS = 37'
WRITE(0,*)  'Florianopolis  Flo = 38'
WRITE(0,*)  'DeAar          DeA = 39'
WRITE(0,*)  'Syowa          Syo = 40'
WRITE(0,*)  'GeorgVNeumayer GvN = 41'
WRITE(0,*)  'SouthPole      SPo = 42'  


WRITE(0,*)'Enter index of station; 99 is for all 42 stations in one go'
WRITE(0,*)'Enter index of station; 0 for McClatchey_on_60levels'
read(*,8900) ISTAT
  if (ISTAT.eq.99) then
    ISTART=1
    ISEND=42
    LMCCLAT=.FALSE.
  end if   

else if (NWSET.EQ.2) then
  KLEV=60
  NSTAT=22
  KLON=1
  KIDIA=1
  KFDIA=1
  DO IS=1,NSTAT
    CSTA(IS)=CSTAV(IS)
    SLAT(IS)=SLATV(IS)
    SLON(IS)=SLONV(IS)
    IDLAT(IS)=IDLATV(IS)
    IDLON(IS)=IDLONV(IS)
    HEIGHT(IS)=HEIGHTV(IS)
    SHIFT(IS)=SHIFTV(IS)
    CIDS(IS)=CIDSV(IS)
    DPR(IS)=DPRV(IS)
  END DO  
               
WRITE(0,*)'Station names and indices'
WRITE(0,*)  'C1 =  1'
WRITE(0,*)  'E1 =  2'
WRITE(0,*)  'E2 =  3'
WRITE(0,*)  'E3 =  4'
WRITE(0,*)  'E4 =  5'
WRITE(0,*)  'E5 =  6'
WRITE(0,*)  'E6 =  7'
WRITE(0,*)  'E7 =  8'
WRITE(0,*)  'E8 =  9'
WRITE(0,*)  'E9 = 10'
WRITE(0,*)  'E10= 11'
WRITE(0,*)  'E11= 12'
WRITE(0,*)  'E12= 13'
WRITE(0,*)  'E13= 14'
WRITE(0,*)  'E15= 15'  
WRITE(0,*)  'E16= 16'
WRITE(0,*)  'E18= 17'
WRITE(0,*)  'E19= 18'
WRITE(0,*)  'E20= 19'
WRITE(0,*)  'E22= 20'
WRITE(0,*)  'E24= 21'
WRITE(0,*)  'E25= 22'

WRITE(0,*)'Enter index of station; 99 is for all 22 stations in one go'
WRITE(0,*)'Enter index of station; 0 for McClatchey_on_60levels'
read(*,8900) ISTAT
  if (ISTAT.eq.99) then
    ISTART=1
    ISEND=22
    LMCCLAT=.FALSE.
  end if   

else
  WRITE(0,*) 'Station ID is not consistent: STOP!!!!'
  STOP
end if  



8900 format(I2)

if (ISTAT.EQ.0.and.NWSET.eq.0) then
   ISTART=1
   ISEND=5
   LMCCLAT=.TRUE.
   WRITE(0,*) 'Enter Emphasis on tropo- (0,2) or strato- sphere (1,3)?'
   WRITE(0,*) '0,1: clear-sky and total LW; 2,4: clear-sky and total SW'
   read (*,8900) ITROSTRA
   LCLOUD=.TRUE.
   LFORCLD=.TRUE.
else if (ISTAT.ne.0 .and. ISTAT.ne.99) then
   ISTART=ISTAT
   ISEND=ISTAT
   LMCCLAT=.FALSE.
end if      

IWRT(90)=0
IWRT(91)=0
IWRT(92)=0
IWRT(93)=0
IWRT(94)=0
IWRT(95)=0
IWRT(96)=0
!WRITE(0,*)'What to store 0=No, 1=Yes'
!WRITE(0,*) '90_Input_Atmospheres'
!read(*,8902) IWRT(90)
!WRITE(0,*) '91_Output_radiation_fluxes'
!read(*,8902) IWRT(91)
!WRITE(0,*) '92_Output_radiative_heating'
!read(*,8902) IWRT(92)
!WRITE(0,*) '93_Output_SurfDownFluxes'
!read(*,8902) IWRT(93)
!WRITE(0,*) '94_Output_BoxToASurfFluxesNTau'
!read(*,8902) IWRT(94)
!WRITE(0,*) '95_Output_BoxCloud'
!read(*,8902) IWRT(95)
!WRITE(0,*) '96_Output_SurfDownSW'
!read(*,8902) IWRT(96)

IWRT(90)=1
IWRT(91)=1
IWRT(92)=1
IWRT(93)=1
IWRT(94)=0
IWRT(95)=0
IWRT(96)=0

WRITE(0,*)'Land=1 or ocean=0'
print *  ,'Land=1 or ocean=0'
read(*,8902) ILSM
8902 format(I1)
ZLSM=FLOAT(ILSM)


!WRITE(0,*)'Start DATE YYYYMMDD; irrelevant for McClatchey atmospheres'
!read(*,8903) NINDAT
8903 format(I8)
!WRITE(0,*)'Start time sec I5'
!read(*,8904) NSSSSS
8904 format(I5)
NINDAT=19990401
NSSSSS=43200


!ILAY=60    
!WRITE(0,*)'Enter NLAY number of levels'
!WRITE(0,*)' McClatchey, ARM_SGP, ARM_BSRN_SURFRAD are 60 levels'
!WRITE(0,*)' ATEX=64 ; BOMEX=83 ; GATE_A,B,C=46 ; OPEN_CELLS=63'
!read (*,9011) ILAY
ILAY=KLEV
      
!KIDIA=JP_IDIA
!KFDIA=JP_FDIA
!KLON =JP_LON
KTDIA=JP_TDIA
!KLEV =JP_LEV
KMODE=JP_MODE
KAER =JP_AER
KLW  =JP_LW
KBOX =0
NBOX =1
    
WRITE(0,*)'Enter KSW number of SW spectral intervals: 2, 4 or 6'
read (*,9011) KSW
!KSW=4
PRINT *,'KSW number of SW spectral intervals: ',KSW

WRITE(0,*)'Enter NAER aerosols  no=0 Tanre=1 GADS=2'
read (*,9011) NAER
if (NAER.EQ.2) THEN
  LNEWAER=.TRUE.
END IF
print *,'LNEWAER= ',LNEWAER

!WRITE(0,*)'Enter Box or not? KBOX=1=box, KBOX=0, no box-type computat.'
!read (*,9011) KBOX
KBOX=0

if (KBOX.EQ.1) then
  WRITE(0,*) 'How many boxes? 10 to 100 is fine: NB Not over 100'
  read (*,8798) NBOX
  if (NBOX.LT.10) then
    write(unit=CBOX,fmt=8796) NBOX
8796 format('00',I1)        
   else if (NBOX.GE.10 .AND. NBOX.LT.100) then
    write(unit=CBOX,fmt=8797) NBOX
8797 format('0',I2)        
   else if (NBOX.GE.100) then
    write(unit=CBOX,fmt=8798) NBOX
8798 format(I3)
   end if        
else
   NBOX=1
   CBOX='000'
end if    

     
!- call to set-up routines for various coefficients
!  N.B. These are independent from the vertical resolution

!- Load Ozone climatology      
    
!WRITE(0,*) 'Enter NOZOCL'
!WRITE(0,*) 'if NOZOCL =-1 tests the vertical quadrature '
!WRITE(0,*) '                                  (no absorber)'
!WRITE(0,*) 'if NOZOCL =0 whatever is read for O3 as input profile'
!WRITE(0,*) '   NOZOCL =1 old ECMWF O3 and climatol. aerosols'
!WRITE(0,*) '   NOZOCL =2 Fortuin-Langematz O3 climatology + aero'
!WRITE(0,*) '   NOZOCL =3 old ECMWF O3 and no aerosol '
!WRITE(0,*) '   NOZOCL =4 Fortuin-Langematz O3 and no aerosol '
!read (*,9012) NOZOCL
NOZOCL=2
WRITE(0,*) '   NOZOCL =2 Fortuin-Langematz O3 climatology + aero'
9012 format(I2)   

!WRITE(0,*) 'Enter IO3ONLY  =0 Oper.   1=Only O3 absorption in UV-Vis.'
!read (*,9012) IO3ONLY
IO3ONLY=0
!if (IO3ONLY.EQ.1) then
!  LO3ONLY=.TRUE.
!else
  LO3ONLY=.FALSE.  
!end if  
 
 
 
 
  
!WRITE(0,*) 'Enter ILWRAD =0 Morcrette, 1991 operational before 20000627'
!WRITE(0,*) '             =1 Mlawer et al., 1997 now ECMWF-operational'
!WRITE(0,*) 'Enter ILWRAD =2 Morcrette, 1991 original as in ERA-15'
!read (*,9012) ILWRAD
ILWRAD=1
if (ILWRAD.EQ.0) then
  LRRTM=.FALSE.
  WRITE(0,*) 'Enter NOVLP =1 MRN, 2=MAX, 3=RAN, 4=Hogan'
  read (*,9012) NOVLP
  
else if (ILWRAD.EQ.1) then
  LRRTM=.TRUE.
  NOVLP=1
  
else if (ILWRAD.EQ.2) then
  LRRTM=.FALSE.
  NOVLP=1
  NRADIP=0
  NRADLP=0
  NICEOPT=0
  NLIQOPT=0
end if  
print *,'NOVLP= ',NOVLP


!WRITE(0,*) '0=clear-sky  1=Cloudy Computations'
!read (*,9012) ILCLOUD
ILCLOUD=1
if (ILCLOUD.eq.1) then
  LCLOUD=.TRUE.

  if (ILWRAD.eq.0) then
    WRITE(0,*) 'Optical Properties M91/G98 IOPTPROP=0'
    WRITE(0,*) 'SW Water: Fouquart, 1987; SW Ice: Ebert-Curry, 1992'
    WRITE(0,*) 'for M91/G98 LW Water        : Smith-Shi, 1992'
    WRITE(0,*) '            LW Ice          : Ebert-Curry, 1992'
    IOPTPROP=0
    LOWASYF=.FALSE.
    LOIFUEC=.FALSE.
    LOWHSSS=.FALSE.
    NICEOPT=1
    NLIQOPT=0
    NRADIP=2
    NRADLP=1
  else if (ILWRAD.eq.2) then
    WRITE(0,*) 'Optical Properties for M91 IOPTPROP=0 as in ERA-15'
    WRITE(0,*) 'SW Water: Fouquart, 1987; SW Ice: Ebert-Curry, 1992'
    WRITE(0,*) 'for M91     LW Water AND Ice: Smith-Shi, 1992'
    IOPTPROP=0
    LOWASYF=.FALSE.
    LOIFUEC=.FALSE.
    LOWHSSS=.FALSE.
    NICEOPT=0
    NLIQOPT=0
    NRADIP=0
    NRADLP=0
  else if (ILWRAD.eq.1) then  
    WRITE(0,*) 'for RRTM enter NICEOPT and NLIQOPT'
    WRITE(0,*) 'Enter NLIQOPT =0  Water LW: Smith-Shi, 1992; SW: Fouquart, 1987'
    WRITE(0,*) '              =1  Water LW: Savijarvi, 1997; SW: Slingo  , 1989'
    WRITE(0,*) '              =2  Water LW: Lindner,Li,2000; SW: Slingo  , 1989'
    read (*,9012) NLIQOPT
    WRITE(0,*) 'Enter NICEOPT =0  Ice LW: Smith,Shi  , 1992; SW: Ebert-Curry, 1992'
    WRITE(0,*) '              =1  Ice LW: Ebert,Curry, 1992; SW: Ebert-Curry, 1992'
    WRITE(0,*) '              =2  Ice LW: Fu,Liou    , 1993; SW: Fu,Liou    , 1993'
    WRITE(0,*) '              =3  Ice LW: Fu et al.  , 1998; SW: Fu         , 1996'
    read (*,9012) NICEOPT
  end if  
  
  WRITE(0,*) 'Enter ICEWAT  =0 Liquid Water    1=Ice    2=both '
  read (*,9012) ICEWAT
  if (LMCCLAT) then
    ICLDMIX=ICEWAT
  end if        
  

  WRITE(0,*) 'Enter IRADLP =0 effective radius - liquid as f(Pressure) '
  WRITE(0,*) '             =1 fixed 10 microns over land, 13 over ocean'
  WRITE(0,*) '             =2 computed from LWC Martin et al, 1994'
  read (*,9012) IRADLP
  NRADLP=IRADLP
  if (IRADLP.le.1) then
    LRADLP=.FALSE.
  else
    LRADLP=.TRUE.
  end if    

  WRITE(0,*) 'Enter IRADIP =0 fixed effective radius - ice 40 microns'
  WRITE(0,*) '             =1   f(T)   40 - 130 microns'
  WRITE(0,*) '             =2   f(T)   30 -  60 microns Jakob-Klein'
  WRITE(0,*) '             =3   f(T,IWC)  Sun-Rikus, 99'
  read (*,9012) IRADIP
  NRADIP=IRADIP

  CRMINICE='000'   
  IF (IRADIP.EQ.3) THEN
  
    WRITE(0,*) 'Enter minimum diameter for ice particles: f4.0 '
    read(*,9013) RMINICE
9013 format(f4.0)
    MINICE=INT(RMINICE)
    IF (MINICE.LT.10) THEN
      WRITE(UNIT=CRMINICE,FMT=9014) MINICE
9014  format('00',I1)
    ELSE IF (MINICE.LT.100) THEN
      WRITE(UNIT=CRMINICE,FMT=9015) MINICE
9015  format('0',I2)
    ELSE IF (MINICE.GE.100) THEN
      WRITE(UNIT=CRMINICE,FMT=9016) MINICE
9016  format(I3)
    END IF
    
  END IF    

!  WRITE(0,*) 'Enter INHOMF  =0 cloud tau is taken as is'
!  WRITE(0,*) 'Enter INHOMF  =1 some provision is made for inhomogeneity'
!  read (*,9012) INHOMF
  INHOMF=0
  if (INHOMF.EQ.0) then
    LINHOM=.FALSE.
  else if (INHOMF.EQ.1) then
    LINHOM=.TRUE.
    WRITE(0,*) 'Enter NHOWINH =1 MT 0.7 factor'
    WRITE(0,*) 'Enter NHOWINH =2 exp(-(sig/tau)^2)'
    WRITE(0,*) 'Enter NHOWINH =3 Cairns, 2000'
    read (*,9012) NHOWINH
  end if      

else
  LCLOUD=.FALSE.
end if    

!WRITE(0,*) 'Enter ICORSFPR =0 nothing done'
!WRITE(0,*) 'Enter ICORSFPR =1 corrected for Delta P'
!WRITE(0,*) 'Enter ICORSFPR =2 corrected for Delta P, ajusted T'
!WRITE(0,*) 'Enter ICORSFPR =3 corrected for Delta P, ajusted T and q'
!read(*,9012) ICORSFPR
ICORSFPR=0

INOVLP=3*NOVLP+ICORSFPR 
 
IOP=NICEOPT*1000000+NLIQOPT*100000+ICEWAT*10000+INHOMF*1000+IRADLP*100+IRADIP*10+INOVLP
if (IOP.lt.10) then
  write (unit=COP,fmt=8696) IOP
8696 format('000000',I1)  
else if (IOP.ge.10.and.IOP.lt.100) then
  write (unit=COP,fmt=8697) IOP
8697 format('00000',I2)  
else if (IOP.ge.100.and.IOP.lt.1000) then
  write (unit=COP,fmt=8698) IOP
8698 format('0000',I3)  
else if (IOP.ge.1000.and.IOP.lt.10000) then
  write (unit=COP,fmt=8699) IOP
8699 format('000',I4)  
else if (IOP.ge.10000.and.IOP.lt.100000) then
  write (unit=COP,fmt=8700) IOP
8700 format('00',I5)  
else if (IOP.ge.100000.and.IOP.lt.1000000) then
  write (unit=COP,fmt=8701) IOP
8701 format('0',I6)  
else if (IOP.ge.1000000) then
  write (unit=COP,fmt=8702) IOP
8702 format(I7)  
end if



IF (NWSET.GE.1) THEN      
  WRITE(0,*) 'How many days? (each with 24 one-hour time-steps)'
  WRITE(0,*) ' will in fact process as many profiles as provided in input file'
  read (*,8910) NDAYS
8910 format(I2)
  NTOT=NDAYS*24
ELSE
  NTOT=1
END IF

      
print *,' KIDIA, KFDIA, KLON ',KIDIA, KFDIA, KLON
print *,' KLEV, KMODE, NAER  ',KLEV, KMODE, NAER
print *,' KLW, KSW           ',KLW, KSW

!--      
LERAD6H=.TRUE.
LERADHS=.TRUE.
LONEWSW=.TRUE.
LEVOIGT=.FALSE.
!      
NFLUX  =6
NMODE  =0
NRAD   =1
NRADFR =-3
NRADPFR=36
NRADPLA=15
NRINT  =4
NRADF2C=1
NRADC2F=1
NUAER  = 31
NTRAER = 19
NLW    = 6
NTSW   = 6
NSW    = KSW
NPRINT = 1

!--

IGEO=0
LLGEOSE =.FALSE.
LLGEOSW =.FALSE.
LLGMS   =.FALSE.
LLINDSA =.FALSE.
LLMTO   =.FALSE.

!--

REPSC  = 1.E-12
REPSCO = 1.E-12
REPSCQ = 1.E-12
REPSCT = 1.E-12
REPSCW = 1.E-12
REPLOG = 1.E-12
NOUT   = 6

!--


     
!- call to set-up routines for various coefficients
!  N.B. These are independent from the vertical resolution

CALL SURDI
CALL SULWN
CALL SUOLW

CALL SURRTAB
CALL SURRTPK 
CALL SURRTRF
CALL SURRTFTR

CALL RRTM_KGB1
CALL RRTM_KGB2
CALL RRTM_KGB3
CALL RRTM_KGB4
CALL RRTM_KGB5
CALL RRTM_KGB6
CALL RRTM_KGB7
CALL RRTM_KGB8
CALL RRTM_KGB9
CALL RRTM_KGB10
CALL RRTM_KGB11
CALL RRTM_KGB12
CALL RRTM_KGB13
CALL RRTM_KGB14
CALL RRTM_KGB15
CALL RRTM_KGB16

! *** mji ***
! Initialization routine for RRTM
! Reduce absorption coefficient data from 256 to 140 g-points

CALL RRTM_INIT_140GP

WRITE(0,*) ' RRTM_INIT_140GP IS INITIALIZED '

write(unit=CONFIGUR,FMT=8913) ILWRAD,KSW,NAER,CBOX,COP
8913 format('Configuration_LW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7)
open (97,file=CONFIGUR)
!
!     ------------------------------------------------------------------
!
!          2.  DEFINE THE ATMOSPHERIC CASE

if (LMCCLAT) then
  ISTART=1
  ISEND=5
  CID(1)='TRO'
  SLAT(1)=0.
  SLON(1)=0.
  CID(2)='MLS'
  SLAT(2)=40.
  SLON(2)=0.
  CID(3)='MLW'
  SLAT(3)=40.
  SLON(3)=0.
  CID(4)='SAS'
  SLAT(4)=75.
  SLON(4)=0.
  CID(5)='SAW'
  SLAT(5)=75.
  SLON(5)=0.
  LFORCLD=.true.
else if (NWSET.GE.1) THEN
  do IS=ISTART,ISEND
    CID(IS)=CIDS(IS)
  end do  
else if (NWSET.eq.-1) then
  ISTART=1
  IEND=1
  CID(1)=CIRCC
  SLAT(1)=0.
  SLON(1)=0.
  KIDIA=1
  KFDIA=1
end if

DO ISTAT = ISTART , ISEND
  

ZLATP=SLAT(ISTAT)
ZLONP=SLON(ISTAT)
CSTAT=CID(ISTAT)
print *,'CID(',ISTAT,')=',CID(ISTAT),' Lat:',SLAT(ISTAT),' Lon:', &
 & SLON(ISTAT),' Height:',HEIGHT(ISTAT)
IF (ZLONP.LT.0.) THEN
  ZLONP=360.-ZLONP
END IF   

if (NWSET.GE.0) then   
  write (unit=CPROF,fmt=8901) CID(ISTAT)
8901 format(A3,'_TQ_CFIL_SurfToA_19990331_19990531')
else if (NWSET.EQ.-1) then
  print *,'ICRCCM3 input file ',CPROF,' with ',KLEV,' levels'
end if  


INPROF=11
INSURF=12
print *,'CPROF=',CPROF
open(INPROF,file=CPROF)


!- basic constants      
CALL SUCST ( KULOUT, NINDAT, NSSSSS, 1 )
if (LDEBUG) print *,RTT,RG,RD
     
!- call to set-up routines for various coefficients
!  N.B. These are independent from the vertical resolution

!- Load Ozone climatology      
       
8911 format(A132)
8912 format(15X,2F7.2)


DO JL=KIDIA,KFDIA
 PLAT5(JL)= ZLATP*RPI/180.
 PLON5(JL)= ZLONP*RPI/180.
END DO 


write(unit=CFILOUT0,FMT=8914) CSTAT,ILWRAD,KSW,NAER,CBOX,COP
write(unit=CFILOUT1,FMT=8915) CSTAT,ILWRAD,KSW,NAER,CBOX,COP
write(unit=CFILOUT2,FMT=8916) CSTAT,ILWRAD,KSW,NAER,CBOX,COP
write(unit=CFILOUT3,FMT=8917) CSTAT,ILWRAD,KSW,NAER,CBOX,COP,CRMINICE
write(unit=CFILOUT4,FMT=8918) CSTAT,ILWRAD,KSW,NAER,CBOX,COP
write(unit=CFILOUT5,FMT=8919) CSTAT,ILWRAD,KSW,NAER,CBOX,COP
write(unit=CFILOUT6,FMT=8920) CSTAT,ILWRAD,KSW,NAER,CBOX,COP
8914 format(A3,'_Input_Atmospheres_LW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7)
8915 format(A3,'_Output_radiation_fluxes_LW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7)
8916 format(A3,'_Output_radiative_heating_LW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7)
8917 format(A3,'_Output_SurfDownFluxes_LW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7,'_',A3)
8918 format(A3,'_Output_BoxToASurfFluxesNTau_LW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7)
8919 format(A3,'_Output_BoxCloud_LW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7)
8920 format(A3,'_Output_SurfDownSW',I1,'_SW',I1,'_AER',I1,'_BX',A3,'_OP',A7)


     
if (IWRT(90).eq.1) open (90,file=CFILOUT0)
if (IWRT(91).eq.1) open (91,file=CFILOUT1)
if (IWRT(92).eq.1) open (92,file=CFILOUT2)
if (IWRT(93).eq.1) open (93,file=CFILOUT3)
if (IWRT(94).eq.1) open (94,file=CFILOUT4)
if (IWRT(95).eq.1) open (95,file=CFILOUT5)
if (IWRT(96).eq.1) open (96,file=CFILOUT6)

DO IDFIL=90,96
  if (IWRT(IDFIL).EQ.1) then
    write (97,8921) ZLATP,ZLONP,PLAT5(KIDIA),PLON5(KIDIA)
8921 format(1x,2F10.3,2f12.6)
    write (97,8922) CPROF
!    write (97,8922) CSURF
8922 format(A80)
    write (97,8923) NOZOCL,ILCLOUD,IOPTPROP,INHOMF,ILWRAD,NAER,KSW,KBOX,NBOX,IRADLP
8923 format('OzoClim:',I2,' Cloud:',I2,' ClOptProp:',I2,' Inhom:',I2 &
    &,' LWscheme:',I2,' NAER:',I2,' KSW:',I2,' KBOX:',I2,' NBOX=',I3 &
    &,' RADLP:',I2)
  end if
END DO  



ICALC=0
ZPHOUR=-0.5

1000 CONTINUE
 ICALC=ICALC+1
 ZPHOUR=ZPHOUR+1.
 
!print *,'Input data' 
!-------------------------- ICRCCM-3 DATA SETS ---------------------
IF (NWSET.EQ.-1) THEN
  do JJJ=1,4
    read(INPROF,5000,END=9999) CBLABLA
    print 5000, CBLABLA
5000 format(A132) 
  enddo


  DO JK=1,KLEV
    read(INPROF,5001) ZAPH(JK),ZATH(JK),ZAOH(JK)
5001 format(28x,f9.3,8x,f8.3,4x,e14.7)

    read(INPROF,5002) ZACFL(JK),ZACFI(JK),ZACF(JK) &
  &              ,ZACQL(JK),ZACQI(JK),ZACQ(JK) &
  &              ,ZSCQL(JK),ZSCQI(JK),ZSCQ(JK) &
  &              ,ZAQ(JK)
!-- 80X up to layer#
!  +16X       mid-layer(km)
!  +16X       density 
!  +16X       cloud fraction (liquid)    <------ ZACFL
!  +16X       cloud fraction (ice)       <------ ZACFI
!  +16x       cloud fraction (both)      <------ ZACF
!  +16X       accum cld.frac. down (liq)
!  +16X       accum cld.frac. down (ice)
!  +16X       accum cld.frac. down (both)
!  +16X       accum cld.frac. up (liq)
!  +16X       accum cld.frac. up (ice)
!  +16X       accum cld.frac. up (both)
!  +16X       mixing ratio cloud liquid  <------ ZACQL
!  +16X       mixing ratio cloud ice     <------ ZACQI
!  +16X       mixing ratio cloud (both)  <------ ZACQ
!  +16X       std mr liq                 <------ ZSCQL
!  +16X       std mr ice                 <------ ZSCQI
!  +16X       std mr both                <------ ZSCQ
!  +16X       log(mr liq) 
!  +16X       log(mr ice) 
!  +16X       log(mr both) 
!  +16X       mixing ratio clear         <------ 
!  +16X       mixing ratio within-cloud  <------
!  +16X       mixing ratio total         <------ ZAQ
 
5002 format(80x,2(16x),3(f16.3),6(16x),3(f16.6),3(f16.6),5(16X),f16.6)

    print 5003, JK,ZAPH(JK),ZATH(JK),ZAOH(JK),ZAQ(JK) &
  &                   ,ZACFL(JK),ZACFI(JK),ZACF(JK) &
  &                   ,ZACQL(JK),ZACQI(JK),ZACQ(JK) &
  &                   ,ZSCQL(JK),ZSCQI(JK),ZSCQ(JK)
5003 format(1x,I3,2f10.3,e14.7,10f16.6)  
    ZAPH(JK)=ZAPH(JK)*100.  
  END DO
  
  JK=KLEV+1
  read(INPROF,5001) ZAPH(JK),ZATH(JK),ZAOH(JK)
  print 5003, JK,ZAPH(JK),ZATH(JK),ZAOH(JK)
  ZAPH(JK)=ZAPH(JK)*100.  
  
  DO JK=1,KLEV
    ZAQ(JK)=ZAQ(JK)*1.E-3    
    ZAO(JK)=0.5*(ZAOH(JK)+ZAOH(JK+1))
    ZAT(JK)=0.5*(ZATH(JK)+ZATH(JK+1))
    ZAP(JK)=0.5*(ZAPH(JK)+ZAPH(JK+1))
  END DO  
  DO JK=1,KLEV+1
    DO JL=KIDIA,KFDIA
      PAPH5(JL,JK)=ZAPH(JK)
      PTH5(JL,JK)=ZATH(JK)
    END DO
  END DO

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PT5(JL,JK)=ZAT(JK)
      PQ5(JL,JK)=ZAQ(JK)
      POZON5(JL,JK)=ZAO(JK)
      PAP5(JL,JK)=ZAP(JK)
      PCLFR5(JL,JK)=0.
      PSQLW5(JL,JK)=1.
      PSQIW5(JL,JK)=1.
      PQLWP5(JL,JK)=0.
      PQIWP5(JL,JK)=0.
      PRLVRI5(JL,JK)=0.
      PRLVRL5(JL,JK)=0.
      PRH5(JL,JK)=0.
    
      IF (ICEWAT.EQ.0) THEN
        PQLWP5(JL,JK)=ZACQL(JK)*1.E-03
        PQIWP5(JL,JK)=0.
        PCLFR5(JL,JK)=ZACFL(JK)
      ELSE IF (ICEWAT.EQ.1) THEN
        PQIWP5(JL,JK)=ZACQI(JK)*1.E-03
        PQLWP5(JL,JK)=0.
        PCLFR5(JL,JK)=ZACFI(JK)
      ELSE IF (ICEWAT.EQ.2) THEN
        PQLWP5(JL,JK)=ZACQL(JK)*1.E-03
        PQIWP5(JL,JK)=ZACQI(JK)*1.E-03
        PCLFR5(JL,JK)=ZACF (JK)
      ELSE IF (ICEWAT.EQ.-1) THEN
        PQLWP5(JL,JK)=ZACQ (JK)*1.E-03
        PQIWP5(JL,JK)=0.
        PCLFR5(JL,JK)=ZACF (JK)
      END IF
        
      IF (PCLFR5(JL,JK) .GT. 0. .AND. ICEWAT.NE.-1 ) THEN    
        IF (PQLWP5(JL,JK) .GT. 0.) THEN      
          PRLVRL5(JL,JK)= (ZSCQL(JK)/ZACQL(JK))**2
          PSQLW5(JL,JK)=EXP(- (ZSCQL(JK)/ZACQL(JK))**2 )
        END IF
        IF (PQIWP5(JL,JK) .GT. 0.) THEN  
          PRLVRI5(JL,JK)= (ZSCQI(JK)/ZACQI(JK))**2 
          PSQIW5(JL,JK)=EXP(- (ZSCQI(JK)/ZACQI(JK))**2 )
        END IF  
      END IF  
    
      IF (PCLFR5(JL,JK) .GT. 0. .AND. ICEWAT.EQ.-1 ) THEN    
        IF (PQLWP5(JL,JK) .GT. 0.) THEN      
          PRLVRL5(JL,JK)= (ZSCQ (JK)/ZACQ (JK))**2 
          PSQLW5(JL,JK)=EXP(- (ZSCQ (JK)/ZACQ (JK))**2 )
        END IF
      END IF  
    
    END DO
  END DO
  DO JL=KIDIA,KFDIA
    PTS5(JL)=ZATH(KLEV+1)
  END DO  
  
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
!---- old-style saturation function for compatibility with previous code  
      ZPP5 = PAP5(JL,JK)/100._JPRB
      ZTT5 = PT5(JL,JK)
      ZRR5 = 1._JPRB
      CALL QSAT  ( ZDDQQ, ZQQ5, ZEE, ZPP5, ZTT5, ZRR5 )
      PQS5(JL,JK) = ZQQ5
      PRH5(JL,JK) = PQ5(JL,JK)/PQS5(JL,JK)
    END DO
  END DO  
    
!- basic constants      
  CALL SUCST ( KULOUT, NINDAT, NSSSSS, 0 )

  ZALBD=0.2
  ZALBP=0.2
  ZTIM=0.
  print *,ZALBD,ZALBP,NINDAT,NSSSSS,ZTIM



!------------------------ ECMWF-TYPE DATA SETS ---------------------------------

ELSE IF (NWSET.GE.0) THEN

  DO JL=KIDIA,KFDIA

    if (NWSET.eq.0 .or. NWSET.eq.2) then
      read (INPROF,8952,END=9999) (VAR2D(IV),IV= 1, 3),NINDATSF,IFTIM
      read (INPROF,8952         ) (VAR2D(IV),IV= 4, 6)
8952  format(14X,3e13.6,1X,I8,I3) 
    else if (NWSET.eq.1) then
      read (INPROF,8951,END=9999) (VAR2D(IV),IV= 1, 3),NINDATSF,IFTIM
      read (INPROF,8951         ) (VAR2D(IV),IV= 4, 6)
8951 format(14X,3e13.6,5(13X),1X,I8,I3) 
    end if

!--- attention: initial NINDAT and IRTIM correspond to time reference in FC
! must be changed to real time.
!
    IRTIM=IFTIM
    if (IFTIM.GT.12) then
      ZTIM=FLOAT(IFTIM)-12.5
      IF (NINDATSF.EQ.19990331) then
        NINDAT=19990401
      ELSE IF (NINDATSF.EQ.19990430) then
        NINDAT=19990501
      ELSE IF (NINDATSF.EQ.19990531) then
        NINDAT=19990601
      ELSE    
        NINDAT=NINDATSF+1
      END IF  
    else
      ZTIM=FLOAT(IFTIM)+11.5
      NINDAT=NINDATSF
    end if    
    NSSSSS=INT(ZTIM*3600)

!- basic constants      
    CALL SUCST ( KULOUT, NINDAT, NSSSSS, 0 )
    if (LDEBUG) print *,RTT,RG,RD

    PAPH0(JL,KLEV+1)=VAR2D(1)
    PTS5(JL)=VAR2D(2)
    ZALBD=VAR2D(3)
    ZALBP=VAR2D(3)
  END DO
 
!print *,'Just as profiles are read '      
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      read (INPROF,8953         ) &
    & PAP0(JL,JK),PT0(JL,JK),PQ0(JL,JK),PCLFR5(JL,JK) &
    &,PQIWP5(JL,JK),PQLWP5(JL,JK),POZON5(JL,JK)
   8953 format(14X,8(1x,E13.6))
      PRH5(JL,JK)=0.

!_-- for Anton_IMET type file    
!    read (INPROF,8953,END=9999) &
!    & NINDAT,ZTIM,PAP5(JL,JK),PT5(JL,JK),PQ5(JL,JK) &
!    &,PQLWP5(JL,JK),PQIWP5(JL,JK),PCLFR5(JL,JK),PRH5(JL,JK)
!    8953 format(15X,I8,1x,f5.2,4X,(1x,E12.5),13x,13x,6(1x,E12.5),13x)
!_-- for ARM-SGP type file    
!!  8953 format(1X,I8,1x,f5.2,1X,7(1x,E13.6))
      IF (NDUMP.LE.3) then
        print 8953, &
        & PAP0(JL,JK),PT0(JL,JK),PQ0(JL,JK) &
        &,PQLWP5(JL,JK),PQIWP5(JL,JK),PCLFR5(JL,JK),PRH5(JL,JK),POZON5(JL,JK)
      END IF  
    
!      POZON5(JL,JK)=0.
    
    END DO
  END DO  


  DO jnu=1,KSW
    DO jl=KIDIA,KFDIA
      PALBD5(JL,JNU)=ZALBD
      PALBP5(JL,JNU)=ZALBP
    ENDDO 
  ENDDO 

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PAPH0(JL,JK  )=AAM(JK  )+BBM(JK  )*PAPH0(JL,KLEV+1)
      PAPH0(JL,JK+1)=AAM(JK+1)+BBM(JK+1)*PAPH0(JL,KLEV+1)
      PSQLW5(JL,JK)=1.
      PSQIW5(JL,JK)=1.
      PRLVRI5(JL,JK)=0.
      PRLVRL5(JL,JK)=0.
    END DO
  END DO 
  
         
!-- tentative correction for model-observation difference in orography

!-- compute the model relative humidity
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
!---- old-style saturation function for compatibility with previous code  
      ZPP5 = PAP0(JL,JK)/100._JPRB
      ZTT5 = PT0(JL,JK)
      ZRR5 = 1._JPRB
      CALL QSAT  ( ZDDQQ, ZQQ5, ZEE, ZPP5, ZTT5, ZRR5 )
      PQS0(JL,JK) = ZQQ5
      PRH0(JL,JK) = PQ0(JL,JK)/PQS0(JL,JK)
    END DO
  END DO    

!JL=KIDIA
!print 8771,PAPH0(JL,KLEV+1),(PQ0(JL,JK),PT0(JL,JK),PQS0(JL,JK),PRH0(JL,JK),JK=KLEV,KLEV-3,-1)
8771 FORMAT(1x,'Before: ',F10.2,16E12.5)

IF (ICORSFPR.EQ.0) THEN
!-- profiles are used as read
  DO JL=KIDIA,KFDIA
    PAPH5(JL,KLEV+1)=PAPH0(JL,KLEV+1)
  END DO
  
!-- half- and full-level pressures, full-level T, Q, Qs, RH

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PAPH5(JL,JK  )=AAM(JK  )+BBM(JK  )*PAPH5(JL,KLEV+1)
      PAPH5(JL,JK+1)=AAM(JK+1)+BBM(JK+1)*PAPH5(JL,KLEV+1)
      PAP5(JL,JK)=0.5*(PAPH5(JL,JK)+PAPH5(JL,JK+1))
      PT5(JL,JK) =PT0(JL,JK)
      PQ5(JL,JK) =PQ0(JL,JK)
      PQS5(JL,JK)=PQS0(JL,JK)
      PRH5(JL,JK)=PRH0(JL,JK)
    END DO
  END DO    
  
ELSE IF (ICORSFPR.EQ.1) THEN
!-- only pressure structure is adjusted
  DO JL=KIDIA,KFDIA
    PAPH5(JL,KLEV+1)=PAPH0(JL,KLEV+1)+DPR(ISTAT)*100.
  END DO
     
!-- half- and full-level pressures, full-level T, Q, Qs, RH

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PAPH5(JL,JK  )=AAM(JK  )+BBM(JK  )*PAPH5(JL,KLEV+1)
      PAPH5(JL,JK+1)=AAM(JK+1)+BBM(JK+1)*PAPH5(JL,KLEV+1)
      PAP5(JL,JK)=0.5*(PAPH5(JL,JK)+PAPH5(JL,JK+1))
      PT5(JL,JK) =PT0(JL,JK)
      PQ5(JL,JK) =PQ0(JL,JK)
      PQS5(JL,JK)=PQS0(JL,JK)
      PRH5(JL,JK)=PRH0(JL,JK)
    END DO
  END DO    
  
ELSE IF (ICORSFPR.EQ.2) THEN
!-- pressure and temperature structures are adjusted
  DO JL=KIDIA,KFDIA
    PAPH5(JL,KLEV+1)=PAPH0(JL,KLEV+1)+DPR(ISTAT)*100.
  END DO
     
!-- half- and full-level pressures

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PAPH5(JL,JK  )=AAM(JK  )+BBM(JK  )*PAPH5(JL,KLEV+1)
      PAPH5(JL,JK+1)=AAM(JK+1)+BBM(JK+1)*PAPH5(JL,KLEV+1)
      PAP5(JL,JK)=0.5*(PAPH5(JL,JK)+PAPH5(JL,JK+1))
      PQ5(JL,JK) =PQ0(JL,JK)
      PQS5(JL,JK)=PQS0(JL,JK)
      PRH5(JL,JK)=PRH0(JL,JK)
    END DO
  END DO 
  
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      PTH5(JL,JK)=(PT0(JL,JK-1)*PAP0(JL,JK-1)&
       &*(PAP0(JL,JK)-PAPH5(JL,JK))&
       &+PT0(JL,JK)*PAP0(JL,JK)*(PAPH5(JL,JK)-PAP0(JL,JK-1)))&
       &*(1./(PAPH5(JL,JK)*(PAP0(JL,JK)-PAP0(JL,JK-1))))
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    PTH5(JL,1)= PT0(JL,1)-PAP0(JL,1)*(PT0(JL,1)-PTH5(JL,2))&
      &/(PAP0(JL,1)-PAPH5(JL,2)) 
    PTH5(JL,KLEV+1)=PT0(JL,KLEV)&
      &            +(PAPH5(JL,KLEV+1)-PAP0(JL,KLEV))&
      &            *(PT0(JL,KLEV)-PTH5(JL,KLEV))&
      &            /(PAP0(JL,KLEV)-PAPH5(JL,KLEV))
  ENDDO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PT5(JL,JK)=0.5*(PTH5(JL,JK)+PTH5(JL,JK+1))
    END DO
  END DO    
    
    
ELSE IF (ICORSFPR.EQ.3) THEN
!-- pressure, temperature and humidity structures are adjusted
  DO JL=KIDIA,KFDIA
    PAPH5(JL,KLEV+1)=PAPH0(JL,KLEV+1)+DPR(ISTAT)*100.
  END DO
     
!-- half- and full-level pressures

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PAPH5(JL,JK  )=AAM(JK  )+BBM(JK  )*PAPH5(JL,KLEV+1)
      PAPH5(JL,JK+1)=AAM(JK+1)+BBM(JK+1)*PAPH5(JL,KLEV+1)
      PAP5(JL,JK)=0.5*(PAPH5(JL,JK)+PAPH5(JL,JK+1))
    END DO
  END DO    
  
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA
      PTH5(JL,JK)=(PT0(JL,JK-1)*PAP0(JL,JK-1)&
       &*(PAP0(JL,JK)-PAPH5(JL,JK))&
       &+PT0(JL,JK)*PAP0(JL,JK)*(PAPH5(JL,JK)-PAP0(JL,JK-1)))&
       &*(1./(PAPH5(JL,JK)*(PAP0(JL,JK)-PAP0(JL,JK-1))))
    ENDDO
  ENDDO
  DO JL=KIDIA,KFDIA
    PTH5(JL,1)= PT0(JL,1)-PAP0(JL,1)*(PT0(JL,1)-PTH5(JL,2))&
      &/(PAP0(JL,1)-PAPH5(JL,2)) 
    PTH5(JL,KLEV+1)=PT0(JL,KLEV)&
      &            +(PAPH5(JL,KLEV+1)-PAP0(JL,KLEV))&
      &            *(PT0(JL,KLEV)-PTH5(JL,KLEV))&
      &            /(PAP0(JL,KLEV)-PAPH5(JL,KLEV))
  ENDDO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      PT5(JL,JK)=0.5*(PTH5(JL,JK)+PTH5(JL,JK+1))
    END DO
  END DO    
  
!-- compute the model relative humidity from the new pressure and temperature
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
!---- old-style saturation function for compatibility with previous code  
      ZPP5 = PAP5(JL,JK)/100._JPRB
      ZTT5 = PT5(JL,JK)
      ZRR5 = 1._JPRB
      CALL QSAT  ( ZDDQQ, ZQQ5, ZEE, ZPP5, ZTT5, ZRR5 )
      PQS5(JL,JK) = ZQQ5
      PQ5(JL,JK) = PRH0(JL,JK)*PQS5(JL,JK)
    END DO
  END DO    

END IF


!JL=KIDIA 
!print 8772,PAPH5(JL,KLEV+1),(PQ5(JL,JK),PT5(JL,JK),PQS5(JL,JK),PRH5(JL,JK),JK=KLEV,KLEV-4,-1)
8772 FORMAT(1x,'After : ',F10.2,20E12.5)

!-- half-level pressures

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PAPH5(JL,JK)=AAM(JK)+BBM(JK)*PAPH5(JL,KLEV+1)
  END DO
END DO    
     
!-- full-level pressures

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PAP5(JL,JK)=0.5*(PAPH5(JL,JK)+PAPH5(JL,JK+1))
  END DO
END DO    

     
!-- half-level temperatures

DO JK=2,klev
  DO JL=KIDIA,KFDIA
    PTH5(JL,JK)=&
    &      (  PT5(JL,JK-1)*PAP5(JL,JK-1)*(PAP5(JL,JK)-PAPH5(JL,JK))   &                          
    &      + PT5(JL,JK)  *PAP5(JL,JK)  *(PAPH5(JL,JK)-PAP5(JL,JK-1)) )&    
    &      /(PAPH5(JL,JK)*(PAP5(JL,JK)-PAP5(JL,JK-1)))
  ENDDO     
ENDDO  

DO jl=KIDIA,KFDIA        
  PTH5(JL,1)= PT5(JL,1)-PAP5(JL,1)*(PT5(JL,1)-PTH5(JL,2))         &
  &                                /(PAP5(JL,1)-PAPH5(JL,2))  
!- There is no temperature discontinuity between surface and air 
! just above      
  PTH5(JL,klev+1)= PTS5(JL)
!- There is a temperature discontinuity between surface and air 
! just above. In that case, the air temperature is extrapolated from
! full- and half-level temperatures in the lowest layer      
!        PTH5(JL,klev+1)= PT5(JL,klev) &
!     &                +(PAPH5(JL,klev+1)-PAP5(JL,klev)) &
!     &                *(PT5(JL,klev)-PTH5(JL,klev)) &
!     &                /(PAP5(JL,klev)-PAPH5(JL,klev))
ENDDO 






END IF






DO jnu=1,KSW
  DO jl=KIDIA,KFDIA
    PALBD5(JL,JNU)=ZALBD
    PALBP5(JL,JNU)=ZALBP
  ENDDO 
ENDDO 


IF (NDUMP.LE.3) THEN
  print *,'After PTH and PAPH are created'
  JL=KIDIA
  DO JK=1,KLEV
    print 8956,PAPH5(JL,JK),PTH5(JL,JK),PAP5(JL,JK),PT5(JL,JK),PQ5(JL,JK) &
     &,PQLWP5(JL,JK),PQIWP5(JL,JK),PCLFR5(JL,JK),PRH5(JL,JK),POZON5(JL,JK)
8956 format(1x,2(F10.2,1x,F8.2),6(1x,E13.6))
  END DO
  print 8957,PAPH5(JL,KLEV+1),PTH5(JL,KLEV+1),ZALBD,ZALBP,NINDAT,ZTIM
8957 format(1x,1(F10.2,1x,F8.2),2(1x,f8.4),1x,I8,1x,F5.2) 
END IF  
      


!NSSSSS=INT(ZTIM)*3600         
!print 9022,NINDAT,NSSSSS
9022 format(1x,'Date: ',I8,' Time in seconds: ',I6)
 
MONTH=NMM(NINDAT)
IDAY =NDD(NINDAT)
IYR  =NAA(NINDAT)
RTIMTR=RTIME(IYR,MONTH,IDAY,NSSSSS)
print 9023,CSTAT,NINDAT,NSSSSS,IYR,MONTH,IDAY,RTIMTR
9023 FORMAT(1x,A3,1x,'Time: ',1X,I9,I6,I5,I3,I3,E20.12) 
if (IWRT(90).EQ.1) write(90,9023) CSTAT,NINDAT,NSSSSS,IYR,MONTH,IDAY,RTIMTR
if (IWRT(91).EQ.1) write(91,9023) CSTAT,NINDAT,NSSSSS,IYR,MONTH,IDAY,RTIMTR
if (IWRT(92).EQ.1) write(92,9023) CSTAT,NINDAT,NSSSSS,IYR,MONTH,IDAY,RTIMTR
!write(93,9023) CSTAT,NINDAT,NSSSSS,IYR,MONTH,IDAY,RTIMTR

RHGMT=ZTIM*3600.

                                           
ZTETA=RTETA(RTIMTR)
RDEASO=RRS(ZTETA)
RDECLI=RDS(ZTETA)
REQTIM=RET(ZTETA)
RSOVR =REQTIM+RHGMT
RWSOVR=RSOVR*_TWO_*RPI/RDAY
RIP0=RI0*REA*REA/(RDEASO*RDEASO)
RII0 = RIP0


RCODEC=COS(RDECLI)
RSIDEC=SIN(RDECLI)
RCOVSR=COS(RWSOVR)
RSIVSR=SIN(RWSOVR)
!print *,'RC/S-DEC, RC/S-VSR:',RCODEC,RSIDEC,RCOVSR,RSIVSR
 
 
CALL SURDI
CALL SULWN
CALL SUSWN    ( NTSW, KSW )
CALL SUAERL
CALL SUAERH
CALL SUAERSN  ( NTSW, KSW )
CALL SUCLOPN  ( NTSW, KSW, KLEV )


IF (NOVLP == 4) THEN

  RKBOL=1.380658
  RNAVO=6.0221367
  R=RNAVO*RKBOL
  RG=9.80665
  RMD=28.9644

  ZAZH(KLEV+1)=0.
  ZAZ(1)=100000.
  JL=KIDIA
  DO jk=KLEV,2,-1
    ZTVIR = PT5(JL,jk)/(1.-0.608*PQ5(jl,jk))
    ZFACT = LOG(PAPH5(jl,jk+1))-LOG(PAPH5(jl,jk))
    ZAZH(jk) = ZAZH(jk+1) + R * ZTVIR/(RMD*RG)*ZFACT
    ZAZ(jk) = 0.5*(ZAZH(jk+1)+ZAZH(jk))*1000.
!    print 9024,jk,ZAZ(jk)
9024 format(1x,I3,f12.2)    
  END DO

  CALL SUOVLP ( KLEV , ZAZ )
  
!  PRINT *,'RAOVLP,RBOVLP ',RAOVLP,RBOVLP
!  DO JK=1,KLEV
!    PRINT 9025,jk,RA1OVLP(JK)
9025 format(1x,I3,f12.9)    
!  END DO  
  
END IF  

!print *,'After all the set-up routines'

      
!     ------------------------------------------------------------------
!
       NPRINT = 1
! 
!!- basic constants 
           
RPI=3.1415926535898
RTT=273.16
REA=149597870000.
RTWAT=RTT
RTICE=RTT-23.
RDAY=86400.
RG=9.80665
RMD=28.9644
RKBOL=1.380658E-23
RNAVO=6.0221367E+23
R=RNAVO*RKBOL
RD=1000.*R/RMD
RCPD=3.5*RD
RI0=1370.
RSIGMA=5.67E-08
!      
RCDAY   = RDAY * RG / RCPD
DIFF   = 1.66
R10E   = 0.4342945
!--
PRII05 = RII0
PCCO25 = 360.E-06*44./29.
ZEPAER = 1.E-12

      
      
IF (NOZOCL.EQ.-1) THEN
  RCH4  =1.E-18
  RN2O  =1.E-18
  RO3   =1.E-18
  RCFC11=1.E-18
  RCFC12=1.E-18
END IF   
                                      
IMINUT=INT(FLOAT(NSSSSS)/60.)
!- Fortuin-Langematz O3 climatology      
CALL SUECOZC ( NINDAT , IMINUT )
!- ECMWF Geleyn O3 climatology      
  ZTHETOZ=RTETA(RTIMTR)
  ZANGOZC=REL(ZTHETOZ)-1.7535
CALL SUECOZO ( ZANGOZC )

CALL SUECAEBC
CALL SUECAEOR
CALL SUECAESD
CALL SUECAESS
CALL SUECAESU
CALL SUECAEC ( NINDAT, IMINUT )
      
       
ZDEGRAD = RPI / 180.

IF (NOZOCL.EQ.-1) THEN
  PCCO25=1.E-18
END IF   

IF (LMCCLAT .OR. NWSET.LE.0) THEN
WRITE(0,*)'Enter ZMU0 the cosine of the solar zenith angle'
  read (*,9026) ZZMU0
  print 9027,ZZMU0
END IF
9026 format(F6.4)      
9027 format(1x,'All profiles computed for Mu0= ',F6.4)
!print*,'Computations done for Mu0 going from 0.05 to 1 by 0.05'

!**======+++++++======++++++======++++++=====++++++====+++++  

IF (NDUMP.LE.3) THEN
  PRINT 9800
   9800 format(1x,'----- Geographical coordinates -----')
END IF 
DO JL=KIDIA,KFDIA
  PGELAM5(JL)=PLON5(JL)
  PGEMU5(JL)=SIN(PLAT5(JL))
  PCLON5(JL)=COS(PLON5(JL))
  PSLON5(JL)=SIN(PLON5(JL))
  IF (NDUMP.LE.3) THEN
    print 9801,JL,PLAT5(JL),PLON5(JL),PGEMU5(JL),PGELAM5(JL)&
    &   ,PCLON5(JL),PSLON5(JL)
9801 FORMAT(1X,'Lat lon (radians): ',I3,2F10.3&
    &   ,'  SinLat ',F10.6,'  LonRad ',F10.6&
    &   ,'  CosLon ',F10.6,'  SinLon ',F10.6)
  END IF 
ENDDO 
      
  DO jk=1,klev
    zeta(jk)=pap5(1,jk)/paph5(1,klev+1)
  ENDDO 
  DO jk=1,klev+1
    zetah(jk)=paph5(1,jk)/paph5(1,klev+1)
  ENDDO 
  IF (NDUMP.LT.2) THEN
    DO jk=1,klev
      print 9802,jk,zetah(jk),zeta(jk)
    ENDDO 
    jk=klev+1
    print 9802,jk,zetah(jk)
9802   format(1x,'Eta levels ',i2,2f12.9)
  END IF
 
                
!- call to other set-up routines for various coefficients
!  N.B. These are dependent on the vertical resolution
CALL SUCLD ( KLEV  , ZETA )
      if (NDUMP.LT.2) print *,'After SUCLD'

! print *,'RCAEROS ',RCAEROS      
                    
CALL SUAERV&
 & ( KLEV   ,ZETAH&
 & , CVDAES ,CVDAEL ,CVDAEU ,CVDAED&
 & , RCTRBGA,RCVOBGA,RCSTBGA,RCAEOPS,RCAEOPL,RCAEOPU&
 & , RCAEOPD,RCTRPT ,RCAEADK,RCAEADM,RCAEROS& 
 & )
 
if (NDUMP.LT.2) print *,'After SUAERV'

if (LTIMTST) then
  DO jl=KIDIA+1,KFDIA
    DO JK = 1 , klev
      pap5(jl,jk)=pap5(1,jk)
      paph5(jl,jk)=paph5(1,jk)
      pt5(jl,jk)=pt5(1,jk)
      pq5(jl,jk)=pq5(1,jk)
      pozon5(jl,jk)=pozon5(1,jk)
    ENDDO 
    paph5(jl,klev+1)=paph5(1,klev+1)
    pts5(jl)=pts5(1)
  ENDDO 
9900   FORMAT (A40)
9901   format (3x,2f10.2,f10.5,2e18.10)
9902   format (13X,f10.2,f10.5)
END IF 
                                                              
DO jnu=1,KSW
  DO jl=KIDIA,KFDIA
    PALBD5(JL,JNU)=ZALBD
    PALBP5(JL,JNU)=ZALBP
  ENDDO 
ENDDO 
!print *,'After PALBD/P '
                    
DO jl=KIDIA,KFDIA
  PEMIS5(JL)=ZEMIS
  PEMIW5(JL)=ZEMIW
  PLSM5(JL)=ZLSM
  
  PMU05(JL)=MAX( RSIDEC*PGEMU5(JL)&
   &-RCODEC*RCOVSR*SQRT(_ONE_-PGEMU5(JL)**2)*COS(PGELAM5(JL))&
   &+RCODEC*RSIVSR*SQRT(_ONE_-PGEMU5(JL)**2)*SIN(PGELAM5(JL))&
   &,_ZERO_)
   
  IF (LMCCLAT .OR. NWSET.LE.0) THEN
    PMU05(JL)=ZZMU0
  END IF   
  
  PNBAS5(JL)=1.
  PNTOP5(JL)=1.
ENDDO 
!print *,'After PEMIS/W '
      
DO jk=1,klev
  DO jl=KIDIA,KFDIA
    if (LFORCLD) then
      PCLFR5(JL,JK)=CLOUD(JK)
      PQIWP5(JL,JK)=0.
      PQLWP5(JL,JK)=0.
      PSQIW5(JL,JK)=1.
      PSQLW5(JL,JK)=1.
!    else
!      PCLFR5(JL,JK)=CF(JK)
!      PQIWP5(JL,JK)=CI(JK)
!      PQLWP5(JL,JK)=CL(JK)
    end if   
    PDP5(JL,JK)=PAPH5(JL,JK+1)-PAPH5(JL,JK) 
    IF (.NOT.LMCCLAT) THEN
      POZON5(JL,JK)=POZON5(JL,JK)*PDP5(JL,JK)
    END IF  
    PQS5(JL,JK)=0.
    PQRAIN5(JL,JK)=0.
    PRAINT5(JL,JK)=0. 
  ENDDO 
ENDDO 
!print *,'After PDP, POZON '

      
       
!- derive the aerosols and ozone distribution from climatology              



CALL  RADACA& 
 &( 1     , KLON   , KLON  , 1     , KLEV&
 &, PAPH5 , PGELAM5, PGEMU5, PCLON5, PSLON5, PTH5&
 &, ZAER5 , ZOZON5&
 &  )
!print *,'After RADACA'

 
IF (.NOT.LMCCLAT .AND. (NOZOCL.EQ.1 .OR. NOZOCL.EQ.3)) THEN
  IF (NDUMP.LT.2) then 
    print *,'Ozone profiles for the different atmospheres'
  END IF  
  DO JK=1,KLEV
    IF (NDUMP.LT.2) then 
      PRINT 9905,JK,(ZOZON5(JL,JK),JL=KIDIA,KFDIA)
    END IF  
    DO JL=KIDIA,KFDIA
      POZON5(JL,JK)=ZOZON5(JL,JK)
    ENDDO       
  ENDDO          
END IF  
 
IF (.NOT.LMCCLAT .AND. (NOZOCL.EQ.-1 .OR. NAER.EQ.0)) THEN
  DO jk=1,klev
    DO jaer=1,KAER
      DO jl=KIDIA,KFDIA
        PAER5(JL,JAER,JK)=ZEPAER
      ENDDO 
    ENDDO 
  ENDDO 
  DO jk=1,klev
    DO jl=KIDIA,KFDIA
      POZON5(JL,JK)=0.
    ENDDO 
  ENDDO 
END IF   

IF ( (LMCCLAT  .OR. NWSET.LE.0)  .AND. NAER.EQ.0) THEN
  DO jk=1,klev
    DO jaer=1,KAER
      DO jl=KIDIA,KFDIA
        PAER5(JL,JAER,JK)=ZEPAER
      ENDDO 
    ENDDO 
  ENDDO 
ELSE  
  DO jk=1,klev
    DO jaer=1,KAER
      DO jl=KIDIA,KFDIA
        PAER5(JL,JAER,JK)=ZAER5(JL,JAER,JK)
      ENDDO 
    ENDDO 
  ENDDO 
END IF  

IF (NDUMP.LT.2) then 
  print *,'Aerosol profiles for the 1st atmosphere'
  DO JK=1,KLEV
    PRINT 9905,JK,(PAER5(KIDIA,JA,JK),JA=1,KAER)
  ENDDO          
  print *,' just after RADACA'
END IF  

IF (.NOT.LMCCLAT .AND. (NOZOCL.EQ.2 .OR. NOZOCL.EQ.4) ) THEN                                       
  IF (NDUMP.LT.2) then 
    print 9904
9904    format(1x,'Entering RADOZC: FORTUIN LANGEMATZ O3 CLIMAT.') 
  END IF      
  
  CALL RADOZC &
   &( 1 , KLON , KLON  , 1, KLEV&
   &, 1    , KLON  , 0  &
   &, PAPH5, PGEMU5&
   &, ZOZON5&
   & )
   
  IF (NDUMP.LT.2) then 
    print *,'Ozone profiles for the 1st atmosphere'
  END IF      
  DO JK=1,KLEV
    IF (NDUMP.LT.2) then 
      PRINT 9905,JK,(ZOZON5(JL,JK),JL=KIDIA,KIDIA)
9905        FORMAT(1X,I3,6E12.5)
    END IF      
    DO JL=KIDIA,KFDIA
      POZON5(JL,JK)=ZOZON5(JL,JK)
    ENDDO          
  ENDDO          
END IF

IF (NOZOCL.GT.2 .OR. NAER.EQ.0) THEN
  DO jk=1,klev
    DO jaer=1,KAER
      DO jl=KIDIA,KFDIA
        PAER5(JL,JAER,JK)=ZEPAER
      ENDDO 
    ENDDO 
  ENDDO 
ELSE  
!- VOLCANIC AEROSOL SET TO epsilon IN ABSENCE OF ERUPTION
  DO jk=1,klev
    DO jl=KIDIA,KFDIA
      PAER5(JL,5,JK)=ZEPAER
    ENDDO 
  ENDDO 
END IF  
      
      
!-- SECURITY CHECK ON AEROSOL AMOUNTS
DO JK=1,KLEV
  DO JAER=1,KAER
    DO JL=KIDIA,KFDIA
      PAER5(JL,JAER,JK)=MAX(ZEPAER,PAER5(JL,JAER,JK))
    ENDDO 
  ENDDO      
ENDDO 
      
           
!-- generates pseudo-clouds (THIS IS VERY CRUDE AND SHOULD BE REPLACED 
! IN ANY REAL APPLICATION)

if (LMCCLAT) then
                                                     
!IFLAG=2
!CALL SATUR&
! & (KIDIA, KFDIA, KLON, KTDIA, KLEV, PAP5, PT5, PQS5, IFLAG )
 
DO JK=1,klev
  DO JL=KIDIA,KFDIA
    if (JL.eq.1) then
      PCLFR5(JL,JK)=CLOUD1(JK)
    else if (JL.eq.2) then
      PCLFR5(JL,JK)=CLOUD2(JK)
    else if (JL.eq.3) then
      PCLFR5(JL,JK)=CLOUD3(JK)
    else if (JL.eq.4) then
      PCLFR5(JL,JK)=CLOUD4(JK)
    else if (JL.eq.5) then
      PCLFR5(JL,JK)=CLOUD5(JK)
    end if
  
  
!---- old-style saturation function for compatibility with previous code  
     ZPP5 = PAP5(JL,JK)/100._JPRB
     ZTT5 = PT5(JL,JK)
     ZRR5 = 1._JPRB
!!     print 9006,ZPP5,ZTT5,ZRR5
!9006 format(1x,'before: P, T, R ',3F10.3)
     CALL QSAT  ( ZDDQQ, ZQQ5, ZEE, ZPP5, ZTT5, ZRR5 )
     PQS5(JL,JK) = ZQQ5
!--------------------------------------------------------------     
     
     
    if (LCLOUD.and.LFORCLD) then
! - water clouds only 
!
   
      IF (ICLDMIX.EQ.0) THEN
        PQLWP5(JL,JK) = RGAMMAS * PQS5(JL,JK)*PCLFR5(JL,JK)
        PQIWP5(JL,JK) = 0.
!         
! - ice clouds only
!
      ELSE IF (ICLDMIX.EQ.1) THEN
        PQIWP5(JL,JK) = RGAMMAS * PQS5(JL,JK)*PCLFR5(JL,JK)
        PQLWP5(JL,JK) = 0.
!         
! - mixed clouds  (mixing depends on local T)
!        
      ELSE IF (ICLDMIX.EQ.2) THEN
        ZCONDW=  RGAMMAS * PQS5(JL,JK)*PCLFR5(JL,JK)
        PTARG=PT5(JL,JK)
        ZFRACT = MIN(1.,((MAX(RTICE,MIN(RTWAT,PTARG))-RTICE)&
     &               /(RTWAT-RTICE))**2) 

        PQLWP5(JL,JK) = ZFRACT* ZCONDW
        PQIWP5(JL,JK) = (1.-ZFRACT)* ZCONDW  
      END IF   
    end if   
  ENDDO 
ENDDO 

END IF

!-- compute total cloud cover
  DO JL=KIDIA,KFDIA
    PTCC5(JL)=1.-PCLFR5(JL,1)
  END DO
  DO JK=2,KLEV
    DO JL=KIDIA,KFDIA  
      PTCC5(JL)=PTCC5(JL)*(1.-MAX(PCLFR5(JL,JK),PCLFR5(JL,JK-1))) &
      & /(1.-MIN(PCLFR5(JL,JK-1),1.-REPCLC))
    END DO
  END DO
  DO JL=KIDIA,KFDIA
    PTCC5(JL)=1.-PTCC5(JL)
  END DO
  
  IF (ICEWAT.EQ.0) THEN
!-- clouds are liquid water only  
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PQLWP5(JL,JK)=PQLWP5(JL,JK)+PQIWP5(JL,JK)
        PQIWP5(JL,JK)=0.
      END DO
    END DO    
  ELSE IF (ICEWAT.EQ.1) THEN
!-- clouds are ice water only  
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PQIWP5(JL,JK)=PQLWP5(JL,JK)+PQIWP5(JL,JK)
        PQLWP5(JL,JK)=0.
      END DO
    END DO 
  END IF     
      
      



IF (NDUMP.LE.2) THEN
  print 9098
9098 format(1x,100('*'))      
  print *,' Driver just before entering RADLSW for reference'
  print 9099
  print *,' KIDIA, KFDIA, KLON, KTDIA ',KIDIA, KFDIA, KLON, KTDIA
  print *,' KLEV, KMODE, KAER  ',KLEV, KMODE, KAER
  print *,' KLW, KSW           ',KLW, KSW
  print 9100,PRII05,PCCO25
9100 format(1x,' RI0, CO2: ',f10.3,2x,e13.7)  
  DO JL=KIDIA,KFDIA
    print 9099,JL
9099   format(1x,10('*'),I3,1x,10('*'))      
    print 9101,(PALBD5(JL,JSW),JSW=1,KSW)                        
    print 9102,(PALBP5(JL,JSW),JSW=1,KSW)   
    print 9103, PMU05(JL)
9101   FORMAT(1x,'Surf.Alb.Diffuse: ',15F10.4)                          
9102   FORMAT(1x,'Surf.Alb.Direct : ',15F10.4)                          
9103   FORMAT(1x,'Cos. Solar Zenith Angle: ',1F10.6) 
    print *,'        PAPH     PTH    PAP      PT      PDP         &
    &PQ   PCLFR    PQIWP       PQLWP       POZON       PQS'          
    DO jk=1,klev
      print 9104,jk,PAPH5(JL,JK),PTH5(JL,JK),PAP5(JL,JK),PT5(JL,JK)&
      &            ,PDP5(JL,JK)& 
      &            ,PQ5(JL,JK),PCLFR5(JL,JK),PQIWP5(JL,JK),PQLWP5(JL,JK)&
      &            ,POZON5(JL,JK),PQS5(JL,JK)
9104     format(1x,i3,f9.1,f8.2,f9.1,f8.2,f9.1,e12.5,f7.4,4e12.5)
    ENDDO 
    jk=klev+1
    print 9104,jk,PAPH5(JL,JK),PTH5(JL,JK)
    print 9105,PTS5(JL)
!    write (91,9108) jk,PAPH5(JL,JK),PTH5(JL,JK),PTS5(JL),PEMIS5(JL)&
!      & ,PEMIW5(JL),PALBD5(JL,1),PALBP5(JL,1)
9105   FORMAT(13X,f8.2)  
    print *,'Aerosols 1=Continental   2=Maritime   3=Desert     4=Urba&
    &n     5=Volcanic  6=Stratos.Bckgnd'
    DO jk=1,klev
      print 9106,jk,(PAER5(JL,JA,JK),JA=1,KAER)
9106     FORMAT(1x,'Aeros. ',I3,6E13.5)
    ENDDO   
    print 9107,PEMIS5(jl),PEMIW5(jl)
9107   FORMAT(1x,'Surf.Emis: ',2F8.4)      
  ENDDO 
  print 9098
END IF
    
IF (IWRT(90).EQ.1) THEN    
  DO JL=KIDIA,KFDIA    
  DO jk=1,klev
   write (90,9108) jk,PAPH5(JL,JK),PTH5(JL,JK),PAP5(JL,JK),PT5(JL,JK)&
    &            ,PDP5(JL,JK)& 
    &            ,PQ5(JL,JK),PCLFR5(JL,JK),PQIWP5(JL,JK),PQLWP5(JL,JK)&
    &            ,POZON5(JL,JK),PQS5(JL,JK),(PAER5(JL,JAER,JK),JAER=1,6)&
    &            ,PSQLW5(JL,JK),PSQIW5(JL,JK),PRLVRL5(JL,JK),PRLVRI5(JL,JK)
9108 format(1x,i3,24(1X,E13.6))
  ENDDO 
  jk=klev+1
  write (90,9108) jk,PAPH5(JL,JK),PTH5(JL,JK),PTS5(JL),PEMIS5(JL)&
   & ,PEMIW5(JL),PALBD5(JL,1),PALBP5(JL,1),PMU05(JL)
  END DO
END IF


! 
if (NDUMP.LE.3) print *,'Just before calling RADLSW'

CALL RADLSW &
 &( KIDIA , KFDIA , KLON  , KTDIA , KLEV   , KMODE , KAER, KBOX, NBOX &
 &, NDUMP , ILWRAD&
 &, PRII05&
 &, PAER5 , PALBD5, PALBP5, PAPH5 , PAP5&
 &, PCCO25, PCLFR5, PDP5  , PEMIS5, PEMIW5 , PLSM5 , PMU05, POZON5&
 &, PQ5   , PQIWP5, PQLWP5, PSQIW5, PSQLW5 , PQS5  , PQRAIN5, PRAINT5&
 &, PRLVRI5,PRLVRL5,PTH5  , PT5   , PTS5   , PNBAS5, PNTOP5&
!
 &, PEMIT5, PFCT5 , PFLT5 , PFCS5 , PFLS5  , PFRSOD5, PSUDU5, PUVDF5, PPARF5&
 &, PFDCT5, PFUCT5, PFDLT5, PFULT5, PFDCS5 , PFUCS5, PFDLS5, PFULS5&
 &, ASWBOX, OLRBOX, SLWBOX, SSWBOX, TAUBOX , CLDBOX &
 &)

!     ----------------------------------------------------------------

      DO JL=KIDIA,KFDIA
        IF (NDUMP.LE.3) THEN
          PRINT *,'Net Fluxes LW & SW, Total and Clear for Atm: ',JL &
          &,' i.e. Mu0= ',PMU05(JL)
        END IF  
        IF (IWRT(91).EQ.1) THEN
          WRITE (91,9991)
        END IF  
9991    format('Net Fluxes LW & SW, Total and Clear for Atm: ',I2,&              
        &,' i.e. Mu0= ',E13.6)
        DO JK=1,KLEV+1
          ZSIGMA=PAPH5(JL,JK)/PAPH5(JL,KLEV+1)
          IF (NDUMP.LE.3) THEN
            print 9992,JL,PFLT5(JL,JK),PFCT5(JL,JK),PFLS5(JL,JK) &
            & ,PFCS5(JL,JK),ZSIGMA
          END IF  
9992      format(1x,I3,4F10.3,F10.6)
          IF (IWRT(91).EQ.1) THEN
            WRITE(91,9994) PFLT5(JL,JK),PFCT5(JL,JK),PFLS5(JL,JK) &
            & ,PFCS5(JL,JK) &
            & ,PFDCT5(JL,JK),PFUCT5(JL,JK),PFDLT5(JL,JK),PFULT5(JL,JK)&
            & ,PFDCS5(JL,JK),PFUCS5(JL,JK),PFDLS5(JL,JK),PFULS5(JL,JK)&
            & ,ZSIGMA 
          END IF  
9993      FORMAT(4X,20(1X,E13.6))          
9994      FORMAT(4X,12(1X,F8.2),1X,E13.6)          
        END DO 
        
        
        IF (NDUMP.LE.3) THEN
          PRINT *,'Corresponding Heating Rates K/day  for Atm: ',JL
        END IF  
        IF (IWRT(92).EQ.1) THEN
          WRITE (92,9995)
        END IF  
9995    format('Corresponding Heating Rates K/day  for Atm: ',I2)
        DO JK=1,KLEV
          ZSIGMA=PAP5(JL,JK)/PAPH5(JL,KLEV+1)
          ZLSIGMA=LOG10(PAP5(JL,JK)/100.)
          ZDFLWT=PFLT5(JL,JK)-PFLT5(JL,JK+1)
          ZDFLWC=PFCT5(JL,JK)-PFCT5(JL,JK+1)
          ZDFSWT=PFLS5(JL,JK)-PFLS5(JL,JK+1)
          ZDFSWC=PFCS5(JL,JK)-PFCS5(JL,JK+1)
          ZDPGCP=RCDAY/(PAPH5(JL,JK+1)-PAPH5(JL,JK))
          ZCOOLT=ZDFLWT*ZDPGCP
          ZCOOLC=ZDFLWC*ZDPGCP
          ZHEATT=ZDFSWT*ZDPGCP
          ZHEATC=ZDFSWC*ZDPGCP
          IF (NDUMP.LE.3) THEN
            print 9996,JK,ZCOOLT,ZCOOLC,ZHEATT,ZHEATC,ZSIGMA 
!            & ,ZDPGCP,ZDFLWT,ZDFLWC,ZDFSWT,ZDFSWC
9996        format(1x,I3,4F10.3,F10.6,5F10.3)
!            print 9996,JK,ZCOOLC,ZSIGMA
!9996        format(1x,I3,1F10.3,F10.6,5F10.3)
          END IF  
          IF (IWRT(92).EQ.1) THEN
!            WRITE(92,9997) ZCOOLT,ZCOOLC,ZHEATT,ZHEATC,ZLSIGMA,ZSIGMA
!            IF (ITROSTRA.EQ.0) THEN
!!              WRITE(92,9997) ZCOOLC,ZSIGMA
!              WRITE(92,9997) ZCOOLT,ZSIGMA
!            ELSE IF (ITROSTRA.EQ.1) THEN
!              WRITE(92,9997) ZCOOLC,ZLSIGMA
!            ELSE
!              WRITE(92,9997) ZCOOLT,ZSIGMA,ZHEATT,ZSIGMA,ZLSIGMA
!            END IF     
            IF (ITROSTRA.EQ.0) THEN
              WRITE(92,9997) ZSIGMA,ZCOOLC,ZCOOLT
            ELSE IF (ITROSTRA.EQ.1) THEN
              WRITE(92,9997) ZLSIGMA,ZCOOLC,ZCOOLT
            ELSE IF (ITROSTRA.EQ.2) THEN
              WRITE(92,9997) ZSIGMA,ZHEATC,ZHEATT
            ELSE IF (ITROSTRA.EQ.3) THEN
              WRITE(92,9997) ZLSIGMA,ZHEATC,ZHEATT
            ELSE
              WRITE(92,9997) ZSIGMA,ZCOOLC,ZHEATC,ZCOOLT,ZHEATT
            END IF     
          END IF  
!9997      FORMAT(4X,4(1X,F8.4),2(1X,E13.6))
!9997      FORMAT(4X,1(1X,F8.4),1(1X,E13.6),1(1X,F8.4),2(1X,E13.6))
9997      FORMAT(4X,E13.6,5(1X,F8.4))
        END DO        
        
        JK=KLEV+1
!        WRITE(93,9998) ZPHOUR &
!          & ,PFLT5(JL,JK) ,PFCT5(JL,JK) ,PFLS5(JL,JK) ,PFCS5(JL,JK) &
!          & ,PFDCT5(JL,JK),PFUCT5(JL,JK),PFDLT5(JL,JK),PFULT5(JL,JK)&
!          & ,PFDCS5(JL,JK),PFUCS5(JL,JK),PFDLS5(JL,JK),PFULS5(JL,JK) 
        WRITE(93,9998) ZPHOUR &
          & ,-PFDLT5(JL,JK),PFULT5(JL,JK),PFDLS5(JL,JK),PFULS5(JL,JK)&
          & ,-PFDCT5(JL,JK),PFDCS5(JL,JK),PFULT5(JL,1),PFULS5(JL,1)&
          & , PFUCT5(JL,1),PFUCS5(JL,1), PTCC5(JL), PMU05(JL)*RII0 &
          & , PSUDU5(JL),PUVDF5(JL), PPARF5(JL)
!          & ,PFLT5(JL,JK) ,PFCT5(JL,JK) ,PFLS5(JL,JK) ,PFCS5(JL,JK) &
!          & ,PFDCT5(JL,JK),PFUCT5(JL,JK),PFDLT5(JL,JK),PFULT5(JL,JK)&
!          & ,PFDCS5(JL,JK),PFUCS5(JL,JK),PFDLS5(JL,JK),PFULS5(JL,JK) 
9998      FORMAT(1X,F10.2,15(1X,F8.3))

IF (IWRT(96).EQ.1) THEN
        if (PSUDU5(JL).GE.120.) THEN
          ZSUDUR=1.
        else
          ZSUDUR=0.
        end if    
        WRITE(96,9998) ZPHOUR &
          & ,PFDLS5(JL,JK),PFULS5(JL,JK),PFDCS5(JL,JK),PUVDF5(JL) &
          & ,PPARF5(JL),PSUDU5(JL),ZSUDUR,PTCC5(JL)
END IF          

          
ZMFDLT=ZMFDLT-PFDLT5(JL,JK)
ZMFULT=ZMFULT+PFULT5(JL,JK)
ZMFDST=ZMFDST+PFDLS5(JL,JK)
ZMFUST=ZMFUST+PFULS5(JL,JK)
ZMFDLC=ZMFDLC-PFDCT5(JL,JK)
ZMFDSC=ZMFDSC+PFDCS5(JL,JK)
ZMFDSI=ZMFDSI+PMU05(JL)*RII0
ITIMES=ITIMES+1

IF (PFDCS5(JL,JK).GT.0.) THEN
  ZMFDSTD=ZMFDSTD+PFDLS5(JL,JK)
  ZMFUSTD=ZMFUSTD+PFULS5(JL,JK)
  ZMFDSCD=ZMFDSCD+PFDCS5(JL,JK)
  ITIMSW=ITIMSW+1
END IF

IF (PFDLT5(JL,JK).NE.PFDCT5(JL,JK)) THEN
  ZMFDLTC=ZMFDLTC-PFDLT5(JL,JK)
  ITIMLWC=ITIMLWC+1
END IF  
        
IF (PFDCS5(JL,JK).GT.0. .AND. (PFDCS5(JL,JK).NE.PFDLS5(JL,JK)) ) THEN
  ZMFDSTC=ZMFDSTC+PFDLS5(JL,JK)
  ITIMSWC=ITIMSWC+1
END IF  

IF (IWRT(94).EQ.1) THEN
  DO ICBOX=1,NBOX
    WRITE(94,1111) ZPHOUR,TAUBOX(JL,ICBOX),ASWBOX(JL,ICBOX) &
    &  , OLRBOX(JL,ICBOX),SLWBOX(JL,ICBOX),SSWBOX(JL,ICBOX)
1111 format(1x,f10.2,5(1x,F10.4))
  END DO
END IF

IF (IWRT(95).EQ.1) THEN
  WRITE(95,1112) ZPHOUR
1112 format(1x,f10.2,50('-'))
  DO JK=1,KLEV
    WRITE(95,1113) (CLDBOX(JL,ICBOX,JK),ICBOX=1,NBOX)
1113 format(1x,100f2.0)
  END DO
END IF

      
      END DO    
      
!      print *,'About to go back to 1000'
      GO TO 1000         
     
!     ----------------------------------------------------------------
 
9999  CONTINUE

ZMFDLT=ZMFDLT/ITIMES
ZMFULT=ZMFULT/ITIMES
ZMFDST=ZMFDST/ITIMES
ZMFUST=ZMFUST/ITIMES
ZMFDLC=ZMFDLC/ITIMES
ZMFDSC=ZMFDSC/ITIMES
ZMFDSI=ZMFDSI/ITIMES

ZMFDSTD=ZMFDSTD/ITIMSW
ZMFUSTD=ZMFUSTD/ITIMSW
ZMFDSCD=ZMFDSCD/ITIMSW

ZMFDLTC=ZMFDLTC/ITIMLWC
ZMFDSTC=ZMFDSTC/ITIMSWC

WRITE(93,9998) 99999.99 &
  & ,ZMFDLT,ZMFULT,ZMFDST,ZMFUST,ZMFDLC,ZMFDSC &
  & ,ZMFDSTD,ZMFUSTD,ZMFDSCD,ZMFDLTC,ZMFDSTC,ZMFDSI

   
END DO
!-- end of loop on ISTAT, the station index      
           

STOP

END PROGRAM RAD1DRIV
