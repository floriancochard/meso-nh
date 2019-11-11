!OPTIONS XOPT(HSFUN)
SUBROUTINE SWU &
  &( KIDIA, KFDIA , KLON  , KLEV &
  &, PSCT , PCARDI, PCLDSW, PPMB , PPSOL, PRMU0, PTAVE, PWV &
  &, PAKI , PCLD  , PCLEAR, PDSIG, PFACT, PRMU , PSEC , PUD &
  &)

!**** *SWU* - SHORTWAVE RADIATION, ABSORBER AMOUNTS

!     PURPOSE.
!     --------
!           COMPUTES THE ABSORBER AMOUNTS USED IN SHORTWAVE RADIATION
!     CALCULATIONS

!**   INTERFACE.
!     ----------
!          *SWU* IS CALLED BY *SW*


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES ABSORBER AMOUNTS WITH TEMPERATURE AND PRESSURE
!     SCALING.

!     EXTERNALS.
!     ----------

!          *SWTT*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14

!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE YOECLD   , ONLY : REPSEC
USE OYOERAD   , ONLY : NOVLP    ,NSW
USE YOERDU   , ONLY : REPSCQ
USE OYOESW    , ONLY : RPDH1    ,RPDU1    ,RPNH     ,RPNU     ,&
            &RTDH2O   ,RTDUMG   ,RTH2O    ,RTUMG
USE YOEOVLP  , ONLY : RA1OVLP
!
USE MODI_SWTT1
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON

!     DUMMY REAL SCALARS
REAL_B :: PCARDI
REAL_B :: PSCT



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PCLDSW(KLON,KLEV), PPMB(KLON,KLEV+1), PPSOL(KLON)&
  &,  PRMU0(KLON)      , PTAVE(KLON,KLEV) , PWV(KLON,KLEV)

REAL_B :: PAKI(KLON,2,NSW)&
  &,  PCLD(KLON,KLEV)  , PCLEAR(KLON)&
  &,  PDSIG(KLON,KLEV) , PFACT(KLON)      , PRMU(KLON)&
  &,  PSEC(KLON)       , PUD(KLON,5,KLEV+1)
  
INTEGER_M :: INUIR  

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

INTEGER_M :: IIND(2)
REAL_B :: ZC1J(KLON,KLEV+1),ZCLEAR(KLON),ZCLOUD(KLON)&
  &,  ZN175(KLON), ZN190(KLON), ZO175(KLON)&
  &,  ZO190(KLON), ZSIGN(KLON)&
  &,  ZR(KLON,2) , ZSIGO(KLON), ZUD(KLON,2)

!     LOCAL INTEGER SCALARS
INTEGER_M :: JA, JK, JKL, JKLP1, JKP1, JL, JNU

!     LOCAL REAL SCALARS
REAL_B :: ZDSCO2, ZDSH2O, ZFPPW, ZRTH, ZRTU, ZWH2O, ZALPHA1


!     ------------------------------------------------------------------

!*         1.     COMPUTES AMOUNTS OF ABSORBERS
!                 -----------------------------


IIND(1)=1
IIND(2)=2


!*         1.1    INITIALIZES QUANTITIES
!                 ----------------------


DO JL = KIDIA,KFDIA
  PUD(JL,1,KLEV+1)=_ZERO_
  PUD(JL,2,KLEV+1)=_ZERO_
  PUD(JL,3,KLEV+1)=_ZERO_
  PUD(JL,4,KLEV+1)=_ZERO_
  PUD(JL,5,KLEV+1)=_ZERO_
  PFACT(JL)= PRMU0(JL) * PSCT
! MODIF
!- already accounted for in RADINT      
  PRMU(JL)=SQRT(1224.* PRMU0(JL) * PRMU0(JL) + 1.) / 35.
! PRMU(JL)=PRMU0(JL)
  PSEC(JL)=_ONE_/PRMU(JL)
  ZC1J(JL,KLEV+1)=_ZERO_
ENDDO

!*          1.3    AMOUNTS OF ABSORBERS
!                  --------------------


DO JL= KIDIA,KFDIA
  ZUD(JL,1) = _ZERO_
  ZUD(JL,2) = _ZERO_
  ZO175(JL) = PPSOL(JL)** RPDU1
  ZO190(JL) = PPSOL(JL)** RPDH1
  ZSIGO(JL) = PPSOL(JL)
  ZCLEAR(JL)=_ONE_
  ZCLOUD(JL)=_ZERO_
ENDDO

DO JK = 1 , KLEV
  JKP1 = JK + 1
  JKL = KLEV+1 - JK
  JKLP1 = JKL+1
  DO JL = KIDIA,KFDIA
    ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
    ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
    ZWH2O = MAX (PWV(JL,JKL) , REPSCQ )
    ZSIGN(JL) = 100._JPRB * PPMB(JL,JKP1)
    PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN(JL))/PPSOL(JL)
    ZN175(JL) = ZSIGN(JL) ** RPDU1
    ZN190(JL) = ZSIGN(JL) ** RPDH1
    ZDSCO2 = ZO175(JL) - ZN175(JL)
    ZDSH2O = ZO190(JL) - ZN190(JL)
    PUD(JL,1,JK) = RPNH * ZDSH2O * ZWH2O  * ZRTH
    PUD(JL,2,JK) = RPNU * ZDSCO2 * PCARDI * ZRTU
    ZFPPW=1.6078_JPRB*ZWH2O/(_ONE_+0.608_JPRB*ZWH2O)
    PUD(JL,4,JK)=PUD(JL,1,JK)*ZFPPW
    PUD(JL,5,JK)=PUD(JL,1,JK)*(_ONE_-ZFPPW)
    ZUD(JL,1) = ZUD(JL,1) + PUD(JL,1,JK)
    ZUD(JL,2) = ZUD(JL,2) + PUD(JL,2,JK)
    ZSIGO(JL) = ZSIGN(JL)
    ZO175(JL) = ZN175(JL)
    ZO190(JL) = ZN190(JL)

!++MODIF_MESONH
    IF ((NOVLP == 1).OR.(NOVLP==6).OR.(NOVLP==8)) THEN
!--MODIF_MESONH
      ZCLEAR(JL)=ZCLEAR(JL)&
       &*(_ONE_-MAX(PCLDSW(JL,JKL),ZCLOUD(JL)))&
       &/(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPSEC))
      ZC1J(JL,JKL)= _ONE_ - ZCLEAR(JL)
      ZCLOUD(JL) = PCLDSW(JL,JKL)
    ELSEIF ((NOVLP == 2).OR.(NOVLP==7)) THEN
      ZCLOUD(JL) = MAX(PCLDSW(JL,JKL),ZCLOUD(JL))
      ZC1J(JL,JKL) = ZCLOUD(JL)
!++MODIF_MESONH
    ELSEIF ((NOVLP == 3).OR.(NOVLP==5)) THEN
!--MODIF_MESONH
      ZCLEAR(JL) = ZCLEAR(JL)*(_ONE_-PCLDSW(JL,JKL))
      ZCLOUD(JL) = _ONE_ - ZCLEAR(JL)
      ZC1J(JL,JKL) = ZCLOUD(JL)
    ELSEIF (NOVLP == 4) THEN
!** Hogan & Illingworth (2001)      
      ZALPHA1=RA1OVLP(KLEV+1-JK)
      ZCLEAR(JL)=ZCLEAR(JL)*( &
        & ZALPHA1*(_ONE_-MAX(PCLDSW(JL,JKL),ZCLOUD(JL))) &
        &        /(_ONE_-MIN(ZCLOUD(JL),_ONE_-REPSEC)) &
        & +(_ONE_-ZALPHA1)*(_ONE_-PCLDSW(JL,JKL)) )
      ZC1J(JL,JKL) = _ONE_ - ZCLEAR(JL) 
      ZCLOUD(JL) = PCLDSW(JL,JKL)
    ENDIF
  ENDDO
ENDDO
DO JL=KIDIA,KFDIA
  PCLEAR(JL)=_ONE_-ZC1J(JL,1)
ENDDO
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    IF (PCLEAR(JL) < _ONE_) THEN
      PCLD(JL,JK)=PCLDSW(JL,JK)/(_ONE_-PCLEAR(JL))
    ELSE
      PCLD(JL,JK)=_ZERO_
    ENDIF
  ENDDO
ENDDO


!*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
!                 -----------------------------------------------


DO JA = 1,2
  DO JL = KIDIA,KFDIA
    ZUD(JL,JA) = ZUD(JL,JA) * PSEC(JL)
  ENDDO
ENDDO

IF (NSW.LE.4) THEN
  INUIR=2
ELSE IF (NSW.EQ.6) THEN
  INUIR=4
END IF     


DO JNU= INUIR,NSW

  CALL SWTT1 ( KIDIA,KFDIA,KLON, JNU, 2, IIND &
   &, ZUD &
   &, ZR                            )

  DO JA = 1,2
    DO JL = KIDIA,KFDIA
      PAKI(JL,JA,JNU) = -LOG( ZR(JL,JA) ) / ZUD(JL,JA)
    ENDDO
  ENDDO
ENDDO


!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SWU

