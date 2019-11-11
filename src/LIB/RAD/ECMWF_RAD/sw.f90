SUBROUTINE SW &
 &( KIDIA, KFDIA , KLON  , KLEV , KAER &
 &, PSCT , PCARDI, PPSOL , PALBD, PALBP , PWV, PQS &
 &, PRMU0, PCG   , PCLDSW, PDP  , POMEGA, POZ, PPMB &
 &, PTAU , PTAVE , PAER &
 &, PHEAT, PFDOWN, PFUP  &
 &, PCEAT, PCDOWN, PCUP  &
 &, PFDNN, PFDNV , PFUPN, PFUPV &
 &, PCDNN, PCDNV , PCUPN, PCUPV &
 &, PSUDU, PUVDF , PPARF, PDIFFS, PDIRFS &
 &, ODUST, PPIZA_DST, PCGA_DST, PTAUREL_DST  )

!**** *SW* - COMPUTES THE SHORTWAVE RADIATION FLUXES.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW* IS CALLED FROM *RADLSW*


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
!          2. COMPUTES FLUXES IN U.V./VISIBLE  SPECTRAL INTERVAL (SW1S)
!          3. COMPUTES FLUXES IN NEAR-INFRARED SPECTRAL INTERVAL (SWNI)

!     EXTERNALS.
!     ----------

!          *SWU*, *SW1S*, *SWNI*

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
!        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
!        95-12-07   J.-J. MORCRETTE  Near-Infrared in nsw-1 Intervals
!        990128     JJMorcrette      sunshine duration
!        99-05-25   JJMorcrette      Revised aerosols
!        00-12-18   JJMorcrette      6 spectral intervals

!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE OYOERAD   , ONLY : NSW
USE YOERDU   , ONLY : RCDAY
!++MODIF_MESONH
USE OYOERAD   , ONLY : NOVLP
!--MODIF_MESONH
USE MODI_SWU
USE MODI_SW1S
USE MODI_SWNI
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KAER
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


REAL_B :: PPSOL(KLON), PAER(KLON,6,KLEV),PRMU0(KLON)&
  &,  PWV(KLON,KLEV),PQS(KLON,KLEV)

REAL_B :: PALBD(KLON,NSW)      , PALBP(KLON,NSW)&
  &,  PCG(KLON,NSW,KLEV)   , PCLDSW(KLON,KLEV)&
  &,  PDP(KLON,KLEV)  &
  &,  POMEGA(KLON,NSW,KLEV), POZ(KLON,KLEV)&
  &,  PPMB(KLON,KLEV+1)&
  &,  PTAU(KLON,NSW,KLEV)  , PTAVE(KLON,KLEV)

REAL_B :: PHEAT(KLON,KLEV), PFDOWN(KLON,KLEV+1), PFUP(KLON,KLEV+1),&
     &PFUPV(KLON), PFUPN(KLON), PFDNV(KLON), PFDNN(KLON)&
  &,  PCEAT(KLON,KLEV), PCDOWN(KLON,KLEV+1), PCUP(KLON,KLEV+1)&
  &,  PCUPV(KLON), PCUPN(KLON), PCDNV(KLON), PCDNN(KLON)&
  &,  PSUDU(KLON), PUVDF(KLON), PPARF(KLON), PDIFFS(KLON,NSW)&
  &,  PDIRFS(KLON,NSW)

!++MODIF_MESONH
LOGICAL :: ODUST                       ! flag for dust
REAL_B  :: PPIZA_DST(KLON,KLEV,NSW)    !wvl dependent ssa dust
REAL_B  :: PCGA_DST(KLON,KLEV,NSW)     !wvl dependent assym.fact dust
REAL_B  :: PTAUREL_DST(KLON,KLEV,NSW)  !wvl dependent tau/tau_{550} for dust
!--MODIF_MESONH

!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZAKI(KLON,2,NSW)&
  &,  ZCLD(KLON,KLEV)    , ZCLEAR(KLON) &
  &,  ZDSIG(KLON,KLEV)   , ZFACT(KLON)&
  &,  ZFD(KLON,KLEV+1)   , ZCD(KLON,KLEV+1)&
  &,  ZCDOWN(KLON,KLEV+1), ZCDNIR(KLON,KLEV+1), ZCDUVS(KLON,KLEV+1)&
  &,  ZFDOWN(KLON,KLEV+1), ZFDNIR(KLON,KLEV+1), ZFDUVS(KLON,KLEV+1)&
  &,  ZFU(KLON,KLEV+1)   , ZCU(KLON,KLEV+1)&
  &,  ZCUP(KLON,KLEV+1)  , ZCUNIR(KLON,KLEV+1), ZCUUVS(KLON,KLEV+1)&
  &,  ZFUP(KLON,KLEV+1)  , ZFUNIR(KLON,KLEV+1), ZFUUVS(KLON,KLEV+1)&
  &,  ZRMU(KLON)         , ZSEC(KLON)         &
  &,  ZSUDU1(KLON)       , ZSUDU2(KLON)       &
  &,  ZSUDU1T(KLON)      , ZSUDU2T(KLON)      &
  &,  ZUD(KLON,5,KLEV+1) , ZDIFF(KLON,KLEV)  ,ZDIRF(KLON,KLEV)&
  &,  ZDIFF2(KLON,KLEV)  , ZDIRF2(KLON,KLEV)

!     LOCAL INTEGER SCALARS
INTEGER_M :: INU, JK, JKL, JL, JNU, INUVS, INUIR

!     LOCAL REAL SCALARS
REAL_B :: ZDCNET, ZDFNET

!     ------------------------------------------------------------------

!*         1.     ABSORBER AMOUNTS AND OTHER USEFUL QUANTITIES
!                 --------------------------------------------


CALL SWU ( KIDIA,KFDIA ,KLON  ,KLEV &
         &, PSCT ,PCARDI,PCLDSW,PPMB ,PPSOL &
         &, PRMU0,PTAVE ,PWV &
         &, ZAKI ,ZCLD  ,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD )
         

!     ------------------------------------------------------------------

!*         2.     INTERVAL (0.185/0.25-0.68 MICRON): U.V. AND VISIBLE
!                 ---------------------------------------------------

IF (NSW.LE.4) THEN
  INUVS=1
  INUIR=2
ELSE IF (NSW.EQ.6) THEN
  INUVS=1
  INUIR=4
END IF     

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    ZFD(JL,JK) =_ZERO_
    ZFU(JL,JK) =_ZERO_
    ZCD(JL,JK) =_ZERO_
    ZCU(JL,JK) =_ZERO_
    ZSUDU1T(JL)=_ZERO_
    PUVDF(JL)  =_ZERO_
    PPARF(JL)  =_ZERO_
  ENDDO
ENDDO

DO JNU = INUVS , INUIR-1
   
   !++MODIF_MESONH
     CALL SW1S &
           &( KIDIA , KFDIA, KLON , KLEV , KAER  , JNU &
           &,  PAER , PALBD , PALBP, PCG  , ZCLD , ZCLEAR &
           &,  ZDSIG, POMEGA, POZ  , ZRMU , ZSEC , PTAU  , ZUD  &
           &,  ZFDUVS,ZFUUVS, ZCDUVS,ZCUUVS, ZSUDU1, ZDIFF,ZDIRF &
           &,  ODUST,PPIZA_DST(:,:,JNU) &       ! SSA for this wavelength
           &,  PCGA_DST(:,:,JNU)   &            ! GCA for this wavelengt
           &,  PTAUREL_DST(:,:,JNU) )           ! TAUREL for this wavelength
   !--MODIF_MESONH

  DO JL=KIDIA,KFDIA
    PDIFFS(JL,JNU)=ZDIFF(JL,1)*ZFACT(JL)
    PDIRFS(JL,JNU)=ZDIRF(JL,1)*ZFACT(JL)
  ENDDO
  DO JK = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZFD(JL,JK)=ZFD(JL,JK)+ZFDUVS(JL,JK)
      ZFU(JL,JK)=ZFU(JL,JK)+ZFUUVS(JL,JK)
      ZCD(JL,JK)=ZCD(JL,JK)+ZCDUVS(JL,JK)
      ZCU(JL,JK)=ZCU(JL,JK)+ZCUUVS(JL,JK)
    ENDDO
  ENDDO
  DO JL = KIDIA,KFDIA
    ZSUDU1T(JL)=ZSUDU1T(JL)+ZSUDU1(JL)
  ENDDO
  
  IF (NSW.EQ.6) THEN
    IF (JNU.LT.INUIR-1) THEN
      DO JL=KIDIA,KFDIA
        PUVDF(JL)=PUVDF(JL)+ZFDUVS(JL,1)
      END DO
    ELSE     
      DO JL=KIDIA,KFDIA
        PPARF(JL)=PPARF(JL)+ZFDUVS(JL,1)
      END DO
    END IF
  END IF    
  
ENDDO
!     ------------------------------------------------------------------

!*         3.     INTERVAL (0.68-4.00 MICRON): NEAR-INFRARED
!                 ------------------------------------------


DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    ZFDOWN(JL,JK)=_ZERO_
    ZFUP  (JL,JK)=_ZERO_
    ZCDOWN(JL,JK)=_ZERO_
    ZCUP  (JL,JK)=_ZERO_
    ZSUDU2T(JL)  =_ZERO_
  ENDDO
ENDDO

! PRINT*,"sw  PDIFFS,PDIRFS"
! 
! PRINT*,PDIFFS
! 
! PRINT*,PDIRFS

DO JNU = INUIR , NSW

   !++MODIF_MESONH
      CALL SWNI &
           &(  KIDIA ,KFDIA , KLON , KLEV , KAER , JNU &
           &,  PAER  ,ZAKI  , PALBD, PALBP, PCG  , ZCLD, ZCLEAR &
           &,  ZDSIG ,POMEGA, POZ  , ZRMU , ZSEC , PTAU, ZUD      &
           &,  PWV   ,PQS &
           &,  ZFDNIR,ZFUNIR,ZCDNIR,ZCUNIR,ZSUDU2,ZDIFF2,ZDIRF2 &
           &,  ODUST,PPIZA_DST(:,:,JNU)  &
           &,  PCGA_DST(:,:,JNU)    &
           &,  PTAUREL_DST(:,:,JNU) &
           &)
    !--MODIF_MESONH
    
    

  DO JL=KIDIA,KFDIA
!     PRINT*,JL,JNU
!     PRINT*,"SW"
!     PRINT*,ZDIFF2(JL,1)
    
    PDIFFS(JL,JNU)=ZDIFF2(JL,1)*ZFACT(JL)
    PDIRFS(JL,JNU)=ZDIRF2(JL,1)*ZFACT(JL)
  ENDDO
  DO JK = 1 , KLEV+1
    DO JL = KIDIA,KFDIA
      ZFDOWN(JL,JK)=ZFDOWN(JL,JK)+ZFDNIR(JL,JK)
      ZFUP  (JL,JK)=ZFUP  (JL,JK)+ZFUNIR(JL,JK)
      ZCDOWN(JL,JK)=ZCDOWN(JL,JK)+ZCDNIR(JL,JK)
      ZCUP  (JL,JK)=ZCUP  (JL,JK)+ZCUNIR(JL,JK)
    ENDDO
  ENDDO
  DO JL = KIDIA,KFDIA
    ZSUDU2T(JL)=ZSUDU2T(JL)+ZSUDU2(JL)
  ENDDO
ENDDO
! 
! PRINT*,"sw  PDIFFS,PDIRFS"
! 
! PRINT*,PDIFFS
! 
! PRINT*,PDIRFS
! 
! pause

!     ------------------------------------------------------------------

!*         4.     FILL THE DIAGNOSTIC ARRAYS
!                 --------------------------


DO JL = KIDIA,KFDIA
  PFDNN(JL)=ZFDOWN(JL,1)*ZFACT(JL)
  PFDNV(JL)=ZFD(JL,1)*ZFACT(JL)
  PFUPN(JL)=ZFUP(JL,KLEV+1)*ZFACT(JL)
  PFUPV(JL)=ZFU(JL,KLEV+1)*ZFACT(JL)

  PCDNN(JL)=ZCDOWN(JL,1)*ZFACT(JL)
  PCDNV(JL)=ZCD(JL,1)*ZFACT(JL)
  PCUPN(JL)=ZCUP(JL,KLEV+1)*ZFACT(JL)
  PCUPV(JL)=ZCU(JL,KLEV+1)*ZFACT(JL)

  PSUDU(JL)=(ZSUDU1T(JL)+ZSUDU2T(JL))*ZFACT(JL) 
  PUVDF(JL)=PUVDF(JL)*ZFACT(JL)
  PPARF(JL)=PPARF(JL)*ZFACT(JL)
ENDDO

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PFUP(JL,JK)   = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
    PFDOWN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
    PCUP(JL,JK)   = (ZCUP(JL,JK)   + ZCU(JL,JK)) * ZFACT(JL)
    PCDOWN(JL,JK) = (ZCDOWN(JL,JK) + ZCD(JL,JK)) * ZFACT(JL)
  ENDDO
ENDDO

DO JKL = 1 , KLEV
  JK = KLEV+1 - JKL
  DO JL = KIDIA,KFDIA
    ZDFNET = PFUP(JL,JK+1) - PFDOWN(JL,JK+1)-PFUP(JL,JK  ) + PFDOWN(JL,JK  )
    PHEAT(JL,JK) = RCDAY * ZDFNET / PDP(JL,JKL)
    ZDCNET = PCUP(JL,JK+1) - PCDOWN(JL,JK+1)-PCUP(JL,JK  ) + PCDOWN(JL,JK  )
    PCEAT(JL,JK) = RCDAY * ZDCNET / PDP(JL,JKL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SW



