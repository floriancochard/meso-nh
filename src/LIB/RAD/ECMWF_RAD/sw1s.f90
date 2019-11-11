SUBROUTINE SW1S &
 &( KIDIA , KFDIA , KLON , KLEV , KAER , KNU &
 &, PAER  , PALBD , PALBP, PCG  , PCLD , PCLEAR &
 &, PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  &
 &, PFD   , PFU   , PCD  , PCU  , PSUDU1,PDIFF,PDIRF &
!++MODIF_MESONH
 &, ODUST,PPIZA_DST,PCGA_DST,PTAUREL_DST  &
!--MODIF_MESONH
 &)

!**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
!     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

!**   INTERFACE.
!     ----------

!          *SW1S* IS CALLED FROM *SW*.


!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!          1. COMPUTES QUANTITIES FOR THE CLEAR-SKY FRACTION OF THE
!     COLUMN
!          2. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
!     CONTINUUM SCATTERING
!          3. MULTIPLY BY OZONE TRANSMISSION FUNCTION

!     EXTERNALS.
!     ----------

!          *SWCLR*, *SWR*, *SWTT*, *SWUVO3*

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
!        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
!        96-01-15   J.-J. MORCRETTE    SW in nsw SPECTRAL INTERVALS 
!        990128     JJMorcrette        sunshine duration
!        99-05-25   JJMorcrette        Revised aerosols
!        00-12-18   JJMorcrette        6 spectral intervals

!     ------------------------------------------------------------------


#include "tsmbkind.h"

USE OYOESW    , ONLY : RRAY     ,RSUN
!++MODIF_MESONH
USE OYOERAD   , ONLY : NSW , NOVLP
!--MODIF_MESONH
USE MODI_SWCLR
USE MODI_SWR
USE MODI_SWTT1
USE MODI_SWUVO3
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KAER
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KNU



!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

REAL_B :: PAER(KLON,6,KLEV)&
  &,  PALBD(KLON,NSW)      , PALBP(KLON,NSW)&
  &,  PCG(KLON,NSW,KLEV)   , PCLD(KLON,KLEV) &
  &,  PCLEAR(KLON)&
  &,  PDSIG(KLON,KLEV)&
  &,  POMEGA(KLON,NSW,KLEV), POZ(KLON,KLEV)&
  &,  PRMU(KLON)           , PSEC(KLON)&
  &,  PTAU(KLON,NSW,KLEV)  , PUD(KLON,5,KLEV+1)

REAL_B :: PFD(KLON,KLEV+1) , PFU(KLON,KLEV+1)&
  &,  PCD(KLON,KLEV+1)     , PCU(KLON,KLEV+1)&
  &,  PSUDU1(KLON)         , PDIFF(KLON,KLEV)&
  &,  PDIRF(KLON,KLEV)

!++MODIF_MESONH
LOGICAL          :: ODUST                   ! flag for DUST
REAL_B  :: PPIZA_DST(KLON,KLEV)    !wvl dependent ssa dust
REAL_B  :: PCGA_DST(KLON,KLEV)     !wvl dependent assym.fact dust
REAL_B  :: PTAUREL_DST(KLON,KLEV)  !wvl dependent tau/tau_{550} for dust
!--MODIF_MESONH


!     ------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

INTEGER_M :: IIND6(6), IIND4(4)

REAL_B :: ZCGAZ(KLON,KLEV)&
  &,  ZDIFF(KLON)        , ZDIRF(KLON)        &
  &,  ZDIFT(KLON)        , ZDIRT(KLON)        &
  &,  ZPIZAZ(KLON,KLEV)&
  &,  ZRAYL(KLON), ZRAY1(KLON,KLEV+1), ZRAY2(KLON,KLEV+1)&
  &,  ZREFZ(KLON,2,KLEV+1)&
  &,  ZRJ(KLON,6,KLEV+1), ZRJ0(KLON,6,KLEV+1)&
  &,  ZRK(KLON,6,KLEV+1), ZRK0(KLON,6,KLEV+1)&
  &,  ZRMUE(KLON,KLEV+1), ZRMU0(KLON,KLEV+1)&
  &,  ZR6(KLON,6)       , ZR4(KLON,4)&
  &,  ZTAUAZ(KLON,KLEV)&
  &,  ZTRA1(KLON,KLEV+1), ZTRA2(KLON,KLEV+1)&
  &,  ZTRCLD(KLON)      , ZTRCLR(KLON)&
  &,  ZW6(KLON,6)       , ZW4(KLON,4), ZO(KLON,2) ,ZT(KLON,2) 

REAL_B :: ZTA1(KLON), ZTO1(KLON)
REAL_B :: ZCLDIR   
  
!     LOCAL INTEGER SCALARS
INTEGER_M :: IKL, IKM1, JAJ, JK, JL




!     ------------------------------------------------------------------

!*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
!                 ----------------------- ------------------


!*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
!                 -----------------------------------------

DO JL = KIDIA,KFDIA
  ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)&
   &* (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)&
   &* (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
ENDDO


!     ------------------------------------------------------------------

!*         2.    CONTINUUM SCATTERING CALCULATIONS
!                ---------------------------------


!*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
!                --------------------------------

!++MODIF_MESONH
   CALL SWCLR &
        &( KIDIA  , KFDIA , KLON  , KLEV , KAER , KNU &
        &, PAER   , PALBP , PDSIG , ZRAYL, PSEC &
        &, ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
        &, ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2, ZTRCLR &
        &, ODUST  , PPIZA_DST,PCGA_DST  &
       &, PTAUREL_DST )
!--MODIF_MESONH

!*         2.2   CLOUDY FRACTION OF THE COLUMN
!                -----------------------------


CALL SWR &
  &( KIDIA ,KFDIA ,KLON  ,KLEV  , KNU &
  &, PALBD ,PCG   ,PCLD  ,POMEGA, PSEC , PTAU &
  &, ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 , ZREFZ, ZRJ  ,ZRK , ZRMUE &
  &, ZTAUAZ,ZTRA1 ,ZTRA2 ,ZTRCLD &
  &)

!     ------------------------------------------------------------------

!*         3.    OZONE ABSORPTION
!                ----------------

IF (NSW <= 4) THEN

!*         3.1   TWO OR FOUR SPECTRAL INTERVALS
!                ------------------------------

  IIND6(1)=1
  IIND6(2)=2
  IIND6(3)=3
  IIND6(4)=1
  IIND6(5)=2
  IIND6(6)=3


!*         3.1.1  DOWNWARD FLUXES
!                 ---------------


  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW6(JL,1)=_ZERO_
    ZW6(JL,2)=_ZERO_
    ZW6(JL,3)=_ZERO_
    ZW6(JL,4)=_ZERO_
    ZW6(JL,5)=_ZERO_
    ZW6(JL,6)=_ZERO_
    PFD(JL,KLEV+1)=((_ONE_-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
     &+ PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
  ENDDO
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZW6(JL,1)=ZW6(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW6(JL,2)=ZW6(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
      ZW6(JL,3)=ZW6(JL,3)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW6(JL,4)=ZW6(JL,4)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW6(JL,5)=ZW6(JL,5)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
      ZW6(JL,6)=ZW6(JL,6)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ENDDO
    
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 6 &
      &, IIND6 &
      &, ZW6  &
      &, ZR6                          )

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR6(JL,1)*ZR6(JL,2)*ZR6(JL,3)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR6(JL,4)*ZR6(JL,5)*ZR6(JL,6)*ZRJ0(JL,JAJ,IKL)
      PDIFF(JL,IKL) = ZDIFF(JL)*RSUN(KNU)*(_ONE_-PCLEAR(JL))
      PDIRF(JL,IKL) = ZDIRF(JL) * RSUN(KNU)*PCLEAR(JL)
      PFD(JL,IKL) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
       &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    ZDIFT(JL) = ZR6(JL,1)*ZR6(JL,2)*ZR6(JL,3)*ZTRCLD(JL)
    ZDIRT(JL) = ZR6(JL,4)*ZR6(JL,5)*ZR6(JL,6)*ZTRCLR(JL)
    PSUDU1(JL) = ((_ONE_-PCLEAR(JL)) * ZDIFT(JL)&
     &+PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)
  ENDDO


!*         3.1.2  UPWARD FLUXES
!                 -------------


  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((_ONE_-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
     &+ PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
     &* RSUN(KNU)
    PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW6(JL,1)=ZW6(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW6(JL,2)=ZW6(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW6(JL,3)=ZW6(JL,3)+POZ(JL,  IKM1)*1.66_JPRB
      ZW6(JL,4)=ZW6(JL,4)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW6(JL,5)=ZW6(JL,5)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW6(JL,6)=ZW6(JL,6)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO
    
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 6 &
      &, IIND6 &
      &, ZW6  &
      &, ZR6                          )
  
    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR6(JL,1)*ZR6(JL,2)*ZR6(JL,3)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR6(JL,4)*ZR6(JL,5)*ZR6(JL,6)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
       &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO




ELSE IF (NSW == 6) THEN

!*         3.2   SIX SPECTRAL INTERVALS
!                ----------------------

  IIND4(1)=1
  IIND4(2)=2
  IIND4(3)=1
  IIND4(4)=2


!*         3.2,1  DOWNWARD FLUXES
!                 ---------------

  JAJ = 2

  DO JL = KIDIA,KFDIA
    ZW4(JL,1)=_ZERO_
    ZW4(JL,2)=_ZERO_
    ZW4(JL,3)=_ZERO_
    ZW4(JL,4)=_ZERO_
  
    ZO(JL,1)=_ZERO_
    ZO(JL,2)=_ZERO_
    PFD(JL,KLEV+1)=((_ONE_-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)& ! TOA flux
      &+ PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)
    PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)        ! TOA flux CS
  ENDDO
  
  ! Quentin
  DO JL = KIDIA,KFDIA
     ZTA1(JL)=_ZERO_
     ZTO1(JL)=_ZERO_
  ENDDO
  ! Quentin
  
  DO JK = 1 , KLEV
    IKL = KLEV+1-JK
    DO JL = KIDIA,KFDIA
      ZW4(JL,1)=ZW4(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW4(JL,2)=ZW4(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
      ZW4(JL,3)=ZW4(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW4(JL,4)=ZW4(JL,4)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
    
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
    ENDDO
 
    ! transmission fucntion for all absorbers
    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 4 &
      &, IIND4 &
      &, ZW4  &
      &, ZR4  &
      & )
      ! ZR4 transmission fucntion

    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, 2 &
      &, ZO  &
      &, ZT  &
      & )
      ! ZT transmission function

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR4(JL,1)*ZR4(JL,2)*ZT(JL,1)*ZRJ(JL,JAJ,IKL) ! multiplication of absorber contributions for clouds
      ZDIRF(JL) = ZR4(JL,3)*ZR4(JL,4)*ZT(JL,2)*ZRJ0(JL,JAJ,IKL) ! flux in clear sky part
    ! PDIFF(JL,IKL) = ZDIFF(JL) * RSUN(KNU)*(_ONE_-PCLEAR(JL))
    ! PDIRF(JL,IKL) = ZDIRF(JL) * RSUN(KNU)*PCLEAR(JL)
      PFD(JL,IKL) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&  ! total downward flux
        &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)             ! total downward clear-sky
      
    ! Quentin
      ZTA1(JL) = ZTA1(JL) + ZTAUAZ(JL,IKL)            ! aerosol + rayleigh OD
      ZTO1(JL) = PTAU(JL,KNU,IKL)*(1.-(POMEGA(JL,KNU,IKL)* &  ! cloud OD
      &           PCG(JL,KNU,IKL)*PCG(JL,KNU,IKL))) + ZTO1(JL)
      ZCLDIR = ZDIRF(JL)/ZRJ0(JL,JAJ,1)*EXP(-ZTA1(JL)/PRMU(JL))   ! remaining direct in clear-sky (otherwise diffuse)
      PDIRF(JL,IKL) = ((_ONE_-PCLEAR(JL))*ZCLDIR*EXP(-ZTO1(JL)/PRMU(JL))+& ! some direct through cloud
      &              PCLEAR(JL)*ZCLDIR) * RSUN(KNU)
      PDIRF(JL,IKL) = MIN(PFD(JL,IKL),PDIRF(JL,IKL))
      PDIFF(JL,IKL) = PFD(JL,IKL) - PDIRF(JL,IKL)
    ! Quentin
    
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    ZDIFT(JL) = ZR4(JL,1)*ZR4(JL,2)*ZT(JL,1)*ZTRCLD(JL)     ! true components with corrected cloudiness
    ZDIRT(JL) = ZR4(JL,3)*ZR4(JL,4)*ZT(JL,2)*ZTRCLR(JL)
    PSUDU1(JL) = ((_ONE_-PCLEAR(JL)) * ZDIFT(JL)&           ! not used by ECMWF_VERSION_2
      &+PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)
  ENDDO

!*         3.2.2  UPWARD FLUXES
!                 -------------


  DO JL = KIDIA,KFDIA
    PFU(JL,1) = ((_ONE_-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)&
      &+ PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))&
      &* RSUN(KNU)
    PCU(JL,1) = ZDIRF(JL) * PALBP(JL,KNU) * RSUN(KNU)
  ENDDO

  DO JK = 2 , KLEV+1
    IKM1=JK-1
    DO JL = KIDIA,KFDIA
      ZW4(JL,1)=ZW4(JL,1)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW4(JL,2)=ZW4(JL,2)+PUD(JL,2,IKM1)*1.66_JPRB
      ZW4(JL,3)=ZW4(JL,3)+PUD(JL,1,IKM1)*1.66_JPRB
      ZW4(JL,4)=ZW4(JL,4)+PUD(JL,2,IKM1)*1.66_JPRB
      
      ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKM1)*1.66_JPRB
      ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKM1)*1.66_JPRB
    ENDDO

    CALL SWTT1 ( KIDIA, KFDIA, KLON, KNU, 4 &
      &, IIND4 &
      &, ZW4  &
      &, ZR4  &
      & )

    CALL SWUVO3 ( KIDIA, KFDIA, KLON, KNU, 2 &
      &, ZO  &
      &, ZT  &
      & )

    DO JL = KIDIA,KFDIA
      ZDIFF(JL) = ZR4(JL,1)*ZR4(JL,2)*ZT(JL,1)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR4(JL,3)*ZR4(JL,4)*ZT(JL,2)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((_ONE_-PCLEAR(JL)) * ZDIFF(JL)&
        &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
      PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
    ENDDO
  ENDDO
  
END IF  

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SW1S



