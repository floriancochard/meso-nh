!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2d 2006/02/03 10:00:52
!-----------------------------------------------------------------
SUBROUTINE RRTM_RTRN1A_140GP (KLEV,ISTART,IEND,ICLDLYR,CLDFRAC,TAUCLD,ABSS1 &
  &, OD,TAUSF1,CLFNET,CLHTR,FNET,HTR,TOTDFLUC,TOTDFLUX,TOTUFLUC,TOTUFLUX &
  &, TAVEL,PZ,TZ,TBOUND,PFRAC,SEMISS,SEMISLW,IREFLECT)

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714
!     Speed-up by D.Salmond, ECMWF, 9907
!     Bug-fix by M.J. Iacono, AER, Inc., 9911
!     Bug-fix by JJMorcrette, ECMWF, 991209 (RAT1, RAT2 initialization)
!     Speed-up by D. Salmond, ECMWF, 9912
!     Bug-fix by JJMorcrette, ECMWF, 0005 (extrapolation T<160K)
!     Speed-up by D. Salmond, ECMWF, 000515

!-* This program calculates the upward fluxes, downward fluxes,
!   and heating rates for an arbitrary atmosphere.  The input to
!   this program is the atmospheric profile and all Planck function
!   information.  First-order "numerical" quadrature is used for the 
!   angle integration, i.e. only one exponential is computed per layer
!   per g-value per band.  Cloud overlap is treated with a generalized
!   maximum/random method in which adjacent cloud layers are treated
!   with maximum overlap, and non-adjacent cloud groups are treated
!   with random overlap.  For adjacent cloud layers, cloud information
!   is carried from the previous two layers.


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPGPT   ,JPLAY
USE OYOERRTAB , ONLY : BPADE
USE OYOERRTWN , ONLY : TOTPLNK  ,DELWAVE
USE OYOERRTFTR, ONLY : NGB

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV
INTEGER_M :: ISTART
INTEGER_M :: IEND

INTEGER_M :: ICLDLYR(JPLAY)      ! Cloud indicator
REAL_B :: CLDFRAC(JPLAY)         ! Cloud fraction
REAL_B :: TAUCLD(JPLAY,JPBAND)   ! Spectral optical thickness
REAL_B :: ABSS1 (JPGPT*JPLAY)
REAL_B :: OD    (JPGPT,JPLAY)
REAL_B :: TAUSF1(JPGPT*JPLAY)
REAL_B :: CLFNET  (0:JPLAY)
REAL_B :: CLHTR   (0:JPLAY)
REAL_B :: FNET    (0:JPLAY)
REAL_B :: HTR     (0:JPLAY)
REAL_B :: TOTDFLUC(0:JPLAY)
REAL_B :: TOTDFLUX(0:JPLAY)
REAL_B :: TOTUFLUC(0:JPLAY)
REAL_B :: TOTUFLUX(0:JPLAY)

!- from PROFILE             
REAL_B :: TAVEL(JPLAY)
REAL_B :: PZ(0:JPLAY)
REAL_B :: TZ(0:JPLAY)
REAL_B :: TBOUND

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

!- from SURFACE             
REAL_B :: SEMISS(JPBAND)
REAL_B :: SEMISLW
INTEGER_M :: IREFLECT

INTEGER_M :: INDLAY(JPLAY),INDLEV(0:JPLAY)

REAL_B :: BBU1(JPGPT*JPLAY),BBUTOT1(JPGPT*JPLAY)
REAL_B :: TLAYFRAC(JPLAY),TLEVFRAC(0:JPLAY)
REAL_B :: BGLEV(JPGPT)
!-- DS_000515
REAL_B :: PLVL(0:JPLAY,JPBAND+1),PLAY(0:JPLAY,JPBAND+1),WTNUM(3)
!-- DS_000515
REAL_B :: ODCLD(JPBAND,JPLAY),EFCLFR1(JPBAND,JPLAY)
REAL_B :: ODCLDNW(JPGPT,JPLAY)
REAL_B :: SEMIS(JPGPT),RADUEMIT(JPGPT)

REAL_B :: RADCLRU1(JPGPT) ,RADCLRD1(JPGPT)
REAL_B :: RADLU1(JPGPT)   ,RADLD1(JPGPT)
REAL_B :: ABSCLD1(JPBAND,JPLAY)
!-- DS_000515
REAL_B :: TRNCLD(JPLAY,JPBAND+1)
!-- DS_000515
REAL_B :: ABSCLDNW(JPGPT,JPLAY)
REAL_B :: ATOT1(JPGPT*JPLAY)

REAL_B :: SURFEMIS(JPBAND),PLNKEMIT(JPBAND)

! dimension of arrays required for cloud overlap calculations

REAL_B :: clrradu(jpgpt),cldradu(jpgpt),oldcld(jpgpt)
REAL_B :: oldclr(jpgpt),rad(jpgpt),faccld1(jplay+1),faccld2(jplay+1)
REAL_B :: facclr1(jplay+1),facclr2(jplay+1)
REAL_B :: faccmb1(jplay+1),faccmb2(jplay+1)
REAL_B :: faccld1d(0:jplay),faccld2d(0:jplay),facclr1d(0:jplay)
REAL_B :: facclr2d(0:jplay),faccmb1d(0:jplay),faccmb2d(0:jplay)
REAL_B :: clrradd(jpgpt),cldradd(jpgpt)
INTEGER_M :: istcld(jplay+1),istcldd(0:jplay)
!******

REAL_B :: ZPLVL(JPGPT+1,JPLAY)  ,ZPLAY(JPGPT+1,JPLAY)
REAL_B :: ZTRNCLD(JPGPT+1,JPLAY),ZTAUCLD(JPGPT+1,JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IBAND, ICLDDN, IENT, INDBOUND, INDEX, IPR, LAY, LEV, NBI

!     LOCAL REAL SCALARS
REAL_B :: BBD, BBDTOT, BGLAY, CLDSRC, DBDTLAY, DBDTLEV,&
          &DELBGDN, DELBGUP, DRAD1, DRADCL1, FACTOT1, &
          &FMAX, FMIN, GASSRC, ODSM, PLANKBND, RADCLD, &
          &RADCLU, RADD, RADMOD, RADU, RAT1, RAT2, SUMPL, &
          &SUMPLEM, TBNDFRAC, TRNS, TTOT, URAD1, URADCL1

!--------------------------------------------------------------------------
! Input
!  JPLAY                 ! Maximum number of model layers
!  JPGPT                 ! Total number of g-point subintervals
!  JPBAND                ! Number of longwave spectral bands
!  SECANG                ! Diffusivity angle
!  WTNUM                 ! Weight for radiance to flux conversion
!  KLEV                  ! Number of model layers
!  PAVEL(JPLAY)          ! Mid-layer pressures (hPa)
!  PZ(0:JPLAY)           ! Interface pressures (hPa)
!  TAVEL(JPLAY)          ! Mid-layer temperatures (K)
!  TZ(0:JPLAY)           ! Interface temperatures (K)
!  TBOUND                ! Surface temperature
!  CLDFRAC(JPLAY)        ! Layer cloud fraction
!  TAUCLD(JPLAY,JPBAND)  ! Layer cloud optical thickness
!  ITR
!  PFRAC(JPGPT,JPLAY)    ! Planck function fractions
!  ICLDLYR(JPLAY)        ! Flag for cloudy layers
!  ICLD                  ! Flag for cloudy column
!  IREFLECT              ! Flag for specular reflection
!  SEMISS(JPBAND)        ! Surface spectral emissivity
!  BPADE                 ! Pade constant
!  OD                    ! Clear-sky optical thickness
!  TAUSF1                ! 
!  ABSS1                 !  
!
! Local
!  ABSS(JPGPT*JPLAY)     !
!  ABSCLD(JPLAY)         !
!  ATOT(JPGPT*JPLAY)     !
!  ODCLR(JPGPT,JPLAY)    ! 
!  ODCLD(JPBAND,JPLAY)   !
!  EFCLFR1(JPBAND,JPLAY) ! Effective cloud fraction
!  RADLU(JPGPT)          ! Upward radiance
!  URAD                  ! Spectrally summed upward radiance
!  RADCLRU(JPGPT)        ! Clear-sky upward radiance
!  CLRURAD               ! Spectrally summed clear-sky upward radiance
!  RADLD(JPGPT)          ! Downward radiance
!  DRAD                  ! Spectrally summed downward radiance
!  RADCLRD(JPGPT)        ! Clear-sky downward radiance
!  CLRDRAD               ! Spectrally summed clear-sky downward radiance
!
! Output
!  TOTUFLUX(0:JPLAY)     ! Upward longwave flux
!  TOTDFLUX(0:JPLAY)     ! Downward longwave flux
!  TOTUFLUC(0:JPLAY)     ! Clear-sky upward longwave flux
!  TOTDFLUC(0:JPLAY)     ! Clear-sky downward longwave flux
!
! Maximum/Random cloud overlap variables
! for upward radiaitve transfer
!  FACCLR2  fraction of clear radiance from previous layer that needs to 
!           be switched to cloudy stream
!  FACCLR1  fraction of the radiance that had been switched in the previous
!           layer from cloudy to clear that needs to be switched back to
!           cloudy in the current layer
!  FACCLD2  fraction of cloudy radiance from previous layer that needs to 
!           be switched to clear stream
!           be switched to cloudy stream
!  FACCLD1  fraction of the radiance that had been switched in the previous
!           layer from clear to cloudy that needs to be switched back to
!           clear in the current layer
! for downward radiaitve transfer
!  FACCLR2D fraction of clear radiance from previous layer that needs to 
!           be switched to cloudy stream
!  FACCLR1D fraction of the radiance that had been switched in the previous
!           layer from cloudy to clear that needs to be switched back to
!           cloudy in the current layer
!  FACCLD2D fraction of cloudy radiance from previous layer that needs to 
!           be switched to clear stream
!           be switched to cloudy stream
!  FACCLD1D fraction of the radiance that had been switched in the previous
!           layer from clear to cloudy that needs to be switched back to
!           clear in the current layer
!
!--------------------------------------------------------------------------
!JUAN
  TRNCLD  = _ZERO_
!JUAN
WTNUM(1)=_HALF_
WTNUM(2)=_ZERO_
WTNUM(3)=_ZERO_

!-start JJM_000511
IF (TBOUND < 339._JPRB .AND. TBOUND >= 160._JPRB ) THEN
  INDBOUND = TBOUND - 159._JPRB
  TBNDFRAC = TBOUND - INT(TBOUND)
ELSE IF (TBOUND >= 339._JPRB ) THEN
  INDBOUND = 180
  TBNDFRAC = TBOUND - 339._JPRB
ELSE IF (TBOUND < 160._JPRB ) THEN
  INDBOUND = 1
  TBNDFRAC = TBOUND - 160._JPRB
ENDIF  
!-end JJM_000511
  
DO LAY = 0, KLEV
  TOTUFLUC(LAY) = _ZERO_
  TOTDFLUC(LAY) = _ZERO_
  TOTUFLUX(LAY) = _ZERO_
  TOTDFLUX(LAY) = _ZERO_
!-start JJM_000511
  IF (TZ(LAY) < 339._JPRB .AND. TZ(LAY) >= 160._JPRB ) THEN
    INDLEV(LAY) = TZ(LAY) - 159._JPRB
    TLEVFRAC(LAY) = TZ(LAY) - INT(TZ(LAY))
  ELSE IF (TZ(LAY) >= 339._JPRB ) THEN
    INDLEV(LAY) = 180
    TLEVFRAC(LAY) = TZ(LAY) - 339._JPRB
  ELSE IF (TZ(LAY) < 160._JPRB ) THEN
    INDLEV(LAY) = 1
    TLEVFRAC(LAY) = TZ(LAY) - 160._JPRB
  ENDIF    
!-end JJM_000511
ENDDO

!_start_jjm 991209
DO LEV=0,KLEV
  FACCLD1(LEV+1) = _ZERO_
  FACCLD2(LEV+1) = _ZERO_
  FACCLR1(LEV+1) = _ZERO_
  FACCLR2(LEV+1) = _ZERO_
  FACCMB1(LEV+1) = _ZERO_
  FACCMB2(LEV+1) = _ZERO_
  FACCLD1D(LEV) = _ZERO_
  FACCLD2D(LEV) = _ZERO_
  FACCLR1D(LEV) = _ZERO_
  FACCLR2D(LEV) = _ZERO_
  FACCMB1D(LEV) = _ZERO_
  FACCMB2D(LEV) = _ZERO_
END DO  
RAT1 = _ZERO_
RAT2 = _ZERO_
!_end_jjm 991209



SUMPL   = _ZERO_
SUMPLEM = _ZERO_

ISTCLD(1) = 1
ISTCLDD(KLEV) = 1

DO LEV = 1, KLEV
!-- DS_000515
!-start JJM_000511
  IF (TAVEL(LEV) < 339._JPRB .AND. TAVEL(LEV) >= 160._JPRB ) THEN
    INDLAY(LEV) = TAVEL(LEV) - 159._JPRB
    TLAYFRAC(LEV) = TAVEL(LEV) - INT(TAVEL(LEV))
  ELSE IF (TAVEL(LEV) >= 339._JPRB ) THEN
    INDLAY(LEV) = 180
    TLAYFRAC(LEV) = TAVEL(LEV) - 339._JPRB
  ELSE IF (TAVEL(LEV) < 160._JPRB ) THEN
    INDLAY(LEV) = 1
    TLAYFRAC(LEV) = TAVEL(LEV) - 160._JPRB
  ENDIF  
!-end JJM_000511
END DO
!-- DS_000515

!-- DS_000515
!OCL SCALAR

DO LEV = 1, KLEV
  IF (ICLDLYR(LEV) == 1) THEN

!mji    
    ISTCLD(LEV+1) = 0
    IF (LEV  ==  KLEV) THEN
      FACCLD1(LEV+1) = _ZERO_
      FACCLD2(LEV+1) = _ZERO_
      FACCLR1(LEV+1) = _ZERO_
      FACCLR2(LEV+1) = _ZERO_
!-- DS_000515      
!      FACCMB1(LEV+1) = _ZERO_
!      FACCMB2(LEV+1) = _ZERO_
!mji      ISTCLD(LEV+1) = _ZERO_
    ELSEIF (CLDFRAC(LEV+1)  >=  CLDFRAC(LEV)) THEN
      FACCLD1(LEV+1) = _ZERO_
      FACCLD2(LEV+1) = _ZERO_
      IF (ISTCLD(LEV)  ==  1) THEN
!mji        ISTCLD(LEV+1) = 0
        FACCLR1(LEV+1) = _ZERO_
!mji        
        FACCLR2(LEV+1) = _ZERO_
        IF (CLDFRAC(LEV) < _ONE_) THEN
        FACCLR2(LEV+1) = (CLDFRAC(LEV+1)-CLDFRAC(LEV))/&
         &(_ONE_-CLDFRAC(LEV))
        END IF   
      ELSE
        FMAX = MAX(CLDFRAC(LEV),CLDFRAC(LEV-1))
!mji
        IF (CLDFRAC(LEV+1)  >  FMAX) THEN
          FACCLR1(LEV+1) = RAT2
          FACCLR2(LEV+1) = (CLDFRAC(LEV+1)-FMAX)/(_ONE_-FMAX)
!mji          
        ELSE IF (CLDFRAC(LEV+1) < FMAX) THEN
          FACCLR1(LEV+1) = (CLDFRAC(LEV+1)-CLDFRAC(LEV))/&
           &(CLDFRAC(LEV-1)-CLDFRAC(LEV))
          FACCLR2(LEV+1) = _ZERO_
!mji
        ELSE
          FACCLR1(LEV+1) = RAT2  
          FACCLR2(LEV+1) = _ZERO_
        ENDIF
      ENDIF
      IF (FACCLR1(LEV+1) > _ZERO_ .OR. FACCLR2(LEV+1) > _ZERO_) THEN
        RAT1 = _ONE_
        RAT2 = _ZERO_
      ENDIF
    ELSE
      FACCLR1(LEV+1) = _ZERO_
      FACCLR2(LEV+1) = _ZERO_
      IF (ISTCLD(LEV)  ==  1) THEN
!mji        ISTCLD(LEV+1) = 0
        FACCLD1(LEV+1) = _ZERO_
        FACCLD2(LEV+1) = (CLDFRAC(LEV)-CLDFRAC(LEV+1))/CLDFRAC(LEV)
      ELSE
        FMIN = MIN(CLDFRAC(LEV),CLDFRAC(LEV-1))
        IF (CLDFRAC(LEV+1)  <=  FMIN) THEN
          FACCLD1(LEV+1) = RAT1
          FACCLD2(LEV+1) = (FMIN-CLDFRAC(LEV+1))/FMIN
        ELSE
          FACCLD1(LEV+1) = (CLDFRAC(LEV)-CLDFRAC(LEV+1))/&
           &(CLDFRAC(LEV)-FMIN)
          FACCLD2(LEV+1) = _ZERO_
        ENDIF
      ENDIF
      IF (FACCLD1(LEV+1) > _ZERO_ .OR. FACCLD2(LEV+1) > _ZERO_) THEN
        RAT1 = _ZERO_
        RAT2 = _ONE_
      ENDIF
    ENDIF
!fcc
    IF (LEV == 1) THEN
      FACCMB1(LEV+1) = 0.
      FACCMB2(LEV+1) = FACCLD1(LEV+1) * FACCLR2(LEV)
    ELSE
      FACCMB1(LEV+1) = FACCLR1(LEV+1) * FACCLD2(LEV) *CLDFRAC(LEV-1)
      FACCMB2(LEV+1) = FACCLD1(LEV+1) * FACCLR2(LEV) *&
       &(_ONE_ - CLDFRAC(LEV-1)) 
    ENDIF
!end fcc
  ELSE
!-- DS_000515
    ISTCLD(LEV+1) = 1
  ENDIF
ENDDO

!_start_jjm 991209
RAT1 = _ZERO_
RAT2 = _ZERO_
!_end_jjm 991209

!-- DS_000515
!OCL SCALAR

DO LEV = KLEV, 1, -1
  IF (ICLDLYR(LEV) == 1) THEN
!mji
    ISTCLDD(LEV-1) = 0  
    IF (LEV  ==  1) THEN
      FACCLD1D(LEV-1) = _ZERO_
      FACCLD2D(LEV-1) = _ZERO_
      FACCLR1D(LEV-1) = _ZERO_
      FACCLR2D(LEV-1) = _ZERO_
      FACCMB1D(LEV-1) = _ZERO_
      FACCMB2D(LEV-1) = _ZERO_
!mji      ISTCLDD(LEV-1) = _ZERO_
    ELSEIF (CLDFRAC(LEV-1)  >=  CLDFRAC(LEV)) THEN
      FACCLD1D(LEV-1) = _ZERO_
      FACCLD2D(LEV-1) = _ZERO_
      IF (ISTCLDD(LEV)  ==  1) THEN
!mji        ISTCLDD(LEV-1) = 0
        FACCLR1D(LEV-1) = _ZERO_
        FACCLR2D(LEV-1) = _ZERO_
        IF (CLDFRAC(LEV) < _ONE_) THEN
          FACCLR2D(LEV-1) = (CLDFRAC(LEV-1)-CLDFRAC(LEV))/&
           &(_ONE_-CLDFRAC(LEV))
        END IF
      ELSE
        FMAX = MAX(CLDFRAC(LEV),CLDFRAC(LEV+1))
!mji
        IF (CLDFRAC(LEV-1)  >  FMAX) THEN
          FACCLR1D(LEV-1) = RAT2
          FACCLR2D(LEV-1) = (CLDFRAC(LEV-1)-FMAX)/(_ONE_-FMAX)
!mji
        ELSE IF (CLDFRAC(LEV-1) < FMAX) THEN
          FACCLR1D(LEV-1) = (CLDFRAC(LEV-1)-CLDFRAC(LEV))/&
           &(CLDFRAC(LEV+1)-CLDFRAC(LEV))
          FACCLR2D(LEV-1) = _ZERO_
!mji
        ELSE          
          FACCLR1D(LEV-1) = RAT2
          FACCLR2D(LEV-1) = _ZERO_
        ENDIF
      ENDIF
      IF (FACCLR1D(LEV-1) > _ZERO_ .OR. FACCLR2D(LEV-1) > _ZERO_)THEN
        RAT1 = _ONE_
        RAT2 = _ZERO_
      ENDIF
    ELSE
      FACCLR1D(LEV-1) = _ZERO_
      FACCLR2D(LEV-1) = _ZERO_
      IF (ISTCLDD(LEV)  ==  1) THEN
!mji        ISTCLDD(LEV-1) = 0
        FACCLD1D(LEV-1) = _ZERO_
        FACCLD2D(LEV-1) = (CLDFRAC(LEV)-CLDFRAC(LEV-1))/CLDFRAC(LEV)
      ELSE
        FMIN = MIN(CLDFRAC(LEV),CLDFRAC(LEV+1))
        IF (CLDFRAC(LEV-1)  <=  FMIN) THEN
          FACCLD1D(LEV-1) = RAT1
          FACCLD2D(LEV-1) = (FMIN-CLDFRAC(LEV-1))/FMIN
        ELSE
          FACCLD1D(LEV-1) = (CLDFRAC(LEV)-CLDFRAC(LEV-1))/&
           &(CLDFRAC(LEV)-FMIN)
          FACCLD2D(LEV-1) = _ZERO_
        ENDIF
      ENDIF
      IF (FACCLD1D(LEV-1) > _ZERO_ .OR. FACCLD2D(LEV-1) > _ZERO_)THEN
        RAT1 = _ZERO_
        RAT2 = _ONE_
      ENDIF
    ENDIF
    FACCMB1D(LEV-1) = FACCLR1D(LEV-1) * FACCLD2D(LEV) *CLDFRAC(LEV+1)
    FACCMB2D(LEV-1) = FACCLD1D(LEV-1) * FACCLR2D(LEV) *&
     &(_ONE_ - CLDFRAC(LEV+1))
  ELSE
    ISTCLDD(LEV-1) = 1
  ENDIF
ENDDO


!- Loop over frequency bands.

DO IBAND = ISTART, IEND
  DBDTLEV = TOTPLNK(INDBOUND+1,IBAND)-TOTPLNK(INDBOUND,IBAND)
  PLANKBND = DELWAVE(IBAND) * (TOTPLNK(INDBOUND,IBAND) + TBNDFRAC * DBDTLEV)
  DBDTLEV = TOTPLNK(INDLEV(0)+1,IBAND) -TOTPLNK(INDLEV(0),IBAND)
!-- DS_000515
  PLVL(0,IBAND) = DELWAVE(IBAND)&
   &* (TOTPLNK(INDLEV(0),IBAND) + TLEVFRAC(0)*DBDTLEV)

  SURFEMIS(IBAND) = SEMISS(IBAND)
  PLNKEMIT(IBAND) = SURFEMIS(IBAND) * PLANKBND
  SUMPLEM  = SUMPLEM + PLNKEMIT(IBAND)
  SUMPL    = SUMPL   + PLANKBND
!--DS
ENDDO
!---

!-- DS_000515
DO IBAND = ISTART, IEND
  DO LEV = 1, KLEV
!----              
!- Calculate the integrated Planck functions for at the
!  level and layer temperatures.
!  Compute cloud transmittance for cloudy layers.
    DBDTLEV = TOTPLNK(INDLEV(LEV)+1,IBAND) - TOTPLNK(INDLEV(LEV),IBAND)
    DBDTLAY = TOTPLNK(INDLAY(LEV)+1,IBAND) - TOTPLNK(INDLAY(LEV),IBAND)
!-- DS_000515
    PLAY(LEV,IBAND) = DELWAVE(IBAND)&
     &*(TOTPLNK(INDLAY(LEV),IBAND)+TLAYFRAC(LEV)*DBDTLAY)
    PLVL(LEV,IBAND) = DELWAVE(IBAND)&
     &*(TOTPLNK(INDLEV(LEV),IBAND)+TLEVFRAC(LEV)*DBDTLEV)
    IF (ICLDLYR(LEV) > 0) THEN
      TRNCLD(LEV,IBAND) = EXP(-TAUCLD(LEV,IBAND))
    ENDIF
!-- DS_000515
  ENDDO

ENDDO

SEMISLW = SUMPLEM / SUMPL

!--DS
DO IPR = 1, JPGPT
  NBI = NGB(IPR)
  DO LEV =  1 , KLEV
!-- DS_000515
    ZPLAY(IPR,LEV) = PLAY(LEV,NBI)
    ZPLVL(IPR,LEV) = PLVL(LEV-1,NBI)
    ZTAUCLD(IPR,LEV) = TAUCLD(LEV,NBI)
    ZTRNCLD(IPR,LEV) = TRNCLD(LEV,NBI)
!-- DS_000515
  ENDDO
ENDDO
!----      

!- For cloudy layers, set cloud parameters for radiative transfer.
DO LEV = 1, KLEV
  IF (ICLDLYR(LEV) > 0) THEN
    DO IPR = 1, JPGPT
!--DS          
!            NBI = NGB(IPR)
      ODCLDNW(IPR,LEV) = ZTAUCLD(IPR,LEV)
      ABSCLDNW(IPR,LEV) = _ONE_ - ZTRNCLD(IPR,LEV)
!----            
!            EFCLFRNW(IPR,LEV) = ABSCLDNW(IPR,LEV) * CLDFRAC(LEV)
    ENDDO
  ENDIF
ENDDO

!- Initialize for radiative transfer.
DO IPR = 1, JPGPT
  RADCLRD1(IPR) = _ZERO_
  RADLD1(IPR)   = _ZERO_
  NBI = NGB(IPR)
  SEMIS(IPR) = SURFEMIS(NBI)
  RADUEMIT(IPR) = PFRAC(IPR,1) * PLNKEMIT(NBI)
!-- DS_000515
  BGLEV(IPR) = PFRAC(IPR,KLEV) * PLVL(KLEV,NBI)
ENDDO

!- Downward radiative transfer.
!  *** DRAD1 holds summed radiance for total sky stream
!  *** DRADCL1 holds summed radiance for clear sky stream

ICLDDN = 0
DO LEV = KLEV, 1, -1
  DRAD1   = _ZERO_
  DRADCL1 = _ZERO_

  IF (ICLDLYR(LEV) == 1) THEN

!  *** Cloudy layer
    ICLDDN = 1
    IENT = JPGPT * (LEV-1)
    DO IPR = 1, JPGPT
      INDEX = IENT + IPR
!--DS            
!            NBI = NGB(IPR)
      BGLAY = PFRAC(IPR,LEV) * ZPLAY(IPR,LEV)
!----            
      DELBGUP     = BGLEV(IPR) - BGLAY
      BBU1(INDEX) = BGLAY + TAUSF1(INDEX) * DELBGUP
!--DS            
      BGLEV(IPR) = PFRAC(IPR,LEV) * ZPLVL(IPR,LEV)
!----            
      DELBGDN = BGLEV(IPR) - BGLAY
      BBD = BGLAY + TAUSF1(INDEX) * DELBGDN
!- total-sky downward flux          
      ODSM = OD(IPR,LEV) + ODCLDNW(IPR,LEV)
      FACTOT1 = ODSM / (BPADE + ODSM)
      BBUTOT1(INDEX) = BGLAY + FACTOT1 * DELBGUP
      ATOT1(INDEX) = ABSS1(INDEX) + ABSCLDNW(IPR,LEV)&
       &- ABSS1(INDEX) * ABSCLDNW(IPR,LEV)
      BBDTOT = BGLAY + FACTOT1 * DELBGDN
      GASSRC = BBD * ABSS1(INDEX)
!***
      IF (ISTCLDD(LEV)  ==  1) THEN
        CLDRADD(IPR) = CLDFRAC(LEV) * RADLD1(IPR)
        CLRRADD(IPR) = RADLD1(IPR) - CLDRADD(IPR)
        OLDCLD(IPR) = CLDRADD(IPR)
        OLDCLR(IPR) = CLRRADD(IPR)
        RAD(IPR) = _ZERO_
      ENDIF
      TTOT = _ONE_ - ATOT1(INDEX)
      CLDSRC = BBDTOT * ATOT1(INDEX)
      
! Separate RT equations for clear and cloudy streams      
      CLDRADD(IPR) = CLDRADD(IPR) * TTOT + CLDFRAC(LEV) * CLDSRC
      CLRRADD(IPR) = CLRRADD(IPR) * (_ONE_-ABSS1(INDEX)) +&
       &(_ONE_ - CLDFRAC(LEV)) * GASSRC

!  Total sky downward radiance
      RADLD1(IPR) = CLDRADD(IPR) + CLRRADD(IPR)
      DRAD1 = DRAD1 + RADLD1(IPR)
      
!  Clear-sky downward radiance          
      RADCLRD1(IPR) = RADCLRD1(IPR)+(BBD-RADCLRD1(IPR))*ABSS1(INDEX)
      DRADCL1 = DRADCL1 + RADCLRD1(IPR)

!* Code to account for maximum/random overlap:
!   Performs RT on the radiance most recently switched between clear and
!   cloudy streams
      RADMOD = RAD(IPR) * (FACCLR1D(LEV-1) * (_ONE_-ABSS1(INDEX)) +&
       &FACCLD1D(LEV-1) *  TTOT) - &
       &FACCMB1D(LEV-1) * GASSRC + &
       &FACCMB2D(LEV-1) * CLDSRC
       
!   Computes what the clear and cloudy streams would have been had no
!   radiance been switched       
      OLDCLD(IPR) = CLDRADD(IPR) - RADMOD
      OLDCLR(IPR) = CLRRADD(IPR) + RADMOD
      
!   Computes the radiance to be switched between clear and cloudy.      
      RAD(IPR) = -RADMOD + FACCLR2D(LEV-1)*OLDCLR(IPR) -&
       &FACCLD2D(LEV-1)*OLDCLD(IPR)
      CLDRADD(IPR) = CLDRADD(IPR) + RAD(IPR)
      CLRRADD(IPR) = CLRRADD(IPR) - RAD(IPR)
!***

    ENDDO

  ELSE

!  *** Clear layer
!  *** DRAD1 holds summed radiance for total sky stream
!  *** DRADCL1 holds summed radiance for clear sky stream

    IENT = JPGPT * (LEV-1)
    IF (ICLDDN == 1) THEN
      DO IPR = 1, JPGPT
        INDEX = IENT + IPR
!--DS         
!           NBI = NGB(IPR)
        BGLAY = PFRAC(IPR,LEV) * ZPLAY(IPR,LEV)
!----            
        DELBGUP     = BGLEV(IPR) - BGLAY
        BBU1(INDEX) = BGLAY + TAUSF1(INDEX) * DELBGUP
!--DS            
        BGLEV(IPR) = PFRAC(IPR,LEV) * ZPLVL(IPR,LEV)
!----                      
        DELBGDN = BGLEV(IPR) - BGLAY
        BBD = BGLAY + TAUSF1(INDEX) * DELBGDN
        
!- total-sky downward radiance
        RADLD1(IPR) = RADLD1(IPR)+(BBD-RADLD1(IPR))*ABSS1(INDEX)
        DRAD1 = DRAD1 + RADLD1(IPR)
        
!- clear-sky downward radiance
!-  Set clear sky stream to total sky stream as long as layers
!-  remain clear.  Streams diverge when a cloud is reached.
        RADCLRD1(IPR) = RADCLRD1(IPR)+(BBD-RADCLRD1(IPR))*ABSS1(INDEX)
        DRADCL1 = DRADCL1 + RADCLRD1(IPR)
      ENDDO
            
    ELSE
        
      DO IPR = 1, JPGPT
        INDEX = IENT + IPR
!--DS         
!           NBI = NGB(IPR)
        BGLAY = PFRAC(IPR,LEV) * ZPLAY(IPR,LEV)
!----            
        DELBGUP     = BGLEV(IPR) - BGLAY
        BBU1(INDEX) = BGLAY + TAUSF1(INDEX) * DELBGUP
!--DS            
        BGLEV(IPR) = PFRAC(IPR,LEV) * ZPLVL(IPR,LEV)
!----                      
        DELBGDN = BGLEV(IPR) - BGLAY
        BBD = BGLAY + TAUSF1(INDEX) * DELBGDN
!- total-sky downward flux          
        RADLD1(IPR) = RADLD1(IPR)+(BBD-RADLD1(IPR))*ABSS1(INDEX)
        DRAD1 = DRAD1 + RADLD1(IPR)
!- clear-sky downward flux          
!-  Set clear sky stream to total sky stream as long as layers
!-  remain clear.  Streams diverge when a cloud is reached.
        RADCLRD1(IPR) = RADLD1(IPR)
      ENDDO
      DRADCL1 = DRAD1
    ENDIF
    
  ENDIF

  TOTDFLUC(LEV-1) = DRADCL1 * WTNUM(1)
  TOTDFLUX(LEV-1) = DRAD1   * WTNUM(1)

ENDDO





! Spectral reflectivity and reflectance
! Includes the contribution of spectrally varying longwave emissivity 
! and reflection from the surface to the upward radiative transfer.
! Note: Spectral and Lambertian reflections are identical for the one
! angle flux integration used here.

URAD1   = _ZERO_
URADCL1 = _ZERO_

!start JJM_000511
!IF (IREFLECT  ==  0) THEN
!- Lambertian reflection.
  DO IPR = 1, JPGPT
! Clear-sky radiance
!    RADCLD = _TWO_ * (RADCLRD1(IPR) * WTNUM(1) )
    RADCLD = RADCLRD1(IPR)
    RADCLRU1(IPR) = RADUEMIT(IPR) + (_ONE_ - SEMIS(IPR)) * RADCLD
    URADCL1 = URADCL1 + RADCLRU1(IPR)

! Total sky radiance
!    RADD = _TWO_ * (RADLD1(IPR) * WTNUM(1) )
    RADD = RADLD1(IPR)
    RADLU1(IPR) = RADUEMIT(IPR) + (_ONE_ - SEMIS(IPR)) * RADD
    URAD1 = URAD1 + RADLU1(IPR)
  ENDDO
  TOTUFLUC(0) = URADCL1 * _HALF_
  TOTUFLUX(0) = URAD1 * _HALF_
!ELSE
!!- Specular reflection.
!  DO IPR = 1, JPGPT
!    RADCLU = RADUEMIT(IPR)
!    RADCLRU1(IPR) = RADCLU + (_ONE_ - SEMIS(IPR)) * RADCLRD1(IPR)
!    URADCL1 = URADCL1 + RADCLRU1(IPR)
!
!    RADU = RADUEMIT(IPR)
!    RADLU1(IPR) = RADU + (_ONE_ - SEMIS(IPR)) * RADLD1(IPR)
!    URAD1 = URAD1 + RADLU1(IPR)
!  ENDDO
!  TOTUFLUC(0) = URADCL1 * WTNUM(1)
!  TOTUFLUX(0) = URAD1   * WTNUM(1)
!ENDIF

!- Upward radiative transfer.
!- *** URAD1 holds the summed radiance for total sky stream
!- *** URADCL1 holds the summed radiance for clear sky stream
DO LEV = 1, KLEV
  URAD1   = _ZERO_
  URADCL1 = _ZERO_

! Check flag for cloud in current layer
  IF (ICLDLYR(LEV) == 1) THEN

!- *** Cloudy layer
    IENT = JPGPT * (LEV-1)
    DO IPR = 1, JPGPT
      INDEX = IENT + IPR
!- total-sky upward flux          
      GASSRC = BBU1(INDEX) * ABSS1(INDEX)
!
!- If first cloudy layer in sequence, split up radiance into clear and
!    cloudy streams depending on cloud fraction
      IF (ISTCLD(LEV)  ==  1) THEN
        CLDRADU(IPR) = CLDFRAC(LEV) * RADLU1(IPR)
        CLRRADU(IPR) = RADLU1(IPR) - CLDRADU(IPR)
        OLDCLD(IPR) = CLDRADU(IPR)
        OLDCLR(IPR) = CLRRADU(IPR)
        RAD(IPR) = _ZERO_
      ENDIF
      TTOT = _ONE_ - ATOT1(INDEX)
      TRNS = _ONE_ - ABSS1(INDEX)
      CLDSRC = BBUTOT1(INDEX) * ATOT1(INDEX)
!      
!- Separate RT equations for clear and cloudy streams      
      CLDRADU(IPR) = CLDRADU(IPR) * TTOT + CLDFRAC(LEV) * CLDSRC
      CLRRADU(IPR) = CLRRADU(IPR) * TRNS +(_ONE_ - CLDFRAC(LEV)) * GASSRC
!***

!- total sky upward flux
      RADLU1(IPR) = CLDRADU(IPR) + CLRRADU(IPR)
      URAD1 = URAD1 + RADLU1(IPR)
      
!- clear-sky upward flux
      RADCLRU1(IPR) = RADCLRU1(IPR) + (BBU1(INDEX)-RADCLRU1(IPR))&
       &*ABSS1(INDEX)
      URADCL1 = URADCL1 + RADCLRU1(IPR)

!* Code to account for maximum/random overlap:
!   Performs RT on the radiance most recently switched between clear and
!   cloudy streams
      RADMOD = RAD(IPR) * (FACCLR1(LEV+1) * TRNS +&
       &FACCLD1(LEV+1) *  TTOT) - &
       &FACCMB1(LEV+1) * GASSRC + &
       &FACCMB2(LEV+1) * CLDSRC
       
!   Computes what the clear and cloudy streams would have been had no
!   radiance been switched       
      OLDCLD(IPR) = CLDRADU(IPR) - RADMOD
      OLDCLR(IPR) = CLRRADU(IPR) + RADMOD
      
!   Computes the radiance to be switched between clear and cloudy.      
      RAD(IPR) = -RADMOD + FACCLR2(LEV+1)*OLDCLR(IPR) -&
       &FACCLD2(LEV+1)*OLDCLD(IPR)
      CLDRADU(IPR) = CLDRADU(IPR) + RAD(IPR)
      CLRRADU(IPR) = CLRRADU(IPR) - RAD(IPR)
!***
    ENDDO

  ELSE

!- *** Clear layer
    IENT = JPGPT * (LEV-1)
    DO IPR = 1, JPGPT
      INDEX = IENT + IPR
!- total-sky upward flux          
      RADLU1(IPR) = RADLU1(IPR)+(BBU1(INDEX)-RADLU1(IPR))*ABSS1(INDEX)
      URAD1 = URAD1 + RADLU1(IPR)
!- clear-sky upward flux
!   Upward clear and total sky streams must be separate because surface
!   reflectance is different for each.
      RADCLRU1(IPR) = RADCLRU1(IPR)+(BBU1(INDEX)-RADCLRU1(IPR))*ABSS1(INDEX)
      URADCL1 = URADCL1 + RADCLRU1(IPR)
    ENDDO

  ENDIF

  TOTUFLUC(LEV) = URADCL1 * WTNUM(1)
  TOTUFLUX(LEV) = URAD1   * WTNUM(1)

ENDDO


!* Convert radiances to fluxes and heating rates for total and clear sky.
! ** NB: moved to calling routine
!      TOTUFLUC(0) = TOTUFLUC(0) * FLUXFAC
!      TOTDFLUC(0) = TOTDFLUC(0) * FLUXFAC
!      TOTUFLUX(0) = TOTUFLUX(0) * FLUXFAC
!      TOTDFLUX(0) = TOTDFLUX(0) * FLUXFAC

!      CLFNET(0) = TOTUFLUC(0) - TOTDFLUC(0)
!      FNET(0)   = TOTUFLUX(0) - TOTDFLUX(0)
!      DO LEV = 1, KLEV
!        TOTUFLUC(LEV) = TOTUFLUC(LEV) * FLUXFAC
!        TOTDFLUC(LEV) = TOTDFLUC(LEV) * FLUXFAC
!        CLFNET(LEV) = TOTUFLUC(LEV) - TOTDFLUC(LEV)

!        TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
!        TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
!        FNET(LEV) = TOTUFLUX(LEV) - TOTDFLUX(LEV)
!        L = LEV - 1

!- Calculate Heating Rates.
!        CLHTR(L)=HEATFAC*(CLFNET(L)-CLFNET(LEV))/(PZ(L)-PZ(LEV)) 
!        HTR(L)  =HEATFAC*(FNET(L)  -FNET(LEV))  /(PZ(L)-PZ(LEV)) 
!      END DO
!      CLHTR(KLEV) = 0.0
!      HTR(KLEV)   = 0.0

RETURN
END SUBROUTINE RRTM_RTRN1A_140GP
