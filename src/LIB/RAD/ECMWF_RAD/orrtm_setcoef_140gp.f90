!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:41
!-----------------------------------------------------------------
SUBROUTINE ORRTM_SETCOEF_140GP (KLEV,COLDRY,WKL &
 &, FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1 &
 &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,COLO2,CO2MULT &
 &, LAYTROP,LAYSWTCH,LAYLOW,PAVEL,TAVEL,SELFFAC,SELFFRAC,INDSELF)

!     Reformatted for F90 by JJMorcrette, ECMWF, 980714

!     Purpose:  For a given atmosphere, calculate the indices and
!     fractions related to the pressure and temperature interpolations.
!     Also calculate the values of the integrated Planck functions 
!     for each band at the level and layer temperatures.

#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY     ,JPBAND    ,JPGPT   ,JPINPX
USE OYOERRTRF , ONLY : PREF      ,PREFLOG   ,TREF

IMPLICIT NONE

REAL_B :: COLDRY(JPLAY)
REAL_B :: WKL(JPINPX,JPLAY)

!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

!- from INTFAC      
REAL_B :: FAC00(JPLAY)
REAL_B :: FAC01(JPLAY)
REAL_B :: FAC10(JPLAY)
REAL_B :: FAC11(JPLAY)
REAL_B :: FORFAC(JPLAY)

!- from INTIND
INTEGER_M :: JP(JPLAY)
INTEGER_M :: JT(JPLAY)
INTEGER_M :: JT1(JPLAY)

!- from PROFDATA             
REAL_B :: COLH2O(JPLAY)
REAL_B :: COLCO2(JPLAY)
REAL_B :: COLO3 (JPLAY)
REAL_B :: COLN2O(JPLAY)
REAL_B :: COLCH4(JPLAY)
REAL_B :: COLO2 (JPLAY)
REAL_B :: CO2MULT(JPLAY)
INTEGER_M :: LAYTROP
INTEGER_M :: LAYSWTCH
INTEGER_M :: LAYLOW

!- from PROFILE             
REAL_B :: PAVEL(JPLAY)
REAL_B :: TAVEL(JPLAY)

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)


!     LOCAL INTEGER SCALARS
INTEGER_M :: JP1, LAY

!     LOCAL REAL SCALARS
REAL_B :: CO2REG, COMPFP, FACTOR, FP, FT, FT1, PLOG, SCALEFAC, STPFAC, WATER


!#include "yoeratm.h"    

STPFAC = 296._JPRB/1013._JPRB

LAYTROP  = 0
LAYSWTCH = 0
LAYLOW   = 0
DO LAY = 1, KLEV
!        Find the two reference pressures on either side of the
!        layer pressure.  Store them in JP and JP1.  Store in FP the
!        fraction of the difference (in ln(pressure)) between these
!        two values that the layer pressure lies.
  PLOG = LOG(PAVEL(LAY))
  JP(LAY) = INT(36._JPRB - 5*(PLOG+0.04_JPRB))
  IF (JP(LAY)  <  1) THEN
    JP(LAY) = 1
  ELSEIF (JP(LAY)  >  58) THEN
    JP(LAY) = 58
  ENDIF
  JP1 = JP(LAY) + 1
  FP = 5._JPRB * (PREFLOG(JP(LAY)) - PLOG)

!        Determine, for each reference pressure (JP and JP1), which
!        reference temperature (these are different for each  
!        reference pressure) is nearest the layer temperature but does
!        not exceed it.  Store these indices in JT and JT1, resp.
!        Store in FT (resp. FT1) the fraction of the way between JT
!        (JT1) and the next highest reference temperature that the 
!        layer temperature falls.
  JT(LAY) = INT(3._JPRB + (TAVEL(LAY)-TREF(JP(LAY)))/15._JPRB)
  IF (JT(LAY)  <  1) THEN
    JT(LAY) = 1
  ELSEIF (JT(LAY)  >  4) THEN
    JT(LAY) = 4
  ENDIF
  FT = ((TAVEL(LAY)-TREF(JP(LAY)))/15._JPRB) - REAL(JT(LAY)-3)
  JT1(LAY) = INT(3._JPRB + (TAVEL(LAY)-TREF(JP1))/15._JPRB)
  IF (JT1(LAY)  <  1) THEN
    JT1(LAY) = 1
  ELSEIF (JT1(LAY)  >  4) THEN
    JT1(LAY) = 4
  ENDIF
  FT1 = ((TAVEL(LAY)-TREF(JP1))/15._JPRB) - REAL(JT1(LAY)-3)

  WATER = WKL(1,LAY)/COLDRY(LAY)
  SCALEFAC = PAVEL(LAY) * STPFAC / TAVEL(LAY)

!        If the pressure is less than ~100mb, perform a different
!        set of species interpolations.
!         IF (PLOG .LE. 4.56) GO TO 5300
!--------------------------------------         
  IF (PLOG  >  4.56_JPRB) THEN
    LAYTROP =  LAYTROP + 1
!        For one band, the "switch" occurs at ~300 mb. 
    IF (PLOG  >=  5.76_JPRB) LAYSWTCH = LAYSWTCH + 1
    IF (PLOG  >=  6.62_JPRB) LAYLOW = LAYLOW + 1

    FORFAC(LAY) = SCALEFAC / (_ONE_+WATER)

!        Set up factors needed to separately include the water vapor
!        self-continuum in the calculation of absorption coefficient.
!C           SELFFAC(LAY) = WATER * SCALEFAC / (1.+WATER)
    SELFFAC(LAY) = WATER * FORFAC(LAY)
    FACTOR = (TAVEL(LAY)-188.0_JPRB)/7.2_JPRB
    INDSELF(LAY) = MIN(9, MAX(1, INT(FACTOR)-7))
    SELFFRAC(LAY) = FACTOR - REAL(INDSELF(LAY) + 7)

!        Calculate needed column amounts.
    COLH2O(LAY) = 1.E-20_JPRB * WKL(1,LAY)
    COLCO2(LAY) = 1.E-20_JPRB * WKL(2,LAY)
    COLO3(LAY)  = 1.E-20_JPRB * WKL(3,LAY)
    COLN2O(LAY) = 1.E-20_JPRB * WKL(4,LAY)
    COLCH4(LAY) = 1.E-20_JPRB * WKL(6,LAY)
    COLO2(LAY)  = 1.E-20_JPRB * WKL(7,LAY)
    IF (COLCO2(LAY)  ==  _ZERO_) COLCO2(LAY) = 1.E-32_JPRB * COLDRY(LAY)
    IF (COLN2O(LAY)  ==  _ZERO_) COLN2O(LAY) = 1.E-32_JPRB * COLDRY(LAY)
    IF (COLCH4(LAY)  ==  _ZERO_) COLCH4(LAY) = 1.E-32_JPRB * COLDRY(LAY)
!        Using E = 1334.2 cm-1.
    CO2REG = 3.55E-24_JPRB * COLDRY(LAY)
    CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) *&
     &272.63_JPRB*EXP(-1919.4_JPRB/TAVEL(LAY))/(8.7604E-4_JPRB*TAVEL(LAY))
!         GO TO 5400
!------------------
  ELSE
!        Above LAYTROP.
! 5300    CONTINUE

!        Calculate needed column amounts.
    FORFAC(LAY) = SCALEFAC / (_ONE_+WATER)

    COLH2O(LAY) = 1.E-20_JPRB * WKL(1,LAY)
    COLCO2(LAY) = 1.E-20_JPRB * WKL(2,LAY)
    COLO3(LAY)  = 1.E-20_JPRB * WKL(3,LAY)
    COLN2O(LAY) = 1.E-20_JPRB * WKL(4,LAY)
    COLCH4(LAY) = 1.E-20_JPRB * WKL(6,LAY)
    COLO2(LAY)  = 1.E-20_JPRB * WKL(7,LAY)
    IF (COLCO2(LAY)  ==  _ZERO_) COLCO2(LAY) = 1.E-32_JPRB * COLDRY(LAY)
    IF (COLN2O(LAY)  ==  _ZERO_) COLN2O(LAY) = 1.E-32_JPRB * COLDRY(LAY)
    IF (COLCH4(LAY)  ==  _ZERO_) COLCH4(LAY) = 1.E-32_JPRB * COLDRY(LAY)
    CO2REG = 3.55E-24_JPRB * COLDRY(LAY)
    CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) *&
     &272.63_JPRB*EXP(-1919.4_JPRB/TAVEL(LAY))/(8.7604E-4_JPRB*TAVEL(LAY))
!----------------     
  ENDIF
! 5400    CONTINUE

!        We have now isolated the layer ln pressure and temperature,
!        between two reference pressures and two reference temperatures 
!        (for each reference pressure).  We multiply the pressure 
!        fraction FP with the appropriate temperature fractions to get 
!        the factors that will be needed for the interpolation that yields
!        the optical depths (performed in routines TAUGBn for band n).

  COMPFP = _ONE_ - FP
  FAC10(LAY) = COMPFP * FT
  FAC00(LAY) = COMPFP * (_ONE_ - FT)
  FAC11(LAY) = FP * FT1
  FAC01(LAY) = FP * (_ONE_ - FT1)

ENDDO

! MT 981104 
!-- Set LAYLOW for profiles with surface pressure less than 750 hPa. 
IF (LAYLOW == 0) LAYLOW=1

RETURN
END SUBROUTINE ORRTM_SETCOEF_140GP
