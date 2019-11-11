!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:41
!-----------------------------------------------------------------
!***************************************************************************
!                                                                          *
!                RRTM :  RAPID RADIATIVE TRANSFER MODEL                    *
!                                                                          *
!             ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                 *
!                        840 MEMORIAL DRIVE                                *
!                        CAMBRIDGE, MA 02139                               *
!                                                                          *
!                           ELI J. MLAWER                                  *
!                         STEVEN J. TAUBMAN~                               *
!                         SHEPARD A. CLOUGH                                *
!                                                                          *
!                        ~currently at GFDL                                *
!                                                                          *
!                       email:  mlawer@aer.com                             *
!                                                                          *
!        The authors wish to acknowledge the contributions of the          *
!        following people:  Patrick D. Brown, Michael J. Iacono,           *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
!                                                                          *
!***************************************************************************
!     Reformatted for F90 by JJMorcrette, ECMWF, 980714                    * 
!                                                                          *
!***************************************************************************
! *** mji ***
! *** This version of RRTM has been altered to interface with either
!     the ECMWF numerical weather prediction model or the ECMWF column 
!     radiation model (ECRT) package. 

!     Revised, April, 1997;  Michael J. Iacono, AER, Inc.
!          - initial implementation of RRTM in ECRT code
!     Revised, June, 1999;  Michael J. Iacono and Eli J. Mlawer, AER, Inc.
!          - to implement generalized maximum/random cloud overlap

SUBROUTINE RRTM_RRTM_140GP &
 &( KIDIA , KFDIA , KLON , KLEV &
 &, PAER  , PAPH  , PAP &
 &, PTS   , PTH   , PT &
 &, ZEMIS , ZEMIW &
 &, PQ    , PCCO2 , POZN &
 &, PCLDF , PTAUCLD &
 &, PEMIT , PFLUX , PFLUC, PTCLEAR &
 &)

! *** This program is the driver for RRTM, the AER rapid model.  
!     For each atmosphere the user wishes to analyze, this routine
!     a) calls ECRTATM to read in the atmospheric profile 
!     b) calls SETCOEF to calculate various quantities needed for 
!        the radiative transfer algorithm
!     c) calls RTRN to do the radiative transfer calculation for
!        clear or cloudy sky
!     d) writes out the upward, downward, and net flux for each
!        level and the heating rate for each layer


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPBAND   ,JPG      ,JPXSEC   ,JPGPT    ,JPLAY    ,&
            &JPINPX
USE OYOERRTWN , ONLY : WAVENUM1 ,WAVENUM2 ,DELWAVE  ,NG       ,NSPA     ,NSPB
!
!USE MODI_ORRTM_ECRT_140GP
!USE MODI_ORRTM_SETCOEF_140GP
!USE MODI_ORRTM_GASABS1A_140GP
!USE MODI_RRTM_RTRN1A_140GP
!
!------------------------------Arguments--------------------------------

! Input arguments


IMPLICIT NONE
INTEGER_M :: kidia                 ! First atmosphere index
INTEGER_M :: kfdia                 ! Last atmosphere index
INTEGER_M :: klon                  ! Number of atmospheres (longitudes)
INTEGER_M :: klev                  ! Number of atmospheric layers
REAL_B :: paer(klon,6,klev)        ! Aerosol optical thickness
REAL_B :: pap(klon,klev)           ! Layer pressures (Pa)
REAL_B :: paph(klon,klev+1)        ! Interface pressures (Pa)
REAL_B :: pts(klon)                ! Surface temperature (K)
REAL_B :: pt(klon,klev)            ! Layer temperature (K)
REAL_B :: zemis(klon)              ! Non-window surface emissivity
REAL_B :: zemiw(klon)              ! Window surface emissivity
REAL_B :: pth(klon,klev+1)         ! Interface temperatures (K)
REAL_B :: pq(klon,klev)            ! H2O specific humidity (mmr)
REAL_B :: pozn(klon,klev)          ! O3 mass mixing ratio
REAL_B :: pcco2                    ! CO2 mass mixing ratio
REAL_B :: rch4                     ! CH4 mass mixing ratio
REAL_B :: rn2o                     ! N2O mass mixing ratio
REAL_B :: rcfc11                   ! CFC11 mass mixing ratio
REAL_B :: rcfc12                   ! CFC12 mass mixing ratio
REAL_B :: pcldf(klon,klev)         ! Cloud fraction
REAL_B :: ptaucld(klon,klev,JPBAND)! Cloud optical depth
REAL_B :: PFLUX(klon,2,klev+1)     ! LW total sky flux (1=up, 2=down)
REAL_B :: PFLUC(klon,2,klev+1)     ! LW clear sky flux (1=up, 2=down)
REAL_B :: PEMIT(klon)              ! Surface LW emissivity
REAL_B :: PTCLEAR(klon)            ! clear-sky fraction of column

INTEGER_M :: ICLDLYR(JPLAY)        ! Cloud indicator
REAL_B :: CLDFRAC(JPLAY)           ! Cloud fraction
REAL_B :: TAUCLD(JPLAY,JPBAND)     ! Spectral optical thickness

REAL_B :: ABSS1 (JPGPT*JPLAY)
REAL_B :: ATR1  (JPGPT,JPLAY)
EQUIVALENCE (ABSS1(1),ATR1(1,1))

REAL_B :: OD    (JPGPT,JPLAY)

REAL_B :: TAUSF1(JPGPT*JPLAY)
REAL_B :: TF1   (JPGPT,JPLAY)
EQUIVALENCE (TAUSF1(1),TF1(1,1))

REAL_B :: COLDRY(JPLAY)
REAL_B :: WKL(JPINPX,JPLAY)

REAL_B :: WX(JPXSEC,JPLAY)         ! Amount of trace gases

REAL_B :: CLFNET  (0:JPLAY)
REAL_B :: CLHTR   (0:JPLAY)
REAL_B :: FNET    (0:JPLAY)
REAL_B :: HTR     (0:JPLAY)
REAL_B :: TOTDFLUC(0:JPLAY)
REAL_B :: TOTDFLUX(0:JPLAY)
REAL_B :: TOTUFLUC(0:JPLAY)
REAL_B :: TOTUFLUX(0:JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: i, icld, iplon, K
INTEGER_M :: ISTART
INTEGER_M :: IEND

!     LOCAL REAL SCALARS
REAL_B :: FLUXFAC, HEATFAC, PI, ZEPSEC, ZTCLEAR

!- from AER
REAL_B :: TAUAERL(JPLAY,JPBAND)

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

!- from PRECISE             
REAL_B :: ONEMINUS

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
REAL_B :: PZ(0:JPLAY)
REAL_B :: TZ(0:JPLAY)
REAL_B :: TBOUND
INTEGER_M :: NLAYERS

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

!- from SURFACE             
REAL_B :: SEMISS(JPBAND)
REAL_B :: SEMISLW
INTEGER_M :: IREFLECT


!     HEATFAC is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)

ZEPSEC = 1.E-06_JPRB
ONEMINUS = _ONE_ - ZEPSEC
PI = _TWO_*ASIN(_ONE_)
FLUXFAC = PI * 2.D4
HEATFAC = 8.4391_JPRB

! *** mji ***
! For use with ECRT, this loop is over atmospheres (or longitudes)
DO iplon = kidia,kfdia

! *** mji ***
!- Prepare atmospheric profile from ECRT for use in RRTM, and define
!  other RRTM input parameters.  Arrays are passed back through the
!  existing RRTM commons and arrays.
  ZTCLEAR=_ONE_

!  print *,'before RRTM_ECRT_140GP'

  CALL ORRTM_ECRT_140GP &
   &( iplon, klon , klev, icld &
   &, paer , paph , pap &
   &, pts  , pth  , pt &
   &, zemis, zemiw &
   &, pq   , pcco2, pozn, pcldf, ptaucld, ztclear &
   &, CLDFRAC,TAUCLD,COLDRY,WKL,WX &
   &, TAUAERL,PAVEL,TAVEL,PZ,TZ,TBOUND,NLAYERS,SEMISS,IREFLECT)

  PTCLEAR(iplon)=ztclear

  ISTART = 1
  IEND   = 16

!  Calculate information needed by the radiative transfer routine
!  that is specific to this atmosphere, especially some of the 
!  coefficients and indices needed to compute the optical depths
!  by interpolating data from stored reference atmospheres. 

!  print *,'before RRTM_SETCOEF_140GP'
  
    CALL ORRTM_SETCOEF_140GP (KLEV,COLDRY,WKL &
   &, FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1 &
   &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,COLO2,CO2MULT &
   &, LAYTROP,LAYSWTCH,LAYLOW,PAVEL,TAVEL,SELFFAC,SELFFRAC,INDSELF)

! print *,'before RRTM_GASABS1A_140GP'

  CALL ORRTM_GASABS1A_140GP (KLEV,ATR1,OD,TF1,COLDRY,WX &
   &, TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS &
   &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,COLO2,CO2MULT &
   &, LAYTROP,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!- Call the radiative transfer routine.

! *** mji ***
!  Check for cloud in column.  Use ECRT threshold set as flag icld in
!  routine ECRTATM.  If icld=1 then column is cloudy, otherwise it is
!  clear.  Also, set up flag array, icldlyr, for use in radiative
!  transfer.  Set icldlyr to one for each layer with non-zero cloud
!  fraction.

  DO K = 1, KLEV
    IF (ICLD == 1.AND.CLDFRAC(K) > ZEPSEC) THEN
      ICLDLYR(K) = 1
    ELSE
      ICLDLYR(K) = 0
    ENDIF
  ENDDO

!  Clear and cloudy parts of column are treated together in RTRN.
!  Clear radiative transfer is done for clear layers and cloudy radiative
!  transfer is done for cloudy layers as identified by icldlyr.
  
! print *,'before RRTM_RTRN1A_140GP'
 
  CALL RRTM_RTRN1A_140GP (KLEV,ISTART,IEND,ICLDLYR,CLDFRAC,TAUCLD,ABSS1 &
   &, OD,TAUSF1,CLFNET,CLHTR,FNET,HTR,TOTDFLUC,TOTDFLUX,TOTUFLUC,TOTUFLUX &
   &, TAVEL,PZ,TZ,TBOUND,PFRAC,SEMISS,SEMISLW,IREFLECT)

! ***   Pass clear sky and total sky up and down flux profiles to ECRT
!       output arrays (zflux, zfluc). Array indexing from bottom to top 
!       is preserved for ECRT.
!       Invert down flux arrays for consistency with ECRT sign conventions.

  pemit(iplon) = SEMISLW
  DO i = 0, KLEV
    pfluc(iplon,1,i+1) =  TOTUFLUC(i)*FLUXFAC
    pfluc(iplon,2,i+1) = -TOTDFLUC(i)*FLUXFAC
    pflux(iplon,1,i+1) =  TOTUFLUX(i)*FLUXFAC
    pflux(iplon,2,i+1) = -TOTDFLUX(i)*FLUXFAC
  ENDDO
ENDDO

RETURN
END SUBROUTINE RRTM_RRTM_140GP
