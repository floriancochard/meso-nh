!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:41
!-----------------------------------------------------------------
!******************************************************************************
!                                                                             *
!                  Optical depths developed for the                           *
!                                                                             *
!                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
!                                                                             *
!            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
!                        840 MEMORIAL DRIVE                                   *
!                        CAMBRIDGE, MA 02139                                  *
!                                                                             *
!                           ELI J. MLAWER                                     *
!                         STEVEN J. TAUBMAN                                   *
!                         SHEPARD A. CLOUGH                                   *
!                                                                             *
!                       email:  mlawer@aer.com                                *
!                                                                             *
!        The authors wish to acknowledge the contributions of the             *
!        following people:  Patrick D. Brown, Michael J. Iacono,              *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
!                                                                             *
!******************************************************************************
! Modified by:                                                                *
!      JJ Morcrette 980714 ECMWF      for use on ECMWF's Fujitsu VPP770       *
!         Reformatted for F90 by JJMorcrette, ECMWF                           * 
!         - replacing COMMONs by MODULEs                                      *
!         - changing labelled to unlabelled DO loops                          *
!         - creating set-up routines for all block data statements            *
!         - reorganizing the parameter statements                             * 
!         - passing KLEV as argument                                          *
!         - suppressing some equivalencing                                    *
!                                                                             *
!      D Salmond    9907   ECMWF      Speed-up modifications                  *
!      D Salmond    000515 ECMWF      Speed-up modifications                  *
!******************************************************************************
!     TAUMOL                                                                  *
!                                                                             *
!     This file contains the subroutines TAUGBn (where n goes from            *
!     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
!     per g-value and layer for band n.                                       *
!                                                                             *
!  Output:  optical depths (unitless)                                         *
!           fractions needed to compute Planck functions at every layer       *
!               and g-value                                                   *
!                                                                             *
!     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
!     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
!                                                                             *
!  Input                                                                      *
!                                                                             *
!     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
!     COMMON /PRECISE/  ONEMINUS                                              *
!     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
!    &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND                        *
!     COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,                              *
!    &                  COLH2O(MXLAY),COLCO2(MXLAY),                          *
!    &                  COLO3(MXLAY),COLN2O(MXLAY),COLCH4(MXLAY),             *
!    &                  COLO2(MXLAY),CO2MULT(MXLAY)                           *
!     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
!    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
!     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
!     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
!                                                                             *
!     Description:                                                            *
!     NG(IBAND) - number of g-values in band IBAND                            *
!     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
!                   atmospheres that are stored for band IBAND per            *
!                   pressure level and temperature.  Each of these            *
!                   atmospheres has different relative amounts of the         *
!                   key species for the band (i.e. different binary           *
!                   species parameters).                                      *
!     NSPB(IBAND) - same for upper atmosphere                                 *
!     ONEMINUS - since problems are caused in some cases by interpolation     *
!                parameters equal to or greater than 1, for these cases       *
!                these parameters are set to this value, slightly < 1.        *
!     PAVEL - layer pressures (mb)                                            *
!     TAVEL - layer temperatures (degrees K)                                  *
!     PZ - level pressures (mb)                                               *
!     TZ - level temperatures (degrees K)                                     *
!     LAYTROP - layer at which switch is made from one combination of         *
!               key species to another                                        *
!     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
!               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
!               respectively (molecules/cm**2)                                *
!     CO2MULT - for bands in which carbon dioxide is implemented as a         *
!               trace species, this is the factor used to multiply the        *
!               band's average CO2 absorption coefficient to get the added    *
!               contribution to the optical depth relative to 355 ppm.        *
!     FACij(LAY) - for layer LAY, these are factors that are needed to        *
!                  compute the interpolation factors that multiply the        *
!                  appropriate reference k-values.  A value of 0 (1) for      *
!                  i,j indicates that the corresponding factor multiplies     *
!                  reference k-value for the lower (higher) of the two        *
!                  appropriate temperatures, and altitudes, respectively.     *
!     JP - the index of the lower (in altitude) of the two appropriate        *
!          reference pressure levels needed for interpolation                 *
!     JT, JT1 - the indices of the lower of the two appropriate reference     *
!               temperatures needed for interpolation (for pressure           *
!               levels JP and JP+1, respectively)                             *
!     SELFFAC - scale factor needed to water vapor self-continuum, equals     *
!               (water vapor density)/(atmospheric density at 296K and        *
!               1013 mb)                                                      *
!     SELFFRAC - factor needed for temperature interpolation of reference     *
!                water vapor self-continuum data                              *
!     INDSELF - index of the lower of the two appropriate reference           *
!               temperatures needed for the self-continuum interpolation      *
!                                                                             *
!  Data input                                                                 *
!     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
!        (note:  n is the band number)                                        *
!                                                                             *
!     Description:                                                            *
!     KA - k-values for low reference atmospheres (no water vapor             *
!          self-continuum) (units: cm**2/molecule)                            *
!     KB - k-values for high reference atmospheres (all sources)              *
!          (units: cm**2/molecule)                                            *
!     SELFREF - k-values for water vapor self-continuum for reference         *
!               atmospheres (used below LAYTROP)                              *
!               (units: cm**2/molecule)                                       *
!                                                                             *
!     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
!     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
!                                                                             *
!******************************************************************************


SUBROUTINE ORRTM_TAUMOL1 (KLEV,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
  &COLH2O,LAYTROP,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.
!     Revised by Michael J. Iacono, Atmospheric & Environmental Research.

!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
 
! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up


#include "tsmbkind.h"

USE OPARRRTM  , ONLY : JPLAY  ,JPBAND ,JPGPT  ,JPXSEC ,NG1
USE OYOERRTWN , ONLY : NG     ,NSPA   ,NSPB
USE OYOERRTA1 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
            &FORREF   ,KA     ,KB     ,SELFREF 

!#include "yoeratm.h"

!      REAL TAUAER(JPLAY)

IMPLICIT NONE

!  Output
REAL_B :: TAU   (JPGPT,JPLAY)

!     DUMMY INTEGER SCALARS
INTEGER_M :: KLEV

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
INTEGER_M :: LAYTROP

!- from SELF             
REAL_B :: SELFFAC(JPLAY)
REAL_B :: SELFFRAC(JPLAY)
INTEGER_M :: INDSELF(JPLAY)

!- from SP             
REAL_B :: PFRAC(JPGPT,JPLAY)

INTEGER_M :: IND0(JPLAY),IND1(JPLAY),INDS(JPLAY)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IG, LAY

!      EQUIVALENCE (TAUAERL(1,1),TAUAER)

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.  

DO LAY = 1, LAYTROP
  IND0(LAY) = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(1) + 1
  IND1(LAY) = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(1) + 1
  INDS(LAY) = INDSELF(LAY)
ENDDO

DO IG = 1, NG1
  DO LAY = 1, LAYTROP
!-- DS_000515  
    TAU (IG,LAY) = COLH2O(LAY) *&
     &(FAC00(LAY) * ABSA(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSA(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSA(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSA(IND1(LAY)+1,IG) +&
     &SELFFAC(LAY) * (SELFREF(INDS(LAY),IG) + &
     &SELFFRAC(LAY) *&
     &(SELFREF(INDS(LAY)+1,IG) - SELFREF(INDS(LAY),IG)))&
     &+ FORFAC(LAY) * FORREF(IG) ) &
     &+ TAUAERL(LAY,1)
    PFRAC(IG,LAY) = FRACREFA(IG)
  ENDDO
ENDDO

DO LAY = LAYTROP+1, KLEV
  IND0(LAY) = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(1) + 1
  IND1(LAY) = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(1) + 1
ENDDO

!-- JJM000517
DO IG = 1, NG1
  DO LAY = LAYTROP+1, KLEV
!-- JJM000517
    TAU (IG,LAY) = COLH2O(LAY) *&
     &(FAC00(LAY) * ABSB(IND0(LAY)  ,IG) +&
     & FAC10(LAY) * ABSB(IND0(LAY)+1,IG) +&
     & FAC01(LAY) * ABSB(IND1(LAY)  ,IG) +&
     & FAC11(LAY) * ABSB(IND1(LAY)+1,IG)&
     &+ FORFAC(LAY) * FORREF(IG) ) &
     &+ TAUAERL(LAY,1)
    PFRAC(IG,LAY) = FRACREFB(IG)
  ENDDO
ENDDO

RETURN
END SUBROUTINE ORRTM_TAUMOL1
