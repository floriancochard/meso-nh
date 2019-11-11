!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE OYOMLUN


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Logical units used by code

!     NULOUT :   output unit
!     NULNAM :   unit number for namelist
!     NORZUN :   UNIT NO. FOR TEMP HEIGHT ERROR CORRELATION
!     NCMAFL :   UNIT NUMBERS FOR CMA FILES
!     NULCL1 :   unit number for climatological fields (month before)
!     NULCL2 :   unit number for climatological fields (month after)
!     NTRJSH :   unit number for trajectory spectral data          WRTRA
!     NTRJGG :   unit number for trajectory grid-point data        WRTRA
!     NINMSH :   unit number for initial point of the minimization SUVAZX
!     NINMGG :   unit number for initial point of the minimization SUVAZX
!     NINISH :   unit number for initial spectral data             SUSPEC
!     NINIGG :   unit number for initial grid-point data           SUSPEC
!     NFGISH :   unit number for first-guess spectral data
!     NFGIGG :   unit number for first-guess grid-point data
!     NPOSSH :   output unit number (spectral fields)              CPREP1
!     NPOSGG :   output unit number (grid point fields)            CPREP1
!     NSPHMO :   unit number for conversion matrix (to normal mode space)
!     NHOUMO :   unit number for conversion matrix (from normal mode space)
!     NPINMI :   unit number for partially implicit NMI help arrays
!     NFREQ2 :   unit number for freqencies arrays FREQ and FREQ2 (YOMNMI)
!     NSAVSP :   unit number for saving spectral arrays           (CNMI)
!     NTIDE  :   unit number for the LFI file containing the total tendencies
!     NLTYPE :   unit number for LTYPGR, LTYPDB and LTYPTD arrays (YOMNMI)
!     NPDIRL :   unit number for post-processing directory listing
!     NPPPSH :   unit number for post-processed spherical harmonics WRPLPP
!     NPPGG  :   unit number for post-processed grid-point fields   WRPLPP
!     NFCER0 :   unit number for forecast errors on standard error grid
!     NFCERZ :   unit number for horizontally interpolated forecast errors
!                in field form.
!     NFCERH :   unit number for horizontally interpolated forecast errors
!                in line form.
!     NFCEVC :   unit number for compressed error climatology file
!     NFCSPG :   unit number for spectral first guess field

!     NFFUNC :   unit number for Fspalt(1,jlev,jvar) and Fspsur(1,1) with the
!                number of simulations.                             EVCOST
!     NFSPNC :   unit number for Fspalt(jlev,jvar,jn) and Fspsur(1,jn)
!                with the number of simulations.                    EVCOST
!     NFDIT0 :   unit number for global distance at t0 with the
!                number of simulations.                              SIM4D
!     NFSPT0 :   unit number for distance at t0 (by field and wave number)
!                at the begining and at the end of the minimisation COSTRA

!     NULZDID:   unit number for zonal diagnostics of dynamics (ZODIA).
!     NULZDIF:   unit number for zonal diagnostics of physical fluxes (ZODIA).

!     NULDILA:   unit number for dilatation matrix (SUDIL,DILAT,SPDILA)
!     NULCONT:   unit number for contraction matrix (SUDIL,DILAT,SPDILA)
!     NULROTS:   unit number for lower troncature rotation matrix (SUROT,SPORTS)
!     NULROTC:   unit number for upper troncature rotation matrix (SUROT,SPORTS)

!     NULLEG :   unit number for Legendre polynomials
!     NULCO  :   unit number for coupled fields (ICMCO)
!     NPODDH :   unit number for mask diagnostic files (DDH)
!     NULRCF :   unit number for restart control file
!     NULHWF :   unit number for history witness file
!     NBIAS  :   unit number for bias (dig. filt. guess - guess)

!     NEFLS  :   unit number for coupling ALADIN file
!     NEFLSS :   unit number for coupling ALADIN file (initialisation)

!     NULUSR1:   unit numbers for user defined files
!     NULDISP:   unit number for display file
!     NULSTAT:   unit number for status  file

!     NULASE :   unit number for CANARI statistics (forecast error s.d.)
!     NULASS :   unit number for CANARI statistics (analysis error s.d.)

!     NULUSR2
!     NULUSR3
!     NULUSR4
!     NULUSR5
!     NULTMP :   unit numbers for file opened and closed in the same routine

!     NULFPxx    unit numbers for Full-POS output files
!     NSCRTCH:   unit number for Full-POS scratch file (for in-line post-proc.)
!     NULFPOS:   unit number for Full-POS control file (end of post-processing
!                in conf. 001 ; auxilary namelist file in conf. 927)
!     NULDIA :   unit number for error diagnostics
!     NULERR :   unit number for comparison with reference run
!     NULREF :   unit number for storing reference run
!     NULTIM :   unit number for storing timing information for DM run
!     NULRAD :   unit number for writing radiation diagnostics
!     NUO3CH1:   unit number for reading ozone chemistry file 1
!     NUO3CH2:   unit number for reading ozone chemistry file 2
!     NTCSR  :   unit number for fields of radiation coefficients
INTEGER_M :: NULOUT
INTEGER_M :: NULNAM
INTEGER_M :: NCMAFL(26)
INTEGER_M :: NPOSSH
INTEGER_M :: NPOSGG
INTEGER_M :: NSPHMO
INTEGER_M :: NHOUMO
INTEGER_M :: NPINMI
INTEGER_M :: NFREQ2
INTEGER_M :: NLTYPE
INTEGER_M :: NTIDE
INTEGER_M :: NTRJSH
INTEGER_M :: NTRJGG
INTEGER_M :: NINMSH
INTEGER_M :: NINMGG
INTEGER_M :: NINISH
INTEGER_M :: NINIGG
INTEGER_M :: NFGISH
INTEGER_M :: NFGIGG
INTEGER_M :: NPPPSH
INTEGER_M :: NPPGG
INTEGER_M :: NULTMP
INTEGER_M :: NFCER0
INTEGER_M :: NFCERZ
INTEGER_M :: NFCERH
INTEGER_M :: NFCEVC
INTEGER_M :: NFCSPG
INTEGER_M :: NPODDH
INTEGER_M :: NULLEG
INTEGER_M :: NORZUN
INTEGER_M :: NSAVSP
INTEGER_M :: NFFUNC
INTEGER_M :: NFSPNC
INTEGER_M :: NFDIT0
INTEGER_M :: NFSPT0
INTEGER_M :: NULZDID
INTEGER_M :: NULZDIF
INTEGER_M :: NULCL1
INTEGER_M :: NULCL2
INTEGER_M :: NULASE
INTEGER_M :: NULASS
INTEGER_M :: NULDILA
INTEGER_M :: NULCONT
INTEGER_M :: NULROTS
INTEGER_M :: NULROTC
INTEGER_M :: NULRCF
INTEGER_M :: NULHWF
INTEGER_M :: NULUSR1
INTEGER_M :: NULUSR2
INTEGER_M :: NULUSR3
INTEGER_M :: NULUSR4
INTEGER_M :: NULUSR5
INTEGER_M :: NULCO
INTEGER_M :: NEFLS
INTEGER_M :: NEFLSS
INTEGER_M :: NBIAS
INTEGER_M :: NPDIRL
INTEGER_M :: NULDISP
INTEGER_M :: NULSTAT
INTEGER_M :: NULFP01
INTEGER_M :: NULFP02
INTEGER_M :: NULFP03
INTEGER_M :: NULFP04
INTEGER_M :: NULFP05
INTEGER_M :: NULFP06
INTEGER_M :: NULFP07
INTEGER_M :: NULFP08
INTEGER_M :: NULFP09
INTEGER_M :: NULFP10
INTEGER_M :: NULFP11
INTEGER_M :: NULFP12
INTEGER_M :: NULFP13
INTEGER_M :: NULFP14
INTEGER_M :: NULFP15
INTEGER_M :: NULFPOS
INTEGER_M :: NSCRTCH
INTEGER_M :: NULDIA
INTEGER_M :: NULERR
INTEGER_M :: NULREF
INTEGER_M :: NULTIM
INTEGER_M :: NULRAD
INTEGER_M :: NUO3CH1
INTEGER_M :: NUO3CH2
INTEGER_M :: NTCSR
!     ------------------------------------------------------------------
END MODULE OYOMLUN
