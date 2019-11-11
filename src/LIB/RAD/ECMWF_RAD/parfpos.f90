!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE PARFPOS

#include "tsmbkind.h"

USE PARDIM


IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! === basic dimensions for Full POST-PROCESSING ===

!     JPOSDOM : Maximum number of horizontal (sub)domains
!     JPOSFRQ : Maximum number of output frequencies
!     JPOSLEN : Maximum length of a (sub)domain name
!     JPOSLIS : Maximum number of groups of subdomains
!     JPOSDIR : Maximum length of the path (or prefix) for the output files
!     JPOSLE  : Maximum number of eta levels on the output subdomain
!     JPOSGL  : Maximum number of latitude rows of the output gaussian grid
!     JPOS3DF : Maximum number of specific 3D dynamic fields
!     JPOSSCVA: Maximum number of post-processable passive scalars
!     JPOS2DF : Maximum number of specific 2D dynamic fields
!     JPOSPHY : Maximum number of surface fields
!     JPOSCFU : Maximum number of cumulated fluxes
!     JPOSXFU : Maximum number of instantaneous fluxes
!     JPOS3P  : Maximum number of pp. pressure levels
!     JPOS3H  : Maximum number of pp. height (above orography) levels
!     JPOS3TH : Maximum number of pp. potential temperature levels
!     JPOS3PV : Maximum number of pp. potential vorticity levels
!     JPOS3S  : Maximum number of pp. eta levels
!     JPOSVSO : Maximum number of climatologic fields of output format


INTEGER_M, PARAMETER :: JPOSDOM=15
INTEGER_M, PARAMETER :: JPOSFRQ=10
INTEGER_M, PARAMETER :: JPOSLEN=10
INTEGER_M, PARAMETER :: JPOSLIS=10
INTEGER_M, PARAMETER :: JPOSDIR=180
INTEGER_M, PARAMETER :: JPOSLE=JPMXLE
INTEGER_M, PARAMETER :: JPOSGL=JPMXGL
INTEGER_M, PARAMETER :: JPOS3DF=63
INTEGER_M, PARAMETER :: JPOS2DF=15
INTEGER_M, PARAMETER :: JPOSPHY=127
INTEGER_M, PARAMETER :: JPOSCFU=63
INTEGER_M, PARAMETER :: JPOSXFU=63
INTEGER_M, PARAMETER :: JPOS3P=31
INTEGER_M, PARAMETER :: JPOS3H=127
INTEGER_M, PARAMETER :: JPOS3TH=15
INTEGER_M, PARAMETER :: JPOS3PV=15
INTEGER_M, PARAMETER :: JPOS3S=JPOSLE
INTEGER_M, PARAMETER :: JPOSVSO=31
INTEGER_M, PARAMETER :: JPOSSCVA=3
!     ------------------------------------------------------------------
END MODULE PARFPOS
