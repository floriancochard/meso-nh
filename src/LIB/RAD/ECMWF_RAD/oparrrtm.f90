!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:36
!-----------------------------------------------------------------
MODULE OPARRRTM


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!     Parameters relevant to AER's RRTM-LW radiation scheme

!     980714  JJMorcrette
!     021212  P Jabouille increase JPLAY for increase vertical level number
!     ------------------------------------------------------------------

!-- standard tropical INTEGER_M, PARAMETER :: JPLAY  = 64
!-- ATEX              INTEGER_M, PARAMETER :: JPLAY  = 83
!-- BOMEX             INTEGER_M, PARAMETER :: JPLAY  = 84
!-- OPEN_CELLS        INTEGER_M, PARAMETER :: JPLAY  = 63
!__ GATE_A,B,C        INTEGER_M, PARAMETER :: JPLAY  = 46

INTEGER_M, PARAMETER :: JPLAY  = 300

INTEGER_M, PARAMETER :: JPG    = 16
INTEGER_M, PARAMETER :: JPBAND = 16
INTEGER_M, PARAMETER :: JPXSEC = 4
INTEGER_M, PARAMETER :: JPINPX = 35
INTEGER_M, PARAMETER :: JPGPT  = 140

INTEGER_M, PARAMETER :: NG1  = 8
INTEGER_M, PARAMETER :: NG2  = 14
INTEGER_M, PARAMETER :: NG3  = 16
INTEGER_M, PARAMETER :: NG4  = 14
INTEGER_M, PARAMETER :: NG5  = 16
INTEGER_M, PARAMETER :: NG6  = 8
INTEGER_M, PARAMETER :: NG7  = 12
INTEGER_M, PARAMETER :: NG8  = 8
INTEGER_M, PARAMETER :: NG9  = 12
INTEGER_M, PARAMETER :: NG10 = 6
INTEGER_M, PARAMETER :: NG11 = 8
INTEGER_M, PARAMETER :: NG12 = 8
INTEGER_M, PARAMETER :: NG13 = 4
INTEGER_M, PARAMETER :: NG14 = 2
INTEGER_M, PARAMETER :: NG15 = 2
INTEGER_M, PARAMETER :: NG16 = 8

INTEGER_M, PARAMETER :: NGS1  = 8
INTEGER_M, PARAMETER :: NGS2  = 22
INTEGER_M, PARAMETER :: NGS3  = 38
INTEGER_M, PARAMETER :: NGS4  = 52
INTEGER_M, PARAMETER :: NGS5  = 68
INTEGER_M, PARAMETER :: NGS6  = 76
INTEGER_M, PARAMETER :: NGS7  = 88
INTEGER_M, PARAMETER :: NGS8  = 96
INTEGER_M, PARAMETER :: NGS9  = 108
INTEGER_M, PARAMETER :: NGS10 = 114
INTEGER_M, PARAMETER :: NGS11 = 122
INTEGER_M, PARAMETER :: NGS12 = 130
INTEGER_M, PARAMETER :: NGS13 = 134
INTEGER_M, PARAMETER :: NGS14 = 136
INTEGER_M, PARAMETER :: NGS15 = 138

!     ------------------------------------------------------------------
END MODULE OPARRRTM
