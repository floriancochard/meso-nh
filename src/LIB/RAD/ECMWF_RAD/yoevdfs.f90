!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:50
!-----------------------------------------------------------------
MODULE YOEVDFS


#include "tsmbkind.h"

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEVDFS* CONTAINS STABILITY FUNCTION TABLES FOR *VDF...*
!     ------------------------------------------------------------------

INTEGER_M, PARAMETER :: JPRITBL=101
REAL_B :: RITBL(JPRITBL)
REAL_B :: ARITBL(JPRITBL)
REAL_B :: RCHBA
REAL_B :: RCHBB
REAL_B :: RCHBC
REAL_B :: RCHBD
REAL_B :: RCHB23A
REAL_B :: RCHBBCD
REAL_B :: RCHBCD
REAL_B :: RCHETA
REAL_B :: RCHETB
REAL_B :: RCHETC
REAL_B :: RCHBHDL
REAL_B :: RCDHALF
REAL_B :: RCDHPI2
REAL_B :: RIMAX
REAL_B :: DRITBL
REAL_B :: DRI26

!**   ** *YOEDFS* CONTAINS STABILITY FUNCTION TABLES FOR *VDF...*

!     A.C.M. BELJAARS   E.C.M.W.F.       26/03/90.

!      NAME      TYPE        PURPOSE
!      ----      ----        -------

!     *RCHBA*     REAL       *CONSTANT A IN *HOLTSLAG AND *DEBRUIN
!                            FUNCTIONS FOR STABLE SITUATIONS
!     *RCHBB*     REAL       *CONSTANT B IN *HB* FUNCTIONS
!     *RCHBC*     REAL       *CONSTANT C IN *HB* FUNCTIONS
!     *RCHBD*     REAL       *CONSTANT D IN *HB* FUNCTIONS
!     *RCHB23A    REAL       2./3.*A IN *HB* FUNCTIONS
!     *RCHBBCD    REAL       B*C/D IN *HB* FUNCTIONS
!     *RCHBCD     REAL       C/D IN *HB* FUNCTIONS
!     *RCHETA*    REAL       CONSTANT IN THE *HOGSTROM *ELLISON *TURNER
!                            FUNCTIONS FOR STABLY STRATIFIED TURBULENCE 
!     *RCHETB*    REAL       CONSTANT IN THE *HET* FUNCTIONS     
!     *RCHETC*    REAL       CONSTANT IN *HET* FUNCTIONS  
!     *RCHBHDL*   REAL       MAXIM ZNLEV/L FOR STABLE BOUNDARY LAYER  
!     *RCDHALF    REAL       CONSTANT IN *DYER AND *HICKS FORMULAE
!                            FOR UNSTABLE SITUATIONS
!     *RCDHPI2    REAL       PI/2.
!     *RIMAX*     REAL       *MAXIMIM RICHARDSON NUMBER TABULATED
!     *DRITBL*    REAL       *INCREMENT OF THE RICHARDSON NUMBER
!                            BETWEEN TABULATED VALUES.
!     *DRI26*     REAL       DRITBL**2/6.
!     *RITBL*     REAL ARRAY *TABULATED ETA-VALUES (Z/L) AS A FUNCTION
!                            OF THE RICHARDSON NUMBER FOR STABLE CASES.
!     *ARITBL*    REAL ARRAY *SECOND DERIVATIVES OF TABULATED FUNCTION
!                            FOR SPLINE INTERPOLATION.
!     ------------------------------------------------------------------
END MODULE YOEVDFS
