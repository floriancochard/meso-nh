!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:42
!-----------------------------------------------------------------
SUBROUTINE SUCLOP

!**** *SUCLOP*  - INITIALIZE COMMON YOECLOP

!     PURPOSE.
!     --------
!           INITIALIZE YOMCLOP, WITH CLOUD OPTICAL PARAMETERS

!**   INTERFACE.
!     ----------
!        *CALL*  SUCLOP
!     FROM *SUECRAD*

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECLOP

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 92-02-29
!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE YOECLOP  , ONLY : RYFWCA   ,RYFWCB   ,RYFWCC   ,RYFWCD   ,&
            &RYFWCE   ,RYFWCF   ,REBCUA   ,REBCUB   ,REBCUC   ,&
            &REBCUD   ,REBCUE   ,REBCUF   ,REBCUG   ,REBCUH   ,&
            &REFFIA   ,REFFIB   ,RTIW     ,RRIW


IMPLICIT NONE


!      ----------------------------------------------------------------

!* Ice cloud properties - crystal: adapted from Ebert and Curry, 1992

! SW : 2 spectral intervals

REBCUA(1)= 3.448E-03_JPRB
REBCUA(2)= 3.448E-03_JPRB
REBCUB(1)= 2.431_JPRB
REBCUB(2)= 2.431_JPRB
REBCUC(1)= 0.99999_JPRB
REBCUC(2)= 0.975634_JPRB
REBCUD(1)= 0._JPRB
REBCUD(2)= 2.487E-04_JPRB
REBCUE(1)= 0.7661_JPRB
REBCUE(2)= 0.7866_JPRB
REBCUF(1)= 5.851E-04_JPRB
REBCUF(2)= 5.937E-04_JPRB

! LW : spectrally averaged with reference Planck function at 257 K

REBCUG= 1.07677_JPRB
REBCUH= 0.00267_JPRB

! Ice particle Effective Radius as a function of LWC

REFFIA= 40._JPRB
REFFIB=  0._JPRB

!* Water cloud properties - from Fouquart (1987)

! SW : 2 spectral intervals: parameters as a function of Reff

RYFWCA(1)= 0._JPRB
RYFWCA(2)= 0._JPRB
RYFWCB(1)= 1.5_JPRB
RYFWCB(2)= 1.5_JPRB
RYFWCC(1)= 0.9999_JPRB
RYFWCC(2)= 0.9988_JPRB
RYFWCD(1)= 5.000E-04_JPRB
RYFWCD(2)= 2.500E-03_JPRB
RYFWCE(1)= 0.5_JPRB
RYFWCE(2)= 0.05_JPRB
RYFWCF(1)= 0.865_JPRB
RYFWCF(2)= 0.910_JPRB

!* Liquid/Solid water transition

RTIW= 263._JPRB
RRIW= 20._JPRB

!     ------------------------------------------------------------------

END SUBROUTINE SUCLOP
