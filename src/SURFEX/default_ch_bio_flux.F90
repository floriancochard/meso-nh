!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE DEFAULT_CH_BIO_FLUX(OCH_BIO_FLUX,PDAILYPAR,PDAILYTEMP)
!     ########################################################################
!
!!****  *DEFAULT_CH_BIO_FLUX* - routine to set default values for the configuration for CH_BIO_FLUX scheme
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2007 
!!      J.Pianezzej 02/2019 : correction for use of MEGAN
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
LOGICAL,  INTENT(OUT) :: OCH_BIO_FLUX  ! flag for the calculation of biogenic fluxes
REAL,     INTENT(OUT), OPTIONAL :: PDAILYPAR  ! default values for megan PAR   temperature
REAL,     INTENT(OUT), OPTIONAL :: PDAILYTEMP ! default values for megan daily temperature
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_CH_BIO_FLUX',0,ZHOOK_HANDLE)
OCH_BIO_FLUX= .FALSE.
IF (PRESENT(PDAILYPAR))  PDAILYPAR   = 200.
IF (PRESENT(PDAILYTEMP)) PDAILYTEMP  = 293.
IF (LHOOK) CALL DR_HOOK('DEFAULT_CH_BIO_FLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_CH_BIO_FLUX
