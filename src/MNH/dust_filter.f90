!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!!   ##############################
     MODULE MODI_DUST_FILTER
!!   ##############################
!!
INTERFACE
!
SUBROUTINE DUST_FILTER(PSV, PRHODREF)

IMPLICIT NONE

REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF

END SUBROUTINE DUST_FILTER
!!
END INTERFACE
!!
END MODULE MODI_DUST_FILTER
!!
!!   #######################################
     SUBROUTINE DUST_FILTER(PSV, PRHODREF)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!    Pierre TULET (CNRM/GMEI) 
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
! Entry variables:
!
! PRSVS(INOUT)       -Array of moments included in PRSVS
!
!*************************************************************
! Exit variables:
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! ZVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!!   IMPLICIT ARGUMENTS
!
USE MODD_DUST
USE MODD_CSTS_DUST
USE MODD_CST, ONLY : XMNH_TINY
!!
IMPLICIT NONE

!! Declarations d'arguments

REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF

!! Declarations de variables internes

INTEGER :: JN
INTEGER :: JMODEIDX
REAL,    DIMENSION(NMODE_DST*3) :: ZPMIN
REAL,    DIMENSION(NMODE_DST)   :: ZINIRADIUS
REAL,    DIMENSION(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_DST*3)  :: ZM                  ! [aerosol units] local array which goes to output later
REAL,    DIMENSION(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3))  :: ZSIGMA !standard deviation

REAL :: ZRGMIN, ZSIGMIN
REAL :: ZPI, ZRHOP, ZFAC, ZMI
INTEGER,DIMENSION(NMODE_DST) :: NM0, NM3, NM6

PSV(:,:,:,:) =  MAX(PSV(:,:,:,:), XMNH_TINY)


END SUBROUTINE DUST_FILTER
