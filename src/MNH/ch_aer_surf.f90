!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!   ########################
     MODULE MODI_CH_AER_SURF
!!   ########################
!!
INTERFACE
!!
SUBROUTINE CH_AER_SURF(PM, PRG, PSIG, PSURF)

IMPLICIT NONE
REAL,   DIMENSION(:,:), INTENT(IN) :: PM      ! moments
REAL,   DIMENSION(:,:), INTENT(IN) :: PRG     ! radius
REAL,   DIMENSION(:,:), INTENT(IN) :: PSIG    ! dispersion
REAL,   DIMENSION(:,:), INTENT(OUT) :: PSURF  ! aerosol surface
END SUBROUTINE CH_AER_SURF
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_SURF
!!
!!   ############################################################
     SUBROUTINE CH_AER_SURF( PM, PRG, PSIG,  PSURF)
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!      Compute surface aerosol 
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    P. Tulet (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CH_AEROSOL, ONLY :  JPMODE, NM0, NM3, NM6
USE MODD_CST,        ONLY :  XPI 
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:), INTENT(IN) :: PM      ! moments
REAL,   DIMENSION(:,:), INTENT(IN) :: PRG     ! radius
REAL,   DIMENSION(:,:), INTENT(IN) :: PSIG    ! dispersion
REAL,   DIMENSION(:,:), INTENT(OUT) :: PSURF  ! aerosol surface
!
!
!*      0.2    declarations local variables
!
REAL,   DIMENSION(SIZE(PM,1),JPMODE) :: ZSIG0, ZRG0, ZN0
REAL,   DIMENSION(SIZE(PM,1),SIZE(PM,2)) :: ZM
INTEGER :: JN
!
!-------------------------------------------------------------------------------
!
!
DO JN=1, JPMODE
  ZM(:,NM0(JN)) = PM(:,NM0(JN)) ! 1/m3
  ZM(:,NM3(JN)) = PM(:,NM3(JN)) ! 
  ZM(:,NM6(JN)) = PM(:,NM6(JN)) 

  ZN0(:,JN) = ZM(:,NM0(JN)) * 1E-6 ! convert from 1/m3 to 1/cc
  ZRG0(:,JN)= PRG(:,JN) * 1E-6   ! convert from micrometers to meters 
  ZSIG0(:,JN)= PSIG(:,JN)
!
! Surface = Pi * M2 with M2 = N Rg**2 exp (2*(ln(sigma)**2))
!
!surface in m2 / cc (psurf will be convert into SI units in BASIC)
  PSURF(:,JN) =  XPI * ZN0(:,JN) * (ZRG0(:,JN))**2 * EXP(2.*ZSIG0(:,JN)**2)

ENDDO
!
!
END SUBROUTINE CH_AER_SURF
!
