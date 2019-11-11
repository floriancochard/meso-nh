!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #######################
       MODULE MODI_BHMIE_WATER
!      #######################
!
INTERFACE
!
      SUBROUTINE BHMIE_WATER( PWAVELENGTH, PREFINDEX, HDSD, PCONC,      & 
                              KANGLE, PEXTINCTION_COEF, PBACKSCAT_COEF, &
                              KRADIUS, PALPHA, PNU, PLWC, PRADIUS       )
!
REAL,                     INTENT(IN) :: PWAVELENGTH ! EM wavelength
COMPLEX,                  INTENT(IN) :: PREFINDEX   ! Refraction index
CHARACTER(LEN=*),         INTENT(IN) :: HDSD   ! Type de size distribution
REAL,                     INTENT(IN) :: PCONC  ! Particle concentration
INTEGER,                  INTENT(IN) :: KANGLE ! Number of scattering angles
!
REAL,                    INTENT(OUT) :: PEXTINCTION_COEF ! Extinction
REAL,                    INTENT(OUT) :: PBACKSCAT_COEF   ! BackScattering
!
INTEGER, OPTIONAL,        INTENT(IN) :: KRADIUS ! Number of radii
                                                ! to discretize the dsd
REAL, OPTIONAL,           INTENT(IN) :: PALPHA  ! Gamma-dsd shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PNU     ! Gamma-dsd shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PLWC    ! Gamma-dsd liquid water content
REAL, OPTIONAL,           INTENT(IN) :: PRADIUS ! Monodispersed-dsd radius
!
END SUBROUTINE BHMIE_WATER
END INTERFACE
END MODULE MODI_BHMIE_WATER 
!
!     ##################################################################
      SUBROUTINE BHMIE_WATER( PWAVELENGTH, PREFINDEX, HDSD, PCONC,      & 
                              KANGLE, PEXTINCTION_COEF, PBACKSCAT_COEF, &
                              KRADIUS, PALPHA, PNU, PLWC, PRADIUS       )
!     ##################################################################
!
!!****  *BHMIE_WATER* - computes the intgration of the Mie parameters over
!!                      the drop size distributions 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the 
!!
!!**  METHOD
!!    ------
!!      The reflectivities are computed using the n(D) * D**6 formula. The
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XRHOLW               ! Liquid water density
!!      Module MODD_RAIN_ICE_DESCR
!!      Module MODD_RAIN_ICE_PARAM
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      J.-P. Chaboureau * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/04/07
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODI_BHMIE
USE MODI_GAMMA
!USE MODI_GAULAG
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN) :: PWAVELENGTH ! EM wavelength
COMPLEX,                  INTENT(IN) :: PREFINDEX   ! Refraction index
CHARACTER(LEN=*),         INTENT(IN) :: HDSD   ! Type de size distribution
REAL,                     INTENT(IN) :: PCONC  ! Particle concentration
INTEGER,                  INTENT(IN) :: KANGLE ! Number of scattering angles
!
REAL,                    INTENT(OUT) :: PEXTINCTION_COEF ! Extinction
REAL,                    INTENT(OUT) :: PBACKSCAT_COEF   ! BackScattering
!
INTEGER, OPTIONAL,        INTENT(IN) :: KRADIUS ! Number of radii
                                                ! to discretize the dsd
REAL, OPTIONAL,           INTENT(IN) :: PALPHA  ! Gamma-dsd shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PNU     ! Gamma-dsd shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PLWC    ! Gamma-dsd liquid water content
REAL, OPTIONAL,           INTENT(IN) :: PRADIUS ! Monodispersed-dsd radius
!
!*       0.2   Declarations of local variables :
!
INTEGER :: J
!
REAL    :: ZDELTANGLE ! Angle increment
REAL    :: ZNUM1      ! PNU-1.
REAL    :: ZEXT       ! Extinction efficiency
REAL    :: ZBAK       ! BackScattering efficiency
REAL    :: ZSIZE_PARAM ! Size parameter
REAL    :: ZLAMBDA     ! Dsd slope parameter
REAL    :: ZSHAPE_FACTOR ! Shape parameter
REAL    :: ZSURF_FACTOR  ! Surface parameter
!
COMPLEX, DIMENSION(:), ALLOCATABLE ::  ZZS1,ZZS2
REAL,    DIMENSION(:), ALLOCATABLE ::  ZABSCISSI,ZWEIGHTS
!
!-------------------------------------------------------------------------------
!
ZDELTANGLE=0.5E0*XPI/FLOAT(KANGLE-1)
ALLOCATE(ZZS1(2*KANGLE-1))
ALLOCATE(ZZS2(2*KANGLE-1))
PEXTINCTION_COEF = 0.0 
PBACKSCAT_COEF   = 0.0
!
IF( HDSD=="GAMMA" ) THEN
  ALLOCATE(ZABSCISSI(KRADIUS))
  ALLOCATE(ZWEIGHTS(KRADIUS))
  ZNUM1=PNU-1.0E0
  CALL GAULAG(ZABSCISSI,ZWEIGHTS,KRADIUS,ZNUM1)
  ZSHAPE_FACTOR = GAMMA(PNU+3.0/PALPHA)/GAMMA(PNU)
  ZLAMBDA = ((XPI/6.E0)*(PCONC/PLWC)*ZSHAPE_FACTOR)**(1.0/3.0)
  ZSURF_FACTOR = 0.25*XPI*PCONC/(ZLAMBDA)**2
  DO J = 1,KRADIUS
    ZSIZE_PARAM = (XPI/PWAVELENGTH)*                         &
                  ((1.0/ZLAMBDA)*(ZABSCISSI(J))**(1.0/PALPHA))
    CALL BHMIE(ZSIZE_PARAM,PREFINDEX,KANGLE,ZZS1,ZZS2,ZEXT,ZBAK)
    ZWEIGHTS(J) = ZWEIGHTS(J)*ZSURF_FACTOR*(ZABSCISSI(J))**(2.0/PALPHA)
    PEXTINCTION_COEF = PEXTINCTION_COEF+ZWEIGHTS(J)*ZEXT
    PBACKSCAT_COEF   = PBACKSCAT_COEF  +ZWEIGHTS(J)*ZBAK
  END DO
ELSE IF( HDSD=="MONOD" ) THEN
  ZSIZE_PARAM = (2.E0*XPI*PRADIUS/PWAVELENGTH)
  ZSURF_FACTOR = (XPI*PRADIUS**2)*PCONC
  CALL BHMIE(ZSIZE_PARAM,PREFINDEX,KANGLE,ZZS1,ZZS2,ZEXT,ZBAK)
  PEXTINCTION_COEF = ZSURF_FACTOR*ZEXT
  PBACKSCAT_COEF   = ZSURF_FACTOR*ZBAK
END IF
!
END SUBROUTINE BHMIE_WATER
