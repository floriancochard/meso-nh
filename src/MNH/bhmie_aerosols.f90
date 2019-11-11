!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ##########################
       MODULE MODI_BHMIE_AEROSOLS
!      ##########################
!
INTERFACE
!
      SUBROUTINE BHMIE_AEROSOLS( PWAVELENGTH, PREFINDEX_COAT, PREFINDEX_CORE, &
                                 HDSD, PCONC, PFRACVOL_CORE,                  &
                                 PEXTINCTION_COEF, PBACKSCAT_COEF,            &
                                 KRADIUS, PDMODAL, PSIG, PRADIUS              )
!
REAL,                     INTENT(IN) :: PWAVELENGTH ! EM wavelength
COMPLEX,                  INTENT(IN) :: PREFINDEX_COAT ! COAT Refraction index
COMPLEX,                  INTENT(IN) :: PREFINDEX_CORE ! CORE Refraction index
CHARACTER(LEN=*),         INTENT(IN) :: HDSD   ! Type de size distribution
REAL,                     INTENT(IN) :: PCONC  ! Particle concentration
REAL,                     INTENT(IN) :: PFRACVOL_CORE  ! Core fraction volume
!
REAL,                    INTENT(OUT) :: PEXTINCTION_COEF ! Extinction
REAL,                    INTENT(OUT) :: PBACKSCAT_COEF   ! BackScattering
!
INTEGER, OPTIONAL,        INTENT(IN) :: KRADIUS ! Number of radii
                                                ! to discretize the dsd
REAL, OPTIONAL,           INTENT(IN) :: PDMODAL ! Lognormal law shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PSIG    ! Lognormal law shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PRADIUS ! Monodispersed-dsd radius
!
END SUBROUTINE BHMIE_AEROSOLS
END INTERFACE
END MODULE MODI_BHMIE_AEROSOLS 
!
!     ##########################################################################
      SUBROUTINE BHMIE_AEROSOLS( PWAVELENGTH, PREFINDEX_COAT, PREFINDEX_CORE,  &
                                 HDSD, PCONC, PFRACVOL_CORE,                   &
                                 PEXTINCTION_COEF, PBACKSCAT_COEF,             &
                                 KRADIUS, PDMODAL, PSIG, PRADIUS               )
!     ##########################################################################
!
!!****  *BHMIE_AEROSOLS* - computes the intgration of the Mie parameters over
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
USE MODI_BHMIE_BHCOAT
!USE MODI_GAUHER
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                     INTENT(IN) :: PWAVELENGTH ! EM wavelength
COMPLEX,                  INTENT(IN) :: PREFINDEX_COAT ! COAT Refraction index
COMPLEX,                  INTENT(IN) :: PREFINDEX_CORE ! CORE Refraction index
CHARACTER(LEN=*),         INTENT(IN) :: HDSD   ! Type de size distribution
REAL,                     INTENT(IN) :: PCONC  ! Particle concentration
REAL,                     INTENT(IN) :: PFRACVOL_CORE  ! Core fraction volume
!
REAL,                    INTENT(OUT) :: PEXTINCTION_COEF ! Extinction
REAL,                    INTENT(OUT) :: PBACKSCAT_COEF   ! BackScattering
!
INTEGER, OPTIONAL,        INTENT(IN) :: KRADIUS ! Number of radii
                                                ! to discretize the dsd
REAL, OPTIONAL,           INTENT(IN) :: PDMODAL ! Lognormal law shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PSIG    ! Lognormal law shape parameter
REAL, OPTIONAL,           INTENT(IN) :: PRADIUS ! Monodispersed-dsd radius
!
!*       0.2   Declarations of local variables :
!
INTEGER :: J
!
REAL    :: ZEXT             ! Extinction efficiency
REAL    :: ZBAK             ! BackScattering efficiency
REAL    :: ZSIZE_PARAM      ! Size parameter
REAL    :: ZSIZE_PARAM_CORE ! Size parameter of the solid aerosol core
REAL    :: ZLAMBDA          ! Dsd slope parameter
REAL    :: ZSURF_FACTOR     ! Surface parameter
!
REAL,    DIMENSION(:), ALLOCATABLE ::  ZABSCISSI,ZWEIGHTS
!
!-------------------------------------------------------------------------------
!
PEXTINCTION_COEF = 0.0 
PBACKSCAT_COEF   = 0.0
!
IF( HDSD=="LOGNO" ) THEN
  ALLOCATE(ZABSCISSI(KRADIUS))
  ALLOCATE(ZWEIGHTS(KRADIUS))
  CALL GAUHER(ZABSCISSI,ZWEIGHTS,KRADIUS)
  ZLAMBDA=1.0/PDMODAL
  ZSURF_FACTOR = 0.25*XPI*PCONC/(ZLAMBDA)**2
  DO J = 1,KRADIUS
    ZSIZE_PARAM = (XPI/PWAVELENGTH)*                         &
                  (1.0/ZLAMBDA)*EXP(ZABSCISSI(J)*SQRT(2.0)*PSIG)
    ZSIZE_PARAM_CORE = ZSIZE_PARAM*(PFRACVOL_CORE)**(1.0/3.0)
    CALL BHMIE_BHCOAT(ZSIZE_PARAM_CORE,ZSIZE_PARAM,    &
                        PREFINDEX_CORE,PREFINDEX_COAT, &
                                              ZEXT,ZBAK)
    ZWEIGHTS(J) = ZWEIGHTS(J)*ZSURF_FACTOR*EXP(SQRT(8.0)*PSIG*ZABSCISSI(J))
    PEXTINCTION_COEF = PEXTINCTION_COEF+ZWEIGHTS(J)*ZEXT
    PBACKSCAT_COEF   = PBACKSCAT_COEF  +ZWEIGHTS(J)*ZBAK
  END DO
ELSE IF( HDSD=="MONOD" ) THEN
  ZSIZE_PARAM = (2.E0*XPI*PRADIUS/PWAVELENGTH)
  ZSIZE_PARAM_CORE = ZSIZE_PARAM*(PFRACVOL_CORE)**(1.0/3.0)
  CALL BHMIE_BHCOAT(ZSIZE_PARAM_CORE,ZSIZE_PARAM,    &
                      PREFINDEX_CORE,PREFINDEX_COAT, &
                                            ZEXT,ZBAK)
  ZSURF_FACTOR = (XPI*PRADIUS**2)*PCONC
  PEXTINCTION_COEF = ZSURF_FACTOR*ZEXT
  PBACKSCAT_COEF   = ZSURF_FACTOR*ZBAK
END IF
!
END SUBROUTINE BHMIE_AEROSOLS
