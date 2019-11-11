
!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ########################
       MODULE MODI_ELEC_FIELD_n
!      ########################
!
INTERFACE
!
      SUBROUTINE ELEC_FIELD_n (PQ_SOURCE, KTCOUNT, PRELAX, PRHODJ, &
                               PEFIELDU, PEFIELDV, PEFIELDW, PPHIT)
!
INTEGER,                INTENT(IN)  :: KTCOUNT ! counter value of the 
                                               ! model temporal loop
REAL,                   INTENT(IN)  :: PRELAX ! relaxation coefficient for 
                                              ! the Richardson's method
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHODJ ! density of reference state
                                             ! * J
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PQ_SOURCE ! Electric charge
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEFIELDU !  3 components 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEFIELDV !     of the 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEFIELDW ! electric field
REAL, DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: PPHIT ! Electrostatic potential
!
END SUBROUTINE ELEC_FIELD_n
END INTERFACE
END MODULE MODI_ELEC_FIELD_n
!
!     ##################################################################
      SUBROUTINE ELEC_FIELD_n(PQ_SOURCE, KTCOUNT, PRELAX, PRHODJ, &
                              PEFIELDU, PEFIELDV, PEFIELDW, PPHIT)
!     ##################################################################
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY : JPVEXT
USE MODD_CONF, ONLY : LCARTESIAN, LFLAT
USE MODD_CST, ONLY : XCPD
USE MODD_DYN_n, ONLY: XTRIGSX, XTRIGSY, XDXHATM, XDYHATM, NIFAXX, NIFAXY
USE MODD_METRICS_n
USE MODD_ELEC_DESCR, ONLY : CLSOL, NLAPITR_ELEC, XEPSILON, XEPOTFW_TOP
USE MODD_ELEC_n, ONLY : XRHOM_E, XAF_E, XCF_E, XBFY_E, &
                        XBFB_E, XBF_SXP2_YP1_Z_E
!
USE MODI_FLAT_INV
USE MODI_FLAT_INVZ
USE MODI_RICHARDSON
USE MODI_CONJGRAD
USE MODI_CONRESOL
USE MODI_CONRESOLZ
USE MODI_GRADIENT_M
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
USE MODE_ll
USE MODE_FM
!
!
IMPLICIT NONE
!
!*       0.1   declarations of dummy arguments
!
INTEGER,                INTENT(IN)  :: KTCOUNT ! counter value of the 
                                               ! model temporal loop
REAL,                   INTENT(IN)  :: PRELAX ! relaxation coefficient for 
                                              ! the Richardson's method
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHODJ ! density of reference state
                                              ! * J
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PQ_SOURCE ! Electric charge
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEFIELDU !  3 components 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEFIELDV !     of the 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEFIELDW ! electric field
REAL, DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: PPHIT ! Electrostatic potential
!
!
!*       0.2   declarations of local variables 
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZPHIT  ! electrostatic potential
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDV_SOURCE  ! divergence of the sources
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTHETAV  ! virtual potential temperature
!
CHARACTER (LEN=4), DIMENSION(2)             :: ZLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2)             :: ZLBCY    ! y-direction LBC type 
!
INTEGER :: IIB      ! indice I for the first inner mass point along x
INTEGER :: IIE      ! indice I for the last inner mass point along x
INTEGER :: IJB      ! indice J for the first inner mass point along y
INTEGER :: IJE      ! indice J for the last inner mass point along y
INTEGER :: IKB      ! indice K for the first inner mass point along z
INTEGER :: IKE      ! indice K for the last inner mass point along z
INTEGER :: IRESP    ! Return code of FM routines
INTEGER :: ILENG    ! Length of the data field in LFIFM file
INTEGER :: IGRID    ! C-grid indicator in LFIFM file
INTEGER :: ILENCH   ! Length of comment string in LFIFM file
INTEGER :: IIU,IJU,IKU ! array sizes in I,J,K
INTEGER :: JI,JJ,JK    !loop index on the vertical levels
INTEGER :: IINFO_ll  
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!
! 
!------------------------------------------------------------------------------
!
NULLIFY(TZFIELDS_ll)
!
!*       1.    COMPUTE LOOP BOUNDS
!              -------------------
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
CALL GET_DIM_EXT_ll('B',IIU,IJU)
!
IKB = 1 + JPVEXT
IKU = SIZE(PEFIELDU,3)
IKE = IKU - JPVEXT
!
ALLOCATE(ZDV_SOURCE(SIZE(PQ_SOURCE,1),SIZE(PQ_SOURCE,2),SIZE(PQ_SOURCE,3)))
ALLOCATE(ZPHIT(SIZE(PQ_SOURCE,1),SIZE(PQ_SOURCE,2),SIZE(PQ_SOURCE,3)))
ALLOCATE(ZTHETAV(SIZE(PQ_SOURCE,1),SIZE(PQ_SOURCE,2),SIZE(PQ_SOURCE,3)))
!
ZPHIT(:,:,:) = 0.
ZTHETAV(:,:,:) = 1./XCPD ! to compensate for the Cpd*Thetav factor in the
                         ! pressure problem for which the solver was coded
! 
!------------------------------------------------------------------------------
!
!*       2.    BOUNDARIES CONDITIONS
!              ---------------------
!
! Lateral problem is open
! Space charge source is resolved in "ini_elec_field"
ZLBCX = 'OPEN' 
ZLBCY = 'OPEN'
!
! Source term
ZDV_SOURCE(:,:,:) = PQ_SOURCE(:,:,:) / XEPSILON
!
! The non-homogenous Neuman problem is transformed in an homogenous Neuman
! problem in the non-periodic cases !!!
IF (LWEST_ll() ) ZDV_SOURCE(IIB-1,:,:) = 0.
IF (LEAST_ll() ) ZDV_SOURCE(IIE+1,:,:) = 0.
IF (LSOUTH_ll()) ZDV_SOURCE(:,IJB-1,:) = 0.
IF (LNORTH_ll()) ZDV_SOURCE(:,IJE+1,:) = 0.
!
DO JJ = 1, IJU
  DO JI = 1, IIU
    ZDV_SOURCE(JI,JJ,IKE+1) = XEPOTFW_TOP(JI,JJ) ! Filled with E_fw/dZ which is
                                                 ! homogeneous to q/epsilon
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       3.    SOLVE THE ELECTRIC FIELD EQUATION
!              ----------------------------------
!
!*       3.1   Compute the electrostatic potential
!              -----------------------------------
!
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZPHIT)
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZDV_SOURCE)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
IF (LFLAT .AND. LCARTESIAN) THEN
! flat cartesian LHE case -> exact solution 
  IF ( CLSOL /= "ZRESI" ) THEN
    CALL FLAT_INV (ZLBCX, ZLBCY, XDXHATM, XDYHATM,   &
                   XRHOM_E, XAF_E, XBFY_E, XCF_E,    &
                   XTRIGSX, XTRIGSY, NIFAXX, NIFAXY, &
                   ZDV_SOURCE, ZPHIT)
  ELSE
    CALL FLAT_INVZ (ZLBCX, ZLBCY, XDXHATM, XDYHATM,  &
                   XRHOM_E, XAF_E, XBFY_E, XCF_E,    &
                   XTRIGSX, XTRIGSY, NIFAXX, NIFAXY, &
                   ZDV_SOURCE, ZPHIT,                &
                   XBFB_E, XBF_SXP2_YP1_Z_E)
  END IF
ELSE
  SELECT CASE(CLSOL)
  CASE('RICHA')     ! Richardson's method
    CALL RICHARDSON (ZLBCX, ZLBCY, XDXX, XDYY, XDZX, XDZY, XDZZ, &
                     PRHODJ, ZTHETAV, XDXHATM, XDYHATM,          &
                     XRHOM_E, XAF_E, XBFY_E, XCF_E, XTRIGSX, XTRIGSY,    &
                     NIFAXX, NIFAXY, NLAPITR_ELEC, KTCOUNT,      &
                     PRELAX, ZDV_SOURCE, ZPHIT)
! 
  CASE('CGRAD')     ! Conjugate Gradient method
    CALL CONJGRAD (ZLBCX, ZLBCY, XDXX, XDYY, XDZX, XDZY, XDZZ, &
                   PRHODJ,ZTHETAV, XDXHATM, XDYHATM,           &
                   XRHOM_E, XAF_E, XBFY_E, XCF_E, XTRIGSX, XTRIGSY,    &
                   NIFAXX, NIFAXY, NLAPITR_ELEC,               &
                   ZDV_SOURCE, ZPHIT)
  CASE('CRESI')     ! Conjugate Residual method
    CALL CONRESOL (ZLBCX, ZLBCY, XDXX, XDYY, XDZX, XDZY, XDZZ, &
                   PRHODJ,ZTHETAV, XDXHATM, XDYHATM,           &
                   XRHOM_E, XAF_E, XBFY_E, XCF_E, XTRIGSX, XTRIGSY,    &
                   NIFAXX, NIFAXY, NLAPITR_ELEC,               &
                   ZDV_SOURCE, ZPHIT)
  CASE('ZRESI')     ! Conjugate Residual method
    CALL CONRESOLZ (ZLBCX, ZLBCY, XDXX, XDYY, XDZX, XDZY, XDZZ,&
                   PRHODJ,ZTHETAV, XDXHATM, XDYHATM,           &
                   XRHOM_E, XAF_E, XBFY_E, XCF_E, XTRIGSX, XTRIGSY,    &
                   NIFAXX, NIFAXY, NLAPITR_ELEC,               &
                   ZDV_SOURCE, ZPHIT,                          &
                   XBFB_E, XBF_SXP2_YP1_Z_E)
  END SELECT
END IF
!
!
!*       3.2    compute the electric field
!               --------------------------
!              
!  Dirichlet Boundary condition for the electrical potential
ZPHIT(:,:,IKB-1) = -ZPHIT(:,:,IKB)
!
! E =  rhodj * Nabla V  
!
CALL ADD3DFIELD_ll(TZFIELDS_ll,ZPHIT)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
PEFIELDU(:,:,:) = PRHODJ(:,:,:) * GX_M_M(1,IKU,1,ZPHIT,XDXX,XDZZ,XDZX)
PEFIELDV(:,:,:) = PRHODJ(:,:,:) * GY_M_M(1,IKU,1,ZPHIT,XDYY,XDZZ,XDZY)
PEFIELDW(:,:,:) = PRHODJ(:,:,:) * GZ_M_M(1,IKU,1,ZPHIT,XDZZ) 
!
IF (PRESENT(PPHIT)) PPHIT(:,:,:) = - PRHODJ(:,:,:) * ZPHIT(:,:,:) 
!
DEALLOCATE(ZDV_SOURCE)
DEALLOCATE(ZPHIT)
DEALLOCATE(ZTHETAV)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ELEC_FIELD_n
