!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!###################
MODULE MODI_PRESSUREZ
!###################
!
INTERFACE
!
      SUBROUTINE PRESSUREZ(                                                &
      HLBCX,HLBCY,HPRESOPT,KITR,OITRADJ,KTCOUNT,PRELAX,KMI,                &
      PRHODJ,PDXX,PDYY,PDZZ,PDZX,PDZY,PDXHATM,PDYHATM,PRHOT,               &
      PAF,PBF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                           &
      KRR,KRRL,KRRI,PDRYMASST,PREFMASS,PMASS_O_PHI0,                       &
      PTHT,PRT,PRHODREF,PTHVREF,PRVREF,PEXNREF,PLINMASS,                   &
      PRUS,PRVS,PRWS,PPABST,                                               &
      PBFB,                                                                &
      PBF_SXP2_YP1_Z,                                                      &
      PRESIDUAL                                           ) !JUAN Z_SPLITING
!
IMPLICIT NONE
!
CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: HLBCX    ! x-direction LBC type
CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: HLBCY    ! y-direction LBC type
!
CHARACTER (LEN=5), INTENT(IN) :: HPRESOPT        ! choice of the pressure solver
!
INTEGER, INTENT(INOUT) :: KITR                   ! number of iterations for the
                                                 ! pressure solver
LOGICAL, INTENT(IN) :: OITRADJ                   ! switch to adjust or not KITR
INTEGER, INTENT(IN) :: KTCOUNT                   ! counter value of the
                                                 ! model temporal loop
INTEGER, INTENT(IN) :: KMI                       ! Model index
REAL, INTENT(IN)    :: PRELAX                    ! relaxation coefficient for
                                                 ! the Richardson's method
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ     ! density of reference state
                                                 ! * J
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX,PDYY,PDZZ,PDZX,PDZY ! metric coefficients
!
REAL, INTENT(IN) :: PDXHATM                     ! mean grid increment in the x
                                                ! direction
REAL, INTENT(IN) :: PDYHATM                     ! mean grid increment in the y
                                                ! direction
!
REAL, DIMENSION (:), INTENT(IN) :: PRHOT         !  mean of XRHODJ on the plane x y
                                                 !  localized at a mass level
!
REAL, DIMENSION(:), INTENT(IN)     :: PAF,PCF    ! vectors giving the nonvanishing
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF        ! elements of the tri-diag. y-slide
                                                 ! matrix in the pressure eq.
!
                                                 ! arrays of sin or cos values
                                                 ! for the FFT :
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSX        ! - along x
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSY        ! - along y
!
                                                 ! decomposition in prime
                                                 ! numbers for the FFT:
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXX      ! - along x
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXY      ! - along y
!
!
INTEGER,                  INTENT(IN)    :: KRR   ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL  ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI  ! Number of ice water var.
!
REAL,                     INTENT(IN)    :: PDRYMASST   ! Mass of dry air and of
REAL,                     INTENT(IN)    :: PREFMASS    ! the ref. atmosphere
REAL,                     INTENT(IN)    :: PMASS_O_PHI0 !    Mass / Phi0
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT       ! Temperature and water
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT        !  variables at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF    ! dry Density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF     ! Virtual Temperature
                                                       ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVREF      ! mixing ratio of the
                                                       ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF     ! Exner function
                                                       ! of the reference state
REAL,                     INTENT(IN)    :: PLINMASS    ! lineic mass through
                                                       ! open boundaries
!
REAL,       INTENT(INOUT) :: PRUS(:,:,:)         ! source term along x
REAL,       INTENT(INOUT) :: PRVS(:,:,:)         ! source term along y
REAL,       INTENT(INOUT) :: PRWS(:,:,:)         ! source term along z
!
REAL,       INTENT(INOUT)   :: PPABST(:,:,:)        ! pressure(t)
!
!JUAN Z_SPLITING
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBFB       ! elements of the tri-diag b-slide .
                                                 ! matrix in the pressure eq.
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF_SXP2_YP1_Z       ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
REAL, OPTIONAL                     :: PRESIDUAL
!JUAN Z_SPLITING
END SUBROUTINE PRESSUREZ
!
END INTERFACE
!
END MODULE MODI_PRESSUREZ
!     ######################################################################
      SUBROUTINE PRESSUREZ(                                                &
      HLBCX,HLBCY,HPRESOPT,KITR,OITRADJ,KTCOUNT,PRELAX,KMI,                &
      PRHODJ,PDXX,PDYY,PDZZ,PDZX,PDZY,PDXHATM,PDYHATM,PRHOT,               &
      PAF,PBF,PCF,PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,                           &
      KRR,KRRL,KRRI,PDRYMASST,PREFMASS,PMASS_O_PHI0,                       &
      PTHT,PRT,PRHODREF,PTHVREF,PRVREF,PEXNREF,PLINMASS,                   &
      PRUS,PRVS,PRWS,PPABST,                                               &
      PBFB,                                                                &
      PBF_SXP2_YP1_Z,                                                      &
      PRESIDUAL                                           ) !JUAN Z_SPLITING
!     ######################################################################
!
!!****  *PRESSUREZ * - solve the pressure equation and add the pressure term
!!      to the sources
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to solve the pressure equation:
!     with either the conjugate gradient method or the Richardson's method.
!     The pressure gradient is added to the sources in order
!     to nullify the divergence of the momentum* Thetavref*(1+Rvref)
!     at the time t+dt.
!
!!**  METHOD
!!    ------
!!     The divergence of the sources  ( RHS of the pressure equation ) is
!!    computed. The pressure equation is then solved by either CG method,
!!    either Richardson's method, or an exact method. Finally, the pressure
!!    gradient is added to the sources RUS, RVS, RWS.
!!    Finally, the absolute pressure is diagnozed from the total mass
!!    included in the simulation domain.
!!
!!    EXTERNAL
!!    --------
!!      Subroutine MASS_LEAK : assures global non-divergence condition in the
!!                             case of open boundaries
!!      Subroutine FLAT_INV  : solve the pressure equation for the case
!!                             without orography
!!      Subroutine RICHARDSON: solve the pressure equation with the
!!                             Richardson's method
!!      Subroutine CONJGRAD  : solve the pressure equation with the Conjugate
!!                             Gradient algorithm
!!      Function   GX_M_U : compute the gradient along x
!!      Function   GY_M_V : compute the gradient along y
!!      Function   GZ_M_W : compute the gradient along z
!!      Subroutine GDIV     : compute J times the divergence of 1/J times a vector
!!      Function MXM: compute an average in the x direction for a variable
!!      at a mass localization
!!      Function MYM: compute an average in the y direction for a variable
!!      at a mass localization
!!      Function MZM: compute an average in the z direction for a variable
!!      at a mass localization
!!      Subroutine P_ABS   : compute the constant for PABS and therefore, the
!!      absolute pressure function
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONF: model configuration
!!        LFLAT: logical switch for zero orography
!!        L2D  : logical switch for two-dimensional configuration
!!        LCARTESIAN : logical switch for cartesian geometry
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT, JPVEXT: define the number of marginal points out of the
!!        physical domain along horizontal and vertical directions respectively
!!      Module MODD_CST: physical constants
!!        XCPD
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (subroutine PRESSURE) + Book1 (  )
!!
!!    AUTHOR
!!    ------
!!      P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      05/07/94
!!      Modification  03/01/95  (Lafore)  To add the absolute pressure diagnosis
!!      Modification  31/01/95  (Stein)   Copy of the pressure function in the
!!                                        2D case in the two outermost planes
!!      Modification  16/02/95  (Mallet)  Add the call to MASS_LEAK
!!      Modification  16/03/95  (Stein)  change the argument list of the
!!                              gradient and remove R from the historical var.
!!      Modification  30/06/95  (Stein)  Add a test not to compute the absolute
!!                              pressure in the Boussinesq case
!!                    16/10/95 (J. Stein) change the budget calls
!!                    29/01/96 (J. Stein) call iterative resolution for
!!                              non-cartessian geometry
!!                    19/12/96 (J.-P. Pinty) update the budget calls
!!                    14/01/97 (Stein,Lafore) New anelastic equations
!!                    17/12/97 ( Stein )include the case of non-vanishing
!!                              orography at the lbc
!!                    26/03/98 (Stein,Jabouille) fix the value of the corner point
!!                    15/06/98  (D.Lugato, R.Guivarch) Parallelisation
!!                    25/08/99 (J.-P. Pinty) add CRESI option to CPRESOPT
!!                    06/11/02 (V. Masson) update the budget calls
!!                    24/08/2005 (J. escobar) BUG : remove IIE+1, IJE+1 out of bound
!!                                references in parallel run
!!                    08/2010 (V.Masson, C.Lac) Add UPDATE_HALO
!!                    11/2010 (V.Masson, C.Lac) PPABST, must not be cyclic => add temp array
!!                                             to save it before UPDATE_HALO
!!                    07/2011 (J.escobar ) Bypass Bug with ifort11/12 on  HLBCX,HLBCY 
!!                    09/2001 (J.escobar ) reintroduce correctly the GMAXLOC_ll call 
!!                    11/2010 (V.Masson, C.Lac) PPABST must not be cyclic => add temp array
!!                                             to save it before UPDATE_HALO
!!                    02/2013 (J.Escobar ) add a test on abs(err) > 100.O for BG without controle of NAN
!!                    2012    (V.Masson)  Modif update_halo due to CONTRAV
!!                    2014    (C.Lac) correction for 3D run with LBOUSS=.TRUE.
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!   J.escobar : check nb proc versus ZRESI & min(DIMX,DIMY)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  Philippe Wautelet: 22/01/2019: use standard FLUSH statement instead of non standard intrinsics
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ARGSLIST_ll, ONLY: LIST_ll
USE MODD_BUDGET
USE MODD_CST
USE MODD_CONF
USE MODD_DYN_n,       ONLY: LRES, XRES
USE MODD_LUNIT_n,     ONLY: TLUOUT
USE MODD_MPIF
USE MODD_PARAMETERS
USE MODD_REF,         ONLY: LBOUSS
USE MODD_VAR_ll,      ONLY: MPI_PRECISION, NMNH_COMM_WORLD , NPROC
!
USE MODE_IO_ll,       ONLY: CLOSE_ll
USE MODE_ll
USE MODE_MPPDB
USE MODE_MSG
!
USE MODI_BUDGET
USE MODI_CONJGRAD
USE MODI_CONRESOL
USE MODI_CONRESOLZ
USE MODI_FLAT_INV
USE MODI_FLAT_INVZ
USE MODI_GDIV
USE MODI_GRADIENT_M
USE MODI_MASS_LEAK
USE MODI_P_ABS
USE MODI_RICHARDSON
USE MODI_SHUMAN
USE MODI_SUM_ll , ONLY : GMAXLOC_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: HLBCX    ! x-direction LBC type
  CHARACTER (LEN=*), DIMENSION(:), INTENT(IN) :: HLBCY    ! y-direction LBC type
!
CHARACTER (LEN=5), INTENT(IN) :: HPRESOPT        ! choice of the pressure solver
!
INTEGER, INTENT(INOUT) :: KITR                   ! number of iterations for the
                                                 ! pressure solver
LOGICAL, INTENT(IN) :: OITRADJ                   ! switch to adjust or not KITR
INTEGER, INTENT(IN) :: KTCOUNT                   ! counter value of the
                                                 ! model temporal loop
INTEGER, INTENT(IN) :: KMI                       ! Model index
REAL, INTENT(IN)    :: PRELAX                    ! relaxation coefficient for
                                                 ! the Richardson's method
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ     ! density of reference state
                                                 ! * J
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX,PDYY,PDZZ,PDZX,PDZY ! metric coefficients
!
REAL, INTENT(IN) :: PDXHATM                     ! mean grid increment in the x
                                                ! direction
REAL, INTENT(IN) :: PDYHATM                     ! mean grid increment in the y
                                                ! direction
!
REAL, DIMENSION (:), INTENT(IN) :: PRHOT         !  mean of XRHODJ on the plane x y
                                                 !  localized at a mass level
!
REAL, DIMENSION(:), INTENT(IN)     :: PAF,PCF    ! vectors giving the nonvanishing
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF        ! elements of the tri-diag y-slide .
                                                 ! matrix in the pressure eq.
!
                                                 ! arrays of sin or cos values
                                                 ! for the FFT :
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSX        ! - along x
REAL, DIMENSION(:), INTENT(IN) :: PTRIGSY        ! - along y
!
                                                 ! decomposition in prime
                                                 ! numbers for the FFT:
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXX      ! - along x
INTEGER, DIMENSION(19), INTENT(IN) :: KIFAXY      ! - along y
!
INTEGER,                  INTENT(IN)    :: KRR   ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL  ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI  ! Number of ice water var.
!
REAL,                     INTENT(IN)    :: PDRYMASST   ! Mass of dry air and of
REAL,                     INTENT(IN)    :: PREFMASS    ! the ref. atmosphere
REAL,                     INTENT(IN)    :: PMASS_O_PHI0 !    Mass / Phi0
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT       ! Temperature and water
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT        !  variables at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF    ! dry Density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF     ! Virtual Temperature
                                                       ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVREF      ! mixing ratio of the
                                                       ! reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF     ! Exner function
                                                       ! of the reference state
REAL,                     INTENT(IN)    :: PLINMASS    ! lineic mass through
                                                       ! open boundaries
!
REAL,       INTENT(INOUT) :: PRUS(:,:,:)         ! source term along x
REAL,       INTENT(INOUT) :: PRVS(:,:,:)         ! source term along y
REAL,       INTENT(INOUT) :: PRWS(:,:,:)         ! source term along z
!
REAL,       INTENT(INOUT)   :: PPABST(:,:,:)        ! pressure(t)
!
!JUAN Z_SPLITING
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBFB       ! elements of the tri-diag b-slide .
                                                 ! matrix in the pressure eq.
REAL, DIMENSION(:,:,:), INTENT(IN) :: PBF_SXP2_YP1_Z       ! elements of the tri-diag. SXP2_YP1_Z-slide 
                                                   ! matrix in the pressure eq.
REAL, OPTIONAL                     :: PRESIDUAL
!JUAN Z_SPLITING
!
!*       0.2   declarations of local variables
!
!                                                           Metric coefficients:
!
REAL, DIMENSION(SIZE(PPABST,1),SIZE(PPABST,2),SIZE(PPABST,3)) :: ZDV_SOURCE
!                                                   ! divergence of the sources
!
INTEGER :: IIB          ! indice I for the first inner mass point along x
INTEGER :: IIE          ! indice I for the last inner mass point along x
INTEGER :: IJB          ! indice J for the first inner mass point along y
INTEGER :: IJE          ! indice J for the last inner mass point along y
INTEGER :: IKB          ! indice K for the first inner mass point along z
INTEGER :: IKE          ! indice K for the last inner mass point along z
INTEGER :: ILUOUT       ! Logical unit of output listing
INTEGER :: IRESP        ! Return code of FM routines
!
REAL, DIMENSION(SIZE(PPABST,1),SIZE(PPABST,2),SIZE(PPABST,3)) :: ZTHETAV, &
                        ! virtual potential temperature
                                                                 ZPHIT
                        ! MAE + DUR => Exner function perturbation
                        ! LHE       => Exner function perturbation * CPD * THVREF
!
REAL            :: ZRV_OV_RD !  XRV / XRD
REAL                  :: ZMAXVAL, ZMAXRES, ZMAX,ZMAX_ll ! for print
INTEGER, DIMENSION(3) :: IMAXLOC ! purpose
INTEGER         :: JWATER          ! loop index on water species
INTEGER         :: IIU,IJU,IKU     ! array sizes in I,J,K
INTEGER         :: JK              ! loop index on the vertical levels
INTEGER         :: JI,JJ
!
REAL, DIMENSION(SIZE(PDXX,1),SIZE(PDXX,3)) :: ZPABS_S ! local pressure on southern side
REAL, DIMENSION(SIZE(PDXX,1),SIZE(PDXX,3)) :: ZPABS_N ! local pressure on northern side
REAL, DIMENSION(SIZE(PDYY,2),SIZE(PDXX,3)) :: ZPABS_E ! local pressure on eastern side
REAL, DIMENSION(SIZE(PDYY,2),SIZE(PDXX,3)) :: ZPABS_W ! local pressure on western side
INTEGER :: IINFO_ll,KINFO
TYPE(LIST_ll), POINTER :: TZFIELDS_ll, TZFIELDS2_ll  ! list of fields to exchange
!
INTEGER :: IIB_I,IIE_I,IJB_I,IJE_I
INTEGER :: IIMAX_ll,IJMAX_ll
!
!
!------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       1.      PRELIMINARIES
!                -------------
!
ILUOUT = TLUOUT%NLU
!
CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
IF ( ( MIN(IIMAX_ll,IJMAX_ll) < NPROC  ) .AND. ( HPRESOPT /= 'ZRESI' ) ) THEN
   WRITE(UNIT=ILUOUT,FMT=*) 'ERROR IN PRESSUREZ:: YOU WANT TO USE TO MANY PROCESSOR WITHOUT CPRESOPT="ZRESI" '
   WRITE(UNIT=ILUOUT,FMT=*) 'MIN(IIMAX_ll,IJMAX_ll)=',MIN(IIMAX_ll,IJMAX_ll),' < NPROC =', NPROC
   WRITE(UNIT=ILUOUT,FMT=*) 'YOU HAVE TO SET CPRESOPT="ZRESI => JOB ABORTED '
   CALL PRINT_MSG(NVERB_FATAL,'GEN','PRESSUREZ','')
ENDIF
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
CALL GET_DIM_EXT_ll('B',IIU,IJU)
!
IKB= 1+JPVEXT
IKU= SIZE(PPABST,3)
IKE= IKU - JPVEXT
!
ZPABS_S(:,:) = 0.
ZPABS_N(:,:) = 0.
ZPABS_E(:,:) = 0.
ZPABS_W(:,:) = 0.
!
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTE THE LINEIC MASS
!              -----------------------
!
IF ( ANY(HLBCX(:)=='OPEN') .OR. ANY(HLBCY(:)=='OPEN') ) THEN                     
  CALL MASS_LEAK(PDXX,PDYY,HLBCX,HLBCY,PLINMASS,PRHODJ,PRUS,PRVS)
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.    COMPUTE THE FORCING TERM FOR THE PRESSURE EQUATION
!              --------------------------------------------------
!
!
CALL MPPDB_CHECK3D(PRUS,"pressurez 4-before update_halo_ll::PRUS",PRECISION)
CALL MPPDB_CHECK3D(PRVS,"pressurez 4-before update_halo_ll::PRVS",PRECISION)
CALL MPPDB_CHECK3D(PRWS,"pressurez 4-before update_halo_ll::PRWS",PRECISION)
NULLIFY(TZFIELDS_ll)
CALL ADD3DFIELD_ll(TZFIELDS_ll, PRUS)
CALL ADD3DFIELD_ll(TZFIELDS_ll, PRVS)
CALL ADD3DFIELD_ll(TZFIELDS_ll, PRWS)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
CALL MPPDB_CHECK3D(PRUS,"pressurez 4-after update_halo_ll::PRUS",PRECISION)
CALL MPPDB_CHECK3D(PRVS,"pressurez 4-after update_halo_ll::PRVS",PRECISION)
CALL MPPDB_CHECK3D(PRWS,"pressurez 4-after update_halo_ll::PRWS",PRECISION)
!
CALL GDIV(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRUS,PRVS,PRWS,ZDV_SOURCE)
!
! The non-homogenous Neuman problem is transformed in an homogenous Neuman
! problem in the non-periodic cases
IF (HLBCX(1) /= 'CYCL') THEN
  IF (LWEST_ll()) ZDV_SOURCE(IIB-1,:,:) = 0.
  IF (LEAST_ll()) ZDV_SOURCE(IIE+1,:,:) = 0.
ENDIF
!
IF (.NOT. L2D .AND. HLBCY(1) /= 'CYCL') THEN
  IF (LSOUTH_ll()) ZDV_SOURCE(:,IJB-1,:) = 0.
  IF (LNORTH_ll()) ZDV_SOURCE(:,IJE+1,:) = 0.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.    SOLVE THE PRESSURE EQUATION
!              ---------------------------
!
!
!*       5.1   Compute the virtual theta and the pressure perturbation
!              -------------------------------------------------------
!
IF(CEQNSYS=='MAE' .OR. CEQNSYS=='DUR') THEN
  IF(KRR > 0) THEN
  !
  !   compute the ratio : 1 + total water mass / dry air mass
    ZRV_OV_RD = XRV / XRD
    ZTHETAV(:,:,:) = 1. + PRT(:,:,:,1)
    DO JWATER = 2 , 1+KRRL+KRRI
      ZTHETAV(:,:,:) = ZTHETAV(:,:,:) + PRT(:,:,:,JWATER)
    END DO
  !   compute the virtual potential temperature when water is present in any
  !   form
    ZTHETAV(:,:,:) = PTHT(:,:,:) * (1. + PRT(:,:,:,1) * ZRV_OV_RD) / ZTHETAV(:,:,:)
  ELSE
  !   compute the virtual potential temperature when water is absent
    ZTHETAV(:,:,:) = PTHT(:,:,:)
  END IF
  !
  ZPHIT(:,:,:)=(PPABST(:,:,:)/XP00)**(XRD/XCPD)-PEXNREF(:,:,:)
  !
ELSEIF(CEQNSYS=='LHE') THEN
  ZPHIT(:,:,:)= ((PPABST(:,:,:)/XP00)**(XRD/XCPD)-PEXNREF(:,:,:))   &
               * XCPD * PTHVREF(:,:,:)
  !
END IF
!
IF(CEQNSYS=='LHE'.AND. LFLAT .AND. LCARTESIAN) THEN
   ! flat cartesian LHE case -> exact solution
 IF ( HPRESOPT /= "ZRESI" ) THEN
  CALL FLAT_INV(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOT,PAF,PBF,PCF,         &
               PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZDV_SOURCE,ZPHIT)
 ELSE
  CALL FLAT_INVZ(HLBCX,HLBCY,PDXHATM,PDYHATM,PRHOT,PAF,PBF,PCF,         &
               PTRIGSX,PTRIGSY,KIFAXX,KIFAXY,ZDV_SOURCE,ZPHIT,& 
               PBFB,&
               PBF_SXP2_YP1_Z)
 ENDIF
ELSE
  SELECT CASE(HPRESOPT)
  CASE('RICHA')     ! Richardson's method
!
    CALL RICHARDSON(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,ZTHETAV,      &
    PDXHATM,PDYHATM,PRHOT,PAF,PBF,PCF,PTRIGSX,PTRIGSY,                        &
    KIFAXX,KIFAXY,KITR,KTCOUNT,PRELAX,ZDV_SOURCE,ZPHIT)
!
   CASE('CGRAD')     ! Conjugate Gradient method
     CALL CONJGRAD(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,ZTHETAV,       &
     PDXHATM,PDYHATM,PRHOT,PAF,PBF,PCF,PTRIGSX,PTRIGSY,                       &
     KIFAXX,KIFAXY,KITR,ZDV_SOURCE,ZPHIT)
!
  CASE('CRESI')     ! Conjugate Residual method
    CALL CONRESOL(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,ZTHETAV,       &
    PDXHATM,PDYHATM,PRHOT,PAF,PBF,PCF,PTRIGSX,PTRIGSY,                       &
    KIFAXX,KIFAXY,KITR,ZDV_SOURCE,ZPHIT)
!
  CASE('ZRESI')     ! Conjugate Residual method
    CALL CONRESOLZ(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,ZTHETAV,       &
    PDXHATM,PDYHATM,PRHOT,PAF,PBF,PCF,PTRIGSX,PTRIGSY,                       &
    KIFAXX,KIFAXY,KITR,ZDV_SOURCE,ZPHIT,                                     &
    PBFB,&
    PBF_SXP2_YP1_Z) !JUAN Z_SPLITING
  END SELECT
END IF
!
!-------------------------------------------------------------------------------
!
!*       6.    ADD THE PRESSURE GRADIENT TO THE SOURCES
!              ----------------------------------------
!
IF ( HLBCX(1) /= 'CYCL' ) THEN
  IF(LWEST_ll()) ZPHIT(IIB-1,:,IKB-1) = ZPHIT(IIB,:,IKB-1)
  IF(LEAST_ll()) ZPHIT(IIE+1,:,IKB-1) = ZPHIT(IIE,:,IKB-1)
ENDIF
IF ( HLBCY(1) /= 'CYCL' ) THEN
  IF (LSOUTH_ll()) ZPHIT(:,IJB-1,IKB-1) = ZPHIT(:,IJB,IKB-1)
  IF (LNORTH_ll()) ZPHIT(:,IJE+1,IKB-1) = ZPHIT(:,IJE,IKB-1)
ENDIF
!
IF ( L2D ) THEN
  IF (LSOUTH_ll()) ZPHIT(:,IJB-1,:) = ZPHIT(:,IJB,:)
  IF (LNORTH_ll()) ZPHIT(:,IJE+1,:) = ZPHIT(:,IJB,:)
END IF
!
ZDV_SOURCE = GX_M_U(1,IKU,1,ZPHIT,PDXX,PDZZ,PDZX)
!
IF ( HLBCX(1) /= 'CYCL' ) THEN
  IF(LWEST_ll()) THEN
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
    DO JK=2,IKU-1
      ZDV_SOURCE(IIB,:,JK)=                                                    &
       (ZPHIT(IIB,:,JK) - ZPHIT(IIB-1,:,JK) - 0.5 * (                              &
        PDZX(IIB,:,JK)   * (ZPHIT(IIB,:,JK)-ZPHIT(IIB,:,JK-1)) / PDZZ(IIB,:,JK)      &
       +PDZX(IIB,:,JK+1) * (ZPHIT(IIB,:,JK+1)-ZPHIT(IIB,:,JK)) / PDZZ(IIB,:,JK+1)    &
                                              )                              &
       ) / PDXX(IIB,:,JK)
    END DO
  ENDIF
  !
  IF(LEAST_ll()) THEN
    DO JK=2,IKU-1
      ZDV_SOURCE(IIE+1,:,JK)=                                                   &
        (ZPHIT(IIE+1,:,JK) - ZPHIT(IIE+1-1,:,JK) - 0.5 * (                        &
         PDZX(IIE+1,:,JK)   * (ZPHIT(IIE+1-1,:,JK)-ZPHIT(IIE+1-1,:,JK-1))           &
                          / PDZZ(IIE+1-1,:,JK)                                  &
        +PDZX(IIE+1,:,JK+1) * (ZPHIT(IIE+1-1,:,JK+1)-ZPHIT(IIE+1-1,:,JK))           &
                          / PDZZ(IIE+1-1,:,JK+1)                                &
                                                     )                        &
        ) / PDXX(IIE+1,:,JK)
    END DO
  END IF
END IF
!
CALL MPPDB_CHECK3DM("before MXM PRESSUREZ :PRU/V/WS",PRECISION,PRUS,PRVS,PRWS)
IF(CEQNSYS=='MAE' .OR. CEQNSYS=='DUR') THEN
  PRUS = PRUS - MXM(PRHODJ * XCPD * ZTHETAV) * ZDV_SOURCE
  PRWS = PRWS - MZM(1,IKU,1,PRHODJ * XCPD * ZTHETAV) * GZ_M_W(1,IKU,1,ZPHIT,PDZZ)
ELSEIF(CEQNSYS=='LHE') THEN
  PRUS = PRUS - MXM(PRHODJ) * ZDV_SOURCE
  PRWS = PRWS - MZM(1,IKU,1,PRHODJ) * GZ_M_W(1,IKU,1,ZPHIT,PDZZ)
END IF
!
IF(.NOT. L2D) THEN
!
  ZDV_SOURCE = GY_M_V(1,IKU,1,ZPHIT,PDYY,PDZZ,PDZY)
!
  IF ( HLBCY(1) /= 'CYCL' ) THEN
    IF (LSOUTH_ll()) THEN
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
      DO JK=2,IKU-1
        ZDV_SOURCE(:,IJB,JK)=                                                  &
         (ZPHIT(:,IJB,JK) - ZPHIT(:,IJB-1,JK) - 0.5 * (                            &
          PDZY(:,IJB,JK)   * (ZPHIT(:,IJB,JK)-ZPHIT(:,IJB,JK-1)) / PDZZ(:,IJB,JK)    &
         +PDZY(:,IJB,JK+1) * (ZPHIT(:,IJB,JK+1)-ZPHIT(:,IJB,JK)) / PDZZ(:,IJB,JK+1)  &
                                                )                            &
         ) / PDYY(:,IJB,JK)
      END DO
    END IF
    !
    IF (LNORTH_ll()) THEN
      DO JK=2,IKU-1
        ZDV_SOURCE(:,IJE+1,JK)=                                                &
         (ZPHIT(:,IJE+1,JK) - ZPHIT(:,IJE+1-1,JK) - 0.5 * (                      &
          PDZY(:,IJE+1,JK)   * (ZPHIT(:,IJE+1-1,JK)-ZPHIT(:,IJE+1-1,JK-1))         &
                           / PDZZ(:,IJE+1-1,JK)                                &
         +PDZY(:,IJE+1,JK+1) * (ZPHIT(:,IJE+1-1,JK+1)-ZPHIT(:,IJE+1-1,JK))         &
                           / PDZZ(:,IJE+1-1,JK+1)                              &
                                                      )                      &
        ) / PDYY(:,IJE+1,JK)
      END DO
    END IF
  END IF
!
  CALL MPPDB_CHECK3DM("before MYM PRESSUREZ :PRU/V/WS",PRECISION,PRUS,PRVS,PRWS)
  IF(CEQNSYS=='MAE' .OR. CEQNSYS=='DUR') THEN
    PRVS = PRVS - MYM(PRHODJ * XCPD * ZTHETAV) * ZDV_SOURCE
  ELSEIF(CEQNSYS=='LHE') THEN
    PRVS = PRVS - MYM(PRHODJ) * ZDV_SOURCE
  END IF
END IF
!
!! same boundary conditions as in gdiv ... !! (provisory coding)
!! (necessary when NVERB=1)
!!
    PRUS(:,:,IKB-1)=PRUS(:,:,IKB)
    PRUS(:,:,IKE+1)=PRUS(:,:,IKE)
    PRVS(:,:,IKB-1)=PRVS(:,:,IKB)
    PRVS(:,:,IKE+1)=PRVS(:,:,IKE)
!
NULLIFY(TZFIELDS2_ll)
CALL ADD3DFIELD_ll(TZFIELDS2_ll, PRUS)
CALL ADD3DFIELD_ll(TZFIELDS2_ll, PRVS)
CALL ADD3DFIELD_ll(TZFIELDS2_ll, PRWS)
CALL UPDATE_HALO_ll(TZFIELDS2_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS2_ll)
!
!  compute the residual divergence
CALL GDIV(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRUS,PRVS,PRWS,ZDV_SOURCE)
!
IF ( CEQNSYS=='DUR' ) THEN
  IF ( SIZE(PRVREF,1) == 0 ) THEN
    ZDV_SOURCE=ZDV_SOURCE/PRHODJ/XTH00*PRHODREF*PTHVREF
  ELSE
    ZDV_SOURCE=ZDV_SOURCE/PRHODJ/XTH00*PRHODREF*PTHVREF*(1.+PRVREF)
  END IF
ELSEIF( CEQNSYS=='MAE' .OR. CEQNSYS=='LHE' ) THEN
  ZDV_SOURCE=ZDV_SOURCE/PRHODJ*PRHODREF
END IF
!
ZMAXVAL=MAX_ll(ABS(ZDV_SOURCE),IINFO_ll)
!JUANZ
IF (PRESENT(PRESIDUAL)) PRESIDUAL = ZMAXVAL
!JUANZ
IMAXLOC=GMAXLOC_ll( ABS(ZDV_SOURCE) )
!
WRITE(ILUOUT,*) 'residual divergence / 2 DT', ZMAXVAL,     &
                ' located at ',   IMAXLOC
FLUSH(unit=ILUOUT)
IF (ABS(ZMAXVAL) .GT. 100.0 ) THEN
   WRITE(ILUOUT,*) ' pressurez.f90 STOP :: SOMETHING WRONG WITH PRESSURE , ABS(RESIDUAL) > 100.0 '  
   FLUSH(unit=ILUOUT)
   STOP ' pressurez.f90 STOP :: SOMETHING WRONG WITH PRESSURE , ABS(RESIDUAL) > 100.0 '
ENDIF 
! number of iterations adjusted
IF (LRES) THEN
   ZMAXRES = XRES
ELSEIF (LFLAT .AND. LCARTESIAN) THEN
  ZMAXRES = XRES_FLAT_CART
ELSE
  ZMAXRES = XRES_OTHER
END IF
!
IF (OITRADJ) THEN
  IF (ZMAXVAL>10.*ZMAXRES) THEN
    KITR=KITR+2
    WRITE(ILUOUT,*) 'NITR adjusted to ', KITR
  ELSE IF (ZMAXVAL<ZMAXRES) THEN
    KITR=MAX(KITR-1,1)
    WRITE(ILUOUT,*) 'NITR adjusted to ', KITR
  ENDIF
ENDIF
!
!*       7.    STORAGE OF THE FIELDS IN BUDGET ARRAYS
!              --------------------------------------
!
IF (LBUDGET_U) CALL BUDGET (PRUS,1,'PRES_BU_RU')
IF (LBUDGET_V) CALL BUDGET (PRVS,2,'PRES_BU_RV')
IF (LBUDGET_W) CALL BUDGET (PRWS,3,'PRES_BU_RW')
!
!-------------------------------------------------------------------------------
!
!*       8.    ABSOLUTE PRESSURE COMPUTATION
!              -----------------------------
!
ZMAX = MAXVAL(ABS ( PRHODREF(:,:,IKB)-PRHODREF(:,:,IKE)) )
CALL MPI_ALLREDUCE(ZMAX, ZMAX_ll, 1, MPI_PRECISION, MPI_MAX,  &
                   NMNH_COMM_WORLD, KINFO)
!IF (      ABS(PRHODREF(IIB,IJB,IKB)-PRHODREF(IIB,IJB,IKE)) > 1.E-12 &
!  .AND. KTCOUNT >0 ) THEN
IF ((ZMAX_ll > 1.E-12) .AND. KTCOUNT >0 ) THEN 
!IF (  KTCOUNT >0 .AND. .NOT.LBOUSS ) THEN
  CALL P_ABS   ( KRR, KRRL, KRRI, PDRYMASST, PREFMASS, PMASS_O_PHI0, &
                 PTHT, PRT, PRHODJ, PRHODREF, ZTHETAV, PTHVREF,      &
                 PRVREF, PEXNREF,  ZPHIT                             )
!
  IF(CEQNSYS=='MAE' .OR. CEQNSYS=='DUR') THEN
    PPABST(:,:,:)=XP00*(ZPHIT+PEXNREF)**(XCPD/XRD)
  ELSEIF(CEQNSYS=='LHE') THEN
    PPABST(:,:,:)=XP00*(ZPHIT/(XCPD*PTHVREF)+PEXNREF)**(XCPD/XRD)
  ENDIF
!
  IF( HLBCX(1) == 'CYCL' ) THEN
   IF (LWEST_ll()) THEN
       ZPABS_W(:,:)= PPABST(IIB,:,:)
   END IF
!
   IF (LEAST_ll()) THEN
       ZPABS_E(:,:)= PPABST(IIE+1,:,:)
   END IF
!
  END IF
!
  IF( HLBCY(1) == 'CYCL' ) THEN
   IF (LSOUTH_ll()) THEN
      ZPABS_S(:,:)= PPABST(:,IJB,:)
   END IF
!
   IF (LNORTH_ll()) THEN
      ZPABS_N(:,:)= PPABST(:,IJE+1,:)
   END IF
!
  END IF
!
  CALL ADD3DFIELD_ll(TZFIELDS_ll, PPABST)
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
!
  IF( HLBCX(1) == 'CYCL' ) THEN
   IF (LWEST_ll()) THEN
       PPABST(IIB,:,:) = ZPABS_W(:,:)
   END IF
!
   IF (LEAST_ll()) THEN
       PPABST(IIE+1,:,:) = ZPABS_E(:,:)
   END IF
!
  END IF
!
  IF( HLBCY(1) == 'CYCL' ) THEN
   IF (LSOUTH_ll()) THEN
      PPABST(:,IJB,:) = ZPABS_S(:,:)
   END IF
!
   IF (LNORTH_ll()) THEN
      PPABST(:,IJE+1,:) = ZPABS_N(:,:)
   END IF
!
  END IF
!
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PRESSUREZ
