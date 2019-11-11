!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home//MESONH/MNH-V4-6-5/src/SRC_CHIMAQ/ch_solve_ph.f90
!-----------------------------------------------------------------
!!    #######################
      MODULE MODI_CH_SOLVE_PH
!!    #######################
!!
!
INTERFACE
!!
SUBROUTINE CH_SOLVE_PH(KLUOUT,PCONC,PTEMP,PLW,KLW,PPH,KRR)
IMPLICIT NONE
INTEGER, INTENT(IN)                      :: KLUOUT
INTEGER, INTENT(IN)                      :: KLW
REAL,    INTENT(IN),    DIMENSION(:,:)   :: PCONC !molec/cm3
REAL,    INTENT(IN),    DIMENSION(:)     :: PTEMP
REAL,    INTENT(IN),    DIMENSION(:)     :: PLW
REAL,    INTENT(INOUT), DIMENSION(:)     :: PPH
INTEGER, INTENT(IN)                      :: KRR
END SUBROUTINE CH_SOLVE_PH
!!
END INTERFACE
!!
END MODULE MODI_CH_SOLVE_PH
!!
!!    ##########################################################
      SUBROUTINE CH_SOLVE_PH(KLUOUT,PCONC,PTEMP,PLW,KLW,PPH,KRR)
!!    ##########################################################
!!
!!*** *CH_SOLVE_PH  calculate pH value
!!
!!    PURPOSE
!!    -------
!!       The purpose of this subroutine is to calculate the pH value in cloud
!!    water or in rainwater if requested by the user.
!!
!!    METHOD
!!    ------
!!       pH value is obtainded by solving the electoneutrality equation :
!!       [H+]+[NH4+]=[OH-]+[HCO3-]+2[CO3--]+[HSO3-]+2[SO3--]
!!                         +[NO3-]+2[SO4--]+[HCOO-]+SOM(IONS)
!!       i.e search the root of a degree 8 polynomial using Laguerre's method.
!!       Dissociation hypothesis : 
!!            Strong acids: HNO3 = NO3- and H2SO4 = SO4--
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    M. Leriche    *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 04/06/07
!!    J.-P. Pinty 11/07/07 extension to the rain drops
!!    M. Leriche 16/11/07 add sulfuric acid
!!    J.-P. Pinty 11/07/07 add CO3-- and SO3--
!!    M. Leriche 05/06/08 add sum of ions
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODI_CH_PH_POLYROOT
!
USE MODD_CH_M9_n,   ONLY : CNAMES, NEQ, NEQAQ
USE MODD_CST
USE MODD_CH_MNHC_n, ONLY : XCH_PHINIT
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)                      :: KLUOUT
INTEGER, INTENT(IN)                      :: KLW
REAL,    INTENT(IN),    DIMENSION(:,:)   :: PCONC !molec/cm3
REAL,    INTENT(IN),    DIMENSION(:)     :: PTEMP
REAL,    INTENT(IN),    DIMENSION(:)     :: PLW
REAL,    INTENT(INOUT), DIMENSION(:)     :: PPH
INTEGER, INTENT(IN)                      :: KRR
!
!*      0.2    declarations of local variables
!
REAL (KIND(0.0D0)), DIMENSION(KLW) :: K0, KW, K11, K12, K21, K22, K3
                                       !equilibrium constant
REAL (KIND(0.0D0)), DIMENSION(KLW) :: C0, C1, C2, C3, C4, SOM
                                       !chemical species concentrations in M
LOGICAL            :: GPH, GPH_TOT ! test if pH roots are correct
INTEGER            :: ITRUE ! counter
INTEGER            :: JI, JJ
!
REAL,                  DIMENSION(KLW) :: ZFACT
REAL,                  DIMENSION(KLW) :: ZPHINIT ! initial pH value
REAL (KIND(0.0D0)),    DIMENSION(KLW) :: KW2
REAL (KIND(0.0D0)),    DIMENSION(:,:), ALLOCATABLE :: ZCOEFS 
COMPLEX (KIND(0.0D0)), DIMENSION(:,:), ALLOCATABLE :: ZZCOEFS
INTEGER            :: IORDER  ! polynomial order 
LOGICAL            :: GPOLISH ! refine roots finding
COMPLEX (KIND(0.0D0)), DIMENSION(:,:), ALLOCATABLE :: ZZROOTS
REAL               :: ZDIST_PHMIN, ZDIST_PH, ZPH, ZROOT
!
!-------------------------------------------------------------------------------
!
!*       1.     INITIALIZATION
!               --------------
!
ZFACT(:) = 3.3556E-03  -  1./PTEMP(:)
!ACID = 0.  ! acidity term
KW = 1.0e-14*EXP(6716.0*ZFACT(:))   !ionic water product
K0 = 1.7E-05*EXP(4350.0*ZFACT(:))   !NH3 + H2O <-> NH4+ + OH-
K11 = 4.3E-07*EXP(920.0*ZFACT(:))   !CO2 + H2O <-> HCO3- + H+
K12 = 4.7E-11*EXP(1780.0*ZFACT(:))  !HCO3- <-> CO3(2-) + H+
K21 = 1.3E-02*EXP(-1965.0*ZFACT(:)) !SO2 + H2O <-> HSO3- + H+
K22 = 6.4E-08*EXP(-1432.0*ZFACT(:)) !HSO3- <-> SO3(2-) + H+
K3  = 1.8E-04*EXP(150.0*ZFACT(:))   !HCOOH <-> HCOO- + H+
KW2 = KW*KW
C0 = 0.                             !NH3
C1 = 0.                             !CO2
C2 = 0.                             !SO2
C3 = 0.                             !HCOOH = ORA1
C4 = 0.                             !HNO3 + 2 x H2SO4 + HCL = strong acid
SOM = 0.
IORDER = 8 !polynomial order
ALLOCATE(ZCOEFS(KLW,IORDER+1))
!
!-------------------------------------------------------------------------------
!
!*       2.     Find chemical species concentrations
!               ------------------------------------
!
ZFACT(:) = 1.e-3*XAVOGADRO*PLW(:)
SELECT CASE (KRR)
  CASE(2) ! cloud water
    DO JI = 1, NEQAQ/2
      JJ = NEQ-NEQAQ+JI
      IF (TRIM(CNAMES(JJ))=='WC_NH3') C0(:) = PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WC_CO2') C1(:) = PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WC_SO2') C2(:) = PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WC_ORA1') C3(:)= PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WC_HNO3') C4(:)= C4(:)+PCONC(:,JI)/(ZFACT(:))
      IF ((TRIM(CNAMES(JJ))=='WC_SULF') .OR. (TRIM(CNAMES(JJ))=='WC_H2SO4')) &
          C4(:)= C4(:)+2.*PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WC_HCL') C4(:)= C4(:)+PCONC(:,JI)/(ZFACT(:))
      IF (CNAMES(JJ)(1:4)=='WC_A') SOM(:) = SOM(:) + PCONC(:,JI)/(ZFACT(:))
      IF (CNAMES(JJ)(1:4)=='WC_B') SOM(:) = SOM(:) + 2.*PCONC(:,JI)/(ZFACT(:))
    END DO
  CASE(3) ! rain water
    DO JI = 1, NEQAQ/2
      JJ = NEQ-NEQAQ/2+JI
      IF (TRIM(CNAMES(JJ))=='WR_NH3') C0(:) = PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WR_CO2') C1(:) = PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WR_SO2') C2(:) = PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WR_ORA1') C3(:)= PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WR_HNO3') C4(:)= C4(:)+PCONC(:,JI)/(ZFACT(:))
      IF ((TRIM(CNAMES(JJ))=='WR_SULF') .OR. (TRIM(CNAMES(JJ))=='WR_H2SO4')) &
          C4(:)= C4(:)+2.*PCONC(:,JI)/(ZFACT(:))
      IF (TRIM(CNAMES(JJ))=='WR_HCL') C4(:)= C4(:)+PCONC(:,JI)/(ZFACT(:))
      IF (CNAMES(JJ)(1:4)=='WR_A') SOM(:) = SOM(:) + PCONC(:,JI)/(ZFACT(:))
      IF (CNAMES(JJ)(1:4)=='WR_B') SOM(:) = SOM(:) + 2.*PCONC(:,JI)/(ZFACT(:))
    END DO
  CASE DEFAULT
    WRITE(UNIT=KLUOUT,FMT='("CH_SOLVE_PH: bad setting for KRR")')
    RETURN
END SELECT
!
!-------------------------------------------------------------------------------
!
!*       3.     Compute polynomial coefficients
!               -------------------------------
!
ZCOEFS(:,9) = K0
ZCOEFS(:,8) = K0*(C0+K11+K21+K3-C4-SOM)+KW
ZCOEFS(:,7) = K0*(C0*(K11+K21+K3)+K11*(K21+K3-C4-SOM-C1+K12)+K21*(K3-C4-SOM-C2+K22)+K3*(-C4-SOM-C3)-KW)+KW*(K11+K21+K3-C4-SOM)
ZCOEFS(:,6) = K0*(C0*(K11*(K21+K3+K12)+K21*(K3+K22))+K11*(K21*K3+K21*(-C4-SOM-C1+K12-C2+K22)+K3*(-C4-SOM-C1+K12-C3) &
-C4*K12-SOM*K12-KW)+K21*(K3*(-C4-SOM-C2+K22-C3)-C4*K22-SOM*K22-KW)-K3*KW)+KW*(K11*(K21+K3-C4-SOM-C1+K12)+K21*(K3-C4-SOM-C2+K22) &
+K3*(-C4-SOM-C3)-KW)
ZCOEFS(:,5) = K0*(C0*(K11*(K21*(K3+K12+K22)+K3*K12)+K21*K3*K22)+K11*(K21*(K3*(-C4-SOM-C1+K12-C2+K22-C3) &
-C4*K12-C4*K22-SOM*K12-SOM*K22-KW-C1*K22-K12*C2+K12*K22)+K3*(-C4*K12-SOM*K12-KW-K12*C3)-KW*K12) &
+K21*(K3*(-C4*K22-SOM*K22-KW-K22*C3)-KW*K22))+KW*(K11*(K21*(K3-C4-SOM-C1+K12-C2+K22)+K3*(-C4-SOM-C1+K12-C3) &
-C4*K12-SOM*K12-KW)+K21*(K3*(-C4-SOM-C2+K22-C3)-C4*K22-SOM*K22-KW)-K3*KW)
ZCOEFS(:,4) = K0*(C0*(K11*(K21*(K3*(K12+K22)+K12*K22)))+K11*(K21*(K3*(-(K12+K22)*(C4+SOM)-KW-C1*K22-K12*C2+K12*K22-K12*C3-K22*C3) &
-C4*K12*K22-SOM*K12*K22-KW*K12-KW*K22)-K3*KW*K12)-K21*K3*KW*K22)+KW*(K11*(K21*(K3*(-C4-SOM-C1+K12-C2+K22-C3)-(K12+K22)*(C4+SOM) &
-KW-C1*K22-K12*C2+K12*K22)+K3*(-C4*K12-SOM*K12-KW-K12*C3)-KW*K12)+K21*(K3*(-C4*K22-SOM*K22-KW-K22*C3)-KW*K22))
ZCOEFS(:,3) = K0*(C0*K11*K21*K3*K12*K22+K11*(K21*(K3*(-C4*K12*K22-SOM*K12*K22-KW*(K12+K22)-K12*K22*C3)-KW*K12*K22))) &
+KW*(K11*(K21*(K3*(-(K12+K22)*(C4+SOM)-KW-C1*K22+K12*(-C2+K22-C3)-K22*C3)-C4*K12*K22-SOM*K12*K22-KW*(K12+K22))-K3*KW*K12) &
-K21*K3*KW*K22)
ZCOEFS(:,2) = +KW*(-K0*K11*K21*K3*K12*K22+K11*(K21*(K3*(-C4*K12*K22-SOM*K12*K22-KW*K12-KW*K22-K12*K22*C3)-KW*K12*K22)))
ZCOEFS(:,1) = -K11*K21*K3*KW2*K12*K22
!
! include additional terms due to the double acidity
! of an aqueous solution of SO2 and CO2
!
ZCOEFS(:,6) = ZCOEFS(:,6) -2.0*(K0*(K11*C1*K12+K21*C2*K22))
ZCOEFS(:,5) = ZCOEFS(:,5) -2.0*(K0*(K11*(K21*C1*K12+K21*C2*K22+K3*C1*K12)+K21*K3*C2*K22)+KW*(K11*C1*K12+K21*C2*K22))
ZCOEFS(:,4) = ZCOEFS(:,4) -2.0*(K0*(K11*K21*(K3*(C1*K12+C2*K22)+C1*K22*K12)+K11*K21*K12*C2*K22)+KW*(K11*K21*(C1*K12+C2*K22) &
+K3*(K11*C1*K12+K21*C2*K22)))
ZCOEFS(:,3) = ZCOEFS(:,3) -2.0*(K11*K21*(K3*(K0*(C1*K22*K12+K12*C2*K22)+KW*(C1*K12+C2*K22))+KW*(C1*K22*K12+K12*C2*K22)))
ZCOEFS(:,2) = ZCOEFS(:,2) -2.0*(K11*K21*K3*KW*K22*(C1*K12+K12*C2))
!
!-------------------------------------------------------------------------------
!
!*       4.     Search the roots and find pH value
!               ----------------------------------
!
ALLOCATE(ZZCOEFS(KLW,IORDER+1))
ALLOCATE(ZZROOTS(KLW,IORDER))
DO JJ=1,IORDER+1
  ZZCOEFS(:,JJ) = CMPLX(ZCOEFS(:,JJ),0.0)
END DO
GPOLISH=.TRUE.
!
DO JI = 1, KLW
  CALL CH_PH_POLYROOT(ZZCOEFS(JI,:), IORDER, GPOLISH, ZZROOTS(JI,:))
END DO
!
ZPHINIT(:)=PPH(:)
WHERE(ZPHINIT<=0.) ZPHINIT=XCH_PHINIT
PPH(:)  = XCH_PHINIT
GPH_TOT = .true.
ITRUE   = 0
DO JI = 1, KLW
  GPH = .false.
  ZDIST_PHMIN = 100. ! initialize an extremely high value of the pH !
  DO JJ=1,IORDER
    ZROOT = REAL(ZZROOTS(JI,JJ))
    IF( AIMAG(ZZROOTS(JI,JJ))<=1.e-20 .AND. (ZROOT>1.e-12.AND.ZROOT<1.e0) )THEN
      ZPH = - ALOG10(ZROOT)
      ZDIST_PH = ABS(ZPHINIT(JI) - ZPH)
!      ZDIST_PH = ABS(XCH_PHINIT - ZPH)
      IF( ZDIST_PH < ZDIST_PHMIN ) THEN ! find the closest root to initial pH value
        PPH(JI) = ZPH
        ZDIST_PHMIN = ZDIST_PH
        IF( .not.GPH ) ITRUE = ITRUE + 1
        GPH = .true.
      ENDIF
    ENDIF
  END DO
  GPH_TOT = GPH_TOT.AND.GPH
END DO
!
IF( .not.GPH_TOT ) THEN
  WRITE(UNIT=KLUOUT,FMT='("CH_SOLVE_PH: no convergence in the range ", &
                        & "0<pH<12,  Nunber of case =",F6.2," %")')    &
                           100.0*( 1.0-(FLOAT(ITRUE)/FLOAT(KLW)) )
ENDIF
!
DEALLOCATE(ZCOEFS)
DEALLOCATE(ZZCOEFS)
DEALLOCATE(ZZROOTS)
!
END SUBROUTINE CH_SOLVE_PH
