!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 solver 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODI_P_ABS
!     #################
!
INTERFACE
!
      SUBROUTINE P_ABS (KRR, KRRL, KRRI, PDRYMASST, PREFMASS, PMASS_O_PHI0, &
                        PTHT, PRT, PRHODJ, PRHODREF, PTHETAV, PTHVREF,      &
                        PRVREF, PEXNREF, PPHIT )
!  
IMPLICIT NONE
!
INTEGER,                  INTENT(IN)    :: KRR  ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL,                     INTENT(IN)    :: PDRYMASST   ! Mass of dry air and of
REAL,                     INTENT(IN)    :: PREFMASS    ! the ref. atmosphere
                                          !  contained in the simulation domain
REAL,                     INTENT(IN)    :: PMASS_O_PHI0 !    Mass / Phi0 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT        ! Temperature and water
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT         !  variables at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ      ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHETAV     ! virtual potential temp.
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF    ! dry Density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF     ! Virtual Temperature
                                                  ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVREF ! vapor mixing ratio 
                                       ! for the reference state 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF! Exner function of the
                                                  ! reference state
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PPHIT  ! Perturbation of
               ! either the Exner function Pi or Pi * Cpd * THvref
!
!
END SUBROUTINE P_ABS
!
END INTERFACE
!
END MODULE MODI_P_ABS
!     #######################################################################
      SUBROUTINE P_ABS (KRR, KRRL, KRRI, PDRYMASST, PREFMASS, PMASS_O_PHI0, &
		                PTHT, PRT, PRHODJ, PRHODREF, PTHETAV, PTHVREF,      &
                        PRVREF, PEXNREF, PPHIT )
!     #######################################################################
!
!!****  *P_ABS * - routine to compute the absolute Exner pressure deviation PHI
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the absolute Exner
!!      pressure Pi ( or Pi multiplied by Cpd*Thetavref) deviation PHI, 
!!      which is not determined for an anelatic system. 
!!      It also diagnozes the total mass of water Mw.
!!
!!     
!!**  METHOD
!!    ------
!!      The knowledge of the total mass of dry air Md and of water Mw 
!!    (including all water categories), allowed to diagnoze the absolute  
!!    Exner pressure PHI. The equation of state is not anymore linearized.
!!
!!    EXTERNAL
!!    --------
!!      none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST 
!!           XRD,XRV      Gaz constant for dry air Rd and wator vapor Rv
!!           XCPD         Specific heat at constant pressure for dry air Cp
!!           XP00         Reference pressure  
!!
!!      Module MODD_PARAMETERS : contains parameters commun to all models
!!        JPHEXT : Horizontal EXTernal points number (JPHEXT=1 for this version)
!!        JPVEXT : Vertical   EXTernal points number (JPVEXT=1 for this version)
!!      Module MODD_CONF  :
!!        CEQNSYS
!!
!!    REFERENCE
!!    ---------
!!      Book1 and book2 of documentation ( routine P_ABS )
!!
!!    AUTHOR
!!    ------
!!	J.-P. Lafore     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/12/94 
!!      J.P. Lafore 10/02/95   Bug correction in ZMASSGUESS
!!      J. Stein    16/03/95   Remove R from the historical variables
!!      J.P. Lafore 14/01/97   Introduction of 2 anelastic systems:
!!                              Modified Anelastic Equation and one derived 
!!                              from Durran (1989), MAE and DUR respectively
!!                  15/06/98  (D.Lugato, R.Guivarch) Parallelisation
!!      J. Colin       07/13  Add LBOUSS
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS 
!              ------------
!
USE MODD_CST
USE MODD_CONF
USE MODD_PARAMETERS
USE MODD_REF, ONLY : LBOUSS
!
USE MODE_ll
!JUAN
USE MODE_REPRO_SUM
!JUAN
!  
IMPLICIT NONE
!  
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KRR  ! Total number of water var.
INTEGER,                  INTENT(IN)    :: KRRL ! Number of liquid water var.
INTEGER,                  INTENT(IN)    :: KRRI ! Number of ice water var.
!
REAL,                     INTENT(IN)    :: PDRYMASST   ! Mass of dry air and of
REAL,                     INTENT(IN)    :: PREFMASS    ! the ref. atmosphere
                                          !  contained in the simulation domain
REAL,                     INTENT(IN)    :: PMASS_O_PHI0 !    Mass / Phi0 
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT        ! Temperature and water
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT         !  variables at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ      ! dry Density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHETAV     ! virtual potential temp.
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF    ! dry Density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHVREF     ! Virtual Temperature
                                                  ! of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVREF ! vapor mixing ratio 
                                       ! for the reference state 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF! Exner function of the
                                                  ! reference state
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PPHIT  ! Perturbation of
               ! either the Exner function Pi or Pi * Cpd * THvref
!
!
!*       0.2   Declarations of local variables :
!
INTEGER             :: IKU       ! Upper dimension in z direction
INTEGER             :: IIB       ! indice I Beginning in x direction
INTEGER             :: IJB       ! indice J Beginning in y direction
INTEGER             :: IKB       ! indice K Beginning in z direction
INTEGER             :: IIE       ! indice I End       in x direction 
INTEGER             :: IJE       ! indice J End       in y direction 
INTEGER             :: IKE       ! indice K End       in z direction 
INTEGER             :: JI        ! Loop index in x direction
INTEGER             :: JJ        ! Loop index in y direction      
INTEGER             :: JK        ! Loop index in z direction       
REAL     ::  ZP00_O_RD     ! = P00 /  Rd
REAL     ::  ZCVD_O_RD     ! = Cvd /  Rd
REAL     ::   ZRV_O_RD     ! = Rv  /  Rd
REAL     ::  ZCVD_O_RDCPD  ! = Cvd / (Rd * Cpd)
REAL     ::  ZMASS_O_PI    !    Mass / Pi0 
REAL     ::  ZMASSGUESS    ! guess of mass resulting of the pressure function
                                       ! provided by the pressure solveur, to an arbitary constant
REAL     ::  ZWATERMASST   ! Total mass of water Mw
!JUAN16
REAL, ALLOCATABLE, DIMENSION(:,:)     :: ZMASS_O_PI_2D,ZMASSGUESS_2D,ZWATERMASST_2D
!JUAN16
REAL     ::  ZPI0          ! constant to retrieve the absolute Exner pressure
INTEGER  ::  JWATER        ! loop index on the different types of water
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3))       &
         ::  ZRTOT, ZRHOREF, ZWORK
REAL     ::  ZPHI0
!
INTEGER  :: IINFO_ll
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES:
!              ----------------------------------------------
!
IKU = SIZE(PTHT,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
ALLOCATE(ZMASS_O_PI_2D(IIB:IIE,IJB:IJE))
ALLOCATE(ZMASSGUESS_2D(IIB:IIE,IJB:IJE))
ALLOCATE(ZWATERMASST_2D(IIB:IIE,IJB:IJE))
!
ZP00_O_RD = XP00 / XRD
ZCVD_O_RD = (XCPD - XRD) / XRD
!
!-------------------------------------------------------------------------------
!
!
!*       2.     COMPUTES THE ABSOLUTE EXNER FUNCTION (MAE+ DUR) 
!	        -----------------------------------------------
!
!       
!
IF ( CEQNSYS=='DUR' .OR. CEQNSYS=='MAE' ) THEN
!
  IF(KRR > 0) THEN
  !
  !   compute the mixing ratio of the total water (ZRTOT)
    ZRTOT(:,:,:) = PRT(:,:,:,1)
    DO JWATER = 2 , 1+KRRL+KRRI                
      ZRTOT(:,:,:) = ZRTOT(:,:,:) + PRT(:,:,:,JWATER)
    END DO
  ELSE
    ZRTOT(:,:,:) = 0.
  END IF
  !
  ZMASSGUESS_2D  = 0.  
  ZMASS_O_PI_2D  = 0.      
  ZWATERMASST_2D = 0.
!
  IF ( CEQNSYS == 'DUR' ) THEN
    ! compute the Jacobian in ZWORK
    IF ( SIZE(PRVREF,1) == 0 ) THEN
      ZWORK(:,:,:)=  PRHODJ * XTH00  / ( PRHODREF * PTHVREF )
    ELSE
      ZWORK(:,:,:)=PRHODJ * XTH00  &
           / ( PRHODREF * PTHVREF * (1. + PRVREF) )
    END IF
    !
    DO JK = IKB,IKE
      DO JJ = IJB,IJE
        DO JI = IIB,IIE
          ZMASSGUESS_2D(JI,JJ)  = ZMASSGUESS_2D(JI,JJ) +                          &
             (PEXNREF(JI,JJ,JK)+PPHIT(JI,JJ,JK))**ZCVD_O_RD   &
             * ZWORK(JI,JJ,JK) / PTHETAV(JI,JJ,JK)
          ZMASS_O_PI_2D(JI,JJ)  = ZMASS_O_PI_2D(JI,JJ) + ZWORK(JI,JJ,JK) / PTHETAV(JI,JJ,JK)
          ZWATERMASST_2D(JI,JJ) = ZWATERMASST_2D(JI,JJ) +       &
            ZRTOT(JI,JJ,JK) * ZWORK(JI,JJ,JK) * PRHODREF(JI,JJ,JK)
        END DO
      END DO
    END DO
!
  ELSE
    DO JK = IKB,IKE
      DO JJ = IJB,IJE
        DO JI = IIB,IIE
          ZMASSGUESS_2D(JI,JJ)  = ZMASSGUESS_2D(JI,JJ) +                               &
             (PEXNREF(JI,JJ,JK)+PPHIT(JI,JJ,JK))**ZCVD_O_RD        &
            * PRHODJ(JI,JJ,JK) / PRHODREF(JI,JJ,JK)                &
            / PTHETAV(JI,JJ,JK)
          ZMASS_O_PI_2D(JI,JJ)  = ZMASS_O_PI_2D(JI,JJ) +                               &
            PRHODJ(JI,JJ,JK) / PRHODREF(JI,JJ,JK) / PTHETAV(JI,JJ,JK)
          ZWATERMASST_2D(JI,JJ) = ZWATERMASST_2D(JI,JJ) + ZRTOT(JI,JJ,JK) * PRHODJ(JI,JJ,JK) 
        END DO
      END DO
    END DO
  END IF
!
  !
  ZMASSGUESS  = SUM_DD_R2_ll(ZMASSGUESS_2D)
  ZMASS_O_PI  = SUM_DD_R2_ll(ZMASS_O_PI_2D)
  ZWATERMASST = SUM_DD_R2_ll(ZWATERMASST_2D)
  !
  ZMASS_O_PI  = ZMASS_O_PI*ZP00_O_RD*ZCVD_O_RD
  ZPI0 = (PDRYMASST + ZWATERMASST - ZP00_O_RD*ZMASSGUESS ) / ZMASS_O_PI
  PPHIT(:,:,:) = PPHIT(:,:,:) + ZPI0
!
!
  !
  !          Second iteration
  !
  ZMASSGUESS_2D  = 0.
  IF ( CEQNSYS == 'DUR' ) THEN
    DO JK = IKB,IKE
      DO JJ = IJB,IJE
        DO JI = IIB,IIE
          ZMASSGUESS_2D(JI,JJ)  = ZMASSGUESS_2D(JI,JJ) +                               &
           (PEXNREF(JI,JJ,JK)+PPHIT(JI,JJ,JK))**ZCVD_O_RD          &
          * ZWORK(JI,JJ,JK) / PTHETAV(JI,JJ,JK)
        END DO
      END DO
    END DO
  ELSE
    DO JK = IKB,IKE
      DO JJ = IJB,IJE
        DO JI = IIB,IIE
          ZMASSGUESS_2D(JI,JJ)  = ZMASSGUESS_2D(JI,JJ) +                                &
            (PEXNREF(JI,JJ,JK)+PPHIT(JI,JJ,JK))**ZCVD_O_RD          &
           * PRHODJ(JI,JJ,JK) / PRHODREF(JI,JJ,JK) / PTHETAV(JI,JJ,JK)
        END DO
      END DO
    END DO
  END IF
!

  ZMASSGUESS  = SUM_DD_R2_ll(ZMASSGUESS_2D)
  !
  ZPI0 = (PDRYMASST + ZWATERMASST - ZP00_O_RD*ZMASSGUESS ) / ZMASS_O_PI
  PPHIT(:,:,:) = PPHIT(:,:,:) + ZPI0
!
!
ELSEIF( CEQNSYS == 'LHE' ) THEN
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTES THE ABSOLUTE PRESSURE FUNCTION (LHE) 
!	        ---------------------------------------------
!
  !               compute the reference moist density
  !
  ZCVD_O_RDCPD = ZCVD_O_RD / XCPD
  ZCVD_O_RD = (XCPD - XRD) / XRD
  !
  IF (LBOUSS) THEN
    ZRHOREF(:,:,:) = PRHODREF(:,:,:)
  ELSE
    ZRHOREF(:,:,:) = PEXNREF(:,:,:) ** ZCVD_O_RD    &
                  * XP00 / ( XRD * PTHVREF(:,:,:) )
  ENDIF        
  !
  !
  !               compute the virtual potential temperature 
  !
  !
  IF(KRR > 0) THEN
  !
  !   compute the mixing ratio of the total water (ZRRTOT)
    ZRV_O_RD = XRV / XRD
    ZRTOT(:,:,:) = PRT(:,:,:,1)
    DO JWATER = 2 , 1+KRRL+KRRI                
      ZRTOT(:,:,:) = ZRTOT(:,:,:) + PRT(:,:,:,JWATER)
    END DO
  !   compute the virtual potential temperature in ZWORK                 
    ZWORK(:,:,:) = PTHT(:,:,:) * (1. + PRT(:,:,:,1) * ZRV_O_RD)  &
                                / (1. + ZRTOT(:,:,:))
  ELSE
  !   compute the virtual potential temperature when water is absent
    ZWORK(:,:,:)  = PTHT(:,:,:)
    ZRTOT(:,:,:) = 0.
  END IF
  !
  !
  !               compute the absolute pressure function 
  !
  !
  !
  ZMASSGUESS_2D  = 0. 
  ZWATERMASST_2D = 0.
!
  DO JK = IKB,IKE
    DO JJ = IJB,IJE
      DO JI = IIB,IIE
        ZMASSGUESS_2D(JI,JJ)  = ZMASSGUESS_2D(JI,JJ) + ZRHOREF(JI,JJ,JK) /  PTHVREF(JI,JJ,JK) *   &
                     (  ZWORK(JI,JJ,JK)                                       &
                      - ZCVD_O_RDCPD * PPHIT(JI,JJ,JK) / PEXNREF(JI,JJ,JK)    &
                     ) * PRHODJ(JI,JJ,JK) /  PRHODREF(JI,JJ,JK)
        ZWATERMASST_2D(JI,JJ) = ZWATERMASST_2D(JI,JJ) + ZRTOT(JI,JJ,JK) * PRHODJ(JI,JJ,JK)
      END DO
    END DO
  END DO
  !
  ZMASSGUESS  = SUM_DD_R2_ll(ZMASSGUESS_2D)
  ZWATERMASST =  SUM_DD_R2_ll(ZWATERMASST_2D)
  !
  ZPHI0 = (PDRYMASST + ZWATERMASST - 2. * PREFMASS + ZMASSGUESS ) / PMASS_O_PHI0
  PPHIT(:,:,:) = PPHIT(:,:,:) + ZPHI0
  !
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE P_ABS
