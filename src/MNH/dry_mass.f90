!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_DRY_MASS
!     ####################
INTERFACE
      SUBROUTINE DRY_MASS(PTHV,PR,PJ,PPABS,PDRYMASS)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHV      ! virtual potential temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PR        ! water mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PJ        ! jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PPABS     ! absolute pressure
!
REAL,                   INTENT(OUT) :: PDRYMASS  ! dry mass
!
END SUBROUTINE DRY_MASS
END INTERFACE
END MODULE MODI_DRY_MASS
!     ######spl
      SUBROUTINE DRY_MASS(PTHV,PR,PJ,PPABS,PDRYMASS)
!     ##########################################
!
!!****  *DRY_MASS* - computation of the total dry air mass 
!!
!!    PURPOSE
!!    -------
!!    This routine computes the total dry mass in the whole domain from
!!    the virtual potential temperature, the mixing ratio, the local Exner
!!    function at the top of the model and the Jacobian. 
!!
!!**  METHOD
!!    ------
!!    
!!  1 The local Exner function in computed by integration of the hydrostatic
!!    relation from top (PEXNTOP2D) to bottom (routine COMPUTE_EXNER_FROM_TOP).
!!
!!
!!  2 The Exner function at mass level is computed as follows and linearly 
!!    extrapolated for the uppest non-physical level
!!    (routine COMPUTE_EXNER_FROM_TOP).
!!
!!  3 rhod is deduced by the relation:
!!
!!
!!                P / (PI)
!!    rhod=   ----------------
!!            Rd thetav (1+rw)
!!
!!  4 The total dry mass is deduced from rhod and the Jacobian (the integration
!!    is performed on the inner points):
!!               
!!    Md= SUM   rhod J*
!!       i,j,k
!!         
!!
!!    EXTERNAL
!!    --------
!!
!!    routine COMPUTE_EXNER_FROM_TOP : to compute the hydrostatic Exner function
!!    
!!    module MODI_COMPUTE_EXNER_FROM_TOP
!!    SUM3D_ll : distributed function equivalent to SUM
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_GRID1     : contains grid variables for model1
!!         XZZ   : altitude of the w points 
!!      Module MODD_CST       : contains physical constants
!!         XG  : gravity constant
!!         XCPD: specific heat for dry air at constant pressure
!!         XP00: reference pressure
!!         XRD : gas constant for dry air
!!      Module MODD_PARAMETERS
!!         JPVEXT;JPHEXT
!!      Module MODD_FIELD1    : contains the prognostic fields of model1
!!         XDRYMASST : total dry mass at t
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/12/94
!!                  Sept. 21, 1995  (J.Stein and V.Masson) surface pressure
!!                  Jan.  09, 1996  (V. Masson) hydrostatic pressure at mass
!!                                  point
!!                  March 06, 1996  (V. Masson) call to COMPUTE_EXNER_FROM_TOP
!!                  Jan   15, 1997 (Stein,Lafore) Durran anelastic equation 
!!                  Jun   10, 1997 (V. Masson) use absolute pressure
!!                  Jul   10, 1997 (V. Masson) removes use to modules of model 1
!!                  Nov   21, 1997 (V. Masson) use all water species
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!* 0.    DECLARATIONS
!        ------------
!
USE MODD_CONF                   ! declaration modules
USE MODD_LUNIT_n, ONLY : TLUOUT
USE MODD_CST
USE MODD_PARAMETERS
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHV      ! virtual potential temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PR        ! water mixing ratio
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PJ        ! jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PPABS     ! absolute pressure
!
REAL,                   INTENT(OUT) :: PDRYMASS  ! dry mass
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER :: JRR
REAL,   DIMENSION(SIZE(PJ,1),SIZE(PJ,2),SIZE(PJ,3)) :: ZRHOD, ZSUMR
INTEGER :: IINFO_ll       ! return code of parallel routine
!-------------------------------------------------------------------------------
!
!                     
!*       1.   COMPUTATION OF RHOD
!             -------------------
!
!
ZSUMR(:,:,:) = 0.
DO JRR=1,SIZE(PR,4)
  ZSUMR(:,:,:) = ZSUMR(:,:,:) + PR(:,:,:,JRR)
END DO
!
ZRHOD(:,:,:)=PPABS(:,:,:)/(PPABS(:,:,:)/XP00)**(XRD/XCPD) &
            /(XRD*PTHV(:,:,:)*(1.+ZSUMR(:,:,:)))
!
!-------------------------------------------------------------------------------
!
!*       2.   COMPUTATION OF THE TOTAL DRY MASS
!             ---------------------------------
!
PDRYMASS=SUM3D_ll(PJ(:,:,:)*ZRHOD(:,:,:),IINFO_ll)
!
!-------------------------------------------------------------------------------
!
WRITE(TLUOUT%NLU,*) 'Routine DRYMASS completed'
!
END SUBROUTINE DRY_MASS
