!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################################
      MODULE MODI_COMPUTE_EXNER_FROM_GROUND
!     #####################################
INTERFACE COMPUTE_EXNER_FROM_GROUND
            SUBROUTINE COMPUTE_EXNER_FROM_GROUND3D(PTHV,PZFLUX,PEXNSURF2D, &
                                                   PEXNFLUX,PEXNMASS)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHV      ! virtual potential temperature
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX    ! altitude of flux points
REAL, DIMENSION(:,:),   INTENT(IN)  :: PEXNSURF2D! ground Exner function
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNFLUX  ! Exner function at flux points
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNMASS  ! Exner function at mass points
!
END SUBROUTINE COMPUTE_EXNER_FROM_GROUND3D
!
            SUBROUTINE COMPUTE_EXNER_FROM_GROUND1D(PTHV,PZFLUX,PEXNSURF, &
                                                   PEXNFLUX,PEXNMASS)
!
REAL, DIMENSION(:), INTENT(IN)  :: PTHV      ! virtual potential temperature
REAL, DIMENSION(:), INTENT(IN)  :: PZFLUX    ! altitude of flux points
REAL,               INTENT(IN)  :: PEXNSURF  ! ground Exner function
REAL, DIMENSION(:), INTENT(OUT) :: PEXNFLUX  ! Exner function at flux points
REAL, DIMENSION(:), INTENT(OUT) :: PEXNMASS  ! Exner function at mass points
!
END SUBROUTINE COMPUTE_EXNER_FROM_GROUND1D
!
END INTERFACE 
END MODULE MODI_COMPUTE_EXNER_FROM_GROUND
!     ######################################
      MODULE MODI_COMPUTE_EXNER_FROM_GROUND3
!     ######################################
INTERFACE
            SUBROUTINE COMPUTE_EXNER_FROM_GROUND3D(PTHV,PZFLUX,PEXNSURF2D, &
                                                   PEXNFLUX,PEXNMASS)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHV      ! virtual potential temperature
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX    ! altitude of flux points
REAL, DIMENSION(:,:),   INTENT(IN) :: PEXNSURF2D! ground Exner function
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNFLUX  ! Exner function at flux points
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNMASS  ! Exner function at mass points
!
END SUBROUTINE COMPUTE_EXNER_FROM_GROUND3D
!
END INTERFACE
END MODULE MODI_COMPUTE_EXNER_FROM_GROUND3
!     ########################################################################
      SUBROUTINE COMPUTE_EXNER_FROM_GROUND3D(PTHV,PZFLUX,PEXNSURF2D,PEXNFLUX,PEXNMASS)
!     ########################################################################
!
!!****  *COMPUTE_EXNER_FROM_GROUND3D* - computation of hydrostatic 
!!                                      Exner function from ground
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!  1 The local Exner function in computed by integration of the hydrostatic
!!    relation from ground (PEXNSURF2D) to top.
!!
!!    dPI= -g/(Cpd thetav) dz
!!
!!  2 The Exner function at mass level is computed as follows and linearly 
!!    extrapolated for the uppest non-physical level:
!!
!!      ~           PI(k+1)-PI(k)
!!     PI(k) = -----------------------
!!                lnPI(k+1)-lnPI(k)
!!         
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!         LTHINSHELL : logical for thinshell approximation
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_CST       : contains physical constants
!!         XG  : gravity constant
!!         XCPD: specific heat for dry air at constant pressure
!!         XRD : gas constant for dry air
!!      Module MODD_PARAMETERS
!!         JPVEXT,JPHEXT
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
!!      Original    06/03/96
!!                  26/08/96 (V. Masson) thinshell approximation only available
!!                  03/12/02 (P. Jabouille)  add no thinshell condition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!* 0.    DECLARATIONS
!        ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_LUNIT
USE MODD_PARAMETERS
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTHV      ! virtual potential temperature
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PZFLUX    ! altitude of flux points
REAL, DIMENSION(:,:)  , INTENT(IN)  :: PEXNSURF2D! ground Exner function
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNFLUX  ! Exner function at flux points
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNMASS  ! Exner function at mass points
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: IKB,IKU,JK
REAL    :: ZD1      ! switch for thinshell approximation
REAL    :: ZGSCPD   ! = g/Cpd
REAL, DIMENSION(SIZE(PZFLUX,1),SIZE(PZFLUX,2),SIZE(PZFLUX,3)) :: ZZM
                    !altitude of mass points
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATIONS
!              ---------------
!
IKB=JPVEXT+1
IKU=SIZE(PZFLUX,3)
ZGSCPD = XG/XCPD
!
!-------------------------------------------------------------------------------
!
!*       2.   COMPUTATION OF THE EXNER FUNCTION AT FLUX POINTS
!             ------------------------------------------------
!
ZZM=MZF(1,IKU,1,PZFLUX)
PEXNFLUX(:,:,IKB)=PEXNSURF2D(:,:)
IF (LCARTESIAN .OR. LTHINSHELL) THEN
  ZD1=0.
ELSE
  ZD1=1.
ENDIF
DO JK=IKB+1,IKU
  PEXNFLUX(:,:,JK)=(PEXNFLUX(:,:,JK-1)*(1.+ZD1*2./7.*(PZFLUX(:,:,JK-1)-PZFLUX(:,:,JK))/ &
                                       (XRADIUS+ZZM(:,:,JK-1)))+ &
        ZGSCPD/PTHV(:,:,JK-1)*(PZFLUX(:,:,JK-1)-PZFLUX(:,:,JK)))/  &
 (1.-ZD1*2./7.*(PZFLUX(:,:,JK-1)-PZFLUX(:,:,JK))/(XRADIUS+ZZM(:,:,JK-1)))
END DO
DO JK=IKB-1,1,-1
  PEXNFLUX(:,:,JK)=(PEXNFLUX(:,:,JK+1)*(1.+ZD1*2./7.*(PZFLUX(:,:,JK+1)-PZFLUX(:,:,JK))/ &
                                       (XRADIUS+ZZM(:,:,JK)))+ &
        ZGSCPD/PTHV(:,:,JK)*(PZFLUX(:,:,JK+1)-PZFLUX(:,:,JK)))/  &
  (1.-ZD1*2./7.*(PZFLUX(:,:,JK+1)-PZFLUX(:,:,JK))/(XRADIUS+ZZM(:,:,JK)))
END DO
!
!-------------------------------------------------------------------------------
!                     
!*       3.   COMPUTATION OF EXNER FUNCTION AT MASS POINTS
!             --------------------------------------------
!
PEXNMASS(:,:,1:IKU-1)=(PEXNFLUX(:,:,1:IKU-1)-PEXNFLUX(:,:,2:IKU))            &
                     /(LOG(PEXNFLUX(:,:,1:IKU-1))-LOG(PEXNFLUX(:,:,2:IKU))) 
PEXNMASS(:,:,IKU)=  PEXNMASS(:,:,IKU-1)                                      &
                 +( PEXNMASS(:,:,IKU-1)-PEXNMASS(:,:,IKU-2) )                &
                 /( PZFLUX(:,:,IKU-1)-PZFLUX(:,:,IKU-2) )                    &
                 *( PZFLUX(:,:,IKU)-PZFLUX(:,:,IKU-1) )
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COMPUTE_EXNER_FROM_GROUND3D
!     ######################################################################
      SUBROUTINE COMPUTE_EXNER_FROM_GROUND1D(PTHV,PZFLUX,PEXNSURF,PEXNFLUX,PEXNMASS)
!     ######################################################################
!
!!****  *COMPUTE_EXNER_FROM_GROUND1D* - computation of hydrostatic 
!!                                      Exner function from model top
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
!!    IMPLICIT ARGUMENTS
!!    ------------------
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
!!      Original    06/03/96
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!* 0.    DECLARATIONS
!        ------------
!
USE MODI_COMPUTE_EXNER_FROM_GROUND3
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PTHV      ! virtual potential temperature
REAL, DIMENSION(:), INTENT(IN)  :: PZFLUX    ! altitude of flux points
REAL,               INTENT(IN)  :: PEXNSURF  ! ground Exner function
REAL, DIMENSION(:), INTENT(OUT) :: PEXNFLUX  ! Exner function at flux points
REAL, DIMENSION(:), INTENT(OUT) :: PEXNMASS  ! Exner function at mass points
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
REAL, DIMENSION(1,1,SIZE(PZFLUX))  :: ZTHV      ! virtual potential temperature
REAL, DIMENSION(1,1,SIZE(PZFLUX))  :: ZZFLUX    ! altitude of flux points
REAL, DIMENSION(1,1)               :: ZEXNSURF  ! ground Exner function
REAL, DIMENSION(1,1,SIZE(PZFLUX))  :: ZEXNFLUX  ! Exner function at flux points
REAL, DIMENSION(1,1,SIZE(PZFLUX))  :: ZEXNMASS  ! Exner function at mass points                                    
!
!-------------------------------------------------------------------------------
!
ZTHV(1,1,:)=PTHV(:)
ZZFLUX(1,1,:)=PZFLUX(:)
ZEXNSURF(1,1)=PEXNSURF
!
CALL COMPUTE_EXNER_FROM_GROUND3D(ZTHV,ZZFLUX,ZEXNSURF,ZEXNFLUX,ZEXNMASS)
!
PEXNFLUX(:)=ZEXNFLUX(1,1,:)
PEXNMASS(:)=ZEXNMASS(1,1,:)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COMPUTE_EXNER_FROM_GROUND1D
