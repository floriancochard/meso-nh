! $Source$
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_VERT_COORD
!     ######################
!
INTERFACE 
!
      SUBROUTINE VERT_COORD(OSLEVE,PZS,PZSMT,PLEN1,PLEN2,PZHAT,PZZ)
!
LOGICAL,                INTENT(IN) :: OSLEVE! flag for Sleve coordinate
REAL, DIMENSION(:,:),   INTENT(IN) :: PZS   ! fine orography
REAL, DIMENSION(:,:),   INTENT(IN) :: PZSMT ! smooth orography
REAL,                   INTENT(IN) :: PLEN1 ! Decay scale for smooth topography
REAL,                   INTENT(IN) :: PLEN2 ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! Positions z in the cartesian plane
REAL, DIMENSION(:,:,:), INTENT(OUT):: PZZ   ! True altitude of the w grid-point
!
END SUBROUTINE VERT_COORD
!
END INTERFACE
!
END MODULE MODI_VERT_COORD
!
!
!
!     #############################
      SUBROUTINE VERT_COORD(OSLEVE,PZS,PZSMT,PLEN1,PLEN2,PZHAT,PZZ)
!     #############################
!
!!****  *VERT_COORD* computes smoothed orography for SLEVE coordinate
!!
!!    PURPOSE
!!    -------
!
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
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        nov 2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
LOGICAL,                INTENT(IN) :: OSLEVE! flag for Sleve coordinate
REAL, DIMENSION(:,:),   INTENT(IN) :: PZS   ! fine orography
REAL, DIMENSION(:,:),   INTENT(IN) :: PZSMT ! smooth orography
REAL,                   INTENT(IN) :: PLEN1 ! Decay scale for smooth topography
REAL,                   INTENT(IN) :: PLEN2 ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! Positions z in the cartesian plane
REAL, DIMENSION(:,:,:), INTENT(OUT):: PZZ   ! True altitude of the w grid-point
!
!
!*       0.2   declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
IF (OSLEVE) THEN
! Sleve coordinate
  CALL SLEVE_COORD(PZS,PZSMT,PLEN1,PLEN2,PZHAT,PZZ)
ELSE
! Gal Chen coordinate
  CALL GALCHEN_COORD(PZS,PZHAT,PZZ)
END IF
!
!-------------------------------------------------------------------------------
CONTAINS
!
!     #############################
      SUBROUTINE SLEVE_COORD(PZS,PZSMT,PLEN1,PLEN2,PZHAT,PZZ)
!     #############################
!
!!****  *SLEVE_COORD* computes smoothed orography for SLEVE coordinate
!!
!!    PURPOSE
!!    -------
!
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
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	G. Zangler      * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        nov 2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_PARAMETERS, ONLY : JPVEXT
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(IN) :: PZS   ! fine orography
REAL, DIMENSION(:,:),   INTENT(IN) :: PZSMT ! smooth orography
REAL,                   INTENT(IN) :: PLEN1 ! Decay scale for smooth topography
REAL,                   INTENT(IN) :: PLEN2 ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! Positions z in the cartesian plane
REAL, DIMENSION(:,:,:), INTENT(OUT):: PZZ   ! True altitude of the w grid-point
!
!*       0.2   declarations of local variables
!
INTEGER :: IIU        ! number of points in X direction
INTEGER :: IJU        ! number of points in Y direction
INTEGER :: IKU        ! number of points in Z direction
INTEGER :: IKE        ! upper physical point
!
REAL                                     :: ZH      ! model top
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZZSMALL ! small-scale topography deviation (PZS-PZSMT)
!
!-------------------------------------------------------------------------------
!
IIU = SIZE(PZZ,1)
IJU = SIZE(PZZ,2)
IKU = SIZE(PZZ,3)
IKE = IKU - JPVEXT
!
ZH = PZHAT(IKE+1)
!
ZZSMALL(:,:) = PZS(:,:) - PZSMT(:,:)   ! Small-scale topography deviation
!
! Sleve coordinate
PZZ(:,:,:) = SPREAD(SPREAD(PZHAT(:),1,IIU),2,IJU) + &
             SPREAD(PZSMT(:,:),3,IKU) * SINH( (ZH - & 
             SPREAD(SPREAD(PZHAT(:),1,IIU),2,IJU) ) /PLEN1 ) / &
             SINH( ZH /PLEN1 ) + &
             SPREAD(ZZSMALL(:,:),3,IKU) * SINH( (ZH - & 
             SPREAD(SPREAD(PZHAT(:),1,IIU),2,IJU) ) /PLEN2 ) / &
             SINH( ZH /PLEN2 )

! Ensure symmetry of layer depths below/above the true surface level
! This is essential (!) for a correct surface pressure gradient computation over sloping topography
!
PZZ(:,:,1) = 2.*PZZ(:,:,2)-PZZ(:,:,3)
!
!-------------------------------------------------------------------------------
END SUBROUTINE SLEVE_COORD
!
!-------------------------------------------------------------------------------
!
!     #############################
      SUBROUTINE GALCHEN_COORD(PZS,PZHAT,PZZ)
!     #############################
!
!!****  *GALCHEN_COORD* computes smoothed orography for Gal-Chen coordinate
!!
!!    PURPOSE
!!    -------
!
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
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	G. Zangler      * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        nov 2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_PARAMETERS, ONLY : JPVEXT
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(IN) :: PZS   ! fine orography
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! Positions z in the cartesian plane
REAL, DIMENSION(:,:,:), INTENT(OUT):: PZZ   ! True altitude of the w grid-point
!
!*       0.2   declarations of local variables
!
!
INTEGER :: IIU        ! number of points in X direction
INTEGER :: IJU        ! number of points in Y direction
INTEGER :: IKU        ! number of points in Z direction
INTEGER :: IKE        ! upper physical point
!
REAL                                     :: ZH       ! model top
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZCOEF    ! 1-zs/H 
!
!-------------------------------------------------------------------------------
!
IIU = SIZE(PZZ,1)
IJU = SIZE(PZZ,2)
IKU = SIZE(PZZ,3)
IKE = IKU - JPVEXT
!
ZH = PZHAT(IKE+1)
!
ZCOEF(:,:) = 1.-PZS(:,:)/ZH
PZZ(:,:,:) = SPREAD(SPREAD(PZHAT(1:IKU),1,IIU),2,IJU)  &
           * SPREAD(ZCOEF(1:IIU,1:IJU),3,IKU)          &
           + SPREAD(PZS(1:IIU,1:IJU),3,IKU)
!
! This is essential (!) for a correct surface pressure gradient computation over sloping topography
PZZ(:,:,1) = 2.*PZZ(:,:,2)-PZZ(:,:,3)
!
!-------------------------------------------------------------------------------
END SUBROUTINE GALCHEN_COORD
!
!-------------------------------------------------------------------------------
END SUBROUTINE VERT_COORD
