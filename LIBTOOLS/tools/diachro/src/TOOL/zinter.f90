!     ##################
      MODULE MODI_ZINTER
!     ##################
!
INTERFACE ZINTER
      SUBROUTINE ZINTER(PVMNH,PZGMNH,PVZL,PLZL,KKB,PUNDEF,KNIVMOD)
!
REAL,DIMENSION(:,:,:),INTENT(IN) :: PVMNH
REAL,DIMENSION(:,:,:),INTENT(IN) :: PZGMNH 
REAL,DIMENSION(:,:,:),INTENT(OUT):: PVZL 
REAL,DIMENSION(:),INTENT(IN)     :: PLZL
REAL,INTENT(IN)         :: PUNDEF
!
INTEGER,INTENT(IN)      :: KKB
INTEGER,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: KNIVMOD 
!
END SUBROUTINE ZINTER
!
      SUBROUTINE SINTER(PVMNH,PZGMNH,PVZL,PLZL,KKB,PUNDEF,KNIVMOD)
!
REAL,DIMENSION(:,:,:),INTENT(IN) :: PVMNH
REAL,DIMENSION(:,:,:),INTENT(IN) :: PZGMNH 
REAL,DIMENSION(:,:,:),INTENT(OUT):: PVZL 
REAL,DIMENSION(:,:,:),INTENT(IN) :: PLZL
REAL,INTENT(IN)         :: PUNDEF
!
INTEGER,INTENT(IN)      :: KKB
INTEGER,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: KNIVMOD 
!
END SUBROUTINE SINTER
!
END INTERFACE ZINTER
END MODULE MODI_ZINTER
!     ##################
      MODULE MODI_SINTER
!     ##################
!
INTERFACE SINTER
      SUBROUTINE SINTER(PVMNH,PZGMNH,PVZL,PLZL,KKB,PUNDEF,KNIVMOD)
!
REAL,DIMENSION(:,:,:),INTENT(IN) :: PVMNH
REAL,DIMENSION(:,:,:),INTENT(IN) :: PZGMNH 
REAL,DIMENSION(:,:,:),INTENT(OUT):: PVZL 
REAL,DIMENSION(:,:,:),INTENT(IN) :: PLZL
REAL,INTENT(IN)         :: PUNDEF
!
INTEGER,INTENT(IN)      :: KKB
INTEGER,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: KNIVMOD 
!
END SUBROUTINE SINTER
END INTERFACE SINTER
END MODULE MODI_SINTER
!
!------------------------------------------------------------------------------
!
!     ####################################################
      SUBROUTINE SINTER(PVMNH,PZGMNH,PVZL,PLZL,KKB,PUNDEF,KNIVMOD)
!     ####################################################
!
!
!!****  *ZINTER * - routine to linearly interpolate
!!
!!     PURPOSE
!!     -------
!    This routine interpolates an input field on Gal-Chen grid, linearly in 
!    another Z-grid (regular or not).
!
!!**   METHOD
!!     ------
!!
!!
!!     EXTERNAL
!!     --------
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!      None
!!
!!     REFERENCE
!!     ---------
!!      Research manual 2 ECMWF forecast model, 1988, Ref M1.6/3
!!      "adiabatic part", Appendix 6 postprocessing
!!      Section 3.  Vertical interpolation, p. A6.5-6
!!      Section 3.4 Extrapolation, pp. A6.6-7
!!
!!     AUTHOR
!!     ------
!!       P. Mascart     * LA *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original       22/04/96
!!       Modification   11/02/99 Chaboureau - some simplifications
!!-----------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*       0.1  Declaration of arguments 
!
INTEGER,INTENT(IN)           :: KKB  ! 1st level above ground    
REAL,DIMENSION(:,:,:),INTENT(IN) :: PVMNH
!!  PVMNH  = tableau du champ donne au points masse Meso-NH
REAL,DIMENSION(:,:,:),INTENT(IN) :: PZGMNH 
!!  PZGMNH = altitude geopotentiel au point masse Meso-NH
REAL,DIMENSION(:,:,:),INTENT(IN) :: PLZL ! list of the new vertical levels
REAL,DIMENSION(:,:,:),INTENT(OUT):: PVZL ! interpolated output field
REAL,INTENT(IN)         :: PUNDEF  ! undefined value
INTEGER,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: KNIVMOD 
!!                                   first model level above PLZL(:,:,1)
!
!*       0.2  Declaration of local variables
!
INTEGER   :: ILT,ILN  ! number of points in the 1st and 2nd dimensions
INTEGER   :: IKU      ! number of input vertical levels
INTEGER   :: INP      ! number of new vertical levels (1: base ; INP: top)

REAL      :: ZSLOPE
INTEGER   :: JI,JJ,JKZL,JK
INTEGER   :: IKD
!
!------------------------------------------------------------------------------
!
!*       1.   INITIALIZATION
!             --------------
!
ILT=SIZE(PVMNH,1)
ILN=SIZE(PVMNH,2)
IKU=SIZE(PVMNH,3)
INP=SIZE(PVZL,3)
PVZL=PUNDEF
IF (PRESENT (KNIVMOD)) KNIVMOD=KKB
!
print*,'in SINTER ',ILT,ILN,IKU,INP
!------------------------------------------------------------------------------
!
!*       2.   INTERPOLATION
!             -------------
!
OX: DO  JI =1,ILT
  OY:  DO   JJ =1,ILN
    PLEV:  DO   JKZL=1,INP
      !
      !   i) Zones flagging
      !
      IKD=0
      IF(PLZL(JI,JJ,JKZL).GE.PZGMNH(JI,JJ,IKU))       IKD=10*IKU
      DO  JK  =IKU-1,KKB,-1
         IF((PZGMNH(JI,JJ,JK+1).GT.PLZL(JI,JJ,JKZL)).AND.   &
           (PLZL(JI,JJ,JKZL).GE.PZGMNH(JI,JJ,JK)))    IKD=JK
      END DO
      IF(PLZL(JI,JJ,JKZL).LT.PZGMNH(JI,JJ,KKB))       IKD=-10*IKU
      IF(IKD==0) IKD=10*IKU  !! pas propre...
      !
      !   ii) Regular points interpolation
      !
      IF(ABS(IKD).NE.(10*IKU)) THEN
        IF ( PVMNH(JI,JJ,IKD) /= PUNDEF .AND. PVMNH(JI,JJ,IKD+1)/= PUNDEF) THEN
          ZSLOPE=(PLZL(JI,JJ,JKZL)-PZGMNH(JI,JJ,IKD))      &
                 /(PZGMNH(JI,JJ,IKD+1)-PZGMNH(JI,JJ,IKD))
          PVZL(JI,JJ,JKZL)=PVMNH(JI,JJ,IKD)                &
                           +ZSLOPE*(PVMNH(JI,JJ,IKD+1)-PVMNH(JI,JJ,IKD))
          IF (PRESENT (KNIVMOD)) THEN
            KNIVMOD(JI,JJ)=IKD+1
          ENDIF
        ELSE
          PVZL(JI,JJ,JKZL)=PUNDEF
        ENDIF
      ELSE
      !
      !   iii) No extrapolation below the ground and above the top
      !
        PVZL(JI,JJ,JKZL)=PUNDEF
      ENDIF
    END DO PLEV
  END DO OY 
END DO OX
!
END SUBROUTINE SINTER
!
!     ####################################################
      SUBROUTINE ZINTER(PVMNH,PZGMNH,PVZL,PLZL,KKB,PUNDEF,KNIVMOD)
!     ####################################################
!
!
!!****  *ZINTER * - routine to linearly interpolate
!!
!!     PURPOSE
!!     -------
!    This routine interpolates an input field on Gal-Chen grid, linearly in 
!    another Z-grid (regular or not).
!
!!**   METHOD
!!     ------
!!
!!
!!     EXTERNAL
!!     --------
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!      None
!!
!!     REFERENCE
!!     ---------
!!      Research manual 2 ECMWF forecast model, 1988, Ref M1.6/3
!!      "adiabatic part", Appendix 6 postprocessing
!!      Section 3.  Vertical interpolation, p. A6.5-6
!!      Section 3.4 Extrapolation, pp. A6.6-7
!!
!!     AUTHOR
!!     ------
!!       P. Mascart     * LA *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original       22/04/96
!!       Modification   11/02/99 Chaboureau - some simplifications
!!-----------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODI_SINTER
IMPLICIT NONE
!
!*       0.1  Declaration of arguments 
!
REAL,DIMENSION(:,:,:),INTENT(IN) :: PVMNH
REAL,DIMENSION(:,:,:),INTENT(IN) :: PZGMNH 
REAL,DIMENSION(:,:,:),INTENT(OUT):: PVZL 
REAL,DIMENSION(:),INTENT(IN)     :: PLZL
REAL,INTENT(IN)         :: PUNDEF
!
INTEGER,INTENT(IN)      :: KKB
INTEGER,DIMENSION(:,:),INTENT(OUT),OPTIONAL:: KNIVMOD 
!
!*       0.2  Declaration of local variables
!
INTEGER   :: ILT,ILN  ! number of points in the 1st and 2nd dimensions
INTEGER   :: IKU      ! number of input vertical levels
INTEGER   :: INP      ! number of new vertical levels (1: base ; INP: top)
REAL,DIMENSION(:,:,:),ALLOCATABLE :: ZLZL 
!
!------------------------------------------------------------------------------
!
!*       1.   INITIALIZATION
!             --------------
!
ILT=SIZE(PVMNH,1)
ILN=SIZE(PVMNH,2)
INP=SIZE(PVZL,3)
!
ALLOCATE(ZLZL(ILT,ILN,INP))
ZLZL(:,:,:) = SPREAD( SPREAD( PLZL(1:INP),1,ILT ) ,2,ILN )
!
!------------------------------------------------------------------------------
!
!*       2.   INTERPOLATION
!             -------------
!
CALL SINTER(PVMNH,PZGMNH,PVZL,ZLZL,KKB,PUNDEF,KNIVMOD)
!
END SUBROUTINE ZINTER
