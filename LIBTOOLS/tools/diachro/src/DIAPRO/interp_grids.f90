!     ######spl
      SUBROUTINE INTERP_GRIDS(K)
!     ##########################
!
!!****  *INTERP_GRIDS* - 
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!      
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist (former NCAR common)
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       25/10/99
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ALLOC_FORDIACHRO
USE MODD_RESOLVCAR
USE MODD_PVT
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_NMGRID

IMPLICIT NONE
!
!*       0.1   Declaration of arguments and results
!
INTEGER :: K
!
!*       0.2   Local Variables 
!
INTEGER :: IG, IP, IN, II, IJ, IK, IT
INTEGER :: IGRIDIA
REAL,DIMENSION(:,:,:,:),ALLOCATABLE :: ZTEM
!
!-------------------------------------------------------------------------------
!
!*      1.    Preliminary calculations
!             ------------------------
!
IF(LPRESYT)THEN
  IG=1
  IN=1
  IP=1
  II=SIZE(XPRES,1)
  IJ=SIZE(XPRES,2)
  IK=SIZE(XPRES,3)
  IT=1
  ALLOCATE(ZTEM(SIZE(XPRES,1),SIZE(XPRES,2),SIZE(XPRES,3),1))
  IGRIDIA=NMGRID
ELSE
  IG=NGRIDIA(NPROCDIA(NBPROCDIA(K),K))
  IN=NNDIA(1,K)
  IP=NPROCDIA(NBPROCDIA(K),K)
  II=SIZE(XVAR,1)
  IJ=SIZE(XVAR,2)
  IK=SIZE(XVAR,3)
  IT=SIZE(XVAR,4)
  ALLOCATE(ZTEM(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4)))
  IGRIDIA=NGRIDIAM
ENDIF
!
! Mars 20000 Cas ou pas de decalage horizontal mais vertical peut-etre
!
IF(LPRESYT)THEN
  ZTEM(:,:,:,1)=XPRES(:,:,:,NLOOPT,IN,IP)
ELSE
  ZTEM(:,:,:,:)=XVAR(:,:,:,:,IN,IP)
ENDIF
!
! Decalages horizontaux
!
!SELECT CASE(NGRIDIAM)
IF(II /=1 .AND. IJ /=1)THEN
SELECT CASE(IGRIDIA)
  CASE(1)
    SELECT CASE(IG)
      CASE(2,6)
        IF(LPRESYT)THEN
	ZTEM(1:II-1,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(1:II-1,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(II,:,:,:)=2.*ZTEM(II-1,:,:,:)-ZTEM(II-2,:,:,:)
      CASE(3,7)
        IF(LPRESYT)THEN
	ZTEM(:,1:IJ-1,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,1:IJ-1,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
      CASE(5)
        IF(LPRESYT)THEN
	ZTEM(1:II-1,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(1:II-1,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(II,:,:,:)=2.*ZTEM(II-1,:,:,:)-ZTEM(II-2,:,:,:)
	ZTEM(:,1:IJ-1,:,:)=.5*(ZTEM(:,1:IJ-1,:,:)+ZTEM(:,2:IJ,:,:))
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
    END SELECT
  CASE(2)
    SELECT CASE(IG)
      CASE(1,4)
        IF(LPRESYT)THEN
	ZTEM(2:II,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(2:II,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(1,:,:,:)=2.*ZTEM(2,:,:,:)-ZTEM(3,:,:,:)
      CASE(3,7)
        IF(LPRESYT)THEN
	ZTEM(:,1:IJ-1,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,1:IJ-1,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
	ZTEM(2:II,:,:,:)=.5*(ZTEM(1:II-1,:,:,:)+ZTEM(2:II,:,:,:))
        ZTEM(1,:,:,:)=2.*ZTEM(2,:,:,:)-ZTEM(3,:,:,:)
      CASE(5)
        IF(LPRESYT)THEN
	ZTEM(:,1:IJ-1,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,1:IJ-1,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
    END SELECT
  CASE(3)
    SELECT CASE(IG)
      CASE(1,4)
        IF(LPRESYT)THEN
	ZTEM(:,2:IJ,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,2:IJ,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,1,:,:)=2.*ZTEM(:,2,:,:)-ZTEM(:,3,:,:)
      CASE(2,6)
        IF(LPRESYT)THEN
	ZTEM(2:II,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(2:II,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(1,:,:,:)=2.*ZTEM(2,:,:,:)-ZTEM(3,:,:,:)
	ZTEM(:,2:IJ,:,:)=.5*(ZTEM(:,1:IJ-1,:,:)+ZTEM(:,2:IJ,:,:))
        ZTEM(:,1,:,:)=2.*ZTEM(:,2,:,:)-ZTEM(:,3,:,:)
      CASE(5)
        IF(LPRESYT)THEN
	ZTEM(1:II-1,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(1:II-1,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(II,:,:,:)=2.*ZTEM(II-1,:,:,:)-ZTEM(II-2,:,:,:)
    END SELECT
  CASE(4)
    SELECT CASE(IG)
      CASE(3,7)
        IF(LPRESYT)THEN
	ZTEM(:,1:IJ-1,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,1:IJ-1,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
      CASE(2,6)
        IF(LPRESYT)THEN
	ZTEM(1:II-1,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(1:II-1,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(II,:,:,:)=2.*ZTEM(II-1,:,:,:)-ZTEM(II-2,:,:,:)
      CASE(5)
        IF(LPRESYT)THEN
	ZTEM(1:II-1,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(1:II-1,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(II,:,:,:)=2.*ZTEM(II-1,:,:,:)-ZTEM(II-2,:,:,:)
	ZTEM(:,1:IJ-1,:,:)=.5*(ZTEM(:,1:IJ-1,:,:)+ZTEM(:,2:IJ,:,:))
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
    END SELECT
  CASE(5)
    SELECT CASE(IG)
      CASE(1,4)
        IF(LPRESYT)THEN
	ZTEM(2:II,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(2:II,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(1,:,:,:)=2.*ZTEM(2,:,:,:)-ZTEM(3,:,:,:)
	ZTEM(:,2:IJ,:,:)=.5*(ZTEM(:,1:IJ-1,:,:)+ZTEM(:,2:IJ,:,:))
        ZTEM(:,1,:,:)=2.*ZTEM(:,2,:,:)-ZTEM(:,3,:,:)
      CASE(2,6)
        IF(LPRESYT)THEN
	ZTEM(:,2:IJ,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,2:IJ,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,1,:,:)=2.*ZTEM(:,2,:,:)-ZTEM(:,3,:,:)
      CASE(3,7)
        IF(LPRESYT)THEN
	ZTEM(2:II,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(2:II,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(1,:,:,:)=2.*ZTEM(2,:,:,:)-ZTEM(3,:,:,:)
    END SELECT
  CASE(6)
    SELECT CASE(IG)
      CASE(1,4)
        IF(LPRESYT)THEN
	ZTEM(2:II,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(2:II,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(1,:,:,:)=2.*ZTEM(2,:,:,:)-ZTEM(3,:,:,:)
      CASE(3,7)
        IF(LPRESYT)THEN
	ZTEM(:,1:IJ-1,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,1:IJ-1,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
	ZTEM(2:II,:,:,:)=.5*(ZTEM(1:II-1,:,:,:)+ZTEM(2:II,:,:,:))
        ZTEM(1,:,:,:)=2.*ZTEM(2,:,:,:)-ZTEM(3,:,:,:)
      CASE(5)
        IF(LPRESYT)THEN
	ZTEM(:,1:IJ-1,:,1)=.5*(XPRES(:,1:IJ-1,:,NLOOPT,IN,IP)+XPRES(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,1:IJ-1,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,IJ,:,:)=2.*ZTEM(:,IJ-1,:,:)-ZTEM(:,IJ-2,:,:)
    END SELECT
  CASE(7)
    SELECT CASE(IG)
      CASE(1,4)
        IF(LPRESYT)THEN
	ZTEM(:,2:IJ,:,1)=.5*(XVAR(:,1:IJ-1,:,NLOOPT,IN,IP)+XVAR(:,2:IJ,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(:,2:IJ,:,:)=.5*(XVAR(:,1:IJ-1,:,:,IN,IP)+XVAR(:,2:IJ,:,:,IN,IP))
	ENDIF
        ZTEM(:,1,:,:)=2.*ZTEM(:,2,:,:)-ZTEM(:,3,:,:)
      CASE(2,6)
        IF(LPRESYT)THEN
	ZTEM(1:II-1,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(1:II-1,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(II,:,:,:)=2.*ZTEM(II-1,:,:,:)-ZTEM(II-2,:,:,:)
	ZTEM(:,2:IJ,:,:)=.5*(ZTEM(:,1:IJ-1,:,:)+ZTEM(:,2:IJ,:,:))
        ZTEM(:,1,:,:)=2.*ZTEM(:,2,:,:)-ZTEM(:,3,:,:)
      CASE(5)
        IF(LPRESYT)THEN
	ZTEM(1:II-1,:,:,1)=.5*(XPRES(1:II-1,:,:,NLOOPT,IN,IP)+XPRES(2:II,:,:,NLOOPT,IN,IP))
	ELSE
	ZTEM(1:II-1,:,:,:)=.5*(XVAR(1:II-1,:,:,:,IN,IP)+XVAR(2:II,:,:,:,IN,IP))
	ENDIF
        ZTEM(II,:,:,:)=2.*ZTEM(II-1,:,:,:)-ZTEM(II-2,:,:,:)
    END SELECT
END SELECT
ENDIF
!
! Decalages VERTICAUX
!
IF(IK /= 1)THEN
SELECT CASE(NGRIDIAM)
  CASE(1,2,3,5)
    SELECT CASE(IG)
      CASE(4,6,7)
	ZTEM(:,:,1:IK-1,:)=.5*(ZTEM(:,:,1:IK-1,:)+ZTEM(:,:,2:IK,:))
        ZTEM(:,:,IK,:)=2.*ZTEM(:,:,IK-1,:)-ZTEM(:,:,IK-2,:)
    END SELECT
  CASE(4,6,7)
    SELECT CASE(IG)
      CASE(1,2,3,5)
	ZTEM(:,:,2:IK,:)=.5*(ZTEM(:,:,1:IK-1,:)+ZTEM(:,:,2:IK,:))
        ZTEM(:,:,1,:)=2.*ZTEM(:,:,2,:)-ZTEM(:,:,3,:)
    END SELECT
END SELECT
ENDIF

IF(LPRESYT)THEN
  XPRES(:,:,:,NLOOPT,IN,IP)=ZTEM(:,:,:,1)
ELSE
  XVAR(:,:,:,:,IN,IP)=ZTEM(:,:,:,:)
ENDIF

DEALLOCATE(ZTEM)
!
!----------------------------------------------------------------------------
!
!*     3.    EXIT
!            ----
!
RETURN
END SUBROUTINE INTERP_GRIDS
