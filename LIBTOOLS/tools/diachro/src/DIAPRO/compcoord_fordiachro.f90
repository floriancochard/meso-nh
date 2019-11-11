!     ######spl
      SUBROUTINE COMPCOORD_FORDIACHRO(KGRID)
!     ######################################
!
!!****  *COMPCOORD_FORDIACHRO* - Computes gridpoint locations,
!!                     meshsizes and topography 
!!                    for all the possible grids, and true altitude where 
!!                    required.
!!
!!    PURPOSE
!!    -------
!       When called for the first time (KGRID=0), COMPCOORD_FORDIACHRO returns for 
!     the 7 possible grid locations:
!        - XHAT, YHAT, ZHAT values (meters) stored in:
!              XXX(:,1:7), XXY(:,1:7), XXZ(:,1:7)
!        - meshsizes values (meters):
!              XXDXHAT(:,1:7), XXDYHAT(:,1:7)
!        - topography altitudes values (meters):
!              XXZS(:,:,1:7)
!
!       When called subsequently (0<KGRID<8), COMPCOORD_FORDIACHRO returns the true
!     gridpoint altitude (meters) corresponding to the requested KGRID value 
!     in the XZZ(:,:,:) array.
!
!!**  METHOD
!!    ------
!!      Temporary arrays are allocated to store the grid point characteristics
!!    and de-allocated on exit. The 3D gridpoints locations are linearly
!!    interpolated to the expected grid location from their respective
!!    nominal locations. Altitudes are interpolated for the w-grid values,
!!    which are obtained directly from the Gal-Chen Sommerville formula.
!!      For XXX, XXY, XXZ, XXDXHAT, XXDYHAT, XXZS the last index is the grid
!!    selector KGRID ranging from 1 to 7 as follows:
!!    1 -> Mass grid, 
!!    2 -> U grid, 
!!    3 -> V grid, 
!!    4 -> W grid, 
!!    5 -> Vertical vorticity grid, 
!!    6 -> y-component vorticity grid,
!!    7 -> x-component vorticity grid.
!!    all the 7 values are prepared one for all in this subroutine and passed
!!    to the general TRACE environment to be used in the display process.
!!
!!      For the XZZ array the last index is the z direction one, not the grid 
!!    selector one.
!! 
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_COORD      : declares gridpoint coordinates (TRACE use)
!!
!!       XXX,XXY,XXZ   : coordinate values for all the MESO-NH grids
!!       XXDXHAT,XXDYHAT,XXDZHAT:   meshsize values for all the MESO-NH grids
!!       XXZS          : topography values for all the MESO_NH grids
!!
!!      Module MODD_GRID1      : declares grid variables (Model module)
!!
!!       XXHAT, XYHAT : x, y in the conformal or cartesian plane
!!       XZHAT        : Gal-Chen z level
!!       XZS          : topography zs
!!       XZZ          : true gridpoint z altitude
!! 
!!      Module MODD_DIM1       : Contains dimensions
!!
!!         NIMAX,NJMAX,NKMAX :  x, y, and z array dimensions
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!
!!         JPHEXT   : Horizontal external points number
!!         JPVEXT   : Vertical external points number
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
!!     The 7 MESO-NH grid types are defined in:
!!      - Asencio N. et al., 1994, "Le projet de modele non-hydrostatique
!!        commun CNRM-LA, specifications techniques",
!!        Note CNRM/GMME, 26, 139p, (pages 39 to 43).
!!
!!      - Fischer C., 1994, "File structure and content in the Meso-NH
!!        model", Meso-nh internal note, CNRM/GMME,  July 5.
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_COORD
USE MODD_DIM1
USE MODD_CONF
USE MODD_GRID1
USE MODD_PARAMETERS
USE MODD_MEMCV
USE MODD_RESOLVCAR
!
USE MODI_VERT_COORD
!
IMPLICIT NONE
!
!*       0.1   Local variables declarations
!
INTEGER           :: IIU, IJU, IKU

INTEGER           :: IIB, IJB, IKB

INTEGER           :: IIE, IJE, IKE
!
! Calcul des X, Y, Z aux points de masse
REAL,DIMENSION(:),ALLOCATABLE   :: ZXMASS, ZYMASS, ZZMASS 

REAL,DIMENSION(:),ALLOCATABLE   :: ZXTEM, ZYTEM, ZZTEM,  &
			           ZDXTEM, ZDYTEM, ZDZTEM

REAL,DIMENSION(:,:),ALLOCATABLE :: ZZSTEM

REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: ZSCOEF

REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: ZSZ

REAL,SAVE :: ZSH

INTEGER :: JGRIDLOOP,  KGRID,  &
	   OKZXTEM, OKZYTEM, OKZZTEM, OKZDXTEM, OKZDYTEM, OKZDZTEM,  &
	   OKZZSTEM, OKZSCOEF, OKZXMASS, OKZYMASS, OKZZMASS, OKXXZS, OKXZZ, OKZSZ
!
!-------------------------------------------------------------------------------
!
!*       1.    ARRAY DIMENSIONS SETTOING
!              -------------------------
!
if(nverbia > 0)then
if (LEN(CDIRCUR) .LT. 500)THEN
print *,' **COMPCOORD KGRID DIRCUR ',KGRID,CDIRCUR(1:LEN_TRIM(CDIRCUR))
endif
endif
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE == 'SU')THEN
  IKU=1
ENDIF
IIB=1+JPHEXT
IJB=1+JPHEXT
IKB=1+JPVEXT
IIE=IIU-JPHEXT
IJE=IJU-JPHEXT
IKE=IKU-JPVEXT
!
!-------------------------------------------------------------------------------
!
!*       2.   CALCULATIONS PERFORMED FOR THE FIRST CALL
!              ----------------------------------------
!
! Test on KGRID selects processing mode:
! . KGRID=0 --> X, Y, Z + meshsizes + topography computed for ALL the
!   possible grid geometry.
! . 0<KGRID<8 --> true altitude computed for the KGRID gridpoints

IF(KGRID==0)THEN                                   ! if  "KGRID=1" selected

!
!*	 2.1   Array allocation when called for the first time
!
!   1D Arrays
!
  ALLOCATE(ZXTEM(1:IIU),STAT=OKZXTEM) ! RINT *,' OKZXTEM',OKZXTEM
  ALLOCATE(ZYTEM(1:IJU),STAT=OKZYTEM) !PRINT *,' OKZYTEM',OKZYTEM
  ALLOCATE(ZZTEM(1:IKU),STAT=OKZZTEM) !PRINT *,' OKZZTEM',OKZZTEM

  ALLOCATE(ZDXTEM(1:IIU),STAT=OKZDXTEM) !PRINT *,' OKZDXTEM',OKZDXTEM
  ALLOCATE(ZDYTEM(1:IJU),STAT=OKZDYTEM) !PRINT *,' OKZDYTEM',OKZDYTEM
  ALLOCATE(ZDZTEM(1:IKU),STAT=OKZDZTEM) !PRINT *,' OKZDZTEM',OKZDZTEM

  ALLOCATE(ZXMASS(1:IIU),STAT=OKZXMASS) !PRINT *,' OKZXMASS',OKZXMASS
  ALLOCATE(ZYMASS(1:IJU),STAT=OKZYMASS) !PRINT *,' OKZYMASS',OKZYMASS
  ALLOCATE(ZZMASS(1:IKU),STAT=OKZZMASS) !PRINT *,' OKZZMASS',OKZZMASS

!   2D Arrays
!
  ALLOCATE(ZZSTEM(1:IIU,1:IJU),STAT=OKZZSTEM) !PRINT *,' OKZZSTEM',OKZZSTEM
!
  IF(ALLOCATED(ZSCOEF))THEN
    DEALLOCATE(ZSCOEF)
  END IF
    ALLOCATE(ZSCOEF(1:IIU,1:IJU),STAT=OKZSCOEF) !PRINT *,' OKZSCOEF',OKZSCOEF
  IF(ALLOCATED(XXX))THEN
    DEALLOCATE(XXX)
  END IF
    ALLOCATE(XXX(1:IIU,7))
  IF(ALLOCATED(XXY))THEN
    DEALLOCATE(XXY)
  END IF
    ALLOCATE(XXY(1:IJU,7))
  IF(ALLOCATED(XXZ))THEN
    DEALLOCATE(XXZ)
  END IF
    ALLOCATE(XXZ(1:IKU,7))
  IF(ALLOCATED(XXDXHAT))THEN
    DEALLOCATE(XXDXHAT)
  END IF
    ALLOCATE(XXDXHAT(1:IIU,7))
  IF(ALLOCATED(XXDYHAT))THEN
    DEALLOCATE(XXDYHAT)
  END IF
    ALLOCATE(XXDYHAT(1:IJU,7))

!   3D Arrays
!
  IF(ALLOCATED(XXZS))THEN
    DEALLOCATE(XXZS)
  END IF
    ALLOCATE(XXZS(1:IIU,1:IJU,7),STAT=OKXXZS) !PRINT *,' OKXXZS',OKXXZS
  IF(ALLOCATED(ZSZ))THEN
    DEALLOCATE(ZSZ)
  END IF
    ALLOCATE(ZSZ(1:IIU,1:IJU,IKU),STAT=OKZSZ) !PRINT *,' OKZSZ',OKZSZ
!
!*       2.2   Computes true altitudes on the W grid (KGRID=4)
!

IF(CSTORAGE_TYPE /= 'PG' .AND. CSTORAGE_TYPE /='SU')THEN
  print *,' ******* COMPCOORD_FORDIACHRO  ZHAT(IKE+1) ',XZHAT(IKE+1)
  CALL VERT_COORD(LSLEVE,XZS,XZSMT,XLEN1,XLEN2,XZHAT,ZSZ)
ENDIF
!
!*       2.3    Interpolates XHAT, YHAT, ZHAT at mass gridpoints
!

  ZXMASS(1:IIU-1)=.5*(XXHAT(2:IIU)+XXHAT(1:IIU-1))
  ZXMASS(IIU)=2.*ZXMASS(IIU-1)-ZXMASS(IIU-2)
  ZYMASS(1:IJU-1)=.5*(XYHAT(2:IJU)+XYHAT(1:IJU-1))
  ZYMASS(IJU)=2.*ZYMASS(IJU-1)-ZYMASS(IJU-2)
  IF(IKU == 1)THEN
    ZZMASS(1:IKU)=XZHAT(1:IKU)
    !ZZMASS(1:IKU)=.5*(XZHAT(2:IKU)+XZHAT(1:IKU-1))  !!! size(XZHAT)=1 !!!
  ELSE
    ZZMASS(1:IKU-1)=.5*(XZHAT(2:IKU)+XZHAT(1:IKU-1))
    ZZMASS(IKU)=2.*ZZMASS(IKU-1)-ZZMASS(IKU-2)
  ENDIF

!
!*       2.4    Interpolates X, Y, Z, meshsizes, and topography 
!*              for all the KGRID selection  locations
!
  DO JGRIDLOOP=1,7

    SELECT CASE(JGRIDLOOP)
   
      CASE(1)
        ZXTEM(:)=ZXMASS(:)
        ZYTEM(:)=ZYMASS(:)
        ZZTEM(:)=ZZMASS(:)
        ZZSTEM(:,:)=XZS(:,:)
  
      CASE(2)
        ZXTEM(:)=XXHAT(:)
        ZYTEM(:)=ZYMASS(:)
        ZZTEM(:)=ZZMASS(:)
        ZZSTEM(2:IIU,:)=.5*(XZS(2:IIU,:)+XZS(1:IIU-1,:))
        ZZSTEM(1,:)=XZS(1,:)
  
      CASE(3)
        ZXTEM(:)=ZXMASS(:)
        ZYTEM(:)=XYHAT(:)
        ZZTEM(:)=ZZMASS(:)
        ZZSTEM(:,2:IJU)=.5*(XZS(:,2:IJU)+XZS(:,1:IJU-1))
        ZZSTEM(:,1)=XZS(:,1)
   
      CASE(4)
        ZXTEM(:)=ZXMASS(:)
        ZYTEM(:)=ZYMASS(:)
        ZZTEM(:)=XZHAT(:)
        ZZSTEM(:,:)=XZS(:,:)
   
      CASE(5)
        ZXTEM(:)=XXHAT(:)
        ZYTEM(:)=XYHAT(:)
        ZZTEM(:)=ZZMASS(:)
        ZZSTEM(2:IIU,:)=.5*(XZS(2:IIU,:)+XZS(1:IIU-1,:))
        ZZSTEM(1,:)=XZS(1,:)
        ZZSTEM(:,2:IJU)=.5*(ZZSTEM(:,2:IJU)+ZZSTEM(:,1:IJU-1))
        ZZSTEM(:,1)=ZZSTEM(:,2)
  
      CASE(6)
        ZXTEM(:)=XXHAT(:)
        ZYTEM(:)=ZYMASS(:)
        ZZTEM(:)=XZHAT(:)
        ZZSTEM(2:IIU,:)=.5*(XZS(2:IIU,:)+XZS(1:IIU-1,:))
        ZZSTEM(1,:)=XZS(1,:)
  
      CASE(7)
        ZXTEM(:)=ZXMASS(:)
        ZYTEM(:)=XYHAT(:)
        ZZTEM(:)=XZHAT(:)
        ZZSTEM(:,2:IJU)=.5*(XZS(:,2:IJU)+XZS(:,1:IJU-1))
        ZZSTEM(:,1)=XZS(:,1)
  
    END SELECT
   
    ZDXTEM(1:IIU-1)=ZXTEM(2:IIU)-ZXTEM(1:IIU-1)
!
! NOTICE: An extra meshlength is added to the max size of the arrays
! in order to avoid a lot of testing hereafter...
!
    ZDXTEM(IIU)=ZDXTEM(IIU-1)
    ZDYTEM(1:IJU-1)=ZYTEM(2:IJU)-ZYTEM(1:IJU-1)
    ZDYTEM(IJU)=ZDYTEM(IJU-1)
    IF(IKU /= 1)THEN
    ZDZTEM(1:IKU-1)=ZZTEM(2:IKU)-ZZTEM(1:IKU-1)
    ZDZTEM(IKU)=ZDZTEM(IKU-1)
    ENDIF
  
! X, Y, Z as functions of KGRID
    XXX(:,JGRIDLOOP)=ZXTEM
    XXY(:,JGRIDLOOP)=ZYTEM
    XXZ(:,JGRIDLOOP)=ZZTEM
  
! Topography as a function of KGRID
    XXZS(:,:,JGRIDLOOP)=ZZSTEM
   
! Meshsizes as functions of KGRID
    XXDXHAT(:,JGRIDLOOP)=ZDXTEM(:)
    XXDYHAT(:,JGRIDLOOP)=ZDYTEM(:)
  
  ENDDO
   
  DEALLOCATE(ZXMASS,ZYMASS,ZZMASS)
  DEALLOCATE(ZXTEM,ZYTEM,ZZTEM)
  DEALLOCATE(ZDXTEM,ZDYTEM,ZDZTEM)
  DEALLOCATE(ZZSTEM)

!-------------------------------------------------------------------------------
!
!*       3.   CALCULATIONS PERFORMED FOR ALL SUBSEQUENT CALLS
!             -----------------------------------------------
!
ELSE                                            ! else if KGRID =/=1 selected
   
! True altitudes

  SELECT CASE(KGRID)
   
    CASE(1)
      XZZ(:,:,1:IKU-1)=0.5*(ZSZ(:,:,1:IKU-1)+ZSZ(:,:,2:IKU))
      XZZ(:,:,IKU)=2.*XZZ(:,:,IKU-1)-XZZ(:,:,IKU-2)
   
    CASE(2)
      XZZ(:,:,1:IKU-1)=0.5*(ZSZ(:,:,1:IKU-1)+ZSZ(:,:,2:IKU))
      XZZ(:,:,IKU)=2.*XZZ(:,:,IKU-1)-XZZ(:,:,IKU-2)
      XZZ(2:IIU,:,:)=0.5*(XZZ(2:IIU,:,:)+XZZ(1:IIU-1,:,:))
      XZZ(1,:,:)=2*XZZ(2,:,:)-XZZ(3,:,:)

    CASE(3)
      XZZ(:,:,1:IKU-1)=0.5*(ZSZ(:,:,1:IKU-1)+ZSZ(:,:,2:IKU))
      XZZ(:,:,IKU)=2.*XZZ(:,:,IKU-1)-XZZ(:,:,IKU-2)
      XZZ(:,2:IJU,:)=0.5*(XZZ(:,2:IJU,:)+XZZ(:,1:IJU-1,:))
      XZZ(:,1,:)=2*XZZ(:,2,:)-XZZ(:,3,:)

    CASE(4)
      XZZ(:,:,:)=ZSZ(:,:,:)

    CASE(5)
      XZZ(:,:,1:IKU-1)=0.5*(ZSZ(:,:,1:IKU-1)+ZSZ(:,:,2:IKU))
      XZZ(:,:,IKU)=2.*XZZ(:,:,IKU-1)-XZZ(:,:,IKU-2)
      XZZ(2:IIU,:,:)=0.5*(XZZ(2:IIU,:,:)+XZZ(1:IIU-1,:,:))
      XZZ(1,:,:)=2*XZZ(2,:,:)-XZZ(3,:,:)
      XZZ(:,2:IJU,:)=0.5*(XZZ(:,2:IJU,:)+XZZ(:,1:IJU-1,:))
      XZZ(:,1,:)=2*XZZ(:,2,:)-XZZ(:,3,:)

    CASE(6)
      XZZ(2:IIU,:,:)=0.5*(ZSZ(2:IIU,:,:)+ZSZ(1:IIU-1,:,:))
      XZZ(1,:,:)=2*XZZ(2,:,:)-XZZ(3,:,:)

    CASE(7)
      XZZ(:,2:IJU,:)=0.5*(ZSZ(:,2:IJU,:)+ZSZ(:,1:IJU-1,:))
      XZZ(:,1,:)=2*XZZ(:,2,:)-XZZ(:,3,:)

  END SELECT

END IF                                                ! End KGRID selection
!
!---------------------------------------------------------------------------
!
!*       4.   EXIT
!             ----
!
RETURN
END SUBROUTINE COMPCOORD_FORDIACHRO
