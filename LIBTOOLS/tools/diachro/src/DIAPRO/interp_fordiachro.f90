!     ######spl
      MODULE MODI_INTERP_FORDIACHRO
!     #############################
!
INTERFACE
!
SUBROUTINE INTERP_FORDIACHRO(KLREF,KD,KF,PTAB,PTABREF)
REAL,DIMENSION(:,:,:), INTENT(IN)         :: PTAB    
REAL,DIMENSION(SIZE(PTAB,1),SIZE(PTAB,2)) :: PTABREF
INTEGER          :: KLREF
INTEGER          :: KD, KF
END SUBROUTINE INTERP_FORDIACHRO
!
END INTERFACE
!
END MODULE MODI_INTERP_FORDIACHRO
!     ######spl
      SUBROUTINE INTERP_FORDIACHRO(KLREF,KD,KF,PTAB,PTABREF)
!     ######################################################
!
!!****  *INTERP_FORDIACHRO* - Horizontal cross-section interpolation
!!
!!    PURPOSE
!!    -------
!       Interpolates 2D horizontal cross-sections within the Meso-NH 3D
!     arrays. These horizontal sections can be:
!     -> constant model-level sections (no interpolation, only sampling
!        of a particular level);
!     -> constant Z (sea-level  altitude) sections;
!     -> constant P (hydrostatic pressure) sections 
!     -> isentropic (constant potential temperature) 
!                                           sections
!
!!**  METHOD
!!    ------
!!      
!!      Linear interpolation of the model field with
!!    respect to "height"  when required
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist (former NCAR common)
!!         CTYPHOR :  Horizontal cross-section type
!!                    (='K' --> model level section;
!!                     ='Z' --> constant-altitude section;
!!                     ='P' --> isobar section (planned)
!!                     ='T' --> isentrope section (planned))
!!         XSPVAL  : Special value
!!
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist (former PARA common)
!!         Module MODD_DIM1 : contains dimensions of data arrays
!!                               NKMAX       : z array dimension
!!                               NIINF, NISUP: lower and upper bounds of arrays
!!                                             to be plotted in x direction
!!                               NJINF, NJSUP: lower and upper bounds of arrays
!!                                             to be plotted in y direction
!! 
!!      Module MODD_PARAMETERS : Contains array border depths
!!          JPHEXT : Horizontal external points number
!!          JPVEXT : Vertical external points number
!!         
!!      Module MODD_GRID1      : declares grid variables (Model module)
!!          XZZ    : true gridpoint z altitude
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
!!      Original       06/06/94
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODN_NCAR
USE MODN_PARA           !NOTICE: MODN_PARA includes MODD_DIM1
USE MODD_PARAMETERS
USE MODD_MASK3D
USE MODD_GRID1
USE MODD_TYPE_AND_LH
!  07/08/96  !
USE MODD_NMGRID
USE MODD_PT_FOR_CH_FORDIACHRO
!  07/08/96  !
USE MODD_RESOLVCAR

IMPLICIT NONE
INTERFACE
  SUBROUTINE COMPLAT(PLAT)
  REAL,DIMENSION(:,:) :: PLAT
  END SUBROUTINE
END INTERFACE
!
!*       0.1   Declaration of arguments and results
!
REAL,DIMENSION(:,:,:), INTENT(IN)         :: PTAB    !Input arrays where the 
                                                     !horizontal section is cut 
REAL,DIMENSION(SIZE(PTAB,1),SIZE(PTAB,2)) :: PTABREF !Output array containing                                                        !the sampled plane
INTEGER :: KLREF  !Sampled level location:
                  !If CTYPHOR='K'-> model level index given,
                  !If CTYPHOR='Z'-> sea-level altitude given in meters,
                  !If CTYPHOR='P'-> pressure level given in hPa,
                  !If CTYPHOR='T'-> potential temperature level given in K. 
INTEGER :: KD, KF ! K Bounds
!
!*       0.2   Local Variables 
!
INTEGER :: IID,IJD
INTEGER :: II, IJ, IK
INTEGER :: JILOOP, JJLOOP, JKLOOP, IKB, IKE
REAL    :: ZREF, ZDIXEPS, ZXM, ZXP
!  07/08/96  !
INTEGER :: IIUP,IJUP,IKU
INTEGER :: IND1, IND2
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: ZPTH, ZPTHPROV
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZLAT
!  07/08/96  !
!
!-------------------------------------------------------------------------------
!
!*      1.    Preliminary calculations
!             ------------------------
!
IKB=1+JPVEXT
ZDIXEPS=10.*EPSILON(1.)
!  07/08/96  !
IKU=NKMAX+2*JPVEXT
!IKU=SIZE(PTAB,3)
IKE=IKU-JPVEXT
IIUP=NIMAX+2*JPHEXT
!IIUP=SIZE(PTAB,1)
!IJUP=SIZE(PTAB,2)
IJUP=NJMAX+2*JPHEXT
!print *,' IIUP,IJUP,IKU ',IIUP,IJUP,IKU
!print *,' INTERP_FORDIACHRO SIZE(XPRES,1),SIZE(XPRES,2),SIZE(XPRES,3) ',SIZE(XPRES,1),SIZE(XPRES,2),SIZE(XPRES,3)
!  07/08/96  !
IF(LPR)THEN
  IF(ALLOCATED(ZPTH)) DEALLOCATE(ZPTH)
  ALLOCATE(ZPTH(SIZE(XPRES,1),SIZE(XPRES,2),SIZE(XPRES,3)))
ENDIF
IF(LTK .OR. LEV .OR. LSV3)THEN
  IF(ALLOCATED(ZPTH)) DEALLOCATE(ZPTH)
  ALLOCATE(ZPTH(SIZE(XTH,1),SIZE(XTH,2),SIZE(XTH,3)))
ENDIF
if(nverbia > 0)then
print *,' INTERP_FORDIACHRO LPR,LTK,LEV,LSV3 ',LPR,LTK,LEV,LSV3
endif
!
! If not a model level request, convert KLREF to 
! the appropriate variable for interpolation
!
print *,' *** Interp KLREF, XLOOPZ ',KLREF,XLOOPZ
IF(CTYPHOR.EQ.'P')THEN                        ! 'P' requested
!>>>>>>>>>>>>>YET TO BE COMPLETED
!Mars 2000
  IF(LCHREEL)THEN
    ZREF=ALOG10(XLOOPZ*100.)
  ELSE
!Mars 2000
!  07/08/96  !
    ZREF=ALOG10(FLOAT(KLREF)*100.)
!  07/08/96  !
!Mars 2000
  ENDIF
!Mars 2000
ELSE                                          ! 'Z' requested
!Mars 2000
  IF(LCHREEL)THEN
    ZREF=XLOOPZ
  ELSE
!Mars 2000
    ZREF=FLOAT(KLREF)
!Mars 2000
  ENDIF
!Mars 2000
END IF
!
!-------------------------------------------------------------------------------
!
!*      2.    Sampling of the requested horizontal section
!             --------------------------------------------
!
!*      2.1   Sampling of a model level: no interpolation necessary
!
CALL CPSETC('CFT','CONSTANT FIELD - VALUE IS $ZDV$')
IF(CTYPHOR.EQ.'K')THEN
  if(nverbia >0)then
   print *, ' ** INTERP CTYPHOR.EQ. K, KLREF,KD,KF ', KLREF,KD,KF
  endif
  IF(LMSKTOP)THEN
  if(nverbia >0)then
  print *,' INTERP MSKTOP NLOOPT ',NLOOPT
  endif
    DO JILOOP=NIINF,NISUP
      DO JJLOOP=NJINF,NJSUP
	DO JKLOOP=KF,KD,-1
	  IID=JILOOP-NIINF+1
	  IJD=JJLOOP-NJINF+1
	    PTABREF(IID,IJD)=XSPVAL
	    IF(LMASK3(JILOOP,JJLOOP,JKLOOP,NLOOPT))THEN
	      PTABREF(IID,IJD)=PTAB(IID,IJD,JKLOOP-KD+1)
	      EXIT
	    ENDIF
	ENDDO
      ENDDO
    ENDDO
  if(nverbia >0)then
     print *,' ** interp CTYPHOR=K AV RETURN DANS LMSKTOP, LPR=',LPR
  endif
    RETURN
  ELSE
  if(nverbia >0)then
   print *, ' ** INTERP AV IF(KLREF < KD)THEN KLREF,KD,KF ', KLREF,KD,KF
  endif
  IF(KLREF < KD)THEN
! Ajout LKCP Avril 2001 -> prise en compte bilans compresses en K
  IF(LKCP)THEN
    PTABREF(:,:)=PTAB(:,:,1)
  ELSE
  CALL CPSETC('CFT','UNDER IKB or 1st recorded LEVEL')
  PTABREF(:,:)=XSPVAL
  ENDIF
  RETURN
  ELSE IF(KLREF > KF)THEN
  CALL CPSETC('CFT','OVER IKE or last recorded LEVEL')
  PTABREF(:,:)=XSPVAL
  RETURN
  ELSE
  PTABREF(:,:)=PTAB(:,:,KLREF-KD+1)
  IND1=0
  DO JILOOP=1,SIZE(PTABREF,1)
  DO JJLOOP=1,SIZE(PTABREF,2)
    IF(PTABREF(JILOOP,JJLOOP) /= XSPVAL)THEN
      IND1=1
      EXIT
    ENDIF
  ENDDO
  ENDDO
  IF(IND1 == 0)THEN
    CALL CPSETC('CFT','No Value')
    IND1=0
  ENDIF
  RETURN
  ENDIF
  ENDIF
END IF
!  07/08/96  !
IF(CTYPHOR.EQ.'P')THEN
    ZPTH=XPRES(:,:,:,NLOOPT,1,1)
ENDIF
IF(CTYPHOR.EQ.'T' .OR. CTYPHOR.EQ.'E' .OR. CTYPHOR.EQ.'V')THEN
    IF(CTYPHOR.EQ.'E')THEN
      II=SIZE(XTH,1)
      IJ=SIZE(XTH,2)
      ALLOCATE(ZLAT(II,IJ))
      CALL COMPLAT(ZLAT)
      IK=SIZE(XTH,3)
! 7 Mars 2000
      print *,' interpol *** NLOOPT ',NLOOPT
      ZPTH=XTH(:,:,:,NLOOPT,1,1)
      DO JKLOOP=1,IK
	WHERE(ZPTH(:,:,JKLOOP) /= XSPVAL)
	  ZPTH(:,:,JKLOOP)=ZPTH(:,:,JKLOOP)* &
	  SIGN(1.,ZLAT(:,:))
	ENDWHERE
!       WHERE(XTH(:,:,JKLOOP,NLOOPT,1,1) /= XSPVAL)
!         XTH(:,:,JKLOOP,NLOOPT,1,1)=XTH(:,:,JKLOOP,NLOOPT,1,1)* &
!         SIGN(1.,ZLAT(:,:))
!       ENDWHERE
      ENDDO
      DEALLOCATE(ZLAT)
    ELSE
      ZPTH=XTH(:,:,:,NLOOPT,1,1)
    ENDIF
!   ZPTH=XTH(:,:,:,NLOOPT,1,1)
ENDIF
IF(CTYPHOR == 'T' .OR. CTYPHOR == 'P' .OR. CTYPHOR == 'E' .OR. CTYPHOR.EQ.'V')THEN
!Mars 2000
IF(LSV3 .OR. LXYZ)THEN
  IF(ALLOCATED(ZPTHPROV))THEN
    DEALLOCATE(ZPTHPROV)
  ENDIF
  ALLOCATE(ZPTHPROV(SIZE(ZPTH,1),SIZE(ZPTH,2),SIZE(ZPTH,3)))
  ZPTHPROV=XSPVAL
SELECT CASE (NMGRID)
  CASE(1)
  CASE(2)
    WHERE(ZPTH(2:IIUP,:,:) /= XSPVAL .AND. ZPTH(1:IIUP-1,:,:) /= XSPVAL)
      ZPTHPROV(2:IIUP,:,:)=.5*(ZPTH(2:IIUP,:,:) + ZPTH(1:IIUP-1,:,:))
    ENDWHERE
    WHERE(ZPTHPROV(2,:,:) /= XSPVAL .AND. ZPTHPROV(3,:,:) /= XSPVAL)
      ZPTHPROV(1,:,:)=2.*ZPTHPROV(2,:,:)-ZPTHPROV(3,:,:)
    ENDWHERE
    ZPTH=ZPTHPROV
  CASE(3)
    WHERE(ZPTH(:,2:IJUP,:) /= XSPVAL .AND. ZPTH(:,1:IJUP-1,:) /= XSPVAL)
      ZPTHPROV(:,2:IJUP,:)=.5*(ZPTH(:,2:IJUP,:) + ZPTH(:,1:IJUP-1,:))
    ENDWHERE
    WHERE(ZPTHPROV(:,2,:) /= XSPVAL .AND. ZPTHPROV(:,3,:) /= XSPVAL)
      ZPTHPROV(:,1,:)=2.*ZPTHPROV(:,2,:)-ZPTHPROV(:,3,:)
    ENDWHERE
    ZPTH=ZPTHPROV
  CASE(4)
    WHERE(ZPTH(:,:,2:IKU) /= XSPVAL .AND. ZPTH(:,:,1:IKU-1) /= XSPVAL)
      ZPTHPROV(:,:,2:IKU)=.5*(ZPTH(:,:,2:IKU) + ZPTH(:,:,1:IKU-1))
    ENDWHERE
    WHERE(ZPTHPROV(:,:,2) /= XSPVAL .AND. ZPTHPROV(:,:,3) /= XSPVAL)
      ZPTHPROV(:,:,1)=2.*ZPTHPROV(:,:,2)-ZPTHPROV(:,:,3)
    ENDWHERE
    ZPTH=ZPTHPROV
  CASE(5)
    WHERE(ZPTH(2:IIUP,:,:) /= XSPVAL .AND. ZPTH(1:IIUP-1,:,:) /= XSPVAL)
      ZPTHPROV(2:IIUP,:,:)=.5*(ZPTH(2:IIUP,:,:) + ZPTH(1:IIUP-1,:,:))
    ENDWHERE
    WHERE(ZPTHPROV(2,:,:) /= XSPVAL .AND. ZPTHPROV(3,:,:) /= XSPVAL)
      ZPTHPROV(1,:,:)=2.*ZPTHPROV(2,:,:)-ZPTHPROV(3,:,:)
    ENDWHERE
    ZPTH=ZPTHPROV
    ZPTHPROV=XSPVAL
    WHERE(ZPTH(:,2:IJUP,:) /= XSPVAL .AND. ZPTH(:,1:IJUP-1,:) /= XSPVAL)
      ZPTHPROV(:,2:IJUP,:)=.5*(ZPTH(:,2:IJUP,:) + ZPTH(:,1:IJUP-1,:))
    ENDWHERE
    WHERE(ZPTHPROV(:,2,:) /= XSPVAL .AND. ZPTHPROV(:,3,:) /= XSPVAL)
      ZPTHPROV(:,1,:)=2.*ZPTHPROV(:,2,:)-ZPTHPROV(:,3,:)
    ENDWHERE
    ZPTH=ZPTHPROV
  CASE(6)
    WHERE(ZPTH(:,:,2:IKU) /= XSPVAL .AND. ZPTH(:,:,1:IKU-1) /= XSPVAL)
      ZPTHPROV(:,:,2:IKU)=.5*(ZPTH(:,:,2:IKU) + ZPTH(:,:,1:IKU-1))
    ENDWHERE
    WHERE(ZPTHPROV(:,:,2) /= XSPVAL .AND. ZPTHPROV(:,:,3) /= XSPVAL)
      ZPTHPROV(:,:,1)=2.*ZPTHPROV(:,:,2)-ZPTHPROV(:,:,3)
    ENDWHERE
    ZPTH=ZPTHPROV
    ZPTHPROV=XSPVAL
    WHERE(ZPTH(2:IIUP,:,:) /= XSPVAL .AND. ZPTH(1:IIUP-1,:,:) /= XSPVAL)
      ZPTHPROV(2:IIUP,:,:)=.5*(ZPTH(2:IIUP,:,:) + ZPTH(1:IIUP-1,:,:))
    ENDWHERE
    WHERE(ZPTHPROV(2,:,:) /= XSPVAL .AND. ZPTHPROV(3,:,:) /= XSPVAL)
      ZPTHPROV(1,:,:)=2.*ZPTHPROV(2,:,:)-ZPTHPROV(3,:,:)
    ENDWHERE
    ZPTH=ZPTHPROV
  CASE(7)
    WHERE(ZPTH(:,:,2:IKU) /= XSPVAL .AND. ZPTH(:,:,1:IKU-1) /= XSPVAL)
      ZPTHPROV(:,:,2:IKU)=.5*(ZPTH(:,:,2:IKU) + ZPTH(:,:,1:IKU-1))
    ENDWHERE
    WHERE(ZPTHPROV(:,:,2) /= XSPVAL .AND. ZPTHPROV(:,:,3) /= XSPVAL)
      ZPTHPROV(:,:,1)=2.*ZPTHPROV(:,:,2)-ZPTHPROV(:,:,3)
    ENDWHERE
    ZPTH=ZPTHPROV
    ZPTHPROV=XSPVAL
    WHERE(ZPTH(:,2:IJUP,:) /= XSPVAL .AND. ZPTH(:,1:IJUP-1,:) /= XSPVAL)
      ZPTHPROV(:,2:IJUP,:)=.5*(ZPTH(:,2:IJUP,:) + ZPTH(:,1:IJUP-1,:))
    ENDWHERE
    WHERE(ZPTHPROV(:,2,:) /= XSPVAL .AND. ZPTHPROV(:,3,:) /= XSPVAL)
      ZPTHPROV(:,1,:)=2.*ZPTHPROV(:,2,:)-ZPTHPROV(:,3,:)
    ENDWHERE
    ZPTH=ZPTHPROV
END SELECT
DEALLOCATE(ZPTHPROV)

ELSE
!Mars 2000

SELECT CASE (NMGRID)
  CASE(1)
  CASE(2)
    ZPTH(2:IIUP,:,:)=.5*(ZPTH(2:IIUP,:,:) + ZPTH(1:IIUP-1,:,:))
    ZPTH(1,:,:)=2.*ZPTH(2,:,:)-ZPTH(3,:,:)
  CASE(3)
    ZPTH(:,2:IJUP,:)=.5*(ZPTH(:,2:IJUP,:) + ZPTH(:,1:IJUP-1,:))
    ZPTH(:,1,:)=2.*ZPTH(:,2,:)-ZPTH(:,3,:)
  CASE(4)
    ZPTH(:,:,2:IKU)=.5*(ZPTH(:,:,2:IKU) + ZPTH(:,:,1:IKU-1))
    ZPTH(:,:,1)=2.*ZPTH(:,:,2)-ZPTH(:,:,3)
  CASE(5)
    ZPTH(2:IIUP,:,:)=.5*(ZPTH(2:IIUP,:,:) + ZPTH(1:IIUP-1,:,:))
    ZPTH(1,:,:)=2.*ZPTH(2,:,:)-ZPTH(3,:,:)
    ZPTH(:,2:IJUP,:)=.5*(ZPTH(:,2:IJUP,:) + ZPTH(:,1:IJUP-1,:))
    ZPTH(:,1,:)=2.*ZPTH(:,2,:)-ZPTH(:,3,:)
  CASE(6)
    ZPTH(:,:,2:IKU)=.5*(ZPTH(:,:,2:IKU) + ZPTH(:,:,1:IKU-1))
    ZPTH(:,:,1)=2.*ZPTH(:,:,2)-ZPTH(:,:,3)
    ZPTH(2:IIUP,:,:)=.5*(ZPTH(2:IIUP,:,:) + ZPTH(1:IIUP-1,:,:))
    ZPTH(1,:,:)=2.*ZPTH(2,:,:)-ZPTH(3,:,:)
  CASE(7)
    ZPTH(:,:,2:IKU)=.5*(ZPTH(:,:,2:IKU) + ZPTH(:,:,1:IKU-1))
    ZPTH(:,:,1)=2.*ZPTH(:,:,2)-ZPTH(:,:,3)
    ZPTH(:,2:IJUP,:)=.5*(ZPTH(:,2:IJUP,:) + ZPTH(:,1:IJUP-1,:))
    ZPTH(:,1,:)=2.*ZPTH(:,2,:)-ZPTH(:,3,:)
END SELECT

!Mars 2000
ENDIF
!Mars 2000
!IF(CTYPHOR == 'P')print *,' ZPTH AP MISE SUR GRILLE '
ENDIF
!  07/08/96  !
!
!*      2.2   Not a model level request: interpolation necessary
!
DO JILOOP=NIINF,NISUP
  DO JJLOOP=NJINF,NJSUP

    IF((CTYPHOR.EQ.'E' .OR. CTYPHOR.EQ.'V') .AND. LINTERPTOP)THEN
    DO JKLOOP=KF,KD,-1
!
      IID=JILOOP-NIINF+1
      IJD=JJLOOP-NJINF+1
!
!*     2.2.3  Potential vorticity request: prepares PV interpolation
!
!  07/08/96  !
        ZXM=ZPTH(JILOOP,JJLOOP,JKLOOP)
        ZXP=ZPTH(JILOOP,JJLOOP,MIN(KF,JKLOOP+1))
!  07/08/96  !
!
!*     2.3    Selects points within the TRACE display window
!
!
!  07/08/96  !
      PTABREF(IID,IJD)=XSPVAL
!  18/02/2000 Essai pour prise en compte des valeurs speciales
      IF(LSV3 .AND. LXYZ00)THEN
	IF(ZXP == XSPVAL .OR. ZXM  == XSPVAL .OR. (ZXP == XSPVAL .AND. &
	  ZXM  == XSPVAL))THEN
	  if(nverbia == 20)then
	  print *,' ***interp JILOOP JJLOOP JKLOOP ZXP ZXM ',JILOOP,&
	  JJLOOP,JKLOOP,ZXP,ZXM
	  endif
	  CYCLE
	ENDIF
      ENDIF
!  18/02/2000 Essai pour prise en compte des valeurs speciales
      IF((ZXP-ZREF)*(ZREF-ZXM).GE.0.)THEN
        IF(JKLOOP+1 <= IKB .OR. JKLOOP+1 > IKE)THEN
          CYCLE
        ELSE
          GO TO 4
        ENDIF
      ELSE IF(ZXP.GE.ZXM-ZDIXEPS.AND.ZXP.LE.ZXM+ZDIXEPS.AND.  &
      ZREF.GE.ZXM-ZDIXEPS.AND.ZREF.LE.ZXM+ZDIXEPS)THEN
        IF(JKLOOP+1 <= IKB .OR. JKLOOP+1 > IKE)THEN
          CYCLE
        ELSE
          GO TO 4
        ENDIF
      ENDIF
!  07/08/96  !
!
    ENDDO

    ELSE

    DO JKLOOP=KD,KF
!
      IID=JILOOP-NIINF+1
      IJD=JJLOOP-NJINF+1
!
!*     2.2.1  Pressure level request: prepares Log(P) interpolation
!
      IF(CTYPHOR.EQ.'P')THEN
!>>>>>>>>>>>>YET TO BE DEVELOPED
        ZXM=ALOG10(ZPTH(JILOOP,JJLOOP,JKLOOP))
        ZXP=ALOG10(ZPTH(JILOOP,JJLOOP,MIN(KF,JKLOOP+1)))
!
!*     2.2.2  Altitude level request: prepares Z interpolation
!
      ELSE IF (CTYPHOR.EQ.'Z')THEN
        ZXM=XZZ(JILOOP,JJLOOP,JKLOOP)
        ZXP=XZZ(JILOOP,JJLOOP,MIN(KF,JKLOOP+1))
!
!*     2.2.3  Potential temperature request: prepares Theta interpolation
!
      ELSE IF(CTYPHOR.EQ.'T')THEN
!>>>>>>>>>>>>YET TO BE DEVELOPED
!  07/08/96  !
        ZXM=ZPTH(JILOOP,JJLOOP,JKLOOP)
        ZXP=ZPTH(JILOOP,JJLOOP,MIN(KF,JKLOOP+1))
!  07/08/96  !
! Mars 2000 Ajout possibilite de faire interpolation a partir du bas
! pour la vorticite potentielle et SV3
      ELSE IF(CTYPHOR.EQ.'E' .AND. .NOT.LINTERPTOP)THEN
        ZXM=ZPTH(JILOOP,JJLOOP,JKLOOP)
        ZXP=ZPTH(JILOOP,JJLOOP,MIN(KF,JKLOOP+1))
      ELSE IF(CTYPHOR.EQ.'V' .AND. .NOT.LINTERPTOP)THEN
        ZXM=ZPTH(JILOOP,JJLOOP,JKLOOP)
        ZXP=ZPTH(JILOOP,JJLOOP,MIN(KF,JKLOOP+1))
! Mars 2000 
      END IF
!
!*     2.3    Selects points within the TRACE display window
!
!
!  07/08/96  !
      PTABREF(IID,IJD)=XSPVAL
!  23/03/2000 Essai pour prise en compte des valeurs speciales
      IF(LSV3 .AND. LXYZ00)THEN
	IF(ZXP == XSPVAL .OR. ZXM  == XSPVAL .OR. (ZXP == XSPVAL .AND. &
	  ZXM  == XSPVAL))THEN
	  if(nverbia == 20)then
	  print *,' ***interp JILOOP JJLOOP JKLOOP ZXP ZXM ',JILOOP,&
	  JJLOOP,JKLOOP,ZXP,ZXM
	  endif
	  CYCLE
	ENDIF
      ENDIF
!  23/03/2000 Essai pour prise en compte des valeurs speciales
      IF((ZXP-ZREF)*(ZREF-ZXM).GE.0.)THEN
        IF(JKLOOP+1 <= IKB .OR. JKLOOP+1 > IKE)THEN
          CYCLE
        ELSE
          GO TO 4
        ENDIF
      ELSE IF(ZXP.GE.ZXM-ZDIXEPS.AND.ZXP.LE.ZXM+ZDIXEPS.AND.  &
      ZREF.GE.ZXM-ZDIXEPS.AND.ZREF.LE.ZXM+ZDIXEPS)THEN
        IF(JKLOOP+1 <= IKB .OR. JKLOOP+1 > IKE)THEN
          CYCLE
        ELSE
          GO TO 4
        ENDIF
      ENDIF
!  07/08/96  !
!
    ENDDO
    ENDIF
!
!*    2.4    Out of display window: inserts a NCAR special value 
!            to suppress display
!
    PTABREF(IID,IJD)=XSPVAL
    GO TO 5
!
4   CONTINUE
!
!*    2.5   Requested level colocated with a model level: no interpolation
! 
    IF(ZXP==ZXM)THEN
      PTABREF(IID,IJD)=PTAB(IID,IJD,JKLOOP-KD+1)
!       print *,' INTERP_FORDIACHRO ZXM ZXP ',ZXM,ZXP
!     IF(CTYPHOR == 'P')THEN
!       print *,' CAS ZXM=ZXP '
!     ENDIF 
!
!*    2.6   Requested level located between model levels: linear interpolation
!
    ELSE
      SELECT CASE(CTYPHOR)
        CASE('Z')
!         print *,' ZXP - ZXM ',ZXP-ZXM
          PTABREF(IID,IJD)=(PTAB(IID,IJD,JKLOOP-KD+1)*(ZXP-ZREF)+  &
          PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP+1-KD+1))*  &
          (ZREF-ZXM))/MAX(1.E-8,(ZXP-ZXM))
        CASE('T','E','V')
!         print *,' ZXP - ZXM ',ZXP-ZXM
          LTHSTAB=.TRUE.
          IF(JKLOOP+1 > IKB)THEN
            IF(ZXP-ZXM >= 0.)THEN
              LTHSTAB=.TRUE.
            ELSE
              LTHSTAB=.FALSE.
!             print *,' JKLOOP, ZXP, ZXM ',JKLOOP,ZXP,ZXM
          ENDIF
          ENDIF
          PTABREF(IID,IJD)=(PTAB(IID,IJD,JKLOOP-KD+1)*(ZXP-ZREF)+  &
          PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP+1-KD+1))*  &
          (ZREF-ZXM))/(ZXP-ZXM)
        CASE('P')
          PTABREF(IID,IJD)=(PTAB(IID,IJD,JKLOOP-KD+1)*(ZXP-ZREF)+  &
          PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP+1-KD+1))*  &
          (ZREF-ZXM))/MIN(-1.E-8,(ZXP-ZXM))
      END SELECT
!     IF(CTYPHOR == 'P' .AND. IID == 4 .AND. IJD == 8)THEN
!       print *,' IID,IJD,JKLOOP-KD+1,PTAB,ZXP-ZREF ',IID,IJD,JKLOOP-KD+1,PTAB(IID,IJD,JKLOOP-KD+1),ZXP-ZREF
!       print *,' IID,IJD,JKLOOP-KD+1,PTAB,ZXP-ZREF ZXP-ZXM',IID,IJD,JKLOOP-KD+1,PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP-KD+1+1)),ZREF-ZXM,ZXP-ZXM
!     ENDIF
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sept 2000 test suivant supprime
!   IF(LXYZ .OR. LSV3)THEN
    IF(PTAB(IID,IJD,JKLOOP-KD+1) == XSPVAL .AND. &
      PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP+1-KD+1)) == XSPVAL)THEN
      PTABREF(IID,IJD)=XSPVAL
    ELSE IF(PTAB(IID,IJD,JKLOOP-KD+1) /= XSPVAL .AND. &
      PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP+1-KD+1)) == XSPVAL)THEN
      PTABREF(IID,IJD)=XSPVAL
!     PTABREF(IID,IJD)=PTAB(IID,IJD,JKLOOP-KD+1)
    ELSE IF(PTAB(IID,IJD,JKLOOP-KD+1) == XSPVAL .AND. &
      PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP+1-KD+1)) /= XSPVAL)THEN
!     PTABREF(IID,IJD)=PTAB(IID,IJD,MIN(KF-KD+1,JKLOOP+1-KD+1))
      PTABREF(IID,IJD)=XSPVAL
    ENDIF
!   ENDIF
! Sept 2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
5 CONTINUE
!
  ENDDO
ENDDO
!
IND1=0
IND2=0
DO JILOOP=1,SIZE(PTABREF,1)
  DO JJLOOP=1,SIZE(PTABREF,2)
    IF(PTABREF(JILOOP,JJLOOP) /= XSPVAL)THEN
      IND1=1
      EXIT
    ELSE
      IND2=1
    ENDIF
  ENDDO
ENDDO
!print *,' PTABREF 1-8 '
!DO JJLOOP=1,SIZE(PTABREF,2)
!  print *,(PTABREF(JILOOP,JJLOOP),JILOOP=1,8)
!ENDDO
!print *,' PTABREF 9-16 '
!DO JJLOOP=1,SIZE(PTABREF,2)
!  print *,(PTABREF(JILOOP,JJLOOP),JILOOP=9,16)
!ENDDO
!print *,' PTABREF 17-24 '
!DO JJLOOP=1,SIZE(PTABREF,2)
!  print *,(PTABREF(JILOOP,JJLOOP),JILOOP=17,24)
!ENDDO
IF(IND1 == 0 .AND. IND2 /= 0)THEN
  CALL CPSETC('CFT','<IKB or 1st recorded LEVEL or >IKE LEVEL')
  IF(LSV3 .OR. LXYZ)THEN
    CALL CPSETC('CFT','No value')
  ENDIF
ELSE 
  CALL CPSETC('CFT','CONSTANT FIELD - VALUE IS $ZDV$')
ENDIF
if(nverbia > 0)then
  print *,' INTERP_FORDIACHRO end: LPR,LTK,LEV,LSV3 ',LPR,LTK,LEV,LSV3
endif
    
!
!----------------------------------------------------------------------------
!
!*     3.    EXIT
!            ----
!
RETURN
END SUBROUTINE INTERP_FORDIACHRO
