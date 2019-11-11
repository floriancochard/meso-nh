!     ######spl
      SUBROUTINE COUPE_FORDIACHRO(PTABI,PTABO,K)
!     ##########################################
!
!!****  *COUPE_FORDIACHRO* - Vertical cross-section interpolation
!!
!!    PURPOSE
!!    -------
!         Interpolates 2D vertical cross-sections within the Meso-NH 3D
!       arrays. Model fields, iapprpriate gridpoint altitudes as well as 
!       appropriate topography height are interpolated. 
!
!!**  METHOD
!!    ------
!!        The general case of a vertical cross-section along any oblique 
!!      direction with respect to the x-y model axes is considered. Simple
!!      linear interpolation is done.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!       XXX,XXY,XXZ    : coordinate values for all the MESO-NH grids
!!       XXZS           : topography values for all the MESO_NH grids
!!       XDSX, XDSY     : projections on the MESO-NH cartesian axes of the XDS
!!                        oblique abscissa (meters), for all grid locations
!!
!!      Module MODD_GRID1  : declares grid variables (Model module)
!!       XZZ            : true z altitude for the current NMGRID grid location
!!
!!      Module MODN_NCAR   : defines NAM_DIRTRA_POS namelist (form. NCAR common)
!!       XSPVAL         : Special value
!!
!!      Module MODD_CVERT  :  Declares work arrays for vertical cross-sections
!!       XWORKZ         : working array for true altitude storage (all grids)
!!       XWZ            : working array for topography (all grids)
!!
!!      Module MODN_PARA   : Defines NAM_DOMAIN_POS namelist (form. PARA common)
!!       NLMAX          :  Number of points horizontally along
!!                         the vertical section
!!       Module MODD_DIM1 : contains dimensions of data arrays
!!        NKMAX  : z array dimension
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!       JPHEXT         : Horizontal external points number
!!       JPVEXT         : Vertical external points number
!!
!!
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!       NMGRID         : Current MESO-NH grid indicator
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!     NCAR Graphics Technical documentation, UNIX version 3.2,
!!     Scientific computing division, NCAR/UCAR, Boulder, USA.
!!      Volume 1: Fundamentals, Vers. 1, May 1993
!!      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   06/01/59
!-------------------------------------------------------------------------
!
!*     0.   DECLARATIONS
!           ------------
!
USE MODD_COORD
USE MODD_GRID1
USE MODN_NCAR
USE MODD_CVERT
USE MODD_MEMCV
USE MODN_PARA
USE MODD_PARAMETERS
USE MODD_NMGRID 
USE MODD_RESOLVCAR
USE MODD_TYPE_AND_LH 
!
IMPLICIT NONE
!
!*     0.1  Dummy arguments and results
!
REAL,DIMENSION(:,:)      :: PTABI    ! Input data array to be interpolated
REAL,DIMENSION(:)        :: PTABO    ! Returned interpolated 2D array
INTEGER                  ::  K       ! Model level where interpolation is done
!
!*     0.2  Local variables
!
REAL    :: ZCIINF,ZCISUP,ZCJINF,ZCJSUP,ZXX,ZYY
INTEGER :: IMIM1,IMI,IMJM1,IMJ,JILOOP, JI, JJ
INTEGER :: IIU, IJU, IKU, IKB, IKE
!
!------------------------------------------------------------------------------
!
!*     1.   PERFORMING VERTICAL INTERPOLATIONS
!           ----------------------------------
!
!*     1.0  Presetting array extends
!
!print *,' ++++coupe NMGRID DIRCUR ',NMGRID,CDIRCUR(1:LEN_TRIM(CDIRCUR))
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
IKE=IKU-JPVEXT
IKB=1+JPVEXT
! Oct 2000 prise en compte du 2D horizontal -> PH=CH+CV
! Ajout LKCP Avril 2001 -> prise en compte du 3D compresse en K
IF(LCHXY .OR. LKCP)THEN
  IKU=1; IKB=1; IKE=1
  if(nverbia > 0)then
    print *,' **coupe NKL NKH LCHXY ',NKL,NKH,LCHXY
  endif
ENDIF
!
!*     1.1  Scans along-section  oblique x axis 
!
DO JILOOP=1,NLMAX   ! Start of general X- scanning loop
!
!*     1.2  Locates the current gridpoint along the cross-section 
!*          oblique x-axis within the Meso-NH grid
!
  ZXX=XDSX(JILOOP,NMGRID)   ! Collects the model X- and Y- axes projections
  ZYY=XDSY(JILOOP,NMGRID)   ! onto the oblique vertical section plane
                            ! for the current (new) section point
!
!*     1.3   The current section point is located along
!*           the x- axis
!
    DO JI=2,IIU
      IF(ZXX.LE.XXX(JI,NMGRID).AND.ZXX.GE.XXX(JI-1,NMGRID))GO TO 1
    ENDDO

1 CONTINUE

  IMIM1=MAX(1,JI-1)   
  IMI=JI             ! JI is the index of the first model bin to the left
!
!*    1.4    Then, it is located along
!*           the Y- axis
!
    DO JJ=2,IJU
      IF(ZYY.LE.XXY(JJ,NMGRID).AND.ZYY.GE.XXY(JJ-1,NMGRID))GO TO 2
    ENDDO

2 CONTINUE

  IMJM1=MAX(1,JJ-1)
  IMJ=JJ             ! JJ is the index of the first model bin below
!
!*   1.5     Finally the X- and Y- distances between the current section 
!*           point and closest model box to the left-bottom are calculated
!
  IF(IMI==IMIM1)THEN       ! Left wall special case
    ZCIINF=0.
    ZCISUP=0.
  ELSE
    !print *,'XX(IMI IMIM1) ZXX',XX(IMI),XX(IMIM1),ZXX
    ZCIINF=(XXX(IMI,NMGRID)-ZXX)/MAX(1.E-10,(XXX(IMI,NMGRID)-XXX(IMIM1,NMGRID)))
    ZCISUP=(ZXX-XXX(IMIM1,NMGRID))/MAX(1.E-10,(XXX(IMI,NMGRID)-XXX(IMIM1,NMGRID)))
  END IF
  !
  IF(IMJ==IMJM1)THEN      ! Bottom wall special case
    ZCJINF=0.
    ZCJSUP=0.
  ELSE
    !PRINT *,'XY(IMJ IMJM1) ZXY',XY(IMJ),XY(IMJM1),ZYY
    ZCJINF=(XXY(IMJ,NMGRID)-ZYY)/MAX(1.E-10,(XXY(IMJ,NMGRID)-XXY(IMJM1,NMGRID)))
    ZCJSUP=(ZYY-XXY(IMJM1,NMGRID))/MAX(1.E-10,(XXY(IMJ,NMGRID)-XXY(IMJM1,NMGRID)))
  END IF
!
!*   1.6     Computes the interpolated altitude of the 
!*           current section point
!
if(nverbia > 1)then
print *,' ** coupe   AV XWORKZ  K= ',K,' size XWORKZ ',size(XWORKZ,1),size(XWORKZ,2),size(XWORKZ,3), &
' IMIM1,IMJM1,IMJ,IMI ',IMIM1,IMJM1,IMJ,IMI
print *,' ** coupe   AV XWORKZ  ZCIINF,ZCJINF,ZCISUP,ZCJSUP,NMGRID ',ZCIINF,ZCJINF,ZCISUP,ZCJSUP,NMGRID
print *,' ** coupe   AV XWORKZ JILOOP XZZ(IMIM1,IMJM1,K),XZZ(IMIM1,IMJ,K),XZZ(IMI,IMJM1,K),XZZ(IMI,IMJ,K) ',&
JILOOP,XZZ(IMIM1,IMJM1,K),XZZ(IMIM1,IMJ,K),XZZ(IMI,IMJM1,K),XZZ(IMI,IMJ,K)
endif

  XWORKZ(JILOOP,K,NMGRID)=ZCIINF*ZCJINF*XZZ(IMIM1,IMJM1,K)+  &
         ZCIINF*ZCJSUP*XZZ(IMIM1,IMJ,K)+                     &
         ZCISUP*ZCJINF*XZZ(IMI,IMJM1,K)+                     &
         ZCISUP*ZCJSUP*XZZ(IMI,IMJ,K)    

if(nverbia > 1)then
print *,' ** coupe   AP XWORKZ  K= ',K,' size XZZ ',size(XZZ,1),size(XZZ,2),size(XZZ,3),' IMIM1,IMJM1,IMJ,IMI ',IMIM1,IMJM1,IMJ,IMI
endif
!
!*   1.7     Computes the interpolated value of the field for
!*           current section point
! 
!  Modifs for diachro
! Avril 2001 Ajout LKCP -> prise en compte 3D compresse sur K
  IF((K.LT.MAX(NKL,IKB).OR.K.GT.MIN(NKH,IKE)) .AND. .NOT.LKCP)THEN
    PTABO(JILOOP)=XSPVAL
  ELSE
! Ajout pour les PH definis avec _CV__K_  ou _Z_ etc... le 10/3/99
!   IF(LCV .AND. LCH)THEN
! idem (02/04/04) pour les _CV_ classiques (cas obs2mesonh)
! Ds ce cas on ne travaille pas necessairemnet sur les niveaux du modele
! mais sur un plan Z ou PR ou TK ou EV qui peut contenir des valeurs speciales
! XSPVAL. Il faut donc en tenir compte en limite de relief
! A PEAUFINER
! Le calcul des altitudes et du relief est fait mais je ne m'en sers pas
! Il peut etre aberrant . PENSER a les eliminer avec LPRINT et LPRINTXY
      IF((PTABI(IMIM1,IMJM1)==XSPVAL .AND. PTABI(IMIM1,IMJ)==XSPVAL).OR.&
	 (PTABI(IMIM1,IMJM1)==XSPVAL .AND. PTABI(IMI,IMJM1)==XSPVAL).OR.&
	 (PTABI(IMIM1,IMJM1)==XSPVAL .AND. PTABI(IMI,IMJ)==XSPVAL).OR. &
	 (PTABI(IMI,IMJM1)==XSPVAL .AND. PTABI(IMI,IMJ)==XSPVAL).OR. &
	 (PTABI(IMIM1,IMJ)==XSPVAL .AND. PTABI(IMI,IMJ)==XSPVAL).OR. &
	 (PTABI(IMIM1,IMJ)==XSPVAL .AND. PTABI(IMI,IMJM1)==XSPVAL))THEN
	 PTABO(JILOOP)=XSPVAL
      ELSE IF(PTABI(IMIM1,IMJM1)==XSPVAL .AND. PTABI(IMIM1,IMJ)/=XSPVAL.AND.&
	      PTABI(IMI,IMJM1)/=XSPVAL .AND. PTABI(IMI,IMJ)/=XSPVAL)THEN
        PTABO(JILOOP)=                                   &
          ZCIINF*ZCJSUP*PTABI(IMIM1,IMJ)+            &
          ZCISUP*ZCJINF*PTABI(IMI,IMJM1)+            &
          ZCISUP*ZCJSUP*PTABI(IMI,IMJ)    
      ELSE IF(PTABI(IMIM1,IMJM1)/=XSPVAL .AND. PTABI(IMIM1,IMJ)==XSPVAL.AND.&
	      PTABI(IMI,IMJM1)/=XSPVAL .AND. PTABI(IMI,IMJ)/=XSPVAL)THEN
        PTABO(JILOOP)=ZCIINF*ZCJINF*PTABI(IMIM1,IMJM1)+  &
          ZCISUP*ZCJINF*PTABI(IMI,IMJM1)+            &
          ZCISUP*ZCJSUP*PTABI(IMI,IMJ)    
      ELSE IF(PTABI(IMIM1,IMJM1)/=XSPVAL .AND. PTABI(IMIM1,IMJ)/=XSPVAL.AND.&
	      PTABI(IMI,IMJM1)==XSPVAL .AND. PTABI(IMI,IMJ)/=XSPVAL)THEN
        PTABO(JILOOP)=ZCIINF*ZCJINF*PTABI(IMIM1,IMJM1)+  &
          ZCIINF*ZCJSUP*PTABI(IMIM1,IMJ)+            &
          ZCISUP*ZCJSUP*PTABI(IMI,IMJ)    
      ELSE IF(PTABI(IMIM1,IMJM1)/=XSPVAL .AND. PTABI(IMIM1,IMJ)/=XSPVAL.AND.&
	      PTABI(IMI,IMJM1)/=XSPVAL .AND. PTABI(IMI,IMJ)==XSPVAL)THEN
        PTABO(JILOOP)=ZCIINF*ZCJINF*PTABI(IMIM1,IMJM1)+  &
          ZCIINF*ZCJSUP*PTABI(IMIM1,IMJ)+            &
          ZCISUP*ZCJINF*PTABI(IMI,IMJM1)
      ELSE
        PTABO(JILOOP)=ZCIINF*ZCJINF*PTABI(IMIM1,IMJM1)+  &
          ZCIINF*ZCJSUP*PTABI(IMIM1,IMJ)+            &
          ZCISUP*ZCJINF*PTABI(IMI,IMJM1)+            &
          ZCISUP*ZCJSUP*PTABI(IMI,IMJ)    
      ENDIF

!   ELSE
! Cas habituel
!   PTABO(JILOOP)=ZCIINF*ZCJINF*PTABI(IMIM1,IMJM1)+  &
!      ZCIINF*ZCJSUP*PTABI(IMIM1,IMJ)+            &
!      ZCISUP*ZCJINF*PTABI(IMI,IMJM1)+            &
!      ZCISUP*ZCJSUP*PTABI(IMI,IMJ)    
!   ENDIF
  END IF
!
!*   1.8     Computes the interpolated topography height for
!*           current section point
!
if(nverbia > 1)then
  print *,' ** coupe   AV XWZ '
endif

XWZ(JILOOP,NMGRID)=ZCIINF*ZCJINF*XXZS(IMIM1,IMJM1,NMGRID)+  &
          ZCIINF*ZCJSUP*XXZS(IMIM1,IMJ,NMGRID)+             &
          ZCISUP*ZCJINF*XXZS(IMI,IMJM1,NMGRID)+             &
          ZCISUP*ZCJSUP*XXZS(IMI,IMJ,NMGRID)    
!
ENDDO                     ! End of the general X- scanning loop

if(nverbia > 0)then
print *,' >>>> SORTIE COUPE  NMGRID= ',NMGRID,' size(XWZ)' ,size(XWZ,1)
endif
if(nverbia > 1)then
print *,' XWZ ',XWZ(:,NMGRID)
endif
!
RETURN
!------------------------------------------------------------------------
!
!*   2.     EXIT
!           ----
!
END SUBROUTINE COUPE_FORDIACHRO
