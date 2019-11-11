!     ######spl
      MODULE MODI_TRACEH_FORDIACHRO
!     #############################
!
INTERFACE
!
SUBROUTINE TRACEH_FORDIACHRO(KLREF,P3D,KLOOP)
REAL,DIMENSION(:,:,:)  :: P3D
INTEGER          :: KLREF,KLOOP
END SUBROUTINE TRACEH_FORDIACHRO
!
END INTERFACE
!
END MODULE MODI_TRACEH_FORDIACHRO
!     #############################################
      SUBROUTINE TRACEH_FORDIACHRO(KLREF,P3D,KLOOP)
!     #############################################
!
!!****  *TRACEH_FORDIACHRO* - Manager for the horizontal cross-section plots
!!
!!    PURPOSE
!!    -------
!       In the case of horizontal cross-sections, call the interpolation and
!     display routines: - for contour plots
!                       - for vector arrow plots
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!    VALNGRID   : retrieves the NGRID grid number when given the variable name
!!    COMCOORD   : computes true altitudes corresponding to the NGRID value
!!    INTERP     : vertically interpolates horizontal cross-sections
!!    IMAGE      : contour plot manager for horizontal cross-sections
!!    IMAGEV     : vector  plot manager for horizontal cross-sections
!!    READ_ALLVAR: reads any generic variable from the LFIFM file given its name
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_TITLE  : Declares heading variables for the plots (TRACE)
!!         NCONT  :  Current plot number
!!
!!      Module MODD_NMGRID : declares global variable  NMGRID (TRACE)
!!         NMGRID  : Current MESO-NH grid indicator
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist
!!                         (former NCAR common)
!!         NHI        : Extrema detection
!!                       (=0 --> H+L, <0 nothing)
!!         NDOT       : Line style
!!                       (=0|1|1023|65535 --> solid lines;
!!                       <0 --> solid lines for positive values and
!!                       dotted lines(ABS(NDOT))for negative values;
!!                       >0 --> dotted lines(ABS(NDOT)) )
!!         CTYPHOR    : Horizontal cross-section type
!!                       (='K' --> model level section;
!!                        ='Z' --> constant-altitude section;
!!                        ='P' --> isobar section (planned)
!!                        ='T' --> isentrope section (planned)
!!
!!      Module MODD_OUT  : defines various logical units and dimensions
!!         NIMAXT     : x-size of the displayed section of the MESO-NH arrays
!!         NJMAXT     : y-size of the displayed section of the MESO-NH arrays
!!         NKMAXT     : z-size of the displayed section of the MESO-NH arrays
!!
!!      Module MODN_PARA
!!         Module MODD_DIM1 : contains dimensions of data arrays
!!             NKMAX :  z array dimensions  
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!         JPVEXT   : Vertical external points number
!!
!!      Module MODD_SUPER   : defines plot overlay control variables 
!!         LSUPER   : =.TRUE. --> plot overlay is active
!!                    =.FALSE. --> plot overlay is not active
!!         NSUPER   : Rank of the current plot in the overlay
!!                    sequence. The initial plot is rank 1.
!!
!!
!!      Module MODD_ALLVAR  : contains generic variables arrays and structures
!!         XWORK3D  : 3D generic scalar field array
!!         XWORKX3D : 3D generic vector field x-component array
!!         XWORKY3D : 3D generic vector field y-component array
!!         XWORKZ3D : 3D generic vector field z-component array
!!         XWORK2D  : 2D generic scalar field
!!         XWORKX3D : 2D generic vector field x-component array
!!         XWORKY3D : 2D generic vector field y-component array
!!>>>>>DRAGOON
!!>>>>>DRAGOON NOTICE: I don't see why a 2D generic vector should not have
!!>>>>>DRAGOON         a w-component as well. Exemple: a 2D map of the u-w 
!!                     vectors...
!!>>>>>DRAGOON
!!         XT1      : structure defining the name, grid number and unit name
!!                    for a 3D generic scalar field (TRACE derived type X_Y_Z_)
!!         XT2      : structure defining the name, grid number and unit name
!!                    for a 2D generic scalar field (TRACE derived type X_Y_)
!!         XT3      : structure defining the name, grid numbers and unit name
!!                    for a 3D generic 3D vector field 
!!                    (TRACE derived type VX_VY_VZ_)
!!         XT4      : structure defining the name, grid numbers and unit name
!!                    for a 2D generic 2D vector  field 
!!                    (TRACE derived type VX_VY_)
!!         XT5      : structure defining the name, grid number and unit name
!!                    for a 1D generic scalar field (TRACE derived type Z_)
!!
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
!!      Updated   PM   06/12/94
!!      Updated   JD   09/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TITLE
USE MODD_MASK3D
USE MODD_TIT
USE MODD_DEFCV
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO
USE MODD_NMGRID
USE MODN_NCAR
USE MODD_OUT
USE MODD_DIM1
USE MODN_PARA
USE MODD_PARAMETERS
USE MODD_TYPE_AND_LH
USE MODD_SUPER
USE MODD_ALLVAR
USE MODI_INTERP_FORDIACHRO
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_COORD
USE MODI_RESOLV_TIT
USE MODI_RESOLV_TITY
USE MODI_COMPUTEDIR

IMPLICIT NONE
!
!*      0.1   Interfaces declaration
!
INTERFACE
      SUBROUTINE PRECOU_FORDIACHRO(PWORK3D,PTEMCV)
      REAL,DIMENSION(:,:,:)  :: PWORK3D
      REAL,DIMENSION(:,:)    :: PTEMCV
      END SUBROUTINE PRECOU_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE IMAGE_FORDIACHRO(PTAB,KLREF,PTABINT,KNHI,KNDOT,HTEXTE)
      CHARACTER(LEN=*)   :: HTEXTE
      REAL                :: PTABINT
      REAL,DIMENSION(:,:) :: PTAB
      INTEGER :: KNHI, KNDOT, KLREF
      END SUBROUTINE IMAGE_FORDIACHRO
END INTERFACE
INTERFACE 
      SUBROUTINE IMAGEV_FORDIACHRO(PU,PV,KLREF,HTEXTE)
      REAL,DIMENSION(:,:) :: PU,PV
      CHARACTER(LEN=*) :: HTEXTE
      INTEGER :: KLREF
      END SUBROUTINE IMAGEV_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE TRAXY(PTEMX,PTEMY,KLOOP,HTITX,HTITY,PTIMED,PTIMEF)
      INTEGER    :: KLOOP
      REAL,DIMENSION(:)  :: PTEMX, PTEMY
      REAL               :: PTIMED, PTIMEF
      CHARACTER(LEN=*) :: HTITX, HTITY
      END SUBROUTINE TRAXY
END INTERFACE
INTERFACE 
      SUBROUTINE TRAHTRAXY(KLOOP,PTEMCV,HTEXTE)
      INTEGER :: KLOOP
      REAL,DIMENSION(:,:) :: PTEMCV
      CHARACTER(LEN=40)   :: HTEXTE
      END SUBROUTINE TRAHTRAXY
END INTERFACE
!
COMMON/TEMH/XZZX,XZZY,NIIMAX,NIJMAX
#include "big.h"
REAL,DIMENSION(N2DVERTX) :: XZZX
REAL,DIMENSION(N2DVERTX) :: XZZY
INTEGER :: NIIMAX, NIJMAX
!
!
!*       0.15 Dummy arguments
!          
INTEGER          :: KLREF, KLOOP, JU
REAL,DIMENSION(:,:,:)   :: P3D
!
!*       0.2  local variables
!          

REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE   :: ZTEM, ZTEM2
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE     :: ZX
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE     :: ZSTAB, ZSTAB1, ZSTAB2, ZTEMCV,ZSTABM
REAL,DIMENSION(:),ALLOCATABLE,SAVE       :: ZTEMX, ZTEMY
REAL,DIMENSION(:),ALLOCATABLE,SAVE       :: ZZY

REAL             :: ZTIMED, ZTIMEF

INTEGER          :: ITER, JTER, IUB1, IUB2, ISKIP
INTEGER          :: ISTA, IER, INB, IWK
INTEGER          :: IWIU, IWJU
INTEGER          :: ILENT, ILENU, ILENCTIMECS, ILENE
INTEGER          :: IBEGTXT, IENDTXT

CHARACTER(LEN=20) :: YNOM
CHARACTER(LEN=40) :: YTEXTE
CHARACTER(LEN=15) :: YEND
CHARACTER(LEN=8)  :: YCAR8
CHARACTER(LEN=16) :: YTITX, YTITY

!
!-------------------------------------------------------------------------------
!
!*      1.    PRELIMINARY CALCULATION
!             -----------------------
!
if(nverbia > 0)then
 print *,' **entree traceh  LPR,LTK,LEV,LSV3,CTYPHOR ', LPR,LTK,LEV,LSV3,CTYPHOR
endif
NIMAXT=NISUP-NIINF+1
NJMAXT=NJSUP-NJINF+1
NKMAXT=NKMAX+2*JPVEXT
!
!*       1.1    Array allocations
!
IF(ALLOCATED(ZSTAB))THEN
  DEALLOCATE(ZSTAB)
END IF
  ALLOCATE(ZSTAB(NIMAXT,NJMAXT))
!
!*       1.2    NCAR setting
!
!
!*       1.3    Interactive option selection and plot overlay management  
!
!
IWIU=SIZE(P3D,1)
IWJU=SIZE(P3D,2)
if(nverbia >0)then
  print *,' ** Entree traceh KLREF ',KLREF
endif
YNOM=ADJUSTL(CGROUP)
IF(YNOM.EQ.'QUIT')THEN
  CALL GQOPS(ISTA)
  CALL GQACWK(1,IER,INB,IWK)
  IF(ISTA >1 .AND. INB >1)THEN
    CALL GDAWK(2)
    CALL GCLWK(2)
  ENDIF
! CALL FRAME
  CALL NGPICT(1,1)
  CALL CLSGKS
  STOP
END IF

IBEGTXT=1
IENDTXT=LEN(YTEXTE)

IF(NSUPERDIA > 1)THEN
  IF(LMINUS .OR. LPLUS)THEN
		      IF(NBPM > 1)THEN
			DO JU=1,NBPM
			  IF(NUMPM(JU) == 3)THEN
		            LSUPER=.TRUE.
			    EXIT
			  ELSE
		            LSUPER=.FALSE.
			  ENDIF
			ENDDO
		      ELSE
		        LSUPER=.FALSE.
		      ENDIF
  ELSE
    LSUPER=.TRUE.
  ENDIF
ELSE
  LSUPER=.FALSE.
ENDIF
IF(KLOOP == 1)NSUPER=0
XLWIDTH=XLWDEF
!
!
! Selects "model levels" mode
!
!
! Selects altitude mode
!

! If no keyword has been detected so far, TRACE tries to 
! interpret the last entry as a new model level number.
!
!
IF(.NOT.LCN .AND. .NOT.LCNCUM)THEN
IF (CTYPHOR.EQ.'K')THEN
  IF(LMSKTOP)THEN
    KLREF=NKH
  ELSE
  IF(KLREF.GT.NKH.OR.KLREF.LT.NKL)THEN
!  IF(KLREF.GT.NKMAX+2*JPVEXT.OR.KLREF.LT.1)THEN
    print *,' This model level is unknown!'
  END IF
  END IF
END IF
END IF
!
!*       2.     PROCESSING OF THE BASIC SET OF VARIABLES
!               ---------------------------------------------------

! WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A8,2X,A1,''='',I5)')CGROUP,CTYPHOR,KLREF

  YTEXTE(1:LEN(YTEXTE))=' '
  IBEGTXT=1
  ILENT=LEN_TRIM(CTITGAL)
  ILENU=LEN_TRIM(CUNITGAL)

IF(LCN .OR. LCNCUM)THEN
  YTEXTE(1:ILENT)=CTITGAL(1:ILENT)

  ZSTAB(:,:)=P3D(:,:,1)
  CALL COMPCOORD_FORDIACHRO(NMGRID)
  IF(LCN)THEN
    YTEXTE(ILENT+1:ILENT+1)=' '
    YTEXTE(ILENT+2:ILENT+9)=ADJUSTL(CTIMEC(8:15))
  ELSE
    YCAR8(1:LEN(YCAR8))=' '
    YTEXTE(ILENT+1:ILENT+1)=' '
    YTEXTE(ILENT+2:ILENT+9)=ADJUSTL(CTIMECS(8:15))
    ILENT=LEN_TRIM(YTEXTE)
    YTEXTE(ILENT+1:ILENT+1)='-'
    ILENCTIMECS=LEN_TRIM(CTIMECS)
    YCAR8=CTIMECS(ILENCTIMECS-7:ILENCTIMECS)
    YTEXTE(ILENT+2:ILENT+9)=ADJUSTL(YCAR8)
  ENDIF
  CALL IMAGE_FORDIACHRO(ZSTAB,KLREF,XDIAINT,NHI,NDOT,YTEXTE(1: &
  LEN_TRIM(YTEXTE)))

ELSE

  !YTEXTE(ILENT+1:ILENT+1)=' '
  !YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
  !IBEGTXT=ILENT+2+ILENU
  !YTEXTE(IBEGTXT:IBEGTXT+2)=' '
  !IBEGTXT=IBEGTXT+3
  
  YEND(1:LEN(YEND))=' '
  
  IF(LEV .AND. CTYPHOR(1:1)=='E')THEN
    IF(LCHREEL)THEN
      WRITE(YEND,'(A2,''='',F7.1)')'PV',XLOOPZ
      ILENE=LEN_TRIM(YEND)
      !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A2,''='',F7.1)')'PV',XLOOPZ
    ELSE
      WRITE(YEND,'(A2,''='',I5)')'PV',KLREF
      ILENE=LEN_TRIM(YEND)
      !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A2,''='',I5)')'PV',KLREF
    ENDIF
  ELSE IF(LMSKTOP)THEN
    WRITE(YEND,'(A9)')' MSKTOP=T'
    ILENE=LEN_TRIM(YEND)
    !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A9)')' MSKTOP=T'
  ELSE IF(LSV3)THEN
    IF(LXYZ00)THEN
      IF(LCHREEL .AND. CTYPHOR /= 'K')THEN
        WRITE(YEND,'(A5,''='',F7.1)')CGROUPSV3(1:5),XLOOPZ
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A4,''='',F7.1)')CGROUPSV3(1:4),XLOOPZ
      ELSE
        WRITE(YEND,'(A5,''='',I5)')CGROUPSV3(1:5),KLREF
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A4,''='',I5)')CGROUPSV3(1:4),KLREF
!     WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A4,''='',I5)')' Z00',KLREF
      ENDIF
    ELSE
      IF(LCHREEL  .AND. CTYPHOR /= 'K')THEN
        WRITE(YEND,'(A3,''='',F7.1)')'SV3',XLOOPZ
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A3,''='',F7.1)')'SV3',XLOOPZ
      ELSE
        WRITE(YEND,'(A3,''='',I5)')'SV3',KLREF
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A3,''='',I5)')'SV3',KLREF
      ENDIF
    ENDIF
  ELSE
    IF(LXYZ)THEN
      IF(LCHREEL  .AND. CTYPHOR /= 'K')THEN
        WRITE(YEND,'(A1,''='',F7.1,A6)')CTYPHOR,XLOOPZ,' MSK=T'
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A1,''='',F7.1,A6)')CTYPHOR,XLOOPZ,' MSK=T'
      ELSE
        WRITE(YEND,'(A1,''='',I5,A6)')CTYPHOR,KLREF,' MSK=T'
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A1,''='',I5,A6)')CTYPHOR,KLREF,' MSK=T'
      ENDIF
    ELSE
      IF(LCHREEL  .AND. CTYPHOR /= 'K')THEN
        WRITE(YEND,'(A1,''='',F7.1)')CTYPHOR,XLOOPZ
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A1,''='',F7.1)')CTYPHOR,XLOOPZ
      ELSE
        WRITE(YEND,'(A1,''='',I5)')CTYPHOR,KLREF
        ILENE=LEN_TRIM(YEND)
        !WRITE(YTEXTE(IBEGTXT:IENDTXT),'(A1,''='',I5)')CTYPHOR,KLREF
      ENDIF
    ENDIF
  ENDIF
  ! YTEXTE est rempli a partir de la fin
  IBEGTXT=IENDTXT-ILENE+1
  YTEXTE(IBEGTXT:IENDTXT)=TRIM(YEND)
  ! 3 car blancs
  IENDTXT=IBEGTXT-1
  IBEGTXT=IENDTXT-2
  YTEXTE(IBEGTXT:IENDTXT)=' '
  ! unite
  IENDTXT=IBEGTXT-1
  IBEGTXT=IENDTXT-ILENU+1
  YTEXTE(IBEGTXT:IENDTXT)=CUNITGAL(1:ILENU)
  ! 1 car blanc
  IENDTXT=IBEGTXT-1
  IBEGTXT=IENDTXT
  YTEXTE(IBEGTXT:IENDTXT)=' '
  ! titre (tronque eventuellement)
  IENDTXT=IBEGTXT-1
  IBEGTXT=MAX(1,IENDTXT-ILENT+1)
  YTEXTE(IBEGTXT:IENDTXT)=CTITGAL(1:ILENT)
  YTEXTE=ADJUSTL(YTEXTE)
if(nverbia > 0)then
  print*,' ** TRACEH: TIT=',CTITGAL(1:ILENT),' UNIT=',CUNITGAL(1:ILENU),&
         ' TEXTE= ',TRIM(YTEXTE)
endif

  CALL COMPCOORD_FORDIACHRO(NMGRID)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
     LDIRWM .OR. LDIRWT .OR. &
     LSUMVM .OR. LSUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
    ALLOCATE(ZTEM(IWIU,IWJU,SIZE(P3D,3)))
    ALLOCATE(ZSTAB1(IWIU,IWJU))
    ALLOCATE(ZSTAB2(IWIU,IWJU))
    IF(NGRIDIA(1) == 1)THEN
      ZTEM(:,:,:)=P3D(:,:,:)
      print *,' *****TRACEH..  PAS D''INTERPOLATION de U sur la grille de masse, GROUPE: ',CGROUP,' IGRID: ',NGRIDIA(1)
    ELSE
      ZTEM(1:IWIU-1,:,:)=0.5*(P3D(1:IWIU-1,:,:)+P3D(2:IWIU,:,:))
      ZTEM(IWIU,:,:)=2.*ZTEM(IWIU-1,:,:)-ZTEM(IWIU-2,:,:)
    ENDIF
!!!!!!!!!!PROVISOIRE POUR VERIF
!   ZTEM(5,:,:)=10.
    CALL INTERP_FORDIACHRO(KLREF,NKL,NKH,ZTEM,ZSTAB1)
if(nverbia >0)then
  print *,' ** Traceh AP INTERP1 KLREF ',KLREF
endif

! Avril 2000 LCV+LCH+LUMVM ou LUTVT -> PH vecteurs
! Traitement U
    IF(LCV)THEN
! Je remets le plan horizontal demande (ZSTAB1) arbitrairement au niveau 2 ou
! NKL de ZTEM et je fais toutes les operations concernant une coupe verticale
! (sauf pour le 2D horizontal)
      IF(SIZE(ZTEM,3) == 1)THEN
        ZTEM(:,:,1)=ZSTAB1(:,:)
      ELSE
        ZTEM(:,:,MAX(2,NKL))=ZSTAB1(:,:)
      ENDIF
      CALL VERIFLEN_FORDIACHRO
      CALL MEMCV
      IF(ALLOCATED(ZTEMCV))THEN
        DEALLOCATE(ZTEMCV)
      ENDIF
      IF(NVERBIA >0)THEN
        print *,' ** TRACEH av PRECOU NLMAX IKU ',NLMAX,SIZE(ZTEM,3)
      ENDIF
      ALLOCATE(ZTEMCV(NLMAX,1:SIZE(ZTEM,3)))
      IF(ALLOCATED(XTEMCVU))THEN
	DEALLOCATE(XTEMCVU)
      ENDIF
      ALLOCATE(XTEMCVU(NLMAX,1))
      CALL PRECOU_FORDIACHRO(ZTEM,ZTEMCV)
      IF(SIZE(ZTEMCV,2) == 1)THEN
        XTEMCVU(:,1)=ZTEMCV(:,1)
      ELSE
        XTEMCVU(:,1)=ZTEMCV(:,2)
      ENDIF
      DEALLOCATE(ZTEMCV)
    ENDIF
! Avril 2000

    DEALLOCATE(ZTEM)
    IF(NVERBIA > 0)THEN
    print *,' DEALLOCATE(ZTEM) '
    ENDIF

! Traitement V
    ALLOCATE(ZTEM2(IWIU,IWJU,SIZE(P3D,3)))
    ZTEM2(:,:,:)=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
               NJINF-NJL+1:NJSUP-NJL+1, &
             :,NLOOPT,1,1)
    IF(NGRIDIA(1) == 1)THEN
      print *,' *****TRACEH..  PAS D''INTERPOLATION de V sur la grille de masse, GROUPE: ',CGROUP,' IGRID: ',NGRIDIA(1)
    ELSE
      ZTEM2(:,1:IWJU-1,:)=0.5*(ZTEM2(:,1:IWJU-1,:)+ZTEM2(:,2:IWJU,:))
      ZTEM2(:,IWJU,:)=2.*ZTEM2(:,IWJU-1,:)-ZTEM2(:,IWJU-2,:)
    ENDIF
!!!!!!!!!!PROVISOIRE POUR VERIF
!   ZTEM(:,10,:)=10.
    CALL INTERP_FORDIACHRO(KLREF,NKL,NKH,ZTEM2,ZSTAB2)
if(nverbia >0)then
  print *,' ** Traceh AP INTERP2 KLREF ',KLREF
endif
! Avril 2000 LCV+LCH+LUMVM ou LUTVT -> PH vecteurs
    IF(LCV)THEN
! Je remets le plan horizontal demande (ZSTAB2) arbitrairement au niveau 2 ou
! NKL de ZTEM2 et je fais toutes les operations concernant une coupe verticale
! (sauf pour le 2D horizontal)
      IF(SIZE(ZTEM2,3) == 1)THEN
        ZTEM2(:,:,1)=ZSTAB2(:,:)
      ELSE
        ZTEM2(:,:,MAX(2,NKL))= ZSTAB2(:,:)
      ENDIF
      CALL VERIFLEN_FORDIACHRO
      CALL MEMCV
      IF(ALLOCATED(ZTEMCV))THEN
        DEALLOCATE(ZTEMCV)
      ENDIF
      IF(NVERBIA >0)THEN
        print *,' ** TRACEH av PRECOU NLMAX IKU ',NLMAX,SIZE(ZTEM2,3)
      ENDIF
      ALLOCATE(ZTEMCV(NLMAX,1:SIZE(ZTEM2,3)))
      IF(ALLOCATED(XTEMCVV))THEN
	DEALLOCATE(XTEMCVV)
      ENDIF
      ALLOCATE(XTEMCVV(NLMAX,1))
      CALL PRECOU_FORDIACHRO(ZTEM2,ZTEMCV)
! Nov 2001
!     XTEMCVV(:,1)=ZTEMCV(:,2)
      IF(SIZE(ZTEMCV,2) == 1)THEN
        XTEMCVV(:,1)=ZTEMCV(:,1)
      ELSE
        XTEMCVV(:,1)=ZTEMCV(:,2)
      ENDIF
! Nov 2001
      DEALLOCATE(ZTEMCV)
    ENDIF
! Avril 2000
    DEALLOCATE(ZTEM2)
    IF(NVERBIA > 0)THEN
    print *,' DEALLOCATE(ZTEM2) '
    ENDIF

    IF(LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT .OR. LDIRWM .OR. LDIRWT .OR. &
       LMUMVM .OR.LMUTVT)THEN
! Avril 2000 LCV+LCH+LUMVM ou LUTVT -> PH vecteurs
      IF(LCV)THEN
!! Nov 2001
         IF(LMUMVM .OR.LMUTVT)THEN
           IF(ALLOCATED(ZTEMCV))THEN
             DEALLOCATE(ZTEMCV)
           ENDIF
!          print *,' XTEMCVU ',XTEMCVU
!          print *,' XTEMCVV ',XTEMCVV
           ALLOCATE(ZTEMCV(SIZE(XTEMCVU,1),SIZE(XTEMCVU,2)))
           WHERE(XTEMCVV(:,:) == XSPVAL)XTEMCVU=XSPVAL
           WHERE(XTEMCVU(:,:) == XSPVAL)XTEMCVV=XSPVAL
           WHERE(XTEMCVU(:,:) /= XSPVAL)ZTEMCV=XTEMCVU*XTEMCVU
           XTEMCVU=ZTEMCV
           WHERE(XTEMCVV(:,:) /= XSPVAL)ZTEMCV=XTEMCVV*XTEMCVV
           XTEMCVV=ZTEMCV
           WHERE(XTEMCVU(:,:) /= XSPVAL)XTEMCVU=SQRT(XTEMCVU+XTEMCVV)
           CALL TRAHTRAXY(KLOOP,XTEMCVU,YTEXTE)
           DEALLOCATE(ZTEMCV)

         ELSEIF(LDIRWM .OR. LDIRWT)THEN
           IUB1=SIZE(XTEMCVU,1)
           IUB2=SIZE(XTEMCVU,2)
           ISKIP=1
           ITER=IUB1; JTER=IUB2
           IF(ALLOCATED(ZX))THEN
             DEALLOCATE(ZX)
           ENDIF
           IF(ALLOCATED(ZZY))THEN
             DEALLOCATE(ZZY)
           ENDIF
           ALLOCATE(ZX(ITER,1),ZZY(JTER))
          print *,' **traceh av ZX,ZZY '
           ZX(:,1)=XZZX(1:IUB1:ISKIP)
           ZZY=XZZY(1:IUB2:ISKIP)
          print *,' **traceh aP ZX,ZZY ',ZX(1:IUB1,1)
          print *,' **traceh aP ZX,ZZY ',ZZY(1:IUB2)
! Calcul de la direction du vent par DIR.... Retour ds XTEMCVV
           CALL COMPUTEDIR(ITER,JTER,IUB1,IUB2,ISKIP,XTEMCVU,XTEMCVV)
          print *,' **traceh ap computedir , av trahtraxy'
           CALL TRAHTRAXY(KLOOP,XTEMCVV,YTEXTE)

         ENDIF
!! Nov 2001
      ELSE
! Avril 2000
        IF(LMUMVM .OR.LMUTVT)THEN
          ZSTAB(:,:)=SQRT(ZSTAB1(:,:)**2+ZSTAB2(:,:)**2)
          WHERE(ZSTAB1(:,:) == XSPVAL)ZSTAB=XSPVAL
          WHERE(ZSTAB2(:,:) == XSPVAL)ZSTAB=XSPVAL
          CALL IMAGE_FORDIACHRO(ZSTAB,KLREF,XDIAINT,NHI,NDOT,YTEXTE(1:LEN_TRIM(YTEXTE)))
        ELSE IF((LDIRWM.OR.LDIRWT).AND. .NOT. LDIRWIND) THEN
          !! direction par DD....
          print*,'traceh dd ',minval(ZSTAB1),maxval(ZSTAB1),minval(ZSTAB2), &
                                             maxval(ZSTAB2)
          print*,'traceh dd ',minloc(ZSTAB1),maxloc(ZSTAB1),minloc(ZSTAB2), &
                                             maxloc(ZSTAB2)
          IUB1=SIZE(ZSTAB1,1)
          IUB2=SIZE(ZSTAB1,2)
          ISKIP=1
          ITER=IUB1; JTER=IUB2
          XZZX(1:IUB1)=XXX(NIINF:NISUP,NMGRID)
          XZZY(1:IUB2)=XXY(NJINF:NJSUP,NMGRID)
          CALL COMPUTEDIR(ITER,JTER,IUB1,IUB2,ISKIP,ZSTAB1,ZSTAB2)
          print*,'traceh dd ',minval(ZSTAB2),maxval(ZSTAB2)
          print*,'traceh dd ',minloc(ZSTAB2),maxloc(ZSTAB2)
          CALL IMAGE_FORDIACHRO(ZSTAB2,KLREF,XDIAINT,NHI,NDOT,YTEXTE(1:LEN_TRIM(YTEXTE)))
        ELSE
        CALL IMAGEV_FORDIACHRO(ZSTAB1,ZSTAB2,KLREF,YTEXTE)
        ENDIF
      ENDIF
! Avril 2000

    ELSE

      ZSTAB(:,:)=SQRT(ZSTAB1(:,:)**2+ZSTAB2(:,:)**2)
      WHERE(ZSTAB1(:,:) == XSPVAL)ZSTAB=XSPVAL
      WHERE(ZSTAB2(:,:) == XSPVAL)ZSTAB=XSPVAL
      CALL IMAGE_FORDIACHRO(ZSTAB,KLREF,XDIAINT,NHI,NDOT,YTEXTE(1:LEN_TRIM(YTEXTE)))

    ENDIF

    IF(ALLOCATED(ZTEM))DEALLOCATE(ZTEM)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL INTERP_FORDIACHRO(KLREF,NKL,NKH,P3D,ZSTAB)
if(nverbia >0)then
  print *,' ** Traceh AP INTERP3 KLREF ',KLREF
endif
!print *,' ZSTAB'
!print *,ZSTAB
!REGLER LE PB DE L'INTERVALLE

    IF(NJMAXT == 1 .AND. NIMAXT /= 1)THEN         !;;;;;;;;;;;;
      IF(ALLOCATED(ZTEMX))THEN
	DEALLOCATE(ZTEMX)
      ENDIF
      IF(ALLOCATED(ZTEMY))THEN
	DEALLOCATE(ZTEMY)
      ENDIF
      ALLOCATE(ZTEMX(SIZE(ZSTAB,1)))
      ALLOCATE(ZTEMY(SIZE(ZSTAB,1)))
      ZTEMX(:)=XXX(NIINF:NISUP,NMGRID)
! Ajout Nov 99
      ZTEMX(:)=ZTEMX(:)-XXX(NIINF,NMGRID)
! Ajout Nov 99
      ZTEMY(:)=ZSTAB(:,1)
      WHERE(ZTEMY == XSPVAL)
!     WHERE(ZTEMY == 999.)
	ZTEMY=1.E36
      END WHERE
      YTITX(1:LEN(YTITX))=' '
      YTITY(1:LEN(YTITX))=' '
      YTITX='X(M)'
      YTITY=CUNITGAL(1:LEN(CUNITGAL))
      ZTIMED=XTRAJT(NLOOPT,1)
      ZTIMEF=ZTIMED
      CALL TRAXY(ZTEMX,ZTEMY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)
      IF(KLOOP == 1)THEN
	IF(LDATFILE)CALL DATFILE_FORDIACHRO
	CALL RESOLV_TIMES(NLOOPT)
	CALL PLCHHQ(.99,.007,YTEXTE(1:LEN_TRIM(YTEXTE)),.011,0.,+1.)
      ENDIF

    ELSE IF(NIMAXT == 1 .AND. NJMAXT /= 1)THEN    !;;;;;;;;;;;;

      IF(ALLOCATED(ZTEMX))THEN
	DEALLOCATE(ZTEMX)
      ENDIF
      IF(ALLOCATED(ZTEMY))THEN
	DEALLOCATE(ZTEMY)
      ENDIF
      ALLOCATE(ZTEMX(SIZE(ZSTAB,2)))
      ALLOCATE(ZTEMY(SIZE(ZSTAB,2)))
      ZTEMX(:)=XXY(NJINF:NJSUP,NMGRID)
! Ajout Nov 99
      ZTEMX(:)=ZTEMX(:)-XXY(NJINF,NMGRID)
! Ajout Nov 99
      ZTEMY(:)=ZSTAB(1,:)
      WHERE(ZTEMY == XSPVAL)
!     WHERE(ZTEMY == 999.)
	ZTEMY=1.E36
      END WHERE
      YTITX(1:LEN(YTITX))=' '
      YTITY(1:LEN(YTITX))=' '
      YTITX='Y(M)'
      YTITY=CUNITGAL(1:LEN(CUNITGAL))
      ZTIMED=XTRAJT(NLOOPT,1)
      ZTIMEF=ZTIMED
      CALL TRAXY(ZTEMX,ZTEMY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)
      IF(KLOOP == 1)THEN
	IF(LDATFILE)CALL DATFILE_FORDIACHRO
	CALL RESOLV_TIMES(NLOOPT)
	CALL PLCHHQ(.99,.007,YTEXTE(1:LEN_TRIM(YTEXTE)),.011,0.,+1.)
      ENDIF

    ELSE                                          !;;;;;;;;;;;;

! Ajout PH = intersection CV et CH 10/3/99 (Defini avec _cv__k_ (ou _z_ etc))
      IF(LCV)THEN       !......................................

! Je remets le plan horizontal demande (ZSTAB) arbitrairement au niveau NKL
! de P3D(ZWORK3D) et je fais toutes les operations concernant une coupe verticale
! Je recupere le profil dans ZTEMCV(1:NLMAX,NKL)
! J'ai les X ds XDS(1:NLMAX) Penser a mettre les latlon pts extremes
       IF(NVERBIA >0)THEN
	 print *,' ** TRACEH SZ(1,2) de P3D et ZSTAB NKL ',&
	 SIZE(P3D,1),SIZE(P3D,2),SIZE(ZSTAB,1),SIZE(ZSTAB,2),NKL
       ENDIF
       ALLOCATE(ZSTABM(SIZE(ZSTAB,1),SIZE(ZSTAB,2)))
! prise en compte du 2D hor. -> PH Oct 2000)
       IF(SIZE(P3D,3) == 1)THEN
         ZSTABM(:,:)=P3D(:,:,1)
         P3D(:,:,1)=ZSTAB(:,:)
       ELSE
         ZSTABM(:,:)=P3D(:,:,MAX(2,NKL))
         P3D(:,:,MAX(2,NKL))=ZSTAB(:,:)
       ENDIF
       CALL VERIFLEN_FORDIACHRO
       CALL MEMCV
       IF(ALLOCATED(ZTEMCV))THEN
	 DEALLOCATE(ZTEMCV)
       ENDIF
       IF(NVERBIA >0)THEN
	 print *,' ** TRACEH av PRECOU NLMAX IKU ',NLMAX,SIZE(P3D,3)
       ENDIF
       ALLOCATE(ZTEMCV(NLMAX,1:SIZE(P3D,3)))
       CALL PRECOU_FORDIACHRO(P3D,ZTEMCV)
! prise en compte du 2D hor. -> PH Oct 2000)
       IF(SIZE(P3D,3) == 1)THEN
         P3D(:,:,1)=ZSTABM(:,:)
       ELSE
         P3D(:,:,MAX(2,NKL))=ZSTABM(:,:)
       ENDIF
       DEALLOCATE(ZSTABM)
       IF(NVERBIA >0)THEN
	 IF(SIZE(P3D,3) == 1)THEN
	   print *,' ** TRACEH ap PRECOU ZTEMCV(:,NKL)', ZTEMCV(:,1)
	 ELSE
	   print *,' ** TRACEH ap PRECOU ZTEMCV(:,NKL)', ZTEMCV(:,MAX(2,NKL))
	 ENDIF
       ENDIF
!!!!!!!!!!!!! Supprime le 30/11/01
       CALL TRAHTRAXY(KLOOP,ZTEMCV,YTEXTE)
!!!!!!!!!!!!! Supprime le 30/11/01
	
      ELSE              !......................................

        CALL IMAGE_FORDIACHRO(ZSTAB,KLREF,XDIAINT,NHI,NDOT,YTEXTE(1:LEN_TRIM(YTEXTE)))
      ENDIF             !......................................

    ENDIF                                         !;;;;;;;;;;;;

  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDIF
!

! CALL FRAME
!
DEALLOCATE(ZSTAB)
IF(ALLOCATED(ZSTAB1))THEN
  DEALLOCATE(ZSTAB1)
ENDIF
IF(ALLOCATED(ZSTAB2))THEN
  DEALLOCATE(ZSTAB2)
ENDIF
if(nverbia > 0)then
 print *,' **sortie traceh  LPR,LTK,LEV,LSV3 ', LPR,LTK,LEV,LSV3
endif
  RETURN
!------------------------------------------------------------------------------
!
!*     5.       EXIT
!               ----
!
1000 FORMAT(5X,I4,3X,A12)
!
!*     5.1         Heading formats
!
1001 FORMAT('Horiz. profile XDEB=',F6.0,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1002 FORMAT('Horiz. profile XDEB=',F6.0,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1003 FORMAT('Horiz. profile XDEB=',E7.2,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1004 FORMAT('Horiz. profile XDEB=',E6.2,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1010 FORMAT('Horiz. profile IDEB=',I4,' JDEB=',I4,' ANG.=',I3,' NBPTS=',I4)
1011 FORMAT('Horiz. profile XDEB=',F6.0,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1013 FORMAT('Horiz. profile XDEB=',F6.0,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1014 FORMAT('Horiz. profile XDEB=',E7.2,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I4)
1015 FORMAT('Horiz. profile XDEB=',E6.2,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I4)
1018 FORMAT('Horiz. profile IND I,J (BEGIN)-(END)=(',I4,',',I4,')-(',I4,',',I4,')')
1019 FORMAT('Horiz. profile LAT,LON (BEGIN)-(END)=(',F4.1,',',F5.1,')-(',F4.1,',',F5.1,')')
1020 FORMAT('Horiz. profile CONF. COORD.(BEGIN)-(END)=(',F8.0,',',F8.0,')-(',F8.0,',',F8.0,')')
!
END SUBROUTINE  TRACEH_FORDIACHRO
