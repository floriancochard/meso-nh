!     #################################################
      SUBROUTINE COUPEUW_FORDIACHRO(PTABI,PTABO,K,KCOMP)
!     ##################################################
!
!!****  *COUPEUW_FORDIACHRO* - Vertical cross-section interpolation for U and W
!!                  wind components
!!
!!    PURPOSE
!!    -------
!         Interpolates 2D vertical cross-sections within the Meso-NH 3D
!       arrays. U and W model fields, appropriate gridpoint altitudes
!       as well as appropriate topography height are interpolated. 
!
!!**  METHOD
!!    ------
!!        The general case of a vertical cross-section along any oblique 
!!      direction with respect to the x-y model axes is considered. Simple
!!      linear interpolation is done.
!!      (First, wind components were co-located onto mass gridpoint)
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
!!      Original       19/09/95
!!      Updated   PM 
!-------------------------------------------------------------------------
!
!*     0.   DECLARATIONS
!           ------------
!
USE MODD_COORD
USE MODD_MEMCV
USE MODD_GRID1
USE MODN_NCAR
USE MODD_CVERT
USE MODN_PARA
USE MODD_PARAMETERS
USE MODD_NMGRID 
USE MODD_TYPE_AND_LH
USE MODD_RESOLVCAR  
USE MODD_MEMGRIUV
!
IMPLICIT NONE
!
!*     0.1  Dummy arguments and results
!
REAL,DIMENSION(:,:)      :: PTABI    ! Input data array to be interpolated
REAL,DIMENSION(:)        :: PTABO    ! Returned interpolated 1D array
INTEGER                  ::  K       ! Model level where interpolation is done
INTEGER                  ::  KCOMP   ! Code = 1 --> U wind component
                                     !      = 2 --> V       "
                                     !      = 3 --> W       "
!
!*     0.2  Local variables
!
REAL    :: ZCIINF,ZCISUP,ZCJINF,ZCJSUP,ZXX,ZYY
INTEGER :: IMIM1,IMI,IMJM1,IMJ,JILOOP, JI, JJ
INTEGER :: IIU, IJU, IKU, IKB, IKE
INTEGER :: IGRID
!
!------------------------------------------------------------------------------
!
!*     1.   PERFORMING VERTICAL INTERPOLATIONS
!           ----------------------------------
!
!*     1.0  Presetting array extends
!
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
IKE=IKU-JPVEXT
IKB=1+JPVEXT
!
!*     1.1  Scans along-section  oblique x axis 
!
! NOTA :
! L'utilisation explicite et volontaire de la valeur 1 comme dernier indice
! des tableaux presents dans la routine signifie que la representation se fait
! en definitive par rapport a la grille de masse en replacant  les com-
! -posantes du vent dans celle-ci
!
SELECT CASE(KCOMP)
  CASE(1)
    IGRID=2
    IF(NGRIU == 1)THEN
      IF(nverbia >0)then
      print *,' **coupeuw NGRIU,CGROUP ',NGRIU,CGROUP
      endif
      IGRID=1
    ENDIF
  CASE(2)
    IGRID=3
    IF(NGRIV == 1)THEN
      IF(nverbia >0)then
      print *,' **coupeuw NGRIV,CGROUP ',NGRIV,CGROUP
      endif
      IGRID=1
    ENDIF
  CASE(3)
    IGRID=1           ! W components put at mass gridpoints
END SELECT

if(nverbia > 0)then
print *,' **COUPEUW NMGRID DIRCUR ',NMGRID,'  ',CDIRCUR(1:LEN_TRIM(CDIRCUR))
endif
DO JILOOP=1,NLMAX   ! Start of general X- scanning loop
!
!*     1.2  Locates the current gridpoint along the cross-section 
!*          oblique x-axis within the Meso-NH grid
!
  ZXX=XDSX(JILOOP,1)        ! Collects the model X- and Y- axes projections
  ZYY=XDSY(JILOOP,1)        ! onto the oblique vertical section plane
                            ! for the current (new) section point
!
!*     1.3   The current section point is located along
!*           the x- axis
!
    DO JI=2,IIU
      IF(ZXX.LE.XXX(JI,IGRID).AND.ZXX.GE.XXX(JI-1,IGRID))GO TO 1
    ENDDO

1 CONTINUE

  IMIM1=MAX(1,JI-1)   
  IMI=JI             ! JI is the index of the first model bin to the left
!
!*    1.4    Then, it is located along
!*           the Y- axis
!
    DO JJ=2,IJU
      IF(ZYY.LE.XXY(JJ,IGRID).AND.ZYY.GE.XXY(JJ-1,IGRID))GO TO 2
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
    ZCIINF=(XXX(IMI,IGRID)-ZXX)/MAX(1.E-10,(XXX(IMI,IGRID)-XXX(IMIM1,IGRID)))
    ZCISUP=(ZXX-XXX(IMIM1,IGRID))/MAX(1.E-10,(XXX(IMI,IGRID)-XXX(IMIM1,IGRID)))
  END IF
  !
  IF(IMJ==IMJM1)THEN      ! Bottom wall special case
    ZCJINF=0.
    ZCJSUP=0.
  ELSE
    !PRINT *,'XY(IMJ IMJM1) ZXY',XY(IMJ),XY(IMJM1),ZYY
    ZCJINF=(XXY(IMJ,IGRID)-ZYY)/MAX(1.E-10,(XXY(IMJ,IGRID)-XXY(IMJM1,IGRID)))
    ZCJSUP=(ZYY-XXY(IMJM1,IGRID))/MAX(1.E-10,(XXY(IMJ,IGRID)-XXY(IMJM1,IGRID)))
  END IF
!
!*   1.6     Computes the interpolated altitude of the 
!*           current section point
!
  XWORKZ(JILOOP,K,1)=ZCIINF*ZCJINF*XZZ(IMIM1,IMJM1,K)+  &
         ZCIINF*ZCJSUP*XZZ(IMIM1,IMJ,K)+                     &
         ZCISUP*ZCJINF*XZZ(IMI,IMJM1,K)+                     &
         ZCISUP*ZCJSUP*XZZ(IMI,IMJ,K)    
!
!*   1.7     Computes the interpolated value of the field for
!*           current section point
! 

!  Modifs for diachro
  IF(K.LT.MAX(NKL,IKB).OR.K.GT.MIN(NKH,IKE))THEN
! IF(K.LT.IKB.OR.K.GT.IKE)THEN
    PTABO(JILOOP)=XSPVAL
  ELSE
    PTABO(JILOOP)=ZCIINF*ZCJINF*PTABI(IMIM1,IMJM1)+  &
          ZCIINF*ZCJSUP*PTABI(IMIM1,IMJ)+            &
          ZCISUP*ZCJINF*PTABI(IMI,IMJM1)+            &
          ZCISUP*ZCJSUP*PTABI(IMI,IMJ)    
  END IF
!
!*   1.8     Computes the interpolated topography height for
!*           current section point
!
XWZ(JILOOP,1)=ZCIINF*ZCJINF*XXZS(IMIM1,IMJM1,NMGRID)+  &
          ZCIINF*ZCJSUP*XXZS(IMIM1,IMJ,NMGRID)+             &
          ZCISUP*ZCJINF*XXZS(IMI,IMJM1,NMGRID)+             &
          ZCISUP*ZCJSUP*XXZS(IMI,IMJ,NMGRID)    
!
ENDDO                     ! End of the general X- scanning loop
!
RETURN
!------------------------------------------------------------------------
!
!*   2.     EXIT
!           ----
!
END SUBROUTINE COUPEUW_FORDIACHRO
