!     ############################
      SUBROUTINE CONVXY2IJ(HCARIN)
!     ############################
!
!!****  *CONVXY2IJ* - Convertit des coordonnees conformes et coordonnees
!!                    geographiques en indices de grille I,J
!!
!!    PURPOSE
!!    -------
!
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
!!      Module MODD_CONVIJ2XY
!!
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!         XXX : XXHAT coordinate values for all the MESO-NH grids
!!         XXY : XYHAT                      "
!!
!!      Module MODE_GRIDPROJ   
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
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/04/99
!!      Updated   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_GRIDPROJ
USE MODD_COORD
USE MODD_FILES_DIACHRO
USE MODD_CONF
USE MODD_GRID
USE MODD_DIM1
USE MODD_GRID1
USE MODD_ALLOC_FORDIACHRO
USE MODD_RESOLVCAR
! Utilisation du meme module pour les operations inverses
USE MODD_CONVIJ2XY
USE MODD_PARAMETERS
USE MODI_RESOLVXISOLEV
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!
CHARACTER(LEN=800) :: HCARIN
!
!*       0.2   Local variables
!
INTEGER           :: IMGRID, J, I, JM
INTEGER           :: II, IJ, IIM, IJM
INTEGER           :: IIU, IJU, ICONVXY2IJ
REAL              :: ZLAT,ZLON,ZX,ZY
!REAL,DIMENSION(:),ALLOCATABLE :: ZCONVLAT, ZCONVLON
!
REAL,DIMENSION(100) :: ZIJ
CHARACTER(LEN=8)    :: YMGRID
!
!-------------------------------------------------------------------------------
!
!*      1.     
!              ----------------------------
!
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
CALL INI_CST
!
!
!*      1.1    
!
HCARIN=ADJUSTL(HCARIN)
if(nverbia >0)then
  print *,' **CONVXY2IJ HCARIN ',HCARIN
endif
IF(NBFILES == 0)THEN
  print *,' Vous devez ouvrir le fichier pour lequel vous demandez l''information avec _file1_...'
  print *,' puis entrer a nouveau votre directive '
  LPBREAD=.TRUE.
  RETURN
ENDIF
IF (LCARTESIAN) THEN
  print *,' In the cartesian geometry, reference latitude and longitude are the same for the whole domain:'
  print *,XLAT0,XLON0
  LPBREAD=.TRUE.
  RETURN
ENDIF

ICONVXY2IJ=INDEX(HCARIN,'CONVXY2IJ')
ZIJ(:)=9999.
CALL RESOLVXISOLEV(HCARIN(1:LEN_TRIM(HCARIN)),ICONVXY2IJ,ZIJ)
DO J=SIZE(ZIJ,1),1,-1
  IF(ZIJ(J) /= 9999.)THEN
    JM=J
    EXIT
  ENDIF
ENDDO
if(nverbia >0)then
  print *,' convxy2ij: ZIJ ',ZIJ(1:JM)
endif
ALLOCATE(XCONVI(JM/2))
ALLOCATE(XCONVJ(JM/2))
ALLOCATE(XCONVX(JM/2))
ALLOCATE(XCONVY(JM/2))
ALLOCATE(XCONVLAT(JM/2))
ALLOCATE(XCONVLON(JM/2))
!ALLOCATE(ZCONVLAT(JM/2*7))
!ALLOCATE(ZCONVLON(JM/2*7))
J=JM/2
XCONVLAT(1:J)=ZIJ(1:JM-1:2)
XCONVLON(1:J)=ZIJ(2:JM:2)
IF(NVERBIA > 0)THEN
  print *,' convxy2ij: XCONVLAT, XCONVLON'
  print *,XCONVLAT
  print *,XCONVLON
ENDIF
DO IMGRID=1,7
DO I=1,J
ZLAT=XCONVLAT(I)
ZLON=XCONVLON(I)

!
CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZX,ZY)

XCONVX(I)=ZX
XCONVY(I)=ZY

DO II=2,SIZE(XXX,1)
  IF(ZX >= XXX(II-1,IMGRID) .AND. ZX < XXX(II,IMGRID))THEN
    IIM=II-1
    EXIT
  ELSE
    IF(II == SIZE(XXX,1))THEN
      IIM=II
      EXIT
    ENDIF
  ENDIF
ENDDO
IF(IIM == SIZE(XXX,1))THEN
  XCONVI(I)=IIM
ELSE
  XCONVI(I)=IIM+((ZX-XXX(IIM,IMGRID))/(XXX(IIM+1,IMGRID)-XXX(IIM,IMGRID)))
ENDIF

DO IJ=2,SIZE(XXY,1)
  IF(ZY >= XXY(IJ-1,IMGRID) .AND. ZY < XXY(IJ,IMGRID))THEN
    IJM=IJ-1
    EXIT
  ELSE
    IF(IJ == SIZE(XXY,1))THEN
      IJM=IJ
      EXIT
    ENDIF
  ENDIF
ENDDO
IF(IJM == SIZE(XXY,1))THEN
  XCONVJ(I)=IJM
ELSE
  XCONVJ(I)=IJM+((ZY-XXY(IJM,IMGRID))/(XXY(IJM+1,IMGRID)-XXY(IJM,IMGRID)))
ENDIF
!
!IF(I == 1)THEN
! ZCONVLAT(IMGRID*2-1)=ZLAT
! ZCONVLON(IMGRID*2-1)=ZLON
!ELSE
! ZCONVLAT(IMGRID*2)=ZLAT
! ZCONVLON(IMGRID*2)=ZLON
!ENDIF
IF(IMGRID == 1 .AND. I == 1)THEN

print *,' GRILLES *    LAT     *    LON     *      X      *      Y      *   I   *   J  '
print *,'******************************************************************************'
endif
IF(IMGRID == 1)THEN
YMGRID=' 1 et 4 '
ELSE IF(IMGRID == 2)THEN
YMGRID=' 2 et 6 '
ELSE IF(IMGRID == 3)THEN
YMGRID=' 3 et 7 '
ELSE IF(IMGRID == 5)THEN
YMGRID=' 5      '
ENDIF
IF(IMGRID == 1 .OR. IMGRID == 2 .OR. IMGRID == 3 .OR. IMGRID == 5)THEN
print 10,YMGRID,XCONVLAT(I),XCONVLON(I),XCONVX(I),XCONVY(I),XCONVI(I),XCONVJ(I)
print *,'------------------------------------------------------------------------------'
ENDIF
ENDDO
ENDDO
!if (nverbia > 0)then
!DO I=1,J*7
! ZLAT=ZCONVLAT(I)
! ZLON=ZCONVLON(I)
! CALL SM_XYHAT_S(XLATORI,XLONORI,ZLAT,ZLON,ZX,ZY)
! print *,' ZLAT=',ZLAT,' ZLON=',ZLON,' ZX=',ZX,' ZY=',ZY
!ENDDO
!endif
10 FORMAT(1X,A8,' * ',F10.6,' * ',F11.6,'* ',F10.0,'  * ',F10.0,'  *',F6.2,' *',F6.2)
DEALLOCATE(XCONVI)
DEALLOCATE(XCONVJ)
DEALLOCATE(XCONVX)
DEALLOCATE(XCONVY)
DEALLOCATE(XCONVLAT)
DEALLOCATE(XCONVLON)
!DEALLOCATE(ZCONVLAT)
!DEALLOCATE(ZCONVLON)

!
!
!------------------------------------------------------------------------------
!
!*      2.    EXIT
!             ----
!
!
RETURN
END SUBROUTINE CONVXY2IJ
