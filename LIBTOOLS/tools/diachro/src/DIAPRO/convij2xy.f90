!     ######spl
      MODULE  MODI_CONVIJ2XY
!     ######################
!
INTERFACE
!
SUBROUTINE CONVIJ2XY(HCARIN)
CHARACTER(LEN=*) :: HCARIN
END SUBROUTINE CONVIJ2XY
!
END INTERFACE
!
END MODULE MODI_CONVIJ2XY
!     ######spl
      SUBROUTINE CONVIJ2XY(HCARIN)
!     ##################
!
!!****  *CONVIJ2XY* - Convertit des indices de grille I,J en coordonnees
!!                    conformes et coordonnees geographiques
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
!!      Module MODD_CONIJ2XY
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
USE MODD_CONVIJ2XY
USE MODD_PARAMETERS
USE MODI_RESOLVXISOLEV
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!
CHARACTER(LEN=*)  :: HCARIN
!
!*       0.2   Local variables
!
INTEGER           :: JJLOOP,JILOOP ,IMGRID, J, JJ, I, JM
INTEGER           :: IIU, IJU, ICONVIJ2XY, ICONVI, ICONVJ
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
  print *,' **CONVIJ2XY HCARIN ',TRIM(HCARIN)
endif
IF(NBFILES == 0)THEN
  print *,' Vous devez ouvrir le fichier pour lequel vous demandez l''information avec _file1_...'
  print *,' puis entrer a nouveau votre directive '
  LPBREAD=.TRUE.
  RETURN
ENDIF
ICONVIJ2XY=INDEX(HCARIN,'CONVIJ2XY')
ZIJ(:)=9999.
CALL RESOLVXISOLEV(HCARIN(1:LEN_TRIM(HCARIN)),ICONVIJ2XY,ZIJ)
DO J=SIZE(ZIJ,1),1,-1
  IF(ZIJ(J) /= 9999.)THEN
    JM=J
    EXIT
  ENDIF
ENDDO
if(nverbia >0)then
  print *,' ZIJ ',ZIJ(1:JM)
endif
ALLOCATE(XCONVIJ(JM))
ALLOCATE(XCONVI(JM/2))
ALLOCATE(XCONVJ(JM/2))
ALLOCATE(XCONVX(JM/2))
ALLOCATE(XCONVY(JM/2))
ALLOCATE(XCONVLAT(JM/2))
ALLOCATE(XCONVLON(JM/2))
!ALLOCATE(ZCONVLAT(JM/2*7))
!ALLOCATE(ZCONVLON(JM/2*7))
J=JM/2
XCONVIJ(1:JM)=ZIJ(1:JM)
XCONVI(1:J)=XCONVIJ(1:JM-1:2)
XCONVJ(1:J)=XCONVIJ(2:JM:2)
IF(NVERBIA > 0)THEN
  print *,' convij2xy: XCONVIJ,XCONVI,XCONVJ'
  print *,XCONVIJ
  print *,XCONVI,'  ',XCONVJ
ENDIF
!
DO IMGRID=1,7
DO I=1,J
ICONVI=INT(XCONVI(I))
ICONVJ=INT(XCONVJ(I))
XCONVX(I)=XXX(ICONVI,IMGRID)+(XXX(MIN(ICONVI+1,SIZE(XXX,1)),IMGRID)-XXX(ICONVI,IMGRID))*(XCONVI(I)-FLOAT(ICONVI))
XCONVY(I)=XXY(ICONVJ,IMGRID)+(XXY(MIN(ICONVJ+1,SIZE(XXY,1)),IMGRID)-XXY(ICONVJ,IMGRID))*(XCONVJ(I)-FLOAT(ICONVJ))
ZX=XCONVX(I); ZY=XCONVY(I)
IF (.NOT. LCARTESIAN) THEN
  CALL SM_LATLON_S(XLATORI,XLONORI,ZX,ZY,ZLAT,ZLON)
  XCONVLAT(I)=ZLAT
  XCONVLON(I)=ZLON
  !IF(I == 1)THEN
  !  ZCONVLAT(IMGRID*2-1)=ZLAT
  !  ZCONVLON(IMGRID*2-1)=ZLON
  !ELSE
  !  ZCONVLAT(IMGRID*2)=ZLAT
  !  ZCONVLON(IMGRID*2)=ZLON
  !ENDIF
  IF(IMGRID == 1 .AND. I == 1)THEN
print *,' GRILLES *   I   *   J   *      X      *      Y      *    LAT     *   LON  '
print *,'******************************************************************************'
  ENDIF
ELSE
  IF(IMGRID == 1 .AND. I == 1)THEN
print *,' GRILLES *   I   *   J   *      X      *      Y      '
print *,'*******************************************************'
  ENDIF
ENDIF
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
IF (.NOT. LCARTESIAN) THEN
  print 10,YMGRID,XCONVI(I),XCONVJ(I),XCONVX(I),XCONVY(I),XCONVLAT(I),XCONVLON(I)
ELSE
  print 20,YMGRID,XCONVI(I),XCONVJ(I),XCONVX(I),XCONVY(I)
ENDIF  
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
10 FORMAT(1X,A8,' *',F6.2,' *',F6.2,' * ',F10.0,'  * ',F10.0,'  *',F10.6,' *',F11.6)
20 FORMAT(1X,A8,' *',F6.2,' *',F6.2,' * ',F10.0,'  * ',F10.0)
DEALLOCATE(XCONVIJ)
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
END SUBROUTINE CONVIJ2XY
