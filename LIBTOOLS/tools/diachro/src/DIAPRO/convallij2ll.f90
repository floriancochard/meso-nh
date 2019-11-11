!     ######spl
      SUBROUTINE CONVALLIJ2LL(HCARIN)
!     ###############################
!
!!****  *CONVALLIJ2LL* - Convertit des indices de grille I,J en coordonnees
!!                    conformes et coordonnees geographiques
!!                    sur l'ensemble du domaine (points de garde)
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
USE MODD_GRID, ONLY: XLONORI,XLATORI
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
CHARACTER(LEN=2400) :: HCARIN
!
!*       0.2   Local variables
!
INTEGER           :: IMGRID, J, I, JM, INUM, IRESP
INTEGER           :: IEGAL, IG, ITER, II, IDEB, IFIN
INTEGER           :: IIU, IJU, ICONVALLIJ2LL, ICONVI, ICONVJ
REAL,DIMENSION(:,:),ALLOCATABLE :: ZLAT,ZLON
REAL              :: ZX1, ZY1, ZLA, ZLO
!REAL,DIMENSION(:),ALLOCATABLE :: ZCONVLAT, ZCONVLON
!
REAL,DIMENSION(100) :: ZIJ
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
  print *,' **CONVALLIJ2LL HCARIN ',TRIM(HCARIN)
endif
IF(NBFILES == 0)THEN
  print *,' Vous devez ouvrir le fichier pour lequel vous demandez l''information avec _file1_...'
  print *,' puis entrer a nouveau votre directive '
  LPBREAD=.TRUE.
  RETURN
ENDIF
CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
IF(IRESP /= 0)THEN
  CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
  OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
  PRINT '('' CONVALLIJ2LL --> Les valeurs seront mises dans le fichier FICVAL '')'
ENDIF
IF (LCARTESIAN) THEN
  print *,' In the cartesian geometry, reference latitude and longitude are the same for the whole domain:'
  print *,XLAT0,XLON0
  LPBREAD=.TRUE.
  RETURN
ENDIF

ICONVALLIJ2LL=INDEX(HCARIN,'CONVALLIJ2LL')
ZIJ(:)=9999.
IEGAL=INDEX(HCARIN,'=')
IF(IEGAL == 0)THEN
  JM=4
  ZIJ(1)=1.
  ZIJ(2)=2.
  ZIJ(3)=3.
  ZIJ(4)=5.
ELSE
CALL RESOLVXISOLEV(HCARIN(1:LEN_TRIM(HCARIN)),ICONVALLIJ2LL,ZIJ)
DO J=SIZE(ZIJ,1),1,-1
  IF(ZIJ(J) /= 9999.)THEN
    JM=J
    EXIT
  ENDIF
ENDDO
ENDIF
if(nverbia >0)then
  print *,' ZIJ ',ZIJ(1:JM)
endif
ALLOCATE(XCONVI(IIU))
ALLOCATE(XCONVJ(IJU))
ALLOCATE(ZLAT(IIU,IJU))
ALLOCATE(ZLON(IIU,IJU))
DO I=1,IIU
  XCONVI(I)=I
ENDDO
DO J=1,IJU
  XCONVJ(J)=J
ENDDO
IF(NVERBIA > 0)THEN
ENDIF
DO IG=1,JM

IMGRID=ZIJ(IG)
DO J=1,IJU
DO I=1,IIU
ICONVI=INT(XCONVI(I))
ICONVJ=INT(XCONVJ(J))
!IF(I < 5 .AND. J < 5)print *,' IMGRID,ICONVI,ICONVJ ',IMGRID,ICONVI,ICONVJ
ZX1=XXX(ICONVI,IMGRID)+(XXX(MIN(ICONVI+1,SIZE(XXX,1)),IMGRID)-XXX(ICONVI,IMGRID))*(XCONVI(I)-FLOAT(ICONVI))
ZY1=XXY(ICONVJ,IMGRID)+(XXY(MIN(ICONVJ+1,SIZE(XXY,1)),IMGRID)-XXY(ICONVJ,IMGRID))*(XCONVJ(J)-FLOAT(ICONVJ))
CALL SM_LATLON_S(XLATORI,XLONORI,ZX1,ZY1,ZLA,ZLO)
ZLAT(I,J)=ZLA
ZLON(I,J)=ZLO
!IF(I < 5 .AND. J < 5)print *,' ZLA,ZLO ',ZLA,ZLO
ENDDO
ENDDO
ITER=IIU/3
IF(ITER*3 < IIU)ITER=ITER+1
IF(IMGRID == 1 .OR. IMGRID == 4)THEN
WRITE(INUM,*)' FICHIER: ',CFILEDIAS(NUMFILECUR)
WRITE(INUM,*)' GRILLES N: 1  et  4   ITER:',ITER,' CONVERSION I,J -> LAT,LON '
ELSE IF(IMGRID == 2 .OR. IMGRID == 6)THEN
WRITE(INUM,*)' FICHIER: ',CFILEDIAS(NUMFILECUR)
WRITE(INUM,*)' GRILLES N: 2  et  6   ITER:',ITER,' CONVERSION I,J -> LAT,LON '
ELSE IF(IMGRID == 3 .OR. IMGRID == 7)THEN
WRITE(INUM,*)' FICHIER: ',CFILEDIAS(NUMFILECUR)
WRITE(INUM,*)' GRILLES N: 3  et  7   ITER:',ITER,' CONVERSION I,J -> LAT,LON '
ELSE IF(IMGRID == 5)THEN
WRITE(INUM,*)' FICHIER: ',CFILEDIAS(NUMFILECUR)
WRITE(INUM,*)' GRILLE N: 5   ITER:',ITER,' CONVERSION I,J -> LAT,LON '
ENDIF
WRITE(INUM,'(''  niinf'',i4,'' njinf'',i4,'' nisup'',i4,'' njsup'',i4)')LBOUND(ZLAT,1),&
LBOUND(ZLAT,2),IIU,IJU
WRITE(INUM,'(''  NBVAL en I '',i4,'' NBVAL en J '',i4)')IIU,IJU
DO I=1,ITER
  IF(I == 1)THEN
    IDEB=1; IFIN=3
  ELSE
    IDEB=IFIN+1; IFIN=IFIN+3
  IF(I == ITER)THEN
    IFIN=IIU
  ENDIF
  ENDIF
  WRITE(INUM,'(1X,78(1H*))')
  WRITE(INUM,'(''    I-> '',7X,I4,9X,2(9X,I4,9X))')(/(II,II=IDEB,IFIN)/)
  WRITE(INUM,'('' J      '',3X,''Lat  ,  Lon'',6X,2(8X,''Lat  ,  Lon'',6X))')
  WRITE(INUM,'(1X,78(1H*))')
  DO J=IJU,1,-1
!   WRITE(INUM,'(I3,3('' *'',F10.5,'' ,'',F10.5,1X))')J,ZLAT(IDEB,J), &
!   ZLON(IDEB,J),ZLAT(IDEB+1,J),ZLON(IDEB+1,J),ZLAT(IDEB+2,J),ZLON(IDEB+2,J)
    WRITE(INUM,'(I4,3('' *'',F10.5,'' ,'',F10.5,1X))')J,(ZLAT(II,J), &
    ZLON(II,J),II=IDEB,IFIN)
  ENDDO
ENDDO
ENDDO
DEALLOCATE(XCONVI)
DEALLOCATE(XCONVJ)
DEALLOCATE(ZLAT)
DEALLOCATE(ZLON)
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
END SUBROUTINE CONVALLIJ2LL
