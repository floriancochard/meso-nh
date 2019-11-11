!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:./s.ini2lalo.f90, Version:1.5, Date:03/06/05, Last modified:01/10/23
!-----------------------------------------------------------------
!     ######spl
MODULE MODI_INI2LALO
!###################
!
INTERFACE
      SUBROUTINE INI2LALO(PLATLON,KNX,KNY,            &
                          KIDEB,KIFIN,KJDEB,KJFIN,PDLON,PDLAT)
!
REAL, DIMENSION(4), INTENT(INOUT) :: PLATLON ! NSWE target domain bounds (deg)
INTEGER, INTENT(OUT) :: KNX,KNY        ! NUMBER OF TARGET POINTS IN X,Y
INTEGER, INTENT(IN), OPTIONAL :: KIDEB,KIFIN ! limites du
INTEGER, INTENT(IN), OPTIONAL :: KJDEB,KJFIN !zoom eventuel
REAL,    INTENT(OUT),OPTIONAL :: PDLON,PDLAT ! resolutions in LOn-LAt computed
!
END SUBROUTINE INI2LALO
END INTERFACE
END MODULE MODI_INI2LALO
!     ########################################################
      SUBROUTINE INI2LALO(PLATLON,KNX,KNY,            &
                          KIDEB,KIFIN,KJDEB,KJFIN,PDLON,PDLAT)
!     ########################################################
!
!!    PURPOSE
!!    -------
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,        ONLY : XRADIUS,XPI
USE MODD_PARAMETERS, ONLY : JPHEXT
USE MODD_GRID,       ONLY : XLAT0,XLATORI,XLONORI
USE MODD_DIM1,       ONLY : NIMAX,NJMAX
USE MODD_GRID1,      ONLY : XLAT,XMAP,XXHAT,XYHAT
!
USE MODE_GRIDPROJ
!
IMPLICIT NONE
!
!*       0.1   Arguments
!
REAL, DIMENSION(4), INTENT(INOUT) :: PLATLON ! NSWE target domain bounds (millideg)
INTEGER,            INTENT(OUT)   :: KNX,KNY ! nb of target points
INTEGER, INTENT(IN), OPTIONAL :: KIDEB,KIFIN ! limites du
INTEGER, INTENT(IN), OPTIONAL :: KJDEB,KJFIN !zoom eventuel
REAL   , INTENT(OUT),OPTIONAL :: PDLON,PDLAT ! resolutions in LOn-LAt computed
!
!*       0.2   Local variables
!
REAL :: ZLONW,ZLONE,ZLATN,ZLATS  ! LAT/LON rounded to nearest millidegree
REAL :: ZDX,ZDY                  ! increments in LAT/LON
REAL :: ZLATM                    ! extreme latitude of input domain
REAL :: ZLA,ZLO,ZI,ZJ,ZXHAT,ZYHAT
INTEGER, DIMENSION(2) :: IMAP
INTEGER :: II,IJ,IIN
INTEGER :: JX,JY
!
!------------------------------------------------------------------------------
!
!*       1.    CHECK Lat/Lon DOMAIN
!              --------------------
!
ZLATN=PLATLON(1) ; ZLATS=PLATLON(2)
ZLONW=PLATLON(3) ; ZLONE=PLATLON(4) 
!
! round to nearest millidegree, longitudes in (0..360) interval
ZLATN=REAL(NINT(ZLATN))
ZLATS=REAL(NINT(ZLATS))
ZLONW=MOD(ZLONW,360000.) ; ZLONE=MOD(ZLONE,360000.)
ZLONW=REAL(NINT(ZLONW))
ZLONE=REAL(NINT(ZLONE))
PLATLON(1)=ZLATN
PLATLON(2)=ZLATS
PLATLON(3)=ZLONW
PLATLON(4)=ZLONE
!
! check if domain is well-defined
IF(ABS(ZLATN)>90000.) THEN
  PRINT*, 'INI2LALO: Bad N latitude - abort: ZLATN=',ZLATN
  STOP
END IF
IF(ABS(ZLATS)>90000) THEN
  PRINT*, 'INI2LALO: Bad S latitude - abort: ZLATS=',ZLATS
  STOP
END IF
IF(ZLATN<=ZLATS) THEN
  PRINT*, 'Bad latitude interval - abort'
  STOP
END IF
!
! compute optimum resolution
IF (PRESENT(KIDEB)) THEN
  IMAP=MINLOC(XMAP(KIDEB:KIFIN,KJDEB:KJFIN)) 
ELSE
  IMAP=MINLOC(XMAP(:,:)) 
END IF
ZLATM=XLAT(IMAP(1),IMAP(2))
ZDX=(XXHAT(3)-XXHAT(2))*180./XPI&
                       /(XRADIUS*COS(ZLATM*XPI/180.))
ZDY=(XYHAT(3)-XYHAT(2))*180./XPI/XRADIUS
print*, 'INI2LALO: equivalent resolution in lat ',ZLATM,IMAP
print*, '        ',ZDX,ZDY
WRITE(6,'(A,I4,1X,I4,A,F6.1,A)')'INI2LALO: point where map scale is minimum ', &
                               IMAP,' (lat ',ZLATM,')'
PRINT*,'equivalent resolution in lon. and lat.: ' ,ZDX,ZDY
!
! compute number of points and
! move E & S boundaries so that lon/lat increments are in millidegrees
! (GRIB constraint)
KNX=NINT( (ZLONE-ZLONW)/(1000*ZDX) +1)
IF(KNX<0) KNX=NINT( (ZLONE-ZLONW+360000.)/(1000*ZDX) +1)
IF(ZDX/=REAL(NINT(ZDX*1000.))/1000.) THEN  ! need to fix longitude
  ZDX=REAL(NINT(ZDX*1000.))
  ZLONE=ZLONW+(KNX-1)*ZDX
  IF(ZLONE>360000.) ZLONE=ZLONE-360000.
  print*, 'INI2LALO: fixing E longitude to ',ZLONE
  PLATLON(4)=ZLONE
  ZDX=ZDX/1000.
ENDIF
!
KNY=NINT( (ZLATN-ZLATS)/(1000*ZDY) +1)
IF(ZDY/=REAL(NINT(ZDY*1000.))/1000.) THEN  ! need to fix latitude
  ZDY=REAL(NINT(ZDY*1000.))
  ZLATS=ZLATN-(KNY-1)*ZDY
  IF(ABS(ZLATS)>90000.) THEN
    STOP "TOO BIG DOMAIN in LATITUDE"
  END IF
  print*, 'INI2LALO: fixing S latitude to ',ZLATS
  PLATLON(2)=ZLATS
    ZDY=ZDY/1000.
ENDIF
!
IF(PRESENT(PDLON))THEN
  PDLON=ZDX
  PDLAT=ZDY
END IF
!
print*, 'INI2LALO: number of points in lon. and lat. domain:', KNX,KNY
IF (PRESENT(KIDEB)) THEN
  PRINT*, 'number of points of input domain (i,j,i*j): ',(KIFIN-KIDEB+1),&
         (KJFIN-KJDEB+1),(KIFIN-KIDEB+1)*(KJFIN-KJDEB+1)
ELSE
  PRINT*, 'number of points of input domain (i,j,i*j): ',NIMAX,NJMAX,NIMAX*NJMAX
END IF
PRINT*, 'number of points of lon.-lat. domain(x,y,x*y):', KNX,KNY,KNX*KNY
!
! check if target domain is inside file domain
IIN=0
DO JY=1,KNY ; DO JX=1,KNX
  ZLO=MOD(ZLONW/1000.+ZDX*(JX-1),360.)
  ZLA=ZLATN/1000.-ZDY*(JY-1)          ! output has N->S scanning
  CALL SM_XYHAT(XLATORI,XLONORI,ZLA,ZLO,ZXHAT,ZYHAT)
  II=MAX(MIN(COUNT(XXHAT(:)<ZXHAT),NIMAX+JPHEXT),1+JPHEXT)
  IJ=MAX(MIN(COUNT(XYHAT(:)<ZYHAT),NJMAX+JPHEXT),1+JPHEXT)
  ZI=(ZXHAT-XXHAT(II))/(XXHAT(II+1)-XXHAT(II))+FLOAT(II)-1
  ZJ=(ZYHAT-XYHAT(IJ))/(XYHAT(IJ+1)-XYHAT(IJ))+FLOAT(IJ)-1
  !
  IF (      (ZI>=1.) .AND. (ZI<=NIMAX) &
     .AND.  (ZJ>=1.) .AND. (ZJ<=NJMAX) ) THEN
    IIN=IIN+1   ! points inside
  ENDIF
END DO ; END DO
PRINT*, 'number of points of lon.-lat. domain inside input file one:', IIN
!
!
END SUBROUTINE INI2LALO
