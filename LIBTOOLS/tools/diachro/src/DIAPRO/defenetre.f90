!     ######spl
      SUBROUTINE DEFENETRE
!     ####################
!
!!****  *DEFENETRE* - Defines the display window for a cartesian model
!!
!!    PURPOSE
!!    -------
!       Defines the display window in the cartesian case for horizontal
!     cross-sections
!
!!**  METHOD
!!    ------
!!      NCAR routines are called to select a display window 
!!    corresponding to the post-processed section of the model 
!!    arrays (NIINFxNISUP).(NJINFxNJSUP)
!!     
!!
!!    EXTERNAL
!!    --------
!!      SET      : defines NCAR window and viewport in normalized and user
!!                 coordinates
!!      LABMOD   : defines axis label format
!!      GRIDAL   : draws axis divisions and ticks
!!      PERIM    : draws a perimeter box for the current plot
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!       XXX,XXY  : coordinate values for all the MESO-NH grids
!!
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!         NMGRID      : Current MESO-NH grid indicator
!!
!!      Module MODD_DIM1 : contains dimensions of data arrays
!!         NIINF, NISUP : lower and upper bounds of arrays
!!                        to be plotted in x direction
!!         NJINF, NJSUP : lower and upper bounds of arrays
!!                        to be plotted in y direction
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
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   13/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_COORD
USE MODD_NMGRID
USE MODD_RESOLVCAR
USE MODD_DIM1
USE MODD_CTL_AXES_AND_STYL
USE MODN_NCAR
!
IMPLICIT NONE
!
REAL :: ZWL, ZWR, ZWB, ZWT, ZDIFWLR, ZDIFWBT, ZDIFVLR, ZDIFVBT
REAL :: ZXMIN, ZXMAX,ZYMIN, ZYMAX
REAL :: ZVL, ZVR, ZVB, ZVT
REAL :: ZDIFVPTLR, ZDIFVPTBT
REAL :: ZI, ZJ, ZX, ZY, ZSZ, ZPOS, ZCENT
INTEGER :: IPOS, ICONVI, ICONVJ, JCAR, ICOLS, ICOLN
CHARACTER(LEN=10) :: FORMAX, FORMAY
CHARACTER(LEN=20) :: YNOM
CHARACTER(LEN=1)  :: YSYMB
!
!-------------------------------------------------------------------------------
!
!*       1.    DISPLAY WINDOW SETTING AND DRAWING
!              ----------------------------------
!
ZWL=XXX(NIINF,NMGRID)
ZWR=XXX(NISUP,NMGRID)
ZWB=XXY(NJINF,NMGRID)
ZWT=XXY(NJSUP,NMGRID)
ZXMIN=ZWL; ZXMAX=ZWR; ZYMIN=ZWB; ZYMAX=ZWT
!
ZDIFWLR=ZWR-ZWL
ZDIFWBT=ZWT-ZWB
if(nverbia > 0)then
print *,' defenetre ENTREE NMGRID NIINF,NISUP,NJINF,NJSUP,ZDIFWLR,ZDIFWBT'
print *,NMGRID,NIINF,NISUP,NJINF,NJSUP,ZDIFWLR,ZDIFWBT
print *,'ZWL,ZWR,ZWB,ZWT, ',ZWL,ZWR,ZWB,ZWT
endif
!
IF(LVPTUSER)THEN
  ZDIFVPTBT=XVPTT-XVPTB
  ZDIFVPTLR=XVPTR-XVPTL
  IF(ZDIFVPTBT >= ZDIFVPTLR*ZDIFWBT/ZDIFWLR)THEN
    ZDIFVBT=ZDIFVPTLR*ZDIFWBT/ZDIFWLR
    ZVB=XVPTB+ABS(ZDIFVPTBT-ZDIFVBT)/2.
!   XVPTB=XVPTB+ABS(ZDIFVPTBT-ZDIFVBT)/2.
    ZVT=XVPTT-ABS(ZDIFVPTBT-ZDIFVBT)/2.
!   XVPTT=XVPTT-ABS(ZDIFVPTBT-ZDIFVBT)/2.
    ZVL=XVPTL; ZVR=XVPTR
  ELSE
    ZDIFVLR=ZDIFVPTBT*ZDIFWLR/ZDIFWBT
    ZVL=XVPTL+ABS(ZDIFVPTLR-ZDIFVLR)/2.
!   XVPTL=XVPTL+ABS(ZDIFVPTLR-ZDIFVLR)/2.
    ZVR=XVPTR-ABS(ZDIFVPTLR-ZDIFVLR)/2.
!   XVPTR=XVPTR-ABS(ZDIFVPTLR-ZDIFVLR)/2.
    ZVB=XVPTB; ZVT=XVPTT
  ENDIF
if(nverbia > 0)then
print *,'ZVL,ZVR,ZVB,ZVT LVPTUSER=T, ',ZVL,ZVR,ZVB,ZVT
endif
ELSE
  IF(ZDIFWLR.GT.ZDIFWBT)THEN
    ZVL=.1
    ZVR=.90
  ! ZVR=.95
    ZDIFVLR=ZVR-ZVL
    ZDIFVBT=ZDIFVLR/ZDIFWLR*ZDIFWBT
    ZVB=(1.-ZDIFVBT)/2.
    ZVT=1.-ZVB
if(nverbia > 0)then
print *,'ZVL,ZVR,ZVB,ZVT, ',ZVL,ZVR,ZVB,ZVT
endif
  ELSE
    ZVB=.1
    ZVT=.90
  ! ZVT=.95
    ZDIFVBT=ZVT-ZVB
    ZDIFVLR=ZDIFVBT/ZDIFWBT*ZDIFWLR
    ZVL=(1.-ZDIFVLR)/2.
    ZVR=1.-ZVL
  END IF
END IF
!
if(nverbia > 0)then
print *,' defenetre ZVL,ZVR,ZVB,ZVT ',ZVL,ZVR,ZVB,ZVT
endif
!!!!!!!!!!!!!!! Sept 99
IF(LINDAX)THEN
if(nverbia > 0)then
print *, '***********DEFENETRE NIINF ...',NIINF,NISUP,NJINF,NJSUP
endif
CALL SET(ZVL,ZVR,ZVB,ZVT,FLOAT(NIINF),FLOAT(NISUP),FLOAT(NJINF),FLOAT(NJSUP),1)    ! Sets NCAR user coordinates
FORMAX='          '
IF(LFMTAXEX)THEN
  FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
ELSE
  FORMAX='(F5.1)'
ENDIF
FORMAY='          '
IF(LFMTAXEY)THEN
  FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
ELSE
  FORMAY='(F5.1)'
ENDIF

CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) ! Sets axis label formats
!CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0) ! Sets axis label formats
!CALL LABMOD('(F5.1)','(F5.1)',0,0,10,10,0,0,0) ! Sets axis label formats

ELSE

CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)    ! Sets NCAR user coordinates
!                                              ! and normalized coordinates
FORMAX='          '
IF(LFMTAXEX)THEN
  FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
ELSE
  FORMAX='(F8.0)'
ENDIF
FORMAY='          '
IF(LFMTAXEY)THEN
  FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
ELSE
  FORMAY='(F8.0)'
ENDIF
CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0) ! Sets axis label formats
!CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0) ! Sets axis label formats
!CALL LABMOD('(F8.0)','(F8.0)',0,0,10,10,0,0,0) ! Sets axis label formats
ENDIF
!!!!!!!!!!!!!!! Sept 99
!CALL GRIDAL(1,1,1,1,1,1,5,0.,0.)
CALL GASETI('LTY',1)                           ! Labels printed by PLCHHQ
IF(LINDAX)THEN
! Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,0,0,5,0.,0)   ! Draws axis tiks and labels
  ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,0,1,5,0.,0)   ! Draws axis tiks and labels
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,1,0,5,0.,0)   ! Draws axis tiks and labels
  ELSE
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,1,1,5,0.,0)   ! Draws axis tiks and labels
  ENDIF
! Avril 2002
ELSE
! Avril 2002
  IF(LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,0,0,5,0.,0)   ! Draws axis tiks and labels
  ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,0,1,5,0.,0)   ! Draws axis tiks and labels
  ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,1,0,5,0.,0)   ! Draws axis tiks and labels
  ELSE
    CALL GRIDAL(NCHITVXMJ,NCHITVXMN,NCHITVYMJ,NCHITVYMN,1,1,5,0.,0)   ! Draws axis tiks and labels
  ENDIF
! Avril 2002
ENDIF
!CALL GRIDAL(5,0,4,0,1,1,5,0.,0)                ! Draws axis tiks and labels
CALL PERIM(1,0,1,0)                            ! Draws perimeter box
!
!!!!!!!!!!!!!!! Sept 99
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,1)    ! Sets NCAR user coordinates
!!!!!!!!!!!!!!! Sept 99
!!!!!!!!!!!!!!! NOv  R2000

! deplace en juillet 2010 dans bcgrd_fordiachro.f90 par G. TANGUY
!if(nverbia > 0)then
!  print *,' **defenetre NIJCAR ',NIJCAR
!endif
!IF(NIJCAR.GE.1)THEN
!  IF(.NOT.LCOLAREA .AND. .NOT.LCOLINE)THEN
!    call tabcol_fordiachro
!  ENDIF
!! IF(LUMVM .OR. LUTVT .AND. NSUPERDIA == 1)THEN
!!   call tabcol_fordiachro
!! ENDIF
!  DO JCAR=1,NIJCAR
!    ZI=XICAR(JCAR)
!    ZJ=XJCAR(JCAR)
!    print *,' **defenetre ZI,ZJ ',ZI,ZJ
!    YSYMB=CSYMCAR(JCAR)
!    ZPOS=XPOSNOM(JCAR)
!    ICOLS=ICOLSYM(JCAR)
!    ICOLN=ICOLNOM(JCAR)
!    IF(XSZSYM(JCAR) /= 0.)THEN
!      ZSZ=XSZSYM(JCAR)
!      IF(ZSZ == 9999.)ZSZ=.012
!    ELSE
!      ZSZ=.012
!    ENDIF
!    ICONVI=INT(ZI)
!    ICONVJ=INT(ZJ)
!    if(nverbia > 0)then
!    print *,' **defenetre ICONVI, ICONVJ ',ICONVI,ICONVJ
!    endif
!    ZX=XXX(ICONVI,NMGRID)+(XXX(MIN(ICONVI+1,SIZE(XXX,1)),NMGRID)-XXX(ICONVI,NMGRID))*(ZI-FLOAT(ICONVI))
!    ZY=XXY(ICONVJ,NMGRID)+(XXY(MIN(ICONVJ+1,SIZE(XXY,1)),NMGRID)-XXY(ICONVJ,NMGRID))*(ZJ-FLOAT(ICONVJ))
!    if(nverbia > 0)then
!    print *,' **defenetre ZX,ZY ',ZX,ZY
!    endif
!!   CALL SM_XYHAT_S(XLATOR,XLONOR,ZLAT,ZLON,ZU,ZV)
!    CALL PCSETI('OC',ICOLS)
!    IF(YSYMB == '.')THEN
!      CALL NGWSYM('N',8,ZX,ZY,ZSZ,ICOLS,0)
!    ELSE
!      CALL PCSETI('OF',2)
!      CALL PCSETR('OL',1.5)
!      CALL PLCHHQ(ZX,ZY,YSYMB,ZSZ,0.,0.)
!      CALL PCSETI('OF',0)
!      CALL PCSETR('OL',0.)
!    ENDIF
!    CALL PCSETI('OC',1)
!    IF(XSZNOM(JCAR) /= 0.)THEN
!      ZSZ=XSZNOM(JCAR)
!      IF(ZSZ == 9999.)ZSZ=.012
!    ELSE
!      ZSZ=.012
!    ENDIF
!    IPOS=ZPOS
!!   print *,' ZSZ NOM ',ZSZ
!    SELECT CASE(IPOS)
!      CASE(0)
!	ZCENT=-1.
!	ZX=ZX+ZSZ*1.1*(ZXMAX-ZXMIN)
!      CASE(45)
!	ZCENT=-1.
!	ZX=ZX+ZSZ*1.0*(ZXMAX-ZXMIN)
!	ZY=ZY+ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!      CASE(90)
!	ZCENT=0.
!	ZY=ZY+ZSZ*1.5*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!      CASE(135)
!	ZCENT=1.
!	ZX=ZX-ZSZ*1.0*(ZXMAX-ZXMIN)
!	ZY=ZY+ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!      CASE(180)
!	ZCENT=1.
!	ZX=ZX-ZSZ*1.1*(ZXMAX-ZXMIN)
!      CASE(225)
!	ZCENT=1.
!	ZX=ZX-ZSZ*1.0*(ZXMAX-ZXMIN)
!	ZY=ZY-ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!      CASE(270)
!	ZCENT=0.
!	ZY=ZY-ZSZ*1.5*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!      CASE(315)
!	ZCENT=-1.
!	ZX=ZX+ZSZ*1.0*(ZXMAX-ZXMIN)
!	ZY=ZY-ZSZ*1.0*(MAX(ZXMAX-ZXMIN,ZYMAX-ZYMIN))
!    END SELECT 
!    IF(CNOMCAR(JCAR) /= ' ')THEN
!      YNOM=CNOMCAR(JCAR)
!      YNOM=ADJUSTL(YNOM)
!      CALL PCSETI('OF',2)
!      CALL PCSETI('OC',ICOLN)
!      CALL PCSETR('OL',1.5)
!!     CALL GSTXCI(ICOLN)
!!     CALL GSPLCI(ICOLN)
!      CALL PLCHHQ(ZX,ZY,YNOM(1:LEN_TRIM(YNOM)),ZSZ,0.,ZCENT)
!    ENDIF
!    CALL PCSETI('OF',0)
!    CALL PCSETR('OL',0.)
!    CALL PCSETI('OC',1)
!    CALL GSTXCI(1)
!  ENDDO
!ENDIF
!!!!!!!!!!!!!!! NOv  R2000
!-----------------------------------------------------------------------------
!
!*      2.   EXIT
!            ----
!
RETURN
END SUBROUTINE  DEFENETRE
