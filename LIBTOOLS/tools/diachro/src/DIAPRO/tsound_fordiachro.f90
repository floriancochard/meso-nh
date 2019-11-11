!     ######spl
SUBROUTINE TSOUND_FORDIACHRO(PPRES,PPTEMP,PPQV,PPU,PPV,KNN,HEADER,HTEXTE, &
		  OMXRAT, &
                  OMIXRAT,ODOFRAME,OSAMPLEUV)
!##########################################################################
!
!!****  *TSOUND_FORDIACHRO* - Emagram plotting routine
!!
!!    PURPOSE
!!    -------
!                                                                        
!        Plots soundings on a skew-T, log P Thermodynamic diagram
!       All units are in the international system.   
!
!!**  METHOD
!!    ------
!!       A standard sounding background is first drawn, and the current
!!      data are plotted on a skew-T, Log P diagram. Various functions
!!      are defined for scale conversion and moisture calculations.
!!
!!    EXPLICIT ARGUMENTS 
!!    ------------------
!!
!!       PRES      - Pressure array for thermodynamic data (Pascals)
!!       PTEMP     - Temperature array (Kelvin)
!!       PQV       - Water vapour mixing ratio (KG/KG)
!!       PU,PV     - Wind (M/S)
!!       KNN       - Number of data points
!!       HEADER    - 40 Character Header (var. name and misc.)
!!       HTEXTE    - Header with gridpoint location (grid indexes)
!!       OMXRAT    - Logical to control dew point line drawing 
!!       OMIXRAT   - Logical for water vapour variable mode selection
!!       ODOFRAME  - Logical for issuing a FRAME after plotting this emagram
!!       OSAMPLEUV - Logical for wind vector decimation
!!
!!    EXTERNAL
!!    --------
!!      OS   : computes the equivalent potential temperature
!!      TSA  : computes the pseudo-moist adiabat
!!      DEWP : computes the dew point
!!     
!!      Notice: two statement functions, ZFY, ZFX are also defined to 
!!              map the (T,P) points onto the user coordinates, and a 
!!              third one, ZCNP, is converts wind directions to the 
!!              meteorological standard.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!     MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!     NCAR Graphics Technical documentation, UNIX version 3.2,
!!     Scientific computing division, NCAR/UCAR, Boulder, USA.
!!      Volume 1: Fundamentals, Vers. 1, May 1993
!!      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
!!
!!     For thermodynamical functions, see for instance: 
!!      Bluestein H. B., 1992, "Synoptic-Dynamic Meteorology in mid-latitudes"
!!      Volume 1, Priciples of Kinematics and Dynamics, Section 4.3, p. 195,
!!      Oxford University Press.
!!
!!
!!    AUTHOR
!!    ------
!!      - Initial version Peridot TRACE Program, P.Bougeault *Meteo-France*,
!!      modified by R. Benoit (mc2, april 91) for the PYREX Oracle data base.
!!      - Present version J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   10/01/95
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TITLE  
USE MODD_TIT
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_RESOLVCAR
USE MODD_TYPE_AND_LH
USE MODN_NCAR
USE MODD_DIM1
USE MODD_RSISOCOL
USE MODD_PARAMETERS
USE MODI_FMREAD

IMPLICIT NONE
!
!*       0.1   Dummy arguments and results
!
INTEGER           :: KNN                      ! Number of data points
REAL,DIMENSION(:) :: PPRES, PPTEMP, PPQV, PPU, PPV ! Sounding state variables
REAL              :: PP, PT, PY, PA           ! Dummies for definitions
CHARACTER(LEN=*)  :: HEADER    ! Header containing variable name              
CHARACTER(LEN=*)  :: HTEXTE    ! Header containing sounding location
LOGICAL           :: OMXRAT    ! Logical keys pecifying whether moisture data
LOGICAL           :: OMIXRAT   ! are present, and if the moisture variable qv
                               ! contains mixing ratio or dewpoint temperature
LOGICAL           :: ODOFRAME  ! Logical for FRAME after plot control
LOGICAL           :: OSAMPLEUV ! Logical for wind plotting only

!
!*       0.2   Local variables
!
INTEGER,PARAMETER     ::  JPNWK=1000
INTEGER               ::  J, JJ, IK, JJJ, II, ID
INTEGER               ::  INUM, IRESP
INTEGER               ::  INC, IANGU, IENCD, ILEN, INEG  !,IPCK
INTEGER               ::  ILENT, ILEN2, JLOOP2, JLOOPT
INTEGER               ::  IKB, IKE, IKU
INTEGER               ::  IB, IE, IN
INTEGER               ::  IERR, ICOLI
INTEGER,DIMENSION(13) ::  IASF
REAL,DIMENSION(8,2)   ::  ZRAT
REAL,DIMENSION(15,2)  ::  ZTP
REAL,DIMENSION(81)    ::  ZSX, ZSY
REAL,DIMENSION(7)     ::  ZXB, ZYB
REAL,DIMENSION(9,2)   ::  ZPLN
REAL,DIMENSION(162)   ::  ZY45, ZDX, ZDY
REAL,DIMENSION(10)    ::  ZPLV
REAL                  ::  ZINT, ZVL, ZVR, ZVB, ZVT, ZWL, ZWR, ZWB, ZWT
REAL                  ::  ZXPOSTITT1, ZXYPOSTITT1
REAL                  ::  ZXPOSTITB1, ZXYPOSTITB1
REAL                  ::  ZXPOSTITB2, ZXYPOSTITB2
!
! Work vectors ZWORKS1...5 dimensioned to JPNWK=1000
! to receive high resolution souding as well.
!
!REAL,DIMENSION(JPNWK) :: ZWORKS1, ZWORKS2, ZWORKS3, ZWORKS4, ZWORKS5
!
REAL                  :: ZDTR, ZTS, ZTK, ZP, ZT, ZTD
REAL                  :: ZAOS, ZATSA, ZX1, ZX2, ZY1, ZY2, ZYD, ZYPD, ZXPD
REAL                  :: ZTX, ZX, ZY, ZDWPT, ZVSCALE, ZVVMAX, ZXM
REAL                  :: ZDYSMPL, ZYSMPL
REAL                  :: ZHA
REAL                  :: ZFX, ZFY, ZCNP
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: PRES, PTEMP, PQV, PU, PV
!
CHARACTER(LEN=2),DIMENSION(8)  :: YLRAT
CHARACTER(LEN=4)               :: YIT
CHARACTER(LEN=1)               :: YC1, Y1
CHARACTER(LEN=2)               :: YC2
CHARACTER(LEN=16)              :: YTEM
CHARACTER(LEN=80)              :: YTEM80
CHARACTER(LEN=19)              :: YGROUP
!
! Logical keys to activate wind, temperature plotting 
!
LOGICAL                        :: GDOTEMP, GDOUV, GDOUVM
!
! To prevent arrows overcrowding when high resolution data are used,
! a maximum number of arrows is set
!
INTEGER                        :: IMXSMPLUV=50
!
!*       0.3  Interface declarations
!
INTERFACE
  FUNCTION OS(PT,PP)
  REAL,INTENT(IN)                :: PT, PP
  REAL                           :: OS
  END FUNCTION OS
END INTERFACE
INTERFACE
  FUNCTION TSA(POS,PP)
  REAL,INTENT(IN)                :: POS, PP
  REAL                           :: TSA
  END FUNCTION TSA
END INTERFACE
INTERFACE
  FUNCTION DEWP(PQ,PP)
  REAL,INTENT(IN)                :: PQ, PP
  REAL                           :: DEWP
  END FUNCTION DEWP
END INTERFACE
INTERFACE
  SUBROUTINE WTSTR (PX,PY,CH,IS,IO,IC)
    CHARACTER*(*) CH
    REAL,INTENT(INOUT) :: PX,PY
    INTEGER :: IS,IO,IC
  END SUBROUTINE WTSTR
END INTERFACE
INTERFACE
  SUBROUTINE ECHELLE(KLEN,PHA)
    INTEGER, INTENT(OUT) :: KLEN
    REAL,    INTENT(OUT) :: PHA
  END SUBROUTINE ECHELLE
END INTERFACE
INTERFACE
  SUBROUTINE FLECHE(PX,PY,PU,PV,KLEN,PHA)
    INTEGER           :: KLEN
    REAL              :: PX, PY
    REAL              :: PU, PV  
    REAL              :: PHA
  END SUBROUTINE FLECHE
END INTERFACE
INTERFACE
  SUBROUTINE RESOLV_TIT(HTIT,HOUT)
    CHARACTER(LEN=*) :: HTIT, HOUT
  END SUBROUTINE RESOLV_TIT
END INTERFACE
!
!*      0.4   Statement function declarations
!
ZFY(PP) = 132.182-44.061*ALOG10(PP) ! Functions mapping the (T,P) values onto
ZFX(PT,PY) = 0.54*PT+0.90692*PY     ! the defined NCAR user coordinates
ZCNP(PA) = AMOD((450.-PA),360.)     ! Wind direction standardization
!
!------------------------------------------------------------------------------
!
!*      1.  BACKGROUND DATA TABLES SET UP
!           -----------------------------
!
!*      1.1  Defines an emagram  color table 
!

IASF(:)=1
CALL GSASF(IASF)
IF(LINVWB)THEN
CALL GSCR(1,1,0.,0.,0.)
CALL GSCR(1,0,1.,1.,1.)
ELSE
CALL GSCR(1,0,0.,0.,0.)
CALL GSCR(1,1,1.,1.,1.)
ENDIF
CALL GSCR(1,2,1.,0.,0.)
CALL GSCR(1,3,0.,1.,0.)
CALL GSCR(1,62,1.,.625,0.)

IKB=1+JPVEXT
IKU=NKMAX+2*JPVEXT
IKE=IKU-JPVEXT
YTEM80(1:LEN(YTEM80))=' '
CALL PCGETC('FC',Y1)
if(nverbia > 0)then
print *,' **tsou Y1 ',Y1
endif
CALL PCSETC('FC','?')
!
!*      1.2  Parameter checking
!
IF(ALLOCATED(PRES))THEN
  DEALLOCATE(PRES)
ENDIF
IF(ALLOCATED(PTEMP))THEN
  DEALLOCATE(PTEMP)
ENDIF
IF(ALLOCATED(PQV))THEN
  DEALLOCATE(PQV)
ENDIF
IF(ALLOCATED(PU))THEN
  DEALLOCATE(PU)
ENDIF
IF(ALLOCATED(PV))THEN
  DEALLOCATE(PV)
ENDIF
ALLOCATE(PRES(SIZE(PPRES)))
ALLOCATE(PTEMP(SIZE(PPTEMP)))
ALLOCATE(PQV(SIZE(PPQV)))
ALLOCATE(PU(SIZE(PPU)))
ALLOCATE(PV(SIZE(PPV)))
PRES(:)=PPRES(:)
PTEMP(:)=PPTEMP(:)
PQV(:)=PPQV(:)
PU(:)=PPU(:)
PV(:)=PPV(:)
PRINT *,' ********** TSOUND_FORDIACHRO'
IF(nverbia > 0)then
PRINT *,' PRES'
PRINT *,PRES
PRINT *,' PTEMP'
PRINT *,PTEMP
PRINT *,' PQV'
PRINT *,PQV
PRINT *,' PU'
PRINT *,PU
PRINT *,' PV'
PRINT *,PV
endif
PRINT *,' HEADER',HEADER ,'LEN ',LEN(HEADER),' LEN_TRIM ',LEN_TRIM(HEADER)
PRINT *,' HTEXTE',HTEXTE
PRINT *,' OMIXRAT ',OMIXRAT
PRINT *,' ODOFRAME ',ODOFRAME
PRINT *,' OSAMPLEUV ',OSAMPLEUV
IF(KNN.GT.JPNWK)THEN                                            ! if 1
  PRINT *,' Emagram TSOUND_FORDIACHRO... data overflows available arrays!'
  PRINT *,' KNN=',KNN,' when maximum allowed size is ',JPNWK,', return'
  RETURN
ENDIF                                                           ! endif 1
! ------nn <=> nwk -------
INC=KNN
  GDOTEMP=KNN.GT.0
  GDOUV=GDOTEMP
!! ESSAI
  IF(LNOUVRS)THEN
  GDOUV=.FALSE.
  ENDIF
!! ESSAI
  GDOUVM=GDOUV
!
!*     1.3   Data for constant mixing ratio lines 
!
ZRAT(1,1)=13.284
ZRAT(2,1)=8.91
ZRAT(3,1)=5.616
ZRAT(4,1)=1.944
ZRAT(5,1)=-1.782
ZRAT(6,1)=-4.698
ZRAT(7,1)=-9.234
ZRAT(8,1)=-14.796
ZRAT(1,2)=16.283
ZRAT(2,2)=12.125
ZRAT(3,2)=8.94
ZRAT(4,2)=5.45
ZRAT(5,2)=1.865
ZRAT(6,2)=-.858
ZRAT(7,2)=-5.313
ZRAT(8,2)=-10.686
!
YLRAT(1)='20'
YLRAT(2)='12'
YLRAT(3)=' 8'
YLRAT(4)=' 5'
YLRAT(5)=' 3'
YLRAT(6)=' 2'
YLRAT(7)=' 1'
YLRAT(8)='.4'
!                
!*    1.4   Data for constant temperature lines
!
ZTP(1,1)=1000.
ZTP(2,1)=1000.
ZTP(3,1)=1000.
ZTP(4,1)=1000.
ZTP(5,1)=1000.
ZTP(6,1)=1000.
ZTP(7,1)=1000.
ZTP(8,1)=1000.
ZTP(9,1)=855.
ZTP(10,1)=625.
ZTP(11,1)=459.
ZTP(12,1)=337.
ZTP(13,1)=247.
ZTP(14,1)=181.
ZTP(15,1)=132.
ZTP(1,2)=730.
ZTP(2,2)=580.
ZTP(3,2)=500.
ZTP(4,2)=430.
ZTP(5,2)=342.
ZTP(6,2)=251.
ZTP(7,2)=185.
ZTP(8,2)=135.
ZTP(9,2)=100.
ZTP(10,2)=100.
ZTP(11,2)=100.
ZTP(12,2)=100.
ZTP(13,2)=100.
ZTP(14,2)=100.
ZTP(15,2)=100.
!
!*    1.5   Data for constant pressure lines
!
ZPLV(1)=100.
ZPLV(2)=200.
ZPLV(3)=300.
ZPLV(4)=400.
ZPLV(5)=500.
ZPLV(6)=600.
ZPLV(7)=700.
ZPLV(8)=800.
ZPLV(9)=850.
ZPLV(10)=1000.
!                 
!*    1.6   Frame of the emagram plot
!
ZXB(1)= -19.
ZXB(2)=27.1
ZXB(3)=27.1
ZXB(4)=18.6
ZXB(5)=18.6
ZXB(6)=-19.
ZXB(7)=-19.
!               
ZYB(1)=0.
ZYB(2)=0.
ZYB(3)=9.
ZYB(4)=17.53
ZYB(5)=44.061
ZYB(6)=44.061
ZYB(7)=0.
!            
!*    1.7   Initial and final points of the
!*          constant pressure lines
!
!     IPCK = 0
!            
ZPLN(1,1)=-19.
ZPLN(2,1)=-19.
ZPLN(3,1)=-19.
ZPLN(4,1)=-19.
ZPLN(5,1)=-19.
ZPLN(6,1)=-19.
ZPLN(7,1)=-19.
ZPLN(8,1)=-19.
ZPLN(9,1)=-19.
ZPLN(1,2)=18.6
ZPLN(2,2)=18.6
ZPLN(3,2)=18.6
ZPLN(4,2)=18.6
ZPLN(5,2)=22.83
ZPLN(6,2)=26.306
ZPLN(7,2)=27.1
ZPLN(8,2)=27.1
ZPLN(9,2)=27.1
!                 
!*    1.8   Various constants
!
ZDTR = ATAN(1.)/45.
IANGU = 359.
!           
!-----------------------------------------------------------------------------
!
!*    2.    DRAWING THE BACKGROUND OF THE EMAGRAM PLOT
!           ------------------------------------------
!
!*    2.1   Draws outline of skew-T Log P diagram 
!
CALL GSTXCI(62)
CALL GSPLCI(62)
CALL GSFACI(62)                                   ! The NCAR user coordinate
CALL SET(.05,.95,.05,.95,-19.0,27.1,0.0,44.061,1) ! system is here set in
                                                  ! accordance with ZFY, ZFX
                                                  ! statement functions defined
                                                  ! above.
CALL CURVE(ZXB,ZYB,7)
!                       
!*    2.2   Draws satured adiabat. curves
!                              
CALL GSTXCI(2)
CALL GSPLCI(3)
CALL GSFACI(3)
ZTS = 32.
DO JJ = 1,7                                                      ! do 1
! CALL SETUSV ('IN',8000)
  CALL DASHDB(990)
  ZP = 1010.
  ZTK = ZTS+273.16
  ZAOS = OS(ZTK,1000.)
    DO J = 1,81                                                 ! do 2
      ZP = ZP-10.
      ZATSA = TSA(ZAOS,ZP)-273.16
      ZSY(J) = ZFY(ZP)
      ZSX(J) = ZFX(ZATSA,ZSY(J))
    ENDDO                                                       ! enddo 2
  CALL CURVED(ZSX,ZSY,81)
  IENCD = IFIX(ZTS)
  WRITE(YIT,100) IENCD
  YIT=ADJUSTL(YIT)
 100     FORMAT(I2)
  ZTS = ZTS-4.
  ZSY(81) = ZSY(81)+0.6
  CALL WTSTR(ZSX(81),ZSY(81),YIT(1:LEN_TRIM(YIT)),1,IANGU,0)
! CALL WTSTR(ZSX(81),ZSY(81),YIT(1:2),1,IANGU,0)
ENDDO                                                           ! enddo 1
!            
!*    2.3   Draws constant mixing ratio lines
!                       
DO J = 1,8                                                      ! do 1
! CALL SETUSV ('IN',8000)
  CALL DASHDB(29127)
  CALL LINED(ZRAT(J,1),-0.1,ZRAT(J,2),6.824)
  YIT(1:2) = YLRAT(J)
  YIT=ADJUSTL(YIT)
  ZY1=6.42
  CALL WTSTR(ZRAT(J,2),ZY1,YIT(1:LEN_TRIM(YIT)),1,IANGU,0)
! CALL WTSTR(ZRAT(J,2),1.42,YIT(1:2),1,IANGU,0)
! print *,' Mixing ratio lines'
ENDDO                                                           ! enddo 1
!            
!*    2.4   Draws constant temperature lines
!                      
CALL GSTXCI(62)
CALL GSPLCI(62)
CALL GSFACI(62)
ZT = 40.         
DO J = 1,15                                                     ! do 1
! CALL SETUSV('IN',8000)
  ZY1 = ZFY(ZTP(J,1))
  ZY2 = ZFY(ZTP(J,2))
  ZX1 = ZFX(ZT,ZY1)
  ZX2 = ZFX(ZT,ZY2)
  CALL LINE(ZX1,ZY1,ZX2,ZY2)
  IF(ZT.EQ.20.)GO TO 19
  IF(ABS(ZT) > 90)THEN
  ZX2 = ZX2+0.4
  ZY2 = ZY2+.441
  ELSEIF(ZT > -100 .AND. ZT < -30)THEN
  ZX2 = ZX2+0.4
  ZY2 = ZY2+.53
  ELSEIF(ZT > -40 .AND. ZT < 0)THEN
  ZX2 = ZX2+0.76
  ZY2 = ZY2+.453
  ELSE
  ZX2 = ZX2+0.88
! ZX2 = ZX2+0.4
  ZY2 = ZY2+.451
  ENDIF
! ZY2 = ZY2+.441
  IENCD = IFIX(ZT)
  WRITE(YIT,101) IENCD
  YIT=ADJUSTL(YIT)
  101     FORMAT(I4  )
  CALL WTSTR (ZX2,ZY2,YIT(1:LEN_TRIM(YIT)),2,45,0)
! CALL WTSTR (ZX2,ZY2,YIT(1:4),2,45,0)
! print *,' Temperature lines'
    19     ZT = ZT-10.
ENDDO                                                           ! enddo 1
!            
!*   2.5    Draws constant dry adiabat. curves
!                       
CALL GSTXCI(3)
CALL GSPLCI(3)
CALL GSFACI(3)
ZT = 51.
DO J = 1,162                                                    ! do 1
  ZY45(J) = 66.67*(5.7625544-ALOG(ZT+273.16))
  ZT = ZT-1.0
ENDDO                                                           ! enddo 1
ZT = 450.
ZTD = 52.
DO JJ = 1,20                                                     ! do 1
! CALL SETUSV('IN',8000)
  CALL DASHDB(13107)
  ZT = ZT-10.
  IK = 0
  ZYD = 66.67*(ALOG(ZT)-5.7625544)
    DO J = 1,162                                                ! do 2
      ZYPD = ZY45(J)+ZYD
      ZTX = ZTD-J
      IF(ZYPD.GT.44.061)EXIT
      IF(ZYPD.LT.0.0)CYCLE
      ZXPD = ZFX(ZTX,ZYPD)
      IF(ZXPD.LT.-19.0)EXIT
      IF(ZXPD.GT.27.1)CYCLE
      IF(ZXPD.GT.18.6.AND.ZT.GT.350.0)CYCLE
      IK = IK+1
      ZDX(IK) = ZXPD
      ZDY(IK) = ZYPD
    ENDDO                                                       ! enddo 2
  CALL CURVED(ZDX,ZDY,IK)
  IENCD = IFIX(ZT)
  WRITE(YIT,102) IENCD
  102     FORMAT(I3)
  CALL WTSTR(ZDX(IK-3),ZDY(IK-3),YIT(1:3),1,IANGU,0)
!print *,' constant dry adiabat. curves IK YIT ',IK,YIT
ENDDO                                                           ! enddo 1
!
!*     2.6    Draws constant pressure lines
!  
CALL GSTXCI(62)
CALL GSPLCI(62)
DO J = 1,10                                                     ! do 1
! CALL SETUSV('IN',8000)
  ZY1 = ZFY(ZPLV(J))
  IF(J.NE.1.AND.J.NE.10)CALL LINE(ZPLN(J,1),ZY1,ZPLN(J,2),ZY1)
  IENCD = IFIX(ZPLV(J) )
  WRITE(YIT,101) IENCD
  YIT=ADJUSTL(YIT)
  IF(J==10)THEN
    ZX1 = -20.4 
    CALL WTSTR(ZX1,ZY1,YIT(1:LEN_TRIM(YIT)),2,IANGU,0)
!   CALL WTSTR(-20.4,ZY1,YIT(1:4),2,IANGU,0)
  ELSE
    ZX1 = -20.3
    CALL WTSTR(ZX1,ZY1,YIT(1:LEN_TRIM(YIT)),2,IANGU,0)
!   CALL WTSTR(-20.7,ZY1,YIT(1:4),2,IANGU,0)
  ENDIF
! CALL WTSTR(-20.9,ZY1,YIT(1:4),1,IANGU,0)
ENDDO                                                           ! enddo 1
!
!*     2.7    Draws  ticks every 2 degrees at 500 MB
!
!CALL SETUSV('IN',8000)
ZY1 = 13.2627
ZY2 = 13.75
ZT = -52.
DO J = 1,31                                                     ! do 1
  ZT = ZT+2.
  IF(AMOD(ZT,10.).EQ.0.)CYCLE
  ZX1 = ZFX(ZT,ZY1)
  ZX2 = ZFX(ZT,ZY2)
  CALL LINE(ZX1,ZY1,ZX2,ZY2)
ENDDO                                                           ! enddo 1
!     IPCK = 1
!
!----------------------------------------------------------------------------
!
!*     3.     DRAWING THE SOUNDING DATA LINES ON THE SKEW-T-LOGP DIAGRAM
!             ----------------------------------------------------------
!
111 CONTINUE                !------111-------
!
!*     3.1   Plot Temperature and dewpoint curves
!
IANGU = 0.
!
CALL GSTXCI(1)
CALL GSPLCI(1)
CALL GSFACI(1)
!                                           
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
!Mars 2000
! Titre N1 BOTTOM
  ZXPOSTITB1=.002
  ZXYPOSTITB1=.005
  IF(XPOSTITB1 /= 0.)THEN
    ZXPOSTITB1=XPOSTITB1
  ENDIF
  IF(XYPOSTITB1 /= 0.)THEN
    ZXYPOSTITB1=XYPOSTITB1
  ENDIF
  CALL RESOLV_TIT('CTITB1',HEADER(1:100))
  IF(HEADER /= ' ')THEN
    IF(XSZTITB1 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,CLEGEND2(1:LEN_TRIM(CLEGEND2)),XSZTITB1,0.,-1.)
      if(nverbia > 0)then
      print *,' **tsound CLEGEND2 ',CLEGEND2(1:LEN_TRIM(CLEGEND2))
      endif
!     CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,CLEGEND2,XSZTITB1,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HEADER(1:LEN_TRIM(HEADER)),.007,0.,-1.)
      if(nverbia > 0)then
      print *,' **tsound HEADER ',HEADER(1:LEN_TRIM(HEADER))
      endif
!     CALL PLCHHQ(ZXPOSTITB1,ZXYPOSTITB1,HEADER,.007,0.,-1.)
    ENDIF
  ENDIF
!CALL PLCHHQ(0.002,0.005,HEADER,.007,0.,-1.)
!Mars 2000
!CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
CALL SET(.05,.95,.05,.95,-19.0,27.1,0.0,44.061,1)
CALL GSCLIP(0)
!print *,' ap GSCLIP'
CALL PLCHHQ(22.8,-1.,HTEXTE(1:LEN_TRIM(HTEXTE)),.01,0.,1.)
!CALL WTSTR(-19.,-1.,HEADER(1:60),1,IANGU,-1)
!!!!CALL WTSTR(22.8,-1.,HTEXTE(1:LEN_TRIM(HTEXTE)),1,IANGU,1)
!print *,' ap WTSTR '

IF(LRS1 .AND. CTYPE == 'CART')THEN
  ILENT=SIZE(XTRS,2)
  ILEN2=2
ELSE IF(LRS1 .AND. CTYPE == 'RSPL')THEN
  ILENT=SIZE(XTRS,1)
  ILEN2=2
ELSE
  ILENT=1
  ILEN2=1
ENDIF

! Memorisation des tableaux passes en arguments pour les restaurer par la suite
!DO JJJ=1,INC                                                     ! do 1
!  ZWORKS1(JJJ) = PRES(JJJ)
!  ZWORKS2(JJJ) = PTEMP(JJJ)
!  ZWORKS3(JJJ) = PQV(JJJ)
!  ZWORKS4(JJJ) = PU(JJJ)
!  ZWORKS5(JJJ) = PV(JJJ)
!ENDDO

DO JLOOP2=1,ILEN2

  DO JLOOPT=1,ILENT
!print *,' Boucle JLOOPT ',JLOOPT

    IF(JLOOP2 == 2 .OR. (JLOOP2 == 1 .AND. LRS1 .AND. JLOOPT >1))THEN

      IF(CTYPE == 'CART')THEN

        CTIMEC(1:LEN(CTIMEC))=' '
        WRITE(CTIMEC(1:8),'(F8.0)')XTIMRS(JLOOPT)
        CTIMEC(LEN_TRIM(CTIMEC)+1:LEN_TRIM(CTIMEC)+1)='s'
        CTIMEC=ADJUSTL(CTIMEC)
        IF(JLOOP2 == 1)THEN
          YTEM(1:LEN(YTEM))=' '
          YTEM=CTIMEC
          CTIMEC(1:LEN(CTIMEC))=' '
          YTEM=ADJUSTL(YTEM)
          IF(NVERBIA > 0)THEN
          print *,' YTEM ',YTEM
	  ENDIF
          WRITE(CTIMEC(1:1),'(I1)')JLOOPT
          CTIMEC(2:2)=' '
          CTIMEC(1+2:LEN_TRIM(YTEM)+2)=YTEM(1:LEN_TRIM(YTEM))
          IF(NVERBIA > 0)THEN
          print *,' CTIMEC ',CTIMEC
	  ENDIF
        ENDIF
      
      ELSE IF(CTYPE == 'RSPL')THEN

        CTIMECS(1:LEN(CTIMECS))=' '
        WRITE(CTIMECS(1:8),'(F8.0)')XTIMRS2(JLOOPT,1)
	CTIMECS=ADJUSTL(CTIMECS)

	IF(JLOOP2 == 1)THEN
	  YTEM(1:LEN(YTEM))=' '
	  YTEM=CTIMECS(1:LEN_TRIM(CTIMECS))
	  YTEM=ADJUSTL(YTEM)
          CTIMECS(1:LEN(CTIMECS))=' '
	  IF(NNST(JLOOPT) < 10)THEN
	    IN=1
	    WRITE(CTIMECS(1:IN),'(I1)')NNST(JLOOPT)
	  ELSE IF(NNST(JLOOPT) >= 10 .AND. NNST(JLOOPT) < 100)THEN
	    IN=2
	    WRITE(CTIMECS(1:IN),'(I2)')NNST(JLOOPT)
	  ELSE
	    IN=3
	    WRITE(CTIMECS(1:IN),'(I3)')NNST(JLOOPT)
	  ENDIF
	  IN=IN+1
	  CTIMECS(IN:IN)=' '
	  IN=IN+1
	  II=LEN_TRIM(YTEM)
	  CTIMECS(IN:IN+II-1)=YTEM(1:II)
	  IN=IN+II
	  CTIMECS(IN:IN)='-'
	  IN=IN+1
	  YTEM(1:II)=' '
	  WRITE(YTEM(1:8),'(F8.0)')XTIMRS2(JLOOPT,NST(JLOOPT))
	  YTEM=ADJUSTL(YTEM)
	  II=LEN_TRIM(YTEM)
	  CTIMECS(IN:IN+II-1)=YTEM(1:II)
	  IN=IN+II
	  CTIMECS(IN:IN)='s'

        ENDIF
        
      ENDIF

    ENDIF

    IF(JLOOP2 == 1 .AND. JLOOPT == 1)THEN
      IF(LRS1)THEN
! Cas LRS : CTIMEC est charge necessairement dans OPER_PROCESS

	SELECT CASE(CTYPE)

	  CASE('CART')
            CTIMEC(1:LEN(CTIMEC))=' '
            CTIMEC(1:3)='  ('
            WRITE(CTIMEC(4:11),'(F8.0)')XTIMRS(JLOOPT)
            CTIMEC(LEN_TRIM(CTIMEC)+1:LEN_TRIM(CTIMEC)+2)='s)'
	  CASE('RSPL')
            CTIMECS(1:LEN(CTIMECS))=' '
            CTIMECS(1:3)='  ('
            WRITE(CTIMECS(4:11),'(F8.0)')XTIMRS2(JLOOPT,1)
            CTIMECS(LEN_TRIM(CTIMECS)+1:LEN_TRIM(CTIMECS)+1)='-'
	    IN=LEN_TRIM(CTIMECS)+1
	    YTEM(1:LEN(YTEM))=' '
	    WRITE(YTEM(1:8),'(F8.0)')XTIMRS2(JLOOPT,NST(JLOOPT))
	    YTEM=ADJUSTL(YTEM)
	    II=LEN_TRIM(YTEM)
	    CTIMECS(IN:IN+II-1)=YTEM(1:II)
	    IN=IN+II
	    CTIMECS(IN:IN+1)='s)'

	END SELECT
      ENDIF

      II=LEN_TRIM(CLEGEND2)+1
!     print *,' **tsound II,len_trim(header) ',II,LEN_TRIM(HEADER)

      SELECT CASE(CTYPE)
	CASE('CART')
          CLEGEND2(II:II+LEN_TRIM(CTIMEC)-1)=CTIMEC(1:LEN_TRIM(CTIMEC))
	CASE('RSPL')
          CLEGEND2(II:II+LEN_TRIM(CTIMECS)-1)=CTIMECS(1:LEN_TRIM(CTIMECS))
      END SELECT
      if(nverbia > 0)then
      print *,' **tsound len_trim(clegend2),len_trim(header) ',LEN_TRIM(CLEGEND2),LEN_TRIM(HEADER)
      endif

      CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
      CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)

! Mars 2000
! Titre N2 BOTTOM
  ZXPOSTITB2=.002
  ZXYPOSTITB2=.025
  IF(XPOSTITB2 /= 0.)THEN
    ZXPOSTITB2=XPOSTITB2
  ENDIF
  IF(XYPOSTITB2 /= 0.)THEN
    ZXYPOSTITB2=XYPOSTITB2
  ENDIF
  CALL RESOLV_TIT('CTITB2',CLEGEND2)
  IF(CLEGEND2 /= ' ')THEN
    IF(XSZTITB2 /= 0.)THEN
      CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2(1:LEN_TRIM(CLEGEND2)),XSZTITB2,0.,-1.)
!     CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,XSZTITB2,0.,-1.)
    ELSE
      CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2(1:LEN_TRIM(CLEGEND2)),.007,0.,-1.)
!     CALL PLCHHQ(ZXPOSTITB2,ZXYPOSTITB2,CLEGEND2,.007,0.,-1.)
    ENDIF
  ENDIF
!     CALL PLCHHQ(0.002,0.025,CLEGEND2,.007,0.,-1.)
! Mars 2000
      CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
      IF(LDATFILE)CALL DATFILE_FORDIACHRO
!print *,' AP DATFILE2'
    ENDIF

IF (.NOT.GDOTEMP) GO TO 61

IF(LRS1)THEN
  SELECT CASE(CTYPE)
    CASE('CART')
      IB=IKB ; IE=IKE
      PRES(:)=XPRS(IB:IE,JLOOPT)
      PTEMP(:)=XTRS(IB:IE,JLOOPT)
      PQV(:)=XRVRS(IB:IE,JLOOPT)
      PU(:)=XURS(IB:IE,JLOOPT)
      PV(:)=XVRS(IB:IE,JLOOPT)
    CASE('RSPL')
      IB=1 ; IE=NST(JLOOPT)
      IF(ALLOCATED(PRES))THEN
	DEALLOCATE(PRES)
      ENDIF
      IF(ALLOCATED(PTEMP))THEN
	DEALLOCATE(PTEMP)
      ENDIF
      IF(ALLOCATED(PQV))THEN
	DEALLOCATE(PQV)
      ENDIF
      IF(ALLOCATED(PU))THEN
	DEALLOCATE(PU)
      ENDIF
      IF(ALLOCATED(PV))THEN
	DEALLOCATE(PV)
      ENDIF
      ALLOCATE(PRES(IE))
      ALLOCATE(PTEMP(IE))
      ALLOCATE(PQV(IE))
      ALLOCATE(PU(IE))
      ALLOCATE(PV(IE))
      PRES(:)=XPRS(JLOOPT,IB:IE)
      PTEMP(:)=XTRS(JLOOPT,IB:IE)
      PQV(:)=XRVRS(JLOOPT,IB:IE)
      PU(:)=XURS(JLOOPT,IB:IE)
      PV(:)=XVRS(JLOOPT,IB:IE)
      INC=SIZE(PRES)
  END SELECT
ENDIF
!
! Avril 99
!
IF(JLOOP2 == 1)THEN
IF(LPRINT)THEN
  CALL FMLOOK('FICVAL','FICVAL',INUM,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMATTR('FICVAL','FICVAL',INUM,IRESP)
    OPEN(UNIT=INUM,FILE='FICVAL',FORM='FORMATTED')
    PRINT '('' LPRINT=T --> Les valeurs seront mises dans le fichier FICVAL '')'
  ENDIF
  SELECT CASE(CTYPE)
    CASE('CART')
      IF(CGROUP == 'UM' .OR. CGROUP == 'VM' .OR. CGROUP == 'THM' .OR. &
      CGROUP == 'PABSM' .OR. CGROUP == 'RVM')THEN
	YGROUP='THM-PABSM-RVM-UM-VM'
      ELSE
	YGROUP='THT-PABST-RVT-UT-VT'
      ENDIF
      WRITE(INUM,'(''RS  '',''G:'',A19,25X,'' T:'',F8.0,''s'',''   (1-IKU)'')')YGROUP,&
&   XTIMRS(JLOOPT)

      WRITE(INUM,'(A19,20X,A4,6X,''NBVAL '',I5)')YGROUP,CTYPE,SIZE(XTRS,1)
      IF(XIRS /= -999.)THEN
        WRITE(INUM,'(''xirs'',F10.5,'' xjrs'',F10.5)')XIRS,XJRS
      ELSE
        WRITE(INUM,'(''nirs'',I5,'' njrs'',I5,'' (grille 1)'')')NIRS,NJRS
      ENDIF
      WRITE(INUM,'(1X,78(1H*))')
! JUin 2001 Ecriture des dates (Demande G.Jaubert ) si LPRDAT=T
  IF(LPRDAT)THEN
    IF(.NOT.ALLOCATED(XPRDAT))THEN
      print *,'**TSOUND XPRDAT NON ALLOUE.Dates non ecrites ds FICVAL .Prevenir J.Duron'
    ELSE
      WRITE(INUM,'(1X,75(1H*))')
      WRITE(INUM,'(1X,''    Dates courante   *     modele      *   experience    *      segment'')')
      WRITE(INUM,'(1X,'' J   An  M  J  Sec.  * An  M  J  Sec.  * An  M  J  Sec.  * An  M  J  Sec.'')')
      WRITE(INUM,'(1X,75(1H*))')
      DO J=1,SIZE(XPRDAT,2)
        WRITE(INUM,'(1X,I3,1X,3(I4,I3,I3,I6,'' *''),I4,I3,I3,I6)')J,INT(XPRDAT(:,J))
      ENDDO
    ENDIF
  ENDIF
! JUin 2001 Ecriture des dates 
      IF(CGROUP(LEN_TRIM(CGROUP):LEN_TRIM(CGROUP)) == 'M')THEN
        WRITE(INUM,'(5X,''K'',4X,''*  THM_RS *  PABSM  *'',7X,''RVM'',7X,&
        &''*    UM   *    VM'')')
      ELSE
        WRITE(INUM,'(5X,''K'',4X,''*  THT_RS *  PABST  *'',7X,''RVT'',7X,&
        &''*    UT   *    VT'')')
      ENDIF
      WRITE(INUM,'(1X,78(1H*))')
      DO J=SIZE(XTRS,1),1,-1
	IF(J == SIZE(XTRS,1))THEN
	  WRITE(INUM,'(''(IKU)'',I4,'' * '',F7.2,'' * '',F7.0,'' * '',E15.8,'' * '', &
& F7.2,'' * '',F7.2)')J,XTRS(J,JLOOPT),XPRS(J,JLOOPT), &
	  XRVRS(J,JLOOPT),XURS(J,JLOOPT),XVRS(J,JLOOPT)
	ELSE
	  WRITE(INUM,'(5X,I4,'' * '',F7.2,'' * '',F7.0,'' * '',E15.8,'' * '',&
  &  F7.2,'' * '',F7.2)')J,XTRS(J,JLOOPT),XPRS(J,JLOOPT), &
	  XRVRS(J,JLOOPT),XURS(J,JLOOPT),XVRS(J,JLOOPT)
	ENDIF
      ENDDO
      WRITE(INUM,'(1X,78(1H*))')
    CASE('RSPL')
      WRITE(INUM,'(''RS  '',''G:'',A16,28X,'' T:'',F8.0,''s'',''   (1-IK)'')')CGROUP, &
   XTIMRS2(JLOOPT,1)
      WRITE(INUM,'(''NBVAL '',I5)')SIZE(XTRS,2)
      WRITE(INUM,'(1X,78(1H*))')
        WRITE(INUM,'(5X,''K'',4X,''*  THT_RS *  PABST  *'',7X,''RVT'',7X,&
 &      ''*    UT   *    VT'')')
      WRITE(INUM,'(1X,78(1H*))')
      DO J=SIZE(XTRS,2),1,-1
	IF(J == SIZE(XTRS,2))THEN
	  WRITE(INUM,'(''(IK) '',I4,'' * '',F7.2,'' * '',F7.0,'' * '',E15.8,'' * '',&
          &  F7.2,'' * '',F7.2)')XTRS(JLOOPT,J),XPRS(JLOOPT,J), &
		    XRVRS(JLOOPT,J),XURS(JLOOPT,J),XVRS(JLOOPT,J)
	ELSE
	  WRITE(INUM,'(5X,I4,'' * '',F7.2,'' * '',F7.0,'' * '',E15.8,'' * '',&
        &  F7.2,'' * '',F7.2)')XTRS(JLOOPT,J),XPRS(JLOOPT,J), &
		    XRVRS(JLOOPT,J),XURS(JLOOPT,J),XVRS(JLOOPT,J)
	ENDIF
      ENDDO
      WRITE(INUM,'(1X,78(1H*))')
  END SELECT
ENDIF
ENDIF
!
! Avril 99
!
!
!*    3.1.1  Data conversion in mb and g/kg
!
DO JJJ=1,INC                                                     ! do 1
  PRES(JJJ)    = PRES(JJJ) * 1.E-2
  PTEMP(JJJ)   = PTEMP(JJJ)-273.16
    IF (OMIXRAT) THEN                                           ! if 1
      PQV(JJJ) = PQV(JJJ) * 1.E3  ! Mixing ratio used
    ELSE                                                        ! else 1
      PQV(JJJ) = PQV(JJJ)-273.16  ! Dew point used
    ENDIF                                                       ! endif 1
ENDDO                                                           ! enddo 1

IF(JLOOP2 == 1)THEN  !00000000000000

!
!*   3.1.2  Draws the temperature of state line
!
IF(LCOLINE)THEN
  ! 45. = 44.061/.95*.97
!Mars 2000
  IF(ILENT == 1)THEN

    IF(LCOLRSONE)THEN
      CALL GSPLCI(NCOLRSONE)
      CALL GSTXCI(NCOLRSONE)
      CALL GSPMCI(NCOLRSONE)
      CALL GSFACI(NCOLRSONE)
    ENDIF

  ELSE 

    IF(LCOLRS1ONE)THEN
      IF(JLOOPT == 1)THEN
        CALL GSPLCI(NCOLRS1ONE1)
        CALL GSTXCI(NCOLRS1ONE1)
        CALL GSPMCI(NCOLRS1ONE1)
        CALL GSFACI(NCOLRS1ONE1)
      ELSEIF(JLOOPT == 2)THEN
        CALL GSPLCI(NCOLRS1ONE2)
        CALL GSTXCI(NCOLRS1ONE2)
        CALL GSPMCI(NCOLRS1ONE2)
        CALL GSFACI(NCOLRS1ONE2)
      ELSEIF(JLOOPT == 3)THEN
        CALL GSPLCI(NCOLRS1ONE3)
        CALL GSTXCI(NCOLRS1ONE3)
        CALL GSPMCI(NCOLRS1ONE3)
        CALL GSFACI(NCOLRS1ONE3)
      ELSEIF(JLOOPT == 4)THEN
        CALL GSPLCI(NCOLRS1ONE4)
        CALL GSTXCI(NCOLRS1ONE4)
        CALL GSPMCI(NCOLRS1ONE4)
        CALL GSFACI(NCOLRS1ONE4)
      ELSEIF(JLOOPT == 5)THEN
        CALL GSPLCI(NCOLRS1ONE5)
        CALL GSTXCI(NCOLRS1ONE5)
        CALL GSPMCI(NCOLRS1ONE5)
        CALL GSFACI(NCOLRS1ONE5)
      ELSE
      ENDIF

    ELSE
!Mars 2000
      IF(JLOOPT == 2)THEN
        CALL GSPLCI(2)
        CALL GSTXCI(2)
        CALL GSPMCI(2)
        CALL GSFACI(2)
      ELSE IF(JLOOPT == 3)THEN
        CALL GSPLCI(7)
        CALL GSTXCI(7)
        CALL GSPMCI(7)
        CALL GSFACI(7)
      ELSE IF(JLOOPT == 4)THEN
        CALL GSPLCI(5)
        CALL GSTXCI(5)
        CALL GSPMCI(5)
        CALL GSFACI(5)
      ELSE IF(JLOOPT == 5)THEN
        CALL GSPLCI(4)
        CALL GSTXCI(4)
        CALL GSPMCI(4)
        CALL GSFACI(4)
      ELSE IF(JLOOPT == 6)THEN
        CALL GSPLCI(6)
        CALL GSTXCI(6)
        CALL GSPMCI(6)
        CALL GSFACI(6)
      ELSE
        CALL GSPLCI(1)
        CALL GSTXCI(1)
        CALL GSPMCI(1)
        CALL GSFACI(1)
      ENDIF
!Mars 2000
    ENDIF
  ENDIF
!Mars 2000
ENDIF

IF(JLOOPT >1)THEN
  SELECT CASE(CTYPE)
    CASE('CART')
      ZX = .05 +(JLOOPT-2)*(.73/6.)
    CASE('RSPL')
      ZX = .05 +(JLOOPT-2)*(.73/3.)
  END SELECT
  ZY = .985
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
  SELECT CASE(CTYPE)
    CASE('CART')
      if(nverbia > 0)then
      PRINT *,'CTIMEC ',CTIMEC(1:LEN_TRIM(CTIMEC)),' JLOOPT ',JLOOPT,ZX,ZY
      endif
      CALL PLCHHQ(ZX,ZY,CTIMEC(1:LEN_TRIM(CTIMEC)),.008,0.,-1.)
    CASE('RSPL')
      if(nverbia > 0)then
      PRINT *,'CTIMECS ',CTIMECS(1:LEN_TRIM(CTIMECS)),' JLOOPT ',JLOOPT,ZX,ZY
      endif
      CALL PLCHHQ(ZX,ZY,CTIMECS(1:LEN_TRIM(CTIMECS)),.008,0.,-1.)
  END SELECT
  CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
! Mars 2000
ELSE
! IF(LRS)THEN
    CALL GQTXCI(IERR,ICOLI)
    CALL GSPLCI(1)
    CALL GSTXCI(1)
    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
    CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
    CALL RESOLV_TIT('CTITT1',YTEM80)
    IF(YTEM80 /= ' ' .AND. YTEM80 /= 'DEFAULT')THEN
      ZXPOSTITT1=.005; ZXYPOSTITT1=.98
      IF(XPOSTITT1 /= 0.)THEN
        ZXPOSTITT1=XPOSTITT1
      ENDIF
      IF(XYPOSTITT1 /= 0.)THEN
        ZXYPOSTITT1=XYPOSTITT1
      ENDIF
      IF(XSZTITT1 /= 0.)THEN
	CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80(1:LEN_TRIM(YTEM80)),XSZTITT1,0.,-1.)
!       CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80,XSZTITT1,0.,-1.)
      ELSE
	CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80(1:LEN_TRIM(YTEM80)),.012,0.,-1.)
!       CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80,.012,0.,-1.)
      ENDIF

    ENDIF
    CALL GSPLCI(ICOLI)
    CALL GSTXCI(ICOLI)
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
! ENDIF
! Mars 2000
ENDIF

CALL SETUSV ('LW',2000)    ! Heavy line used for the
!CALL SETUSV ('IN',10000)  ! sounding data 
!

DO J = 1,INC                                                     ! do 1
  IF( PRES(J).LT.100. )EXIT
  ZY = ZFY(PRES(J))
  ZX = ZFX(PTEMP(J),ZY)
  IF(J.EQ.1)CALL FRSTPT(ZX,ZY)
  CALL VECTOR(ZX,ZY)
ENDDO                                                           ! enddo 1

CALL SFLUSH
!print *,' AP CALL SFLUSH'
IF(JLOOPT > 1 .AND. .NOT. LCOLINE)THEN
  CALL GSLWSC(1.)
  CALL GSLN(3)     ! Sets dotted line mode
  CALL VECTOR(ZX,ZY+.5*JLOOPT)
  CALL SFLUSH
  CALL GSLN(1)
  SELECT CASE(CTYPE)
    CASE('CART')
      IF(JLOOPT <10)THEN
        WRITE(YC1,'(I1)')JLOOPT
	IN=1
      ELSE
        WRITE(YC2,'(I2)')JLOOPT
	IN=2
      ENDIF
    CASE('RSPL')
      IF(NNST(JLOOPT) <10)THEN
        WRITE(YC1,'(I1)')NNST(JLOOPT)
	IN=1
      ELSE
        WRITE(YC2,'(I2)')NNST(JLOOPT)
	IN=2
      ENDIF
  END SELECT

  IF(IN == 1)THEN
    CALL PLCHHQ(ZX,ZY+.7*JLOOPT,YC1,.008,0.,0.)
  ELSE
    CALL PLCHHQ(ZX,ZY+.7*JLOOPT,YC2,.008,0.,0.)
  ENDIF

  CALL GSLWSC(2.)

ENDIF
!
!*   3.1.3  Draws dewpoint as function of pressure
!
!CALL GSLN(3)     ! Sets dotted line mode
!
IF(OMXRAT)THEN
!
  DO J = 1,INC                                                    ! do 1
    IF(PTEMP(J).LE.-40.)EXIT
    ZY = ZFY(PRES(J))
      IF (OMIXRAT) THEN                 ! Converts mixing ratio to
        ZDWPT = DEWP( PQV(J),PRES(J) )  ! dewpoint temperature
      ELSE
        ZDWPT = PQV(J)                  ! No conversion necessary here
      END IF 
    ZX = ZFX(ZDWPT,ZY)
!   IF(J.EQ.1)CALL FRSTPT(ZX,ZY)
!   CALL VECTOR(ZX,ZY)
    IF(J == 1)THEN
      INEG=0
      CALL FRSTPT(ZX,ZY)
      IF(PQV(J) <= 0.)INEG=1
      IF(PQV(J) >  0.)CALL VECTOR(ZX,ZY)
    ELSE
      IF(PQV(J) <= 0.)THEN
        INEG=1
        CALL FRSTPT(ZX,ZY)
      ELSE
        SELECT CASE(INEG)
          CASE(0)
            CALL VECTOR(ZX,ZY)
          CASE(1)
            CALL FRSTPT(ZX,ZY)
            CALL VECTOR(ZX,ZY)
            INEG=0
        END SELECT
	IF(MOD(J,4) == 0)THEN
	  CALL GSMK(2)
	  CALL GPM(1,ZX,ZY)
	ENDIF
      END IF
    END IF
  ENDDO                                                           ! enddo 1
!
IF(JLOOPT > 1 .AND. .NOT. LCOLINE)THEN
  CALL GSLWSC(1.)
  CALL GSLN(3)
  CALL VECTOR(ZX+1.5,ZY+.7*JLOOPT)
  CALL SFLUSH
  WRITE(YC1,'(I1)')JLOOPT
  CALL GSLN(1)
  CALL PLCHHQ(ZX,ZY+.5*JLOOPT,YC1,.008,0.,0.)
  CALL GSLWSC(2.)
ENDIF
END IF
!
CALL SFLUSH
!print *,' AP CALL SFLUSH2'
IF(LCOLINE)THEN
  IF(JLOOPT == 2)THEN
  ELSE IF(JLOOPT == 3)THEN
  ELSE IF(JLOOPT == 4)THEN
  ELSE IF(JLOOPT == 5)THEN
  ELSE IF(JLOOPT == 6)THEN
  ELSE
  ENDIF
  CALL GSPLCI(1)
  CALL GSPMCI(1)
  CALL GSTXCI(1)
  CALL GSFACI(1)
ENDIF

CALL GSLN(1)  ! Restores solid line 


ENDIF      !00000000000000
!
 61 CONTINUE
!
IF(LRS1 .AND. JLOOP2 == 1 .AND. JLOOPT >1)THEN
  GDOUV=.FALSE.
ELSE
  GDOUV=GDOUVM
ENDIF
!
!
!*     3.2   Plots wind vectors
!
IF(.NOT.GDOUV)GO TO 66
!
!*     3.2.1  Sets arrow scale
!
ZVSCALE=SQRT(PU(1)*PU(1)+PV(1)*PV(1))
! print *,' ZVSCALE ',ZVSCALE
DO JJJ=1,INC                                                    ! do 1
! ZWORKS1(JJJ) = PRES(JJJ)
! ZWORKS2(JJJ) = PTEMP(JJJ)
! ZWORKS3(JJJ) = PQV(JJJ)
! ZWORKS4(JJJ) = PU(JJJ)
! ZWORKS5(JJJ) = PV(JJJ)
  ZVVMAX=SQRT(PU(JJJ)*PU(JJJ)+PV(JJJ)*PV(JJJ))
  IF (ZVVMAX.GT.ZVSCALE) ZVSCALE=ZVVMAX
! print *,' JJJ ZVSCALE ',JJJ,ZVSCALE
!       PRES(JJJ) = PRES(JJJ) * 1.E-2
ENDDO                                                           ! enddo 1
!
if(nverbia >0)then
print *,' AV CALL ECHELLE'
endif
CALL PCSETC('FC',':')
CALL ECHELLE(ILEN,ZHA) ! Sets arrow size
CALL PCSETC('FC','/')
!
if(nverbia >0)then
print *,' AP CALL ECHELLE'
endif
IF(JLOOP2 == 2)THEN
  IF(JLOOPT == 1)THEN
!   print *,' ILENT ',ILENT
    ZINT=(22.5 - (-14.4))/(ILENT-1)
  ENDIF
  ZXM=-14.4+(JLOOPT-1)*ZINT
  SELECT CASE(CTYPE)
    CASE('CART')
      CALL PLCHHQ(ZXM-1.8,43.,CTIMEC(1:LEN_TRIM(CTIMEC)),.009,0.,-1.)
    CASE('RSPL')
      IF(MOD(JLOOPT,2) /= 0)THEN
        CALL PLCHHQ(ZXM-1.8,43.,CTIMECS(1:LEN_TRIM(CTIMECS)),.009,0.,-1.)
      ELSE
        CALL PLCHHQ(ZXM-1.8,42.,CTIMECS(1:LEN_TRIM(CTIMECS)),.009,0.,-1.)
      ENDIF
  END SELECT
ELSE
  ZXM=22.5
  ZINT=1.
ENDIF
if(nverbia >0)then
print *,' ZXM  ZINT ',ZXM,ZINT
endif
CALL LINE(ZXM,0.0,ZXM,44.061)  ! Draws a vertical line for wind display
CALL SFLUSH
!
!!!!!CALL SETUSV('LW',1000)
!
!*    3.2.2  Optional arrow sampling
!
! Only when winds are displayed, computes the distance between
! two adjacent arrows if the arrow number is limited to IMXSMPLUV
!
IF (OSAMPLEUV) THEN                                             ! if 1
  ZDYSMPL=44.061/FLOAT(IMXSMPLUV-1)
ELSE                                                            ! else 1
  ZDYSMPL=0.         
ENDIF                                                           ! endif 1
ZYSMPL=-ZFY(PRES(1))
!
!*    3.3.3  Plots the vectors
!
CALL GSLWSC(2.) ! Sets heavy line
!
#ifdef O2000
CALL VVSETI('CPM',2 )
!CALL VVSETR('AMX',.05 )
!CALL VVSETR('AMN',.005 )
#endif
DO J = 1,INC                                                     ! do 1
!DO J = 1,KNN                                                     ! do 1
  IF( PRES(J).LT.100. )GO TO 66
  ZY1 = ZFY(PRES(J))   ! Locates arrow at the relevant pressure level
  IF(J.GT.1.AND.(OSAMPLEUV.AND.(ZY1-ZYSMPL.LT.ZDYSMPL)))CYCLE
! print *,' ZY1 ',ZY1
! print *,' AVV FLECHE'
  CALL FLECHE(ZXM,ZY1,PU(J),PV(J),ILEN,ZHA)
! print *,' AP FLECHE ZXM,ZY1 ',ZXM,ZY1
  ZYSMPL=ZY1
ENDDO                                                           ! enddo 1
!
 66 CONTINUE
if(nverbia >0)then
print *,' AP 66'
endif
! 
CALL GSLWSC(1.) !Restores initial line width
!
!
!-----------------------------------------------------------------------------
!
!*    4.    NORMAL EXIT 
!           -----------
!
IF (ODOFRAME) CALL FRAME ! FRAME issued if required

  ENDDO          ! Fin DO JLOOPT
  
  IF(LRS1 .AND. JLOOP2 == 1)THEN
    CALL FRAME
    CALL SET(.05,.95,.05,.95,-19.0,27.1,0.0,44.061,1)
    CALL FRSTPT(-19.,0.)
    CALL VECTOR(-19.,44.061)
    CALL VECTOR(27.1,44.061)
    CALL VECTOR(27.1,0.)
    CALL VECTOR(-19.,0.)
    CALL GSCLIP(0)
    CALL PLCHHQ(-19.,-1.,HTEXTE(1:LEN_TRIM(HTEXTE)),.010,0.,-1.)
!!  CALL GSCLIP(1)
    CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!Mars 2000 Altitudes IKB IKE grille de masse
    IF(CTYPE == 'CART')THEN
    ENDIF
!Mars 2000
    CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
    IF(LDATFILE)CALL DATFILE_FORDIACHRO
if(nverbia >0)then
print *,' AP DATFILE'
endif
!Mars 2000
    CALL RESOLV_TIT('CTITT1',YTEM80)
    IF(YTEM80 /= ' ' .AND. YTEM80 /= 'DEFAULT')THEN

      ZXPOSTITT1=.005; ZXYPOSTITT1=.98
      IF(XPOSTITT1 /= 0.)THEN
        ZXPOSTITT1=XPOSTITT1
      ENDIF
      IF(XYPOSTITT1 /= 0.)THEN
        ZXYPOSTITT1=XYPOSTITT1
      ENDIF
      IF(XSZTITT1 /= 0.)THEN
	CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80(1:LEN_TRIM(YTEM80)),XSZTITT1,0.,-1.)
!       CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80,XSZTITT1,0.,-1.)
      ELSE
	CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80(1:LEN_TRIM(YTEM80)),.012,0.,-1.)
!       CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM80,.012,0.,-1.)
      ENDIF

    ENDIF
!Mars 2000
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  ENDIF
!  DO JJJ=1,INC                                                    ! do 1
!    PRES(JJJ)  =  ZWORKS1(JJJ)
!    PTEMP(JJJ) =  ZWORKS2(JJJ)
!    PQV(JJJ)   =  ZWORKS3(JJJ)
!    PU(JJJ)    =  ZWORKS4(JJJ)
!    PV(JJJ)    =  ZWORKS5(JJJ)
!  ENDDO                                                           ! enddo 1
 
ENDDO            ! Fin DO JLOOP2
!
if(nverbia >0)then
print *,' AV RETURN '
endif
!
CALL PCSETC('FC',Y1)
RETURN
!
!-----------------------------------------------------------------------------
!
!*    5.     ARRAY OVERFLOW CONTROL 
!            ----------------------
! Notice: 
! This section has been implemented to conform to
! the former TRACE implentation. It is not called
! in the present TRACE implementation.
!
!*    5.1    Test on T and moisture array sizes
!
      ENTRY TSOUNDTD (PPRES,PPTEMP,PPQV,PPU,PPV,KNN,HEADER, OMIXRAT, ODOFRAME)
!
INC=KNN !00000000 nn <=> nwk 0000000000
!
IF(KNN.GT.JPNWK)THEN
  PRINT *,' Emagram TSOUNDTD: too much data points requested'
  PRINT *,' NN=',KNN,' when maximum allowed is ',JPNWK,', return.'
RETURN
ENDIF
! 
GDOTEMP=.TRUE.
GDOUV=.FALSE.
GO TO 111
!
!*    5.2    Test on wind  array sizes
!
      ENTRY TSOUNDUV (PPRES,PPTEMP,PPQV,PPU,PPV,KNN,HEADER, OMIXRAT, ODOFRAME)
!
INC=KNN  !00000000 nn <=> nwk 0000000000
!
IF(KNN.GT.JPNWK)THEN
  PRINT *,' Emagram TSOUNDUV: too much data points requested'
  PRINT *,' NN=',KNN,' when maximum allowed is ',JPNWK,', return.'
RETURN
ENDIF
! 
GDOTEMP=.FALSE.
GDOUV=.TRUE.
GO TO 111
!
!----------------------------------------------------------------------------
!
!*    6.     EXIT
!            ----
!
END SUBROUTINE TSOUND_FORDIACHRO
