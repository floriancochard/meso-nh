      SUBROUTINE FRAME
      COMMON /GFLASH/MODEF,IOPWKS(100),IOACT(100),NUMOP,IWISSI
C
C  FRAME is designed to effect a break in the picture drawing
C  sequence depending upon whether the workstation type is 
C  MO, or whether it is an OUTPUT or OUTIN workstation.
C  
C  An UPDATE WORKSTATION and CLEAR WORKSTATION is done on all 
C  metafiles and all workstations of type OUTPUT.  For metafiles
C  this inserts an END PICTURE into the metafile.
C  
C  If there are any OUTIN workstations, all of them are updated
C  with an UPDATE WORKSTATION and a pause is done on the OUTIN
C  workstation of most recent creation.  After return from the 
C  pause, a CLEAR WORKSTATION is done on all OUTIN workstations.
C
      INTEGER WKID
      CHARACTER*80 DATREC,STR,ISTR
C
C  First, flush the pen-move buffer.
C
      CALL PLOTIF (0.,0.,2)
C
C  If no workstations are open, return.
C
      CALL GQOPWK (1,IER,NO,ID)
      IF (NO .EQ. 0) RETURN
C
C  Update all workstations.
C
      DO 200 I=1,NO
C
C  Get the workstation ID.
C
        CALL GQOPWK (I,IERR,NO,WKID)
C
C  Get workstation type.
C
        CALL GQWKC (WKID,IER,ICON,ITYPE)
C
C  Get workstation category (0=output; 2=out/in; 4=metafile).
C
        CALL GQWKCA (ITYPE,IER,ICAT)
C
        IF (ICAT .EQ. 4) THEN
C
C  Illegal to call FRAME while a FLASH buffer is open.
C
          IF (MODEF .EQ. 1) THEN
            CALL SETER 
     -    ('FRAME - ILLEGAL TO CALL FRAME WHILE A FLASH BUFFER IS OPEN',      
     -      16,2)
          ENDIF
          CALL GCLRWK(WKID,0)
        ELSE IF (ICAT.EQ.0 .OR. ICAT.EQ.2) THEN
          CALL GUWK(WKID,0)
          IF (ICAT .EQ. 0) THEN
            CALL GCLRWK(WKID,1)
          ENDIF
        ENDIF
  200 CONTINUE
C
C  Pause on the OUTIN workstaton of most recent creation.
C
      DO 100 I=NO,1,-1
        CALL GQOPWK (I,IERR,NO,WKID)
        CALL GQWKC (WKID,IER,ICON,ITYPE)
        CALL GQWKCA (ITYPE,IER,ICAT)
        IF (ICAT.EQ.2) THEN
          ISTR(1:1) = CHAR(0)
          CALL GINST(WKID,1,0,ISTR,1,0.,1279.,0.,1023.,1,1,1,DATREC)       
          CALL GSSTM(WKID,1,0,0)
          CALL GRQST(WKID,1,ISTAT,LOSTR,STR)
          GO TO 110
        ENDIF
  100 CONTINUE
  110 CONTINUE
C
C  Clear all OUTIN worktations.
C
      DO 300 I=1,NO
        CALL GQOPWK (I,IERR,NO,WKID)
        CALL GQWKC (WKID,IER,ICON,ITYPE)
        CALL GQWKCA (ITYPE,IER,ICAT)
        IF (ICAT.EQ.2) THEN
          CALL GCLRWK(WKID,1)
        ENDIF
  300 CONTINUE
      RETURN
      END
C------------------------------------------------------------------------
C
C     ###########################################
      SUBROUTINE CPMPXY(IMAP,XINP,YINP,XOTP,YOTP)
C     ########################################### 
C
C
CC****  *CPMPXY* - Maps compack isocontour points on the Meso-NH coordinate
CC****             sytem verically or horizontally.
CC
CC    PURPOSE
CC    -------
C       Maps compack isocontour points on the Meso-NH coordinate
C    sytem vertically or horizontally. This routine is directly called
C    by the NCAR CPRECT and CPCLDR cotour drawing routines.
C
CC**  METHOD
CC    ------
CC
CC    CPMPXY routine is used within the NCAR Conpack calls to map the contoured
CC   array matrix onto the stretched model cartographic space. 
CC     The plotted data are NOT interpolated onto a regular grid before 
CC   plotting, instead a coordinate stretching technique is used. Basically, 
CC   the contour calculations are made in a "grid index space" where the 
CC   meshsize is uniform and equal to 1 between successive model points (this
CC   corresponds to the x_bar_* and y_bar_* coordinates of the Meso-NH 
CC   technical specification book, page 41). In this "grid index space"
CC   contourlines points are located by two floating-point index coordinates
CC   vaying between 1 and the corresponding array dimension. This "grid index"
CC   coordinates are latter converted back to screen coordinates by CPMPXY to
CC   obtain a correct display.
CC    Using this routine assumes that the NCAR internal "IMAP" parameter
CC   is given the value 4 (arbitrary convention).
CC
CC
CC NOTICE:    CPMPXY and the NCAR graphical utilities are NOT written
CC ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
CC          does not follow the Meso-NH usual rules: it has to be using
CC          a COMMON stack with  static memory allocation of XZZXX and
CC          XZZXY arrays.
CC
CC    EXTERNAL
CC    --------
CC     None
CC
CC    EXPLICIT ARGUMENTS
CC    ------------------
CC
CC       IMAP : Selects the customized mapping, has to be set to 4 (input).
CC       XINP : x-coordinate of the current contour point given as a 
CC              fractionnal grid index (input).
CC       YINP : y-coordinate of the current contour point given as a
CC              fractionnal grid index (input).
CC       XOTP : x-coordinate of the current contour point after re-mapping onto
CC              the true display geometry, given in the NCAR "user coordinate"
CC              system (meters, output)
CC       YOTP : y-coordinate of the current contour point after re-mapping onto
CC              the true display geometry, given in the NCAR "user coordinate"
CC              system (meters, output)
CC
CC       NOTICE: All these dummy arguments are required
CC       ------  by the NCAR CALLS
CC
CC    IMPLICIT ARGUMENTS
CC    ------------------
CC
CC     Common TEMV: Vertical cross-section grid information
CC       ZWORKZ: True altitudes of the current data point iwithin the section
CC               (in meters)
CC       ZZDS  : Abscissa of the section gridpoint along the oblique horizontal
CC               axis of the section (meters)
CC       INX   : Number of datapoint along the section's abscissa
CC       INY   : Number of gridlevel along the section's vertical axis
CC
CC     Common LOGI: Section geometry information flags copied from the 
CC                  fortran-90 MODN_PARA module to be passed to the 
CC                  fortran-77 part of TRACE.
CC       LVERT : copy of LVERTI, .TRUE. if horizontal section activated
CC       LHOR  : copy of LHORIZ, .TRUE. if vertical section activated. 
CC
CC     Common TEMH: Horizontal section grid information
CC       ZZXX  : Meso-NH X coordinate values for the current data points
CC       ZZXY  : Meso-NH Y coordinate values for the current data points
CC       IIMAX : X array dimension
CC       IJMAX : Y array dimension 
CC
CC    REFERENCE
CC    ---------
CC
CC      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
CC       + Book1: Concepts and Fundamentals, to appear in 1994;
CC       + Book2: Technical Reference and Flowcharts, to appear in 1994;
CC       + Book3: Tutorial, November 1994.
CC
CC     NCAR Graphics Technical documentation, UNIX version 3.2,
CC     Scientific computing division, NCAR/UCAR, Boulder, USA.
CC      Volume 1: Fundamentals, Vers. 1, May 1993
CC      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
CC
CC    AUTHOR
CC    ------
CC
CC      J. Duron    * Laboratoire d'Aerologie *
CC
CC    MODIFICATIONS
CC    -------------
CC      Original       01/07/94
CC      Updated   PM   24/01/95
C-------------------------------------------------------------------------------
C
C*     0.   DECLARATIONS
C           ------------
C
C>>>>>>>DRAGOON NOTICE: I ENFORCED "IMPLICIT NONE" IT'S WISE CHECKING...
C
      IMPLICIT NONE
C
C*     0.1  Dummy arguments
C
      INTEGER IMAP
      REAL XINP,YINP
      REAL XOTP,YOTP
C
C*     0.1  Commons 
C
      COMMON/TEMV/ZWORKZ,ZZDS,INX,INY
      COMMON/LOGI/LVERT,LHOR,LPT,LXABS
      COMMON/TEMH/ZZXX,ZZXY,IIMAX,IJMAX
#include "big.h"
C     REAL ZWORKZ(600,300),ZZDS(600),ZZXX(600),ZZXY(300)
c     REAL ZWORKZ(1000,400),ZZDS(1000),ZZXX(1000),ZZXY(400)
      REAL ZWORKZ(N2DVERTX,2500),ZZDS(N2DVERTX)
      REAL ZZXX(N2DVERTX),ZZXY(N2DVERTX)
C     REAL ZWORKZ(200,200),ZZDS(200),ZZXX(200),ZZXY(200)
      LOGICAL  LVERT,LHOR,LPT,LXABS
      INTEGER INX,INY,IIMAX,IJMAX
C
C*    0.2   Local variables
C
c     REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZZZXY
C     DIMENSION ZZZXY(1000,400)
C     DIMENSION ZZZXY(200,200)
C     REAL ZZZXY
      INTEGER LL,JJ,I,J,IX,IY,IXP1,IYP1
      REAL ZDIFX,ZX1,ZX2,ZY,ZDIFY,ZW1,ZW2,ZW3,ZW4,Z1,Z2,ZR
      
C
C------------------------------------------------------------------------------
C
C*    1.   RE-MAPS THE CONTOUR POINTS ONTO THE STRECHED DISPLAY COORDINATES
C          ----------------------------------------------------------------
C
C*    1.1  Stores horizontal section's Y in a 2D workarray
C
c     IF(ALLOCATED(ZZZXY))THEN
c       DEALLOCATE(ZZZXY)
c     ENDIF
c     print *,' MON CPMPXY A MOI',XINP,YINP
C     PRINT *,' In CPMPXY IMAP=',IMAP
      IF(IMAP.EQ.4)THEN
C     PRINT *,' In CPMPXY LMAX=',INX
c     IF(LHOR)THEN
c       ALLOCATE(ZZZXY(1000,400))
c     LL=IIMAX
c     JJ=IJMAX

c     DO 1 I=1,LL
c     DO 2 J=1,JJ
c     DO 1 J=1,JJ
c     DO 2 I=1,LL
c     ZZZXY(I,J)=ZZXY(J)
c2    CONTINUE
c1    CONTINUE
c     ENDIF 
C
C*    1.2  Computes streched X's 
C 
C Nearest gridpoint is located in fractionnal coordinates,
C distance to nearest gridpoint is computed, and converted 
C to Meso NH true location (NCAR user coordinates).
C
      IX=INT(XINP)
C     IF(FLOAT(IX)+.989.LE.XINP)IX=IX+1
      ZDIFX=XINP-FLOAT(IX)
c     print *,' XINP IX ZDIFX LHOR+V ',XINP,IX,ZDIFX,LHOR,LVERT

      IF(LVERT)THEN
      ZX1=ZZDS(MAX(IX,1))
      ZX2=ZZDS(MIN(IX+1,INX))
C     PRINT *,' cpmpxy XINP IX',XINP,IX,' ZX1 2',ZX1,ZX2
      ELSE
      ZX1=ZZXX(MAX(IX,1))
      ZX2=ZZXX(MIN(IX+1,IIMAX))
C     PRINT *,' cpmpxy XINP IX',XINP,IX,' ZX1 2',ZX1,ZX2
      ENDIF
c     PRINT *,' cpmpxy XINP IX',XINP,IX,' ZX1 2',ZX1,ZX2
      XOTP=ZX1+ZDIFX*(ZX2-ZX1)

C
C*    1.3  Computes streched Y's
C
C Same as above, but altitudes are used here, when
C LVERT=.T. Here the four surrounding corners in
C fractional space are located. And a 2D linear
C interpolation is performed to remap onto true
C altitudes and true distances
C
      ZY=YINP
      IY=INT(ZY)
C     IF(FLOAT(IY)+.989.LE.YINP)IY=IY+1
      ZDIFY=ZY-FLOAT(IY)
      
c     print *,' INX,INY ',INX,INY
      IF(LVERT)THEN
c     PRINT *,' cpmpxy YINP IY',YINP,IY
       IXP1=MIN(INX,IX+1)
       IYP1=MIN(INY,IY+1)
       IF(LPT .AND. LXABS)THEN
C Cas LPXT=.T. et LXABSC=.T.
C Cas profil horizontal // X . Permutation volontaire des indices I et J
C car chargement (pour des pbs de place memoire) des temps en I (alors qu'ils
C sont representes en Y) et des valeurs en J alors qu'elles sont representees 
C en abscisse (Chargement dans PVFCT)
C Nota : les X sont eux charges normalement dans ZZDS (de 1 a INX)
C LPT=LPXT
         ZW1=ZWORKZ(IY,IX)
         ZW2=ZWORKZ(IYP1,IX)
         ZW3=ZWORKZ(IY,IXP1)
         ZW4=ZWORKZ(IYP1,IXP1)
       ELSE
         ZW1=ZWORKZ(IX,IY)
         ZW2=ZWORKZ(IX,IYP1)
         ZW3=ZWORKZ(IXP1,IY)
         ZW4=ZWORKZ(IXP1,IYP1)
       ENDIF
       Z1=ZW1+ZDIFY*(ZW2-ZW1)
       Z2=ZW3+ZDIFY*(ZW4-ZW3)
       ZR=Z1+ZDIFX*(Z2-Z1)
      ELSE
       ZW1=ZZXY(MAX(IY,1))
       ZW2=ZZXY(MIN(IY+1,IJMAX))
       ZR=ZW1+ZDIFY*(ZW2-ZW1)
      ENDIF
      YOTP=ZR
c     PRINT *,' xotp,yotp',xotp,yotp
      END IF
      
c     IF(ALLOCATED(ZZZXY))THEN
c       DEALLOCATE(ZZZXY)
c     ENDIF

      RETURN
C
C----------------------------------------------------------------------------
C
C*    2.    EXIT
C           ----
C
      END 
C----------------------------------------------------------------------------
C
C	$Id$
C
C***********************************************************************
C P A C K A G E   E Z M A P   -   I N T R O D U C T I O N
C***********************************************************************
C
C This file contains implementation instructions and the code for the
C package EZMAP.  Banners like the one above delimit the major sections
C of the file.  The code itself is separated into three sections: user-
C level routines, internal routines, and the block data routine which
C determines the default values of internal parameters.  Within each
C section, routines appear in alphabetical order.
C
C***********************************************************************
C P A C K A G E   E Z M A P   -   I M P L E M E N T A T I O N
C***********************************************************************
C
C The EZMAP package is written in FORTRAN-77 and should be relatively
C easy to implement.  The outline data required may be generated by
C running the program
C
C     PROGRAM CONVRT
C       DIMENSION FLIM(4),PNTS(200)
C       REWIND 1
C       REWIND 2
C   1   READ (1,3,END=2) NPTS,IGID,IDLS,IDRS,(FLIM(I),I=1,4)
C       IF (NPTS.GT.1) READ (1,4,END=2) (PNTS(I),I=1,NPTS)
C       WRITE (2) NPTS,IGID,IDLS,IDRS,(FLIM(I),I=1,4),(PNTS(I),I=1,NPTS)
C       GO TO 1
C   2   STOP
C   3   FORMAT (4I4,4F8.3)
C   4   FORMAT (10F8.3)
C     END
C
C with the EZMAP card-image dataset on unit 1.  The output file, on unit
C 2, contains the binary outline data to be used by EZMAP.  The EZMAP
C routine MAPIO (which see) must then be modified to access this file.
C
C***********************************************************************
C T H E   C O D E   -   U S E R - L E V E L   R O U T I N E S
C***********************************************************************
C
      SUBROUTINE MAPDRW
C
C Declare required common blocks.  See MAPBD for descriptions of these
C common blocks and the variables in them.
C
#if defined(NCL511)
      COMMON /MAPCM4/  GRDR,GRID,GRLA,GRLO,GRPO,OTOL,PDRE,PLA1,PLA2,
     +                   PLA3,PLA4,PLB1,PLB2,PLB3,PLB4,PLNO,PLTO,ROTA,
     +                   SRCH,XLOW,XROW,YBOW,YTOW,IDOT,IDSH,IDTL,ILCW,
     +                   ILTS,JPRJ,ELPF,INTF,LBLF,PRMF
      DOUBLE PRECISION GRDR,GRID,GRLA,GRLO,GRPO,OTOL,PDRE,PLA1,PLA2,
     +                   PLA3,PLA4,PLB1,PLB2,PLB3,PLB4,PLNO,PLTO,ROTA,
     +                   SRCH,XLOW,XROW,YBOW,YTOW
      INTEGER          IDOT,IDSH,IDTL,ILCW,ILTS,JPRJ
      LOGICAL          ELPF,INTF,LBLF,PRMF
      SAVE   /MAPCM4/
#else
      COMMON /MAPCM4/ INTF,JPRJ,PHIA,PHIO,ROTA,ILTS,PLA1,PLA2,PLA3,PLA4,
     +                PLB1,PLB2,PLB3,PLB4,PLTR,GRID,IDSH,IDOT,LBLF,PRMF,
     +                ELPF,XLOW,XROW,YBOW,YTOW,IDTL,GRDR,SRCH,ILCW,GRLA,
     +                GRLO,GRPO
      LOGICAL         INTF,LBLF,PRMF,ELPF
      SAVE /MAPCM4/
#endif
      COMMON/EPAISCONT/ZLWCONT
      COMMON/FDC/IFDC

C
C Initialize the package, draw and label the grid, and draw outlines.
C
c     print *,' INTF ',INTF
      IF (INTF) CALL MAPINT
      CALL MAPGRD
      CALL MAPLBL
      CALL GQLWSC(IERR,ZWIDTH)
      CALL GSLWSC(ZLWCONT)
C     CALL GSLWSC(5.)
      IF(IFDC .EQ. 1 .OR. IFDC .EQ. 3)THEN
C     IF(IFDC .NE. 0)THEN
      CALL MPLNDR('Earth..1',3)
c     print *,' MAPDRW AP MPLNDR( IFDC= ',IFDC
      ENDIF
C     CALL MAPLOT
      CALL GSLWSC(ZWIDTH)
C
      RETURN
      END
C     ###############################################
      SUBROUTINE VVUMXY (X,Y,U,V,UVM,XB,YB,XE,YE,IST)
C     ###############################################
C
C
CC****  *VVUMXY* - Maps velocity vectors onto the Meso-NH coordinate system
CC****             for horizontal cross-sections (so far)
CC
CC    PURPOSE
CC    -------
C       Maps velocity vectors onto the Meso-NH coordinate system
C   for horizontal cross-sections. This routine is called directly by 
C   VVINIT and VVECTR NCAR uitilities to draw wind or flux vectors 
C   making allowance for variable mesh sizes. For the time being,
C   only the case of horizontal cross-section is adressed, vertical 
C   cross-sections vectors are not yet implemented. 
C
CC**  METHOD
CC    ------
CC
CC      With the settings used in TRACE (i.e. parameter SET=0, and IMAP=4),
CC   VVUMXY receives arrow locations (X,Y) as grid array indices (values
CC   ranging between 1 and IIMAX or IJMAX), and wind components (U,V) in 
CC   Meso-NH physical units (m/s for winds) from VVINIT or VVECTR. 
CC      First, VVUMXY converts the locations of the vector starting points to 
CC   the Meso-NH  x- and y-  coordinates by using the Meso-NH gridpoint 
CC   locations given in  arrays ZZX and ZZY, and these arrow locations  are 
CC   finally converted to the NCAR normalized device coordinate system by CUFX
CC   or CUFY calls. 
CC      Next, the wind components are converted into arrow lengthes expressed 
CC   in NCAR nomalized device coordinates using the SXDC and SYDC scale
CC   factors (these later being provided automatically by VVINIT). 
CC      Finally VVUMXY returns the vector endpoint coordinates (XE,YE) computed
CC   by adding origin locations and arrow lengthes, both expressed in NCAR 
CC   normalized device coordinates (See NCAR User Guide "Fundamentals", 
CC   Appendix A, p345 section 1).
CC
CC NOTICE:
CC ------
CC
CC   - This calculation assumes that the plotted arrows origins are located on
CC   one of the model grids, and that both wind components  are colocated. The
CC   necessary calculations are done by TRACE. This VVUMXY routine is probably 
CC   not suitable to plot vectors at arbitrary locations between model 
CC   gridpoints.
CC   - Many usefull informations on NCAR vector plots are in form of man pages.
CC   See "man vectors-params" for the description of the tunable parameters
CC   of VVINIT and VVECTR, see "man vvumxy" for the custom mapping of arrows
CC   onto the user coordinate space.
CC   -  Using this routine assumes that the NCAR internal "IMAP" parameter
CC   is given the value 4 (arbitrary convention).
CC   -  VVUMXY and the NCAR graphical utilities are NOT written
CC   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
CC   does not follow the Meso-NH usual rules: it has to be using
CC   COMMON stacks with  static memory allocations.
CC
CC    EXTERNAL
CC    --------
CC
CC     CUFX  : routine to convert a NCAR user coordinate X value into its
CC             NCAR normalized device coordinate equivalent.
CC     CUFY  : routine to convert a NCAR user coordinate Y value into its
CC             NCAR normalized device coordinate equivalent.
CC
CC    EXPLICIT ARGUMENTS
CC    ------------------
CC
CC       X,Y  : (input) position of the vector origin in the grid array index
CC              space (values ranging between 1 and IIMAX or IJMAX, the size
CC              of post-processing section of the Meso-NH arrays),
CC       U,V  : (input) vector components from the U,V arrays for this position
CC       UVM  : (input, not used) magnitude of the U,V components
CC       XB,YB: (output) starting point of the vector in the NCAR normalized 
CC              device coordinate system 
CC       XE,YE: (output) ending point of the vector in the NCAR normalized
CC              device coordinate system
CC       IST  : (output, not used) status results of the mapping: 0 indicates 
CC              success
CC       
CC       NOTICE: All these dummy arguments are required
CC       ------  by the NCAR CALLS
CC
CC    IMPLICIT ARGUMENTS
CC    ------------------
CC     Common VVMAP: Mapping information provided by the NCAR package
CC       IMAP  : Map projection selector, has to be 4 for present TRACE
CC               implementation
CC       SXDC  : X Scale factor to convert physical vector component values to
CC               normalized device coordinate values.
CC       SYDC  : Y Scale factor to convert physical vector component values to
CC               normalized device coordinate values.
CC
CC     Common LOGI: Section geometry information flags copied from the 
CC                  fortran-90 MODN_PARA module to be passed to the 
CC                  fortran-77 part of TRACE (not used so far).
CC       LVERT : copy of LVERTI, .TRUE. if horizontal section activated
CC       LHOR  : copy of LHORIZ, .TRUE. if vertical section activated. 
CC
CC     Common TEMH: Horizontal section grid information
CC       ZZX   : Meso-NH X coordinate values for the current data points
CC       ZZY   : Meso-NH Y coordinate values for the current data points
CC       IIMAX : X array dimension of the postprocessing Meso-NH array section
CC       IJMAX : Y array dimension of the postprocessing Meso-NH array section
CC
CC    REFERENCE
CC    ---------
CC
CC      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
CC       + Book1: Concepts and Fundamentals, to appear in 1994;
CC       + Book2: Technical Reference and Flowcharts, to appear in 1994;
CC       + Book3: Tutorial, November 1994.
CC
CC     NCAR Graphics Technical documentation, UNIX version 3.2,
CC     Scientific computing division, NCAR/UCAR, Boulder, USA.
CC      Volume 1: Fundamentals, Vers. 1, May 1993
CC      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
CC
CC    AUTHOR
CC    ------
CC
CC      J. Duron    * Laboratoire d'Aerologie *
CC
CC    MODIFICATIONS
CC    -------------
CC      Original       01/07/94
CC      Updated   PM   26/01/95
C-------------------------------------------------------------------------------
C
C*     0.   DECLARATIONS
C           ------------
C
C>>>>>>>DRAGOON NOTICE: I ENFORCED "IMPLICIT NONE" IT'S WISE CHECKING...
C
      USE MODD_PVT
C Janvier 2001
      USE MODD_RESOLVCAR
      USE MODN_PARA
C Janvier 2001
      IMPLICIT NONE
C
C*     0.0  Dummy arguments
C
      REAL X, Y,U,V,UVM,XB,YB,XE,YE
      REAL CUFX, CUFY
      INTEGER IST
C
C*     0.1  Commons
C
      COMMON /VVMAP/
     +                IMAP       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                SXDC       ,SYDC       ,NXCT       ,NYCT       ,
     +                RLEN       ,LNLG       ,INVX       ,INVY       ,
     +                ITRT       ,IWCT       ,FW2W       ,FH2H       ,
     +                DVMN       ,DVMX       ,RBIG       ,IBIG
C
      SAVE /VVMAP/
      REAL XVPL,XVPR,YVPB,YVPT,WXMN,WXMX,WYMN,WYMX,XLOV,XHIV,YLOV,YHIV,
     +     SXDC,SYDC,RLEN,FW2W,FH2H,DVMN,DVMX,RBIG
      INTEGER IMAP,NXCT,NYCT,LNLG,INVX,INVY,ITRT,IWCT,IBIG
C
      COMMON/LOGI/LVERT,LHOR,LPT,LXABS
      LOGICAL LVERT,LHOR,LPT,LXABS
C
      COMMON/TEMH/ZZX,ZZY,IIMAX,IJMAX
      COMMON/TEMV/ZWORKZ,ZZDS,INX,INY
#include "big.h"
C     DIMENSION ZZX(200),ZZY(200)
c     DIMENSION ZZX(1000),ZZY(400)
      DIMENSION ZZX(N2DVERTX),ZZY(N2DVERTX)
      REAL ZZX,ZZY
c     REAL ZWORKZ(1000,400),ZZDS(1000)
      REAL ZWORKZ(N2DVERTX,2500),ZZDS(N2DVERTX)
C     REAL ZWORKZ(200,200),ZZDS(200)
      INTEGER IIMAX,IJMAX
      INTEGER INX,INY
      INTEGER ICOLUVG
C Janvier 2001
      INTEGER IER,ICLIP
      REAL ZBID(4)
C Janvier 2001
C
C*     0.2  Local variables
C
      REAL PDTOR,PRTOD,P1XPI,P2XPI,P1D2PI,P5D2PI
C
      INTEGER IX,IY
C
C
C*    0.3   Math constants initialization (not used here)
C
      PARAMETER (PDTOR  = 0.017453292519943,
     +           PRTOD  = 57.2957795130823,
     +           P1XPI  = 3.14159265358979,
     +           P2XPI  = 6.28318530717959,
     +           P1D2PI = 1.57079632679489,
     +           P5D2PI = 7.85398163397448) 
      DATA ICOLUVG/1/
C
C---------------------------------------------------------------------
C
C*    1.    VECTOR ARROW LOCATION AND SCALING
C           ---------------------------------
C
C*    1.1   Converts vector starting point from section array indices
C*          to normalized device coordinates
C
C     print *,' MON VVU....A MOI'
      IF(IMAP.EQ.4)THEN
C	 print *, ' X Y',X,Y,'  SXDC SYDC',SXDC,SYDC
C	 print *, ' X Y',X,Y
C
C NOTICE: It is mandatory to use nearest integer function  NINT here
C
         IX=NINT(X)
         IY=NINT(Y)
C
         IF(LHOR)THEN
	       X=ZZX(IX)
	       Y=ZZY(IY)
         ELSE
C Janvier 2001
           IF(LPV)THEN
	     IF(IX == NPROFILE)THEN
	       X=(ZZDS(1) + ZZDS(NLMAX))/2
	       Y=ZWORKZ(IX,IY)
	     ELSE
	       RETURN
	     ENDIF
	   ELSE
C Janvier 2001
	       X=ZZDS(IX)
	       Y=ZWORKZ(IX,IY)
C Janvier 2001
           ENDIF
	   CALL GQCLIP(IER,ICLIP,ZBID)
	   IF(ICLIP == 0 .AND. (Y > XHMAX .OR. Y < XHMIN))THEN
	     RETURN
	   ENDIF
C Janvier 2001
         ENDIF
C
         XB=CUFX(X)
         YB=CUFY(Y)
C	 PRINT *,' IX IY ',IX,IY,' ZZX(IX)ZZY(IY) ',
C    1         ZZX(IX),ZZY(IY)
 
C        PRINT *,'ZZDS(IX),ZWORKZ(IX,IY) ',ZZDS(IX),ZWORKZ(IX,IY)
C*   1.2   End of vector normalized device coordinate location
C
         XE=XB+U*SXDC
         YE=YB+V*SYDC
C        PRINT *,' XB YB XE YE ',XB,YB,XE,YE
C        PRINT *,' U V SXDC SYDC ',U,V,SXDC,SYDC
      ENDIF
C Essai couleur Mars 2000
      IF(LCOLPVT)THEN
      CALL GSPLCI(NCOL2DUV(IX,IY))
      ELSE
C       IF(NCOLUVG .NE. ICOLUVG)THEN
          CALL GSPLCI(NCOLUVG)
	  ICOLUVG=NCOLUVG
C       ENDIF
      ENDIF
      RETURN
C
C-----------------------------------------------------------------------------
C
C*   2.    EXIT
C          ----
C
      END
C
C	$Id$
C
      SUBROUTINE GERHND(ERRNR,FCTID,ERRFIL)
C
C  ERROR HANDLING
C
      INTEGER ERRNR,FCTID,ERRFIL
C
#if defined(NCL511)
      include 'gkscom-5.1.1.h'
#else
      include 'gkscom.h'
#endif
C
C  Special common blocks containing current error number
C  and file identifier.
C
      COMMON /GKERR1/ ENUM
      COMMON /GKERR2/ FNAME
      INTEGER ENUM
      CHARACTER*6 FNAME
C
C  Record number of error message and maximum number of allowable
C  errors before abort.
C
C  AUGMENTATION VOLONTAIRE DE MAXERR (AVANT = 10)
      DATA MNERR,MAXERR/0,1000/
C
      IF (CUFLAG.EQ.-1 .OR. ERRNR.NE.-109) MNERR = MNERR+1
      IF (MNERR .GT. MAXERR) THEN
        CALL GERLOG(-107,FCTID,ERRFIL)
        STOP
      ENDIF
      ENUM  = ERRNR
      FNAME = GNAM(FCTID+1)
      CALL GERLOG(ERRNR,FCTID,ERRFIL)
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	$Id$
C
C***********************************************************************
C L A B E L B A R   -   I N T R O D U C T I O N
C***********************************************************************
C
C This file contains materials for a package which draws "label bars" -
C horizontal or vertical rectangles divided into boxes (each of which
C is either colored or filled with a pattern), and having labels next
C to it, which serves as a key for a solid-filled plot.
C
C***********************************************************************
C L A B E L B A R   -   I M P L E M E N T A T I O N
C***********************************************************************
C
C LABELBAR is written in standard FORTRAN 77.  No special effort should
C be required to implement it.  It does require various other parts of
C the NCAR Graphics package to have been implemented; in particular, it
C requires the package SOFTFILL, the support routine SETER, and various
C routines from SPPS.
C
C***********************************************************************
C L A B E L B A R   -   U S E R - L E V E L   R O U T I N E S
C***********************************************************************
C
      SUBROUTINE LBLBAR_FORDIACHRO(IHOV,XLEB,XREB,YBEB,YTEB,NBOX,
     +                   WSFB,HSFB,LFIN,
     +                   IFTP,LLBS,NLBS,LBAB)
C
        DIMENSION LFIN(*)
        CHARACTER*(*) LLBS(*)
C
C This routine draws a horizontal or vertical label bar to serve as a
C key for a solid-filled plot.
C
C IHOV is 0 if a horizontal label bar is to be drawn, 1 if a vertical
C label bar is to be drawn.
C
C XLEB is a value between 0 and 1, specifying the position of the left
C edge of the bar.
C
C XREB is a value between 0 and 1, specifying the position of the right
C edge of the bar.
C
C YBEB is a value between 0 and 1, specifying the position of the bottom
C edge of the bar.
C
C YTEB is a value between 0 and 1, specifying the position of the top
C edge of the bar.
C
C ABS(NBOX) is the number of boxes into which the bar is to be divided.
C If NBOX is positive, the boxes will be outlined after being filled;
C if NBOX is negative, this will not be done.
C
C WSFB and HSFB are the width and height, respectively, of each little
C solid-filled box, as fractions of the rectangles resulting from the
C division of the bar into ABS(NBOX) pieces.
C
C LFIN is a list of indices, each of which specifies, in some manner,
C how one of the solid-filled boxes is to be filled.  (For example,
C each may be a color index.)
C
C IFTP specifies the type of solid fill to be used.  If IFTP is zero,
C the routine SFSGFA, in the package SOFTFILL, will be called, with
C an index from LFIN as the value of the argument ICI.  (By default,
C this will result in color fill; the value of the SOFTFILL internal
C parameter 'TY' may be changed to select some other kind of fill by
C SFSGFA.)  If IFTP is non-zero, the user-replaceable routine LBFILL
C will be used to fill the boxes; the default version of this routine
C just does color fill.
C
C LLBS is a list of labels for the solid-filled boxes.
C
C NLBS is the number of labels in the list LLBS.  If NLBS is equal to
C ABS(NBOX)-1, then label I applies to the line separating box I from
C box I+1.  If NLBS is equal to NBOX, then label I applies to box I.  If
C NLBS is equal to ABS(NBOX)+1, then labels 1 and NLBS apply to the left
C and right ends (if IHOV is non-zero, the bottom and top ends) of the
C whole color bar; for values of I not equal to 1 or NLBS, label I
C applies to the line separating box I-1 from box I.
C
C LBAB is a flag having the value 0 if the bar is to be unlabelled, 1
C if the labels are to be below a horizontal bar or to the right of a
C vertical bar, 2 if the labels are to be above a horizontal bar or to
C the left of a vertical bar, 3 if the labels are to be on both sides
C of the bar.
C
C
C Declare the common block where internal parameters are stored.
C
        COMMON /LBCOMN/ ICBL,ICFL,ICLB,WOBL,WOFL,WOLB
        SAVE   /LBCOMN/
        COMMON/GENF/NBCU
C
C Declare the block data routine external to force it to load.
C
        EXTERNAL LBBLDA
C
C Define local arrays to hold X and Y coordinates of boxes.
C
        DIMENSION XCRA(5),YCRA(5)
C
C Define local arrays for use as work arrays by the routine SFSGFA.
C
        DIMENSION RWRK(6),IWRK(8)
C
C Save the current SET parameters and arrange for the use of normalized
C device coordinates.
C
        CALL GETSET (XLVP,XRVP,YBVP,YTVP,XLWD,XRWD,YBWD,YTWD,LNLG)
        CALL    SET (  0.,  1.,  0.,  1.,  0.,  1.,  0.,  1.,   1)
C
C Compute the width and height of each section of the bar and the
C coordinates of the edges of the first solid-filled box.
C
        IF (IHOV.EQ.0) THEN
          WSOB=(XREB-XLEB)/REAL(ABS(NBOX))
          WINC=WSOB
          HSOB=YTEB-YBEB
          HINC=0.
          XLB1=XLEB+.5*(1.-WSFB)*WSOB
          XRB1=XLB1+WSFB*WSOB
          IF (LBAB.EQ.1) THEN
            YBB1=YTEB-HSFB*HSOB
            YTB1=YTEB
          ELSE IF (LBAB.EQ.2) THEN
            YBB1=YBEB
            YTB1=YBEB+HSFB*HSOB
          ELSE
            YBB1=YBEB+.5*(1.-HSFB)*HSOB
            YTB1=YTEB-.5*(1.-HSFB)*HSOB
          END IF
        ELSE
          WSOB=XREB-XLEB
          WINC=0.
          HSOB=(YTEB-YBEB)/REAL(ABS(NBOX))
          HINC=HSOB
          IF (LBAB.EQ.1) THEN
            XLB1=XLEB
            XRB1=XLEB+WSFB*WSOB
          ELSE IF (LBAB.EQ.2) THEN
            XLB1=XREB-WSFB*WSOB
            XRB1=XREB
          ELSE
            XLB1=XLEB+.5*(1.-WSFB)*WSOB
            XRB1=XREB-.5*(1.-WSFB)*WSOB
          END IF
          YBB1=YBEB+.5*(1.-HSFB)*HSOB
          YTB1=YBB1+HSFB*HSOB
        END IF
C
C Draw the bar by filling all of the individual boxes.
C
        CALL GQFACI (IERR,ISFC)
        IF (IERR.NE.0) THEN
          CALL SETER ('LBLBAR - ERROR EXIT FROM GQFACI',1,2)
          STOP
        END IF
C
        IF (ICFL.GE.0) THEN
          CALL GQPLCI (IERR,ISPC)
          IF (IERR.NE.0) THEN
            CALL SETER ('LBLBAR - ERROR EXIT FROM GQPLCI',2,2)
            STOP
          END IF
          CALL GSPLCI (ICFL)
        END IF
C
        IF (WOFL.GT.0.) THEN
          CALL GQLWSC (IERR,STLW)
          IF (IERR.NE.0) THEN
            CALL SETER ('LBLBAR - ERROR EXIT FROM GQLWSC',3,2)
            STOP
          END IF
          CALL GSLWSC (WOFL)
        END IF
C
        DO 101 I=1,ABS(NBOX)
          XCRA(1)=XLB1+REAL(I-1)*WINC
          YCRA(1)=YBB1+REAL(I-1)*HINC
          XCRA(2)=XRB1+REAL(I-1)*WINC
          YCRA(2)=YCRA(1)
          XCRA(3)=XCRA(2)
          YCRA(3)=YTB1+REAL(I-1)*HINC
          XCRA(4)=XCRA(1)
          YCRA(4)=YCRA(3)
          XCRA(5)=XCRA(1)
          YCRA(5)=YCRA(1)
          IF (IFTP.EQ.0) THEN
            CALL SFSGFA (XCRA,YCRA,4,RWRK,6,IWRK,8,LFIN(I))
          ELSE
            CALL LBFILL (IFTP,XCRA,YCRA,5,LFIN(I))
          END IF
  101   CONTINUE
C
        CALL GSFACI (ISFC)
        IF (ICFL.GE.0) CALL GSPLCI (ISPC)
        IF (WOFL.GT.0.) CALL GSLWSC (STLW)
C
C If it is to be done, outline the boxes now.
C
        IF (NBOX.GT.0) THEN
C
          IF (ICBL.GE.0) THEN
            CALL GQPLCI (IERR,ISPC)
            IF (IERR.NE.0) THEN
              CALL SETER ('LBLBAR - ERROR EXIT FROM GQPLCI',4,2)
              STOP
            END IF
            CALL GSPLCI (ICBL)
          END IF
C
          IF (WOBL.GT.0.) THEN
            CALL GQLWSC (IERR,STLW)
            IF (IERR.NE.0) THEN
              CALL SETER ('LBLBAR - ERROR EXIT FROM GQLWSC',5,2)
              STOP
            END IF
            CALL GSLWSC (WOBL)
          END IF
C
          DO 102 I=1,ABS(NBOX)
            XCRA(1)=XLB1+REAL(I-1)*WINC
            YCRA(1)=YBB1+REAL(I-1)*HINC
            XCRA(2)=XRB1+REAL(I-1)*WINC
            YCRA(2)=YCRA(1)
            XCRA(3)=XCRA(2)
            YCRA(3)=YTB1+REAL(I-1)*HINC
            XCRA(4)=XCRA(1)
            YCRA(4)=YCRA(3)
            XCRA(5)=XCRA(1)
            YCRA(5)=YCRA(1)
            IF (IHOV.EQ.0) THEN
              IF (I.EQ.1.OR.WSFB.NE.1.) THEN
                CALL GPL (5,XCRA,YCRA)
              ELSE
                CALL GPL (4,XCRA,YCRA)
              END IF
            ELSE
              IF (I.EQ.1.OR.HSFB.NE.1.) THEN
                CALL GPL (5,XCRA,YCRA)
              ELSE
                CALL GPL (4,XCRA(2),YCRA(2))
              END IF
            END IF
  102     CONTINUE
C
          IF (ICBL.GE.0) CALL GSPLCI (ISPC)
          IF (WOBL.GT.0.) CALL GSLWSC (STLW)

        END IF
C
C If labelling is to be done at all ...
C
        IF (LBAB.NE.0) THEN
C
C ... save the current setting of the PLOTCHAR "text extent" parameter
C and reset it to force computation of "text extent" quantities.
C
          CALL PCGETI ('TE - TEXT EXTENT FLAG',ITEX)
          CALL PCSETI ('TE - TEXT EXTENT FLAG',1)
C
C Find the dimensions of the largest label in the list of labels.
C
          WMAX=0.
          HMAX=0.
C
          DO 104 I=1,NLBS
            NCLB=LEN(LLBS(I))
  103       IF (LLBS(I)(NCLB:NCLB).EQ.' ') THEN
              NCLB=NCLB-1
              IF (NCLB.NE.0) GO TO 103
            END IF
            IF (NCLB.NE.0) THEN
              CALL PLCHHQ (.5,.5,LLBS(I)(1:NCLB),.01,360.,0.)
              CALL PCGETR ('DL - DISTANCE TO LEFT EDGE'  ,DSTL)
              CALL PCGETR ('DR - DISTANCE TO RIGHT EDGE' ,DSTR)
              CALL PCGETR ('DB - DISTANCE TO TOP EDGE'   ,DSTB)
              CALL PCGETR ('DT - DISTANCE TO BOTTOM EDGE',DSTT)
              WMAX=MAX(WMAX,DSTL+DSTR+.02)
              HMAX=MAX(HMAX,DSTB+DSTT+.02)
            END IF
  104     CONTINUE
C
C If the maximum height and width are undefined, quit.
C
          IF (WMAX.LE..02.OR.HMAX.LE..02) GO TO 107
C
C Determine the character width to be used and the resulting offset
C distance to the bottom or top of the label.
C
C         print *,' WSOB ',WSOB
        IF(IHOV /= 0 .AND. NBCU <= 7 .AND. WSOB < .06)WSOB=.06
C         print *,' WSOB MODIFIE ',WSOB
          IF (IHOV.EQ.0) THEN
            HOLA=(1.-HSFB)*HSOB
            IF (LBAB.GE.3) HOLA=HOLA/2.
            WCHR=.01*MIN(WSOB/WMAX,HOLA/HMAX)
            DSTB=(DSTB+.01)*(WCHR/.01)
            DSTT=(DSTT+.01)*(WCHR/.01)
          ELSE
            WOLA=(1.-WSFB)*WSOB
            IF (LBAB.GE.3) WOLA=WOLA/2.
            WCHR=.01*MIN(WOLA/WMAX,HSOB/HMAX)
          END IF
C         print *,' WCHR ',WCHR
C
C Draw the labels.
C
          CALL GQPLCI (IERR,ISCL)
          IF (IERR.NE.0) THEN
            CALL SETER ('LBLBAR - ERROR EXIT FROM GQPLCI',6,2)
            STOP
          END IF
          CALL GQTXCI (IERR,ISCT)
          IF (IERR.NE.0) THEN
            CALL SETER ('LBLBAR - ERROR EXIT FROM GQTXCI',7,2)
            STOP
          END IF
          IF (ICLB.LT.0) THEN
            CALL GSPLCI (ISCT)
          ELSE
            CALL GSPLCI (ICLB)
            CALL GSTXCI (ICLB)
          END IF
          IF (WOLB.GT.0.) THEN
            CALL GQLWSC (IERR,STLW)
            IF (IERR.NE.0) THEN
              CALL SETER ('LBLBAR - ERROR EXIT FROM GQLWSC',8,2)
              STOP
            END IF
            CALL GSLWSC (WOLB)
          END IF
C
          IF (NLBS.LT.ABS(NBOX)) THEN
            XLB1=XLB1+WINC
            YBB1=YBB1+HINC
C           print *,'1 XLB1,YBB1 ',XLB1,YBB1
          ELSE IF (NLBS.EQ.ABS(NBOX)) THEN
            XLB1=XLB1+WSFB*WINC/2.
            YBB1=YBB1+HSFB*HINC/2.
C           print *,'2 XLB1,YBB1 ',XLB1,YBB1
          END IF
C
          DO 106 I=1,NLBS
            NCLB=LEN(LLBS(I))
  105       IF (LLBS(I)(NCLB:NCLB).EQ.' ') THEN
              NCLB=NCLB-1
              IF (NCLB.NE.0) GO TO 105
            END IF
            IF (NCLB.NE.0) THEN
              IF (IHOV.EQ.0) THEN
                IF (LBAB.EQ.1.OR.LBAB.GE.3)
     +            CALL PLCHHQ (XLB1+REAL(I-1)*WSOB,YBB1-DSTT,
     +                            LLBS(I)(1:NCLB),WCHR,0.,0.)
                IF (LBAB.EQ.2.OR.LBAB.GE.3)
     +            CALL PLCHHQ (XLB1+REAL(I-1)*WSOB,YTB1+DSTB,
     +                            LLBS(I)(1:NCLB),WCHR,0.,0.)
              ELSE
C IHOV /= 0 Barre verticale ; LBAB=1 Valeurs a dte ; LBAB=2 Valeurs a g
C JDJDJD
C               IF (LBAB.EQ.1.OR.LBAB.GE.3)
                IF (LBAB.EQ.1)
     +            CALL PLCHHQ (XRB1,YBB1+REAL(I-1)*HSOB,
     +                            LLBS(I)(1:NCLB),WCHR,0.,-1.)
                IF (LBAB.GE.3)
     +            CALL PLCHHQ (XRB1+WCHR,YBB1+REAL(I-1)*HSOB,
     +                            LLBS(I)(1:NCLB),WCHR,0.,-1.)
C JDJDJD
C               IF (LBAB.EQ.2.OR.LBAB.GE.3)
                IF (LBAB.EQ.2)
     +            CALL PLCHHQ (XLB1,YBB1+REAL(I-1)*HSOB,
     +                            LLBS(I)(1:NCLB),WCHR,0.,+1.)
                IF (LBAB.GE.3)
     +            CALL PLCHHQ (XLB1-WCHR,YBB1+REAL(I-1)*HSOB,
     +                            LLBS(I)(1:NCLB),WCHR,0.,+1.)
              END IF
            END IF
  106     CONTINUE
C
          CALL GSPLCI (ISCL)
          IF (ICLB.GE.0) CALL GSTXCI (ISCT)
          IF (WOLB.GT.0.) CALL GSLWSC (STLW)
C
C Restore the original setting of the PLOTCHAR text extent flag.
C
  107     CALL PCSETI ('TE - TEXT EXTENT FLAG',ITEX)
C
        END IF
C
C Restore the original SET parameters.
C
        CALL SET (XLVP,XRVP,YBVP,YTVP,XLWD,XRWD,YBWD,YTWD,LNLG)
C
C Done.
C
        RETURN
C
      END
C     #################################################
      SUBROUTINE SFILL(XWRK,YWRK,NWRK,IAREA,IGRP,NGRPS)
C     #################################################
C
C
CC****  *SFILL* - Performs hatching of plot areas were the
CC                true altitude is lower than the topograpy
CC
CC    PURPOSE
CC    -------
C       When contour plot is drawn, all the locations where the displayed
C     points are below the model topography have to be hatched. SFILL
C     detects these points and perform the hatching.
C
CC**  METHOD
CC    ------
CC
CC      In IMAGE, IMAGEv or IMCOU.., as the contour plots are prepared, the
CC    altitude of the displayed section points are checked to locate points
CC    lower than the local topography. When such points are found they are
CC    marked with a specific "area number" used by SFILL as a mask to 
CC    decide where hatching has to be performed. See the NCAR manual to 
CC    understand how "area numbers" work, this topic is slightly 
CC    involved.. (NCAR contouring tutorial, Vol. 2, pages 12-19, page 120, 
CC    and pages 130-133). 
CC
CC      To summarize, all the lines composing a plot are grouped by "edge 
CC    groups" which may be individually accessed using "group numbers" to
CC    perform specific tasks. For the present purpose only the lines drawn
CC    by CONPACK are important, and they belong to group number 3.
CC      When the contours are computed, CONPACK  assigns "area numbers" to the 
CC    different sub-regions of the plot: typically screen points out of the 
CC    model domain are given a negative area number,  areas between 
CC    isocontours receive area numbers greater than 2, with increasing area 
CC    numbers from the lower contour to the higher one, and TRACE gives an 
CC    area number of 2 to regions under the topography. 
CC      The hatching is therefore performed by scanning the group and area 
CC    numbers to locate the screen points to be hatched, as follows:
CC    - SFILL is called by CONPACK for each contour polygon, with XWRK-YWRK
CC    containing the NWRK points of the current contour, and IAREA-IGRP
CC    containing the corresponding group and area numbers;
CC    - First, the group number is checked to select CONPACK items only,
CC    - Second, the area number is checked to select underground areas,
CC    - If so, the hatching parameters are set (SP=.008, and AN=45 for
CC    slanting hatching) and the SFNORM pattern filling routine is called
CC    to fill the current contour (XWRK-YWRK) with the prescribed pattern.
CC
CC NOTICE:    SFILL and the NCAR graphical utilities are NOT written
CC ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
CC          does not follow the Meso-NH usual rules: it has to be directly 
CC          called by the NCAR CONPACK utility.
CC
CC    EXTERNAL
CC    --------
CC     None
CC
CC    EXPLICIT ARGUMENTS
CC    ------------------
CC
CC       XWRK : x-coordinates (in NCAR fractional system) of the successive 
CC              points forming a given contour enclosing a polygonal area. 
CC       YWRK : y-coordinates (in NCAR fractional system) of the successive
CC              points forming a given contour enclosing a polygonal area.
CC       NWRK : Number of points in XWRK-YWRK to build the contour.
CC       IAREA: Area identifiers for the polygon defined by the XWRK-YWRK and
CC              for each of the NGRPS groups of edges in this plot.
CC       IGRP : Group identifiers for the polygon defined by the XWRK-YWRK and
CC              for each of the NGRPS groups of edges in this plot.
CC       NGRPS: Maximum number of edge groups defined in this plot.
CC
CC       NOTICE: All these dummy arguments are required
CC       ------  by the NCAR CALLS
CC
CC    IMPLICIT ARGUMENTS
CC    ------------------
CC       None
CC
CC    REFERENCE
CC    ---------
CC
CC      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
CC       + Book1: Concepts and Fundamentals, to appear in 1994;
CC       + Book2: Technical Reference and Flowcharts, to appear in 1994;
CC       + Book3: Tutorial, November 1994.
CC
CC     NCAR Graphics Technical documentation, UNIX version 3.2,
CC     Scientific computing division, NCAR/UCAR, Boulder, USA.
CC      Volume 1: Fundamentals, Vers. 1, May 1993
CC      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
CC
CC    AUTHOR
CC    ------
CC
CC      J. Duron    * Laboratoire d'Aerologie *
CC
CC    MODIFICATIONS
CC    -------------
CC      Original       01/07/94
CC      Updated   PM   24/01/95
C-------------------------------------------------------------------------------
C
C*     0.   DECLARATIONS
C
C>>>>>>>DRAGOON NOTICE: I ENFORCED "IMPLICIT NONE" IT'S WISE CHECKING...
C
      IMPLICIT NONE
C
C*     0.1  Dummy arguments
C
      REAL XWRK(*), YWRK(*)
      INTEGER IAREA(*), IGRP(*)
      INTEGER NGRPS,NWRK
C
C*     0.2  Local variables
C
      REAL RSCR(10000)
      INTEGER ISCR(10000)
      INTEGER IA,I,J
C
C------------------------------------------------------------------------------
C
C*     1.    UNDERGROUND AREAS HATCHING
C            --------------------------
C
C*     1.1   Locates CONPACK contour edge lines (group number=3)
C
      DO I=1,NGRPS
C     print *,' IGRP IAREA',IGRP(I),IAREA(I),'  I',I
      IF(IGRP(I).EQ.3)IA=IAREA(I)
      ENDDO
C
C*     1.2   Locates areas with number=2 (underground) and hatches
C
      IF(IA.eq.2)THEN
C       print *,'NWRK ',NWRK,' XWRK YWRK '
	DO J=1,NWRK
C	PRINT *,XWRK(J),YWRK(J)
	ENDDO
      CALL SFSETR('SP',.008)
      CALL SFSETI('AN',45)
      CALL SFSETI('DO',0)
      CALL SFSETI('CH',0)
      CALL  GSMKSC(1.)
      CALL SFNORM(XWRK,YWRK,NWRK,RSCR,10000,ISCR,10000)
      ENDIF
C
C-----------------------------------------------------------------------------
C
C*     2.    EXIT
C            ----
C
      RETURN
      END
C
C     #################################################
      SUBROUTINE SFILLH(XWRK,YWRK,NWRK,IAREA,IGRP,NGRPS)
C     #################################################
C
C
CC****  *SFILLH* - Performs hatching of plot areas were the
CC                true altitude is lower than the topograpy
CC
CC    PURPOSE
CC    -------
C       When contour plot is drawn, all the locations where the displayed
C     points are below the model topography have to be hatched. SFILLH
C     detects these points and perform the hatching.
C
CC**  METHOD
CC    ------
CC
CC      In IMAGE, IMAGEv or IMCOU.., as the contour plots are prepared, the
CC    altitude of the displayed section points are checked to locate points
CC    lower than the local topography. When such points are found they are
CC    marked with a specific "area number" used by SFILLH as a mask to 
CC    decide where hatching has to be performed. See the NCAR manual to 
CC    understand how "area numbers" work, this topic is slightly 
CC    involved.. (NCAR contouring tutorial, Vol. 2, pages 12-19, page 120, 
CC    and pages 130-133). 
CC
CC      To summarize, all the lines composing a plot are grouped by "edge 
CC    groups" which may be individually accessed using "group numbers" to
CC    perform specific tasks. For the present purpose only the lines drawn
CC    by CONPACK are important, and they belong to group number 3.
CC      When the contours are computed, CONPACK  assigns "area numbers" to the 
CC    different sub-regions of the plot: typically screen points out of the 
CC    model domain are given a negative area number,  areas between 
CC    isocontours receive area numbers greater than 2, with increasing area 
CC    numbers from the lower contour to the higher one, and TRACE gives an 
CC    area number of 2 to regions under the topography. 
CC      The hatching is therefore performed by scanning the group and area 
CC    numbers to locate the screen points to be hatched, as follows:
CC    - SFILLH is called by CONPACK for each contour polygon, with XWRK-YWRK
CC    containing the NWRK points of the current contour, and IAREA-IGRP
CC    containing the corresponding group and area numbers;
CC    - First, the group number is checked to select CONPACK items only,
CC    - Second, the area number is checked to select underground areas,
CC    - If so, the hatching parameters are set (SP=.008, and AN=45 for
CC    slanting hatching) and the SFNORM pattern filling routine is called
CC    to fill the current contour (XWRK-YWRK) with the prescribed pattern.
CC
CC NOTICE:    SFILLH and the NCAR graphical utilities are NOT written
CC ------   in Fortran 90, but in Fortran 77.. This sub-section of TRACE
CC          does not follow the Meso-NH usual rules: it has to be directly 
CC          called by the NCAR CONPACK utility.
CC
CC    EXTERNAL
CC    --------
CC     None
CC
CC    EXPLICIT ARGUMENTS
CC    ------------------
CC
CC       XWRK : x-coordinates (in NCAR fractional system) of the successive 
CC              points forming a given contour enclosing a polygonal area. 
CC       YWRK : y-coordinates (in NCAR fractional system) of the successive
CC              points forming a given contour enclosing a polygonal area.
CC       NWRK : Number of points in XWRK-YWRK to build the contour.
CC       IAREA: Area identifiers for the polygon defined by the XWRK-YWRK and
CC              for each of the NGRPS groups of edges in this plot.
CC       IGRP : Group identifiers for the polygon defined by the XWRK-YWRK and
CC              for each of the NGRPS groups of edges in this plot.
CC       NGRPS: Maximum number of edge groups defined in this plot.
CC
CC       NOTICE: All these dummy arguments are required
CC       ------  by the NCAR CALLS
CC
CC    IMPLICIT ARGUMENTS
CC    ------------------
CC       None
CC
CC    REFERENCE
CC    ---------
CC
CC      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
CC       + Book1: Concepts and Fundamentals, to appear in 1994;
CC       + Book2: Technical Reference and Flowcharts, to appear in 1994;
CC       + Book3: Tutorial, November 1994.
CC
CC     NCAR Graphics Technical documentation, UNIX version 3.2,
CC     Scientific computing division, NCAR/UCAR, Boulder, USA.
CC      Volume 1: Fundamentals, Vers. 1, May 1993
CC      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
CC
CC    AUTHOR
CC    ------
CC
CC      J. Duron    * Laboratoire d'Aerologie *
CC
CC    MODIFICATIONS
CC    -------------
CC      Original       01/07/94
CC      Updated   PM   24/01/95
C-------------------------------------------------------------------------------
C
C*     0.   DECLARATIONS
C
C>>>>>>>DRAGOON NOTICE: I ENFORCED "IMPLICIT NONE" IT'S WISE CHECKING...
C
      IMPLICIT NONE
C
C*     0.1  Dummy arguments
C
      REAL XWRK(*), YWRK(*)
      INTEGER IAREA(*), IGRP(*)
      INTEGER NGRPS,NWRK
C
C*     0.2  Commons
C
C
      COMMON/HACHAREA/IHACH(300)
      INTEGER IHACH
C*     0.3  Local variables
C
      REAL RSCR(50000)
      INTEGER ISCR(50000)
      INTEGER IA,I,J,N
C
      REAL ZSP(66)
      INTEGER IND(66),IDO(66),ICH(66),IANG(66)
      INTEGER INDM

C
C------------------------------------------------------------------------------
      DATA ZSP/2*0.,.02,.01,.005,.0025,5*.009,5*.0045,.009,.0045,
     1              .02,.01,.005,.0025,5*.009,5*.0045,.009,.0045,
     2              .00045,.002,.003,.004,.005,.006,.007,.008,.009,
     3              .01,.011,.012,.013,.014,.015,.016,
     4              .001,.002,.003,.004,.005,.006,.007,.008,
     5              .001,.002,.003,.004,.005,.006,.007,.008/
C    5              6*.006/
      DATA IDO/2*0,4*1,11*0,17*1,16*0,16*1/
C     DATA IDO/2*0,4*1,11*0,17*1,16*0,14*1/
      DATA ICH/66*0/
C     DATA ICH/58*0,-1,-2,-3,-4,-5,-1/
      DATA IANG/6*0,45,0,90,-45,-90,135,0,90,-45,-90,135,45,16*0,
     10,135,2*0,45,0,90,-45,-90,45,0,90,-45,-90,135,45,8*135,8*135/
C    14*0,45,0,90,-45,-90,45,0,90,-45,-90,135,45,8*135,5*0,135/
      N=66
      DO I=1,N
	IND(I)=I-1
      ENDDO
C       print *,'NWRK ',NWRK,' XWRK YWRK '
C
C*     1.    UNDERGROUND AREAS HATCHING
C            --------------------------
C
C*     1.1   Locates CONPACK contour edge lines (group number=3)
C
      IA=-5
      DO I=1,NGRPS
C     print *,' IGRP IAREA',IGRP(I),IAREA(I),'  I',I
      IF(IGRP(I).EQ.3)THEN
        IF(IAREA(I) .GT.0)THEN
        IA=IHACH(IAREA(I))
C       print *,' IGRP IAREA',IGRP(I),IAREA(I),'  I',I
C       print *,' IA ',IA
        ENDIF
      ENDIF
      ENDDO
C
C*     1.2   Hatches
C
      IF(IA.GT.0)THEN

C       print *,'NWRK ',NWRK,' XWRK YWRK '
	DO J=1,N
	  IF(IA.EQ.IND(J))THEN
	    INDM=J
C           print *,' SFILLH INDM ',INDM
	  ENDIF
	ENDDO
      IF(INDM .EQ. 1)THEN
C       CALL GSFACI(0)
C       CALL GFA(NWRK,XWRK,YWRK)
      ELSE IF(INDM .EQ. 2)THEN
	CALL GSFACI(1)
	CALL GFA(NWRK,XWRK,YWRK)
      ELSE
      CALL SFSETR('SP',ZSP(INDM))
      CALL SFSETI('AN',IABS(IANG(INDM)))
      CALL SFSETI('DO',IDO(INDM))
      CALL SFSETI('CH',ICH(INDM))
      IF(INDM .GE. 59)CALL GSMKSC(2.)
      CALL SFWRLD(XWRK,YWRK,NWRK,RSCR,50000,ISCR,50000)
      IF(IANG(INDM) .LT. 0)THEN
        CALL SFSETI('AN',IABS(IANG(INDM))+90)
        CALL SFNORM(XWRK,YWRK,NWRK,RSCR,50000,ISCR,50000)
      ENDIF
      ENDIF

      ENDIF
      CALL GSMKSC(1.)
C
C-----------------------------------------------------------------------------
C
C*     2.    EXIT
C            ----
C
      RETURN
      END
C
C
C	$Id$
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LBFILL (IFTP,XCRA,YCRA,NCRA,INDX)
        DIMENSION XCRA(*),YCRA(*)
	INTEGER ISCR(1000)
	REAL    RSCR(1000)
	REAL ZSP(66)
	INTEGER IDO(66),ICH(66),IANG(66)
C
	DATA ZSP/2*0.,.02,.01,.005,.0025,5*.009,5*.0045,.009,.0045,
     1              .02,.01,.005,.0025,5*.009,5*.0045,.009,.0045,
     2              .00045,.002,.003,.004,.005,.006,.007,.008,.009,
     3              .01,.011,.012,.013,.014,.015,.016,
     4              .001,.002,.003,.004,.005,.006,.007,.008,
     5              .001,.002,.003,.004,.005,.006,.007,.008/
C    5              6*.006/
        DATA IDO/2*0,4*1,11*0,17*1,16*0,16*1/
C       DATA IDO/2*0,4*1,11*0,17*1,16*0,14*1/
        DATA ICH/66*0/
C       DATA ICH/58*0,-1,-2,-3,-4,-5,-1/
        DATA IANG/6*0,45,0,90,-45,-90,135,0,90,-45,-90,135,45,16*0,
     10,135,2*0,45,0,90,-45,-90,45,0,90,-45,-90,135,45,8*135,8*135/
C    14*0,45,0,90,-45,-90,45,0,90,-45,-90,135,45,8*135,5*0,135/

C Couleurs
	IF(IFTP.EQ.1)THEN
        CALL GSFACI (INDX)
        CALL GFA (NCRA-1,XCRA,YCRA)
C Hachures et grises
	ELSE
	  IF(INDX.EQ.0)THEN
          ELSE IF(INDX.EQ.1)THEN
C         IF(INDX.EQ.0 .OR. INDX.EQ.1)THEN
            CALL GSFACI (INDX)
            CALL GFA (NCRA-1,XCRA,YCRA)
	  ELSE
	    INDM=INDX+1
	    CALL SFSETR('SP',ZSP(INDM))
            CALL SFSETI('AN',IABS(IANG(INDM)))
	    CALL SFSETI('DO',IDO(INDM))
	    CALL SFSETI('CH',ICH(INDM))
            IF(INDM .GE. 59)CALL GSMKSC(2.)
	    CALL SFNORM(XCRA,YCRA,NCRA,RSCR,1000,ISCR,1000)
	    IF(IANG(INDM) .LT. 0)THEN
	      CALL SFSETI('AN',IABS(IANG(INDM))+90)
	      CALL SFNORM(XCRA,YCRA,NCRA,RSCR,1000,ISCR,1000)
            ENDIF
	    CALL GSMKSC(1.)
	  ENDIF
	ENDIF
        RETURN
      END
C
C Janvier 2001 . Routine importee du Ncar ds package personnel pour
C modif (-> essai de definir une echelle pour les fleches en supprimant
C l'elimination des fleches > ABS(XVHC) ds le cas ou XVHC est <0
C
C       $Id$
C
      SUBROUTINE VVECTR (U,V,P,IAM,VVUDMV,WRK)
C  Janvier 2001
      USE MODD_RESOLVCAR
C  Janvier 2001

C
C Argument dimensions
C
      DIMENSION U(IUD1,*), V(IVD1,*), P(IPD1,*)
C
      DIMENSION WRK(*),IAM(*)
C
      EXTERNAL VVUDMV
C
C Input parameters
C
C U,V    - 2-d arrays holding the component values of a vector field
C P      - A 2-d array containing a scalar data field. The contents
C          of this array may be used to color the vectors 
C IAM    - Area mask array
C VVUDMV - User modifiable masked vector drawing function
C WRK    - work array (currently unused)
C
C Output parameters:
C
C None
C
C PURPOSE                VVECTR draws a representation of a two-
C                        dimensional velocity field by drawing arrows
C                        from each data location.  The length of the
C                        arrow is proportional to the strength of the
C                        field at that location and the direction of
C                        the arrow indicates the direction of the flow
C                        at that location.
C
C ---------------------------------------------------------------------
C
C NOTE:
C Since implicit typing is used for all real and integer variables
C a consistent length convention has been adopted to help clarify the
C significance of the variables encountered in the code for this 
C utility. All local variable and subroutine parameter identifiers 
C are limited to 1,2,or 3 characters. Four character names identify  
C members of common blocks. Five and 6 character variable names 
C denote PARAMETER constants or subroutine or function names.
C
C Declare the VV common blocks.
C
C IPLVLS - Maximum number of color threshold level values
C IPAGMX - Maximum number of area groups allowed in the area map
C
      PARAMETER (IPLVLS = 256, IPAGMX = 64)
C
C
C Integer and real common block variables
C
C
      COMMON /VVCOM/
     +                IUD1       ,IVD1       ,IPD1       ,IXDM       ,
     +                IYDN       ,VLOM       ,VHIM       ,ISET       ,
     +                VRMG       ,VRLN       ,VFRC       ,IXIN       ,
     +                IYIN       ,ISVF       ,UUSV       ,UVSV       ,
     +                UPSV       ,IMSK       ,ICPM       ,UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                RVMN       ,RVMX       ,RDMN       ,RDMX       ,
     +                ISPC       ,RVMD       ,IPLR       ,IVST       ,
     +                IVPO       ,ILBL       ,IDPF       ,IMSG       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
C Arrow size/shape parameters
C
        COMMON / VVARO /
     +                HDSZ       ,HINF       ,HANG       ,IAST       ,
     +                HSIN       ,HCOS       ,FAMN       ,FAMX       ,
     +                UVMG       ,FAIR       ,FAWR       ,FAWF       ,
     +                FAXR       ,FAXF       ,FAYR       ,FAYF       ,
     +                AROX(8)    ,AROY(8)    ,FXSZ       ,FYSZ       ,
     +                FXRF       ,FXMN       ,FYRF       ,FYMN       ,
     +                FWRF       ,FWMN       ,FIRF       ,FIMN       ,
     +                AXMN       ,AXMX       ,AYMN       ,AYMX       ,
     +     	      IACM       ,IAFO       ,WBAD       ,WBTF       ,
     +                WBCF       ,WBDF       ,WBSC
C
C
C Text related parameters
C
        COMMON /VVTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC  ,
     +                FLBS    ,ILBC

C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=36)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CLBT,CILT
C
C Text string parameters
C
      COMMON /VVCHAR/ CSTR,CMNT,CMXT,CZFT,CLBT,CILT
C
      SAVE /VVCOM/, /VVARO/, /VVTXP/, /VVCHAR/
C
C The mapping common block: made available to user mapping routines
C
      COMMON /VVMAP/
     +                IMAP       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                SXDC       ,SYDC       ,NXCT       ,NYCT       ,
     +                RLEN       ,LNLG       ,INVX       ,INVY       ,
     +                ITRT       ,IWCT       ,FW2W       ,FH2H       ,
     +                DVMN       ,DVMX       ,RBIG       ,IBIG
C
      SAVE /VVMAP/
C
C Math constants
C
      PARAMETER (PDTOR  = 0.017453292519943,
     +           PRTOD  = 57.2957795130823,
     +           P1XPI  = 3.14159265358979,
     +           P2XPI  = 6.28318530717959,
     +           P1D2PI = 1.57079632679489,
     +           P5D2PI = 7.85398163397448) 
C
C --------------------------------------------------------------------
C
C Local variable dimensions
C
      PARAMETER (IPLBSZ=10)
      CHARACTER*(IPLBSZ)LBL
      REAL IAR(4)
C
C Local variables
C
C
C The following status and count variables are used to gather
C statistics that are not currently available to the user
C
C IST - Status flag returned from the mapping routine
C ISC - Count of vectors rejected by the mapping routine
C ICT - Count of vector actually plotted
C MXO - Count of vectors rejected because magnitude > maximum
C MNO - Count of vectors rejected because magnitude < minimum
C
C Variables relating to the vector magnitude label
C
C LBL - Character string to hold the vector magnitude label
C NC - Number of characters in the vector magnitude label
C IDP - Local decimal flag for the ENCD routine
C ASH - Scale factor for the vector magnitude label
C
C Zero-field processing and label
C
C IZF - Zero field flag, set TRUE if no vectors are plotted
C XF,YF - fractional length of Zero field string
C IB,IE - beginning and end characters of the string
C W,H   - width and height of the string in fractional coordinates
C XW,YW - position of the string in window coordinates
C
C Vector length adjustment
C
C RAT - Temporary ratio variable
C VA  - adjusted length of current vector
C RA  - ratio of adjusted length to current length
C SMN,SMX - saved value of DVMN and DVMX so they can be restored
C
C Other variables
C
C IOC - the old (saved) color
C IOW - the old (saved) linewidth
C IDA - Do area masking flag
C VMN - The minimum vector size actually plotted (in frac coords)
C VMX - The maximum vector size actually plotted (in frac coords)
C I,J - loop indices for traversing the vector arrays
C K   - loop index for traversing the threshold values
C UI,VI - local copies of the current vector values
C XB,XE,YB,XE - the beginning/ending points of the vector in 
C               the fractional system
C X,Y - mapping of the array indices to a coordinate system
C VLN - length of the current vector in fractional coordinates
C XGV,YGV - X and Y grid value, the scaled distance between each
C           array grid point
C VPL,VPR,VPB,VPT,WDL,WDR,WDB,WDT,ILG - Saved SET call values
C IER,ICL,IAR - Clip query values
C 
C ---------------------------------------------------------------------
C
C Check for valid area map and area group overflow if masking is enabled
C
      IF (IMSK.GT.0) THEN
         IF (IAM(7).GT.IPAGMX) THEN
            CSTR(1:29)='VVECTR - TOO MANY AREA GROUPS'
            CALL SETER (CSTR(1:29),1,1)
            RETURN
         END IF
         IF (IAM(7).LE.0) THEN
            CSTR(1:25)='VVECTR - INVALID AREA MAP'
            CALL SETER (CSTR(1:29),2,1)
            RETURN
         END IF
      END IF
C
C Initialize local variables
C
      NC  = 0
      ICT = 0
      IVC = 0
      ISC = 0
      IZC = 0
      ITH = 0
      MXO = 0
      MNO = 0
      IDA = IMSK
      VMN = RBIG
      VMX = 0.0
      IZF = 1
      SMN=DVMN
      SMX=DVMX
C 
C Save the current color and linewidth, then set the vector
C linewidth. Color must be set on a per vector basis within the 
C main loop. Label text color is set here if a single color is
C specified for all labels. 
C
      CALL GQPLCI(IER,IOC)
      CALL GQTXCI(IER,IOT)
      CALL GQFAIS(IER,IOF)
      CALL GQFACI(IER,IOK)
      CALL GQLWSC(IER,ROW)
      CALL GSLWSC(WDLV)
      IF (ILBC .GE. 0) THEN
         CALL GSTXCI(ILBC)
      END IF
      IF (IAST.NE.0) THEN
         CALL GSFAIS(1)
      END IF
C
C If there are no drawable vectors skip the main loop
C
      IF (UVMX .LE. 0.0) THEN
         IZC=NXCT*NYCT
         DVMX=0.0
         DVMN=0.0
         VMN=0.0
         VMX=0.0
         VFR=0.0
         DRL=0.0
         IAV=0
         GOTO 9800
      END IF
C
C Initialize variables (both local and common block values) that 
C control the mapping between vector magnitude and the realized 
C vector length. 
C
      CALL VVILNS(DRL,VFR,IAV)
C
      IF (DVMX .GT. 2.0*(XVPR - XVPL)) THEN
         CSTR(1:36)='VVECTR - VECTOR NDC LENGTH TOO GREAT'
         CALL SETER (CSTR(1:36),3,1)
         RETURN
      END IF
C
C If using filled arrows initialize the fill arrow data
C For wind barbs initialize data, set up for calling NGDOTS, and
C set the fill color the same as the line color
C
      IF (IAST.EQ.1) THEN
         CALL VVINFA
      ELSE IF (IAST.GE.2) THEN
         CALL NGGETI('CT',ICI)
         CALL NGSETI('CT',1)
         CALL GSFACI(IOC)
         CALL VVINWB
      END IF
C
C Set the scaling for the optional vector labels
C
      IDP = IDPF
      IF (UVMN.NE.0.0 .AND. (ABS(UVMN).LT.0.1 .OR. ABS(UVMN).GE.1.E5))
     +    IDP = 1
      IF (UVMX.NE.0.0 .AND. (ABS(UVMX).LT.0.1 .OR. ABS(UVMX).GE.1.E5))
     +    IDP = 1
      ASH = 1.0
      IF (IDP .NE. 0) ASH =
     +     10.**(3-IFIX(ALOG10(AMAX1(ABS(UVMN),ABS(UVMX)))-500.)-500)
C
C If thinning is in effect, set up the thinning arrays
C
      IV=IXDM*IYDN+1
      IF (RVMD.GT.0.0) THEN
         CALL VVTHIN(U,V,P,WRK(1),WRK(IV))
      END IF
C
C Calculate the grid interval represented by adjacent array
C elements along each axis
C
      XGV=(XHIV-XLOV)/REAL(MAX(1,IXDM-1))
      YGV=(YHIV-YLOV)/REAL(MAX(1,IYDN-1))
C
C Draw the vectors. Note the extra processing if there are special 
C values to consider or the independent scalar array is processed.
C
      DO 201 J=1,IYDN,IYIN
         DO 200 I=1,IXDM,IXIN
C
            UI = U(I,J)
            VI = V(I,J)
C
C If thinning remove thinned out vectors
C
            IF (RVMD.GT.0.0) THEN
               CALL VVTHND(I,J,WRK(1),IS)
               IF (IS.EQ.1) GO TO 194
            END IF
C
C Cull out special values
C
            IF (ISVF .GT. 0) THEN
               IF (UI .EQ. UUSV) THEN
                  IF (ISVF .EQ. 1 .OR. ISVF .EQ. 3) GO TO 199
                  IF (VI .EQ. UVSV .AND. ISVF .EQ. 4) GO TO 199
               ELSE IF (VI .EQ. UVSV) THEN
                  IF (ISVF .EQ. 2 .OR. ISVF .EQ. 3) GO TO 199
               END IF
            END IF
C
C Calculate the vector magnitude or if the polar flag is set
C compute the cartesian component values
C
            IF (IPLR .LE. 0) THEN
               UVMG = SQRT(UI*UI+VI*VI)
            ELSE
               UVMG = ABS(UI)
               IF (IPLR .EQ. 1) VI = PDTOR * VI
               UI = UVMG * COS(VI)
               VI = UVMG * SIN(VI)
            END IF
C
C Bypass vectors that fall outside the user-specified range.
C
            IF (UVMG .LT. UVMN) GO TO 196
C
CCCCCCCCCCCCCSuppression pour voir!!!!!!! -> ca marche 
            IF(LVSUPSCA)THEN
C             IF (UVMG .GT. UVMX) GO TO 197
            ELSE
              IF (UVMG .GT. UVMX) GO TO 197
	    ENDIF
CCCCCCCCCCCCCSuppression pour voir!!!!!!!
C
C Eliminate zero vectors unless using wind barbs
C
            IF (UVMG .EQ. 0.0 .AND. IAST .LT. 2) GO TO 198
C
C If using a scalar array, check for special values in the array, 
C then determine the color to use for the vector
C
            IF (ABS(ICTV) .GE. 2) THEN
C
               IF (ISPC .EQ. 0 .AND. P(I,J) .EQ. UPSV) THEN
                  GO TO 199
               ELSE IF (ISPC .GT. 0 .AND. P(I,J) .EQ. UPSV) THEN
                  IF (IAST .EQ. 0) THEN
                     CALL GSPLCI(ISPC)
                  ELSE IF (IAST .EQ. 1) THEN
                     IF (IACM .EQ. -1 .OR. IACM .GE. 1) THEN
                        CALL GSPLCI(ISPC)
                     END IF
                     IF (IACM .EQ. 0 .OR. ABS(IACM) .GE. 2) THEN
                        CALL GSFACI(ISPC)
                     END IF
                  ELSE
                     CALL GSPLCI(ISPC)
                     CALL GSFACI(ISPC)
                  END IF
                  GO TO 129
               END IF
C
               DO 128 K=1,NLVL,1
                  IF (P(I,J).LE.TVLU(K) .OR. K.EQ.NLVL) THEN
                     IF (IAST .EQ. 0) THEN
                        CALL GSPLCI(ICLR(K))
                     ELSE IF (IAST .EQ. 1) THEN
                        IF (IACM .EQ. -1 .OR. IACM .GE. 1) THEN
                           CALL GSPLCI(ICLR(K))
                        END IF
                        IF (IACM .EQ. 0 .OR. ABS(IACM) .GE. 2) THEN
                           CALL GSFACI(ICLR(K))
                        END IF
                     ELSE
                        CALL GSPLCI(ICLR(K))
                        CALL GSFACI(ICLR(K))
                     END IF
                     IF (ILBC .EQ. -1) THEN
                        CALL GSTXCI(ICLR(K))
                     END IF
                     GO TO 129
                  END IF
 128           CONTINUE
C
 129           CONTINUE
C               
            ELSE IF (ICTV .NE. 0) THEN
C
C If coloring based on vector magnitude, figure out the color
C
               DO 130 K=1,NLVL,1
                  IF (UVMG.LE.TVLU(K) .OR. K.EQ.NLVL) THEN
                     IF (IAST .EQ. 0) THEN
                        CALL GSPLCI(ICLR(K))
                     ELSE IF (IAST .EQ. 1) THEN
                        IF (IACM .EQ. -1 .OR. IACM .GE. 1) THEN
                           CALL GSPLCI(ICLR(K))
                        END IF
                        IF (IACM .EQ. 0 .OR. ABS(IACM) .GE. 2) THEN
                           CALL GSFACI(ICLR(K))
                        END IF
                     ELSE
                        CALL GSPLCI(ICLR(K))
                        CALL GSFACI(ICLR(K))
                     END IF
                     IF (ILBC .EQ. -1) THEN
                        CALL GSTXCI(ICLR(K))
                     END IF
                     GO TO 131
                  END IF
 130           CONTINUE
C
 131           CONTINUE
C
            END IF
C
C Map the vector. If the compatiblity flag is set use the 
C compatibility subroutine.
C
            IF (ICPM .GT. 0) THEN
C
               CALL VVFCPM(I,J,UI,VI,UVMG,XB,YB,XE,YE,IST)
               IF (IST .NE. 0 .AND. IST .NE. -999) GO TO 195
C
            ELSE
C
               X=XLOV+REAL(I-1)*XGV
               Y=YLOV+REAL(J-1)*YGV
               CALL HLUVVMPXY(X,Y,UI,VI,UVMG,XB,YB,XE,YE,IST)
               IF (IST .NE. 0 .AND. IST .NE. -999) GO TO 195
C
            END IF
C
            IF (IAST .GE. 2 .AND. IST .EQ. -999) THEN
               VLN = DVMX
            ELSE
               VLN = SQRT((XE-XB)*(XE-XB)+(YE-YB)*(YE-YB))
               IF (VLN .EQ. 0.0) GO TO 198
C
C Adjust the vector length in proportion to the difference between
C the minimum and maximum display vector magnitudes
C
               IF (IAV.NE.0) THEN
                  VA = VFR+(DVMX - VFR)*(UVMG - UVMN) /(UVMX - UVMN)
                  RA = VA / VLN
                  XE = XB + RA *(XE-XB)
                  YE = YB + RA *(YE-YB)
                  VLN = VA
               END IF
            END IF
C
C Track the minimum/maximum displayed values
C
            IF (UVMG .LT. VMN) VMN=UVMG
            IF (UVMG .GT. VMX) VMX=UVMG
C
C Turn zero field flag off; encode the number if a label is to
C be drawn
C
            IZF = 0
            IF (ILBL .NE. 0) CALL ENCD(UVMG,ASH,LBL,NC,IDP)
C
C Draw the vector
C
            IF (IAST .EQ. 0) THEN
               CALL VVDRAW (XB,YB,XE,YE,VLN,LBL,NC,IAM,VVUDMV,IDA)
            ELSE IF (IAST .EQ. 1) THEN
               CALL VVDRFL (XB,YB,XE,YE,VLN,LBL,NC,IAM,VVUDMV,IDA)
            ELSE
               CALL VVDRWB (XB,YB,XE,YE,VLN,LBL,NC,IAM,VVUDMV,IDA)
            END IF
C
C Statistical data:
C
C Vectors plotted
C
            ICT=ICT + 1
            GOTO 200
C
 194        CONTINUE
C
C Vectors culled out by thinning algorithm
C
            ITH=ITH+1
            GO TO 200
C
 195        CONTINUE
C
C Vectors rejected by mapping routine
C
            ISC=ISC+1
            GO TO 200
C
 196        CONTINUE
C
C Vectors under minimum magnitude
C
            MNO=MNO+1
            GO TO 200
C
 197        CONTINUE
C
C Vectors over maximum magnitude
C
            MXO=MXO + 1
            GO TO 200
C
C Zero length vectors cannot be drawn even if UVMN is 0.0, but
C need to be treated as if they were drawn.
C
 198        CONTINUE
C
            IF (UVMG .LT. VMN) VMN=UVMG
            IZC=IZC + 1
            GO TO 200
C
C Special values
C
 199        CONTINUE
            IVC = IVC+1
C
 200     CONTINUE
 201  CONTINUE
C
C End of main loop.
C
 9800 CONTINUE
C
C Plot statistics
C
      IF (IVST .EQ. 1) THEN
         LUN=I1MACH(2)
         WRITE(LUN,*) 'VVECTR Statistics'
         WRITE(LUN,*) '                    Vectors plotted:',ICT
         WRITE(LUN,*) 'Vectors rejected by mapping routine:',ISC
         WRITE(LUN,*) '    Vectors under minimum magnitude:',MNO
         WRITE(LUN,*) '     Vectors over maximum magnitude:',MXO
         WRITE(LUN,*) '          Other zero length vectors:',IZC
         WRITE(LUN,*) '            Rejected special values:',IVC
         IF (RVMD.GT.0) THEN
            WRITE(LUN,*) '     Vectors below minimum distance:',ITH
         END IF
         WRITE(LUN,*) '   Minimum plotted vector magnitude:',VMN
         WRITE(LUN,*) '   Maximum plotted vector magnitude:',VMX
         IF (ABS(ICTV).GE.2) THEN
            WRITE(LUN,*) '               Minimum scalar value:',PMIN
            WRITE(LUN,*) '               Maximum scalar value:',PMAX
         END IF
         WRITE(LUN,*) ' '
      END IF
C
C Reset attributes
C
      CALL GSPLCI(IOC)
      CALL GSLWSC(ROW)
      CALL GSTXCI(IOT)
      CALL GSFACI(IOK)
      CALL GSFAIS(IOF)
C
C Set the read-only min/max vector sizes to reflect the vectors
C actually drawn
C
      IF (IAV.EQ.0) THEN
         RDMN=VMN*SXDC
      ELSE
         RDMN = VFR+(DVMX - VFR)*(VMN - UVMN) /(UVMX - UVMN)
      END IF
      RDMX=VMX*SXDC
      RVMX=VMX
      RVMN=VMN
C
C If vectors were drawn, write out the vector informational text if 
C called for, else conditionally write the zero field text.
C The size printed out depends on whether absolute or relative
C size mode is in effect.
C 
      IF (IZF .EQ. 0) THEN
C
         IF (CMXT(1:1) .NE. ' ') THEN
            IF (VRMG .GT. 0.0) THEN
               CALL VVARTX(CMXT,IMXP,FMXX,FMXY,FMXS,IMXC,VRMG,DRL)
            ELSE IF (VHIM .LT. 0.0) THEN
               CALL VVARTX(CMXT,IMXP,FMXX,FMXY,FMXS,IMXC,UVMX,DVMX)
            ELSE
               CALL VVARTX(CMXT,IMXP,FMXX,FMXY,FMXS,IMXC,VMX,RDMX)
            ENDIF
         END IF
         IF (CMNT(1:1) .NE. ' ') THEN
            IF (VLOM .LT. 0.0) THEN
               CALL VVARTX(CMNT,IMNP,FMNX,FMNY,FMNS,IMNC,UVMN,DVMN)
            ELSE
               CALL VVARTX(CMNT,IMNP,FMNX,FMNY,FMNS,IMNC,VMN,RDMN)
            END IF
         END IF
C
      ELSE
C
         IF (CZFT(1:1) .NE. ' ') THEN
C
C Turn clipping off and SET to an identity transform
C
            CALL GQCLIP(IER,ICL,IAR)
            CALL GSCLIP(0)
            CALL GETSET(VPL,VPR,VPB,VPT,WDL,WDR,WDB,WDT,ILG)
            CALL SET(0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,1)
C
            XF = XVPL + FZFX * FW2W
            YF = YVPB + FZFY * FH2H
            CALL VVTXLN(CZFT,IPCHSZ,IB,IE)
            CALL VVTXIQ(CZFT(IB:IE),FZFS*FW2W,W,H)
            CALL VVTXPO(IZFP,XF,YF,W,H,XW,YW)
            IF (IZFC .GE. 0) THEN
               CALL GSTXCI(IZFC)
               CALL GSPLCI(IZFC)
            ELSE
               CALL  GSPLCI(IOT)
            END IF
C      
            CALL PLCHHQ(XW,YW,CZFT(IB:IE),FZFS*FW2W,0.0,0.0)
C
            CALL GSTXCI(IOT)
            CALL GSPLCI(IOC)
C
C Restore clipping and the set transformation.
C
            CALL NGSETI('CT',ICI)
            CALL GSCLIP(ICL)
            CALL SET(VPL,VPR,VPB,VPT,WDL,WDR,WDB,WDT,ILG)
C
         END IF
C
      END IF
C
C Restore DVMN and DVMX
C
      DVMN=SMN
      DVMX=SMX
C
C Done
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C $Id$
C
      SUBROUTINE EZXY (XDRA,YDRA,NPTS,LABG)
C
      USE MODD_RESOLVCAR
      USE MODD_TYPE_AND_LH
      USE MODN_NCAR

      REAL XDRA(*),YDRA(*)
C
      CHARACTER*(*) LABG
C
C
C The routine EZXY draws one curve through the points (XDRA(I),YDRA(I)),
C for I = 1, 2, ... NPTS.
C
      CALL AGGETI ('SET .',ISET)
      CALL AGGETI ('FRAM.',IFRA)
      if(nverbia > 0)then
      print *,' EZXY ISET IFRA CTYPE LCOLINE ',ISET,IFRA,CTYPE,LCOLINE
      endif
C
      CALL AGEZSU (2,XDRA,YDRA,NPTS,1,NPTS,LABG,IIVX,IIEX,IIVY,IIEY)
      CALL AGBACK
C
      IF(CTYPE == 'SPXY' .AND. LCOLINE)THEN
	CALL GSLWSC(2.)
	IF(LPHALO .OR. LPHAO)THEN
	  CALL GSPLCI(4)
	ELSEIF(NLOOPN == 1)THEN
	  CALL GSPLCI(3)
	ELSEIF(NLOOPN == 2)THEN
	  CALL GSPLCI(2)
	ELSE
	  CALL GSPLCI(1)
	ENDIF
      ENDIF
      IF (ISET.GE.0) CALL AGCURV (XDRA,1,YDRA,1,NPTS,1)
      IF(CTYPE == 'SPXY' .AND. LCOLINE)THEN
	CALL SFLUSH
	print *,' LSPO,LOSPLO,LSPLO,LPHALO,LPHAO ',LSPO,LOSPLO,LSPLO,LPHALO,LPHAO
	CALL GSLWSC(1.)
        CALL GSPLCI(1)
      ENDIF
C
      IF (IFRA.EQ.1) CALL FRAME
C
      RETURN
C
      END
