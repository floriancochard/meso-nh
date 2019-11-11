CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C AVRIL 2002 
C Ces routines ne sont presentes que pour les streamlines pour
C augmenter la dimension d'1 tableau 750 -> 1500
ccccc Intervention perso dans 2 routines des streamlines (Fin du fichier)
C ce parametre existe aussi ds stinit.f ou je suis intervenue
C Intervention totale ds stumxy.f
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       $Id$
C
      BLOCK DATA STDATA
C
C This routine defines the default values of the Streamline parameters.
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C Old STRMLN interface common blocks
C
      COMMON /STR02/  EXT , SIDE , XLT , YBT
C
      COMMON /STR03/  INITA , INITB , AROWL , ITERP , ITERC , IGFLG
     +             ,  IMSG , UVMSG , ICYC , DISPL , DISPC , CSTOP
C
C ---------------------------------------------------------------------
C
C Initialization of STPAR
C
C IUD1 -- 'UD1' -- First dimension of U
C
      DATA     IUD1 / -1 /
C
C IVD1 -- 'VD1' -- First dimension of V
C
      DATA     IVD1 / -1 /
C
C IPD1 -- 'PD1' -- First dimension of P
C
      DATA     IPD1 / -1 /
C
C IXD1 -- 'XD1' -- Array index for start of data, first dimension
C
      DATA     IXD1 / 1 /
C
C IXDM -- 'XDM' -- Array index for end of data, first dimension
C
      DATA     IXDM / -1 /
C
C IYD1 -- 'YD1' -- Array index for start of data, second dimension
C
      DATA     IYD1 / 1 /
C   
C IYDN -- 'YDN' -- Array index for end of data, second dimension
C
      DATA     IYDN / -1 /
C
C IWKD -- 'WKD' -- Dimension of work array
C
      DATA     IWKD / -1 /
C
C IWKU -- 'WKU' -- Amount of work array actually used (read-only)
C
      DATA     IWKU / 0 /
C
C ISET -- 'SET' -- The Set call flag - Old NSET parameter
C
      DATA     ISET / 1 /
C
C IERR -- 'ERR' -- Error code set by STRMLN (read-only)
C                  -101 - Cyclic flag set for non-cyclic data
C
      DATA     IERR / 0 /
C
C
C IXIN -- 'XIN' -- The X Axis grid increment, must be > 0
C IYIN -- 'YIN' -- The Y Axis grid increment, must be > 0
C
      DATA IXIN / 1 /
      DATA IYIN / 1 /
C
C IXM1 -- (IXDM - 1) (not user accessible)
C IXM2 -- (IXDM - 2) (not user accessible)
C IYM1 -- (IYDN - 1) (not user accessible)
C IYM2 -- (IYDN - 2) (not user accessible)
C
C IMSK -- 'MSK' -- Mask streamlines to an area map: <1 -- no mapping,
C                  >=1 - mapping;
C
      DATA IMSK / 0 /
C
C ICPM -- 'CPM' -- the compatibility mode. If >0 the FX,FY,
C                  functions are used. Additionally, when
C                  used in conjunction with the STRMLN routine, 
C                  has a meaningful range from -4 to +4 inclusive,
C                  where various combinations are allowed to use or
C                  ignore 1) the optional input parameters to
C                  VELVCT, 2) the data in STR01,STR02,STR03,STR04
C                  common, 3) FX, etc routines, as follows:
C
C                  -4: no FX, ignore params, ignore old common data
C                  -3: no FX, ignore params, use old common data
C                  -2: no FX, use params, ignore old common data
C                  -1: no FX, use params, use old common data
C                   0: default, same as -4 if STINIT,STREAM called,
C                      same as +1 if STRMLN or EZSTRM called
C                  +1: FX, use params, use old common data
C                  +2: FX, use params, ignore old common data
C                  +3: FX, ignore params, use old common data
C                  +4: FX, ignore params, ignore old common data
C
C                  FX means using FX,FY
C                  When parameters and common block values are
C                  used they override any values set using the
C                  STSETx routines
C
      DATA ICPM / 0 /
C
C NLVL -- 'NLV' -- number of distinct colors to use for the
C                    independent variable mapping -- cannot exceed
C                    IPLVLS -- default: 16
C                    
      DATA  NLVL /  0 /
C
C IPAI -- 'PAI' -- the current level -- must be set before 
C                   modifying an internal level array value
C
      DATA   IPAI /   1     /
C
C ICTV -- 'CTV' -- compute thresholds flag:
C                  0 -- no vector coloring
C                  < 0: color vectors by magnitude
C                  > 0: color vectors by contents of scalar array P
C                  +-1: number of levels and threshold values already
C                       set
C                  >1,<1: use CTV equally spaced levels
C
      DATA  ICTV /   0     /
C
C WDLV -- 'LWD' -- the width of a streamline
C 
      DATA  WDLV /   1.0   /
C
C UVMN -- 'VMN' -- the minimum displayed vector magnitude, read-only
C UVMX -- 'VMX' -- the maximum displayed vector magnitude, read-only
C PMIN -- 'PMN' -- the minimum scalar array value, read-only
C PMAX -- 'PMX' -- the maximum scalar array value, read-only
C
      DATA UVMN / 0.0 /
      DATA UVMX / 0.0 /
      DATA PMIN / 0.0 /
      DATA PMAX / 0.0 /
C
C ITHN -- 'THN' -- streamline thinning flag
C
      DATA ITHN / 0 /
C
C IPLR -- 'PLR' -- Polar coordinates for UV array flag
C
      DATA IPLR / 0 /
C
C ISST -- 'SST' -- Streamline statistics flag
C
      DATA ISST / 0 /
C
C ICLR -- 'CLR' -- the GKS color index value
C
      DATA  ICLR / IPLVLS * 1 /
C
C TVLU -- 'TVL' -- the list of threshold values
C
      DATA  TVLU / IPLVLS * 0.0 /
C
C End of STPAR intialization
C
C --------------------------------------------------------------------
C
C STTRAN initialization 
C
C User coordinate system to viewport, UV array to user coordinates
C
C UVPS -- 'VPS' -- The viewport mode
C
      DATA UVPS / 0.25 /
C
C UVPL -- 'VPL' -- Viewport left
C
      DATA UVPL / 0.05 /
C
C UVPR -- 'VPR' -- Viewport right
C
      DATA UVPR / 0.95 /
C
C UVPB -- 'VPB' -- Viewport bottom
C
      DATA UVPB / 0.05 /
C
C UVPT -- 'VPT' -- Viewport top
C
      DATA UVPT / 0.95 /
C
C UWDL -- 'WDL' -- Window left
C
      DATA UWDL / 0.0 /
C
C UWDR -- 'WDR' -- Window right
C
      DATA UWDR / 0.0 /
C
C UWDB -- 'WDB' -- Window bottom
C
      DATA UWDB / 0.0 /
C
C UWDT -- 'WDT' -- Window top
C
      DATA UWDT / 0.0 /
C
C UXC1 -- 'XC1' -- minimum X coord
C
      DATA UXC1 / 0.0 /
C
C UXCM -- 'XCM' -- maximum Y coord
C
      DATA UXCM / 0.0 /
C
C UYC1 -- 'YC1' -- minimum Y coord
C
      DATA UYC1 / 0.0 /
C
C UYCN -- 'YCN' -- maximum Y coord
C
      DATA UYCN / 0.0 /
C
C End of STTRAN
C ----------------------------------------------------------------------
C
C STSTRM - Parameters affecting the stream processing algorithm
C
C ISGD -- 'SGD' - Stream starting grid increment (INITA)
C
      DATA ISGD / 2 /
C
C IAGD -- 'AGD' - Arrow placement grid increment (INITB)
C
      DATA IAGD / 2 /
C
C RARL -- 'ARL' - Length of one side of arrow as fraction 
C                 of the viewport width (replaces AROWL)
C
      DATA RARL / 0.012 /
C
C ICKP -- 'CKP' - Check progress after this many iterations (ITERP)
C
      DATA ICKP / 35 /
C
C ICKX -- 'CKX' - Check streamline crossover after this many 
C                 iterations (ITERC). (If negative crossover is 
C                 checked at each entrance to a new grid cell)
C
      DATA ICKX / -99 /
C
C ITRP -- 'TRP' - Interpolaton method (IGFLG)
C                 0 - Use 16 point bessel where possible
C                 non 0 - use bi-linear interpolation everywhere
C
      DATA ITRP / 0 /
C
C ICYK -- 'CYK' - Cyclical data flag (ICYC) If non-zero, instructs
C                 the utility to use cyclic interpolation formulas.
C                 If set and data is non-cyclic the error flag is set.
C
      DATA ICYK / 0 /
C
C RVNL -- 'VNL' - Normalization factor for the differential magnitude.
C                 This controls number of steps in compatibility mode
C                 only when the FX,FY mapping routines are used. See 
C                 parameter 'DFM' for step control when STMPXY and
C                 associated routines are used
C
      DATA RVNL / 0.33 /
C
C ISVF -- 'SVF' - Special value flag  (IMSG)
C                 0 - no special values
C                 non 0 - there may be special values, use only
C                         bi-linear interpolation
      DATA ISVF / 0 /
C
C RUSV -- 'USV' -- The U array special value (UVMSG)
C
      DATA RUSV / 1.0E12 /
C
C RVSV -- 'VSV' -- The V array special value (UVMSG)
C
      DATA RVSV / 1.0E12 /
C
C RNDA -- assigned the NDC value of the arrow size.
C
C ISPC -- 'SPC' -- Special color -- 
C                      < 0: no P special value
C                      = 0: don't draw streamline that has a P spec val
C                      > 0: draw P special values using color SPC
C
      DATA ISPC / -1 /
C
C RPSV -- 'PSV' -- The P array special value
C 
      DATA RPSV / 1.0E12 /
C
C RCDS -- 'CDS' - The critical displacement as a multiple of 'DFM'.
C                 Replaces DISPC. If the streamline has not moved
C                 CDS*DFM units in NDC space after ICKP iterations,
C                 the streamline is terminated
C
      DATA RCDS / 2.0 /
C
C RSSP -- 'SSP' - Stream spacing value as a fraction of the viewport
C                 width; replaces CSTOP. Checked when a new grid box is
C                 entered.
C
      DATA RSSP / 0.015 /
C
C RDFM -- 'DFM' - Differential magnitude as a fraction of the viewport
C                 width. Smaller values result in more steps and a more
C                 accurate approximation of the streamline.
C
      DATA RDFM / 0.02 /
C
C RSMD -- 'SMD' - Streamline minimum distance as a fraction of the 
C                 viewport width.
C
      DATA RSMD / 0.0 /
C
C RAMD -- 'AMD' - Arrow minimum distance as a fraction of the 
C                 viewport width.
C
      DATA RAMD / 0.0 /
C
C IGBS -- 'GBS' - Grid based spacing flag
C
      DATA IGBS / 0 /
C
C End of STSTRM
C --------------------------------------------------------------------
C
C STTXP - Text parameters 
C
C ICCM -- internal - maximum length of character strings
C
      DATA ICSZ / IPCHSZ /
C
C FZFS -- 'ZFS' -- size of text for zero field string as FVPW
C FZFX -- 'ZFX' -- X position of zero field string as FVPW
C FZFY -- 'ZFY' -- Y position of zero field string as FVPW
C IZFP -- 'ZFP' -- zero field string position flag
C IZFC -- 'ZFC' -- color of text for zero field label
C 
      DATA FZFS / 0.033 /
      DATA FZFX / 0.5 /
      DATA FZFY / 0.5 /
      DATA IZFP / 0 /
      DATA IZFC / -1 /
C
C ---------------------------------------------------------------------
C
C Beginning of STCHAR initialization
C
      DATA CZFT / 'ZERO FIELD' /
C
C End of STCHAR initialization
C
C
C ---------------------------------------------------------------------
C
C STMAP initialization
C
C IMAP -- 'MAP' -- the mapping transformation to use
C
      DATA  IMAP / 0 /
C
C ITRT -- 'TRT' -- Transform type flag: 
C                      0  - transform position only
C                      1  - transform position and angle
C                     -1  - transform position, angle, and magnitude
C
      DATA ITRT / 1 /
C
C XVPL,XVPT,YVPB,YVPT -- the viewport values (NDC boundaries)
C
C WXMN,WXMX,WYMN,WYMX -- the window minimum and maximum values
C                        (User coordinate space)
C
C XLOV,XHIV,YLOV,YHIV -- the mapped array endpoint values
C                        (Data coordinate space)
C
C XGDS,YGDS -- size in data coordinates of a grid box
C
C NXCT,NYCT -- number of points in X and Y used for the plot
C
C DFMG -- The magnitude of the diffential increment in NDC space
C
C LNLG -- the log scale mapping flag from SET call
C
C INVX,INVY -- inverse flags for the window boundaries
C
C IWCT - unused
C
C FW2W,FH2H -- fraction of viewport to fraction of viewspace
C
C RBIG,IBIG -- maximum expressible real and integer values
C
C ---------------------------------------------------------------------
C
C STRMLN compatibility common blocks
C
C Beginning of STR02 initialization
C
      DATA EXT  / 0.25 /
      DATA SIDE / 0.90  /
      DATA XLT  / 0.05 /
      DATA YBT  / 0.05 /
C
C End of STR02 initialization
C
C Beginning of STR03 initialization
C
      DATA INITA  / 2 /
      DATA INITB  / 2  /
      DATA AROWL  / 0.33 /
      DATA ITERP  / 35 /
      DATA ITERC  / -99 /
      DATA IGFLG  / 0 /
      DATA ICYC   / 0 /
      DATA IMSG   / 0 /
      DATA UVMSG  / 1.E+36 /
      DATA DISPL  / 0.33 /
      DATA DISPC  / 0.67 /
      DATA CSTOP  / 0.50 /
C
C End of STR03 initialization
C
      END
C
C       $Id$
C
      SUBROUTINE STDRAW  (U,V,UX,VY,IAM,STUMSL)
C
C This routine draws the streamlines.
C
      DIMENSION  U(IUD1,*)             ,V(IVD1,*)
      DIMENSION  UX(IXDM,IYDN)         ,VY(IXDM,IYDN)
      DIMENSION  IAM(*)
      EXTERNAL STUMSL
C
C Input parameters:
C
C U,V    - Vector component arrays
C UX,UY  - Work arrays
C IAM    - Mask array
C STUMSL - User-defined masked streamline drawing routine
C
C The work array has been broken up into two arrays for clarity.  The
C top half of WORK (called UX) will have the normalized (and
C possibly transformed) U components and will be used for book
C keeping.  the lower half of the WORK array (called VY) will
C contain the normalized (and possibly transformed) V components.
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Local declarations
C
C Point and list buffers
C
C The XLS and YLS arrays serve as a circular list. they
C are used to prevent lines from crossing one another.
C
      DIMENSION PX(IPNPTS), PY(IPNPTS)
      DIMENSION XLS(IPLSTL), YLS(IPLSTL)
C
C Parameters:
C
C IPZERO, IPONE, IPTWO - the numbers 0,1,2
C PRZERO - the number 0.0
C PTHREE - the number 3.0
C PSMALL - a small floating point number, large enough to be 
C          detectable by any standard processor
C PMXITR - maximum iteration count for figuring when determining
C          the streamline edge
C
      PARAMETER (IPZERO=0, IPONE=1, IPTWO=2, PRZERO=0.0, PTHREE=3.0)
      PARAMETER (PSMALL=0.000001, PMXITR=32)
C
C Local variables
C
C VSM      - A small value in comparison to the normalized vector mag.
C ISK      - Number of bits to skip in bit routines
C IS1      - ISK + 1
C SSP      - Stream spacing value in fractional (ND) coordinates
C CDS      - Critical displacement in fractional (ND) coordinates
C LCT      - Count of streamlines drawn
C ITO      - Total number of points used to draw all the streamlines
C LCU      - Amount of list currently in use
C LCK      - Current list index
C IDR      - drawing direction 0 + direction 1 - direction
C SGN      - multiplier to change sign based on drawing direction
C IPC      - number of points currently in the point buffer
C ICT      - count of iterations in current streamline
C I,J      - Grid indices
C UIJ,VIJ  - individual vector components
C CVF      - component-wise vector normalizing factor
C LST      - flag indicating the last point in a streamline
C IUX      - integer storage for retrieved bits
C ISV, JSV - saved grid indices where stream starts in + direction
C NBX      - count of grid boxes for current streamline
C LBC      - box checking variable
C X, Y     - current X,Y coordinates (grid coordinates
C DU, DV   - Current normalized interpolated vector components
C XDA, YDA - Current position in data coordinates
C XUS, YUS - Current position in user coordinates
C XND, YND - Current position in NDC space
C XNS, YNS - value of XND and YND saved at the start of the streamline 
C                           and after each progress check
C XN1, YN1 - Previous position in NDC space
C TA       - The tangent angle in NDC space
C DUV      - The differential normalized interpolated vector magnitude
C CSA,SNA  - Cosine and sine of the tangent angle
C XN2,YN2  - The previous previous position in NDC space
C TMG      - Temporary magnitude 
C XT,YT    - Temporary x and y values
C XU1,YU1  - Previous X and Y user coordinate values
C NCT      - Iteration count for determining the streamline edge
C LI       - Index into circular crossover checking list
C IZO      - Zero field flag
C
C --------------------------------------------------------------------
C
C Initialize local variables.
C
C Bit manipulation values
C
c     print *,' ++entree STDRAW'
      VSM = R1MACH(3)*VNML
      ISK = I1MACH(5) - 2
      IS1 = ISK + 1
C
C Stream spacing (setting depends on whether grid relative sizing is
C in effect) and critical displacement
C
      IF (IGBS.EQ.0) THEN
         SSP=RSSP*FW2W
      ELSE
         SSP=RSSP*FW2W/REAL(IXDM)
      END IF
      CDS=RCDS*DFMG
      SMD=RSMD*FW2W
      AMD=RAMD*FW2W
C
C Stream and arrow counters
C
      LCT=0
      ITO=0
      IAC=0
C
C Crossover list variables
C
      LCU = 1
      LCK = 1
      XLS(1) = 0.0
      YLS(1) = 0.0
C
C Current streamline variables
C
      IDR = 0
      SGN = 1.0
      IPC = 0
      ICT = 0
      IUC = 0
      JSV = IYD1
C
C
C Compute the X and Y normalized (and possibly transformed)
C displacement components (UX and VY).
C
      IZO = 1
      DO  40 J=IYD1,IYDN
         DO  30 I=IXD1,IXDM
C
            CALL STMPUV(U(I,J),V(I,J),UIJ,VIJ,IST)
            IF (UIJ.NE.0. .OR. VIJ.NE.0.) THEN
               IZO = 0
               CVF = VNML/SQRT(UIJ*UIJ + VIJ*VIJ)
               UIJ = CVF*UIJ
               VIJ = CVF*VIJ
            END IF
C
C Bookkeeping is done in the least significant bits of the UX array.
C When UIJ is exactly zero this can present some problems.
C To get around this problem, set it to a relatively small number.
C
            IF (UIJ.EQ.0.0) UIJ = VSM
C
C Mask out the least significant two bits as flags for each grid box
C A grid box is any region surrounded by four grid points.
C Flag 1 indicates whether any streamline has previously passed
C through this box.
C Flag 2 indicates whether any directional arrow has already
C appeared in this box.
C Judicious use of these flags prevents overcrowding of
C streamlines and directional arrows.
C
            CALL SBYTES(UIJ,IPZERO,ISK,2,0,1)
C
            IF (MOD(I,ISGD).NE.0 .OR. MOD(J,ISGD).NE.0) THEN
               CALL SBYTES(UIJ,IPONE,IS1,1,0,1)
            END IF
            IF (MOD(I,IAGD).NE.0 .OR. MOD(J,IAGD).NE.0) THEN
               CALL SBYTES(UIJ,IPONE,ISK,1,0,1)
            END IF
C
            UX(I,J) = UIJ
            VY(I,J) = VIJ
C
 30      CONTINUE
 40   CONTINUE
C
C If Zero field bail out
C
      IF (IZO .EQ. 1) THEN
         LCT = 0
         ITO = 0
         GO TO 190
      END IF
C
C
C Start a streamline. Experience has shown that a pleasing picture
C will be produced if new streamlines are started only in grid
C boxes that previously have not had other streamlines pass through
C them. As long as a reasonably dense pattern of available boxes
C is initially prescribed, the order of scanning the grid pts. for
C available boxes is immaterial.
C
 50   CONTINUE
C
C First ensure that the point buffer is clear
C
      IF (IPC.GT.1) CALL STLNSG(PX,PY,IPC,IAM,STUMSL)
C
      LST=0
C
C Find an available box for starting a streamline.
C
      IF (IDR .EQ. 0) THEN
C
         LCT=LCT+1
         ITO = ITO+ICT
         ICT = 0
         DO  70 J=JSV,IYM1
            DO  60 I=IXD1,IXM1
               CALL GBYTES(UX(I,J),IUX,ISK,2,0,1)
               IF (IAND(IUX,IPONE) .EQ. IPZERO) GO TO 80
 60         CONTINUE
 70      CONTINUE
C
C Must be no available boxes for starting a streamline.
C This is the final exit from the streamline drawing loop
C
         GO TO 190
C
 80      CONTINUE
C
C Initialize parameters for starting a streamline.
C Turn the box off for starting a streamline.
C If the special value parameter is turned on, check to see if 
C this box has missing data. If so, find a new starting box.
C
         CALL SBYTES(UX(I,J),IPONE,IS1,1,0,1)
         IF (ISVF .NE. 0) THEN
            CALL STSVCK(U,V,I,J,IST)
            IF (IST .NE. 0) GO TO 50
         END IF
C
         ISV = I
         JSV = J
         IDR = 1
         SGN = +1.0
         IUC = 0
         DST = 0.0
C
      ELSE
C
C Come to here to draw in the opposite direction
C
         IDR = 0
         SGN = -1.
         I = ISV
         J = JSV
         DST = 0.0
         ITO = ITO+ICT
      END IF
C
C Initiate the drawing sequence, resetting counters.
C Start all streamlines in the center of a box.
C Find the initial normalized interpolated vector components.
C
      NBX = 0
      IF (IDR.NE.0) LBC = LCK+1
      IF (LBC.GT.IPLSTL) LBC = 1
      X = FLOAT(I)+0.5
      Y = FLOAT(J)+0.5
      CALL  STDUDV(UX,VY,I,J,X,Y,DU,DV)
      XDA=XLOV+(X-1.0)*XGDS
      YDA=YLOV+(Y-1.0)*YGDS
      DU=DU*SGN
      DV=DV*SGN
C
C Get initial point in the various coordinate systems
C and the tangent angle of the stream. If the compatibility flag
C is positive the FX,FY routines must be used.
C
      IF (ICPM.LE.0) THEN
C
         XDA=XLOV+(X-1.0)*XGDS
         YDA=YLOV+(Y-1.0)*YGDS
         CALL HLUSTMPXY(XDA,YDA,XUS,YUS,IST)
         IF (IST .LT. 0) GO TO 50
         XND=CUFX(XUS)
         YND=CUFY(YUS)
         XN1=XND
         YN1=YND
         CALL HLUSTMPTA(XDA,YDA,XUS,YUS,XND,YND,DU,DV,TA,IST)
         IF (IST .LT. 0) GO TO 50
C
      ELSE
C
         XUS=FX(X,Y)
         IF (XUS.LT.WXMN .OR. XUS.GT.WXMX) GO TO 50 
         YUS=FY(X,Y)
         IF (YUS.LT.WYMN .OR. YUS.GT.WYMX) GO TO 50 
         XND=CUFX(XUS)
         YND=CUFY(YUS)
         TA=ATAN2(DV,DU)
C
      END IF
C
      XNS=XND
      YNS=YND
      ICT=1
      IPC=1
      PX(IPC)=XUS
      PY(IPC)=YUS
C      
C Check grid box directional arrow eligibility
C If a minimum arrow distance is set then the first arrow is not drawn
C
      IF (AMD.LE.0.0) THEN
         CALL GBYTES(UX(I,J),IUX,ISK,2,0,1)
C
         IF (IDR.NE.0 .AND. IAND(IUX,IPTWO).EQ.0) THEN
            IAC=IAC+1
            CALL STARDR(XUS,YUS,XND,YND,TA,IAM,STUMSL,IST)
            IF (IST.EQ.0) THEN
               CALL SBYTES(UX(I,J),IPONE,ISK,1,0,1)
            END IF
C
         END IF
      END IF
C
      ADS = 0.0
C
C Loop to this point until streamline ends
C
 110  CONTINUE
C
C Check to see if the streamline has entered a new grid box.
C
      IF (I.EQ.IFIX(X) .AND. J.EQ.IFIX(Y)) THEN
C
C Must be in same box --  Clear the point buffer if required
C
         IF (IPC .EQ. IPNPTS) THEN
c           print *,' IPC IPNPTS ',IPC,IPNPTS
            CALL STLNSG(PX,PY,IPNPTS,IAM,STUMSL)
            PX(1)=PX(IPNPTS)
            PY(1)=PY(IPNPTS)
            IPC=1
         ENDIF
C
C Determine the interpolated normalized vector at this point
C
         CALL STDUDV (UX,VY,I,J,X,Y,DU,DV)
         IF (DU.EQ.0.0 .AND. DV.EQ.0.0) GO TO 50
C
C Processing diverges depending on the compatibility mode
C
         IF (ICPM .LE. 0) THEN
C
C Get the tangent angle of the streamline at the current point
C in NDC space
C
            CALL HLUSTMPTA(XDA,YDA,XUS,YUS,XND,YND,DU,DV,TA,IST)
            IF (IST.NE.0) GO TO 50
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            IF (XUS.LT.WXMN .OR. XUS.GT.WXMX) GO TO 50
            IF (YUS.LT.WYMN .OR. YUS.GT.WYMX) GO TO 50
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C            
         ELSE
C
C A new point is found in grid space, then transformed into
C user and NDC space. There is no transformation of the tangent
C angle.
            X=X+SGN*DU
            Y=Y+SGN*DV
            XUS=FX(X,Y)
            IF (XUS.LT.WXMN .OR. XUS.GT.WXMX) GO TO 50 
            YUS=FY(X,Y)
            IF (YUS.LT.WYMN .OR. YUS.GT.WYMX) GO TO 50 
            XND=CUFX(XUS)
            YND=CUFY(YUS)
            TA=ATAN2(DV,DU)
C
         END IF
C
C Count the point and add it to the point buffer
C
         ICT=ICT+1
         IPC=IPC+1
         PX(IPC)=XUS
         PY(IPC)=YUS
C
         IF (ICPM.LT.1) THEN
C
            IF (LST .EQ. 1) GO TO 50
C
C The increment in NDC space needs to be proportional to the
C magnitude of the interpolated vector, in order to ensure that
C progress checking works at points of convergence or divergence.
C The square enhances the effectiveness of the technique.
C
            DUV=(DU*DU+DV*DV)/(VNML*VNML)
            CSA=COS(TA)*SGN
            SNA=SIN(TA)*SGN
C
C The current point is adjusted one third of the distance back to
C the previous point. Empirically, in most cases, this seems to
C decrease the inaccuracy resulting from the use of a finite valued
C differential step.
C
            XN2=XN1
            YN2=YN1
            XN1=XND+(XN2-XND)/PTHREE
            YN1=YND+(YN2-YND)/PTHREE
            XND=XN1+CSA*DFMG*DUV
            YND=YN1+SNA*DFMG*DUV
            XD = XND - XN1
            YD = YND - YN1
            DST = DST + SQRT(XD*XD+YD*YD)
C
C If the increment takes the line outside the viewport, find an
C interpolated point on the grid edge. Set a flag indicating
C the end of the stream
C
            IF (XND .LT. XVPL) THEN
               XND = XVPL
               IF (ABS(CSA).GT.0.1) THEN
                  TMG = (XND-XN1)/CSA
                  YND = YN1+SNA*TMG
               ENDIF
               LST = 1
            ELSE IF (XND .GT. XVPR) THEN
               XND = XVPR
               IF (ABS(CSA).GT.0.1) THEN
                  TMG = (XND-XN1)/CSA
                  YND = YN1+SNA*TMG
               ENDIF
               LST = 1
            ELSE IF (YND .LT. YVPB) THEN
               YND = YVPB
               IF (ABS(SNA).GT.0.1) THEN
                  TMG = (YND-YN1)/SNA
                  XND = XN1+CSA*TMG
               END IF
               LST = 1
            ELSE IF (YND .GT. YVPT) THEN
               YND = YVPT
               IF (ABS(SNA).GT.0.1) THEN
                  TMG = (YND-YN1)/SNA
                  XND = XN1+CSA*TMG
               END IF
               LST = 1
            END IF
C
C Now that the new point has been found in NDC space, find its
C coordinates in user, data, and grid space.
C
            XU1=XUS
            YU1=YUS
            XUS=CFUX(XND)
            YUS=CFUY(YND)
C
C Even if the point is within NDC and User boundaries it can still be 
C outside the data area. In this case we use an iterative technique to
C determine the end of the streamline.
C
            CALL HLUSTIMXY(XUS,YUS,XDA,YDA,IST)
            IF (IST.GE.0) THEN
               X=(XDA-XLOV)/XGDS+1.0
               Y=(YDA-YLOV)/YGDS+1.0
            ELSE
               NCT=1
C
C Loop to this point dividing the distance in half at each step
C
 120           CONTINUE
               XT=XU1+(XUS-XU1)/2.0
               YT=YU1+(YUS-YU1)/2.0
               IF (NCT.GE.PMXITR) GO TO 50
               IF (ABS(XUS-XU1).LE.PSMALL .AND. 
     +              ABS(YUS-YU1).LE.PSMALL) THEN
                  XUS=XU1
                  YUS=YU1
                  CALL HLUSTIMXY(XUS,YUS,XDA,YDA,IST)
                  IF (IST.LT.0) GO TO 50
               ELSE
                  CALL HLUSTIMXY(XT,YT,XDA,YDA,IST)
                  NCT=NCT+1
                  IF (IST.LT.0) THEN
                     XUS=XT
                     YUS=YT
                  ELSE
                     XU1=XT
                     YU1=YT
                  END IF
                  GO TO 120
               END IF
C
               XND=CUFX(XUS)
               YND=CUFY(YUS)
               LST=1
            END IF
C
C
C If on the top or right edge of the grid space, decrease the X and/or
C Y value by a small amount so the interpolation routine still works.
C
            IF (IFIX(X).GE.IXDM) X=FLOAT(IXDM)-PSMALL
            IF (IFIX(Y).GE.IYDN) Y=FLOAT(IYDN)-PSMALL
C
         END IF
C
C Check streamline progress every 'ICKP' iterations.
C
         IF (MOD(ICT,ICKP).EQ.0) THEN
            IF (ABS(XND-XNS).LT.CDS 
     +           .AND. ABS(YND-YNS).LT.CDS) THEN
               GO TO 50
            END IF
            XNS=XND
            YNS=YND
         END IF
C
C If the circular list does not need to be checked for
C streamline crossover, return to the top of the main loop.
C
         IF (ICKX.LT.0 .OR. MOD(ICT,ICKX).NE.0) GO TO 110
C
      ELSE
C
C Must have entered a new grid box  check for the following :
C (1) Are the new points on the grid?
C (2) Check for missing data if msg data flag (ISVF) has been set.
C (3) Is this box eligible for a directional arrow?
C (4) Location of this entry versus other streamline entries
C
         I = IFIX(X)
         J = IFIX(Y)
         NBX = NBX+1
C
C Check (1) (Only performed in compatibility mode)
C
         IF (ICPM.GT.0) THEN
            IF (I.LT.IXD1 .OR. I.GT.IXM1 
     +           .OR. J.LT.IYD1 .OR. J.GT.IYM1) THEN
               GO TO  50
            END IF
         END IF
C
C Check (2)
C
         IF (ISVF.NE.0) THEN
            CALL STSVCK(U,V,I,J,IST)
            IF (IST .NE. 0) GO TO 50
         END IF
C
C Check (3) -- postpone actually drawing the arrow until after the 
C crossover check, if crossover detected the arrow will not be drawn.
C
         IDA = 0
         CALL GBYTES(UX(I,J),IUX,ISK,2,0,1)
         IF (IAND(IUX,IPTWO) .EQ. 0) THEN
            IF (DST-ADS .GT. AMD) THEN
               ADS = DST
               IDA = 1
            END IF
         END IF
C
      END IF
C
C Check (4) (performed any time streamline crossover is checked)
C
      DO 140 LI=1,LCU
         IF (ABS(XND-XLS(LI)) .LE. SSP .AND.
     +        ABS(YND-YLS(LI)) .LE. SSP) THEN
            IF (LBC.LE.LCK .AND.
     +           (LI.LT.LBC .OR. LI.GT.LCK)) THEN
               GO TO 50
            ELSE IF (LBC.GT.LCK .AND. 
     +              (LI.LT.LBC .AND. LI.GT.LCK)) THEN
               GO TO 50
            END IF
         END IF
 140  CONTINUE
C
      LCU = MIN0(LCU+1,IPLSTL)
      LCK = LCK+1
c     IF (LCK.GT.IPLSTL)print *,'***attention LCK= ',IPLSTL
      IF (LCK.GT.IPLSTL) LCK = 1
      XLS(LCK) = XND
      YLS(LCK) = YND
      CALL SBYTES(UX(I,J),IPONE,IS1,1,0,1)
      IF (NBX.GE.5) THEN
         LBC = LBC+1
         IF (LBC.GT.IPLSTL) LBC = 1
      END IF
C
      IF (IDA.EQ.1) THEN
         CALL STARDR(XUS,YUS,XND,YND,TA,IAM,STUMSL,IST)
         IAC = IAC + 1
         IF (IST .EQ. 0) THEN
            CALL SBYTES(UX(I,J),IPONE,ISK,1,0,1)
         END IF
         IDA = 0
      END IF

C
C Return to top of drawing loop
C
      GO TO 110
C
C
C Final exit
C
  190 CONTINUE
C
      IF (IZO .EQ. 1) THEN
         CALL STZERO
      END IF
C
C Plot statistics
C
      IF (ISST.EQ.1) THEN
         LUN=I1MACH(2)
         WRITE(LUN,*) 'STREAM Statistics'
         WRITE(LUN,*) '                Streamlines plotted:',LCT
         WRITE(LUN,*) '      Total differential step count:',ITO
         WRITE(LUN,*) ' '
      END IF
C
C Set the workspace used parameter
C
      IWKU = 2*IXDM*IYDN
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
      SUBROUTINE STARDR(XUS,YUS,XND,YND,TA,IAM,STUMSL,IST)
C
C This routine draws the arrow. Calculations are in fractional
C coordinates to ensure uniform arrows irrespective of the 
C mapping in effect.
C A small fraction of the differential change is used to find the
C tangent angle at the current position. Once the angle is known the
C arrow can be drawn at a fixed size independent of the mapping
C routine currently employed.
C
C Input parameters:
C
C XUS,YUS - current position in user space
C XND,YND - current position in NDC space
C TA    - Angle in NDC
C IAM   - Area mask array
C STUMSL - User defined masked streamline drawing routine
C
C Output parameters:
C
C IST - Status code, indicates success or failure
C
      DIMENSION  IAM(*)
      EXTERNAL STUMSL
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Point buffers
C
      DIMENSION AX(3), AY(3)
C
C Local variables
C
C AX, AY   - Arrow head point buffers
C DXW, DYW - Change in X,Y in window coordinates
C XF, YF   - Arrow head position in the fractional system
C DXF,DYF  - Incremental change in the fractional system
C PHI      - Tangent angle
C K        - Loop index and sign factor for each edge of the arrow
C KK       - Index for the arrow head array, within the loop
C D30      - Half the angle of the point of the arrow head (about 30 o)
C XX,YY    - Ends of the arrow in window coordinates
C
C Parameters:
C
C PHFANG - Half the angle of the arrow head (0.5 in radians is 
C          approximately equivalent to 30 degrees)
C PLWFCT - Linewidth factor, arrow size is increased by this 
C          much when the linewidth is greater than 1.0

      PARAMETER (PHFANG=0.5, PLWFCT=0.15)
C
C ---------------------------------------------------------------------
C
c     print *,' ++entree STARDR'
      IST=0
C
      AX(2) = XUS
      AY(2) = YUS
      FLW = 1.0 + PLWFCT*MAX(0.0,WDLV-1.0)
C
      DO 10 K = -1,1,2
C
C K serves as a sign determining factor; KK indexes the point array.
C
         KK=K+2
         D30 = -(P1D2PI-TA)+FLOAT(K)*PHFANG
         XX = +RNDA*FLW*SIN(D30)+XND
         YY = -RNDA*FLW*COS(D30)+YND
         AX(KK) = CFUX(XX)
         AY(KK) = CFUY(YY)
C
 10   CONTINUE
C
      CALL STLNSG(AX,AY,3,IAM,STUMSL)
      
C
C Done
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
      SUBROUTINE STLNSG(X,Y,IPC,IAM,STUMSL)
C
C This routine draws a single streamline segment based on the current
C contents of the point buffers. If masking is in effect the area
C line drawing subroutine, ARDRLN is called. Otherwise CURVE is
C invoked. 
C  
C Input parameters:
C
C X,Y - Point arrays
C IPC - Number of points
C IAM   - Area mask array
C STUMSL - User-defined masked streamline drawing routine
C
      DIMENSION X(IPC), Y(IPC)
      DIMENSION  IAM(*)
      EXTERNAL STUMSL
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
      DIMENSION IAI(IPGRCT),IAG(IPGRCT)
      DIMENSION XO(IPNPTS), YO(IPNPTS)
C
C ---------------------------------------------------------------------
C
c     print *,' ++entree STLNSG'
      IF (IMSK.LT.1) THEN
         CALL CURVE(X,Y,IPC)
         CALL SFLUSH
      ELSE
         CALL ARDRLN(IAM, X, Y, IPC, XO, YO, IPC, 
     +        IAI, IAG, IPGRCT, STUMSL)
      END IF
C
C Done
C 
      RETURN
      END
C
C ---------------------------------------------------------------------
C
      SUBROUTINE STSVCK(U,V,I,J,IST)
C
      DIMENSION  U(IUD1,*), V(IVD1,*)
C
C Checks for special values in the vicinity of I,J
C
C Input parameters
C
C U,V - vector field components array
C I,J - current array position
C
C Output parameters:
C
C IST - status value, 0 if no special values in neighborhood
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C ---------------------------------------------------------------------
C
c     print *,' ++entree STSVCK'
      IST = 0
C
      IF (I.EQ.IXDM .OR. J.EQ.IYDN) THEN
         IF (U(I,J).EQ.RUSV) THEN
            IST = -1
         ELSE IF (V(I,J).EQ.RVSV) THEN
            IST = -1
         END IF
         RETURN
      END IF

      IF (U(I,J).EQ.RUSV) THEN
         IST = -1
      ELSE IF (U(I,J+1).EQ.RUSV) THEN
         IST = -1
      ELSE IF (U(I+1,J).EQ.RUSV) THEN
         IST = -1
      ELSE IF (U(I+1,J+1).EQ.RUSV) THEN
         IST = -1
      ELSE IF (V(I,J).EQ.RVSV) THEN
         IST = -1
      ELSE IF (V(I,J+1).EQ.RVSV) THEN
         IST = -1
      ELSE IF (V(I+1,J).EQ.RVSV) THEN
         IST = -1
      ELSE IF (V(I+1,J+1).EQ.RVSV) THEN
         IST = -1
      END IF
C
C Done
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
      SUBROUTINE STMPUV(UI,VI,UO,VO,IST)
C
C Maps the U,V vector component values
C
C Input parameters:
C
C UI,VI  - Input values of U,V
C
C     Output parameters:
C
C UO,VO  - Output mapped component values
C IST    - Status value
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C Statement functions for field tranformations
C
      FU(X,Y) = X
      FV(X,Y) = Y
C
C ---------------------------------------------------------------------
C
c     print *,' ++entree STMPUV'
      IST = 0
C
C Input array polar mode
C
      IF (IPLR .LT. 1) THEN
         UT=UI
         VT=VI
      ELSE IF (IPLR .EQ. 1) THEN
         UT = UI*COS(PDTOR*VI)
         VT = UI*SIN(PDTOR*VI)
      ELSE IF (IPLR .GT. 1) THEN
         UT = UI*COS(VI)
         VT = UI*SIN(VI)
      END IF
C
C Allow mapping using FU,FV functions
C
      UO = FU(UT,VT)
      VO = FV(UT,VT)
C
C Done
C
      RETURN
      END
C
C ---------------------------------------------------------------------
C
      SUBROUTINE STZERO
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
c     print *,' ++entree STZERO'
      IF (CZFT(1:1) .EQ. ' ') THEN
         RETURN
      END IF
C
      CALL GQPLCI(IER,IOC)
      CALL GQTXCI(IER,IOT)
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
C     Restore clipping and the set transformation.
C     
      CALL GSCLIP(ICL)
      CALL SET(VPL,VPR,VPB,VPT,WDL,WDR,WDB,WDT,ILG)
C
C Done
C
      RETURN
      END



C
C       $Id$
C
      SUBROUTINE STDUDV (UX,VY,I,J,X,Y,DU,DV)
C
C Input parameters:
C
C UX,VY  - the arrays containing normalized vector field data
C I,J    - the current grid indices
C X,Y    - the X,Y position relative to the grid
C
C Output parameters:
C
C DU,DV  - Interpolated value of the vector field components
C          at the specified point 
C
C Interpolation routine to calculate the displacemant components.
C The philosphy here is to utilize as many points as possible
C (within reason) in order to obtain a pleasing and accurate plot.
C Interpolation schemes desired by other users may easily be
C substituted if desired.
C
      DIMENSION UX(IXDM,*), VY(IXDM,*)
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C FDLI  - Double linear interpolation formula
C FBESL - Bessel 16 pt interpolation formula ( most used formula )
C FQUAD - Quadratic interpolation formula
C
      FDLI(Z,Z1,Z2,Z3,DX,DY) = (1.-DX)*((1.-DY)*Z +DY*Z1)
     +                         +     DX *((1.-DY)*Z2+DY*Z3)
      FBESL(Z,ZP1,ZP2,ZM1,DZ)=Z+DZ*(ZP1-Z+0.25*(DZ-1.)*((ZP2-ZP1-Z+ZM1)
     +                        +0.666667*(DZ-0.5)*(ZP2-3.*ZP1+3.*Z-ZM1)))
      FQUAD(Z,ZP1,ZM1,DZ)=Z+0.5*DZ*(ZP1-ZM1+DZ*(ZP1-2.*Z+ZM1))
C
C ---------------------------------------------------------------------
C
c     print *,' ++entree STDUDV'
      DX = X-AINT(X)
      DY = Y-AINT(Y)
      ITF = 1
      IM1 = I-1
      IP2 = I+2
C
C Determine which interpolation formula to use 
C depending on I,J location or the special flags
C
      IF (I.GE.IXDM .OR. J.GE.IYDN) THEN
C
C This branch should never be taken if STDRAW is correct, but is 
C included for safety
C
         RETURN
C
      ELSE IF(ISVF.NE.0 .OR. ITRP.NE.0) THEN
         ITF = 1
      ELSE IF (J.GT.IYD1 .AND. J.LT.IYM1 
     +        .AND. I.GT.IXD1 .AND. I.LT.IXM1) THEN
         ITF = 2
      ELSE IF (J.EQ.IYM1 .AND. I.GT.IXD1 .AND. I.LT.IXM1) THEN
         ITF = 3
      ELSE IF (J.EQ.IYD1) THEN
         ITF = 1
      ELSE IF (ICYK.NE.1) THEN
         IF (I.EQ.IXD1) THEN
            ITF = 1
         ELSE IF (I.EQ.IXM1) THEN
            ITF = 4
         END IF
      ELSE IF (I.EQ.IXD1 .AND. J.LT.IYM1) THEN 
         IM1 = IXM1
         ITF = 2
      ELSE IF (I.EQ.IXM1 .AND. J.LT.IYM1) THEN
         IP2 = IXD1+1
         ITF = 2
      ELSE IF (J.EQ.IYM1 .AND. I.EQ.IXD1) THEN
         IM1 = IXM1
         ITF = 3
      ELSE IF (J.EQ.IYM1 .AND. I.EQ.IXM1) THEN
         IP2 = IXD1+1
         ITF = 3
      END IF
C
      IF (ITF .EQ. 1) THEN
C
C Double linear interpolation formula. This scheme works at all points
C but the resulting streamlines are not as pleasing as those drawn
C by FBESL or FQUAD. Currently this is utilized
C only at certain boundary points or if ITRP is not equal to zero,
C or if special value processing is turned on.
C
         DU = FDLI(UX(I,J),UX(I,J+1),UX(I+1,J),UX(I+1,J+1),DX,DY)
         DV = FDLI(VY(I,J),VY(I,J+1),VY(I+1,J),VY(I+1,J+1),DX,DY)
C
      ELSE IF (ITF .EQ. 2) THEN
C
C 16 point bessel interpolation scheme.
C
         UJM1 = FBESL(UX(I,J-1),UX(I+1,J-1),UX(IP2,J-1),UX(IM1,J-1),DX)
         UJ   = FBESL(UX(I,J),UX(I+1,J),UX(IP2,J),UX(IM1,J),DX)
         UJP1 = FBESL(UX(I,J+1),UX(I+1,J+1),UX(IP2,J+1),UX(IM1,J+1),DX)
         UJP2 = FBESL(UX(I,J+2),UX(I+1,J+2),UX(IP2,J+2),UX(IM1,J+2),DX)
         DU   = FBESL(UJ,UJP1,UJP2,UJM1,DY)
         VJM1 = FBESL(VY(I,J-1),VY(I+1,J-1),VY(IP2,J-1),VY(IM1,J-1),DX)
         VJ   = FBESL(VY(I,J),VY(I+1,J),VY(IP2,J),VY(IM1,J),DX)
         VJP1 = FBESL(VY(I,J+1),VY(I+1,J+1),VY(IP2,J+1),VY(IM1,J+1),DX)
         VJP2 = FBESL(VY(I,J+2),VY(I+1,J+2),VY(IP2,J+2),VY(IM1,J+2),DX)
         DV   = FBESL(VJ,VJP1,VJP2,VJM1,DY)
C
      ELSE IF (ITF .EQ. 3) THEN
C
C 12 point interpolation scheme applicable to one row from top boundary
C
         UJM1 = FBESL(UX(I,J-1),UX(I+1,J-1),UX(IP2,J-1),UX(IM1,J-1),DX)
         UJ   = FBESL(UX(I,J),UX(I+1,J),UX(IP2,J),UX(IM1,J),DX)
         UJP1 = FBESL(UX(I,J+1),UX(I+1,J+1),UX(IP2,J+1),UX(IM1,J+1),DX)
         DU   = FQUAD(UJ,UJP1,UJM1,DY)
         VJM1 = FBESL(VY(I,J-1),VY(I+1,J-1),VY(IP2,J-1),VY(IM1,J-1),DX)
         VJ   = FBESL(VY(I,J),VY(I+1,J),VY(IP2,J),VY(IM1,J),DX)
         VJP1 = FBESL(VY(I,J+1),VY(I+1,J+1),VY(IP2,J+1),VY(IM1,J+1),DX)
         DV   = FQUAD(VJ,VJP1,VJM1,DY)
C
      ELSE IF (ITF .EQ. 4) THEN
C
C 9 point interpolation scheme for use in the non-cyclic case
C at I=IXM1; J > IYD1 and J <= IYM1
C
         UJP1 = FQUAD(UX(I,J+1),UX(I+1,J+1),UX(IM1,J+1),DX)
         UJ   = FQUAD(UX(I,J),UX(I+1,J),UX(IM1,J),DX)
         UJM1 = FQUAD(UX(I,J-1),UX(I+1,J-1),UX(IM1,J-1),DX)
         DU   = FQUAD(UJ,UJP1,UJM1,DY)
         VJP1 = FQUAD(VY(I,J+1),VY(I+1,J+1),VY(IM1,J+1),DX)
         VJ   = FQUAD(VY(I,J),VY(I+1,J),VY(IM1,J),DX)
         VJM1 = FQUAD(VY(I,J-1),VY(I+1,J-1),VY(IM1,J-1),DX)
         DV   = FQUAD(VJ,VJP1,VJM1,DY)
C
      END IF
C
C Done
C
      RETURN
      END
C
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STGETC (CNM,CVL)
C
      CHARACTER*(*) CNM,CVL
C
C This subroutine is called to retrieve the character value of a
C specified parameter.
C
C CNM is the name of the parameter whose value is to be retrieved.
C
C CVL is a character variable in which the desired value is to be
C returned by STGETC.
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Check for a parameter name that is too short.
C
c     print *,' ++entree STGETC'
      IF (LEN(CNM).LT.3) THEN
        CSTR(1:36)='STGETC - PARAMETER NAME TOO SHORT - '
        CSTR(37:36+LEN(CNM))=CNM
        CALL SETER (CSTR(1:36+LEN(CNM)),1,1)
        RETURN
      END IF
C
C Get the proper parameter.
C
      IF (CNM(1:3).EQ.'ZFT'.OR.CNM(1:3).EQ.'zft') THEN
         CALL VVTXLN(CZFT,IPCHSZ,IB,IE)
         CVL=CZFT(IB:IE)
      ELSE
C
         CSTR(1:36)='STGETC - PARAMETER NAME NOT KNOWN - '
         CSTR(37:39)=CNM(1:3)
         CALL SETER (CSTR(1:39),3,1)
         RETURN
C
      END IF
C
C
C Done.
C
      RETURN
C
      END
C
C       $Id$
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STGETR (CNM,RVL)
C
      CHARACTER*(*) CNM
C
C This subroutine is called to retrieve the real value of a specified
C parameter.
C
C CNM is the name of the parameter whose value is to be retrieved.
C
C RVL is a real variable in which the desired value is to be returned
C by STGETR.
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Check for a parameter name that is too short.
C
c     print *,' ++entree STGETR'
      IF (LEN(CNM).LT.3) THEN
        CSTR(1:46)='STGETI OR STGETR - PARAMETER NAME TOO SHORT - '
        CSTR(47:46+LEN(CNM))=CNM
        CALL SETER (CSTR(1:46+LEN(CNM)),1,1)
        RETURN
      END IF
C
C Check for incorrect use of the index parameter.
C
      IF (CNM(1:3).EQ.'CLR'.OR.CNM(1:3).EQ.'clr'
     +    .OR.CNM(1:3).EQ.'TVL'.OR.CNM(1:3).EQ.'tvl') THEN
         IF (IPAI.LT.1.OR.IPAI.GT.NLVL) THEN
            CSTR(1:46)='STGETI OR STGETR - GETTING XXX - PAI INCORRECT'
            CSTR(28:30)=CNM(1:3)
            CALL SETER (CSTR(1:46),2,1)
            RETURN
         END IF
      END IF
C
C Get the appropriate parameter value.
C
C ---------------------------------------------------------------------
C
C Values in STPAR
C
      IF (CNM(1:3).EQ.'UD1'.OR. CNM(1:3).EQ.'ud1') THEN
         RVL=REAL(IUD1)
      ELSE IF (CNM(1:3).EQ.'VD1'.OR. CNM(1:3).EQ.'vd1') THEN
         RVL=REAL(IVD1)
      ELSE IF (CNM(1:3).EQ.'PD1'.OR. CNM(1:3).EQ.'pd1') THEN
         RVL=REAL(IPD1)
      ELSE IF (CNM(1:3).EQ.'XD1'.OR. CNM(1:3).EQ.'xd1') THEN
         RVL=REAL(IXD1)
      ELSE IF (CNM(1:3).EQ.'XDM'.OR. CNM(1:3).EQ.'xdm') THEN
         RVL=REAL(IXDM)
      ELSE IF (CNM(1:3).EQ.'YD1'.OR. CNM(1:3).EQ.'yd1') THEN
         RVL=REAL(IYD1)
      ELSE IF (CNM(1:3).EQ.'YDN'.OR. CNM(1:3).EQ.'ydn') THEN
         RVL=REAL(IYDN)
      ELSE IF (CNM(1:3).EQ.'WKD'.OR.CNM(1:3).EQ.'wkd') THEN
        RVL=REAL(IWKD)
      ELSE IF (CNM(1:3).EQ.'WKU'.OR.CNM(1:3).EQ.'wku') THEN
        RVL=REAL(IWKU)
      ELSE IF (CNM(1:3).EQ.'SET'.OR. CNM(1:3).EQ.'set') THEN
         RVL=REAL(ISET)
      ELSE IF (CNM(1:3).EQ.'ERR'.OR. CNM(1:3).EQ.'err') THEN
         RVL=REAL(IERR)
      ELSE IF (CNM(1:3).EQ.'XIN'.OR.CNM(1:3).EQ.'xin') THEN
        RVL=IXIN
      ELSE IF (CNM(1:3).EQ.'YIN'.OR.CNM(1:3).EQ.'yin') THEN
        RVL=IYIN
      ELSE IF (CNM(1:3).EQ.'MSK'.OR. CNM(1:3).EQ.'msk') THEN
         RVL=REAL(IMSK)
      ELSE IF (CNM(1:3).EQ.'CPM'.OR. CNM(1:3).EQ.'cpm') THEN
         RVL=REAL(ICPM)
      ELSE IF (CNM(1:3).EQ.'NLV'.OR.CNM(1:3).EQ.'nlv') THEN
        RVL=REAL(NLVL)
      ELSE IF (CNM(1:3).EQ.'PAI'.OR.CNM(1:3).EQ.'pai') THEN
        RVL=REAL(IPAI)
      ELSE IF (CNM(1:3).EQ.'CTV'.OR.CNM(1:3).EQ.'ctv') THEN
        RVL=REAL(ICTV)
      ELSE IF (CNM(1:3).EQ.'LWD'.OR.CNM(1:3).EQ.'lwd') THEN
        RVL=WDLV
      ELSE IF (CNM(1:3).EQ.'VMN'.OR.CNM(1:3).EQ.'vmn') THEN
        RVL=UVMN
      ELSE IF (CNM(1:3).EQ.'VMX'.OR.CNM(1:3).EQ.'vmx') THEN
        RVL=UVMX
      ELSE IF (CNM(1:3).EQ.'PMN'.OR.CNM(1:3).EQ.'pmn') THEN
        RVL=PMIN
      ELSE IF (CNM(1:3).EQ.'PMX'.OR.CNM(1:3).EQ.'pmx') THEN
        RVL=PMAX
      ELSE IF (CNM(1:3).EQ.'THN'.OR. CNM(1:3).EQ.'thn') THEN
         RVL=REAL(ITHN)
      ELSE IF (CNM(1:3).EQ.'PLR'.OR. CNM(1:3).EQ.'plr') THEN
         RVL=REAL(IPLR)
      ELSE IF (CNM(1:3).EQ.'SST'.OR. CNM(1:3).EQ.'sst') THEN
         RVL=REAL(ISST)
      ELSE IF (CNM(1:3).EQ.'CLR'.OR.CNM(1:3).EQ.'clr') THEN
         RVL=REAL(ICLR(IPAI))
      ELSE IF (CNM(1:3).EQ.'TVL'.OR.CNM(1:3).EQ.'tvl') THEN
         RVL=TVLU(IPAI)
C
C ---------------------------------------------------------------------
C
C Values in STTRAN
C
      ELSE IF (CNM(1:3).EQ.'VPS'.OR. CNM(1:3).EQ.'vps') THEN
         RVL=REAL(UVPS)
      ELSE IF (CNM(1:3).EQ.'VPL'.OR.CNM(1:3).EQ.'vpl') THEN
         RVL=UVPL
      ELSE IF (CNM(1:3).EQ.'VPR'.OR.CNM(1:3).EQ.'vpr') THEN
         RVL=UVPR
      ELSE IF (CNM(1:3).EQ.'VPB'.OR.CNM(1:3).EQ.'vpb') THEN
         RVL=UVPB
      ELSE IF (CNM(1:3).EQ.'VPT'.OR.CNM(1:3).EQ.'vpt') THEN
         RVL=UVPT
      ELSE IF (CNM(1:3).EQ.'WDL'.OR.CNM(1:3).EQ.'wdl') THEN
         RVL=UWDL
      ELSE IF (CNM(1:3).EQ.'WDR'.OR.CNM(1:3).EQ.'wdr') THEN
         RVL=UWDR
      ELSE IF (CNM(1:3).EQ.'WDB'.OR.CNM(1:3).EQ.'wdb') THEN
         RVL=UWDB
      ELSE IF (CNM(1:3).EQ.'WDT'.OR.CNM(1:3).EQ.'wdt') THEN
         RVL=UWDT
      ELSE IF (CNM(1:3).EQ.'XC1'.OR.CNM(1:3).EQ.'xc1') THEN
         RVL=UXC1
      ELSE IF (CNM(1:3).EQ.'XCM'.OR.CNM(1:3).EQ.'xcm') THEN
         RVL=UXCM
      ELSE IF (CNM(1:3).EQ.'YC1'.OR.CNM(1:3).EQ.'yc1') THEN
         RVL=UYC1
      ELSE IF (CNM(1:3).EQ.'YCN'.OR.CNM(1:3).EQ.'ycn') THEN
         RVL=UYCN
C
C ---------------------------------------------------------------------
C
C Values in STSTRM
C
      ELSE IF (CNM(1:3).EQ.'SGD'.OR. CNM(1:3).EQ.'sgd') THEN
         RVL=REAL(ISGD)
      ELSE IF (CNM(1:3).EQ.'AGD'.OR. CNM(1:3).EQ.'agd') THEN
         RVL=REAL(IAGD)
      ELSE IF (CNM(1:3).EQ.'ARL'.OR. CNM(1:3).EQ.'arl') THEN
         RVL=RARL
      ELSE IF (CNM(1:3).EQ.'CKP'.OR. CNM(1:3).EQ.'ckp') THEN
         RVL=REAL(ICKP)
      ELSE IF (CNM(1:3).EQ.'CKX'.OR. CNM(1:3).EQ.'ckx') THEN
         RVL=REAL(ICKX)
      ELSE IF (CNM(1:3).EQ.'TRP'.OR. CNM(1:3).EQ.'trp') THEN
         RVL=REAL(ITRP)
      ELSE IF (CNM(1:3).EQ.'CYK'.OR. CNM(1:3).EQ.'cyk') THEN
         RVL=REAL(ICYK)
      ELSE IF (CNM(1:3).EQ.'VNL'.OR. CNM(1:3).EQ.'vnl') THEN
         RVL=RVNL
      ELSE IF (CNM(1:3).EQ.'SVF'.OR. CNM(1:3).EQ.'svf') THEN
         RVL=REAL(ISVF)
      ELSE IF (CNM(1:3).EQ.'USV'.OR. CNM(1:3).EQ.'usv') THEN
         RVL=RUSV
      ELSE IF (CNM(1:3).EQ.'VSV'.OR. CNM(1:3).EQ.'vsv') THEN
         RVL=RVSV
      ELSE IF (CNM(1:3).EQ.'PSV'.OR. CNM(1:3).EQ.'psv') THEN
         RVL=RPSV
      ELSE IF (CNM(1:3).EQ.'SPC'.OR. CNM(1:3).EQ.'spc') THEN
         RVL=REAL(ISPC)
      ELSE IF (CNM(1:3).EQ.'CDS'.OR. CNM(1:3).EQ.'cds') THEN
         RVL=RCDS
      ELSE IF (CNM(1:3).EQ.'SSP'.OR. CNM(1:3).EQ.'ssp') THEN
         RVL=RSSP
      ELSE IF (CNM(1:3).EQ.'DFM'.OR. CNM(1:3).EQ.'dfm') THEN
         RVL=RDFM
      ELSE IF (CNM(1:3).EQ.'SMD'.OR. CNM(1:3).EQ.'smd') THEN
         RVL=RSMD
      ELSE IF (CNM(1:3).EQ.'AMD'.OR. CNM(1:3).EQ.'amd') THEN
         RVL=RAMD
      ELSE IF (CNM(1:3).EQ.'GBS'.OR. CNM(1:3).EQ.'gbs') THEN
         RVL=REAL(IGBS)
C
C ---------------------------------------------------------------------
C
C Values in STTXP
C
C character attributes
C
C
      ELSE IF (CNM(1:3).EQ.'ZFS'.OR.CNM(1:3).EQ.'zfs') THEN
         RVL=FZFS
      ELSE IF (CNM(1:3).EQ.'ZFX'.OR.CNM(1:3).EQ.'zfx') THEN
         RVL=FZFX
      ELSE IF (CNM(1:3).EQ.'ZFY'.OR.CNM(1:3).EQ.'zfy') THEN
         RVL=FZFY
      ELSE IF (CNM(1:3).EQ.'ZFP'.OR. CNM(1:3).EQ.'zfp') THEN
         RVL=REAL(IZFP)
      ELSE IF (CNM(1:3).EQ.'ZFC'.OR. CNM(1:3).EQ.'zfc') THEN
         RVL=REAL(IZFC)
C
C ---------------------------------------------------------------------
C
C Values in STMAP
C
      ELSE IF (CNM(1:3).EQ.'MAP'.OR. CNM(1:3).EQ.'map') THEN
         RVL=REAL(IMAP)
      ELSE IF (CNM(1:3).EQ.'TRT'.OR. CNM(1:3).EQ.'trt') THEN
         RVL=REAL(ITRT)
      ELSE IF (CNM(1:3).EQ.'VPL'.OR.CNM(1:3).EQ.'vpl') THEN
         RVL=XVPL
      ELSE IF (CNM(1:3).EQ.'VPR'.OR.CNM(1:3).EQ.'vpr') THEN
         RVL=XVPR
      ELSE IF (CNM(1:3).EQ.'VPB'.OR.CNM(1:3).EQ.'vpb') THEN
         RVL=YVPB
      ELSE IF (CNM(1:3).EQ.'VPT'.OR.CNM(1:3).EQ.'vpt') THEN
         RVL=YVPT
      ELSE IF (CNM(1:3).EQ.'XMN'.OR.CNM(1:3).EQ.'xmn') THEN
         RVL=WXMN
      ELSE IF (CNM(1:3).EQ.'XMX'.OR.CNM(1:3).EQ.'xmx') THEN
         RVL=WXMX
      ELSE IF (CNM(1:3).EQ.'YMN'.OR.CNM(1:3).EQ.'ymn') THEN
         RVL=WYMN
      ELSE IF (CNM(1:3).EQ.'YMX'.OR.CNM(1:3).EQ.'ymx') THEN
         RVL=WYMX
      ELSE IF (CNM(1:3).EQ.'XLV'.OR.CNM(1:3).EQ.'xlv') THEN
         RVL=XLOV
      ELSE IF (CNM(1:3).EQ.'XHV'.OR.CNM(1:3).EQ.'xhv') THEN
         RVL=XHIV
      ELSE IF (CNM(1:3).EQ.'YLV'.OR.CNM(1:3).EQ.'ylv') THEN
         RVL=YLOV
      ELSE IF (CNM(1:3).EQ.'YHV'.OR.CNM(1:3).EQ.'yhv') THEN
         RVL=YHIV
      ELSE IF (CNM(1:3).EQ.'NXC'.OR. CNM(1:3).EQ.'nxc') THEN
         RVL=REAL(NXCT)
      ELSE IF (CNM(1:3).EQ.'NYC'.OR. CNM(1:3).EQ.'nyc') THEN
         RVL=REAL(NYCT)
      ELSE IF (CNM(1:3).EQ.'LLG'.OR. CNM(1:3).EQ.'llg') THEN
         RVL=REAL(LNLG)
      ELSE IF (CNM(1:3).EQ.'IVX'.OR. CNM(1:3).EQ.'ivx') THEN
         RVL=REAL(INVX)
      ELSE IF (CNM(1:3).EQ.'IVY'.OR. CNM(1:3).EQ.'ivy') THEN
         RVL=REAL(INVY)
      ELSE IF (CNM(1:3).EQ.'RBG'.OR. CNM(1:3).EQ.'rbg') THEN
         RVL=REAL(RBIG)
      ELSE IF (CNM(1:3).EQ.'IBG'.OR. CNM(1:3).EQ.'ibg') THEN
         RVL=REAL(IBIG)
C
C ---------------------------------------------------------------------
C
      ELSE
         CSTR(1:46)='STGETI OR STGETR - PARAMETER NAME NOT KNOWN - '
         CSTR(47:49)=CNM(1:3)
         CALL SETER (CSTR(1:49),3,1)
         RETURN
      END IF
C
C Done.
C
      RETURN
C
      END
C
C       $Id$
C
      SUBROUTINE STREAM (U,V,P,IAM,STUMSL,WRK)
C
      DIMENSION  U(IUD1,*), V(IVD1,*), P(IPD1,*), IAM(*), WRK(*)
C
      EXTERNAL STUMSL
C
C Input parameters:
C
C U,V    - arrays containing vector field data
C P      - 2-d scalar data array. (dummy - not implemented yet)
C IAM    - An area map array, may be dummied if 'MSK' is zero
C STUMSL - User modifiable masked drawing function; also may
C          be dummied if 'MSK is zero
C WRK    - workspace 
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C -----------------------------------------------------------------
C
C Check for valid area map and area group overflow if masking is enabled
C
c     print *,' ++entree STREAM'
      IF (IMSK.GT.0) THEN
         IF (IAM(7).GT.IPGRCT) THEN
            CSTR(1:29)='STREAM - TOO MANY AREA GROUPS'
            CALL SETER (CSTR(1:29),1,1)
            RETURN
         END IF
         IF (IAM(7).LE.0) THEN
            CSTR(1:25)='STREAM - INVALID AREA MAP'
            CALL SETER (CSTR(1:29),2,1)
            RETURN
         END IF
      END IF
C
C Save the line color, text color and linewidth.
C Then set up the new linewidth values
C 
      CALL GQPLCI(IER,IOC)
      CALL GQTXCI(IER,IOT)
      CALL GQLWSC(IER,ROW)
      CALL GSLWSC(WDLV)
C
C Calculation of NDC sizing values varies based on whether grid 
C relative sizing is in effect.
C
      IF (IGBS .EQ. 0) THEN
         RNDA=RARL*FW2W
         DFMG=RDFM*FW2W
      ELSE
         RNDA=RARL*FW2W/REAL(IXDM)
         DFMG=RDFM*FW2W/REAL(IXDM)
      END IF
C
C If not using the FX,FY routines, then the vector normalization
C value is fixed. 
C
      IF (ICPM.LT.1) THEN
         VNML=0.3333333
      ELSE
         VNML=RVNL
      END IF
C
C Draw the streamlines.
C Break the work array into two parts.  See STDRAW for further
C comments on this.
C
      CALL STDRAW (U,V,WRK(1),WRK(IXDM*IYDN+1),IAM,STUMSL)
C
C Reset the polyline color, text color, and the linewidth
C
      CALL GSPLCI(IOC)
      CALL GSLWSC(ROW)
      CALL GSTXCI(IOT)
C
      RETURN
      END
C
C --------------------------------------------------------------------
C Original disucussion of the STRMLN algorithm follows:
C
C HISTORY                Written and standardized in November 1973.
C
C                        Converted to FORTRAN 77 and GKS in June, 1984.
C
C
C PORTABILITY            FORTRAN 77
C
C ALGORITHM              Wind components are normalized to the value
C                        of DISPL. The least significant two
C                        bits of the work array are
C                        utilized as flags for each grid box. Flag 1
C                        indicates whether any streamline has
C                        previously passed through this box.  Flag 2
C                        indicates whether a directional arrow has
C                        already appeared in a box. Judicious use
C                        of these flags prevents overcrowding of
C                        streamlines and directional arrows.
C                        Experience indicates that a final pleasing
C                        picture is produced when streamlines are
C                        initiated in the center of a grid box. The
C                        streamlines are drawn in one direction then
C                        in the opposite direction.
C
C REFERENCE              The techniques utilized here are described
C                        in an article by Thomas Whittaker (U. of
C                        Wisconsin) which appeared in the notes and
C                        correspondence section of Monthly Weather
C                        Review, June 1977.
C
C TIMING                 Highly variable
C                          It depends on the complexity of the
C                          flow field and the parameters:  DISPL,
C                          DISPC , CSTOP , INITA , INITB , ITERC ,
C                          and IGFLG. (See below for a discussion
C                          of these parameters.) If all values
C                          are default, then a simple linear
C                          flow field for a 40 x 40 grid will
C                          take about 0.4 seconds on the CRAY1-A;
C                          a fairly complex flow field will take about
C                          1.5 seconds on the CRAY1-A.
C
C
C INTERNAL PARAMETERS
C
C                        NAME     DEFAULT         FUNCTION
C                        ----     -------         --------
C
C                        EXT       0.25   Lengths of the sides of the
C                                         plot are proportional to
C                                         IPTSX and JPTSY except in
C                                         the case when MIN(IPTSX,JPT)
C                                         / MAX(IPTSX,JPTSY) .LT. EXT;
C                                         in that case a square
C                                         graph is plotted.
C
C                        SIDE      0.90   Length of longer edge of
C                                         plot. (See also EXT.)
C
C                        XLT       0.05   Left hand edge of the plot.
C                                         (0.0 = left edge of frame)
C                                         (1.0 = right edge of frame)
C
C                        YBT       0.05   Bottom edge of the plot.
C                                         (0.0 = bottom ; 1.0 = top)
C
C                                         (YBT+SIDE and XLT+SIDE must
C                                         be .LE. 1. )
C
C                        INITA     2      Used to precondition grid
C                                         boxes to be eligible to
C                                         start a streamline.
C                                         For example, a value of 4
C                                         means that every fourth
C                                         grid box is eligible ; a
C                                         value of 2 means that every
C                                         other grid box is eligible.
C                                         (see INITB)
C
C                        INITB     2      Used to precondition grid
C                                         boxes to be eligible for
C                                         direction arrows.
C                                         If the user changes the
C                                         default values of INITA
C                                         and/or INITB, it should
C                                         be done such that
C                                         MOD(INITA,INITB) = 0 .
C                                         For a dense grid try
C                                         INITA=4 and INITB=2 to
C                                         reduce the CPU time.
C
C                        AROWL     0.33   Length of direction arrow.
C                                         For example, 0.33 means
C                                         each directional arrow will
C                                         take up a third of a grid
C                                         box.
C
C                        ITERP     35     Every 'ITERP' iterations
C                                         the streamline progress
C                                         is checked.
C
C                        ITERC     -99    The default value of this
C                                         parameter is such that
C                                         it has no effect on the
C                                         code. When set to some
C                                         positive value, the program
C                                         will check for streamline
C                                         crossover every 'ITERC'
C                                         iterations. (The routine
C                                         currently does this every
C                                         time it enters a new grid
C                                         box.)
C                                         Caution:  When this
C                                         parameter is activated,
C                                         CPU time will increase.
C
C                        IGFLG     0      A value of zero means that
C                                         the sixteen point Bessel
C                                         Interpolation Formula will
C                                         be utilized where possible;
C                                         when near the grid edges,
C                                         quadratic and bi-linear
C                                         interpolation  will be
C                                         used. This mixing of
C                                         interpolation schemes can
C                                         sometimes cause slight
C                                         raggedness near the edges
C                                         of the plot.  If IGFLG.NE.0,
C                                         then only the bilinear
C                                         interpolation formula
C                                         is used; this will generally
C                                         result in slightly faster
C                                         plot times but a less
C                                         pleasing plot.
C
C                        IMSG      0      If zero, then no missing
C                                         U and V components are
C                                         present.
C                                         If .NE. 0, STRMLN will
C                                         utilize the
C                                         bi-linear interpolation
C                                         scheme and terminate if
C                                         any data points are missing.
C
C                        UVMSG     1.E+36 Value assigned to a missing
C                                         point.
C
C                        ICYC      0      Zero means the data are
C                                         non-cyclic in the X
C                                         direction.
C                                         If .NE 0, the
C                                         cyclic interpolation
C                                         formulas will be used.
C                                         (Note:  Even if the data
C                                         are cyclic in X, leaving
C                                         ICYC = 0 will do no harm.)
C
C                        DISPL     0.33   The wind speed is
C                                         normalized to this value.
C                                         (See the discussion below.)
C
C                        DISPC     0.67   The critical displacement.
C                                         If after 'ITERP' iterations
C                                         the streamline has not
C                                         moved this distance, the
C                                         streamline will be
C                                         terminated.
C
C                        CSTOP     0.50   This parameter controls
C                                         the spacing between
C                                         streamlines.  The checking
C                                         is done when a new grid
C                                         box is entered.
C
C DISCUSSION OF          Assume a value of 0.33 for DISPL.  This
C DISPL,DISPC            means that it will take three steps to move
C AND CSTOP              across one grid box if the flow was all in the
C                        X direction. If the flow is zonal, then a
C                        larger value of DISPL is in order.
C                        If the flow is highly turbulent, then
C                        a smaller value is in order.  The smaller
C                        DISPL, the more the CPU time.  A value
C                        of 2 to 4 times DISPL is a reasonable value
C                        for DISPC.  DISPC should always be greater
C                        than DISPL. A value of 0.33 for CSTOP would
C                        mean that a maximum of three stream-
C                        lines will be drawn per grid box. This max
C                        will normally only occur in areas of singular
C                        points.
C
C                                            ***************************
C                                            Any or all of the above
C                                            parameters may be changed
C                                            by utilizing common blocks
C                                            STR02 and/or STR03
C                                            ***************************
C
C                        UXSML               A number which is small
C                                            compared to the average
C                                            normalized u component.
C                                            Set automatically.
C
C                        NCHK      750       This parameter is located
C                                            in STDRAW. It specifies the
C                                            length of the circular
C                                            lists  used for checking
C                                            for STRMLN crossovers.
C                                            For most plots this number
C                                            may be reduced to 500
C                                            or less and the plots will
C                                            not be altered.
C
C                        ISKIP               Number of bits to be
C                                            skipped to get to the
C                                            least two significant bits
C                                            in a floating point number.
C                                            The default value is set to
C                                            I1MACH(5) - 2 . This value
C                                            may have to be changed
C                                            depending on the target
C                                            computer; see subroutine
C                                            STDRAW.
C
C --------------------------------------------------------------------
C
C       $Id$
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STRSET
C
C This subroutine may be called to reset all variables which have
C default values to those values.
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Reset individual parameters.
C
C Common block STPAR
C
c     print *,' ++entree STRSET'
      IUD1 = -1
      IVD1 = -1
      IPD1 = -1
      IXD1 = 1
      IXDM = -1
      IYD1 = 1
      IYDN = -1
      IWKD = -1
      IWKU = 0
      ISET = 1
      IERR = 0
      IXIN = 1
      IYIN = 1
      IMSK = 0
      ICPM = 0
      NLVL = 0
      IPAI = 1
      ICTV = 0
      WDLV = 1.0
      UVMN = 0.0
      UVMX = 0.0
      PMIN = 0.0
      PMAX = 0.0
      ITHN = 0
      IMAP = 0
      IPLR = 0
      ISST = 0
C
C Parameter arrays
C
      DO 101 I=1,IPLVLS,1
         ICLR(I) = 1
         TVLU(I) = 0.0
 101  CONTINUE
C
C
C ---------------------------------------------------------------------
C
C STTRAN
C
      UVPS = 0.25
      UVPL = 0.05
      UVPR = 0.95
      UVPB = 0.05
      UVPT = 0.95
      UWDL = 0.0
      UWDR = 0.0
      UWDB = 0.0
      UWDT = 0.0
      UXC1 = 0.0
      UXCM = 0.0
      UYC1 = 0.0
      UYCN = 0.0
C
C ---------------------------------------------------------------------
C
C STSTRM
C
      ISGD = 2
      IAGD = 2
      RARL = 0.012
      ICKP = 35
      ICKX = -99
      ITRP = 0
      ICYK = 0
      RVNL = 0.33
      ISVF = 0
      RUSV = 1.0E12
      RVSV = 1.0E12
      RPSV = 1.0E12
      ISPC = -1
      RCDS = 2.0
      RSSP = 0.015
      RDFM = 0.02
      RSMD = 0.0
      RAMD = 0.0
      IGBS = 0
C
C ---------------------------------------------------------------------
C
C
      FZFS = 0.033
      FZFX = 0.5
      FZFY = 0.5
      IZFP = 0
      IZFC = -1
C
C ---------------------------------------------------------------------
C
C STCHAR values
C
      CZFT = 'ZERO FIELD'
C
C ---------------------------------------------------------------------
C
C STMAP values
C
      IMAP = 0
      ITRT = 1
      IBIG = I1MACH(9)
      RBIG = R1MACH(2)
C
C ---------------------------------------------------------------------
C
C Done
C
      RETURN
C
      END
C
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STSETC (CNM,CVL)
C
      CHARACTER*(*) CNM,CVL
C
C This subroutine is called to give a specified character value to a
C specified parameter.
C
C CNM is the name of the parameter whose value is to be set.
C
C CVL is a character variable containing the new value of the
C parameter.
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Check for a parameter name that is too short.
C
c     print *,' ++entree STSETC'
      IF (LEN(CNM).LT.3) THEN
        CSTR(1:36)='STSETC - PARAMETER NAME TOO SHORT - '
        CSTR(37:36+LEN(CNM))=CNM
        CALL SETER (CSTR(1:36+LEN(CNM)),1,1)
        RETURN
      END IF
C
C Set the proper parameter.
C
      IF (CNM(1:3).EQ.'ZFT'.OR.CNM(1:3).EQ.'zft') THEN
         CZFT=CVL
      ELSE
C
         CSTR(1:36)='STSETC - PARAMETER NAME NOT KNOWN - '
         CSTR(37:39)=CNM(1:3)
         CALL SETER (CSTR(1:39),3,1)
         RETURN
C
      END IF
C
C Done.
C
      RETURN
C
      END
C
C       $Id$
C
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STSETR (CNM,RVL)
C
      CHARACTER*(*) CNM
C
C This subroutine is called to set the real value of a specified
C parameter.
C
C CNM is the name of the parameter whose value is to be set.
C
C RVL is a real variable containing the new value of the parameter.
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
      PARAMETER (IPNPTS = 256, IPLSTL = 1500, IPGRCT = 64)
c     PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Check for a parameter name that is too short.
C
c     print *,' ++entree STSETR'
      IF (LEN(CNM).LT.3) THEN
        CSTR(1:46)='STSETI OR STSETR - PARAMETER NAME TOO SHORT - '
        CSTR(47:46+LEN(CNM))=CNM
        CALL SETER (CSTR(1:46+LEN(CNM)),1,1)
        RETURN
      END IF
C
C Check for incorrect use of the index parameter.
C
      IF (CNM(1:3).EQ.'CLR'.OR.CNM(1:3).EQ.'clr'
     +    .OR.CNM(1:3).EQ.'TVL'.OR.CNM(1:3).EQ.'tvl') THEN
         IF (IPAI.LT.1.OR.IPAI.GT.IPLVLS) THEN
            CSTR(1:46)='STSETI OR STSETR - SETTING XXX - PAI INCORRECT'
            CSTR(28:30)=CNM(1:3)
            CALL SETER (CSTR(1:46),2,1)
            RETURN
         END IF
      END IF
C
C Set the appropriate parameter value.
C
C ---------------------------------------------------------------------
C
C Values in STPAR
C
      IF (CNM(1:3).EQ.'UD1'.OR. CNM(1:3).EQ.'ud1') THEN
         IUD1=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'VD1'.OR. CNM(1:3).EQ.'vd1') THEN
         IVD1=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'PD1'.OR. CNM(1:3).EQ.'pd1') THEN
         IPD1=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'XD1'.OR. CNM(1:3).EQ.'xd1') THEN
         IXD1=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'XDM'.OR. CNM(1:3).EQ.'xdm') THEN
         IXDM=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'YD1'.OR. CNM(1:3).EQ.'yd1') THEN
         IYD1=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'YDN'.OR. CNM(1:3).EQ.'ydn') THEN
         IYDN=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'WKD'.OR.CNM(1:3).EQ.'wkd') THEN
         IWKD=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'WKU'.OR.CNM(1:3).EQ.'wku') THEN
         IWKU=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'SET'.OR. CNM(1:3).EQ.'set') THEN
         ISET=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'ERR'.OR. CNM(1:3).EQ.'err') THEN
         IERR=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'XIN'.OR. CNM(1:3).EQ.'xin') THEN
         IXIN=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'YIN'.OR. CNM(1:3).EQ.'yin') THEN
         IYIN=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'MSK'.OR. CNM(1:3).EQ.'msk') THEN
         IMSK=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'CPM'.OR. CNM(1:3).EQ.'cpm') THEN
         ICPM=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'NLV'.OR.CNM(1:3).EQ.'nlv') THEN
         NLVL=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'PAI'.OR.CNM(1:3).EQ.'pai') THEN
         IF (RVL .LT. 1.0 .OR. RVL .GT. IPLVLS) GO TO 9800
         IPAI=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'CTV'.OR.CNM(1:3).EQ.'ctv') THEN
         ICTV=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'LWD'.OR.CNM(1:3).EQ.'lwd') THEN
         IF (RVL .LE. 0.0) GO TO 9800
         WDLV=RVL
C
C UVMN,UVMX, PMIN, PMAX are read-only
C
      ELSE IF (CNM(1:3).EQ.'THN'.OR. CNM(1:3).EQ.'thn') THEN
         ITHN=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'PLR'.OR. CNM(1:3).EQ.'plr') THEN
         IPLR=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'SST'.OR. CNM(1:3).EQ.'sst') THEN
         ISST=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'CLR'.OR.CNM(1:3).EQ.'clr') THEN
         ICLR(IPAI)=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'TVL'.OR.CNM(1:3).EQ.'tvl') THEN
         TVLU(IPAI)=RVL
C
C ---------------------------------------------------------------------
C
C Values in STTRAN
C
      ELSE IF (CNM(1:3).EQ.'VPS'.OR. CNM(1:3).EQ.'vps') THEN
         UVPS=RVL
      ELSE IF (CNM(1:3).EQ.'VPL'.OR.CNM(1:3).EQ.'vpl') THEN
         UVPL=MIN(1.0,MAX(0.0,RVL))
      ELSE IF (CNM(1:3).EQ.'VPR'.OR.CNM(1:3).EQ.'vpr') THEN
         UVPR=MIN(1.0,MAX(0.0,RVL))
      ELSE IF (CNM(1:3).EQ.'VPB'.OR.CNM(1:3).EQ.'vpb') THEN
         UVPB=MIN(1.0,MAX(0.0,RVL))
      ELSE IF (CNM(1:3).EQ.'VPT'.OR.CNM(1:3).EQ.'vpt') THEN
         UVPT=MIN(1.0,MAX(0.0,RVL))
      ELSE IF (CNM(1:3).EQ.'WDL'.OR.CNM(1:3).EQ.'wdl') THEN
         UWDL=RVL
      ELSE IF (CNM(1:3).EQ.'WDR'.OR.CNM(1:3).EQ.'wdr') THEN
         UWDR=RVL
      ELSE IF (CNM(1:3).EQ.'WDB'.OR.CNM(1:3).EQ.'wdb') THEN
         UWDB=RVL
      ELSE IF (CNM(1:3).EQ.'WDT'.OR.CNM(1:3).EQ.'wdt') THEN
         UWDT=RVL
      ELSE IF (CNM(1:3).EQ.'XC1'.OR.CNM(1:3).EQ.'xc1') THEN
         UXC1=RVL
      ELSE IF (CNM(1:3).EQ.'XCM'.OR.CNM(1:3).EQ.'xcm') THEN
         UXCM=RVL
      ELSE IF (CNM(1:3).EQ.'YC1'.OR.CNM(1:3).EQ.'yc1') THEN
         UYC1=RVL
      ELSE IF (CNM(1:3).EQ.'YCN'.OR.CNM(1:3).EQ.'ycn') THEN
         UYCN=RVL
C
C ---------------------------------------------------------------------
C
C Values in STSTRM
C
      ELSE IF (CNM(1:3).EQ.'SGD'.OR. CNM(1:3).EQ.'sgd') THEN
         ISGD=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'AGD'.OR. CNM(1:3).EQ.'agd') THEN
         IAGD=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'ARL'.OR. CNM(1:3).EQ.'arl') THEN
         RARL=RVL
      ELSE IF (CNM(1:3).EQ.'CKP'.OR. CNM(1:3).EQ.'ckp') THEN
         ICKP=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'CKX'.OR. CNM(1:3).EQ.'ckx') THEN
         ICKX=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'TRP'.OR. CNM(1:3).EQ.'trp') THEN
         ITRP=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'CYK'.OR. CNM(1:3).EQ.'cyk') THEN
         ICYK=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'VNL'.OR. CNM(1:3).EQ.'vnl') THEN
         RVNL=RVL
      ELSE IF (CNM(1:3).EQ.'SVF'.OR. CNM(1:3).EQ.'svf') THEN
         ISVF=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'USV'.OR. CNM(1:3).EQ.'usv') THEN
         RUSV=RVL
      ELSE IF (CNM(1:3).EQ.'VSV'.OR. CNM(1:3).EQ.'vsv') THEN
         RVSV=RVL
      ELSE IF (CNM(1:3).EQ.'PSV'.OR. CNM(1:3).EQ.'psv') THEN
         RPSV=RVL
      ELSE IF (CNM(1:3).EQ.'SPC'.OR. CNM(1:3).EQ.'spc') THEN
         ISPC=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'CDS'.OR. CNM(1:3).EQ.'cds') THEN
         RCDS=RVL
      ELSE IF (CNM(1:3).EQ.'SSP'.OR. CNM(1:3).EQ.'ssp') THEN
         RSSP=RVL
      ELSE IF (CNM(1:3).EQ.'DFM'.OR. CNM(1:3).EQ.'dfm') THEN
         RDFM=RVL
      ELSE IF (CNM(1:3).EQ.'SMD'.OR. CNM(1:3).EQ.'smd') THEN
         RSMD=RVL
      ELSE IF (CNM(1:3).EQ.'AMD'.OR. CNM(1:3).EQ.'amd') THEN
         RAMD=RVL
      ELSE IF (CNM(1:3).EQ.'GBS'.OR. CNM(1:3).EQ.'gbs') THEN
         IGBS=INT(RVL)
C
C This parameter is special in that it causes RSSP,RDFM, and RARL
C to be reset.
C
        IF (IGBS .EQ. 0) THEN
           RARL = 0.012
           RDFM = 0.02
           RSSP = 0.015
        ELSE
           RARL = 0.33
           RDFM = 0.33
           RSSP = 0.5
        END IF
C
C ---------------------------------------------------------------------
C
C Values in STTXP
C
C Character attributes
C
C
      ELSE IF (CNM(1:3).EQ.'ZFS'.OR.CNM(1:3).EQ.'zfs') THEN
         FZFS=RVL
      ELSE IF (CNM(1:3).EQ.'ZFX'.OR.CNM(1:3).EQ.'zfx') THEN
         FZFX=RVL
      ELSE IF (CNM(1:3).EQ.'ZFY'.OR.CNM(1:3).EQ.'zfy') THEN
         FZFY=RVL
      ELSE IF (CNM(1:3).EQ.'ZFP'.OR. CNM(1:3).EQ.'zfp') THEN
         IZFP=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'ZFC'.OR. CNM(1:3).EQ.'zfc') THEN
         IZFC=INT(RVL)
C
C ---------------------------------------------------------------------
C
C Values in STMAP
C
      ELSE IF (CNM(1:3).EQ.'MAP'.OR. CNM(1:3).EQ.'map') THEN
         IMAP=INT(RVL)
      ELSE IF (CNM(1:3).EQ.'TRT'.OR. CNM(1:3).EQ.'trt') THEN
         ITRT=INT(RVL)
C
C ---------------------------------------------------------------------
C
      ELSE
        CSTR(1:46)='STSETI OR STSETR - PARAMETER NAME NOT KNOWN - '
        CSTR(47:49)=CNM(1:3)
        CALL SETER (CSTR(1:49),3,1)
        RETURN
      END IF
C
      GOTO 9900
C
 9800 CONTINUE
C
      CSTR(1:50)='STSETI OR STSETR - PARAMETER VALUE OUT OF RANGE - '
      CSTR(51:53)=CNM(1:3)
      CALL SETER (CSTR(1:53),3,1)
      RETURN
C      
 9900 CONTINUE
C
C Done.
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       $Id$
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STINIT (U,LU,V,LV,P,LP,M,N,WRK,LW)

      USE MODD_RESOLVCAR
C
C Argument dimensions.
C
      DIMENSION       U(LU,N)    ,V(LV,N)    ,P(LP,N)
      DIMENSION       WRK(LW)
C
C Input parameters
C
C U,V   - 2-d arrays holding the component values of a vector field
C LU,LV - The first dimensions of the U and V arrays, respectively
C ----------------
C P     - A 2-d array containing a scalar data field. The contents
C         of this array may be used to color the streamlines. 
C LP    - The first dimension of the P array
C NOTE:
C Coloring by means of the P scalar data field is not yet
C implemented
C ----------------
C M     - The first data dimension (must be less than or equal to
C         MIN(LU,LV) (or MIN(LU,LV,LP) if the P array is used
C WRK   - an internally used work array
C LW    - dimension of the work array (must be at least 2*M*N) 
C
C Output parameters:
C
C None
C
C Force the block data routine, which sets default variables, to load. 
C
      EXTERNAL STDATA
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
C Declare the ST common blocks.
C
      PARAMETER (IPLVLS = 64)
C
C Integer and real common block variables
C
C
      COMMON / STPAR /
     +                IUD1       ,IVD1       ,IPD1       ,
     +                IXD1       ,IXDM       ,IYD1       ,IYDN       ,
     +                IXM1       ,IYM1       ,IXM2       ,IYM2       ,
     +                IWKD       ,IWKU       ,ISET       ,IERR       ,
     +                IXIN       ,IYIN       ,IMSK       ,ICPM       ,
     +                NLVL       ,IPAI       ,ICTV       ,WDLV       ,
     +                UVMN       ,UVMX       ,PMIN       ,PMAX       ,
     +                ITHN       ,IPLR       ,ISST       ,
     +                ICLR(IPLVLS)           ,TVLU(IPLVLS)
C
      COMMON / STTRAN /
     +                UVPS       ,
     +                UVPL       ,UVPR       ,UVPB       ,UVPT       ,
     +                UWDL       ,UWDR       ,UWDB       ,UWDT       ,
     +                UXC1       ,UXCM       ,UYC1       ,UYCN 
C
C Stream algorithm parameters
C
      COMMON / STSTRM /
     +                ISGD       ,IAGD       ,RARL       ,ICKP       ,
     +                ICKX       ,ITRP       ,ICYK       ,RVNL       ,
     +                ISVF       ,RUSV       ,RVSV       ,RNDA       ,
     +                ISPC       ,RPSV       ,RCDS       ,RSSP       ,
     +                RDFM       ,RSMD       ,RAMD       ,IGBS
C
C Text related parameters
C Note: graphical text output is not yet implemented for the
C       Streamline utility.
C
      COMMON / STTXP /
     +                FCWM    ,ICSZ    ,
     +                FMNS    ,FMNX    ,FMNY    ,IMNP    ,IMNC  ,
     +                FMXS    ,FMXX    ,FMXY    ,IMXP    ,IMXC  ,
     +                FZFS    ,FZFX    ,FZFY    ,IZFP    ,IZFC  ,
     +                FILS    ,FILX    ,FILY    ,IILP    ,IILC 
C
C Character variable declartions
C
      CHARACTER*160 CSTR
      PARAMETER (IPCHSZ=80)
      CHARACTER*(IPCHSZ)  CMNT,CMXT,CZFT,CILT
C
C Text string parameters
C
      COMMON / STCHAR / CSTR,CMNT,CMXT,CZFT,CILT
C
      SAVE /STPAR/, /STTRAN/, /STSTRM/, /STTXP/, /STCHAR/
C
C Internal buffer lengths
C
C IPNPTS - Number of points in the point buffer -- not less than 3
C IPLSTL - Streamline-crossover-check circular list length
C IPGRCT - Number of groups supported for area masking
C
c     PARAMETER (IPNPTS = 256, IPLSTL = 15000, IPGRCT = 64)
      PARAMETER (IPNPTS = 256, IPLSTL = 750, IPGRCT = 64)
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
C Write the array sizes into the common block
C 
      IUD1=LU
      IVD1=LV
      IPD1=LP
      IWKD=LW
C
C Error if M > LU or M > LV?
C
      IF (LU.LT.M .OR. LV.LT.M) THEN
         CSTR(1:45)='STINIT - U AND/OR V ARRAY DIMENSIONS EXCEEDED'
         CALL SETER (CSTR(1:45),1,1)
         RETURN
      END IF
      IXDM=MIN(M,LU,LV)
      IYDN=N
      IXD1=1
ccJD
      IYD1=1
      if(nverbia > 0)then
      print *,' **stinit AV ISGD,IYD1,IYDN NSEUIL ',ISGD,
     1IYD1,IYDN,NSEUIL
      endif
      IF(ISGD == 1)THEN
      IYD1=NSEUIL
      ELSE
      IYD1=1
      IF(NSGD==2)IYDN=NSEUIL+1
      ENDIF
c     IYD1=1
      if(nverbia > 0)then
      print *,' **stinit ISGD,IYD1,IYDN ',ISGD,IYD1,IYDN,NSEUIL
      endif
ccJD
      IXM1=IXDM-1
      IXM2=IXDM-2
      IYM1=IYDN-1
      IYM2=IYDN-2
      IF (LW .LT. 2*IXDM*IYDN) THEN
         CSTR(1:37)='STINIT - WRK ARRAY DIMENSION EXCEEDED'
         CALL SETER (CSTR(1:37),2,1)
         RETURN
      END IF
C
C Initialize and transfer some arguments to local variables.
C
      IBIG = I1MACH(9)
      RBIG = R1MACH(2)
C
C Decide what the range of values in X and Y should be.
C
      IF (UXC1.EQ.UXCM) THEN
        XLOV=1.
        XHIV=REAL(IXDM)
      ELSE
        XLOV=UXC1
        XHIV=UXCM
      END IF
C
      IF (UYC1.EQ.UYCN) THEN
        YLOV=1.
        YHIV=REAL(IYDN)
      ELSE
        YLOV=UYC1
        YHIV=UYCN
      END IF
C
      IXIN = MAX(IXIN,1)
      IYIN = MAX(IYIN,1)
C
      NXCT = IXDM/IXIN
      NYCT = IYDN/IYIN
C
C If the user has done a SET call, retrieve the arguments; if he hasn't
C done a SET call, do it for him.
C
      IF (ISET .EQ .0) THEN
C
        CALL GETSET (XVPL,XVPR,YVPB,YVPT,WXMN,WXMX,WYMN,WYMX,LNLG)
C
      ELSE
C
        LNLG=1
C
        IF (UWDL.EQ.UWDR) THEN
          WXMN=XLOV
          WXMX=XHIV
        ELSE
          WXMN=UWDL
          WXMX=UWDR
        END IF
C
        IF (UWDB.EQ.UWDT) THEN
          WYMN=YLOV
          WYMX=YHIV
        ELSE
          WYMN=UWDB
          WYMX=UWDT
        END IF
C
C Determine the viewport based on the setting of the viewport
C shape and viewport extent parameters
C
        IF (UVPS.LT.0.) THEN
          AR=ABS(UVPS)
        ELSE IF (UVPS.EQ.0.) THEN
          AR=(UVPR-UVPL)/(UVPT-UVPB)
        ELSE IF (UVPS.LE.1.) THEN
          AR=ABS((WXMX-WXMN)/(WYMX-WYMN))
          IF (MIN(AR,1./AR).LT.UVPS) AR=(UVPR-UVPL)/(UVPT-UVPB)
        ELSE
          AR=ABS((WXMX-WXMN)/(WYMX-WYMN))
          IF (MAX(AR,1./AR).GT.UVPS) AR=1.
        END IF
C
        IF (AR.LT.(UVPR-UVPL)/(UVPT-UVPB)) THEN
          XVPL=.5*(UVPL+UVPR)-.5*(UVPT-UVPB)*AR
          XVPR=.5*(UVPL+UVPR)+.5*(UVPT-UVPB)*AR
          YVPB=UVPB
          YVPT=UVPT
        ELSE
          XVPL=UVPL
          XVPR=UVPR
          YVPB=.5*(UVPB+UVPT)-.5*(UVPR-UVPL)/AR
          YVPT=.5*(UVPB+UVPT)+.5*(UVPR-UVPL)/AR
        END IF
C
        CALL SET (XVPL,XVPR,YVPB,YVPT,WXMN,WXMX,WYMN,WYMX,LNLG)
C
      END IF
C
C Calculate fraction of VP width to fractional size factor.
C Calculate fraction of VP height to fractional size factor.
C These are for convenience.
C
      FW2W = XVPR - XVPL
      FH2H = YVPT - YVPB
C
C Swap window rectangle if it is inverted, but keep track
C This makes it easier to exclude out-of-bounds points in the
C projection mapping routines
C
      INVX=0
      INVY=0
      IF (WXMN .GT. WXMX) THEN
         T=WXMN
         WXMN=WXMX
         WXMX=T
         INVX=1
      END IF
      IF (WYMN .GT. WYMX) THEN
         T=WYMN
         WYMN=WYMX
         WYMX=T
         INVY=1
      END IF
C
C If cyclic data specified check to ensure the cyclic condition exists.
C The error flag is set if necessary within STCYCL
C
      IF (ICYK.NE.0) CALL STCYCL(U,V)
C
C Calculate the grid size
C
      XGDS=(XHIV-XLOV)/(REAL(NXCT)-1.0)
      YGDS=(YHIV-YLOV)/(REAL(NYCT)-1.0)
C
C Done.
C
      RETURN
C
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C       $Id$
C
C-----------------------------------------------------------------------
C
      SUBROUTINE STUMXY(XDA,YDA,XUS,YUS,IST)
C
C User modifiable routine for mapping data coordinate space to
C user space
C
C
C Input parameters:
C
C XDA,YDA - Point in data coordinate space
C
C Output parameters:
C
C XUS,YUS - Point in user coordinate space
C IST     - Status code indicating success or failure
C
C --------------------------------------------------------------------
      USE MODD_RESOLVCAR
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
C
C*     0.1  Commons
C
      COMMON/TEMV/ZWORKZ,ZZDS,INX,INY
      COMMON/LOGI/LVERT,LHOR,LPT,LXABS
      COMMON/TEMH/ZZXX,ZZXY,IIMAX,IJMAX
      SAVE /TEMH/
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
C -------------------------------------------------------------
c     IMPLICIT NONE
C
C*     0.1  Dummy arguments
C
#include "big.h"
      INTEGER IST 
      REAL XDA,YDA
      REAL XUS,YUS
      REAL ZZXX(N2DVERTX),ZZXY(N2DVERTX)
c     REAL ZZXX(4000),ZZXY(400)
cc    REAL ZZXX(1000),ZZXY(400)
      REAL ZWORKZ(N2DVERTX,2500),ZZDS(N2DVERTX)
c     REAL ZWORKZ(4000,400),ZZDS(4000)
cc    REAL ZWORKZ(1000,400),ZZDS(1000)
      LOGICAL  LVERT,LHOR,LPT,LXABS
      INTEGER INX,INY,IIMAX,IJMAX
C
C*    0.2   Local variables
C
      INTEGER LL,JJ,I,J,IX,IY,IXP1,IYP1
      REAL ZDIFX,ZX1,ZX2,ZY,ZDIFY,ZW1,ZW2,ZW3,ZW4,Z1,Z2,ZR
ccc Avec Interpol en Z
c      INTEGER IPASZ
c      REAL ZPASZ
c      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZW
c      IF(ALLOCATED(ZW))DEALLOCATE(ZW)
c      IPASZ=100
c      ALLOCATE(ZW(IPASZ))
c      ZW(1)=0.
c      ZPASZ=MAXVAL(ZWORKZ) /(IPASZ-1)
c      print *,' IPASZ ZPASZ MAXVAL(XZWORKZ) ',IPASZ,ZPASZ,
c     1MAXVAL(ZWORKZ)
c      DO J=2,IPASZ
c        ZW(J)=ZW(J-1)+ZPASZ
c      ENDDO
c      print *,' **stum LVERT LHOR ',LVERT,LHOR
ccc Avec Interpol en Z
    

C
C Identity transformation
C
      IST=0
c     XUS=XDA
c     YUS=YDA

C*    1.2  Computes streched X's
      IF(IMAP == 4)THEN

      IX=INT(XDA)

      IF(LVERT)THEN
      IF(IX < 1 .OR. IX > INX)THEN
       print *,' **stumxy  AV IX= XDA bizarre IX INX ',XDA,IX,INX
       IX=NINT(XDA)
       print *,' **stumxy  AP IX= XDA bizarre IX ',XDA,IX
      ENDIF

      ELSE

      IF(IX < 1 .OR. IX > IIMAX)THEN
       print *,' **stumxy  AV IX= XDA bizarre IX IIMAX ',XDA,IX,IIMAX
       IX=NINT(XDA)
       print *,' **stumxy  AP IX=XDA bizarre IX ',XDA,IX
      ENDIF
      ENDIF
C     IF(FLOAT(IX)+.989.LE.XDA)IX=IX+1
      ZDIFX=XDA-FLOAT(IX)
c     print *,' XDA IX ZDIFX LHOR+V ',XDA,IX,ZDIFX,LHOR,LVERT

      IF(LVERT)THEN
      ZX1=ZZDS(MAX(IX,1))
      ZX2=ZZDS(MIN(IX+1,INX))
c     PRINT *,' cpmpxy XDA IX',XDA,IX,' ZX1 2',ZX1,ZX2,' XUS ',
c    1XUS
      ELSE

      ZX1=ZZXX(MAX(IX,1))
      ZX2=ZZXX(MIN(IX+1,IIMAX))
      ENDIF
      IF(LVERT)THEN
c     PRINT *,' cpmpxy XDA IX',XDA,IX,' ZX1 2',ZX1,ZX2,' ZDIFX ',
c    1ZDIFX,' INX ',INX
      ELSE
c     PRINT *,' cpmpxy XDA IX',XDA,IX,' ZX1 2',ZX1,ZX2,' ZDIFX ',
c    1ZDIFX,' IIMAX ',IIMAX
      ENDIF
      XUS=ZX1+ZDIFX*(ZX2-ZX1)
c     PRINT *,' cpmpxy XDA IX',XDA,IX,' ZX1 2',ZX1,ZX2,' XUS ',
c    1XUS

C*    1.3  Computes streched Y's

      ZY=YDA
      IY=INT(ZY)
C     IF(FLOAT(IY)+.989.LE.YDA)IY=IY+1
      ZDIFY=ZY-FLOAT(IY)

      IF(LVERT)THEN
c     PRINT *,' cpmpxy YINP IY',YINP,IY
       IXP1=MIN(INX,IX+1)
       IF(LINTERPOLSTR)THEN

ccc Avec Interpol en Z
       IYP1=MIN(NZSTR,IY+1)
       ZW1=XZSTR(MAX(IY,1))
       ZW2=XZSTR(MIN(IYP1,NZSTR))
       ZR=ZW1+ZDIFY*(ZW2-ZW1)
       if(nverbia > 0)then
       print *,' **stum** YDA,IY,ZW1,ZW2,ZR,NZSTR ',
     1YDA,IY,ZW1,ZW2,ZR,NZSTR
       endif
ccc Avec Interpol en Z
ccc SANS Interpol en Z
       ELSE
       IYP1=MIN(INY,IY+1)
       ZW1=ZWORKZ(IX,IY)
       ZW2=ZWORKZ(IX,IYP1)
       ZW3=ZWORKZ(IXP1,IY)
       ZW4=ZWORKZ(IXP1,IYP1)
       Z1=ZW1+ZDIFY*(ZW2-ZW1)
       Z2=ZW3+ZDIFY*(ZW4-ZW3)
       ZR=Z1+ZDIFX*(Z2-Z1)
       if(nverbia > 0)then
       print *,' **stum** YDA,IY,ZW1,ZW2,ZW3,ZW4,Z1,Z2,ZR,INY ',
     1YDA,IY,ZW1,ZW2,ZW3,ZW4,Z1,Z2,ZR,INY
       endif
       ENDIF
ccc SANS Interpol en Z

      ELSE

       ZW1=ZZXY(MAX(IY,1))
       ZW2=ZZXY(MIN(IY+1,IJMAX))
       ZR=ZW1+ZDIFY*(ZW2-ZW1)
      ENDIF
c     PRINT *,' cpmpxy YDA IY',YDA,IY,' ZW1 2',ZW1,ZW2,' ZDIFY ',
c    1ZDIFY,' IJMAX ',IJMAX
       YUS=ZR
       if(nverbia > 0)then
       print *,' ***stumxy... xda,yda,xus,yus ',XDA,YDA,XUS,YUS 
       endif
       ENDIF

C
C Done.
C
      RETURN
C
      END
C
C ---------------------------------------------------------------------
C
C ---------------------------------------------------------------------
C
      SUBROUTINE STUIXY(XUS,YUS,XDA,YDA,IST)
C
C User modifiable routine for inversely transforming
C a point in user coordinate space to data space
C
C Input parameters:
C
C XUS,YUS - Point in user coordinate space
C
C Output parameters:
C
C XDA,YDA - Point in data coordinate space
C IST     - Status code indicating success or failure
C
C --------------------------------------------------------------------
      USE MODN_NCAR
      USE MODD_RESOLVCAR
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
C
C*     0.1  Commons
C
      COMMON/TEMV/ZWORKZ,ZZDS,INX,INY
      COMMON/LOGI/LVERT,LHOR,LPT,LXABS
      COMMON/TEMH/ZZXX,ZZXY,IIMAX,IJMAX
      SAVE /TEMH/
      

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
C ---------------------------------------------------------------------
C
c     IMPLICIT NONE
C
C*     0.1  Dummy arguments
C
#include "big.h"
      INTEGER IST
      REAL XDA,YDA
      REAL XUS,YUS
      REAL ZZXX(N2DVERTX),ZZXY(N2DVERTX)
c     REAL ZZXX(4000),ZZXY(400)
cc    REAL ZZXX(1000),ZZXY(400)
      REAL ZWORKZ(N2DVERTX,2500),ZZDS(N2DVERTX)
c     REAL ZWORKZ(4000,400),ZZDS(4000)
cc    REAL ZWORKZ(1000,400),ZZDS(1000)
      LOGICAL  LVERT,LHOR,LPT,LXABS
      INTEGER INX,INY, IIMAX,IJMAX
      INTEGER IVM,IVM1,IVM2
   
C
C*    0.2   Local variables
C
      INTEGER LL,JJ,I,J,IX,IY,IXP1,IYP1
      REAL ZDIFX,ZX1,ZX2,ZY,ZDIFY,ZW1,ZW2,ZW3,ZW4,Z1,Z2,ZR
      LOGICAL GOK

ccc Avec Interpol en Z
c      INTEGER IPASZ
c      REAL ZPASZ
c      REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZW
c      IF(ALLOCATED(ZW))DEALLOCATE(ZW)
c      IPASZ=100
c      ALLOCATE(ZW(IPASZ))
c      ZW(1)=0.
c      ZPASZ=MAXVAL(ZWORKZ) /(IPASZ-1)
c      print *,' IPASZ ZPASZ MAXVAL(XZWORKZ) ',IPASZ,ZPASZ,
c     1MAXVAL(ZWORKZ)
c      DO J=2,IPASZ
c        ZW(J)=ZW(J-1)+ZPASZ
c      ENDDO
c      print *,' **stum LVERT LHOR ',LVERT,LHOR
ccc Avec Interpol en Z


      IF(IMAP == 4)THEN
      IST=0
c     XDA=XUS
c     YDA=YUS

      IF(LVERT)THEN
      DO I=1,INX
      IVM=0
      IF(XUS == ZZDS(I))THEN
      XDA=I
      IVM=1
      IVM1=I
      GOK=.TRUE.
      EXIT
      ELSEIF(XUS >= ZZDS(MAX(I,1)) .AND. 
     1       XUS < ZZDS(MIN(I+1,INX)))THEN
      ZDIFX=XUS-ZZDS(MAX(I,1))
      ZX1=ZZDS(MAX(I,1))
      ZX2=ZZDS(MIN(I+1,INX))
      XDA=I+ZDIFX/(ZX2-ZX1)
      IVM=2
      IVM1=MAX(I,1)
      IVM2=MIN(I+1,INX)
      GOK=.TRUE.
      EXIT
      ELSE
c     GOK=.FALSE.
      CYCLE
c     IST=-2
      ENDIF
      IF(I == INX)THEN
      XDA=XSPVAL
      IST=-2
      ELSE
      ENDIF
      ENDDO

      ELSE

      DO I=1,IIMAX
      IF(XUS == ZZXX(I))THEN
      XDA=I
      GOK=.TRUE.
      EXIT
      ELSEIF(XUS >= ZZXX(MAX(I,1)) .AND. 
     1       XUS < ZZXX(MIN(I+1,IIMAX)))THEN
      ZDIFX=XUS-ZZXX(MAX(I,1))
      ZX1=ZZXX(MAX(I,1))
      ZX2=ZZXX(MIN(I+1,IIMAX))
      XDA=I+ZDIFX/(ZX2-ZX1)
      GOK=.TRUE.
      EXIT
      ELSE
c     GOK=.FALSE.
      CYCLE
c     IST=-2
      ENDIF
      ENDDO

      ENDIF

      IF(LVERT)THEN

      IF(LINTERPOLSTR)THEN
ccc Avec Interpol en Z
      DO J=1,NZSTR-1
      IF(YUS == XZSTR(J))THEN
      YDA=J
      GOK=.TRUE.
      EXIT
      ELSEIF(YUS >= XZSTR(MAX(J,1)) .AND. 
     1       YUS < XZSTR(MIN(J+1,NZSTR)))THEN
      ZDIFY=YUS-XZSTR(MAX(J,1))
      ZW1=XZSTR(MAX(J,1))
      ZW2=XZSTR(MIN(J+1,NZSTR))
      YDA=J+ZDIFY/(ZW2-ZW1)
      GOK=.TRUE.
      EXIT
      ELSE
c     GOK=.FALSE.
      CYCLE
      ENDIF
      IF(J == NZSTR-1)THEN
      IST=-2
      ENDIF
      ENDDO
ccc Avec Interpol en Z

      ELSE

ccc SANS Interpol en Z
       IF(IVM == 0)THEN
         YDA=XSPVAL
         IST=-2
         RETURN
       ELSEIF(IVM == 1)THEN
         DO J=2,INY-1
           IF(YUS < ZWORKZ(IVM1,2))THEN
             YDA=XSPVAL
             XDA=XSPVAL
             RETURN
           ELSEIF(YUS == ZWORKZ(IVM1,J))THEN
             YDA=J
             EXIT
c          ELSEIF(YUS >= ZWORKZ(IVM1,J) .AND.
           ELSEIF(YUS >= ZWORKZ(IVM1,MAX(J,2)) .AND.
     1         YUS <   ZWORKZ(IVM1,MIN(J+1,INY)))THEN
             ZW1=ZWORKZ(IVM1,MAX(J,2))
             ZW2=ZWORKZ(IVM1,MIN(J+1,INY))
             ZDIFY=YUS-ZW1
             IF(ZW2 /= ZW1)THEN
               YDA=J+ZDIFY/(ZW2-ZW1)
             EXIT
             ELSE
               YDA=J
             EXIT
             ENDIF
           ENDIF
         ENDDO
       ELSEIF(IVM == 2)THEN
       DO J=2,INY-1
         ZW1=ZWORKZ(IVM1,MAX(J,2))
         ZW2=ZWORKZ(IVM1,MIN(J+1,INY))
         ZW3=ZWORKZ(IVM2,MAX(J,2))
         ZW4=ZWORKZ(IVM2,MIN(J+1,INY))
         IF(ZX2 /= ZX1)THEN
           ZW5=ZW1+ZDIFX/(ZX2-ZX1)*(ZW3-ZW1)
         ELSE
           ZW5=ZW1
         ENDIF
         IF(J == 2)THEN
           ZW5M=ZW5
         ENDIF
         IF(ZX2 /= ZX1)THEN
           ZW6=ZW2+ZDIFX/(ZX2-ZX1)*(ZW4-ZW2)
         ELSE
           ZW6=ZW2
         ENDIF
         IF(YUS < ZW5M)THEN
           YDA=XSPVAL
           XDA=XSPVAL
           if(nverbia >0)then
           print *,' stui*** YUS < ZW5M ',YUS,ZW5M
           endif
           RETURN
         ELSEIF(YUS >= ZW5 .AND. YUS < ZW6)THEN
         ZDIFY=YUS-ZW5
         IF(ZW6 /= ZW5)THEN
           YDA=J+ZDIFY/(ZW6-ZW5)
           EXIT
         ELSE
           YDA=J
           EXIT
         ENDIF
         ENDIF
       ENDDO
       ELSE
       YDA=XSPVAL
           if(nverbia >0)then
           print *,' stui*** YUS  en dehors cas prevus ',YUS
           endif
      IST=-2
       RETURN
       ENDIF
ccc SANS Interpol en Z


       ENDIF

      ELSE

      DO J=1,IJMAX
      IF(YUS == ZZXY(J))THEN
      YDA=J
      GOK=.TRUE.
      EXIT
      ELSEIF(YUS >= ZZXY(MAX(J,1)) .AND. 
     1       YUS < ZZXY(MIN(J+1,IJMAX)))THEN
      ZDIFY=YUS-ZZXY(MAX(J,1))
      ZW1=ZZXY(MAX(J,1))
      ZW2=ZZXY(MIN(J+1,IJMAX))
      YDA=J+ZDIFY/(ZW2-ZW1)
      GOK=.TRUE.
      EXIT
      ELSE
c     GOK=.FALSE.
      CYCLE
c     IST=-2
      ENDIF
      ENDDO

      ENDIF
c     print *,' +++STUIXY(XUS,YUS,XDA,YDA,IST) ',XUS,YUS,XDA,YDA,IST
C
C Done
C
      ENDIF
      if(nverbia >0)then
      print *,' +++STUIXY(XUS,YUS,XDA,YDA,IST) ',XUS,YUS,XDA,YDA,IST
      endif
      RETURN
      END
C
C ---------------------------------------------------------------------
C
      SUBROUTINE STUMTA(XDA,YDA,XUS,YUS,XND,YND,DU,DV,TA,IST)
C
C User modifiable routine for mapping a tangent angle in data space to 
C normalized device coordinate space.
C
C Input parameters:
C
C XDA,YDA - Point in data coordinate space
C XUS,YUS - Point in user coordinate space
C XND,YND - Point in NDC space
C DU,DV   - Differential vector components in data space
C
C Output parameters:
C
C TA      - Streamline tangent angle in NDC space
C IST     - Status code indicating success or failure
C
C --------------------------------------------------------------------
C
C The mapping common block: made available to user mapping routines
C
      COMMON /STMAP/
     +                IMAP       ,LNLG       ,INVX       ,INVY       ,
     +                XLOV       ,XHIV       ,YLOV       ,YHIV       ,
     +                WXMN       ,WXMX       ,WYMN       ,WYMX       ,
     +                XVPL       ,XVPR       ,YVPB       ,YVPT       ,
     +                XGDS       ,YGDS       ,NXCT       ,NYCT       ,
     +                ITRT       ,FW2W       ,FH2H       ,
     +                DFMG       ,VNML       ,RBIG       ,IBIG
C
      SAVE /STMAP/
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
C ---------------------------------------------------------------------
C
      IF(IMAP == 4)THEN
      IST=0
      TA=ATAN2(DV,DU)
c     print *,' +++++++STUMTA XDA,YDA,XUS,YUS,XND,YND,DU,DV,TA ',
c    1XDA,YDA,XUS,YUS,XND,YND,DU,DV,TA
      ENDIF
C
C Done.
C
      RETURN
C
      END




