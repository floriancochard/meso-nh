!     ######spl
      MODULE MODI_READ_DIMGRIDREF_FM2DIA
!     #####################################
!
INTERFACE
!
SUBROUTINE READ_DIMGRIDREF_FM2DIA(K,HNAMFILE,HLUOUT)
INTEGER :: K
CHARACTER(LEN=*) :: HNAMFILE, HLUOUT
END SUBROUTINE READ_DIMGRIDREF_FM2DIA
!
END INTERFACE
!
END MODULE MODI_READ_DIMGRIDREF_FM2DIA
!     #######################################################
      SUBROUTINE READ_DIMGRIDREF_FM2DIA(K,HNAMFILE,HLUOUT)
!     #######################################################
!
!!****  *READ_DIMGRIDREF_FM2DIA* - Lecture et ecriture des parametres
!!         "intouchables" et des profils 1D de l'etat de reference
!! 
!!
!!    PURPOSE
!!    -------
! 
!
!!**  METHOD
!!    ------
!       Lecture des dimensions par appel a SET_GRID
!          "        parametres de grilles par appel a SET_GRID
!          "        des 3 var. de l'etat de ref. 
!      Ecriture de toutes ces informations dans le fichier diachronique
!                  par appel a WRITE_DIMGRIDREF
!!      
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHORS
!!    -------
!!    J. Duron      *Lab. Aerologie* 
!!
!!    Copyright 1994,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/01/96 
!!      Modification 291196 CSTORAGE_TYPE forced to 'PG' (temp.)
!!      Modification 01/2003 suppression de l appel a SET_REF_FORDIACHRO
!           (=SET_REF modifie en supprimant toute la partie calculs inutile)
!!      Modification 12/2003 appel a SET_GRID remplace par SET_LIGHT_GRID
!!      Modification 09/2004 lecture de MASDEV pour masdev4_6
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIM1  ! NIMAX,NJMAX,NKMAX, NIINF,NISUP, NJINF,NJSUP
USE MODD_GRID  ! XLON0,XLAT0,XBETA, XRPK,XLONORI,XLATORI
USE MODD_GRID1 ! LSLEVE,XLEN1,XLEN2
USE MODD_PARAMETERS, ONLY: JPHEXT,JPVEXT
USE MODD_CONF, ONLY: CCONF,CSTORAGE_TYPE,LCARTESIAN,NMASDEV,NBUGFIX,L1D,L2D,LPACK
USE MODD_PARAM1, ONLY: CSURF
USE MODD_TIME
USE MODD_TIME1
!
USE MODD_DIACHRO, ONLY: CMY_NAME_DIA, CDAD_NAME_DIA
USE MODD_OUT_DIA, ONLY : NLUOUTD
USE MODD_REA_LFI
!
USE MODI_SET_DIM
USE MODI_SET_LIGHT_GRID
USE MODI_FMREAD
!
!*       0.1   Dummy arguments
!

INTEGER           :: K

CHARACTER(LEN=*)  :: HNAMFILE
CHARACTER(LEN=*)  :: HLUOUT
!
!*       0.2   Local variables declarations
!
!
INTEGER           :: JJ, J
INTEGER           :: IIU, IJU, IKU ! Upper bounds in x, y, z directions
INTEGER           :: IIB, IJB, IKB ! Begining useful area in x, y, z directions
INTEGER           :: IIE, IJE, IKE ! End useful area in x, y, z directions
!
REAL              :: ZLAT,ZLON ! Emagram soundings gridpoint location 
                               ! latitude and longitude (decimal degrees)
REAL              :: ZX,ZY     ! Emagram soundings gridpoint location 
                               ! cartesian east and north coordinates (meters)
!
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: ZJ ! Jacobian
!
!-------------------------------------------------------------------------------
!
!*       1.    Preseting the general FM2DIACHRO environment
!              ---------------------------------------
!
!*	 1.1   Sets default values
!
CCONF='POSTP'
!
!*	 1.6   Reads the LFIFM file initial section (i.e. Array dimensions)
!
NIINF=0 ; NISUP=0 ; NJINF=0 ; NJSUP=0
!
CALL SET_DIM(HNAMFILE,HLUOUT,NIINF,NISUP,NJINF,NJSUP,NIMAX,NJMAX,NKMAX)
!
CMY_NAME_DIA(1:LEN(CMY_NAME_DIA))=' '
CRECFM='MY_NAME'
NLENG=28
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,CMY_NAME_DIA,NGRID,NLENCH,CCOMMENT,NRESP)
!
CDAD_NAME_DIA(1:LEN(CDAD_NAME_DIA))=' '
CRECFM='DAD_NAME'
NLENG=28
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,CDAD_NAME_DIA,NGRID,NLENCH,CCOMMENT,NRESP)
print *,'CMY_name CDAD_name ',CMY_NAME_DIA,CDAD_NAME_DIA
!
CRECFM='SURF'
NLENG=4
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,CSURF,NGRID,NLENCH,CCOMMENT,NRESP)
!
!  Reads the geometry configuration selector
!
CRECFM='CARTESIAN'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,LCARTESIAN,NGRID,NLENCH,CCOMMENT,NRESP)
!
CRECFM='THINSHELL'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,LTHINSHELL,NGRID,NLENCH,CCOMMENT,NRESP)
!
CRECFM='STORAGE_TYPE'
NLENG=2
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,CSTORAGE_TYPE,NGRID,NLENCH,CCOMMENT,NRESP)
IF(NRESP /= 0) CSTORAGE_TYPE='MT'
!
CRECFM='L1D'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,L1D,NGRID,NLENCH,CCOMMENT,NRESP)
!
CRECFM='L2D'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,L2D,NGRID,NLENCH,CCOMMENT,NRESP)
!
CRECFM='PACK'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,LPACK,NGRID,NLENCH,CCOMMENT,NRESP)
!
!  Reads the MesoNH version
!
CRECFM='MASDEV'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,NMASDEV,NGRID,NLENCH,CCOMMENT,NRESP)
IF (NRESP /=0 ) NMASDEV=43
!
CRECFM='BUGFIX'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,NBUGFIX,NGRID,NLENCH,CCOMMENT,NRESP)
IF (NRESP /=0 ) NBUGFIX=0
!
CRECFM='JPHEXT'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,JPHEXT,NGRID,NLENCH,CCOMMENT,NRESP)
IF (NRESP /=0 ) JPHEXT=1
!*	 1.7   Allocates the first bunch of input arrays
!
!*       1.7.1  Local variables :
!
IIU=NIMAX+2*JPHEXT ; IJU=NJMAX+2*JPHEXT ; IKU=NKMAX+2*JPVEXT
!
print *,' READ_DIMGRIDREF_FM2DIA CSTORAGE_TYPE=',CSTORAGE_TYPE
IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE=='SU')THEN
  IKU=1
  LCARTESIAN=.FALSE.
  NKMAX=1
ENDIF
!
IIB=1+JPHEXT ; IIE=IIU-JPHEXT
IJB=1+JPHEXT ; IJE=IJU-JPHEXT
IKB=1+JPVEXT ; IKE=IKU-JPVEXT
WRITE(NLUOUTD,*) 'MAIN: IIB, IJB, IKB=',IIB,IJB,IKB
WRITE(NLUOUTD,*) 'MAIN: IIE, IJE, IKE=',IIE,IJE,IKE
WRITE(NLUOUTD,*) 'MAIN: IIU, IJU, IKU=',IIU,IJU,IKU
!
!
IF(K == 1)THEN ! premier fichier
  ALLOCATE(ZJ(IIU,IJU,IKU))
  !
  !*       1.7.2  Grid variables (MODD_GRID1 module):
  !
  ALLOCATE(XXHAT(IIU),XYHAT(IJU),XZHAT(IKU))
  ALLOCATE(XMAP(IIU,IJU))
  ALLOCATE(XLAT(IIU,IJU))
  ALLOCATE(XLON(IIU,IJU))
  ALLOCATE(XDXHAT(IIU),XDYHAT(IJU))
  ALLOCATE(XZS(IIU,IJU),XZSMT(IIU,IJU))
  ALLOCATE(XZZ(IIU,IJU,IKU))
  !
  XXHAT=0. ; XYHAT=0. ; XZHAT=0. ; XMAP=0. ; XLAT=0. ; XLON=0.
  XDXHAT=0. ; XDYHAT=0. ; XZS=0. ; XZZ=0.
  !
ENDIF
!
!*	 1.8   Reads the last section of the LFIFM file
! 
! Notice: The whole XXHAT, XYHAT arrays have to be set here
!         to make provision for any grid selector choice 
!
NIINF=1 ; NISUP=IIU
NJINF=1 ; NJSUP=IJU
!
CALL SET_LIGHT_GRID(1,HNAMFILE,HLUOUT, &
                    IIU,IJU,IKU,NIMAX,NJMAX,         &
                    XLONORI,XLATORI,         &
                    XLON,XLAT,XXHAT,XYHAT,   &
                    XDXHAT,XDYHAT,XMAP,      &
                    XZS,XZZ,XZHAT,LSLEVE,XLEN1,XLEN2,XZSMT,&
                    ZJ,                                    &
                    TDTMOD,TDTCUR                          )
!
IF(CSTORAGE_TYPE == 'PG')THEN
  IKU=1
  LCARTESIAN=.FALSE.
  NKMAX=1
  TDTMOD%TIME=0.
  TDTCUR%TIME=0.
  TDTEXP%TIME=0.
  TDTSEG%TIME=0.
  TDTMOD%TDATE%YEAR=0.
  TDTMOD%TDATE%MONTH=0.
  TDTMOD%TDATE%DAY=0.
  TDTCUR%TDATE%YEAR=0.
  TDTCUR%TDATE%MONTH=0.
  TDTCUR%TDATE%DAY=0.
  TDTEXP%TDATE%YEAR=0.
  TDTEXP%TDATE%MONTH=0.
  TDTEXP%TDATE%DAY=0.
  TDTSEG%TDATE%YEAR=0.
  TDTSEG%TDATE%MONTH=0.
  TDTSEG%TDATE%DAY=0.
ELSE IF(CSTORAGE_TYPE == 'SU')THEN
  IKU=1
  LCARTESIAN=.FALSE.
  NKMAX=1
  TDTMOD%TIME= TDTCUR%TIME
  TDTEXP%TIME= TDTCUR%TIME
  TDTSEG%TIME= TDTCUR%TIME
  TDTMOD%TDATE%YEAR= TDTCUR%TDATE%YEAR
  TDTMOD%TDATE%MONTH= TDTCUR%TDATE%MONTH
  TDTMOD%TDATE%DAY= TDTCUR%TDATE%DAY
  TDTEXP%TDATE%YEAR= TDTCUR%TDATE%YEAR
  TDTEXP%TDATE%MONTH= TDTCUR%TDATE%MONTH
  TDTEXP%TDATE%DAY= TDTCUR%TDATE%DAY
  TDTSEG%TDATE%YEAR= TDTCUR%TDATE%YEAR
  TDTSEG%TDATE%MONTH= TDTCUR%TDATE%MONTH
  TDTSEG%TDATE%DAY= TDTCUR%TDATE%DAY
ENDIF
!
!-------------------------------------------------------------------------------
!
RETURN

END SUBROUTINE READ_DIMGRIDREF_FM2DIA
