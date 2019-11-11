!     ######spl
      MODULE MODI_READ_AND_WRITE_DIMGRIDREF
!     #####################################
!
INTERFACE
!
SUBROUTINE READ_AND_WRITE_DIMGRIDREF(K,HNAMFILE,HLUOUT)
INTEGER :: K
CHARACTER(LEN=*) :: HNAMFILE, HLUOUT
END SUBROUTINE READ_AND_WRITE_DIMGRIDREF
!
END INTERFACE
!
END MODULE MODI_READ_AND_WRITE_DIMGRIDREF
!     #######################################################
      SUBROUTINE READ_AND_WRITE_DIMGRIDREF(K,HNAMFILE,HLUOUT)
!     #######################################################
!
!!****  *READ_AND_WRITE_DIMGRIDREF* - Lecture et ecriture des parametres
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIACHRO, ONLY: CMY_NAME_DIA, CDAD_NAME_DIA
USE MODD_DIM1  ! NIMAX,NJMAX,NKMAX, NIINF,NISUP, NJINF,NJSUP
USE MODD_DIMGRID_FORDIACHRO, ONLY: NNBF
USE MODD_GRID  ! XLON0,XLAT0, XBETA,XRPK
USE MODD_GRID1
USE MODD_OUT_DIA, ONLY : NLUOUTD
USE MODD_OUT1
USE MODD_PARAMETERS
USE MODD_DYN , ONLY: XSEGLEN
USE MODD_DYN1, ONLY: XTSTEP
USE MODD_CONF, ONLY: CCONF,CSTORAGE_TYPE,LCARTESIAN
USE MODD_TIME
USE MODD_TIME1
USE MODD_REF  ! XRHODREFZ,XTHVREFZ,XEXNTOP
USE MODD_REA_LFI
!
USE MODI_SET_DIM
USE MODI_SET_GRID
USE MODI_WRITE_DIMGRIDREF
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
REAL,DIMENSION(:), ALLOCATABLE,SAVE  :: IIMAX, IJMAX, IKMAX, ITIMECUR
REAL,DIMENSION(:), ALLOCATABLE,SAVE  :: ZLON0, ZRPK, ZLONOR, ZLATOR, ZLAT0, &
                                        ZBETA
LOGICAL,DIMENSION(:), ALLOCATABLE,SAVE :: OCARTESIAN
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
!  Reads the geometry configuration selector
CRECFM='THINSHELL'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,LTHINSHELL,NGRID,NLENCH,CCOMMENT,NRESP)
!
CRECFM='CARTESIAN'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,LCARTESIAN,NGRID,NLENCH,CCOMMENT,NRESP)
!
CRECFM='STORAGE_TYPE'
NLENG=2
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,CSTORAGE_TYPE,NGRID,NLENCH,CCOMMENT,NRESP)
IF(NRESP /= 0) CSTORAGE_TYPE='MT'
!
!
!*	 1.7   Allocates the first bunch of input arrays
!
!*       1.7.1  Local variables :
!
IIU=NIMAX+2*JPHEXT ; IJU=NJMAX+2*JPHEXT ; IKU=NKMAX+2*JPVEXT
!
print *,' READ_AND_WRITE_DIMGRIDREF ENTREE CSTORAGE_TYPE ',CSTORAGE_TYPE
IF(CSTORAGE_TYPE == 'PG')THEN
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
  ALLOCATE(XZS(IIU,IJU))
  ALLOCATE(XZZ(IIU,IJU,IKU))
  !
  !*       1.7.3  Reference state variables (MODD_REF1 module):
  !
  ALLOCATE(XRHODREFZ(IKU),XTHVREFZ(IKU))
  !
  XXHAT=0. ; XYHAT=0. ; XZHAT=0. ; XMAP=0. ; XLAT=0. ; XLON=0.
  XDXHAT=0. ; XDYHAT=0. ; XZS=0. ; XZZ=0.
  XRHODREFZ=0. ; XTHVREFZ=0.; XEXNTOP=0.
  ALLOCATE(IIMAX(NNBF),IJMAX(NNBF),IKMAX(NNBF),ITIMECUR(NNBF))
  ALLOCATE(ZLON0(NNBF),ZRPK(NNBF),ZLONOR(NNBF),ZLATOR(NNBF),ZLAT0(NNBF),ZBETA(NNBF))
  ALLOCATE(OCARTESIAN(NNBF))
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
! Artifice pour eviter 1 plantage dans SET_GRID
XTSTEP=50.
XSEGLEN=500.
!
CALL SET_GRID(1,HNAMFILE,HLUOUT,IIU,IJU,IKU,NIINF,NISUP,NJINF,NJSUP,XTSTEP,&
              XSEGLEN, XOUT1,XOUT2,XOUT3,XOUT4,XOUT5,XOUT6,XOUT7,XOUT8,    &
                       XOUT9,XOUT10,XOUT11,XOUT12,XOUT13,XOUT14,XOUT15,    &
                       XOUT16,XOUT17,XOUT18,XOUT19,XOUT20,                 &
              XLONOR,XLATOR,XLON,XLAT,XXHAT,XYHAT,                         &
              XDXHAT,XDYHAT,XMAP,XZS,XZZ,XZHAT,                            &
              ZJ,                                                          &
              TDTMOD,TDTCUR,NSTOP,NOUT_TIMES,NOUT_NUMB                     )
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
ENDIF
!
!*       1.9   read 3 variables of ref. state without orography (SET_REF)
!
CRECFM='STORAGE_TYPE'
NLENG=2
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,CSTORAGE_TYPE,NGRID,NLENCH,CCOMMENT,NRESP)
!
CRECFM='RHOREFZ'
NLENG=SIZE(XRHODREFZ)
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,XRHODREFZ,NGRID,NLENCH,CCOMMENT,NRESP)
IF(NRESP == -47)THEN
  print *,' XRHODREFZ ABSENT dans le fichier ',TRIM(HNAMFILE),': MIS a 0. '
  XRHODREFZ(:)=0.
ENDIF
!
CRECFM='THVREFZ'
NLENG=SIZE(XTHVREFZ)
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,XTHVREFZ,NGRID,NLENCH,CCOMMENT,NRESP)
IF(NRESP == -47)THEN
  print *,' XTHVREFZ ABSENT dans le fichier ',TRIM(HNAMFILE),': MIS a 0. '
  XTHVREFZ(:)=0.
ENDIF
!
CRECFM='EXNTOP'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,XEXNTOP,NGRID,NLENCH,CCOMMENT,NRESP)
IF(NRESP == -47)THEN
  print *,' XEXNTOP ABSENT dans le fichier ',TRIM(HNAMFILE),': MIS a 0. '
  XEXNTOP=0.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    WRITING OR CHECKING DIM., GRID., REF. VARIABLES
!              -----------------------------------------------
!
IIMAX(K)=NIMAX ; IJMAX(K)=NJMAX ; IKMAX(K)=NKMAX
ITIMECUR(K)=TDTCUR%TIME
!
ZLON0(K)=XLON0   ; ZLAT0(K)=XLAT0
ZLONOR(K)=XLONOR ; ZLATOR(K)=XLATOR
ZRPK(K)=XRPK     ; ZBETA(K)=XBETA
!
OCARTESIAN(K)=LCARTESIAN
!
!
IF(K == 1)THEN  ! premier fichier
  !
  CALL WRITE_DIMGRIDREF
  !
ENDIF
!
IF(K > 1)THEN   ! fichiers suivants
  !
  IF(IIMAX(K) /= IIMAX(1))THEN
    PRINT *,' K IIMAX(K) IIMAX(1) ',K,IIMAX(K),IIMAX(1)
  ENDIF
  IF(IJMAX(K) /= IJMAX(1))THEN
    PRINT *,' K IJMAX(K) IJMAX(1) ',K,IJMAX(K),IJMAX(1)
  ENDIF
  IF(IKMAX(K) /= IKMAX(1))THEN
    PRINT *,' K IKMAX(K) IKMAX(1) ',K,IKMAX(K),IKMAX(1)
  ENDIF
  IF(ITIMECUR(K) /= ITIMECUR(1))THEN
    PRINT *,' K ITIMECUR(K) ITIMECUR(1) ',K,ITIMECUR(K),ITIMECUR(1)
  ENDIF
  !
  IF(ZLON0(K) /= ZLON0(1))THEN
    PRINT *,' K ZLON0(K) ZLON0(1) ',K,ZLON0(K),ZLON0(1)
  ENDIF
  IF(ZRPK(K) /= ZRPK(1))THEN
    PRINT *,' K ZRPK(K) ZRPK(1) ',K,ZRPK(K),ZRPK(1)
  ENDIF
  IF(ZLONOR(K) /= ZLONOR(1))THEN
    PRINT *,' K ZLONOR(K) ZLONOR(1) ',K,ZLONOR(K),ZLONOR(1)
  ENDIF
  IF(ZLATOR(K) /= ZLATOR(1))THEN
    PRINT *,' K ZLATOR(K) ZLATOR(1) ',K,ZLATOR(K),ZLATOR(1)
  ENDIF
  IF(ZLAT0(K) /= ZLAT0(1))THEN
    PRINT *,' K ZLAT0(K) ZLAT0(1) ',K,ZLAT0(K),ZLAT0(1)
  ENDIF
  IF(ZBETA(K) /= ZBETA(1))THEN
    PRINT *,' K ZBETA(K) ZBETA(1) ',K,ZBETA(K),ZBETA(1)
  ENDIF
  !
  IF((OCARTESIAN(K) .AND..NOT. OCARTESIAN(1)) .OR. &
     (.NOT. OCARTESIAN(K) .AND. OCARTESIAN(1)))THEN
    PRINT *,' K OCARTESIAN(K) OCARTESIAN(1) ',K,OCARTESIAN(K),OCARTESIAN(1)
  ENDIF
  !
ENDIF
!------------------------------------------------------------------------------
!
!*      4.    EPILOG
!             ------
!
IF(K == NNBF)THEN  ! dernier fichier
  DEALLOCATE(IIMAX,IJMAX,IKMAX,ITIMECUR)
  DEALLOCATE(ZLON0,ZRPK,ZLONOR,ZLATOR,ZLAT0,ZBETA)
  DEALLOCATE(OCARTESIAN)
END IF
!
RETURN

END SUBROUTINE READ_AND_WRITE_DIMGRIDREF
