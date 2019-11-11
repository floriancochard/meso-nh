!     ######spl
      MODULE MODI_READ_DIMGRIDREF
!     ###########################
!
INTERFACE
!
SUBROUTINE READ_DIMGRIDREF(K,HNAMFILE,HLUOUT)
INTEGER :: K
CHARACTER(LEN=*) :: HNAMFILE, HLUOUT
END SUBROUTINE READ_DIMGRIDREF
!
END INTERFACE
!
END MODULE MODI_READ_DIMGRIDREF
!     ######spl
      SUBROUTINE READ_DIMGRIDREF(K,HNAMFILE,HLUOUT)
!     #############################################
!
!!****  *READ_DIMGRIDREF* - 
!! 
!!
!!    PURPOSE
!!    -------
! 
!
!!**  METHOD
!!    ------
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
!!      Modification 01/2003 suppression de l appel a SET_REF_FORDIACHRO
!           (=SET_REF modifie en supprimant toute la partie calculs inutile)
!!      Modification 12/2003 appel a SET_GRID remplace par SET_LIGHT_GRID
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF, ONLY: CCONF,CSTORAGE_TYPE,LCARTESIAN,LTHINSHELL
USE MODD_DIM1, ONLY: NIMAX,NJMAX,NKMAX, NIINF,NISUP,NJINF,NJSUP
USE MODD_GRID  ! XLONORI,XLATORI
USE MODD_GRID1, ONLY: XLON,XLAT,XXHAT,XYHAT,&
                      XDXHAT,XDYHAT,XMAP,XZS,XZZ,XZHAT,&
                      LSLEVE,XLEN1,XLEN2,XZSMT
USE MODD_PARAMETERS, ONLY: JPHEXT,JPVEXT
USE MODD_TIME
USE MODD_TIME1
!
USE MODD_REA_LFI
USE MODD_RESOLVCAR, ONLY: NVERBIA
!
USE MODI_SET_DIM
USE MODI_SET_LIGHT_GRID
USE MODI_FMREAD
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!
INTEGER           :: K
!
CHARACTER(LEN=*)  :: HNAMFILE
CHARACTER(LEN=*)  :: HLUOUT
!
!*       0.2   Local variables declarations
!
INTEGER           :: IIU, IJU, IKU ! Upper bounds in x, y, z directions
INTEGER           :: IIB, IJB, IKB ! Begining useful area in x, y, z directions
INTEGER           :: IIE, IJE, IKE ! End useful area in x, y, z directions
!
INTEGER,SAVE      :: IIINF, IISUP, IJINF, IJSUP
!
!REAL              :: ZLAT,ZLON ! Emagram soundings gridpoint location 
                               ! latitude and longitude (decimal degrees)
!REAL              :: ZX,ZY     ! Emagram soundings gridpoint location 
                               ! cartesian east and north coordinates (meters)
REAL,DIMENSION(:,:,:),ALLOCATABLE :: ZJ ! Jacobian
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
IIINF=NIINF; IISUP=NISUP; IJINF=NJINF; IJSUP=NJSUP
NIINF=0    ; NISUP=0    ; NJINF=0    ; NJSUP=0
NIMAX=0
CALL FMREAD(HNAMFILE,'IMAX',HLUOUT,1,NIMAX,NGRID,NLENCH,CCOMMENT,NRESP)
IF(NRESP /= 0)THEN
  NIMAX=0
  print *,' Absence d''entete dans ce fichier '
  RETURN
ENDIF
if(nverbia>=5) print *,'Av SET_DIM NIMAX=',NIMAX
CALL SET_DIM(HNAMFILE,HLUOUT,NIINF,NISUP,NJINF,NJSUP,NIMAX,NJMAX,NKMAX)
if(nverbia>=5) print *,'Ap SET_DIM NIMAX=',NIMAX
!
!  Reads the geometry configuration selector
!
CRECFM='CARTESIAN'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,LCARTESIAN,NGRID,NLENCH,CCOMMENT,NRESP)
if(nverbia>=5)print *,' LCARTESIAN=', LCARTESIAN

CRECFM='THINSHELL'
NLENG=1
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,LTHINSHELL,NGRID,NLENCH,CCOMMENT,NRESP)
if(nverbia>=5)print *,' LTHINSHELL=', LTHINSHELL

CRECFM='STORAGE_TYPE'
NLENG=2
CALL FMREAD(HNAMFILE,CRECFM,HLUOUT,NLENG,CSTORAGE_TYPE,NGRID,NLENCH,CCOMMENT,NRESP)
IF(NRESP /= 0) CSTORAGE_TYPE='MT'
print *,' CSTORAGE_TYPE =',CSTORAGE_TYPE
!
!*	 1.7   Allocates the first bunch of input arrays
!
!
!*       1.7.1  Local variables :
!
IIU=NIMAX+2*JPHEXT ; IJU=NJMAX+2*JPHEXT ; IKU=NKMAX+2*JPVEXT
!
IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE =='SU') IKU=1
!
IIB=1+JPHEXT ; IIE=IIU-JPHEXT
IJB=1+JPHEXT ; IJE=IJU-JPHEXT
IKB=1+JPVEXT ; IKE=IKU-JPVEXT
if(nverbia>=3) print*,'* in READ_DIMGRIDREF'
print*,' IIB, IJB, IKB= ',IIB,IJB,IKB
print*,' IIE, IJE, IKE= ',IIE,IJE,IKE
print*,' IIU, IJU, IKU= ',IIU,IJU,IKU
!
!*       1.7.2  Grid variables (MODD_GRID1 module):
!
IF(ALLOCATED(XXHAT)) DEALLOCATE(XXHAT)
IF(ALLOCATED(XYHAT)) DEALLOCATE(XYHAT)
IF(ALLOCATED(XZHAT)) DEALLOCATE(XZHAT)
IF(ALLOCATED(XMAP))  DEALLOCATE(XMAP)
IF(ALLOCATED(XLAT))  DEALLOCATE(XLAT)
IF(ALLOCATED(XLON))  DEALLOCATE(XLON)
IF(ALLOCATED(XDXHAT))DEALLOCATE(XDXHAT)
IF(ALLOCATED(XDYHAT))DEALLOCATE(XDYHAT)
IF(ALLOCATED(XZS))   DEALLOCATE(XZS)
IF(ALLOCATED(XZSMT)) DEALLOCATE(XZSMT)
IF(ALLOCATED(XZZ))   DEALLOCATE(XZZ)
ALLOCATE(XXHAT(IIU),XYHAT(IJU),XZHAT(IKU))
ALLOCATE(XMAP(IIU,IJU))
ALLOCATE(XLAT(IIU,IJU))
ALLOCATE(XLON(IIU,IJU))
ALLOCATE(XDXHAT(IIU),XDYHAT(IJU))
ALLOCATE(XZS(IIU,IJU),XZSMT(IIU,IJU))
ALLOCATE(XZZ(IIU,IJU,IKU))
!
!*	 1.8   Reads the last section of the LFIFM file
! 
! Notice: The whole XXHAT, XYHAT arrays have to be set here
!         to make provision for any grid selector choice 
!
NIINF=1 ; NISUP=IIU
NJINF=1 ; NJSUP=IJU
!
ALLOCATE(ZJ(IIU,IJU,IKU))
CALL SET_LIGHT_GRID(1,HNAMFILE,HLUOUT, &
                    IIU,IJU,IKU,NIMAX,NJMAX,   &
                    XLONORI,XLATORI,           &
                    XLON,XLAT,XXHAT,XYHAT,             &
                    XDXHAT,XDYHAT,XMAP,        &
                    XZS,XZZ,XZHAT,LSLEVE,XLEN1,XLEN2,XZSMT,  &
                    ZJ,                                &
                    TDTMOD,TDTCUR )
!
DEALLOCATE(ZJ)
IF(IIINF /= 0 .AND. IISUP /=0 .AND. IJINF /=0 .AND. IJSUP /=0)THEN
  NIINF=IIINF; NISUP=IISUP; NJINF=IJINF; NJSUP=IJSUP
ENDIF
!
!------------------------------------------------------------------------------
!
!*      4.    EPILOGUE
!             --------
RETURN

END SUBROUTINE READ_DIMGRIDREF
