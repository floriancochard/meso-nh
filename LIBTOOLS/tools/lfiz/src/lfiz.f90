PROGRAM LFIZ
#ifdef NAGf95
  USE F90_UNIX
#endif

IMPLICIT NONE 

#ifndef NAGf95
INTEGER :: IARGC
! CRAY specific
INTEGER :: arglen
!!!!!!!!!!!!!!!!!
#endif
INTEGER :: inarg
CHARACTER(LEN=50) :: yexe


INTEGER, PARAMETER :: FM_FIELD_SIZE = 16
INTEGER, PARAMETER :: ISRCLU  = 11
INTEGER, PARAMETER :: IDESTLU = 12
INTEGER :: JPHEXT
INTEGER :: iverb
INTEGER :: inap ! nb d'articles prevus (utile a la creation)
INTEGER :: inaf ! nb d'articles presents dans un fichier existant
INTEGER :: inafdest

CHARACTER(LEN=128) :: filename,DESTFNAME
INTEGER :: JI,JJ
INTEGER :: IRESP
CHARACTER(LEN=FM_FIELD_SIZE),DIMENSION(:),ALLOCATABLE :: yrecfm
INTEGER,                     DIMENSION(:),ALLOCATABLE :: ileng
INTEGER(KIND=8),             DIMENSION(:),ALLOCATABLE :: iwork

INTEGER :: ilengs
INTEGER :: ipos
INTEGER :: sizemax

INTEGER            :: IGRID
INTEGER            :: ICOMLEN,ICH
CHARACTER(LEN=100) :: COMMENT
INTEGER :: I2DSIZE,I3DSIZE,DATASIZE
INTEGER :: IDIMX,IDIMY,IDIMZ
LOGICAL :: GUSEDIM
INTEGER :: CPT
INTEGER :: LFICOMP
INTEGER :: NEWSIZE
INTEGER :: searchndx
INTEGER :: INDDATIM
INARG = IARGC()

#if defined(F90HP)
#define HPINCR 1
#else
#define HPINCR 0
#endif

#if defined(FUJI) || defined(NAGf95) || defined(NEC) || defined(HP) || defined(pgf) || defined(G95) || defined(GFORTRAN)
  CALL GETARG(0+HPINCR,yexe)
  IF (LEN_TRIM(yexe) == 0) THEN
    PRINT *, 'FATAL ERROR : Activer la macro -DF90HP dans le Makefile et recompiler'
    STOP
  END IF
#else
  CALL PXFGETARG(0,yexe,arglen,iresp)
#endif
!  PRINT *,yexe, ' avec ',INARG,' arguments.'
  IF (INARG == 1) THEN 
#if defined(FUJI) || defined(NAGf95) || defined(NEC) || defined(HP) || defined(pgf) || defined(G95)|| defined(GFORTRAN)
     CALL GETARG(1+HPINCR,filename)
#else
     CALL PXFGETARG(1,filename,arglen,iresp)
#endif
  ELSE 
     PRINT *,'Usage : ', TRIM(yexe), ' [fichier lfi]'
     STOP
  END IF

searchndx = INDEX(TRIM(filename),".lfi",.TRUE.)
IF (searchndx /= 0 .AND. (LEN_TRIM(filename)-searchndx) == 3) THEN
  DESTFNAME=filename(1:searchndx)//'Z.lfi'
ELSE
  PRINT *,'ERROR : extension invalide'
  STOP
END IF

iverb = 0 ! verbosity level


IDIMX = 0
IDIMY = 0
IDIMZ = 0
GUSEDIM = .FALSE.

CALL LFIOUV(IRESP,ISRCLU,.TRUE.,filename,'OLD',.FALSE.&
            & ,.FALSE.,iverb,inap,inaf)

CALL FMREADLFIN1(ISRCLU,'LFI_COMPRESSED',LFICOMP,iresp)
IF (iresp == 0) THEN
  SELECT CASE (LFICOMP)
    CASE(1)
      PRINT *,TRIM(filename),' : already compressed'
    CASE(0)
      PRINT *,'Data are in 32bits real format'
    CASE default
      PRINT *,'File in an unknown compression mode'
    END SELECT
    CALL LFIFER(IRESP,ISRCLU,'KEEP')
    STOP 9
END IF

CALL FMREADLFIN1(ISRCLU,'JPHEXT',JPHEXT,iresp)
IF (iresp /= 0) JPHEXT = 1

! First check if IMAX,JMAX,KMAX exist in LFI file
! to handle 3D, 2D variables -> update IDIMX,IDIMY,IDIMZ
CALL FMREADLFIN1(ISRCLU,'IMAX',IDIMX,iresp)
IF (iresp == 0) IDIMX = IDIMX+2*JPHEXT  ! IMAX + 2*JPHEXT
!
CALL FMREADLFIN1(ISRCLU,'JMAX',IDIMY,iresp)
IF (iresp == 0) IDIMY = IDIMY+2*JPHEXT  ! JMAX + 2*JPHEXT
!
CALL FMREADLFIN1(ISRCLU,'KMAX',IDIMZ,iresp)
IF (iresp == 0) IDIMZ = IDIMZ+2  ! KMAX + 2*JPVEXT

I2DSIZE = IDIMX*IDIMY 
I3DSIZE = IDIMX*IDIMY*IDIMZ
 
GUSEDIM = (I2DSIZE > 0)
IF (GUSEDIM) THEN
  PRINT *,'MESONH 3D, 2D articles DIMENSIONS used :'
  PRINT *,'DIMX =',IDIMX
  PRINT *,'DIMY =',IDIMY
  PRINT *,'DIMZ =',IDIMZ ! IDIMZ may be equal to 0 (PGD files)
ELSE
  PRINT *,'Can''t find IMAX or JMAX variables in the file : Compression ABORTED'
  CALL LFIFER(IRESP,ISRCLU,'KEEP')
  STOP 
END IF


PRINT *,'compressed file : ',DESTFNAME
CALL LFIOUV(IRESP,IDESTLU,.TRUE.,DESTFNAME,'NEW'&
     & ,.FALSE.,.FALSE.,iverb,inaf+1,inafdest)

CALL LFIPOS(IRESP,ISRCLU)
ALLOCATE(yrecfm(inaf))
ALLOCATE(ileng(inaf))
yrecfm(:) = ''
sizemax=0
DO ji=1,inaf
  CALL LFICAS(IRESP,ISRCLU,yrecfm(ji),ileng(ji),ipos,.TRUE.)
  IF (ileng(ji) > sizemax) sizemax=ileng(ji)
END DO
PRINT *,' Nombre total d''articles dans fichier source :', inaf
PRINT *,'sizemax =',sizemax
ALLOCATE(IWORK(sizemax))

CPT=0
DO JI=1,inaf
  CALL LFILEC(IRESP,ISRCLU,yrecfm(JI),iwork,ileng(JI))
  IGRID = IWORK(1)
  ICOMLEN = IWORK(2)
  IF (ICOMLEN > LEN(COMMENT)) THEN
    PRINT *,'ERROR : COMMENT string is too small'
    STOP
  END IF
  
  COMMENT = ''
  DO JJ=1,ICOMLEN
    ICH = iwork(2+JJ)
    COMMENT(JJ:JJ) = CHAR(ICH)
  END DO
  DATASIZE=ileng(JI)-ICOMLEN-2

!  IF (DATASIZE == I2DSIZE .OR. DATASIZE == I3DSIZE) THEN
  !IF (MODULO(DATASIZE,I2DSIZE) == 0) THEN

  INDDATIM=INDEX(yrecfm(JI),'.DATIM')
  IF ((MODULO(DATASIZE,I2DSIZE) == 0).AND. (TRIM(yrecfm(ji))/='ZS').AND.&
  (INDDATIM == 0))THEN
    CPT=CPT+1

!    PRINT *,'GRID=',IGRID
!    PRINT *,'COMMENT = ',TRIM(COMMENT)
!    PRINT *,'Taille data = ',DATASIZE
    PRINT *,'***** compression de ',JI,': ',TRIM(yrecfm(JI))
    CALL COMPRESS_FIELD(IWORK(3+ICOMLEN),IDIMX,IDIMY,DATASIZE,NEWSIZE)
!    NEWSIZE=DATASIZE
    PRINT *,'***** ARTICLE compressé ',JI,': ',TRIM(yrecfm(JI)),', taille=',DATASIZE,',comp=',NEWSIZE
    ileng(JI) = NEWSIZE+ICOMLEN+2
  ELSE
    PRINT *,'ARTICLE ',JI,': ',TRIM(yrecfm(JI)),', taille =',ileng(JI)
  END IF
  CALL LFIECR(iresp,IDESTLU,yrecfm(JI),iwork,ileng(JI))  
END DO

IF (CPT > 0) THEN
  ! ADD a new article to TAG the compressed file
  IWORK(1) = 0
  COMMENT = "Compressed articles"
  ICOMLEN = LEN_TRIM(COMMENT)
  IWORK(2) = ICOMLEN
  DO JJ=1,ICOMLEN
    IWORK(2+JJ)=ICHAR(COMMENT(JJ:JJ))
  END DO
  ILENGS = 3+ICOMLEN
  IWORK(ILENGS) = 1
  CALL LFIECR(iresp,IDESTLU,'LFI_COMPRESSED',iwork,ilengs) 
END IF


PRINT *,' Nombre total d''articles      :', inaf
PRINT *,' Nombre d''articles compresses :', CPT
PRINT *,'sizemax =',sizemax
CALL LFIFER(IRESP,ISRCLU,'KEEP')
CALL LFIFER(IRESP,IDESTLU,'KEEP')

CONTAINS 

SUBROUTINE FMREADLFIN1(klu,hrecfm,kval,kresp)
INTEGER, INTENT(IN)         :: klu ! logical fortran unit au lfi file
CHARACTER(LEN=*),INTENT(IN) :: hrecfm ! article name to be read
INTEGER, INTENT(OUT)        :: kval ! integer value for hrecfm article
INTEGER, INTENT(OUT)        :: kresp! return code null if OK
!
INTEGER(KIND=8),DIMENSION(:),ALLOCATABLE::iwork
INTEGER :: iresp,ilenga,iposex,icomlen
!
CALL LFINFO(iresp,klu,hrecfm,ilenga,iposex)
IF (iresp /=0 .OR. ilenga == 0) THEN
  kresp = -1
  kval = 0
ELSE
  ALLOCATE(IWORK(ilenga))
  CALL LFILEC(iresp,klu,hrecfm,iwork,ilenga)
  icomlen = iwork(2)
  kval = iwork(3+icomlen)
  kresp = iresp
  DEALLOCATE(IWORK)
END IF
END SUBROUTINE FMREADLFIN1

END PROGRAM LFIZ
