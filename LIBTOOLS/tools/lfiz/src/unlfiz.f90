PROGRAM UNLFIZ
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
INTEGER :: iverb
INTEGER :: inap ! nb d'articles prevus (utile a la creation)
INTEGER :: inaf ! nb d'articles presents dans un fichier existant
INTEGER :: inafdest

CHARACTER(LEN=128) :: filename,DESTFNAME
INTEGER :: JI,JJ
INTEGER :: IRESP
CHARACTER(LEN=FM_FIELD_SIZE),DIMENSION(:),ALLOCATABLE :: yrecfm
INTEGER,                     DIMENSION(:),ALLOCATABLE :: ileng
INTEGER(KIND=8),             DIMENSION(:),ALLOCATABLE :: iwork,iworknew

INTEGER :: ilengs
INTEGER :: ipos
INTEGER :: sizemax

INTEGER            :: ICOMLEN
CHARACTER(LEN=100) :: COMMENT
INTEGER :: DATASIZE,NEWSIZE
INTEGER :: IDIMX,IDIMY,IDIMZ
LOGICAL :: GUSEDIM
INTEGER :: CPT
INTEGER :: LFICOMP
INTEGER :: searchndx
INTEGER :: ITYPCOD
INTEGER :: ITOTAL,ITOTALMAX

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
#if defined(FUJI) || defined(NAGf95) || defined(NEC) || defined(HP) || defined(pgf) || defined(G95) || defined(GFORTRAN)
     CALL GETARG(1+HPINCR,filename)
#else
     CALL PXFGETARG(1,filename,arglen,iresp)
#endif
  ELSE 
     PRINT *,'Usage : ', TRIM(yexe), ' [fichier lfi]'
     STOP
  END IF


searchndx = INDEX(TRIM(filename),".Z.lfi",.TRUE.)
IF (searchndx /= 0 .AND. (LEN_TRIM(filename)-searchndx) == 5) THEN
  PRINT *,'Extension fichier compresse trouvee'
  DESTFNAME=filename(1:searchndx)//'lfi'
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
IF (iresp /= 0 .OR. LFICOMP /= 1) THEN
  PRINT *, 'File ',TRIM(filename),' doesn''t need to be decompressed'
  CALL LFIFER(IRESP,ISRCLU,'KEEP')
  STOP 9
END IF  

PRINT *,'Uncompressed (but 32 bits REAL precision) file : ',DESTFNAME
CALL LFIOUV(IRESP,IDESTLU,.TRUE.,DESTFNAME,'NEW'&
     & ,.FALSE.,.FALSE.,iverb,inaf,inafdest)

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
ITOTALMAX=sizemax
ALLOCATE(IWORKNEW(ITOTALMAX))

CPT=0
DO JI=1,inaf
  CALL LFILEC(IRESP,ISRCLU,yrecfm(JI),iwork,ileng(JI))
  ICOMLEN  = IWORK(2)
  DATASIZE = ileng(JI)-ICOMLEN-2

  CALL GET_COMPHEADER(IWORK(3+ICOMLEN),DATASIZE,NEWSIZE,ITYPCOD)
  IF (NEWSIZE >= 0) THEN 
    
    CPT=CPT+1
    ITOTAL = NEWSIZE+2+ICOMLEN
    PRINT *,'***** ARTICLE compressé ',JI,': ',TRIM(yrecfm(JI)),', taille=',DATASIZE,',decomp=',NEWSIZE
    ! compressed data found
    IF (ITOTALMAX < ITOTAL) THEN
      ITOTALMAX = ITOTAL
      DEALLOCATE(IWORKNEW)
      ALLOCATE(IWORKNEW(ITOTALMAX))
    END IF
    IWORKNEW(1:2+ICOMLEN) = IWORK(1:2+ICOMLEN)
    CALL DECOMPRESS_FIELD(IWORKNEW(3+ICOMLEN),NEWSIZE,IWORK(3+ICOMLEN),DATASIZE,ITYPCOD)
    CALL LFIECR(iresp,IDESTLU,yrecfm(JI),IWORKNEW,ITOTAL)
  ELSE
    PRINT *,'ARTICLE ',JI,': ',TRIM(yrecfm(JI)),', taille =',ileng(JI)
    CALL LFIECR(iresp,IDESTLU,yrecfm(JI),IWORK,ileng(JI))
  END IF
END DO

IF (CPT > 0) THEN
  ! ADD a new article to TAG the compressed file
  IWORK(1) = 0
  COMMENT = "UnCompressed articles"
  ICOMLEN = LEN_TRIM(COMMENT)
  IWORK(2) = ICOMLEN
  DO JJ=1,ICOMLEN
    IWORK(2+JJ)=ICHAR(COMMENT(JJ:JJ))
  END DO
  ILENGS = 3+ICOMLEN
  IWORK(ILENGS) = 2
  CALL LFIECR(iresp,IDESTLU,'LFI_COMPRESSED',iwork,ilengs) 
END IF


PRINT *,' Nombre total d''articles      :', inaf
PRINT *,' Nombre d''articles decompresses :', CPT
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

END PROGRAM UNLFIZ
