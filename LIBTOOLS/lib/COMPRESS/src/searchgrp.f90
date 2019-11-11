!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
MODULE MODE_SEARCHGRP
IMPLICIT NONE 
TYPE SOP_t
  INTEGER :: NBGRP
  INTEGER,DIMENSION(:),POINTER :: IBEG
  INTEGER,DIMENSION(:),POINTER :: IEND
  INTEGER,DIMENSION(:),POINTER :: VALMIN
  INTEGER,DIMENSION(:),POINTER :: VALMAX
END TYPE SOP_t

INTEGER,EXTERNAL :: FMINBITS_IN_WORD

! Private variables
INTEGER,SAVE,                           PRIVATE :: IGRP
INTEGER,DIMENSION(:),ALLOCATABLE,TARGET,PRIVATE :: IBEG,IEND,VALMAX,VALMIN
INTEGER,PARAMETER,                      PRIVATE :: MAINSEUIL=8
INTEGER,SAVE,                           PRIVATE :: IGRPMAX
INTEGER,SAVE,                           PRIVATE :: ICOUNT
INTEGER,DIMENSION(16),PARAMETER,      PRIVATE :: MINELT=(/4,4,4,4,5,5,6,6,7,8,9,11,13,17,26,51/)

! Private routines
PRIVATE :: RECSEARCH_GRP

CONTAINS 
SUBROUTINE INI_SOPDATA(SOPDATA)
TYPE(SOP_t), INTENT(OUT) :: SOPDATA

SOPDATA%NBGRP = 0
NULLIFY(SOPDATA%IBEG)
NULLIFY(SOPDATA%IEND)
NULLIFY(SOPDATA%VALMIN)
NULLIFY(SOPDATA%VALMAX)

END SUBROUTINE INI_SOPDATA

SUBROUTINE RECSEARCH(KTAB,SOPDATA)
INTEGER,DIMENSION(:) :: KTAB
TYPE(SOP_t), INTENT(OUT) :: SOPDATA

INTEGER :: NELT
INTEGER :: GELT,BGELT

IF (ALLOCATED(IBEG)) THEN
  DEALLOCATE(IBEG,IEND,VALMAX,VALMIN)
END IF

NELT=SIZE(KTAB)
ALLOCATE(IBEG(NELT),IEND(NELT),VALMAX(NELT),VALMIN(NELT))
ICOUNT = 0
IGRP   = 0
IGRPMAX = NELT
CALL RECSEARCH_GRP(1,NELT,KTAB,MAINSEUIL)
GELT = MAXVAL(IEND(1:IGRP)-IBEG(1:IGRP)+1)
BGELT = FMINBITS_IN_WORD(GELT)

#ifdef DEBUG
PRINT *,'Routine RECSEARCH_GRP appelee',ICOUNT,'fois.'
PRINT *,'Nbre de groupes =',IGRP
PRINT *,'Nbre maxi d''elements dans groupes',GELT
PRINT *,'Nbre de bits pour coder le nombre d''elements:',BGELT
#endif

SOPDATA%NBGRP=IGRP
SOPDATA%IBEG=>IBEG
SOPDATA%IEND=>IEND
SOPDATA%VALMIN=>VALMIN
SOPDATA%VALMAX=>VALMAX

END SUBROUTINE RECSEARCH

RECURSIVE SUBROUTINE RECSEARCH_GRP(IND1,IND2,ITAB,ISEUIL)
INTEGER,             INTENT(IN) :: IND1,IND2,ISEUIL
INTEGER,DIMENSION(:),INTENT(IN) :: ITAB

INTEGER :: II
INTEGER :: IMAX,IMIN
INTEGER :: IVAL
INTEGER :: nbitcod
INTEGER :: tmpidx1,tmpidx2

ICOUNT=ICOUNT+1

IF (IGRP == 0) THEN
  IMIN = MINVAL(ITAB(IND1:IND2))
  IMAX = MAXVAL(ITAB(IND1:IND2))
  IGRP = 1
  VALMIN(IGRP) = IMIN
  VALMAX(IGRP) = IMAX
  IBEG(IGRP) = IND1
  IEND(IGRP) = IND2
ELSE
  IMIN = VALMIN(IGRP)
  IMAX = VALMAX(IGRP)
END IF

IF (IMAX > IMIN) THEN

  IBEG(IGRP) = IND1
  IEND(IGRP) = IND1
  VALMIN(IGRP) = ITAB(IND1)
  VALMAX(IGRP) = ITAB(IND1)
  
  DO II=IND1,IND2-1
    IVAL = ITAB(II+1)
    IMAX=MAX(VALMAX(IGRP),IVAL)
    IMIN=MIN(VALMIN(IGRP),IVAL)
    IF ((IMAX-IMIN)<(2**ISEUIL)) THEN
      ! II+1 belong to group IGRP
      IEND(IGRP) = II+1
      VALMIN(IGRP) = IMIN
      VALMAX(IGRP) = IMAX
    ELSE
      ! Search the created group
      nbitcod=FMINBITS_IN_WORD(VALMAX(IGRP)-VALMIN(IGRP))
#ifdef DEBUG
      PRINT *,'F:(IGRP,IBEG,IEND,MAX,MIN,nbitcod)=',IGRP,',',IBEG(IGRP),',',IEND(IGRP),',',VALMAX(IGRP),',',VALMIN(IGRP),',',nbitcod
#endif      
      IF (IEND(IGRP)-IBEG(IGRP)>MINELT(nbitcod+1)) THEN
        IF (nbitcod > 0) THEN
          tmpidx1=IBEG(IGRP)
          tmpidx2=IEND(IGRP)
#ifdef DEBUG
          PRINT *,'Appel 1 RECSEARCH_GRP (first,last,seuil):',tmpidx1,tmpidx2,nbitcod/2
#endif
          CALL RECSEARCH_GRP(tmpidx1,tmpidx2,ITAB,nbitcod/2)
        END IF
      ELSE
        IF (IGRP > 1) THEN
          nbitcod=FMINBITS_IN_WORD(VALMAX(IGRP-1)-VALMIN(IGRP-1))
          IMIN=MIN(VALMIN(IGRP-1),VALMIN(IGRP))
          IMAX=MAX(VALMAX(IGRP-1),VALMAX(IGRP))
          IF (IEND(IGRP-1)-IBEG(IGRP-1)<=MINELT(nbitcod+1)) THEN
            IF ((IMAX-IMIN) < 2**15) THEN 
            ! concat IGRP-1 and IGRP
              IEND(IGRP-1) = IEND(IGRP)
              VALMIN(IGRP-1) = IMIN
              VALMAX(IGRP-1) = IMAX
              IGRP = IGRP-1
            END IF
          ELSE
            IF (FMINBITS_IN_WORD(IMAX-IMIN) <= nbitcod) THEN
              ! concat IGRP-1 and IGRP
              IEND(IGRP-1) = IEND(IGRP)
              VALMIN(IGRP-1) = IMIN
              VALMAX(IGRP-1) = IMAX
              IGRP = IGRP-1
            END IF
          END IF
        END IF
      END IF
      ! New group is created
      IGRP = IGRP+1
      IF (IGRP>IGRPMAX) THEN
        PRINT *,'ERROR max number of group exceeded !'
        STOP
      END IF
      IBEG(IGRP) = II+1
      IEND(IGRP) = II+1
      VALMIN(IGRP) = IVAL
      VALMAX(IGRP) = IVAL
    END IF
  END DO
#ifdef DEBUG
  PRINT *,'L:',IGRP,':',VALMAX(IGRP)-VALMIN(IGRP),FMINBITS_IN_WORD(VALMAX(IGRP)-VALMIN(IGRP))
#endif
  nbitcod = FMINBITS_IN_WORD(VALMAX(IGRP)-VALMIN(IGRP))
  IF (IEND(IGRP)-IBEG(IGRP)>= MINELT(nbitcod+1)) THEN
    IF (nbitcod > 0) THEN
      tmpidx1=IBEG(IGRP)
      tmpidx2=IEND(IGRP)
#ifdef DEBUG
      PRINT *,'Appel 2 RECSEARCH_GRP (first,last,seuil):',tmpidx1,tmpidx2,nbitcod/2
#endif
      CALL RECSEARCH_GRP(tmpidx1,tmpidx2,ITAB,nbitcod/2)
    END IF
  END IF
END IF
    
END SUBROUTINE RECSEARCH_GRP

END MODULE MODE_SEARCHGRP

SUBROUTINE INVERTCOL(ITAB,KX,KY)
IMPLICIT NONE 
INTEGER,                  INTENT(IN)   :: KX,KY
INTEGER,DIMENSION(KX,KY), INTENT(INOUT)::ITAB

ITAB(:,2:KY:2) = ITAB(KX:1:-1,2:KY:2)

END SUBROUTINE INVERTCOL

