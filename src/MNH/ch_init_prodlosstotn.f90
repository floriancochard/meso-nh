!!    ########################### 
      MODULE MODI_CH_INIT_PRODLOSSTOT_n
!!    ########################### 
!!
INTERFACE
!!
SUBROUTINE CH_INIT_PRODLOSSTOT_n(KLUOUT)
!!
IMPLICIT NONE
!!
INTEGER, INTENT(IN)  :: KLUOUT   ! output listing channel
!!
!!
END SUBROUTINE CH_INIT_PRODLOSSTOT_n
!!
END INTERFACE
!!
END MODULE MODI_CH_INIT_PRODLOSSTOT_n
!!
!!    ########################################
      SUBROUTINE CH_INIT_PRODLOSSTOT_n(KLUOUT)
!!    ########################################
!!
!!****  *CH_INIT_PRODLOSSTOT_n* - prepare arrays for chemical prod/loss terms computation
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to search for which species the integrated
!     production/destruction is desired (DIAG1.nam) and to allocate needed arrays
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH
!!          
!!    AUTHOR
!!    ------
!!	    F. Brosse *Laboratoire d'Aerologie UPS-CNRS*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    October 2016                   
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
USE MODD_CH_MNHC_n,  ONLY: LUSECHAQ, CSPEC_PRODLOSS
USE MODD_CH_M9_n,    ONLY: CNAMES, NEQ
USE MODD_CH_PRODLOSSTOT_n, ONLY: NEQ_PLT, CNAMES_PRODLOSST, NIND_SPEC
USE MODD_DIAG_FLAG, ONLY : LCHEMDIAG, CSPEC_DIAG
IMPLICIT NONE 

INTEGER, INTENT(IN)  :: KLUOUT   ! output listing channel
!local variables
CHARACTER(LEN(CSPEC_PRODLOSS)+2) :: YWORKSTR, YS ! YWORKSTR is CCSPEC surrounded by ','
CHARACTER(LEN=NEQ*34)   :: YCHECKSTR
CHARACTER(LEN=32)       :: YSPEC_NAME
LOGICAL                 :: GCHECKFAILED, GFOUND
INTEGER                 :: I
INTEGER                 :: J
INTEGER                 :: ILEN

NEQ_PLT = 0
ALLOCATE(CNAMES_PRODLOSST(NEQ))
ALLOCATE(NIND_SPEC(NEQ))

IF(LCHEMDIAG.AND.LEN(TRIM(CSPEC_DIAG))/=0) CSPEC_PRODLOSS=CSPEC_DIAG

YWORKSTR = CLEANSTR(CSPEC_PRODLOSS)
WRITE (KLUOUT,*) 'YWORKSTR=',TRIM(YWORKSTR)
IF (YWORKSTR /= '') THEN
   YCHECKSTR=','
   DO I=1,NEQ
      ! Create check string
      YCHECKSTR=TRIM(YCHECKSTR)//TRIM(CNAMES(I))//','

      GFOUND = .FALSE.
      IF (INDEX(YWORKSTR,','//TRIM(CNAMES(I))//',') /= 0) THEN
         GFOUND = .TRUE.
      ELSE
         IF (LUSECHAQ) THEN
            IF ((INDEX(CNAMES(I),'WC_') == 1 .OR. INDEX(CNAMES(I),'WR_') == 1)) THEN
               IF (INDEX(YWORKSTR,','//TRIM(CNAMES(I)(4:))//',') /= 0) THEN
                  GFOUND = .TRUE.
               END IF
            END IF
           
         END IF
      END IF
      IF (GFOUND) THEN
         NEQ_PLT = NEQ_PLT + 1
         CNAMES_PRODLOSST(NEQ_PLT) = TRIM(CNAMES(I))
         NIND_SPEC(NEQ_PLT)  = I
      END IF
   END DO
  
   
   !WRITE(KLUOUT, *) 'YCHECKSTR=',TRIM(YCHECKSTR)
   ! Compare inputs with CNAMES check string
   GCHECKFAILED = .FALSE.
   YS = YWORKSTR(2:)
   DO
      I = SCAN(YS,',')
      IF (I == 0) EXIT
      YSPEC_NAME=YS(:I-1)
      YS = YS(I+1:)
      IF (INDEX(YCHECKSTR,','//TRIM(YSPEC_NAME)//',') == 0) THEN
         WRITE (KLUOUT,*) TRIM(YSPEC_NAME),' not found in CNAMES'
         GCHECKFAILED = .TRUE.
      END IF
   END DO
   IF (GCHECKFAILED) THEN
      WRITE(KLUOUT,*) 'Wrong (misspelled) CSPEC encountered...ABORTING !'
      CALL ABORT
      STOP 
   END IF
ELSE
   DEALLOCATE(CNAMES_PRODLOSST)
   DEALLOCATE(NIND_SPEC)
   ALLOCATE(CNAMES_PRODLOSST(0))
   ALLOCATE(NIND_SPEC(0))
END IF

IF (NEQ_PLT > 0) THEN
   PRINT *,'     CNAMES_PRODLOSSTOT    ISPEC'
   DO I=1,NEQ_PLT
      PRINT "(I2,'      ',A10,'     ',I3)", I, CNAMES_PRODLOSST(I), NIND_SPEC(I)
   END DO
ELSE
   PRINT *, "Nothing to do !"
END IF


CONTAINS 

FUNCTION CLEANSTR(ins) RESULT(outs)
CHARACTER(LEN=*), INTENT(IN) :: ins
CHARACTER(LEN=LEN(ins)+2)    :: outs ! add ',' around input string
INTEGER :: i, j

outs = ','
j = 2
DO i=1, LEN_TRIM(ins)
   IF (ins(i:i) == ' ') CYCLE                             ! remove spaces
   IF (ins(i:i) == ',' .AND. outs(j-1:j-1) == ',') CYCLE  ! unique ','
   outs(j:j) = ins(i:i)
   j = j+1
END DO
IF (outs(j-1:j-1) /= ',') outs(j:j) = ','
IF (outs == ',') outs = ''

END FUNCTION CLEANSTR

END SUBROUTINE CH_INIT_PRODLOSSTOT_n


