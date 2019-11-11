MODULE MODE_MNH_TIMING
!
! Modification :
! J.ESCOBAR 13/11/2008 : change (2) in (:) for bug in IBM-SP6 compiler
!

INTEGER     :: NLUOUT_TIMING

CONTAINS

SUBROUTINE SECOND_MNH2(XT)
!
USE modd_mpif
!
REAL,DIMENSION(2)           :: XT
!
CALL CPU_TIME(XT(1))
XT(2) = MPI_Wtime()
!ZT = MPI_WTICK()
!print*,"MPI_WTICK()=",ZT,XT
END SUBROUTINE SECOND_MNH2

!JUAN
      SUBROUTINE SET_ILUOUT_TIMING(KLUOUT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: KLUOUT
        NLUOUT_TIMING = KLUOUT
      END SUBROUTINE SET_ILUOUT_TIMING

      SUBROUTINE TIMING_SEPARATOR(HSEP)
        IMPLICIT NONE
        CHARACTER :: HSEP
        INTEGER   :: J  
        WRITE(NLUOUT_TIMING,FMT= "('|',120(A1),'|')" ) ( HSEP , J=1,120 )
      END SUBROUTINE TIMING_SEPARATOR
      !
      SUBROUTINE TIMING_LEGEND()
        CALL  TIMING_SEPARATOR('-')
        WRITE(NLUOUT_TIMING,FMT="( '|     CPUTIM/ELAPSE                     |&
         &|   SUM(PROC)   |  MEAN(PROC)   |  MIN(PROC)    |  MAX(PROC)    | PERCENTAGE %  |')" ) 
        CALL  TIMING_SEPARATOR('-')
      END SUBROUTINE TIMING_LEGEND
!     ########################################
      SUBROUTINE TIME_HEADER_ll(KMI)
!     ########################################
        IMPLICIT NONE 
        INTEGER :: KMI
        CALL  TIMING_SEPARATOR('-')
        CALL  TIMING_SEPARATOR(' ')
        WRITE(NLUOUT_TIMING,FMT= "('|',32X,'COMPUTING TIME ANALYSIS in MODEL',I0,55X,'|')" ) KMI
        CALL  TIMING_SEPARATOR(' ')
        CALL  TIMING_LEGEND()  
      END SUBROUTINE TIME_HEADER_ll
!     ########################################
      SUBROUTINE TIME_STAT_ll(PRES, PSUM, HPRINT, HSEP)
!     ########################################
!
!*       0.    DECLARATIONS
!
  USE MODD_MPIF
  USE MODD_VAR_ll, ONLY : MPI_PRECISION, NPROC
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
! 
  REAL,DIMENSION(:), INTENT(IN)         :: PRES ! (1)=CPU & (2)=ELAPSE Proccessors Timing
!
  REAL,DIMENSION(:), INTENT(INOUT)      :: PSUM ! (1)=SUM(CPU) & (2)=SUM(ELAPSE) Timing
!
  CHARACTER(len=*),  INTENT(IN),OPTIONAL :: HPRINT
  CHARACTER       ,  INTENT(IN),OPTIONAL :: HSEP
!
!*       0.2   Declarations of local variables :
!
  INTEGER,PARAMETER       :: NSTAT=5
  REAL,DIMENSION(2,NSTAT) :: ZSTAT ! (1)=Sum(proc),(2)=Sum/Nproc,(3)=Min(proc),(4)=Max(proc),(5)=Purcent(1)
  INTEGER                 :: INFO
  CHARACTER(len=30)       :: VIDE = ""
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
INFO = -1
! 1.1 Sum(Proc)
  CALL MPI_ALLREDUCE(PRES, ZSTAT(1,1), 2, MPI_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, INFO)
! 1.2 Sum/Proc
  ZSTAT(:,2) = ZSTAT(:,1 ) / NPROC
! 1.3 Min(Proc)
  CALL MPI_ALLREDUCE(PRES, ZSTAT(1,3), 2, MPI_PRECISION, &
                     MPI_MIN, MPI_COMM_WORLD, INFO)
! 1.4 Max(Proc)
  CALL MPI_ALLREDUCE(PRES, ZSTAT(1,4), 2, MPI_PRECISION, &
                     MPI_MAX, MPI_COMM_WORLD, INFO)
!
  IF (.NOT.PRESENT(HPRINT)) THEN
   ! return total sum 
   PSUM = ZSTAT(:,1)
  !

  ELSEIF ( ZSTAT(1,1) > 0.0 ) THEN
   ! use Psum , for print stat & pourcent
   ! Purcent
   ZSTAT(:,5) = 100.0 * ZSTAT(:,1) / PSUM(:)
   ! print stat
   !
   IF (PRESENT(HSEP)) CALL  TIMING_SEPARATOR(HSEP)
   WRITE(NLUOUT_TIMING,FMT= "('|',A30,'| CPUTIM ||',5(F15.3,'|'))" ) HPRINT//VIDE,ZSTAT(1,:)
   WRITE(NLUOUT_TIMING,FMT= "('|',A30,'| ELAPSE ||',5(F15.3,'|'))" ) HPRINT//VIDE,ZSTAT(2,:)

   !
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE TIME_STAT_ll
!JUAN
    END MODULE MODE_MNH_TIMING
