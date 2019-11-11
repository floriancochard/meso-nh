!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_MNH_TIMING
!
! Modification :
! J.ESCOBAR 13/11/2008 : change (2) in (:) for bug in IBM-SP6 compiler
! J.Escobar 1/09/2011  : reduce 'timing' format
! J.Escobar 12/02/2013 : tribulle to slow on large BG partition , inhib it by a early return in the code 
! Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!

INTEGER     :: NLUOUT_TIMING

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SECOND_MNH2(XT)
!
USE modd_mpif
!
REAL*8,DIMENSION(2)           :: XT
!
CALL CPU_TIME(XT(1))
XT(2) = MPI_Wtime()
END SUBROUTINE SECOND_MNH2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!JUAN
      SUBROUTINE SET_ILUOUT_TIMING(KLUOUT)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: KLUOUT
        NLUOUT_TIMING = KLUOUT
      END SUBROUTINE SET_ILUOUT_TIMING

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE TIMING_SEPARATOR(HSEP)
        IMPLICIT NONE
        CHARACTER :: HSEP
        INTEGER   :: J  
        WRITE(NLUOUT_TIMING,FMT= "('|',100(A1),'|')" ) ( HSEP , J=1,100 )
      END SUBROUTINE TIMING_SEPARATOR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE TIMING_LEGEND()
        CALL  TIMING_SEPARATOR('-')
        WRITE(NLUOUT_TIMING,FMT="( '|     CPUTIM/ELAPSE                     |&
         &|   SUM(PROC)   |MEAN(PROC)| MIN(PROC | MAX(PROC)| PERCENT %|')" ) 
        CALL  TIMING_SEPARATOR('-')
      END SUBROUTINE TIMING_LEGEND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ########################################
      SUBROUTINE TIME_HEADER_ll(KMI)
!     ########################################
        IMPLICIT NONE 
        INTEGER :: KMI
        CALL  TIMING_SEPARATOR('-')
        CALL  TIMING_SEPARATOR(' ')
        WRITE(NLUOUT_TIMING,FMT= "('|',32X,'COMPUTING TIME ANALYSIS in MODEL',I0,35X,'|')" ) KMI
        CALL  TIMING_SEPARATOR(' ')
        CALL  TIMING_LEGEND()  
      END SUBROUTINE TIME_HEADER_ll


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ########################################
      SUBROUTINE TIME_STAT_ll(PRES, PSUM, HPRINT, HSEP,HFULL)
!     ########################################
!
!*       0.    DECLARATIONS
!
  USE MODD_MPIF
  USE MODD_VAR_ll, ONLY : MPI_PRECISION, NPROC, IP
  !JUANZ
  USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
  !JUANZ
  !
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
! 
  REAL*8,DIMENSION(:), INTENT(IN)         :: PRES ! (1)=CPU & (2)=ELAPSE Proccessors Timing
!
  REAL*8,DIMENSION(:), INTENT(INOUT)      :: PSUM ! (1)=SUM(CPU) & (2)=SUM(ELAPSE) Timing
!
  CHARACTER(len=*),  INTENT(IN),OPTIONAL :: HPRINT
  CHARACTER       ,  INTENT(IN),OPTIONAL :: HSEP
  CHARACTER(len=*),  INTENT(IN),OPTIONAL :: HFULL
!
!*       0.2   Declarations of local variables :
!
  INTEGER,PARAMETER         :: NSTAT=5
  REAL*8,DIMENSION(2,NSTAT) :: ZSTAT ! (1)=Sum(proc),(2)=Sum/Nproc,(3)=Min(proc),(4)=Max(proc),(5)=Purcent(1)
  INTEGER                   :: INFO,IROOT,JP
  CHARACTER(len=30)         :: VIDE = ""
  CHARACTER(len=30)        :: FILE = ""
  INTEGER                   :: IC

  REAL*8,DIMENSION(2,NPROC)   :: ZSTAT_ALL 
  INTEGER, DIMENSION(NPROC) :: IND
  INTEGER :: ILU
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
INFO = -1
! 1.1 Sum(Proc)
  CALL MPI_ALLREDUCE(PRES, ZSTAT(:,1), 2, MPI_REAL8, &
                     MPI_SUM, NMNH_COMM_WORLD, INFO)
! 1.2 Sum/Proc
  ZSTAT(:,2) = ZSTAT(:,1 ) / NPROC
! 1.3 Min(Proc)
  CALL MPI_ALLREDUCE(PRES, ZSTAT(:,3), 2, MPI_REAL8, &
                     MPI_MIN, NMNH_COMM_WORLD, INFO)
! 1.4 Max(Proc)
  CALL MPI_ALLREDUCE(PRES, ZSTAT(:,4), 2, MPI_REAL8, &
                     MPI_MAX, NMNH_COMM_WORLD, INFO)


  IF (.NOT.PRESENT(HPRINT)) THEN
   ! return total sum 
   PSUM = ZSTAT(:,1)
  !

  ELSEIF ( ZSTAT(1,1) > 0.0 ) THEN
   ! use Psum , for print stat & pourcent
   ! Purcent
     WHERE ( PSUM /= 0.0 )
        ZSTAT(:,5) = 100.0 * ZSTAT(:,1) / PSUM(:)
     ELSEWHERE
        ZSTAT(:,5) = 0.0
     END WHERE
   ! print stat
   !
   IF (PRESENT(HSEP)) CALL  TIMING_SEPARATOR(HSEP)
   WRITE(NLUOUT_TIMING,FMT= "('|',A30,'| CPUTIM ||',F15.3,'|',4(F10.3,'|'),F7.3,'|')" ) HPRINT//VIDE,ZSTAT(1,:)
   WRITE(NLUOUT_TIMING,FMT= "('|',A30,'| ELAPSE ||',F15.3,'|',4(F10.3,'|'),F7.3,'|')" ) HPRINT//VIDE,ZSTAT(2,:)

   IF (PRESENT(HFULL)) THEN
      ! gather all data
      !CALL  TIMING_SEPARATOR(HSEP)
      IROOT = 0
      CALL MPI_GATHER(PRES(:),2,MPI_REAL8,ZSTAT_ALL(:,1),2,MPI_REAL8,&
           IROOT,NMNH_COMM_WORLD, INFO)
      IF (IP.EQ.1) THEN
         FILE = trim(adjustl(HPRINT))
         DO IC=1,LEN_trim(FILE)
            SELECT CASE( ICHAR(FILE(IC:IC)) )
               CASE(ICHAR(' '),ICHAR('='),ICHAR('>'))
                 FILE(IC:IC) = '_'
            END SELECT

         END DO
         OPEN (newunit=ILU,file=FILE)
         WRITE(unit=ILU,FMT= "(10('#'),A30,10('#'))" ) HPRINT//VIDE
         call  triabulle(ZSTAT_ALL(2,:),ind)
         DO JP=1,NPROC
            WRITE(unit=ILU,FMT= "(5(I8,' ',F15.3,' '))" ) &
                 ind(JP),ZSTAT_ALL(2,ind(JP))
         END DO
         CLOSE (unit=ILU)
      END IF
   END IF
   !
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE TIME_STAT_ll
!JUAN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine triabulle(vec,ind)
implicit none

real*8   , intent(inout)   :: vec(:)
integer, intent(out)  :: ind(:)


integer               :: n 

logical :: a
integer :: i
integer :: mem

n = size(vec)
a = .true.
do i=1,n 
   ind(i) = i 
enddo

return
!JUAN TO SLOW ON BG !!!

do while (a)
   a=.false.
   do i=1,n-1
      if (vec(ind(i))>vec(ind(i+1))) then
         mem = ind(i)
         ind(i) = ind(i+1)
         ind(i+1) = mem
         a=.true.
      endif
   enddo
enddo

end subroutine triabulle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END MODULE MODE_MNH_TIMING
