!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!Correction :
!  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-----------------------------------------------------------------
MODULE modd_repro_sum
  TYPE DOUBLE_DOUBLE
     SEQUENCE
     REAL :: R , E
  END TYPE DOUBLE_DOUBLE
END MODULE modd_repro_sum

MODULE mode_repro_sum

  USE MODD_MPIF 
  USE modd_repro_sum

  IMPLICIT NONE

  LOGICAL :: FIRST_CALL_DD     = .TRUE.
  INTEGER :: MNH_DOUBLE_DOUBLE 
  INTEGER :: MNH_SUM_DD  

!!$  INTERFACE ADD   
!!$     MODULE PROCEDURE RPDD 
!!$  END INTERFACE
!!$
!!$  INTERFACE OPERATOR(+)
!!$     MODULE PROCEDURE ADD_DD_R
!!$  END INTERFACE
!!$
!!$  INTERFACE OPERATOR(-)
!!$     MODULE PROCEDURE SUB_DD_R
!!$  END INTERFACE
!!$
!!$  INTERFACE ASSIGNMENT(=)
!!$     MODULE PROCEDURE DD_FROM_R , R_FROM_DD
!!$  END INTERFACE

CONTAINS

  SUBROUTINE INIT_DD(KINFO)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: KINFO ! MPI return status
    REAL                 :: REAL_DEFAULT
    INTEGER,PARAMETER    :: REAL_KIND=KIND(REAL_DEFAULT)
    INTEGER              :: MNH_MPI_REAL

    !
    ! find the default type of REAL for MESONH on MPI
    !
     IF (REAL_KIND .EQ. 4 ) THEN
        MNH_MPI_REAL = MPI_REAL4
     ELSE
        MNH_MPI_REAL = MPI_REAL8
    END IF
    !
    ! define the double-double for MPI
    !
    CALL MPI_TYPE_CONTIGUOUS(2, MNH_MPI_REAL ,MNH_DOUBLE_DOUBLE , KINFO)
    CALL MPI_TYPE_COMMIT(MNH_DOUBLE_DOUBLE , KINFO)
    !
    ! define the double-double sum = MNH_SUM_DD  for MPI 
    !
    CALL MPI_OP_CREATE(DDPDD, .TRUE., MNH_SUM_DD, KINFO)
    FIRST_CALL_DD = .FALSE.
    !
  END SUBROUTINE INIT_DD
  
  PURE SUBROUTINE DDPDD (dda, ddb, len, itype)
    !----------------------------------------------------------------------
    ! 
    ! Purpose: 
    ! Modification of original codes written by David H. Bailey    
    ! This subroutine computes ddb(i) = dda(i)+ddb(i)
    ! for use with MPI_*_REDUCE
    ! 
    !----------------------------------------------------------------------
    !
    ! Arguments
    !
    INTEGER, INTENT(in)                :: len       ! array length
    TYPE(DOUBLE_DOUBLE), INTENT(in)    :: dda(len)  ! input
    TYPE(DOUBLE_DOUBLE), INTENT(inout) :: ddb(len)  ! result
    INTEGER, INTENT(in)                :: itype     ! unused
    !
    ! Local workspace
    !
    REAL e, t1, t2
    INTEGER i
    !
    !-----------------------------------------------------------------------
    !
    DO i = 1, len
       !
       !   Compute dda + ddb using Knuth's trick.
       !
       t1 = dda(i)%R + ddb(i)%R
       e  = t1 - dda(i)%R
       t2 = ((ddb(i)%R - e) + (dda(i)%R - (t1 - e))) &
            + dda(i)%E + ddb(i)%E
       !
       !   The result is t1 + t2, after normalization.
       !
       ddb(i)%R = t1 + t2
       ddb(i)%E = t2 - ((t1 + t2) - t1) 
    ENDDO

    RETURN

  END SUBROUTINE DDPDD

  ELEMENTAL SUBROUTINE RPDD (a, ddb)
    !----------------------------------------------------------------------
    ! 
    ! Purpose: 
    ! Modification of original codes written by David H. Bailey    
    ! This subroutine computes ddb = a +ddb
    ! Could be inlined by compiler <=> elemental function
    ! 
    !----------------------------------------------------------------------
    !
    ! Arguments
    !
    REAL   , INTENT(in)                :: a     ! input
    TYPE(DOUBLE_DOUBLE), INTENT(inout) :: ddb   ! result
    !
    ! Local workspace
    !
    REAL e, t1, t2
    !
    !-----------------------------------------------------------------------
    !
    !   Compute dda + ddb using Knuth's trick.
    t1 = a  + ddb%R
    e  = t1 - a
    t2 = ((ddb%R - e) + (a - (t1 - e))) &
         + ddb%E
    !
    !   The result is t1 + t2, after normalization.
    ddb%R = t1 + t2
    ddb%E = t2 - ((t1 + t2) - t1) 

    RETURN

  END SUBROUTINE RPDD

  FUNCTION SUM_DD_DD1 (dda) RESULT (ddb)
    !----------------------------------------------------------------------
    ! 
    ! Purpose: 
    ! Modification of original codes written by David H. Bailey    
    ! This subroutine computes ddb = sum (dda(:))
    ! 
    !----------------------------------------------------------------------
    !
    ! Arguments
    !
    TYPE(DOUBLE_DOUBLE), DIMENSION(:) , INTENT(in)    :: dda   ! input
    TYPE(DOUBLE_DOUBLE)                               :: ddb   ! result
    !
    ! Local workspace
    !
    !REAL e, t1, t2
    REAL               ,DIMENSION(SIZE(dda,1))      :: e, t1, t2
    TYPE(DOUBLE_DOUBLE),ALLOCATABLE, DIMENSION(:)   :: ddc,ddd
    INTEGER i,j
    INTEGER  nx , lnxb2 , nxb2 , ipas, ipasm
    !
    !-----------------------------------------------------------------------
    !
    nx    =  SIZE(dda,1)
    lnxb2 =  nint(0.5+log(1.0*nx)/log(2.0))
    nxb2  =  2**lnxb2
    ALLOCATE(ddc(nxb2),ddd(nxb2))
    !print*,"nx,nxb2",nx,nxb2

    ddb%R = 0.0
    ddb%E = 0.0

!!$    ddc(1:nx)%R = dda(1:nx)%R
!!$    ddc(1:nx)%E = dda(1:nx)%E
!!$    !
!!$    ! copie directly the contribution corresponding
!!$    ! to ddd(i) = dda(i) + 0 when not a power of 2 size
!!$    !
!!$    ipas  = 2**(lnxb2-1)
!!$    ipasm = nx - ipas
!!$    ddd(ipasm+1:ipas)%R = dda(ipasm+1:ipas)%R
!!$    ddd(ipasm+1:ipas)%E = dda(ipasm+1:ipas)%E
!!$    DO j=lnxb2-1,0,-1
!!$       ! 
!!$       ! test for nx not power of 2 
!!$       ! 
!!$       ipas  = 2**j
!!$       ipasm = min(ipas,nx-ipas) 
!!$
!!$       DO i = 1, ipasm
!!$          !
!!$          !   Compute dda + ddb using Knuth's trick.
!!$          !
!!$          t1(i) = ddc(i)%R + ddc(i+ipas)%R
!!$          e(i)  = t1(i) - ddc(i)%R
!!$          t2(i) = ((ddc(i+ipas)%R - e(i)) + (ddc(i)%R - (t1(i) - e(i)))) &
!!$               + ddc(i)%E + ddc(i+ipas)%E
!!$          !
!!$          !   The result is t1 + t2, after normalization.
!!$          !
!!$          ddd(i)%R = t1(i) + t2(i)
!!$          ddd(i)%E = t2(i) - ((t1(i) + t2(i)) - t1(i)) 
!!$       ENDDO
!!$       ddc(1:ipas)%R = ddd(1:ipas)%R
!!$       ddc(1:ipas)%E = ddd(1:ipas)%E
!!$    END DO
!!$    ddb = ddc(1)

    DO i = 1, SIZE(dda,1)
       !
       !   Compute dda + ddb using Knuth's trick.
       !
       t1(i) = dda(i)%R + ddb%R
       e(i)  = t1(i) - dda(i)%R
       t2(i) = ((ddb%R - e(i)) + (dda(i)%R - (t1(i) - e(i)))) &
            + dda(i)%E + ddb%E
       !
       !   The result is t1 + t2, after normalization.
       !
       ddb%R = t1(i) + t2(i)
       ddb%E = t2(i) - ((t1(i) + t2(i)) - t1(i)) 
    ENDDO

  END FUNCTION SUM_DD_DD1

  FUNCTION SUM_DD_R2_ll (a) RESULT(c)
    !----------------------------------------------------------------------
    ! 
    ! Purpose: 
    ! Modification of original codes written by David H. Bailey    
    ! This subroutine computes ddc = ddb + a
    ! Could be inlined by compiler <=> elemental function
    ! 
    !----------------------------------------------------------------------
    USE MODE_ll , ONLY : REDUCESUM_ll
    !
    ! Arguments
    !
    REAL                                     :: c     ! result
    REAL,DIMENSION(:,:), INTENT(in)          :: a     ! input
    !
    ! Local workspace
    !
    TYPE(DOUBLE_DOUBLE)                      :: ddc  
    TYPE(DOUBLE_DOUBLE),DIMENSION(SIZE(a,1)) :: ddb   
    REAL               ,DIMENSION(SIZE(a,1)) :: e, t1, t2
    INTEGER                                  :: i,j
    INTEGER                                  :: IINFO_ll 
    !
    !-----------------------------------------------------------------------
    !
    !   Compute dda + ddb using Knuth's trick.
    ddb%R = 0.0
    ddb%E = 0.0
    DO j=1,SIZE(a,2)
       DO i=1,SIZE(a,1)
          t1(i) = a(i,j)  + ddb(i)%R
          e(i)  = t1(i) - a(i,j)
          t2(i) = ((ddb(i)%R - e(i)) + (a(i,j) - (t1(i) - e(i)))) &
               + ddb(i)%E
          !
          !   The result is t1 + t2, after normalization.
          ddb(i)%R = t1(i) + t2(i)
          ddb(i)%E = t2(i) - ((t1(i) + t2(i)) - t1(i)) 
       END DO
    END DO
    ddc = SUM_DD_DD1(ddb)
    CALL REDUCESUM_ll(ddc,IINFO_ll)
    c = ddc%R
  END FUNCTION SUM_DD_R2_LL

 FUNCTION SUM_DD_R2_R1_ll (a) RESULT(c)
    !----------------------------------------------------------------------
    ! 
    ! Purpose: 
    ! Modification of original codes written by David H. Bailey    
    ! This subroutine computes c(1:n2) = sum_dd(a(:,1:n2)) on all processors
    ! 
    !----------------------------------------------------------------------
    USE MODE_ll , ONLY : REDUCESUM_ll
    !
    ! Arguments
    !
    REAL,DIMENSION      (:,:), INTENT(in)          :: a     ! input
    REAL,DIMENSION(SIZE(a,2))                      :: c     ! result
    !
    ! Local workspace
    !
    TYPE(DOUBLE_DOUBLE),DIMENSION(SIZE(a,2)) :: ddb   
    REAL               ,DIMENSION(SIZE(a,2)) :: e, t1, t2
    INTEGER                                  :: i,j
    INTEGER                                  :: IINFO_ll 
    !
    !-----------------------------------------------------------------------
    !
    !   Compute dda + ddb using Knuth's trick.
    ddb(:)%R = 0.0
    ddb(:)%E = 0.0
    DO j=1,SIZE(a,2)
       DO i=1,SIZE(a,1)
          t1(j) = a(i,j)  + ddb(j)%R
          e(j)  = t1(j) - a(i,j)
          t2(j) = ((ddb(j)%R - e(j)) + (a(i,j) - (t1(j) - e(j)))) &
               + ddb(j)%E
          !
          !   The result is t1 + t2, after normalization.
          ddb(j)%R = t1(j) + t2(j)
          ddb(j)%E = t2(j) - ((t1(j) + t2(j)) - t1(j)) 
       END DO
    END DO
    CALL REDUCESUM_ll(ddb,IINFO_ll)
    c(:) = ddb(:)%R
  END FUNCTION SUM_DD_R2_R1_ll

 FUNCTION SUM_DD_R1_ll (a) RESULT(c)
    !----------------------------------------------------------------------
    ! 
    ! Purpose: 
    ! Modification of original codes written by David H. Bailey    
    ! This subroutine computes c = sum_dd(a(:)) on all processors
    ! 
    !----------------------------------------------------------------------
    USE MODE_ll , ONLY : REDUCESUM_ll
    !
    ! Arguments
    !
    REAL,DIMENSION        (:), INTENT(in)          :: a     ! input
    REAL                                           :: c     ! result
    !
    ! Local workspace
    !
    TYPE(DOUBLE_DOUBLE)                      :: ddb   
    REAL                                     :: e, t1, t2
    INTEGER                                  :: i,j
    INTEGER                                  :: IINFO_ll 
    !
    !-----------------------------------------------------------------------
    !
    !   Compute dda + ddb using Knuth's trick.
    ddb%R = 0.0
    ddb%E = 0.0
    DO i=1,SIZE(a,1)
       t1 = a(i)  + ddb%R
       e  = t1 - a(i)
       t2 = ((ddb%R - e) + (a(i) - (t1 - e))) &
            + ddb%E
       !
       !   The result is t1 + t2, after normalization.
       ddb%R = t1 + t2
       ddb%E = t2 - ((t1 + t2) - t1) 
    END DO
    CALL REDUCESUM_ll(ddb,IINFO_ll)
    c = ddb%R
  END FUNCTION SUM_DD_R1_ll

!!$  ELEMENTAL FUNCTION ADD_DD_R (ddb, a) RESULT(ddc)
!!$    !----------------------------------------------------------------------
!!$    ! 
!!$    ! Purpose: 
!!$    ! Modification of original codes written by David H. Bailey    
!!$    ! This subroutine computes ddc = ddb + a
!!$    ! Could be inlined by compiler <=> elemental function
!!$    ! 
!!$    !----------------------------------------------------------------------
!!$    !
!!$    ! Arguments
!!$    !
!!$    TYPE(DOUBLE_DOUBLE)                :: ddc   ! result
!!$    TYPE(DOUBLE_DOUBLE), INTENT(in)    :: ddb   ! input
!!$    REAL               , INTENT(in)    :: a     ! input
!!$    !
!!$    ! Local workspace
!!$    !
!!$    REAL e, t1, t2
!!$    !
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    !   Compute dda + ddb using Knuth's trick.
!!$    t1 = a  + ddb%R
!!$    e  = t1 - a
!!$    t2 = ((ddb%R - e) + (a - (t1 - e))) &
!!$         + ddb%E
!!$    !
!!$    !   The result is t1 + t2, after normalization.
!!$    ddc%R = t1 + t2
!!$    ddc%E = t2 - ((t1 + t2) - t1) 
!!$
!!$  END FUNCTION ADD_DD_R

!!$  FUNCTION SUM_DD_R1 (a) RESULT(ddc)
!!$    !----------------------------------------------------------------------
!!$    ! 
!!$    ! Purpose: 
!!$    ! Modification of original codes written by David H. Bailey    
!!$    ! This subroutine computes ddc = sum( a(:))
!!$    ! 
!!$    !----------------------------------------------------------------------
!!$    !
!!$    ! Arguments
!!$    !
!!$    TYPE(DOUBLE_DOUBLE)                :: ddc   ! result
!!$    REAL,DIMENSION(:)  , INTENT(in)    :: a     ! input
!!$    !
!!$    ! Local workspace
!!$    !
!!$    REAL e, t1, t2
!!$    INTEGER i
!!$    !
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    !   Compute dda + ddb using Knuth's trick.
!!$    ddc%R = 0.0
!!$    ddc%E = 0.0
!!$    DO i=1,SIZE(a)
!!$       t1 = a(i)  + ddc%R
!!$       e  = t1 - a(i)
!!$       t2 = ((ddc%R - e) + (a(i) - (t1 - e))) &
!!$            + ddc%E
!!$       !
!!$       !   The result is t1 + t2, after normalization.
!!$       ddc%R = t1 + t2
!!$       ddc%E = t2 - ((t1 + t2) - t1) 
!!$    END DO
!!$  END FUNCTION SUM_DD_R1
  
!!$  ELEMENTAL FUNCTION SUB_DD_R (ddb, a) RESULT(ddc)
!!$    !----------------------------------------------------------------------
!!$    ! 
!!$    ! Purpose: 
!!$    ! Modification of original codes written by David H. Bailey    
!!$    ! This subroutine computes ddc = ddb - a
!!$    ! 
!!$    !----------------------------------------------------------------------
!!$    !
!!$    ! Arguments
!!$    !
!!$    TYPE(DOUBLE_DOUBLE)                :: ddc   ! result
!!$    TYPE(DOUBLE_DOUBLE), INTENT(in)    :: ddb   ! input
!!$    REAL               , INTENT(in)    :: a     ! input
!!$    !
!!$    ! Local workspace
!!$    !
!!$    REAL e, t1, t2
!!$    !
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    !   Compute dda + ddb using Knuth's trick.
!!$    t1 = - a  + ddb%R
!!$    e  = t1 + a
!!$    t2 = ((ddb%R - e) + ( - a - (t1 - e))) &
!!$         + ddb%E
!!$    !
!!$    !   The result is t1 + t2, after normalization.
!!$    ddc%R = t1 + t2
!!$    ddc%E = t2 - ((t1 + t2) - t1) 
!!$
!!$  END FUNCTION SUB_DD_R

!!$  ELEMENTAL SUBROUTINE DD_FROM_R (ddb, a) 
!!$    !----------------------------------------------------------------------
!!$    ! 
!!$    ! Purpose: 
!!$    ! Modification of original codes written by David H. Bailey    
!!$    ! This subroutine computes ddb = a
!!$    ! Could be inlined by compiler <=> elemental function
!!$    ! 
!!$    !----------------------------------------------------------------------
!!$    !
!!$    ! Arguments
!!$    !
!!$    TYPE(DOUBLE_DOUBLE), INTENT(out)   :: ddb   ! input
!!$    REAL               , INTENT(in)    :: a     ! input
!!$    !
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    ddb%R = a
!!$    ddB%E = 0.0 
!!$    
!!$  END SUBROUTINE DD_FROM_R

!!$  ELEMENTAL SUBROUTINE R_FROM_DD (a, ddb) 
!!$    !----------------------------------------------------------------------
!!$    ! 
!!$    ! Purpose: 
!!$    ! Modification of original codes written by David H. Bailey    
!!$    ! This subroutine computes a = ddb
!!$    ! Could be inlined by compiler <=> elemental function
!!$    ! 
!!$    !----------------------------------------------------------------------
!!$    !
!!$    ! Arguments
!!$    !
!!$    REAL               , INTENT(out)    :: a     ! input
!!$    TYPE(DOUBLE_DOUBLE), INTENT(in)     :: ddb   ! input
!!$    !
!!$    !-----------------------------------------------------------------------
!!$    !
!!$    a =  ddb%R
!!$    
!!$  END SUBROUTINE R_FROM_DD

END MODULE mode_repro_sum
