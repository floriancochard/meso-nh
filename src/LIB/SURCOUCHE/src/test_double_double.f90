!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
PROGRAM TEST_DOUBLE_DOUBLE
  ! This code calculates the summation of an array of real numbers
  ! distributed on multiple processors using double-double precision.

  USE MODD_MPIF
  USE mode_repro_sum

  IMPLICIT NONE

  !INCLUDE 'mpif.h'
  INTEGER , PARAMETER :: n =1024
  REAL array(n)
  INTEGER myPE, totPEs, stat(MPI_STATUS_SIZE), ierr
  INTEGER MPI_SUMDD2, itype ,i 
  
  TYPE(DOUBLE_DOUBLE) local_sum, global_sum
  CALL MPI_INIT (ierr)
  CALL MPI_COMM_RANK (MPI_COMM_WORLD, myPE, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, totPEs, ierr)

  array = (myPE+1)*n

  ! operator MPI_SUMDD is created based on an external function CPC.
  !CALL MPI_OP_CREATE (CPC, .TRUE., MPI_SUMDD2, ierr)
  CALL MPI_OP_CREATE (DDPDD, .TRUE., MPI_SUMDD2, ierr)

  ! assume array(n) is the local part of a global distributed array.

  local_sum = DOUBLE_DOUBLE ( 0.0 , 0.0 )

  DO i = 1, n
     CALL RPDD (array(i) , local_sum )
     !local_sum = local_sum + array(i)
  ENDDO
  ! add all local_sums on each PE to PE0 with MPI_SUMDD.
  ! global_sum is a complex number, represents final (sum, error).

 global_sum  =  DOUBLE_DOUBLE ( 0.0 , 0.0 )

  CALL MPI_REDUCE (local_sum, global_sum, 1, MNH_DOUBLE_DOUBLE, MPI_SUMDD2, &
       & 0, MPI_COMM_WORLD, ierr)

  if ( myPE .eq. 0 ) print*,"global_sum =" , global_sum
  CALL MPI_FINALIZE(ierr)

END PROGRAM TEST_DOUBLE_DOUBLE
