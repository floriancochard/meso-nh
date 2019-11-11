!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##################
      MODULE MODI_SUM_ll
!     ##################
!
INTERFACE
!
!!     #######################################################
       FUNCTION EXTRACT_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                            KXEND, KYEND, KZEND )
!!     #######################################################
!
  REAL, DIMENSION(:,:,:), POINTER :: EXTRACT_ll ! Result
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
       END FUNCTION EXTRACT_ll
!
!!     ###########################################################
       FUNCTION SUM1D_ll( PFIELD, KDIR, KINFO, KXOR, KYOR, KZOR, &
                           KXEND, KYEND, KZEND )
!!     ###########################################################
!
  REAL, DIMENSION(:,:), POINTER :: SUM1D_ll ! Result
!
  INTEGER, INTENT(IN) :: KDIR ! Summation direction (1, 2 or 3)
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
       END FUNCTION SUM1D_ll
!
!!     ###################################################################
       FUNCTION SUM2D_ll( PFIELD, KDIR1, KDIR2, KINFO, KXOR, KYOR, KZOR, &
                          KXEND, KYEND, KZEND )
!!     ###################################################################
!
  REAL, DIMENSION(:), POINTER :: SUM2D_ll ! Result
!
  INTEGER, INTENT(IN) :: KDIR1, KDIR2 ! Summation directions (1, 2 or 3)
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
       END FUNCTION SUM2D_ll
!
!!     #####################################################
       FUNCTION SUM3D_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                          KXEND, KYEND, KZEND )
!!     #####################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
       END FUNCTION SUM3D_ll
!
!!     #######################################################################
       FUNCTION SUM_1DFIELD_ll( PFIELD, HDIR, KOR, KEND, KERR ) RESULT( ZSUM )
!!     #######################################################################
!
  REAL, DIMENSION(:), INTENT(IN) :: PFIELD ! 1d Field
!
  CHARACTER(LEN=1) :: HDIR               ! direction of the 1D field
!
  INTEGER, OPTIONAL, INTENT(OUT) :: KERR ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KOR, KEND ! Coordinates of the region
!
  REAL :: ZSUM ! result
!
       END FUNCTION SUM_1DFIELD_ll
!
!!     ########################################################
       REAL FUNCTION MAX_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                             KXEND, KYEND, KZEND )
!!     ########################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
       END FUNCTION MAX_ll
!
!!     ########################################################
       REAL FUNCTION MIN_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                             KXEND, KYEND, KZEND )
!!     ########################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
       END FUNCTION MIN_ll
!
!!     ###########################################
       FUNCTION SUMMASK_ll( PFIELD, OMASK, KINFO )
!!     ###########################################
!
  REAL, DIMENSION(:), POINTER :: SUMMASK_ll ! Result
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASK ! 2d Mask
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
       END FUNCTION SUMMASK_ll
!
!!     ####################################################
       REAL FUNCTION SUMMASKCOMP_ll( PFIELD, OMASK, KINFO )
!!     ####################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASK ! 2d Mask
! 
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
       END FUNCTION SUMMASKCOMP_ll
!
!!     #############################################
       SUBROUTINE SUM_DIM1_ll( PFIELD, PRES, KINFO )
!!     #############################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  REAL, DIMENSION(:, :), INTENT(OUT) :: PRES
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE SUM_DIM1_ll
!
!!     #############################################
       SUBROUTINE SUM_DIM2_ll( PFIELD, PRES, KINFO )
!!     #############################################
!
  REAL, DIMENSION(:), INTENT(IN) :: PFIELD
  REAL, DIMENSION(:), INTENT(OUT) :: PRES
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE SUM_DIM2_ll
!
END INTERFACE
!
INTERFACE GMAXLOC_ll
!
!!     #######################################################
       FUNCTION GMAXLOC3D_ll( PARRAY, MASK ) RESULT( KMAXLOC )
!!     #######################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN)              :: PARRAY   ! input array in
                                    ! which the maximum is to be found
!
  LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: MASK     ! mask 
  INTEGER, DIMENSION(3)                           :: KMAXLOC  ! indices
          ! of the maximum value on the whole domain (global coordinates)
!
       END FUNCTION GMAXLOC3D_ll
!
!!     ##############################################################
       FUNCTION GMAXLOC2D_ll( PARRAY, KDIMS, MASK ) RESULT( KMAXLOC )
!!     ##############################################################
!
  REAL, DIMENSION(:,:), INTENT(IN)              :: PARRAY   ! input array in
                                   ! which the maximum is to be found
  LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK     ! mask 
  INTEGER, DIMENSION(2)                         :: KMAXLOC  ! indices
          ! of the maximum value on the whole domain (global coordinates)
  INTEGER, DIMENSION(2), OPTIONAL               :: KDIMS
!
       END FUNCTION GMAXLOC2D_ll
!
!!     ##############################################################
       FUNCTION GMAXLOC1D_ll( PARRAY, KDIMS, MASK ) RESULT( KMAXLOC )
!!     ##############################################################
!
  REAL, DIMENSION(:), INTENT(IN)              :: PARRAY   ! input array in
                                 ! which the maximum is to be found
  LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: MASK     ! mask 
  INTEGER                                     :: KMAXLOC  ! indice
          ! of the maximum value on the whole domain (global coordinates)
  INTEGER,  OPTIONAL                          :: KDIMS
!
       END FUNCTION GMAXLOC1D_ll
!
END INTERFACE
!
INTERFACE GMINLOC_ll
!
!!     #######################################################
       FUNCTION GMINLOC3D_ll( PARRAY, MASK ) RESULT( KMINLOC )
!!     #######################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN)              :: PARRAY   ! input array in
                                 ! which the minimum is to be found
  LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: MASK     ! mask 
  INTEGER, DIMENSION(3)                           :: KMINLOC  ! indices
           ! of the minimum value on the whole domain (global coordinates)
!
       END FUNCTION GMINLOC3D_ll
!
!!     ##############################################################
       FUNCTION GMINLOC2D_ll( PARRAY, KDIMS, MASK ) RESULT (KMINLOC )
!!     ##############################################################
!
  REAL, DIMENSION(:,:), INTENT(IN)              :: PARRAY   ! input array in
                                 ! which the minimum is to be found
  LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK     ! mask 
  INTEGER, DIMENSION(2)                         :: KMINLOC  ! indices
           ! of the minimum value on the whole domain (global coordinates)
  INTEGER, DIMENSION(2), OPTIONAL               :: KDIMS
!
       END FUNCTION GMINLOC2D_ll
!
!!     ##############################################################
       FUNCTION GMINLOC1D_ll( PARRAY, KDIMS, MASK ) RESULT( KMINLOC )
!!     ##############################################################
!
  REAL, DIMENSION(:), INTENT(IN)              :: PARRAY   ! input array in
                                 ! which the minimum is to be found
  LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: MASK     ! mask 
  INTEGER                                     :: KMINLOC  ! indice
           ! of the minimum value on the whole domain (global coordinates)
  INTEGER,  OPTIONAL                       :: KDIMS
!
       END FUNCTION GMINLOC1D_ll
!
END INTERFACE
!
INTERFACE REDUCESUM_ll
!
!!     ###########################################
       SUBROUTINE REDUCE_SUM_0DD_ll( PRES, KINFO )
!!     ###########################################
!
  USE MODD_REPRO_SUM  , ONLY          :  DOUBLE_DOUBLE
!
  TYPE(DOUBLE_DOUBLE) , INTENT(INOUT) :: PRES ! sum
!
  INTEGER             , INTENT(OUT)   :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_0DD_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_0D_ll( PRES, KINFO )
!!     ##########################################
!
  REAL, INTENT(INOUT) :: PRES ! sum 
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_0D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_1DD_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODD_REPRO_SUM  , ONLY                     :  DOUBLE_DOUBLE
!
  TYPE(DOUBLE_DOUBLE),DIMENSION(:),INTENT(INOUT) :: PRES ! sum
!
  INTEGER                   , INTENT(OUT)  :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_1DD_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_1D_ll( PRES, KINFO )
!!     ##########################################
!
  REAL,DIMENSION(:),INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_1D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_2D_ll( PRES, KINFO )
!!     ##########################################
!
  REAL,DIMENSION(:,:),INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_2D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_3D_ll( PRES, KINFO )
!!     ##########################################
!
  REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_3D_ll
!
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I0D_ll( PRES, KINFO )
!!     ##########################################
!
  INTEGER, INTENT(INOUT) :: PRES ! sum 
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_I0D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I1D_ll( PRES, KINFO )
!!     ##########################################
!
  INTEGER,DIMENSION(:),INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_I1D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I2D_ll( PRES, KINFO )
!!     ##########################################
!
  INTEGER,DIMENSION(:,:),INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_I2D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I3D_ll( PRES, KINFO )
!!     ##########################################
!
  INTEGER, DIMENSION(:,:,:),INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
       END SUBROUTINE REDUCE_SUM_I3D_ll
!
END INTERFACE
!
END MODULE MODI_SUM_ll
