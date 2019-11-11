!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!Correction :
!  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-----------------------------------------------------------------

!     ##################
      MODULE MODE_SUM_ll
!     ##################
!
!!    Purpose
!!    -------
!
!     The purpose of this module is to provide functions for the summations
!     of the budget computations, for the search of minimum and maximum.
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!     SUBROUTINES : SUM_DIM1_ll, SUM_DIM2_ll,
!                   REDUCESUM_ll (REDUCE_SUM_0D_ll, REDUCE_SUM_1D_ll,
!                                 REDUCE_SUM_2D_ll, REDUCE_SUM_3D_ll)
!     FUNCTIONS   : EXTRACT_ll,
!                   SUM1D_ll, SUM2D_ll, SUM3D_ll, SUM_1DFIELD_ll,
!                   MAX_ll, MIN_ll,
!                   SUMMASK_ll, SUMMASKCOMP_ll
!
!!    Reference
!!    ---------
!
!!    Authors
!!    -------
!
!     R. Guivarch  * CERFACS *
!     P. Kloos (CERFACS/CNRM)
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NPROC -
!       IP -
!       JPHALO -
!       MPI_PRECISION - 
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!------------------------------------------------------------------------------
!
   USE MODD_MPIF
   !JUANZ
   USE MODD_VAR_ll, ONLY : NMNH_COMM_WORLD
   !JUANZ
!
!  INCLUDE 'mpif.h'
!
  CONTAINS
!
!     #######################################################
      FUNCTION EXTRACT_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                           KXEND, KYEND, KZEND )
!     #######################################################
!
!!****  *EXTRACT_ll* - function to extract a subset of a 3D field
!                      contained in a geographic region given by its first point
!                      and its last point. If these two points are omitted,
!                      the whole physical domain is considered.
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the intersection of its sub-domain
!     with the geographic region.
!     If this intersection is not empty, the processor fills
!     an intermediate buffer of the dimension of the result with its elements.
!     Then we make a summation of all the intermediate buffers (MPI reduction)
!     and obtain the result.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       MPI_PRECISION - 
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!
!!    Modifications
!!    -------------
!     Original 03/07/98
!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), POINTER :: EXTRACT_ll ! Result
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, &  ! Coordinates
                                    KXEND, KYEND, KZEND ! of the region
!
!*        0.2   declarations of local variables
!
  INTEGER :: IXOR, IYOR, IZOR, & ! Real coordinates of the region
             IXEND, IYEND, IZEND
!
  REAL, ALLOCATABLE, TARGET :: ZBUF(:,:,:), ZBUFD(:,:,:) ! Intermediate buffers
!
  INTEGER :: IIE, IJE, IIB, IJB ! Displacements in PFIELD
!
  INTEGER :: IIB_BUF,IJB_BUF,IIE_BUF,IJE_BUF ! Displacements in the buffers
!
  INTEGER :: IDIM ! Dimension of the result array
!
  INTEGER :: IWEST, IEAST, & ! increments for the boundaries
             ISOUTH, INORTH  !   they are set to 0 if the subdomain
                             !   is not at a boundary
                             !   they are set to -1 or 1 if the subdomain is at
                             ! the corresponding boundary
!
  INTEGER :: IIMAX, IJMAX ! global dimensions of the physical model
!
!-------------------------------------------------------------------------------
!
  KINFO = 0
!
!*        1.    INITIALISATIONS :
!               -----------------
!
!*        1.1   Fill in the increments for the boundaries
!
  IWEST = 0
  IEAST = 0
  ISOUTH = 0
  INORTH = 0
!
  IF (LWEST_ll()) IWEST = -1
  IF (LEAST_ll()) IEAST = 1 
  IF (LSOUTH_ll()) ISOUTH = -1
  IF (LNORTH_ll()) INORTH = 1
!
!*        1.2   Get the model dimensions
!
  CALL GET_GLOBALDIMS_ll(IIMAX, IJMAX)
!
!*        1.3   Fill The Coordinates :
!
  IF(PRESENT(KXOR) .AND. PRESENT(KYOR) .AND. PRESENT(KZOR) .AND. &
     PRESENT(KXEND) .AND. PRESENT(KYEND) .AND. PRESENT(KZEND) ) THEN
!
    IXOR = KXOR
    IYOR = KYOR
    IZOR = KZOR
    IXEND = KXEND
    IYEND = KYEND
    IZEND = KZEND
!
  ELSEIF(.NOT.PRESENT(KXOR) .AND. .NOT.PRESENT(KYOR) &
        .AND. .NOT.PRESENT(KZOR) .AND. .NOT.PRESENT(KXEND) &
        .AND. .NOT.PRESENT(KYEND) .AND. .NOT.PRESENT(KZEND) ) THEN
!
    IXOR = 1 + JPHEXT
    IYOR = 1 + JPHEXT
    IZOR = 1 + JPVEXT
    IXEND = IIMAX + JPHEXT
    IYEND = IJMAX + JPHEXT
    IZEND = SIZE(PFIELD,3) - JPVEXT
!
  ELSE
!
! Error
!
    ALLOCATE(ZBUF(1,1,1))
    ZBUF = 0
!
    KINFO = -1
!
    EXTRACT_ll => ZBUF
!
    DEALLOCATE(ZBUF)
    RETURN
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.    THE COORDINATES ARE NOT IN THE GLOBAL DOMAIN -> ERROR :
!               -----------------------------------------------------
!
  IF((IXOR < 1 ) .OR. (IYOR < 1 ) .OR. &
     (IXEND > IIMAX + 2*JPHEXT) .OR. (IYEND > IJMAX + 2*JPHEXT)) THEN
!
    ALLOCATE(ZBUF(1,1,1))
    ZBUF = 0
!
!   Error
!
    KINFO = -1
!
    EXTRACT_ll => ZBUF
!
    DEALLOCATE(ZBUF)
    RETURN
!
!-------------------------------------------------------------------------------
!
!*        3.    THE COORDINATES ARE IN THE GLOBAL DOMAIN :
!               ----------------------------------------
!
  ELSE
!
!*        3.1   Compute the intersection between the zone specified by the user
!              and the physical splitting extended by the boundaries
!
    IIB = MAX( TCRRT_COMDATA%TSPLIT_B%NXORP+IWEST, IXOR )
    IJB = MAX( TCRRT_COMDATA%TSPLIT_B%NYORP+ISOUTH, IYOR )
!
    IIE = MIN( TCRRT_COMDATA%TSPLIT_B%NXENDP+IEAST, IXEND )
    IJE = MIN( TCRRT_COMDATA%TSPLIT_B%NYENDP+INORTH, IYEND )
!
!*        3.2   Compute the dimension and allocate the buffers
!
    IDIM = (IXEND-IXOR+1)*(IYEND-IYOR+1)*(IZEND-IZOR+1)
!
    ALLOCATE(ZBUF(IXEND-IXOR+1, IYEND-IYOR+1, IZEND-IZOR+1))
    ALLOCATE(ZBUFD(IXEND-IXOR+1, IYEND-IYOR+1, IZEND-IZOR+1))
!
    ZBUF = 0
!
!*        3.3   The intersection is not empty
!
    IF((IIB <= IIE) .AND. (IJB <= IJE) ) THEN
!
!*        3.3.1 compute the displacements in the buffers
!
      IIB_BUF = IIB - IXOR + 1
      IJB_BUF = IJB - IYOR + 1
      IIE_BUF = IIE - IXOR + 1
      IJE_BUF = IJE - IYOR + 1
!
!*        3.3.2 switch the displacements in PFIELD to local indices
!
      IIB = IIB - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJB = IJB - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
      IIE = IIE - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJE = IJE - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
!
!*        3.3.3 extract the elements
!
      ZBUF(IIB_BUF:IIE_BUF, IJB_BUF:IJE_BUF, 1:IZEND-IZOR+1) = &
                                            PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND)
!
    ENDIF
!
!*        3.4    Summation with all the processors
!
    CALL MPI_ALLREDUCE(ZBUF, ZBUFD, IDIM, MPI_PRECISION, &
                       MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
!*        3.5    Return the result
!
    EXTRACT_ll => ZBUFD
!
!*        3.6    Deallocate the buffers
!
    DEALLOCATE(ZBUF)
    DEALLOCATE(ZBUFD)
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION EXTRACT_ll
!
!     ###########################################################
      FUNCTION SUM1D_ll( PFIELD, KDIR, KINFO, KXOR, KYOR, KZOR, &
                         KXEND, KYEND, KZEND )
!     ###########################################################
!
!!****  *SUM1D_ll* - function to perform a summation according a direction KDIR
!                    of a 3D field contained in a geographic region given
!                    by its first point and its last point.
!                    If these two points are omitted,
!                    the whole physical domain is considered.
! 
!                    The result is a 2D array.
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the intersection of its sub-domain
!     with the geographic region.
!     If this intersection is not empty, the processor performs
!     the summation and fills an intermediate buffer of the dimension
!     of the result. Then we make a summation of all
!     the intermediate buffers (MPI reduction) and obtain the result.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       MPI_PRECISION - 
!
!!    Implicit Arguments
!!    ------------------
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!
!!    Modifications
!!    -------------
!     Original 03/07/98
! 
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
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
!*        0.2   declarations of local variables
!
  INTEGER :: IXOR, IYOR, IZOR, & ! Real coordinates of the region
             IXEND, IYEND, IZEND
!
  REAL, ALLOCATABLE, TARGET :: ZBUF(:,:), ZBUFD(:,:) ! Intermediate buffers
!
  INTEGER :: IIE, IJE, IIB, IJB ! Displacements in PFIELD
!
  INTEGER :: IIB_BUF,IJB_BUF,IIE_BUF,IJE_BUF ! Displacements in the buffers
!
  INTEGER :: IDIM ! Dimension of the result array
!
  INTEGER :: IWEST, IEAST, & ! increments for the boundaries
             ISOUTH, INORTH  !   they are set to 0 if the subdomain
                             !   is not at a boundary
                             !   they are set to -1 or 1 if the subdomain is at
                             !   the corresponding boundary
!
  INTEGER :: IIMAX, IJMAX ! global dimensions of the physical model
!
!-------------------------------------------------------------------------------
!
  KINFO = 0
!
!*        1.    INITIALISATIONS :
!               -----------------
!*        1.1   Fill in the increments for the boundaries
!
  IWEST = 0
  IEAST = 0
  ISOUTH = 0
  INORTH = 0
!
  IF (LWEST_ll()) IWEST = -1
  IF (LEAST_ll()) IEAST = 1 
  IF (LSOUTH_ll()) ISOUTH = -1
  IF (LNORTH_ll()) INORTH = 1
!
!*        1.2   Get the global coordinates
!
  CALL GET_GLOBALDIMS_ll(IIMAX, IJMAX)
!
!*        1.3   Fill The Coordinates :
!
  IF(PRESENT(KXOR) .AND. PRESENT(KYOR) .AND. PRESENT(KZOR) .AND. &
     PRESENT(KXEND) .AND. PRESENT(KYEND) .AND. PRESENT(KZEND) ) THEN
!   
    IXOR = KXOR
    IYOR = KYOR
    IZOR = KZOR
    IXEND = KXEND
    IYEND = KYEND
    IZEND = KZEND
!
  ELSEIF(.NOT.PRESENT(KXOR) .AND. .NOT.PRESENT(KYOR) &
       .AND. .NOT.PRESENT(KZOR) .AND. .NOT.PRESENT(KXEND) &
       .AND. .NOT.PRESENT(KYEND) .AND. .NOT.PRESENT(KZEND) ) THEN
!
    IXOR = 1 + JPHEXT
    IYOR = 1 + JPHEXT
    IZOR = 1 + JPVEXT 
    IXEND = IIMAX + JPHEXT
    IYEND = IJMAX + JPHEXT
    IZEND = SIZE(PFIELD,3)-JPVEXT
!
  ELSE
!
! Error
!
    ALLOCATE(ZBUF(1,1))
    ZBUF = 0
!
    KINFO = -1
!
    SUM1D_ll => ZBUF
!
    DEALLOCATE(ZBUF)
    RETURN
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.    THE COORDINATES ARE NOT IN THE GLOBAL DOMAIN -> ERROR :
!               -----------------------------------------------------
!
  IF((IXOR < JPHEXT) .OR. (IYOR < JPHEXT) .OR. &
     (IXEND > IIMAX + 2*JPHEXT) .OR. (IYEND > IJMAX + 2*JPHEXT)) THEN
!
    ALLOCATE(ZBUF(1,1))
    ZBUF = 0
!
!   Error
!
    KINFO = -1
!
    SUM1D_ll => ZBUF
!
    DEALLOCATE(ZBUF)
    RETURN
!
!-------------------------------------------------------------------------------
!
!*        3.    THE COORDINATES ARE IN THE GLOBAL DOMAIN :
!               ----------------------------------------
!
  ELSE
!
!*        3.1   Compute the intersection
!
    IIB = MAX( TCRRT_COMDATA%TSPLIT_B%NXORP+IWEST, IXOR )
    IJB = MAX( TCRRT_COMDATA%TSPLIT_B%NYORP+ISOUTH, IYOR )
!
    IIE = MIN( TCRRT_COMDATA%TSPLIT_B%NXENDP+IEAST, IXEND )
    IJE = MIN( TCRRT_COMDATA%TSPLIT_B%NYENDP+INORTH, IYEND )
!
!*        3.2   Compute the dimension and allocate the buffers
!               according the direction
!
    IF(KDIR .EQ. 1) THEN
!
!*        3.2.1 the direction is 1
!
      ALLOCATE(ZBUF(IYEND-IYOR+1, IZEND-IZOR+1))
      ALLOCATE(ZBUFD(IYEND-IYOR+1, IZEND-IZOR+1))
      IDIM = (IYEND-IYOR+1)*(IZEND-IZOR+1)
!
    ELSEIF(KDIR .EQ. 2) THEN
!
!*        3.2.2 the direction is 2
!
      ALLOCATE(ZBUF(IXEND-IXOR+1, IZEND-IZOR+1))
      ALLOCATE(ZBUFD(IXEND-IXOR+1, IZEND-IZOR+1))
      IDIM = (IXEND-IXOR+1)*(IZEND-IZOR+1)
!
    ELSEIF(KDIR .EQ. 3) THEN
!
!*        3.2.3 the direction is 3
!
      ALLOCATE(ZBUF(IXEND-IXOR+1, IYEND-IYOR+1))
      ALLOCATE(ZBUFD(IXEND-IXOR+1, IYEND-IYOR+1))
      IDIM = (IXEND-IXOR+1)*(IYEND-IYOR+1)
!
    ELSE
!
!*        3.2.4 the direction is not 1, 2 or 3 -> Error
!
      ALLOCATE(ZBUF(1,1))
      ZBUF = 0
!
      KINFO = -1
!
      SUM1D_ll => ZBUF
!
      DEALLOCATE(ZBUF)
      RETURN
!
    ENDIF
!
    ZBUF = 0
    ZBUFD = 0
!
!*        3.3   The intersection is not empty
!
    IF((IIB <= IIE) .AND. (IJB <= IJE) ) THEN
!
!*        3.3.1 compute the displacements in the buffers
!
      IIB_BUF = IIB - IXOR + 1
      IJB_BUF = IJB - IYOR + 1
      IIE_BUF = IIE - IXOR + 1
      IJE_BUF = IJE - IYOR + 1
!
!*        3.3.2 switch the displacements in PFIELD to local indices 
!
      IIB = IIB - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJB = IJB - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
      IIE = IIE - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJE = IJE - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
!
!*        3.3.4 perform the summations
!
      IF(KDIR .EQ. 1) THEN
!
        ZBUF(IJB_BUF:IJE_BUF, 1:IZEND-IZOR+1) = &
                                    SUM(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND), 1)
!
      ELSEIF(KDIR .EQ. 2) THEN
!
        ZBUF(IIB_BUF:IIE_BUF, 1:IZEND-IZOR+1) = &
                                    SUM(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND), 2)
!
      ELSEIF(KDIR .EQ. 3) THEN
!
        ZBUF(IIB_BUF:IIE_BUF, IJB_BUF:IJE_BUF) = &
                                    SUM(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND), 3)
!
      ENDIF
!
    ENDIF
!
!*        3.4    Summation with all the processors
!
    CALL MPI_ALLREDUCE(ZBUF, ZBUFD, IDIM, MPI_PRECISION, &
                       MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
!*        3.5    Return the result
!
    SUM1D_ll => ZBUFD
!
!*        3.6    Deallocate the buffers
!
    DEALLOCATE(ZBUF)
    DEALLOCATE(ZBUFD)
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SUM1D_ll
!
!     ###################################################################
      FUNCTION SUM2D_ll( PFIELD, KDIR1, KDIR2, KINFO, KXOR, KYOR, KZOR, &
                         KXEND, KYEND, KZEND )
!     ###################################################################
!
!!****  *SUM2D_ll* - function to perform a summation according two directions
!                    KDIR1 and KDIR2  of a 3D field contained in a geographic
!                    region given by its first point and its last point.
!                    If these two points are omitted,
!                    the whole physical domain is considered.
! 
!                    The result is a 1D array.
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the intersection of its sub-domain
!     with the geographic region.
!     If this intersection is not empty, the processor performs
!     the summation and fills an intermediate buffer of the dimension
!     of the result. Then we make a summation of all
!     the intermediate buffers (MPI reduction) and obtain the result.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       MPI_PRECISION - 
!
!!    Implicit Arguments
!!    ------------------
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!
!!    Modifications
!!    -------------
!     Original 03/07/98
!     J.ESCOBAR 9/03/2008 : BUG remove deallocate(ZBUF...) or result is NULL !!!
! 
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
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
!*        0.2   declarations of local variables
!
  INTEGER :: IXOR, IYOR, IZOR, & ! Real coordinates of the region
             IXEND, IYEND, IZEND
!
  REAL, POINTER :: ZBUF(:), ZBUFD(:) ! Intermediate buffers
!
  INTEGER :: IIE, IJE, IIB, IJB ! Displacements in PFIELD
!
  INTEGER :: IIB_BUF,IJB_BUF,IIE_BUF,IJE_BUF ! Displacements in the buffers
!
  INTEGER :: IDIM ! Dimension of the result array
!
  INTEGER :: IWEST, IEAST, & ! increments for the boundaries
             ISOUTH, INORTH  !   they are set to 0 if the subdomain
                             !   is not at a boundary
                             !   they are set to -1 or 1 if the subdomain is at
                             !   the corresponding boundary
!
  INTEGER :: IIMAX, IJMAX ! global dimensions of the physical model
!
!-------------------------------------------------------------------------------
!
  KINFO = 0
!
!*        1.    INITIALISATIONS :
!               -----------------
!*        1.1   Fill in the increments for the boundaries
!
  IWEST = 0
  IEAST = 0
  ISOUTH = 0
  INORTH = 0
!
  IF (LWEST_ll()) IWEST = -1
  IF (LEAST_ll()) IEAST = 1
  IF (LSOUTH_ll()) ISOUTH = -1
  IF (LNORTH_ll()) INORTH = 1
!
!*        1.2   Get the dimensions of the current model
!
  CALL GET_GLOBALDIMS_ll(IIMAX, IJMAX)
!
!*        1.3   Fill The Coordinates :
!
  IF(PRESENT(KXOR) .AND. PRESENT(KYOR) .AND. PRESENT(KZOR) .AND. &
     PRESENT(KXEND) .AND. PRESENT(KYEND) .AND. PRESENT(KZEND) ) THEN
!   
    IXOR = KXOR
    IYOR = KYOR
    IZOR = KZOR
    IXEND = KXEND
    IYEND = KYEND
    IZEND = KZEND
!
  ELSEIF(.NOT.PRESENT(KXOR) .AND. .NOT.PRESENT(KYOR) &
        .AND. .NOT.PRESENT(KZOR) .AND. .NOT.PRESENT(KXEND) &
        .AND. .NOT.PRESENT(KYEND) .AND. .NOT.PRESENT(KZEND) ) THEN
!
    IXOR = 1 + JPHEXT
    IYOR = 1 + JPHEXT
    IZOR = 1 + JPVEXT 
    IXEND = IIMAX + JPHEXT
    IYEND = IJMAX + JPHEXT
    IZEND = SIZE(PFIELD,3)-JPVEXT
!
  ELSE
!
! Error
!
    ALLOCATE(ZBUF(1))
    ZBUF = 0
!
    KINFO = -1
!
    SUM2D_ll => ZBUF
!
!JUAN    DEALLOCATE(ZBUF)
    RETURN
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.    THE COORDINATES ARE NOT IN THE GLOBAL DOMAIN -> ERROR :
!               -----------------------------------------------------
!
  IF((IXOR < 1 ) .OR. (IYOR < 1 ) .OR. &
     (IXEND > IIMAX + 2*JPHEXT) .OR. (IYEND > IJMAX + 2*JPHEXT)) THEN
!
    ALLOCATE(ZBUF(1))
    ZBUF = 0
!
!   Error
!
    KINFO = -1
!
    SUM2D_ll => ZBUF
!JUAN    DEALLOCATE(ZBUF)
    RETURN
!
!-------------------------------------------------------------------------------
!
!*        3.    THE COORDINATES ARE IN THE GLOBAL DOMAIN :
!               ----------------------------------------
!
  ELSE
!
!*        3.1   Compute the intersection
!
    IIB = MAX( TCRRT_COMDATA%TSPLIT_B%NXORP+IWEST, IXOR )
    IJB = MAX( TCRRT_COMDATA%TSPLIT_B%NYORP+ISOUTH, IYOR )
!
    IIE = MIN( TCRRT_COMDATA%TSPLIT_B%NXENDP+IEAST, IXEND )
    IJE = MIN( TCRRT_COMDATA%TSPLIT_B%NYENDP+INORTH, IYEND )
!
!*        3.2   Compute the dimension and allocate the buffers
!               according the directions
!
    IF ((KDIR1 .EQ. 1) .AND. (KDIR2 .EQ. 2) ) THEN
!
!*        3.2.1 the directions are 1 and 2
!
      ALLOCATE(ZBUF(IZEND-IZOR+1))
      ALLOCATE(ZBUFD(IZEND-IZOR+1))
      IDIM = (IZEND-IZOR+1)
!
    ELSEIF( (KDIR1 .EQ. 2) .AND. (KDIR2 .EQ. 3) ) THEN
!
!*        3.2.2 the directions are 2 and 3
!
      ALLOCATE(ZBUF(IXEND-IXOR+1))
      ALLOCATE(ZBUFD(IXEND-IXOR+1))
      IDIM = (IXEND-IXOR+1)
!
!*        3.2.3 the directions are 1 and 3
!
    ELSEIF( (KDIR1 .EQ. 1) .AND. (KDIR2 .EQ. 3) ) THEN
!
      ALLOCATE(ZBUF(IYEND-IYOR+1))
      ALLOCATE(ZBUFD(IYEND-IYOR+1))
      IDIM = (IYEND-IYOR+1)
!
    ELSE
!
!*        3.2.4 Error of directions
!
      ALLOCATE(ZBUF(1))
      ZBUF = 0
!
      KINFO = -1
!
      SUM2D_ll => ZBUF
!JUAN      DEALLOCATE(ZBUF)
      RETURN
!
    ENDIF
!
    ZBUF = 0
    ZBUFD = 0
!
!*        3.3   The intersection is not empty
!
    IF((IIB <= IIE) .AND. (IJB <= IJE) ) THEN
!
!*        3.3.1 compute the displacements in the buffers
!
      IIB_BUF = IIB - IXOR + 1
      IJB_BUF = IJB - IYOR + 1
      IIE_BUF = IIE - IXOR + 1
      IJE_BUF = IJE - IYOR + 1
!
!*        3.3.2 switch the displacements in PFIELD to local indices
!
      IIB = IIB - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJB = IJB - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
      IIE = IIE - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJE = IJE - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
!
!*        3.3.4 perform the summations
!
      IF ((KDIR1 .EQ. 1) .AND. (KDIR2 .EQ. 2) ) THEN
!
        ZBUF(1:IZEND-IZOR+1) = &
                            SUM(SUM(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND), 2), 1)
!
      ELSEIF( (KDIR1 .EQ. 2) .AND. (KDIR2 .EQ. 3) ) THEN
!
        ZBUF(IIB_BUF:IIE_BUF) = &
                            SUM(SUM(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND), 3), 2)
!
      ELSEIF( (KDIR1 .EQ. 1) .AND. (KDIR2 .EQ. 3) ) THEN
!
        ZBUF(IJB_BUF:IJE_BUF) = &
                            SUM(SUM(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND), 3), 1)
!
      ENDIF
!
    ENDIF
!
!*        3.4    Summation with all the processors
!
    CALL MPI_ALLREDUCE(ZBUF, ZBUFD, IDIM, MPI_PRECISION, &
                       MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
!*        3.5    Return the result
!
    SUM2D_ll => ZBUFD
!
!*        3.6    Deallocate the buffers
!
    DEALLOCATE(ZBUF)
!    DEALLOCATE(ZBUFD)
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SUM2D_ll
!
!     ##########################################################
      REAL FUNCTION SUM3D_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                              KXEND, KYEND, KZEND )
!     ##########################################################
!
!!****  *SUM3D_ll* - function to perform a summation in the three directions
!                    of a 3D field contained in a geographic region given
!                    by its first point and its last point.
!                    If these two points are omitted,
!                    the whole physical domain is considered.
! 
!                    The result is a scalar
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the intersection of its sub-domain
!     with the geographic region.
!     If this intersection is not empty, the processor performs
!     the summation and compute a local sum.
!     Then we make a summation of all these local sums (MPI reduction)
!     and obtain the result.
!
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       MPI_PRECISION - 
!
!!    Implicit Arguments
!!    ------------------
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!
!!    Modifications
!!    -------------
!     Original 03/07/98
!     Modif : V.Masson , 08/2010 : for reproducibility
! 
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  USE MODE_TOOLS_ll
!
!JUAN
!!$  USE MODE_GATHER_ll
 USE MODE_REPRO_SUM 
!JUAN
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
!*        0.2   declarations of local variables
!
  REAL :: ZSUM3D ! Intermediate Sum
!
  INTEGER :: IXOR, IYOR, IZOR, & ! Real coordinates of the region
             IXEND, IYEND, IZEND
!
  INTEGER :: IIE, IJE, IIB, IJB ! Displacements in PFIELD
!
  INTEGER :: IWEST, IEAST, & ! increments for the boundaries
             ISOUTH, INORTH  !   they are set to 0 if the subdomain
                             !   is not at a boundary
                             !   they are set to -1 or 1 if the subdomain is at
                             !   the corresponding boundary
!
  INTEGER :: IIMAX, IJMAX ! global dimensions of the physical model
!
REAL, DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2)) :: ZSUM
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSUM_ll 
!-------------------------------------------------------------------------------
!
  KINFO = 0
!
  ZSUM3D = 0
!
!*        1.    INITIALISATIONS :
!               -----------------
!*        1.1   Fill in the increments for the boundaries
!
  IWEST = 0
  IEAST = 0
  ISOUTH = 0
  INORTH = 0
!
  IF (LWEST_ll()) IWEST = -1
  IF (LEAST_ll()) IEAST = 1
  IF (LSOUTH_ll()) ISOUTH = -1
  IF (LNORTH_ll()) INORTH = 1
!
!*        1.2   Get the dimensions of the current model
!
  CALL GET_GLOBALDIMS_ll(IIMAX,IJMAX)
!
!*        1.3   Fill The Coordinates :
!
  IF(PRESENT(KXOR) .AND. PRESENT(KYOR) .AND. PRESENT(KZOR) .AND. &
     PRESENT(KXEND) .AND. PRESENT(KYEND) .AND. PRESENT(KZEND) ) THEN
!   
    IXOR = KXOR
    IYOR = KYOR
    IZOR = KZOR
    IXEND = KXEND
    IYEND = KYEND
    IZEND = KZEND
!
  ELSEIF(.NOT.PRESENT(KXOR) .AND. .NOT.PRESENT(KYOR) &
        .AND. .NOT.PRESENT(KZOR) .AND. .NOT.PRESENT(KXEND) &
        .AND. .NOT.PRESENT(KYEND) .AND. .NOT.PRESENT(KZEND) ) THEN
!
    IXOR = 1 + JPHEXT
    IYOR = 1 + JPHEXT
    IZOR = 1 + JPVEXT 
    IXEND = IIMAX + JPHEXT
    IYEND = IJMAX + JPHEXT
    IZEND = SIZE(PFIELD,3)-JPVEXT
!
  ELSE
!
! Error
!
    KINFO = -1
!
    SUM3D_ll = 0
!
    RETURN
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.    THE COORDINATES ARE NOT IN THE GLOBAL DOMAIN -> ERROR :
!               -----------------------------------------------------
!
  IF((IXOR < 1 ) .OR. (IYOR < 1 ) .OR. &
     (IXEND > IIMAX + 2*JPHEXT) .OR. (IYEND > IJMAX + 2*JPHEXT)) THEN
!
    SUM3D_ll = 0
!
!   Error
!
    KINFO = -1
    RETURN
!
!-------------------------------------------------------------------------------
!
!*        3.    THE COORDINATES ARE IN THE GLOBAL DOMAIN :
!               ----------------------------------------
!
  ELSE
!
!*        3.1   Compute the intersection
!
    IIB = MAX( TCRRT_COMDATA%TSPLIT_B%NXORP+IWEST, IXOR )
    IJB = MAX( TCRRT_COMDATA%TSPLIT_B%NYORP+ISOUTH, IYOR )
!
    IIE = MIN( TCRRT_COMDATA%TSPLIT_B%NXENDP+IEAST, IXEND )
    IJE = MIN( TCRRT_COMDATA%TSPLIT_B%NYENDP+INORTH, IYEND )
!
!*        3.2   The intersection is not empty
!
    ZSUM(:,:) = 0.
    IF((IIB <= IIE) .AND. (IJB <= IJE) ) THEN
!
!*        3.2.1 compute the displacements in PFIELD
!
      IIB = IIB - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJB = IJB - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
      IIE = IIE - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJE = IJE - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
!
!*        3.2.2 perform the summations
!
!!$      ZSUM3D = SUM(SUM(SUM(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND), 3), 2), 1)
       ! vertical sum only, on each processor
       ZSUM(:,:) = 0.
       ZSUM(IIB:IIE, IJB:IJE) = SUM(PFIELD(IIB:IIE, IJB:IJE,IZOR:IZEND),3)
!
    ENDIF
!
!*        3.4    Summation with all the processors
!
!!$    CALL MPI_ALLREDUCE(ZSUM3D, SUM3D_ll, 1, MPI_PRECISION, &
!!$                       MPI_SUM, NMNH_COMM_WORLD, KINFO)
!!$!
!!$! gathers the total 2D field
!!$
!!$    ALLOCATE(ZSUM_ll(IIMAX+2,IJMAX+2))
!!$    CALL GATHERALL_FIELD_ll('XY',ZSUM,ZSUM_ll,KINFO)
!!$
!!$!* computes the sum
!!$    SUM3D_ll= SUM(ZSUM_ll(:,:))
!!$    DEALLOCATE(ZSUM_ll)
!JUAN    
    SUM3D_ll= SUM_DD_R2_ll(ZSUM(:,:))
!JUAN    
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SUM3D_ll
!
!     ###################################################################
      FUNCTION SUM_1DFIELD_ll(PFIELD, HDIR, KOR, KEND, KERR) RESULT(ZSUM)
!     ###################################################################
!
!!****  *SUM_1DFIELD_ll* - this function calculates the sum of the PFIELD 
!                          one-dimensional field according to the HDIR direction
!                          The KOR and KEND arguments can be used to specify
!                          the boundaries of the sum.
!                          The dimension of PFIELD is supposed to be
!                          the dimension of the extended subdomain
!                          in the HDIR direction.
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!       The sizes and positions of each local array are gathered on each proc.
!     Then the MPI_ALLGATHERV routine is called to gather the global 1D field.
!     Then the sum of the 1D global field is calculated.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NPROC -
!       MPI_PRECISION - 
!
!!    Reference
!!    ---------
!     User Guide for the MesoNH Parallel Package
!     L. Giraud, R. Guivarch, P. Kloos, D. Lugato (CERFACS)
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    MODIFICATIONS
!!    -------------
!     Original 16/09/98
!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll
  USE MODD_VAR_ll, ONLY : IP, TCRRT_COMDATA, TCRRT_PROCONF, NPROC, MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY :  LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
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
!*        0.2   declarations of local variables
!
  REAL, DIMENSION(:), ALLOCATABLE :: ZGLOBFIELD ! global field
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
  INTEGER, DIMENSION(2) :: IMAX ! global dimensions of the current model 
  INTEGER, DIMENSION(2) :: IOR, IEND, IORP, IENDP, IORE, IENDE
  INTEGER :: IDIR ! direction of the 1D field
  INTEGER :: IERR
  INTEGER :: ISIZE
  INTEGER, DIMENSION(NPROC) :: IDISPL, ISIZES
  INTEGER :: IB, IE
  INTEGER :: IGB, IGE
  INTEGER :: IWEST, IEAST, INORTH, ISOUTH
!
!-------------------------------------------------------------------------------
!
!*        1.    INITIALISATIONS
!               ---------------
!
!*        1.1   Get current splitting
!
  IWEST = 0
  IEAST = 0
  INORTH = 0
  ISOUTH = 0
  IF (LWEST_ll())  IWEST=-1
  IF (LEAST_ll())  IEAST=1
  IF (LNORTH_ll()) INORTH=1
  IF (LSOUTH_ll()) ISOUTH=-1
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
  IOR(1)  = TZSPLIT%NXORP+IWEST
  IOR(2)  = TZSPLIT%NYORP+ISOUTH
  IEND(1) = TZSPLIT%NXENDP+IEAST
  IEND(2) = TZSPLIT%NYENDP+INORTH
!
  IORP(1)  = TZSPLIT%NXORP
  IORP(2)  = TZSPLIT%NYORP
  IENDP(1) = TZSPLIT%NXENDP
  IENDP(2) = TZSPLIT%NYENDP
!
  IORE(1)  = TZSPLIT%NXORE
  IORE(2)  = TZSPLIT%NYORE
  IENDE(1) = TZSPLIT%NXENDE
  IENDE(2) = TZSPLIT%NYENDE
!
!*        1.2   Get the dimensions of the current model
!
  CALL GET_GLOBALDIMS_ll(IMAX(1), IMAX(2))
!
!*        1.3   Set the direction (1 for X and 2 for Y)
!
  IDIR = ICHAR(HDIR) - ICHAR('X') + 1
!
!
!*        1.4   Set the boundaries
!
  IF (PRESENT(KOR) .AND. PRESENT(KEND)) THEN
!
!*      Test the ranges
!
    IF((KOR < 1 ) .OR. (KEND > IMAX(IDIR)+2*JPHEXT)) THEN
!
!   Error
!
      KERR = -1
      RETURN
!
    ENDIF
!
    IB = MAX(KOR, IOR(IDIR)) 
    IE = MIN(KEND, IEND(IDIR)) 
!
  ELSE ! default : physical zone
!
    IB = IORP(IDIR)
    IE = IENDP(IDIR)
!
  ENDIF
!
!*        1.5  Get global begin and end
!
  CALL MPI_ALLREDUCE(IB, IGB, 1, MPI_INTEGER, MPI_MIN, NMNH_COMM_WORLD, IERR)
  CALL MPI_ALLREDUCE(IE, IGE, 1, MPI_INTEGER, MPI_MAX, NMNH_COMM_WORLD, IERR)
!
!-------------------------------------------------------------------------------
!
!*        2.   GET THE GLOBAL FIELD
!              --------------------
!
!*        2.1  Have the sizes and global positions known by all procs
!
  ISIZE = IE - IB + 1
  CALL MPI_ALLGATHER( (/ ISIZE /) , 1, MPI_INTEGER, ISIZES, 1, MPI_INTEGER, &
                     NMNH_COMM_WORLD, IERR)
!
  CALL MPI_ALLGATHER( (/ IB-IGB /), 1, MPI_INTEGER, IDISPL, 1, MPI_INTEGER, &
                     NMNH_COMM_WORLD, IERR)
! 
!*        2.2  Get the global field
!
  ALLOCATE(ZGLOBFIELD(IGE-IGB+1))
  CALL MPI_ALLGATHERV(PFIELD(IB-IORE(IDIR)+1:), ISIZE, &
                      MPI_PRECISION, ZGLOBFIELD, ISIZES, IDISPL, &
                      MPI_PRECISION, NMNH_COMM_WORLD, IERR)
!
!-------------------------------------------------------------------------------
!
!*        3.   CALCULATE THE SUM
!              -----------------
!
!OCL SCALAR
  ZSUM = SUM(ZGLOBFIELD(1:IGE-IGB+1))
!
  DEALLOCATE(ZGLOBFIELD)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SUM_1DFIELD_ll
!
!     ########################################################
      REAL FUNCTION MAX_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                            KXEND, KYEND, KZEND )
!     ########################################################
!
!!****  *MAX_ll* - function to find the maximum value in a 3D field
!                  contained in a geographic region given by its first point
!                  and its last point. If these two points are omitted,
!                  the whole physical domain is considered.
! 
!                  The result is a scalar
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the intersection of its sub-domain
!     with the geographic region.
!     If this intersection is not empty, the processor searchs
!     the maximum in this intersection.
!     Then we make a MPI reduction with all the local maxima.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NPROC -
!       MPI_PRECISION - 
!
!!    Reference
!!    ---------
!     User Guide for the MesoNH Parallel Package
!     L. Giraud, R. Guivarch, P. Kloos, D. Lugato (CERFACS)
!
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!
!!    MODIFICATIONS
!!    -------------
!     Original 03/07/98
!!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
!*       0.2   declarations of local variables
!
  REAL :: ZMAX ! Intermediate Max
!
  INTEGER :: IXOR, IYOR, IZOR, & ! Real coordinates of the region
             IXEND, IYEND, IZEND
!
  INTEGER :: IIE, IJE, IIB, IJB ! Displacements in PFIELD
!
  INTEGER :: IWEST, IEAST, & ! increments for the boundaries
             ISOUTH, INORTH  !   they are set to 0 if the subdomain
                             !   is not at a boundary
                             !   they are set to -1 or 1 if the subdomain is at
                             !   the corresponding boundary
!
  INTEGER :: IIMAX, IJMAX ! global dimensions of the physical model
!
!-------------------------------------------------------------------------------
!
  KINFO = 0
!
  ZMAX = 0
!
!*        1.    INITIALISATIONS :
!               -----------------
!*        1.1   Fill in the increments for the boundaries
!
  IWEST = 0
  IEAST = 0
  ISOUTH = 0
  INORTH = 0
!
  IF (LWEST_ll()) IWEST = -1
  IF (LEAST_ll()) IEAST = 1
  IF (LSOUTH_ll()) ISOUTH = -1
  IF (LNORTH_ll()) INORTH = 1
!
!*        1.2   Get the dimensions of the current model
!
  CALL GET_GLOBALDIMS_ll(IIMAX, IJMAX)
!
!*        1.3   Fill The Coordinates :
!
  IF(PRESENT(KXOR) .AND. PRESENT(KYOR) .AND. PRESENT(KZOR) .AND. &
     PRESENT(KXEND) .AND. PRESENT(KYEND) .AND. PRESENT(KZEND) ) THEN
!   
    IXOR = KXOR
    IYOR = KYOR
    IZOR = KZOR
    IXEND = KXEND
    IYEND = KYEND
    IZEND = KZEND
!
  ELSEIF(.NOT.PRESENT(KXOR) .AND. .NOT.PRESENT(KYOR) &
       .AND. .NOT.PRESENT(KZOR) .AND. .NOT.PRESENT(KXEND) &
       .AND. .NOT.PRESENT(KYEND) .AND. .NOT.PRESENT(KZEND) ) THEN
!
    IXOR = 1 + JPHEXT
    IYOR = 1 + JPHEXT
    IZOR = 1 + JPVEXT 
    IXEND = IIMAX + JPHEXT
    IYEND = IJMAX + JPHEXT
    IZEND = SIZE(PFIELD,3)-JPVEXT
!
  ELSE
!
!   Error
!
    KINFO = -1
!
    MAX_ll = 0
!
    RETURN
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.    THE COORDINATES ARE NOT IN THE GLOBAL DOMAIN -> ERROR :
!               -----------------------------------------------------
!
  IF((IXOR < 1 ) .OR. (IYOR < 1 ) .OR. &
     (IXEND > IIMAX + 2*JPHEXT) .OR. (IYEND > IJMAX + 2*JPHEXT)) THEN
!
    MAX_ll = 0
!
!   Error
!
    KINFO = -1
    RETURN
!
!-------------------------------------------------------------------------------
!
!*        3.    THE COORDINATES ARE IN THE GLOBAL DOMAIN :
!               ----------------------------------------
!
  ELSE
!
!*        3.1   Compute the intersection
!
    IIB = MAX( TCRRT_COMDATA%TSPLIT_B%NXORP+IWEST, IXOR )
    IJB = MAX( TCRRT_COMDATA%TSPLIT_B%NYORP+ISOUTH, IYOR )
!
    IIE = MIN( TCRRT_COMDATA%TSPLIT_B%NXENDP+IEAST, IXEND )
    IJE = MIN( TCRRT_COMDATA%TSPLIT_B%NYENDP+INORTH, IYEND )
!
!*        3.2   The intersection is not empty
!
    IF((IIB <= IIE) .AND. (IJB <= IJE) ) THEN
!
!*        3.2.1 compute the displacements in PFIELD
!
      IIB = IIB - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJB = IJB - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
      IIE = IIE - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJE = IJE - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
!
!*        3.2.2 perform the search
!

      ZMAX = MAXVAL(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND))
!
    ENDIF
!
!*        3.4    Reduction with all the processors
!
    CALL MPI_ALLREDUCE(ZMAX, MAX_ll, 1, MPI_PRECISION, &
                       MPI_MAX, NMNH_COMM_WORLD, KINFO)
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION MAX_ll
!
!     ########################################################
      REAL FUNCTION MIN_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                            KXEND, KYEND, KZEND )
!     ########################################################
!
!!****  *MIN_ll* - function to find the minimum value in a 3D field
!                  contained in a geographic region given by its first point
!                  and its last point. If these two points are omitted,
!                  the whole physical domain is considered.
! 
!                  the result is a scalar
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the intersection of its sub-domain
!     with the geographic region.
!     If this intersection is not empty, the processor searchs the minimum
!     in this intersection.
!     Then we make a MPI reduction with all the local minima.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPVEXT
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NPROC -
!       MPI_PRECISION - 
!
!!    Reference
!!    ---------
!     User Guide for the MesoNH Parallel Package
!     L. Giraud, R. Guivarch, P. Kloos, D. Lugato (CERFACS)
!
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!
!!    MODIFICATIONS
!!    -------------
!     Original 03/07/98
!!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
  INTEGER, OPTIONAL, INTENT(IN) :: KXOR, KYOR, KZOR, & ! Coordinates
                                   KXEND, KYEND, KZEND ! of the region
!
!*        0.2   declarations of local variables
!
  REAL :: ZMIN ! Intermediate Min
!
  INTEGER :: IXOR, IYOR, IZOR, & ! Real coordinates of the region
             IXEND, IYEND, IZEND
!
  INTEGER :: IIE, IJE, IIB, IJB ! Displacements in PFIELD
!
  INTEGER :: IWEST, IEAST, & ! increments for the boundaries
             ISOUTH, INORTH  !   they are set to 0 if the subdomain
                             !   is not at a boundary
                             !   they are set to -1 or 1 if the subdomain is at
                             !   the corresponding boundary
!
  INTEGER :: IIMAX, IJMAX ! global dimensions of the physical model
!
!-------------------------------------------------------------------------------
!
  KINFO = 0
!
  ZMIN = 0
!
!*        1.    INITIALISATIONS :
!              -----------------
!*        1.1   Fill in the increments for the boundaries
!
  IWEST = 0
  IEAST = 0
  ISOUTH = 0
  INORTH = 0
!
  IF (LWEST_ll()) IWEST = -1
  IF (LEAST_ll()) IEAST = 1
  IF (LSOUTH_ll()) ISOUTH = -1
  IF (LNORTH_ll()) INORTH = 1
!
!*        1.2   Get the dimensions of the current model
!
  CALL GET_GLOBALDIMS_ll(IIMAX, IJMAX)
!
!*        1.3   Fill The Coordinates :
!
  IF(PRESENT(KXOR) .AND. PRESENT(KYOR) .AND. PRESENT(KZOR) .AND. &
     PRESENT(KXEND) .AND. PRESENT(KYEND) .AND. PRESENT(KZEND) ) THEN
!   
    IXOR = KXOR
    IYOR = KYOR
    IZOR = KZOR
    IXEND = KXEND
    IYEND = KYEND
    IZEND = KZEND
!
  ELSEIF(.NOT.PRESENT(KXOR) .AND. .NOT.PRESENT(KYOR) &
       .AND. .NOT.PRESENT(KZOR) .AND. .NOT.PRESENT(KXEND) &
       .AND. .NOT.PRESENT(KYEND) .AND. .NOT.PRESENT(KZEND) ) THEN
!
    IXOR = 1 + JPHEXT
    IYOR = 1 + JPHEXT
    IZOR = 1 + JPVEXT 
    IXEND = IIMAX + JPHEXT
    IYEND = IJMAX + JPHEXT
    IZEND = SIZE(PFIELD,3)-JPVEXT
!
  ELSE
!
! Error
!
    KINFO = -1
!
    MIN_ll = 0
!
    RETURN
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.    THE COORDINATES ARE NOT IN THE GLOBAL DOMAIN -> ERROR :
!               -----------------------------------------------------
!
  IF((IXOR < 1 ) .OR. (IYOR < 1 ) .OR. &
     (IXEND > IIMAX + 2*JPHEXT) .OR. (IYEND > IJMAX + 2*JPHEXT)) THEN
!
    MIN_ll = 0
!
!   Error
!
    KINFO = -1
    RETURN
!
!-------------------------------------------------------------------------------
!
!*        3.    THE COORDINATES ARE IN THE GLOBAL DOMAIN :
!               ----------------------------------------
!
  ELSE
!
!         3.1   Compute the intersection
!
    IIB = MAX( TCRRT_COMDATA%TSPLIT_B%NXORP+IWEST, IXOR )
    IJB = MAX( TCRRT_COMDATA%TSPLIT_B%NYORP+ISOUTH, IYOR )
!
    IIE = MIN( TCRRT_COMDATA%TSPLIT_B%NXENDP+IEAST, IXEND )
    IJE = MIN( TCRRT_COMDATA%TSPLIT_B%NYENDP+INORTH, IYEND )
!
!         3.2   The intersection is not empty
!
    IF((IIB <= IIE) .AND. (IJB <= IJE) ) THEN
!
!         3.2.1 compute the displacements in PFIELD
!
      IIB = IIB - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJB = IJB - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
      IIE = IIE - TCRRT_COMDATA%TSPLIT_B%NXORE + 1
      IJE = IJE - TCRRT_COMDATA%TSPLIT_B%NYORE + 1
!
!         3.2.2 perform the search
!
      ZMIN = MINVAL(PFIELD(IIB:IIE, IJB:IJE, IZOR:IZEND))
!
    ENDIF
!
!        3.4    Reduction with all the processors
!
    CALL MPI_ALLREDUCE(ZMIN, MIN_ll, 1, MPI_PRECISION, &
                       MPI_MIN, NMNH_COMM_WORLD, KINFO)
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION MIN_ll
!
!     #########################################
      FUNCTION SUMMASK_ll(PFIELD, OMASK, KINFO)
!     #########################################
!
!!****  *SUMMASK_ll* function to calculate the sum of the 3D field PFIELD
!                    on the whole domain according the directions x and y
!                    for the points that pass the 2D-horizontal mask OMASK.
! 
!                    The result is a 1D array of dimension SIZE(PFIELD,3).
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the summation in its subdomain
!     then we make a MPI reduction with the local result.
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       MPI_PRECISION - 
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 03/07/98
!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION
!
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
  REAL, DIMENSION(:), POINTER :: SUMMASK_ll ! Result
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASK ! 2d Mask
!
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
!*        0.2   declarations of local variables
!
  REAL, DIMENSION(SIZE(PFIELD,3)), TARGET :: ZBUF, ZBUFD ! Intermediate buffers
!
  INTEGER :: JK ! Loop control 
!
!-------------------------------------------------------------------------------
!
!*        1.    PERFOMS THE LOCAL SUM :
!               ---------------------
!
  DO JK = 1, SIZE(PFIELD,3)
    ZBUF(JK) = SUM(SUM(PFIELD(:,:,JK),1,OMASK(:,:)),1)
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       2.    REDUCTION WITH ALL THE PROCESSORS :
!              ---------------------------------
!
  CALL MPI_ALLREDUCE(ZBUF, ZBUFD, SIZE(PFIELD,3), &
                     MPI_PRECISION, MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
!-------------------------------------------------------------------------------
!
!        3.    RETURN RESULT :
!              -------------
!
  SUMMASK_ll => ZBUFD
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SUMMASK_ll
!
!     ##################################################
      REAL FUNCTION SUMMASKCOMP_ll(PFIELD, OMASK, KINFO)
!     ##################################################
!
!!****  *SUMMASKCOMP_ll* function to calculate the sum of the 3D field PFIELD
!                        on the whole domain according the three directions
!                        for the points that pass the 2D-horizontal mask OMASK.
! 
!                        The result is a scalar
! 
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     Each processor computes the summation in its subdomain
!     then we make a MPI reduction with the local result.
! 
!     External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       MPI_PRECISION - 
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 03/07/98
!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : MPI_PRECISION
!
  IMPLICIT NONE
!
!*        0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD ! 3d Field
!
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASK ! 2d Mask
! 
  INTEGER, INTENT(OUT) :: KINFO ! Returned Info
!
!*        0.2   declarations of local variables
!
  REAL :: ZSUM ! Intermediate Sum
!
  REAL, DIMENSION(SIZE(PFIELD,3)) :: ZBUF ! Intermediate buffer
!
  INTEGER :: JK ! Loop control
!
!-------------------------------------------------------------------------------
!
!*        1.    PERFOMS THE LOCAL SUM :
!               ---------------------
!
  DO JK = 1, SIZE(PFIELD,3)
    ZBUF(JK) = SUM(SUM(PFIELD(:,:,JK),1,OMASK(:,:)),1)
  ENDDO
!
  ZSUM = SUM(ZBUF,1)
!
!-------------------------------------------------------------------------------
!
!*        2.    REDUCTION WITH ALL THE PROCESSORS :
!               ---------------------------------
!
  CALL MPI_ALLREDUCE(ZSUM, SUMMASKCOMP_ll, 1, MPI_PRECISION, &
                     MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION SUMMASKCOMP_ll
!
!     ###########################################
      SUBROUTINE SUM_DIM2_ll(PFIELD, PRES, KINFO)
!     ###########################################
!
!!****  *SUM_DIM2_ll*-
!
!!    Purpose
!!    -------
!     The PFIELD argument is a local 1D array according the x-direction,
!     result of local summations in y-direction.
!     The purpose of this routine is to merge all the local sum arrays in a
!     global one PRES.
! 
!!    Method
!!    ------
!       Each processor fills its part of PRES array in an intermediate buffer
!     ZBUF with its local PFIELD, then we reduce the buffer in the result PRES.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       IP -
!       MPI_PRECISION -
!       JPHALO -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll
!
  USE MODD_VAR_ll, ONLY : IP, TCRRT_COMDATA, TCRRT_PROCONF, JPHALO, &
                          MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LWEST_ll, LEAST_ll
!
  IMPLICIT NONE
!
!*        0.1   Declarations of dummy arguments :
!
  REAL, DIMENSION(:), INTENT(IN) :: PFIELD
  REAL, DIMENSION(:), INTENT(OUT) :: PRES
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*        0.2   Declarations of local variables :
!
  INTEGER :: IIB, IIE, IIBG, IIEG ! local and global displacements
!
  REAL, DIMENSION(SIZE(PRES)) :: ZBUF
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT ! Intermediate model splitting
!
!-------------------------------------------------------------------------------
!
!*        1.    GET THE CURRENT SPLITTING CONFIGURATION
!               AND COMPUTE DISPLACEMENTS
!               ---------------------------------------
!
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
!
  IF (LWEST_ll()) THEN
    IIB  = 1
    IIBG = TZSPLIT%NXORE 
  ELSE
    IIB  = 1+JPHALO
    IIBG = TZSPLIT%NXORP
  ENDIF
!
  IF (LEAST_ll()) THEN
    IIE  = SIZE(PFIELD)
    IIEG = TZSPLIT%NXENDE
  ELSE
    IIE  = SIZE(PFIELD)-JPHALO
    IIEG = TZSPLIT%NXENDP
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.    FILL THE INTERMEDIATE BUFFER
!               ----------------------------
!
  ZBUF = 0
  ZBUF(IIBG:IIEG) = PFIELD(IIB:IIE)
!
!-----------------------------------------------------------------
!
!*        3.    MERGE LOCAL SUMS
!               ----------------
!
  CALL MPI_ALLREDUCE(ZBUF, PRES, SIZE(PRES), MPI_PRECISION, MPI_SUM, &
                     NMNH_COMM_WORLD,KINFO)
!
!-----------------------------------------------------------------
!
      END SUBROUTINE SUM_DIM2_ll
!
!     #####################################################
      SUBROUTINE SUM_DIM1_DD_ll(PFIELD, PRES, KDIM, KINFO)
!     #####################################################
!
!!****  *SUM_DIM1_DD_ll*-
!
!!    Purpose
!!    -------
!     The PFIELD argument is a local 1D array according the y-direction,
!     result of local summations in x-direction.
!     The purpose of this routine is to merge all the local sum arrays in a
!     global one PRES.
! 
!!    Method
!!    ------
!     Each processor fills its part of PRES array in an intermediate buffer
!     ZBUF with its local PFIELD, then we reduce the buffer in the result PRES.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LEAST_ll, LWEST_ll, LNORTH_ll, LSOUTH_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       IP -
!       MPI_PRECISION -
!       JPHALO -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!              05/99 : P. Jabouille - N. Gicquel
!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll
!
  USE MODD_VAR_ll, ONLY : IP, TCRRT_COMDATA, TCRRT_PROCONF, JPHALO, &
                          MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LWEST_ll, LEAST_ll, LNORTH_ll, LSOUTH_ll
!
 USE MODE_REPRO_SUM
!
!*       0.    DECLARATIONS
!              ------------
!
  IMPLICIT NONE
!
!
!*       0.    DECLARATIONS
!              ------------
!
!
!*       0.1   Declarations of dummy arguments :
!
  REAL, DIMENSION(:,:), INTENT(IN)   :: PFIELD
  REAL, DIMENSION(:), INTENT(OUT)    :: PRES
  INTEGER              , INTENT(IN)  :: KDIM
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  INTEGER :: IXB, IXE, IXBG, IXEG, IJB, IJE, IJBG, IJEG ! local and global displacements
!
   TYPE(DOUBLE_DOUBLE), DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2)) :: ZBUF
   TYPE(DOUBLE_DOUBLE), DIMENSION(SIZE(PRES))                    :: ZBUF1D_ll
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT ! Intermediate model splitting
!
  INTEGER :: JI,JJ
!
!-----------------------------------------------------------------
!
!*       1.    Get the current splitting configuration and compute displacements
!
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
!
  ZBUF%R      = 0.0   ; ZBUF%E      = 0.0
  ZBUF1D_ll%R = 0.0   ; ZBUF1D_ll%E = 0.0
!
  IF (KDIM.EQ.1) THEN
     IF (LSOUTH_ll()) THEN
        IJB  = 1
        IJBG = TZSPLIT%NYORE
     ELSE
        IJB  = 1+JPHALO
        IJBG = TZSPLIT%NYORP
     ENDIF
!
     IF (LNORTH_ll()) THEN
        IJE  = SIZE(PFIELD, 2)
        IJEG = TZSPLIT%NYENDE
     ELSE
        IJE  = SIZE(PFIELD, 2)-JPHALO
        IJEG = TZSPLIT%NYENDP
     ENDIF
!
!-----------------------------------------------------------------
!
!*       2.    Fill the intermediate buffer
!
     ZBUF(:,IJB:IJE)%R = PFIELD(:,IJB:IJE)
     DO JJ=IJBG,IJEG
        ZBUF1D_ll(JJ) = SUM_DD_DD1 (ZBUF(:,JJ-IJBG+IJB))
     END DO
!
!-----------------------------------------------------------------
!
!*       3.    Merge local sums
!
    CALL  REDUCE_SUM_1DD_ll(ZBUF1D_ll, KINFO)
    PRES = ZBUF1D_ll%R
!
!-----------------------------------------------------------------
  ELSE
     IF (KDIM.EQ.2) THEN
        IF (LWEST_ll()) THEN
           IXB  = 1
           IXBG = TZSPLIT%NXORE
        ELSE
           IXB  = 1+JPHALO
           IXBG = TZSPLIT%NXORP
        ENDIF
!
        IF (LEAST_ll()) THEN
           IXE  = SIZE(PFIELD, 1)
           IXEG = TZSPLIT%NXENDE
        ELSE
           IXE  = SIZE(PFIELD, 1)-JPHALO
           IXEG = TZSPLIT%NXENDP
        ENDIF
!
!-----------------------------------------------------------------
!
!*       2.    Fill the intermediate buffer
!
        ZBUF(IXB:IXE,:)%R = PFIELD(IXB:IXE,:)
        DO JI=IXBG,IXEG
           ZBUF1D_ll(JI) = SUM_DD_DD1 (ZBUF(JI-IXBG+IXB,:))
        END DO
!
!-----------------------------------------------------------------
!
!*       2.    Merge local sums
!
        CALL  REDUCE_SUM_1DD_ll(ZBUF1D_ll, KINFO)     
        PRES = ZBUF1D_ll%R
     ENDIF
  ENDIF
  
      END SUBROUTINE SUM_DIM1_DD_ll
!     #####################################################
      SUBROUTINE SUM_DIM1_ll(PFIELD, PRES, KINFO)
!     #####################################################
!
!!****  *SUM_DIM1_ll*-
!
!!    Purpose
!!    -------
!     The PFIELD argument is a local 1D array according the y-direction,
!     result of local summations in x-direction.
!     The purpose of this routine is to merge all the local sum arrays in a
!     global one PRES.
! 
!!    Method
!!    ------
!     Each processor fills its part of PRES array in an intermediate buffer
!     ZBUF with its local PFIELD, then we reduce the buffer in the result PRES.
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       LEAST_ll, LWEST_ll, LNORTH_ll, LSOUTH_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       IP -
!       MPI_PRECISION -
!       JPHALO -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!              05/99 : P. Jabouille - N. Gicquel
!
!-------------------------------------------------------------------------------
!
!*        0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll
!
  USE MODD_VAR_ll, ONLY : IP, TCRRT_COMDATA, TCRRT_PROCONF, JPHALO, &
                          MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : LWEST_ll, LEAST_ll, LNORTH_ll, LSOUTH_ll
!
!*       0.    DECLARATIONS
!              ------------
!
  IMPLICIT NONE
!
!
!*       0.    DECLARATIONS
!              ------------
!
!
!*       0.1   Declarations of dummy arguments :
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  REAL, DIMENSION(:, :), INTENT(OUT) :: PRES
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  INTEGER :: IXB, IXE, IXBG, IXEG, IJB, IJE, IJBG, IJEG ! local and global displacements
!
  REAL, DIMENSION(SIZE(PRES,1),SIZE(PRES,2)) :: ZBUF
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT ! Intermediate model splitting
!
!-----------------------------------------------------------------
!
!*       1.    Get the current splitting configuration and compute displacements
!
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IP)
!
  IF (SIZE(PFIELD, 1).EQ.1) THEN
     IF (LSOUTH_ll()) THEN
        IJB  = 1
        IJBG = TZSPLIT%NYORE
     ELSE
        IJB  = 1+JPHALO
        IJBG = TZSPLIT%NYORP
     ENDIF
!
     IF (LNORTH_ll()) THEN
        IJE  = SIZE(PFIELD, 2)
        IJEG = TZSPLIT%NYENDE
     ELSE
        IJE  = SIZE(PFIELD, 2)-JPHALO
        IJEG = TZSPLIT%NYENDP
     ENDIF
!
!-----------------------------------------------------------------
!
!*       2.    Fill the intermediate buffer
!
     ZBUF = 0
     ZBUF(IJBG:IJEG, :) = PFIELD(1, IJB:IJE, :)
!
!-----------------------------------------------------------------
!
!*       2.    Merge local sums
!
     CALL MPI_ALLREDUCE(ZBUF, PRES, SIZE(PRES, 1) * SIZE(PRES, 2), MPI_PRECISION, MPI_SUM, NMNH_COMM_WORLD,KINFO)
!
!-----------------------------------------------------------------
  ELSE
     IF (SIZE(PFIELD, 2).EQ.1) THEN
        IF (LWEST_ll()) THEN
           IXB  = 1
           IXBG = TZSPLIT%NXORE
        ELSE
           IXB  = 1+JPHALO
           IXBG = TZSPLIT%NXORP
        ENDIF
!
        IF (LEAST_ll()) THEN
           IXE  = SIZE(PFIELD, 1)
           IXEG = TZSPLIT%NXENDE
        ELSE
           IXE  = SIZE(PFIELD, 1)-JPHALO
           IXEG = TZSPLIT%NXENDP
        ENDIF
!
!-----------------------------------------------------------------
!
!*       2.    Fill the intermediate buffer
!
        ZBUF = 0
        ZBUF(IXBG:IXEG, :) = PFIELD(IXB:IXE, 1, :)
!
!-----------------------------------------------------------------
!
!*       2.    Merge local sums
!
        CALL MPI_ALLREDUCE(ZBUF, PRES, SIZE(PRES, 1) * SIZE(PRES,2), MPI_PRECISION, MPI_SUM, NMNH_COMM_WORLD,KINFO)       
     ENDIF
  ENDIF
  
      END SUBROUTINE SUM_DIM1_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_0DD_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_0DD_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the scalar argument PRES on processors.
!
!     REDUCE_SUM_0Q_ll is the routine for scalar REAL*16 argument
!     of the generic routine REDUCESUM_ll.
!
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_0Q_ll, each processor
!     computes its local sum PRES; in REDUCE_SUM_0Q_ll
!     we reduce this values and return the global sum
!     in the PRES variable  REAL*16.
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       MPI_PRECISION -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!     R. Guivarch 09/07/98 Same argument PRES INOUT
!
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
 USE MODE_REPRO_SUM 
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  TYPE(DOUBLE_DOUBLE), INTENT(INOUT) :: PRES ! sum 
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  TYPE(DOUBLE_DOUBLE)  :: ZRES ! sum 
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  IF (FIRST_CALL_DD) CALL INIT_DD(KINFO)
  ZRES%R = 0.0 ;  ZRES%E = 0.0
  CALL MPI_ALLREDUCE(PRES, ZRES, 1, MNH_DOUBLE_DOUBLE , &
                     MNH_SUM_DD, NMNH_COMM_WORLD, KINFO)

  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE REDUCE_SUM_0DD_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_0D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_0D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the scalar argument PRES on processors.
!
!     REDUCE_SUM_0D_ll is the routine for scalar argument
!     of the generic routine REDUCESUM_ll.
!
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_0D_ll, each processor
!     computes its local sum PRES; in REDUCE_SUM_0D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       MPI_PRECISION -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!     R. Guivarch 09/07/98 Same argument PRES INOUT
!
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : MPI_PRECISION
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  REAL, INTENT(INOUT) :: PRES ! sum 
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  REAL :: ZRES ! Intermediate result
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  CALL MPI_ALLREDUCE(PRES, ZRES, 1, MPI_PRECISION, &
                     MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_0D_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_1DD_ll(PRES, KINFO)
!     ########################################
!

!!    Author
!!    ------
!     J.Escobar 22/10/2010
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
        USE MODE_REPRO_SUM 

  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
   TYPE(DOUBLE_DOUBLE), DIMENSION(:), INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  TYPE(DOUBLE_DOUBLE), DIMENSION(SIZE(PRES,1)) :: ZRES ! Intermediate sum

!
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  IF (FIRST_CALL_DD) CALL INIT_DD(KINFO)
  ZRES%R = 0.0 ;  ZRES%E = 0.0
  CALL MPI_ALLREDUCE(PRES, ZRES, SIZE(PRES), MNH_DOUBLE_DOUBLE , &
                     MNH_SUM_DD, NMNH_COMM_WORLD, KINFO)
PRES = ZRES
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE REDUCE_SUM_1DD_ll
!

!     ########################################
      SUBROUTINE REDUCE_SUM_1D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_1D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the each entry of the one-dimensional vector PRES
!     on all processors.
! 
!     REDUCE_SUM_1D_ll is the routine for 1D argument
!     of the generic routine REDUCESUM_ll.
! 
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_1D_ll, each processor
!     computes its local 1D sum PRES; in REDUCE_SUM_1D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       MPI_PRECISION -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!     R. Guivarch 09/07/98 Same argument PRES INOUT
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : MPI_PRECISION
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  REAL, DIMENSION(:), INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  REAL, DIMENSION(SIZE(PRES,1)) :: ZRES ! Intermediate sum
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  CALL MPI_ALLREDUCE(PRES, ZRES, SIZE(PRES,1), MPI_PRECISION, &
                     MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_1D_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_2D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_2D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the each entry of the two-dimensional vector PRES
!     on all processors.
! 
!     REDUCE_SUM_2D_ll is the routine for 2D argument
!     of the generic routine REDUCESUM_ll.
! 
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_2D_ll, each processor
!     computes its local 2D sum PRES; in REDUCE_SUM_2D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       MPI_PRECISION -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!     R. Guivarch 09/07/98 Same argument PRES INOUT
!
!-------------------------------------------------------------------------------
! 
!*       0.    DECLARATIONS 
!
  USE MODD_VAR_ll, ONLY : MPI_PRECISION
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  REAL, DIMENSION(:,:), INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
! 
  REAL, DIMENSION(SIZE(PRES,1),SIZE(PRES,2)) :: ZRES ! Intermediate sum
!
  INTEGER :: IDIM 
! 
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  IDIM = SIZE(PRES,1) * SIZE(PRES,2)
!
  CALL MPI_ALLREDUCE(PRES, ZRES, IDIM, MPI_PRECISION, MPI_SUM, &
                     NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_2D_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_3D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_3D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the each entry of the three-dimensional vector PRES
!     on all processors.
! 
!     REDUCE_SUM_3D_ll is the routine for 3D argument
!     of the generic routine REDUCESUM_ll.
! 
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_3D_ll, each processor
!     computes its local 3D sum PRES; in REDUCE_SUM_3D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       MPI_PRECISION -
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 27/06/98
!     R. Guivarch 09/07/98 Same argument PRES INOUT
!
!-------------------------------------------------------------------------------
! 
!*       0.    DECLARATIONS 
!
  USE MODD_VAR_ll, ONLY : MPI_PRECISION
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
! 
  REAL, DIMENSION(SIZE(PRES,1),SIZE(PRES,2),SIZE(PRES,3)) :: ZRES ! Intermediate
                                                                  ! sum
!
  INTEGER :: IDIM
! 
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  IDIM = SIZE(PRES,1) * SIZE(PRES,2) * SIZE(PRES,3)
!
  CALL MPI_ALLREDUCE(PRES, ZRES, IDIM, MPI_PRECISION, MPI_SUM, &
                     NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_3D_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_I0D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_I0D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the scalar argument PRES on processors.
!
!     REDUCE_SUM_I0D_ll is the routine for integer scalar argument
!     of the generic routine REDUCESUM_ll.
!
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_I0D_ll, each processor
!     computes its local sum PRES; in REDUCE_SUM_I0D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Author
!!    ------
!     D. Gazen  * L.A. *
!
!!    Modifications
!!    -------------
!     Original 4/09/2000
!
!-----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  INTEGER, INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  INTEGER :: ZRES ! Intermediate result
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  CALL MPI_ALLREDUCE(PRES, ZRES, 1, MPI_INTEGER, &
                     MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_I0D_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_I1D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_I1D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the each entry of the one-dimensional vector PRES
!     on all processors.
!
!     REDUCE_SUM_I1D_ll is the routine for integer 1D argument
!     of the generic routine REDUCESUM_ll.
!
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_I1D_ll, each processor
!     computes its local 1D sum PRES; in REDUCE_SUM_I1D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Author
!!    ------
!     D. Gazen  * L.A. *
!
!!    Modifications
!!    -------------
!     Original 4/09/2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  INTEGER, DIMENSION(:), INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  INTEGER, DIMENSION(SIZE(PRES,1)) :: ZRES ! Intermediate sum
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  CALL MPI_ALLREDUCE(PRES, ZRES, SIZE(PRES,1), MPI_INTEGER, &
                     MPI_SUM, NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_I1D_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_I2D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_2D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the each entry of the two-dimensional vector PRES
!     on all processors.
!
!     REDUCE_SUM_I2D_ll is the routine for integer 2D argument
!     of the generic routine REDUCESUM_ll.
!
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_I2D_ll, each processor
!     computes its local 2D sum PRES; in REDUCE_SUM_I2D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Author
!!    ------
!     D. Gazen  * L.A. *
!
!!    Modifications
!!    -------------
!     Original 4/09/2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  INTEGER, DIMENSION(:,:), INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  INTEGER, DIMENSION(SIZE(PRES,1),SIZE(PRES,2)) :: ZRES ! Intermediate sum
!
  INTEGER :: IDIM
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  IDIM = SIZE(PRES,1) * SIZE(PRES,2)
!
  CALL MPI_ALLREDUCE(PRES, ZRES, IDIM, MPI_INTEGER, MPI_SUM, &
                     NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_I2D_ll
!
!     ########################################
      SUBROUTINE REDUCE_SUM_I3D_ll(PRES, KINFO)
!     ########################################
!
!!****  *REDUCE_SUM_I3D_ll*-
!
!!    Purpose
!!    -------
!     This routine calculates the sum of the values
!     of the each entry of the three-dimensional vector PRES
!     on all processors.
!
!     REDUCE_SUM_I3D_ll is the routine for 3D argument
!     of the generic routine REDUCESUM_ll.
!
!!    Method
!!    ------
!     Before the call to REDUCE_SUM_I3D_ll, each processor
!     computes its local 3D sum PRES; in REDUCE_SUM_I3D_ll
!     we reduce this values and return the global sum
!     in the PRES variable.
!
!!    External
!!    --------
!
!!    Author
!!    ------
!     D. Gazen  * L.A. *
!
!!    Modifications
!!    -------------
!     Original 4/09/2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : MPI_PRECISION
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
  INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: PRES ! sum
!
  INTEGER, INTENT(OUT) :: KINFO ! MPI return status
!
!*       0.2   Declarations of local variables :
!
  INTEGER, DIMENSION(SIZE(PRES,1),SIZE(PRES,2),SIZE(PRES,3)) :: ZRES ! Intermediate
                                                                  ! sum
!
  INTEGER :: IDIM
!
!-------------------------------------------------------------------------------
!
!*       1. CALL THE MPI_ALLREDUCE ROUTINE
!           ------------------------------
!
  IDIM = SIZE(PRES,1) * SIZE(PRES,2) * SIZE(PRES,3)
!
  CALL MPI_ALLREDUCE(PRES, ZRES, IDIM, MPI_INTEGER, MPI_SUM, &
                     NMNH_COMM_WORLD, KINFO)
!
  PRES = ZRES
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REDUCE_SUM_I3D_ll
!
!-------------------------------------------------------------------------------
!
END MODULE MODE_SUM_ll
