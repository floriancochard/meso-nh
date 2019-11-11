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
!-----------------------------------------------------------------

!     ###################
      MODULE MODE_SUM2_ll
!     ###################
!
!!    Purpose
!!    -------
!
!     The purpose of this module is to provide functions for the summations
!     of the budget computations, for the search of minimum and maximum.
!
!!    Functions Of The User Interface
!!    ------------------------------
!
!     FUNCTIONS : GMAXLOC_ll (GMAXLOC3D_ll, GMAXLOC2D_ll, GMAXLOC1D_ll),
!                 GMINLOC_ll (GMINLOC3D_ll, GMINLOC2D_ll, GMINLOC1D_ll)
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
!
!     ###################################################
      FUNCTION GMAXLOC3D_ll(PARRAY, MASK) RESULT(KMAXLOC)
!     ###################################################
!
!!****  *GMAXLOC3D_ll*-
!
!!    Purpose
!!    -------
!     This function calculates the summations of the budget computations
!     for the search of maximum im a three-dimensionnal array
!
!     GMAXLOC3D_ll is the function for 3D argument
!     of the generic function GMAXLOC_ll
!
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
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
!*       0.    DECLARATIONS
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
!
  IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!
  REAL, DIMENSION(:,:,:), INTENT(IN)    :: PARRAY   ! local value
  LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: MASK ! logical mask 
  INTEGER, DIMENSION(3)                 :: KMAXLOC  ! indices
                                                  ! of the
                                                  ! maximum value
                                                  ! on the whole domain
                                                  ! (global coordinates)
!
!*       0.2   Declarations of local variables
!
  INTEGER, DIMENSION(3) :: IMAXLOC
  REAL :: ZMAX
  INTEGER, DIMENSION(3) :: IDIMS
!juan     : local physical domaine
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !
!
!-------------------------------------------------------------------------------
!
!*       1.    CALCULATE THE VALUE AND LOCATION OF THE MAXIMUM
!              ON THE CURRENT PROCESS
!              -----------------------------------------------
!
!
!juan : local physical domaine
!
IIB = 1              + JPHEXT
IIE = SIZE(PARRAY,1) - JPHEXT
IJB = 1              + JPHEXT
IJE = SIZE(PARRAY,2) - JPHEXT
IKB = 1              + JPVEXT
IKE = SIZE(PARRAY,3) - JPHEXT
!
IF (PRESENT(MASK)) THEN
IMAXLOC = MAXLOC(PARRAY(IIB:IIE,IJB:IJE,IKB:IKE), MASK(IIB:IIE,IJB:IJE,IKB:IKE))
ZMAX    = MAXVAL(PARRAY(IIB:IIE,IJB:IJE,IKB:IKE), MASK=MASK(IIB:IIE,IJB:IJE,IKB:IKE))
ELSE
IMAXLOC = MAXLOC(PARRAY(IIB:IIE,IJB:IJE,IKB:IKE))
ZMAX    = MAXVAL(PARRAY(IIB:IIE,IJB:IJE,IKB:IKE))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    FIND OUT ON WHICH PROCESS THE MAXIMUM IS
!
  IDIMS=(/1,2,3/)
  KMAXLOC = UTIL_GMAXLOC_ll(ZMAX, IMAXLOC, IDIMS)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION GMAXLOC3D_ll
!
!     ###########################################################
      FUNCTION GMAXLOC2D_ll(PARRAY, KDIMS, MASK ) RESULT(KMAXLOC)
!     ###########################################################
!
!!****  *GMAXLOC2D_ll*-
!
!!    Purpose
!!    -------
!     This function calculates the summations of the budget computations
!     for the search of maximum im a two-dimensionnal array
!
!     GMAXLOC2D_ll is the function for 2D argument
!     of the generic function GMAXLOC_ll
!
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
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
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!
  REAL, DIMENSION(:,:), INTENT(IN)    :: PARRAY   ! local value
  LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK ! logical mask
  INTEGER, DIMENSION(2)                 :: KMAXLOC  ! indices
                                                    ! of the
                                                    ! maximum value
                                                    ! on the whole domain
                                                    ! (global coordinates)
  INTEGER, DIMENSION(2), OPTIONAL :: KDIMS
!
!*       0.2   Declarations of local variables
!
  INTEGER, DIMENSION(2) :: IMAXLOC
  INTEGER, DIMENSION(2) :: IDIMS
  REAL :: ZMAX
!
!-------------------------------------------------------------------------------
!
!*       1.  SET DIRECTIONS IN WHICH THE MAXIMUM HAS TO BE CALCULATED
!            (DEFAULT : x,y)
!            --------------------------------------------------------
!
  IF (PRESENT(KDIMS)) THEN
    IDIMS=KDIMS
  ELSE
    IDIMS = (/1,2/)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.  CALCULATE THE VALUE AND LOCATION OF THE MAXIMUM
!            ON THE CURRENT PROCESS
!            -----------------------------------------------
!
  IF (PRESENT(MASK)) THEN
    IMAXLOC=MAXLOC(PARRAY, MASK)
    ZMAX = MAXVAL(PARRAY, MASK=MASK)
  ELSE
    IMAXLOC=MAXLOC(PARRAY)
    ZMAX = MAXVAL(PARRAY)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.  FIND OUT ON WHICH PROCESS THE MAXIMUM IS
!            ----------------------------------------
!
  KMAXLOC = UTIL_GMAXLOC_ll(ZMAX, IMAXLOC, KDIMS)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION GMAXLOC2D_ll
!
!     ##########################################################
      FUNCTION GMAXLOC1D_ll(PARRAY, KDIMS, MASK) RESULT(KMAXLOC)
!     ##########################################################
!
!!****  *GMAXLOC1D_ll*-
!
!!    Purpose
!!    -------
!     This function calculates the summations of the budget computations
!     for the search of maximum im a one-dimensionnal array
!
!     GMAXLOC1D_ll is the function for 1D argument
!     of the generic function GMAXLOC_ll
!
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
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
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!
  REAL, DIMENSION(:), INTENT(IN)              :: PARRAY ! input array
                                                        ! in which the maximum
                                                        ! is to be found
  LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: MASK     ! mask
  INTEGER                                     :: KMAXLOC  ! indice
                                                  ! of the
                                                  ! maximum value
                                                  ! on the whole domain
                                                  ! (global coordinates)
  INTEGER, OPTIONAL                           :: KDIMS
!
!*       0.2   Declarations of local variables
!
  INTEGER, DIMENSION(1) :: IMAXLOC, IDIMS, IGMAXLOC
  REAL    :: ZMAX
!
!-------------------------------------------------------------------------------
!
!*       1.  SET DIRECTIONS IN WHICH THE MAXIMUM HAS TO BE CALCULATED
!            (DEFAULT : x)
!            --------------------------------------------------------
!
  IF (PRESENT(KDIMS)) THEN
    IDIMS(1)=KDIMS
  ELSE
    IDIMS(1)=1 
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.  CALCULATE THE VALUE AND LOCATION OF THE MAXIMUM
!            ON THE CURRENT PROCESS
!            -----------------------------------------------
!
  IF (PRESENT(MASK)) THEN
    IMAXLOC=MAXLOC(PARRAY, MASK)
    ZMAX = MAXVAL(PARRAY, MASK=MASK)
  ELSE
    IMAXLOC=MAXLOC(PARRAY)
    ZMAX = MAXVAL(PARRAY)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.  FIND OUT ON WHICH PROCESS THE MAXIMUM IS
!            ----------------------------------------
!
!
  IGMAXLOC = UTIL_GMAXLOC_ll(ZMAX, IMAXLOC, IDIMS)
  KMAXLOC = IMAXLOC(1)
!
!-------------------------------------------------------------------------------
!
      END FUNCTION GMAXLOC1D_ll
!
!
!     ###################################################
      FUNCTION GMINLOC3D_ll(PARRAY, MASK) RESULT(KMINLOC)
!     ###################################################
!
!!****  *GMINLOC3D_ll*-
!
!!    Purpose
!!    -------
!     This function calculates the summations of the budget computations
!     for the search of minimum im a three-dimensionnal array
!
!     GMINLOC3D_ll is the function for 3D argument
!     of the generic function GMINLOC_ll
!
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
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
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!
  REAL, DIMENSION(:,:,:), INTENT(IN)    :: PARRAY   ! local value
  LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: MASK ! logical mask 
  INTEGER, DIMENSION(3)                 :: KMINLOC  ! indices
                                                    ! of the
                                                    ! maximum value
                                                    ! on the whole domain
                                                    ! (global coordinates)
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.  CALCULATE THE GLOBAL POSITION OF THE MINIMUM
!            USING THE GMAXLOC_ll FUNCTION
!            --------------------------------------------
!
  IF (PRESENT(MASK)) THEN
  !  KMINLOC = GMAXLOC_ll(-PARRAY, MASK)
    KMINLOC = GMAXLOC3D_ll(-PARRAY, MASK)
  ELSE
  !  KMINLOC = GMAXLOC_ll(-PARRAY)
    KMINLOC = GMAXLOC3D_ll(-PARRAY)
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION GMINLOC3D_ll
!
!     ###########################################################
      FUNCTION GMINLOC2D_ll(PARRAY, KDIMS, MASK ) RESULT(KMINLOC)
!     ###########################################################
!
!!****  *GMINLOC2D_ll*-
!
!!    Purpose
!!    -------
!     This function calculates the summations of the budget computations
!     for the search of minimum im a two-dimensionnal array
!
!     GMINLOC2D_ll is the function for 2D argument
!     of the generic function GMINLOC_ll
!
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
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
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!
  REAL, DIMENSION(:,:), INTENT(IN)    :: PARRAY   ! local value
  LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK ! logical mask
  INTEGER, DIMENSION(2)                 :: KMINLOC  ! indices
                                                    ! of the
                                                    ! maximum value
                                                    ! on the whole domain
                                                    ! (global coordinates)
  INTEGER, DIMENSION(2), OPTIONAL :: KDIMS
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.  CALCULATE THE GLOBAL POSITION OF THE MINIMUM
!            USING THE GMAXLOC_ll FUNCTION
!            --------------------------------------------
!
  IF (PRESENT(MASK) .AND. .NOT.PRESENT(KDIMS)) THEN
  !  KMINLOC = GMAXLOC_ll(-PARRAY, MASK=MASK)
    KMINLOC = GMAXLOC2D_ll(-PARRAY, MASK=MASK)
  ELSEIF(PRESENT(KDIMS) .AND. .NOT.PRESENT(MASK)) THEN
  !  KMINLOC = GMAXLOC_ll(-PARRAY, KDIMS=KDIMS)
    KMINLOC = GMAXLOC2D_ll(-PARRAY, KDIMS=KDIMS)
  ELSEIF(PRESENT(KDIMS) .AND. PRESENT(MASK)) THEN
  !  KMINLOC = GMAXLOC_ll(-PARRAY, KDIMS, MASK)
    KMINLOC = GMAXLOC2D_ll(-PARRAY, KDIMS, MASK)
  ELSE
  !  KMINLOC = GMAXLOC_ll(-PARRAY)
    KMINLOC = GMAXLOC2D_ll(-PARRAY)
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION GMINLOC2D_ll
!
!     ##########################################################
      FUNCTION GMINLOC1D_ll(PARRAY, KDIMS, MASK) RESULT(KMINLOC)
!     ##########################################################
!
!!****  *GMINLOC1D_ll*-
!
!!    Purpose
!!    -------
!     This function calculates the summations of the budget computations
!     for the search of minimum im a one-dimensionnal array
!
!     GMINLOC1D_ll is the function for 1D argument
!     of the generic function GMINLOC_ll
!
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
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
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!
  REAL, DIMENSION(:), INTENT(IN)              :: PARRAY   ! input array
                                                          ! in which the maximum
                                                          ! is to be found
  LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: MASK     ! mask
  INTEGER                                     :: KMINLOC  ! indice
                                                  ! of the
                                                  ! maximum value
                                                  ! on the whole domain
                                                  ! (global coordinates)
  INTEGER, OPTIONAL                           :: KDIMS
!
!*       0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.  CALCULATE THE GLOBAL POSITION OF THE MINIMUM
!            USING THE GMAXLOC_ll FUNCTION
!            --------------------------------------------
!
  IF (PRESENT(MASK) .AND. .NOT.PRESENT(KDIMS)) THEN
  !  KMINLOC = GMAXLOC_ll(-PARRAY, MASK=MASK)
    KMINLOC = GMAXLOC1D_ll(-PARRAY, MASK=MASK)
  ELSEIF(PRESENT(KDIMS) .AND. .NOT.PRESENT(MASK)) THEN
  !  KMINLOC = GMAXLOC_ll(-PARRAY, KDIMS=KDIMS)
    KMINLOC = GMAXLOC1D_ll(-PARRAY, KDIMS=KDIMS)
  ELSEIF(PRESENT(KDIMS) .AND. PRESENT(MASK)) THEN
  !  KMINLOC = GMAXLOC_ll(-PARRAY, KDIMS, MASK)
    KMINLOC = GMAXLOC1D_ll(-PARRAY, KDIMS, MASK)
  ELSE
  !  KMINLOC = GMAXLOC_ll(-PARRAY)
    KMINLOC = GMAXLOC1D_ll(-PARRAY)
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION GMINLOC1D_ll
!
!     #############################################################
      FUNCTION LMAXLOC_ll(PVALUE, KLOCALMAX, KPROC) RESULT(KMAXLOC)
!     #############################################################
!
!!****  *LMAXLOC_ll* - returns the local coordinates of the maximum value
!                      of a set of real numbers distributed on the procs
!
!!    PURPOSE
!!    -------
!       This function takes as arguments a real number and its local position
!     on the subdomain. It returns the local position of the maximum of
!     the real numbers on all the subdomains, and the number of the subdomain
!     corresponding to the maximum.
!
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_VAR_ll
!       IP -
!       MPI_2PRECISION
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
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : IP, MPI_2PRECISION
!
  IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments and result
!
  REAL, INTENT(IN)                    :: PVALUE     ! local value
  INTEGER, DIMENSION(:), INTENT(IN)   :: KLOCALMAX  ! indices of the max. 
                                                    ! on the local proc.
  INTEGER, DIMENSION(SIZE(KLOCALMAX)) :: KMAXLOC    ! indices
                                                    ! of the
                                                    ! maximum value
                                                    ! on the whole domain
                                                    ! (local coordinates)
  INTEGER, INTENT(OUT)                :: KPROC      ! number of the process
                                                    ! with the max value
!
!*       0.2   Declarations of local variables
!
  INTEGER :: IPROCMAX
  INTEGER :: ISIZE, INFO_ll
  INTEGER, DIMENSION(SIZE(KLOCALMAX)) :: IMAXLOC
  INTEGER :: IERR ! return status
  REAL, DIMENSION(2) :: ZBUFIN, ZBUFOUT
!
!-------------------------------------------------------------------------------
!
!*       1.    REDUCTION
!              ---------
!
  ZBUFIN (1) = PVALUE
  ZBUFIN (2) = IP
  CALL MPI_ALLREDUCE(ZBUFIN, ZBUFOUT, 1, MPI_2PRECISION, MPI_MAXLOC, &
                     NMNH_COMM_WORLD, IERR)
!
!
!-------------------------------------------------------------------------------
!
!*       2.  BROADCAST THE COORDINATES OF THE MAXIMUM VALUE
!            ----------------------------------------------
!
  IPROCMAX = ZBUFOUT(2)
  IMAXLOC = KLOCALMAX
  ISIZE=SIZE(KLOCALMAX)
  CALL MPI_BCAST(IMAXLOC, ISIZE, MPI_INTEGER, IPROCMAX-1, &
                 NMNH_COMM_WORLD, INFO_ll)
!
  KPROC=IPROCMAX
  KMAXLOC=IMAXLOC
!
!-------------------------------------------------------------------------------
!
      END FUNCTION LMAXLOC_ll 
!
!     ################################################################
      FUNCTION UTIL_GMAXLOC_ll(PVALUE, KLOCALMAX, KDIM) RESULT(KMAXLOC)
!     ################################################################
!
!!****  *GMAXLOC_ll* - returns the global coordinates of the maximum value
!                      of a set of real numbers distributed on the procs
!
!!    Purpose
!!    -------
!     This function takes as arguments a real number and its local position
!     on the subdomain. It returns the global position of the maximum of
!     the real numbers on all the subdomains, and the number of the subdomain
!     corresponding to the maximum.
!     The KDIM array can be used to specify which dimensions are used.
! 
!!    Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       type MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!       TCRRT_PROCONF - Current configuration for current model
!       IP -
!       MPI_2PRECISION
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
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : IP, MPI_2PRECISION, TCRRT_PROCONF
!
  USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll
!
  IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!              ------------------------------------
!
  REAL, INTENT(IN)                    :: PVALUE     ! local value
  INTEGER, DIMENSION(:), INTENT(IN)   :: KLOCALMAX  ! indices of the max.
                                                    ! on the local proc.
  INTEGER, DIMENSION(SIZE(KLOCALMAX)) :: KMAXLOC    ! indices
                                                    ! of the
                                                    ! maximum value
                                                    ! on the whole domain
                                                    ! (global coordinates)
!
  INTEGER, DIMENSION(:), OPTIONAL :: KDIM ! used dimensions 
                                         ! (e.g. (2,3) for Y and Z dimensions)
!
!*       0.2   Declarations of local variables
!
  INTEGER :: INFO_ll             ! return status
  INTEGER :: ISIZE               ! number of coordinates
  INTEGER :: IPROCMAX            ! number of the proc. with the maximum value
  INTEGER, DIMENSION(3) :: IBEGIN ! coordinates of the subdomain
  INTEGER, DIMENSION(SIZE(KLOCALMAX)) :: IMAXLOC
  REAL, DIMENSION(2) :: ZBUFIN, ZBUFOUT
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT ! to pick up the coordinates
                                              ! of the subdomain 
!
!-------------------------------------------------------------------------------
!
!*       1.  CALCULATE THE RANK OF THE PROCESS CONTAINING THE MAXIMUM VALUE
!            --------------------------------------------------------------
!
  ISIZE = SIZE(KLOCALMAX)
  ZBUFIN (1) = PVALUE
  ZBUFIN (2) = IP
! 
  CALL MPI_ALLREDUCE(ZBUFIN, ZBUFOUT, 1, MPI_2PRECISION, MPI_MAXLOC, &
                     NMNH_COMM_WORLD, INFO_ll)
!
!-------------------------------------------------------------------------------
!
!*       2.  BROADCAST THE COORDINATES OF THE MAXIMUM VALUE
!            ----------------------------------------------
!
  IPROCMAX = ZBUFOUT(2)
  IMAXLOC = KLOCALMAX
  CALL MPI_BCAST(IMAXLOC, ISIZE, MPI_INTEGER, IPROCMAX-1, &
                 NMNH_COMM_WORLD, INFO_ll)
!
!-------------------------------------------------------------------------------
!
!*       3.  CONVERT GLOBAL COORDINATES INTO LOCAL COORDINATES
!            -------------------------------------------------
!
! we point to the structure which contains information
! on the subdomain containing
! the maximum point
!
  TZSPLIT => TCRRT_PROCONF%TSPLITS_B(IPROCMAX)
!
  IBEGIN(1) = TZSPLIT%NXORE - 1
  IBEGIN(2) = TZSPLIT%NYORE - 1
  IBEGIN(3) = 0
!
  KMAXLOC(:) = IMAXLOC(:) + IBEGIN(KDIM(:))
!
!-------------------------------------------------------------------------------
!
      END FUNCTION UTIL_GMAXLOC_ll
!
END MODULE MODE_SUM2_ll
