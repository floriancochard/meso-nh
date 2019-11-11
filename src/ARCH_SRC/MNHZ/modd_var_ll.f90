!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!      ##################
       MODULE MODD_VAR_ll
!      ##################
!
!!****  *MODD_VAR_ll* - declaration of global parallel variables
!!
!!    Purpose
!!    -------
!       The purpose of this declarative module is to declare the
!       parallel global variables
!
!!   Reference
!!    ---------
!
!!    Authors
!!    -------
!
!     R. Guivarch               * CERFACS - ENSEEIHT *
!     Ph. Kloos                 * CERFACS - CNRM *
!     N. Gicquel                * CNRM *
!
!!    Implicit Arguments
!!    ------------------
!
!     MODD_STRUCTURE_ll
!        most of the structured types defined in this module
!
!!    Modifications
!!    -------------
!
!    Original 04/05/99

!-------------------------------------------------------------------------------
!
  USE MODD_STRUCTURE_ll
! 
  IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
! Current configuration for current model
!
  TYPE(PROCONF_ll), POINTER       :: TCRRT_PROCONF
!
! Current communication data structure for current model
! and local processor
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TCRRT_COMDATA
!
! Number of the local processor
!
  INTEGER :: IP
!
! Number of processors
!
  INTEGER :: NPROC
!
! 'BSPLITTING' or 'XSPLITTING' or 'YSPLITTING'
!
  CHARACTER(LEN=10) :: YSPLITTING = "BSPLITTING"
!
!JUAN
INTEGER,SAVE      :: NZ_PROC = 1  ! Number of proc to use in the Z splitting 
!JUAN
!
! i/o output filename & unit
!
  CHARACTER(LEN=*),PARAMETER :: YOUTPUTFILE="surcouche_listing"
!
  INTEGER :: NIOUNIT
!
! mpi communicators
!
  INTEGER :: NHALO_COM
  INTEGER :: NTRANS_COM
  INTEGER :: NHALO2_COM
  INTEGER :: NGRID_COM
!
! maximum size of a halo
!
  INTEGER :: NMAXSIZEHALO
!
! width of an interior halo ; should be equal to JPHEXT
!
  INTEGER :: JPHALO 
!
! dimensions of the extended domain
! DIMX = NIMAX + 2*JPHEXT ...
!
  INTEGER :: DIMX,DIMY,DIMZ
!
! MPI_PRECISION, MPI_2PRECISION
!
  INTEGER :: MPI_PRECISION
  INTEGER :: MPI_2PRECISION
!
  INTEGER, PARAMETER :: NTMAX = 100
!
! maximum size of 2D and 3D buffers
!
  INTEGER :: NBUFFERSIZE_2D
  INTEGER :: NBUFFERSIZE_3D
!
! buffer sizes
!
  INTEGER  :: NCOMBUFFSIZE, NCOMBUFFSIZE1, NCOMBUFFSIZE2
!
! variable to define message tag
!
  INTEGER, PARAMETER :: NMAXTAG = 5000
  INTEGER, PARAMETER :: NNEXTTAG = 50
!
  INTEGER, PARAMETER :: NMODULO_MSSGTAG = 10
!
END MODULE MODD_VAR_ll
