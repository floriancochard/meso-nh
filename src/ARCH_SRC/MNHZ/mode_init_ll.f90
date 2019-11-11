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
      MODULE MODE_INIT_ll
!     ###################
!!
!!    Purpose
!!    -------
!
!     The purpose of this module is the implementation of the initialisation
!     of parallel data structure
!
!!    Routines Of The User Interface
!!    ------------------------------
! 
!     SUBROUTINES : SET_SPLITTING_ll, SET_LBX_ll, SET_LBY_ll
!                   SET_DIM_ll, SET_JP_ll, SET_XRATIO_ll, SET_YRATIO_ll
!                   SET_DAD_ll, SET_XOR_ll, SET_XEND_ll, SET_YOR_ll,
!                   SET_YEND_ll, SET_DAD0_ll,
! 
!                   INI_PARA_ll, END_PARA_ll
!                   
!!    Implicit Arguments
!!    ------------------
! 
!!    Reference
!!    ---------
! 
!      User Interface for Meso-NH parallel package
!      Ph. Kloos, L. Giraud, R. Guivarch, D. Lugato
!
!!    Authors
!!    -------
! 
!     R. Guivarch               * CERFACS - ENSEEIHT *
!     Ph. Kloos                 * CERFACS - CNRM *
!     N. Gicquel                * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!
!       Original     May 19, 1998
!       Juan     19/08/2005: distinction Halo NORD/SUD & EST/WEST
!
!-------------------------------------------------------------------------------
!
  USE MODD_MPIF
!
  IMPLICIT NONE
!
!  INCLUDE 'mpif.h'
!
  CONTAINS
!
!     #######################################
      SUBROUTINE SET_SPLITTING_ll(HSPLITTING)
!     #######################################
!
!!****  *SET_SPLITTING_ll* - 
!
!!    Purpose
!!    -------
!       Set the variable YSPLITTING with HSPLITTING
!
!-------------------------------------------------------------------------------
!
  USE MODD_VAR_ll, ONLY : YSPLITTING
!
  IMPLICIT NONE
!
  CHARACTER(LEN=10) :: HSPLITTING
!
!-------------------------------------------------------------------------------
!
  YSPLITTING = HSPLITTING
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_SPLITTING_ll
!
!     ################################
      SUBROUTINE SET_LBX_ll(KLBX, KMI)
!     ################################
!
!!****  *SET_LBX_ll *-
!
!!    Purpose
!!    -------
!       Set the variable CLBCX(KMI,:) with KLBX
!
!-------------------------------------------------------------------------------
!
  USE MODD_PARAMETERS_ll, ONLY : JPMODELMAX
  USE MODD_DIM_ll, ONLY : CLBCX
!
  IMPLICIT NONE
!
  CHARACTER(LEN=*) :: KLBX
  INTEGER :: KMI
!
!-------------------------------------------------------------------------------
!
  IF (KMI.LE.JPMODELMAX) THEN
    CLBCX(KMI, :) = KLBX
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LBX_ll
!
!     ################################
      SUBROUTINE SET_LBY_ll(KLBY, KMI)
!     ################################
!
!!****  *SET_LBY_ll *-
!
!!    Purpose
!!    -------
!       Set the variable CLBCY(KMI,:) with KLBY
!
!-------------------------------------------------------------------------------
!
  USE MODD_PARAMETERS_ll, ONLY : JPMODELMAX
  USE MODD_DIM_ll, ONLY : CLBCY
!
  IMPLICIT NONE
!
  CHARACTER(LEN=*) :: KLBY
  INTEGER :: KMI
!
!-------------------------------------------------------------------------------
!
  IF (KMI.LE.JPMODELMAX) THEN
    CLBCY(KMI, :) = KLBY
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LBY_ll
!
!     #################################
      SUBROUTINE SET_DIM_ll(KX, KY, KZ)
!     #################################
!
!!****  *SET_DIM_ll *-
!
!!    Purpose
!!    -------
!       Set the variable CLBCY(KMI,:) with KLBY
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NIMAX_ll, NJMAX_ll, NKMAX_ll
!
  IMPLICIT NONE
!
  INTEGER :: KX,KY,KZ
!
!-------------------------------------------------------------------------------
!
  NIMAX_ll = KX
  NJMAX_ll = KY
  NKMAX_ll = KZ
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_DIM_ll
!
!     #####################################################
      SUBROUTINE SET_JP_ll(KMODELMAX, KHEXT, KVEXT, KPHALO)
!     #####################################################
!
!!****  *SET_JP_ll *-
!
!!    Purpose
!!    -------
!       Set the halo variables and alloacte arrays of MODD_DIM_ll
!
!-------------------------------------------------------------------------------
!
  USE MODD_PARAMETERS_ll, ONLY : JPMODELMAX, JPHEXT, JPVEXT
  USE MODD_DIM_ll, ONLY : NDXRATIO_ALL, NDYRATIO_ALL, &
                          NXOR_ALL, NYOR_ALL, NXEND_ALL, NYEND_ALL, &
                          NDAD, CLBCX, CLBCY
  USE MODD_VAR_ll, ONLY : JPHALO
!
  IMPLICIT NONE
!
  INTEGER :: KMODELMAX, KHEXT, KVEXT, KPHALO
!
!-------------------------------------------------------------------------------
!
  JPMODELMAX = KMODELMAX
  JPHEXT = KHEXT
  JPVEXT = KVEXT 
  JPHALO = KPHALO
!
! Allocate arrays declared in MODD_DIM_ll
!
  ALLOCATE(NDXRATIO_ALL(JPMODELMAX), NDYRATIO_ALL(JPMODELMAX))
  ALLOCATE(NXOR_ALL(JPMODELMAX), NYOR_ALL(JPMODELMAX)) 
  ALLOCATE(NXEND_ALL(JPMODELMAX), NYEND_ALL(JPMODELMAX)) 
  ALLOCATE(NDAD(JPMODELMAX)) 
  ALLOCATE(CLBCX(JPMODELMAX, 2), CLBCY(JPMODELMAX, 2)) 
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_JP_ll
!
!     ######################################
      SUBROUTINE SET_XRATIO_ll(KXRATIO, KMI)
!     ######################################
!
!!****  *SET_XRATIO_ll *-
!
!!    Purpose
!!    -------
!       Set the variable NDXRATIO_ALL(KMI) with KXRATIO
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NDXRATIO_ALL
!
  IMPLICIT NONE
!
  INTEGER :: KXRATIO, KMI
!
!-------------------------------------------------------------------------------
!
  NDXRATIO_ALL(KMI) = KXRATIO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_XRATIO_ll
!
!     ######################################
      SUBROUTINE SET_YRATIO_ll(KYRATIO, KMI)
!     ######################################
!
!!****  *SET_YRATIO_ll *-
!
!!    Purpose
!!    -------
!       Set the variable NDYRATIO_ALL(KMI) with KYRATIO
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NDYRATIO_ALL
!
  IMPLICIT NONE
!
  INTEGER :: KYRATIO, KMI
!
!-------------------------------------------------------------------------------
!
  NDYRATIO_ALL(KMI) = KYRATIO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_YRATIO_ll
!
!     ################################
      SUBROUTINE SET_DAD_ll(KDAD, KMI)
!     ################################
!
!!****  *SET_DAD_ll* -
!
!!    Purpose
!!    -------
!       Set the variable NDAD(KMI) with KDAD
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NDAD
!
  IMPLICIT NONE
!
  INTEGER :: KDAD, KMI
!
!-------------------------------------------------------------------------------
!
  NDAD(KMI) = KDAD
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_DAD_ll
!
!     ################################
      SUBROUTINE SET_XOR_ll(KXOR, KMI)
!     ################################
!
!!****  *SET_XOR_ll* -
!
!!    Purpose
!!    -------
!       Set the variable NXOR_ALL(KMI) with KXOR
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NXOR_ALL
!
  IMPLICIT NONE
!
  INTEGER :: KXOR, KMI
!
!-------------------------------------------------------------------------------
!
  NXOR_ALL(KMI) = KXOR
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_XOR_ll
!
!     ##################################
      SUBROUTINE SET_XEND_ll(KXEND, KMI)
!     ##################################
!
!!****  *SET_XEND_ll* -
!
!!    Purpose
!!    -------
!       Set the variable NXEND_ALL(KMI) with KXEND
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NXEND_ALL
!
  IMPLICIT NONE
!
  INTEGER :: KXEND, KMI
!
!-------------------------------------------------------------------------------
!
  NXEND_ALL(KMI) = KXEND
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_XEND_ll
!
!     ################################
      SUBROUTINE SET_YOR_ll(KYOR, KMI)
!     ################################
!
!!****  *SET_YOR_ll* -
!
!!    Purpose
!!    -------
!       Set the variable NYOR_ALL(KMI) with KYOR
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NYOR_ALL
!
  IMPLICIT NONE
!
  INTEGER :: KYOR, KMI
!
!-------------------------------------------------------------------------------
!
  NYOR_ALL(KMI) = KYOR
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_YOR_ll
!
!     ##################################
      SUBROUTINE SET_YEND_ll(KYEND, KMI)
!     ##################################
!
!!****  *SET_YEND_ll* -
!
!!    Purpose
!!    -------
!       Set the variable NYEND_ALL(KMI) with KYEND
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NYEND_ALL
!
  IMPLICIT NONE
!
  INTEGER :: KYEND, KMI
!
!-------------------------------------------------------------------------------
!
  NYEND_ALL(KMI) = KYEND
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_YEND_ll
!
!     ########################
      SUBROUTINE SET_DAD0_ll()
!     ########################
!
!!****  *SET_DAD0_ll* -
!
!!    Purpose
!!    -------
!       fill the array NDAD with 0
!
!-------------------------------------------------------------------------------
!
  USE MODD_DIM_ll, ONLY : NDAD
!
  IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
  NDAD(:) = 0
!
      END SUBROUTINE SET_DAD0_ll
!
!     ################################
      SUBROUTINE INI_PARA_ll(KINFO_ll)
!     ################################
!
!!****  *INI_PARA_ll* - routine to initialize the parallel variables
!!
!!    Purpose
!!    -------
!     the purpose of the routine is to fill the structured type variables
!     TCRRT_PROCONF and TCRRT_COMDATA
!
!!**  Method
!!    ------
!!
!!    External
!!    --------
!     Module MODE_SPLITTING_ll
!      SPLIT2
!
!     Module MODE_CONSTRUCT_ll
!       INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS,
!       CONSTRUCT_HALO1, CONSTRUCT_HALO2,
!       CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY,
!       COMPUTE_HALO_MAX, COMPUTE_TRANS_MAX
!
!     Module MODE_NEST_ll
!       INI_CHILD
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_DIM_ll
!       JPHEXT - Horizontal External points number
!       NDXRATIO_ALL, NDYRATIO_ALL, NXOR_ALL, NYOR_ALL,
!       NXEND_ALL, NYEND_ALL,...
! 
!     Module MODD_PARALLEL
!       TCRRT_PROCONF - Current configuration for current model
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
! 
!     Reference
!!    ---------
! 
!!    AUTHOR
!!    ------
!       R. Guivarch
! 
!!    MODIFICATIONS
!!    -------------
!     Original 01/05/98
!     R. Guivarch 01/01/98  Grid-Nesting
!     R. Guivarch 29/11/99  x and y splitting -> YSPLITTING 
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll
  USE MODD_PARAMETERS_ll
  USE MODD_STRUCTURE_ll
  USE MODD_VAR_ll
!
  USE MODE_SPLITTING_ll, ONLY : SPLIT2
!
  USE MODE_CONSTRUCT_ll, ONLY : INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS, &
                                CONSTRUCT_HALO1, CONSTRUCT_HALO2, &
                                CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY, &
                                COMPUTE_HALO_MAX, COMPUTE_TRANS_MAX
!
  USE MODE_NEST_ll, ONLY : INI_CHILD
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(OUT) :: KINFO_ll
!
!*       0.2   declarations of local variables
!
#if defined(MNH_MPI_BSEND)
  INTEGER  ,PARAMETER                      :: MPI_BUFFER_SIZE = 40000000
  CHARACTER,SAVE,ALLOCATABLE,DIMENSION(:)  :: MPI_BUFFER
!JUAN
  LOGICAL,SAVE                             :: GFIRSTCALL = .TRUE.
!JUAN
#endif

  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP ! intermediate zone
!
  TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
  TYPE(PROCONF_ll), POINTER :: TZPROCONF
  INTEGER :: JMODEL
  INTEGER     :: IRESP
  LOGICAL     :: GISINIT
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZE MPI :
!              --------------
!
  KINFO_ll = 0
  CALL MPI_INITIALIZED(GISINIT, KINFO_ll)
  IF (.NOT. GISINIT) THEN
     CALL MPI_INIT(KINFO_ll)
  END IF
!
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, IP, KINFO_ll)
!
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, KINFO_ll)
!
  CALL MPI_COMM_DUP(MPI_COMM_WORLD, NHALO_COM, KINFO_ll)
!
  CALL MPI_COMM_DUP(MPI_COMM_WORLD, NHALO2_COM, KINFO_ll)
!
  CALL MPI_COMM_DUP(MPI_COMM_WORLD, NTRANS_COM, KINFO_ll)
!
  CALL MPI_COMM_DUP(MPI_COMM_WORLD, NGRID_COM, KINFO_ll)
!
  IP = IP + 1
!
  MPI_PRECISION = MPI_DOUBLE_PRECISION
  MPI_2PRECISION = MPI_2DOUBLE_PRECISION
!
!-------------------------------------------------------------------------------
!
!*       2.    SET OUTPUT FILE :
!              ---------------

!  CALL OPEN_ll(UNIT=NIOUNIT,FILE=YOUTPUTFILE,ACTION='write',form&
!       &='FORMATTED',MODE=SPECIFIC,IOSTAT=IRESP)
!
!-------------------------------------------------------------------------------
!
!*       3.    ALLOCATION :
!              ----------
!
#if defined(MNH_MPI_BSEND)
  IF (GFIRSTCALL) THEN
  ALLOCATE(MPI_BUFFER(MPI_BUFFER_SIZE))
  CALL MPI_BUFFER_ATTACH(MPI_BUFFER,MPI_BUFFER_SIZE,KINFO_ll)
  GFIRSTCALL = .FALSE.
  ENDIF
#endif

  ALLOCATE(TZDZP(NPROC))
!
  ALLOCATE(TCRRT_PROCONF)
  CALL ALLOC(TCRRT_COMDATA)
  ALLOCATE(TCRRT_PROCONF%TSPLITS_B(NPROC))
  ALLOCATE(TCRRT_PROCONF%TSPLITS_X(NPROC))
  ALLOCATE(TCRRT_PROCONF%TSPLITS_Y(NPROC))
  ALLOCATE(TCRRT_PROCONF%TBOUND(NPROC))
  NULLIFY(TCRRT_PROCONF%TPARENT)
  NULLIFY(TCRRT_COMDATA%TPARENT)
  NULLIFY(TCRRT_PROCONF%TCHILDREN)
  NULLIFY(TCRRT_COMDATA%TCHILDREN)
!
!-------------------------------------------------------------------------------
!
!*       4.    SPLITTING OF THE DOMAIN :
!              -----------------------
!
  DIMX = NIMAX_ll + 2*JPHEXT
  DIMY = NJMAX_ll + 2*JPHEXT
  DIMZ = NKMAX_ll + 2*JPVEXT
!
  TCRRT_PROCONF%NUMBER = 1
!
  CALL SPLIT2(NIMAX_ll,NJMAX_ll,NKMAX_ll,NPROC,TZDZP,YSPLITTING)
!    
!-------------------------------------------------------------------------------
!
!*       5.    INITIALIZATION OF TCRRT_PROCONF :
!              -------------------------------
!
  CALL INI_PZ(TCRRT_PROCONF,TZDZP)
!
  CALL INI_BOUNDARIES(TCRRT_PROCONF)
!
  CALL INI_EZ(TCRRT_PROCONF)
!
  CALL INI_TRANS(TCRRT_PROCONF)
!
!-------------------------------------------------------------------------------
!
!*       6.    INITIALIZATION OF TCRRT_COMDATA :
!              -------------------------------
!
!*       6.1    Model Number
!
  TCRRT_COMDATA%NUMBER = 1
!
!*       6.2    Pointer from TCRRT_COMDATA to TCRRT_PROCONF for 2Way splitting
!
  TCRRT_COMDATA%TSPLIT_B => TCRRT_PROCONF%TSPLITS_B(IP)

  TZSPLIT => TCRRT_COMDATA%TSPLIT_B
!
!
!*       6.3   Pointer from TCRRT_COMDATA to TCRRT_PROCONF
!        for x-slices splitting
  TCRRT_COMDATA%TSPLIT_X => TCRRT_PROCONF%TSPLITS_X(IP)
!
  TZSPLIT => TCRRT_COMDATA%TSPLIT_X
!
!
!*       6.4   Pointer from TCRRT_COMDATA to TCRRT_PROCONF
!              for y-slices splitting
!
  TCRRT_COMDATA%TSPLIT_Y => TCRRT_PROCONF%TSPLITS_Y(IP)
!
  TZSPLIT => TCRRT_COMDATA%TSPLIT_Y
!
!
!*       6.5   Construction of HALO1 communication data
!
  CALL CONSTRUCT_HALO1(TCRRT_COMDATA, TCRRT_PROCONF)
  CALL CONSTRUCT_HALO2(TCRRT_COMDATA, TCRRT_PROCONF)
!
!
!*       6.6   Construction of 1D communication data
!
  ALLOCATE(TCRRT_COMDATA%HALO1DX)
  ALLOCATE(TCRRT_COMDATA%HALO1DX%NSEND_WEST(NPROC))
  ALLOCATE(TCRRT_COMDATA%HALO1DX%NSEND_EAST(NPROC))
  CALL CONSTRUCT_1DX(TCRRT_COMDATA, TCRRT_PROCONF)
!
  ALLOCATE(TCRRT_COMDATA%HALO1DY)
  ALLOCATE(TCRRT_COMDATA%HALO1DY%NSEND_SOUTH(NPROC))
  ALLOCATE(TCRRT_COMDATA%HALO1DY%NSEND_NORTH(NPROC))
  CALL CONSTRUCT_1DY(TCRRT_COMDATA, TCRRT_PROCONF)
!
!
!*       6.7   Construction of Transposition communication data
!
  CALL CONSTRUCT_TRANS(TCRRT_COMDATA, TCRRT_PROCONF)
!
!
!-------------------------------------------------------------------------------
!
!        7.    GRID NESTING :
!              ------------
!
  NULLIFY(TCRRT_PROCONF%TCHILDREN)
  NULLIFY(TCRRT_COMDATA%TCHILDREN)
  NULLIFY(TCRRT_COMDATA%TP2C_DATA)
!
  DO JMODEL = 1, JPMODELMAX
!
    IF( NDAD(JMODEL) .EQ. TCRRT_PROCONF%NUMBER ) THEN
      CALL INI_CHILD(TCRRT_PROCONF, TCRRT_COMDATA, JMODEL)
    ENDIF
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
  TZPROCONF => TCRRT_PROCONF
!
  CALL COMPUTE_TRANS_MAX(NBUFFERSIZE_3D, TCRRT_COMDATA)
  NCOMBUFFSIZE1 = NBUFFERSIZE_3D
  !JUAN NCOMBUFFSIZE1 = 10000000

!
  CALL COMPUTE_HALO_MAX(NMAXSIZEHALO, TCRRT_COMDATA)
!
!NAG4.0 boom avec le 50 lorsqu'on active les scalaires 
!  NBUFFERSIZE_2D = 50*NMAXSIZEHALO
  NBUFFERSIZE_2D = 150*NMAXSIZEHALO
!NAG4.0
  NCOMBUFFSIZE2 = NBUFFERSIZE_2D
!
  DEALLOCATE(TZDZP)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_PARA_ll
!
!     ##################################
      SUBROUTINE END_PARA_ll( KINFO_ll )
!     ##################################
!
!!****  *END_PARA_ll* - routine to finalize the parallel session
!
!!    Purpose
!!    -------
!     the purpose of the routine is to end the parallel session
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_DIM_ll
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
!     Original 01/06/98
!     R. Guivarch 15/09/99  deallocation of grid-nesting arrays
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll
!  USE MODD_STRUCTURE_ll
!  USE MODD_VAR_ll, ONLY : NIOUNIT, YOUTPUTFILE
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(OUT) :: KINFO_ll

!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    CALL TO MPI_FINALIZE
!
!  CALL CLOSE_ll(YOUTPUTFILE)

  CALL MPI_FINALIZE(KINFO_ll)
!
!-------------------------------------------------------------------------------
!
!*       2.    DEALLOCATION
!
  DEALLOCATE(NDXRATIO_ALL, NDYRATIO_ALL)
  DEALLOCATE(NXOR_ALL, NYOR_ALL) 
  DEALLOCATE(NXEND_ALL, NYEND_ALL) 
  DEALLOCATE(NDAD) 
  DEALLOCATE(CLBCX, CLBCY) 
!
END SUBROUTINE END_PARA_ll
!
END MODULE MODE_INIT_ll
