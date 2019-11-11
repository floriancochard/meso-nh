!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!     ########################
      MODULE MODE_EXCHANGE_ll
!     ########################
!
!!    Purpose
!!    -------
!
!     The purpose of this module is the implementation of communication routines
!         UPDATE_HALO_ll, UPDATE_1DHALO_ll,
!         REMAP_2WAY_X, REMAP_X_2WAY, REMAP_X_Y, REMAP_Y_X
! 
!!    Routines Of The User Interface
!!    ------------------------------
! 
!     SUBROUTINES : UPDATE_HALO_ll, UPDATE_1DHALO_ll, REMAP_2WAY_X_ll,
!                   REMAP_X_2WAY_ll, REMAP_X_Y_ll, REMAP_Y_X_ll
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!        type LIST_ll, LIST1D_ll
!
!     Module MODD_STRUCTURE_ll
!        type CRSPD_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!       IP - Number of local processor=subdomain
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       MPI_PRECISION - mpi precision
!       JPHALO - size of the halo
!       NCOMBUFFSIZE1 - buffer sizs
!       NHALO_COM - mpi communicator
!       NTRANS_COM - mpi communicator
!       NNEXTTAG, NMAXTAG - variable to define message tag
!
!     Module MODD_DIM_ll
!       CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!       CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
!       NKMAX_ll - maximum vertical dimension
!
!     Module MODD_PARAMETERS_ll
!       JPVEXT - vertical halo size
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
!       Original     May 19, 1998
!       R. Guivarch June 29, 1998 MPI_PRECISION
!       N. Gicquel, P. Kloos - October 01, 1998 - COPY_CRSPD, 
!                 COPY_ZONE, COPY_CRSPD_TRANS, COPY_ZONE_TRANS
! 
!-------------------------------------------------------------------------------
!
!     Include parameters used by MPI (already in MODE_IO_ll)
!
   USE MODD_MPIF
!  INCLUDE 'mpif.h'
!
!* 
!  
  CONTAINS
!
!     ########################################
      SUBROUTINE UPDATE_HALO_ll(TPLIST, KINFO)
!     ########################################
!
!!****  *UPDATE_HALO_ll* - routine to update halo
! 
!!    Purpose
!!    -------
!       This routine updates the halo with the values computed by the 
!     neighbor subdomains. The fields to be updated are in the
!     TPLIST list of fields. Before UPDATE_HALO_ll is called, TPLIST
!     has been filled with the fields to be communicated
!
!!**  Method
!!    ------
!       We treat first the zones the processor sends or received
!     from the others processors and then the zones it sents or
!     received from itself.
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       SEND_RECV_CRSPD, COPY_CRSPD
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_ARGSLIST_ll
!        type LIST_ll
!
!     Module MODD_VAR_ll
!       NHALO_COM - mpi communicator
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch               * CERFACS - ENSEEIHT *
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
!     Juan     19/08/2005: distinction Halo NORD/SUD & EST/WEST
!                        + modification INTENT -> INTENT(INOUT)
!     J.Escobar 13/11/2008: correction on size of buffer(IBUFFSIZE) in MPI_RECV
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_VAR_ll, ONLY : NHALO_COM, TCRRT_COMDATA
!  USE MODE_ARGSLIST_ll
!  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll
!
!*       0.1   declarations of arguments
!
  TYPE(LIST_ll), POINTER :: TPLIST ! pointer to the list of fields to be updated
  INTEGER                :: KINFO  ! return status
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
!               -------------------------------------------------------------
!
  CALL SEND_RECV_CRSPD(TCRRT_COMDATA%TSEND_HALO1, TCRRT_COMDATA%TRECV_HALO1, &
                       TPLIST, TPLIST, NHALO_COM, KINFO)
!
!*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
!               ------------------------------------------------------------
!
  CALL COPY_CRSPD(TCRRT_COMDATA%TSEND_HALO1, TCRRT_COMDATA%TRECV_HALO1, &
                  TPLIST, TPLIST, KINFO)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_HALO_ll
!
!     ##########################################
      SUBROUTINE UPDATE_1DHALO_ll(TPLIST, KINFO)
!     ##########################################
!
!
!!**** *UPDATE_1DHALO_ll* routine to update the halo of a list of one
!!                        dimensional fields.
! 
!!    Purpose
!!    -------
! 
!!    Method
!!    ------
!       For each field in the list TPLIST the UPDATE_1DHALOX_ll or the
!     UPDATE_1DHALOY_ll is called, depending on the orientation of the field
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       UPDATE_1DHALOX_ll, UPDATE_1DHALOY_ll
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST1D_ll
!
!!    Reference
!!    ---------
!     User Interface for the Mesonh Parallel Package
!
!!    Author
!!    ------
!     P. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original 12 October 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_ARGSLIST_ll, ONLY : LIST1D_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LIST1D_ll), POINTER :: TPLIST
  INTEGER, INTENT(INOUT) :: KINFO
!
!*       0.2   declarations of local variables
!
  TYPE(LIST1D_ll), POINTER :: TZLIST
!
!
!-------------------------------------------------------------------------------
!
!*       1.   GO THROUGH THE TPLIST LIST AND CALl THE APPROPRIATE ROUTINES
!             ------------------------------------------------------------
! 
  TZLIST => TPLIST
  DO WHILE(ASSOCIATED(TZLIST))
    IF (TZLIST%CDIR == "X") THEN
      CALL UPDATE_1DHALOX_ll(TZLIST%ARRAY1D, KINFO)
    ELSEIF (TZLIST%CDIR == "Y") THEN
      CALL UPDATE_1DHALOY_ll(TZLIST%ARRAY1D, KINFO)
    ELSE
      KINFO = -1
      WRITE(*,*) 'Error UPDATE_1DHALO_ll : bad CDIR value in list'
    ENDIF
    TZLIST => TZLIST%NEXT
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_1DHALO_ll
!
!     ###########################################
      SUBROUTINE UPDATE_1DHALOX_ll(PFIELD, KINFO)
!     ###########################################
!
!!**** *UPDATE_1DHALO_ll* routine to update the halo of a one dimensional
!                         field.
! 
!!    Purpose
!!    -------
! 
!!    Method
!!    ------
!       The data calculated by the CONSTRUCT_1DX routine (called by INI_PARA_ll)
!     are used to know from which proc the data will be sent and received
!
!!    Reference
!!    ---------
!!    User Interface for the Mesonh Parallel Package
!
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       GET_INDICE_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!       MPI_PRECISION - mpi precision
!       JPHALO - size of the halo
!
!     Module MODD_DIM_ll
!       CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!       CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
!
!!    Author
!!    ------
!     P. Kloos
!
!!    Modifications
!!    -------------
!     Original 12 October 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : CLBCX, CLBCY
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, MPI_PRECISION, JPHALO
!
  USE MODE_TOOLS_ll, ONLY : GET_INDICE_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!  INCLUDE 'mpif.h'
!
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:), POINTER :: PFIELD
  INTEGER, INTENT(OUT) :: KINFO
!
!*       0.2   declarations of local variables
!
  INTEGER :: IXOR, IYOR, IXEND, IYEND
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
  INTEGER :: J
  INTEGER :: ICURMODEL

! JUAN
#if defined (MNH_MPI_ISEND)
INTEGER,PARAMETER                                     :: MPI_MAX_REQ = 1024
INTEGER,SAVE,DIMENSION(MPI_MAX_REQ)                   :: REQ_TAB
INTEGER,SAVE,DIMENSION(MPI_STATUS_SIZE,MPI_MAX_REQ)   :: STATUS_TAB
INTEGER                                               :: NB_REQ
#endif
! JUAN
!
!-------------------------------------------------------------------------------
!
  ICURMODEL = TCRRT_COMDATA%NUMBER
!
!*       1.   GET DIMENSIONS
!             --------------
!
  CALL GET_INDICE_ll(IXOR, IYOR, IXEND, IYEND)
!
!-------------------------------------------------------------------------------
!
!*       2.   SEND THE HALO
!             -------------
!
!JUAN
#if defined(MNH_MPI_ISEND)
  NB_REQ = 0
#endif
!JUAN
  DO J=1, TCRRT_COMDATA%HALO1DX%NBSEND_WEST
!
#if defined(MNH_MPI_BSEND)
    CALL MPI_BSEND(PFIELD(IXOR), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DX%NSEND_WEST(J)-1, &
                  1, MPI_COMM_WORLD, KINFO)
#else
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_ISEND(PFIELD(IXOR), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DX%NSEND_WEST(J)-1, &
                  1, MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO)
#else
    CALL MPI_SEND(PFIELD(IXOR), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DX%NSEND_WEST(J)-1, &
                  1, MPI_COMM_WORLD, KINFO)
#endif
!JUAN
#endif
!
  ENDDO
!
  DO J=1, TCRRT_COMDATA%HALO1DX%NBSEND_EAST
!
#if defined(MNH_MPI_BSEND)
    CALL MPI_BSEND(PFIELD(IXEND-JPHALO+1), JPHALO, MPI_PRECISION, &
                 TCRRT_COMDATA%HALO1DX%NSEND_EAST(J)-1, &
                 2, MPI_COMM_WORLD, KINFO) 
#else
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_ISEND(PFIELD(IXEND-JPHALO+1), JPHALO, MPI_PRECISION, &
                 TCRRT_COMDATA%HALO1DX%NSEND_EAST(J)-1, &
                 2, MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO) 
#else
    CALL MPI_SEND(PFIELD(IXEND-JPHALO+1), JPHALO, MPI_PRECISION, &
                 TCRRT_COMDATA%HALO1DX%NSEND_EAST(J)-1, &
                 2, MPI_COMM_WORLD, KINFO) 
#endif
!JUAN
#endif

!
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       3.   RECEIVE THE HALO
!             ----------------
!
  IF (.NOT. LEAST_ll() .OR. CLBCX(ICURMODEL, 2) == "CYCL") THEN
!
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_IRECV(PFIELD(IXEND+1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DX%NRECV_EAST-1, 1, &
                  MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO)
#else
    CALL MPI_RECV(PFIELD(IXEND+1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DX%NRECV_EAST-1, 1, &
                  MPI_COMM_WORLD, ISTATUS, KINFO)
#endif
!JUAN
!
  ENDIF
!
  IF (.NOT. LWEST_ll() .OR. CLBCX(ICURMODEL, 1) == "CYCL") THEN
!
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_IRECV(PFIELD(1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DX%NRECV_WEST-1, 2, &
                  MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO)
#else
    CALL MPI_RECV(PFIELD(1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DX%NRECV_WEST-1, 2, &
                  MPI_COMM_WORLD, ISTATUS, KINFO)
#endif
!JUAN
!
  ENDIF

!JUAN
#if defined(MNH_MPI_ISEND)
  CALL MPI_WAITALL(NB_REQ,REQ_TAB,STATUS_TAB,KINFO) 
#endif
!JUAN

!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_1DHALOX_ll
!
!     ###########################################
      SUBROUTINE UPDATE_1DHALOY_ll(PFIELD, KINFO)
!     ###########################################
!
!!**** *UPDATE_1DHALOY_ll* routine to update the halo of a one dimensional
!                         field.
!
!!    Purpose
!!    -------
!
!!    Method
!!    ------
!       The data calculated by the CONSTRUCT_1DX routine (called by INI_PARA_ll)
!     are used to know from which proc the data will be sent and received
!!
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       GET_INDICE_ll, LNORTH_ll, LSOUTH_ll
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!       MPI_PRECISION - mpi precision
!       JPHALO - size of the halo
!
!     Module MODD_DIM_ll
!       CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!       CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : CLBCX, CLBCY
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, JPHALO, MPI_PRECISION
!
  USE MODE_TOOLS_ll, ONLY : GET_INDICE_ll, LNORTH_ll, LSOUTH_ll
!
  IMPLICIT NONE
!
!  INCLUDE 'mpif.h'
!
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:), POINTER :: PFIELD
  INTEGER, INTENT(INOUT) :: KINFO
!
!*       0.2   declarations of local variables
!
  INTEGER :: IXOR, IYOR, IXEND, IYEND
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
  INTEGER :: J
!
  INTEGER :: ICURMODEL
!
! JUAN
#if defined (MNH_MPI_ISEND)
INTEGER,PARAMETER                                     :: MPI_MAX_REQ = 1024
INTEGER,SAVE,DIMENSION(MPI_MAX_REQ)                   :: REQ_TAB
INTEGER,SAVE,DIMENSION(MPI_STATUS_SIZE,MPI_MAX_REQ)   :: STATUS_TAB
INTEGER                                               :: NB_REQ
#endif
! JUAN
!-------------------------------------------------------------------------------

  ICURMODEL = TCRRT_COMDATA%NUMBER
!
!*       1.   GET DIMENSIONS
!             --------------
!
  CALL GET_INDICE_ll(IXOR, IYOR, IXEND, IYEND)
!
!-------------------------------------------------------------------------------
!
!*       2.   SEND THE HALO
!             -------------
!
!JUAN
#if defined(MNH_MPI_ISEND)
  NB_REQ = 0
#endif
  DO J=1, TCRRT_COMDATA%HALO1DY%NBSEND_SOUTH
!
#if defined(MNH_MPI_BSEND)
    CALL MPI_BSEND(PFIELD(IYOR), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NSEND_SOUTH(J)-1, &
                  1, MPI_COMM_WORLD, KINFO)
#else
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_ISEND(PFIELD(IYOR), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NSEND_SOUTH(J)-1, &
                  1, MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO)
#else
    CALL MPI_SEND(PFIELD(IYOR), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NSEND_SOUTH(J)-1, &
                  1, MPI_COMM_WORLD, KINFO)
#endif
!JUAN
#endif
!
  ENDDO
!
  DO J=1, TCRRT_COMDATA%HALO1DY%NBSEND_NORTH
!
#if defined(MNH_MPI_BSEND)
    CALL MPI_BSEND(PFIELD(IYEND-JPHALO+1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NSEND_NORTH(J)-1, &
                  2, MPI_COMM_WORLD, KINFO) 
#else
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_ISEND(PFIELD(IYEND-JPHALO+1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NSEND_NORTH(J)-1, &
                  2, MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO)
#else
    CALL MPI_SEND(PFIELD(IYEND-JPHALO+1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NSEND_NORTH(J)-1, &
                  2, MPI_COMM_WORLD, KINFO)
#endif
!JUAN 

#endif
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       3.   RECEIVE THE HALO
!             ----------------
!
  IF (.NOT. LNORTH_ll() .OR. CLBCY(ICURMODEL, 2) == "CYCL") THEN
!
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_IRECV(PFIELD(IYEND+1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NRECV_NORTH-1, 1, &
                  MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO)
#else
    CALL MPI_RECV(PFIELD(IYEND+1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NRECV_NORTH-1, 1, &
                  MPI_COMM_WORLD, ISTATUS, KINFO)
#endif
!JUAN
!
  ENDIF
!
  IF (.NOT. LSOUTH_ll() .OR. CLBCY(ICURMODEL, 1) == "CYCL") THEN
!
!JUAN
#if defined(MNH_MPI_ISEND)
    NB_REQ = NB_REQ + 1
    CALL MPI_IRECV(PFIELD(1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NRECV_SOUTH-1, 2, &
                  MPI_COMM_WORLD, REQ_TAB(NB_REQ), KINFO)
#else
    CALL MPI_RECV(PFIELD(1), JPHALO, MPI_PRECISION, &
                  TCRRT_COMDATA%HALO1DY%NRECV_SOUTH-1, 2, &
                  MPI_COMM_WORLD, ISTATUS, KINFO)
#endif
!JUAN
!
  ENDIF
!
!JUAN
#if defined(MNH_MPI_ISEND)
  CALL MPI_WAITALL(NB_REQ,REQ_TAB,STATUS_TAB,KINFO) 
#endif
!JUAN
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_1DHALOY_ll
!
!     ######################################################
      SUBROUTINE REMAP_2WAY_X_ll(PFIELDIN, PFIELDOUT, KINFO)
!     ######################################################
!
!!**** *REMAP_2WAY_X_ll* -
!!
!!    Purpose
!!    -------
!       This routine remaps a field PFIELDIN according to a 2way splitting
!       to a field PFIELDOUT according to a x-slices splitting
!
!!**  Method
!!    ------
!!      First the processors send and receive the subsets of their local domain
!!    according to 2way and x-slices splittings from neighboring
!!    processors/subdomains different from themselves,
!!    and then they treat the zones they send or receive from themselves.
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       SEND_RECV_FIELD, COPY_CRSPD_TRANS
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch               * CERFACS - ENSEEIHT *
!!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
!-------------------------------------------------------------------------------
!
!*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
!               -------------------------------------------------------------
!
  CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_TRANS_BX, &
                       TCRRT_COMDATA%TRECV_TRANS_BX, &
                       PFIELDIN, PFIELDOUT, KINFO)
!
!-------------------------------------------------------------------------------
!
!*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
!               ------------------------------------------------------------
!
  CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_TRANS_BX, &
                        TCRRT_COMDATA%TRECV_TRANS_BX, &
                        PFIELDIN, PFIELDOUT, KINFO)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REMAP_2WAY_X_ll
!
!     ######################################################
      SUBROUTINE REMAP_X_2WAY_ll(PFIELDIN, PFIELDOUT, KINFO)
!     ######################################################
!
!!****  *REMAP_X_2WAY_ll* -
!!
!!    Purpose
!!    -------
!       This routine remaps a field PFIELDIN according to a x-slices splitting
!     to a field PFIELDOUT according to a 2way splitting
!
!!**  Method
!!    ------
!       First the processors send and receive the subsets of their local domain
!     according to x-slices and y-slices splittings from neighboring
!     processors/subdomains different from themselves,
!     and then they treat the zones they send or receive from themselves.
!
!
!!    External
!!    --------
!     Module MODE_ARGSLIST_ll
!       ADD3DFIELD_ll, CLEANLIST_ll
!
!     Module MODE_EXCHANGE_ll
!       SEND_RECV_FIELD, COPY_CRSPD_TRANS
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch               * CERFACS - ENSEEIHT *
!!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  USE MODE_ARGSLIST_ll, ONLY : ADD3DFIELD_ll, CLEANLIST_ll
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
!*       0.2   declarations of local variables
!
  TYPE(LIST_ll), POINTER :: TZLIST
!
!-------------------------------------------------------------------------------
!
!*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
!               -------------------------------------------------------------
!
  CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_TRANS_BX, &
                       TCRRT_COMDATA%TSEND_TRANS_BX, &
                       PFIELDIN, PFIELDOUT, KINFO)
!
!-------------------------------------------------------------------------------
!
!*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
!               ------------------------------------------------------------
!
  CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_TRANS_BX, &
                        TCRRT_COMDATA%TSEND_TRANS_BX, &
                        PFIELDIN, PFIELDOUT, KINFO)
!
!*       3.    UPDATE HALO :
!              -----------
!
  NULLIFY(TZLIST)
  CALL ADD3DFIELD_ll(TZLIST,PFIELDOUT)
!
  CALL UPDATE_HALO_ll(TZLIST, KINFO)
  CALL CLEANLIST_ll(TZLIST)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REMAP_X_2WAY_ll
!
!     ###################################################
      SUBROUTINE REMAP_X_Y_ll(PFIELDIN, PFIELDOUT, KINFO)
!     ###################################################
!
!!****  *REMAP_X_Y_ll* -
! 
!!    Purpose
!     -------
!       This routine remaps a field PFIELDIN according to a x-slices splitting
!     to a field PFIELDOUT according to a y-slices splitting
!
!!**  Method
!!    ------
!       First the processors send and receive the subsets of their local domain
!     according to y-slices and x-slices splittings from neighboring
!     processors/subdomains different from themselves,
!     and then they treat the zones they send or receive from themselves.
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       SEND_RECV_FIELD, COPY_CRSPD_TRANS
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch               * CERFACS - ENSEEIHT *
!!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
!-------------------------------------------------------------------------------
!
!*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
!               -------------------------------------------------------------
!
  CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_TRANS_XY, &
                       TCRRT_COMDATA%TRECV_TRANS_XY, &
                       PFIELDIN, PFIELDOUT, KINFO)
!
!-------------------------------------------------------------------------------
!
!*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
!               ------------------------------------------------------------
!
  CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_TRANS_XY, &
                        TCRRT_COMDATA%TRECV_TRANS_XY, &
                        PFIELDIN, PFIELDOUT, KINFO)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REMAP_X_Y_ll
!
!     ###################################################
      SUBROUTINE REMAP_Y_X_ll(PFIELDIN, PFIELDOUT, KINFO)
!     ###################################################
!
!!****  *REMAP_Y_X_ll* -
!
!!    Purpose
!!    -------
!       This routine remaps a field PFIELDIN according to a y-slices splitting
!     to a field PFIELDOUT according to a x-slices splitting
!
!!**  Method
!!    ------
!       First the processors send and receive the subsets of their local domain
!     according to x-slices and 2-way splittings from neighboring
!     processors/subdomains different from themselves,
!     and then they treat the zones they send or receive from themselves.
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       SEND_RECV_FIELD, COPY_CRSPD_TRANS
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch               * CERFACS - ENSEEIHT *
!!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
  INTEGER :: KINFO ! return status
!
!-------------------------------------------------------------------------------
!
!*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
!               -------------------------------------------------------------
!
  CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_TRANS_XY, &
                       TCRRT_COMDATA%TSEND_TRANS_XY, &
                       PFIELDIN, PFIELDOUT, KINFO)
!
!-------------------------------------------------------------------------------
!
!*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
!               ------------------------------------------------------------
!
  CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_TRANS_XY, &
                        TCRRT_COMDATA%TSEND_TRANS_XY, &
                        PFIELDIN, PFIELDOUT, KINFO)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE REMAP_Y_X_ll
!
!     ###################################################
      SUBROUTINE COPY_CRSPD(TPSEND_CRSPD, TPRECV_CRSPD, &
                            TPSENDLIST, TPRECVLIST, KINFO)
!     ####################################################
!
!!****  *COPY_CRSPD* -
!
!!    Purpose
!!    -------
!     copy the zones the a process sends to itself
!     (instead of sending them via MPI)
!
!!**  Method
!!    ------
!     we go over all the zones of the TPSEND_CRSPD variable and find
!     out those whose NUMBER equals IP.
!     In this case the recipient of a ZONE is the same as the sender.
!     To find out where the zone has to be copied we go over the
!     TPRECV_CRSPD variable.
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       COPY_ZONE
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!        type CRSPD_ll
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll
  USE MODD_VAR_ll, ONLY : IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!

  TYPE(CRSPD_ll), POINTER :: TPSEND_CRSPD  ! CRSPD to be sent
  TYPE(CRSPD_ll), POINTER :: TPRECV_CRSPD  ! CRSPD to be received
  TYPE(LIST_ll), POINTER  :: TPSENDLIST, & ! pointer to the list of fields
                             TPRECVLIST    ! to be updated
  INTEGER                 :: KINFO         ! return status
!
!*       0.2   declarations of local variables
!
  TYPE(CRSPD_ll), POINTER :: TZSEND_CRSPD, TZRECV_CRSPD
  INTEGER :: JSEND
  INTEGER :: IMSSGTAG
!
!-------------------------------------------------------------------------------
!
!*       1.    GO OVER THE TPSEND_CRSPD LIST OF ZONES
!              --------------------------------------
!
  IF (.NOT.ASSOCIATED(TPSEND_CRSPD)) THEN
    RETURN
  ENDIF
!
  TZSEND_CRSPD => TPSEND_CRSPD
  DO JSEND = 1, TPSEND_CRSPD%NCARD
!
!*       1.1   Test whether a zone is sent to the same proc
!
    IF (TZSEND_CRSPD%TELT%NUMBER == IP) THEN
!
!*       1.2   If so, go over the TPRECV_CRSPD list of zones
!*             and test whether the zone to be received corresponds to
!*             the zone to be sent
!
      TZRECV_CRSPD => TPRECV_CRSPD
      DO WHILE (ASSOCIATED(TZRECV_CRSPD))
        IF (TZRECV_CRSPD%TELT%NUMBER == IP &
           .AND. TZRECV_CRSPD%TELT%MSSGTAG == TZSEND_CRSPD%TELT%MSSGTAG) THEN
!
!*       1.2.1 If so, copy the zone
!
            CALL COPY_ZONE(TZSEND_CRSPD%TELT, TZRECV_CRSPD%TELT, &
                           TPSENDLIST, TPRECVLIST, KINFO)
!
        ENDIF
!
        TZRECV_CRSPD => TZRECV_CRSPD%TNEXT
!
      ENDDO
!
    ENDIF
!
    TZSEND_CRSPD => TZSEND_CRSPD%TNEXT
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COPY_CRSPD
!
!     ###################################################################
      SUBROUTINE COPY_ZONE(TPSEND, TPRECV, TPSENDLIST, TPRECVLIST, KINFO)
!     ###################################################################
!
!!****  *COPY_ZONE* -
!
!!    Purpose
!!    -------
!       This routine copies the values of the fields in the TPLIST
!     list of fields, from the zone discribed by TPSEND into
!     the TPRECV zone.
!!
!!**  Method
!!    ------
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!        type ZONE_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!  USE MODE_ARGSLIST_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!

  TYPE(ZONE_ll) :: TPSEND                 ! ZONE_ll to be sent
  TYPE(ZONE_ll) :: TPRECV                 ! ZONE_ll to be received
  TYPE(LIST_ll), POINTER :: TPSENDLIST, & ! pointer to the list 
                            TPRECVLIST    ! of fields to be updated
!
  INTEGER                :: KINFO       ! return status
!
!*       0.2   declarations of local variables
!
  TYPE(LIST_ll), POINTER :: TZSENDLIST, TZRECVLIST
  INTEGER :: IIBS, IIES, IJBS, IJES, IIBR, IIER, IJBR, IJER
!
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISE DIMENSIONS
!              ---------------------
!
  IIBS = TPSEND%NXOR
  IIES = TPSEND%NXEND
  IJBS = TPSEND%NYOR
  IJES = TPSEND%NYEND
!
  IIBR = TPRECV%NXOR
  IIER = TPRECV%NXEND
  IJBR = TPRECV%NYOR
  IJER = TPRECV%NYEND
!
!-------------------------------------------------------------------------------
!
!*       2.    COPY THE TPSEND ZONE INTO THE TPRECV ZONE
!              -----------------------------------------
!
!*       2.1   Go over the TPLIST list of fields
!
  TZSENDLIST => TPSENDLIST
  TZRECVLIST => TPRECVLIST
  DO WHILE (ASSOCIATED(TZSENDLIST) .AND. ASSOCIATED(TZRECVLIST))
    IF (TZSENDLIST%L2D) THEN
      TZRECVLIST%ARRAY2D(IIBR:IIER,IJBR:IJER) = &
                                        TZSENDLIST%ARRAY2D(IIBS:IIES,IJBS:IJES)
    ELSEIF(TZSENDLIST%L3D) THEN
      TZRECVLIST%ARRAY3D(IIBR:IIER,IJBR:IJER,:) = &
                                      TZSENDLIST%ARRAY3D(IIBS:IIES,IJBS:IJES,:)
    ENDIF
    TZSENDLIST => TZSENDLIST%NEXT
    TZRECVLIST => TZRECVLIST%NEXT
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COPY_ZONE
!
!     #######################################################################
      SUBROUTINE COPY_CRSPD_TRANS(TPSEND, TPRECV, PFIELDIN, PFIELDOUT, KINFO)
!     #######################################################################
!
!!****  *COPY_CRSPD_TRANS* -
!
!!    Purpose
!!    -------
!       copy the zones a process sends to itself, instead of sending them
!!    via MPI.
!!
!!**  Method
!!    ------
!     we go over all the zones of the TPSEND variable and find
!     out those whose NUMBER equals IP.
!     In this case the recipient of a ZONE is the same as the sender.
!     To find out where the zone has to be copied we go over the
!     TPRECV variable.
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       COPY_ZONE_TRANS
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!        type CRSPD_ll
!
!     Module MODD_VAR_ll
!       IP - Number of local processor=subdomain
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll
  USE MODD_VAR_ll, ONLY : IP
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPSEND     ! CRSPD to be sent
  TYPE(CRSPD_ll), POINTER :: TPRECV     ! CRSPD to be received
  REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PFIELDIN
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT
  INTEGER                 :: KINFO       ! return status
!
!*       0.2   declarations of local variables
!
  TYPE(CRSPD_ll), POINTER :: TZSEND, TZRECV
  INTEGER :: JSEND
!
!-------------------------------------------------------------------------------
!
!*       1.    GO OVER THE TPSEND_CRSPD LIST OF ZONES
!              --------------------------------------
!
  IF (.NOT.ASSOCIATED(TPSEND)) THEN
    RETURN
  ENDIF
!
  TZSEND => TPSEND
  DO JSEND = 1, TPSEND%NCARD
!
!*       1.1   Test whether a zone is sent to the same proc
!
    IF (TZSEND%TELT%NUMBER == IP) THEN
!
!*       1.2   If so, go over the TPRECV list of zones
!*             and test whether the zone to be received corresponds to
!*             the zone to be sent
!
      TZRECV => TPRECV
      DO WHILE (ASSOCIATED(TZRECV))
        IF (TZRECV%TELT%NUMBER == IP &
           .AND. TZRECV%TELT%MSSGTAG == TZSEND%TELT%MSSGTAG) THEN
! 
!*       1.2.1 If so, copy the zone
!
          CALL COPY_ZONE_TRANS(TZSEND%TELT, TZRECV%TELT, PFIELDIN, PFIELDOUT, &
                               KINFO)
!
        ENDIF
!
        TZRECV => TZRECV%TNEXT
!
      ENDDO
!
    ENDIF
!
    TZSEND => TZSEND%TNEXT
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COPY_CRSPD_TRANS
!
!     ######################################################################
      SUBROUTINE COPY_ZONE_TRANS(TPSEND, TPRECV, PFIELDIN, PFIELDOUT, KINFO)
!     ######################################################################
!
!!****  *COPY_ZONE_TRANS* -
!
!!    Purpose
!!    -------
!       this routine copies the values of the PFIELDIN field situated
!     in the TPSEND zone, into the PFIELDOUT field at the TPRECV zone
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!       type ZONE_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!

  TYPE(ZONE_ll) :: TPSEND     ! ZONE_ll to be sent
  TYPE(ZONE_ll) :: TPRECV     ! ZONE_ll to be received
  REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PFIELDIN
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT
  INTEGER                :: KINFO      ! return status
!
!*       0.2   declarations of local variables
!
  TYPE(LIST_ll), POINTER :: TZLIST
  INTEGER :: IIBS, IIES, IJBS, IJES, IKBS, IKES, IIBR, IIER, &
             IJBR, IJER, IKBR, IKER
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISE DIMENSIONS
!              ---------------------
!
  IIBS = TPSEND%NXOR
  IIES = TPSEND%NXEND
  IJBS = TPSEND%NYOR
  IJES = TPSEND%NYEND
  IKBS = TPSEND%NZOR
  IKES = TPSEND%NZEND
!
  IIBR = TPRECV%NXOR
  IIER = TPRECV%NXEND
  IJBR = TPRECV%NYOR
  IJER = TPRECV%NYEND
  IKBR = TPRECV%NZOR
  IKER = TPRECV%NZEND
!
!*       2.    COPY THE VALUES OF PFIELDIN SITUATED IN THE TPSEND ZONE
!*             IN THE ENTRIES OF PFIELDOUT DEFINED BY TPRECV
!              -------------------------------------------------------
!
  PFIELDOUT(IIBR:IIER,IJBR:IJER,IKBR:IKER) = &
                                        PFIELDIN(IIBS:IIES,IJBS:IJES,IKBS:IKES)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COPY_ZONE_TRANS
!
!     ########################################################
      SUBROUTINE FILLIN_BUFFER(TPFIELD, TPZONE, PBUFFER, KINC)
!     ########################################################
!
!!****  *FILLIN_BUFFER* -
!
!!    Purpose
!!    -------
!       this routine fills the PBUFFER buffer with the values of TPFIELD field
!     situated in the TPZONE zone
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!       type ZONE_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LIST_ll), POINTER :: TPFIELD
  TYPE(ZONE_ll) :: TPZONE
  REAL, DIMENSION(:), INTENT(INOUT) :: PBUFFER
  INTEGER :: KINC
!
!*       0.2   declarations of local variables
!
  INTEGER :: JI,JJ,JK
!
! Local index for single vector loop use
  INTEGER :: JKMAX,JJMIN,JJMAX,JIMIN,JIMAX
  INTEGER :: JKOR,JKEND,JJOR,JJEND,JIOR,JIEND

  INTEGER :: JIJ,JIJMAX
  INTEGER :: JIBOX,JJBOX
!
!-------------------------------------------------------------------------------
!
  IF (.NOT.ASSOCIATED(TPFIELD)) THEN
      RETURN
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       1.    FILL IN THE BUFFER
!              ------------------
!
  IF (TPFIELD%L2D) THEN
!
      DO JJ=TPZONE%NYOR, TPZONE%NYEND
        DO JI=TPZONE%NXOR, TPZONE%NXEND
!
          KINC = KINC + 1
          PBUFFER(KINC) = TPFIELD%ARRAY2D(JI,JJ)
!
        ENDDO
      ENDDO
!
  ELSEIF(TPFIELD%L3D) THEN
!
!        Z axe
      JKMAX = SIZE(TPFIELD%ARRAY3D,3)
      JKOR  = TPZONE%NZOR
      JKEND = TPZONE%NZEND
!        Y axe
      JJMAX = SIZE(TPFIELD%ARRAY3D,2)
      JJOR  = TPZONE%NYOR
      JJEND = TPZONE%NYEND
!        X axe
      JIMAX = SIZE(TPFIELD%ARRAY3D,1)
      JIOR  = TPZONE%NXOR
      JIEND = TPZONE%NXEND
!
      JIBOX  = JIEND - JIOR + 1
      JJBOX  = JJEND - JJOR + 1
      JIJMAX = JIBOX*JJBOX
!      JINC = 0
        DO JK=1, SIZE(TPFIELD%ARRAY3D,3)
!CDIR NODEP
!OCL NOVREC
           DO JIJ=1,JIJMAX
              JI = JIOR + MOD( JIJ -1   , JIBOX )
              JJ = JJOR +    ( JIJ -1 ) / JIBOX
              PBUFFER(JIJ+KINC) = TPFIELD%ARRAY3D(JI, JJ, JK)
           ENDDO
           KINC = KINC + JIJMAX   
        ENDDO
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE FILLIN_BUFFER
!
!     #########################################################
      SUBROUTINE FILLOUT_BUFFER(TPFIELD, TPZONE, PBUFFER, KINC)
!     #########################################################
!
!!****  *FILLOUT_BUFFER* -
!
!!    Purpose
!!    -------
!       this routine fills out the PBUFFER buffer in the field TPFIELD
!     according to the TPZONE zone
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!       type ZONE_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LIST_ll), POINTER :: TPFIELD
  TYPE(ZONE_ll) :: TPZONE
  REAL, DIMENSION(:), INTENT(INOUT) :: PBUFFER
  INTEGER :: KINC
!
!*       0.2   declarations of local variables
!
  INTEGER :: JI,JJ,JK
!
!JUAN
  INTEGER :: JKMAX,JJMIN,JJMAX,JIMIN,JIMAX
  INTEGER :: JKOR,JKEND,JJOR,JJEND,JIOR,JIEND

  INTEGER :: JIJ,JIJMAX
  INTEGER :: JIBOX,JJBOX
!JUAN
!-------------------------------------------------------------------------------
!
  IF (.NOT.ASSOCIATED(TPFIELD)) THEN
      RETURN
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       1.    FILL OUT THE BUFFER
!              -------------------
!
  IF (TPFIELD%L2D) THEN
!
      DO JJ=TPZONE%NYOR, TPZONE%NYEND
        DO JI=TPZONE%NXOR, TPZONE%NXEND
!
          KINC = KINC + 1
          TPFIELD%ARRAY2D(JI,JJ) = PBUFFER(KINC) 
!
        ENDDO
      ENDDO
!
  ELSEIF(TPFIELD%L3D) THEN
!
!    Z axe
     JKMAX = SIZE(TPFIELD%ARRAY3D,3)
     JKOR  = TPZONE%NZOR
     JKEND = TPZONE%NZEND
!    Y axe
     JJMAX = SIZE(TPFIELD%ARRAY3D,2)
     JJOR  = TPZONE%NYOR
     JJEND = TPZONE%NYEND
!    X axe
     JIMAX = SIZE(TPFIELD%ARRAY3D,1)
     JIOR  = TPZONE%NXOR
     JIEND = TPZONE%NXEND 
!
     JIBOX  = JIEND - JIOR + 1
     JJBOX  = JJEND - JJOR + 1
     JIJMAX = JIBOX*JJBOX
!
     DO JK = 1, SIZE(TPFIELD%ARRAY3D,3)
!CDIR NODEP
!OCL NOVREC
        DO JIJ=1,JIJMAX
           JI = JIOR + MOD( JIJ -1   , JIBOX )
           JJ = JJOR +    ( JIJ -1 ) / JIBOX
           TPFIELD%ARRAY3D(JI, JJ, JK) = PBUFFER(JIJ+KINC)
        ENDDO
        KINC = KINC + JIJMAX
     ENDDO
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE FILLOUT_BUFFER
!
!     ##########################################################
      SUBROUTINE FILLOUT_BUFFERS(TPFIELD, TPZONE, PBUFFER, KINC)
!     ##########################################################
!
!!****  *FILLOUT_BUFFERS* -
!
!!    Purpose
!!    -------
!       this routine fills out the PBUFFER buffer in the fields of TPFIELD
!     fields list according to the TPZONE zone
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       FILLOUT_BUFFER
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!       type ZONE_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LIST_ll), POINTER :: TPFIELD, TZFIELD
  TYPE(ZONE_ll) :: TPZONE
  REAL, DIMENSION(:) :: PBUFFER
  INTEGER :: KINC
!
!-------------------------------------------------------------------------------
!
!*       1.    GO OVER THE TPFIELD LIST OF FIELDS
!              ----------------------------------
!
  TZFIELD => TPFIELD
  DO WHILE (ASSOCIATED(TZFIELD))
     CALL FILLOUT_BUFFER(TZFIELD, TPZONE, PBUFFER, KINC)
     TZFIELD => TZFIELD%NEXT
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE FILLOUT_BUFFERS
!
!     #########################################################
      SUBROUTINE FILLIN_BUFFERS(TPFIELD, TPZONE, PBUFFER, KINC)
!     #########################################################
!
!!****  *FILLIN_BUFFERS* -
!
!!    Purpose
!!    -------
!       this routine fills the PBUFFER buffer with the values of the fields
!     of TPFIELD field list situated in the TPZONE zone
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       FILLIN_BUFFER
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_STRUCTURE_ll
!       type ZONE_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel                * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!!     1 october 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LIST_ll), POINTER :: TPFIELD, TZFIELD
  TYPE(ZONE_ll) :: TPZONE
  REAL, DIMENSION(:) :: PBUFFER
  INTEGER :: KINC
!
!-------------------------------------------------------------------------------!
!*       1.    GO OVER THE TPFIELD LIST OF FIELDS
!              ----------------------------------
!
  KINC = 0 
  TZFIELD => TPFIELD
  DO WHILE (ASSOCIATED(TZFIELD))
     CALL FILLIN_BUFFER(TZFIELD, TPZONE, PBUFFER, KINC)
     TZFIELD => TZFIELD%NEXT 
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE FILLIN_BUFFERS
!
!     ######################################################
      SUBROUTINE SEND_RECV_FIELD(TPCRSPDSEND, TPCRSPDRECV, &
                                 PFIELDIN, PFIELDOUT, KINFO)
!     ######################################################
!
!!****  *SEND_RECV_FIELD* -
!
!!    Purpose
!!    -------
!     This routine sends the data of the PFIELDIN field 
!     to the correspondants of the TPCRSPDSEND list 
!     and receives the data of the PFIELDOUT field
!     from the correspondants  of the TPCRSPDRECV list
!
!!**  Method
!!    ------
!
!     This routine is based on the following rule : 
!        one sent for one received (if it is possible);
!     The algorithm is the following :
!
!     while (there is some messages to send)
!       OR  (there is some messages to receive)
!     do
!       if (there is some messages to send)
!         send the first of the list
!       end if
!       if (there is some messages to receive)
!         try to receive
!       end if
!     done
!
!     The receptions are non-blocking and don't follow necessarly
!     the order of the TPCRSPDRECV list
!
!!    External
!!    -------
!     Module MODE_TOOLS_ll
!       GET_MAX_SIZE
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type CRSPD_ll, ZONE_ll
!
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_VAR_ll
!       IP - Number of local processor=subdomain
!       NCOMBUFFSIZE1 - buffer size
!       NTRANS_COM - mpi communicator
!       MPI_PRECISION - mpi precision
!       NNEXTTAG, NMAXTAG - variable to define message tag
!
!     Module MODD_PARAMETERS_ll
!       JPVEXT - vertical halo size
!
!     Module MODD_DIM_ll
!       NKMAX_ll - maximum vertical dimension
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel               * CERFACS - CNRM *
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll, ZONE_ll
  USE MODD_VAR_ll, ONLY : NCOMBUFFSIZE1, IP, NTRANS_COM, MPI_PRECISION, &
                          NNEXTTAG, NMAXTAG
  USE MODD_PARAMETERS_ll, ONLY : JPVEXT
  USE MODD_DIM_ll, ONLY : NKMAX_ll
!
  USE MODE_TOOLS_ll, ONLY : GET_MAX_SIZE
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPCRSPDSEND, TPCRSPDRECV
  REAL, DIMENSION(:, :, :), INTENT(INOUT) :: PFIELDIN
  REAL, DIMENSION(:, :, :), INTENT(INOUT) :: PFIELDOUT
!
  INTEGER :: KINFO
!
!*       0.2   declarations of local variables
!
  INTEGER :: JINC, JI, JJ, JK  ! Loop and counter variables
  INTEGER :: FOUND, KERROR
! 
  TYPE(CRSPD_ll), POINTER :: TZMAILSEND, TZMAILRECV
  TYPE(ZONE_ll), POINTER :: TZZONESEND, TZZONERECV
!
! JUAN
#if defined (MNH_MPI_ISEND)
  REAL, DIMENSION (:,:), ALLOCATABLE :: TZBUFFER
#else
  REAL, DIMENSION (:), ALLOCATABLE :: TZBUFFER
#endif
! JUAN
!
  INTEGER IRECVSTATUS(MPI_STATUS_SIZE) ! Status of completed receive request
! 
  LOGICAL :: IRECVFLAG
  INTEGER :: IMSGTAG, ISENDERPROC
!
  INTEGER :: IRECVNB, ISENDNB ! Total numbers of receive and send to do
  INTEGER :: IRECVDONE
!
  INTEGER, SAVE :: ITAGOFFSET = 0
!
  INTEGER :: IMAXSIZESEND, IMAXSIZERECV, IBUFFSIZE
!
!JUAN
  INTEGER :: JKMAX,JJMIN,JJMAX,JIMIN,JIMAX
  INTEGER :: JKOR,JKEND,JJOR,JJEND,JIOR,JIEND

  INTEGER :: JIJ,JIJMAX
  INTEGER :: JIBOX,JJBOX
!JUAN
! JUAN
#if defined (MNH_MPI_ISEND)
INTEGER,PARAMETER                                     :: MPI_MAX_REQ = 1024
INTEGER,SAVE,DIMENSION(MPI_MAX_REQ)                   :: REQ_TAB
INTEGER,SAVE,DIMENSION(MPI_STATUS_SIZE,MPI_MAX_REQ)   :: STATUS_TAB
INTEGER                                               :: NB_REQ,NFIRST_REQ_RECV
#endif
! JUAN
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATIONS
!              ---------------
!
  IRECVNB = 0
  ISENDNB = 0
  IRECVDONE = 0
!
!*       1.1   computation of the buffer'size
!
  IF (.NOT.ASSOCIATED(TPCRSPDSEND)) THEN
    ISENDNB = 0
  ELSE
    ISENDNB = TPCRSPDSEND%NCARDDIF
    IMAXSIZESEND = GET_MAX_SIZE(TPCRSPDSEND)
  ENDIF
!
  IF (.NOT.ASSOCIATED(TPCRSPDRECV)) THEN
    IRECVNB = 0
  ELSE
    IRECVNB = TPCRSPDRECV%NCARDDIF
    IMAXSIZERECV = GET_MAX_SIZE(TPCRSPDRECV)
  ENDIF
!
  IBUFFSIZE = IMAXSIZESEND 
  IF (IMAXSIZERECV > IBUFFSIZE) IBUFFSIZE = IMAXSIZERECV
!
  IBUFFSIZE = IBUFFSIZE * (NKMAX_ll + 2 * JPVEXT)
!
! JUAN
#if defined (MNH_MPI_ISEND)
  ALLOCATE(TZBUFFER(IBUFFSIZE,ISENDNB+IRECVNB))
  TZBUFFER(:,:) = 0.0 
#else
  ALLOCATE(TZBUFFER(IBUFFSIZE))
  TZBUFFER(:) = 0.0 
#endif
! JUAN
!
  TZMAILRECV => TPCRSPDRECV
  TZMAILSEND => TPCRSPDSEND
!
!NZJUAN  CALL MPI_BARRIER(NTRANS_COM, KERROR)
!
!-------------------------------------------------------------------------------
!
!*       2.    MAIN LOOP
!              ---------
!
! JUAN
#if defined (MNH_MPI_ISEND)
  NB_REQ = 0
#endif
! JUAN
  DO WHILE (ASSOCIATED(TZMAILSEND))
!
!*       2.1  if there is still something to send
!
     IF (ASSOCIATED(TZMAILSEND)) THEN		
        TZZONESEND => TZMAILSEND%TELT
        IF (TZZONESEND%NUMBER /= IP) THEN
! JUAN
#if defined (MNH_MPI_ISEND)
           NB_REQ = NB_REQ + 1
#endif
! JUAN 
!          Z axe
           JKMAX = SIZE(PFIELDIN,3)
           JKOR  = TZZONESEND%NZOR
           JKEND = TZZONESEND%NZEND
!          Y axe
           JJMAX = SIZE(PFIELDIN,2)
           JJOR  = TZZONESEND%NYOR
           JJEND = TZZONESEND%NYEND
!          X axe
           JIMAX = SIZE(PFIELDIN,1)
           JIOR  = TZZONESEND%NXOR
           JIEND = TZZONESEND%NXEND   

           JINC = 0
           JIBOX  = JIEND - JIOR + 1
           JJBOX  = JJEND - JJOR + 1
           JIJMAX = JIBOX*JJBOX

           DO JK=TZZONESEND%NZOR, TZZONESEND%NZEND
!CDIR NODEP
!OCL NOVREC
              DO JIJ=1,JIJMAX
                 JI = JIOR + MOD( JIJ -1   , JIBOX )
                 JJ = JJOR +    ( JIJ -1 ) / JIBOX
!JUAN
#if defined (MNH_MPI_ISEND)
                 TZBUFFER(JIJ+JINC,NB_REQ) = PFIELDIN(JI, JJ, JK)
#else
                 TZBUFFER(JIJ+JINC) = PFIELDIN(JI, JJ, JK)
#endif
!JUAN
              ENDDO
              JINC = JINC + JIJMAX   
           ENDDO
!
#if defined(MNH_MPI_BSEND)
           CALL MPI_BSEND(TZBUFFER, JINC, MPI_PRECISION, TZZONESEND%NUMBER - 1, &
                TZZONESEND%MSSGTAG + ITAGOFFSET, NTRANS_COM, KERROR)
#else
!JUAN
#if defined(MNH_MPI_ISEND)
           CALL MPI_ISEND(TZBUFFER(1,NB_REQ), JINC, MPI_PRECISION, TZZONESEND%NUMBER - 1, &
                TZZONESEND%MSSGTAG + ITAGOFFSET, NTRANS_COM, REQ_TAB(NB_REQ), KERROR)
#else
           CALL MPI_SEND(TZBUFFER, JINC, MPI_PRECISION, TZZONESEND%NUMBER - 1, &
                TZZONESEND%MSSGTAG + ITAGOFFSET, NTRANS_COM, KERROR)
#endif
!JUAN
#endif
        ENDIF
        TZMAILSEND => TZMAILSEND%TNEXT
     ENDIF
!
  ENDDO

!NZJUAN  CALL MPI_BARRIER(NTRANS_COM, KERROR)

! JUAN
#if defined (MNH_MPI_ISEND)
  NFIRST_REQ_RECV = NB_REQ
#endif
! JUAN 

  DO WHILE (ASSOCIATED(TZMAILRECV))
     TZZONERECV => TZMAILRECV%TELT 
     IF (TZZONERECV%NUMBER == IP) THEN
        TZMAILRECV =>    TZMAILRECV%TNEXT
     ELSE
! JUAN
#if defined (MNH_MPI_ISEND)
        NB_REQ = NB_REQ + 1
!JUAN NZ   CALL MPI_IRECV(TZBUFFER(1,NB_REQ), NCOMBUFFSIZE1, MPI_PRECISION, &
        CALL MPI_IRECV(TZBUFFER(1,NB_REQ), IBUFFSIZE, MPI_PRECISION, &
             TZZONERECV%NUMBER-1, TZZONERECV%MSSGTAG + ITAGOFFSET, &
             NTRANS_COM, REQ_TAB(NB_REQ), KERROR)
#else
!JUAN NZ        CALL MPI_RECV(TZBUFFER, NCOMBUFFSIZE1, MPI_PRECISION, TZZONERECV%NUMBER-1, &
        CALL MPI_RECV(TZBUFFER, IBUFFSIZE, MPI_PRECISION, TZZONERECV%NUMBER-1, &
             TZZONERECV%MSSGTAG + ITAGOFFSET, NTRANS_COM, IRECVSTATUS, KERROR)
#endif
! JUAN  
        TZMAILRECV =>    TZMAILRECV%TNEXT
!JUAN
#if defined(MNH_MPI_ISEND)
     ENDIF
  ENDDO
 

  CALL MPI_WAITALL(NB_REQ,REQ_TAB,STATUS_TAB,KINFO) 


  TZMAILRECV => TPCRSPDRECV
  NB_REQ = NFIRST_REQ_RECV

  DO WHILE (ASSOCIATED(TZMAILRECV))
     TZZONERECV => TZMAILRECV%TELT 
     IF (TZZONERECV%NUMBER == IP) THEN
        TZMAILRECV =>    TZMAILRECV%TNEXT
     ELSE
        NB_REQ = NB_REQ + 1 
        TZMAILRECV =>    TZMAILRECV%TNEXT 
#endif
!JUAN 
!       Z axe
        JKMAX = SIZE(PFIELDOUT,3)
        JKOR  = TZZONERECV%NZOR
        JKEND = TZZONERECV%NZEND
!       Y axe
        JJMAX = SIZE(PFIELDOUT,2)
        JJOR  = TZZONERECV%NYOR
        JJEND = TZZONERECV%NYEND
!       X axe
        JIMAX = SIZE(PFIELDOUT,1)
        JIOR  = TZZONERECV%NXOR
        JIEND = TZZONERECV%NXEND  
!        
        JINC = 0
        JIBOX  = JIEND - JIOR + 1
        JJBOX  = JJEND - JJOR + 1
        JIJMAX = JIBOX*JJBOX
        DO JK = TZZONERECV%NZOR, TZZONERECV%NZEND
!CDIR NODEP
!OCL NOVREC
           DO JIJ=1,JIJMAX
              JI = JIOR + MOD( JIJ -1   , JIBOX )
              JJ = JJOR +    ( JIJ -1 ) / JIBOX
! JUAN 
#if defined (MNH_MPI_ISEND)
              PFIELDOUT(JI, JJ, JK) = TZBUFFER(JIJ+JINC,NB_REQ)
#else
              PFIELDOUT(JI, JJ, JK) = TZBUFFER(JIJ+JINC)
#endif
! JUAN
           ENDDO
           JINC = JINC + JIJMAX
        ENDDO
!
     ENDIF

!
  ENDDO
!
!NZJUAN  CALL MPI_BARRIER(NTRANS_COM, KERROR)
!
  ITAGOFFSET = MOD((ITAGOFFSET + NNEXTTAG), NMAXTAG)
!
  DEALLOCATE(TZBUFFER)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SEND_RECV_FIELD
!
!     ##############################################################
      SUBROUTINE SEND_RECV_CRSPD(TPCRSPDSEND, TPCRSPDRECV, &
                                 TPFIELDLISTSEND, TPFIELDLISTRECV, &
                                 KMPI_COMM, KINFO, KBARRIER)
!     ##############################################################
!
!!****  *SEND_RECV_CRSPD*-
!
!!    Purpose
!!    -------
!       This routine sends the data of the TPFIELDLISTSEND list 
!       to the correspondants of the TPCRSPDSEND list 
!       and receives the data of the TPFIELDLISTRECV list
!       from the correspondants  of the TPCRSPDRECV list
!
!!**  Method
!!    ------
!
!     This routine is based on the following rule : 
!        one sent for one received (if it is possible);
!     The algorithm is the following :
!
!     while (there is some messages to send)
!       OR  (there is some messages to receive)
!     do
!       if (there is some messages to send)
!         send the first of the list
!       end if
!       if (there is some messages to receive)
!         try to receive
!       end if
!     done
!
!     The receptions are non-blocking and don't follow necessarly
!     the order of the TPCRSPDRECV list
!
!!    External
!!    -------
!     Module MODE_TOOLS_ll
!       GET_MAX_SIZE
!
!     Module MODE_EXCHANGE_ll
!       FILLIN_BUFFERS, FILLOUT_BUFFERS
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type CRSPD_ll, ZONE_ll
!
!     Module MODD_ARGSLIST_ll
!       type LIST_ll
!
!     Module MODD_VAR_ll
!       IP - Number of local processor=subdomain
!       NCOMBUFFSIZE1 - buffer size
!       NTRANS_COM - mpi communicator
!       MPI_PRECISION - mpi precision
!       NNEXTTAG, NMAXTAG - variable to define message tag
!
!     Module MODD_PARAMETERS_ll
!       JPVEXT - vertical halo size
!
!     Module MODD_DIM_ll
!       NKMAX_ll - maximum vertical dimension
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     N. Gicquel               * CERFACS - CNRM *
!
!-------------------------------------------------------------------------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll, ZONE_ll
  USE MODD_VAR_ll, ONLY : NCOMBUFFSIZE1, IP, MPI_PRECISION, NNEXTTAG, NMAXTAG
  USE MODD_DIM_ll, ONLY : NKMAX_ll
  USE MODD_PARAMETERS_ll, ONLY : JPVEXT
!
  USE MODE_TOOLS_ll, ONLY : GET_MAX_SIZE
!
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPCRSPDSEND, TPCRSPDRECV
  TYPE(LIST_ll), POINTER :: TPFIELDLISTSEND, TPFIELDLISTRECV
  INTEGER :: KMPI_COMM
  INTEGER :: KINFO
  INTEGER, OPTIONAL :: KBARRIER
!
!*       0.2   declarations of local variables
!
  INTEGER :: JINC, JI, JJ, JK  ! Loop and counter variables
  INTEGER :: FOUND, KERROR
!
#if defined (MNH_MPI_ISEND)
  REAL, DIMENSION (:,:), ALLOCATABLE :: TZBUFFER
#else      
  REAL, DIMENSION(:), ALLOCATABLE :: TZBUFFER  ! Buffers for info received
#endif
! JUAN
!
  INTEGER IRECVSTATUS(MPI_STATUS_SIZE) ! Status of completed receive request
!
  LOGICAL :: IRECVFLAG, ISENDFLAG
  INTEGER :: IMSGTAG, ISENDERPROC
!
  INTEGER :: IRECVNB, ISENDNB ! Total numbers of receive and send to do
  INTEGER :: IRECVDONE ! RECEIVE COMPLETED (receive and treated)
  INTEGER :: IMAXSIZESEND, IMAXSIZERECV, IBUFFSIZE
!
  TYPE(CRSPD_ll), POINTER :: TPMAILSEND, TPMAILRECV
  TYPE(ZONE_ll), POINTER :: TZZONESEND, TZZONERECV
  TYPE(LIST_ll), POINTER :: TZFIELDLISTSEND, TZFIELDLISTRECV
  INTEGER, SAVE :: ITAGOFFSET = 0
! JUAN
#if defined (MNH_MPI_ISEND)
INTEGER,PARAMETER                                     :: MPI_MAX_REQ = 1024
INTEGER,SAVE,DIMENSION(MPI_MAX_REQ)                   :: REQ_TAB
INTEGER,SAVE,DIMENSION(MPI_STATUS_SIZE,MPI_MAX_REQ)   :: STATUS_TAB
INTEGER                                               :: NB_REQ,NFIRST_REQ_RECV
#endif
! JUAN
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALISATIONS
!              ---------------
!
  IF ((.NOT.ASSOCIATED(TPFIELDLISTSEND)) &
    .OR.(.NOT.ASSOCIATED(TPFIELDLISTRECV))) THEN
    RETURN
  ENDIF
!
  IF (.NOT.ASSOCIATED(TPCRSPDSEND)) THEN
    ISENDNB = 0
    IMAXSIZESEND = 0
  ELSE
    ISENDNB = TPCRSPDSEND%NCARDDIF
    IMAXSIZESEND = GET_MAX_SIZE(TPCRSPDSEND) * TPFIELDLISTSEND%NCARD
  ENDIF
!
  IF (.NOT.ASSOCIATED(TPCRSPDRECV)) THEN
    IRECVNB = 0
    IMAXSIZERECV = 0
  ELSE
    IRECVNB = TPCRSPDRECV%NCARDDIF
    IMAXSIZERECV = GET_MAX_SIZE(TPCRSPDRECV) * TPFIELDLISTRECV%NCARD
  ENDIF
!
  IBUFFSIZE = IMAXSIZESEND 
  IF (IMAXSIZERECV > IBUFFSIZE) IBUFFSIZE = IMAXSIZERECV
!
  IBUFFSIZE = IBUFFSIZE * (NKMAX_ll + 2 * JPVEXT) 
! JUAN
#if defined (MNH_MPI_ISEND)
  ALLOCATE(TZBUFFER(IBUFFSIZE,ISENDNB+IRECVNB))
  TZBUFFER(:,:) = 0.0
#else
  ALLOCATE(TZBUFFER(IBUFFSIZE))
  TZBUFFER(:) = 0.0
#endif
! JUAN
! 
  TZFIELDLISTSEND => TPFIELDLISTSEND
  TZFIELDLISTRECV => TPFIELDLISTRECV
!
  IRECVDONE = 0
  TPMAILRECV => TPCRSPDRECV
  TPMAILSEND => TPCRSPDSEND
!
  IF (.NOT.PRESENT(KBARRIER)) THEN
!NZJUAN    CALL MPI_BARRIER(KMPI_COMM, KERROR)
  ELSE
!NZJUAN    IF (KBARRIER.EQ.1)  CALL MPI_BARRIER(KMPI_COMM, KERROR)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    MAIN LOOP
!              ---------
! JUAN
#if defined (MNH_MPI_ISEND)
  NB_REQ = 0
#endif
! JUAN
  DO WHILE (ASSOCIATED(TPMAILSEND))
!
!*       2.1  if there is still something to send
!
    IF (ASSOCIATED(TPMAILSEND)) THEN
  ! Build the send buffer
       TZZONESEND => TPMAILSEND%TELT 
       IF (TZZONESEND%NUMBER /= IP) THEN 
          JINC = 0
! JUAN
#if defined (MNH_MPI_ISEND)
          NB_REQ = NB_REQ + 1
          CALL FILLIN_BUFFERS(TZFIELDLISTSEND, TZZONESEND, TZBUFFER(:,NB_REQ), JINC)
#else
          CALL FILLIN_BUFFERS(TZFIELDLISTSEND, TZZONESEND, TZBUFFER, JINC)
#endif
! JUAN 
#if defined(MNH_MPI_BSEND)
          CALL MPI_BSEND(TZBUFFER, JINC, MPI_PRECISION, TZZONESEND%NUMBER - 1,  &
               TZZONESEND%MSSGTAG + ITAGOFFSET, KMPI_COMM,  KERROR)
#else
! JUAN
#if defined (MNH_MPI_ISEND)
          CALL MPI_ISEND(TZBUFFER(1,NB_REQ), JINC, MPI_PRECISION, TZZONESEND%NUMBER - 1,  &
               TZZONESEND%MSSGTAG + ITAGOFFSET, KMPI_COMM, REQ_TAB(NB_REQ), KERROR)
#else
          CALL MPI_SEND(TZBUFFER, JINC, MPI_PRECISION, TZZONESEND%NUMBER - 1,  &
               TZZONESEND%MSSGTAG + ITAGOFFSET, KMPI_COMM, KERROR)
#endif
! JUAN
#endif
      ENDIF
      TPMAILSEND => TPMAILSEND%TNEXT
    ENDIF
  ENDDO

!NZJUAN  CALL MPI_BARRIER(KMPI_COMM, KERROR)

! JUAN
#if defined (MNH_MPI_ISEND)
  NFIRST_REQ_RECV = NB_REQ
#endif
! JUAN 

  DO WHILE (ASSOCIATED(TPMAILRECV)) 
     IF (TPMAILRECV%TELT%NUMBER == IP) THEN
        TPMAILRECV => TPMAILRECV%TNEXT
     ELSE
#if defined (MNH_MPI_ISEND)
        NB_REQ = NB_REQ + 1
        CALL MPI_IRECV(TZBUFFER(1,NB_REQ), IBUFFSIZE, MPI_PRECISION, &
             TPMAILRECV%TELT%NUMBER -1 , &
             TPMAILRECV%TELT%MSSGTAG + ITAGOFFSET, &
             KMPI_COMM, REQ_TAB(NB_REQ), KERROR)
#else
        CALL MPI_RECV(TZBUFFER, IBUFFSIZE, MPI_PRECISION, &
             TPMAILRECV%TELT%NUMBER -1 , &
             TPMAILRECV%TELT%MSSGTAG + ITAGOFFSET, &
             KMPI_COMM, IRECVSTATUS, KERROR)
        JINC = 0 
        CALL FILLOUT_BUFFERS(TZFIELDLISTRECV, TPMAILRECV%TELT, TZBUFFER, JINC)
#endif
! JUAN
        TPMAILRECV => TPMAILRECV%TNEXT
        !
     ENDIF
     
     !
  ENDDO

! JUAN
#if defined (MNH_MPI_ISEND)
  CALL MPI_WAITALL(NB_REQ,REQ_TAB,STATUS_TAB,KINFO) 

  TPMAILRECV => TPCRSPDRECV
  NB_REQ = NFIRST_REQ_RECV

  DO WHILE (ASSOCIATED(TPMAILRECV)) 
     IF (TPMAILRECV%TELT%NUMBER == IP) THEN
        TPMAILRECV => TPMAILRECV%TNEXT
     ELSE 
        !
        NB_REQ = NB_REQ + 1
        JINC = 0 
        CALL FILLOUT_BUFFERS(TZFIELDLISTRECV, TPMAILRECV%TELT, TZBUFFER(:,NB_REQ), JINC)
        TPMAILRECV => TPMAILRECV%TNEXT
        !
     ENDIF
     
     !
  ENDDO
#endif
!JUAN 
!
  DEALLOCATE(TZBUFFER)
!  
  IF (.NOT.PRESENT(KBARRIER)) THEN 
!NZJUAN    CALL MPI_BARRIER(KMPI_COMM, KERROR)
  ELSE
  IF (KBARRIER.EQ.1)  CALL MPI_BARRIER(KMPI_COMM, KERROR)
  ENDIF
!
  ITAGOFFSET = MOD((ITAGOFFSET + NNEXTTAG), NMAXTAG)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SEND_RECV_CRSPD
!
END MODULE MODE_EXCHANGE_ll
