!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------

!     #################
      MODULE MODE_LB_ll
!     #################
! 
!!    Purpose
!!    -------
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!     SUBROUTINES : SET_LBFIELD_ll (SET_LB2DFIELD_ll, SET_LB3DFIELD_ll)
!                   UNSET_LBFIELD_ll,
!                   LB_FORCING_ll,
!                   SET_LBSIZEX_ll, SET_LBSIZEY_ll
!                   INIT_LB_ll
!
!!    Reference
!!    --------- 
!
!     User Interface for Meso-NH parallel package
!     Ph. Kloos, L. Giraud, R. Guivarch, D. Lugato
!
!!    Authors
!!    -------
!
!     R. Guivarch    * CERFACS *
!     Ph. Kloos      * CNRM - CERFACS *
!     N. Gicquel     * CNRM - CERFACS *
!
!!   Implicit Arguments
!!   ------------------
!
!     Module MODD_DIM_ll
!       NDXRATIO_ALL, NDYRATIO_ALL - Ratio for all models
!       NXOR_ALL, NXEND_ALL, NYOR_ALL, NYEND_ALL, NDAD
!
!     Module MODD_STRUCTURE_ll
!       types LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll,
!             PARENT2CHILD_DATA_ll, LPROCONF_ll, PROC_COM_DATA_ll,
!             PROCONF_ll, ZONE_ll, CRSPD_ll
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPMODELMAX, NMAXRIM
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NGRID_COM - mpi communicator
!       NPROC, DIMX, DIMY, IP, YSPLITTING
!
!------------------------------------------------------------------------------
!  
  USE MODD_STRUCTURE_ll
!
  CONTAINS
!
!     #############################################################
      SUBROUTINE SET_LB2DFIELD_ll(P2DFIELD, PTFIELD, KFINELBSIZE, &
                                  HSIDE, KMODEL)
!     #############################################################
!
!!****  *SET_LB2DFIELD_ll*- routine to fill the parallel data lists
!                           of forcing
!!    Purpose
!!    -------
!     this routine adds the 2D field P2DFIELD (resp. PTFIELD)
!     to the forcing exchange parent (resp. child) list
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!     Module MODE_ARGSLIST_ll
!       ADD2DFIELD_ll
!
!     Module MODE_NEST_ll
!       LBFINE2COARSE
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_DIM_ll
!       NDXRATIO_ALL, NDYRATIO_ALL - Ratio for all models
!
!     Module MODD_STRUCTURE_ll
!       types LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : NDXRATIO_ALL, NDYRATIO_ALL
  USE MODD_STRUCTURE_ll, ONLY : LCRSPD_ll, LPROC_COM_DATA_ll, &
                                LPARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  USE MODE_ARGSLIST_ll, ONLY : ADD2DFIELD_ll
  USE MODE_NEST_ll, ONLY : LBFINE2COARSE
!  
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
  INTEGER, INTENT(IN) :: KFINELBSIZE, KMODEL
  CHARACTER(LEN=*), INTENT(IN) :: HSIDE
!
!*       0.2   declarations of local variables
!
  INTEGER :: ICOARSE
  TYPE(LCRSPD_ll), POINTER :: TZPAR, TZCHILD
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZLCOMDATA
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZP2CDATA
!
!-------------------------------------------------------------------------------
!
!*       1.    Initialisation
!
!*       1.1   Checks
!
  IF (.NOT.ASSOCIATED(TCRRT_COMDATA%TCHILDREN) &
    & .OR. .NOT.ASSOCIATED(TCRRT_COMDATA%TP2C_DATA)) THEN
    WRITE(*,*) 'Problem in set_lbfield_ll'
    WRITE(*,*) 'The current model has no child'
    RETURN
  ENDIF
!
!*       1.    Calculate the coarse size needed for the interpolation
!
  IF (HSIDE == "NORTH" .OR. HSIDE == "SOUTH") THEN
    ICOARSE = LBFINE2COARSE(NDYRATIO_ALL(KMODEL), KFINELBSIZE)
  ELSE
    ICOARSE = LBFINE2COARSE(NDXRATIO_ALL(KMODEL), KFINELBSIZE)
  ENDIF
!-------------------------------------------------------------------------------
!
!*       2.    Point to the appropriate communication structures
!
!*       2.1   Point to the child communication data structure
!
  TZLCOMDATA => TCRRT_COMDATA%TCHILDREN
  DO WHILE (TZLCOMDATA%TELT%NUMBER /= KMODEL)
    TZLCOMDATA => TZLCOMDATA%TNEXT
  ENDDO
  IF (.NOT.ASSOCIATED(TZLCOMDATA)) THEN
    WRITE(*,*) 'Error SET_LBFIELD_ll : ', KMODEL, &
               ' is not a child of the current model'
    STOP
  ENDIF
!
!*       2.2   Point to the parent2child data structure
!
  TZP2CDATA => TCRRT_COMDATA%TP2C_DATA
  DO WHILE (TZP2CDATA%TELT%NUMBER /= KMODEL)
    TZP2CDATA => TZP2CDATA%TNEXT
  ENDDO
  IF (.NOT.ASSOCIATED(TZP2CDATA)) THEN
    WRITE(*,*) 'Error SET_LBFIELD_ll : ', KMODEL, &
               ' is not a child of the current model'
    STOP
  ENDIF
!
!*       2.3   Point to the appropriate side
!
  SELECT CASE(HSIDE)
    CASE("NORTH")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBYN
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBYN
    CASE("SOUTH")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBYS
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBYS
    CASE("WEST")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBXW
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBXW
    CASE("EAST")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBXE
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBXE
  END SELECT
!
!-------------------------------------------------------------------------------
!
!*       3.    Point to the LCRSPD corresponding to 
!*             the coarse size needed for the interpolation
!*             and put P2DFIELD and PTFIELD in the list of fields
!
  DO WHILE (ASSOCIATED(TZCHILD))
    IF (TZCHILD%NWIDTH == ICOARSE) THEN
      CALL ADD2DFIELD_ll(TZCHILD%TLIST, PTFIELD)
      EXIT
    ENDIF
    TZCHILD => TZCHILD%TNEXT
  ENDDO
!
  DO WHILE (ASSOCIATED(TZPAR)) 
    IF (ICOARSE == TZPAR%NWIDTH) THEN
      CALL ADD2DFIELD_ll(TZPAR%TLIST, P2DFIELD)
      EXIT
    ENDIF
    TZPAR => TZPAR%TNEXT
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LB2DFIELD_ll
!
!     #############################################################
      SUBROUTINE SET_LB3DFIELD_ll(P3DFIELD, PTFIELD, KFINELBSIZE, &
                                  HSIDE, KMODEL)
!     #############################################################
!
!!****  *SET_LB3DFIELD_ll*- routine to fill the parallel data lists
!                           of forcing
!!    Purpose
!!    -------
!     this routine adds the 3D field P3DFIELD (resp. PTFIELD)
!     to the forcing exchange parent (resp. child) list
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!     Module MODE_ARGSLIST_ll
!       ADD3DFIELD_ll
!
!     Module MODE_NEST_ll
!       LBFINE2COARSE
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_DIM_ll
!       NDXRATIO_ALL, NDYRATIO_ALL - Ratio for all models
!
!     Module MODD_STRUCTURE_ll
!       types LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     P. Kloos
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : NDXRATIO_ALL, NDYRATIO_ALL
  USE MODD_STRUCTURE_ll, ONLY : LCRSPD_ll, LPROC_COM_DATA_ll, &
                                LPARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  USE MODE_ARGSLIST_ll, ONLY : ADD3DFIELD_ll
  USE MODE_NEST_ll, ONLY : LBFINE2COARSE
!
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
  INTEGER, INTENT(IN) :: KFINELBSIZE, KMODEL
  CHARACTER(LEN=*), INTENT(IN) :: HSIDE
!
!*       0.2   declarations of local variables
!
  INTEGER :: ICOARSE
  TYPE(LCRSPD_ll), POINTER :: TZPAR, TZCHILD
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZLCOMDATA
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZP2CDATA
!
!-------------------------------------------------------------------------------
!
!*       1.    Initialisation
!
!*       1.1   Checks
!
  IF (.NOT.ASSOCIATED(TCRRT_COMDATA%TCHILDREN) &
    & .OR. .NOT.ASSOCIATED(TCRRT_COMDATA%TP2C_DATA)) THEN
    WRITE(*,*) 'Problem in set_lbfield_ll'
    WRITE(*,*) 'The current model has no child'
    RETURN
  ENDIF
!
!*       1.    Calculate the coarse size needed for the interpolation
!
  IF (HSIDE == "NORTH" .OR. HSIDE == "SOUTH") THEN
    ICOARSE = LBFINE2COARSE(NDYRATIO_ALL(KMODEL), KFINELBSIZE)
  ELSE
    ICOARSE = LBFINE2COARSE(NDXRATIO_ALL(KMODEL), KFINELBSIZE)
  ENDIF
!-------------------------------------------------------------------------------
!
!*       2.    Point to the appropriate communication structures
!
!*       2.1   Point to the child communication data structure
!
  TZLCOMDATA => TCRRT_COMDATA%TCHILDREN
  DO WHILE (TZLCOMDATA%TELT%NUMBER /= KMODEL)
    TZLCOMDATA => TZLCOMDATA%TNEXT
  ENDDO
  IF (.NOT.ASSOCIATED(TZLCOMDATA)) THEN
    WRITE(*,*) 'Error SET_LBFIELD_ll : ', KMODEL, &
               ' is not a child of the current model'
    STOP
  ENDIF
!
!*       2.2   Point to the parent2child data structure
!
  TZP2CDATA => TCRRT_COMDATA%TP2C_DATA
  DO WHILE (TZP2CDATA%TELT%NUMBER /= KMODEL)
    TZP2CDATA => TZP2CDATA%TNEXT
  ENDDO
  IF (.NOT.ASSOCIATED(TZP2CDATA)) THEN
    WRITE(*,*) 'Error SET_LBFIELD_ll : ', KMODEL, &
               ' is not a child of the current model'
    STOP
  ENDIF
!
!*       2.3   Point to the appropriate side
!
  SELECT CASE(HSIDE)
    CASE("NORTH")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBYN
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBYN
    CASE("SOUTH")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBYS
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBYS
    CASE("WEST")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBXW
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBXW
    CASE("EAST")
      TZCHILD => TZLCOMDATA%TELT%TRECV_LBXE
      TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LBXE
  END SELECT
!
!-------------------------------------------------------------------------------
!
!*       3.    Point to the LCRSPD corresponding to 
!*             the coarse size needed for the interpolation
!*             and put P3DFIELD and PTFIELD in the list of fields
!
  DO WHILE (ASSOCIATED(TZCHILD))
    IF (TZCHILD%NWIDTH == ICOARSE) THEN
      CALL ADD3DFIELD_ll(TZCHILD%TLIST, PTFIELD)
      EXIT
    ENDIF
    TZCHILD => TZCHILD%TNEXT
  ENDDO
!
  DO WHILE (ASSOCIATED(TZPAR)) 
    IF (ICOARSE == TZPAR%NWIDTH) THEN
      CALL ADD3DFIELD_ll(TZPAR%TLIST, P3DFIELD)
      EXIT
    ENDIF
    TZPAR => TZPAR%TNEXT
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LB3DFIELD_ll
!
!     #############################
      SUBROUTINE UNSET_LBFIELD_ll()
!     #############################
!
!!****  *UNSET_LBFIELD_ll*- routine to clean the parallel data lists
!                           of forcing
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!     Module MODE_CONSTRUCT_ll
!       CLEANLIST_LCRSPD
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       types LPARENT2CHILD_DATA_ll, PARENT2CHILD_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     P. Kloos
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!!
  USE MODD_STRUCTURE_ll, ONLY : LPARENT2CHILD_DATA_ll, PARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
  USE MODE_CONSTRUCT_ll, ONLY : CLEANLIST_LCRSPD
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!*       0.2   declarations of local variables
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZP2C
  TYPE(PARENT2CHILD_DATA_ll), POINTER :: TZP2CDATA
!
!-------------------------------------------------------------------------------
!
  CALL CLEANLIST_LCRSPD(TCRRT_COMDATA%TRECV_LBXW)
  CALL CLEANLIST_LCRSPD(TCRRT_COMDATA%TRECV_LBXE)
  CALL CLEANLIST_LCRSPD(TCRRT_COMDATA%TRECV_LBYS)
  CALL CLEANLIST_LCRSPD(TCRRT_COMDATA%TRECV_LBYN)
!
  IF (TCRRT_COMDATA%NUMBER > 1) THEN
    TZP2C => TCRRT_COMDATA%TPARENT%TP2C_DATA
!
! Find the parent's PARENT2CHILD_DATA_ll data structure
! associated to the current model
!
    DO WHILE(ASSOCIATED(TZP2C))
      IF (TZP2C%TELT%NUMBER == TCRRT_COMDATA%NUMBER) THEN
        EXIT
      ENDIF
      TZP2C => TZP2C%TNEXT
    ENDDO
!
    TZP2CDATA => TZP2C%TELT
!
! Clean the lists corresponding to the parent's send for the LB
!
    CALL CLEANLIST_LCRSPD(TZP2CDATA%TSEND_1WAY_LBXW)
    CALL CLEANLIST_LCRSPD(TZP2CDATA%TSEND_1WAY_LBXE)
    CALL CLEANLIST_LCRSPD(TZP2CDATA%TSEND_1WAY_LBYS)
    CALL CLEANLIST_LCRSPD(TZP2CDATA%TSEND_1WAY_LBYN)
!
  ELSE
!
    WRITE(*,*) 'Problem in UNSET_LBFIELD'
    WRITE(*,*) 'The current model is 1'
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UNSET_LBFIELD_ll
!
!     #########################################
      SUBROUTINE LB_FORCING_ll( KCHILD, KINFO )
!     #########################################
!
!!****  *LB_FORCING_ll* - routine to do the forcing
!
!!    Purpose
!!    -------
!     This routine fill the fields of the child model with the values
!     of the parent model
!
!     DO NOT SWITCH TO THE CHILD MODEL
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!     Module MODE_EXCHANGE_ll
!       SEND_RECV_CRSPD, COPY_CRSPD
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       types LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll,
!             LPROCONF_ll, PROC_COM_DATA_ll,
!             PARENT2CHILD_DATA_ll, LCRSPD_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NGRID_COM - mpi communicator
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
!     Original 14/02/00
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll, &
                                LPROCONF_ll, PROC_COM_DATA_ll,           &
                                PARENT2CHILD_DATA_ll, LCRSPD_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, TCRRT_PROCONF, NGRID_COM
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_CRSPD, COPY_CRSPD
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN) :: KCHILD
!
  INTEGER, INTENT(OUT) :: KINFO
!
!*       0.2   declarations of local variables
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZCURRENT ! intermediate variables
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZCHILDCURRENT
  TYPE(LPROCONF_ll), POINTER :: TZPROCONFCURRENT
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TZCHILD_COMDATA    ! child
  TYPE(PARENT2CHILD_DATA_ll), POINTER :: TZP2C_DATA
!
  TYPE(LIST_ll), POINTER :: TZPARENTFIELD
  TYPE(LIST_ll), POINTER :: TZCHILDFIELD
!
  TYPE(LCRSPD_ll), POINTER :: TZCHILDLCRSPD, TZPARLCRSPD
!
!-------------------------------------------------------------------------------
!
!*        1.   POINT TO THE CHILD COMDATA, PROCONF
!*             POINT TO PARENT2CHILDDATA
!              -----------------------------------
!
  KINFO = -1
  TZCURRENT => TCRRT_COMDATA%TCHILDREN
  TZCHILDCURRENT => TCRRT_COMDATA%TP2C_DATA
  TZPROCONFCURRENT => TCRRT_PROCONF%TCHILDREN
!
  DO WHILE(ASSOCIATED(TZCURRENT))
    IF(TZCURRENT%TELT%NUMBER .EQ. KCHILD) THEN
      KINFO = 0
      EXIT
    ENDIF
    TZCURRENT => TZCURRENT%TNEXT
    TZCHILDCURRENT => TZCHILDCURRENT%TNEXT
    TZPROCONFCURRENT => TZPROCONFCURRENT%TNEXT
  ENDDO
!
  IF(.NOT.ASSOCIATED(TZCURRENT)) RETURN
!
  TZCHILD_COMDATA => TZCURRENT%TELT
  TZP2C_DATA => TZCHILDCURRENT%TELT
!
!-------------------------------------------------------------------------------
!
!*        2.   SEND AND RECEIVE THE LB FIELDS
!              ------------------------------
!
!*        2.1  West side
!
  TZCHILDLCRSPD => TZCHILD_COMDATA%TRECV_LBXW
  TZPARLCRSPD => TZP2C_DATA%TSEND_1WAY_LBXW
!
  DO WHILE(ASSOCIATED(TZCHILDLCRSPD) .OR. ASSOCIATED(TZPARLCRSPD))
!
    CALL SEND_RECV_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                         TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                         NGRID_COM, KINFO, 0)  ! 0 Prevent CALL TO MPI BARRIER
                                               ! (No symetric communications) 
!
    CALL COPY_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                    TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                    KINFO)
!
    TZPARLCRSPD => TZPARLCRSPD%TNEXT
    TZCHILDLCRSPD => TZCHILDLCRSPD%TNEXT
!
  ENDDO
!
!*        2.2  East side
!
  TZCHILDLCRSPD => TZCHILD_COMDATA%TRECV_LBXE
  TZPARLCRSPD => TZP2C_DATA%TSEND_1WAY_LBXE
!
  DO WHILE(ASSOCIATED(TZCHILDLCRSPD) .OR. ASSOCIATED(TZPARLCRSPD))
!
    CALL SEND_RECV_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                         TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                         NGRID_COM, KINFO, 0)  ! 0 Prevent CALL TO MPI BARRIER
                                               ! (No symetric communications)
!
    CALL COPY_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                    TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                    KINFO)
!
    TZPARLCRSPD => TZPARLCRSPD%TNEXT
    TZCHILDLCRSPD => TZCHILDLCRSPD%TNEXT
!
  ENDDO
!
!*        2.3  South side
!
  TZCHILDLCRSPD => TZCHILD_COMDATA%TRECV_LBYS
  TZPARLCRSPD => TZP2C_DATA%TSEND_1WAY_LBYS
!
  DO WHILE(ASSOCIATED(TZCHILDLCRSPD) .OR. ASSOCIATED(TZPARLCRSPD))
!
    CALL SEND_RECV_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                         TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                         NGRID_COM, KINFO, 0)
!
    CALL COPY_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                    TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                    KINFO)
!
    TZPARLCRSPD => TZPARLCRSPD%TNEXT
    TZCHILDLCRSPD => TZCHILDLCRSPD%TNEXT
!
  ENDDO
!
!*        2.4  North side
!
  TZCHILDLCRSPD => TZCHILD_COMDATA%TRECV_LBYN
  TZPARLCRSPD => TZP2C_DATA%TSEND_1WAY_LBYN
!
  DO WHILE(ASSOCIATED(TZCHILDLCRSPD) .OR. ASSOCIATED(TZPARLCRSPD))
!
    CALL SEND_RECV_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                         TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                         NGRID_COM, KINFO, 0)
!
    CALL COPY_CRSPD(TZPARLCRSPD%TCRSPD, TZCHILDLCRSPD%TCRSPD, &
                    TZPARLCRSPD%TLIST, TZCHILDLCRSPD%TLIST, &
                    KINFO)
!
    TZPARLCRSPD => TZPARLCRSPD%TNEXT
    TZCHILDLCRSPD => TZCHILDLCRSPD%TNEXT
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE LB_FORCING_ll
!
!     ##########################################
      SUBROUTINE SET_LBSIZEX_ll(KNBRIM, KRIMTAB)
!     ##########################################
!
!!****  *SET_LBSIZEX_ll*- routine to initialize the array of rim 
!                         in x-direction of the current  child model
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     N. Gicquel      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!*       0.1  declarations of arguments
!
  INTEGER :: KNBRIM
  INTEGER, DIMENSION(:) :: KRIMTAB
!
!-------------------------------------------------------------------------------
!
  TCRRT_COMDATA%NDIMRIMLBX = KNBRIM
  TCRRT_COMDATA%NRIMLBX(1:KNBRIM) = KRIMTAB(1:KNBRIM)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LBSIZEX_ll
!
!     ##########################################
      SUBROUTINE SET_LBSIZEY_ll(KNBRIM, KRIMTAB)
!     ##########################################
!
!!****  *SET_LBSIZEY_ll*- routine to initialize the array of rim 
!                         in y-direction of the current  child model
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     N. Gicquel      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!*       0.1  declarations of arguments
!
  INTEGER :: KNBRIM
  INTEGER, DIMENSION(:) :: KRIMTAB
!
!-------------------------------------------------------------------------------
!
  TCRRT_COMDATA%NDIMRIMLBY = KNBRIM
  TCRRT_COMDATA%NRIMLBY(1:KNBRIM) = KRIMTAB(1:KNBRIM)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LBSIZEY_ll
!
!     ##################################################################
      RECURSIVE SUBROUTINE INIT_LB_MODEL (TPPROCONF, TPCOMDATA, KMODEL ) 
!     ##################################################################
!
!!****  *INIT_LB_MODEL* - routine to construct recursively the data structures
!                          for LB FIELDS of all models
! 
!!    Purpose
!!    -------
!     this routine fills all the data structures (splitting data structures,
!     communications data structures) for LB FIELDS
!
!!**  Method
!!    ------
! 
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll, EXTRACT_ZONE
!
!     Module MODE_SPLITTING_ll
!       SPLIT2
!
!     Module MODE_NEST_ll
!       CHILDS_PT_OF_VIEW, LBFINE2COARSE, PARENTS_PT_OF_VIEW 
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       types PROC_COM_DATA_ll, PROCONF_ll, LPROCONF_ll,
!             PARENT2CHILD_DATA_ll, LPARENT2CHILD_DATA_ll,
!             LPROC_COM_DATA_ll, ZONE_ll, CRSPD_ll, LCRSPD_ll
!
!     Module MODD_PARAMETERS_ll
!       JPHEXT, JPMODELMAX, NMAXRIM
!
!     Module MODD_VAR_ll
!       NPROC, DIMX, DIMY, IP, YSPLITTING, NIOUNIT
!
!     Module MODDMODD_DIM_ll
!       NDXRATIO_ALL, NXOR_ALL, NXEND_ALL
!       NDYRATIO_ALL, NYOR_ALL, NYEND_ALL, NDAD
! 
!!    Reference
!!    ---------
!     Notes sur les routines multi-modeles de la surcouche, R. Guivarch
! 
!!    Author
!!    ------
!     N. Gicquel      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!     Philippe 12/01/2018: renamed dimension variables NKMAX_ll->NKMAX_TMP_ll to prevent mix-up with modd_dimn
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : NDXRATIO_ALL, NXOR_ALL, NXEND_ALL, &
                          NDYRATIO_ALL, NYOR_ALL, NYEND_ALL, NDAD,&
                          NKMAX_TMP_ll
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPMODELMAX, NMAXRIM
  USE MODD_STRUCTURE_ll, ONLY : PROC_COM_DATA_ll, PROCONF_ll, LPROCONF_ll,   &
                                PARENT2CHILD_DATA_ll, LPARENT2CHILD_DATA_ll, &
                                LPROC_COM_DATA_ll, ZONE_ll, CRSPD_ll,        &
                                LCRSPD_ll
  USE MODD_VAR_ll, ONLY : NPROC, DIMX, DIMY, IP, YSPLITTING, NIOUNIT
!
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll, &
                            EXTRACT_ZONE
  USE MODE_SPLITTING_ll, ONLY : SPLIT2
  USE MODE_NEST_ll, ONLY : CHILDS_PT_OF_VIEW, LBFINE2COARSE, PARENTS_PT_OF_VIEW
!
!  USE MODE_TEST_ll, ONLY : DISPLAY_ZONE, DISPLAY_SPLITTING, &
!                           DISPLAY_CRSPD, DISPLAY_LCRSPD
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data structure
  TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
  INTEGER :: KMODEL
!
!*       0.2   declarations of local variables
!
  TYPE(LPROCONF_ll), POINTER :: TZLCHILD_PROCONF
!
  TYPE(PROCONF_ll), POINTER :: TZCHILD_PROCONF ! variable for the child
                                               ! splitting data
!
  TYPE(PARENT2CHILD_DATA_ll), POINTER :: TZP2C_DATA ! variable for parent
                                                    ! to child communications
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZLP2C_DATA
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZLCHILD_COMDATA
  TYPE(PROC_COM_DATA_ll), POINTER :: TZCHILD_COMDATA ! variable for the child
                                                     ! communication data
                                                     ! (parent to child, halo,
                                                     ! transposition)
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZCOARSE ! Coarse grid splitting
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZFINE, &! Fine grid splitting
             				      TZFINEPB ! tzfine + boundary halo
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZINTER ! Intermediate zone array
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZP ! Physical zone splitting
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZPB ! Physical + boundary halo
                                                   ! splitting
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZLB ! LB points of the
                                                   ! coarse grid to be
						   ! send to the child
  TYPE(ZONE_ll) :: TZCHILDLB
!
  TYPE(LCRSPD_ll), POINTER :: TZLCRSPD, TZLTMPCRSPD, & ! temporary lists
                              TZLCRSPDS, TZLCRSPDR     ! of correspondants
!
  TYPE(CRSPD_ll), POINTER :: TZCRSPDS, TZCRSPDR
!
  INTEGER :: J, JMODEL, JI ! loop controls
!
  INTEGER :: ICOARSELBX, ICOARSELBY ! size of the LB zones
  INTEGER :: IPREV ! used for the LB zones
!
  INTEGER :: MESSTAG
!
  INTEGER :: IWLBMSG, IELBMSG, INLBMSG, ISLBMSG ! message tags for LB
!
  INTEGER :: IDIMRIMLBX, IDIMRIMLBY ! number of sizes of fine RIM zones
  INTEGER, DIMENSION(NMAXRIM) :: IRIMLBX, IRIMLBY ! sizes of the fine RIM zones
  INTEGER :: IERR
  INTEGER :: ICARD
  LOGICAL :: GINTER
!
!-------------------------------------------------------------------------------
!
  IWLBMSG = 210
  IELBMSG = 220
  ISLBMSG = 230
  INLBMSG = 240
!
!-------------------------------------------------------------------------------
!
!*        1.   GET PROCONF AND PROCCOMDATA OF THE CHILDREN
!              -------------------------------------------
!
  TZLCHILD_PROCONF => TPPROCONF%TCHILDREN
!
  DO WHILE (TZLCHILD_PROCONF%TELT%NUMBER /= KMODEL)
     TZLCHILD_PROCONF => TZLCHILD_PROCONF%TNEXT
  ENDDO
!
  TZCHILD_PROCONF => TZLCHILD_PROCONF%TELT
!
  TZLCHILD_COMDATA => TPCOMDATA%TCHILDREN
!
  DO WHILE (TZLCHILD_COMDATA%TELT%NUMBER /= KMODEL)
     TZLCHILD_COMDATA => TZLCHILD_COMDATA%TNEXT
  ENDDO
!
  TZCHILD_COMDATA => TZLCHILD_COMDATA%TELT
!
!-------------------------------------------------------------------------------
!
!*        2.   ALLOCATE ARRAYS FOR TEMPORARY COMPUTATIONS
!              ------------------------------------------
!
  ALLOCATE(TZLB(NPROC))
  DIMX = NDXRATIO_ALL(KMODEL) * (NXEND_ALL(KMODEL)-NXOR_ALL(KMODEL) - 1) &
         + 2*JPHEXT
  DIMY = NDYRATIO_ALL(KMODEL) * (NYEND_ALL(KMODEL)-NYOR_ALL(KMODEL) - 1) &
         + 2*JPHEXT
!
!-------------------------------------------------------------------------------
!
!*       3.    CONSTRUCTION OF USEFUL SPLITTINGS :
!              ---------------------------------
  ALLOCATE(TZCOARSE(NPROC))
  CALL SPLIT2(NXEND_ALL(KMODEL)-NXOR_ALL(KMODEL)-1, &
              NYEND_ALL(KMODEL)-NYOR_ALL(KMODEL)-1, &
              NKMAX_TMP_ll,                         &
              NPROC,TZCOARSE, YSPLITTING)
!
!  WRITE(NIOUNIT,*) 'TZCOARSE',NPROC,KMODEL
!  CALL DISPLAY_SPLITTING(TZCOARSE,NPROC,1)
!
!        3.1   TZFINE splitting of fine grid of the child model
!                     (physical zones)
!
  ALLOCATE(TZFINE(NPROC))
  DO J = 1, NPROC
!
    TZFINE(J)%NUMBER = TZCOARSE(J)%NUMBER
    TZFINE(J)%NXOR = (TZCOARSE(J)%NXOR - 2) * NDXRATIO_ALL(KMODEL) + 2
    TZFINE(J)%NYOR = (TZCOARSE(J)%NYOR - 2) * NDYRATIO_ALL(KMODEL) + 2
    TZFINE(J)%NXEND = (TZCOARSE(J)%NXEND - 1) * NDXRATIO_ALL(KMODEL) + 1
    TZFINE(J)%NYEND = (TZCOARSE(J)%NYEND - 1) * NDYRATIO_ALL(KMODEL) + 1
!
  ENDDO
!
!  WRITE(NIOUNIT,*) 'TZFINE'
!  CALL DISPLAY_SPLITTING(TZFINE,NPROC,1)
!
!        3.2   TZFINEPB : splitting of the fine grid of the child model
!                         (physical zones + boundary halo)
!
  ALLOCATE(TZFINEPB(NPROC))
  TZFINEPB(1:NPROC) = TZFINE(1:NPROC)
  DO J=1, NPROC
!
    IF (LWEST_ll(J)) TZFINEPB(J)%NXOR = TZFINE(J)%NXOR - JPHEXT
    IF (LEAST_ll(J)) TZFINEPB(J)%NXEND = TZFINE(J)%NXEND + JPHEXT
    IF (LSOUTH_ll(J)) TZFINEPB(J)%NYOR = TZFINE(J)%NYOR - JPHEXT
    IF (LNORTH_ll(J)) TZFINEPB(J)%NYEND = TZFINE(J)%NYEND + JPHEXT
!
  ENDDO
!
!  WRITE(NIOUNIT,*) 'TZFINEPB'
!  CALL DISPLAY_SPLITTING(TZFINEPB,NPROC,1)
!
!        3.3   TZPB splitting of the physical zone + the boundary halo
!              TZP splitting of the physical zone
!
  ALLOCATE(TZP(NPROC), TZPB(NPROC))
!
  CALL EXTRACT_ZONE( TPPROCONF%TSPLITS_B, TZP, TZPB )
!
!  WRITE(NIOUNIT,*) 'TZP'
!  CALL DISPLAY_SPLITTING(TZP,NPROC,1)
!
  DO J = 1, NPROC
!
    IF(.NOT.LNORTH_ll(J)) TZPB(J)%NYEND = TZP(J)%NYEND
    IF(.NOT.LSOUTH_ll(J)) TZPB(J)%NYOR = TZP(J)%NYOR
    IF(.NOT.LWEST_ll(J))  TZPB(J)%NXOR = TZP(J)%NXOR
    IF(.NOT.LEAST_ll(J))  TZPB(J)%NXEND = TZP(J)%NXEND
!
  ENDDO
!
!  WRITE(NIOUNIT,*) 'TZPB'
!  CALL DISPLAY_SPLITTING(TZPB,NPROC,1)
!
  TZLP2C_DATA => TPCOMDATA%TP2C_DATA
!
  IF (.NOT.ASSOCIATED(TZLP2C_DATA)) RETURN
!
! Search the parent to child corresponding to KMODEL
!
  DO WHILE((TZLP2C_DATA%TELT%NUMBER /= KMODEL).AND.ASSOCIATED(TZLP2C_DATA))
     TZLP2C_DATA => TZLP2C_DATA%TNEXT
  ENDDO
!
  IF (.NOT.ASSOCIATED(TZLP2C_DATA)) RETURN
!
  TZP2C_DATA => TZLP2C_DATA%TELT
!   
!-------------------------------------------------------------------------------
!
!*        4.   CONSTRUCTION OF TZP2C_DATA%TSEND_1WAY_LBXW
!              ------------------------------------------
!
  ICARD = 0
  IPREV = 0
!  GINTER = .FALSE.
!
!  WRITE(NIOUNIT,*) '****** OUEST ******'
!  WRITE(NIOUNIT,*) 'NDIMRIMLBX = ', TZCHILD_COMDATA%NDIMRIMLBX
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBX
!
!         4.1  Calculate the width of the coarse zone to be interpolated
!
    ICOARSELBX = LBFINE2COARSE(NDXRATIO_ALL(KMODEL), &
                               TZCHILD_COMDATA%NRIMLBX(JI))
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBX(JI)', &
!                     TZCHILD_COMDATA%NRIMLBX(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBX =',ICOARSELBX
!
    IF (ICOARSELBX /= IPREV) THEN
!
!         4.2  Set up LB zones on each child proc
!
      CALL ALLOC(TZLTMPCRSPD)
      DO J = 1, NPROC
        IF (TZCHILD_COMDATA%NRIMLBX(JI)+JPHEXT >= TZFINEPB(J)%NXOR) THEN
          TZLB(J) = ZONE_ll(0,0,0,0,0,0,0,0)
          TZLB(J)%NUMBER = TZCOARSE(J)%NUMBER
          TZLB(J)%NXOR = NXOR_ALL(KMODEL) + TZCOARSE(J)%NXOR - 3
          TZLB(J)%NYOR = NYOR_ALL(KMODEL) + TZCOARSE(J)%NYOR - 3
          TZLB(J)%NXEND = MIN(NXOR_ALL(KMODEL) + TZCOARSE(J)%NXEND + 1, &
                            NXOR_ALL(KMODEL) + ICOARSELBX - 2)
          TZLB(J)%NYEND = NYOR_ALL(KMODEL) + TZCOARSE(J)%NYEND + 1
        ELSE
          TZLB(J) = ZONE_ll(0,0,0,0,0,0,0,0)
        ENDIF
!        CALL DISPLAY_ZONE(TZLB(J),1)
      ENDDO
!
      MESSTAG = IWLBMSG + ICOARSELBX
      CALL PARENTS_PT_OF_VIEW(TZLB, TZPB(IP), TPPROCONF, MESSTAG, &
                              TZLTMPCRSPD, GINTER)
!
      IF (GINTER) THEN
         TZLTMPCRSPD%NWIDTH = ICOARSELBX
         IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBXW)) THEN
            TZLCRSPD%TNEXT => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ELSE
            TZP2C_DATA%TSEND_1WAY_LBXW => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ENDIF
!
!         4.3  Increase number of elements in the list
!
        ICARD = ICARD + 1
      ELSE
         DEALLOCATE(TZLTMPCRSPD)
      ENDIF
      IPREV = ICOARSELBX
    ENDIF
  ENDDO
!
  IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBXW)) THEN
    TZP2C_DATA%TSEND_1WAY_LBXW%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZP2C_DATA%TSEND_1WAY_LBXW)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        5.   CONSTRUCTION OF TZP2C_DATA%TSEND_1WAY_LBXE
!              ------------------------------------------
!
  ICARD = 0
  IPREV = 0
!  GINTER = .FALSE.
!
!  WRITE(NIOUNIT,*) '****** EST ******'
!  WRITE(NIOUNIT,*) 'NDIMRIMLBX = ', TZCHILD_COMDATA%NDIMRIMLBX
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBX
!
!         5.1  Calculate the width of the coarse zone to be interpolated
!
    ICOARSELBX = LBFINE2COARSE(NDXRATIO_ALL(KMODEL), &
                 TZCHILD_COMDATA%NRIMLBX(JI))
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBX(JI)', &
!                     TZCHILD_COMDATA%NRIMLBX(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBX =',ICOARSELBX
!
!
    IF (ICOARSELBX /= IPREV) THEN
!
!         5.2  Set up the list of zones to be communicated
!
      CALL ALLOC(TZLTMPCRSPD)
      DO J = 1, NPROC
        IF (DIMX - TZCHILD_COMDATA%NRIMLBX(JI) <= TZFINEPB(J)%NXEND) THEN
          TZLB(J)%NUMBER = TZCOARSE(J)%NUMBER
          TZLB(J)%NXOR = MAX(NXOR_ALL(KMODEL) + TZCOARSE(J)%NXOR - 3, &
                             NXEND_ALL(KMODEL) - ICOARSELBX + 2)
          TZLB(J)%NYOR = NYOR_ALL(KMODEL) + TZCOARSE(J)%NYOR - 3
          TZLB(J)%NXEND = NXOR_ALL(KMODEL) + TZCOARSE(J)%NXEND + 1
          TZLB(J)%NYEND = NYOR_ALL(KMODEL) + TZCOARSE(J)%NYEND + 1
        ELSE
          TZLB(J) = ZONE_ll(0,0,0,0,0,0,0,0)
        ENDIF
!        CALL DISPLAY_ZONE(TZLB(J),1)
      ENDDO
!
      MESSTAG = IELBMSG + ICOARSELBX
      CALL PARENTS_PT_OF_VIEW(TZLB, TZPB(IP), TPPROCONF, MESSTAG, &
                              TZLTMPCRSPD, GINTER)
!
      IF (GINTER) THEN
         TZLTMPCRSPD%NWIDTH = ICOARSELBX
         IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBXE)) THEN
            TZLCRSPD%TNEXT => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ELSE
            TZP2C_DATA%TSEND_1WAY_LBXE => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ENDIF
!
!         5.3  Increase number of elements in the list
!
         ICARD = ICARD + 1
      ELSE
         DEALLOCATE(TZLTMPCRSPD)

      ENDIF
      IPREV = ICOARSELBX
    ENDIF
  ENDDO
!
  IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBXE)) THEN
    TZP2C_DATA%TSEND_1WAY_LBXE%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZP2C_DATA%TSEND_1WAY_LBXE)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        6.   CONSTRUCTION OF TZP2C_DATA%TSEND_1WAY_LBYS
!              ------------------------------------------
!
  ICARD = 0
  IPREV = 0
!  GINTER = .FALSE.
!
!  WRITE(NIOUNIT,*) '****** SUD ******'
!  WRITE(NIOUNIT,*) 'NDIMRIMLBX = ', TZCHILD_COMDATA%NDIMRIMLBX
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBY
!
!         6.1  Calculate the width of the coarse zone to be interpolated
!
    ICOARSELBY = LBFINE2COARSE(NDYRATIO_ALL(KMODEL), &
                               TZCHILD_COMDATA%NRIMLBY(JI))
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBY(JI)', &
!                     TZCHILD_COMDATA%NRIMLBY(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBY =',ICOARSELBY
!
    IF (ICOARSELBY /= IPREV) THEN
!
!         6.2  Set up the list of zones to be communicated
!
      CALL ALLOC(TZLTMPCRSPD)
      DO J = 1, NPROC
        IF (TZCHILD_COMDATA%NRIMLBY(JI)+JPHEXT >= TZFINEPB(J)%NYOR) THEN
          TZLB(J)%NUMBER = TZCOARSE(J)%NUMBER
          TZLB(J)%NXOR = NXOR_ALL(KMODEL) + TZCOARSE(J)%NXOR - 3
          TZLB(J)%NYOR = NYOR_ALL(KMODEL) + TZCOARSE(J)%NYOR - 3
          TZLB(J)%NXEND = NXOR_ALL(KMODEL) + TZCOARSE(J)%NXEND + 1
          TZLB(J)%NYEND = MIN(NYOR_ALL(KMODEL) + TZCOARSE(J)%NYEND + 1, &
                              NYOR_ALL(KMODEL) + ICOARSELBY - 2)
        ELSE
          TZLB(J) = ZONE_ll(0,0,0,0,0,0,0,0)
        ENDIF
!        CALL DISPLAY_ZONE(TZLB(J),1)
      ENDDO
!
      MESSTAG = ISLBMSG + ICOARSELBY
      CALL PARENTS_PT_OF_VIEW(TZLB, TZPB(IP), TPPROCONF, MESSTAG, &
                              TZLTMPCRSPD, GINTER)
!
      IF (GINTER) THEN
         TZLTMPCRSPD%NWIDTH = ICOARSELBY
         IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBYS)) THEN
            TZLCRSPD%TNEXT => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ELSE
            TZP2C_DATA%TSEND_1WAY_LBYS => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ENDIF
!
!         6.3  Increase number of elements in the list
!
         ICARD = ICARD + 1
      ELSE
         DEALLOCATE(TZLTMPCRSPD)
      ENDIF
      IPREV = ICOARSELBY
    ENDIF
  ENDDO
!
  IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBYS)) THEN
    TZP2C_DATA%TSEND_1WAY_LBYS%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZP2C_DATA%TSEND_1WAY_LBYS)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        7.   CONSTRUCTION OF TZP2C_DATA%TSEND_1WAY_LBYN
!              ------------------------------------------
!
  ICARD = 0
  IPREV = 0
!  GINTER = .FALSE.
!
!  WRITE(NIOUNIT,*) '****** NORD ******'
!  WRITE(NIOUNIT,*) 'NDIMRIMLBX = ', TZCHILD_COMDATA%NDIMRIMLBX
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBY
!
!         7.1  Calculate the width of the coarse zone to be interpolated
!
    ICOARSELBY = LBFINE2COARSE(NDYRATIO_ALL(KMODEL), &
                               TZCHILD_COMDATA%NRIMLBY(JI))
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBY(JI)', &
!                     TZCHILD_COMDATA%NRIMLBY(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBY =',ICOARSELBY
!
    IF (ICOARSELBY /= IPREV) THEN
!
!         7.2  Set up the list of zones to be communicated
!
      CALL ALLOC(TZLTMPCRSPD)
      DO J = 1, NPROC
        IF (DIMY - TZCHILD_COMDATA%NRIMLBY(JI) <= TZFINEPB(J)%NYEND) THEN
          TZLB(J)%NUMBER = TZCOARSE(J)%NUMBER
          TZLB(J)%NXOR = NXOR_ALL(KMODEL) + TZCOARSE(J)%NXOR - 3
          TZLB(J)%NYOR = MAX(NYOR_ALL(KMODEL) + TZCOARSE(J)%NYOR - 3, &
                             NYEND_ALL(KMODEL) - ICOARSELBY + 2)
          TZLB(J)%NXEND = NXOR_ALL(KMODEL) + TZCOARSE(J)%NXEND + 1
          TZLB(J)%NYEND = NYOR_ALL(KMODEL) + TZCOARSE(J)%NYEND + 1
        ELSE
          TZLB(J) = ZONE_ll(0,0,0,0,0,0,0,0)
        ENDIF
!        CALL DISPLAY_ZONE(TZLB(J),1)
      ENDDO
!
      MESSTAG = INLBMSG + ICOARSELBY
      CALL PARENTS_PT_OF_VIEW(TZLB, TZPB(IP), TPPROCONF, MESSTAG, &
                              TZLTMPCRSPD, GINTER)
!
      IF (GINTER) THEN
         TZLTMPCRSPD%NWIDTH = ICOARSELBY
         IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBYN)) THEN
            TZLCRSPD%TNEXT => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ELSE
            TZP2C_DATA%TSEND_1WAY_LBYN => TZLTMPCRSPD
            TZLCRSPD => TZLTMPCRSPD
         ENDIF
!
!         7.3  Increase number of elements in the list
!
        ICARD = ICARD + 1
      ENDIF
      IPREV = ICOARSELBY
    ENDIF
  ENDDO
!
  IF (ASSOCIATED(TZP2C_DATA%TSEND_1WAY_LBYN)) THEN
    TZP2C_DATA%TSEND_1WAY_LBYN%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZP2C_DATA%TSEND_1WAY_LBYN)
  ENDIF
  NULLIFY(TZLCRSPD)
!
!-------------------------------------------------------------------------------
!
!*        8.   CONSTRUCTION OF THE PART OF TPCOMDATA CONCERNING DATA
!              FOR LB COMMUNICATIONS FOR THE CHILD KMODEL (RECEPTION)
!              ------------------------------------------------------
!
!         8.1 West side
!
  IPREV = 0
  ICARD = 0
!
!  WRITE(NIOUNIT,*) '****** OUEST RECEPTION ******'
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBX
!
!         8.1.1 Calculate the width of the coarse zone to be interpolated
!
      ICOARSELBX = LBFINE2COARSE(NDXRATIO_ALL(KMODEL), &
                                 TZCHILD_COMDATA%NRIMLBX(JI))
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBX(JI)', &
!                     TZCHILD_COMDATA%NRIMLBX(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBX =',ICOARSELBX
!
      IF (ICOARSELBX /= IPREV) THEN
!
!         8.1.2 Set up the list of zones to be communicated
!
         IF (TZCHILD_COMDATA%NRIMLBX(JI)+JPHEXT >= TZFINEPB(IP)%NXOR) THEN
            IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBXW)) THEN
               CALL ALLOC(TZLCRSPD%TNEXT)
               TZLCRSPD => TZLCRSPD%TNEXT
            ELSE
               CALL ALLOC(TZCHILD_COMDATA%TRECV_LBXW)
               TZLCRSPD => TZCHILD_COMDATA%TRECV_LBXW
            ENDIF
!
            TZCHILDLB%NUMBER = TZCOARSE(IP)%NUMBER
            TZCHILDLB%NXOR = NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXOR - 3
            TZCHILDLB%NYOR = NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYOR - 3
            TZCHILDLB%NXEND = MIN(NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXEND + 1, &
                                  NXOR_ALL(KMODEL) + ICOARSELBX - 2)
            TZCHILDLB%NYEND = NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYEND + 1
            TZLCRSPD%NWIDTH = ICOARSELBX
!
            MESSTAG = IWLBMSG + ICOARSELBX
            CALL CHILDS_PT_OF_VIEW(TZCHILDLB, TZPB, TZCHILDLB, MESSTAG, &
                                   TZLCRSPD%TCRSPD)
            ICARD = ICARD + 1
!            CALL DISPLAY_ZONE(TZCHILDLB,1)
         ENDIF
      ENDIF
      IPREV = ICOARSELBX
   ENDDO
!
   IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBXW)) THEN
      TZCHILD_COMDATA%TRECV_LBXW%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZCHILD_COMDATA%TRECV_LBXW)
   ENDIF
!
!         8.2  East side
!
  IPREV = 0
  ICARD = 0
!
!  WRITE(NIOUNIT,*) '****** EST RECEPTION ******'
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBX
!
!         8.2.1 Calculate the width of the coarse zone to be interpolated
!
      ICOARSELBX = LBFINE2COARSE(NDXRATIO_ALL(KMODEL), &
                                 TZCHILD_COMDATA%NRIMLBX(JI))
!
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBX(JI)', &
!                     TZCHILD_COMDATA%NRIMLBX(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBX =',ICOARSELBX
!
      IF (ICOARSELBX /= IPREV) THEN
!
!         8.2.2  Set up the list of zones to be communicated
!
         IF (DIMX - TZCHILD_COMDATA%NRIMLBX(JI) <= TZFINEPB(IP)%NXEND) THEN
            IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBXE)) THEN
               CALL ALLOC(TZLCRSPD%TNEXT)
               TZLCRSPD => TZLCRSPD%TNEXT
            ELSE
               CALL ALLOC(TZCHILD_COMDATA%TRECV_LBXE)
               TZLCRSPD => TZCHILD_COMDATA%TRECV_LBXE
            ENDIF
!
            TZCHILDLB%NUMBER = TZCOARSE(IP)%NUMBER
            TZCHILDLB%NXOR = MAX(NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXOR - 3, &
                                 NXEND_ALL(KMODEL) - ICOARSELBX + 2)
            TZCHILDLB%NYOR = NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYOR - 3
            TZCHILDLB%NXEND = NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXEND + 1
            TZCHILDLB%NYEND = NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYEND + 1
            TZLCRSPD%NWIDTH = ICOARSELBX
            MESSTAG = IELBMSG + ICOARSELBX
            CALL CHILDS_PT_OF_VIEW(TZCHILDLB, TZPB, TZCHILDLB, MESSTAG, &
                                   TZLCRSPD%TCRSPD)
            ICARD = ICARD + 1
!            CALL DISPLAY_ZONE(TZCHILDLB,1)
         ENDIF
      ENDIF
      IPREV = ICOARSELBX
  ENDDO
!
  IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBXE)) THEN
      TZCHILD_COMDATA%TRECV_LBXE%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZCHILD_COMDATA%TRECV_LBXE)
  ENDIF
!
!         8.3  South side
!
  IPREV = 0
  ICARD = 0
!
!  WRITE(NIOUNIT,*) '****** SUD RECEPTION ******'
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBY
!
!         8.3.1 Calculate the width of the coarse zone to be interpolated
!
      ICOARSELBY = LBFINE2COARSE(NDYRATIO_ALL(KMODEL), &
                                 TZCHILD_COMDATA%NRIMLBY(JI))
!
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBY(JI)', &
!                     TZCHILD_COMDATA%NRIMLBY(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBY =',ICOARSELBY
!
      IF (ICOARSELBY /= IPREV) THEN
!
!         8.3.2 Set up the list of zones to be communicated
!
         IF (TZCHILD_COMDATA%NRIMLBY(JI)+JPHEXT >= TZFINEPB(IP)%NYOR) THEN
            IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBYS)) THEN
               CALL ALLOC(TZLCRSPD%TNEXT)
               TZLCRSPD => TZLCRSPD%TNEXT
            ELSE
               CALL ALLOC(TZCHILD_COMDATA%TRECV_LBYS)
               TZLCRSPD => TZCHILD_COMDATA%TRECV_LBYS
            ENDIF
!
            TZCHILDLB%NUMBER = TZCOARSE(IP)%NUMBER
            TZCHILDLB%NXOR = NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXOR - 3
            TZCHILDLB%NYOR = NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYOR - 3
            TZCHILDLB%NXEND = NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXEND + 1
            TZCHILDLB%NYEND = MIN(NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYEND + 1, &
                                  NYOR_ALL(KMODEL) + ICOARSELBY - 2)
            TZLCRSPD%NWIDTH = ICOARSELBY
            MESSTAG = ICOARSELBY + ISLBMSG
            CALL CHILDS_PT_OF_VIEW(TZCHILDLB, TZPB, TZCHILDLB, MESSTAG, &
                                   TZLCRSPD%TCRSPD)
            ICARD = ICARD + 1
!            CALL DISPLAY_ZONE(TZCHILDLB,1)
         ENDIF
      ENDIF
      IPREV = ICOARSELBY
  ENDDO
!
  IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBYS)) THEN
      TZCHILD_COMDATA%TRECV_LBYS%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZCHILD_COMDATA%TRECV_LBYS)
  ENDIF
!
!         8.4  North side
!
  IPREV = 0
  ICARD = 0
!
!  WRITE(NIOUNIT,*) '****** NORD RECEPTION ******'
!
  DO JI=1, TZCHILD_COMDATA%NDIMRIMLBY
!
!         8.4.1 Calculate the width of the coarse zone to be interpolated
!
      ICOARSELBY = LBFINE2COARSE(NDYRATIO_ALL(KMODEL), &
                                 TZCHILD_COMDATA%NRIMLBY(JI))
!
!    WRITE(NIOUNIT,*) JI, 'TZCHILD_COMDATA%NRIMLBY(JI)', &
!                     TZCHILD_COMDATA%NRIMLBY(JI)
!    WRITE(NIOUNIT,*) JI, 'ICOARSELBY =',ICOARSELBY
!
      IF (ICOARSELBY /= IPREV) THEN
!
!         8.4.2 Set up the list of zones to be communicated
!
         IF (DIMY -  TZCHILD_COMDATA%NRIMLBY(JI) <= TZFINEPB(IP)%NYEND) THEN
            IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBYN)) THEN
               CALL ALLOC(TZLCRSPD%TNEXT)
               TZLCRSPD => TZLCRSPD%TNEXT
            ELSE
               CALL ALLOC(TZCHILD_COMDATA%TRECV_LBYN)
               TZLCRSPD => TZCHILD_COMDATA%TRECV_LBYN
            ENDIF
!
            TZCHILDLB%NUMBER = TZCOARSE(IP)%NUMBER
            TZCHILDLB%NXOR = NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXOR - 3
            TZCHILDLB%NYOR = MAX(NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYOR - 3, &
                                 NYEND_ALL(KMODEL) - ICOARSELBY + 2)
            TZCHILDLB%NXEND = NXOR_ALL(KMODEL) + TZCOARSE(IP)%NXEND + 1
            TZCHILDLB%NYEND = NYOR_ALL(KMODEL) + TZCOARSE(IP)%NYEND + 1
            TZLCRSPD%NWIDTH = ICOARSELBY
            MESSTAG = INLBMSG + ICOARSELBY
            CALL CHILDS_PT_OF_VIEW(TZCHILDLB, TZPB, TZCHILDLB, MESSTAG, &
                                   TZLCRSPD%TCRSPD)
            ICARD = ICARD + 1
!            CALL DISPLAY_ZONE(TZCHILDLB,1)
         ENDIF
      ENDIF
      IPREV = ICOARSELBY
  ENDDO
!
  IF (ASSOCIATED(TZCHILD_COMDATA%TRECV_LBYN)) THEN
    TZCHILD_COMDATA%TRECV_LBYN%NCARD = ICARD
!    CALL DISPLAY_LCRSPD(TZCHILD_COMDATA%TRECV_LBYN)
  ENDIF
!
  DEALLOCATE(TZLB)
  DEALLOCATE(TZFINE)
  DEALLOCATE(TZFINEPB)
  DEALLOCATE(TZCOARSE)
  DEALLOCATE(TZP)
  DEALLOCATE(TZPB)
!
!-------------------------------------------------------------------------------
!
!*        9.   INITIALIZATION FOR THIS CHILDREN (MODEL KMODEL)
!              -----------------------------------------------
!
  DO JMODEL = 1, JPMODELMAX
!
    IF( NDAD(JMODEL) .EQ. TZCHILD_PROCONF%NUMBER ) THEN
      CALL INIT_LB_MODEL(TZCHILD_PROCONF, TZCHILD_COMDATA, JMODEL )
    ENDIF
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INIT_LB_MODEL
!
!    ########################
      SUBROUTINE INIT_LB_ll()
!    ########################
!
!!****  *INIT_LB_ll* - routine to construct the data structures
!                      for LB FIELDS of all descendants models of
!                      current model ;
!                      called with model number 1, this routine
!                      constructs the LB-structures for all models
! 
!!    Purpose
!!    -------
!     this routine fills all the data structures (splitting data structures,
!     communications data structures) for LB FIELDS
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
!       NDAD
!
!     Module MODD_PARAMETERS_ll
!       JPMODELMAX
!
!     Module MODD_VAR_ll
!       TCRRT_PROCONF, TCRRT_COMDATA
!
!!    Reference
!!    ---------
!     Notes sur les routines multi-modeles de la surcouche, R. Guivarch
! 
!!    Author
!!    ------
!     N. Gicquel      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : NDAD
  USE MODD_PARAMETERS_ll, ONLY : JPMODELMAX
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, TCRRT_COMDATA
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
!
!*       0.2   declarations of local variables
!
  INTEGER :: JMODEL, IINFO
!
!-------------------------------------------------------------------------------
!
  DO JMODEL = 1, JPMODELMAX
!
    IF( NDAD(JMODEL) .EQ. TCRRT_PROCONF%NUMBER ) THEN
      CALL INIT_LB_MODEL(TCRRT_PROCONF, TCRRT_COMDATA, JMODEL )
    ENDIF
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INIT_LB_ll
!
!
!
      SUBROUTINE SET_LB_FIELD_ll( HLBTYPE, PFIELD, PLBXFIELD, PLBYFIELD, IIB, IJB,&
                     IIE, IJE, SHIFTWEST, SHIFTEAST, SHIFTSOUTH, SHIFTNORTH )
!     #######################################################################
!
!!****  *SET_LB_FIELD_ll * - subroutine to copy the values associated with the
!!                           Lateral Boundaries to the corresoponding LB field
!!
!!    AUTHOR
!!    ------
!!
!!       M. Moge     * LA, CNRS *
!!
!!      Original     28/11/14
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
  USE MODD_CONF
!  USE MODD_DIM_n
  USE MODD_DYN_n
  USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
!  USE MODE_ll
  USE MODE_IO_ll
  USE MODE_MPPDB
  USE MODE_DISTRIB_LB
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*),INTENT(IN) :: HLBTYPE ! LB type : 'LB','LBU'
  REAL, DIMENSION(:,:,:), INTENT(IN)  :: PFIELD      ! field on the whole domain (or subdomain)
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLBXFIELD    ! LB field - X direction
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLBYFIELD    ! LB field - Y direction
  !beginning and end of the local physical subdomain
  INTEGER, INTENT(IN)   :: IIB            ! indice I Beginning in x direction
  INTEGER, INTENT(IN)   :: IJB            ! indice J Beginning in y direction
  INTEGER, INTENT(IN)   :: IIE            ! indice I End       in x direction
  INTEGER, INTENT(IN)   :: IJE            ! indice J End       in y direction
  INTEGER, INTENT(IN)   :: SHIFTWEST, SHIFTEAST, SHIFTSOUTH, SHIFTNORTH ! shifting applied to the indices copied from PFIELD in each direction
                                                                        ! it is used for LBXUM et LBXVM
                                                                        ! I do not know why...
  !
  ! LOCAL VARIABLES
  CHARACTER(4) :: YLBTYPEX ! LB type : 'LBX','LBXU'
  CHARACTER(4) :: YLBTYPEY ! LB type : 'LBY','LBYV'
  ! local indices for the intersection of the local subdomain and the LB zone
  INTEGER             :: IIB_LOCLB           ! indice I Beginning in x direction
  INTEGER             :: IJB_LOCLB           ! indice J Beginning in y direction
  INTEGER             :: IIE_LOCLB           ! indice I End       in x direction
  INTEGER             :: IJE_LOCLB           ! indice J End       in y direction
  ! global indices for the intersection of the local subdomain and the LB zone
  INTEGER             :: IIB_GLBLB           ! indice I Beginning in x direction
  INTEGER             :: IJB_GLBLB           ! indice J Beginning in y direction
  INTEGER             :: IIE_GLBLB           ! indice I End       in x direction
  INTEGER             :: IJE_GLBLB           ! indice J End       in y direction
  INTEGER             :: LOCLBSIZEE, LOCLBSIZEW, LOCLBSIZEN, LOCLBSIZES ! size of the local portion of the LB zone in each direction (East, West, North, South)
  INTEGER             :: GLBLBBEGIN,GLBLBEND
  !
  ! SET LB TYPE
  IF ( HLBTYPE == 'LB' ) THEN
    YLBTYPEX = 'LBX'
    YLBTYPEY = 'LBY'
  ELSE IF ( HLBTYPE == 'LBU' ) THEN
    YLBTYPEX = 'LBXU'
    YLBTYPEY = 'LBYV'
  ELSE
    WRITE(*,*) "ERROR: from SET_LB_FIELD_ll, UNKNOWN LB TYPE", HLBTYPE
    CALL ABORT
  ENDIF
!
! get the local indices of the West-East LB arrays for the local subdomain
  CALL GET_DISTRIB_LB(YLBTYPEX,ISP,'LOC','WRITE',NRIMX,IIB_LOCLB,IIE_LOCLB,IJB_LOCLB,IJE_LOCLB)
! and the corresponding indices for the LB global arrays
  CALL GET_DISTRIB_LB(YLBTYPEX,ISP,'FM','WRITE',NRIMX,IIB_GLBLB,IIE_GLBLB,IJB_GLBLB,IJE_GLBLB)
  IF ( IIE_LOCLB-IIB_LOCLB /= IIE_GLBLB-IIB_GLBLB ) THEN
    WRITE(*,*) "ERROR: from SET_LB_FIELD_ll, West-East IIE_LOCLB-IIB_LOCLB =",&
        IIE_LOCLB-IIB_LOCLB, " /= IIE_GLBLB-IIB_GLBLB =", IIE_GLBLB-IIB_GLBLB
    CALL ABORT
  ENDIF
  LOCLBSIZEW = 0
  LOCLBSIZEE = 0
  IF ( IIB_LOCLB /= 0 ) THEN  ! if the LB zone of the local subdomain is non-empty
    ! WARNING : The size of the local portion of the LB zone can be less than NRIMX
    ! Example : if the size of the subdomain is 4 and NRIMX=6, the LB zone will be divided between 2 processes
    !           and LOCLBSIZEW will be 5 on the first process, and 2 on the second process
    IF ( IIB_GLBLB <= NRIMX+JPHEXT .AND. IIE_GLBLB >= NRIMX+JPHEXT+1 ) THEN ! the local west and east LB zones are both non empty
      LOCLBSIZEW = NRIMX+JPHEXT-IIB_GLBLB
      PLBXFIELD(IIB_LOCLB:IIB_LOCLB+LOCLBSIZEW,:,:)  = PFIELD(IIB_GLBLB+SHIFTWEST:IIB_GLBLB+SHIFTWEST+LOCLBSIZEW,:,:)
      PLBXFIELD(IIE_LOCLB-LOCLBSIZEW:IIE_LOCLB,:,:)  = PFIELD(IIE+JPHEXT-LOCLBSIZEW+SHIFTEAST:IIE+JPHEXT+SHIFTEAST,:,:)
    ELSE IF ( IIB_GLBLB <= NRIMX+JPHEXT ) THEN  ! the local west LB zone only is non empty
      LOCLBSIZEW = NRIMX+JPHEXT-IIB_GLBLB
      PLBXFIELD(IIB_LOCLB:IIE_LOCLB,:,:)  = PFIELD(IIB_GLBLB+SHIFTWEST:IIE_GLBLB+SHIFTWEST,:,:)
    ELSE IF ( IIB_GLBLB >= NRIMX+JPHEXT+1 ) THEN  ! the local east LB zone only is non empty
!      LOCLBSIZEE = IIE_LOCLB-IIB_LOCLB
      GLBLBBEGIN = IIE+JPHEXT-(2*NRIMX+2*JPHEXT-IIB_GLBLB)+SHIFTEAST
      GLBLBEND = IIE+JPHEXT-(2*NRIMX+2*JPHEXT-IIE_GLBLB)+SHIFTEAST
      PLBXFIELD(IIB_LOCLB:IIE_LOCLB,:,:)  = PFIELD(GLBLBBEGIN:GLBLBEND,:,:)
!      PLBXFIELD(NRIMX+1+IIB_LOCLB:NRIMX+1+IIE_LOCLB,:,:)  = PFIELD(GLBLBBEGIN:GLBLBEND,:,:)
    ELSE
      WRITE(*,*) "ERROR: from SET_LB_FIELD_ll, This type of partition is not allowed !"
      CALL ABORT
    ENDIF
  ENDIF !( IIB_LOCLB /= 0 )
!
!*       5.9.1.8  Y-direction variables
!
  IF( .NOT. L2D ) THEN
    LOCLBSIZES = 0
    LOCLBSIZEN = 0
  ! get the local indices of the South-North LB arrays for the local subdomain
    CALL GET_DISTRIB_LB(YLBTYPEY,ISP,'LOC','WRITE',NRIMY,IIB_LOCLB,IIE_LOCLB,IJB_LOCLB,IJE_LOCLB)
  ! and the corresponding indices for the LB global arrays
    CALL GET_DISTRIB_LB(YLBTYPEY,ISP,'FM','WRITE',NRIMY,IIB_GLBLB,IIE_GLBLB,IJB_GLBLB,IJE_GLBLB)
    IF ( IJE_LOCLB-IJB_LOCLB /= IJE_GLBLB-IJB_GLBLB ) THEN
      WRITE(*,*) "ERROR: from SET_LB_FIELD_ll, South-North IJE_LOCLB-IJB_LOCLB =",&
           IJE_LOCLB-IJB_LOCLB, " /= IJE_GLBLB-IJB_GLBLB =", IJE_GLBLB-IJB_GLBLB
      CALL ABORT
    ENDIF
    IF ( IJB_LOCLB /= 0 ) THEN  ! if the LB zone of the local subdomain is non-empty
      IF ( IJB_GLBLB <= NRIMY+JPHEXT .AND. IJE_GLBLB >= NRIMY+JPHEXT+1 ) THEN ! the local south and north LB zones are non empty
        LOCLBSIZES = NRIMY+JPHEXT-IJB_GLBLB
        PLBYFIELD(:,IJB_LOCLB:IJB_LOCLB+LOCLBSIZES,:)  = PFIELD(:,IJB_GLBLB+SHIFTSOUTH:IJB_GLBLB+LOCLBSIZES+SHIFTSOUTH,:)
        PLBYFIELD(:,IJE_LOCLB-LOCLBSIZES:IJE_LOCLB,:)  = PFIELD(:,IJE+JPHEXT-LOCLBSIZES+SHIFTNORTH:IJE+JPHEXT+SHIFTNORTH,:)
      ELSE IF ( IJB_GLBLB <= NRIMY+JPHEXT ) THEN  ! the local south LB zone only is non empty
        LOCLBSIZES = NRIMY+JPHEXT-IJB_GLBLB
        PLBYFIELD(:,IJB_LOCLB:IJE_LOCLB,:)  = PFIELD(:,IJB_GLBLB+SHIFTSOUTH:IJE_GLBLB+SHIFTSOUTH,:)
      ELSE IF ( IJB_GLBLB >= NRIMY+JPHEXT+1 ) THEN  ! the local north LB zone only is non empty
        GLBLBBEGIN = IJE+JPHEXT-(2*NRIMY+2*JPHEXT-IJB_GLBLB)+SHIFTNORTH
        GLBLBEND = IJE+JPHEXT-(2*NRIMY+2*JPHEXT-IJE_GLBLB)+SHIFTNORTH
        PLBYFIELD(:,IJB_LOCLB:IJE_LOCLB,:)  = PFIELD(:,GLBLBBEGIN:GLBLBEND,:)
!        PLBYFIELD(:,NRIMY+1+IJB_LOCLB:NRIMY+1+IJE_LOCLB,:)  = PFIELD(:,GLBLBBEGIN:GLBLBEND,:)
      ELSE
        WRITE(*,*) "ERROR: from SET_LB_FIELD_ll, This type of partition is not allowed !"
        CALL ABORT
      ENDIF

    ENDIF !( IJB_LOCLB /= 0 )
  ENDIF !( .NOT. L2D )
!
      END SUBROUTINE SET_LB_FIELD_ll
!
!
!
      FUNCTION GET_LOCAL_LB_SIZE_X_ll( KRIMX  ) RESULT( LBSIZEX )
!     #######################################################################
!
!!****  *GET_LOCAL_LB_SIZE_X_ll * - get the local LB size in X direction,
!!       i.e. the size of the array containing the local portion of the LB zone
!!
!!    AUTHOR
!!    ------
!!
!!       M. Moge     * LA, CNRS *
!!
!!      Original     01/12/14
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
  USE MODE_ll
  !
  IMPLICIT NONE
  !

  INTEGER, INTENT(IN) :: KRIMX               ! global LB size in X direction (input)
  INTEGER             :: LBSIZEX             ! local LB size in X direction (output)
                                             ! Size of the array containing the local portion of the LB zone
  LBSIZEX = 0
  IF( LWEST_ll() ) THEN
    LBSIZEX = LBSIZEX + KRIMX+1
  ENDIF
  IF( LEAST_ll() ) THEN
    LBSIZEX = LBSIZEX + KRIMX+1
  ENDIF
!
      END FUNCTION GET_LOCAL_LB_SIZE_X_ll
!
!
!
      FUNCTION GET_LOCAL_LB_SIZE_Y_ll( KRIMY  ) RESULT( LBSIZEY )
!     #######################################################################
!
!!****  *GET_LOCAL_LB_SIZE_Y_ll * - get the local LB size in Y direction,
!!       i.e. the size of the array containing the local portion of the LB zone
!!
!!    AUTHOR
!!    ------
!!
!!       M. Moge     * LA, CNRS *
!!
!!      Original     01/12/14
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
  USE MODE_ll
  !
  IMPLICIT NONE
  !

  INTEGER, INTENT(IN) :: KRIMY               ! global LB size in Y direction (input)
  INTEGER             :: LBSIZEY             ! local LB size in Y direction (output)
                                             ! Size of the array containing the local portion of the LB zone
  LBSIZEY = 0
  IF( LSOUTH_ll() ) THEN
    LBSIZEY = LBSIZEY + KRIMY+1
  ENDIF
  IF( LNORTH_ll() ) THEN
    LBSIZEY = LBSIZEY + KRIMY+1
  ENDIF
!
      END FUNCTION GET_LOCAL_LB_SIZE_Y_ll
!
END MODULE MODE_LB_ll
