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

!     #################
      MODULE MODE_LS_ll
!     #################
! 
!!    Purpose
!!    -------
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!     SUBROUTINES : SET_LSFIELD_1WAY_ll
!                       (SET_LS2DFIELD_1WAY_ll, SET_LS2DFIELD_1WAY_ll),
!                   UNSET_LSFIELD_1WAY_ll,
!                   SET_LSFIELD_2WAY_ll
!                       (SET_LS2DFIELD_2WAY_ll, SET_LS3DFIELD_2WAY_ll),
!                   UNSET_LSFIELD_2WAY_ll,
!                   LS_FORCING_ll, LS_FEEDBACK_ll
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
!
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       type LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll,
!            PARENT2CHILD_DATA_ll, LPROCONF_ll, PROC_COM_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NGRID_COM - mpi communicator
!       NBUFFERSIZE_3D - buffer size
!
!------------------------------------------------------------------------------
!
  USE MODD_STRUCTURE_ll
!
  CONTAINS
!
!     ###########################################################
      SUBROUTINE SET_LS2DFIELD_1WAY_ll(P2DFIELD, PTFIELD, KMODEL)
!     ###########################################################
!
!!****  *SET_LS2DFIELD_1WAY_ll*- routine to fill the parallel data lists
!                                of forcing
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
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll
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
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 11 fev. 2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LCRSPD_ll, LPROC_COM_DATA_ll, &
                                LPARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  USE MODE_ARGSLIST_ll, ONLY : ADD2DFIELD_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
  INTEGER, INTENT(IN) :: KMODEL
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
    WRITE(*,*) 'Problem in SET_LS2DFIELD_1WAY_ll'
    WRITE(*,*) 'The current model has no child'
    RETURN
  ENDIF
!
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
    WRITE(*,*) 'Error SET_LS2DFIELD_1WAY_ll : ', KMODEL, &
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
    WRITE(*,*) 'Error SET_LS2DFIELD_1WAY_ll : ', KMODEL, &
               ' is not a child of the current model'
    STOP
  ENDIF
!
  TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LS
  TZCHILD => TZLCOMDATA%TELT%TRECV_1WAY_LS
!
!-------------------------------------------------------------------------------
!
!*       3.    Put P2DFIELD and PTFIELD in the list of fields
!
   CALL ADD2DFIELD_ll(TZCHILD%TLIST, PTFIELD)
!
   CALL ADD2DFIELD_ll(TZPAR%TLIST, P2DFIELD)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LS2DFIELD_1WAY_ll
!
!     ###########################################################
      SUBROUTINE SET_LS3DFIELD_1WAY_ll(P3DFIELD, PTFIELD, KMODEL)
!     ###########################################################
!
!!****  *SET_LS3DFIELD_1WAY_ll*- routine to fill the parallel data lists
!                                of forcing
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
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll
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
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 11 fev. 2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LCRSPD_ll, LPROC_COM_DATA_ll, &
                                LPARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  USE MODE_ARGSLIST_ll, ONLY : ADD3DFIELD_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
  INTEGER, INTENT(IN) :: KMODEL
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
    WRITE(*,*) 'Problem in SET_LS3DFIELD_1WAY_ll'
    WRITE(*,*) 'The current model has no child'
    RETURN
  ENDIF
!
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
    WRITE(*,*) 'Error SET_LS3DFIELD_1WAY_ll : ', KMODEL, &
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
    WRITE(*,*) 'Error SET_LS3DFIELD_1WAY_ll : ', KMODEL, &
               ' is not a child of the current model'
    STOP
  ENDIF
!
  TZPAR => TZP2CDATA%TELT%TSEND_1WAY_LS
  TZCHILD => TZLCOMDATA%TELT%TRECV_1WAY_LS
!
!-------------------------------------------------------------------------------
!
!*       3.    Put P3DFIELD and PTFIELD in the list of fields
!
   CALL ADD3DFIELD_ll(TZCHILD%TLIST, PTFIELD)
!
   CALL ADD3DFIELD_ll(TZPAR%TLIST, P3DFIELD)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LS3DFIELD_1WAY_ll
!
!     ##################################
      SUBROUTINE UNSET_LSFIELD_1WAY_ll()
!     ##################################
!
!!****  *UNSET_LSFIELD_1WAY_ll*- routine to clean the parallel data lists
!                                of forcing
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
!       type LPARENT2CHILD_DATA_ll, PARENT2CHILD_DATA_ll
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
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 11 fev. 2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
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
  IF (TCRRT_COMDATA%NUMBER > 1) THEN
!
    CALL CLEANLIST_LCRSPD(TCRRT_COMDATA%TRECV_1WAY_LS)
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
! Clean the lists corresponding to the parent's send for the LS
!
    CALL CLEANLIST_LCRSPD(TZP2CDATA%TSEND_1WAY_LS)
!
  ELSE
!
    WRITE(*,*) 'Problem in UNSET_LSFIELD_1WAY_ll'
    WRITE(*,*) 'The current model is 1'
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UNSET_LSFIELD_1WAY_ll
!
!     ###################################################
      SUBROUTINE SET_LS2DFIELD_2WAY_ll(P2DFIELD, PTFIELD)
!     ###################################################
!
!!****  *SET_LS2DFIELD_2WAY_ll*- routine to fill the parallel data lists
!                                of feed back
!!    Purpose
!!    -------
!     this routine adds the 2D field P2DFIELD (resp. PTFIELD)
!     to the feed back exchange parent (resp. child) list
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!     Module MODE_ARGSLIST_ll
!       ADD2DFIELD_ll
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll
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
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 11 fev. 2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LCRSPD_ll, LPROC_COM_DATA_ll, &
                                LPARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  USE MODE_ARGSLIST_ll, ONLY : ADD2DFIELD_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
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
!*       1.    Point to the appropriate communication structures
!
!*       1.1   Point to the child communication data structure
!
  TZP2CDATA => TCRRT_COMDATA%TPARENT%TP2C_DATA
  DO WHILE (TZP2CDATA%TELT%NUMBER /= TCRRT_COMDATA%NUMBER)
    TZP2CDATA => TZP2CDATA%TNEXT
  ENDDO
!
  TZPAR => TZP2CDATA%TELT%TRECV_2WAY_LS
  TZCHILD => TCRRT_COMDATA%TSEND_2WAY_LS
!
!-------------------------------------------------------------------------------
!
!*       2.   Put P2DFIELD and PTFIELD in the list of fields
!
   CALL ADD2DFIELD_ll(TZCHILD%TLIST, PTFIELD)
!
   CALL ADD2DFIELD_ll(TZPAR%TLIST, P2DFIELD)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LS2DFIELD_2WAY_ll
!
!     ###################################################
      SUBROUTINE SET_LS3DFIELD_2WAY_ll(P3DFIELD, PTFIELD)
!     ###################################################
!
!!****  *SET_LS3DFIELD_2WAY_ll*- routine to fill the parallel data lists
!                                of feed back
!!    Purpose
!!    -------
!     this routine adds the 3D field P3DFIELD (resp. PTFIELD)
!     to the feed back exchange parent (resp. child) list
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!     Module MODE_ARGSLIST_ll
!       ADD3DFIELD_ll
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       type LCRSPD_ll, LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll
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
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 11 fev. 2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LCRSPD_ll, LPROC_COM_DATA_ll, &
                                LPARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  USE MODE_ARGSLIST_ll, ONLY : ADD3DFIELD_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
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
!*       1.    Point to the appropriate communication structures
!
!*       1.1   Point to the child communication data structure
!
  TZP2CDATA => TCRRT_COMDATA%TPARENT%TP2C_DATA
  DO WHILE (TZP2CDATA%TELT%NUMBER /= TCRRT_COMDATA%NUMBER)
    TZP2CDATA => TZP2CDATA%TNEXT
  ENDDO
!
  TZPAR => TZP2CDATA%TELT%TRECV_2WAY_LS
  TZCHILD => TCRRT_COMDATA%TSEND_2WAY_LS
!
!-------------------------------------------------------------------------------
!
!*       2.   Put P3DFIELD and PTFIELD in the list of fields
!
   CALL ADD3DFIELD_ll(TZCHILD%TLIST, PTFIELD)
!
   CALL ADD3DFIELD_ll(TZPAR%TLIST, P3DFIELD)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_LS3DFIELD_2WAY_ll
!
!     ########################################
      SUBROUTINE UNSET_LSFIELD_2WAY_ll(KMODEL)
!     ########################################
!
!!****  *UNSET_LSFIELD_2WAY_ll*- routine to clean the parallel data lists
!                                of forcing
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
!       type LPARENT2CHILD_DATA_ll, PARENT2CHILD_DATA_ll
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
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 11 fev. 2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LPARENT2CHILD_DATA_ll, PARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
  USE MODE_CONSTRUCT_ll, ONLY : CLEANLIST_LCRSPD
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN) :: KMODEL
!
!*       0.2   declarations of local variables
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZP2CDATA
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZLCOMDATA
!
!-------------------------------------------------------------------------------
!
!*       1.1   Point to the child communication data structure
!
  TZLCOMDATA => TCRRT_COMDATA%TCHILDREN
  DO WHILE (TZLCOMDATA%TELT%NUMBER /= KMODEL)
    TZLCOMDATA => TZLCOMDATA%TNEXT
  ENDDO
  IF (.NOT.ASSOCIATED(TZLCOMDATA)) THEN
    WRITE(*,*) 'Error UNSET_LSFIELD_2WAY_ll : ', KMODEL, &
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
    WRITE(*,*) 'Error UNSET_LSFIELD_2WAY_ll : ', KMODEL, &
               ' is not a child of the current model'
    STOP
  ENDIF
!
  CALL CLEANLIST_LCRSPD(TZLCOMDATA%TELT%TSEND_2WAY_LS)
  CALL CLEANLIST_LCRSPD(TZP2CDATA%TELT%TRECV_2WAY_LS)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UNSET_LSFIELD_2WAY_ll
!
!     #########################################
      SUBROUTINE LS_FORCING_ll( KCHILD, KINFO, OEXTRAPOL, OCYCLIC_EXTRAPOL )
!     #########################################
!!
!!****  *LS_FORCING_ll* - routine to do the forcing
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
!     Module MODE_NEST_ll
!       GO_TOMODEL_ll
!
!!    Implicit Arguments
!!    ------------------
!     Module MODD_STRUCTURE_ll
!       types LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll,
!             LPROCONF_ll, PROC_COM_DATA_ll,
!             PARENT2CHILD_DATA_ll
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
!     Original 11 fev. 2000
!         24/02/2015 (M.Moge) calling EXTRAPOL_ON_PSEUDO_HALO for cyclic cases where the child grid is the whole father grid
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LPROC_COM_DATA_ll, LPARENT2CHILD_DATA_ll, &
                                LPROCONF_ll, PROC_COM_DATA_ll,           &
                                PARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, TCRRT_PROCONF, NGRID_COM 
!
  USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_CRSPD, COPY_CRSPD
  USE MODE_NEST_ll, ONLY : GO_TOMODEL_ll
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODE_EXTRAPOL, ONLY : EXTRAPOL_ON_PSEUDO_HALO
  USE MODE_MODELN_HANDLER, ONLY : GOTO_MODEL
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN) :: KCHILD 
  INTEGER, INTENT(OUT) :: KINFO
  LOGICAL, OPTIONAL, INTENT(IN) :: OEXTRAPOL   !if TRUE, call EXTRAPOL_ON_PSEUDO_HALO
  LOGICAL, OPTIONAL, INTENT(IN) :: OCYCLIC_EXTRAPOL   !pass to EXTRAPOL_ON_PSEUDO_HALO, perform a cyclic extrapolation if TRUE
!
!
!*       0.2   declarations of local variables
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZCURRENT ! intermediate variables
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZCHILDCURRENT
  TYPE(LPROCONF_ll), POINTER :: TZPROCONFCURRENT
  TYPE(PROC_COM_DATA_ll), POINTER :: TZCHILD_COMDATA    ! child
  TYPE(PARENT2CHILD_DATA_ll), POINTER :: TZP2C_DATA
  INTEGER :: KINITIALMODEL, KINFO2
  TYPE(LIST_ll), POINTER :: TZLISTCURRENT
!
!-------------------------------------------------------------------------------
!
!*       1.    SEARCH THE KCHILD MODEL
!              -----------------------
!
  KINITIALMODEL = TCRRT_COMDATA%NUMBER 
  CALL GO_TOMODEL_ll(KCHILD, KINFO)
!
  IF (KINFO == -1) THEN
    RETURN
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*        2.   GO BACK ON ITS FATHER :
!              IT HAS THE STRUCTURES FOR LS COMMUNICATION
!              ------------------------------------------
!
!         2.1  Test if model has a parent (if not, it's the model 1)
!
  IF (.NOT.ASSOCIATED(TCRRT_COMDATA%TPARENT)) THEN
     KINFO = -2
     RETURN
  ENDIF
!
  TCRRT_COMDATA => TCRRT_COMDATA%TPARENT
  TCRRT_PROCONF => TCRRT_PROCONF%TPARENT
!  
!-------------------------------------------------------------------------------
!
!*        3.   POINT TO THE CHILD COMDATA, PROCONF
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
!*        4.   SEND THE LS FIELDS
!              ------------------
!
  CALL SEND_RECV_CRSPD(TZP2C_DATA%TSEND_1WAY_LS%TCRSPD, &
                       TZCHILD_COMDATA%TRECV_1WAY_LS%TCRSPD, &
                       TZP2C_DATA%TSEND_1WAY_LS%TLIST, &
                       TZCHILD_COMDATA%TRECV_1WAY_LS%TLIST, NGRID_COM, KINFO)
!
  CALL COPY_CRSPD(TZP2C_DATA%TSEND_1WAY_LS%TCRSPD, &
                  TZCHILD_COMDATA%TRECV_1WAY_LS%TCRSPD, &
                  TZP2C_DATA%TSEND_1WAY_LS%TLIST, &
                  TZCHILD_COMDATA%TRECV_1WAY_LS%TLIST, KINFO)
!
  CALL GO_TOMODEL_ll(KINITIALMODEL, KINFO2)
!
!  CALL GO_TOMODEL_ll(KCHILD, KINFO2)
!  CALL GOTO_MODEL(KCHILD)
  IF ( PRESENT(OEXTRAPOL) ) THEN
  IF ( OEXTRAPOL ) THEN
    TZLISTCURRENT => TZCHILD_COMDATA%TRECV_1WAY_LS%TLIST
    DO WHILE(ASSOCIATED(TZLISTCURRENT))
      IF( ASSOCIATED(TZLISTCURRENT%ARRAY3D) )THEN
        IF ( PRESENT(OCYCLIC_EXTRAPOL) ) THEN
          CALL EXTRAPOL_ON_PSEUDO_HALO(TZLISTCURRENT%ARRAY3D,OCYCLIC_EXTRAPOL)
        ELSE
          CALL EXTRAPOL_ON_PSEUDO_HALO(TZLISTCURRENT%ARRAY3D)
        ENDIF
      ENDIF
      IF( ASSOCIATED(TZLISTCURRENT%ARRAY2D) )THEN
        IF ( PRESENT(OCYCLIC_EXTRAPOL) ) THEN
          CALL EXTRAPOL_ON_PSEUDO_HALO(TZLISTCURRENT%ARRAY2D,OCYCLIC_EXTRAPOL)
        ELSE
          CALL EXTRAPOL_ON_PSEUDO_HALO(TZLISTCURRENT%ARRAY2D)
        ENDIF
      ENDIF
!      IF( ASSOCIATED(TZLISTCURRENT%ARRAY1D) )THEN
!        IF ( PRESENT(OCYCLIC_EXTRAPOL) ) THEN
!        CALL EXTRAPOL_ON_PSEUDO_HALO(TZLISTCURRENT%ARRAY1D,OCYCLIC_EXTRAPOL)
!        ELSE
!        CALL EXTRAPOL_ON_PSEUDO_HALO(TZLISTCURRENT%ARRAY1D)
!        ENDIF
!      ENDIF
      TZLISTCURRENT => TZLISTCURRENT%NEXT
    ENDDO
  ENDIF
  ENDIF
!  CALL GO_TOMODEL_ll(KINITIALMODEL, KINFO2)
!  CALL GOTO_MODEL(KINITIALMODEL)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE LS_FORCING_ll
!
!     ###################################
      SUBROUTINE LS_FEEDBACK_ll( KINFO )
!     ###################################
!
!!****  *LS_FEEDBACK_ll* - routine to do the feed back
!
!!    Purpose
!!    -------
!     This routine fill the fields of the parent model with the values
!     of the child model
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
!       types LPARENT2CHILD_DATA_ll, PROC_COM_DATA_ll,
!             PARENT2CHILD_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NGRID_COM - mpi communicator
!       NBUFFERSIZE_3D - buffer size
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
!     Original 11 fev. 2000
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LPARENT2CHILD_DATA_ll, PROC_COM_DATA_ll, &
                                PARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, TCRRT_PROCONF, NGRID_COM
!
  USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_CRSPD, COPY_CRSPD
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(OUT) :: KINFO
!
!*       0.2   declarations of local variables
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZCHILDCURRENT
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TZPARENT_COMDATA    ! child
  TYPE(PARENT2CHILD_DATA_ll), POINTER :: TZP2C_DATA
!
!-------------------------------------------------------------------------------
!
  KINFO = -1
!
!*        1.   POINT TO THE PARENT COMDATA, PROCONF
!              ------------------------------------
!
  TZPARENT_COMDATA => TCRRT_COMDATA%TPARENT
  TZCHILDCURRENT => TZPARENT_COMDATA%TP2C_DATA
!
!-------------------------------------------------------------------------------
!
!*        2.   FIND THE PROCONF STRUCTURE OF THE CURRENT CHILD
!              -----------------------------------------------
!
  DO WHILE(ASSOCIATED(TZCHILDCURRENT))
    IF(TZCHILDCURRENT%TELT%NUMBER .EQ. TCRRT_COMDATA%NUMBER) THEN
      EXIT
    ENDIF
    TZCHILDCURRENT => TZCHILDCURRENT%TNEXT
  ENDDO
!
  TZP2C_DATA => TZCHILDCURRENT%TELT
!
!-------------------------------------------------------------------------------
!
!*        3.   SEND AND RECEIVE THE LS FIELDS
!              ------------------------------
!
  CALL SEND_RECV_CRSPD(TCRRT_COMDATA%TSEND_2WAY_LS%TCRSPD, &
                       TZP2C_DATA%TRECV_2WAY_LS%TCRSPD, &
                       TCRRT_COMDATA%TSEND_2WAY_LS%TLIST, &
                       TZP2C_DATA%TRECV_2WAY_LS%TLIST, NGRID_COM, KINFO)
!
  CALL COPY_CRSPD(TCRRT_COMDATA%TSEND_2WAY_LS%TCRSPD, &
                  TZP2C_DATA%TRECV_2WAY_LS%TCRSPD, &
                  TCRRT_COMDATA%TSEND_2WAY_LS%TLIST, &
                  TZP2C_DATA%TRECV_2WAY_LS%TLIST, KINFO)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE LS_FEEDBACK_ll
!
END MODULE MODE_LS_ll
