!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
! Philippe 12/01/2018: renamed dimension variable NKMAX_ll->NKMAX_TMP_ll to prevent mix-up with modd_dimn
!-----------------------------------------------------------------

!     ###################
      MODULE MODE_NEST_ll
!     ###################
! 
!!    Purpose
!!    -------
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!     SUBROUTINES : GO_TOMODEL_ll, GET_CHILD_DIM_ll, GET_FEEDBACK_COORD_ll
!                   GET_MODEL_NUMBER_ll
!     FUNCTIONS   : LBFINE2COARSE
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
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_DIM_ll
!       NDXRATIO_ALL, NDYRATIO_ALL - Ratio for all models
!       NXOR_ALL, NXEND_ALL, NYOR_ALL, NYEND_ALL, NDAD
!
!     Module MODD_STRUCTURE_ll
!       types PARENT2CHILD_DATA_ll, LPARENT2CHILD_DATA_ll
!             PROC_COM_DATA_ll, LPROC_COM_DATA_ll
!             PROCONF_ll, LPROCONF_ll
!             CRSPD_ll, LCRSPD_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!       NPROC, DIMX, DIMY, IP, YSPLITTING
!  
!------------------------------------------------------------------------------
!  
  USE MODD_STRUCTURE_ll
!
  CONTAINS
!
!     #######################################
      SUBROUTINE ADD_PROCONF( TPHEAD, TPELT )
!     #######################################
!
!!****  *ADD_PROCONF* - routine to add a child at the end of a LPROCONF_ll
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     if the list is void, we create the list and put the element
!     as the first element else we add the element at the end of the list
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       type LPROCONF_ll
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
!     Original 22/07/98
!     Juan     19/08/2005: distinction Halo NORD/SUD & EST/WEST
!     J. Escobar 27/06/2011 correction for gridnesting with different SHAPE 
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : LPROCONF_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LPROCONF_ll), POINTER :: TPHEAD ! head of the list
!
  TYPE(PROCONF_ll), POINTER :: TPELT ! element to be added
!
!
!*       0.2   declarations of local variables
!
  TYPE(LPROCONF_ll), POINTER :: TZCURRENT, TZNEW ! intermediate variables
!
!-------------------------------------------------------------------------------
!
!*       1.    ADD THE ELEMENT TPELT :
!              ---------------------
!
  IF(.NOT.ASSOCIATED(TPHEAD)) THEN
!
! first element of the list
!
    ALLOCATE(TPHEAD)
    TPHEAD%TELT = TPELT
    NULLIFY( TPHEAD%TNEXT )
    TPHEAD%NCARD = 1
!
  ELSE
!
! others elements
!
    TZCURRENT => TPHEAD
!
!   Go to the end of the list
    DO WHILE(ASSOCIATED(TZCURRENT%TNEXT))
      TZCURRENT => TZCURRENT%TNEXT
    ENDDO
!
!   Add the element
    ALLOCATE(TZNEW)
    TZNEW%NCARD = TPHEAD%NCARD + 1
    TZNEW%TELT = TPELT
    NULLIFY(TZNEW%TNEXT)
!
    TZCURRENT%TNEXT => TZNEW
!
    TPHEAD%NCARD = TPHEAD%NCARD + 1
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE ADD_PROCONF
!
!     #####################################
      SUBROUTINE ADD_CHILD( TPHEAD, TPELT )
!     #####################################
!
!!****  *ADD_CHILD* - routine to add a child at the end of a LPROCONF_ll
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     if the list is void, we create the list and put the element
!     as the first element else we add the element at the end of the list
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       types PARENT2CHILD_DATA_ll, LPARENT2CHILD_DATA_ll
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
!     Original 22/07/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : PARENT2CHILD_DATA_ll, LPARENT2CHILD_DATA_ll
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TPHEAD ! head of the list
!
  TYPE(PARENT2CHILD_DATA_ll), INTENT(IN) :: TPELT ! element to be added
!
!
!*       0.2   declarations of local variables
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZCURRENT, TZNEW ! intermediate
                                                           ! variables
!
!-------------------------------------------------------------------------------
!
!*       1.    ADD THE ELEMENT TPELT :
!              ---------------------
!
  IF(.NOT.ASSOCIATED(TPHEAD)) THEN
!
! first element of the list
!
    ALLOCATE(TPHEAD)
    TPHEAD%TELT = TPELT
    NULLIFY( TPHEAD%TNEXT )
    TPHEAD%NCARD = 1
!
  ELSE
!
! others elements
!
    TZCURRENT => TPHEAD
!
!   Go to the end of the list
    DO WHILE(ASSOCIATED(TZCURRENT%TNEXT))
      TZCURRENT => TZCURRENT%TNEXT
    ENDDO
!
!   Add the element
    ALLOCATE(TZNEW)
    TZNEW%TELT = TPELT
    TZNEW%NCARD = TPHEAD%NCARD + 1
    NULLIFY(TZNEW%TNEXT)
!
    TZCURRENT%TNEXT => TZNEW
!
    TPHEAD%NCARD = TPHEAD%NCARD + 1
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE ADD_CHILD
!
!     #######################################
      SUBROUTINE ADD_COMDATA( TPHEAD, TPELT )
!     #######################################
!
!!****  *ADD_COMDATA* - routine to add a child at the end of a LPROC_COM_DATA_ll
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!     if the list is void, we create the list and put the element
!     as the first element else we add the element at the end of the list
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       types PROC_COM_DATA_ll, LPROC_COM_DATA_ll
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
!     Original 22/07/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : PROC_COM_DATA_ll, LPROC_COM_DATA_ll
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TPHEAD ! head of the list
!
  TYPE(PROC_COM_DATA_ll), INTENT(IN) :: TPELT ! element to be added
!
!
!*       0.2   declarations of local variables
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZCURRENT, TZNEW ! intermediate variables
!
!-------------------------------------------------------------------------------
!
!*       1.    ADD THE ELEMENT TPELT :
!              ---------------------
!
  IF(.NOT.ASSOCIATED(TPHEAD)) THEN
!
! first element of the list
!
    ALLOCATE(TPHEAD)
    TPHEAD%TELT = TPELT
    NULLIFY( TPHEAD%TNEXT )
    TPHEAD%NCARD = 1
!
  ELSE
!
! others elements
!
    TZCURRENT => TPHEAD
!
!   Go to the end of the list
    DO WHILE(ASSOCIATED(TZCURRENT%TNEXT))
      TZCURRENT => TZCURRENT%TNEXT
    ENDDO
!
!   Add the element
    ALLOCATE(TZNEW)
    TZNEW%TELT = TPELT
    TZNEW%NCARD = TPHEAD%NCARD + 1
    NULLIFY(TZNEW%TNEXT)
!
    TZCURRENT%TNEXT => TZNEW
!
    TPHEAD%NCARD = TPHEAD%NCARD + 1
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE ADD_COMDATA
!
!!    #######################################
      SUBROUTINE GO_TOMODEL_ll(KMODEL, KINFO)
!     #######################################
!
!!****  *GO_TOMODEL_ll* -Switch to the model KMODEL
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       types PROCONF_ll, PROC_COM_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     N. Gicquel     * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : PROCONF_ll, PROC_COM_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, TCRRT_PROCONF
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER :: KMODEL, KINFO
!
!*       0.2   declarations of local variables
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TZCURRENT_COM
  TYPE(PROCONF_ll), POINTER :: TZCURRENT_CONF
  INTEGER :: IFOUND
!
!-------------------------------------------------------------------------------
!
!         1.   GO BACK TO THE ROOT MODEL
!              -------------------------
!
  TZCURRENT_COM => TCRRT_COMDATA
  TZCURRENT_CONF => TCRRT_PROCONF
!
  DO WHILE(ASSOCIATED(TZCURRENT_COM%TPARENT))
    TZCURRENT_COM => TZCURRENT_COM%TPARENT
    TZCURRENT_CONF => TZCURRENT_CONF%TPARENT
  ENDDO
!
  IFOUND = 0
!
!-------------------------------------------------------------------------------
!
!         2.   SEARCH MODEL KMODEL
!              -------------------
!
  CALL SWITCH_TOMODEL(TZCURRENT_COM, TZCURRENT_CONF, KMODEL, IFOUND)
!
  IF (IFOUND == 0) THEN
    KINFO = -1 ! Model not found. Search error !!
  ELSE
    KINFO = 0
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GO_TOMODEL_ll
!
!     ####################################################################
      RECURSIVE SUBROUTINE SWITCH_TOMODEL(TPCURRENT_COM, TPCURRENT_CONF, &
                                          KMODEL, KFOUND)
!     ####################################################################
!
!!****  *SWITCH_TOMODEL* -
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       types PROCONF_ll, PROC_COM_DATA_ll, LPROC_COM_DATA_ll, LPROCONF_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!       TCRRT_PROCONF - Current configuration for current model
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     N. Gicquel     * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 02 fev. 1999
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : PROC_COM_DATA_ll, PROCONF_ll, &
                                LPROC_COM_DATA_ll, LPROCONF_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, TCRRT_PROCONF
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCURRENT_COM
  TYPE(PROCONF_ll), POINTER :: TPCURRENT_CONF
!
  INTEGER, INTENT(IN) :: KMODEL
  INTEGER :: KFOUND
!
!*       0.2   declarations of local variables
!
  INTEGER :: KSTART
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZCURRENT_COM ! intermediate variables
  TYPE(LPROCONF_ll), POINTER :: TZCURRENT_CONF ! intermediate variables
  TYPE(PROC_COM_DATA_ll), POINTER :: TZCHILD_COM
  TYPE(PROCONF_ll), POINTER ::  TZCHILD_CONF 
!
!-------------------------------------------------------------------------------
!
!         1.   IF THE CURRENT MODEL THE SEARCHED ONE?
!              --------------------------------------
!
  IF (TPCURRENT_COM%NUMBER == KMODEL) THEN
    KFOUND = 1
    TCRRT_COMDATA => TPCURRENT_COM
    TCRRT_PROCONF => TPCURRENT_CONF
    RETURN
  ENDIF
!
!-------------------------------------------------------------------------------
!
!         2.   NO => EXPLORE ALL THE CHILDREN MODEL OF THE CURRENT MODEL
!              ---------------------------------------------------------
!
  TZCURRENT_COM => TPCURRENT_COM%TCHILDREN
  TZCURRENT_CONF => TPCURRENT_CONF%TCHILDREN
  DO WHILE(ASSOCIATED(TZCURRENT_COM))
!
    IF (TZCURRENT_COM%TELT%NUMBER == KMODEL) THEN
      KFOUND = 1
      TCRRT_COMDATA => TZCURRENT_COM%TELT
      TCRRT_PROCONF => TZCURRENT_CONF%TELT
      RETURN
    ENDIF
!
    TZCHILD_COM => TZCURRENT_COM%TELT
    TZCHILD_CONF  => TZCURRENT_CONF%TELT
!
    CALL SWITCH_TOMODEL(TZCHILD_COM, TZCHILD_CONF, KMODEL, KFOUND)
!
    IF (KFOUND /= 1) THEN
      TZCURRENT_COM => TZCURRENT_COM%TNEXT
      TZCURRENT_CONF => TZCURRENT_CONF%TNEXT
    ELSE 
      RETURN
    ENDIF
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SWITCH_TOMODEL
!
!     ###############################################
      SUBROUTINE GET_MODEL_NUMBER_ll( KMODEL_NUMBER )
!     ###############################################
!
!!****  *GET_MODEL_NUMBER_ll* -
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
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
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch
!!
!!    Modifications
!!    -------------
!     Original 25/02/00
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER :: KMODEL_NUMBER
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
  KMODEL_NUMBER =TCRRT_COMDATA%NUMBER
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_MODEL_NUMBER_ll
!
!     ####################################################
      SUBROUTINE GET_CHILD_DIM_ll( KCHILD, KX, KY, KINFO )
!     ####################################################
!
!!****  *GET_CHILD_DIM* - routine to get the dimensions of model
!                         KCHILD (dimensions in the coarse grid
!                         reference of parent model)
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!     Module MODD_STRUCTURE_ll
!       type LPROC_COM_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!
!!    Author
!!    ------
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 30/07/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll
!
  USE MODD_STRUCTURE_ll, ONLY : LPROC_COM_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN) :: KCHILD
!
  INTEGER, INTENT(OUT) :: KX, KY
!
  INTEGER, INTENT(OUT) :: KINFO
!
!
!*       0.2   declarations of local variables
!
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZCURRENT ! intermediate variables
!
!-------------------------------------------------------------------------------
!
  KINFO = -1
  TZCURRENT => TCRRT_COMDATA%TCHILDREN
!
  DO WHILE(ASSOCIATED(TZCURRENT))
    IF(TZCURRENT%TELT%NUMBER .EQ. KCHILD) THEN
      KINFO = 0
      EXIT
    ENDIF
    TZCURRENT => TZCURRENT%TNEXT
  ENDDO
!
  IF(.NOT.ASSOCIATED(TZCURRENT)) RETURN
!
  KX = TZCURRENT%TELT%NLSDIMX
  KY = TZCURRENT%TELT%NLSDIMY
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_CHILD_DIM_ll
!
!     ###################################################################
      SUBROUTINE GET_FEEDBACK_COORD_ll( KXOR, KYOR, KXEND, KYEND, KINFO )
!     ###################################################################
!
!!****  *GET_FEEDBACK_COORD_ll* - routine to get the coordinates of model
!                                 KCHILD (dimensions in the coarse grid
!                                 reference of parent model)
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!       types CRSPD_ll, LPARENT2CHILD_DATA_ll
!
!     Module MODD_VAR_ll
!       TCRRT_COMDATA - Current communication data structure for current model
!                       and local processor
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch
!
!!    Modifications
!!    -------------
!     Original 06/10/99
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll, LPARENT2CHILD_DATA_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
!
  INTEGER, INTENT(OUT) :: KINFO
!
!
!*       0.2   declarations of local variables
!
!
  TYPE(LPARENT2CHILD_DATA_ll), POINTER :: TZCURRENT ! intermediate variables
  TYPE(CRSPD_ll), POINTER :: TZFEEDBACK ! intermediate variables
!
!-------------------------------------------------------------------------------
!
  KINFO = -1
!
  TZCURRENT => TCRRT_COMDATA%TPARENT%TP2C_DATA
!
  DO WHILE(ASSOCIATED(TZCURRENT) &
    .AND.(TZCURRENT%TELT%NUMBER.NE.TCRRT_COMDATA%NUMBER))
    TZCURRENT => TZCURRENT%TNEXT
  ENDDO
!
  IF(.NOT.ASSOCIATED(TZCURRENT)) RETURN
!
  TZFEEDBACK => TZCURRENT%TELT%TFEEDBACK_COORD
!
  IF(.NOT.ASSOCIATED(TZFEEDBACK)) RETURN
  KINFO = 0
!
  KXOR = TZFEEDBACK%TELT%NXOR
  KYOR = TZFEEDBACK%TELT%NYOR
  KXEND = TZFEEDBACK%TELT%NXEND
  KYEND = TZFEEDBACK%TELT%NYEND
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE GET_FEEDBACK_COORD_ll
!
!     #################################################################
      FUNCTION LBFINE2COARSE(KRATIO, KLBSIZE) RESULT(KCOARSE)
!     #################################################################
!
!
!!****  *LBFINE2COARSE* - This function computes the number of points
!                         on the coarse grid useful to interpolate
!                         the LB zone
!
!!    Purpose
!!    -------
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Gicquel     * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 09/10/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.0   declaration of result
!
  INTEGER :: KCOARSE
!
!*       0.1   declarations of arguments
!
  INTEGER :: KRATIO, & ! ratio between fine and coarse grid
             KLBSIZE   ! size of LB on the fine grid
!
!*       0.2   declarations of local variables  
!
  INTEGER :: I0
!
!!------------------------------------------------------------------------------
!
!*       1.    Calculate the number of points needed to interpolate
!              the LB zone
!              -----------
!
  I0 = CEILING(REAL(KRATIO)/2.) + 1 
!
  KCOARSE = 5 + FLOOR(REAL(KLBSIZE-I0)/REAL(KRATIO))
!
  IF (KRATIO == 1 .AND. KLBSIZE == 1) KCOARSE = 4
!
!!------------------------------------------------------------------------------
!
      END FUNCTION LBFINE2COARSE
!
!     ################################################
      SUBROUTINE CHILDS_PT_OF_VIEW(TPRECVZONE, TPPE, &
                                   TPORIG, KRECVTAG, &
                                   TPRECVCRSPD)
!     ################################################
!
!!****  *CHILDS_PT_OF_VIEW*-
!
!!    Purpose
!!    -------
!     Given the zone to be sent and received by the current
!     child processor (TPRECVZONE), this routine initializes the
!     communication data structure TPRECVCRSPD, communication
!     that the current processor has to make as a child
!!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       INTERSECTION, ADD_ZONE
!
!   Implicit Arguments
!!   ------------------
!
!     Module MODD_STRUCTURE_ll
!       types ZONE_ll, CRSPD_ll
!
!     Module MODD_VAR_ll
!       NPROC
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 22/07/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, CRSPD_ll
  USE MODD_VAR_ll, ONLY : NPROC
  USE MODE_TOOLS_ll, ONLY : INTERSECTION, ADD_ZONE
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(ZONE_ll), INTENT(IN) :: TPRECVZONE ! zone to be received
!
  TYPE(ZONE_ll), DIMENSION(:),  INTENT(IN) :: TPPE   ! physical splitting
                                                     ! +boundary halo
!
  INTEGER, INTENT(IN) :: KRECVTAG ! MPI message tags
!
  TYPE(ZONE_ll), INTENT(IN) :: TPORIG ! coordinates of the child's array
                                      ! in the parent's model
!
  TYPE(CRSPD_ll), POINTER :: TPRECVCRSPD
!
!*       0.2   declarations of local variables
!
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZINTER 
  INTEGER :: JI
!
!-------------------------------------------------------------------------------
!
!*        1.    Construction of the 1way communication data
!               --------------------------------------------
!
  CALL INTERSECTION(TPPE, NPROC, TPRECVZONE, TZINTER)
!
  DO JI = 1, NPROC
    IF (TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = KRECVTAG
      TZINTER(JI)%NXOR = TZINTER(JI)%NXOR - TPORIG%NXOR + 1
      TZINTER(JI)%NYOR = TZINTER(JI)%NYOR - TPORIG%NYOR + 1
      TZINTER(JI)%NXEND = TZINTER(JI)%NXEND  - TPORIG%NXOR + 1
      TZINTER(JI)%NYEND = TZINTER(JI)%NYEND - TPORIG%NYOR + 1
      CALL ADD_ZONE( TPRECVCRSPD, TZINTER(JI) )
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE CHILDS_PT_OF_VIEW
!
!     #######################################################
      SUBROUTINE PARENTS_PT_OF_VIEW(TPSENDZONES, TPPE,      &
                                    TPPROCONF, KSENDTAG,    &
                                    TPSENDLCRSPD, OINTER)
!     #######################################################
!
!!****  *PARENTS_PT_OF_VIEW*-
!
!!    Purpose
!!    -------
!     Given the zones to be sent and received to each child processor
!     (TPSENDZONES), this routine initializes the communication
!     data structure TPSENDCRSPD, communication that the current
!     processor has to make as a parent
!!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!       INTERSECTION, ADD_ZONE, GLOBAL2LOCAL
!
!   Implicit Arguments
!!   ------------------
!
!     Module MODD_STRUCTURE_ll
!       types ZONE_ll, CRSPD_ll
!
!     Module MODD_VAR_ll
!       NPROC
!
!!    Reference
!!    ---------
!     Notes sur les routines multi-modeles de la surcouche, R. Guivarch
!
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!     Ph. Kloos      * CNRM - CERFACS *
!
!!    Modifications
!!    -------------
!     Original 22/07/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, PROCONF_ll, LCRSPD_ll
  USE MODD_VAR_ll, ONLY : NPROC
  USE MODE_TOOLS_ll, ONLY : ADD_ZONE, INTERSECTION, GLOBAL2LOCAL
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(ZONE_ll), DIMENSION(:), INTENT(IN) :: TPSENDZONES ! zones to be sent
!
  TYPE(ZONE_ll),   INTENT(IN) :: TPPE   ! physical splitting+boundary halo
!
  TYPE(PROCONF_ll), POINTER :: TPPROCONF
!
  INTEGER, INTENT(IN) :: KSENDTAG ! MPI message tags
!
  TYPE(LCRSPD_ll), POINTER :: TPSENDLCRSPD
!
  LOGICAL, INTENT(OUT) :: OINTER
!
!*       0.2   declarations of local variables
!
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZINTER 
  INTEGER :: JI
!
  TYPE(ZONE_ll), POINTER :: TZZONE
  TYPE(CRSPD_ll), POINTER :: TZCRSPD
!
!-------------------------------------------------------------------------------
!
!*        1.    Construction of the 1way communication data
!               --------------------------------------------
!
!*        1.1   Compute the intersection between the subdomain
!               of the current proc (coarse mesh)
!               and the set of LB zones corresponding to 
!               each child proc (coarse mesh to be interpolated)
!
  CALL INTERSECTION(TPSENDZONES, NPROC, TPPE, TZINTER)
!
!*        1.2   Add the intersection to TPSENDLCRSPD%TCRSPD structure
!
  OINTER = .FALSE.
  DO JI = 1, NPROC
    IF (TZINTER(JI)%NUMBER.NE.0) THEN
      IF (.NOT. OINTER) THEN
        NULLIFY(TPSENDLCRSPD%TCRSPD)
        OINTER = .TRUE.
      ENDIF
      TZINTER(JI)%MSSGTAG = KSENDTAG
      CALL ADD_ZONE( TPSENDLCRSPD%TCRSPD, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       2.    Switch to local coordinates
!              ---------------------------
!
  IF (OINTER) THEN
    CALL GLOBAL2LOCAL(TPPROCONF, TPSENDLCRSPD%TCRSPD)
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE PARENTS_PT_OF_VIEW
!
!     #########################################################
      RECURSIVE SUBROUTINE INI_CHILD( TPPROCONF, TPCOMDATA, K )
!     #########################################################
!
!!****  *INI_CHILD* - routine to construct recursively the data structures
!                     of all the model
!
!!    Purpose
!!    -------
!     this routine fills all the data structures (splitting data structures,
!     communications data structures) for grid-nesting
!
!!**  Method
!!    ------
!
!!    External
!!    --------
!
!     Module MODE_SPLITTING_ll
!       SPLIT2      
!
!     Module MODE_TOOLS_ll 
!       LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll,
!       INTERSECTION, GLOBAL2LOCAL, ADD_ZONE, EXTRACT_ZONE
!
!     Module MODE_CONSTRUCT_ll
!       INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS
!       CONSTRUCT_HALO1, CONSTRUCT_HALO2,
!       CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_DIM_ll
!       NDXRATIO_ALL, NDYRATIO_ALL - Ratio for all models
!       NXOR_ALL, NXEND_ALL, NYOR_ALL, NYEND_ALL, NDAD
!
!     Module MODD_STRUCTURE_ll
!       types ZONE_ll, PROCONF_ll, PROC_COM_DATA_ll
!             PARENT2CHILD_DATA_ll, CRSPD_ll, LCRSPD_ll
!
!     Module MODD_VAR_ll
!       NPROC, DIMX, DIMY, IP, YSPLITTING
!
!!    Reference
!!    ---------
!     Notes sur les routines multi-modeles de la surcouche, R. Guivarch
!
!!    Author
!!    ------
!     R. Guivarch    * CERFACS *
!
!!    Modifications
!!    -------------
!     Original 22/07/98
!     R. Guivarch 29/11/99  x and y splitting -> YSPLITTING 
!     J. Escobar  24/09/2013 : temp patch for problem of gridnesting with different SHAPE
!     M.Moge      10/02/2015 construct halo extended (needed for an interpolation in SPAWNING)
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_DIM_ll, ONLY : NDXRATIO_ALL, NXOR_ALL, NXEND_ALL, &
                          NDYRATIO_ALL, NYOR_ALL, NYEND_ALL, NDAD, &
                          NKMAX_TMP_ll
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT, JPMODELMAX
!
  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, PROCONF_ll, PROC_COM_DATA_ll, &
                                PARENT2CHILD_DATA_ll, CRSPD_ll, LCRSPD_ll
!
  USE MODD_VAR_ll, ONLY : NPROC, DIMX, DIMY, DIMZ, IP, YSPLITTING
!JUANZ
  USE MODD_CONFZ, ONLY  : NZ_VERB,NZ_PROC
!JUANZ
!
  USE MODE_TOOLS_ll, ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll, &
        INTERSECTION, GLOBAL2LOCAL, ADD_ZONE, EXTRACT_ZONE
!
  USE MODE_CONSTRUCT_ll, ONLY : INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS, &
                                CONSTRUCT_HALO1, CONSTRUCT_HALO2, CONSTRUCT_HALO_EXTENDED, &
                                CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY
!
  USE MODE_SPLITTING_ll , ONLY : SPLIT2, def_splitting2
  USE MODE_TOOLSZ_ll    , ONLY : ini_pzz, ini_boundariesz, ini_ezz, construct_transz
  USE MODE_TOOLSZ_ll    , ONLY : splitz
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN) :: K ! Number of the model
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data structure
  TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
!
!*       0.2   declarations of local variables
!
  TYPE(PROCONF_ll), POINTER :: TZCHILD_PROCONF ! variable for the child
                                               ! splitting data
!
  TYPE(PARENT2CHILD_DATA_ll), POINTER :: TZP2C_DATA ! variable for parent
                                                    ! to child communications
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TZCHILD_COMDATA ! variable for the child
                                                     ! communication data
                                                     ! (parent to child, halo,
                                                     ! transposition)
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZCOARSE ! Coarse grid splitting
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZFINE   ! Fine grid splitting
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZINTER ! Intermediate zone array
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZP ! Physical zone splitting
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZPB ! Physical + boundary halo
                                                   ! splitting
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:)  :: TZSEND, TZRECV ! Points of
!
  TYPE(ZONE_ll) :: TZFEEDBACK
                              ! the coarse grid
                              ! to be sent or to be received
                              ! to (from) the child
!
  TYPE(LCRSPD_ll), POINTER :: TZLCRSPD, & ! temporary list of correspondants
                              TZLCRSPDS, TZLCRSPDR
!
  TYPE(CRSPD_ll), POINTER :: TZCRSPDS, TZCRSPDR, TZPTR
!
  INTEGER :: J, JMODEL, JI ! loop controls
!
  INTEGER :: IERR
  INTEGER :: ICARD
  INTEGER :: IXMIN, IYMIN, IXMAX, IYMAX

!JUAN Z_SPLITTING
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP1_YP2_Z ! intermediate Full Z = B splitting  without halo zone
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SX_YP2_ZP1 ! intermediate Full X     splitting zone
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_Y_ZP1 ! intermediate Full Y     splitting zone
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_YP1_Z ! intermediate Full Z = B transposed splitting  without halo zone
  INTEGER :: JX_DOMAINS,JY_DOMAINS
  LOGICAL :: LPREM
  INTEGER :: P1,P2
!JUAN Z_SPLITTING
 INTEGER :: IXSIZE,IYSIZE,IXSIZE_COARSE,IYSIZE_COARSE ! size of child without halo
!
!-------------------------------------------------------------------------------
!
!
!*       1.    ALLOCATION AND NUMBER INITIALIZATION :
!              ------------------------------------
!
  ALLOCATE(TZCHILD_PROCONF)
  CALL ALLOC(TZP2C_DATA)
  CALL ALLOC(TZCHILD_COMDATA)
!
  ALLOCATE(TZINTER(NPROC))
!JUAN Z_SPLITTING
  ALLOCATE(TZDZP_SXP1_YP2_Z(NPROC))
  ALLOCATE(TZDZP_SXP2_YP1_Z(NPROC))
  ALLOCATE(TZDZP_SX_YP2_ZP1(NPROC))
  ALLOCATE(TZDZP_SXP2_Y_ZP1(NPROC))
!JUAN Z_SPLITTING
!
  NULLIFY(TZCRSPDS)
  NULLIFY(TZCRSPDR)
!
!        1.1   Number of the child model
!
  TZCHILD_PROCONF%NUMBER = K
  TZP2C_DATA%NUMBER = K
  TZCHILD_COMDATA%NUMBER = K
!
!        1.2   TZCHILD_PROCONF%TPARENT and TZCHILD_COMDATA%TPARENT
!              point to the parent structure
!
  TZCHILD_PROCONF%TPARENT => TPPROCONF
  TZCHILD_COMDATA%TPARENT => TPCOMDATA
!
!        1.3   Dimensions of the child model
!
  IXSIZE_COARSE=NXEND_ALL(K)-NXOR_ALL(K) -2*JPHEXT +1
  IYSIZE_COARSE=NYEND_ALL(K)-NYOR_ALL(K) -2*JPHEXT +1
  IXSIZE = NDXRATIO_ALL(K) * IXSIZE_COARSE
  IYSIZE = NDYRATIO_ALL(K) * IYSIZE_COARSE
  DIMX = IXSIZE + 2*JPHEXT
  DIMY = IYSIZE + 2*JPHEXT
!JUAN Z_SPLITTING
  DIMZ = NKMAX_TMP_ll + 2*JPVEXT
!JUAN Z_SPLITTING
!
!*       2.    CONSTRUCTION OF USEFUL SPLITTINGS :
!              ---------------------------------
!
!        2.1   TZCOARSE splitting of the coarse grid including
!              the child model zone
!
  ALLOCATE(TZCOARSE(NPROC))
!!$  CALL SPLIT2(NXEND_ALL(K)-NXOR_ALL(K)-1, NYEND_ALL(K)-NYOR_ALL(K)-1, &
!!$!JUAN
!!$       NKMAX_TMP_ll,&
!!$!JUAN
!!$       NPROC,TZCOARSE, YSPLITTING)
    !
    ! find the B splitting, dimension without halo for FFT
    !
    CALL DEF_SPLITTING2(JX_DOMAINS,JY_DOMAINS,IXSIZE,IYSIZE,NPROC,LPREM)
    !
    P1 = MIN(DIMZ,JX_DOMAINS)
    !JUAN PATCH NESTING DIFFERENT SHAPE
    P1 = NZ_PROC
    !JUAN PATCH NESTING DIFFERENT SHAPE
    P2 = NPROC / P1
    IF ( NZ_VERB .GE. 5 ) THEN
       IF ( IP .EQ. 1 )THEN
          print*," INI_CHILDZ:: NZ_PROC   =",NZ_PROC
          print*," INI_CHILDZ:: JX_DOMAINS=",JX_DOMAINS
          print*," INI_CHILDZ:: JY_DOMAINS=",JY_DOMAINS
          print*
          !
          print*," INI_CHILDZ:: P1=MIN(DIMZ,MAX(JX_DOMAINS,JY_DOMAINS))=",  P1
          !
          print*," INI_CHILDZ:: P2=NPROC/P1/                           =",  P2
       END IF
    END IF
    !
    ! find Z splitting wihout halo in X&Y , with halo in Z
    !
    CALL SPLIT2(IXSIZE_COARSE,IYSIZE_COARSE, NKMAX_TMP_ll, NPROC,TZCOARSE, YSPLITTING,P1,P2)

    CALL SPLITZ(IXSIZE_COARSE,IYSIZE_COARSE,DIMZ,NPROC,TZDZP_SXP1_YP2_Z,'P1P2SPLITT', 1 ,P1,P2)
    CALL SPLITZ(IXSIZE_COARSE,IYSIZE_COARSE,DIMZ,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING', P1,P1,P2)
    CALL SPLITZ(IXSIZE_COARSE,IYSIZE_COARSE,DIMZ,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING', P1,P1,P2)
    CALL SPLITZ(IXSIZE_COARSE,IYSIZE_COARSE,DIMZ,NPROC,TZDZP_SXP2_YP1_Z,'P2P1SPLITT', 1 ,P1,P2)

    CALL COARSE_TO_FINE(TZDZP_SXP1_YP2_Z)
    CALL COARSE_TO_FINE(TZDZP_SX_YP2_ZP1)
    CALL COARSE_TO_FINE(TZDZP_SXP2_Y_ZP1)
    CALL COARSE_TO_FINE(TZDZP_SXP2_YP1_Z)

!!$    CALL SPLITZ(DIMX-2*JPHEXT,DIMY-2*JPHEXT,DIMZ,NPROC,TZDZP_SXP1_YP2_Z,'P1P2SPLITT', 1 ,P1,P2)
!!$    CALL SPLITZ(DIMX-2*JPHEXT,DIMY-2*JPHEXT,DIMZ,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING', P1,P1,P2)
!!$    CALL SPLITZ(DIMX-2*JPHEXT,DIMY-2*JPHEXT,DIMZ,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING', P1,P1,P2)
!!$    CALL SPLITZ(DIMX-2*JPHEXT,DIMY-2*JPHEXT,DIMZ,NPROC,TZDZP_SXP2_YP1_Z,'P2P1SPLITT', 1 ,P1,P2)


!
!        2.2   TZFINE splitting of fine grid of the child model
!                     (physical zones)
!
  ALLOCATE(TZFINE(NPROC))
  DO J = 1, NPROC
!
    TZFINE(J)%NUMBER = TZCOARSE(J)%NUMBER

    TZFINE(J)%NXOR  = (TZCOARSE(J)%NXOR  - JPHEXT -1 ) * NDXRATIO_ALL(K) + JPHEXT +1
    TZFINE(J)%NYOR  = (TZCOARSE(J)%NYOR  - JPHEXT -1 ) * NDYRATIO_ALL(K) + JPHEXT +1
    TZFINE(J)%NXEND = (TZCOARSE(J)%NXEND - JPHEXT    ) * NDXRATIO_ALL(K) + JPHEXT
    TZFINE(J)%NYEND = (TZCOARSE(J)%NYEND - JPHEXT    ) * NDYRATIO_ALL(K) + JPHEXT

!!$    TZFINE(J)%NXOR = (TZCOARSE(J)%NXOR - 2) * NDXRATIO_ALL(K) + 2
!!$    TZFINE(J)%NYOR = (TZCOARSE(J)%NYOR - 2) * NDYRATIO_ALL(K) + 2
!!$    TZFINE(J)%NXEND = (TZCOARSE(J)%NXEND - 1) * NDXRATIO_ALL(K) + 1
!!$    TZFINE(J)%NYEND = (TZCOARSE(J)%NYEND - 1) * NDYRATIO_ALL(K) + 1

!JUAN Z_SPLITTING
    TZFINE(J)%NZOR  = TZCOARSE(J)%NZOR
    TZFINE(J)%NZEND = TZCOARSE(J)%NZEND
!JUAN Z_SPLITTING
!
  ENDDO
!
!        2.3   TZPB splitting of the physical zone + the boundary halo
!              TZP splitting of the physical zone
!
  ALLOCATE(TZP(NPROC), TZPB(NPROC))
!
  CALL EXTRACT_ZONE( TPPROCONF%TSPLITS_B, TZP, TZPB )
!
  DO J = 1, NPROC
!
!!$    IF(.NOT.LNORTH_ll(J)) TZPB(J)%NYEND = TZP(J)%NYEND
!!$    IF(.NOT.LSOUTH_ll(J)) TZPB(J)%NYOR = TZP(J)%NYOR
!!$    IF(.NOT.LWEST_ll(J))  TZPB(J)%NXOR = TZP(J)%NXOR
!!$    IF(.NOT.LEAST_ll(J))  TZPB(J)%NXEND = TZP(J)%NXEND

  IF (.NOT.TPPROCONF%TBOUND(J)%NORTH)  TZPB(J)%NYEND = TZP(J)%NYEND
  IF (.NOT.TPPROCONF%TBOUND(J)%SOUTH)  TZPB(J)%NYOR  = TZP(J)%NYOR
  IF (.NOT.TPPROCONF%TBOUND(J)%WEST)   TZPB(J)%NXOR  = TZP(J)%NXOR
  IF (.NOT.TPPROCONF%TBOUND(J)%EAST)   TZPB(J)%NXEND = TZP(J)%NXEND

!
  ENDDO
!
!        2.4   TZSEND  points of the coarse grid to be sent to the child
!              for each processor
!
  ALLOCATE(TZSEND(NPROC))
!
  DO J = 1, NPROC
!
    TZSEND(J)%NUMBER = TZCOARSE(J)%NUMBER

    TZSEND(J)%NXOR = NXOR_ALL(K) + TZCOARSE(J)%NXOR    -1  -JPHEXT -1 ! - 3
    TZSEND(J)%NYOR = NYOR_ALL(K) + TZCOARSE(J)%NYOR    -1  -JPHEXT -1 ! - 3
    TZSEND(J)%NXEND = NXOR_ALL(K) + TZCOARSE(J)%NXEND  -1  +JPHEXT +1 ! + 1
    TZSEND(J)%NYEND = NYOR_ALL(K) + TZCOARSE(J)%NYEND  -1  +JPHEXT +1 ! + 1
!!$
!!$    TZSEND(J)%NXOR = NXOR_ALL(K) + TZCOARSE(J)%NXOR - 3
!!$    TZSEND(J)%NYOR = NYOR_ALL(K) + TZCOARSE(J)%NYOR - 3
!!$    TZSEND(J)%NXEND = NXOR_ALL(K) + TZCOARSE(J)%NXEND + 1
!!$    TZSEND(J)%NYEND = NYOR_ALL(K) + TZCOARSE(J)%NYEND + 1
!
  ENDDO
!
  TZCHILD_COMDATA%NLSDIMX = TZCOARSE(IP)%NXEND - TZCOARSE(IP)%NXOR + 1 +2*(JPHEXT+1) ! + 5
  TZCHILD_COMDATA%NLSDIMY = TZCOARSE(IP)%NYEND - TZCOARSE(IP)%NYOR + 1 +2*(JPHEXT+1) ! + 5
!!$
!!$  TZCHILD_COMDATA%NLSDIMX = TZCOARSE(IP)%NXEND - TZCOARSE(IP)%NXOR + 5
!!$  TZCHILD_COMDATA%NLSDIMY = TZCOARSE(IP)%NYEND - TZCOARSE(IP)%NYOR + 5
!
!        2.5   TZRECV points of the coarse grid to be received from the child
!              for each processor
!
  ALLOCATE(TZRECV(NPROC))
!
  DO J = 1, NPROC
!
    TZRECV(J)%NUMBER = TZCOARSE(J)%NUMBER
    TZRECV(J)%NXOR = NXOR_ALL(K) + TZCOARSE(J)%NXOR - 1
    TZRECV(J)%NYOR = NYOR_ALL(K) + TZCOARSE(J)%NYOR - 1
    TZRECV(J)%NXEND = NXOR_ALL(K) + TZCOARSE(J)%NXEND - 1
    TZRECV(J)%NYEND = NYOR_ALL(K) + TZCOARSE(J)%NYEND - 1
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       3.    PARENT'S POINT OF VIEW :
!              ----------------------
!
!        3.1   Construction of TZP2C_DATA%TSEND_1WAY_LS, the variable of type
!              CRSPD_ll which contains data for 1way communication
!              between the parent and its child
!
!        3.1.1 Computation of the zones the parent sends
!
  CALL INTERSECTION( TZSEND, NPROC, TZPB(IP), TZINTER )
!
  CALL ALLOC(TZP2C_DATA%TSEND_1WAY_LS)
!
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 200
      CALL ADD_ZONE( TZP2C_DATA%TSEND_1WAY_LS%TCRSPD, TZINTER(JI) )
    ENDIF
  ENDDO
!
!        3.1.2 Switch to local coordinates
!
  CALL GLOBAL2LOCAL(TPPROCONF, TZP2C_DATA%TSEND_1WAY_LS%TCRSPD)
!
!        3.2   Construction of TZP2C_DATA%TRECV_2WAY_LS%TCRSPD and 
!              TZP2C_DATA%TFEEDBACK_COORD, the variables of type
!              CRSPD_ll which contain data for 2way communication
!              between the parent and its child
!
!        3.2.1 Computation of the zones the parent receives
!
  CALL INTERSECTION( TZRECV, NPROC, TZP(IP), TZINTER )
!
  CALL ALLOC(TZP2C_DATA%TRECV_2WAY_LS)
!
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 300
      CALL ADD_ZONE( TZP2C_DATA%TRECV_2WAY_LS%TCRSPD, TZINTER(JI) )
    ENDIF
  ENDDO
!
  TZPTR => TZP2C_DATA%TRECV_2WAY_LS%TCRSPD
!
!        3.2.2 Search of the extremum of these zones
!
  IF(ASSOCIATED(TZPTR)) THEN

    IXMIN = TZPTR%TELT%NXOR
    IYMIN = TZPTR%TELT%NYOR
    IXMAX = TZPTR%TELT%NXEND
    IYMAX = TZPTR%TELT%NYEND

    DO WHILE(ASSOCIATED(TZPTR))
      IF(TZPTR%TELT%NXOR < IXMIN) IXMIN = TZPTR%TELT%NXOR
      IF(TZPTR%TELT%NYOR < IYMIN) IYMIN = TZPTR%TELT%NYOR
      IF(TZPTR%TELT%NXEND > IXMAX) IXMAX = TZPTR%TELT%NXEND
      IF(TZPTR%TELT%NYEND > IYMAX) IYMAX = TZPTR%TELT%NYEND
      TZPTR => TZPTR%TNEXT
    ENDDO
!
    TZFEEDBACK%NXOR = IXMIN
    TZFEEDBACK%NYOR = IYMIN
    TZFEEDBACK%NXEND = IXMAX
    TZFEEDBACK%NYEND = IYMAX
!
!
!        3.2.3 Memorization of these extremum in TZP2C_DATA%TFEEDBACK_COORD
!
    CALL ADD_ZONE( TZP2C_DATA%TFEEDBACK_COORD, TZFEEDBACK)
!
!        3.2.4 Modification of TZP2C_DATA%TRECV_2WAY_LS%TCRSPD
!		(translation with (-IXMIN+1,-IYMIN+1) vector)
!		in order to haves indices between 1 and IXMAX or IYMAX
!
    TZPTR => TZP2C_DATA%TRECV_2WAY_LS%TCRSPD
    DO WHILE(ASSOCIATED(TZPTR))
      TZPTR%TELT%NXOR = TZPTR%TELT%NXOR - IXMIN + 1
      TZPTR%TELT%NYOR = TZPTR%TELT%NYOR - IYMIN + 1
      TZPTR%TELT%NXEND = TZPTR%TELT%NXEND - IXMIN + 1
      TZPTR%TELT%NYEND = TZPTR%TELT%NYEND - IYMIN + 1
      TZPTR => TZPTR%TNEXT
    ENDDO
!
!        3.2.5 Switch to local coordinates
!
    CALL GLOBAL2LOCAL(TPPROCONF, TZP2C_DATA%TFEEDBACK_COORD)
!
  ENDIF
!
!        3.5   Set TZCHILD_COMDATA%NRIMLBX and TZCHILD_COMDATA%NRIMLBY,
!              the arrays which contain the sizes of the LB zones
!
!  CALL SET_LB_SIZES(K, TZCHILD_COMDATA%NDIMRIMLBX, TZCHILD_COMDATA%NDIMRIMLBY, &
!                    TZCHILD_COMDATA%NRIMLBX, TZCHILD_COMDATA%NRIMLBY)
!
!        3.6   Construction of TZP2C_DATA%TSEND_1WAY_LBX(W/E) and
!              TZP2C_DATA%TSEND_1WAY_LBY(S/N), the variables of type
!              LCRSPD_ll which contains the regions to be sent to the child
!              model procs to build the LB fields in the x direction
!
!-------------------------------------------------------------------------------
!
!*       4.    CHILD'S POINT OF VIEW :
!              ---------------------
!
!        4.1   Construction of TZCHILD_PROCONF, variable which contains
!              all the splittings of the child model
!
  ALLOCATE(TZCHILD_PROCONF%TSPLITS_B(NPROC))
  ALLOCATE(TZCHILD_PROCONF%TSPLITS_X(NPROC))
  ALLOCATE(TZCHILD_PROCONF%TSPLITS_Y(NPROC))
!JUAN Z_SPLITTING
  ALLOCATE(TZCHILD_PROCONF%TSPLITS_SXP1_YP2_Z(NPROC))
  ALLOCATE(TZCHILD_PROCONF%TSPLITS_SXP2_YP1_Z(NPROC))
  ALLOCATE(TZCHILD_PROCONF%TSPLITS_SX_YP2_ZP1(NPROC))
  ALLOCATE(TZCHILD_PROCONF%TSPLITS_SXP2_Y_ZP1(NPROC))
!JUAN Z_SPLITTING
  ALLOCATE(TZCHILD_PROCONF%TBOUND(NPROC))
!
!
  CALL INI_PZ(TZCHILD_PROCONF, TZFINE)
!JUAN Z_SPLITTING
  CALL INI_PZZ(TZCHILD_PROCONF%TSPLITS_SXP1_YP2_Z,TZDZP_SXP1_YP2_Z)
  CALL INI_PZZ(TZCHILD_PROCONF%TSPLITS_SXP2_YP1_Z,TZDZP_SXP2_YP1_Z)
  CALL INI_PZZ(TZCHILD_PROCONF%TSPLITS_SX_YP2_ZP1,TZDZP_SX_YP2_ZP1)
  CALL INI_PZZ(TZCHILD_PROCONF%TSPLITS_SXP2_Y_ZP1,TZDZP_SXP2_Y_ZP1)
!JUAN Z_SPLITTING
!
  CALL INI_BOUNDARIES(TZCHILD_PROCONF)
!JUAN Z_SPLITTING
  CALL INI_BOUNDARIESZ(TZCHILD_PROCONF)
!JUAN Z_SPLITTING
!
  CALL INI_EZ(TZCHILD_PROCONF)
!JUAN Z_SPLITTING
  CALL INI_EZZ(TZCHILD_PROCONF)
!JUAN Z_SPLITTING
!
  CALL INI_TRANS(TZCHILD_PROCONF)
!
!        4.2   Construction of TZCHILD_COMDATA (LS communications)
!
!        4.2.1 Construction of TZCHILD_COMDATA%TRECV_1WAY_LS%TCRSPD,
!              the variable of type CRSPD_ll which contains
!              data for 1WAY_LS communication between the child and its parent
!
  CALL INTERSECTION(TZPB, NPROC, TZSEND(IP), TZINTER)
!
  CALL ALLOC(TZCHILD_COMDATA%TRECV_1WAY_LS)
!
!  NULLIFY(TZCHILD_COMDATA%TRECV_1WAY_LS%TCRSPD)
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 200
      TZINTER(JI)%NXOR = TZINTER(JI)%NXOR - TZSEND(IP)%NXOR + 1
      TZINTER(JI)%NYOR = TZINTER(JI)%NYOR - TZSEND(IP)%NYOR + 1
      TZINTER(JI)%NXEND = TZINTER(JI)%NXEND  - TZSEND(IP)%NXOR + 1
      TZINTER(JI)%NYEND = TZINTER(JI)%NYEND - TZSEND(IP)%NYOR + 1
      CALL ADD_ZONE( TZCHILD_COMDATA%TRECV_1WAY_LS%TCRSPD, TZINTER(JI) )
    ENDIF
  ENDDO
!
!        4.2.2 Construction of TZCHILD_COMDATA%TSEND_2WAY_LS%TCRSPD,
!              the variable of type CRSPD_ll which contains
!              data for 2WAY_LS communication between the child and its parent
!
  CALL INTERSECTION( TZP, NPROC, TZRECV(IP), TZINTER )
!
  CALL ALLOC(TZCHILD_COMDATA%TSEND_2WAY_LS)
!
!  NULLIFY(TZCHILD_COMDATA%TSEND_2WAY_LS%TCRSPD)
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 300
      TZINTER(JI)%NXOR = TZINTER(JI)%NXOR - TZSEND(IP)%NXOR + 1
      TZINTER(JI)%NYOR = TZINTER(JI)%NYOR - TZSEND(IP)%NYOR + 1
      TZINTER(JI)%NXEND = TZINTER(JI)%NXEND  - TZSEND(IP)%NXOR + 1
      TZINTER(JI)%NYEND = TZINTER(JI)%NYEND - TZSEND(IP)%NYOR + 1
      CALL ADD_ZONE( TZCHILD_COMDATA%TSEND_2WAY_LS%TCRSPD, TZINTER(JI) )
    ENDIF
  ENDDO
!
!        4.4 Construction of the part of TZCHILD_COMDATA concerning data
!        for halo and transposition communications
!
  TZCHILD_COMDATA%TSPLIT_B => TZCHILD_PROCONF%TSPLITS_B(IP)
  TZCHILD_COMDATA%TSPLIT_X => TZCHILD_PROCONF%TSPLITS_X(IP)
  TZCHILD_COMDATA%TSPLIT_Y => TZCHILD_PROCONF%TSPLITS_Y(IP)
!JUAN Z_SPLITTING
  TZCHILD_COMDATA%TSPLIT_SXP1_YP2_Z => TZCHILD_PROCONF%TSPLITS_SXP1_YP2_Z(IP)
  TZCHILD_COMDATA%TSPLIT_SXP2_YP1_Z => TZCHILD_PROCONF%TSPLITS_SXP2_YP1_Z(IP)
  TZCHILD_COMDATA%TSPLIT_SX_YP2_ZP1 => TZCHILD_PROCONF%TSPLITS_SX_YP2_ZP1(IP)
  TZCHILD_COMDATA%TSPLIT_SXP2_Y_ZP1 => TZCHILD_PROCONF%TSPLITS_SXP2_Y_ZP1(IP)
!JUAN Z_SPLITTING

  CALL CONSTRUCT_HALO1(TZCHILD_COMDATA, TZCHILD_PROCONF)
  CALL CONSTRUCT_HALO2(TZCHILD_COMDATA, TZCHILD_PROCONF)
  CALL CONSTRUCT_HALO_EXTENDED(TZCHILD_COMDATA, TZCHILD_PROCONF, JPHEXT+1)
!
  CALL CONSTRUCT_TRANS(TZCHILD_COMDATA, TZCHILD_PROCONF)
!JUAN Z_SPLITTING
  CALL CONSTRUCT_TRANSZ(TZCHILD_COMDATA, TZCHILD_PROCONF)
!JUAN Z_SPLITTING
!
  ALLOCATE(TZCHILD_COMDATA%HALO1DX)
  ALLOCATE(TZCHILD_COMDATA%HALO1DX%NSEND_WEST(NPROC))
  ALLOCATE(TZCHILD_COMDATA%HALO1DX%NSEND_EAST(NPROC))
  CALL CONSTRUCT_1DX(TZCHILD_COMDATA, TZCHILD_PROCONF)
!
  ALLOCATE(TZCHILD_COMDATA%HALO1DY)
  ALLOCATE(TZCHILD_COMDATA%HALO1DY%NSEND_SOUTH(NPROC))
  ALLOCATE(TZCHILD_COMDATA%HALO1DY%NSEND_NORTH(NPROC))
  CALL CONSTRUCT_1DY(TZCHILD_COMDATA, TZCHILD_PROCONF)
!
!        4.5   Grid nesting for the child's children (recursivity)
!
  NULLIFY(TZCHILD_PROCONF%TCHILDREN)
!  NULLIFY(TZCHILD_COMDATA%TCHILDREN)
!  NULLIFY(TZCHILD_COMDATA%TP2C_DATA)
!
  DO JMODEL = 1, JPMODELMAX
    IF( NDAD(JMODEL) .EQ. TZCHILD_PROCONF%NUMBER ) THEN
      CALL INI_CHILD(TZCHILD_PROCONF, TZCHILD_COMDATA, JMODEL)
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
!        5.    UPDATE OF THE PARENT'S STRUCTURE :
!              --------------------------------
!
!        5.1   We add the variable TZCHILD_PROCONF to the list
!              of PROCONF_ll children of the parent
!
  CALL ADD_PROCONF(TPPROCONF%TCHILDREN, TZCHILD_PROCONF)
!
!        5.2  We add the variable TZP2C_DATA to the list
!             of PARENT2CHILD_DATA_ll children of the parent
!
  CALL ADD_CHILD(TPCOMDATA%TP2C_DATA, TZP2C_DATA)
!
!        5.3  We add the variable TZCHILD_COMDATA to the list
!             of PROC_COM_DATA_ll children of the parent
!
  CALL ADD_COMDATA(TPCOMDATA%TCHILDREN,TZCHILD_COMDATA)
!
!-------------------------------------------------------------------------------
!
CONTAINS
  SUBROUTINE COARSE_TO_FINE(TZ)

    IMPLICIT NONE
    
    TYPE(ZONE_ll), DIMENSION(:) :: TZ   ! grid splitting to transform from coarse (father) resolution/grid
                                        ! to fien ( son ) resolution/grid    

    INTEGER :: J

    DO J = 1, NPROC
       !
       TZ(J)%NUMBER = TZ(J)%NUMBER

       TZ(J)%NXOR = (TZ(J)%NXOR - JPHEXT -1) * NDXRATIO_ALL(K) + JPHEXT +1 ! -/+2
       TZ(J)%NYOR = (TZ(J)%NYOR - JPHEXT -1) * NDYRATIO_ALL(K) + JPHEXT +1 ! -/+2
       TZ(J)%NXEND = (TZ(J)%NXEND - JPHEXT) * NDXRATIO_ALL(K) + JPHEXT       ! -/+1
       TZ(J)%NYEND = (TZ(J)%NYEND - JPHEXT) * NDYRATIO_ALL(K) + JPHEXT       ! -/+1
!!$
!!$       TZ(J)%NXOR = (TZ(J)%NXOR - 2) * NDXRATIO_ALL(K) + 2
!!$       TZ(J)%NYOR = (TZ(J)%NYOR - 2) * NDYRATIO_ALL(K) + 2
!!$       TZ(J)%NXEND = (TZ(J)%NXEND - 1) * NDXRATIO_ALL(K) + 1
!!$       TZ(J)%NYEND = (TZ(J)%NYEND - 1) * NDYRATIO_ALL(K) + 1

       !JUAN Z_SPLITTING
       TZ(J)%NZOR  = TZ(J)%NZOR
       TZ(J)%NZEND = TZ(J)%NZEND
       !JUAN Z_SPLITTING
       !
    ENDDO

  END SUBROUTINE COARSE_TO_FINE

      END SUBROUTINE INI_CHILD
!
END MODULE MODE_NEST_ll
