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
      MODULE MODE_CONSTRUCT_ll
!     ########################
!
!!    Purpose
!!    -------
!
!     The purpose of this module is to provide subroutines and functions for
!     for the initialization of parallel data variables
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!     None
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
!     R. Guivarch, D. Lugato    * CERFACS *
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original 01/05/98
!     Juan 19/08/2005: distinction Halo NORD/SUD & EST/WEST
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types PROCONF_ll, ZONE_ll, MODELSPLITTING_ll, CRSPD_ll,
!              PROC_COM_DATA_ll, LPROC_COM_DATA_ll, LCRSPD_ll,
!              LOCALISATION_ll
!
!     Module MODD_PARAMETERS_ll
!        JPHEXT - Horizontal External points number
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
!        JPHALO - size of the halo
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!        DIMX, DIMY, DIMZ - dimensions of the extended domain
!
!     Module MODD_DIM_ll
!        CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!        CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
!
!------------------------------------------------------------------------------
!
  USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll, ZONE_ll, MODELSPLITTING_ll, &
                                 CRSPD_ll, PROC_COM_DATA_ll, LCRSPD_ll,  &
                                 LPROC_COM_DATA_ll, LOCALISATION_ll
!
  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
!
  USE MODD_VAR_ll, ONLY        : IP, NPROC, JPHALO, TCRRT_COMDATA,       &
                                 DIMX, DIMY, DIMZ 
!
  USE MODD_DIM_ll, ONLY        : CLBCX, CLBCY
!
  USE MODD_MPIF
!
!  INCLUDE 'mpif.h'
!
  CONTAINS
!
!     #####################################
      SUBROUTINE INI_PZ( TPPROCONF, TPPZS )
!     #####################################
!
!!****  *INI_PZ* - routine to initialize the physical 2way splitting
!                  of the TPPROCONF variable
! 
!!    Purpose
!!    -------
!     the purpose of this routine is to fill the arguments of the
!     variable TPPROCONF$TSPLITS_B concerning physical subdomains
!     with a given splitting TPPZS
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
!        types PROCONF_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!        NPROC - Number of processors
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
!     R. Guivarch 01/08/98 arguments for grid-nesting
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : PROCONF_ll, ZONE_ll
!  USE MODD_VAR_ll, ONLY       : NPROC
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(PROCONF_ll), POINTER              :: TPPROCONF ! splitting data structure
  TYPE(ZONE_ll), DIMENSION(:), INTENT(IN):: TPPZS     ! Physical Zone Splitting
!
!
!*       0.2   declarations of local variables
!
  INTEGER :: J ! loop control variable
!
!-------------------------------------------------------------------------------
!
!*       1.    FILL TPPROCONF%TSPLITS_B FOR EACH J :
!              ---------------------------------------
!
  DO J = 1, NPROC
!
    TPPROCONF%TSPLITS_B(J)%NUMBER = TPPZS(J)%NUMBER
!
    TPPROCONF%TSPLITS_B(J)%NXORP  = TPPZS(J)%NXOR
    TPPROCONF%TSPLITS_B(J)%NXENDP = TPPZS(J)%NXEND
    TPPROCONF%TSPLITS_B(J)%NDIMXP = TPPZS(J)%NXEND - TPPZS(J)%NXOR + 1
!
    TPPROCONF%TSPLITS_B(J)%NYORP  = TPPZS(J)%NYOR
    TPPROCONF%TSPLITS_B(J)%NYENDP = TPPZS(J)%NYEND
    TPPROCONF%TSPLITS_B(J)%NDIMYP = TPPZS(J)%NYEND - TPPZS(J)%NYOR + 1
!
    TPPROCONF%TSPLITS_B(J)%NZORP  = TPPZS(J)%NZOR
    TPPROCONF%TSPLITS_B(J)%NZENDP = TPPZS(J)%NZEND
    TPPROCONF%TSPLITS_B(J)%NDIMZP = TPPZS(J)%NZEND - TPPZS(J)%NZOR + 1
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INI_PZ
!
!     ##############################
      SUBROUTINE INI_EZ( TPPROCONF )
!     ##############################
!
!!****  *INI_PZ* - routine to initialize the extended 2way splitting
!                  of the TPPROCONF variable
! 
!!    Purpose
!!    -------
!     the purpose of this routine is to fill the arguments of the
!     variable TPPROCONF$TSPLITS_B concerning extended subdomains
!
!!**  Method
!!    ------
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!        LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types PROCONF_ll MODELSPLITTING_ll
!
!     Module MODD_VAR_ll
!        JPHALO - size of the halo
!        NPROC - Number of processors
! 
!     Module MODD_PARAMETERS_ll
!        JPHEXT - Horizontal External points number
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
!     R. Guivarch 01/08/98 arguments for grid-nesting
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll, MODELSPLITTING_ll
!  USE MODD_VAR_ll, ONLY        : NPROC, JPHALO
   USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
!
  USE MODE_TOOLS_ll, ONLY      : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
!
!*       0.2   declarations of local variables
!
  INTEGER                          :: J  ! loop control variable
  TYPE(MODELSPLITTING_ll), POINTER :: TB ! Intermediate MODELSPLITTING_ll
                                         ! variable
!
!-------------------------------------------------------------------------------
!
!*       1.    FILL TPPROCONF%TSPLITS_B FOR EACH J :
!              ---------------------------------------
!
  DO J = 1, NPROC
!
    TB => TPPROCONF%TSPLITS_B(J)
    IF(LWEST_ll(J)) THEN
      TB%NXORE = TB%NXORP - JPHEXT
    ELSE
      TB%NXORE = TB%NXORP - JPHALO
    ENDIF
!
    IF(LSOUTH_ll(J)) THEN
      TB%NYORE = TB%NYORP - JPHEXT
    ELSE
      TB%NYORE = TB%NYORP - JPHALO
    ENDIF
!
    IF(LEAST_ll(J)) THEN
      TB%NXENDE = TB%NXENDP + JPHEXT
    ELSE
      TB%NXENDE = TB%NXENDP + JPHALO
    ENDIF
!
    IF(LNORTH_ll(J)) THEN
      TB%NYENDE = TB%NYENDP + JPHEXT
    ELSE
      TB%NYENDE = TB%NYENDP + JPHALO
    ENDIF
!
    !JUAN
    TB%NZORE  = TB%NZORP  - JPVEXT
    TB%NZENDE = TB%NZENDP + JPVEXT
    !JUAN
!
    TB%NDIMXE = TB%NXENDE - TB%NXORE + 1
    TB%NDIMYE = TB%NYENDE - TB%NYORE + 1
    TB%NDIMZE = TB%NZENDE - TB%NZORE + 1
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INI_EZ
!
!     ######################################################
      SUBROUTINE INI_CYCLIC( TPSENDHALO, TPRECVHALO, &
                             TPSENDBOUNDX, TPRECVBOUNDX, &
                             TPSENDBOUNDY, TPRECVBOUNDY, &
                             TPSENDBOUNDXY, TPRECVBOUNDXY, &
                             TPPZS ,TPEZS, KH )
!     ######################################################
!
!!****  *INI_CYCLIC* - routine to updates the correspondants'lists
!                      in case of cyclic conditions
! 
!!    Purpose
!!    -------
!     The purpose of this routine is to complete the correspondants'lists
!     for the processors on the boundaries of the domain with the informations
!     on cyclic conditions
!
!!**  Method
!!    ------
!     for the different cases of cyclic conditions and the different boundaries,
!     we list the processors and construct the intersections to determine
!     which zones to exchange
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!        ADD_ZONE, LWEST_ll, LSOUTH_ll, LEAST_ll, LNORTH_ll
! 
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types CRSPD_ll, ZONE_ll, PROC_COM_DATA_ll
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
! 
!     Module MODD_DIM_ll
!        CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!        CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
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
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll, ZONE_ll, PROC_COM_DATA_ll
!  USE MODD_VAR_ll, ONLY       : TCRRT_COMDATA, NPROC, IP, DIMZ
!  USE MODD_DIM_ll, ONLY       : CLBCX, CLBCY
!
  USE MODE_TOOLS_ll, ONLY     : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll, &
                                ADD_ZONE
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPSENDHALO ! Correspondants
                                        ! for send commmunications
  TYPE(CRSPD_ll), POINTER :: TPRECVHALO ! Correspondants
                                        ! for recv commmunications

  TYPE(CRSPD_ll), POINTER :: TPSENDBOUNDX, TPRECVBOUNDX, &
                             TPSENDBOUNDY, TPRECVBOUNDY, &
                             TPSENDBOUNDXY, TPRECVBOUNDXY
!
  TYPE(ZONE_ll), DIMENSION(:), INTENT(IN) :: TPPZS ! Physical Zone splitting
                                                   ! of the global domain
  TYPE(ZONE_ll), DIMENSION(:), INTENT(IN) :: TPEZS ! Extended Zone splitting
                                                   ! of the global domain
!
  INTEGER, INTENT(IN)                     :: KH    ! size of the external
                                                   ! halo (=JPHEXT)
!
!
!*       0.2   declarations of local variables
!
  INTEGER       :: J ! loop control variable
  INTEGER       :: IXMIN, IXMAX, IYMIN, IYMAX ! intermediate variables
                                        !for intersections
!
  TYPE(ZONE_ll) :: TZZONE_INTER !Intermediate zone variable
!
  INTEGER       :: ICURMODEL
!
!-------------------------------------------------------------------------------
!
  ICURMODEL = TCRRT_COMDATA%NUMBER
!
!*       1.    WEST BOUNDARY :
!              -------------
!
  IF(CLBCX(ICURMODEL, 1) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the WEST boundary
!
    IF(LWEST_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the EAST ones
        IF(LEAST_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the extended subdomain of processor j
!         to determine the send zone (for processor IP)
          IYMIN = MAX( TPPZS(IP)%NYOR, TPEZS(J)%NYOR)
          IYMAX = MIN( TPPZS(IP)%NYEND, TPEZS(J)%NYEND)
!
          IF(IYMIN <= IYMAX) THEN
!           if the intersection is not void we add the zone to
!           the send correspondant
            TZZONE_INTER = ZONE_ll( J, 4, &
                                 TPPZS(IP)%NXOR,&
                                 TPPZS(IP)%NXOR + KH - 1,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ )

!!$            TZZONE_INTER = ZONE_ll( J, 4, &
!!$                                 TPPZS(IP)%NXOR, IYMIN, 1, &
!!$                                 TPPZS(IP)%NXOR + KH - 1, IYMAX, DIMZ )

            CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPSENDBOUNDX, TZZONE_INTER)
            CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
          ENDIF
!
!         compute the intersection between the extended subdomain
!         of processor IP with the physical subdomain of processor j
!         to determine the receive zone (for processor IP)
          IYMIN = MAX( TPEZS(IP)%NYOR, TPPZS(J)%NYOR)
          IYMAX = MIN( TPEZS(IP)%NYEND, TPPZS(J)%NYEND)
!
          IF(IYMIN <= IYMAX) THEN
!           if the intersection is not void we add the zone to
!           the receive correspondant
            TZZONE_INTER = ZONE_ll( J, 8, &
                                 TPEZS(IP)%NXOR,&
                                 TPEZS(IP)%NXOR + KH - 1,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, 8, &
!!$                                 TPEZS(IP)%NXOR, IYMIN, 1, &
!!$                                 TPEZS(IP)%NXOR + KH - 1, IYMAX, DIMZ)

            CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPRECVBOUNDX, TZZONE_INTER)
            CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          ENDIF
!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    EAST BOUNDARY :
!              -------------
!
  IF(CLBCX(ICURMODEL, 2) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the EAST boundary
!
    IF(LEAST_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the WEST ones
        IF(LWEST_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the extended subdomain of processor j
!         to determine the send zone (for processor IP)
          IYMIN = MAX( TPPZS(IP)%NYOR, TPEZS(J)%NYOR)
          IYMAX = MIN( TPPZS(IP)%NYEND, TPEZS(J)%NYEND)
!
          IF(IYMIN <= IYMAX) THEN
!           if the intersection is not void we add the zone to
!           the send correspondant

            TZZONE_INTER = ZONE_ll( J, 8, &
                                 TPPZS(IP)%NXEND - KH + 1,&
                                 TPPZS(IP)%NXEND,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ )

!!$            TZZONE_INTER = ZONE_ll( J, 8, &
!!$                                 TPPZS(IP)%NXEND - KH + 1, IYMIN, 1, &
!!$                                 TPPZS(IP)%NXEND, IYMAX, DIMZ)

            CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPSENDBOUNDX, TZZONE_INTER)
            CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
          ENDIF
!
!         compute the intersection between the extended subdomain
!         of processor IP with the physical subdomain of processor j
!         to determine the receive zone (for processor IP)
          IYMIN = MAX( TPEZS(IP)%NYOR, TPPZS(J)%NYOR)
          IYMAX = MIN( TPEZS(IP)%NYEND, TPPZS(J)%NYEND)
!
          IF(IYMIN <= IYMAX) THEN
!           if the intersection is not void we add the zone to
!           the receive correspondant
            TZZONE_INTER = ZONE_ll( J, 4, &
                                 TPEZS(IP)%NXEND - KH + 1,&
                                 TPEZS(IP)%NXEND,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, 4, &
!!$                                 TPEZS(IP)%NXEND - KH + 1, IYMIN, 1, &
!!$                                 TPEZS(IP)%NXEND, IYMAX, DIMZ)

            CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPRECVBOUNDX, TZZONE_INTER)
            CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          ENDIF

!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.    SOUTH BOUNDARY :
!              --------------
!
  IF(CLBCY(ICURMODEL, 1) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the SOUTH boundary
!
    IF(LSOUTH_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the NORTH ones
        IF(LNORTH_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the extended subdomain of processor j
!         to determine the send zone (for processor IP)
          IXMIN = MAX( TPPZS(IP)%NXOR, TPEZS(J)%NXOR)
          IXMAX = MIN( TPPZS(IP)%NXEND, TPEZS(J)%NXEND)
!
          IF(IXMIN <= IXMAX) THEN
!           if the intersection is not void we add the zone to
!           the send correspondant
            TZZONE_INTER = ZONE_ll( J, 6, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPPZS(IP)%NYOR,&
                                 TPPZS(IP)%NYOR + KH - 1,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, 6, &
!!$                                 IXMIN, TPPZS(IP)%NYOR, 1, &
!!$                                 IXMAX, TPPZS(IP)%NYOR + KH - 1, DIMZ)

            CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPSENDBOUNDY, TZZONE_INTER)
            CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
          ENDIF
!
!         compute the intersection between the extended subdomain
!         of processor IP with the physical subdomain of processor j
!         to determine the receive zone (for processor IP)
          IXMIN = MAX( TPEZS(IP)%NXOR, TPPZS(J)%NXOR)
          IXMAX = MIN( TPEZS(IP)%NXEND, TPPZS(J)%NXEND)
!
          IF(IXMIN <= IXMAX) THEN
!           if the intersection is not void we add the zone to
!           the receive correspondant
            TZZONE_INTER = ZONE_ll( J, 2, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPEZS(IP)%NYOR,&
                                 TPEZS(IP)%NYOR + KH - 1,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, 2, &
!!$                                 IXMIN, TPEZS(IP)%NYOR, 1, &
!!$                                 IXMAX, TPEZS(IP)%NYOR + KH - 1, DIMZ)

            CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPRECVBOUNDY, TZZONE_INTER)
            CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          ENDIF
!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.    NORTH BOUNDARY :
!              --------------
!
  IF(CLBCY(ICURMODEL, 2) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the NORTH boundary
!
    IF(LNORTH_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the SOUTH ones
        IF(LSOUTH_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the extended subdomain of processor j
!         to determine the send zone (for processor IP)
          IXMIN = MAX( TPPZS(IP)%NXOR, TPEZS(J)%NXOR)
          IXMAX = MIN( TPPZS(IP)%NXEND, TPEZS(J)%NXEND)
!
          IF(IXMIN <= IXMAX) THEN
!           if the intersection is not void we add the zone to
!           the send correspondant
            TZZONE_INTER = ZONE_ll( J, 2, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPPZS(IP)%NYEND - KH + 1,&
                                 TPPZS(IP)%NYEND,&
                                 1, &
                                 DIMZ )

!!$            TZZONE_INTER = ZONE_ll( J, 2, &
!!$                                 IXMIN, TPPZS(IP)%NYEND - KH + 1, 1, &
!!$                                 IXMAX, TPPZS(IP)%NYEND, DIMZ )

            CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPSENDBOUNDY, TZZONE_INTER)
            CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
          ENDIF
!
!         compute the intersection between the extended subdomain
!         of processor IP with the physical subdomain of processor j
!         to determine the receive zone (for processor IP)
          IXMIN = MAX( TPEZS(IP)%NXOR, TPPZS(J)%NXOR)
          IXMAX = MIN( TPEZS(IP)%NXEND, TPPZS(J)%NXEND)
!
          IF(IXMIN <= IXMAX) THEN
!           if the intersection is not void we add the zone to
!           the receive correspondant
            TZZONE_INTER = ZONE_ll( J, 6, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPEZS(IP)%NYEND - KH + 1,&
                                 TPEZS(IP)%NYEND,&
                                 1, &
                                 DIMZ )

!!$            TZZONE_INTER = ZONE_ll( J, 6, &
!!$                                 IXMIN, TPEZS(IP)%NYEND - KH + 1, 1, &
!!$                                 IXMAX, TPEZS(IP)%NYEND, DIMZ )

            CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
            CALL ADD_ZONE( TPRECVBOUNDY, TZZONE_INTER)
            CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          ENDIF
!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.    SOUTH-EAST CORNER :
!              -----------------
!
  IF((CLBCX(ICURMODEL, 2).EQ. 'CYCL') &
    .AND. (CLBCY(ICURMODEL, 1).EQ. 'CYCL')) THEN
!
!    the local subdomain (processor) is on the SOUTH-EAST corner
!
    IF( LEAST_ll() .AND. LSOUTH_ll() ) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the NORTH_WEST one
        IF( LNORTH_ll(J) .AND. LWEST_ll(J) ) THEN
!
!         we add the zone to the send correspondant
          TZZONE_INTER = ZONE_ll( J, 7, &
                     TPPZS(IP)%NXEND - KH + 1,&
                     TPPZS(IP)%NXEND,&
                     TPPZS(IP)%NYOR,&
                     TPPZS(IP)%NYOR + KH - 1,&
                     1, &
                     DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 7, &
!!$                     TPPZS(IP)%NXEND - KH + 1, TPPZS(IP)%NYOR, 1, &
!!$                     TPPZS(IP)%NXEND, TPPZS(IP)%NYOR + KH - 1, DIMZ )

          CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
!
!         we add the zone to the receive correspondant
          TZZONE_INTER = ZONE_ll( J, 3, &
                     TPEZS(IP)%NXEND - KH + 1,&
                     TPEZS(IP)%NXEND,&
                     TPEZS(IP)%NYOR,&
                     TPEZS(IP)%NYOR + KH - 1,&
                     1, &
                     DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 3, &
!!$                     TPEZS(IP)%NXEND - KH + 1, TPEZS(IP)%NYOR, 1, &
!!$                     TPEZS(IP)%NXEND, TPEZS(IP)%NYOR + KH - 1, DIMZ )

          CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          EXIT
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       6.    NORTH-WEST CORNER :
!              -----------------
!
  IF((CLBCX(ICURMODEL, 1).EQ. 'CYCL') &
    .AND. (CLBCY(ICURMODEL, 2).EQ. 'CYCL')) THEN
!
!    the local subdomain (processor) is on the NORTH-WEST corner
!
    IF( LNORTH_ll() .AND. LWEST_ll() ) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the SOUTH-EAST one
        IF( LEAST_ll(J) .AND. LSOUTH_ll(J) ) THEN
!
!         we add the zone to the send correspondant
          TZZONE_INTER = ZONE_ll( J, 3, &
                     TPPZS(IP)%NXOR,&
                     TPPZS(IP)%NXOR + KH - 1,&
                     TPPZS(IP)%NYEND - KH + 1,&
                     TPPZS(IP)%NYEND,&
                     1, &
                     DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 3, &
!!$                     TPPZS(IP)%NXOR, TPPZS(IP)%NYEND - KH + 1, 1, &
!!$                     TPPZS(IP)%NXOR + KH - 1, TPPZS(IP)%NYEND, DIMZ )

          CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
!
!         we add the zone to the receive correspondant
          TZZONE_INTER = ZONE_ll( J, 7, &
                     TPEZS(IP)%NXOR,&
                     TPEZS(IP)%NXOR + KH - 1,&
                     TPEZS(IP)%NYEND - KH + 1,&
                     TPEZS(IP)%NYEND,&
                     1, &
                     DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 7, &
!!$                     TPEZS(IP)%NXOR, TPEZS(IP)%NYEND - KH + 1, 1, &
!!$                     TPEZS(IP)%NXOR + KH - 1, TPEZS(IP)%NYEND, DIMZ )

          CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          EXIT
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       7.    NORTH-EAST CORNER :
!              -----------------
!
  IF((CLBCX(ICURMODEL, 2).EQ. 'CYCL') &
    .AND. (CLBCY(ICURMODEL, 2).EQ. 'CYCL')) THEN
!
!    the local subdomain (processor) is on the NORTH-EAST corner
!
    IF( LNORTH_ll() .AND. LEAST_ll() ) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the SOUTH-WEST one
        IF( LSOUTH_ll(J) .AND. LWEST_ll(J) ) THEN
!
!         we add the zone to the send correspondant
          TZZONE_INTER = ZONE_ll( J, 9, &
            TPPZS(IP)%NXEND - KH + 1,&
            TPPZS(IP)%NXEND,&
            TPPZS(IP)%NYEND - KH + 1,&
            TPPZS(IP)%NYEND,&
            1, &
            DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 9, &
!!$            TPPZS(IP)%NXEND - KH + 1, TPPZS(IP)%NYEND - KH + 1, 1, &
!!$            TPPZS(IP)%NXEND, TPPZS(IP)%NYEND, DIMZ )

          CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
!
!         we add the zone to the receive correspondant
          TZZONE_INTER = ZONE_ll( J, 5, &
             TPEZS(IP)%NXEND - KH + 1,&
             TPEZS(IP)%NXEND,&
             TPEZS(IP)%NYEND - KH + 1,&
             TPEZS(IP)%NYEND,&
             1, &
             DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 5, &
!!$             TPEZS(IP)%NXEND - KH + 1, TPEZS(IP)%NYEND - KH + 1, 1, &
!!$             TPEZS(IP)%NXEND, TPEZS(IP)%NYEND, DIMZ )

          CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          EXIT
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       8.    SOUTH-WEST CORNER :
!              -----------------
!
  IF((CLBCX(ICURMODEL, 1).EQ. 'CYCL') &
    .AND. (CLBCY(ICURMODEL, 1).EQ. 'CYCL')) THEN
!
!    the local subdomain (processor) is on the SOUTH-WEST corner
!
    IF( LSOUTH_ll() .AND. LWEST_ll() ) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the NORTH-EAST one
        IF( LNORTH_ll(J) .AND. LEAST_ll(J) ) THEN
!
!         we add the zone to the send correspondant
          TZZONE_INTER = ZONE_ll( J, 5, &
            TPPZS(IP)%NXOR,&
            TPPZS(IP)%NXOR + KH - 1,&
            TPPZS(IP)%NYOR,&
            TPPZS(IP)%NYOR + KH - 1,&
            1, &
            DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 5, &
!!$            TPPZS(IP)%NXOR, TPPZS(IP)%NYOR, 1, &
!!$            TPPZS(IP)%NXOR + KH - 1, TPPZS(IP)%NYOR + KH - 1, DIMZ )

          CALL ADD_ZONE( TPSENDHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPSENDBOUNDXY, TZZONE_INTER)
!
!         we add the zone to the receive correspondant
          TZZONE_INTER = ZONE_ll( J, 9, &
            TPEZS(IP)%NXOR,&
            TPEZS(IP)%NXOR + KH - 1,&
            TPEZS(IP)%NYOR,&
            TPEZS(IP)%NYOR + KH - 1,&
            1, &
            DIMZ )

!!$          TZZONE_INTER = ZONE_ll( J, 9, &
!!$            TPEZS(IP)%NXOR, TPEZS(IP)%NYOR, 1, &
!!$            TPEZS(IP)%NXOR + KH - 1, TPEZS(IP)%NYOR + KH - 1, DIMZ )

          CALL ADD_ZONE( TPRECVHALO, TZZONE_INTER )
          CALL ADD_ZONE( TPRECVBOUNDXY, TZZONE_INTER)
          EXIT
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INI_CYCLIC
!
!     ##################################################
      SUBROUTINE CONSTRUCT_HALO1( TPCOMDATA, TPPROCONF )
!     ##################################################
!
!!****  *CONSTRUCT_HALO1* - routine to construct the halo1 correspondants
! 
!!    Purpose
!!    -------
!     the purpose of the routine is to fill the structured type variable
!     TPCOMDATA with informations concerning the communications of
!     halo of size 1
!
!!**  Method
!!    ------
!     we compute for the local processor,
!      - intersections between extended zones of the global domain
!        and local physical zone to find the send correspondant
!        of the local processor
!      - intersections between physical zones of the global domain
!        and local extended zone to find the receive correspondant
!        of the local processor
! 
!     we complete these correspondants in case of cyclic conditions
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!        ADD_ZONE, INTERSECTION, GLOBAL2LOCAL, EXTRACT_ZONE
!        LWEST_ll, LSOUTH_ll, LEAST_ll, LNORTH_ll
!
!     Module MODE_CONSTRUCT_ll
!        INI_CYCLIC
! 
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!
!     Module MODD_PARAMETERS_ll
!        JPHEXT - Horizontal External points number
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!
!     Module MODD_DIM_ll
!        CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!        CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
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
!     R. Guivarch 01/08/98 arguments for grid-nesting
!       ??          ??       zones intersects itself
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY  : ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
!  USE MODD_VAR_ll, ONLY        : NPROC, IP, TCRRT_COMDATA
!  USE MODD_DIM_ll, ONLY        : CLBCX, CLBCY
!
  USE MODE_TOOLS_ll, ONLY      : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll, &
                                 INTERSECTION, GLOBAL2LOCAL, ADD_ZONE,     &
                                 EXTRACT_ZONE
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data structure
  TYPE(PROCONF_ll), POINTER       :: TPPROCONF ! splitting data structure
!
!*       0.2   declarations of local variables
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZPZS ! Physical zone splitting
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZEZS ! Extended zone splitting
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZINTER ! Intermediate zone
!
  INTEGER                                  :: JI ! loop control variable

  INTEGER                                  :: ICURMODEL
  INTEGER                                  :: ISHIFTS, ISHIFTN,   &
                                              ISHIFTE, ISHIFTW
  INTEGER                                  :: ISHIFTSI, ISHIFTNI, &
                                              ISHIFTEI, ISHIFTWI
  INTEGER                                  :: IS, IE, IW ,IN
!
!-------------------------------------------------------------------------------
!
!*       1.    ALLOCATE OF THE LOCAL VARIABLES :
!              -------------------------------
!
  ALLOCATE( TZPZS(NPROC), TZEZS(NPROC), TZINTER(NPROC) )
!
!-------------------------------------------------------------------------------
!
!*       2.    EXTRACTION OF PHYSICAL AND EXTENDED 2WAY SPLITTING :
!              --------------------------------------------------
!
  CALL EXTRACT_ZONE(TPPROCONF%TSPLITS_B, TZPZS, TZEZS )
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTATION OF INTERSECTION BETWEEN LOCAL PHYSICAL ZONE
!*             AND EXTENDED SPLITTING -> SEND CORRESPONDANT :
!              --------------------------------------------
!
  CALL INTERSECTION( TZEZS, NPROC, TZPZS(IP), TZINTER )
!
  ICURMODEL = TCRRT_COMDATA%NUMBER
!
  ISHIFTS = 0
  ISHIFTW = 0
  ISHIFTN = 0
  ISHIFTE = 0
!
  IF (LSOUTH_ll()) ISHIFTS = 1
  IF (LWEST_ll()) ISHIFTW = 1
  IF (LNORTH_ll()) ISHIFTN = 1
  IF (LEAST_ll()) ISHIFTE = 1
!
  IF ((ISHIFTS.NE.0).OR.(ISHIFTW.NE.0).OR.(ISHIFTN.NE.0).OR. &
      (ISHIFTE.NE.0)) THEN
!
    DO JI = 1, NPROC
!
!     if intersection not void and intersected zone is zone itself
!
      IF ((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
        ISHIFTSI = 2
        ISHIFTWI = 2
        ISHIFTNI = 2
        ISHIFTEI = 2
        IF (LSOUTH_ll(JI)) ISHIFTSI = 1
        IF (LWEST_ll(JI)) ISHIFTWI = 1
        IF (LNORTH_ll(JI)) ISHIFTNI = 1
        IF (LEAST_ll(JI)) ISHIFTEI = 1

        IS = 0
        IN = 0
        IW = 0
        IE = 0
!
!     if intersected zone is on a border too
!
        IF ((ISHIFTS == ISHIFTSI).AND.(CLBCX(ICURMODEL, 1) /= 'CYCL')) THEN
          IS = -1
        ENDIF
!
        IF ((ISHIFTN == ISHIFTNI).AND.(CLBCX(ICURMODEL, 2) /= 'CYCL')) THEN
          IN = 1
        ENDIF
!
        IF ((ISHIFTW == ISHIFTWI).AND.(CLBCY(ICURMODEL, 1) /= 'CYCL')) THEN
          IW = -1
        ENDIF
!
        IF ((ISHIFTE == ISHIFTEI).AND.(CLBCY(ICURMODEL, 2) /= 'CYCL')) THEN
          IE = 1
        ENDIF
!
        TZINTER(JI) = ZONE_ll(&
             TZINTER(JI)%NUMBER         ,&
             TZINTER(JI)%MSSGTAG        ,&
             TZINTER(JI)%NXOR + IW      ,&
             TZINTER(JI)%NXEND + IE     ,&
             TZINTER(JI)%NYOR + IS      ,&
             TZINTER(JI)%NYEND + IN     ,&
             TZINTER(JI)%NZOR           ,&
             TZINTER(JI)%NZEND           )

!!$        TZINTER(JI) = ZONE_ll(TZINTER(JI)%NUMBER, TZINTER(JI)%MSSGTAG, &
!!$                              TZINTER(JI)%NXOR + IW, TZINTER(JI)%NYOR + IS, &
!!$                              TZINTER(JI)%NZOR, TZINTER(JI)%NXEND + IE, &
!!$                              TZINTER(JI)%NYEND + IN, TZINTER(JI)%NZEND)

!!$!
      ENDIF
!
    ENDDO
!
  ENDIF
!
  NULLIFY(TPCOMDATA%TSEND_HALO1)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = 1
      CALL ADD_ZONE( TPCOMDATA%TSEND_HALO1, TZINTER(JI) )
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4.    COMPUTATION OF INTERSECTION BETWEEN LOCAL EXTENDED ZONE
!              AND PHYSICAL SPLITTING -> RECV CORRESPONDANT :
!              --------------------------------------------
!
  CALL INTERSECTION( TZPZS, NPROC, TZEZS(IP), TZINTER )
!
  IF ((ISHIFTS.NE.0).OR.(ISHIFTW.NE.0).OR.(ISHIFTN.NE.0).OR. &
      (ISHIFTE.NE.0)) THEN
!
    DO JI = 1, NPROC
!
!     if intersection not void and intersected zone is zone itself
!
      IF ((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
        ISHIFTSI = 2
        ISHIFTWI = 2
        ISHIFTNI = 2
        ISHIFTEI = 2
        IF (LSOUTH_ll(JI)) ISHIFTSI = 1
        IF (LWEST_ll(JI)) ISHIFTWI = 1
        IF (LNORTH_ll(JI)) ISHIFTNI = 1
        IF (LEAST_ll(JI)) ISHIFTEI = 1
!
        IS = 0
        IN = 0
        IW = 0
        IE = 0
!
!     if intersected zone is on a border too
!
        IF ((ISHIFTS == ISHIFTSI).AND.(CLBCX(ICURMODEL, 1) /= 'CYCL')) THEN
          IS = -1
        ENDIF
!
        IF ((ISHIFTN == ISHIFTNI).AND.(CLBCX(ICURMODEL, 2) /= 'CYCL')) THEN
          IN = 1
        ENDIF
!
        IF ((ISHIFTW == ISHIFTWI).AND.(CLBCY(ICURMODEL, 1) /= 'CYCL')) THEN
          IW = -1
        ENDIF
!
        IF ((ISHIFTE == ISHIFTEI).AND.(CLBCY(ICURMODEL, 2) /= 'CYCL')) THEN
          IE = 1
        ENDIF
!
        TZINTER(JI) = ZONE_ll(TZINTER(JI)%NUMBER, TZINTER(JI)%MSSGTAG, &
                              TZINTER(JI)%NXOR + IW,&
                              TZINTER(JI)%NXEND + IE, &
                              TZINTER(JI)%NYOR + IS, &
                              TZINTER(JI)%NYEND + IN,&
                              TZINTER(JI)%NZOR,&
                              TZINTER(JI)%NZEND)

!!$       TZINTER(JI) = ZONE_ll(TZINTER(JI)%NUMBER, TZINTER(JI)%MSSGTAG, &
!!$                              TZINTER(JI)%NXOR + IW, TZINTER(JI)%NYOR + IS, &
!!$                              TZINTER(JI)%NZOR, TZINTER(JI)%NXEND + IE, &
!!$                              TZINTER(JI)%NYEND + IN, TZINTER(JI)%NZEND)
!!$!
      ENDIF
!
    ENDDO
!
  ENDIF
!
  NULLIFY(TPCOMDATA%TRECV_HALO1)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = 1
      CALL ADD_ZONE( TPCOMDATA%TRECV_HALO1, TZINTER(JI) )
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       5.    MODIFICATIONS IN CASE OF CYCLIC CONDITIONS :
!              ------------------------------------------
!
  NULLIFY(TPCOMDATA%TSEND_BOUNDX)
  NULLIFY(TPCOMDATA%TRECV_BOUNDX)
  NULLIFY(TPCOMDATA%TSEND_BOUNDY)
  NULLIFY(TPCOMDATA%TRECV_BOUNDY)
  NULLIFY(TPCOMDATA%TSEND_BOUNDXY)
  NULLIFY(TPCOMDATA%TRECV_BOUNDXY)
!
  CALL INI_CYCLIC( TPCOMDATA%TSEND_HALO1, &
                   TPCOMDATA%TRECV_HALO1, &
                   TPCOMDATA%TSEND_BOUNDX, &
                   TPCOMDATA%TRECV_BOUNDX, &
                   TPCOMDATA%TSEND_BOUNDY, &
                   TPCOMDATA%TRECV_BOUNDY, &
                   TPCOMDATA%TSEND_BOUNDXY, &
                   TPCOMDATA%TRECV_BOUNDXY, &
                   TZPZS ,TZEZS, JPHEXT )
!
!
!-------------------------------------------------------------------------------
!
!*       6.    SWITCH FROM GLOBAL COORDINATES TO LOCAL COORDINATES :
!              ---------------------------------------------------
!
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TSEND_HALO1)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TRECV_HALO1)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TSEND_BOUNDX)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TRECV_BOUNDX)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TSEND_BOUNDY)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TRECV_BOUNDY)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TSEND_BOUNDXY)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TRECV_BOUNDXY)
!
!
!-------------------------------------------------------------------------------
!
!*       7.    DEALLOCATION OF LOCAL VARIABLES :
!              -------------------------------
!
  DEALLOCATE( TZPZS, TZEZS, TZINTER )
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE CONSTRUCT_HALO1
!
!     ################################################
      SUBROUTINE CONSTRUCT_1DX( TPCOMDATA, TPPROCONF )
!     ################################################
!
!!**** *CONSTRUCT_1DX* - routine to compute the communication
!                       data for the 1D fields in the x direction
! 
!!    Purpose
!!    -------
!       This routine is used to compute the communication data
!     for the 1D halo updates
!
!!**  Method
!!    ------
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!        INTERSECTION, EXTRACT_ZONE, LWEST_ll, LEAST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!
!     Module MODD_PARAMETERS_ll
!        JPHEXT - Horizontal External points number
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!        DIMX - dimensions of the extended domain
!
!     Module MODD_DIM_ll
!        CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
! 
!!    Modifications
!!    -------------
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY  : ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!  USE MODD_DIM_ll, ONLY        : CLBCX
!  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
!  USE MODD_VAR_ll, ONLY        : TCRRT_COMDATA, DIMX, NPROC, IP
!
  USE MODE_TOOLS_ll, ONLY      : INTERSECTION, EXTRACT_ZONE, LEAST_ll, LWEST_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data
!
  TYPE(PROCONF_ll), POINTER       :: TPPROCONF ! proc configuration data
!
!*       0.2   declarations of local variables
!
  INTEGER                         :: J
  INTEGER                         :: ICOUNTWEST, ICOUNTEAST
  INTEGER, DIMENSION(NPROC)       :: IWEST, IEAST
!
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZPZS, TZEZS, TZBOTTOMTAB, TZINTER
  TYPE(ZONE_ll)                   :: TZBOTTOM
!
  INTEGER                         :: ICURMODEL
!
!-------------------------------------------------------------------------------
!
  ICURMODEL = TCRRT_COMDATA%NUMBER
!
!*       1.   COMPUTE WHICH POINTS HAVE TO BE RECEIVED
!             ----------------------------------------
!
!*       1.1  Extract physical and extended zones
!
  CALL EXTRACT_ZONE(TPPROCONF%TSPLITS_B, TZPZS, TZEZS)
!
!*       1.2  Build a bottom zone
!
  TZBOTTOM%NUMBER = 0
  TZBOTTOM%MSSGTAG = 0
  TZBOTTOM%NXOR = TZEZS(IP)%NXOR
  TZBOTTOM%NYOR = TZPZS(IP)%NYOR
  TZBOTTOM%NXEND = TZEZS(IP)%NXEND
  TZBOTTOM%NYEND = TZPZS(IP)%NYOR
!
!*       1.3  Intersect the bottom zone with the physical zone splitting
!
  CALL INTERSECTION(TZPZS, NPROC, TZBOTTOM, TZINTER)
  TZINTER(IP) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
!
!*       1.4  Set the information in TPCOMDATA structure
!
  DO J = 1,NPROC
    IF (TZINTER(J)%NXOR /= 0) THEN
!
!* Test if the intersection is on the west of the subdomain
!
      IF (TZINTER(J)%NXOR == TZEZS(IP)%NXOR) THEN
        TPCOMDATA%HALO1DX%NRECV_WEST = J
      ENDIF
!
!* Test if the intersection is on the east of the subdomain
!
      IF (TZINTER(J)%NXEND == TZEZS(IP)%NXEND) THEN
        TPCOMDATA%HALO1DX%NRECV_EAST = J
      ENDIF
!
    ENDIF
  ENDDO
!
!* Test whether the subdomain is on the west boundary
!
  IF (LWEST_ll()) THEN
    IF (CLBCX(ICURMODEL, 1) == 'CYCL') THEN
      TZBOTTOM%NXOR = DIMX - 2*JPHEXT + 1
      TZBOTTOM%NXEND = DIMX - JPHEXT
      CALL INTERSECTION(TZPZS, NPROC, TZBOTTOM, TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          TPCOMDATA%HALO1DX%NRECV_WEST = J
        ENDIF
      ENDDO
    ELSE
      TPCOMDATA%HALO1DX%NRECV_WEST = MPI_PROC_NULL 
    ENDIF
  ENDIF
!
!* Test whether the subdomain is on the east boundary
!
  IF (LEAST_ll()) THEN
    IF (CLBCX(ICURMODEL, 2) == 'CYCL') THEN
      TZBOTTOM%NXOR = 1 + JPHEXT
      TZBOTTOM%NXEND = 2*JPHEXT
      CALL INTERSECTION(TZPZS, NPROC, TZBOTTOM, TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          TPCOMDATA%HALO1DX%NRECV_EAST = J
        ENDIF
      ENDDO
    ELSE
      TPCOMDATA%HALO1DX%NRECV_EAST = MPI_PROC_NULL 
    ENDIF
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.   COMPUTE WHICH POINTS HAVE TO BE SENT
!             ------------------------------------
!
!*       2.1  Build bottom zone
!
  DO J=1,NPROC
    TZBOTTOMTAB(J)%NUMBER = 0
    TZBOTTOMTAB(J)%MSSGTAG = 0
    TZBOTTOMTAB(J)%NXOR = TZEZS(J)%NXOR
    TZBOTTOMTAB(J)%NYOR = TZPZS(J)%NYOR
    TZBOTTOMTAB(J)%NXEND = TZEZS(J)%NXEND
    TZBOTTOMTAB(J)%NYEND = TZPZS(J)%NYOR
  ENDDO
!
  TZBOTTOMTAB(IP) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
!
!*       2.2  Intersect the bottom zones with the zone of the current proc
!
  CALL INTERSECTION(TZBOTTOMTAB, NPROC, TZPZS(IP), TZINTER)
!
!*       2.3  Set the information in TPCOMDATA structure
!
  ICOUNTWEST = 0
  ICOUNTEAST = 0
!
  DO J = 1,NPROC
    IF (TZINTER(J)%NXOR /= 0) THEN
!
!* Test if the intersection is on the west of the subdomain
!
      IF (TZINTER(J)%NXOR == TZPZS(IP)%NXOR) THEN
        ICOUNTWEST = ICOUNTWEST + 1
        IWEST(ICOUNTWEST) = J
      ENDIF
!
!* Test if the intersection is on the east of the subdomain
!
      IF (TZINTER(J)%NXEND == TZPZS(IP)%NXEND) THEN
        ICOUNTEAST = ICOUNTEAST + 1
        IEAST(ICOUNTEAST) = J
      ENDIF
!
    ENDIF
  ENDDO
!
!* Test whether the subdomain is on the west boundary
!
  IF (LWEST_ll()) THEN
    IF (CLBCX(ICURMODEL, 1) == 'CYCL') THEN
      DO J=1, NPROC
        IF (LEAST_ll(J)) THEN
          TZBOTTOMTAB(J)%NXOR = 1 + JPHEXT
          TZBOTTOMTAB(J)%NXEND = 2*JPHEXT
        ELSE
          TZBOTTOMTAB(J) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
        ENDIF
      ENDDO
!
      CALL INTERSECTION(TZBOTTOMTAB, NPROC, TZPZS(IP), TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          ICOUNTWEST = ICOUNTWEST + 1
          IWEST(ICOUNTWEST) = J
        ENDIF
      ENDDO
!
    ELSE
!
      ICOUNTWEST = 0
!      IWEST(1) = MPI_PROC_NULL
!
    ENDIF
  ENDIF
!
!* Test whether the subdomain is on the east boundary
!
  IF (LEAST_ll()) THEN
    IF (CLBCX(ICURMODEL, 2) == 'CYCL') THEN
      DO J=1,NPROC
        IF (LWEST_ll(J)) THEN
          TZBOTTOMTAB(J)%NXOR = DIMX - 2*JPHEXT + 1
          TZBOTTOMTAB(J)%NXEND = DIMX - JPHEXT
        ELSE
          TZBOTTOMTAB(J) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
        ENDIF
      ENDDO
!
      CALL INTERSECTION(TZBOTTOMTAB, NPROC, TZPZS(IP), TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          ICOUNTEAST = ICOUNTEAST + 1
          IEAST(ICOUNTEAST) = J
        ENDIF
      ENDDO
!
    ELSE
!
      ICOUNTEAST = 0
!      IEAST(1) = MPI_PROC_NULL
!
    ENDIF
  ENDIF
!
  IF (LWEST_ll() .AND. LEAST_ll() .AND. CLBCX(ICURMODEL, 1)=='CYCL') THEN
    ICOUNTWEST = 1
    ICOUNTEAST = 1
!
    IWEST(1) = IP
    IEAST(1) = IP
    TPCOMDATA%HALO1DX%NRECV_WEST = IP
    TPCOMDATA%HALO1DX%NRECV_EAST = IP
  ENDIF
!
  TPCOMDATA%HALO1DX%NSEND_WEST(1:ICOUNTWEST) = IWEST(1:ICOUNTWEST)
  TPCOMDATA%HALO1DX%NSEND_EAST(1:ICOUNTEAST) = IEAST(1:ICOUNTEAST)
!
  TPCOMDATA%HALO1DX%NBSEND_WEST = ICOUNTWEST
  TPCOMDATA%HALO1DX%NBSEND_EAST = ICOUNTEAST
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE CONSTRUCT_1DX
!
!     ################################################
      SUBROUTINE CONSTRUCT_1DY( TPCOMDATA, TPPROCONF )
!     ################################################
!
!!**** *CONSTRUCT_1D* - routine to compute the communication
!                       data for the 1D fields
! 
!!    Purpose
!!    -------
! 
!       This routine is used to compute the communication data
!     for the 1D halo updates
! 
!!**  Method
!!    ------
! 
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!        INTERSECTION, EXTRACT_ZONE, LNORTH_ll, LSOUTH_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!
!     Module MODD_PARAMETERS_ll
!        JPHEXT - Horizontal External points number
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!        DIMY - dimensions of the extended domain
!
!     Module MODD_DIM_ll
!        CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
! 
!!    Modifications
!!    -------------
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY  : ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
!  USE MODD_DIM_ll, ONLY        : CLBCY
!  USE MODD_VAR_ll, ONLY        : TCRRT_COMDATA, NPROC, IP, DIMY
!
  USE MODE_TOOLS_ll, ONLY      : INTERSECTION, EXTRACT_ZONE, &
                                 LNORTH_ll, LSOUTH_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data
!
  TYPE(PROCONF_ll), POINTER       :: TPPROCONF ! proc configuration data
!
!*       0.2   declarations of local variables
!
  INTEGER                         :: J
  INTEGER                         :: ICOUNTSOUTH, ICOUNTNORTH
  INTEGER, DIMENSION(NPROC)       :: ISOUTH, INORTH
!
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZPZS, TZEZS, TZLEFTTAB, TZINTER
  TYPE(ZONE_ll)                   :: TZLEFT
!
  INTEGER                         :: ICURMODEL
!
!-------------------------------------------------------------------------------
!
  ICURMODEL = TCRRT_COMDATA%NUMBER
!
!*       1.   COMPUTE WHICH POINTS HAVE TO BE RECEIVED
!             ----------------------------------------
!
!*       1.1  Extract physical and extended zones
!
  CALL EXTRACT_ZONE(TPPROCONF%TSPLITS_B, TZPZS, TZEZS)
!
!*       1.2  Build a bottom zone
!
  TZLEFT%NUMBER = 0
  TZLEFT%MSSGTAG = 0
  TZLEFT%NXOR = TZPZS(IP)%NXOR
  TZLEFT%NYOR = TZEZS(IP)%NYOR
  TZLEFT%NXEND = TZPZS(IP)%NXOR
  TZLEFT%NYEND = TZEZS(IP)%NYEND
!
!*       1.3  Intersect the bottom zone with the physical zone splitting
!
  CALL INTERSECTION(TZPZS, NPROC, TZLEFT, TZINTER)
  TZINTER(IP) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
!
!*       1.4  Set the information in TPCOMDATA structure
!
  DO J = 1,NPROC
    IF (TZINTER(J)%NXOR /= 0) THEN
!
!* Test if the intersection is on the south of the subdomain
!
      IF (TZINTER(J)%NYOR == TZEZS(IP)%NYOR) THEN
        TPCOMDATA%HALO1DY%NRECV_SOUTH = J
      ENDIF
!
!* Test if the intersection is on the north of the subdomain
!
      IF (TZINTER(J)%NYEND == TZEZS(IP)%NYEND) THEN
        TPCOMDATA%HALO1DY%NRECV_NORTH = J
      ENDIF
!
    ENDIF
  ENDDO
!
!* Test whether the subdomain is on the south boundary
!
  IF (LSOUTH_ll()) THEN
    IF (CLBCY(ICURMODEL, 1) == 'CYCL') THEN
      TZLEFT%NYOR = DIMY - 2*JPHEXT + 1
      TZLEFT%NYEND = DIMY - JPHEXT
      CALL INTERSECTION(TZPZS, NPROC, TZLEFT, TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          TPCOMDATA%HALO1DY%NRECV_SOUTH = J
        ENDIF
      ENDDO
    ELSE
      TPCOMDATA%HALO1DY%NRECV_SOUTH = MPI_PROC_NULL 
    ENDIF
  ENDIF
!
!* Test whether the subdomain is on the north boundary
!
  IF (LNORTH_ll()) THEN
    IF (CLBCY(ICURMODEL, 2) == 'CYCL') THEN
      TZLEFT%NYOR = 1 + JPHEXT
      TZLEFT%NYEND = 2*JPHEXT
      CALL INTERSECTION(TZPZS, NPROC, TZLEFT, TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          TPCOMDATA%HALO1DY%NRECV_NORTH = J
        ENDIF
      ENDDO
    ELSE
      TPCOMDATA%HALO1DY%NRECV_NORTH = MPI_PROC_NULL 
    ENDIF
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.   COMPUTE WHICH POINTS HAVE TO BE SENT
!             ------------------------------------
!
!*       2.1  Build left zone
!
  DO J=1,NPROC
    TZLEFTTAB(J)%NUMBER = 0
    TZLEFTTAB(J)%MSSGTAG = 0
    TZLEFTTAB(J)%NXOR = TZPZS(J)%NXOR
    TZLEFTTAB(J)%NYOR = TZEZS(J)%NYOR
    TZLEFTTAB(J)%NXEND = TZPZS(J)%NXOR
    TZLEFTTAB(J)%NYEND = TZEZS(J)%NYEND
  ENDDO
!
!
  TZLEFTTAB(IP) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
!
!*       2.2  Intersect the bottom zones with the zone of the current proc
!
  CALL INTERSECTION(TZLEFTTAB, NPROC, TZPZS(IP), TZINTER)
!
!*       2.3  Set the information in TPCOMDATA structure
!
  ICOUNTSOUTH = 0
  ICOUNTNORTH = 0
!
  DO J = 1,NPROC
    IF (TZINTER(J)%NXOR /= 0) THEN
!
!* Test if the intersection is on the south of the subdomain
!
      IF (TZINTER(J)%NYOR == TZPZS(IP)%NYOR) THEN
        ICOUNTSOUTH = ICOUNTSOUTH + 1
        ISOUTH(ICOUNTSOUTH) = J
      ENDIF
!
!* Test if the intersection is on the north of the subdomain
!
      IF (TZINTER(J)%NYEND == TZPZS(IP)%NYEND) THEN
        ICOUNTNORTH = ICOUNTNORTH + 1
        INORTH(ICOUNTNORTH) = J
      ENDIF
!
    ENDIF
  ENDDO
!
!* Test whether the subdomain is on the south boundary
!
  IF (LSOUTH_ll()) THEN
    IF (CLBCY(ICURMODEL, 1) == 'CYCL') THEN
      DO J=1, NPROC
        IF (LNORTH_ll(J)) THEN
          TZLEFTTAB(J)%NYOR = 1 + JPHEXT
          TZLEFTTAB(J)%NYEND = 2*JPHEXT
        ELSE
          TZLEFTTAB(J) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
        ENDIF
      ENDDO
!
      CALL INTERSECTION(TZLEFTTAB, NPROC, TZPZS(IP), TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          ICOUNTSOUTH = ICOUNTSOUTH + 1
          ISOUTH(ICOUNTSOUTH) = J
        ENDIF
      ENDDO
!
    ELSE
!
      ICOUNTSOUTH = 0
!      ISOUTH(1) = MPI_PROC_NULL
!
    ENDIF
  ENDIF
!
!* Test whether the subdomain is on the north boundary
!
  IF (LNORTH_ll()) THEN
    IF (CLBCY(ICURMODEL, 2) == 'CYCL') THEN
      DO J=1,NPROC
        IF (LSOUTH_ll(J)) THEN
          TZLEFTTAB(J)%NYOR = DIMY - 2*JPHEXT + 1
          TZLEFTTAB(J)%NYEND = DIMY - JPHEXT
        ELSE
          TZLEFTTAB(J) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
        ENDIF
      ENDDO
!
      CALL INTERSECTION(TZLEFTTAB, NPROC, TZPZS(IP), TZINTER)
      DO J=1, NPROC
        IF (TZINTER(J)%NXOR /= 0) THEN
          ICOUNTNORTH = ICOUNTNORTH + 1
          INORTH(ICOUNTNORTH) = J
        ENDIF
      ENDDO
!
    ELSE
!
      ICOUNTNORTH = 0
!      INORTH(1) = MPI_PROC_NULL
!
    ENDIF
  ENDIF
!
  IF (LNORTH_ll() .AND. LSOUTH_ll() .AND. CLBCY(ICURMODEL, 1)=='CYCL') THEN
!
    ICOUNTNORTH = 1
    ICOUNTSOUTH = 1
!
    ISOUTH(1) = IP
    INORTH(1) = IP
    TPCOMDATA%HALO1DY%NRECV_SOUTH = IP
    TPCOMDATA%HALO1DY%NRECV_NORTH = IP
!
  ENDIF
!
  TPCOMDATA%HALO1DY%NSEND_SOUTH(1:ICOUNTSOUTH) = ISOUTH(1:ICOUNTSOUTH)
  TPCOMDATA%HALO1DY%NSEND_NORTH(1:ICOUNTNORTH) = INORTH(1:ICOUNTNORTH)
!
  TPCOMDATA%HALO1DY%NBSEND_SOUTH = ICOUNTSOUTH
  TPCOMDATA%HALO1DY%NBSEND_NORTH = ICOUNTNORTH
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE CONSTRUCT_1DY
!
!     #################################################
      SUBROUTINE COMPUTE_CRSPD_MAX( KMAXSIZE, TPINTER )
!     #################################################
!
!
!!**** *COMPUTE_CRSPD_MAX* - routine to compute the maximum size
!                            of the zones included in a list of zones
!
!!    Purpose
!!    -------
!     Computes recursively the maximum size of the zones included
!     in the TZINTER list of zones
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
!        types CRSPD_ll
!
!     Module MODD_VAR_ll
!        DIMZ - dimensions of the extended domain
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
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll
!  USE MODD_VAR_ll, ONLY       : DIMZ
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER :: TPINTER
  INTEGER                 :: KMAXSIZE
!
!*       0.2   declarations of local variables
!
  TYPE(CRSPD_ll), POINTER :: TZINTER
  INTEGER                 :: ISIZE
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTES THE MAXIMUM SIZE OF THE ZONES
!
  TZINTER => TPINTER
!
  DO WHILE ( ASSOCIATED(TZINTER) )
    ISIZE = (TZINTER%TELT%NYEND - TZINTER%TELT%NYOR + 1) * &
            (TZINTER%TELT%NXEND - TZINTER%TELT%NXOR + 1) * DIMZ
!
    IF(ISIZE .GT. KMAXSIZE) KMAXSIZE = ISIZE
!
    TZINTER => TZINTER%TNEXT
!
  ENDDO
!
  NULLIFY(TZINTER)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COMPUTE_CRSPD_MAX
!
!     ############################################################
      RECURSIVE SUBROUTINE COMPUTE_HALO_MAX( KMAXSIZE, TPCOMDATA )
!     ############################################################
!
!!****  *COMPUTE_HALO_MAX* - routine to compute the maximum size
!                            of the HALO zones of all models
!
!!    Purpose
!!    -------
!     The purpose of this routine is to compute recursively the maximum
!     size of the HALO zones fot the local processor
!
!!**  Method
!!    ------
!     we list all the send and receive correspondants and computes for
!     each element the number of points for the model described by
!     TPCOMDATA and switch to the children model (for grid-nesting)
! 
!!    ExternaL
!!    --------
!
!     Module MODE_CONSTRUCT_ll
!        COMPUTE_CRSPD_MAX
!
!!    Implicit Arguments
!!    ------------------
! 
!     Module MODD_STRUCTURE_ll
!        types PROC_COM_DATA_ll, LPROC_COM_DATA_ll
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
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : PROC_COM_DATA_ll, LPROC_COM_DATA_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data structure
!
  INTEGER, INTENT(INOUT)            :: KMAXSIZE ! size max !JUAN 
!
!*       0.2   declarations of local variables
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZCHILD ! (for grid nesting)
  TYPE(PROC_COM_DATA_ll), POINTER  :: TZCOMDATA ! communications data structure
!
  INTEGER                          :: ISIZE ! Intermediate Size variable
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE KMAXSIZE FOR CURRENT MODEL :
!              ----------------------------------
!
!*       1.1   in the send correspondant (halo1)
!
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_HALO1)
!
!*       1.2   in the receive correspondant (halo1)
!
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_HALO1)
!
!*       1.3   in the send correspondant (halo2)
!
!  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_HALO2)
!
!*       1.4   in the receive correspondant (halo2)
!
!  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_HALO2)
!
!-------------------------------------------------------------------------------
!
!*       2.     GRID NESTING :
!               ------------
!
  TZCHILD => TPCOMDATA%TCHILDREN
!
  DO WHILE(ASSOCIATED(TZCHILD))
!
!    describe all the children and call recursively COMPUTE_HALO_MAX
!
    TZCOMDATA => TZCHILD%TELT
    CALL COMPUTE_HALO_MAX(KMAXSIZE, TZCOMDATA)
    TZCHILD => TZCHILD%TNEXT
!
  ENDDO
!
  NULLIFY(TZCHILD)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COMPUTE_HALO_MAX
!
!     #############################################################
      RECURSIVE SUBROUTINE COMPUTE_TRANS_MAX( KMAXSIZE, TPCOMDATA )
!     #############################################################
!
!!****  *COMPUTE_HALO1_MAX* - routine to compute the maximum size
!                             of the HALO1 zones of all models
!
!!    Purpose
!!    -------
!     The purpose of this routine is to compute recursively the maximum
!     size of the TRANS zones fot the local processor
!
!!**  Method
!!    ------
!     we list all the send and receive correspondants and computes for
!     each element the number of points for the model described by
!     TPCOMDATA and switch to the children model (for grid-nesting)
! 
!!    External
!!    --------
! 
!     Module MODE_CONSTRUCT_ll
!        COMPUTE_CRSPD_MAX
!
!!    Implicit ArgumentS
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types PROC_COM_DATA_ll, LPROC_COM_DATA_ll
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
!     R. Guivarch 01/08/98  argument TPCOMDATA instead of TCCRT_COMDATA
!     P. Kloos    20/10/98  use subroutine COMPUTE_CRSPD_MAX
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : PROC_COM_DATA_ll, LPROC_COM_DATA_ll
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data structure
!
  INTEGER, INTENT(INOUT)            :: KMAXSIZE ! size max !JUAN 
!
!*       0.2   declarations of local variables
!
  TYPE(LPROC_COM_DATA_ll), POINTER :: TZCHILD ! (for grid nesting)
  TYPE(PROC_COM_DATA_ll), POINTER  :: TZCOMDATA ! communications data structure
!
  INTEGER                          :: ISIZE ! Intermediate Size variable
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE KMAXSIZE FOR CURRENT MODEL :
!              ----------------------------------
!
!*       1.1   in the send_trans_bx correspondant
!
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_TRANS_BX)
!
!*       1.2   in the recv_trans_bx correspondant
!
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_TRANS_BX)
!
!*       1.3   in the send_trans_xy correspondant
!
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_TRANS_XY)
!
!*       1.4   in the recv_trans_xy correspondant
!
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_TRANS_XY)
!
!JUAN ZSPLITTING
!
!*       2.1   in the correspondant
!
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_SXP2_YP1_Z_SX_YP2_ZP1) !JUAN BUG? 
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_SXP2_YP1_Z_SX_YP2_ZP1) !JUAN BUG? 
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_SX_YP2_ZP1_SXP2_Y_ZP1)
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_SX_YP2_ZP1_SXP2_Y_ZP1)
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_SXP2_Y_ZP1_SXP2_YP1_Z)
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_SXP2_Y_ZP1_SXP2_YP1_Z)
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_B_SX_YP2_ZP1)
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_B_SX_YP2_ZP1)
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TSEND_SXP2_Y_ZP1_B)
  CALL COMPUTE_CRSPD_MAX(KMAXSIZE,TPCOMDATA%TRECV_SXP2_Y_ZP1_B)
!JUAN ZSPLITTING
!
!-------------------------------------------------------------------------------
!
!*       2.     GRID NESTING :
!               ------------
!
  TZCHILD => TPCOMDATA%TCHILDREN
!
  DO WHILE (ASSOCIATED(TZCHILD))
!
!    describe all the children and call recursively COMPUTE_HALO_MAX
!
    TZCOMDATA => TZCHILD%TELT
    CALL COMPUTE_TRANS_MAX(KMAXSIZE, TZCOMDATA)
    TZCHILD => TZCHILD%TNEXT
!
  ENDDO
!
  NULLIFY(TZCHILD)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE COMPUTE_TRANS_MAX
!
!     #################################
      SUBROUTINE INI_TRANS( TPPROCONF )
!     #################################
!
!!****  *INI_TRANS* - routine to initialize the transposition splitting
!                     of the global domain
!
!!    Purpose
!!    -------
!     the purpose of this routine is to split the global domain
!     in x-slice splitting and y-slice splitting ; these data are put
!     in the variable TPPROCONF
!
!!**  Method
!!    ------
!     for each direction, we divise the dimension by the number of processors.
!     each processor receives the result of this division and the rest is
!     ditributed between the first processors ( one line (colonne) per
!     processor)
!
!!    External
!!    --------
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types PROCONF_ll, MODELSPLITTING_ll
!
!     Module MODD_PARAMETERS_ll
!        JPHEXT - Horizontal External points number
!
!     Module MODD_VAR_ll
!        NPROC - Number of processors
!        DIMX, DIMY - dimensions of the extended domain
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
!     R. Guivarch 01/08/98 arguments for grid-nesting
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll, MODELSPLITTING_ll
!  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
!  USE MODD_VAR_ll, ONLY        : NPROC, DIMX, DIMY
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(PROCONF_ll), POINTER :: TPPROCONF ! splitting data structure
!
!*       0.2   declarations of local variables
!
  INTEGER                                        :: IELEM, ICHUNK, IRESTE
!
  INTEGER                                        :: J ! loop control variable
!
  INTEGER                                        :: IDEB
!
  TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TX, TY ! Intermediate
                                                 ! MODELSPLITTING_ll variables
!
!-------------------------------------------------------------------------------
!
!*       1.    TX, TY POINT TO THE TRANSPOSITION SPLITTING IN TPPROCONF :
!              ------------------------------------------------------------
!
  TX => TPPROCONF%TSPLITS_X
  TY => TPPROCONF%TSPLITS_Y
!
!-------------------------------------------------------------------------------
!
!*       2.    KNOWN VALUES :
!              --------------
  DO J = 1, NPROC
    TX(J)%NUMBER = J
    TY(J)%NUMBER = J
!
    TX(J)%NXORP = 1 + JPHEXT
    TX(J)%NXENDP = DIMX - JPHEXT
!
    TY(J)%NYORP = 1
    TY(J)%NYENDP = DIMY
!
    TX(J)%NXORE = 0
    TX(J)%NYORE = 0
    TX(J)%NXENDE = 0
    TX(J)%NYENDE = 0
!
    TX(J)%NDIMXE = 0
    TX(J)%NDIMYE = 0
!
    TY(J)%NXORE = 0
    TY(J)%NYORE = 0
    TY(J)%NXENDE = 0
    TY(J)%NYENDE = 0
!
    TY(J)%NDIMXE = 0
    TY(J)%NDIMYE = 0
!
  ENDDO
!
!-------------------------------------------------------------------------------
!
!        3.    X-SLICES SPLITTING :
!              ------------------
!
  IELEM=DIMY-2*JPHEXT
!
! result of the integer division
  ICHUNK=IELEM/NPROC
! rest
  IRESTE=IELEM-ICHUNK*NPROC
!
  TX(1)%NYORP = 1
!
  TX(1)%NYENDP = ICHUNK
  IF(IRESTE.NE.0) TX(1)%NYENDP = TX(1)%NYENDP + 1
!
  DO J = 2 , IRESTE
    TX(J)%NYORP = TX(J-1)%NYENDP + 1
    TX(J)%NYENDP = TX(J)%NYORP + ICHUNK
  ENDDO
!
  IF(IRESTE.EQ.0) THEN
    IDEB = 2
  ELSE
    IDEB = IRESTE + 1
  ENDIF
!
  DO J = IDEB, NPROC
    TX(J)%NYORP = TX(J-1)%NYENDP + 1
    TX(J)%NYENDP = TX(J)%NYORP + ICHUNK - 1
  ENDDO
!
  TX(:)%NYORP = TX(:)%NYORP + JPHEXT
  TX(:)%NYENDP = TX(:)%NYENDP + JPHEXT
  TX(:)%NDIMXP = DIMX - 2*JPHEXT
  TX(:)%NDIMYP = TX(:)%NYENDP - TX(:)%NYORP + 1
!
!-------------------------------------------------------------------------------
!
!        4.    Y-SLICES SPLITTING :
!              ------------------
!
  IELEM=DIMX
  ICHUNK=IELEM/NPROC
  IRESTE=IELEM-ICHUNK*NPROC
  TY(1)%NXORP = 1
!
  TY(1)%NXENDP = ICHUNK
  IF(IRESTE.NE.0)  TY(1)%NXENDP = TY(1)%NXENDP + 1
!
  DO J = 2 , IRESTE
    TY(J)%NXORP = TY(J-1)%NXENDP + 1
    TY(J)%NXENDP = TY(J)%NXORP + ICHUNK
  ENDDO
!
  IF(IRESTE.EQ.0) THEN
    IDEB = 2
  ELSE
    IDEB = IRESTE + 1
  ENDIF
!
  DO J = IDEB , NPROC
    TY(J)%NXORP = TY(J-1)%NXENDP + 1
    TY(J)%NXENDP = TY(J)%NXORP + ICHUNK - 1
  ENDDO
!
  TY(:)%NDIMXP = TY(:)%NXENDP - TY(:)%NXORP + 1
  TY(:)%NDIMYP = DIMY
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INI_TRANS
!
!     ##################################################
      SUBROUTINE CONSTRUCT_TRANS( TPCOMDATA, TPPROCONF )
!     ##################################################
!
!!****  *CONSTRUCT_TRANS* - routine to construct
!                           the transposition correspondants
!
!!    Purpose
!!    -------
!     the purpose of the routine is to fill the structured type variable
!     TPCOMDATA with informations concerning the communications
!     during transposition
!
!!**  Method
!!    ------
!     we compute for the local processor,
!      - intersections between zones of the global domain in x-slice splitting
!        and local zone in 2way splitting to find the send correspondant
!        of the local processor for transposition 2way / x-slices
! 
!      - intersections between physical zones of the global domain
!        in 2way splitting and local zone in x-slices splitting
!        to find the receive correspondant of the local processor
!        for transposition 2way / x-slices
! 
!      - intersections between zones of the global domain in y-slice splitting
!        and local zone in x-slices splitting to find the send correspondant
!        of the local processor for transposition x-slices / y-slices
! 
!      - intersections between zones of the global domain
!        in x-slices splitting and local zone in y-slices splitting
!        to find the receive correspondant of the local processor
!        for transposition x-slices / y-slices
! 
!!    External
!!    --------
!     Module MODE_TOOLS_ll
!        ADD_ZONE, INTERSECTION, GLOBAL2LOCAL, EXTRACT_ZONE, G2LX
!
!     Module MODE_CONSTRUCT_ll
!        RESIZE_TRANS
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
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
!     R. Guivarch 01/08/98 arguments for grid-nesting
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
!  USE MODD_VAR_ll, ONLY       : IP, NPROC
!
  USE MODE_TOOLS_ll, ONLY     : INTERSECTION, G2LX, GLOBAL2LOCAL, ADD_ZONE, &
                                EXTRACT_ZONE
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data
  TYPE(PROCONF_ll), POINTER       :: TPPROCONF ! splitting data
!
!*       0.2   declarations of local variables
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZPZS ! 2way Splitting
                                                    ! of the global domain
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZTRANSXZS ! x-slices Splitting
                                                         ! of the global domain
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZTRANSYZS ! y-slices Splitting
                                                         ! of the global domain
!
  TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZINTER ! Intermediate zone
!
  INTEGER                                  :: JI ! loop control variable
!
!-------------------------------------------------------------------------------
!
!*       1.    ALLOCATE OF THE LOCAL VARIABLES :
!              -------------------------------
!
  ALLOCATE( TZPZS(NPROC), TZTRANSXZS(NPROC), &
            TZTRANSYZS(NPROC), TZINTER(NPROC) )
!
!-------------------------------------------------------------------------------
!
!*       2.    CONSTRUCTION OF TRANSPOSITION 2WAY -> X-SLICES CORRESPONDANTS :
!              -------------------------------------------------------------
!
!*       2.1   extraction of physical 2way splitting
!              from TPPROCONF structure
!
  CALL EXTRACT_ZONE(TPPROCONF%TSPLITS_B, TZPZS, TZINTER)
!
!*       2.2   extraction of x-slices splitting
!              from TPPROCONF structure
!
  CALL EXTRACT_ZONE(TPPROCONF%TSPLITS_X, TZTRANSXZS, TZINTER)
!
!*       2.3   computation of intersection between local 2way zone
!              and x-slices splitting -> send correspondants
!
  CALL INTERSECTION( TZTRANSXZS, NPROC, TZPZS(IP), TZINTER )
!
  NULLIFY(TPCOMDATA%TSEND_TRANS_BX)
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 100
      CALL ADD_ZONE( TPCOMDATA%TSEND_TRANS_BX, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       2.4   computation of intersection between local x-slices zone
!              and 2way splitting -> recv correspondants
!
  CALL INTERSECTION( TZPZS, NPROC, TZTRANSXZS(IP), TZINTER )
!
  NULLIFY(TPCOMDATA%TRECV_TRANS_BX)
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 100
      CALL ADD_ZONE( TPCOMDATA%TRECV_TRANS_BX, TZINTER(JI) )
    ENDIF
  ENDDO
!
!
!*       2.5   Switch to local coordinates
!
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TSEND_TRANS_BX)
  CALL G2LX(TPPROCONF%TSPLITS_X(IP),TPCOMDATA%TRECV_TRANS_BX)
!

  CALL RESIZE_TRANS(TPPROCONF)
!
!-------------------------------------------------------------------------------
!
!*       3.    CONSTRUCTION OF
!              TRANSPOSITION X-SLICES -> Y-SLICES CORRESPONDANTS :
!              -------------------------------------------------
!
!*       3.1   extraction of x-slices splitting
!
  CALL EXTRACT_ZONE(TPPROCONF%TSPLITS_X, TZTRANSXZS, TZINTER)
!
!*       3.2   extraction of y-slices splitting
!
  CALL EXTRACT_ZONE(TPPROCONF%TSPLITS_Y, TZTRANSYZS, TZINTER)
!
!*       3.3   computation of intersection between local x-slices zone
!              and y-slices splitting -> send correspondants
!
  CALL INTERSECTION( TZTRANSYZS, NPROC, TZTRANSXZS(IP), TZINTER )
!
  NULLIFY(TPCOMDATA%TSEND_TRANS_XY)
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 1
      CALL ADD_ZONE( TPCOMDATA%TSEND_TRANS_XY, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       3.4   computation of intersection between local y-slices zone
!              and x-slices splitting -> recv correspondants
!
  CALL INTERSECTION( TZTRANSXZS, NPROC, TZTRANSYZS(IP), TZINTER )
!
  NULLIFY(TPCOMDATA%TRECV_TRANS_XY)
  DO JI = 1, NPROC
    IF(TZINTER(JI)%NUMBER.NE.0) THEN
      TZINTER(JI)%MSSGTAG = 1
      CALL ADD_ZONE( TPCOMDATA%TRECV_TRANS_XY, TZINTER(JI) )
    ENDIF
  ENDDO
!
!
!*       3.5   Switch to local coordinates
!
  CALL G2LX(TPPROCONF%TSPLITS_X(IP),TPCOMDATA%TSEND_TRANS_XY)
  CALL G2LX(TPPROCONF%TSPLITS_Y(IP),TPCOMDATA%TRECV_TRANS_XY)
!
!
!-------------------------------------------------------------------------------
!
!*       4.    DEALLOCATION OF LOCAL VARIABLES :
!              -------------------------------
!
  DEALLOCATE(TZPZS, TZTRANSXZS, TZTRANSYZS, TZINTER)
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE CONSTRUCT_TRANS
!
!     ####################################
      SUBROUTINE RESIZE_TRANS( TPPROCONF )
!     ####################################
!
!!****  *RESIZE_TRANS* - routine to resize the x-slices domains
!                        for transpositions
!
!!    Purpose
!!    -------
!     the purpose of this routine is to resize the x-slices domains
!     for the FFT
!
!!**  Method
!!    ------
!     all x-slices'zones are increased by 2*JPHEXT in x-direction
!     the last x-slices'zone is increased by 2*JPHEXT in y-direction
!
!!    External
!!    --------
!
!     Module MODD_STRUCTURE_ll
!        types PROCONF_ll
!
!     Module MODD_PARAMETERS_ll
!        JPHEXT - Horizontal External points number
!
!     Module MODD_VAR_ll
!        NPROC - Number of processors
!
!!    Reference
!!    ---------
!!
!!    Author
!!    ------
!!    R. Guivarch
!!
!!    Modifications
!!    -------------
!!    Original 01/05/98
!!    R. Guivarch 01/08/98 arguments for grid-nesting
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll
!  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
!  USE MODD_VAR_ll, ONLY        : NPROC
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    RESIZE X-SLICE DOMAINS :
!              ----------------------
!
  TPPROCONF%TSPLITS_X(:)%NXORP = &
    TPPROCONF%TSPLITS_X(:)%NXORP - JPHEXT
!
  TPPROCONF%TSPLITS_X(:)%NYORP = &
    TPPROCONF%TSPLITS_X(:)%NYORP - JPHEXT
!
  TPPROCONF%TSPLITS_X(:)%NXENDP = &
    TPPROCONF%TSPLITS_X(:)%NXENDP + JPHEXT
!
  TPPROCONF%TSPLITS_X(:)%NYENDP = &
    TPPROCONF%TSPLITS_X(:)%NYENDP - JPHEXT
!
  TPPROCONF%TSPLITS_X(:)%NDIMXP = &
    TPPROCONF%TSPLITS_X(:)%NDIMXP + 2*JPHEXT
!
  TPPROCONF%TSPLITS_X(NPROC)%NYENDP = &
    TPPROCONF%TSPLITS_X(NPROC)%NYENDP + 2*JPHEXT
!
  TPPROCONF%TSPLITS_X(NPROC)%NDIMYP = &
    TPPROCONF%TSPLITS_X(NPROC)%NDIMYP + 2*JPHEXT
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE RESIZE_TRANS
!
!     ######################################
      SUBROUTINE INI_BOUNDARIES( TPPROCONF )
        !     ######################################
        !
        !!****  *INI_BOUNDARIES* - routine to compute the localisation variable TBOUND
        !
        !!    Purpose
        !!    -------
        !     the purpose of this routine is to fill the localisation variable TBOUND
        !     of TPPROCONF
        !
        !!**  Method
        !!    ------
        !     for each subdomain, in function of its origin point and end point
        !     we determine if the subdomain is located on a boundarie and fill
        !     the variable TBOUND
        !
        !!    Implicit Arguments
        !!    ------------------
        !
        !     Module MODD_STRUCTURE_ll
        !        types PROCONF_ll, MODELSPLITTING_ll
        !
        !     Module MODD_PARAMETERS_ll
        !        JPHEXT - Horizontal External points number
        !
        !     Module MODD_VAR_ll
        !        NPROC - Number of processors
        !        DIMX, DIMY - dimensions of the extended domain
        !
        !!    Reference
        !!    ---------
        !!
        !!    Author
        !!    ------
        !!    R. Guivarch
        !!
        !!    Modifications
        !!    -------------
        !!    Original 01/05/98
        !!    R. Guivarch 01/08/98 arguments for grid-nesting
        !!
        !-------------------------------------------------------------------------------
        !
        !*       0.    DECLARATIONS
        !
        !  USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll, MODELSPLITTING_ll, LOCALISATION_ll
        !  USE MODD_PARAMETERS_ll, ONLY : JPHEXT
        !  USE MODD_VAR_ll, ONLY        : NPROC, DIMX, DIMY
        !
        IMPLICIT NONE
        !
        !
        !*       0.1   declarations of arguments
        !
        TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
        !
        !*       0.2   declarations of local variables
        !
        TYPE(MODELSPLITTING_ll) :: TZSPLIT ! Intermediate MODELSPLITTING_ll variable
        !
        INTEGER                 :: J ! loop control variable
        !
        !-------------------------------------------------------------------------------
        !
        !*       1.    FILL TBOUND OF EACH PROCESSOR :
        !              ------------------------------
        !
        DO J = 1, NPROC
           !
           TZSPLIT = TPPROCONF%TSPLITS_B(J)
           !
           TPPROCONF%TBOUND(J) = &
                !JUAN
                LOCALISATION_ll( .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.)
           !JUAN
           !
           IF( TZSPLIT%NXORP .EQ. (1+JPHEXT)) &
                TPPROCONF%TBOUND(J)%WEST = .TRUE.
           !
           IF( TZSPLIT%NXENDP  .EQ. (DIMX-JPHEXT)) &
                TPPROCONF%TBOUND(J)%EAST = .TRUE.
           !
           IF( TZSPLIT%NYORP .EQ. (1+JPHEXT)) &
                TPPROCONF%TBOUND(J)%SOUTH = .TRUE.
           !
           IF( TZSPLIT%NYENDP  .EQ. (DIMY-JPHEXT)) &
                TPPROCONF%TBOUND(J)%NORTH = .TRUE.
           !
        ENDDO
        !
        !-------------------------------------------------------------------------------
        !
      END SUBROUTINE INI_BOUNDARIES
!
!     ##################################################
      SUBROUTINE CONSTRUCT_HALO2( TPCOMDATA, TPPROCONF )
!     ##################################################
!
!!****  *CONSTRUCT_HALO2l* - routine to set the list of zones for
!                            second layer halo
!!    Purpose
!!    -------
!       The purpose of this routine is to set the list of zones
!     that will be exchanged during the update of the second
!     layer of the halo.
!
!!**  Method
!!    ------
!       Four arrays of type ZONE_ll are used to model the
!     second layer of the halo that will be received during the update.
!     They are WESTEXTHALO2(:), EASTEXTHALO2(:), SOUTHEXTHALO2(:),
!     NORTHEXTHALO2(:). The subscripts correspond to the subdomain number.
!       To set the zones that will be received during the update,
!     WESTEXTHALO2(IP), EASTEXTHALO2(IP), SOUTHEXTHALO2(IP), and
!     NORTHEXTHALO2(IP) are intersected with the Physical Zone
!     Splitting of the domain, where IP is the number of the current
!     processor/subdomain (see MODE_PARALLEL for the definition of
!     the Physical Zone Splitting).
!       To set the zones that will be sent during the update,
!     WESTEXTHALO2(:), EASTEXTHALO2(:), SOUTHEXTHALO2(:), and
!     NORTHEXTHALO2(:) are intersected with the physical zone
!     of the current processor/subdomain.
!
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!        ADD_ZONE, INTERSECTION, GLOBAL2LOCAL, EXTRACT_ZONE
!
!     Module MODE_CONSTRUCT_ll
!        INI_CYCLIC2
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types PROC_COM_DATA_ll, PROCONF_ll, ZONE_ll
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
!
!!    Reference
!!    ---------
!!
!!    Authors
!!    -------
!
!     R. Guivarch               * CERFACS - ENSEEIHT *
!     Ph. Kloos                 * CERFACS - CNRM *
!!
!!    Modifications
!!    -------------
!     Original 25/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : PROC_COM_DATA_ll, PROCONF_ll, ZONE_ll
!  USE MODD_VAR_ll, ONLY       : NPROC, IP
!
  USE MODE_TOOLS_ll, ONLY : ADD_ZONE, INTERSECTION, GLOBAL2LOCAL, EXTRACT_ZONE
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data
  TYPE(PROCONF_ll), POINTER       :: TPPROCONF ! splitting data
!
!*       0.2   declarations of local variables
!
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZPZS ! Physical zone splitting
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZEZS ! Extended zone splitting
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZWESTEXTHALO2  ! west second layer of
                                                     ! the external halo
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZSOUTHEXTHALO2 ! south second layer of
                                                     ! the external halo
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZNORTHEXTHALO2 ! north second layer of
                                                     ! the external halo
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZEASTEXTHALO2  ! east second layer of
                                                     ! the external halo
!
  TYPE(ZONE_ll), DIMENSION(NPROC) :: TZINTER ! Intermediate zone
!
  INTEGER                         :: ISENDSOUTHNORTH, ISENDEASTWEST, &
                                     ISENDNORTHSOUTH, ISENDWESTEAST
!
  INTEGER                         :: JI ! loop control variable
!
!-------------------------------------------------------------------------------
!
!*       1.    SET MPI MESSAGE TAGS :
!              --------------------
!
  ISENDSOUTHNORTH = 2
  ISENDEASTWEST   = 4
  ISENDNORTHSOUTH = 6
  ISENDWESTEAST   = 8
!
!-------------------------------------------------------------------------------
!
!*       1.    EXTRACTION OF PHYSICAL AND EXTENDED 2WAY SPLITTING :
!              --------------------------------------------------
!
  CALL EXTRACT_ZONE( TPPROCONF%TSPLITS_B, TZPZS, TZEZS )
!
!-------------------------------------------------------------------------------
!
!*       2.    CONSTRUCTION  OF EXTERNAL HALOS :
!              ---------------------------------
!
  CALL SET_EXT_HALO2(TZPZS, TZWESTEXTHALO2, TZSOUTHEXTHALO2, TZNORTHEXTHALO2, &
                     TZEASTEXTHALO2)
!
!-------------------------------------------------------------------------------
!
!*       3.    CONSTRUCTION OF RECEIVE CORRESPONDANT :
!              ---------------------------------------
!
  NULLIFY(TPCOMDATA%TRECV_HALO2)
!
!*       3.1   Find out the zones to be received from the south
!
  CALL INTERSECTION(TZPZS, NPROC, TZSOUTHEXTHALO2(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDSOUTHNORTH
      CALL ADD_ZONE( TPCOMDATA%TRECV_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       3.2   Find out the zones to be received from the east
!
  CALL INTERSECTION(TZPZS, NPROC, TZEASTEXTHALO2(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDEASTWEST
      CALL ADD_ZONE( TPCOMDATA%TRECV_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       3.3   Find out the zones to be received from the north
!
  CALL INTERSECTION(TZPZS, NPROC, TZNORTHEXTHALO2(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDNORTHSOUTH
      CALL ADD_ZONE( TPCOMDATA%TRECV_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       3.4   Find out the zones to be received from the west
!
  CALL INTERSECTION(TZPZS, NPROC, TZWESTEXTHALO2(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDWESTEAST
      CALL ADD_ZONE( TPCOMDATA%TRECV_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4.    CONSTRUCTION OF SEND CORRESPONDANT :
!              ------------------------------------
!
  NULLIFY(TPCOMDATA%TSEND_HALO2)
!
!*       4.1   Find out the zones to be sent to the north
!
  CALL INTERSECTION(TZSOUTHEXTHALO2, NPROC, TZPZS(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDSOUTHNORTH
      CALL ADD_ZONE( TPCOMDATA%TSEND_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       4.2   Find out the zones to be sent to the west
!
  CALL INTERSECTION(TZEASTEXTHALO2, NPROC, TZPZS(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDEASTWEST
      CALL ADD_ZONE( TPCOMDATA%TSEND_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       4.3   Find out the zones to be sent to the south
!
  CALL INTERSECTION(TZNORTHEXTHALO2, NPROC, TZPZS(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDNORTHSOUTH
      CALL ADD_ZONE( TPCOMDATA%TSEND_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!*       4.4   Find out the zones to be sent to the east
!
  CALL INTERSECTION(TZWESTEXTHALO2, NPROC, TZPZS(IP), TZINTER)
  DO JI = 1, NPROC
    IF((TZINTER(JI)%NUMBER.NE.0).AND.(TZINTER(JI)%NUMBER.NE.IP)) THEN
      TZINTER(JI)%MSSGTAG = ISENDWESTEAST
      CALL ADD_ZONE( TPCOMDATA%TSEND_HALO2, TZINTER(JI) )
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
!*       5.    MODIFICATIONS IN CASE OF CYCLIC CONDITIONS :
!              ------------------------------------------
!
  CALL INI_CYCLIC2( TPCOMDATA%TSEND_HALO2, &
                    TPCOMDATA%TRECV_HALO2, &
                    TZPZS )
!
!
!-------------------------------------------------------------------------------
!
!*       6.    SWITCH FROM GLOBAL COORDINATES TO LOCAL COORDINATES :
!              ---------------------------------------------------
!
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TSEND_HALO2)
  CALL GLOBAL2LOCAL(TPPROCONF, TPCOMDATA%TRECV_HALO2)
!
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE CONSTRUCT_HALO2
!
!     #############################################################
      SUBROUTINE SET_EXT_HALO2( TPPZS, TPWESTHALO2, TPSOUTHHALO2, &
                                TPNORTHHALO2, TPEASTHALO2 )
!     #############################################################
!
!!****  *SET_EXT_HALO2* - initialization of location variable
!                         for second layer halo
!
!!    Purpose
!!    -------
!     The purpose of this routine is to set the values of the
!     TPWESTEXTHALO2, TPSOUTHEXTHALO2, TPNORTHEXTHALO2 and TPEASTEXTHALO2
!     variables. This variables contain the location of the second layer
!     of the halo, for each subdomain
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        type ZONE_ll
!
!!    Reference
!!    ---------
!
!!    Author
!!    ------
!
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original 25/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(ZONE_ll), DIMENSION(:) :: TPPZS            ! Physical Zone Splitting
  TYPE(ZONE_ll), DIMENSION(:) :: TPWESTHALO2   ! zones modeling
  TYPE(ZONE_ll), DIMENSION(:) :: TPSOUTHHALO2  ! each side of
  TYPE(ZONE_ll), DIMENSION(:) :: TPNORTHHALO2  ! the second layer of the
  TYPE(ZONE_ll), DIMENSION(:) :: TPEASTHALO2   ! halo
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    SET TPWESTHALO2
!
  TPWESTHALO2(:)%MSSGTAG= 4
  TPWESTHALO2(:)%NUMBER = TPPZS(:)%NUMBER
  TPWESTHALO2(:)%NXOR   = TPPZS(:)%NXOR - 2
  TPWESTHALO2(:)%NXEND  = TPPZS(:)%NXOR - 2
  TPWESTHALO2(:)%NYOR   = TPPZS(:)%NYOR
  TPWESTHALO2(:)%NYEND  = TPPZS(:)%NYEND
!
!-------------------------------------------------------------------------------
!
!*       2.    SET TPEASTHALO2
!
  TPEASTHALO2(:)%NUMBER = TPPZS(:)%NUMBER
  TPEASTHALO2(:)%MSSGTAG= 8
  TPEASTHALO2(:)%NXOR   = TPPZS(:)%NXEND+2
  TPEASTHALO2(:)%NXEND  = TPPZS(:)%NXEND+2
  TPEASTHALO2(:)%NYOR   = TPPZS(:)%NYOR
  TPEASTHALO2(:)%NYEND  = TPPZS(:)%NYEND
!
!-------------------------------------------------------------------------------
!
!*       3.    SET TPNORTHHALO2
!
  TPNORTHHALO2(:)%NUMBER = TPPZS(:)%NUMBER
  TPNORTHHALO2(:)%MSSGTAG= 2
  TPNORTHHALO2(:)%NXOR   = TPPZS(:)%NXOR
  TPNORTHHALO2(:)%NXEND  = TPPZS(:)%NXEND
  TPNORTHHALO2(:)%NYOR   = TPPZS(:)%NYEND + 2
  TPNORTHHALO2(:)%NYEND  = TPPZS(:)%NYEND + 2
!
!-------------------------------------------------------------------------------
!
!*       4.    SET TPSOUTHHALO2
!
  TPSOUTHHALO2(:)%NUMBER = TPPZS(:)%NUMBER
  TPSOUTHHALO2(:)%MSSGTAG= 6
  TPSOUTHHALO2(:)%NXOR   = TPPZS(:)%NXOR
  TPSOUTHHALO2(:)%NXEND  = TPPZS(:)%NXEND
  TPSOUTHHALO2(:)%NYOR   = TPPZS(:)%NYOR - 2
  TPSOUTHHALO2(:)%NYEND  = TPPZS(:)%NYOR - 2
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE SET_EXT_HALO2
!
!     ###############################################
      SUBROUTINE INI_CYCLIC2( TPSEND, TPRECV, TPPZS )
!     ###############################################
!
!!****  *INI_CYCLIC2* - routine to updates the correspondants'lists
!                       in case of cyclic conditions
!!
!!    Purpose
!!    -------
!     The purpose of this routine is to complete the correspondants'lists
!     for the processors on the boundaries of the domain with the informations
!     on cyclic conditions
!
!!**  Method
!!    ------
!     For the different cases of cyclic conditions and the different boundaries,
!     we list the processors and construct the intersections to determine
!     which zones to exchange
!
!!    External
!!    --------
!
!     Module MODE_TOOLS_ll
!        LWEST_ll, LSOUTH_ll, LEAST_ll, LNORTH_ll
!        ADD_ZONE
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        types CRSPD_ll, ZONE_ll, PROC_COM_DATA_ll
!
!     Module MODD_VAR_ll
!        IP - Number of local processor=subdomain
!        NPROC - Number of processors
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!
!     Module MODD_DIM_ll
!        CLBCX - X-direction LBC type at left(1) and right(2) boundaries
!        CLBCY - Y-direction LBC type at left(1) and right(2) boundaries
!
!!    Reference
!!    ---------
!
!!    Authors
!!    -------
!
!     R. Guivarch               * CERFACS - ENSEEIHT *
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original 28/05/98
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!  USE MODD_STRUCTURE_ll, ONLY : CRSPD_ll, ZONE_ll, PROC_COM_DATA_ll
!  USE MODD_VAR_ll, ONLY       : TCRRT_COMDATA, DIMZ, NPROC, IP
!  USE MODD_DIM_ll, ONLY       : CLBCX, CLBCY
!
  USE MODE_TOOLS_ll, ONLY     : LWEST_ll, LSOUTH_ll, LEAST_ll, LNORTH_ll, &
                                ADD_ZONE
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(CRSPD_ll), POINTER                 :: TPSEND ! Correspondants
                                                    ! for send commmunications
  TYPE(CRSPD_ll), POINTER                 :: TPRECV ! Correspondants
                                                    ! for revc commmunications
!
  TYPE(ZONE_ll), DIMENSION(:), INTENT(IN) :: TPPZS  ! Physical Zone splitting
                                                    ! of the global domain
!
!*       0.2   declarations of local variables
!
  INTEGER       :: J ! loop control variable
  INTEGER       :: IXMIN, IXMAX, IYMIN, IYMAX ! intermediate variables
                                              ! for intersections
  INTEGER       :: ISENDWESTEAST, ISENDEASTWEST, &  ! MPI message tags
                   ISENDNORTHSOUTH, ISENDSOUTHNORTH
!
  TYPE(ZONE_ll) :: TZZONE_INTER !Intermediate zone variable
!
  INTEGER       :: ICURMODEL
!
!-------------------------------------------------------------------------------
!
  ICURMODEL = TCRRT_COMDATA%NUMBER
!
!*       1.    SET MPI MESSAGE TAGS :
!              --------------------
!
  ISENDWESTEAST   = 14
  ISENDEASTWEST   = 18
  ISENDNORTHSOUTH = 12
  ISENDSOUTHNORTH = 16
!
!-------------------------------------------------------------------------------
!
!*       2.    WEST BOUNDARY :
!              -------------
!
  IF(CLBCX(ICURMODEL, 1) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the WEST boundary
!
    IF(LWEST_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the EAST ones
        IF(LEAST_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the physical subdomain of processor j
!         to determine the send and receive zones (for processor IP)
          IYMIN = MAX( TPPZS(IP)%NYOR, TPPZS(J)%NYOR)
          IYMAX = MIN( TPPZS(IP)%NYEND, TPPZS(J)%NYEND)
!
          IF(IYMIN <= IYMAX) THEN        ! the intersection is not void
!
!           add the adequate zone to the send correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDWESTEAST, &
                                 TPPZS(IP)%NXOR+1,&
                                 TPPZS(IP)%NXOR+1,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ )

!!$            TZZONE_INTER = ZONE_ll( J, ISENDWESTEAST, &
!!$                                 TPPZS(IP)%NXOR+1, IYMIN, 1, &
!!$                                 TPPZS(IP)%NXOR+1, IYMAX, DIMZ )

            CALL ADD_ZONE( TPSEND, TZZONE_INTER )
!
!           add the adequate zone to the receive correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDEASTWEST, &
                                 TPPZS(IP)%NXOR-2,&
                                 TPPZS(IP)%NXOR-2,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, ISENDEASTWEST, &
!!$                                 TPPZS(IP)%NXOR-2, IYMIN, 1, &
!!$                                 TPPZS(IP)%NXOR-2, IYMAX, DIMZ)

            CALL ADD_ZONE( TPRECV, TZZONE_INTER )
          ENDIF
!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.    EAST BOUNDARY :
!              -------------
!
  IF(CLBCX(ICURMODEL, 2) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the EAST boundary
!
    IF(LEAST_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the WEST ones
        IF(LWEST_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the physical subdomain of processor j
!         to determine the send and receive zones (for processor IP)
!
          IYMIN = MAX( TPPZS(IP)%NYOR, TPPZS(J)%NYOR)
          IYMAX = MIN( TPPZS(IP)%NYEND, TPPZS(J)%NYEND)
!
          IF(IYMIN <= IYMAX) THEN     ! the intersection is not void
!
!           add the adequate zone to the send correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDEASTWEST, &
                                 TPPZS(IP)%NXEND - 1,&
                                 TPPZS(IP)%NXEND - 1,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ)
!!$
!!$            TZZONE_INTER = ZONE_ll( J, ISENDEASTWEST, &
!!$                                 TPPZS(IP)%NXEND - 1, IYMIN, 1, &
!!$                                 TPPZS(IP)%NXEND - 1, IYMAX, DIMZ)

            CALL ADD_ZONE( TPSEND, TZZONE_INTER )
!
!           add the adequate zone to the receive correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDWESTEAST, &
                                 TPPZS(IP)%NXEND + 2,&
                                 TPPZS(IP)%NXEND + 2,&
                                 IYMIN,&
                                 IYMAX,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, ISENDWESTEAST, &
!!$                                 TPPZS(IP)%NXEND + 2, IYMIN, 1, &
!!$                                 TPPZS(IP)%NXEND + 2, IYMAX, DIMZ)

            CALL ADD_ZONE( TPRECV, TZZONE_INTER )
          ENDIF
!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.    SOUTH BOUNDARY :
!              --------------
!
  IF(CLBCY(ICURMODEL, 1) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the SOUTH boundary
!
    IF(LSOUTH_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the NORTH ones
!
        IF(LNORTH_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the extended subdomain of processor j
!         to determine the send zone (for processor IP)
!
          IXMIN = MAX( TPPZS(IP)%NXOR, TPPZS(J)%NXOR)
          IXMAX = MIN( TPPZS(IP)%NXEND, TPPZS(J)%NXEND)
!
          IF(IXMIN <= IXMAX) THEN ! the intersection is not void
!
!           add the adequate zone to the send correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDSOUTHNORTH, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPPZS(IP)%NYOR + 1,&
                                 TPPZS(IP)%NYOR + 1,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, ISENDSOUTHNORTH, &
!!$                                 IXMIN, TPPZS(IP)%NYOR + 1, 1, &
!!$                                 IXMAX, TPPZS(IP)%NYOR + 1, DIMZ)

            CALL ADD_ZONE( TPSEND, TZZONE_INTER )
!
!           add the adequate zone to the receive correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDNORTHSOUTH, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPPZS(IP)%NYOR - 2,&
                                 TPPZS(IP)%NYOR - 2,&
                                 1, &
                                 DIMZ)

!!$            TZZONE_INTER = ZONE_ll( J, ISENDNORTHSOUTH, &
!!$                                 IXMIN, TPPZS(IP)%NYOR - 2, 1, &
!!$                                 IXMAX, TPPZS(IP)%NYOR - 2, DIMZ)

            CALL ADD_ZONE( TPRECV, TZZONE_INTER )
          ENDIF
!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.    NORTH BOUNDARY :
!              --------------
!
  IF(CLBCY(ICURMODEL, 2) .EQ. 'CYCL') THEN
!
!    the local subdomain (processor) is on the NORTH boundary
!
    IF(LNORTH_ll()) THEN
!
      DO J = 1, NPROC
!
!       list all the subdomains and select the SOUTH ones
        IF(LSOUTH_ll(J)) THEN
!
!         compute the intersection between the physical subdomain
!         of processor IP with the physical subdomain of processor j
!         to determine the send and receive zones (for processor IP)
!
          IXMIN = MAX( TPPZS(IP)%NXOR, TPPZS(J)%NXOR)
          IXMAX = MIN( TPPZS(IP)%NXEND, TPPZS(J)%NXEND)
!
          IF(IXMIN <= IXMAX) THEN ! the intersection is not void
!
!           add the adequate zone to the send correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDNORTHSOUTH, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPPZS(IP)%NYEND - 1,&
                                 TPPZS(IP)%NYEND - 1,&
                                 1, &
                                 DIMZ )

!!$            TZZONE_INTER = ZONE_ll( J, ISENDNORTHSOUTH, &
!!$                                 IXMIN, TPPZS(IP)%NYEND - 1, 1, &
!!$                                 IXMAX, TPPZS(IP)%NYEND - 1, DIMZ )

            CALL ADD_ZONE( TPSEND, TZZONE_INTER )
!
!           add the adequate zone to the receive correspondant
!
            TZZONE_INTER = ZONE_ll( J, ISENDSOUTHNORTH, &
                                 IXMIN,&
                                 IXMAX,&
                                 TPPZS(IP)%NYEND + 2,&
                                 TPPZS(IP)%NYEND + 2,&
                                 1, &
                                 DIMZ )

!!$            TZZONE_INTER = ZONE_ll( J, ISENDSOUTHNORTH, &
!!$                                 IXMIN, TPPZS(IP)%NYEND + 2, 1, &
!!$                                 IXMAX, TPPZS(IP)%NYEND + 2, DIMZ )

            CALL ADD_ZONE( TPRECV, TZZONE_INTER )
          ENDIF
!
        ENDIF
!
      ENDDO
!
    ENDIF
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE INI_CYCLIC2
!
!!     #######################################
       SUBROUTINE CLEANLIST_LCRSPD( TPLCRSPD )
!!     #######################################
!
!!****  *CLEANLIST_LCRSPD* - routine to destroy a list of CRSPD
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
!     Module MODE_ARGSLIST_ll
!        CLEANLIST_ll
!
!!    Implicit Arguments
!!    ------------------
!
!     Module MODD_STRUCTURE_ll
!        type CRSPD_ll
!
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!!
!!    Modifications
!!    -------------
!     Original 01/05/98
!
!------------------------------------------------------------------------------
!
!  USE MODD_STRUCTURE_ll, ONLY : LCRSPD_ll
!
  USE MODE_ARGSLIST_ll, ONLY  : CLEANLIST_ll
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(LCRSPD_ll), POINTER :: TPLCRSPD
!
!*       0.2   declarations of local variables
!
  TYPE(LCRSPD_ll), POINTER :: TZLCRSPD
!
!------------------------------------------------------------------------------
!
!*       1.    Dealloacte one by one the elements of TPLIST
!              --------------------------------------------
!
  TZLCRSPD => TPLCRSPD
  DO WHILE (ASSOCIATED(TZLCRSPD))
    CALL CLEANLIST_ll(TZLCRSPD%TLIST)
    TZLCRSPD => TZLCRSPD%TNEXT
  ENDDO
!
!------------------------------------------------------------------------------
!
       END SUBROUTINE CLEANLIST_LCRSPD
!
END MODULE MODE_CONSTRUCT_ll
