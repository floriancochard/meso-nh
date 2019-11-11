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

!     #########################
      MODULE MODE_BOUNDARIES_ll
!     #########################
!
!!****  *MODE_BOUNDARIES_ll* -
!
!!    Purpose
!!    -------
!
!!    Routines Of The User Interface
!!    ------------------------------
!
!     SUBROUTINE : UPDATE_BOUNDARIES_ll
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
!     N. Gicquel                * CNRM *
!
!!    Implicit Arguments
!!    ------------------
!       Module MODD_VAR_ll :
!         TCRRT_COMDATA : information about the communications on the current
!                         model
!         TCRRT_COMDATA%TSEND_BOUNDX  : list of regions to be sent
!         TCRRT_COMDATA%TSEND_BOUNDY  : when updating the boundaries
!         TCRRT_COMDATA%TSEND_BOUNDXY : (type CRSPD_ll)
!
!       Module MODD_ARGSLIST_ll :
!         LIST_ll : type for a list of fields
!
!-------------------------------------------------------------------------------
!
  CONTAINS
!
!     ###########################################################
      SUBROUTINE UPDATE_BOUNDARIES_ll( HDIRECTION, TPLIST, KINFO)
!     ###########################################################
!
!!****  *UPDATE_BOUNDARIES_ll* -
!
!!    Purpose
!!    -------
!     This routine updates the boundaries with the values computed by the
!     neighbor subdomains only for processors on the physical global
!     domain's boundaries.
!     The fields to be updated are in the
!     TPLIST list of fields. Before UPDATE_BOUNDARIES_ll is called, TPLIST
!     has been filled with the fields to be communicated
!
!!**  Method
!!    ------
!
!     According to the selected direction, each processor sends and
!     receives the fields on their boundaries
! 
!!    Externals
!!    ---------
!
!     Module MODE_EXCHANGE_ll :
!        SEND_RECV_CRSPD, COPY_CRSPD
!
!!    Implicit Arguments
!!    ------------------
!       Module MODD_VAR_ll :
!         TCRRT_COMDATA : information about the communications on the current
!                         model
!         TCRRT_COMDATA%TSEND_BOUNDX  : list of regions to be sent
!         TCRRT_COMDATA%TSEND_BOUNDY  : when updating the boundaries
!         TCRRT_COMDATA%TSEND_BOUNDXY : (type CRSPD_ll)
!
!       Module MODD_ARGSLIST_ll :
!         LIST_ll : type for a list of fields
!
!!    Authors
!!    -------
!
!     R. Guivarch               * CERFACS - ENSEEIHT *
!     Ph. Kloos                 * CERFACS - CNRM *
!     N. Gicquel                * CNRM *
!
!!    Modifications
!!    -------------
!       Original    June 15, 1998
!       N. Gicquel  October 30, 1998 COPY_CRSPD
!
!-------------------------------------------------------------------------------
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA, NHALO_COM
!
  USE MODE_EXCHANGE_ll, ONlY : SEND_RECV_CRSPD, COPY_CRSPD
!
  IMPLICIT NONE
!
!*       0.    DECLARATIONS
!
!
!*       0.1   declarations of arguments
!
  CHARACTER*2, INTENT(IN) :: HDIRECTION
  TYPE(LIST_ll), POINTER  :: TPLIST ! pointer to the list of fields
                                    ! to be updated
  INTEGER                 :: KINFO  ! return status
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    IF TPLIST IS NOT EMPTY
!              ----------------------
  IF(ASSOCIATED(TPLIST)) THEN
!
!*       1.1   XX Boundaries
!
    IF (HDIRECTION.EQ.'XX') THEN
      CALL SEND_RECV_CRSPD(TCRRT_COMDATA%TSEND_BOUNDX, &
                           TCRRT_COMDATA%TRECV_BOUNDX, &
                           TPLIST, TPLIST, NHALO_COM, KINFO)

      CALL COPY_CRSPD(TCRRT_COMDATA%TSEND_BOUNDX, &
                      TCRRT_COMDATA%TRECV_BOUNDX, &
                      TPLIST,TPLIST,KINFO)
    ENDIF
!
!*       1.2   YY Boundaries
!
    IF (HDIRECTION.EQ.'YY') THEN
      CALL SEND_RECV_CRSPD(TCRRT_COMDATA%TSEND_BOUNDY, &
                           TCRRT_COMDATA%TRECV_BOUNDY, &
                           TPLIST, TPLIST, NHALO_COM, KINFO)

      CALL COPY_CRSPD(TCRRT_COMDATA%TSEND_BOUNDY, &
                      TCRRT_COMDATA%TRECV_BOUNDY, TPLIST, &
                      TPLIST,KINFO)

    ENDIF
!
!*       1.3   XY Boundaries
!
   IF (HDIRECTION.EQ.'XY') THEN
     CALL SEND_RECV_CRSPD(TCRRT_COMDATA%TSEND_BOUNDXY, &
                          TCRRT_COMDATA%TRECV_BOUNDXY, &
                          TPLIST, TPLIST, NHALO_COM, KINFO)

     CALL COPY_CRSPD(TCRRT_COMDATA%TSEND_BOUNDXY, &
                     TCRRT_COMDATA%TRECV_BOUNDXY, TPLIST, &
                     TPLIST,KINFO)
!
   ENDIF
 ELSE
   KINFO = -1
 ENDIF
!
!-------------------------------------------------------------------------------
!
      END SUBROUTINE UPDATE_BOUNDARIES_ll
!
END MODULE MODE_BOUNDARIES_ll
