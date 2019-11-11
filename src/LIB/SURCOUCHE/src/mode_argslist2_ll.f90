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

!!     ########################
       MODULE MODE_ARGSLIST2_ll
!!     ########################
!
!!****  *MODD_ARGSLIST_ll* - declaration of lists type 
!
!!     Routines Of The User Interface
!!     ------------------------------
!
!      SUBROUTINES : ADD_FIELD2_ll, DEL_FIELD2_ll
! 
!!     Purpose
!!     -------
!      This module manages a list of "halos2". A "halo2" is a variable of 
!      type HALO2_ll. The type HALO2_ll contains four 2D array, 
!      one for each boundary of the subdomain (see module MODD_STRUCTURE2).
!      Halos may be added (routine ADD_FIELD2_ll) or deleted 
!      (routine DEL_FIELD2_ll) from the list. The list can then be used 
!      to update the second layer of the halo.
! 
!!     Reference
!!     ---------
! 
!      User interface for the Meso-NH parallel package
! 
!!     Authors
!!     -------
! 
!      Ph. Kloos                 * CERFACS - CNRM *
!
!!     Implicit Arguments
!!     ------------------
!
!      Module MODD_STRUCTURE2 : contains type HALO2_ll
!
!      Module MODD_ARGSLIST : contains type HALO2LIST_ll, lists of halos2
!
!
!!     Modifications
!!     -------------
!      Original    May 19, 1998
!
!-------------------------------------------------------------------------------
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2_ll, HALO2LIST_ll
!
  CONTAINS
!
!
!!     ###############################################
       SUBROUTINE ADD_FIELD2_ll(TPLIST_ll, TPHALO2_ll)
!!     ###############################################
!
!!****  *ADD_FIELD2_ll* - 
!
!!     Purpose
!!     -------
!      This routine is used to add a halo (TPHALO2_ll) to a list 
!      of halos (HALO2LIST_ll)
!
!!     Reference
!!     ---------
!
!      User interface for the Meso-NH parallel package
!
!!     Implicit Arguments
!!     ------------------
!
!      Module MODD_ARGSLIST
!         HALO2_ll, HALO2LIST_ll - structure and list of halo2
!
!!     Author
!!     ------
!
!      Ph. Kloos                 * CERFACS - CNRM *
!
!!     Modifications
!!     -------------
!      Original    May 19, 1998
!
!-------------------------------------------------------------------------------
!
  IMPLICIT NONE
!
!*       0.    DECLARATIONS
!
!
!*       0.1   declarations of arguments
!
  TYPE(HALO2LIST_ll), POINTER :: TPLIST_ll  ! list of HALO2
  TYPE(HALO2_ll), TARGET      :: TPHALO2_ll ! HALO2 to be added
!
!*       0.2   declarations of local variables
!
  TYPE(HALO2LIST_ll), POINTER :: TZLIST_ll
!
!-------------------------------------------------------------------------------
!
!*       1.    In case TPLIST_ll has not been already used
!              -------------------------------------------
!
  IF (.NOT.ASSOCIATED(TPLIST_ll)) THEN ! TPLIST_ll is empty
!
    ALLOCATE(TPLIST_ll)
    TPLIST_ll%HALO2 => TPHALO2_ll
    NULLIFY(TPLIST_ll%NEXT)
    TPLIST_ll%NCARD = 1
!
  ELSE 
!
!
!*       2.    TPLIST already contains fields;
!*             add the PFIELD field at the end of the list
!              -------------------------------------------
!
    TZLIST_ll => TPLIST_ll
    DO WHILE(ASSOCIATED(TZLIST_ll%NEXT))
      TZLIST_ll => TZLIST_ll%NEXT
    ENDDO
!
    ALLOCATE(TZLIST_ll%NEXT)
    TZLIST_ll => TZLIST_ll%NEXT
    TZLIST_ll%HALO2 => TPHALO2_ll
    NULLIFY(TZLIST_ll%NEXT)
!
    TPLIST_ll%NCARD = TPLIST_ll%NCARD + 1
!
  ENDIF
!
!-------------------------------------------------------------------------------
!
       END SUBROUTINE ADD_FIELD2_ll
!
!!     ######################################################
       SUBROUTINE DEL_FIELD2_ll(TPLIST_ll, TPHALO2_ll, KINFO)
!!     ######################################################
!
!!****  *DEL_FIELD2_ll* -
!
!!     Purpose
!!     -------
!      This routine deletes the TPHALO2_ll halo from the list
!      of halos (HALO2LIST_ll)
!
!!     Reference
!!     ---------
!
!      User interface for the Meso-NH parallel package
!
!!     Implicit Arguments
!!     ------------------
!
!      Module MODD_ARGSLIST
!         HALO2_ll, HALO2LIST_ll - structure and list of halo2
!
!!     Author
!!     ------
!
!      Ph. Kloos                 * CERFACS - CNRM *
!
!!     Modifications
!!     -------------
!      Original    May 19, 1998
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  TYPE(HALO2LIST_ll), POINTER :: TPLIST_ll
  TYPE(HALO2_ll), TARGET :: TPHALO2_ll
  INTEGER :: KINFO ! return status : 0 if PFIELD has been found in TPLIST
                   !                 1 otherwise
!
!
!*       0.2   declarations of local variables
!
  TYPE(HALO2LIST_ll), POINTER :: TZLIST, TZTMP
!
!-------------------------------------------------------------------------------
!
!*       1.    Check presence of TPHALO2_ll in the first element of the list
!              -------------------------------------------------------------
!
  KINFO=1
  DO WHILE (ASSOCIATED(TPLIST_ll%HALO2, TPHALO2_ll)) 
    IF (ASSOCIATED(TPLIST_ll%NEXT)) THEN
      TZTMP => TPLIST_ll%NEXT
      TZTMP%NCARD = TPLIST_ll%NCARD - 1
    ENDIF
    TPLIST_ll => TPLIST_ll%NEXT
    KINFO = 0
  END DO
!
!-------------------------------------------------------------------------------
!
!*       2.    Check presence of TPHALO2_ll in other elements of TPLIST_ll
!*             There may be several occurences of TPHALO2_ll
!              --------------------------------------------------------------
!
  TZLIST => TPLIST_ll
  DO WHILE(ASSOCIATED(TZLIST))
    TZTMP => TZLIST
    TZLIST => TZLIST%NEXT
    IF (ASSOCIATED(TZLIST)) THEN
      IF (ASSOCIATED(TZLIST%HALO2, TPHALO2_ll)) THEN
        TZTMP%NEXT => TZLIST%NEXT
        TZLIST => TZTMP
        KINFO = 0
      ENDIF
    ENDIF
  ENDDO
!
!-------------------------------------------------------------------------------
!
       END SUBROUTINE DEL_FIELD2_ll
!
END MODULE MODE_ARGSLIST2_ll
