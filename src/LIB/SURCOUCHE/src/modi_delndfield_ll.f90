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

!     ########################
      MODULE MODI_DELnDFIELD_ll
!     #########################
!
INTERFACE
!
!!     #################################################
       SUBROUTINE DEL1DFIELD_ll( TPLIST, PFIELD, KINFO )
!!     #################################################
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER     :: TPLIST ! list of fields
  REAL, DIMENSION(:), TARGET :: PFIELD ! field to be deleted
                                       ! from the list of fields
  INTEGER, INTENT(OUT)       :: KINFO  ! return status :
                                       !   0 if PFIELD has been found
                                       !   1 otherwise.
!
       END SUBROUTINE DEL1DFIELD_ll
!
!!     #################################################
       SUBROUTINE DEL2DFIELD_ll( TPLIST, PFIELD, KINFO )
!!     #################################################
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER     :: TPLIST ! list of fields
  REAL, DIMENSION(:,:), TARGET :: PFIELD ! field to be deleted
                                       ! from the list of fields
  INTEGER, INTENT(OUT)       :: KINFO  ! return status :
                                       !   0 if PFIELD has been found
                                       !   1 otherwise.
!
       END SUBROUTINE DEL2DFIELD_ll
!
!!     #################################################
       SUBROUTINE DEL3DFIELD_ll( TPLIST, PFIELD, KINFO )
!!     #################################################
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER     :: TPLIST ! list of fields
  REAL, DIMENSION(:,:,:), TARGET :: PFIELD ! field to be deleted
                                       ! from the list of fields
  INTEGER, INTENT(OUT)       :: KINFO  ! return status :
                                       !   0 if PFIELD has been found
                                       !   1 otherwise.
!
       END SUBROUTINE DEL3DFIELD_ll
!
!!     #################################
       SUBROUTINE CLEANLIST_ll( TPLIST )
!!     #################################
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER     :: TPLIST ! list of fields
!
       END SUBROUTINE CLEANLIST_ll
!
END INTERFACE
!
END MODULE MODI_DELnDFIELD_ll
