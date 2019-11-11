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

!     ###########################
      MODULE MODI_ADDDELFIELD2_ll
!     ###########################
!
INTERFACE
!
!!     #################################################
       SUBROUTINE ADD_FIELD2_ll( TPLIST_ll, TPHALO2_ll )
!!     #################################################
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2_ll, HALO2LIST_ll
!
  TYPE(HALO2LIST_ll), POINTER :: TPLIST_ll  ! list of HALO2
  TYPE(HALO2_ll), TARGET      :: TPHALO2_ll ! HALO2 to be added
!
       END SUBROUTINE ADD_FIELD2_ll
!
!!     ########################################################
       SUBROUTINE DEL_FIELD2_ll( TPLIST_ll, TPHALO2_ll, KINFO )
!!     ########################################################
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2_ll, HALO2LIST_ll
!
  TYPE(HALO2LIST_ll), POINTER :: TPLIST_ll ! list of fields
  TYPE(HALO2_ll), TARGET      :: TPHALO2_ll! field to be deleted
                                           ! from the list of fields
  INTEGER                     :: KINFO     ! return status :
                                           !   0 if PFIELD has been found
                                           !   1 otherwise
!
       END SUBROUTINE DEL_FIELD2_ll
!
END INTERFACE
!
END MODULE MODI_ADDDELFIELD2_ll

