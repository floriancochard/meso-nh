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
      MODULE MODI_ADDnDFIELD_ll
!     #########################
!
INTERFACE
!
!!     ################################################
       SUBROUTINE ADD1DFIELD_ll( HDIR, TPLIST, PFIELD )
!!     ################################################
!
  USE MODD_ARGSLIST_ll, ONLY : LIST1D_ll
!
  CHARACTER(LEN=1), INTENT(IN) :: HDIR   ! direction of the field
                                         ! ("X" or "Y")
  TYPE(LIST1D_ll), POINTER     :: TPLIST ! list of fields
  REAL, DIMENSION(:), TARGET   :: PFIELD ! field to be added 
                                         ! to the list of fields
!
       END SUBROUTINE ADD1DFIELD_ll
!
!!     ##########################################
       SUBROUTINE ADD2DFIELD_ll( TPLIST, PFIELD )
!!     ##########################################
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER       :: TPLIST ! list of fields
  REAL, DIMENSION(:,:), TARGET :: PFIELD ! field to be added
                                         ! to the list of fields
!
       END SUBROUTINE ADD2DFIELD_ll
!
!!     ##########################################
       SUBROUTINE ADD3DFIELD_ll( TPLIST, PFIELD )
!!     ##########################################
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER         :: TPLIST   ! list of fields
  REAL, DIMENSION(:,:,:), TARGET :: PFIELD   ! field to be added
                                             ! to the list of fields
!
       END SUBROUTINE ADD3DFIELD_ll
!
END INTERFACE
!
END MODULE MODI_ADDnDFIELD_ll
