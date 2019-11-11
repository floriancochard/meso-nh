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

!     #######################
      MODULE MODI_LOCATION_ll
!     #######################
!
INTERFACE
!
!!     ###########################################
       LOGICAL FUNCTION LNORTH_ll( K, HSPLITTING )
!!     ###########################################
!
  INTEGER, INTENT(IN), OPTIONAL     :: K ! number of the subdomain
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
!
       END FUNCTION LNORTH_ll
!
!!     ##########################################
       LOGICAL FUNCTION LWEST_ll( K, HSPLITTING )
!!     ##########################################
!
  INTEGER, INTENT(IN), OPTIONAL     :: K ! number of the subdomain
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
!
       END FUNCTION LWEST_ll
!
!!     ###########################################
       LOGICAL FUNCTION LSOUTH_ll( K, HSPLITTING )
!!     ###########################################
!
  INTEGER, INTENT(IN), OPTIONAL     :: K ! number of the subdomain
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
!
       END FUNCTION LSOUTH_ll
!
!!     ##########################################
       LOGICAL FUNCTION LEAST_ll( K, HSPLITTING )
!!     ##########################################
!
  INTEGER, INTENT(IN), OPTIONAL     :: K ! number of the subdomain
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
!
       END FUNCTION LEAST_ll
!
END INTERFACE
!
END MODULE MODI_LOCATION_ll
