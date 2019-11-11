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

!     ######spl
      MODULE MODI_SCATTER
!     #####################
!
INTERFACE
      SUBROUTINE SCATTER(P1,P2)
!
REAL, DIMENSION(:,:),  INTENT(IN)    :: P1
REAL, DIMENSION(:,:),  INTENT(OUT)   :: P2
!
END SUBROUTINE SCATTER
!
END INTERFACE
!
END MODULE MODI_SCATTER
