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

!     ###################
      MODULE MODI_NEST_ll
!     ###################
!
INTERFACE
!
!     ###############################################
      SUBROUTINE GET_MODEL_NUMBER_ll( KMODEL_NUMBER )
!     ###############################################
!
  INTEGER :: KMODEL_NUMBER
!
      END SUBROUTINE GET_MODEL_NUMBER_ll
!
!     ####################################################
      SUBROUTINE GET_CHILD_DIM_ll( KCHILD, KX, KY, KINFO )
!     ####################################################
!
  INTEGER, INTENT(IN) :: KCHILD
!
  INTEGER, INTENT(OUT) :: KX, KY
!
  INTEGER, INTENT(OUT) :: KINFO
!
       END SUBROUTINE GET_CHILD_DIM_ll
!
!     ###################################################################
      SUBROUTINE GET_FEEDBACK_COORD_ll( KXOR, KYOR, KXEND, KYEND, KINFO )
!     ###################################################################
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
!
  INTEGER, INTENT(OUT) :: KINFO
!
       END SUBROUTINE GET_FEEDBACK_COORD_ll
!
!     ##################################
      SUBROUTINE UNSET_LSFIELD_1WAY_ll()
!     ##################################
!
      END SUBROUTINE UNSET_LSFIELD_1WAY_ll
!
!     ##########################################
      SUBROUTINE UNSET_LSFIELD_2WAY_ll( KMODEL )
!     ##########################################
!
  INTEGER, INTENT(IN) :: KMODEL
!
      END SUBROUTINE UNSET_LSFIELD_2WAY_ll
!
!     #########################################
      SUBROUTINE LS_FORCING_ll( KCHILD, KINFO, OEXTRAPOL, OCYCLIC_EXTRAPOL )
!     #########################################
!
  INTEGER, INTENT(IN) :: KCHILD 
  INTEGER, INTENT(OUT) :: KINFO
  LOGICAL, OPTIONAL, INTENT(IN) :: OEXTRAPOL
  LOGICAL, OPTIONAL, INTENT(IN) :: OCYCLIC_EXTRAPOL
!
      END SUBROUTINE LS_FORCING_ll
!
!     ###################################
      SUBROUTINE LS_FEEDBACK_ll( KINFO  )
!     ###################################
!
  INTEGER, INTENT(OUT) :: KINFO
!
      END SUBROUTINE LS_FEEDBACK_ll
!

!     #############################
      SUBROUTINE UNSET_LBFIELD_ll()
!     #############################
!
      END SUBROUTINE UNSET_LBFIELD_ll
!
!     #########################################
      SUBROUTINE LB_FORCING_ll( KCHILD, KINFO )
!     #########################################
!
  INTEGER, INTENT(IN) :: KCHILD
!
  INTEGER, INTENT(OUT) :: KINFO
!
       END SUBROUTINE LB_FORCING_ll
!
!!     ###########################################################
       FUNCTION LBFINE2COARSE( KRATIO, KLBSIZE ) RESULT( KCOARSE )
!!     ###########################################################
!
  IMPLICIT NONE
!
  INTEGER :: KCOARSE
!
  INTEGER :: KRATIO, KLBSIZE
!
      END FUNCTION LBFINE2COARSE
!
!     #########################################
      SUBROUTINE GO_TOMODEL_ll( KMODEL, KINFO )
!     #########################################
!
INTEGER :: KMODEL, KINFO
!
      END SUBROUTINE GO_TOMODEL_ll
!
END INTERFACE
!
INTERFACE SET_LSFIELD_1WAY_ll
!
!     #############################################################
      SUBROUTINE SET_LS2DFIELD_1WAY_ll( P2DFIELD, PTFIELD, KMODEL )
!     #############################################################
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
  INTEGER, INTENT(IN) :: KMODEL
!
      END SUBROUTINE SET_LS2DFIELD_1WAY_ll
!
!     #############################################################
      SUBROUTINE SET_LS3DFIELD_1WAY_ll( P3DFIELD, PTFIELD, KMODEL )
!     #############################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
  INTEGER, INTENT(IN) :: KMODEL
!
      END SUBROUTINE SET_LS3DFIELD_1WAY_ll
!
END INTERFACE
!
INTERFACE SET_LSFIELD_2WAY_ll
!
!     #####################################################
      SUBROUTINE SET_LS2DFIELD_2WAY_ll( P2DFIELD, PTFIELD )
!     #####################################################
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
!
      END SUBROUTINE SET_LS2DFIELD_2WAY_ll
!
!     #####################################################
      SUBROUTINE SET_LS3DFIELD_2WAY_ll( P3DFIELD, PTFIELD )
!     #####################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
!
      END SUBROUTINE SET_LS3DFIELD_2WAY_ll
!
END INTERFACE
!
INTERFACE SET_LBFIELD_ll
!
!     ##############################################################
      SUBROUTINE SET_LB2DFIELD_ll( P2DFIELD, PTFIELD, KFINELBSIZE, &
                                   HSIDE, KMODEL )
!     ##############################################################
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
!
  INTEGER, INTENT(IN) :: KFINELBSIZE, KMODEL
!
  CHARACTER(LEN=*), INTENT(IN) :: HSIDE
!
      END SUBROUTINE SET_LB2DFIELD_ll
!
!     ##############################################################
      SUBROUTINE SET_LB3DFIELD_ll( P3DFIELD, PTFIELD, KFINELBSIZE, &
                                   HSIDE, KMODEL )
!     ##############################################################
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
!
  INTEGER, INTENT(IN) :: KFINELBSIZE, KMODEL
!
  CHARACTER(LEN=*), INTENT(IN) :: HSIDE
!
      END SUBROUTINE SET_LB3DFIELD_ll
!
END INTERFACE
!
END MODULE MODI_NEST_ll
