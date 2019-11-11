!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/goto_model_surfex_mnh.f90
!-----------------------------------------------------------------
!#######################
MODULE MODI_GOTO_MODEL_SURFEX_MNH
  !#######################
  !
  INTERFACE
    !     ###############################
          SUBROUTINE GOTO_MODEL_SURFEX_MNH(U,KMI, KINFO_ll)
    !     ###############################
    !!
    !!    PURPOSE
    !!    -------
    !!
    !!    Initializes local sizes in SURFEX module MODD_SURF_ATM_n for model KMI
    !!    and calls GOTO_MODEL(KMI)
    !!             GOTO_MODEL_ll(KMI, KINFO_ll)
    !!
    !!    METHOD
    !!    ------
    !!
    !!    EXTERNAL
    !!    --------
    !!
    !!
    !!    IMPLICIT ARGUMENTS
    !!    ------------------
    !!
    !!
    !!    REFERENCE
    !!    ---------
    !!
    !!    AUTHOR
    !!    ------
    !!
    !!    M. Moge                   LA - CNRS
    !!
    !!    MODIFICATION
    !!    ------------
    !!
    !!    Original      08/2015
    !----------------------------------------------------------------------------
    !
    !*    0.     DECLARATION
    !            -----------
    !
    USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
    !
    IMPLICIT NONE
    !
    !*    0.1    Declaration of dummy arguments
    !            ------------------------------
    !
    TYPE(SURF_ATM_t), INTENT(INOUT) :: U
    INTEGER,                         INTENT(IN)    :: KMI    !model id
    INTEGER,                         INTENT(OUT)    :: KINFO_ll
          END SUBROUTINE GOTO_MODEL_SURFEX_MNH
  !
  END INTERFACE
  !
END MODULE MODI_GOTO_MODEL_SURFEX_MNH
!     ###############################
      SUBROUTINE GOTO_MODEL_SURFEX_MNH(U,KMI, KINFO_ll)
!     ###############################
!!
!!    PURPOSE
!!    -------
!!
!!    Initializes local sizes in SURFEX module MODD_SURF_ATM_n for model KMI
!!
!!    METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    M. Moge                   LA - CNRS
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original      08/2015
!!  06/2016     (G.Delautier) phasage surfex 8
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
!USE MODE_ll
USE MODE_TOOLS_ll, ONLY : GET_GLOBALDIMS_ll, GET_DIM_PHYS_ll
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT
USE MODD_VAR_ll, ONLY : YSPLITTING
USE MODE_MODELN_HANDLER
!
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
INTEGER,                         INTENT(IN)    :: KMI    !model id
INTEGER,                         INTENT(OUT)    :: KINFO_ll
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!INTEGER :: IINFO_ll ! return code of // routines
INTEGER :: IMI ! return code of // routines
CHARACTER*1 :: HSPLIT
!
!------------------------------------------------------------------------------
!
!IMI = GET_CURRENT_MODEL_INDEX()
!IF ( KMI /= IMI ) THEN
!  WRITE(*,*) "ERROR IN PGD_GOTO_MODEL_SURFEX_MNH : current MNH model ", IMI, " is different from target SURFEX model ", KMI
!  CALL ABORT
!ENDIF
!! il faudrait faire ce test mais cette GOTO_SURFEX ne marche pas, donc on laisse tomber et on bricole
!
!------------------------------------------------------------------------------
!
!CALL GOTO_SURFEX(KMI, LKFROM)  ! cette routine plante, donc on ne touche pas aux modeles
!CALL GO_TOMODEL_ll(KMI,KINFO_ll)
!CALL GOTO_MODEL(KMI)
!
CALL GET_GLOBALDIMS_ll(U%NIMAX_SURF_ll,U%NJMAX_SURF_ll)
U%NDIM_FULL = U%NIMAX_SURF_ll*U%NJMAX_SURF_ll
!
IF ( YSPLITTING == "BSPLITTING" ) THEN
  HSPLIT = 'B'
ELSE IF ( YSPLITTING == "XSPLITTING" ) THEN
  HSPLIT = 'X'
ELSE IF ( YSPLITTING == "YSPLITTING" ) THEN
  HSPLIT = 'Y'
ELSE
  HSPLIT = ''
ENDIF
CALL GET_DIM_PHYS_ll(HSPLIT,U%NIMAX_SURF_LOC,U%NJMAX_SURF_LOC)
U%NSIZE_FULL = U%NIMAX_SURF_LOC*U%NJMAX_SURF_LOC
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GOTO_MODEL_SURFEX_MNH
