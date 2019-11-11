!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_NEST_FIELD_n
!     ######################
!
INTERFACE
!
      SUBROUTINE NEST_FIELD_n(HFIELD)
!
CHARACTER(LEN=6),   INTENT(IN)  :: HFIELD ! name of the field to nest
!
END SUBROUTINE NEST_FIELD_n
!
END INTERFACE
!
END MODULE MODI_NEST_FIELD_n
!
!
!
!     #############################
      SUBROUTINE NEST_FIELD_n(HFIELD)
!     #############################
!
!!****  *NEST_FIELD$n* - make a field coherent between model $n and its sons
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        27/09/96
!!                      30/07/97 (Masson) group MODI_FILL_SONFIELDn
!!                      04/08/97 (Masson) correction of land fraction
!!                      15/03/99 (Masson) new cover types
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CONF
USE MODD_NESTING
USE MODD_PARAMETERS
USE MODD_NEST_PGD_n
USE MODD_GRID_n
USE MODD_DIM_n
!
USE MODI_FILL_SONFIELD_n
!
USE MODE_MODELN_HANDLER
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: HFIELD ! name of the field to nest
!
!
!*       0.2   declarations of local variables
!
INTEGER :: JMI   ! model index
INTEGER :: JLSON ! index in XNESTFIELD array
INTEGER :: IMI   ! current model index
REAL, DIMENSION(:,:,:,:), POINTER :: DPTR_XNESTFIELD ! dummy pointer to correct ifort bug
!
!-------------------------------------------------------------------------------
!
!*       1.    Allocations
!              -----------
!
SELECT CASE(HFIELD)
!
!* modification of orography
!
  CASE ('ZS    ')
    IMI = GET_CURRENT_MODEL_INDEX()
    ALLOCATE ( XNESTFIELD(NIMAX+2*JPHEXT,NJMAX+2*JPHEXT,1+COUNT(NDAD(:)==IMI),1))
!
!
  CASE DEFAULT
    RETURN
END SELECT
!
XNESTFIELD(:,:,:,:) = 0.
!
!-------------------------------------------------------------------------------
!
!*       2.    Retrieve information from son(s)
!              --------------------------------
!
DO JLSON=1,SIZE(NSON)
  DPTR_XNESTFIELD=>XNESTFIELD ! correct an ifort bug with following call
  CALL FILL_SONFIELD_n(NSON(JLSON),HFIELD,DPTR_XNESTFIELD,JLSON)
END DO
!
!-------------------------------------------------------------------------------
!
!*       3.    Update father
!              -------------
!
SELECT CASE(HFIELD)
!
!* modification of orography
!
  CASE ('ZS    ')
!
    XZS(:,:)=SUM(XNESTFIELD(:,:,:,1)*NNESTMASK(:,:,:),DIM=3)
!
  CASE DEFAULT
    RETURN
END SELECT
!
!-------------------------------------------------------------------------------
!
!*       4.    Deallocations
!              -------------
!
DEALLOCATE ( XNESTFIELD )
!-------------------------------------------------------------------------------
!
END SUBROUTINE NEST_FIELD_n
