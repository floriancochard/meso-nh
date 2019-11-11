!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_TOTAL_DMASS
!     #######################
!
INTERFACE
!
      SUBROUTINE TOTAL_DMASS(PJ,PRHOD,PDRYMASS)
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! Jacobian             
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOD     ! dry density
!
REAL,                   INTENT(OUT) :: PDRYMASS ! Mass of dry air Md
                                          !  contained in the simulation domai
!
END SUBROUTINE TOTAL_DMASS
!
END INTERFACE
!
END MODULE MODI_TOTAL_DMASS
!
!
!
!
!     #########################################################################
      SUBROUTINE TOTAL_DMASS(PJ,PRHOD,PDRYMASS)
!     #########################################################################
!
!!****  *TOTAL_DMASS* - routine to set reference state for anelastic approximation
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set the total mass of dry air Md
!
!!**  METHOD
!!    ------
!!      
!!      The total mass of dry air Md is taken as the one corresponding
!!    to the reference atmosphere.
!!      This assumption is made for a simulation preparation of an ideal case.
!!    For a real case, the total mass of dry air Md, can be supplied by the
!!    LS coupling fields.
!!      
!!    EXTERNAL
!!    --------   
!!      none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS : contains declaration of parameter variables
!!       
!!         JPHEXT   : Horizontal external points number
!!         JPVEXT   : Vertical external points number
!!
!!      Module MODD_CONF   : contains configuration variables
!!      
!!       NVERB  : verbosity level
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine TOTAL_DMASS)
!!      
!!
!!    AUTHOR
!!    ------
!!	J.P. Lafore     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        06/02/95  
!!      J.-P. Pinty     16/12/99  Parallelization
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_CONF,    ONLY: NVERB
USE MODD_LUNIT_n, ONLY: TLUOUT
!
USE MODE_ll
!  
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! Jacobian             
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHOD     ! dry density
!
REAL,                   INTENT(OUT) :: PDRYMASS ! Mass of dry air Md
                                          !  contained in the simulation domain                        
!
!*       0.2   declarations of local variables
!
INTEGER             :: IINFO_ll   ! Return code of parallel routine
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE TOTAL MASS OF DRY AIR   
!	        ----------------------------------
!
PDRYMASS=SUM3D_ll(PJ(:,:,:)*PRHOD(:,:,:),IINFO_ll)
!
!-------------------------------------------------------------------------------
!
!*       3.    PRINT ON OUTPUT-LISTING
!              -----------------------
!
IF(NVERB >= 5) THEN  
  WRITE(TLUOUT%NLU,*) 'TOTAL_DMASS:    M= ',PDRYMASS
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TOTAL_DMASS
