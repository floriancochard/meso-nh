!MNH_LIC Copyright 2000-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #################
MODULE MODI_SPEC_VER_INT
!      #################
!
INTERFACE
!
      SUBROUTINE SPEC_VER_INT(KMI,PA_MNH, PA_SPEC)

INTEGER,                   INTENT(IN)  :: KMI 
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA_MNH
!
REAL,    DIMENSION(:,:,:), INTENT(OUT) :: PA_SPEC
!
END SUBROUTINE SPEC_VER_INT
!
END INTERFACE
!
END MODULE MODI_SPEC_VER_INT
!
!     ######################################
      SUBROUTINE SPEC_VER_INT(KMI, PA_MNH, PA_SPEC)
!     ######################################
!
!
!!****  *SPEC_VER_INT* interpolates a MESONH field
!!                    on the LES output levels
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_LES
USE MODD_PARAMETERS
!
USE MODE_ll
USE MODE_GATHER_ll
USE MODI_VER_INTERP_LIN
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
INTEGER,                   INTENT(IN)  :: KMI 
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA_MNH
!
REAL,    DIMENSION(:,:,:), INTENT(OUT) :: PA_SPEC
!
!
!       0.2  declaration of local variables
!
INTEGER :: IRESP
INTEGER :: JK     ! vertical loop counter
!
REAL, DIMENSION(SIZE(PA_MNH,1),SIZE(PA_MNH,2),SIZE(PA_SPEC,3)) :: ZA
REAL, DIMENSION(:,:,:), ALLOCATABLE                            :: ZA_ll
!
INTEGER :: IIMAX_ll ! total physical domain size (x)
INTEGER :: IJMAX_ll ! total physical domain size (y)
!
REAL, DIMENSION(SIZE(PA_SPEC,3)) :: ZA_MEAN
INTEGER :: JI, JJ ! loop counter
!-------------------------------------------------------------------------------
!
IF (CSPECTRA_LEVEL_TYPE=='N') THEN
  PA_SPEC(:,:,:) = XUNDEF
  RETURN
ELSE IF (CSPECTRA_LEVEL_TYPE=='K') THEN
  DO JK = 1, NSPECTRA_K
    ZA(:,:,JK) = PA_MNH(:,:,NSPECTRA_LEVELS(JK))
  END DO
ELSE IF (CSPECTRA_LEVEL_TYPE=='Z') THEN
  PA_SPEC = VER_INTERP_LIN(PA_MNH,NKLIN_CURRENT_SPEC,XCOEFLIN_CURRENT_SPEC)
  !
  WHERE(NKLIN_CURRENT_SPEC<2)
    PA_SPEC = XUNDEF
  END WHERE
ELSE
  PRINT*, '-------> STOP in SPEC_VER_INT <----------'
 !callabortstop
CALL ABORT
  STOP
END IF
!
!-------------------------------------------------------------------------------
!
! ONE PROCESSOR ONLY
! ------------------
!
CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
ALLOCATE(ZA_ll(IIMAX_ll+2*JPHEXT,IJMAX_ll+2*JPHEXT,NSPECTRA_K))
!
!ZA_ll(:,:,:) = ZA(:,:,:)
CALL GATHERALL_FIELD_ll('XY',ZA,ZA_ll,IRESP)
!
!-------------------------------------------------------------------------------
!
PA_SPEC(:,:,:) = ZA_ll(NLESn_IINF(KMI):NLESn_ISUP(KMI),NLESn_JINF(KMI):NLESn_JSUP(KMI),:)
!
DO JK=1,SIZE(PA_SPEC,3)
  ZA_MEAN(JK) = SUM(PA_SPEC(:,:,JK)) / SIZE(PA_SPEC(:,:,JK))
END DO
DO JK=1,SIZE(PA_SPEC,3)
  DO JJ=1,SIZE(PA_SPEC,2)
    DO JI=1,SIZE(PA_SPEC,1)
       PA_SPEC(JI,JJ,JK) = PA_SPEC(JI,JJ,JK) - ZA_MEAN(JK)
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SPEC_VER_INT
