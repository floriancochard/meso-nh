!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 les 2006/05/18 13:07:25
!-----------------------------------------------------------------
!      ###################
MODULE MODI_LES_MEAN_1PROC
!      ###################
!
INTERFACE LES_MEAN_1PROC
!
      SUBROUTINE LES_MEAN_1PROC_2D(PA, OMASK,                   &
                                   PA_MEAN, KAVG_PTS, KUND_PTS  )
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PA_MEAN
INTEGER,                 INTENT(OUT) :: KAVG_PTS
INTEGER,                 INTENT(OUT) :: KUND_PTS
!
END SUBROUTINE LES_MEAN_1PROC_2D
!
      SUBROUTINE LES_MEAN_1PROC_3D(PA, OMASK,                   &
                                   PA_MEAN, KAVG_PTS, KUND_PTS  )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),     INTENT(OUT) :: PA_MEAN
INTEGER, DIMENSION(:),     INTENT(OUT) :: KAVG_PTS
INTEGER, DIMENSION(:),     INTENT(OUT) :: KUND_PTS
!
END SUBROUTINE LES_MEAN_1PROC_3D
!
      SUBROUTINE LES_MEAN_1PROC_3DM(PA, OMASK,                   &
                                    PA_MEAN, KAVG_PTS, KUND_PTS  )

REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),     INTENT(OUT) :: PA_MEAN
INTEGER, DIMENSION(:),     INTENT(OUT) :: KAVG_PTS
INTEGER, DIMENSION(:),     INTENT(OUT) :: KUND_PTS
!
END SUBROUTINE LES_MEAN_1PROC_3DM
!
END INTERFACE
!
END MODULE MODI_LES_MEAN_1PROC
!
!     ###########################################################
      SUBROUTINE LES_MEAN_1PROC_2D(PA, OMASK,                   &
                                   PA_MEAN, KAVG_PTS, KUND_PTS  )
!     ###########################################################
!
!
!!****  *LES_MEAN_1PROC* computes the average of one field on one processor
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
!!      V. Masson        06/11/02 optimization
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:), INTENT(IN)  :: OMASK
!
REAL,                    INTENT(OUT) :: PA_MEAN
INTEGER,                 INTENT(OUT) :: KAVG_PTS
INTEGER,                 INTENT(OUT) :: KUND_PTS
!
!       0.2  declaration of local variables
!
INTEGER :: IIB, IIE, IJB, IJE ! physical domain boundary
INTEGER :: JI, JJ             ! loop counters
!
!-------------------------------------------------------------------------------
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
!
PA_MEAN  = 0.
KAVG_PTS = 0
KUND_PTS = 0
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IF ( OMASK(JI,JJ) ) THEN
      IF ( PA(JI,JJ)/=XUNDEF ) THEN
        PA_MEAN  = PA_MEAN  + PA(JI,JJ)
        KAVG_PTS = KAVG_PTS + 1
      ELSE
        KUND_PTS = KUND_PTS + 1
      END IF
    END IF
  END DO
END DO
!
IF ( KAVG_PTS > 0 ) THEN
  PA_MEAN = PA_MEAN / KAVG_PTS
ELSE
  PA_MEAN = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_1PROC_2D
!
!     ###########################################################
      SUBROUTINE LES_MEAN_1PROC_3D(PA, OMASK,                   &
                                   PA_MEAN, KAVG_PTS, KUND_PTS  )
!     ###########################################################
!
!
!!****  *LES_MEAN_1PROC* computes the average of one field on one processor
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
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),     INTENT(OUT) :: PA_MEAN
INTEGER, DIMENSION(:),     INTENT(OUT) :: KAVG_PTS
INTEGER, DIMENSION(:),     INTENT(OUT) :: KUND_PTS
!
!
!       0.2  declaration of local variables
!
INTEGER :: IIB, IIE, IJB, IJE ! physical domain boundary
INTEGER :: JI, JJ             ! loop counters
!
INTEGER :: IK                 ! number of vertical levels
INTEGER :: JK                 ! loop counter
!-------------------------------------------------------------------------------
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
IK=SIZE(PA,3)
!
PA_MEAN  = 0.
KAVG_PTS = 0
KUND_PTS = 0
!
DO JK=1,IK
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IF ( OMASK(JI,JJ,JK) ) THEN
        IF ( PA(JI,JJ,JK)/=XUNDEF ) THEN
          PA_MEAN (JK) = PA_MEAN (JK) + PA(JI,JJ,JK)
          KAVG_PTS(JK) = KAVG_PTS(JK) + 1
        ELSE
          KUND_PTS(JK) = KUND_PTS(JK) + 1
        END IF
      END IF
    END DO
  END DO
END DO
!
WHERE ( KAVG_PTS(:) > 0 )
  PA_MEAN(:) = PA_MEAN(:) / KAVG_PTS(:)
ELSEWHERE
  PA_MEAN(:) = XUNDEF
END WHERE
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_1PROC_3D
!
!     ###########################################################
      SUBROUTINE LES_MEAN_1PROC_3DM(PA, OMASK,                   &
                                    PA_MEAN, KAVG_PTS, KUND_PTS  )
!     ###########################################################
!
!
!!****  *LES_MEAN_1PROC* computes the average of one field on one processor
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
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
!
REAL,    DIMENSION(:),     INTENT(OUT) :: PA_MEAN
INTEGER, DIMENSION(:),     INTENT(OUT) :: KAVG_PTS
INTEGER, DIMENSION(:),     INTENT(OUT) :: KUND_PTS
!
!
!       0.2  declaration of local variables
!
INTEGER :: IIB, IIE, IJB, IJE ! physical domain boundary
INTEGER :: JI, JJ             ! loop counters
!
INTEGER :: IK                 ! number of vertical levels
INTEGER :: JK                 ! loop counter
!-------------------------------------------------------------------------------
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
IK=SIZE(PA,3)
!
PA_MEAN  = 0.
KAVG_PTS = 0
KUND_PTS = 0
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IF ( OMASK(JI,JJ) ) THEN
      DO JK=1,IK
        IF ( PA(JI,JJ,JK)/=XUNDEF ) THEN
          PA_MEAN (JK) = PA_MEAN (JK) + PA(JI,JJ,JK)
          KAVG_PTS(JK) = KAVG_PTS(JK) + 1
        ELSE
          KUND_PTS(JK) = KUND_PTS(JK) + 1
        END IF
      END DO
    END IF
  END DO
END DO
!
WHERE ( KAVG_PTS(:) > 0 )
  PA_MEAN(:) = PA_MEAN(:) / KAVG_PTS(:)
ELSEWHERE
  PA_MEAN(:) = XUNDEF
END WHERE
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_MEAN_1PROC_3DM

