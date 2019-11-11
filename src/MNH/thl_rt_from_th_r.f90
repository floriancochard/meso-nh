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
!    ############################
     MODULE MODI_THL_RT_FROM_TH_R  
!    ############################ 
!
INTERFACE
!
      SUBROUTINE THL_RT_FROM_TH_R(OUSERV, OUSERC, OUSERR,             &
                                  OUSERI, OUSERS, OUSERG, OUSERH,     &
                                  PL_O_EXN_CP,                        &
                                  PTH, PR, PTHL, PRT                  )
!
LOGICAL, INTENT(IN) :: OUSERV  ! flag to use water vapor
LOGICAL, INTENT(IN) :: OUSERC  ! flag to use water cloud
LOGICAL, INTENT(IN) :: OUSERR  ! flag to use rain
LOGICAL, INTENT(IN) :: OUSERI  ! flag to use ice cloud
LOGICAL, INTENT(IN) :: OUSERS  ! flag to use snow
LOGICAL, INTENT(IN) :: OUSERG  ! flag to use graupel
LOGICAL, INTENT(IN) :: OUSERH  ! flag to use hail
REAL, DIMENSION(:,:,:), INTENT(IN) :: PL_O_EXN_CP ! L/(Exn*Cp)
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTH    ! theta (or theta tendency) to
                                             ! transform into th_l (or tendency)
REAL, DIMENSION(:,:,:,:),INTENT(IN):: PR     ! mixing ratios (or tendency) to
                                             ! transform into total water
                                             ! (or its tendency)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PTHL   ! th_l (or th_l tendency)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PRT    ! total water (or its tendency)
!-------------------------------------------------------------------------------
!
END SUBROUTINE THL_RT_FROM_TH_R
!
END INTERFACE
!
END MODULE MODI_THL_RT_FROM_TH_R
!
!     #################################################################
      SUBROUTINE THL_RT_FROM_TH_R(OUSERV, OUSERC, OUSERR,             &
                                  OUSERI, OUSERS, OUSERG, OUSERH,     &
                                  PL_O_EXN_CP,                        &
                                  PTH, PR, PTHL, PRT                  )
!     #################################################################
!
!
!!****  *THL_RT_FROM_TH_R* - computes the conservative variables
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
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
!!      V. Masson               * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         20/09/02
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_PARAMETERS
USE MODD_CST
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
LOGICAL, INTENT(IN) :: OUSERV  ! flag to use water vapor
LOGICAL, INTENT(IN) :: OUSERC  ! flag to use water cloud
LOGICAL, INTENT(IN) :: OUSERR  ! flag to use rain
LOGICAL, INTENT(IN) :: OUSERI  ! flag to use ice cloud
LOGICAL, INTENT(IN) :: OUSERS  ! flag to use snow
LOGICAL, INTENT(IN) :: OUSERG  ! flag to use graupel
LOGICAL, INTENT(IN) :: OUSERH  ! flag to use hail
REAL, DIMENSION(:,:,:), INTENT(IN) :: PL_O_EXN_CP ! L/(Exn*Cp)
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTH    ! theta (or theta tendency) to
                                             ! transform into th_l (or tendency)
REAL, DIMENSION(:,:,:,:),INTENT(IN):: PR     ! mixing ratios (or tendency) to
                                             ! transform into total water
                                             ! (or its tendency)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PTHL   ! th_l (or th_l tendency)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PRT    ! total water (or its tendency)
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER                                                    :: IRR
! water variable counter
!----------------------------------------------------------------------------
!
!*      1.2  preliminary calculations
!            ------------------------
!
!
!* latent heat of vaporization at air temperature
! 
IRR=2
!
IF (OUSERC) THEN
  PRT(:,:,:) = PR(:,:,:,1) + PR(:,:,:,2)
  IF (OUSERR) THEN
    IRR=IRR+1
    PRT(:,:,:) = PRT(:,:,:) +  PR(:,:,:,IRR)
  END IF
  IF (OUSERI) THEN
    IRR=IRR+1
    PRT(:,:,:) = PRT(:,:,:) +  PR(:,:,:,IRR)
  END IF
  IF (OUSERS) THEN
    IRR=IRR+1
    PRT(:,:,:) = PRT(:,:,:) +  PR(:,:,:,IRR)
  END IF
  IF (OUSERG) THEN
    IRR=IRR+1
    PRT(:,:,:) = PRT(:,:,:) +  PR(:,:,:,IRR)
  END IF
  IF (OUSERH) THEN
    IRR=IRR+1
    PRT(:,:,:) = PRT(:,:,:) +  PR(:,:,:,IRR)
  END IF
  !
  !* liquid potential temperature
  !
  PTHL(:,:,:) = PTH(:,:,:) - PR(:,:,:,2)*PL_O_EXN_CP
  !
ELSE
  PRT(:,:,:) = 0.
  PTHL(:,:,:) = PTH(:,:,:)
  IF (OUSERV) PRT(:,:,:) = PR(:,:,:,1)
END IF

END SUBROUTINE THL_RT_FROM_TH_R
