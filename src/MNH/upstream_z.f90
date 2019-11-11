!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 forcing 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_UPSTREAM_Z 
!     ######################
!
INTERFACE
!
FUNCTION UPSTREAM_Z(PA,PRWCT)  RESULT(PUPSTRZ)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! advected tendency
                                                            ! at flux height
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRWCT  !  w-component 
                                                            ! contravariant 
                                                            !   velocity     
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PUPSTRZ! vertical upstream
                                                            ! advection 
END FUNCTION UPSTREAM_Z
!
END INTERFACE
!
END MODULE MODI_UPSTREAM_Z
!
!     ##############################################
      FUNCTION UPSTREAM_Z(PA,PRWCT)  RESULT(PUPSTRZ)
!     ##############################################
!
!!****  *UPSTREAM_Z* -  vertical upstream operator
!!
!!    PURPOSE
!!    -------
!!      The purpose of this function  is to compute a vertically advection
!!    (ascending or subsiding motion) with an upstream differencing scheme. 
!!
!!**  METHOD
!!    ------ 
!!      The upstream advection PUPSTRZ(i,:,:) is selected automatically
!!    according to the sign of PRWCT:
!!        - if PRWCT>0  the advective tendency PA is actually valid at the 
!!    advected variable location
!!        - if PRWCT<0  the advective tendency PA is shifted downward to be
!!    valid at the advected variable location
!!    For the upper boundary condition, in case of subsidence, the gradient
!!    from level IKE-1 is copied onto level IKE, which means that the upper
!!    boundary satisfies a "constant gradient" condition.
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      NONE
!!
!!
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty  * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/03/96 
!!      P.Jabouille 15/10/98 set default values at the first and last vertical levels
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! advected tendency
                                                            ! at flux height
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRWCT  !  w-component 
                                                            ! contravariant 
                                                            !   velocity     
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PUPSTRZ! vertical upstream
                                                            ! advection
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: ZA     ! work array
!
INTEGER :: IKB, IKE
!
!-------------------------------------------------------------------------------
!
IKB=1+JPVEXT
IKE=SIZE(PRWCT,3) - JPVEXT
!
!*       1.0    DEFINITION OF UPSTREAM_Z
!               ------------------------
!
ZA(:,:,:) = EOSHIFT(PA(:,:,:),SHIFT=1,DIM=3)
ZA(:,:,IKE) = ZA(:,:,IKE-1)
!
PUPSTRZ(:,:,IKB:IKE) = (0.5+SIGN(0.5,PRWCT(:,:,IKB:IKE)))  &
                                         *PA(:,:,IKB:IKE)  &
                     + (0.5-SIGN(0.5,PRWCT(:,:,IKB:IKE)))  &
                                         *ZA(:,:,IKB:IKE) 
PUPSTRZ(:,:,IKB-1)=999.
PUPSTRZ(:,:,IKE+1)=999.
!
!-------------------------------------------------------------------------------
!
END FUNCTION UPSTREAM_Z
