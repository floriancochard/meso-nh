!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##################
      MODULE MODI_FLUX 
!     ##################
!
INTERFACE
!
FUNCTION FXM(PA,PRUCT)  RESULT(PFXM)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRUCT  ! u-component 
                                                            ! contravariant velocity     
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PFXM   ! result at flux
                                                            ! side
END FUNCTION FXM
!
FUNCTION FYM(PA,PRVCT)  RESULT(PFYM)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRVCT  ! v-component 
                                                            ! contravariant velocity     
                                                            ! antidiffusive vel.     
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PFYM   ! result at flux
                                                            ! side
END FUNCTION FYM
!
FUNCTION FZM(PA,PRWCT)  RESULT(PFZM)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRWCT  ! w-component 
                                                            ! contravariant velocity     
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PFZM   ! result at flux
                                                            ! side
END FUNCTION FZM
!
END INTERFACE
!
END MODULE MODI_FLUX
!
!     #####################################
      FUNCTION FXM(PA,PRUCT)  RESULT(PFXM)
!     #####################################
!
!!****  *FXM* -  MPDATA flux operator : flux operator in x direction
!!                                      for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a flux operator 
!     along the x direction (I index) for a field PA localized at a mass
!     point. This flux operator is used to calculate the advection of a quantity 
!     according to the upstream scheme. The result is localized at a x-flux point (u point).
!
!!**  METHOD
!!    ------ 
!!      The result PFXM(i,:,:) is defined by 
!!    (PA(i-1,:,:)*MAX(0.,PRUCT(i,:,:)) + PA(i,:,:)*MIN(0.,PRUCT(i,:,:))
!!    At i=1, PFXM(1,:,:) are replaced by the values of PFXM,
!!    which are the right values in the x-cyclic case. 
!!    
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
!!      MPDATA advection scheme (book1)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Vila-Guerau       * Meteo France *
!!	J.-P. Lafore         * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/09/95 
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
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRUCT  ! u-component  
                                                            ! contravariant velocity      
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PFXM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! Size of the array in the x direction
!
!-------------------------------------------------------------------------------
!
!*       0.3    DEFINITION OF FXM
!              ------------------
!
IIU = SIZE(PA,1)
!
DO JI=2,IIU
  PFXM(JI,:,:) =  PA(JI-1,:,:)*AMAX1(0.,PRUCT(JI,:,:))    & 
                 +PA(JI,:,:)  *AMIN1(0.,PRUCT(JI,:,:)) 
END DO
!
PFXM(1,:,:)    = PFXM(IIU-2*JPHEXT+1,:,:) 
!
!-------------------------------------------------------------------------------
!
END FUNCTION FXM
!
!
!     #####################################
      FUNCTION FYM(PA,PRVCT)  RESULT(PFYM)
!     #####################################
!
!!****  *FYM* -  MPDATA flux operator : flux operator in y direction
!!                                      for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a flux operator 
!     along the y direction (J index) for a field PA localized at a mass
!     point. This flux operator is used to calculate the advection of a quantity 
!     according to the upwind scheme. The result is localized at a y-flux point (v point).
!
!!**  METHOD
!!    ------ 
!!      The result PFYM(:,j,:) is defined by 
!!    (PA(:,j-1,:)*MAX(0.,PRVCT(:,j,:) + PA(:,j,:)*MIN(0.,PRVCT(:,j,:))
!!    At j=1, PFYM(:,j,:) are replaced by the values of PFYM,
!!    which are the right values in the y-cyclic case. 
!!    
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
!!      MPDATA advection scheme (book1)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Vila-Guerau       * Meteo France *
!!	J.-P. Lafore         * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/09/95 
!-------------------------------------------------------------------------------
!
!*       1.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       1.1   Declarations of argument and result
!              ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRVCT  ! v-component 
                                                            ! contravariant velocity     
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PFYM   ! result at flux
                                                            ! side
!
!*       1.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!-------------------------------------------------------------------------------
!
!*       1.3    DEFINITION OF FYM
!              ------------------
!
IJU = SIZE(PA,2)
!
DO JJ=2,IJU
  PFYM(:,JJ,:) =  PA(:,JJ-1,:)*AMAX1(0.,PRVCT(:,JJ,:))     &
                 +PA(:,JJ,:)  *AMIN1(0.,PRVCT(:,JJ,:)) 
END DO
!
PFYM(:,1,:)    = PFYM(:,IJU-2*JPHEXT+1,:)  
!
!-------------------------------------------------------------------------------
!
END FUNCTION FYM
!
!
!     #####################################
      FUNCTION FZM(PA,PRWCT)  RESULT(PFZM)
!     #####################################
!
!!****  *FZM* -  MPDATA flux operator : flux operator in z direction
!!                                      for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a flux operator 
!     along the z direction (K index) for a field PA localized at a mass
!     point. This flux operator is used to calculate the advection of a quantity 
!     according to the upwind scheme. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------ 
!!      The result PFZM(:,:,k) is defined by 
!!    (PA(:,:,k-1)*MAX(0.,PRWCT(:,:,k) + PA(:,:,k)*MIN(0.,PRWCT(:,;,k))
!!    At k=1, PFZM(:,:,k) is defined by -999. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      MPDATA advection scheme (book1)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Vila-Guerau       * Meteo France *
!!	J.-P. Lafore         * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/09/95 
!-------------------------------------------------------------------------------
!
!*       2.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       2.1   Declarations of argument and result
!              ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PRWCT  ! w-component 
                                                            ! contravariant velocity     
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PFZM   ! result at flux
                                                            ! side
!
!*       2.2    Declarations of local variables
!               -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
INTEGER :: IKU            ! Size of the array in the z direction
!
!-------------------------------------------------------------------------------
!
!*       2.3    DEFINITION OF FZM
!                ------------------
!
IKU = SIZE(PA,3)
!
DO JK=2,IKU
  PFZM(:,:,JK) =  PA(:,:,JK-1)*AMAX1(0.,PRWCT(:,:,JK))      &
                 +PA(:,:,JK)  *AMIN1(0.,PRWCT(:,:,JK)) 
END DO
!
PFZM(:,:,1)    =  -999. 
!
!-------------------------------------------------------------------------------
!
END FUNCTION FZM
