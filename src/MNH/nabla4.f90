!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 operators 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_NABLA4
!     ##################
!
INTERFACE
!
FUNCTION DX4(PA)  RESULT(PDX4)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA   ! field variable
REAL, DIMENSION(SIZE(PA,1)-4,SIZE(PA,2),SIZE(PA,3)) :: PDX4 ! result
END FUNCTION DX4
!
FUNCTION DY4(PA)  RESULT(PDY4)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA   ! field variable
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)-4,SIZE(PA,3)) :: PDY4 ! result
END FUNCTION DY4
!
FUNCTION DX4_2(PA)  RESULT(PDX4_2)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! field variable
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDX4_2 ! result
END FUNCTION DX4_2
!
FUNCTION DY4_2(PA)  RESULT(PDY4_2)
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! field variable
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDY4_2 ! result
END FUNCTION DY4_2
!
END INTERFACE
!
END MODULE MODI_NABLA4
!     ######spl
      FUNCTION DX4(PA)  RESULT(PDX4)
!     ###############################
!
!!****  *DX4* -  Nabla_4 operator in X direction applied to any variable
!!               with circular shift on the boundaries
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute a fourth-order finite
!     difference along the X direction (I index) for a field PA. The
!     result PDX4 is localized on the same grid of the field PA.
!
!!**  METHOD
!!    ------ 
!!        The result PDX4(i,:,:) is defined by PA(i-2,:,:)-4*PA(i-1,:,:)
!!        +6*PA(i,:,:)-4*PA(i+1,:,:)-PA(i+2,:,:).
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
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty and J.-P. Lafore      * LA & Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/10/94 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA   ! field variable
REAL, DIMENSION(SIZE(PA,1)-4,SIZE(PA,2),SIZE(PA,3)) :: PDX4 ! result
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in X direction
INTEGER :: IIU            ! Upper bound in X direction of PA 
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DX4
!              ------------------
!
IIU = SIZE(PA,1)
!
PDX4(1:IIU-4,:,:) = PA(1:IIU-4,:,:) + PA(5:IIU,:,:)            &
          - 4.0*( PA(2:IIU-3,:,:) + PA(4:IIU-1,:,:) )          &
          + 6.0*PA(3:IIU-2,:,:)
!
!-------------------------------------------------------------------------------
!
END FUNCTION DX4
!     ######spl
      FUNCTION DX4_2(PA)  RESULT(PDX4_2)
!     ##################################
!
!!****  *DX4_2* -  Nabla_4 operator in X direction applied to any variable
!!                 in case of open conditions which degenerates to a nabla_2
!!                 operator on the boundaries
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a forth-order finite
!     difference along the X direction (I index) for a field PA and a
!     second order finite difference on the boundaries. The result PDX4_2
!     is localized on the same grid of the field PA.
!
!!**  METHOD
!!    ------ 
!!        The result PDX4_2(i,:,:) is defined by PA(i-2,:,:)-4*PA(i-1,:,:)
!!        +6*PA(i,:,:)-4*PA(i+1,:,:)-PA(i+2,:,:) everywhere except on the
!!        boundaries where is reduced to PA(i-1,:,:)-2*PA(i,:,:)+PA(i+1,:,:).
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
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/10/94 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDX4_2   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in X direction
INTEGER :: IIU            ! Upper bound in X direction of PA 
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DX4_2
!              ------------------
!
IIU = SIZE(PA,1)
!
PDX4_2(3:IIU-2,:,:) = PA(1:IIU-4,:,:) + PA(5:IIU,:,:)  &
            - 4.0*( PA(2:IIU-3,:,:) + PA(4:IIU-1,:,:) ) + 6.0*PA(3:IIU-2,:,:)
!
PDX4_2(2,:,:)     = PA(1,:,:)     - 2.*PA(2,:,:)     + PA(3,:,:)
PDX4_2(IIU-1,:,:) = PA(IIU-2,:,:) - 2.*PA(IIU-1,:,:) + PA(IIU,:,:)
!
PDX4_2(1,:,:)   = 0.0
PDX4_2(IIU,:,:) = 0.0
!
!-------------------------------------------------------------------------------
!
END FUNCTION DX4_2
!     ######spl
      FUNCTION DY4(PA)  RESULT(PDY4)
!     ###############################
!
!!****  *DY4* -  Nabla_4 operator in Y direction applied to any variable
!!               with circular shift on the boundaries
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute a fourth-order finite
!     difference along the Y direction (J index) for a field PA. The
!     result PDY4 is localized on the same grid of the field PA.
!
!!**  METHOD
!!    ------ 
!!        The result PDY4(:,j,:) is defined by PA(:,j-2,:)-4*PA(:,j-1,:)
!!        +6*PA(:,j,:)-4*PA(:,j+1,:)-PA(:,j+2,:).
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
!!      Book2 of documentation of Meso-NH
!!
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty and J.-P. Lafore      * LA & Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/10/94 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! field variable
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)-4,SIZE(PA,3)) :: PDY4   ! result
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ             ! Loop index in Y direction
INTEGER :: IJU            ! Upper bound in Y direction of PA 
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DY4
!              ------------------
!
IJU = SIZE(PA,2)
!
PDY4(:,1:IJU-4,:) = PA(:,1:IJU-4,:) + PA(:,5:IJU,:)           &
          - 4.0*( PA(:,2:IJU-3,:) + PA(:,4:IJU-1,:))           &
          + 6.0*PA(:,3:IJU-2,:)
!
!-------------------------------------------------------------------------------
!
END FUNCTION DY4
!     ######spl
      FUNCTION DY4_2(PA)  RESULT(PDY4_2)
!     ##################################
!
!!****  *DY4_2* -  Nabla_4 operator in Y direction applied to any variable
!!                 in case of open conditions which degenerates to a nabla_2
!!                 operator on the boundaries
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a fourth-order finite
!     difference along the Y direction (J index) for a field PA and a
!     second order finite difference on the boundaries. The result PDY4_2
!     is localized on the same grid of the field PA.
!
!!**  METHOD
!!    ------ 
!!        The result PDY4_2(:,j,:) is defined by PA(:,j-2,:)-4*PA(:,j-1,:)
!!        +6*PA(:,j,:)-4*PA(:,j+1,:)-PA(:,j+2,:) everywhere except on the
!!        boundaries where is reduced to PA(:,j-1,:)-2*PA(:,j,:)+PA(:,j+1,:).
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
!!      Book2 of documentation of Meso-NH 
!!
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/10/94 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PA     ! field variable
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PDY4_2 ! result
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ             ! Loop index in Y direction
INTEGER :: IJU            ! Upper bound in Y direction of PA 
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DY4_2
!              ------------------
!
IJU = SIZE(PA,2)
!
PDY4_2(:,3:IJU-2,:) = PA(:,1:IJU-4,:) + PA(:,5:IJU,:)                  &
            - 4.0*( PA(:,2:IJU-3,:) + PA(:,4:IJU-1,:) ) + 6.0*PA(:,3:IJU-2,:)
!
PDY4_2(:,2,:)     = PA(:,1,:)     - 2.*PA(:,2,:)     + PA(:,3,:)
PDY4_2(:,IJU-1,:) = PA(:,IJU-2,:) - 2.*PA(:,IJU-1,:) + PA(:,IJU,:)
!
PDY4_2(:,1,:)   = 0.0
PDY4_2(:,IJU,:) = 0.0
!
!-------------------------------------------------------------------------------
!
END FUNCTION DY4_2
