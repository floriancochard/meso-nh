!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 surfex 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_UNPACK_1D_2D
!     ##########################
INTERFACE UNPACK_1D_2D
      SUBROUTINE UNPACK_1D_2D_FROMI2D(KM,K1D,K2D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
INTEGER, DIMENSION(:),   INTENT(IN) :: K1D
INTEGER, DIMENSION(:,:), INTENT(OUT):: K2D
!
END SUBROUTINE UNPACK_1D_2D_FROMI2D
!
      SUBROUTINE UNPACK_1D_2D_FROML2D(KM,O1D,O2D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
LOGICAL, DIMENSION(:),   INTENT(IN) :: O1D
LOGICAL, DIMENSION(:,:), INTENT(OUT):: O2D
!
END SUBROUTINE UNPACK_1D_2D_FROML2D
!
      SUBROUTINE UNPACK_1D_2D_FROM2D(KM,P1D,P2D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:),   INTENT(IN) :: P1D
REAL, DIMENSION(:,:), INTENT(OUT):: P2D
!
END SUBROUTINE UNPACK_1D_2D_FROM2D
!
      SUBROUTINE UNPACK_1D_2D_FROM3D(KM,P1D,P2D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:),   INTENT(IN) :: P1D
REAL, DIMENSION(:,:,:), INTENT(OUT):: P2D
!
END SUBROUTINE UNPACK_1D_2D_FROM3D
!
      SUBROUTINE UNPACK_1D_2D_FROM4D(KM,P1D,P2D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:,:),   INTENT(IN) :: P1D
REAL, DIMENSION(:,:,:,:), INTENT(OUT):: P2D
!
END SUBROUTINE UNPACK_1D_2D_FROM4D
!
END INTERFACE UNPACK_1D_2D
!
END MODULE MODI_UNPACK_1D_2D
!
!     ##############################################
      SUBROUTINE UNPACK_1D_2D_FROMI2D(KM,K1D,K2D)
!     ##############################################
!
!!****  *UNPACK_1D_2D* - put a 1D field into a 2D field
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/01/03
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,   ONLY : NUNDEF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
INTEGER, DIMENSION(:),   INTENT(IN) :: K1D
INTEGER, DIMENSION(:,:), INTENT(OUT):: K2D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER, DIMENSION(:), ALLOCATABLE :: I1D ! 1D field reshaped into a vector
!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(I1D(SIZE(K2D)))
!
I1D(:) = NUNDEF
!
DO JI=1,SIZE(K1D)
  I1D(KM(JI)) = K1D(JI) 
ENDDO
!
K2D=RESHAPE(I1D, (/ SIZE(K2D,1), SIZE(K2D,2) /) )
DEALLOCATE(I1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE UNPACK_1D_2D_FROMI2D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE UNPACK_1D_2D_FROML2D(KM,O1D,O2D)
!     ##############################################
!
!!****  *UNPACK_1D_2D* - put a 1D field into a 2D field
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/01/03
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,   ONLY : NUNDEF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
LOGICAL, DIMENSION(:),   INTENT(IN) :: O1D
LOGICAL, DIMENSION(:,:), INTENT(OUT):: O2D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
LOGICAL, DIMENSION(:), ALLOCATABLE :: G1D ! 1D field reshaped into a vector
!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(G1D(SIZE(O2D)))
!
G1D(:) = .FALSE.
!
DO JI=1,SIZE(O1D)
  G1D(KM(JI)) = O1D(JI) 
ENDDO
!
O2D=RESHAPE(G1D, (/ SIZE(O2D,1), SIZE(O2D,2) /) )
DEALLOCATE(G1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE UNPACK_1D_2D_FROML2D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE UNPACK_1D_2D_FROM2D(KM,P1D,P2D)
!     ##############################################
!
!!****  *UNPACK_1D_2D* - put a 1D field into a 2D field
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/01/03
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:),   INTENT(IN) :: P1D
REAL, DIMENSION(:,:), INTENT(OUT):: P2D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:), ALLOCATABLE :: Z1D ! 1D field reshaped into a vector
!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(Z1D(SIZE(P2D)))
!
Z1D(:) = XUNDEF
!
DO JI=1,SIZE(P1D)
  Z1D(KM(JI)) = P1D(JI) 
ENDDO
!
P2D=RESHAPE(Z1D, (/ SIZE(P2D,1), SIZE(P2D,2) /) )
DEALLOCATE(Z1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE UNPACK_1D_2D_FROM2D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE UNPACK_1D_2D_FROM3D(KM,P1D,P2D)
!     ##############################################
!
!!****  *UNPACK_1D_2D* - put a 1D field into a 2D field
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/01/03
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),     INTENT(IN) :: KM
REAL, DIMENSION(:,:),   INTENT(IN) :: P1D
REAL, DIMENSION(:,:,:), INTENT(OUT):: P2D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:,:), ALLOCATABLE :: Z1D ! 3D field reshaped (2 first components)
!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(Z1D(SIZE(P2D,1)*SIZE(P2D,2), SIZE(P2D,3)))
!
Z1D(:,:) = XUNDEF

DO JI=1,SIZE(P1D,1)
  Z1D(KM(JI),:) = P1D(JI,:)
ENDDO
!
P2D(:,:,:)=RESHAPE(Z1D, (/ SIZE(P2D,1), SIZE(P2D,2), SIZE(P2D,3) /) )
DEALLOCATE(Z1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE UNPACK_1D_2D_FROM3D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE UNPACK_1D_2D_FROM4D(KM,P1D,P2D)
!     ##############################################
!
!!****  *UNPACK_1D_2D* - extract the defined data from a 2D field into a 1D field
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/01/03
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,   ONLY : XUNDEF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),       INTENT(IN) :: KM
REAL, DIMENSION(:,:,:),   INTENT(IN) :: P1D
REAL, DIMENSION(:,:,:,:), INTENT(OUT):: P2D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: Z1D ! 3D field reshaped (3 first components)
!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(Z1D(SIZE(P2D,1)*SIZE(P2D,2), SIZE(P2D,3), SIZE(P2D,4)))
!
Z1D(:,:,:) = XUNDEF

DO JI=1,SIZE(P1D)
  Z1D(KM(JI),:,:) = P1D(JI,:,:)
ENDDO
!
P2D(:,:,:,:)=RESHAPE(Z1D, (/ SIZE(P2D,1), SIZE(P2D,2), SIZE(P2D,3), SIZE(P2D,4) /) )
DEALLOCATE(Z1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE UNPACK_1D_2D_FROM4D
