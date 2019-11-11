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
      MODULE MODI_PACK_2D_1D
!     ##########################
INTERFACE PACK_2D_1D
      SUBROUTINE PACK_2D_1D_FROMI2D(KM,K2D,K1D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
INTEGER, DIMENSION(:,:), INTENT(IN) :: K2D
INTEGER, DIMENSION(:),   INTENT(OUT):: K1D
!
END SUBROUTINE PACK_2D_1D_FROMI2D

      SUBROUTINE PACK_2D_1D_FROML2D(KM,K2D,K1D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
LOGICAL, DIMENSION(:,:), INTENT(IN) :: K2D
LOGICAL, DIMENSION(:),   INTENT(OUT):: K1D
!
END SUBROUTINE PACK_2D_1D_FROML2D
!
      SUBROUTINE PACK_2D_1D_FROM2D(KM,P2D,P1D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:), INTENT(IN) :: P2D
REAL, DIMENSION(:),   INTENT(OUT):: P1D
!
END SUBROUTINE PACK_2D_1D_FROM2D
!
      SUBROUTINE PACK_2D_1D_FROM3D(KM,P2D,P1D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:,:), INTENT(IN) :: P2D
REAL, DIMENSION(:,:),   INTENT(OUT):: P1D
!
END SUBROUTINE PACK_2D_1D_FROM3D
!
      SUBROUTINE PACK_2D_1D_FROM4D(KM,P2D,P1D)

INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: P2D
REAL, DIMENSION(:,:,:),   INTENT(OUT):: P1D
!
END SUBROUTINE PACK_2D_1D_FROM4D
!
END INTERFACE PACK_2D_1D
!
END MODULE MODI_PACK_2D_1D
!
!     ##############################################
      SUBROUTINE PACK_2D_1D_FROMI2D(KM,K2D,K1D)
!     ##############################################
!
!!****  *PACK_2D_1D* - extract the defined data from a 2D field into a 1D field
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
!!      J.Escobar   22/09/2012 :  PGI BUG on bound/reshape
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
INTEGER, DIMENSION(:,:), INTENT(IN) :: K2D
INTEGER, DIMENSION(:),   INTENT(OUT):: K1D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER, DIMENSION(:)  , ALLOCATABLE :: I1D ! 2D field reshaped into a vector
INTEGER, DIMENSION(:,:), ALLOCATABLE :: I2D ! PGI BUG on bound/reshape

!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(I1D(SIZE(K2D)))
ALLOCATE(I2D(SIZE(K2D,1),(SIZE(K2D,2))))
I2D = K2D
I1D=RESHAPE(I2D , (/ SIZE(K2D) /) )
!
DO JI=1,SIZE(K1D,1)
  K1D(JI) = I1D(KM(JI)) 
ENDDO
!
DEALLOCATE(I1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PACK_2D_1D_FROMI2D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE PACK_2D_1D_FROML2D(KM,O2D,O1D)
!     ##############################################
!
!!****  *PACK_2D_1D* - extract the defined data from a 2D field into a 1D field
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
!!      J.Escobar   22/09/2012 :  PGI BUG on bound/reshape
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
LOGICAL, DIMENSION(:,:), INTENT(IN) :: O2D
LOGICAL, DIMENSION(:),   INTENT(OUT):: O1D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
LOGICAL, DIMENSION(:), ALLOCATABLE :: G1D ! 2D field reshaped into a vector
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: G2D ! PGI BUG on bound/reshape
!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(G1D(SIZE(O2D)))
ALLOCATE(G2D(SIZE(O2D,1),SIZE(O2D,2)))
G2D = O2D
G1D=RESHAPE(G2D, (/ SIZE(O2D) /) )
!
DO JI=1,SIZE(O1D,1)
  O1D(JI) = G1D(KM(JI)) 
ENDDO
!
DEALLOCATE(G1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PACK_2D_1D_FROML2D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE PACK_2D_1D_FROM2D(KM,P2D,P1D)
!     ##############################################
!
!!****  *PACK_2D_1D* - extract the defined data from a 2D field into a 1D field
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
!!      J.Escobar   22/09/2012 :  PGI BUG on bound/reshape
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:), INTENT(IN) :: P2D
REAL, DIMENSION(:),   INTENT(OUT):: P1D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:)  , ALLOCATABLE :: Z1D ! 2D field reshaped into a vector
REAL, DIMENSION(:,:), ALLOCATABLE :: Z2D  ! PGI BUG on bound/reshape
!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(Z1D(SIZE(P2D)))
ALLOCATE(Z2D(SIZE(P2D,1),SIZE(P2D,2)))
Z2D = P2D
Z1D=RESHAPE(Z2D  , (/ SIZE(P2D) /) ) ! PGI BUG in bound => generate local  tmp array
DO JI=1,SIZE(P1D,1)
  P1D(JI) = Z1D(KM(JI)) 
ENDDO
!
DEALLOCATE(Z1D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PACK_2D_1D_FROM2D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE PACK_2D_1D_FROM3D(KM,P3D,P2D)
!     ##############################################
!
!!****  *PACK_2D_1D* - extract the defined data from a 3D field into a 2D field
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
!!      J.Escobar   22/09/2012 :  PGI BUG on bound/reshape
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:,:), INTENT(IN) :: P3D
REAL, DIMENSION(:,:),   INTENT(OUT):: P2D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:,:)  , ALLOCATABLE :: Z2D ! 3D field reshaped (2 first components)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: Z3D ! PGI BUG on bound/reshape

!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(Z2D(SIZE(P3D,1)*SIZE(P3D,2), SIZE(P3D,3)))
ALLOCATE(Z3D(SIZE(P3D,1),SIZE(P3D,2), SIZE(P3D,3)))
Z3D = P3D
Z2D(:,:)=RESHAPE(Z3D, (/ SIZE(P3D,1)*SIZE(P3D,2), SIZE(P3D,3) /) )
!
DO JI=1,SIZE(P2D,1)
  P2D(JI,:) = Z2D(KM(JI),:) 
ENDDO
!
DEALLOCATE(Z2D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PACK_2D_1D_FROM3D
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     ##############################################
      SUBROUTINE PACK_2D_1D_FROM4D(KM,P4D,P3D)
!     ##############################################
!
!!****  *PACK_4D_3D* - extract the defined data from a 4D field into a 3D field
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
!!      J.Escobar   22/09/2012 :  PGI BUG on bound/reshape
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, DIMENSION(:),   INTENT(IN) :: KM
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: P4D
REAL, DIMENSION(:,:,:),   INTENT(OUT):: P3D
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(:,:,:)  , ALLOCATABLE :: Z3D ! 3D field reshaped (3 first components)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: Z4D !  PGI BUG on bound/reshape

!
INTEGER :: JI ! loop counter
!
!-------------------------------------------------------------------------------
!
ALLOCATE(Z3D(SIZE(P4D,1)*SIZE(P4D,2), SIZE(P4D,3), SIZE(P4D,4)))
ALLOCATE(Z4D(SIZE(P4D,1),SIZE(P4D,2), SIZE(P4D,3), SIZE(P4D,4)))
Z4D = P4D
Z3D(:,:,:)=RESHAPE(Z4D, (/ SIZE(P4D,1)*SIZE(P4D,2), SIZE(P4D,3), SIZE(P4D,4) /) )
!
DO JI=1,SIZE(P3D,1)
  P3D(JI,:,:) = Z3D(KM(JI),:,:)
ENDDO
!
DEALLOCATE(Z3D)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PACK_2D_1D_FROM4D
