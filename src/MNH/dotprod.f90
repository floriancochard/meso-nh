!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 solver 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_DOTPROD
!     ###################
!
INTERFACE
!
      FUNCTION DOTPROD(PA,PB,HLBCX,HLBCY)  RESULT(PDOTPROD)
!  
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA, PB     ! input vectors
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type 
!
REAL :: PDOTPROD                                     ! dot product
!
END FUNCTION DOTPROD
!
END INTERFACE
!
END MODULE MODI_DOTPROD
!
!
!
!     #####################################################
      FUNCTION DOTPROD(PA,PB,HLBCX,HLBCY)  RESULT(PDOTPROD)
!     #####################################################
!
!!****  *DOTPROD* -  compute the dot product of two vectors
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute dot product of the vectors 
!     stored in the arrays PA, PB. The elements of PA and PB are localized at
!     mass points.
!
!!**  METHOD
!!    ------ 
!!      The scalar product DOTPROD  of 2 vectors  A and B is defined by : 
!!        DOTPROD = SUM( A(i,j,k)* B(i,j,k) )
!!    The bounds for the summation depend on the l.b.c.
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT, JPVEXT: define the number of marginal points out of the 
!!        physical domain along horizontal and vertical directions respectively
!!     Module MODD_CONF: model configurations
!!        L2D: logical switch for 2D model version
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (function DOTPROD)
!!
!!    AUTHOR
!!    ------
!!	P. Hereil and J. Stein      * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/07/94 
!!      J.-P. Pinty 12/11/99 Parallelization
!!
!-------------------------------------------------------------------------------
!
!*      0.     DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CONF
!
USE MODE_ll
!JUAN
USE MODE_REPRO_SUM
!JUAN
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments and result
!             ------------------------------------
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA, PB     ! input vectors
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type 
!
REAL :: PDOTPROD                                     ! dot product
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JI,JJ                         ! loop indexes
!
INTEGER :: IIB          ! indice I for the first inner mass point along x
INTEGER :: IIE          ! indice I for the last inner mass point along x
INTEGER :: IJB          ! indice J for the first inner mass point along y
INTEGER :: IJE          ! indice J for the last inner mass point along y
INTEGER :: IKB          ! indice K for the first inner mass point along z
INTEGER :: IKE          ! indice K for the last inner mass point along z
!
INTEGER :: ILBXB,ILBYB,ILBXE,ILBYE          ! loop indices depending on the
                                            ! lateral boundary conditions
!
INTEGER  :: IINFO_ll
!JUAN16
REAL, ALLOCATABLE, DIMENSION(:,:)     :: ZDOTPROD
!JUAN16
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE LOOP BOUNDS
!              -------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
IKB=1+JPVEXT
IKE=SIZE(PA,3) - JPVEXT
!
IF(HLBCX(1)/='CYCL' .AND. LWEST_ll()) THEN
  ILBXB = IIB-1    ! non cyclic condition at the physical boundary
ELSE
  ILBXB = IIB
ENDIF
!
IF(HLBCX(2)/='CYCL' .AND. LEAST_ll()) THEN
  ILBXE = IIE+1    ! non cyclic condition at the physical boundary
ELSE
  ILBXE = IIE
ENDIF
!
ILBYB = IJB
ILBYE = IJE
!
IF (.NOT.L2D) THEN        ! 3d version
  IF(HLBCY(1)/='CYCL' .AND. LSOUTH_ll()) THEN
    ILBYB = IJB-1    ! non cyclic condition at the physical boundary
  ELSE
    ILBYB = IJB
  ENDIF
!
  IF(HLBCY(2)/='CYCL' .AND. LNORTH_ll()) THEN
    ILBYE = IJE+1    ! non cyclic condition at the physical boundary
  ELSE
    ILBYE = IJE
  ENDIF
ELSE                      ! 2d version
  ILBYB = IJB
  ILBYE = IJB
ENDIF
!
!*       2.    COMPUTE THE DOT PRODUCT 
!              -----------------------
!
!JUAN16
ALLOCATE(ZDOTPROD(ILBXB:ILBXE,ILBYB:ILBYE))
ZDOTPROD    = 0.
DO JK = IKB-1,IKE+1
   DO JJ = ILBYB,ILBYE
      DO JI = ILBXB,ILBXE
         ZDOTPROD(JI,JJ) = ZDOTPROD(JI,JJ) + PA(JI,JJ,JK) * PB(JI,JJ,JK)
      END DO
   END DO
END DO
PDOTPROD = SUM_DD_R2_ll(ZDOTPROD)
!JUAN16
!
!-------------------------------------------------------------------------------
!
END FUNCTION DOTPROD
