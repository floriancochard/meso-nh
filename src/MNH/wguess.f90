!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODI_WGUESS
!     ##################
INTERFACE
      SUBROUTINE WGUESS(PRHODJU,PRHODJV,PZZ,PDXX,PDYY,PDZZ,PDZX,PDZY,PRHODJW)
!
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODJU ! rhodJU on the MESO-NH grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODJV ! rhodJV on the MESO-NH grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZZ     ! height of w points
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZX    ! metric coefficient dzx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZY    ! metric coefficient dzy
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PRHODJW ! rhodJw on the MESO-NH grid
!
END SUBROUTINE WGUESS
END INTERFACE
END MODULE MODI_WGUESS
!
!     #######################################################################
      SUBROUTINE WGUESS(PRHODJU,PRHODJV,PZZ,PDXX,PDYY,PDZZ,PDZX,PDZY,PRHODJW)
!     #######################################################################
!
!!****  *WGUESS* - compute the first guess of w
!!
!!    PURPOSE
!!    -------
!!    This routine computes a value of w in order to have a first guess of w
!!    for the anelastic correction routine. The bottom boundary condition is
!!    verified, not the top one.
!!
!!**  METHOD
!!    ------
!!
!!    1 the vertical contravariant component Wc of the momentum is computed 
!!      from the uncompressible form of the continuity equation:
!!                     _                _                _
!!        d(rhodJ Uc)/dx + d(rhodJ Vc)/dy + d(rhodJ Wc)/dz = 0.
!!        with Wc=0 at ground level
!!
!!    2 the final value of w ( vertical catesian component) is deduced from:
!!
!!         Wc=1/dzz * (w-Udzx/dxx-Vdzy/dyy)
!!
!!      CAUTION: the values of rhoJw at JI=IJU-1, JJ=IJU-1,JK=IKU-1 or JK=IKB 
!!      are duplicated on the points JI=IJU, JJ=IJU JK=IKU and JK=IKB-1 to JK=1
!!      respectively.
!!
!!      
!!
!!    EXTERNAL
!!    --------
!!
!!    DXF,DYF,MXF,MYF,MZM             : Shuman operators
!!    Module MODI_SHUMAN              : interface for Shuman operators
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB   : verbosity level for output-listing
!!      Module MODD_PARAMETERS
!!         JPVEXT
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/12/94
!!                  15/05/96  spread the residual divergence on the whole domain
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_SHUMAN ! interface module
USE MODD_CONF   ! declaration modules
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODJU ! rhodJU on the MESO-NH grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRHODJV ! rhodJV on the MESO-NH grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZZ     ! height of w points
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZX    ! metric coefficient dzx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZY    ! metric coefficient dzy
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PRHODJW ! rhodJw on the MESO-NH grid
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER :: IIU,IJU
INTEGER :: IKB,IKU
INTEGER :: JK
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) ::  & !
                                         ZRHODJWC             ! rhoJ Wc
REAL, DIMENSION(SIZE(PDZZ,1),SIZE(PDZZ,2)) ::  ZLAMBDA 
                  ! characteristic length  for the weight function
!-------------------------------------------------------------------------------
!
IIU=SIZE(PDZZ,1)
IJU=SIZE(PDZZ,2)
IKB=JPVEXT+1
IKU=SIZE(PDZZ,3)
!
!*       1.    INTEGRATION OF THE CONTRAVARIANT W
!              ----------------------------------
!
ZRHODJWC(:,:,IKB)=0.
DO JK=IKB,IKU-1
  ZRHODJWC(:,:,JK+1:JK+1)=ZRHODJWC(:,:,JK:JK)                     &
                         -DXF(PRHODJU(:,:,JK:JK)/PDXX(:,:,JK:JK)) &
                         -DYF(PRHODJV(:,:,JK:JK)/PDYY(:,:,JK:JK))
END DO
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTATION OF rhoJw
!              --------------------
!
!*       2.1   General case
!              ------------
!
PRHODJW= PDZZ*ZRHODJWC + MXF(PDZX*MZM(1,IKU,1,PRHODJU/PDXX)) &
                       + MYF(PDZY*MZM(1,IKU,1,PRHODJV/PDYY))
!
!*       2.2   Copies on boundaries
!              --------------------
!
PRHODJW( 1 , : ,:)=PRHODJW(  2  ,  :  ,:)
PRHODJW(IIU, : ,:)=PRHODJW(IIU-1,  :  ,:)
PRHODJW( : , 1 ,:)=PRHODJW(  :  ,  2  ,:)
PRHODJW( : ,IJU,:)=PRHODJW(  :  ,IJU-1,:)
!
!*       2.3   Apply a weight function
!              -----------------------
ZLAMBDA(:,:)= (PZZ(:,:,IKU)-PZZ(:,:,IKB)) / 10
DO JK=IKB,IKU
  PRHODJW(:,:,JK) = PRHODJW(:,:,JK) *           &
     ( 1. - EXP( (PZZ(:,:,JK)-PZZ(:,:,IKU)) / ZLAMBDA(:,:) ) )
END DO
!
DO JK=1,IKB-1
 PRHODJW(:,:,JK)=PRHODJW(:,:,IKB)
END DO
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE WGUESS
