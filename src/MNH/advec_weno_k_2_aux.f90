!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##############################
      MODULE MODI_ADVEC_WENO_K_2_AUX
!     ##############################
!
INTERFACE
!
      SUBROUTINE ADVEC_WENO_K_2_UX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_WENO_K_2_UX
!
!                    ----------------------------
!
      SUBROUTINE ADVEC_WENO_K_2_MX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_WENO_K_2_MX
!
!                     ---------------------------
!
      SUBROUTINE ADVEC_WENO_K_2_VY(HLBCY,PSRC, PRVCT, PR, TPHALO2)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_WENO_K_2_VY
!
!                  ------------------------------
!
      SUBROUTINE ADVEC_WENO_K_2_MY(HLBCY,PSRC, PRVCT, PR, TPHALO2)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_WENO_K_2_MY
!
!                     -------------------------------
!
FUNCTION WENO_K_2_WZ(PSRC, PRWCT) RESULT(PR)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on W grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
END FUNCTION WENO_K_2_WZ
!
!                      ------------------------------
!
FUNCTION WENO_K_2_MZ(PSRC, PRWCT) RESULT(PR)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on W grid
!
! output source term
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
END FUNCTION WENO_K_2_MZ
!
END INTERFACE
!
END MODULE MODI_ADVEC_WENO_K_2_AUX
!
!-----------------------------------------------------------------------------
!
!     ############################################################
      SUBROUTINE ADVEC_WENO_K_2_UX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
!     ############################################################
!!
!!**** Computes PRUCT * PUT. Upstream fluxes of U in X direction.  
!!     Input PUT is on U Grid 'ie' (i,j,k) based on UGRID reference
!!     Output PR is on mass Grid 'ie' (i+1/2,j,k) based on UGRID reference
!!              
!!    AUTHOR
!!    ------
!!    F. Visentin   *CNRS/LA*               
!!
!!    MODIFICATIONS
!!    -------------
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
! output source term
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER:: IW,IE,IWF,IEF   ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFPOS1, ZFPOS2
!
! intermediate reconstruction fluxes for negative wind case
! we need only one since ZFNEG2 = ZFPOS2
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFNEG1, ZFNEG2
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBPOS1, ZBPOS2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBNEG1, ZBNEG2
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2
!
! standard weights
!
REAL, PARAMETER :: ZGAMMA1 = 1./3.
REAL, PARAMETER :: ZGAMMA2 = 2./3.
!
REAL, PARAMETER :: ZEPS = 1.0E-15
!
!-----------------------------------------------------------------------------
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!-------------------------------------------------------------------------------
!
!*       0.4.   INITIALIZE THE FIELD 
!               ---------------------
!
PR(:,:,:) = 0.0
!
ZFPOS1  = 0.0
ZFPOS2  = 0.0
ZFNEG1  = 0.0
ZFNEG2  = 0.0
ZBPOS1  = 0.0
ZBPOS2  = 0.0
ZBNEG1  = 0.0
ZBNEG2  = 0.0
ZOMP1   = 0.0
ZOMP2   = 0.0
ZOMN1   = 0.0
ZOMN2   = 0.0
!
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE X DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
  IW=IIB
  IE=IIE
!
! r: many left cells in regard to 'i' cell for each stencil
!
! intermediate fluxes at the mass point on Ugrid u(i+1/2,j,k) for positive wind
! (r=1 for the first stencil ZFPOS1, r=0 for the second ZFPOS2)
!
   ZFPOS1(IW:IE+1,:,:) = 0.5 * (3.0*PSRC(IW:IE+1,:,:) - PSRC(IW-1:IE,:,:))
   ZFPOS1(IW-1,   :,:) = 0.5 * (3.0*PSRC(IW-1,   :,:) - TPHALO2%WEST(:,:))
!
   ZFPOS2(IW-1:IE,:,:) = 0.5 * (PSRC(IW-1:IE,:,:) + PSRC(IW:IE+1,:,:))
   ZFPOS2(IE+1,   :,:) = 0.5 * (PSRC(IE+1,   :,:) + TPHALO2%EAST(:,:))
!
! intermediate flux at the mass point on Ugrid (i+1/2,j,k) for negative wind 
! case (from the right to the left)
! (r=0 for the second stencil ZFNEG2=ZFPOS2, r=-1 for the first ZFNEG1)  
!
  ZFNEG1(IW-1:IE-1,:,:) = 0.5 * (3.0*PSRC(IW:IE,:,:) - PSRC(IW+1:IE+1,:,:))
  ZFNEG1(IE,   :,:) = 0.5 * (3.0*PSRC(IE+1,   :,:) - TPHALO2%EAST(:,:))
  ZFNEG2(IW-1:IE,:,:) = 0.5 * (PSRC(IW-1:IE,:,:) + PSRC(IW:IE+1,:,:))
  ZFNEG2(IE+1,   :,:) = 0.5 * (PSRC(IE+1,   :,:) + TPHALO2%EAST(:,:))
!
! smoothness indicators for positive wind case
!
  ZBPOS1(IW:IE+1,:,:) = (PSRC(IW:IE+1,:,:) - PSRC(IW-1:IE,:,:))**2
  ZBPOS1(IW-1,   :,:) = (PSRC(IW-1,   :,:) - TPHALO2%WEST(:,:))**2
!
  ZBPOS2(IW-1:IE,:,:) = (PSRC(IW:IE+1,:,:) - PSRC(IW-1:IE,:,:))**2
  ZBPOS2(IE+1,   :,:) = (TPHALO2%EAST(:,:) - PSRC(IE+1,   :,:))**2
!
! smoothness indicators for negative wind case
!       
  ZBNEG1(IW-1:IE-1,:,:) = (PSRC(IW:IE,:,:)   - PSRC(IW+1:IE+1,:,:))**2
  ZBNEG1(IE,   :,:)     = (PSRC(IE+1,   :,:) - TPHALO2%EAST(:,:))**2
  ZBNEG2(IW-1:IE,:,:)   = (PSRC(IW-1:IE,:,:) - PSRC(IW:IE+1,:,:))**2
  ZBNEG2(IE+1,   :,:)   = (PSRC(IE+1,   :,:) - TPHALO2%EAST(:,:))**2
!
! WENO weights
!
  ZOMP1 = ZGAMMA1 / (ZEPS + ZBPOS1)**2
  ZOMP2 = ZGAMMA2 / (ZEPS + ZBPOS2)**2
  ZOMN1 = ZGAMMA1 / (ZEPS + ZBNEG1)**2
  ZOMN2 = ZGAMMA2 / (ZEPS + ZBNEG2)**2
!
! WENO fluxes
!
  PR = (ZOMN2/(ZOMN1+ZOMN2) * ZFNEG2 +                           &
       (ZOMN1/(ZOMN1+ZOMN2) * ZFNEG1)) * (0.5-SIGN(0.5,PRUCT)) + &
       (ZOMP2/(ZOMP1+ZOMP2) * ZFPOS2 +                           &
       (ZOMP1/(ZOMP1+ZOMP2) * ZFPOS1)) * (0.5+SIGN(0.5,PRUCT))
!
!
!       OPEN, WALL, NEST CASE IN THE X DIRECTION
!
CASE ('OPEN','WALL','NEST')
!
  IW=IIB
  IE=IIE
!
!       USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSICAL BORDER
!
  IF(LWEST_ll()) THEN
    PR(IW-1,:,:) = PSRC(IW-1,:,:) * (0.5+SIGN(0.5,PRUCT(IW-1,:,:))) + PSRC(IW,:,:) * (0.5-SIGN(0.5,PRUCT(IW-1,:,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE 
    ZFPOS1(IW-1,:,:) = 0.5 * (3.0*PSRC(IW-1,:,:) - TPHALO2%WEST(:,:))
    ZFPOS2(IW-1,:,:) = 0.5 * (PSRC(IW-1,    :,:) + PSRC(IW,:,:))
    ZBPOS1(IW-1,:,:) = (PSRC(IW-1,:,:) - TPHALO2%WEST(:,:))**2
    ZBPOS2(IW-1,:,:) = (PSRC(IW,  :,:) - PSRC(IW-1,:,:))**2
!
    ZFNEG1(IW-1,:,:) = 0.5 * (3.0*PSRC(IW,:,:) - PSRC(IW+1,:,:))
    ZFNEG2(IW-1,:,:) = 0.5 * (PSRC(IW-1,  :,:) + PSRC(IW,  :,:))
    ZBNEG1(IW-1,:,:) = (PSRC(IW,  :,:) - PSRC(IW+1,:,:))**2
    ZBNEG2(IW-1,:,:) = (PSRC(IW-1,:,:) - PSRC(IW,  :,:))**2
!
    ZOMP1(IW-1,:,:) = ZGAMMA1 / (ZEPS + ZBPOS1(IW-1,:,:))**2
    ZOMP2(IW-1,:,:) = ZGAMMA2 / (ZEPS + ZBPOS2(IW-1,:,:))**2
    ZOMN1(IW-1,:,:) = ZGAMMA1 / (ZEPS + ZBNEG1(IW-1,:,:))**2
    ZOMN2(IW-1,:,:) = ZGAMMA2 / (ZEPS + ZBNEG2(IW-1,:,:))**2
!
    PR(IW-1,:,:) = (ZOMN2(IW-1,:,:)/(ZOMN1(IW-1,:,:)+ZOMN2(IW-1,:,:)) * ZFNEG2(IW-1,:,:) +   &
                   (ZOMN1(IW-1,:,:)/(ZOMN1(IW-1,:,:)+ZOMN2(IW-1,:,:)) * ZFNEG1(IW-1,:,:))) * (0.5-SIGN(0.5,PRUCT(IW-1,:,:))) + &
                   (ZOMP2(IW-1,:,:)/(ZOMP1(IW-1,:,:)+ZOMP2(IW-1,:,:)) * ZFPOS2(IW-1,:,:) +  &
                   (ZOMP1(IW-1,:,:)/(ZOMP1(IW-1,:,:)+ZOMP2(IW-1,:,:)) * ZFPOS1(IW-1,:,:))) * (0.5+SIGN(0.5,PRUCT(IW-1,:,:)))
!
  ENDIF
!
  IF(LEAST_ll()) THEN
    PR(IE,:,:) = PSRC(IE,:,:) * (0.5+SIGN(0.5,PRUCT(IE,:,:))) + PSRC(IE+1,:,:) * (0.5-SIGN(0.5,PRUCT(IE,:,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    ZFPOS1(IE,:,:) = 0.5 * (3.0*PSRC(IE,:,:) - PSRC(IE-1,:,:))
    ZFPOS2(IE,:,:) = 0.5 * (PSRC(IE,    :,:) + PSRC(IE+1,:,:)) 
    ZBPOS1(IE,:,:) = (PSRC(IE,:,:) - PSRC(IE-1,:,:))**2
    ZBPOS2(IE,:,:) = (PSRC(IE+1,:,:) - PSRC(IE,:,:))**2
!
    ZFNEG1(IE,:,:) = 0.5 * (3.0*PSRC(IE+1,:,:) - TPHALO2%EAST(:,:))
    ZFNEG2(IE,:,:) = 0.5 * (PSRC(IE,:,:) + PSRC(IE+1,:,:))
    ZBNEG1(IE,:,:) = (PSRC(IE+1,:,:) - TPHALO2%EAST(:,:))**2
    ZBNEG2(IE,:,:) = (PSRC(IE,  :,:) - PSRC(IE+1,:,:))**2
!
    ZOMP1(IE,:,:) = ZGAMMA1 / (ZEPS + ZBPOS1(IE,:,:))**2
    ZOMP2(IE,:,:) = ZGAMMA2 / (ZEPS + ZBPOS2(IE,:,:))**2
    ZOMN1(IE,:,:) = ZGAMMA1 / (ZEPS + ZBNEG1(IE,:,:))**2
    ZOMN2(IE,:,:) = ZGAMMA2 / (ZEPS + ZBNEG2(IE,:,:))**2
!
    PR(IE,:,:) = (ZOMN2(IE,:,:)/(ZOMN1(IE,:,:)+ZOMN2(IE,:,:)) * ZFNEG2(IE,:,:) +  &
       (ZOMN1(IE,:,:)/(ZOMN1(IE,:,:)+ZOMN2(IE,:,:)) * ZFNEG1(IE,:,:))) * (0.5-SIGN(0.5,PRUCT(IE,:,:))) + &
                 (ZOMP2(IE,:,:)/(ZOMP1(IE,:,:)+ZOMP2(IE,:,:)) * ZFPOS2(IE,:,:) +  &
       (ZOMP1(IE,:,:)/(ZOMP1(IE,:,:)+ZOMP2(IE,:,:)) * ZFPOS1(IE,:,:))) * (0.5+SIGN(0.5,PRUCT(IE,:,:)))
!
  ENDIF
!
!      USE A THIRD ORDER UPSTREAM WENO SCHEME ELSEWHERE 
!
  ZFPOS1(IW:IE-1,:,:) = 0.5 * (3.0*PSRC(IW:IE-1,:,:) - PSRC(IW-1:IE-2,:,:))
  ZFPOS2(IW:IE-1,:,:) = 0.5 * (PSRC(IW:IE-1,    :,:) + PSRC(IW+1:IE,  :,:))
  ZBPOS1(IW:IE-1,:,:) = (PSRC(IW:IE-1,:,:) - PSRC(IW-1:IE-2,:,:))**2
  ZBPOS2(IW:IE-1,:,:) = (PSRC(IW+1:IE,:,:) - PSRC(IW:IE-1,  :,:))**2
!
  ZFNEG1(IW:IE-1,:,:) = 0.5 * (3.0*PSRC(IW+1:IE,:,:) - PSRC(IW+2:IE+1,:,:))
  ZFNEG2(IW:IE-1,:,:) = 0.5 * (PSRC(IW:IE-1,    :,:) + PSRC(IW+1:IE,  :,:))
  ZBNEG1(IW:IE-1,:,:) = (PSRC(IW+1:IE,:,:) - PSRC(IW+2:IE+1,:,:))**2
  ZBNEG2(IW:IE-1,:,:)   = (PSRC(IW:IE-1,:,:) - PSRC(IW+1:IE,:,:))**2
!
  ZOMP1(IW:IE-1,:,:) = ZGAMMA1 / (ZEPS + ZBPOS1(IW:IE-1,:,:))**2
  ZOMP2(IW:IE-1,:,:) = ZGAMMA2 / (ZEPS + ZBPOS2(IW:IE-1,:,:))**2
  ZOMN1(IW:IE-1,:,:) = ZGAMMA1 / (ZEPS + ZBNEG1(IW:IE-1,:,:))**2
  ZOMN2(IW:IE-1,:,:) = ZGAMMA2 / (ZEPS + ZBNEG2(IW:IE-1,:,:))**2
!
    PR(IW:IE-1,:,:) = (ZOMN2(IW:IE-1,:,:)/(ZOMN1(IW:IE-1,:,:)+ZOMN2(IW:IE-1,:,:)) * ZFNEG2(IW:IE-1,:,:) +   &
       (ZOMN1(IW:IE-1,:,:)/(ZOMN1(IW:IE-1,:,:)+ZOMN2(IW:IE-1,:,:)) * ZFNEG1(IW:IE-1,:,:))) * (0.5-SIGN(0.5,PRUCT(IW:IE-1,:,:))) + &
       (ZOMP2(IW:IE-1,:,:)/(ZOMP1(IW:IE-1,:,:)+ZOMP2(IW:IE-1,:,:)) * ZFPOS2(IW:IE-1,:,:) +   &
       (ZOMP1(IW:IE-1,:,:)/(ZOMP1(IW:IE-1,:,:)+ZOMP2(IW:IE-1,:,:)) * ZFPOS1(IW:IE-1,:,:))) * (0.5+SIGN(0.5,PRUCT(IW:IE-1,:,:)))
!
END SELECT
!
PR = PR * PRUCT
CALL GET_HALO(PR)
!
END SUBROUTINE ADVEC_WENO_K_2_UX
!
!------------------------------------------------------------------------------
!
!     ############################################################
      SUBROUTINE ADVEC_WENO_K_2_MX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
!     ############################################################
!!
!!**** Computes PRUCT * PWT (or PRUCT * PVT). Upstream fluxes of W (or V) 
!!     variables in X direction.  
!!     Input PWT is on W Grid 'ie' (i,j,k) based on WGRID reference
!!     Output PR is on mass Grid 'ie' (i-1/2,j,k) based on WGRID reference  
!!
!!    AUTHOR
!!    ------
!!    F. Visentin   *CNRS/LA*                
!!
!!    MODIFICATIONS
!!    -------------
!!
!------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
! output source term
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER::  IW,IE   ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFPOS1, ZFPOS2
!
! intermediate reconstruction fluxes for negative wind case
! we need only one since ZFNEG2 = ZFPOS2
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFNEG1, ZFNEG2
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBPOS1, ZBPOS2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBNEG1, ZBNEG2
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2
!
! standard weights
!
REAL, PARAMETER :: ZGAMMA1 = 1./3.
REAL, PARAMETER :: ZGAMMA2 = 2./3.
!
REAL, PARAMETER :: ZEPS = 1.0E-15
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!-----------------------------------------------------------------------------
!
!*       0.4.   INITIALIZE THE FIELD 
!               ---------------------
!
PR(:,:,:) = 0.0
!
ZFPOS1 = 0.0
ZFPOS2 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0 
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0 
!
!------------------------------------------------------------------------------
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE X DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
  IW=IIB
  IE=IIE
!
! intermediate fluxes for positive wind case
!
  ZFPOS1(IW+1:IE+1,:,:) = 0.5 * (3.0*PSRC(IW:IE,:,:) - PSRC(IW-1:IE-1,:,:))
  ZFPOS1(IW,       :,:) = 0.5 * (3.0*PSRC(IW-1, :,:) - TPHALO2%WEST(:,:))
!!  ZFPOS1(IW-1,     :,:) = - 999.
!
  ZFPOS2(IW:IE+1,:,:) = 0.5 * (PSRC(IW-1:IE,:,:) + PSRC(IW:IE+1,:,:))
  ZFPOS2(IW-1,   :,:) = 0.5 * (TPHALO2%WEST(:,:) + PSRC(IW-1,   :,:))
!
! intermediate flux for negative wind case
!
  ZFNEG1(IW-1:IE,:,:) = 0.5 * (3.0*PSRC(IW-1:IE,:,:) - PSRC(IW:IE+1,:,:))
  ZFNEG1(IE+1,   :,:) = 0.5 * (3.0*PSRC(IE+1,   :,:) - TPHALO2%EAST(:,:))
!
  ZFNEG2(IW:IE+1,:,:) = 0.5 * (PSRC(IW:IE+1,:,:) + PSRC(IW-1:IE,:,:))
  ZFNEG2(IW-1,       :,:) = 0.5 * (PSRC(IW-1, :,:) + TPHALO2%WEST(:,:))
! 
! smoothness indicators for positive wind case
!
  ZBPOS1(IW+1:IE+1,:,:) = (PSRC(IW:IE,:,:) - PSRC(IW-1:IE-1,:,:))**2
  ZBPOS1(IW,       :,:) = (PSRC(IW-1, :,:) - TPHALO2%WEST(:,:))**2
!!  ZBPOS1(IW-1,     :,:) = - 999.
!
  ZBPOS2(IW:IE+1,:,:) = (PSRC(IW:IE+1,:,:) - PSRC(IW-1:IE,:,:))**2
  ZBPOS2(IW-1,   :,:) = (PSRC(IW-1,   :,:) - TPHALO2%WEST(:,:))**2
!
! smoothness indicators for negative wind case
!       
  ZBNEG1(IW-1:IE,:,:) = (PSRC(IW-1:IE,:,:) - PSRC(IW:IE+1,:,:))**2
  ZBNEG1(IE+1,   :,:) = (PSRC(IE+1,   :,:) - TPHALO2%EAST(:,:))**2
!
  ZBNEG2(IW:IE+1,:,:) = (PSRC(IW-1:IE,:,:) - PSRC(IW:IE+1,:,:))**2
  ZBNEG2(IW-1,   :,:) = (TPHALO2%WEST(:,:) - PSRC(IW-1,:,:))**2
!
! WENO weights
!
  ZOMP1 = ZGAMMA1 / (ZEPS + ZBPOS1)**2
  ZOMP2 = ZGAMMA2 / (ZEPS + ZBPOS2)**2
  ZOMN1 = ZGAMMA1 / (ZEPS + ZBNEG1)**2
  ZOMN2 = ZGAMMA2 / (ZEPS + ZBNEG2)**2
!
! WENO fluxes
!
  PR = (ZOMP2/(ZOMP1+ZOMP2) * ZFPOS2 +                            &
       (ZOMP1/(ZOMP1+ZOMP2) * ZFPOS1)) * (0.5+SIGN(0.5,PRUCT )) + &
       (ZOMN2/(ZOMN1+ZOMN2) * ZFNEG2 +                            &
       (ZOMN1/(ZOMN1+ZOMN2) * ZFNEG1)) * (0.5-SIGN(0.5,PRUCT ))
!
!
!       OPEN, WALL, NEST CASE IN THE X DIRECTION
!
CASE ('OPEN','WALL','NEST')
!
  IW=IIB
  IE=IIE
!
!       USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSICAL BORDER
!
  IF(LWEST_ll()) THEN
    PR(IW,:,:) = PSRC(IW-1,:,:) * (0.5+SIGN(0.5,PRUCT(IW,:,:))) + PSRC(IW,:,:) * (0.5-SIGN(0.5,PRUCT(IW,:,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    ZFPOS1(IW,:,:) = 0.5 * (3.0*PSRC(IW-1, :,:) - TPHALO2%WEST(:,:))
    ZFPOS2(IW,:,:) = 0.5 * (PSRC(IW-1,     :,:) + PSRC(IW,     :,:))
    ZBPOS1(IW,:,:) = (PSRC(IW-1,:,:) - TPHALO2%WEST(:,:))**2
    ZBPOS2(IW,:,:) = (PSRC(IW,  :,:) - PSRC(IW-1,:,:))**2
!
    ZFNEG1(IW,:,:) = 0.5 * (3.0*PSRC(IW,:,:) - PSRC(IW+1,:,:))
    ZFNEG2(IW,:,:) = 0.5 * (PSRC(IW,    :,:) + PSRC(IW-1,:,:))
    ZBNEG1(IW,:,:) = (PSRC(IW,:,:) - PSRC(IW+1,:,:))**2
    ZBNEG2(IW,:,:) = (PSRC(IW-1,:,:) - PSRC(IW,:,:))**2
!
    ZOMP1(IW,:,:) = ZGAMMA1 / (ZEPS + ZBPOS1(IW,:,:))**2
    ZOMP2(IW,:,:) = ZGAMMA2 / (ZEPS + ZBPOS2(IW,:,:))**2
    ZOMN1(IW,:,:) = ZGAMMA1 / (ZEPS + ZBNEG1(IW,:,:))**2
    ZOMN2(IW,:,:) = ZGAMMA2 / (ZEPS + ZBNEG2(IW,:,:))**2
!
    PR(IW,:,:) = (ZOMP2(IW,:,:)/(ZOMP1(IW,:,:)+ZOMP2(IW,:,:)) * ZFPOS2(IW,:,:) +  &
       (ZOMP1(IW,:,:)/(ZOMP1(IW,:,:)+ZOMP2(IW,:,:)) * ZFPOS1(IW,:,:))) * (0.5+SIGN(0.5,PRUCT(IW,:,:))) + &
       (ZOMN2(IW,:,:)/(ZOMN1(IW,:,:)+ZOMN2(IW,:,:)) * ZFNEG2(IW,:,:) +   &
       (ZOMN1(IW,:,:)/(ZOMN1(IW,:,:)+ZOMN2(IW,:,:)) * ZFNEG1(IW,:,:))) * (0.5-SIGN(0.5,PRUCT(IW,:,:)))
!
  ENDIF
!
  IF(LEAST_ll()) THEN
    PR(IE+1,:,:) = PSRC(IE,:,:) * (0.5+SIGN(0.5,PRUCT(IE+1,:,:))) + PSRC(IE+1,:,:) * (0.5-SIGN(0.5,PRUCT(IE+1,:,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    ZFPOS1(IE+1,:,:) = 0.5 * (3.0*PSRC(IE,:,:) - PSRC(IE-1,:,:))
    ZFPOS2(IE+1,:,:) = 0.5 * (PSRC(IE,    :,:) + PSRC(IE+1,:,:))
    ZBPOS1(IE+1,:,:) = (PSRC(IE,:,:) - PSRC(IE-1,:,:))**2
    ZBPOS2(IE+1,:,:) = (PSRC(IE+1,:,:) - PSRC(IE,:,:))**2
!
    ZFNEG1(IE+1,:,:) = 0.5 * (3.0*PSRC(IE+1,:,:) - TPHALO2%EAST(:,:))
    ZFNEG2(IE+1,:,:) = 0.5 * (PSRC(IE+1,    :,:) + PSRC(IE,:,:))
    ZBNEG1(IE+1,:,:) = (PSRC(IE+1,:,:) - TPHALO2%EAST(:,:))**2
    ZBNEG2(IE+1,:,:) = (PSRC(IE,  :,:) - PSRC(IE+1,:,:))**2
!
    ZOMP1(IE+1,:,:) = ZGAMMA1 / (ZEPS + ZBPOS1(IE+1,:,:))**2
    ZOMP2(IE+1,:,:) = ZGAMMA2 / (ZEPS + ZBPOS2(IE+1,:,:))**2
    ZOMN1(IE+1,:,:) = ZGAMMA1 / (ZEPS + ZBNEG1(IE+1,:,:))**2
    ZOMN2(IE+1,:,:) = ZGAMMA2 / (ZEPS + ZBNEG2(IE+1,:,:))**2
!
    PR(IE+1,:,:) = (ZOMP2(IE+1,:,:)/(ZOMP1(IE+1,:,:)+ZOMP2(IE+1,:,:)) * ZFPOS2(IE+1,:,:) +  &
       (ZOMP1(IE+1,:,:)/(ZOMP1(IE+1,:,:)+ZOMP2(IE+1,:,:)) * ZFPOS1(IE+1,:,:))) * (0.5+SIGN(0.5,PRUCT(IE+1,:,:))) + &
       (ZOMN2(IE+1,:,:)/(ZOMN1(IE+1,:,:)+ZOMN2(IE+1,:,:)) * ZFNEG2(IE+1,:,:) +   &
       (ZOMN1(IE+1,:,:)/(ZOMN1(IE+1,:,:)+ZOMN2(IE+1,:,:)) * ZFNEG1(IE+1,:,:))) * (0.5-SIGN(0.5,PRUCT(IE+1,:,:)))
!
  ENDIF
!
!      USE A THIRD ORDER UPSTREAM WENO SCHEME ELSEWHERE 
!
  ZFPOS1(IW+1:IE,:,:) = 0.5 * (3.0*PSRC(IW:IE-1,:,:) - PSRC(IW-1:IE-2,:,:))
  ZFPOS2(IW+1:IE,:,:) = 0.5 * (PSRC(IW:IE-1,    :,:) + PSRC(IW+1:IE,  :,:))
  ZBPOS1(IW+1:IE,:,:) = (PSRC(IW:IE-1,:,:) - PSRC(IW-1:IE-2,:,:))**2
  ZBPOS2(IW+1:IE,:,:) = (PSRC(IW+1:IE,:,:) - PSRC(IW:IE-1,:,:))**2
!
  ZFNEG1(IW+1:IE,:,:) = 0.5 * (3.0*PSRC(IW+1:IE,:,:) - PSRC(IW+2:IE+1,:,:))
  ZFNEG2(IW+1:IE,:,:) = 0.5 * (PSRC(IW+1:IE,    :,:) + PSRC(IW:IE-1,  :,:))
  ZBNEG1(IW+1:IE,:,:) = (PSRC(IW+1:IE,:,:) - PSRC(IW+2:IE+1,:,:))**2
  ZBNEG2(IW+1:IE,:,:) = (PSRC(IW:IE-1,:,:) - PSRC(IW+1:IE,:,:))**2
!
  ZOMP1(IW+1:IE,:,:) = ZGAMMA1 / (ZEPS + ZBPOS1(IW+1:IE,:,:))**2
  ZOMP2(IW+1:IE,:,:) = ZGAMMA2 / (ZEPS + ZBPOS2(IW+1:IE,:,:))**2
  ZOMN1(IW+1:IE,:,:) = ZGAMMA1 / (ZEPS + ZBNEG1(IW+1:IE,:,:))**2
  ZOMN2(IW+1:IE,:,:) = ZGAMMA2 / (ZEPS + ZBNEG2(IW+1:IE,:,:))**2
!
  PR(IW+1:IE,:,:) = (ZOMP2(IW+1:IE,:,:)/(ZOMP1(IW+1:IE,:,:)+ZOMP2(IW+1:IE,:,:)) * ZFPOS2(IW+1:IE,:,:) + &
       (ZOMP1(IW+1:IE,:,:)/(ZOMP1(IW+1:IE,:,:)+ZOMP2(IW+1:IE,:,:)) * ZFPOS1(IW+1:IE,:,:))) * (0.5+SIGN(0.5,PRUCT(IW+1:IE,:,:))) + &
       (ZOMN2(IW+1:IE,:,:)/(ZOMN1(IW+1:IE,:,:)+ZOMN2(IW+1:IE,:,:)) * ZFNEG2(IW+1:IE,:,:) +              &
       (ZOMN1(IW+1:IE,:,:)/(ZOMN1(IW+1:IE,:,:)+ZOMN2(IW+1:IE,:,:)) * ZFNEG1(IW+1:IE,:,:))) * (0.5-SIGN(0.5,PRUCT(IW+1:IE,:,:)))
!
END SELECT
!
PR = PR * PRUCT
CALL GET_HALO(PR)
!
END SUBROUTINE ADVEC_WENO_K_2_MX
!
!-------------------------------------------------------------------------------
!
!     ############################################################
      SUBROUTINE ADVEC_WENO_K_2_MY(HLBCY,PSRC, PRVCT, PR, TPHALO2)
!     ############################################################
!!
!!****  Computes PRVCT * PUT (or PRVCT * PWT). Upstream fluxes of U (or W) 
!!      variables in Y direction.  
!!      Input PUT is on U Grid 'ie' (i,j,k) based on UGRID reference
!!      Output PR is on mass Grid 'ie' (i,j-1/2,k) based on UGRID reference 
!!
!!    AUTHOR
!!    ------
!!    F. Visentin   *CNRS/LA*                 
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
! output source term
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2 ! halo2 for the field at t
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER::  IS,IN      ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFPOS1, ZFPOS2
!
! intermediate reconstruction fluxes for negative wind case
! we need only one since ZFNEG2 = ZFPOS2
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFNEG1, ZFNEG2
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBPOS1, ZBPOS2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBNEG1, ZBNEG2
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2
!
! standard weights
!
REAL, PARAMETER :: ZGAMMA1 = 1./3.
REAL, PARAMETER :: ZGAMMA2 = 2./3.
!
REAL, PARAMETER :: ZEPS = 1.0E-15
!
!-----------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!---------------------------------------------------------------------------
!
!*       0.4.   INITIALIZE THE FIELD 
!               ---------------------
!
PR(:,:,:) = 0.0
!
ZFPOS1 = 0.0
ZFPOS2 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0 
!
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCY(1) ) ! 
!
!*       1.1    CYCLIC CASE IN THE Y DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCY(1) == HLBCY(2)
!
  IS=IJB
  IN=IJE
!
! intermediate fluxes for positive wind case
!
  ZFPOS1(:,IS+1:IN+1,:) = 0.5 * (3.0*PSRC(:,IS:IN,:) - PSRC(:,IS-1:IN-1,:))
  ZFPOS1(:,IS,       :) = 0.5 * (3.0*PSRC(:,IS-1, :) - TPHALO2%SOUTH(:,:))
!!  ZFPOS1(:,IS-1,     :) = - 999.
!
  ZFPOS2(:,IS:IN+1,:) = 0.5 * (PSRC(:,IS-1:IN,:) + PSRC(:,IS:IN+1,:))
  ZFPOS2(:,IS-1,   :) = 0.5 * (TPHALO2%SOUTH(:,:) + PSRC(:,IS-1,   :))
!
  ZFNEG1(:,IS-1:IN,:) = 0.5 * (3.0*PSRC(:,IS-1:IN,:) - PSRC(:,IS:IN+1,:))
  ZFNEG1(:,IN+1,   :) = 0.5 * (3.0*PSRC(:,IN+1,   :) - TPHALO2%NORTH(:,:))
!
  ZFNEG2(:,IS:IN+1,:) = 0.5 * (PSRC(:,IS:IN+1,:) + PSRC(:,IS-1:IN,:))
  ZFNEG2(:,IS-1,   :) = 0.5 * (PSRC(:,IS-1,   :) + TPHALO2%SOUTH(:,:))
!
! smoothness indicators for positive wind case
!
  ZBPOS1(:,IS+1:IN+1,:) = (PSRC(:,IS:IN,:) - PSRC(:,IS-1:IN-1,:))**2
  ZBPOS1(:,IS,       :) = (PSRC(:,IS-1,   :) - TPHALO2%SOUTH(:,:))**2
!!  ZBPOS1(:,IS-1,     :) = - 999. 
!
  ZBPOS2(:,IS:IN+1,:) = (PSRC(:,IS:IN+1,:) - PSRC(:,IS-1:IN,:))**2
  ZBPOS2(:,IS-1,   :) = (PSRC(:,IS-1,   :) - TPHALO2%SOUTH(:,:))**2
!
! smoothness indicators for negative wind case
!
  ZBNEG1(:,IS-1:IN,:) = (PSRC(:,IS-1:IN,:) - PSRC(:,IS:IN+1,:))**2
  ZBNEG1(:,IN+1,   :) = (PSRC(:,IN+1,   :) - TPHALO2%NORTH(:,:))**2
!
  ZBNEG2(:,IS:IN+1,:) = (PSRC(:,IS-1:IN,:) - PSRC(:,IS:IN+1,:))**2
  ZBNEG2(:,IS-1,   :) = (TPHALO2%SOUTH(:,:) - PSRC(:,IS-1,:))**2
!
! WENO weights
!
  ZOMP1 = ZGAMMA1 / (ZEPS + ZBPOS1)**2
  ZOMP2 = ZGAMMA2 / (ZEPS + ZBPOS2)**2
  ZOMN1 = ZGAMMA1 / (ZEPS + ZBNEG1)**2
  ZOMN2 = ZGAMMA2 / (ZEPS + ZBNEG2)**2
!
! WENO fluxes
!
  PR = (ZOMP2/(ZOMP1+ZOMP2) * ZFPOS2 +                           &
       (ZOMP1/(ZOMP1+ZOMP2) * ZFPOS1)) * (0.5+SIGN(0.5,PRVCT)) + &
       (ZOMN2/(ZOMN1+ZOMN2) * ZFNEG2 +                           &
       (ZOMN1/(ZOMN1+ZOMN2) * ZFNEG1)) * (0.5-SIGN(0.5,PRVCT))
!
!
!       OPEN, WALL, NEST CASE IN THE Y DIRECTION
!
CASE ('OPEN','WALL','NEST')
!
  IS=IJB
  IN=IJE
!
!       USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSICAL BORDER
!
  IF(LSOUTH_ll()) THEN
    PR(:,IS,:) = PSRC(:,IS-1,:) * (0.5+SIGN(0.5,PRVCT(:,IS,:))) + PSRC(:,IS,:) * (0.5-SIGN(0.5,PRVCT(:,IS,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    ZFPOS1(:,IS,:) = 0.5 * (3.0*PSRC(:,IS-1,:) - TPHALO2%SOUTH(:,:))
    ZFPOS2(:,IS,:) = 0.5 * (PSRC(:,IS-1,:) + PSRC(:,IS,:))
    ZBPOS1(:,IS,:) = (PSRC(:,IS-1,:) - TPHALO2%SOUTH(:,:))**2
    ZBPOS2(:,IS,:) = (PSRC(:,IS,  :) - PSRC(:,IS-1,:))**2
!
    ZFNEG1(:,IS,:) = 0.5 * (3.0*PSRC(:,IS,:) - PSRC(:,IS+1,:))
    ZFNEG2(:,IS,:) = 0.5 * (PSRC(:,IS,    :) + PSRC(:,IS-1,:))
    ZBNEG1(:,IS,:) = (PSRC(:,IS,  :) - PSRC(:,IS+1,:))**2
    ZBNEG2(:,IS,:) = (PSRC(:,IS-1,:) - PSRC(:,IS,  :))**2
!
    ZOMP1(:,IS,:) = ZGAMMA1 / (ZEPS + ZBPOS1(:,IS,:))**2
    ZOMP2(:,IS,:) = ZGAMMA2 / (ZEPS + ZBPOS2(:,IS,:))**2
    ZOMN1(:,IS,:) = ZGAMMA1 / (ZEPS + ZBNEG1(:,IS,:))**2
    ZOMN2(:,IS,:) = ZGAMMA2 / (ZEPS + ZBNEG2(:,IS,:))**2
!
    PR(:,IS,:) = (ZOMP2(:,IS,:)/(ZOMP1(:,IS,:)+ZOMP2(:,IS,:)) * ZFPOS2(:,IS,:) + &
       (ZOMP1(:,IS,:)/(ZOMP1(:,IS,:)+ZOMP2(:,IS,:)) * ZFPOS1(:,IS,:))) * (0.5+SIGN(0.5,PRVCT(:,IS,:))) + &
       (ZOMN2(:,IS,:)/(ZOMN1(:,IS,:)+ZOMN2(:,IS,:)) * ZFNEG2(:,IS,:) +  &
       (ZOMN1(:,IS,:)/(ZOMN1(:,IS,:)+ZOMN2(:,IS,:)) * ZFNEG1(:,IS,:))) * (0.5-SIGN(0.5,PRVCT(:,IS,:)))
!
  ENDIF
!
  IF(LNORTH_ll()) THEN
    PR(:,IN+1,:) = PSRC(:,IN,:) * (0.5+SIGN(0.5,PRVCT(:,IN+1,:))) + PSRC(:,IN+1,:) * (0.5-SIGN(0.5,PRVCT(:,IN+1,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    ZFPOS1(:,IN+1,:) = 0.5 * (3.0*PSRC(:,IN,:) - PSRC(:,IN-1,:))
    ZFPOS2(:,IN+1,:) = 0.5 * (PSRC(:,IN,    :) + PSRC(:,IN+1,:))
    ZBPOS1(:,IN+1,:) = (PSRC(:,IN,:) - PSRC(:,IN-1,:))**2
    ZBPOS2(:,IN+1,:) = (PSRC(:,IN+1,:) - PSRC(:,IN,:))**2
!
    ZFNEG1(:,IN+1,:) = 0.5 * (3.0*PSRC(:,IN+1,:) - TPHALO2%NORTH(:,:))
    ZFNEG2(:,IN+1,:) = 0.5 * (PSRC(:,IN+1,    :) + PSRC(:,IN,:))
    ZBNEG1(:,IN+1,:) = (PSRC(:,IN+1,:) - TPHALO2%NORTH(:,:))**2
    ZBNEG2(:,IN+1,:) = (PSRC(:,IN,  :) - PSRC(:,IN+1,:))**2
!
    ZOMP1(:,IN+1,:) = ZGAMMA1 / (ZEPS + ZBPOS1(:,IN+1,:))**2
    ZOMP2(:,IN+1,:) = ZGAMMA2 / (ZEPS + ZBPOS2(:,IN+1,:))**2
    ZOMN1(:,IN+1,:) = ZGAMMA1 / (ZEPS + ZBNEG1(:,IN+1,:))**2
    ZOMN2(:,IN+1,:) = ZGAMMA2 / (ZEPS + ZBNEG2(:,IN+1,:))**2
!
    PR(:,IN+1,:) = (ZOMP2(:,IN+1,:)/(ZOMP1(:,IN+1,:)+ZOMP2(:,IN+1,:)) * ZFPOS2(:,IN+1,:) +   &
       (ZOMP1(:,IN+1,:)/(ZOMP1(:,IN+1,:)+ZOMP2(:,IN+1,:)) * ZFPOS1(:,IN+1,:))) * (0.5+SIGN(0.5,PRVCT(:,IN+1,:))) + &
       (ZOMN2(:,IN+1,:)/(ZOMN1(:,IN+1,:)+ZOMN2(:,IN+1,:)) * ZFNEG2(:,IN+1,:) +     &
       (ZOMN1(:,IN+1,:)/(ZOMN1(:,IN+1,:)+ZOMN2(:,IN+1,:)) * ZFNEG1(:,IN+1,:))) * (0.5-SIGN(0.5,PRVCT(:,IN+1,:)))
!
  ENDIF
!
!      USE A THIRD ORDER UPSTREAM WENO SCHEME ELSEWHERE 
!
  ZFPOS1(:,IS+1:IN,:) = 0.5 * (3.0*PSRC(:,IS:IN-1,:) - PSRC(:,IS-1:IN-2,:))
  ZFPOS2(:,IS+1:IN,:) = 0.5 * (PSRC(:,IS:IN-1,    :) + PSRC(:,IS+1:IN,  :))
  ZBPOS1(:,IS+1:IN,:) = (PSRC(:,IS:IN-1,:) - PSRC(:,IS-1:IN-2,:))**2
  ZBPOS2(:,IS+1:IN,:) = (PSRC(:,IS+1:IN,:) - PSRC(:,IS:IN-1,  :))**2
!
  ZFNEG1(:,IS+1:IN,:) = 0.5 * (3.0*PSRC(:,IS+1:IN,:) - PSRC(:,IS+2:IN+1,:))
  ZFNEG2(:,IS+1:IN,:) = 0.5 * (PSRC(:,IS+1:IN,    :) + PSRC(:,IS:IN-1,  :))
  ZBNEG1(:,IS+1:IN,:) = (PSRC(:,IS+1:IN,:) - PSRC(:,IS+2:IN+1,:))**2
  ZBNEG2(:,IS+1:IN,:) = (PSRC(:,IS:IN-1,:) - PSRC(:,IS+1:IN,:))**2
!
  ZOMP1(:,IS+1:IN,:) = ZGAMMA1 / (ZEPS + ZBPOS1(:,IS+1:IN,:))**2
  ZOMP2(:,IS+1:IN,:) = ZGAMMA2 / (ZEPS + ZBPOS2(:,IS+1:IN,:))**2
  ZOMN1(:,IS+1:IN,:) = ZGAMMA1 / (ZEPS + ZBNEG1(:,IS+1:IN,:))**2
  ZOMN2(:,IS+1:IN,:) = ZGAMMA2 / (ZEPS + ZBNEG2(:,IS+1:IN,:))**2
!
  PR(:,IS+1:IN,:) = (ZOMP2(:,IS+1:IN,:)/(ZOMP1(:,IS+1:IN,:)+ZOMP2(:,IS+1:IN,:)) * ZFPOS2(:,IS+1:IN,:) +  &
       (ZOMP1(:,IS+1:IN,:)/(ZOMP1(:,IS+1:IN,:)+ZOMP2(:,IS+1:IN,:)) * ZFPOS1(:,IS+1:IN,:))) * (0.5+SIGN(0.5,PRVCT(:,IS+1:IN,:))) + &
       (ZOMN2(:,IS+1:IN,:)/(ZOMN1(:,IS+1:IN,:)+ZOMN2(:,IS+1:IN,:)) * ZFNEG2(:,IS+1:IN,:) + &
       (ZOMN1(:,IS+1:IN,:)/(ZOMN1(:,IS+1:IN,:)+ZOMN2(:,IS+1:IN,:)) * ZFNEG1(:,IS+1:IN,:))) * (0.5-SIGN(0.5,PRVCT(:,IS+1:IN,:)))
!
END SELECT
!
PR = PR * PRVCT
CALL GET_HALO(PR)
!
END SUBROUTINE ADVEC_WENO_K_2_MY
!-------------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE ADVEC_WENO_K_2_VY(HLBCY, PSRC, PRVCT, PR, TPHALO2)
!     #############################################################
!!
!!**** Computes PRVCT * PVT. Upstream fluxes of V in Y direction.  
!!     Input PVT is on V Grid 'ie' (i,j,k) based on VGRID reference
!!     Output PR is on mass Grid 'ie' (i,j+1/2,k) based on VGRID reference
!!
!!    AUTHOR
!!    ------
!!    F. Visentin   *CNRS/LA*
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER::  IS,IN      ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFPOS1, ZFPOS2
!
! intermediate reconstruction fluxes for negative wind case
! we need only one since ZFNEG2 = ZFPOS2
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFNEG1, ZFNEG2
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBPOS1, ZBPOS2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBNEG1, ZBNEG2
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2
!
! standard weights
!
REAL, PARAMETER :: ZGAMMA1 = 1./3.
REAL, PARAMETER :: ZGAMMA2 = 2./3.
!
REAL, PARAMETER :: ZEPS = 1.0E-15
!
!----------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!--------------------------------------------------------------------------
!
!*       0.4.   INITIALIZE THE FIELD 
!               ---------------------
!
PR(:,:,:) = 0.0
!
ZFPOS1 = 0.0
ZFPOS2 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0
!
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCY(1) ) ! Y direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE Y DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
  IS=IJB
  IN=IJE
!
! intermediate fluxes for positive wind case
!
  ZFPOS1(:,IS:IN+1,:) = 0.5 * (3.0*PSRC(:,IS:IN+1,:) - PSRC(:,IS-1:IN,:))
  ZFPOS1(:,IS-1,   :) = 0.5 * (3.0*PSRC(:,IS-1,   :) - TPHALO2%SOUTH(:,:))
!
  ZFPOS2(:,IS-1:IN,:) = 0.5 * (PSRC(:,IS-1:IN,:) + PSRC(:,IS:IN+1,:))
  ZFPOS2(:,IN+1,   :) = 0.5 * (PSRC(:,IN+1,   :) + TPHALO2%NORTH(:,:))
!
! intermediate flux for negative wind case
!
  ZFNEG1(:,IS-1:IN-1,:) = 0.5 * (3.0*PSRC(:,IS:IN,:) - PSRC(:,IS+1:IN+1,:))
  ZFNEG1(:,IN,   :) = 0.5 * (3.0*PSRC(:,IN+1,   :) - TPHALO2%NORTH(:,:))
!
  ZFNEG2(:,IS-1:IN,:) = 0.5 * (PSRC(:,IS-1:IN,:) + PSRC(:,IS:IN+1,:))
  ZFNEG2(:,IN+1,   :) = 0.5 * (PSRC(:,IN+1,   :) + TPHALO2%NORTH(:,:))
!
! smoothness indicators for positive wind case
!
  ZBPOS1(:,IS:IN+1,:) = (PSRC(:,IS:IN+1,:) - PSRC(:,IS-1:IN,:))**2
  ZBPOS1(:,IS-1,   :) = (PSRC(:,IS-1,   :) - TPHALO2%SOUTH(:,:))**2
!
  ZBPOS2(:,IS-1:IN,:) = (PSRC(:,IS:IN+1,:) - PSRC(:,IS-1:IN,:))**2
  ZBPOS2(:,IN+1,   :) = (TPHALO2%NORTH(:,:) - PSRC(:,IN+1,     :))**2
!
! smoothness indicators for negative wind case
!
  ZBNEG1(:,IS-1:IN-1,:) = (PSRC(:,IS:IN,:) - PSRC(:,IS+1:IN+1,:))**2
  ZBNEG1(:,IN,       :) = (PSRC(:,IN+1, :) - TPHALO2%NORTH(:,:))**2
!
  ZBNEG2(:,IS-1:IN,:) = (PSRC(:,IS-1:IN,:) - PSRC(:,IS:IN+1,:))**2
  ZBNEG2(:,IN+1,   :) = (PSRC(:,IN+1,   :) - TPHALO2%NORTH(:,:))**2 
!
! WENO weights
!
  ZOMP1 = ZGAMMA1 / (ZEPS + ZBPOS1)**2
  ZOMP2 = ZGAMMA2 / (ZEPS + ZBPOS2)**2
  ZOMN1 = ZGAMMA1 / (ZEPS + ZBNEG1)**2
  ZOMN2 = ZGAMMA2 / (ZEPS + ZBNEG2)**2
!
  PR = (ZOMP2/(ZOMP1+ZOMP2) * ZFPOS2 +                           &
       (ZOMP1/(ZOMP1+ZOMP2) * ZFPOS1)) * (0.5+SIGN(0.5,PRVCT)) + &
       (ZOMN2/(ZOMN1+ZOMN2) * ZFNEG2 +                           &
       (ZOMN1/(ZOMN1+ZOMN2) * ZFNEG1)) * (0.5-SIGN(0.5,PRVCT))
!
!
!       OPEN, WALL, NEST CASE IN THE Y DIRECTION
!
CASE ('OPEN','WALL','NEST')
!
  IS=IJB
  IN=IJE
!
!       USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSICAL BORDER
!
  IF(LSOUTH_ll()) THEN
    PR(:,IS-1,:) = PSRC(:,IS-1,:) * (0.5+SIGN(0.5,PRVCT(:,IS-1,:))) + PSRC(:,IS,:) * (0.5-SIGN(0.5,PRVCT(:,IS-1,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    ZFPOS1(:,IS-1,:) = 0.5 * (3.0*PSRC(:,IS-1,:) - TPHALO2%SOUTH(:,:))
    ZFPOS2(:,IS-1,:) = 0.5 * (PSRC(:,IS-1,    :) + PSRC(:,IS,:))
    ZBPOS1(:,IS-1,:) = (PSRC(:,IS-1,:) - TPHALO2%SOUTH(:,:))**2
    ZBPOS2(:,IS-1,:) = (PSRC(:,IS,  :) - PSRC(:,IS-1,:))**2
!
    ZFNEG1(:,IS-1,:) = 0.5 * (3.0*PSRC(:,IS,:) - PSRC(:,IS+1,:))
    ZFNEG2(:,IS-1,:) = 0.5 * (PSRC(:,IS-1,  :) + PSRC(:,IS,:))
    ZBNEG1(:,IS-1,:) = (PSRC(:,IS,:) - PSRC(:,IS+1,:))**2
    ZBNEG2(:,IS-1,:) = (PSRC(:,IS-1,:) - PSRC(:,IS,:))**2
!
    ZOMP1(:,IS-1,:) = ZGAMMA1 / (ZEPS + ZBPOS1(:,IS-1,:))**2
    ZOMP2(:,IS-1,:) = ZGAMMA2 / (ZEPS + ZBPOS2(:,IS-1,:))**2
    ZOMN1(:,IS-1,:) = ZGAMMA1 / (ZEPS + ZBNEG1(:,IS-1,:))**2
    ZOMN2(:,IS-1,:) = ZGAMMA2 / (ZEPS + ZBNEG2(:,IS-1,:))**2
!
    PR(:,IS-1,:) = (ZOMP2(:,IS-1,:)/(ZOMP1(:,IS-1,:)+ZOMP2(:,IS-1,:)) * ZFPOS2(:,IS-1,:) +   &
       (ZOMP1(:,IS-1,:)/(ZOMP1(:,IS-1,:)+ZOMP2(:,IS-1,:)) * ZFPOS1(:,IS-1,:))) * (0.5+SIGN(0.5,PRVCT(:,IS-1,:))) + &
       (ZOMN2(:,IS-1,:)/(ZOMN1(:,IS-1,:)+ZOMN2(:,IS-1,:)) * ZFNEG2(:,IS-1,:) +        &
       (ZOMN1(:,IS-1,:)/(ZOMN1(:,IS-1,:)+ZOMN2(:,IS-1,:)) * ZFNEG1(:,IS-1,:))) * (0.5-SIGN(0.5,PRVCT(:,IS-1,:)))
!
  ENDIF
!
  IF(LNORTH_ll()) THEN
    PR(:,IN,:) = PSRC(:,IN,:) * (0.5+SIGN(0.5,PRVCT(:,IN,:))) + PSRC(:,IN+1,:) * (0.5-SIGN(0.5,PRVCT(:,IN,:)))
!
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    ZFPOS1(:,IN,:) = 0.5 * (3.0*PSRC(:,IN,:) - PSRC(:,IN-1,:))
    ZFPOS2(:,IN,:) = 0.5 * (PSRC(:,IN,    :) + PSRC(:,IN+1,:))
    ZBPOS1(:,IN,:) = (PSRC(:,IN,  :) - PSRC(:,IN-1,:))**2
    ZBPOS2(:,IN,:) = (PSRC(:,IN+1,:) - PSRC(:,IN,  :))**2
!
    ZFNEG1(:,IN,:) = 0.5 * (3.0*PSRC(:,IN+1,:) - TPHALO2%NORTH(:,:))
    ZFNEG2(:,IN,:) = 0.5 * (PSRC(:,IN,      :) + PSRC(:,IN+1,:))
    ZBNEG1(:,IN,:) = (PSRC(:,IN+1,:) - TPHALO2%NORTH(:,:))**2
    ZBNEG2(:,IN,:) = (PSRC(:,IN,  :) - PSRC(:,IN+1,:))**2
!
    ZOMP1(:,IN,:) = ZGAMMA1 / (ZEPS + ZBPOS1(:,IN,:))**2
    ZOMP2(:,IN,:) = ZGAMMA2 / (ZEPS + ZBPOS2(:,IN,:))**2
    ZOMN1(:,IN,:) = ZGAMMA1 / (ZEPS + ZBNEG1(:,IN,:))**2
    ZOMN2(:,IN,:) = ZGAMMA2 / (ZEPS + ZBNEG2(:,IN,:))**2
!
    PR(:,IN,:) = (ZOMP2(:,IN,:)/(ZOMP1(:,IN,:)+ZOMP2(:,IN,:)) * ZFPOS2(:,IN,:) + &
       (ZOMP1(:,IN,:)/(ZOMP1(:,IN,:)+ZOMP2(:,IN,:)) * ZFPOS1(:,IN,:))) * (0.5+SIGN(0.5,PRVCT(:,IN,:))) + &
       (ZOMN2(:,IN,:)/(ZOMN1(:,IN,:)+ZOMN2(:,IN,:)) * ZFNEG2(:,IN,:) +  &
       (ZOMN1(:,IN,:)/(ZOMN1(:,IN,:)+ZOMN2(:,IN,:)) * ZFNEG1(:,IN,:))) * (0.5-SIGN(0.5,PRVCT(:,IN,:)))
!
  ENDIF
!
!      USE A THIRD ORDER UPSTREAM WENO SCHEME ELSEWHERE 
!
  ZFPOS1(:,IS:IN-1,:) = 0.5 * (3.0*PSRC(:,IS:IN-1,:) - PSRC(:,IS-1:IN-2,:))
  ZFPOS2(:,IS:IN-1,:) = 0.5 * (PSRC(:,IS:IN-1,    :) + PSRC(:,IS+1:IN,  :))
  ZBPOS1(:,IS:IN-1,:) = (PSRC(:,IS:IN-1,:) - PSRC(:,IS-1:IN-2,:))**2
  ZBPOS2(:,IS:IN-1,:) = (PSRC(:,IS+1:IN,:) - PSRC(:,IS:IN-1,  :))**2
!  
  ZFNEG1(:,IS:IN-1,:) = 0.5 * (3.0*PSRC(:,IS+1:IN,:) - PSRC(:,IS+2:IN+1,:))  
  ZFNEG2(:,IS:IN-1,:) = 0.5 * (PSRC(:,IS:IN-1,    :) + PSRC(:,IS+1:IN,  :))
  ZBNEG1(:,IS:IN-1,:) = (PSRC(:,IS+1:IN,:) - PSRC(:,IS+2:IN+1,:))**2
  ZBNEG2(:,IS:IN-1,:) = (PSRC(:,IS:IN-1,:) - PSRC(:,IS+1:IN,  :))**2
!
  ZOMP1(:,IS:IN-1,:) = ZGAMMA1 / (ZEPS + ZBPOS1(:,IS:IN-1,:))**2
  ZOMP2(:,IS:IN-1,:) = ZGAMMA2 / (ZEPS + ZBPOS2(:,IS:IN-1,:))**2
  ZOMN1(:,IS:IN-1,:) = ZGAMMA1 / (ZEPS + ZBNEG1(:,IS:IN-1,:))**2
  ZOMN2(:,IS:IN-1,:) = ZGAMMA2 / (ZEPS + ZBNEG2(:,IS:IN-1,:))**2
!
  PR(:,IS:IN-1,:) = (ZOMP2(:,IS:IN-1,:)/(ZOMP1(:,IS:IN-1,:)+ZOMP2(:,IS:IN-1,:)) * ZFPOS2(:,IS:IN-1,:) + &
       (ZOMP1(:,IS:IN-1,:)/(ZOMP1(:,IS:IN-1,:)+ZOMP2(:,IS:IN-1,:)) * ZFPOS1(:,IS:IN-1,:))) * (0.5+SIGN(0.5,PRVCT(:,IS:IN-1,:))) + &
       (ZOMN2(:,IS:IN-1,:)/(ZOMN1(:,IS:IN-1,:)+ZOMN2(:,IS:IN-1,:)) * ZFNEG2(:,IS:IN-1,:) + &
       (ZOMN1(:,IS:IN-1,:)/(ZOMN1(:,IS:IN-1,:)+ZOMN2(:,IS:IN-1,:)) * ZFNEG1(:,IS:IN-1,:))) * (0.5-SIGN(0.5,PRVCT(:,IS:IN-1,:)))
!
END SELECT
!
PR = PR * PRVCT
CALL GET_HALO(PR)
!
END SUBROUTINE ADVEC_WENO_K_2_VY
!
!-------------------------------------------------------------------------------
!
!     ############################################
      FUNCTION WENO_K_2_WZ(PSRC, PRWCT) RESULT(PR)
!     ############################################
!!
!!* Computes PRWCT * PWT. Upstream fluxes of W in Z direction.  
!!  Input PWT is on W Grid 'ie' (i,j,k) based on WGRID reference
!!  Output PR is on mass Grid 'ie' (i,j,k+1/2) based on WGRID reference
!!
!!    AUTHOR
!!    ------
!!    F. Visentin   *CNRS/LA*
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_CONF
USE MODD_PARAMETERS,ONLY: JPVEXT
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on W grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IB    ! Begining useful area in x,y,z directions
INTEGER :: IT    ! End useful area in x,y,z directions
!
! WENO-related variables:
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFPOS1, ZFPOS2
!
! intermediate reconstruction fluxes for negative wind case
! we need only one since ZFNEG2 = ZFPOS2
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFNEG1, ZFNEG2
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBPOS1, ZBPOS2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBNEG1, ZBNEG2
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2
!
! standard weights
!
REAL, PARAMETER :: ZGAMMA1 = 1./3.
REAL, PARAMETER :: ZGAMMA2 = 2./3.
!
REAL, PARAMETER :: ZEPS = 1.0E-15
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
IB = 1 + JPVEXT
IT = SIZE(PSRC,3) - JPVEXT
!
PR(:,:,:) = 0.0
!
ZFPOS1 = 0.0
ZFPOS2 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0 
!
! intermediate fluxes at the mass point on Wgrid w(i,j,k+1/2) for positive 
! wind case (L. to the R.)
! (r=1 for the first stencil ZFPOS1, r=0 for the second ZFPOS2)
!
ZFPOS1(:,:,IB:IT-1) = 0.5 * (3.0*PSRC(:,:,IB:IT-1) - PSRC(:,:,IB-1:IT-2))
ZFPOS2(:,:,IB:IT-1) = 0.5 * (PSRC(:,:,IB:IT-1) + PSRC(:,:,IB+1:IT))
!
! intermediate flux at the mass point on Wgrid w(i,j,k+1/2) for negative 
! wind case (R. to the L.)
! (r=-1 for the first stencil ZFNEG1, r=0 for the second ZFNEG2=ZFPOS2)
!
ZFNEG1(:,:,IB-1:IT-1) = 0.5 * (3.0*PSRC(:,:,IB:IT) - PSRC(:,:,IB+1:IT+1))
ZFNEG2(:,:,IB-1:IT) = 0.5 * (PSRC(:,:,IB-1:IT) + PSRC(:,:,IB:IT+1))
!
! smoothness indicators for positive wind case
!
ZBPOS1(:,:,IB:IT-1) = (PSRC(:,:,IB:IT-1) - PSRC(:,:,IB-1:IT-2))**2
ZBPOS2(:,:,IB:IT-1) = (PSRC(:,:,IB+1:IT) - PSRC(:,:,IB:IT-1))**2
!
! smoothness indicators for negative wind case
!
ZBNEG1(:,:,IB-1:IT-1) = (PSRC(:,:,IB:IT) - PSRC(:,:,IB+1:IT+1))**2
ZBNEG2(:,:,IB-1:IT) = (PSRC(:,:,IB-1:IT) - PSRC(:,:,IB:IT+1))**2
!
! WENO weights
!
ZOMP1 = ZGAMMA1 / (ZEPS + ZBPOS1)**2
ZOMP2 = ZGAMMA2 / (ZEPS + ZBPOS2)**2
ZOMN1 = ZGAMMA1 / (ZEPS + ZBNEG1)**2
ZOMN2 = ZGAMMA2 / (ZEPS + ZBNEG2)**2
!
! WENO fluxes
!
PR(:,:,IB:IT-1) = (ZOMP2(:,:,IB:IT-1)/(ZOMP1(:,:,IB:IT-1)+ZOMP2(:,:,IB:IT-1))* &
                                                         ZFPOS2(:,:,IB:IT-1) + &
                  (ZOMP1(:,:,IB:IT-1)/(ZOMP1(:,:,IB:IT-1)+ZOMP2(:,:,IB:IT-1))* &
                 ZFPOS1(:,:,IB:IT-1))) * (0.5+SIGN(0.5,PRWCT(:,:,IB:IT-1) )) + &
                  (ZOMN2(:,:,IB:IT-1)/(ZOMN1(:,:,IB:IT-1)+ZOMN2(:,:,IB:IT-1))* &
                                                         ZFNEG2(:,:,IB:IT-1) + &
                  (ZOMN1(:,:,IB:IT-1)/(ZOMN1(:,:,IB:IT-1)+ZOMN2(:,:,IB:IT-1))* &
                 ZFNEG1(:,:,IB:IT-1))) * (0.5-SIGN(0.5,PRWCT(:,:,IB:IT-1) ))
!
PR(:,:,IB-1) = PSRC(:,:,IB-1) * (0.5+SIGN(0.5,PRWCT(:,:,IB-1) )) + &
               PSRC(:,:,IB)   * (0.5-SIGN(0.5,PRWCT(:,:,IB-1) ))
PR(:,:,IT)   = PSRC(:,:,IT)   * (0.5+SIGN(0.5,PRWCT(:,:,IT) ))   + &
               PSRC(:,:,IT+1) * (0.5-SIGN(0.5,PRWCT(:,:,IT) ))
PR(:,:,IT+1) = -999.
!
PR = PR * PRWCT
CALL GET_HALO(PR)
!
END FUNCTION WENO_K_2_WZ
!
!-----------------------------------------------------------------------------
!
!     ############################################
      FUNCTION WENO_K_2_MZ(PSRC, PRWCT) RESULT(PR)
!     ############################################
!!
!!* Computes PRWCT * PUT (or PRWCT * PVT). Upstream fluxes of U (or V) 
!!  variables in Z direction.  
!!  Input PUT is on U Grid 'ie' (i,j,k) based on UGRID reference                
!!  Output PR is on mass Grid 'ie' (i,j,k-1/2) based on UGRID reference
!!
!!    AUTHOR
!!    ------
!!    F. Visentin   *CNRS/LA*
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_CONF
USE MODD_PARAMETERS,ONLY: JPVEXT
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on W grid
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IB    ! Begining useful area in x,y,z directions
INTEGER :: IT    ! End useful area in x,y,z directions
!
! WENO-related variables:
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFPOS1, ZFPOS2
!
! intermediate reconstruction fluxes for negative wind case
! we need only one since ZFNEG2 = ZFPOS2
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZFNEG1, ZFNEG2
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBPOS1, ZBPOS2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZBNEG1, ZBNEG2
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2
!
! standard weights
!
REAL, PARAMETER :: ZGAMMA1 = 1./3.
REAL, PARAMETER :: ZGAMMA2 = 2./3.
!
REAL, PARAMETER :: ZEPS = 1.0E-15
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
IB = 1 + JPVEXT
IT = SIZE(PSRC,3) - JPVEXT
!
PR(:,:,:) = 0.0
!
ZFPOS1 = 0.0
ZFPOS2 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0 
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0
!
! intermediate fluxes at the flux point on the Wgrid u(i,j,k-1/2) for 
! positive wind case
!
ZFPOS1(:,:,IB+1:IT) = 0.5 * (3.0*PSRC(:,:,IB:IT-1) - PSRC(:,:,IB-1:IT-2))
ZFPOS2(:,:,IB+1:IT) = 0.5 * (PSRC(:,:,IB:IT-1) + PSRC(:,:,IB+1:IT))
!
! intermediate flux at the flux point on the Wgrid u(i,j,k-1/2) for 
! negative wind case
!
ZFNEG1(:,:,IB+1:IT) = 0.5 * (3.0*PSRC(:,:,IB+1:IT) - PSRC(:,:,IB+2:IT+1))
ZFNEG2(:,:,IB+1:IT) = 0.5 * (PSRC(:,:,IB:IT-1) + PSRC(:,:,IB+1:IT))
!
! smoothness indicators for positive wind case
!
ZBPOS1(:,:,IB+1:IT) = (PSRC(:,:,IB:IT-1) - PSRC(:,:,IB-1:IT-2))**2
ZBPOS2(:,:,IB+1:IT) = (PSRC(:,:,IB+1:IT) - PSRC(:,:,IB:IT-1))**2
!
! smoothness indicators for negative wind case
!
ZBNEG1(:,:,IB+1:IT) = (PSRC(:,:,IB+1:IT) - PSRC(:,:,IB+2:IT+1))**2
ZBNEG2(:,:,IB+1:IT) = (PSRC(:,:,IB:IT-1) - PSRC(:,:,IB+1:IT))**2
!
! WENO weights
!
ZOMP1(:,:,IB+1:IT) = ZGAMMA1 / (ZEPS + ZBPOS1(:,:,IB+1:IT))**2
ZOMP2(:,:,IB+1:IT) = ZGAMMA2 / (ZEPS + ZBPOS2(:,:,IB+1:IT))**2
ZOMN1(:,:,IB+1:IT) = ZGAMMA1 / (ZEPS + ZBNEG1(:,:,IB+1:IT))**2
ZOMN2(:,:,IB+1:IT) = ZGAMMA2 / (ZEPS + ZBNEG2(:,:,IB+1:IT))**2
!
PR(:,:,IB+1:IT) = (ZOMP2(:,:,IB+1:IT)/(ZOMP1(:,:,IB+1:IT)+ZOMP2(:,:,IB+1:IT))* &
                                                         ZFPOS2(:,:,IB+1:IT) + &
                  (ZOMP1(:,:,IB+1:IT)/(ZOMP1(:,:,IB+1:IT)+ZOMP2(:,:,IB+1:IT))* &
                 ZFPOS1(:,:,IB+1:IT))) * (0.5+SIGN(0.5,PRWCT(:,:,IB+1:IT) )) + &
                  (ZOMN2(:,:,IB+1:IT)/(ZOMN1(:,:,IB+1:IT)+ZOMN2(:,:,IB+1:IT))* &
                                                         ZFNEG2(:,:,IB+1:IT) + &
                  (ZOMN1(:,:,IB+1:IT)/(ZOMN1(:,:,IB+1:IT)+ZOMN2(:,:,IB+1:IT))* &
                 ZFNEG1(:,:,IB+1:IT))) * (0.5-SIGN(0.5,PRWCT(:,:,IB+1:IT) ))
!
PR(:,:,IB)   = PSRC(:,:,IB-1) * (0.5+SIGN(0.5,PRWCT(:,:,IB) ))   + &
               PSRC(:,:,IB)   * (0.5-SIGN(0.5,PRWCT(:,:,IB) ))
PR(:,:,IT+1) = PSRC(:,:,IT)   * (0.5+SIGN(0.5,PRWCT(:,:,IT+1) )) + &
               PSRC(:,:,IT+1) * (0.5-SIGN(0.5,PRWCT(:,:,IT+1) ))
!
PR = PR * PRWCT
CALL GET_HALO(PR)
!
END FUNCTION WENO_K_2_MZ
