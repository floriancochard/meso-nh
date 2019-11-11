!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##############################
      MODULE MODI_ADVEC_WENO_K_3_AUX
!     ##############################
!
INTERFACE
!
      SUBROUTINE ADVEC_WENO_K_3_UX(HLBCX,PSRC, PRUCT, PR)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
!
END SUBROUTINE ADVEC_WENO_K_3_UX
!
!---------------------------------------------------------------------------------
!
      SUBROUTINE ADVEC_WENO_K_3_MX(HLBCX,PSRC, PRUCT, PR)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRUCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
!
END SUBROUTINE ADVEC_WENO_K_3_MX
!
!---------------------------------------------------------------------------------
!
      SUBROUTINE ADVEC_WENO_K_3_VY(HLBCY,PSRC, PRVCT, PR)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
!
END SUBROUTINE ADVEC_WENO_K_3_VY
!
!---------------------------------------------------------------------------------
!
      SUBROUTINE ADVEC_WENO_K_3_MY(HLBCY,PSRC, PRVCT, PR)
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY  ! Y direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on U grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRVCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
!
END SUBROUTINE ADVEC_WENO_K_3_MY
!
!---------------------------------------------------------------------------------
!
FUNCTION WENO_K_3_WZ(PSRC, PRWCT) RESULT(PR)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on W grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
END FUNCTION WENO_K_3_WZ
!
!---------------------------------------------------------------------------------
!
FUNCTION WENO_K_3_MZ(PSRC, PRWCT) RESULT(PR)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on W grid
!
! output source term
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
END FUNCTION WENO_K_3_MZ
!
END INTERFACE
!
END MODULE MODI_ADVEC_WENO_K_3_AUX
!
!-----------------------------------------------------------------------------
!
!     ############################################################
      SUBROUTINE ADVEC_WENO_K_3_UX(HLBCX,PSRC, PRUCT, PR)
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
!!    T. Lunet 02/10/2014:  Correction of periodic boudary conditions
!!       Change of structure in order to adapt WENO to NHALOK
!!       Suppression of second layer HALO pointers
!!       Complete code documentation
!!      J.Escobar : 25/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 02/10/2015 : correction on CYCL/OPEN boundaries
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
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
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER :: IW,IE      ! Physical boundary index
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFPOS1, ZFPOS2, ZFPOS3
!
! intermediate reconstruction fluxes for negative wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFNEG1, ZFNEG2, ZFNEG3
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBPOS1, ZBPOS2, ZBPOS3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBNEG1, ZBNEG2, ZBNEG3
!
! WENO non-normalized weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2, ZOMP3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2, ZOMN3
!
! EPSILON for weno weights calculation
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
!*       0.4.   INITIALIZE THE FIELD 
!               ---------------------
!
PR(:,:,:) = 0.0
!
ZFPOS1  = 0.0
ZFPOS2  = 0.0
ZFPOS3  = 0.0
ZFNEG1  = 0.0
ZFNEG2  = 0.0
ZFNEG3  = 0.0
ZBPOS1  = 0.0
ZBPOS2  = 0.0
ZBPOS3  = 0.0
ZBNEG1  = 0.0
ZBNEG2  = 0.0
ZBNEG3  = 0.0
ZOMP1   = 0.0
ZOMP2   = 0.0
ZOMP3   = 0.0
ZOMN1   = 0.0
ZOMN2   = 0.0
ZOMN3   = 0.0 
!
!-------------------------------------------------------------------------------
!*       1.1.   Interior Fluxes 
!              ---------------------
IW=IIB
IE=IIE
!
!-------------------------------------------------------------------------------
! Flux calculation in the physical domain far enough from the boundary 
! WENO scheme order 5, IW+1 -> IE-2
! Computation at the mass point on Ugrid u(i+1/2,j,k)
!-------------------------------------------------------------------------------
!
! ----- Positive fluxes -----
!
! First positive stencil, needs indices i-2, i-1, i
ZFPOS1(IW+1:IE-2,:,:) = 1./6.   * (2.0*PSRC(IW-1:IE-4,:,:) - 7.0*PSRC(IW:IE-3,:,:) + 11.0*PSRC(IW+1:IE-2,:,:)) ! Flux
ZBPOS1(IW+1:IE-2,:,:) = 13./12. * (    PSRC(IW-1:IE-4,:,:) - 2.0*PSRC(IW:IE-3,:,:) +      PSRC(IW+1:IE-2,:,:))**2 & 
                       + 1./4    * (    PSRC(IW-1:IE-4,:,:) - 4.0*PSRC(IW:IE-3,:,:) + 3.0* PSRC(IW+1:IE-2,:,:))**2 ! Smoothness indicator
ZOMP1(IW+1:IE-2,:,:)  = 1./10. /  (ZEPS + ZBPOS1(IW+1:IE-2,:,:))**2
!
! Second positive stencil, needs indices i-1, i, i+1
ZFPOS2(IW+1:IE-2,:,:) = 1./6.  * (-1.0*PSRC(IW:IE-3,:,:) + 5.0*PSRC(IW+1:IE-2,:,:) + 2.0*PSRC(IW+2:IE-1,:,:))! Flux
ZBPOS2(IW+1:IE-2,:,:) = 13./12 * (     PSRC(IW:IE-3,:,:) - 2.0*PSRC(IW+1:IE-2,:,:) +     PSRC(IW+2:IE-1,:,:))**2 &
                         + 1./4   * (  PSRC(IW:IE-3,:,:) - PSRC(IW+2:IE-1,:,:))**2! Smoothness indicator
ZOMP2(IW+1:IE-2,:,:)  = 3./5. /  (ZEPS + ZBPOS2(IW+1:IE-2,:,:))**2
!
! Third positive stencil, needs indices i, i+1, i+2
ZFPOS3(IW+1:IE-2,:,:) = 1./6   * (2.0*PSRC(IW+1:IE-2,:,:) + 5.0*PSRC(IW+2:IE-1,:,:) - PSRC(IW+3:IE,:,:))! Flux
ZBPOS3(IW+1:IE-2,:,:) = 13./12 * (PSRC(IW+1:IE-2,:,:) - 2.0*PSRC(IW+2:IE-1,:,:) + PSRC(IW+3:IE,:,:))**2 &
+ 1./4   * (3.0*PSRC(IW+1:IE-2,:,:) - 4.0*PSRC(IW+2:IE-1,:,:) + PSRC(IW+3:IE,:,:))**2! Smoothness indicator
ZOMP3(IW+1:IE-2,:,:)  = 3./10. / (ZEPS + ZBPOS3(IW+1:IE-2,:,:))**2
!
! ----- Negative fluxes ----- 
!
! First negative stencil, needs indices i+1, i+2, i+3
ZFNEG1(IW+1:IE-2,:,:) = 1./6.   * (11.0*PSRC(IW+2:IE-1,:,:) - 7.0*PSRC(IW+3:IE,:,:) + 2.0*PSRC(IW+4:IE+1,:,:)) ! Flux
ZBNEG1(IW+1:IE-2,:,:) = 13./12. * (     PSRC(IW+2:IE-1,:,:) - 2.0*PSRC(IW+3:IE,:,:) +     PSRC(IW+4:IE+1,:,:))**2 & 
                        + 1./4   * (3.0* PSRC(IW+2:IE-1,:,:) - 4.0*PSRC(IW+3:IE,:,:) +     PSRC(IW+4:IE+1,:,:))**2 ! Smoothness indicator
ZOMN1(IW+1:IE-2,:,:)  = 1./10. /  (ZEPS + ZBNEG1(IW+1:IE-2,:,:))**2
!
! Second negative stencil, needs indices i, i+1, i+2
ZFNEG2(IW+1:IE-2,:,:) = 1./6.  * (2.0*PSRC(IW+1:IE-2,:,:) + 5.0*PSRC(IW+2:IE-1,:,:) - 1.0*PSRC(IW+3:IE,:,:))! Flux
ZBNEG2(IW+1:IE-2,:,:) = 13./12 * (    PSRC(IW+1:IE-2,:,:) - 2.0*PSRC(IW+2:IE-1,:,:) +     PSRC(IW+3:IE,:,:))**2 &
                         + 1./4   * (    PSRC(IW+1:IE-2,:,:) -  PSRC(IW+3:IE,:,:))**2! Smoothness indicator
ZOMN2(IW+1:IE-2,:,:)  = 3./5. /  (ZEPS + ZBNEG2(IW+1:IE-2,:,:))**2! Non-normalized weight
!
! Third negative stencil, needs indices i-1, i, i+1
ZFNEG3(IW+1:IE-2,:,:) = 1./6   * (-1.0*PSRC(IW:IE-3,:,:) + 5.0*PSRC(IW+1:IE-2,:,:) + 2.0*PSRC(IW+2:IE-1,:,:))! Flux
ZBNEG3(IW+1:IE-2,:,:) = 13./12 * ( PSRC(IW:IE-3,:,:) - 2.0*PSRC(IW+1:IE-2,:,:) +    PSRC(IW+2:IE-1,:,:))**2 &
                          + 1./4   * (     PSRC(IW:IE-3,:,:) - 4.0*PSRC(IW+1:IE-2,:,:) + 3.0*PSRC(IW+2:IE-1,:,:))**2! Smoothness indicator
ZOMN3(IW+1:IE-2,:,:)  = 3./10. / (ZEPS + ZBNEG3(IW+1:IE-2,:,:))**2! Non-normalized weight
!
!
! ----- Total flux -----
!
PR(IW+1:IE-2,:,:) = (ZOMP1(IW+1:IE-2,:,:)/(ZOMP1(IW+1:IE-2,:,:)+ZOMP2(IW+1:IE-2,:,:)+ZOMP3(IW+1:IE-2,:,:)) &
                    * ZFPOS1(IW+1:IE-2,:,:) + &
                      ZOMP2(IW+1:IE-2,:,:)/(ZOMP1(IW+1:IE-2,:,:)+ZOMP2(IW+1:IE-2,:,:)+ZOMP3(IW+1:IE-2,:,:)) &
                    * ZFPOS2(IW+1:IE-2,:,:) + & 
                      ZOMP3(IW+1:IE-2,:,:)/(ZOMP1(IW+1:IE-2,:,:)+ZOMP2(IW+1:IE-2,:,:)+ZOMP3(IW+1:IE-2,:,:)) &
                    * ZFPOS3(IW+1:IE-2,:,:)) &
                    * (0.5+SIGN(0.5,PRUCT(IW+1:IE-2,:,:))) &
                  + (ZOMN1(IW+1:IE-2,:,:)/(ZOMN1(IW+1:IE-2,:,:)+ZOMN2(IW+1:IE-2,:,:)+ZOMN3(IW+1:IE-2,:,:)) &
                   * ZFNEG1(IW+1:IE-2,:,:)  &
                     + ZOMN2(IW+1:IE-2,:,:)/(ZOMN1(IW+1:IE-2,:,:)+ZOMN2(IW+1:IE-2,:,:)+ZOMN3(IW+1:IE-2,:,:)) &
                   * ZFNEG2(IW+1:IE-2,:,:)  &
                     + ZOMN3(IW+1:IE-2,:,:)/(ZOMN1(IW+1:IE-2,:,:)+ZOMN2(IW+1:IE-2,:,:)+ZOMN3(IW+1:IE-2,:,:)) &
                    * ZFNEG3(IW+1:IE-2,:,:))  & 
                    * (0.5-SIGN(0.5,PRUCT(IW+1:IE-2,:,:)))
!
!-------------------------------------------------------------------------------
!*       1.2.   West border
!               ---------------------
!
!!IF( LWEST_ll() .AND. .FALSE. ) THEN 
IF( LWEST_ll() ) THEN
!-----------------------------------------------------------------------------
! West border is physical -- IW,IW-1
!-----------------------------------------------------------------------------
SELECT CASE (HLBCX(1)) ! X direction LBC type on left side
! 
CASE ('CYCL')
!---------------------------------------------------------------------------
! Periodic boundary condition
!---------------------------------------------------------------------------
!
IF( LEAST_ll() .AND. .FALSE. ) THEN! East boundary is physical (monoproc)
!
! First positive stencil, needs indices i-2, i-1, i 
ZFPOS1(IW,:,:)  = 1./6. * (2.0*PSRC(IE,:,:)   - 7.0*PSRC(IW-1,:,:) + 11.0*PSRC(IW,:,:))! Flux IW
ZFPOS1(IW-1,:,:) = 1./6. * (2.0*PSRC(IE-1,:,:) - 7.0*PSRC(IE,:,:)   + 11.0*PSRC(IW-1,:,:))! Flux IW-1
ZBPOS1(IW,:,:)  = 13./12. * (PSRC(IE,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 & 
 + 1./4    * (PSRC(IE,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2 ! Smoothness indicator IW
ZBPOS1(IW-1,:,:) = 13./12. * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) +     PSRC(IW-1,:,:))**2 & 
 + 1./4    * (PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + 3.0*PSRC(IW-1,:,:))**2 ! Smoothness indicator IW-1
ZOMP1(IW-1:IW,:,:)  = 1./10. / (ZEPS + ZBPOS1(IW-1:IW,:,:))**2! Non-normalized weight IW,IW-1
!
! Second positive stencil, needs indices i-1, i, i+1
ZFPOS2(IW,:,:)   = 1./6.* (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:))! Flux IW
ZFPOS2(IW-1,:,:) = 1./6.* (-1.0*PSRC(IE,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW-1
ZBPOS2(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
 + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2! Smoothness indicator IW
ZBPOS2(IW-1,:,:) = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IW-1,:,:) + PSRC(IW,:,:))**2 &
 + 1./4   * (PSRC(IE,:,:) -                      PSRC(IW,:,:))**2! Smoothness indicator IW-1
ZOMP2(IW-1:IW,:,:)  = 3./5. / (ZEPS + ZBPOS2(IW-1:IW,:,:))**2! Non-normalized weight IW,IW-1
!
! Third negative stencil, needs indices i-1, i, i+1
ZFNEG3(IW,:,:)   = 1./6 * (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:))! Flux IW
ZFNEG3(IW-1,:,:) = 1./6 * (-1.0*PSRC(IE,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:))! Flux IW-1
ZBNEG3(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) +     PSRC(IW+1,:,:))**2 &
   + 1./4   * (PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + 3.0*PSRC(IW+1,:,:))**2! Smoothness indicator IW
ZBNEG3(IW-1,:,:) = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 &
   + 1./4   * (PSRC(IE,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2! Smoothness indicator IW-1
ZOMN3(IW-1:IW,:,:) = 3./10. / (ZEPS + ZBNEG3(IW-1:IW,:,:))**2! Non-normalized weight IW,IW-1
! 
ELSEIF(IW>3) THEN! East boundary is proc border, with minimum 3 HALO points on west side
!
! First positive stencil, needs indices i-2, i-1, i 
ZFPOS1(IW,:,:)   = 1./6. * (2.0*PSRC(IW-2,:,:)   - 7.0*PSRC(IW-1,:,:) + 11.0*PSRC(IW,:,:))! Flux IW
ZFPOS1(IW-1,:,:) = 1./6. * (2.0*PSRC(IW-3,:,:) - 7.0*PSRC(IW-2,:,:)   + 11.0*PSRC(IW-1,:,:))! Flux IW-1
ZBPOS1(IW,:,:)   = 13./12. * (PSRC(IW-1,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 & 
 + 1./4    * (PSRC(IW-1,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2 ! Smoothness indicator IW
ZBPOS1(IW-1,:,:) = 13./12. * (PSRC(IW-3,:,:) - 2.0*PSRC(IW-2,:,:) +     PSRC(IW-1,:,:))**2 & 
 + 1./4    * (PSRC(IW-3,:,:) - 4.0*PSRC(IW-2,:,:) + 3.0*PSRC(IW-1,:,:))**2 ! Smoothness indicator IW-1
ZOMP1(IW-1:IW,:,:)  = 1./10. / (ZEPS + ZBPOS1(IW-1:IW,:,:))**2! Non-normalized weight IW,IW-1
!
! Second positive stencil, needs indices i-1, i, i+1
ZFPOS2(IW,:,:)   = 1./6.* (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:))! Flux IW
ZFPOS2(IW-1,:,:) = 1./6.* (-1.0*PSRC(IW-2,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:))! Flux IW-1
ZBPOS2(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
 + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2! Smoothness indicator IW
ZBPOS2(IW-1,:,:) = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) + PSRC(IW,:,:))**2 &
 + 1./4   * (PSRC(IW-2,:,:) -                      PSRC(IW,:,:))**2! Smoothness indicator IW-1
ZOMP2(IW-1:IW,:,:)  = 3./5. / (ZEPS + ZBPOS2(IW-1:IW,:,:))**2! Non-normalized weight IW,IW-1
!
! Third negative stencil, needs indices i-1, i, i+1
ZFNEG3(IW,:,:)   = 1./6 * (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:))! Flux IW
 ZFNEG3(IW-1,:,:) = 1./6 * (-1.0*PSRC(IW-2,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW-1
 ZBNEG3(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) +     PSRC(IW+1,:,:))**2 &
        + 1./4   * (PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + 3.0*PSRC(IW+1,:,:))**2 ! Smoothness indicator IW
 ZBNEG3(IW-1,:,:) = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 &
        + 1./4   * (PSRC(IW-2,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2 ! Smoothness indicator IW-1
 ZOMN3(IW-1:IW,:,:) = 3./10. / (ZEPS + ZBNEG3(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ELSE ! East boundary is proc border, with NHALO < 3 on west side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on west side'
  CALL ABORT  
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
!
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IW,:,:)   = 1./6 * (2.0*PSRC(IW,:,:) + 5.0*PSRC(IW+1,:,:) - PSRC(IW+2,:,:)) ! Flux IW
 ZFPOS3(IW-1,:,:) = 1./6 * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:) - PSRC(IW+1,:,:)) ! Flux IW-1
 ZBPOS3(IW,:,:)   = 13./12 * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 ! Smoothness indicator IW
 ZBPOS3(IW-1,:,:) = 13./12 * (    PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 ! Smoothness indicator IW-1
 ZOMP3(IW-1:IW,:,:)  = 3./10. / (ZEPS + ZBPOS3(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1 
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(IW,:,:)   = 1./6. * (11.0*PSRC(IW+1,:,:) - 7.0*PSRC(IW+2,:,:) + 2.0*PSRC(IW+3,:,:)) ! Flux IW
 ZFNEG1(IW-1,:,:) = 1./6. * (11.0*PSRC(IW,:,:)   - 7.0*PSRC(IW+1,:,:) + 2.0*PSRC(IW+2,:,:)) ! Flux IW-1
 ZBNEG1(IW,:,:)   = 13./12. * (    PSRC(IW+1,:,:) - 2.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW+1,:,:) - 4.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2  ! Smoothness indicator IW
 ZBNEG1(IW-1,:,:) = 13./12. * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2  ! Smoothness indicator IW-1
 ZOMN1(IW-1:IW,:,:) = 1./10. / (ZEPS + ZBNEG1(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(IW,:,:)   = 1./6. * (2.0*PSRC(IW,:,:)   + 5.0*PSRC(IW+1,:,:) - 1.0*PSRC(IW+2,:,:)) ! Flux IW
 ZFNEG2(IW-1,:,:) = 1./6. * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   - 1.0*PSRC(IW+1,:,:)) ! Flux IW-1
 ZBNEG2(IW,:,:)   = 13./12 * (PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (PSRC(IW,:,:) -                      PSRC(IW+2,:,:))**2 ! Smoothness indicator IW
 ZBNEG2(IW-1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW-1
 ZOMN2(IW-1:IW,:,:) = 3./5. / (ZEPS + ZBNEG2(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ! ----- Total flux -----
! 
 PR(IW-1:IW,:,:) = ( ZOMP1(IW-1:IW,:,:)/(ZOMP1(IW-1:IW,:,:)+ZOMP2(IW-1:IW,:,:)+ZOMP3(IW-1:IW,:,:)) * ZFPOS1(IW-1:IW,:,:) &
                       + ZOMP2(IW-1:IW,:,:)/(ZOMP1(IW-1:IW,:,:)+ZOMP2(IW-1:IW,:,:)+ZOMP3(IW-1:IW,:,:)) * ZFPOS2(IW-1:IW,:,:) & 
                       + ZOMP3(IW-1:IW,:,:)/(ZOMP1(IW-1:IW,:,:)+ZOMP2(IW-1:IW,:,:)+ZOMP3(IW-1:IW,:,:)) * ZFPOS3(IW-1:IW,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IW-1:IW,:,:))) &
                    + ( ZOMN1(IW-1:IW,:,:)/(ZOMN1(IW-1:IW,:,:)+ZOMN2(IW-1:IW,:,:)+ZOMN3(IW-1:IW,:,:)) * ZFNEG1(IW-1:IW,:,:) &
                       + ZOMN2(IW-1:IW,:,:)/(ZOMN1(IW-1:IW,:,:)+ZOMN2(IW-1:IW,:,:)+ZOMN3(IW-1:IW,:,:)) * ZFNEG2(IW-1:IW,:,:) &
                       + ZOMN3(IW-1:IW,:,:)/(ZOMN1(IW-1:IW,:,:)+ZOMN2(IW-1:IW,:,:)+ZOMN3(IW-1:IW,:,:)) * ZFNEG3(IW-1:IW,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IW-1:IW,:,:)))
!
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IW-1
    PR(IW-1,:,:) = PSRC(IW-1,:,:) * (0.5+SIGN(0.5,PRUCT(IW-1,:,:))) + &
                   PSRC(IW,:,:) * &
                   (0.5-SIGN(0.5,PRUCT(IW-1,:,:)))
!
!   ! WENO scheme order 3, IW
    ZFPOS1(IW,:,:) = 0.5 * (3.0*PSRC(IW,:,:) - PSRC(IW-1,:,:)) ! First positive flux
    ZFPOS2(IW,:,:) = 0.5 * ( PSRC(IW,:,:) + PSRC(IW+1,:,:)) ! Second positive flux
    ZBPOS1(IW,:,:) = (PSRC(IW,:,:)   - PSRC(IW-1,:,:))**2 ! First positive smoothness indicator
    ZBPOS2(IW,:,:) = (PSRC(IW+1,:,:) - PSRC(IW,:,:))**2  ! Second positive smoothness indicator
!
    ZFNEG1(IW,:,:) = 0.5 * (3.0*PSRC(IW+1,:,:) - PSRC(IW+2,:,:)) ! First negative flux
    ZFNEG2(IW,:,:) = 0.5 * ( PSRC(IW,:,:)   + PSRC(IW+1,:,:)) ! Second negative flux
    ZBNEG1(IW,:,:) = (PSRC(IW+1,:,:) - PSRC(IW+2,:,:))**2 ! First negative smoothness indicator
    ZBNEG2(IW,:,:) = (PSRC(IW,:,:)   - PSRC(IW+1,:,:))**2 ! Second negative smoothness indicator
!
    ZOMP1(IW,:,:) = 1./3. / (ZEPS + ZBPOS1(IW,:,:))**2 ! First positive non-normalized weight
    ZOMP2(IW,:,:) = 2./3. / (ZEPS + ZBPOS2(IW,:,:))**2 ! Second positive non-normalized weight
    ZOMN1(IW,:,:) = 1./3. / (ZEPS + ZBNEG1(IW,:,:))**2 ! First negative non-normalized weight
    ZOMN2(IW,:,:) = 2./3. / (ZEPS + ZBNEG2(IW,:,:))**2 ! Second negative non-normalized weight
! 
    PR(IW,:,:) = (ZOMN2(IW,:,:)/(ZOMN1(IW,:,:)+ZOMN2(IW,:,:)) * ZFNEG2(IW,:,:) + &
             (ZOMN1(IW,:,:)/(ZOMN1(IW,:,:)+ZOMN2(IW,:,:)) * ZFNEG1(IW,:,:))) &
           *(0.5-SIGN(0.5,PRUCT(IW,:,:))) + &
     (ZOMP2(IW,:,:)/(ZOMP1(IW,:,:)+ZOMP2(IW,:,:)) * ZFPOS2(IW,:,:) + &
     (ZOMP1(IW,:,:)/(ZOMP1(IW,:,:)+ZOMP2(IW,:,:)) * ZFPOS1(IW,:,:))) &
     *(0.5+SIGN(0.5,PRUCT(IW,:,:)))  ! Total flux
! 
 END SELECT ! SELECT CASE (HLBCX(1)) ! X direction LBC type on left side
!
ELSE
 !-----------------------------------------------------------------------------
 ! West border is proc border -- IW,IW-1
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/west-int not parallelisable with NHALO < 3' 
 CALL ABORT 
 STOP ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
!
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(IW,:,:)   = 1./6. * (2.0*PSRC(IW-2,:,:) - 7.0*PSRC(IW-1,:,:) + 11.0*PSRC(IW,:,:)) ! Flux IW
 ZFPOS1(IW-1,:,:) = 1./6. * (2.0*PSRC(IW-3,:,:) - 7.0*PSRC(IW-2,:,:) + 11.0*PSRC(IW-1,:,:)) ! Flux IW-1
 ZBPOS1(IW,:,:)   = 13./12. * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 & 
        + 1./4    * (PSRC(IW-2,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2   ! Smoothness indicator IW
 ZBPOS1(IW-1,:,:) = 13./12. * (PSRC(IW-3,:,:) - 2.0*PSRC(IW-2,:,:) +     PSRC(IW-1,:,:))**2 & 
        + 1./4    * (PSRC(IW-3,:,:) - 4.0*PSRC(IW-2,:,:) + 3.0*PSRC(IW-1,:,:))**2  ! Smoothness indicator IW-1
 ZOMP1(IW-1:IW,:,:)  = 1./10. / (ZEPS + ZBPOS1(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(IW,:,:)   = 1./6.* (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW
 ZFPOS2(IW-1,:,:) = 1./6.* (-1.0*PSRC(IW-2,:,:) + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW-1
 ZBPOS2(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW
 ZBPOS2(IW-1,:,:) = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (PSRC(IW-2,:,:) -                      PSRC(IW,:,:))**2 ! Smoothness indicator IW-1
 ZOMP2(IW-1:IW,:,:)  = 3./5. / (ZEPS + ZBPOS2(IW-1:IW,:,:))**2  ! Non-normalized weight IW,IW-1
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IW,:,:)   = 1./6 * (2.0*PSRC(IW,:,:) + 5.0*PSRC(IW+1,:,:) - PSRC(IW+2,:,:)) ! Flux IW
 ZFPOS3(IW-1,:,:) = 1./6 * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:) - PSRC(IW+1,:,:)) ! Flux IW-1
 ZBPOS3(IW,:,:)   = 13./12 * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 ! Smoothness indicator IW
 ZBPOS3(IW-1,:,:) = 13./12 * (    PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 ! Smoothness indicator IW-1
 ZOMP3(IW-1:IW,:,:)  = 3./10. / (ZEPS + ZBPOS3(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(IW,:,:)   = 1./6. * (11.0*PSRC(IW+1,:,:) - 7.0*PSRC(IW+2,:,:) + 2.0*PSRC(IW+3,:,:)) ! Flux IW
 ZFNEG1(IW-1,:,:) = 1./6. * (11.0*PSRC(IW,:,:)   - 7.0*PSRC(IW+1,:,:) + 2.0*PSRC(IW+2,:,:)) ! Flux IW-1
 ZBNEG1(IW,:,:)   = 13./12. * (    PSRC(IW+1,:,:) - 2.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW+1,:,:) - 4.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2  ! Smoothness indicator IW
 ZBNEG1(IW-1,:,:) = 13./12. * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2  ! Smoothness indicator IW-1
 ZOMN1(IW-1:IW,:,:) = 1./10. / (ZEPS + ZBNEG1(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(IW,:,:)   = 1./6. * (2.0*PSRC(IW,:,:)   + 5.0*PSRC(IW+1,:,:) - 1.0*PSRC(IW+2,:,:)) ! Flux IW
 ZFNEG2(IW-1,:,:) = 1./6. * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   - 1.0*PSRC(IW+1,:,:)) ! Flux IW-1
 ZBNEG2(IW,:,:)   = 13./12 * (PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (PSRC(IW,:,:) -                      PSRC(IW+2,:,:))**2 ! Smoothness indicator IW
 ZBNEG2(IW-1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW-1
 ZOMN2(IW-1:IW,:,:) = 3./5. / (ZEPS + ZBNEG2(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(IW,:,:)   = 1./6 * (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW
 ZFNEG3(IW-1,:,:) = 1./6 * (-1.0*PSRC(IW-2,:,:) + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW-1
 ZBNEG3(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) +     PSRC(IW+1,:,:))**2 &
        + 1./4   * (PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + 3.0*PSRC(IW+1,:,:))**2 ! Smoothness indicator IW
 ZBNEG3(IW-1,:,:) = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 &
        + 1./4   * (PSRC(IW-2,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2 ! Smoothness indicator IW-1
 ZOMN3(IW-1:IW,:,:) = 3./10. / (ZEPS + ZBNEG3(IW-1:IW,:,:))**2 ! Non-normalized weight IW,IW-1
!
 ! ----- Total flux -----
! 
 PR(IW-1:IW,:,:) = ( ZOMP1(IW-1:IW,:,:)/(ZOMP1(IW-1:IW,:,:)+ZOMP2(IW-1:IW,:,:)+ZOMP3(IW-1:IW,:,:)) * ZFPOS1(IW-1:IW,:,:) &
                       + ZOMP2(IW-1:IW,:,:)/(ZOMP1(IW-1:IW,:,:)+ZOMP2(IW-1:IW,:,:)+ZOMP3(IW-1:IW,:,:)) * ZFPOS2(IW-1:IW,:,:) & 
                       + ZOMP3(IW-1:IW,:,:)/(ZOMP1(IW-1:IW,:,:)+ZOMP2(IW-1:IW,:,:)+ZOMP3(IW-1:IW,:,:)) * ZFPOS3(IW-1:IW,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IW-1:IW,:,:))) &
                    + ( ZOMN1(IW-1:IW,:,:)/(ZOMN1(IW-1:IW,:,:)+ZOMN2(IW-1:IW,:,:)+ZOMN3(IW-1:IW,:,:)) * ZFNEG1(IW-1:IW,:,:) &
                       + ZOMN2(IW-1:IW,:,:)/(ZOMN1(IW-1:IW,:,:)+ZOMN2(IW-1:IW,:,:)+ZOMN3(IW-1:IW,:,:)) * ZFNEG2(IW-1:IW,:,:) &
                       + ZOMN3(IW-1:IW,:,:)/(ZOMN1(IW-1:IW,:,:)+ZOMN2(IW-1:IW,:,:)+ZOMN3(IW-1:IW,:,:)) * ZFNEG3(IW-1:IW,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IW-1:IW,:,:)))
!
 END IF ! NHALO
!
END IF ! IF(LWEST_ll()) 
!
!-------------------------------------------------------------------------------
!*       1.3.   East border
!               ---------------------
!
!! IF( LEAST_ll() .AND. .FALSE. ) THEN 
IF( LEAST_ll() ) THEN 

 !-----------------------------------------------------------------------------
 ! East border is physical -- IE-1,IE
 !-----------------------------------------------------------------------------
 SELECT CASE (HLBCX(2)) ! X direction LBC type on right side
! 
 CASE ('CYCL')
 !---------------------------------------------------------------------------
 ! Periodic boundary condition
 !---------------------------------------------------------------------------
! 
 IF (LWEST_ll() .AND. .FALSE. ) THEN  ! West boundary is physical (monoproc)
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IE-1,:,:) = 1./6 * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:) - PSRC(IE+1,:,:)) ! Flux IE-1
 ZFPOS3(IE,:,:)   = 1./6 * (2.0*PSRC(IE,:,:) + 5.0*PSRC(IE+1,:,:) - PSRC(IW,:,:)) ! Flux IE
 ZBPOS3(IE-1,:,:) = 13./12 * (    PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 ! Smoothness indicator IE-1
 ZBPOS3(IE,:,:)   = 13./12 * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2  ! Smoothness indicator IE
 ZOMP3(IE-1:IE,:,:)  = 3./10. / (ZEPS + ZBPOS3(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE 
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(IE-1,:,:) = 1./6. * (11.0*PSRC(IE,:,:)   - 7.0*PSRC(IE+1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IE-1
 ZFNEG1(IE,:,:)   = 1./6. * (11.0*PSRC(IE+1,:,:) - 7.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IE
 ZBNEG1(IE-1,:,:) = 13./12. * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2  ! Smoothness indicator IE-1
 ZBNEG1(IE,:,:)   = 13./12. * (    PSRC(IE+1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE+1,:,:) - 4.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2  ! Smoothness indicator IE
 ZOMN1(IE-1:IE,:,:) = 1./10. / (ZEPS + ZBNEG1(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(IE-1,:,:) = 1./6. * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   - 1.0*PSRC(IE+1,:,:)) ! Flux IE-1
 ZFNEG2(IE,:,:)   = 1./6. * (2.0*PSRC(IE,:,:)   + 5.0*PSRC(IE+1,:,:) - 1.0*PSRC(IW,:,:)) ! Flux IE
 ZBNEG2(IE-1,:,:)  = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                    PSRC(IE+1,:,:))**2 ! Smoothness indicator IE-1
 ZBNEG2(IE,:,:)   = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (PSRC(IE,:,:) -                      PSRC(IW,:,:))**2  ! Smoothness indicator IE
 ZOMN2(IE-1:IE,:,:) = 3./5. / (ZEPS + ZBNEG2(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ELSEIF(IE<=SIZE(PSRC,1)-3) THEN ! West boundary is proc border, with minimum 3 HALO points on east side
!
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IE-1,:,:) = 1./6 * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:) - PSRC(IE+1,:,:)) ! Flux IE-1
 ZFPOS3(IE,:,:)   = 1./6 * (2.0*PSRC(IE,:,:) + 5.0*PSRC(IE+1,:,:) - PSRC(IE+2,:,:)) ! Flux IE
 ZBPOS3(IE-1,:,:) = 13./12 * (    PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 ! Smoothness indicator IE-1
 ZBPOS3(IE,:,:)   = 13./12 * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE
 ZOMP3(IE-1:IE,:,:)  = 3./10. / (ZEPS + ZBPOS3(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(IE-1,:,:) = 1./6. * (11.0*PSRC(IE,:,:)   - 7.0*PSRC(IE+1,:,:) + 2.0*PSRC(IE+2,:,:)) ! Flux IE-1
 ZFNEG1(IE,:,:)   = 1./6. * (11.0*PSRC(IE+1,:,:) - 7.0*PSRC(IE+2,:,:)   + 2.0*PSRC(IE+3,:,:)) ! Flux IE
 ZBNEG1(IE-1,:,:) = 13./12. * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE-1
 ZBNEG1(IE,:,:)   = 13./12. * (    PSRC(IE+1,:,:) - 2.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE+1,:,:) - 4.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2  ! Smoothness indicator IE
 ZOMN1(IE-1:IE,:,:) = 1./10. / (ZEPS + ZBNEG1(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(IE-1,:,:) = 1./6. * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   - 1.0*PSRC(IE+1,:,:)) ! Flux IE-1
 ZFNEG2(IE,:,:)   = 1./6. * (2.0*PSRC(IE,:,:)   + 5.0*PSRC(IE+1,:,:) - 1.0*PSRC(IE+2,:,:)) ! Flux IE
 ZBNEG2(IE-1,:,:)  = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                    PSRC(IE+1,:,:))**2 ! Smoothness indicator IE-1
 ZBNEG2(IE,:,:)   = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (PSRC(IE,:,:) -                      PSRC(IE+2,:,:))**2  ! Smoothness indicator IE
 ZOMN2(IE-1:IE,:,:) = 3./5. / (ZEPS + ZBNEG2(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ELSE ! West boundary is proc border, with NHALO < 3 on east side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on east side'
  CALL ABORT
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(IE-1,:,:) = 1./6. * (2.0*PSRC(IE-3,:,:) - 7.0*PSRC(IE-2,:,:) + 11.0*PSRC(IE-1,:,:)) ! Flux IE-1
 ZFPOS1(IE,:,:)   = 1./6. * (2.0*PSRC(IE-2,:,:) - 7.0*PSRC(IE-1,:,:) + 11.0*PSRC(IE,:,:)) ! Flux IE
 ZBPOS1(IE-1,:,:) = 13./12. * (PSRC(IE-3,:,:) - 2.0*PSRC(IE-2,:,:) +     PSRC(IE-1,:,:))**2 & 
        + 1./4    * (PSRC(IE-3,:,:) - 4.0*PSRC(IE-2,:,:) + 3.0*PSRC(IE-1,:,:))**2  ! Smoothness indicator IE-1
 ZBPOS1(IE,:,:)   = 13./12. * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 & 
        + 1./4    * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2  ! Smoothness indicator IE
 ZOMP1(IE-1:IE,:,:)  = 1./10. / (ZEPS + ZBPOS1(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(IE-1,:,:) = 1./6. * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE-1
 ZFPOS2(IE,:,:)   = 1./6. * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE
 ZBPOS2(IE-1,:,:) = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) + PSRC(IE,:,:))**2 &
      + 1./4   * (PSRC(IE-2,:,:) -                      PSRC(IE,:,:))**2 ! Smoothness indicator IE-1
 ZBPOS2(IE,:,:)   = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                   PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZOMP2(IE-1:IE,:,:)  = 3./5. / (ZEPS + ZBPOS2(IE-1:IE,:,:))**2  ! Non-normalized weight IE-1,IE
! 
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(IE-1,:,:) = 1./6 * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE-1
 ZFNEG3(IE,:,:)   = 1./6 * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE
 ZBNEG3(IE-1,:,:) = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 &
        + 1./4   * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2 ! Smoothness indicator IE-1
 ZBNEG3(IE,:,:)   = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) +     PSRC(IE+1,:,:))**2 &
        + 1./4   * (PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + 3.0*PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZOMN3(IE-1:IE,:,:) = 3./10. / (ZEPS + ZBNEG3(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! ----- Total flux -----
! 
 PR(IE-1:IE,:,:) = ( ZOMP1(IE-1:IE,:,:)/(ZOMP1(IE-1:IE,:,:)+ZOMP2(IE-1:IE,:,:)+ZOMP3(IE-1:IE,:,:)) * ZFPOS1(IE-1:IE,:,:) &
                       + ZOMP2(IE-1:IE,:,:)/(ZOMP1(IE-1:IE,:,:)+ZOMP2(IE-1:IE,:,:)+ZOMP3(IE-1:IE,:,:)) * ZFPOS2(IE-1:IE,:,:) & 
                       + ZOMP3(IE-1:IE,:,:)/(ZOMP1(IE-1:IE,:,:)+ZOMP2(IE-1:IE,:,:)+ZOMP3(IE-1:IE,:,:)) * ZFPOS3(IE-1:IE,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IE-1:IE,:,:))) &
                    + ( ZOMN1(IE-1:IE,:,:)/(ZOMN1(IE-1:IE,:,:)+ZOMN2(IE-1:IE,:,:)+ZOMN3(IE-1:IE,:,:)) * ZFNEG1(IE-1:IE,:,:) &
                       + ZOMN2(IE-1:IE,:,:)/(ZOMN1(IE-1:IE,:,:)+ZOMN2(IE-1:IE,:,:)+ZOMN3(IE-1:IE,:,:)) * ZFNEG2(IE-1:IE,:,:) &
                       + ZOMN3(IE-1:IE,:,:)/(ZOMN1(IE-1:IE,:,:)+ZOMN2(IE-1:IE,:,:)+ZOMN3(IE-1:IE,:,:)) * ZFNEG3(IE-1:IE,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IE-1:IE,:,:)))
!
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IE
    PR(IE,:,:) = PSRC(IE,:,:) * (0.5+SIGN(0.5,PRUCT(IE,:,:))) + &
                 PSRC(IE+1,:,:) * &
                 (0.5-SIGN(0.5,PRUCT(IE,:,:)))
!
!   ! WENO scheme order 3, IE-1
    ZFPOS1(IE-1,:,:) = 0.5 * (3.0*PSRC(IE-1,:,:) - PSRC(IE-2,:,:)) ! First positive flux
    ZFPOS2(IE-1,:,:) = 0.5 * ( PSRC(IE-1,:,:) + PSRC(IE,:,:)) ! Second positive flux
    ZBPOS1(IE-1,:,:) = (PSRC(IE-1,:,:) - PSRC(IE-2,:,:))**2 ! First positive smoothness indicator
    ZBPOS2(IE-1,:,:) = (PSRC(IE,:,:)   - PSRC(IE-1,:,:))**2 ! Second positive smoothness indicator
!
    ZFNEG1(IE-1,:,:) = 0.5 * (3.0*PSRC(IE,:,:)   - PSRC(IE+1,:,:)) ! First negative flux
    ZFNEG2(IE-1,:,:) = 0.5 * ( PSRC(IE-1,:,:) + PSRC(IE,:,:)) ! Second negative flux
    ZBNEG1(IE-1,:,:) = (PSRC(IE,:,:)   - PSRC(IE+1,:,:))**2 ! First negative smoothness indicator
    ZBNEG2(IE-1,:,:) = (PSRC(IE-1,:,:) - PSRC(IE,:,:))**2  ! Second negative smoothness indicator
!
    ZOMP1(IE-1,:,:) = 1./3. / (ZEPS + ZBPOS1(IE-1,:,:))**2 ! First positive non-normalized weight
    ZOMP2(IE-1,:,:) = 2./3. / (ZEPS + ZBPOS2(IE-1,:,:))**2 ! Second positive non-normalized weight
    ZOMN1(IE-1,:,:) = 1./3. / (ZEPS + ZBNEG1(IE-1,:,:))**2 ! First negative non-normalized weight
    ZOMN2(IE-1,:,:) = 2./3. / (ZEPS + ZBNEG2(IE-1,:,:))**2 ! Second negative non-normalized weight
! 
    PR(IE-1,:,:) = (ZOMN2(IE-1,:,:)/(ZOMN1(IE-1,:,:)+ZOMN2(IE-1,:,:))*ZFNEG2(IE-1,:,:) + &
      (ZOMN1(IE-1,:,:)/(ZOMN1(IE-1,:,:)+ZOMN2(IE-1,:,:))*ZFNEG1(IE-1,:,:))) &
      *(0.5-SIGN(0.5,PRUCT(IE-1,:,:))) + &
                 (ZOMP2(IE-1,:,:)/(ZOMP1(IE-1,:,:)+ZOMP2(IE-1,:,:))*ZFPOS2(IE-1,:,:) + &
                   (ZOMP1(IE-1,:,:)/(ZOMP1(IE-1,:,:)+ZOMP2(IE-1,:,:))*ZFPOS1(IE-1,:,:))) &
      * (0.5+SIGN(0.5,PRUCT(IE-1,:,:)))  ! Total flux
! 
 END SELECT ! SELECT CASE (HLBCX(2)) ! X direction LBC type on right side
!
ELSE
 !-----------------------------------------------------------------------------
 ! East border is proc border -- IE-1,IE
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/east-int not parallelisable with NHALO < 3' 
 CALL ABORT
 STOP  ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >= 3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
! 
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(IE-1,:,:) = 1./6. * (2.0*PSRC(IE-3,:,:) - 7.0*PSRC(IE-2,:,:) + 11.0*PSRC(IE-1,:,:)) ! Flux IE-1
 ZFPOS1(IE,:,:)   = 1./6. * (2.0*PSRC(IE-2,:,:) - 7.0*PSRC(IE-1,:,:) + 11.0*PSRC(IE,:,:)) ! Flux IE
 ZBPOS1(IE-1,:,:) = 13./12. * (PSRC(IE-3,:,:) - 2.0*PSRC(IE-2,:,:) +     PSRC(IE-1,:,:))**2 & 
        + 1./4    * (PSRC(IE-3,:,:) - 4.0*PSRC(IE-2,:,:) + 3.0*PSRC(IE-1,:,:))**2  ! Smoothness indicator IE-1
 ZBPOS1(IE,:,:)   = 13./12. * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 & 
        + 1./4    * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2  ! Smoothness indicator IE
 ZOMP1(IE-1:IE,:,:)  = 1./10. / (ZEPS + ZBPOS1(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(IE-1,:,:) = 1./6. * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE-1
 ZFPOS2(IE,:,:)   = 1./6. * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE
 ZBPOS2(IE-1,:,:) = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) + PSRC(IE,:,:))**2 &
      + 1./4   * (PSRC(IE-2,:,:) -                      PSRC(IE,:,:))**2 ! Smoothness indicator IE-1
 ZBPOS2(IE,:,:)   = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                   PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZOMP2(IE-1:IE,:,:)  = 3./5. / (ZEPS + ZBPOS2(IE-1:IE,:,:))**2  ! Non-normalized weight IE-1,IE
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IE-1,:,:) = 1./6 * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:) - PSRC(IE+1,:,:)) ! Flux IE-1
 ZFPOS3(IE,:,:)   = 1./6 * (2.0*PSRC(IE,:,:) + 5.0*PSRC(IE+1,:,:) - PSRC(IE+2,:,:)) ! Flux IE
 ZBPOS3(IE-1,:,:) = 13./12 * (    PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 ! Smoothness indicator IE-1
 ZBPOS3(IE,:,:)   = 13./12 * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE
 ZOMP3(IE-1:IE,:,:)  = 3./10. / (ZEPS + ZBPOS3(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(IE-1,:,:) = 1./6. * (11.0*PSRC(IE,:,:)   - 7.0*PSRC(IE+1,:,:) + 2.0*PSRC(IE+2,:,:)) ! Flux IE-1
 ZFNEG1(IE,:,:)   = 1./6. * (11.0*PSRC(IE+1,:,:) - 7.0*PSRC(IE+2,:,:) + 2.0*PSRC(IE+3,:,:)) ! Flux IE
 ZBNEG1(IE-1,:,:) = 13./12. * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE-1
 ZBNEG1(IE,:,:)   = 13./12. * (    PSRC(IE+1,:,:) - 2.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE+1,:,:) - 4.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2  ! Smoothness indicator IE
 ZOMN1(IE-1:IE,:,:) = 1./10. / (ZEPS + ZBNEG1(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(IE-1,:,:) = 1./6. * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   - 1.0*PSRC(IE+1,:,:)) ! Flux IE-1
 ZFNEG2(IE,:,:)   = 1./6. * (2.0*PSRC(IE,:,:)   + 5.0*PSRC(IE+1,:,:) - 1.0*PSRC(IE+2,:,:)) ! Flux IE
 ZBNEG2(IE-1,:,:)  = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                    PSRC(IE+1,:,:))**2 ! Smoothness indicator IE-1
 ZBNEG2(IE,:,:)   = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (PSRC(IE,:,:) -                      PSRC(IE+2,:,:))**2 ! Smoothness indicator IE
 ZOMN2(IE-1:IE,:,:) = 3./5. / (ZEPS + ZBNEG2(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(IE-1,:,:) = 1./6 * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE-1
 ZFNEG3(IE,:,:)   = 1./6 * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE
 ZBNEG3(IE-1,:,:) = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 &
        + 1./4   * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2 ! Smoothness indicator IE-1
 ZBNEG3(IE,:,:)   = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) +     PSRC(IE+1,:,:))**2 &
        + 1./4   * (PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + 3.0*PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZOMN3(IE-1:IE,:,:) = 3./10. / (ZEPS + ZBNEG3(IE-1:IE,:,:))**2 ! Non-normalized weight IE-1,IE
!
 ! ----- Total flux -----
! 
 PR(IE-1:IE,:,:) = ( ZOMP1(IE-1:IE,:,:)/(ZOMP1(IE-1:IE,:,:)+ZOMP2(IE-1:IE,:,:)+ZOMP3(IE-1:IE,:,:)) * ZFPOS1(IE-1:IE,:,:) &
                       + ZOMP2(IE-1:IE,:,:)/(ZOMP1(IE-1:IE,:,:)+ZOMP2(IE-1:IE,:,:)+ZOMP3(IE-1:IE,:,:)) * ZFPOS2(IE-1:IE,:,:) & 
                       + ZOMP3(IE-1:IE,:,:)/(ZOMP1(IE-1:IE,:,:)+ZOMP2(IE-1:IE,:,:)+ZOMP3(IE-1:IE,:,:)) * ZFPOS3(IE-1:IE,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IE-1:IE,:,:))) &
                    + ( ZOMN1(IE-1:IE,:,:)/(ZOMN1(IE-1:IE,:,:)+ZOMN2(IE-1:IE,:,:)+ZOMN3(IE-1:IE,:,:)) * ZFNEG1(IE-1:IE,:,:) &
                       + ZOMN2(IE-1:IE,:,:)/(ZOMN1(IE-1:IE,:,:)+ZOMN2(IE-1:IE,:,:)+ZOMN3(IE-1:IE,:,:)) * ZFNEG2(IE-1:IE,:,:) &
                       + ZOMN3(IE-1:IE,:,:)/(ZOMN1(IE-1:IE,:,:)+ZOMN2(IE-1:IE,:,:)+ZOMN3(IE-1:IE,:,:)) * ZFNEG3(IE-1:IE,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IE-1:IE,:,:)))
! 
 END IF ! NHALO
!
END IF ! IF(LWEST_ll()) 
!-------------------------------------------------------------------------------
!
PR = PR * PRUCT ! Add contravariant flux
!
END SUBROUTINE ADVEC_WENO_K_3_UX
!
!------------------------------------------------------------------------------
!
!     ############################################################
      SUBROUTINE ADVEC_WENO_K_3_MX(HLBCX,PSRC, PRUCT, PR)
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
!! T. Lunet 02/10/2014:  Correction of periodic boudary conditions
!!       Change of structure in order to adapt WENO to NHALOK
!!       Suppression of second layer HALO pointers
!!       Complete code documentation
!!
!------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
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
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFPOS1, ZFPOS2, ZFPOS3
!
! intermediate reconstruction fluxes for negative wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFNEG1, ZFNEG2, ZFNEG3
!
! smoothness indicators for positive wind case
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBPOS1, ZBPOS2, ZBPOS3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBNEG1, ZBNEG2, ZBNEG3
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2, ZOMP3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2, ZOMN3
!
! EPSILON for weno weights calculation
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
ZFPOS3 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZFNEG3 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBPOS3 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZBNEG3 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMP3  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0
ZOMN3  = 0.0 
!
!-------------------------------------------------------------------------------
!*       1.1.   Interior Fluxes 
!              ---------------------
IW=IIB
IE=IIE
!
!-------------------------------------------------------------------------------
! Flux calculation in the physical domain far enough from the boundary 
! WENO scheme order 5, IW+2 -> IE-1
! Computation at the mass point on Ugrid u(i-1/2,j,k)
!-------------------------------------------------------------------------------
!
! ----- Positive fluxes -----
!
! First positive stencil, needs indices i-3, i-2, i-1
ZFPOS1(IW+2:IE-1,:,:) = 1./6.   * (2.0*PSRC(IW-1:IE-4,:,:) - 7.0*PSRC(IW:IE-3,:,:) + 11.0*PSRC(IW+1:IE-2,:,:))   ! Flux
ZBPOS1(IW+2:IE-1,:,:) = 13./12. * (    PSRC(IW-1:IE-4,:,:) - 2.0*PSRC(IW:IE-3,:,:) +      PSRC(IW+1:IE-2,:,:))**2 & 
      + 1./4    * (    PSRC(IW-1:IE-4,:,:) - 4.0*PSRC(IW:IE-3,:,:) + 3.0* PSRC(IW+1:IE-2,:,:))**2  ! Smoothness indicator
ZOMP1(IW+2:IE-1,:,:)  = 1./10. /  (ZEPS + ZBPOS1(IW+2:IE-1,:,:))**2             ! Non-normalized weight
!
! Second positive stencil, needs indices i-2, i-1, i
ZFPOS2(IW+2:IE-1,:,:) = 1./6.  * (-1.0*PSRC(IW:IE-3,:,:) + 5.0*PSRC(IW+1:IE-2,:,:) + 2.0*PSRC(IW+2:IE-1,:,:))  ! Flux
ZBPOS2(IW+2:IE-1,:,:) = 13./12 * (     PSRC(IW:IE-3,:,:) - 2.0*PSRC(IW+1:IE-2,:,:) +     PSRC(IW+2:IE-1,:,:))**2 &
      + 1./4   * (     PSRC(IW:IE-3,:,:) -                               PSRC(IW+2:IE-1,:,:))**2 ! Smoothness indicator
ZOMP2(IW+2:IE-1,:,:)  = 3./5. /  (ZEPS + ZBPOS2(IW+2:IE-1,:,:))**2             ! Non-normalized weight
!
! Third positive stencil, needs indices i-1, i, i+1
ZFPOS3(IW+2:IE-1,:,:) = 1./6   * (2.0*PSRC(IW+1:IE-2,:,:) + 5.0*PSRC(IW+2:IE-1,:,:) - PSRC(IW+3:IE,:,:))  ! Flux
ZBPOS3(IW+2:IE-1,:,:) = 13./12 * ( PSRC(IW+1:IE-2,:,:) - 2.0*PSRC(IW+2:IE-1,:,:) + PSRC(IW+3:IE,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW+1:IE-2,:,:) - 4.0*PSRC(IW+2:IE-1,:,:) + PSRC(IW+3:IE,:,:))**2 ! Smoothness indicator
ZOMP3(IW+2:IE-1,:,:)  = 3./10. / (ZEPS + ZBPOS3(IW+2:IE-1,:,:))**2           ! Non-normalized weight
!
! ----- Negative fluxes ----- 
!
! First negative stencil, needs indices i, i+1, i+2
ZFNEG1(IW+2:IE-1,:,:) = 1./6.   * (11.0*PSRC(IW+2:IE-1,:,:) - 7.0*PSRC(IW+3:IE,:,:) + 2.0*PSRC(IW+4:IE+1,:,:))   ! Flux
ZBNEG1(IW+2:IE-1,:,:) = 13./12. * (     PSRC(IW+2:IE-1,:,:) - 2.0*PSRC(IW+3:IE,:,:) +     PSRC(IW+4:IE+1,:,:))**2 & 
      + 1./4    * (3.0* PSRC(IW+2:IE-1,:,:) - 4.0*PSRC(IW+3:IE,:,:) +     PSRC(IW+4:IE+1,:,:))**2  ! Smoothness indicator
ZOMN1(IW+2:IE-1,:,:)  = 1./10. /  (ZEPS + ZBNEG1(IW+2:IE-1,:,:))**2             ! Non-normalized weight
!
! Second negative stencil, needs indices i-1, i, i+1
ZFNEG2(IW+2:IE-1,:,:) = 1./6.  * (2.0*PSRC(IW+1:IE-2,:,:) + 5.0*PSRC(IW+2:IE-1,:,:) - 1.0*PSRC(IW+3:IE,:,:))  ! Flux
ZBNEG2(IW+2:IE-1,:,:) = 13./12 * (    PSRC(IW+1:IE-2,:,:) - 2.0*PSRC(IW+2:IE-1,:,:) +     PSRC(IW+3:IE,:,:))**2 &
      + 1./4   * (    PSRC(IW+1:IE-2,:,:) -                               PSRC(IW+3:IE,:,:))**2 ! Smoothness indicator
ZOMN2(IW+2:IE-1,:,:)  = 3./5. /  (ZEPS + ZBNEG2(IW+2:IE-1,:,:))**2            ! Non-normalized weight
!
! Third negative stencil, needs indices i-2, i-1, i
ZFNEG3(IW+2:IE-1,:,:) = 1./6   * (-1.0*PSRC(IW:IE-3,:,:) + 5.0*PSRC(IW+1:IE-2,:,:) + 2.0*PSRC(IW+2:IE-1,:,:))  ! Flux
ZBNEG3(IW+2:IE-1,:,:) = 13./12 * (  PSRC(IW:IE-3,:,:) - 2.0*PSRC(IW+1:IE-2,:,:) +     PSRC(IW+2:IE-1,:,:))**2 &
      + 1./4   * (     PSRC(IW:IE-3,:,:) - 4.0*PSRC(IW+1:IE-2,:,:) + 3.0*PSRC(IW+2:IE-1,:,:))**2 ! Smoothness indicator
ZOMN3(IW+2:IE-1,:,:)  = 3./10. / (ZEPS + ZBNEG3(IW+2:IE-1,:,:))**2             ! Non-normalized weight
!
!
! ----- Total flux -----
!
PR(IW+2:IE-1,:,:) = (ZOMP1(IW+2:IE-1,:,:)/(ZOMP1(IW+2:IE-1,:,:)+ZOMP2(IW+2:IE-1,:,:)+ZOMP3(IW+2:IE-1,:,:)) &
           * ZFPOS1(IW+2:IE-1,:,:) + &
                      ZOMP2(IW+2:IE-1,:,:)/(ZOMP1(IW+2:IE-1,:,:)+ZOMP2(IW+2:IE-1,:,:)+ZOMP3(IW+2:IE-1,:,:)) &
           * ZFPOS2(IW+2:IE-1,:,:) + & 
                      ZOMP3(IW+2:IE-1,:,:)/(ZOMP1(IW+2:IE-1,:,:)+ZOMP2(IW+2:IE-1,:,:)+ZOMP3(IW+2:IE-1,:,:)) &
           * ZFPOS3(IW+2:IE-1,:,:)) &
                    * (0.5+SIGN(0.5,PRUCT(IW+2:IE-1,:,:))) &
                  + (ZOMN1(IW+2:IE-1,:,:)/(ZOMN1(IW+2:IE-1,:,:)+ZOMN2(IW+2:IE-1,:,:)+ZOMN3(IW+2:IE-1,:,:)) &
           * ZFNEG1(IW+2:IE-1,:,:)  &
                     + ZOMN2(IW+2:IE-1,:,:)/(ZOMN1(IW+2:IE-1,:,:)+ZOMN2(IW+2:IE-1,:,:)+ZOMN3(IW+2:IE-1,:,:)) &
           * ZFNEG2(IW+2:IE-1,:,:)  &
                     + ZOMN3(IW+2:IE-1,:,:)/(ZOMN1(IW+2:IE-1,:,:)+ZOMN2(IW+2:IE-1,:,:)+ZOMN3(IW+2:IE-1,:,:)) &
           * ZFNEG3(IW+2:IE-1,:,:))  & 
                    * (0.5-SIGN(0.5,PRUCT(IW+2:IE-1,:,:)))
!
!-------------------------------------------------------------------------------
!*       1.2.   West border
!               ---------------------
!
!! IF( LWEST_ll()  .AND. .FALSE. ) THEN 
IF( LWEST_ll() ) THEN 
 !-----------------------------------------------------------------------------
 ! West border is physical -- IW+1,IW
 !-----------------------------------------------------------------------------
 SELECT CASE (HLBCX(1)) ! X direction LBC type on left side
! 
 CASE ('CYCL')
 !---------------------------------------------------------------------------
 ! Periodic boundary condition
 !---------------------------------------------------------------------------
!
 IF(LEAST_ll()  .AND. .FALSE. ) THEN  ! East border is physical
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(IW+1,:,:) = 1./6. * (2.0*PSRC(IE,:,:)   - 7.0*PSRC(IW-1,:,:) + 11.0*PSRC(IW,:,:)) ! Flux IW+1
 ZFPOS1(IW,:,:)   = 1./6. * (2.0*PSRC(IE-1,:,:) - 7.0*PSRC(IE,:,:)   + 11.0*PSRC(IW-1,:,:)) ! Flux IW
 ZBPOS1(IW+1,:,:) = 13./12. * (PSRC(IE,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 & 
        + 1./4    * (PSRC(IE,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2   ! Smoothness indicator IW+1
 ZBPOS1(IW,:,:) = 13./12. * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) +     PSRC(IW-1,:,:))**2 & 
      + 1./4    * (PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + 3.0*PSRC(IW-1,:,:))**2  ! Smoothness indicator IW
 ZOMP1(IW:IW+1,:,:)  = 1./10. / (ZEPS + ZBPOS1(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(IW+1,:,:) = 1./6.* (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW+1
 ZFPOS2(IW,:,:)   = 1./6.* (-1.0*PSRC(IE,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW
 ZBPOS2(IW+1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW+1
 ZBPOS2(IW,:,:)   = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IW-1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (PSRC(IE,:,:) -                      PSRC(IW,:,:))**2  ! Smoothness indicator IW
 ZOMP2(IW:IW+1,:,:)  = 3./5. / (ZEPS + ZBPOS2(IW:IW+1,:,:))**2  ! Non-normalized weight IW+1,IW
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(IW+1,:,:) = 1./6 * (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW+1
 ZFNEG3(IW,:,:)   = 1./6 * (-1.0*PSRC(IE,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW
 ZBNEG3(IW+1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) +     PSRC(IW+1,:,:))**2 &
        + 1./4   * (PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + 3.0*PSRC(IW+1,:,:))**2 ! Smoothness indicator IW+1
 ZBNEG3(IW,:,:)   = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 &
        + 1./4   * (PSRC(IE,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2 ! Smoothness indicator IW
 ZOMN3(IW:IW+1,:,:) = 3./10. / (ZEPS + ZBNEG3(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ELSEIF(IW>3) THEN ! East boundary is proc border, with minimum 3 HALO points on west side
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(IW+1,:,:) = 1./6. * (2.0*PSRC(IW-2,:,:)   - 7.0*PSRC(IW-1,:,:) + 11.0*PSRC(IW,:,:)) ! Flux IW+1
 ZFPOS1(IW,:,:)   = 1./6. * (2.0*PSRC(IW-3,:,:) - 7.0*PSRC(IW-2,:,:)   + 11.0*PSRC(IW-1,:,:)) ! Flux IW
 ZBPOS1(IW+1,:,:) = 13./12. * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 & 
        + 1./4    * (PSRC(IW-2,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2   ! Smoothness indicator IW+1
 ZBPOS1(IW,:,:) = 13./12. * (PSRC(IW-3,:,:) - 2.0*PSRC(IW-2,:,:) +     PSRC(IW-1,:,:))**2 & 
      + 1./4    * (PSRC(IW-3,:,:) - 4.0*PSRC(IW-2,:,:) + 3.0*PSRC(IW-1,:,:))**2  ! Smoothness indicator IW
 ZOMP1(IW:IW+1,:,:)  = 1./10. / (ZEPS + ZBPOS1(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(IW+1,:,:) = 1./6.* (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW+1
 ZFPOS2(IW,:,:)   = 1./6.* (-1.0*PSRC(IW-2,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW
 ZBPOS2(IW+1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW+1
 ZBPOS2(IW,:,:)   = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (PSRC(IW-2,:,:) -                      PSRC(IW,:,:))**2  ! Smoothness indicator IW
 ZOMP2(IW:IW+1,:,:)  = 3./5. / (ZEPS + ZBPOS2(IW:IW+1,:,:))**2  ! Non-normalized weight IW+1,IW
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(IW+1,:,:) = 1./6 * (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW+1
 ZFNEG3(IW,:,:)   = 1./6 * (-1.0*PSRC(IW-2,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW
 ZBNEG3(IW+1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) +     PSRC(IW+1,:,:))**2 &
        + 1./4   * (PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + 3.0*PSRC(IW+1,:,:))**2 ! Smoothness indicator IW+1
 ZBNEG3(IW,:,:)   = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 &
        + 1./4   * (PSRC(IW-2,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2 ! Smoothness indicator IW
 ZOMN3(IW:IW+1,:,:) = 3./10. / (ZEPS + ZBNEG3(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ELSE ! East boundary is proc border, with NHALO < 3 on west side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on west side'
  CALL ABORT
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
! 
 ! Third positive stencil, needs indices i-1, i, i+1
 ZFPOS3(IW+1,:,:) = 1./6 * (2.0*PSRC(IW,:,:) + 5.0*PSRC(IW+1,:,:) - PSRC(IW+2,:,:)) ! Flux IW+1
 ZFPOS3(IW,:,:)   = 1./6 * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:) - PSRC(IW+1,:,:)) ! Flux IW
 ZBPOS3(IW+1,:,:) = 13./12 * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 ! Smoothness indicator IW+1
 ZBPOS3(IW,:,:)  = 13./12 * (    PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 ! Smoothness indicator IW
 ZOMP3(IW:IW+1,:,:)  = 3./10. / (ZEPS + ZBPOS3(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(IW+1,:,:) = 1./6. * (11.0*PSRC(IW+1,:,:) - 7.0*PSRC(IW+2,:,:) + 2.0*PSRC(IW+3,:,:)) ! Flux IW+1
 ZFNEG1(IW,:,:)   = 1./6. * (11.0*PSRC(IW,:,:)   - 7.0*PSRC(IW+1,:,:) + 2.0*PSRC(IW+2,:,:)) ! Flux IW
 ZBNEG1(IW+1,:,:) = 13./12. * (    PSRC(IW+1,:,:) - 2.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW+1,:,:) - 4.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2  ! Smoothness indicator IW+1
 ZBNEG1(IW,:,:)   = 13./12. * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2  ! Smoothness indicator IW
 ZOMN1(IW:IW+1,:,:) = 1./10. / (ZEPS + ZBNEG1(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(IW+1,:,:) = 1./6. * (2.0*PSRC(IW,:,:)   + 5.0*PSRC(IW+1,:,:) - 1.0*PSRC(IW+2,:,:)) ! Flux IW+1
 ZFNEG2(IW,:,:)   = 1./6. * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   - 1.0*PSRC(IW+1,:,:)) ! Flux IW
 ZBNEG2(IW+1,:,:) = 13./12 * (PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (PSRC(IW,:,:) -                      PSRC(IW+2,:,:))**2 ! Smoothness indicator IW+1
 ZBNEG2(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW
 ZOMN2(IW:IW+1,:,:) = 3./5. / (ZEPS + ZBNEG2(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! ----- Total flux -----
! 
 PR(IW:IW+1,:,:) = ( ZOMP1(IW:IW+1,:,:)/(ZOMP1(IW:IW+1,:,:)+ZOMP2(IW:IW+1,:,:)+ZOMP3(IW:IW+1,:,:)) * ZFPOS1(IW:IW+1,:,:) &
                       + ZOMP2(IW:IW+1,:,:)/(ZOMP1(IW:IW+1,:,:)+ZOMP2(IW:IW+1,:,:)+ZOMP3(IW:IW+1,:,:)) * ZFPOS2(IW:IW+1,:,:) & 
                       + ZOMP3(IW:IW+1,:,:)/(ZOMP1(IW:IW+1,:,:)+ZOMP2(IW:IW+1,:,:)+ZOMP3(IW:IW+1,:,:)) * ZFPOS3(IW:IW+1,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IW:IW+1,:,:))) &
                    + ( ZOMN1(IW:IW+1,:,:)/(ZOMN1(IW:IW+1,:,:)+ZOMN2(IW:IW+1,:,:)+ZOMN3(IW:IW+1,:,:)) * ZFNEG1(IW:IW+1,:,:) &
                       + ZOMN2(IW:IW+1,:,:)/(ZOMN1(IW:IW+1,:,:)+ZOMN2(IW:IW+1,:,:)+ZOMN3(IW:IW+1,:,:)) * ZFNEG2(IW:IW+1,:,:) &
                       + ZOMN3(IW:IW+1,:,:)/(ZOMN1(IW:IW+1,:,:)+ZOMN2(IW:IW+1,:,:)+ZOMN3(IW:IW+1,:,:)) * ZFNEG3(IW:IW+1,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IW:IW+1,:,:)))
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IW
    PR(IW,:,:) = PSRC(IW-1,:,:) * (0.5+SIGN(0.5,PRUCT(IW,:,:))) + &
                 PSRC(IW,:,:)  * (0.5-SIGN(0.5,PRUCT(IW,:,:)))
!
!   ! WENO scheme order 3, IW+1
    ZFPOS1(IW+1,:,:) = 0.5 * (3.0*PSRC(IW,:,:) - PSRC(IW-1,:,:)) ! First positive flux
    ZFPOS2(IW+1,:,:) = 0.5 * ( PSRC(IW,:,:) + PSRC(IW+1,:,:)) ! Second positive flux
    ZBPOS1(IW+1,:,:) = (PSRC(IW,:,:)   - PSRC(IW-1,:,:))**2 ! First positive smoothness indicator
    ZBPOS2(IW+1,:,:) = (PSRC(IW+1,:,:) - PSRC(IW,:,:))**2  ! Second positive smoothness indicator
!
    ZFNEG1(IW+1,:,:) = 0.5 * (3.0*PSRC(IW+1,:,:) - PSRC(IW+2,:,:)) ! First negative flux
    ZFNEG2(IW+1,:,:) = 0.5 * ( PSRC(IW,:,:)   + PSRC(IW+1,:,:)) ! Second negative flux
    ZBNEG1(IW+1,:,:) = (PSRC(IW+1,:,:) - PSRC(IW+2,:,:))**2 ! First negative smoothness indicator
    ZBNEG2(IW+1,:,:) = (PSRC(IW,:,:)   - PSRC(IW+1,:,:))**2 ! Second negative smoothness indicator
!
    ZOMP1(IW+1,:,:) = 1./3. / (ZEPS + ZBPOS1(IW+1,:,:))**2 ! First positive non-normalized weight
    ZOMP2(IW+1,:,:) = 2./3. / (ZEPS + ZBPOS2(IW+1,:,:))**2 ! Second positive non-normalized weight
    ZOMN1(IW+1,:,:) = 1./3. / (ZEPS + ZBNEG1(IW+1,:,:))**2 ! First negative non-normalized weight
    ZOMN2(IW+1,:,:) = 2./3. / (ZEPS + ZBNEG2(IW+1,:,:))**2 ! Second negative non-normalized weight
! 
    PR(IW+1,:,:) = (ZOMN2(IW+1,:,:)/(ZOMN1(IW+1,:,:)+ZOMN2(IW+1,:,:)) * ZFNEG2(IW+1,:,:) + &
            (ZOMN1(IW+1,:,:)/(ZOMN1(IW+1,:,:)+ZOMN2(IW+1,:,:)) * ZFNEG1(IW+1,:,:))) &
          *(0.5-SIGN(0.5,PRUCT(IW+1,:,:))) + &
      (ZOMP2(IW+1,:,:)/(ZOMP1(IW+1,:,:)+ZOMP2(IW+1,:,:)) * ZFPOS2(IW+1,:,:) + &
      (ZOMP1(IW+1,:,:)/(ZOMP1(IW+1,:,:)+ZOMP2(IW+1,:,:)) * ZFPOS1(IW+1,:,:))) &
      *(0.5+SIGN(0.5,PRUCT(IW+1,:,:)))  ! Total flux
!
 END SELECT ! SELECT CASE (HLBCX(1)) ! X direction LBC type on left side
!
ELSE
 !-----------------------------------------------------------------------------
 ! West border is proc border -- IW,IW-1
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/west-int not parallelisable with NHALO < 3' 
 CALL ABORT
 STOP  ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
!
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(IW+1,:,:) = 1./6. * (2.0*PSRC(IW-2,:,:)   - 7.0*PSRC(IW-1,:,:) + 11.0*PSRC(IW,:,:)) ! Flux IW+1
 ZFPOS1(IW,:,:)   = 1./6. * (2.0*PSRC(IW-3,:,:) - 7.0*PSRC(IW-2,:,:)   + 11.0*PSRC(IW-1,:,:)) ! Flux IW
 ZBPOS1(IW+1,:,:) = 13./12. * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 & 
        + 1./4    * (PSRC(IW-2,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2   ! Smoothness indicator IW+1
 ZBPOS1(IW,:,:) = 13./12. * (PSRC(IW-3,:,:) - 2.0*PSRC(IW-2,:,:) +     PSRC(IW-1,:,:))**2 & 
      + 1./4    * (PSRC(IW-3,:,:) - 4.0*PSRC(IW-2,:,:) + 3.0*PSRC(IW-1,:,:))**2  ! Smoothness indicator IW
 ZOMP1(IW:IW+1,:,:)  = 1./10. / (ZEPS + ZBPOS1(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(IW+1,:,:) = 1./6.* (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW+1
 ZFPOS2(IW,:,:)   = 1./6.* (-1.0*PSRC(IW-2,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW
 ZBPOS2(IW+1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW+1
 ZBPOS2(IW,:,:)   = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (PSRC(IW-2,:,:) -                      PSRC(IW,:,:))**2  ! Smoothness indicator IW
 ZOMP2(IW:IW+1,:,:)  = 3./5. / (ZEPS + ZBPOS2(IW:IW+1,:,:))**2  ! Non-normalized weight IW+1,IW
! 
 ! Third positive stencil, needs indices i-1, i, i+1
 ZFPOS3(IW+1,:,:) = 1./6 * (2.0*PSRC(IW,:,:) + 5.0*PSRC(IW+1,:,:) - PSRC(IW+2,:,:)) ! Flux IW+1
 ZFPOS3(IW,:,:)   = 1./6 * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:) - PSRC(IW+1,:,:)) ! Flux IW
 ZBPOS3(IW+1,:,:) = 13./12 * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 ! Smoothness indicator IW+1
 ZBPOS3(IW,:,:)  = 13./12 * (    PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 ! Smoothness indicator IW
 ZOMP3(IW:IW+1,:,:)  = 3./10. / (ZEPS + ZBPOS3(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(IW+1,:,:) = 1./6. * (11.0*PSRC(IW+1,:,:) - 7.0*PSRC(IW+2,:,:) + 2.0*PSRC(IW+3,:,:)) ! Flux IW+1
 ZFNEG1(IW,:,:)   = 1./6. * (11.0*PSRC(IW,:,:)   - 7.0*PSRC(IW+1,:,:) + 2.0*PSRC(IW+2,:,:)) ! Flux IW
 ZBNEG1(IW+1,:,:) = 13./12. * (    PSRC(IW+1,:,:) - 2.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW+1,:,:) - 4.0*PSRC(IW+2,:,:) + PSRC(IW+3,:,:))**2  ! Smoothness indicator IW+1
 ZBNEG1(IW,:,:)   = 13./12. * (    PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IW,:,:) - 4.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2  ! Smoothness indicator IW
 ZOMN1(IW:IW+1,:,:) = 1./10. / (ZEPS + ZBNEG1(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(IW+1,:,:) = 1./6. * (2.0*PSRC(IW,:,:)   + 5.0*PSRC(IW+1,:,:) - 1.0*PSRC(IW+2,:,:)) ! Flux IW+1
 ZFNEG2(IW,:,:)   = 1./6. * (2.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   - 1.0*PSRC(IW+1,:,:)) ! Flux IW
 ZBNEG2(IW+1,:,:) = 13./12 * (PSRC(IW,:,:) - 2.0*PSRC(IW+1,:,:) + PSRC(IW+2,:,:))**2 &
      + 1./4   * (PSRC(IW,:,:) -                      PSRC(IW+2,:,:))**2 ! Smoothness indicator IW+1
 ZBNEG2(IW,:,:)   = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 &
      + 1./4   * (PSRC(IW-1,:,:) -                    PSRC(IW+1,:,:))**2 ! Smoothness indicator IW
 ZOMN2(IW:IW+1,:,:) = 3./5. / (ZEPS + ZBNEG2(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(IW+1,:,:) = 1./6 * (-1.0*PSRC(IW-1,:,:) + 5.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IW+1
 ZFNEG3(IW,:,:)   = 1./6 * (-1.0*PSRC(IW-2,:,:)   + 5.0*PSRC(IW-1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IW
 ZBNEG3(IW+1,:,:) = 13./12 * (PSRC(IW-1,:,:) - 2.0*PSRC(IW,:,:) +     PSRC(IW+1,:,:))**2 &
        + 1./4   * (PSRC(IW-1,:,:) - 4.0*PSRC(IW,:,:) + 3.0*PSRC(IW+1,:,:))**2 ! Smoothness indicator IW+1
 ZBNEG3(IW,:,:)   = 13./12 * (PSRC(IW-2,:,:) - 2.0*PSRC(IW-1,:,:) +     PSRC(IW,:,:))**2 &
        + 1./4   * (PSRC(IW-2,:,:) - 4.0*PSRC(IW-1,:,:) + 3.0*PSRC(IW,:,:))**2 ! Smoothness indicator IW
 ZOMN3(IW:IW+1,:,:) = 3./10. / (ZEPS + ZBNEG3(IW:IW+1,:,:))**2 ! Non-normalized weight IW+1,IW
!
 ! ----- Total flux -----
! 
 PR(IW:IW+1,:,:) = ( ZOMP1(IW:IW+1,:,:)/(ZOMP1(IW:IW+1,:,:)+ZOMP2(IW:IW+1,:,:)+ZOMP3(IW:IW+1,:,:)) * ZFPOS1(IW:IW+1,:,:) &
                       + ZOMP2(IW:IW+1,:,:)/(ZOMP1(IW:IW+1,:,:)+ZOMP2(IW:IW+1,:,:)+ZOMP3(IW:IW+1,:,:)) * ZFPOS2(IW:IW+1,:,:) & 
                       + ZOMP3(IW:IW+1,:,:)/(ZOMP1(IW:IW+1,:,:)+ZOMP2(IW:IW+1,:,:)+ZOMP3(IW:IW+1,:,:)) * ZFPOS3(IW:IW+1,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IW:IW+1,:,:))) &
                    + ( ZOMN1(IW:IW+1,:,:)/(ZOMN1(IW:IW+1,:,:)+ZOMN2(IW:IW+1,:,:)+ZOMN3(IW:IW+1,:,:)) * ZFNEG1(IW:IW+1,:,:) &
                       + ZOMN2(IW:IW+1,:,:)/(ZOMN1(IW:IW+1,:,:)+ZOMN2(IW:IW+1,:,:)+ZOMN3(IW:IW+1,:,:)) * ZFNEG2(IW:IW+1,:,:) &
                       + ZOMN3(IW:IW+1,:,:)/(ZOMN1(IW:IW+1,:,:)+ZOMN2(IW:IW+1,:,:)+ZOMN3(IW:IW+1,:,:)) * ZFNEG3(IW:IW+1,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IW:IW+1,:,:)))
!
 END IF ! NHALO
!
END IF ! IF(LWEST_ll()) 
!
!-------------------------------------------------------------------------------
!*       1.3.   East border
!               ---------------------
!
!! IF(LEAST_ll()  .AND. .FALSE. ) THEN 
IF(LEAST_ll() ) THEN 
 !-----------------------------------------------------------------------------
 ! East border is physical -- IE,IE+1
 !-----------------------------------------------------------------------------
 SELECT CASE (HLBCX(2)) ! X direction LBC type on right side
! 
 CASE ('CYCL')
 !---------------------------------------------------------------------------
 ! Periodic boundary condition
 !---------------------------------------------------------------------------
!
 IF(LWEST_ll()  .AND. .FALSE. ) THEN  ! West border is physical 
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IE,:,:)   = 1./6 * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:) - PSRC(IE+1,:,:)) ! Flux IE
 ZFPOS3(IE+1,:,:) = 1./6 * (2.0*PSRC(IE,:,:) + 5.0*PSRC(IE+1,:,:) - PSRC(IW,:,:)) ! Flux IE+1
 ZBPOS3(IE,:,:)   = 13./12 * (    PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZBPOS3(IE+1,:,:) = 13./12 * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2  ! Smoothness indicator IE+1
 ZOMP3(IE:IE+1,:,:) = 3./10. / (ZEPS + ZBPOS3(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(IE,:,:)   = 1./6. * (11.0*PSRC(IE,:,:)   - 7.0*PSRC(IE+1,:,:) + 2.0*PSRC(IW,:,:)) ! Flux IE
 ZFNEG1(IE+1,:,:) = 1./6. * (11.0*PSRC(IE+1,:,:) - 7.0*PSRC(IW,:,:)   + 2.0*PSRC(IW+1,:,:)) ! Flux IE+1
 ZBNEG1(IE,:,:)   = 13./12. * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2  ! Smoothness indicator IE
 ZBNEG1(IE+1,:,:) = 13./12. * (    PSRC(IE+1,:,:) - 2.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE+1,:,:) - 4.0*PSRC(IW,:,:) + PSRC(IW+1,:,:))**2  ! Smoothness indicator IE+1
 ZOMN1(IE:IE+1,:,:) = 1./10. / (ZEPS + ZBNEG1(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(IE,:,:)   = 1./6. * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   - 1.0*PSRC(IE+1,:,:)) ! Flux IE
 ZFNEG2(IE+1,:,:) = 1./6. * (2.0*PSRC(IE,:,:)   + 5.0*PSRC(IE+1,:,:) - 1.0*PSRC(IW,:,:)) ! Flux IE+1
 ZBNEG2(IE,:,:)   = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                    PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZBNEG2(IE+1,:,:) = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IW,:,:))**2 &
      + 1./4   * (PSRC(IE,:,:) -                      PSRC(IW,:,:))**2  ! Smoothness indicator IE+1
 ZOMN2(IE:IE+1,:,:) = 3./5. / (ZEPS + ZBNEG2(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ELSEIF(IE<=SIZE(PSRC,1)-3) THEN ! West boundary is proc border, with minimum 3 HALO points on east side
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IE,:,:)   = 1./6 * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:) - PSRC(IE+1,:,:)) ! Flux IE
 ZFPOS3(IE+1,:,:) = 1./6 * (2.0*PSRC(IE,:,:) + 5.0*PSRC(IE+1,:,:) - PSRC(IE+2,:,:)) ! Flux IE+1
 ZBPOS3(IE,:,:)   = 13./12 * (    PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZBPOS3(IE+1,:,:) = 13./12 * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE+1
 ZOMP3(IE:IE+1,:,:) = 3./10. / (ZEPS + ZBPOS3(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(IE,:,:)   = 1./6. * (11.0*PSRC(IE,:,:)   - 7.0*PSRC(IE+1,:,:) + 2.0*PSRC(IE+2,:,:)) ! Flux IE
 ZFNEG1(IE+1,:,:) = 1./6. * (11.0*PSRC(IE+1,:,:) - 7.0*PSRC(IE+2,:,:)   + 2.0*PSRC(IE+3,:,:)) ! Flux IE+1
 ZBNEG1(IE,:,:)   = 13./12. * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE
 ZBNEG1(IE+1,:,:) = 13./12. * (    PSRC(IE+1,:,:) - 2.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE+1,:,:) - 4.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2  ! Smoothness indicator IE+1
 ZOMN1(IE:IE+1,:,:) = 1./10. / (ZEPS + ZBNEG1(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(IE,:,:)   = 1./6. * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   - 1.0*PSRC(IE+1,:,:)) ! Flux IE
 ZFNEG2(IE+1,:,:) = 1./6. * (2.0*PSRC(IE,:,:)   + 5.0*PSRC(IE+1,:,:) - 1.0*PSRC(IE+2,:,:)) ! Flux IE+1
 ZBNEG2(IE,:,:)   = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                    PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZBNEG2(IE+1,:,:) = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (PSRC(IE,:,:) -                      PSRC(IE+2,:,:))**2  ! Smoothness indicator IE+1
 ZOMN2(IE:IE+1,:,:) = 3./5. / (ZEPS + ZBNEG2(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ELSE ! West boundary is proc border, with NHALO < 3 on east side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on east side'
  CALL ABORT
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(IE,:,:)   = 1./6. * (2.0*PSRC(IE-3,:,:) - 7.0*PSRC(IE-2,:,:) + 11.0*PSRC(IE-1,:,:)) ! Flux IE
 ZFPOS1(IE+1,:,:) = 1./6. * (2.0*PSRC(IE-2,:,:) - 7.0*PSRC(IE-1,:,:) + 11.0*PSRC(IE,:,:)) ! Flux IE+1
 ZBPOS1(IE,:,:)   = 13./12. * (PSRC(IE-3,:,:) - 2.0*PSRC(IE-2,:,:) +     PSRC(IE-1,:,:))**2 & 
        + 1./4    * (PSRC(IE-3,:,:) - 4.0*PSRC(IE-2,:,:) + 3.0*PSRC(IE-1,:,:))**2  ! Smoothness indicator IE
 ZBPOS1(IE+1,:,:) = 13./12. * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 & 
        + 1./4    * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2  ! Smoothness indicator IE+1
 ZOMP1(IE:IE+1,:,:) = 1./10. / (ZEPS + ZBPOS1(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(IE,:,:)   = 1./6. * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE
 ZFPOS2(IE+1,:,:) = 1./6. * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE+1
 ZBPOS2(IE,:,:)   = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) + PSRC(IE,:,:))**2 &
      + 1./4   * (PSRC(IE-2,:,:) -                      PSRC(IE,:,:))**2 ! Smoothness indicator IE
 ZBPOS2(IE+1,:,:) = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                   PSRC(IE+1,:,:))**2 ! Smoothness indicator IE+1
 ZOMP2(IE:IE+1,:,:)  = 3./5. / (ZEPS + ZBPOS2(IE:IE+1,:,:))**2  ! Non-normalized weight IE,IE+1
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(IE,:,:)   = 1./6 * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE
 ZFNEG3(IE+1,:,:) = 1./6 * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE+1
 ZBNEG3(IE,:,:)   = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 &
        + 1./4   * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2 ! Smoothness indicator IE
 ZBNEG3(IE+1,:,:) = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) +     PSRC(IE+1,:,:))**2 &
        + 1./4   * (PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + 3.0*PSRC(IE+1,:,:))**2 ! Smoothness indicator IE+1
 ZOMN3(IE:IE+1,:,:) = 3./10. / (ZEPS + ZBNEG3(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! ----- Total flux -----
! 
 PR(IE:IE+1,:,:) = ( ZOMP1(IE:IE+1,:,:)/(ZOMP1(IE:IE+1,:,:)+ZOMP2(IE:IE+1,:,:)+ZOMP3(IE:IE+1,:,:)) * ZFPOS1(IE:IE+1,:,:) &
                       + ZOMP2(IE:IE+1,:,:)/(ZOMP1(IE:IE+1,:,:)+ZOMP2(IE:IE+1,:,:)+ZOMP3(IE:IE+1,:,:)) * ZFPOS2(IE:IE+1,:,:) & 
                       + ZOMP3(IE:IE+1,:,:)/(ZOMP1(IE:IE+1,:,:)+ZOMP2(IE:IE+1,:,:)+ZOMP3(IE:IE+1,:,:)) * ZFPOS3(IE:IE+1,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IE:IE+1,:,:))) &
                    + ( ZOMN1(IE:IE+1,:,:)/(ZOMN1(IE:IE+1,:,:)+ZOMN2(IE:IE+1,:,:)+ZOMN3(IE:IE+1,:,:)) * ZFNEG1(IE:IE+1,:,:) &
                       + ZOMN2(IE:IE+1,:,:)/(ZOMN1(IE:IE+1,:,:)+ZOMN2(IE:IE+1,:,:)+ZOMN3(IE:IE+1,:,:)) * ZFNEG2(IE:IE+1,:,:) &
                       + ZOMN3(IE:IE+1,:,:)/(ZOMN1(IE:IE+1,:,:)+ZOMN2(IE:IE+1,:,:)+ZOMN3(IE:IE+1,:,:)) * ZFNEG3(IE:IE+1,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IE:IE+1,:,:)))
!
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IE+1
    PR(IE+1,:,:) = PSRC(IE,:,:)  * (0.5+SIGN(0.5,PRUCT(IE+1,:,:))) + &
                   PSRC(IE+1,:,:) * (0.5-SIGN(0.5,PRUCT(IE+1,:,:)))
!
!   ! WENO scheme order 3, IE
    ZFPOS1(IE,:,:) = 0.5 * (3.0*PSRC(IE-1,:,:) - PSRC(IE-2,:,:)) ! First positive flux
    ZFPOS2(IE,:,:) = 0.5 * ( PSRC(IE-1,:,:) + PSRC(IE,:,:)) ! Second positive flux
    ZBPOS1(IE,:,:) = (PSRC(IE-1,:,:) - PSRC(IE-2,:,:))**2 ! First positive smoothness indicator
    ZBPOS2(IE,:,:) = (PSRC(IE,:,:)   - PSRC(IE-1,:,:))**2 ! Second positive smoothness indicator
!
    ZFNEG1(IE,:,:) = 0.5 * (3.0*PSRC(IE,:,:)   - PSRC(IE+1,:,:)) ! First negative flux
    ZFNEG2(IE,:,:) = 0.5 * ( PSRC(IE-1,:,:) + PSRC(IE,:,:)) ! Second negative flux
    ZBNEG1(IE,:,:) = (PSRC(IE,:,:)   - PSRC(IE+1,:,:))**2 ! First negative smoothness indicator
    ZBNEG2(IE,:,:) = (PSRC(IE-1,:,:) - PSRC(IE,:,:))**2  ! Second negative smoothness indicator
!
    ZOMP1(IE,:,:) = 1./3. / (ZEPS + ZBPOS1(IE,:,:))**2 ! First positive non-normalized weight
    ZOMP2(IE,:,:) = 2./3. / (ZEPS + ZBPOS2(IE,:,:))**2 ! Second positive non-normalized weight
    ZOMN1(IE,:,:) = 1./3. / (ZEPS + ZBNEG1(IE,:,:))**2 ! First negative non-normalized weight
    ZOMN2(IE,:,:) = 2./3. / (ZEPS + ZBNEG2(IE,:,:))**2 ! Second negative non-normalized weight
! 
    PR(IE,:,:) = (ZOMN2(IE,:,:)/(ZOMN1(IE,:,:)+ZOMN2(IE,:,:))*ZFNEG2(IE,:,:) + &
     (ZOMN1(IE,:,:)/(ZOMN1(IE,:,:)+ZOMN2(IE,:,:))*ZFNEG1(IE,:,:))) &
     *(0.5-SIGN(0.5,PRUCT(IE,:,:))) + &
                 (ZOMP2(IE,:,:)/(ZOMP1(IE,:,:)+ZOMP2(IE,:,:))*ZFPOS2(IE,:,:) + &
                 (ZOMP1(IE,:,:)/(ZOMP1(IE,:,:)+ZOMP2(IE,:,:))*ZFPOS1(IE,:,:))) &
     *(0.5+SIGN(0.5,PRUCT(IE,:,:)))  ! Total flux
! 
 END SELECT ! SELECT CASE (HLBCX(2)) ! X direction LBC type on right side
!
ELSE
 !-----------------------------------------------------------------------------
 ! East border is proc border -- IE,IE+1
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/east-int not parallelisable with NHALO < 3' 
 CALL ABORT
 STOP  ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >= 3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
! 
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(IE,:,:)   = 1./6. * (2.0*PSRC(IE-3,:,:) - 7.0*PSRC(IE-2,:,:) + 11.0*PSRC(IE-1,:,:)) ! Flux IE
 ZFPOS1(IE+1,:,:) = 1./6. * (2.0*PSRC(IE-2,:,:) - 7.0*PSRC(IE-1,:,:) + 11.0*PSRC(IE,:,:)) ! Flux IE+1
 ZBPOS1(IE,:,:)   = 13./12. * (PSRC(IE-3,:,:) - 2.0*PSRC(IE-2,:,:) +     PSRC(IE-1,:,:))**2 & 
        + 1./4    * (PSRC(IE-3,:,:) - 4.0*PSRC(IE-2,:,:) + 3.0*PSRC(IE-1,:,:))**2  ! Smoothness indicator IE
 ZBPOS1(IE+1,:,:) = 13./12. * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 & 
        + 1./4    * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2  ! Smoothness indicator IE+1
 ZOMP1(IE:IE+1,:,:) = 1./10. / (ZEPS + ZBPOS1(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(IE,:,:)   = 1./6. * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE
 ZFPOS2(IE+1,:,:) = 1./6. * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE+1
 ZBPOS2(IE,:,:)   = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) + PSRC(IE,:,:))**2 &
      + 1./4   * (PSRC(IE-2,:,:) -                      PSRC(IE,:,:))**2 ! Smoothness indicator IE
 ZBPOS2(IE+1,:,:) = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                   PSRC(IE+1,:,:))**2 ! Smoothness indicator IE+1
 ZOMP2(IE:IE+1,:,:)  = 3./5. / (ZEPS + ZBPOS2(IE:IE+1,:,:))**2  ! Non-normalized weight IE,IE+1
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(IE,:,:)   = 1./6 * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:) - PSRC(IE+1,:,:)) ! Flux IE
 ZFPOS3(IE+1,:,:) = 1./6 * (2.0*PSRC(IE,:,:) + 5.0*PSRC(IE+1,:,:) - PSRC(IE+2,:,:)) ! Flux IE+1
 ZBPOS3(IE,:,:)   = 13./12 * (    PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZBPOS3(IE+1,:,:) = 13./12 * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE+1
 ZOMP3(IE:IE+1,:,:) = 3./10. / (ZEPS + ZBPOS3(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(IE,:,:)   = 1./6. * (11.0*PSRC(IE,:,:)   - 7.0*PSRC(IE+1,:,:) + 2.0*PSRC(IE+2,:,:)) ! Flux IE
 ZFNEG1(IE+1,:,:) = 1./6. * (11.0*PSRC(IE+1,:,:) - 7.0*PSRC(IE+2,:,:)   + 2.0*PSRC(IE+3,:,:)) ! Flux IE+1
 ZBNEG1(IE,:,:)   = 13./12. * (    PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE,:,:) - 4.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2  ! Smoothness indicator IE
 ZBNEG1(IE+1,:,:) = 13./12. * (    PSRC(IE+1,:,:) - 2.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2 & 
      + 1./4    * (3.0*PSRC(IE+1,:,:) - 4.0*PSRC(IE+2,:,:) + PSRC(IE+3,:,:))**2  ! Smoothness indicator IE+1
 ZOMN1(IE:IE+1,:,:) = 1./10. / (ZEPS + ZBNEG1(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(IE,:,:)   = 1./6. * (2.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   - 1.0*PSRC(IE+1,:,:)) ! Flux IE
 ZFNEG2(IE+1,:,:) = 1./6. * (2.0*PSRC(IE,:,:)   + 5.0*PSRC(IE+1,:,:) - 1.0*PSRC(IE+2,:,:)) ! Flux IE+1
 ZBNEG2(IE,:,:)   = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) + PSRC(IE+1,:,:))**2 &
      + 1./4   * (PSRC(IE-1,:,:) -                    PSRC(IE+1,:,:))**2 ! Smoothness indicator IE
 ZBNEG2(IE+1,:,:) = 13./12 * (PSRC(IE,:,:) - 2.0*PSRC(IE+1,:,:) + PSRC(IE+2,:,:))**2 &
      + 1./4   * (PSRC(IE,:,:) -                      PSRC(IE+2,:,:))**2  ! Smoothness indicator IE+1
 ZOMN2(IE:IE+1,:,:) = 3./5. / (ZEPS + ZBNEG2(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(IE,:,:)   = 1./6 * (-1.0*PSRC(IE-2,:,:) + 5.0*PSRC(IE-1,:,:) + 2.0*PSRC(IE,:,:)) ! Flux IE
 ZFNEG3(IE+1,:,:) = 1./6 * (-1.0*PSRC(IE-1,:,:) + 5.0*PSRC(IE,:,:)   + 2.0*PSRC(IE+1,:,:)) ! Flux IE+1
 ZBNEG3(IE,:,:)   = 13./12 * (PSRC(IE-2,:,:) - 2.0*PSRC(IE-1,:,:) +     PSRC(IE,:,:))**2 &
        + 1./4   * (PSRC(IE-2,:,:) - 4.0*PSRC(IE-1,:,:) + 3.0*PSRC(IE,:,:))**2 ! Smoothness indicator IE
 ZBNEG3(IE+1,:,:) = 13./12 * (PSRC(IE-1,:,:) - 2.0*PSRC(IE,:,:) +     PSRC(IE+1,:,:))**2 &
        + 1./4   * (PSRC(IE-1,:,:) - 4.0*PSRC(IE,:,:) + 3.0*PSRC(IE+1,:,:))**2 ! Smoothness indicator IE+1
 ZOMN3(IE:IE+1,:,:) = 3./10. / (ZEPS + ZBNEG3(IE:IE+1,:,:))**2 ! Non-normalized weight IE,IE+1
!
 ! ----- Total flux -----
! 
 PR(IE:IE+1,:,:) = ( ZOMP1(IE:IE+1,:,:)/(ZOMP1(IE:IE+1,:,:)+ZOMP2(IE:IE+1,:,:)+ZOMP3(IE:IE+1,:,:)) * ZFPOS1(IE:IE+1,:,:) &
                       + ZOMP2(IE:IE+1,:,:)/(ZOMP1(IE:IE+1,:,:)+ZOMP2(IE:IE+1,:,:)+ZOMP3(IE:IE+1,:,:)) * ZFPOS2(IE:IE+1,:,:) & 
                       + ZOMP3(IE:IE+1,:,:)/(ZOMP1(IE:IE+1,:,:)+ZOMP2(IE:IE+1,:,:)+ZOMP3(IE:IE+1,:,:)) * ZFPOS3(IE:IE+1,:,:)) &
                      * (0.5+SIGN(0.5,PRUCT(IE:IE+1,:,:))) &
                    + ( ZOMN1(IE:IE+1,:,:)/(ZOMN1(IE:IE+1,:,:)+ZOMN2(IE:IE+1,:,:)+ZOMN3(IE:IE+1,:,:)) * ZFNEG1(IE:IE+1,:,:) &
                       + ZOMN2(IE:IE+1,:,:)/(ZOMN1(IE:IE+1,:,:)+ZOMN2(IE:IE+1,:,:)+ZOMN3(IE:IE+1,:,:)) * ZFNEG2(IE:IE+1,:,:) &
                       + ZOMN3(IE:IE+1,:,:)/(ZOMN1(IE:IE+1,:,:)+ZOMN2(IE:IE+1,:,:)+ZOMN3(IE:IE+1,:,:)) * ZFNEG3(IE:IE+1,:,:)) &
                      * (0.5-SIGN(0.5,PRUCT(IE:IE+1,:,:)))
!
 END IF ! NHALO
!
END IF ! IF(LWEST_ll()) 
!-------------------------------------------------------------------------------
!
PR = PR * PRUCT ! Add contravariant flux
!
END SUBROUTINE ADVEC_WENO_K_3_MX
!
!-------------------------------------------------------------------------------
!
!     ########################################################################
      SUBROUTINE ADVEC_WENO_K_3_MY(HLBCY,PSRC, PRVCT, PR)
!     ########################################################################
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
!! T. Lunet 02/10/2014:  Correction of periodic boudary conditions
!!       Change of structure in order to adapt WENO to NHALOK
!!       Suppression of second layer HALO pointers
!!       Complete code documentation
!!      J.Escobar : 02/10/2015 : correction on CYCL/OPEN boundaries
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
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
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER::  IS,IN      ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFPOS1, ZFPOS2, ZFPOS3
!
! intermediate reconstruction fluxes for negative wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFNEG1, ZFNEG2, ZFNEG3
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBPOS1, ZBPOS2, ZBPOS3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBNEG1, ZBNEG2, ZBNEG3
!
! WENO weights 
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZOMP1, ZOMP2, ZOMP3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZOMN1, ZOMN2, ZOMN3
!
! EPSILON for weno weights calculation
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
ZFPOS3 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZFNEG3 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBPOS3 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZBNEG3 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMP3  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0
ZOMN3  = 0.0 
!
!-------------------------------------------------------------------------------
!*       1.1.   Interior Fluxes 
!              ---------------------
IS=IJB
IN=IJE
!
!-------------------------------------------------------------------------------
! Flux calculation in the physical domain far enough from the boundary 
! WENO scheme order 5, IS+2 -> IN-1
! Computation at the mass point on Vgrid v(i-1/2,j,k)
!-------------------------------------------------------------------------------
!
! ----- Positive fluxes -----
!
! First positive stencil, needs indices i-3, i-2, i-1
ZFPOS1(:,IS+2:IN-1,:) = 1./6.   * (2.0*PSRC(:,IS-1:IN-4,:)- 7.0*PSRC(:,IS:IN-3,:)+ 11.0*PSRC(:,IS+1:IN-2,:))   ! Flux
ZBPOS1(:,IS+2:IN-1,:) = 13./12. * (    PSRC(:,IS-1:IN-4,:)- 2.0*PSRC(:,IS:IN-3,:)+      PSRC(:,IS+1:IN-2,:))**2 &
      + 1./4    * (    PSRC(:,IS-1:IN-4,:)- 4.0*PSRC(:,IS:IN-3,:)+ 3.0* PSRC(:,IS+1:IN-2,:))**2  ! Smoothness indicator
ZOMP1(:,IS+2:IN-1,:)  = 1./10. /  (ZEPS + ZBPOS1(:,IS+2:IN-1,:))**2             ! Non-normalized weight
!
! Second positive stencil, needs indices i-2, i-1, i
ZFPOS2(:,IS+2:IN-1,:) = 1./6.  * (-1.0*PSRC(:,IS:IN-3,:) + 5.0*PSRC(:,IS+1:IN-2,:) + 2.0*PSRC(:,IS+2:IN-1,:))  ! Flux
ZBPOS2(:,IS+2:IN-1,:) = 13./12 * (     PSRC(:,IS:IN-3,:) - 2.0*PSRC(:,IS+1:IN-2,:) +     PSRC(:,IS+2:IN-1,:))**2 &
      + 1./4   * (     PSRC(:,IS:IN-3,:) -                               PSRC(:,IS+2:IN-1,:))**2 ! Smoothness indicator
ZOMP2(:,IS+2:IN-1,:)  = 3./5. /  (ZEPS + ZBPOS2(:,IS+2:IN-1,:))**2             ! Non-normalized weight
!
! Third positive stencil, needs indices i-1, i, i+1
ZFPOS3(:,IS+2:IN-1,:) = 1./6   * (2.0*PSRC(:,IS+1:IN-2,:) + 5.0*PSRC(:,IS+2:IN-1,:) - PSRC(:,IS+3:IN,:))  ! Flux
ZBPOS3(:,IS+2:IN-1,:) = 13./12 * ( PSRC(:,IS+1:IN-2,:) - 2.0*PSRC(:,IS+2:IN-1,:) + PSRC(:,IS+3:IN,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS+1:IN-2,:) - 4.0*PSRC(:,IS+2:IN-1,:) + PSRC(:,IS+3:IN,:))**2 ! Smoothness indicator
ZOMP3(:,IS+2:IN-1,:)  = 3./10. / (ZEPS + ZBPOS3(:,IS+2:IN-1,:))**2           ! Non-normalized weight
!
! ----- Negative fluxes ----- 
!
! First negative stencil, needs indices i, i+1, i+2
ZFNEG1(:,IS+2:IN-1,:) = 1./6.   * (11.0*PSRC(:,IS+2:IN-1,:) - 7.0*PSRC(:,IS+3:IN,:) + 2.0*PSRC(:,IS+4:IN+1,:))   ! Flux
ZBNEG1(:,IS+2:IN-1,:) = 13./12. * (     PSRC(:,IS+2:IN-1,:) - 2.0*PSRC(:,IS+3:IN,:) +     PSRC(:,IS+4:IN+1,:))**2 &
      + 1./4    * (3.0* PSRC(:,IS+2:IN-1,:) - 4.0*PSRC(:,IS+3:IN,:) +     PSRC(:,IS+4:IN+1,:))**2  ! Smoothness indicator
ZOMN1(:,IS+2:IN-1,:)  = 1./10. /  (ZEPS + ZBNEG1(:,IS+2:IN-1,:))**2             ! Non-normalized weight
!
! Second negative stencil, needs indices i-1, i, i+1
ZFNEG2(:,IS+2:IN-1,:) = 1./6.  * (2.0*PSRC(:,IS+1:IN-2,:) + 5.0*PSRC(:,IS+2:IN-1,:) - 1.0*PSRC(:,IS+3:IN,:))  ! Flux
ZBNEG2(:,IS+2:IN-1,:) = 13./12 * (    PSRC(:,IS+1:IN-2,:) - 2.0*PSRC(:,IS+2:IN-1,:) +     PSRC(:,IS+3:IN,:))**2 &
      + 1./4   * (    PSRC(:,IS+1:IN-2,:) -                               PSRC(:,IS+3:IN,:))**2 ! Smoothness indicator
ZOMN2(:,IS+2:IN-1,:)  = 3./5. /  (ZEPS + ZBNEG2(:,IS+2:IN-1,:))**2            ! Non-normalized weight
!
! Third negative stencil, needs indices i-2, i-1, i
ZFNEG3(:,IS+2:IN-1,:) = 1./6   * (-1.0*PSRC(:,IS:IN-3,:) + 5.0*PSRC(:,IS+1:IN-2,:) + 2.0*PSRC(:,IS+2:IN-1,:))  ! Flux
ZBNEG3(:,IS+2:IN-1,:) = 13./12 * (  PSRC(:,IS:IN-3,:) - 2.0*PSRC(:,IS+1:IN-2,:) +     PSRC(:,IS+2:IN-1,:))**2 &
      + 1./4   * (     PSRC(:,IS:IN-3,:) - 4.0*PSRC(:,IS+1:IN-2,:) + 3.0*PSRC(:,IS+2:IN-1,:))**2 ! Smoothness indicator
ZOMN3(:,IS+2:IN-1,:)  = 3./10. / (ZEPS + ZBNEG3(:,IS+2:IN-1,:))**2             ! Non-normalized weight
!
! ----- Total flux -----
!
PR(:,IS+2:IN-1,:) = (ZOMP1(:,IS+2:IN-1,:)/(ZOMP1(:,IS+2:IN-1,:)+ZOMP2(:,IS+2:IN-1,:)+ZOMP3(:,IS+2:IN-1,:)) &
           * ZFPOS1(:,IS+2:IN-1,:) + &
                      ZOMP2(:,IS+2:IN-1,:)/(ZOMP1(:,IS+2:IN-1,:)+ZOMP2(:,IS+2:IN-1,:)+ZOMP3(:,IS+2:IN-1,:)) &
           * ZFPOS2(:,IS+2:IN-1,:) + & 
                      ZOMP3(:,IS+2:IN-1,:)/(ZOMP1(:,IS+2:IN-1,:)+ZOMP2(:,IS+2:IN-1,:)+ZOMP3(:,IS+2:IN-1,:)) &
           * ZFPOS3(:,IS+2:IN-1,:)) &
                    * (0.5+SIGN(0.5,PRVCT(:,IS+2:IN-1,:))) &
                  + (ZOMN1(:,IS+2:IN-1,:)/(ZOMN1(:,IS+2:IN-1,:)+ZOMN2(:,IS+2:IN-1,:)+ZOMN3(:,IS+2:IN-1,:)) &
           * ZFNEG1(:,IS+2:IN-1,:)  &
                     + ZOMN2(:,IS+2:IN-1,:)/(ZOMN1(:,IS+2:IN-1,:)+ZOMN2(:,IS+2:IN-1,:)+ZOMN3(:,IS+2:IN-1,:)) &
           * ZFNEG2(:,IS+2:IN-1,:)  &
                     + ZOMN3(:,IS+2:IN-1,:)/(ZOMN1(:,IS+2:IN-1,:)+ZOMN2(:,IS+2:IN-1,:)+ZOMN3(:,IS+2:IN-1,:)) &
           * ZFNEG3(:,IS+2:IN-1,:))  & 
                    * (0.5-SIGN(0.5,PRVCT(:,IS+2:IN-1,:)))
!
!-------------------------------------------------------------------------------
!*       1.2.   South border
!               ---------------------
!
!! IF(LSOUTH_ll()  .AND. .FALSE. ) THEN 
IF(LSOUTH_ll()) THEN 
 !-----------------------------------------------------------------------------
 ! South border is physical -- IS+1,IS
 !-----------------------------------------------------------------------------
 SELECT CASE (HLBCY(1)) ! Y direction LBC type on south side
! 
 CASE ('CYCL')
 !---------------------------------------------------------------------------
 ! Periodic boundary condition
 !---------------------------------------------------------------------------
!
 IF(LNORTH_ll()  .AND. .FALSE. ) THEN ! North border is physical 
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(:,IS+1,:) = 1./6. * (2.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IS-1,:) + 11.0*PSRC(:,IS,:)) ! Flux IS+1
 ZFPOS1(:,IS,:)   = 1./6. * (2.0*PSRC(:,IN-1,:) - 7.0*PSRC(:,IN,:)   + 11.0*PSRC(:,IS-1,:)) ! Flux IS
 ZBPOS1(:,IS+1,:) = 13./12. * (PSRC(:,IN,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 & 
        + 1./4    * (PSRC(:,IN,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2   ! Smoothness indicator IS+1
 ZBPOS1(:,IS,:) = 13./12. * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) +     PSRC(:,IS-1,:))**2 & 
      + 1./4    * (PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + 3.0*PSRC(:,IS-1,:))**2  ! Smoothness indicator IS
 ZOMP1(:,IS:IS+1,:)  = 1./10. / (ZEPS + ZBPOS1(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(:,IS+1,:) = 1./6.* (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS+1
 ZFPOS2(:,IS,:)   = 1./6.* (-1.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS
 ZBPOS2(:,IS+1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS+1
 ZBPOS2(:,IS,:)   = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IS-1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IS,:))**2  ! Smoothness indicator IS
 ZOMP2(:,IS:IS+1,:)  = 3./5. / (ZEPS + ZBPOS2(:,IS:IS+1,:))**2  ! Non-normalized weight IS+1,IS
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(:,IS+1,:) = 1./6 * (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS+1
 ZFNEG3(:,IS,:)   = 1./6 * (-1.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS
 ZBNEG3(:,IS+1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) +     PSRC(:,IS+1,:))**2 &
        + 1./4   * (PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + 3.0*PSRC(:,IS+1,:))**2 ! Smoothness indicator IS+1
 ZBNEG3(:,IS,:)   = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 &
        + 1./4   * (PSRC(:,IN,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2 ! Smoothness indicator IS
 ZOMN3(:,IS:IS+1,:) = 3./10. / (ZEPS + ZBNEG3(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
! 
 ELSEIF(IS>3) THEN ! North boundary is proc border, with minimum 3 HALO points on sounth side
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(:,IS+1,:) = 1./6. * (2.0*PSRC(:,IS-2,:)   - 7.0*PSRC(:,IS-1,:) + 11.0*PSRC(:,IS,:)) ! Flux IS+1
 ZFPOS1(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS-3,:) - 7.0*PSRC(:,IS-2,:)   + 11.0*PSRC(:,IS-1,:)) ! Flux IS
 ZBPOS1(:,IS+1,:) = 13./12. * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 & 
        + 1./4    * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2   ! Smoothness indicator IS+1
 ZBPOS1(:,IS,:) = 13./12. * (PSRC(:,IS-3,:) - 2.0*PSRC(:,IS-2,:) +     PSRC(:,IS-1,:))**2 & 
      + 1./4    * (PSRC(:,IS-3,:) - 4.0*PSRC(:,IS-2,:) + 3.0*PSRC(:,IS-1,:))**2  ! Smoothness indicator IS
 ZOMP1(:,IS:IS+1,:)  = 1./10. / (ZEPS + ZBPOS1(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(:,IS+1,:) = 1./6.* (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS+1
 ZFPOS2(:,IS,:)   = 1./6.* (-1.0*PSRC(:,IS-2,:) + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS
 ZBPOS2(:,IS+1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS+1
 ZBPOS2(:,IS,:)   = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IS-2,:) -                      PSRC(:,IS,:))**2  ! Smoothness indicator IS
 ZOMP2(:,IS:IS+1,:)  = 3./5. / (ZEPS + ZBPOS2(:,IS:IS+1,:))**2  ! Non-normalized weight IS+1,IS
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(:,IS+1,:) = 1./6 * (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS+1
 ZFNEG3(:,IS,:)   = 1./6 * (-1.0*PSRC(:,IS-2,:) + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS
 ZBNEG3(:,IS+1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) +     PSRC(:,IS+1,:))**2 &
        + 1./4   * (PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + 3.0*PSRC(:,IS+1,:))**2 ! Smoothness indicator IS+1
 ZBNEG3(:,IS,:)   = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 &
        + 1./4   * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2 ! Smoothness indicator IS
 ZOMN3(:,IS:IS+1,:) = 3./10. / (ZEPS + ZBNEG3(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
! 
 ELSE ! North boundary is proc border, with NHALO < 3 on south side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on south side'
  CALL ABORT
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
! 
 ! Third positive stencil, needs indices i-1, i, i+1
 ZFPOS3(:,IS+1,:) = 1./6 * (2.0*PSRC(:,IS,:) + 5.0*PSRC(:,IS+1,:) - PSRC(:,IS+2,:)) ! Flux IS+1
 ZFPOS3(:,IS,:)   = 1./6 * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:) - PSRC(:,IS+1,:)) ! Flux IS
 ZBPOS3(:,IS+1,:) = 13./12 * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 ! Smoothness indicator IS+1
 ZBPOS3(:,IS,:)  = 13./12 * (    PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZOMP3(:,IS:IS+1,:)  = 3./10. / (ZEPS + ZBPOS3(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(:,IS+1,:) = 1./6. * (11.0*PSRC(:,IS+1,:) - 7.0*PSRC(:,IS+2,:) + 2.0*PSRC(:,IS+3,:)) ! Flux IS+1
 ZFNEG1(:,IS,:)   = 1./6. * (11.0*PSRC(:,IS,:)   - 7.0*PSRC(:,IS+1,:) + 2.0*PSRC(:,IS+2,:)) ! Flux IS
 ZBNEG1(:,IS+1,:) = 13./12. * (    PSRC(:,IS+1,:) - 2.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS+1,:) - 4.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2  ! Smoothness indicator IS+1
 ZBNEG1(:,IS,:)   = 13./12. * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2  ! Smoothness indicator IS
 ZOMN1(:,IS:IS+1,:) = 1./10. / (ZEPS + ZBNEG1(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(:,IS+1,:) = 1./6. * (2.0*PSRC(:,IS,:)   + 5.0*PSRC(:,IS+1,:) - 1.0*PSRC(:,IS+2,:)) ! Flux IS+1
 ZFNEG2(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   - 1.0*PSRC(:,IS+1,:)) ! Flux IS
 ZBNEG2(:,IS+1,:) = 13./12 * (PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (PSRC(:,IS,:) -                      PSRC(:,IS+2,:))**2 ! Smoothness indicator IS+1
 ZBNEG2(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZOMN2(:,IS:IS+1,:) = 3./5. / (ZEPS + ZBNEG2(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! ----- Total flux -----
! 
 PR(:,IS:IS+1,:) = ( ZOMP1(:,IS:IS+1,:)/(ZOMP1(:,IS:IS+1,:)+ZOMP2(:,IS:IS+1,:)+ZOMP3(:,IS:IS+1,:)) * ZFPOS1(:,IS:IS+1,:) &
                       + ZOMP2(:,IS:IS+1,:)/(ZOMP1(:,IS:IS+1,:)+ZOMP2(:,IS:IS+1,:)+ZOMP3(:,IS:IS+1,:)) * ZFPOS2(:,IS:IS+1,:) & 
                       + ZOMP3(:,IS:IS+1,:)/(ZOMP1(:,IS:IS+1,:)+ZOMP2(:,IS:IS+1,:)+ZOMP3(:,IS:IS+1,:)) * ZFPOS3(:,IS:IS+1,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IS:IS+1,:))) &
                    + ( ZOMN1(:,IS:IS+1,:)/(ZOMN1(:,IS:IS+1,:)+ZOMN2(:,IS:IS+1,:)+ZOMN3(:,IS:IS+1,:)) * ZFNEG1(:,IS:IS+1,:) &
                       + ZOMN2(:,IS:IS+1,:)/(ZOMN1(:,IS:IS+1,:)+ZOMN2(:,IS:IS+1,:)+ZOMN3(:,IS:IS+1,:)) * ZFNEG2(:,IS:IS+1,:) &
                       + ZOMN3(:,IS:IS+1,:)/(ZOMN1(:,IS:IS+1,:)+ZOMN2(:,IS:IS+1,:)+ZOMN3(:,IS:IS+1,:)) * ZFNEG3(:,IS:IS+1,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IS:IS+1,:)))
!
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IS
    PR(:,IS,:) = PSRC(:,IS-1,:) * (0.5+SIGN(0.5,PRVCT(:,IS,:))) + &
                 PSRC(:,IS,:)  * (0.5-SIGN(0.5,PRVCT(:,IS,:)))
!
!   ! WENO scheme order 3, IS+1
    ZFPOS1(:,IS+1,:) = 0.5 * (3.0*PSRC(:,IS,:) - PSRC(:,IS-1,:)) ! First positive flux
    ZFPOS2(:,IS+1,:) = 0.5 * ( PSRC(:,IS,:) + PSRC(:,IS+1,:)) ! Second positive flux
    ZBPOS1(:,IS+1,:) = (PSRC(:,IS,:)   - PSRC(:,IS-1,:))**2 ! First positive smoothness indicator
    ZBPOS2(:,IS+1,:) = (PSRC(:,IS+1,:) - PSRC(:,IS,:))**2  ! Second positive smoothness indicator
!
    ZFNEG1(:,IS+1,:) = 0.5 * (3.0*PSRC(:,IS+1,:) - PSRC(:,IS+2,:)) ! First negative flux
    ZFNEG2(:,IS+1,:) = 0.5 * ( PSRC(:,IS+1,:) + PSRC(:,IS,:)) ! Second negative flux
    ZBNEG1(:,IS+1,:) = (PSRC(:,IS+1,:) - PSRC(:,IS+2,:))**2 ! First negative smoothness indicator
    ZBNEG2(:,IS+1,:) = (PSRC(:,IS,:)   - PSRC(:,IS+1,:))**2 ! Second negative smoothness indicator
!
    ZOMP1(:,IS+1,:) = 1./3. / (ZEPS + ZBPOS1(:,IS+1,:))**2 ! First positive non-normalized weight
    ZOMP2(:,IS+1,:) = 2./3. / (ZEPS + ZBPOS2(:,IS+1,:))**2 ! Second positive non-normalized weight
    ZOMN1(:,IS+1,:) = 1./3. / (ZEPS + ZBNEG1(:,IS+1,:))**2 ! First negative non-normalized weight
    ZOMN2(:,IS+1,:) = 2./3. / (ZEPS + ZBNEG2(:,IS+1,:))**2 ! Second negative non-normalized weight
! 
    PR(:,IS+1,:) = (ZOMP2(:,IS+1,:)/(ZOMP1(:,IS+1,:)+ZOMP2(:,IS+1,:)) * ZFPOS2(:,IS+1,:) + &
      (ZOMP1(:,IS+1,:)/(ZOMP1(:,IS+1,:)+ZOMP2(:,IS+1,:)) * ZFPOS1(:,IS+1,:))) &
      *(0.5+SIGN(0.5,PRVCT(:,IS+1,:))) + &
      (ZOMN2(:,IS+1,:)/(ZOMN1(:,IS+1,:)+ZOMN2(:,IS+1,:)) * ZFNEG2(:,IS+1,:) + &
            (ZOMN1(:,IS+1,:)/(ZOMN1(:,IS+1,:)+ZOMN2(:,IS+1,:)) * ZFNEG1(:,IS+1,:))) &
          *(0.5-SIGN(0.5,PRVCT(:,IS+1,:)))  ! Total flux
!
 END SELECT ! SELECT CASE (HLBCY(1)) ! Y direction LBC type on south side
!
ELSE
 !-----------------------------------------------------------------------------
 ! South border is proc border -- IS,IS-1
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/south-int not parallelisable with NHALO < 3' 
 CALL ABORT
 STOP  ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
!
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(:,IS+1,:) = 1./6. * (2.0*PSRC(:,IS-2,:)   - 7.0*PSRC(:,IS-1,:) + 11.0*PSRC(:,IS,:)) ! Flux IS+1
 ZFPOS1(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS-3,:) - 7.0*PSRC(:,IS-2,:)   + 11.0*PSRC(:,IS-1,:)) ! Flux IS
 ZBPOS1(:,IS+1,:) = 13./12. * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 & 
        + 1./4    * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2   ! Smoothness indicator IS+1
 ZBPOS1(:,IS,:) = 13./12. * (PSRC(:,IS-3,:) - 2.0*PSRC(:,IS-2,:) +     PSRC(:,IS-1,:))**2 & 
      + 1./4    * (PSRC(:,IS-3,:) - 4.0*PSRC(:,IS-2,:) + 3.0*PSRC(:,IS-1,:))**2  ! Smoothness indicator IS
 ZOMP1(:,IS:IS+1,:)  = 1./10. / (ZEPS + ZBPOS1(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(:,IS+1,:) = 1./6.* (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS+1
 ZFPOS2(:,IS,:)   = 1./6.* (-1.0*PSRC(:,IS-2,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS
 ZBPOS2(:,IS+1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS+1
 ZBPOS2(:,IS,:)   = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IS-2,:) -                      PSRC(:,IS,:))**2  ! Smoothness indicator IS
 ZOMP2(:,IS:IS+1,:)  = 3./5. / (ZEPS + ZBPOS2(:,IS:IS+1,:))**2  ! Non-normalized weight IS+1,IS
! 
 ! Third positive stencil, needs indices i-1, i, i+1
 ZFPOS3(:,IS+1,:) = 1./6 * (2.0*PSRC(:,IS,:) + 5.0*PSRC(:,IS+1,:) - PSRC(:,IS+2,:)) ! Flux IS+1
 ZFPOS3(:,IS,:)   = 1./6 * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:) - PSRC(:,IS+1,:)) ! Flux IS
 ZBPOS3(:,IS+1,:) = 13./12 * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 ! Smoothness indicator IS+1
 ZBPOS3(:,IS,:)  = 13./12 * (    PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZOMP3(:,IS:IS+1,:)  = 3./10. / (ZEPS + ZBPOS3(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(:,IS+1,:) = 1./6. * (11.0*PSRC(:,IS+1,:) - 7.0*PSRC(:,IS+2,:) + 2.0*PSRC(:,IS+3,:)) ! Flux IS+1
 ZFNEG1(:,IS,:)   = 1./6. * (11.0*PSRC(:,IS,:)   - 7.0*PSRC(:,IS+1,:) + 2.0*PSRC(:,IS+2,:)) ! Flux IS
 ZBNEG1(:,IS+1,:) = 13./12. * (    PSRC(:,IS+1,:) - 2.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS+1,:) - 4.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2  ! Smoothness indicator IS+1
 ZBNEG1(:,IS,:)   = 13./12. * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2  ! Smoothness indicator IS
 ZOMN1(:,IS:IS+1,:) = 1./10. / (ZEPS + ZBNEG1(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(:,IS+1,:) = 1./6. * (2.0*PSRC(:,IS,:)   + 5.0*PSRC(:,IS+1,:) - 1.0*PSRC(:,IS+2,:)) ! Flux IS+1
 ZFNEG2(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   - 1.0*PSRC(:,IS+1,:)) ! Flux IS
 ZBNEG2(:,IS+1,:) = 13./12 * (PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (PSRC(:,IS,:) -                      PSRC(:,IS+2,:))**2 ! Smoothness indicator IS+1
 ZBNEG2(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZOMN2(:,IS:IS+1,:) = 3./5. / (ZEPS + ZBNEG2(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(:,IS+1,:) = 1./6 * (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS+1
 ZFNEG3(:,IS,:)   = 1./6 * (-1.0*PSRC(:,IS-2,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS
 ZBNEG3(:,IS+1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) +     PSRC(:,IS+1,:))**2 &
        + 1./4   * (PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + 3.0*PSRC(:,IS+1,:))**2 ! Smoothness indicator IS+1
 ZBNEG3(:,IS,:)   = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 &
        + 1./4   * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2 ! Smoothness indicator IS
 ZOMN3(:,IS:IS+1,:) = 3./10. / (ZEPS + ZBNEG3(:,IS:IS+1,:))**2 ! Non-normalized weight IS+1,IS
!
 ! ----- Total flux -----
! 
 PR(:,IS:IS+1,:) = ( ZOMP1(:,IS:IS+1,:)/(ZOMP1(:,IS:IS+1,:)+ZOMP2(:,IS:IS+1,:)+ZOMP3(:,IS:IS+1,:)) * ZFPOS1(:,IS:IS+1,:) &
                       + ZOMP2(:,IS:IS+1,:)/(ZOMP1(:,IS:IS+1,:)+ZOMP2(:,IS:IS+1,:)+ZOMP3(:,IS:IS+1,:)) * ZFPOS2(:,IS:IS+1,:) & 
                       + ZOMP3(:,IS:IS+1,:)/(ZOMP1(:,IS:IS+1,:)+ZOMP2(:,IS:IS+1,:)+ZOMP3(:,IS:IS+1,:)) * ZFPOS3(:,IS:IS+1,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IS:IS+1,:))) &
                    + ( ZOMN1(:,IS:IS+1,:)/(ZOMN1(:,IS:IS+1,:)+ZOMN2(:,IS:IS+1,:)+ZOMN3(:,IS:IS+1,:)) * ZFNEG1(:,IS:IS+1,:) &
                       + ZOMN2(:,IS:IS+1,:)/(ZOMN1(:,IS:IS+1,:)+ZOMN2(:,IS:IS+1,:)+ZOMN3(:,IS:IS+1,:)) * ZFNEG2(:,IS:IS+1,:) &
                       + ZOMN3(:,IS:IS+1,:)/(ZOMN1(:,IS:IS+1,:)+ZOMN2(:,IS:IS+1,:)+ZOMN3(:,IS:IS+1,:)) * ZFNEG3(:,IS:IS+1,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IS:IS+1,:)))
!
 END IF ! NHALO
!
END IF ! IF(LSOUTH_ll()) 
!
!-------------------------------------------------------------------------------
!*       1.3.   East border
!               ---------------------
!
!! IF(LNORTH_ll()  .AND. .FALSE. ) THEN 
IF( LNORTH_ll() ) THEN 
 !-----------------------------------------------------------------------------
 ! North border is physical -- IN,IN+1
 !-----------------------------------------------------------------------------
 SELECT CASE (HLBCY(2)) ! Y direction LBC type on north side
! 
 CASE ('CYCL')
 !---------------------------------------------------------------------------
 ! Periodic boundary condition
 !---------------------------------------------------------------------------
!
 IF (LSOUTH_ll()  .AND. .FALSE. ) THEN ! South border is physical 
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IN,:)   = 1./6 * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:) - PSRC(:,IN+1,:)) ! Flux IN
 ZFPOS3(:,IN+1,:) = 1./6 * (2.0*PSRC(:,IN,:) + 5.0*PSRC(:,IN+1,:) - PSRC(:,IS,:)) ! Flux IN+1
 ZBPOS3(:,IN,:)   = 13./12 * (    PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZBPOS3(:,IN+1,:) = 13./12 * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2  ! Smoothness indicator IN+1
 ZOMP3(:,IN:IN+1,:) = 3./10. / (ZEPS + ZBPOS3(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(:,IN,:)   = 1./6. * (11.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IN+1,:) + 2.0*PSRC(:,IS,:)) ! Flux IN
 ZFNEG1(:,IN+1,:) = 1./6. * (11.0*PSRC(:,IN+1,:) - 7.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IN+1
 ZBNEG1(:,IN,:)   = 13./12. * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2  ! Smoothness indicator IN
 ZBNEG1(:,IN+1,:) = 13./12. * (    PSRC(:,IN+1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN+1,:) - 4.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2  ! Smoothness indicator IN+1
 ZOMN1(:,IN:IN+1,:) = 1./10. / (ZEPS + ZBNEG1(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   - 1.0*PSRC(:,IN+1,:)) ! Flux IN
 ZFNEG2(:,IN+1,:) = 1./6. * (2.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IN+1,:) - 1.0*PSRC(:,IS,:)) ! Flux IN+1
 ZBNEG2(:,IN,:)   = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                    PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZBNEG2(:,IN+1,:) = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IS,:))**2  ! Smoothness indicator IN+1
 ZOMN2(:,IN:IN+1,:) = 3./5. / (ZEPS + ZBNEG2(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ELSEIF(IN<=SIZE(PSRC,2)-3) THEN ! South boundary is proc border, with minimum 3 HALO points on north side
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IN,:)   = 1./6 * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:) - PSRC(:,IN+1,:)) ! Flux IN
 ZFPOS3(:,IN+1,:) = 1./6 * (2.0*PSRC(:,IN,:) + 5.0*PSRC(:,IN+1,:) - PSRC(:,IN+2,:)) ! Flux IN+1
 ZBPOS3(:,IN,:)   = 13./12 * (    PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZBPOS3(:,IN+1,:) = 13./12 * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2  ! Smoothness indicator IN+1
 ZOMP3(:,IN:IN+1,:) = 3./10. / (ZEPS + ZBPOS3(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(:,IN,:)   = 1./6. * (11.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IN+1,:) + 2.0*PSRC(:,IN+2,:)) ! Flux IN
 ZFNEG1(:,IN+1,:) = 1./6. * (11.0*PSRC(:,IN+1,:) - 7.0*PSRC(:,IN+2,:)   + 2.0*PSRC(:,IN+3,:)) ! Flux IN+1
 ZBNEG1(:,IN,:)   = 13./12. * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2  ! Smoothness indicator IN
 ZBNEG1(:,IN+1,:) = 13./12. * (    PSRC(:,IN+1,:) - 2.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN+1,:) - 4.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2  ! Smoothness indicator IN+1
 ZOMN1(:,IN:IN+1,:) = 1./10. / (ZEPS + ZBNEG1(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   - 1.0*PSRC(:,IN+1,:)) ! Flux IN
 ZFNEG2(:,IN+1,:) = 1./6. * (2.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IN+1,:) - 1.0*PSRC(:,IN+2,:)) ! Flux IN+1
 ZBNEG2(:,IN,:)   = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                    PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZBNEG2(:,IN+1,:) = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IN+2,:))**2  ! Smoothness indicator IN+1
 ZOMN2(:,IN:IN+1,:) = 3./5. / (ZEPS + ZBNEG2(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ELSE ! South boundary is proc border, with NHALO < 3 on south side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on south side'
  CALL ABORT
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN-3,:) - 7.0*PSRC(:,IN-2,:) + 11.0*PSRC(:,IN-1,:)) ! Flux IN
 ZFPOS1(:,IN+1,:) = 1./6. * (2.0*PSRC(:,IN-2,:) - 7.0*PSRC(:,IN-1,:) + 11.0*PSRC(:,IN,:)) ! Flux IN+1
 ZBPOS1(:,IN,:)   = 13./12. * (PSRC(:,IN-3,:) - 2.0*PSRC(:,IN-2,:) +     PSRC(:,IN-1,:))**2 & 
        + 1./4    * (PSRC(:,IN-3,:) - 4.0*PSRC(:,IN-2,:) + 3.0*PSRC(:,IN-1,:))**2  ! Smoothness indicator IN
 ZBPOS1(:,IN+1,:) = 13./12. * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 & 
        + 1./4    * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2  ! Smoothness indicator IN+1
 ZOMP1(:,IN:IN+1,:) = 1./10. / (ZEPS + ZBPOS1(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(:,IN,:)   = 1./6. * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN
 ZFPOS2(:,IN+1,:) = 1./6. * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN+1
 ZBPOS2(:,IN,:)   = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) + PSRC(:,IN,:))**2 &
      + 1./4   * (PSRC(:,IN-2,:) -                      PSRC(:,IN,:))**2 ! Smoothness indicator IN
 ZBPOS2(:,IN+1,:) = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                   PSRC(:,IN+1,:))**2 ! Smoothness indicator IN+1
 ZOMP2(:,IN:IN+1,:)  = 3./5. / (ZEPS + ZBPOS2(:,IN:IN+1,:))**2  ! Non-normalized weight IN,IN+1
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(:,IN,:)   = 1./6 * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN
 ZFNEG3(:,IN+1,:) = 1./6 * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN+1
 ZBNEG3(:,IN,:)   = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 &
        + 1./4   * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2 ! Smoothness indicator IN
 ZBNEG3(:,IN+1,:) = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) +     PSRC(:,IN+1,:))**2 &
        + 1./4   * (PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + 3.0*PSRC(:,IN+1,:))**2 ! Smoothness indicator IN+1
 ZOMN3(:,IN:IN+1,:) = 3./10. / (ZEPS + ZBNEG3(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! ----- Total flux -----
! 
 PR(:,IN:IN+1,:) = ( ZOMP1(:,IN:IN+1,:)/(ZOMP1(:,IN:IN+1,:)+ZOMP2(:,IN:IN+1,:)+ZOMP3(:,IN:IN+1,:)) * ZFPOS1(:,IN:IN+1,:) &
                       + ZOMP2(:,IN:IN+1,:)/(ZOMP1(:,IN:IN+1,:)+ZOMP2(:,IN:IN+1,:)+ZOMP3(:,IN:IN+1,:)) * ZFPOS2(:,IN:IN+1,:) & 
                       + ZOMP3(:,IN:IN+1,:)/(ZOMP1(:,IN:IN+1,:)+ZOMP2(:,IN:IN+1,:)+ZOMP3(:,IN:IN+1,:)) * ZFPOS3(:,IN:IN+1,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IN:IN+1,:))) &
                    + ( ZOMN1(:,IN:IN+1,:)/(ZOMN1(:,IN:IN+1,:)+ZOMN2(:,IN:IN+1,:)+ZOMN3(:,IN:IN+1,:)) * ZFNEG1(:,IN:IN+1,:) &
                       + ZOMN2(:,IN:IN+1,:)/(ZOMN1(:,IN:IN+1,:)+ZOMN2(:,IN:IN+1,:)+ZOMN3(:,IN:IN+1,:)) * ZFNEG2(:,IN:IN+1,:) &
                       + ZOMN3(:,IN:IN+1,:)/(ZOMN1(:,IN:IN+1,:)+ZOMN2(:,IN:IN+1,:)+ZOMN3(:,IN:IN+1,:)) * ZFNEG3(:,IN:IN+1,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IN:IN+1,:)))
!
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IN+1
    PR(:,IN+1,:) = PSRC(:,IN,:)  * (0.5+SIGN(0.5,PRVCT(:,IN+1,:))) + &
                   PSRC(:,IN+1,:) * (0.5-SIGN(0.5,PRVCT(:,IN+1,:)))
!
!   ! WENO scheme order 3, IN
    ZFPOS1(:,IN,:) = 0.5 * (3.0*PSRC(:,IN-1,:) - PSRC(:,IN-2,:)) ! First positive flux
    ZFPOS2(:,IN,:) = 0.5 * ( PSRC(:,IN-1,:) + PSRC(:,IN,:)) ! Second positive flux
    ZBPOS1(:,IN,:) = (PSRC(:,IN-1,:) - PSRC(:,IN-2,:))**2 ! First positive smoothness indicator
    ZBPOS2(:,IN,:) = (PSRC(:,IN,:)   - PSRC(:,IN-1,:))**2 ! Second positive smoothness indicator
!
    ZFNEG1(:,IN,:) = 0.5 * (3.0*PSRC(:,IN,:) - PSRC(:,IN+1,:)) ! First negative flux
    ZFNEG2(:,IN,:) = 0.5 * ( PSRC(:,IN,:) + PSRC(:,IN-1,:)) ! Second negative flux
    ZBNEG1(:,IN,:) = (PSRC(:,IN,:)   - PSRC(:,IN+1,:))**2 ! First negative smoothness indicator
    ZBNEG2(:,IN,:) = (PSRC(:,IN-1,:) - PSRC(:,IN,:))**2  ! Second negative smoothness indicator
!
    ZOMP1(:,IN,:) = 1./3. / (ZEPS + ZBPOS1(:,IN,:))**2 ! First positive non-normalized weight
    ZOMP2(:,IN,:) = 2./3. / (ZEPS + ZBPOS2(:,IN,:))**2 ! Second positive non-normalized weight
    ZOMN1(:,IN,:) = 1./3. / (ZEPS + ZBNEG1(:,IN,:))**2 ! First negative non-normalized weight
    ZOMN2(:,IN,:) = 2./3. / (ZEPS + ZBNEG2(:,IN,:))**2 ! Second negative non-normalized weight
! 
    PR(:,IN,:) = (ZOMP2(:,IN,:)/(ZOMP1(:,IN,:)+ZOMP2(:,IN,:))*ZFPOS2(:,IN,:) + &
                 (ZOMP1(:,IN,:)/(ZOMP1(:,IN,:)+ZOMP2(:,IN,:))*ZFPOS1(:,IN,:))) &
     *(0.5+SIGN(0.5,PRVCT(:,IN,:))) + &
     (ZOMN2(:,IN,:)/(ZOMN1(:,IN,:)+ZOMN2(:,IN,:))*ZFNEG2(:,IN,:) + &
     (ZOMN1(:,IN,:)/(ZOMN1(:,IN,:)+ZOMN2(:,IN,:))*ZFNEG1(:,IN,:))) &
     *(0.5-SIGN(0.5,PRVCT(:,IN,:)))  ! Total flux
!
 END SELECT ! SELECT CASE (HLBCX(2)) ! X direction LBC type on right side
!
ELSE
 !-----------------------------------------------------------------------------
 ! North border is proc border -- IN,IN+1
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/north-int not parallelisable with NHALO < 3' 
 CALL ABORT
 STOP  ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >= 3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
! 
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-3, i-2, i-1 
 ZFPOS1(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN-3,:) - 7.0*PSRC(:,IN-2,:) + 11.0*PSRC(:,IN-1,:)) ! Flux IN
 ZFPOS1(:,IN+1,:) = 1./6. * (2.0*PSRC(:,IN-2,:) - 7.0*PSRC(:,IN-1,:) + 11.0*PSRC(:,IN,:)) ! Flux IN+1
 ZBPOS1(:,IN,:)   = 13./12. * (PSRC(:,IN-3,:) - 2.0*PSRC(:,IN-2,:) +     PSRC(:,IN-1,:))**2 & 
        + 1./4    * (PSRC(:,IN-3,:) - 4.0*PSRC(:,IN-2,:) + 3.0*PSRC(:,IN-1,:))**2  ! Smoothness indicator IN
 ZBPOS1(:,IN+1,:) = 13./12. * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 & 
        + 1./4    * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2  ! Smoothness indicator IN+1
 ZOMP1(:,IN:IN+1,:) = 1./10. / (ZEPS + ZBPOS1(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! Second positive stencil, needs indices i-2, i-1, i
 ZFPOS2(:,IN,:)   = 1./6. * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN
 ZFPOS2(:,IN+1,:) = 1./6. * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN+1
 ZBPOS2(:,IN,:)   = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) + PSRC(:,IN,:))**2 &
      + 1./4   * (PSRC(:,IN-2,:) -                      PSRC(:,IN,:))**2 ! Smoothness indicator IN
 ZBPOS2(:,IN+1,:) = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                   PSRC(:,IN+1,:))**2 ! Smoothness indicator IN+1
 ZOMP2(:,IN:IN+1,:)  = 3./5. / (ZEPS + ZBPOS2(:,IN:IN+1,:))**2  ! Non-normalized weight IN,IN+1
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IN,:)   = 1./6 * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:) - PSRC(:,IN+1,:)) ! Flux IN
 ZFPOS3(:,IN+1,:) = 1./6 * (2.0*PSRC(:,IN,:) + 5.0*PSRC(:,IN+1,:) - PSRC(:,IN+2,:)) ! Flux IN+1
 ZBPOS3(:,IN,:)   = 13./12 * (    PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZBPOS3(:,IN+1,:) = 13./12 * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2  ! Smoothness indicator IN+1
 ZOMP3(:,IN:IN+1,:) = 3./10. / (ZEPS + ZBPOS3(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i, i+1, i+2
 ZFNEG1(:,IN,:)   = 1./6. * (11.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IN+1,:) + 2.0*PSRC(:,IN+2,:)) ! Flux IN
 ZFNEG1(:,IN+1,:) = 1./6. * (11.0*PSRC(:,IN+1,:) - 7.0*PSRC(:,IN+2,:)   + 2.0*PSRC(:,IN+3,:)) ! Flux IN+1
 ZBNEG1(:,IN,:)   = 13./12. * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2  ! Smoothness indicator IN
 ZBNEG1(:,IN+1,:) = 13./12. * (    PSRC(:,IN+1,:) - 2.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN+1,:) - 4.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2  ! Smoothness indicator IN+1
 ZOMN1(:,IN:IN+1,:) = 1./10. / (ZEPS + ZBNEG1(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! Second negative stencil, needs indices i-1, i, i+1
 ZFNEG2(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   - 1.0*PSRC(:,IN+1,:)) ! Flux IN
 ZFNEG2(:,IN+1,:) = 1./6. * (2.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IN+1,:) - 1.0*PSRC(:,IN+2,:)) ! Flux IN+1
 ZBNEG2(:,IN,:)   = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                    PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZBNEG2(:,IN+1,:) = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IN+2,:))**2  ! Smoothness indicator IN+1
 ZOMN2(:,IN:IN+1,:) = 3./5. / (ZEPS + ZBNEG2(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! Third negative stencil, needs indices i-2, i-1, i
 ZFNEG3(:,IN,:)   = 1./6 * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN
 ZFNEG3(:,IN+1,:) = 1./6 * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN+1
 ZBNEG3(:,IN,:)   = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 &
        + 1./4   * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2 ! Smoothness indicator IN
 ZBNEG3(:,IN+1,:) = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) +     PSRC(:,IN+1,:))**2 &
        + 1./4   * (PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + 3.0*PSRC(:,IN+1,:))**2 ! Smoothness indicator IN+1
 ZOMN3(:,IN:IN+1,:) = 3./10. / (ZEPS + ZBNEG3(:,IN:IN+1,:))**2 ! Non-normalized weight IN,IN+1
!
 ! ----- Total flux -----
! 
 PR(:,IN:IN+1,:) = ( ZOMP1(:,IN:IN+1,:)/(ZOMP1(:,IN:IN+1,:)+ZOMP2(:,IN:IN+1,:)+ZOMP3(:,IN:IN+1,:)) * ZFPOS1(:,IN:IN+1,:) &
                       + ZOMP2(:,IN:IN+1,:)/(ZOMP1(:,IN:IN+1,:)+ZOMP2(:,IN:IN+1,:)+ZOMP3(:,IN:IN+1,:)) * ZFPOS2(:,IN:IN+1,:) & 
                       + ZOMP3(:,IN:IN+1,:)/(ZOMP1(:,IN:IN+1,:)+ZOMP2(:,IN:IN+1,:)+ZOMP3(:,IN:IN+1,:)) * ZFPOS3(:,IN:IN+1,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IN:IN+1,:))) &
                    + ( ZOMN1(:,IN:IN+1,:)/(ZOMN1(:,IN:IN+1,:)+ZOMN2(:,IN:IN+1,:)+ZOMN3(:,IN:IN+1,:)) * ZFNEG1(:,IN:IN+1,:) &
                       + ZOMN2(:,IN:IN+1,:)/(ZOMN1(:,IN:IN+1,:)+ZOMN2(:,IN:IN+1,:)+ZOMN3(:,IN:IN+1,:)) * ZFNEG2(:,IN:IN+1,:) &
                       + ZOMN3(:,IN:IN+1,:)/(ZOMN1(:,IN:IN+1,:)+ZOMN2(:,IN:IN+1,:)+ZOMN3(:,IN:IN+1,:)) * ZFNEG3(:,IN:IN+1,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IN:IN+1,:)))
!
 END IF ! NHALO
!
END IF ! IF(LNORTH_ll()) 
!-------------------------------------------------------------------------------
!
PR = PR * PRVCT ! Add contravariant flux
!
END SUBROUTINE ADVEC_WENO_K_3_MY
!
!-------------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE ADVEC_WENO_K_3_VY(HLBCY, PSRC, PRVCT, PR)
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
!! T. Lunet 02/10/2014:  Correction of periodic boudary conditions
!!       Change of structure in order to adapt WENO to NHALOK
!!       Suppression of second layer HALO pointers
!!       Complete code documentation
!!      J.Escobar : 02/10/2015 : correction on CYCL/OPEN boundaries
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_LUNIT
USE MODD_CONF
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
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER::  IS,IN    ! Physical boundary index
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFPOS1, ZFPOS2, ZFPOS3
!
! intermediate reconstruction fluxes for negative wind case
!
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFNEG1, ZFNEG2, ZFNEG3
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBPOS1, ZBPOS2, ZBPOS3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBNEG1, ZBNEG2, ZBNEG3
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMP1, ZOMP2, ZOMP3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: ZOMN1, ZOMN2, ZOMN3
!
! EPSILON for weno weights calculation
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
ZFPOS3 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZFNEG3 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBPOS3 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZBNEG3 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMP3  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0
ZOMN3  = 0.0
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!*       1.1.   Interior Fluxes 
!              ---------------------
IS=IJB
IN=IJE
!
!-------------------------------------------------------------------------------
! Flux calculation in the physical domain far enough from the boundary 
! WENO scheme order 5, IS+1 -> IN-2
! Computation at the mass point on Vgrid v(i+1/2,j,k)
!-------------------------------------------------------------------------------
!
! ----- Positive fluxes -----
!
! First positive stencil, needs indices i-2, i-1, i
ZFPOS1(:,IS+1:IN-2,:) = 1./6.   * (2.0*PSRC(:,IS-1:IN-4,:) - 7.0*PSRC(:,IS:IN-3,:) + 11.0*PSRC(:,IS+1:IN-2,:))   ! Flux
ZBPOS1(:,IS+1:IN-2,:) = 13./12. * (    PSRC(:,IS-1:IN-4,:) - 2.0*PSRC(:,IS:IN-3,:) +      PSRC(:,IS+1:IN-2,:))**2 & 
      + 1./4    * (    PSRC(:,IS-1:IN-4,:) - 4.0*PSRC(:,IS:IN-3,:) + 3.0* PSRC(:,IS+1:IN-2,:))**2  ! Smoothness indicator
ZOMP1(:,IS+1:IN-2,:)  = 1./10. /  (ZEPS + ZBPOS1(:,IS+1:IN-2,:))**2             ! Non-normalized weight
!
! Second positive stencil, needs indices i-1, i, i+1
ZFPOS2(:,IS+1:IN-2,:) = 1./6.  * (-1.0*PSRC(:,IS:IN-3,:) + 5.0*PSRC(:,IS+1:IN-2,:) + 2.0*PSRC(:,IS+2:IN-1,:))  ! Flux
ZBPOS2(:,IS+1:IN-2,:) = 13./12 * (     PSRC(:,IS:IN-3,:) - 2.0*PSRC(:,IS+1:IN-2,:) +     PSRC(:,IS+2:IN-1,:))**2 &
      + 1./4   * (     PSRC(:,IS:IN-3,:) -                               PSRC(:,IS+2:IN-1,:))**2 ! Smoothness indicator
ZOMP2(:,IS+1:IN-2,:)  = 3./5. /  (ZEPS + ZBPOS2(:,IS+1:IN-2,:))**2             ! Non-normalized weight
!
! Third positive stencil, needs indices i, i+1, i+2
ZFPOS3(:,IS+1:IN-2,:) = 1./6   * (2.0*PSRC(:,IS+1:IN-2,:) + 5.0*PSRC(:,IS+2:IN-1,:) - PSRC(:,IS+3:IN,:))  ! Flux
ZBPOS3(:,IS+1:IN-2,:) = 13./12 * ( PSRC(:,IS+1:IN-2,:) - 2.0*PSRC(:,IS+2:IN-1,:) + PSRC(:,IS+3:IN,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS+1:IN-2,:) - 4.0*PSRC(:,IS+2:IN-1,:) + PSRC(:,IS+3:IN,:))**2 ! Smoothness indicator
ZOMP3(:,IS+1:IN-2,:)  = 3./10. / (ZEPS + ZBPOS3(:,IS+1:IN-2,:))**2           ! Non-normalized weight
!
! ----- Negative fluxes ----- 
!
! First negative stencil, needs indices i+1, i+2, i+3
ZFNEG1(:,IS+1:IN-2,:) = 1./6.   * (11.0*PSRC(:,IS+2:IN-1,:) - 7.0*PSRC(:,IS+3:IN,:) + 2.0*PSRC(:,IS+4:IN+1,:))   ! Flux
ZBNEG1(:,IS+1:IN-2,:) = 13./12. * (     PSRC(:,IS+2:IN-1,:) - 2.0*PSRC(:,IS+3:IN,:) +     PSRC(:,IS+4:IN+1,:))**2 & 
      + 1./4    * (3.0* PSRC(:,IS+2:IN-1,:) - 4.0*PSRC(:,IS+3:IN,:) +     PSRC(:,IS+4:IN+1,:))**2  ! Smoothness indicator
ZOMN1(:,IS+1:IN-2,:)  = 1./10. /  (ZEPS + ZBNEG1(:,IS+1:IN-2,:))**2             ! Non-normalized weight
!
! Second negative stencil, needs indices i, i+1, i+2
ZFNEG2(:,IS+1:IN-2,:) = 1./6.  * (2.0*PSRC(:,IS+1:IN-2,:) + 5.0*PSRC(:,IS+2:IN-1,:) - 1.0*PSRC(:,IS+3:IN,:))  ! Flux
ZBNEG2(:,IS+1:IN-2,:) = 13./12 * (    PSRC(:,IS+1:IN-2,:) - 2.0*PSRC(:,IS+2:IN-1,:) +     PSRC(:,IS+3:IN,:))**2 &
      + 1./4   * (    PSRC(:,IS+1:IN-2,:) -                               PSRC(:,IS+3:IN,:))**2 ! Smoothness indicator
ZOMN2(:,IS+1:IN-2,:)  = 3./5. /  (ZEPS + ZBNEG2(:,IS+1:IN-2,:))**2            ! Non-normalized weight
!
! Third negative stencil, needs indices i-1, i, i+1
ZFNEG3(:,IS+1:IN-2,:) = 1./6   * (-1.0*PSRC(:,IS:IN-3,:) + 5.0*PSRC(:,IS+1:IN-2,:) + 2.0*PSRC(:,IS+2:IN-1,:))  ! Flux
ZBNEG3(:,IS+1:IN-2,:) = 13./12 * (  PSRC(:,IS:IN-3,:) - 2.0*PSRC(:,IS+1:IN-2,:) +     PSRC(:,IS+2:IN-1,:))**2 &
      + 1./4   * (     PSRC(:,IS:IN-3,:) - 4.0*PSRC(:,IS+1:IN-2,:) + 3.0*PSRC(:,IS+2:IN-1,:))**2 ! Smoothness indicator
ZOMN3(:,IS+1:IN-2,:)  = 3./10. / (ZEPS + ZBNEG3(:,IS+1:IN-2,:))**2             ! Non-normalized weight
!
!
! ----- Total flux -----
!
PR(:,IS+1:IN-2,:) = (ZOMP1(:,IS+1:IN-2,:)/(ZOMP1(:,IS+1:IN-2,:)+ZOMP2(:,IS+1:IN-2,:)+ZOMP3(:,IS+1:IN-2,:)) &
           * ZFPOS1(:,IS+1:IN-2,:) + &
                      ZOMP2(:,IS+1:IN-2,:)/(ZOMP1(:,IS+1:IN-2,:)+ZOMP2(:,IS+1:IN-2,:)+ZOMP3(:,IS+1:IN-2,:)) &
           * ZFPOS2(:,IS+1:IN-2,:) + & 
                      ZOMP3(:,IS+1:IN-2,:)/(ZOMP1(:,IS+1:IN-2,:)+ZOMP2(:,IS+1:IN-2,:)+ZOMP3(:,IS+1:IN-2,:)) &
           * ZFPOS3(:,IS+1:IN-2,:)) &
                    * (0.5+SIGN(0.5,PRVCT(:,IS+1:IN-2,:))) &
                  + (ZOMN1(:,IS+1:IN-2,:)/(ZOMN1(:,IS+1:IN-2,:)+ZOMN2(:,IS+1:IN-2,:)+ZOMN3(:,IS+1:IN-2,:)) &
           * ZFNEG1(:,IS+1:IN-2,:)  &
                     + ZOMN2(:,IS+1:IN-2,:)/(ZOMN1(:,IS+1:IN-2,:)+ZOMN2(:,IS+1:IN-2,:)+ZOMN3(:,IS+1:IN-2,:)) &
           * ZFNEG2(:,IS+1:IN-2,:)  &
                     + ZOMN3(:,IS+1:IN-2,:)/(ZOMN1(:,IS+1:IN-2,:)+ZOMN2(:,IS+1:IN-2,:)+ZOMN3(:,IS+1:IN-2,:)) &
           * ZFNEG3(:,IS+1:IN-2,:))  & 
                    * (0.5-SIGN(0.5,PRVCT(:,IS+1:IN-2,:)))
!
!-------------------------------------------------------------------------------
!*       1.2.   South border
!               ---------------------
!
!! IF(LSOUTH_ll()  .AND. .FALSE. ) THEN 
IF(LSOUTH_ll() ) THEN 
 !-----------------------------------------------------------------------------
 ! South border is physical -- IS,IS-1
 !-----------------------------------------------------------------------------
 SELECT CASE (HLBCY(1)) ! Y direction LBC type on south side
! 
 CASE ('CYCL')
 !---------------------------------------------------------------------------
 ! Periodic boundary condition
 !---------------------------------------------------------------------------
!
 IF(LNORTH_ll()  .AND. .FALSE. ) THEN ! North border is physical 
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(:,IS,:)   = 1./6. * (2.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IS-1,:) + 11.0*PSRC(:,IS,:)) ! Flux IS
 ZFPOS1(:,IS-1,:) = 1./6. * (2.0*PSRC(:,IN-1,:) - 7.0*PSRC(:,IN,:)   + 11.0*PSRC(:,IS-1,:)) ! Flux IS-1
 ZBPOS1(:,IS,:)   = 13./12. * (PSRC(:,IN,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 & 
        + 1./4    * (PSRC(:,IN,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2   ! Smoothness indicator IS
 ZBPOS1(:,IS-1,:) = 13./12. * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) +     PSRC(:,IS-1,:))**2 & 
        + 1./4    * (PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + 3.0*PSRC(:,IS-1,:))**2  ! Smoothness indicator IS-1
 ZOMP1(:,IS-1:IS,:)  = 1./10. / (ZEPS + ZBPOS1(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(:,IS,:)   = 1./6.* (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS
 ZFPOS2(:,IS-1,:) = 1./6.* (-1.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS-1
 ZBPOS2(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZBPOS2(:,IS-1,:) = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IS-1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IS,:))**2  ! Smoothness indicator IS-1
 ZOMP2(:,IS-1:IS,:)  = 3./5. / (ZEPS + ZBPOS2(:,IS-1:IS,:))**2  ! Non-normalized weight IS,IS-1
!
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(:,IS,:)   = 1./6 * (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS
 ZFNEG3(:,IS-1,:) = 1./6 * (-1.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS-1
 ZBNEG3(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) +     PSRC(:,IS+1,:))**2 &
        + 1./4   * (PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + 3.0*PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZBNEG3(:,IS-1,:) = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 &
        + 1./4   * (PSRC(:,IN,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2 ! Smoothness indicator IS-1
 ZOMN3(:,IS-1:IS,:) = 3./10. / (ZEPS + ZBNEG3(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
! 
 ELSEIF(IS>3) THEN ! North boundary is proc border, with minimum 3 HALO points on south side
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS-2,:)   - 7.0*PSRC(:,IS-1,:) + 11.0*PSRC(:,IS,:)) ! Flux IS
 ZFPOS1(:,IS-1,:) = 1./6. * (2.0*PSRC(:,IS+3,:) - 7.0*PSRC(:,IS-2,:)   + 11.0*PSRC(:,IS-1,:)) ! Flux IS-1
 ZBPOS1(:,IS,:)   = 13./12. * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 & 
        + 1./4    * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2   ! Smoothness indicator IS
 ZBPOS1(:,IS-1,:) = 13./12. * (PSRC(:,IS-3,:) - 2.0*PSRC(:,IS-2,:) +     PSRC(:,IS-1,:))**2 & 
        + 1./4    * (PSRC(:,IS-3,:) - 4.0*PSRC(:,IS-2,:) + 3.0*PSRC(:,IS-1,:))**2  ! Smoothness indicator IS-1
 ZOMP1(:,IS-1:IS,:)  = 1./10. / (ZEPS + ZBPOS1(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(:,IS,:)   = 1./6.* (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS
 ZFPOS2(:,IS-1,:) = 1./6.* (-1.0*PSRC(:,IS-2,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS-1
 ZBPOS2(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZBPOS2(:,IS-1,:) = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IS-2,:) -                      PSRC(:,IS,:))**2  ! Smoothness indicator IS-1
 ZOMP2(:,IS-1:IS,:)  = 3./5. / (ZEPS + ZBPOS2(:,IS-1:IS,:))**2  ! Non-normalized weight IS,IS-1
!
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(:,IS,:)   = 1./6 * (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS
 ZFNEG3(:,IS-1,:) = 1./6 * (-1.0*PSRC(:,IS-2,:)   + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS-1
 ZBNEG3(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) +     PSRC(:,IS+1,:))**2 &
        + 1./4   * (PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + 3.0*PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZBNEG3(:,IS-1,:) = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 &
        + 1./4   * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2 ! Smoothness indicator IS-1
 ZOMN3(:,IS-1:IS,:) = 3./10. / (ZEPS + ZBNEG3(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
! 
 ELSE ! North boundary is proc border, with NHALO < 3 on south side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on south side'
  CALL ABORT
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IS,:)   = 1./6 * (2.0*PSRC(:,IS,:) + 5.0*PSRC(:,IS+1,:) - PSRC(:,IS+2,:)) ! Flux IS
 ZFPOS3(:,IS-1,:) = 1./6 * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:) - PSRC(:,IS+1,:)) ! Flux IS-1
 ZBPOS3(:,IS,:)   = 13./12 * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 ! Smoothness indicator IS
 ZBPOS3(:,IS-1,:) = 13./12 * (    PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 ! Smoothness indicator IS-1
 ZOMP3(:,IS-1:IS,:)  = 3./10. / (ZEPS + ZBPOS3(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(:,IS,:)   = 1./6. * (11.0*PSRC(:,IS+1,:) - 7.0*PSRC(:,IS+2,:) + 2.0*PSRC(:,IS+3,:)) ! Flux IS
 ZFNEG1(:,IS-1,:) = 1./6. * (11.0*PSRC(:,IS,:)   - 7.0*PSRC(:,IS+1,:) + 2.0*PSRC(:,IS+2,:)) ! Flux IS-1
 ZBNEG1(:,IS,:)   = 13./12. * (    PSRC(:,IS+1,:) - 2.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS+1,:) - 4.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2  ! Smoothness indicator IS
 ZBNEG1(:,IS-1,:) = 13./12. * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2  ! Smoothness indicator IS-1
 ZOMN1(:,IS-1:IS,:) = 1./10. / (ZEPS + ZBNEG1(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS,:)   + 5.0*PSRC(:,IS+1,:) - 1.0*PSRC(:,IS+2,:)) ! Flux IS
 ZFNEG2(:,IS-1,:) = 1./6. * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   - 1.0*PSRC(:,IS+1,:)) ! Flux IS-1
 ZBNEG2(:,IS,:)   = 13./12 * (PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (PSRC(:,IS,:) -                      PSRC(:,IS+2,:))**2 ! Smoothness indicator IS
 ZBNEG2(:,IS-1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS-1
 ZOMN2(:,IS-1:IS,:) = 3./5. / (ZEPS + ZBNEG2(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! ----- Total flux -----
! 
 PR(:,IS-1:IS,:) = ( ZOMP1(:,IS-1:IS,:)/(ZOMP1(:,IS-1:IS,:)+ZOMP2(:,IS-1:IS,:)+ZOMP3(:,IS-1:IS,:)) * ZFPOS1(:,IS-1:IS,:) &
                       + ZOMP2(:,IS-1:IS,:)/(ZOMP1(:,IS-1:IS,:)+ZOMP2(:,IS-1:IS,:)+ZOMP3(:,IS-1:IS,:)) * ZFPOS2(:,IS-1:IS,:) & 
                       + ZOMP3(:,IS-1:IS,:)/(ZOMP1(:,IS-1:IS,:)+ZOMP2(:,IS-1:IS,:)+ZOMP3(:,IS-1:IS,:)) * ZFPOS3(:,IS-1:IS,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IS-1:IS,:))) &
                    + ( ZOMN1(:,IS-1:IS,:)/(ZOMN1(:,IS-1:IS,:)+ZOMN2(:,IS-1:IS,:)+ZOMN3(:,IS-1:IS,:)) * ZFNEG1(:,IS-1:IS,:) &
                       + ZOMN2(:,IS-1:IS,:)/(ZOMN1(:,IS-1:IS,:)+ZOMN2(:,IS-1:IS,:)+ZOMN3(:,IS-1:IS,:)) * ZFNEG2(:,IS-1:IS,:) &
                       + ZOMN3(:,IS-1:IS,:)/(ZOMN1(:,IS-1:IS,:)+ZOMN2(:,IS-1:IS,:)+ZOMN3(:,IS-1:IS,:)) * ZFNEG3(:,IS-1:IS,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IS-1:IS,:)))
!
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IS-1
    PR(:,IS-1,:) = PSRC(:,IS-1,:) * (0.5+SIGN(0.5,PRVCT(:,IS-1,:))) + &
                   PSRC(:,IS,:) * &
                   (0.5-SIGN(0.5,PRVCT(:,IS-1,:)))
!
!   ! WENO scheme order 3, IS
    ZFPOS1(:,IS,:) = 0.5 * (3.0*PSRC(:,IS,:) - PSRC(:,IS-1,:)) ! First positive flux
    ZFPOS2(:,IS,:) = 0.5 * ( PSRC(:,IS,:) + PSRC(:,IS+1,:)) ! Second positive flux
    ZBPOS1(:,IS,:) = (PSRC(:,IS,:)   - PSRC(:,IS-1,:))**2 ! First positive smoothness indicator
    ZBPOS2(:,IS,:) = (PSRC(:,IS+1,:) - PSRC(:,IS,:))**2  ! Second positive smoothness indicator
!
    ZFNEG1(:,IS,:) = 0.5 * (3.0*PSRC(:,IS+1,:) - PSRC(:,IS+2,:)) ! First negative flux
    ZFNEG2(:,IS,:) = 0.5 * ( PSRC(:,IS,:)   + PSRC(:,IS+1,:)) ! Second negative flux
    ZBNEG1(:,IS,:) = (PSRC(:,IS+1,:) - PSRC(:,IS+2,:))**2 ! First negative smoothness indicator
    ZBNEG2(:,IS,:) = (PSRC(:,IS,:)   - PSRC(:,IS+1,:))**2 ! Second negative smoothness indicator
!
    ZOMP1(:,IS,:) = 1./3. / (ZEPS + ZBPOS1(:,IS,:))**2 ! First positive non-normalized weight
    ZOMP2(:,IS,:) = 2./3. / (ZEPS + ZBPOS2(:,IS,:))**2 ! Second positive non-normalized weight
    ZOMN1(:,IS,:) = 1./3. / (ZEPS + ZBNEG1(:,IS,:))**2 ! First negative non-normalized weight
    ZOMN2(:,IS,:) = 2./3. / (ZEPS + ZBNEG2(:,IS,:))**2 ! Second negative non-normalized weight
! 
    PR(:,IS,:) = (ZOMN2(:,IS,:)/(ZOMN1(:,IS,:)+ZOMN2(:,IS,:)) * ZFNEG2(:,IS,:) + &
             (ZOMN1(:,IS,:)/(ZOMN1(:,IS,:)+ZOMN2(:,IS,:)) * ZFNEG1(:,IS,:))) &
           *(0.5-SIGN(0.5,PRVCT(:,IS,:))) + &
     (ZOMP2(:,IS,:)/(ZOMP1(:,IS,:)+ZOMP2(:,IS,:)) * ZFPOS2(:,IS,:) + &
     (ZOMP1(:,IS,:)/(ZOMP1(:,IS,:)+ZOMP2(:,IS,:)) * ZFPOS1(:,IS,:))) &
     *(0.5+SIGN(0.5,PRVCT(:,IS,:)))  ! Total flux
! 
 END SELECT ! SELECT CASE (HLBCY(1)) ! Y direction LBC type on south side
!
ELSE
 !-----------------------------------------------------------------------------
 ! South border is proc border -- IS,IS-1
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/south-int not parallelisable with NHALO < 3' 
 CALL ABORT
 STOP  ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
!
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS-2,:) - 7.0*PSRC(:,IS-1,:) + 11.0*PSRC(:,IS,:)) ! Flux IS
 ZFPOS1(:,IS-1,:) = 1./6. * (2.0*PSRC(:,IS-3,:) - 7.0*PSRC(:,IS-2,:) + 11.0*PSRC(:,IS-1,:)) ! Flux IS-1
 ZBPOS1(:,IS,:)   = 13./12. * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 & 
        + 1./4    * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2   ! Smoothness indicator IS
 ZBPOS1(:,IS-1,:) = 13./12. * (PSRC(:,IS-3,:) - 2.0*PSRC(:,IS-2,:) +     PSRC(:,IS-1,:))**2 & 
        + 1./4    * (PSRC(:,IS-3,:) - 4.0*PSRC(:,IS-2,:) + 3.0*PSRC(:,IS-1,:))**2  ! Smoothness indicator IS-1
 ZOMP1(:,IS-1:IS,:)  = 1./10. / (ZEPS + ZBPOS1(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(:,IS,:)   = 1./6.* (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS
 ZFPOS2(:,IS-1,:) = 1./6.* (-1.0*PSRC(:,IS-2,:) + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS-1
 ZBPOS2(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZBPOS2(:,IS-1,:) = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IS-2,:) -                      PSRC(:,IS,:))**2 ! Smoothness indicator IS-1
 ZOMP2(:,IS-1:IS,:)  = 3./5. / (ZEPS + ZBPOS2(:,IS-1:IS,:))**2  ! Non-normalized weight IS,IS-1
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IS,:)   = 1./6 * (2.0*PSRC(:,IS,:) + 5.0*PSRC(:,IS+1,:) - PSRC(:,IS+2,:)) ! Flux IS
 ZFPOS3(:,IS-1,:) = 1./6 * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:) - PSRC(:,IS+1,:)) ! Flux IS-1
 ZBPOS3(:,IS,:)   = 13./12 * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 ! Smoothness indicator IS
 ZBPOS3(:,IS-1,:) = 13./12 * (    PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 ! Smoothness indicator IS-1
 ZOMP3(:,IS-1:IS,:)  = 3./10. / (ZEPS + ZBPOS3(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(:,IS,:)   = 1./6. * (11.0*PSRC(:,IS+1,:) - 7.0*PSRC(:,IS+2,:) + 2.0*PSRC(:,IS+3,:)) ! Flux IS
 ZFNEG1(:,IS-1,:) = 1./6. * (11.0*PSRC(:,IS,:)   - 7.0*PSRC(:,IS+1,:) + 2.0*PSRC(:,IS+2,:)) ! Flux IS-1
 ZBNEG1(:,IS,:)   = 13./12. * (    PSRC(:,IS+1,:) - 2.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS+1,:) - 4.0*PSRC(:,IS+2,:) + PSRC(:,IS+3,:))**2  ! Smoothness indicator IS
 ZBNEG1(:,IS-1,:) = 13./12. * (    PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IS,:) - 4.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2  ! Smoothness indicator IS-1
 ZOMN1(:,IS-1:IS,:) = 1./10. / (ZEPS + ZBNEG1(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(:,IS,:)   = 1./6. * (2.0*PSRC(:,IS,:)   + 5.0*PSRC(:,IS+1,:) - 1.0*PSRC(:,IS+2,:)) ! Flux IS
 ZFNEG2(:,IS-1,:) = 1./6. * (2.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   - 1.0*PSRC(:,IS+1,:)) ! Flux IS-1
 ZBNEG2(:,IS,:)   = 13./12 * (PSRC(:,IS,:) - 2.0*PSRC(:,IS+1,:) + PSRC(:,IS+2,:))**2 &
      + 1./4   * (PSRC(:,IS,:) -                      PSRC(:,IS+2,:))**2 ! Smoothness indicator IS
 ZBNEG2(:,IS-1,:) = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 &
      + 1./4   * (PSRC(:,IS-1,:) -                    PSRC(:,IS+1,:))**2 ! Smoothness indicator IS-1
 ZOMN2(:,IS-1:IS,:) = 3./5. / (ZEPS + ZBNEG2(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(:,IS,:)   = 1./6 * (-1.0*PSRC(:,IS-1,:) + 5.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IS
 ZFNEG3(:,IS-1,:) = 1./6 * (-1.0*PSRC(:,IS-2,:) + 5.0*PSRC(:,IS-1,:) + 2.0*PSRC(:,IS,:)) ! Flux IS-1
 ZBNEG3(:,IS,:)   = 13./12 * (PSRC(:,IS-1,:) - 2.0*PSRC(:,IS,:) +     PSRC(:,IS+1,:))**2 &
        + 1./4   * (PSRC(:,IS-1,:) - 4.0*PSRC(:,IS,:) + 3.0*PSRC(:,IS+1,:))**2 ! Smoothness indicator IS
 ZBNEG3(:,IS-1,:) = 13./12 * (PSRC(:,IS-2,:) - 2.0*PSRC(:,IS-1,:) +     PSRC(:,IS,:))**2 &
        + 1./4   * (PSRC(:,IS-2,:) - 4.0*PSRC(:,IS-1,:) + 3.0*PSRC(:,IS,:))**2 ! Smoothness indicator IS-1
 ZOMN3(:,IS-1:IS,:) = 3./10. / (ZEPS + ZBNEG3(:,IS-1:IS,:))**2 ! Non-normalized weight IS,IS-1
!
 ! ----- Total flux -----
! 
 PR(:,IS-1:IS,:) = ( ZOMP1(:,IS-1:IS,:)/(ZOMP1(:,IS-1:IS,:)+ZOMP2(:,IS-1:IS,:)+ZOMP3(:,IS-1:IS,:)) * ZFPOS1(:,IS-1:IS,:) &
                       + ZOMP2(:,IS-1:IS,:)/(ZOMP1(:,IS-1:IS,:)+ZOMP2(:,IS-1:IS,:)+ZOMP3(:,IS-1:IS,:)) * ZFPOS2(:,IS-1:IS,:) & 
                       + ZOMP3(:,IS-1:IS,:)/(ZOMP1(:,IS-1:IS,:)+ZOMP2(:,IS-1:IS,:)+ZOMP3(:,IS-1:IS,:)) * ZFPOS3(:,IS-1:IS,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IS-1:IS,:))) &
                    + ( ZOMN1(:,IS-1:IS,:)/(ZOMN1(:,IS-1:IS,:)+ZOMN2(:,IS-1:IS,:)+ZOMN3(:,IS-1:IS,:)) * ZFNEG1(:,IS-1:IS,:) &
                       + ZOMN2(:,IS-1:IS,:)/(ZOMN1(:,IS-1:IS,:)+ZOMN2(:,IS-1:IS,:)+ZOMN3(:,IS-1:IS,:)) * ZFNEG2(:,IS-1:IS,:) &
                       + ZOMN3(:,IS-1:IS,:)/(ZOMN1(:,IS-1:IS,:)+ZOMN2(:,IS-1:IS,:)+ZOMN3(:,IS-1:IS,:)) * ZFNEG3(:,IS-1:IS,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IS-1:IS,:)))
!
 END IF ! NHALO
!
END IF ! IF(LSOUTH_ll()) 
!
!-------------------------------------------------------------------------------
!*       1.3.   North border
!               ---------------------
!
!! IF(LNORTH_ll()  .AND. .FALSE. ) THEN 
IF(LNORTH_ll()) THEN 
 !-----------------------------------------------------------------------------
 ! North border is physical -- IN-1,IN
 !-----------------------------------------------------------------------------
 SELECT CASE (HLBCY(2)) ! Y direction LBC type on north side
! 
 CASE ('CYCL')
 !---------------------------------------------------------------------------
 ! Periodic boundary condition
 !---------------------------------------------------------------------------
! 
 IF(LSOUTH_ll()  .AND. .FALSE. ) THEN  ! South border is physical 
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IN-1,:) = 1./6 * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:) - PSRC(:,IN+1,:)) ! Flux IN-1
 ZFPOS3(:,IN,:)   = 1./6 * (2.0*PSRC(:,IN,:) + 5.0*PSRC(:,IN+1,:) - PSRC(:,IS,:)) ! Flux IN
 ZBPOS3(:,IN-1,:) = 13./12 * (    PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 ! Smoothness indicator IN-1
 ZBPOS3(:,IN,:)   = 13./12 * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2  ! Smoothness indicator IN
 ZOMP3(:,IN-1:IN,:)  = 3./10. / (ZEPS + ZBPOS3(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(:,IN-1,:) = 1./6. * (11.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IN+1,:) + 2.0*PSRC(:,IS,:)) ! Flux IN-1
 ZFNEG1(:,IN,:)   = 1./6. * (11.0*PSRC(:,IN+1,:) - 7.0*PSRC(:,IS,:)   + 2.0*PSRC(:,IS+1,:)) ! Flux IN
 ZBNEG1(:,IN-1,:) = 13./12. * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2  ! Smoothness indicator IN-1
 ZBNEG1(:,IN,:)   = 13./12. * (    PSRC(:,IN+1,:) - 2.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN+1,:) - 4.0*PSRC(:,IS,:) + PSRC(:,IS+1,:))**2  ! Smoothness indicator IN
 ZOMN1(:,IN-1:IN,:) = 1./10. / (ZEPS + ZBNEG1(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(:,IN-1,:) = 1./6. * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   - 1.0*PSRC(:,IN+1,:)) ! Flux IN-1
 ZFNEG2(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IN+1,:) - 1.0*PSRC(:,IS,:)) ! Flux IN
 ZBNEG2(:,IN-1,:)  = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                    PSRC(:,IN+1,:))**2 ! Smoothness indicator IN-1
 ZBNEG2(:,IN,:)   = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IS,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IS,:))**2  ! Smoothness indicator IN
 ZOMN2(:,IN-1:IN,:) = 3./5. / (ZEPS + ZBNEG2(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ELSEIF(IN<=SIZE(PSRC,2)-3) THEN ! South boundary is proc border, with minimum 3 HALO points on north side
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IN-1,:) = 1./6 * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:) - PSRC(:,IN+1,:)) ! Flux IN-1
 ZFPOS3(:,IN,:)   = 1./6 * (2.0*PSRC(:,IN,:) + 5.0*PSRC(:,IN+1,:) - PSRC(:,IN+2,:)) ! Flux IN
 ZBPOS3(:,IN-1,:) = 13./12 * (    PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 ! Smoothness indicator IN-1
 ZBPOS3(:,IN,:)   = 13./12 * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 ! Smoothness indicator IN
 ZOMP3(:,IN-1:IN,:)  = 3./10. / (ZEPS + ZBPOS3(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(:,IN-1,:) = 1./6. * (11.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IN+1,:) + 2.0*PSRC(:,IN+2,:)) ! Flux IN-1
 ZFNEG1(:,IN,:)   = 1./6. * (11.0*PSRC(:,IN+1,:) - 7.0*PSRC(:,IN+2,:)   + 2.0*PSRC(:,IN+3,:)) ! Flux IN
 ZBNEG1(:,IN-1,:) = 13./12. * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2  ! Smoothness indicator IN-1
 ZBNEG1(:,IN,:)   = 13./12. * (    PSRC(:,IN+1,:) - 2.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN+1,:) - 4.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2  ! Smoothness indicator IN
 ZOMN1(:,IN-1:IN,:) = 1./10. / (ZEPS + ZBNEG1(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(:,IN-1,:) = 1./6. * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   - 1.0*PSRC(:,IN+1,:)) ! Flux IN-1
 ZFNEG2(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IN+1,:) - 1.0*PSRC(:,IN+2,:)) ! Flux IN
 ZBNEG2(:,IN-1,:)  = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                    PSRC(:,IN+1,:))**2 ! Smoothness indicator IN-1
 ZBNEG2(:,IN,:)   = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IN+2,:))**2 ! Smoothness indicator IN
 ZOMN2(:,IN-1:IN,:) = 3./5. / (ZEPS + ZBNEG2(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ELSE ! South boundary is proc border, with NHALO < 3 on north side
  PRINT *,'ERROR : WENO5/CYCL fluxes calculation needs JPHEXT (&NHALO) >= 3 on north side'
  CALL ABORT
  STOP  ' Error in advec_weno_k_3_aux.f90 '
 ENDIF
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(:,IN-1,:) = 1./6. * (2.0*PSRC(:,IN-3,:) - 7.0*PSRC(:,IN-2,:) + 11.0*PSRC(:,IN-1,:)) ! Flux IN-1
 ZFPOS1(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN-2,:) - 7.0*PSRC(:,IN-1,:) + 11.0*PSRC(:,IN,:)) ! Flux IN
 ZBPOS1(:,IN-1,:) = 13./12. * (PSRC(:,IN-3,:) - 2.0*PSRC(:,IN-2,:) +     PSRC(:,IN-1,:))**2 & 
        + 1./4    * (PSRC(:,IN-3,:) - 4.0*PSRC(:,IN-2,:) + 3.0*PSRC(:,IN-1,:))**2  ! Smoothness indicator IN-1
 ZBPOS1(:,IN,:)   = 13./12. * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 & 
        + 1./4    * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2  ! Smoothness indicator IN
 ZOMP1(:,IN-1:IN,:)  = 1./10. / (ZEPS + ZBPOS1(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(:,IN-1,:) = 1./6. * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN-1
 ZFPOS2(:,IN,:)   = 1./6. * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN
 ZBPOS2(:,IN-1,:) = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) + PSRC(:,IN,:))**2 &
      + 1./4   * (PSRC(:,IN-2,:) -                      PSRC(:,IN,:))**2 ! Smoothness indicator IN-1
 ZBPOS2(:,IN,:)   = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                   PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZOMP2(:,IN-1:IN,:)  = 3./5. / (ZEPS + ZBPOS2(:,IN-1:IN,:))**2  ! Non-normalized weight IN-1,IN
!
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(:,IN-1,:) = 1./6 * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN-1
 ZFNEG3(:,IN,:)   = 1./6 * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN
 ZBNEG3(:,IN-1,:) = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 &
        + 1./4   * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2 ! Smoothness indicator IN-1
 ZBNEG3(:,IN,:)   = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) +     PSRC(:,IN+1,:))**2 &
        + 1./4   * (PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + 3.0*PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZOMN3(:,IN-1:IN,:) = 3./10. / (ZEPS + ZBNEG3(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! ----- Total flux -----
! 
 PR(:,IN-1:IN,:) = ( ZOMP1(:,IN-1:IN,:)/(ZOMP1(:,IN-1:IN,:)+ZOMP2(:,IN-1:IN,:)+ZOMP3(:,IN-1:IN,:)) * ZFPOS1(:,IN-1:IN,:) &
                       + ZOMP2(:,IN-1:IN,:)/(ZOMP1(:,IN-1:IN,:)+ZOMP2(:,IN-1:IN,:)+ZOMP3(:,IN-1:IN,:)) * ZFPOS2(:,IN-1:IN,:) & 
                       + ZOMP3(:,IN-1:IN,:)/(ZOMP1(:,IN-1:IN,:)+ZOMP2(:,IN-1:IN,:)+ZOMP3(:,IN-1:IN,:)) * ZFPOS3(:,IN-1:IN,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IN-1:IN,:))) &
                    + ( ZOMN1(:,IN-1:IN,:)/(ZOMN1(:,IN-1:IN,:)+ZOMN2(:,IN-1:IN,:)+ZOMN3(:,IN-1:IN,:)) * ZFNEG1(:,IN-1:IN,:) &
                       + ZOMN2(:,IN-1:IN,:)/(ZOMN1(:,IN-1:IN,:)+ZOMN2(:,IN-1:IN,:)+ZOMN3(:,IN-1:IN,:)) * ZFNEG2(:,IN-1:IN,:) &
                       + ZOMN3(:,IN-1:IN,:)/(ZOMN1(:,IN-1:IN,:)+ZOMN2(:,IN-1:IN,:)+ZOMN3(:,IN-1:IN,:)) * ZFNEG3(:,IN-1:IN,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IN-1:IN,:)))
!
!
 CASE ('OPEN','WALL','NEST') 
 !---------------------------------------------------------------------------
 ! Open, or Wall, or Nest boundary condition => WENO order reduction
 !---------------------------------------------------------------------------
!
 ! WENO scheme order 1, IN
    PR(:,IN,:) = PSRC(:,IN,:) * (0.5+SIGN(0.5,PRVCT(:,IN,:))) + &
                 PSRC(:,IN+1,:) * &
                 (0.5-SIGN(0.5,PRVCT(:,IN,:)))
!
!   ! WENO scheme order 3, IN-1
    ZFPOS1(:,IN-1,:) = 0.5 * (3.0*PSRC(:,IN-1,:) - PSRC(:,IN-2,:)) ! First positive flux
    ZFPOS2(:,IN-1,:) = 0.5 * (    PSRC(:,IN-1,:) + PSRC(:,IN,:)) ! Second positive flux
    ZBPOS1(:,IN-1,:) = (PSRC(:,IN-1,:) - PSRC(:,IN-2,:))**2 ! First positive smoothness indicator
    ZBPOS2(:,IN-1,:) = (PSRC(:,IN,:)   - PSRC(:,IN-1,:))**2 ! Second positive smoothness indicator
!
    ZFNEG1(:,IN-1,:) = 0.5 * (3.0*PSRC(:,IN,:)   - PSRC(:,IN+1,:)) ! First negative flux
    ZFNEG2(:,IN-1,:) = 0.5 * (    PSRC(:,IN-1,:) + PSRC(:,IN,:)) ! Second negative flux
    ZBNEG1(:,IN-1,:) = (PSRC(:,IN,:)   - PSRC(:,IN+1,:))**2 ! First negative smoothness indicator
    ZBNEG2(:,IN-1,:) = (PSRC(:,IN-1,:) - PSRC(:,IN,:))**2  ! Second negative smoothness indicator
!
    ZOMP1(:,IN-1,:) = 1./3. / (ZEPS + ZBPOS1(:,IN-1,:))**2 ! First positive non-normalized weight
    ZOMP2(:,IN-1,:) = 2./3. / (ZEPS + ZBPOS2(:,IN-1,:))**2 ! Second positive non-normalized weight
    ZOMN1(:,IN-1,:) = 1./3. / (ZEPS + ZBNEG1(:,IN-1,:))**2 ! First negative non-normalized weight
    ZOMN2(:,IN-1,:) = 2./3. / (ZEPS + ZBNEG2(:,IN-1,:))**2 ! Second negative non-normalized weight
! 
    PR(:,IN-1,:) = (ZOMN2(:,IN-1,:)/(ZOMN1(:,IN-1,:)+ZOMN2(:,IN-1,:))*ZFNEG2(:,IN-1,:) + &
      (ZOMN1(:,IN-1,:)/(ZOMN1(:,IN-1,:)+ZOMN2(:,IN-1,:))*ZFNEG1(:,IN-1,:))) &
      *(0.5-SIGN(0.5,PRVCT(:,IN-1,:))) + &
                 (ZOMP2(:,IN-1,:)/(ZOMP1(:,IN-1,:)+ZOMP2(:,IN-1,:))*ZFPOS2(:,IN-1,:) + &
                   (ZOMP1(:,IN-1,:)/(ZOMP1(:,IN-1,:)+ZOMP2(:,IN-1,:))*ZFPOS1(:,IN-1,:))) &
      * (0.5+SIGN(0.5,PRVCT(:,IN-1,:)))  ! Total flux
! 
 END SELECT ! SELECT CASE (HLBCY(2)) ! Y direction LBC type on north side
!
ELSE
 !-----------------------------------------------------------------------------
 ! North border is proc border -- IN-1,IN
 !-----------------------------------------------------------------------------
!
 IF (NHALO<3) THEN
 PRINT *,'ERROR : WENO5/north-int not parallelisable with NHALO < 3' 
 CALL ABORT
 STOP  ' Error in advec_weno_k_3_aux.f90 '
 ELSEIF (NHALO>=3) THEN
 !---------------------------------------------------------------------------
 ! NHALO >= 3 => WENO5 for all boundary points
 !---------------------------------------------------------------------------
! 
 ! ----- Positive fluxes -----
!
 ! First positive stencil, needs indices i-2, i-1, i 
 ZFPOS1(:,IN-1,:) = 1./6. * (2.0*PSRC(:,IN-3,:) - 7.0*PSRC(:,IN-2,:) + 11.0*PSRC(:,IN-1,:)) ! Flux IN-1
 ZFPOS1(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN-2,:) - 7.0*PSRC(:,IN-1,:) + 11.0*PSRC(:,IN,:)) ! Flux IN
 ZBPOS1(:,IN-1,:) = 13./12. * (PSRC(:,IN-3,:) - 2.0*PSRC(:,IN-2,:) +     PSRC(:,IN-1,:))**2 & 
        + 1./4    * (PSRC(:,IN-3,:) - 4.0*PSRC(:,IN-2,:) + 3.0*PSRC(:,IN-1,:))**2  ! Smoothness indicator IN-1
 ZBPOS1(:,IN,:)   = 13./12. * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 & 
        + 1./4    * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2  ! Smoothness indicator IN
 ZOMP1(:,IN-1:IN,:)  = 1./10. / (ZEPS + ZBPOS1(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! Second positive stencil, needs indices i-1, i, i+1
 ZFPOS2(:,IN-1,:) = 1./6. * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN-1
 ZFPOS2(:,IN,:)   = 1./6. * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN
 ZBPOS2(:,IN-1,:) = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) + PSRC(:,IN,:))**2 &
      + 1./4   * (PSRC(:,IN-2,:) -                      PSRC(:,IN,:))**2 ! Smoothness indicator IN-1
 ZBPOS2(:,IN,:)   = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                   PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZOMP2(:,IN-1:IN,:)  = 3./5. / (ZEPS + ZBPOS2(:,IN-1:IN,:))**2  ! Non-normalized weight IN-1,IN
! 
 ! Third positive stencil, needs indices i, i+1, i+2
 ZFPOS3(:,IN-1,:) = 1./6 * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:) - PSRC(:,IN+1,:)) ! Flux IN-1
 ZFPOS3(:,IN,:)   = 1./6 * (2.0*PSRC(:,IN,:) + 5.0*PSRC(:,IN+1,:) - PSRC(:,IN+2,:)) ! Flux IN
 ZBPOS3(:,IN-1,:) = 13./12 * (    PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 ! Smoothness indicator IN-1
 ZBPOS3(:,IN,:)   = 13./12 * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2  ! Smoothness indicator IN
 ZOMP3(:,IN-1:IN,:)  = 3./10. / (ZEPS + ZBPOS3(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! ----- Negative fluxes ----- 
!
 ! First negative stencil, needs indices i+1, i+2, i+3
 ZFNEG1(:,IN-1,:) = 1./6. * (11.0*PSRC(:,IN,:)   - 7.0*PSRC(:,IN+1,:) + 2.0*PSRC(:,IN+2,:)) ! Flux IN-1
 ZFNEG1(:,IN,:)   = 1./6. * (11.0*PSRC(:,IN+1,:) - 7.0*PSRC(:,IN+2,:) + 2.0*PSRC(:,IN+3,:)) ! Flux IN
 ZBNEG1(:,IN-1,:) = 13./12. * (    PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN,:) - 4.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2  ! Smoothness indicator IN-1
 ZBNEG1(:,IN,:)   = 13./12. * (    PSRC(:,IN+1,:) - 2.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2 & 
      + 1./4    * (3.0*PSRC(:,IN+1,:) - 4.0*PSRC(:,IN+2,:) + PSRC(:,IN+3,:))**2  ! Smoothness indicator IN
 ZOMN1(:,IN-1:IN,:) = 1./10. / (ZEPS + ZBNEG1(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! Second negative stencil, needs indices i, i+1, i+2
 ZFNEG2(:,IN-1,:) = 1./6. * (2.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   - 1.0*PSRC(:,IN+1,:)) ! Flux IN-1
 ZFNEG2(:,IN,:)   = 1./6. * (2.0*PSRC(:,IN,:)   + 5.0*PSRC(:,IN+1,:) - 1.0*PSRC(:,IN+2,:)) ! Flux IN
 ZBNEG2(:,IN-1,:)  = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) + PSRC(:,IN+1,:))**2 &
      + 1./4   * (PSRC(:,IN-1,:) -                    PSRC(:,IN+1,:))**2 ! Smoothness indicator IN-1
 ZBNEG2(:,IN,:)   = 13./12 * (PSRC(:,IN,:) - 2.0*PSRC(:,IN+1,:) + PSRC(:,IN+2,:))**2 &
      + 1./4   * (PSRC(:,IN,:) -                      PSRC(:,IN+2,:))**2 ! Smoothness indicator IN
 ZOMN2(:,IN-1:IN,:) = 3./5. / (ZEPS + ZBNEG2(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! Third negative stencil, needs indices i-1, i, i+1
 ZFNEG3(:,IN-1,:) = 1./6 * (-1.0*PSRC(:,IN-2,:) + 5.0*PSRC(:,IN-1,:) + 2.0*PSRC(:,IN,:)) ! Flux IN-1
 ZFNEG3(:,IN,:)   = 1./6 * (-1.0*PSRC(:,IN-1,:) + 5.0*PSRC(:,IN,:)   + 2.0*PSRC(:,IN+1,:)) ! Flux IN
 ZBNEG3(:,IN-1,:) = 13./12 * (PSRC(:,IN-2,:) - 2.0*PSRC(:,IN-1,:) +     PSRC(:,IN,:))**2 &
        + 1./4   * (PSRC(:,IN-2,:) - 4.0*PSRC(:,IN-1,:) + 3.0*PSRC(:,IN,:))**2 ! Smoothness indicator IN-1
 ZBNEG3(:,IN,:)   = 13./12 * (PSRC(:,IN-1,:) - 2.0*PSRC(:,IN,:) +     PSRC(:,IN+1,:))**2 &
        + 1./4   * (PSRC(:,IN-1,:) - 4.0*PSRC(:,IN,:) + 3.0*PSRC(:,IN+1,:))**2 ! Smoothness indicator IN
 ZOMN3(:,IN-1:IN,:) = 3./10. / (ZEPS + ZBNEG3(:,IN-1:IN,:))**2 ! Non-normalized weight IN-1,IN
!
 ! ----- Total flux -----
! 
 PR(:,IN-1:IN,:) = ( ZOMP1(:,IN-1:IN,:) /(ZOMP1(:,IN-1:IN,:) +ZOMP2(:,IN-1:IN,:) +ZOMP3(:,IN-1:IN,:)) * ZFPOS1(:,IN-1:IN,:) &
                       + ZOMP2(:,IN-1:IN,:) /(ZOMP1(:,IN-1:IN,:) +ZOMP2(:,IN-1:IN,:) +ZOMP3(:,IN-1:IN,:)) * ZFPOS2(:,IN-1:IN,:) & 
                       + ZOMP3(:,IN-1:IN,:) /(ZOMP1(:,IN-1:IN,:) +ZOMP2(:,IN-1:IN,:) +ZOMP3(:,IN-1:IN,:)) * ZFPOS3(:,IN-1:IN,:)) &
                      * (0.5+SIGN(0.5,PRVCT(:,IN-1:IN,:))) &
                    + ( ZOMN1(:,IN-1:IN,:)/(ZOMN1(:,IN-1:IN,:)+ZOMN2(:,IN-1:IN,:)+ZOMN3(:,IN-1:IN,:)) * ZFNEG1(:,IN-1:IN,:) &
                       + ZOMN2(:,IN-1:IN,:)/(ZOMN1(:,IN-1:IN,:)+ZOMN2(:,IN-1:IN,:)+ZOMN3(:,IN-1:IN,:)) * ZFNEG2(:,IN-1:IN,:) &
                       + ZOMN3(:,IN-1:IN,:)/(ZOMN1(:,IN-1:IN,:)+ZOMN2(:,IN-1:IN,:)+ZOMN3(:,IN-1:IN,:)) * ZFNEG3(:,IN-1:IN,:)) &
                      * (0.5-SIGN(0.5,PRVCT(:,IN-1:IN,:)))
! 
 END IF ! NHALO
!
END IF ! IF(LNORTH_ll()) 
!-------------------------------------------------------------------------------
!
PR = PR * PRVCT ! Add contravariant flux
!
END SUBROUTINE ADVEC_WENO_K_3_VY
!
!-------------------------------------------------------------------------------
!
!     ############################################
      FUNCTION WENO_K_3_WZ(PSRC, PRWCT) RESULT(PR)
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
!! T. Lunet 02/10/2014:  Complete code documentation
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_CONF
USE MODD_PARAMETERS,ONLY: JPVEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on W grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
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
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFPOS1, ZFPOS2, ZFPOS3
!
! intermediate reconstruction fluxes for negative wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFNEG1, ZFNEG2, ZFNEG3
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBPOS1, ZBPOS2, ZBPOS3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBNEG1, ZBNEG2, ZBNEG3
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZOMP1, ZOMP2, ZOMP3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZOMN1, ZOMN2, ZOMN3
!
! EPSILON for weno weights calculation
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
ZFPOS3 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZFNEG3 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBPOS3 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZBNEG3 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMP3  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0
ZOMN3  = 0.0 
!
!-------------------------------------------------------------------------------
!*       1.1.   Interior Fluxes 
!              ---------------------
!
!-------------------------------------------------------------------------------
! Flux calculation in the physical domain far enough from the boundary 
! WENO scheme order 5, IB+1 -> IT-2
! Computation at the mass point on Wgrid v(i+1/2,j,k)
!-------------------------------------------------------------------------------
!
! ----- Positive fluxes -----
!
! First positive stencil, needs indices i-2, i-1, i
ZFPOS1(:,:,IB+1:IT-2) = 1./6.   * (2.0*PSRC(:,:,IB-1:IT-4) - 7.0*PSRC(:,:,IB:IT-3) + 11.0*PSRC(:,:,IB+1:IT-2))   ! Flux
ZBPOS1(:,:,IB+1:IT-2) = 13./12. * (    PSRC(:,:,IB-1:IT-4) - 2.0*PSRC(:,:,IB:IT-3) +      PSRC(:,:,IB+1:IT-2))**2 & 
      + 1./4    * (    PSRC(:,:,IB-1:IT-4) - 4.0*PSRC(:,:,IB:IT-3) + 3.0* PSRC(:,:,IB+1:IT-2))**2  ! Smoothness indicator
ZOMP1(:,:,IB+1:IT-2)  = 1./10. /  (ZEPS + ZBPOS1(:,:,IB+1:IT-2))**2             ! Non-normalized weight
!
! Second positive stencil, needs indices i-1, i, i+1
ZFPOS2(:,:,IB+1:IT-2) = 1./6.  * (-1.0*PSRC(:,:,IB:IT-3) + 5.0*PSRC(:,:,IB+1:IT-2) + 2.0*PSRC(:,:,IB+2:IT-1))  ! Flux
ZBPOS2(:,:,IB+1:IT-2) = 13./12 * (     PSRC(:,:,IB:IT-3) - 2.0*PSRC(:,:,IB+1:IT-2) +     PSRC(:,:,IB+2:IT-1))**2 &
      + 1./4   * (     PSRC(:,:,IB:IT-3) -                               PSRC(:,:,IB+2:IT-1))**2 ! Smoothness indicator
ZOMP2(:,:,IB+1:IT-2)  = 3./5. /  (ZEPS + ZBPOS2(:,:,IB+1:IT-2))**2             ! Non-normalized weight
!
! Third positive stencil, needs indices i, i+1, i+2
ZFPOS3(:,:,IB+1:IT-2) = 1./6   * (2.0*PSRC(:,:,IB+1:IT-2) + 5.0*PSRC(:,:,IB+2:IT-1) - PSRC(:,:,IB+3:IT))  ! Flux
ZBPOS3(:,:,IB+1:IT-2) = 13./12 * ( PSRC(:,:,IB+1:IT-2) - 2.0*PSRC(:,:,IB+2:IT-1) + PSRC(:,:,IB+3:IT))**2 &
      + 1./4   * (3.0*PSRC(:,:,IB+1:IT-2) - 4.0*PSRC(:,:,IB+2:IT-1) + PSRC(:,:,IB+3:IT))**2 ! Smoothness indicator
ZOMP3(:,:,IB+1:IT-2)  = 3./10. / (ZEPS + ZBPOS3(:,:,IB+1:IT-2))**2           ! Non-normalized weight
!
! ----- Negative fluxes ----- 
!
! First negative stencil, needs indices i+1, i+2, i+3
ZFNEG1(:,:,IB+1:IT-2) = 1./6.   * (11.0*PSRC(:,:,IB+2:IT-1) - 7.0*PSRC(:,:,IB+3:IT) + 2.0*PSRC(:,:,IB+4:IT+1))   ! Flux
ZBNEG1(:,:,IB+1:IT-2) = 13./12. * (     PSRC(:,:,IB+2:IT-1) - 2.0*PSRC(:,:,IB+3:IT) +     PSRC(:,:,IB+4:IT+1))**2 & 
      + 1./4    * (3.0* PSRC(:,:,IB+2:IT-1) - 4.0*PSRC(:,:,IB+3:IT) +     PSRC(:,:,IB+4:IT+1))**2  ! Smoothness indicator
ZOMN1(:,:,IB+1:IT-2)  = 1./10. /  (ZEPS + ZBNEG1(:,:,IB+1:IT-2))**2             ! Non-normalized weight
!
! Second negative stencil, needs indices i, i+1, i+2
ZFNEG2(:,:,IB+1:IT-2) = 1./6.  * (2.0*PSRC(:,:,IB+1:IT-2) + 5.0*PSRC(:,:,IB+2:IT-1) - 1.0*PSRC(:,:,IB+3:IT))  ! Flux
ZBNEG2(:,:,IB+1:IT-2) = 13./12 * (    PSRC(:,:,IB+1:IT-2) - 2.0*PSRC(:,:,IB+2:IT-1) +     PSRC(:,:,IB+3:IT))**2 &
      + 1./4   * (    PSRC(:,:,IB+1:IT-2) -                               PSRC(:,:,IB+3:IT))**2 ! Smoothness indicator
ZOMN2(:,:,IB+1:IT-2)  = 3./5. /  (ZEPS + ZBNEG2(:,:,IB+1:IT-2))**2            ! Non-normalized weight
!
! Third negative stencil, needs indices i-1, i, i+1
ZFNEG3(:,:,IB+1:IT-2) = 1./6   * (-1.0*PSRC(:,:,IB:IT-3) + 5.0*PSRC(:,:,IB+1:IT-2) + 2.0*PSRC(:,:,IB+2:IT-1))  ! Flux
ZBNEG3(:,:,IB+1:IT-2) = 13./12 * (  PSRC(:,:,IB:IT-3) - 2.0*PSRC(:,:,IB+1:IT-2) +     PSRC(:,:,IB+2:IT-1))**2 &
      + 1./4   * (     PSRC(:,:,IB:IT-3) - 4.0*PSRC(:,:,IB+1:IT-2) + 3.0*PSRC(:,:,IB+2:IT-1))**2 ! Smoothness indicator
ZOMN3(:,:,IB+1:IT-2)  = 3./10. / (ZEPS + ZBNEG3(:,:,IB+1:IT-2))**2             ! Non-normalized weight
!
!
! ----- Total flux -----
!
PR(:,:,IB+1:IT-2) = (ZOMP1(:,:,IB+1:IT-2)/(ZOMP1(:,:,IB+1:IT-2)+ZOMP2(:,:,IB+1:IT-2)+ZOMP3(:,:,IB+1:IT-2)) &
           * ZFPOS1(:,:,IB+1:IT-2) + &
                      ZOMP2(:,:,IB+1:IT-2)/(ZOMP1(:,:,IB+1:IT-2)+ZOMP2(:,:,IB+1:IT-2)+ZOMP3(:,:,IB+1:IT-2)) &
           * ZFPOS2(:,:,IB+1:IT-2) + & 
                      ZOMP3(:,:,IB+1:IT-2)/(ZOMP1(:,:,IB+1:IT-2)+ZOMP2(:,:,IB+1:IT-2)+ZOMP3(:,:,IB+1:IT-2)) &
           * ZFPOS3(:,:,IB+1:IT-2)) &
                    * (0.5+SIGN(0.5,PRWCT(:,:,IB+1:IT-2))) &
                  + (ZOMN1(:,:,IB+1:IT-2)/(ZOMN1(:,:,IB+1:IT-2)+ZOMN2(:,:,IB+1:IT-2)+ZOMN3(:,:,IB+1:IT-2)) &
           * ZFNEG1(:,:,IB+1:IT-2)  &
                     + ZOMN2(:,:,IB+1:IT-2)/(ZOMN1(:,:,IB+1:IT-2)+ZOMN2(:,:,IB+1:IT-2)+ZOMN3(:,:,IB+1:IT-2)) &
           * ZFNEG2(:,:,IB+1:IT-2)  &
                     + ZOMN3(:,:,IB+1:IT-2)/(ZOMN1(:,:,IB+1:IT-2)+ZOMN2(:,:,IB+1:IT-2)+ZOMN3(:,:,IB+1:IT-2)) &
           * ZFNEG3(:,:,IB+1:IT-2))  & 
                    * (0.5-SIGN(0.5,PRWCT(:,:,IB+1:IT-2)))
!
!-------------------------------------------------------------------------------
!*       1.2.   Bottom border
!               ---------------------
!---------------------------------------------------------------------------
! WENO order reduction
!---------------------------------------------------------------------------
!
! WENO scheme order 1, IB-1
PR(:,:,IB-1) =  PSRC(:,:,IB-1) * (0.5+SIGN(0.5,PRWCT(:,:,IB-1) )) &
              + PSRC(:,:,IB  ) * (0.5-SIGN(0.5,PRWCT(:,:,IB-1) ))
!
!   ! WENO scheme order 3, IB
ZFPOS1(:,:,IB) = 0.5 * (3.0*PSRC(:,:,IB) - PSRC(:,:,IB-1)) ! First positive flux
ZFPOS2(:,:,IB) = 0.5 * ( PSRC(:,:,IB) + PSRC(:,:,IB+1)) ! Second positive flux
ZBPOS1(:,:,IB) = (PSRC(:,:,IB)   - PSRC(:,:,IB-1))**2 ! First positive smoothness indicator
ZBPOS2(:,:,IB) = (PSRC(:,:,IB+1) - PSRC(:,:,IB))**2  ! Second positive smoothness indicator
!
ZFNEG1(:,:,IB) = 0.5 * (3.0*PSRC(:,:,IB+1) - PSRC(:,:,IB+2)) ! First negative flux
ZFNEG2(:,:,IB) = 0.5 * ( PSRC(:,:,IB)   + PSRC(:,:,IB+1)) ! Second negative flux
ZBNEG1(:,:,IB) = (PSRC(:,:,IB+1) - PSRC(:,:,IB+2))**2 ! First negative smoothness indicator
ZBNEG2(:,:,IB) = (PSRC(:,:,IB)   - PSRC(:,:,IB+1))**2 ! Second negative smoothness indicator
!
ZOMP1(:,:,IB) = 1./3. / (ZEPS + ZBPOS1(:,:,IB))**2 ! First positive non-normalized weight
ZOMP2(:,:,IB) = 2./3. / (ZEPS + ZBPOS2(:,:,IB))**2 ! Second positive non-normalized weight
ZOMN1(:,:,IB) = 1./3. / (ZEPS + ZBNEG1(:,:,IB))**2 ! First negative non-normalized weight
ZOMN2(:,:,IB) = 2./3. / (ZEPS + ZBNEG2(:,:,IB))**2 ! Second negative non-normalized weight
! 
PR(:,:,IB) = (ZOMN2(:,:,IB)/(ZOMN1(:,:,IB)+ZOMN2(:,:,IB)) * ZFNEG2(:,:,IB) + &
         (ZOMN1(:,:,IB)/(ZOMN1(:,:,IB)+ZOMN2(:,:,IB)) * ZFNEG1(:,:,IB))) &
       *(0.5-SIGN(0.5,PRWCT(:,:,IB))) + &
    (ZOMP2(:,:,IB)/(ZOMP1(:,:,IB)+ZOMP2(:,:,IB)) * ZFPOS2(:,:,IB) + &
    (ZOMP1(:,:,IB)/(ZOMP1(:,:,IB)+ZOMP2(:,:,IB)) * ZFPOS1(:,:,IB))) &
    *(0.5+SIGN(0.5,PRWCT(:,:,IB)))  ! Total flux
! 
!-------------------------------------------------------------------------------
!*       1.3.   Top border
!               ---------------------
!
!---------------------------------------------------------------------------
! Open, or Wall, or Nest boundary condition => WENO order reduction
!---------------------------------------------------------------------------
!
! WENO scheme order 1, IT
PR(:,:,IT) = PSRC(:,:,IT  ) * (0.5+SIGN(0.5,PRWCT(:,:,IT) )) &
           + PSRC(:,:,IT+1) * (0.5-SIGN(0.5,PRWCT(:,:,IT) ))
!
!   ! WENO scheme order 3, IT-1
ZFPOS1(:,:,IT-1) = 0.5 * (3.0*PSRC(:,:,IT-1) - PSRC(:,:,IT-2)) ! First positive flux
ZFPOS2(:,:,IT-1) = 0.5 * ( PSRC(:,:,IT-1) + PSRC(:,:,IT)) ! Second positive flux
ZBPOS1(:,:,IT-1) = (PSRC(:,:,IT-1) - PSRC(:,:,IT-2))**2 ! First positive smoothness indicator
ZBPOS2(:,:,IT-1) = (PSRC(:,:,IT)   - PSRC(:,:,IT-1))**2 ! Second positive smoothness indicator
!
ZFNEG1(:,:,IT-1) = 0.5 * (3.0*PSRC(:,:,IT)   - PSRC(:,:,IT+1)) ! First negative flux
ZFNEG2(:,:,IT-1) = 0.5 * ( PSRC(:,:,IT-1) + PSRC(:,:,IT)) ! Second negative flux
ZBNEG1(:,:,IT-1) = (PSRC(:,:,IT)   - PSRC(:,:,IT+1))**2 ! First negative smoothness indicator
ZBNEG2(:,:,IT-1) = (PSRC(:,:,IT-1) - PSRC(:,:,IT))**2  ! Second negative smoothness indicator
!
ZOMP1(:,:,IT-1) = 1./3. / (ZEPS + ZBPOS1(:,:,IT-1))**2 ! First positive non-normalized weight
ZOMP2(:,:,IT-1) = 2./3. / (ZEPS + ZBPOS2(:,:,IT-1))**2 ! Second positive non-normalized weight
ZOMN1(:,:,IT-1) = 1./3. / (ZEPS + ZBNEG1(:,:,IT-1))**2 ! First negative non-normalized weight
ZOMN2(:,:,IT-1) = 2./3. / (ZEPS + ZBNEG2(:,:,IT-1))**2 ! Second negative non-normalized weight
! 
PR(:,:,IT-1) = (ZOMN2(:,:,IT-1)/(ZOMN1(:,:,IT-1)+ZOMN2(:,:,IT-1))*ZFNEG2(:,:,IT-1) + &
     (ZOMN1(:,:,IT-1)/(ZOMN1(:,:,IT-1)+ZOMN2(:,:,IT-1))*ZFNEG1(:,:,IT-1))) &
     *(0.5-SIGN(0.5,PRWCT(:,:,IT-1))) + &
             (ZOMP2(:,:,IT-1)/(ZOMP1(:,:,IT-1)+ZOMP2(:,:,IT-1))*ZFPOS2(:,:,IT-1) + &
               (ZOMP1(:,:,IT-1)/(ZOMP1(:,:,IT-1)+ZOMP2(:,:,IT-1))*ZFPOS1(:,:,IT-1))) &
     * (0.5+SIGN(0.5,PRWCT(:,:,IT-1)))  ! Total flux
!
PR = PR * PRWCT ! Add contravariant flux
!
END FUNCTION WENO_K_3_WZ
!
!-----------------------------------------------------------------------------
!
!     ########################################################################
      FUNCTION WENO_K_3_MZ(PSRC, PRWCT) RESULT(PR)
!     ########################################################################
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
!! T. Lunet 02/10/2014:  Complete code documentation
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
USE MODD_CONF
USE MODD_PARAMETERS,ONLY: JPVEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX  ! X direction LBC type
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on W grid
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
!
! intermediate reconstruction fluxes for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFPOS1, ZFPOS2, ZFPOS3
!
! intermediate reconstruction fluxes for negative wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZFNEG1, ZFNEG2, ZFNEG3
!
! smoothness indicators for positive wind case
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBPOS1, ZBPOS2, ZBPOS3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZBNEG1, ZBNEG2, ZBNEG3
!
! WENO weights
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZOMP1, ZOMP2, ZOMP3
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)):: ZOMN1, ZOMN2, ZOMN3
!
! EPSILON for weno weights calculation
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
ZFPOS3 = 0.0
ZFNEG1 = 0.0
ZFNEG2 = 0.0
ZFNEG3 = 0.0
ZBPOS1 = 0.0
ZBPOS2 = 0.0
ZBPOS3 = 0.0
ZBNEG1 = 0.0
ZBNEG2 = 0.0
ZBNEG3 = 0.0
ZOMP1  = 0.0
ZOMP2  = 0.0
ZOMP3  = 0.0
ZOMN1  = 0.0
ZOMN2  = 0.0
ZOMN3  = 0.0 
!
!-------------------------------------------------------------------------------
!*       1.1.   Interior Fluxes 
!              ---------------------
!
!-------------------------------------------------------------------------------
! Flux calculation in the physical domain far enough from the boundary 
! WENO scheme order 5, IB+2 -> IT-1
! Computation at the mass point on Wgrid v(i-1/2,j,k)
!-------------------------------------------------------------------------------
!
! ----- Positive fluxes -----
!
! First positive stencil, needs indices i-3, i-2, i-1
ZFPOS1(:,:,IB+2:IT-1) = 1./6.   * (2.0*PSRC(:,:,IB-1:IT-4) - 7.0*PSRC(:,:,IB:IT-3) + 11.0*PSRC(:,:,IB+1:IT-2))   ! Flux
ZBPOS1(:,:,IB+2:IT-1) = 13./12. * (    PSRC(:,:,IB-1:IT-4) - 2.0*PSRC(:,:,IB:IT-3) +      PSRC(:,:,IB+1:IT-2))**2 & 
      + 1./4    * (    PSRC(:,:,IB-1:IT-4) - 4.0*PSRC(:,:,IB:IT-3) + 3.0* PSRC(:,:,IB+1:IT-2))**2  ! Smoothness indicator
ZOMP1(:,:,IB+2:IT-1)  = 1./10. /  (ZEPS + ZBPOS1(:,:,IB+2:IT-1))**2             ! Non-normalized weight
!
! Second positive stencil, needs indices i-2, i-1, i
ZFPOS2(:,:,IB+2:IT-1) = 1./6.  * (-1.0*PSRC(:,:,IB:IT-3) + 5.0*PSRC(:,:,IB+1:IT-2) + 2.0*PSRC(:,:,IB+2:IT-1))  ! Flux
ZBPOS2(:,:,IB+2:IT-1) = 13./12 * (     PSRC(:,:,IB:IT-3) - 2.0*PSRC(:,:,IB+1:IT-2) +     PSRC(:,:,IB+2:IT-1))**2 &
      + 1./4   * (     PSRC(:,:,IB:IT-3) -                               PSRC(:,:,IB+2:IT-1))**2 ! Smoothness indicator
ZOMP2(:,:,IB+2:IT-1)  = 3./5. /  (ZEPS + ZBPOS2(:,:,IB+2:IT-1))**2             ! Non-normalized weight
!
! Third positive stencil, needs indices i-1, i, i+1
ZFPOS3(:,:,IB+2:IT-1) = 1./6   * (2.0*PSRC(:,:,IB+1:IT-2) + 5.0*PSRC(:,:,IB+2:IT-1) - PSRC(:,:,IB+3:IT))  ! Flux
ZBPOS3(:,:,IB+2:IT-1) = 13./12 * ( PSRC(:,:,IB+1:IT-2) - 2.0*PSRC(:,:,IB+2:IT-1) + PSRC(:,:,IB+3:IT))**2 &
      + 1./4   * (3.0*PSRC(:,:,IB+1:IT-2) - 4.0*PSRC(:,:,IB+2:IT-1) + PSRC(:,:,IB+3:IT))**2 ! Smoothness indicator
ZOMP3(:,:,IB+2:IT-1)  = 3./10. / (ZEPS + ZBPOS3(:,:,IB+2:IT-1))**2           ! Non-normalized weight
!
! ----- Negative fluxes ----- 
!
! First negative stencil, needs indices i, i+1, i+2
ZFNEG1(:,:,IB+2:IT-1) = 1./6.   * (11.0*PSRC(:,:,IB+2:IT-1) - 7.0*PSRC(:,:,IB+3:IT) + 2.0*PSRC(:,:,IB+4:IT+1))   ! Flux
ZBNEG1(:,:,IB+2:IT-1) = 13./12. * (     PSRC(:,:,IB+2:IT-1) - 2.0*PSRC(:,:,IB+3:IT) +     PSRC(:,:,IB+4:IT+1))**2 & 
      + 1./4    * (3.0* PSRC(:,:,IB+2:IT-1) - 4.0*PSRC(:,:,IB+3:IT) +     PSRC(:,:,IB+4:IT+1))**2  ! Smoothness indicator
ZOMN1(:,:,IB+2:IT-1)  = 1./10. /  (ZEPS + ZBNEG1(:,:,IB+2:IT-1))**2             ! Non-normalized weight
!
! Second negative stencil, needs indices i-1, i, i+1
ZFNEG2(:,:,IB+2:IT-1) = 1./6.  * (2.0*PSRC(:,:,IB+1:IT-2) + 5.0*PSRC(:,:,IB+2:IT-1) - 1.0*PSRC(:,:,IB+3:IT))  ! Flux
ZBNEG2(:,:,IB+2:IT-1) = 13./12 * (    PSRC(:,:,IB+1:IT-2) - 2.0*PSRC(:,:,IB+2:IT-1) +     PSRC(:,:,IB+3:IT))**2 &
      + 1./4   * (    PSRC(:,:,IB+1:IT-2) -                               PSRC(:,:,IB+3:IT))**2 ! Smoothness indicator
ZOMN2(:,:,IB+2:IT-1)  = 3./5. /  (ZEPS + ZBNEG2(:,:,IB+2:IT-1))**2            ! Non-normalized weight
!
! Third negative stencil, needs indices i-2, i-1, i
ZFNEG3(:,:,IB+2:IT-1) = 1./6   * (-1.0*PSRC(:,:,IB:IT-3) + 5.0*PSRC(:,:,IB+1:IT-2) + 2.0*PSRC(:,:,IB+2:IT-1))  ! Flux
ZBNEG3(:,:,IB+2:IT-1) = 13./12 * (  PSRC(:,:,IB:IT-3) - 2.0*PSRC(:,:,IB+1:IT-2) +     PSRC(:,:,IB+2:IT-1))**2 &
      + 1./4   * (     PSRC(:,:,IB:IT-3) - 4.0*PSRC(:,:,IB+1:IT-2) + 3.0*PSRC(:,:,IB+2:IT-1))**2 ! Smoothness indicator
ZOMN3(:,:,IB+2:IT-1)  = 3./10. / (ZEPS + ZBNEG3(:,:,IB+2:IT-1))**2             ! Non-normalized weight
!
!
! ----- Total flux -----
!
PR(:,:,IB+2:IT-1) = (ZOMP1(:,:,IB+2:IT-1)/(ZOMP1(:,:,IB+2:IT-1)+ZOMP2(:,:,IB+2:IT-1)+ZOMP3(:,:,IB+2:IT-1)) &
           * ZFPOS1(:,:,IB+2:IT-1) + &
                      ZOMP2(:,:,IB+2:IT-1)/(ZOMP1(:,:,IB+2:IT-1)+ZOMP2(:,:,IB+2:IT-1)+ZOMP3(:,:,IB+2:IT-1)) &
           * ZFPOS2(:,:,IB+2:IT-1) + & 
                      ZOMP3(:,:,IB+2:IT-1)/(ZOMP1(:,:,IB+2:IT-1)+ZOMP2(:,:,IB+2:IT-1)+ZOMP3(:,:,IB+2:IT-1)) &
           * ZFPOS3(:,:,IB+2:IT-1)) &
                    * (0.5+SIGN(0.5,PRWCT(:,:,IB+2:IT-1))) &
                  + (ZOMN1(:,:,IB+2:IT-1)/(ZOMN1(:,:,IB+2:IT-1)+ZOMN2(:,:,IB+2:IT-1)+ZOMN3(:,:,IB+2:IT-1)) &
           * ZFNEG1(:,:,IB+2:IT-1)  &
                     + ZOMN2(:,:,IB+2:IT-1)/(ZOMN1(:,:,IB+2:IT-1)+ZOMN2(:,:,IB+2:IT-1)+ZOMN3(:,:,IB+2:IT-1)) &
           * ZFNEG2(:,:,IB+2:IT-1)  &
                     + ZOMN3(:,:,IB+2:IT-1)/(ZOMN1(:,:,IB+2:IT-1)+ZOMN2(:,:,IB+2:IT-1)+ZOMN3(:,:,IB+2:IT-1)) &
           * ZFNEG3(:,:,IB+2:IT-1))  & 
                    * (0.5-SIGN(0.5,PRWCT(:,:,IB+2:IT-1)))
!
!-------------------------------------------------------------------------------
!*       1.2.   Bottom border
!               ---------------------
!---------------------------------------------------------------------------
! WENO order reduction
!---------------------------------------------------------------------------
!
! WENO scheme order 1, IB
PR(:,:,IB) = PSRC(:,:,IB-1) * (0.5+SIGN(0.5,PRWCT(:,:,IB) )) &
           + PSRC(:,:,IB  ) * (0.5-SIGN(0.5,PRWCT(:,:,IB) ))
!
! WENO scheme order 3, IB+1
ZFPOS1(:,:,IB+1) = 0.5 * (3.0*PSRC(:,:,IB) - PSRC(:,:,IB-1)) ! First positive flux
ZFPOS2(:,:,IB+1) = 0.5 * ( PSRC(:,:,IB) + PSRC(:,:,IB+1)) ! Second positive flux
ZBPOS1(:,:,IB+1) = (PSRC(:,:,IB)   - PSRC(:,:,IB-1))**2 ! First positive smoothness indicator
ZBPOS2(:,:,IB+1) = (PSRC(:,:,IB+1) - PSRC(:,:,IB))**2  ! Second positive smoothness indicator
!
ZFNEG1(:,:,IB+1) = 0.5 * (3.0*PSRC(:,:,IB+1) - PSRC(:,:,IB+2)) ! First negative flux
ZFNEG2(:,:,IB+1) = 0.5 * ( PSRC(:,:,IB)   + PSRC(:,:,IB+1)) ! Second negative flux
ZBNEG1(:,:,IB+1) = (PSRC(:,:,IB+1) - PSRC(:,:,IB+2))**2 ! First negative smoothness indicator
ZBNEG2(:,:,IB+1) = (PSRC(:,:,IB)   - PSRC(:,:,IB+1))**2 ! Second negative smoothness indicator
!
ZOMP1(:,:,IB+1) = 1./3. / (ZEPS + ZBPOS1(:,:,IB+1))**2 ! First positive non-normalized weight
ZOMP2(:,:,IB+1) = 2./3. / (ZEPS + ZBPOS2(:,:,IB+1))**2 ! Second positive non-normalized weight
ZOMN1(:,:,IB+1) = 1./3. / (ZEPS + ZBNEG1(:,:,IB+1))**2 ! First negative non-normalized weight
ZOMN2(:,:,IB+1) = 2./3. / (ZEPS + ZBNEG2(:,:,IB+1))**2 ! Second negative non-normalized weight
! 
PR(:,:,IB+1) = (ZOMP2(:,:,IB+1)/(ZOMP1(:,:,IB+1)+ZOMP2(:,:,IB+1)) * ZFPOS2(:,:,IB+1)  &
              + ZOMP1(:,:,IB+1)/(ZOMP1(:,:,IB+1)+ZOMP2(:,:,IB+1)) * ZFPOS1(:,:,IB+1)) * &
    (0.5+SIGN(0.5,PRWCT(:,:,IB+1) ))           &
             + (ZOMN2(:,:,IB+1)/(ZOMN1(:,:,IB+1)+ZOMN2(:,:,IB+1)) * ZFNEG2(:,:,IB+1)    &
              + ZOMN1(:,:,IB+1)/(ZOMN1(:,:,IB+1)+ZOMN2(:,:,IB+1)) * ZFNEG1(:,:,IB+1)) * &
    (0.5-SIGN(0.5,PRWCT(:,:,IB+1) )) ! Total flux
! 
!-------------------------------------------------------------------------------
!*       1.3.   Top border
!               ---------------------
!
!---------------------------------------------------------------------------
! Open, or Wall, or Nest boundary condition => WENO order reduction
!---------------------------------------------------------------------------
!
! WENO scheme order 1, IT+1
PR(:,:,IT+1) = PSRC(:,:,IT  ) * (0.5+SIGN(0.5,PRWCT(:,:,IT+1) )) &
             + PSRC(:,:,IT+1) * (0.5-SIGN(0.5,PRWCT(:,:,IT+1) ))
!
! WENO scheme order 3, IT
ZFPOS1(:,:,IT) = 0.5 * (3.0*PSRC(:,:,IT-1) - PSRC(:,:,IT-2)) ! First positive flux
ZFPOS2(:,:,IT) = 0.5 * ( PSRC(:,:,IT-1) + PSRC(:,:,IT)) ! Second positive flux
ZBPOS1(:,:,IT) = (PSRC(:,:,IT-1) - PSRC(:,:,IT-2))**2 ! First positive smoothness indicator
ZBPOS2(:,:,IT) = (PSRC(:,:,IT)   - PSRC(:,:,IT-1))**2 ! Second positive smoothness indicator
!
ZFNEG1(:,:,IT) = 0.5 * (3.0*PSRC(:,:,IT)   - PSRC(:,:,IT+1)) ! First negative flux
ZFNEG2(:,:,IT) = 0.5 * ( PSRC(:,:,IT-1) + PSRC(:,:,IT)) ! Second negative flux
ZBNEG1(:,:,IT) = (PSRC(:,:,IT)   - PSRC(:,:,IT+1))**2 ! First negative smoothness indicator
ZBNEG2(:,:,IT) = (PSRC(:,:,IT-1) - PSRC(:,:,IT))**2  ! Second negative smoothness indicator
!
ZOMP1(:,:,IT) = 1./3. / (ZEPS + ZBPOS1(:,:,IT))**2 ! First positive non-normalized weight
ZOMP2(:,:,IT) = 2./3. / (ZEPS + ZBPOS2(:,:,IT))**2 ! Second positive non-normalized weight
ZOMN1(:,:,IT) = 1./3. / (ZEPS + ZBNEG1(:,:,IT))**2 ! First negative non-normalized weight
ZOMN2(:,:,IT) = 2./3. / (ZEPS + ZBNEG2(:,:,IT))**2 ! Second negative non-normalized weight
! 
PR(:,:,IT) = (ZOMP2(:,:,IT)/(ZOMP1(:,:,IT)+ZOMP2(:,:,IT)) * ZFPOS2(:,:,IT)     &
            + ZOMP1(:,:,IT)/(ZOMP1(:,:,IT)+ZOMP2(:,:,IT)) * ZFPOS1(:,:,IT)) *  &
            (0.5+SIGN(0.5,PRWCT(:,:,IT) ))                                     &
           + (ZOMN2(:,:,IT)/(ZOMN1(:,:,IT)+ZOMN2(:,:,IT)) * ZFNEG2(:,:,IT)     &
            + ZOMN1(:,:,IT)/(ZOMN1(:,:,IT)+ZOMN2(:,:,IT)) * ZFNEG1(:,:,IT)) *  &
            (0.5-SIGN(0.5,PRWCT(:,:,IT) ))  ! Total flux
!
PR = PR * PRWCT ! Add contravariant flux
!
END FUNCTION WENO_K_3_MZ
