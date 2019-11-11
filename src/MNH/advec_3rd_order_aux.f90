!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ###############################
      MODULE MODI_ADVEC_3RD_ORDER_AUX
!     ###############################
!!    AUTHOR
!!    ------
!!
!! Correction :	
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
INTERFACE
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ADVEC_3RD_ORDER_UX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_3RD_ORDER_UX
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ADVEC_3RD_ORDER_MX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_3RD_ORDER_MX
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ADVEC_3RD_ORDER_VY(HLBCY,PSRC, PRVCT, PR, TPHALO2)
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_3RD_ORDER_VY
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE ADVEC_3RD_ORDER_MY(HLBCY,PSRC, PRVCT, PR, TPHALO2)
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
END SUBROUTINE ADVEC_3RD_ORDER_MY
!
!------------------------------------------------------------------------
!
      FUNCTION UP3_WZ(PSRC, PRWCT) RESULT(PR)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on W grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
! output source term
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
END FUNCTION UP3_WZ
!
!-------------------------------------------------------------------------------
!
      FUNCTION UP3_MZ(PSRC, PRWCT) RESULT(PR)
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on W grid
!
! output source term
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
END FUNCTION UP3_MZ
!
END INTERFACE
!
END MODULE MODI_ADVEC_3RD_ORDER_AUX
!
!-------------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE ADVEC_3RD_ORDER_UX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
!     #############################################################
!!
!!****  ADVEC_3RD_ORDER_UX - 3rd order upstream fluxes of U in X direction
!!              input variable PSRC is on U grid, and output PR is on mass grid
!!
!!    AUTHOR
!!    ------
!!      C.Lac            * CNRM/GMME *               
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
USE MODE_ll
!
USE MODD_LUNIT
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
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
!-------------------------------------------------------------------------------
!
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
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE X DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
  IW=IIB+1
  IE=IIE
!
  IWF=IW-1
  IEF=IE-1
!
  PR(IWF:IEF,:,:) = 1./6. * ( (2.*PSRC(IW:IE,:,:) + 5.*PSRC(IW-1:IE-1,:,:) -   &
                   PSRC(IW-2:IE-2,:,:)) * (0.5+SIGN(0.5,PRUCT(IW-1:IE-1,:,:))) &
                            + (5.*PSRC(IW:IE,:,:) + 2.*PSRC(IW-1:IE-1,:,:) -   &
                   PSRC(IW+1:IE+1,:,:)) * (0.5-SIGN(0.5,PRUCT(IW-1:IE-1,:,:))) )
!
  PR(IEF+1,:,:) = 1./6. * ( (2.*PSRC(IE+1,:,:) + 5.*PSRC(IE,:,:) -             &
                               PSRC(IE-1,:,:)) * (0.5+SIGN(0.5,PRUCT(IE,:,:))) &
                          + (5.*PSRC(IE+1,:,:) + 2.*PSRC(IE,:,:) -             &
                            TPHALO2%EAST(:,:)) * (0.5-SIGN(0.5,PRUCT(IE,:,:))))
!
  PR(IWF-1,:,:) = 1./6. * ( (2.*PSRC(IW-1,:,:) + 5.*PSRC(IW-2,:,:) -           &
                          TPHALO2%WEST(:,:)) * (0.5+SIGN(0.5,PRUCT(IW-2,:,:))) &
                          + (5.*PSRC(IW-1,:,:) + 2.*PSRC(IW-2,:,:) -           &
                               PSRC(IW,:,:)) * (0.5-SIGN(0.5,PRUCT(IW-2,:,:))) )
!
!       OPEN, WALL, NEST CASE IN THE X DIRECTION
!
CASE ('OPEN','WALL','NEST')
!
!       USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSICAL BORDER
!
  IF (LWEST_ll()) THEN
      IW=IIB+2          ! special case of C grid
  ELSE
!!$    IF(NHALO == 1) THEN
      IW=IIB+1
!!$    ELSE
!!$      IW=IIB
!!$    ENDIF
  ENDIF
!!$  IF (LEAST_ll() .OR. NHALO == 1) THEN
  IF (LEAST_ll()) THEN
    IE=IIE
  ELSE
    IE=IIE
  END IF
!
  IWF=IW-1
  IEF=IE-1
!
  IF(LWEST_ll()) THEN
    PR(IWF-1,:,:) = PSRC(IW-2,:,:) * (0.5+SIGN(0.5,PRUCT(IW-2,:,:))) &
                  + PSRC(IW-1,:,:) * (0.5-SIGN(0.5,PRUCT(IW-2,:,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(IWF-1,:,:) = 1./6. * ( (2.*PSRC(IW-1,:,:) + 5.*PSRC(IW-2,:,:) -         &
                          TPHALO2%WEST(:,:)) * (0.5+SIGN(0.5,PRUCT(IW-2,:,:))) &
                          + (5.*PSRC(IW-1,:,:) + 2.*PSRC(IW-2,:,:) -           &
                               PSRC(IW,:,:)) * (0.5-SIGN(0.5,PRUCT(IW-2,:,:))) )
  ENDIF
!
  IF(LEAST_ll()) THEN
    PR(IEF+1,:,:) = PSRC(IE,:,:)   * (0.5+SIGN(0.5,PRUCT(IE,:,:))) &
                  + PSRC(IE+1,:,:) * (0.5-SIGN(0.5,PRUCT(IE,:,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(IEF+1,:,:) = 1./6. * ( (2.*PSRC(IE+1,:,:) + 5.*PSRC(IE,:,:) -           &
                               PSRC(IE-1,:,:)) * (0.5+SIGN(0.5,PRUCT(IE,:,:))) &
                          + (5.*PSRC(IE+1,:,:) + 2.*PSRC(IE,:,:) -             &
                            TPHALO2%EAST(:,:)) * (0.5-SIGN(0.5,PRUCT(IE,:,:))))
  ENDIF
!
!      USE A THIRD ORDER UPSTREAM SCHEME ELSEWHERE 
!
  PR(IWF:IEF,:,:) = 1./6. * ( (2.*PSRC(IW:IE,:,:) + 5.*PSRC(IW-1:IE-1,:,:) -   &
                   PSRC(IW-2:IE-2,:,:)) * (0.5+SIGN(0.5,PRUCT(IW-1:IE-1,:,:))) &
                            + (5.*PSRC(IW:IE,:,:) + 2.*PSRC(IW-1:IE-1,:,:) -   &
                   PSRC(IW+1:IE+1,:,:)) * (0.5-SIGN(0.5,PRUCT(IW-1:IE-1,:,:))) )
!
END SELECT
!
PR = PR * PRUCT
!
END SUBROUTINE ADVEC_3RD_ORDER_UX 
!
!-------------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE ADVEC_3RD_ORDER_MX(HLBCX,PSRC, PRUCT, PR, TPHALO2)
!     #############################################################
!!
!!**** ADVEC_3RD_ORDER_MX - 3rd order upstream fluxes of variable in X direction
!!     Input variable PSRC is on MASS grid, and output PR is on U grid
!!
!!    AUTHOR
!!    ------
!!      C.Lac            * CNRM/GMME *               
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER:: IW,IE,IWF,IEF   ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!-------------------------------------------------------------------------------
!
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
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCX(1) ) ! X direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE X DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
  IW=IIB+1
  IE=IIE
!
  IWF=IW
  IEF=IE
!
  PR(IWF:IEF,:,:) = 1./6. * ( (2.*PSRC(IW:IE,:,:) + 5.*PSRC(IW-1:IE-1,:,:) -   &
                       PSRC(IW-2:IE-2,:,:)) * (0.5+SIGN(0.5,PRUCT(IW:IE,:,:))) &
                            + (5.*PSRC(IW:IE,:,:) + 2.*PSRC(IW-1:IE-1,:,:) -   &
                       PSRC(IW+1:IE+1,:,:)) * (0.5-SIGN(0.5,PRUCT(IW:IE,:,:))) )
!
  PR(IWF-1,:,:) = 1./6. * ( (2.*PSRC(IW-1,:,:) + 5.*PSRC(IW-2,:,:) -           &
                          TPHALO2%WEST(:,:)) * (0.5+SIGN(0.5,PRUCT(IW-1,:,:))) &
                          + (5.*PSRC(IW-1,:,:) + 2.*PSRC(IW-2,:,:) -           &
                               PSRC(IW,:,:)) * (0.5-SIGN(0.5,PRUCT(IW-1,:,:))) )
!
  PR(IEF+1,:,:) = 1./6. * ( (2.*PSRC(IE+1,:,:) + 5.*PSRC(IE,:,:) -             &
                             PSRC(IE-1,:,:)) * (0.5+SIGN(0.5,PRUCT(IE+1,:,:))) &
                          + (5.*PSRC(IE+1,:,:) + 2.*PSRC(IE,:,:) -             &
                               TPHALO2%EAST) * (0.5-SIGN(0.5,PRUCT(IE+1,:,:))) )
!
!    OPEN, WALL, NEST CASE IN THE X DIRECTION
!
CASE ('OPEN','WALL','NEST')
!
!    USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSCIAL BORDER 
!
  IF (LWEST_ll()) THEN
    IW=IIB+1
  ELSE
!!$    IF(NHALO == 1) THEN
      IW=IIB+1
!!$    ELSE
!!$      IW=IIB
!!$    ENDIF
  ENDIF
!!$  IF (LEAST_ll() .OR. NHALO == 1) THEN
  IF (LEAST_ll()) THEN
    IE=IIE
  ELSE
    IE=IIE
  END IF  
!
  IWF=IW
  IEF=IE
!
  IF(LWEST_ll()) THEN
    PR(IWF-1,:,:) = PSRC(IW-2,:,:) * (0.5+SIGN(0.5,PRUCT(IW-1,:,:))) &
                  + PSRC(IW-1,:,:) * (0.5-SIGN(0.5,PRUCT(IW-1,:,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(IWF-1,:,:) = 1./6. * ( (2.*PSRC(IW-1,:,:) + 5.*PSRC(IW-2,:,:) -         &
                          TPHALO2%WEST(:,:)) * (0.5+SIGN(0.5,PRUCT(IW-1,:,:))) &
                          + (5.*PSRC(IW-1,:,:) + 2.*PSRC(IW-2,:,:) -           &
                               PSRC(IW,:,:)) * (0.5-SIGN(0.5,PRUCT(IW-1,:,:))) )
  ENDIF
!
  IF(LEAST_ll()) THEN
    PR(IEF+1,:,:) = PSRC(IE,:,:)   * (0.5+SIGN(0.5,PRUCT(IE+1,:,:))) &
                  + PSRC(IE+1,:,:) * (0.5-SIGN(0.5,PRUCT(IE+1,:,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(IEF+1,:,:) = 1./6. * ( (2.*PSRC(IE+1,:,:) + 5.*PSRC(IE,:,:) -           &
                             PSRC(IE-1,:,:)) * (0.5+SIGN(0.5,PRUCT(IE+1,:,:))) &
                          + (5.*PSRC(IE+1,:,:) + 2.*PSRC(IE,:,:) -             &
                               TPHALO2%EAST) * (0.5-SIGN(0.5,PRUCT(IE+1,:,:))) )
  ENDIF
!
!    USE A THIRD ORDER UPSTREAM SCHEME ELSEWHERE
!
  PR(IWF:IEF,:,:) = 1./6. * ( (2.*PSRC(IW:IE,:,:) + 5.*PSRC(IW-1:IE-1,:,:) -   &
                       PSRC(IW-2:IE-2,:,:)) * (0.5+SIGN(0.5,PRUCT(IW:IE,:,:))) &
                            + (5.*PSRC(IW:IE,:,:) + 2.*PSRC(IW-1:IE-1,:,:) -   &
                       PSRC(IW+1:IE+1,:,:)) * (0.5-SIGN(0.5,PRUCT(IW:IE,:,:))) )
!
END SELECT
!
PR = PR * PRUCT
!
END SUBROUTINE ADVEC_3RD_ORDER_MX
!
!-------------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE ADVEC_3RD_ORDER_VY(HLBCY,PSRC, PRVCT, PR, TPHALO2)
!     #############################################################
!!
!!****  ADVEC_3RD_ORDER_VY - 3rd order upstream fluxes of V in Y direction
!!      Input variable PSRC is on V grid, and output PR is on MASS grid
!!
!!    AUTHOR
!!    ------
!!      C.Lac            * CNRM/GMME *               
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER:: IS,IN,ISF,INF   ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!-------------------------------------------------------------------------------
!
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
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCY(1) ) ! 
!
!*       1.1    CYCLIC CASE IN THE Y DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCY(1) == HLBCY(2)
!
  IS=IJB+1
  IN=IJE
!
  ISF=IS-1
  INF=IN-1
!
  PR(:,ISF:INF,:) = 1./6. * ( (2.*PSRC(:,IS:IN,:) + 5.*PSRC(:,IS-1:IN-1,:) -   &
                   PSRC(:,IS-2:IN-2,:)) * (0.5+SIGN(0.5,PRVCT(:,IS-1:IN-1,:))) &
                            + (5.*PSRC(:,IS:IN,:) + 2.*PSRC(:,IS-1:IN-1,:) -   &
                   PSRC(:,IS+1:IN+1,:)) * (0.5-SIGN(0.5,PRVCT(:,IS-1:IN-1,:))) )
!
  PR(:,ISF-1,:) = 1./6. * ( (2.*PSRC(:,IS-1,:) + 5.*PSRC(:,IS-2,:) -           &
                         TPHALO2%SOUTH(:,:)) * (0.5+SIGN(0.5,PRVCT(:,IS-2,:))) &
                          + (5.*PSRC(:,IS-1,:) + 2.*PSRC(:,IS-2,:) -           &
                               PSRC(:,IS,:)) * (0.5-SIGN(0.5,PRVCT(:,IS-2,:))) )
!
  PR(:,INF+1,:) = 1./6. * ( (2.*PSRC(:,IN+1,:) + 5.*PSRC(:,IN,:) -             &
                               PSRC(:,IN-1,:)) * (0.5+SIGN(0.5,PRVCT(:,IN,:))) &
                          + (5.*PSRC(:,IN+1,:) + 2.*PSRC(:,IN,:) -             &
                           TPHALO2%NORTH(:,:)) * (0.5-SIGN(0.5,PRVCT(:,IN,:))) )
!
!       OPEN, WALL, NEST CASES IN THE Y DIRECTION 
!
CASE ('OPEN','WALL','NEST')
!
!       USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSICAL BORDER 
!
  IF (LSOUTH_ll()) THEN
    IS=IJB+2
  ELSE
!!$    IF(NHALO == 1) THEN
      IS=IJB+1
!!$    ELSE
!!$      IS=IJB
!!$    ENDIF
  ENDIF
!!$  IF (LNORTH_ll() .OR. NHALO == 1) THEN
  IF (LNORTH_ll()) THEN
    IN=IJE
  ELSE
    IN=IJE
  END IF
!
  ISF=IS-1
  INF=IN-1
!
  IF(LSOUTH_ll()) THEN
    PR(:,ISF-1,:) = PSRC(:,IS-2,:) * (0.5+SIGN(0.5,PRVCT(:,IS-2,:))) &
                  + PSRC(:,IS-1,:) * (0.5-SIGN(0.5,PRVCT(:,IS-2,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(:,ISF-1,:) = 1./6. * ( (2.*PSRC(:,IS-1,:) + 5.*PSRC(:,IS-2,:) -         &
                         TPHALO2%SOUTH(:,:)) * (0.5+SIGN(0.5,PRVCT(:,IS-2,:))) &
                          + (5.*PSRC(:,IS-1,:) + 2.*PSRC(:,IS-2,:) -           &
                               PSRC(:,IS,:)) * (0.5-SIGN(0.5,PRVCT(:,IS-2,:))) )
  ENDIF
!
  IF(LNORTH_ll()) THEN
    PR(:,INF+1,:) = PSRC(:,IN,:)   * (0.5+SIGN(0.5,PRVCT(:,IN,:))) &
                  + PSRC(:,IN+1,:) * (0.5-SIGN(0.5,PRVCT(:,IN,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(:,INF+1,:) = 1./6. * ( (2.*PSRC(:,IN+1,:) + 5.*PSRC(:,IN,:) -           &
                               PSRC(:,IN-1,:)) * (0.5+SIGN(0.5,PRVCT(:,IN,:))) &
                          + (5.*PSRC(:,IN+1,:) + 2.*PSRC(:,IN,:) -             &
                           TPHALO2%NORTH(:,:)) * (0.5-SIGN(0.5,PRVCT(:,IN,:))) )
  ENDIF
!
!       USE A 3RD ORDER UPSTREAM SCHEME ELSEWHERE
!
  PR(:,ISF:INF,:) = 1./6. * ( (2.*PSRC(:,IS:IN,:) + 5.*PSRC(:,IS-1:IN-1,:) -   &
                   PSRC(:,IS-2:IN-2,:)) * (0.5+SIGN(0.5,PRVCT(:,IS-1:IN-1,:))) &
                            + (5.*PSRC(:,IS:IN,:) + 2.*PSRC(:,IS-1:IN-1,:) -   &
                   PSRC(:,IS+1:IN+1,:)) * (0.5-SIGN(0.5,PRVCT(:,IS-1:IN-1,:))) )
!
END SELECT
!
PR = PR * PRVCT
!
END SUBROUTINE ADVEC_3RD_ORDER_VY
!
!-------------------------------------------------------------------------------
!
!     ##############################################################
      SUBROUTINE ADVEC_3RD_ORDER_MY(HLBCY, PSRC, PRVCT, PR, TPHALO2)
!     ##############################################################
!!
!!**** ADVEC_3RD_ORDER_MY - 3rd order upstream fluxes of variable in Y direction
!!     Input variable PSRC is on MASS grid, and output PR is on V grid
!!
!!    AUTHOR
!!    ------
!!      C.Lac            * CNRM/GMME *               
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
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PR
TYPE(HALO2_ll), OPTIONAL, POINTER :: TPHALO2      ! halo2 for the field at t
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IJB    ! Begining useful area in x,y,z directions
INTEGER :: IIE,IJE    ! End useful area in x,y,z directions
INTEGER:: IS,IN,ISF,INF   ! Coordinate of third order diffusion area
!
INTEGER:: ILUOUT,IRESP   ! for prints
!
!-------------------------------------------------------------------------------
!
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
!-------------------------------------------------------------------------------
!
SELECT CASE ( HLBCY(1) ) ! Y direction LBC type: (1) for left side
!
!*       1.1    CYCLIC CASE IN THE Y DIRECTION:
!
CASE ('CYCL')          ! In that case one must have HLBCX(1) == HLBCX(2)
!
  IS=IJB+1
  IN=IJE
!
  ISF=IS
  INF=IN
!
  PR(:,ISF:INF,:) = 1./6. * ( (2.*PSRC(:,IS:IN,:) + 5.*PSRC(:,IS-1:IN-1,:) -   &
                       PSRC(:,IS-2:IN-2,:)) * (0.5+SIGN(0.5,PRVCT(:,IS:IN,:))) &
                            + (5.*PSRC(:,IS:IN,:) + 2.*PSRC(:,IS-1:IN-1,:) -   &
                       PSRC(:,IS+1:IN+1,:)) * (0.5-SIGN(0.5,PRVCT(:,IS:IN,:))) )
!
  PR(:,ISF-1,:) = 1./6. * ( (2.*PSRC(:,IS-1,:) + 5.*PSRC(:,IS-2,:) -           &
                         TPHALO2%SOUTH(:,:)) * (0.5+SIGN(0.5,PRVCT(:,IS-1,:))) &
                          + (5.*PSRC(:,IS-1,:) + 2.*PSRC(:,IS-2,:) -           &
                               PSRC(:,IS,:)) * (0.5-SIGN(0.5,PRVCT(:,IS-1,:))) )
!
  PR(:,INF+1,:) = 1./6. * ( (2.*PSRC(:,IN+1,:) + 5.*PSRC(:,IN,:) -             &
                             PSRC(:,IN-1,:)) * (0.5+SIGN(0.5,PRVCT(:,IN+1,:))) &
                          + (5.*PSRC(:,IN+1,:) + 2.*PSRC(:,IN,:) -             &
                         TPHALO2%NORTH(:,:)) * (0.5-SIGN(0.5,PRVCT(:,IN+1,:))) )
!
!       OPEN, WALL, NEST CASES IN THE Y DIRECTION 
!
CASE ('OPEN','WALL','NEST')
!
!       USE A FIRST ORDER UPSTREAM SCHEME AT THE PHYSICAL BORDER 
!
  IF (LSOUTH_ll()) THEN
    IS=IJB+1
  ELSE
!!$    IF(NHALO == 1) THEN
      IS=IJB+1
!!$    ELSE
!!$      IS=IJB
!!$    ENDIF
  ENDIF
!!$  IF (LNORTH_ll() .OR. NHALO == 1) THEN
  IF (LNORTH_ll()) THEN
    IN=IJE
  ELSE
    IN=IJE
  END IF
!
  ISF=IS
  INF=IN
!
  IF(LSOUTH_ll()) THEN
    PR(:,ISF-1,:) = PSRC(:,IS-2,:) * (0.5+SIGN(0.5,PRVCT(:,IS-1,:))) &
                  + PSRC(:,IS-1,:) * (0.5-SIGN(0.5,PRVCT(:,IS-1,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(:,ISF-1,:) = 1./6. * ( (2.*PSRC(:,IS-1,:) + 5.*PSRC(:,IS-2,:) -         &
                         TPHALO2%SOUTH(:,:)) * (0.5+SIGN(0.5,PRVCT(:,IS-1,:))) &
                          + (5.*PSRC(:,IS-1,:) + 2.*PSRC(:,IS-2,:) -           &
                               PSRC(:,IS,:)) * (0.5-SIGN(0.5,PRVCT(:,IS-1,:))) )
  END IF
!
  IF(LNORTH_ll()) THEN
    PR(:,INF+1,:) = PSRC(:,IN,:)   * (0.5+SIGN(0.5,PRVCT(:,IN+1,:))) &
                  + PSRC(:,IN+1,:) * (0.5-SIGN(0.5,PRVCT(:,IN+1,:)))
!!$  ELSEIF (NHALO == 1) THEN
  ELSE
    PR(:,INF+1,:) = 1./6. * ( (2.*PSRC(:,IN+1,:) + 5.*PSRC(:,IN,:) -           &
                             PSRC(:,IN-1,:)) * (0.5+SIGN(0.5,PRVCT(:,IN+1,:))) &
                          + (5.*PSRC(:,IN+1,:) + 2.*PSRC(:,IN,:) -             &
                         TPHALO2%NORTH(:,:)) * (0.5-SIGN(0.5,PRVCT(:,IN+1,:))) )
  END IF
!
!       USE A THIRD ORDER UPSTREAM SCHEME ELSEWHERE 
!
   PR(:,ISF:INF,:) = 1./6. * ( (2.*PSRC(:,IS:IN,:) + 5.*PSRC(:,IS-1:IN-1,:) -  &
                       PSRC(:,IS-2:IN-2,:)) * (0.5+SIGN(0.5,PRVCT(:,IS:IN,:))) &
                            + (5.*PSRC(:,IS:IN,:) + 2.*PSRC(:,IS-1:IN-1,:) -   &
                       PSRC(:,IS+1:IN+1,:)) * (0.5-SIGN(0.5,PRVCT(:,IS:IN,:))) )
!
END SELECT
!
PR = PR * PRVCT
!
END SUBROUTINE ADVEC_3RD_ORDER_MY
!
!-------------------------------------------------------------------------------
!
!     #######################################
      FUNCTION UP3_WZ(PSRC, PRWCT) RESULT(PR)
!     #######################################
!!
!!****  UP3_WZ - upstream fluxes of W in Z direction
!!              input variable PSRC is on W grid, and output PR is on MASS grid
!!
!!    AUTHOR
!!    ------
!!      C.Lac            * CNRM/GMME *               
!!
!!    MODIFICATIONS
!!    -------------
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on W grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on MASS GRID
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB    ! Begining useful area in x,y,z directions
INTEGER :: IKE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
! upstream flux on mass points
!
PR(:,:,IKB:IKE-1) = 1./6. * ( (2.*PSRC(:,:,IKB+1:IKE) + 5.*PSRC(:,:,IKB:IKE-1)-&
                 PSRC(:,:,IKB-1:IKE-2)) * (0.5+SIGN(0.5,PRWCT(:,:,IKB:IKE-1))) &
                            + (5.*PSRC(:,:,IKB+1:IKE) + 2.*PSRC(:,:,IKB:IKE-1)-&
                 PSRC(:,:,IKB+2:IKE+1)) * (0.5-SIGN(0.5,PRWCT(:,:,IKB:IKE-1))) )
!
PR(:,:,IKB-1) = PSRC(:,:,IKB-1) * (0.5+SIGN(0.5,PRWCT(:,:,IKB-1))) &
              + PSRC(:,:,IKB  ) * (0.5-SIGN(0.5,PRWCT(:,:,IKB-1)))
PR(:,:,IKE  ) = PSRC(:,:,IKE  ) * (0.5+SIGN(0.5,PRWCT(:,:,IKE  ))) &
              + PSRC(:,:,IKE+1) * (0.5-SIGN(0.5,PRWCT(:,:,IKE  )))
PR(:,:,IKE+1) = -999.                  
!
PR = PR * PRWCT
!
END FUNCTION UP3_WZ
!
!-------------------------------------------------------------------------------
!
!     #######################################
      FUNCTION UP3_MZ(PSRC, PRWCT) RESULT(PR)
!     #######################################
!!
!!****  UP3_MZ - upstream fluxes of variable in Z direction
!!      input variable PSRC is on MASS grid, and output PR is on W grid
!!
!!    AUTHOR
!!    ------
!!      C.Lac            * CNRM/GMME *               
!!
!!    MODIFICATIONS
!!    -------------
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PSRC  ! variable on MASS grid at t
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRWCT ! contrav. comp. on W grid
!
! output source term
!
REAL, DIMENSION(SIZE(PSRC,1),SIZE(PSRC,2),SIZE(PSRC,3)) :: PR
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB    ! Begining useful area in x,y,z directions
INTEGER :: IKE    ! End useful area in x,y,z directions
!
!-------------------------------------------------------------------------------
!
!*       0.3.     COMPUTES THE DOMAIN DIMENSIONS
!                 ------------------------------
!
IKB = 1 + JPVEXT
IKE = SIZE(PSRC,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
! upstream flux on mass points
!
PR(:,:,IKB+1:IKE) = 1./6. * ( (2.*PSRC(:,:,IKB+1:IKE) + 5.*PSRC(:,:,IKB:IKE-1)-&
                 PSRC(:,:,IKB-1:IKE-2)) * (0.5+SIGN(0.5,PRWCT(:,:,IKB+1:IKE))) &
                            + (5.*PSRC(:,:,IKB+1:IKE) + 2.*PSRC(:,:,IKB:IKE-1)-&
                 PSRC(:,:,IKB+2:IKE+1)) * (0.5-SIGN(0.5,PRWCT(:,:,IKB+1:IKE))) )
!
PR(:,:,IKB  ) = PSRC(:,:,IKB-1) * (0.5+SIGN(0.5,PRWCT(:,:,IKB  ))) &
              + PSRC(:,:,IKB  ) * (0.5-SIGN(0.5,PRWCT(:,:,IKB  )))
PR(:,:,IKE+1) = PSRC(:,:,IKE  ) * (0.5+SIGN(0.5,PRWCT(:,:,IKE+1))) &
              + PSRC(:,:,IKE+1) * (0.5-SIGN(0.5,PRWCT(:,:,IKE+1)))
PR(:,:,IKB-1) = -999.                  
!
PR = PR * PRWCT
!
END FUNCTION UP3_MZ
