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
!     ################
      MODULE MODI_GDIV
!     ################
!
INTERFACE
!
      SUBROUTINE GDIV(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PU,PV,PW,PGDIV)
!  
IMPLICIT NONE
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type
! 
                                                 ! Metric coefficients:       
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! d*xx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! d*yy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX      ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY      ! d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ      ! d*zz
!
                                                 ! Field components
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PU        ! along x             
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PV        ! along y
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PW        ! along z
!
REAL, DIMENSION(:,:,:), INTENT(OUT)    :: PGDIV     ! divergence at 
                                                    ! a mass point
!
END SUBROUTINE GDIV 
!
END INTERFACE
!
END MODULE MODI_GDIV
!
!     ####################################################################
      SUBROUTINE GDIV(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PU,PV,PW,PGDIV)
!     ####################################################################
!
!!****  *GDIV * - Compute J times the divergence of 1/J times a vector defined 
!!      by its cartesian components 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute J times the divergence of 
!     1/J times a vector which cartesian components are (U, V, W). The result
!     is localized at a mass point: 
!    
!                    GDIV = dxf (UC) + dyf (VC) + dzf (WC)
!       
!     where UC, VC, WC are the contravariant components of the vector.
!     The array is completed outside the physical domain by the value of the
!     normal component at the boundary. 
!
!!**  METHOD
!!    ------
!!      First, we compute the contravariant components by using 
!!    the suboutine CONTRAV (The metric coefficients are dummy arguments). Then 
!!    we use the Shuman finite difference operators DXF, DYF, DZF to compute 
!!    the divergence. The result is localized at a mass point.
!!
!!    EXTERNAL
!!    --------
!!      SUBROUTINE CONTRAV : compute the contavariant components 
!!    Shuman operators :
!!      FUNCTION DXF : compute finite difference along x for a variable 
!!    localized at a flux side
!!      FUNCTION DYF : compute finite difference along y for a variable 
!!    localized at a flux side
!!      FUNCTION DZF : compute finite difference along z for a variable 
!!    localized at a flux side
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT, JPVEXT: define the number of marginal points out of the 
!!        physical domain along horizontal and vertical directions respectively
!!     Module MODD_CONF: model configurations
!!        L2D: logical for 2D model version
!!      Module MODI_SHUMAN : interface for the Shuman operators
!!      Module MODI_CONTRAV : interface for the contravariant components
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine GDIV)
!!      
!!
!!    AUTHOR
!!    ------
!!	P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/07/94 
!!                  17/07/97  ( J. Stein and V. Masson) initialize the corner
!!                             verticals
!!                  30/09/97  ( J. Stein ) bug correction for the case of
!!                             non-vanishing oro. at the open lbc
!!     Modification 15/06/98  (D.Lugato, R.Guivarch) Parallelisation
!!                  22/08/02  (P Jabouille) simplification of parallel coding
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CONF
USE MODI_SHUMAN
USE MODI_CONTRAV
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
! 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX    ! x-direction LBC type 
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY    ! y-direction LBC type
! 
                                                 ! Metric coefficients:       
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! d*xx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! d*yy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX      ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY      ! d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ      ! d*zz
!
                                                 ! Field components
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PU        ! along x             
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PV        ! along y
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PW        ! along z
!
REAL, DIMENSION(:,:,:), INTENT(OUT)     :: PGDIV             ! divergence at 
                                                             ! a mass point
!
!*       0.2   declarations of local variables
!
                                         ! Contravariant components along:
REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2),SIZE(PU,3)) :: ZUC      ! x
REAL, DIMENSION(SIZE(PV,1),SIZE(PV,2),SIZE(PV,3)) :: ZVC      ! y 
REAL, DIMENSION(SIZE(PW,1),SIZE(PW,2),SIZE(PW,3)) :: ZWC, &   ! z 
                                                     Z1,Z2,Z3 !work arrays
!
INTEGER :: IIB          ! indice I for the first inner mass point along x
INTEGER :: IIE          ! indice I for the last inner mass point along x
INTEGER :: IJB          ! indice J for the first inner mass point along y
INTEGER :: IJE          ! indice J for the last inner mass point along y
INTEGER :: IKB          ! indice K for the first inner mass point along z
INTEGER :: IKE          ! indice K for the last inner mass point along z
!
INTEGER :: JI,JJ,JK                         ! loop indexes
!
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE LOOP BOUNDS
!              -------------------
!
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PU,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE CONTRAVARIANT COMPONENTS
!              ------------------------------------
!
!*       2.1   prepare the boundary conditions
!
!
    PU(:,:,IKB-1)=PU(:,:,IKB)
    PU(:,:,IKE+1)=PU(:,:,IKE)
    PV(:,:,IKB-1)=PV(:,:,IKB)
    PV(:,:,IKE+1)=PV(:,:,IKE)

! 
!
!*       2.1   compute the contravariant components
!
CALL CONTRAV(HLBCX,HLBCY,PU,PV,PW,PDXX,PDYY,PDZZ,PDZX,PDZY,ZUC,ZVC,ZWC,4)
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTE THE DIVERGENCE 
!              ----------------------
!
PGDIV=0. !usefull for the four corners and halo zones
!
Z1(IIB:IIE,:,:)=ZUC(IIB+1:IIE+1,:,:)-ZUC(IIB:IIE,:,:)
Z2(:,IJB:IJE,:)=ZVC(:,IJB+1:IJE+1,:)-ZVC(:,IJB:IJE,:)
Z3(:,:,IKB:IKE)=ZWC(:,:,IKB+1:IKE+1)-ZWC(:,:,IKB:IKE)
!
PGDIV(IIB:IIE,IJB:IJE,IKB:IKE)= Z1(IIB:IIE,IJB:IJE,IKB:IKE) +  &
                                Z2(IIB:IIE,IJB:IJE,IKB:IKE) +  &
                                Z3(IIB:IIE,IJB:IJE,IKB:IKE) 
                               ! only the divergences computed 
                               ! in the inner mass points are meaningful  
!
!-------------------------------------------------------------------------------
!
!*       4.    SET DIVERGENCE AT THE OUTER POINTS
!              ----------------------------------
!
!*       4.1   set divergence at the upper and lower boundary
!
!   we set the divergence equal to the vertical contravariant component above
!   and under the physical domain
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    PGDIV(JI,JJ,IKB-1)=ZWC(JI,JJ,IKB)
    PGDIV(JI,JJ,IKE+1)=ZWC(JI,JJ,IKE+1)
  END DO
END DO
!
!*       4.2   set divergence at the lateral boundaries
!
!   we set the divergence equal to the horizontal contravariant component at 
!   the right and the left of the physical domain in both horizontal directions
!   for non-periodic cases
!
IF(HLBCX(1) /= 'CYCL' .AND. LWEST_ll()) THEN
  DO JK=IKB,IKE
    DO JJ=IJB,IJE
      PGDIV(IIB-1,JJ,JK)=ZUC(IIB,JJ,JK)
    END DO
  END DO
END IF
!
IF(HLBCX(2) /= 'CYCL' .AND. LEAST_ll()) THEN
  DO JK=IKB,IKE
    DO JJ=IJB,IJE
      PGDIV(IIE+1,JJ,JK)=ZUC(IIE+1,JJ,JK)
    END DO
  END DO
END IF
!
!
IF (.NOT. L2D .AND. HLBCY(1) /= 'CYCL' .AND. LSOUTH_ll()) THEN
  DO JK=IKB,IKE
    DO JI=IIB,IIE
      PGDIV(JI,IJB-1,JK)=ZVC(JI,IJB,JK)
    END DO
  END DO
END IF
!
IF (.NOT. L2D .AND. HLBCY(2) /= 'CYCL' .AND. LNORTH_ll()) THEN
  DO JK=IKB,IKE
    DO JI=IIB,IIE
      PGDIV(JI,IJE+1,JK)=ZVC(JI,IJE+1,JK)
    END DO
  END DO
END IF
!
!*       4.3   set divergence at the corner points  
!
! it is the following of the condition of copy the horizontal component
! under the bottom of the model
!
IF(HLBCX(1) /= 'CYCL' .AND. LWEST_ll()) THEN
  PGDIV(IIB-1,IJB:IJE,IKB-1)=PGDIV(IIB-1,IJB:IJE,IKB)
  PGDIV(IIB-1,IJB:IJE,IKE+1)=PGDIV(IIB-1,IJB:IJE,IKE)
END IF
!
IF (HLBCX(2) /= 'CYCL' .AND. LEAST_ll()) THEN
  PGDIV(IIE+1,IJB:IJE,IKB-1)=PGDIV(IIE+1,IJB:IJE,IKB)
  PGDIV(IIE+1,IJB:IJE,IKE+1)=PGDIV(IIE+1,IJB:IJE,IKE)
END IF
!
IF (.NOT. L2D .AND. HLBCY(1) /= 'CYCL' .AND. LSOUTH_ll()) THEN
  PGDIV(IIB:IIE,IJB-1,IKB-1)=PGDIV(IIB:IIE,IJB-1,IKB)
  PGDIV(IIB:IIE,IJB-1,IKE+1)=PGDIV(IIB:IIE,IJB-1,IKE)
END IF
!
IF (.NOT. L2D .AND. HLBCY(2) /= 'CYCL' .AND. LNORTH_ll()) THEN
  PGDIV(IIB:IIE,IJE+1,IKB-1)=PGDIV(IIB:IIE,IJE+1,IKB)
  PGDIV(IIB:IIE,IJE+1,IKE+1)=PGDIV(IIB:IIE,IJE+1,IKE)
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GDIV
