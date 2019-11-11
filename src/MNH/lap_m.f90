!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #################
      MODULE MODI_LAP_M
!     #################
!
INTERFACE
!
      FUNCTION LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PY)  &
               RESULT(PLAP_M)
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHODJ    ! density_reference * J
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY        ! field components
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PLAP_M ! final divergence 
!
END FUNCTION LAP_M
!
END INTERFACE
!
END MODULE MODI_LAP_M
!
!
!
!     #########################################################################
      FUNCTION LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PY)  &
               RESULT(PLAP_M)
!     #########################################################################
!
!!****  *LAP_M * - compute the Laplacian of a field PY 
!!
!!    PURPOSE
!!    -------
!       This function computes laplacian of a scalar field PY
!     localized at mass points, with bottom topography.
!     The result is localized at a mass point and defined by: 
!                                   (         ( GX_M_U (PY) )  )
!                    PLAP_M  = GDIV ( rho * J ( GX_M_V (PY) )  )
!                                   (         ( GX_M_W (PY) )  )  
!
!     Where GX_M_.. are the cartesian components of the gradient of PY and
!     GDIV is the operator acting on a vector AA: 
!                   
!                   GDIV ( AA ) = J * divergence (1/J  AA  ) 
!     
!!**  METHOD
!!    ------
!!      First, we compute the gradients along x, y , z of the P field. The 
!!    result is multiplied by rhod.
!!      Then, the pseudo-divergence ( J * DIV (1/J o ) ) is computed by the 
!!    subroutine GDIV. The result is localized at a mass point.
!!
!!    EXTERNAL
!!    --------
!!      Function GX_M_U : compute the gradient along x 
!!      Function GY_M_V : compute the gradient along y 
!!      Function GZ_M_W : compute the gradient along z 
!!      FUNCTION MXM: compute an average in the x direction for a variable  
!!      at a mass localization
!!      FUNCTION MYM: compute an average in the y direction for a variable  
!!      at a mass localization
!!      FUNCTION MZM: compute an average in the z direction for a variable  
!!      at a mass localization
!!      Subroutine GDIV : compute J times the divergence of 1/J times a vector
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS: JPHEXT, JPVEXT
!!      Module MODD_CONF: L2D
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/02/07 simplified from QLAP function, T.Maric 
!!      Modification   
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
USE MODD_CONF
!USE MODD_CST
USE MODI_GDIV
!USE MODI_GDIV_M
USE MODI_GRADIENT_M
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
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
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PRHODJ    ! density of reference * J
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY        ! field components
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PLAP_M ! final divergence 
!
!*       0.2   declarations of local variables
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: ZU ! rho*J*gradient along x
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: ZV ! rho*J*gradient along y
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: ZW ! rho*J*gradient along z
!
INTEGER                          :: IIU,IJU,IKU         ! I,J,K array sizes
INTEGER                          :: JK,JJ,JI            ! vertical loop index
!-------------------------------------------------------------------------------
!
!
!*       1.    COMPUTE THE GRADIENT COMPONENTS
!              -------------------------------
!
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=SIZE(PY,3)
!
ZU = GX_M_U(1,IKU,1,PY,PDXX,PDZZ,PDZX)
!
IF ( HLBCX(1) /= 'CYCL' .AND. LWEST_ll() ) THEN
  DO JK=2,IKU-1
    DO JJ=1,IJU
      ZU(2,JJ,JK)=  (PY(2,JJ,JK) - PY(1,JJ,JK) - 0.5 * (                  &
       PDZX(2,JJ,JK)   * (PY(2,JJ,JK)-PY(2,JJ,JK-1)) / PDZZ(2,JJ,JK)      &
      +PDZX(2,JJ,JK+1) * (PY(2,JJ,JK+1)-PY(2,JJ,JK)) / PDZZ(2,JJ,JK+1)    &  
                                      )    ) / PDXX(2,JJ,JK)
    END DO
  END DO
END IF
!
IF ( HLBCX(1) /= 'CYCL' .AND. LEAST_ll() ) THEN
  DO JK=2,IKU-1
    DO JJ=1,IJU
      ZU(IIU,JJ,JK)=  (PY(IIU,JJ,JK) - PY(IIU-1,JJ,JK) - 0.5 * (                    &
        PDZX(IIU,JJ,JK)   * (PY(IIU-1,JJ,JK)-PY(IIU-1,JJ,JK-1)) / PDZZ(IIU-1,JJ,JK)  &
       +PDZX(IIU,JJ,JK+1) * (PY(IIU-1,JJ,JK+1)-PY(IIU-1,JJ,JK)) / PDZZ(IIU-1,JJ,JK+1)&  
                                            ) ) / PDXX(IIU,JJ,JK)
    END DO
  END DO
END IF
!
IF(.NOT. L2D) THEN 
!
  ZV = GY_M_V(1,IKU,1,PY,PDYY,PDZZ,PDZY)
!
  IF ( HLBCY(1) /= 'CYCL' .AND. LSOUTH_ll() ) THEN 
    DO JK=2,IKU-1
      DO JI=1,IIU
        ZV(JI,2,JK)=   (PY(JI,2,JK) - PY(JI,1,JK) - 0.5 * (                  &
          PDZY(JI,2,JK)   * (PY(JI,2,JK)-PY(JI,2,JK-1)) / PDZZ(JI,2,JK)      &
         +PDZY(JI,2,JK+1) * (PY(JI,2,JK+1)-PY(JI,2,JK)) / PDZZ(JI,2,JK+1)    &  
                                            )   ) / PDYY(JI,2,JK) 
      END DO
    END DO
  END IF
  IF ( HLBCY(1) /= 'CYCL' .AND. LNORTH_ll() ) THEN
!
    DO JK=2,IKU-1
      DO JI=1,IIU
        ZV(JI,IJU,JK)=    (PY(JI,IJU,JK) - PY(JI,IJU-1,JK) - 0.5 * (                  &
          PDZY(JI,IJU,JK)   * (PY(JI,IJU-1,JK)-PY(JI,IJU-1,JK-1)) / PDZZ(JI,IJU-1,JK)  &
         +PDZY(JI,IJU,JK+1) * (PY(JI,IJU-1,JK+1)-PY(JI,IJU-1,JK)) / PDZZ(JI,IJU-1,JK+1)&
                                                )  ) / PDYY(JI,IJU,JK) 
      END DO
    END DO
  END IF
!
ELSE
  ZV=0.
ENDIF
!
ZU = MXM(PRHODJ) * ZU
!
IF(.NOT. L2D) THEN 
   ZV = MYM(PRHODJ) * ZV
ENDIF
!
ZW = MZM(1,IKU,1,PRHODJ) * GZ_M_W(1,IKU,1,PY,PDZZ)
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE DIVERGENCE  
!              ----------------------
!
CALL GDIV(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,ZU,ZV,ZW,PLAP_M)    
!
PLAP_M(:,:,1) = 0.
!
!-------------------------------------------------------------------------------
!
END FUNCTION LAP_M
