!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ######spl
MODULE MODI_BEAMBROAD

  INTERFACE 
  SUBROUTINE BEAMBROAD(PDISCR,PSINGPT,PX_H,PX_V,PW_H,PW_V,OMASK)
    REAL,DIMENSION(:,:,:,:,:,:),INTENT(INOUT)  :: PDISCR  
    REAL,DIMENSION(:,:,:,:),  INTENT(OUT) :: PSINGPT
    
    REAL, DIMENSION(:),           INTENT(IN) :: PX_H ! Gaussian horizontal nodes
    REAL, DIMENSION(:),           INTENT(IN) :: PX_V ! Gaussian vertical nodes
    REAL, DIMENSION(:),           INTENT(IN) :: PW_H ! Gaussian horizontal weights
    REAL, DIMENSION(:),           INTENT(IN) :: PW_V ! Gaussian vertical weights
    LOGICAL,                      INTENT(IN) :: OMASK ! compute bins located after partial mask
  END SUBROUTINE BEAMBROAD
  END INTERFACE
  
END MODULE MODI_BEAMBROAD
!
!     ##############################################################
      SUBROUTINE BEAMBROAD(PDISCR,PSINGPT,PX_H,PX_V,PW_H,PW_V,OMASK)
!     ##############################################################
!
!!****  *BEAMBROAD * - takes into account beam broadening with range
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute variables on the beam ray from
!!    discretized pinpoint rays.
!!
!!**  METHOD
!!    ------
!!      Book2 of documentation ( routine RADAR_SIMULATOR )
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!      Module MODD_PARAMETERS
!!      Module MODD_RADAR 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine RADAR_SIMULATOR )
!!
!!    AUTHOR
!!    ------
!!      O. Caumont       * Météo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/04/2008
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST       , ONLY: XPI
USE MODD_PARAMETERS, ONLY: XUNDEF
USE MODD_RADAR     , ONLY: LQUAD,NBELEV

IMPLICIT NONE
    
REAL,DIMENSION(:,:,:,:,:,:),INTENT(INOUT)  :: PDISCR  
REAL,DIMENSION(:,:,:,:),  INTENT(OUT) :: PSINGPT

REAL, DIMENSION(:),           INTENT(IN) :: PX_H ! Gaussian horizontal nodes
REAL, DIMENSION(:),           INTENT(IN) :: PX_V ! Gaussian vertical nodes
REAL, DIMENSION(:),           INTENT(IN) :: PW_H ! Gaussian horizontal weights
REAL, DIMENSION(:),           INTENT(IN) :: PW_V ! Gaussian vertical weights
LOGICAL,                      INTENT(IN) :: OMASK ! compute bins located after partial mask

INTEGER  :: JI,JL,JEL,JAZ,JH,JV ! Loop control variables
INTEGER  :: IEL
INTEGER  :: INBRAD,INPTS_H,INPTS_V ! sizes of the arrays
REAL :: ZVTEMP
    
!
!*       1.     INITIALIZATION 
!   	        --------------
INBRAD=SIZE(PDISCR,1)
INPTS_H=SIZE(PDISCR,5)
INPTS_V=SIZE(PDISCR,6)
PSINGPT(:,:,:,:)=0.

!
!*       2.     CALCULATIONS 
!   	        --------------
DO JI=1,INBRAD  
   IEL=NBELEV(JI)
   DO JEL=1,IEL  
      DO JAZ=1,SIZE(PDISCR,3) 
         DO JL=1,SIZE(PDISCR,4)
            DO JH=1,INPTS_H
               ZVTEMP=0.
               DO JV=1,INPTS_V  ! Loop on Jv
                  ! if previously underground on this beam
                  IF(JL > 1) THEN
                     IF(PDISCR(JI,JEL,JAZ,JL-1,JH,JV)==-XUNDEF.AND..NOT.OMASK) &
                     PDISCR(JI,JEL,JAZ,JL,JH,JV)=-XUNDEF
                  END IF
                  IF(PDISCR(JI,JEL,JAZ,JL,JH,JV) /= -XUNDEF.AND.ZVTEMP /= -XUNDEF) THEN
                     ! Quadrature on vertical reflectivities 
                     IF(LQUAD) THEN
                        ZVTEMP=ZVTEMP+PDISCR(JI,JEL,JAZ,JL,JH,JV)*PW_V(ABS((2*JV-INPTS_V-1)/2)+1) &
                             *EXP(-2.*LOG(2.)*PX_V(ABS((2*JV-INPTS_V-1)/2)+1)**2)
                     ELSE
                        ZVTEMP=ZVTEMP+PDISCR(JI,JEL,JAZ,JL,JH,JV)*PW_V(ABS((2*JV-INPTS_V-1)/2)+1)
                     END IF
                  ELSE
                     ZVTEMP=-XUNDEF
                  END IF
               END DO ! End loop on JV
               
               IF(ZVTEMP /= -XUNDEF .AND. PSINGPT(JI,JEL,JAZ,JL) /= -XUNDEF) THEN
                  IF(LQUAD) THEN
                     PSINGPT(JI,JEL,JAZ,JL)=PSINGPT(JI,JEL,JAZ,JL)+ZVTEMP*PW_H(ABS((2*JH-INPTS_H-1)/2)+1) &
                          *EXP(-2.*LOG(2.)*PX_H(ABS((2*JH-INPTS_H-1)/2)+1)**2)
                  ELSE
                     PSINGPT(JI,JEL,JAZ,JL)=PSINGPT(JI,JEL,JAZ,JL)+ZVTEMP*PW_H(ABS((2*JH-INPTS_H-1)/2)+1)
                  END IF
               ELSE 
                  PSINGPT(JI,JEL,JAZ,JL)=-XUNDEF
               END IF
            END DO ! End loop on JH 
                
            IF(PSINGPT(JI,JEL,JAZ,JL) /= -XUNDEF) THEN
               IF(LQUAD) THEN
                  PSINGPT(JI,JEL,JAZ,JL)=PSINGPT(JI,JEL,JAZ,JL)*2.*LOG(2.)/XPI
               ELSE
                  PSINGPT(JI,JEL,JAZ,JL)=PSINGPT(JI,JEL,JAZ,JL)/XPI! ELSE REMAINS -XUNDEF
               END IF
            END IF
            
         END DO
      END DO
   END DO
END DO

END SUBROUTINE BEAMBROAD
