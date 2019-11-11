!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!#####################
MODULE MODI_ADV_BOUNDARIES
!#####################
!
INTERFACE
!
      SUBROUTINE ADV_BOUNDARIES ( HLBCX,HLBCY,PFIELD,PFIELDI,HFIELD )    
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)      :: HLBCX,HLBCY   ! X and Y-direc. LBC type
REAL, DIMENSION(:,:,:),   INTENT(INOUT)         :: PFIELD                        
REAL, DIMENSION(:,:,:),   INTENT(IN), OPTIONAL  :: PFIELDI                        
CHARACTER(LEN=1),         INTENT(IN), OPTIONAL  :: HFIELD  ! Field type
!
END SUBROUTINE ADV_BOUNDARIES
!
END INTERFACE
!

END MODULE MODI_ADV_BOUNDARIES
!
!
!     ####################################################################
      SUBROUTINE ADV_BOUNDARIES ( HLBCX,HLBCY,PFIELD,PFIELDI,HFIELD )
!     ####################################################################
!
!!****  *ADV_BOUNDARIES* - routine to prepare the top and bottom Boundary Conditions 
!!
!!
!!    AUTHOR
!!    ------
!!   V.Masson
!! Correction :	
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!   C.Lac     : 10/16 : top BC for W 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!         
USE MODD_PARAMETERS
USE MODE_ll
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)      :: HLBCX,HLBCY   ! X and Y-direc. LBC type
REAL, DIMENSION(:,:,:),   INTENT(INOUT)         :: PFIELD                        
REAL, DIMENSION(:,:,:),   INTENT(IN), OPTIONAL  :: PFIELDI                        
CHARACTER(LEN=1),         INTENT(IN), OPTIONAL  :: HFIELD  ! Field type
!
!
!*       0.2   declarations of local variables
!
INTEGER             :: IKB       ! indice K Beginning in z direction
INTEGER             :: IKE       ! indice K End       in z direction 
INTEGER             :: IIU, IJU  ! Index End in X and Y directions
INTEGER             :: IIB,IIE,IJB,IJE ! interior domaine bound
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES:
!              ----------------------------------------------
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PFIELD,3) - JPVEXT
IIU=SIZE(PFIELD,1)
IJU=SIZE(PFIELD,2)
!
IF (SIZE(PFIELD)==0) RETURN
!-------------------------------------------------------------------------------
!
!*       2.    UPPER AND LOWER BC FILLING:   
!              ---------------------------
!
!*       2.1    COMPUTE THE FIELD EXTRAPOLATIONS AT THE GROUND
!
!
   IF (PRESENT(HFIELD) .AND. PRESENT(PFIELDI)) THEN
     IF (HFIELD=='W') &
     PFIELD  (:,:,IKB  )   = PFIELDI (:,:,IKB) 
   END IF
!
   PFIELD  (:,:,IKB-1)   = PFIELD  (:,:,IKB) 

!
!*       2.2    COMPUTE THE FIELD EXTRAPOLATIONS AT THE TOP
!
  PFIELD  (:,:,IKE+1)   = PFIELD  (:,:,IKE) 
!
!
!*       3.    LATERAL BC FILLING                                
!              ---------------------------
!
IF( PRESENT(PFIELDI) )  THEN
  IF (HLBCX(1)=='OPEN' .AND. LWEST_ll()) THEN
     PFIELD(:IIB-1,:,:) = PFIELDI(:IIB-1,:,:) ! 1
     IF (PRESENT(HFIELD)) THEN
       IF (HFIELD=='U') &
       PFIELD(:IIB,:,:) = PFIELDI(:IIB,:,:)   ! 2
     END IF
  END IF
  IF (HLBCX(2)=='OPEN' .AND. LEAST_ll()) THEN
     PFIELD(IIE+1:,:,:) = PFIELDI(IIE+1:,:,:) ! IIU
  END IF
  IF (HLBCY(1)=='OPEN' .AND. LSOUTH_ll()) THEN
     PFIELD(:,:IJB-1,:) = PFIELDI(:,:IJB-1,:) ! 1
     IF (PRESENT(HFIELD)) THEN
       IF (HFIELD=='V') &
       PFIELD(:,:IJB,:) = PFIELDI(:,:IJB,:) ! 2
     END IF
  END IF
  IF (HLBCY(2)=='OPEN' .AND. LNORTH_ll()) THEN
     PFIELD(:,IJE+1:,:) = PFIELDI(:,IJE+1:,:) ! IJU
  END IF
END IF
!
!*       4. TOP BC for W
!
 IF (PRESENT(HFIELD)) THEN
   IF (HFIELD=='W') PFIELD(:,:,IKE+1) = 0.
 END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADV_BOUNDARIES
