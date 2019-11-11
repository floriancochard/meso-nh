!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ######spl
      MODULE MODI_MEAN_Z
!     ###################### 
!
INTERFACE
!
!
FUNCTION MW_Z(PA,PRHODJ)      RESULT(PMW_Z)
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PMW_Z ! result mass point
!
END FUNCTION MW_Z
!
END INTERFACE 
!
END MODULE MODI_MEAN_Z
!-----------------------------------------------------------------
!     #######################################################
      FUNCTION MW_Z(PA,PRHODJ)      RESULT(PMW_Z)
!     #######################################################
!
!!****  *MW_Z* - compute the vertical average weighted mw_z= SUM(PA*RHO)/ (sum(RHO))
!
!!    AUTHOR
!!    ------
!!      P.Peyrille CNRM
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/02/04 
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN
USE MODD_CONF
USE MODD_PARAMETERS ! (JPVEXT)

IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: PMW_Z ! result mass point
!
!
!*       0.2   declaration of local variables
REAL,DIMENSION(SIZE(PA,1),SIZE(PA,2))     :: ZSOM,ZSOMRHO
INTEGER :: JK,IKB,IKE,IKU
INTEGER             :: IRESP   ! Return code of FM routines
INTEGER :: ZSTOP
!              NONE
!
!----------------------------------------------------------------------------

!
!*       1.    DEFINITION of GX_M_M
!              --------------------
!
IKU=SIZE(PA,3) ! nb points sur Z
IKB=1+JPVEXT
IKE=IKU-JPVEXT

ZSOM(:,:)=0.
ZSOMRHO(:,:)=0.
DO JK=IKB,IKE
  ZSOM(:,:) = ZSOM(:,:) + PA(:,:,JK) * PRHODJ(:,:,JK)
  ZSOMRHO(:,:) = ZSOMRHO(:,:) + PRHODJ(:,:,JK)
ENDDO
!
PMW_Z(:,:)= ZSOM(:,:)/ZSOMRHO(:,:)
!
!----------------------------------------------------------------------------
!
END FUNCTION MW_Z
