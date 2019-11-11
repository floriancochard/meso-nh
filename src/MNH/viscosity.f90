!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
! 
!     #####################
      MODULE MODI_VISCOSITY
!     #####################
!
INTERFACE
!
    SUBROUTINE VISCOSITY(HLBCX, HLBCY, KRR, KSV, PNU, PPRANDTL,          &
                        OVISC_UVW, OVISC_TH, OVISC_SV, OVISC_R,         &
                        ODRAG,  &
                        PUT, PVT, PWT, PTHT, PRT, PSVT,                 &
                        PRHODJ, PDXX, PDYY, PDZZ, PDZX, PDZY,           &
                        PRUS, PRVS, PRWS, PRTHS, PRRS, PRSVS,PDRAG      )
!
     IMPLICIT NONE
!
!*       0.1   Declarations of arguments:
!
! X and Y lateral boundary conditions
     CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY 
!
     INTEGER, INTENT(IN) :: KRR      ! number of moist variables
     INTEGER, INTENT(IN) :: KSV      ! number of scalar variables
!
     REAL, INTENT(IN) :: PNU         ! viscosity coefficient
     REAL, INTENT(IN) :: PPRANDTL    ! Parandtl number
!
!
! logical switches
     LOGICAL, INTENT(IN) :: OVISC_UVW ! momentum
     LOGICAL, INTENT(IN) :: OVISC_TH  ! theta
     LOGICAL, INTENT(IN) :: OVISC_SV  ! scalar tracer
     LOGICAL, INTENT(IN) :: OVISC_R   ! moisture
     LOGICAL, INTENT(IN) :: ODRAG     ! noslip/freeslip 
!
!
! input variables at time t
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PUT
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PVT
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PWT
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT
     REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT
     REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PSVT
!
!
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ ! rho_ref * Jacobian
!
      REAL, DIMENSION(:,:), INTENT(IN) :: PDRAG ! Array -1/1 defining where the no-slipcondition is applied
! metric coefficients
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDYY
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZX
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZY
!
! output source terms
     REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRUS, PRVS, PRWS
     REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRTHS
     REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS, PRSVS
!
   END SUBROUTINE VISCOSITY
!
END INTERFACE
!
END MODULE MODI_VISCOSITY
!
!-------------------------------------------------------------------------------
!
SUBROUTINE VISCOSITY(HLBCX, HLBCY, KRR, KSV, PNU, PPRANDTL,          &
                     OVISC_UVW, OVISC_TH, OVISC_SV, OVISC_R,         &
                     ODRAG,        &
                     PUT, PVT, PWT, PTHT, PRT, PSVT,                 &
                     PRHODJ, PDXX, PDYY, PDZZ, PDZX, PDZY,           &
                     PRUS, PRVS, PRWS, PRTHS, PRRS, PRSVS,PDRAG      )
!
!    IMPLICIT ARGUMENTS
!    ------------------ 
!      Module MODD_PARAMETERS: JPHEXT, JPVEXT
!      Module MODD_CONF: L2D
!
!!    AUTHOR
!!    ------
!!      Jeanne Colin                      * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      01/18 (C.Lac) Add budgets
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODI_LAP_M
  USE MODI_SHUMAN  
  USE MODD_PARAMETERS
  USE MODD_CONF
  USE MODD_VISCOSITY
  USE MODD_DRAG_n
  USE MODD_BUDGET
  USE MODE_ll
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
  USE MODE_FM
  USE MODI_BUDGET
!
!-------------------------------------------------------------------------------
!
  IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
! X and Y lateral boundary conditions
  CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY 
!
  INTEGER, INTENT(IN) :: KRR      ! number of moist variables
  INTEGER, INTENT(IN) :: KSV      ! number of scalar variables
!
  REAL, INTENT(IN) :: PPRANDTL    ! Parandtl number
  REAL, INTENT(IN) :: PNU         ! viscous diffusion rate 
!
! logical switches
     LOGICAL, INTENT(IN) :: OVISC_UVW ! momentum
     LOGICAL, INTENT(IN) :: OVISC_TH  ! theta
     LOGICAL, INTENT(IN) :: OVISC_SV  ! scalar tracer
     LOGICAL, INTENT(IN) :: OVISC_R   ! moisture
     LOGICAL, INTENT(IN) :: ODRAG     ! noslip/freeslip
!
! input variables at time t
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PUT
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PVT
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PWT
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT
     REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRT
     REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PSVT
!
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODJ ! rho_ref * Jacobian
!
!
REAL, DIMENSION(:,:), INTENT(IN) :: PDRAG ! Array -1/1 defining where the no-slip condition is applied

     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDYY
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZX
     REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZY
!
! output source terms
     REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRUS, PRVS, PRWS
     REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRTHS
     REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS, PRSVS
!
!
!*       0.2   Declarations of local variables
!
  INTEGER :: IK ! counter
  INTEGER :: IKB, IKE
!
  REAL, DIMENSION(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)) :: ZTMP ! temp storage
  REAL, DIMENSION(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)) :: ZLAPu ! temp storage
  REAL, DIMENSION(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)) :: ZLAPv ! temp storage
  REAL, DIMENSION(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)) :: ZY1 ! temp storage
  REAL, DIMENSION(SIZE(PWT,1),SIZE(PWT,2),SIZE(PWT,3)) :: ZY2 ! temp storage
!
!
INTEGER                          :: IIU,IJU,IKU         ! I,J,K array sizes
!
INTEGER                          :: JI,JJ,JK  ! I loop index
INTEGER :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!
!
!-------------------------------------------------------------------------------
!
IIU=SIZE(PWT,1)
IJU=SIZE(PWT,2)
IKU=SIZE(PWT,3)

!*       1.    Viscous forcing for potential temperature
!	       -----------------------------------------
!
!
IF (OVISC_TH) THEN
!
!
      PRTHS = PRTHS + PNU/PPRANDTL * &
              LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PTHT)
!
!
END IF
!
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'VISC_BU_RU')
!
!-------------------------------------------------------------------------------
!
!*       2.    Viscous forcing for moisture
!	       ----------------------------
!
IF (OVISC_R .AND. (SIZE(PRT,1) > 0)) THEN
!
!
     DO IK = 1, KRR
        PRRS(:,:,:,IK) = PRRS(:,:,:,IK) + PNU/PPRANDTL * &
             LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PRT(:,:,:,IK))
     END DO
!
!
END IF
!
IF (LBUDGET_RV) CALL BUDGET (PRRS(:,:,:,1),6,'VISC_BU_RRV')
IF (LBUDGET_RC) CALL BUDGET (PRRS(:,:,:,2),7,'VISC_BU_RRC')
IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:,3),8,'VISC_BU_RRR')
IF (LBUDGET_RI) CALL BUDGET (PRRS(:,:,:,4),9,'VISC_BU_RRI')
IF (LBUDGET_RS) CALL BUDGET (PRRS(:,:,:,5),10,'VISC_BU_RRS')
IF (LBUDGET_RG) CALL BUDGET (PRRS(:,:,:,6),11,'VISC_BU_RRG')
IF (LBUDGET_RH) CALL BUDGET (PRRS(:,:,:,7),12,'VISC_BU_RRH')
!
!-------------------------------------------------------------------------------
!
!*       3.    Viscous forcing for passive scalars
!	       -----------------------------------
!
IF (OVISC_SV .AND. (SIZE(PSVT,1) > 0)) THEN
!
!
      DO IK = 1, KSV
         PRSVS(:,:,:,IK) = PRSVS(:,:,:,IK) + PNU/PPRANDTL * &
              LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,PSVT(:,:,:,IK))
      END DO
!
END IF
!
IF (LBUDGET_SV) THEN
  DO  IK = 1, KSV
    CALL BUDGET (PRSVS(:,:,:,IK), 12+IK, 'VISC_BU_RSV')
  END DO
END IF
!

!-------------------------------------------------------------------------------
!
!*       4.    Viscous forcing for momentum
!	       ----------------------------
!
IF (OVISC_UVW) THEN
!
!*       4.1   U - component
!
!
      ZY1 = MXF(PUT)
      IF (ODRAG) THEN
         ZY1(:,:,1) = PDRAG * ZY1(:,:,2)
      ENDIF
!
! 
       ZLAPu = LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,   &
                   PDZY,PDZZ,PRHODJ,ZY1)
!! Update halo to compute the source term
 NULLIFY(TZFIELDS_ll)
 CALL ADD3DFIELD_ll(TZFIELDS_ll,ZLAPu)
 CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
 CALL CLEANLIST_ll(TZFIELDS_ll)
!
 PRUS = PRUS + MXM(PNU*ZLAPu)
!
!*       4.2   V - component
!              -------------
!
  IF (.NOT. L2D) THEN

      ZY2 = MYF(PVT) 
      IF (ODRAG) THEN
        ZY2(:,:,1) = PDRAG * ZY2(:,:,2)
      ENDIF
!
      ZLAPv =  LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,   &
                     PDZY,PDZZ,PRHODJ,ZY2)
!! Update halo to compute the source term
!
 NULLIFY(TZFIELDS_ll)
 CALL ADD3DFIELD_ll(TZFIELDS_ll,ZLAPv)
 CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
 CALL CLEANLIST_ll(TZFIELDS_ll)
!
 PRVS = PRVS + MYM(PNU*ZLAPv)

ENDIF 

!
!*       4.3   W - component
!              -------------
!
   IKB = JPVEXT + 1
   IKE = SIZE(PWT,3) - JPVEXT

   ZTMP = MZF(1,IKU,1,PWT)
!
   IF (ODRAG) THEN
         WHERE (PDRAG==-1)
         ZTMP(:,:,IKB) = 0.
         ENDWHERE
   ENDIF
!
   DO IK = 1,JPVEXT
      ZTMP(:,:,IK) = ZTMP(:,:,IKB)
      ZTMP(:,:,IKE+IK) = ZTMP(:,:,IKE)
   END DO
!
   ZTMP = MZM(1,IKU,1, PNU * &
          LAP_M(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,PRHODJ,ZTMP) )
!
   DO IK = 1,JPVEXT
      ZTMP(:,:,IK) = ZTMP(:,:,IKB)
      ZTMP(:,:,IKE+IK) = ZTMP(:,:,IKE) 
   END DO
   PRWS = PRWS + ZTMP
!
!!! Debug provisoire dans le cas ou le noslip est applique jusqu'au bord de
!sortie de flux en OPEN
  IF ( LWEST_ll().AND.(ODRAG)) THEN
    IF ( MINVAL(PDRAG(IIU,:))== -1) THEN
              DO JK=1,IKU
                WHERE(PDRAG(IIU,:)== -1)
            PRUS(IIU,:,JK) = PRUS(IIU-1,:,JK)
            PRVS(IIU,:,JK) = PRVS(IIU-1,:,JK)
            PRWS(IIU,:,JK) = PRWS(IIU-1,:,JK)
                ENDWHERE
            END DO
   ENDIF
  ENDIF
END IF
!
IF (LBUDGET_U) CALL BUDGET (PRUS,1,'VISC_BU_RU')
IF (LBUDGET_V) CALL BUDGET (PRVS,2,'VISC_BU_RV')
IF (LBUDGET_W) CALL BUDGET (PRWS,2,'VISC_BU_RW')
!
END SUBROUTINE VISCOSITY
