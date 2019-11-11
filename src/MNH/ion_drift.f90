
!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #####################
      MODULE MODI_ION_DRIFT
!     #####################

INTERFACE
!
      SUBROUTINE ION_DRIFT(PDRIFTP, PDRIFTM, PSVT, HLBCX, HLBCY)
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)    :: HLBCX,HLBCY
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PDRIFTP, PDRIFTM
REAL, DIMENSION(:,:,:,:),       INTENT(INOUT) :: PSVT
!
END SUBROUTINE ION_DRIFT
END INTERFACE
END MODULE MODI_ION_DRIFT
!
!     ##########################################################
      SUBROUTINE ION_DRIFT(PDRIFTP, PDRIFTM, PSVT, HLBCX, HLBCY)
!     ##########################################################
!
!!    PURPOSE
!!    -------
!!
!!    Compute the ion drift defined as
!!        -/+ Div ( mu(+/-) N(+/-) E )
!!        where mu is the ion mobility (m2/V.s), N is the ion mixing ratio (1/kg)
!!              and E is the electric field (V/m)
!!           +/- mu(+/-) E  = drift motion
!!
!!    AUTHOR
!!    ------
!!          M. Chong      01/2010
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!         
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_METRICS_n, ONLY : XDXX, XDYY, XDZX, XDZY, XDZZ
USE MODD_NSV,  ONLY: NSV_ELECBEG, NSV_ELECEND, XSVMIN
USE MODD_ELEC_n, ONLY : XCION_POS_FW, XCION_NEG_FW, &
                        XMOBIL_POS, XMOBIL_NEG,     &
                        XEFIELDU, XEFIELDV, XEFIELDW
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
USE MODE_ll
USE MODE_ELEC_ll
!
USE MODI_SHUMAN
USE MODI_GDIV
USE MODI_ION_BOUND4DRIFT
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)    :: HLBCX,HLBCY
REAL, DIMENSION(:,:,:),         INTENT(INOUT) :: PDRIFTP, PDRIFTM
REAL, DIMENSION(:,:,:,:),       INTENT(INOUT) :: PSVT
!
!
!*       0.2   declarations of local variables
!
INTEGER :: IIB, IIE  ! index of first and last inner mass points along x
INTEGER :: IJB, IJE  ! index of first and last inner mass points along y
INTEGER :: IKB, IKE  ! index of first and last inner mass points along z
INTEGER :: IKU
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZDRIFTX
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZDRIFTY
REAL, DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZDRIFTZ
!
INTEGER :: IINFO_ll      ! return code of parallel routine
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!
NULLIFY(TZFIELDS_ll)
!
!
!------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES
!              ----------------------------------------------
!
! beginning and end indexes of the physical subdomain
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PSVT,3) - JPVEXT
IKU = SIZE(PSVT,3)
!
!
!-------------------------------------------------------------------------------
!
!*      2.     UPDATE BOUNDARY CONDITION FOR IONS ACCORDING TO THE DRIFT MOTION
!              ----------------------------------------------------------------
!
IF (LWEST_ll() ) THEN
  XEFIELDU(IIB-1,:,:) = XEFIELDU(IIB,:,:)
ENDIF
!
IF (LEAST_ll() ) THEN
  XEFIELDU(IIE+1,:,:) = XEFIELDU(IIE,:,:)
ENDIF
!
IF (LSOUTH_ll()) THEN
  XEFIELDV(:,IJB-1,:) = XEFIELDV(:,IJB,:)
ENDIF   
!
IF (LNORTH_ll()) THEN
  XEFIELDV(:,IJE+1,:) = XEFIELDV(:,IJE,:)
ENDIF   
!
CALL ADD3DFIELD_ll(TZFIELDS_ll,XEFIELDU)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XEFIELDV)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
!  specify lateral boundary ion mixing ratio
CALL ION_BOUND4DRIFT (HLBCX,HLBCY,XEFIELDU,XEFIELDV,PSVT)
!
CALL ADD3DFIELD_ll(TZFIELDS_ll,PSVT(:,:,:,NSV_ELECBEG))
CALL ADD3DFIELD_ll(TZFIELDS_ll,PSVT(:,:,:,NSV_ELECEND))
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
! specify upper boundary ion mixing ratio
WHERE (XEFIELDW(:,:,IKE+1) .GE. 0.)   ! Out(In)flow for positive (negative) ions 
  PSVT (:,:,IKE+1,NSV_ELECBEG) = MAX(2. * PSVT (:,:,IKE,NSV_ELECBEG) -  &
                                  PSVT (:,:,IKE-1,NSV_ELECBEG),XSVMIN(NSV_ELECBEG))
  PSVT (:,:,IKE+1,NSV_ELECEND) = XCION_NEG_FW(:,:,IKE+1)
ELSE WHERE      ! In(Out)flow for positive (negative) ions
  PSVT (:,:,IKE+1,NSV_ELECBEG) = XCION_POS_FW(:,:,IKE+1)
  PSVT (:,:,IKE+1,NSV_ELECEND) = MAX(2.* PSVT (:,:,IKE,NSV_ELECEND) -  &
                                  PSVT (:,:,IKE-1,NSV_ELECEND),XSVMIN(NSV_ELECEND))
END WHERE  
!
XEFIELDW(:,:,IKB-1) = XEFIELDW(:,:,IKB)
XEFIELDW(:,:,IKE+1) = XEFIELDW(:,:,IKE)
!
!-------------------------------------------------------------------------------
!
!*      3.     DRIFT MOTION TAKEN AS THE DIVERGENCE OF THE DRIFT FLUXES
!              --------------------------------------------------------
!
!*      3.1  positive ion source (drifting along E)
!
! x-component of div term
ZDRIFTX(:,:,:) = PSVT(:,:,:,NSV_ELECBEG) * XMOBIL_POS(:,:,:)
ZDRIFTX(:,:,:) = ZDRIFTX(:,:,:) * XEFIELDU(:,:,:)
ZDRIFTX(:,:,:) = -MXM(ZDRIFTX(:,:,:)) ! Put components at flux sides
!
! y-component of div term
ZDRIFTY(:,:,:) = PSVT(:,:,:,NSV_ELECBEG) * XMOBIL_POS(:,:,:)
ZDRIFTY(:,:,:) = ZDRIFTY(:,:,:) * XEFIELDV(:,:,:)
ZDRIFTY(:,:,:) = -MYM(ZDRIFTY(:,:,:)) ! Put components at flux sides
!
! z-component of div term
ZDRIFTZ(:,:,:) = PSVT(:,:,:,NSV_ELECBEG) * XMOBIL_POS(:,:,:)
ZDRIFTZ(:,:,:) = ZDRIFTZ(:,:,:) * XEFIELDW(:,:,:)
ZDRIFTZ(:,:,:) = -MZM(1,IKU,1,ZDRIFTZ(:,:,:)) ! Put components at flux sides
!
IF (LWEST_ll( ))  ZDRIFTX(IIB-1,:,:) = ZDRIFTX(IIB,:,:)
IF (LEAST_ll( ))  ZDRIFTX(IIE+1,:,:) = ZDRIFTX(IIE,:,:)
IF (LSOUTH_ll( )) ZDRIFTY(:,IJB-1,:) = ZDRIFTY(:,IJB,:)
IF (LNORTH_ll( )) ZDRIFTY(:,IJE+1,:) = ZDRIFTY(:,IJE,:)
ZDRIFTZ(:,:,IKB-1) = ZDRIFTZ(:,:,IKB)
ZDRIFTZ(:,:,IKE+1) = ZDRIFTZ(:,:,IKE)
!
CALL GDIV(HLBCX,HLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ,ZDRIFTX,ZDRIFTY,ZDRIFTZ,PDRIFTP)
!
!*      3.2    negative ion source (drifting counter E)
!
! x-component of div term
ZDRIFTX(:,:,:) = PSVT(:,:,:,NSV_ELECEND) * XMOBIL_NEG(:,:,:)
ZDRIFTX(:,:,:) = ZDRIFTX(:,:,:) * XEFIELDU(:,:,:)
ZDRIFTX(:,:,:) = +MXM(ZDRIFTX(:,:,:)) ! Put components at flux sides
!
! y-component of div term
ZDRIFTY(:,:,:) = PSVT(:,:,:,NSV_ELECEND) * XMOBIL_NEG(:,:,:)
ZDRIFTY(:,:,:) = ZDRIFTY(:,:,:) * XEFIELDV(:,:,:)
ZDRIFTY(:,:,:) = +MYM(ZDRIFTY(:,:,:)) ! Put components at flux sides
!
! z-component of div term
ZDRIFTZ(:,:,:) = PSVT(:,:,:,NSV_ELECEND) * XMOBIL_NEG(:,:,:)
ZDRIFTZ(:,:,:) = ZDRIFTZ(:,:,:) * XEFIELDW(:,:,:)
ZDRIFTZ(:,:,:) = +MZM(1,IKU,1,ZDRIFTZ(:,:,:)) ! Put components at flux sides

!
IF (LWEST_ll( ))  ZDRIFTX(IIB-1,:,:) = ZDRIFTX(IIB,:,:)
IF (LEAST_ll( ))  ZDRIFTX(IIE+1,:,:) = ZDRIFTX(IIE,:,:)
IF (LSOUTH_ll( )) ZDRIFTY(:,IJB-1,:) = ZDRIFTY(:,IJB,:)
IF (LNORTH_ll( )) ZDRIFTY(:,IJE+1,:) = ZDRIFTY(:,IJE,:)
ZDRIFTZ(:,:,IKB-1) = ZDRIFTZ(:,:,IKB)
ZDRIFTZ(:,:,IKE+1) = ZDRIFTZ(:,:,IKE)
!
CALL GDIV(HLBCX,HLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ,ZDRIFTX,ZDRIFTY,ZDRIFTZ,PDRIFTM)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ION_DRIFT
