!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##########################
      MODULE MODI_INI_FIELD_ELEC
!     ##########################
!
INTERFACE
!
      SUBROUTINE INI_FIELD_ELEC (PDXX, PDYY, PDZZ, PDZX, PDZY, PZZ)
!
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDXX     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDYY     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDZZ     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDZX     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDZY     ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PZZ      ! vertical grid
!
END SUBROUTINE INI_FIELD_ELEC
END INTERFACE
END MODULE MODI_INI_FIELD_ELEC
!
!     ############################################################
      SUBROUTINE INI_FIELD_ELEC(PDXX, PDYY, PDZZ, PDZX, PDZY, PZZ)
!     ############################################################
!
!
!!****  *INI_FIELD_ELEC* - routine to initialize the electric field 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize the variables
!     of the electric field computation 
!
!!**  METHOD
!!    ------
!!      The initialization of the scheme is performed as follows :
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty    * Laboratoire d'Aérologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    29/11/02
!!      C. Barthe   10/11/09   phasage en version 4.8.1
!!      M. Chong    26/01/10   Small ions parameters 
!!                           + Fair weather field from Helsdon-Farley
!!                             (JGR, 1987, 5661-5675)
!!      J.-P. Pinty 01/07/12   Add a non-homogeneous Neuman fair-weather 
!!                             boundary condition at the top
!!
!!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_CST
USE MODD_REF_n, ONLY: XRHODREF
USE MODD_ELEC_DESCR
USE MODD_ELEC_n        
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
USE MODI_GDIV
USE MODI_SHUMAN
!
USE MODE_ll
USE MODE_FM
!
IMPLICIT NONE
!
!*	0.1	Declaration of dummy arguments
!
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDXX  ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDYY  ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDZZ  ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDZX  ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PDZY  ! Metric coefficients
REAL, DIMENSION(:,:,:),  INTENT(IN) ::  PZZ   ! vertical grid
!
!*	0.2	Declaration of local variables
!
! 
CHARACTER(LEN=4), DIMENSION(2) :: ZLBCX  ! x-direction LBC type 
CHARACTER(LEN=4), DIMENSION(2) :: ZLBCY  ! y-direction LBC type 
!
INTEGER :: JK     ! loop over the vertical levels
INTEGER :: IINFO_ll ! 
INTEGER :: IKB,IKE,IKU      ! Indices for the first and last point along vertical
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZZMASS, ZWORK, ZWORK1, ZWORK2
!
TYPE(LIST_ll),POINTER  :: TZFIELDS_ll ! list of fields to exchange
!
!
!------------------------------------------------------------------------------
!
!*      1.    INITIALIZATIONS
!             ---------------
! 
IKB = 1 + JPVEXT
IKE = SIZE(PZZ,3) - JPVEXT
IKU = SIZE(PZZ,3)
ZLBCX = 'OPEN'  ! forced LBC
ZLBCY = 'OPEN'  ! forced LBC
!
NULLIFY(TZFIELDS_ll)
!
! Allocations
!
ALLOCATE( XEFIELDU(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( XEFIELDV(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( XEFIELDW(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( XESOURCEFW(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( ZZMASS(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )  !Alt at mass point
ALLOCATE( ZWORK(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( ZWORK1(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
IF( .NOT. LCOSMIC_APPROX ) THEN
  ALLOCATE( ZWORK2(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
END IF
ALLOCATE( XIONSOURCEFW(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( XCION_POS_FW(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( XCION_NEG_FW(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( XMOBIL_POS(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
ALLOCATE( XMOBIL_NEG(SIZE(PDZZ,1),SIZE(PDZZ,2),SIZE(PDZZ,3)) )
!
!++ jpp
ALLOCATE( XEPOTFW_TOP(SIZE(PDZZ,1),SIZE(PDZZ,2)) )
!-- jpp
!
!------------------------------------------------------------------------------
!
!*       2.    FAIR WEATHER ELECTRIC FIELD
!              ---------------------------  
!
! The vertical component of the electric field is given by : E=E_0 * exp(k_e*z) 
! where E_0 = -100 V m^-1 and k_e = -292e-6 m^-1 
! We define the electric field as : E = - rhodJ * Nabla V
!
!  Helsdon-Farley: E=E_0 (b1 exp(-a1 z) + b2 exp(-a2 z) + b3 exp(-a3 z)
XEFIELDU(:,:,:) = 0.
XEFIELDV(:,:,:) = 0. 
!
! Initialization of Fair Weather Electric Field  at W-point
IF( .NOT. LFW_HELFA ) THEN
  XEFIELDW(:,:,:) = XE_0 * EXP(XKEF * PZZ(:,:,:))
ELSE
  XEFIELDW(:,:,:) = XE0_HF * (XB1_HF*EXP(-XA1_HF*PZZ(:,:,:)) &
                            + XB2_HF*EXP(-XA2_HF*PZZ(:,:,:)) &
                            + XB3_HF*EXP(-XA3_HF*PZZ(:,:,:)))
END IF
!++ jpp
XEPOTFW_TOP(:,:) = XEFIELDW(:,:,SIZE(PDZZ,3)) ! used in the top boundary
                                              ! condition when inverting the
                                              ! Gauss equation in V (here EPOT)
!-- jpp
XEFIELDW(:,:,SIZE(PDZZ,3)) = 2. * XEFIELDW(:,:,SIZE(PDZZ,3)-1) -  &
                                  XEFIELDW(:,:,SIZE(PDZZ,3)-2)

! Computing the mobility of small positive (negative) ions at Mass-point
ZZMASS = MZF(1,IKU,1, PZZ )   ! altitude at mass point

DO JK = 2,SIZE(PZZ,3)-1
  XMOBIL_POS(:,:,JK) = XF_POS * EXP( XEXPMOB* ZZMASS(:,:,JK) )
  XMOBIL_NEG(:,:,JK) = XF_NEG * EXP( XEXPMOB* ZZMASS(:,:,JK) )
END DO

XMOBIL_POS(:,:,1) = 2.0 * XMOBIL_POS(:,:,2) - XMOBIL_POS(:,:,3)
XMOBIL_POS(:,:,SIZE(PDZZ,3)) = 2. * XMOBIL_POS(:,:,SIZE(PDZZ,3)-1) - &
                                    XMOBIL_POS(:,:,SIZE(PDZZ,3)-2)
XMOBIL_NEG(:,:,1) = 2.0*XMOBIL_NEG(:,:,2) - XMOBIL_NEG(:,:,3)
XMOBIL_NEG(:,:,SIZE(PDZZ,3)) = 2. * XMOBIL_NEG(:,:,SIZE(PDZZ,3)-1) - &
                                     XMOBIL_NEG(:,:,SIZE(PDZZ,3)-2)
!
! Initial number concentrations of small positive (negative) free ions
IF( .NOT. LFW_HELFA ) THEN
  ZWORK(:,:,:) = XE_0 * EXP(XKEF * ZZMASS(:,:,:))
  ZWORK1(:,:,:) = XE_0 * XKEF * EXP(XKEF * ZZMASS(:,:,:))
  IF(.NOT. LCOSMIC_APPROX) THEN
    ZWORK2(:,:,:) = XE_0 * XKEF * XKEF * EXP(XKEF * ZZMASS(:,:,:))
  END IF
ELSE
  ZWORK(:,:,:)= XE0_HF * (XB1_HF*EXP(-XA1_HF*ZZMASS(:,:,:)) &
                        + XB2_HF*EXP(-XA2_HF*ZZMASS(:,:,:)) &
                        + XB3_HF*EXP(-XA3_HF*ZZMASS(:,:,:)))
  ZWORK1(:,:,:)= XE0_HF * (-XB1_HF*XA1_HF*EXP(-XA1_HF*ZZMASS(:,:,:)) &
                           -XB2_HF*XA2_HF*EXP(-XA2_HF*ZZMASS(:,:,:)) &
                           -XB3_HF*XA3_HF*EXP(-XA3_HF*ZZMASS(:,:,:)))
  IF(.NOT. LCOSMIC_APPROX) THEN
    ZWORK2(:,:,:)= XE0_HF * (XB1_HF*XA1_HF*XA1_HF*EXP(-XA1_HF*ZZMASS(:,:,:)) &
                            +XB2_HF*XA2_HF*XA2_HF*EXP(-XA2_HF*ZZMASS(:,:,:)) &
                            +XB3_HF*XA3_HF*XA3_HF*EXP(-XA3_HF*ZZMASS(:,:,:)))
  END IF
END IF
!
XCION_POS_FW(:,:,:) = (XMOBIL_NEG(:,:,:) * XEPSILON * ZWORK1(:,:,:) + &
                       XJCURR_FW / ZWORK(:,:,:)) /                     &
                      (XECHARGE * (XMOBIL_POS(:,:,:) + XMOBIL_NEG(:,:,:)))
XCION_NEG_FW(:,:,:) = XCION_POS_FW - XEPSILON * ZWORK1(:,:,:) / XECHARGE
XCION_POS_FW(:,:,SIZE(PDZZ,3)) = 2. * XCION_POS_FW(:,:,SIZE(PDZZ,3)-1) - &
                                      XCION_POS_FW(:,:,SIZE(PDZZ,3)-2)
XCION_NEG_FW(:,:,SIZE(PDZZ,3)) = 2. * XCION_NEG_FW(:,:,SIZE(PDZZ,3)-1) - &
                                       XCION_NEG_FW(:,:,SIZE(PDZZ,3)-2)
!
WHERE(XCION_NEG_FW < 0.) XCION_NEG_FW = 0.
!
! Computing the ion source from cosmic rays
XIONSOURCEFW(:,:,:) = XIONCOMB * XCION_POS_FW(:,:,:) * XCION_NEG_FW(:,:,:)
!
IF ( .NOT. LCOSMIC_APPROX ) THEN
  XIONSOURCEFW(:,:,:) = XIONSOURCEFW(:,:,:) +                               &
                        XMOBIL_POS(:,:,:) * XMOBIL_NEG(:,:,:) * XEPSILON * &
                       (XEXPMOB * ZWORK(:,:,:) * ZWORK1(:,:,:) +            &
                        ZWORK1(:,:,:) * ZWORK1(:,:,:)        +              &
                        ZWORK(:,:,:) * ZWORK2(:,:,:)) /                     &
                       (XECHARGE * (XMOBIL_POS(:,:,:) + XMOBIL_NEG(:,:,:)))

  XIONSOURCEFW(:,:,1) = 0.
  XIONSOURCEFW(:,:,SIZE(PDZZ,3)) = 2. * XIONSOURCEFW(:,:,SIZE(PDZZ,3)-1) -  &
                                        XIONSOURCEFW(:,:,SIZE(PDZZ,3)-2)
END IF
!
!  Transform ion concentration into ion mixing ratio (Number/kg of air)

XCION_POS_FW(:,:,:)  = XCION_POS_FW(:,:,:)  / XRHODREF(:,:,:)
XCION_NEG_FW(:,:,:) = XCION_NEG_FW(:,:,:) / XRHODREF(:,:,:)
XCION_POS_FW(:,:,SIZE(PDZZ,3)) = 2. * XCION_POS_FW(:,:,SIZE(PDZZ,3)-1) - &
                                      XCION_POS_FW(:,:,SIZE(PDZZ,3)-2)
XCION_NEG_FW(:,:,SIZE(PDZZ,3)) = 2. * XCION_NEG_FW(:,:,SIZE(PDZZ,3)-1) - &
                                       XCION_NEG_FW(:,:,SIZE(PDZZ,3)-2)
!
XEFIELDW(:,:,1) = 0.            ! Electric field null in a conductor
XEFIELDW(:,:,SIZE(PDZZ,3)) = 0. ! either the ground or the ionosphere!
!
!
!------------------------------------------------------------------------------
!
!*       3.    FAIR WEATHER SPACE CHARGE
!              -------------------------
!
CALL ADD3DFIELD_ll(TZFIELDS_ll,XEFIELDU)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XEFIELDV)
CALL ADD3DFIELD_ll(TZFIELDS_ll,XEFIELDW)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
CALL GDIV (ZLBCX, ZLBCY,                 &
           PDXX, PDYY, PDZX, PDZY, PDZZ, &
           XEFIELDU,XEFIELDV,XEFIELDW,   &
           XESOURCEFW                    )
!
XESOURCEFW(:,:,:) = XESOURCEFW(:,:,:) * XEPSILON  ! Nabla E * epsilon  = + rho
                                                  ! C / m^3
XESOURCEFW(:,:,SIZE(PZZ,3)) = XESOURCEFW(:,:,SIZE(PZZ,3)-1)
!
DEALLOCATE(ZZMASS)
DEALLOCATE(ZWORK)
DEALLOCATE(ZWORK1)
IF( .NOT. LCOSMIC_APPROX) THEN
  DEALLOCATE(ZWORK2)
END IF
!
!
!------------------------------------------------------------------------------
!
END SUBROUTINE INI_FIELD_ELEC
