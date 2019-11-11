!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ######################################
     MODULE MODI_COMPUTE_FUNCTION_THERMO
!    ######################################
!
INTERFACE
      
!     #################################################################
      SUBROUTINE COMPUTE_FUNCTION_THERMO( KRR,                            &
                                       PTH, PR, PEXN, PPABS,                 &
                                       PT,PLVOCPEXN,PLSOCPEXN,               &
                                       PAMOIST,PATHETA                       )
!     #################################################################

!*               1.1  Declaration of Arguments
!

INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.

REAL, DIMENSION(:,:,:), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:,:)  , INTENT(IN) :: PPABS,PEXN    ! pressure, Exner funct.

REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PT      ! temperature
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PLVOCPEXN,PLSOCPEXN     ! L/(cp*Pi)

REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL  ::  PAMOIST,PATHETA
!
END SUBROUTINE COMPUTE_FUNCTION_THERMO

END INTERFACE
!
END MODULE MODI_COMPUTE_FUNCTION_THERMO

                                       
!     #################################################################
      SUBROUTINE COMPUTE_FUNCTION_THERMO( KRR,                            &
                                       PTH, PR, PEXN, PPABS,                 &
                                       PT,PLVOCPEXN,PLSOCPEXN,               &
                                       PAMOIST,PATHETA                       )
!     #################################################################
!
!!
!!****  *COMPUTE_FUNCTION_THERMO* - 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!     
!!     JP Pinty      *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!     Original   24/02/03
!!     Externalisation of computations done in TURB and MF_TURB (Malardel and Pergaud, fev. 2007)
!!     Optimization : V.Masson, 09/2010
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.

REAL, DIMENSION(:,:,:), INTENT(IN)   :: PTH      ! theta
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PR       ! water species
REAL, DIMENSION(:,:,:)  , INTENT(IN) :: PPABS,PEXN    ! pressure, Exner funct.

REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PT      ! temperature
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PLVOCPEXN,PLSOCPEXN     ! L/(cp*Pi)

REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL  ::  PAMOIST,PATHETA
!
!-------------------------------------------------------------------------------
! 
!*       0.2   Declarations of local variables
!
REAL                :: ZEPS         ! XMV / XMD
REAL, DIMENSION(SIZE(PTH,1),SIZE(PTH,2),SIZE(PTH,3)) ::     &
          ZCP,                        &  ! Cp 
          ZE,                         &  ! Saturation mixing ratio
          ZDEDT,                      &  ! Saturation mixing ratio derivative
          ZFRAC_ICE,                  &  ! Ice fraction
          ZAMOIST_W,                  &  ! Coefficients for s = f (Thetal,Rnp)
          ZATHETA_W,                  &  !
          ZAMOIST_I,                  &  !
          ZATHETA_I                      !

INTEGER             :: JRR
INTEGER             :: IRRL, IRRI
!
!-------------------------------------------------------------------------------
!
  IRRI=0
  IRRL=0
  IF (KRR>=2) IRRL=MIN(KRR-1,2)
  IF (KRR>=4) IRRI=KRR-3
!
  ZEPS = XMV / XMD

  PLVOCPEXN(:,:,:) = 0.
  PLSOCPEXN(:,:,:) = 0.
  ZFRAC_ICE(:,:,:) = 0.0

!
!*       Cph
!
ZCP=XCPD

IF (KRR > 0) ZCP(:,:,:) = ZCP(:,:,:) + XCPV * PR(:,:,:,1)

DO JRR = 2,1+IRRL  ! loop on the liquid components  
   ZCP(:,:,:)  = ZCP(:,:,:) + XCL * PR(:,:,:,JRR)
END DO

DO JRR = 2+IRRL,1+IRRL+IRRI ! loop on the solid components   
  ZCP(:,:,:)  = ZCP(:,:,:)  + XCI * PR(:,:,:,JRR)
END DO

!*      Temperature
!
PT(:,:,:) =  PTH(:,:,:) * PEXN(:,:,:)
!
!
!! Liquid water
!
IF ( IRRL >= 1 ) THEN 
!
!*       Lv/Cph 
!
  PLVOCPEXN(:,:,:) = (XLVTT + (XCPV-XCL) *  (PT(:,:,:)-XTT) ) / ZCP(:,:,:)
!
 IF (PRESENT(PAMOIST)) THEN
!
!*      Saturation vapor pressure with respect to water
!
  ZE(:,:,:) =  EXP( XALPW - XBETAW/PT(:,:,:) - XGAMW*ALOG( PT(:,:,:) ) )
!
!*      Saturation  mixing ratio with respect to water
!
  ZE(:,:,:) =  ZE(:,:,:) * ZEPS / ( PPABS(:,:,:) - ZE(:,:,:) )
!
!*      Compute the saturation mixing ratio derivative (rvs')
!
  ZDEDT(:,:,:) = ( XBETAW / PT(:,:,:)  - XGAMW ) / PT(:,:,:)   &
                 * ZE(:,:,:) * ( 1. + ZE(:,:,:) / ZEPS )
!
!*      Compute Amoist
!
  PAMOIST(:,:,:) = 0.
  PATHETA(:,:,:) = 0.
  ZAMOIST_W(:,:,:)=  0.5 / ( 1.0 + ZDEDT(:,:,:) * PLVOCPEXN(:,:,:) )
!
!*      Compute Atheta
!
  ZATHETA_W(:,:,:)= ZAMOIST_W(:,:,:) * PEXN(:,:,:) *                          &
        ( ( ZE(:,:,:) - PR(:,:,:,1) ) * PLVOCPEXN(:,:,:) /                    &
          ( 1. + ZDEDT(:,:,:) * PLVOCPEXN(:,:,:) )           *              &
          (                                                             &
           ZE(:,:,:) * (1. + ZE(:,:,:)/ZEPS)                                &
                        * ( -2.*XBETAW/PT(:,:,:) + XGAMW ) / PT(:,:,:)**2   &
          +ZDEDT(:,:,:) * (1. + 2. * ZE(:,:,:)/ZEPS)                        &
                        * ( XBETAW/PT(:,:,:) - XGAMW ) / PT(:,:,:)          &
          )                                                             &
         - ZDEDT(:,:,:)                                                   &
        )
 END IF
!
!! Solid water
!
  IF ( IRRI >= 1 ) THEN 

!
!*       Fraction of ice
!
    WHERE(PR(:,:,:,2)+PR(:,:,:,4)>0.0)
      ZFRAC_ICE(:,:,:) = PR(:,:,:,4) / ( PR(:,:,:,2)+PR(:,:,:,4) )
    END WHERE
!
!*       Ls/Cph 
!
    PLSOCPEXN(:,:,:) = (XLSTT + (XCPV-XCI) *  (PT(:,:,:)-XTT) ) / ZCP(:,:,:)
!
   IF (PRESENT(PAMOIST)) THEN
!
!*      Saturation vapor pressure with respect to ice
!
    ZE(:,:,:) =  EXP( XALPI - XBETAI/PT(:,:,:) - XGAMI*ALOG( PT(:,:,:) ) )
!
!*      Saturation  mixing ratio with respect to ice
!
    ZE(:,:,:) =  ZE(:,:,:) * ZEPS / ( PPABS(:,:,:) - ZE(:,:,:) )
!
!*      Compute the saturation mixing ratio derivative (rvs')
!
    ZDEDT(:,:,:) = ( XBETAI / PT(:,:,:)  - XGAMI ) / PT(:,:,:)   &
                   * ZE(:,:,:) * ( 1. + ZE(:,:,:) / ZEPS )
!
!*      Compute Amoist
!
    ZAMOIST_I(:,:,:)=  0.5 / ( 1.0 + ZDEDT(:,:,:) * PLSOCPEXN(:,:,:) )
!
!*      Compute Atheta
!
    ZATHETA_I(:,:,:)= ZAMOIST_I(:,:,:) * PEXN(:,:,:) *                        &
        ( ( ZE(:,:,:) - PR(:,:,:,1) ) * PLSOCPEXN(:,:,:) /                    &
          ( 1. + ZDEDT(:,:,:) * PLSOCPEXN(:,:,:) )           *              &
          (                                                             &
           ZE(:,:,:) * (1. + ZE(:,:,:)/ZEPS)                                &
                        * ( -2.*XBETAI/PT(:,:,:) + XGAMI ) / PT(:,:,:)**2   &
          +ZDEDT(:,:,:) * (1. + 2. * ZE(:,:,:)/ZEPS)                        &
                        * ( XBETAI/PT(:,:,:) - XGAMI ) / PT(:,:,:)          &
          )                                                             &
         - ZDEDT(:,:,:)                                                   &
        )


   ENDIF
  ENDIF

 IF (PRESENT(PAMOIST)) THEN
  PAMOIST(:,:,:) = (1.0-ZFRAC_ICE(:,:,:))*ZAMOIST_W(:,:,:) &
                         +ZFRAC_ICE(:,:,:) *ZAMOIST_I(:,:,:)
  PATHETA(:,:,:) = (1.0-ZFRAC_ICE(:,:,:))*ZATHETA_W(:,:,:) &
                         +ZFRAC_ICE(:,:,:) *ZATHETA_I(:,:,:)
 ENDIF

!
!*      Lv/Cph/Exner and  Ls/Cph/Exner
!
  PLVOCPEXN(:,:,:) = PLVOCPEXN(:,:,:) / PEXN(:,:,:)
  PLSOCPEXN(:,:,:) = PLSOCPEXN(:,:,:) / PEXN(:,:,:)
!
ENDIF

END SUBROUTINE COMPUTE_FUNCTION_THERMO
