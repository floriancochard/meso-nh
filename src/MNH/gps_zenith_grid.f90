!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!##########################################
MODULE MODI_GPS_ZENITH_GRID
!##########################################
INTERFACE
      SUBROUTINE GPS_ZENITH_GRID(PRT,PTEMP,PPABST,PZTD,PZHD,PZWD)
!
!*       0.1   Declarations of dummy arguments :
!
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PRT      ! Water vapor mixing ratio
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PTEMP    ! Air temperature
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PPABST   ! Absolute pressure
!
REAL,  DIMENSION(:,:),   INTENT(OUT)     :: PZTD    ! Zenithal Total Delay
REAL,  DIMENSION(:,:),   INTENT(OUT)     :: PZHD    ! Zenithal Hydrostatic Delay
REAL,  DIMENSION(:,:),   INTENT(OUT)     :: PZWD    ! Zenithal Wet Delay
!
END SUBROUTINE GPS_ZENITH_GRID
END INTERFACE
END MODULE MODI_GPS_ZENITH_GRID
!
!     ###################################################################
      SUBROUTINE GPS_ZENITH_GRID(PRT,PTEMP,PPABST,PZTD,PZHD,PZWD)
!     ###################################################################
!
!!****  *GPS_ZENITH_GRID * - computes  GPS zenithal delays parameters
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the synthetic ZTD, ZHD, ZWD
!!
!!**  METHOD
!!    ------
!!      The delays are computed using the refractivity formula and its different
!!    approximation. The delay can be separated in two parts.
!!    The ZTD is the sum of the ZHD, including all component of the atmosphere, and 
!!    the ZWD, including the contribution of the non-dipolar moment of the molecule
!!    of water vapor. 
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!      Brenot et al, 2006, JGR
!!
!!    AUTHORS
!!    -------
!!      H. Brenot      * Laboratoire de Geophysique Interne et Tectonophysique *
!!                                   /  * Meteo France *      
!!       &  V. Ducrocq    * Meteo France *    
!!       &  A. Walpersdorf    * Laboratoire de Geophysique Interne et Tectonophysique *    
!!     
!!    MODIFICATIONS
!!    -------------
!!      Original    18/11/04
!!      Modified    4/12/2007
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_GR_FIELD_n
USE MODD_DIAG_FLAG
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_GRID_n
USE MODE_GRIDPROJ
USE MODE_ll
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PRT        ! Water vapour mixing ratio
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PTEMP      ! Air temperature
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PPABST     ! Absolute pressure
!
REAL,  DIMENSION(:,:),   INTENT(OUT)        :: PZTD    ! Zenithal Total Delay
REAL,  DIMENSION(:,:),   INTENT(OUT)        :: PZHD    ! Zenithal Hydrostatic Delay
REAL,  DIMENSION(:,:),   INTENT(OUT)        :: PZWD    ! Zenithal Wet Delay
!
!*       0.2   Declarations of local variables :
!
!-------- Calculation parameters --------
REAL ::  ZK1,ZK2,ZK3            ! k1, k2 and K3 atmospheric refractivity constants
REAL  :: ZRDSRV                 ! XRD/XRV
!-------- Loop parameters ---------------
INTEGER :: IIB,IIE          ! Loop limits for coordinate X
INTEGER :: IJB,IJE          ! Loop limits for coordinate Y
INTEGER :: IKB,IKE          ! Loop limits for coordinate Z
INTEGER :: JK               ! Loop variables of control
INTEGER :: IIU,IJU,IKU      ! Loop variables of model 
REAL,  DIMENSION(:),ALLOCATABLE         :: ZXHATM,ZYHATM   ! mass-point positions 
REAL,  DIMENSION(:,:,:),ALLOCATABLE     :: ZZHATM   ! mass level altitude  
!-------- Physical parameters for the integration ----------------------------
REAL, DIMENSION(:,:,:),ALLOCATABLE ::  ZE              !  Partial pressure of water vapor
REAL, DIMENSION(:,:,:),ALLOCATABLE ::  ZTV             !  Virtual temperature
!-------- External model contributions ---------------------------------------
REAL, DIMENSION(:,:),ALLOCATABLE   ::  ZZHDX           !  Zenith hydrostatic delay external (of the model contribution)
!-------- Top model parameters -----------------------------------------------
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPTOP            !  pressure of top of model
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTVTOP           !  Tv of top of model
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZETOP            !  water vapor partial pressure of top of model
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTEMPTOP         !  absolute temperature of top of model
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZG0              !  gravity of ground (Vedel approx.)
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZR0              !  radius of ground (Vedel approx.)
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZGTOP            !  gravity of top of model using Vedel approx.
!  
!-------------------------------------------------------------------------------
!
!*       1.     INTIALIZE DIMENSIONS AND ALLOCATE ARRAYS
!               ----------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IIU = SIZE (PTEMP,1)
IJU = SIZE (PTEMP,2)
IKU = SIZE (PTEMP,3)
IKB = JPVEXT + 1
IKE = IKU - JPVEXT
!
ALLOCATE(ZXHATM(IIU))
ALLOCATE(ZYHATM(IJU))
ALLOCATE(ZZHATM(IIU,IJU,IKU))
ALLOCATE(ZE(IIU,IJU,IKU))
ALLOCATE(ZTV(IIU,IJU,IKU))
ALLOCATE(ZGTOP(IIU,IJU))
ALLOCATE(ZPTOP(IIU,IJU))
ALLOCATE(ZTVTOP(IIU,IJU))
ALLOCATE(ZTEMPTOP(IIU,IJU))
!
ALLOCATE(ZZHDX(IIU,IJU))
!
ALLOCATE(ZG0(IIU,IJU))
ALLOCATE(ZR0(IIU,IJU))
!
PZTD(:,:) = 0.    ! Zenithal Total Delay
PZHD(:,:) = 0.    ! Zenithal Hydrostatic Delay
PZWD(:,:) = 0.    ! Zenithal Wet Delay
!
ZZHDX(:,:)    = 0.
!-------------------------------------------------------------------------------
!
!*       3.    REFRACTIVITY  COEFFICIENTS AND OTHER CONSTANTS
!              ---------------------------------------------
!
! Refractivity coeficients
! Bevis et al. (1994)
ZK1 = 0.776       ! K/Pa
ZK2 = 0.704       ! K/Pa
ZK3 = 3739.       ! K2/Pa
ZRDSRV=XRD/XRV
!
!-------------------------------------------------------------------------------!
!*       4.    AUXILLARY VARIABLES
!              ------------------- 
!
! 
ZXHATM(1:IIU-1) = 0.5*(XXHAT(1:IIU-1)+XXHAT(2:IIU))
ZXHATM(IIU)     = 2.*XXHAT(IIU)-ZXHATM(IIU-1)
ZYHATM(1:IJU-1) = 0.5*(XYHAT(1:IJU-1)+XYHAT(2:IJU))
ZYHATM(IJU)     = 2.*XYHAT(IJU)-ZYHATM(IJU-1)  
ZZHATM(:,:,1:IKU-1)=0.5*(XZZ(:,:,1:IKU-1)+XZZ(:,:,2:IKU))
ZZHATM(:,:,IKU)    = 2.*XZZ(:,:,IKU) -ZZHATM(:,:,IKU-1)
! 
! Vapor pressure in Pa : ZE
ZE(:,:,:) = PPABST(:,:,:) *  PRT(:,:,:)  / ( ZRDSRV + PRT(:,:,:) )
! Virtual  temperature
ZTV(:,:,:) = PTEMP(:,:,:) * ( 1. + PRT(:,:,:) / ZRDSRV        )    &
     / ( 1. +  PRT(:,:,:)  )
ZTVTOP(:,:)= ( ZTV(:,:,IKE) - ZTV(:,:,IKE-1) ) *  ( XZZ(:,:,IKE+1) - ZZHATM(:,:,IKE) )  &
     /(ZZHATM(:,:,IKE) - ZZHATM(:,:,IKE-1))      +  ZTV(:,:,IKE)
! for extrapolation above model top : gtop, Ptop et Ttop
ZG0(:,:) = 9.780356 * ( 1. + 5.2885E-3 * (SIN(XLAT(:,:)*XPI/180.))**2 - &
                        5.9E-6 * ( SIN(2.*XLAT(:,:)*XPI/180.))**2 )
ZR0(:,:) = 6378.1E3  /  SQRT ( ( SIN(XLAT(:,:)*XPI/180.)*6378.1/6356.6 )**2 + &
          ( COS(XLAT(:,:)*XPI/180.))**2 )
ZGTOP(:,:) = ZG0(:,:) * ( ZR0(:,:) / (  ZR0(:,:) + XZZ(:,:,IKE+1)  ) )**2
ZPTOP(:,:) = PPABST(:,:,IKE)*EXP(XG*0.5*(XZZ(:,:,IKE)-XZZ(:,:,IKE+1))/(XRD*ZTVTOP(:,:)))
ZTEMPTOP(:,:) = (( PTEMP(:,:,IKE) - PTEMP(:,:,IKE-1) )  * 0.5 * ( XZZ(:,:,IKE+1) - XZZ(:,:,IKE) ) &
            /   ( XZZ(:,:,IKE) - XZZ(:,:,IKE-1) ))   +  PTEMP(:,:,IKE)
! 
!-------------------------------------------------------------------------------!
!*       5.    ZENITH DELAYS
!              -------------  
!
!        5.1  Integrated  Model Contribution 
!
!            5.1.1  ZHD 
!
! level contribution for integrated formulation
DO JK=IKB,IKE
   PZHD(:,:) = PZHD(:,:) + ( 1.E-6 * ZK1 * PPABST(:,:,JK) * &
                     ( XZZ(:,:,JK+1) - XZZ(:,:,JK) ) / ZTV(:,:,JK) )
END DO
!
!            5.1.2  ZWD 
! 
DO JK=IKB,IKE 
  PZWD(:,:) = PZWD(:,:) + (1.E-6 * ( (ZK2-ZRDSRV*ZK1) + ( ZK3/PTEMP(:,:,JK) ) ) * &
        ZE(:,:,JK)* ( XZZ(:,:,JK+1) - XZZ(:,:,JK) ) / PTEMP(:,:,JK))
END DO
!
!        5.2  External Model Contribution  (hydrostatic delay only, wet delay negligeable)
!
ZZHDX(:,:) = 1.E-6 * ZK1 * ZPTOP(:,:) * XRD * ( 1. + 2. * XRD * ZTEMPTOP(:,:) &
             / ( ( XRADIUS + XZZ(:,:,IKE+1) ) * ZGTOP(:,:) )  +  2. * ( XRD * ZTEMPTOP(:,:) &
            / ( (XRADIUS + XZZ(:,:,IKE+1)) * ZGTOP(:,:) ))**2 ) / ZGTOP(:,:)
    
PZHD(:,:) = PZHD(:,:) + ZZHDX(:,:)
!
!        5.3  Total Zenith delay 
!
PZTD(:,:) = PZHD(:,:) + PZWD(:,:) 
!  
!
DEALLOCATE(ZE)
DEALLOCATE(ZTV)
DEALLOCATE(ZGTOP)
DEALLOCATE(ZPTOP)
DEALLOCATE(ZTVTOP)
DEALLOCATE(ZTEMPTOP)
DEALLOCATE(ZZHDX)
! 
END SUBROUTINE GPS_ZENITH_GRID
