!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!##########################################
MODULE MODI_GPS_ZENITH
!##########################################
INTERFACE
      SUBROUTINE GPS_ZENITH(HFGRI,PRM,PTEMP,PPABSM,PTSRAD,PZTD,PZHD,PZWD)
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),            INTENT(IN)  :: HFGRI    ! ASCII file name
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PRM      ! Water vapor mixing ratio
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PTEMP    ! Air temperature
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PPABSM   ! Absolute pressure
REAL,  DIMENSION(:,:),       INTENT(IN)  :: PTSRAD   ! Surface temperature
!
REAL,  DIMENSION(:,:),   INTENT(OUT)     :: PZTD    ! Zenithal Total Delay
REAL,  DIMENSION(:,:),   INTENT(OUT)     :: PZHD    ! Zenithal Hydrostatic Delay
REAL,  DIMENSION(:,:),   INTENT(OUT)     :: PZWD    ! Zenithal Wet Delay
!
END SUBROUTINE GPS_ZENITH
END INTERFACE
END MODULE MODI_GPS_ZENITH
!
!     ###################################################################
      SUBROUTINE GPS_ZENITH(HFGRI,PRM,PTEMP,PPABSM,PTSRAD,PZTD,PZHD,PZWD)
!     ###################################################################
!
!!****  *GPS_ZENITH * - computes  GPS zenithal delays parameters
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
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_CST
USE MODD_DIAG_FLAG
USE MODD_GR_FIELD_n
USE MODD_GRID,             ONLY: XLONORI,XLATORI
USE MODD_GRID_n
USE MODE_GRIDPROJ
USE MODD_IO_ll,            ONLY: TFILEDATA
USE MODD_LUNIT
USE MODD_PARAMETERS
!
USE MODE_FM,               ONLY: IO_FILE_CLOSE_ll,IO_FILE_OPEN_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST
USE MODE_TOOLS_ll,         ONLY: LWEST_ll, LEAST_ll, LNORTH_ll, LSOUTH_ll
!
USE MODI_INTERPOL_STATION
!
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*),            INTENT(IN)  :: HFGRI      ! ASCII file name
!
REAL,  DIMENSION(:,:,:),   INTENT(IN)    :: PRM        ! Water vapour mixing ratio
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PTEMP      ! Air temperature
REAL,  DIMENSION(:,:,:),     INTENT(IN)  :: PPABSM     ! Absolute pressure
REAL,  DIMENSION(:,:),       INTENT(IN)  :: PTSRAD     ! Surface temperature
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
INTEGER :: JI,JJ,JK               ! Loop variables of control
INTEGER :: IIU,IJU,IKU      ! Loop variables of model 
INTEGER :: ILUOUT0, IRESP   ! file unit and return code for output
INTEGER :: JL               !
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
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRTTOP           !  mix ratio of top of model
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTEMPTOP         !  absolute temperature of top of model
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZG0              !  gravity of ground (Vedel approx.)
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZR0              !  radius of ground (Vedel approx.)
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZGTOP            !  gravity of top of model using Vedel approx.
      
!--------- Parameters to assess GPS sites delays ----------------------------- 

INTEGER :: ISTATIONS,IFGRI
REAL :: ZZM_STAT,ZTM_STAT,ZTV_STAT,ZPM_STAT,ZEM_STAT,ZRV_STAT
REAL:: ZPTOP_STAT,ZGTOP_STAT,ZTEMPTOP_STAT
INTEGER, DIMENSION(:), ALLOCATABLE :: II,IJ  
LOGICAL, DIMENSION(:), ALLOCATABLE :: GSTATION
INTEGER :: II1,IJ1,IX,IY,IBOT
REAL :: ZI,ZJ,ZST_PROC
REAL, DIMENSION(:), ALLOCATABLE ::ZXHAT_GPS, ZYHAT_GPS
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZ_STA
REAL, DIMENSION(:),ALLOCATABLE :: ZGPS_ZTD, ZGPS_ZHD, ZGPS_ZWD
REAL, DIMENSION(:),ALLOCATABLE :: ZE_STA, ZP_STA, ZT_STA, ZTV_STA
!--------- Misc --------------------------------------------------------------
TYPE(TFILEDATA),POINTER :: TZFILE
!-------------------------------------------------------------------------------
!
!*       1.     INTIALIZE DIMENSIONS AND ALLOCATE ARRAYS
!               ----------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
TZFILE => NULL()
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
ZE(:,:,:) = PPABSM(:,:,:) *  PRM(:,:,:)  / ( ZRDSRV + PRM(:,:,:) )
! Virtual  temperature
ZTV(:,:,:) = PTEMP(:,:,:) * ( 1. + PRM(:,:,:) / ZRDSRV        )    &
     / ( 1. +  PRM(:,:,:)  )
ZTVTOP(:,:)= ( ZTV(:,:,IKE) - ZTV(:,:,IKE-1) ) *  ( XZZ(:,:,IKE+1) - ZZHATM(:,:,IKE) )  &
     /(ZZHATM(:,:,IKE) - ZZHATM(:,:,IKE-1))      +  ZTV(:,:,IKE)
! for extrapolation above model top : gtop, Ptop et Ttop
ZG0(:,:) = 9.780356 * ( 1. + 5.2885E-3 * (SIN(XLAT(:,:)*XPI/180.))**2 - &
                        5.9E-6 * ( SIN(2.*XLAT(:,:)*XPI/180.))**2 )
ZR0(:,:) = 6378.1E3  /  SQRT ( ( SIN(XLAT(:,:)*XPI/180.)*6378.1/6356.6 )**2 + &
          ( COS(XLAT(:,:)*XPI/180.))**2 )
ZGTOP(:,:) = ZG0(:,:) * ( ZR0(:,:) / (  ZR0(:,:) + XZZ(:,:,IKE+1)  ) )**2
ZPTOP(:,:) = PPABSM(:,:,IKE)*EXP(XG*0.5*(XZZ(:,:,IKE)-XZZ(:,:,IKE+1))/(XRD*ZTVTOP(:,:)))
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
   PZHD(:,:) = PZHD(:,:) + ( 1.E-6 * ZK1 * PPABSM(:,:,JK) * &
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
!-------------------------------------------------------------------------------!
!*       6.     STATION VALUES  
!              --------------
!
!        6.1 Some initializations
!
! 
ISTATIONS=COUNT(XLAT_GPS(:) < XUNDEF)
ISTATIONS=MIN(ISTATIONS,COUNT(XLON_GPS(:)< XUNDEF))
ISTATIONS=MIN(ISTATIONS,COUNT(XZS_GPS(:)> -999.0))
PRINT *,'Number of GPS STATIONS ', ISTATIONS
!
IF (ISTATIONS >0 ) THEN 
!
  CALL IO_FILE_ADD2LIST(TZFILE,HFGRI,'GPS','WRITE')
  CALL IO_FILE_OPEN_ll(TZFILE)
  IFGRI = TZFILE%NLU
  PRINT *,'File ',TRIM(HFGRI),' opened with unit= ',IFGRI,' IRESP= ',IRESP
  WRITE(IFGRI,*,IOSTAT=IRESP) 'Number of STATIONS', ISTATIONS
! 
  ALLOCATE(ZGPS_ZTD(ISTATIONS))
  ALLOCATE(ZGPS_ZHD(ISTATIONS))
  ALLOCATE(ZGPS_ZWD(ISTATIONS))
  ALLOCATE(ZE_STA(IKU),ZP_STA(IKU))
  ALLOCATE(ZT_STA(IKU),ZTV_STA(IKU))
  ALLOCATE(ZZ_STA(IKU,ISTATIONS))
  ALLOCATE(II(ISTATIONS),IJ(ISTATIONS))
  ALLOCATE(ZXHAT_GPS(ISTATIONS),ZYHAT_GPS(ISTATIONS))
  ALLOCATE(GSTATION(ISTATIONS))
!
  ZGPS_ZTD(:) = 0. 
  ZGPS_ZWD(:) = 0.
  ZGPS_ZHD(:) = 0.
!
!
  DO JL=1,ISTATIONS ! LOOP over the GPS stations 
!
!        6.2 Compute positions inside the Meso-NH domain 
!
    GSTATION(JL)=.TRUE.
    CALL SM_XYHAT(XLATORI,XLONORI,   &
       XLAT_GPS(JL),XLON_GPS(JL),ZXHAT_GPS(JL),ZYHAT_GPS(JL))
!
    II(JL)=COUNT(ZXHATM(:)<=ZXHAT_GPS(JL))
    IX=COUNT(XXHAT(:)<=ZXHAT_GPS(JL))
    IF (IX<IIB .AND. LWEST_ll()) THEN
     ! station outside the MESO-NH domain 
      GSTATION(JL)=.FALSE.
    ENDIF
    IF (IX>IIE .AND. LEAST_ll()) THEN
     ! station outside the MESO-NH domain 
      GSTATION(JL)=.FALSE.
    ENDIF
    IJ(JL)=COUNT(ZYHATM(:)<=ZYHAT_GPS(JL))
    IY=COUNT(XYHAT(:)<=ZYHAT_GPS(JL))
    IF (IY<IJB .AND. LSOUTH_ll()) THEN
     ! stations outside MESO-NH domain 
      GSTATION(JL)=.FALSE.
    ENDIF
    IF (IY>IJE .AND. LNORTH_ll()) THEN
     ! stations outside the MESO-NH domain
      GSTATION(JL)=.FALSE.
    ENDIF
    IF (.NOT.GSTATION(JL)  ) THEN
      WRITE(IFGRI,*) 'NO DATA, GPS station ',CNAM_GPS(JL),' not in the domain Lat:',XLAT_GPS(JL),' Lon:',XLON_GPS(JL)
    END IF
!     
! Position of the station according to processors
!      
    ZST_PROC=0.
    IF (IX>=IIB .AND. IX <=IIE .AND. IY >=IJB .AND. IY <=IJE)  ZST_PROC=1.
!
    IF ( ZST_PROC >0. .AND. GSTATION(JL)) THEN 
      II1=II(JL)
      IJ1=IJ(JL)
!     interpolate Z at station position and check that the difference between model relief and station altitude is weaker than XDIFORO
      CALL INTERPOL_STATION(XZZ(:,:,:),ZXHATM,ZYHATM,II1,IJ1,ZXHAT_GPS(JL),ZYHAT_GPS(JL),ZZ_STA(:,JL))
      IF ( ABS( ZZ_STA(IKB,JL)-XZS_GPS(JL)) > XDIFFORO ) THEN
        WRITE(IFGRI,*) 'NO DATA, Difference between the model orography and the GPS station height too large for ',CNAM_GPS(JL)
        GSTATION(JL)=.FALSE.
        CYCLE
      END IF
!
!        6.3 Interpolate to the station positions 
!
! interpolate model variables to obs point
      CALL INTERPOL_STATION(PPABSM(:,:,:),ZXHATM,ZYHATM,II1,IJ1,ZXHAT_GPS(JL),ZYHAT_GPS(JL),ZP_STA(:))
      CALL INTERPOL_STATION(ZE(:,:,:),ZXHATM,ZYHATM,II1,IJ1,ZXHAT_GPS(JL),ZYHAT_GPS(JL),ZE_STA(:))
      CALL INTERPOL_STATION(PTEMP(:,:,:),ZXHATM,ZYHATM,II1,IJ1,ZXHAT_GPS(JL),ZYHAT_GPS(JL),ZT_STA(:))
      CALL INTERPOL_STATION(ZTV(:,:,:),ZXHATM,ZYHATM,II1,IJ1,ZXHAT_GPS(JL),ZYHAT_GPS(JL),ZTV_STA(:))
!  
!            6.3.1 For stations above model orography
!  
      IF (ZZ_STA(IKB,JL) < XZS_GPS(JL)) THEN ! station above the model orography                           
        DO JK=IKB+1,IKE
          IF ( ZZ_STA(JK,JL) <   XZS_GPS(JL)) THEN ! whole layer to remove
           IBOT=JK
          ELSE      ! partial layer to add 
            ZGPS_ZHD(JL)=ZGPS_ZHD(JL)+ ( 1.E-6 * ZK1 * ZP_STA(JK-1) * ( ZZ_STA(JK,JL) - XZS_GPS(JL) ) / ZTV_STA(JK-1)) 
            ZGPS_ZWD(JL)=ZGPS_ZWD(JL)+ ( 1.E-6 *  ( (ZK2-ZRDSRV*ZK1) + ( ZK3/ZT_STA(JK-1) ) ) * &
                 ZE_STA(JK-1)* ( ZZ_STA(JK,JL) - XZS_GPS(JL)  ) / ZT_STA(JK-1) ) 
            IBOT=JK
            EXIT
          END IF
        END DO
!        
!            6.3.2 For station below the model orography
!       
      ELSE  ! station below the model orography 
        IBOT=IKB
! Extrapolate variables below the model orography: 
        ZZM_STAT=0.5*(XZS_GPS(JL)+ZZ_STA(IKB,JL))
! assume constant temperature gradient
        ZTM_STAT=ZT_STA(IKB) + ( (ZZM_STAT-0.5*(ZZ_STA(IKB,JL)+ZZ_STA(IKB+1,JL)))*&
             ( PTEMP(II(JL),IJ(JL),IKB)- PTEMP(II(JL),IJ(JL),IKB+1))/&
             (ZZHATM(II(JL),IJ(JL),IKB)-ZZHATM(II(JL),IJ(JL),IKB+1)))
! assume constant virtual temperature gradient
        ZTV_STAT=ZTV_STA(IKB) + ( (ZZM_STAT-0.5*(ZZ_STA(IKB,JL)+ZZ_STA(IKB+1,JL)))*&
             ( ZTV(II(JL),IJ(JL),IKB)- ZTV(II(JL),IJ(JL),IKB+1))/&
             (ZZHATM(II(JL),IJ(JL),IKB)-ZZHATM(II(JL),IJ(JL),IKB+1)))
! assume hydrostatic law
        ZPM_STAT = ZP_STA(IKB) * EXP(XG *(ZZM_STAT-0.5*(ZZ_STA(IKB,JL)+ZZ_STA(IKB+1,JL)))&
             /(XRD* 0.5 *(ZTV_STAT+ZTV_STA(IKB))))
        
! assume constant rvn for Vapor pressure
        CALL INTERPOL_STATION( PRM(:,:,IKB),ZXHATM,ZYHATM,II(JL),IJ(JL),ZXHAT_GPS(JL), &
             ZYHAT_GPS(JL),ZRV_STAT)
        ZEM_STAT =  ZPM_STAT *  ZRV_STAT  / ( ZRDSRV + ZRV_STAT )
! add contribution below the model orography        
        ZGPS_ZHD(JL) = ZGPS_ZHD(JL) + ( 1.E-6 * ZK1 * ZPM_STAT * &
             ( ZZ_STA(IKB,JL) - XZS_GPS(JL) ) / ZTV_STAT )
        ZGPS_ZWD (JL)=ZGPS_ZWD(JL) + ( 1.E-6 * ( (ZK2-ZRDSRV*ZK1) + ( ZK3/ZTM_STAT ) )&
             * ZEM_STAT* (  ZZ_STA(IKB,JL) - XZS_GPS(JL) ) / ZTM_STAT  )
      ENDIF
!  
!            6.3.3 Add Full model layer contributions
!  
      DO JK=IBOT,IKE
        ZGPS_ZHD(JL)=ZGPS_ZHD(JL)+  ( 1.E-6 * ZK1 * ZP_STA(JK) * &
             ( ZZ_STA(JK+1,JL) - ZZ_STA(JK,JL) ) / ZTV_STA(JK) )
!             
        ZGPS_ZWD(JL)=ZGPS_ZWD(JL)+ (1.E-6 * ( (ZK2-ZRDSRV*ZK1) + ( ZK3/ZT_STA(JK) ) ) * &
             ZE_STA(JK)* ( ZZ_STA(JK+1,JL) - ZZ_STA(JK,JL) ) / ZT_STA(JK))
      END DO
!
!            6.3.4 Add external  contribution for ZHD
!
      CALL INTERPOL_STATION( ZPTOP(:,:),ZXHATM,ZYHATM,II(JL),IJ(JL),ZXHAT_GPS(JL), &
           ZYHAT_GPS(JL),ZPTOP_STAT)
      CALL INTERPOL_STATION( ZGTOP(:,:),ZXHATM,ZYHATM,II(JL),IJ(JL),ZXHAT_GPS(JL), &
           ZYHAT_GPS(JL),ZGTOP_STAT)
      CALL INTERPOL_STATION( ZTEMPTOP(:,:),ZXHATM,ZYHATM,II(JL),IJ(JL),ZXHAT_GPS(JL), &
           ZYHAT_GPS(JL),ZTEMPTOP_STAT)
      ZGPS_ZHD(JL)=ZGPS_ZHD(JL)+ ( 1.E-6 * ZK1 * ZPTOP_STAT * XRD * ( 1. + 2. * XRD * ZTEMPTOP_STAT &
           / ( ( XRADIUS + ZZ_STA(IKE+1,JL) ) * ZGTOP_STAT )  +  2. * ( XRD * ZTEMPTOP_STAT &
           / ( (XRADIUS + ZZ_STA(IKE+1,JL)) * ZGTOP_STAT ))**2 ) / ZGTOP_STAT)
!  
!            6.3.5 Final ZTD
!
      ZGPS_ZTD(JL) =  ZGPS_ZHD(JL)+ ZGPS_ZWD(JL)
      WRITE(UNIT=IFGRI,FMT=1000,IOSTAT=IRESP) CNAM_GPS(JL), ZGPS_ZHD(JL), &
                                              ZGPS_ZWD(JL), ZGPS_ZTD(JL)
    ENDIF ! station included in the domain and proc
!    
  END DO  ! end loop over the GPS Stations                             
!-------------------------------------------------------------------------------
!
!*       7.    DEALLOCATE ALL THE ARRAYS 
!
  1000 FORMAT('STATION ',A10,' ZHD: ',F8.5,' ZWD: ',F8.5,' ZTD: ',F8.5)
! 
  CALL IO_FILE_CLOSE_ll(TZFILE,IRESP)
  PRINT *,'File ',TRIM(HFGRI),' closed, IRESP= ',IRESP
!
  DEALLOCATE(ZXHATM)
  DEALLOCATE(ZYHATM)
  DEALLOCATE(ZGPS_ZTD)
  DEALLOCATE(ZGPS_ZHD)
  DEALLOCATE(ZGPS_ZWD)
  DEALLOCATE(ZE_STA)
  DEALLOCATE(ZP_STA)
  DEALLOCATE(ZT_STA)
  DEALLOCATE(ZTV_STA)
  DEALLOCATE(ZZ_STA)
END IF ! ISTATIONS test
!
DEALLOCATE(ZE)
DEALLOCATE(ZTV)
DEALLOCATE(ZGTOP)
DEALLOCATE(ZPTOP)
DEALLOCATE(ZTVTOP)
DEALLOCATE(ZTEMPTOP)
DEALLOCATE(ZZHDX)
! 
END SUBROUTINE GPS_ZENITH
