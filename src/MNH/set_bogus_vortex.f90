!MNH_LIC Copyright 2001-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ############################
      MODULE MODI_SET_BOGUS_VORTEX
!     ############################
INTERFACE
!
SUBROUTINE SET_BOGUS_VORTEX(PUT,PVT,PTHT)
!
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PUT  ! Pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PVT  ! Pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PTHT ! Dry potential temperature
!
END SUBROUTINE SET_BOGUS_VORTEX
!
END INTERFACE
!
END MODULE MODI_SET_BOGUS_VORTEX
!
!
!     #########################################
      SUBROUTINE SET_BOGUS_VORTEX(PUT,PVT,PTHT)
!     #########################################
!
!!**  *SET_BOGUS_VORTEX* - add a bogus vortex in a filtered large-scale
!!                         analyse and balance it using thermal wind balance
!!                         relation to prescrib the mass field.
!!
!!    PURPOSE
!!    -------
!      The purpose of this routine is to add a pseudo analytical bogus vortex
!      from Doppler radar observations
!
!!**  METHOD
!!    ------
!!     The bogus vortex is retrieved from the EVTD methodology
!!     (Roux and Marks 1996, Roux et al. 2003), and then the tangential wind
!!     is adjusted to a modified version of the Holland's analytical wind profil
!!     (Holland 1980). 
!!     
!!
!!    EXTERNAL
!!    --------
!!      FMATTR
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    
!!    none
!!
!!    REFERENCE
!!    ---------
!!   This routine is not available in Meso-NH documentation yet.
!!
!!    AUTHOR
!!    ------
!!	O. Nuissier         * L.A *
!!      R. Rogers           * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!      Modification   01/02/08 (D.Barbary) Add convergence angle in bogus
!!                                          and use modd_hurr_param for Holland's parameters
!!                     20/02/08 (D.Barbary) Change condition of ZRADBOGMAX
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODE_GRIDPROJ
USE MODE_MSG
!
USE MODD_HURR_CONF,  ONLY: XLATBOG,XLONBOG,XVTMAXSURF,XRADWINDSURF, &
                           XANGCONV0,XANGCONV1000,XANGCONV2000, &
			   XB_0, XMAX
USE MODD_PARAMETERS, ONLY: XUNDEF,JPVEXT
USE MODD_LUNIT,      ONLY: TLUOUT0
USE MODD_CST,        ONLY: XPI,XOMEGA 
USE MODD_CONF,       ONLY: NVERB
USE MODD_GRID,       ONLY: XLONORI,XLATORI
USE MODD_GRID_n,     ONLY: XXHAT,XYHAT,XZHAT,XLAT
USE MODD_REF,        ONLY: XTHVREFZ
!
USE MODI_SHUMAN
USE MODI_HOLLAND_VT   ! Module which compute Vt(r,z)
USE MODI_THERM_WIND_BAL
!
IMPLICIT NONE
!
!
!*       0.1   Dummy arguments
!
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PUT  ! Pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PVT  ! Pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PTHT ! Dry potential temperature
!
!*       0.2   Local variables
!
INTEGER                                                   :: IRET,ILUOUT0
INTEGER                                                   :: IREF_MAX
!
REAL,DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZTH_BOG
REAL,DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZU_BOG
REAL,DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZV_BOG
REAL,DIMENSION(:,:,:), ALLOCATABLE :: ZZHAT3D,ZZHATM ! altitude of mass points
REAL,DIMENSION(:)    , ALLOCATABLE :: ZZ_REF,ZTH_REF,ZZETA
REAL,DIMENSION(:)    , ALLOCATABLE :: ZVT_GRID,ZVR_GRID,ZTH_GRID
REAL,DIMENSION(:)    , ALLOCATABLE :: ZANGCONV ! Angle de convergence
					       ! en fonction de l altitude
!
REAL             :: ZRADBOGMAX ! Maximum for the bogus vortex
REAL             :: ZXK,ZYK
REAL             :: ZI,ZJ      ! Fractionnal indexes of the first guess
REAL             :: ZR,ZPHI
REAL             :: ZRRB,ZRRB_MAX
REAL             :: ZTH_REF_MEAN
REAL             :: ZXHAT,ZYHAT
REAL             :: ZCORIO ! Coriolis parameter at (i,j)
REAL             :: ZRADSDG,ZRADIM_IB
REAL             :: ZVT, ZVT_0 ! Tangentiel Wind at Z=0 and R
!
INTEGER          :: JI,JJ,JK     ! Loop indexes along the x,y,z directions
INTEGER          :: II,IJ        ! Indexes of the point of the 1st guess
INTEGER          :: IIBOG,IJBOG
INTEGER          :: IIU,IJU,IKU,IKB,IKE
INTEGER          :: IRESP   ! Return code of FM-routines
CHARACTER(LEN=100) :: YMSG
!
!-------------------------------------------------------------------------------
!
!*       1. INITIALIZATION
!           --------------
!
ILUOUT0 = TLUOUT0%NLU
!
WRITE(ILUOUT0,'(A)')' Begin of SET_BOGUS_VORTEX routine'
!
IIU = SIZE(PUT,1)
IJU = SIZE(PUT,2)
IKU = SIZE(PUT,3)
IKB = 1+JPVEXT
IKE = IKU-JPVEXT
!
ZRADSDG = XPI / 180.
!
IF (XLATBOG == XUNDEF .OR. XLONBOG == XUNDEF) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_BOGUS_VORTEX','You do not specify a first guess position of the bogus (XLATBOG, XLONBOG)')
END IF
IF (XVTMAXSURF == XUNDEF) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_BOGUS_VORTEX','You do not specify a maximum for the tangential wind (XVTMAXSURF,m/s)')
END IF
IF (XRADWINDSURF == XUNDEF) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN',&
                 'SET_BOGUS_VORTEX','You do not specify a radius for the maximum tangential wind (XRADWINDSURF,m)')
END IF
!
ZU_BOG(:,:,:) = 0.
ZV_BOG(:,:,:) = 0.
ZTH_BOG(:,:,:)= 0.
!
!-------------------------------------------------------------------------------
!
!*     2. LOCALIZING ON THE MESO-NH GRID THE POSITION OF THE NEW
!                      SPECIFIED BOGUS VORTEX
!         -------------------------------------------------------
!
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'Localizing the position of the bogus vortex'
!
CALL SM_XYHAT(XLATORI,XLONORI,XLATBOG,XLONBOG,ZXHAT,ZYHAT)
II=MAX(MIN(COUNT(XXHAT(:)<ZXHAT),IIU-1),1)
IJ=MAX(MIN(COUNT(XYHAT(:)<ZYHAT),IJU-1),1)
ZI=(ZXHAT-XXHAT(II))/(XXHAT(II+1)-XXHAT(II))+FLOAT(II)
ZJ=(ZYHAT-XYHAT(IJ))/(XYHAT(IJ+1)-XYHAT(IJ))+FLOAT(IJ)
IIBOG = INT(ZI)
IJBOG = INT(ZJ)
IF (NVERB>=5) WRITE(ILUOUT0,'(A,I3,A,I3)')' equivalent indexes in the Meso-NH grid: I= ',IIBOG,' J= ',IJBOG
!
IF ( (IIBOG<1) .OR. (IIBOG>IIU+1) .OR. &
     (IJBOG<1) .OR. (IJBOG>IJU+1)      ) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_BOGUS_VORTEX','The first guess position of the vortex is not in the Meso-NH domain')
END IF
!
!
!-------------------------------------------------------------------------------
!
!*	 3. DETERMINATION OF RADIAL EXTENSION OF THE VORTEX
!         -----------------------------------------------
!        In order to assume a coherent representation of the vortex
!        in the different domains, we select a criterion on Vt which
!        determREAL1-BretB2.nam ine the maximum radius of the vortex
!        (see also Vt expression in holland_vt_calc.f90)
!        Here, the radius ZRmax is such as <---> Vt(ZRmax,0)=0.5 m /s
!        and Vt(ZR+1)>=Vt(ZR) if ZR>=ZRVmax
!
ZVT=XVTMAXSURF
ZVT_0=ZVT
ZRADBOGMAX=XRADWINDSURF
! Use Coriolis parameter at bogus localisation at fisrt approximation
ZCORIO = 2. *XOMEGA *SIN(XLAT(IIBOG,IJBOG)*ZRADSDG)
DO WHILE ( (ZVT > 0.5 ) .AND. (ZVT_0 >= ZVT) )
  ZVT_0=ZVT
  ZRADBOGMAX=ZRADBOGMAX+XRADWINDSURF
  ZRADIM_IB = (ZRADBOGMAX / XRADWINDSURF) ** XB_0
  ZVT = SQRT ( XVTMAXSURF**2 * EXP (1. - 1. / ZRADIM_IB ) / ZRADIM_IB + &
        (ZCORIO * ZRADBOGMAX*1000. / 2.)**2 ) - ABS(ZCORIO * ZRADBOGMAX*1000. / 2.)
END DO

!
IF (NVERB>=10) THEN
  WRITE(ILUOUT0,'(A,F5.1,A,F6.2,A)') &
               'For a maximum surface tangential wind of ',XVTMAXSURF,&
               'm/s and a surface radius of maximum wind of ',XRADWINDSURF,'km'
  WRITE(ILUOUT0,'(A,F7.2,A)') &
               ' a vortex with a radial extension of ',ZRADBOGMAX,'km is taken into account'
END IF
!
ZRADBOGMAX=ZRADBOGMAX*1000.  ! conversion from km to m
!
!-------------------------------------------------------------------------------
!
!*     4. COMPUTES FOR EACH MESONH GRID POINT AN ANALYTICAL EXPRESSION
!         OF TANGENTIAL WIND Vt AND DEDUCES THE TEMPERATURE FIELD
!          (OR THE MASS FIELD) WHICH BALANCES THE VORTEX
!
ALLOCATE(ZZHAT3D(1,1,IKU),ZZHATM(1,1,IKU)) ! to compute altitude of mass points
ZZHAT3D(1,1,:) = XZHAT(:)
ZZHATM = MZF(1,IKU,1,ZZHAT3D)
DEALLOCATE(ZZHAT3D)
!
! Definition de l angle de convergence
! Forme lineaire de 0-1km puis 1-2km puis angle constant (nul)
! d'après litterature
ALLOCATE(ZANGCONV(IKU))
ZANGCONV(:)=0.
DO JK = 1, IKU
  IF (ZZHATM(1,1,JK).LT.0.) ZANGCONV(JK)=XANGCONV0
  IF ((ZZHATM(1,1,JK).GE.0.).AND.(ZZHATM(1,1,JK).LT.1000.))         &
    ZANGCONV(JK)=XANGCONV0 + &
        (XANGCONV1000-XANGCONV0)/(1000.-0.)*(ZZHATM(1,1,JK)-0)
  IF ((ZZHATM(1,1,JK).GE.1000.).AND.(ZZHATM(1,1,JK).LT.2000.))       &
    ZANGCONV(JK)=XANGCONV1000 +              &
        (XANGCONV2000-XANGCONV1000)/(2000.-1000.)*(ZZHATM(1,1,JK)-1000)
  IF (ZZHATM(1,1,JK).GE.2000.)    &
    ZANGCONV(JK)=XANGCONV2000
END DO
! Conversion en radian
ZANGCONV(:) = ZANGCONV(:) * ZRADSDG

!
! standard atmosphere for the thermal wind balance
!is taken as the 3D MesoNH reference atmosphere (XTHVREF)
IREF_MAX=-1
DO JK = IKB+1, IKE
  IF((ZZHATM(1,1,JK+1).GT.XMAX).AND.(XMAX.GE.ZZHATM(1,1,JK))) IREF_MAX=JK-1
END DO
IF (IREF_MAX==-1) THEN
  WRITE(YMSG,'(A,F4.1,A)')'The maximum height for the bogus vortex (',XMAX,' km) is not in the vertical grid'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_BOGUS_VORTEX',YMSG)
END IF
ZTH_REF_MEAN = SUM(XTHVREFZ(IKB:IREF_MAX)) / (IREF_MAX-IKB+1)
IF (NVERB>=10) THEN
  WRITE(ILUOUT0,'(A,I2,A,F7.1,A)')'Reference atmosphere is taken until level ',&
       IREF_MAX,' (above ',XMAX,'m)'
  WRITE(ILUOUT0,'(A,F6.2,A)')' mean reference theta(K)= ',ZTH_REF_MEAN
END IF
!
! the bogus vortex does existe below XMAX
ALLOCATE(ZTH_REF(IREF_MAX),ZZ_REF(IREF_MAX),ZZETA(IREF_MAX))
ZTH_REF(:) = ZTH_REF_MEAN
ZZ_REF (:) = 0.
ZZETA  (:) = 0.
ZTH_REF(IKB:IREF_MAX) = XTHVREFZ(IKB:IREF_MAX)
ZZ_REF (IKB:IREF_MAX) = ZZHATM(1,1,IKB:IREF_MAX)
ZZETA  (IKB:IREF_MAX) = ZZ_REF(IKB:IREF_MAX) &
                        *ZTH_REF(IKB:IREF_MAX)/ZTH_REF_MEAN
DEALLOCATE(ZZHATM)
!
ALLOCATE(ZVT_GRID(IREF_MAX),ZVR_GRID(IREF_MAX),ZTH_GRID(IREF_MAX))
!
DO JJ = 1, IJU
  DO JI = 1, IIU
    ZXK = XXHAT(JI) - XXHAT(IIBOG)
    ZYK = XYHAT(JJ) - XYHAT(IJBOG)
    ZR  = SQRT(ZXK*ZXK + ZYK*ZYK)
    ZCORIO = 2. *XOMEGA *SIN(XLAT(JI,JJ)*ZRADSDG)
    !
    ! the bogus vortex does existe inside ZRADBOGMAX
    IF ( ZR < ZRADBOGMAX ) THEN
      IF ( ZR>0 ) THEN
        ZPHI = ATAN2(ZYK,ZXK)
        IF (ZPHI < 0.) ZPHI = ZPHI + 2*XPI
      !
      ! a) Vt obeys to the Holland's formulation 
        ! integration en ZETA 
        CALL HOLLAND_VT(ZZETA(IREF_MAX),ZZETA,ZR,ZCORIO,ZVT_GRID)
	IF (ZR < XRADWINDSURF*1000.) THEN
      ! Vr exist only for R>XRADWINDSURF
	  ZVR_GRID(1:IREF_MAX) = 0.
        ELSE
          ZVR_GRID(1:IREF_MAX) = ZVT_GRID(1:IREF_MAX) *    &
                                 TAN( ZANGCONV(1:IREF_MAX) )
        END IF
      ! Reverse the signe of tangentiel wind when it is in southern hemisphere
        IF (XLATBOG<0.) THEN
	  ZVT_GRID = -ZVT_GRID
        END IF

      !
      ! Compute u and v from Vt
        ZU_BOG(JI,JJ,1:IREF_MAX) = - ZVT_GRID(1:IREF_MAX) * SIN(ZPHI) - &
				   ZVR_GRID(1:IREF_MAX) * COS(ZPHI)
        ZV_BOG(JI,JJ,1:IREF_MAX) =   ZVT_GRID(1:IREF_MAX) * COS(ZPHI) - &
				   ZVR_GRID(1:IREF_MAX) * SIN(ZPHI)

      ENDIF
      !
      ! b) while mass field is obtained from thermal wind relation
      ! integration en ZETA 
      CALL THERM_WIND_BAL(ZZETA(IREF_MAX),ZCORIO,ZR,ZRADBOGMAX, &
                          ZZETA,ZTH_REF, &
                          ZTH_GRID)
      ZTH_BOG(JI,JJ,IKB:IREF_MAX) = ZTH_GRID(IKB:IREF_MAX)
      !
    END IF        
  END DO
END DO
!
!-------------------------------------------------------------------------------
!    
!*     5. SUPERPOSE THE SPECIFIED VORTEX ONTO THE ENVIRONMENTAL FIELDS
!         ----------------------------------------------------------
!
PUT(:,:,:) = PUT(:,:,:) + ZU_BOG(:,:,:)
PVT(:,:,:) = PVT(:,:,:) + ZV_BOG(:,:,:)
PTHT(:,:,:)= PTHT(:,:,:)+ ZTH_BOG(:,:,:)
!
WRITE(ILUOUT0,'(A)')' Routine SET_BOGUS_VORTEX completed'
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_BOGUS_VORTEX 
