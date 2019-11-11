!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_CSTN
!     ####################
!
INTERFACE
!
SUBROUTINE SET_CSTN(TPFILE,TPEXPREFILE,HFUNU,HFUNV,KILOC,KJLOC,OBOUSS,OPV_PERT,ORMV_BL,PJ,OSHIFT,PCORIOZ) 
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
TYPE(TFILEDATA),        INTENT(IN)  :: TPFILE ! outpput data file
TYPE(TFILEDATA),        INTENT(IN)  :: TPEXPREFILE ! input data file
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNU  ! type of variation of U
                                              ! in y direction
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNV  ! type of variation of V
                                              ! in x direction
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
LOGICAL,                INTENT(IN)  :: OPV_PERT! logical switch for PV inversion
LOGICAL,                INTENT(IN)  :: ORMV_BL! logical switch for remouve boundary layer
REAL, DIMENSION(:,:,:), INTENT(IN) :: PJ ! jacobien 
LOGICAL,                INTENT(IN)  :: OSHIFT ! logical switch for vertical shift
!
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCORIOZ ! Coriolis parameter
                                               ! (exceptionnaly 3D array)
!
END SUBROUTINE SET_CSTN
!
END INTERFACE
!
END MODULE MODI_SET_CSTN
!
!
!
!     ##################################################################################
      SUBROUTINE SET_CSTN(TPFILE,TPEXPREFILE,HFUNU,HFUNV,KILOC,KJLOC,OBOUSS,OPV_PERT,ORMV_BL,PJ,OSHIFT,PCORIOZ) 
!     ##################################################################################
!
!!****  *SET_CSTN * - routine to initialize mass and wind fields from a Nv=cste 
!!                    profile
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is  to initialize the mass (theta,r,
!     thetavrefz,rhorefz) and the wind fields on model grid from a 
!     vertical profile ( ILEVEL - 1 layers  of Nv=cste) located at point 
!     (KILOC,KJLOC) : 
!     
!                      
!                     |------ Z(ILEVEL),U(ILEVEL),V(ILEVEL),Hu(ILEVEL)
!                     |
!                     | Nv(ILEVEL-1)
!                     |
!                     |------ Z(ILEVEL-1),U(ILEVEL-1),V(ILEVEL-1),Hu(ILEVEL-1)
!                     |  .
!                     |  .
!                     |  .
!                     |  .
!                     |------ Z(2),U(2),V(2),Hu(2)
!                     |
!                     |Nv(1)
!                     |
!                     |
!  Pground Thvground  |------ Z(1), U(1),V(1),Hu(1) ---------------ground level
!
!               (KILOC,KJLOC)
!
!
!
!        The free-format part of EXPRE file contains the vertical profile.The data 
!     are stored in the following order :
!
!         - number of  data levels ( variable ILEVEL)
!         - Thetav at ground
!         - Pressure at ground
!         - height of the ILEVEL levels
!         - U-wind component of the ILEVEL levels
!         - V-wind component of the ILEVEL levels
!         - Relative humidity of the ILEVEL levels
!         - Moist Brunt Vaisala frequency of the ILEVEL-1 layers
!
!
!!**  METHOD
!!    ------
!!      Firstly, the vertical profile is read, then the virtual potential 
!!    temperature Thetav is retrieved from the integrated form of the Brunt 
!!    Vaisala frequency definition :
!!    
!!      thetav(k+1) = thetav(k) exp( Nv(k) * Nv(k) (z(k+1)-z(k))/ g)
!!
!!      Then, Thetav and the humidity are interpolated on a vertical grid which
!!    is which is a mixed grid calaculated with VERT_COORD from the vertical levels of MNH
!!    grid and with a constant ororgraphy equal to the altitude of the vertical
!!    profile (ZZGROUND) (It permits to keep low levels information with a
!!    shifting function (as in PREP_REAL_CASE))
!!    For the vapor mixing ratio computation, the pressure is first determined
!!    from the hydrostatic relation (done by PRESS_HEIGHT) and the virtual
!!    temperature is computed from the virtual potential  temperature :
!!                            P    Rd/Cpd  
!!          Tv    = Thetav ( ---- )        
!!                           P00                             
!!
!!      The vapor mixing ratio is determined by an iterative procedure from
!!    pressure,humidity and virtual temperature (done by SM_PMR_HU)
!!
!!      Then, the horizontal structures of the 3D mass fields are deduced
!!    in SET_MASS 
!!
!!    EXTERNAL
!!    --------
!!      Module MODE_THERMO : contains thermodynamic routines
!!         SM_PMR_HU : to compute vapor mixing ratio from pressure, virtual
!!                    temperature and relative humidity
!!      PRESS_HEIGHT : to compute pressure from height and potentail virtual
!!                     temperature
!!
!!      SET_MASS : to compute mass and wind fields on the 3D-model grid
!!
!!      Module MODI_PRESS_HEIGHT : interface for function PRESS_HEIGHT
!!      Module MODI_SET_MASS     : interface for subroutine SET_MASS
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST     : contains physical constants
!!        XG    : Gravity constant 
!!        XP00  : reference pressure
!!        XRD   : Gas constant  for dry air
!!        XCPD : Specific heat for dry air at constant pressure
!!
!!      Module MODD_LUNIT1  : contains logical unit names 
!!        CLUOUT : name of output-listing
!!
!!      Module MODD_CONF    : contains configuration variables for all models. 
!!        NVERB : verbosity level for output-listing
!!
!!      Module MODD_GRID1  : contains grid variables
!!        XZHAT : height of w-levels of vertical model grid without orography
!!  
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (routine SET_CSTN)
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/08/94
!!      change the call to SET_MASS                    5/11/94  (J.Stein)     
!!      change the computation of Rv in function of Hu 25/2/95  (J.Stein)  
!!      J.Stein     30/01/96  use the RS ground pressure to initialize the
!!                            hydrostatic pressure computation
!!      G. Tanguy   26/10/10  change the interpolation of the RS : we use now a
!!                            mixed grid (PREP_REAL_CASE method)
!!      V.Masson    12/08/13  Parallelization of the initilization profile
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_GRID_n
USE MODD_IO_ll,      ONLY : TFILEDATA
USE MODD_LUNIT_n,    ONLY: CLUOUT, TLUOUT
USE MODD_PARAMETERS, ONLY: JPHEXT
!
USE MODE_FM
USE MODE_THERMO
USE MODE_ll
USE MODE_MPPDB
!
USE MODI_PRESS_HEIGHT
USE MODI_SET_MASS
USE MODI_SHUMAN
USE MODI_VERT_COORD
!
IMPLICIT NONE
!  
!*       0.1   Declarations of arguments :
!
TYPE(TFILEDATA),        INTENT(IN)  :: TPFILE ! outpput data file
TYPE(TFILEDATA),        INTENT(IN)  :: TPEXPREFILE ! input data file
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNU  ! type of variation of U
                                              ! in y direction
CHARACTER(LEN=*),       INTENT(IN)  :: HFUNV  ! type of variation of V
                                              ! in x direction
INTEGER,                INTENT(IN)  :: KILOC  ! I Localisation of vertical profile
INTEGER,                INTENT(IN)  :: KJLOC  ! J Localisation of vertical profile
LOGICAL,                INTENT(IN)  :: OBOUSS ! logical switch for Boussinesq version
LOGICAL,                INTENT(IN)  :: OPV_PERT! logical switch for PV inversion
LOGICAL,                INTENT(IN)  :: ORMV_BL! logical switch for remouve boundary layer
REAL, DIMENSION(:,:,:), INTENT(IN) :: PJ ! jacobien 
LOGICAL,                INTENT(IN)  :: OSHIFT ! logical switch for vertical shift
!
REAL, DIMENSION(:,:,:), INTENT(OUT), OPTIONAL :: PCORIOZ ! Coriolis parameter
                                               ! (exceptionnaly 3D array)
!
!*       0.2   Declarations of local variables :
!
!
! fields and data on the sounding levels
!
INTEGER                         :: ILUPRE,IRESP ! logical unit number of the 
                                                ! EXPRE and FM return code
INTEGER                         :: ILUOUT    ! Logical unit number for
                                             ! output-listing   
INTEGER                         :: ILEVEL   ! number of levels 
INTEGER                         :: ILAYER   ! number of layers
REAL                            :: ZPGROUND ! pressure at the ground level
REAL, DIMENSION(:), ALLOCATABLE :: ZHEIGHT  ! Height of levels
REAL, DIMENSION(:), ALLOCATABLE :: ZU,ZV    ! wind components
REAL, DIMENSION(:), ALLOCATABLE :: ZHU      ! relative humidity
REAL, DIMENSION(:), ALLOCATABLE :: ZNV      ! Moist Brunt Vaisala frequency 
                                            ! for  layers
REAL, DIMENSION(:), ALLOCATABLE :: ZTHV     ! Thetav
INTEGER                         :: JK,JKLEV ! Loop indexes
!
!  variables on the grid without orography
!
REAL, DIMENSION(SIZE(XZHAT))    :: ZZHATM  ! Height of mass model grid levels 
                                           ! without orography 
REAL, DIMENSION(SIZE(XZHAT))    :: ZTHVM   ! Virtual potential Temperature
                                           ! at mass model grid levels
REAL, DIMENSION(SIZE(XZHAT))    :: ZTVM    ! Virtual Temperature at mass model
                                           ! grid levels
REAL, DIMENSION(SIZE(XZHAT))    :: ZUW,ZVW ! Wind at w model grid levels
REAL, DIMENSION(SIZE(XZHAT))    :: ZHUM    ! humidity at mass model
                                           ! grid levels 
REAL, DIMENSION(SIZE(XZHAT))    :: ZMRM    ! mixing ratio at mass model
                                           ! grid levels 
REAL, DIMENSION(SIZE(XZHAT))    :: ZPM     ! pressure at mass model
                                           ! grid levels 
REAL                            :: ZRDSCPD ! Rd/Cpd
REAL                            :: ZDZSDH,ZDZ1SDH,ZDZ2SDH ! interpolation
                                                          ! working arrays
REAL                            :: ZEXNGRDM, ZPGRDM       ! Exner function and 
                                          ! pressure at z = -delta z /2
INTEGER                         :: IKU    ! Upper bound in z direction of model
                                          ! arrays   
!
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT))   ::ZZS_LS
REAL,DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) ::ZZFLUX_MX,ZZMASS_MX ! mixed grid
REAL, DIMENSION(SIZE(XZHAT))    :: ZZFLUX_PROFILE ! altitude of flux points on the initialization columns
REAL, DIMENSION(SIZE(XZHAT))    :: ZZMASS_PROFILE ! altitude of mass points on the initialization columns
!
INTEGER         :: IIB, IIE, IJB, IJE
INTEGER         :: IXOR_ll, IYOR_ll
INTEGER         :: IINFO_ll
LOGICAL         :: GPROFILE_IN_PROC   ! T : initialization profile is in current processor
!
!-------------------------------------------------------------------------------
!
!*	 1.     PROLOGUE : RETRIEVE LOGICAL UNIT NUMBERS AND INITIALIZE SOME
!                          CONSTANTS
!	        ------------------------------------------------------------
!
ILUPRE = TPEXPREFILE%NLU
ILUOUT = TLUOUT%NLU
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
CALL GET_OR_ll('B',IXOR_ll,IYOR_ll)
!
ZRDSCPD = XRD / XCPD
!-------------------------------------------------------------------------------
!
!*	 2.     READ DATA
!	        ---------
!
!*       2.1    Read levels number of  and allocate memory
!
READ(ILUPRE,*) ILEVEL
ILAYER = ILEVEL -1
ALLOCATE(ZNV(ILAYER))
ALLOCATE(ZHEIGHT(ILEVEL),ZU(ILEVEL),ZV(ILEVEL),ZHU(ILEVEL))
ALLOCATE(ZTHV(ILEVEL))
!
!*       2.2    Read data
!
READ(ILUPRE,*) ZTHV(1)       ! The first level is at the ground 
READ(ILUPRE,*) ZPGROUND
READ(ILUPRE,*) ZHEIGHT
READ(ILUPRE,*) ZU
READ(ILUPRE,*) ZV
READ(ILUPRE,*) ZHU
READ(ILUPRE,*) ZNV
!
!-------------------------------------------------------------------------------
!
!*	 3.     COMPUTE THETAV 
!	        --------------
!
DO JK = 2,ILEVEL
  ZTHV(JK) = ZTHV(JK-1) * EXP((ZNV(JK-1)**2) * (ZHEIGHT(JK)-ZHEIGHT(JK-1))/XG)
END DO
!
!
!-------------------------------------------------------------------------------
!
!*	 4.     INTERPOLATE ON THE  VERTICAL MIXED MODEL GRID 
!	        ---------------------------------------------------------
!
IKU=SIZE(XZHAT)
!
!*       4.1    Compute mixed grid
!
IF (PRESENT(PCORIOZ)) THEN 
! LGEOSBAL=T (no shift allowed, MNH grid without ororgraphy)        
  ZZS_LS(:,:)=0
ELSE
ZZS_LS(:,:)=ZHEIGHT(1)
ENDIF
CALL VERT_COORD(LSLEVE,ZZS_LS,ZZS_LS,XLEN1,XLEN2,XZHAT,ZZFLUX_MX)
ZZMASS_MX(:,:,:)=MZF(1,IKU,1,ZZFLUX_MX)
ZZMASS_MX(:,:,IKU)=1.5*ZZFLUX_MX(:,:,IKU)-0.5*ZZFLUX_MX(:,:,IKU-1)
!
CALL MPPDB_CHECK3D(ZZMASS_MX,"SET_CSTN::ZZMASS_MX",PRECISION)
!
!* vertical grid at initialization profile location
GPROFILE_IN_PROC=(KILOC+JPHEXT-IXOR_ll+1>=IIB .AND. KILOC+JPHEXT-IXOR_ll+1<=IIE ) &
         & .AND. (KJLOC+JPHEXT-IYOR_ll+1>=IJB .AND. KJLOC+JPHEXT-IYOR_ll+1<=IJE)
!
IF (GPROFILE_IN_PROC) THEN
  ZZMASS_PROFILE(:) = ZZMASS_MX(KILOC+JPHEXT-IXOR_ll+1,KJLOC+JPHEXT-IYOR_ll+1,:)
  ZZFLUX_PROFILE(:) = ZZFLUX_MX(KILOC+JPHEXT-IXOR_ll+1,KJLOC+JPHEXT-IYOR_ll+1,:)
ELSE
  ZZMASS_PROFILE(:) = 0.
  ZZFLUX_PROFILE(:) = 0.
END IF
DO JK = 1,IKU
  CALL REDUCESUM_ll(ZZMASS_PROFILE(JK), IINFO_ll)
  CALL REDUCESUM_ll(ZZFLUX_PROFILE(JK), IINFO_ll)
END DO
!
!*       4.2    Interpolate and extrapolate U and V on w-mixed grid levels :
!      
DO JK = 1,IKU
  IF (ZZFLUX_PROFILE(JK) <= ZHEIGHT(1)) THEN        ! extrapolation below the first level
    ZDZSDH  = (ZZFLUX_PROFILE(JK) - ZHEIGHT(1)) / (ZHEIGHT(2) - ZHEIGHT(1))
    ZUW(JK) = ZU(1) + (ZU(2) - ZU(1)) * ZDZSDH
    ZVW(JK) = ZV(1) + (ZV(2) - ZV(1)) * ZDZSDH 
  ELSE IF (ZZFLUX_PROFILE(JK) > ZHEIGHT(ILEVEL) ) THEN    ! extrapolation above the last
    ZDZSDH  = (ZZFLUX_PROFILE(JK) - ZHEIGHT(ILEVEL))  &   ! level               
            / (ZHEIGHT(ILEVEL) - ZHEIGHT(ILEVEL-1))
    ZUW(JK) = ZU(ILEVEL) + (ZU(ILEVEL) -ZU(ILEVEL -1)) * ZDZSDH
    ZVW(JK) = ZV(ILEVEL) + (ZV(ILEVEL) -ZV(ILEVEL -1)) * ZDZSDH
  ELSE        ! interpolation between the first and last levels
    DO JKLEV = 1,ILEVEL-1
      IF ( (ZZFLUX_PROFILE(JK) > ZHEIGHT(JKLEV)).AND.(ZZFLUX_PROFILE(JK) <= ZHEIGHT(JKLEV+1))) THEN 
        ZDZ1SDH = (ZZFLUX_PROFILE(JK) - ZHEIGHT(JKLEV))                 &
              / (ZHEIGHT(JKLEV+1)-ZHEIGHT(JKLEV)) 
        ZDZ2SDH = 1.- ZDZ1SDH 
        ZUW(JK) = (ZU(JKLEV) * ZDZ2SDH) + (ZU(JKLEV+1) *ZDZ1SDH)
        ZVW(JK) = (ZV(JKLEV) * ZDZ2SDH) + (ZV(JKLEV+1) *ZDZ1SDH)
      END IF 
    END DO
  END IF
END DO
!
!
!*       4.3    Interpolate and extrapolate Thetav and Hu on mass mixed grid levels
!                  
DO JK = 1,IKU
  IF (ZZMASS_PROFILE(JK) <= ZHEIGHT(1)) THEN             ! extrapolation below the first
    ZDZSDH  = (ZZMASS_PROFILE(JK) - ZHEIGHT(1)) / (ZHEIGHT(2) - ZHEIGHT(1)) !  level
    ZTHVM(JK) = ZTHV(1) * EXP((ZNV(1)**2) * (ZZMASS_PROFILE(JK) - ZHEIGHT(1))/XG)
    ZHUM(JK)  = ZHU(1) + (ZHU(2) - ZHU(1)) * ZDZSDH 
  ELSE IF (ZZMASS_PROFILE(JK) > ZHEIGHT(ILEVEL) ) THEN   ! extrapolation above the last
    ZDZSDH  = (ZZMASS_PROFILE(JK) - ZHEIGHT(ILEVEL))          &        ! level
            / (ZHEIGHT(ILEVEL) - ZHEIGHT(ILEVEL-1))
    ZTHVM(JK) = ZTHV(ILEVEL) + (ZTHV(ILEVEL) -ZTHV(ILEVEL -1)) * ZDZSDH
    ZHUM(JK) = ZHU(ILEVEL) + (ZHU(ILEVEL) -ZHU(ILEVEL -1)) * ZDZSDH
  ELSE             ! interpolation between the first and last levels
    DO JKLEV = 1,ILEVEL-1
      IF ( (ZZMASS_PROFILE(JK) > ZHEIGHT(JKLEV)).AND.                  &
           (ZZMASS_PROFILE(JK) <= ZHEIGHT(JKLEV+1)) ) THEN 
!
!     a linear interpolation is used for the humidity field and the 
!     logarithmic variation law (Nv = cst) is used for THV 
        ZDZ1SDH = (ZZMASS_PROFILE(JK) - ZHEIGHT(JKLEV))                 &
              / (ZHEIGHT(JKLEV+1)-ZHEIGHT(JKLEV)) 
        ZDZ2SDH = 1. -ZDZ1SDH
        ZHUM(JK)  = (ZHU(JKLEV) * ZDZ2SDH) + (ZHU(JKLEV+1) *ZDZ1SDH)
        ZTHVM(JK) = ZTHV(JKLEV) * EXP         &
        ( (ZNV(JKLEV)**2) * ( ZZMASS_PROFILE(JK) - ZHEIGHT(JKLEV) ) /XG )
      END IF 
    END DO
  END IF
END DO
!
!
!*       4.3    Compute Mixing ratio
!
! determines the pressure under the ground
ZEXNGRDM= ( ZPGROUND / XP00) ** ZRDSCPD  &
          - XG/XCPD  / (0.5*(ZTHV(1)+ZTHVM(1))) * (ZZMASS_PROFILE(1) - ZHEIGHT(1))
ZPGRDM = XP00 * ZEXNGRDM ** (1./ZRDSCPD)
ZPM(:)  = PRESS_HEIGHT(ZZMASS_PROFILE(:),ZTHVM,ZPGRDM,ZTHVM(1),ZZMASS_PROFILE(1)) ! compute P
ZTVM(:) = ZTHVM(:) * (ZPM(:) / XP00) ** ZRDSCPD                ! compute Tv
ZMRM(:) = SM_PMR_HU(ZPM(:),ZTVM(:),ZHUM(:),      &
                              SPREAD(ZMRM(:),2,1))             ! compute vapor
                                                               ! mixing ratio  
!-------------------------------------------------------------------------------
!
!*	 5.     COMPUTE FIELDS ON THE MODEL GRID (WITH OROGRAPHY)
!	        -------------------------------------------------
!
IF (PRESENT(PCORIOZ)) THEN
    CALL SET_MASS(TPFILE,GPROFILE_IN_PROC, ZZFLUX_PROFILE,                       &
                  KILOC+JPHEXT,KJLOC+JPHEXT,ZZS_LS,ZZMASS_MX,ZZFLUX_MX,ZPGROUND, &
                  ZTHVM,ZMRM,ZUW,ZVW,OSHIFT,OBOUSS,PJ,HFUNU,HFUNV,PCORIOZ=PCORIOZ)
ELSE
    CALL SET_MASS(TPFILE,GPROFILE_IN_PROC, ZZFLUX_PROFILE,                       &
                  KILOC+JPHEXT,KJLOC+JPHEXT,ZZS_LS,ZZMASS_MX,ZZFLUX_MX,ZPGROUND, &
                  ZTHVM,ZMRM,ZUW,ZVW,OSHIFT,OBOUSS,PJ,HFUNU,HFUNV)
ENDIF
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_CSTN
