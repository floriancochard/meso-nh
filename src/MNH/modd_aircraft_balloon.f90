!MNH_LIC Copyright 2000-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ############################
      MODULE MODD_AIRCRAFT_BALLOON
!     ############################
!
!!****  *MODD_AIRCRAFT_BALLOON* - declaration of balloons
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!      the different balloons types.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!
!!    AUTHOR
!!    ------
!! P. Jabouille   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/05/00
!!              Apr,19, 2001 (G.Jaubert) add CVBALL type
!!              March, 2013 : O.Caumont, C.Lac : add vertical profiles
!!              Oct,2016 : G.DELAUTIER LIMA
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_TYPE_DATE
!
TYPE FLYER
!
!
!* general information
!
CHARACTER(LEN=3)              :: MODEL  ! type of model used for each balloon/aircraft
                                        ! 'FIX' : NMODEL used during the run
                                        ! 'MOB' : change od model depends of the 
                                        !         balloon/aircraft location
INTEGER                       :: NMODEL ! model number for each balloon/aircraft
CHARACTER(LEN=6)              :: TYPE   ! flyer type:
                                        ! 'RADIOS' : radiosounding balloon
                                        ! 'ISODEN' : iso-density balloon
                                        ! 'AIRCRA' : aircraft
                                        ! 'CVBALL' : Constant Volume balloon
CHARACTER(LEN=10)             :: TITLE  ! title or name for the balloon/aircraft
TYPE(DATE_TIME)               :: LAUNCH ! launch/takeoff date and time
LOGICAL                       :: CRASH  ! occurence of crash
LOGICAL                       :: FLY    ! occurence of flying
!
!* storage monitoring
!
REAL                          :: T_CUR  ! current time since last storage
INTEGER                       :: N_CUR  ! current step of storage
REAL                          :: STEP   ! storage time step
!
!* balloon dynamical characteristics
!
REAL                          :: LAT    ! latitude of launch
REAL                          :: LON    ! lontitude of launch
REAL                          :: XLAUNCH! X coordinate of launch
REAL                          :: YLAUNCH! Y coordinate of launch
REAL                          :: ALT    ! altitude of launch (if 'RADIOS' or 'ISODEN' or 'CVBALL')
REAL                          :: WASCENT! ascent vertical speed (if 'RADIOS')
REAL                          :: RHO    ! density of launch (if 'ISODEN')
REAL                          :: PRES   ! pressure of launch (if 'ISODEN')
REAL                          :: DIAMETER! apparent diameter of the balloon (m) (if 'CVBALL')
REAL                          :: AERODRAG! aerodynamic drag coefficient of the balloon (if 'CVBALL')
REAL                          :: INDDRAG! induced drag coefficient (i.e. air shifted by the balloon) (if 'CVBALL')
REAL                          :: VOLUME ! volume of the balloon (m3) (if 'CVBALL')
REAL                          :: MASS   ! mass of the balloon (kg) (if 'CVBALL')
!
!* aircraft flight definition
!
INTEGER                       :: SEG      ! number of aircraft flight segments
INTEGER                       :: SEGCURN  ! current flight segment number
REAL                          :: SEGCURT  ! current flight segment time spent
REAL, DIMENSION(:),   POINTER :: SEGLAT  => NULL() ! latitude of flight segment extremities  (LEG+1)
REAL, DIMENSION(:),   POINTER :: SEGLON  => NULL() ! longitude of flight segment extremities (LEG+1)
REAL, DIMENSION(:),   POINTER :: SEGX    => NULL() ! X of flight segment extremities         (LEG+1)
REAL, DIMENSION(:),   POINTER :: SEGY    => NULL() ! Y of flight segment extremities         (LEG+1)
REAL, DIMENSION(:),   POINTER :: SEGP    => NULL() ! pressure of flight segment extremities  (LEG+1)
REAL, DIMENSION(:),   POINTER :: SEGZ    => NULL() ! altitude of flight segment extremities  (LEG+1)
REAL, DIMENSION(:),   POINTER :: SEGTIME => NULL() ! duration of flight segments             (LEG  )
!
!* aircraft altitude type definition
!
LOGICAL                       :: ALTDEF   ! TRUE == altitude given in pressure
!
!* current position of the balloon/aircraft
!
REAL                          :: X_CUR    ! current x
REAL                          :: Y_CUR    ! current y
REAL                          :: Z_CUR    ! current z (if 'RADIOS' or 'AIRCRA' and 'ALTDEF' = T)
REAL                          :: P_CUR    ! current p (if 'AIRCRA' and 'ALTDEF' = F)
!
!* data records
!
REAL, DIMENSION(:),    POINTER :: TIME      => NULL() ! t(n)  (n: recording instants)
REAL, DIMENSION(:),    POINTER :: X         => NULL() ! X(n)
REAL, DIMENSION(:),    POINTER :: Y         => NULL() ! Y(n)
REAL, DIMENSION(:),    POINTER :: Z         => NULL() ! Z(n)
REAL, DIMENSION(:),    POINTER :: XLON      => NULL() ! longitude(n)
REAL, DIMENSION(:),    POINTER :: YLAT      => NULL() ! latitude (n)
REAL, DIMENSION(:),    POINTER :: ZON       => NULL() ! zonal wind(n)
REAL, DIMENSION(:),    POINTER :: MER       => NULL() ! meridian wind(n)
REAL, DIMENSION(:),    POINTER :: W         => NULL() ! w(n)  (air vertical speed)
REAL, DIMENSION(:),    POINTER :: P         => NULL() ! p(n)
REAL, DIMENSION(:),    POINTER :: TKE       => NULL() ! tke(n)
REAL, DIMENSION(:),    POINTER :: TKE_DISS  => NULL() ! tke dissipation rate
REAL, DIMENSION(:),    POINTER :: TH        => NULL() ! th(n)
REAL, DIMENSION(:,:),  POINTER :: R         => NULL() ! r*(n)
REAL, DIMENSION(:,:),  POINTER :: SV        => NULL() ! Sv*(n)
REAL, DIMENSION(:,:),  POINTER :: RTZ       => NULL() ! tot hydrometeor mixing ratio
REAL, DIMENSION(:,:,:),POINTER :: RZ        => NULL() ! water vapour mixing ratio
REAL, DIMENSION(:,:),  POINTER :: FFZ       => NULL() ! horizontal wind
REAL, DIMENSION(:,:),  POINTER :: IWCZ      => NULL() ! ice water content
REAL, DIMENSION(:,:),  POINTER :: LWCZ      => NULL() ! liquid water content
REAL, DIMENSION(:,:),  POINTER :: CIZ       => NULL() ! Ice concentration
REAL, DIMENSION(:,:),  POINTER :: CCZ       => NULL() ! Cloud concentration (LIMA)
REAL, DIMENSION(:,:),  POINTER :: CRZ       => NULL() ! Rain concentration (LIMA)
REAL, DIMENSION(:,:),  POINTER :: CRARE     => NULL() ! cloud radar reflectivity
REAL, DIMENSION(:,:),  POINTER :: CRARE_ATT => NULL() ! attenuated (= more realistic) cloud radar reflectivity
REAL, DIMENSION(:,:),  POINTER :: WZ        => NULL() ! vertical profile of vertical velocity
REAL, DIMENSION(:,:),  POINTER :: ZZ        => NULL() ! vertical profile of mass point altitude (above sea)
REAL, DIMENSION(:,:),  POINTER :: AER       => NULL() ! Extinction at 550 nm
REAL, DIMENSION(:,:),  POINTER :: DST_WL    => NULL() ! Extinction by wavelength
REAL, DIMENSION(:),    POINTER :: ZS        => NULL() ! zs(n)
REAL, DIMENSION(:),    POINTER :: TSRAD     => NULL() ! Ts(n)
REAL, DIMENSION(:,:),  POINTER :: DATIME    => NULL() ! record for diachro
!
REAL, DIMENSION(:)  ,   POINTER :: THW_FLUX => NULL() ! thw_flux(n)
REAL, DIMENSION(:)  ,   POINTER :: RCW_FLUX => NULL() ! rcw_flux(n)
REAL, DIMENSION(:,:),   POINTER :: SVW_FLUX => NULL() ! psw_flux(n)
END TYPE FLYER
REAL :: XLAM_CRAD ! cloud radar wavelength (m)
!
!-------------------------------------------------------------------------------------------
!
LOGICAL     :: LFLYER    ! flag to use aircraft/balloons
!
TYPE(FLYER) :: TBALLOON1 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON2 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON3 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON4 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON5 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON6 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON7 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON8 ! characteristics and records of a balloon
TYPE(FLYER) :: TBALLOON9 ! characteristics and records of a balloon
!
TYPE(FLYER) :: TAIRCRAFT1 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT2 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT3 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT4 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT5 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT6 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT7 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT8 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT9 ! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT10! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT11! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT12! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT13! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT14! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT15! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT16! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT17! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT18! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT19! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT20! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT21! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT22! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT23! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT24! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT25! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT26! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT27! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT28! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT29! characteristics and records of an aircraft
TYPE(FLYER) :: TAIRCRAFT30! characteristics and records of an aircraft
!
END MODULE MODD_AIRCRAFT_BALLOON
