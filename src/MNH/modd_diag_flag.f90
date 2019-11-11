!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/modd_diag_flag.f90,v $ $Revision: 1.2.4.1.2.2.10.2.2.2.2.2 $
! masdev4_8 modd 2008/06/30 15:13:13
!-----------------------------------------------------------------
!     #####################
      MODULE MODD_DIAG_FLAG
!     ######################
!
!!****  *MODD_DIAG_FLAG* - declaration of flags related to the diagnostic
!!                          fields computed in diag program
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       16/11/98     
!!       J.-P. Pinty   29/11/02 add C3R5, ICE2, ICE4, ELEC
!!       J.-P. Chaboureau 15/04/03  add LRAD_SUBG_COND
!!       L. Leriche 21/04/06 add aqueous phase chemistry LCHAQDIAG
!!       D.Ricard 2015 : add LMOIST_ES
!!     C.Lac 10/2016  Add visibility diagnostic
!!  10/2016     (F Brosse) Add prod/loss terms computation for chemistry  
!! 10/2017      (G.Delautier) New boundary layer height : replace LBLTOP by CBLTOP 
!!       T. Dauhut     10/2017 add parallel 3D clustering
!!       J.-P. Chaboureau 01/2018 add altitude interpolation
!!       J.-P. Chaboureau 01/2018 add coarse graining
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER :: NVERSION_RAD ! version of radar diagnostics routine 
                        ! 1 : RADAR_RAIN_ICE from JPP (original version)
                        ! 2 : RADAR_SIMULATOR from OC 
!
CHARACTER(LEN=6)  :: CISO       ! PABSM,POVOM,THM for isoPR,EV,TK_level display
LOGICAL     :: LVAR_RS              ! UM,VM,WM,RVM
LOGICAL     :: LVAR_LS              ! LSUM,LSVM,LSWM,LSRVM
INTEGER     :: NCONV_KF             ! Convective scheme
INTEGER     :: NRAD_3D              ! Radiative scheme
CHARACTER(LEN=44)  :: CRAD_SAT      ! GOES-E,GOES-W,GMS,INDSAT,METEOSAT
LOGICAL     :: LRAD_SUBG_COND       ! to activate subgrid condensation
                                    !scheme in the radiatif transfer code
!rttov and satellites variables
INTEGER, DIMENSION(4,10) :: NRTTOVINFO  ! 4-column table where c1=pltform,
! c2=sat c3=sensor c4=choise !for calculations (tb, jacobian, adjoint)
LOGICAL     :: LVAR_TURB
LOGICAL     :: LTURBFLX             ! Turbulent fluxes
LOGICAL     :: LTURBDIAG            ! Diagnostics of turbulent quantities
LOGICAL     :: LMFFLX               ! Mass-Flux fluxes
REAL        :: XDTSTEP              ! time step when initial file ('TT')
LOGICAL     :: LVAR_MRW
LOGICAL     :: LVAR_MRSV
LOGICAL     :: LVAR_FRC
LOGICAL     :: LTPZH
LOGICAL     :: LMOIST_V
LOGICAL     :: LMOIST_E
LOGICAL     :: LMOIST_ES
LOGICAL     :: LMOIST_S1 ! Diagnostics for Thetas 
LOGICAL     :: LMOIST_S2 ! 
LOGICAL     :: LMOIST_L  !
LOGICAL     :: LCOREF
LOGICAL     :: LVORT, LDIV 
LOGICAL     :: LMEAN_POVO
REAL, DIMENSION(2) :: XMEAN_POVO
LOGICAL     :: LGEO        ! Geostrophic wind (m/s)
LOGICAL     :: LAGEO       ! Ageostrophic wind (m/s)
LOGICAL     :: LWIND_ZM    ! Zonal and meridien components of wind (m/s)
LOGICAL     :: LWIND_CONTRAV ! Contravariant component of wind (m/s)
LOGICAL     :: LMSLP       ! Mean Sea Level Pressure (hPa)
LOGICAL     :: LTHW
LOGICAL     :: LCLD_COV
LOGICAL     :: LHU_FLX
LOGICAL     :: LVAR_PR
LOGICAL     :: LTOTAL_PR
LOGICAL     :: LMEAN_PR
REAL, DIMENSION(2) :: XMEAN_PR
INTEGER     :: NCAPE       ! CAPE, DCAPE, CIN, CAPEMAX, CINMAX
LOGICAL     :: LBV_FR
LOGICAL     :: LRADAR
CHARACTER (LEN=5)  :: CBLTOP
LOGICAL     :: LVISI
LOGICAL     :: LTRAJ       ! to compute trajectories
LOGICAL     :: LCHEMDIAG = .FALSE.  ! flag for chemistry
CHARACTER (LEN=4) :: CAERDIAG  ! aerosols optical thickness type
LOGICAL     :: LCHAQDIAG   ! flag for aqueous phase chemistry
REAL, DIMENSION(10)  :: XCHEMLAT,XCHEMLON ! positions of vertical profiles written by routine write_ts1d
CHARACTER (LEN=1024) :: CSPEC_BU_DIAG
CHARACTER (LEN=1024) :: CSPEC_DIAG
LOGICAL     :: LAIRCRAFT_BALLOON    ! aircraft and balloon trajectories
INTEGER     :: NTIME_AIRCRAFT_BALLOON ! time in seconds of trajectories computing
REAL        :: XSTEP_AIRCRAFT_BALLOON ! minimum time step for trajectories calculations (s)
REAL, DIMENSION(9) :: XLAT_BALLOON  ! initial latitudes of the balloons
                                    !(at file time minus NTIME_AIRCRAFT_BALLOON/2)
REAL, DIMENSION(9) :: XLON_BALLOON  ! initial longitudes of the balloons
REAL, DIMENSION(9) :: XALT_BALLOON  ! initial altitude of the balloons (m)
LOGICAL     :: LC2R2
LOGICAL     :: LC3R5
LOGICAL     :: LELECDIAG            ! flag for atmospheric electricity
!
INTEGER                :: NGPS  ! GPS=0 : ZTD, GPS=1 : ZTD,ZHD,ZWD
CHARACTER (LEN=10), DIMENSION(50) :: CNAM_GPS ! name of the GPS stations
REAL, DIMENSION(50)    :: XLAT_GPS            ! latitude of the GPS stations
REAL, DIMENSION(50)    :: XLON_GPS            ! longitude of the GPS stations
REAL, DIMENSION(50)    :: XZS_GPS             ! height of the GPS stations
REAL                   :: XDIFFORO            ! maximum difference between model
                                !orography and station height accepted 
!
LOGICAL, DIMENSION(100) :: LDIAG
REAL,    DIMENSION(100) :: XDIAG
!
LOGICAL           :: LLIDAR ! flag for lidar computation
CHARACTER (LEN=5) :: CVIEW_LIDAR ! 'NADIR' or 'ZENIT'
REAL              :: XALT_LIDAR  ! altitude of the lidar source 
REAL              :: XWVL_LIDAR  ! wavelength of the lidar source 
!
LOGICAL           :: LISOPR ! flag to write on isobaric level
REAL,DIMENSION(10):: XISOPR ! list of level for isobaric interpolation
LOGICAL           :: LISOTH ! flag to write on isentropic level
REAL,DIMENSION(10):: XISOTH ! list of level for isentropic interpolation
LOGICAL           :: LISOAL ! flag to write on altitude level
REAL,DIMENSION(99):: XISOAL ! list of level for altitude interpolation
!
LOGICAL           :: LCOARSE   ! flag for coarse graining
INTEGER           :: NDXCOARSE ! gridpoint number for coarse graining
!
LOGICAL           :: LCLSTR ! flag for 3D clustering
LOGICAL           :: LBOTUP ! to propagate clustering from bottom to top, otherwise top to bottom
CHARACTER(LEN=8)  :: CFIELD ! field on which clustering is applied, could be 'W' or 'CLOUD'
REAL              :: XTHRES ! threshold value to detect the 3D structures

END MODULE MODD_DIAG_FLAG
