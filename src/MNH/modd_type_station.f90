!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_TYPE_STATION
!     ############################
!
!!****  *MODD_STATION* - declaration of stations
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!      the different stations types.
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
!!	P. Tulet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/01/02
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
TYPE STATION
!
!
!* general information
!
!
!* storage monitoring
!
REAL                          :: T_CUR  ! current time since last storage
INTEGER                       :: N_CUR  ! current step of storage
REAL                          :: STEP   ! storage time step
!
!* data records
!
CHARACTER(LEN=8),DIMENSION(:), POINTER   :: NAME=>NULL()   ! station name
CHARACTER(LEN=8),DIMENSION(:), POINTER   :: TYPE=>NULL()   ! station type
REAL, DIMENSION(:),   POINTER :: TIME=>NULL()     ! t(n)  (n: recording instants)
LOGICAL, DIMENSION(:),  POINTER :: ERROR=>NULL()  ! 
REAL, DIMENSION(:),     POINTER :: X=>NULL()      ! X(n)
REAL, DIMENSION(:),     POINTER :: Y=>NULL()      ! Y(n)
REAL, DIMENSION(:),     POINTER :: Z=>NULL()      ! Z(n)
REAL, DIMENSION(:),     POINTER :: LON=>NULL()    ! longitude(n)
REAL, DIMENSION(:),     POINTER :: LAT=>NULL()    ! latitude (n)
REAL, DIMENSION(:,:),   POINTER :: ZON=>NULL()    ! zonal wind(n)
REAL, DIMENSION(:,:),   POINTER :: MER=>NULL()    ! meridian wind(n)
REAL, DIMENSION(:,:),   POINTER :: W=>NULL()      ! w(n)  (air vertical speed)
REAL, DIMENSION(:,:),   POINTER :: P=>NULL()      ! p(n)
REAL, DIMENSION(:,:),   POINTER :: TKE=>NULL()    ! tke(n)
REAL, DIMENSION(:,:),   POINTER :: TH=>NULL()     ! th(n)
REAL, DIMENSION(:,:,:), POINTER :: R=>NULL()      ! r*(n)
REAL, DIMENSION(:,:,:), POINTER :: SV=>NULL()     ! Sv*(n)
REAL, DIMENSION(:),     POINTER :: ZS=>NULL()     ! zs(n)
REAL, DIMENSION(:,:),   POINTER :: TSRAD=>NULL()  ! Ts(n)
REAL, DIMENSION(:,:),   POINTER :: DATIME=>NULL() ! record for diachro
!
REAL, DIMENSION(:,:),   POINTER :: T2M=>NULL()    ! 
REAL, DIMENSION(:,:),   POINTER :: Q2M=>NULL()    ! 
REAL, DIMENSION(:,:),   POINTER :: HU2M=>NULL()   ! 
REAL, DIMENSION(:,:),   POINTER :: ZON10M=>NULL() ! 
REAL, DIMENSION(:,:),   POINTER :: MER10M=>NULL() ! 
REAL, DIMENSION(:,:),   POINTER :: RN=>NULL()     ! 
REAL, DIMENSION(:,:),   POINTER :: H=>NULL()      ! 
REAL, DIMENSION(:,:),   POINTER :: LE=>NULL()     !
REAL, DIMENSION(:,:),   POINTER :: LEI=>NULL()    ! 
REAL, DIMENSION(:,:),   POINTER :: GFLUX=>NULL()  !
REAL, DIMENSION(:,:),   POINTER :: SWD=>NULL()     ! 
REAL, DIMENSION(:,:),   POINTER :: SWU=>NULL()     ! 
REAL, DIMENSION(:,:),   POINTER :: LWD=>NULL()     !
REAL, DIMENSION(:,:),   POINTER :: LWU=>NULL()     !
REAL, DIMENSION(:,:),   POINTER :: SWDIR=>NULL()  ! 
REAL, DIMENSION(:,:),   POINTER :: SWDIFF=>NULL() !
REAL, DIMENSION(:,:),   POINTER :: DSTAOD=>NULL() ! Dust Aerosol Optical Depth
REAL, DIMENSION(:,:),   POINTER :: SFCO2=>NULL()  ! CO2 surface flux
!
INTEGER, DIMENSION(:),  POINTER :: K=>NULL()      ! Model level for altitude
                                                  ! comparisons
INTEGER, DIMENSION(:),  POINTER :: I=>NULL()      ! i index  (n)
INTEGER, DIMENSION(:),  POINTER :: J=>NULL()      ! j index  (n)

END TYPE STATION
!
END MODULE MODD_TYPE_STATION
