!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_LIGHT_GRID
!     ####################
!
INTERFACE
!
      SUBROUTINE SET_LIGHT_GRID(KMI,HINIFILE,HLUOUT,                          &
                          KIU,KJU,KKU,KIMAX_ll,KJMAX_ll,                      &
                          PLONORI,PLATORI,PLON,PLAT,                          &
                          PXHAT,PYHAT,PDXHAT,PDYHAT, PMAP,                    &
                          PZS,PZZ,PZHAT,OSLEVE,PLEN1,PLEN2,PZSMT,             &
                          PJ,                                                 &
                          TPDTMOD,TPDTCUR         )
!
USE MODD_TYPE_DATE
!
INTEGER,                INTENT(IN)  :: KMI       ! Model index 
CHARACTER (LEN=*),      INTENT(IN)  :: HINIFILE  ! Name of the initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
INTEGER,                INTENT(IN)  :: KIU       ! Upper dimension in x direction
                                                 ! for sub-domain arrays  
INTEGER,                INTENT(IN)  :: KJU       ! Upper dimension in y direction
                                                 ! for sub-domain arrays 
INTEGER,                INTENT(IN)  :: KKU       ! Upper dimension in z direction
                                                 ! for domain arrays 
INTEGER,               INTENT(IN)   :: KIMAX_ll  !  Dimensions  in x direction 
                                                 ! of the physical domain,
INTEGER,               INTENT(IN)   :: KJMAX_ll  !  Dimensions  in y direction 
                                                 ! of the physical domain,
!  
REAL,                   INTENT(OUT) :: PLONORI    ! Longitude  of the
                                                  ! Origine point of the  
                                                  ! conformal projection
REAL,                   INTENT(OUT) :: PLATORI    ! Latitude of the
                                                  ! Origine point of the
                                                  ! conformal projection
REAL, DIMENSION(:,:),   INTENT(OUT) :: PLON,PLAT ! Longitude and latitude  
REAL, DIMENSION(:),     INTENT(OUT) :: PXHAT     ! Position x in the conformal
                                                 ! plane or on the cartesian plane
REAL, DIMENSION(:),     INTENT(OUT) :: PYHAT     ! Position y in the conformal
                                                 ! plane or on the cartesian plane
REAL, DIMENSION(:),     INTENT(OUT) :: PDXHAT    ! horizontal stretching in x
REAL, DIMENSION(:),     INTENT(OUT) :: PDYHAT    ! horizontal stretching in y
REAL, DIMENSION(:,:),   INTENT(OUT) :: PMAP      ! Map factor
!
REAL, DIMENSION(:,:),   INTENT(OUT) :: PZS       ! orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PZZ       ! Height z                                           
REAL, DIMENSION(:),     INTENT(OUT) :: PZHAT     ! Height  level   
LOGICAL,                INTENT(OUT) :: OSLEVE    ! flag for SLEVE coordinate
REAL,                   INTENT(OUT) :: PLEN1     ! Decay scale for smooth topography
REAL,                   INTENT(OUT) :: PLEN2     ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:),   INTENT(OUT) :: PZSMT     ! smooth-orography
!
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTMOD   ! date and time of the model
                                                 ! beginning
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTCUR   ! Current date and time 
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PJ        ! Jacobian 
!
END SUBROUTINE SET_LIGHT_GRID
!
END INTERFACE
!
END MODULE MODI_SET_LIGHT_GRID
!
!
!
!
!
!     #########################################################################
      SUBROUTINE SET_LIGHT_GRID(KMI,HINIFILE,HLUOUT,                          &
                          KIU,KJU,KKU,KIMAX_ll,KJMAX_ll,                      &
                          PLONORI,PLATORI,PLON,PLAT,                          &
                          PXHAT,PYHAT,PDXHAT,PDYHAT, PMAP,                    &
                          PZS,PZZ,PZHAT,OSLEVE,PLEN1,PLEN2,PZSMT,             &
                          PJ,                                                 &
                          TPDTMOD,TPDTCUR         )
!     #########################################################################
!
!!****  *SET_LIGHT_GRID* - routine to set grid variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set spatio-temporal grid variables
!
!!**  METHOD
!!    ------
!!
!!      The spatial grid variables are read in initial file : 
!!        * The reference latitude (XLAT0), the reference longitude (XLON0) and
!!      the projection parameter (XPRPK) are read if spherical geometry is used.
!!      (LCARTESIAN=.FALSE.) and only at the first call (by INI_MODEL1,i.e. KMI=1),
!!     since it is the same for all nested models.
!!        * The rotation angle (XBETA) is read only at the first call for the
!!     same reason. 
!!        * The latitude and longitude of the origine points (XLATOR and XLONOR)
!!     are read for a spherical geometry (LCARTESIAN=.FALSE.).
!!        * The horizontal positions (PXHAT and PYHAT) are always read. 
!!
!!      The temporal grid variables are read in initial file : 
!!        * The number of time steps for the current segment depends on the time step
!!     PTSTEP and on the segment length PSEGLEN plus one time step of the first
!!     model for all models. 
!!        * The time of the beginning of experiment (TDTEXP of type DATE_TIME) 
!!     is read only at the first call  by INI_MODEL1 (KMI=1), 
!!     since it is the same for all nested models.
!!        * The times of the  beginning of model (TPDTMOD of type DATE_TIME),
!!     of beginning of segment (TPDTSEG  of type DATE_TIME) are read for
!!     all models
!!
!!      Then, the other spatial grid variables are deduced :
!!        * If Cartesian geometry (LCARTESIAN=.TRUE.), SM_GRIDCART computes 
!!      the horizontal stretchings (PDXHAT and PDYHAT) the height (PZZ) and the 
!!      Jacobian (PJ).
!!        * if Spherical geometry (LCARTESIAN=.FALSE.), SM_GRIDPROJ computes 
!!      the horizontal stretchings (PDXHAT and PDYHAT) the height (PZZ), the 
!!      Jacobian (PJ), the map factor (PMAP), the latitude (PLAT) and the 
!!      longitude (PLON).    
!!
!!      and  the other temporal  grid variables are deduced :
!!        The current time (TPDTCUR of type DATE_TIME) is set equal to the time
!!    of beginning of segment.
!!
!!     IF verbose option (NVERB >=5), the time is printed on output-listing
!!    EXTERNAL
!!    --------   
!!      FMREAD      : to read data in LFIFM file 
!!      FMLOOK      : to retrieve a logical unit number 
!!
!!      Module MODE_GRIDPROJ : contains conformal projection routines 
!!        SM_GRIDPROJ : to compute some grid variables in case of conformal
!!                       projection
!!        SM_LATLON   : to compute latitude and longitude, giving the 
!!                      positions on the grid
!!      Module MODE_GRIDCART : contains  cartesian geometry routines 
!!        SM_GRIDCART : to compute some grid_variables in case of cartesian
!!                       geometry 
!!      Module MODE_TIME : contains SM_PRINT_TIME routine
!!                         and uses module MODD_TIME (for definition
!!                         of types DATE_TIME and DATE
!!       ZS_BOUNDARY   : replace the orography outside the fine-mesh model by
!!                       the large-scale orography of the DAD model
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!      Module MODD_CONF       : contains declaration of configuration variables
!!                              for all models
!!         CCONF      : Configuration for all models ( START, RESTART or POST)
!!         LCARTESIAN :  Logical for cartesian geometry 
!!                       .TRUE.  = cartesian geometry 
!!         NVERB      : Level of informations on output-listing
!!                          0 for minimum  prints
!!                          5 for intermediate level of prints
!!                         10 for maximum  prints 
!!         CSTORAGE_TYPE : type of stored informations ( 2 or one instant)
!! 
!!
!!      Module MODD_GRID       : contains spatial  grid variables for all model
!!
!!         XLON0 : Reference longitude for the conformal projection
!!         XLAT0 : Reference latitude  
!!         XBETA : Rotation angle 
!!         XRPK  : Projection parameter for the conformal projection
!!
!!      Module MODE_TIME      : uses module MODD_TIME (contains temporal grid
!!                            variables for all model
!!                  TDTEXP : Date and time for the experiment beginning
!!                  TDTSEG : Date and time for the segment beginning
!! 
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine SET_LIGHT_GRID)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    /06/94 
!!      J. STEIN    02/01/95  correct the TPDTCUR initialization 
!!      J. STEIN    26/01/95  read TPDTCUR in the FM-file 
!!      J. STEIN    16/03/95  bug in the TPDTCUR reading
!!      J. STEIN    16/04/95  another bug in the TPDTCUR initialization
!!      J. STEIN    03/01/96  change the temporal grid 
!!      J. STEIN P.JABOUILLE 30/04/96 add the storage-type reading
!!      J. STEIN    25/05/96  read RPK only in the non-cartesian case
!!      J.P. LAFORE 03/07/97  gridnesting implementation
!!      V. DUCROCQ   13/08/98  //
!!      J. STEIN    01/02/99  change the orography at the boundary for the
!!                            grid-nesting lbc
!!     V.MASSON 12/10/00 read of the orography in all cases, even if LFLAT=T
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
!
USE MODD_CONF
USE MODD_GRID
USE MODD_TIME
!
USE MODE_GRIDCART
USE MODE_GRIDPROJ
!USE MODE_ll
!USE MODI_GATHER_ll  !!!! a mettre dans mode_ll
!
!USE MODE_FMREAD
USE MODI_FMREAD
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!  
INTEGER,                INTENT(IN)  :: KMI       ! Model index 
CHARACTER (LEN=*),      INTENT(IN)  :: HINIFILE  ! Name of the initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
INTEGER,                INTENT(IN)  :: KIU       ! Upper dimension in x direction
                                                 ! for sub-domain arrays  
INTEGER,                INTENT(IN)  :: KJU       ! Upper dimension in y direction
                                                 ! for sub-domain arrays 
INTEGER,                INTENT(IN)  :: KKU       ! Upper dimension in z direction
                                                 ! for domain arrays 
INTEGER,               INTENT(IN)   :: KIMAX_ll  !  Dimensions  in x direction 
                                                 ! of the physical domain,
INTEGER,               INTENT(IN)   :: KJMAX_ll  !  Dimensions  in y direction 
                                                 ! of the physical domain,
!
REAL,                   INTENT(OUT) :: PLONORI    ! Longitude  of the
                                                  ! Origine point of the  
                                                  ! conformal projection
REAL,                   INTENT(OUT) :: PLATORI    ! Latitude of the
                                                  ! Origine point of the
                                                  ! conformal projection
REAL, DIMENSION(:,:),   INTENT(OUT) :: PLON,PLAT ! Longitude and latitude  
REAL, DIMENSION(:),     INTENT(OUT) :: PXHAT     ! Position x in the conformal
                                                 ! plane or on the cartesian plane
REAL, DIMENSION(:),     INTENT(OUT) :: PYHAT     ! Position y in the conformal
                                                 ! plane or on the cartesian plane
REAL, DIMENSION(:),     INTENT(OUT) :: PDXHAT    ! horizontal stretching in x
REAL, DIMENSION(:),     INTENT(OUT) :: PDYHAT    ! horizontal stretching in y
REAL, DIMENSION(:,:),   INTENT(OUT) :: PMAP      ! Map factor
!
REAL, DIMENSION(:,:),   INTENT(OUT) :: PZS       ! orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PZZ       ! Height z                                           
REAL, DIMENSION(:),     INTENT(OUT) :: PZHAT     ! Height  level   
LOGICAL,                INTENT(OUT) :: OSLEVE    ! flag for SLEVE coordinate
REAL,                   INTENT(OUT) :: PLEN1     ! Decay scale for smooth topography
REAL,                   INTENT(OUT) :: PLEN2     ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:),   INTENT(OUT) :: PZSMT     ! smooth-orography
!
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTMOD   ! date and time of the model
                                                 ! beginning
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTCUR   ! Current date and time 
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PJ        ! Jacobian 
!  
!*       0.2   declarations of local variables
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
REAL                   :: ZXHATM,ZYHATM    ! coordinates of mass point 
REAL                   :: ZLONORI,ZLATORI  ! lon/lat of mass point (x=0,y=0)
INTEGER                :: ILENG,IGRID,ILENCH,IRESP  !   File 
CHARACTER (LEN=16)     :: YRECFM              ! management
CHARACTER (LEN=100)    :: YCOMMENT            ! variables  
!CHARACTER (LEN=2)      :: YDIR                !
INTEGER, DIMENSION(3)  :: ITDATE           ! date array
INTEGER                :: IMASDEV                   ! masdev of the file
LOGICAL                :: GSLEVE    ! local flag for SLEVE coordinate
!
!-------------------------------------------------------------------------------
!
YRECFM='MASDEV' 
!YDIR='--'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,IMASDEV,IGRID,ILENCH,YCOMMENT,IRESP)
IF (IRESP /=0) IMASDEV=43
!
!*       1.    READ GRID  VARIABLES IN INITIAL FILE
!              ------------------------------------
!
!*       1.1   Spatial grid
!
IF (KMI == 1) THEN
  YRECFM='STORAGE_TYPE' 
  !YDIR='--'
  ILENG=2
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,CSTORAGE_TYPE,IGRID,ILENCH,YCOMMENT,IRESP)
  IF (IRESP /= 0) CSTORAGE_TYPE='MT'
  !
  YRECFM='LON0'     ! this parameter is also useful in the cartesian to
  !YDIR='--'        ! compute the sun position for the radiation scheme
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XLON0,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LAT0'     ! this parameter is also useful in the cartesian to 
  !YDIR='--'        ! compute the Coriolis parameter
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XLAT0,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='BETA'     ! this parameter is also useful in the cartesian to 
  !YDIR='--'           ! rotate the simulatin domain
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XBETA,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
YRECFM='XHAT'
!YDIR='XX'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,KIU,PXHAT,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='YHAT'
!YDIR='YY'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,KJU,PYHAT,IGRID,ILENCH,YCOMMENT,IRESP)
!
IF (.NOT.LCARTESIAN) THEN
  YRECFM='RPK'
  !YDIR='--'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XRPK,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LONORI'
  !YDIR='--'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLONORI,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LATORI'
  !YDIR='--'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLATORI,IGRID,ILENCH,YCOMMENT,IRESP)
!  compute  PLATORI,PLONORI i.e. latitude and longitude of
!  coordinates x=0, y=0 of the grid.
  IF (IMASDEV<=45) THEN
!! compute  PLATOR,PLONOR of each sub-domain
!! i.e. latitude and longitude of mass point (1,1)
  !IF (NPROC > 1) THEN
  !  ALLOCATE(ZXHAT_ll(KIMAX_ll+ 2 * JPHEXT),ZYHAT_ll(KJMAX_ll+2 * JPHEXT))
  !  CALL GATHERALL_FIELD_ll('XX',PXHAT,ZXHAT_ll,IRESP) !//
  !  CALL GATHERALL_FIELD_ll('YY',PYHAT,ZYHAT_ll,IRESP) !//
  !  ZXHATM =0.5 * (PXHAT(1)+PXHAT(2))
  !  ZYHATM =0.5 * (PYHAT(1)+PYHAT(2))
  !  CALL SM_LATLON(ZXHAT_ll,ZYHAT_ll,PLATOR_ll,PLONOR_ll,ZXHATM,ZYHATM,&
  !       PLATOR,PLONOR)
  !  DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
  !ELSE
  ! PLATOR = PLATOR_ll
  ! PLONOR = PLONOR_ll
  !END IF 
  YRECFM='LONOR'
  !YDIR='--'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLONORI,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LATOR'
  !YDIR='--'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLATORI,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  ZXHATM = - 0.5 * (PXHAT(1)+PXHAT(2))
  ZYHATM = - 0.5 * (PYHAT(1)+PYHAT(2))
  CALL SM_LATLON(PLATORI,PLONORI,ZXHATM,ZYHATM,&
                   ZLATORI,ZLONORI)
  PLATORI = ZLATORI
  PLONORI = ZLONORI
  END IF
  !
END IF
!
YRECFM='ZS'
!YDIR='XY'
ILENG=KIU*KJU
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PZS,IGRID,ILENCH,YCOMMENT,IRESP)
IF (IRESP /= 0)THEN
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG/3,PZS(:,2),IGRID,ILENCH,YCOMMENT,IRESP)
  IF(IRESP == 0)THEN
    PZS(:,1)=PZS(:,2)
    PZS(:,3)=PZS(:,2)
  ELSE
    PZS(:,:) = 0.
  ENDIF
ENDIF
!
YRECFM='ZHAT'
!YDIR='--'
ILENG=KKU
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PZHAT,IGRID,ILENCH,YCOMMENT,IRESP)
!
!CALL DEFAULT_SLEVE(OSLEVE,PLEN1,PLEN2)
OSLEVE=.FALSE.
PLEN1=7500.
PLEN2=2500.
!
IF (IMASDEV<=46) THEN
  PZSMT  = PZS
  OSLEVE = .FALSE.
ELSE
  YRECFM='SLEVE'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,GSLEVE,IGRID,ILENCH,YCOMMENT,IRESP)
  IF (IRESP ==0) OSLEVE=GSLEVE
  !
  YRECFM='ZSMT'
  ILENG=KIU*KJU
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PZSMT,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
IF (OSLEVE) THEN
  YRECFM='LEN1'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLEN1,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LEN2'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLEN2,IGRID,ILENCH,YCOMMENT,IRESP)
  print *,'set_light_grid: SLEVE=',OSLEVE,PLEN1,PLEN2
END IF
!
!*       1.2   Temporal grid
!
IF (KMI == 1) THEN
  YRECFM='DTEXP%TDATE' 
  !YDIR='--'
  ILENG=3
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
  TDTEXP%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3))  
  YRECFM='DTEXP%TIME'
  !YDIR='--'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TDTEXP%TIME,IGRID,ILENCH,           &
             YCOMMENT,IRESP)
END IF 
!   
YRECFM='DTCUR%TDATE' 
!YDIR='--'
ILENG=3
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
TPDTCUR%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3)) 
YRECFM='DTCUR%TIME'
!YDIR='--'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TPDTCUR%TIME,IGRID,ILENCH,           &
            YCOMMENT,IRESP) 
!
YRECFM='DTMOD%TDATE' 
!YDIR='--'
ILENG=3
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
TPDTMOD%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3)) 
YRECFM='DTMOD%TIME'
!YDIR='--'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TPDTMOD%TIME,IGRID,ILENCH,           &
            YCOMMENT,IRESP)
!
YRECFM='DTSEG%TDATE' 
!YDIR='--'
ILENG=3
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
TDTSEG%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3)) 
YRECFM='DTSEG%TIME'
!YDIR='--'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TDTSEG%TIME,IGRID,ILENCH,           &
            YCOMMENT,IRESP)
!
!-------------------------------------------------------------------------------
!
!*       2.    SET OTHER GRID VARIABLES 
!              ------------------------
!
!*       2.1    Spatial grid
! 
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(HLUOUT,PXHAT,PYHAT,PZHAT,PZS,OSLEVE,PLEN1,PLEN2,PZSMT,PDXHAT,PDYHAT,PZZ,PJ) 
ELSE
  CALL SM_GRIDPROJ(HLUOUT,PXHAT,PYHAT,PZHAT,PZS,OSLEVE,PLEN1,PLEN2,PZSMT,PLATORI,PLONORI, &
                   PMAP,PLAT,PLON,PDXHAT,PDYHAT,PZZ,PJ)
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_LIGHT_GRID
