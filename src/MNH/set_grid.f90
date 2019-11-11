!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_GRID
!     ####################
!
INTERFACE
!
      SUBROUTINE SET_GRID(KMI,TPINIFILE,                                      &
                          KIU,KJU,KKU,KIMAX_ll,KJMAX_ll,                      &
                          PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,    &
                          PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,    &
                          KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,            &
                          HLBCX,HLBCY,                                        &
                          PTSTEP,PSEGLEN,                                     &
                          PLONORI,PLATORI,PLON,PLAT,                          &
                          PXHAT,PYHAT,PDXHAT,PDYHAT, PMAP,                    &
                          PZS,PZZ,PZHAT,PZTOP,OSLEVE,PLEN1,PLEN2,PZSMT,       &
                          PJ,                                                 &
                          TPDTMOD,TPDTCUR,KSTOP,                              &
                          KBAK_NUMB,KOUT_NUMB,TPBACKUPN,TPOUTPUTN             )
!
USE MODD_TYPE_DATE
USE MODD_IO_ll, ONLY: TFILEDATA,TOUTBAK
!
INTEGER,                INTENT(IN)  :: KMI       ! Model index
TYPE(TFILEDATA),        INTENT(IN)  :: TPINIFILE !Initial file
INTEGER,                INTENT(IN)  :: KIU       ! Upper dimension in x direction
                                                 ! for sub-domain arrays
INTEGER,                INTENT(IN)  :: KJU       ! Upper dimension in y direction
                                                 ! for sub-domain arrays
INTEGER,                INTENT(IN)  :: KKU       ! Upper dimension in z direction
                                                 ! for domain arrays
INTEGER,                INTENT(IN)  :: KIMAX_ll  !  Dimensions  in x direction
                                                 ! of the physical domain,
INTEGER,                INTENT(IN)  :: KJMAX_ll  !  Dimensions  in y direction
                                                 ! of the physical domain,
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
!
REAL,                   INTENT(IN)  :: PTSTEP    ! time step of model KMI
REAL,                   INTENT(INOUT) :: PSEGLEN ! segment duration (in seconds)
!
REAL,                   INTENT(OUT) :: PLONORI   ! Longitude  of the
                                                 ! Origine point of
                                                 ! conformal projection
REAL,                   INTENT(OUT) :: PLATORI   ! Latitude of the
                                                 ! Origine point of
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
REAL,                   INTENT(OUT) :: PZTOP     ! Model top
LOGICAL,                INTENT(OUT) :: OSLEVE    ! flag for SLEVE coordinate
REAL,                   INTENT(OUT) :: PLEN1     ! Decay scale for smooth topography
REAL,                   INTENT(OUT) :: PLEN2     ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:),   INTENT(OUT) :: PZSMT     ! smooth-orography
!
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTMOD   ! date and time of the model
                                                 ! beginning
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTCUR   ! Current date and time
INTEGER,                INTENT(OUT) :: KSTOP     ! number of time steps for
                                                 ! current segment
INTEGER,POINTER,        INTENT(OUT) :: KBAK_NUMB ! number of backups
INTEGER,POINTER,        INTENT(OUT) :: KOUT_NUMB ! number of outputs
TYPE(TOUTBAK),DIMENSION(:),POINTER,INTENT(OUT) :: TPBACKUPN ! List of backups
TYPE(TOUTBAK),DIMENSION(:),POINTER,INTENT(OUT) :: TPOUTPUTN ! List of outputs
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PJ        ! Jacobian
!
END SUBROUTINE SET_GRID
!
END INTERFACE
!
END MODULE MODI_SET_GRID
!
!
!
!
!
!     #########################################################################
      SUBROUTINE SET_GRID(KMI,TPINIFILE,                                      &
                          KIU,KJU,KKU,KIMAX_ll,KJMAX_ll,                      &
                          PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,    &
                          PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,    &
                          KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,            &
                          HLBCX,HLBCY,                                        &
                          PTSTEP,PSEGLEN,                                     &
                          PLONORI,PLATORI,PLON,PLAT,                          &
                          PXHAT,PYHAT,PDXHAT,PDYHAT, PMAP,                    &
                          PZS,PZZ,PZHAT,PZTOP,OSLEVE,PLEN1,PLEN2,PZSMT,       &
                          PJ,                                                 &
                          TPDTMOD,TPDTCUR,KSTOP,                              &
                          KBAK_NUMB,KOUT_NUMB,TPBACKUPN,TPOUTPUTN             )
!     #########################################################################
!
!!****  *SET_GRID* - routine to set grid variables
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
!!        * The latitude and longitude of the origine points (PLATORI and PLONORI)
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
!!      Book2 of documentation (routine SET_GRID)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/06/94
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
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_BUDGET
USE MODD_CONF
USE MODD_CONF_n
USE MODD_DYN
USE MODD_GRID
USE MODD_IO_ll,   ONLY: TFILEDATA,TOUTBAK
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_OUT_n,   ONLY: OUT_MODEL
USE MODD_PARAMETERS
USE MODD_NESTING
!
USE MODE_FIELD,   ONLY: TFIELDDATA,TFIELDLIST,FIND_FIELD_ID_FROM_MNHNAME
USE MODE_FM
USE MODE_FMREAD
USE MODE_GATHER_ll
USE MODE_GRIDCART
USE MODE_GRIDPROJ
USE MODE_IO_MANAGE_STRUCT
USE MODE_ll
USE MODE_TIME
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
INTEGER,                INTENT(IN)  :: KMI       ! Model index
TYPE(TFILEDATA),        INTENT(IN)  :: TPINIFILE !Initial file
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
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the inner model domain, relative to outer model
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
!
REAL,                   INTENT(IN)  :: PTSTEP    ! time step of model KMI
REAL,                   INTENT(INOUT) :: PSEGLEN ! segment duration (in seconds)
!
REAL,                   INTENT(OUT) :: PLONORI   ! Longitude  of the
                                                 ! Origine point of the
                                                 ! conformal projection
REAL,                   INTENT(OUT) :: PLATORI   ! Latitude of the
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
REAL,                   INTENT(OUT) :: PZTOP     ! Model top
LOGICAL,                INTENT(OUT) :: OSLEVE    ! flag for SLEVE coordinate
REAL,                   INTENT(OUT) :: PLEN1     ! Decay scale for smooth topography
REAL,                   INTENT(OUT) :: PLEN2     ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:),   INTENT(OUT) :: PZSMT     ! smooth-orography
!
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTMOD   ! date and time of the model
                                                 ! beginning
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTCUR   ! Current date and time
INTEGER,                INTENT(OUT) :: KSTOP     ! number of time steps for
                                                 ! current segment
INTEGER,POINTER,        INTENT(OUT) :: KBAK_NUMB ! number of backups
INTEGER,POINTER,        INTENT(OUT) :: KOUT_NUMB ! number of outputs
TYPE(TOUTBAK),DIMENSION(:),POINTER,INTENT(OUT) :: TPBACKUPN ! List of backups
TYPE(TOUTBAK),DIMENSION(:),POINTER,INTENT(OUT) :: TPOUTPUTN ! List of outputs
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PJ        ! Jacobian
!
!*       0.2   declarations of local variables
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
REAL                         :: ZXHATM,ZYHATM    ! coordinates of mass point
REAL                         :: ZLATORI, ZLONORI ! lat and lon of left-bottom point
INTEGER                :: IID,IRESP
CHARACTER (LEN=40)     :: YTITLE                    ! Title for date print
INTEGER                :: ILUOUT                    ! Logical unit number for
                                                    ! output-listing
INTEGER                :: IIUP,IJUP ,ISUP=1         ! size  of working
                                                    ! window arrays,
                                                    ! supp. time steps
!
TYPE(TFIELDDATA)       :: TZFIELD
!-------------------------------------------------------------------------------
!
!*       1.    READ GRID  VARIABLES IN INITIAL FILE
!              ------------------------------------
!
!*       1.1   Spatial grid
!
CALL IO_READ_FIELD(TPINIFILE,'STORAGE_TYPE',CSTORAGE_TYPE,IRESP)
IF (IRESP /= 0) CSTORAGE_TYPE='TT'
!
IF (KMI == 1) THEN
  ! this parameter is also useful in the cartesian to
  ! compute the sun position for the radiation scheme
  CALL IO_READ_FIELD(TPINIFILE,'LON0',XLON0)
  !
  ! this parameter is also useful in the cartesian to
  ! compute the Coriolis parameter
  CALL IO_READ_FIELD(TPINIFILE,'LAT0',XLAT0)
  !
  ! this parameter is also useful in the cartesian to
  ! rotate the simulatin domain
  CALL IO_READ_FIELD(TPINIFILE,'BETA',XBETA)
END IF
!
CALL IO_READ_FIELD(TPINIFILE,'XHAT',PXHAT)
CALL IO_READ_FIELD(TPINIFILE,'YHAT',PYHAT)
!
IF (.NOT.LCARTESIAN) THEN
  CALL IO_READ_FIELD(TPINIFILE,'RPK',XRPK)
  !
  IF ( (TPINIFILE%NMNHVERSION(1)==4 .AND. TPINIFILE%NMNHVERSION(2)>5) .OR. TPINIFILE%NMNHVERSION(1)>4 ) THEN
    CALL IO_READ_FIELD(TPINIFILE,'LONORI',PLONORI)
    CALL IO_READ_FIELD(TPINIFILE,'LATORI',PLATORI)
  !
  ELSE                     
    CALL FIND_FIELD_ID_FROM_MNHNAME('LONORI',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CMNHNAME = 'LONOR'
    CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PLONORI)
    !
    CALL FIND_FIELD_ID_FROM_MNHNAME('LATORI',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%CMNHNAME = 'LATOR'
    CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PLATORI)
    !
    ALLOCATE(ZXHAT_ll(KIMAX_ll+ 2 * JPHEXT),ZYHAT_ll(KJMAX_ll+2 * JPHEXT))
    CALL GATHERALL_FIELD_ll('XX',PXHAT,ZXHAT_ll,IRESP) !//
    CALL GATHERALL_FIELD_ll('YY',PYHAT,ZYHAT_ll,IRESP) !//
    ZXHATM = - 0.5 * (ZXHAT_ll(1)+ZXHAT_ll(2))
    ZYHATM = - 0.5 * (ZYHAT_ll(1)+ZYHAT_ll(2))
    CALL SM_LATLON(PLATORI,PLONORI,ZXHATM,ZYHATM,ZLATORI,ZLONORI)
    DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
    PLATORI = ZLATORI
    PLONORI = ZLONORI
  END IF
  !
END IF

CALL IO_READ_FIELD(TPINIFILE,'ZS',PZS)
CALL IO_READ_FIELD(TPINIFILE,'ZHAT',PZHAT)
CALL IO_READ_FIELD(TPINIFILE,'ZTOP',PZTOP)
!
CALL DEFAULT_SLEVE(OSLEVE,PLEN1,PLEN2)
!
IF ( TPINIFILE%NMNHVERSION(1)<4 .OR. (TPINIFILE%NMNHVERSION(1)==4 .AND. TPINIFILE%NMNHVERSION(2)<=6) ) THEN
  PZSMT  = PZS
  OSLEVE = .FALSE.
ELSE
  CALL IO_READ_FIELD(TPINIFILE,'ZSMT',PZSMT)
  CALL IO_READ_FIELD(TPINIFILE,'SLEVE',OSLEVE)
END IF
!
IF (OSLEVE) THEN
  CALL IO_READ_FIELD(TPINIFILE,'LEN1',PLEN1)
  CALL IO_READ_FIELD(TPINIFILE,'LEN2',PLEN2)
END IF
!
!*       1.2   Temporal grid
!
CALL IO_READ_FIELD(TPINIFILE,'DTMOD',TPDTMOD)
CALL IO_READ_FIELD(TPINIFILE,'DTCUR',TPDTCUR)
!
IF (KMI == 1) THEN
CALL IO_READ_FIELD(TPINIFILE,'DTEXP',TDTEXP)
END IF
!
CALL IO_READ_FIELD(TPINIFILE,'DTSEG',TDTSEG)
!
!-------------------------------------------------------------------------------
!
!*       2.    SET OTHER GRID VARIABLES
!              ------------------------
!
!
!*       2.1    Spatial grid
!
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(PXHAT,PYHAT,PZHAT,PZS,OSLEVE,PLEN1,PLEN2,PZSMT,PDXHAT,PDYHAT,PZZ,PJ)
ELSE
  CALL SM_GRIDPROJ(PXHAT,PYHAT,PZHAT,PZS,OSLEVE,PLEN1,PLEN2,PZSMT,PLATORI,PLONORI, &
                   PMAP,PLAT,PLON,PDXHAT,PDYHAT,PZZ,PJ)
END IF
!
!*       2.2    Temporal grid - segment length
!
IF (CPROGRAM /= 'DIAG  ') TDTSEG = TPDTCUR
ISUP = 1     ! 1 corresponds to a last timestep
   ! to obtain the prognostic and diagnostic fields all along this timestep
!
IF ( KMI == 1) PSEGLEN = PSEGLEN + PTSTEP*ISUP ! needed for the gridnesting case to get
                                               ! the same PSEGLEN for all nested models
KSTOP = NINT(PSEGLEN/PTSTEP)
!
!
!*       2.3    Temporal grid - outputs managment
!
! The output/backups times have been read only by model 1
IF (CPROGRAM == 'MESONH' .AND. KMI == 1) CALL IO_PREPARE_BAKOUT_STRUCT(ISUP,PTSTEP,PSEGLEN)
!
KBAK_NUMB => OUT_MODEL(KMI)%NBAK_NUMB
KOUT_NUMB => OUT_MODEL(KMI)%NOUT_NUMB
TPBACKUPN => OUT_MODEL(KMI)%TBACKUPN
TPOUTPUTN => OUT_MODEL(KMI)%TOUTPUTN
!
!-------------------------------------------------------------------------------
!
!*       3.    PRINT ON OUTPUT-LISTING
!              -----------------------
!
ILUOUT = TLUOUT%NLU
IF  (NVERB >= 10) THEN
  IIUP = SIZE(PXHAT)
  IJUP = SIZE(PYHAT)
  IF(LCARTESIAN) THEN
    WRITE(ILUOUT,FMT=*) 'SET_GRID : No map projection used.'
  ELSE
    IF (XRPK == 1.) THEN
      WRITE(ILUOUT,FMT=*) 'SET_GRID : Polar stereo used.'
    ELSE IF (XRPK == 0.) THEN
      WRITE(ILUOUT,FMT=*) 'SET_GRID : Mercator used.'
    ELSE
      WRITE(ILUOUT,FMT=*) 'SET_GRID : Lambert used, cone factor=',XRPK
      WRITE(ILUOUT,FMT=*) ' SET_GRID : LON0 = ',XLON0,' LAT0 = ',XLAT0, &
       ' RPK = ',XRPK,' BETA = ',XBETA,' LONORI = ',PLONORI,' LATORI = ',PLATORI
    END IF
  END IF
  WRITE(ILUOUT,FMT=*) ' SET_GRID : Some PXHAT values:'
  WRITE(ILUOUT,FMT=*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,FMT=*) PXHAT(1),PXHAT(IIUP/2),PXHAT(IIUP)
!
  WRITE(ILUOUT,FMT=*) ' SET_GRID : Some PYHAT values:'
  WRITE(ILUOUT,FMT=*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,FMT=*) PYHAT(1),PYHAT(IJUP/2),PYHAT(IJUP)
!
  WRITE(ILUOUT,FMT=*) ' SET_GRID : Some PZHAT values:'
  WRITE(ILUOUT,FMT=*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,FMT=*) PZHAT(1),PZHAT(KKU/2),PZHAT(KKU)
!
  WRITE(ILUOUT,FMT=*) ' SET_GRID : Some PZS values:'
  WRITE(ILUOUT,FMT=*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,FMT=*) PZS(1,1),PZS(IIUP/2,IJUP/2),PZS(IIUP,IJUP)
!
  YTITLE='CURRENT DATE AND TIME'
  CALL SM_PRINT_TIME(TPDTCUR,TLUOUT,YTITLE)
END IF
IF (NVERB >= 5) THEN
  YTITLE='DATE AND TIME OF EXPERIMENT BEGINNING'
  CALL SM_PRINT_TIME(TDTEXP,TLUOUT,YTITLE)
  YTITLE='DATE AND TIME OF MODEL BEGINNING'
  CALL SM_PRINT_TIME(TPDTMOD,TLUOUT,YTITLE)
END IF
YTITLE='DATE AND TIME OF SEGMENT BEGINNING'
CALL SM_PRINT_TIME(TDTSEG,TLUOUT,YTITLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_GRID
