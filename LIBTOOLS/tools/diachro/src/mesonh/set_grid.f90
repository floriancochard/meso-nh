!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/init/s.set_grid.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_GRID
!     ####################
!
INTERFACE
!
      SUBROUTINE SET_GRID(KMI,HINIFILE,HLUOUT,                                &
                          KIU,KJU,KKU,KIINF,KISUP,KJINF,KJSUP,                &
                          PTSTEP,PSEGLEN,                                     &
                          POUT1,POUT2,POUT3,POUT4,POUT5,POUT6,POUT7,POUT8,    &
                          POUT9,POUT10,POUT11,POUT12,POUT13,POUT14,POUT15,    &
                          POUT16,POUT17,POUT18,POUT19,POUT20,                 &
                          PLONOR,PLATOR,PLON,PLAT,                            &
                          PXHAT,PYHAT,PDXHAT,PDYHAT, PMAP,                    &
                          PZS,PZZ,PZHAT,                                      &
                          PJ,                                                 &
                          TPDTMOD,TPDTCUR,KSTOP,KOUT_TIMES,KOUT_NUMB)
!
USE MODE_TIME
!
INTEGER,                INTENT(IN)  :: KMI       ! Model index 
CHARACTER (LEN=*),      INTENT(IN)  :: HINIFILE  ! Name of the initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
INTEGER,                INTENT(IN)  :: KIU       ! Upper dimension in x direction
                                                 ! for arrays in initial file  
INTEGER,                INTENT(IN)  :: KJU       ! Upper dimension in y direction
                                                 ! for arrays in initial file  
INTEGER,                INTENT(IN)  :: KKU       ! Upper dimension in z direction
                                                 ! for arrays in initial file
INTEGER,                INTENT(IN)  :: KIINF,KISUP   
                                                 ! Lower and upper  dimensions
                                                 ! in x direction for working 
                                                 ! window  
INTEGER,                INTENT(IN)  :: KJINF,KJSUP 
                                                 ! Lower and upper dimensions
                                                 !  in y direction for working
                                                 ! window
REAL,                   INTENT(IN)  :: PTSTEP    ! time step of model KMI
REAL,                   INTENT(INOUT) :: PSEGLEN ! segment duration (in seconds)
REAL, INTENT(INOUT)  ::  POUT1,POUT2,POUT3,POUT4,POUT5,POUT6,POUT7,POUT8
REAL, INTENT(INOUT)  ::  POUT9,POUT10,POUT11,POUT12,POUT13,POUT14,POUT15
REAL, INTENT(INOUT)  ::  POUT16,POUT17,POUT18,POUT19,POUT20
! increments in seconds from the beginning of the segment to the 
! instant where the n-th fields output on FM-files is realized
!  
REAL,                   INTENT(OUT) :: PLONOR    ! Longitude  of the
                                                 ! Origine point for the 
                                                 ! conformal projection
REAL,                   INTENT(OUT) :: PLATOR    ! Latitude of the
                                                 ! Origine point for the 
                                                 ! conformal projectio
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
!
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTMOD   ! date and time of the model
                                                 ! beginning
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTCUR   ! Current date and time 
INTEGER,                INTENT(OUT) :: KSTOP     ! number of time steps for
                                                 ! current segment 
INTEGER, DIMENSION(:), INTENT(OUT)  :: KOUT_TIMES ! list of the values
               ! of the temporal index in the temporal model loop where fields
               !  outputs on FM-files are realized
INTEGER,                INTENT(OUT) :: KOUT_NUMB ! number of outputs
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
!     #########################################################################
      SUBROUTINE SET_GRID(KMI,HINIFILE,HLUOUT,                                &
                          KIU,KJU,KKU,KIINF,KISUP,KJINF,KJSUP,                &
                          PTSTEP,PSEGLEN,                                     &
                          POUT1,POUT2,POUT3,POUT4,POUT5,POUT6,POUT7,POUT8,    &
                          POUT9,POUT10,POUT11,POUT12,POUT13,POUT14,POUT15,    &
                          POUT16,POUT17,POUT18,POUT19,POUT20,                 &
                          PLONOR,PLATOR,PLON,PLAT,                            &
                          PXHAT,PYHAT,PDXHAT,PDYHAT, PMAP,                    &
                          PZS,PZZ,PZHAT,                                      &
                          PJ,                                                 &
                          TPDTMOD,TPDTCUR,KSTOP,KOUT_TIMES,KOUT_NUMB)
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
!!        * The latitude and longitude of the origine points (XLATOR and XLONOR)
!!     are read for a spherical geometry (LCARTESIAN=.FALSE.).
!!        * The horizontal positions (PXHAT and PYHAT) are always read. 
!!        * The orography (PZS) is set equal to zero if zero orography is needed 
!!     (LFLAT=.TRUE.), else it is  read in initial file.
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
!!         LFLAT      : Logical for zero ororography
!!                       .TRUE.  = no orography (zs=0.)
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
!!	V. Ducrocq       * Meteo France *
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_GRID
USE MODD_BUDGET
USE MODD_DYN
USE MODD_FMOUT
USE MODD_NESTING
!
USE MODE_GRIDCART
USE MODE_GRIDPROJ
USE MODE_TIME
!
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
                                                 ! for arrays in initial file  
INTEGER,                INTENT(IN)  :: KJU       ! Upper dimension in y direction
                                                 ! for arrays in initial file  
INTEGER,                INTENT(IN)  :: KKU       ! Upper dimension in z direction
                                                 ! for arrays in initial file
INTEGER,                INTENT(IN)  :: KIINF,KISUP   
                                                 ! Lower and upper  dimensions
                                                 ! in x direction for working 
                                                 ! window  
INTEGER,                INTENT(IN)  :: KJINF,KJSUP 
                                                 ! Lower and upper dimensions
                                                 !  in y direction for working
                                                 ! window
REAL,                   INTENT(IN)  :: PTSTEP    ! time step of model KMI
REAL,                   INTENT(INOUT) :: PSEGLEN ! segment duration (in seconds)
REAL, INTENT(INOUT)  ::  POUT1,POUT2,POUT3,POUT4,POUT5,POUT6,POUT7,POUT8
REAL, INTENT(INOUT)  ::  POUT9,POUT10,POUT11,POUT12,POUT13,POUT14,POUT15
REAL, INTENT(INOUT)  ::  POUT16,POUT17,POUT18,POUT19,POUT20
! increments in seconds from the beginning of the segment to the 
! instant where the n-th fields output on FM-files is realized
!  
REAL,                   INTENT(OUT) :: PLONOR    ! Longitude  of the
                                                 ! Origine point for the 
                                                 ! conformal projection
REAL,                   INTENT(OUT) :: PLATOR    ! Latitude of the
                                                 ! Origine point for the 
                                                 ! conformal projectio
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
!
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTMOD   ! date and time of the model
                                                 ! beginning
TYPE (DATE_TIME),       INTENT(OUT) :: TPDTCUR   ! Current date and time 
INTEGER,                INTENT(OUT) :: KSTOP     ! number of time steps for
                                                 ! current segment 
INTEGER, DIMENSION(:), INTENT(OUT)  :: KOUT_TIMES ! list of the values
               ! of the temporal index in the temporal model loop where fields
               !  outputs on FM-files are realized
INTEGER,                INTENT(OUT) :: KOUT_NUMB ! number of outputs
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PJ        ! Jacobian 
!  
!*       0.2   declarations of local variables
!
REAL, DIMENSION(KIU)         :: Z1DI              ! 1D array (x direction) used
                                                  ! to read data in inital file
REAL, DIMENSION(KJU)         :: Z1DJ              ! 1D array (y direction) used
                                                  ! to read data in inital file
REAL                         :: ZXHATM,ZYHATM     ! coordinates of mass point 
                                                  ! (KIINF,KJINF)
REAL                         :: ZLATORNEW,ZLONORNEW ! geographical coordinates 
                                                  ! of mass point (KIINF,KJINF)
REAL, DIMENSION(KIU,KJU)     :: Z2D               ! 2D array (x,y directions) used
                                                  ! to read data in inital file
INTEGER                      :: I2D               ! size of 2D arrays
INTEGER                :: ILENG,IGRID,ILENCH,IRESP  !   File 
CHARACTER (LEN=16)     :: YRECFM                    ! management
CHARACTER (LEN=100)    :: YCOMMENT                  ! variables  
INTEGER, DIMENSION(3)  :: ITDATE           ! date array
CHARACTER (LEN=40)     :: YTITLE                    ! Title for date print 
INTEGER                :: ILUOUT                    ! Logical unit number for
                                                    ! output-listing
INTEGER                :: JKLOOP,JOUT               ! Loop index
INTEGER                :: IIUP,IJUP ,ISUP=1         ! size  of working 
                                                    ! window arrays, 
                                                    ! supp. time steps
INTEGER, DIMENSION(2)  :: ISTORAGE_TYPE             ! integer values of the 
                                                    ! ASCII codes for CSTORAGE_TYPE 
!
!-------------------------------------------------------------------------------
!
!*       1.    READ GRID  VARIABLES IN INITIAL FILE
!              ------------------------------------
!
!*       1.1   Spatial grid
!
IF (KMI == 1) THEN
  YRECFM='STORAGE_TYPE' 
  ILENG=2
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ISTORAGE_TYPE,IGRID,ILENCH,YCOMMENT,IRESP)
  IF (IRESP == 0) THEN
    CSTORAGE_TYPE(1:1)=ACHAR(ISTORAGE_TYPE(1))
    CSTORAGE_TYPE(2:2)=ACHAR(ISTORAGE_TYPE(2))
  ELSE
    CSTORAGE_TYPE='MT'
  END IF
  !
  YRECFM='LON0'     ! this parameter is also useful in the cartesian to
  ILENG=1           ! compute the sun position for the radiation scheme
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XLON0,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LAT0'     ! this parameter is also useful in the cartesian to 
  ILENG=1           ! compute the Coriolis parameter
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XLAT0,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='BETA'     ! this parameter is also useful in the cartesian to 
  ILENG=1           ! rotate the simulatin domain
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XBETA,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
!
IF (.NOT.LCARTESIAN) THEN
  !
  YRECFM='RPK'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,XRPK,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LONOR'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLONOR,IGRID,ILENCH,YCOMMENT,IRESP)
  !
  YRECFM='LATOR'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,PLATOR,IGRID,ILENCH,YCOMMENT,IRESP)
END IF
YRECFM='XHAT'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,KIU,Z1DI,IGRID,ILENCH,YCOMMENT,IRESP)
PXHAT(:)=Z1DI(KIINF:KISUP)
YRECFM='YHAT'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,KJU,Z1DJ,IGRID,ILENCH,YCOMMENT,IRESP)
PYHAT(:)=Z1DJ(KJINF:KJSUP) 
!
! in case of postprocessing working window, compute new PLATOR,PLONOR
!  i.e. latitude and longitude of mass point (KIINF,KJINF)
IF (.NOT.LCARTESIAN) THEN
  IF ((KIINF /= 1).OR.(KJINF /= 1)) THEN
    ZXHATM =0.5 * (PXHAT(1)+PXHAT(2))
    ZYHATM =0.5 * (PYHAT(1)+PYHAT(2))
    CALL SM_LATLON(Z1DI,Z1DJ,PLATOR,PLONOR,ZXHATM,ZYHATM,ZLATORNEW,ZLONORNEW)
    PLATOR = ZLATORNEW
    PLONOR = ZLONORNEW
  END IF 
END IF
!
IF (LFLAT) THEN
  PZS(:,:) = 0.
ELSE
  YRECFM='ZS'
  I2D=KIU*KJU 
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,I2D,Z2D,IGRID,ILENCH,YCOMMENT,IRESP)
  IF(IRESP /= 0)THEN
    CALL FMREAD(HINIFILE,YRECFM,HLUOUT,I2D/3,Z2D(:,2),IGRID,ILENCH,YCOMMENT,IRESP)
    IF(IRESP == 0)THEN
      Z2D(:,1)=Z2D(:,2)
      Z2D(:,3)=Z2D(:,2)
      PZS(:,:) = Z2D(KIINF:KISUP,KJINF:KJSUP)
    ELSE
      PZS(:,:) = 0.
    ENDIF
  ELSE
!   print *,' SET_GRID KIINF,KISUP,KJINF,KJSUP ',KIINF,KISUP,KJINF,KJSUP
!   print *,' SET_GRID size Z2D et PZS ',size(Z2D,1),size(Z2D,2),size(PZS,1),size(PZS,2)
    PZS(:,:) = Z2D(KIINF:KISUP,KJINF:KJSUP)
  ENDIF
END IF
YRECFM='ZHAT'
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,KKU,PZHAT,IGRID,ILENCH,YCOMMENT,IRESP)
!
!*       1.2   Temporal grid
!
IF (KMI == 1) THEN
  YRECFM='DTEXP%TDATE' 
  ILENG=3
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
  TDTEXP%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3))  
  YRECFM='DTEXP%TIME'
  ILENG=1
  CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TDTEXP%TIME,IGRID,ILENCH,           &
             YCOMMENT,IRESP)
END IF 
!   
YRECFM='DTCUR%TDATE' 
ILENG=3
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
TPDTCUR%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3)) 
YRECFM='DTCUR%TIME'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TPDTCUR%TIME,IGRID,ILENCH,           &
            YCOMMENT,IRESP) 
!
YRECFM='DTMOD%TDATE' 
ILENG=3
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
TPDTMOD%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3)) 
YRECFM='DTMOD%TIME'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TPDTMOD%TIME,IGRID,ILENCH,           &
            YCOMMENT,IRESP)
!
YRECFM='DTSEG%TDATE' 
ILENG=3
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,ITDATE,IGRID,ILENCH,YCOMMENT,IRESP)
TDTSEG%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3)) 
YRECFM='DTSEG%TIME'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,TDTSEG%TIME,IGRID,ILENCH,           &
            YCOMMENT,IRESP)
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
  CALL SM_GRIDCART(HLUOUT,PXHAT,PYHAT,PZHAT,PZS,PDXHAT,PDYHAT,PZZ,PJ) 
ELSE
  CALL SM_GRIDPROJ(HLUOUT,PXHAT,PYHAT,PZHAT,PZS,PLATOR,PLONOR, &
                   PMAP,PLAT,PLON,PDXHAT,PDYHAT,PZZ,PJ)
END IF    
!
!*       2.2    Temporal grid - segment length
!
TDTSEG = TPDTCUR 
ISUP = 1     ! 1 corresponds to a last timestep 
   ! to obtain the prognostic and diagnostic fields all along this timestep
!
KOUT_TIMES(:) = -999
!
IF ( KMI == 1) PSEGLEN = PSEGLEN + PTSTEP*ISUP ! needed for the gridnesting case to get
                                               ! the same PSEGLEN for all nested models
KSTOP = NINT(PSEGLEN/PTSTEP)
!
!
!*       2.3    Temporal grid - outputs managment 
!
!*       2.3.1  a) synchronization between nested models through XFMOUT arrays (MODD_FMOUT)
!
DO JOUT = 1,20 
  IF (XFMOUT(KMI,JOUT) /= -999.) THEN
    XFMOUT(KMI,JOUT) = NINT(XFMOUT(KMI,JOUT)/PTSTEP) * PTSTEP
    DO JKLOOP = KMI,JPMODELMAX
      XFMOUT(JKLOOP,JOUT) = XFMOUT(KMI,JOUT)
    END DO
  END IF
END DO
!
!*       2.3.2  b) back to original XOUT variables (MODD_OUTn)
!
POUT1  = XFMOUT(KMI,1 )
POUT2  = XFMOUT(KMI,2 )
POUT3  = XFMOUT(KMI,3 )
POUT4  = XFMOUT(KMI,4 )
POUT5  = XFMOUT(KMI,5 )
POUT6  = XFMOUT(KMI,6 )
POUT7  = XFMOUT(KMI,7 )
POUT8  = XFMOUT(KMI,8 )
POUT9  = XFMOUT(KMI,9 )
POUT10 = XFMOUT(KMI,10)
POUT11 = XFMOUT(KMI,11)
POUT12 = XFMOUT(KMI,12)
POUT13 = XFMOUT(KMI,13)
POUT14 = XFMOUT(KMI,14)
POUT15 = XFMOUT(KMI,15)
POUT16 = XFMOUT(KMI,16)
POUT17 = XFMOUT(KMI,17)
POUT18 = XFMOUT(KMI,18)
POUT19 = XFMOUT(KMI,19)
POUT20 = XFMOUT(KMI,20)
!
KOUT_NUMB =0 
!
IF(POUT1 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT1/PTSTEP) + 1
END IF
!
IF(POUT2 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT2/PTSTEP) + 1
END IF
!
IF(POUT3 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT3/PTSTEP) + 1
END IF
!
IF(POUT4 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT4/PTSTEP) + 1
END IF
!
IF(POUT5 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT5/PTSTEP) + 1
END IF
!
IF(POUT6 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT6/PTSTEP) + 1
END IF
!
IF(POUT7 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT7/PTSTEP) + 1
END IF
!
IF(POUT8 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT8/PTSTEP) + 1
END IF
!
IF(POUT9 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT9/PTSTEP) + 1
END IF
!
IF(POUT10 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT10/PTSTEP) + 1
END IF
!
IF(POUT11 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT11/PTSTEP) + 1
END IF
!
IF(POUT12 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT12/PTSTEP) + 1
END IF
!
IF(POUT13 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT13/PTSTEP) + 1
END IF
!
IF(POUT14 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT14/PTSTEP) + 1
END IF
!
IF(POUT15 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT15/PTSTEP) + 1
END IF
!
IF(POUT16 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT16/PTSTEP) + 1
END IF
!
IF(POUT17 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT17/PTSTEP) + 1
END IF
!
IF(POUT18 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT18/PTSTEP) + 1
END IF
!
IF(POUT19 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT19/PTSTEP) + 1
END IF
!
IF(POUT20 /= -999.) THEN
  KOUT_NUMB = KOUT_NUMB + 1
  KOUT_TIMES(KOUT_NUMB) = NINT(POUT20/PTSTEP) + 1
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       3.    PRINT ON OUTPUT-LISTING 
!              -----------------------
!
CALL FMLOOK(HLUOUT,HLUOUT,ILUOUT,IRESP)
IF  (NVERB >= 10) THEN
  IIUP = SIZE(PXHAT)
  IJUP = SIZE(PYHAT) 
  WRITE(ILUOUT,*) ' SET_GRID : XLON0 = ', XLON0,' XLAT0 = ',XLAT0, &
       ' XRPK = ',XRPK,' XBETA = ',XBETA,' PLONOR = ',PLONOR,      &
       ' PLATOR = ' , PLATOR
  IF(LCARTESIAN) THEN
    WRITE(ILUOUT,*) 'SET_GRID : No map projection used.'
  ELSE
    IF (XRPK == 1.) THEN
      WRITE(ILUOUT,*) 'SET_GRID : Polar stereo used.'
    ELSE IF (XRPK == 0.) THEN
      WRITE(ILUOUT,*) 'SET_GRID : Mercator used.'
    ELSE
      WRITE(ILUOUT,*) 'SET_GRID : Lambert used, cone factor=',XRPK 
    END IF
  END IF
  WRITE(ILUOUT,*) ' SET_GRID : Some PXHAT values:'
  WRITE(ILUOUT,*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,*) PXHAT(1),PXHAT(IIUP/2),PXHAT(IIUP) 
! 
  WRITE(ILUOUT,*) ' SET_GRID : Some PYHAT values:'
  WRITE(ILUOUT,*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,*) PYHAT(1),PYHAT(IJUP/2),PYHAT(IJUP)  
!
  WRITE(ILUOUT,*) ' SET_GRID : Some PZHAT values:'
  WRITE(ILUOUT,*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,*) PZHAT(1),PZHAT(KKU/2),PZHAT(KKU) 
! 
  WRITE(ILUOUT,*) ' SET_GRID : Some PZS values:'
  WRITE(ILUOUT,*) ' I= 1        I=IIU/2       I=IIU'
  WRITE(ILUOUT,*) PZS(1,1),PZS(IIUP/2,IJUP/2),PZS(IIUP,IJUP)  
!
  YTITLE='CURRENT DATE AND TIME'
  CALL SM_PRINT_TIME(TPDTCUR,HLUOUT,YTITLE)
END IF
IF (NVERB >= 5) THEN
  YTITLE='DATE AND TIME OF EXPERIMENT BEGINNING'
  CALL SM_PRINT_TIME(TDTEXP,HLUOUT,YTITLE)
  YTITLE='DATE AND TIME OF MODEL BEGINNING'
  CALL SM_PRINT_TIME(TPDTMOD,HLUOUT,YTITLE)
END IF
YTITLE='DATE AND TIME OF SEGMENT BEGINNING'
CALL SM_PRINT_TIME(TDTSEG,HLUOUT,YTITLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_GRID
