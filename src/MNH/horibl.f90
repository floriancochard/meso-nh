!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODI_HORIBL
!     ##################
!
INTERFACE
    SUBROUTINE HORIBL(PILA1,PILO1,PILA2,PILO2,KINLA,KINLO,KILEN,PARIN, &
                      KOLEN,PXOUT,PYOUT,PAROUT,ODVECT,                 &
                      PTIME,OMINMAX,                                    &
                      KLSMIN,KLSMOUT                           )
!
REAL,                      INTENT(IN)  :: PILA1   ! Lat. (y) of first input point
REAL,                      INTENT(IN)  :: PILO1   ! Lon. (x) of first input point
REAL,                      INTENT(IN)  :: PILA2   ! Lat. (y) of last input point
REAL,                      INTENT(IN)  :: PILO2   ! Lon. (x) of last input point
INTEGER,                   INTENT(IN)  :: KINLA   ! Number of parallels
INTEGER, DIMENSION(KINLA), INTENT(IN)  :: KINLO   ! Nb. of points on each parallel
INTEGER,                   INTENT(IN)  :: KILEN   ! size of input arrays
REAL,    DIMENSION(KILEN), INTENT(IN)  :: PARIN   ! input array
INTEGER,                   INTENT(IN)  :: KOLEN   ! size of output array
REAL,    DIMENSION(KOLEN), INTENT(IN)  :: PXOUT   ! X (lon.) of output points
REAL,    DIMENSION(KOLEN), INTENT(IN)  :: PYOUT   ! Y (lat.) of output points
REAL,    DIMENSION(KOLEN), INTENT(OUT) :: PAROUT  ! output array
LOGICAL,                   INTENT(IN)  :: ODVECT  ! data is vectorial (True/False)
REAL,                      INTENT(INOUT) :: PTIME ! time spent in routine
LOGICAL,                   INTENT(IN)    :: OMINMAX ! TRUE pour borner les valeurs aux min et max du champ d entree 
INTEGER, DIMENSION(KILEN), INTENT(IN), OPTIONAL  :: KLSMIN  ! input land/sea mask
INTEGER, DIMENSION(KOLEN), INTENT(IN), OPTIONAL  :: KLSMOUT ! output land/sea mask
!
END SUBROUTINE HORIBL
!
END INTERFACE
!
END MODULE MODI_HORIBL
!
!
!   ###########################################################################
    SUBROUTINE HORIBL(PILA1,PILO1,PILA2,PILO2,KINLA,KINLO,KILEN,PARIN, &
                      KOLEN,PXOUT,PYOUT,PAROUT,ODVECT,                 &
                      PTIME,OMINMAX,                                    &
                      KLSMIN,KLSMOUT                           )
!   ###########################################################################
!
!!****  *HORIBL* - horitontal bilinear interpolation
!!
!!    PURPOSE
!!    -------
!!
!!    Interpolates a field, supports masks.
!!
!!    METHOD
!!    ------
!!
!!    This routine performs a bilinear interpolation based on the 12 surrounding
!!    points. It begins with an interpolation along the latitudes (with third order
!!    polynoms interpolation with 4 points and linear interpolation for 2 points)
!!    and then a second along the longitude (third order polynoms interpolation).
!!    Two interpolations are performed : first along the parallels then between the
!!    four resulting points.
!!
!!    The disposition of the points is the following :
!!
!!
!!            N         1   2
!!
!!            ^     3   4   5   6
!!            |           x
!!            |     7   8   9  10
!!            |
!!                     11  12
!!            S
!!              W ---------------> E
!!
!!   Note : the name 'south', 'north', may not be exact if the last data point is
!!     to the south of first (delta latitude < 0). This does not affect computations.
!!
!!   The formula used to compute the weight is :
!!        (Lon   - Lon.i) . (Lon   - Lon.i) . (Lon   - Lon.i)
!!   Wi = ---------------------------------------------------
!!        (Lon.i - Lon.j) . (Lon.i - Lon.k) . (Lon.i - Lon.l)
!!   Where j,k,l are the other points of the line.
!!
!!   When masks are used, points with different types than the output points are
!!   not taken in account (in the formula, the corresponding coefficient is set
!!   to 1). If no points of the same nature are available, the interpolation is
!!   performed anyway with the 12 points. It is the task of the calling program
!!   to react to this situation.
!!
!!   When the inputs parameters define a circular map (or global), the inputs data
!!   are extended. The value of the parameter ODVECT is used to know if the datas
!!   are vectorial or scalar (this affects the sign of extended values).
!!
!!   EXTERNAL
!!   --------
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
!!   REFERENCE
!!   ---------
!!
!!   This routine is based on the one used by the software FULL-POS from Meteo France.
!!   More informations may be found in 'Book 1'
!!
!!   AUTHOR
!!   ------
!!
!!   J.Pettre & V.Bousquet
!!
!!   MODIFICATIONS
!!   -------------
!!
!!   Original       07/01/1999
!!                  21/04/1999 (V. Masson) set correct prefixes and bug in
!!                             a logical definition
!!                  21/04/1999 (V. Masson) bug in north and south poles
!!                             extension for input map land-sea mask
!!                  27/05/1999 (V. Masson) bug in 'grib south pole'
!!                             extrapolation (number of point per parallel)
!!                  27/05/1999 (V. Masson) bug in 'grib pole' extrapolation
!!                             extra latitudes are now computed symetrically
!!                             to the poles.
!!                  16/06/2010 (G. Tanguy) bug in 'grib north pole"
!!                              extrapolation (tabular ZARIN not totaly filled)
!!                  12/10/2012 J.Escobar & F.Tocquer , interface mismatch, remove OPTIONAL for OMINMAX
!!                  2011       (P.Peyrille) 2D problem in the formulation
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!------------------------------------------------------------------------------
!
!
!*      0. DECLARATIONS
!       ---------------
!
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
!
USE MODD_LUNIT
!
USE MODD_PARAMETERS,ONLY : XUNDEF
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!*      0.1. Declaration of arguments
!         
REAL,                      INTENT(IN)  :: PILA1   ! Lat. (y) of first input point
REAL,                      INTENT(IN)  :: PILO1   ! Lon. (x) of first input point
REAL,                      INTENT(IN)  :: PILA2   ! Lat. (y) of last input point
REAL,                      INTENT(IN)  :: PILO2   ! Lon. (x) of last input point
INTEGER,                   INTENT(IN)  :: KINLA   ! Number of parallels
INTEGER, DIMENSION(KINLA), INTENT(IN)  :: KINLO   ! Number of point along a parallel
INTEGER,                   INTENT(IN)  :: KILEN   ! size of input arrays
REAL,    DIMENSION(KILEN), INTENT(IN)  :: PARIN   ! input array
INTEGER,                   INTENT(IN)  :: KOLEN   ! size of output array
REAL,    DIMENSION(KOLEN), INTENT(IN)  :: PXOUT   ! X (lon.) of output points
REAL,    DIMENSION(KOLEN), INTENT(IN)  :: PYOUT   ! Y (lat.) of output points
REAL,    DIMENSION(KOLEN), INTENT(OUT) :: PAROUT  ! output array
LOGICAL,                   INTENT(IN)  :: ODVECT  ! data is vectorial (True/False)
REAL,                      INTENT(INOUT) :: PTIME ! time spent in routine
LOGICAL,                   INTENT(IN)  :: OMINMAX ! TRUE pour borner les valeurs aux min et max du champ d entree 
INTEGER, DIMENSION(KILEN), INTENT(IN), OPTIONAL  :: KLSMIN  ! input land/sea mask
INTEGER, DIMENSION(KOLEN), INTENT(IN), OPTIONAL  :: KLSMOUT ! output land/sea mask

!
!*      0.2. Declaration of local variables
!            
 ! Variables used to perform the interpolation
REAL                               :: ZOLA     ! Latitude of the output point
REAL                               :: ZOLO     ! Longitude of the output point
REAL                               :: ZIDLA    ! Delta latitude
REAL                               :: ZIDLO    ! Delta longitude
INTEGER, DIMENSION(:), ALLOCATABLE :: IOFS     ! Offset of each parallel in the array
  ! Number of the surrounding latitudes
INTEGER                            :: IOSS,IOS,ION,IONN
  ! Posiiton in the array of the twelwe surrounding points
INTEGER                            :: IP1,IP2,IP3,IP4,IP5,IP6,IP7,IP8,IP9,IP10, &
                                      IP11,IP12
  ! Latitudes and longitudes of the surrounding points
REAL                               :: ZLANN,ZLAN,ZLAS,ZLASS
REAL                               :: ZLOP1,ZLOP2,ZLOP3,ZLOP4 ,ZLOP5 ,ZLOP6,    &
                                      ZLOP7,ZLOP8,ZLOP9,ZLOP10,ZLOP11,ZLOP12
  ! Weights of the latitudes and of the points
REAL                               :: ZWNN,ZWN,ZWS,ZWSS
REAL                               :: ZW1,ZW2,ZW3,ZW4,ZW5,ZW6,ZW7,ZW8,ZW9,ZW10, &
                                      ZW11,ZW12
  ! Land/sea mask coefficient for each point : 0 -> point not taken in account,
  !                                            1 -> point taken in account
REAL                               :: ZLSM1,ZLSM2 ,ZLSM3 ,ZLSM4 ,ZLSM5 ,ZLSM6,ZLSM7,ZLSM8, &
                                      ZLSM9,ZLSM10,ZLSM11,ZLSM12,ZLSMNN,ZLSMN,ZLSMS,ZLSMSS,&
                                      ZLSMTOT
 ! Variables implied in the extension procedure
REAL                               :: ZILO1     ! Longitude of the first data point
REAL                               :: ZILO2     ! Longitude of the last data point
LOGICAL                            :: GGLOBLON  ! True if the map is circular
LOGICAL                            :: GGLOBN    ! True if the map has the north pole
LOGICAL                            :: GGLOBS    ! True if the map has the south pole
INTEGER                            :: IBIGSIZE  ! Size of the extended map
INTEGER                            :: IMIDDLE   ! Used for extensions around the poles
INTEGER                            :: IOFFSET1  ! Offset in map
INTEGER                            :: IOFFSET2  ! Offset in map
REAL                               :: ZSOUTHPOLE! south pole latitude (-90 or  90)
REAL                               :: ZNORTHPOLE! north pole latitude ( 90 or -90)
REAL,    DIMENSION(:), ALLOCATABLE :: ZARIN     ! Extended input datas
INTEGER, DIMENSION(:), ALLOCATABLE :: ILSMIN    ! Extended land/sea mask
INTEGER, DIMENSION(:), ALLOCATABLE :: IINLO     ! Extended KINLO
INTEGER                            :: IINLA     ! Number of parallel
REAL                               :: ZVECT     ! -1 if input is vectorial
LOGICAL                            :: LDLSM     ! Specify if land/sea mask is present or not
 ! Variables used for the out put listing
INTEGER                            :: ILUOUT0   ! Logical unit number
 ! Loop counters
INTEGER                            :: JOPOS     ! Output position
INTEGER                            :: JIPOS     ! Input position
INTEGER                            :: JLOOP1    ! Dummy counter
!
REAL                               :: ZTIME1    ! CPU time start counter
REAL                               :: ZTIME2    ! CPU time end   counter
!
!------------------------------------------------------------------------------
REAL                               :: ZMAX      ! Max of 12 surrounding values
REAL                               :: ZMIN      ! Min of 12 surrounding values
INTEGER                            :: JLOOP2    ! Dummy counter
INTEGER,    DIMENSION(12)          :: IP        ! Array for IPn 
INTEGER                            :: IRESP     ! Return code of FM-routines
!------------------------------------------------------------------------------
!
CALL SECOND_MNH(ZTIME1)
!
ILUOUT0 = TLUOUT0%NLU
!
!------------------------------------------------------------------------------
!
!*     1. DETERMINATION  of the latitude of the poles (depending of the latitude
!         -------------                                 of the first data point)
!
IF (PILA1>0.) THEN
  ZSOUTHPOLE= 90.
  ZNORTHPOLE=-90.
ELSE
  ZSOUTHPOLE=-90.
  ZNORTHPOLE= 90.
END IF
!
!------------------------------------------------------------------------------
!
!*     2. EXTEND DATA GRID
!         ----------------
  ! Land / Sea mask
LDLSM = .FALSE.
IF (PRESENT(KLSMIN) .AND. PRESENT(KLSMOUT)) LDLSM = .TRUE.
!
!*    2.1 Alias input data
!
ZILO1 = PILO1
ZILO2 = PILO2
ZVECT = 1.
IF (ODVECT) ZVECT=-1.
!
!*   2.2 Center input domain in order to have Lo1 < Lo 2
!
IF (ZILO2 < 0.)    ZILO2 = ZILO2 + 360.
IF (ZILO1 < 0.)    ZILO1 = ZILO1 + 360.
IF (ZILO2 < ZILO1) ZILO1 = ZILO1 - 360.
!
!*   2.3 Extend one point (needed for reduced grids)
!
! Longitude coordinate of points are found by :
!                      i
!  Lon(i) = Lon1 + ------------- . (Lon2 - Lon1)
!                   Npts(Lat)-1
! Where i goes from 0 to Npts(Lat)-1. The result of this is that the last point of 
! each parallel is located at Lon2. This is not the case for reduced grid where the 
! position of the last point depends upon the number of points of the parallel. For
! reduced grid, the right formula to use is the following :
!                       i
!  Lon(i) = Lon1 + ----------- . (Lon2' - Lon1)
!                   Npts(Lat)
! Where Lon2' = Lon1 + 2.PI.
!
!                                              Lon2 - Lon1
! This can be generalized with Lon2' = Lon2 + -------------
!                                              Nptsmax - 1
!
JOPOS = MAXVAL(KINLO(1:KINLA))
ZILO2 = ZILO1 + (ZILO2 - ZILO1) * JOPOS / (JOPOS - 1.)
!
!* 2.4 Test if the input is global or partially global
!
! Note that we must have a global map to make extension around the poles
GGLOBN   = .FALSE.
GGLOBS   = .FALSE.
GGLOBLON = .FALSE.
IF (ZILO2-360.>ZILO1-1.E-3) GGLOBLON = .TRUE.
ZIDLA = (PILA2 - PILA1) / (KINLA - 1)
IF ((PILA1-ZIDLA>= 90.) .OR. (PILA1-ZIDLA<=-90.)) GGLOBS=GGLOBLON
IF ((PILA2+ZIDLA>= 90.) .OR. (PILA2+ZIDLA<=-90.)) GGLOBN=GGLOBLON
! Aladin case (input PILA2, PILO2 are in meters) no extension
IF ( PILA2 > 100. ) THEN
  GGLOBN   = .FALSE.
  GGLOBS   = .FALSE.
  GGLOBLON = .FALSE.
END IF
!
!* 2.5  Compute the size of the resulting map
!
IBIGSIZE = KILEN
IF (GGLOBS  ) IBIGSIZE=IBIGSIZE+(4+KINLO(    1))+(4+KINLO(      2))
IF (GGLOBN  ) IBIGSIZE=IBIGSIZE+(4+KINLO(KINLA))+(4+KINLO(KINLA-1))
IF (GGLOBLON) IBIGSIZE=IBIGSIZE+ 4*KINLA
!
!* 2.6 Compute the resulting map
!
ALLOCATE (ZARIN(IBIGSIZE))
ALLOCATE (ILSMIN(IBIGSIZE))
!
! 2.6.1 Compute the longitude extension
!
! This is a basic copy of the data. If extension is possible, the first and last
! two lines are copied twice this way :
!
!    /---------------\
!    |               |
!   [.] [.] [....   ...] [.] [.] 
!        |            |
!        \------------/
!
! A point represent a data.
!
JIPOS = 1
JOPOS = 1
IF (GGLOBS) JOPOS=JOPOS+(4+KINLO(1))+(4+KINLO(2))
IF (GGLOBLON) THEN
  DO JLOOP1 = 1, KINLA
    ZARIN(JOPOS  ) = PARIN(JIPOS+KINLO(JLOOP1)-2)
    ZARIN(JOPOS+1) = PARIN(JIPOS+KINLO(JLOOP1)-1)
    ZARIN(JOPOS+2:JOPOS+2+KINLO(JLOOP1)-1) = PARIN(JIPOS:JIPOS+KINLO(JLOOP1)-1)
    ZARIN(JOPOS+2+KINLO(JLOOP1)  ) = PARIN(JIPOS  )
    ZARIN(JOPOS+2+KINLO(JLOOP1)+1) = PARIN(JIPOS+1)
    IF (LDLSM) THEN
      ILSMIN(JOPOS  ) = KLSMIN(JIPOS+KINLO(JLOOP1)-2)
      ILSMIN(JOPOS+1) = KLSMIN(JIPOS+KINLO(JLOOP1)-1)
      ILSMIN(JOPOS+2:JOPOS+2+KINLO(JLOOP1)-1) = KLSMIN(JIPOS:JIPOS+KINLO(JLOOP1)-1)
      ILSMIN(JOPOS+2+KINLO(JLOOP1)  ) = KLSMIN(JIPOS  )
      ILSMIN(JOPOS+2+KINLO(JLOOP1)+1) = KLSMIN(JIPOS+1)
    END IF
    JIPOS = JIPOS + KINLO(JLOOP1)
    JOPOS = JOPOS + KINLO(JLOOP1) + 4
  END DO
ELSE
  ZARIN(JOPOS:JOPOS+KILEN-1) = PARIN(JIPOS:JIPOS+KILEN-1)
  IF (LDLSM) THEN
    ILSMIN(JOPOS:JOPOS+KILEN-1) = KLSMIN(JIPOS:JIPOS+KILEN-1)
  END IF
END IF
!
! 2.6.2 Compute the south pole extension
!
! Pole extension is performed by copying the first half datas to the last half 
! datas of the extension parallel :
!
!  [.] [.] [....] [....] [.] [.]
!                  ||||
!            /-------/
!           ||||
!  [.] [.] [....] [....] [.] [.]
!
IF (GGLOBS) THEN ! South pole (south meaning begining of the grib)
  IOFFSET1 = 4 + KINLO(2)
  IOFFSET2 = IOFFSET1 + 4 + KINLO(1)
  IMIDDLE = (KINLO(1)+4) / 2
  ZARIN(IOFFSET1+1:IOFFSET1+IMIDDLE) = &
    ZVECT*ZARIN(IOFFSET2+1+IMIDDLE-2:IOFFSET2+2*IMIDDLE-2)
  ZARIN(IOFFSET1+IMIDDLE+1:IOFFSET1+KINLO(1)+4) = &
    ZVECT*ZARIN(IOFFSET2+1+2:IOFFSET2+KINLO(1)+4-IMIDDLE+2)
  IF (LDLSM) THEN
    ILSMIN(IOFFSET1+1:IOFFSET1+IMIDDLE) = &
      ILSMIN(IOFFSET2+1+IMIDDLE-2:IOFFSET2+2*IMIDDLE-2)
    ILSMIN(IOFFSET1+IMIDDLE+1:IOFFSET1+KINLO(1)+4) = &
      ILSMIN(IOFFSET2+1+2:IOFFSET2+KINLO(1)+4-IMIDDLE+2)
  END IF
  IOFFSET2 = IOFFSET2 + 4 + KINLO(1)
  IMIDDLE = (KINLO(2)+4) / 2
  ZARIN(1:IMIDDLE) = ZVECT*ZARIN(IOFFSET2+1+IMIDDLE-2:IOFFSET2+2*IMIDDLE-2)
  ZARIN(IMIDDLE+1:KINLO(2)+4) = &
    ZVECT*ZARIN(IOFFSET2+1+2:IOFFSET2+KINLO(2)+4-IMIDDLE+2)
  IF (LDLSM) THEN
    ILSMIN(1:IMIDDLE) = ILSMIN(IOFFSET2+1+IMIDDLE:IOFFSET2+2*IMIDDLE)
    ILSMIN(IMIDDLE+1:KINLO(2)+4) = ILSMIN(IOFFSET2+1+2:IOFFSET2+KINLO(2)+4-IMIDDLE+2)
  END IF
END IF
!
! 2.6.3 Compute the north pole extension
!
IF (GGLOBN) THEN ! North pole (north meaning end of the grib)
  IOFFSET1 = IBIGSIZE - (4+KINLO(KINLA-1)) - (4+KINLO(KINLA))
  IOFFSET2 = IOFFSET1 - (4+KINLO(KINLA))
  IMIDDLE = (KINLO(KINLA)+4) / 2
  ZARIN(IOFFSET1+1:IOFFSET1+IMIDDLE) = &
    ZVECT*ZARIN(IOFFSET2+1+IMIDDLE-2:IOFFSET2+2*IMIDDLE-2)
  ZARIN(IOFFSET1+IMIDDLE+1:IOFFSET1+KINLO(KINLA)+4) = &
    ZVECT*ZARIN(IOFFSET2+1+2:IOFFSET2+KINLO(KINLA)+4-IMIDDLE+2)
  IF (LDLSM) THEN
    ILSMIN(IOFFSET1+1:IOFFSET1+IMIDDLE) = &
      ILSMIN(IOFFSET2+1+IMIDDLE-2:IOFFSET2+2*IMIDDLE-2)
    ILSMIN(IOFFSET1+IMIDDLE+1:IOFFSET1+KINLO(KINLA)+4) = &
      ILSMIN(IOFFSET2+1+2:IOFFSET2+KINLO(KINLA)+4-IMIDDLE+2)
  END IF
  IOFFSET1 = IOFFSET1 + (4+KINLO(KINLA))
  IOFFSET2 = IOFFSET2 - (4+KINLO(KINLA-1))
  IMIDDLE = (KINLO(KINLA-1)+4) / 2
  ZARIN(IOFFSET1+1:IOFFSET1+IMIDDLE) = &
    ZVECT*ZARIN(IOFFSET2+1+IMIDDLE-2:IOFFSET2+2*IMIDDLE-2)
  ZARIN(IOFFSET1+IMIDDLE+1:IOFFSET1+KINLO(KINLA-1)+4) = &
    ZVECT*ZARIN(IOFFSET2+1+2:IOFFSET2+KINLO(KINLA-1)+4-IMIDDLE+2)
  IF (LDLSM) THEN
    ILSMIN(IOFFSET1+1:IOFFSET1+IMIDDLE) = &
      ILSMIN(IOFFSET2+1+IMIDDLE-2:IOFFSET2+2*IMIDDLE-2)
    ILSMIN(IOFFSET1+IMIDDLE+1:IOFFSET1+KINLO(KINLA-1)+4) = &
      ILSMIN(IOFFSET2+1+2:IOFFSET2+KINLO(KINLA-1)+4-IMIDDLE+2)
  END IF
END IF
!
!*  2.7  Compute the resulting parameters of the map
!
IINLA = KINLA
IF (GGLOBS) IINLA = IINLA + 2
IF (GGLOBN) IINLA = IINLA + 2
!
ALLOCATE (IINLO(IINLA))
IOFFSET1 = 0
IF (GGLOBS) THEN
  IINLO(IOFFSET1+1) = KINLO(2)
  IINLO(IOFFSET1+2) = KINLO(1)
  IOFFSET1 = IOFFSET1 + 2
END IF
IINLO(IOFFSET1+1:IOFFSET1+KINLA) = KINLO(1:KINLA)
IOFFSET1 = IOFFSET1 + KINLA
IF (GGLOBN) THEN
  IINLO(IOFFSET1+1) = KINLO(KINLA)
  IINLO(IOFFSET1+2) = KINLO(KINLA-1)
  IOFFSET1 = IOFFSET1 + 2
END IF
!
!*  2.8  Compute Offset array used to acces the lines
!
ALLOCATE (IOFS(IINLA))
IOFS(1) = 1
IF (GGLOBLON) IOFS(1)=IOFS(1)+2
DO JLOOP1=2, IINLA
  IOFS(JLOOP1) = IOFS(JLOOP1-1) + IINLO(JLOOP1-1)
  IF (GGLOBLON) IOFS(JLOOP1) = IOFS(JLOOP1) + 4
END DO
!
!------------------------------------------------------------------------------
!
!*     3.   LOOP OVER ALL THE POINTS OF THE OUTPUT GRID
!           -------------------------------------------
JOPOS = 0
DO JLOOP1 = 1, KOLEN
  JOPOS = JOPOS + 1
  ZOLA  = PYOUT(JOPOS)
  ZOLO  = PXOUT(JOPOS)
  ! PPeyrille modification
! IF (ZOLO < ZILO1) ZOLO = ZOLO + 360.
! IF (ZOLO > ZILO2) ZOLO = ZOLO - 360.
  IF ((ZOLO-ZILO1) < -0.1E-06) ZOLO = ZOLO + 360.
  IF ((ZOLO -ZILO2) > 0.1E-06) ZOLO = ZOLO - 360.  
!
!* 3.1 Locate the 12 input points around the interpolated output point
!*
!*            N         1   2
!*
!*            ^     3   4   5   6
!*            |           x
!*            |     7   8   9  10
!*            |
!*                     11  12
!*            S
!*              W ---------------> E
!*
!* Note : the name 'south', 'north', may not be exact if the point 2 is
!*   to the south of point 1 (IDLA < 0). This does not affect computation.
!
    ! 3.1.1. find positions of latitudes
  IOS = NINT( (ZOLA-PILA1)/ZIDLA - 0.5) ! because of the zero
  ZLAS = PILA1 + IOS * ZIDLA
  IOS  = IOS + 1
  IF (GGLOBS)  IOS  = IOS + 2
  IOSS = IOS - 1
  ION  = IOS + 1
  IONN = ION + 1
  ZLASS = ZLAS - ZIDLA
  ZLAN  = ZLAS + ZIDLA
  ZLANN = ZLAN + ZIDLA
      !
      ! extra latitudes are computed symetrically compared to the poles
      !
  IF (GGLOBS .AND. IOS==2) THEN
    ZLASS = 2. * ZSOUTHPOLE - ZLANN
    ZLAS  = 2. * ZSOUTHPOLE - ZLAN
  END IF
  IF (GGLOBS .AND. IOS==3) THEN
    ZLASS = 2. * ZSOUTHPOLE - ZLAS
  END IF
  IF (GGLOBN .AND. IOS==IINLA-2) THEN
    ZLANN = 2. * ZNORTHPOLE - ZLASS
    ZLAN  = 2. * ZNORTHPOLE - ZLAS
  END IF
  IF (GGLOBN .AND. IOS==IINLA-3) THEN
     ZLANN = 2. * ZNORTHPOLE - ZLAN
  END IF
!
  IF ((IOSS<1).OR.(IONN>IINLA).OR. &
     (IOSS<1).OR.(IONN>IINLA)) THEN
!callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','HORIBL','input domain is smaller than output one')
 END IF
!
      ! 3.1.2. northern
  ZIDLO = (ZILO2 - ZILO1) / (IINLO(IONN))
  IP1   = INT((ZOLO - ZILO1) / ZIDLO)
  IP2   = IP1  + 1
  ZLOP1 = ZILO1 + IP1 * ZIDLO
  ZLOP2 = ZLOP1 + ZIDLO
!
      ! 3.1.3. north
  ZIDLO = (ZILO2 - ZILO1) / (IINLO(ION ))
  IP4   = INT((ZOLO - ZILO1) / ZIDLO)
  IP3   = IP4  - 1
  IP5   = IP4  + 1
  IP6   = IP5  + 1
  ZLOP4 = ZILO1 + IP4 * ZIDLO
  ZLOP3 = ZLOP4 - ZIDLO
  ZLOP5 = ZLOP4 + ZIDLO
  ZLOP6 = ZLOP5 + ZIDLO
!
      ! 3.1.4. south
  ZIDLO = (ZILO2 - ZILO1) / (IINLO(IOS ))
  IP8   = INT((ZOLO - ZILO1) / ZIDLO)
  IP7   = IP8  - 1
  IP9   = IP8  + 1
  IP10  = IP9  + 1
  ZLOP8 = ZILO1 + IP8 * ZIDLO
  ZLOP7 = ZLOP8 - ZIDLO
  ZLOP9 = ZLOP8 + ZIDLO
  ZLOP10= ZLOP9 + ZIDLO
!
      ! 3.1.5. southern
  ZIDLO = (ZILO2 - ZILO1) / (IINLO(IOSS))
  IP11  = INT((ZOLO - ZILO1) / ZIDLO)
  IP12  = IP11 + 1
  ZLOP11= ZILO1 + IP11* ZIDLO
  ZLOP12= ZLOP11+ ZIDLO
!
      ! 3.1.6. check position of points
  IF (GGLOBLON) THEN
    IF ((IP1 <-2) .OR. (IP2 >IINLO(IONN)+1) .OR. &
        (IP3 <-2) .OR. (IP6 >IINLO(ION )+1) .OR. &
        (IP7 <-2) .OR. (IP10>IINLO(IOS )+1) .OR. &
      (IP11<-2) .OR. (IP12>IINLO(IOSS)+1)) THEN
!callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','HORIBL','input domain is smaller than output one - longitude global')
    END IF
  ELSE
    IF ((IP1 <0) .OR. (IP2 >IINLO(IONN)-1) .OR. &
        (IP3 <0) .OR. (IP6 >IINLO(ION )-1) .OR. &
        (IP7 <0) .OR. (IP10>IINLO(IOS )-1) .OR. &
        (IP11<0) .OR. (IP12>IINLO(IOSS)-1)) THEN
!callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','HORIBL','input domain is smaller than output one - longitude local')
    END IF
  END IF
!
      ! 3.1.7. add parallel offset
  IP1 =IP1 + IOFS(IONN)
  IP2 =IP2 + IOFS(IONN)
  IP3 =IP3 + IOFS(ION )
  IP4 =IP4 + IOFS(ION )
  IP5 =IP5 + IOFS(ION )
  IP6 =IP6 + IOFS(ION )
  IP7 =IP7 + IOFS(IOS )
  IP8 =IP8 + IOFS(IOS )
  IP9 =IP9 + IOFS(IOS )
  IP10=IP10+ IOFS(IOS )
  IP11=IP11+ IOFS(IOSS)
  IP12=IP12+ IOFS(IOSS)
!
!*  3.2 Land / Sea mask
!
  ZLSM1  = 1.
  ZLSM2  = 1.
  ZLSM3  = 1.
  ZLSM4  = 1.
  ZLSM5  = 1.
  ZLSM6  = 1.
  ZLSM7  = 1.
  ZLSM8  = 1.
  ZLSM9  = 1.
  ZLSM10 = 1.
  ZLSM11 = 1.
  ZLSM12 = 1.
  ZLSMNN = 1.
  ZLSMN  = 1.
  ZLSMS  = 1.
  ZLSMSS = 1.
  IF (LDLSM) THEN
    IF (ILSMIN(IP1 ).NE.KLSMOUT(JOPOS)) ZLSM1  = 0.
    IF (ILSMIN(IP2 ).NE.KLSMOUT(JOPOS)) ZLSM2  = 0.
    IF (ILSMIN(IP3 ).NE.KLSMOUT(JOPOS)) ZLSM3  = 0.
    IF (ILSMIN(IP4 ).NE.KLSMOUT(JOPOS)) ZLSM4  = 0.
    IF (ILSMIN(IP5 ).NE.KLSMOUT(JOPOS)) ZLSM5  = 0.
    IF (ILSMIN(IP6 ).NE.KLSMOUT(JOPOS)) ZLSM6  = 0.
    IF (ILSMIN(IP7 ).NE.KLSMOUT(JOPOS)) ZLSM7  = 0.
    IF (ILSMIN(IP8 ).NE.KLSMOUT(JOPOS)) ZLSM8  = 0.
    IF (ILSMIN(IP9 ).NE.KLSMOUT(JOPOS)) ZLSM9  = 0.
    IF (ILSMIN(IP10).NE.KLSMOUT(JOPOS)) ZLSM10 = 0.
    IF (ILSMIN(IP11).NE.KLSMOUT(JOPOS)) ZLSM11 = 0.
    IF (ILSMIN(IP12).NE.KLSMOUT(JOPOS)) ZLSM12 = 0.
    ZLSMNN = MIN(ZLSM1 +ZLSM2,1.)
    ZLSMN  = MIN(ZLSM3 +ZLSM4 +ZLSM5 +ZLSM6,1.)
    ZLSMS  = MIN(ZLSM7 +ZLSM8 +ZLSM9 +ZLSM10,1.)
    ZLSMSS = MIN(ZLSM11+ZLSM12,1.)
    ZLSMTOT = MIN(ZLSMNN+ZLSMN+ZLSMS+ZLSMSS,1.)
    IF (ZLSMNN < 1.E-3) THEN
      ZLSM1 = 1.
      ZLSM2 = 1.
    END IF
    IF (ZLSMN  < 1.E-3) THEN
      ZLSM3 = 1.
      ZLSM4 = 1.
      ZLSM5 = 1.
      ZLSM6 = 1.
    END IF
    IF (ZLSMS  < 1.E-3) THEN
      ZLSM7 = 1.
      ZLSM8 = 1.
      ZLSM9 = 1.
      ZLSM10= 1.
    END IF
    IF (ZLSMSS < 1.E-3) THEN
      ZLSM11= 1.
      ZLSM12= 1.
    END IF
    IF (ZLSMTOT < 1.E-3) THEN
      ZLSMNN = 1.
      ZLSMN  = 1.
      ZLSMS  = 1.
      ZLSMSS = 1.
    END IF
  ENDIF
!
!*  3.3 Weight of points
!
      ! 3.3.1 northern
  ZW1  = ZLSM1 * (1.+ZLSM2 *(ZOLO -ZLOP1 )/(ZLOP1 -ZLOP2 ))
  ZW2  = 1. - ZW1
  ZWNN = ZLSMNN* (1.+ZLSMN *(ZOLA -ZLANN)/(ZLANN-ZLAN )) &
               * (1.+ZLSMS *(ZOLA -ZLANN)/(ZLANN-ZLAS )) &
               * (1.+ZLSMSS*(ZOLA -ZLANN)/(ZLANN-ZLASS))
!
      ! 3.3.2. north
  ZW3  = ZLSM3 * (1.+ZLSM4 *(ZOLO -ZLOP3 )/(ZLOP3 -ZLOP4 )) &
               * (1.+ZLSM5 *(ZOLO -ZLOP3 )/(ZLOP3 -ZLOP5 )) &
               * (1.+ZLSM6 *(ZOLO -ZLOP3 )/(ZLOP3 -ZLOP6 )) 
  ZW4  = ZLSM4 * (1.+ZLSM3 *(ZOLO -ZLOP4 )/(ZLOP4 -ZLOP3 )) &
               * (1.+ZLSM5 *(ZOLO -ZLOP4 )/(ZLOP4 -ZLOP5 )) &
               * (1.+ZLSM6 *(ZOLO -ZLOP4 )/(ZLOP4 -ZLOP6 ))
  ZW5  = ZLSM5 * (1.+ZLSM3 *(ZOLO -ZLOP5 )/(ZLOP5 -ZLOP3 )) &
               * (1.+ZLSM4 *(ZOLO -ZLOP5 )/(ZLOP5 -ZLOP4 )) &
               * (1.+ZLSM6 *(ZOLO -ZLOP5 )/(ZLOP5 -ZLOP6 ))
  ZW6 = 1. - ZW3 - ZW4 - ZW5
  ZWN  = ZLSMN * (1.+ZLSMNN*(ZOLA -ZLAN )/(ZLAN -ZLANN)) &
               * (1.+ZLSMS *(ZOLA -ZLAN )/(ZLAN -ZLAS )) &
               * (1.+ZLSMSS*(ZOLA -ZLAN )/(ZLAN -ZLASS))
!
      ! 3.3.3. south
  ZW7  = ZLSM7 * (1.+ZLSM8 *(ZOLO -ZLOP7 )/(ZLOP7 -ZLOP8 )) &
               * (1.+ZLSM9 *(ZOLO -ZLOP7 )/(ZLOP7 -ZLOP9 )) &
               * (1.+ZLSM10*(ZOLO -ZLOP7 )/(ZLOP7 -ZLOP10))
  ZW8  = ZLSM8 * (1.+ZLSM7 *(ZOLO -ZLOP8 )/(ZLOP8 -ZLOP7 )) &
               * (1.+ZLSM9 *(ZOLO -ZLOP8 )/(ZLOP8 -ZLOP9 )) &
               * (1.+ZLSM10*(ZOLO -ZLOP8 )/(ZLOP8 -ZLOP10))
  ZW9  = ZLSM9 * (1.+ZLSM7 *(ZOLO -ZLOP9 )/(ZLOP9 -ZLOP7 )) &
               * (1.+ZLSM8 *(ZOLO -ZLOP9 )/(ZLOP9 -ZLOP8 )) &
               * (1.+ZLSM10*(ZOLO -ZLOP9 )/(ZLOP9 -ZLOP10))
  ZW10 = 1. - ZW7 - ZW8 - ZW9
  ZWS  = ZLSMS * (1.+ZLSMNN*(ZOLA -ZLAS )/(ZLAS -ZLANN)) &
               * (1.+ZLSMN *(ZOLA -ZLAS )/(ZLAS -ZLAN )) &
               * (1.+ZLSMSS*(ZOLA -ZLAS )/(ZLAS -ZLASS))
!
      ! 3.3.4. southern
  ZW11 = ZLSM11* (1.+ZLSM12*(ZOLO -ZLOP11)/(ZLOP11-ZLOP12))
  ZW12 = 1. - ZW11
  ZWSS = 1. - ZWNN - ZWN - ZWS
!
      ! 3.3.5. longitude weight x latitude weight
  ZW1  = ZW1  * ZWNN
  ZW2  = ZW2  * ZWNN
  ZW3  = ZW3  * ZWN
  ZW4  = ZW4  * ZWN
  ZW5  = ZW5  * ZWN
  ZW6  = ZW6  * ZWN
  ZW7  = ZW7  * ZWS
  ZW8  = ZW8  * ZWS
  ZW9  = ZW9  * ZWS
  ZW10 = ZW10 * ZWS
  ZW11 = ZW11 * ZWSS
  ZW12 = ZW12 * ZWSS
!
  PAROUT (JOPOS) = ZW1  * ZARIN(IP1 ) + &
                   ZW2  * ZARIN(IP2 ) + &
                   ZW3  * ZARIN(IP3 ) + &
                   ZW4  * ZARIN(IP4 ) + &
                   ZW5  * ZARIN(IP5 ) + &
                   ZW6  * ZARIN(IP6 ) + &
                   ZW7  * ZARIN(IP7 ) + &
                   ZW8  * ZARIN(IP8 ) + &
                   ZW9  * ZARIN(IP9 ) + &
                   ZW10 * ZARIN(IP10) + &
                   ZW11 * ZARIN(IP11) + &
                   ZW12 * ZARIN(IP12)        
!
! For surface fields, the interpoalted value is bounded 
! by the min max values of the initial field
!
IF (PRESENT(KLSMIN) .OR. OMINMAX ) THEN
!
  IP(1)=IP1
  IP(2)=IP2
  IP(3)=IP3
  IP(4)=IP4
  IP(5)=IP5
  IP(6)=IP6
  IP(7)=IP7
  IP(8)=IP8
  IP(9)=IP9
  IP(10)=IP10
  IP(11)=IP11
  IP(12)=IP12
  ZMIN=XUNDEF
  ZMAX=XUNDEF
  DO JLOOP2=1,12
    IF (ZARIN(IP(JLOOP2))==XUNDEF) CYCLE
    IF ((ZMAX==XUNDEF)) THEN
      ZMAX=ZARIN(IP(JLOOP2))
      ZMIN=ZARIN(IP(JLOOP2))
    ELSE
      ZMAX=MAX(ZMAX,ZARIN(IP(JLOOP2)))
      ZMIN=MIN(ZMIN,ZARIN(IP(JLOOP2)))
    ENDIF
  END DO
  PAROUT(JOPOS) = MAX(MIN(PAROUT(JOPOS),ZMAX),ZMIN)
ENDIF

END DO
!
DEALLOCATE (IINLO)
DEALLOCATE (ZARIN)
DEALLOCATE (ILSMIN)
DEALLOCATE (IOFS)
!
!
CALL SECOND_MNH(ZTIME2)
PTIME = PTIME + ZTIME2 - ZTIME1
!
RETURN
!
END SUBROUTINE HORIBL
