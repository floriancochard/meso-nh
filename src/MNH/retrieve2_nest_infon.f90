!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################################
      MODULE MODI_RETRIEVE2_NEST_INFO_n
!     ################################
!
INTERFACE 
!
      SUBROUTINE RETRIEVE2_NEST_INFO_n(KMI,KDAD,KXOR_C_ll,KYOR_C_ll,KXSIZE,KYSIZE,KDXRATIO,KDYRATIO)
!
INTEGER,INTENT(IN)  :: KMI      ! son model index
INTEGER,INTENT(IN)  :: KDAD     ! dad model index
INTEGER,INTENT(OUT) :: KXOR_C_ll     ! position of pgd model origine points
INTEGER,INTENT(OUT) :: KYOR_C_ll     ! according to father domain
INTEGER,INTENT(OUT) :: KXSIZE   ! number of grid meshes in father grid to be
INTEGER,INTENT(OUT) :: KYSIZE   ! covered by the pgd domain
INTEGER,INTENT(OUT) :: KDXRATIO ! resolution ratio between father grid
INTEGER,INTENT(OUT) :: KDYRATIO ! and son grid
!
END SUBROUTINE RETRIEVE2_NEST_INFO_n
!
END INTERFACE
!
END MODULE MODI_RETRIEVE2_NEST_INFO_n
!
!
!
!     ###############################################################
      SUBROUTINE RETRIEVE2_NEST_INFO_n(KMI,KDAD,KXOR_C_ll,KYOR_C_ll,KXSIZE,KYSIZE, &
                                           KDXRATIO,KDYRATIO)
!     ###############################################################
!
!!****  *RETRIEVE2_NEST_INFO_n* - routine to test coherence between grid
!!                                of model KMI and of current model IDAD.
!!                                
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_GRID : contains projection definition
!!        XLAT0
!!        XLON0
!!        XRPK
!!        XBETA
!!      Module MODD_PGDGRID : contains domain definition
!!        XPGDLATOR
!!        XPGDLONOR
!!        XPGDXHAT
!!        XPGDYHAT
!!      Module MODD_PGDDIM : contains domain size
!!        NPGDIMAX
!!        NPGDJMAX
!!        XLATORI
!!        XLONORI
!!      Module MODD_GRID_n :
!!        XXHAT
!!        XYHAT 
!!      Module MODD_DIM_n :
!!        NIMAX, NJMAX
!!      Module MODD_PARAMETERS :
!!        JPHEXT
!!      Module MODD_LUNIT :
!!        CLUOUT
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        25/09/96
!!                      22/09/99 PGD modules for dad, and _n module for son
!!      J Stein         04/07/01 add cartesian case
!!      M.Faivre            2014
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 01/06/2016 : Bug in type of ZBUF INTEGER => REAL & use MPI_PRECISION for r4/R8 compatibility
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CONF
USE MODD_DIM_ll,       ONLY: NXOR_ALL, NXEND_ALL, NYOR_ALL, NYEND_ALL, NIMAX_TMP_ll, NJMAX_TMP_ll
USE MODD_DIM_n,        ONLY: NIMAX, NJMAX
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IO_ll,        ONLY: ISNPROC, ISP
USE MODD_LUNIT,        ONLY: TLUOUT0
USE MODD_MPIF
USE MODD_PARAMETERS
USE MODD_PGDDIM
USE MODD_PGDGRID
USE MODD_STRUCTURE_ll, ONLY: ZONE_ll
USE MODD_VAR_ll,       ONLY: YSPLITTING, NMNH_COMM_WORLD, MPI_PRECISION
!
USE MODE_GRIDPROJ
USE MODE_MODELN_HANDLER 
USE MODE_MPPDB
USE MODE_NEST_ll,      ONLY: GO_TOMODEL_ll
USE MODE_SPLITTING_ll, ONLY: SPLIT2
USE MODE_TOOLS_ll,     ONLY: LEAST_ll, LWEST_ll, LNORTH_ll, LSOUTH_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,INTENT(IN)  :: KMI      ! son model index
INTEGER,INTENT(IN)  :: KDAD     ! dad model index
INTEGER,INTENT(OUT) :: KXOR_C_ll     ! position of pgd model origine points
INTEGER,INTENT(OUT) :: KYOR_C_ll     ! according to father (next refered as 1) domain
INTEGER,INTENT(OUT) :: KXSIZE   ! number of grid meshes in model 1 to be
INTEGER,INTENT(OUT) :: KYSIZE   ! covered by the pgd domain
INTEGER,INTENT(OUT) :: KDXRATIO ! resolution ratio between grid 1
INTEGER,INTENT(OUT) :: KDYRATIO ! and its son (next refered as 2) grid
!
!
!*       0.2   declarations of local variables
!
INTEGER              :: ILUOUT
INTEGER              :: IIU           ! relatively to model 1
INTEGER              :: IJU           ! relatively to model 1
INTEGER              :: IIUGLB           ! relatively to model 1
INTEGER              :: IJUGLB           ! relatively to model 1
INTEGER              :: IPGDIU        ! relatively to model 2
INTEGER              :: IPGDJU        ! relatively to model 2
REAL                 :: ZLAT2         ! geographical coordinates of the first
REAL                 :: ZLON2         ! local physical flux point of model 2
REAL                 :: ZLAT2GLB      ! geographical coordinates of the first
REAL                 :: ZLON2GLB      ! global physical flux point of model 2
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPGDLAT1 ! geographical coordinates of all
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPGDLON1 ! the local flux points of model 1
!
INTEGER,DIMENSION(2) :: IXY1          ! first point relatively to model 1
                                      ! corresponding to local physical domain 2 (local coords)
INTEGER,DIMENSION(2) :: IXY1GLB          ! first point relatively to model 1
                                      ! corresponding to global physical domain 2 (global coords)
INTEGER,DIMENSION(1) :: IX2,IY2       ! point relatively to model 2 corresponding
                                      ! to second physical point of model 1
INTEGER,DIMENSION(1) :: IXSUP1,IYSUP1 ! last point relatively to model 1
                                      ! corresponding to physical domain 2
REAL :: IXSUPCOORD1,IYSUPCOORD1  ! coordinates of the last point relatively to model 1
                     ! corresponding to physical domain 2
!
REAL                 :: ZEPS = 1.E-6  ! a small number
!
INTEGER     :: JI,JJ         ! loop controls relatively to model 2
INTEGER     :: JIBOX,JJBOX   ! grid mesh relatively to model 1
INTEGER     :: IINFO_ll
INTEGER     :: IROOTBUF
INTEGER     :: IROOT
INTEGER     :: IPROC
INTEGER     :: IXOR_F, IYOR_F    ! origin of local father subdomain (global coord)
INTEGER     :: IXEND_F, IYEND_F    ! end of local father subdomain (global coord)
!INTEGER     :: IXOR_C, IYOR_C    ! origin of local father subdomain (global coord)
!INTEGER     :: IXEND_C, IYEND_C    ! end of local father subdomain (global coord)
!INTEGER     :: IIMAX_C_ll, IJMAX_C_ll   ! global dimensions of child model
INTEGER     :: II
REAL        :: ZSENDBUF, ZRECVBUF
REAL        :: ZCOEF         ! ponderation coefficient for linear interpolation
REAL, DIMENSION(:), ALLOCATABLE :: ZXHAT, ZYHAT ! coordinates of model 2
!                            ! recomputed from coordinates of model 1 and ratios
REAL, DIMENSION(:), ALLOCATABLE :: ZPGDXHAT, ZPGDYHAT ! as XPGDXHAT and XPGDYHAT
!                                                     ! with one more point
REAL :: ZERROR_X,ZERROR_Y
REAL :: ZPGDXHATIXY1,ZPGDYHATIXY1         ! value of XPGDXHAT and XPGDYHAT at origin point of son model
REAL :: ZPGDXHATIXY1_1,ZPGDYHATIXY1_1     ! value of XPGDXHAT and XPGDYHAT at the next points in X and Y direction respectively
REAL :: ZXHATFIRSTENTRY_C,ZYHATFIRSTENTRY_C     ! value of XXHAT and XYHAT at the first physical point of son model
REAL :: ZXHATLASTENTRY_C,ZYHATLASTENTRY_C     ! value of XXHAT and XYHAT at the last physical point of son model
REAL :: ZPGDXHATIXY2,ZPGDYHATIXY2         ! value of XPGDXHAT and XPGDYHAT at end point of son model
REAL :: ZPGDXHATIXY2_1,ZPGDYHATIXY2_1     ! value of XPGDXHAT and XPGDYHAT at the next points in X and Y direction respectively
TYPE(ZONE_ll), DIMENSION(:), ALLOCATABLE :: TZSPLITTING
INTEGER, DIMENSION(2) :: IOR_C     ! position of pgd model origin points according to father (refered as model 1) domain / 0 if not on local father subdomain
!TYPE(LIST1D_ll), POINTER :: TZFIELDS_ll  ! list of fields to exchange
!
! variables needed for asynchronous communications
!INTEGER,PARAMETER                                     :: MPI_MAX_REQ = 1024
!INTEGER,SAVE,DIMENSION(MPI_MAX_REQ)                   :: REQ_TAB
!INTEGER                                               :: NB_REQ
!
!-------------------------------------------------------------------------------
! Current model is DAD model
!
! get splitting of father model
ALLOCATE(TZSPLITTING(ISNPROC))
CALL SPLIT2 ( NIMAX_TMP_ll, NJMAX_TMP_ll, 1, ISNPROC, TZSPLITTING, YSPLITTING )
IXOR_F = TZSPLITTING(ISP)%NXOR-JPHEXT
IYOR_F = TZSPLITTING(ISP)%NYOR-JPHEXT
IXEND_F = TZSPLITTING(ISP)%NXEND-JPHEXT
IYEND_F = TZSPLITTING(ISP)%NYEND-JPHEXT
!
! go to son model
CALL GOTO_MODEL(KMI)
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!! get global dims of son model
!IIMAX_C_ll = NXEND_ALL(KMI) - NXOR_ALL(KMI) - JPHEXT  !c'est bizarre mais on l'a init comme ca car sinon get_globaldims_ll donne un resultat faux...
!IJMAX_C_ll = NYEND_ALL(KMI) - NYOR_ALL(KMI) - JPHEXT  !c'est bizarre mais on l'a init comme ca car sinon get_globaldims_ll donne un resultat faux...
DEALLOCATE(TZSPLITTING)
!
ILUOUT = TLUOUT0%NLU
!20131008 adapt calculation NPGDIMAX AND NPGDSJMAX and IIU,IJU from retrieve1 !
!IIU=NPGDIMAX+2*JPHEXT
!IJU=NPGDJMAX+2*JPHEXT
IF ( CPROGRAM == 'REAL ' ) THEN
!IF ( CPROGRAM == 'REAL ' .OR. CPROGRAM == 'NESPGD' ) THEN
!20131009 adapt all changes from retrieve1
  XPGDLATOR=XLATORI
  XPGDLONOR=XLONORI
  NPGDIMAX =NIMAX
  NPGDJMAX =NJMAX
  IF (ALLOCATED(XPGDXHAT)) DEALLOCATE(XPGDXHAT)
  IF (ALLOCATED(XPGDYHAT)) DEALLOCATE(XPGDYHAT)
  ALLOCATE(XPGDXHAT(SIZE(XXHAT)))
  ALLOCATE(XPGDYHAT(SIZE(XYHAT)))
  XPGDXHAT(:)=XXHAT(:)
  XPGDYHAT(:)=XYHAT(:)
ELSE
!JUAN correction pour PREP_NEST_PGD 4/04/2014
!!$NPGDIMAX =NIMAX_TMP_ll
!!$NPGDJMAX =NJMAX_TMP_ll
ENDIF
!
!20131008 : now compute IIU & IJU
IIU=NPGDIMAX+2*JPHEXT
IJU=NPGDJMAX+2*JPHEXT
!
!
!*      1.    KXOR_C_ll,KYOR_C_ll
!             ---------
!
IF(.NOT.LCARTESIAN) THEN
!
!*      1.1   latitude and longitude of first local flux point (model2)
!             ---------------------------------------------------
!
  CALL SM_LATLON(XLATORI,XLONORI,                 &
                 XXHAT(JPHEXT+1),XYHAT(JPHEXT+1), &
                 ZLAT2,ZLON2)

!
!*      1.2   latitude and longitude of all local flux points (model1)
!             --------------------------------------------------
!
  ALLOCATE(ZPGDLAT1(IIU,IJU))
  ALLOCATE(ZPGDLON1(IIU,IJU))
  CALL SM_LATLON(XPGDLATOR,XPGDLONOR,                                 &
                 SPREAD(XPGDXHAT(:),2,IJU),SPREAD(XPGDYHAT(:),1,IIU), &
                 ZPGDLAT1(:,:),ZPGDLON1(:,:))
  CALL MPPDB_CHECK2D(ZPGDLAT1,"retrieve2_nest_info:ZPGDLAT1",PRECISION)
  CALL MPPDB_CHECK2D(ZPGDLON1,"retrieve2_nest_info:ZPGDLON1",PRECISION)
ENDIF  
!
!*      1.3   KXOR_C_ll, KYOR_C_ll - origin (global) of son model in father grid
!
!
! get origin of the intersection of father subdomain and son subdomain (in local coordinates)
! we do not differenciate case LCARTESIAN and the other cases
!
IXY1(1:1)=MINLOC(ABS(XPGDXHAT(:)-XXHAT(JPHEXT+1)))
IXY1(2:2)=MINLOC(ABS(XPGDYHAT(:)-XYHAT(JPHEXT+1)))
! check if there is an intersection
IF ( IXY1(1) == SIZE(XPGDXHAT) ) THEN
  IF ( XPGDXHAT(SIZE(XPGDXHAT)) <  XXHAT(JPHEXT+1) ) THEN
    ! there is no intersection - son subdomain is west of father subdomain
    IXY1(1) = 0
  ENDIF
ELSE IF ( IXY1(1) == 1 ) THEN
  IF ( XPGDXHAT(1) >  XXHAT(SIZE(XXHAT)-JPHEXT) ) THEN
    ! there is no intersection - son subdomain is east of father subdomain
    IXY1(1) = 0
  ENDIF
ENDIF
IF ( IXY1(2) == SIZE(XPGDYHAT) ) THEN
  IF ( XPGDYHAT(SIZE(XPGDYHAT)) <  XYHAT(JPHEXT+1) ) THEN
    ! there is no intersection - son subdomain is north of father subdomain
    IXY1(2) = 0
  ENDIF
ELSE IF ( IXY1(2) == 1 ) THEN
  IF ( XPGDYHAT(1) >  XYHAT(SIZE(XYHAT)-JPHEXT) ) THEN
    ! there is no intersection - son subdomain is south of father subdomain
    IXY1(2) = 0
  ENDIF
ENDIF
!
! Get the indices of the origin of global son model in father model (global coordinates) : KXOR_C_ll, KYOR_C_ll
!
! get the value of XXHAT and XYHAT at the origin of global son model
  ZXHATFIRSTENTRY_C = XXHAT(JPHEXT+1)
  ZYHATFIRSTENTRY_C = XYHAT(JPHEXT+1)
  CALL MPI_ALLREDUCE(XXHAT(JPHEXT+1), ZXHATFIRSTENTRY_C, 1,MPI_PRECISION, MPI_MIN, NMNH_COMM_WORLD, IINFO_ll)
  CALL MPI_ALLREDUCE(XYHAT(JPHEXT+1), ZYHATFIRSTENTRY_C, 1,MPI_PRECISION, MPI_MIN, NMNH_COMM_WORLD, IINFO_ll)
! get the latitude and longitude ZLAT2 and ZLON2 at the origin of global son model
  ZLAT2GLB = ZLAT2
  ZLON2GLB = ZLON2
  CALL MPI_ALLREDUCE(ZLAT2, ZLAT2GLB, 1,MPI_PRECISION, MPI_MIN, NMNH_COMM_WORLD, IINFO_ll)
  CALL MPI_ALLREDUCE(ZLON2, ZLON2GLB, 1,MPI_PRECISION, MPI_MIN, NMNH_COMM_WORLD, IINFO_ll)

  ! identify the process that own the origin of global son model, and communicate the global indices of the origin to all processes
  IF ( ZXHATFIRSTENTRY_C > XPGDXHAT(JPHEXT+1) .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) .AND. &
       ZYHATFIRSTENTRY_C > XPGDYHAT(JPHEXT+1) .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    IOR_C(1:1) = MINLOC(ABS(XPGDXHAT(:)-ZXHATFIRSTENTRY_C))
    IOR_C(2:2) = MINLOC(ABS(XPGDYHAT(:)-ZYHATFIRSTENTRY_C))
    IOR_C(1:1) = IOR_C(1:1) + IXOR_F - 1 - JPHEXT
    IOR_C(2:2) = IOR_C(2:2) + IYOR_F - 1 - JPHEXT
    ! we do some tests....
!    IF (LCARTESIAN ) THEN
      ZERROR_X=MINVAL(ABS(XPGDXHAT(:)-ZXHATFIRSTENTRY_C))
      ZERROR_Y=MINVAL(ABS(XPGDYHAT(:)-ZYHATFIRSTENTRY_C))
      IF ( ZERROR_X+ZERROR_Y > ZEPS ) THEN
	WRITE(ILUOUT,*) 'the first physical flux point of model ',KDAD,' does not correspond'
	WRITE(ILUOUT,*) 'to any of its father.'
	WRITE(ILUOUT,*) 'error on x and y : ', ZERROR_X,ZERROR_Y
    !callabortstop
    !CALL ABORT
    !    STOP
      END IF
!    ELSE
!      IF (MINVAL(ABS(ZPGDLAT1(:,:)-ZLAT2)+ABS(ZPGDLON1(:,:)-ZLON2))>ZEPS) THEN
!	WRITE(ILUOUT,*) 'the first physical flux point of model ',KDAD,' does not correspond'
!	WRITE(ILUOUT,*) 'to any of its father.'
!	WRITE(ILUOUT,*) 'sum of error on latitude and longitude: ', &
!		      MINVAL(ABS(ZPGDLAT1(:,:)-ZLAT2)+ABS(ZPGDLON1(:,:)-ZLON2))
    !callabortstop
    !CALL ABORT
    !    STOP
!      END IF
!    END IF
  ELSE
    IOR_C(1:1)=0
    IOR_C(2:2)=0
  ENDIF
  CALL MPI_ALLREDUCE(IOR_C(1:1), KXOR_C_ll, 1,MPI_INTEGER, MPI_SUM, NMNH_COMM_WORLD, IINFO_ll)
  CALL MPI_ALLREDUCE(IOR_C(2:2), KYOR_C_ll, 1,MPI_INTEGER, MPI_SUM, NMNH_COMM_WORLD, IINFO_ll)
!
!*      1.4   modify coordinates
! so that XXHAT(JPEXT+1) and XYHAT(JPEXT+1) correspond to the coordinates of the closest father grid points east (resp. north) of XXHAT(JPEXT+1) and XYHAT(JPEXT+1)
!             ------------------
!
! we need to do communications :
! each process must get the value of XPGDXHAT at the origin of its local son subdomain
!
!
! 1.4.1- Identify the process that owns the origin of local son model
!    we do not know the size of son domain in the father grid nor the global index of the origin of the local son subdmain,
!    so it is tricky.
!    We use the coordinates of the origin of local son model : XXHAT(JPHEXT+1) and XYHAT(JPHEXT+1)
! 1.4.2- communicate the values of XPGDXHAT and XPGDYHAT at the origin of local son model
DO IPROC = 0,ISNPROC-1  !loop on all processes
  ! XXHAT(JPHEXT+1), XYHAT(JPHEXT+1)  is the first physical entry of local son subdomain
  ZXHATFIRSTENTRY_C = XXHAT(JPHEXT+1)
  ZYHATFIRSTENTRY_C = XYHAT(JPHEXT+1)
  ! broadcast XXHAT(JPHEXT+1) and find which process' father subdomain contains the coords of the first physical entry of local son subdomain
  CALL MPI_BCAST( ZXHATFIRSTENTRY_C, 1, MPI_PRECISION, IPROC, NMNH_COMM_WORLD, IINFO_ll )
  CALL MPI_BCAST( ZYHATFIRSTENTRY_C, 1, MPI_PRECISION, IPROC, NMNH_COMM_WORLD, IINFO_ll )
  !
  ! communicating the value of XPGDXHAT (X direction) at the origin of local son subdomain
  IF (  IPROC == ISP-1 .AND. ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
    .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
    .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
    .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! in this case, the local father subdomain contains the first physical point of local son subdomain
    ZPGDXHATIXY1 = XPGDXHAT(IXY1(1))
    ZPGDYHATIXY1 = XPGDYHAT(IXY1(2))
  ELSE IF ( ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
  .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
  .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
  .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! the local father subdomain of current process contains the first physical point of local son subdomain of IPROC
    ! search for the first father physical grid point east and north of (not strictly) the first physical point of local son subdomain
    II=SIZE(XPGDXHAT)-JPHEXT
    DO WHILE ( XPGDXHAT(II) > ZXHATFIRSTENTRY_C )
      II=II-1
    END DO
    ! the index of the first physical point of the local son subdomain of IPROC is II on the current process
    ! send XPGDXHAT(II) to process IPROC
    ZSENDBUF = XPGDXHAT(II)
    CALL MPI_SEND( ZSENDBUF,1,MPI_PRECISION,IPROC,ISP+II,NMNH_COMM_WORLD,IINFO_ll )
  ELSE IF ( IPROC == ISP-1 ) THEN
    CALL MPI_RECV( ZRECVBUF,1,MPI_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,NMNH_COMM_WORLD,MPI_STATUS_IGNORE,IINFO_ll )
    ZPGDXHATIXY1 = ZRECVBUF
  ELSE
    ! the other processes do nothing...
  ENDIF
  !
  ! communicating the value of XPGDYHAT (Y direction) at the origin of local son subdomain
  IF (  IPROC == ISP-1 .AND. ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
    .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
    .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
    .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! in this case, the local father subdomain contains the first physical point of local son subdomain
    ZPGDXHATIXY1 = XPGDXHAT(IXY1(1))
    ZPGDYHATIXY1 = XPGDYHAT(IXY1(2))
  ELSE IF ( ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
  .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
  .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
  .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! the local father subdomain of current process contains the first physical point of local son subdomain
    ! search for the first father physical grid point east and north of (not strictly) the first physical point of local son subdomain
    II=SIZE(XPGDYHAT)-JPHEXT
    DO WHILE ( XPGDYHAT(II) > ZYHATFIRSTENTRY_C )
      II=II-1
    END DO
    ! the index of the first physical point of the local son subdomain is II on the current process
    ! send XPGDYHAT(II) to process IPROC
    ZSENDBUF = XPGDYHAT(II)
    CALL MPI_SEND( ZSENDBUF,1,MPI_PRECISION,IPROC,ISP+II,NMNH_COMM_WORLD,IINFO_ll )
  ELSE IF ( IPROC == ISP-1 ) THEN
    CALL MPI_RECV( ZRECVBUF,1,MPI_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,NMNH_COMM_WORLD,MPI_STATUS_IGNORE,IINFO_ll )
    ZPGDYHATIXY1 = ZRECVBUF
  ELSE
    ! the other processes do nothing...
  ENDIF
  ! REMARK :
  ! I have to do synchronous communications since the receiving process does not know the rank
  ! of the sending process, nor the tag of the message
ENDDO
!
! 1.4.3- communicate the values of XPGDXHAT (resp. XPGDYHAT) at the next point east (resp. north) of the origin of local son model
!     (same as for 1.4.2)
!
DO IPROC = 0,ISNPROC-1  !loop on all processes
  ZXHATFIRSTENTRY_C = XXHAT(JPHEXT+1)
  ZYHATFIRSTENTRY_C = XYHAT(JPHEXT+1)
  ! broadcast XXHAT(JPHEXT+1) and find which process' father subdomain contains the coords of the first physical entry of local son subdomain
  CALL MPI_BCAST( ZXHATFIRSTENTRY_C, 1, MPI_PRECISION, IPROC, NMNH_COMM_WORLD, IINFO_ll )
  CALL MPI_BCAST( ZYHATFIRSTENTRY_C, 1, MPI_PRECISION, IPROC, NMNH_COMM_WORLD, IINFO_ll )
  !
  ! communicating the value of XPGDXHAT (X direction) at the origin of local son subdomain
  IF (  IPROC == ISP-1 .AND. ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
    .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
    .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
    .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! in this case, the local father subdomain contains the first physical point of local son subdomain
    ZPGDXHATIXY1_1 = XPGDXHAT(IXY1(1)+1)
    ZPGDYHATIXY1_1 = XPGDYHAT(IXY1(2)+1)
  ELSE IF ( ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
  .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
  .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
  .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! the local father subdomain of current process contains the first physical point of local son subdomain
    ! search for the first father physical grid point east and north of (not strictly) the first physical point of local son subdomain
    II=SIZE(XPGDXHAT)-JPHEXT
    DO WHILE ( XPGDXHAT(II) > ZXHATFIRSTENTRY_C )
      II=II-1
    END DO
    ! the index of the first physical point of the local son subdomain is II on the current process
    ! XPGDXHAT(II+1) is also defined on current process since HALO is at least 1
    ! send XPGDXHAT(II+1) to process IPROC
    ZSENDBUF = XPGDXHAT(II+1)
    CALL MPI_SEND( ZSENDBUF,1,MPI_PRECISION,IPROC,ISP+II+1,NMNH_COMM_WORLD,IINFO_ll )
  ELSE IF ( IPROC == ISP-1 ) THEN
    CALL MPI_RECV( ZRECVBUF,1,MPI_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,NMNH_COMM_WORLD,MPI_STATUS_IGNORE,IINFO_ll )
    ZPGDXHATIXY1_1 = ZRECVBUF
  ELSE
    ! the other processes do nothing...
  ENDIF
  !
  ! communicating the value of XPGDYHAT (Y direction) at the origin of local son subdomain
  IF (  IPROC == ISP-1 .AND. ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
    .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
    .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
    .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! in this case, the local father subdomain contains the first physical point of local son subdomain
    ZPGDXHATIXY1_1 = XPGDXHAT(IXY1(1)+1)
    ZPGDYHATIXY1_1 = XPGDYHAT(IXY1(2)+1)
  ELSE IF ( ZXHATFIRSTENTRY_C >= XPGDXHAT(JPHEXT+1) &
  .AND. ZXHATFIRSTENTRY_C <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) &
  .AND. ZYHATFIRSTENTRY_C >= XPGDYHAT(JPHEXT+1) &
  .AND. ZYHATFIRSTENTRY_C <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
    ! the local father subdomain of current process contains the first physical point of local son subdomain
    ! search for the first father physical grid point east and north of (not strictly) the first physical point of local son subdomain
    II=SIZE(XPGDYHAT)-JPHEXT
    DO WHILE ( XPGDYHAT(II) > ZYHATFIRSTENTRY_C )
      II=II-1
    END DO
    ! the index of the first physical point of the local son subdomain is II on the current process
    ! XPGDYHAT(II+1) is also defined on current process since HALO is at least 1
    ! send XPGDYHAT(II+1) to process IPROC
    ZSENDBUF = XPGDYHAT(II+1)
    CALL MPI_SEND( ZSENDBUF,1,MPI_PRECISION,IPROC,ISP+II+1,NMNH_COMM_WORLD,IINFO_ll )
  ELSE IF ( IPROC == ISP-1 ) THEN
    CALL MPI_RECV( ZRECVBUF,1,MPI_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,NMNH_COMM_WORLD,MPI_STATUS_IGNORE,IINFO_ll )
    ZPGDYHATIXY1_1 = ZRECVBUF
  ELSE
    ! the other processes do nothing...
  ENDIF
ENDDO
!
! 1.4.4 - modify coordinates so that XXHAT(JPEXT+1) and XYHAT(JPEXT+1) correspond to the coordinates of the closest father grid points east (resp. north) of XXHAT(JPEXT+1) and XYHAT(JPEXT+1)
!
XXHAT(:) = XXHAT(:) + ZPGDXHATIXY1-XXHAT(JPHEXT+1)
XYHAT(:) = XYHAT(:) + ZPGDYHATIXY1-XYHAT(JPHEXT+1)
!XXHAT(:) = XXHAT(:) + XPGDXHAT(IXY1(1))-XXHAT(JPHEXT+1)
!XYHAT(:) = XYHAT(:) + XPGDYHAT(IXY1(2))-XYHAT(JPHEXT+1)
! 
!-------------------------------------------------------------------------------
!
!*      2.    KDXRATIO, KDYRATIO
!             ------------------
!
IX2(:)=MINLOC(ABS(ZPGDXHATIXY1_1-XXHAT(:)))
IY2(:)=MINLOC(ABS(ZPGDYHATIXY1_1-XYHAT(:)))
!
KDXRATIO=IX2(1)-JPHEXT-1
KDYRATIO=IY2(1)-JPHEXT-1
!
!-------------------------------------------------------------------------------
!
!*      3.    KXSIZE,KYSIZE
!             -------------

! 3.1- Identify the process that owns the end of local son model
!    we do not know the size of son domain in the father grid nor the global index of the end of the local son subdmain,
!    so it is tricky.
!    We use the coordinates of the origin of local son model : XXHAT(JPHEXT+1) and XYHAT(JPHEXT+1)
! 3.2- communicate the values of XPGDXHAT and XPGDYHAT at the point just past the end of local son model
! WARNING: we assume JPHEXT >= 1
DO IPROC = 0,ISNPROC-1  !loop on all processes
  ZXHATLASTENTRY_C = XXHAT(SIZE(XXHAT)-JPHEXT)
  ZYHATLASTENTRY_C = XYHAT(SIZE(XYHAT)-JPHEXT)
  ! broadcast XXHAT(SIZE(XXHAT)-JPHEXT) and find which process' father subdomain contains the coords of the last physical entry of local son subdomain
  CALL MPI_BCAST( ZXHATLASTENTRY_C, 1, MPI_PRECISION, IPROC, NMNH_COMM_WORLD, IINFO_ll )
  CALL MPI_BCAST( ZYHATLASTENTRY_C, 1, MPI_PRECISION, IPROC, NMNH_COMM_WORLD, IINFO_ll )
  !
  ! communicating the value of XPGDXHAT (X direction) at the origin of local son subdomain
  IF (  IPROC == ISP-1 .AND. ZXHATLASTENTRY_C >= XPGDXHAT(JPHEXT+1) &
    .AND. ZXHATLASTENTRY_C < XPGDXHAT(SIZE(XPGDXHAT)) &
    .AND. ZYHATLASTENTRY_C >= XPGDYHAT(JPHEXT+1) &
    .AND. ZYHATLASTENTRY_C < XPGDYHAT(SIZE(XPGDYHAT)) ) THEN
    ! the local father subdomain of current process contains the last physical point of local son subdomain
    ! search for the last father physical grid point west and south of (not strictly) the last physical point of local son subdomain
    II=SIZE(XPGDXHAT)-JPHEXT
    DO WHILE ( XPGDXHAT(II) > ZXHATLASTENTRY_C )
      II=II-1
    END DO
    ! the index of the last physical point of the local son subdomain is II on the current process
    ! XPGDYHAT(II+1) is also defined on current process since HALO is at least 1
    ZPGDXHATIXY2_1 = XPGDXHAT(II)
  ELSE IF ( ZXHATLASTENTRY_C >= XPGDXHAT(JPHEXT+1) &
  .AND. ZXHATLASTENTRY_C < XPGDXHAT(SIZE(XPGDXHAT)) &
  .AND. ZYHATLASTENTRY_C >= XPGDYHAT(JPHEXT+1) &
  .AND. ZYHATLASTENTRY_C < XPGDYHAT(SIZE(XPGDYHAT)) ) THEN
    ! the local father subdomain of current process contains the last physical point of local son subdomain
    ! search for the last father physical grid point west and south of (not strictly) the last physical point of local son subdomain
    II=SIZE(XPGDXHAT)-JPHEXT
    DO WHILE ( XPGDXHAT(II) > ZXHATLASTENTRY_C )
      II=II-1
    END DO
    ! the index of the last physical point of the local son subdomain is II on the current process
    ! send XPGDXHAT(II) to process IPROC
    ! XPGDYHAT(II+1) is also defined on current process since HALO is at least 1
    ZSENDBUF = XPGDXHAT(II)
    CALL MPI_SEND( ZSENDBUF,1,MPI_PRECISION,IPROC,ISP+II,NMNH_COMM_WORLD,IINFO_ll )
  ELSE IF ( IPROC == ISP-1 ) THEN
    CALL MPI_RECV( ZRECVBUF,1,MPI_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,NMNH_COMM_WORLD,MPI_STATUS_IGNORE,IINFO_ll )
    ZPGDXHATIXY2_1 = ZRECVBUF
  ELSE
    ! the other processes do nothing...
  ENDIF
  !
  ! communicating the value of XPGDYHAT (Y direction) at the origin of local son subdomain
  IF (  IPROC == ISP-1 .AND. ZXHATLASTENTRY_C >= XPGDXHAT(JPHEXT+1) &
    .AND. ZXHATLASTENTRY_C < XPGDXHAT(SIZE(XPGDXHAT)) &
    .AND. ZYHATLASTENTRY_C >= XPGDYHAT(JPHEXT+1) &
    .AND. ZYHATLASTENTRY_C < XPGDYHAT(SIZE(XPGDYHAT)) ) THEN
    ! the local father subdomain of current process contains the last physical point of local son subdomain
    ! search for the last father physical grid point west and south of (not strictly) the last physical point of local son subdomain
    II=SIZE(XPGDYHAT)-JPHEXT
    DO WHILE ( XPGDYHAT(II) > ZYHATLASTENTRY_C )
      II=II-1
    END DO
    ! the index of the last physical point of the local son subdomain is II on the current process
    ! send XPGDYHAT(II) to process IPROC
    ZPGDYHATIXY2_1 = XPGDYHAT(II)
  ELSE IF ( ZXHATLASTENTRY_C >= XPGDXHAT(JPHEXT+1) &
  .AND. ZXHATLASTENTRY_C < XPGDXHAT(SIZE(XPGDXHAT)) &
  .AND. ZYHATLASTENTRY_C >= XPGDYHAT(JPHEXT+1) &
  .AND. ZYHATLASTENTRY_C < XPGDYHAT(SIZE(XPGDYHAT)) ) THEN
    ! the local father subdomain of current process contains the last physical point of local son subdomain
    ! search for the last father physical grid point west and south of (not strictly) the last physical point of local son subdomain
    II=SIZE(XPGDYHAT)-JPHEXT
    DO WHILE ( XPGDYHAT(II) > ZYHATLASTENTRY_C )
      II=II-1
    END DO
    ! the index of the last physical point of the local son subdomain is II on the current process
    ! send XPGDYHAT(II) to process IPROC
    ZSENDBUF = XPGDYHAT(II)
    CALL MPI_SEND( ZSENDBUF,1,MPI_PRECISION,IPROC,ISP+II,NMNH_COMM_WORLD,IINFO_ll )
  ELSE IF ( IPROC == ISP-1 ) THEN
    CALL MPI_RECV( ZRECVBUF,1,MPI_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,NMNH_COMM_WORLD,MPI_STATUS_IGNORE,IINFO_ll )
    ZPGDYHATIXY2_1 = ZRECVBUF
  ELSE
    ! the other processes do nothing...
  ENDIF
ENDDO
  ! REMARK :
  ! I have to do synchronous communications since the receiving process does not know the rank
  ! of the sending process, nor the tag of the message
  ! For the same reason (tag unknown to receiving process),
  ! I cannot send/recv XPGDXHAT(II) and XPGDYHAT(II) at the same time

! 3.3 - now we have the coordinates (ZPGDXHATIXY2_1, ZPGDYHATIXY2_1) of the point in father grid just right+north of the LOCAL son subdomain
!       We compute the coordinates of the last point in father grid of the GLOBAL son subdomain
CALL MPI_ALLREDUCE(ZPGDXHATIXY2_1, IXSUPCOORD1, 1,MPI_PRECISION, MPI_MAX, NMNH_COMM_WORLD, IINFO_ll)
CALL MPI_ALLREDUCE(ZPGDYHATIXY2_1, IYSUPCOORD1, 1,MPI_PRECISION, MPI_MAX, NMNH_COMM_WORLD, IINFO_ll)

!     we compute the index of this point in local father grid
IF ( IXSUPCOORD1 >= XPGDXHAT(1+JPHEXT) .AND. IXSUPCOORD1 <= XPGDXHAT(SIZE(XPGDXHAT)-JPHEXT) .AND. &
   IYSUPCOORD1 >= XPGDYHAT(1+JPHEXT) .AND. IYSUPCOORD1 <= XPGDYHAT(SIZE(XPGDYHAT)-JPHEXT) ) THEN
  ! the point in father grid just right+north of the local son subdomain is in local subdomain
  ! compute the local index in X (resp. Y) direction of this point
  IXSUP1(:)=1
  DO WHILE( XPGDXHAT(IXSUP1(1)) < IXSUPCOORD1 )
    IXSUP1(:)=IXSUP1(:)+1
  ENDDO
  IYSUP1(:)=1
  DO WHILE( XPGDYHAT(IYSUP1(1)) < IYSUPCOORD1 )
    IYSUP1(:)=IYSUP1(:)+1
  ENDDO
  ! switch to global coordinates
  IXSUP1(:) = IXSUP1(:) + IXOR_F - 1
  IYSUP1(:) = IYSUP1(:) + IYOR_F - 1
ELSE
  IXSUP1(:)=0
  IYSUP1(:)=0
ENDIF
CALL MPI_ALLREDUCE(IXSUP1(1), KXSIZE, 1,MPI_INTEGER, MPI_MAX, NMNH_COMM_WORLD, IINFO_ll)
CALL MPI_ALLREDUCE(IYSUP1(1), KYSIZE, 1,MPI_INTEGER, MPI_MAX, NMNH_COMM_WORLD, IINFO_ll)
IXSUP1(1) = KXSIZE
IYSUP1(1) = KYSIZE
!
! compute the global size of son model in the father grid
KXSIZE=IXSUP1(1)-(KXOR_C_ll+JPHEXT)+1
KYSIZE=IYSUP1(1)-(KYOR_C_ll+JPHEXT)+1
!
! some more tests
!
CALL MPI_ALLREDUCE(IIU-2*JPHEXT, IIUGLB, 1,MPI_INTEGER, MPI_SUM, NMNH_COMM_WORLD, IINFO_ll)
CALL MPI_ALLREDUCE(IJU-2*JPHEXT, IJUGLB, 1,MPI_INTEGER, MPI_SUM, NMNH_COMM_WORLD, IINFO_ll)
IIUGLB = IIUGLB + 2*JPHEXT
IJUGLB = IJUGLB + 2*JPHEXT
IF (     KXOR_C_ll<1 .OR. KXOR_C_ll+KXSIZE+2*JPHEXT>IIUGLB      &
    .OR. KYOR_C_ll<1 .OR. KYOR_C_ll+KYSIZE+2*JPHEXT>IJUGLB) THEN
  WRITE(ILUOUT,*) 'KXEND or KYEND (last point used in domain',KMI,') outside of the domain'
  WRITE(ILUOUT,*) 'KXEND= ', KXOR_C_ll+KXSIZE+2*JPHEXT-1, 'KYEND= ', KYOR_C_ll+KYSIZE+2*JPHEXT-1
 !callabortstop
!CALL ABORT
!  STOP
END IF
!-------------------------------------------------------------------------------
!
!*      4.    Tests on coordinate arrays
!             --------------------------
!
ALLOCATE(ZXHAT(NIMAX+2*JPHEXT))
ALLOCATE(ZYHAT(NJMAX+2*JPHEXT))
!
IPGDIU = NPGDIMAX+2*JPHEXT
IPGDJU = NPGDJMAX+2*JPHEXT
!
ALLOCATE(ZPGDXHAT(0:IPGDIU+1))
ALLOCATE(ZPGDYHAT(0:IPGDJU+1))
!
! it is too complicated to test on the HALO
! it would require communications to determine the neighbouring processes
! and updating the extra halo points we added in ZPGDXHAT / ZPGDYHAT
!
ZPGDXHAT(1:IPGDIU) = XPGDXHAT(:)
ZPGDYHAT(1:IPGDJU) = XPGDYHAT(:)
!IF ( LEAST_ll() ) THEN
ZPGDXHAT(IPGDIU+1) = 2.* XPGDXHAT(IPGDIU) - XPGDXHAT(IPGDIU-1)
!ENDIF
!IF ( LNORTH_ll() ) THEN
ZPGDYHAT(IPGDJU+1) = 2.* XPGDYHAT(IPGDJU) - XPGDYHAT(IPGDJU-1)
!ENDIF
!IF ( LWEST_ll() ) THEN
ZPGDXHAT(0)        = 2.* XPGDXHAT(1) - XPGDXHAT(2)
!ENDIF
!IF ( LSOUTH_ll() ) THEN
ZPGDYHAT(0)        = 2.* XPGDYHAT(1) - XPGDYHAT(2)
!ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! JE COMMENTE TOUTE LA PARTIE 4 CAR IL S'AGIT SEULEMENT DE TESTS,
!!! ET POUR LES FAIRE CORRECTEMENT IL FAUT FAIRE DES COMMUNICATIONS : C'EST INUTILE DE LE FAIRE SYSTEMATIQUEMENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
DO JI=1,NIMAX+2*JPHEXT
  JIBOX=(JI+KDXRATIO-1-JPHEXT)/KDXRATIO + KXOR_C_ll
  ZCOEF= FLOAT(MOD(JI+KDXRATIO-1-JPHEXT,KDXRATIO))/FLOAT(KDXRATIO)
  ZXHAT(JI)=(1.-ZCOEF)*ZPGDXHAT(JIBOX+JPHEXT-1)+ZCOEF*ZPGDXHAT(JIBOX+JPHEXT) ! +1
END DO
!
DO JJ=1,NJMAX+2*JPHEXT
  JJBOX=(JJ+KDYRATIO-1-JPHEXT)/KDYRATIO + KYOR_C_ll
  ZCOEF= FLOAT(MOD(JJ+KDYRATIO-1-JPHEXT,KDYRATIO))/FLOAT(KDYRATIO)
  ZYHAT(JJ)=(1.-ZCOEF)*ZPGDYHAT(JJBOX+JPHEXT-1)+ZCOEF*ZPGDYHAT(JJBOX+JPHEXT) ! +1
END DO
!
IF (     ANY(ABS(XXHAT(:)-ZXHAT(:))>ZEPS)            &
    .OR. ANY(ABS(XYHAT(:)-ZYHAT(:))>ZEPS) ) THEN
  WRITE(ILUOUT,*) 'XHAT or YHAT functions are incoherent'
  WRITE(ILUOUT,*) ' '
  DO JI=1,NIMAX+2*JPHEXT
    WRITE(ILUOUT,*) '  XXHAT(',JI,')  = ', XXHAT(JI) , &
                    '  ZXHAT(',JI,')  = ', ZXHAT(JI)
  END DO
  WRITE(ILUOUT,*) ' '
  DO JJ=1,NJMAX+2*JPHEXT
    WRITE(ILUOUT,*) '  XYHAT(',JJ,')  = ', XYHAT(JJ) , &
                    '  ZYHAT(',JJ,')  = ', ZYHAT(JJ)
  END DO
 !callabortstop
!CALL ABORT
!  STOP
END IF
#endif
!
DEALLOCATE(ZXHAT)
DEALLOCATE(ZYHAT)
!
DEALLOCATE(ZPGDXHAT)
DEALLOCATE(ZPGDYHAT)
!
IF (.NOT. LCARTESIAN) THEN
  DEALLOCATE(ZPGDLAT1)
  DEALLOCATE(ZPGDLON1)
ENDIF
  
!-------------------------------------------------------------------------------
!
IF (KDAD > 0) CALL GOTO_MODEL(KDAD)
!
END SUBROUTINE RETRIEVE2_NEST_INFO_n
