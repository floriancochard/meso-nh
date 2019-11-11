!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------

#ifdef MNH_MPI_DOUBLE_PRECISION
#define MNH_MPI_REAL MPI_DOUBLE_PRECISION
#define MNH_MPI_2REAL MPI_2DOUBLE_PRECISION
#else
#define MNH_MPI_REAL MPI_REAL
#define MNH_MPI_2REAL MPI_2REAL
#endif

!       ########################
MODULE MODE_SPLITTINGZ_ll
  !     ########################
  !
  !!    Purpose
  !!    -------
  !
  !     The purpose of this module is to provide subroutines 
  !     for the splitting of the domain in X Y & Z direction 
  !
  !!    Routines Of The User Interface
  !!    ------------------------------
  !
  !     None
  !
  !!    Reference
  !!    ---------
  !
  !     User Interface for Meso-NH parallel package
  !     Ph. Kloos, L. Giraud, R. Guivarch, D. Lugato
  !
  !!    Authors
  !!    -------
  !
  !     J. ESCOBAR                 * LA *
  ! Modifications:
  !  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
  USE MODD_MPIF
  !
  USE MODE_SPLITTING_ll
  !
  LOGICAL,SAVE :: LINI_PARAZ = .FALSE. !Useful to know if a call to INI_PARAZ_ll or INI_PARAZ_CHILD_ll has already be done
  !
CONTAINS
  !       ################################
  SUBROUTINE INI_PARAZ_ll(KINFO_ll)
    !     ################################
    !
    !!****  *INI_PARA_ll* - routine to initialize the parallel variables
    !!
    !!    Purpose
    !!    -------
    !     the purpose of the routine is to fill the structured type variables
    !     TCRRT_PROCONF and TCRRT_COMDATA
    !
    !!**  Method
    !!    ------
    !!
    !!    External
    !!    --------
    !     Module MODE_SPLITTING_ll
    !      SPLIT2
    !
    !     Module MODE_CONSTRUCT_ll
    !       INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS,
    !       CONSTRUCT_HALO1, CONSTRUCT_HALO2, CONSTRUCT_HALO_EXTENDED,
    !       CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY,
    !       COMPUTE_HALO_MAX, COMPUTE_TRANS_MAX
    !
    !     Module MODE_NEST_ll
    !       INI_CHILD
    !
    !!    Implicit Arguments
    !!    ------------------
    !     Module MODD_DIM_ll
    !       JPHEXT - Horizontal External points number
    !       NDXRATIO_ALL, NDYRATIO_ALL, NXOR_ALL, NYOR_ALL,
    !       NXEND_ALL, NYEND_ALL,...
    ! 
    !     Module MODD_PARALLEL
    !       TCRRT_PROCONF - Current configuration for current model
    !       TCRRT_COMDATA - Current communication data structure for current model
    !                       and local processor
    ! 
    !     Reference
    !!    ---------
    ! 
    !!    AUTHOR
    !!    ------
    !       R. Guivarch
    ! 
    !!    MODIFICATIONS
    !!    -------------
    !     Original 01/05/98
    !     R. Guivarch 01/01/98  Grid-Nesting
    !     R. Guivarch 29/11/99  x and y splitting -> YSPLITTING 
    !     J. Escobar  24/09/2013 : temp patch for problem of gridnesting with different SHAPE
    !     M.Moge      10/02/2015 construct halo extended (needed for an interpolation in SPAWNING)
    !     J. Escobar 5/06/2018 : add cpp key MNH_USE_MPI_STATUSES_IGNORE for use of true MPI_STATUSES_IGNORE
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_DIM_ll
    USE MODD_PARAMETERS_ll
    USE MODD_STRUCTURE_ll
    USE MODD_VAR_ll
    !
    USE MODE_SPLITTING_ll, ONLY : SPLIT2
    !
    USE MODE_CONSTRUCT_ll, ONLY : INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS, &
         CONSTRUCT_HALO1, CONSTRUCT_HALO2, CONSTRUCT_HALO_EXTENDED, &
         CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY, &
         COMPUTE_HALO_MAX, COMPUTE_TRANS_MAX
    !
    USE MODE_NEST_ll, ONLY : INI_CHILD

    USE MODE_TOOLSZ_ll, ONLY : SPLITZ, ini_pzz, ini_boundariesz, ini_ezz, construct_transz
    !
    !JUANZ
    USE  MODE_MNH_WORLD , ONLY :  INIT_NMNH_COMM_WORLD
    USE  MODD_CONFZ     , ONLY :  NZ_VERB,NZ_PROC,MPI_BUFFER_SIZE,LMNH_MPI_ALLTOALLV_REMAP,NZ_SPLITTING
    !JUANZ
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(OUT) :: KINFO_ll
    !
    !*       0.2   declarations of local variables
    !
    !INTEGER  ,PARAMETER                      :: MPI_BUFFER_SIZE = 140000000
    CHARACTER,SAVE,ALLOCATABLE,DIMENSION(:)  :: MPI_BUFFER
    !JUAN
    LOGICAL,SAVE                             :: GFIRSTCALL = .TRUE.
    !JUAN

    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP ! intermediate zone
    !
    TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
    TYPE(PROCONF_ll), POINTER :: TZPROCONF
    INTEGER :: JMODEL
    INTEGER     :: IRESP
    LOGICAL     :: GISINIT
    !
    !JUAN
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP1_YP2_Z ! intermediate Full Z = B splitting  without halo zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SX_YP2_ZP1 ! intermediate Full X     splitting zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_Y_ZP1 ! intermediate Full Y     splitting zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_YP1_Z ! intermediate Full Z = B transposed splitting  without halo zone

    INTEGER :: JX_DOMAINS,JY_DOMAINS
    LOGICAL :: LPREM
    INTEGER :: P1,P2
    !JUANZ
    INTEGER :: P1P2(2), P1P2COORD(2) , IROW , ICOL, NROW, NCOL
    LOGICAL :: Lperiodic(2), remain_dims(2) , Lreorder
    INTEGER :: JI
    !JUANZ
    !JUAN
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    INITIALIZE MPI :
    !              --------------
    !
    KINFO_ll = 0
    CALL INIT_NMNH_COMM_WORLD(KINFO_ll)
    !
    CALL MPI_COMM_DUP(NMNH_COMM_WORLD, NHALO_COM, KINFO_ll)
    !
    CALL MPI_COMM_DUP(NMNH_COMM_WORLD, NHALO2_COM, KINFO_ll)
    !
    CALL MPI_COMM_DUP(NMNH_COMM_WORLD, NTRANS_COM, KINFO_ll)
    !
    CALL MPI_COMM_DUP(NMNH_COMM_WORLD, NGRID_COM, KINFO_ll)
    !
    MPI_PRECISION  = MNH_MPI_REAL
    MPI_2PRECISION = MNH_MPI_2REAL
    !
    ! For bug with intelmpi+ilp64+i8 declare MNH_STATUSES_IGNORE
    !
#ifndef MNH_USE_MPI_STATUSES_IGNORE
    ALLOCATE(MNH_STATUSES_IGNORE(MPI_STATUS_SIZE,NPROC*2))
#endif
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.    SET OUTPUT FILE :
    !              ---------------

    !  CALL OPEN_ll(UNIT=NIOUNIT,FILE=YOUTPUTFILE,ACTION='write',form&
    !       &='FORMATTED',MODE=SPECIFIC,IOSTAT=IRESP)
    !
    !-------------------------------------------------------------------------------
    !
    !*       3.    ALLOCATION :
    !              ----------
    !
    IF (GFIRSTCALL) THEN
       IF (LMNH_MPI_ALLTOALLV_REMAP) NZ_SPLITTING = 14
       IF (NZ_VERB .GE. 5 ) THEN
          IF ( IP .EQ. 1 ) THEN
             print*," INI_PARAZ_ll : MPI_BUFFER_SIZE in MB =",MPI_BUFFER_SIZE
             IF (LMNH_MPI_ALLTOALLV_REMAP)  print*," INI_PARAZ_ll : USING MPI_ALLTOALLV FOR REMAPPING "
          END IF
       END IF
       ALLOCATE(MPI_BUFFER(MPI_BUFFER_SIZE*1000000))
       CALL MPI_BUFFER_ATTACH(MPI_BUFFER,MPI_BUFFER_SIZE*1000000,KINFO_ll)
       GFIRSTCALL = .FALSE.
    ENDIF

    ALLOCATE(TZDZP(NPROC))
    !JUAN
    ALLOCATE(TZDZP_SXP1_YP2_Z(NPROC))
    ALLOCATE(TZDZP_SXP2_YP1_Z(NPROC))
    ALLOCATE(TZDZP_SX_YP2_ZP1(NPROC))
    ALLOCATE(TZDZP_SXP2_Y_ZP1(NPROC))
    !JUAN
    !
    ALLOCATE(TCRRT_PROCONF)
    CALL ALLOC(TCRRT_COMDATA)
    ALLOCATE(TCRRT_PROCONF%TSPLITS_B(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_X(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_Y(NPROC))
    !JUAN
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1(NPROC))
    !JUAN
    ALLOCATE(TCRRT_PROCONF%TBOUND(NPROC))
    NULLIFY(TCRRT_PROCONF%TPARENT)
    NULLIFY(TCRRT_COMDATA%TPARENT)
    NULLIFY(TCRRT_PROCONF%TCHILDREN)
    NULLIFY(TCRRT_COMDATA%TCHILDREN)
    !
    !-------------------------------------------------------------------------------
    !
    !*       4.    SPLITTING OF THE DOMAIN :
    !              -----------------------
    !
    DIMX = NIMAX_TMP_ll + 2*JPHEXT
    DIMY = NJMAX_TMP_ll + 2*JPHEXT
    DIMZ = NKMAX_TMP_ll + 2*JPVEXT
    !
    TCRRT_PROCONF%NUMBER = 1
    !

    !JUAN CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SXP2_YP1_Z,'BSPLITTING',NZ_PROC)
!!$    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SXP2_YP1_Z,'BSPLITTING',1)
!!$    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING',NZ_PROC)
!!$    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING',NZ_PROC)
    ! Add halo directly in Z direction 



    !
    ! find the B spltting
    !
    CALL DEF_SPLITTING2(JX_DOMAINS,JY_DOMAINS,NIMAX_TMP_ll,NJMAX_TMP_ll,NPROC,LPREM)
    !
    P1 = JX_DOMAINS
    IF (DIMZ .NE. 3 )    P1 = MIN(DIMZ,JX_DOMAINS)
    IF (NZ_PROC .GT. 0 ) P1 = NZ_PROC
    P2 = NPROC / P1
    !JUAN PATCH NESTING DIFFERENT SHAPE
    NZ_PROC = P1
    IF (NZ_VERB .GE. 5 ) THEN
       IF ( IP .EQ. 1 )THEN
          print*," INI_PARAZ_ll:: NZ_PROC   =",NZ_PROC
          print*," INI_PARAZ_ll:: JX_DOMAINS=",JX_DOMAINS
          print*," INI_PARAZ_ll:: JY_DOMAINS=",JY_DOMAINS
          print*
          !
          print*," INI_PARAZ_ll:: P1=MIN(NZ_PROC,DIMZ) > 0 .OR. MIN(DIMZ,MAX(JX_DOMAINS,JY_DOMAINS))=",  P1
          !
          print*," INI_PARAZ_ll:: P2=NPROC/P1/                                            =",  P2
       END IF
    END IF
    NP1 = P1
    NP2 = P2
    !
    !JUANZ
    P1P2(1) = NP2
    P1P2(2) = NP1
    Lperiodic(1) = .false.
    Lperiodic(2) = .false.
    Lreorder=.false.
    ! creating cartesian processor grid
    call MPI_Cart_create(NMNH_COMM_WORLD,2,P1P2,Lperiodic,Lreorder,NMNH_P1P2_WORLD,KINFO_ll)
    ! Obtaining process ids with in the cartesian grid
    call MPI_Cart_coords(NMNH_P1P2_WORLD,IP-1,2,P1P2COORD,KINFO_ll)
   
    ! using cart comworld create east-west(row) sub comworld
    remain_dims(1) = .false.
    remain_dims(2) = .true.
    call MPI_Cart_sub(NMNH_P1P2_WORLD,remain_dims,NMNH_ROW_WORLD,KINFO_ll)
    CALL MPI_COMM_RANK(NMNH_ROW_WORLD, IROW, KINFO_ll)
    CALL MPI_COMM_SIZE(NMNH_ROW_WORLD, NROW, KINFO_ll)

    ! using cart comworld create north-south(column) sub comworld
    remain_dims(1) = .true.
    remain_dims(2) = .false.
    call MPI_Cart_sub(NMNH_P1P2_WORLD,remain_dims,NMNH_COL_WORLD,KINFO_ll)
    CALL MPI_COMM_RANK(NMNH_COL_WORLD, ICOL, KINFO_ll)
    CALL MPI_COMM_SIZE(NMNH_COL_WORLD, NCOL, KINFO_ll)
    !JUANZ 

    CALL SPLIT2(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP,YSPLITTING,P1,P2)

    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,DIMZ,NPROC,TZDZP_SXP1_YP2_Z,'P1P2SPLITT', 1 ,P1,P2)
    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,DIMZ,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING', P1,P1,P2)
    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,DIMZ,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING', P1,P1,P2)
    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,DIMZ,NPROC,TZDZP_SXP2_YP1_Z,'P2P1SPLITT', 1 ,P1,P2)


    !JUAN    CALL SPLIT2(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP,YSPLITTING)


    !    
    !-------------------------------------------------------------------------------
    !
    !*       5.    INITIALIZATION OF TCRRT_PROCONF :
    !              -------------------------------
    !
    CALL INI_PZ(TCRRT_PROCONF,TZDZP)
    !JUAN
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z,TZDZP_SXP1_YP2_Z)
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z,TZDZP_SXP2_YP1_Z)
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1,TZDZP_SX_YP2_ZP1)
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1,TZDZP_SXP2_Y_ZP1)
    !JUAN
    !
    CALL INI_BOUNDARIES(TCRRT_PROCONF)
    !JUAN
    CALL INI_BOUNDARIESZ(TCRRT_PROCONF)
    !JUAN
    !
    CALL INI_EZ(TCRRT_PROCONF)
    !JUAN
    CALL INI_EZZ(TCRRT_PROCONF)
    !JUAN
    !
    CALL INI_TRANS(TCRRT_PROCONF)
    !
    !-------------------------------------------------------------------------------
    !
    !*       6.    INITIALIZATION OF TCRRT_COMDATA :
    !              -------------------------------
    !
    !*       6.1    Model Number
    !
    TCRRT_COMDATA%NUMBER = 1
    !
    !*       6.2    Pointer from TCRRT_COMDATA to TCRRT_PROCONF for 2Way splitting
    !
    TCRRT_COMDATA%TSPLIT_B => TCRRT_PROCONF%TSPLITS_B(IP)

    !TZSPLIT => TCRRT_COMDATA%TSPLIT_B
    !
    !
    !*       6.3   Pointer from TCRRT_COMDATA to TCRRT_PROCONF
    !        for x-slices splitting

    TCRRT_COMDATA%TSPLIT_X => TCRRT_PROCONF%TSPLITS_X(IP)
    !
    !TZSPLIT => TCRRT_COMDATA%TSPLIT_X
    !
    !
    !*       6.4   Pointer from TCRRT_COMDATA to TCRRT_PROCONF
    !              for y-slices splitting
    !
    TCRRT_COMDATA%TSPLIT_Y => TCRRT_PROCONF%TSPLITS_Y(IP)
    !
    !TZSPLIT => TCRRT_COMDATA%TSPLIT_Y
    !
    !JUAN
    DO JI=1, NPROC
       IF ( TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z(JI)%NUMBER .EQ. IP ) THEN 
          TCRRT_COMDATA%TSPLIT_SXP1_YP2_Z => TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z(JI)
       ENDIF
       IF ( TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z(JI)%NUMBER .EQ. IP ) THEN 
          TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z => TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z(JI)
       ENDIF
       IF (  TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1(JI)%NUMBER .EQ. IP ) THEN 
          TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1 => TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1(JI)
       ENDIF
       IF ( TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1(JI)%NUMBER .EQ. IP ) THEN
          TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1 => TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1(JI)
       END IF
    END DO
    !JUAN
    !
    !*       6.5   Construction of HALO1 communication data
    !
    CALL CONSTRUCT_HALO1(TCRRT_COMDATA, TCRRT_PROCONF)
    CALL CONSTRUCT_HALO2(TCRRT_COMDATA, TCRRT_PROCONF)
    CALL CONSTRUCT_HALO_EXTENDED(TCRRT_COMDATA, TCRRT_PROCONF, JPHEXT+1)
    !
    !
    !*       6.6   Construction of 1D communication data
    !
    ALLOCATE(TCRRT_COMDATA%HALO1DX)
    ALLOCATE(TCRRT_COMDATA%HALO1DX%NSEND_WEST(NPROC))
    ALLOCATE(TCRRT_COMDATA%HALO1DX%NSEND_EAST(NPROC))
    CALL CONSTRUCT_1DX(TCRRT_COMDATA, TCRRT_PROCONF)
    !
    ALLOCATE(TCRRT_COMDATA%HALO1DY)
    ALLOCATE(TCRRT_COMDATA%HALO1DY%NSEND_SOUTH(NPROC))
    ALLOCATE(TCRRT_COMDATA%HALO1DY%NSEND_NORTH(NPROC))
    CALL CONSTRUCT_1DY(TCRRT_COMDATA, TCRRT_PROCONF)
    !
    !
    !*       6.7   Construction of Transposition communication data
    !
    CALL CONSTRUCT_TRANS(TCRRT_COMDATA, TCRRT_PROCONF)
    CALL CONSTRUCT_TRANSZ(TCRRT_COMDATA, TCRRT_PROCONF)
    !
    !
    !-------------------------------------------------------------------------------
    !
    !        7.    GRID NESTING :
    !              ------------
    !
    NULLIFY(TCRRT_PROCONF%TCHILDREN)
    NULLIFY(TCRRT_COMDATA%TCHILDREN)
    NULLIFY(TCRRT_COMDATA%TP2C_DATA)
    !
    DO JMODEL = 1, JPMODELMAX
       !
       IF( NDAD(JMODEL) .EQ. TCRRT_PROCONF%NUMBER ) THEN
          CALL INI_CHILD(TCRRT_PROCONF, TCRRT_COMDATA, JMODEL)
       ENDIF
       !
    ENDDO
    !
    !-------------------------------------------------------------------------------
    !
    TZPROCONF => TCRRT_PROCONF
    !
    CALL COMPUTE_TRANS_MAX(NBUFFERSIZE_3D, TCRRT_COMDATA)
    IF (NZ_VERB .GE. 5 ) THEN
       IF (IP.EQ.1) print*,"INI_PARAZ_ll::COMPUTE_TRANS_MAX(NBUFFERSIZE_3D, TCRRT_COMDATA)=",NBUFFERSIZE_3D
    END IF
    !JUAN NCOMBUFFSIZE1 = NBUFFERSIZE_3D
    !NCOMBUFFSIZE1 = NBUFFERSIZE_3D*2
    NCOMBUFFSIZE1 = NBUFFERSIZE_3D
    !JUAN NCOMBUFFSIZE1 = 10000000
    !
    CALL COMPUTE_HALO_MAX(NMAXSIZEHALO, TCRRT_COMDATA)
    !
    !NAG4.0 boom avec le 50 lorsqu'on active les scalaires 
    !  NBUFFERSIZE_2D = 50*NMAXSIZEHALO
    NBUFFERSIZE_2D = 150*NMAXSIZEHALO
    !NAG4.0
    NCOMBUFFSIZE2 = NBUFFERSIZE_2D
    !
    DEALLOCATE(TZDZP)
    !
    LINI_PARAZ = .TRUE.
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE INI_PARAZ_ll
  !
  !       ################################
  SUBROUTINE INI_PARAZ_CHILD_ll(KINFO_ll)
    !     ################################
    !
    !!****  *INI_PARAZ_CHILD_ll* - routine to initialize the parallel variables for a child model
    !!                             constructed from a father model in PREP_PGD.
    !!                             Should be called after INI_PARAZ_ll on the father model
    !!                             Similar to INI_PARAZ_ll and INI_CHILD
    !!
    !!    Purpose
    !!    -------
    !     the purpose of the routine is to fill the structured type variables
    !     TCRRT_PROCONF and TCRRT_COMDATA
    !
    !!**  Method
    !!    ------
    !!
    !!    External
    !!    --------
    !     Module MODE_SPLITTING_ll
    !      SPLIT2
    !
    !     Module MODE_CONSTRUCT_ll
    !       INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS,
    !       CONSTRUCT_HALO1, CONSTRUCT_HALO2, CONSTRUCT_HALO_EXTENDED,
    !       CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY,
    !       COMPUTE_HALO_MAX, COMPUTE_TRANS_MAX
    !
    !!    Implicit Arguments
    !!    ------------------
    !     Module MODD_DIM_ll
    !       JPHEXT - Horizontal External points number
    !       NDXRATIO_ALL, NDYRATIO_ALL, NXOR_ALL, NYOR_ALL,
    !       NXEND_ALL, NYEND_ALL,...
    ! 
    !     Module MODD_PARALLEL
    !       TCRRT_PROCONF - Current configuration for current model
    !       TCRRT_COMDATA - Current communication data structure for current model
    !                       and local processor
    ! 
    !     Reference
    !!    ---------
    ! 
    !!    AUTHOR
    !!    ------
    !       M. Moge
    ! 
    !!    MODIFICATIONS
    !!    -------------
    !     Original 21/07/15
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_DIM_ll
    USE MODD_PARAMETERS_ll
    USE MODD_STRUCTURE_ll
    USE MODD_VAR_ll
    !
    USE MODE_SPLITTING_ll, ONLY : SPLIT2
    !
    USE MODE_CONSTRUCT_ll, ONLY : INI_PZ, INI_EZ, INI_BOUNDARIES, INI_TRANS, &
         CONSTRUCT_HALO1, CONSTRUCT_HALO2, CONSTRUCT_HALO_EXTENDED, &
         CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY, &
         COMPUTE_HALO_MAX, COMPUTE_TRANS_MAX
    !
    USE MODE_TOOLSZ_ll, ONLY : SPLITZ, ini_pzz, ini_boundariesz, ini_ezz, construct_transz
    !
    !JUANZ
    USE  MODE_MNH_WORLD , ONLY :  INIT_NMNH_COMM_WORLD
    USE  MODD_CONFZ     , ONLY :  NZ_VERB,NZ_PROC,MPI_BUFFER_SIZE,LMNH_MPI_ALLTOALLV_REMAP,NZ_SPLITTING
    !JUANZ
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(OUT) :: KINFO_ll
    !
    !*       0.2   declarations of local variables
    !
    !INTEGER  ,PARAMETER                      :: MPI_BUFFER_SIZE = 140000000
    CHARACTER,SAVE,ALLOCATABLE,DIMENSION(:)  :: MPI_BUFFER
    !JUAN
    LOGICAL,SAVE                             :: GFIRSTCALL = .TRUE.
    !JUAN

    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP ! intermediate zone
    !
    TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
    TYPE(PROCONF_ll), POINTER :: TZPROCONF
    INTEGER :: JMODEL
    INTEGER     :: IRESP
    LOGICAL     :: GISINIT
    !
    !JUAN
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP1_YP2_Z ! intermediate Full Z = B splitting  without halo zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SX_YP2_ZP1 ! intermediate Full X     splitting zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_Y_ZP1 ! intermediate Full Y     splitting zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_YP1_Z ! intermediate Full Z = B transposed splitting  without halo zone

    INTEGER :: JX_DOMAINS,JY_DOMAINS
    LOGICAL :: LPREM
    INTEGER :: P1,P2
    !JUANZ
    INTEGER :: P1P2(2), P1P2COORD(2) , IROW , ICOL, NROW, NCOL
    LOGICAL :: Lperiodic(2), remain_dims(2) , Lreorder
    INTEGER :: JI
    INTEGER :: IXSIZE_ll	! global sizes of son domain in father grid
    INTEGER :: IYSIZE_ll
    !JUANZ
    !JUAN
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    INITIALIZE MPI :
    !              --------------
    !
    ! MPI should already be initialized
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.    SET OUTPUT FILE :
    !              ---------------

    !  CALL OPEN_ll(UNIT=NIOUNIT,FILE=YOUTPUTFILE,ACTION='write',form&
    !       &='FORMATTED',MODE=SPECIFIC,IOSTAT=IRESP)
    !
    !-------------------------------------------------------------------------------
    !
    !*       3.    ALLOCATION :
    !              ----------
    !
    ! buffer has already been alloacated in the call to INI_PARAZ_ll on the father model

    ALLOCATE(TZDZP(NPROC))
    !JUAN
    ALLOCATE(TZDZP_SXP1_YP2_Z(NPROC))
    ALLOCATE(TZDZP_SXP2_YP1_Z(NPROC))
    ALLOCATE(TZDZP_SX_YP2_ZP1(NPROC))
    ALLOCATE(TZDZP_SXP2_Y_ZP1(NPROC))
    !JUAN
    !
    ALLOCATE(TCRRT_PROCONF)
    CALL ALLOC(TCRRT_COMDATA)
    ALLOCATE(TCRRT_PROCONF%TSPLITS_B(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_X(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_Y(NPROC))
    !JUAN
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1(NPROC))
    ALLOCATE(TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1(NPROC))
    !JUAN
    ALLOCATE(TCRRT_PROCONF%TBOUND(NPROC))
    NULLIFY(TCRRT_PROCONF%TPARENT)
    NULLIFY(TCRRT_COMDATA%TPARENT)
    NULLIFY(TCRRT_PROCONF%TCHILDREN)
    NULLIFY(TCRRT_COMDATA%TCHILDREN)
    !
    !-------------------------------------------------------------------------------
    !
    !*       4.    SPLITTING OF THE DOMAIN :
    !              -----------------------
    !
    IXSIZE_ll = NIMAX_TMP_ll/NDXRATIO_ALL(1)
    IYSIZE_ll = NJMAX_TMP_ll/NDYRATIO_ALL(1)
    DIMX = IXSIZE_ll*NDXRATIO_ALL(1) + 2*JPHEXT
    DIMY = IYSIZE_ll*NDYRATIO_ALL(1) + 2*JPHEXT
    DIMZ = NKMAX_TMP_ll + 2*JPVEXT
    !
    TCRRT_PROCONF%NUMBER = 1
    !

    !JUAN CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SXP2_YP1_Z,'BSPLITTING',NZ_PROC)
!!$    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SXP2_YP1_Z,'BSPLITTING',1)
!!$    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING',NZ_PROC)
!!$    CALL SPLITZ(NIMAX_TMP_ll,NJMAX_TMP_ll,NKMAX_TMP_ll,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING',NZ_PROC)
    ! Add halo directly in Z direction 



    !
    ! find the B spltting
    !
    CALL DEF_SPLITTING2(JX_DOMAINS,JY_DOMAINS,IXSIZE_ll,IYSIZE_ll,NPROC,LPREM)
    !
    P1 = JX_DOMAINS
    IF (DIMZ .NE. 3 )    P1 = MIN(DIMZ,JX_DOMAINS)
    IF (NZ_PROC .GT. 0 ) P1 = NZ_PROC
    P2 = NPROC / P1
    !JUAN PATCH NESTING DIFFERENT SHAPE
    NZ_PROC = P1
    IF (NZ_VERB .GE. 5 ) THEN
       IF ( IP .EQ. 1 )THEN
          print*," INI_PARAZ_ll:: NZ_PROC   =",NZ_PROC
          print*," INI_PARAZ_ll:: JX_DOMAINS=",JX_DOMAINS
          print*," INI_PARAZ_ll:: JY_DOMAINS=",JY_DOMAINS
          print*
          !
          print*," INI_PARAZ_ll:: P1=MIN(NZ_PROC,DIMZ) > 0 .OR. MIN(DIMZ,MAX(JX_DOMAINS,JY_DOMAINS))=",  P1
          !
          print*," INI_PARAZ_ll:: P2=NPROC/P1/                                            =",  P2
       END IF
    END IF
    NP1 = P1
    NP2 = P2
    !
    !JUANZ
    P1P2(1) = NP2
    P1P2(2) = NP1
    Lperiodic(1) = .false.
    Lperiodic(2) = .false.
    Lreorder=.false.
    ! creating cartesian processor grid
    call MPI_Cart_create(NMNH_COMM_WORLD,2,P1P2,Lperiodic,Lreorder,NMNH_P1P2_WORLD,KINFO_ll)
    ! Obtaining process ids with in the cartesian grid
    call MPI_Cart_coords(NMNH_P1P2_WORLD,IP-1,2,P1P2COORD,KINFO_ll)
   
    ! using cart comworld create east-west(row) sub comworld
    remain_dims(1) = .false.
    remain_dims(2) = .true.
    call MPI_Cart_sub(NMNH_P1P2_WORLD,remain_dims,NMNH_ROW_WORLD,KINFO_ll)
    CALL MPI_COMM_RANK(NMNH_ROW_WORLD, IROW, KINFO_ll)
    CALL MPI_COMM_SIZE(NMNH_ROW_WORLD, NROW, KINFO_ll)

    ! using cart comworld create north-south(column) sub comworld
    remain_dims(1) = .true.
    remain_dims(2) = .false.
    call MPI_Cart_sub(NMNH_P1P2_WORLD,remain_dims,NMNH_COL_WORLD,KINFO_ll)
    CALL MPI_COMM_RANK(NMNH_COL_WORLD, ICOL, KINFO_ll)
    CALL MPI_COMM_SIZE(NMNH_COL_WORLD, NCOL, KINFO_ll)
    !JUANZ 

    
    ! split the child model according to the father grid elements (coarse)
    CALL SPLIT2(IXSIZE_ll,IYSIZE_ll,NKMAX_TMP_ll,NPROC,TZDZP,YSPLITTING,P1,P2)
    CALL SPLITZ(IXSIZE_ll,IYSIZE_ll,DIMZ,NPROC,TZDZP_SXP1_YP2_Z,'P1P2SPLITT', 1 ,P1,P2)
    CALL SPLITZ(IXSIZE_ll,IYSIZE_ll,DIMZ,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING', P1,P1,P2)
    CALL SPLITZ(IXSIZE_ll,IYSIZE_ll,DIMZ,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING', P1,P1,P2)
    CALL SPLITZ(IXSIZE_ll,IYSIZE_ll,DIMZ,NPROC,TZDZP_SXP2_YP1_Z,'P2P1SPLITT', 1 ,P1,P2)

    ! 'convert' the splitting from coarse (father) to fine (son) grid using NDXRATIO_ALL(1), NDYRATIO_ALL(1)
    CALL COARSE_TO_FINE(TZDZP)
    CALL COARSE_TO_FINE(TZDZP_SXP1_YP2_Z)
    CALL COARSE_TO_FINE(TZDZP_SX_YP2_ZP1)
    CALL COARSE_TO_FINE(TZDZP_SXP2_Y_ZP1)
    CALL COARSE_TO_FINE(TZDZP_SXP2_YP1_Z)

    !    
    !-------------------------------------------------------------------------------
    !
    !*       5.    INITIALIZATION OF TCRRT_PROCONF :
    !              -------------------------------
    !
    CALL INI_PZ(TCRRT_PROCONF,TZDZP)
    !JUAN
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z,TZDZP_SXP1_YP2_Z)
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z,TZDZP_SXP2_YP1_Z)
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1,TZDZP_SX_YP2_ZP1)
    CALL INI_PZZ(TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1,TZDZP_SXP2_Y_ZP1)
    !JUAN
    !
    CALL INI_BOUNDARIES(TCRRT_PROCONF)
    !JUAN
    CALL INI_BOUNDARIESZ(TCRRT_PROCONF)
    !JUAN
    !
    CALL INI_EZ(TCRRT_PROCONF)
    !JUAN
    CALL INI_EZZ(TCRRT_PROCONF)
    !JUAN
    !
    CALL INI_TRANS(TCRRT_PROCONF)
    !
    !-------------------------------------------------------------------------------
    !
    !*       6.    INITIALIZATION OF TCRRT_COMDATA :
    !              -------------------------------
    !
    !*       6.1    Model Number
    !
    TCRRT_COMDATA%NUMBER = 1
    !
    !*       6.2    Pointer from TCRRT_COMDATA to TCRRT_PROCONF for 2Way splitting
    !
    TCRRT_COMDATA%TSPLIT_B => TCRRT_PROCONF%TSPLITS_B(IP)

    !TZSPLIT => TCRRT_COMDATA%TSPLIT_B
    !
    !
    !*       6.3   Pointer from TCRRT_COMDATA to TCRRT_PROCONF
    !        for x-slices splitting

    TCRRT_COMDATA%TSPLIT_X => TCRRT_PROCONF%TSPLITS_X(IP)
    !
    !TZSPLIT => TCRRT_COMDATA%TSPLIT_X
    !
    !
    !*       6.4   Pointer from TCRRT_COMDATA to TCRRT_PROCONF
    !              for y-slices splitting
    !
    TCRRT_COMDATA%TSPLIT_Y => TCRRT_PROCONF%TSPLITS_Y(IP)
    !
    !TZSPLIT => TCRRT_COMDATA%TSPLIT_Y
    !
    !JUAN
    DO JI=1, NPROC
       IF ( TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z(JI)%NUMBER .EQ. IP ) THEN 
          TCRRT_COMDATA%TSPLIT_SXP1_YP2_Z => TCRRT_PROCONF%TSPLITS_SXP1_YP2_Z(JI)
       ENDIF
       IF ( TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z(JI)%NUMBER .EQ. IP ) THEN 
          TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z => TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z(JI)
       ENDIF
       IF (  TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1(JI)%NUMBER .EQ. IP ) THEN 
          TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1 => TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1(JI)
       ENDIF
       IF ( TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1(JI)%NUMBER .EQ. IP ) THEN
          TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1 => TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1(JI)
       END IF
    END DO
    !JUAN
    !
    !*       6.5   Construction of HALO1 communication data
    !
    CALL CONSTRUCT_HALO1(TCRRT_COMDATA, TCRRT_PROCONF)
    CALL CONSTRUCT_HALO2(TCRRT_COMDATA, TCRRT_PROCONF)
    CALL CONSTRUCT_HALO_EXTENDED(TCRRT_COMDATA, TCRRT_PROCONF, JPHEXT+1)
    !
    !
    !*       6.6   Construction of 1D communication data
    !
    ALLOCATE(TCRRT_COMDATA%HALO1DX)
    ALLOCATE(TCRRT_COMDATA%HALO1DX%NSEND_WEST(NPROC))
    ALLOCATE(TCRRT_COMDATA%HALO1DX%NSEND_EAST(NPROC))
    CALL CONSTRUCT_1DX(TCRRT_COMDATA, TCRRT_PROCONF)
    !
    ALLOCATE(TCRRT_COMDATA%HALO1DY)
    ALLOCATE(TCRRT_COMDATA%HALO1DY%NSEND_SOUTH(NPROC))
    ALLOCATE(TCRRT_COMDATA%HALO1DY%NSEND_NORTH(NPROC))
    CALL CONSTRUCT_1DY(TCRRT_COMDATA, TCRRT_PROCONF)
    !
    !
    !*       6.7   Construction of Transposition communication data
    !
    CALL CONSTRUCT_TRANS(TCRRT_COMDATA, TCRRT_PROCONF)
    CALL CONSTRUCT_TRANSZ(TCRRT_COMDATA, TCRRT_PROCONF)
    !
    !
    !-------------------------------------------------------------------------------
    !
    !        7.    GRID NESTING :
    !              ------------
    !
    ! No grid nesting in this case : We are initializing a child domain directly in PREP_PGD, 
    ! after having called INI_PARAZ_ll on father grid alone
    !
    NULLIFY(TCRRT_PROCONF%TCHILDREN)
    NULLIFY(TCRRT_COMDATA%TCHILDREN)
    NULLIFY(TCRRT_COMDATA%TP2C_DATA)
    !
    !-------------------------------------------------------------------------------
    !
    TZPROCONF => TCRRT_PROCONF
    !
    CALL COMPUTE_TRANS_MAX(NBUFFERSIZE_3D, TCRRT_COMDATA)
    IF (NZ_VERB .GE. 5 ) THEN
       IF (IP.EQ.1) print*,"INI_PARAZ_ll::COMPUTE_TRANS_MAX(NBUFFERSIZE_3D, TCRRT_COMDATA)=",NBUFFERSIZE_3D
    END IF
    !JUAN NCOMBUFFSIZE1 = NBUFFERSIZE_3D
    !NCOMBUFFSIZE1 = NBUFFERSIZE_3D*2
    NCOMBUFFSIZE1 = NBUFFERSIZE_3D
    !JUAN NCOMBUFFSIZE1 = 10000000
    !
    CALL COMPUTE_HALO_MAX(NMAXSIZEHALO, TCRRT_COMDATA)
    !
    !NAG4.0 boom avec le 50 lorsqu'on active les scalaires 
    !  NBUFFERSIZE_2D = 50*NMAXSIZEHALO
    NBUFFERSIZE_2D = 150*NMAXSIZEHALO
    !NAG4.0
    NCOMBUFFSIZE2 = NBUFFERSIZE_2D
    !
    DEALLOCATE(TZDZP)
    !
    LINI_PARAZ = .TRUE.
    !-------------------------------------------------------------------------------
    !
    CONTAINS
      SUBROUTINE COARSE_TO_FINE(TZ)

	IMPLICIT NONE
	
	TYPE(ZONE_ll), DIMENSION(:) :: TZ   ! grid splitting to transform from coarse (father) resolution/grid
					    ! to fien ( son ) resolution/grid    

	INTEGER :: J

	DO J = 1, NPROC
	  !
	  TZ(J)%NUMBER = TZ(J)%NUMBER
	  TZ(J)%NXOR  = (TZ(J)%NXOR - JPHEXT -1 ) * NDXRATIO_ALL(1) + JPHEXT +1 
	  TZ(J)%NYOR  = (TZ(J)%NYOR - JPHEXT -1 ) * NDYRATIO_ALL(1) + JPHEXT +1 
	  TZ(J)%NXEND = (TZ(J)%NXEND - JPHEXT   ) * NDXRATIO_ALL(1) + JPHEXT
	  TZ(J)%NYEND = (TZ(J)%NYEND - JPHEXT   ) * NDYRATIO_ALL(1) + JPHEXT
	  !JUAN Z_SPLITTING
	  TZ(J)%NZOR  = TZ(J)%NZOR
	  TZ(J)%NZEND = TZ(J)%NZEND
	  !JUAN Z_SPLITTING
	  !
	ENDDO

      END SUBROUTINE COARSE_TO_FINE
  
  END SUBROUTINE INI_PARAZ_CHILD_ll
  !
  !     #######################################
!!$  SUBROUTINE SET_NZ_PROC_ll(KZ_PROC)
!!$    !     #######################################
!!$    !
!!$    !!****  *SET_NZ_PROC_ll* - 
!!$    !
!!$    !!    Purpose
!!$    !!    -------
!!$    !       Set the number of Proc in Z splitting
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$    USE MODD_VAR_ll, ONLY : NZ_PROC
!!$    !
!!$    IMPLICIT NONE
!!$    !
!!$    INTEGER :: KZ_PROC
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$    NZ_PROC  = KZ_PROC
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$  END SUBROUTINE SET_NZ_PROC_ll

  !     #################################################
  SUBROUTINE GET_DIM_EXTZ_ll( HSPLIT, KXDIM, KYDIM, KZDIM )
    !     #################################################
    !
    !!****  *GET_DIM_EXTZ_ll* - returns the dimensions of the extended z-2way subdomain
    !                   or of the zx-slices subdomain or of the zy-slices
    !                   subdomain of the local processor
    ! 
    !!    Purpose
    !!    -------
    !     the Purpose of this routine is to give subdomain dimension
    !
    !!**  Method
    !!    ------
    !     if HSPLIT='SXP1_YP2_Z', the dimensions of the extended 2way subdomain are returned
    !     if HSPLIT='SXP2_YP1_Z', the dimensions of the extended 2way subdomain are returned
    !     if HSPLIT='SX_YP2_ZP1', the dimensions of x-slices subdomain are returned
    !     if HSPLIT='SXP2_Y_ZP1', the dimensions of y-slices subdomain are returned
    ! 
    !!    External
    !!    --------
    ! 
    !!    Implicit Arguments
    !!    ------------------
    !     Module MODD_VAR_ll
    !        TCRRT_COMDATA - Current communication data structure for current model
    !                        and local processor
    ! 
    !!    Reference
    !!    ---------
    ! 
    !!    Author
    !!    ------
    !     R. Guivarch
    ! 
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    ! 
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    CHARACTER(len=*), INTENT(IN) :: HSPLIT
    !
    INTEGER, INTENT(OUT) :: KXDIM, KYDIM, KZDIM
    !
    !
    !*       0.2   declarations of local variables
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    Return the dimensions
    !
    IF( HSPLIT .EQ. 'SXP2_YP1_Z' ) THEN
       KXDIM = TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z%NDIMXE
       KYDIM = TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z%NDIMYE
       KZDIM = TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z%NDIMZE
    ELSEIF ( HSPLIT .EQ. 'SXP1_YP2_Z' ) THEN
       KXDIM = TCRRT_COMDATA%TSPLIT_SXP1_YP2_Z%NDIMXE
       KYDIM = TCRRT_COMDATA%TSPLIT_SXP1_YP2_Z%NDIMYE
       KZDIM = TCRRT_COMDATA%TSPLIT_SXP1_YP2_Z%NDIMZE
    ELSEIF ( HSPLIT .EQ. 'SX_YP2_ZP1' ) THEN
       KXDIM = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NDIMXE
       KYDIM = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NDIMYE
       KZDIM = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NDIMZE
    ELSE
       KXDIM = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NDIMXE
       KYDIM = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NDIMYE
       KZDIM = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NDIMZE
    ENDIF
!!$    ELSEIF ( HSPLIT .EQ. 'SX_YP2_ZP1' ) THEN
!!$       KXDIM = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NDIMXP
!!$       KYDIM = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NDIMYP
!!$       KZDIM = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NDIMZP
!!$    ELSE
!!$       KXDIM = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NDIMXP
!!$       KYDIM = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NDIMYP
!!$       KZDIM = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NDIMZP
!!$    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE GET_DIM_EXTZ_ll
  !
  !
  !     ######################################################
  SUBROUTINE REMAP_B_SX_YP2_ZP1_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_B_SX_YP2_ZP1_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in B splitting to PFIELDOUT in SX_YP2_ZP1 splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER :: KINFO ! return status
    !
    !*       0.2   declarations of local var
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_B_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_B_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_B_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_B_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    !
    !-------------------------------------------------------------------------------
    !
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TSEND_BOX_B_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_BOX_B_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
  END SUBROUTINE REMAP_B_SX_YP2_ZP1_ll
  !
  !     ######################################################
  SUBROUTINE REMAP_SXP1_YP2_Z_SX_YP2_ZP1_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SXP1_YP2_Z_SX_YP2_ZP1_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in B splitting to PFIELDOUT in SX_YP2_ZP1 splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER                               :: KINFO ! return status
    !
    !*       0.2   declarations of local var
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_SXP1_YP2_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_SXP1_YP2_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_SXP1_YP2_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_SXP1_YP2_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    !
    !-------------------------------------------------------------------------------
    !
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TSEND_BOX_SXP1_YP2_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_BOX_SXP1_YP2_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
  END SUBROUTINE REMAP_SXP1_YP2_Z_SX_YP2_ZP1_ll

  !     ######################################################
  SUBROUTINE REMAP_SX_YP2_ZP1_B_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SX_YP2_ZP1_B_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in B splitting to PFIELDOUT in SX_YP2_ZP1 splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER :: KINFO ! return status
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_B_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_B_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_B_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_B_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TRECV_BOX_B_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_BOX_B_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SX_YP2_ZP1_B_ll
    !
  !     ######################################################
  SUBROUTINE REMAP_SX_YP2_ZP1_SXP1_YP2_Z_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SX_YP2_ZP1_SXP1_YP2_Z_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in B splitting to PFIELDOUT in SX_YP2_ZP1 splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER :: KINFO ! return status
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_SXP1_YP2_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_SXP1_YP2_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_SXP1_YP2_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_SXP1_YP2_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TRECV_BOX_SXP1_YP2_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_BOX_SXP1_YP2_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SX_YP2_ZP1_SXP1_YP2_Z_ll
  !
  !     ######################################################
  SUBROUTINE REMAP_SX_YP2_ZP1_SXP2_Y_ZP1_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SX_YP2_ZP1_SXP2_Y_ZP1_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in SX_YP2_ZP1splitting to PFIELDOUT in SXP2_Y_ZP1splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER                               :: KINFO ! return status
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_SX_YP2_ZP1_SXP2_Y_ZP1, &
         TCRRT_COMDATA%TRECV_SX_YP2_ZP1_SXP2_Y_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_SX_YP2_ZP1_SXP2_Y_ZP1, &
         TCRRT_COMDATA%TRECV_SX_YP2_ZP1_SXP2_Y_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    !
    !-------------------------------------------------------------------------------
    !
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TSEND_BOX_SX_YP2_ZP1_SXP2_Y_ZP1, &
         TCRRT_COMDATA%TRECV_BOX_SX_YP2_ZP1_SXP2_Y_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SX_YP2_ZP1_SXP2_Y_ZP1_ll

  !     ######################################################
  SUBROUTINE REMAP_SXP2_Y_ZP1_SX_YP2_ZP1_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SXP2_Y_ZP1_SX_YP2_ZP1_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in SXP2_Y_ZP1 splitting to PFIELDOUT in SX_YP2_ZP1 splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER :: KINFO ! return status
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_SX_YP2_ZP1_SXP2_Y_ZP1, &
         TCRRT_COMDATA%TSEND_SX_YP2_ZP1_SXP2_Y_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_SX_YP2_ZP1_SXP2_Y_ZP1, &
         TCRRT_COMDATA%TSEND_SX_YP2_ZP1_SXP2_Y_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    !
    !-------------------------------------------------------------------------------
    !
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TRECV_BOX_SX_YP2_ZP1_SXP2_Y_ZP1, &
         TCRRT_COMDATA%TSEND_BOX_SX_YP2_ZP1_SXP2_Y_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SXP2_Y_ZP1_SX_YP2_ZP1_ll

  !     ######################################################
  SUBROUTINE REMAP_SXP2_Y_ZP1_SXP2_YP1_Z_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SXP2_Y_ZP1_SXP2_YP1_Z_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in SXP2_Y_ZP1 splitting to PFIELDOUT in SXP2_YP1_Z splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER :: KINFO ! return status
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_SXP2_YP1_Z, &
         TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_SXP2_YP1_Z, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_SXP2_YP1_Z, &
         TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_SXP2_YP1_Z, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TSEND_BOX_SXP2_Y_ZP1_SXP2_YP1_Z, &
         TCRRT_COMDATA%TRECV_BOX_SXP2_Y_ZP1_SXP2_YP1_Z, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SXP2_Y_ZP1_SXP2_YP1_Z_ll

  !     ######################################################
!!$  SUBROUTINE REMAP_SXP2_Y_ZP1_B_ll(PFIELDIN, PFIELDOUT, KINFO)
!!$    !     ######################################################
!!$    !
!!$    !!**** *REMAP_SXP2_Y_ZP1_B_ll* -
!!$    !!
!!$    !!    Purpose
!!$    !!    -------
!!$    !       remaps PFIELDIN from in SXP2_Y_ZP1 splitting to PFIELDOUT in B splitting
!!$    !
!!$    ! 
!!$    !!    Author
!!$    !!    ------
!!$    !     J. ESCOBAR             * LA - CNRS *
!!$    !!
!!$    !!    Modifications
!!$    !!    -------------
!!$    !     Original 24/01/2008
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$    !*       0.    DECLARATIONS
!!$    !
!!$    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
!!$    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
!!$    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
!!$    !
!!$    IMPLICIT NONE
!!$    !
!!$    !*       0.1   declarations of arguments
!!$    !
!!$    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
!!$    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
!!$    INTEGER :: KINFO ! return status
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
!!$    !               -------------------------------------------------------------
!!$    !
!!$    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
!!$    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_B, &
!!$         TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_B, &
!!$         PFIELDIN, PFIELDOUT, KINFO)
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
!!$    !               ------------------------------------------------------------
!!$    !
!!$    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_B, &
!!$         TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_B, &
!!$         PFIELDIN, PFIELDOUT, KINFO)
!!$    ELSE
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$    CALL ALL_SEND_RECV(TCRRT_COMDATA%TSEND_BOX_SXP2_Y_ZP1_B, &
!!$         TCRRT_COMDATA%TRECV_BOX_SXP2_Y_ZP1_B, &
!!$         PFIELDIN, PFIELDOUT, KINFO)
!!$    ENDIF
!!$    !
!!$    !-------------------------------------------------------------------------------
!!$    !
!!$  END SUBROUTINE REMAP_SXP2_Y_ZP1_B_ll

  !     ######################################################
  SUBROUTINE REMAP_SXP2_YP1_Z_SXP2_Y_ZP1_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SXP2_YP1_Z_SXP2_Y_ZP1_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in B splitting to PFIELDOUT in SXP2_Y_ZP1 splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER :: KINFO ! return status
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_SXP2_YP1_Z, &
         TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_SXP2_YP1_Z, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_SXP2_YP1_Z, &
         TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_SXP2_YP1_Z, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TRECV_BOX_SXP2_Y_ZP1_SXP2_YP1_Z, &
         TCRRT_COMDATA%TSEND_BOX_SXP2_Y_ZP1_SXP2_YP1_Z, &
         PFIELDIN, PFIELDOUT, KINFO)
    END IF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SXP2_YP1_Z_SXP2_Y_ZP1_ll

  !     ######################################################
  SUBROUTINE REMAP_B_SXP2_Y_ZP1_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_B_SXP2_Y_ZP1_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in B splitting to PFIELDOUT in SXP2_Y_ZP1 splitting
    !
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR             * LA - CNRS *
    !!
    !!    Modifications
    !!    -------------
    !     Original 24/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODE_EXCHANGE_ll, ONLY : SEND_RECV_FIELD, COPY_CRSPD_TRANS
    USE MODD_VAR_ll     , ONLY : TCRRT_COMDATA
    USE MODD_CONFZ      , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
    INTEGER :: KINFO ! return status
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.     UPDATE THE ZONES NOT SENT OR RECEIVED BY THE PROCESSOR ITSELF 
    !               -------------------------------------------------------------
    !
    IF ( IAND(NZ_SPLITTING,4) .EQ. 0 ) THEN 
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_B, &
         TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_B, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_B, &
         TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_B, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    !
    !-------------------------------------------------------------------------------
    !
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TRECV_BOX_SXP2_Y_ZP1_B, &
         TCRRT_COMDATA%TSEND_BOX_SXP2_Y_ZP1_B, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_B_SXP2_Y_ZP1_ll



!     #######################################
      LOGICAL FUNCTION LWESTZ_ll(K,HSPLITTING)
!     #######################################
!
!!****  *LWEST_ll* - function which returns the position on to the boundaries
!                  of the subdomain K according to the splitting
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to offer a transparent way to obtain
!     the position of a subdomain
!
!!**  Method
!!    ------
!     if the argument HSPLITTING is omitted the 2Way splitting is considered
!     if the argument K is omitted, the local subdomain is considered
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_PROCONF - Current configuration for current model
!        IP - Number of local processor=subdomain
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN), OPTIONAL :: K ! number of the subdomain
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting

!!
!*       0.2   declarations of local variables
!
  INTEGER :: IT ! number of the tested subdomain
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF THE RESULT :
!              -------------------------
!
  IF( PRESENT(K) ) THEN
     IT = K
  ELSE
     IT = IP
  ENDIF
!
  LWESTZ_ll = .FALSE.
!
  IF(.NOT.PRESENT(HSPLITTING)) THEN
    LWESTZ_ll = TCRRT_PROCONF%TBOUND_SXP2_YP1_Z(IT)%WEST
  ELSEIF(HSPLITTING .EQ. 'SXP2_YP1_Z') THEN
    LWESTZ_ll = TCRRT_PROCONF%TBOUND_SXP2_YP1_Z(IT)%WEST
  ELSEIF(HSPLITTING .EQ. 'SXP1_YP2_Z') THEN
    LWESTZ_ll = TCRRT_PROCONF%TBOUND_SXP1_YP2_Z(IT)%WEST
  ELSEIF(HSPLITTING .EQ. 'SX_YP2_ZP1') THEN
    LWESTZ_ll = TCRRT_PROCONF%TBOUND_SX_YP2_ZP1(IT)%WEST
  ELSEIF(HSPLITTING .EQ. 'SXP2_Y_ZP1') THEN
    LWESTZ_ll = TCRRT_PROCONF%TBOUND_SXP2_Y_ZP1(IT)%WEST
  ENDIF
!
!-------------------------------------------------------------------------------
!
      END FUNCTION LWESTZ_ll

!
!     ########################################
      LOGICAL FUNCTION LSOUTHZ_ll(K,HSPLITTING)
!     ########################################
!
!!****  *LSOUTH_ll* - function which returns the position on to the boundaries
!!                 of the subdomain K according to the splitting
!!
!!    Purpose
!!    -------
!     the Purpose of this routine is to offer a transparent way to obtain
!     the position of a subdomain
!
!!**  Method
!!    ------
!     if the argument HSPLITTING is omitted the 2Way splitting is considered
!     if the argument K is omitted, the local subdomain is considered
! 
!!    External
!!    --------
! 
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_PROCONF - Current configuration for current model
!        IP - Number of local processor=subdomain
! 
!!    Reference
!!    ---------
! 
!!    Author
!!    ------
!     R. Guivarch
! 
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_PROCONF, IP
!
  IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
  INTEGER, INTENT(IN), OPTIONAL :: K ! number of the subdomain
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: HSPLITTING ! kind of splitting
 
!!
!*       0.2   declarations of local variables
!
  INTEGER :: IT ! number of the tested subdomain
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF THE RESULT :
!              -------------------------
!
  IF( PRESENT(K) ) THEN
     IT = K
  ELSE
     IT = IP
  ENDIF
!
  LSOUTHZ_ll = .FALSE.
!
  IF(.NOT.PRESENT(HSPLITTING)) THEN
    LSOUTHZ_ll = TCRRT_PROCONF%TBOUND_SXP2_YP1_Z(IT)%SOUTH
  ELSEIF(HSPLITTING .EQ. 'SXP2_YP1_Z') THEN
    LSOUTHZ_ll = TCRRT_PROCONF%TBOUND_SXP2_YP1_Z(IT)%SOUTH
  ELSEIF(HSPLITTING .EQ. 'SXP1_YP2_Z') THEN
    LSOUTHZ_ll = TCRRT_PROCONF%TBOUND_SXP1_YP2_Z(IT)%SOUTH
  ELSEIF(HSPLITTING .EQ. 'SX_YP2_ZP1') THEN
    LSOUTHZ_ll = TCRRT_PROCONF%TBOUND_SX_YP2_ZP1(IT)%SOUTH
  ELSEIF(HSPLITTING .EQ. 'SXP2_Y_ZP1') THEN
    LSOUTHZ_ll = TCRRT_PROCONF%TBOUND_SXP2_Y_ZP1(IT)%SOUTH
  ENDIF
!
!-------------------------------------------------------------------------------
!
END FUNCTION LSOUTHZ_ll

    SUBROUTINE ALL_SEND_RECV(TSEND_BOX_FROM,TRECV_BOX_TO, &
         PFIELDIN, PFIELDOUT, KINFO)
      !
      USE MODD_STRUCTURE_ll , ONLY : BOX_ll
      USE MODD_VAR_ll       , ONLY : MPI_PRECISION
      !JUANZ
      !USE MODD_MPIF         , ONLY : MPI_COMM_WORLD
      USE MODD_VAR_ll        , ONLY : NMNH_COMM_WORLD
      !JUANZ
      IMPLICIT NONE
      !
      ! Argument
      !
      TYPE(BOX_ll)                , POINTER :: TSEND_BOX_FROM,TRECV_BOX_TO
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDIN ! field to be sent
      REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT ! reception field
      INTEGER                               :: KINFO ! return status
      !
      ! local var
      !
      INTEGER                               :: JB,JI,JJ,JK, JCNT ! loop
      REAL, DIMENSION(:), ALLOCATABLE       :: ZSEND ! buffer to be sent
      REAL, DIMENSION(:), ALLOCATABLE       :: ZRECV ! buffer to be recv      
      !
      ALLOCATE(ZSEND(TSEND_BOX_FROM%NSIZE))
      ALLOCATE(ZRECV(TRECV_BOX_TO%NSIZE))
      !
      JCNT = 0
      DO JB = 1, TSEND_BOX_FROM%NBOX
         IF  ( TSEND_BOX_FROM%NCNT(JB) .NE. 0 ) THEN
            DO JK=TSEND_BOX_FROM%NZOR(JB),TSEND_BOX_FROM%NZEND(JB)
               DO JJ=TSEND_BOX_FROM%NYOR(JB),TSEND_BOX_FROM%NYEND(JB)
                  DO JI=TSEND_BOX_FROM%NXOR(JB),TSEND_BOX_FROM%NXEND(JB)
                     JCNT =  JCNT + 1
                     ZSEND(JCNT) = PFIELDIN(JI,JJ,JK)
                  END DO
               END DO
            END DO
         END IF
      END DO
      !
      CALL mpi_alltoallv(ZSEND,TSEND_BOX_FROM%NCNT,TSEND_BOX_FROM%NSTRT,MPI_PRECISION,&
                         ZRECV,TRECV_BOX_TO%NCNT  ,TRECV_BOX_TO%NSTRT  ,MPI_PRECISION,&
                         TSEND_BOX_FROM%NCOM,KINFO)
      !
      JCNT = 0
      PFIELDOUT = 0.0
      DO JB = 1, TRECV_BOX_TO%NBOX
         IF  ( TRECV_BOX_TO%NCNT(JB) .NE. 0 ) THEN
            DO JK=TRECV_BOX_TO%NZOR(JB),TRECV_BOX_TO%NZEND(JB)
               DO JJ=TRECV_BOX_TO%NYOR(JB),TRECV_BOX_TO%NYEND(JB)
                  DO JI=TRECV_BOX_TO%NXOR(JB),TRECV_BOX_TO%NXEND(JB)
                     JCNT =  JCNT + 1
                     PFIELDOUT(JI,JJ,JK) =  ZRECV(JCNT)
                  END DO
               END DO
            END DO
         END IF
      END DO
      !
      DEALLOCATE(ZSEND,ZRECV)
    END SUBROUTINE ALL_SEND_RECV

!     ##########################################
      SUBROUTINE GET_ORZ_ll( HSPLIT, KXOR, KYOR )
!     ##########################################
!
!!****  *GET_ORZ_ll* - returns the origin'coordinates of the extended
!                     2way subdomain or of the x-slices subdomain
!                     or of the y-slices
!                     subdomain of the local processor (global indices)
! 
!!    Purpose
!!    -------
!     the Purpose of this routine is to give subdomain origin
!
!!**  Method
!!    ------
!     if HSPLIT = 'SXP2_YP1_Z', the origin of the extended subdomain are returned
!     if HSPLIT = 'SX_YP2_ZP1', the origin of x-slices subdomain are returned
!     if HSPLIT = 'SXP2_Y_ZP1', the origin of y-slices subdomain are returned
! 
!!    External
!!    --------
!     
!!    Implicit Arguments
!!    ------------------
!     Module MODD_VAR_ll
!        TCRRT_COMDATA - Current communication data structure for current model
!                        and local processor
!     
!!    Reference
!!    ---------
!     
!!    Author
!!    ------
!     R. Guivarch
!     
!!    Modifications
!!    -------------
!     Original 01/05/98
! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
  USE MODD_VAR_ll, ONLY : TCRRT_COMDATA
!
  IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
  CHARACTER(len=*), INTENT(IN) :: HSPLIT
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!
!*       1.    
!
  IF( HSPLIT .EQ. 'SXP2_YP1_Z' ) THEN
    KXOR = TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z%NXORE
    KYOR = TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z%NYORE
  ELSEIF ( HSPLIT .EQ. 'SX_YP2_ZP1' ) THEN
    KXOR = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NXORP
    KYOR = TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1%NYORP
  ELSE
    KXOR = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NXORP
    KYOR = TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1%NYORP
  ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_ORZ_ll

END MODULE MODE_SPLITTINGZ_ll
