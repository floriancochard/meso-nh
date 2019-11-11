!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

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
  USE MODD_MPIF
  !
  USE MODE_SPLITTING_ll
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
    !       CONSTRUCT_HALO1, CONSTRUCT_HALO2,
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
         CONSTRUCT_HALO1, CONSTRUCT_HALO2, &
         CONSTRUCT_TRANS, CONSTRUCT_1DX, CONSTRUCT_1DY, &
         COMPUTE_HALO_MAX, COMPUTE_TRANS_MAX
    !
    USE MODE_NEST_ll, ONLY : INI_CHILD
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(OUT) :: KINFO_ll
    !
    !*       0.2   declarations of local variables
    !
#if defined(MNH_MPI_BSEND)
    INTEGER  ,PARAMETER                      :: MPI_BUFFER_SIZE = 40000000
    CHARACTER,SAVE,ALLOCATABLE,DIMENSION(:)  :: MPI_BUFFER
    !JUAN
    LOGICAL,SAVE                             :: GFIRSTCALL = .TRUE.
    !JUAN
#endif

    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP ! intermediate zone
    !
    TYPE(MODELSPLITTING_ll), POINTER :: TZSPLIT
    TYPE(PROCONF_ll), POINTER :: TZPROCONF
    INTEGER :: JMODEL
    INTEGER     :: IRESP
    LOGICAL     :: GISINIT
    !
    !JUAN
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_YP1_Z ! intermediate Z_2WAY splitting zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SX_YP2_ZP1 ! intermediate Z_X    splitting zone
    TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TZDZP_SXP2_Y_ZP1 ! intermediate Z_Y    splitting zone

    INTEGER :: JX_DOMAINS,JY_DOMAINS
    LOGICAL :: LPREM
    INTEGER :: P1,P2
    !JUAN
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    INITIALIZE MPI :
    !              --------------
    !
    KINFO_ll = 0
    CALL MPI_INITIALIZED(GISINIT, KINFO_ll)
    IF (.NOT. GISINIT) THEN
       CALL MPI_INIT(KINFO_ll)
    END IF
    !
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, IP, KINFO_ll)
    !
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, KINFO_ll)
    !
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, NHALO_COM, KINFO_ll)
    !
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, NHALO2_COM, KINFO_ll)
    !
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, NTRANS_COM, KINFO_ll)
    !
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, NGRID_COM, KINFO_ll)
    !
    IP = IP + 1
    !
    MPI_PRECISION = MPI_DOUBLE_PRECISION
    MPI_2PRECISION = MPI_2DOUBLE_PRECISION
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
#if defined(MNH_MPI_BSEND)
    IF (GFIRSTCALL) THEN
       ALLOCATE(MPI_BUFFER(MPI_BUFFER_SIZE))
       CALL MPI_BUFFER_ATTACH(MPI_BUFFER,MPI_BUFFER_SIZE,KINFO_ll)
       GFIRSTCALL = .FALSE.
    ENDIF
#endif

    ALLOCATE(TZDZP(NPROC))
    !JUAN
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
    DIMX = NIMAX_ll + 2*JPHEXT
    DIMY = NJMAX_ll + 2*JPHEXT
    DIMZ = NKMAX_ll + 2*JPVEXT
    !
    TCRRT_PROCONF%NUMBER = 1
    !

    !JUAN CALL SPLITZ(NIMAX_ll,NJMAX_ll,NKMAX_ll,NPROC,TZDZP_SXP2_YP1_Z,'BSPLITTING',NZ_PROC)
!!$    CALL SPLITZ(NIMAX_ll,NJMAX_ll,NKMAX_ll,NPROC,TZDZP_SXP2_YP1_Z,'BSPLITTING',1)
!!$    CALL SPLITZ(NIMAX_ll,NJMAX_ll,NKMAX_ll,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING',NZ_PROC)
!!$    CALL SPLITZ(NIMAX_ll,NJMAX_ll,NKMAX_ll,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING',NZ_PROC)
    ! Add halo directly in Z direction 

    CALL SPLIT2(NIMAX_ll,NJMAX_ll,NKMAX_ll,NPROC,TZDZP,YSPLITTING)

    !
    ! find the B spltting
    !
    CALL DEF_SPLITTING2(JX_DOMAINS,JY_DOMAINS,NIMAX_ll,NJMAX_ll,NPROC,LPREM)
    !
    P1 = MIN(DIMZ,JX_DOMAINS)
    P2 = NPROC / P1
    IF ( IP .EQ. 1 )THEN
       print*," INI_PARAZ_ll:: NZ_PROC   =",NZ_PROC
       print*," INI_PARAZ_ll:: JX_DOMAINS=",JX_DOMAINS
       print*," INI_PARAZ_ll:: JY_DOMAINS=",JY_DOMAINS
       print*
       !
       print*," INI_PARAZ_ll:: P1=MIN(DIMZ,MAX(JX_DOMAINS,JY_DOMAINS))=",  P1
       !
       print*," INI_PARAZ_ll:: P2=NPROC/P1/                           =",  P2
    END IF
    !

    CALL SPLITZ(NIMAX_ll,NJMAX_ll,DIMZ,NPROC,TZDZP_SX_YP2_ZP1,'YSPLITTING',P1)
    CALL SPLITZ(NIMAX_ll,NJMAX_ll,DIMZ,NPROC,TZDZP_SXP2_Y_ZP1,'XSPLITTING',P1)
    CALL SPLITZ(NIMAX_ll,NJMAX_ll,DIMZ,NPROC,TZDZP_SXP2_YP1_Z,'BSPLITTING',1)


!JUAN    CALL SPLIT2(NIMAX_ll,NJMAX_ll,NKMAX_ll,NPROC,TZDZP,YSPLITTING)


    !    
    !-------------------------------------------------------------------------------
    !
    !*       5.    INITIALIZATION OF TCRRT_PROCONF :
    !              -------------------------------
    !
    CALL INI_PZ(TCRRT_PROCONF,TZDZP)
    !JUAN
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
    TCRRT_COMDATA%TSPLIT_SXP2_YP1_Z => TCRRT_PROCONF%TSPLITS_SXP2_YP1_Z(IP)
    TCRRT_COMDATA%TSPLIT_SX_YP2_ZP1 => TCRRT_PROCONF%TSPLITS_SX_YP2_ZP1(IP)
    TCRRT_COMDATA%TSPLIT_SXP2_Y_ZP1 => TCRRT_PROCONF%TSPLITS_SXP2_Y_ZP1(IP)
    !JUAN
    !
    !*       6.5   Construction of HALO1 communication data
    !
    CALL CONSTRUCT_HALO1(TCRRT_COMDATA, TCRRT_PROCONF)
    CALL CONSTRUCT_HALO2(TCRRT_COMDATA, TCRRT_PROCONF)
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
    IF (IP.EQ.1) print*,"INI_PARAZ_ll::COMPUTE_TRANS_MAX(NBUFFERSIZE_3D, TCRRT_COMDATA)=",NBUFFERSIZE_3D
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
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE INI_PARAZ_ll
  !
  !     #################################################################
  SUBROUTINE SPLITZ(X_DIM,Y_DIM,Z_DIM,NB_PROC,TPROC,HSPLITTING,KZ_PROC)
    !     #################################################################
    !
    !!****  *SPLITZ* - routine which splits a domain in NB_PROC sub-domains
    !                     by using DEFINE_SPLITTING2.
    !
    !!    Purpose
    !!    -------
    !     this routine fills the fields of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !     Module MODD_VAR_ll
    !       JPHEXT
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !     R.Guivarch 29/11/99 : x and y splitting : HSPLITTING
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
    !JUAN
    USE MODD_VAR_ll, ONLY : IP
    USE MODD_CONF  , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
    !JUAN
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(IN) :: NB_PROC,X_DIM,Y_DIM,Z_DIM
    CHARACTER*10, INTENT(IN) :: HSPLITTING ! kind of splitting
    TYPE(ZONE_LL), INTENT(OUT), DIMENSION(NB_PROC),TARGET  :: TPROC
    !
    !JUAN
    INTEGER, INTENT(IN) :: KZ_PROC ! number of proc for Z splitting

    !JUAN
    !
    !*       0.2   declarations of local variables
    !
    INTEGER                 :: X_DOMAINS,Y_DOMAINS,Z_DOMAINS,X_DOMAINS_NEW
    LOGICAL                 :: PREM
    INTEGER                 :: IK
    INTEGER                 :: IDOM 
    TYPE(ZONE_LL), POINTER  :: TP
    INTEGER                 :: NB_PROC_XY
    !
    !        0. CHECK NB_PROC/NZ_PROC
    IF ( MOD(NB_PROC,KZ_PROC) .NE. 0 ) THEN
       PRINT*
       WRITE(*,1000) NB_PROC, KZ_PROC
       PRINT*
1000   FORMAT("MODE_SPLITTINGZ::SPLITZ --> NB_PROC=", I4 ," NOT DIVISIBLE BY KZ_PROC=", I4)
       STOP
    ENDIF
    !
    !   Splitting in Z possible so
    !
    NB_PROC_XY = NB_PROC / KZ_PROC
    Z_DOMAINS  = KZ_PROC
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    FIND THE SPLITTING in XY & Z 
    !
    !  CALL DEF_SPLITTINGZ(X_DOMAINS,Y_DOMAINS,Z_DOMAINS,X_DIM,Y_DIM,Z_DIM,NB_PROC,KZ_PROC,PREM)
    CALL DEF_SPLITTING2(X_DOMAINS,Y_DOMAINS,X_DIM,Y_DIM,NB_PROC_XY,PREM)
    !
    ! transpose the data X & Y
    X_DOMAINS_NEW = Y_DOMAINS
    Y_DOMAINS     = X_DOMAINS
    X_DOMAINS     = X_DOMAINS_NEW
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    FILL THE FIELDS OF TPROC
    !
    IF(HSPLITTING.EQ."BSPLITTING") THEN
       IF ((PREM).AND.(NB_PROC_XY.GT.2)) THEN
          !
          !   split x direction only on NB_PROC_XY - 1 processors
          !   and on reducted x-size = X_DIM - X_DIM/NB_PROC_XY -1
          !
          CALL DEF_SPLITTING2(X_DOMAINS,Y_DOMAINS,X_DIM - X_DIM/NB_PROC_XY -1,  &
               Y_DIM,NB_PROC_XY-1,PREM)
          !
          !     the last Z processor slide hold last all Y dim 
          !
          IDOM  = NB_PROC - KZ_PROC
          DO IK = 1,KZ_PROC
             IDOM = IDOM +  1
             TP => TPROC(IDOM)
             TP%NUMBER = IDOM
             ! X coordonate
             TP%NXOR  = X_DIM - X_DIM/NB_PROC_XY
             TP%NXEND = X_DIM
             ! Y coordonate
             TP%NYOR  = 1
             TP%NYEND = Y_DIM
             ! Z coordonate
             CALL SLIDE_COORD(Z_DIM,KZ_PROC,IK,TP%NZOR,TP%NZEND)
          ENDDO
          !
          ! cartesian splitting with NB_PROC-1
          !
          CALL CARTESIANZ(TPROC,NB_PROC_XY-1,X_DIM-X_DIM/NB_PROC_XY-1,   &
               Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS)
          !
       ELSE ! XY & Z splitting --> general case
          !
          CALL CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS)

       END IF
    ELSEIF(HSPLITTING.EQ."XSPLITTING") THEN
       !
       X_DOMAINS=NB_PROC_XY
       Y_DOMAINS=1
       CALL CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS)
       !
    ELSE
       !
       X_DOMAINS=1
       Y_DOMAINS=NB_PROC_XY
       CALL CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS)
       !
    END IF
    !
    !*       3.     shift from physical to extended domain
    !
    IF ( IAND(NZ_SPLITTING,2) > 0 ) THEN 
       IF ( IP .EQ. 1 )THEN
          PRINT*,"********************************************************************"
          PRINT*,"******************** HSPLITTING=",HSPLITTING," *************************"
          PRINT*,"********************************************************************"
          PRINT*,"NB_PROC=",NB_PROC," KZ_PROC=",KZ_PROC
          PRINT*,"X_DIM=", X_DIM ,", X_DOMAINS=",X_DOMAINS
          PRINT*,"Y_DIM=", Y_DIM ,", Y_DOMAINS=",Y_DOMAINS
          PRINT*,"Z_DIM=", Z_DIM ,", Z_DOMAINS=",Z_DOMAINS
          
          DO IK = 1,NB_PROC,MAX(1,NB_PROC-1)
             PRINT*," ============== NPROC=",IK,"========================"
             PRINT*,"NXOR=",TPROC(IK)%NXOR," NXEND=",TPROC(IK)%NXEND," TAILLE=",1+TPROC(IK)%NXEND-TPROC(IK)%NXOR
             PRINT*,"NYOR=",TPROC(IK)%NYOR," NYEND=",TPROC(IK)%NYEND," TAILLE=",1+TPROC(IK)%NYEND-TPROC(IK)%NYOR
             PRINT*,"NZOR=",TPROC(IK)%NZOR," NZEND=",TPROC(IK)%NZEND," TAILLE=",1+TPROC(IK)%NZEND-TPROC(IK)%NZOR
          END DO
       ENDIF
    END IF
    !    STOP
    !
    ! Add 'Halo points' to global coordonne in X & Y direction
    !
    TPROC(:)%NXOR  = TPROC(:)%NXOR  + JPHEXT
    TPROC(:)%NXEND = TPROC(:)%NXEND + JPHEXT
    !
    TPROC(:)%NYOR  = TPROC(:)%NYOR  + JPHEXT
    TPROC(:)%NYEND = TPROC(:)%NYEND + JPHEXT
    !
    ! In Z direction Halo already intégrated in Z dimension 
    !
!!$    TPROC(:)%NZOR  = TPROC(:)%NZOR  + JPVEXT
!!$    TPROC(:)%NZEND = TPROC(:)%NZEND + JPVEXT
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE SPLITZ

  !     #######################################
  SUBROUTINE SET_NZ_PROC_ll(KZ_PROC)
    !     #######################################
    !
    !!****  *SET_NZ_PROC_ll* - 
    !
    !!    Purpose
    !!    -------
    !       Set the number of Proc in Z splitting
    !
    !-------------------------------------------------------------------------------
    !
    USE MODD_VAR_ll, ONLY : NZ_PROC
    !
    IMPLICIT NONE
    !
    INTEGER :: KZ_PROC
    !
    !-------------------------------------------------------------------------------
    !
    NZ_PROC  = KZ_PROC
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE SET_NZ_PROC_ll

  !###################################################################
  SUBROUTINE CARTESIANZ(TPROC,NB_PROC,X_DIM,Y_DIM,Z_DIM,X_DOMAINS,Y_DOMAINS,Z_DOMAINS)
    !###################################################################
    !
    !!****  *CARTESIAN* - routine which splits a domain if NB_PROC 
    !  		      is not a prime number
    !
    !!    Purpose
    !!    -------
    !     this routine fills the elements of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     J. ESCOBAR  * LA *
    !
    !!    Modifications
    !!    -------------
    !     Original 19/11/2007
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(ZONE_LL), INTENT(INOUT), DIMENSION(:),TARGET  :: TPROC                         ! split structur to file
    INTEGER      , INTENT(IN)                   :: NB_PROC                       ! number of processors
    INTEGER      , INTENT(IN)                   :: X_DIM,Y_DIM,Z_DIM             ! global data_grid dimension
    INTEGER      , INTENT(IN)                   :: X_DOMAINS,Y_DOMAINS,Z_DOMAINS ! processors_grid dimension
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: IDOM,II,IJ,IK
    !JUAN
    TYPE(ZONE_LL), POINTER  :: TP
    !JUAN
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    COMPUTE THE AVERAGE DIMENSION
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    FILL THE FIELDS OF TPROC
    ! 

    IDOM = 0
    DO IK=1,Z_DOMAINS
       DO IJ=1,Y_DOMAINS
          DO II=1,X_DOMAINS
             !
             !       file processor number
             !
             IDOM = IDOM + 1
             TP => TPROC(IDOM)
             TP%NUMBER = IDOM
             !
             !       compute/file x coordonate
             !
             CALL SLIDE_COORD(X_DIM,X_DOMAINS,II,TP%NXOR,TP%NXEND)
             !
             !       compute/file y coordonate
             !
             CALL SLIDE_COORD(Y_DIM,Y_DOMAINS,IJ,TP%NYOR,TP%NYEND)
             !
             !       compute/file Z coordonate
             !
             CALL SLIDE_COORD(Z_DIM,Z_DOMAINS,IK,TP%NZOR,TP%NZEND)
             !
          END DO
          !
       END DO

    END DO
    ! 
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE CARTESIANZ

  SUBROUTINE SLIDE_COORD(KDIM_DATA,KDIM_PROC,THIS_PROC,KOR,KEND)

    !!    Purpose
    !
    !   Compute for the processor=THIS_PROC the origine/end of slide in decomposing
    !   an array of data of dimension=KDIM_DATA on KDIM_PROC
    !
    !!    Author
    !!    ------
    !     J. ESCOBAR    * LA *

    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(IN)  :: KDIM_DATA ! dimension of data to split
    INTEGER, INTENT(IN)  :: KDIM_PROC ! numbers of processor to use in splitting
    INTEGER, INTENT(IN)  :: THIS_PROC ! processor id from 1..NB_PROC
    INTEGER, INTENT(OUT) :: KOR,KEND  ! Origine/End coordonate
    !
    !*       0.2   declarations of local variables
    !
    INTEGER             :: IDIM_SLIDE ! slide dimension ( without rest/delta )
    INTEGER             :: IREST      ! number of point in surabondance to distribut 
    INTEGER             :: IDELTAOR,IDELTAEND     ! offset in origine to apply 

    IDIM_SLIDE   = KDIM_DATA/KDIM_PROC
    IREST        = MOD(KDIM_DATA,KDIM_PROC)
    IDELTAOR     = MIN(IREST,THIS_PROC-1)
    IDELTAEND    = MIN(IREST,THIS_PROC)

    KOR   = ( THIS_PROC - 1 ) * IDIM_SLIDE + 1 + IDELTAOR
    KEND  =   THIS_PROC       * IDIM_SLIDE     + IDELTAEND

  END SUBROUTINE SLIDE_COORD


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
  !     ########################################
  SUBROUTINE INI_PZZ( TP_TSPLITS, TPPZS)
    !     ########################################
    !
    !!****  *INI_PZ* - routine to initialize the physical 2way splitting
    !                  of the TPPROCONF variable
    ! 
    !!    Purpose
    !!    -------
    !     the purpose of this routine is to fill the arguments of the
    !     variable TPPROCONF$TSPLITS_B concerning physical subdomains
    !     with a given splitting TPPZS
    !
    !!**  Method
    !!    ------
    ! 
    !!    External
    !!    --------
    ! 
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !        types PROCONF_ll, ZONE_ll
    !
    !     Module MODD_VAR_ll
    !        NPROC - Number of processors
    !
    !!    Reference
    !!    ---------
    ! 
    !!    Author
    !!    ------
    !     R. Guivarch               * CERFACS - ENSEEIHT *
    !!
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !     R. Guivarch 01/08/98 arguments for grid-nesting
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : PROCONF_ll, ZONE_ll, MODELSPLITTING_ll
    USE MODD_VAR_ll, ONLY       : NPROC
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(MODELSPLITTING_ll), DIMENSION(:), INTENT(INOUT):: TP_TSPLITS     ! Physical Zone Splitting
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(IN):: TPPZS     ! Physical Zone Splitting
    !
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: J ! loop control variable
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    FILL TPPROCONF%TSPLITS_B FOR EACH J :
    !              ---------------------------------------
    !
    DO J = 1, NPROC
       !
       TP_TSPLITS(J)%NUMBER = TPPZS(J)%NUMBER
       !
       TP_TSPLITS(J)%NXORP  = TPPZS(J)%NXOR
       TP_TSPLITS(J)%NXENDP = TPPZS(J)%NXEND
       TP_TSPLITS(J)%NDIMXP = TPPZS(J)%NXEND - TPPZS(J)%NXOR + 1
       !
       TP_TSPLITS(J)%NYORP  = TPPZS(J)%NYOR
       TP_TSPLITS(J)%NYENDP = TPPZS(J)%NYEND
       TP_TSPLITS(J)%NDIMYP = TPPZS(J)%NYEND - TPPZS(J)%NYOR + 1
       !
       TP_TSPLITS(J)%NZORP  = TPPZS(J)%NZOR
       TP_TSPLITS(J)%NZENDP = TPPZS(J)%NZEND
       TP_TSPLITS(J)%NDIMZP = TPPZS(J)%NZEND - TPPZS(J)%NZOR + 1
       !
    ENDDO
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE INI_PZZ
  !
  !##############################
  SUBROUTINE INI_EZZ( TPPROCONF )
    !############################
    !
    !!    Author
    !!    ------
    !     J. ESCOBAR              * LA - CNRS *
    ! 
    !!    Modifications
    !!    -------------
    !     Original 22/01/2008
    ! 
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
    !
    !-------------------------------------------------------------------------------
    !
    CALL INI_EZZZ( TPPROCONF%TBOUND_SXP2_YP1_Z, TPPROCONF%TSPLITS_SXP2_YP1_Z)
    CALL INI_EZZZ( TPPROCONF%TBOUND_SX_YP2_ZP1, TPPROCONF%TSPLITS_SX_YP2_ZP1)
    CALL INI_EZZZ( TPPROCONF%TBOUND_SXP2_Y_ZP1, TPPROCONF%TSPLITS_SXP2_Y_ZP1)
    !-------------------------------------------------------------------------------
    !
  CONTAINS
    SUBROUTINE INI_EZZZ(TP_TBOUNDS,TP_TSPLITS)
      !
      !*       0.    DECLARATIONS
      !
      USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT ! Horizontal/Vertical External points number
      USE MODD_STRUCTURE_ll , ONLY : MODELSPLITTING_ll, LOCALISATION_ll
      USE MODD_VAR_ll       , ONLY : NPROC,  JPHALO
      !
      IMPLICIT NONE
      !
      !*       0.1   declarations of arguments
      !
      TYPE(LOCALISATION_ll)   , DIMENSION(:), POINTER :: TP_TBOUNDS
      TYPE(MODELSPLITTING_ll) , DIMENSION(:), POINTER :: TP_TSPLITS
      !
      !*       0.2   declarations of local variables
      !
      INTEGER                 :: J ! loop control variable
      INTEGER                 :: IFLAG
      IFLAG=0
      DO J = 1, NPROC
         !
         ! Origine indexe 
         !
         IF(TP_TBOUNDS(J)%WEST) THEN
            TP_TSPLITS%NXORE = TP_TSPLITS%NXORP - JPHEXT*IFLAG
         ELSE
            TP_TSPLITS%NXORE = TP_TSPLITS%NXORP - JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%SOUTH) THEN
            TP_TSPLITS%NYORE = TP_TSPLITS%NYORP - JPHEXT*IFLAG
         ELSE
            TP_TSPLITS%NYORE = TP_TSPLITS%NYORP - JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%BOTTOM) THEN
            TP_TSPLITS%NZORE = TP_TSPLITS%NZORP - JPVEXT*IFLAG
         ELSE
            TP_TSPLITS%NZORE = TP_TSPLITS%NZORP - JPHALO*IFLAG
         ENDIF
         !
         ! End indexe
         !
         IF(TP_TBOUNDS(J)%EAST) THEN
            TP_TSPLITS%NXENDE = TP_TSPLITS%NXENDP + JPHEXT*IFLAG
         ELSE
            TP_TSPLITS%NXENDE = TP_TSPLITS%NXENDP + JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%NORTH) THEN
            TP_TSPLITS%NYENDE = TP_TSPLITS%NYENDP + JPHEXT*IFLAG
         ELSE
            TP_TSPLITS%NYENDE = TP_TSPLITS%NYENDP + JPHALO*IFLAG
         ENDIF
         !
         IF(TP_TBOUNDS(J)%TOP) THEN
            TP_TSPLITS%NZENDE = TP_TSPLITS%NZENDP + JPVEXT*IFLAG
         ELSE
            TP_TSPLITS%NZENDE = TP_TSPLITS%NZENDP + JPHALO*IFLAG
         ENDIF
         !
         !
         ! Size indexe
         !
         TP_TSPLITS%NDIMXE = TP_TSPLITS%NXENDE - TP_TSPLITS%NXORE + 1
         TP_TSPLITS%NDIMYE = TP_TSPLITS%NYENDE - TP_TSPLITS%NYORE + 1
         TP_TSPLITS%NDIMZE = TP_TSPLITS%NZENDE - TP_TSPLITS%NZORE + 1
         !
      ENDDO
    END SUBROUTINE INI_EZZZ
  END SUBROUTINE INI_EZZ

  !     #######################################
  SUBROUTINE INI_BOUNDARIESZ( TPPROCONF )
    !####################################
    !
    !!****  *INI_BOUNDARIESZ* - routine to compute the localisation variable TBOUND
    !
    !!
    !!    Author
    !!    ------
    !!    J. ESCOBAR
    !!
    !!    Modifications
    !!    -------------
    !!    Original 22/01/2008
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY  : PROCONF_ll
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PROCONF_ll), POINTER :: TPPROCONF       ! splitting data structure
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    FILL TBOUND OF EACH PROCESSOR :
    !              ------------------------------
    !
    CALL INI_BOUNZ( TPPROCONF%TBOUND_SXP2_YP1_Z, TPPROCONF%TSPLITS_SXP2_YP1_Z)
    CALL INI_BOUNZ( TPPROCONF%TBOUND_SX_YP2_ZP1, TPPROCONF%TSPLITS_SX_YP2_ZP1)
    CALL INI_BOUNZ( TPPROCONF%TBOUND_SXP2_Y_ZP1, TPPROCONF%TSPLITS_SXP2_Y_ZP1)
    !
    !-------------------------------------------------------------------------------
    !
  CONTAINS
    SUBROUTINE INI_BOUNZ(TP_TBOUNDS,TP_TSPLITS)
      !
      !*       0.    DECLARATIONS
      !
      USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT ! Horizontal/Vertical External points number
      USE MODD_STRUCTURE_ll , ONLY : MODELSPLITTING_ll, LOCALISATION_ll
      USE MODD_VAR_ll       , ONLY : NPROC, DIMX, DIMY, DIMZ
      !
      IMPLICIT NONE
      !
      !*       0.1   declarations of arguments
      !
      TYPE(LOCALISATION_ll)   , DIMENSION(:), POINTER :: TP_TBOUNDS
      TYPE(MODELSPLITTING_ll) , DIMENSION(:), POINTER :: TP_TSPLITS
      !
      !*       0.2   declarations of local variables
      !
      INTEGER                 :: J ! loop control variable
      !
      ALLOCATE(TP_TBOUNDS(NPROC))
      !
      DO J = 1, NPROC
         !
         !
         TP_TBOUNDS(J) = &
              LOCALISATION_ll( .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.)
         !
         IF( TP_TSPLITS(J)%NXORP .EQ. (1+JPHEXT)) &
              TP_TBOUNDS(J)%WEST = .TRUE.
         !
         IF( TP_TSPLITS(J)%NXENDP  .EQ. (DIMX-JPHEXT)) &
              TP_TBOUNDS(J)%EAST = .TRUE.
         !
         IF( TP_TSPLITS(J)%NYORP .EQ. (1+JPHEXT)) &
              TP_TBOUNDS(J)%SOUTH = .TRUE.
         !
         IF( TP_TSPLITS(J)%NYENDP  .EQ. (DIMY-JPHEXT)) &
              TP_TBOUNDS(J)%NORTH = .TRUE.
         !
         IF( TP_TSPLITS(J)%NZORP .EQ. (1+JPVEXT)) &
              TP_TBOUNDS(J)%BOTTOM = .TRUE.
         !
         IF( TP_TSPLITS(J)%NZENDP  .EQ. (DIMZ-JPVEXT)) &
              TP_TBOUNDS(J)%TOP = .TRUE.
      ENDDO
    END SUBROUTINE INI_BOUNZ
  END SUBROUTINE INI_BOUNDARIESZ

  !
  !     ##################################################
  SUBROUTINE CONSTRUCT_TRANSZ(TPCOMDATA, TPPROCONF )
    !     ##################################################
    !
    !!****  *CONSTRUCT_TRANS* - routine to construct
    !                           the transposition correspondants
    !
    !!    Purpose
    !!    -------
    !     the purpose of the routine is to fill the structured type variable
    !     TPCOMDATA with informations concerning the communications
    !     during transposition
    !
    !!**  Method
    !!    ------
    !     we compute for the local processor,
    !      - intersections between zones of the global domain in x-slice splitting
    !        and local zone in 2way splitting to find the send correspondant
    !        of the local processor for transposition 2way / x-slices
    ! 
    !      - intersections between physical zones of the global domain
    !        in 2way splitting and local zone in x-slices splitting
    !        to find the receive correspondant of the local processor
    !        for transposition 2way / x-slices
    ! 
    !      - intersections between zones of the global domain in y-slice splitting
    !        and local zone in x-slices splitting to find the send correspondant
    !        of the local processor for transposition x-slices / y-slices
    ! 
    !      - intersections between zones of the global domain
    !        in x-slices splitting and local zone in y-slices splitting
    !        to find the receive correspondant of the local processor
    !        for transposition x-slices / y-slices
    ! 
    !!    External
    !!    --------
    !     Module MODE_TOOLS_ll
    !        ADD_ZONE, INTERSECTION, GLOBAL2LOCAL, EXTRACT_ZONE, G2LX
    !
    !     Module MODE_CONSTRUCT_ll
    !        RESIZE_TRANS
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !        types ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
    !
    !     Module MODD_VAR_ll
    !        IP - Number of local processor=subdomain
    !        NPROC - Number of processors
    !
    !!    Reference
    !!    ---------
    ! 
    !!    Author
    !!    ------
    !     R. Guivarch               * CERFACS - ENSEEIHT *
    ! 
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !     R. Guivarch 01/08/98 arguments for grid-nesting
    !!
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, PROC_COM_DATA_ll, PROCONF_ll
    USE MODD_VAR_ll, ONLY       : IP, NPROC
    !
    USE MODE_TOOLS_ll, ONLY     : INTERSECTION, G2LX, GLOBAL2LOCAL, ADD_ZONE, &
         EXTRACT_ZONE
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PROC_COM_DATA_ll), POINTER :: TPCOMDATA ! communications data
    TYPE(PROCONF_ll), POINTER       :: TPPROCONF ! splitting data
    !
    !*       0.2   declarations of local variables
    !
    !
    !-------------------------------------------------------------------------------
    !
    !        1.1    CONSTRUCTION OF TRANSPOSITION B -> SX_YP2_ZP1-SLICES CORRESPONDANTS :
    !              -------------------------------------------------------------
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,1000,TPCOMDATA%TSEND_B_SX_YP2_ZP1,&
         TPCOMDATA%TSEND_BOX_B_SX_YP2_ZP1,&
         TPPROCONF%TSPLITS_B,TPPROCONF%TSPLITS_SX_YP2_ZP1)
    !CALL G2LXZ(TPPROCONF%TSPLITS_B(IP),TPCOMDATA%TSEND_B_SX_YP2_ZP1)
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,1000,TPCOMDATA%TRECV_B_SX_YP2_ZP1,&
         TPCOMDATA%TRECV_BOX_B_SX_YP2_ZP1,&     
         TPPROCONF%TSPLITS_SX_YP2_ZP1,TPPROCONF%TSPLITS_B)
    !-------------------------------------------------------------------------------
    !
    !*       1.2   CONSTRUCTION OF
    !              TRANSPOSITION SXP2_Y_ZP1-SLICES -> B-SLICES CORRESPONDANTS :
    !              -------------------------------------------------
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,2000,TPCOMDATA%TSEND_SXP2_Y_ZP1_B,&
         TPCOMDATA%TSEND_BOX_SXP2_Y_ZP1_B,&
         TPPROCONF%TSPLITS_SXP2_Y_ZP1,TPPROCONF%TSPLITS_B)
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,2000,TPCOMDATA%TRECV_SXP2_Y_ZP1_B,&
         TPCOMDATA%TRECV_BOX_SXP2_Y_ZP1_B,&
         TPPROCONF%TSPLITS_B,TPPROCONF%TSPLITS_SXP2_Y_ZP1)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.    CONSTRUCTION OF TRANSPOSITION SXP2_YP1_Z -> SX_YP2_ZP1-SLICES CORRESPONDANTS :
    !              -------------------------------------------------------------
    !
!!$    CALL CONSTRUCT_TRANSZZ(NPROC,3000,TPCOMDATA%TSEND_SXP2_YP1_Z_SX_YP2_ZP1,&
!!$         TPCOMDATA%TSEND_BOX_,&
!!$         TPPROCONF%TSPLITS_SXP2_YP1_Z,TPPROCONF%TSPLITS_SX_YP2_ZP1)
!!$    !
!!$    CALL CONSTRUCT_TRANSZZ(NPROC,3000,TPCOMDATA%TRECV_SXP2_YP1_Z_SX_YP2_ZP1,&
!!$         TPCOMDATA%TRECV_BOX_,&
!!$         TPPROCONF%TSPLITS_SX_YP2_ZP1,TPPROCONF%TSPLITS_SXP2_YP1_Z)
    !
    !-------------------------------------------------------------------------------
    !
    !*       3.    CONSTRUCTION OF
    !              TRANSPOSITION SX_YP2_ZP1-SLICES -> SXP2_Y_ZP1-SLICES CORRESPONDANTS :
    !              -------------------------------------------------
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,4000,TPCOMDATA%TSEND_SX_YP2_ZP1_SXP2_Y_ZP1,&
         TPCOMDATA%TSEND_BOX_SX_YP2_ZP1_SXP2_Y_ZP1,&
         TPPROCONF%TSPLITS_SX_YP2_ZP1,TPPROCONF%TSPLITS_SXP2_Y_ZP1)
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,4000,TPCOMDATA%TRECV_SX_YP2_ZP1_SXP2_Y_ZP1,&
         TPCOMDATA%TRECV_BOX_SX_YP2_ZP1_SXP2_Y_ZP1,&
         TPPROCONF%TSPLITS_SXP2_Y_ZP1,TPPROCONF%TSPLITS_SX_YP2_ZP1)
    !
    !-------------------------------------------------------------------------------
    !
    !*       4.    CONSTRUCTION OF
    !              TRANSPOSITION SXP2_Y_ZP1-SLICES -> SXP2_YP1_Z-SLICES CORRESPONDANTS :
    !              -------------------------------------------------
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,5000,TPCOMDATA%TSEND_SXP2_Y_ZP1_SXP2_YP1_Z,&
         TPCOMDATA%TSEND_BOX_SXP2_Y_ZP1_SXP2_YP1_Z,&
         TPPROCONF%TSPLITS_SXP2_Y_ZP1,TPPROCONF%TSPLITS_SXP2_YP1_Z)
    !
    CALL CONSTRUCT_TRANSZZ(NPROC,5000,TPCOMDATA%TRECV_SXP2_Y_ZP1_SXP2_YP1_Z,&
         TPCOMDATA%TRECV_BOX_SXP2_Y_ZP1_SXP2_YP1_Z,&
         TPPROCONF%TSPLITS_SXP2_YP1_Z,TPPROCONF%TSPLITS_SXP2_Y_ZP1)
    !
    !*       3.1   extraction of x-slices splitting
    !
    !
    !
    !*       3.5   Switch to local coordinates
    !
    !CALL G2LX(TPPROCONF%TSPLITS_X(IP),TPCOMDATA%TSEND_TRANS_XY)
    !CALL G2LX(TPPROCONF%TSPLITS_Y(IP),TPCOMDATA%TRECV_TRANS_XY)
    !
    !
    !-------------------------------------------------------------------------------
    !
    !*       4.    DEALLOCATION OF LOCAL VARIABLES :
    !              -------------------------------
    !    
    !
    !-------------------------------------------------------------------------------
    !
  CONTAINS
    SUBROUTINE CONSTRUCT_TRANSZZ(KPROC,KMESG,TP_TRANS_FROM_TO,&
         TBOX_FROM_TO,&
         TP_TSPLITS_FROM,TP_TSPLITS_TO)
      !
      USE MODD_STRUCTURE_ll , ONLY : CRSPD_ll, MODELSPLITTING_ll, BOX_ll
      USE MODD_VAR_ll       , ONLY : IP
      IMPLICIT NONE
      !
      ! Argument
      !
      INTEGER                                         :: KPROC,KMESG
      TYPE(CRSPD_ll)                        , POINTER :: TP_TRANS_FROM_TO
      TYPE(BOX_ll)                          , POINTER :: TBOX_FROM_TO
      TYPE(MODELSPLITTING_ll),  DIMENSION(:), POINTER :: TP_TSPLITS_FROM,TP_TSPLITS_TO
      !
      ! Local Variable
      !
      !TYPE(ZONE_ll), ALLOCATABLE, DIMENSION(:) :: TP_ZONE_FROM,TP_ZONE_TO,TZINTER
      TYPE(ZONE_ll)          , DIMENSION(KPROC) :: TP_ZONE_FROM,TP_ZONE_TO,TZINTER
      INTEGER                                   :: JI
      !-------------------------------------------------------------------------------
      !
      !*       1.    ALLOCATE OF THE LOCAL VARIABLES :
      !              -------------------------------
      !
      !ALLOCATE( TZPZS(KPROC), TZTRANSXZS(KPROC), &
      !     TZTRANSYZS(KPROC), TZINTER(KPROC) )
      !
      !-------------------------------------------------------------------------------
      !
      !*       2.    CONSTRUCTION OF TRANSPOSITION 2WAY -> X-SLICES CORRESPONDANTS :
      !              -------------------------------------------------------------
      !
      !*       2.1   extraction of physical 2way splitting
      !              from TPPROCONF structure
      !
      CALL EXTRACT_ZONEZ(TP_TSPLITS_FROM, TP_ZONE_FROM, TZINTER)
      !
      !*       2.2   extraction of x-slices splitting
      !              from TPPROCONF structure
      !
      CALL EXTRACT_ZONEZ(TP_TSPLITS_TO, TP_ZONE_TO, TZINTER)
      !
      !*       2.3   computation of intersection between local 2way zone
      !              and x-slices splitting -> send correspondants
      !
      CALL INTERSECTIONZ(TP_ZONE_TO, KPROC, TP_ZONE_FROM(IP), TZINTER,&
                         TBOX_FROM_TO)
      !
      NULLIFY(TP_TRANS_FROM_TO)
      DO JI = 1, KPROC
         IF(TZINTER(JI)%NUMBER.NE.0) THEN
            !TZINTER(JI)%MSSGTAG = KMESG 
            TZINTER(JI)%MSSGTAG = KMESG + JI*0000 + TZINTER(JI)%NUMBER*0
            CALL ADD_ZONE(TP_TRANS_FROM_TO, TZINTER(JI) )
         ENDIF
      ENDDO
      !
      CALL G2LXZ(TP_TSPLITS_FROM(IP),TP_TRANS_FROM_TO,TBOX_FROM_TO)
      !
    END SUBROUTINE CONSTRUCT_TRANSZZ
  END SUBROUTINE CONSTRUCT_TRANSZ
  !
  !     ####################################################
  SUBROUTINE INTERSECTIONZ( TPSPLIT, K, TPZONE, TPRES,&
                            TBOX_FROM_TO)
    !     ####################################################
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR
    !!
    !!    Modifications
    !!    -------------
    !     Original 23/01/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_ll, BOX_ll, ALLOC_BOX_LL
    USE MODD_VAR_ll, ONLY : DIMZ
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(IN) :: TPSPLIT ! Splitting of the domain
    !
    INTEGER, INTENT(IN) :: K ! Number of elements of TPSPLIT
    !
    TYPE(ZONE_ll), INTENT(IN) :: TPZONE ! Zone to be splitted
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(OUT) :: TPRES ! Splitting of the zone
    !
     TYPE(BOX_ll)              , POINTER     :: TBOX_FROM_TO
    !*       0.2   declarations of local variables
    !
    INTEGER :: J      ! loop control variable
    INTEGER :: JCNTS  ! Siez of the current box
    INTEGER :: JSTRT  ! displacement in the pack box 
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    LIST AND COMPUTE INTERSECTION BETWEEN TPSPLIT(J) AND TPZONE :
    !              -----------------------------------------------------------
    !
    CALL ALLOC_BOX_ll(TBOX_FROM_TO,K) 
    !
    JSTRT = 0
    !
    DO J = 1, K
       ! 
       !   Which subdomain is the owner of TPSPLIT(J)
       TPRES(J)%NUMBER = TPSPLIT(J)%NUMBER
       !
       !   Computation of the origin coordinate
       TPRES(J)%NXOR = MAX( TPZONE%NXOR, TPSPLIT(J)%NXOR )
       TPRES(J)%NYOR = MAX( TPZONE%NYOR, TPSPLIT(J)%NYOR )
       TPRES(J)%NZOR = MAX( TPZONE%NZOR, TPSPLIT(J)%NZOR )
       !
       !   Computation of the last coordinate
       TPRES(J)%NXEND = MIN( TPZONE%NXEND, TPSPLIT(J)%NXEND )
       TPRES(J)%NYEND = MIN( TPZONE%NYEND, TPSPLIT(J)%NYEND )
       TPRES(J)%NZEND = MIN( TPZONE%NZEND, TPSPLIT(J)%NZEND )
       !
       !   for z-direction all the domain is considered
!!$    TPRES(J)%NZOR = 1
!!$    TPRES(J)%NZEND = DIMZ
       ! 
       !   if the intersection is void, the result is nullified
       IF((TPRES(J)%NXOR > TPRES(J)%NXEND) .OR. &
            (TPRES(J)%NYOR > TPRES(J)%NYEND) .OR. &
            (TPRES(J)%NZOR > TPRES(J)%NZEND) ) THEN
          
          TPRES(J) = ZONE_ll ( 0, 0, 0, 0, 0, 0, 0, 0 )
          JCNTS    = 0
       ELSE                    
          JCNTS = ( TPRES(J)%NXEND - TPRES(J)%NXOR + 1 ) * &
                  ( TPRES(J)%NYEND - TPRES(J)%NYOR + 1 ) * &     
                  ( TPRES(J)%NZEND - TPRES(J)%NZOR + 1 )
       END IF
       ! 
       ! Fill the box_ll for mpialltoallv
       !
       !
       TBOX_FROM_TO%NXOR(J)  = TPRES(J)%NXOR
       TBOX_FROM_TO%NXEND(J) = TPRES(J)%NXEND
       !
       TBOX_FROM_TO%NYOR(J)  = TPRES(J)%NYOR
       TBOX_FROM_TO%NYEND(J) = TPRES(J)%NYEND
       !
       TBOX_FROM_TO%NZOR(J)  = TPRES(J)%NZOR
       TBOX_FROM_TO%NZEND(J) = TPRES(J)%NZEND
       !
       TBOX_FROM_TO%NCNT(J)  = JCNTS
       TBOX_FROM_TO%NSTRT(J) = JSTRT
       JSTRT = JSTRT + JCNTS ! next box start pointer
       !
    ENDDO
    TBOX_FROM_TO%NSIZE = JSTRT
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE INTERSECTIONZ

  !     ######################################################
  SUBROUTINE REMAP_SXP2_YP1_Z_SX_YP2_ZP1_ll(PFIELDIN, PFIELDOUT, KINFO)
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
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_SXP2_YP1_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_SXP2_YP1_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_SXP2_YP1_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TRECV_SXP2_YP1_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SXP2_YP1_Z_SX_YP2_ZP1_ll

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
    USE MODD_CONF       , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
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

  !     ######################################################
  SUBROUTINE REMAP_SX_YP2_ZP1_SXP2_YP1_Z_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SX_YP2_ZP1_B_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in B splitting to PFIELDOUT in SX_YP2_ZP1 splitting
    !vi run
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
    USE MODD_CONF       , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
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
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TRECV_SXP2_YP1_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_SXP2_YP1_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TRECV_SXP2_YP1_Z_SX_YP2_ZP1, &
         TCRRT_COMDATA%TSEND_SXP2_YP1_Z_SX_YP2_ZP1, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    !
    !-------------------------------------------------------------------------------
    !
    !
    !-------------------------------------------------------------------------------
    !
!!$    CALL ALL_SEND_RECV(TCRRT_COMDATA%TSEND_BOX_SXP2_YP1_Z_SX_YP2_ZP1, &
!!$         TCRRT_COMDATA%TRECV_BOX_SXP2_YP1_Z_SX_YP2_ZP1, &
!!$         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
  END SUBROUTINE REMAP_SX_YP2_ZP1_SXP2_YP1_Z_ll

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
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SX_YP2_ZP1_B_ll

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
    USE MODD_CONF       , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
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
    USE MODD_CONF       , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
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
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SXP2_Y_ZP1_SXP2_YP1_Z_ll

  !     ######################################################
  SUBROUTINE REMAP_SXP2_Y_ZP1_B_ll(PFIELDIN, PFIELDOUT, KINFO)
    !     ######################################################
    !
    !!**** *REMAP_SXP2_Y_ZP1_B_ll* -
    !!
    !!    Purpose
    !!    -------
    !       remaps PFIELDIN from in SXP2_Y_ZP1 splitting to PFIELDOUT in B splitting
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
    USE MODD_CONF       , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
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
    CALL SEND_RECV_FIELD(TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_B, &
         TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_B, &
         PFIELDIN, PFIELDOUT, KINFO)
    !
    !-------------------------------------------------------------------------------
    !
    !*       2.     UPDATE THE ZONES THE PROCESSOR SENDS OR RECEIVED FROM ITSELF
    !               ------------------------------------------------------------
    !
    CALL COPY_CRSPD_TRANS(TCRRT_COMDATA%TSEND_SXP2_Y_ZP1_B, &
         TCRRT_COMDATA%TRECV_SXP2_Y_ZP1_B, &
         PFIELDIN, PFIELDOUT, KINFO)
    ELSE
    !
    !-------------------------------------------------------------------------------
    !
    CALL ALL_SEND_RECV(TCRRT_COMDATA%TSEND_BOX_SXP2_Y_ZP1_B, &
         TCRRT_COMDATA%TRECV_BOX_SXP2_Y_ZP1_B, &
         PFIELDIN, PFIELDOUT, KINFO)
    ENDIF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE REMAP_SXP2_Y_ZP1_B_ll

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
    USE MODD_CONF       , ONLY : NZ_SPLITTING ! for debug IZ=1=flat_inv;  IZ=2=flat_invz ;  IZ=1+2=the two
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

  !     #################################################
  SUBROUTINE EXTRACT_ZONEZ( TPSPLITS, TPPZS, TPEZS )
    !     #################################################
    !
    !!****  *EXTRACT_ZONEZ* - routine to construct two splittings variables
    !!                       from a MODELSPLITTING_ll variable
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR   LA - CNRS
    ! 
    !!    Modifications
    !!    -------------
    !     Original 28/01/2008
    ! 
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll, ZONE_ll
    USE MODD_VAR_ll, ONLY : NPROC
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(MODELSPLITTING_ll), DIMENSION(:), POINTER :: TPSPLITS
    !
    TYPE(ZONE_ll), DIMENSION(:), INTENT(OUT) :: TPPZS, TPEZS
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: J ! loop control variable
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    FILL TPPZS AND TPEZS FOR EACH J :
    !              -------------------------------
    !
    DO J = 1, NPROC
       !
       TPPZS(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
       TPEZS(J) = ZONE_ll( 0, 0, 0, 0, 0, 0, 0, 0 )
       !
       TPPZS(J)%NUMBER = TPSPLITS(J)%NUMBER
       TPPZS(J)%NXOR   = TPSPLITS(J)%NXORP
       TPPZS(J)%NYOR   = TPSPLITS(J)%NYORP
       TPPZS(J)%NXEND  = TPSPLITS(J)%NXENDP
       TPPZS(J)%NYEND  = TPSPLITS(J)%NYENDP
       !JUAN
!!$       TPPZS(J)%NZOR   = TPSPLITS(J)%NZORP
!!$       TPPZS(J)%NZEND  = TPSPLITS(J)%NZENDP
       ! JUAN :: For Z splitting extended domain is needed in Z slide
       TPPZS(J)%NZOR   = TPSPLITS(J)%NZORE
       TPPZS(J)%NZEND  = TPSPLITS(J)%NZENDE
       !JUAN
       !
       TPEZS(J)%NUMBER = TPSPLITS(J)%NUMBER
       TPEZS(J)%NXOR   = TPSPLITS(J)%NXORE
       TPEZS(J)%NYOR   = TPSPLITS(J)%NYORE
       TPEZS(J)%NXEND  = TPSPLITS(J)%NXENDE
       TPEZS(J)%NYEND  = TPSPLITS(J)%NYENDE
       !JUAN
       TPEZS(J)%NZOR   = TPSPLITS(J)%NZORE
       TPEZS(J)%NZEND  = TPSPLITS(J)%NZENDE
       !JUAN
    ENDDO
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE EXTRACT_ZONEZ

  !     ################################
  SUBROUTINE G2LXZ(TPSPLIT,TPCRSPD,TBOX)
    !     ################################
    !
    !!****  *G2LXZ* - routine to switch from global coordinates to local ones
    ! 
    !
    !!    Reference
    !!    ---------
    ! 
    !!    Author
    !!    ------
    !     J. ESCOBAR
    !
    !!    Modifications
    !!    -------------
    !     Original 01/02/2008
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : MODELSPLITTING_ll, CRSPD_ll, ZONE_ll, BOX_ll
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(MODELSPLITTING_ll), INTENT(IN) :: TPSPLIT ! x-slices or y-slices
    ! splitting
    !
    TYPE(CRSPD_ll), POINTER :: TPCRSPD ! CRSPD_ll to be switch
    TYPE(BOX_ll)  , POINTER :: TBOX
    !
    !*       0.2   declarations of local variables
    !
    ! intermediate variables to describe the list and the elements
    TYPE(ZONE_ll), POINTER :: TZZONE
    TYPE(CRSPD_ll), POINTER :: TZCRSPD
    !
    !-------------------------------------------------------------------------------
    !
    !*       1.    SWITCH :
    !              ------
    ! we list the variable TPCRSPD of type CRSPD_ll
    TZCRSPD => TPCRSPD
    DO WHILE (ASSOCIATED(TZCRSPD))
       TZZONE => TZCRSPD%TELT
       !
       !   we substract the origin of the local subdomain (physical subdomain)
!!$       TZZONE%NXOR = TZZONE%NXOR - TPSPLIT%NXORP + 1
!!$       TZZONE%NXEND = TZZONE%NXEND - TPSPLIT%NXORP + 1
!!$       TZZONE%NYOR = TZZONE%NYOR - TPSPLIT%NYORP + 1
!!$       TZZONE%NYEND = TZZONE%NYEND - TPSPLIT%NYORP + 1
!!$       !JUAN
!!$       TZZONE%NZOR = TZZONE%NZOR - TPSPLIT%NZORP + 1
!!$       TZZONE%NZEND = TZZONE%NZEND - TPSPLIT%NZORP + 1
       !JUAN
       !JUAN
       !   we substract the origin of the local subdomain (extended subdomain)
       TZZONE%NXOR = TZZONE%NXOR - TPSPLIT%NXORE + 1
       TZZONE%NXEND = TZZONE%NXEND - TPSPLIT%NXORE + 1
       !
       TZZONE%NYOR = TZZONE%NYOR - TPSPLIT%NYORE + 1
       TZZONE%NYEND = TZZONE%NYEND - TPSPLIT%NYORE + 1
       !JUAN
       TZZONE%NZOR = TZZONE%NZOR - TPSPLIT%NZORE + 1
       TZZONE%NZEND = TZZONE%NZEND - TPSPLIT%NZORE + 1
       !JUAN
       !
       TZCRSPD => TZCRSPD%TNEXT
    ENDDO
    TBOX%NXOR(:)  =  TBOX%NXOR(:)  -  TPSPLIT%NXORE + 1
    TBOX%NXEND(:) =  TBOX%NXEND(:) -  TPSPLIT%NXORE + 1
    !
    TBOX%NYOR(:)  =  TBOX%NYOR(:)  -  TPSPLIT%NYORE + 1
    TBOX%NYEND(:) =  TBOX%NYEND(:) -  TPSPLIT%NYORE + 1
    !
    TBOX%NZOR(:)  =  TBOX%NZOR(:)  -  TPSPLIT%NZORE + 1
    TBOX%NZEND(:) =  TBOX%NZEND(:) -  TPSPLIT%NZORE + 1
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE G2LXZ

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
      USE MODD_MPIF         , ONLY : MPI_COMM_WORLD
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
                         MPI_COMM_WORLD,KINFO)
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
