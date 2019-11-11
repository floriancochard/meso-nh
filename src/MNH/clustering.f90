!    ################
     MODULE MODI_CLUSTERING
!    ################
!
INTERFACE
!
      SUBROUTINE CLUSTERING(OBOTTOMUP,OCLOUD,PFIELD,KCLUSTERIDT,KCLUSTERLVL,PCLUSTERSEC)
!
LOGICAL,                   INTENT(IN)  :: OBOTTOMUP
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OCLOUD 
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PFIELD 
INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: KCLUSTERIDT
INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: KCLUSTERLVL
REAL,    DIMENSION(:,:,:), INTENT(OUT) :: PCLUSTERSEC
!
END SUBROUTINE CLUSTERING
!
END INTERFACE
!
END MODULE MODI_CLUSTERING
!
!     #######################
      SUBROUTINE CLUSTERING(OBOTTOMUP,OCLOUD,PFIELD,KCLUSTERIDT,KCLUSTERLVL,PCLUSTERSEC)
!     #######################
      !
!!    PURPOSE
!!    -------
!!    
!!    Identify structures as 3D objects made of connected points.
!!
!!**  METHOD
!!    ------
!!
!!    Uses a 3D mask as input. Outputs 3D fields containing:
!!      IDT the identitdy number of each structure
!!      LVL the level from which each structure is identified
!!      SEC the section of each structure at current level
!!
!!    Both IDT and LVL are necessary to identify unequivocally the structures.
!!
!!    The identification is first led level by level.
!!    The connection between the levels is done after.
!!
!!
!!    AUTHOR
!!    ------
!!      T.    Dauhut                     * LA *
!!      J.-P. Chaboureau                 * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/06/13
!!      Modified    05/08/14    T. Dauhut   v75
!!      Modified    10/10/14    T. Dauhut   toward 3D version
!!      Modified    13/11/14    T. Dauhut   adding property field analyse
!!      Modified    13/06/17    T. Dauhut   to start volume scan from top
!!      Modified    04/10/17    T. Dauhut   to be added to next MNH versions
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
 USE MODD_MPIF , ONLY : MPI_INTEGER
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_VAR_ll, ONLY : MPI_PRECISION, NPROC, IP, NMNH_COMM_WORLD
USE MODD_DYN_n, ONLY : XDXHATM, XDYHATM
!USE&
!      MPI
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
LOGICAL,                   INTENT(IN)  :: OBOTTOMUP
LOGICAL, DIMENSION(:,:,:), INTENT(IN)  :: OCLOUD 
REAL   , DIMENSION(:,:,:), INTENT(IN)  :: PFIELD 
INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: KCLUSTERIDT
INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: KCLUSTERLVL
REAL,    DIMENSION(:,:,:), INTENT(OUT) :: PCLUSTERSEC
!
!*       0.2   declarations of local variables
!
INTEGER, PARAMETER :: JPUNDEF = 999999            ! NaN value outside clusters
INTEGER            :: IKBEG,IKEND,IKINC           ! for the loop on the levels
INTEGER            :: JI,JJ,JC,JK,JH,JDI,JDJ      ! loop counters
INTEGER            :: IMINI, IMIX, IMUX           ! local cluster index
INTEGER            :: IIB,IJB,IIE,IJE,IKE,IXOR,IYOR  ! physical bounds of local domain
INTEGER            :: IIU, IJU              ! dimensions of local domain+halo
INTEGER            :: NIMAX_ll, NJMAX_ll    ! dimensions of global domain
INTEGER            :: ICLUSMAX              ! maximum nbr of clusters in local domain
INTEGER            :: ICPT, ICHANGES, ICPR  ! local counter
INTEGER            :: ITEMPIDT,ITEMPLVL,ITEMPAFL ! temporary cluster ID, baselevel, affluent section
REAL               :: ZMESHAREA
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ICLUSIND ! 1D table for local indices
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ICLUSIZE ! 1D table for local cluster size in pixels
INTEGER, ALLOCATABLE, DIMENSION(:,:)  ::  IMAPIND  ! local 2D field pointing toward local indices
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ICLUSTERIDT ! 1D table for cluster ID (need LVL to be unequivocal)
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ICLUSTERLVL ! 1D table for cluster base level
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  IAFFLUSECTN ! 1D table for cluster affluent section
INTEGER, DIMENSION(SIZE(KCLUSTERIDT,1),SIZE(KCLUSTERIDT,2)) :: IMAPIDT_ll,IMAPLVL_ll,IMAPAFL_ll
REAL,    DIMENSION(SIZE(KCLUSTERIDT,1),SIZE(KCLUSTERIDT,2)) :: ZMAPIDT_ll,ZMAPLVL_ll,ZMAPAFL_ll
REAL,    ALLOCATABLE, DIMENSION(:)    ::  ZCLUSSUMFLD ! 1D table for cluster local sum of field values
!
!*       0.3 declarations of variables for MPI exchanges
!
INTEGER                 :: ITOTNBR                   ! counter to create final tables (.GLBLIST..)
INTEGER                 :: IRANK, IROOT, INFO        ! for MPI subroutines
INTEGER                 :: IINFO_ll                  ! (MNH-MPI) return code of parallel routine
TYPE(LIST_ll), POINTER  :: TZFIELDS_ll               ! (MNH-MPI) list of fields to exchange
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ICLUSNBR   ! 1D table for nbr of clusters in local domains
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  IPROCDPL   ! 1D table for procs starting indices in global tables
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ILOCIND    ! old indices point toward new ones after optimization
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ILOCLISTIDT,ILOCLISTIDT2 ! non-redundant list of local cluster IDs
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ILOCLISTLVL,ILOCLISTLVL2 ! their corresponding base level
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  ILOCLISTSEC,ILOCLISTSEC2 ! their corresponding sections
REAL   , ALLOCATABLE, DIMENSION(:)    ::  ZLOCLISTFLD,ZLOCLISTFLD2 ! their corresponding field average
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  IGLBLISTIDT ! concatenation of local domain cluster ID lists
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  IGLBLISTLVL ! their corresponding base level
INTEGER, ALLOCATABLE, DIMENSION(:)    ::  IGLBLISTSEC ! their corresponding sections
REAL   , ALLOCATABLE, DIMENSION(:)    ::  ZGLBLISTFLD ! their corresponding field average
!
!-------------------------------------------------------------------------------
!
!*       1.     INITIALISATION
!               --------------
!
INFO  = -1
IROOT =  0 
CALL MPI_COMM_RANK(NMNH_COMM_WORLD, IRANK, INFO)  ! get the rank of the current processor
!
IKE = SIZE(OCLOUD,3)
!
IF (OBOTTOMUP) THEN
  IKBEG=1       ! scans the volume from the surface
  IKEND=IKE     ! up to the top
  IKINC=1       ! increment to next level
ELSE
  IKBEG=IKE     ! scans the volume from the top
  IKEND=1       ! down to the surface
  IKINC=-1      ! increment to next level
END IF
!
DO JK=IKBEG,IKEND,IKINC
!PRINT *,IRANK,'initialisation, level =',JK
!
NULLIFY(TZFIELDS_ll)
!
CALL GET_GLOBALDIMS_ll (NIMAX_ll,NJMAX_ll)
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_OR_ll('B',IXOR,IYOR)
!
!-------------------------------------------------------------------------------
!
!*       2.     IDENTIFY CLUSTER PIXELS THANKS TO OCLOUD MASK
!               ---------------------------------------------
!PRINT *,IRANK,'identification'
!
!*--------------------------------------------------------------------*C
!*       2.1    identify clusters by uninterrupted zonal segments
!* 1st label in IMAPIND  -- & --  1st measure of size in ICLUSIZE
!*--------------------------------------------------------------------*C
!
ALLOCATE(IMAPIND(IIU,IJU))
IMAPIND(:,:) = JPUNDEF      ! The pixels first point toward no local indices.
!
ICLUSMAX=IIU*IJU
ALLOCATE(ICLUSIZE(ICLUSMAX))
ICLUSIZE(:) = 0             ! The local cluster size is initiated to 0.
!
ICPT = 0
!
DO JJ=IJB,IJE
  IF (OCLOUD(IIB,JJ,JK)) THEN
    ICPT = ICPT + 1         ! at begining of a row, pixel points toward a new local index
    IMAPIND(IIB,JJ) = ICPT
    ICLUSIZE(ICPT)   = ICLUSIZE(ICPT) + 1
  END IF
  DO JI=IIB+1,IIE
    IF (OCLOUD(JI,JJ,JK)) THEN
      IF (.NOT.OCLOUD(JI-1,JJ,JK)) THEN
        ICPT = ICPT + 1     ! keep previous local index only if previous pixel mask was also true
      END IF
      IMAPIND(JI,JJ) = ICPT
      ICLUSIZE(ICPT)   = ICLUSIZE(ICPT) + 1
    END IF
  END DO
END DO
!PRINT *,'IP=',IP,' ICPT=',ICPT,' COUNT(OCLOUD)=',COUNT(OCLOUD)
!IF (ICPT.GT.ICLUSMAX) PRINT *,"ICLUSMAX=",ICLUSMAX," LT  ICPT=",ICPT
!
! IMAPIND(JI,JJ) = JC means that (JI,JJ) pixel points toward JC as local index
! IMAPIND is JPUNDEF where OCLOUD is FALSE, it is in [1,ICLUSMAX] elsewhere
! IMAPIND have the same value along uninterrupted zonal rows
! IMAPIND will not change before step 4. 
! 
! ICLUSIZE(JC) is equal to the number of pixels pointing toward JC as local index
! ICLUSIZE will not change before step 3.
! sum(ICLUSIZE) is a constant and it is equal to num(OCLOUD) in the local domain
!
!*--------------------------------------------------------------------*C
!*       2.2    GATHER CONNECTED CLUSTERS
!*  create ICLUSIND table which will contain Local  Cluster Indices
!*   those Indices will further point toward Global Cluster ID in ICLUSTER...
!*--------------------------------------------------------------------*C
!
ALLOCATE(ICLUSIND(ICLUSMAX))
ICLUSIND = (/ (JC, JC=1,ICLUSMAX) /)
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
!  DO JI=IIB,IIE-1
    IF (OCLOUD(JI,JJ,JK)) THEN
!*--------------------------------------------------------------------*C
!* MINIMUM OF THE REGIONS THAT BORDER THIS ONE
      IMINI = IMAPIND(JI,JJ)
      IMINI = MIN(IMINI,MIN(IMAPIND(JI-1,JJ  ),IMAPIND(JI  ,JJ-1)))
!     IMINI = MIN(IMINI,MIN(IMAPIND(JI-1,JJ-1),IMAPIND(JI+1,JJ-1)))
!*--------------------------------------------------------------------*C
!* SEE WHETHER THERE ARE TWO CHAINS
!*--------------------------------------------------------------------*C
      IF (IMINI.GT.ICLUSIND(IMAPIND(JI,JJ))) THEN
!*--------------------------------------------------------------------*C
!* LOOK FOR THE MINIMUM OF THE CHAIN STARTING WITH KCLUSTER(I,J)
!*--------------------------------------------------------------------*C
        IMIX = ICLUSIND(IMAPIND(JI,JJ))
        DO
          IF ( ICLUSIND(IMIX) .NE. IMIX ) THEN
            IMIX = ICLUSIND(IMIX)
          ELSE
            EXIT
          END IF
        END DO
!*--------------------------------------------------------------------*C
!* LOOK FOR THE MINIMUM OF THE CHAIN STARTING WITH MINI
!*--------------------------------------------------------------------*C
        IMUX = ICLUSIND(IMINI)
        DO
          IF ( ICLUSIND(IMUX) .NE. IMUX ) THEN
            IMUX = ICLUSIND(IMUX)
          ELSE
            EXIT
          END IF
        END DO
!*--------------------------------------------------------------------*C
!* USE THE ABSOLUTE MINIMUM OF THE TWO CHAINS TO CONNECT THEM
!*--------------------------------------------------------------------*C
        IF (IMUX.GT.IMIX) THEN
          ICLUSIND(IMUX) = IMIX
        ELSEIF (IMIX.GT.IMUX) THEN
          ICLUSIND(IMIX) = IMUX
        END IF
!*--------------------------------------------------------------------*C
!* SINGLE CHAIN
!*--------------------------------------------------------------------*C
      ELSE
        ICLUSIND(IMAPIND(JI,JJ)) = IMINI
      END IF
    END IF
  END DO
END DO
!
! since now, each Local Index JC that belongs to a pixel in the mask
!    points (via ICLUSIND) toward itself -OR- toward a lower Local Index 
!                          that belongs to a connected pixel in the mask
!
!
!---------------------------------------------------------------------------------
!
!*       3.     OPTIMISE ICLUSIND TABLE -&- GATHER LOCAL CLUSTER SIZES IN ICLUSIZE
!               ------------------------------------------------------------------
!PRINT *,IRANK,'optimisation'
!
DO JC=1,ICLUSMAX
  ICLUSIND(JC) = ICLUSIND(ICLUSIND(JC))
  IF (ICLUSIND(JC).LT.JC) THEN
    ICLUSIZE(ICLUSIND(JC)) = ICLUSIZE(ICLUSIND(JC)) + ICLUSIZE(JC)
    ICLUSIZE(JC) = 0
  END IF
END DO
!
! since now, each Local Index JC that belong to a pixel in the mask
!    points (via ICLUSIND) toward the LOWEST Local Index
!        that belongs to a *connected* pixel in the mask
!
! ICLUSIZE(JC) gives now the number of pixels
!    that point (via IMAPIND and ICLUSIND) toward JC as Local Index
! sum(ICLUSIZE) is still equal to num(OCLOUD) in the local domain
!
!----------------------------------------------------------------------------------
!
!*       4.     GIVE TO CLUSTER PIXEL LOWEST LOCAL INDEX AND THEN GLOBAL CLUSTER ID
!               -------------------------------------------------------------------
!PRINT *,IRANK,'naming'
!
!*--------------------------------------------------------------------*C
!* 2nd label in IMAPIND  -- & --  1st label in in ICLUSTER... and IMAP..._ll
!*--------------------------------------------------------------------*C
!
ALLOCATE(ICLUSTERIDT(ICLUSMAX))
ALLOCATE(ICLUSTERLVL(ICLUSMAX))
ALLOCATE(IAFFLUSECTN(ICLUSMAX))
ALLOCATE(ZCLUSSUMFLD(ICLUSMAX))
ICLUSTERIDT(:)   = 0       ! local indices first indicate no global cluster ID
ICLUSTERLVL(:)   = IKE+1   ! and cluster base level is setup above highest level
IAFFLUSECTN(:)   = 0       ! no affluent from under so affluent section = 0
ZCLUSSUMFLD(:)   = 0.      ! sum of field values inside the clusters is set up to 0.
!
!        4.1    CHOOSING CLUSTER FUTURE NAMES & WRITING IT AND PROPERTIES IN ICLUSTER...
!              (UPDATE 1D TABLES)
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IF (OCLOUD(JI,JJ,JK)) THEN
    ! using optimised ICLUSIND to connect pixel to lowest local index pixel
    ! belonging to the same cluster in this horizontal level
      IMAPIND(JI,JJ) = ICLUSIND(IMAPIND(JI,JJ))           ! lowest local index
    ! adding field value of current pixel to its cluster sum
      ZCLUSSUMFLD(IMAPIND(JI,JJ)) = ZCLUSSUMFLD(IMAPIND(JI,JJ))+PFIELD(JI,JJ,JK)
    ! in case that pixel belongs to an already named cluster:
      ITEMPIDT = ICLUSTERIDT(IMAPIND(JI,JJ))                ! potential IDnumber
      ITEMPLVL = ICLUSTERLVL(IMAPIND(JI,JJ))                ! potential baselvl
      ITEMPAFL = IAFFLUSECTN(IMAPIND(JI,JJ))                ! potential affluent section
    ! checking whether there is an identified cluster at same position in previous level, an "affluent"
    ! (previous level is JK-IKINC)
    ! in this case, we settup future name in consequence
      IF (JK.NE.IKBEG) THEN
      ! 
      DO JDI=-1,1 ! cross-shape sieve of conseidered underlying pixels
      DO JDJ=-1,1 ! JI +/- 1 and JJ +/- 1 and (JI,JJ)
      IF ((ABS(JDI)+ABS(JDJ)).LT.2) THEN
      !
      IF (OCLOUD(JI+JDI,JJ+JDJ,JK-IKINC)) THEN
        ! 1st priority given to biggest affluent in section        
        IF (PCLUSTERSEC(JI+JDI,JJ+JDJ,JK-IKINC).GT.ITEMPAFL) THEN
          ITEMPIDT = KCLUSTERIDT(JI+JDI,JJ+JDJ,JK-IKINC)
          ITEMPLVL = KCLUSTERLVL(JI+JDI,JJ+JDJ,JK-IKINC)
          ITEMPAFL = PCLUSTERSEC(JI+JDI,JJ+JDJ,JK-IKINC)
        ELSEIF (PCLUSTERSEC(JI+JDI,JJ+JDJ,JK-IKINC).EQ.ITEMPAFL) THEN
          ! 2nd priority given to affluent with lowest base
          IF (KCLUSTERLVL(JI+JDI,JJ+JDJ,JK-IKINC).LT.ITEMPLVL) THEN
            ITEMPIDT = KCLUSTERIDT(JI+JDI,JJ+JDJ,JK-IKINC)
            ITEMPLVL = KCLUSTERLVL(JI+JDI,JJ+JDJ,JK-IKINC)
          ELSEIF (KCLUSTERLVL(JI+JDI,JJ+JDJ,JK-IKINC).EQ.ITEMPLVL) THEN
            ! 3rd priority given to affluent with lowest IDnumber
            IF (ITEMPIDT.EQ.0) THEN
              ITEMPIDT = KCLUSTERIDT(JI+JDI,JJ+JDJ,JK-IKINC)
            ELSE
              ITEMPIDT = MIN(KCLUSTERIDT(JI+JDI,JJ+JDJ,JK-IKINC),ITEMPIDT)
            END IF
          END IF
        END IF
      END IF
      !
      END IF
      END DO
      END DO
      !
      END IF
      !    
      IF (ITEMPIDT.EQ.0) THEN
        ICLUSTERIDT(IMAPIND(JI,JJ)) = NIMAX_ll*(IYOR+JJ-1) + (IXOR+JI)
        ICLUSTERLVL(IMAPIND(JI,JJ)) = JK
        IAFFLUSECTN(IMAPIND(JI,JJ)) = 0
      ELSE
        ICLUSTERIDT(IMAPIND(JI,JJ)) = ITEMPIDT
        ICLUSTERLVL(IMAPIND(JI,JJ)) = ITEMPLVL
        IAFFLUSECTN(IMAPIND(JI,JJ)) = ITEMPAFL
      END IF
      !
    END IF ! OCLOUD(JI,JJ) is TRUE
  END DO
END DO
!
DEALLOCATE(ICLUSIND)
!
! now, ICLUSTER... has a name for each cluster in the current level
! the same as the cluster under if there is one, a new name otherwise
! in case of several possibilities, we linked current level cluster prioritarly with:
! 1. cluster under with the bigest section
! 2. cluster under with the lowest base level
! 3. cluster under with the lowest ID number
!
!*       4.2    GIVING ID, BASE LEVEL and AFFLUENT SECTION in IMAP..._ll 
!              (COLORATION OF 2D TABLES) 
!
IMAPIDT_ll(:,:) = 0       ! pixels fist have no global cluster ID
IMAPLVL_ll(:,:) = IKE+1   ! and cluster base level is setup above highest level
IMAPAFL_ll(:,:) = 0       ! no affluent from under so affluent section = 0
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IF (OCLOUD(JI,JJ,JK)) THEN
      IMAPIDT_ll(JI,JJ) = ICLUSTERIDT(IMAPIND(JI,JJ)) ! global cluster ID
      IMAPLVL_ll(JI,JJ) = ICLUSTERLVL(IMAPIND(JI,JJ)) ! cluster baselevel
      IMAPAFL_ll(JI,JJ) = IAFFLUSECTN(IMAPIND(JI,JJ)) ! affluent section
    END IF
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       5.     UPDATE HALO AND UNIFY ADJACENT CLUSTERS GIVING COMMON ID*
!               --------------------------------------------------------
!               *ID is defined unequivocally by (IDT number, base LVL)
!
!              (MPI COMMUNICATIONS FOR 2D FIELDS:
!               IDT, LVL, AFL useful for clustering)
!
!PRINT *,IRANK,'halo-update'
!
DO JH=1,9999            ! update halo until there is no more changes anywhere
!PRINT *,'compteur premiere  boucle infinie =',JH
!
ZMAPIDT_ll=IMAPIDT_ll
ZMAPLVL_ll=IMAPLVL_ll
ZMAPAFL_ll=IMAPAFL_ll
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZMAPIDT_ll)
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZMAPLVL_ll)
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZMAPAFL_ll)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
IMAPIDT_ll=NINT(ZMAPIDT_ll)
IMAPLVL_ll=NINT(ZMAPLVL_ll)
IMAPAFL_ll=NINT(ZMAPAFL_ll)
ICHANGES = 0    ! counter of changes
!
DO JJ=IJB,IJE   ! check Western and Eastern borders of local domain
!
IF ( (IMAPIDT_ll(IIB,JJ).NE.0).AND.(IMAPIDT_ll(1,JJ).NE.0) ) THEN
  IF ( IMAPAFL_ll(1,JJ).GT.IAFFLUSECTN(IMAPIND(IIB,JJ)) ) THEN
  ! if halo pixel has a bigger affluent than current pixel
    ICLUSTERIDT(IMAPIND(IIB,JJ)) = IMAPIDT_ll(1,JJ)
    ICLUSTERLVL(IMAPIND(IIB,JJ)) = IMAPLVL_ll(1,JJ)
    IAFFLUSECTN(IMAPIND(IIB,JJ)) = IMAPAFL_ll(1,JJ)
    ICHANGES = ICHANGES + 1
  ELSEIF (IMAPAFL_ll(1,JJ).EQ.IAFFLUSECTN(IMAPIND(IIB,JJ))) THEN
    IF (IMAPLVL_ll(1,JJ).LT.ICLUSTERLVL(IMAPIND(IIB,JJ))) THEN
    ! if halo pixel belongs to a cluster with a lower base
      ICLUSTERIDT(IMAPIND(IIB,JJ)) = IMAPIDT_ll(1,JJ)
      ICLUSTERLVL(IMAPIND(IIB,JJ)) = IMAPLVL_ll(1,JJ)
      ICHANGES = ICHANGES + 1
    ELSEIF (IMAPLVL_ll(1,JJ).EQ.ICLUSTERLVL(IMAPIND(IIB,JJ))) THEN
      IF (IMAPIDT_ll(1,JJ).LT.ICLUSTERIDT(IMAPIND(IIB,JJ))) THEN
      ! if halo pixel belongs to a cluster with a lower ID number
      ICLUSTERIDT(IMAPIND(IIB,JJ)) = IMAPIDT_ll(1,JJ)
      ICHANGES = ICHANGES + 1
      END IF
    END IF
  END IF
END IF
!
IF ( (IMAPIDT_ll(IIE,JJ).NE.0).AND.(IMAPIDT_ll(IIU,JJ).NE.0) ) THEN
  IF ( IMAPAFL_ll(IIU,JJ).GT.IAFFLUSECTN(IMAPIND(IIE,JJ)) ) THEN
  ! if halo pixel has a bigger affluent than current pixel
    ICLUSTERIDT(IMAPIND(IIE,JJ)) = IMAPIDT_ll(IIU,JJ)
    ICLUSTERLVL(IMAPIND(IIE,JJ)) = IMAPLVL_ll(IIU,JJ)
    IAFFLUSECTN(IMAPIND(IIE,JJ)) = IMAPAFL_ll(IIU,JJ)
    ICHANGES = ICHANGES + 1
  ELSEIF (IMAPAFL_ll(IIU,JJ).EQ.IAFFLUSECTN(IMAPIND(IIE,JJ))) THEN
    IF (IMAPLVL_ll(IIU,JJ).LT.ICLUSTERLVL(IMAPIND(IIE,JJ))) THEN
    ! if halo pixel belongs to a cluster with a lower base
      ICLUSTERIDT(IMAPIND(IIE,JJ)) = IMAPIDT_ll(IIU,JJ)
      ICLUSTERLVL(IMAPIND(IIE,JJ)) = IMAPLVL_ll(IIU,JJ)
      ICHANGES = ICHANGES + 1
    ELSEIF (IMAPLVL_ll(IIU,JJ).EQ.ICLUSTERLVL(IMAPIND(IIE,JJ))) THEN
      IF (IMAPIDT_ll(IIU,JJ).LT.ICLUSTERIDT(IMAPIND(IIE,JJ))) THEN
      ! if halo pixel belongs to a cluster with a lower ID number
      ICLUSTERIDT(IMAPIND(IIE,JJ)) = IMAPIDT_ll(IIU,JJ)
      ICHANGES = ICHANGES + 1
      END IF
    END IF
  END IF
END IF
!
END DO          ! two borders checked
!
DO JI=IIB,IIE   ! check Southern and Northern borders of local domain
!
IF ( (IMAPIDT_ll(JI,IJB).NE.0).AND.(IMAPIDT_ll(JI,1).NE.0) ) THEN
  IF ( IMAPAFL_ll(JI,1).GT.IAFFLUSECTN(IMAPIND(JI,IJB)) ) THEN
  ! if halo pixel has a bigger affluent than current pixel
    ICLUSTERIDT(IMAPIND(JI,IJB)) = IMAPIDT_ll(JI,1)
    ICLUSTERLVL(IMAPIND(JI,IJB)) = IMAPLVL_ll(JI,1)
    IAFFLUSECTN(IMAPIND(JI,IJB)) = IMAPAFL_ll(JI,1)
    ICHANGES = ICHANGES + 1
  ELSEIF (IMAPAFL_ll(JI,1).EQ.IAFFLUSECTN(IMAPIND(JI,IJB))) THEN
    IF (IMAPLVL_ll(JI,1).LT.ICLUSTERLVL(IMAPIND(JI,IJB))) THEN
    ! if halo pixel belongs to a cluster with a lower base
      ICLUSTERIDT(IMAPIND(JI,IJB)) = IMAPIDT_ll(JI,1)
      ICLUSTERLVL(IMAPIND(JI,IJB)) = IMAPLVL_ll(JI,1)
      ICHANGES = ICHANGES + 1
    ELSEIF (IMAPLVL_ll(JI,1).EQ.ICLUSTERLVL(IMAPIND(JI,IJB))) THEN
      IF (IMAPIDT_ll(JI,1).LT.ICLUSTERIDT(IMAPIND(JI,IJB))) THEN
      ! if halo pixel belongs to a cluster with a lower ID number
      ICLUSTERIDT(IMAPIND(JI,IJB)) = IMAPIDT_ll(JI,1)
      ICHANGES = ICHANGES + 1
      END IF
    END IF
  END IF
END IF
!
IF ( (IMAPIDT_ll(JI,IJE).NE.0).AND.(IMAPIDT_ll(JI,IJU).NE.0) ) THEN
  IF ( IMAPAFL_ll(JI,IJU).GT.IAFFLUSECTN(IMAPIND(JI,IJE)) ) THEN
  ! if halo pixel has a bigger affluent than current pixel
    ICLUSTERIDT(IMAPIND(JI,IJE)) = IMAPIDT_ll(JI,IJU)
    ICLUSTERLVL(IMAPIND(JI,IJE)) = IMAPLVL_ll(JI,IJU)
    IAFFLUSECTN(IMAPIND(JI,IJE)) = IMAPAFL_ll(JI,IJU)
    ICHANGES = ICHANGES + 1
  ELSEIF (IMAPAFL_ll(JI,IJU).EQ.IAFFLUSECTN(IMAPIND(JI,IJE))) THEN
    IF (IMAPLVL_ll(JI,IJU).LT.ICLUSTERLVL(IMAPIND(JI,IJE))) THEN
    ! if halo pixel belongs to a cluster with a lower base
      ICLUSTERIDT(IMAPIND(JI,IJE)) = IMAPIDT_ll(JI,IJU)
      ICLUSTERLVL(IMAPIND(JI,IJE)) = IMAPLVL_ll(JI,IJU)
      ICHANGES = ICHANGES + 1
    ELSEIF (IMAPLVL_ll(JI,IJU).EQ.ICLUSTERLVL(IMAPIND(JI,IJE))) THEN
      IF (IMAPIDT_ll(JI,IJU).LT.ICLUSTERIDT(IMAPIND(JI,IJE))) THEN
      ! if halo pixel belongs to a cluster with a lower ID number
      ICLUSTERIDT(IMAPIND(JI,IJE)) = IMAPIDT_ll(JI,IJU)
      ICHANGES = ICHANGES + 1
      END IF
    END IF
  END IF
END IF
!
END DO          ! all 4 borders checked
!
IF (ICHANGES.GT.0) THEN  ! in case of changes due to halo, rename clusters
  DO JJ=IJB,IJE
  DO JI=IIB,IIE
  IF (OCLOUD(JI,JJ,JK)) THEN
     IMAPIDT_ll(JI,JJ) = ICLUSTERIDT(IMAPIND(JI,JJ))
     IMAPLVL_ll(JI,JJ) = ICLUSTERLVL(IMAPIND(JI,JJ))
     IMAPAFL_ll(JI,JJ) = IAFFLUSECTN(IMAPIND(JI,JJ))
  ELSE
     IMAPIDT_ll(JI,JJ) = 0
     IMAPLVL_ll(JI,JJ) = IKE+1
     IMAPAFL_ll(JI,JJ) = 0
  END IF
  END DO
  END DO
END IF
!
CALL REDUCESUM_ll(ICHANGES,IINFO_ll)
IF (ICHANGES.EQ.0) EXIT 
!
END DO ! infinite loop until there is no more changes anywhere
!
KCLUSTERIDT(:,:,JK) = IMAPIDT_ll
KCLUSTERLVL(:,:,JK) = IMAPLVL_ll
!
!
CALL CLEANLIST_ll(TZFIELDS_ll)
DEALLOCATE(IAFFLUSECTN)
!
!-------------------------------------------------------------------------------
!
!*       6.     COMPUTE GLOBAL CLUSTER SECTIONS AND FIELD AVERAGES
!               --------------------------------------------------
!
!              (MPI COMMUNICATIONS FOR 1D TABLES
!               useful for computing properties distributions
!               within cluster population)
!
!PRINT *,IRANK,'measuring'
!
!*       6.1    SUM UP LOCAL CLUSTER INFO IN 1D TABLES
!PRINT *,IRANK,'build non-redundant 1D table'
!
ICPR = COUNT(ICLUSTERIDT.NE.0) ! number of cluster parts in the local domain, may be redundant
!
ALLOCATE(ILOCLISTIDT(ICPR))  ! list of cluster IDs       in local domain  
ALLOCATE(ILOCLISTLVL(ICPR))  ! list of cluster baselevel in local domain  
ALLOCATE(ILOCLISTSEC(ICPR))  ! list of cluster sections  in local domain
ALLOCATE(ZLOCLISTFLD(ICPR))  ! list of cluster field sum in local domain
ALLOCATE(ILOCIND(ICLUSMAX))  ! pointer from old indices (1..iclusmax) toward new non-redundant tables
ILOCLISTIDT = 0
ILOCLISTLVL = 0
ILOCLISTSEC = 0
ZLOCLISTFLD = 0.
ILOCIND     = 0
!
! create ILOCLIST... and ILOCLISTSEC from ICLUSTER... and ICLUSIZE to get rid of all 0 elements
! and to sum the sizes of clusters appearing in different parts inside the local domain
! the sum of field values inside each cluster is continuously updated
!
ICPT=0                ! counter to know last index of ILOCLIST... that has been filled
IF (ICPR.GT.0) THEN
DO JC=1,ICLUSMAX
  IF (ICLUSTERIDT(JC).NE.0) THEN
    JJ=1        ! to check whether ILOCLIST... already contains ICLUSTER...(JC)
    DO
      IF (ILOCLISTIDT(JJ).EQ.ICLUSTERIDT(JC)) THEN
      IF (ILOCLISTLVL(JJ).EQ.ICLUSTERLVL(JC)) THEN
        ZLOCLISTFLD(JJ) = ZLOCLISTFLD(JJ) + ZCLUSSUMFLD(JC)
        ILOCLISTSEC(JJ) = ILOCLISTSEC(JJ) + ICLUSIZE(JC)
        ILOCIND(JC)     = JJ
        EXIT
      END IF
      END IF
      IF (JJ.EQ.ICPT+1) THEN
        ILOCLISTIDT(JJ) = ICLUSTERIDT(JC)
        ILOCLISTLVL(JJ) = ICLUSTERLVL(JC)
        ZLOCLISTFLD(JJ) = ZCLUSSUMFLD(JC)
        ILOCLISTSEC(JJ) = ICLUSIZE(JC)
        ILOCIND(JC)     = JJ
        ICPT = ICPT+1
        EXIT
      END IF
      JJ = JJ+1
    END DO
  END IF
END DO
END IF
!
DEALLOCATE(ZCLUSSUMFLD)
DEALLOCATE(ICLUSIZE)
DEALLOCATE(ICLUSTERIDT)
DEALLOCATE(ICLUSTERLVL)
!
! normally, ICPT = COUNT(ILOCLISTID.NE.0)
! number of cluster IDs in local domain, *non-redundant*
!
ALLOCATE(ILOCLISTIDT2(ICPT))
ALLOCATE(ILOCLISTLVL2(ICPT))
ALLOCATE(ILOCLISTSEC2(ICPT))
ALLOCATE(ZLOCLISTFLD2(ICPT))
ILOCLISTIDT2 = ILOCLISTIDT(1:ICPT)  ! JUST TO REMOVE ALL 0 AT THE END
ILOCLISTLVL2 = ILOCLISTLVL(1:ICPT)
ILOCLISTSEC2 = ILOCLISTSEC(1:ICPT)
ZLOCLISTFLD2 = ZLOCLISTFLD(1:ICPT)  ! /ILOCLISTSEC(1:ICPT) ! from sum compute average
DEALLOCATE(ILOCLISTIDT)
DEALLOCATE(ILOCLISTLVL)
DEALLOCATE(ILOCLISTSEC)
DEALLOCATE(ZLOCLISTFLD)
!
!*       6.2    CONCATENATE ALL LOCAL INFORMATION
!PRINT *,IRANK,'concatenate non-redundant 1D table'
!
ALLOCATE(ICLUSNBR(NPROC))   ! number of clusters in each proc's local domain
CALL MPI_ALLGATHER(ICPT, 1, MPI_INTEGER, ICLUSNBR, 1, MPI_INTEGER, NMNH_COMM_WORLD, INFO)
! each processor knows now how many clusters appear in all other processor domains
!
!PRINT *,IRANK,'build IPROCDPL'
ALLOCATE(IPROCDPL(NPROC))
IPROCDPL(1)=0
IF (NPROC.GT.1) THEN
DO JC=2,NPROC
  IPROCDPL(JC)=IPROCDPL(JC-1)+ICLUSNBR(JC-1)
END DO
END IF
! IPROCDPL(JC) = nbr of clusters contained by processors whose rank < JC (i.e. cumulative)
!
ITOTNBR = SUM(ICLUSNBR)  ! number of clusters ! with redundance from one local domain to the other !
!
!
ALLOCATE(IGLBLISTIDT(ITOTNBR)) ! concatenated list of cluster IDs        of all local domains
ALLOCATE(IGLBLISTLVL(ITOTNBR)) ! concatenated list of cluster baselevels of all local domains
ALLOCATE(IGLBLISTSEC(ITOTNBR)) ! their corresponding LOCAL sections   !!
ALLOCATE(ZGLBLISTFLD(ITOTNBR)) ! their corresponding LOCAL field sum
IGLBLISTIDT=0
IGLBLISTLVL=0
IGLBLISTSEC=0
ZGLBLISTFLD=0.
!
!PRINT *,IRANK,'call all-gatherv'
CALL MPI_ALLGATHERV(ILOCLISTIDT2, ICPT, MPI_INTEGER, IGLBLISTIDT, ICLUSNBR, IPROCDPL, MPI_INTEGER, &
                    NMNH_COMM_WORLD, INFO)
CALL MPI_ALLGATHERV(ILOCLISTLVL2, ICPT, MPI_INTEGER, IGLBLISTLVL, ICLUSNBR, IPROCDPL, MPI_INTEGER, &
                    NMNH_COMM_WORLD, INFO)
CALL MPI_ALLGATHERV(ILOCLISTSEC2, ICPT, MPI_INTEGER, IGLBLISTSEC, ICLUSNBR, IPROCDPL, MPI_INTEGER, &
                    NMNH_COMM_WORLD, INFO)
CALL MPI_ALLGATHERV(ZLOCLISTFLD2, ICPT, MPI_PRECISION, ZGLBLISTFLD, ICLUSNBR, IPROCDPL, MPI_PRECISION, &
                    NMNH_COMM_WORLD, INFO)
!
!*       6.3    EACH PROC COMPUTES GLOBAL SECTIONS AND FIELD AVERAGE OF ITS CLUSTERS
!PRINT *,IRANK,'global sizes of its own clusters' 
!
! update cluster sections in ILOCLISTSEC2 in order to have GLOBAL sections
DO JC=1,ICPT
  DO JJ=1,ITOTNBR
    ! TO BE SURE NOT TO SUM ITS OWN SIZE TO THE OTHER DOMAIN PART SIZES
    IF ( (JJ.LE.IPROCDPL(IRANK+1)) .OR. (JJ.GT.(IPROCDPL(IRANK+1)+ICLUSNBR(IRANK+1))) ) THEN
    IF ( (IGLBLISTIDT(JJ).EQ.ILOCLISTIDT2(JC)) .AND.  (IGLBLISTLVL(JJ).EQ.ILOCLISTLVL2(JC)) ) THEN
    ! if IGLBLIST...(JJ) correspond to ILOCLIST...2(JC) then sum sections and field avg
    ILOCLISTSEC2(JC)=ILOCLISTSEC2(JC)+IGLBLISTSEC(JJ)
    ZLOCLISTFLD2(JC)=ZLOCLISTFLD2(JC)+ZGLBLISTFLD(JJ)
    END IF
    END IF
  END DO
END DO
!
DO JC=1,ICPT
  ZLOCLISTFLD2(JC)=ZLOCLISTFLD2(JC)/ILOCLISTSEC2(JC) ! from sum toward average
END DO
!
!PRINT *,'ILOCLISTID2 = ',ILOCLISTID2
!PRINT *,'ILOCLISTSZ2 = ',ILOCLISTSZ2
!
! since now, ILOCLISTSEC2 contains global sizes     of local clusters
!       and, ZLOCLISTFLD2 contains global field avg of local clusters
!
!*       6.4    FILL PCLUSTERSEC TABLE WITH GLOBAL SIZES
!PRINT *,IRANK,'filling PCLUSTERSEC'
  !
  PCLUSTERSEC(:,:,JK) = 0. ! pixels first belong to no clusters, so clustersize is zero
  ZMESHAREA = XDXHATM*XDYHATM/10**6         ! in km2
    DO JJ=IJB,IJE          
    DO JI=IIB,IIE          
    IF (OCLOUD(JI,JJ,JK)) THEN
    PCLUSTERSEC(JI,JJ,JK) = ILOCLISTSEC2(ILOCIND(IMAPIND(JI,JJ)))*ZMESHAREA
    END IF
    END DO
    END DO
!
!* DEALLOCATIONS
!
DEALLOCATE(ILOCLISTIDT2)
DEALLOCATE(ILOCLISTLVL2)
DEALLOCATE(ILOCLISTSEC2)
DEALLOCATE(ZLOCLISTFLD2)
!
DEALLOCATE(ILOCIND,IMAPIND)
DEALLOCATE(ICLUSNBR,IPROCDPL)
DEALLOCATE(IGLBLISTIDT,IGLBLISTLVL,IGLBLISTSEC,ZGLBLISTFLD)
!
END DO ! JK
!
!
END SUBROUTINE CLUSTERING
