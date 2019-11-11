!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_MPPDB
!
!       Modifs :
!!      J.Escobar 23/10/2012: correct CHECK_LB & format print output 
!!      M.Moge 05/02/2015: MPPDB_CHECK_SURFEX2D and MPPDB_CHECK_SURFEX3D + bug fix in MPPDB_CHECK2D and MPPDB_CHECK3D (call MPI_AllReduce at the beginning)
!  J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!  G.Delautier : 23/06/2016 : surfex v8
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!  Philippe Wautelet: 22/01/2019: use standard FLUSH statement instead of non standard intrinsics
!  Philippe Wautelet: 22/01/2019: use sleep_c subroutine instead of non-standard call system
!
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use modi_tools_c

  IMPLICIT NONE


  INTEGER            ,PARAMETER    :: NSLEEP=30 !Sleep duration to wait for (in seconds)
  INTEGER            ,PARAMETER    :: chlg=256
  CHARACTER(LEN=chlg),PARAMETER    :: MPPDB_CONF = "mppdb.nam"

  LOGICAL                          :: MPPDB_INITIALIZED = .FALSE.
  LOGICAL                          :: MPPDB_DEBUG = .FALSE.

  CHARACTER(LEN=chlg)              :: MPPDB_EXEC  = "PREP_REAL_CASE" 
  CHARACTER(LEN=chlg)              :: MPPDB_HOST  = "localhost" 
  CHARACTER(LEN=chlg)              :: MPPDB_WDIR  = "."
  INTEGER                          :: MPPDB_NBSON = 1

  INTEGER                          :: MPPDB_INTER_COMM,MPPDB_INTRA_COMM
  INTEGER                          :: MPPDB_IRANK_WORLD,MPPDB_IRANK_INTRA
  INTEGER                          :: MPPDB_NBPROC_WORLD,MPPDB_NBPROC_INTRA

  LOGICAL                          :: MPPDB_FATHER_WORLD = .FALSE.

  REAL                             :: PRECISION = 1e-8 * 0.0
  LOGICAL                          :: MPPDB_CHECK_LB = .FALSE.

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_INIT()
#ifdef MNH_SP4
    !pas de mpi_spawn sur IBM-SP ni MPI_ARGV_NULL etc ...
    RETURN           
#else
    USE MODD_MPIF
    !JUANZ
    USE  MODE_MNH_WORLD , ONLY :  INIT_NMNH_COMM_WORLD
    USE  MODD_VAR_ll    , ONLY :  NMNH_COMM_WORLD
    !JUANZ
    IMPLICIT NONE


    INTEGER                         :: IUNIT
    INTEGER                         :: IERR
    INTEGER                         :: IRANK_WORLD,IRANK_INTRA
    INTEGER                         :: NBPROC_WORLD,NBPROC_INTRA

    LOGICAL                         :: GISINIT
    LOGICAL                         :: drapeau

    INTEGER                         :: INFO_SPAWN
    INTEGER                         :: RANK_FATHER = 0
    INTEGER,ALLOCATABLE             :: info_error(:)
    CHARACTER(LEN=40)               :: chaine
    LOGICAL                         :: isset



    NAMELIST /NAM_MPPDB/ MPPDB_DEBUG,MPPDB_EXEC,MPPDB_HOST,MPPDB_NBSON,MPPDB_WDIR,MPPDB_CHECK_LB

    !NMNH_COMM_WORLD = MPI_COMM_WORLD

    ! If already initialized, nothing to do
    IF (MPPDB_INITIALIZED) RETURN
    !
    MPPDB_INITIALIZED = .TRUE.
    !
    ! Init MPI
    !
    CALL MPI_INITIALIZED(GISINIT, ierr)
    IF (.NOT. GISINIT) THEN
       !CALL MPI_INIT(ierr) 
       CALL  INIT_NMNH_COMM_WORLD(ierr)
    END IF
    !
    ! Get me rank in the my world
    !
    CALL MPI_COMM_RANK(NMNH_COMM_WORLD, MPPDB_IRANK_WORLD , ierr)
    CALL MPI_COMM_SIZE(NMNH_COMM_WORLD, MPPDB_NBPROC_WORLD, ierr)
    !
    !... Have I a father ?
    !
    CALL MPI_COMM_GET_PARENT(MPPDB_INTER_COMM, ierr)
    IF (MPPDB_INTER_COMM == MPI_COMM_NULL) THEN 
       !
       ! NO Father !
       MPPDB_FATHER_WORLD = .TRUE.  
       CALL MPI_BARRIER(NMNH_COMM_WORLD,ierr) 
       !
       ! if no config file , inactive MPPDB routines
       !
       OPEN(newunit=IUNIT,file=MPPDB_CONF,STATUS='OLD',iostat=IERR)
       IF (IERR.NE.0 ) THEN
          ! PRINT*,"IOSTAT=",IERR
          !
          !no config file => exit from MPPDB not actived
          !
          MPPDB_INITIALIZED = .FALSE.
          RETURN          
       END IF
       !
       ! read the parameter
       !
       READ(unit=IUNIT,NML=NAM_MPPDB)
       CLOSE(unit=IUNIT)
       IF ( .NOT. MPPDB_DEBUG ) THEN
          ! why don't when MMPDB to work
          MPPDB_INITIALIZED = .FALSE.
          RETURN
       ENDIF
       !
       IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          !-------------------------------------------------------------------------!
          !                                                                         !
          ! NO Father & rank=0 <=> I'm the root Father so Init all and Clone        !
          !                                                                         !
          !-------------------------------------------------------------------------!
          IF (MPPDB_DEBUG) THEN
             WRITE(*,NML=NAM_MPPDB)
          ENDIF
          !
          ! Create the info contecte for the son
          !
          ! host
          !
          CALL MPI_INFO_CREATE (INFO_SPAWN , ierr)
          !CALL MPI_INFO_SET    (INFO_SPAWN , "host", MPPDB_HOST , ierr)
          !CALL MPI_INFO_GET    (INFO_SPAWN , "host", 40, chaine, isset ,ierr)
          !IF (MPPDB_DEBUG) PRINT*,"MPPDB_INIT:: FATHER ::INFO_SPAWN , host=",isset,chaine
          !IF (ierr.NE.0) STOP 'MPPDB_INIT :: PB MPI_INFO_SET "host" '
          !
          ! working directory
          !
          CALL MPI_INFO_SET    (INFO_SPAWN , "wdir", MPPDB_WDIR , ierr)
          CALL MPI_INFO_GET    (INFO_SPAWN , "wdir", 40, chaine, isset ,ierr)
          IF (MPPDB_DEBUG) PRINT*,"MPPDB_INIT:: FATHER :: INFO_SPAWN , wdir=",isset,chaine
          IF (ierr.NE.0) STOP 'MPPDB_INIT:: PB MPI_INFO_SET  "wdir" '
          !
       ELSE
          ! other father only do nothing but participate 
          INFO_SPAWN =  MPI_INFO_NULL
          !
       END IF !  IF (MPPDB_IRANK_WORLD.EQ.0)
       !
       ! clone the son
       !
       ALLOCATE(info_error(MPPDB_NBSON))
       CALL MPI_BARRIER(NMNH_COMM_WORLD,ierr) 
       !
       CALL MPI_COMM_SPAWN(MPPDB_EXEC, MPI_ARGV_NULL,MPPDB_NBSON,INFO_SPAWN, &
            RANK_FATHER, NMNH_COMM_WORLD,MPPDB_INTER_COMM ,info_error, ierr)
       !
       CALL MPI_BARRIER(NMNH_COMM_WORLD,ierr) 
       !
       DEALLOCATE(info_error)
       !
       ! merge the communicator
       !
       drapeau=.FALSE.
       CALL MPI_INTERCOMM_MERGE(MPPDB_INTER_COMM, drapeau, MPPDB_INTRA_COMM, ierr)
       !
       !... Numbre of processus in MPPDB_INTRA_COMM.
       CALL MPI_COMM_SIZE(MPPDB_INTRA_COMM, mppdb_nbproc_intra, ierr)
       !         
       !... My rank in MPPDB_INTRA_COMM
       CALL MPI_COMM_RANK(MPPDB_INTRA_COMM, mppdb_irank_intra, ierr)
       IF (MPPDB_IRANK_WORLD.EQ.0) THEN 
          !  I'm the first father 
       IF (MPPDB_DEBUG) print*,"MPPDB_INIT :: FIRST FATHER mppdb_irank_intra=", mppdb_irank_intra &
            ,"mppdb_nbproc_intra=",mppdb_nbproc_intra
       flush(unit=OUTPUT_UNIT)
       endif
       !
       ! Wait the sons
       !
       CALL MPI_BARRIER  ( MPPDB_INTRA_COMM , ierr )
       ! WAIT FOR TOTALVIEW IF NEEDED 
       call sleep_c(NSLEEP)
       !
    ELSE ! (MPPDB_INTER_COMM <> MPI_COMM_NULL)
       !-------------------------------------------------------------------------!
       !                                                                         !
       !  I've a father <=> I'm a son                                            !
       !                                                                         !
       !-------------------------------------------------------------------------!
       !
       CALL MPI_BARRIER(NMNH_COMM_WORLD,ierr) 
       !
       ! merge the communicator
       !
       drapeau=.TRUE.
       CALL MPI_INTERCOMM_MERGE(MPPDB_INTER_COMM, drapeau, MPPDB_INTRA_COMM, ierr)
       !
       !... Numbre of processus in MPPDB_INTRA_COMM.
       CALL MPI_COMM_SIZE(MPPDB_INTRA_COMM, mppdb_nbproc_intra, ierr)
       !
       !... My rang in MPPDB_INTRA_COMM
       CALL MPI_COMM_RANK(MPPDB_INTRA_COMM, mppdb_irank_intra, ierr)
       !
       !
       ! Wait the FATHER's
       !
       CALL MPI_BARRIER  ( MPPDB_INTRA_COMM , ierr )
       !
       ! WAIT FOR TOTALVIEW IF NEEDED 
       call sleep_c(NSLEEP)
       !
       MPPDB_DEBUG = .TRUE.
       IF (MPPDB_DEBUG) write(200,*) "MPPDB_INIT :: FIRST SON mppdb_irank_intra=", mppdb_irank_intra &
            ,"MPPDB_IRANK_WORLD=",MPPDB_IRANK_WORLD
       !
       IF (MPPDB_IRANK_WORLD.EQ.0) THEN 
          !  I'm the first son 
          MPPDB_DEBUG = .TRUE.
          IF (MPPDB_DEBUG) print*,"MPPDB_INIT :: FIRSTSON mppdb_irank_intra=", mppdb_irank_intra
       END IF
    END IF !  IF (MPPDB_INTER_COMM == MPI_COMM_NULL)

    CALL MPI_BARRIER  ( MPPDB_INTRA_COMM , ierr )

#endif
  END SUBROUTINE MPPDB_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_BARRIER()
#ifdef MNH_SP4
    !pas de mpi_spawn sur IBM-SP ni MPI_ARGV_NULL etc ...
    RETURN           
#else
    IMPLICIT NONE
    INTEGER      :: IERR
    !
    ! synchronize all father & sons
    !
    IF ( .NOT. MPPDB_INITIALIZED ) RETURN 
    !
    CALL MPI_BARRIER(MPPDB_INTRA_COMM,IERR)
    !
#endif
  END SUBROUTINE MPPDB_BARRIER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_CHECK3D(PTAB,MESSAGE,PRECISION)

    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_VAR_ll    , ONLY : MPI_PRECISION
    USE MODD_MPIF      , ONLY : MPI_INTEGER, MPI_STATUS_IGNORE, MPI_SUM
    USE MODE_GATHER_ll

    IMPLICIT NONE


    REAL, DIMENSION(:,:,:)               :: PTAB
    CHARACTER(len=*)                     :: MESSAGE
    REAL                                 :: PRECISION

    !
    ! local var
    !
    REAL,ALLOCATABLE,TARGET, DIMENSION(:,:,:)   :: TAB_ll,TAB_SON_ll,TAB_SAVE_ll
    INTEGER                              :: IIMAX_ll,IJMAX_ll
    INTEGER                              :: IIU,IJU,IIU_ll,IJU_ll,IKU_ll
    INTEGER                              :: IINFO_ll

    INTEGER,PARAMETER                    :: ITAG1  = 12345 , ITAG2 = 123456

    INTEGER                              :: I_FIRST_SON
    INTEGER                              :: I_FIRST_FATHER
    REAL                                 :: MAX_DIFF , MAX_VAL
    INTEGER                              :: IIB_ll,IIE_ll,IJB_ll,IJE_ll
    INTEGER                              :: IGLBSIZEPTAB

    REAL,POINTER, DIMENSION(:,:,:)   :: TAB_INTERIOR_ll ! for easy debug
    INTEGER                              :: IK
    INTEGER                              :: KSIZEBUF

    INTEGER                              :: IIU_SON_ll,IJU_SON_ll,IKU_SON_ll
    INTEGER                              :: IIB_SON_ll,IIE_SON_ll,IJB_SON_ll,IJE_SON_ll
    INTEGER                              :: IHEXT_SON_ll , IDIFF_HEXT

#ifdef MNH_SP4
    !pas de mpi_spawn sur IBM-SP ni MPI_ARGV_NULL etc ...
    RETURN           
#else
    IF ( ( .NOT. MPPDB_INITIALIZED ) ) RETURN
    !get the global size of PTAB
    CALL MPI_ALLREDUCE(SIZE(PTAB), IGLBSIZEPTAB, 1,MPI_INTEGER, MPI_SUM, MPPDB_INTER_COMM, IINFO_ll)
    IF ( IGLBSIZEPTAB == 0 ) RETURN
    !
    CALL MPPDB_BARRIER()
    !
    IF(MPPDB_FATHER_WORLD) THEN
       !
       ! Reconstruct the whole PTAB in TAB_ll
       !
       CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
       IIU_ll = IIMAX_ll+2*JPHEXT
       IJU_ll = IJMAX_ll+2*JPHEXT
       IKU_ll = SIZE(PTAB,3)
       ALLOCATE(TAB_ll(IIU_ll,IJU_ll,IKU_ll))
       ALLOCATE(TAB_SAVE_ll(IIU_ll,IJU_ll,IKU_ll))
       CALL GATHERALL_FIELD_ll('XY',PTAB,TAB_ll,IINFO_ll)

       IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          !
          ! I'm the first FATHER => recieve the correct globale ARRAY from first son
          !
          !
          ! the first son , is the next processus after this 'world' so
          !
          I_FIRST_SON = MPPDB_NBPROC_WORLD         
          !
          ! recieve JPHEXT from son if different
          !
          CALL MPI_RECV(IHEXT_SON_ll,1,MPI_INTEGER,I_FIRST_SON, &
               ITAG1, MPPDB_INTRA_COMM,MPI_STATUS_IGNORE, IINFO_ll)

          !IHEXT_SON_ll = JPHEXT

          IIU_SON_ll = IIMAX_ll+2*IHEXT_SON_ll
          IJU_SON_ll = IJMAX_ll+2*IHEXT_SON_ll
          IKU_SON_ll = SIZE(PTAB,3)
    
          ALLOCATE(TAB_SON_ll(IIU_SON_ll,IJU_SON_ll,IKU_SON_ll))
          !
          CALL MPI_RECV(TAB_SON_ll,SIZE(TAB_SON_ll),MPI_PRECISION,I_FIRST_SON, &
               ITAG2, MPPDB_INTRA_COMM,MPI_STATUS_IGNORE, IINFO_ll)
          !

          !

          IF (MPPDB_CHECK_LB) THEN
             IDIFF_HEXT = MIN(JPHEXT,IHEXT_SON_ll)
          ELSE
             IDIFF_HEXT = 0
          ENDIF
          IIB_ll   = 1 + JPHEXT    ; IJB_ll = 1 + JPHEXT
          IIE_ll   = IIU_ll-JPHEXT ; IJE_ll = IJU_ll-JPHEXT
          
          IIB_SON_ll   = 1 + IHEXT_SON_ll    ; IJB_SON_ll = 1 + IHEXT_SON_ll
          IIE_SON_ll   = IIU_SON_ll-IHEXT_SON_ll ; IJE_SON_ll = IJU_SON_ll-IHEXT_SON_ll
          !
          TAB_SAVE_ll = TAB_ll
          TAB_ll      = 0.0
          TAB_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll)   &
            = ABS ( TAB_SAVE_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll) & 
            -       TAB_SON_ll(IIB_SON_ll-IDIFF_HEXT:IIE_SON_ll+IDIFF_HEXT,IJB_SON_ll-IDIFF_HEXT:IJE_SON_ll+IDIFF_HEXT &
                              ,1:IKU_SON_ll) )

          MAX_VAL  = MAXVAL( ABS (TAB_SON_ll(IIB_SON_ll-IDIFF_HEXT:IIE_SON_ll+IDIFF_HEXT,&
                                             IJB_SON_ll-IDIFF_HEXT:IJE_SON_ll+IDIFF_HEXT,1:IKU_SON_ll) ) )
          IF ( MAX_VAL .EQ. 0.0 ) MAX_VAL = 1.0
          MAX_DIFF=MAXVAL(TAB_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll)/MAX_VAL)
          TAB_INTERIOR_ll=> TAB_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll)
          !
          IF (MAX_DIFF > PRECISION ) THEN
             write(6, '(" MPPDB_CHECK3D :: PB MPPDB_CHECK3D =",A40," ERROR=",e15.8," MAXVAL=",e15.8)' ) MESSAGE,MAX_DIFF , MAX_VAL
          ELSE
             write(6, '(" MPPDB_CHECK3D :: OK MPPDB_CHECK3D =",A40," ERROR=",e15.8," MAXVAL=",e15.8)' ) MESSAGE,MAX_DIFF , MAX_VAL
          END IF
          flush(unit=OUTPUT_UNIT)
          !
          DEALLOCATE(TAB_ll,TAB_SON_ll)
          !
       END IF
    ELSE
       !
       ! Reconstruct the all PTAB in TAB_ll
       !
       CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
       IIU_ll = IIMAX_ll+2*JPHEXT
       IJU_ll = IJMAX_ll+2*JPHEXT
       IKU_ll = SIZE(PTAB,3)
       ALLOCATE(TAB_ll(IIU_ll,IJU_ll,IKU_ll))
       CALL GATHERALL_FIELD_ll('XY',PTAB,TAB_ll,IINFO_ll)
       !
       ! SON WORLD 
       !
       IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          !
          ! first son --> send the good array to the first father
          !
          I_FIRST_FATHER = 0
          IHEXT_SON_ll = JPHEXT
          CALL MPI_BSEND(IHEXT_SON_ll,1,MPI_INTEGER,I_FIRST_FATHER, &
               ITAG1, MPPDB_INTRA_COMM, IINFO_ll)

          CALL MPI_BSEND(TAB_ll,SIZE(TAB_ll),MPI_PRECISION,I_FIRST_FATHER, &
               ITAG2, MPPDB_INTRA_COMM, IINFO_ll)
 
       END IF
    END IF

    CALL MPPDB_BARRIER()
#endif
  END SUBROUTINE MPPDB_CHECK3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_CHECK3DM(MESSAGE,PRECISION &
                           ,PTAB1,PTAB2,PTAB3,PTAB4,PTAB5,PTAB6,PTAB7,PTAB8,PTAB9,PTAB10 &
                           ,PTAB11,PTAB12,PTAB13,PTAB14,PTAB15,PTAB16,PTAB17,PTAB18,PTAB19,PTAB20 &
                           )

    IMPLICIT NONE

    CHARACTER(len=*)                     :: MESSAGE
    REAL                                 :: PRECISION
    REAL, DIMENSION(:,:,:), OPTIONAL     :: PTAB1,PTAB2,PTAB3,PTAB4,PTAB5,PTAB6,PTAB7,PTAB8,PTAB9,PTAB10
    REAL, DIMENSION(:,:,:), OPTIONAL     :: PTAB11,PTAB12,PTAB13,PTAB14,PTAB15,PTAB16,PTAB17,PTAB18,PTAB19,PTAB20

    IF (PRESENT(PTAB1)) CALL MPPDB_CHECK3D(PTAB1,MESSAGE//"::PTAB1",PRECISION)
    IF (PRESENT(PTAB2)) CALL MPPDB_CHECK3D(PTAB2,MESSAGE//"::PTAB2",PRECISION)
    IF (PRESENT(PTAB3)) CALL MPPDB_CHECK3D(PTAB3,MESSAGE//"::PTAB3",PRECISION)
    IF (PRESENT(PTAB4)) CALL MPPDB_CHECK3D(PTAB4,MESSAGE//"::PTAB4",PRECISION)
    IF (PRESENT(PTAB5)) CALL MPPDB_CHECK3D(PTAB5,MESSAGE//"::PTAB5",PRECISION)
    IF (PRESENT(PTAB6)) CALL MPPDB_CHECK3D(PTAB6,MESSAGE//"::PTAB6",PRECISION)
    IF (PRESENT(PTAB7)) CALL MPPDB_CHECK3D(PTAB7,MESSAGE//"::PTAB7",PRECISION)
    IF (PRESENT(PTAB8)) CALL MPPDB_CHECK3D(PTAB8,MESSAGE//"::PTAB8",PRECISION)
    IF (PRESENT(PTAB9)) CALL MPPDB_CHECK3D(PTAB9,MESSAGE//"::PTAB9",PRECISION)
    IF (PRESENT(PTAB10)) CALL MPPDB_CHECK3D(PTAB10,MESSAGE//"::PTAB10",PRECISION)
    IF (PRESENT(PTAB11)) CALL MPPDB_CHECK3D(PTAB11,MESSAGE//"::PTAB11",PRECISION)
    IF (PRESENT(PTAB12)) CALL MPPDB_CHECK3D(PTAB12,MESSAGE//"::PTAB12",PRECISION)
    IF (PRESENT(PTAB13)) CALL MPPDB_CHECK3D(PTAB13,MESSAGE//"::PTAB13",PRECISION)
    IF (PRESENT(PTAB14)) CALL MPPDB_CHECK3D(PTAB14,MESSAGE//"::PTAB14",PRECISION)
    IF (PRESENT(PTAB15)) CALL MPPDB_CHECK3D(PTAB15,MESSAGE//"::PTAB15",PRECISION)
    IF (PRESENT(PTAB16)) CALL MPPDB_CHECK3D(PTAB16,MESSAGE//"::PTAB16",PRECISION)
    IF (PRESENT(PTAB17)) CALL MPPDB_CHECK3D(PTAB17,MESSAGE//"::PTAB17",PRECISION)
    IF (PRESENT(PTAB18)) CALL MPPDB_CHECK3D(PTAB18,MESSAGE//"::PTAB18",PRECISION)
    IF (PRESENT(PTAB19)) CALL MPPDB_CHECK3D(PTAB19,MESSAGE//"::PTAB19",PRECISION)
    IF (PRESENT(PTAB20)) CALL MPPDB_CHECK3D(PTAB20,MESSAGE//"::PTAB20",PRECISION)

  END SUBROUTINE MPPDB_CHECK3DM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_CHECK2D(PTAB,MESSAGE,PRECISION)

    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODE_GATHER_ll
    USE MODD_VAR_ll    , ONLY : MPI_PRECISION
    USE MODD_MPIF      , ONLY : MPI_INTEGER, MPI_STATUS_IGNORE, MPI_SUM

    USE  MODD_VAR_ll    , ONLY :  NMNH_COMM_WORLD

    IMPLICIT NONE

    REAL, DIMENSION(:,:)                 :: PTAB
    CHARACTER(len=*)                     :: MESSAGE
    REAL                                 :: PRECISION

    !
    ! local var
    !
    REAL,ALLOCATABLE,TARGET, DIMENSION(:,:)     :: TAB_ll,TAB_SON_ll,TAB_SAVE_ll
    INTEGER                              :: IIMAX_ll,IJMAX_ll
    INTEGER                              :: IIU,IJU,IIU_ll,IJU_ll
    INTEGER                              :: IINFO_ll

    INTEGER,PARAMETER                    :: ITAG = 12345

    INTEGER                              :: I_FIRST_SON
    INTEGER                              :: I_FIRST_FATHER
    REAL                                 :: MAX_DIFF , MAX_VAL
    INTEGER                              :: IIB_ll,IIE_ll,IJB_ll,IJE_ll

    REAL,POINTER, DIMENSION(:,:)   :: TAB_INTERIOR_ll ! for easy debug
    INTEGER                              :: IGLBSIZEPTAB

    INTEGER                              :: IIU_SON_ll,IJU_SON_ll
    INTEGER                              :: IIB_SON_ll,IIE_SON_ll,IJB_SON_ll,IJE_SON_ll
    INTEGER                              :: IHEXT_SON_ll , IDIFF_HEXT

#ifdef MNH_SP4
    !pas de mpi_spawn sur IBM-SP ni MPI_ARGV_NULL etc ...
    RETURN           
#else
    IF ( ( .NOT. MPPDB_INITIALIZED ) ) RETURN
    CALL MPI_ALLREDUCE(SIZE(PTAB), IGLBSIZEPTAB, 1,MPI_INTEGER, MPI_SUM, MPPDB_INTRA_COMM, IINFO_ll)
    IF ( IGLBSIZEPTAB == 0 ) RETURN

    CALL MPPDB_BARRIER()

    IF(MPPDB_FATHER_WORLD) THEN
       !
       ! Reconstruct the all PTAB in TAB_ll
       !
       CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
       IIU_ll = IIMAX_ll+2*JPHEXT
       IJU_ll = IJMAX_ll+2*JPHEXT
       ALLOCATE(TAB_ll(IIU_ll,IJU_ll))
       ALLOCATE(TAB_SAVE_ll(IIU_ll,IJU_ll))
       CALL GATHERALL_FIELD_ll('XY',PTAB,TAB_ll,IINFO_ll)

       IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          !
          ! I'm the first FATHER => recieve the correct globale ARRAY from first son
          !
          !
          ! the first son , is the next processus after this 'world' so
          !
          I_FIRST_SON = MPPDB_NBPROC_WORLD
          !
          ! recieve JPHEXT from son if different
          !
          CALL MPI_RECV(IHEXT_SON_ll,1,MPI_INTEGER,I_FIRST_SON, &
               ITAG, MPPDB_INTRA_COMM,MPI_STATUS_IGNORE, IINFO_ll)

          IIU_SON_ll = IIMAX_ll+2*IHEXT_SON_ll
          IJU_SON_ll = IJMAX_ll+2*IHEXT_SON_ll
          
          ALLOCATE(TAB_SON_ll(IIU_SON_ll,IJU_SON_ll))

          !
          CALL MPI_RECV(TAB_SON_ll,SIZE(TAB_SON_ll),MPI_PRECISION,I_FIRST_SON, &
               ITAG, MPPDB_INTRA_COMM,MPI_STATUS_IGNORE, IINFO_ll)
          !
 
          IF (MPPDB_CHECK_LB) THEN
             IDIFF_HEXT = MIN(JPHEXT,IHEXT_SON_ll)
          ELSE
             IDIFF_HEXT = 0
          ENDIF

          IIB_ll   = 1 + JPHEXT    ; IJB_ll = 1 + JPHEXT
          IIE_ll   = IIU_ll-JPHEXT ; IJE_ll = IJU_ll-JPHEXT
          IIB_SON_ll   = 1 + IHEXT_SON_ll    ; IJB_SON_ll = 1 + IHEXT_SON_ll
          IIE_SON_ll   = IIU_SON_ll-IHEXT_SON_ll ; IJE_SON_ll = IJU_SON_ll-IHEXT_SON_ll
       
          !
          TAB_SAVE_ll = TAB_ll
          TAB_ll      = 0.0
          TAB_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT)   &
            = ABS ( TAB_SAVE_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT) & 
            -       TAB_SON_ll(IIB_SON_ll-IDIFF_HEXT:IIE_SON_ll+IDIFF_HEXT,IJB_SON_ll-IDIFF_HEXT:IJE_SON_ll+IDIFF_HEXT) )         

          MAX_VAL  = MAXVAL( ABS (TAB_SON_ll(IIB_SON_ll-IDIFF_HEXT:IIE_SON_ll+IDIFF_HEXT, &
                                             IJB_SON_ll-IDIFF_HEXT:IJE_SON_ll+IDIFF_HEXT) ) ) 
          IF ( MAX_VAL .EQ. 0.0 ) MAX_VAL = 1.0
          MAX_DIFF = MAXVAL( TAB_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IIB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT) / MAX_VAL )
          TAB_INTERIOR_ll => TAB_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IIB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT)
          IF (MAX_DIFF > PRECISION ) THEN
             write(6, '(" MPPDB_CHECK2D :: PB MPPDB_CHECK2D =",A40," ERROR=",e15.8," MAXVAL=",e15.8)' ) MESSAGE,MAX_DIFF , MAX_VAL
          ELSE
             write(6, '(" MPPDB_CHECK2D :: OK MPPDB_CHECK2D =",A40," ERROR=",e15.8," MAXVAL=",e15.8)' ) MESSAGE,MAX_DIFF , MAX_VAL
          END IF
          flush(unit=OUTPUT_UNIT)
          !
          DEALLOCATE(TAB_ll,TAB_SON_ll)
          !
       END IF
    ELSE
       !
       ! SON WORLD 
       !
       !
       ! Reconstruct the all PTAB in TAB_ll
       !
       CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
       IIU_ll = IIMAX_ll+2*JPHEXT
       IJU_ll = IJMAX_ll+2*JPHEXT
       ALLOCATE(TAB_ll(IIU_ll,IJU_ll))
       CALL GATHERALL_FIELD_ll('XY',PTAB,TAB_ll,IINFO_ll)

       IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          !
          ! first son --> send the good array to the first father
          !
          I_FIRST_FATHER = 0
          CALL MPI_BSEND(JPHEXT,1,MPI_INTEGER,I_FIRST_FATHER, &
               ITAG, MPPDB_INTRA_COMM, IINFO_ll)
          CALL MPI_BSEND(TAB_ll,SIZE(TAB_ll),MPI_PRECISION,I_FIRST_FATHER, &
               ITAG, MPPDB_INTRA_COMM, IINFO_ll)
       END IF
    END IF

    CALL MPPDB_BARRIER()

#endif

  END SUBROUTINE MPPDB_CHECK2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_CHECKLB(PLB,MESSAGE,PRECISION,HLBTYPE,KRIM)

    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    USE MODD_VAR_ll    , ONLY : MPI_PRECISION ,  NMNH_COMM_WORLD
    USE MODD_IO_ll,        ONLY : ISP,ISNPROC,GSMONOPROC,LPACK,L2D
    USE MODD_MPIF      , ONLY   : MPI_INTEGER, MPI_STATUS_IGNORE

    USE MODE_DISTRIB_LB
    USE MODE_TOOLS_ll,     ONLY : GET_GLOBALDIMS_ll
    IMPLICIT NONE

    REAL, DIMENSION(:,:,:) , TARGET      :: PLB
    CHARACTER(len=*)                     :: MESSAGE
    REAL                                 :: PRECISION
    CHARACTER(LEN=*),       INTENT(IN)   ::HLBTYPE! 'LBX','LBXU','LBY' or 'LBYV'
    INTEGER,                INTENT(IN) ::KRIM  ! size of the LB area

    !
    ! local var
    !
    REAL,ALLOCATABLE, DIMENSION(:,:,:)       :: TAB_ll,TAB_SON_ll,TAB_SAVE_ll
    REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: Z3D
    REAL,DIMENSION(:,:,:), POINTER           :: TX3DP,TAB_INTERIOR_ll
    INTEGER                              :: IIMAX_ll,IJMAX_ll
    INTEGER                              :: IIU,IJU,IIU_ll,IJU_ll,IKU_ll
    INTEGER                              :: IINFO_ll

    INTEGER,PARAMETER                    :: ITAG = 12345

    INTEGER                              :: I_FIRST_SON
    INTEGER                              :: I_FIRST_FATHER
    REAL                                 :: MAX_DIFF , MAX_VAL
    INTEGER                              :: IIB_ll,IIE_ll,IJB_ll,IJE_ll
    INTEGER                              :: JI
    INTEGER                              :: IIB,IIE,IJB,IJE

    INTEGER                              :: IIU_SON_ll,IJU_SON_ll,IKU_SON_ll
    INTEGER                              :: IIB_SON_ll,IIE_SON_ll,IJB_SON_ll,IJE_SON_ll
    INTEGER                              :: IHEXT_SON_ll , IDIFF_HEXT , IRIM_ll , IRIM_SON_ll

#ifdef MNH_SP4
    !pas de mpi_spawn sur IBM-SP ni MPI_ARGV_NULL etc ...
    RETURN           
#else
    IF ( .NOT. MPPDB_INITIALIZED ) RETURN 
    !
    CALL MPPDB_BARRIER()
    !
    IF(MPPDB_FATHER_WORLD) THEN
       !
       ! Reconstruct the all PLB in TAB_ll
       !
       CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
       IIU_ll = IIMAX_ll+2*JPHEXT
       IJU_ll = IJMAX_ll+2*JPHEXT
       IKU_ll = SIZE(PLB,3)
       IRIM_ll = (KRIM+JPHEXT)*2

       IF (HLBTYPE == 'LBX' .OR. HLBTYPE == 'LBXU') THEN 
          IIU_ll = IRIM_ll
       ELSE
          IJU_ll = IRIM_ll
       END IF
 
       IF (MPPDB_IRANK_WORLD.EQ.0)  THEN
          ! I/O proc case
          ALLOCATE(Z3D(IIU_ll,IJU_ll,SIZE(PLB,3)))
          DO JI = 1,ISNPROC
             CALL GET_DISTRIB_LB(HLBTYPE,JI,'FM','WRITE',KRIM,IIB,IIE,IJB,IJE)
             IF (IIB /= 0) THEN
                TX3DP=>Z3D(IIB:IIE,IJB:IJE,:)
                IF (ISP /= JI) THEN
                   CALL MPI_RECV(TX3DP,SIZE(TX3DP),MPI_PRECISION,JI-1 &
                        ,99,NMNH_COMM_WORLD,MPI_STATUS_IGNORE,IINFO_ll) 
                ELSE
                   CALL GET_DISTRIB_LB(HLBTYPE,JI,'LOC','WRITE',KRIM,IIB,IIE,IJB,IJE)
                   TX3DP = PLB(IIB:IIE,IJB:IJE,:)
                END IF
             END IF
          END DO

          TX3DP=>Z3D

       ELSE
          ! Other processors
          CALL GET_DISTRIB_LB(HLBTYPE,ISP,'LOC','WRITE',KRIM,IIB,IIE,IJB,IJE)
          IF (IIB /= 0) THEN
             TX3DP=>PLB(IIB:IIE,IJB:IJE,:)
             CALL MPI_BSEND(TX3DP,SIZE(TX3DP),MPI_PRECISION,0,99,NMNH_COMM_WORLD,IINFO_ll)
          END IF
       END IF

       IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          !
          ! I'm the first FATHER => recieve the correct globale ARRAY from first son
          !
          !
          ! the first son , is the next processus after this 'world' so
          !
          I_FIRST_SON = MPPDB_NBPROC_WORLD
          !
          ! recieve JPHEXT from son if different
          !
          CALL MPI_RECV(IHEXT_SON_ll,1,MPI_INTEGER,I_FIRST_SON, &
               ITAG, MPPDB_INTRA_COMM,MPI_STATUS_IGNORE, IINFO_ll)

          IIU_SON_ll = IIMAX_ll+2*IHEXT_SON_ll
          IJU_SON_ll = IJMAX_ll+2*IHEXT_SON_ll
          IKU_SON_ll = SIZE(PLB,3)
          IRIM_SON_ll = (KRIM+IHEXT_SON_ll)*2
          !
          IF (HLBTYPE == 'LBX' .OR. HLBTYPE == 'LBXU') THEN 
             IIU_SON_ll = IRIM_SON_ll
          ELSE
             IJU_SON_ll = IRIM_SON_ll
          END IF          
          !
          ALLOCATE(TAB_SON_ll(IIU_SON_ll,IJU_SON_ll,IKU_SON_ll))
          !
          CALL MPI_RECV(TAB_SON_ll,SIZE(TAB_SON_ll),MPI_PRECISION,I_FIRST_SON, &
               ITAG, MPPDB_INTRA_COMM,MPI_STATUS_IGNORE, IINFO_ll)
          !
          IDIFF_HEXT = MIN(JPHEXT,IHEXT_SON_ll)
          !
          ALLOCATE(TAB_SAVE_ll(SIZE(Z3D,1),SIZE(Z3D,2),SIZE(Z3D,3)))
          !
          IF (HLBTYPE == 'LBX' .OR. HLBTYPE == 'LBXU') THEN 

          ELSE
          END IF
          IIB_ll   = 1 + JPHEXT    ; IJB_ll = 1 + JPHEXT
          IIE_ll   = IIU_ll-JPHEXT ; IJE_ll = IJU_ll-JPHEXT
          
          IIB_SON_ll   = 1 + IHEXT_SON_ll    ; IJB_SON_ll = 1 + IHEXT_SON_ll
          IIE_SON_ll   = IIU_SON_ll-IHEXT_SON_ll ; IJE_SON_ll = IJU_SON_ll-IHEXT_SON_ll
          !
          TAB_SAVE_ll = Z3D
          Z3D      = 0.0
          Z3D(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll)   &
            = ABS ( TAB_SAVE_ll(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll) & 
            -       TAB_SON_ll(IIB_SON_ll-IDIFF_HEXT:IIE_SON_ll+IDIFF_HEXT,IJB_SON_ll-IDIFF_HEXT:IJE_SON_ll+IDIFF_HEXT &
                              ,1:IKU_SON_ll) )
          !
          MAX_VAL  = MAXVAL( ABS (TAB_SON_ll(IIB_SON_ll-IDIFF_HEXT:IIE_SON_ll+IDIFF_HEXT,&
                                             IJB_SON_ll-IDIFF_HEXT:IJE_SON_ll+IDIFF_HEXT,1:IKU_SON_ll) ) )
          IF ( MAX_VAL .EQ. 0.0 ) MAX_VAL = 1.0
          !
          MAX_DIFF=MAXVAL(Z3D(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll)/MAX_VAL)
          TAB_INTERIOR_ll=> Z3D(IIB_ll-IDIFF_HEXT:IIE_ll+IDIFF_HEXT,IJB_ll-IDIFF_HEXT:IJE_ll+IDIFF_HEXT,1:IKU_ll)
          !
          IF (MAX_DIFF > PRECISION ) THEN
             print*," MPPDB_CHECKLB :: PB MPPDB_CHECKLB =", MESSAGE ," ERROR=",MAX_DIFF , MAX_VAL
          ELSE
             print*," MPPDB_CHECKLB :: OK MPPDB_CHECKLB =", MESSAGE ," ERROR=",MAX_DIFF , MAX_VAL
          END IF
          flush(unit=OUTPUT_UNIT)
          !
          DEALLOCATE(TAB_SON_ll)
          !
       END IF
    ELSE
       !
       ! SON WORLD 
       !
       IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          !
          ! first son --> send the good array to the first father
          !
          I_FIRST_FATHER = 0
          IHEXT_SON_ll = JPHEXT
          CALL MPI_BSEND(IHEXT_SON_ll,1,MPI_INTEGER,I_FIRST_FATHER, &
               ITAG, MPPDB_INTRA_COMM, IINFO_ll)
          CALL MPI_BSEND(PLB,SIZE(PLB),MPI_PRECISION,I_FIRST_FATHER, &
               ITAG, MPPDB_INTRA_COMM, IINFO_ll)
       END IF
    END IF

    CALL MPPDB_BARRIER()
#endif
  END SUBROUTINE MPPDB_CHECKLB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_CHECK_SURFEX2D(PTAB,MESSAGE,PRECISION,KLUOUT,HTYPE,KIU,KJU)

    USE MODD_PARAMETERS, ONLY : JPHEXT
    USE MODD_VAR_ll    , ONLY : MPI_PRECISION
    USE MODI_GET_1D_MASK
    USE MODI_UNPACK_SAME_RANK
    USE MODI_GET_SURF_MASK_n
    USE MODD_IO_SURF_MNH, ONLY : NHALO
    USE MODD_MNH_SURFEX_n

    IMPLICIT NONE

    REAL, DIMENSION(:), INTENT(IN)         :: PTAB
    CHARACTER(len=*), INTENT(IN)           :: MESSAGE
    REAL, INTENT(IN)                       :: PRECISION
    CHARACTER(LEN=*), INTENT(IN),OPTIONAL  :: HTYPE    ! 'WATER', 'NATURE', 'TOWN', 'SEA', 'FULL' (default is 'FULL')
    INTEGER, INTENT(IN)                    :: KLUOUT   ! output listing logical unit
    INTEGER, INTENT(IN),OPTIONAL           :: KIU    ! size of local subdomain in X direction, useful in case where GET_INDICE_ll does not give the sire of the desired model (e.g. in pgd2)
    INTEGER, INTENT(IN),OPTIONAL           :: KJU    ! size of local subdomain in Y direction
    !
    ! local var
    !
    REAL,ALLOCATABLE, DIMENSION(:)       :: PTAB_UNPACKED
    REAL,ALLOCATABLE, DIMENSION(:,:)     :: ZFIELD2D
    INTEGER                              :: IIU,IJU
    INTEGER                              :: KXOR, KYOR, KXEND, KYEND  ! origin and end of the local physical subdomain
    INTEGER                              :: II,IJ
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: KMASK
    INTEGER                              :: KSIZE
    INTEGER                              :: KSIZE_FULL
    !
    IF ( ( .NOT. MPPDB_INITIALIZED ) ) RETURN
    !
!    IF ( SIZE(PTAB) == 0 ) THEN
!      ALLOCATE(ZFIELD2D(0,0))
!      RETURN
    !
    ! Get the dimensions of the subdomain
    !
    IF ( PRESENT(KIU) .AND. PRESENT(KJU) ) THEN
      IIU = KIU+2*JPHEXT
      IJU = KJU+2*JPHEXT
      KSIZE_FULL = KIU*KJU
    ELSE
      CALL GET_INDICE_ll( KXOR, KYOR, KXEND, KYEND )
      IIU = KXEND-KXOR+1+2*JPHEXT
      IJU = KYEND-KYOR+1+2*JPHEXT
      KSIZE_FULL = (KXEND-KXOR+1)*(KYEND-KYOR+1)
      IF ( PRESENT(HTYPE) .AND. KSIZE_FULL /= SIZE(YSURF_CUR%U%XCOVER,1) .AND. NHALO /= JPHEXT ) THEN
        !IIU = KXEND-KXOR+1+2*JPHEXT+2*NHALO
        !IJU = KYEND-KYOR+1+2*JPHEXT+2*NHALO
        KSIZE_FULL = (KXEND-KXOR+1+2*NHALO) * (KYEND-KYOR+1+2*NHALO)
      ENDIF
    ENDIF
    !
    ! Unpack PTAB
    !
    IF(PRESENT(HTYPE)) THEN
      KSIZE = SIZE( PTAB, 1 )
      ALLOCATE( KMASK(KSIZE) )
      ALLOCATE( PTAB_UNPACKED(KSIZE_FULL) )
      CALL GET_SURF_MASK_n(YSURF_CUR%DTCO,YSURF_CUR%U,HTYPE,KSIZE,KMASK,KSIZE_FULL,KLUOUT)
      CALL UNPACK_SAME_RANK( KMASK, PTAB, PTAB_UNPACKED )
    ELSE
      KSIZE = KSIZE_FULL
      ALLOCATE( PTAB_UNPACKED(KSIZE) )
      PTAB_UNPACKED(:) = PTAB(:)
    ENDIF
    !
    ! Redimension PTAB into a 2D field
    !
    ALLOCATE(ZFIELD2D(IIU,IJU))
    ZFIELD2D = 0.
    DO IJ=1+JPHEXT,IJU-JPHEXT
      DO II=1+JPHEXT,IIU-JPHEXT
         IF(PRESENT(HTYPE)) THEN
            ZFIELD2D(II,IJ) = PTAB_UNPACKED((IJ-JPHEXT-1+NHALO)*(KXEND-KXOR+1+2*NHALO)+II-JPHEXT+NHALO)
         ELSE
            ZFIELD2D(II,IJ) = PTAB_UNPACKED((IJ-JPHEXT-1)*(KXEND-KXOR+1)+II-JPHEXT)
         END IF
      ENDDO
    ENDDO
    !
    ! Call MPPDB_CHECK2D on ZFIELD3D
    !
    IF (MPPDB_IRANK_WORLD.EQ.0) THEN
      write(6,*) ' MPPDB_CHECK_SURFEX2D :'
    ENDIF
    CALL MPPDB_CHECK2D(ZFIELD2D,MESSAGE,PRECISION)

    IF (ALLOCATED(KMASK)) DEALLOCATE( KMASK )
    IF (ALLOCATED(PTAB_UNPACKED)) DEALLOCATE( PTAB_UNPACKED )
    IF (ALLOCATED(ZFIELD2D)) DEALLOCATE( ZFIELD2D )
    !
  END SUBROUTINE MPPDB_CHECK_SURFEX2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE MPPDB_CHECK_SURFEX3D(PTAB,MESSAGE,PRECISION,KLUOUT,HTYPE,KZSIZE)

    USE MODD_PARAMETERS, ONLY : JPHEXT
    USE MODD_VAR_ll    , ONLY : MPI_PRECISION
    USE MODI_GET_1D_MASK
    USE MODI_UNPACK_SAME_RANK
    USE MODI_GET_SURF_MASK_n
    USE MODD_IO_SURF_MNH, ONLY : NHALO
    USE MODD_CONFZ     , ONLY : MPI_BUFFER_SIZE
    USE MODD_MPIF      , ONLY : MPI_INTEGER, MPI_STATUS_IGNORE, MPI_SUM
    USE MODD_MNH_SURFEX_n
!
    IMPLICIT NONE
!
    REAL, DIMENSION(:,:)               :: PTAB
    CHARACTER(len=*)                   :: MESSAGE
    REAL                               :: PRECISION
    CHARACTER(LEN=*), INTENT(IN),OPTIONAL  :: HTYPE    ! 'WATER', 'NATURE', 'TOWN', 'SEA', 'FULL' (default is 'FULL')
    INTEGER, INTENT(IN)                    :: KLUOUT   ! output listing logical unit
    INTEGER, INTENT(IN),OPTIONAL           :: KZSIZE   ! size of Z-dimension. Necessary if PTAB is of size 0 on one process
    !
    ! local var
    !
    REAL,ALLOCATABLE, DIMENSION(:,:)       :: PTAB_UNPACKED
    REAL,ALLOCATABLE, DIMENSION(:,:,:)   :: ZFIELD3D
    INTEGER                              :: IIU,IJU,IKU
    INTEGER                              :: KXOR, KYOR, KXEND, KYEND  ! origin and end of the local physical subdomain
    INTEGER                              :: II,IJ,IK
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: KMASK
    INTEGER                              :: KSIZE
    INTEGER                              :: KSIZEBUF
    INTEGER                              :: KSIZE_FULL
    INTEGER                              :: IGLBSIZEPTAB
    INTEGER                              :: INBSLICES
    INTEGER                              :: IINFO_ll
    !
    IF ( ( .NOT. MPPDB_INITIALIZED ) ) RETURN
    CALL MPI_ALLREDUCE(SIZE(PTAB), IGLBSIZEPTAB, 1,MPI_INTEGER, MPI_SUM, MPPDB_INTRA_COMM, IINFO_ll)
    IF ( IGLBSIZEPTAB == 0 ) RETURN
    !
    IF ( SIZE(PTAB) == 0 ) THEN   !if the local size of the field is 0, we need to define ZFIELD3D filled with default value 1e20
      CALL GET_INDICE_ll( KXOR, KYOR, KXEND, KYEND )
      IIU = KXEND-KXOR+1+2*JPHEXT
      IJU = KYEND-KYOR+1+2*JPHEXT
      IKU = KZSIZE
      ALLOCATE(ZFIELD3D(IIU,IJU,IKU))
      ZFIELD3D = 1.E20
    ELSE
      !
      ! Get the dimensions of the subdomain
      !
      CALL GET_INDICE_ll( KXOR, KYOR, KXEND, KYEND )
      IIU = KXEND-KXOR+1+2*JPHEXT
      IJU = KYEND-KYOR+1+2*JPHEXT
      IKU = SIZE(PTAB,2)
      KSIZE_FULL = (KXEND-KXOR+1)*(KYEND-KYOR+1)
      IF ( PRESENT(HTYPE) .AND. KSIZE_FULL /= SIZE(YSURF_CUR%U%XCOVER,1) .AND. NHALO /= JPHEXT ) THEN
        KSIZE_FULL = (KXEND-KXOR+1+2*NHALO) * (KYEND-KYOR+1+2*NHALO)
      ENDIF
      !
      ! Unpack PTAB
      !
      IF(PRESENT(HTYPE)) THEN
        KSIZE = SIZE( PTAB, 1 )
        ALLOCATE( KMASK(KSIZE) )
        ALLOCATE( PTAB_UNPACKED(KSIZE_FULL,IKU) )
        CALL GET_SURF_MASK_n(YSURF_CUR%DTCO,YSURF_CUR%U,HTYPE,KSIZE,KMASK,KSIZE_FULL,KLUOUT)
        DO II=1,IKU
          CALL UNPACK_SAME_RANK( KMASK, PTAB(:,II), PTAB_UNPACKED(:,II) )
        ENDDO
      ELSE
        KSIZE = KSIZE_FULL
        ALLOCATE( PTAB_UNPACKED(KSIZE,IKU) )
        PTAB_UNPACKED(:,:) = PTAB(:,:)
      ENDIF
      !
      ! Redimension PTAB into a 2D field
      !
      ALLOCATE(ZFIELD3D(IIU,IJU,IKU))
      ZFIELD3D = 0.
      DO IJ=1+JPHEXT,IJU-JPHEXT
        DO II=1+JPHEXT,IIU-JPHEXT
          !ZFIELD3D(II,IJ,:) = PTAB_UNPACKED((IJ-JPHEXT-1)*(KXEND-KXOR+1)+II-JPHEXT,:)
         ZFIELD3D(II,IJ,:) = PTAB_UNPACKED((IJ-JPHEXT-1+NHALO)*(KXEND-KXOR+1+2*NHALO)+II-JPHEXT+NHALO,:)
        ENDDO
      ENDDO
    ENDIF
    !
    ! Call MPPDB_CHECK3D on ZFIELD3D
    !
    ! pour eviter de communiquer des tableaux trop grands qui ne passent pas en memoire,
    ! on "decoupe" le champ en morceaux de taille inferieure a MPI_BUFFER_SIZE*1000000/8
    !ATTENTION : en fait ça ne suffit pas, il faut prendre une limite plus petite
    !je choisi arbitrairement 52*102*102 comme limite a la taille globale du champ
!    IF ( SIZE(ZFIELD3D) > MPI_BUFFER_SIZE*1000000/8 ) THEN
!      KSIZEBUF = SIZE(ZFIELD3D,3)*8/MPI_BUFFER_SIZE*1000000
!    IF ( SIZE(ZFIELD3D) > 52*102*102 ) THEN
    IF ( IGLBSIZEPTAB > MPI_BUFFER_SIZE*1000000/16 ) THEN
      INBSLICES = IGLBSIZEPTAB/(MPI_BUFFER_SIZE*1000000/16)
      IF (SIZE(ZFIELD3D,3) >= INBSLICES ) THEN
        KSIZEBUF = SIZE(ZFIELD3D,3)/INBSLICES
      ELSE
        write(6,*) ' MPPDB_CHECK_SURFEX3D : field \"',MESSAGE,'\" is too large to be checked with MPPDB. No checking was done...'
      ENDIF
!    IF ( IGLBSIZEPTAB > 52*102*102 ) THEN
!      INBSLICES = 52
!      IF (SIZE(ZFIELD3D,3) >=52 ) THEN
!        KSIZEBUF = SIZE(ZFIELD3D,3)/52
!      ELSE
!        KSIZEBUF=1
!      ENDIF
      DO IK=1,INBSLICES
        IF (MPPDB_IRANK_WORLD.EQ.0) THEN
          IF ( KSIZEBUF*INBSLICES==SIZE(ZFIELD3D,3) ) THEN
            write(6,*) ' MPPDB_CHECK_SURFEX3D part ',IK,'/',INBSLICES,' :'
          ELSE
            write(6,*) ' MPPDB_CHECK_SURFEX3D part ',IK,'/',INBSLICES+1,' :'
          ENDIF
        ENDIF
        CALL MPPDB_CHECK3D(ZFIELD3D(:,:,(IK-1)*KSIZEBUF+1:IK*KSIZEBUF),MESSAGE,PRECISION)
      ENDDO
        IF ( KSIZEBUF*INBSLICES==SIZE(ZFIELD3D,3) ) THEN
        ELSE
          IF (MPPDB_IRANK_WORLD.EQ.0) THEN
            write(6,*) ' MPPDB_CHECK_SURFEX3D part ',IK,'/',INBSLICES+1,' :'
          ENDIF
        CALL MPPDB_CHECK3D(ZFIELD3D(:,:,KSIZEBUF*INBSLICES+1:),MESSAGE,PRECISION)
        ENDIF
    ELSE
      IF (MPPDB_IRANK_WORLD.EQ.0) THEN
        write(6,*) ' MPPDB_CHECK_SURFEX3D :'
      ENDIF
      CALL MPPDB_CHECK3D(ZFIELD3D,MESSAGE,PRECISION)
    ENDIF
    IF (ALLOCATED(KMASK)) DEALLOCATE( KMASK )
    IF (ALLOCATED(PTAB_UNPACKED)) DEALLOCATE( PTAB_UNPACKED )
    IF (ALLOCATED(ZFIELD3D)) DEALLOCATE( ZFIELD3D )
    !
  END SUBROUTINE MPPDB_CHECK_SURFEX3D


END MODULE MODE_MPPDB


