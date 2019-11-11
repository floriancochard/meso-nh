!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    MODIFICATIONS
!!    -------------
!!
!!  J.Escobar 3/12/2014 : typo form -> from
!!  Philippe 03/10/2017: set IP and NPROC in INIT_NMNH_COMM_WORLD
!!  Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!!
MODULE MODE_MNH_WORLD
  IMPLICIT NONE
  CHARACTER(len=*), parameter   :: conf_mnh_world="conf_mnh_world.nam"
  CHARACTER(len=10)             :: cmapping = "DEFAULT"

  NAMELIST /NAM_CONF_MNH_WORLD/ cmapping

CONTAINS

  SUBROUTINE  INIT_NMNH_COMM_WORLD(KINFO_ll)
    !JUANZ
    USE MODD_MPIF  , ONLY : MPI_COMM_WORLD, MPI_CHARACTER
#ifdef MNH_GA
    USE MODD_MPIF  , ONLY :  MPI_THREAD_MULTIPLE
#endif
    USE MODD_VAR_ll, ONLY : IP, NPROC, NMNH_COMM_WORLD , LMNH_ISINIT 
    !JUANZ
    IMPLICIT NONE

    INTEGER :: KINFO_ll

    INTEGER :: IRANK,IPROC
    LOGICAL :: GISINIT

    !JUANZ
    INTEGER                         :: myrank_key,new_rank,new_size
    INTEGER                         :: COLOR 
    INTEGER                         :: IX_DIM,I
    INTEGER, DIMENSION(:) , POINTER :: IX_COORD,IY_COORD
    INTEGER                         :: IERR,IROOT
#ifdef MNH_GA
    INTEGER                         :: REQUIRED=MPI_THREAD_MULTIPLE,PROVIDED 
#endif
    !JUANZ
    INTEGER :: ILU

    !
    KINFO_ll = 0
    CALL MPI_INITIALIZED(GISINIT, KINFO_ll)
    IF (.NOT. GISINIT) THEN
#ifdef MNH_GA
       CALL MPI_INIT_thread(REQUIRED,PROVIDED,KINFO_ll)
#else
       CALL MPI_INIT(KINFO_ll)
#endif
       !
       CALL MPI_COMM_RANK(MPI_COMM_WORLD, IRANK, KINFO_ll)
       !
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD, IPROC, KINFO_ll)
       IERR = 1
       !
       ! Read namelist config file
       !
       IF ( irank .EQ. 0 ) THEN
          OPEN(newunit=ILU,form="formatted",file=conf_mnh_world,STATUS='OLD',iostat=IERR)
          ! Read IO parameter
          IF (IERR.EQ.0) THEN
             READ(unit=ilu,NML=NAM_CONF_MNH_WORLD)
             WRITE(*,NAM_CONF_MNH_WORLD)
          ENDIF
          CLOSE(unit=ILU)
       ENDIF
       iroot = 0
       ! Brodcast mapping
       CALL MPI_BCAST(cmapping,10,MPI_CHARACTER,iroot,MPI_COMM_WORLD,KINFO_ll)

       SELECT CASE (cmapping)
          !
!       CASE ("HILBERT2D")
!          !
!          ! 2D hilbert
!          !
!          IX_DIM = 2** INT (0.5 + log(1.0*IPROC)/log(2.0)/2.0)
!          CALL HILBERT_REMAP(IX_DIM,IX_COORD,IY_COORD)
!          myrank_key =  IY_COORD(IRANK+1) + IX_COORD(IRANK+1) * IX_DIM
!          IF (IRANK .EQ. 0 ) THEN
!             print *,"IX_DIM=",IX_DIM
!             OPEN(newunit=ILU,form="formatted",file="hilbert_2D")
!             DO  I=1,IPROC
!                write(unit=ILU,fmt=*) I, IX_COORD(I),IY_COORD(I) , &
!                     1+ IY_COORD(I) + IX_COORD(I) * IX_DIM 
!             END DO
!             CLOSE(unit=ILU)
!          END IF
!          color = 1
!          call MPI_Comm_split( MPI_COMM_WORLD, color, myrank_key,NMNH_COMM_WORLD , KINFO_ll )
!          call MPI_Comm_rank ( NMNH_COMM_WORLD, new_rank, KINFO_ll)
!          call MPI_Comm_size ( NMNH_COMM_WORLD, new_size, KINFO_ll )
          !
       CASE default 
          !
          ! mpi_comm_world mapping
          !
          !myrank_key = irank
          NMNH_COMM_WORLD = MPI_COMM_WORLD
          !        
       END SELECT
       !
       !JUANZ create new/remapped communicator 
       !
    END IF
    IF (.NOT. LMNH_ISINIT ) THEN
       LMNH_ISINIT = .TRUE.
       !
       CALL MPI_COMM_RANK(NMNH_COMM_WORLD, IP, KINFO_ll)
       IP = IP + 1
       !
       CALL MPI_COMM_SIZE(NMNH_COMM_WORLD, NPROC, KINFO_ll)
       !
       IF ( IP .EQ. 1 ) THEN
          PRINT*,"hello MNH world from IP=",IP," nproc=",NPROC
       END IF
    END IF

  END SUBROUTINE INIT_NMNH_COMM_WORLD

!  RECURSIVE SUBROUTINE HILBERT_REMAP(KN,KX_COORD,KY_COORD)
!
!    IMPLICIT NONE
!    INTEGER, INTENT(IN)                            :: KN   ! >= sqrt(Nproc number)
!    INTEGER, DIMENSION(:) , POINTER, INTENT(INOUT)   :: KX_COORD,KY_COORD
!
!    !local variable
!    INTEGER, DIMENSION(:) , POINTER                :: IX_COORD_N2,IY_COORD_N2
!    INTEGER                                        :: IX_SIZE_N2 ,IY_SIZE_N2
!    INTEGER                                        :: IB, IE
!
!    IF ( KN .EQ. 2 ) THEN
!       ALLOCATE(KX_COORD(4),KY_COORD(4))
!       KX_COORD(:) = (/ 0 ,0 ,1 ,1 /)
!       KY_COORD(:) = (/ 0 ,1 ,1 ,0 /)
!    ELSE
!       CALL  HILBERT_REMAP(KN/2,IX_COORD_N2,IY_COORD_N2)
!       IX_SIZE_N2 = SIZE(IX_COORD_N2)
!       IY_SIZE_N2 = SIZE(IY_COORD_N2)
!       ALLOCATE(KX_COORD(2*IX_SIZE_N2+2*IY_SIZE_N2))
!       ALLOCATE(KY_COORD(2*IX_SIZE_N2+2*IY_SIZE_N2))
!       !
!       ! xl =  [ y // x // N/2 + x // N-1-y ] 
!       !
!       IB = 1      ; IE =      IY_SIZE_N2 ; KX_COORD(IB:IE) =   IY_COORD_N2
!       IB = IE + 1 ; IE = IE + IX_SIZE_N2 ; KX_COORD(IB:IE) =   IX_COORD_N2
!       IB = IE + 1 ; IE = IE + IX_SIZE_N2 ; KX_COORD(IB:IE) =   IX_COORD_N2 + KN/2
!       IB = IE + 1 ; IE = IE + IY_SIZE_N2 ; KX_COORD(IB:IE) = - IY_COORD_N2 + KN-1
!       !
!       ! yl =  [ x // N/2+y // N/2 + y // N/2-1-X ] 
!       !
!       IB = 1      ; IE =      IX_SIZE_N2 ; KY_COORD(IB:IE) =   IX_COORD_N2
!       IB = IE + 1 ; IE = IE + IY_SIZE_N2 ; KY_COORD(IB:IE) =   IY_COORD_N2 + KN/2
!       IB = IE + 1 ; IE = IE + IY_SIZE_N2 ; KY_COORD(IB:IE) =   IY_COORD_N2 + KN/2
!       IB = IE + 1 ; IE = IE + IX_SIZE_N2 ; KY_COORD(IB:IE) = - IX_COORD_N2 + KN/2 -1 
!       
!       DEALLOCATE(IX_COORD_N2,IY_COORD_N2)
!
!    END IF
!
!
!  END SUBROUTINE HILBERT_REMAP


END MODULE MODE_MNH_WORLD
