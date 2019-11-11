!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

#ifdef MNH_MPI_DOUBLE_PRECISION
#define MPI_FLOAT MPI_DOUBLE_PRECISION
#else
#define MPI_FLOAT MPI_REAL
#endif

MODULE MODE_FMREAD
!
!Correction :
!  J.Escobar : 22/08/2005 : BUG : manque un "GOTO 1000" si champs
!              lue non trouvé !!!
!
USE MODD_MPIF
IMPLICIT NONE 

PRIVATE

INTERFACE FMREAD
  MODULE PROCEDURE FMREADX0_ll,FMREADX1_ll,FMREADX2_ll,FMREADX3_ll,&
       & FMREADX4_ll,FMREADX5_ll,FMREADX6_ll,&
       & FMREADN0_ll,FMREADN1_ll,FMREADN2_ll,&
       & FMREADL0_ll,FMREADL1_ll,FMREADC0_ll,FMREADT0_ll
END INTERFACE
!

PUBLIC FMREAD_LB,FMREAD,FMREADX0_ll,FMREADX1_ll,FMREADX2_ll,FMREADX3_ll,&
       & FMREADX4_ll,FMREADX5_ll,FMREADX6_ll,&
       & FMREADN0_ll,FMREADN1_ll,FMREADN2_ll,&
       & FMREADL0_ll,FMREADL1_ll,FMREADC0_ll,FMREADT0_ll

!INCLUDE 'mpif.h'

CONTAINS 
SUBROUTINE FM_READ_ERR(HFUNC,HFILEM,HFIPRI,HRECFM,HDIR,KRESP)
USE MODE_FM, ONLY : FMLOOK_ll

CHARACTER(LEN=*) :: HFUNC 
CHARACTER(LEN=*) :: HFILEM
CHARACTER(LEN=*) :: HFIPRI
CHARACTER(LEN=*) :: HRECFM
CHARACTER(LEN=*) :: HDIR
INTEGER          :: KRESP

INTEGER          :: ILUPRI
INTEGER          :: IRESP

CALL FMLOOK_ll(HFIPRI,HFIPRI,ILUPRI,IRESP)
WRITE (ILUPRI,*) ' exit from ',HFUNC, ' with RESP:',KRESP
WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
WRITE (ILUPRI,*) '   | HRECFM = ',HRECFM
WRITE (ILUPRI,*) '   | HDIR  = ',HDIR

END SUBROUTINE FM_READ_ERR


SUBROUTINE BCAST_HEADER(TPFD,TPFMH)
USE MODE_FD_ll, ONLY : FD_ll
USE MODD_FM
TYPE(FD_ll),     POINTER    :: TPFD
TYPE(FMHEADER), INTENT(IN) :: TPFMH

INTEGER :: ierr 

CALL MPI_BCAST(TPFMH%GRID,1,MPI_INTEGER,TPFD%OWNER-1,TPFD%COMM,IERR)
CALL MPI_BCAST(TPFMH%COMLEN,1,MPI_INTEGER,TPFD%OWNER-1,TPFD%COMM,IERR)
CALL MPI_BCAST(TPFMH%COMMENT,TPFMH%COMLEN,MPI_CHARACTER,TPFD%OWNER-1,TPFD%COMM,IERR)

END SUBROUTINE BCAST_HEADER

SUBROUTINE FMREADX0_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
!
!*      0.    DECLARATIONS
!             ------------
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),             INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),             INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),             INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),             INTENT(IN) ::HDIR   ! field form
REAL,                         INTENT(INOUT)::PFIELD ! array containing the data field 
INTEGER,                      INTENT(INOUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                      INTENT(INOUT)::KLENCH ! length of comment string
CHARACTER(LEN=*),             INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                      INTENT(INOUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
!----------------------------------------------------------------
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
TYPE(FMHEADER)               :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,1,PFIELD,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,1,PFIELD,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    CALL MPI_BCAST(PFIELD,1,MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADX0_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
KRESP = IRESP
RETURN
    
END SUBROUTINE FMREADX0_ll

SUBROUTINE FMREADX1_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll
!
!*      0.    DECLARATIONS
!             ------------
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),        INTENT(IN) ::HRECFM   ! name of the article to read
CHARACTER(LEN=*),        INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),        INTENT(IN) ::HDIR     ! Field form
REAL,DIMENSION(:),TARGET,INTENT(INOUT)::PFIELD   ! array containing the data field 
INTEGER,                 INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                 INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),        INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                 INTENT(INOUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
!----------------------------------------------------------------
CHARACTER(LEN=JPFINL)     :: YFNLFI
INTEGER                   :: IERR
TYPE(FD_ll), POINTER      :: TZFD
INTEGER                   :: IRESP
REAL,DIMENSION(:),POINTER :: ZFIELDP
LOGICAL                   :: GALLOC
TYPE(FMHEADER)            :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    ELSE
      ALLOCATE(ZFIELDP(0))
      GALLOC = .TRUE.
    END IF
      
    !  
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    IF (HDIR /= 'XX' .AND. HDIR /='YY') THEN
      ! Broadcast Field
      CALL MPI_BCAST(PFIELD,SIZE(PFIELD),MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
    ELSE 
      !Scatter Field
      CALL SCATTER_XXFIELD(HDIR,ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM) 
    END IF
  END IF !(GSMONOPROC)
  
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADX1_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF

IF (GALLOC) DEALLOCATE (ZFIELDP)
KRESP = IRESP
RETURN
!------------------------------------------------------------------
END SUBROUTINE FMREADX1_ll

SUBROUTINE FMREADX2_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll

CHARACTER(LEN=*),           INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),           INTENT(IN) ::HRECFM   ! name of the article to read
CHARACTER(LEN=*),           INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),           INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:,:),TARGET, INTENT(INOUT)::PFIELD   ! array containing the data field
INTEGER,                    INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                    INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),           INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                   INTENT(INOUT)::KRESP     ! return-code
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
REAL,DIMENSION(:,:), POINTER :: ZFIELDP
LOGICAL                      :: GALLOC
TYPE(FMHEADER)               :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'

!------------------------------------------------------------------

TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:)=SPREAD(SPREAD(PFIELD(2,2),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
      ZFIELDP=>PFIELD(:,2:2)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:)=SPREAD(PFIELD(:,2),DIM=2,NCOPIES=3)
    ELSE
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF
    IF (IRESP /= 0) GOTO 1000
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      ! I/O processor case
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC) 
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    ELSE
      ALLOCATE(ZFIELDP(0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      ! XX or YY Scatter Field
      CALL SCATTER_XXFIELD(HDIR,ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM) 
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        ! 2D compact case
        CALL SCATTER_XXFIELD('XX',ZFIELDP(:,1),PFIELD(:,2),TZFD%OWNER,TZFD%COMM)
        PFIELD(:,:) = SPREAD(PFIELD(:,2),DIM=2,NCOPIES=3)
      ELSE
        ! XY Scatter Field
        CALL SCATTER_XYFIELD(ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM)
      END IF
    ELSE
      CALL MPI_BCAST(PFIELD,SIZE(PFIELD),MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
    END IF
  END IF !(GSMONOPROC)
  
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADX2_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
IF (GALLOC) DEALLOCATE (ZFIELDP)
KRESP = IRESP
RETURN
!------------------------------------------------------------------


END SUBROUTINE FMREADX2_ll

SUBROUTINE FMREADX3_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll

CHARACTER(LEN=*),             INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),             INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),             INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),             INTENT(IN) ::HDIR   ! field form
REAL, DIMENSION(:,:,:),TARGET,INTENT(INOUT)::PFIELD ! array containing the data field
INTEGER,                      INTENT(INOUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                      INTENT(INOUT)::KLENCH ! length of comment string
CHARACTER(LEN=*),             INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                      INTENT(INOUT)::KRESP    ! return-code
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                    :: YFNLFI
INTEGER                                  :: IERR
TYPE(FD_ll), POINTER                     :: TZFD
INTEGER                                  :: IRESP
REAL,DIMENSION(:,:,:),POINTER            :: ZFIELDP
LOGICAL                                  :: GALLOC
TYPE(FMHEADER)                           :: TZFMH
!JUAN
INTEGER                                  :: JK
CHARACTER(LEN=LEN(HRECFM))               :: YK,YRECZSLIDE
REAL,DIMENSION(:,:),POINTER              :: ZSLIDE_ll,ZSLIDE
!JUAN
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP  = 0
YFNLFI = TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D  .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2,:)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:,:)=SPREAD(SPREAD(PFIELD(2,2,:),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
      ALLOCATE (ZFIELDP(SIZE(PFIELD,1),1,SIZE(PFIELD,3)))
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:,:)=SPREAD(ZFIELDP(:,1,:),DIM=2,NCOPIES=3)
    ELSE
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF
    IF (IRESP /= 0) GOTO 1000
  ELSE ! multiprocessor execution
!
!JUAN BG Z SLIDE 
!
    IF (ISP == TZFD%OWNER)  THEN
!JUAN      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
      ZSLIDE => PFIELD(:,:,1)
      CALL ALLOCBUFFER_ll(ZSLIDE_ll,ZSLIDE,HDIR,GALLOC)
!JUAN      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
!JUAN           & ,IRESP)
    ELSE 
!JUAN      ALLOCATE(ZFIELDP(0,0,0))
      ALLOCATE(ZSLIDE_ll(0,0))
      GALLOC = .TRUE. 
    END IF

Z_SLIDE: DO JK=1,SIZE(PFIELD,3)
!JUAN
    WRITE(YK,'(I4.4)')  JK
    YRECZSLIDE = TRIM(HRECFM)//YK

    IF (ISP == TZFD%OWNER)  THEN
!JUAN      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
!JUAN           & ,IRESP)
      CALL FM_READ_ll(TZFD%FLU,YRECZSLIDE,.TRUE.,SIZE(ZSLIDE_ll),ZSLIDE_ll,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
!JUAN    !
!JUAN    CALL BCAST_HEADER(TZFD,TZFMH)
!JUAN    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      ! XX or YY Scatter Field
       STOP " XX NON PREVU SUR BG POUR LE MOMENT "
      CALL SCATTER_XXFIELD(HDIR,ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM) 
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        ! 2D compact case
       STOP " L2D NON PREVU SUR BG POUR LE MOMENT "
        CALL SCATTER_XXFIELD('XX',ZFIELDP(:,1,:),PFIELD(:,2,:),TZFD%OWNER,TZFD%COMM)
        PFIELD(:,:,:) = SPREAD(PFIELD(:,2,:),DIM=2,NCOPIES=3)
      ELSE
        ! XY Scatter Field
        ZSLIDE => PFIELD(:,:,JK)
!JUAN   CALL SCATTER_XYFIELD(ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM)
        CALL SCATTER_XYFIELD(ZSLIDE_ll,ZSLIDE,TZFD%OWNER,TZFD%COMM)
      END IF
    ELSE
      ! Broadcast Field
       STOP "  Broadcast Field NON PREVU SUR BG POUR LE MOMENT "
      CALL MPI_BCAST(PFIELD,SIZE(PFIELD),MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
    END IF
 END DO Z_SLIDE
 !
 CALL BCAST_HEADER(TZFD,TZFMH)
 !
!JUAN BG Z SLIDE  
  END IF !(GSMONOPROC) 
  
  KGRID    = TZFMH%GRID
  KLENCH   = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61          
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADX3_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
!JUAN IF (GALLOC) DEALLOCATE (ZFIELDP)
IF (GALLOC) DEALLOCATE (ZSLIDE_ll)
KRESP = IRESP
RETURN
!------------------------------------------------------------------
END SUBROUTINE FMREADX3_ll

SUBROUTINE FMREADX4_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll

CHARACTER(LEN=*),              INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),              INTENT(IN) ::HRECFM   ! name of the article to read
CHARACTER(LEN=*),              INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),              INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:,:,:,:),TARGET,INTENT(INOUT)::PFIELD   ! array containing the data field
INTEGER,                       INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                       INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),              INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                       INTENT(INOUT)::KRESP  ! return-code if
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)           :: YFNLFI
INTEGER                         :: IERR
TYPE(FD_ll), POINTER            :: TZFD
INTEGER                         :: IRESP
REAL,DIMENSION(:,:,:,:),POINTER :: ZFIELDP
LOGICAL                         :: GALLOC
TYPE(FMHEADER)                  :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2,:,:)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:,:,:)=SPREAD(SPREAD(PFIELD(2,2,:,:),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
      ZFIELDP=>PFIELD(:,2:2,:,:)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:,:,:)=SPREAD(PFIELD(:,2,:,:),DIM=2,NCOPIES=3)
    ELSE
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF
    IF (IRESP /= 0) GOTO 1000
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      ! XX or YY Scatter Field
      CALL SCATTER_XXFIELD(HDIR,ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM) 
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        ! 2D compact case
        CALL SCATTER_XXFIELD('XX',ZFIELDP(:,1,:,:),PFIELD(:,2,:,:),TZFD%OWNER,TZFD%COMM)
        PFIELD(:,:,:,:) = SPREAD(PFIELD(:,2,:,:),DIM=2,NCOPIES=3)
      ELSE
        ! XY Scatter Field
        CALL SCATTER_XYFIELD(ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM)
      END IF
    ELSE
      ! Broadcast Field
      CALL MPI_BCAST(PFIELD,SIZE(PFIELD),MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
    END IF
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADX4_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF

IF (GALLOC) DEALLOCATE (ZFIELDP)
KRESP = IRESP
RETURN
!------------------------------------------------------------------
END SUBROUTINE FMREADX4_ll

SUBROUTINE FMREADX5_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll

CHARACTER(LEN=*),                INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),                INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),                INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),                INTENT(IN) ::HDIR   ! field form
REAL,DIMENSION(:,:,:,:,:),TARGET,INTENT(INOUT)::PFIELD ! array containing the data field
INTEGER,                         INTENT(INOUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                         INTENT(INOUT)::KLENCH ! length of comment string
CHARACTER(LEN=*),                INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                         INTENT(INOUT)::KRESP  ! return-code
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)             :: YFNLFI
INTEGER                           :: IERR
TYPE(FD_ll), POINTER              :: TZFD
INTEGER                           :: IRESP
REAL,DIMENSION(:,:,:,:,:),POINTER :: ZFIELDP
LOGICAL                           :: GALLOC
TYPE(FMHEADER)                    :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2,:,:,:)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:,:,:,:)=SPREAD(SPREAD(PFIELD(2,2,:,:,:),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
      ZFIELDP=>PFIELD(:,2:2,:,:,:)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
      PFIELD(:,:,:,:,:)=SPREAD(PFIELD(:,2,:,:,:),DIM=2,NCOPIES=3)
    ELSE
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF  
    IF (IRESP /= 0) GOTO 1000
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      ! XX or YY Scatter Field
      CALL SCATTER_XXFIELD(HDIR,ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM) 
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        ! 2D compact case
        CALL SCATTER_XXFIELD('XX',ZFIELDP(:,1,:,:,:),PFIELD(:,2,:,:,:),&
             & TZFD%OWNER,TZFD%COMM)
        PFIELD(:,:,:,:,:) = SPREAD(PFIELD(:,2,:,:,:),DIM=2,NCOPIES=3)
      ELSE
        ! XY Scatter Field
        CALL SCATTER_XYFIELD(ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM)
      END IF
    ELSE
      ! Broadcast Field
      CALL MPI_BCAST(PFIELD,SIZE(PFIELD),MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
    END IF
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADX5_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
IF (GALLOC) DEALLOCATE (ZFIELDP)
KRESP = IRESP
RETURN
!------------------------------------------------------------------
END SUBROUTINE FMREADX5_ll

SUBROUTINE FMREADX6_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll

CHARACTER(LEN=*),                  INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),                  INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),                  INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),                  INTENT(IN) ::HDIR   ! field form
REAL,DIMENSION(:,:,:,:,:,:),TARGET,INTENT(INOUT)::PFIELD ! array containing the data field
INTEGER,                           INTENT(INOUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                           INTENT(INOUT)::KLENCH ! length of comment string
CHARACTER(LEN=*),                  INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                           INTENT(INOUT)::KRESP  ! return-code
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)               :: YFNLFI
INTEGER                             :: IERR
TYPE(FD_ll), POINTER                :: TZFD
INTEGER                             :: IRESP
REAL,DIMENSION(:,:,:,:,:,:),POINTER :: ZFIELDP
LOGICAL                             :: GALLOC
TYPE(FMHEADER)                      :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      ! XX or YY Scatter Field
      CALL SCATTER_XXFIELD(HDIR,ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM) 
    ELSE IF (HDIR == 'XY') THEN
      ! XY Scatter Field
      CALL SCATTER_XYFIELD(ZFIELDP,PFIELD,TZFD%OWNER,TZFD%COMM)
    ELSE
      ! Broadcast Field
      CALL MPI_BCAST(PFIELD,SIZE(PFIELD),MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
    END IF
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADX6_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
IF (GALLOC) DEALLOCATE (ZFIELDP)
KRESP = IRESP
RETURN
!------------------------------------------------------------------
END SUBROUTINE FMREADX6_ll

SUBROUTINE FMREADN0_ll(HFILEM,HRECFM,HFIPRI,HDIR,KFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL

!
!*      0.    DECLARATIONS
!             ------------
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),          INTENT(IN) ::HRECFM   ! name of the article to read
CHARACTER(LEN=*),          INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),          INTENT(IN) ::HDIR     ! Field form
INTEGER,                   INTENT(INOUT)::KFIELD ! array containing the data field     
INTEGER,                   INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),          INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                   INTENT(INOUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
TYPE(FMHEADER)               :: TZFMH

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!  
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,1,KFIELD,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,1,KFIELD,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)       
    !
    CALL MPI_BCAST(KFIELD,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADN0_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
KRESP = IRESP
RETURN

END SUBROUTINE FMREADN0_ll

SUBROUTINE FMREADN1_ll(HFILEM,HRECFM,HFIPRI,HDIR,KFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll

!*      0.    DECLARATIONS
!             ------------
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),           INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),           INTENT(IN) ::HRECFM   ! name of the article to read
CHARACTER(LEN=*),           INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),           INTENT(IN) ::HDIR     ! Field form
INTEGER,DIMENSION(:),TARGET,INTENT(INOUT)::KFIELD ! array containing the data field     
INTEGER,                    INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                    INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),           INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                    INTENT(INOUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)            :: YFNLFI
INTEGER                          :: IERR
TYPE(FD_ll), POINTER             :: TZFD
INTEGER                          :: IRESP
INTEGER,DIMENSION(:),POINTER     :: IFIELDP
LOGICAL                          :: GALLOC
TYPE(FMHEADER)                   :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(KFIELD),KFIELD,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(IFIELDP,KFIELD,HDIR,GALLOC)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP,TZFMH&
           & ,IRESP)
    ELSE
      ALLOCATE(IFIELDP(0))
      GALLOC = .TRUE.
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    IF (HDIR /= 'XX' .AND. HDIR /='YY') THEN
      ! Broadcast Field
      CALL MPI_BCAST(KFIELD,SIZE(KFIELD),MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    ELSE 
      !Scatter Field
      CALL SCATTER_XXFIELD(HDIR,IFIELDP,KFIELD,TZFD%OWNER,TZFD%COMM) 
    END IF
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADN1_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
IF (GALLOC) DEALLOCATE (IFIELDP)
KRESP = IRESP
RETURN
  
END SUBROUTINE FMREADN1_ll

SUBROUTINE FMREADN2_ll(HFILEM,HRECFM,HFIPRI,HDIR,KFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_SCATTER_ll
USE MODE_ALLOCBUFFER_ll

CHARACTER(LEN=*),              INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),              INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),              INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),              INTENT(IN) ::HDIR   ! field form
INTEGER, DIMENSION(:,:),TARGET,INTENT(INOUT)::KFIELD ! array containing the data field
INTEGER,                       INTENT(INOUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                       INTENT(INOUT)::KLENCH ! length of comment string
CHARACTER(LEN=*),              INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                       INTENT(INOUT)::KRESP  ! return-code
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)          :: YFNLFI
INTEGER                        :: IERR
TYPE(FD_ll), POINTER           :: TZFD
INTEGER                        :: IRESP
INTEGER,DIMENSION(:,:),POINTER :: IFIELDP
LOGICAL                        :: GALLOC
TYPE(FMHEADER)                 :: TZFMH
!
!*      1.1   THE NAME OF LFIFM
!
GALLOC = .FALSE.
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(KFIELD,1)==3 .AND. SIZE(KFIELD,2)==3) THEN 
      IFIELDP=>KFIELD(2:2,2:2)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP,TZFMH,IRESP)
      KFIELD(:,:)=SPREAD(SPREAD(KFIELD(2,2),DIM=1,NCOPIES=3),DIM=2,NCOPIES=3)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(KFIELD,2)==3) THEN
      IFIELDP=>KFIELD(:,2:2)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP,TZFMH,IRESP)
      KFIELD(:,:)=SPREAD(KFIELD(:,2),DIM=2,NCOPIES=3)
    ELSE
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(KFIELD),KFIELD,TZFMH,IRESP)
    END IF
    IF (IRESP /= 0) GOTO 1000
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(IFIELDP,KFIELD,HDIR,GALLOC)
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP&
           & ,TZFMH,IRESP)
    ELSE
      ALLOCATE(IFIELDP(0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      ! XX or YY Scatter Field
      CALL SCATTER_XXFIELD(HDIR,IFIELDP,KFIELD,TZFD%OWNER,TZFD&
           & %COMM) 
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        ! 2D compact case
        CALL SCATTER_XXFIELD('XX',IFIELDP(:,1),KFIELD(:,2),TZFD%OWNER,TZFD%COMM)
        KFIELD(:,:) = SPREAD(KFIELD(:,2),DIM=2,NCOPIES=3)
      ELSE
        ! XY Scatter Field
        CALL SCATTER_XYFIELD(IFIELDP,KFIELD,TZFD%OWNER,TZFD%COMM)
      END IF
    ELSE
      ! Broadcast Field
      IF (ISP == TZFD%OWNER) KFIELD = IFIELDP
      CALL MPI_BCAST(KFIELD,SIZE(KFIELD),MPI_INTEGER,TZFD%OWNER-1&
           & ,TZFD%COMM,IERR)
    END IF
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADN2_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
!
IF (GALLOC) DEALLOCATE (IFIELDP)
KRESP = IRESP
RETURN
!------------------------------------------------------------------
END SUBROUTINE FMREADN2_ll


SUBROUTINE FMREADL0_ll(HFILEM,HRECFM,HFIPRI,HDIR,OFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL

!*      0.    DECLARATIONS
!             ------------
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),          INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),          INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),          INTENT(IN) ::HDIR   ! field form
LOGICAL,                   INTENT(INOUT)::OFIELD ! array containing the data field
INTEGER,                   INTENT(INOUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(INOUT)::KLENCH ! length of comment string
CHARACTER(LEN=*),          INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                   INTENT(INOUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
INTEGER                      :: IFIELD
TYPE(FMHEADER)               :: TZFMH

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,1,IFIELD,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,1,IFIELD,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    CALL MPI_BCAST(IFIELD,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,&
         & IERR)
  END IF
  IF (IFIELD==1) THEN
    OFIELD=.TRUE.
  ELSE
    OFIELD=.FALSE.
  END IF
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADL0_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
KRESP = IRESP
RETURN

END SUBROUTINE FMREADL0_ll

SUBROUTINE FMREADL1_ll(HFILEM,HRECFM,HFIPRI,HDIR,OFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
!
!*      0.    DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN) ::HFILEM  ! FM-file name
CHARACTER(LEN=*),          INTENT(IN) ::HRECFM  ! name of the article to read
CHARACTER(LEN=*),          INTENT(IN) ::HFIPRI  ! output file for error messages
CHARACTER(LEN=*),          INTENT(IN) ::HDIR    ! Field form
LOGICAL, DIMENSION(:),     INTENT(INOUT)::OFIELD  ! array containing the data field
INTEGER,                   INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),          INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                   INTENT(INOUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!

CHARACTER(LEN=JPFINL)            :: YFNLFI
INTEGER                          :: IERR
TYPE(FD_ll), POINTER             :: TZFD
INTEGER                          :: IRESP
INTEGER, DIMENSION(SIZE(OFIELD)) :: IFIELD
TYPE(FMHEADER)                   :: TZFMH

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELD),IFIELD,TZFMH&
         & ,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELD),IFIELD,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    CALL MPI_BCAST(IFIELD,SIZE(IFIELD),MPI_INTEGER,TZFD%OWNER-1,TZFD&
       & %COMM,IERR)
  END IF
  WHERE (IFIELD==1)
    OFIELD=.TRUE.
  ELSEWHERE
    OFIELD=.FALSE.
  END WHERE
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADL1_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
KRESP = IRESP
RETURN

END SUBROUTINE FMREADL1_ll

SUBROUTINE FMREADC0_ll(HFILEM,HRECFM,HFIPRI,HDIR,HFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
!
!*      0.    DECLARATIONS
!             ------------
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),          INTENT(IN) ::HRECFM   ! name of the article to read
CHARACTER(LEN=*),          INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),          INTENT(IN) ::HDIR     ! Field form
CHARACTER(LEN=*),          INTENT(INOUT)::HFIELD   ! array containing the data field    
INTEGER,                   INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),          INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                   INTENT(INOUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)             :: YFNLFI
INTEGER                           :: IERR
TYPE(FD_ll), POINTER              :: TZFD
INTEGER                           :: IRESP
INTEGER                           :: JLOOP
INTEGER, DIMENSION(:), ALLOCATABLE:: IFIELD
INTEGER                           :: ILENG
TYPE(FMHEADER)                    :: TZFMH

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
ILENG=LEN(HFIELD)
ALLOCATE(IFIELD(ILENG))
!
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN  ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,ILENG,IFIELD,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE ! parallel execution
    IF (ISP == TZFD%OWNER)  THEN
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.FALSE.,ILENG,IFIELD,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !
    CALL BCAST_HEADER(TZFD,TZFMH)
    !
    CALL MPI_BCAST(IFIELD,ILENG,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,&
       & IERR)
  END IF ! parallel execution
  !
  DO JLOOP=1,ILENG
    HFIELD(JLOOP:JLOOP)=ACHAR(IFIELD(JLOOP))
  END DO
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADC0_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
IF (ALLOCATED(IFIELD)) DEALLOCATE(IFIELD)
KRESP = IRESP
RETURN

END SUBROUTINE FMREADC0_ll

SUBROUTINE FMREADT0_ll(HFILEM,HRECFM,HFIPRI,HDIR,TFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
!*      0.    DECLARATIONS
!             ------------
!
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC 
USE MODD_TYPE_DATE
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),          INTENT(IN) ::HRECFM   ! name of the article to read
CHARACTER(LEN=*),          INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),          INTENT(IN) ::HDIR     ! Field form
TYPE (DATE_TIME),          INTENT(INOUT)::TFIELD ! array containing the data field
INTEGER,                   INTENT(INOUT)::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(INOUT)::KLENCH   ! length of comment string
CHARACTER(LEN=*),          INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,                   INTENT(INOUT)::KRESP    ! return-code
!
!
!*      0.2   Declarations of local variables
!
!-------------------------------------------------------------------------------


CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
INTEGER,DIMENSION(3)         :: ITDATE
REAL                         :: ZTIME
TYPE(FMHEADER)               :: TZFMH

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0

YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'

TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    CALL FM_READ_ll(TZFD%FLU,TRIM(HRECFM)//'%TDATE',.FALSE.,3,ITDATE&
         & ,TZFMH,IRESP)
    CALL FM_READ_ll(TZFD%FLU,TRIM(HRECFM)//'%TIME',.TRUE.,1,ZTIME&
         & ,TZFMH,IRESP)
    IF (IRESP /= 0) GOTO 1000
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL FM_READ_ll(TZFD%FLU,TRIM(HRECFM)//'%TDATE',.FALSE.,3,ITDATE&
           & ,TZFMH,IRESP)
      CALL FM_READ_ll(TZFD%FLU,TRIM(HRECFM)//'%TIME',.TRUE.,1,ZTIME&
           & ,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    ! Last header is significant
    CALL BCAST_HEADER(TZFD,TZFMH)      
    !
    CALL MPI_BCAST(ITDATE,3,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    CALL MPI_BCAST(ZTIME,1,MPI_FLOAT,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
  TFIELD%TDATE = DATE(ITDATE(1),ITDATE(2),ITDATE(3))
  TFIELD%TIME = ZTIME
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREADT0_ll",HFILEM,HFIPRI,HRECFM,HDIR,IRESP)
ENDIF
KRESP = IRESP
RETURN

END SUBROUTINE FMREADT0_ll

SUBROUTINE FMREAD_LB(HFILEM,HRECFM,HFIPRI,HLBTYPE,PLB,KRIM,KL3D,&
     & KGRID,KLENCH,HCOMMENT,KRESP)
USE MODD_FM
USE MODD_IO_ll,        ONLY : ISP,ISNPROC,GSMONOPROC,LPACK,L2D 
USE MODD_PARAMETERS_ll,ONLY : JPHEXT
USE MODE_DISTRIB_LB
USE MODE_TOOLS_ll,     ONLY : GET_GLOBALDIMS_ll
USE MODE_FD_ll,        ONLY : GETFD,JPFINL,FD_LL

CHARACTER(LEN=*),     INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*),     INTENT(IN) ::HRECFM   ! name of the article to be written
CHARACTER(LEN=*),     INTENT(IN) ::HFIPRI   ! file for prints
CHARACTER(LEN=*),     INTENT(IN) ::HLBTYPE  ! 'LBX','LBXU','LBY' or 'LBYV'
REAL, DIMENSION(:,:,:),TARGET, INTENT(INOUT)::PLB ! array containing the LB field
INTEGER,              INTENT(IN) :: KRIM  ! size of the LB area
INTEGER,              INTENT(IN) :: KL3D  ! size of the LB array in FM
INTEGER,              INTENT(INOUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,              INTENT(INOUT)::KLENCH ! length of comment string
CHARACTER(LEN=*),     INTENT(INOUT)::HCOMMENT ! comment string
INTEGER,              INTENT(INOUT)::KRESP  ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: Z3D
REAL,DIMENSION(:,:,:), POINTER           :: TX3DP
TYPE(FMHEADER)               :: TZFMH
INTEGER :: IIMAX_ll,IJMAX_ll
INTEGER :: IIB,IIE,IJB,IJE
INTEGER :: JI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: STATUS

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    IF (HLBTYPE == 'LBX' .OR. HLBTYPE == 'LBXU') THEN 
      ALLOCATE(Z3D(KL3D,SIZE(PLB,2),SIZE(PLB,3)))
      Z3D = 0.0
      IF (LPACK .AND. L2D) THEN
        TX3DP=>Z3D(:,2:2,:)
        CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(TX3DP),TX3DP,TZFMH,IRESP)
        Z3D(:,:,:) = SPREAD(Z3D(:,2,:),DIM=2,NCOPIES=3)
      ELSE
        CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(Z3D),Z3D,TZFMH,IRESP)
      END IF
      PLB(1:KRIM+1,:,:)          = Z3D(1:KRIM+1,:,:)
      PLB(KRIM+2:2*(KRIM+1),:,:) = Z3D(KL3D-KRIM:KL3D,:,:)
    ELSE !(HLBTYPE == 'LBY' .OR. HLBTYPE == 'LBYV') 
      ALLOCATE(Z3D(SIZE(PLB,1),KL3D,SIZE(PLB,3)))
      Z3D = 0.0
      CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(Z3D),Z3D,TZFMH,IRESP)
      PLB(:,1:KRIM+1,:)          = Z3D(:,1:KRIM+1,:)
      PLB(:,KRIM+2:2*(KRIM+1),:) = Z3D(:,KL3D-KRIM:KL3D,:)
    END IF
    IF (IRESP /= 0) GOTO 1000
  ELSE                 ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
      IF (HLBTYPE == 'LBX' .OR. HLBTYPE == 'LBXU') THEN 
        ALLOCATE(Z3D(KL3D,IJMAX_ll+2*JPHEXT,SIZE(PLB,3)))
        Z3D = 0.0
        IF (LPACK .AND. L2D) THEN
          TX3DP=>Z3D(:,2:2,:)
          CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(TX3DP),TX3DP,TZFMH,IRESP)
          Z3D(:,:,:) = SPREAD(Z3D(:,2,:),DIM=2,NCOPIES=3)
        ELSE
          CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(Z3D),Z3D,TZFMH,IRESP)
        END IF
        ! erase gap in LB field
        Z3D(KRIM+2:2*(KRIM+1),:,:) = Z3D(KL3D-KRIM:KL3D,:,:)
      ELSE !(HLBTYPE == 'LBY' .OR. HLBTYPE == 'LBYV') 
        ALLOCATE(Z3D(IIMAX_ll+2*JPHEXT,KL3D,SIZE(PLB,3)))
        Z3D = 0.0
        CALL FM_READ_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(Z3D),Z3D,TZFMH,IRESP)
        ! erase gap in LB field
        Z3D(:,KRIM+2:2*(KRIM+1),:) = Z3D(:,KL3D-KRIM:KL3D,:)
      END IF
    END IF
    !  
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
    IF (IRESP /= 0) GOTO 1000
    !  
    CALL BCAST_HEADER(TZFD,TZFMH)
    ! 
    IF (ISP == TZFD%OWNER)  THEN
      DO JI = 1,ISNPROC
        CALL GET_DISTRIB_LB(HLBTYPE,JI,'FM','READ',KRIM,IIB,IIE,IJB,IJE)
        IF (IIB /= 0) THEN
          TX3DP=>Z3D(IIB:IIE,IJB:IJE,:)
          IF (ISP /= JI) THEN 
            CALL MPI_SEND(TX3DP,SIZE(TX3DP),MPI_FLOAT,JI-1,99,TZFD%COMM,IERR)
          ELSE
            CALL GET_DISTRIB_LB(HLBTYPE,JI,'LOC','READ',KRIM,IIB,IIE,IJB,IJE)
            PLB(IIB:IIE,IJB:IJE,:) = TX3DP(:,:,:)
          END IF
        END IF
      END DO
    ELSE
      CALL GET_DISTRIB_LB(HLBTYPE,ISP,'LOC','READ',KRIM,IIB,IIE,IJB,IJE)
      IF (IIB /= 0) THEN
        TX3DP=>PLB(IIB:IIE,IJB:IJE,:)
        CALL MPI_RECV(TX3DP,SIZE(TX3DP),MPI_FLOAT,TZFD%OWNER-1,99,TZFD%COMM,STATUS,IERR)
      END IF
    END IF
  END IF !(GSMONOPROC)
  KGRID  = TZFMH%GRID
  KLENCH = TZFMH%COMLEN
  HCOMMENT = TZFMH%COMMENT(1:TZFMH%COMLEN)
ELSE 
  IRESP = -61          
END IF
!----------------------------------------------------------------
1000 CONTINUE
!! Error handler
IF (IRESP.NE.0) THEN
  CALL FM_READ_ERR("FMREAD_LB",HFILEM,HFIPRI,HRECFM,HLBTYPE,IRESP)
ENDIF
!
IF (ALLOCATED(Z3D)) DEALLOCATE (Z3D)
KRESP = IRESP
!
RETURN
END SUBROUTINE FMREAD_LB

END MODULE MODE_FMREAD


!
