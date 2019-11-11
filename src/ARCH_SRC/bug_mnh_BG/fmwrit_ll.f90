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

MODULE MODE_FMWRIT

USE MODD_MPIF

IMPLICIT NONE 

PRIVATE

INTERFACE FMWRIT
  MODULE PROCEDURE FMWRITX0_ll,FMWRITX1_ll,FMWRITX2_ll,FMWRITX3_ll,&
       & FMWRITX4_ll,FMWRITX5_ll,FMWRITX6_ll,&
       & FMWRITN0_ll,FMWRITN1_ll,FMWRITN2_ll,&
       & FMWRITL0_ll,FMWRITL1_ll,FMWRITC0_ll,FMWRITT0_ll
END INTERFACE

INTERFACE FMWRITBOX
  MODULE PROCEDURE FMWRITBOXX2_ll,FMWRITBOXX3_ll,FMWRITBOXX4_ll,&
       & FMWRITBOXX5_ll,FMWRITBOXX6_ll
END INTERFACE

PUBLIC FMWRIT_LB,FMWRITBOX,FMWRIT,FMWRITX0_ll,FMWRITX1_ll,FMWRITX2_ll,FMWRITX3_ll,&
       & FMWRITX4_ll,FMWRITX5_ll,FMWRITX6_ll,FMWRITN0_ll,FMWRITN1_ll,FMWRITN2_ll,&
       & FMWRITL0_ll,FMWRITL1_ll,FMWRITC0_ll,FMWRITT0_ll,FMWRITBOXX2_ll,&
       & FMWRITBOXX3_ll,FMWRITBOXX4_ll,FMWRITBOXX5_ll,FMWRITBOXX6_ll

!INCLUDE 'mpif.h'

CONTAINS 

SUBROUTINE FM_WRIT_ERR(HFUNC,HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH&
     & ,KRESP)
USE MODE_FM, ONLY : FMLOOK_ll

CHARACTER(LEN=*) :: HFUNC 
CHARACTER(LEN=*) :: HFILEM
CHARACTER(LEN=*) :: HFIPRI
CHARACTER(LEN=*) :: HRECFM
CHARACTER(LEN=*) :: HDIR
INTEGER          :: KGRID
INTEGER          :: KLENCH
INTEGER          :: KRESP

INTEGER          :: ILUPRI
INTEGER          :: IRESP

CALL FMLOOK_ll(HFIPRI,HFIPRI,ILUPRI,IRESP)
WRITE (ILUPRI,*) ' exit from ',HFUNC,' with RESP:',KRESP
WRITE (ILUPRI,*) '   | HFILEM = ',HFILEM
WRITE (ILUPRI,*) '   | HRECFM = ',HRECFM
WRITE (ILUPRI,*) '   | HDIR  = ',HDIR
WRITE (ILUPRI,*) '   | KGRID  = ',KGRID
WRITE (ILUPRI,*) '   | KLENCH = ',KLENCH

END SUBROUTINE FM_WRIT_ERR



SUBROUTINE FMWRITX0_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
!
!*      0.    DECLARATIONS
!             ------------
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN) ::HFILEM  ! FM-file name
CHARACTER(LEN=*),        INTENT(IN) ::HRECFM  ! name of the article to write
CHARACTER(LEN=*),        INTENT(IN) ::HFIPRI  ! output file for error messages
CHARACTER(LEN=*),        INTENT(IN) ::HDIR    ! field form
REAL,                    INTENT(IN) ::PFIELD  ! array containing the data field
INTEGER,                 INTENT(IN) ::KGRID   ! C-grid indicator (u,v,w,T)
INTEGER,                 INTENT(IN) ::KLENCH  ! length of comment string
CHARACTER(LEN=*),        INTENT(IN) ::HCOMMENT! comment string
INTEGER,                 INTENT(OUT)::KRESP   ! return-code 
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
!
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,1,PFIELD,TZFMH,IRESP)    
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,1,PFIELD,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITX0_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
KRESP = IRESP
END SUBROUTINE FMWRITX0_ll
  
SUBROUTINE FMWRITX1_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!
!*      0.    DECLARATIONS
!             ------------
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),        INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),        INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),        INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                 INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                 INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),        INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                 INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
!----------------------------------------------------------------
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
TYPE(FMHEADER)               :: TZFMH
REAL,DIMENSION(:),POINTER    :: ZFIELDP
LOGICAL                      :: GALLOC
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------    
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
    ELSE
      ALLOCATE(ZFIELDP(0))
      GALLOC = .TRUE.
    END IF
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      CALL GATHER_XXFIELD(HDIR,PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITX1_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITX1_ll
  
SUBROUTINE FMWRITX2_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),          INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),          INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),          INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:,:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),          INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                   INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                  :: YFNLFI
INTEGER                                :: IERR
TYPE(FD_ll), POINTER                   :: TZFD
INTEGER                                :: IRESP
REAL,DIMENSION(:,:),POINTER            :: ZFIELDP
TYPE(FMHEADER)                         :: TZFMH
LOGICAL                                :: GALLOC
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
      ZFIELDP=>PFIELD(:,2:2)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
    ELSE
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      ! I/O processor case
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
    ELSE
      ALLOCATE(ZFIELDP(0,0))
      GALLOC = .TRUE.
    END IF
    !    
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      CALL GATHER_XXFIELD(HDIR,PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        CALL GATHER_XXFIELD('XX',PFIELD(:,2),ZFIELDP(:,1),TZFD%OWNER,TZFD%COMM)
      ELSE
        CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
      END IF
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
  END IF
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITX2_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITX2_ll
  
SUBROUTINE FMWRITX3_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),            INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),            INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),            INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),            INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:,:,:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                     INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                     INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),            INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                     INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                    :: YFNLFI
INTEGER                                  :: IERR
TYPE(FD_ll), POINTER                     :: TZFD
INTEGER                                  :: IRESP
REAL,DIMENSION(:,:,:),POINTER            :: ZFIELDP
TYPE(FMHEADER)                           :: TZFMH
LOGICAL                                  :: GALLOC
!JUAN
INTEGER                                  :: JK
CHARACTER(LEN=LEN(HRECFM))               :: YK,YRECZSLIDE
REAL,DIMENSION(:,:),POINTER              :: ZSLIDE_ll,ZSLIDE
!JUAN
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2,:)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
      ZFIELDP=>PFIELD(:,2:2,:)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
    ELSE
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF
  ELSE ! multiprocessor execution
!
!JUAN BG Z SLIDE 
!
    IF (ISP == TZFD%OWNER)  THEN
!JUAN      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
      ZSLIDE => PFIELD(:,:,1)
      CALL ALLOCBUFFER_ll(ZSLIDE_ll,ZSLIDE,HDIR,GALLOC)
    ELSE
!JUAN       ALLOCATE(ZFIELDP(0,0,0))
      ALLOCATE(ZSLIDE_ll(0,0))
      GALLOC = .TRUE.
    END IF
    !
Z_SLIDE: DO JK=1,SIZE(PFIELD,3)
!JUAN
    WRITE(YK,'(I4.4)')  JK
    YRECZSLIDE = TRIM(HRECFM)//YK
!JUAN
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
       STOP " XX NON PREVU SUR BG POUR LE MOMENT "
      CALL GATHER_XXFIELD(HDIR,PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
       STOP " L2D NON PREVU SUR BG POUR LE MOMENT "
        CALL GATHER_XXFIELD('XX',PFIELD(:,2,:),ZFIELDP(:,1,:),TZFD%OWNER,TZFD%COMM)
      ELSE
!JUAN         IF ( JK == 1 ) CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
        ZSLIDE => PFIELD(:,:,JK)
        CALL GATHER_XYFIELD(ZSLIDE,ZSLIDE_ll,TZFD%OWNER,TZFD%COMM)
      END IF
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
!JUAN        IF ( JK == 1 ) CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
!JUAN            & ,IRESP)
      CALL FM_WRIT_ll(TZFD%FLU,YRECZSLIDE,.TRUE.,SIZE(ZSLIDE_ll),ZSLIDE_ll,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
 END DO Z_SLIDE
!JUAN BG Z SLIDE  
  END IF ! multiprocessor execution
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITX3_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
!JUAN IF (GALLOC) DEALLOCATE(ZFIELDP)
IF (GALLOC) DEALLOCATE(ZSLIDE_ll)
KRESP = IRESP
END SUBROUTINE FMWRITX3_ll
  
SUBROUTINE FMWRITX4_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),              INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),              INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),              INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),              INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:,:,:,:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                       INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                       INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),              INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                       INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                    :: YFNLFI
INTEGER                                  :: IERR
TYPE(FD_ll), POINTER                     :: TZFD
INTEGER                                  :: IRESP
REAL,DIMENSION(:,:,:,:),POINTER          :: ZFIELDP
TYPE(FMHEADER)                           :: TZFMH
LOGICAL                                  :: GALLOC
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2,:,:)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D  .AND. SIZE(PFIELD,2)==3) THEN
      ZFIELDP=>PFIELD(:,2:2,:,:)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
    ELSE
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      CALL GATHER_XXFIELD(HDIR,PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        CALL GATHER_XXFIELD('XX',PFIELD(:,2,:,:),ZFIELDP(:,1,:,:),TZFD%OWNER,TZFD%COMM)
      ELSE
        CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
      END IF
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF ! multiprocessor execution
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITX4_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITX4_ll

SUBROUTINE FMWRITX5_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),                INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),                INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),                INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),                INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:,:,:,:,:),TARGET,INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                         INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                         INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),                INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                         INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                    :: YFNLFI
INTEGER                                  :: IERR
TYPE(FD_ll), POINTER                     :: TZFD
INTEGER                                  :: IRESP
REAL,DIMENSION(:,:,:,:,:),POINTER        :: ZFIELDP
TYPE(FMHEADER)                           :: TZFMH
LOGICAL                                  :: GALLOC
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(PFIELD,1)==3 .AND. SIZE(PFIELD,2)==3) THEN 
      ZFIELDP=>PFIELD(2:2,2:2,:,:,:)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(PFIELD,2)==3) THEN
      ZFIELDP=>PFIELD(:,2:2,:,:,:)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
    ELSE
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
    END IF
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      CALL GATHER_XXFIELD(HDIR,PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        CALL GATHER_XXFIELD('XX',PFIELD(:,2,:,:,:),ZFIELDP(:,1,:,:,:),&
             & TZFD%OWNER,TZFD%COMM)
      ELSE
        CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
      END IF
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF ! multiprocessor execution
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITX5_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITX5_ll

SUBROUTINE FMWRITX6_ll(HFILEM,HRECFM,HFIPRI,HDIR,PFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),                INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),                INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),                INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),                INTENT(IN) ::HDIR     ! field form
REAL,DIMENSION(:,:,:,:,:,:),TARGET,INTENT(IN) ::PFIELD ! array containing the data field
INTEGER,                         INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                         INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),                INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                         INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                    :: YFNLFI
INTEGER                                  :: IERR
TYPE(FD_ll), POINTER                     :: TZFD
INTEGER                                  :: IRESP
REAL,DIMENSION(:,:,:,:,:,:),POINTER        :: ZFIELDP
TYPE(FMHEADER)                           :: TZFMH
LOGICAL                                  :: GALLOC
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PFIELD),PFIELD,TZFMH,IRESP)
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(ZFIELDP,PFIELD,HDIR,GALLOC)
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      CALL GATHER_XXFIELD(HDIR,PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
    ELSE IF (HDIR == 'XY') THEN
      CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM)
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF ! multiprocessor execution
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITX6_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITX6_ll

SUBROUTINE FMWRITN0_ll(HFILEM,HRECFM,HFIPRI,HDIR,KFIELD,KGRID,&
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
CHARACTER(LEN=*),   INTENT(IN) ::HFILEM  ! FM-file name
CHARACTER(LEN=*),   INTENT(IN) ::HRECFM  ! name of the article to read
CHARACTER(LEN=*),   INTENT(IN) ::HFIPRI  ! output file for error messages
CHARACTER(LEN=*),   INTENT(IN) ::HDIR    ! field form
INTEGER,            INTENT(IN) ::KFIELD  ! array containing the data field
INTEGER,            INTENT(IN) ::KGRID   ! C-grid indicator (u,v,w,T)
INTEGER,            INTENT(IN) ::KLENCH  ! length of comment string
CHARACTER(LEN=*),   INTENT(IN) ::HCOMMENT! comment string
INTEGER,            INTENT(OUT)::KRESP   ! return-code
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
TYPE(FMHEADER)               :: TZFMH

!----------------------------------------------------------------
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,1,KFIELD,TZFMH,IRESP)
  ELSE 
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,1,KFIELD,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF ! multiprocessor execution
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITN0_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH&
       & ,IRESP)
END IF
KRESP = IRESP
END SUBROUTINE FMWRITN0_ll

SUBROUTINE FMWRITN1_ll(HFILEM,HRECFM,HFIPRI,HDIR,KFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!*      0.    DECLARATIONS
!             ------------
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),           INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),           INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),           INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),           INTENT(IN) ::HDIR     ! field form
INTEGER,DIMENSION(:),TARGET,INTENT(IN) ::KFIELD   ! array containing the data field
INTEGER,                    INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                    INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),           INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                    INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
TYPE(FMHEADER)               :: TZFMH
INTEGER,DIMENSION(:),POINTER :: IFIELDP
LOGICAL                      :: GALLOC
!----------------------------------------------------------------
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(KFIELD),KFIELD,TZFMH,IRESP)
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(IFIELDP,KFIELD,HDIR,GALLOC)
    ELSE
      ALLOCATE(IFIELDP(0))
      GALLOC = .TRUE.
    END IF
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      CALL GATHER_XXFIELD(HDIR,KFIELD,IFIELDP,TZFD%OWNER,TZFD%COMM)
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITN1_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH&
       & ,IRESP)
END IF
IF (GALLOC) DEALLOCATE(IFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITN1_ll

SUBROUTINE FMWRITN2_ll(HFILEM,HRECFM,HFIPRI,HDIR,KFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC,LPACK,L1D,L2D
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_ALLOCBUFFER_ll
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),             INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),             INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),             INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),             INTENT(IN) ::HDIR     ! field form
INTEGER,DIMENSION(:,:),TARGET,INTENT(IN) ::KFIELD ! array containing the data field
INTEGER,                      INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
INTEGER,                      INTENT(IN) ::KLENCH   ! length of comment string
CHARACTER(LEN=*),             INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                      INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                    :: YFNLFI
INTEGER                                  :: IERR
TYPE(FD_ll), POINTER                     :: TZFD
INTEGER                                  :: IRESP
INTEGER,DIMENSION(:,:),POINTER           :: IFIELDP
TYPE(FMHEADER)                           :: TZFMH
LOGICAL                                  :: GALLOC

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
!    IF (LPACK .AND. L1D .AND. HDIR=='XY') THEN 
    IF (LPACK .AND. L1D .AND. SIZE(KFIELD,1)==3 .AND. SIZE(KFIELD,2)==3) THEN 
      IFIELDP=>KFIELD(2:2,2:2)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP,TZFMH,IRESP)
!    ELSE IF (LPACK .AND. L2D .AND. HDIR=='XY') THEN
    ELSE IF (LPACK .AND. L2D .AND. SIZE(KFIELD,2)==3) THEN
      IFIELDP=>KFIELD(:,2:2)
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP,TZFMH,IRESP)
    ELSE
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(KFIELD),KFIELD,TZFMH,IRESP)
    END IF
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      CALL ALLOCBUFFER_ll(IFIELDP,KFIELD,HDIR,GALLOC)
    ELSE
      ALLOCATE(IFIELDP(0,0))
      GALLOC = .TRUE.
    END IF
    !
    IF (HDIR == 'XX' .OR. HDIR =='YY') THEN
      CALL GATHER_XXFIELD(HDIR,KFIELD,IFIELDP,TZFD%OWNER,TZFD%COMM)
    ELSE IF (HDIR == 'XY') THEN
      IF (LPACK .AND. L2D) THEN
        CALL GATHER_XXFIELD('XX',KFIELD(:,2),IFIELDP(:,1),TZFD%OWNER,TZFD%COMM)
      ELSE
        CALL GATHER_XYFIELD(KFIELD,IFIELDP,TZFD%OWNER,TZFD%COMM)
      END IF
    END IF
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELDP),IFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF

ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITN2_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH,IRESP)
END IF
IF (GALLOC) DEALLOCATE(IFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITN2_ll
  
  
SUBROUTINE FMWRITL0_ll(HFILEM,HRECFM,HFIPRI,HDIR,OFIELD,KGRID,&
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
CHARACTER(LEN=*), INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*), INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*), INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*), INTENT(IN) ::HDIR   ! field form
LOGICAL,          INTENT(IN) ::OFIELD ! array containing the data field
INTEGER,          INTENT(IN)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,          INTENT(IN)::KLENCH ! length of comment string
CHARACTER(LEN=*), INTENT(IN)::HCOMMENT ! comment string
INTEGER,          INTENT(OUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
INTEGER                      :: IFIELD
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
TYPE(FMHEADER)               :: TZFMH

!----------------------------------------------------------------
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
IF (OFIELD) THEN
  IFIELD=1
ELSE
  IFIELD=0
END IF
!----------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,1,IFIELD,TZFMH,IRESP)
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,1,IFIELD,TZFMH,IRESP)
    END IF
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITL0_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH&
       & ,IRESP)
END IF
KRESP = IRESP
END SUBROUTINE FMWRITL0_ll

SUBROUTINE FMWRITL1_ll(HFILEM,HRECFM,HFIPRI,HDIR,OFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL

!*      0.    DECLARATIONS
!             ------------
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),    INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),    INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),    INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),    INTENT(IN) ::HDIR   ! field form
LOGICAL,DIMENSION(:),INTENT(IN) ::OFIELD ! array containing the data field
INTEGER,             INTENT(IN)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,             INTENT(IN)::KLENCH ! length of comment string
CHARACTER(LEN=*),    INTENT(IN)::HCOMMENT ! comment string
INTEGER,             INTENT(OUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
INTEGER, DIMENSION(SIZE(OFIELD)) :: IFIELD
CHARACTER(LEN=JPFINL)            :: YFNLFI
INTEGER                          :: IERR
TYPE(FD_ll), POINTER             :: TZFD
INTEGER                          :: IRESP
TYPE(FMHEADER)                   :: TZFMH

!----------------------------------------------------------------
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
WHERE (OFIELD)
  IFIELD=1
ELSEWHERE
  IFIELD=0
END WHERE
!----------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMLEN=KLENCH
    TZFMH%COMMENT=HCOMMENT
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELD),IFIELD,TZFMH,IRESP)
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,SIZE(IFIELD),IFIELD,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITL1_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH&
       & ,IRESP)
END IF
KRESP = IRESP
END SUBROUTINE FMWRITL1_ll

SUBROUTINE FMWRITC0_ll(HFILEM,HRECFM,HFIPRI,HDIR,HFIELD,KGRID,&
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
CHARACTER(LEN=*),  INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),  INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),  INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),  INTENT(IN) ::HDIR   ! field form
CHARACTER(LEN=*),  INTENT(IN) ::HFIELD ! array containing the data field
INTEGER,           INTENT(IN)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,           INTENT(IN)::KLENCH ! length of comment string
CHARACTER(LEN=*),  INTENT(IN)::HCOMMENT ! comment string
INTEGER,           INTENT(OUT)::KRESP    ! return-code
!
!*      0.2   Declarations of local variables
!
INTEGER                          :: JLOOP
INTEGER,DIMENSION(:),ALLOCATABLE :: IFIELD
INTEGER                          :: ILENG
CHARACTER(LEN=JPFINL)            :: YFNLFI
INTEGER                          :: IERR
TYPE(FD_ll), POINTER             :: TZFD
INTEGER                          :: IRESP
TYPE(FMHEADER)                   :: TZFMH

!----------------------------------------------------------------
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
ILENG=LEN(HFIELD)
!
IF (ILENG==0) THEN
  ILENG=1
  ALLOCATE(IFIELD(1))
  IFIELD(1)=IACHAR(' ')
ELSE
  ALLOCATE(IFIELD(ILENG))
  DO JLOOP=1,ILENG
    IFIELD(JLOOP)=IACHAR(HFIELD(JLOOP:JLOOP))
  END DO
END IF
!----------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN  ! sequential execution
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,ILENG,IFIELD,TZFMH,KRESP)
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.FALSE.,ILENG,IFIELD,TZFMH,KRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITC0_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH&
       & ,IRESP)
END IF
IF (ALLOCATED(IFIELD)) DEALLOCATE(IFIELD)
KRESP = IRESP
END SUBROUTINE FMWRITC0_ll

SUBROUTINE FMWRITT0_ll(HFILEM,HRECFM,HFIPRI,HDIR,TFIELD,KGRID,&
     KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_TYPE_DATE
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),    INTENT(IN) ::HFILEM ! FM-file name
CHARACTER(LEN=*),    INTENT(IN) ::HRECFM ! name of the article to read
CHARACTER(LEN=*),    INTENT(IN) ::HFIPRI ! output file for error messages
CHARACTER(LEN=*),    INTENT(IN) ::HDIR   ! field form
TYPE (DATE_TIME),    INTENT(IN) ::TFIELD ! array containing the data field
INTEGER,             INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,             INTENT(IN) ::KLENCH ! length of comment string
CHARACTER(LEN=*),    INTENT(IN) ::HCOMMENT ! comment string
INTEGER,             INTENT(OUT)::KRESP    ! return-code
!--------------------------------------------------------------------
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)        :: YFNLFI
INTEGER                      :: IERR
TYPE(FD_ll), POINTER         :: TZFD
INTEGER                      :: IRESP
TYPE(FMHEADER)               :: TZFMH
INTEGER, DIMENSION(3)        :: ITDATE    ! date array
!
!-------------------------------------------------------------------------------
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
ITDATE(1)=TFIELD%TDATE%YEAR
ITDATE(2)=TFIELD%TDATE%MONTH
ITDATE(3)=TFIELD%TDATE%DAY
!-------------------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID=KGRID
    TZFMH%COMMENT='YYYYMMDD'
    TZFMH%COMLEN=LEN_TRIM(TZFMH%COMMENT)
    CALL FM_WRIT_ll(TZFD%FLU,TRIM(HRECFM)//'%TDATE',.FALSE.,3,ITDATE&
         & ,TZFMH,IRESP)
    TZFMH%COMMENT='SECONDS'
    TZFMH%COMLEN=LEN_TRIM(TZFMH%COMMENT)
    CALL FM_WRIT_ll(TZFD%FLU,TRIM(HRECFM)//'%TIME',.TRUE.,1,TFIELD%TIME&
         & ,TZFMH,IRESP)
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID=KGRID
      TZFMH%COMMENT='YYYYMMDD'
      TZFMH%COMLEN=LEN_TRIM(TZFMH%COMMENT)
      CALL FM_WRIT_ll(TZFD%FLU,TRIM(HRECFM)//'%TDATE',.FALSE.,3,ITDATE&
           & ,TZFMH,IRESP)
      TZFMH%COMMENT='SECONDS'
      TZFMH%COMLEN=LEN_TRIM(TZFMH%COMMENT)
      CALL FM_WRIT_ll(TZFD%FLU,TRIM(HRECFM)//'%TIME',.TRUE.,1,TFIELD%TIME&
           & ,TZFMH,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD%COMM,IERR)
  END IF
ELSE 
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITT0_ll",HFILEM,HFIPRI,HRECFM,HDIR,KGRID,KLENCH&
       & ,IRESP)
END IF
KRESP = IRESP
END SUBROUTINE FMWRITT0_ll

SUBROUTINE FMWRIT_LB(HFILEM,HRECFM,HFIPRI,HLBTYPE,PLB,KRIM,KL3D,&
     & KGRID,KLENCH,HCOMMENT,KRESP)
USE MODD_IO_ll,        ONLY : ISP,ISNPROC,GSMONOPROC,LPACK,L2D
USE MODD_PARAMETERS_ll,ONLY : JPHEXT
USE MODD_FM
USE MODE_DISTRIB_LB
USE MODE_TOOLS_ll,     ONLY : GET_GLOBALDIMS_ll
USE MODE_FD_ll,        ONLY : GETFD,JPFINL,FD_LL
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),       INTENT(IN) ::HFILEM ! file name
CHARACTER(LEN=*),       INTENT(IN) ::HRECFM ! name of the article to be written
CHARACTER(LEN=*),       INTENT(IN) ::HFIPRI ! file for prints in FM
CHARACTER(LEN=*),       INTENT(IN) ::HLBTYPE! 'LBX','LBXU','LBY' or 'LBYV'
REAL,DIMENSION(:,:,:),TARGET,INTENT(IN) ::PLB ! array containing the LB field
INTEGER,                INTENT(IN) ::KRIM  ! size of the LB area
INTEGER,                INTENT(IN) ::KL3D  ! size of the LB array in FM
INTEGER,                INTENT(IN) ::KGRID ! C-grid indicator (u,v,w,T)
INTEGER,                INTENT(IN) ::KLENCH ! length of comment string
CHARACTER(LEN=*),       INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                INTENT(OUT)::KRESP  ! return-code
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)                    :: YFNLFI
INTEGER                                  :: IERR
TYPE(FD_ll), POINTER                     :: TZFD
INTEGER                                  :: IRESP
REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: Z3D
REAL,DIMENSION(:,:,:), POINTER           :: TX3DP
TYPE(FMHEADER)                           :: TZFMH
INTEGER                                  :: IIMAX_ll,IJMAX_ll
INTEGER                                  :: JI
INTEGER :: IIB,IIE,IJB,IJE
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: STATUS
!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
IF (KL3D /= 2*(KRIM+1)) THEN
  IRESP = -30
  GOTO 1000
END IF
!
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN  ! sequential execution
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      IF (LPACK .AND. L2D) THEN
        TX3DP=>PLB(:,2:2,:)
        CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(TX3DP),TX3DP,TZFMH,IRESP)
      ELSE
        CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(PLB),PLB,TZFMH,IRESP)
      END IF
  ELSE
    IF (ISP == TZFD%OWNER)  THEN
      ! I/O proc case
      CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
      IF (HLBTYPE == 'LBX' .OR. HLBTYPE == 'LBXU') THEN 
        ALLOCATE(Z3D((KRIM+1)*2,IJMAX_ll+2*JPHEXT,SIZE(PLB,3)))
      ELSE ! HLBTYPE == 'LBY' .OR. HLBTYPE == 'LBYV' 
        ALLOCATE(Z3D(IIMAX_ll+2*JPHEXT,(KRIM+1)*2,SIZE(PLB,3)))
      END IF
      DO JI = 1,ISNPROC
        CALL GET_DISTRIB_LB(HLBTYPE,JI,'FM','WRITE',KRIM,IIB,IIE,IJB,IJE)
        IF (IIB /= 0) THEN
          TX3DP=>Z3D(IIB:IIE,IJB:IJE,:)
          IF (ISP /= JI) THEN
            CALL MPI_RECV(TX3DP,SIZE(TX3DP),MPI_FLOAT,JI-1,99,TZFD%COMM,STATUS,IERR) 
          ELSE
            CALL GET_DISTRIB_LB(HLBTYPE,JI,'LOC','WRITE',KRIM,IIB,IIE,IJB,IJE)
            TX3DP = PLB(IIB:IIE,IJB:IJE,:)
          END IF
        END IF
      END DO
      TZFMH%GRID=KGRID
      TZFMH%COMLEN=KLENCH
      TZFMH%COMMENT=HCOMMENT
      IF (LPACK .AND. L2D) THEN
        TX3DP=>Z3D(:,2:2,:)
      ELSE
        TX3DP=>Z3D
      END IF
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(TX3DP),TX3DP,TZFMH,IRESP)
    ELSE
      ! Other processors
      CALL GET_DISTRIB_LB(HLBTYPE,ISP,'LOC','WRITE',KRIM,IIB,IIE,IJB,IJE)
      IF (IIB /= 0) THEN
        TX3DP=>PLB(IIB:IIE,IJB:IJE,:)
        CALL MPI_SEND(TX3DP,SIZE(TX3DP),MPI_FLOAT,TZFD%OWNER-1,99,TZFD%COMM,IERR)
      END IF
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
  END IF !(GSMONOPROC)
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
1000 CONTINUE
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRIT_LB",HFILEM,HFIPRI,HRECFM,HLBTYPE,KGRID,KLENCH,IRESP)
END IF
!
IF (ALLOCATED(Z3D)) DEALLOCATE(Z3D)
KRESP = IRESP
END SUBROUTINE FMWRIT_LB

SUBROUTINE FMWRITBOXX2_ll(HFILEM,HRECFM,HFIPRI,HBUDGET,PFIELD,KGRID,&
         HCOMMENT,KXOBOX,KXEBOX,KYOBOX,KYEBOX,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),            INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),            INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),            INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),            INTENT(IN) ::HBUDGET  ! 'BUDGET' (budget)  or 'OTHER' (MesoNH field)
REAL,DIMENSION(:,:),TARGET,  INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                     INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
CHARACTER(LEN=*),            INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                     INTENT(IN) ::KXOBOX   ! 
INTEGER,                     INTENT(IN) ::KXEBOX   ! Global coordinates of the box
INTEGER,                     INTENT(IN) ::KYOBOX   ! 
INTEGER,                     INTENT(IN) ::KYEBOX   ! 
INTEGER,                     INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)            :: YFNLFI
INTEGER                          :: IERR
TYPE(FD_ll), POINTER             :: TZFD
INTEGER                          :: IRESP
REAL,DIMENSION(:,:),POINTER      :: ZFIELDP
TYPE(FMHEADER)                   :: TZFMH
LOGICAL                          :: GALLOC

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID    = KGRID
    TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
    TZFMH%COMMENT = HCOMMENT
    IF (HBUDGET /= 'BUDGET') THEN
      ! take the sub-section of PFIELD defined by the box
      ZFIELDP=>PFIELD(KXOBOX:KXEBOX,KYOBOX:KYEBOX)
    ELSE
      ! take the field as a budget
      ZFIELDP=>PFIELD
    END IF
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      ! Allocate the box
      ALLOCATE(ZFIELDP(KXEBOX-KXOBOX+1,KYEBOX-KYOBOX+1))
      GALLOC = .TRUE.
    ELSE
      ALLOCATE(ZFIELDP(0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM,&
         & KXOBOX,KXEBOX,KYOBOX,KYEBOX,HBUDGET)
    !
    IF (ISP == TZFD%OWNER) THEN
      TZFMH%GRID    = KGRID
      TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
      TZFMH%COMMENT = HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
  END IF ! multiprocessor execution
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITBOXX2_ll",HFILEM,HFIPRI,HRECFM,'XY',KGRID,LEN(HCOMMENT),IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITBOXX2_ll

SUBROUTINE FMWRITBOXX3_ll(HFILEM,HRECFM,HFIPRI,HBUDGET,PFIELD,KGRID,&
         HCOMMENT,KXOBOX,KXEBOX,KYOBOX,KYEBOX,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),            INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),            INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),            INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),            INTENT(IN) ::HBUDGET  ! 'BUDGET' (budget)  or 'OTHER' (MesoNH field)
REAL,DIMENSION(:,:,:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                     INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
CHARACTER(LEN=*),            INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                     INTENT(IN) ::KXOBOX   ! 
INTEGER,                     INTENT(IN) ::KXEBOX   ! Global coordinates of the box
INTEGER,                     INTENT(IN) ::KYOBOX   ! 
INTEGER,                     INTENT(IN) ::KYEBOX   ! 
INTEGER,                     INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)               :: YFNLFI
INTEGER                             :: IERR
TYPE(FD_ll), POINTER                :: TZFD
INTEGER                             :: IRESP
REAL,DIMENSION(:,:,:),POINTER       :: ZFIELDP
TYPE(FMHEADER)                      :: TZFMH
LOGICAL                             :: GALLOC

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID    = KGRID
    TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
    TZFMH%COMMENT = HCOMMENT
    IF (HBUDGET /= 'BUDGET') THEN
      ! take the sub-section of PFIELD defined by the box
      ZFIELDP=>PFIELD(KXOBOX:KXEBOX,KYOBOX:KYEBOX,:)
    ELSE
      ! take the field as a budget
      ZFIELDP=>PFIELD
    END IF
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      ! Allocate the box
      ALLOCATE(ZFIELDP(KXEBOX-KXOBOX+1,KYEBOX-KYOBOX+1,SIZE(PFIELD,3)))
      GALLOC = .TRUE.
    ELSE
      ALLOCATE(ZFIELDP(0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM,&
         & KXOBOX,KXEBOX,KYOBOX,KYEBOX,HBUDGET)
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID    = KGRID
      TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
      TZFMH%COMMENT = HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
  END IF ! multiprocessor execution
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITBOXX3_ll",HFILEM,HFIPRI,HRECFM,'XY',KGRID,LEN(HCOMMENT),IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITBOXX3_ll

SUBROUTINE FMWRITBOXX4_ll(HFILEM,HRECFM,HFIPRI,HBUDGET,PFIELD,KGRID,&
         HCOMMENT,KXOBOX,KXEBOX,KYOBOX,KYEBOX,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),              INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),              INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),              INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),              INTENT(IN) ::HBUDGET  ! 'BUDGET' (budget)  or 'OTHER' (MesoNH field)
REAL,DIMENSION(:,:,:,:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                       INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
CHARACTER(LEN=*),              INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                       INTENT(IN) ::KXOBOX   ! 
INTEGER,                       INTENT(IN) ::KXEBOX   ! Global coordinates of the box
INTEGER,                       INTENT(IN) ::KYOBOX   ! 
INTEGER,                       INTENT(IN) ::KYEBOX   ! 
INTEGER,                       INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)               :: YFNLFI
INTEGER                             :: IERR
TYPE(FD_ll), POINTER                :: TZFD
INTEGER                             :: IRESP
REAL,DIMENSION(:,:,:,:),POINTER     :: ZFIELDP
TYPE(FMHEADER)                      :: TZFMH
LOGICAL                             :: GALLOC

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID    = KGRID
    TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
    TZFMH%COMMENT = HCOMMENT
    IF (HBUDGET /= 'BUDGET') THEN
      ! take the sub-section of PFIELD defined by the box
      ZFIELDP=>PFIELD(KXOBOX:KXEBOX,KYOBOX:KYEBOX,:,:)
    ELSE
      ! take the field as a budget
      ZFIELDP=>PFIELD
    END IF
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      ! Allocate the box
      ALLOCATE(ZFIELDP(KXEBOX-KXOBOX+1,KYEBOX-KYOBOX+1,SIZE(PFIELD,3),SIZE(PFIELD,4)))
      GALLOC = .TRUE.
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM,&
         & KXOBOX,KXEBOX,KYOBOX,KYEBOX,HBUDGET)
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID    = KGRID
      TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
      TZFMH%COMMENT = HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
  END IF ! multiprocessor execution
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITBOXX4_ll",HFILEM,HFIPRI,HRECFM,'XY',KGRID,LEN(HCOMMENT),IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITBOXX4_ll

SUBROUTINE FMWRITBOXX5_ll(HFILEM,HRECFM,HFIPRI,HBUDGET,PFIELD,KGRID,&
         HCOMMENT,KXOBOX,KXEBOX,KYOBOX,KYEBOX,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),              INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),              INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),              INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),              INTENT(IN) ::HBUDGET  ! 'BUDGET' (budget)  or 'OTHER' (MesoNH field)
REAL,DIMENSION(:,:,:,:,:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                       INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
CHARACTER(LEN=*),              INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                       INTENT(IN) ::KXOBOX   ! 
INTEGER,                       INTENT(IN) ::KXEBOX   ! Global coordinates of the box
INTEGER,                       INTENT(IN) ::KYOBOX   ! 
INTEGER,                       INTENT(IN) ::KYEBOX   ! 
INTEGER,                       INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)               :: YFNLFI
INTEGER                             :: IERR
TYPE(FD_ll), POINTER                :: TZFD
INTEGER                             :: IRESP
REAL,DIMENSION(:,:,:,:,:),POINTER   :: ZFIELDP
TYPE(FMHEADER)                      :: TZFMH
LOGICAL                             :: GALLOC

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID    = KGRID
    TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
    TZFMH%COMMENT = HCOMMENT
    IF (HBUDGET /= 'BUDGET') THEN
      ! take the sub-section of PFIELD defined by the box
      ZFIELDP=>PFIELD(KXOBOX:KXEBOX,KYOBOX:KYEBOX,:,:,:)
    ELSE
      ! take the field as a budget
      ZFIELDP=>PFIELD
    END IF
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      ! Allocate the box
      ALLOCATE(ZFIELDP(KXEBOX-KXOBOX+1,KYEBOX-KYOBOX+1,SIZE(PFIELD,3),&
           & SIZE(PFIELD,4),SIZE(PFIELD,5)))
      GALLOC = .TRUE.
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM,&
         & KXOBOX,KXEBOX,KYOBOX,KYEBOX,HBUDGET)
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID    = KGRID
      TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
      TZFMH%COMMENT = HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
  END IF ! multiprocessor execution
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITBOXX5_ll",HFILEM,HFIPRI,HRECFM,'XY',KGRID,LEN(HCOMMENT),IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITBOXX5_ll

SUBROUTINE FMWRITBOXX6_ll(HFILEM,HRECFM,HFIPRI,HBUDGET,PFIELD,KGRID,&
         HCOMMENT,KXOBOX,KXEBOX,KYOBOX,KYEBOX,KRESP)
USE MODD_IO_ll, ONLY : ISP,GSMONOPROC
USE MODD_FM
USE MODE_FD_ll, ONLY : GETFD,JPFINL,FD_LL
USE MODE_GATHER_ll
!
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*),              INTENT(IN) ::HFILEM   ! FM-file name
CHARACTER(LEN=*),              INTENT(IN) ::HRECFM   ! name of the article to write
CHARACTER(LEN=*),              INTENT(IN) ::HFIPRI   ! output file for error messages
CHARACTER(LEN=*),              INTENT(IN) ::HBUDGET  ! 'BUDGET' (budget)  or 'OTHER' (MesoNH field)
REAL,DIMENSION(:,:,:,:,:,:),TARGET,INTENT(IN) ::PFIELD   ! array containing the data field
INTEGER,                       INTENT(IN) ::KGRID    ! C-grid indicator (u,v,w,T)
CHARACTER(LEN=*),              INTENT(IN) ::HCOMMENT ! comment string
INTEGER,                       INTENT(IN) ::KXOBOX   ! 
INTEGER,                       INTENT(IN) ::KXEBOX   ! Global coordinates of the box
INTEGER,                       INTENT(IN) ::KYOBOX   ! 
INTEGER,                       INTENT(IN) ::KYEBOX   ! 
INTEGER,                       INTENT(OUT)::KRESP    ! return-code 
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=JPFINL)               :: YFNLFI
INTEGER                             :: IERR
TYPE(FD_ll), POINTER                :: TZFD
INTEGER                             :: IRESP
REAL,DIMENSION(:,:,:,:,:,:),POINTER :: ZFIELDP
TYPE(FMHEADER)                      :: TZFMH
LOGICAL                             :: GALLOC

!
!*      1.1   THE NAME OF LFIFM
!
IRESP = 0
GALLOC = .FALSE.
YFNLFI=TRIM(ADJUSTL(HFILEM))//'.lfi'
!------------------------------------------------------------------
TZFD=>GETFD(YFNLFI)
IF (ASSOCIATED(TZFD)) THEN
  IF (GSMONOPROC) THEN ! sequential execution
    TZFMH%GRID    = KGRID
    TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
    TZFMH%COMMENT = HCOMMENT
    IF (HBUDGET /= 'BUDGET') THEN
      ! take the sub-section of PFIELD defined by the box
      ZFIELDP=>PFIELD(KXOBOX:KXEBOX,KYOBOX:KYEBOX,:,:,:,:)
    ELSE
      ! take the field as a budget
      ZFIELDP=>PFIELD
    END IF
    CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH,IRESP)
  ELSE ! multiprocessor execution
    IF (ISP == TZFD%OWNER)  THEN
      ! Allocate the box
      ALLOCATE(ZFIELDP(KXEBOX-KXOBOX+1,KYEBOX-KYOBOX+1,SIZE(PFIELD,3),&
           & SIZE(PFIELD,4),SIZE(PFIELD,5),SIZE(PFIELD,6)))
      GALLOC = .TRUE.
    ELSE
      ALLOCATE(ZFIELDP(0,0,0,0,0,0))
      GALLOC = .TRUE.
    END IF
    !
    CALL GATHER_XYFIELD(PFIELD,ZFIELDP,TZFD%OWNER,TZFD%COMM,&
         & KXOBOX,KXEBOX,KYOBOX,KYEBOX,HBUDGET)
    !
    IF (ISP == TZFD%OWNER)  THEN
      TZFMH%GRID    = KGRID
      TZFMH%COMLEN  = LEN_TRIM(HCOMMENT)
      TZFMH%COMMENT = HCOMMENT
      CALL FM_WRIT_ll(TZFD%FLU,HRECFM,.TRUE.,SIZE(ZFIELDP),ZFIELDP,TZFMH&
           & ,IRESP)
    END IF
    !
    CALL MPI_BCAST(IRESP,1,MPI_INTEGER,TZFD%OWNER-1,TZFD&
         & %COMM,IERR)
  END IF ! multiprocessor execution
ELSE
  IRESP = -61
END IF
!----------------------------------------------------------------
IF (IRESP.NE.0) THEN
  CALL FM_WRIT_ERR("FMWRITBOXX6_ll",HFILEM,HFIPRI,HRECFM,'XY',KGRID,LEN(HCOMMENT),IRESP)
END IF
IF (GALLOC) DEALLOCATE(ZFIELDP)
KRESP = IRESP
END SUBROUTINE FMWRITBOXX6_ll
  
END MODULE MODE_FMWRIT


