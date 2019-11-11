!     ########################
      PROGRAM COMPUTE_VER_GRID
!     ########################
!
!!****  *COMPUTE_VER_GRID* - compute the vertigal grid from data in namelist 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/04/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef NAGf95
USE F90_UNIX  ! for GETARG
#endif
!
USE MODE_POS
!
IMPLICIT NONE
!
!*       0.1   Declaration of local variables
!              ------------------------------
!
! Variables for input line 
CHARACTER(LEN=100) :: yexe
integer :: ilenexe
#ifndef NAGf95
INTEGER :: iargc
! CRAY specific
INTEGER :: arglen
!!!!!!!!!!!!!!!!!
#endif
INTEGER :: inarg
!
INTEGER, PARAMETER :: JPVEXT = 1      ! Vertical External points number
!
CHARACTER(LEN=28)              :: YNAM1  ! name of the namelist file
INTEGER                        :: INAM1
CHARACTER(LEN=18)              :: YLUOUT0    ! Name of output_listing file
INTEGER                        :: ILUOUT0
INTEGER                        :: IRESP
LOGICAL                        :: GFOUND
!
INTEGER :: JK         ! vertical loop control
INTEGER :: IKB        ! first inner vertical point index
INTEGER :: IKE        ! last inner vertical point index
INTEGER :: IKU        ! upper vertical point index
INTEGER :: NKMAX      ! Dimensions in z direction
! namelist NAM_VER_GRID in PRE_REAL
CHARACTER(LEN=6) :: YZGRID_TYPE ! type of input vertical grid
REAL :: ZDZGRD        ! vertical mesh length near the ground
REAL :: ZDZTOP        ! vertical mesh length near the top 
REAL :: ZSTRGRD       ! stretching value near the ground
REAL :: ZSTRTOP       ! stretching value near the top of the model
REAL :: ZZMAX_STRGRD  ! maximum height under which the stretching is equal to
                      ! ZSTRGRD
! namelist NAM_GRIDn_PRE in PRE_IDEAL
CHARACTER(LEN=6) :: CZGRID_TYPE ! type of input vertical grid
REAL :: XDZGRD        ! vertical mesh length near the ground
REAL :: XDZTOP        ! vertical mesh length near the top 
REAL :: XSTRGRD       ! stretching value near the ground
REAL :: XSTRTOP       ! stretching value near the top of the model
REAL :: XZMAX_STRGRD  ! maximum height under which the stretching is equal to
REAL :: XLATOR,XLONOR ! latitude and longitude of the Origine point
REAL :: XLATCEN,XLONCEN ! latitude and longitude of the center of the domain
REAL :: XDELTAX,XDELTAY ! horizontal mesh lengths  
REAL :: XHMAX ! Maximum height for orography
REAL :: NEXPX,NEXPY     ! Exponents for  orography in case of CZS='SINE'
REAL :: XAX, XAY        ! Widths for orography in case CZS='BELL'
INTEGER :: NIZS , NJZS  ! Localization of the center in case CZS ='BELL' 
! namelist NAM_DIMn_PRE in PRE_IDEAL
INTEGER :: NIMAX, NJMAX ! Dimensions in x,y directions
!
REAL :: ZSTRETCH      ! running stretching value
LOGICAL :: LTHINSHELL ! thinshell approximation
REAL, DIMENSION(:), ALLOCATABLE :: ZZHAT ! height level without orography
REAL, DIMENSION(:), ALLOCATABLE :: ZSTRETCHING ! stretching between two 
!                                              ! consecutive vertical levels
!
!*       0.3   Declaration of namelists
!              ------------------------
!
! in PRE_REAL1.nam
NAMELIST/NAM_VER_GRID/ LTHINSHELL,NKMAX, &
                     YZGRID_TYPE,ZDZGRD,ZDZTOP,ZZMAX_STRGRD,ZSTRGRD,ZSTRTOP, &
                     NIMAX,NJMAX, &
                     XLONOR,XLATOR,XLATCEN,XLONCEN,XDELTAX,XDELTAY, &
                     XHMAX,NEXPX,NEXPY,XAX,XAY,NIZS,NJZS
! in PRE_IDEA1.nam
NAMELIST/NAM_DIMN_PRE/ NIMAX,NJMAX,NKMAX
NAMELIST/NAM_GRIDN_PRE/ CZGRID_TYPE,XLONOR,XLATOR,XLONCEN,XLATCEN,XDELTAX, &
                       XDELTAY,XDZGRD,XDZTOP,XZMAX_STRGRD,XSTRGRD,XSTRTOP, &
                       XHMAX,NEXPX,NEXPY,XAX,XAY,NIZS,NJZS
!
NAMELIST/NAM_VER_OUT/ YZGRID_TYPE,NKMAX,ZDZGRD,ZDZTOP, &
                      ZZMAX_STRGRD,ZSTRGRD,ZSTRTOP
!-------------------------------------------------------------------------------
!
!*       1.    SET DEFAULT VALUES
!              ------------------
!
NKMAX=0
ZDZGRD=300. ; XDZGRD=300.
ZDZTOP=300. ; XDZTOP=300.
ZZMAX_STRGRD=0. ; XZMAX_STRGRD=0.
ZSTRGRD=0. ; XSTRGRD=0.
ZSTRTOP=0. ; XSTRTOP=0.
YZGRID_TYPE='FUNCTN' ; CZGRID_TYPE='FUNCTN'
!
YLUOUT0='OUTPUT_VER_GRID'
YNAM1='VER_GRID1.nam'
!
!-------------------------------------------------------------------------------
!
!*       2.    RETRIEVE THE NAME OF THE NAMELIST FILE
!              --------------------------------------
!
inarg = iargc()
#if defined(F90HP)
#define HPINCR 1
#else
#define HPINCR 0
#endif
!
#if defined(FUJI) || defined(NAGf95) || defined(NEC) || defined(HP) || defined(pgf) || defined(G95) || defined(GFORTRAN)
CALL GETARG(1+HPINCR,yexe)
IF (LEN_TRIM(yexe) == 0) THEN
  PRINT *, 'FATAL ERROR : Activer la macro -DF90HP dans le Makefile et recompiler'
  STOP
END IF
#else
CALL PXFGETARG(1,yexe,arglen,iresp)
#endif
YNAM1=TRIM(yexe)
PRINT *,'Input file is ',YNAM1
!
!-------------------------------------------------------------------------------
!
!*       3.    OPENNING OF THE FILES
!              ---------------------
!
!CALL FMATTR(YLUOUT0,YLUOUT0,ILUOUT0,IRESP)
ILUOUT0=20
OPEN(ILUOUT0,FILE=YLUOUT0)
!
!CALL FMATTR(YNAM1,YLUOUT0,INAM1,IRESP)
INAM1=21
OPEN(INAM1,FILE=YNAM1,STATUS='OLD',iostat=iresp)
IF (IRESP==0) THEN
  PRINT *,'Opening namelist file ',YNAM1
ELSE
  STOP 'ERROR in opening namelist file'
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.    READING OF THE DATA
!*       1.2   Vertical grid value
!              -------------------
!
CALL POSNAM(INAM1,'NAM_VER_GRID',GFOUND)
IF (GFOUND) THEN
  READ(INAM1,NAM_VER_GRID) 
  PRINT *, '  namelist NAM_VER_GRID read'
ENDIF
!
IF(NKMAX==0) THEN
  CALL POSNAM(INAM1,'NAM_GRIDN_PRE',GFOUND)
  IF (GFOUND) THEN 
    READ(INAM1,NAM_GRIDN_PRE) 
    PRINT *, '  namelist NAM_GRIDN_PRE read'
  ENDIF
  CALL POSNAM(INAM1,'NAM_DIMN_PRE',GFOUND)
  IF (GFOUND) THEN 
    READ(INAM1,NAM_DIMN_PRE) 
    PRINT *, '  namelist NAM_DIMN_PRE read'
  ENDIF
  IRESP=-1
ENDIF
!
IF (NKMAX==0) THEN
  CLOSE(INAM1)
  STOP 'Bad initialization of vertical parameters'
ENDIF
!
IF (IRESP==-1) THEN   ! PRE_IDEA1.nam case
  YZGRID_TYPE=CZGRID_TYPE 
  ZDZGRD=XDZGRD
  ZDZTOP=XDZTOP
  ZZMAX_STRGRD=XZMAX_STRGRD
  ZSTRGRD=XSTRGRD
  ZSTRTOP=XSTRTOP
ENDIF
!-------------------------------------------------------------------------------
!
!*       5.    COMPUTATION OF VERTICAL STRETCHING :
!              ----------------------------------
!
IKB=JPVEXT+1
IKE=NKMAX+JPVEXT
IKU=NKMAX+2*JPVEXT
!
IF (.NOT. ALLOCATED(ZZHAT)) ALLOCATE(ZZHAT(IKU))
!
IF (YZGRID_TYPE=='FUNCTN') THEN
!
  IF (ABS(ZDZTOP-ZDZGRD) < 1.E-10) THEN
    ZZHAT(:) = (/ (FLOAT(JK-IKB)*ZDZGRD, JK=1,IKU) /)
!
  ELSE
    IF (ZDZGRD>ZDZTOP) THEN
      WRITE(ILUOUT0,*) 'ZDZGRD MUST BE SMALLER THAN OR EQUAL TO ZDZTOP'
      WRITE(ILUOUT0,*) 'CHANGE THESE PARAMETERS AND TRY AGAIN'
      WRITE(ILUOUT0,*) 'ZDZGRD =', ZDZGRD,'  ZDZTOP =', ZDZTOP
      STOP
    END IF 
!
    ZZHAT(IKB-1)=-ZDZGRD
    ZZHAT(IKB)= 0.
    ZZHAT(IKB+1)=ZDZGRD
    DO JK=IKB+2,IKU
      IF ( ZZHAT(JK-1) < ZZMAX_STRGRD - 1.E-10 ) THEN
        ZSTRETCH=ZSTRGRD/100.
      ELSE
        ZSTRETCH=ZSTRTOP/100.
      END IF
!
      ZZHAT(JK)=ZZHAT(JK-1)+(ZZHAT(JK-1)-ZZHAT(JK-2))*(1.+ZSTRETCH)
!
      IF ( ZZHAT(JK)-ZZHAT(JK-1) > ZDZTOP ) THEN
        ZZHAT(JK)=ZZHAT(JK-1)+ZDZTOP
      END IF
    END DO
!
  END IF
!
END IF
!-------------------------------------------------------------------------------
!
!*       6.    MANUALLY SPECIFIED LEVELS :
!              -------------------------
!
IF (YZGRID_TYPE=='MANUAL') THEN
!
  CALL POSKEY(INAM1,ILUOUT0,'ZHAT')
  READ(INAM1,*) (ZZHAT(JK), JK=JPVEXT+1,NKMAX+JPVEXT+1)
!
  DO JK=JPVEXT,1,-1
    ZZHAT(JK)=ZZHAT(JK+1) - (ZZHAT(JPVEXT+2)-ZZHAT(JPVEXT+1))
  END DO
  DO JK=NKMAX+JPVEXT+2,IKU
    ZZHAT(JK)=ZZHAT(JK-1) + (ZZHAT(NKMAX+JPVEXT+1)-ZZHAT(NKMAX+JPVEXT))
  END DO
!
END IF
!
!-------------------------------------------------------------------------------
!
!*       7.    TEST ON STRETCHING :
!              ------------------
!
WRITE(ILUOUT0,nml=NAM_VER_OUT) 
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,1) 1,ZZHAT(1)
WRITE(ILUOUT0,1) 2,ZZHAT(2)
ALLOCATE(ZSTRETCHING(IKU))
DO JK=3,IKU
  ZSTRETCHING(JK)=(ZZHAT(JK)-ZZHAT(JK-1))/(ZZHAT(JK-1)-ZZHAT(JK-2))-1.
  IF ( ABS(ZSTRETCHING(JK) ) > 0.20 + 1.E-10 ) THEN
     WRITE(ILUOUT0,4) JK,ZZHAT(JK),100.*ZSTRETCHING(JK)
  ELSE IF ( ABS(ZSTRETCHING(JK) ) > 0.07 ) THEN
     WRITE(ILUOUT0,3) JK,ZZHAT(JK),100.*ZSTRETCHING(JK)
  ELSE
     WRITE(ILUOUT0,2) JK,ZZHAT(JK),100.*ZSTRETCHING(JK)
  ENDIF
ENDDO
IF ( ANY(ABS(ZSTRETCHING(3:) ) > 0.20 + 1.E-10 ) ) THEN
  WRITE(ILUOUT0,*)
  WRITE(ILUOUT0,*) '   +-------------------------------------+'
  WRITE(ILUOUT0,*) '   | STRETCHING TOO HIGH (MORE THAN 20%) |'
  WRITE(ILUOUT0,*) '   +-------------------------------------+'
  WRITE(ILUOUT0,*)
  STOP
END IF
WRITE(ILUOUT0,*)
!
1 FORMAT('ZHAT(',I3,')=',F18.12)
2 FORMAT('ZHAT(',I3,')=',F18.12,'  (+',F6.2,' %)')
3 FORMAT('ZHAT(',I3,')=',F18.12,'  (+',F6.2,' %) WARNING: high stretching')
4 FORMAT('ZHAT(',I3,')=',F18.12,'  (+',F6.2,' %) ERROR  : stretching too high')
!
DEALLOCATE(ZSTRETCHING)
!-------------------------------------------------------------------------------
!
PRINT *, 'VERGRID completed'
PRINT *, '=> output grid and stretching in file ', YLUOUT0
!
!-------------------------------------------------------------------------------
!
!*       8.    CLOSING OF THE FILES
!              --------------------
!
CLOSE(INAM1)
!CALL FMFREE(YNAM1,YLUOUT0,IRESP)
CLOSE(ILUOUT0)
!CALL FMFREE(YLUOUT0,YLUOUT0,IRESP)
!
!-------------------------------------------------------------------------------
!
END PROGRAM COMPUTE_VER_GRID
