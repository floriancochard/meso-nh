!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_READ_VER_GRID
!     #########################
INTERFACE
      SUBROUTINE READ_VER_GRID(TPPRE_REAL1,PZHAT,OSLEVE,PLEN1,PLEN2)
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
TYPE(TFILEDATA),POINTER,      INTENT(IN) :: TPPRE_REAL1! namelist file
REAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: PZHAT      ! vertival grid of input fmfile
LOGICAL,            OPTIONAL, INTENT(IN) :: OSLEVE     ! flag for SLEVE coordinate
REAL,               OPTIONAL, INTENT(IN) :: PLEN1      ! Decay scale for smooth topography
REAL,               OPTIONAL, INTENT(IN) :: PLEN2      ! Decay scale for small-scale topography deviation
!
END SUBROUTINE READ_VER_GRID
END INTERFACE
END MODULE MODI_READ_VER_GRID
!
!     ##############################################################
      SUBROUTINE READ_VER_GRID(TPPRE_REAL1,PZHAT,OSLEVE,PLEN1,PLEN2)
!     ##############################################################
!
!!****  *READ_VER_GRID* - reads namelist data in file PRE_REAL1 for the 
!!                        initialization of the vertical grid for real cases.
!! 
!!    PURPOSE
!!    -------
!!    This routine initializes the vertical grid definition variables 
!!    aand configuration variables defined in the namelist file
!!    given by the MESO-NH user, and stores these variables in the adequates
!!    modules.                                                        _
!!    In particular, this routine initializes the vertical coordinate z array.
!!
!!**  METHOD
!!    ------
!!
!!    The routine reads NKMAX, the number of mass levels.
!!                              ^
!!    The routine allocates the z array at size NKMAX+2*JPVEXT
!!                                          ^
!!    The routine initializes the values of z as follows:
!!
!!    IF YZGRID_TYPE =='SAMEGR'
!!
!!      the vertical grid is the same as in the input mesonh file
!! 
!!    IF YZGRID_TYPE =='FUNCTN'
!!
!!      if a stretching value is given in namelist, 
!!      then     the stretching is constant from the low grid mesh value to 
!!               the high one
!!      else     The grid is stretched along the z direction, the mesh varies 
!!               from XDZGRD near the ground to XDZTOP near the top and the 
!!               weigthing function is a TANH function characterized by its 
!!               center and width above and under this center.
!!
!!    IF YZGRID_TYPE =='MANUAL'
!!
!!      The NKMAX+1 levels are read in * format after the namelists, from 
!!      ground level to rigid top level
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_GRID1
!!         XZHAT
!!      Module MODD_DIM1
!!         NKMAX
!!      Module MODD_PARAMETERS
!!         JPVEXT
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
!!      Original    16/12/94
!!                  Sept 21, 1995 (V.Masson J.Stein) change the namelist  
!!                  Aug. 10, 1996 (V. Masson) add LTHINSHELL in namelist
!!                  Oct, 25, 1996 (V.Masson) deallocations
!!                  Oct. 10, 2001 (I.Mallet) allow namelists in different orders
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF           ! declaration modules
USE MODD_DIM_n, NKMAX_n=>NKMAX
USE MODD_GRID_n, LSLEVE_n=>LSLEVE, XLEN1_n=>XLEN1, XLEN2_n=>XLEN2
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_LUNIT
USE MODD_PARAMETERS
!
USE MODE_FM
USE MODE_MSG
USE MODE_POS
!
USE MODI_DEFAULT_SLEVE
!
USE MODN_BLANK
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
TYPE(TFILEDATA),POINTER,      INTENT(IN) :: TPPRE_REAL1! namelist file
LOGICAL,            OPTIONAL, INTENT(IN) :: OSLEVE     ! flag for SLEVE coordinate
REAL,               OPTIONAL, INTENT(IN) :: PLEN1      ! Decay scale for smooth topography
REAL,               OPTIONAL, INTENT(IN) :: PLEN2      ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:), OPTIONAL, INTENT(IN) :: PZHAT      ! vertival grid of input fmfile
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER :: IPRE_REAL1 ! logical unit for namelist file HPRE_REAL1
INTEGER :: ILUOUT0    ! logical unit for output listing file
INTEGER :: IRESP      ! return-code if problems eraised
LOGICAL :: GFOUND     !    "        when reading namelist
INTEGER :: JK         ! vertical loop control
INTEGER :: IKB        ! first inner vertical point index
INTEGER :: IKU        ! upper vertical point index
CHARACTER(LEN=6) :: YZGRID_TYPE ! type of input vertical grid
REAL :: ZDZGRD        ! vertical mesh length near the ground
REAL :: ZDZTOP        ! vertical mesh length near the top 
REAL :: ZSTRGRD       ! stretching value near the ground
REAL :: ZSTRTOP       ! stretching value near the top of the model
REAL :: ZZMAX_STRGRD  ! maximum height under which the stretching is equal to
                      ! ZSTRGRD
REAL :: ZSTRETCH      ! running stretching value
LOGICAL :: GTHINSHELL ! thinshell approximation
INTEGER :: NKMAX ! re-declaration because of namelist
LOGICAL :: LSLEVE ! re-declaration because of namelist
REAL    :: XLEN2 ! re-declaration because of namelist
REAL    :: XLEN1 ! re-declaration because of namelist
REAL, DIMENSION(:), ALLOCATABLE :: ZSTRETCHING ! stretching between two 
!                                              ! consecutive vertical levels
!
!*       0.3   Declaration of namelists
!              ------------------------
!
NAMELIST/NAM_VER_GRID/ LTHINSHELL,NKMAX,YZGRID_TYPE,ZDZGRD,ZDZTOP,ZZMAX_STRGRD,ZSTRGRD,ZSTRTOP,&
                       LSLEVE,XLEN1,XLEN2
!-------------------------------------------------------------------------------
!
!*       1.    READING OF THE NAMELIST IN FILE HPRE_REAL1 :  
!              ------------------------------------------
ILUOUT0 = TLUOUT0%NLU
IPRE_REAL1 = TPPRE_REAL1%NLU
!
!*       1.1   Vertical grid default value
!              ---------------------------
!
IF (PRESENT(PZHAT)) THEN
  YZGRID_TYPE='SAMEGR'
ELSE
  YZGRID_TYPE='FUNCTN'
END IF
!
GTHINSHELL=LTHINSHELL
NKMAX_n=0
ZDZGRD=300.
ZDZTOP=300.
ZZMAX_STRGRD=0.
ZSTRGRD=0.
ZSTRTOP=0.
CALL DEFAULT_SLEVE(LSLEVE,XLEN1,XLEN2)
!
!*       1.2   Vertical grid value
!              -------------------
!
NKMAX=NKMAX_n
CALL POSNAM(IPRE_REAL1,'NAM_VER_GRID',GFOUND,ILUOUT0)
IF (GFOUND) READ(UNIT=IPRE_REAL1,NML=NAM_VER_GRID) 
NKMAX_n  = NKMAX
LSLEVE_n = LSLEVE
XLEN1_n  = XLEN1
XLEN2_n  = XLEN2
!
IF (CPROGRAM=='REAL  ') THEN
  IF (ASSOCIATED (XZHAT) ) DEALLOCATE(XZHAT)
  CALL POSNAM(IPRE_REAL1,'NAM_BLANK',GFOUND,ILUOUT0)
  IF (GFOUND) READ(UNIT=IPRE_REAL1,NML=NAM_BLANK)
END IF 
!
IKB=JPVEXT+1
IKU=NKMAX_n+2*JPVEXT
!
SELECT CASE(YZGRID_TYPE)
!-------------------------------------------------------------------------------
!
!*       2.    KEEPING OF THE SAME GRID :
!              ------------------------
!
CASE('SAMEGR')
  IF (PRESENT(PZHAT) .AND. PRESENT(OSLEVE) .AND. PRESENT(PLEN1) .AND.  PRESENT(PLEN2)) THEN
    IF (NKMAX_n==0)  NKMAX_n=SIZE(PZHAT)-2*JPVEXT
    ALLOCATE(XZHAT(NKMAX_n+2*JPVEXT))

    IF ( (NKMAX_n+2*JPVEXT) > SIZE(PZHAT) ) THEN
      WRITE(ILUOUT0,*) 'ERROR IN READ_VER_GRID :'
      WRITE(ILUOUT0,*) '  YOU WANT TO KEEP THE SAME VERTICAL GRID, BUT YOU ASK'
      WRITE(ILUOUT0,*) '  FOR MORE LEVELS THAN IN INPUT FM FILE.'
      WRITE(ILUOUT0,*) '  FM file: ',SIZE(PZHAT)-2*JPVEXT
      WRITE(ILUOUT0,*) '  NKMAX  : ',NKMAX_n
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_VER_GRID','')
    END IF

    XZHAT(:)=PZHAT(1:NKMAX_n+2*JPVEXT)
    LTHINSHELL = GTHINSHELL
    LSLEVE_n     = OSLEVE
    XLEN1_n      = PLEN1
    XLEN2_n      = PLEN2

    IF ( (NKMAX_n+2*JPVEXT) == SIZE(PZHAT) ) THEN
      WRITE(ILUOUT0,*) 'same vertical grid kept.'
    ELSE
      WRITE(ILUOUT0,*) NKMAX_n,' first levels in vertical grid kept.'
    END IF
    WRITE(ILUOUT0,*) 'LTHINSHELL set to : ',LTHINSHELL
    WRITE(ILUOUT0,*) 'LSLEVE     set to : ',LSLEVE_n
    IF (LSLEVE_n) THEN
      WRITE(ILUOUT0,*) 'XLEN1      set to : ',XLEN1_n
      WRITE(ILUOUT0,*) 'XLEN2      set to : ',XLEN2_n
    END IF

  ELSE
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_VER_GRID','the vertical grid can be kept only if the input file is a MESONH file')
  END IF
!
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTATION OF VERTICAL STRETCHING :
!              ----------------------------------
!
CASE('FUNCTN')
!
  IF (.NOT. ASSOCIATED(XZHAT)) ALLOCATE(XZHAT(IKU))
!
  IF (ABS(ZDZTOP-ZDZGRD) < 1.E-10) THEN
    XZHAT(:) = (/ (FLOAT(JK-IKB)*ZDZGRD, JK=1,IKU) /)
!
  ELSE
    IF (ZDZGRD>ZDZTOP) THEN
      WRITE(ILUOUT0,*) 'ZDZGRD MUST BE SMALLER THAN OR EQUAL TO ZDZTOP'
      WRITE(ILUOUT0,*) 'CHANGE THESE PARAMETERS AND TRY AGAIN'
      WRITE(ILUOUT0,*) 'ZDZGRD =', ZDZGRD,'  ZDZTOP =', ZDZTOP
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_VER_GRID','')
    END IF 
!
    XZHAT(IKB-1)=-ZDZGRD
    XZHAT(IKB)= 0.
    XZHAT(IKB+1)=ZDZGRD
    DO JK=IKB+2,IKU
      IF ( XZHAT(JK-1) < ZZMAX_STRGRD - 1.E-10 ) THEN
        ZSTRETCH=ZSTRGRD/100.
      ELSE
        ZSTRETCH=ZSTRTOP/100.
      END IF
!
      XZHAT(JK)=XZHAT(JK-1)+(XZHAT(JK-1)-XZHAT(JK-2))*(1.+ZSTRETCH)
!
      IF ( XZHAT(JK)-XZHAT(JK-1) > ZDZTOP ) THEN
        XZHAT(JK)=XZHAT(JK-1)+ZDZTOP
      END IF
    END DO
!
  END IF
!
!-------------------------------------------------------------------------------
!
!*       4.    MANUALLY SPECIFIED LEVELS :
!              -------------------------
!
CASE('MANUAL')
!
  IF (.NOT. ASSOCIATED(XZHAT)) ALLOCATE(XZHAT(IKU))
!
  WRITE(ILUOUT0,FMT=*) 'YZGRID_TYPE="MANUAL", ATTEMPT TO READ VECTOR XZHAT(2,NKU)'
  CALL POSKEY(IPRE_REAL1,ILUOUT0,'ZHAT')
  READ(IPRE_REAL1,*) (XZHAT(JK), JK=JPVEXT+1,NKMAX_n+JPVEXT+1)
!
  DO JK=JPVEXT,1,-1
    XZHAT(JK)=XZHAT(JK+1) - (XZHAT(JPVEXT+2)-XZHAT(JPVEXT+1))
  END DO
  DO JK=NKMAX_n+JPVEXT+2,IKU
    XZHAT(JK)=XZHAT(JK-1) + (XZHAT(NKMAX_n+JPVEXT+1)-XZHAT(NKMAX_n+JPVEXT))
  END DO
!
!-------------------------------------------------------------------------------
!
CASE DEFAULT
  WRITE(ILUOUT0,FMT=*) 'UNKNOWN TYPE OF VERTICAL GRID GENERATION'
  WRITE(ILUOUT0,FMT=*) ' YZGRID_TYPE =',YZGRID_TYPE
  WRITE(ILUOUT0,FMT=*) '-> JOB ABORTED'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','READ_VER_GRID','')
END SELECT
!
!Set model top
XZTOP = XZHAT(IKU)
!
!-------------------------------------------------------------------------------
!
!*       5.    TEST ON STRETCHING :
!              ------------------
!
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,1) 1,XZHAT(1)
WRITE(ILUOUT0,1) 2,XZHAT(2)
ALLOCATE(ZSTRETCHING(IKU))
DO JK=3,IKU
  ZSTRETCHING(JK)=(XZHAT(JK)-XZHAT(JK-1))/(XZHAT(JK-1)-XZHAT(JK-2))-1.
  IF ( ABS(ZSTRETCHING(JK) ) > 0.20 + 1.E-10 ) THEN
     WRITE(ILUOUT0,4) JK,XZHAT(JK),100.*ZSTRETCHING(JK)
  ELSE IF ( ABS(ZSTRETCHING(JK) ) > 0.07 ) THEN
     WRITE(ILUOUT0,3) JK,XZHAT(JK),100.*ZSTRETCHING(JK)
  ELSE
     WRITE(ILUOUT0,2) JK,XZHAT(JK),100.*ZSTRETCHING(JK)
  ENDIF
ENDDO

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
WRITE(ILUOUT0,*) 'Routine READ_VER_GRID completed'
!
END SUBROUTINE READ_VER_GRID
