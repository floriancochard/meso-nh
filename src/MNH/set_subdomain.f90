!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_SET_SUBDOMAIN
!     #########################
INTERFACE
      SUBROUTINE SET_SUBDOMAIN(TPNMLFILE,TPATMFILE,KXOR_DAD,KYOR_DAD, &
                               KXOR,KYOR,KDXRATIO,KDYRATIO            )
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
TYPE(TFILEDATA),   INTENT(IN) :: TPNMLFILE ! namelist file
TYPE(TFILEDATA),   INTENT(IN) :: TPATMFILE ! atmospheric MNH file
INTEGER,           INTENT(OUT):: KXOR_DAD  ! XOR compared to Dad file, if any
INTEGER,           INTENT(OUT):: KYOR_DAD  ! YOR compared to Dad file, if any
INTEGER,           INTENT(OUT):: KXOR      ! XOR given or computed
INTEGER,           INTENT(OUT):: KYOR      ! YOR given or computed
INTEGER,           INTENT(OUT):: KDXRATIO  ! DXRATIO compared to Dad file, if any
INTEGER,           INTENT(OUT):: KDYRATIO  ! DYRATIO compared to Dad file, if any
!
END SUBROUTINE SET_SUBDOMAIN
END INTERFACE
END MODULE MODI_SET_SUBDOMAIN
!
!     #################################################################
      SUBROUTINE SET_SUBDOMAIN(TPNMLFILE,TPATMFILE,KXOR_DAD,KYOR_DAD, &
                               KXOR,KYOR,KDXRATIO,KDYRATIO            )
!     #################################################################
!
!!****  *SET_SUBDOMAIN* - computes the horizontal MESO-NH domain
!! 
!!    PURPOSE
!!    -------
!!    This routine computes the boundaries of the MESO-NH domain inside the
!!    Aladin domain. This subroutine initializes also XXHAT and XYHAT.
!!
!!**  METHOD
!!    ------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         TLUOUT0 : output-listing file
!!      Module MODD_PARAMETERS
!!         JPHEXT
!!      Module MODD_GRID1
!!         XLONOR
!!         XLATOR
!!         XXHAT
!!         XYHAT
!!      Module MODD_DIM1
!!         NIMAX
!!         NJMAX
!!      Module MODD_GRID   : contains grid definition
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
!!      Original     16/12/94
!!      Modification 30/01/96 Flux point as center of the domain (V. Masson)
!!      Modification 02/10/96 coordinates are copied from pgd file
!!                   25/10/96 Add deallocations (V. Masson)
!!                   15/03/99 Possibility to use IXOR and IYOR (V. Masson)
!!                   15/10/01 Allow namelists in different orders (I. Mallet)
!!                   30/10/03 No more used in masdev4_6 version (V. Masson)
!!                   30/11/04 Used again in masdev4_7 version (V. Masson)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_GRIDPROJ       ! executive module
USE MODE_POS
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
!
USE MODD_CONF           ! declaration modules
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_LUNIT
USE MODD_GRID
USE MODD_GRID_n
USE MODD_DIM_n, ONLY: NIMAX_n=>NIMAX,NJMAX_n=>NJMAX
USE MODD_PGDGRID
USE MODD_PGDDIM
USE MODD_PARAMETERS
!
USE MODE_FMREAD
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
TYPE(TFILEDATA),   INTENT(IN) :: TPNMLFILE ! namelist file
TYPE(TFILEDATA),   INTENT(IN) :: TPATMFILE ! atmospheric MNH file
INTEGER,           INTENT(OUT):: KXOR_DAD  ! XOR compared to Dad file, if any
INTEGER,           INTENT(OUT):: KYOR_DAD  ! YOR compared to Dad file, if any
INTEGER,           INTENT(OUT):: KXOR      ! XOR given or computed
INTEGER,           INTENT(OUT):: KYOR      ! YOR given or computed
INTEGER,           INTENT(OUT):: KDXRATIO  ! DXRATIO compared to Dad file, if any
INTEGER,           INTENT(OUT):: KDYRATIO  ! DYRATIO compared to Dad file, if any
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
CHARACTER(LEN=28) :: YDADFILE ! name of atmospheric file father (if any)
INTEGER :: ILUNAM     ! logical unit for namelist file HNAMELIST
INTEGER :: ILUOUT0    ! logical unit for listing file
INTEGER :: IRESP      ! return-code if problems araised
INTEGER :: IIU,IJU    ! number of points in x,y directions of the large grid
INTEGER :: IIINF      ! left point index
INTEGER :: IISUP      ! right point index
INTEGER :: IJINF      ! bottom point index
INTEGER :: IJSUP      ! top point index
                      ! Declare namelist variables
                               ! horizontal position (i,j) of the ORigin
INTEGER :: NXOR,NYOR           ! of the model 2 domain, relative to model 1
INTEGER :: NIMAX,NJMAX         ! Dimensions respectively in x,  y directions
REAL, DIMENSION(:),   ALLOCATABLE :: ZXHAT      ! = XXHAT(:)
REAL, DIMENSION(:),   ALLOCATABLE :: ZYHAT      ! = XYHAT(:)
!
LOGICAL :: GFOUND     ! return code when searching namelist
!
!*       0.3   Declaration of namelists
!              ------------------------
NAMELIST/NAM_MESONH_DOM/ NIMAX,NJMAX,NXOR,NYOR
!                                                    ! definition of MESO-NH
!                                                    ! horizontal domain
!-------------------------------------------------------------------------------
ILUOUT0 = TLUOUT0%NLU
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF THE LATITUDE AND LONGITUDE OF LARGE DOMAIN POINTS
!              ----------------------------------------------------------------
!
ALLOCATE(ZXHAT(SIZE(XXHAT)))
ALLOCATE(ZYHAT(SIZE(XYHAT)))
ZXHAT(:) = XXHAT(:)
ZYHAT(:) = XYHAT(:)
!
IIU=SIZE(ZXHAT)
IJU=SIZE(ZYHAT)
!
KXOR_DAD=1
KYOR_DAD=1
KDXRATIO=1
KDYRATIO=1
!-------------------------------------------------------------------------------
!
!*       2.    DOMAIN DEFINITION
!              -----------------
!
!*       2.1   Default values
!              --------------
!
NXOR=NUNDEF
NYOR=NUNDEF
!
!
!*       2.2   Reading of the namelist file (or given as dummy argument)
!              ---------------------------- 
!
ILUNAM = TPNMLFILE%NLU
CALL POSNAM(ILUNAM,'NAM_MESONH_DOM',GFOUND,ILUOUT0)  
IF (GFOUND) THEN
  NIMAX=NIMAX_n ; NJMAX=NJMAX_n
  READ(UNIT=ILUNAM,NML=NAM_MESONH_DOM)
  NIMAX_n=NIMAX ; NJMAX_n=NJMAX
END IF
!
!*       2.3   Default values if none was given in namelist: domain is centered
!              ---------------------------- ---------------
!
IF ( NXOR==NUNDEF .OR. NYOR==NUNDEF ) THEN
  IF ( NXOR/=NUNDEF ) THEN
    WRITE(ILUOUT0,*) 'bottom-left domain X index is given, but not its Y index'
    WRITE(ILUOUT0,*) 'NXOR = ',NXOR
    WRITE(ILUOUT0,*) 'Please also define NYOR or choose an other domain initialization mode'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_SUBDOMAIN','bottom-left domain X index is given, but not its Y index')
  END IF
  IF ( NYOR/=NUNDEF ) THEN
    WRITE(ILUOUT0,*) 'bottom-left domain Y index is given, but not its X index'
    WRITE(ILUOUT0,*) 'NYOR = ',NYOR
    WRITE(ILUOUT0,*) 'Please also define NXOR or choose an other domain initialization mode'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_SUBDOMAIN','bottom-left domain Y index is given, but not its X index')
  END IF
  NXOR= (IIU - NIMAX)/2
  NYOR= (IJU - NJMAX)/2
END IF
!
!-------------------------------------------------------------------------------
!
WRITE(ILUOUT0,*) 'given or computed NIMAX = ',NIMAX
WRITE(ILUOUT0,*) 'given or computed NJMAX = ',NJMAX
WRITE(ILUOUT0,*) 'given or computed NXOR  = ',NXOR
WRITE(ILUOUT0,*) 'given or computed NYOR  = ',NYOR
!
!-------------------------------------------------------------------------------
!
!*       4.    COHERENCE WITH GRID MESH INDICES
!              --------------------------------
!
!*       4.1   TEST if FATHER of atmospheric MNH file exists:
!              ---------------------------------------------
!
CALL IO_READ_FIELD(TPATMFILE,'DAD_NAME',YDADFILE,IRESP)
IF ( IRESP /= 0  ) YDADFILE='                          '
!
IF (LEN_TRIM(YDADFILE)/=0) THEN
  CALL IO_READ_FIELD(TPATMFILE,'DXRATIO',KDXRATIO,IRESP)
  IF ( IRESP /= 0 .OR. KDXRATIO == 0 ) THEN
    KDXRATIO=1
  END IF
  !
  CALL IO_READ_FIELD(TPATMFILE,'DYRATIO',KDYRATIO,IRESP)
  IF ( IRESP /= 0 .OR. KDYRATIO == 0 ) THEN
    KDYRATIO=1
  END IF
  !
  CALL IO_READ_FIELD(TPATMFILE,'XOR',KXOR,IRESP)
  IF ( IRESP /= 0 ) KXOR_DAD=1
  !
  CALL IO_READ_FIELD(TPATMFILE,'YOR',KYOR,IRESP)
  IF ( IRESP /= 0 ) KYOR_DAD=1
END IF
!
!*       4.2   CONTROL of new domain position relatively to father mesh size
!              -------------------------------------------------------------
!
IF ( MOD(NIMAX,KDXRATIO) /= 0 ) THEN
  WRITE(ILUOUT0,*) '**************************************'
  WRITE(ILUOUT0,*) 'In order to be included in father grid'
  WRITE(ILUOUT0,*) 'NIMAX must be divisible by: ', KDXRATIO
  WRITE(ILUOUT0,*) 'Your NIMAX is :', NIMAX
  WRITE(ILUOUT0,*) '**************************************'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_SUBDOMAIN','')
END IF
!
IF ( MOD(NJMAX,KDYRATIO) /= 0 ) THEN
  WRITE(ILUOUT0,*) '**************************************'
  WRITE(ILUOUT0,*) 'In order to be included in father grid'
  WRITE(ILUOUT0,*) 'NJMAX must be divisible by: ', KDYRATIO
  WRITE(ILUOUT0,*) 'Your NJMAX is :', NJMAX
  WRITE(ILUOUT0,*) '**************************************'
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_SUBDOMAIN','')
END IF
!
IF ( MOD(NXOR-1,KDXRATIO) /= 0 ) THEN
  NXOR = KDXRATIO * ((NXOR-1)/KDXRATIO) + 1
  WRITE(ILUOUT0,*) '**************************************'
  WRITE(ILUOUT0,*) 'In order to be included in father grid'
  WRITE(ILUOUT0,*) 'NXOR is set to : ', NXOR
  WRITE(ILUOUT0,*) '**************************************'
END IF
KXOR=NXOR
!* new left coordinate compared to father grid, if any
KXOR_DAD = KXOR_DAD + (NXOR-1)/KDXRATIO
!
IF ( MOD(NYOR-1,KDYRATIO) /= 0 ) THEN
  NYOR = KDYRATIO * ((NYOR-1)/KDYRATIO) + 1
  WRITE(ILUOUT0,*) '**************************************'
  WRITE(ILUOUT0,*) 'In order to be included in father grid'
  WRITE(ILUOUT0,*) 'NYOR is set to : ', NYOR
  WRITE(ILUOUT0,*) '**************************************'
END IF
KYOR=NYOR
!* new bottom coordinate compared to father grid, if any
KYOR_DAD = KYOR_DAD + (NYOR-1)/KDYRATIO
!
!-------------------------------------------------------------------------------
!
!*       5.    COMPUTATION OF THE BOUNDARIES OF THE MESO-NH DOMAIN 
!              ---------------------------------------------------
!
IIINF=NXOR
IJINF=NYOR
IISUP=NXOR  + NIMAX    + 2*JPHEXT - 1
IJSUP=NYOR  + NJMAX    + 2*JPHEXT - 1
!
!-------------------------------------------------------------------------------
!
!*       6.    TEST ON THE BOUNDARIES
!              ----------------------
!
WRITE(ILUOUT0,*) 'IINF=', IIINF, ' JINF=', IJINF
WRITE(ILUOUT0,*) 'ISUP=', IISUP, ' JSUP=', IJSUP
!
IF (   (IIINF<1)              &
  .OR. (IJINF<1)              &
  .OR. (IISUP>IIU)            &
  .OR. (IJSUP>IJU) )          &
THEN
  WRITE(ILUOUT0,*)
  WRITE(ILUOUT0,*) '  +----------------------------------------------------------+'
  WRITE(ILUOUT0,*) '  | INPUT PGD DOMAIN TOO SMALL OR OUTPUT DOMAIN NOT WELL CENTERED |'
  WRITE(ILUOUT0,*) '  +----------------------------------------------------------+'
  WRITE(ILUOUT0,*)
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_SUBDOMAIN','input PGD domain too small or output domain not well centered')
ENDIF
!
!-------------------------------------------------------------------------------
!                        
!*       8.    COMPUTATION OF xhat AND yhat
!              ----------------------------
!
DEALLOCATE(XXHAT)
DEALLOCATE(XYHAT)
ALLOCATE(XXHAT(IISUP-IIINF+1))
XXHAT(:)=ZXHAT(IIINF:IISUP)
ALLOCATE(XYHAT(IJSUP-IJINF+1))
XYHAT(:)=ZYHAT(IJINF:IJSUP)
DEALLOCATE(ZXHAT)
DEALLOCATE(ZYHAT)
!-------------------------------------------------------------------------------
!                        
WRITE(ILUOUT0,*) 'Routine SET_SUBDOMAIN completed'
!
END SUBROUTINE SET_SUBDOMAIN
