!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!###################
MODULE MODI_CHECK_ZS
!###################
INTERFACE
      SUBROUTINE CHECK_ZS(HFMFILE,HDAD_NAME,KIINF,KJINF)
!
CHARACTER(LEN=*),    INTENT(IN)    :: HFMFILE   ! name of the Mesonh input file
CHARACTER(LEN=*),    INTENT(INOUT) :: HDAD_NAME ! true name of the Mesonh input file
INTEGER,             INTENT(IN)    :: KIINF     ! indexes of new truncated
INTEGER,             INTENT(IN)    :: KJINF     ! domain, compared to the old
                                                ! one (input pgd)
!
END SUBROUTINE CHECK_ZS
END INTERFACE
END MODULE MODI_CHECK_ZS
!     ##################################################
      SUBROUTINE CHECK_ZS(HFMFILE,HDAD_NAME,KIINF,KJINF)
!     ##################################################
!
!!****  *CHECK_ZS* - checks coherence between large scale and fine orographies
!!                   for nesting purposes
!! 
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!  1 Resolution ratios during previous spawning are read in FM file.
!!  2 The 2 orographies are averaged on the grid with coarser resolution.
!!  3 If the 2 orographies on coarse grids are identical, then DAD_NAME
!!    is kept; if not, it is initialized to ' ', and nesting wont be
!!    allowed between the output file and its father.
!!
!!    EXTERNAL
!!    --------
!!
!!    function FMREAD
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     : contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_GRID1 
!!         XZS
!!      Module MODD_DIM1 
!!         NIMAX,NJMAX
!!      Module MODD_PARAMETERS
!!         JPHEXT
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
!!      Original    24/09/96
!!                  20/05/98  (V. Masson and J. Stein) include the case where
!!                             the domain is reduced 
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_DIM_n
USE MODD_GRID_n
USE MODD_IO_ll,            ONLY: TFILEDATA
USE MODD_LUNIT,            ONLY: TLUOUT0
USE MODD_NESTING
USE MODD_PARAMETERS
!
USE MODE_FMREAD
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_FIND_BYNAME
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
CHARACTER(LEN=*),    INTENT(IN)    :: HFMFILE   ! name of the Mesonh input file
CHARACTER(LEN=*),    INTENT(INOUT) :: HDAD_NAME ! true name of the Mesonh input file
INTEGER,             INTENT(IN)    :: KIINF     ! indexes of new truncated
INTEGER,             INTENT(IN)    :: KJINF     ! domain, compared to the old
                                                ! one (input pgd)
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER             :: IRESP                ! return-code if problems occured
INTEGER             :: ILUOUT0              ! logical unit for file CLUOUT0
!
INTEGER             :: IDXRATIO = 0         ! aspect ratios during previous
INTEGER             :: IDYRATIO = 0         ! spawning (if any)
INTEGER             :: IIMAX, IJMAX         ! dimensions printed in input file
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS    ! zs in input file
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS1   ! averaged orographies
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS2   !
INTEGER                           :: IIMAXC ! inner dimensions of coarse
INTEGER                           :: IJMAXC ! definition arrays
INTEGER                           :: JI,JJ  ! loop counters
TYPE(TFILEDATA),POINTER           :: TZFMFILE
!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
!-------------------------------------------------------------------------------
!
!*            1. Reading of aspect ratios and dimensions
!                ---------------------------------------
!
CALL IO_FILE_FIND_BYNAME(TRIM(HFMFILE),TZFMFILE,IRESP)
!
CALL IO_READ_FIELD(TZFMFILE,'DXRATIO',IDXRATIO,IRESP)
IF ( IRESP /= 0 .OR. IDXRATIO == 0 ) THEN
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'resolution ratio in x direction not present in fmfile; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
  HDAD_NAME='                            '
  RETURN
END IF
!
CALL IO_READ_FIELD(TZFMFILE,'DYRATIO',IDYRATIO,IRESP)
IF ( IRESP /= 0 .OR. IDYRATIO == 0 ) THEN
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'resolution ratio in y direction not present in fmfile; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
  HDAD_NAME='                            '
  RETURN 
END IF
!
CALL IO_READ_FIELD(TZFMFILE,'XOR',NXOR_ALL(1),IRESP)
IF ( IRESP /= 0 ) THEN
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'position XOR not present in fmfile; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
  HDAD_NAME='                            '
  RETURN
END IF
!
CALL IO_READ_FIELD(TZFMFILE,'YOR',NYOR_ALL(1),IRESP)
IF ( IRESP /= 0 ) THEN
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'resolution YOR not present in fmfile; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
  HDAD_NAME='                            '
  RETURN 
END IF
!
CALL IO_READ_FIELD(TZFMFILE,'IMAX',IIMAX)
CALL IO_READ_FIELD(TZFMFILE,'JMAX',IJMAX)
!
ALLOCATE(ZZS(IIMAX+2*JPHEXT,IJMAX+2*JPHEXT))
CALL IO_READ_FIELD(TZFMFILE,'ZS',ZZS)
!
!*            2. Allocate coarse arrays
!                ----------------------
!
IIMAXC = NIMAX/IDXRATIO
IJMAXC = NJMAX/IDYRATIO
!
IF ( MOD(NIMAX,IDXRATIO) /= 0 .OR. MOD(NJMAX,IDYRATIO) /= 0 ) THEN
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'Your truncated domain does not contain an integer number'
  WRITE (ILUOUT0,*) 'of grid meshes of its father domain; no nesting allowed'
  WRITE (ILUOUT0,*) 'NIMAX=',NIMAX,' IDXRATIO=',IDXRATIO
  WRITE (ILUOUT0,*) 'NJMAX=',NJMAX,' IDYRATIO=',IDYRATIO
  WRITE (ILUOUT0,*) '********************************************************'
  HDAD_NAME='                         '
  RETURN
END IF
!
IF ( MOD(KIINF-1,IDXRATIO) /= 0 .OR. MOD(KJINF-1,IDYRATIO) /= 0 ) THEN
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'Your truncated domain does not start on an grid mesh'
  WRITE (ILUOUT0,*) 'border of its father domain; no nesting allowed'
  WRITE (ILUOUT0,*) 'KIINF-1=',KIINF-1,' IDXRATIO=',IDXRATIO
  WRITE (ILUOUT0,*) 'KIINF-1=',KJINF-1,' IDYRATIO=',IDYRATIO
  WRITE (ILUOUT0,*) '********************************************************'
  HDAD_NAME='                         '
  RETURN
END IF
!
ALLOCATE(ZZS1(IIMAXC,IJMAXC))
ALLOCATE(ZZS2(IIMAXC,IJMAXC))
!
NXOR_ALL(1)=NXOR_ALL(1) + (KIINF-JPHEXT)/IDXRATIO
NYOR_ALL(1)=NYOR_ALL(1) + (KJINF-JPHEXT)/IDYRATIO
!
!*            3. Average the orographies
!                -----------------------
!
DO JI=1,IIMAXC
  DO JJ=1,IJMAXC
    ZZS1(JI,JJ) = SUM(ZZS( (JI-1)*IDXRATIO+KIINF+JPHEXT:JI*IDXRATIO+JPHEXT+KIINF-1,   &
                           (JJ-1)*IDYRATIO+KJINF+JPHEXT:JJ*IDYRATIO+JPHEXT+KJINF-1) ) &
                / (IDXRATIO*IDYRATIO)
    ZZS2(JI,JJ) = SUM(XZS( (JI-1)*IDXRATIO+1+JPHEXT:JI*IDXRATIO+JPHEXT,   &
                           (JJ-1)*IDYRATIO+1+JPHEXT:JJ*IDYRATIO+JPHEXT) ) &
                / (IDXRATIO*IDYRATIO)
  END DO
END DO
!
!*            4. Check the orographies
!                ----------------------
!
IF ( ANY(ABS(ZZS2(:,:)-ZZS1(:,:))>1.E-2) ) THEN
  HDAD_NAME='                         '
  WRITE (ILUOUT0,*) '********************************************************'
  WRITE (ILUOUT0,*) 'Maximum difference between averaged orographies: ', &
                    MAXVAL(ABS(ZZS2(:,:)-ZZS1(:,:)))
  WRITE (ILUOUT0,*) 'Averaged orographies are different; no nesting allowed'
  WRITE (ILUOUT0,*) '********************************************************'
END IF
!
!-------------------------------------------------------------------------------
!
IF (LSLEVE) THEN
!
  CALL IO_READ_FIELD(TZFMFILE,'ZSMT',ZZS)
!
!*            5. Average the smooth orographies
!                ------------------------------
!
  DO JI=1,IIMAXC
    DO JJ=1,IJMAXC
      ZZS1(JI,JJ) = SUM(ZZS( (JI-1)*IDXRATIO+KIINF+JPHEXT:JI*IDXRATIO+JPHEXT+KIINF-1,   &
                             (JJ-1)*IDYRATIO+KJINF+JPHEXT:JJ*IDYRATIO+JPHEXT+KJINF-1) ) &
                  / (IDXRATIO*IDYRATIO)
      ZZS2(JI,JJ) = SUM(XZSMT( (JI-1)*IDXRATIO+1+JPHEXT:JI*IDXRATIO+JPHEXT,   &
                               (JJ-1)*IDYRATIO+1+JPHEXT:JJ*IDYRATIO+JPHEXT) ) &
                  / (IDXRATIO*IDYRATIO)
    END DO
  END DO
!
!*            6. Check the smooth orographies
!                ----------------------------
!
  IF ( ANY(ABS(ZZS2(:,:)-ZZS1(:,:))>1.E-2) ) THEN
    HDAD_NAME='                         '
    WRITE (ILUOUT0,*) '********************************************************'
    WRITE (ILUOUT0,*) 'Maximum difference between averaged SMOOTH orographies: ', &
                      MAXVAL(ABS(ZZS2(:,:)-ZZS1(:,:)))
    WRITE (ILUOUT0,*) 'Averaged orographies are different; no nesting allowed'
    WRITE (ILUOUT0,*) '********************************************************'
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
IF (LEN_TRIM(HDAD_NAME)>0) THEN
  NDXRATIO_ALL(1)=IDXRATIO
  NDYRATIO_ALL(1)=IDYRATIO
END IF
!
!-------------------------------------------------------------------------------
!
WRITE (ILUOUT0,*) 'Routine CHECK_ZS completed'
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CHECK_ZS
