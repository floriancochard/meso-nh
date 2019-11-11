!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!     ###############
      SUBROUTINE INIT_MNH
!     ###############
!
!!****  *INIT_MNH * - monitor to initialize the variables of the model
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize all the variables
!     used in the  model temporal loop or in the post-processings
!
!!**  METHOD
!!    ------
!!      This initialization  is separated in three parts :
!!         1. A part common to all models where :
!!            - The output-listing file common to all models is opened.
!!            - The physical constants are initialized.
!!            - The other constants for all models are initialized.
!!         2. The treatment of descriptor files model by model :
!!             The DESFM and EXSEG files are read and the EXSEG file is updated
!!         3. The sequential initialization of  nested models :
!!             The initial data fields are read in different files for each
!!             model and  variables which are not in these initial files are
!!             deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      INI_CST    : to initialize physical constants
!!      INI_CTURB  : to initialize for all models the constants used in the
!!                   turbulence scheme
!!      INI_SEG_n  : to read and update descriptor files  
!!      INI_SIZE   : to initialize the sizes of the different models
!!      INI_MODEL  : to initialize each nested model
!!      INI_PARA_ll: to build the ll data structures
!!      GO_TOMODEL : displace the ll lists to the right nested model
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : JPMODELMAX
!!
!!      Module MODD_CONF       : NMODEL,NVERB
!!
!!      Module MODD_LUNIT      : CLUOUT0
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INIT_MNH)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/06/94
!!      J.Stein     05/01/95  add ini_cturb
!!      J.P. Lafore 18/08/95  Time STEP change
!!      J.P. Lafore 22/07/96  ZTSTEP_ALL introduction for nesting
!!      V. Ducrocq  7/08/98   //
!!      P. Jabouille 7/07/99  split ini_modeln in 2 parts+ cleaning
!!      V. Masson   15/03/99  call to ini_data_cover
!!      P.Jabouille 15/07/99  special initialisation for spawning
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_LUNIT
USE MODD_LBC_n, ONLY: CLBCX,CLBCY   ! only for spawning purpose
USE MODD_DYN_n, ONLY: CPRESOPT,NITR ! only for spawning purpose
!
USE MODE_IO_ll
USE MODE_ll
USE MODE_FM
!
USE MODI_INI_CST
USE MODI_INI_CTURB
USE MODI_INI_SEG_n
USE MODI_INI_MODEL_n
USE MODI_INI_SIZE_n
USE MODI_INI_SIZE_SPAWN
USE MODI_RESET_EXSEG
USE MODE_MODELN_HANDLER
!
!JUAN
USE MODE_SPLITTINGZ_ll
!JUAN
!
IMPLICIT NONE
!
!*       0.1   Local variables
!
INTEGER :: JMI                                        !  Loop index
CHARACTER (LEN=16), DIMENSION(JPMODELMAX) :: YLUOUT   ! Name for output-listing
                                                      ! of nested models
CHARACTER (LEN=28), DIMENSION(JPMODELMAX) :: YINIFILE ! names of
                                                      ! the initial files
INTEGER  :: ILUOUT0,IRESP                             ! Logical unit number for
                                                      ! output-listing common
                                                      ! to all models and return
                                                      ! code of file management
REAL, DIMENSION(JPMODELMAX)            :: ZTSTEP_OLD  ! OLD Time STEP (DESFM)
REAL, DIMENSION(JPMODELMAX)            :: ZTSTEP_ALL  ! Time STEP of ALL models
INTEGER                                :: IINFO_ll    ! return code of // routines
!
! Dummy pointers needed to correct an ifort Bug
CHARACTER(LEN=4), DIMENSION(:), POINTER :: DPTR_CLBCX,DPTR_CLBCY

!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION  COMMON TO ALL MODELS
!              ------------------------------------
!
!*       1.1   initialize // E/S and open  output-listing file
!
!
IF (CPROGRAM/='REAL  ') THEN
  CLUOUT0 = 'OUTPUT_LISTING0'
  CALL OPEN_ll(UNIT=ILUOUT0,FILE=CLUOUT0,IOSTAT=IRESP,FORM='FORMATTED', &
               ACTION='WRITE',MODE=GLOBAL)
ELSE
  CALL FMLOOK_ll(CLUOUT0,CLUOUT0,ILUOUT0,IRESP)
END IF
!
WRITE(UNIT=ILUOUT0,FMT="(50('*'),/,'*',48X,'*',/,                  &
                  &  7('*'),10X, ' MESO-NH MODEL ',10X,8('*'),/,   &
                  & '*',48X,'*',/,                                 &
                  & 7('*'),12X,' CNRM - LA ',12X,8('*'),/,         &
                  & '*',48X,'*',/,  50('*'))")
!
!
!*      1.2   initialize physical constants
!
CALL INI_CST
!
!
!*       1.3    initialize other constants
!
CALL INI_CTURB                      ! constants for the turbulence scheme
!
!
!*       1.4    initialize surface parameters for each cover type
!
!
!-------------------------------------------------------------------------------
!
!*       2.    READ AND UPDATE DESCRIPTOR FILES
!              --------------------------------
!
DO JMI=1,JPMODELMAX
  CALL GOTO_MODEL(JMI)
  CALL INI_SEG_n(JMI,YLUOUT(JMI),YINIFILE(JMI),ZTSTEP_OLD(JMI),ZTSTEP_ALL)
  IF (JMI.EQ.NMODEL) EXIT
END DO
!
IF (CPROGRAM=='DIAG') CALL RESET_EXSEG(YLUOUT(1))
!
!-------------------------------------------------------------------------------
!
!
!*       3.    INITIALIZE EACH MODEL SIZES AND DEPENDENCY
!              ------------------------------------------
!
DO JMI=1,NMODEL
  CALL GOTO_MODEL(JMI)
  CALL INI_SIZE_n(JMI,YLUOUT(JMI),YINIFILE(JMI))
END DO
!
IF (CPROGRAM=='SPAWN ') THEN 
  DPTR_CLBCX=>CLBCX
  DPTR_CLBCY=>CLBCY
  CALL INI_SIZE_SPAWN(DPTR_CLBCX,DPTR_CLBCY,CPRESOPT,NITR,YINIFILE(1))
END IF
!
!   INITIALIZE data structures of ComLib
!
!JUAN CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!
!-------------------------------------------------------------------------------
!
!
!*       4.    INITIALIZE EACH MODEL
!              ---------------------
!
DO JMI=1,NMODEL
  CALL GO_TOMODEL_ll(JMI,IINFO_ll)
  CALL GOTO_MODEL(JMI)
  CALL INI_MODEL_n(JMI,ZTSTEP_OLD(JMI),YLUOUT(JMI),YINIFILE(JMI))
END DO
!
!-------------------------------------------------------------------------------
!
!*       5.    WRITE MESSAGE ON OUTPUT-LISTING
!              -------------------------------
!
IF (NVERB >= 5) THEN
  WRITE(UNIT=ILUOUT0,FMT="(50('*'),/,'*',48X,'*',/,                &
                & '*',10X,' INITIALIZATION TERMINATED',10X,'*',/,  &
                & '*',48X,'*',/,50('*'))")
END IF
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_MNH
