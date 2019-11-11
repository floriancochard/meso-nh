!MNH_LIC Copyright 1998-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_INI_LS
!     ######################
!
INTERFACE 
!
      SUBROUTINE INI_LS(TPINIFILE,HGETRVM,OLSOURCE,              &
           PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,              &
           PDRYMASSS,                                            &
           PLSUMM,PLSVMM,PLSWMM,PLSTHMM,PLSRVMM,PLSZWSMM, PDRYMASST,PLENG, &
           OSTEADY_DMASS)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),        INTENT(IN)    :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN)    :: HGETRVM   ! GET indicator
LOGICAL,                INTENT(IN)    :: OLSOURCE  ! Switch for the source term
! Larger Scale fields (source if OLSOURCE=T,  fields at time t-dt if OLSOURCE=F) :
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLSUM,PLSVM,PLSWM  ! Wind
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLSTHM,PLSRVM      ! Mass
REAL, DIMENSION(:,:),            INTENT(INOUT) :: PLSZWSM                ! sea wave
!if OLSOURCE=T : 
REAL,                   INTENT(INOUT), OPTIONAL :: PDRYMASSS             ! Md source
!Large Scale  fields at time t-dt (if OLSOURCE=T) : 
REAL, DIMENSION(:,:,:), INTENT(IN),    OPTIONAL :: PLSUMM,PLSVMM,PLSWMM  ! Wind
REAL, DIMENSION(:,:,:), INTENT(IN),    OPTIONAL :: PLSTHMM,PLSRVMM       ! Mass
REAL, DIMENSION(:,:),   INTENT(IN),    OPTIONAL  :: PLSZWSMM             ! Sea wave
REAL,                   INTENT(IN),    OPTIONAL :: PDRYMASST             ! Md(t)
REAL,                   INTENT(IN),    OPTIONAL :: PLENG                 ! Interpolation length
LOGICAL,                INTENT(IN),    OPTIONAL :: OSTEADY_DMASS         ! Md evolution logical switch
!
END SUBROUTINE INI_LS
!
END INTERFACE
!
END MODULE MODI_INI_LS
!
!
!     ############################################################
      SUBROUTINE INI_LS(TPINIFILE,HGETRVM,OLSOURCE,              &
           PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,              &
           PDRYMASSS,                                            &
           PLSUMM,PLSVMM,PLSWMM,PLSTHMM,PLSRVMM,PLSZWSMM,PDRYMASST,PLENG, &
           OSTEADY_DMASS)
!     ############################################################
!
!!****  *INI_LS* - routine to initialize  Larger Scale fields
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to read Larger Scale fields
!!    and to initialize  larger scale fields sources (if OLSOURCE=T)
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!      FMREAD   : to read data in LFIFM file
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CONF   : NVERB
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine INI_LS)
!!      
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        22/09/98
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 14/03/2019: correct ZWS when variable not present in file
!!
!! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CONF
USE MODD_FIELD_n, only: XZWS_DEFAULT
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_TIME ! for type DATE_TIME
!
USE MODE_FMREAD
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
TYPE(TFILEDATA),        INTENT(IN)    :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN)    :: HGETRVM   ! GET indicator
LOGICAL,                INTENT(IN)    :: OLSOURCE  ! Switch for the source term
! Larger Scale fields (source if OLSOURCE=T,  fields at time t-dt if OLSOURCE=F) :
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLSUM,PLSVM,PLSWM  ! Wind
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLSTHM,PLSRVM      ! Mass
REAL, DIMENSION(:,:),            INTENT(INOUT) :: PLSZWSM                ! sea wave
!if OLSOURCE=T : 
REAL,                   INTENT(INOUT), OPTIONAL :: PDRYMASSS             ! Md source
!Large Scale  fields at time t-dt (if OLSOURCE=T) : 
REAL, DIMENSION(:,:,:), INTENT(IN),    OPTIONAL :: PLSUMM,PLSVMM,PLSWMM  ! Wind
REAL, DIMENSION(:,:,:), INTENT(IN),    OPTIONAL :: PLSTHMM,PLSRVMM       ! Mass
REAL, DIMENSION(:,:),   INTENT(IN),    OPTIONAL :: PLSZWSMM              ! Sea wave
REAL,                   INTENT(IN),    OPTIONAL :: PDRYMASST             ! Md(t)
REAL,                   INTENT(IN),    OPTIONAL :: PLENG                 ! Interpolation length
LOGICAL,                INTENT(IN),    OPTIONAL :: OSTEADY_DMASS         ! Md evolution logical switch
!
!
!*       0.2   declarations of local variables
!
CHARACTER(LEN=15) :: YVAL
INTEGER           :: IRESP
!
!-------------------------------------------------------------------------------
!
!
!
!*       1.    SOME INITIALIZATIONS
!              --------------------
!
!NONE
!
!
!-------------------------------------------------------------------------------
!
!*       2.    READ LARGE SCALE FIELDS
!              -----------------------
!
CALL IO_READ_FIELD(TPINIFILE,'LSUM', PLSUM)
CALL IO_READ_FIELD(TPINIFILE,'LSVM', PLSVM)
CALL IO_READ_FIELD(TPINIFILE,'LSWM', PLSWM)
CALL IO_READ_FIELD(TPINIFILE,'LSTHM',PLSTHM)
CALL IO_READ_FIELD(TPINIFILE,'ZWS',PLSZWSM, IRESP)
!If the field ZWS is not in the file, set its value to XZWS_DEFAULT
!ZWS is present in files since MesoNH 5.4.2
IF ( IRESP/=0 ) THEN
  WRITE (YVAL,'( E15.8 )') XZWS_DEFAULT
  CALL PRINT_MSG(NVERB_WARNING,'IO','INI_LS','ZWS not found in file: using default value: '//TRIM(YVAL)//' m')
  PLSZWSM(:,:) = XZWS_DEFAULT
END IF
!
IF (HGETRVM == 'READ') THEN         ! LS-vapor                                    
  CALL IO_READ_FIELD(TPINIFILE,'LSRVM',PLSRVM)
ENDIF
!
IF (PRESENT(PLSZWSMM)) THEN
    PLSZWSM(:,:)= (PLSZWSM(:,:) - PLSZWSMM(:,:)) / PLENG
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE LARGE SCALE SOURCES
!              -------------------------------
!
! IN case of initialization of large scale source terms (OLSOURCE=T) :
! xxxM are large scale source terms 
! xxxMM are large scale fields at time t -dt 
IF (OLSOURCE) THEN 
  IF (PRESENT(PLSUMM).AND.PRESENT(PLSVMM).AND.PRESENT(PLSWMM).AND.PRESENT(PLSTHMM)) THEN 
    PLSUM(:,:,:) = (PLSUM(:,:,:) - PLSUMM(:,:,:))   / PLENG
    PLSVM(:,:,:) = (PLSVM(:,:,:) - PLSVMM(:,:,:))   / PLENG
    PLSWM(:,:,:) = (PLSWM(:,:,:) - PLSWMM(:,:,:))   / PLENG
    PLSTHM(:,:,:)= (PLSTHM(:,:,:) - PLSTHMM(:,:,:)) / PLENG
  ELSE
    CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LS','missing argument')
  ENDIF
!  LS-vapor  
  IF (HGETRVM == 'READ') THEN        
    IF (PRESENT(PLSRVMM))   THEN                           
      PLSRVM(:,:,:)  = (PLSRVM(:,:,:) - PLSRVMM(:,:,:)) / PLENG
    ELSE
      CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LS','missing argument (PLSRVMM)')
    ENDIF
  ENDIF
! Dry mass
   IF(.NOT. OSTEADY_DMASS) THEN
     IF (PRESENT(PDRYMASSS).AND.PRESENT(PDRYMASST)) THEN
       CALL IO_READ_FIELD(TPINIFILE,'DRYMASST',PDRYMASSS)
       PDRYMASSS   = (PDRYMASSS - PDRYMASST) / PLENG
     ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE INI_LS
