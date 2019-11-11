!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #########################
MODULE MODI_INI_POSPROFILER_n
!      #########################
!
INTERFACE
!
      SUBROUTINE INI_POSPROFILER_n(PTSTEP, TPDTSEG, PSEGLEN, &
                                   KRR, KSV, OUSETKE,        &
                                   PLATOR, PLONOR            )
!
USE MODD_TYPE_DATE
!
REAL,               INTENT(IN) :: PTSTEP  ! time step
TYPE(DATE_TIME),    INTENT(IN) :: TPDTSEG ! segment date and time
REAL,               INTENT(IN) :: PSEGLEN ! segment length
INTEGER,            INTENT(IN) :: KRR     ! number of moist variables
INTEGER,            INTENT(IN) :: KSV     ! number of scalar variables
LOGICAL,            INTENT(IN) :: OUSETKE ! flag to use tke
REAL,               INTENT(IN) :: PLATOR  ! latitude of origine point
REAL,               INTENT(IN) :: PLONOR  ! longitude of origine point
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_POSPROFILER_n
!
END INTERFACE
!
END MODULE MODI_INI_POSPROFILER_n
!
!     ########################################################
      SUBROUTINE INI_POSPROFILER_n(PTSTEP, TPDTSEG, PSEGLEN, &
                                   KRR, KSV, OUSETKE,        &
                                   PLATOR, PLONOR            )
!     ########################################################
!
!
!!****  *INI_POSPROFILER_n* - 
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
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
!!    AUTHOR
!!    ------
!!      Valery Masson             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     P. Tulet 15/01/2002  
!!     C.Lac 10/2016  Add visibility diagnostic
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CONF
USE MODD_DYN_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_LUNIT_n,      ONLY: TLUOUT
USE MODD_PARAMETERS
USE MODD_PROFILER_n
USE MODD_RADIATIONS_n, ONLY: NAER
USE MODD_TYPE_PROFILER
USE MODD_TYPE_DATE
!
USE MODE_GRIDPROJ
USE MODE_IO_ll
USE MODE_ll
USE MODE_MSG
!
USE MODI_INI_PROFILER_N
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
REAL,               INTENT(IN) :: PTSTEP  ! time step
TYPE(DATE_TIME),    INTENT(IN) :: TPDTSEG ! segment date and time
REAL,               INTENT(IN) :: PSEGLEN ! segment length
INTEGER,            INTENT(IN) :: KRR     ! number of moist variables
INTEGER,            INTENT(IN) :: KSV     ! number of scalar variables
LOGICAL,            INTENT(IN) :: OUSETKE ! flag to use tke
REAL,               INTENT(IN) :: PLATOR  ! latitude of origine point
REAL,               INTENT(IN) :: PLONOR  ! longitude of origine point
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER :: ISTORE   ! number of storage instants
INTEGER :: ILUOUT   ! logical unit
INTEGER :: IKU      !
!
!----------------------------------------------------------------------------
ILUOUT = TLUOUT%NLU
!----------------------------------------------------------------------------
!
!*      1.   Default values
!            --------------
IKU =  SIZE(XZZ,3)     ! nombre de niveaux verticaux
!
CALL DEFAULT_PROFILER_n(TPROFILER)
!
!
!*      3.   Stations initialization
!            -----------------------
!
CALL INI_PROFILER_n
LPROFILER = (NUMBPROFILER>0)
!
!----------------------------------------------------------------------------
!
!*      4.   Allocations of storage arrays
!            -----------------------------
!
IF(NUMBPROFILER>0) THEN
  CALL ALLOCATE_PROFILER_n(TPROFILER)
  CALL INI_INTERP_PROFILER_n(TPROFILER)
END IF
!----------------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE DEFAULT_PROFILER_n(TPROFILER)
!
TYPE(PROFILER), INTENT(INOUT) :: TPROFILER
!
NUMBPROFILER     = 0
TPROFILER%T_CUR  = XUNDEF
TPROFILER%N_CUR  = 0
TPROFILER%STEP   = XTSTEP 
!
END SUBROUTINE DEFAULT_PROFILER_n
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ALLOCATE_PROFILER_n(TPROFILER)
!
TYPE(PROFILER), INTENT(INOUT) :: TPROFILER
!
ISTORE = INT ( (PSEGLEN-XTSTEP) / TPROFILER%STEP ) + 1
!
ALLOCATE(TPROFILER%TIME  (ISTORE))
ALLOCATE(TPROFILER%ERROR (NUMBPROFILER))
ALLOCATE(TPROFILER%X     (NUMBPROFILER))
ALLOCATE(TPROFILER%Y     (NUMBPROFILER))
ALLOCATE(TPROFILER%ZON   (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%MER   (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%FF    (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%DD    (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%W     (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%P     (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%ZZ    (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%TH    (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%THV   (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%RHOD  (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%VISI  (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%VISIKUN(ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%RARE  (ISTORE,IKU,NUMBPROFILER))
ALLOCATE(TPROFILER%R     (ISTORE,IKU,NUMBPROFILER,KRR))
ALLOCATE(TPROFILER%SV    (ISTORE,IKU,NUMBPROFILER,KSV))
ALLOCATE(TPROFILER%AER   (ISTORE,IKU,NUMBPROFILER,NAER))
IF (OUSETKE) THEN
  ALLOCATE(TPROFILER%TKE (ISTORE,IKU,NUMBPROFILER))
ELSE
  ALLOCATE(TPROFILER%TKE (0,IKU,0))
END IF
ALLOCATE(TPROFILER%DATIME(16,ISTORE))
ALLOCATE(TPROFILER%T2M     (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%Q2M     (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%HU2M    (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%ZON10M  (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%MER10M  (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%RN      (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%H       (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%LE      (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%LEI     (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%GFLUX   (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%LWD     (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%LWU     (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%SWD     (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%SWU     (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%IWV   (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%ZTD   (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%ZWD   (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%ZHD   (ISTORE,NUMBPROFILER))
ALLOCATE(TPROFILER%TKE_DISS(ISTORE,IKU,NUMBPROFILER))
!
!
TPROFILER%ERROR= .FALSE.
TPROFILER%TIME = XUNDEF
TPROFILER%ZON  = XUNDEF
TPROFILER%MER  = XUNDEF
TPROFILER%FF   = XUNDEF
TPROFILER%DD   = XUNDEF
TPROFILER%W    = XUNDEF
TPROFILER%P    = XUNDEF
TPROFILER%ZZ   = XUNDEF
TPROFILER%TH   = XUNDEF
TPROFILER%THV  = XUNDEF
TPROFILER%RHOD = XUNDEF
TPROFILER%VISI = XUNDEF
TPROFILER%VISIKUN = XUNDEF
TPROFILER%RARE = XUNDEF
TPROFILER%IWV  = XUNDEF
TPROFILER%ZTD  = XUNDEF
TPROFILER%ZWD  = XUNDEF
TPROFILER%ZHD  = XUNDEF
TPROFILER%R    = XUNDEF
TPROFILER%SV   = XUNDEF
TPROFILER%AER  = XUNDEF
TPROFILER%TKE  = XUNDEF
TPROFILER%T2M      = XUNDEF
TPROFILER%Q2M      = XUNDEF
TPROFILER%HU2M     = XUNDEF
TPROFILER%ZON10M   = XUNDEF
TPROFILER%MER10M   = XUNDEF
TPROFILER%RN       = XUNDEF
TPROFILER%H        = XUNDEF
TPROFILER%LE       = XUNDEF
TPROFILER%GFLUX    = XUNDEF
TPROFILER%LEI      = XUNDEF
TPROFILER%LWD      = XUNDEF
TPROFILER%LWU      = XUNDEF
TPROFILER%SWD      = XUNDEF
TPROFILER%SWU      = XUNDEF
TPROFILER%TKE_DISS = XUNDEF
!
END SUBROUTINE ALLOCATE_PROFILER_n
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE INI_INTERP_PROFILER_n(TPROFILER)
!
TYPE(PROFILER), INTENT(INOUT) :: TPROFILER
INTEGER :: III
INTEGER :: IIU, IJU
!
DO III=1,NUMBPROFILER
  CALL GET_DIM_EXT_ll ('B',IIU,IJU)
  CALL SM_XYHAT(PLATOR,PLONOR,                          &
                TPROFILER%LAT(III), TPROFILER%LON(III), &
                TPROFILER%X(III),   TPROFILER%Y(III)    )
ENDDO
!
IF ( ANY(TPROFILER%LAT(:)==XUNDEF) .OR. ANY(TPROFILER%LON(:)==XUNDEF) ) THEN
  WRITE(ILUOUT,*) 'Error in station position '
  WRITE(ILUOUT,*) 'either LATitude or LONgitude segment'
  WRITE(ILUOUT,*) 'definiton is not complete.'
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_INTERP_PROFILER_n','')
END IF
!
TPROFILER%STEP  = MAX ( PTSTEP, TPROFILER%STEP )
!
!
END SUBROUTINE INI_INTERP_PROFILER_n
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
END SUBROUTINE INI_POSPROFILER_n
