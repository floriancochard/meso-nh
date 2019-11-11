!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #########################
MODULE MODI_INI_SURFSTATION_n
!      #########################
!
INTERFACE
!
      SUBROUTINE INI_SURFSTATION_n(PTSTEP, TPDTSEG, PSEGLEN, &
                                   KRR, KSV, OUSETKE,        &
                                   PLATOR, PLONOR            )
!
USE MODD_TYPE_DATE
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
END SUBROUTINE INI_SURFSTATION_n
!
END INTERFACE
!
END MODULE MODI_INI_SURFSTATION_n
!
!     ########################################################
      SUBROUTINE INI_SURFSTATION_n(PTSTEP, TPDTSEG, PSEGLEN, &
                                   KRR, KSV, OUSETKE,        &
                                   PLATOR, PLONOR            )
!     ########################################################
!
!
!!****  *INI_SURFSTATION_n* - 
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
!!     A. Lemonsu 19/11/2002 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CONF
USE MODD_DYN_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAMETERS
USE MODD_STATION_n
USE MODD_TYPE_DATE
!
USE MODE_GRIDPROJ
USE MODE_IO_ll
USE MODE_ll
USE MODE_MSG
!
USE MODI_INI_STATION_N
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
INTEGER :: ISTORE ! number of storage instants
INTEGER :: ILUOUT ! logical unit
!
!----------------------------------------------------------------------------
ILUOUT = TLUOUT%NLU
!----------------------------------------------------------------------------
!
!*      1.   Default values
!            --------------
!
CALL DEFAULT_STATION_n(TSTATION)
!
!
!*      3.   Stations initialization
!            -----------------------
!
CALL INI_STATION_n
LSTATION = (NUMBSTAT>0)
!
!----------------------------------------------------------------------------
!
!*      4.   Allocations of storage arrays
!            -----------------------------
!
IF(NUMBSTAT>0) THEN
  CALL ALLOCATE_STATION_n(TSTATION)
  CALL INI_INTERP_STATION_n(TSTATION)
ENDIF
!----------------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------------
SUBROUTINE DEFAULT_STATION_n(TSTATION)
!
TYPE(STATION), INTENT(INOUT) :: TSTATION
!
NUMBSTAT   = 0
!
TSTATION%T_CUR  = XUNDEF
TSTATION%N_CUR  = 0
TSTATION%STEP   = XTSTEP  
!
END SUBROUTINE DEFAULT_STATION_n
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE ALLOCATE_STATION_n(TSTATION)
!
TYPE(STATION), INTENT(INOUT) :: TSTATION   ! 
!
ISTORE = INT ( (PSEGLEN-XTSTEP) / TSTATION%STEP ) + 1
!
!
!
ALLOCATE(TSTATION%TIME(ISTORE))
ALLOCATE(TSTATION%ERROR (NUMBSTAT))
ALLOCATE(TSTATION%X   (NUMBSTAT))
ALLOCATE(TSTATION%Y   (NUMBSTAT))
ALLOCATE(TSTATION%SV  (ISTORE,NUMBSTAT,KSV))
ALLOCATE(TSTATION%TSRAD (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%ZS  (NUMBSTAT))
ALLOCATE(TSTATION%DATIME(16,ISTORE))
ALLOCATE(TSTATION%ZON   (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%MER   (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%W     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%P     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%TH    (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%R     (ISTORE,NUMBSTAT,KRR))
IF (OUSETKE) THEN
  ALLOCATE(TSTATION%TKE (ISTORE,NUMBSTAT))
ELSE
  ALLOCATE(TSTATION%TKE (0,0))
END IF
ALLOCATE(TSTATION%T2M     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%Q2M     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%HU2M    (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%ZON10M  (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%MER10M  (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%RN      (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%H       (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%LE      (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%GFLUX   (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%LEI     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%SWD     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%SWU     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%LWD     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%LWU     (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%DSTAOD  (ISTORE,NUMBSTAT))
ALLOCATE(TSTATION%SFCO2   (ISTORE,NUMBSTAT))
!
TSTATION%ERROR = .FALSE.
TSTATION%TIME  = XUNDEF
TSTATION%ZON   = XUNDEF
TSTATION%MER   = XUNDEF
TSTATION%W     = XUNDEF
TSTATION%P     = XUNDEF
TSTATION%TH    = XUNDEF
TSTATION%R     = XUNDEF
TSTATION%SV    = XUNDEF
TSTATION%TKE   = XUNDEF
TSTATION%TSRAD = XUNDEF
TSTATION%ZS    = XUNDEF
TSTATION%T2M     = XUNDEF
TSTATION%Q2M     = XUNDEF
TSTATION%HU2M    = XUNDEF
TSTATION%ZON10M  = XUNDEF
TSTATION%MER10M  = XUNDEF
TSTATION%RN      = XUNDEF
TSTATION%H       = XUNDEF
TSTATION%LE      = XUNDEF
TSTATION%GFLUX   = XUNDEF
TSTATION%LEI     = XUNDEF
TSTATION%SWD     = XUNDEF
TSTATION%SWU     = XUNDEF
TSTATION%LWD     = XUNDEF
TSTATION%LWU     = XUNDEF
TSTATION%DSTAOD  = XUNDEF
TSTATION%SFCO2   = XUNDEF
!
END SUBROUTINE ALLOCATE_STATION_n
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE INI_INTERP_STATION_n(TSTATION)
!
TYPE(STATION), INTENT(INOUT) :: TSTATION   ! 
INTEGER :: JII                             ! 
INTEGER :: IIU, IJU                        !   
!
IF ( ALL(TSTATION%LAT(:)/=XUNDEF) .AND. ALL(TSTATION%LON(:)/=XUNDEF) ) THEN
 LSTATLAT = .TRUE.
 DO JII=1,NUMBSTAT
   CALL GET_DIM_EXT_ll ('B',IIU,IJU)
   CALL SM_XYHAT(PLATOR,PLONOR,                        &
                 TSTATION%LAT(JII), TSTATION%LON(JII), &
                 TSTATION%X(JII),   TSTATION%Y(JII)    )
 ENDDO
ELSE
 LSTATLAT = .FALSE.
 DO JII=1,NUMBSTAT
   TSTATION%X(JII) = XXHAT(TSTATION%I(JII))
   TSTATION%Y(JII) = XYHAT(TSTATION%I(JII))
   CALL GET_DIM_EXT_ll ('B',IIU,IJU)
   CALL SM_LATLON(PLATOR,PLONOR,                       &
                 TSTATION%X(JII),   TSTATION%Y(JII),   &
                 TSTATION%LAT(JII), TSTATION%LON(JII)  )
 ENDDO
END IF
!
IF ( ANY(TSTATION%LAT(:)==XUNDEF) .OR. ANY(TSTATION%LON(:)==XUNDEF) ) THEN
  WRITE(ILUOUT,*) 'Error in station position '
  WRITE(ILUOUT,*) 'either LATitude or LONgitude segment'
  WRITE(ILUOUT,*) 'or I and J segment'
  WRITE(ILUOUT,*) 'definition is not complete.'
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SURFSTATION_n','')
END IF
!
TSTATION%STEP  = MAX ( PTSTEP, TSTATION%STEP )
!
!
END SUBROUTINE INI_INTERP_STATION_n
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
END SUBROUTINE INI_SURFSTATION_n
