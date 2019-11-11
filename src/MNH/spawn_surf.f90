!MNH_LIC Copyright 2004-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################## 
MODULE MODI_SPAWN_SURF
!########################
!
INTERFACE
!
      SUBROUTINE SPAWN_SURF (HINIFILE, HINIFILEPGD, TPOUTDATAFILE, OSPAWN_SURF)
!
USE MODD_IO_ll, ONLY: TFILEDATA

!
CHARACTER (LEN=*),      INTENT(IN) :: HINIFILE     ! Input file
CHARACTER (LEN=*),      INTENT(IN) :: HINIFILEPGD
TYPE(TFILEDATA),POINTER,INTENT(IN) :: TPOUTDATAFILE
LOGICAL,                INTENT(IN) :: OSPAWN_SURF  ! flag to spawn surface fields
!
END SUBROUTINE SPAWN_SURF
!
END INTERFACE
!
END MODULE MODI_SPAWN_SURF
!
!
!     #########################################################################
      SUBROUTINE SPAWN_SURF (HINIFILE, HINIFILEPGD, TPOUTDATAFILE, OSPAWN_SURF)
!     #########################################################################
!
!!****  *SPAWN_SURF * - subroutine to  call spawning of surface fields
!!
!!    PURPOSE
!!    -------
!!
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson       * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     01/2004  
!!  06/2016     (G.Delautier) phasage surfex 8
!!     P.Wautelet 08/07/2016 : removed MNH_NCWRIT define
!!  01/2018      (G.Delautier) SURFEX 8.1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CONF,         ONLY : NVERB
USE MODD_GRID_n,       ONLY : XZS
USE MODD_IO_ll,        ONLY : TFILEDATA
USE MODD_IO_SURF_MNH,  ONLY : COUTFILE
USE MODD_LUNIT,        ONLY : TPGDFILE, TOUTDATAFILE
USE MODD_MNH_SURFEX_n
USE MODD_NESTING,      ONLY : CMY_NAME, CDAD_NAME
USE MODD_PARAM_n,      ONLY : CSURF
USE MODD_TIME_n,       ONLY : TDTCUR
!
USE MODE_ll
USE MODE_FMWRIT
USE MODE_IO_ll
USE MODE_MODELN_HANDLER
USE MODE_PREP_CTL,     ONLY : PREP_CTL
!
USE MODI_INIT_PGD_SURF_ATM
USE MODI_MNHPUT_ZS_n
USE MODI_PREP_SURF_ATM
USE MODI_WRITE_DIAG_SURF_ATM_N
USE MODI_WRITE_HGRIDn
USE MODI_WRITE_PGD_SURF_ATM_N
USE MODI_WRITE_SURF_ATM_N
USE MODI_ZOOM_PGD_SURF_ATM
!
IMPLICIT NONE
!
!*       0.1.1  Declarations of global variables not declared in the modules :
!
!*       0.1.2  Declarations of dummy arguments :
!
CHARACTER (LEN=*),      INTENT(IN) :: HINIFILE     ! Input file
CHARACTER (LEN=*),      INTENT(IN) :: HINIFILEPGD
TYPE(TFILEDATA),POINTER,INTENT(IN) :: TPOUTDATAFILE
LOGICAL,                INTENT(IN) :: OSPAWN_SURF  ! flag to spawn surface fields
!
!*       0.1.3  Declarations of local variables :
!
INTEGER         :: IINFO_ll
TYPE (PREP_CTL) :: YLCTL
!  
!-------------------------------------------------------------------------------
!
!*       7.      Surface variables :
!
IF (CSURF=='EXTE') THEN
  IF (OSPAWN_SURF) THEN
    TPGDFILE     => TPOUTDATAFILE !Corresponding to file with CNAME = CMY_NAME(2)
    TOUTDATAFILE => TPOUTDATAFILE !Corresponding to file with CNAME = CMY_NAME(2)
    COUTFILE   = CMY_NAME(2)
    !* spawn PGD fields
    CALL GOTO_SURFEX(1)
    !JUAN REALZ
    !CALL GO_TOMODEL_ll(1,IINFO_ll)
    !CALL GOTO_MODEL(1)
    !JUAN REALZ
    CALL ZOOM_PGD_SURF_ATM(YSURF_CUR,'MESONH',HINIFILEPGD,'MESONH',TPGDFILE%CNAME,'MESONH')
    CALL MNHPUT_ZS_n
    !* writing of physiographic fields in the file
    ALLOCATE(YSURF_CUR%DUO%CSELECT(0))
    CALL WRITE_PGD_SURF_ATM_n(YSURF_CUR,'MESONH')
    !* rereading of physiographic fields and definition of prognostic fields
    CALL INIT_PGD_SURF_ATM(YSURF_CUR,'MESONH','PRE',HINIFILE,'MESONH',      &
                           TDTCUR%TDATE%YEAR, TDTCUR%TDATE%MONTH, &
                           TDTCUR%TDATE%DAY, TDTCUR%TIME          )
    CALL PREP_SURF_ATM(YSURF_CUR,'MESONH',HINIFILE,'MESONH',HINIFILEPGD,'MESONH',YLCTL)
    !* writing of all surface fields
    CALL WRITE_SURF_ATM_n(YSURF_CUR,'MESONH','ALL',.FALSE.)
    CALL WRITE_DIAG_SURF_ATM_n(YSURF_CUR,'MESONH','ALL')
    !
  ELSE
    CSURF='EXRM'
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SPAWN_SURF
