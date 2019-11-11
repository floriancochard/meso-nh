!MNH_LIC Copyright 1997-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_SELECT_AVG_PGD
!     #############################
INTERFACE
      SUBROUTINE SELECT_AVG_PGD(HFIELD_NAME,HATYPE,HAREA)
!
CHARACTER(LEN=*), INTENT(IN)  :: HFIELD_NAME
CHARACTER(LEN=3), INTENT(OUT) :: HATYPE
CHARACTER(LEN=3), INTENT(OUT) :: HAREA
!
END SUBROUTINE SELECT_AVG_PGD
END INTERFACE
END MODULE MODI_SELECT_AVG_PGD
!
!     ###################################################
      SUBROUTINE SELECT_AVG_PGD(HFIELD_NAME,HATYPE,HAREA)
!     ###################################################
!
!!****  
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
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
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    15/12/97
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_LUNIT
USE MODE_MSG
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=*), INTENT(IN)  :: HFIELD_NAME
CHARACTER(LEN=3), INTENT(OUT) :: HATYPE
CHARACTER(LEN=3), INTENT(OUT) :: HAREA
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!
CHARACTER(LEN=20) :: YFIELD
!
INTEGER :: ILUOUT0, IRESP
!-------------------------------------------------------------------------------
!
YFIELD='                    '
YFIELD=HFIELD_NAME//YFIELD
!
SELECT CASE (YFIELD)
!
!*      1.     Vegetation parameters
!              ---------------------
!
!
  CASE('VEG                     ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('LAI                     ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('RSMIN                   ')
  HATYPE = 'INV'
  HAREA  = 'NAT'
!
  CASE('GAMMA                   ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('RGL                     ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('CV                      ')
  HATYPE = 'INV'
  HAREA  = 'NAT'
!
  CASE('DG2                     ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('DG3                     ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('Z0VEG                   ')
  HATYPE = 'LOG'
  HAREA  = 'NAT'
!
  CASE('Z0HVEG                  ')
  HATYPE = 'LOG'
  HAREA  = 'NAT'
!
  CASE('ALBNIR_ECO              ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('ALBVIS_ECO              ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('EMIS_ECO                ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('VEGTYPE                 ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('GMES                    ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('BSLAI                   ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('LAIMIN                  ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('SEFOLD                  ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
  CASE('H_VEG                   ')
  HATYPE = 'ARI'
  HAREA  = 'NAT'
!
!
!-------------------------------------------------------------------------------
!
!*      2.     Town parameters
!              ---------------
!
  CASE('Z0_TOWN                 ')
  HATYPE = 'LOG'
  HAREA  = 'TWN'
!
  CASE('ALBNIR_ROOF             ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('ALBVIS_ROOF             ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('EMIS_ROOF               ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('HC_ROOF                 ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('TC_ROOF                 ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('D_ROOF                  ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('ALBNIR_ROAD             ')
  HATYPE = 'ARI'
  HAREA  = 'STR'
!
  CASE('ALBVIS_ROAD             ')
  HATYPE = 'ARI'
  HAREA  = 'STR'
!
  CASE('EMIS_ROAD               ')
  HATYPE = 'ARI'
  HAREA  = 'STR'
!
  CASE('HC_ROAD                 ')
  HATYPE = 'ARI'
  HAREA  = 'STR'
!
  CASE('TC_ROAD                 ')
  HATYPE = 'ARI'
  HAREA  = 'STR'
!
  CASE('D_ROAD                  ')
  HATYPE = 'ARI'
  HAREA  = 'STR'
!
  CASE('ALBNIR_WALL             ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('ALBVIS_WALL             ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('EMIS_WALL               ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('HC_WALL                 ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('TC_WALL                 ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('D_WALL                  ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('BLD                     ')
  HATYPE = 'ARI'
  HAREA  = 'TWN'
!
  CASE('BLD_HEIGHT              ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('BLD_HL_RATIO            ')
  HATYPE = 'ARI'
  HAREA  = 'BLD'
!
  CASE('CAN_HW_RATIO            ')
  HATYPE = 'ARI'
  HAREA  = 'TWN'
!
  CASE('H_TRAFFIC               ')
  HATYPE = 'ARI'
  HAREA  = 'TWN'
!
  CASE('LE_TRAFFIC              ')
  HATYPE = 'ARI'
  HAREA  = 'TWN'
!
  CASE('H_INDUSTRY              ')
  HATYPE = 'ARI'
  HAREA  = 'TWN'
!
  CASE('LE_INDUSTRY             ')
  HATYPE = 'ARI'
  HAREA  = 'TWN'
!
!-------------------------------------------------------------------------------
!
  CASE DEFAULT
    ILUOUT0 = TLUOUT0%NLU
    WRITE(ILUOUT0) '*****************************************************************************'
    WRITE(ILUOUT0) ' '
    WRITE(ILUOUT0) 'Error in field name specification for a personal initialization:'
    WRITE(ILUOUT0) 'Field to be initialized by user is: ',YFIELD
    WRITE(ILUOUT0) ' '
    WRITE(ILUOUT0) 'MESONH available fields are: '
    WRITE(ILUOUT0) 'VEG          ','LAI          ','RSMIN      ','GAMMA   ','RGL     ','CV '
    WRITE(ILUOUT0) 'DG2          ','DG3          ','Z0VEG      ','Z0HVEG '
    WRITE(ILUOUT0) 'ALBNIR_ECO   ','ALBVIS_ECO   ','EMIS_ECO   '
    WRITE(ILUOUT0) 'VEGTYPE      ','GMES         ','BSLAI      ','LAIMIN  ','SEFOLD  ','H_VEG '
    WRITE(ILUOUT0) 'Z0_TOWN      ','BLD          ','BLD_HEIGHT '
    WRITE(ILUOUT0) 'BLD_HL_RATIO ','CAN_HW_RATIO '
    WRITE(ILUOUT0) 'ALBNIR_ROOF  ','ALBVIS_ROOF  ','EMIS_ROOF  ','HC_ROOF ','TC_ROOF ','D_ROOF '
    WRITE(ILUOUT0) 'ALBNIR_ROAD  ','ALBVIS_ROAD  ','EMIS_ROAD  ','HC_ROAD ','TC_ROAD ','D_ROAD '
    WRITE(ILUOUT0) 'ALBNIR_WALL  ','ALBVIS_WALL  ','EMIS_WALL  ','HC_WALL ','TC_WALL ','D_WALL '
    WRITE(ILUOUT0) 'H_TRAFFIC    ','LE_TRAFFIC   ','H_INDUSTRY ','LE_INDUSTRY '
    WRITE(ILUOUT0) ' '
    WRITE(ILUOUT0) 'Check the name if the field you want to initialize is one of these,'
    WRITE(ILUOUT0) 'or define your field as a dummy field if it is a new field'
    WRITE(ILUOUT0) ' '
    WRITE(ILUOUT0) '*****************************************************************************'
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SELECT_AVG_PGD','')
!
END SELECT
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SELECT_AVG_PGD
