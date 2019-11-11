!MNH_LIC Copyright 1997-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_SELECT_STD_PGD
!     ##########################
INTERFACE
      SUBROUTINE SELECT_STD_PGD(HFIELD_NAME,PFIELD)
!
CHARACTER(LEN=*),     INTENT(IN) :: HFIELD_NAME ! pgd field name
REAL, DIMENSION(:,:), INTENT(IN) :: PFIELD      ! pgd field
!
END SUBROUTINE SELECT_STD_PGD
END INTERFACE
END MODULE MODI_SELECT_STD_PGD
!
!     ######################################
      SUBROUTINE SELECT_STD_PGD(HFIELD_NAME,PFIELD)
!     ######################################
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
!!    F.Solmon    06/00 patch approach : Rq
!!                      value of surface variable are atributed to NPT_USER class 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_LUNIT
USE MODD_PGDFIELDS
!
!
USE MODD_GROUND_PAR
!
!
USE MODI_PGD_INDEX
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=*),     INTENT(IN) :: HFIELD_NAME ! pgd field name
REAL, DIMENSION(:,:), INTENT(IN) :: PFIELD      ! pgd field
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!
CHARACTER(LEN=20) :: YFIELD
!
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
  LNOCLASS_PGD0(PGD_INDEX('VEG'))=.TRUE.
  XPGDVEG(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('LAI                     ')
  LNOCLASS_PGD0(PGD_INDEX('LAI'))=.TRUE.
  XPGDLAI(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('RSMIN                   ')
  LNOCLASS_PGD0(PGD_INDEX('RSMIN'))=.TRUE.
  XPGDRSMIN(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('GAMMA                   ')
  LNOCLASS_PGD0(PGD_INDEX('GAMMA'))=.TRUE.
  XPGDGAMMA(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('RGL                     ')
  LNOCLASS_PGD0(PGD_INDEX('RGL'))=.TRUE.
  XPGDRGL(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('CV                      ')
  LNOCLASS_PGD0(PGD_INDEX('CV'))=.TRUE.
  XPGDCV(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('DG2                     ')
  LNOCLASS_PGD0(PGD_INDEX('DG2'))=.TRUE.
  XPGDDG(:,:,2,NPT_USER) = PFIELD (:,:)

!
  CASE('DG3                     ')
  LNOCLASS_PGD0(PGD_INDEX('DG3'))=.TRUE.
  XPGDDG(:,:,3,NPT_USER) = PFIELD (:,:)

!
  CASE('Z0VEG                   ')
  LNOCLASS_PGD0(PGD_INDEX('Z0VEG'))=.TRUE.
  XPGDZ0VEG(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('Z0HVEG                  ')
  LNOCLASS_PGD0(PGD_INDEX('Z0HVEG'))=.TRUE.
  XPGDZ0HVEG(:,:,NPT_USER) = PFIELD (:,:)

  CASE('ALBNIR_ECO              ')
  LNOCLASS_PGD0(PGD_INDEX('ALBNIR_ECO'))=.TRUE.
  XPGDALBNIR_ECO(:,:,NPT_USER) = PFIELD (:,:)
!
  CASE('ALBVIS_ECO              ')
  LNOCLASS_PGD0(PGD_INDEX('ALBVIS_ECO'))=.TRUE.
  XPGDALBVIS_ECO(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('EMIS_ECO                ')
  LNOCLASS_PGD0(PGD_INDEX('EMIS_ECO'))=.TRUE.
  XPGDEMIS_ECO(:,:,NPT_USER) = PFIELD (:,:)

!
  CASE('GMES                    ')
  LNOCLASS_PGD0(PGD_INDEX('GMES'))=.TRUE.
  XPGDGMES(:,:,NPT_USER) =PFIELD (:,:)

!
  CASE('BSLAI                   ')
  LNOCLASS_PGD0(PGD_INDEX('BSLAI'))=.TRUE.
  XPGDBSLAI(:,:,NPT_USER) =PFIELD (:,:)

!
  CASE('LAIMIN                  ')
  LNOCLASS_PGD0(PGD_INDEX('LAIMIN'))=.TRUE.
  XPGDLAIMIN(:,:,NPT_USER) =PFIELD (:,:)

!
  CASE('SEFOLD                  ')
  LNOCLASS_PGD0(PGD_INDEX('SEFOLD'))=.TRUE.
  XPGDSEFOLD(:,:,NPT_USER) =PFIELD (:,:)

!
  CASE('H_TREE                  ')
  LNOCLASS_PGD0(PGD_INDEX('H_TREE'))=.TRUE.
  XPGDH_TREE(:,:,NPT_USER) =PFIELD (:,:)

!
!-------------------------------------------------------------------------------
!
!*      2.     Town parameters
!              ---------------
!
  CASE('Z0_TOWN                 ')
  LNOCLASS_PGD0(PGD_INDEX('Z0_TOWN'))=.TRUE.
  XPGDZ0_TOWN(:,:) = PFIELD (:,:)
!
  CASE('ALBNIR_ROOF             ')
  LNOCLASS_PGD0(PGD_INDEX('ALBNIR_ROOF'))=.TRUE.
  XPGDALBNIR_ROOF(:,:) = PFIELD (:,:)
!
  CASE('ALBVIS_ROOF             ')
  LNOCLASS_PGD0(PGD_INDEX('ALBVIS_ROOF'))=.TRUE.
  XPGDALBVIS_ROOF(:,:) = PFIELD (:,:)
!
  CASE('EMIS_ROOF               ')
  LNOCLASS_PGD0(PGD_INDEX('EMIS_ROOF'))=.TRUE.
  XPGDEMIS_ROOF(:,:) = PFIELD (:,:)
!
  CASE('HC_ROOF                 ')
  LNOCLASS_PGD0(PGD_INDEX('HC_ROOF'))=.TRUE.
  XPGDHC_ROOF(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDHC_ROOF,3))
!
  CASE('TC_ROOF                 ')
  LNOCLASS_PGD0(PGD_INDEX('TC_ROOF'))=.TRUE.
  XPGDTC_ROOF(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDTC_ROOF,3))
!
  CASE('D_ROOF                  ')
  LNOCLASS_PGD0(PGD_INDEX('D_ROOF'))=.TRUE.
  XPGDD_ROOF(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDD_ROOF,3))
!
  CASE('ALBNIR_ROAD             ')
  LNOCLASS_PGD0(PGD_INDEX('ALBNIR_ROAD'))=.TRUE.
  XPGDALBNIR_ROAD(:,:) = PFIELD (:,:)
!
  CASE('ALBVIS_ROAD             ')
  LNOCLASS_PGD0(PGD_INDEX('ALBVIS_ROAD'))=.TRUE.
  XPGDALBVIS_ROAD(:,:) = PFIELD (:,:)
!
  CASE('EMIS_ROAD               ')
  LNOCLASS_PGD0(PGD_INDEX('EMIS_ROAD'))=.TRUE.
  XPGDEMIS_ROAD(:,:) = PFIELD (:,:)
!
  CASE('HC_ROAD                 ')
  LNOCLASS_PGD0(PGD_INDEX('HC_ROAD'))=.TRUE.
  XPGDHC_ROAD(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDHC_ROAD,3))
!
  CASE('TC_ROAD                 ')
  LNOCLASS_PGD0(PGD_INDEX('TC_ROAD'))=.TRUE.
  XPGDTC_ROAD(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDTC_ROAD,3))
!
  CASE('D_ROAD                  ')
  LNOCLASS_PGD0(PGD_INDEX('D_ROAD'))=.TRUE.
  XPGDD_ROAD(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDD_ROAD,3))
!
  CASE('ALBNIR_WALL             ')
  LNOCLASS_PGD0(PGD_INDEX('ALBNIR_WALL'))=.TRUE.
  XPGDALBNIR_WALL(:,:) = PFIELD (:,:)
!
  CASE('ALBVIS_WALL             ')
  LNOCLASS_PGD0(PGD_INDEX('ALBVIS_WALL'))=.TRUE.
  XPGDALBVIS_WALL(:,:) = PFIELD (:,:)
!
  CASE('EMIS_WALL               ')
  LNOCLASS_PGD0(PGD_INDEX('EMIS_WALL'))=.TRUE.
  XPGDEMIS_WALL(:,:) = PFIELD (:,:)
!
  CASE('HC_WALL                 ')
  LNOCLASS_PGD0(PGD_INDEX('HC_WALL'))=.TRUE.
  XPGDHC_WALL(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDHC_WALL,3))
!
  CASE('TC_WALL                 ')
  LNOCLASS_PGD0(PGD_INDEX('TC_WALL'))=.TRUE.
  XPGDTC_WALL(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDTC_WALL,3))
!
  CASE('D_WALL                  ')
  LNOCLASS_PGD0(PGD_INDEX('D_WALL'))=.TRUE.
  XPGDD_WALL(:,:,:) = SPREAD(PFIELD (:,:),3,SIZE(XPGDD_WALL,3))
!
  CASE('BLD                     ')
  LNOCLASS_PGD0(PGD_INDEX('BLD'))=.TRUE.
  XPGDBLD(:,:) = PFIELD (:,:)
!
  CASE('BLD_HEIGHT              ')
  LNOCLASS_PGD0(PGD_INDEX('BLD_HEIGHT'))=.TRUE.
  XPGDBLD_HEIGHT(:,:) = PFIELD (:,:)
!
  CASE('BLD_HL_RATIO            ')
  LNOCLASS_PGD0(PGD_INDEX('BLD_HL_RATIO'))=.TRUE.
  XPGDBLD_HL_RATIO(:,:) = PFIELD (:,:)
!
  CASE('CAN_HW_RATIO            ')
  LNOCLASS_PGD0(PGD_INDEX('CAN_HW_RATIO'))=.TRUE.
  XPGDCAN_HW_RATIO(:,:) = PFIELD (:,:)
!
  CASE('H_TRAFFIC               ')
  LNOCLASS_PGD0(PGD_INDEX('H_TRAFFIC'))=.TRUE.
  XPGDH_TRAFFIC(:,:) = PFIELD (:,:)
!
  CASE('LE_TRAFFIC              ')
  LNOCLASS_PGD0(PGD_INDEX('LE_TRAFFIC'))=.TRUE.
  XPGDLE_TRAFFIC(:,:) = PFIELD (:,:)
!
  CASE('H_INDUSTRY               ')
  LNOCLASS_PGD0(PGD_INDEX('H_INDUSTRY'))=.TRUE.
  XPGDH_INDUSTRY(:,:) = PFIELD (:,:)
!
  CASE('LE_INDUSTRY              ')
  LNOCLASS_PGD0(PGD_INDEX('LE_INDUSTRY'))=.TRUE.
  XPGDLE_INDUSTRY(:,:) = PFIELD (:,:)
!
!-------------------------------------------------------------------------------
!
  CASE DEFAULT

  PRINT*, ' '
  PRINT*, 'The field ',HFIELD_NAME, ' is not yet a standard PGD field.'
  PRINT*, 'IT WILL NOT BE SAVED ON THE PGD FILE.'
  PRINT*, ' '

END SELECT
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SELECT_STD_PGD
