!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #####################
      MODULE MODI_PGD_INDEX
!     #####################
INTERFACE
      FUNCTION PGD_INDEX(HFIELD) RESULT(KINDEX)
!
CHARACTER(LEN=*),       INTENT(IN)    :: HFIELD  ! name of PGD field
INTEGER                               :: KINDEX  ! index of this field
!
END FUNCTION PGD_INDEX
END INTERFACE
END MODULE MODI_PGD_INDEX
!
!     #########################################
      FUNCTION PGD_INDEX(HFIELD) RESULT(KINDEX)
!     #########################################
!
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
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=*),       INTENT(IN)    :: HFIELD  ! name of PGD field
INTEGER                               :: KINDEX  ! index of this field
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
CHARACTER(LEN=20) :: YFIELD
!
!-------------------------------------------------------------------------------
!
YFIELD='                    '
YFIELD=HFIELD//YFIELD
!
SELECT CASE (YFIELD)
!
!*      1.     Vegetation parameters
!              ---------------------
!
!
  CASE('VEG                     ')
  KINDEX=1
!
  CASE('LAI                     ')
  KINDEX=2
!
  CASE('RSMIN                   ')
  KINDEX=3
!
  CASE('GAMMA                   ')
  KINDEX=4
!
  CASE('RGL                     ')
  KINDEX=5
!
  CASE('CV                      ')
  KINDEX=6
!
  CASE('DG2                     ')
  KINDEX=7
!
  CASE('Z0VEG                   ')
  KINDEX=8
!
  CASE('Z0HVEG                  ')
  KINDEX=9
!
  CASE('ALBNIR_ECO              ')
  KINDEX=10
!
  CASE('ALBVIS_ECO              ')
  KINDEX=11
!
  CASE('EMIS_VEG                ')
  KINDEX=12
!
  CASE('DG3                     ')
  KINDEX=13
!
  CASE('GMES                    ')
  KINDEX=15
!
  CASE('BSLAI                   ')
  KINDEX=16
!
  CASE('LAIMIN                   ')
  KINDEX=17
!
  CASE('SEFOLD                  ')
  KINDEX=18
!
  CASE('H_VEG                   ')
  KINDEX=19

!
!-------------------------------------------------------------------------------
!
!*      2.     Town parameters
!              ---------------
!
  CASE('Z0_TOWN                 ')
  KINDEX=20
!
  CASE('ALBNIR_ROOF             ')
  KINDEX=21
!
  CASE('ALBVIS_ROOF             ')
  KINDEX=22
!
  CASE('EMIS_ROOF               ')
  KINDEX=23
!
  CASE('HC_ROOF                 ')
  KINDEX=24
!
  CASE('ALBNIR_ROAD             ')
  KINDEX=25
!
  CASE('ALBVIS_ROAD             ')
  KINDEX=26
!
  CASE('EMIS_ROAD               ')
  KINDEX=27
!
  CASE('HC_ROAD                 ')
  KINDEX=28
!
  CASE('ALBNIR_WALL             ')
  KINDEX=29
!
  CASE('ALBVIS_WALL             ')
  KINDEX=30
!
  CASE('EMIS_WALL               ')
  KINDEX=31
!
  CASE('HC_WALL                 ')
  KINDEX=32
!
  CASE('BLD                     ')
  KINDEX=33
!
  CASE('BLD_HEIGHT              ')
  KINDEX=34
!
  CASE('BLD_HL_RATIO            ')
  KINDEX=35
!
  CASE('CAN_HW_RATIO            ')
  KINDEX=36
!
  CASE('TC_ROOF                 ')
  KINDEX=37
!
  CASE('TC_ROAD                 ')
  KINDEX=38
!
  CASE('TC_WALL                 ')
  KINDEX=39
!
  CASE('D_ROOF                  ')
  KINDEX=40
!
  CASE('D_ROAD                  ')
  KINDEX=41
!
  CASE('D_WALL                  ')
  KINDEX=42
!
  CASE('H_TRAFFIC               ')
  KINDEX=43
!
  CASE('LE_TRAFFIC              ')
  KINDEX=44
!
  CASE('H_INDUSTRY              ')
  KINDEX=45
!
  CASE('LE_INDUSTRY             ')
  KINDEX=46
!
!-------------------------------------------------------------------------------
!
  CASE DEFAULT
  KINDEX=100
!
END SELECT
!-------------------------------------------------------------------------------
!
END FUNCTION PGD_INDEX
