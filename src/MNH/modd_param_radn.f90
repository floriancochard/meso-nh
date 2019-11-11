!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/11/23 17:28:26
!-----------------------------------------------------------------
!     ######################## 
      MODULE MODD_PARAM_RAD_n
!     ########################
!
!!****  *MODD_PARAM_RAD$n* - declaration of the control parameters for
!!                           calling the radiation scheme
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to define the set of space
!!    and time control parameters for the radiation computations.
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PARAM_RAD$n)
!!
!!    AUTHOR
!!    ------
!!     J.-P. Pinty   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      26/02/95
!!      J.Stein       2/1/96    add the control parameter for the split of 
!!                              the memory  
!!      J.-P. Pinty   19/11/96    redefine the SPLIT parameter
!!      J.Stein       17/7/99   add the LRAD_DIAG switch
!!      F.Solmon      15/03/02  add the control parameter for aerosol and cloud radiative 
!!                              properties. Remove the NSPOT option.
!!      B.Aouizerats  07/11     add aerosol optical properties CAOP
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE PARAM_RAD_t
!
  REAL   :: XDTRAD        ! Interval of time between two full radiation
                             ! computations.
  REAL   :: XDTRAD_CLONLY ! Interval of time between two radiation
                             ! computations for the cloudy columns only.
  REAL  :: XFUDG          ! subgrid cloud inhomogeneity factor

  LOGICAL:: LCLEAR_SKY    ! .TRUE. the radiation computations are made for
                             ! one mean clear-sky but the whole cloudy columns
                             ! .FALSE. the radiation computations are performed
                             ! for each individual vertical column.
  LOGICAL:: LINIRAD       ! logical switch to initialize or read in the FM
                             ! file the radiative tendencies, surface fluxes
                             ! and temporal informations
  INTEGER:: NRAD_COLNBR   ! Maximal number of air columns called in a single
                             ! call of the ECMWF_radiation subroutine
  INTEGER:: NRAD_DIAG     ! index for the number of fields in the output file
!
!
  CHARACTER (LEN=4) :: CAER  ! aerosols optical thickness climatology
  CHARACTER (LEN=4) :: CAOP  ! aerosols optical properties                       
  CHARACTER (LEN=4) :: CLW  ! LW code 
  CHARACTER (LEN=4) :: CEFRADL ! parametrisation of liquid effective radius
  CHARACTER (LEN=4) :: CEFRADI ! parametrisation of ice  effective radius
  CHARACTER (LEN=4) :: COPWLW ! parametrisation of LW optical properties for cloud water 
  CHARACTER (LEN=4) :: COPILW ! parametrisation of LW optical properties for cloud ice
  CHARACTER (LEN=4) :: COPWSW ! parametrisation of SW optical properties for cloud water 
  CHARACTER (LEN=4) :: COPISW ! parametrisation of SW optical properties for cloud ice
  LOGICAL           :: LAERO_FT ! logical switch for a temporal interpolation of
                                ! aerosol and ozone distribution
  LOGICAL           :: LFIX_DAT ! logical switch to fix the date to a constant 
                                ! perpetual day ( diurnal cycle is conserved)
!-------------------------------------------------------------------------------
!







END TYPE PARAM_RAD_t

TYPE(PARAM_RAD_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: PARAM_RAD_MODEL

REAL, POINTER :: XDTRAD=>NULL()
REAL, POINTER :: XDTRAD_CLONLY=>NULL()
REAL, POINTER :: XFUDG=>NULL()
LOGICAL, POINTER :: LCLEAR_SKY=>NULL()
LOGICAL, POINTER :: LINIRAD=>NULL()
INTEGER, POINTER :: NRAD_COLNBR=>NULL()
INTEGER, POINTER :: NRAD_DIAG=>NULL()
CHARACTER (LEN=4), POINTER :: CAER=>NULL()
CHARACTER (LEN=4), POINTER :: CAOP=>NULL()
CHARACTER (LEN=4), POINTER :: CLW=>NULL()
CHARACTER (LEN=4), POINTER :: CEFRADL=>NULL()
CHARACTER (LEN=4), POINTER :: CEFRADI=>NULL()
CHARACTER (LEN=4), POINTER :: COPWLW=>NULL()
CHARACTER (LEN=4), POINTER :: COPILW=>NULL()
CHARACTER (LEN=4), POINTER :: COPWSW=>NULL()
CHARACTER (LEN=4), POINTER :: COPISW=>NULL()
LOGICAL, POINTER :: LAERO_FT=>NULL()
LOGICAL, POINTER :: LFIX_DAT=>NULL()

CONTAINS

SUBROUTINE PARAM_RAD_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
XDTRAD=>PARAM_RAD_MODEL(KTO)%XDTRAD
XDTRAD_CLONLY=>PARAM_RAD_MODEL(KTO)%XDTRAD_CLONLY
XFUDG=>PARAM_RAD_MODEL(KTO)%XFUDG
LCLEAR_SKY=>PARAM_RAD_MODEL(KTO)%LCLEAR_SKY
LINIRAD=>PARAM_RAD_MODEL(KTO)%LINIRAD
NRAD_COLNBR=>PARAM_RAD_MODEL(KTO)%NRAD_COLNBR
NRAD_DIAG=>PARAM_RAD_MODEL(KTO)%NRAD_DIAG
CAER=>PARAM_RAD_MODEL(KTO)%CAER
CAOP=>PARAM_RAD_MODEL(KTO)%CAOP
CLW=>PARAM_RAD_MODEL(KTO)%CLW
CEFRADL=>PARAM_RAD_MODEL(KTO)%CEFRADL
CEFRADI=>PARAM_RAD_MODEL(KTO)%CEFRADI
COPWLW=>PARAM_RAD_MODEL(KTO)%COPWLW
COPILW=>PARAM_RAD_MODEL(KTO)%COPILW
COPWSW=>PARAM_RAD_MODEL(KTO)%COPWSW
COPISW=>PARAM_RAD_MODEL(KTO)%COPISW
LAERO_FT=>PARAM_RAD_MODEL(KTO)%LAERO_FT
LFIX_DAT=>PARAM_RAD_MODEL(KTO)%LFIX_DAT

END SUBROUTINE PARAM_RAD_GOTO_MODEL

END MODULE MODD_PARAM_RAD_n
