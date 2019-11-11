!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 12:43:28
!-----------------------------------------------------------------
!!    ############################ 
      MODULE MODD_SUB_ELEC_n
!!    ############################ 
!
!!****  *MODD_SUB_ELEC_n*- declaration for electric flashes 
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!       
!!    AUTHOR
!!    ------
!!	C.Barthe, C.Lac  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE SUB_ELEC_t
!
  LOGICAL :: GEFIRSTCALL  = .TRUE.
!
  INTEGER , DIMENSION(:), POINTER :: ISFLASH_NUMBER=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: ISNB_FLASH=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: ISCELL_NUMBER=>NULL() ! Cell number of the lightning flash
  INTEGER , DIMENSION(:), POINTER :: ISNBSEG=>NULL() ! Number of flash segments
  INTEGER , DIMENSION(:), POINTER :: ISTCOUNT_NUMBER=>NULL() ! Temporal loop number of the flash
  INTEGER , DIMENSION(:), POINTER :: ISTYPE=>NULL() ! flash type :IC, CGN or CGP
  REAL    , DIMENSION(:), POINTER :: ZXMASS=>NULL() ! Coord. at mass points
  REAL    , DIMENSION(:), POINTER :: ZYMASS=>NULL() ! Coord. at mass points
  REAL    , DIMENSION(:,:,:), POINTER :: ZZMASS=>NULL() ! Coord. at mass points
  REAL    , DIMENSION(:,:,:), POINTER :: ZPRES_COEF=>NULL() ! Pressure effect for E
  REAL    , DIMENSION(:,:,:), POINTER :: ZSCOORD_SEG=>NULL() ! Global coordinates of segments
  REAL    , DIMENSION(:), POINTER :: ZSEM_TRIG=>NULL() ! Electric field module at the triggering pt
  REAL    , DIMENSION(:), POINTER :: ZSNEUT_POS=>NULL() ! Positive charge
                             ! neutralized at each segment
  REAL    , DIMENSION(:), POINTER :: ZSNEUT_NEG=>NULL()    ! Negative charge       
                             ! neutralized at each segment
END TYPE SUB_ELEC_t

TYPE(SUB_ELEC_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: SUB_ELEC_MODEL

  LOGICAL, POINTER :: GEFIRSTCALL=>NULL()
  INTEGER , DIMENSION(:), POINTER :: ISFLASH_NUMBER=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: ISNB_FLASH=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: ISCELL_NUMBER=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: ISNBSEG=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: ISTCOUNT_NUMBER=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: ISTYPE=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZXMASS=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZYMASS=>NULL() 
  REAL    , DIMENSION(:,:,:), POINTER :: ZZMASS=>NULL() 
  REAL    , DIMENSION(:,:,:), POINTER :: ZPRES_COEF=>NULL() 
  REAL    , DIMENSION(:,:,:), POINTER :: ZSCOORD_SEG=>NULL()
  REAL    , DIMENSION(:), POINTER :: ZSEM_TRIG=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZSNEUT_POS=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZSNEUT_NEG=>NULL() 

CONTAINS

SUBROUTINE SUB_ELEC_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
SUB_ELEC_MODEL(KFROM)%ISFLASH_NUMBER=>ISFLASH_NUMBER
SUB_ELEC_MODEL(KFROM)%ISNB_FLASH=>ISNB_FLASH
SUB_ELEC_MODEL(KFROM)%ISCELL_NUMBER=>ISCELL_NUMBER
SUB_ELEC_MODEL(KFROM)%ISNBSEG=>ISNBSEG            
SUB_ELEC_MODEL(KFROM)%ISTCOUNT_NUMBER=>ISTCOUNT_NUMBER
SUB_ELEC_MODEL(KFROM)%ISTYPE=>ISTYPE           
SUB_ELEC_MODEL(KFROM)%ZXMASS=>ZXMASS
SUB_ELEC_MODEL(KFROM)%ZYMASS=>ZYMASS
SUB_ELEC_MODEL(KFROM)%ZZMASS=>ZZMASS
SUB_ELEC_MODEL(KFROM)%ZPRES_COEF=>ZPRES_COEF
SUB_ELEC_MODEL(KFROM)%ZSCOORD_SEG=>ZSCOORD_SEG
SUB_ELEC_MODEL(KFROM)%ZSEM_TRIG=>ZSEM_TRIG
SUB_ELEC_MODEL(KFROM)%ZSNEUT_POS=>ZSNEUT_POS
SUB_ELEC_MODEL(KFROM)%ZSNEUT_NEG=>ZSNEUT_NEG   
!
! Current model is set to model KTO
GEFIRSTCALL=>SUB_ELEC_MODEL(KTO)%GEFIRSTCALL
ISFLASH_NUMBER=>SUB_ELEC_MODEL(KTO)%ISFLASH_NUMBER
ISNB_FLASH=>SUB_ELEC_MODEL(KTO)%ISNB_FLASH 
ISCELL_NUMBER=>SUB_ELEC_MODEL(KTO)%ISCELL_NUMBER 
ISNBSEG=>SUB_ELEC_MODEL(KTO)%ISNBSEG        
ISTCOUNT_NUMBER=>SUB_ELEC_MODEL(KTO)%ISTCOUNT_NUMBER 
ISTYPE=>SUB_ELEC_MODEL(KTO)%ISTYPE              
ZXMASS=>SUB_ELEC_MODEL(KTO)%ZXMASS 
ZYMASS=>SUB_ELEC_MODEL(KTO)%ZYMASS 
ZZMASS=>SUB_ELEC_MODEL(KTO)%ZZMASS 
ZPRES_COEF=>SUB_ELEC_MODEL(KTO)%ZPRES_COEF 
ZSCOORD_SEG=>SUB_ELEC_MODEL(KTO)%ZSCOORD_SEG
ZSEM_TRIG=>SUB_ELEC_MODEL(KTO)%ZSEM_TRIG
ZSNEUT_POS=>SUB_ELEC_MODEL(KTO)%ZSNEUT_POS
ZSNEUT_NEG=>SUB_ELEC_MODEL(KTO)%ZSNEUT_NEG   

END SUBROUTINE SUB_ELEC_GOTO_MODEL

END MODULE MODD_SUB_ELEC_n
