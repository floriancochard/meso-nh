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
      MODULE MODD_SUB_PASPOL_n
!!    ############################ 
!
!!****  *MODD_SUB_PASPOL_n*- declaration of passive pollutants
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!      the kinetics of the releases for passive pollutants
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
!!	C. Lac     *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/08
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE SUB_PASPOL_t
!
  LOGICAL :: GPPFIRSTCALL  = .TRUE.
!
  LOGICAL , DIMENSION(:), POINTER :: GBEGEMIS=>NULL()! flag to initialize of emission
  INTEGER , DIMENSION(:), POINTER :: IDEBYY=>NULL() ! Begin of release YYYY
  INTEGER , DIMENSION(:), POINTER :: IDEBMM=>NULL() ! Begin of release MM  
  INTEGER , DIMENSION(:), POINTER :: IDEBDD=>NULL() ! Begin of release DD  
  REAL    , DIMENSION(:), POINTER :: ZDEBSS=>NULL() ! Begin of release SS
  REAL    , DIMENSION(:), POINTER :: ZCHRO2=>NULL() ! Time of increase (sec)
  REAL    , DIMENSION(:), POINTER :: ZCHRO3=>NULL() ! Time of decrease(sec)
  REAL    , DIMENSION(:), POINTER :: ZCHRO4=>NULL() ! End of release (sec)
  INTEGER , DIMENSION(:), POINTER :: IPIGI=>NULL()  ! Index on X of the release
  INTEGER , DIMENSION(:), POINTER :: IPIGJ=>NULL()  ! Index on Y of the release
  REAL    , DIMENSION(:), POINTER :: ZRATEINC=>NULL() ! Rate of release in the
                                                      ! first stage
  REAL    , DIMENSION(:), POINTER :: ZRATECST=>NULL() ! Constant rate of release
  REAL    , DIMENSION(:), POINTER :: ZRATEDEC=>NULL() ! Rate of release in the
                                                      ! last stage
  REAL    , DIMENSION(:), POINTER :: ZQTOT=>NULL()    ! Total mass released   
  CHARACTER(LEN=3), DIMENSION(:), POINTER :: CNIVO=>NULL() ! Type of emission :
                                                           ! SRF ou ALT
  REAL    , DIMENSION(:,:,:), POINTER :: ZDIFHORZ=>NULL()  ! Horizontal
                                        ! diffusion of the release
  REAL    , DIMENSION(:,:), POINTER ::  ZREPVERT=>NULL()  ! Column of emission
END TYPE SUB_PASPOL_t

TYPE(SUB_PASPOL_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: SUB_PASPOL_MODEL

  LOGICAL, POINTER :: GPPFIRSTCALL=>NULL()
  LOGICAL,  DIMENSION(:), POINTER :: GBEGEMIS=>NULL()
  INTEGER , DIMENSION(:), POINTER :: IDEBYY=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: IDEBMM=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: IDEBDD=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZDEBSS=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZCHRO2=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZCHRO3=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZCHRO4=>NULL() 
  INTEGER , DIMENSION(:), POINTER :: IPIGI=>NULL()  
  INTEGER , DIMENSION(:), POINTER :: IPIGJ=>NULL()  
  REAL    , DIMENSION(:), POINTER :: ZRATEINC=>NULL()
  REAL    , DIMENSION(:), POINTER :: ZRATECST=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZRATEDEC=>NULL() 
  REAL    , DIMENSION(:), POINTER :: ZQTOT=>NULL() 
  CHARACTER(LEN=3), DIMENSION(:), POINTER :: CNIVO=>NULL() 
  REAL    , DIMENSION(:,:,:), POINTER :: ZDIFHORZ=>NULL()  
  REAL    , DIMENSION(:,:), POINTER ::  ZREPVERT=>NULL()

CONTAINS

SUBROUTINE SUB_PASPOL_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
SUB_PASPOL_MODEL(KFROM)%GBEGEMIS=>GBEGEMIS
SUB_PASPOL_MODEL(KFROM)%IDEBYY=>IDEBYY
SUB_PASPOL_MODEL(KFROM)%IDEBMM=>IDEBMM
SUB_PASPOL_MODEL(KFROM)%IDEBDD=>IDEBDD
SUB_PASPOL_MODEL(KFROM)%ZDEBSS=>ZDEBSS
SUB_PASPOL_MODEL(KFROM)%ZCHRO2=>ZCHRO2
SUB_PASPOL_MODEL(KFROM)%ZCHRO3=>ZCHRO3
SUB_PASPOL_MODEL(KFROM)%ZCHRO4=>ZCHRO4
SUB_PASPOL_MODEL(KFROM)%IPIGI=>IPIGI 
SUB_PASPOL_MODEL(KFROM)%IPIGJ=>IPIGJ 
SUB_PASPOL_MODEL(KFROM)%ZRATEINC=>ZRATEINC
SUB_PASPOL_MODEL(KFROM)%ZRATECST=>ZRATECST
SUB_PASPOL_MODEL(KFROM)%ZRATEDEC=>ZRATEDEC
SUB_PASPOL_MODEL(KFROM)%ZQTOT=>ZQTOT   
SUB_PASPOL_MODEL(KFROM)%CNIVO=>CNIVO         
SUB_PASPOL_MODEL(KFROM)%ZDIFHORZ=>ZDIFHORZ      
SUB_PASPOL_MODEL(KFROM)%ZREPVERT=>ZREPVERT      
!
! Current model is set to model KTO
GPPFIRSTCALL=>SUB_PASPOL_MODEL(KTO)%GPPFIRSTCALL
GBEGEMIS=>SUB_PASPOL_MODEL(KTO)%GBEGEMIS   
IDEBYY=>SUB_PASPOL_MODEL(KTO)%IDEBYY
IDEBMM=>SUB_PASPOL_MODEL(KTO)%IDEBMM 
IDEBDD=>SUB_PASPOL_MODEL(KTO)%IDEBDD 
ZDEBSS=>SUB_PASPOL_MODEL(KTO)%ZDEBSS 
ZCHRO2=>SUB_PASPOL_MODEL(KTO)%ZCHRO2 
ZCHRO3=>SUB_PASPOL_MODEL(KTO)%ZCHRO3 
ZCHRO4=>SUB_PASPOL_MODEL(KTO)%ZCHRO4 
IPIGI=>SUB_PASPOL_MODEL(KTO)%IPIGI  
IPIGJ=>SUB_PASPOL_MODEL(KTO)%IPIGJ  
ZRATEINC=>SUB_PASPOL_MODEL(KTO)%ZRATEINC
ZRATECST=>SUB_PASPOL_MODEL(KTO)%ZRATECST
ZRATEDEC=>SUB_PASPOL_MODEL(KTO)%ZRATEDEC
ZQTOT=>SUB_PASPOL_MODEL(KTO)%ZQTOT   
CNIVO=>SUB_PASPOL_MODEL(KTO)%CNIVO 
ZDIFHORZ=>SUB_PASPOL_MODEL(KTO)%ZDIFHORZ  
ZREPVERT=>SUB_PASPOL_MODEL(KTO)%ZREPVERT  

END SUBROUTINE SUB_PASPOL_GOTO_MODEL

END MODULE MODD_SUB_PASPOL_n
