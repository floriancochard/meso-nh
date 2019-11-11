!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_CONF_n
!     #################
!
!!****  *MODD_CONF$n* - declaration of configuration variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the variables
!     which concern the configuration of the model. For exemple, 
!     the  type of moist variables. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_CONFn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!       
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94                      
!!      J.-P. Pinty 11/04/96  include the ice concentration
!!      D. Gazen    22/01/01  move NSV to MODD_NSV module
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CONF_t
  LOGICAL :: LUSERV  ! Logical to use rv
  LOGICAL :: LUSERC  ! Logical to use rc
  LOGICAL :: LUSERR  ! Logical to use rr
  LOGICAL :: LUSERI  ! Logical to use ri
  LOGICAL :: LUSERS  ! Logical to use rs
  LOGICAL :: LUSERG  ! Logical to use rg
  LOGICAL :: LUSERH  ! Logical to use rh
  INTEGER :: NRR     ! Total number of water variables
  INTEGER :: NRRL    ! Number of liquid water variables
  INTEGER :: NRRI    ! Number of solid water variables
!
  INTEGER :: IDX_RVT = -1 ! Position in array for rv
  INTEGER :: IDX_RCT = -1 ! Position in array for rc
  INTEGER :: IDX_RRT = -1 ! Position in array for rr
  INTEGER :: IDX_RIT = -1 ! Position in array for ri
  INTEGER :: IDX_RST = -1 ! Position in array for rs
  INTEGER :: IDX_RGT = -1 ! Position in array for rg
  INTEGER :: IDX_RHT = -1 ! Position in array for rh
!
  CHARACTER (LEN=2) :: CSTORAGE_TYPE ! storage type for the informations 
                                 ! written in the FM files ( 'TT' if the MesoNH 
                                 ! prognostic fields are at the same instant;
                                 ! 'MT' if they are taken at two instants in
                                 ! succession; 'PG' for PGD files informations )
  LOGICAL :: LUSECI  ! Logical to use Ci
!
END TYPE CONF_t

TYPE(CONF_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CONF_MODEL

LOGICAL, POINTER :: LUSERV=>NULL()
LOGICAL, POINTER :: LUSERC=>NULL()
LOGICAL, POINTER :: LUSERR=>NULL()
LOGICAL, POINTER :: LUSERI=>NULL()
LOGICAL, POINTER :: LUSERS=>NULL()
LOGICAL, POINTER :: LUSERG=>NULL()
LOGICAL, POINTER :: LUSERH=>NULL()
INTEGER, POINTER :: NRR=>NULL()
INTEGER, POINTER :: NRRL=>NULL()
INTEGER, POINTER :: NRRI=>NULL()
INTEGER, POINTER :: IDX_RVT=>NULL()
INTEGER, POINTER :: IDX_RCT=>NULL()
INTEGER, POINTER :: IDX_RRT=>NULL()
INTEGER, POINTER :: IDX_RIT=>NULL()
INTEGER, POINTER :: IDX_RST=>NULL()
INTEGER, POINTER :: IDX_RGT=>NULL()
INTEGER, POINTER :: IDX_RHT=>NULL()
LOGICAL, POINTER :: LUSECI=>NULL()
CHARACTER (LEN=2),POINTER :: CSTORAGE_TYPE=>NULL()

CONTAINS

SUBROUTINE CONF_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
LUSERV=>CONF_MODEL(KTO)%LUSERV
LUSERC=>CONF_MODEL(KTO)%LUSERC
LUSERR=>CONF_MODEL(KTO)%LUSERR
LUSERI=>CONF_MODEL(KTO)%LUSERI
LUSERS=>CONF_MODEL(KTO)%LUSERS
LUSERG=>CONF_MODEL(KTO)%LUSERG
LUSERH=>CONF_MODEL(KTO)%LUSERH
NRR=>CONF_MODEL(KTO)%NRR
NRRL=>CONF_MODEL(KTO)%NRRL
NRRI=>CONF_MODEL(KTO)%NRRI
IDX_RVT=>CONF_MODEL(KTO)%IDX_RVT
IDX_RCT=>CONF_MODEL(KTO)%IDX_RCT
IDX_RRT=>CONF_MODEL(KTO)%IDX_RRT
IDX_RIT=>CONF_MODEL(KTO)%IDX_RIT
IDX_RST=>CONF_MODEL(KTO)%IDX_RST
IDX_RGT=>CONF_MODEL(KTO)%IDX_RGT
IDX_RHT=>CONF_MODEL(KTO)%IDX_RHT
LUSECI=>CONF_MODEL(KTO)%LUSECI
CSTORAGE_TYPE=>CONF_MODEL(KTO)%CSTORAGE_TYPE

END SUBROUTINE CONF_GOTO_MODEL

END MODULE MODD_CONF_n
