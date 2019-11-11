!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODD_DEEP_CONVECTION_n
!     #############################
!
!!****  *MODD_DEEP_CONVECTION$n* - Contains convective tendencies 
!!
!!    PURPOSE
!!    -------
!!      Contains global convective tendencies and convective counter
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_DEEP_CONVECTION$n)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!      Modif       11/12/98 : add diagnostic variables
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

INTEGER, DIMENSION(:,:),POINTER :: NCOUNTCONV=>NULL()  ! convective counter used
                                                       ! to refresh the
                                                       ! convective tendencies
REAL, DIMENSION(:,:,:) ,POINTER :: XDTHCONV=>NULL()    ! convective TH tendency (K/s)
REAL, DIMENSION(:,:,:) ,POINTER :: XDRVCONV=>NULL()    ! convective r_v tendency (1/s)
REAL, DIMENSION(:,:,:) ,POINTER :: XDRCCONV=>NULL()    ! convective r_c tendency (1/s)
REAL, DIMENSION(:,:,:) ,POINTER :: XDRICONV=>NULL()    ! convective r_i tendency (1/s)
REAL, DIMENSION(:,:)   ,POINTER :: XPRCONV=>NULL()     ! total precipitation rate (m/s)
REAL, DIMENSION(:,:)   ,POINTER :: XPRSCONV=>NULL()    ! solid precipitation rate (m/s)
REAL, DIMENSION(:,:)   ,POINTER :: XPACCONV=>NULL()    ! accumulated convective
                                                       ! precipitation (m)
REAL, DIMENSION(:,:,:,:),POINTER ::XDSVCONV=>NULL()    ! tracer tendencies (1/s)
!diagnostic variables
REAL, DIMENSION(:,:,:) ,POINTER :: XUMFCONV=>NULL()    ! updraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:) ,POINTER :: XDMFCONV=>NULL()    ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:) ,POINTER :: XMFCONV=>NULL()     ! convective mass flux (kg/s m2)
REAL, DIMENSION(:,:,:) ,POINTER :: XPRLFLXCONV=>NULL() ! liquid precip flux (m/s)
REAL, DIMENSION(:,:,:) ,POINTER :: XPRSFLXCONV=>NULL() ! solid  precip flux (m/s)
REAL, DIMENSION(:,:)   ,POINTER :: XCAPE=>NULL()       ! CAPE (J/kg)
INTEGER,DIMENSION(:,:), POINTER :: NCLTOPCONV=>NULL()  ! convective cloud top level
INTEGER,DIMENSION(:,:), POINTER :: NCLBASCONV=>NULL()  ! convective cloud base level

REAL, DIMENSION(:,:) ,  POINTER :: XIC_RATE=>NULL()          ! IC lightning frequency
REAL, DIMENSION(:,:) ,  POINTER :: XCG_RATE=>NULL()          ! CG lightning frequency
REAL, DIMENSION(:,:) ,  POINTER :: XIC_TOTAL_NUMBER=>NULL()  ! Total number of IC
REAL, DIMENSION(:,:) ,  POINTER :: XCG_TOTAL_NUMBER=>NULL()  ! Total number of CG

CONTAINS

SUBROUTINE DEEP_CONVECTION_GOTO_MODEL(KFROM, KTO)
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
END SUBROUTINE DEEP_CONVECTION_GOTO_MODEL

END MODULE MODD_DEEP_CONVECTION_n
