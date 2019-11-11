!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_REF_n
!     #################
!
!!****  *MODD_REF$n* - declaration of reference state  variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the variables
!     which concern the reference state, used for the anelastic 
!     approximation. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_REFn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94                      
!!      Modification   03/01/95  (Lafore)  To add the reference mass variables                      
!!      Modification   09/02/95  (Masson)  To add LINMASS
!!      Modification   25/07/97  (Stein)   To add XRVREF
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE REF_t
!          
!  REAL, DIMENSION(:,:,:), POINTER :: XRHODREF=>NULL()! dry density for ref. state 
!  REAL, DIMENSION(:,:,:), POINTER :: XTHVREF=>NULL() ! Thetav for ref. state 
  REAL, DIMENSION(:,:,:), POINTER :: XRVREF=>NULL()  ! mixing Ratio for Vapor computed
                                                    ! for the REFerence state
  REAL, DIMENSION(:,:,:), POINTER :: XEXNREF=>NULL() ! Exner function for ref. state 
!
  REAL, DIMENSION(:,:,:), POINTER :: XRHODJ=>NULL()  ! ( rhod J ) = dry density 
              ! for reference state * Jacobian of the GCS transformation.
  REAL                          ::  XREFMASS     ! Total mass of the ref.
                                                    ! state
  REAL                          ::  XMASS_O_PHI0 ! Mass / Phi0
  REAL                          ::  XLINMASS     ! Lineic mass for open
                                                    ! boundaries
! 
END TYPE REF_t

TYPE(REF_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: REF_MODEL

REAL, DIMENSION(:,:,:), POINTER :: XRHODREF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XTHVREF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRVREF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XEXNREF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRHODJ=>NULL()
REAL, POINTER :: XREFMASS=>NULL()
REAL, POINTER :: XMASS_O_PHI0=>NULL()
REAL, POINTER :: XLINMASS=>NULL()

CONTAINS

SUBROUTINE REF_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!REF_MODEL(KFROM)%XRHODREF=>XRHODREF !Done in FIELDLIST_GOTO_MODEL
!REF_MODEL(KFROM)%XTHVREF=>XTHVREF !Done in FIELDLIST_GOTO_MODEL
REF_MODEL(KFROM)%XRVREF=>XRVREF
REF_MODEL(KFROM)%XEXNREF=>XEXNREF
REF_MODEL(KFROM)%XRHODJ=>XRHODJ
!
! Current model is set to model KTO
!XRHODREF=>REF_MODEL(KTO)%XRHODREF !Done in FIELDLIST_GOTO_MODEL
!XTHVREF=>REF_MODEL(KTO)%XTHVREF !Done in FIELDLIST_GOTO_MODEL
XRVREF=>REF_MODEL(KTO)%XRVREF
XEXNREF=>REF_MODEL(KTO)%XEXNREF
XRHODJ=>REF_MODEL(KTO)%XRHODJ
XREFMASS=>REF_MODEL(KTO)%XREFMASS
XMASS_O_PHI0=>REF_MODEL(KTO)%XMASS_O_PHI0
XLINMASS=>REF_MODEL(KTO)%XLINMASS

END SUBROUTINE REF_GOTO_MODEL

END MODULE MODD_REF_n
