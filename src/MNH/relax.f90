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
!     ##################
      MODULE MODI_RELAX
!     ##################
!
INTERFACE
!
!
FUNCTION RELAX(PA,KB)  RESULT(PRELAX)
REAL, INTENT(IN)           :: PA   ! normalised distance to the inner boundary
INTEGER, INTENT(IN)        :: KB   ! maximum number of points in the rim zone
REAL                       :: PRELAX   ! result
END FUNCTION RELAX
!
END INTERFACE
!
END MODULE MODI_RELAX
!
!
!     #####################################
      FUNCTION RELAX(PA,KB)  RESULT(PRELAX)
!     #####################################
!
!!****  *RELAX* -  to compute lateral relaxation coefficient
!!
!!    PURPOSE
!!    -------
!!       The purpose of this function  is to compute the lateral relaxation
!!      coefficients which depend of the distance to the closest boundary.
!!      Distance greater than 1 corresponds to the corner areas. A
!!      threshold of 2. is applied in these areas.
!!
!!**  METHOD
!!    ------ 
!!      The profile variation of coefficients is specified by the user in this
!!    function.
!!      Note that the result is non-dimensioned.
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (RELAX function)
!!      Test of lateral boundary relaxation scheme in a barotropic model:
!!        Kallberg, 1977, ECMWF
!!      Le modele de prevision numerique Peridot, note de travail EERM n161
!!        juillet 1986, pg 30-31
!!
!!
!!    AUTHOR
!!    ------
!!      I. Mallet       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/03/96
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
REAL, INTENT(IN)           :: PA   ! normalised distance to the inner boundary
INTEGER, INTENT(IN)        :: KB   ! maximum number of points in the rim zone
REAL                       :: PRELAX   ! result
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL :: ZNBR            ! number of mesh-sizes to the lateral boundary
REAL :: ZFCT            ! Peridot coefficient
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION
!              ----------
!
!  previous formulation (sinus square)
!  PRELAX = SIN(PA*0.5*XPI)**2
!
!          Peridot profile
!
ZFCT = 0.45339
ZNBR = FLOAT(KB)*(1.-PA)
PRELAX = MIN(2.,ZFCT**ZNBR) 
!
!-------------------------------------------------------------------------------
!
END FUNCTION RELAX
