!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!     #########################
      MODULE MODD_PARAMETERS_ll
!     #########################
!
!!****  *MODD_PARAMETERS_ll* - declaration of parameter variables
!!                             communication layer
!
!!
!!    Purpose
!!    -------
!       The purpose of this declarative module is to specify  the variables 
!     which have the PARAMETER attribute   
!
!!   Reference
!!    ---------
!
!!    Authors
!!    -------
!
!     R. Guivarch               * CERFACS - ENSEEIHT *
!     Ph. Kloos                 * CERFACS - CNRM *
!     N. Gicquel                * CNRM *
!
!!    Implicit Arguments
!!    ------------------
!
!     None
!
!!    Modifications
!!    -------------
!
!    Original 04/05/99

!-------------------------------------------------------------------------------
!
!
  INTEGER :: JPHEXT        ! Horizontal External points number
  INTEGER :: JPVEXT        ! Vertical External points number
  INTEGER :: JPMODELMAX    ! Maximum allowed number of nested models 
!
  INTEGER, PARAMETER :: NMAXRIM = 10 ! maximum number of different RIM sizes
!
END MODULE MODD_PARAMETERS_ll
