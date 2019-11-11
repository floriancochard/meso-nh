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

!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/d10/free/mesonh/sources/modd/s.modd_type_date.f90, Version:1.6, Date:01/01/15, Last modified:00/12/21
!-----------------------------------------------------------------
!     #################
      MODULE MODD_TYPE_DATE
!     #################
!
!!****  *MODD_TYPE_DATE* - declaration of temporal types
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to define
!      the time types. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE 
!!
!!    REFERENCE
!!    --------- 
!!      Book2 of documentation of Meso-NH (module MODD_TYPE_DATE)
!!       
!!    AUTHOR
!!    ------
!!	P. Jabouille   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/08/97                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
TYPE DATE
INTEGER :: YEAR
INTEGER :: MONTH
INTEGER :: DAY
END TYPE DATE
!
TYPE DATE_TIME
TYPE (DATE) :: TDATE
REAL :: TIME
END TYPE DATE_TIME 
!
END MODULE MODD_TYPE_DATE
