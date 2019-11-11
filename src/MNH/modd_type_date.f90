!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
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
