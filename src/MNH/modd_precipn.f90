!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_PRECIP_n
!     ####################
!
!!****  *MODD_PRECIP$n* - declaration of precipitating fields
!!
!!    PURPOSE
!!    -------
!       Stores the INstantaneous and ACcumulated PRecipitating fields of 
!!      resolved clouds
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_PRECIPn)
!!          
!!    AUTHOR
!!    ------
!!      J.-P. Pinty *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!       J.-P. Pinty   29/11/02 add C3R5, ICE2, ICE4
!!                    10/2016 (C.Lac) Add droplet deposition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

REAL, DIMENSION(:,:), POINTER :: XINPRC=>NULL(), XACPRC=>NULL(), & ! Instant and cumul of ground
                                 XINDEP=>NULL(), XACDEP=>NULL(), & ! precipitation fields of Cloud, Rain,
                                 XINPRR=>NULL(), XACPRR=>NULL(), & ! Snow, Graupel and Hail
                                 XINPRS=>NULL(), XACPRS=>NULL(), &
                                 XINPRG=>NULL(), XACPRG=>NULL(), &
                                 XINPRH=>NULL(), XACPRH=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XINPRR3D=>NULL(), & !Instant precip 3D
                                   XEVAP3D =>NULL()    ! Evap profile 3D

CONTAINS

SUBROUTINE PRECIP_GOTO_MODEL(KFROM, KTO)
!
INTEGER, INTENT(IN) :: KFROM, KTO
!
END SUBROUTINE PRECIP_GOTO_MODEL

END MODULE MODD_PRECIP_n
