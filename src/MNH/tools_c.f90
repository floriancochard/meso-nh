!MNH_LIC Copyright 2018-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------

!#############
module modi_tools_c
!#############
!
!!    Purpose
!!    -------
!
!     The Purpose of this module is to provide interfaces
!     to C functions
!
!    Authors
!    -------
!
!     P. Wautelet 04/12/2018
!
  use, intrinsic :: iso_c_binding

  implicit none

  interface
    subroutine sleep_c(ksec) bind(c, name="sleep")
      import C_INT
      integer(kind=C_INT), VALUE :: ksec
    end subroutine sleep_c
  end interface

end module modi_tools_c
