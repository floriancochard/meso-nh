!MNH_LIC Copyright 2018-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!  Creation:
!    P. Wautelet : 14/12/2018
!-----------------------------------------------------------------
module mode_io_tools_lfi

use modd_io_ll, only: tfiledata

implicit none

private

public :: io_prepare_verbosity_lfi

contains

subroutine io_prepare_verbosity_lfi(tpfile, kmelev, ostats)
  type(tfiledata),       intent(in)  :: tpfile
  integer(kind=LFI_INT), intent(out) :: kmelev
  logical,               intent(out) :: ostats

  select case (tpfile%nlfiverb)
    case(:2)
      ostats = .false.
      kmelev = 0
    case(3:6)
      ostats = .false.
      kmelev = 1
    case(7:9)
      ostats = .false.
      kmelev = 2
    case(10:)
      ostats = .true.
      kmelev = 2
  end select

end subroutine io_prepare_verbosity_lfi


end module mode_io_tools_lfi
