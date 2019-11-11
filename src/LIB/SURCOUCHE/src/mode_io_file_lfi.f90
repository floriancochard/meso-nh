!MNH_LIC Copyright 2018-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!  Author: P. Wautelet 14/12/2018
!
!  Remarks: some of the code comes from mode_fm.f90 and mode_io.f90
!           (was duplicated in the 2 files)
!
!  Modifications:
!     Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!                                    + move IOFREEFLU and IONEWFLU to mode_io_file_lfi.f90
!                                    + move management of NNCID and NLFIFLU to the nc4 and lfi subroutines
!
!-----------------------------------------------------------------
module mode_io_file_lfi

use modd_io_ll,  only: tfiledata
use modd_netcdf, only: idcdf_kind

use mode_msg

implicit none

private

public :: io_create_file_lfi, io_close_file_lfi, io_open_file_lfi

integer, parameter :: JPRESERVED_UNIT   = 11
integer, parameter :: JPMAX_UNIT_NUMBER = JPRESERVED_UNIT + 300

logical,save :: galloc(JPRESERVED_UNIT:JPMAX_UNIT_NUMBER) = .false.

contains

subroutine io_create_file_lfi(tpfile, kstatus)
  use mode_io_tools,            only: io_construct_filename
  use mode_io_tools_lfi,        only: io_prepare_verbosity_lfi
  use mode_io_tools_mnhversion, only: io_set_mnhversion

  type(tfiledata), intent(inout) :: tpfile
  integer,         intent(inout) :: kstatus

  character(len=:), allocatable :: yfilem        ! name of the file
  character(len=:), allocatable :: yforstatus    ! Status for open of a file (for LFI) ('OLD','NEW','UNKNOWN','SCRATCH','REPLACE')
  integer(kind=LFI_INT)         :: iresou, inumbr
  integer(kind=LFI_INT)         :: imelev, inprar
  integer(kind=LFI_INT)         :: ininar        ! Number of articles present in LFI file
  logical                       :: gnewfi
  logical                       :: gnamfi, gfater, gstats

  call print_msg(NVERB_DEBUG,'IO','io_create_file_lfi','called for '//trim(tpfile%cname))

  kstatus = 0

  if (tpfile%lmaster) then
    call io_construct_filename(tpfile, yfilem)

    iresou = 0
    if ( tpfile%nlfiflu /= -1 ) call print_msg(NVERB_ERROR,'IO', &
                                               'io_create_file_lfi','file '//trim(yfilem)//'.lfi has already a unit number')
    tpfile%nlfiflu = ionewflu()
    gnamfi = .true.
    yforstatus = 'REPLACE'
    gfater = .true.

    call io_prepare_verbosity_lfi(tpfile, imelev, gstats)

    inumbr = tpfile%nlfiflu
    inprar = tpfile%nlfinprar
    call lfiouv(iresou, inumbr, gnamfi, trim(yfilem)//'.lfi', yforstatus, gfater, gstats, imelev, inprar, ininar)

    tpfile%nlfininar = ininar

    if (iresou/=0) kstatus = int(iresou, kind=kind(kstatus))

    !test if file is newly defined
    gnewfi = (ininar==0) .or. (imelev<2)
    if (.not.gnewfi) then
      call print_msg(NVERB_INFO,'IO','io_create_file_lfi','file '//trim(yfilem)//'.lfi previously created with LFI')
    endif
  end if
  call io_set_mnhversion(tpfile)
end subroutine io_create_file_lfi


subroutine io_close_file_lfi(tpfile, kstatus)
  type(tfiledata),   intent(inout)  :: tpfile
  integer, optional, intent(out)    :: kstatus

  character(len=*), parameter :: YSTATUS = 'KEEP'

  integer(kind=LFI_INT) :: istatus

  call print_msg(NVERB_DEBUG,'IO','io_close_file_lfi','called for '//trim(tpfile%cname))

  istatus = 0

  if (tpfile%lmaster) then
    if ( tpfile%nlfiflu /= -1 ) then
      call lfifer(istatus, tpfile%nlfiflu, YSTATUS)
      call iofreeflu(int(tpfile%nlfiflu))
      tpfile%nlfiflu = -1
    else
      istatus = -1
      call print_msg(NVERB_WARNING, 'IO', 'io_close_file_lfi', 'file '//trim(tpfile%cname)//'.lfi is not opened')
    end if
  end if

  if (present(kstatus)) kstatus = int(istatus,kind=kind(kstatus))
end subroutine io_close_file_lfi


subroutine io_open_file_lfi(tpfile, kstatus)
  use mode_io_tools,            only: io_construct_filename
  use mode_io_tools_lfi,        only: io_prepare_verbosity_lfi
  use mode_io_tools_mnhversion, only: io_get_mnhversion

  type(tfiledata), intent(inout) :: tpfile
  integer,         intent(inout) :: kstatus

  character(len=:),allocatable :: yfilem        ! name of the file
  character(len=:),allocatable :: yforstatus    ! Status for open of a file (for LFI) ('OLD','NEW','UNKNOWN','SCRATCH','REPLACE')
  integer                      :: istatus
  integer(kind=LFI_INT)        :: iresou, inumbr
  integer(kind=LFI_INT)        :: imelev, inprar
  integer(kind=LFI_INT)        :: ininar        ! Number of articles present in LFI file
  logical                      :: gnewfi
  logical                      :: gnamfi, gfater, gstats

  call print_msg(NVERB_DEBUG,'IO','io_open_file_lfi','called for '//trim(tpfile%cname))

  kstatus = 0

  if (tpfile%lmaster) then
    call io_construct_filename(tpfile, yfilem)

    iresou = 0
    if ( tpfile%nlfiflu /= -1 ) call print_msg(NVERB_ERROR,'IO', &
                                               'io_open_file_lfi','file '//trim(yfilem)//'.lfi has already a unit number')
    tpfile%nlfiflu = ionewflu()
    gnamfi = .true.
    yforstatus = 'OLD'
    gfater = .true.

    call io_prepare_verbosity_lfi(tpfile, imelev, gstats)

    inumbr = tpfile%nlfiflu
    inprar = tpfile%nlfinprar

    call lfiouv(iresou, inumbr, gnamfi, trim(yfilem)//'.lfi', yforstatus, gfater, gstats, imelev, inprar, ininar)

    tpfile%nlfininar = ininar

    if (iresou/=0) kstatus = int(iresou, kind=kind(kstatus))
  end if
  call io_get_mnhversion(tpfile)
end subroutine io_open_file_lfi


function ionewflu()
  use modd_io_ll, only: nnullunit

  integer :: ionewflu

  integer :: ji
  integer :: ios
  logical :: gexists, gopened, gfound

  gfound = .false.

  do ji = JPRESERVED_UNIT, JPMAX_UNIT_NUMBER
    if ( galloc(ji) ) cycle
    inquire(unit=ji, exist=gexists, opened=gopened, iostat=ios)
    if (gexists .and. .not. gopened .and. ios == 0) then
      ionewflu   = ji
      gfound     = .true.
      galloc(ji) = .true.
      exit
    end if
  end do

  if (.not. gfound) then
    call print_msg(NVERB_ERROR,'IO','ionewflu','wrong unit number')
    ionewflu = nnullunit !/dev/null Fortran unit
  end if
end function ionewflu


subroutine iofreeflu(koflu)
  integer :: koflu

  if ( (koflu >= JPRESERVED_UNIT) .and. (koflu <= JPMAX_UNIT_NUMBER) ) then
    galloc(koflu) = .false.
  else
    call print_msg(NVERB_ERROR,'IO','iofreeflu','wrong unit number')
  end if
end subroutine iofreeflu


end module mode_io_file_lfi
