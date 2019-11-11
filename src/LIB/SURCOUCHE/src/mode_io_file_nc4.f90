!MNH_LIC Copyright 2018-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!  Author: P. Wautelet 13/12/2018
!
!  Remarks: some of the code comes from mode_fm.f90 and mode_io.f90
!           (was duplicated in the 2 files)
!
!  Modifications:
!     Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!                                    + move IOFREEFLU and IONEWFLU to mode_io_file_lfi.f90
!                                    + move management of NNCID and NLFIFLU to the nc4 and lfi subroutines
!     Philippe Wautelet: 10/01/2019: replace handle_err by io_handle_err_nc4 for better netCDF error messages
!  P. Wautelet 07/03/2019: bugfix: io_set_mnhversion must be called by all the processes
!
!-----------------------------------------------------------------
#if defined(MNH_IOCDF4)
module mode_io_file_nc4

use modd_io_ll,  only: tfiledata
use modd_netcdf, only: IDCDF_KIND

use mode_io_tools_nc4, only: io_handle_err_nc4, io_set_knowndims_nc4, newiocdf
use mode_msg

use NETCDF,      only: NF90_CLOBBER, NF90_GLOBAL, NF90_NETCDF4, NF90_NOERR, NF90_NOWRITE,  &
                       NF90_CLOSE, NF90_CREATE, NF90_GET_ATT, NF90_INQUIRE, NF90_INQUIRE_ATTRIBUTE, &
                       NF90_OPEN, NF90_PUT_ATT, NF90_STRERROR

implicit none

private

public :: io_create_file_nc4, io_close_file_nc4, io_open_file_nc4

contains

subroutine io_create_file_nc4(tpfile,hprogram_orig)
  use mode_io_tools,            only: io_construct_filename
  use mode_io_tools_mnhversion, only: io_set_mnhversion

  type(tfiledata),           intent(inout) :: tpfile
  character(len=*),optional, intent(in)    :: hprogram_orig !to emulate a file coming from this program

  character(len=:),allocatable :: yfilem  ! name of the file
  integer(kind=IDCDF_KIND)     :: istatus

  call print_msg(NVERB_DEBUG,'IO','io_create_file_nc4','called for '//trim(tpfile%cname))

  if (tpfile%lmaster) then
    call io_construct_filename(tpfile, yfilem)

    tpfile%tncdims => newiocdf()
    istatus = NF90_CREATE(adjustl(trim(yfilem))//".nc", ior(NF90_CLOBBER,NF90_NETCDF4), tpfile%nncid)
    if (istatus /= NF90_NOERR) then
      call print_msg(NVERB_FATAL,'IO','io_create_file_nc4','NF90_CREATE for '//trim(yfilem)//'.nc: '//NF90_STRERROR(istatus))
    end if
    call io_set_not_cleanly_closed_nc4(tpfile)
    call io_set_knowndims_nc4(tpfile, hprogram_orig=hprogram_orig)
  end if
  call io_set_mnhversion(tpfile)
end subroutine io_create_file_nc4


subroutine io_close_file_nc4(tpfile,kstatus)
  use mode_io_tools_nc4, only: cleaniocdf

  type(tfiledata),                    intent(inout) :: tpfile
  integer(kind=IDCDF_KIND), optional, intent(out)   :: kstatus

  integer(kind=IDCDF_KIND) :: istatus

  call print_msg(NVERB_DEBUG,'IO','io_close_file_nc4','called for '//trim(tpfile%cname))

  istatus = 0

  if (tpfile%lmaster ) then
    if (tpfile%nncid == -1) then
      call print_msg(NVERB_WARNING, 'IO', 'io_close_file_nc4', 'file '//trim(tpfile%cname)//'.nc is not opened')
    else
      if (trim(tpfile%cmode) == 'WRITE') call io_set_cleanly_closed_nc4(tpfile)
      istatus = NF90_CLOSE(tpfile%nncid)
      if (istatus /= NF90_NOERR) then
        call print_msg(NVERB_WARNING, 'IO', 'io_close_file_nc4', 'NF90_CLOSE error: '//trim(NF90_STRERROR(istatus)))
      end if
      tpfile%nncid = -1
      if (associated(tpfile%tncdims)) call cleaniocdf(tpfile%tncdims)
    end if
  end if

  if (present(kstatus)) kstatus = istatus
end subroutine io_close_file_nc4


subroutine io_open_file_nc4(tpfile)
  use mode_io_tools,            only: io_construct_filename
  use mode_io_tools_mnhversion, only: io_get_mnhversion

  type(tfiledata), intent(inout) :: tpfile

  character(len=:),allocatable :: yfilem  ! name of the file
  integer(kind=IDCDF_KIND)     :: istatus

  call print_msg(NVERB_DEBUG,'IO','io_open_file_nc4','called for '//trim(tpfile%cname))

  if (tpfile%lmaster) then
    call io_construct_filename(tpfile, yfilem)

    tpfile%tncdims => newiocdf()
    istatus = NF90_OPEN(adjustl(trim(yfilem))//".nc", NF90_NOWRITE, tpfile%nncid)
    if (istatus /= NF90_NOERR) then
      call print_msg(NVERB_FATAL, 'IO', 'io_open_file_nc4', 'NF90_OPEN for '//trim(yfilem)//'.nc: '//NF90_STRERROR(istatus))
    end if

    istatus = NF90_INQUIRE(tpfile%nncid, nvariables=tpfile%nncnar)
    if (istatus /= NF90_NOERR) then
      call print_msg(NVERB_FATAL,'IO','io_open_file_nc4','NF90_INQUIRE for '//trim(yfilem)//'.nc: '//NF90_STRERROR(istatus))
    end if
  end if

  if (trim(tpfile%cmode) == 'READ') then
    call io_get_mnhversion(tpfile)
    if (tpfile%lmaster) call io_check_cleanly_closed_nc4(tpfile)
  end if

end subroutine io_open_file_nc4


subroutine io_check_cleanly_closed_nc4(tpfile)
  type(tfiledata), intent(in) :: tpfile

  character(len=:), allocatable :: yclean
  integer(kind=IDCDF_KIND) :: ilen, istatus
  integer, dimension(3) :: imnhversion

  call print_msg(NVERB_DEBUG,'IO','io_check_cleanly_closed_nc4','called for '//trim(tpfile%cname))

  imnhversion = tpfile%nmnhversion
  if ( imnhversion(1)<5                                                 .OR. &
      (imnhversion(1)==5 .AND. imnhversion(2)<4)                        .OR. &
      (imnhversion(1)==5 .AND. imnhversion(2)==4 .AND. imnhversion(3)<2)     ) then
    call print_msg(NVERB_DEBUG,'IO','io_check_cleanly_closed_nc4', &
                   'file '//trim(tpfile%cname)//' is too old (before MNH 5.4.2) to check if cleanly closed')
    return
  end if

  istatus = NF90_INQUIRE_ATTRIBUTE(tpfile%nncid, NF90_GLOBAL, 'MNH_cleanly_closed', len = ilen)
  if (istatus /= NF90_NOERR) then
    call print_msg(NVERB_ERROR,'IO','io_check_cleanly_closed_nc4', &
                   'MNH_cleanly_closed attribute not found in file '//trim(tpfile%cname))
  else
    allocate( character(len=ilen) :: yclean )
    istatus = NF90_GET_ATT(tpfile%nncid, NF90_GLOBAL, 'MNH_cleanly_closed', yclean)
    if (istatus /= NF90_NOERR) then
      call print_msg(NVERB_WARNING,'IO','io_check_cleanly_closed_nc4', &
                    'MNH_cleanly_closed attribute not found in file '//trim(tpfile%cname))
    else
      if (yclean == 'yes') then
        call print_msg(NVERB_DEBUG,'IO','io_check_cleanly_closed_nc4', &
                      'file '//trim(tpfile%cname)//' was cleanly closed before opening')
      else if (yclean == 'no') then
        call print_msg(NVERB_ERROR,'IO','io_check_cleanly_closed_nc4', &
                      'file '//trim(tpfile%cname)//' was not cleanly closed before opening')
      else
        call print_msg(NVERB_ERROR,'IO','io_check_cleanly_closed_nc4', &
                      'invalid MNH_cleanly_closed attribute for file '//trim(tpfile%cname))
      end if
    end if
  end if
end subroutine io_check_cleanly_closed_nc4


subroutine io_set_cleanly_closed_nc4(tpfile)
  type(tfiledata), intent(in) :: tpfile

  integer(kind=IDCDF_KIND) :: istatus

  call print_msg(NVERB_DEBUG,'IO','io_set_cleanly_closed_nc4','called for '//trim(tpfile%cname))

  istatus = NF90_PUT_ATT(tpfile%nncid, NF90_GLOBAL, 'MNH_cleanly_closed', 'yes')
  if (istatus /= NF90_NOERR) call io_handle_err_nc4(istatus,'io_set_cleanly_closed_nc4','NF90_PUT_ATT','MNH_cleanly_closed')
end subroutine io_set_cleanly_closed_nc4


subroutine io_set_not_cleanly_closed_nc4(tpfile)
  type(tfiledata), intent(in) :: tpfile

  integer(kind=IDCDF_KIND) :: istatus

  call print_msg(NVERB_DEBUG,'IO','io_set_not_cleanly_closed_nc4','called for '//trim(tpfile%cname))

  istatus = NF90_PUT_ATT(tpfile%nncid, NF90_GLOBAL, 'MNH_cleanly_closed', 'no')
  if (istatus /= NF90_NOERR) call io_handle_err_nc4(istatus,'io_set_not_cleanly_closed_nc4','NF90_PUT_ATT','MNH_cleanly_closed')
end subroutine io_set_not_cleanly_closed_nc4

end module mode_io_file_nc4
#else
!
! External dummy subroutines
!
subroutine io_create_file_nc4(a, b)
use mode_msg
integer :: a, b
CALL PRINT_MSG(NVERB_ERROR,'IO','io_create_file_nc4','empty call. Compile with -DMNH_IOCDF4 flag to enable NetCDF')
end subroutine io_create_file_nc4
!
subroutine io_close_file_nc4(a)
use mode_msg
integer :: a
CALL PRINT_MSG(NVERB_ERROR,'IO','io_close_file_nc4','empty call. Compile with -DMNH_IOCDF4 flag to enable NetCDF')
end subroutine io_close_file_nc4
!
subroutine io_open_file_nc4(a)
use mode_msg
integer :: a
CALL PRINT_MSG(NVERB_ERROR,'IO','io_open_file_nc4','empty call. Compile with -DMNH_IOCDF4 flag to enable NetCDF')
end subroutine io_open_file_nc4
!
#endif
