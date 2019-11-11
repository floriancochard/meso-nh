!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!  Modifications:
!    P. Wautelet : 13/12/2018 : extracted from mode_io.f90
!    P. Wautelet : 14/12/2018 : added io_construct_filename
!-----------------------------------------------------------------
module mode_io_tools

use modd_io_ll, only: tfiledata

implicit none

private

public :: io_file, io_rank, io_construct_filename

contains

  FUNCTION io_file(k,nb_proc_io)
    !
    ! return the file number where to write the K level of data
    !
    IMPLICIT NONE
    INTEGER(kind=MNH_MPI_RANK_KIND)                   :: k,nb_proc_io
    INTEGER(kind=MNH_MPI_RANK_KIND)                   :: io_file

    io_file = MOD ((k-1) , nb_proc_io )

  END FUNCTION io_file

  FUNCTION IO_RANK(IFILE,nb_proc,nb_proc_io,offset_rank)
    !
    ! return the proc number which must write the 'IFILE' file
    !
    IMPLICIT NONE
    INTEGER(kind=MNH_MPI_RANK_KIND)                  :: IFILE,nb_proc,nb_proc_io
    INTEGER(kind=MNH_MPI_RANK_KIND),OPTIONAL         :: offset_rank

    INTEGER(kind=MNH_MPI_RANK_KIND)                  :: IO_RANK

    INTEGER(kind=MNH_MPI_RANK_KIND)                  :: ipas,irest

    ipas  =        nb_proc / nb_proc_io
    irest =  MOD ( nb_proc , nb_proc_io )

    IF  (ipas /= 0 ) THEN
       IO_RANK=ipas * IFILE + MIN(IFILE , irest )
    ELSE
       IO_RANK=MOD(IFILE , nb_proc )
    ENDIF

    !
    ! optional rank to shift for read test
    !
    IF (PRESENT(offset_rank)) THEN
       IF ( offset_rank .GT.0 ) IO_RANK=MOD(IO_RANK+offset_rank,nb_proc)
       IF ( offset_rank .LT.0 ) IO_RANK=MOD(nb_proc-IO_RANK+offset_rank,nb_proc)
    ENDIF

  END FUNCTION IO_RANK


subroutine io_construct_filename(tpfile,hfilem)
  type(tfiledata),               intent(inout) :: tpfile
  character(len=:), allocatable, intent(out)   :: hfilem

  if (allocated(tpfile%cdirname)) then
    if(len_trim(tpfile%cdirname)>0) then
      hfilem = trim(tpfile%cdirname)//'/'//trim(tpfile%cname)
    else
      hfilem = trim(tpfile%cname)
    end if
  else
    hfilem = trim(tpfile%cname)
  end if

end subroutine io_construct_filename


end module mode_io_tools



module mode_io_tools_mnhversion

use modd_io_ll, only: tfiledata

use mode_msg

implicit none

private

public :: io_get_mnhversion, io_set_mnhversion

contains

  subroutine io_get_mnhversion(tpfile)
  !Compare MNHVERSION of file with current version and store it in file metadata
    use modd_conf,   only: nmnhversion
    use modd_io_ll,  only: tfiledata
    use mode_field,  only: tfielddata,typeint
    use mode_fmread, only: io_read_field

    type(tfiledata), intent(inout) :: tpfile

    character(len=12)       :: ymnhversion_file,ymnhversion_curr
    integer :: imasdev,ibugfix
    integer :: iresp
    integer,dimension(3)    :: imnhversion
    type(tfielddata)        :: tzfield

    call print_msg(NVERB_DEBUG,'IO','io_get_mnhversion','called for '//trim(tpfile%cname))

    if ( trim(tpfile%cmode) /= 'READ' ) &
      call print_msg(NVERB_FATAL,'IO','io_get_mnhversion',trim(tpfile%cname)// 'not opened in read mode')

    imnhversion(:) = 0
    !use tzfield because tfieldlist could be not initialised
    tzfield%cmnhname   = 'MNHVERSION'
    tzfield%cstdname   = ''
    tzfield%clongname  = 'MesoNH version'
    tzfield%cunits     = ''
    tzfield%cdir       = '--'
    tzfield%ccomment   = ''
    tzfield%ngrid      = 0
    tzfield%ntype      = TYPEINT
    tzfield%ndims      = 1
    tzfield%ltimedep   = .false.
    call io_read_field(tpfile,tzfield,imnhversion,iresp)
    if (iresp/=0) then
      tzfield%cmnhname   = 'MASDEV'
      tzfield%clongname  = 'MesoNH version (without bugfix)'
      tzfield%ndims      = 0
      call io_read_field(tpfile,tzfield,imasdev,iresp)
      if (iresp/=0) then
        call print_msg(NVERB_WARNING,'IO','io_get_mnhversion','unknown MASDEV version for '//trim(tpfile%cname))
      else
        if (imasdev<100) then
          imnhversion(1)=imasdev/10
          imnhversion(2)=mod(imasdev,10)
        else !for example for mnh 4.10
          imnhversion(1)=imasdev/100
          imnhversion(2)=mod(imasdev,100)
        end if
      end if
      !
      tzfield%cmnhname   = 'BUGFIX'
      tzfield%clongname  = 'MesoNH bugfix number'
      call io_read_field(tpfile,tzfield,ibugfix,iresp)
      if (iresp/=0) then
        call print_msg(NVERB_WARNING,'IO','io_get_mnhversion','unknown BUGFIX version for '//trim(tpfile%cname))
      else
        imnhversion(3)=ibugfix
      end if
    end if
    !
    write(ymnhversion_file,"( I0,'.',I0,'.',I0 )" ) imnhversion(1),imnhversion(2),imnhversion(3)
    write(ymnhversion_curr,"( I0,'.',I0,'.',I0 )" ) nmnhversion(1),nmnhversion(2),nmnhversion(3)
    !
    if ( imnhversion(1)==0 .and. imnhversion(2)==0 .and. imnhversion(3)==0 ) then
      call print_msg(NVERB_WARNING,'IO','io_get_mnhversion','file '//trim(tpfile%cname)//&
                    ' was written with an unknown version of MesoNH')
      else if (  imnhversion(1)< nmnhversion(1) .or. &
              (imnhversion(1)==nmnhversion(1) .and. imnhversion(2)< nmnhversion(2)) .or. &
              (imnhversion(1)==nmnhversion(1) .and. imnhversion(2)==nmnhversion(2) .and. imnhversion(3)<nmnhversion(3)) ) then
      call print_msg(NVERB_WARNING,'IO','io_get_mnhversion','file '//trim(tpfile%cname)//&
                      ' was written with an older version of MesoNH ('//trim(ymnhversion_file)//&
                      ' instead of '//trim(ymnhversion_curr)//')')
      else if (  imnhversion(1)> nmnhversion(1) .or. &
                (imnhversion(1)==nmnhversion(1) .and. imnhversion(2)> nmnhversion(2)) .or. &
                (imnhversion(1)==nmnhversion(1) .and. imnhversion(2)==nmnhversion(2) .and. imnhversion(3)>nmnhversion(3)) ) then
        call print_msg(NVERB_WARNING,'IO','io_get_mnhversion','file '//trim(tpfile%cname)//&
                      ' was written with a more recent version of MesoNH ('//trim(ymnhversion_file)//&
                      ' instead of '//trim(ymnhversion_curr)//')')
      else
        call print_msg(NVERB_DEBUG,'IO','io_get_mnhversion','file '//trim(tpfile%cname)//&
                      ' was written with the same version of MesoNH ('//trim(ymnhversion_curr)//')')
      end if
      !
      tpfile%nmnhversion(:) = imnhversion(:)
  end subroutine io_get_mnhversion


  subroutine io_set_mnhversion(tpfile)
    use modd_conf,  only: nmnhversion
    use modd_io_ll, only: tfiledata

    type(tfiledata), intent(inout) :: tpfile

    call print_msg(NVERB_DEBUG,'IO','io_set_mnhversion','called for '//trim(tpfile%cname))

    if ( trim(tpfile%cmode) /= 'WRITE' ) &
      call print_msg(NVERB_FATAL,'IO','io_set_mnhversion',trim(tpfile%cname)// 'not opened in write mode')

    tpfile%nmnhversion(:) = nmnhversion(:)
  end subroutine io_set_mnhversion

end module mode_io_tools_mnhversion
