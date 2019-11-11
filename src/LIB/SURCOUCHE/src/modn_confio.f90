!MNH_LIC Copyright 2014-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODN_CONFIO
!     ##################
!
!!****  *MODN_CONFIO* - declaration of namelist NAM_CONFIO
!!
!!    PURPOSE
!!    -------
!     Define I/O configuration variables that can be set with the NAM_CONFIO namelist
!!    /!\ These variables must be transmitted to the SURCOUCHE library via the
!!    SET_CONFIO_ll subroutine before the FIRST call to IO_FILE_OPEN_ll
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      D.Gazen L.A.
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    31/03/2014
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  Philippe Wautelet: 21/01/2019: add LIO_ALLOW_NO_BACKUP and LIO_NO_WRITE to modd_io_ll to allow to disable writes (for bench purposes)
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_IO_ll, ONLY : LVERB_OUTLST, LVERB_STDOUT, LVERB_ALLPRC, &
                       NIO_VERB, NIO_ABORT_LEVEL, NGEN_VERB, NGEN_ABORT_LEVEL, &
                       CIO_DIR, LIO_ALLOW_NO_BACKUP, LIO_NO_WRITE
!
IMPLICIT NONE
!
LOGICAL,SAVE :: LCDF4    = .FALSE. ! TRUE : enable NetCDF4 Input/Output
LOGICAL,SAVE :: LLFIOUT  = .FALSE. ! TRUE : add LFI output when NetCDF4 I/O is enabled (debug)
LOGICAL,SAVE :: LLFIREAD = .FALSE. ! TRUE : enable LFI reading (disable NetCDF4 reading)
                                   !        when NetCDF4 I/O is enabled (debug)

NAMELIST/NAM_CONFIO/LCDF4, LLFIOUT, LLFIREAD,                 &
                    LVERB_OUTLST, LVERB_STDOUT, LVERB_ALLPRC, &
                    NIO_VERB,  NIO_ABORT_LEVEL,               &
                    NGEN_VERB, NGEN_ABORT_LEVEL, CIO_DIR,     &
                    LIO_ALLOW_NO_BACKUP, LIO_NO_WRITE
!
END MODULE MODN_CONFIO

