!MNH_LIC Copyright 2001-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_WRITE_BALLOON_n
!     ###########################
!
INTERFACE
!
SUBROUTINE WRITE_BALLOON_n(TPFILE)
USE MODD_IO_ll, ONLY: TFILEDATA
!
IMPLICIT NONE
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! File characteristics
!
END SUBROUTINE WRITE_BALLOON_n
!
END INTERFACE
!
END MODULE MODI_WRITE_BALLOON_n
!
!
!     ###################################
      SUBROUTINE WRITE_BALLOON_n(TPFILE)
!     ###################################
!
!!****  *WRITE_BALLOON_n* - routine to write balloon records in a LFIFM file
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      FMWRIT     : FM-routine to write a record
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_AIRCRAFT_BALLOON_n : contains balloon and aircraft variables
!!      Module MODD_GRID_n : contains spatial grid variables
!!      Module MODD_LUNIT_n   : contains logical unit variables
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!  	G.Jaubert   *Meteo France* 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/06/01 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_AIRCRAFT_BALLOON
USE MODD_GRID,  ONLY: XLONORI,XLATORI
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LUNIT_n
!
USE MODE_GRIDPROJ
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! File characteristics
!
!*       0.2   Declarations of local variables
!
!
IF (TBALLOON1%FLY) CALL WRITE_LFI_BALLOON(TBALLOON1)
IF (TBALLOON2%FLY) CALL WRITE_LFI_BALLOON(TBALLOON2)
IF (TBALLOON3%FLY) CALL WRITE_LFI_BALLOON(TBALLOON3)
IF (TBALLOON4%FLY) CALL WRITE_LFI_BALLOON(TBALLOON4)
IF (TBALLOON5%FLY) CALL WRITE_LFI_BALLOON(TBALLOON5)
IF (TBALLOON6%FLY) CALL WRITE_LFI_BALLOON(TBALLOON6)
IF (TBALLOON7%FLY) CALL WRITE_LFI_BALLOON(TBALLOON7)
IF (TBALLOON8%FLY) CALL WRITE_LFI_BALLOON(TBALLOON8)
IF (TBALLOON9%FLY) CALL WRITE_LFI_BALLOON(TBALLOON9)
!
!
CONTAINS
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE WRITE_LFI_BALLOON(TPFLYER)
!
USE MODE_FIELD, ONLY : TFIELDDATA, TYPEREAL
USE MODE_FMWRIT
!
TYPE(FLYER),        INTENT(IN)       :: TPFLYER
!
!
!*       0.2   Declarations of local variables
!
REAL               :: ZLAT          ! latitude of the balloon
REAL               :: ZLON          ! longitude of the balloon
TYPE(TFIELDDATA)   :: TZFIELD
!
!
CALL SM_LATLON(XLATORI,XLONORI,  &
     TPFLYER%X_CUR,TPFLYER%Y_CUR,ZLAT,ZLON)
!
!
TZFIELD%CMNHNAME   = TRIM(TPFLYER%TITLE)//'LAT'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
TZFIELD%CUNITS     = 'degree'
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = ''
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZLAT)
!
TZFIELD%CMNHNAME   = TRIM(TPFLYER%TITLE)//'LON'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
TZFIELD%CUNITS     = 'degree'
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = ''
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZLON)
!
TZFIELD%CMNHNAME   = TRIM(TPFLYER%TITLE)//'ALT'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
TZFIELD%CUNITS     = 'm'
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = ''
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPFILE,TZFIELD,TPFLYER%Z_CUR)
!
TZFIELD%CMNHNAME   = TRIM(TPFLYER%TITLE)//'WASCENT'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
TZFIELD%CUNITS     = 'm s-1'
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = ''
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPFILE,TZFIELD,TPFLYER%WASCENT)
!
TZFIELD%CMNHNAME   = TRIM(TPFLYER%TITLE)//'RHO'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
TZFIELD%CUNITS     = 'kg m-3'
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = ''
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .TRUE.
CALL IO_WRITE_FIELD(TPFILE,TZFIELD,TPFLYER%RHO)
!
!
!
END SUBROUTINE WRITE_LFI_BALLOON
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE WRITE_BALLOON_n
