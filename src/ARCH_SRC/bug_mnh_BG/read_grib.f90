!#################################
MODULE MODI_READ_GRIB
!#################################
  INTERFACE
    SUBROUTINE READ_GRIB(HGRIB,KLUOUT,KPARAM,KRET,          &
                         KSEC1,KSEC2,PSEC2,KNI,PDOUT,       &
                         KLTYPE,KLEV1,KLEV2                 )
 CHARACTER(LEN=*),                   INTENT(IN)    :: HGRIB  ! Grib file name
 INTEGER,                            INTENT(IN)    :: KLUOUT ! logical unit of output listing
 INTEGER,                            INTENT(IN)    :: KPARAM ! Parameter to read
 INTEGER,                            INTENT(OUT)   :: KRET   ! Result
 INTEGER, OPTIONAL, DIMENSION(1024), INTENT(OUT)   :: KSEC1  ! These are parameters to pass
 INTEGER, OPTIONAL, DIMENSION(1024), INTENT(OUT)   :: KSEC2  ! to each call (used to avoid
 REAL,    OPTIONAL, DIMENSION(512),  INTENT(OUT)   :: PSEC2  ! making calculations twice)
 INTEGER, OPTIONAL,                  INTENT(IN)    :: KNI    ! Number of input points
 REAL,    OPTIONAL, DIMENSION(:),    POINTER       :: PDOUT  ! Output datas
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLTYPE ! type of level (Grib code table 3)
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLEV1  ! level definition
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLEV2  ! level definition
    END SUBROUTINE READ_GRIB
  END INTERFACE
END MODULE MODI_READ_GRIB
!
!     ##########################################################################
    SUBROUTINE READ_GRIB(HGRIB,KLUOUT,KPARAM,KRET,          &
                         KSEC1,KSEC2,PSEC2,KNI,PDOUT,       &
                         KLTYPE,KLEV1,KLEV2                 )
!     ##########################################################################
!!****  *READ_GRIB* - Read a grib field & performs interpolation (optional)
!!
!!    PURPOSE
!!    -------
!!
!!    Searchs & reads a field in a grib file. Returns the field or an interpolated
!!    one.
!!
!!    METHOD
!!    ------
!!
!!   The field to read is defined by :
!!    . its number (required),
!!    . its level type (may be set to -1 to accept any type),
!!    . its level id 1 (may be set to -1 to accept any value),
!!    . its level id 2 (may be set to -1 to accept any value).
!!   If '-1' values are used, the routine returns the founded values.
!!
!!
!!   EXTERNAL
!!   --------
!!   
!!   subroutine PBOPEN        : open a grib file
!!   subroutine PBGRIB        : read datas from a grib file
!!   subroutine PBCLOSE       : close a gribfile
!!   subroutine GRIBEX        : decodes grib data
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
!!   REFERENCE
!!   ---------
!!
!!   This routine is based on books describing the Grib file format :
!!     'Encoding and Decoding Grib data', John D.Chambers(ECMWF), October 1995
!!     'Accessing GRIB and BUFR data', John D.Chambers(ECMWF), May 1994
!!     'A guide to Grib' Edition 1, John D.Stackpole(NOAA), March 1994
!!
!!   AUTHOR
!!   ------
!!
!!   V.Bousquet
!!
!!   MODIFICATIONS
!!   -------------
!!
!!   Original       07/01/1999
!!
!!   Modification 06/2003 (V. Masson) simplification for externalization
!!   Modification 11/2005 (I. Mallet) increase JPACK to read CEP files 
!!                                    from cy29r1 (24bits coding)
!!   Modification 03/2006 (I. Mallet) increase JPACK to read CEP files 
!!-----------------------------------------------------------------------------------
!
! 0. DECLARATIONS
! ---------------
!
!
IMPLICIT NONE
!
! 0.1. Declaration of arguments
! -----------------------------
 CHARACTER(LEN=*),                   INTENT(IN)    :: HGRIB  ! Grib file name
 INTEGER,                            INTENT(IN)    :: KLUOUT ! logical unit of output listing
 INTEGER,                            INTENT(IN)    :: KPARAM ! Parameter to read
 INTEGER,                            INTENT(OUT)   :: KRET   ! Result
 INTEGER, OPTIONAL, DIMENSION(1024), INTENT(OUT)   :: KSEC1  ! These are parameters to pass
 INTEGER, OPTIONAL, DIMENSION(1024), INTENT(OUT)   :: KSEC2  ! to each call (used to avoid
 REAL,    OPTIONAL, DIMENSION(512),  INTENT(OUT)   :: PSEC2  ! making calculations twice)
 INTEGER, OPTIONAL,                  INTENT(IN)    :: KNI    ! Number of input points
 REAL,    OPTIONAL, DIMENSION(:),    POINTER       :: PDOUT  ! Output datas
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLTYPE ! type of level (Grib code table 3)
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLEV1  ! level definition
 INTEGER, OPTIONAL,                  INTENT(INOUT) :: KLEV2  ! level definition
!
!-----------------------------------------------------------------------------------
!
END SUBROUTINE READ_GRIB
