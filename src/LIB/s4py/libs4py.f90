!MNH_LIC Copyright 2014-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
SUBROUTINE WLFIOUV(KRETURNCODE, CDFILE, CDSTATE, KNUMER)
! ** PURPOSE
!    Open a LFI file
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    CDFILE: path to file to open
!    CDSTATE: state of file ('NEW', 'OLD', 'UNKNOWN', 'SCRATCH')
!    KNUMER: logical unit number associated to file
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!    8 nov 2018, S. Riette: Meso-NH version
!  P. Wautelet 21/02/2019: add copyright notice + use INT64 for 64-bits integers
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
CHARACTER(LEN=*), INTENT(IN) :: CDFILE
CHARACTER(LEN=*), INTENT(IN) :: CDSTATE
INTEGER(KIND=INT64), INTENT(OUT) :: KNUMER
!
! II. Local variables declaration
INTEGER, PARAMETER :: JPMAXLOGICALUNITNUMBER=5000
INTEGER(KIND=LFI_INT) :: IRETURNCODE
LOGICAL :: LLEXISTS, LLOPEN
INTEGER(KIND=LFI_INT) :: IRECORDNUMBER
INTEGER(KIND=LFI_INT) :: INUMER
!
! III. File opening
!
! III.a Search for an available logical unit
INUMER=0
LLEXISTS=.FALSE.
LLOPEN=.TRUE.
IRETURNCODE=0
DO WHILE(INUMER.LT.JPMAXLOGICALUNITNUMBER .AND. (LLOPEN .OR. .NOT. LLEXISTS))
  INUMER=INUMER+1
  INQUIRE(UNIT=INUMER, EXIST=LLEXISTS, OPENED=LLOPEN)
ENDDO
IF(LLOPEN .OR. .NOT. LLEXISTS) THEN
  IRETURNCODE=-999
ENDIF
!
#ifdef __GFORTRAN__
! III.b (Re)-init of libgfortran to enable big_endian file reading
!**** *** ** * only gfortran will work with this * ** *** ****
CALL INIT_GFORTRAN_BIG_ENDIAN()
#endif
!
! III.c LFI file opening
CALL LFIOUV(IRETURNCODE, INUMER, .TRUE., CDFILE, CDSTATE, .FALSE.,&
           &.FALSE., INT(0, LFI_INT), INT(1, LFI_INT), IRECORDNUMBER)
IF(IRETURNCODE/=0)THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFIOUV', '')
ENDIF
!
#ifdef __GFORTRAN__
! III.d (Re)-init of libgfortran to enable native endianess file reading
!**** *** ** * only gfortran will work with this * ** *** ****
CALL INIT_GFORTRAN_NATIVE_ENDIAN()
#endif
!
KNUMER=INT(INUMER, 8)
KRETURNCODE=INT(IRETURNCODE,8)

END SUBROUTINE WLFIOUV

!______________________________________________________________________

SUBROUTINE WLFINAF(KRETURNCODE, KNUMER, KNALDO, KNTROU, KNARES, KNAMAX)
! ** PURPOSE
!    Wrapper to LFINAF
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KNUMER: logical unit number associated to file
!    KNALDO: Number of actual logical data records (holes excluded)
!    KNTROU: Number of logical records which are holes
!    KNARES: Number of logical records which can be written in the reserved part of index (holes included)
!    KNAMAX: Maximum number of logical records which one can write on logical unit
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!    8 nov 2018, S. Riette: Meso-NH version
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: KNUMER
INTEGER(KIND=INT64), INTENT(OUT) :: KNALDO, KNTROU, KNARES, KNAMAX
!
! II. Local variables declaration
INTEGER(KIND=LFI_INT) :: IRETURNCODE
INTEGER(KIND=LFI_INT) :: INUMER
INTEGER(KIND=LFI_INT) :: INALDO, INTROU, INARES, INAMAX
!
! III. LFINAF call
INUMER=INT(KNUMER, KIND(INUMER))
CALL LFINAF(IRETURNCODE, INUMER, INALDO, INTROU, INARES, INAMAX)
IF(IRETURNCODE/=0)THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFINAF', '')
ENDIF
KRETURNCODE=INT(IRETURNCODE,8)
KNALDO=INT(INALDO, 8)
KNTROU=INT(INTROU, 8)
KNARES=INT(INARES, 8)
KNAMAX=INT(INAMAX, 8)
!
END SUBROUTINE WLFINAF

!______________________________________________________________________

SUBROUTINE WLFIPOS(KRETURNCODE, KNUMER)
! ** PURPOSE
!    Wrapper to LFIPOS
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KNUMER: logical unit number associated to file
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!    8 nov 2018, S. Riette: Meso-NH version
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: KNUMER
!
! II. Local variables declaration
INTEGER(KIND=LFI_INT) :: IRETURNCODE
INTEGER(KIND=LFI_INT) :: INUMER
!
! III. LFIPOS call
INUMER=INT(KNUMER, KIND(INUMER))
CALL LFIPOS(IRETURNCODE, INUMER)
IF(IRETURNCODE/=0)THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFIPOS', '')
ENDIF
KRETURNCODE=INT(IRETURNCODE,8)
!
END SUBROUTINE WLFIPOS

!______________________________________________________________________

SUBROUTINE WLFICAS(KRETURNCODE, KNUMER, CDNOMA, KLONG, KPOSEX, LDAVAN)
! ** PURPOSE
!    Wrapper to LFICAS
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KNUMER: logical unit number associated to file
!    CDNOMA: name of next record
!    KLONG: length of next record
!    KPOSEX: position in file of the first word of next record
!    LDAVAN: true if one must move forward the pointer
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!                             use of true logical instead of integer
!    8 nov 2018, S. Riette: Meso-NH version
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: KNUMER
CHARACTER(LEN=16), INTENT(OUT) :: CDNOMA
INTEGER(KIND=INT64), INTENT(OUT) :: KLONG, KPOSEX
LOGICAL, INTENT(IN) :: LDAVAN
!
! II. Local variables declaration
INTEGER(KIND=LFI_INT) :: IRETURNCODE
INTEGER(KIND=LFI_INT) :: INUMER
INTEGER(KIND=LFI_INT) :: ILONG, IPOSEX
!
! III. LFICAS call
INUMER=INT(KNUMER, KIND(INUMER))
CALL LFICAS(IRETURNCODE, INUMER, CDNOMA, ILONG, IPOSEX, LDAVAN)
IF(IRETURNCODE/=0)THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFICAS', '')
ENDIF
KRETURNCODE=INT(IRETURNCODE,8)
KLONG=INT(ILONG, 8)
KPOSEX=INT(IPOSEX, 8)
!
END SUBROUTINE WLFICAS

!______________________________________________________________________

SUBROUTINE WLFINFO(KRETURNCODE, KNUMER, CDNOMA, KLONG, KPOSEX)
! ** PURPOSE
!    Wrapper to LFINFO
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KNUMER: logical unit number associated to file
!    CDNOMA: name of record
!    KLONG: length of record
!    KPOSEX: position in file of the first word of next record

! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!    8 nov 2018, S. Riette: Meso-NH version
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: KNUMER
CHARACTER(LEN=16), INTENT(IN) :: CDNOMA
INTEGER(KIND=INT64), INTENT(OUT) :: KLONG, KPOSEX
!
! II. Local variables declaration
INTEGER(KIND=LFI_INT) :: IRETURNCODE
INTEGER(KIND=LFI_INT) :: INUMER
INTEGER(KIND=LFI_INT) :: ILONG, IPOSEX
!
! III. LFINFO call
INUMER=INT(KNUMER, KIND(INUMER))
CALL LFINFO(IRETURNCODE, INUMER, CDNOMA, ILONG, IPOSEX)
IF(IRETURNCODE/=0)THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFINFO', '')
ENDIF
KRETURNCODE=INT(IRETURNCODE,8)
KLONG=INT(ILONG, 8)
KPOSEX=INT(IPOSEX, 8)
!
END SUBROUTINE WLFINFO

!______________________________________________________________________

SUBROUTINE WLFILEC(KRETURNCODE, KNUMER, CDNOMA, KLONG, LDABORT, KTAB)
! ** PURPOSE
!    Wrapper to LFILEC
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KNUMER: logical unit number associated to file
!    CDNOMA: name of record
!    KLONG: length of record
!    LDABORT: must we raise an exception on error -21 ?
!    KTAB: integer array read

! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!                             use of true logical instead of integer
!    8 nov 2018, S. Riette: Meso-NH version
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: KNUMER
CHARACTER(LEN=16), INTENT(IN) :: CDNOMA
INTEGER(KIND=INT64), INTENT(IN) :: KLONG
LOGICAL, INTENT(IN) :: LDABORT
INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(KLONG) :: KTAB
!
! II. Local variables declaration
INTEGER(KIND=LFI_INT) :: IRETURNCODE
INTEGER(KIND=LFI_INT) :: INUMER, ILONG
INTEGER(KIND=LFI_INT) :: ITOTLONG, IPOSEX
INTEGER(KIND=INT64), ALLOCATABLE :: KTABTOT(:)
!
! III. LFILEC call
INUMER=INT(KNUMER, KIND(INUMER))
ILONG=INT(KLONG, KIND(ILONG))
!
!Because NERFAG cannot be changed easily, we read the entire article
!even if only a part is needed, otherwise NERFAG=2 would be sufficient
CALL WLFINFO(IRETURNCODE, INUMER, CDNOMA, ITOTLONG, IPOSEX)
IF(ILONG .LT. ITOTLONG) THEN
  ALLOCATE(KTABTOT(ITOTLONG))
  CALL LFILEC(IRETURNCODE, INUMER, CDNOMA, KTABTOT, ITOTLONG)
ELSE
  CALL LFILEC(IRETURNCODE, INUMER, CDNOMA, KTAB, ILONG)
ENDIF
IF (IRETURNCODE/=0 .AND. .NOT. (IRETURNCODE==-21 .AND. .NOT. LDABORT)) THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFILEC', '')
  KRETURNCODE=INT(IRETURNCODE,8)
ELSE
  KRETURNCODE=INT(0,8)
ENDIF
IF(ILONG .LT. ITOTLONG) THEN
  KTAB(:)=KTABTOT(1:ILONG)
  DEALLOCATE(KTABTOT)
ENDIF
!
END SUBROUTINE WLFILEC

!_________________________________________________________________________________________________

SUBROUTINE WGET_COMPHEADER(KSIZE, KDATA, KLONG, KLONU, KTYPECOMP)
! ** PURPOSE
!    Wrapper to GET_COMPHEADER
!
! ** DUMMY ARGUMENTS
!    KSIZE: Size of KDATA
!    KDATA: (part of) integer array read from record
!    KLONG: length of compressed data
!    KLONU: length of uncompressed data
!    KTYPECOMP: type of compression
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(IN) :: KSIZE
INTEGER(KIND=INT64), INTENT(IN), DIMENSION(KSIZE) :: KDATA
INTEGER(KIND=INT64), INTENT(IN) :: KLONG
INTEGER(KIND=INT64), INTENT(OUT) :: KLONU
INTEGER(KIND=INT64), INTENT(OUT) :: KTYPECOMP
!
! II. Local variables declaration
INTEGER :: ILONG
INTEGER :: ILONU, ITYPECOMP
!
! III. GET_COMPHEADER call
#ifdef MNH_COMPRESS
ILONG=KLONG
CALL GET_COMPHEADER(KDATA, ILONG, ILONU, ITYPECOMP)
KLONU=INT(ILONU, 8)
KTYPECOMP=INT(ITYPECOMP, 8)
#else
print*, "Error: code was compiled without COMPRESS support, please define MNH_COMPRESS"
KLONU=INT(-1, 8)
KTYPECOMP=INT(-1, 8)
#endif
!
END SUBROUTINE WGET_COMPHEADER

!_________________________________________________________________________________________________

SUBROUTINE WCOMPRESS_FIELD(KTAB, KX, KY, KSIZEDECOMP, KSIZECOMP)
! ** PURPOSE
!    Wrapper to COMPRESS_FIELD
!
! ** DUMMY ARGUMENTS
!    KTAB: decompressed integer array (IN)
!          compressed data integer array (OUT)
!    KX, KY: x and y dimensions
!    KSIZEDECOMP: size of decompressed data
!    KSIZECOMP: size of compressed integer array
!
! ** AUTHOR
!    16 July 2015, S. Riette
!
! ** MODIFICATIONS
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(IN) :: KX, KY, KSIZEDECOMP
INTEGER(KIND=INT64), INTENT(INOUT), DIMENSION(KSIZEDECOMP) :: KTAB
INTEGER(KIND=INT64), INTENT(OUT) :: KSIZECOMP
!
! II. Local variables declaration
INTEGER :: ISIZEDECOMP, ISIZECOMP, IX, IY
!
! III. COMPRESS_FIELD call
#ifdef MNH_COMPRESS
ISIZEDECOMP=KSIZEDECOMP
IX=KX
IY=KY
CALL COMPRESS_FIELD(KTAB, IX, IY, ISIZEDECOMP, ISIZECOMP)
KSIZECOMP=ISIZECOMP
#else
print*, "Error: code was compiled without COMPRESS support, please define MNH_COMPRESS"
KSIZECOMP=INT(-1, 8)
#endif
!
END SUBROUTINE WCOMPRESS_FIELD

!_________________________________________________________________________________________________

SUBROUTINE WDECOMPRESS_FIELD(KSIZE, KCOMP, KTYPECOMP, KLDECOMP, KDECOMP)
! ** PURPOSE
!    Wrapper to DECOMPRESS_FIELD
!
! ** DUMMY ARGUMENTS
!    KSIZE: size of KCOMP
!    KCOMP: compressed integer array
!    KTYPECOMP: type of compression
!    KDECOMP: decompressed data integer array
!    KLDECOMP: length of decompressed data
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(IN) :: KSIZE
INTEGER(KIND=INT64), INTENT(IN), DIMENSION(KSIZE) :: KCOMP
INTEGER(KIND=INT64), INTENT(IN) :: KTYPECOMP
INTEGER(KIND=INT64), INTENT(IN) :: KLDECOMP
INTEGER(KIND=INT64), INTENT(OUT), DIMENSION(KLDECOMP) :: KDECOMP
!
! II. Local variables declaration
INTEGER :: ITYPECOMP, ILDECOMP
!
! III. DECOMPRESS_FIELD call
#ifdef MNH_COMPRESS
ILDECOMP=KLDECOMP
ITYPECOMP=KTYPECOMP
CALL DECOMPRESS_FIELD(KDECOMP, ILDECOMP, KCOMP, SIZE(KCOMP,1), ITYPECOMP)
#else
print*, "Error: code was compiled without COMPRESS support, please define MNH_COMPRESS"
KDECOMP(:)=-1
#endif
!
END SUBROUTINE WDECOMPRESS_FIELD

!_________________________________________________________________________________________________

SUBROUTINE WLFIFER(KRETURNCODE, KNUMER, CDSTTC)
! ** PURPOSE
!    Close a LFI file
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KNUMER: logical unit number associated to file
!    CDSTTC: close status ('KEEP', 'SCRATCH', 'DELETE')
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!    8 nov 2018, S. Riette: Meso-NH version
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: KNUMER
CHARACTER(LEN=7), INTENT(IN) :: CDSTTC
!
! II.  Local variables declaration
INTEGER(KIND=LFI_INT) :: IRETURNCODE
INTEGER(KIND=LFI_INT) :: INUMER
!
! III. LFIFER call
INUMER=INT(KNUMER, KIND(INUMER))
CALL LFIFER(IRETURNCODE, INUMER, CDSTTC)
IF(IRETURNCODE/=0)THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFIFER', '')
ENDIF
KRETURNCODE=INT(IRETURNCODE,8)
!
END SUBROUTINE WLFIFER

!_________________________________________________________________________________________________

SUBROUTINE WLFIECR(KRETURNCODE, KNUMER, CDNOMA, KSIZE, KTAB)
! ** PURPOSE
!    Wrapper to LFIECR
!
! ** DUMMY ARGUMENTS
!    KRETURNCODE: error code
!    KNUMER: logical unit number associated to file
!    CDNOMA: name of field to write
!    KSIZE: Size of KTAB
!    KTAB: integer array to write
!
! ** AUTHOR
!    9 April 2014, S. Riette
!
! ** MODIFICATIONS
!    26 sept 2014, S. Riette: use 64bits LFI subroutines
!    8 nov 2018, S. Riette: Meso-NH version
!
! I. Dummy arguments declaration
use iso_fortran_env, only: INT64
IMPLICIT NONE
INTEGER(KIND=INT64), INTENT(OUT) :: KRETURNCODE
INTEGER(KIND=INT64), INTENT(IN) :: KNUMER
CHARACTER(LEN=16), INTENT(IN) :: CDNOMA
INTEGER(KIND=INT64), INTENT(IN) :: KSIZE
INTEGER(KIND=INT64), INTENT(IN), DIMENSION(KSIZE) :: KTAB
!
! II.  Local variables declaration
INTEGER(KIND=LFI_INT) :: IRETURNCODE
INTEGER(KIND=LFI_INT) :: INUMER
INTEGER(KIND=LFI_INT) :: ISIZE
!
! III. LFIECR call
INUMER=INT(KNUMER, KIND(INUMER))
ISIZE=INT(KSIZE, KIND(ISIZE))
CALL LFIECR(IRETURNCODE, INUMER, CDNOMA, KTAB, ISIZE)
IF(IRETURNCODE/=0)THEN
  CALL LFIENG(INUMER, INT(0, LFI_INT), IRETURNCODE, .FALSE., '', 'LFIECR', '')
ENDIF
KRETURNCODE=INT(IRETURNCODE,8)
!
END SUBROUTINE WLFIECR

!_________________________________________________________________________________________________
