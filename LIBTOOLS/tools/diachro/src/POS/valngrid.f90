!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.valngrid.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE VALNGRID(HCAR)
!     #########################
!
!!****  *VALNGRID* - Selects the NGRID value (alias KGRID or IGRID
!!                   or NMGRID)
!!
!!    PURPOSE
!!    -------
!       Given only the name of a variable, returns the corresponding 
!      NGRID value, and calculates the true altitude array for this
!      grid location.
!
!!**  METHOD
!!    ------
!!     
!!      The name is given as a character string, the NGRID value is found
!!     by searching the LFIFM record for this string. 
!!      Next, the relevant altitude array is built by a call to COMPCOORD.
!!
!!    EXTERNAL
!!    --------
!!      COMPCOORD : computes the true sea-level altitude corresponding to the
!!     current NGRID selection.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!         NMGRID    : Current MESO-NH grid indicator
!!
!!      Module MODD_OUT       : Defines a log. unit for printing
!!         CNAMFILE : filename prefix of the FM files to be processed
!!         NLUOUT   : Logical unit number for printed output
!!
!!      Module MODD_DIM1       : Contains dimensions
!!         NIMAX,NJMAX,NKMAX :  x, y, and z array dimensions
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!         JPHEXT   : Horizontal external points number
!!         JPVEXT   : Vertical external points number
!!
!!      Module MODD_LUNIT1     : Declares names and log. unit of files
!!         CLUOUT   : Name of output_listing file
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   01/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_NMGRID
USE MODD_DIM1
USE MODD_LUNIT1
USE MODD_PARAMETERS
USE MODD_OUT
USE MODI_FMREAD

IMPLICIT NONE
!
!*       0.1   Declarations and dummy arguments
!
CHARACTER(LEN=*)  :: HCAR       ! name of the requested variable, as a string                  
!
!*       0.2   Local variables
!
INTEGER             :: I3D                        ! size of 3D   arrays
INTEGER             :: IIU, IJU, IKU              ! array sises

INTEGER             :: ILENG, IGRID,ILENCH,IRESP  !   File 
CHARACTER (LEN=16)  :: YRECFM                     ! management
CHARACTER (LEN=100) :: YCOMMENT                   ! variables   
!
CHARACTER(LEN=10) :: YCAR                         ! work array

REAL, DIMENSION(:,:,:),ALLOCATABLE,SAVE :: ZS3D   ! 3D array used to read  data
                                                  ! in initial file 
!
!*       0.3   String justified left to avoid trouble
!
!WRITE(NLUOUT,*)' HCAR ',HCAR
YCAR=ADJUSTL(HCAR)
!
!-------------------------------------------------------------------------------
!
!*       1.    SETS ARRAY SIZES AND ALLOCATIONS
!              --------------------------------
!
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT

IF(.NOT.ALLOCATED(ZS3D))THEN
  ALLOCATE(ZS3D(IIU,IJU,IKU))
END IF

I3D=IIU*IJU*IKU
!
!-------------------------------------------------------------------------------
!
!*       2.   SEARCHES THE LFIFM FILE FOR THE STRING
!               AND COMPUTES APPROPRIATE TRUE-ALTITUDES
!             -------------------------------- --------
!
YRECFM = YCAR
CALL FMREAD(CNAMFILE,YRECFM,CLUOUT,I3D,ZS3D,IGRID,ILENCH,YCOMMENT,IRESP)
IF(IRESP.EQ.-47)THEN
  NMGRID=1
ELSE
  NMGRID=IGRID
END IF
CALL COMPCOORD(NMGRID)
DEALLOCATE(ZS3D)
!
!-------------------------------------------------------------------------------!
!*       3.   EXIT
!             ----
RETURN
!
END SUBROUTINE VALNGRID
