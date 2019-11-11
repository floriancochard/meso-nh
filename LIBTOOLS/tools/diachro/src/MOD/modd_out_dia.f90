!     ######spl
      MODULE  MODD_OUT_DIA
!     ####################
!
!!****  *MODD_OUT_DIA* - defines a logical unit number for printed outputs
!!
!!    PURPOSE
!!    -------
!       So far, this declarative module is a garbage box containing items
!     fitting nothing else...
!       Content:
!         - logical unit number for the printed output,
!         - size of the matrix section to be displayed in the MESO-NH
!           field arrays,
!         - indexes to locate the displayed MESO-NH column in the 
!           the "radio-sounding" mode,
!         - filename prefix of the LFI file to be displayed.
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_FIELD1_CV2D), to appear in 1994 
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original        01/02/96
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


INTEGER                 :: NLUOUTD                !  Logical unit number for
!INTEGER,DIMENSION(50)  :: NLUOUTD                !  Logical unit number for
                                                 !  printed outputs


CHARACTER(LEN=16)                :: CLUOUTD      !  Names for printed outputs
!CHARACTER(LEN=16),DIMENSION(50)  :: CLUOUTD      !  Names for printed outputs

CHARACTER(LEN=32),DIMENSION(100)  :: CLFIFMD      !  Full names of the ".lfi"
CHARACTER(LEN=32),DIMENSION(100)  :: CDESFMD      !  and ".des" files to be
                                                 !  processed

CHARACTER(LEN=28),DIMENSION(100) :: CNAMFILED      !  Filename prefix of the files
                                                 !  to be processed
!
END MODULE MODD_OUT_DIA
