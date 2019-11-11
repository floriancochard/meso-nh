!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_out.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE  MODD_OUT
!     #################
!
!!****  *MODD_OUT* - defines a logical unit number for printed outputs
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
!!     original        02/06/94
!!     updated   PM    21/11/94
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


INTEGER           :: NLUOUT                      !  Logical unit number for
                                                 !  printed outputs

INTEGER           :: NIMAXT, NJMAXT, NKMAXT      !  Size of the displayed 
                                                 !  section of the MESO-NH
                                                 !  field arrays
INTEGER           :: NNAMRS                      !  =0  --> RS
                                                 !  =1  =/= RS

INTEGER           :: NIRS,  NJRS                 !  Grid indexes to locate
                                                 !  a "radiosounding" point

INTEGER :: NRRM     ! Total number of water variables at time t
INTEGER :: NRRT     ! Total number of water variables at time t

CHARACTER(LEN=20) :: CNAMRS                      !  Contains 'RS' in case of RS
                                                 !           something else in
                                                 !           the others cases
CHARACTER(LEN=32) :: CLFIFM, CDESFM              !  Full names of the ".lfi"
                                                 !  and ".des" files to be
                                                 !  processed

CHARACTER(LEN=28) :: CNAMFILE                    !  Filename prefix of the files
                                                 !  to be processed
!
END MODULE MODD_OUT
