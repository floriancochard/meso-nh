!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_title.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE MODD_TITLE
!     #################
!
!!****  *MODD_TITLE* - Declares heading variables for the plots
!!
!!    PURPOSE
!!    -------
!       This declarative module defines a character variable containing
!     the heading title of the current plot, and the rank of the current
!     plot from the start of the TRACE session.
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_TITLE), to appear in 1994 
!!       
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/06/94                      
!!      Updated  PM 22/11/94
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
INTEGER           :: NCONT    ! Current plot number

CHARACTER(LEN=110) :: CLEGEND  ! Current plot heading title
CHARACTER(LEN=100) :: CLEGEND2 ! Current plot heading title

END MODULE MODD_TITLE
