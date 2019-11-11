!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_cvert.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_CVERT
!     ###################
!
!!****  *MODD_CVERT* - Declares work arrays for vertical cross-sections
!!
!!    PURPOSE
!!    -------
!       For vertical cross-sections only, this declarative module declares 
!     the arrays containing the sea-level altitudes and the model topography 
!     of the oblique cross-section points.     
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_CVERT) 
!!       
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/06/94                      
!!      Updated   PM   17/11/94  
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XWORKZ ! Sea-level altitude array 
                                                 ! (meters)
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: XWZ    ! Topography array (meters)

END MODULE MODD_CVERT
