!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_super.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE MODD_SUPER
!     #################
!
!!****  *MODD_SUPER* - This declaration module defines variables controlling
!!                     the overlay of several successive plots over either
!!                     an horizontal map or a vertical cross-section.
!!                     
!!    PURPOSE
!!    -------
!       To control the possibility of plot overlay, two global variables are
!     defined. LSUPER is a logical specifying if the overlay option is
!     activated for the current plot. NSUPER is an integer giving the rank
!     of the current plot in the overlay sequence. The first plot of the
!     sequence, i.e. the background plot, is given the rank NSUPER=1.
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!     The principle of overlay handling in TRACE is detailed in:
!!     - Book3 of the TRACE volume of the Meso-NH user manual.
!!
!!     The technicalities are found in:
!!     - Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_SUPER)
!!
!!       
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      23/11/94                      
!!      Updated  PM   24/11/94
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
LOGICAL,SAVE   :: LSUPER  ! =.T. --> plot overlay is active
                          ! =.F. --> plot overlay is not active
!
INTEGER,SAVE   :: NSUPER  ! Rank of the current plot in the overlay
                          ! sequence. The initial plot is rank 1.
!
END MODULE MODD_SUPER
