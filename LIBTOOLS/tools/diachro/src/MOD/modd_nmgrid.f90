!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_nmgrid.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_NMGRID
!     ##################
!
!!****  *MODD_NMGRID* - Global variable  NMGRID declaration
!!
!!    PURPOSE
!!    -------
!!      This declarative module defines the NMGRID global variable, which
!!    stores the value of the grid indicator for the current displayed field
!!    (local alias IGRID or KGRID). The grid indicator is the number of the
!!    grid where the displayed variable is located in the MESO-NH model. Seven
!!    different grids are used, so far. See Book-1 for grid definitions.
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_FIELD1_CV2D), to appear in 1994 
!!
!!     The 7 MESO-NH grid types are defined in:
!!      
!!      - Asencio N. et al., 1994, "Le projet de modele non-hydrostatique
!!        commun CNRM-LA, specifications techniques", 
!!        Note CNRM/GMME, 26, 139p, (pages 39 to 43).
!!
!!      - Fischer C., 1994, "File structure and content in the Meso-NH 
!!        model", Meso-nh internal note, CNRM/GMME,  July 5.
!!       
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
INTEGER           :: NMGRID    ! Current MESO-NH grid indicator
!
END MODULE MODD_NMGRID
