!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_coord.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE MODD_COORD
!     #################
!
!!****  *MODD_COORD* - Declaration of the TRACE arrays giving the gridpoint
!!                     coordinates for the 7 MESO-NH grid types.
!!
!!    PURPOSE
!!    -------
!       This declarative module defines a set of arrays containing:
!
!     - XHAT, YHAT, ZHAT coordinate values for all the MESO-NH grids 
!           --> XXX(:,:) XXY(:,:) XXZ(:,:)
!     - Meshsize values for all the available MESO-NH grids:
!           --> XXDXHAT(:,:) XXDYHAT(:,:) XXDZHAT(:,:)
!     - Oblique meshsize values along the abscissa (horizontal) axis 
!       of the oblique vertical cross-sections (for all the MESO-NH grids):
!           --> XXDS(:,:)
!     - Oblique abscissa values for the gridpoints along the horizontal 
!       axis of the oblique vertical cross-sections(for all the MESO-NH grids):
!           --> XDS(:,:)
!     - X- and Y- projections on the MESO-NH axes directions for the gridpoints 
!       along the oblique vertical cross-sections (for all the MESO-NH grids): 
!           --> XDSX(:,:) XDSY(:,:)
!     - Interpolated topography for all the available MESO_NH grids:
!           --> XXZS(:,:,:)
!     
!     In all the forecoming arrays, the last index is the grid indicator,
!     NGRID, i. e. the number of the grid where the displayed variable is 
!     located. Seven grids are available so far, see the MESO-NH Book-1 
!     for definitions. The local name for this grid indicator may be IGRID,
!     or NMGRID, or KGRID according to the context. 
!
!     Lengthes are given in meters.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
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
IMPLICIT NONE
!
! XXDXHAT, XXDYHAT, Mesh size arrays (meters), last index is the grid indicator 
! XXDZHAT        
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XXDXHAT, XXDYHAT, XXDZHAT

! XXX, XXY,   Values of XHAT, YHAT, ZHAT (meters) for the different grids, 
! XXZ         as given by NGRID (second index, grid indicator) 
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XXX, XXY, XXZ

! XXDS        Mesh size (meters) along the horizontal axis of an oblique  
!             vertical cross-section, for all the grids (given by NGRID, 
!             second index)
! XDS         Abscissa array along the horizontal axis of an oblique vertical
!             cross-section (meters), for all the grids (given by NGRID, 
!             second index)
! XDSX, XDSY  Projections on the MESO-NH cartesian axes of the XDS oblique
!             abscissa (meters), for all the grids (given by NGRID, second
!             index)
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XXDS, XDS, XDSX, XDSY

! XXZS        Terrain topography (meters) interpolated at the NGRID gridpoint
!             location (for all the possible grids, given by the third index)
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: XXZS
				   

END MODULE MODD_COORD
