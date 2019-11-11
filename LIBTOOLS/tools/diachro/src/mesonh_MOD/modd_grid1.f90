!-----------------------------------------------------------------
!     ##################
      MODULE MODD_GRID1
!     ##################
!
!!****  *MODD_GRID1* - declaration of grid variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the variables
!     describing the grid. 
!    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_GRIDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94                      
!!      J. Stein    15/11/95  add the slope angle
!!      V. Ducrocq   13/08/98  // : add XLATOR_ll and XLONOR_ll       
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
REAL,SAVE :: XLONOR,XLATOR  ! Longitude and latitude of the Origine point
                            !  for the conformal projection of the sub-domain (//)
REAL,SAVE :: XLONOR_ll,XLATOR_ll  ! Longitude and latitude of the Origine point
                            !  for the conformal projection of the domain 
REAL,SAVE, DIMENSION(:,:), ALLOCATABLE :: XLON,XLAT ! Longitude and latitude  
!
REAL,SAVE, DIMENSION(:),   ALLOCATABLE :: XXHAT   ! Position x in the 
                                         ! conformal or cartesian plane
REAL,SAVE, DIMENSION(:),   ALLOCATABLE :: XYHAT   ! Position y in the 
                                         ! conformal or cartesian plane
REAL,SAVE, DIMENSION(:),   ALLOCATABLE :: XDXHAT  ! horizontal stretching in x
REAL,SAVE, DIMENSION(:),   ALLOCATABLE :: XDYHAT  ! horizontal stretching in y
REAL,SAVE, DIMENSION(:,:), ALLOCATABLE :: XMAP    ! Map factor 
!
REAL,SAVE, DIMENSION(:,:),   ALLOCATABLE :: XZS   ! orography
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XZZ   ! height z 
REAL,SAVE, DIMENSION(:),     ALLOCATABLE :: XZHAT ! height level without orography
!
REAL, DIMENSION(:,:)  , ALLOCATABLE :: XDIRCOSXW,XDIRCOSYW,XDIRCOSZW 
                                               ! director cosinus of the normal 
                                               ! to the ground surface 
!  
REAL,SAVE, DIMENSION(:,:),  ALLOCATABLE  ::  XCOSSLOPE  ! cosinus of the angle
                                 ! between i and the slope vector
REAL,SAVE, DIMENSION(:,:),  ALLOCATABLE  ::  XSINSLOPE  ! sinus of the angle
                                 ! between i and the slope vector
!
!* quantities for SLEVE vertical coordinate
LOGICAL,SAVE                             :: LSLEVE    ! Logical for SLEVE coordinate
REAL,SAVE                                :: XLEN1     ! Decay scale for smooth topography
REAL,SAVE                                :: XLEN2     ! Decay scale for small-scale topography deviation
REAL,SAVE, DIMENSION(:,:),   ALLOCATABLE :: XZSMT   ! smooth orography for SLEVE coordinate
!
END MODULE MODD_GRID1
