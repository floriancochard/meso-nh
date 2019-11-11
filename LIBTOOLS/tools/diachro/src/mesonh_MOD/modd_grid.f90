!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_grid.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE MODD_GRID
!     #################
!
!!****  *MODD_GRID* - declaration of grid variables for all models
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the variables
!     describing the grid for all models. 
!    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_GRID)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94                      
!!      V. Masson   nov 2004  : add XLATORI and XLONORI    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
REAL,SAVE :: XLON0,XLAT0    ! Reference longitude and latitude 
                            !  for the conformal projection
REAL,SAVE :: XBETA,XRPK     ! Rotation angle and projection parameter
                            !  for the conformal projection
REAL,SAVE :: XLONORI,XLATORI ! Longitude and latitude of the point
                             ! of coordinates x=0, y=0
                             ! for the conformal projection
!  
END MODULE MODD_GRID
