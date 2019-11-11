!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_DIM1
!     ##################
!
!!****  *MODD_DIM1* - declaration of dimensions
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the dimensions 
!     of the data arrays.   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DIMn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94     
!!      Modifications 13/08/98 (V. Ducrocq) // NIINF .. NJSUP are no more used in the init part                
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
INTEGER,SAVE :: NIMAX,NJMAX,NKMAX  !  Dimensions respectively  in x , 
                              ! y ,  z directions of the physical sub-domain.
INTEGER,SAVE :: NIMAX_ll,NJMAX_ll  !  Dimensions respectively  in x and y
                                   ! directions of the physical domain
INTEGER,SAVE :: NIINF, NISUP       !  Lower bound and upper bound of the arrays 
                                   ! in x direction 
INTEGER,SAVE :: NJINF, NJSUP       !  Lower bound and upper bound of the arrays 
                                   ! in y direction
!
END MODULE MODD_DIM1
