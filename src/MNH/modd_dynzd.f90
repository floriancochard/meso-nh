!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 14:00:36
!-----------------------------------------------------------------
!     #################
      MODULE MODD_DYNZD
!     #################
!
!!****  *MODD_DYN - declaration of dynamic control variables 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the dynamic
!     control variables.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DYNn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/05/94                      
!!      Modifications 16/11/94   (Lafore+Pinty)  For NUM_DIFF
!!      Modifications 06/01/95   (Lafore)        For LSTEADY_DMASS
!!      Modifications 28/07/96   (Masson)        Supress LSTEADY_DMASS
!!      Modifications 15/03/98   (Stein)         Add LHO_RELAX for each variables 
!!      Modifications 22/01/01   (Gazen)         Add LHORELAX_SVC2R2, _SVCHEM, _SVLG
!!      Modifications 29/11/02   (Pinty)         Add  LHORELAX_SVC1R3, _SVELEC
!!      Modifications 03/11/04   (Zängl)         Add fields for truly horizontal diffusion
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
IMPLICIT NONE
!

! Additional variables needed for truly horizontal diffusion (G. Zängl)

LOGICAL, SAVE   :: LZDIFFU     ! Logical switch if modified diffusion is used

TYPE TYPE_ZDIFFU_HALO2
!  
!  data and all is dimension
!
   REAL, DIMENSION(:,:,:), POINTER :: XZZ=>NULL()   ! height with halo2

!  Interpolation coefficients
   REAL, DIMENSION(:,:,:), POINTER :: XRKIP1=>NULL(),XRKIP2=>NULL(),&
                                      XRKIM1=>NULL(),XRKIM2=>NULL(),&
                                      XRKJP1=>NULL(),XRKJP2=>NULL(),&
                                      XRKJM1=>NULL(),XRKJM2=>NULL()
!
!  Reduction factors for diffusion coefficients near the surface
   REAL    , DIMENSION(:,:)   , POINTER  :: XREDFACI=>NULL(),XREDFACJ=>NULL()
!
! lowest model level for which z-diffusion can be calculated
! on all grid points without intersecting the topography
   INTEGER                               :: NZDLB
!
! lowest model level for which z-diffusion can be calculated
! on the local grid point without intersecting the topography
   INTEGER , DIMENSION(:,:)   , POINTER  :: NZDI=>NULL(),NZDJ=>NULL()
!
END TYPE TYPE_ZDIFFU_HALO2

END MODULE MODD_DYNZD
