!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######spl
      MODULE MODD_PGDSSOWORK
!     ######################
!
!!****  *MODD_PGDSSOWORK* - declaration of work arrays and variables
!!                          for SSO computation
!!
!!    PURPOSE
!!    -------  
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/05/97                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!
!
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZSSQO   ! mean of orography in a
!                                               ! SSO subgrid square from
!                                               ! ZMAXSSQ averaged values
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GSSQO   ! presence of data in a SSO
!                                               ! subgrid square
INTEGER :: ISSQOX                               ! number of SSO subgrid squares
INTEGER :: ISSQOY                               ! in each direction in grid mesh
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_PGDSSOWORK
