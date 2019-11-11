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
!     ##############
      MODULE MODD_LG
!     ##############
!
!!****  *MODD_RAIN_LG* - declaration of variables used to manage lagrangian variables
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_LG)
!!          
!!    AUTHOR
!!    ------
!!	P. Jabouille  *CNRM*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/06/2001
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
CHARACTER(LEN=10),DIMENSION(3),PARAMETER :: CLGNAMES=(/'LGX','LGY','LGZ'/)
                                       ! basenames of the lagrangian articles stored
                                       ! in the binary files
REAL,PARAMETER :: XLG1MIN=-1.E+9, XLG2MIN=-1.E+9, XLG3MIN=0.
                                       ! minimum values for lagrangian variables
REAL,PARAMETER :: XLGSPEED=10.         ! advective velocity 
!
END MODULE MODD_LG
!
!
