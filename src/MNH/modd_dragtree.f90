!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!!
!!    #####################
      MODULE MODD_DRAGTREE
!!    #####################
!!
!!*** *MODD_DRAGTREE*
!!
!!    PURPOSE
!!    -------
!       Declaration to take into account tree drag in Meso-NH              
!              instead of SURFEX. 
!!
!!**  AUTHOR
!!    ------
!!    C.Lac                   *CNRM*
!
!!    MODIFICATIONS
!!    -------------
!!    Original 30/06/11
!!    06/16  (C.Lac) Add droplet deposition
!-----------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
LOGICAL    ::     LDRAGTREE    ! flag used to  take into account tree drag in 
!                              ! the atmospheric model instead of SURFEX.
LOGICAL    ::     LDEPOTREE    ! flag for droplet deposition on trees
!
REAL       ::     XVDEPOTREE   ! Droplet deposition velocity
!
END MODULE MODD_DRAGTREE
