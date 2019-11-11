!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------

!!    ######################
      SUBROUTINE INIT_SALT
!!    ######################
!     PURPOSE
!!     -------
!!
!!     initialization of variables for the sea salt scheme
!!
!!     METHOD
!!     ------
!!
!!
!!     REFERENCE
!!     ---------
!!     none
!!
!!
!!     AUTHOR
!!     ------
!!     Marine Claeys (CNRM)
!!
!!     MODIFICATIONS
!!     -------------
!!
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!
USE MODD_SALT
!
IMPLICIT NONE

IF(NMODE_SLT == 5) THEN

!JPSALTORDER = (/5, 4, 3, 2, 1 /)
!Initial dry number median radius (um) from Ova et al., 2014
XINIRADIUS_SLT=  (/0.009, 0.021, 0.045, 0.115, 0.415/)
!Initial, standard deviation from  Ova et al., 2014
XINISIG_SLT =  (/ 1.37, 1.5, 1.42, 1.53, 1.85 /)
!Minimum allowed number concentration for any mode (#/m3)
XN0MIN_SLT  = (/1. , 1., 1., 1., 1. /)


ELSE IF ( NMODE_SLT == 3) THEN

! Set the order of the loops sorted by importance
!This means that if a user choses 1 mode it will have characteristics of mode 2
!2 modes will be mode 2 & 3, whereas 3 modes will modes 1, 2 and 3
!JPSALTORDER = (/3, 2, 1, 4, 5/)
! 
   !Initial dry number median radius (um) from Vignati et al., 2001
   ! XINIRADIUS_SLT= (/0.2, 2., 12./)
   !Initial, standard deviation from Vignati et al., 2001
   ! XINISIG_SLT =  (/1.9, 2., 3./)
   !Minimum allowed number concentration for any mode (#/m3)
   ! XN0MIN_SLT  = (/1.e1 , 1. , 1.e-4 /)


!Pour 3 modes Schultz
!Initial dry number median radius (um) from Schultz et al., 2004
 XINIRADIUS_SLT= 0.5*(/0.28, 2.25, 15.32, 0., 0. /)
!Initial, standard deviation from  Schultz et al., 2004
 XINISIG_SLT =  (/1.9, 2., 2., 0., 0./)
!Minimum allowed number concentration for any mode (#/m3)
 XN0MIN_SLT  = (/1.e1 , 1. , 1.e-4, 0., 0. /)
!
END IF


END SUBROUTINE INIT_SALT
