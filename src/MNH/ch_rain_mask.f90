!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home//MESONH/MNH-V4-6-5/src/SRC_CHIMAQ/ch_cloud_mask.f90
!-----------------------------------------------------------------
!!    #########################
      MODULE MODI_CH_RAIN_MASK
!!    #########################
!
!
INTERFACE
!
SUBROUTINE CH_RAIN_MASK(KMI,KVECNPT,PLWR,OLWR,KLWR)
IMPLICIT NONE
INTEGER, INTENT(IN)                  :: KMI
INTEGER, INTENT(IN)                  :: KVECNPT
REAL,    INTENT(IN),    DIMENSION(KVECNPT)          :: PLWR
LOGICAL, INTENT(INOUT), DIMENSION(KVECNPT)          :: OLWR
INTEGER, INTENT(INOUT)                              :: KLWR !number of points 
                                                            !passing cloud mask
END SUBROUTINE CH_RAIN_MASK
!
END INTERFACE
!
END MODULE MODI_CH_RAIN_MASK
!
!!    ####################################################
      SUBROUTINE CH_RAIN_MASK(KMI,KVECNPT,PLWR,OLWR,KLWR)
!!    ####################################################
!!
!!*** *CH_RAIN_MASK cloud mask for aqueous phase chemistry
!!
!!    PURPOSE
!!    -------
!!       The purpose of this subroutine is to select the points in the
!!    vector (computed by ch_monitorn) where the liquid water content of cloud 
!!    or rain is greater than a thresold defined by the user. 
!!
!!    METHOD
!!    ------
!!       uses a logical where LWR > thresold
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    M. Leriche    *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 4/06/07
!!
!!    EXTERNAL
!!    --------
!     none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
USE MODD_CH_MNHC_n, ONLY : XRTMIN_AQ
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)                  :: KMI
INTEGER, INTENT(IN)                  :: KVECNPT
REAL,    INTENT(IN),    DIMENSION(KVECNPT)          :: PLWR
LOGICAL, INTENT(INOUT), DIMENSION(KVECNPT)          :: OLWR
INTEGER, INTENT(INOUT)                              :: KLWR !number of points 
                                                            !passing cloud mask
!
!*      0.2    declarations of local variables
!
INTEGER :: JL
!
!-------------------------------------------------------------------------------
!
!*       1.     INITIALIZATION
!               --------------
!
OLWR(:) = .FALSE.
!
!*      2.    SELECT POINTS - RAIN MASK
!             --------------------------
!
OLWR(:) = PLWR(:)>XRTMIN_AQ
KLWR = COUNT(OLWR(:))
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CH_RAIN_MASK
