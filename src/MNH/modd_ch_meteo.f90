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
!!    ###################### 
      MODULE MODD_CH_METEO
!!    ######################
!!
!!*** *MODD_CH_METEO*
!!
!!    PURPOSE
!!    -------
!       contains the meteovariables that will be updated by ch_update_meteo
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 20/04/99
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
SAVE
!
INTEGER                           :: NMETEORECS 
		    ! number of records
!
INTEGER                           :: NMETEORECACT = 1
		    ! actual record (used in temporal interpolation)
!
REAL, DIMENSION(:),   ALLOCATABLE :: XMETEOTIME
                    ! the time of the individual records
!
REAL, DIMENSION(:,:), ALLOCATABLE :: XMETEODATA 
		    ! the meteodata (the first dimension is the number of
		    ! elements NMETEOVARS, the second dimension is the number
		    ! of records)
!
END MODULE MODD_CH_METEO
