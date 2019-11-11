!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!    ########################### 
      MODULE MODI_CH_UPDATE_METEO
!!    ########################### 
!
INTERFACE
SUBROUTINE CH_UPDATE_METEO(TPM,PTIME)
USE MODD_CH_M9_n, ONLY: METEOTRANSTYPE
IMPLICIT NONE
TYPE(METEOTRANSTYPE), INTENT(INOUT) :: TPM
REAL,                 INTENT(IN)    :: PTIME ! current simulation time
END SUBROUTINE CH_UPDATE_METEO
END INTERFACE
END MODULE MODI_CH_UPDATE_METEO
!!
!!    ##################################### 
      SUBROUTINE CH_UPDATE_METEO(TPM,PTIME)
!!    #####################################
!!
!!***  *CH_UPDATE_METEO*
!!
!!    PURPOSE
!!    -------
!!    update a set of meteo variables
!!
!!**  METHOD
!!    ------
!!    interpolate NMETEOVARS variables in time
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 21/04/95
!!    27/07/96 (K. Suhre) restructured
!!    20/04/99 (K. Suhre) meteo variables are now interpolated in time
!!    01/12/03 (D. Gazen)   change Chemical scheme interface!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9_n,    ONLY: NMETEOVARS, METEOTRANSTYPE
USE MODD_CH_METEO,   ONLY: NMETEORECS, XMETEOTIME, XMETEODATA, NMETEORECACT
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
TYPE(METEOTRANSTYPE), INTENT(INOUT) :: TPM  ! the meteo variables
REAL,                 INTENT(IN)    :: PTIME ! current simulation time
!
!*       0.2  declaration of local variables
!     ----------------
REAL :: ZALPHA ! interpolation weight
!
!------------------------------------------------------------------------------
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
IF (PTIME .LE. XMETEOTIME(1)) THEN
!
! take first record
  TPM%XMETEOVAR(:) = XMETEODATA(:,1) 
!
ELSE IF (PTIME .GE. XMETEOTIME(NMETEORECS)) THEN
!
! take last record
  TPM%XMETEOVAR(1:NMETEOVARS) = XMETEODATA(1:NMETEOVARS,NMETEORECS) 
!
ELSE
!
! interpolate meteo variables in time
  IF (PTIME .GE. XMETEOTIME(NMETEORECACT+1)) NMETEORECACT = NMETEORECACT+1
!
  ZALPHA = (PTIME                      - XMETEOTIME(NMETEORECACT)) &
	 / (XMETEOTIME(NMETEORECACT+1) - XMETEOTIME(NMETEORECACT))
!
  TPM%XMETEOVAR(1:NMETEOVARS) = ZALPHA * XMETEODATA(1:NMETEOVARS,NMETEORECACT+1) &
  	           + (1.-ZALPHA) * XMETEODATA(1:NMETEOVARS,NMETEORECACT)
!
END IF
!
END SUBROUTINE CH_UPDATE_METEO
