!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!   ########################
     MODULE MODI_CH_AER_EQM_CORMASS
!!   ########################
!!
INTERFACE
!!
SUBROUTINE CH_AER_EQM_CORMASS(PSVT)
IMPLICIT NONE
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)   :: PSVT
END SUBROUTINE CH_AER_EQM_CORMASS
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_EQM_CORMASS
!!
!!
!!   ############################################################
     SUBROUTINE CH_AER_EQM_CORMASS(PSVT) 
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!    Realise la conservation de la masse 
!!    Filtre les valeurs des moments 0 et 6 inferieures aux valeurs 
!!    minimales (Rg et SIG) introduites en nameliste 
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!    none
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CST, ONLY :  XMNH_TINY   
!!
IMPLICIT NONE
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)   :: PSVT
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
!
!*      0.2    declarations local variables
!
!-------------------------------------------------------------------------------

!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               ------------------------------------
!
  PSVT(:,:,:,:) = MAX(PSVT(:,:,:,:),XMNH_TINY)

!
END SUBROUTINE CH_AER_EQM_CORMASS
