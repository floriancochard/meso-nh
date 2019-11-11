!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 profiler 2006/10/24 10:07:27
!-----------------------------------------------------------------
!      #########################
MODULE MODI_END_DIAG_IN_RUN
!      #########################
!
INTERFACE
!
SUBROUTINE END_DIAG_IN_RUN
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE END_DIAG_IN_RUN
!
END INTERFACE
!
END MODULE MODI_END_DIAG_IN_RUN
!
!     ####################
SUBROUTINE END_DIAG_IN_RUN
!     ####################
!
!
!!****  *END_DIAG_IN_RUN* - 
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Valery Masson             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original 11/2003
!!
!!                   02/2018 Q.Libois ECRAD
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_PARAMETERS, ONLY : XUNDEF
USE MODD_DIAG_IN_RUN
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!-------------------------------------------------------------------------------
!
DEALLOCATE(XCURRENT_RN    )! net radiation
DEALLOCATE(XCURRENT_H     )! sensible heat flux
DEALLOCATE(XCURRENT_LE    )! latent heat flux
DEALLOCATE(XCURRENT_LEI   )! Solid latent heat flux
DEALLOCATE(XCURRENT_GFLUX )! ground flux
DEALLOCATE(XCURRENT_LWD   )! incoming longwave at the surface
DEALLOCATE(XCURRENT_LWU   )! outcoming longwave at the surface
DEALLOCATE(XCURRENT_SWD   )! incoming Shortwave at the surface
DEALLOCATE(XCURRENT_SWU   )! outcoming Shortwave at the surface
IF(ALLOCATED(XCURRENT_SWDIR)) DEALLOCATE(XCURRENT_SWDIR )! incoming Shortwave direct at the surface
IF(ALLOCATED(XCURRENT_SWDIFF))DEALLOCATE(XCURRENT_SWDIFF)! incoming Shortwave diffuse at the surface
DEALLOCATE(XCURRENT_T2M   )! temperature at 2m
DEALLOCATE(XCURRENT_Q2M   )! humidity at 2m
DEALLOCATE(XCURRENT_HU2M  )! humidity at 2m
DEALLOCATE(XCURRENT_ZON10M)! zonal wind at 10m
DEALLOCATE(XCURRENT_MER10M)! meridian wind at 10m
DEALLOCATE(XCURRENT_DSTAOD)! dust aerosol optical depth
DEALLOCATE(XCURRENT_SFCO2   ) ! CO2 Surface flux
DEALLOCATE(XCURRENT_TKE_DISS) ! Tke dissipation rate
DEALLOCATE(XCURRENT_SLTAOD)   ! Salt aerosol optical depth
DEALLOCATE(XCURRENT_ZWS   )   ! Significant height of waves

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE END_DIAG_IN_RUN
