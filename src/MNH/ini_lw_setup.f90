!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_INI_LW_SETUP
!     ##########################
!
INTERFACE
!
    SUBROUTINE INI_LW_SETUP(HRAD,KLWB_MNH,PLW_BANDS)
!
CHARACTER(LEN=*), INTENT(IN)   :: HRAD      ! type of radiation scheme
INTEGER,          INTENT(IN)   :: KLWB_MNH  ! number of SW band
REAL,             DIMENSION(:) :: PLW_BANDS ! wavelength in the middle of each band
!
END SUBROUTINE INI_LW_SETUP
!
END INTERFACE
!
END MODULE MODI_INI_LW_SETUP
!
!
!   #######################################################################
    SUBROUTINE INI_LW_SETUP(HRAD,KLWB_MNH,PLW_BANDS)
!   #######################################################################
!
!!****  *INI_LW_SETUP * - initialisation for RRTM radiation LW bands
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!! 
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/03/03
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=*), INTENT(IN)   :: HRAD      ! type of radiation scheme
INTEGER,          INTENT(IN)   :: KLWB_MNH  ! number of SW band
REAL,             DIMENSION(:) :: PLW_BANDS ! wavelength in the middle of each band
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
SELECT CASE (HRAD)
  CASE ('ECMW','ECRA')

! number of LW band used in MESONH between surface and radiation
! beware of spectral match with RRTM bands which might not be monotically increasing

       PLW_BANDS(1) = 520E-6
       PLW_BANDS(2) = 30E-6
       PLW_BANDS(3) = 17.9E-6
       PLW_BANDS(4) = 15.1E-6
       PLW_BANDS(5) = 13.2E-6
       PLW_BANDS(6) = 11.2E-6
       PLW_BANDS(7) = 9.73E-6
       PLW_BANDS(8) = 8.87E-6
       PLW_BANDS(9) = 7.83E-6
       PLW_BANDS(10) = 6.98E-6
       PLW_BANDS(11) = 6.16E-6
       PLW_BANDS(12) = 5.18E-6
       PLW_BANDS(13) = 4.62E-6
       PLW_BANDS(14) = 4.32E-6
       PLW_BANDS(15) = 4.02E-6      
       PLW_BANDS(16) = 3.59E-6
  CASE ('FIXE','TOPA','NONE')

! number of SW band used in MESONH between surface and radiation

    IF (KLWB_MNH  == 1) THEN
       PLW_BANDS = 12E-6
    ELSEIF (KLWB_MNH  == 16) THEN
       PLW_BANDS(1) = 520E-6
       PLW_BANDS(2) = 30E-6
       PLW_BANDS(3) = 17.9E-6
       PLW_BANDS(4) = 15.1E-6
       PLW_BANDS(5) = 13.2E-6
       PLW_BANDS(6) = 11.2E-6
       PLW_BANDS(7) = 9.73E-6
       PLW_BANDS(8) = 8.87E-6
       PLW_BANDS(9) = 7.83E-6
       PLW_BANDS(10) = 6.98E-6
       PLW_BANDS(11) = 6.16E-6
       PLW_BANDS(12) = 5.18E-6
       PLW_BANDS(13) = 4.62E-6
       PLW_BANDS(14) = 4.32E-6
       PLW_BANDS(15) = 4.02E-6      
       PLW_BANDS(16) = 3.59E-6
    ELSE
!callabortstop
       CALL ABORT
       STOP     
    ENDIF

!
END SELECT
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_LW_SETUP
