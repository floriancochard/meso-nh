!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 surfex 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_INI_SW_SETUP
!     ##########################
!
INTERFACE
!
    SUBROUTINE INI_SW_SETUP(HRAD,KSWB_MNH,PSW_BANDS)
!
CHARACTER(LEN=*), INTENT(IN)   :: HRAD      ! type of radiation scheme
INTEGER,          INTENT(IN)   :: KSWB_MNH  ! number of SW band
REAL,             DIMENSION(:) :: PSW_BANDS ! wavelength in the middle of each band
!
END SUBROUTINE INI_SW_SETUP
!
END INTERFACE
!
END MODULE MODI_INI_SW_SETUP
!
!
!   #######################################################################
    SUBROUTINE INI_SW_SETUP(HRAD,KSWB_MNH,PSW_BANDS)
!   #######################################################################
!
!!****  *INI_SW_SETUP * - initialisation for ECMWF radiation SW bands
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
!!      modification : 01/09/03  Y. Seity, KSWB_MNH=6
!!                   02/2018 Q.Libois ECRAD
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
INTEGER,          INTENT(IN)   :: KSWB_MNH  ! number of SW band
REAL,             DIMENSION(:) :: PSW_BANDS ! wavelength in the middle of each band
!
!*       0.2   declarations of local variables
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
SELECT CASE (HRAD)
  CASE ('ECMW')

! number of SW band used in MESONH between surface and radiation

    IF (KSWB_MNH  == 1) THEN
       PSW_BANDS = 0.7E-6
    ELSEIF (KSWB_MNH  == 6) THEN
       PSW_BANDS(1) = 0.2175E-6
       PSW_BANDS(2) = 0.345E-6
       PSW_BANDS(3) = 0.565E-6
       PSW_BANDS(4) = 0.94E-6
       PSW_BANDS(5) = 1.785E-6
       PSW_BANDS(6) = 3.19E-6
    ELSE
       !callabortstop
       CALL ABORT
       STOP     
    ENDIF
    
  CASE ('ECRA') 
       PSW_BANDS(1) = 3.462E-6
       PSW_BANDS(2) = 2.788E-6
       PSW_BANDS(3) = 2.325E-6
       PSW_BANDS(4) = 2.046E-6
       PSW_BANDS(5) = 1.784E-6
       PSW_BANDS(6) = 1.462E-6
       PSW_BANDS(7) = 1.270E-6
       PSW_BANDS(8) = 1.010E-6
       PSW_BANDS(9) = 0.702E-6
       PSW_BANDS(10) = 0.533E-6
       PSW_BANDS(11) = 0.393E-6
       PSW_BANDS(12) = 0.304E-6
       PSW_BANDS(13) = 0.231E-6
       PSW_BANDS(14) = 8.021E-6       

  CASE ('FIXE','TOPA','NONE')

! number of SW band used in MESONH between surface and radiation

    IF (KSWB_MNH  == 1) THEN
       PSW_BANDS = 0.7E-6
    ELSEIF (KSWB_MNH  == 6) THEN
       PSW_BANDS(1) = 0.2175E-6
       PSW_BANDS(2) = 0.345E-6
       PSW_BANDS(3) = 0.565E-6
       PSW_BANDS(4) = 0.94E-6
       PSW_BANDS(5) = 1.785E-6
       PSW_BANDS(6) = 3.19E-6
    ELSE
!callabortstop
CALL ABORT
       STOP     
    ENDIF

!
END SELECT
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_SW_SETUP
