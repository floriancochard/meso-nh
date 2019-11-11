!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 profiler 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #######################
      SUBROUTINE INI_PROFILER_n
!     #######################
!
!
!!****  *INI_PROFILER_n* - user initializes the station location
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
!!   Must be defined (for each aircraft):
!!   ---------------
!!
!!  No default exist for these variables.
!!  ************************************
!!
!!  1) Number of stations
!!  2) the model in which these stations are
!!     if NOT initialized, the stations are NOT used.
!!
!!  3) the (LAT, LON, ALT) latitude,longitude and altitude of the station location.
!!  4) the station name
!!
!!
!!
!!   Can be defined  (for each aircraft):
!!   --------------
!!
!!
!!  9) the time step for data storage.
!!    default is 60s
!!
!! 10) the name or title describing the balloon (8 characters)
!!     default is the balloon type (6 characters) + the balloon numbers (2 characters)
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
!!      Pierre Tulet             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original 15/01/2002
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_PROFILER_n
USE MODD_PARAMETERS
!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
!----------------------------------------------------------------------------
!
!*      1.   Nameliste 
!            ---------
NUMBPROFILER             = 0
!
IF (NUMBPROFILER > 0) THEN
ALLOCATE(TPROFILER%LAT (NUMBPROFILER))
ALLOCATE(TPROFILER%LON (NUMBPROFILER))
ALLOCATE(TPROFILER%ALT (NUMBPROFILER))
ALLOCATE(TPROFILER%NAME(NUMBPROFILER))
ALLOCATE(TPROFILER%TYPE(NUMBPROFILER))
!
TPROFILER%LON = XUNDEF
TPROFILER%LAT = XUNDEF
TPROFILER%ALT = XUNDEF
TPROFILER%NAME = "    "
TPROFILER%TYPE = "         "
!
TPROFILER%STEP = 900.
!
!* location (latitude, longitude, altitude)
!
!
TPROFILER%LAT                = (/ 43.3000, 43.3300, 43.6200, 43.6550, 43.3400, &
                                  43.3000, 43.3500, 43.8128, 44.1711, 44.1689, &
                                  44.0833, 43.6200, 43.7164, 43.5333, 43.4833, &
                                  43.4800, 43.4856,                            &
                                  43.5423, 43.3000, 43.5000, 43.2568, 43.5000, &
                                  43.3333, 43.9070, 43.5430, 43.5300, 43.5000, &
                                  43.4300,                                     &
                                  43.2780, 43.3277, 43.3396, 43.3230, 43.3559, &
                                  43.3060, 43.3724, 43.2568, 43.3834, 43.3417, &
                                  43.3165, 43.3873, 43.2833, 43.2900, 43.2842, &
                                  43.3597, 43.2788, 43.4999 /)
!
TPROFILER%LON                = (/ 5.3790, 4.8200, 5.2000, 6.0986, 5.4100, &
                                  5.3833, 5.4000, 5.1506, 5.2875, 5.0692, &
                                  5.0500, 5.4000, 5.7653, 5.0667, 5.2833, &
                                  5.3200, 5.3405,                         &
                                  4.9010, 5.4000, 4.9300, 5.4054, 5.3700, &
                                  5.4167, 4.8988, 5.0740, 5.0700, 5.3700, &
                                  5.2300,                                 &
                                  5.5762, 5.4782, 5.5873, 5.3668, 5.4550, &
                                  5.3951, 5.4842, 5.4054, 5.4253, 5.4113, &
                                  5.4369, 5.3893, 5.5114, 5.3788, 5.4611, &
                                  5.3994, 5.3538, 5.3674 /)
TPROFILER%ALT                = (/ 2.,2.,2.,2.,2.,&
                                  2.,2.,2.,2.,2.,&
                                  2.,2.,2.,2.,2.,&
                                  2.,2.,         &
                                  2.,2.,2.,2.,2.,&
                                  2.,2.,2.,2.,2.,&
                                  2.,            &
                                  2.,2.,2.,2.,2.,&
                                  2.,2.,2.,2.,2.,&
                                  2.,2.,2.,2.,2.,&
                                  2.,2.,2.       /)                          
!
TPROFILER%NAME               = (/ 'CAAM    ', 'CRAU    ', 'BARD    ', 'MONT    ', 'IUTF    ', &
                                  'OBSF    ', 'VALF    ', 'LUBE    ', 'VENT    ', 'DENT    ', &
                                  'CARP    ', 'DUPA    ', 'VINO    ', 'CHAM    ', 'REA1    ', &
                                  'REA2    ', 'REA3    ',                                     &
                                  'SOD_4M  ', 'UHF_4M  ', 'VHF_STM ', 'SOD_CO  ', 'UHF_DEG ', &
                                  'SOD_EC  ', 'SOD_E1  ', 'SOD_E2  ', 'UHF_ED  ', 'VHF_LS  ', &
                                  'UHF_DEP ',                                                 &
                                  'AGRI    ', 'ALLA    ', 'BOYE    ', 'CANE    ', 'CHAT    ', &
                                  'CINQ    ', 'COTE    ', 'DR12    ', 'ETOI    ', 'IUT1    ', &
                                  'JULI    ', 'ONYX    ', 'PENN    ', 'TRIB    ', 'VALB    ', &
                                  'VALL    ', 'MARS    ', 'AIRB    ' /)                                               
!
TPROFILER%TYPE               = (/ 'flux    ', 'flux    ', 'flux    ', 'flux    ', 'flux    ', &
                                  'flux    ', 'flux    ', 'flux    ', 'flux    ', 'flux    ', &
                                  'flux    ', 'flux    ', 'flux    ', 'flux    ', 'flux    ', &
                                  'flux    ', 'flux    ',                                     &
                                  'sodar   ', 'uhf     ', 'vhf     ', 'sodar   ', 'uhf     ', &
                                  'sodar   ', 'sodar   ', 'sodar   ', 'uhf     ', 'vhf     ', &
                                  'uhf     ',                                                 &
                                  'gps     ', 'gps     ', 'gps     ', 'gps     ', 'gps     ', &
                                  'gps     ', 'gps     ', 'gps     ', 'gps     ', 'gps     ', &
                                  'gps     ', 'gps     ', 'gps     ', 'gps     ', 'gps     ', &
                                  'gps     ', 'gps     ', 'gps     ' /)                     
!
!----------------------------------------------------------------------------
ENDIF
!
END SUBROUTINE INI_PROFILER_n
