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
      SUBROUTINE INI_STATION_n
!     #######################
!
!
!!****  *INI_STATION_n* - user initializes the station location
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
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_STATION_n
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
NUMBSTAT             = 0
!
IF (NUMBSTAT > 0) THEN
ALLOCATE  (TSTATION%LAT(NUMBSTAT))
ALLOCATE  (TSTATION%LON(NUMBSTAT))
ALLOCATE  (TSTATION%I(NUMBSTAT))
ALLOCATE  (TSTATION%J(NUMBSTAT))
ALLOCATE  (TSTATION%Z(NUMBSTAT))
ALLOCATE  (TSTATION%K(NUMBSTAT))
ALLOCATE  (TSTATION%NAME(NUMBSTAT))
ALLOCATE  (TSTATION%TYPE(NUMBSTAT))
!
TSTATION%LON  = XUNDEF
TSTATION%LAT  = XUNDEF
TSTATION%Z    = XUNDEF
TSTATION%K    = XUNDEF
TSTATION%I    = XUNDEF
TSTATION%J    = XUNDEF
TSTATION%NAME = "        "
TSTATION%TYPE = "        "
!
TSTATION%STEP = 10.
!
!* location (latitude, longitude, altitude)
!
!***************************************************************
! * Horizontal location
! You have to choose between (TSTATION%LAT,TSTATION%LON) 
! or  (TSTATION%I,TSTATION%J) for all the stations 
! if both are defined it will choose (TSTATION%LAT,TSTATION%LON)
!***************************************************************
!
!TSTATION%LAT              = (/ 45.0  /) 
!TSTATION%LON              = (/ 4.5 /)
TSTATION%I                = (/ 25  /) 
TSTATION%J                = (/ 20 /)
!
!***************************************************************
! * Vertical location
! You have to choose between TSTATION%K and TSTATION%Z
! for all the stations
! if both are defined it will choose TSTATION%K 
!***************************************************************
!TSTATION%Z                  = (/ 10., 500. /) 
!
TSTATION%K   = (/ 10 /)              
!
!***************************************************************
!* station name
!***************************************************************
TSTATION%NAME        = (/ 'BIDON'  /)
!***************************************************************
!* station type
!***************************************************************
TSTATION%TYPE        = (/ 'sol     '/)
!
!----------------------------------------------------------------------------
ENDIF
!
END SUBROUTINE INI_STATION_n
