!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 balloon 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################
      SUBROUTINE INI_BALLOON
!     ######################
!
!
!!****  *INI_BALLOON* - user initializes the balloon characteristics
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
!!    For constant volume Balloon, horizontal advection using horizontal wind
!!        vertical spped of the balloon calculated using the balloon equation
!!        (Koffi et AL 2000, JAS vol 57 P.2007-2021)
!!
!!   Must be defined (for each balloon):
!!   ---------------
!!
!!  No default exist for these variables.
!!  ************************************
!!
!!  1) the model in which the balloon will evolve
!!     if NOT initialized, the balloon is NOT used.
!!  1.1) the possibility to switch from a model to its dad or kid
!!       'FIX' : NMODEL used during the run
!!       'MOB' : best resolution model used. NMODEL=1 is used at the beginning
!!
!!  2) the type of balloon
!!
!!     'RADIOS' for radiosounding balloon
!!     'ISODEN' for iso-density balloon
!!     'CVBALL' for constant volume Balloon
!!
!!  3) the launching date and time
!!
!!  4) the latitude of the launching site
!!
!!  5) the longitude of the launching site
!!
!!  6) the altitude of the launching site (for 'RADIOS')
!!
!!                      OR
!!
!!     the altitude OR pressure of balloon at start of the leveled flight
!!     (for 'ISODEN'). In this case, the density of this level will be computed,
!!     and the balloon will evolve at this density level.
!!
!!
!!
!!   Can be defined  (for each balloon):
!!   --------------
!!
!!  7) the ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!!     default is 5m/s
!!
!!  8) the time step for data storage.
!!    default is 60s
!!
!!  9) the name or title describing the balloon (8 characters)
!!     default is the balloon type (6 characters) + the balloon numbers (2 characters)
!!
!!  10) for 'CVBALL' the aerodynamic drag coefficient of the balloon
!!
!!  11) for 'CVBALL' the induced drag coefficient (i.e. air shifted by the balloon)
!!
!!  12) for 'CVBALL' the volume of the balloon
!!
!!  13) for 'CVBALL' the mass of the balloon
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
!!     Original 15/05/2000
!!              Apr,19, 2001 (G.Jaubert) add CVBALL type and switch in models
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_AIRCRAFT_BALLOON
USE MODD_CST
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
!*      1.   Balloon number 1
!            ----------------
!
!* model number
!
TBALLOON1%NMODEL             = 0
TBALLOON1%MODEL             = 'MOB'
!
!* balloon type
!
TBALLOON1%TYPE               = 'CVBALL'
!
!* balloon name
!
TBALLOON1%TITLE              = 'CVB1MOBI'
!
!* launching date and time
!
TBALLOON1%LAUNCH%TDATE%YEAR  =  1999
TBALLOON1%LAUNCH%TDATE%MONTH =    09
TBALLOON1%LAUNCH%TDATE%DAY   =    19
TBALLOON1%LAUNCH%TIME        = 32460.
!
!* latitude and longitude of launching site (decimal degree)
!
TBALLOON1%LAT                = 45.800
TBALLOON1%LON                =  8.629
!
!* altitude of the launching site for 'RADIOS'
!* altitude or pressure of the flight level for 'ISODEN'
!
!TBALLOON1%ALT                =   3959.
TBALLOON1%PRES               = 98450.
!
!* time step for data storage  (s)
!
TBALLOON1%STEP               = 20.
!
!* ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!
TBALLOON1%WASCENT            = 0.
!
!* aerodynamic drag coefficient of the balloon (for 'CVBALL')
!* induced drag coefficient (i.e. air shifted by the balloon) (for 'CVBALL')
!* volume of the balloon (m3) (for 'CVBALL')
!* mass of the balloon (kg) (for 'CVBALL')
!
TBALLOON1%AERODRAG           = 0.44
TBALLOON1%INDDRAG           = 0.014
TBALLOON1%VOLUME           = 3.040
TBALLOON1%MASS           = 2.4516
TBALLOON1%DIAMETER           = ((3.*TBALLOON1%VOLUME)/(4.*XPI))**(1./3.)
!
!----------------------------------------------------------------------------
!
!*      2.   Balloon number 2
!            ----------------
!
!* model number
!
TBALLOON2%NMODEL             = 0
TBALLOON2%MODEL             = 'MOB'
!
!* balloon type
!
TBALLOON2%TYPE               = 'CVBALL'
!
!* balloon name
!
TBALLOON2%TITLE              = 'CVB2MOBI'
!
!* launching date and time
!
TBALLOON2%LAUNCH%TDATE%YEAR  =  1999
TBALLOON2%LAUNCH%TDATE%MONTH =    09
TBALLOON2%LAUNCH%TDATE%DAY   =    19
TBALLOON2%LAUNCH%TIME        = 39660.
!
!* latitude and longitude of launching site (decimal degree)
!
TBALLOON2%LAT                = 45.800
TBALLOON2%LON                =  8.630
!
!* altitude of the launching site for 'RADIOS'
!* altitude or pressure of the flight level for 'ISODEN'
!
!TBALLOON2%ALT                =   3959.
TBALLOON2%PRES               = 98490.
!
!* time step for data storage  (s)
!
TBALLOON2%STEP               = 20.
!
!* ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!
TBALLOON2%WASCENT            = 0.
!
!* aerodynamic drag coefficient of the balloon (for 'CVBALL')
!* induced drag coefficient (i.e. air shifted by the balloon) (for 'CVBALL')
!* volume of the balloon (m3) (for 'CVBALL')
!* mass of the balloon (kg) (for 'CVBALL')
!
TBALLOON2%AERODRAG           = 0.44
TBALLOON2%INDDRAG           = 0.014
TBALLOON2%VOLUME           = 3.040
TBALLOON2%MASS           = 2.58087
TBALLOON2%DIAMETER           = ((3.*TBALLOON2%VOLUME)/(4.*XPI))**(1./3.)
!
!-------------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
!*      3.   Balloon number 3
!            ----------------
!
!* model number
!
TBALLOON3%NMODEL             = 0
TBALLOON3%MODEL             = 'MOB'
!
!* balloon type
!
TBALLOON3%TYPE               = 'RADIOS'
!
!* balloon name
!
TBALLOON3%TITLE              = 'RSMASE19'
!
!* launching date and time
!
TBALLOON3%LAUNCH%TDATE%YEAR  =  1999
TBALLOON3%LAUNCH%TDATE%MONTH =    09
TBALLOON3%LAUNCH%TDATE%DAY   =    19
TBALLOON3%LAUNCH%TIME        = 68400.
!
!* latitude and longitude of launching site (decimal degree)
!
TBALLOON3%LAT                = 46.810
TBALLOON3%LON                =  9.39 
!
!* altitude of the launching site for 'RADIOS'
!* altitude or pressure of the flight level for 'ISODEN'
!
TBALLOON3%ALT                =   865. 
!TBALLOON3%PRES               = 62360.
!
!* time step for data storage  (s)
!
TBALLOON3%STEP               = 20.
!
!* ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!
TBALLOON3%WASCENT            = 4.85
!
!* aerodynamic drag coefficient of the balloon (for 'CVBALL')
!* induced drag coefficient (i.e. air shifted by the balloon) (for 'CVBALL')
!* volume of the balloon (m3) (for 'CVBALL')
!* mass of the balloon (kg) (for 'CVBALL')
!
TBALLOON3%AERODRAG           = 0.44
TBALLOON3%INDDRAG           = 0.014
TBALLOON3%VOLUME           = 3.040
TBALLOON3%MASS           = 2.4516
TBALLOON3%DIAMETER           = ((3.*TBALLOON3%VOLUME)/(4.*XPI))**(1./3.)
!
!
!----------------------------------------------------------------------------
!
!*      4.   Balloon number 4
!            ----------------
!
!* model number
!
TBALLOON4%NMODEL             = 0
TBALLOON4%MODEL             = 'FIX'
!
!* balloon type
!
TBALLOON4%TYPE               = 'CVBALL'
!
!* balloon name
!
TBALLOON4%TITLE              = 'CVB1ACVB'
!
!* launching date and time
!
TBALLOON4%LAUNCH%TDATE%YEAR  =  1999
TBALLOON4%LAUNCH%TDATE%MONTH =    09
TBALLOON4%LAUNCH%TDATE%DAY   =    19
TBALLOON4%LAUNCH%TIME        = 32460.
!
!* latitude and longitude of launching site (decimal degree)
!
TBALLOON4%LAT                = 45.922
TBALLOON4%LON                =  8.646
!
!* altitude of the launching site for 'RADIOS'
!* altitude or pressure of the flight level for 'ISODEN'
!
TBALLOON4%ALT                =   3959.
!TBALLOON4%PRES               = 62360.
!
!* time step for data storage  (s)
!
TBALLOON4%STEP               = 20.
!
!* ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!
!TBALLOON4%WASCENT            = 2.55
!
!* aerodynamic drag coefficient of the balloon (for 'CVBALL')
!* induced drag coefficient (i.e. air shifted by the balloon) (for 'CVBALL')
!* volume of the balloon (m3) (for 'CVBALL')
!* mass of the balloon (kg) (for 'CVBALL')
!
TBALLOON4%AERODRAG           = 0.44
TBALLOON4%INDDRAG           = 0.014
TBALLOON4%VOLUME           = 3.040
TBALLOON4%MASS           = 2.4516
TBALLOON4%DIAMETER           = ((3.*TBALLOON4%VOLUME)/(4.*XPI))**(1./3.)
!
!----------------------------------------------------------------------------
!
!*      5.   Balloon number 5
!            ----------------
!
!* model number
!
TBALLOON5%NMODEL             = 0
TBALLOON5%MODEL             = 'FIX'
!
!* balloon type
!
TBALLOON5%TYPE               = 'CVBALL'
!
!* balloon name
!
TBALLOON5%TITLE              = 'CVB1DEPA'
!
!* launching date and time
!
TBALLOON5%LAUNCH%TDATE%YEAR  =  1999
TBALLOON5%LAUNCH%TDATE%MONTH =    09
TBALLOON5%LAUNCH%TDATE%DAY   =    19
TBALLOON5%LAUNCH%TIME        = 32435.
!
!* latitude and longitude of launching site (decimal degree)
!
TBALLOON5%LAT                = 45.800
TBALLOON5%LON                =  8.630
!
!* altitude of the launching site for 'RADIOS'
!* altitude or pressure of the flight level for 'ISODEN'
!
TBALLOON5%ALT                =    340.
!TBALLOON5%PRES               = 62360.
!
!* time step for data storage  (s)
!
TBALLOON5%STEP               = 20.
!
!* ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!
!TBALLOON5%WASCENT            = 2.55
!
!* aerodynamic drag coefficient of the balloon (for 'CVBALL')
!* induced drag coefficient (i.e. air shifted by the balloon) (for 'CVBALL')
!* volume of the balloon (m3) (for 'CVBALL')
!* mass of the balloon (kg) (for 'CVBALL')
!
TBALLOON5%AERODRAG           = 0.44
TBALLOON5%INDDRAG           = 0.014
TBALLOON5%VOLUME           = 3.040
TBALLOON5%MASS           = 2.4516
TBALLOON5%DIAMETER           = ((3.*TBALLOON5%VOLUME)/(4.*XPI))**(1./3.)
!
!----------------------------------------------------------------------------
!
!*      6.   Balloon number 6
!            ----------------
!
!* model number
!
TBALLOON6%NMODEL             = 0
TBALLOON6%MODEL             = 'FIX'
!
!* balloon type
!
TBALLOON6%TYPE               = 'CVBALL'
!
!* balloon name
!
TBALLOON6%TITLE              = 'CVB1RCVB'
!
!* launching date and time
!
TBALLOON6%LAUNCH%TDATE%YEAR  =  1999
TBALLOON6%LAUNCH%TDATE%MONTH =    09
TBALLOON6%LAUNCH%TDATE%DAY   =    19
TBALLOON6%LAUNCH%TIME        = 32460.
!
!* latitude and longitude of launching site (decimal degree)
!
TBALLOON6%LAT                = 45.922
TBALLOON6%LON                =  8.646
!
!* altitude of the launching site for 'RADIOS'
!* altitude or pressure of the flight level for 'ISODEN'
!
!TBALLOON6%ALT                =   3959.
!TBALLOON6%PRES               = 62360.
!
!* time step for data storage  (s)
!
TBALLOON6%STEP               = 20.
!
!* ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!
!TBALLOON6%WASCENT            = 2.55
!
!* aerodynamic drag coefficient of the balloon (for 'CVBALL')
!* induced drag coefficient (i.e. air shifted by the balloon) (for 'CVBALL')
!* volume of the balloon (m3) (for 'CVBALL')
!* mass of the balloon (kg) (for 'CVBALL')
!
TBALLOON6%AERODRAG           = 0.44
TBALLOON6%INDDRAG           = 0.014
TBALLOON6%VOLUME           = 3.040
TBALLOON6%MASS           = 2.4516
TBALLOON6%DIAMETER           = ((3.*TBALLOON6%VOLUME)/(4.*XPI))**(1./3.)
!
!----------------------------------------------------------------------------
!
!*      7.   Balloon number 7
!            ----------------
!
!* model number
!
TBALLOON7%NMODEL             = 0
TBALLOON7%MODEL             = 'FIX'
!
!* balloon type
!
TBALLOON7%TYPE               = 'CVBALL'
!
!* balloon name
!
TBALLOON7%TITLE              = 'CVB1PISO'
!
!* launching date and time
!
TBALLOON7%LAUNCH%TDATE%YEAR  =  1999
TBALLOON7%LAUNCH%TDATE%MONTH =    09
TBALLOON7%LAUNCH%TDATE%DAY   =    19
TBALLOON7%LAUNCH%TIME        = 32460.
!
!* latitude and longitude of launching site (decimal degree)
!
TBALLOON7%LAT                = 45.922
TBALLOON7%LON                =  8.646
!
!* altitude of the launching site for 'RADIOS'
!* altitude or pressure of the flight level for 'ISODEN'
!
!TBALLOON7%ALT                =   3959.
TBALLOON7%PRES               = 62360.
!
!* time step for data storage  (s)
!
TBALLOON7%STEP               = 20.
!
!* ascentional vertical speed of the ballon (in calm air) (for 'RADIOS')
!
!TBALLOON7%WASCENT            = 2.55
!
!* aerodynamic drag coefficient of the balloon (for 'CVBALL')
!* induced drag coefficient (i.e. air shifted by the balloon) (for 'CVBALL')
!* volume of the balloon (m3) (for 'CVBALL')
!* mass of the balloon (kg) (for 'CVBALL')
!
TBALLOON7%AERODRAG           = 0.44
TBALLOON7%INDDRAG           = 0.014
TBALLOON7%VOLUME           = 3.040
TBALLOON7%MASS           = 2.4516
TBALLOON7%DIAMETER           = ((3.*TBALLOON7%VOLUME)/(4.*XPI))**(1./3.)
!
!----------------------------------------------------------------------------
!
END SUBROUTINE INI_BALLOON
