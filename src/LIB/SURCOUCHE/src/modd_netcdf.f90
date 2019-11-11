!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------
MODULE MODD_NETCDF
IMPLICIT NONE 

INTEGER,PARAMETER :: IDCDF_KIND = SELECTED_INT_KIND(8)

TYPE IOCDF
   TYPE(DIMCDF), POINTER :: DIM_NI      => NULL()
   TYPE(DIMCDF), POINTER :: DIM_NJ      => NULL()
   TYPE(DIMCDF), POINTER :: DIM_LEVEL   => NULL()
   TYPE(DIMCDF), POINTER :: DIM_NI_U    => NULL()
   TYPE(DIMCDF), POINTER :: DIM_NJ_U    => NULL()
   TYPE(DIMCDF), POINTER :: DIM_NI_V    => NULL()
   TYPE(DIMCDF), POINTER :: DIM_NJ_V    => NULL()
   TYPE(DIMCDF), POINTER :: DIM_LEVEL_W => NULL()
   TYPE(DIMCDF), POINTER :: DIMTIME     => NULL()
   TYPE(DIMCDF), POINTER :: DIMSTR      => NULL()
   TYPE(DIMCDF), POINTER :: DIMLIST     => NULL()
END TYPE IOCDF

TYPE DIMCDF
   CHARACTER(LEN=32)        :: NAME = ''
   INTEGER(KIND=IDCDF_KIND) :: LEN  = 0
   INTEGER(KIND=IDCDF_KIND) :: ID   = -1
   TYPE(DIMCDF), POINTER    :: NEXT => NULL()
END TYPE DIMCDF

TYPE TPTR2DIMCDF
  TYPE(DIMCDF),POINTER :: TDIM
END TYPE TPTR2DIMCDF

TYPE(DIMCDF),TARGET :: TDIM_DUMMY !Dummy dimension

END MODULE MODD_NETCDF
