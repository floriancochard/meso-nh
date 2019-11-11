!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
  SUBROUTINE GET_NB_PROCIO_WRITE_MNH( KNB_PROCIO, KRESP )
!
!!****  *GET_NB_PROCIO_WRITE_MNH* - gets the number of processes used for Output of file MODD_IO_SURF_MNH::COUTFILE
!!                        
!!
!!    PURPOSE
!!    -------
!!      call GET_NB_PROCIO_WRITE_MNH from SURFEX to get the number of processes used 
!!      for Output of file MODD_IO_SURF_MNH::COUTFILE in MESO-NH (defined by user in namelist)
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
!!
!!    AUTHOR
!!    ------
!!	M. Moge   *LA - UPS*  08/01/2016
!!      J. escobar 19/04/2016 : bypass , For pb IOZ/NETCDF , pretende alway 2 ( > 1 ) I/O processors for homogenus PGD files	
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!
IMPLICIT NONE
!
!*      0.    DECLARATIONS
!             ------------
!
!*      0.1   Declarations of arguments
!
INTEGER, INTENT(OUT) :: KNB_PROCIO ! number of processes used for IO
INTEGER, INTENT(OUT) :: KRESP      ! return-code
!
!*      0.2   Declarations of local variables
!
!NONE
!----------------------------------------------------------------
KNB_PROCIO = 2
KRESP = 0
!
  END SUBROUTINE GET_NB_PROCIO_WRITE_MNH
