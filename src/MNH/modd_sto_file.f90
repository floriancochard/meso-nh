!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------
!     ######spl
MODULE MODD_STO_FILE
!################
CHARACTER (LEN=28), SAVE :: CFILES(100)       ! names of the files to be treated
INTEGER           , SAVE :: NSTART_SUPP(100)  ! supplementary starts 
                                              ! for the lagrangian trajectories 
!
END MODULE MODD_STO_FILE
