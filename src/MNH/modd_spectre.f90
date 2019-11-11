!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #################
      MODULE MODD_SPECTRE
!     #################


LOGICAL,SAVE :: LSPECTRE_U    ! computation of spectra for UT
LOGICAL,SAVE :: LSPECTRE_V    ! computation of spectra for VT
LOGICAL,SAVE :: LSPECTRE_W    ! computation of spectra for WT
LOGICAL,SAVE :: LSPECTRE_TH   ! computation of spectra for THT
LOGICAL,SAVE :: LSPECTRE_RV   ! computation of spectra for RVT
LOGICAL,SAVE :: LSPECTRE_LSU  ! computation of spectra for LSUM
LOGICAL,SAVE :: LSPECTRE_LSV  ! computation of spectra for LSVM
LOGICAL,SAVE :: LSPECTRE_LSW  ! computation of spectra for LSWM
LOGICAL,SAVE :: LSPECTRE_LSTH ! computation of spectra for LSTHM
LOGICAL,SAVE :: LSPECTRE_LSRV ! computation of spectra for LSRVM
LOGICAL,SAVE :: LSMOOTH       ! switch to smooth spectra
LOGICAL,SAVE :: LZOOM         ! switch to calculate spectra on a zoom of the domain
INTEGER,SAVE :: NITOT         ! Dimensions in x direction
INTEGER,SAVE :: NJTOT         ! Dimensions in y direction
INTEGER,SAVE :: NXDEB         ! horizontal position (i) of the ORigin of the zoom
INTEGER,SAVE :: NYDEB         ! horizontal position (j) of the ORigin of the zoom
CHARACTER(LEN=6):: CTYPEFILE  ! type of the INPUT FM-file
LOGICAL,SAVE :: LSTAT         ! flag to have some statistics 
END MODULE MODD_SPECTRE
