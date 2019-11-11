!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
MODULE modd_wp
  !! necessary for compiling Gelato, but not included in LIB/GELATO for reasons relevant to 
  !! Gelato build process
  !! Modi J.Escobar 24/03/2017 : auto detection of size of real for possible compilation in 4/8bytes
   REAL , PRIVATE             :: REAL_DEF_WP
   INTEGER, PUBLIC, PARAMETER ::   wp =  KIND(REAL_DEF_WP) ! SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
END MODULE modd_wp
