//MNH_LIC Copyright 2019-2019 CNRS, Meteo-France and Universite Paul Sabatier
//MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
//MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
//MNH_LIC for details. version 1.
//-----------------------------------------------------------------
#ifdef __GFORTRAN__
/* Philippe Marguinaud idea */

void init_gfortran_big_endian_(){
  _gfortran_set_convert (2);
}
void init_gfortran_native_endian_(){
  _gfortran_set_convert (0);
}
#endif
