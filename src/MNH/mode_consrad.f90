!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/29 12:06:27
!-----------------------------------------------------------------
!        ###################
         MODULE MODE_CONSRAD
!        ###################
!
!!
!!    PURPOSE
!!    -------
!!      Initialisation of the constante use for the computation of the
!       cloud optical properties not define in the ECMWF radiative scheme
!!
!!    AUTHOR
!!    ------
!!      
!!
!!    MODIFICATIONS
!!    -------------
!!     
IMPLICIT NONE
CONTAINS
SUBROUTINE INI_CONS_PROP_OP()
USE MODD_CRAD,ONLY: XLWSAVI,XLWLILI,XLWC2R2
IMPLICIT NONE
!Longwave parametrisation make with the size
!distribution hypothesis of the 2 moment 
!scheme of Meso-NH  
XLWC2R2( 1,1:4)=(/-1.59254, 0.727208  , -0.764753, 0.139631/)
XLWC2R2( 2,1:4)=(/ -1.8151  , 0.745122  , -0.683449, 0.117338/)
XLWC2R2( 3,1:4)=(/ -1.61771  , 0.704344  , -0.639996, 0.0802414/)
XLWC2R2( 4,1:4)=(/ -1.4645  , 0.726893  , -0.685239, 0.0747236 /)
XLWC2R2( 5,1:4)=(/ -1.58243  , 0.50894  , -0.367452, -0.0189502/)
XLWC2R2( 6,1:4)=(/ -1.10935  , 0.944325  , -1.14812, 0.204628 /)
XLWC2R2( 7,1:4)=(/ -1.30808  , 0.696888  , -0.635871, 0.0290943/)
XLWC2R2( 8,1:4)=(/ -1.33627  , 0.531407  , -0.407361, -0.0412926/)
XLWC2R2( 9,1:4)=(/ -1.33626  , 0.431992  , -0.292512, -0.0723066/)
XLWC2R2( 10,1:4)=(/ -1.29659  , 0.418514  , -0.320223, -0.0570699/)
XLWC2R2( 11,1:4)=(/ -1.0167  , 0.593283  , -0.834146, 0.128019/)
XLWC2R2( 12,1:4)=(/ -0.658514  , 0.495732  , -1.14801, 0.2716/)
XLWC2R2( 13,1:4)=(/ -0.637296  , 0.627436  , -1.32186, 0.324686/)
XLWC2R2( 14,1:4)=(/ -0.77736  , 1.01068  , -1.61091, 0.392262/)
XLWC2R2( 15,1:4)=(/ -1.27301  , 1.76796  , -1.87305, 0.390043/)
XLWC2R2( 16,1:4)=(/ -1.43287  , 0.691154  , -0.116589, -0.229455/)
!Savijarvi (1997) parameterisation
XLWSAVI( 1,1:2)=(/ 0.15,0.005/)
XLWSAVI( 2,1:2)=(/0.15,0.005/)
XLWSAVI( 3,1:2)=(/0.15,0.005/)
XLWSAVI( 4,1:2)=(/0.15,0.005/)
XLWSAVI( 5,1:2)=(/0.15,0.005/)
XLWSAVI( 6,1:2)=(/0.15,0.005/)
XLWSAVI( 7,1:2)=(/0.15,0.005/)
XLWSAVI( 8,1:2)=(/0.15,0.005/)
XLWSAVI( 9,1:2)=(/0.12,0.003/)
XLWSAVI( 10,1:2)=(/0.12,0.003/)
XLWSAVI( 11,1:2)=(/0.25,0.07/)
XLWSAVI( 12,1:2)=(/ 0.40,0.09/)
XLWSAVI( 13,1:2)=(/0.40,0.09/)
XLWSAVI( 14,1:2)=(/0.40,0.09/)
XLWSAVI( 15,1:2)=(/0.40,0.09/)
XLWSAVI( 16,1:2)=(/0.40,0.09/)
!Correction of the Linder and Li (2000)
!parameterisation implement
! in the ECMWF scheme
XLWLILI( 1,1:5)=(/ 8.8116E-02 , -1.2857E-03,   0.81658,-3.9428,   4.6652/)
XLWLILI( 2,1:5)=(/4.1307E-04  ,-5.9631E-05 ,  2.4275 ,-9.0838 ,  9.6069/)
XLWLILI( 3,1:5)=(/-5.7709E-02 ,  9.9071E-04,   3.1118,-9.554,   9.0189/)
XLWLILI( 4,1:5)=(/-5.3069E-02 ,  9.9992E-04,   2.8045,-7.2836,   6.2573/)
XLWLILI( 5,1:5)=(/-2.3627E-02 ,  5.5291E-04,   2.1785,-5.4664,   4.7379/)
XLWLILI( 6,1:5)=(/2.9022E-02  ,-3.9657E-04 ,  1.4902 ,-5.0777,   5.217/)
XLWLILI( 7,1:5)=(/-2.4901E-02 ,  1.6195E-04,  2.9375,-11.437,  12.273/)
XLWLILI( 8,1:5)=(/-0.14269   ,2.2282E-03   ,4.6478,-16.369 , 16.533/)
XLWLILI( 9,1:5)=(/-0.20398   ,3.4708E-03   ,5.2858,-16.603 , 15.392/)
XLWLILI(10,1:5)=(/ -0.18318  , 3.308E-03  , 4.612,-11.55,   8.7086/)
XLWLILI(11,1:5)=(/-0.2042   ,3.72E-03   ,4.8566 ,-11.972 ,  8.6344/)
XLWLILI(12,1:5)=(/-0.14037   ,2.8058E-03   ,3.4969 ,-3.377 , -2.3541/)
XLWLILI(13,1:5)=(/ -0.14037  , 2.8058E-03  , 3.4969,-3.377,  -2.3541/)
XLWLILI(14,1:5)=(/ -0.14037  , 2.8058E-03  , 3.4969,-3.377,  -2.3541/)
XLWLILI(15,1:5)=(/ -0.14037  , 2.8058E-03  , 3.4969,-3.377,  -2.3541/)
XLWLILI(16,1:5)=(/ -0.14037  , 2.8058E-03  , 3.4969, -3.377,  -2.3541/)
END SUBROUTINE INI_CONS_PROP_OP
END MODULE MODE_CONSRAD
