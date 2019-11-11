/*************************************************************************

Include file : unifacparma.h

Purpose: Unifac parameters (replace file read)

Include dependencies:  Included in Unidriver.c

Notes: 5 Type A SOA + butandioic acid + H2O at 25 C
       (Data of Lyman, Reehl, Rosenblatt, 1990)

       NMOL and NFUNC need to match DIMMOL and DIMFUN in unidriver.c

       Parameters inputted as in the original fortran input file, therefore
       a transpose is needed in the unidriver C program for matrices A and NU
       in order to pass them properly into the Fortran unifac routine.

revision History:  1. Developed by Betty Pun, AER, December, 1999 
	              under CARB funding

**************************************************************************/

                                                                               

#ifndef UNIPARM_H
#define UNIPARM_H

/* no. of molecules */  	
int NMOL = 7;

/* no. of functional groups */
int NFUNC = 11;	

/* Z = 10 is a fixed parameter in Unifac */
double Z = 10.0;

/* original file input has temperature, 
but temperature is in main input file now */

/* group volume parameters */
/* dimension of RG is the same as NFUNC */
double RG[DIMFUN] = {0.9011, 0.6744, 0.4469, 1.1167, 0.8886, 0.6605, 1.00, 0.92, 1.6724, 0.998, 1.3013};

/* group surface area parameters */
/* dimension of QG is the same as NFUNC */
double QG[DIMFUN] = {0.8480, 0.5400, 0.2280, 0.8670, 0.6760, 0.4850, 1.20, 1.40, 1.4880, 0.948, 1.2240};

/* no. of groups in each molecule*/
int NU[DIMMOL][DIMFUN] = {
{0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2}, 
{1, 0, 0, 0, 2, 0, 0, 0, 0, 1, 2},
{2, 0, 0, 0, 1, 1, 1, 0, 0, 2, 0},
{1, 2, 1, 1, 0, 0, 1, 0, 1, 0, 1},
{0, 6, 1, 0, 0, 0, 1, 0, 1, 1, 0},
{0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2},
{0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0}};
  
/* no. of groups in each molecule*/
double A[DIMFUN][DIMFUN] = {
{0.0,      0.0,      0.0,      -200.00,  -200.00,  -200.00,  986.500,  1318.00,  476.400,  677.000,  663.500},
{0.0,      0.0,      0.0,      -200.00,  -200.00,  -200.00,  986.500,  1318.00,  476.400,  677.000,  663.500},
{0.0,      0.0,      0.0,      -200.00,  -200.00,  -200.00,  986.500,  1318.00,  476.400,  677.000,  663.500},
{2520.00,  2520.00,  2520.00,  0.0,      0.0,      0.0,      693.900,  634.200,  524.500,  0.0,      730.400},
{2520.00,  2520.00,  2520.00,  0.0,      0.0,      0.0,      693.900,  634.200,  524.500,  0.0,      730.400},
{2520.00,  2520.00,  2520.00,  0.0,      0.0,      0.0,      693.900,  634.200,  524.500,  0.0,      730.400},
{156.400,  156.400,  156.400,  869.400,  869.400,  869.400,  0.0,      353.500,  84.0000,  441.800,  119.000},
{300.000,  300.000,  300.000,  692.700,  692.700,  692.700,  -229.10,  0.0,      -195.40,  -257.30,  -14.090},
{26.7600,  26.7600,  26.7600,  -82.920,  -82.920,  -82.920,  164.500,  472.500,  0.0,      -37.360,  669.400},
{505.700,  505.700,  505.700,  0.0,      0.0,      0.0,      -404.80,  232.700,  128.000,  0.0,      0.0},
{315.300,  315.300,  315.311,  349.200,  349.200,  349.200,  -151.00,  -66.170,  -297.80,  0.0,      0.0}};

#endif

/********************END unifacparam.h**********************************/


