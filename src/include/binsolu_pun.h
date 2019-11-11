/**************************************************************************
Include file: binsolu.h  

Purpose: concentration of solute (micromol solute / microgram water)
	 in a binary solution at RH = Aw = 0, 0.1, 0.2, 0.3 etc

Include dependencies: used in zsr.c

Revision history: Developed by Betty Pun, December 1999, under CARB funding
                  
***********************************************************************/

#ifndef BINARYSOLU_H
#define BINARYSOLU_H

/* graduation of RH scale in molalbin definition, 
starting with RH = 0 (first row), ending with 1 (last row) */
double RHgrad = 0.1;            

/* binary solution molality (umol/ug water) corresponding to RH 
   at 10% graduations --> 11 rows */   
double molalbin[11][NAMOL] = {
{555.600, 555.600, 555.600, 555.600, 555.600, 555.600},
{0.49396, 0.45742, 0.86272, 0.55832, 0.92947, 0.54700},
{0.22574, 0.21651, 0.40057, 0.26300, 0.42461, 0.25189},
{0.13569, 0.13549, 0.24605, 0.16412, 0.25726, 0.15283},
{0.09018, 0.09439, 0.16855, 0.11418, 0.17392, 0.10277},
{0.06232, 0.06918, 0.12160, 0.08375, 0.12435, 0.07213},
{0.04319, 0.05184, 0.08988, 0.06308, 0.09149, 0.05104},
{0.02890, 0.03891, 0.06671, 0.04780, 0.06809, 0.03516},
{0.01740, 0.02848, 0.04860, 0.03574, 0.05035, 0.02212},
{0.00755, 0.01921, 0.03329, 0.02543, 0.03601, 0.00997},
{0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000}};

#endif

/************************END BINSOLU.H *********************************/



