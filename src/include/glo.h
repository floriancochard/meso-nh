/* BKP April 2000 global variables for the OA module */

/* Control Variables */

int zsrflag = 1;   /* flag to use binary solution and zsr if = 1 */
                   /* zsrflag = 0 calls unifac and solves implicit
		      equation for a.c. water = RH using newt1 */ 
     
int Newtflag = 0;      /* flag to use NEWT in Type B module 
			  and type A module with absorption when = 1, 
			  if = 0; don't use NEWT */

int saturationflag = 0;  /* saturationflag = 1 means to 
			    use saturation to determine particulate-phase
                            concentration when inorganic particle is dry. 
			    If = 0, use absorption */


/* General Variables */

double temperature;
double PAOM;            /* as input to Type B, PAOM defined as 
			 non-volatile organic (other OC) */
double VPCrit = 1.0e-8; /* Criterion for setting intital particle 
			   conc, compared to VP in torr or VP 
			   rated by PAOM */
                      

/* Type A variables */

double totA[NAMOL+1], acHP, RH, LWC;   

double g[NAMOL+1];     /* results of gas phase concentrations */
double negcharge;      /* total mole per m3 air of -ve charges */
double KA[NAMOL+1];    /* dry Type A absorption partition constants */

int naero;          /* no. of species in aerosol phase: molecule + ions */
int aidx[NAAERO+1];  /* list of molecular or ionic species present and equations to use) */
int aidxmol[NAMOL+1];/* list of molecular species present (equations to use in absorption) */
int anrerrflag;      /* numerical recipe error flag for Type A */

/* Type A parameters */

double Critsol = 100.0;
double LWCTOL = 0.01;  

int NK[NAMOL+1] = {0, 3, 3, 1, 2, 1, 3};  /* no. of eq. relationship for 
					    each solute, in the order 
					    listed in amina.c */
/* original parameters
double K[NAAERO+1] = {0.0, 0.0245, 3.95e-4, 7.70e-6, 22.01, 3.95e-4, 7.70e-6,
		  1.47e-4, 0.0489, 6.52e-4, 0.0196, 7.34e-3, 3.95e-4, 7.70e-6};
*/

/* below are used in the sensitivity analysis, H1 is updated based on 2
   explicit compounds rather than 3 */
double K[NAAERO+1] = {0.0, 3.87e-5, 3.95e-4, 7.70e-6, 22.01, 3.95e-4, 7.70e-6,
		  1.47e-4, 0.0489, 6.52e-4, 0.0196, 7.34e-3, 3.95e-4, 7.70e-6};
       /* partition parameters H and K
          H is in units of (aq. microgram / m3 air)/(microgram water / m3 air)
                           / g. microgram / m3 )
          estimated based on Suzuki et al., 1992
          K is in units of mol/kg water (same as {H+}) with
          concentrations of molecules in ions in the same mass-based units
          K's of malic acid and glyoxalic acid used respectively for
          compounds that dissociate twice and once */
double MW[NAAERO+1] = {18.0, 104.0, 103.0, 102.0, 184.0, 183.0, 182.0,
                      154.0, 186.0, 185.0, 186.0, 118.0, 117.0, 116.0 };
       /* Molecular weights in the order defined in amain.c */
double DRH[NAMOL+1] = {0.0, 0.79, 0.0, 0.0, 0.0, 0.0, 0.0};
       /* deliquescence humidities for molecules */
double VP[NAMOL+1] = {0.0, 462.22, 0.0127, 1.86e5, 1.11, 71.62, 50.16};    
       /* vapor pressure in units of mass (microgram) per m3 air */    
double VPAtorr[NAMOL+1] = {0.0,8.26e-5,1.28e-9,2.24e-2,1.1e-7,7.16e-6,7.9e-6};
       /* vapor pressure in torr */    


/* Type B variables */

double cb[NBSP+1];      /* total Type B condensable concentrations */
double KB[NBSP+1];      /* Type B partition constants */
int aidxb[NBSP+1];      /* index for Type B equations to solve */    
int bnrerrflag;         /* numerical recipe error flag for type b */  

/* Type B parameters */

double MWB[NBSP+1] = {0.0, 197.0,  164.0, 181.0, 301.0, 170.0};
double VPB[NBSP+1] = {0.0, 3.e-10, 3.e-6, 7.e-6, 5.e-8, 3.e-8};   /* torr */

/* the following are Type B absorbing medium parameters used in 
   subroutine typeb */

double xaom[(NBSPAOM+1)] = {0., 0.4, 0.05, 0.15, 0.12, 0.28};
/* mole fraction of non-volatile organics */

/* double fom = 0.1; */

double MWaom = 280;                                                        

                                                                               
