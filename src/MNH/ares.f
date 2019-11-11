!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
c/////////////////////////////////////////////////////////////////////////////
            
C    Calculate the aerosol chemical speciation and water content.

      SUBROUTINE RPMARES ( SO4, HNO3, NO3, NH3, NH4, RH, TEMP,
     &                     ASO4, ANO3, AH2O, ANH4, GNH3, GNO3) 

C-----------------------------------------------------------------------
C
C Description:
C
C   ARES calculates the chemical composition of a sulfate/nitrate/
C   ammonium/water aerosol based on equilibrium thermodynamics.
C
C   This code considers two regimes depending upon the molar ratio 
C   of ammonium to sulfate. 
C
C   For values of this ratio less than 2,the code solves a cubic for 
C   hydrogen ion molality, HPLUS,  and if enough ammonium and liquid
C   water are present calculates the dissolved nitric acid. For molal
C   ionic strengths greater than 50, nitrate is assumed not to be present. 
C   
C   For values of the molar ratio of 2 or greater, all sulfate is assumed
C   to be ammonium sulfate and a calculation is made for the presence of
C   ammonium nitrate.
C
C   The Pitzer multicomponent approach is used in subroutine ACTCOF to
C   obtain the activity coefficients. Abandoned -7/30/97 FSB 

c   The Bromley method of calculating the activity coefficients is s used
c    in this version

c   The calculation of liquid water
C   is done in subroutine water. Details for both calculations are given
C   in the respective subroutines.
C
C   Based upon MARS due to 
C   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld, 
C   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
C
C   and SCAPE due to 
C   Kim, Seinfeld, and Saxeena, Aerosol Ceience and Technology,
C   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
C
C NOTE: All concentrations supplied to this subroutine are TOTAL
C       over gas and aerosol phases
C
C Parameters:
C 
C  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
C  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
C  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
C  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
C  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
C  RH    : Fractional relative humidity (IN)
C  TEMP  : Temperature in Kelvin (IN)
C  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
C  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
C  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT) 
C  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
C  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
C  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
C  NITR  : Number of iterations for obtaining activity coefficients  (OUT) 
C  NR    : Number of real roots to the cubic in the low ammonia case (OUT)
C 
C Revision History:
C      Who       When        Detailed description of changes
C   ---------   --------  -------------------------------------------
C   S.Roselle   11/10/87  Received the first version of the MARS code
C   S.Roselle   12/30/87  Restructured code
C   S.Roselle   2/12/88   Made correction to compute liquid-phase 
C                         concentration of H2O2.  
C   S.Roselle   5/26/88   Made correction as advised by SAI, for 
C                         computing H+ concentration.
C   S.Roselle   3/1/89    Modified to operate with EM2
C   S.Roselle   5/19/89   Changed the maximum ionic strength from 
C                         100 to 20, for numerical stability.
C   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
C                         using equations for nitrate budget.
C   F.Binkowski 6/18/91   New ammonia poor case which 
C                         omits letovicite.
C   F.Binkowski 7/25/91   Rearranged entire code, restructured
C                         ammonia poor case.
C   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
C                         as SO4--
C   F.Binkowski 12/6/91   Changed the ammonia defficient case so that 
C                         there is only neutralized sulfate (ammonium
C                         sulfate) and sulfuric acid.
C   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement 
C                          with the Cohen et al. (1987)  maximum molality
C                          of 36.2 in Table III.( J. Phys Chem (91) page
C                          4569, and Table IV p 4587.)
C   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
C                         possibility for denomenator becoming zero; 
C                         this involved solving for HPLUS first.
C                         Note that for a relative humidity
C                          less than 50%, the model assumes that there is no 
C                          aerosol nitrate.
C   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)  
C                          Redid logic as follows
C                         1. Water algorithm now follows Spann & Richardson
C                         2. Pitzer Multicomponent method used
C                         3. Multicomponent practical osmotic coefficient 
C                            use to close iterations.
C                         4. The model now assumes that for a water
C                            mass fraction WFRAC less than 50% there is
C                            no aerosol nitrate.
C   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor 
C                         case, and changed the WFRAC criterion to 40%.
C                         For ammonium to sulfate ratio less than 1.0 
C                         all ammonium is aerosol and no nitrate aerosol
C                         exists.
C   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
C                         allow gas-phase ammonia to exist. 
C   F.Binkowski 7/26/95   Changed equilibrium constants to values from 
C                         Kim et al. (1993) 
C   F.Binkowski 6/27/96   Changed to new water format
c   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent 
c                         activity coefficients. The binary activity coefficients 
c                         are the same as the previous version
c   F.Binkowski 8/1/97    Chenged minimum sulfate from 0.0 to 1.0e-6 i.e.
c                         1 picogram per cubic meter
c   P.Wautelet 13/02/2018 use ifdef MNH_REAL to prevent problems with intrinsics on Blue Gene/Q
C
C-----------------------------------------------------------------------
*   K. Suhre    6/6/98    
*
*   note different values of the following constants,
*   the latter ones are used in subroutine awater ..
*
*           MWSO4   = 96.0576
*           mwso4   = 96.0636
*
*           MWNH4   = 18.03858
*           mwnh4   = 18.0985
*
*           MWNO3   = 62.0049
*           mwno3   = 62.0649
*
* --
*
*   increased precision by
*   changing  TOLER2 = 0.001  to   TOLER2 = 0.00001
*
C-----------------------------------------------------------------------
 
      IMPLICIT NONE

C...........INCLUDES and their descriptions

ccc      INCLUDE SUBST_CONST          ! constants

C...........PARAMETERS and their descriptions:

      REAL        MWNACL           ! molecular weight for NaCl
      PARAMETER ( MWNACL = 58.44277 )

      REAL        MWNO3            ! molecular weight for NO3
      PARAMETER ( MWNO3  = 62.0049 ) 

      REAL        MWHNO3           ! molecular weight for HNO3
      PARAMETER ( MWHNO3 = 63.01287 )       

      REAL        MWSO4            ! molecular weight for SO4
      PARAMETER ( MWSO4 = 96.0576 )

      REAL        MWHSO4           ! molecular weight for HSO4
      PARAMETER ( MWHSO4 = MWSO4 + 1.0080 ) 

      REAL        MH2SO4           ! molecular weight for H2SO4
      PARAMETER ( MH2SO4 = 98.07354 ) 

      REAL        MWNH3            ! molecular weight for NH3
      PARAMETER ( MWNH3 = 17.03061 ) 

      REAL        MWNH4            ! molecular weight for NH4
      PARAMETER ( MWNH4 = 18.03858 )

      REAL        MWORG            ! molecular weight for Organic Species
      PARAMETER ( MWORG = 16.0 )

      REAL        MWCL             ! molecular weight for Chloride  
      PARAMETER ( MWCL = 35.453 )

      REAL        MWAIR            ! molecular weight for AIR
      PARAMETER ( MWAIR = 28.964 )

      REAL        MWLCT            ! molecular weight for Letovicite
      PARAMETER ( MWLCT = 3.0 * MWNH4 + 2.0 * MWSO4 + 1.0080 )

      REAL        MWAS             ! molecular weight for Ammonium Sulfate
      PARAMETER ( MWAS = 2.0 * MWNH4 + MWSO4 )

      REAL        MWABS            ! molecular weight for Ammonium Bisulfate 
      PARAMETER ( MWABS = MWNH4 + MWSO4 + 1.0080 )

C...........ARGUMENTS and their descriptions

      REAL        SO4              ! Total sulfate in micrograms / m**3 
ciamodels3
      REAL        HNO3             ! Total nitric acid in micrograms / m**3
      REAL        NO3              ! Total nitrate in micrograms / m**3
      REAL        NH3              ! Total ammonia in micrograms / m**3
      REAL        NH4              ! Total ammonium in micrograms / m**3
      REAL        RH               ! Fractional relative humidity 
      REAL        TEMP             ! Temperature in Kelvin 
      REAL        ASO4             ! Aerosol sulfate in micrograms / m**3 
      REAL        ANO3             ! Aerosol nitrate in micrograms / m**3
      REAL        AH2O             ! Aerosol liquid water content water in micrograms / m**3
      REAL        ANH4             ! Aerosol ammonium in micrograms / m**3
      REAL        GNO3             ! Gas-phase nitric acid in micrograms / m**3
      REAL        GNH3             ! Gas-phase ammonia in micrograms / m**3

C...........SCRATCH LOCAL VARIABLES and their descriptions:
       
      INTEGER     IRH              ! Index set to percent relative humidity  
      INTEGER     NITR             ! Number of iterations for activity coefficients
      INTEGER     NNN              ! Loop index for iterations 
      INTEGER     NR               ! Number of roots to cubic equation for HPLUS

      REAL        A0               ! Coefficients and roots of 
      REAL        A1               ! Coefficients and roots of 
      REAL        A2               ! Coefficients and roots of 
      REAL        AA               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      REAL        BAL              ! internal variables ( high ammonia case)
      REAL        BB               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      REAL        BHAT             ! Variables used for ammonia solubility 
      REAL        CC               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      REAL        CONVT            ! Factor for conversion of units  
      REAL        DD               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      REAL        DISC             ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      REAL        EROR             ! Relative error used for convergence test  
      REAL        FNH3             ! "Free ammonia concentration", that which exceeds TWOSO4       
      REAL        GAMAAB           ! Activity Coefficient for (NH4+, HSO4-)GAMS( 2,3 )
      REAL        GAMAAN           ! Activity coefficient for (NH4+, NO3-) GAMS( 2,2 )
      REAL        GAMAHAT          ! Variables used for ammonia solubility 
      REAL        GAMANA           ! Activity coefficient for (H+ ,NO3-)   GAMS( 1,2 )
      REAL        GAMAS1           ! Activity coefficient for (2H+, SO4--) GAMS( 1,1 )
      REAL        GAMAS2           ! Activity coefficient for (H+, HSO4-)  GAMS( 1,3 )
      REAL        GAMOLD           ! used for convergence of iteration
      REAL        GASQD            ! internal variables ( high ammonia case)
      REAL        HPLUS            ! Hydrogen ion (low ammonia case) (moles / kg water)
      REAL        K1A              ! Equilibrium constant for ammoniua to ammonium
      REAL        K2SA             ! Equilibrium constant for sulfate-bisulfate (aqueous)
      REAL        K3               ! Dissociation constant for ammonium nitrate 
      REAL        KAN              ! Equilibrium constant for ammonium nitrate (aqueous)
      REAL        KHAT             ! Variables used for ammonia solubility 
      REAL        KNA              ! Equilibrium constant for nitric acid (aqueous)   
      REAL        KPH              ! Henry's Law Constant for ammonia       
      REAL        KW               ! Equilibrium constant for water dissociation             
      REAL        KW2              ! Internal variable using KAN 
      REAL        MAN              ! Nitrate (high ammonia case) (moles / kg water)
      REAL        MAS              ! Sulfate (high ammonia case) (moles / kg water)
      REAL        MHSO4            ! Bisulfate (low ammonia case) (moles / kg water)
      REAL        MNA              ! Nitrate (low ammonia case) (moles / kg water)
      REAL        MNH4             ! Ammonium (moles / kg water)
      REAL        MOLNU            ! Total number of moles of all ions
      REAL        MSO4             ! Sulfate (low ammonia case) (moles / kg water)
      REAL        PHIBAR           ! Practical osmotic coefficient      
      REAL        PHIOLD           ! Previous value of practical osmotic coefficient used for convergence of iteration
      REAL        RATIO            ! Molar ratio of ammonium to sulfate
      REAL        RK2SA            ! Internal variable using K2SA
      REAL        RKNA             ! Internal variables using KNA
      REAL        RKNWET           ! Internal variables using KNA
      REAL        RR1
      REAL        RR2
      REAL        STION            ! Ionic strength
      REAL        T1               ! Internal variables for temperature corrections
      REAL        T2               ! Internal variables for temperature corrections
      REAL        T21              ! Internal variables of convenience (low ammonia case)
      REAL        T221             ! Internal variables of convenience (low ammonia case)
      REAL        T3               ! Internal variables for temperature corrections
      REAL        T4               ! Internal variables for temperature corrections
      REAL        T6               ! Internal variables for temperature corrections
      REAL        TNH4             ! Total ammonia and ammonium in micromoles / meter ** 3
      REAL        TNO3             ! Total nitrate in micromoles / meter ** 3
      REAL        TOLER1           ! Tolerances for convergence test 
      REAL        TOLER2           ! Tolerances for convergence test 
      REAL        TSO4             ! Total sulfate in micromoles / meter ** 3
      REAL        TWOSO4           ! 2.0 * TSO4  (high ammonia case) (moles / kg water)
      REAL        WFRAC            ! Water mass fraction 
      REAL        WH2O             ! Aerosol liquid water content (internally)
                                   ! micrograms / meter **3 on output
                                   ! internally it is 10 ** (-6) kg (water) / meter ** 3
                                   ! the conversion factor (1000 g = 1 kg) is applied 
                                   ! for AH2O output
      REAL        WSQD             ! internal variables ( high ammonia case)
      REAL        XNO3             ! Nitrate aerosol concentration in micromoles / meter ** 3
      REAL        XXQ              ! Variable used in quadratic solution
      REAL        YNH4             ! Ammonium aerosol concentration in micromoles / meter** 3
      REAL        ZH2O             ! Water variable saved in case ionic strength too high.
      REAL        ZSO4             ! Total sulfate molality - mso4 + mhso4 (low ammonia case) (moles / kg water)

      REAL        CAT( 2 )         ! Array for cations (1, H+); (2, NH4+) (moles / kg water)
      REAL        AN ( 3 )         ! Array for anions (1, SO4--); (2, NO3-); (3, HSO4-)  (moles / kg water) 
      REAL        CRUTES( 3 )      ! Coefficients and roots of 
      REAL        GAMS( 2, 3 )     ! Array of activity coefficients 
      REAL        MINSO4           ! Minimum value of sulfate laerosol concentration
       PARAMETER( MINSO4 = 1.0E-6 / MWSO4 ) 
      REAL        FLOOR
       PARAMETER( FLOOR = 1.0E-30) ! minimum concentration       

C-----------------------------------------------------------------------
C  begin body of subroutine RPMARES
                                                                         
C...convert into micromoles/m**3
 
ccc      WRITE( 10, * ) 'SO4, NO3, NH3 ', SO4, NO3, NH3
Ciamodels3 merge NH3/NH4 , HNO3,NO3 here
      TSO4 = MAX( 0.0, SO4 / MWSO4  )
      TNO3 = MAX( 0.0, (NO3 / MWNO3 + HNO3 / MWHNO3) )
      TNH4 = MAX( 0.0, (NH3 / MWNH3 + NH4 / MWNH4)  )
ccc      WRITE( 10, * ) 'TSO4, TNO3, TNH3, RH ', TSO4, TNO3, TNH3, RH
 
C...now set humidity index IRH as a percent

      IRH = NINT( 100.0 * RH )

C...Check for valid IRH

      IRH = MAX(  1, IRH )
      IRH = MIN( 99, IRH )
ccc      WRITE(10,*)'RH,IRH ',RH,IRH

C...Specify the equilibrium constants at  correct
C...  temperature.  Also change units from ATM to MICROMOLE/M**3 (for KAN,
C...  KPH, and K3 )
C...  Values from Kim et al. (1993) except as noted.
 
      CONVT = 1.0 / ( 0.082 * TEMP ) 
      T6 = 0.082E-9 *  TEMP
      T1 = 298.0 / TEMP
      T2 = ALOG( T1 )
      T3 = T1 - 1.0
      T4 = 1.0 + T2 - T1
      KNA  = 2.511E+06 *  EXP(  29.17 * T3 + 16.83 * T4 ) * T6
      K1A  = 1.805E-05 *  EXP(  -1.50 * T3 + 26.92 * T4 )
      K2SA = 1.015E-02 *  EXP(   8.85 * T3 + 25.14 * T4 )
      KW   = 1.010E-14 *  EXP( -22.52 * T3 + 26.92 * T4 )
      KPH  = 57.639    *  EXP(  13.79 * T3 - 5.39  * T4 ) * T6
ccc      K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6  
      KHAT =  KPH * K1A / KW  
      KAN  =  KNA * KHAT  

C...Compute temperature dependent equilibrium constant for NH4NO3
C...  ( from Mozurkewich, 1993)

      K3 = EXP( 118.87  - 24084.0 / TEMP -  6.025  * ALOG( TEMP ) )

C...Convert to (micromoles/m**3) **2

      K3 = K3 * CONVT * CONVT

      WH2O   = 0.0
      STION  = 0.0
      AH2O   = 0.0
      MAS    = 0.0
      MAN    = 0.0
      HPLUS  = 0.0
      TOLER1 = 0.00001
      TOLER2 = 0.00001
      NITR   = 0
      NR     = 0
      RATIO  = 0.0
      GAMAAN = 1.0
      GAMOLD = 1.0

C...set the ratio according to the amount of sulfate and nitrate

      IF ( TSO4 .GT. MINSO4 ) THEN
        RATIO = TNH4 / TSO4

C...If there is no sulfate and no nitrate, there can be no ammonium
C...  under the current paradigm. Organics are ignored in this version.

      ELSE 
      
       IF ( TNO3 .EQ. 0.0 ) THEN

C *** If there is very little sulfate and no nitrate set concentrations
c      to a very small value and return.        
          ASO4 = MAX(FLOOR, ASO4)
          ANO3 = MAX(FLOOR, ANO3 )          
          WH2O = 0.0
          AH2O = 0.0
          GNH3 = MAX(FLOOR,GNH3)
          GNO3 = MAX(FLOOR,GNO3)
          RETURN
       END IF
       
C...For the case of no sulfate and nonzero nitrate, set ratio to 5
C...  to send the code to the high ammonia case

        RATIO = 5.0
      END IF 

 
C....................................
C......... High Ammonia Case ........
C....................................
 
      IF ( RATIO .GT. 2.0 ) THEN
 
        GAMAAN = 0.1

C...Set up twice the sulfate for future use.

        TWOSO4 = 2.0 * TSO4
        XNO3 = 0.0            
        YNH4 = TWOSO4

C...Treat different regimes of relative humidity 

C...ZSR relationship is used to set water levels. Units are
C...  10**(-6) kg water/ (cubic meter of air)
C...  start with ammomium sulfate solution without nitrate

      CALL awater(IRH,TSO4,YNH4,TNO3,AH2O) !**** note TNO3
        WH2O = 1.0E-3 * AH2O  
        ASO4 = TSO4 * MWSO4
        ANO3 = 0.0
        ANH4 = YNH4 * MWNH4
        WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )
ccc        IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water       
        IF ( WFRAC .LT. 0.2 ) THEN
 
C..."dry" ammonium sulfate and ammonium nitrate
C...  compute free ammonia

          FNH3 = TNH4 - TWOSO4
          CC = TNO3 * FNH3 - K3

C...check for not enough to support aerosol      

          IF ( CC .LE. 0.0 ) THEN
            XNO3 = 0.0
          ELSE
            AA = 1.0
            BB = -( TNO3 + FNH3 ) 
            DISC = BB * BB - 4.0 * CC

C...Check for complex roots of the quadratic
C...  set nitrate to zero and RETURN if complex roots are found

            IF ( DISC .LT. 0.0 ) THEN
              XNO3 = 0.0
              AH2O = 1000.0 * WH2O
              YNH4 = TWOSO4
              GNO3 = TNO3 * MWHNO3
              GNH3 = ( TNH4 - YNH4 ) * MWNH3
              ASO4 = TSO4 * MWSO4
              ANO3 = 0.0
              ANH4 = YNH4 * MWNH4
              RETURN
            END IF

C...to get here, BB .lt. 0.0, CC .gt. 0.0 always      

            DD = SQRT( DISC )
            XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )

C...Since both roots are positive, select smaller root.      

            XNO3 = MIN( XXQ / AA, CC / XXQ )
          
          END IF
          AH2O = 1000.0 * WH2O
          YNH4 = 2.0 * TSO4 + XNO3
          GNO3 = ( TNO3 - XNO3 ) * MWHNO3
          GNH3 = ( TNH4 - YNH4 ) * MWNH3
          ASO4 = TSO4 * MWSO4
          ANO3 = XNO3 * MWNO3
          ANH4 = YNH4 * MWNH4
          RETURN

        END IF

C...liquid phase containing completely neutralized sulfate and
C...  some nitrate.  Solve for composition and quantity.
 
        MAS = TSO4 / WH2O
        MAN = 0.0
        XNO3 = 0.0
        YNH4 = TWOSO4
        PHIOLD = 1.0

C...Start loop for iteration
 
C...The assumption here is that all sulfate is ammonium sulfate,
C...  and is supersaturated at lower relative humidities.
 
        DO 1501 NNN = 1, 150 
          NITR = NNN
          GASQD = GAMAAN * GAMAAN
          WSQD = WH2O * WH2O
          KW2 = KAN * WSQD / GASQD
          AA = 1.0 - KW2
          BB = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
          CC = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

C...This is a quadratic for XNO3 [MICROMOLES / M**3] of nitrate in solution.

          DISC = BB * BB - 4.0 * AA * CC

C...Check for complex roots, if so set nitrate to zero and RETURN 

          IF ( DISC .LT. 0.0 ) THEN
            XNO3 = 0.0
            AH2O = 1000.0 * WH2O
            YNH4 = TWOSO4
            GNO3 = TNO3 * MWHNO3
            GNH3 = ( TNH4 - YNH4 ) * MWNH3
            ASO4 = TSO4 * MWSO4
            ANO3 = 0.0
            ANH4 = YNH4 * MWNH4
ccc            WRITE( 10, * ) ' COMPLEX ROOTS '
            RETURN
          END IF

          DD = SQRT( DISC )
          XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )
          RR1 = XXQ / AA
          RR2 = CC / XXQ

C...choose minimum positve root         

          IF ( ( RR1 * RR2 ) .LT. 0.0 ) THEN
            XNO3 = MAX( RR1, RR2 )
          ELSE 
            XNO3 = MIN( RR1, RR2 )
          END IF

          XNO3 = MIN( XNO3, TNO3 )

C...This version assumes no solid sulfate forms (supersaturated ) 
C...  Now update water

          CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

C...ZSR relationship is used to set water levels. Units are
C...  10**(-6) kg water/ (cubic meter of air)
C...  The conversion from micromoles to moles is done by the units of WH2O.

          WH2O = 1.0E-3 * AH2O

C...Ionic balance determines the ammonium in solution.

          MAN = XNO3 / WH2O
          MAS = TSO4 / WH2O
          MNH4 = 2.0 * MAS + MAN
          YNH4 = MNH4 * WH2O

C...MAS, MAN and MNH4 are the aqueous concentrations of sulfate, nitrate,
C...  and ammonium in molal units (moles/(kg water) ).

          STION = 3.0 * MAS + MAN
          CAT( 1 ) = 0.0
          CAT( 2 ) = MNH4 
          AN ( 1 ) = MAS
          AN ( 2 ) = MAN
          AN ( 3 ) = 0.0
          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )
          GAMAAN = GAMS( 2, 2 )

C...Use GAMAAN for convergence control

          EROR = ABS( GAMOLD - GAMAAN ) / GAMOLD
          GAMOLD = GAMAAN

C...Check to see if we have a solution

          IF ( EROR .LE. TOLER1 ) THEN 
ccc            WRITE( 11, * ) RH, STION, GAMS( 1, 1 ),GAMS( 1, 2 ), GAMS( 1, 3 ),
ccc     &      GAMS( 2, 1 ), GAMS( 2, 2 ), GAMS( 2, 3 ), PHIBAR

            ASO4 = TSO4 * MWSO4
            ANO3 = XNO3 * MWNO3
            ANH4 = YNH4 * MWNH4
            GNO3 = ( TNO3 - XNO3 ) * MWHNO3
            GNH3 = ( TNH4 - YNH4 ) * MWNH3
            AH2O = 1000.0 * WH2O
            RETURN
          END IF

1501    CONTINUE

C...If after NITR iterations no solution is found, then:

        ASO4 = TSO4 * MWSO4
        ANO3 = 0.0
        YNH4 = TWOSO4
        ANH4 = YNH4 * MWNH4
        CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )
        GNO3 = TNO3 * MWHNO3
        GNH3 = ( TNH4 - YNH4 ) * MWNH3
        RETURN

      ELSE
       
C......................................
C......... Low Ammonia Case ...........
C......................................
      
C...coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)

C...All cases covered by this logic
 
        WH2O = 0.0
        CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
        WH2O = 1.0E-3 * AH2O
        ZH2O = AH2O

C...convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
C...  per cubic meter of air (1000 g = 1 kg)
      
        ASO4 = TSO4 * MWSO4
        ANH4 = TNH4 * MWNH4
        ANO3 = 0.0
        GNO3 = TNO3 * MWHNO3
        GNH3 = 0.0 

C...Check for zero water.      

        IF ( WH2O .EQ. 0.0 ) RETURN
        ZSO4 = TSO4 / WH2O 

C...ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4      

ccc         IF ( ZSO4 .GT. 11.0 ) THEN

C...do not solve for aerosol nitrate for total sulfate molality
C...  greater than 11.0 because the model parameters break down
C...  greater than  9.0 because the model parameters break down

        IF ( ZSO4 .GT. 9.0 ) THEN   ! 18 June 97
          RETURN
        END IF

C...First solve with activity coeffs of 1.0, then iterate.

        PHIOLD = 1.0
        GAMANA = 1.0
        GAMAS1 = 1.0
        GAMAS2 = 1.0
        GAMAAB = 1.0
        GAMOLD = 1.0

C...All ammonia is considered to be aerosol ammonium. 

        MNH4 = TNH4 / WH2O

C...MNH4 is the molality of ammonium ion.

        YNH4 = TNH4
      
C...loop for iteration
 
        DO 1601 NNN = 1, 150
          NITR = NNN

C...set up equilibrium constants including activities
C...  solve the system for hplus first then sulfate & nitrate

          RK2SA = K2SA * GAMAS2 * GAMAS2 / ( GAMAS1 * GAMAS1 * GAMAS1 )
          RKNA = KNA / ( GAMANA * GAMANA )
          RKNWET = RKNA * WH2O       
          T21  = ZSO4 - MNH4
          T221 = ZSO4 + T21

C...set up coefficients for cubic       

          A2 = RK2SA + RKNWET - T21
          A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET )
     &       - RK2SA * ZSO4 - RKNA * TNO3
          A0 = - (T21 * RK2SA * RKNWET
     &       + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )
          CALL CUBIC ( A2, A1, A0, NR, CRUTES )
       
C...Code assumes the smallest positive root is in CRUTES(1)
 
          HPLUS = CRUTES( 1 )
          BAL = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0
          MSO4 = RK2SA * ZSO4 / ( HPLUS + RK2SA )   ! molality of sulfate ion
          MHSO4 = ZSO4 - MSO4                       ! molality of bisulfate ion
          MNA = RKNA * TNO3 / ( HPLUS + RKNWET )    ! molality of nitrate ion
          MNA = MAX( 0.0, MNA )
          MNA = MIN( MNA, TNO3 / WH2O )
          XNO3 = MNA * WH2O
          ANO3 = MNA * WH2O * MWNO3
          GNO3 = ( TNO3 - XNO3 ) * MWHNO3
        
C...Calculate ionic strength      

          STION = 0.5 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0 * MSO4 )
          
C...Update water

          CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

C...Convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
C...  per cubic meter of air (1000 g = 1 kg)                       

          WH2O = 1.0E-3 * AH2O 
          CAT( 1 ) = HPLUS
          CAT( 2 ) = MNH4
          AN ( 1 ) = MSO4
          AN ( 2 ) = MNA
          AN ( 3 ) = MHSO4

          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )

          GAMANA = GAMS( 1, 2 )
          GAMAS1 = GAMS( 1, 1 )
          GAMAS2 = GAMS( 1, 3 )
          GAMAAN = GAMS( 2, 2 )

          GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
          BHAT = KHAT * GAMAHAT 
ccc          EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
ccc          PHIOLD = PHIBAR
          EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD 
          GAMOLD = GAMAHAT   

C...write out molalities and activity coefficient
C...  and return with good solution

          IF ( EROR .LE. TOLER2 ) THEN
ccc            WRITE(12,*) RH, STION,HPLUS,ZSO4,MSO4,MHSO4,MNH4,MNA
ccc            WRITE(11,*) RH, STION, GAMS(1,1),GAMS(1,2),GAMS(1,3),
ccc     &                  GAMS(2,1),GAMS(2,2),GAMS(2,3), PHIBAR
            RETURN
          END IF

1601    CONTINUE     

C...after NITR iterations, failure to solve the system, no ANO3

        GNO3 = TNO3 * MWHNO3
        ANO3 = 0.0
        CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )      
        RETURN
            
      END IF   ! ratio .gt. 2.0

      END ! end RPMares
C ///////////////////////////////////////////////////
C ///////////////////////////////////////////////////
       SUBROUTINE ACTCOF ( CAT, AN, GAMA, MOLNU, PHIMULT )

C-----------------------------------------------------------------------
C
C DESCRIPTION:
C
C  This subroutine computes the activity coefficients of (2NH4+,SO4--),
C  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
C  multicomponent solution, using Bromley's model and Pitzer's method.
C
C REFERENCES:
C
C   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
C     in aqueous solutions.  AIChE J. 19, 313-320.
C
C   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of 
C     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.
C
C   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures 
C     of strong acids over saline solutions - I HNO3, 
C     Atmos. Environ. (22): 91-100
C
C   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
C     and mean activity and osmotic coefficients of 0-100% nitric acid 
C     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380
C
C   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
C     general equilibrium model for inorganic multicomponent atmospheric
C     aerosols.  Atmos. Environ. 21(11), 2453-2466.
C


C
C ARGUMENT DESCRIPTION:
C
C     CAT(1) : conc. of H+    (moles/kg)
C     CAT(2) : conc. of NH4+  (moles/kg)
C     AN(1)  : conc. of SO4-- (moles/kg)
C     AN(2)  : conc. of NO3-  (moles/kg)
C     AN(3)  : conc. of HSO4- (moles/kg)
C     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
C     GAMA(2,2)    :  "    "     "       "       "    "  (NH4+,NO3-)
C     GAMA(2,3)    :  "    "     "       "       "    "  (NH4+. HSO4-)
C     GAMA(1,1)    :  "    "     "       "       "    "  (2H+,SO4--)
C     GAMA(1,2)    :  "    "     "       "       "    "  (H+,NO3-)
C     GAMA(1,3)    :  "    "     "       "       "    "  (H+,HSO4-)
C     MOLNU   : the total number of moles of all ions.
C     PHIMULT : the multicomponent paractical osmotic coefficient.
C
C REVISION HISTORY:
C      Who       When        Detailed description of changes
C   ---------   --------  -------------------------------------------
C   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
C                         new routine using a method described by Pilinis
C                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-2466.
C   S.Roselle   7/30/97   Modified for use in Models-3
C   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

C...........INCLUDES and their descriptions

c      INCLUDE SUBST_XSTAT     ! M3EXIT status codes
C....................................................................

      INTEGER    XSTAT0       ! Normal, successful completion
      PARAMETER (XSTAT0 = 0)
      INTEGER    XSTAT1       ! File I/O error
      PARAMETER (XSTAT1 = 1)
      INTEGER    XSTAT2       ! Execution error
      PARAMETER (XSTAT2 = 2)
      INTEGER    XSTAT3       ! Special  error
      PARAMETER (XSTAT3 = 3)

      CHARACTER*120 XMSG

C...........PARAMETERS and their descriptions:

      INTEGER      NCAT                 ! number of cations
      PARAMETER  ( NCAT = 2 )

      INTEGER      NAN                  ! number of anions
      PARAMETER  ( NAN = 3 )

C...........ARGUMENTS and their descriptions

      REAL         MOLNU                ! tot # moles of all ions
      REAL         PHIMULT              ! multicomponent paractical osmotic coef
      REAL         CAT( NCAT )          ! cation conc in moles/kg (input)
      REAL         AN ( NAN )           ! anion conc in moles/kg (input)
      REAL         GAMA( NCAT, NAN )    ! mean molal ionic activity coefs

C...........SCRATCH LOCAL VARIABLES and their descriptions:

      CHARACTER*16 PNAME            ! driver program name
      SAVE         PNAME

      INTEGER      IAN                  ! anion indX
      INTEGER      ICAT                 ! cation indX

      REAL         FGAMA                ! 
      REAL         I                    ! ionic strength 
      REAL         R                    ! 
      REAL         S                    ! 
      REAL         TA                   ! 
      REAL         TB                   ! 
      REAL         TC                   ! 
      REAL         TEXPV                ! 
      REAL         TRM                  ! 
      REAL         TWOI                 ! 2*ionic strength
      REAL         TWOSRI               ! 2*sqrt of ionic strength
      REAL         ZBAR                 ! 
      REAL         ZBAR2                ! 
      REAL         ZOT1                 ! 
      REAL         SRI                  ! square root of ionic strength 
      REAL         F2( NCAT )           ! 
      REAL         F1( NAN )            ! 
      REAL         ZP( NCAT )           ! absolute value of charges of cation
      REAL         ZM( NAN )            ! absolute value of charges of anion
      REAL         BGAMA ( NCAT, NAN )  ! 
      REAL         X     ( NCAT, NAN )  ! 
      REAL         M     ( NCAT, NAN )  ! molality of each electrolyte
      REAL         LGAMA0( NCAT, NAN )  ! binary activity coefficients
      REAL         Y     ( NAN, NCAT )  ! 
      REAL         BETA0 ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         BETA1 ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         CGAMA ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         V1    ( NCAT, NAN )  ! number of cations in electrolyte formula
      REAL         V2    ( NCAT, NAN )  ! number of anions in electrolyte formula

      DATA         ZP / 1.0, 1.0 /
      DATA         ZM / 2.0, 1.0, 1.0 /
      DATA          XMSG / ' ' /
      DATA         PNAME / 'ACTCOF' /

C *** Sources for the coefficients BETA0, BETA1, CGAMA:
 
C *** (1,1);(1,3)  - Clegg & Brimblecombe (1988)
C *** (2,3)        - Pilinis & Seinfeld (1987), cgama different 
C *** (1,2)        - Clegg & Brimblecombe (1990)
C *** (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)                                 
      
c *** now set the basic constants, BETA0, BETA1, CGAMA 

      DATA BETA0(1,1) /2.98E-2/,    BETA1(1,1) / 0.0/, 
     &     CGAMA(1,1) / 4.38E-2/                                ! 2H+SO4-
     
      DATA BETA0(1,2) / 1.2556E-1/,   BETA1(1,2) / 2.8778E-1/, 
     &     CGAMA(1,2) / -5.59E-3/                               ! HNO3
     
      DATA BETA0(1,3) / 2.0651E-1/,   BETA1(1,3) / 5.556E-1/,
     &     CGAMA(1,3) /0.0/                                     ! H+HSO4-
     
      DATA BETA0(2,1) /4.6465E-2/,   BETA1(2,1) /-0.54196/, 
     &     CGAMA(2,1) /-1.2683E-3/                              ! (NH4)2SO4
     
      DATA BETA0(2,2) /-7.26224E-3/, BETA1(2,2) /-1.168858/, 
     &     CGAMA(2,2) /3.51217E-5/                              ! NH4NO3
     
      DATA BETA0(2,3) / 4.494E-2/,    BETA1(2,3) / 2.3594E-1/,
     &     CGAMA(2,3) /-2.962E-3/                               ! NH4HSO4

      DATA V1(1,1), V2(1,1) / 2.0, 1.0 /     ! 2H+SO4-
      DATA V1(2,1), V2(2,1) / 2.0, 1.0 /     ! (NH4)2SO4
      DATA V1(1,2), V2(1,2) / 1.0, 1.0 /     ! HNO3 
      DATA V1(2,2), V2(2,2) / 1.0, 1.0 /     ! NH4NO3
      DATA V1(1,3), V2(1,3) / 1.0, 1.0 /     ! H+HSO4-
      DATA V1(2,3), V2(2,3) / 1.0, 1.0 /     ! NH4HSO4

C-----------------------------------------------------------------------
C  begin body of subroutine ACTCOF

C...compute ionic strength

      I = 0.0

      DO ICAT = 1, NCAT
        I = I + CAT( ICAT ) * ZP( ICAT ) * ZP( ICAT )
      END DO

      DO IAN = 1, NAN
        I = I + AN( IAN ) * ZM( IAN ) * ZM( IAN )
      END DO

      I = 0.5 * I

C...check for problems in the ionic strength

      IF ( I .EQ. 0.0 ) THEN

        DO IAN = 1, NAN
          DO ICAT = 1, NCAT
            GAMA( ICAT, IAN ) = 0.0
          END DO
        END DO
        
c        XMSG = 'Ionic strength is zero...returning zero activities'
c        CALL M3WARN ( PNAME, 0, 0, XMSG )
        write(*,*)  'Ionic strength is zero...returning zero activities'
        RETURN

      ELSE IF ( I .LT. 0.0 ) THEN
c        XMSG = 'Ionic strength below zero...negative concentrations'
c        CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT1 )
CALL ABORT
       STOP 'Ionic strength below zero...negative concentrations'
      END IF

C...compute some essential expressions

      SRI    = SQRT( I )
      TWOSRI = 2.0 * SRI
      TWOI   = 2.0 * I
      TEXPV  = 1.0 - EXP( -TWOSRI ) * ( 1.0 + TWOSRI - TWOI )
      R      = 1.0 + 0.75 * I
      S      = 1.0 + 1.5  * I
      ZOT1   = 0.511 * SRI / ( 1.0 + SRI )

C...Compute binary activity coeffs

      FGAMA = -0.392 * ( ( SRI / ( 1.0 + 1.2 * SRI )
     &      + ( 2.0 / 1.2 ) * ALOG( 1.0 + 1.2 * SRI ) ) )

      DO ICAT = 1, NCAT
        DO IAN = 1, NAN

        BGAMA( ICAT, IAN ) = 2.0 * BETA0( ICAT, IAN )
     &                    + ( 2.0 * BETA1( ICAT, IAN ) / ( 4.0 * I ) )
     &                    * TEXPV

C...compute the molality of each electrolyte for given ionic strength

          M( ICAT, IAN ) = ( CAT( ICAT )**V1( ICAT, IAN )
     &                   *   AN( IAN )**V2( ICAT, IAN ) )**( 1.0
     &                   / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) ) )

C...calculate the binary activity coefficients

          LGAMA0( ICAT, IAN ) = ( ZP( ICAT ) * ZM( IAN ) * FGAMA
     &                      + M( ICAT, IAN )
     &                      * ( 2.0 * V1( ICAT, IAN ) * V2( ICAT, IAN )
     &                      / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )
     &                      * BGAMA( ICAT, IAN ) )
     &                      + M( ICAT, IAN ) * M( ICAT, IAN )
     &                      * ( 2.0 * ( V1( ICAT, IAN )
     &                      * V2( ICAT, IAN ) )**1.5
     &                      / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )
     &                      * CGAMA( ICAT, IAN ) ) ) / 2.302585093

        END DO
      END DO

C...prepare variables for computing the multicomponent activity coeffs

      DO IAN = 1, NAN
        DO ICAT = 1, NCAT
          ZBAR = ( ZP( ICAT ) + ZM( IAN ) ) * 0.5
          ZBAR2 = ZBAR * ZBAR
          Y( IAN, ICAT ) = ZBAR2 * AN( IAN ) / I
          X( ICAT, IAN ) = ZBAR2 * CAT( ICAT ) / I
        END DO
      END DO

      DO IAN = 1, NAN
        F1( IAN ) = 0.0
        DO ICAT = 1, NCAT
          F1( IAN ) = F1( IAN ) + X( ICAT, IAN ) * LGAMA0( ICAT, IAN )
     &              + ZOT1 * ZP( ICAT ) * ZM( IAN ) * X( ICAT, IAN )
        END DO
      END DO

      DO ICAT = 1, NCAT
        F2( ICAT ) = 0.0
        DO IAN = 1, NAN
          F2( ICAT ) = F2( ICAT ) + Y( IAN, ICAT ) * LGAMA0( ICAT, IAN )
     &               + ZOT1 * ZP( ICAT ) * ZM( IAN ) * Y( IAN, ICAT )
        END DO
      END DO

C...now calculate the multicomponent activity coefficients

      DO IAN = 1, NAN
        DO ICAT = 1, NCAT

          TA = -ZOT1 * ZP( ICAT ) * ZM( IAN )
          TB = ZP( ICAT ) * ZM( IAN ) / ( ZP( ICAT ) + ZM( IAN ) )
          TC = ( F2( ICAT ) / ZP( ICAT ) + F1( IAN ) / ZM( IAN ) )
          TRM = TA + TB * TC

          IF ( TRM .GT. 30.0 ) THEN
            GAMA( ICAT, IAN ) = 1.0E+30
c            XMSG = 'Multicomponent activity coefficient is extremely large'
c            CALL M3WARN ( PNAME, 0, 0, XMSG )
CALL ABORT
           STOP 'Multicomponent activity coefficient is extremely large'
          ELSE
            GAMA( ICAT, IAN ) = 10.0**TRM
          END IF

        END DO
      END DO
 
      RETURN
      END ! actcof 

C ///////////////////////////////////


*------------------------------------------------------------------
* subroutine  to find the roots of a cubic equation / 3rd order polynomial
* formulae can be found in numer. recip.  on page 145
*   kiran  developed  this version on 25/4/1990
*   dr. francis binkowski modified the routine on 6/24/91, 8/7/97
C ***
C234567
      subroutine cubic(a2,a1,a0,nr,crutes)
      implicit none
      integer  nr      
      real a2,a1,a0,crutes(3)
      real*8 qq,rr,a2sq,theta, sqrt3, one3rd
      real*8 dum1,dum2,part1,part2,part3,rrsq,phi,yy1,yy2,yy3
      real*8 costh, sinth
      data sqrt3/1.732050808/, one3rd/0.333333333/
	  a2sq=a2*a2
	  qq=(a2sq-3.*a1)/9.
	  rr=( a2*(2.*a2sq - 9.*a1) + 27.*a0 )/54.
* CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
      dum1=qq*qq*qq 
      rrsq=rr*rr
      dum2=dum1 - rrsq
      if(dum2.ge.0.) then
* NOW WE HAVE THREE REAL ROOTS
         phi=sqrt(dum1)
         if(abs(phi).lt.1.e-20) then 
C           write(10,*) ' cubic phi small, phi = ',phi
            crutes(1) = 0.0 
            crutes(2) = 0.0
            crutes(3) = 0.0
            nr = 0            
		   stop 
         end if
         theta=acos(rr/phi)/3.0
         costh = cos(theta)
         sinth = sin(theta)
c *** use trig identities to simplify the expressions 
c *** binkowski's modification
         part1=sqrt(qq)
         yy1=part1*costh
         yy2=yy1-a2/3.0
         yy3=sqrt3*part1*sinth
         crutes(3) = -2.0*yy1 - a2/3.0
         crutes(2) = yy2 + yy3
         crutes(1) = yy2 - yy3
C *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
         if(crutes(1) .lt. 0.0) crutes(1) = 1.0e9
         if(crutes(2) .lt. 0.0) crutes(2) =1.0e9
         if(crutes(3) .lt. 0.0) crutes(3) = 1.0e9
c *** put smallest positive root in crutes(1)
         crutes(1)=min( crutes(1),crutes(2),crutes(3))
         nr=3
      else  ! dum IS NEGATIVE
C     NOW HERE WE HAVE ONLY ONE REAL ROOT
         part1=sqrt(rrsq-dum1)
         part2=abs(rr)
         part3=(part1+part2)**one3rd
         crutes(1) = 
#if (MNH_REAL == 8)
     &        -sign(1.0E0,rr) * ( part3 + (qq/part3) ) - a2/3. 
#else
     &        -sign(1.0D0,rr) * ( part3 + (qq/part3) ) - a2/3. 
#endif
         crutes(2)=0.
         crutes(3)=0.
      nr=1
      end if
	  return
	  end ! cubic	
c //////////////////////////  


c
c
c##############################################################################
      subroutine awater(irhx,mso4,mnh4,mno3, wh2o)
c NOTE!!! wh2o is returned in micrograms / cubic meter
c         mso4,mnh4,mno3 are in microMOLES / cubic meter      
c
c  This  version uses polynomials rather than tables, and uses empirical
c polynomials for the mass fraction of solute (mfs) as a function of water activity
c   where:
c     
c            mfs = ms / ( ms + mw) 
c             ms is the mass of solute
c             mw is the mass of water.
c
c  Define y = mw/ ms
c
c  then  mfs = 1 / (1 + y)
c
c    y can then be obtained from the values of mfs as
c
c             y = (1 - mfs) / mfs  
c
c
c     the aerosol is assumed to be in a metastable state if the rh is
c     is below the rh of deliquescence, but above the rh of crystallization.
c
c     ZSR interpolation is used for sulfates with x ( the molar ratio of
c     ammonium to sulfate in eh range 0 <= x <= 2, by sections.
c     section 1: 0 <= x < 1
c     section 2: 1 <= x < 1.5
c     section 3: 1.5 <= x < 2.0
c     section 4: 2 <= x 
c     In sections 1 through 3, only the sulfates can affect the amount of water
c     on the particles.
c     In section 4, we have fully neutralized sulfate, and extra ammonium which
c     allows more nitrate to be present. Thus, the ammount of water is calculated
c     using ZSR for ammonium sulfate and ammonium nitrate. Crystallization is
c     assumed to occur in sections 2,3,and 4. See detailed discussion below. 
c
c

c definitions:
c     mso4, mnh4, and mno3 are the number of micromoles/(cubic meter of air) 
c      for sulfate, ammonium, and nitrate respectively
c     irhx is the relative humidity (%)
c     wh2o is the returned water amount in micrograms / cubic meter of air
c     x is the molar ratio of ammonium to sulfate 
c     y0,y1,y1.5, y2 are the water contents in mass of water/mass of solute
c     for pure aqueous solutions with x equal 1, 1.5, and 2 respectively.
c     y3 is the value of the mass ratio of water to solute for
c     a pure ammonium nitrate  solution.
c
c
coded by Dr. Francis S. Binkowski, 4/8/96.
c
      implicit none
      integer irhx, irh 
      real mso4,mnh4,mno3     
      real tso4,tnh4,tno3, wh2o, x
      real  aw, awc
      real  poly4, poly6     
      real mfs0,mfs1,mfs15, mfs2
      real c0(4), c1(4), c15(4), c2(4)
      real y, y0,y1,y15,y2,y3, y40, y140, y1540, yc
      real kSO4(6),kNO3(6), mfsSO4,mfsNO3

c      
c
      real mwso4, mwnh4, mwno3, mw2, mwano3

c *** molecular weights:      
      parameter(
     &           mwso4   = 96.0636,
     &           mwnh4   = 18.0985,
     &           mwno3   = 62.0649,
     &           mw2     =  mwso4 + 2.0 * mwnh4,                
     &           mwano3  = mwno3 + mwnh4 )
 
c     The polynomials use data for aw as a function of mfs from Tang and
c     Munkelwitz, JGR 99: 18801-18808, 1994.
c     The polynomials were fit to Tang's values of water activity as a 
c     function of mfs. 
     
c *** coefficients of polynomials fit to Tang and Munkelwitz data 
c     now give mfs as a function of water activity. 
  
      data c1/0.9995178, -0.7952896, 0.99683673, -1.143874/
      data c15/1.697092,-4.045936, 5.833688, -3.463783/
      data c2/2.085067, -6.024139, 8.967967, -5.002934/

C *** the following coefficients are a fit to the data in Table 1 of
c     Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975 
c      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
C *** New data fit to data from 
c       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975 
c       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
c       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
       data c0/ 0.798079, -1.574367, 2.536686, -1.735297 / 
      
       
c *** polynomials for ammonium nitrate and ammonium sulfate are from:       
c     Chan et al.1992, Atmospheric Environment (26A): 1661-1673.

      data kNO3/0.2906, 6.83665, -26.9093,
     &          46.6983, -38.803, 11.8837/  
      data kSO4/ 2.27515, -11.147, 36.3369, 
     &        -64.2134, 56.8341, -20.0953/
      

c *** check range of per cent relative humidity 
       irh = irhx
       irh = max(1,irh)
       irh = min(irh,100)
       aw  = float(irh) / 100.0 ! water activity = fractional relative humidity
       tso4 = max( mso4 , 0.0 )
       tnh4 = max( mnh4 , 0.0 )
       tno3 = max( mno3 , 0.0 )
       x = 0.0
c *** if there is non-zero sulfate calculate the molar ratio      
       if (tso4 .gt. 0.0 ) then
         x = tnh4 / tso4
       else
c *** otherwise check for non-zero nitrate and ammonium       
         if ( tno3 .gt. 0.0 . and. tnh4 .gt. 0.0 ) x = 10.0
       end if
c
c 

c *** begin screen on x for calculating wh2o      
       if ( x .lt. 1.0 ) then
c 
          mfs0 = poly4(c0,aw)
          mfs1 = poly4(c1,aw)
          y0 = (1.0 - mfs0 ) / mfs0
          y1 = (1.0 - mfs1 ) / mfs1
          y = (1.0 - x) * y0 + x * y1       
             
c       
       else if ( x .lt. 1.5) then
c
         if ( irh .ge. 40 ) then         
            mfs1  = poly4(c1,aw)
            mfs15 = poly4(c15,aw)
            y1  = (1.0 - mfs1 ) / mfs1
            y15 = (1.0 - mfs15) / mfs15
            y = 2.0 * ( y1 * (1.5 - x) + y15 *( x - 1.0) ) 
         else
c *** set up for crystalization

c *** Crystallization is done as follows:
c      For 1.5 <= x, crystallization is assumed to occur at rh = 0.4
c      For x <= 1.0, crystallization is assumed to occur at an rh < 0.01,
c      and since the code does not allow ar rh < 0.01, crystallization
c      is assumed not to occur in this range.
c      For 1.0 <= x <= 1.5 the crystallization curve is a straignt line
c      from a value of y15 at rh = 0.4 to a value of zero at y1. From
C      point B to point A in the diagram. 
c      The algorithm does a double interpolation to calculate the amount of 
c      water. 
c
c        y1(0.40)               y15(0.40)
c         +                     + Point B 
c                              
c
c
c
c         +--------------------+
c       x=1                   x=1.5
c      Point A
c
c
                
            awc = 0.80 * (x - 1.0) ! rh along the crystallization curve.             
            y = 0.0 
             if ( aw .ge. awc ) then ! interpolate using crystalization curve
               mfs1  = poly4(c1,0.40)
               mfs15 = poly4(c15,0.40)
               y140  = (1.0 - mfs1 ) / mfs1
               y1540 = (1.0 - mfs15) / mfs15
               y40 = 2.0 * ( y140 * (1.5 - x) + y1540 *( x - 1.0) )
               yc = 2.0 * y1540 * (x -1.0) ! y along crystallization curve
               y = y40 - (y40 - yc) * (0.40-aw) / (0.40 - awc)
            end if ! end of checking for aw
          end if ! end of checking on irh

       else if( x .lt. 1.9999) then
c
           y= 0.0
           if( irh .ge. 40) then 
             mfs15 = poly4(c15,aw)
             mfs2  = poly4(c2,aw)
             y15 = (1.0 - mfs15) / mfs15
             y2  = (1.0 - mfs2) / mfs2
             y = 2.0 * (y15 * (2.0 - x) + y2 * (x - 1.5) )
           end if ! end of check for crystallization
c
c
c       
c   
      else ! 1.9999 < x 
        
c regime where ammonium sulfate and ammonium nitrate are in solution.
c      
c *** following cf&s for both ammonium sulfate and ammonium nitrate
c *** check for crystallization here. their data indicate a 40% value
c     is appropriate.          
            y2 = 0.0
            y3 = 0.0
            if ( irh .ge. 40) then
              mfsSO4 = poly6(kSO4,aw)
              mfsNO3 = poly6(kNO3,aw)         
              y2 = (1.0 - mfsSO4) / mfsSO4
              y3 = (1.0 - mfsNO3) / mfsNO3   
                         
            end if
c       
       end if ! end of checking on x 
c
c *** now set up output of wh2o

c      wh2o units are micrograms (liquid water) / cubic meter of air       
c
       if ( x .lt. 1.9999) then
       
         wh2o =  y * (tso4 * mwso4 + mwnh4 * tnh4)
         
       else
       
c *** this is the case that all the sulfate is ammonium sulfate
c     and the excess ammonium forms ammonum nitrate
        
        wh2o =   y2 * tso4 * mw2 + y3 * tno3 * mwano3 
       
       end if 
c
       return
       end

c23456789012345678901234567890123456789012345678901234567890123456789012     
      
      function poly4(A,X)
      real poly4
      real A(4), X     
       poly4 = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) )))
      return
      end

      function poly6(A,X)
      real poly6
      real A(6), X 
      poly6 = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) + 
     &           X * ( A(5) + X * (A(6)  )))))
       return
      end      ! awater
c //////////////////////////////////////////////////////////////////
