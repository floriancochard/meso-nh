C         Include file : unifacparam.h
C
C         Purpose: Unifac parameters
C
C         Included in Unidrivb.f
C
CNotes: 8 primary compounds: 
C             nC29
C             C4 diacid
C             naphthalene - 2,6-diacid
C             benzo(ghi) perylene
C             biomarker compound
C             4-carboxybenzoic acid
C             C18 acid
C             1,5-ditertbutyl-4,8-di(dimethylpropyl)-decalin
C       + 10 secondary compounds at 25C
C       (Data of Lyman, Reehl, Rosenblatt, 1990)
C
C       NMOL and NFUNC need to match DIMMOL and DIMFUN in unidrivb.f
C
C revision History:  1. Developed by Betty Pun, AER, December, 1999 
C                 under CARB funding
C                   2. Modified by Betty Pun, AER, December, 1999 
C                      under CARB funding for Type B module
C                   3. updated for new POA mix and proper structures
C                      and for use in fortran, Rob Griffin, CIT 1/00
C
          REAL RG(16), QG(16), A(16,16)
          INTEGER NU(18,16)

C
C
C    no. of possible molecules      
c values from Rob          NMOL = 18
C orilam have 10 SOA + 1 POA
         NMOL = 11
C
C
C   no. of functional groups
          NFUNC = 16
C
C
      DIMFUN = 16
C          max dimension of unifac parameter arrays for functional groups 
C          = NFUNC in unifacparamb.h
C 
      DIMMOL = 18
C          max dimension of unifac parameter arrays for functional groups 
C          = NMOL in unifacparamb.h
C    Z = 10 is a fixed parameter in Unifac
C
          Z = 10.0
C
         TINY = 1.0E-08
C
C       order for functional groups: CH3, CH2, CH, C, C=C, C=CH
C       aroC, pahC, tol, ethylbenz, OH, phen, ketone, alde, COOH, NO2
C
C     group volume parameters
C     dimension of RG is the same as NFUNC
          DATA RG  /0.9011, 0.6744, 0.4469, 0.2195, 0.6605, 0.8886,
     +              0.5313, 0.3562, 1.2663, 1.0396, 1.0000, 0.8952,
     +              1.6724, 0.9980, 1.3013, 1.4199/
C
C    group surface area parameters
C    dimension of QG is the same as NFUNC
          DATA QG  /0.8480, 0.5400, 0.2280, 0.0, 0.4850, 0.6760,
     +              0.4000, 0.1200, 0.9680, 0.6600, 1.2000, 0.6800,
     +              1.4880, 0.9480, 1.2440, 1.1040/
C
C   no. of groups in each molecule (# CH3 in each 1-18, then
C    # CH2 in each 1-18, then # CH in each 1-18, etc.)
C
c     start with CH3
c
      DATA  NU  /2, 0, 0, 0, 8, 0, 1, 12, 0, 0, 0, 2, 3,
     +      0, 1, 2, 2, 2,
c        CH2
     +      27, 2, 0, 0, 11, 0, 16, 6, 0, 0, 0, 12, 3, 0,
     +      0, 0, 0, 2,
c        CH  
     +      0, 0, 0, 0, 6, 0, 0, 6, 0, 0, 0, 1, 2, 0, 0, 0,
     +      2, 3,
c        c
     +      0, 0, 0, 0, 5, 0, 0, 4, 0, 0, 0, 0, 1, 0, 0, 0,
     +      0, 0,
c        C=C
     +      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     +      0, 0,
c         C=H
     +      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1,
     +      1, 0,
C       aroC
     +      0, 0, 8, 12, 0, 6, 0, 0, 2, 3, 7, 0, 0, 0, 0, 0,
     +      0, 0,
C       pahC
     +      0, 0, 2, 10, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0,
     +      0, 0,
c        tol
     +      0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0,
     +      0, 0,
c         EB
     +      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     +      0, 0,
c       OH
     +      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1,
     +      1, 1,
c      phen
     +      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     +      0, 0,
c      ket
     +      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     +      1, 1,
c     ald
     +      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2,
     +      0, 1,
c      COOH
     +      0, 2, 2, 0, 0, 2, 1, 0, 1, 1, 0, 0, 0, 2, 2, 0,
     +      1, 0,
c      no2
     +      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0,
     +      0, 0/
C
c   CH3 with 1-16 then CH2 with 1-16 then CH with 1,16
c   then C with 1-16, etc.
c
c    start with CH3
c
C  interaction parameters group-group
      DATA A / 0.0, 0.0, 0.0, 0.0, 2520.00, 2520.00, -11.1200, 
     +       -11.1200, -69.700, -69.700, 156.400, 10000.00,
     +       26.7600, 505.700, 315.300, 5541.00,
c       CH2
     +       0.0, 0.0, 0.0,
     +       0.0, 2520.00, 2520.00, -11.1200, -11.1200, -69.700,
     +       -69.700, 156.400, 10000.00, 26.7600, 505.700,
     +       315.300, 5541.00,
c       CH
     +       0.0, 0.0, 0.0, 0.0, 2520.00,
     +       2520.00, -11.1200, -11.1200, -69.700, -69.700,
     +       156.400, 10000.00, 26.7600, 505.700, 315.300,
     +       5541.00,
c        C
     +       0.0, 0.0, 0.0, 0.0, 2520.00, 2520.00,
     +       -11.1200, -11.1200, -69.700, -69.700, 156.400,
     +       10000.00, 26.7600, 505.700, 315.300, 5541.00, 
c        C=C
     +       -200.00, -200.00, -200.00, -200.00, 0.0, 0.0,
     +       -97.7800, -97.7800, -269.700, -269.700, 8694.00,
     +       732.200, -82.9200, 0.0, 349.200, 0.0,
c        C=CH
     +       -200.00,
     +       -200.00, -200.00, -200.00, 0.0, 0.0, -97.7800,
     +       -97.7800, -269.700, -269.700, 8694.00, 732.200,
     +       -82.9200, 0.0, 349.200, 0.0,
c        aro
     +       61.1300, 61.1300,
     +       61.1300, 61.1300, 340.700, 340.700, 0.0, 0.0,
     +       -146.800, -146.800, 89.600, 270.200, 140.100, 0.0,
     +       62.3200, 1824.00,
c        pah
     +       61.1300, 61.1300, 61.1300, 61.1300,
     +       340.700, 340.700, 0.0, 0.0, -146.800, -146.800, 89.600,
     +       270.200, 140.100, 0.0, 62.3200, 1824.00,
c        tol
     +       76.500,
     +       76.500, 76.500, 76.500, 4102.00, 4102.00, 167.00,
     +       167.00, 0.0, 0.0, 25.8200, 10000.00, 365.800, 0.0,
     +       268.200, -127.800,
c        EB
     +       76.500, 76.500, 76.500, 76.500,
     +       4102.00, 4102.00, 167.00, 167.00, 0.0, 0.0, 25.8200,
     +       10000.00, 365.800, 0.0, 268.200, -127.800,
c        OH
     +       986.500, 986.500, 986.500, 986.500, 693.900, 693.900,
     +       636.100, 636.100, 803.200, 803.200, 0.0, -274.500, 
     +       164.500, -404.800, -151.00, 561.600,
c       phenol
     +       912.200, 912.200,
     +       912.200, 912.200, 926.300, 926.300, 1174.00, 1174.00,
     +       674.300, 674.300, -442.100, 0.0, -246.800, 0.0, 0.0,
     +       0.0, 
c        ket
     +       476.400, 476.400, 476.400, 476.400, 524.500,
     +       524.500, 25.7700, 25.7700, -52.100, -52.100, 84.00,
     +       -158.800, 0.0, 128.00, -297.800, 0.0,
c        ald
     +       677.00, 677.00,
     +       677.00, 677.00, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 441.800,
     +       0.0, -37.3600, 0.0, 0.0, 0.0,
c       COOH
     +       663.500, 663.500, 663.500,
     +       663.500, 730.400, 730.400, 537.400, 537.400, 603.800,
     +       603.800, 119.00, 0.0, 669.400, 0.0, 0.0, 0.0,
c       C-NO2
     +       543.00,
     +       543.00, 543.00, 543.00, 0.0, 0.0, 194.900, 194.900, 
     +       4448.00, 4448.00, 157.100, 0.0, 0.0, 0.0, 0.0, 0.0/
C
