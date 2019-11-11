!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!####################
MODULE MODI_GAMMA_INC_LOW
!####################
!
INTERFACE
!
FUNCTION GAMMA_INC_LOW(PA,PX)  RESULT(PGAMMA_INC_LOW)
REAL, INTENT(IN)                                  :: PA
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PGAMMA_INC_LOW
END FUNCTION GAMMA_INC_LOW
!
END INTERFACE
!
END MODULE MODI_GAMMA_INC_LOW
!     #############################################
      FUNCTION GAMMA_INC_LOW(PA,PX)  RESULT(PGAMMA_INC_LOW)
!     #############################################
!     
!
!!****  *GAMMA_INC_LOW * -  Generalized gamma  function  
!!                   
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the generalized 
!!   lower incomplete Gamma function of its argument.
!!
!!                             /X
!!                             |
!!    GAMMA_INC_LOW(A,X)= ---- | Z**(A-1) EXP(-Z) dZ with A >0
!!                             | 
!!                             /0
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODULE MODI_GAMMA : computation of the Gamma function
!!
!!    REFERENCE
!!    ---------
!!    U. Blahak : Efficient approximation of the incomplete gamma function for
!!                use in cloud model applications, GMD, 2010  
!!
!!
!!    AUTHOR
!!    ------
!!	   V. Vionnet (CNRM/GMME)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     20/09/10
!
!*       0. DECLARATIONS
!           ------------
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     :: PA
REAL, INTENT(IN)                     :: PX
REAL                                 :: PGAMMA_INC_LOW
!
!*       0.2 declarations of local variables
!
REAL                                 :: ZP(6), ZQ(4), ZR(4), ZS(5)
REAL                                 :: ZC(4)
REAL                                 :: ZWORK
!
!*       0.3 initializations of local variables
!
ZP(1) = 9.4368392235E-3
ZP(2) = -1.0782666481E-4
ZP(3) = -5.8969657295E-6
ZP(4) = 2.8939523781E-7
ZP(5) = 1.0043326298E-1
ZP(6) = 5.5637848465E-1

ZQ(1) = 1.1464706419E-1
ZQ(2) = 2.6963429121
ZQ(3) = -2.9647038257
ZQ(4) = 2.1080724954

ZR(1) = 0.0
ZR(2) = 1.1428716184
ZR(3) = -6.6981186438E-3
ZR(4) = 1.0480765092E-4

ZS(1) = 1.0356711153
ZS(2) = 2.3423452308
ZS(3) = -3.6174503174E-1
ZS(4) = -3.1376557650
ZS(5) = 2.9092306039
!
!*       1 Compute coefficients
!
IF( (PX.LT.0.0).OR.(PA.LE.0.0) ) THEN
  PRINT *,' BAD ARGUMENTS IN GAMMA_INC_LOW'
!callabortstop
CALL ABORT
  STOP
END IF
!
!
ZC(1) = 1.+ZP(1)*PA+ZP(2)*PA**2+ZP(3)*PA**3+ZP(4)*PA**4+ZP(5)*(EXP(-ZP(6)*PA)-1)
!
ZC(2) = ZQ(1) + ZQ(2)/PA + ZQ(3)/PA**2 + ZQ(4)/PA**3
!
ZC(3) = ZR(1)+ZR(2)*PA+ZR(3)*PA**2+ZR(4)*PA**3
!
ZC(4) = ZS(1) + ZS(2)/PA + ZS(3)/PA**2 + ZS(4)/PA**3 + ZS(5)/PA**4
!
!*       2 Compute final results
!
ZWORK = 0.5+0.5*TANH(ZC(2)*(PX-ZC(3)))

PGAMMA_INC_LOW = EXP(-PX)* PX**PA * (1./PA +ZC(1)*PX/(PA*(PA+1.))+(ZC(1)*PX)**2/(PA*(PA+1.)*(PA+2.))) &
                            * (1.-ZWORK) + GAMMA(PA)*ZWORK*(1.-ZC(4)**(-PX)) 
RETURN
!
END FUNCTION GAMMA_INC_LOW
