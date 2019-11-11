!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ################
      MODULE MODI_TROE
!     ################

INTERFACE TROE
   MODULE PROCEDURE TROE8, TROE9
END INTERFACE TROE

CONTAINS 
!
!========================================================================
!!    ##################################################################
      FUNCTION TROE8(PCOEF, PKO, PNEXP, PKINF, PMEXP, PM, PT, KVECNPT)
!!    ##################################################################
!!
!!*** *TROE*
!!
!!    PURPOSE
!!    -------
!     this function implements the TROE reaction rate for RACM
!!
!!    REFERENCE
!!    ---------
!!    Stockwell et al., JGR, 1997
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original 27/01/98
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: TROE8
REAL,                     INTENT(IN) :: PCOEF,PKO, PNEXP, PKINF, PMEXP
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT
!!
!!    LOCAL VARIABLES
!!    ---------------
REAL, DIMENSION(KVECNPT) :: ZKOTM, ZKINFT, ZFACT
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
!*        1. THE EXPRESSION
!         -----------------
!
TROE8(:) = TROE9(PCOEF, PKO, PNEXP, PKINF, PMEXP, 0.6, PM, PT, KVECNPT)
!
END FUNCTION TROE8
!
!!    ###################################################################
      FUNCTION TROE9(PCOEF, PKO, PNEXP, PKINF, PMEXP, PF, PM, PT, KVECNPT) 
!!    ###################################################################
!!
!!*** *TROE*
!!
!!    PURPOSE
!!    -------
!     this function implements the TROE reaction rate for RACM
!!
!!    REFERENCE
!!    ---------
!!    Stockwell et al., JGR, 1997
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original 27/01/98
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask 
REAL, DIMENSION(KVECNPT)             :: TROE9
REAL,                     INTENT(IN) :: PCOEF, PKO, PNEXP, PKINF, PMEXP, PF 
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT 
!!
!!    LOCAL VARIABLES
!!    ---------------
REAL, DIMENSION(KVECNPT) :: ZKOTM, ZKINFT, ZFACT
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
!*        1. THE EXPRESSION
!         -----------------
!
ZKOTM(:) = PM(:) * PKO * ( (PT(:)/300.)**(-PNEXP) ) 
ZKINFT(:)= PKINF * ( (PT(:)/300.)**(-PMEXP) ) 
ZFACT(:) = PF**(1./(1.+ALOG10(ZKOTM(:)/ZKINFT(:))**2 )) 
TROE9(:) = PCOEF*(ZKOTM(:)/(1.+ZKOTM(:)/ZKINFT(:)))*ZFACT(:) 
!
END FUNCTION TROE9
!
END MODULE MODI_TROE
!

!     ######################
      MODULE MODI_TROE_EQUIL
!     ######################
INTERFACE TROE_EQUIL
   MODULE PROCEDURE TROE_EQUIL9, TROE_EQUIL10
END INTERFACE TROE_EQUIL

CONTAINS 

!!    #########################################################################
      FUNCTION TROE_EQUIL9(PKO, PNEXP, PKINF, PMEXP, PAFACT, PB, PM, PT, KVECNPT)
!!    #########################################################################
!!
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL,DIMENSION(KVECNPT)              :: TROE_EQUIL9
REAL,                     INTENT(IN) :: PKO, PNEXP, PKINF, PMEXP, PAFACT, PB
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
!*        1. THE EXPRESSION
!         -----------------
!
TROE_EQUIL9(:) = TROE_EQUIL10(PKO, PNEXP, PKINF, PMEXP, PAFACT, PB, 0.6, PM, PT, KVECNPT) 

END FUNCTION TROE_EQUIL9
!
!!    ############################################################################## 
      FUNCTION TROE_EQUIL10(PKO, PNEXP, PKINF, PMEXP, PAFACT, PB, PF, PM, PT, KVECNPT) 
!!    ############################################################################## 
!! 
!!*** *TROE_EQUIL* 
!! 
!!    PURPOSE 
!!    ------- 
!     this function implements the TROE_EQUIL reaction rate for RACM 
!!    Formulation is from RACM and PKO, PNEXP, PKINF, PMEXP, PF  
!!    factors are issued from CACM, 
!!    whereas PAFACT, PB factors are from RACM. 
!! 
!!    REFERENCE 
!!    --------- 
!!    Stockwell et al., JGR, 1997 
!! 
!!    AUTHOR 
!!    ------ 
!!    Karsten Suhre (LA) 
!!     
!!    MODIFICATIONS 
!!    ------------- 
!!    Original 27/01/98 
!!    Tulet P. 08/11/04  Update to CACM 
!! 
!!------------------------------------------------------------------------------ 
!! 
!!    EXTERNAL 
!!    -------- 
!!    none 
!! 
!!    IMPLICIT ARGUMENTS 
!!    ------------------ 
!!    none 
!! 
!!    EXPLICIT ARGUMENTS 
!!    ------------------ 
IMPLICIT NONE 
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask 
REAL,DIMENSION(KVECNPT)              :: TROE_EQUIL10
REAL,                     INTENT(IN) :: PKO, PNEXP, PKINF, PMEXP, PAFACT, PB, PF 
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT 
!! 
!!    LOCAL VARIABLES 
!!    --------------- 
REAL, DIMENSION(KVECNPT) :: ZKOTM, ZKINFT, ZFACT 
!! 
!------------------------------------------------------------------------------ 
!! 
!!    EXECUTABLE STATEMENTS 
!!    --------------------- 
! 
!*        1. THE EXPRESSION 
!         ----------------- 
! 
ZKOTM(:)         = PM(:) * PKO * ( (PT(:)/300.)**(-PNEXP) ) 
ZKINFT(:)        = PKINF * ( (PT(:)/300.)**(-PMEXP) ) 
ZFACT(:)         = PF**(1./(1.+ALOG10(ZKOTM(:)/ZKINFT(:))**2 )) 
TROE_EQUIL10(:) = PAFACT*exp(-PB/PT(:))*(ZKOTM(:) & 
               /(1.+ZKOTM(:)/ZKINFT(:)))*ZFACT(:) 
!
END FUNCTION TROE_EQUIL10

END MODULE MODI_TROE_EQUIL

!
!======================================================================== 
! 
!     ###################### 
      MODULE MODI_TROE_KA 
!     ###################### 
INTERFACE 
FUNCTION TROE_KA(PK1, PK1EXP, PK2, PK2EXP, PK3, PK3EXP, PM, PT,KVECNPT) 
IMPLICIT NONE 
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask 
REAL,DIMENSION(KVECNPT)              :: TROE_KA 
REAL,                     INTENT(IN) :: PK1, PK1EXP, PK2, PK2EXP, PK3, PK3EXP  
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT 
END FUNCTION TROE_KA 
END INTERFACE 
END MODULE MODI_TROE_KA 
! 
!!    #######################################################################
      FUNCTION TROE_KA(PK1, PK1EXP, PK2, PK2EXP, PK3, PK3EXP, PM, PT,KVECNPT) 
!!    #######################################################################
!! 
!!*** *TROE_EQUIL* 
!! 
!!    PURPOSE 
!!    ------- 
!     this function implements rate constant Calculation for CACM 
!! 
!!    REFERENCE 
!!    --------- 
!!    Griffin et al., JGR, 2002 
!! 
!!    AUTHOR 
!!    ------ 
!!    Tulet (LA) 
!!     
!!    MODIFICATIONS 
!!    ------------- 
!!    Original 08/11/04 
!! 
!!------------------------------------------------------------------------------ 
!! 
!!    EXTERNAL 
!!    -------- 
!!    none 
!! 
!!    IMPLICIT ARGUMENTS 
!!    ------------------ 
!!    none 
!! 
!!    EXPLICIT ARGUMENTS 
!!    ------------------ 
IMPLICIT NONE 
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask 
REAL,DIMENSION(KVECNPT)              :: TROE_KA 
REAL,                     INTENT(IN) :: PK1, PK1EXP, PK2, PK2EXP, PK3, PK3EXP  
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT 
!! 
!!    LOCAL VARIABLES 
!!    --------------- 
REAL, DIMENSION(KVECNPT) :: ZK1, ZK2, ZK3 
!! 
!------------------------------------------------------------------------------ 
!! 
!!    EXECUTABLE STATEMENTS 
!!    --------------------- 
! 
!*        1. THE EXPRESSION 
!         ----------------- 
! 
ZK1(:)        = PK1* EXP(PK1EXP / PT(:)) 
ZK2(:)        = PK2* EXP(PK2EXP / PT(:)) 
ZK3(:)        = PK3* EXP(PK3EXP / PT(:)) 
TROE_KA(:)    = ZK1(:) + ZK3(:)*PM(:)*(1 + ZK3(:)*PM(:) / ZK2(:)) 
! 
END FUNCTION TROE_KA 
!
!======================================================================== 
! 
!     ###################### 
      MODULE MODI_TROE_KB 
!     ###################### 
INTERFACE 
FUNCTION TROE_KB(PK1, PK1EXP, PK2, PK2EXP, PM, PT,KVECNPT) 
IMPLICIT NONE 
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask 
REAL,DIMENSION(KVECNPT)              :: TROE_KB 
REAL,                     INTENT(IN) :: PK1, PK1EXP, PK2, PK2EXP 
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT 
END FUNCTION TROE_KB 
END INTERFACE 
END MODULE MODI_TROE_KB 
! 
!!    ########################################################## 
      FUNCTION TROE_KB(PK1, PK1EXP, PK2, PK2EXP, PM, PT,KVECNPT) 
!!    ########################################################## 
!! 
!!*** *TROE_EQUIL* 
!! 
!!    PURPOSE 
!!    ------- 
!     this function implements rate constant Calculation for CACM 
!! 
!!    REFERENCE 
!!    --------- 
!!    Griffin et al., JGR, 2002 
!! 
!!    AUTHOR 
!!    ------ 
!!    Tulet (LA) 
!!     
!!    MODIFICATIONS 
!!    ------------- 
!!    Original 08/11/04 
!! 
!!------------------------------------------------------------------------------ 
!! 
!!    EXTERNAL 
!!    -------- 
!!    none 
!! 
!!    IMPLICIT ARGUMENTS 
!!    ------------------ 
!!    none 
!! 
!!    EXPLICIT ARGUMENTS 
!!    ------------------ 
IMPLICIT NONE 
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask 
REAL,DIMENSION(KVECNPT)              :: TROE_KB 
REAL,                     INTENT(IN) :: PK1, PK1EXP, PK2, PK2EXP 
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PM, PT 
!! 
!!    LOCAL VARIABLES 
!!    --------------- 
REAL, DIMENSION(KVECNPT) :: ZK1, ZK2 
!! 
!------------------------------------------------------------------------------ 
!! 
!!    EXECUTABLE STATEMENTS 
!!    --------------------- 
! 
!*        1. THE EXPRESSION 
!         ----------------- 
! 
ZK1(:)        = PK1* EXP(PK1EXP / PT(:)) 
ZK2(:)        = PK2* EXP(PK2EXP / PT(:)) 
TROE_KB(:)    = ZK1(:) + ZK2(:)*PM(:) 
! 
END FUNCTION TROE_KB 
!
!========================================================================
!
!     ##############
      MODULE MODI_KT
!     ##############
INTERFACE
FUNCTION KT(PALPHA, PMM, PT, PRAD, KVECNPT)
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: KT
REAL,                     INTENT(IN) :: PALPHA,PMM
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT, PRAD
END FUNCTION KT
END INTERFACE
END MODULE MODI_KT
!
!
!!    ###########################################
      FUNCTION KT(PALPHA, PMM, PT, PRAD, KVECNPT)
!!    ###########################################
!!
!!*** *KT*
!!
!!    PURPOSE
!!    -------
!     this function implements the mass transfer reaction rate
!     for exchange between gas and aqueous phase
!!
!!    REFERENCE
!!    ---------
!!    Schwartz, 1986
!!
!!    AUTHORS
!!    ------
!!    Céline Mari & Maud Leriche (LA)
!! 
!!    MODIFICATIONS
!!    -------------
!!    Original 22/02/2007
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: KT
REAL,                     INTENT(IN) :: PALPHA,PMM
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT, PRAD
!!
!!    LOCAL VARIABLES
!!    ---------------
REAL, DIMENSION(KVECNPT)             :: ZDG, ZV
!
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
!*        1. THE EXPRESSION
!         -----------------
!
ZV(:) = SQRT(8.*PT(:)*8.3144/(PMM*1.e-3)/3.1415926535898)
ZDG(:) = 1.e-7 * ZV(:) /3.    !gas-phase diffusion
KT(:) = 1./(PRAD(:)*PRAD(:)/(3.*ZDG(:))+(4.*PRAD(:)/(3.*ZV(:)*PALPHA)))
!
END FUNCTION KT
!
!========================================================================
!
!     #################
      MODULE MODI_HENRY
!     #################
INTERFACE
FUNCTION HENRY(PH, PDEPT, PT, KVECNPT)
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: HENRY
REAL,                     INTENT(IN) :: PH,PDEPT
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT
END FUNCTION HENRY
END INTERFACE
END MODULE MODI_HENRY
!
!
!!    ###########################################
      FUNCTION HENRY(PH, PDEPT, PT, KVECNPT)
!!    ###########################################
!!
!!*** *HENRY*
!!
!!    PURPOSE
!!    -------
!     this function computes the henry's law constant at a given temperature
!!
!!    AUTHOR
!!    ------
!!    Maud Leriche (LA)
!! 
!!    MODIFICATIONS
!!    -------------
!!    Original 23/02/2007
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: HENRY
REAL,                     INTENT(IN) :: PH,PDEPT
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
!*        1. THE EXPRESSION
!         -----------------
!
HENRY(:) = PH*EXP(-PDEPT*((1./PT(:))-(1./298.15)))
!
END FUNCTION HENRY
!
!========================================================================
!
!     #################
      MODULE MODI_HEFFA
!     #################
INTERFACE
FUNCTION HEFFA(PH, PDEPT, PK1, PDEPT1, PK2, PDEPT2, PPH, PT, KVECNPT)
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: HEFFA
REAL,                     INTENT(IN) :: PH,PDEPT,PK1,PDEPT1,PK2,PDEPT2
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT,PPH
END FUNCTION HEFFA
END INTERFACE
END MODULE MODI_HEFFA
!
!
!!    #####################################################################
      FUNCTION HEFFA(PH, PDEPT, PK1, PDEPT1, PK2, PDEPT2, PPH, PT, KVECNPT)
!!    #####################################################################
!!
!!*** *HEFFA*
!!
!!    PURPOSE
!!    -------
!     this function computes the effective henry's law constant
!     at a given temperature for acid
!!
!!    AUTHOR
!!    ------
!!    Maud Leriche (LA)
!! 
!!    MODIFICATIONS
!!    -------------
!!    Original 23/02/2007
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: HEFFA
REAL,                     INTENT(IN) :: PH,PDEPT,PK1,PDEPT1,PK2,PDEPT2
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT,PPH
!!
!!    LOCAL VARIABLES
!!    ---------------
REAL, DIMENSION(KVECNPT)             :: ZHPLUS, ZH, ZK1, ZK2
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
!*        0. Concentration of H+
!         -------------------------
!
ZHPLUS(:) = 10.**(-PPH(:))
!
!*        1. Temperature dependency
!         -------------------------
!
ZH(:) = PH*EXP(-PDEPT*((1./PT(:))-(1./298.15)))
ZK1(:) = PK1*EXP(-PDEPT1*((1./PT(:))-(1./298.15)))
ZK2(:) = PK2*EXP(-PDEPT2*((1./PT(:))-(1./298.15)))
!
!*        2. THE EXPRESSION
!         -----------------
!
HEFFA(:) = ZH(:)*( 1.0 + (ZK1(:)/ZHPLUS(:))*(1.0+(ZK2(:)/(ZHPLUS(:)))) )
!
END FUNCTION HEFFA
!
!========================================================================
!
!     #################
      MODULE MODI_HEFFB
!     #################
INTERFACE
FUNCTION HEFFB(PH, PDEPT, PK1, PDEPT1, PPH, PT, KVECNPT)
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: HEFFB
REAL,                     INTENT(IN) :: PH,PDEPT,PK1,PDEPT1
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT,PPH
END FUNCTION HEFFB
END INTERFACE
END MODULE MODI_HEFFB
!
!
!!    ########################################################
      FUNCTION HEFFB(PH, PDEPT, PK1, PDEPT1, PPH, PT, KVECNPT)
!!    ########################################################
!!
!!*** *HEFFB*
!!
!!    PURPOSE
!!    -------
!     this function computes the effective henry's law constant
!     at a given temperature for base
!!
!!    AUTHOR
!!    ------
!!    Maud Leriche (LA)
!! 
!!    MODIFICATIONS
!!    -------------
!!    Original 23/02/2007
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER,                  INTENT(IN) :: KVECNPT ! no. of points in vector mask
REAL, DIMENSION(KVECNPT)             :: HEFFB
REAL,                     INTENT(IN) :: PH,PDEPT,PK1,PDEPT1
REAL, DIMENSION(KVECNPT), INTENT(IN) :: PT,PPH
!!
!!    LOCAL VARIABLES
!!    ---------------
REAL, DIMENSION(KVECNPT)             :: ZKH2O, ZHPLUS, ZOHM, ZH, ZK1
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
!*        0. Concentration of OH-
!         -------------------------
!
ZHPLUS(:) = 10.**(-PPH(:))
ZKH2O = 1.e-14*EXP(-6716*((1./PT(:))-(1./298.15)))
ZOHM(:) = ZKH2O(:)/ZHPLUS(:)
!
!*        1. Temperature dependency
!         -------------------------
!
ZH(:) = PH*EXP(-PDEPT*((1./PT(:))-(1./298.15)))
ZK1(:) = PK1*EXP(-PDEPT1*((1./PT(:))-(1./298.15)))
!
!*        2. THE EXPRESSION
!         -----------------
!
HEFFB(:) = ZH(:)*( 1.0 + ZK1(:)/ZOHM(:) )
!
END FUNCTION HEFFB
!


