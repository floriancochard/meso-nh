!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!   #####################
     MODULE MODI_CH_NNARES
!!   #####################
!!
INTERFACE
!!
SUBROUTINE CH_NNARES (PAER,PRH, PDENAIR, PPRESSURE, PTEMP, PRC)
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(INOUT) :: PAER
REAL, DIMENSION(:),   INTENT(IN)    :: PRH, PDENAIR, PPRESSURE, PTEMP, PRC

END SUBROUTINE CH_NNARES
!!
END INTERFACE
!!
END MODULE MODI_CH_NNARES
!!
!!    ##########################################################
      SUBROUTINE CH_NNARES (PAER,PRH, PDENAIR, PPRESSURE, PTEMP, PRC)
!!    ##########################################################
!!
!!    PURPOSE
!!    -------
!!
!!    Calculate the aerosol chemical speciation and water content.
!!    THIS IS NEURAL NET VERSION OF THE ORIGINAL ARES CODE,
!!    GENERATED AND TRAINED USING NNFIT
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    V. Crassier & K. Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!
!*****************************************************************
! Parameters of ARES:
!
!  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
!  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
!  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
!  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
!  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
!  RH    : Fractional relative humidity (IN)
!  TEMP  : Temperature in Kelvin (IN)
!  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
!  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT)
!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
!
!***************************************************************
!!
!!   EXTERNAL
!!   -------
!!
USE MODD_CH_AEROSOL
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
IMPLICIT NONE
!...........ARGUMENTS and their descriptions
REAL, DIMENSION(:,:), INTENT(INOUT) :: PAER
REAL, DIMENSION(:),   INTENT(IN)    :: PRH, PDENAIR, PPRESSURE, PTEMP, PRC
!
! PAER(:,1) :: H2SO4    in micrograms / m**3
! PAER(:,2) :: NH3(g)   in micrograms / m**3
! PAER(:,3) :: HNO3(g)  in micrograms / m**3
! PAER(:,4) :: H2O(a)   in micrograms / m**3
! PAER(:,5) :: NO3(a)   in micrograms / m**3
! PAER(:,6) :: NH4(a)   in micrograms / m**3
!
!...........PARAMETERS and their descriptions:

REAL, PARAMETER ::  ZMWH2O = 18.0           ! molecular weight for water      
REAL, PARAMETER ::  ZMWNO3 = 62.0049        ! molecular weight for NO3
REAL, PARAMETER ::  ZMWHNO3 = 63.01287      ! molecular weight for HNO3
REAL, PARAMETER ::  ZMH2SO4 = 98.07354      ! molecular weight for H2SO4
REAL, PARAMETER ::  ZMWNH3 = 17.03061       ! molecular weight for NH3
REAL, PARAMETER ::  ZMWNH4 = 18.03858       ! molecular weight for NH4

!-----------------------------------------------------------------------------
!     .. the independant variables of the NN are:
!     1) mole fraction of total ammonium           X = (NH4 + NH3) / (SO4)
!     2) mole fraction of total nitrate            Y = (NO3 + HNO3) / (SO4 + NO3 +  HNO3)
!     3) fractional relative humidity              R = RH
!     4) temperature                               T = TEMP
!
!     .. the output variables of the NN are:
!     1) mole fraction of water             A = AH2O / (AH20+ASO4+ANO3+ANH4)
!     2) mole fraction of gas phase ammonia B = GNH3 / (GNH3 + ANH4)
!     3) mole fraction of gas phase nitrate C = GNO3 / (GNO3 + ANO3)
!
REAL, DIMENSION(SIZE(PAER,1)) :: TOTAMM,TOTNIT,TOTSUL
!     .. input ..
REAL, DIMENSION(SIZE(PAER,1)) :: X, Y
!
!     .. output ..
REAL, DIMENSION(SIZE(PAER,1)) :: A, B, C
!
!-----------------------------------------------------------------------------
! print *, "PAER(1:5,1) H2SO4 ", PAER(1:5,1)
! print *, "PAER(1:5,2) NH3   ", PAER(1:5,2)
! print *, "PAER(1:5,3) HNO3  ", PAER(1:5,3)

!     convert to NN input variables X, Y, R, T
      TOTAMM(:) = PAER(:,2)/ZMWNH3 + PAER(:,6)/ZMWNH4
      TOTNIT(:) = PAER(:,3)/ZMWHNO3 + PAER(:,5)/ZMWNO3
      TOTSUL(:) = PAER(:,1)/ZMH2SO4

      X(:) = TOTAMM(:)/ TOTSUL(:)
      Y(:) = TOTNIT(:) / (TOTNIT(:) + TOTSUL(:))


      
!   call the neural net
      CALL NN ( X, Y, PRH, PTEMP, A, B, C, SIZE(PAER,1))

!     compute return variables in mole concentrations
      PAER(:,1) = TOTSUL(:) 
      PAER(:,2) = B(:) * TOTAMM(:) 
      PAER(:,3) = C(:) * TOTNIT(:) 
      PAER(:,5) = (1-C(:)) * TOTNIT(:) 
      PAER(:,6) = (1.-B(:)) * TOTAMM(:) 
      PAER(:,4) = (1./(1.-A(:))) * (PAER(:,1) + PAER(:,5) + PAER(:,6))

!     convert return variables to mass concentration
      PAER(:,1)=PAER(:,1)*ZMH2SO4
      PAER(:,2)=PAER(:,2)*ZMWNH3
      PAER(:,3)=PAER(:,3)*ZMWHNO3
      PAER(:,4)=PAER(:,4)*ZMWH2O
      PAER(:,5)=PAER(:,5)*ZMWNO3
      PAER(:,6)=PAER(:,6)*ZMWNH4

END SUBROUTINE CH_NNARES


      SUBROUTINE NN( X, Y, PRH, TEMP, A, B, C, KVECTN )
!     .. call the neural net ..
USE MODD_CH_AEROSOL
      IMPLICIT NONE

!     .. parameters
INTEGER, INTENT(IN) :: KVECTN
REAL, DIMENSION(KVECTN), INTENT(INOUT) :: X, Y, A, B, C
REAL, DIMENSION(KVECTN), INTENT(IN)    :: PRH, TEMP

REAL, DIMENSION(KVECTN) :: ZRH, ZTEMP
REAL, DIMENSION(KVECTN) :: X1,X2
REAL, DIMENSION(KVECTN) :: A1,A2,B2,C1,C2
REAL,DIMENSION(KVECTN,100) :: XIN,YOUT 


X(:)=MAX(0.1,X(:))
X(:)=MIN(3.999,X(:))
Y(:)=MAX(0.05,Y(:))
Y(:)=MIN(0.99,Y(:))
ZRH(:)=MAX(0.5,PRH(:))
ZRH(:)=MIN(0.9,PRH(:))
ZTEMP(:)=MAX(275.,TEMP(:))
ZTEMP(:)=MIN(315.,TEMP(:))

X2(:)=REAL(INT(X(:)/2.))
X1(:)=REAL(INT(2.-X(:)/2.))

XIN(:,1) = X(:)
XIN(:,2) = Y(:)
XIN(:,3) = ZRH(:)
XIN(:,4) = ZTEMP(:)

CALL FCTRESO(I1IA,J1JA,K1KA,W1IJA,W1JKA,X1MAXA,X1MINA,X1MODA,XIN,YOUT,KVECTN) 
A1(:) = (YOUT(:,1))
CALL FCTRESO(I1IC,J1JC,K1KC,W1IJC,W1JKC,X1MAXC,X1MINC,X1MODC,XIN,YOUT,KVECTN)
C1(:) = (YOUT(:,1))
CALL FCTRESO(I2IA,J2JA,K2KA,W2IJA,W2JKA,X2MAXA,X2MINA,X2MODA,XIN,YOUT,KVECTN) 
A2(:) = (YOUT(:,1))
CALL FCTRESO(I2IB,J2JB,K2KB,W2IJB,W2JKB,X2MAXB,X2MINB,X2MODB,XIN,YOUT,KVECTN) 
B2(:) = exp(YOUT(:,1))
CALL FCTRESO(I2IC,J2JC,K2KC,W2IJC,W2JKC,X2MAXC,X2MINC,X2MODC,XIN,YOUT,KVECTN) 
C2(:) = (YOUT(:,1))

A(:)=X1(:)*A1(:)+X2(:)*A2(:)
B(:)=X2(:)*B2(:)
C(:)=X1(:)*C1(:)+X2(:)*C2(:)

RETURN
END SUBROUTINE NN 

SUBROUTINE FCTRESO(II,JJ,KK,WIJ,WJK,XMAX,XMIN,XMOD,X,Y,KVECTN) 
IMPLICIT NONE
INTEGER, INTENT(IN) :: KVECTN
REAL, DIMENSION(100,100), INTENT(INOUT) :: WIJ,WJK
REAL, DIMENSION(2,100), INTENT(INOUT) :: XMIN,XMAX,XMOD
INTEGER, INTENT(INOUT) :: II,JJ,KK
REAL, DIMENSION(KVECTN,100), INTENT(INOUT) :: X,Y

REAL, DIMENSION(KVECTN,100) :: U,H,S
INTEGER :: I, J, K

!     .. LES EQUATIONS DU RESEAU 
DO I=1,II 
 U(:,I)=(X(:,I)-XMIN(1,I))/(XMAX(1,I)-XMIN(1,I)) 
ENDDO
U(:,II+1)=1. 

DO J=1,JJ 
 H(:,J)=0. 
 DO I=1,II+1 
   H(:,J)=H(:,J)+WIJ(I,J)*U(:,I) 
 ENDDO
 H(:,J)=1./(1.+EXP(-H(:,J))) 
ENDDO
H(:,JJ+1)=1. 

DO K=1,KK 
 S(:,K)=0. 
 DO J=1,JJ+1 
  S(:,K)=S(:,K)+WJK(J,K)*H(:,J) 
 ENDDO
 S(:,K)=1./(1.+EXP(-S(:,K))) 

 Y(:,K)=XMIN(2,K)+S(:,K)*(XMAX(2,K)-XMIN(2,K)) 
ENDDO

RETURN 
END SUBROUTINE FCTRESO
