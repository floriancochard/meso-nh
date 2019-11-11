!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ######spl
      MODULE MODE_FGAU
!     ####################
!
!!****  *MODE_FGAU* -  module routines  
!!
!!    PURPOSE
!!    -------
!!      Compute some Gaussian quadrature nodes and weights and the factorial of an
!!    integer. --- they might need some optimization, idea: look at newer Numerical
!!    Recipes  ---
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST   
!!        XPI
!!
!!    REFERENCE
!!    ---------
!!      Press, W. H., B. P. Flannery, S. A. Teukolsky et W. T. Vetterling, 1986: 
!!    Numerical Recipes: The Art of Scientific Computing. Cambridge University 
!!    Press, 818 pp.
!!
!!    AUTHOR
!!    ------
!!      O. Caumont & V. Ducrocq  * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   26/03/2004    
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS

!              ------------
USE MODD_CST, ONLY: XPI
!
!-------------------------------------------------------------------------------
!
CONTAINS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       2.   ROUTINE GAULEG
!             -------------------
!-------------------------------------------------------------------------------
!     #########################
      SUBROUTINE GAULEG(N,X,W)
!     ########################
! computes positive Gauss-Legendre nodes and weights over ]-1,1[
! adapted from Press et al. 1986
    IMPLICIT NONE
    REAL,PARAMETER :: EPS=1E-6
    INTEGER,INTENT(IN) :: N
    REAL,DIMENSION((N+1)/2),INTENT(OUT) :: X,W
    INTEGER :: I,J
    REAL :: P1,P2,P3,PP,Z,Z1
    
    DO I=1,(N+1)/2
       Z=COS(XPI*(I-.25)/(N+.5))
       Z1=Z+1.
       
       DO WHILE(ABS(Z-Z1) > EPS)
          P1=1.
          P2=0.
          DO J=1,N
             P3=P2
             P2=P1
             P1=((2.*J-1.)*Z*P2-(J-1.)*P3)/J
          END DO
          
          PP=N*(Z*P1-P2)/(Z*Z-1.)
          Z1=Z
          Z=Z1-P1/PP
       END DO
       X((N+1)/2-I+1)=Z
       W((N+1)/2-I+1)=2./(1-Z*Z)/PP/PP
    END DO
    DO I=1,(N+1)/2
    END DO
  END SUBROUTINE GAULEG

!-------------------------------------------------------------------------------
!
!*       2.   ROUTINE GAULAG
!             -------------------
!-------------------------------------------------------------------------------
!     #########################
      SUBROUTINE GAULAG(N,X,W)
!     #########################
! returns nodes and weights of Gauss-Laguerre quadrature. (From Press et al.)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL, DIMENSION(N), INTENT(OUT) :: X,W
    REAL, PARAMETER :: EPS=3E-13

    REAL :: ALF=0.
    INTEGER :: ITS,J
    INTEGER,PARAMETER :: MAXIT=10
    REAL :: ANU
    REAL,PARAMETER :: C1=9.084064E-1,C2=5.214976E-2,C3=2.579930E-3,C4=3.986126e-3
    REAL,DIMENSION(N) :: RHS,R2,R3,THETA,P1,P2,P3,PP,Z,Z1
    LOGICAL,DIMENSION(N) :: UNFINISHED

!    WRITE(0,*) 'IN GAULAG'
    ANU=4.*N+2.*ALF+2.
    RHS=ARTH(4.*N-1.,-4.,N)*XPI/ANU
    R3=RHS**(1./3.)
    R2=R3**2
    THETA=R3*(C1+R2*(C2+R2*(C3+R2*C4)))
    Z=ANU*COS(THETA)**2
    UNFINISHED=.TRUE.
    DO ITS=1,MAXIT
       WHERE(UNFINISHED)
          P1=1.
          P2=0.
       END WHERE
       DO J=1,N
          WHERE(UNFINISHED)
             P3=P2
             P2=P1
             P1=((2.*J-1.+ALF-Z)*P2-(J-1.+ALF)*P3)/J
          END WHERE
       END DO
       WHERE(UNFINISHED)
          PP=(N*P1-(N+ALF)*P2)/Z
          Z1=Z
          Z=Z1-P1/PP
          UNFINISHED=(ABS(Z-Z1)>EPS*Z)
       END WHERE
       IF(.NOT.ANY(UNFINISHED)) EXIT
    END DO
    IF(ITS==MAXIT+1) PRINT*,'ERROR IN GAULAG: TOO MANY ITERATIONS'
    X=Z
    W=-EXP(GAMMLN(ALF+N)-GAMMLN(REAL(N)))/(PP*N*P2)
    
  END SUBROUTINE GAULAG
!-------------------------------------------------------------------------------
!
!*       3.   ROUTINE GAUHER
!             -------------------
!-------------------------------------------------------------------------------
!     #########################
      SUBROUTINE GAUHER(N,X2,W)
!     #########################
!   returns POSITIVE nodes and weights of Gauss-Hermite quadrature.
    IMPLICIT NONE
    ! N : ordre du polynôme de Hermite
    ! X2 : abscisses POSITIVES de la quadrature
    ! W : poids correspondants de la quadrature
    INTEGER,INTENT(IN) :: N
    REAL,DIMENSION((N+1)/2),INTENT(OUT) :: X2,W
    REAL :: EPS=1e-2
    REAL :: PX,DPX,X,Y
    INTEGER,DIMENSION(N+1) :: P0,P1,P2
    REAL,DIMENSION((N+1)/2) :: X1
    
    INTEGER :: I,J,K
    
    IF(N>=15) THEN
       PRINT*,'SUBROUTINE GAUHER FAILS TO CONVERGE FOR N>=15. ANYWAY, THIS NUMBER IS TOO HIGH.'
       PRINT*,'PLEASE TAKE A SMALLER NUMBER OF POINTS OR MODIFY THIS SUBROUTINE.'
       STOP
    END IF

    P0(:)=0
    P1(:)=0
    P2(:)=0
    
    P0(1)=1 ! N=0 H0(x)=1
    P1(1)=0 ! N=1 H1(x)=2x
    P1(2)=2 
    X1(1)=1.
    IF(N==1) THEN
       X2(1)=0.
       W(1)=SQRT(XPI)
       RETURN
    END IF
    DO I=1,N-1
       DO J=0,I+1
          ! CALCUL DU POLY DE HERMITE DE DEGRÉ I+1
          IF(J==0) THEN
             P2(J+1)=-I*P0(J+1)
          ELSE
             P2(J+1)=P1(J)-I*P0(J+1)
          END IF
          P2(J+1)=2*P2(J+1)
       END DO
       P0(:)=P1(:)
       P1(:)=P2(:)
       ! CALCUL DES I+1 ZÉROS
       DO J=1,(I+2)/2
          IF(J==(I+2)/2) THEN
             X=4.*X1((I+1)/2)+1.
          ELSE IF(J==1) THEN
             IF(MOD(I+1,2)==0) THEN
                X=X1(2)/2.
             ELSE
                X=0.
             END IF
          ELSE
             IF(MOD(I+1,2)==0) THEN
                X=(X1(J+1)+X1(J))/2.
             ELSE
                X=(X1(J)+X1(J-1))/2.
             END IF
          END IF
          ! MÉTHODE DE NEWTON
          DO
             ! CALCUL DE P(X) ET P'(X)
             PX=P2(I+2)*X+P2(I+1)
             DPX=P2(I+2)
             DO K=I,1,-1
                DPX=PX+DPX*X
                PX=P2(K)+PX*X
             END DO
             
             Y=X-PX/DPX
             IF(ABS(Y-X)<EPS) EXIT
             X=Y
          END DO
          X2(J)=Y
       END DO
       X1(:)=X2(:)
    END DO
    ! CALCUL DES POIDS
    DO I=1,(N+1)/2
       ! CALCUL DE P'(X(I))
       DPX=N*P2(N+1)
       DO K=N,2,-1
          DPX=(K-1)*P2(K)+DPX*X2(I)
       END DO
       W(I)=2.**(N+1)*Factorial(N)*SQRT(XPI)/DPX**2
    END DO
    RETURN
  END SUBROUTINE GAUHER
!-------------------------------------------------------------------------------
!
!*       4.   FUNCTION Factorial
!             -------------------
!-------------------------------------------------------------------------------
!   ###########################################
    RECURSIVE FUNCTION Factorial(n)  RESULT(Fact)
!   ##########################################
!   returns n!
    IMPLICIT NONE
    INTEGER :: Fact
    INTEGER, INTENT(IN) :: n
    
    IF (n == 0) THEN
       Fact = 1
    ELSE
       Fact = n * Factorial(n-1)
    END IF
    
  END FUNCTION Factorial
!

  FUNCTION ARTH(FIRST,INCREMENT,N)
    REAL,INTENT(IN) :: FIRST,INCREMENT
    INTEGER,INTENT(IN) :: N
    REAL,DIMENSION(N) :: ARTH
    INTEGER :: K
    
    DO K=1,N
       ARTH(K)=FIRST+INCREMENT*(K-1)
    END DO
  END FUNCTION ARTH


  FUNCTION gammln(xx)
    IMPLICIT NONE
    REAL, INTENT(IN) :: xx
    REAL :: gammln
    REAL :: tmp,x
    REAL :: stp = 2.5066282746310005
    REAL, DIMENSION(6) :: coef = (/76.18009172947146,& 
         -86.50532032941677,24.01409824083091,& 
         -1.231739572450155,0.1208650973866179e-2,&
         -0.5395239384953e-5/)
    x=xx
    tmp=x+5.5
    tmp=(x+0.5)*log(tmp)-tmp
    gammln=tmp+log(stp*(1.000000000190015+&
         sum(coef(:)/arth(x+1.,1.,size(coef))))/x)
  END FUNCTION gammln

END MODULE MODE_FGAU
