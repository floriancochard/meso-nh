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
      MODULE MODE_FSCATTER
!     ####################
!
!!****  *MODE_FSCATTER* -  module routines for RADAR_SCATTERING 
!!
!!    PURPOSE
!!    -------
!!       Compute the Pth moment order of the generalized gamma law, dielectric
!!    functions for water and ice, and backscattering and extinction 
!!    efficiencies from Mie theory for a homogeneous or coated sphere.
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE          
!!
!!    REFERENCE
!!    ---------
!!      Bohren, C.F., and D. R. Huffman, 1983 : Absorption and Scattering of 
!!    Light by Small Particles. John Wiley & Sons.
!!
!!
!!    AUTHOR
!!    ------
!!      O. Caumont & V. Ducrocq * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  26/03/2004  
!!      27/05/2014 (O. Caumont) Added Maxwell Garnett equation
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS

!              ------------
!
!-------------------------------------------------------------------------------
!
CONTAINS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       1.   FUNCTION MOMG
!             -------------------
!-------------------------------------------------------------------------------
!     #################################################
       FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
!     #################################################
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA
!
  IMPLICIT NONE
!
  REAL,INTENT(IN)  :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL,INTENT(IN)  :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL,INTENT(IN)  :: PP     ! order of the moment
  REAL             :: PMOMG  ! result: moment of order ZP
!
!------------------------------------------------------------------------------
!
!
  PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
  END FUNCTION MOMG
!
!-------------------------------------------------------------------------------
!
!*       2.   FUNCTION QEPSW
!             -------------------
!-------------------------------------------------------------------------------
!     ###########################################
      FUNCTION QEPSW(PTEMP,PFREQ) RESULT(PQEPSW)
!     ###########################################
   ! water complex dielectric function (Liebe et al., 1991)
    ! electromagnetic fields in exp(-i*omega*t), i.e. Im(epsw)>=0
    ! in  : ptemp=temperature in K
    !       pfreq=frequency in Hz
    ! out : eps=epsilon
    
    IMPLICIT NONE
    REAL,INTENT(IN) :: PFREQ
    REAL,INTENT(IN) :: PTEMP
    REAL :: ZTHETA,ZEZ,ZEINF,ZF
    COMPLEX :: PQEPSW

    ZTHETA=1-300./PTEMP
    ZEZ=77.66-103.3*ZTHETA
    ZEINF=0.066*ZEZ
    ZF=(20.27+146.5*ZTHETA+314.*ZTHETA**2)*1.E9
    
    PQEPSW=ZEINF+(ZEZ-ZEINF)/(1.-(0.,1.)*PFREQ/ZF)

  END FUNCTION QEPSW
!
!-------------------------------------------------------------------------------
!
!*       3.   FUNCTION QEPSI
!             -------------------
!-------------------------------------------------------------------------------
!   ##########################################
    FUNCTION QEPSI(PTEMP,PFREQ) RESULT(PQEPSI)
!   ##########################################
    ! ice complex dielectric function (Hufford, 1991)
    ! electromagnetic fields in exp(-i*omega*t), i.e. Im(epsi)>=0
    ! in  : ptemp=temperature in K
    !       pfreq=frequency in Hz
    ! out : eps=epsilon
    
    IMPLICIT NONE
    REAL,INTENT(IN) :: PFREQ
    REAL,INTENT(IN) :: PTEMP
    REAL :: ZTHETA
    COMPLEX :: PQEPSI

    ZTHETA=1-300./PTEMP
    
    PQEPSI=3.15+(0.,1.)*((50.4-62*ZTHETA)*1.E5*EXP(22.1*ZTHETA)/PFREQ+ &
         ((0.502+0.131*ZTHETA)*1.E-13/(1.-ZTHETA)+ &
         0.542E-15*((1-ZTHETA)/(.0073-ZTHETA))**2)*PFREQ)
!
  END FUNCTION QEPSI
!
!-------------------------------------------------------------------------------
!
!*       4.   SUBROUTINE BHMIE
!             -----------------
!-------------------------------------------------------------------------------
!   ##########################################
    SUBROUTINE BHMIE(X,REFREL,QEXT,QSCA,QBACK)
!   ##########################################
    IMPLICIT NONE
    ! Arguments:
    REAL,INTENT(IN) :: X
    REAL,INTENT(OUT) :: QEXT,QSCA,QBACK
    COMPLEX,INTENT(IN) :: REFREL
    ! Local variables:
    INTEGER :: N,NSTOP
    REAL :: CHI0X,CHI1X,CHIX,DX,CHIPX
    COMPLEX :: Y,XBACK,AN,BN,DY
!
! Modification (C.Lac) 04/2014 : exclude very small values of x 
!
    !         -----------------------------------------------------------
    !              del is the inner sphere convergence criterion
    !         -----------------------------------------------------------
    Y = refrel*x
    NSTOP = X+4.*X**0.3333+2.0
    !         -----------------------------------------------------------
    !              series terminated after nstop terms
    !         -----------------------------------------------------------
    dx = COS(x)/SIN(x)
    DY  = COS(Y)/SIN(Y)
    chi0x = -sin(x)
    chi1x = cos(x)
    qsca = 0.0
    qext = 0.0
    xback = (0.0,0.0)
    n=1
    IF (x <= 1.E-07) THEN
     QEXT = 0.
     QSCA = 0.
     QBACK = 0.
     RETURN
    ELSE
     do while(n<=nstop)
   ! DO n=1,nstop 
       DX  = 1.0/(n/x-dx) - n/x       
       DY  = 1.0/(n/y-dy) - n/y
       chix = (2.0*n-1.)*chi1x/x - chi0x
       chi0x=chi1x
       chi1x=chix
       CHIPX  = CHI0X  - N*CHI1X/X
       an = 1./(1.-(0.,1.)*(CHI1X*DX-CHIPX)*(chi1x*dy-REFREL*chipx)/(dy-refrel*dx))
       bn = 1./(1.-(0.,1.)*(CHI1X*DX-CHIPX)*(REFREL*chi1x*dy-chipx)/(REFREL*dy-dx))
       qsca = qsca + (2.0*n+1.0)*(ABS(an)**2+ABS(bn)**2)
       xback = xback + (2.0*n+1.0)*(-1.)**n*(an-bn)
       qext = qext + (2.0*n + 1.0)*(an+bn)
       n=n+1
     END DO
     QSCA = (2.0/(X*X))*qsca
     QEXT = (2.0/(X*X))*qext
     qback = (1.0/(X*X))*ABS(XBACK)**2
    !qext=4.*x*AIMAG((REFREL**2-1.)/(REFREL**2+2.)*&
    !     (1+X**2/15.*(REFREL**2-1.)/(REFREL**2+2.)*(REFREL**4+27.*REFREL**2+38.)/(2.*REFREL**2+3.)))&
    !     +8./3.*X**4*REAL(((REFREL**2-1.)/(REFREL**2+2.))**2)
    !qback=4.*X**4*ABS((REFREL**2-1)/(REFREL**2+2))**2
    END IF
    RETURN
  END SUBROUTINE BHMIE
!
!-------------------------------------------------------------------------------
!
!*       5.   SUBROUTINE BHCOAT
!             -----------------
!-------------------------------------------------------------------------------
!   ######################################################
    SUBROUTINE BHCOAT(X,Y,RFREL1,RFREL2,QEXT,QSCA,QBACK)
!   ######################################################
    IMPLICIT NONE
! PARAMETERS
    REAL,PARAMETER :: DEL=1.e-8
    COMPLEX,PARAMETER :: II=(0.e0,1.e0)
    ! Arguments:
    REAL,INTENT(IN) :: X,Y
    REAL,INTENT(OUT) :: QEXT,QSCA,QBACK
    COMPLEX,INTENT(IN) :: RFREL1,RFREL2
    ! Local variables:
    LOGICAL :: IFLAG
    INTEGER :: N,NSTOP
    REAL :: CHI0Y,CHI1Y,CHIY,D0Y,D1Y
    COMPLEX :: CHIPX2,CHIPY,CHIPY2,CHIY2,CHIX2,CHI0X2,CHI0Y2,CHI1X2,CHI1Y2
    COMPLEX :: D1X1,D1X2,D1Y2,D0X1,D0X2,D0Y2,BRACK,REFREL,X1,X2,Y2
    COMPLEX :: XBACK,AN,BN,DNBAR,GNBAR,CRACK
    !***********************************************************************
    !
    ! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back for coated sphere.
    ! All bessel functions computed by upward recurrence.
    ! Input:
    !        X = 2*PI*RCORE*REFMED/WAVEL
    !        Y = 2*PI*RMANT*REFMED/WAVEL
    !        RFREL1 = REFCOR/REFMED
    !        RFREL2 = REFMAN/REFMED 
    ! where  REFCOR = complex refr.index of core)
    !        REFMAN = complex refr.index of mantle)
    !        REFMED = real refr.index of medium)
    !        RCORE = radius of core
    !        RMANT = radius of mantle
    !        WAVEL = wavelength of light in ambient medium
    !
    ! Routine BHCOAT is taken from Bohren & Huffman (1983)
    ! Obtained from C.L.Joseph
    !
    ! History:
    ! 92/11/24 (BTD) Explicit declaration of all variables
    !***********************************************************************
    !
    !         -----------------------------------------------------------
    !              del is the inner sphere convergence criterion
    !         -----------------------------------------------------------
    x1 = rfrel1*x
    x2 = rfrel2*x
    y2 = rfrel2*y
    Nstop = y + 4.*y**0.3333 + 2.0
    refrel = rfrel2/rfrel1
    !         -----------------------------------------------------------
    !              series terminated after nstop terms
    !         -----------------------------------------------------------
    d0x1 = COS(x1)/SIN(x1)
    d0x2 = COS(x2)/SIN(x2)
    D0Y  = COS(Y)/SIN(Y)
    d0y2 = COS(y2)/SIN(y2)
    chi0y = -sin(y)
    chi1y = cos(y)
    chi0y2 = -SIN(y2)
    chi1y2 = COS(y2)
    chi0x2 = -SIN(x2)
    chi1x2 = COS(x2)
    qsca = 0.0
    qext = 0.0
    xback = (0.0,0.0)
    iflag = .FALSE.
    n=1
    DO while(n<=nstop) 
       chiy = (2.0*n-1.)*chi1y/y - chi0y
       D1Y  = 1.0/(n/y-d0y) - n/y
       d1y2 = 1.0/(n/y2-d0y2) - n/y2
!       print*,'step1',n,nstop
       IF(.not.iflag) THEN
          d1x1 = 1.0/(n/x1-d0x1) - n/x1
          d1x2 = 1.0/(n/x2-d0x2) - n/x2
          chix2 = (2.0*n - 1.0)*chi1x2/x2 - chi0x2
          chiy2 = (2.0*n - 1.0)*chi1y2/y2 - chi0y2
          chipx2 = chi1x2 - n*chix2/x2
          chipy2 = chi1y2 - n*chiy2/y2
          CHIPY  = CHI1Y  - N*CHIY/Y
          brack=(refrel*d1x1-d1x2)/(refrel*d1x1*chix2-chipx2)/(chix2*d1x2-chipx2) &
               *(chiy2*d1y2-chipy2)
          crack=(refrel*d1x2-d1x1)/(refrel*chipx2-d1x1*chix2)/(chix2*d1x2-chipx2) &
               *(chiy2*d1y2-chipy2)
       END IF
       IF(ABS(brack*chipy2)<=del*ABS(d1y2).AND.ABS(brack*chiy2)<=del &
            .AND.ABS(crack*chipy2)<=del*ABS(d1y2).AND.ABS(crack*chiy2)<=del&
            .AND..NOT.iflag) THEN
          brack = (0.,0.)
          crack = (0.,0.)
          iflag = .true.
       ELSE
          dnbar = (d1y2 - brack*chipy2)/(1.0-brack*chiy2)
          gnbar = (d1y2 - crack*chipy2)/(1.0-crack*chiy2)
          an = 1./(1.-II*(chiy*d1y-chipy)*(CHIY*dnbar-rfrel2*chipy)/(dnbar-rfrel2*d1y))
          bn = 1./(1.-II*(chiy*d1y-chipy)*(rfrel2*CHIY*gnbar-chipy)/(rfrel2*gnbar-d1y))
          qsca = qsca + (2.0*n+1.0)*(ABS(an)**2+ABS(bn)**2)
          xback = xback + (2.0*n+1.0)*(-1.)**n*(an-bn)
          qext = qext + (2.0*n + 1.0)*(an+bn)
          chi0y = chi1y
          chi1y = chiy
          chi0x2 = chi1x2
          chi1x2 = chix2
          chi0y2 = chi1y2
          chi1y2 = chiy2
          d0x1 = d1x1
          d0x2 = d1x2
          D0Y  = D1Y
          d0y2 = d1y2
          n=n+1
       END IF
    END DO
    QSCA = (2.0/(y*y))*qsca
    QEXT = (2.0/(y*y))*qext
    qback = (1.0/(y*y))*ABS(XBACK)**2
    !qback=4.*Y**4*ABS((RFREL1**2-1)/(RFREL1**2+2))**2
    RETURN
  END SUBROUTINE BHCOAT
!
!-------------------------------------------------------------------------------
!
!*       6.   FUNCTION MG
!             -------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   ##########################################
    FUNCTION MG(QEPSINC,QEPSMAT,PF) RESULT(PQEPS)
!   ##########################################
    ! Maxwell Garnett (1904) equation for dielectric function of effective medium (sphere inclusions in a matrix)
    
    IMPLICIT NONE
    COMPLEX, INTENT(IN)  :: QEPSINC ! dielectric function of inclusions
    COMPLEX, INTENT(IN)  :: QEPSMAT ! dielectric function of matrix
    REAL,    INTENT(IN)  :: PF      ! volume fraction of the inclusions
    COMPLEX              :: PQEPS   ! dielectric function of effective medium
    COMPLEX              :: QEPSF

    QEPSF=(QEPSINC-QEPSMAT)/(QEPSINC+2.*QEPSMAT)
    PQEPS=(1.+2.*PF*QEPSF)/(1.-PF*QEPSF)*QEPSMAT
!
  END FUNCTION MG
!
!-------------------------------------------------------------------------------
END MODULE MODE_FSCATTER
