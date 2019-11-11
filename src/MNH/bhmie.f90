!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!      #################
       MODULE MODI_BHMIE
!      #################
!
INTERFACE
      SUBROUTINE BHMIE(PSIZE_PARAM,PPREFR,KNANG,PPS1,PPS2,PQEXT,PQBAK)
!
REAL,                  INTENT(IN)  :: PSIZE_PARAM
COMPLEX,               INTENT(IN)  :: PPREFR
INTEGER,               INTENT(IN)  :: KNANG
COMPLEX, DIMENSION(:), INTENT(OUT) :: PPS1
COMPLEX, DIMENSION(:), INTENT(OUT) :: PPS2
REAL,                  INTENT(OUT) :: PQEXT,PQBAK
!
END SUBROUTINE BHMIE
END INTERFACE
END MODULE MODI_BHMIE 
!
!     ################################################################
      SUBROUTINE BHMIE(PSIZE_PARAM,PPREFR,KNANG,PPS1,PPS2,PQEXT,PQBAK)
!     ################################################################
!
!!***********************************************************************
!! Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
!!    to calculate scattering and absorption by a homogenous isotropic
!!    sphere.
!! Given:
!!    X = 2*pi*a/lambda
!!    REFREL = (complex refr. index of sphere)/(real index of medium)
!!    NANG = number of angles between 0 and 90 degrees
!!           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
!!           if called with NANG<2, will set NANG=2 and will compute
!!           scattering for theta=0,90,180.
!! Returns:
!!    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
!!                                scatt. E perp. to scatt. plane)
!!    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
!!                                scatt. E parr. to scatt. plane)
!!    PQEXT = C_ext/pi*a**2 = efficiency factor for extinction
!!    PQBAK = (dC_sca/domega)/pi*a**2
!!          = backscattering efficiency [NB: this is (1/4*pi) smaller
!!            than the "radar backscattering efficiency"; see Bohren &
!!            Huffman 1983 pp. 120-123]
!!
!! Original program taken from Bohren and Huffman (1983), Appendix A
!! Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
!! in order to compute <cos(theta)>
!! 91/05/07 (BTD): Modified to allow NANG=1
!! 91/08/15 (BTD): Corrected error (failure to initialize P)
!! 91/08/15 (BTD): Modified to enhance vectorizability.
!! 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
!! 91/08/15 (BTD): Changed definition of PQBAK.
!! 92/01/08 (BTD): Converted to full double precision and double complex
!!                 eliminated 2 unneed lines of code
!!                 eliminated redundant variables (e.g. APSI,APSI0)
!!                 renamed RN -> EN = double precision N
!!                 Note that DOUBLE COMPLEX and DCMPLX are not part
!!                 of f77 standard, so this version may not be fully
!!                 portable.  In event that portable version is
!!                 needed, use src/bhmie_f77.f
!! 93/06/01 (BTD): Changed AMAX1 to generic function MAX
!! 22/01/2019 (P.Wautelet): correct kind of complex datatype
!!***********************************************************************
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                  INTENT(IN)  :: PSIZE_PARAM
COMPLEX,               INTENT(IN)  :: PPREFR
INTEGER,               INTENT(IN)  :: KNANG
COMPLEX, DIMENSION(:), INTENT(OUT) :: PPS1
COMPLEX, DIMENSION(:), INTENT(OUT) :: PPS2
REAL,                  INTENT(OUT) :: PQEXT,PQBAK
!
!*       0.1   Declarations of local variables :
!
REAL (KIND(0.0D0)),    DIMENSION(:), ALLOCATABLE :: ZPI,ZPI0,ZPI1
REAL (KIND(0.0D0)),    DIMENSION(:), ALLOCATABLE :: ZAMU,ZTAU
COMPLEX (KIND(0.0D0)), DIMENSION(:), ALLOCATABLE :: ZZD
!
INTEGER :: J,JJ,JJJ
INTEGER :: IEN
INTEGER :: ISTOP,INMX
REAL (KIND(0.0D0))    :: ZCHI,ZCHI0,ZCHI1
REAL (KIND(0.0D0))    :: ZPSI,ZPSI0,ZPSI1
REAL (KIND(0.0D0))    :: ZDELTA_ANGLE,ZEN,ZFN,ZONE,ZSIZE_PARAM_STOP,ZYMOD
COMPLEX (KIND(0.0D0)) :: ZZAN,ZZAN1,ZZBN,ZZBN1,ZZXI,ZZXI1,ZZEN,ZZY
!
!***********************************************************************
!
ZZY = PSIZE_PARAM*PPREFR
ZYMOD = ABS(ZZY)
!
!*** Series expansion terminated after ISTOP terms
!    Logarithmic derivatives calculated from INMX on down
!
ZSIZE_PARAM_STOP = PSIZE_PARAM+4.*PSIZE_PARAM**0.3333+2.
INMX  = INT(MAX(ZSIZE_PARAM_STOP,ZYMOD))+15
ISTOP = INT(ZSIZE_PARAM_STOP)
!
! BTD experiment 91/1/15: add one more term to series and compare results
!      INMX=AMAX1(ZSIZE_PARAM_STOP,ZYMOD)+16
! test: compute 7001 wavelengths between .0001 and 1000 micron
! for a=1.0micron SiC grain.  When INMX increased by 1, only a single
! computed number changed (out of 4*7001) and it only changed by 1/8387
! conclusion: we are indeed retaining enough terms in series!
!
ALLOCATE(ZTAU(KNANG))
ALLOCATE(ZAMU(KNANG))
ZDELTA_ANGLE = 0.5*XPI/FLOAT(KNANG-1)
DO J = 1,KNANG
  ZAMU(J) = COS( ZDELTA_ANGLE*FLOAT(J-1) )
ENDDO
!
ALLOCATE(ZPI(KNANG))
ZPI(:) = 0.
ALLOCATE(ZPI0(KNANG))
ALLOCATE(ZPI1(KNANG))
ZPI0(:) = 0.
ZPI1(:) = 1.
!
PPS1(:) = (0.,0.)
PPS2(:) = (0.,0.)
!
!*** Logarithmic derivative D(J) calculated by downward recurrence
!    beginning with initial value (0.,0.) at J=INMX
!
ALLOCATE(ZZD(INMX))
ZZD(INMX) = (0.,0.)
!
DO J = 1,INMX-1
  IEN = INMX-J+1
  ZZEN = FLOAT(IEN)/ZZY
  ZZD(INMX-J) = ZZEN-(1.0/(ZZD(IEN)+ZZEN))
ENDDO
!
!*** Riccati-Bessel functions with real argument X
!    calculated by upward recurrence
!
ZPSI0 = COS(PSIZE_PARAM)
ZPSI1 = SIN(PSIZE_PARAM)
ZCHI0 =-SIN(PSIZE_PARAM)
ZCHI1 = COS(PSIZE_PARAM)
ZZXI1 = CMPLX(ZPSI1,-ZCHI1)
ZONE = -1.
!
ZZAN1 = CMPLX(0.0,0.0)
ZZBN1 = CMPLX(0.0,0.0)
DO J = 1,ISTOP
  ZEN = FLOAT(J)
  ZFN = (2.0*ZEN+1.0)/(ZEN*(ZEN+1.0))
!
! for given N, ZPSI  = psi_n        ZCHI  = chi_n
!              ZPSI1 = psi_{n-1}    ZCHI1 = chi_{n-1}
!              ZPSI0 = psi_{n-2}    ZCHI0 = chi_{n-2}
! Calculate psi_n and chi_n
!
  ZPSI = (2.0*ZEN-1.0)*ZPSI1/PSIZE_PARAM-ZPSI0
  ZCHI = (2.0*ZEN-1.0)*ZCHI1/PSIZE_PARAM-ZCHI0
  ZZXI = CMPLX(ZPSI,-ZCHI)
!
!*** Compute AN and BN:
!
  ZZAN = ((ZZD(J)/PPREFR+ZEN/PSIZE_PARAM)*ZPSI-ZPSI1)/ &
         ((ZZD(J)/PPREFR+ZEN/PSIZE_PARAM)*ZZXI-ZZXI1)
  ZZBN = ((PPREFR*ZZD(J)+ZEN/PSIZE_PARAM)*ZPSI-ZPSI1)/ &
         ((PPREFR*ZZD(J)+ZEN/PSIZE_PARAM)*ZZXI-ZZXI1)
!
!*** Store previous values of AN and BN for use
!    in computation of g=<cos(theta)>
!
  ZZAN1 = ZZAN
  ZZBN1 = ZZBN
!
!*** Now calculate scattering intensity pattern
!    First do angles from 0 to 90
!
  DO JJ = 1,KNANG
    ZPI(JJ) = ZPI1(JJ)
    ZTAU(JJ) = ZEN*ZAMU(JJ)*ZPI(JJ) - (ZEN+1.)*ZPI0(JJ)
    PPS1(JJ) = PPS1(JJ) + CMPLX(ZFN*(ZZAN*ZPI(JJ)+ZZBN*ZTAU(JJ)),kind=kind(PPS1(1)))
    PPS2(JJ) = PPS2(JJ) + CMPLX(ZFN*(ZZAN*ZTAU(JJ)+ZZBN*ZPI(JJ)),kind=kind(PPS2(1)))
  ENDDO
!
!*** Now do angles greater than 90 using PI and TAU from
!    angles less than 90.
!    P=1 for N=1,3,...; P=-1 for N=2,4,...
!
  ZONE = -ZONE
  DO JJ = 1,KNANG-1
    JJJ = 2*KNANG-JJ
    PPS1(JJJ) = PPS1(JJJ) + CMPLX(ZFN*ZONE*(ZZAN*ZPI(JJ)-ZZBN*ZTAU(JJ)),kind=kind(PPS1(1)))
    PPS2(JJJ) = PPS2(JJJ) + CMPLX(ZFN*ZONE*(ZZBN*ZPI(JJ)-ZZAN*ZTAU(JJ)),kind=kind(PPS2(1)))
  ENDDO
  ZPSI0 = ZPSI1
  ZPSI1 = ZPSI
  ZCHI0 = ZCHI1
  ZCHI1 = ZCHI
  ZZXI1 = CMPLX(ZPSI1,-ZCHI1)
!
!*** Compute pi_n for next value of n
!    For each angle J, compute pi_n+1
!    from PI = pi_n , ZPI0 = pi_n-1
!
  ZPI1(:) = ((2.0*ZEN+1.0)*ZAMU(:)*ZPI(:)-(ZEN+1.0)*ZPI0(:))/ZEN
  ZPI0(:) = ZPI(:)
ENDDO
!
!*** Have summed sufficient terms.
!    Now compute PQEXT,PQBAK
!
      PQEXT=(4./(PSIZE_PARAM*PSIZE_PARAM))*REAL(PPS1(1))
      PQBAK=(2.*ABS(PPS1(2*KNANG-1))/PSIZE_PARAM)**2
!
!*** Deallocations
!
DEALLOCATE(ZTAU)
DEALLOCATE(ZAMU)
DEALLOCATE(ZPI)
DEALLOCATE(ZPI0)
DEALLOCATE(ZPI1)
DEALLOCATE(ZZD)
!
RETURN
END SUBROUTINE BHMIE
