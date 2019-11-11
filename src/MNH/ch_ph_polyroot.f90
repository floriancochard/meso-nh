!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##########################
      MODULE MODI_CH_PH_POLYROOT
!     ##########################
INTERFACE
      SUBROUTINE CH_PH_POLYROOT(PPCOEF, KORDER, OPOLISH, PPALL_ROOTS)
!
COMPLEX (KIND(0.0D0)), DIMENSION(:), INTENT(IN)    :: PPCOEF  ! Polynomial coef.
INTEGER,                             INTENT(IN)    :: KORDER  ! Polynomial deg.
LOGICAL,                             INTENT(IN)    :: OPOLISH ! A better 
                                                              ! accuracy if TRUE
COMPLEX (KIND(0.0D0)), DIMENSION(:), INTENT(INOUT) :: PPALL_ROOTS ! Complex 
                                                                  ! roots
!
END SUBROUTINE
END INTERFACE
END MODULE MODI_CH_PH_POLYROOT
!
!     ###############################################################
      SUBROUTINE CH_PH_POLYROOT(PPCOEF, KORDER, OPOLISH, PPALL_ROOTS)
!     ###############################################################
!!    REFERENCE
!!    ---------
!!      Numerical Recipes in Fortran: Press, Teukolsky, Vetterling and Flannery
!!      pp 365-368.
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty *LA*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/07
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!
COMPLEX (KIND(0.0D0)), DIMENSION(:), INTENT(IN)    :: PPCOEF  ! Polynomial coef.
INTEGER,                             INTENT(IN)    :: KORDER  ! Polynomial deg.
LOGICAL,                             INTENT(IN)    :: OPOLISH ! A better 
                                                              ! accuracy if TRUE
COMPLEX (KIND(0.0D0)), DIMENSION(:), INTENT(INOUT) :: PPALL_ROOTS ! Complex 
                                                                  ! roots
!
!*       0.2   Declaration of local variables
!
INTEGER :: JJ, JL                        ! loop counters
INTEGER :: IITER                         ! number of iteration in 
                                         ! Laguerre routine
COMPLEX (KIND(0.0D0)), DIMENSION(SIZE(PPCOEF)) :: ZZDEFLATED_COEF
COMPLEX (KIND(0.0D0))                          :: ZROOT, ZB, ZC   ! ancillaries
REAL    :: ZEPS=1.0E-6
!
!-------------------------------------------------------------------------------
!
ZZDEFLATED_COEF(:) = PPCOEF(:)
!
!  First estimate of the roots
!
DO JJ=KORDER,1,-1
  ZROOT = CMPLX(0.0,0.0)
  CALL LAGUERRE(ZZDEFLATED_COEF, JJ, ZROOT, IITER)
  IF( ABS(AIMAG(ZROOT))<=2.0*ZEPS**(2*ABS(REAL(ZROOT))) ) THEN
    ZROOT = CMPLX(REAL(ZROOT),0.0)
  END IF
  PPALL_ROOTS(JJ) = ZROOT
  ZB = ZZDEFLATED_COEF(JJ+1)
  DO JL=JJ,1,-1
    ZC= ZZDEFLATED_COEF(JL)
    ZZDEFLATED_COEF(JL) = ZB
    ZB = ZROOT*ZB+ZC
  END DO
END DO
!
!  Further increase of the accuracy of the roots
!
IF (OPOLISH) THEN
  DO JJ=1,KORDER
    CALL LAGUERRE(PPCOEF, KORDER, PPALL_ROOTS(JJ), IITER)
  END DO
END IF
!
DO JJ=2,KORDER
  ZROOT = PPALL_ROOTS(JJ)
  DO JL=JJ-1,1,-1
    IF(REAL(PPALL_ROOTS(JL)) <= REAL(ZROOT)) GO TO 10
    PPALL_ROOTS(JL+1) = PPALL_ROOTS(JL)
  END DO
  JL = 0
10  PPALL_ROOTS(JL+1) = ZROOT
END DO
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE LAGUERRE(PA, IM, PX, IITS)
!
  COMPLEX (KIND(0.0D0)), DIMENSION(:), INTENT(IN)    :: PA
  INTEGER                            , INTENT(IN)    :: IM
  COMPLEX (KIND(0.0D0))              , INTENT(INOUT) :: PX
  INTEGER                            , INTENT(OUT)   :: IITS  
!
  COMPLEX (KIND(0.0D0)) :: ZZDX,ZZX1,ZZB,ZZD,ZZF,ZZG,ZZH,ZZSQ,ZZGP,ZZGM,ZZG2
  REAL, DIMENSION(8)    :: ZFRAC=(/0.50,0.25,0.75,0.13,0.38,0.62,0.88,1.00/)
  REAL(KIND(0.0D0))     :: ZABX,ZABP,ZABM,ZERR
  INTEGER               :: JJ, JITER ! iterative loops
  INTEGER               :: IMR=8 ! dimension of ZFRAC
  INTEGER               :: IMT=10
  INTEGER               :: IMAXIT
  REAL                  :: ZEPSS=2.0E-7
!
!
!
  IMAXIT= IMT*IMR
  LOOP:  DO JITER=1,IMAXIT
    IITS = JITER
    ZZB  = PA(IM+1)
    ZERR = ABS(ZZB)
    ZZD  = CMPLX(0.0,0.0)
    ZZF  = CMPLX(0.0,0.0)
    ZABX = ABS(PX)
    DO JJ=IM,1,-1
      ZZF = PX*ZZF+ZZD
      ZZD = PX*ZZD+ZZB
      ZZB = PX*ZZB+PA(JJ)
      ZERR = ABS(ZZB)+ZABX*ZERR
    END DO
!
    ZERR = ZEPSS*ZERR
    IF(ABS(ZZB) <= ZERR) THEN
      EXIT LOOP
    ELSE
      ZZG  = ZZD/ZZB
      ZZG2 = ZZG*ZZG
      ZZH  = ZZG2 - 2.0*(ZZF/ZZB)
      ZZSQ = SQRT( FLOAT(IM-1)*(FLOAT(IM)*ZZH-ZZG2) ) 
      ZZGP = ZZG + ZZSQ
      ZZGM = ZZG - ZZSQ
!
      ZABP = ABS(ZZGP)
      ZABM = ABS(ZZGM)
      IF(ZABP < ZABM) THEN
        ZZGP = ZZGM
      END IF
      IF(MAX(ZABP,ZABM) > 0.0) THEN
        ZZDX = FLOAT(IM)/ZZGP
        ELSE
        ZZDX = EXP(CMPLX(LOG(1.0+ZABX),FLOAT(JITER)))
      END IF 
    END IF
    ZZX1 = PX-ZZDX
    IF(PX == ZZX1) EXIT LOOP ! Convergence is reached
    IF(MOD(JITER,IMT) /= 0) THEN
      PX = ZZX1
      ELSE
      PX = PX - ZZDX*ZFRAC(JITER/IMT)
    END IF
  END DO LOOP
  END SUBROUTINE LAGUERRE
END SUBROUTINE CH_PH_POLYROOT
