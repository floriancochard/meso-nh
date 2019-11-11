!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/10/18 12:20:58
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_CH_AER_WETDEP_n
!!   ##############################
!!
INTERFACE
!
SUBROUTINE CH_AER_WETDEP_n(PDTMONITOR,PSVT,PCWETDEP,&
                 PRHODREF,PSEDA)

IMPLICIT NONE

REAL,                       INTENT(IN)    :: PDTMONITOR
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PSVT, PCWETDEP
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSEDA


END SUBROUTINE CH_AER_WETDEP_n
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_WETDEP_n
!!
!!   #######################################
     SUBROUTINE CH_AER_WETDEP_n(PDTMONITOR,PSVT,PCWETDEP,&
                       PRHODREF, PSEDA)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!    Pierre TULET (LACy)
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!    M.Leriche 2015 correction bug
!!
! Entry variables:
!
! PM(IN)       -Array of moments
!
!************************************************************
!!
!!   IMPLICIT ARGUMENTS
!
USE MODD_CH_AEROSOL
USE MODD_NSV
USE MODD_CH_AERO_n, ONLY : GSEDFIX
USE MODE_AERO_PSD
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL,                       INTENT(IN)    :: PDTMONITOR
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PSVT, PCWETDEP
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSEDA
!
!*       0.2   Declarations of local variables
!
INTEGER ::  JN
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE*3) :: ZPM, ZPMOLD
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE)   :: ZRG, ZSIG, ZRGN, ZSIGN
REAL    :: ZINIRADIUSI, ZINIRADIUSJ, ZRGMIN
INTEGER             :: ILU  ! indice K End       in z direction
REAL,    DIMENSION(JPIN) :: ZPMIN

!-------------------------------------------------------------------------------
!
ILU = SIZE(PSVT,3) -1

IF (CRGUNIT=="MASS") THEN
    ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
    ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
    ZINIRADIUSI = XINIRADIUSI
    ZINIRADIUSJ = XINIRADIUSJ
END IF
!
!Get minimum values possible
ZPMIN(1) = XN0IMIN
ZRGMIN =  ZINIRADIUSI
ZPMIN(2) = ZPMIN(1) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGIMIN)**2)
ZPMIN(3) = ZPMIN(1) * (ZRGMIN**6)*EXP(18. * LOG(XSIGIMIN)**2)
ZPMIN(4) = XN0JMIN
ZRGMIN =  ZINIRADIUSJ
ZPMIN(5) = ZPMIN(4) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGJMIN)**2)
ZPMIN(6) = ZPMIN(4) * (ZRGMIN**6)*EXP(18. * LOG(XSIGJMIN)**2)
!
CALL  PPP2AERO(PSVT, PRHODREF, PSIG3D=ZSIG, PRG3D=ZRG, PM3D=ZPMOLD)
CALL  PPP2AERO(PCWETDEP, PRHODREF, PSIG3D=ZSIGN, PRG3D=ZRGN, PM3D=ZPM)
!
DO JN=1,JPMODE
WHERE (ZPM(:,:,:,NM3(JN)) .LT. ZPMIN(NM3(JN)))
   ZPM(:,:,:,NM3(JN)) = ZPMOLD(:,:,:,NM3(JN))
END WHERE
! Calcul pour maintenir le rayon fixe pour la sedimentation :
 IF (LRGFIX) THEN
    ZPM(:,:,:,NM0(JN)) = ZPM(:,:,:,NM3(JN)) / &
                        (ZRG(:,:,:,JN)**3)*EXP(4.5 * LOG(ZSIG(:,:,:,JN))**2)
 ENDIF

! Calcul pour maintenir la dispersion fixe pour la sedimentation :
! sinon Rg augmente lors de la reconstruction d'une loi log-normale
! a partir des nouveau Mk (sigma diminue plus vite que Rg)

! calcul de M6 en conservant sigma
ZPM(:,:,:, NM6(JN)) = ZPM(:,:,:,NM0(JN)) &
       * ( (ZPM(:,:,:,NM3(JN))/ZPM(:,:,:,NM0(JN)))**(1./3.) &
       * exp(-(3./2.)*LOG(ZSIG(:,:,:,JN))**2))**6 &
       * exp(18.*LOG(ZSIG(:,:,:,JN))**2)

  IF ((GSEDFIX).AND.&
      (((JN .EQ. 1).AND. (LVARSIGI)).OR.&
       ((JN .EQ. 2).AND. (LVARSIGJ)))) THEN
! calcul de M6 en conservant Rg
    ZPM(:,:,:,NM6(JN)) = ZPM(:,:,:,NM3(JN)) ** 4 / &
                     (ZRG(:,:,:,JN)**6 * ZPM(:,:,:,NM0(JN))**3)

  END IF
END DO
!
IF (GSEDFIX) THEN
  GSEDFIX = .FALSE.
ELSE
  GSEDFIX = .TRUE.
END IF

DO JN=1,JPMODE
  WHERE ((ZPM(:,:,:,NM0(JN)) .LT. ZPMIN(NM0(JN))).OR.&
         (ZPM(:,:,:,NM3(JN)) .LT. ZPMIN(NM3(JN))).OR.&
         (ZPM(:,:,:,NM6(JN)) .LT. ZPMIN(NM6(JN))))
    ZPM(:,:,:,NM0(JN)) = ZPMOLD(:,:,:,NM0(JN))
    ZPM(:,:,:,NM3(JN)) = ZPMOLD(:,:,:,NM3(JN))
    ZPM(:,:,:,NM6(JN)) = ZPMOLD(:,:,:,NM6(JN))
  END WHERE
ENDDO

DO JN=1,JPIN
 PSEDA(:,:,:,JN) = PSEDA(:,:,:,JN) + & 
                  (ZPM(:,:,:,JN) - ZPMOLD(:,:,:,JN)) / PDTMONITOR
END DO
!
END SUBROUTINE CH_AER_WETDEP_n
