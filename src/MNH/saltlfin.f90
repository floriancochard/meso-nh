!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/10/19 17:13:51
!-----------------------------------------------------------------
!     ########################
      MODULE MODI_SALTLFI_n
!     ########################
!
INTERFACE
!
!++cb++24/10/16
!SUBROUTINE SALTLFI_n(PSV, PRHODREF)
SUBROUTINE SALTLFI_n(PSV, PRHODREF, PZZ)
IMPLICIT NONE
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PZZ

END SUBROUTINE SALTLFI_n
!
END INTERFACE
!
END MODULE MODI_SALTLFI_n
!
!
!     ############################################################
!      SUBROUTINE SALTLFI_n(PSV, PRHODREF)
      SUBROUTINE SALTLFI_n(PSV, PRHODREF, PZZ)
!     ############################################################
!
!!    PURPOSE
!!    -------
!!    Realise l'équilibre des moments à partir du sigma et du diametre moyen
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    2014    P.Tulet  modif calcul ZM
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes

!!    EXTERNAL
!!    --------
!!    None
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SALT
USE MODD_NSV
!++cb++24/10/16
!USE MODD_GRID_n, ONLY: XZZ
!--cb--
USE MODD_CSTS_SALT
USE MODD_CST, ONLY :    &
       XPI              & !Definition of pi
      ,XBOLTZ           & ! Boltzman constant 
      ,XAVOGADRO        & ![molec/mol] avogadros number
      ,XG               & ! Gravity constant
      ,XP00             & ! Reference pressure
      ,XMD              & ![kg/mol] molar weight of air
      ,XRD              & ! Gaz constant for dry air
      ,XCPD               !  Cpd (dry air)
! 
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:,:,:),    INTENT(INOUT) :: PSV
REAL,   DIMENSION(:,:,:),      INTENT(IN) :: PRHODREF
REAL,   DIMENSION(:,:,:),      INTENT(IN) :: PZZ
!
!
!*      0.2    declarations local variables
!
REAL   :: ZDEN2MOL, ZRHOI, ZMI, ZFAC, ZRGMIN
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZCTOTA
REAL,DIMENSION(:,:,:,:), ALLOCATABLE  :: ZM
REAL,DIMENSION(:,:,:),   ALLOCATABLE  :: ZSIGMA
INTEGER,DIMENSION(:),    ALLOCATABLE  :: IM0, IM3, IM6
REAL,DIMENSION(:),       ALLOCATABLE  :: ZMMIN
REAL,DIMENSION(:),       ALLOCATABLE  :: ZINIRADIUS, ZINISIGMA
REAL,DIMENSION(:,:),     ALLOCATABLE  :: ZSEA
INTEGER :: IKU
!+Marine
INTEGER :: IMOMENTS
!-Marine
INTEGER :: JI, JJ, JN, JK  ! loop counter
INTEGER :: IMODEIDX  ! index mode
REAL, PARAMETER  :: ZN_SALT=1E4 ! multiplcative factor for X0MIN
REAL, PARAMETER  :: ZCLM=800. ! Marine Salt layer (m)
REAL    :: ZN_SALTN
!
!-------------------------------------------------------------------------------
!
!*       1.     TRANSFER FROM GAS TO AEROSOL MODULE
!               -----------------------------------
!
!        1.1    initialisation 
!
IKU=SIZE(PSV,3)
!+ Marine
!
ALLOCATE (IM0(NMODE_SLT))
ALLOCATE (IM3(NMODE_SLT))
ALLOCATE (IM6(NMODE_SLT))
ALLOCATE (ZCTOTA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT))
ALLOCATE (ZM(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT*3))
ALLOCATE (ZSIGMA(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3)))
ALLOCATE (ZINIRADIUS(NMODE_SLT))
ALLOCATE (ZINISIGMA(NMODE_SLT))
ALLOCATE (ZMMIN(NMODE_SLT*3))
ALLOCATE (ZSEA(SIZE(PSV,1), SIZE(PSV,2)))

ZSEA(:,:) = 0.
!++cb++20/10/16
!WHERE ((XZZ(:,:,1) .LT. 0.1).AND.(XZZ(:,:,1) .GE. 0.)) 
!  ZSEA(:,:) = 1.
!END WHERE
!++cb++24/10/16
!WHERE (XZZ(:,:,1) .LE. 0.01) 
WHERE (PZZ(:,:,1) .LE. 0.01)
!--cb--
  ZSEA(:,:) = 1.
END WHERE
!--cb--
!
!
!+ Marine
DO JN = 1, NMODE_SLT
  IM0(JN) = 1+(JN-1)*3
  IM3(JN) = 2+(JN-1)*3
  IM6(JN) = 3+(JN-1)*3
  !
  !Get the sea salt mode we are talking about, MODE 2 is treated first, then mode 3, then 1
  !This index is only needed to get the right radius out of the XINIRADIUS_SLT array and the
  !right XINISIG_SLT out of the XINISIG_SLT-array
  IMODEIDX = JPSALTORDER(JN)
  !
  !Convert initial mass median radius to number median radius
  IF (CRGUNITS=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
  ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
  END IF
  ZINISIGMA(JN)  = XINISIG_SLT(IMODEIDX)
  !
  ZMMIN(IM0(JN)) = XN0MIN_SLT(IMODEIDX)
! ZRGMIN   = XCOEFRADMIN * ZINIRADIUS(JN)
  ZRGMIN   = ZINIRADIUS(JN)
  ZMMIN(IM3(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**3)*EXP(4.5 * LOG(ZINISIGMA(JN))**2) 
  ZMMIN(IM6(JN)) = XN0MIN_SLT(IMODEIDX) * (ZRGMIN**6)*EXP(18. * LOG(ZINISIGMA(JN))**2)
ENDDO
!
!
!XDENSITY_SALT est fixé dans modd_csts_salt.f90
ZRHOI = XDENSITY_SALT 
ZMI   = XMOLARWEIGHT_SALT 
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
ZFAC=(4./3.)*XPI*ZRHOI*1.e-9

!
DO JN=1,NMODE_SLT

!*       1.1    calculate moment 0 from sea salt number by m3
!
! initial vertical profil of sea salt  and convert in  #/m3
!+Marine : (reprendre XN0MIN_SLT de modd_salt.f90).
! Pas plus simple de fixer une dimension à ZN_SALT qui dépend de JN pour ne pas
! avoir à rappeler le schéma d'émission?
  IF(NMODE_SLT == 5)THEN
    IF (JN == 1) ZN_SALTN = XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT 
    IF (JN == 2) ZN_SALTN = XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT
    IF (JN == 3) ZN_SALTN = XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT
    IF (JN == 4) ZN_SALTN = XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT
    IF (JN == 5) ZN_SALTN = XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT
  ELSE 
    IF (JN == 1) ZN_SALTN =  XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT 
    IF (JN == 2) ZN_SALTN =  XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT
    IF (JN == 3) ZN_SALTN =  XN0MIN_SLT(JPSALTORDER(JN)) *  ZN_SALT
  END IF
!-Marine
  DO JK=1, SIZE(PSV,3) 
    DO JJ=1, SIZE(PSV,2) 
      DO JI=1, SIZE(PSV,1) 
!++cb++24/10/16
!        IF (XZZ(JI,JJ,JK) .LT. 600.) THEN
        IF (PZZ(JI,JJ,JK) .LT. 600.) THEN
          ZM(JI,JJ,JK,IM0(JN)) = ZN_SALTN 
!        ELSE IF ((XZZ(JI,JJ,JK) .GE. 600.).AND.(XZZ(JI,JJ,JK) .LT. 1000.)) THEN 
        ELSE IF ((PZZ(JI,JJ,JK) .GE. 600.).AND.(PZZ(JI,JJ,JK) .LT. 1000.)) THEN
!          ZM(JI,JJ,JK,IM0(JN)) = ZN_SALTN - ZN_SALTN*(1.-1E-3)*(XZZ(JI,JJ,JK)-600.) / 400.
          ZM(JI,JJ,JK,IM0(JN)) = ZN_SALTN - &
                                 ZN_SALTN * (1.-1E-3) * (PZZ(JI,JJ,JK)-600.) / 400.
!        ELSE IF (XZZ(JI,JJ,JK) .GE. 1000.) THEN
        ELSE IF (PZZ(JI,JJ,JK) .GE. 1000.) THEN
          ZM(JI,JJ,JK,IM0(JN)) = ZN_SALTN * 1E-3
!--cb--
        END IF
      END DO
    END DO
    ! Over continent value of the free troposphere
    WHERE (ZSEA(:,:) == 0.)
      ZM(:,:,JK,IM0(JN)) = ZN_SALTN *1E-3
    END WHERE
    WHERE ((ZSEA(:,:) .GT. 0.).AND.(ZSEA(:,:) .LT. 1.))
      ZM(:,:,JK,IM0(JN)) = ZM(:,:,JK,IM0(JN))-(ZM(:,:,JK,IM0(JN)) -ZN_SALTN *1E-3) * &
                           (1. - ZSEA(:,:))
    END WHERE
  END DO

  ZM(:,:,:,IM0(JN)) = MAX(ZMMIN(IM0(JN)), ZM(:,:,:,IM0(JN)))
!
!*       1.2    calculate moment 3 from m0,  RG and SIG 
!
  ZM(:,:,:,IM3(JN)) = ZM(:,:,:,IM0(JN)) * &
              (ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(ZINISIGMA(JN))**2) 
  ZM(:,:,:,IM3(JN)) = MAX(ZMMIN(IM3(JN)), ZM(:,:,:,IM3(JN)))
!
!*       1.3    calculate moment 6 from m0,  RG and SIG 
!
  ZM(:,:,:,IM6(JN))= ZM(:,:,:,IM0(JN)) * ((ZINIRADIUS(JN)**6)*&
                        EXP(18. * (LOG(ZINISIGMA(JN)))**2))
  ZM(:,:,:,IM6(JN)) = MAX(ZMMIN(IM6(JN)), ZM(:,:,:,IM6(JN)))
!
!*       1.4    output concentration
!+ Marine
!  PSV(:,:,:,1+(JN-1)*3) = ZM(:,:,:,IM0(JN)) * XMD / (XAVOGADRO*PRHODREF(:,:,:))
!  PSV(:,:,:,2+(JN-1)*3) = ZM(:,:,:,IM3(JN)) * XMD*XPI * 4./3.  / &
!                           (ZMI*PRHODREF(:,:,:)*(1.d0/ZRHOI)*XM3TOUM3_SALT)
!
!  PSV(:,:,:,3+(JN-1)*3) = ZM(:,:,:,IM6(JN)) *  XMD / (XAVOGADRO*PRHODREF(:,:,:)*1.d-6)
!
!++cb++20/10/16
  IMOMENTS = INT(NSV_SLTEND - NSV_SLTBEG + 1) / NMODE_SLT
!--cb--

  IF (IMOMENTS == 3) THEN
    PSV(:,:,:,1+(JN-1)*3) = ZM(:,:,:,IM0(JN)) * XMD / (XAVOGADRO*PRHODREF(:,:,:))
    PSV(:,:,:,2+(JN-1)*3) = ZM(:,:,:,IM3(JN)) * XMD*XPI * 4./3.  / &
                            (ZMI*PRHODREF(:,:,:)*(1.d0/ZRHOI)*XM3TOUM3_SALT)

    PSV(:,:,:,3+(JN-1)*3) = ZM(:,:,:,IM6(JN)) * XMD / (XAVOGADRO*PRHODREF(:,:,:)*1.d-6)
  ELSE IF (IMOMENTS == 2) THEN
    PSV(:,:,:,1+(JN-1)*2) = ZM(:,:,:,IM0(JN)) * XMD / (XAVOGADRO*PRHODREF(:,:,:))
    PSV(:,:,:,2+(JN-1)*2) = ZM(:,:,:,IM3(JN)) * XMD*XPI * 4./3.  / &
                            (ZMI*PRHODREF(:,:,:)*(1.d0/ZRHOI)*XM3TOUM3_SALT)
  ELSE 
    PSV(:,:,:,JN) = ZM(:,:,:,IM3(JN)) * XMD*XPI * 4./3.  / &
                            (ZMI*PRHODREF(:,:,:)*(1.d0/ZRHOI)*XM3TOUM3_SALT)
  END IF
!
END DO
!
DEALLOCATE(ZSEA)
DEALLOCATE(ZMMIN)
DEALLOCATE(ZINISIGMA)
DEALLOCATE(ZINIRADIUS)
DEALLOCATE(ZSIGMA)
DEALLOCATE(ZM)
DEALLOCATE(ZCTOTA)
DEALLOCATE(IM6)
DEALLOCATE(IM3)
DEALLOCATE(IM0)
!
!
END SUBROUTINE SALTLFI_n
