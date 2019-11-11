!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/ch_aer_solv.f90,v $ $Revision: 1.1.2.1.2.1.16.2.2.1.2.1 $ $Date: 2015/12/01 15:26:23 $
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!!   #######################
     MODULE MODI_CH_AER_SOLV
!!   #######################
!!
INTERFACE
!!
SUBROUTINE CH_AER_SOLV(PM, PSIG0, PRG0, PN0,PCTOTG, PCTOTA, PCCTOT, &
                         PDMINTRA,PDMINTER,PDMCOND, PSEDA,PDT,&
                         POM, PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PTIME,PSOLORG)
IMPLICIT NONE
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL, DIMENSION(:,:),   INTENT(INOUT) :: POM
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PDMINTRA
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PDMINTER
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PDMCOND
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL, INTENT(IN)                      :: PDT, PTIME
REAL, DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC


END SUBROUTINE CH_AER_SOLV
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_SOLV
!!
!!   ##############################################################################
     SUBROUTINE CH_AER_SOLV(PM, PSIG0, PRG0, PN0,PCTOTG, PCTOTA, PCCTOT, &
                              PDMINTRA,PDMINTER,PDMCOND,PSEDA, PDT, POM, &
                              PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PTIME,PSOLORG)
!!   ##############################################################################
!!
!!   PURPOSE
!!   -------
!!   Time variable solver of the modal aerosol equations
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    P. Tulet for nesting
!!    P. Tulet organic condensation
!!    P. Tulet thermodynamic equilibrium for each mode
!!    P. Tulet add third mode
!!    M. Leriche 2015 correction bug
!!    M. Leriche 08/16 suppress moments index declaration already in modd_aerosol
!!    M. Leriche 08/16 add an other particular case for the M0 resolution to
!!               avoid a division by zero (when ZK = 1)
!!
!!    EXTERNAL
!!    --------
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CH_AEROSOL
USE MODD_CST, ONLY : XMNH_TINY
USE MODI_CH_AER_MINERAL
USE MODI_CH_AER_ORGANIC
USE MODI_CH_AER_MPMPO
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL, DIMENSION(:,:),   INTENT(INOUT) :: POM
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PDMINTRA
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PDMINTER
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PDMCOND
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL, INTENT(IN)                      :: PDT, PTIME
REAL, DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
!
!*       0.2   declarations of local variables
!
INTEGER :: JI,JJ,JK, JN, IDT
REAL, DIMENSION(SIZE(PM,1)) :: ZSUM
REAL, DIMENSION(SIZE(PM,1),JPMODE) :: ZNEWM, ZMASK
REAL, DIMENSION(SIZE(PM,1)) :: ZSIGMA

REAL, DIMENSION(SIZE(PM,1)) :: ZA,ZB,ZC,ZD
REAL, DIMENSION(SIZE(PM,1)) :: ZCONST1,ZCONST2
REAL, DIMENSION(SIZE(PM,1)) :: Z0,ZK,ZKEXP

REAL, SAVE,  DIMENSION(JPMODE*3)  :: ZPMIN
REAL, SAVE,  DIMENSION(JPMODE)   :: ZRATIOBC, ZRATIOOC
REAL    :: ZINIRADIUSI, ZINIRADIUSJ
REAL, SAVE,  DIMENSION(JPMODE)  :: ZRGMIN, ZRGMAX
LOGICAL, SAVE               :: GPHYSLIM = .TRUE. ! flag
!
!-------------------------------------------------------------------------------
!
  !
  IF (CRGUNIT=="MASS") THEN
    ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
    ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
  ELSE
    ZINIRADIUSI = XINIRADIUSI 
    ZINIRADIUSJ = XINIRADIUSJ
  END IF
  !
  ZPMIN(1) = XN0IMIN
  ZPMIN(2) = ZPMIN(1) * (ZINIRADIUSI**3)*EXP(4.5 * LOG(XSIGIMIN)**2) 
  ZPMIN(3) = ZPMIN(1) * (ZINIRADIUSI**6)*EXP(18. * LOG(XSIGIMIN)**2)
  !
  ZPMIN(4) = XN0JMIN
  ZPMIN(5) = ZPMIN(4) * (ZINIRADIUSJ**3)*EXP(4.5 * LOG(XSIGJMIN)**2) 
  ZPMIN(6) = ZPMIN(4) * (ZINIRADIUSJ**6)*EXP(18. * LOG(XSIGJMIN)**2)
  !
!write(*,*)
!write(*,*) '******************************************'
!write(*,*) '         Debut Solveur Aerosol            '
!write(*,*) '******************************************'
!write(*,*) 
!write(*,*) 'Pas de temps:',PDT,'s'

!*****************************************************************
!*****************************************************************
! SOLVEUR DE lA PARTIE MICROPHYSIQUE
!*****************************************************************
!*****************************************************************
!
DO JI=1,JPMODE

!*************************************************************
! Resolution du moment d'ordre 0: pour cela il faut resoudre
! une equation differentielle du type dY/dt=-AY^2-BY+C
!*************************************************************


! Pour la resolution plusieurs cas particuliers seront traites
ZA(:)=0.
ZB(:)=0.
ZC(:)=0.
ZA(:)=PDMINTRA(:,NM0(JI))
ZB(:)=PDMINTER(:,NM0(JI))
ZC(:)=PDMCOND(:,NM0(JI))


DO JK=1,SIZE(PM,1)
 IF ((ZB(JK) == 0. .AND. ZC(JK)/PM(JK,NM0(JI)) <= 1.e-10).OR. &
     (ZC(JK) <= 1.e-10 .AND. ZB(JK)/ZA(JK) <= 1.e-3))  THEN
! type dY/dt=-AY^2
   Z0(JK)=PM(JK,NM0(JI)) 
   PM(JK,NM0(JI))=Z0(JK)/(1.+ZA(JK)*Z0(JK)*PDT)
 ELSE
   ZCONST1(JK)=ZB(JK)/(2.*ZA(JK))
   Z0(JK)=PM(JK,NM0(JI))+ZCONST1(JK)
   IF (((ZB(JK)**2+4.*ZA(JK)*ZC(JK))) < 0.) THEN
     ZD(JK)=SQRT(ABS(ZB(JK)**2+4.*ZA(JK)*ZC(JK)))
     PM(JK,NM0(JI))=-ZCONST1(JK)+ZD(JK)*TAN(ATAN(Z0(JK)/ZD(JK))-ZA(JK)*ZD(JK)*PDT)
   ELSE
     ZD(JK)=SQRT(ZB(JK)**2+4.*ZA(JK)*ZC(JK))
     ZCONST2(JK)=ZD(JK)/(2.*ABS(ZA(JK)))
     ZKEXP(JK)=EXP(-2.*ZA(JK)*ZCONST2(JK)*PDT)
     ZK(JK)=(Z0(JK)-ZCONST2(JK))/(Z0(JK)+ZCONST2(JK))*ZKEXP(JK)
     PM(JK,NM0(JI))=-ZCONST1(JK)+ZCONST2(JK)*(1.+ZK(JK))/(1.-ZK(JK))
   ENDIF
 ENDIF
ENDDO

  ! Sedimentation for particules number

PM(:,NM0(JI))= PM(:,NM0(JI)) + PSEDA(:,NM0(JI)) * PDT
PM(:,NM0(JI))= MAX(PM(:,NM0(JI)), XMNH_TINY )

 
!*************************************************************
! Resolution du moment d'ordre 3
!*************************************************************

PM(:,NM3(JI))=PM(:,NM3(JI))+ &
              (PDMINTRA(:,NM3(JI))+PDMINTER(:,NM3(JI))+PDMCOND(:,NM3(JI))+&
               PSEDA(:,NM3(JI)))*PDT
PM(:,NM3(JI))= MAX(PM(:,NM3(JI)), XMNH_TINY)

!*************************************************************
! Resolution du moment d'ordre 6
!*************************************************************

PM(:,NM6(JI))=PM(:,NM6(JI))+ (PM(:,NM0(JI))**2*PDMINTRA(:,NM6(JI))+&
             PDMINTER(:,NM6(JI))+PDMCOND(:,NM6(JI)) + PSEDA(:,NM6(JI)) )*PDT

PM(:,NM6(JI))= MAX(PM(:,NM6(JI)), XMNH_TINY)
ENDDO


!*****************************************************************
!*****************************************************************
! SOLVEUR DE L'EQUILIBRE CHIMIQUE (MARS sera utilise)
!*****************************************************************
!*****************************************************************

!******************************************************************
! Calcul de la variation de concentration des differents
! composes pour trouver le nouveau moment d'ordre 3
!******************************************************************

DO JI=1,JPMODE

! Coagulation intermodale 
!-------------------------

DO JJ=1,NSP+NCARB+NSOA

 PCTOTA(:,JJ,JI)=PCTOTA(:,JJ,JI) &
                +(PCCTOT(:,JJ,1)*PDMINTER(:,NM3(JI)) + PCCTOT(:,JJ,JI)* PDMINTRA(:,NM3(JI))) &
                *XFAC(JJ)*PDT

! Sedimentation
!--------------
 PCTOTA(:,JJ,JI)=  PCTOTA(:,JJ,JI) + PCCTOT(:,JJ,JI)*PSEDA(:,NM3(JI))*XFAC(JJ)*PDT
 PCTOTA(:,JJ,JI)=  MAX(PCTOTA(:,JJ,JI), XMNH_TINY)


ENDDO
ENDDO

! H2SO4 Condensation + Nucleation 
!---------------------------------

 PCTOTA(:,JP_AER_SO4,1)=PCTOTA(:,JP_AER_SO4,1) &
                        +PDMCOND(:,NM3(1))*XFAC(JP_AER_SO4)*PDT
 PCTOTA(:,JP_AER_SO4,2)=PCTOTA(:,JP_AER_SO4,2) &
                        +PDMCOND(:,NM3(2))*XFAC(JP_AER_SO4)*PDT
!
!*************************************************************
! Calcul de la fraction massique entre les modes
!*************************************************************
ZSUM (:) = 0.
DO JI=1,JPMODE
  DO JJ=1,NSP+NCARB+NSOA
    ZSUM (:) = ZSUM (:) + PCTOTA(:,JJ,JI)
  ENDDO
ENDDO
POM(:,:) = 0.
DO JI=1,JPMODE
  DO JJ=1,NSP+NCARB+NSOA
    POM(:,JI)  =  POM(:,JI) + PCTOTA(:,JJ,JI) / ZSUM (:) 
  ENDDO
ENDDO
  

! Equilibre mineraux
!-------------------

IDT = INT(MAX(5.*PDT,1.))
IF ((PDT .GT. 0.).AND.( MOD(INT(PTIME) , IDT) .EQ. 0)) THEN
!IF (PDT .GT. 0.) THEN
  CALL CH_AER_MINERAL(PCTOTG, PCTOTA,PRV, PDENAIR, PPRESSURE, PTEMP, PRC, POM,&
                      PCCTOT,PSIG0, PRG0, PDT)

! Equilibre Organiques
!---------------------

  IF (NSOA .EQ. 10)   CALL CH_AER_ORGANIC(PCTOTG, PCTOTA,PRV, PDENAIR, &
                                          PPRESSURE, PTEMP,&
                                          PRC, POM, PCCTOT,PSIG0, PRG0, PDT, PSOLORG)
END IF

! Mass need to be positive
PCTOTA(:,:,:)= MAX (PCTOTA(:,:,:),0.)
PCTOTG(:,:)= MAX (PCTOTG(:,:),0.)

DO JI=1,JPMODE
  ZSUM(:)=0.
  DO JJ=1,NSP+NCARB+NSOA
   ZSUM(:)=ZSUM(:)+PCTOTA(:,JJ,JI)/XRHOI(JJ)
  ENDDO

  DO JJ=1,NSP+NCARB+NSOA
   PCCTOT(:,JJ,JI)=PCTOTA(:,JJ,JI)/XRHOI(JJ)/ZSUM(:)

  ENDDO
ENDDO
 
!******************************************************************************
! Calcul des nouveaux moments d'ordre 3 et 6
! Le moment d'ordre 3 est recalcule a partir de la composition de chaque mode
! Le moment d'ordre 6 est calcule pour garder sigma constant pendant l'equilibre chimique
!******************************************************************************
DO JN=1,JPMODE
!
  
  IF (JN .EQ. 1) THEN

    IF (LVARSIGI) THEN ! variable dispersion for mode 1

      ZSIGMA(:)=PM(:,NM3(JN))**2./(PM(:,NM0(JN))*PM(:,NM6(JN)))
      ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
      ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
      ZSIGMA(:)= LOG(ZSIGMA(:))
      ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
      WHERE (ZSIGMA(:) > XSIGIMAX)
      ZSIGMA(:) =  XSIGIMAX
      END WHERE
      WHERE (ZSIGMA(:) < XSIGIMIN)
      ZSIGMA(:) =  XSIGIMIN
      END WHERE

    ELSE ! fixed dispersion for mode 1
      ZSIGMA(:) = XINISIGI
    END IF
  END IF
!
  IF (JN .EQ. 2) THEN

    IF (LVARSIGJ) THEN ! variable dispersion for mode 2

      ZSIGMA(:)=PM(:,NM3(JN))**2./(PM(:,NM0(JN))*PM(:,NM6(JN)))
      ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
      ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
      ZSIGMA(:)= LOG(ZSIGMA(:))
      ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
      WHERE (ZSIGMA(:) > XSIGJMAX)
      ZSIGMA(:) =  XSIGJMAX
      END WHERE
      WHERE (ZSIGMA(:) < XSIGJMIN)
      ZSIGMA(:) =  XSIGJMIN
      END WHERE

    ELSE ! fixed dispersion for mode 2
      ZSIGMA(:) = XINISIGJ
    END IF
  END IF

 PSIG0(:,JN) = LOG(ZSIGMA(:))

 PN0(:,JN) = PM(:,NM0(JN))

END DO

DO JN=1,JPMODE
! Calcul du nouveau moment d'ordre 3
  ZNEWM(:,JN)=0.
  DO JJ=1,NSP+NCARB+NSOA
    PCTOTA(:,JJ,JN) = MAX(PCTOTA(:,JJ,JN),0.)
    ZNEWM(:,JN)=ZNEWM(:,JN)+PCTOTA(:,JJ,JN)/XFAC(JJ)
  ENDDO
  PM(:,NM3(JN))=ZNEWM(:,JN)
END DO

DO JN=1,JPMODE
  PM(:,NM6(JN)) = PM(:,NM0(JN)) &
          * ( (PM(:,NM3(JN))/PM(:,NM0(JN)))**(1./3.) * exp(-(3./2.)*PSIG0(:,JN)**2))**6 &
          * exp(18.*PSIG0(:,JN)**2)

  PRG0(:,JN)= (PM(:,NM3(JN))**4/(PM(:,NM6(JN)) * PM(:,NM0(JN))**3))**(1./6.)
ENDDO
!*************************************************************
! Blindages pour valeurs inferieurs au mininmum accepte
!*************************************************************
! ratio selon ch_aer_reallfin.f90 ((modifiable)
ZRATIOBC(1)  = 5.84E-3 / (5.84E-3+2.336E-2)
ZRATIOOC(1)  = 2.336E-2 / (5.84E-3+2.336E-2)
ZRATIOBC(2)  = 1.46E-3 / (1.46E-3+5.84E-3)
ZRATIOOC(2)  = 5.84E-3 / (1.46E-3+5.84E-3)
ZRGMIN(1) = ZINIRADIUSI / XCOEFRADIMIN
ZRGMIN(2) = ZINIRADIUSJ / XCOEFRADJMIN
ZRGMAX(1) = XCOEFRADIMAX * ZINIRADIUSI
ZRGMAX(2) = XCOEFRADJMAX * ZINIRADIUSJ

DO JN=1,JPMODE
  ZMASK(:,JN) = 1.
  WHERE ((PM(:,NM0(JN)) .LT. ZPMIN(NM0(JN))).OR.&
         (PM(:,NM3(JN)) .LT. ZPMIN(NM3(JN))).OR.&
         (PM(:,NM6(JN)) .LT. ZPMIN(NM6(JN))))
!         (PM(:,NM6(JN)) .LE. ZPMIN(NM6(JN))).OR.&
!         (PRG0(:,JN)) .LE. ZRGMIN(JN).OR.&
!         (PRG0(:,JN)) .GT. ZRGMAX(JN))

    PM(:,NM0(JN)) = ZPMIN(NM0(JN))
    PM(:,NM3(JN)) = ZPMIN(NM3(JN))
    PM(:,NM6(JN)) = ZPMIN(NM6(JN))

    ZMASK(:,JN)  = 0.
  END WHERE
  DO JJ=1,NSP+NCARB+NSOA
    PCTOTA(:,JJ,JN) = PCTOTA(:,JJ,JN) * ZMASK(:,JN)
  ENDDO
  WHERE (ZMASK(:,JN) == 0.)
    PCTOTA(:,JP_AER_BC,JN) = ZRATIOBC(JN) * ZPMIN(NM3(JN)) * XFAC(JP_AER_BC)
    PCTOTA(:,JP_AER_OC,JN) = ZRATIOOC(JN) * ZPMIN(NM3(JN)) * XFAC(JP_AER_OC)
  END WHERE

ENDDO
!
!
END SUBROUTINE CH_AER_SOLV
