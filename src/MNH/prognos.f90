!MNH_LIC Copyright 2012-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #######################
      MODULE MODI_PROGNOS
!     #######################
!
INTERFACE
!
SUBROUTINE PROGNOS(PDT,PDZ,PLV,PCPH,PPRES,PRHOD,PRR,PTT,PRV,PRC,PS0,PCN,PCL)
!
REAL,                     INTENT(IN)    :: PDT
REAL, DIMENSION(:),       INTENT(IN)    :: PPRES
REAL, DIMENSION(:),       INTENT(IN)    :: PDZ
REAL, DIMENSION(:),       INTENT(IN)    :: PLV
REAL, DIMENSION(:),       INTENT(IN)    :: PCPH
REAL, DIMENSION(:),       INTENT(IN)    :: PRHOD
REAL, DIMENSION(:),       INTENT(IN)    :: PRR
REAL, DIMENSION(:),       INTENT(INOUT) :: PTT
REAL, DIMENSION(:),       INTENT(INOUT) :: PRV
REAL, DIMENSION(:),       INTENT(INOUT) :: PRC
REAL, DIMENSION(:),       INTENT(INOUT) :: PS0
REAL, DIMENSION(:),       INTENT(INOUT) :: PCN
REAL, DIMENSION(:),       INTENT(INOUT) :: PCL
!
END SUBROUTINE PROGNOS
!
END INTERFACE
!
END MODULE MODI_PROGNOS
!
!     ###################################################################################
      SUBROUTINE PROGNOS(PDT,PDZ,PLV,PCPH,PPRES,PRHOD,PRR,PTT,PRV,PRC,PS0,PCN,PCL)
!     ###################################################################################
!
!!****  * -  compute pseudo-prognostic of supersaturation according to Thouron
!                                                                     et al. 2012
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!
!!    REFERENCE
!!    ---------
!!
!!      Thouron, O., J.-L. Brenguier, and F. Burnet, Supersaturation calculation
!!      in large eddy simulation models for prediction of the droplet number
!!      concentration, Geosci. Model Dev., 5, 761-772, 2012.
!!
!!    AUTHOR
!!    ------
!!
!!      O. Thouron   * CNRM Meteo-France* : 
!!
!!    MODIFICATIONS
!!    -------------
!!     2014 G.Delautier : remplace MODD_RAIN_C2R2_PARAM par MODD_RAIN_C2R2_KHKO_PARAM
!!     2015 M.Mazoyer and O.Thouron : Physical tunings
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CST
USE MODD_PARAM_C2R2
USE MODD_RAIN_C2R2_KHKO_PARAM
!
USE MODE_IO_ll
USE MODE_MSG
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
REAL,                     INTENT(IN)    :: PDT
REAL, DIMENSION(:),       INTENT(IN)    :: PPRES
REAL, DIMENSION(:),       INTENT(IN)    :: PDZ
REAL, DIMENSION(:),       INTENT(IN)    :: PLV
REAL, DIMENSION(:),       INTENT(IN)    :: PCPH
REAL, DIMENSION(:),       INTENT(IN)    :: PRHOD
REAL, DIMENSION(:),       INTENT(IN)    :: PRR
REAL, DIMENSION(:),       INTENT(INOUT) :: PTT
REAL, DIMENSION(:),       INTENT(INOUT) :: PRV
REAL, DIMENSION(:),       INTENT(INOUT) :: PRC
REAL, DIMENSION(:),       INTENT(INOUT) :: PS0
REAL, DIMENSION(:),       INTENT(INOUT) :: PCN
REAL, DIMENSION(:),       INTENT(INOUT) :: PCL
!
!
!*       0.2   Declarations of local variables :
!
!
REAL, DIMENSION(SIZE(PRHOD,1))   ::  ZZW1,ZZW2,ZDZRC2,ZDZRC,ZCPH
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZA1,ZA2,ZB,ZC,ZG
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZLV,ZTT1,ZRT,ZTL,ZTT1_TEMP,ZTT2_TEMP
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZRMOY,ZRVSAT1,ZRVSAT2
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZVEC2  ! Work vectors forinterpolations
INTEGER, DIMENSION(SIZE(PRHOD,1)):: IVEC2   ! Vectors of indices for interpolations
INTEGER :: J1,J2
REAL,DIMENSION(SIZE(PS0,1))      ::MEM_PS0,ADJU2
REAL::AER_RAD
REAL, DIMENSION(SIZE(PRHOD,1))   :: ZFLAG_ACT   !Flag for activation
!
INTEGER                          :: IRESP      ! Return code of FM routines
INTEGER                          :: ILUOUT     ! Logical unit of output listing
CHARACTER(LEN=100)               :: YMSG
!
!minimum radius of cloud droplet
AER_RAD=1.0E-6
!
ZFLAG_ACT(:)=0.0
!ACTIVATION
ZVEC2(:) =0.0
IVEC2(:) =0.0
!
DO J1 = 1,4
  WHERE (PS0(:).GT.0.0)
   ZVEC2(:) = MAX( 1.00001, MIN( FLOAT(NHYP)-0.00001,      &
                   XHYPINTP1*LOG(PS0(:))+XHYPINTP2 ) )
   IVEC2(:) = INT( ZVEC2(:) )
   ZVEC2(:) = ZVEC2(:) - FLOAT( IVEC2(:) )
  END WHERE
END DO
ZZW1(:) =0.0
WHERE (PS0(:).GT.0.0)
ZZW1(:) =   XHYPF12( IVEC2(:)+1 )* ZVEC2(:)      &
             - XHYPF12( IVEC2(:)   )*(ZVEC2(:) - 1.0)
END WHERE
!
! the CCN spectra formula uses ZSMAX in percent
!
ZZW2(:)=0.0
IF (XCONC_CCN > 0) THEN
WHERE (PS0(:).GT.0.0)
   ZZW2(:) = MIN( XCONC_CCN,XCHEN * (100.0*PS0(:))**XKHEN * ZZW1(:) ) 
END WHERE
ELSE
WHERE (PS0(:).GT.0.0)
   ZZW2(:) = XCHEN * (100.0*PS0(:))**XKHEN * ZZW1(:) 
END WHERE
ENDIF
ZZW2(:)=MAX( (ZZW2(:)-PCN(:)),0.0 )
!
WHERE (ZZW2(:).LT.1.0)   !Non physique d'activer moins d'une particule
ZZW2(:)=0
END WHERE
!
!
!FLAG ACTIVE A TRUE (1.0) si on active pas
DO J2=1,SIZE(PRC,1)
 IF (ZZW2(J2).EQ.0.0) THEN
 ZFLAG_ACT(J2)=1.0
 ENDIF
ENDDO
!
! Mean radius
ZRMOY(:)=0.0
DO J2=1,SIZE(PRC,1)
 IF (PRC(J2).NE.0.0) THEN
 ZRMOY(J2)=(MOMG(XALPHAC,XNUC,3.0)*4.0*XPI*PCL(J2)*XRHOLW/&
        (3.0*PRC(J2)*PRHOD(J2)))**(1.0/3.0)
 ZRMOY(J2)=(PCL(J2)*MOMG(XALPHAC,XNUC,1.0)/ZRMOY(J2))
 ENDIF
 ZRMOY(J2)=ZRMOY(J2)+(ZZW2(J2)*AER_RAD)
ENDDO
!
PCL(:)= PCL(:)+ZZW2(:)
PCN(:)= PCN(:)+ZZW2(:)
!
!CALCUL DE A1 => Estimation de (drs/dt)f
!T(=à determiner) avant forcage; T'(=PTT) apres forcage
!Calcul de ZTT1: calculé en inversant S0(T)jusqu'à T:
! l'erreur faite sur cette inversion est supérieur à la précision
! recherchée, on applique à rs(T') pour cxalculer le DT=T'-T qui
! correspond à la variation rs(T')-rs(T). Permet de recuperer une valeur
! correcte de DT et donc de determiner T comme T=T'-DT         
!ZRVSAT1=rs(T)
!
ZRVSAT1(:)=PRV(:)/(PS0(:)+1.0)
!ZTT1<--es(T) de rs(T)
ZTT1_TEMP(:)=PPRES(:)*((((XMV / XMD)/ZRVSAT1(:))+1.0)**(-1D0))
!ZTT1<--T de es(T)
ZTT1_TEMP(:)=LOG(ZTT1_TEMP(:)/610.8)
ZTT1_TEMP(:)=(31.25*ZTT1_TEMP(:) -17.5688*273.15)/(ZTT1_TEMP(:) - 17.5688)
!es(T')
ZZW1(:)=EXP(XALPW-XBETAW/PTT(:)-XGAMW*LOG(PTT(:)))
!ZRVSAT2=rs(T')
ZRVSAT2(:)=(XMV / XMD)*ZZW1(:)/(PPRES(:)-ZZW1(:))
!ZTT2<--es(T') de rs(T')
ZTT2_TEMP(:)=PPRES(:)*((((XMV / XMD)/ZRVSAT2(:))+1.0)**(-1D0))
!ZTT2<--T' de es(T')
IF (MINVAL(ZTT2_TEMP).LT.0.0) THEN
  WRITE(YMSG,*) 'ZTT2_TEMP',MINVAL(ZTT2_TEMP),MINLOC(ZTT2_TEMP)
  CALL PRINT_MSG(NVERB_FATAL,'GEN','PROGNOS',YMSG)
ENDIF
!
ZTT2_TEMP(:)=LOG(ZZW1(:)/610.8)
ZTT2_TEMP(:)=(31.25*ZTT2_TEMP(:) -17.5688*273.15)/(ZTT2_TEMP(:) - 17.5688)
!ZTT1=T'-DT
ZTT1(:)=PTT(:)-(ZTT2_TEMP(:)-ZTT1_TEMP(:))
!Lv(T)
ZLV(:) = XLVTT+(XCPV-XCL)*(ZTT1(:)-XTT)                
!
ZA1(:)=-(((PS0(:)+1.0)**2.0)/PRV(:))*(ZRVSAT2(:)-(PRV(:)/(PS0(:)+1.0)))/PDT
!G
ZG(:)= 1.0/(XRHOLW*((XRV*ZTT1(:)/(XDIVA*EXP(XALPW-(XBETAW/ZTT1(:))-(XGAMW*LOG(ZTT1(:))))))     &
+((ZLV(:)/(XTHCO*ZTT1(:)))*((ZLV(:)/(ZTT1(:)*XRV))-1.0))))
!
ZC(:)=4.0*XPI*(XRHOLW/PRHOD(:))*ZG(:)
ZDZRC(:)=0.0
ZDZRC(:)=ZC(:)*PS0(:)*ZRMOY(:)
MEM_PS0(:)=PS0(:)
!CALCUL DE B => Estimation de (drs/dT)ce
!T(=PTT) avant condensation; T'(=à determiner) apres condensation
!Lv(T),Cph(T)
ZLV(:) = XLVTT+(XCPV-XCL)*(PTT(:)-XTT)                
ZCPH(:)= XCPD+XCPV*PRV(:)+XCL*(PRC(:)+PRR(:))
!T'=T+(DT)ce
ZTT1(:)=PTT(:)+(ZDZRC(:)*PDT*ZLV(:)/ZCPH(:))
!es(T')
ZZW1(:)=EXP(XALPW-XBETAW/PTT(:)-XGAMW*LOG(PTT(:)))
!rs(T')
ZZW1(:)=(XMV / XMD)*ZZW1(:)/(PPRES(:)-ZZW1(:))
!es(Tcond)
ZZW2(:)=EXP(XALPW-XBETAW/ZTT1(:)-XGAMW*LOG(ZTT1(:)))
!rs(Tcond)
ZZW2(:)=(XMV / XMD)*ZZW2(:)/(PPRES(:)-ZZW2(:))
!
WHERE (ZTT1(:).NE.PTT(:))
 ZB(:)=(ZLV(:)/ZCPH(:))*((ZZW2(:)-ZZW1(:))/(ZTT1(:)-PTT(:)))
ELSEWHERE
 ZB(:)=0.0
 ZDZRC(:)=0.0
ENDWHERE
!Calcul de S+dS
PS0(:)=PS0(:)+((ZA1(:)-(((ZB(:)*(PS0(:)+1.0)+1.0)*ZDZRC(:))/ZRVSAT1(:)))*PDT)
!
!Ajustement tel que rv=(s+1)*rvs
ZTL(:)=PTT(:)-(PLV(:)/PCPH(:))*PRC(:)
ZRT(:)=PRC(:)+PRV(:)
ZDZRC2(:)=PRC(:)
DO J2=1,SIZE(ZDZRC,1)
  IF ((ZDZRC(J2).NE.0.0).OR.(ZDZRC2(J2).NE.0.0)) THEN 
    DO J1=1,5
     ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
     ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
     ZZW1(J2)=EXP(XALPW-XBETAW/PTT(J2)-XGAMW*LOG(PTT(J2)))
     ZRVSAT1(J2)=(XMV / XMD)*ZZW1(J2)/(PPRES(J2)-ZZW1(J2))
     PRV(J2)=MIN(ZRT(J2),(PS0(J2)+1.0)*ZRVSAT1(J2))
     PRC(J2)=MAX(ZRT(J2)-PRV(J2),0.0)
     PTT(J2)=0.5*PTT(J2)+0.5*(ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2)))
    ENDDO
    ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
    ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
    PTT(J2)=ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2))
 ENDIF
ENDDO
ADJU2(:)=0.0
!
!Correction dans les mailles où ds a été surestimée
ZDZRC2(:)=PRC(:)-ZDZRC2(:)
WHERE ((MEM_PS0(:).LE.0.0).AND.(PS0(:).GT.0.0).AND.(ZDZRC2(:).LT.0.0))
  PS0(:)=0.0
  ADJU2(:)=1.0
ENDWHERE
!
WHERE ((MEM_PS0(:).GE.0.0).AND.(PS0(:).LT.0.0).AND.(ZDZRC2(:).GT.0.0))
  PS0(:)=0.0
  ADJU2(:)=1.0
ENDWHERE
!
DO J2=1,SIZE(ADJU2,1)
  IF (ADJU2(J2)==1) THEN 
   DO J1=1,5
    ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
    ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
    ZZW1(J2)=EXP(XALPW-XBETAW/PTT(J2)-XGAMW*LOG(PTT(J2)))
    ZRVSAT1(J2)=(XMV / XMD)*ZZW1(J2)/(PPRES(J2)-ZZW1(J2))
    PRV(J2)=MIN(ZRT(J2),(PS0(J2)+1.0)*ZRVSAT1(J2))
    PRC(J2)=MAX(ZRT(J2)-PRV(J2),0.0)
    PTT(J2)=0.5*PTT(J2)+0.5*(ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2)))
   ENDDO
   ZLV(J2) = XLVTT+(XCPV-XCL)*(PTT(J2)-XTT)                
   ZCPH(J2)=XCPD+XCPV*PRV(J2)+XCL*(PRC(J2)+PRR(J2))
   PTT(J2)=ZTL(J2)+(ZLV(J2)*PRC(J2)/ZCPH(J2))
  ENDIF
ENDDO
!
!Elimination de l'eau liquide dans les mailles où le rayon des gouttelettes est
!inférieur à AER_RAD
ZRMOY(:)=0.0
DO J2=1,SIZE(PRC,1)
 IF (PRC(J2).NE.0.0) THEN
  ZRMOY(J2)=(MOMG(XALPHAC,XNUC,3.0)*4.0*XPI*PCL(J2)*XRHOLW/&
        (3.0*PRC(J2)*PRHOD(J2)))**(1.0/3.0)
  ZRMOY(J2)=MOMG(XALPHAC,XNUC,1.0)/ZRMOY(J2)
IF ((ZFLAG_ACT(J2).EQ.1.0).AND.(MEM_PS0(J2).LT.0.0).AND.(ZRMOY(J2).LT.AER_RAD)) THEN
   PTT(J2)=ZTL(J2)
   PRV(J2)=ZRT(J2)
   PRC(J2)=0.0
  ENDIF
 ENDIF
ENDDO
!
!Calcul de S au regard de T et rv en fin de pas de temps
ZZW1=EXP(XALPW-XBETAW/PTT(:)-XGAMW*LOG(PTT(:)))
 !rvsat
ZRVSAT1(:)=(XMV / XMD)*ZZW1(:)/(PPRES-ZZW1(:))
!
WHERE (PRC(:)==0.0D0)
 PS0(:)=(PRV(:)/ZRVSAT1(:))-1D0
ENDWHERE
!
CONTAINS
!
FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
USE MODI_GAMMA
IMPLICIT NONE
REAL     :: PALPHA ! first shape parameter of the DIMENSIONnal distribution
REAL     :: PNU    ! second shape parameter of the DIMENSIONnal distribution
REAL     :: PP     ! order of the moment
REAL     :: PMOMG  ! result: moment of order ZP
PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
END FUNCTION MOMG
!
END SUBROUTINE PROGNOS
