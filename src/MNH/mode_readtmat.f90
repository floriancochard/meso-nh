!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ######spl
      MODULE MODE_READTMAT
!     ####################
!
!!****  *MODE_READTMAT* -  module routines  
!!
!!    PURPOSE
!!    -------
!!      Reads coefficients RE IM S11,S22 in Tmat lookup table, build from Tmat_20140121.f
!!      in rep /home/augros/Programmes/T-Matrice/Clotilde/ClotildeV3
!!      which is a modification of Mishchenko code
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    
!!
!!    AUTHOR
!!    ------
!!      C. Augros   * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   30/10/2012    
!!    
!!      C. Augros 15/11/2012 :  if ZLAM, D, T or ZELEV are out of the table limits
!!      S11, S22... are set to 0 => this is realistic only in case of small D, below Dmin
!!      but in reality it happens also for large D (> 8mm) => needs to be corrected
!!
!!      ----- MY_MODIF 4 ------
!!      C. Augros 7/12/2012 : quadrilinear interpolation instead of the ponderation by 1/d**2
!!      as done in ReadResuTmatV3.f90 (Programmes/T-Matrice/Clotilde/ClotildeV3/ReadResuTmatV3)
!!
!!      ----- MY_MODIF 5 ------
!!      C. Augros 12/12/2012 
!!	    SUBROUTINE READTMAT(PLAM_RAD,PELEV_RAD,PT,Dm,RES11,IMS11,&
!!      RES22,IMS22,RES11f,IMS11f,RES22f,IMS22f) instead of
!!	    SUBROUTINE READTMAT(PLAM_RAD,PELEV_RAD,PT,Dm,S11,S22,S11f,S22f)
!!
!!      ----- MY_MODIF 6 ------
!!      C. Augros 13/12/2012 
!!	    add of :
!!      SUBROUTINE READTMATINT(PLAM_RAD,PELEV_RAD,PT,M,PS11_CARRE,PS22_CARRE,&
!!      S22S11,PRE_S22FMS11F,PIM_S22FT,PIM_S11FT)
!!      read the 2nd Tmat table with the PS11_CARRE, PS22_CARRE... integrated
!!      over all diameters (table made by MakeTmatIntegree.f90)
!!      as done in ReadTmatInt (Programmes/T-Matrice/Clotilde/ClotildeV3/
!!      TmatIntegree/ReadTmatInt/)
!!
!!      C. Augros 01/2013
!!      Corrections dans subroutine READTMAT et READTMATINT des erreurs de lecture 
!!      des tables (remplacement de ZLAM_INF, ZTC_INF... par PLAM_MIN, PTC_MIN)
!!      
!!      ----- MY_MODIF 12 ------
!!      C. Augros 13/12/2012 
!!
!!	SUBROUTINE CALC_KTMAT(PLAM_RAD,PELEV_RAD,PT,M,PLAM_MIN,PLAM_MAX,PLAM_STEP,&
!!    	PELEV_MIN,PELEV_MAX,PELEV_STEP,PTC_MIN,PTC_MAX,PTC_STEP,&
!!    	KTMAT,PLAM_RED,PELEV_RED,PTC_RED,PM_RED)!!
!!	=> calcul des positions dans la table Tmat (ktmat) des coefficients à interpoler
!!	et des variables réduites qui traduisent la position entre 0 et 1 de ZLAM, ZELEV,
!!	ZTC et M par rapport aux bornes sup et inf
!!
!!      ------------ MY_MODIF 13 ------
!!      C. Augros 6/02/2014
!!      modification CALC_KTMAT et INTERPOL pour tenir compte de la colonne supplémentaire ZFW (wet graupel)
!!      uniquement pour les tables Tmat 06
!!      suppresion de la routine READTMATINT

!!      C. Augros 3/06/2014
!       Correction signe < dans mode_readtmat pour condition T° wet graupel
!       et modif dans calc_ktmat pour initialiser ZFW à 0 ou à max si proche des min et max

!!    ------------- MY_MODIF 14 ------
!     C. Augros 3/06/2014
!     modif seuils sur PR_RAY, PS_RAY...: on remplace par une condition sur M (hydromet content)> minM=10-7
!     pour être cohérent avec la valeur min dans les tables Tmat => evite plein de calculs pour rien
!
!     C. Augros 11/06/2014
!     Correction du bug dans CALC_KTMAT (mode_readtmat): INB_FW et pas ZFW !!!
!     => correction du pb dans la bande brillante

!!    ------------- MY_MODIF 15 ------
!     C. Augros 11/06/2014 (pas d'impact sur mode_readtmat)
!     => dans radar_scat: calcul fonction diélectique dans la glace avec maxwell garnett (et pas Smith 84)
!     si NDIFF=7 => calcul avec Mie pour la glace primaire
!     => dans write_lfifm1: ecriture des sorties sur une seule colonne F12.5 (on supprime az, porte)

!!    ------------- MY_MODIF 16 ------
!     C. Augros 22/08/2014
!     On ne conserve que les routines CALC_KTMAT et INTERPOL (qui sont appelées dans radar_scat si NDIFF=7)
!     Suppression de READTMAT (avec la nouvelle version de radarscat, les tables tmat
!     intégrées sont lues au début du code, pas besoin de les relire à chaque fois)
!-------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS

!              ------------
USE MODD_RADAR, ONLY:XVALGROUND
USE MODD_CST, ONLY: XPI
USE MODD_PARAMETERS, ONLY:XUNDEF,NUNDEF
!
!-------------------------------------------------------------------------------
!
CONTAINS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       1.   SUBROUTINE CALC_KTMAT
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   #################################################
    SUBROUTINE CALC_KTMAT(PELEV_RAD,PT,PFW,PM,&
    PELEV_MIN,PELEV_MAX,PELEV_STEP,PTC_MIN,PTC_MAX,PTC_STEP,&
    PFW_MIN,PFW_MAX,PFW_STEP,PEXPM_MIN,PEXPM_MAX,PEXPM_STEP,&
    KTMAT,PELEV_RED,PTC_RED,PFW_RED,PM_RED)
!   #################################################

IMPLICIT NONE
! Arguments entree et sortie
REAL,INTENT(IN) :: PELEV_RAD,PT,PFW,PM,PEXPM_MIN,PEXPM_MAX,PEXPM_STEP,&
PELEV_MIN,PELEV_MAX,PELEV_STEP,PTC_MIN,PTC_MAX,PTC_STEP,PFW_MIN,PFW_MAX,PFW_STEP
INTEGER,DIMENSION(16),INTENT(OUT):: KTMAT
REAL,INTENT(OUT) :: PELEV_RED,PTC_RED,PFW_RED,PM_RED
    
!************* Declarations ***************
REAL :: ZELEV,ZELEV_INF,ZELEV_SUP
REAL :: ZTC,ZTC_INF,ZTC_SUP
REAL :: ZFW_INF,ZFW_SUP,ZFW
REAL :: ZEXPM,ZEXPM_INF,ZEXPM_SUP,ZEXPM_RED,ZM_INF,ZM_SUP
INTEGER :: IELEV,ITC,IFW,IEXPM
INTEGER :: INB_ELEV,INB_TC,INB_FW,INB_M
INTEGER :: IELEVS,ITCS,IEXPMS,IFWS
!**********************************************
!*****           Parametres !!!          ******
!**********************************************
!Conversion de l'elevation en degre (a partir de la valeur en radian)
ZELEV=PELEV_RAD*180./XPI
!Conversion de la temperature de °K en °C
ZTC=PT-273.15 
!Hydromet content
ZEXPM=LOG10(PM)  
!Liquid water fraction
ZFW=PFW


!On verifie que  ZELEV, ZTC, ZFW et M sont compris dans les bornes min et max

IF (ABS(ZELEV-PELEV_MIN) < PELEV_STEP/10) ZELEV=PELEV_MIN
IF (ABS(ZELEV-PELEV_MAX) < PELEV_STEP/10) ZELEV=PELEV_MAX
IF (ABS(ZTC-PTC_MIN) < PTC_STEP/10) ZTC=PTC_MIN
IF (ABS(ZTC-PTC_MAX) < PTC_STEP/10) ZTC=PTC_MAX
IF (ABS(ZFW-PFW_MIN) < PFW_STEP/10) ZFW=PFW_MIN
IF (ABS(ZFW-PFW_MAX) < PFW_STEP/10) ZFW=PFW_MAX
IF (ABS(ZEXPM-PEXPM_MIN) < PEXPM_STEP/10) ZEXPM=PEXPM_MIN
IF (ABS(ZEXPM-PEXPM_MAX) < PEXPM_STEP/10) ZEXPM=PEXPM_MAX

IF ((ZELEV >=PELEV_MIN).AND. (ZELEV<=PELEV_MAX) .AND.&
    (ZTC >=PTC_MIN)    .AND. (ZTC<=PTC_MAX)   .AND.(ZFW >=PFW_MIN)    .AND. (ZFW<=PFW_MAX)     .AND.&
    (ZEXPM >=PEXPM_MIN).AND. (ZEXPM<=PEXPM_MAX)) THEN
    
  !Recherche dans le fichier de la position des valeurs encadrant les
  !valeurs données ci-dessus
    !------- ZELEV ------------------
  IELEV=floor((ZELEV-PELEV_MIN)/PELEV_STEP)
  ZELEV_INF=PELEV_MIN+IELEV*PELEV_STEP
  IF (ZELEV==ZELEV_INF) THEN
    IELEVS=IELEV
  ELSE
    IELEVS=IELEV+1
  ENDIF
  ZELEV_SUP=PELEV_MIN+IELEVS*PELEV_STEP
  INB_ELEV=nint((PELEV_MAX-PELEV_MIN)/PELEV_STEP)+1
  !WRITE(0,*) "IELEV,IELEVS,ZELEV_INF,ZELEV_SUP,INB_ELEV : ",IELEV,IELEVS,ZELEV_INF,ZELEV_SUP,INB_ELEV
  !------- ZTC ------------------
  ITC=floor((ZTC-PTC_MIN)/PTC_STEP)
  ZTC_INF=PTC_MIN+ITC*PTC_STEP
  IF (ZTC==ZTC_INF) THEN
    ITCS=ITC
  ELSE
    ITCS=ITC+1
  ENDIF
  ZTC_SUP=PTC_MIN+(ITCS)*PTC_STEP
  INB_TC=nint((PTC_MAX-PTC_MIN)/PTC_STEP)+1
  !WRITE(0,*) "ITC,ITCS,ZTC_INF,ZTC_SUP,INB_TC : ",ITC,ITCS,ZTC_INF,ZTC_SUP,INB_TC

  !------- ZFW ------------------
  IFW=floor((ZFW-PFW_MIN)/PFW_STEP)
  ZFW_INF=PFW_MIN+IFW*PFW_STEP
  IF (ZFW==ZFW_INF) THEN
    IFWS=IFW
  ELSE
    IFWS=IFW+1
  ENDIF
  ZFW_SUP=PFW_MIN+(IFWS)*PFW_STEP
  INB_FW=nint((PFW_MAX-PFW_MIN)/PFW_STEP)+1
  !WRITE(0,*) "IFW,IFWS,ZFW_INF,ZFW_SUP,INB_FW : ",IFW,IFWS,ZFW_INF,ZFW_SUP,INB_FW

  !------- PM ------------------
  IEXPM=floor((ZEXPM-PEXPM_MIN)/PEXPM_STEP)
  ZEXPM_INF=PEXPM_MIN+IEXPM*PEXPM_STEP
  IF (ZEXPM==ZEXPM_INF) THEN
    IEXPMS=IEXPM
  ELSE
    IEXPMS=IEXPM+1
  ENDIF
  ZEXPM_SUP=PEXPM_MIN+IEXPMS*PEXPM_STEP
  INB_M=nint((PEXPM_MAX-PEXPM_MIN)/PEXPM_STEP)+1
  ZM_INF=10**ZEXPM_INF
  ZM_SUP=10**ZEXPM_SUP
  !WRITE(0,*) "IEXPM,IEXPMS,ZEXPM_INF,ZEXPM_SUP,INB_M,ZM_INF,ZM_SUP : ",&
  !IEXPM,IEXPMS,ZEXPM_INF,ZEXPM_SUP,INB_M,ZM_INF,ZM_SUP
  !WRITE(0,*) "  "

  !-- Calcul des variables reduites (comprises entre 0 et 1)
  ! pour l'interpolation linaire
  IF (ZELEV_SUP .NE. ZELEV_INF) THEN 
    PELEV_RED=(ZELEV-ZELEV_INF)/(ZELEV_SUP-ZELEV_INF)
  ELSE
    PELEV_RED=0
  ENDIF
  IF (ZTC_SUP .NE. ZTC_INF) THEN 
    PTC_RED=(ZTC-ZTC_INF)/(ZTC_SUP-ZTC_INF)
  ELSE
    PTC_RED=0
  ENDIF
  IF (ZFW_SUP .NE. ZFW_INF) THEN 
    PFW_RED=(ZFW-ZFW_INF)/(ZFW_SUP-ZFW_INF)
  ELSE
    PFW_RED=0
  ENDIF
  IF (ZEXPM_SUP .NE. ZEXPM_INF) THEN 
    PM_RED=(PM-ZM_INF)/(ZM_SUP-ZM_INF)
    ZEXPM_RED=(ZEXPM-ZEXPM_INF)/(ZEXPM_SUP-ZEXPM_INF)
  ELSE
    PM_RED=0
  ENDIF
  KTMAT(1)=ITC*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFW*INB_M+IEXPM+1
  KTMAT(2)=ITC*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFW*INB_M+IEXPMS+1
  KTMAT(3)=ITC*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFWS*INB_M+IEXPM+1
  KTMAT(4)=ITC*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFWS*INB_M+IEXPMS+1
  KTMAT(5)=ITC*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFW*INB_M+IEXPM+1
  KTMAT(6)=ITC*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFW*INB_M+IEXPMS+1
  KTMAT(7)=ITC*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFWS*INB_M+IEXPM+1
  KTMAT(8)=ITC*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFWS*INB_M+IEXPMS+1
  KTMAT(9)=ITCS*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFW*INB_M+IEXPM+1
  KTMAT(10)=ITCS*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFW*INB_M+IEXPMS+1
  KTMAT(11)=ITCS*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFWS*INB_M+IEXPM+1
  KTMAT(12)=ITCS*INB_ELEV*INB_FW*INB_M+IELEV*INB_FW*INB_M+IFWS*INB_M+IEXPMS+1
  KTMAT(13)=ITCS*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFW*INB_M+IEXPM+1
  KTMAT(14)=ITCS*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFW*INB_M+IEXPMS+1
  KTMAT(15)=ITCS*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFWS*INB_M+IEXPM+1
  KTMAT(16)=ITCS*INB_ELEV*INB_FW*INB_M+IELEVS*INB_FW*INB_M+IFWS*INB_M+IEXPMS+1
ELSE
!  WRITE(0,*) "ZM, ZTC, ZELEV ou en dehors des bornes:"
!  WRITE(0,*) ",ZELEV,ZTC,ZEXPM, ZFW : ",ZELEV,ZTC,ZEXPM, ZFW
!  WRITE(0,*) "PELEV_MIN,PELEV_STEP,PELEV_MAX",PELEV_MIN,PELEV_STEP,PELEV_MAX
!  WRITE(0,*) "PTC_MIN,PTC_STEP,PTC_MAX",PTC_MIN,PTC_STEP,PTC_MAX
!  WRITE(0,*) "PFW_MIN,PFW_STEP,PFW_MAX",PFW_MIN,PFW_STEP,PFW_MAX
!  WRITE(0,*) "PEXPM_MIN,PEXPM_STEP,PEXPM_MAX",PEXPM_MIN,PEXPM_STEP,PEXPM_MAX
!  WRITE(0,*) "--------------------------------"
!  IF ((ZELEV >=PELEV_MIN).AND. (ZELEV<=PELEV_MAX)) THEN
!    WRITE(0,*) "ok ZELEV :",ZELEV
!  ELSE
!    WRITE(0,*) "Nok ZELEV :",ZELEV
!  ENDIF
!  IF ((ZTC >=PTC_MIN).AND. (ZTC<=PTC_MAX)) THEN
!    WRITE(0,*) "ok ZTC :",ZTC
!  ELSE
!    WRITE(0,*) "Nok ZTC :",ZTC
!  ENDIF
!  IF ((ZFW >=PFW_MIN).AND. (ZFW<=PFW_MAX)) THEN
!    WRITE(0,*) "ok ZFW :",ZFW
!  ELSE
!    WRITE(0,*) "Nok ZFW :",ZFW
!  ENDIF
!  IF ((ZEXPM >=PEXPM_MIN).AND. (ZEXPM<=PEXPM_MAX)) THEN
!    WRITE(0,*) "ok ZEXPM :",ZEXPM
!  ELSE
!    WRITE(0,*) "Nok ZEXPM :",ZEXPM
!  ENDIF
  KTMAT(:)=-NUNDEF
  PTC_RED=-XUNDEF
  PELEV_RED=-XUNDEF
  PFW_RED=-XUNDEF
  PM_RED=-XUNDEF
ENDIF

RETURN
END SUBROUTINE CALC_KTMAT

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       2.   SUBROUTINE INTERPOL
!         --------------------------------------
!-------------------------------------------------------------------------------
!   #################################################
  SUBROUTINE INTERPOL(PELEV_RED,PTC_RED,PFW_RED,PM_RED,PMAT_COEF,&
  PS11_CARRE,PS22_CARRE,PRE_S22S11,PIM_S22S11,PRE_S22FMS11F,PIM_S22FT,PIM_S11FT)
!   #################################################

IMPLICIT NONE
! Arguments entree et sortie
REAL,INTENT(IN) :: PELEV_RED,PTC_RED,PFW_RED,PM_RED
REAL,DIMENSION(7,16),INTENT(IN) :: PMAT_COEF !matrice contenant tous les coef interpolés: RES11, RES22...
REAL,INTENT(OUT) :: PS11_CARRE,PS22_CARRE,PRE_S22S11,PIM_S22S11,PRE_S22FMS11F,PIM_S22FT,PIM_S11FT 
INTEGER:: JCOEF
REAL,DIMENSION(7) :: ZVECT_COEF !vecteur contenant tous les coef interpolés: RES11, RES22...

!------------------- Interpolation --------------------------
!on calcule la distance entre le point ZELEV,T,M avec les 
!bornes inf et sup
! si on ne se trouve pas exactement sur une des bornes
!IF ((ZELEV_SUP/=ZELEV_INF) .OR. (ZTC_SUP/=ZTC_INF) .OR. (ZM_SUP/=ZM_INF)) THEN
!	WRITE(0,*) "IF ( (ZELEV_SUP/=ZELEV_INF) .OR. (ZTC_SUP/=ZTC_INF) .OR. (ZM_SUP/=ZM_INF))"
!ENDIF

   !--- Interpolation linéaire ---
DO JCOEF=1,7,1
  ZVECT_COEF(JCOEF)= & 
    (1-PM_RED)*(1-PFW_RED)*(1-PELEV_RED)*(1-PTC_RED)*PMAT_COEF(JCOEF,1)+&
    PM_RED*(1-PFW_RED)*(1-PELEV_RED)*(1-PTC_RED)*PMAT_COEF(JCOEF,2)+&
    (1-PM_RED)*PFW_RED*(1-PELEV_RED)*(1-PTC_RED)*PMAT_COEF(JCOEF,3)+&
    PM_RED*PFW_RED*(1-PELEV_RED)*(1-PTC_RED)*PMAT_COEF(JCOEF,4)+&
    (1-PM_RED)*(1-PFW_RED)*PELEV_RED*(1-PTC_RED)*PMAT_COEF(JCOEF,5)+&
    PM_RED*(1-PFW_RED)*PELEV_RED*(1-PTC_RED)*PMAT_COEF(JCOEF,6)+&
    (1-PM_RED)*PFW_RED*PELEV_RED*(1-PTC_RED)*PMAT_COEF(JCOEF,7)+&
    PM_RED*PFW_RED*PELEV_RED*(1-PTC_RED)*PMAT_COEF(JCOEF,8)+&

    (1-PM_RED)*(1-PFW_RED)*(1-PELEV_RED)*PTC_RED*PMAT_COEF(JCOEF,9)+&
    PM_RED*(1-PFW_RED)*(1-PELEV_RED)*PTC_RED*PMAT_COEF(JCOEF,10)+&
    (1-PM_RED)*PFW_RED*(1-PELEV_RED)*PTC_RED*PMAT_COEF(JCOEF,11)+&
    PM_RED*PFW_RED*(1-PELEV_RED)*PTC_RED*PMAT_COEF(JCOEF,12)+&
    (1-PM_RED)*(1-PFW_RED)*PELEV_RED*PTC_RED*PMAT_COEF(JCOEF,13)+&
    PM_RED*(1-PFW_RED)*PELEV_RED*PTC_RED*PMAT_COEF(JCOEF,14)+&
    (1-PM_RED)*PFW_RED*PELEV_RED*PTC_RED*PMAT_COEF(JCOEF,15)+&
    PM_RED*PFW_RED*PELEV_RED*PTC_RED*PMAT_COEF(JCOEF,16)

ENDDO!--- fin interpolation lineaire -----------------------

PS11_CARRE=ZVECT_COEF(1)
PS22_CARRE=ZVECT_COEF(2)
PRE_S22S11=ZVECT_COEF(3)
PIM_S22S11=ZVECT_COEF(4)
PRE_S22FMS11F=ZVECT_COEF(5)
PIM_S22FT=ZVECT_COEF(6)
PIM_S11FT=ZVECT_COEF(7)

RETURN

  END SUBROUTINE INTERPOL
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       3.   SUBROUTINE CALC_KTMAT_LIMA  : idem CALC_KTMAT mais pour adapter 
!             aux concentrations pour l'espace pluie pour LIMA
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   #################################################
    SUBROUTINE CALC_KTMAT_LIMA(PELEV_RAD,PT,PCC,PM,&
    PELEV_MIN,PELEV_MAX,PELEV_STEP,PTC_MIN,PTC_MAX,PTC_STEP,&
    PEXPCC_MIN,PEXPCC_MAX,PEXPCC_STEP,PEXPM_MIN,PEXPM_MAX,PEXPM_STEP,&
    KTMAT,PELEV_RED,PTC_RED,PCC_RED,PM_RED)
!   #################################################

IMPLICIT NONE
! Arguments entree et sortie
REAL,INTENT(IN) :: PELEV_RAD,PT,PCC,PM,PEXPM_MIN,PEXPM_MAX,PEXPM_STEP,&
PELEV_MIN,PELEV_MAX,PELEV_STEP,PTC_MIN,PTC_MAX,PTC_STEP,PEXPCC_MIN,PEXPCC_MAX,PEXPCC_STEP
INTEGER,DIMENSION(16),INTENT(OUT):: KTMAT
REAL,INTENT(OUT) :: PELEV_RED,PTC_RED,PCC_RED,PM_RED
    
!************* Declarations ***************
REAL :: ZELEV,ZELEV_INF,ZELEV_SUP
REAL :: ZTC,ZTC_INF,ZTC_SUP
REAL :: ZEXPCC_INF,ZEXPCC_SUP,ZEXPCC,ZCC_INF,ZCC_SUP,ZEXPCC_RED
REAL :: ZEXPM,ZEXPM_INF,ZEXPM_SUP,ZEXPM_RED,ZM_INF,ZM_SUP
INTEGER :: IELEV,ITC,IEXPCC,IEXPM
INTEGER :: INB_ELEV,INB_TC,INB_CC,INB_M
INTEGER :: IELEVS,ITCS,IEXPMS,IEXPCCS
!**********************************************
!*****           Parametres !!!          ******
!**********************************************

!Conversion de l'elevation en degre (a partir de la valeur en radian)
ZELEV=PELEV_RAD*180./XPI
!Conversion de la temperature de °K en °C
ZTC=PT-273.15 
!Hydromet content
ZEXPM=LOG10(PM)  
!Concentration
ZEXPCC=LOG10(PCC)


!On verifie que  ZELEV, ZTC, ZCC et M sont compris dans les bornes min et max
IF (ABS(ZELEV-PELEV_MIN) < PELEV_STEP/10) ZELEV=PELEV_MIN
IF (ABS(ZELEV-PELEV_MAX) < PELEV_STEP/10) ZELEV=PELEV_MAX
IF (ABS(ZTC-PTC_MIN) < PTC_STEP/10) ZTC=PTC_MIN
IF (ABS(ZTC-PTC_MAX) < PTC_STEP/10) ZTC=PTC_MAX
IF (ABS(ZEXPCC-PEXPCC_MIN) < PEXPCC_STEP/10) ZEXPCC=PEXPCC_MIN
IF (ABS(ZEXPCC-PEXPCC_MAX) < PEXPCC_STEP/10) ZEXPCC=PEXPCC_MAX
IF (ABS(ZEXPM-PEXPM_MIN) < PEXPM_STEP/10) ZEXPM=PEXPM_MIN
IF (ABS(ZEXPM-PEXPM_MAX) < PEXPM_STEP/10) ZEXPM=PEXPM_MAX

IF ((ZELEV >=PELEV_MIN).AND. (ZELEV<=PELEV_MAX) .AND.&
    (ZTC >=PTC_MIN)    .AND. (ZTC<=PTC_MAX)   .AND.(ZEXPCC >=PEXPCC_MIN)    .AND. (ZEXPCC<=PEXPCC_MAX)     .AND.&
    (ZEXPM >=PEXPM_MIN).AND. (ZEXPM<=PEXPM_MAX)) THEN
    
  !Recherche dans le fichier de la position des valeurs encadrant les
  !valeurs données ci-dessus
  !------- ZELEV ------------------
  IELEV=floor((ZELEV-PELEV_MIN)/PELEV_STEP)
  ZELEV_INF=PELEV_MIN+IELEV*PELEV_STEP
  IF (ZELEV==ZELEV_INF) THEN
    IELEVS=IELEV
  ELSE
    IELEVS=IELEV+1
  ENDIF
  ZELEV_SUP=PELEV_MIN+IELEVS*PELEV_STEP
  INB_ELEV=nint((PELEV_MAX-PELEV_MIN)/PELEV_STEP)+1
  !WRITE(0,*) "IELEV,IELEVS,ZELEV_INF,ZELEV_SUP,INB_ELEV : ",IELEV,IELEVS,ZELEV_INF,ZELEV_SUP,INB_ELEV
  !------- ZTC ------------------
  ITC=floor((ZTC-PTC_MIN)/PTC_STEP)
  ZTC_INF=PTC_MIN+ITC*PTC_STEP
  IF (ZTC==ZTC_INF) THEN
    ITCS=ITC
  ELSE
    ITCS=ITC+1
  ENDIF
  ZTC_SUP=PTC_MIN+(ITCS)*PTC_STEP
  INB_TC=nint((PTC_MAX-PTC_MIN)/PTC_STEP)+1
  !WRITE(0,*) "ITC,ITCS,ZTC_INF,ZTC_SUP,INB_TC : ",ITC,ITCS,ZTC_INF,ZTC_SUP,INB_TC

  !------- ZCC ------------------
  IEXPCC=floor((ZEXPCC-PEXPCC_MIN)/PEXPCC_STEP)
  ZEXPCC_INF=PEXPCC_MIN+IEXPCC*PEXPCC_STEP
  IF (ZEXPCC==ZEXPCC_INF) THEN
    IEXPCCS=IEXPCC
  ELSE
    IEXPCCS=IEXPCC+1
  ENDIF
  ZEXPCC_SUP=PEXPCC_MIN+(IEXPCCS)*PEXPCC_STEP
  INB_CC=nint((PEXPCC_MAX-PEXPCC_MIN)/PEXPCC_STEP)+1
  ZCC_INF=10**ZEXPCC_INF
  ZCC_SUP=10**ZEXPCC_SUP
  !------- PM ------------------
  IEXPM=floor((ZEXPM-PEXPM_MIN)/PEXPM_STEP)
  ZEXPM_INF=PEXPM_MIN+IEXPM*PEXPM_STEP
  IF (ZEXPM==ZEXPM_INF) THEN
    IEXPMS=IEXPM
  ELSE
    IEXPMS=IEXPM+1
  ENDIF
  ZEXPM_SUP=PEXPM_MIN+IEXPMS*PEXPM_STEP
  INB_M=nint((PEXPM_MAX-PEXPM_MIN)/PEXPM_STEP)+1
  ZM_INF=10**ZEXPM_INF
  ZM_SUP=10**ZEXPM_SUP
  !WRITE(0,*) "IEXPM,IEXPMS,ZEXPM_INF,ZEXPM_SUP,INB_M,ZM_INF,ZM_SUP : ",&
  !IEXPM,IEXPMS,ZEXPM_INF,ZEXPM_SUP,INB_M,ZM_INF,ZM_SUP
  !WRITE(0,*) "  "

  !-- Calcul des variables reduites (comprises entre 0 et 1)
  ! pour l'interpolation linaire
  IF (ZELEV_SUP .NE. ZELEV_INF) THEN 
    PELEV_RED=(ZELEV-ZELEV_INF)/(ZELEV_SUP-ZELEV_INF)
  ELSE
    PELEV_RED=0
  ENDIF
  IF (ZTC_SUP .NE. ZTC_INF) THEN 
    PTC_RED=(ZTC-ZTC_INF)/(ZTC_SUP-ZTC_INF)
  ELSE
    PTC_RED=0
  ENDIF
  IF (ZEXPCC_SUP .NE. ZEXPCC_INF) THEN 
    PCC_RED=(PCC-ZCC_INF)/(ZCC_SUP-ZCC_INF)
    ZEXPCC_RED=(ZEXPCC-ZEXPCC_INF)/(ZEXPCC_SUP-ZEXPCC_INF)
  ELSE
    PCC_RED=0
  ENDIF
  IF (ZEXPM_SUP .NE. ZEXPM_INF) THEN 
    PM_RED=(PM-ZM_INF)/(ZM_SUP-ZM_INF)
    ZEXPM_RED=(ZEXPM-ZEXPM_INF)/(ZEXPM_SUP-ZEXPM_INF)
  ELSE
    PM_RED=0
  ENDIF
  KTMAT(1)=ITC*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCC*INB_M+IEXPM+1
  KTMAT(2)=ITC*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCC*INB_M+IEXPMS+1
  KTMAT(3)=ITC*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCCS*INB_M+IEXPM+1
  KTMAT(4)=ITC*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCCS*INB_M+IEXPMS+1
  KTMAT(5)=ITC*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCC*INB_M+IEXPM+1
  KTMAT(6)=ITC*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCC*INB_M+IEXPMS+1
  KTMAT(7)=ITC*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCCS*INB_M+IEXPM+1
  KTMAT(8)=ITC*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCCS*INB_M+IEXPMS+1
  KTMAT(9)=ITCS*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCC*INB_M+IEXPM+1
  KTMAT(10)=ITCS*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCC*INB_M+IEXPMS+1
  KTMAT(11)=ITCS*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCCS*INB_M+IEXPM+1
  KTMAT(12)=ITCS*INB_ELEV*INB_CC*INB_M+IELEV*INB_CC*INB_M+IEXPCCS*INB_M+IEXPMS+1
  KTMAT(13)=ITCS*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCC*INB_M+IEXPM+1
  KTMAT(14)=ITCS*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCC*INB_M+IEXPMS+1
  KTMAT(15)=ITCS*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCCS*INB_M+IEXPM+1
  KTMAT(16)=ITCS*INB_ELEV*INB_CC*INB_M+IELEVS*INB_CC*INB_M+IEXPCCS*INB_M+IEXPMS+1
ELSE
!  WRITE(0,*) "ZM, ZTC, ZELEV en dehors des bornes:"
!  WRITE(0,*) "ZELEV,ZTC,ZEXPM, ZEXPCC : ",ZELEV,ZTC,ZEXPM, ZEXPCC
!  WRITE(0,*) "PELEV_MIN,PELEV_STEP,PELEV_MAX",PELEV_MIN,PELEV_STEP,PELEV_MAX
!  WRITE(0,*) "PTC_MIN,PTC_STEP,PTC_MAX",PTC_MIN,PTC_STEP,PTC_MAX
!  WRITE(0,*) "PEXPCC_MIN,PEXPCC_STEP,PEXPCC_MAX",PEXPCC_MIN,PEXPCC_STEP,PEXPCC_MAX
!  WRITE(0,*) "PEXPM_MIN,PEXPM_STEP,PEXPM_MAX",PEXPM_MIN,PEXPM_STEP,PEXPM_MAX
!  WRITE(0,*) "--------------------------------"
!  IF ((ZELEV >=PELEV_MIN).AND. (ZELEV<=PELEV_MAX)) THEN
!    WRITE(0,*) "ok ZELEV :",ZELEV
!  ELSE
!    WRITE(0,*) "Nok ZELEV :",ZELEV
!  ENDIF
!  IF ((ZTC >=PTC_MIN).AND. (ZTC<=PTC_MAX)) THEN
!    WRITE(0,*) "ok ZTC :",ZTC
!  ELSE
!    WRITE(0,*) "Nok ZTC :",ZTC
!  ENDIF
!  IF ((ZEXPCC >=PEXPCC_MIN).AND. (ZEXPCC<=PEXPCC_MAX)) THEN
!    WRITE(0,*) "ok ZEXPCC :",ZEXPCC
!  ELSE
!    WRITE(0,*) "Nok ZEXPCC :",ZEXPCC
!  ENDIF
!  IF ((ZEXPM >=PEXPM_MIN).AND. (ZEXPM<=PEXPM_MAX)) THEN
!    WRITE(0,*) "ok ZEXPM :",ZEXPM
!  ELSE
!    WRITE(0,*) "Nok ZEXPM :",ZEXPM
!  ENDIF
  KTMAT(:)=-NUNDEF
  PTC_RED=-XUNDEF
  PELEV_RED=-XUNDEF
  PCC_RED=-XUNDEF
  PM_RED=-XUNDEF
ENDIF

RETURN
END SUBROUTINE CALC_KTMAT_LIMA



END MODULE MODE_READTMAT
