!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/05/24 18:05:53
!-----------------------------------------------------------------
!     ############################################################
      MODULE MODI_PINTER
!     ############################################################
!
INTERFACE
      SUBROUTINE PINTER(PVMNH,PPMNH,PZGMNH,PTMNH,PVPL,PLPL,      &
                        KIU,KJU,KKU,KKB,KNP,HOPINT,HNAMV)
!
INTEGER,INTENT(IN)           :: KIU,KJU,KKU,KKB,KNP
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PVMNH
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PPMNH
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PZGMNH 
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PTMNH
REAL,DIMENSION(KIU,KJU,KNP),INTENT(OUT):: PVPL 
REAL,DIMENSION(KIU,KJU,KNP),INTENT(IN) :: PLPL
!
CHARACTER (LEN=3),INTENT(IN) :: HOPINT
CHARACTER (LEN=4),INTENT(IN) :: HNAMV
!
END SUBROUTINE PINTER
END INTERFACE
END MODULE MODI_PINTER
!
!------------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE PINTER(PVMNH,PPMNH,PZGMNH,PTMNH,PVPL,PLPL,       &
                        KIU,KJU,KKU,KKB,KNP,HOPINT,HNAMV)
!     #############################################################
!
!      PURPOSE
!!     -------
!  Interpolation verticale lineaire en Log(P) ou P
!  suivant HOPTINT, avec extrapolation style ECMWF.
!
!!**   METHOD
!!     ------
!!  Interpolation verticale lineaire en Log(P) ou P
!!  suivant OPTINT, avec extrapolation style ECMWF.
!!
!!  PVMNH  = tableau du champ donne au points masse Meso-NH
!!  PPMNH  = pressions au points masse Meso-NH
!!  PZGMNH = altitude geopotentiel au point masse Meso-NH
!!  PTMNH  = temperature au point masse Meso-NH
!!  PVPL   = tableau du champ interpole aux niveaux de pression
!!  PLPL   = liste des niveaux de pression demandes
!!  KIU    = nombre de points Meso-NH sur OX
!!  KJU    = nombre de points Meso-NH sur OY
!!  KKU    = nombre de niveaux Meso-NH (KKU=sommet, KKB=sol)
!!  KKB    = indice du niveau de reference CLS (1er niv. au dessus du sol)
!!  KNP    = nombre de niveaux de pression demandes (1=base, KNP=sommet)
!!  HOPINT = otion type d'interpolation ('LOG', ou 'LIN')
!!  HNAMV  = nom symbolique variable traite
!!          ('ZMGP','T...','RHU.','UB..','VB..','OMT.')
!!
!!
!!     EXTERNAL
!!     --------
!!      MODD_CST        XRD, XG
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!      None
!!
!!     REFERENCE
!!     ---------
!!      Research manual 2 ECMWF forecast model, 1988, Ref M1.6/3
!!      "adiabatic part", Appendix 6 postprocessing
!!      Section 3.  Vertical interpolation, p. A6.5-6
!!      Section 3.4 Extrapolation, pp. A6.6-7
!!
!!     AUTHOR
!!     ------
!!       P. Mascart     * LA *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original       22/04/96
!!       Chaboureau     26/04/04 add 2 dimension to PLPL
!!
!! ----------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_CST        ! XRD, XG
!
IMPLICIT NONE
!
!*       0.1  Declaration of arguments 
!
INTEGER,INTENT(IN)           :: KIU,KJU,KKU,KKB,KNP
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PVMNH
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PPMNH
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PZGMNH 
REAL,DIMENSION(KIU,KJU,KKU),INTENT(IN) :: PTMNH
REAL,DIMENSION(KIU,KJU,KNP),INTENT(OUT):: PVPL 
REAL,DIMENSION(KIU,KJU,KNP),INTENT(IN) :: PLPL
!
CHARACTER (LEN=3),INTENT(IN) :: HOPINT
CHARACTER (LEN=4),INTENT(IN) :: HNAMV
!
!*       0.2  Declaration of local variables
!
REAL      :: ZSLOPE,ZFP,ZTSTAR,ZTNUL,ZALF,ZZSOL,ZTPLAT,ZTPRNUL
INTEGER   :: JI,JJ,JKPL,JK
INTEGER   :: IKD
!
!         
!------------------------------------------------------------------------------
!
!*       1.   INTERPOLATION DES POINTS 
!             ------------------------
!
OX: DO  JI =1,KIU
  OY:  DO   JJ =1,KJU
    PLEV:  DO   JKPL=1,KNP
      !
      !   i) Reperage des zones
      !
      IKD=0
      IF(PLPL(JI,JJ,JKPL).LE.PPMNH(JI,JJ,KKU))             IKD=10*KKU
      DO  JK  =KKU-1,KKB,-1
         IF((PPMNH(JI,JJ,JK+1).LT.PLPL(JI,JJ,JKPL)).AND.   &
           (PLPL(JI,JJ,JKPL).LE.PPMNH(JI,JJ,JK)))          IKD=JK
      END DO
      IF(PLPL(JI,JJ,JKPL).GE.PPMNH(JI,JJ,KKB))             IKD=-10*KKU
      !
      !   ii) Interpolation des points reguliers
      !
      IF(ABS(IKD).NE.(10*KKU)) THEN
	IF(HOPINT.EQ.'LOG') THEN
	   ZSLOPE=LOG(PLPL(JI,JJ,JKPL)/PPMNH(JI,JJ,IKD))       &
                 /LOG(PPMNH(JI,JJ,IKD+1)/PPMNH(JI,JJ,IKD))
	ELSE
	   ZSLOPE=(PLPL(JI,JJ,JKPL)-PPMNH(JI,JJ,IKD))          &
                 /(PPMNH(JI,JJ,IKD+1)-PPMNH(JI,JJ,IKD))
	ENDIF
	PVPL(JI,JJ,JKPL)=PVMNH(JI,JJ,IKD)                &
                        +ZSLOPE*(PVMNH(JI,JJ,IKD+1)-PVMNH(JI,JJ,IKD))
      ENDIF
      !
      !  iii) Extrapolation base et sommet (algorithme ECMWF!)
      !
      !  Extrapolation au dessus du sommet
      !
      IF(IKD.EQ.10*KKU) THEN 
	IF(HNAMV.EQ.'ZMGP'.OR.HNAMV.EQ.'UB..'.OR.HNAMV.EQ.'VB..') THEN
          ! ZMGP ou UB ou VB
	  PVPL(JI,JJ,JKPL)=PVMNH(JI,JJ,KKU)+(PLPL(JI,JJ,JKPL)-PPMNH(JI,JJ,KKU))*  &
                          (PVMNH(JI,JJ,KKU-1)-PVMNH(JI,JJ,KKU))/            &
			  (PPMNH(JI,JJ,KKU-1)-PPMNH(JI,JJ,KKU))
	ELSE 
          ! autres (T ou OMT ou RHU)
	  PVPL(JI,JJ,JKPL)=PVMNH(JI,JJ,KKU) 
	ENDIF
      ENDIF
      !
      !   Extrapolation au dessous du sol
      !
      IF(IKD.EQ.-10*KKU) THEN
	IF(HNAMV.EQ.'ZMGP') THEN
           ! ZMGP
	   ZTSTAR=PTMNH(JI,JJ,KKB)
	   ZTNUL=ZTSTAR+.0065*PZGMNH(JI,JJ,KKB)
	   ZALF=.0065*XRD/XG
	   IF(ZTNUL.GE.290.5.AND.ZTSTAR.LT.290.5) THEN
	      ZALF=XRD*(290.5-ZTSTAR)/(PZGMNH(JI,JJ,KKB)*XG)
	   ELSE IF (ZTNUL.GT.290.5.AND.ZTSTAR.GT.290.5) THEN
	      ZALF=0.
	      ZTSTAR=.5*(ZTSTAR+290.5)
	   ENDIF
	   IF(ZTSTAR.LT.255) THEN
	      ZTSTAR=.5*(ZTSTAR+255.)
	   ENDIF
	   ZFP=ZALF*LOG(PLPL(JI,JJ,JKPL)/PPMNH(JI,JJ,KKB))
	   PVPL(JI,JJ,JKPL)=PVMNH(JI,JJ,KKB)-XRD*(ZTSTAR/XG)*     &
                            LOG(PLPL(JI,JJ,JKPL)/PPMNH(JI,JJ,KKB))*     &
                            (1.+ZFP/2.+(ZFP**2)/6.)
	ELSE IF(HNAMV.EQ.'T...') THEN
           ! T
	   ZTSTAR=PTMNH(JI,JJ,KKB)
	   IF(PPMNH(JI,JJ,KKB).GE.PLPL(JI,JJ,JKPL)) THEN
	      PVPL(JI,JJ,JKPL)=                                   &
              ((PPMNH(JI,JJ,KKB)-PLPL(JI,JJ,JKPL))*PTMNH(JI,JJ,KKB+1)   &
              +(PLPL(JI,JJ,JKPL)-PPMNH(JI,JJ,KKB+1))*ZTSTAR)            &
              /(PPMNH(JI,JJ,KKB)-PPMNH(JI,JJ,KKB+1))
	   ELSE IF(PLPL(JI,JJ,JKPL).GT.PPMNH(JI,JJ,KKB)) THEN
	      ZZSOL=PZGMNH(JI,JJ,KKB)
	      IF(ZZSOL.LT.2000.) THEN
		 ZALF=.0065*(XRD/XG)
	      ELSE
		 ZTNUL =ZTSTAR+.0065*ZZSOL
		 ZTPLAT=MIN(ZTNUL,298.)
		 IF(ZZSOL.GT.2500.) THEN
		    ZTPRNUL=ZTPLAT
		 ELSE 
		    ZTPRNUL=.002*((2500.-ZZSOL)*ZTNUL             &
		    +(ZZSOL-2000.)*ZTPLAT)
		 ENDIF
		 IF(ZTPRNUL.LT.ZTSTAR) THEN
		    ZALF=0.
		 ELSE
		    ZALF=XRD*(ZTPRNUL-ZTSTAR)/(PZGMNH(JI,JJ,KKB)*XG)
		 ENDIF
	      ENDIF
	      ZFP=ZALF*LOG(PLPL(JI,JJ,JKPL)/PPMNH(JI,JJ,KKB))
	      PVPL(JI,JJ,JKPL)=ZTSTAR*(1.+ZFP+(ZFP**2)/2.+(ZFP**3)/6.)
	   ENDIF
	ELSE
          ! autres (UB ou VB ou OMT ou RHU ...)
	  PVPL(JI,JJ,JKPL)=PVMNH(JI,JJ,KKB) 
	ENDIF
      ENDIF
    END DO PLEV
  END DO OY 
END DO OX
!
!
END SUBROUTINE PINTER
!------------------------------------------------------------------------------

