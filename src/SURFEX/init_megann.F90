!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!    ###############################
SUBROUTINE INIT_MEGAN_n(IO, S, K, NP, MSF, MGN, PLAT, HSV, PMEGAN_FIELDS)
!    ###############################
!!
!!***  *BVOCEM*
!! 
!!    PURPOSE
!!    -------
!!    Calculate the biogenic emission fluxes upon the MEGAN code
!!    http://lar.wsu.edu/megan/
!!
!!    METHOD
!!    ------
!!
!!
!!    AUTHOR
!!    ------
!!    P. Tulet (LACy)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original: 25/10/14
!!    Modified: 06/2017, J. Pianezze, adaptation for SurfEx v8.0
!!    Modified: 06/2018, P. Tulet,  add PFT and LAI
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!
USE MODD_MEGAN_SURF_FIELDS_n, ONLY : MEGAN_SURF_FIELDS_t
USE MODD_MEGAN_n, ONLY : MEGAN_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_P_t, ISBA_K_t, ISBA_NP_t
!
USE MODD_DATA_COVER_PAR, ONLY : NVT_C4, NVT_TRBE, NVT_TRBD, NVT_TEBE, &
                    NVT_TEBD, NVT_TENE, NVT_BOBD, NVT_BONE, NVT_BOND, &
                    NVT_BOGR, NVT_SHRB, NVT_GRAS, NVT_TROG, NVT_C3,   &
                    NVT_NO, NVT_ROCK, NVT_SNOW, NVT_IRR, NVT_PARK
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_DATA_COVER, ONLY : XDATA_LAI
!
USE MODI_VEGTYPE_TO_PATCH
#ifdef MNH_MEGAN
USE MODD_MEGAN
USE MODI_INIT_MGN2MECH 
#endif
USE MODI_ABOR1_SFX
!
!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_NP_t), INTENT(INOUT) :: NP
TYPE(MEGAN_SURF_FIELDS_t), INTENT(INOUT) :: MSF
TYPE(MEGAN_t), INTENT(INOUT) :: MGN
!
!*       0.1  declaration of arguments
!
REAL, DIMENSION(:), INTENT(IN)     :: PLAT      ! Lat of the grid cell
CHARACTER(LEN=6), DIMENSION(:),INTENT(IN) :: HSV ! name of all scalar variables
REAL, DIMENSION(:,:),INTENT(IN)     :: PMEGAN_FIELDS     ! EF factors
!
!*       0.1 Declaration of local variables
#ifdef MNH_MEGAN
!
INTEGER :: JI, JSV, JP
!
INTEGER:: IP_TRBE, IP_TRBD, IP_TEBE, IP_TEBD, IP_TENE, &
          IP_BOBD, IP_BONE, IP_BOND, IP_SHRB
!
REAL, DIMENSION(SIZE(K%XCLAY,1),IO%NPATCH) :: ZH_TREE
REAL,DIMENSION(SIZE(K%XCLAY,1)) :: ZSILT
REAL,DIMENSION(SIZE(K%XCLAY,1)) :: ZLAI
!
!IF (.NOT.IO%LTR_ML) THEN
!  CALL ABOR1_SFX('INIT_MEGANN: FATAL ERROR PUT LTR_ML = T in NAM_ISBA (PREP_PGD step)')
!END IF
!
ALLOCATE(MGN%XPFT   (N_MGN_PFT,SIZE(K%XCLAY,1)))
ALLOCATE(MGN%XEF    (N_MGN_SPC,SIZE(K%XCLAY,1)))
ALLOCATE(MGN%NSLTYP (SIZE(K%XCLAY,1)))
ALLOCATE(MGN%XBIOFLX(SIZE(K%XCLAY,1)))
MGN%XBIOFLX(:) = 0.
!
! Prepare the mechanism conversion between Megan and MesoNH
MGN%CMECHANISM = "RELACS2" ! scheme default in MesoNH
!
DO JSV=1,SIZE(HSV)
  IF (TRIM(HSV(JSV))=="DIEN") MGN%CMECHANISM = "RACM"
  IF (TRIM(HSV(JSV))=="ALKA") MGN%CMECHANISM = "RELACS"
  IF (TRIM(HSV(JSV))=="ALKA") MGN%CMECHANISM = "RELACS"
  IF (TRIM(HSV(JSV))=="OLEH") MGN%CMECHANISM = "CACM"
  IF (TRIM(HSV(JSV))=="URG7") MGN%CMECHANISM = "RELACS2"
END DO
!
IF (TRIM(MGN%CMECHANISM)=="RACM"    .OR.TRIM(MGN%CMECHANISM)=="RADM2".OR.TRIM(MGN%CMECHANISM)=="SAPRCII" .OR.&
    TRIM(MGN%CMECHANISM)=="SAPRC99" .OR.TRIM(MGN%CMECHANISM)=="CBMZ" .OR.TRIM(MGN%CMECHANISM)=="SAPRC99X".OR.&
    TRIM(MGN%CMECHANISM)=="SAPRC99Q".OR.TRIM(MGN%CMECHANISM)=="CB05" .OR.TRIM(MGN%CMECHANISM)=="CB6"     .OR.&
    TRIM(MGN%CMECHANISM)=="SOAX") THEN
  MGN%CMECHANISM2 = MGN%CMECHANISM
ELSE
  MGN%CMECHANISM2 = "SAPRC99" ! megan default
END IF
!
MGN%LCONVERSION = .TRUE.
!
CALL INIT_MGN2MECH(MGN%CMECHANISM2, MGN%LCONVERSION, MGN%CVNAME3D, MGN%CMECH_SPC, MGN%NSPMH_MAP, &
                   MGN%NMECH_MAP, MGN%XCONV_FAC, MGN%XMECH_MWT, MGN%NVARS3D, MGN%N_SCON_SPC      )
!
DO JSV=1,SIZE(HSV)
  IF (TRIM(HSV(JSV)) == "NO")     MGN%NNO    = JSV  ! ReLACS
  IF (TRIM(HSV(JSV)) == "ALD")    MGN%NALD   = JSV  ! ReLACS
  IF (TRIM(HSV(JSV)) == "BIO")    MGN%NBIO   = JSV  ! ReLACS
  IF (TRIM(HSV(JSV)) == "ALKA")   MGN%NALKA  = JSV  ! ReLACS
  IF (TRIM(HSV(JSV)) == "ALKE")   MGN%NALKE  = JSV  ! ReLACS
  IF (TRIM(HSV(JSV)) == "ARO")    MGN%NARO   = JSV  ! ReLACS
  IF (TRIM(HSV(JSV)) == "CARBO")  MGN%NCARBO = JSV  ! ReLACS
  !
  IF (TRIM(HSV(JSV)) == "ETHE")  MGN%NETHE = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "OLEL")  MGN%NOLEL = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "OLEH")  MGN%NOLEH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ALKL")  MGN%NALKL = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ALKM")  MGN%NALKM = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ALKH")  MGN%NALKH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "AROH")  MGN%NAROH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "AROL")  MGN%NAROL = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "AROO")  MGN%NAROO = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "AROL")  MGN%NAROL = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ARAL")  MGN%NARAL = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ARAC")  MGN%NARAC = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "PAH")   MGN%NPAH  = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ALD2")  MGN%NALD2 = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "KETL")  MGN%NKETL = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "KETH")  MGN%NKETH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "MEOH")  MGN%NMEOH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ETOH")  MGN%NETOH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ALCH")  MGN%NALCH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ISOP")  MGN%NISOP = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "BIOL")  MGN%NBIOL = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "BIOH")  MGN%NBIOH = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "MTBE")  MGN%NMTBE = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "MVK")   MGN%NMVK  = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "MCR")   MGN%NMCR  = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "MGLY")  MGN%NMGLY = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ACID")  MGN%NACID = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ORA1")  MGN%NORA1 = JSV  ! ReLACS2 or CACM
  IF (TRIM(HSV(JSV)) == "ORA2")  MGN%NORA2 = JSV  ! ReLACS2 or CACM
END DO
!
DO JSV=1,SIZE(MGN%CVNAME3D) ! megan species (racm family)
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ISO")  MGN%NISO  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "CH4")  MGN%NCH4  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ETH")  MGN%NETH  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "HC3")  MGN%NHC3  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "HC5")  MGN%NHC5  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "HC8")  MGN%NHC8  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "OL2")  MGN%NOL2  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "OLI")  MGN%NOLI  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "OLT")  MGN%NOLT  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ALD")  MGN%NALD  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "KET")  MGN%NKET  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "TOL")  MGN%NTOL  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "HCHO") MGN%NHCHO = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ORA1") MGN%NORA1 = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ORA2") MGN%NORA2 = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "API")  MGN%NAPI  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "LIM")  MGN%NLIM  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "CO")   MGN%NCO   = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "SO2")  MGN%NSO2  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "NO")   MGN%NNO   = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "HNO3") MGN%NHNO3 = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "NO2")  MGN%NNO2  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "NR")   MGN%NNR   = JSV
END DO
!
DO JSV=1,SIZE(MGN%CVNAME3D) ! megan species (soax family)
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ISP")  MGN%NISP  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "TRP")  MGN%NTRP  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "XYLA") MGN%NXYLA = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "CG5")  MGN%NCG5  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "SQT")  MGN%NSQT  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "TOLA") MGN%NTOLA = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "CG6")  MGN%NCG6  = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "CG4")  MGN%NCG4  = JSV
END DO
!
DO JSV=1,SIZE(MGN%CVNAME3D) !megan species (saprc family)
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ISOPRENE") MGN%NISOPRENE = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "TRP1")     MGN%NTRP1     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ACET")     MGN%NACET     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "MEK")      MGN%NMEK      = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "HCOOH")    MGN%NHCOOH    = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "CCO_OH")   MGN%NCCO_OH   = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "CCHO")     MGN%NCCHO     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "RCHO")     MGN%NRCHO     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "RCO_OH")   MGN%NRCO_OH   = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "BALD")     MGN%NBALD     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ETHENE")   MGN%NETHENE   = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ALK4")     MGN%NALK4     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ALK5")     MGN%NALK5     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ARO1")     MGN%NARO1     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "ARO2")     MGN%NARO2     = JSV
  IF (TRIM(MGN%CVNAME3D(JSV)) == "OLE1")     MGN%NOLE1     = JSV
END DO
!
! Compute soil USDA type 
!
! CLAY : CLAY >= 0.40 SILT < 0.40 SAND < 0.45
! SANDY CLAY : CLAY >= 0.36 SAND >= 0.45 
! SILTY CLAY : CLAY >= 0.40 SILT >= 0.40 
! SILT : SILT >= 0.8 CLAY < 0.12
! SAND : SAND >= 0.3*CLAY + 0.87
! SANDY CLAY LOAM : CLAY >= 0.28  CLAY < 0.36 SAND >= 0.45 | CLAY >= 0.20 CLAY < 0.28 SILT < 0.28
! SILTY CLAY LOAM : CLAY >= 0.28 CLAY < 0.40 SAND < 0.20
! CLAY LOAM : CLAY >= 0.28 CLAY < 0.40 SAND >= 0.20 SAND < 0.45
! SILT LOAM : SILT >= 0.8 CLAY >= 0.12 | SILT >= 0.5 SILT < 0.8 CLAY < 0.28
! LOAMY SAND : SAND >= CLAY + 0.7 SAND < 0.3*CLAY + 0.87
! SANDY LOAM : SAND >= 0.52  CLAY < 0.20 | SAND >= (0.5 - CLAY)  CLAY < 0.07
! LOAM : CLAY >= 0.20 CLAY < 0.28 SILT >= 0.28 SILT < 0.5 | SAND >= (0.5 - CLAY)  CLAY < 0.20
!
ZSILT(:) = 1. - K%XCLAY(:,1) - K%XSAND(:,1)
!
WHERE (ZSILT(:) <= 0.) ZSILT(:) = 0.0
!
DO JI = 1, SIZE(K%XCLAY,1)

  IF ( K%XCLAY(JI,1)>=0.28 ) THEN
    IF ( K%XSAND(JI,1)>=0.45 ) THEN
      IF (K%XCLAY(JI,1)>=0.36 ) THEN     ! Sandy Clay 
        MGN%NSLTYP(JI) = 9
      ELSE                            ! Sandy Clay Loam 
        MGN%NSLTYP(JI) = 6
      ENDIF
    ELSEIF ( K%XCLAY(JI,1)>=0.40 ) THEN
      IF ( ZSILT(JI)>=0.40 ) THEN    ! Silty Clay 
        MGN%NSLTYP(JI) = 10
      ELSE                            ! Clay 
        MGN%NSLTYP(JI) = 11
      ENDIF
    ELSEIF (K%XSAND(JI,1)>=0.20 ) THEN  ! Clay Loam 
      MGN%NSLTYP(JI) = 8
    ELSE                             ! Silty Clay Loam 
      MGN%NSLTYP(JI) = 7
    ENDIF
  ENDIF
  !
  IF ( ZSILT(JI)>=0.8 .AND. K%XCLAY(JI,1)<0.12 ) THEN ! Silt 
    MGN%NSLTYP(JI) = 12
  ELSEIF ( K%XCLAY(JI,1)<0.28 ) THEN    ! ( clay est forcément < 0.28 )
    IF ( ZSILT(JI) >= 0.5 ) THEN      ! Silt Loam 
      MGN%NSLTYP(JI) = 4
    ELSEIF ( K%XCLAY(JI,1)>=0.20 ) THEN
      IF ( ZSILT(JI)>=0.28 ) THEN   ! Loam  
        MGN%NSLTYP(JI) = 5
      ELSE                           ! Sandy Clay Loam 
        MGN%NSLTYP(JI) = 6
      ENDIF
    ENDIF
  ENDIF
  !
  IF ( K%XSAND(JI,1)>=(0.3*K%XCLAY(JI,1) + 0.87) ) THEN   ! Sand 
    MGN%NSLTYP(JI) = 1
  ELSEIF ( K%XSAND(JI,1)>=(K%XCLAY(JI,1) + 0.7) ) THEN    ! Loamy Sand
    MGN%NSLTYP(JI) = 2
  ELSEIF ( K%XSAND(JI,1)>=0.52 .AND. K%XCLAY(JI,1)<0.20 ) THEN ! Sandy Loam            
    MGN%NSLTYP(JI) = 3
  ELSEIF ( K%XSAND(JI,1)>=(0.5 - K%XCLAY(JI,1)) ) THEN
    IF ( K%XCLAY(JI,1)<0.07 ) THEN                   ! Sandy Loam
      MGN%NSLTYP(JI) = 3
    ELSEIF ( K%XCLAY(JI,1)<0.20 ) THEN               ! Loam
      MGN%NSLTYP(JI) = 5
    ENDIF
  ENDIF
  !
ENDDO
!
! Passage des type de végétation isba/vegtype avec ceux de Megan
!
IP_TRBE = VEGTYPE_TO_PATCH(NVT_TRBE, IO%NPATCH)
IP_TRBD = VEGTYPE_TO_PATCH(NVT_TRBD, IO%NPATCH)
IP_TEBE = VEGTYPE_TO_PATCH(NVT_TEBE, IO%NPATCH)
IP_TEBD = VEGTYPE_TO_PATCH(NVT_TEBD, IO%NPATCH)
IP_TENE = VEGTYPE_TO_PATCH(NVT_TENE, IO%NPATCH)
IP_BOBD = VEGTYPE_TO_PATCH(NVT_BOBD, IO%NPATCH)
IP_BONE = VEGTYPE_TO_PATCH(NVT_BONE, IO%NPATCH)
IP_BOND = VEGTYPE_TO_PATCH(NVT_BOND, IO%NPATCH)
IP_SHRB = VEGTYPE_TO_PATCH(NVT_SHRB, IO%NPATCH)
!
MGN%XPFT(:,:) = 0.
!
ZH_TREE(:,:) = XUNDEF
DO JP = 1,IO%NPATCH
  DO JI = 1,NP%AL(JP)%NSIZE_P
    ZH_TREE(NP%AL(JP)%NR_P(JI),JP) = NP%AL(JP)%XH_TREE(JI)
  ENDDO
ENDDO
!
! 1 Needleleaf evergreen temperate tree
! -------------------------------------
! utilisation de la classe NVT_CONI pour 30 < LAT < 60 
WHERE ((PLAT(:) .GE. 30.) .AND. (PLAT(:) .LT. 60.))
  MGN%XPFT(1,:) = S%XVEGTYPE(:,NVT_TENE)
END WHERE
WHERE ((PLAT(:) .LE. -30.) .AND. (PLAT(:) .GT. -60.))
  MGN%XPFT(1,:) = S%XVEGTYPE(:,NVT_TENE)
END WHERE
!
! 2 Needleleaf evergreen boreal tree
! -------------------------------------
!utilisation de la classe NVT_CONI  pour LAT > 60
WHERE ((PLAT(:) .GE. 60.) .OR. (PLAT(:) .LE. -60.))
  MGN%XPFT(2,:) = S%XVEGTYPE(:,NVT_BONE)
END WHERE
!
!3 Needleleaf deciduous boreal tree
! -------------------------------------
!utilisation de la classe NVT_TREE pour LAT > 60
WHERE ((PLAT(:) .GE. 60.) .OR. (PLAT(:) .LE. -60.))
  MGN%XPFT(3,:) = S%XVEGTYPE(:,NVT_BOND)
END WHERE
!
!4 Broadleaf evergreen tropical tree
! -------------------------------------
!utilisation de la classe NVT_EVER pour -30 < LAT < 30
! et une hauteur d'arbre supérieur à 3 m
WHERE (((PLAT(:) .GE. -30.) .AND. (PLAT(:) .LE. 30.)).AND.&
       (ZH_TREE(:,IP_TRBE) .GE. 3.).AND.(ZH_TREE(:,IP_TRBE) .NE. XUNDEF))
MGN%XPFT(4,:) = S%XVEGTYPE(:,NVT_TRBE)
END WHERE
!
!5 Broadleaf evergreen temperate tree
! -------------------------------------
! utilisation de la classe NVT_EVER pour 30 < LAT < 60 
! et une hauteur d'arbre supérieur à 3 m. 
WHERE (((PLAT(:) .GE. 30.) .AND. (PLAT(:) .LT. 60.)).AND.&
       (ZH_TREE(:,IP_TEBE) .GE. 3.).AND.(ZH_TREE(:,IP_TEBE) .NE. XUNDEF))
MGN%XPFT(5,:) = S%XVEGTYPE(:,NVT_TEBE)
END WHERE
WHERE (((PLAT(:) .LE. -30.) .AND. (PLAT(:) .GT. -60.)).AND.&
       (ZH_TREE(:,IP_TEBE) .GE. 3.).AND.(ZH_TREE(:,IP_TEBE) .NE. XUNDEF))
MGN%XPFT(5,:) = S%XVEGTYPE(:,NVT_TEBE)
END WHERE
!
!6 Broadleaf deciduous tropical tree
! -------------------------------------
!utilisation de la classe NVT_TREE pour -30 < LAT < 30
WHERE ((PLAT(:) .GE. -30.) .AND. (PLAT(:) .LE. 30.))
MGN%XPFT(6,:) = S%XVEGTYPE(:,NVT_TRBD)
END WHERE
!
!7 Broadleaf deciduous temperate tree
! -------------------------------------
!utilisation de la classe NVT_TREE pour 30 < LAT < 60
! en utilisant une hauteur d'arbre supérieur à 3 m
WHERE (((PLAT(:) .GE. 30.) .AND. (PLAT(:) .LT. 60.)).AND.&
       (ZH_TREE(:,IP_TEBD) .GE. 3.).AND.(ZH_TREE(:,IP_TEBD) .NE. XUNDEF))
MGN%XPFT(7,:) = S%XVEGTYPE(:,NVT_TEBD)
END WHERE
WHERE (((PLAT(:) .LE. -30.) .AND. (PLAT(:) .GT. -60.)).AND.&
       (ZH_TREE(:,IP_TEBD) .GE. 3.).AND.(ZH_TREE(:,IP_TEBD) .NE. XUNDEF))
MGN%XPFT(7,:) = S%XVEGTYPE(:,NVT_TEBD)
END WHERE
!
!8 Broadleaf deciduous boreal tree
! -------------------------------------
!utilisation de la classe NVT_TREE pour LAT > 60
WHERE (((PLAT(:) .GE. 60.) .OR. (PLAT(:) .LE. -60.)).AND.&
       (ZH_TREE(:,IP_BOBD) .GE. 3.).AND.(ZH_TREE(:,IP_BOBD) .NE. XUNDEF))
MGN%XPFT(8,:) = S%XVEGTYPE(:,NVT_BOBD)
END WHERE
!
!9 Broadleaf evergreen shrub
! -------------------------------------
!utilisation de la classe NVT_EVER pour une hauteur d'arbre inférieure à 3 m
WHERE (ZH_TREE(:,IP_SHRB) .LT. 3.)
MGN%XPFT(9,:) = S%XVEGTYPE(:,NVT_SHRB)
END WHERE
!
!10 Broadleaf deciduous temperate shrub
! -------------------------------------
!utilisation de la classe NVT_TREE pour une hauteur d'arbre inférieure à 3 m
! et  pour 30 < LAT < 60
WHERE ((ZH_TREE(:,IP_SHRB) .LT. 3.) .AND. (ZH_TREE(:,IP_SHRB).NE. XUNDEF) .AND. &
       ((PLAT(:) .GE. 30.) .AND. (PLAT(:) .LT. 60.)))
MGN%XPFT(10,:) = S%XVEGTYPE(:,NVT_SHRB)
END WHERE
WHERE ((ZH_TREE(:,IP_SHRB) .LT. 3.) .AND. (ZH_TREE(:,IP_SHRB).NE. XUNDEF) .AND. &
       ((PLAT(:) .LE. -30.) .AND. (PLAT(:) .GT. -60.)))
MGN%XPFT(10,:) = S%XVEGTYPE(:,NVT_SHRB)
END WHERE
!
!11 Broadleaf deciduous boreal_shrub
! -------------------------------------
!utilisation de la classe NVT_TREE pour une hauteur d'arbre inférieure à 3 m
! et  pour  LAT > 60
WHERE ((ZH_TREE(:,IP_SHRB) .LT. 3.) .AND. (ZH_TREE(:,IP_SHRB).NE. XUNDEF) .AND. &
       ((PLAT(:) .GE. 60.) .OR. (PLAT(:) .LE. -60.)))
MGN%XPFT(11,:) = S%XVEGTYPE(:,NVT_SHRB)
END WHERE
!
!12 C3 arctic grass
! -------------------------------------
!utilisation de la classe NVT_GRAS + NVT_PARK pour LAT > 60
WHERE ((PLAT(:) .GE. 60.) .OR. (PLAT(:) .LE. -60.))
MGN%XPFT(12,:) = S%XVEGTYPE(:,NVT_GRAS) + S%XVEGTYPE(:,NVT_PARK)
ELSEWHERE
!
!13 C3 non-arctic grass
! -------------------------------------
!utilisation de la classe NVT_GRAS + NVT_PARK ailleur
MGN%XPFT(13,:) = S%XVEGTYPE(:,NVT_GRAS) +  S%XVEGTYPE(:,NVT_PARK)
END WHERE
!
!14 C4 grass
! -------------------------------------
! utilisation de la classe NVT_TROG 
MGN%XPFT(14,:) = S%XVEGTYPE(:,NVT_TROG) 
!
!15 Corn
! -------------------------------------
! utilisation de la classe NVT_C4
MGN%XPFT(15,:) = S%XVEGTYPE(:,NVT_C4) 
!
!16 Wheat
! -------------------------------------
! utilisation de la classe NVT_C3
MGN%XPFT(16,:) = S%XVEGTYPE(:,NVT_C3)
!
! Emission factor
MGN%XEF(:,:) = 0.
!
! Default values
! 1: ISOP isoprene
MGN%XEF(1,:) = 6000.
! 2: MYRC myrcene
MGN%XEF(2,:) = 20.
! 3: SABI sabinene
MGN%XEF(3,:) = 300.
! 4: LIMO limonene 
MGN%XEF(4,:) = 80.
! 5: A_3CAR carene_3
MGN%XEF(5,:) = 40.
! 6: OCIM ocimene_t_b
MGN%XEF(6,:) = 40.
! 7: BPIN pinene_b
MGN%XEF(7,:) = 125.
! 8: APIN pinene_a
MGN%XEF(8,:) = 300.
! 9: OMTP A_2met_styrene + cymene_p + cymene_o + phellandrene_a + thujene_a + terpinene_a 
!          + terpinene_g + terpinolene + phellandrene_b + camphene + bornene + fenchene_a
!          + ocimene_al + ....
! 10: FARN
! 11: BCAR
! 12: OSQT
! 13: MBO
! 14: MEOH
! 15: ACTO
! 16: CO
! 17: NO
!MGN%XEF(17,:,1) = 30.
! 18: BIDER
! 19: STRESS
! 20: OTHER
! Values from the megan maps fields (to be introduced at the PREP_PGD step - nameliste PRE_PGD1.nam)
DO JSV=1, MSF%NMEGAN_NBR
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFISOP")   MGN%XEF(1,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFMYRC")   MGN%XEF(2,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFSABI")   MGN%XEF(3,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFLIMO")   MGN%XEF(4,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFCARE")   MGN%XEF(5,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFOCIM")   MGN%XEF(6,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFBPIN")   MGN%XEF(7,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFAPIN")   MGN%XEF(8,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFOMTP")   MGN%XEF(9,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFFARN")   MGN%XEF(10,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFBCAR")   MGN%XEF(11,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFOSQT")   MGN%XEF(12,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFMBO")    MGN%XEF(13,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFMEOH")   MGN%XEF(14,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFACTO")   MGN%XEF(15,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFCO")     MGN%XEF(16,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFNO")     MGN%XEF(17,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFBIDER")  MGN%XEF(18,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFSTRESS") MGN%XEF(19,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "EFOTHER")  MGN%XEF(20,:) = PMEGAN_FIELDS(:,JSV)
!  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "LAI")      PLAI(:,1)     = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT1")     MGN%XPFT(1,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT2")     MGN%XPFT(2,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT3")     MGN%XPFT(3,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT4")     MGN%XPFT(4,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT5")     MGN%XPFT(5,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT6")     MGN%XPFT(6,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT7")     MGN%XPFT(7,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT8")     MGN%XPFT(8,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT9")     MGN%XPFT(9,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT10")    MGN%XPFT(10,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT11")    MGN%XPFT(11,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT12")    MGN%XPFT(12,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT13")    MGN%XPFT(13,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT14")    MGN%XPFT(14,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT15")    MGN%XPFT(15,:) = PMEGAN_FIELDS(:,JSV)
  IF (TRIM(MSF%CMEGAN_NAME(JSV)) == "PFT16")    MGN%XPFT(16,:) = PMEGAN_FIELDS(:,JSV)
END DO

#endif
!
!---------------------------------------------------------------------------
!
END SUBROUTINE INIT_MEGAN_n

