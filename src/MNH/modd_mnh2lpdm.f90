!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
        MODULE MODD_MNH2LPDM
!	##################
!-----------------------------------------------------------------------------
!*	MODD_MNH2S2 : BORDEL MIMIESQUE, NO COMMENT.
!
!	Auteur      : Michel Bouzom, DP/SERV/ENV
!	Creation    : 18.01.2002
!-----------------------------------------------------------------------------
!
!
!*	FICHIERS ET NAMELISTS.
!	----------------------
!
INTEGER,PARAMETER :: JPMNHMAX = 30   ! Nombre max de  fichiers FM.
INTEGER,SAVE      :: NBMNH           ! Nombre reel de fichiers FM.
INTEGER,SAVE      :: IGRILLE         ! numero de la grille Meso-nh
!
CHARACTER(LEN=50), SAVE, DIMENSION(JPMNHMAX) :: CFMTO  ! Nom des fichiers METEO.
CHARACTER(LEN=50), SAVE                      :: CFGRI  ! Nom du fichier GRILLE.
CHARACTER(LEN=50), SAVE                      :: CFDAT  ! Nom du fichier DATE.
CHARACTER(LEN=28), SAVE, DIMENSION(JPMNHMAX) :: CFMNH  ! Noms des FM.
!
CHARACTER*10, SAVE :: CTURBPARAM="ISOTROPE"   ! ISOTROPE ou HANNA
!
NAMELIST/NAM_FIC/ CFMTO,CFGRI,CFDAT,CFMNH
NAMELIST/NAM_TURB/ CTURBPARAM
!
!
!
!*	MESO-NH.
!	--------
!
!*	Dimensions et indices utiles Meso-NH.
!
INTEGER, SAVE :: NIU,NJU,NKU             ! Taille des tableaux de travail.
INTEGER, SAVE :: NIB,NIE,NJB,NJE,NKB,NKE ! Bornes du domaine physiques.
!
!*	Champs Meso-NH a extraire.
!
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XUT,XVT,XWT  ! Vent.
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XTHT         ! Tempe potentielle.
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XTKET  ! m2/s2
REAL, DIMENSION(:,:),  ALLOCATABLE :: XINRT ! Taux de precipitation en mm/h
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XLM  ! m
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XDISSIP  ! J/Kg
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XRMVT ! Rapport de melange en vapeur d eau
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XRMCT ! Rapport de melange en eau liquide
                                            ! nuageuse
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XRMRT ! Rapport de melange en eau liquide
                                            ! pluie
REAL, DIMENSION(:,:),ALLOCATABLE :: XSFU, XSFV ! flux cinematiques w'u' et w'v'
REAL, DIMENSION(:,:,:),ALLOCATABLE :: XWPTHP ! flux de chaleur sensible 
                                              ! w'theta'
REAL, DIMENSION(:,:),ALLOCATABLE :: XZ0
!
!*	VARIABLES POUR FICHIER METEO.
!	-----------------------------
!
INTEGER, SAVE :: NSIMAX,NSJMAX          ! Dimensions du domaine a extraire.
INTEGER, SAVE :: NSIB,NSIE,NSJB,NSJE    ! Bornes du domaine a extraire.
!
REAL,    SAVE :: XPASXM,   XPASYM       ! Pas selon X et Y en metres.
REAL,    SAVE :: ZMAILLE                ! Max de XPASXM XPASYM.
REAL,    SAVE :: XPTSOLON, XPTSOLAT     ! Lon et Lat du point sud-ouest.
REAL,    SAVE :: XPTX0LON, XPTY0LAT     ! Lon et Lat du point x=0, y=0
INTEGER, SAVE :: NMDLAA,NMDLMM,NMDLJJ   ! Date du modele.
INTEGER, SAVE :: NMDLHH,NMDLMN,NMDLSS   ! Heure du modele.
!
REAL, SAVE, DIMENSION(:),     ALLOCATABLE :: XSHAUT ! niveaux hauteurs LPDM  [m]
REAL, SAVE, DIMENSION(:,:),   ALLOCATABLE :: XSREL    ! Relief   [m]
REAL, SAVE, DIMENSION(:,:),   ALLOCATABLE :: XSCORIOZ  ! Coef. coriolis
REAL, SAVE, DIMENSION(:,:),   ALLOCATABLE :: XSZ0     ! Longueur de rugosite [m]
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSU      
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSV     
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSW   
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSTH   
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSTKE
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSLM
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSDISSIP
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: XSWPTHP
REAL, DIMENSION(:,:),   ALLOCATABLE :: XSSFU
REAL, DIMENSION(:,:),   ALLOCATABLE :: XSSFV
REAL, DIMENSION(:,:),   ALLOCATABLE :: XSINRT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSRMV
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSRMC
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSRMR
REAL, DIMENSION(:,:),   ALLOCATABLE :: XSHMIX ! hauteur de melange en m
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSSIGU ! variance du vent composante x
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSSIGW ! variance du vent composante y
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSTIMEU ! temps lagrangien Tlx
REAL, DIMENSION(:,:,:), ALLOCATABLE :: XSTIMEW ! temps lagrangien Tly
REAL, SAVE , DIMENSION(:,:,:), ALLOCATABLE :: XSPHI   ! Hauteur geopotentielle
REAL, DIMENSION(:,:),   ALLOCATABLE :: XSUSTAR ! Vitesse de friction
REAL, DIMENSION(:,:),   ALLOCATABLE :: XSWSTAR ! Vitesse verticale convective
REAL, DIMENSION(:,:),   ALLOCATABLE :: XSLMO   ! Longueur de Monin Obukov
REAL, DIMENSION(:), ALLOCATABLE :: XSTHETAV ! Tempe virtuelle potentielle
!
REAL               :: XRVSRD                 ! Rv/Rd.
REAL               :: XTHSOL
!
END MODULE MODD_MNH2LPDM
