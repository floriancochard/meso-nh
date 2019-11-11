!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_AEROZON
!     ##########################
!
INTERFACE
    SUBROUTINE AEROZON(PPABST,PTHT,PTSRAD,PLAT,PLON,TPDTCUR,TPDTEXP,  &
         KDLON,KFLEV,HAER,KAER,KSTATM,                             &
         PSINDEL,PCOSDEL,PTSIDER,PCORSOL,                          &
         PSTATM,POZON, PAER)
!
USE MODD_TIME
!
CHARACTER (LEN=*),      INTENT(IN) :: HAER       ! aerosol optical thickness climatology
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPABST! pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT  !Temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PTSRAD ! surface radiative temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PLAT, PLON ! arrays of latitude-longitude
!
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date and time

!
INTEGER, INTENT(IN) :: KDLON   ! number of columns where the radiation
                                ! calculations are performed
INTEGER, INTENT(IN) :: KFLEV   ! number of vertical levels where the radiation
                                ! calculations are performed
INTEGER, INTENT(IN) :: KAER    ! number of AERosol classes
INTEGER, INTENT(IN) :: KSTATM  ! index of the STAndard ATMosphere level just
                                ! above the model top
REAL, INTENT(OUT)    :: PSINDEL ! sine   of the solar declination angle
REAL, INTENT(OUT)    :: PCOSDEL ! cosine of the solar declination angle
REAL, INTENT(OUT)    :: PTSIDER ! sideral decimal time correction
REAL, INTENT(OUT)    :: PCORSOL ! daily solar constant correction
!
REAL, DIMENSION(:,:),   INTENT(IN) :: PSTATM   ! working standard atmosphere
!
REAL, DIMENSION(:,:,:),   POINTER :: POZON      ! ozone mixing ratio ( from climato.)
REAL, DIMENSION(:,:,:,:), POINTER :: PAER       ! aerosols optical thickness (from climato)
!
END SUBROUTINE AEROZON
!
END INTERFACE
!
END MODULE MODI_AEROZON
!
!  #################################################################
    SUBROUTINE AEROZON(PPABST,PTHT,PTSRAD,PLAT,PLON,TPDTCUR,TPDTEXP,   &
         KDLON,KFLEV,HAER,KAER,KSTATM,                             &
         PSINDEL,PCOSDEL,PTSIDER,PCORSOL,                          &
         PSTATM,POZON, PAER)
!   ###############################################################
!
!!****  *AEROZON * - initialisation for ECMWF radiation scheme in the MesoNH framework
!!
!!    PURPOSE
!!    -------
!!      Basically the same purpose as ini_radiations but only for aerosols and
!!      ozon distribution. 
!!
!!**  METHOD
!!    ------
!!      The formulae of astronomy are taken from Paltridge and Platt (1976).
!!    The latitudinal and seasonal dependences for the choice of the standard
!!    atmosphere is straightforward.
!!     The set up of ECMWF radiation package is mainly done in the "ini_radconf" routine
!!    in function of defined options ( eg use of RRTM for LW).   
!!     Ozone field and Aerosols optical thickness are defined through specific 
!!    climatologies. For aerosols distribution set up, several option are allowed : 
!!      1) in function of underlying cover types ( land, sea, desert, town) 
!!      2) using the ECMWF global distribution of Tanre (1984) interpolated on the domain  
!!   The carbon dioxide concentration is homogeneously distributed at 360 ppm. Note that 
!!   other active gaz are directly initialised in ECMWF Package
!!   Finally, the radiative tendency and the short and long waves surface
!!    fluxes are read or set to 0 and also the instants at which full or 
!!    partial radiations' call has been performed. 
!!
!!    EXTERNAL
!!    --------
!!      
!!      Subroutine INI_STAND_ATM : Initialization of 5 standard atmosphere     
!!      Subroutine AEROZON   : Initialization of ECMWF radiation package constant
!!   CAUTION: the following routine do not use the MNH norm   
!!      Subroutine  SUECOZC         : ozone climatology loading
!!      Subroutine  RADOZC          : ozone climatology interpolation
!!      Subroutine  INI_HOR_AERCLIM : aerosol climatology : horizontal distribution on domain
!!      Subroutine  SUAERV          : aerosol climatology : vertical distribution
!!      Subroutine  RADAER          : interpolation on domain 
!!     these routine use ECMWF specific module "yo===="
!!   END CAUTION  
!!      GET_DIM_EXT_ll : get extended sub-domain sizes
!!      GET_INDICE_ll  : get physical sub-domain bounds
!!      GET_GLOBALDIMS_ll : get physical global domain sizes
!!      REDUCESUM_ll   : sum into a scalar variable
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!! 
!!     Use of "yo===" ECMWF radiation package specific modules for climatologies 
!!         initialisation           
!!     Module MODD_CONF : define model configuration
!!         LCARTESIAN : knowledge of PLAT and PLON
!!         L2D        : bidimensionnal case
!!         LFLAT      : flat terrain
!!         CCONF      : temporal configuration START or RESTART
!!      Module MODD_CST : define constants
!!         XPI : pi
!!         XMD : molar mass of dry air
!!         XDAY: length of the day
!!      Module MODD_GRID : define the projection parameters
!!         XLAT0 : reference latitude
!!         XLON0 : reference longitude
!!      Module MODD_PARAMETERS : parameter variables
!!         JPVEXT   : Vertical EXTernal point number
!!      Module MODD_STAND_ATM : contains 5 standard atmospheres
!!         XSTROATM : standard tropical atmosphere
!!         XSMLSATM : standard mid-latitudes summer atmosphere
!!         XSMLWATM : standard mid-latitudes winter atmosphere
!!         XSPOSATM : standard polar summer atmosphere
!!         XSPOWATM : standard polar winter atmosphere
!!
!!    REFERENCE
!!    ---------
!!      ECMWF IFS radiation documentation
!!      Book2 of documentation ( routine AEROZON )
!!      Paltridge, G.W. and Platt, C.M.R. (1976) in "Radiative Processes
!!        in Meteorology and Climatology", Elsevier, New-York
!!
!!    AUTHOR
!!    ------
!!  	P. Peyrille
!!
!!    MODIFICATIONS
!!    -------------
!!      (P.Peyrille) 20/07/04 : add LFIX_DAT to have perpetual day
!!      J.Escobar 30/03/2017  : Management of compilation of ECMWF_RAD in REAL*8 with MNH_REAL=R4
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!ECMWF radiation scheme specific modules 
!
USE PARKIND1 , ONLY : JPRB
USE YOEAERD  , ONLY : RCAEOPS  ,RCAEOPL  ,RCAEOPU  ,RCAEOPD  ,RCTRBGA  ,&
     RCVOBGA  ,RCSTBGA  ,RCTRPT   ,RCAEADM  ,RCAEROS  ,&
     RCAEADK
!
!MESO-NH modules
!
USE MODD_CONF
USE MODD_CST
USE MODD_GRID
USE MODD_PARAMETERS
USE MODD_STAND_ATM
USE MODD_TIME
USE MODD_GROUND_PAR
USE MODD_PARAM_RAD_n,  ONLY: LFIX_DAT
!
USE MODE_ll
USE MODE_FM
USE MODE_FMREAD
!
USE MODI_SHUMAN
USE MODI_INI_RADCONF
USE MODI_INI_HOR_AERCLIM
USE MODI_SUECOZC
USE MODI_RADOZC
USE MODI_SUAERV
USE MODI_RADAER
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=*),      INTENT(IN) :: HAER       ! aerosol optical thickness climatology
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPABST! pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT  !Temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PTSRAD ! surface radiative temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PLAT, PLON ! arrays of latitude-longitude
!
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date and time

!
INTEGER, INTENT(IN) :: KDLON   ! number of columns where the radiation
                                ! calculations are performed
INTEGER, INTENT(IN) :: KFLEV   ! number of vertical levels where the radiation
                                ! calculations are performed
INTEGER, INTENT(IN) :: KAER    ! number of AERosol classes
INTEGER, INTENT(IN) :: KSTATM  ! index of the STAndard ATMosphere level just
                                ! above the model top
REAL, INTENT(OUT)    :: PSINDEL ! sine   of the solar declination angle
REAL, INTENT(OUT)    :: PCOSDEL ! cosine of the solar declination angle
REAL, INTENT(OUT)    :: PTSIDER ! sideral decimal time correction
REAL, INTENT(OUT)    :: PCORSOL ! daily solar constant correction
!
REAL, DIMENSION(:,:),   INTENT(IN) :: PSTATM   ! working standard atmosphere
!
REAL, DIMENSION(:,:,:),   POINTER :: POZON      ! ozone mixing ratio ( from climato.)
REAL, DIMENSION(:,:,:,:), POINTER :: PAER       ! aerosols optical thickness (from climato)
!
                         ! last radiation call only for the cloudy verticals
!
!
INTEGER :: JI, JJ, JK, JK1, JKRAD,IIJ,JL 
!
INTEGER :: IIB           ! I index value of the first inner mass point
INTEGER :: IJB           ! J index value of the first inner mass point
INTEGER :: IKB           ! K index value of the first inner mass point
INTEGER :: IIE           ! I index value of the last inner mass point
INTEGER :: IJE           ! J index value of the last inner mass point
INTEGER :: IKE           ! K index value of the last inner mass point
INTEGER :: IIU           ! array size for the first  index
INTEGER :: IJU           ! array size for the second index
INTEGER :: IKU           ! array size for the third  index
INTEGER :: IKUP          ! 
!
INTEGER, DIMENSION(0:11) :: IBIS, INOBIS ! Cumulative number of days per month
                                         ! for bissextile and regular years
REAL :: ZDATE         ! Julian day of the year
REAL :: ZAD           ! Angular Julian day of the year
REAL :: ZDECSOL       ! Daily solar declination angle 
REAL :: ZA1, ZA2      ! Ancillary variables
!
!
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZEXNT ! Exner function
!
! Variables for aerosols and ozone climatologies set up

REAL, DIMENSION (:,:),   ALLOCATABLE  :: ZPAVE, ZWORK_GRID
REAL, DIMENSION (:),   ALLOCATABLE  :: ZAESEA, ZAELAN, ZAEURB, ZAEDES
!
REAL(KIND=JPRB), DIMENSION (:,:),   ALLOCATABLE  :: ZT_HL
REAL(KIND=JPRB), DIMENSION (:),   ALLOCATABLE  :: ZAESEA_RAD, ZAELAN_RAD, ZAEURB_RAD, ZAEDES_RAD
REAL(KIND=JPRB), DIMENSION (:,:), ALLOCATABLE  :: ZCVDAES, ZCVDAEL, ZCVDAEU, ZCVDAED
REAL(KIND=JPRB), DIMENSION (:,:), ALLOCATABLE  :: ZPRES_HL,ZOZON,ZETAH
REAL(KIND=JPRB), DIMENSION (:),   ALLOCATABLE  :: ZGEMU 
INTEGER :: ZYMD, ZHOURS   ! date for climatology initialisation
INTEGER :: JKCEP,JK_NH
REAL(KIND=JPRB), DIMENSION (:,:,:), ALLOCATABLE  :: ZAER
REAL, DIMENSION(:),      ALLOCATABLE  :: ZAECOV_SEA, ZAECOV_URB, ZAECOV_LAN, ZAECOV_DES
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       0.1  INITIALIZATIONS
!
!
!*       0.2  COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
!
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
IKU = SIZE(PPABST,3)
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
! size of global physical domain
!
!*       1.    COMPUTES THE DAY OF THE YEAR
!              ----------------------------
!
INOBIS(:) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
IBIS(0) = INOBIS(0)
DO JI=1,11
  IBIS(JI) = INOBIS(JI)+1
END DO
!
IF ( LFIX_DAT ) THEN 
  IF( MOD(TPDTEXP%TDATE%YEAR,4).EQ.0 ) THEN
    ZDATE = FLOAT(TPDTEXP%TDATE%DAY +   IBIS(TPDTEXP%TDATE%MONTH-1)) - 1
    ZAD = 2.0*XPI*ZDATE/366.0
  ELSE
    ZDATE = FLOAT(TPDTEXP%TDATE%DAY + INOBIS(TPDTEXP%TDATE%MONTH-1)) - 1
    ZAD = 2.0*XPI*ZDATE/365.0
  END IF
ELSE
  IF( MOD(TPDTCUR%TDATE%YEAR,4).EQ.0 ) THEN
    ZDATE = FLOAT(TPDTCUR%TDATE%DAY +   IBIS(TPDTCUR%TDATE%MONTH-1)) - 1
    ZAD = 2.0*XPI*ZDATE/366.0
  ELSE
    ZDATE = FLOAT(TPDTCUR%TDATE%DAY + INOBIS(TPDTCUR%TDATE%MONTH-1)) - 1
    ZAD = 2.0*XPI*ZDATE/365.0
  END IF
END IF 
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTES THE SOLAR DECLINATION ANGLE
!	        ------------------------------------
!
ZDECSOL = 0.006918-0.399912*COS(ZAD)   +0.070257*SIN(ZAD)    &
         -0.006758*COS(2.*ZAD)+0.000907*SIN(2.*ZAD) &
         -0.002697*COS(3.*ZAD)+0.00148 *SIN(3.*ZAD)
PSINDEL = SIN(ZDECSOL)
PCOSDEL = COS(ZDECSOL)
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTES THE SIDERAL HOUR CORRECTION
!	        ------------------------------------
!
ZA1 = (1.00554*ZDATE- 6.28306)*(XPI/180.0)
ZA2 = (1.93946*ZDATE+23.35089)*(XPI/180.0)
PTSIDER = (7.67825*SIN(ZA1)+10.09176*SIN(ZA2)) / 60.0
!
!-------------------------------------------------------------------------------
!
!*       4.     COMPUTES THE DAILY SOLAR CONSTANT CORRECTION
!	        --------------------------------------------
!
PCORSOL = 1.00011+0.034221*COS(ZAD)   +0.001280*SIN(ZAD)    &
                 +0.000719*COS(2.*ZAD)+0.000077*SIN(2.*ZAD)
!
!-------------------------------------------------------------------------------
!
!*       8.     INITIALIZE RADIATIVELY ACTIVE COMPOUNDS (3D FIELDS) 
!	        ------------------------------------------------------ 
!
!*       8.1   set up for grid dependant quantitites (from initial state) 
! 
ALLOCATE (ZPRES_HL(KDLON,KFLEV+1))
ALLOCATE (ZPAVE(KDLON,KFLEV))
ALLOCATE (ZETAH(KDLON,KFLEV+1))
ALLOCATE (ZT_HL(KDLON,KFLEV+1))
ALLOCATE (ZGEMU(KDLON))
!
ZEXNT(:,:,:)= ( PPABST(:,:,:)/XP00 ) ** (XRD/XCPD)
!
DO JK=IKB,IKE+1
  JKRAD = JK-JPVEXT
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZPRES_HL(IIJ,JKRAD) = XP00 * (0.5*(ZEXNT(JI,JJ,JK)+ZEXNT(JI,JJ,JK-1)))**(XCPD/XRD)
      ZT_HL(IIJ,JKRAD) = 0.5*(PTHT(JI,JJ,JK)*ZEXNT(JI,JJ,JK)+PTHT(JI,JJ,JK-1)*ZEXNT(JI,JJ,JK-1))
    END DO
  END DO
END DO
!
!  Surface temperature at the first level
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    ZT_HL(IIJ,1) = PTSRAD(JI,JJ)
  END DO
END DO
!
!  Standard atmosphere extension
!* begining at ikup+1 level allows to use a model domain higher than 50km
IKUP   = IKE-JPVEXT+1
!
DO JK=IKUP+1,KFLEV+1
  JK1 = (KSTATM-1)+(JK-IKUP)
  ZPRES_HL(:,JK) = PSTATM(JK1,2)*100.0
  ZT_HL(:,JK) = PSTATM(JK1,3)
END DO
!
! vertical grid inversion for compatibility with ECMWF routines
ALLOCATE (ZWORK_GRID(SIZE(ZPRES_HL,1),KFLEV+1))
!
!half level pressure
ZWORK_GRID(:,:)=ZPRES_HL(:,:)
DO JKRAD=1, KFLEV+1
  JK1=(KFLEV+1)+1-JKRAD
  ZPRES_HL(:,JKRAD) = ZWORK_GRID(:,JK1)
END DO
!
!half level temperature
ZWORK_GRID(:,:)=ZT_HL(:,:)
DO  JKRAD=1, KFLEV+1
  JK1=(KFLEV+1)+1-JKRAD
  ZT_HL(:,JKRAD)=ZWORK_GRID(:,JK1)
END DO
DEALLOCATE(ZWORK_GRID)
!
DO JKRAD=1,KFLEV
  ZPAVE(:,JKRAD)=0.5*(ZPRES_HL(:,JKRAD)+ZPRES_HL(:,JKRAD+1))
END DO
!
!coo geographique 
!
IF(.NOT.LCARTESIAN) THEN
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
      ZGEMU(IIJ) = SIN ( PLAT(JI,JJ) * XPI/180)  
    END DO
  END DO
ELSE
  ZGEMU(:) = SIN ( XLAT0 * XPI/180) 
END IF
!
!*     8.2     OZONE climatology 
!     -----------------------------
!
ALLOCATE (ZOZON(KDLON,KFLEV))
!
IF (LFIX_DAT ) THEN 
  ZYMD = TPDTEXP%TDATE%YEAR * 1E4 + TPDTEXP%TDATE%MONTH * 1E2 + TPDTEXP%TDATE%DAY
  ZHOURS = INT(TPDTCUR%TIME / 60.)
ELSE
  ZYMD = TPDTCUR%TDATE%YEAR * 1E4 + TPDTCUR%TDATE%MONTH * 1E2 + TPDTCUR%TDATE%DAY
  ZHOURS = INT(TPDTCUR%TIME / 60.)
END IF 
!
! Fortuin langematz climatology loading
CALL SUECOZC ( ZYMD , ZHOURS )
!
! DOESN'T WORK WITH A ROTATED OR STRETCHED GRID
! Interpolation on the simulation domain
CALL RADOZC ( 1 , KDLON, KDLON , 1, KFLEV, 1 ,&
     KDLON ,0,ZPRES_HL, ZGEMU,       &
     ZOZON                         )
! 
!     8.3 time interpolation of AEROSOLS climatogy 
!     ----------------------------------------------
!
ALLOCATE (ZAER(KDLON, KFLEV,KAER))
IF(HAER /= 'NONE') THEN
!
!     8.3.1 horizontal ditributions   
!
  ALLOCATE (ZAESEA(KDLON))
  ALLOCATE (ZAELAN(KDLON))
  ALLOCATE (ZAEURB(KDLON))
  ALLOCATE (ZAEDES(KDLON))

  ALLOCATE (ZAESEA_RAD(KDLON)) 
  ALLOCATE (ZAELAN_RAD(KDLON))
  ALLOCATE (ZAEURB_RAD(KDLON))
  ALLOCATE (ZAEDES_RAD(KDLON))
!
! AEROSOLS ECMWF climatologies
!
  IF ( HAER == 'TEGE' ) THEN
    CALL INI_HOR_AERCLIM (HAER,IIB,IIE,IJB,IJE,KDLON,ZYMD,ZHOURS, &
                          PLAT,PLON,ZAESEA,ZAELAN,ZAEURB,ZAEDES )
  END IF
!
!
!     8.3.2 vertical ditributions (standard profiles derived from Tanre)
!
  ALLOCATE (ZCVDAES(KDLON,KFLEV+1))
  ALLOCATE (ZCVDAEL(KDLON,KFLEV+1))
  ALLOCATE (ZCVDAEU(KDLON,KFLEV+1))
  ALLOCATE (ZCVDAED(KDLON,KFLEV+1))
  DO JL=1,KDLON
    ZETAH(JL,:)=ZPRES_HL(JL,:)/101300. ! reference pressure for normalisation
  END DO
  WHERE (ZETAH (:,:) > 1.)
    ZETAH(:,:)=1.
  END WHERE
!
! set up of vertical ditribution parameters
  CALL SUAERV ( KDLON, KFLEV   , ZETAH, &
        ZCVDAES ,ZCVDAEL ,ZCVDAEU ,ZCVDAED, &
        RCTRBGA,RCVOBGA,RCSTBGA,RCAEOPS,RCAEOPL,RCAEOPU,&
        RCAEOPD,RCTRPT ,RCAEADK,RCAEADM,RCAEROS )
!  
! modification of initial ECMWF maximum optical thickness 
! for aerosols classes in case of HAER=SURF
! note : these variables belongs to yoeaerd module   
!
! final aerosol profiles on mnh grid
!
  ZAESEA_RAD = ZAESEA ; ZAELAN_RAD = ZAELAN ; ZAEURB_RAD = ZAEURB ; ZAEDES_RAD = ZAEDES
  CALL RADAER (1, KDLON, KDLON, 1, KFLEV, ZPRES_HL,ZT_HL, &
       ZCVDAES ,ZCVDAEL ,ZCVDAEU ,ZCVDAED, &
       ZAESEA_RAD, ZAELAN_RAD, ZAEURB_RAD, ZAEDES_RAD, &
       ZAER )
!
!!- VOLCANIC AEROSOL SET TO epsilon IN ABSENCE OF ERUPTION 
  ZAER(:,:,5) = 1.E-12
!
  DEALLOCATE (ZCVDAES)
  DEALLOCATE (ZCVDAEL)
  DEALLOCATE (ZCVDAEU)
  DEALLOCATE (ZCVDAED)
  DEALLOCATE (ZAESEA)
  DEALLOCATE (ZAELAN)
  DEALLOCATE (ZAEURB)
  DEALLOCATE (ZAEDES)

  DEALLOCATE (ZAESEA_RAD)
  DEALLOCATE (ZAELAN_RAD)
  DEALLOCATE (ZAEURB_RAD)
  DEALLOCATE (ZAEDES_RAD)
ELSE
  ZAER(:,:,:)= 1E-12
END IF
!
!*       8.4   Adaptation on mnh domain
!       --------------------------------
!
!
POZON (:,:,:)=0.
PAER (:,:,:,:)=0.
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    DO JKRAD=1,KFLEV      
      JK1 = KFLEV +1 -JKRAD
      PAER  (JI,JJ,JKRAD,:) = ZAER (IIJ,JK1,:)
      POZON (JI,JJ,JKRAD) = ZOZON (IIJ,JK1)
    END DO
  END DO
END DO
!
!
!
DEALLOCATE (ZPRES_HL) 
DEALLOCATE (ZPAVE)
DEALLOCATE (ZT_HL)
DEALLOCATE (ZETAH)
DEALLOCATE (ZGEMU)
DEALLOCATE (ZOZON)
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AEROZON









