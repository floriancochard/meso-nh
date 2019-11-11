!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_INI_RADIATIONS_ECMWF
!     ##########################
!
INTERFACE
!
    SUBROUTINE INI_RADIATIONS_ECMWF(HINIFILE,HLUOUT,                      &
         PZHAT,PPABST,PTHT,PTSRAD,PLAT,PLON,TPDTCUR,TPDTEXP,              &
         HLW,KDLON,KFLEV,KFLUX,KRAD,KSWB,HAER,KAER,KSTATM,                &
         PSTATM,PSEA,PTOWN,PBARE,POZON, PAER,PDST_WL, OSUBG_COND                    )
!
USE MODD_TYPE_DATE
!
CHARACTER (LEN=*),      INTENT(IN)  :: HINIFILE  ! Name of the initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
CHARACTER (LEN=*),      INTENT(IN) :: HAER       ! aerosol optical thickness climatology
CHARACTER (LEN=4),      INTENT(IN) :: HLW        ! LW scheme used
!
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! height level without orography
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPABST! pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT  !Temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PTSRAD ! surface radiative temperature
REAL, DIMENSION (:,:),  INTENT(IN) :: PSEA   ! sea fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PTOWN  ! town fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PBARE  ! bare soil fraction
REAL, DIMENSION(:,:),   INTENT(IN) :: PLAT, PLON ! arrays of latitude-longitude
!
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
                                                 !  which must be communicated in INIT
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date&time. Ajout PP
!
INTEGER, INTENT(OUT) :: KDLON   ! number of columns where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLEV   ! number of vertical levels where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLUX   ! number of top and ground fluxes in the output
INTEGER, INTENT(OUT) :: KRAD    ! number of satellite radiances to synthesize
INTEGER, INTENT(OUT) :: KAER    ! number of AERosol classes
INTEGER, INTENT(OUT) :: KSWB    ! number of SW band
INTEGER, INTENT(OUT) :: KSTATM  ! index of the STAndard ATMosphere level just
                                ! above the model top
!
REAL, DIMENSION(:,:),   INTENT(OUT) :: PSTATM   ! working standard atmosphere
!
REAL, DIMENSION(:,:,:),   POINTER :: POZON      ! ozone mixing ratio ( from climato.)
REAL, DIMENSION(:,:,:,:), POINTER :: PAER       ! aerosols optical thickness (from climato)
REAL, DIMENSION(:,:,:,:), POINTER :: PDST_WL    ! aerosols optical thickness (from climato)
LOGICAL, INTENT(IN)  ::   OSUBG_COND ! Switch for sub-grid condensation
!
END SUBROUTINE INI_RADIATIONS_ECMWF
!
END INTERFACE
!
END MODULE MODI_INI_RADIATIONS_ECMWF
!
!
!   #######################################################################
    SUBROUTINE INI_RADIATIONS_ECMWF(HINIFILE,HLUOUT,                      &
         PZHAT,PPABST,PTHT,PTSRAD,PLAT,PLON,TPDTCUR,TPDTEXP,              &
         HLW,KDLON,KFLEV,KFLUX,KRAD,KSWB,HAER,KAER,KSTATM,                &
         PSTATM,PSEA,PTOWN,PBARE,POZON, PAER, PDST_WL,OSUBG_COND                    )
!   #######################################################################
!
!!****  *INI_RADIATIONS * - initialisation for ECMWF radiation scheme in the MesoNH framework
!!
!!    PURPOSE
!!    -------
!!      Firstly, the purpose of this routine is to compute the solar declination angle,
!!    the daily sideral hour and solar constant corrections and the slope 
!!    angles. Then the closest standard atmosphere is selected by checking 
!!    the mean latitude and the seasonal extrema. It will allow to complement the
!!    atmosphere above the top of the domain for radiation calculation.
!!    Secondly the initialisation of ECMWF radiation package is performed,  
!!     as well as the set up of radiatively active component.
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
!!      Subroutine INI_RADCONF   : Initialization of ECMWF radiation package constant
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
!!      Book2 of documentation ( routine INI_RADIATIONS )
!!      Paltridge, G.W. and Platt, C.M.R. (1976) in "Radiative Processes
!!        in Meteorology and Climatology", Elsevier, New-York
!!
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/02/95
!!      Revised     12/09/95
!!       (J.Stein)  02/10/95 add the RESTA configuration
!!       (J.Stein)  01/04/96 bug correction in the slope computations + rm the
!!                           sunset and sunrise computation + bug in ZAD 
!!      (J.Stein and V. Masson)  27/08/96  bug correction in the slope
!!                                         computation
!!      (J.-P.Pinty)  21/11/96 add ABS to the ZLATMEAN test in the sub-tropics
!!      (V.Masson)  18/08/97 call to fmread directly with dates and strings
!!      (V.Masson)  07/11/97 test on pressure levels
!!      (N. Asencio)  11/08/98 add parallel code
!!      (P. Jabouille) 25/05/99 replace 'DTRAD_CLONLY' by 'DTRAD_CLLY' (size too long)
!!
!!      (F. Solmon 03/03/02) : MAJOR MODIFICATIONS, updtating of MNH radiation scheme from
!!                            ECMWF last version 
!!      (V. Masson) 01/2004  split of the routine in two (externalization 
!!                           of the surface)
!!      (A.Grini) 07/2005 add dust
!!      (M.Tomasini P.Peyrille)  06/2012  to set date to a perpetual day if LFIX_DAT=T
!!      (V. Masson)          replaces cover fractions by sea/town/bare soil fractions
!!      J.Escobar 30/03/2017  : Management of compilation of ECMWF_RAD in REAL*8 with MNH_REAL=R4
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!ECMWF radiation scheme specific modules 
!
USE PARKIND1,  ONLY : JPRB
USE YOEAERD  , ONLY : RCAEOPS  ,RCAEOPL  ,RCAEOPU  ,RCAEOPD  ,RCTRBGA  ,&
     RCVOBGA  ,RCSTBGA  ,RCTRPT   ,RCAEADM  ,RCAEROS  ,&
     RCAEADK
USE OYOMCST , ONLY : RTT
USE YOETHF , ONLY : RTICE
!
!MESO-NH modules
!
USE MODD_TYPE_DATE
USE MODD_CONF
USE MODD_CST
USE MODD_GRID
USE MODD_PARAMETERS
USE MODD_STAND_ATM
USE MODD_PARAM_RAD_n,  ONLY: LFIX_DAT
!
USE MODE_ll
USE MODE_FM
!
USE MODI_INI_RADCONF
USE MODI_INI_HOR_AERCLIM
!
USE MODE_DUSTOPT
USE MODE_SALTOPT
USE MODE_CONSRAD
USE MODE_REPRO_SUM
USE MODI_INI_STAND_ATM
USE MODI_SUECOZC
USE MODI_RADOZC
USE MODI_SUAERV
USE MODI_RADAER
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=*),      INTENT(IN)  :: HINIFILE  ! Name of the initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
CHARACTER (LEN=*),      INTENT(IN) :: HAER       ! aerosol optical thickness climatology
CHARACTER (LEN=4),      INTENT(IN) :: HLW        ! LW scheme used
!
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! height level without orography
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPABST! pressure
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTHT  !Temperature
REAL, DIMENSION(:,:),   INTENT(IN) :: PTSRAD ! surface radiative temperature
REAL, DIMENSION (:,:),  INTENT(IN) :: PSEA   ! sea fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PTOWN  ! town fraction
REAL, DIMENSION (:,:),  INTENT(IN) :: PBARE  ! bare soil fraction
REAL, DIMENSION(:,:),   INTENT(IN) :: PLAT, PLON ! arrays of latitude-longitude
!
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date&time. Ajout PP
!
INTEGER, INTENT(OUT) :: KDLON   ! number of columns where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLEV   ! number of vertical levels where the radiation
                                ! calculations are performed
INTEGER, INTENT(OUT) :: KFLUX   ! number of top and ground fluxes in the output
INTEGER, INTENT(OUT) :: KRAD    ! number of satellite radiances to synthesize
INTEGER, INTENT(OUT) :: KAER    ! number of AERosol classes
INTEGER, INTENT(OUT) :: KSWB    ! number of SW band
INTEGER, INTENT(OUT) :: KSTATM  ! index of the STAndard ATMosphere level just
                                ! above the model top
!
REAL, DIMENSION(:,:),   INTENT(OUT) :: PSTATM   ! working standard atmosphere
!
REAL, DIMENSION(:,:,:),   POINTER :: POZON      ! ozone mixing ratio ( from climato.)
REAL, DIMENSION(:,:,:,:), POINTER :: PAER       ! aerosols optical thickness (from climato)
REAL, DIMENSION(:,:,:,:), POINTER :: PDST_WL    ! aerosols optical thickness (from climato)
LOGICAL, INTENT(IN)  ::   OSUBG_COND ! Switch for sub-grid condensation
!
!
!*       0.2   declarations of local variables
!
LOGICAL :: GSUMMER ! .T. when SUMMERtime
LOGICAL :: GWINTER ! .T. when WINTERtime
LOGICAL :: GSEASON ! .T. when SUMMERtime in the northern hemisphere or
                   !     when WINTERtime in the southern hemisphere
!
INTEGER :: JI, JJ, JK, JK1, JKRAD,IIJ,JL  ! loop index
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
REAL :: ZLATMEAN      ! MEAN LATitude in the domain
REAL :: ZLAT_TROPICAL ! TROPIQUE LATitude
REAL :: ZLAT_POLAR    ! POLAR circle LATitude
!
REAL, DIMENSION(:,:),ALLOCATABLE :: ZLON          ! longitude
REAL, DIMENSION(SIZE(PSTATM,1)) :: ZZSTAT ! half level altitudes of standard atm.
!
INTEGER :: IINFO_ll                   ! return code of parallel routine
INTEGER :: IIMAX_ll,IJMAX_ll          ! Number of points of
                                      ! Global physical domain
                                      ! in the x and y directions
!
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZEXNT ! Exner function
!
! Variables for aerosols and ozone climatologies set up
REAL, DIMENSION (:),     ALLOCATABLE  :: ZAESEA, ZAELAN, ZAEURB, ZAEDES
REAL(KIND=JPRB), DIMENSION (:), ALLOCATABLE  :: ZAESEA_RAD, ZAELAN_RAD, ZAEURB_RAD, ZAEDES_RAD
LOGICAL, DIMENSION (:,:),ALLOCATABLE  :: GAFRICA, GASIA, GAUSTRALIA
REAL, DIMENSION (:,:),   ALLOCATABLE  :: ZDESERT ! desert fraction
REAL, DIMENSION (:,:),   ALLOCATABLE  :: ZPAVE, ZWORK_GRID
REAL(KIND=JPRB), DIMENSION (:,:,:), ALLOCATABLE  :: ZAER
REAL(KIND=JPRB), DIMENSION (:,:),   ALLOCATABLE  :: ZPRES_HL,ZT_HL,ZOZON
REAL(KIND=JPRB), DIMENSION (:,:),   ALLOCATABLE  :: ZCVDAES, ZCVDAEL, ZCVDAEU, ZCVDAED,ZETAH
REAL(KIND=JPRB), DIMENSION (:),     ALLOCATABLE  :: ZGEMU 
REAL, DIMENSION(:),      ALLOCATABLE  :: ZAECOV_SEA, ZAECOV_URB, ZAECOV_LAN, ZAECOV_DES
INTEGER :: ZYMD, ZHOURS   ! date for climatology initialisation
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       0.1  INITIALIZATIONS
!
!*       0.2  COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
!
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
IKU = SIZE(PTHT,3)
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
! size of global physical domain
!
CALL GET_GLOBALDIMS_ll ( IIMAX_ll,IJMAX_ll )
!
!-------------------------------------------------------------------------------
!
!*       6.     INITIALIZES THE STANDARD ATMOSPHERES AND CHOOSES THE CLOSEST ONE
!	        ----------------------------------------------------------------
!
CALL INI_STAND_ATM
!
!  latitude dependency
!
IF(.NOT.LCARTESIAN) THEN
  ! . global sum
  ZLATMEAN = SUM_DD_R2_ll( PLAT(IIB:IIE,IJB:IJE))
  ZLATMEAN = ZLATMEAN /FLOAT(IIMAX_ll*IJMAX_ll)
ELSE 
  ZLATMEAN = XLAT0
ENDIF
!
!  summer/winter dependency
!
IF (LFIX_DAT ) THEN           ! Ajout PP 
   IF( (TPDTEXP%TDATE%MONTH.GE.4).AND.(TPDTEXP%TDATE%MONTH.LE.9) ) THEN
     GSUMMER = .TRUE.
     GWINTER = .FALSE.
   ELSE
     GSUMMER = .FALSE.
     GWINTER = .TRUE.
   END IF
ELSE
   IF( (TPDTCUR%TDATE%MONTH.GE.4).AND.(TPDTCUR%TDATE%MONTH.LE.9) ) THEN
     GSUMMER = .TRUE.
     GWINTER = .FALSE.
   ELSE
     GSUMMER = .FALSE.
     GWINTER = .TRUE.
   END IF
END IF
!
GSEASON = (ZLATMEAN.GE.0.0.AND.GSUMMER).OR.(ZLATMEAN.LE.0.0.AND.GWINTER)
ZLAT_TROPICAL = 23.5
ZLAT_POLAR    = 66.5
!
IF( ABS(ZLATMEAN).LE.ZLAT_TROPICAL ) THEN
  PSTATM(:,:) = XSTROATM(:,:)
ELSE IF( ABS(ZLATMEAN).LE.ZLAT_POLAR ) THEN
  IF( GSEASON ) THEN
    PSTATM(:,:) = XSMLSATM(:,:)
  ELSE
    PSTATM(:,:) = XSMLWATM(:,:)
  END IF
ELSE IF( GSEASON ) THEN
  PSTATM(:,:) = XSPOSATM(:,:)
ELSE
  PSTATM(:,:) = XSPOWATM(:,:)
END IF
!
!
!  search for KSTATM the level number in the standard atmosphere which
!                    is both immediately above the highest model level
!                    and the lowest model pressure level
!
! securities:
! ----------
!
!       z (ike+1)        < 0.5 * z(kstatm-1) + 0.5 * z(kstatm)
! and   pressure (ike+1) > pressure (kstatm)
!
! The first condition insures a correct succession in temperature profiles
! The second one insures a monotoneous pressure profile, with a security
! margin of one half model level.
!
ZZSTAT(1)  = 0.
ZZSTAT(2:) = 0.5* ( PSTATM(:SIZE(PSTATM,1)-1,1) + PSTATM(2:,1) )
!
KSTATM=MAX( 1 + COUNT(        PZHAT     (IKE+1)  >= 1000.0*ZZSTAT(:)   ) ,  &
            1 + COUNT( MIN_ll(PPABST,IINFO_ll,1,1,IKE+1,                    &
                              IIMAX_ll+2*JPHEXT,IJMAX_ll+2*JPHEXT,IKE+1)    &
                                                 <=  100.0*PSTATM(:,2) ) )
!
KSTATM=MAX(2,KSTATM)
!
! if KSTATM is greater than the dimension of the standard atmosphere profile
! array upper boundary, then the MESO-NH profile will be used alone
! in the radiation code.
!
!  computes KFLEV the number of vertical levels used in 
!                 the radiative calculations:
!  - in the IKE-JPVEXT lowest levels the model thermodynamical structure is used
!  - in the upper SIZE(PSTATM,1)-KSTATM levels the thermal structure is given
!    by the standard atmosphere
!
KFLEV = IKE -JPVEXT + SIZE(PSTATM,1) - KSTATM + 1
!
!-------------------------------------------------------------------------------
!
!*       7.     INITIALIZES THE ECMWF RADIATION PACKAGE 
!	        ------------------------------------------------
!
KDLON = (IIE-IIB+1)*(IJE-IJB+1) ! number of column
KFLUX = 6
KRAD  = 0
KAER  = 6
KSWB  = 6 ! number of SW band 
!
!       7.1    Initialization of the ECMWF physical constants,
!              and grid-independant quantities
!       
CALL INI_RADCONF (HLW,KSWB,OSUBG_COND)
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
ZT_HL(:,:) = 273.
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
WHERE (ZT_HL(:,:) == 0.)
  ZT_HL(:,:) = 273.
ENDWHERE
!
!  Standard atmosphere extension
!* begining at ikup+1 level allows to use a model domain higher than 50km
!
IKUP   = IKE-JPVEXT+1
!
DO JK=IKUP+1,KFLEV+1
  JK1 = (KSTATM-1)+(JK-IKUP)
  ZPRES_HL(:,JK) = PSTATM(JK1,2)*100.0
  ZT_HL(:,JK) = PSTATM(JK1,3)
END DO
!
! vertical grid inversion for compatibility with ECMWF routines
!
ALLOCATE (ZWORK_GRID(SIZE(ZPRES_HL,1),KFLEV+1))
!
!half level pressure
!
ZWORK_GRID(:,:)=ZPRES_HL(:,:)
DO JKRAD=1, KFLEV+1
  JK1=(KFLEV+1)+1-JKRAD
  ZPRES_HL(:,JKRAD) = ZWORK_GRID(:,JK1)
END DO
!
!half level temperature
!
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
!
ALLOCATE (ZOZON(KDLON,KFLEV))
!
IF (LFIX_DAT ) THEN     ! Ajout PP 
 ZYMD = TPDTEXP%TDATE%YEAR * 1E4 + TPDTEXP%TDATE%MONTH * 1E2 + TPDTEXP%TDATE%DAY
 ZHOURS = INT(TPDTEXP%TIME / 60.)
ELSE
 ZYMD = TPDTCUR%TDATE%YEAR * 1E4 + TPDTCUR%TDATE%MONTH * 1E2 + TPDTCUR%TDATE%DAY
 ZHOURS = INT(TPDTCUR%TIME / 60.)
END IF
!
! Fortuin langematz climatology loading
!
CALL SUECOZC ( ZYMD , ZHOURS )
!
! Interpolation on the simulation domain
!
CALL RADOZC ( 1 , KDLON, KDLON , 1, KFLEV, 1 ,&
     KDLON ,0,ZPRES_HL, ZGEMU,       &
     ZOZON                         )
!
!*    8.3         AEROSOLS climatogy 
!
ALLOCATE (ZAER(KDLON, KFLEV,KAER))
!
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
  IF ( HAER =='TANR' .OR. HAER == 'TEGE' ) THEN
    CALL INI_HOR_AERCLIM (HAER,IIB,IIE,IJB,IJE,KDLON,ZYMD,ZHOURS, &
                          PLAT,PLON,ZAESEA,ZAELAN,ZAEURB,ZAEDES )
  END IF
!
! AEROSOLS from SURFACE FRACTIONS
!
  IF( HAER =='SURF') THEN
!
    !* deserts are only considered over Africa, southern Asia, Australia
    !* longitude between -180 and 180 for geographical tests
    !  Only bare soil fractions larger than 0.5 are supposed to contribute to
    !  desert aerosols
    !
    ALLOCATE(ZDESERT   (IIU,IJU))
    ZDESERT(:,:) = 0.

    IF (.NOT.LCARTESIAN) THEN
      !* deserts are only considered over Africa, southern Asia, Australia
      ALLOCATE(ZLON      (IIU,IJU))
      ALLOCATE(GAFRICA   (IIU,IJU))
      ALLOCATE(GASIA     (IIU,IJU))
      ALLOCATE(GAUSTRALIA(IIU,IJU))    
      !* longitude between -180 and 180 for geographical tests
      ZLON = PLON(:,:) - NINT(PLON/360.)*360.
      GAFRICA   (:,:) = PLAT(:,:) > -36.086389  .AND. PLAT(:,:) < 36.010556 &
                  .AND. ZLON(:,:) > -73.18      .AND. ZLON(:,:) < 34.158611
      GASIA     (:,:) = PLAT(:,:) > 4.358056    .AND. PLAT(:,:) < 55.335278 &
                  .AND. ZLON(:,:) > -123.157778 .AND. ZLON(:,:) <-34.285556
      GAUSTRALIA(:,:) = PLAT(:,:) > -39.561389  .AND. PLAT(:,:) < -10.251667 &
                  .AND. ZLON(:,:) > -155.041944 .AND. ZLON(:,:) < -111.405556
      !
      !  Only bare soil fractions larger than 0.5 are supposed to contribute to
      !  desert aerosols
      !
      WHERE (GAFRICA(:,:) .OR. GASIA(:,:) .OR. GAUSTRALIA(:,:)) &
      ZDESERT(:,:) = MAX( 2.*(PBARE(:,:)-0.5) , 0.)
      !
      !
    ELSE
      !
      ZDESERT(:,:) = MAX( 2.*(PBARE(:,:)-0.5) , 0.)
      !
    ENDIF
    !
    !* fills sea, town, desert and land surface covers for aerosols distributions
    DO JJ=IJB,IJE
      DO JI=IIB,IIE
        IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
        ZAESEA(IIJ) = PSEA(JI,JJ)
        ZAEURB(IIJ) = PTOWN(JI,JJ)
        ZAEDES(IIJ) = ZDESERT(JI,JJ)
        ZAELAN(IIJ) = MAX( 1.- ZAESEA(IIJ) - ZAEURB(IIJ) - ZAEDES(IIJ) , 0.)
      END DO
    END DO
    IF (ALLOCATED(ZLON)) DEALLOCATE(ZLON)
    IF (ALLOCATED(GAFRICA))     DEALLOCATE(GAFRICA)
    IF (ALLOCATED(GASIA))     DEALLOCATE(GASIA)
    IF (ALLOCATED(GAUSTRALIA))     DEALLOCATE(GAUSTRALIA)
    IF (ALLOCATED(ZDESERT))     DEALLOCATE(ZDESERT)

  END IF
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
  IF ( HAER =='SURF') THEN
     RCAEOPS = 0.001 ! Sea instead 0.02 to agree with TEGEN
     RCAEOPL = 0.05  ! Land (continental)
     RCAEOPU = 0.3   ! Urban zone
     RCAEOPD = 0.5   ! Desert
  END IF
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
!* 8.4              Adaptation on mnh domain
!
ALLOCATE (POZON(IIU,IJU,KFLEV))
ALLOCATE (PAER(IIU,IJU,KFLEV,KAER))
ALLOCATE (PDST_WL(IIU,IJU,KFLEV,KSWB))
!
POZON (:,:,:)=0.
PAER (:,:,:,:)=0.
PDST_WL (:,:,:,:)=0.
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    IIJ = 1 + (JI-IIB) + (IIE-IIB+1)*(JJ-IJB)
    DO JKRAD=1,KFLEV      
      JK1 = KFLEV +1 -JKRAD
      POZON (JI,JJ,JKRAD) = ZOZON (IIJ,JK1)
      PAER  (JI,JJ,JKRAD,:) = ZAER (IIJ,JK1,:)
    END DO
  END DO
END DO
!
!        Read in look up tables of dust optical properties
!No arguments, all look up tables are defined in module
!modd_dust_opt_lkt
CALL DUST_OPT_LKT_SET1()
CALL DUST_OPT_LKT_SET2()
CALL DUST_OPT_LKT_SET3()
CALL DUST_OPT_LKT_SET4()
CALL DUST_OPT_LKT_SET5()
CALL DUST_OPT_LKT_SET6()
CALL DUST_OPT_LKT_SET7()
CALL DUST_OPT_LKT_SET8()
CALL DUST_OPT_LKT_SET9()
CALL DUST_OPT_LKT_SET10()

CALL SALT_OPT_LKT_SET1()
CALL SALT_OPT_LKT_SET2()
CALL SALT_OPT_LKT_SET3()
CALL SALT_OPT_LKT_SET4()
CALL SALT_OPT_LKT_SET5()
CALL SALT_OPT_LKT_SET6()
CALL SALT_OPT_LKT_SET7()
CALL SALT_OPT_LKT_SET8()
CALL SALT_OPT_LKT_SET9()
CALL SALT_OPT_LKT_SET10()


!
CALL INI_CONS_PROP_OP
DEALLOCATE (ZPRES_HL) 
DEALLOCATE (ZPAVE)
DEALLOCATE (ZT_HL)
DEALLOCATE (ZETAH)
DEALLOCATE (ZGEMU)
DEALLOCATE (ZOZON)
DEALLOCATE (ZAER)
!
!-------------------------------------------------------------------------------
!
RTICE=RTT-23.
!
END SUBROUTINE INI_RADIATIONS_ECMWF
