!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_INI_RADIATIONS
!     ##########################
!
INTERFACE
!
    SUBROUTINE INI_RADIATIONS(TPINIFILE,HLUOUT,OINIRAD,TPDTCUR,TPDTEXP,&
         PZZ,PDXX,PDYY,                                                &
         PSINDEL,PCOSDEL,PTSIDER,PCORSOL,PSLOPANG,PSLOPAZI,            &
         PDTHRAD,PDIRFLASWD,PSCAFLASWD,                                &
         PFLALWD,PDIRSRFSWD,KCLEARCOL_TM1,                             &
         PZENITH, PAZIM, TPDTRAD_FULL,TPDTRAD_CLONLY,TPINITHALO2D_ll,  &
         PRADEFF,PSWU,PSWD,PLWU,PLWD,PDTHRADSW,PDTHRADLW               )
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_IO_ll,       ONLY : TFILEDATA
USE MODD_TYPE_DATE
!
TYPE(TFILEDATA),        INTENT(IN)  :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
LOGICAL,                INTENT(IN)  :: OINIRAD   ! switch to initialize or read
                                                 ! the radiation informations
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date and time
                                                 ! Ajout PP
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ        ! height z
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX,PDYY  ! metric coefficients
REAL,                   INTENT(OUT):: PSINDEL    ! sine   of the solar declination angle
REAL,                   INTENT(OUT):: PCOSDEL    ! cosine of the solar declination angle
REAL,                   INTENT(OUT):: PTSIDER    ! sideral decimal time correction
REAL,                   INTENT(OUT):: PCORSOL    ! daily solar constant correction
REAL, DIMENSION(:,:),   INTENT(OUT):: PSLOPANG   ! slope angle
REAL, DIMENSION(:,:),   INTENT(OUT):: PSLOPAZI   ! azimuthal slope angle
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDTHRAD    ! radiative tendency
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDIRFLASWD ! Direct    Short Wave flux at the flat surface
REAL, DIMENSION(:,:,:), INTENT(OUT):: PSCAFLASWD ! Scaterred Short Wave flux at the flat surface
REAL, DIMENSION(:,:),   INTENT(OUT):: PFLALWD    ! Long Wave flux at the surface
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDIRSRFSWD ! Direct    Short Wave flux at the surface
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCLEARCOL_TM1 ! trace of cloud col at the previous
                                                    ! radiation step
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PZENITH        ! Solar zenithal angle
REAL, DIMENSION(:,:), INTENT(INOUT):: PAZIM          ! Solar azimuthal angle
TYPE (DATE_TIME),       INTENT(OUT):: TPDTRAD_FULL   ! date and time of the 
                                                     ! last full radiation call
TYPE (DATE_TIME),       INTENT(OUT):: TPDTRAD_CLONLY ! date and time of the 
                         ! last radiation call only for the cloudy verticals
TYPE(LIST_ll), POINTER             :: TPINITHALO2D_ll! pointer for the list of fields
                                                     !  which must be communicated in INIT
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PSWU ! upward SW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PSWD ! downward SW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PLWU ! upward LW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PLWD ! downward LW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PDTHRADSW !  dthrad sw
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PDTHRADLW !  dthrad lw
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PRADEFF ! effective radius
!
END SUBROUTINE INI_RADIATIONS
!
END INTERFACE
!
END MODULE MODI_INI_RADIATIONS
!
!
!   ####################################################################
    SUBROUTINE INI_RADIATIONS(TPINIFILE,HLUOUT,OINIRAD,TPDTCUR,TPDTEXP,&
         PZZ,PDXX,PDYY,                                                &
         PSINDEL,PCOSDEL,PTSIDER,PCORSOL,PSLOPANG,PSLOPAZI,            &
         PDTHRAD,PDIRFLASWD,PSCAFLASWD,                                &
         PFLALWD,PDIRSRFSWD,KCLEARCOL_TM1,                             &
         PZENITH, PAZIM, TPDTRAD_FULL,TPDTRAD_CLONLY,TPINITHALO2D_ll,  &
         PRADEFF,PSWU,PSWD,PLWU,PLWD,PDTHRADSW,PDTHRADLW               )
!   ####################################################################
!
!!****  *INI_RADIATION_TIME * - initialisation for radiation scheme in the MesoNH framework
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!! 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	V. Masson        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2003  split of ini_radiations. externalization of ISBA
!!      O.Thouron   06/2008  Add diagnostics
!!      P. Peyrille, M. Tomasini 06/2012:  if LFIX_DAT=T, TDTCUR is replaced by
!!                               TDTEXP to  have a perpetual day ie. the diurnal cycle is retained 
!!                               but the day stays the same during the whole run 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!MESO-NH modules
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_CST,         ONLY : XPI
USE MODD_CONF,        ONLY : LFLAT, L2D
USE MODD_IO_ll,       ONLY : TFILEDATA
USE MODD_LES
USE MODD_PARAMETERS,  ONLY : JPVEXT
USE MODD_PARAM_RAD_n, ONLY : LFIX_DAT
USE MODD_TYPE_DATE
!
USE MODE_FMREAD
USE MODE_ll
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(TFILEDATA),        INTENT(IN)  :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN)  :: HLUOUT    ! name for output-listing
                                                 !  of nested models
LOGICAL,                INTENT(IN)  :: OINIRAD   ! switch to initialize or read
                                                 ! the radiation informations
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR    ! Current date and time
TYPE (DATE_TIME),       INTENT(IN) :: TPDTEXP    ! Current date and time
                                                 ! Ajout PP
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ        ! height z
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX,PDYY  ! metric coefficients
REAL,                   INTENT(OUT):: PSINDEL    ! sine   of the solar declination angle
REAL,                   INTENT(OUT):: PCOSDEL    ! cosine of the solar declination angle
REAL,                   INTENT(OUT):: PTSIDER    ! sideral decimal time correction
REAL,                   INTENT(OUT):: PCORSOL    ! daily solar constant correction
REAL, DIMENSION(:,:),   INTENT(OUT):: PSLOPANG   ! slope angle
REAL, DIMENSION(:,:),   INTENT(OUT):: PSLOPAZI   ! azimuthal slope angle
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDTHRAD    ! radiative tendency
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDIRFLASWD ! Direct    Short Wave flux at the flat surface
REAL, DIMENSION(:,:,:), INTENT(OUT):: PSCAFLASWD ! Scaterred Short Wave flux at the flat surface
REAL, DIMENSION(:,:),   INTENT(OUT):: PFLALWD    ! Long Wave flux at the surface
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDIRSRFSWD ! Direct    Short Wave flux at the surface
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCLEARCOL_TM1 ! trace of cloud col at the previous
                                                    ! radiation step
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PZENITH        ! Solar zenithal angle
REAL, DIMENSION(:,:), INTENT(INOUT):: PAZIM          ! Solar azimuthal angle
TYPE (DATE_TIME),       INTENT(OUT):: TPDTRAD_FULL   ! date and time of the 
                                                     ! last full radiation call
TYPE (DATE_TIME),       INTENT(OUT):: TPDTRAD_CLONLY ! date and time of the 
TYPE(LIST_ll), POINTER             :: TPINITHALO2D_ll! pointer for the list of fields
                                                     !  which must be communicated in INIT
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PSWU ! upward SW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PSWD ! downward SW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PLWU ! upward LW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PLWD ! downward LW Flux 
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PDTHRADSW !  dthrad sw
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PDTHRADLW !  dthrad lw
REAL, DIMENSION(:,:,:),     INTENT(OUT) :: PRADEFF ! effective radius
!
!*       0.2   declarations of local variables
!
INTEGER, DIMENSION(0:11) :: IBIS, INOBIS ! Cumulative number of days per month
                                         ! for bissextile and regular years
REAL :: ZDATE         ! Julian day of the year
REAL :: ZAD           ! Angular Julian day of the year
REAL :: ZDECSOL       ! Daily solar declination angle 
REAL :: ZA1, ZA2      ! Ancillary variables
!
INTEGER :: JI         ! loop counter
INTEGER :: IIB        ! I index value of the first inner mass point
INTEGER :: IJB        ! J index value of the first inner mass point
INTEGER :: IKB        ! K index value of the first inner mass point
INTEGER :: IIE        ! I index value of the last inner mass point
INTEGER :: IJE        ! J index value of the last inner mass point

REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),1) :: ZSLOPX, ZSLOPY ! Terrain slopes in
                                                             ! x and y directions
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
!
!*       1.    COMPUTES THE DAY OF THE YEAR
!              ----------------------------
!
INOBIS(:) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
IBIS(0) = INOBIS(0)
DO JI=1,11
  IBIS(JI) = INOBIS(JI)+1
END DO
IF ( LFIX_DAT ) THEN   ! Ajout PP 
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
!*       5.     COMPUTES THE SLOPE ANGLE AND THE AZIMUTHAL SLOPE ANGLE
!	        ------------------------------------------------------
!
!
IF(LFLAT) THEN
  PSLOPANG(:,:) = 0.0
  PSLOPAZI(:,:) = -0.5*XPI 
ELSE 
  !  . local computation
  ZSLOPX(:,:,:) = MXF( DXM(PZZ(:,:,IKB:IKB))/PDXX(:,:,IKB:IKB) )
  !
  IF(L2D) THEN
    PSLOPANG(:,:) = ATAN(ABS(ZSLOPX(:,:,1)))
    PSLOPAZI(:,:) = -0.5*XPI
  ELSE
    !    . local computation
    ZSLOPY(:,:,:) = MYF( DYM(PZZ(:,:,IKB:IKB))/PDYY(:,:,IKB:IKB) ) 
    PSLOPANG(:,:) = ATAN(SQRT(ZSLOPX(:,:,1)**2+ZSLOPY(:,:,1)**2))
    PSLOPAZI(:,:) = 0.0
    PSLOPAZI(:,:) = - 0.5*XPI +                       &
         ATAN2( ZSLOPY(:,:,1), ZSLOPX(:,:,1) + SIGN(1.E-30,ZSLOPX(:,:,1)) )
  END IF
END IF
!
!   5.2 Update halo of PSLOPANG and PSLOPAZI at the end of ini_modeln
!
CALL ADD2DFIELD_ll (TPINITHALO2D_ll,PSLOPANG)
CALL ADD2DFIELD_ll (TPINITHALO2D_ll,PSLOPAZI)
!
!-------------------------------------------------------------------------------
!
!*        9.    INITIALIZE TIME FOR THE RADIATION CALL
!	            --------------------------------------
!
PSWU(:,:,:)   = 0.
PSWD(:,:,:)   = 0.
PLWU(:,:,:)   = 0.
PLWD(:,:,:)   = 0.
PDTHRADSW(:,:,:)   = 0.
PDTHRADLW(:,:,:)   = 0.
PRADEFF(:,:,:)   = 0.
!
IF ( OINIRAD ) THEN
  IF (LFIX_DAT ) THEN                      ! Ajout PP
     TPDTRAD_FULL   = TPDTEXP             ! Ajout PP
     TPDTRAD_CLONLY = TPDTEXP             ! Ajout PP
  ELSE                                    ! Ajout PP
     TPDTRAD_FULL     = TPDTCUR
     TPDTRAD_CLONLY   = TPDTCUR
  END IF
  PDTHRAD(:,:,:)   = 0.
  PDIRFLASWD(:,:,:)= 0.
  PSCAFLASWD(:,:,:)= 0.
  PFLALWD(:,:)     = 0.
  PDIRSRFSWD(:,:,:)= 0.
  KCLEARCOL_TM1    = 0
ELSE
  CALL IO_READ_FIELD(TPINIFILE,'DTRAD_FULL',  TPDTRAD_FULL)
  CALL IO_READ_FIELD(TPINIFILE,'DTRAD_CLLY',  TPDTRAD_CLONLY)
  CALL IO_READ_FIELD(TPINIFILE,'DTHRAD',      PDTHRAD)
  CALL IO_READ_FIELD(TPINIFILE,'FLALWD',      PFLALWD)
  CALL IO_READ_FIELD(TPINIFILE,'DIRFLASWD',   PDIRFLASWD)
  CALL IO_READ_FIELD(TPINIFILE,'SCAFLASWD',   PSCAFLASWD)
  CALL IO_READ_FIELD(TPINIFILE,'DIRSRFSWD',   PDIRSRFSWD)
  CALL IO_READ_FIELD(TPINIFILE,'CLEARCOL_TM1',KCLEARCOL_TM1)
  CALL IO_READ_FIELD(TPINIFILE,'ZENITH',      PZENITH)
  CALL IO_READ_FIELD(TPINIFILE,'AZIM',        PAZIM)
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_RADIATIONS
