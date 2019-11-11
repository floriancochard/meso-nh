!MNH_LIC Copyright 2002-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###########################
MODULE MODI_WRITE_PROFILER_n
!      ###########################
!
INTERFACE
!
      SUBROUTINE WRITE_PROFILER_n(TPDIAFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE ! diachronic file to write
!
END SUBROUTINE WRITE_PROFILER_n
!
END INTERFACE
!
END MODULE MODI_WRITE_PROFILER_n
!
!     ##########################################
      SUBROUTINE WRITE_PROFILER_n(TPDIAFILE)
!     ##########################################
!
!
!!****  *WRITE_PROFILER* - write the balloon and aircraft trajectories and records
!!                      in the diachronic file
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
!!
!!
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Pierre TULET             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original 15/02/2002
!!     2016 : G.DELAUTIER : LIMA
!!              Oct, 2016 (C.Lac) Add visibility diagnostics for fog
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  J. Escobar : 16/08/2018: From Pierre & Maud , correction use CNAMES(JSV-NSV_CHEMBEG+1)
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LUNIT
USE MODD_PARAMETERS
!
USE MODD_TYPE_PROFILER
USE MODD_PROFILER_n
USE MODD_CH_M9_n,         ONLY: CNAMES
USE MODD_CH_AEROSOL,      ONLY: CAERONAMES, LORILAM, JPMODE
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_ICE_C1R3_DESCR,  ONLY: C1R3NAMES
USE MODD_ELEC_DESCR,      ONLY: CELECNAMES
USE MODD_LG,              ONLY: CLGNAMES
USE MODD_DUST,            ONLY: CDUSTNAMES, LDUST, NMODE_DST
USE MODD_SALT,            ONLY: CSALTNAMES, LSALT
USE MODD_NSV
USE MODD_RADIATIONS_n,     ONLY:NAER
USE MODD_DIAG_IN_RUN
!
USE MODE_DUST_PSD
USE MODE_AERO_PSD
!
USE MODI_WRITE_DIACHRO
USE MODD_PARAM_LIMA_WARM, ONLY: CLIMA_WARM_NAMES, CAERO_MASS
USE MODD_PARAM_LIMA_COLD, ONLY: CLIMA_COLD_NAMES
USE MODD_PARAM_LIMA     , ONLY: NINDICE_CCN_IMM,NMOD_CCN,NMOD_IFN,NMOD_IMM
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE ! diachronic file to write
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER     ::  II  ! loop
!
!----------------------------------------------------------------------------
!
DO II=1,NUMBPROFILER
  CALL PROFILER_DIACHRO_n(TPROFILER, II)
ENDDO
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE PROFILER_DIACHRO_n(TPROFILER,II)
!
TYPE(PROFILER),     INTENT(IN)       :: TPROFILER
INTEGER,            INTENT(IN)       :: II
!
!*      0.2  declaration of local variables for diachro
!
REAL,    DIMENSION(:,:,:,:,:,:),  ALLOCATABLE :: ZWORK6   ! contains temporal serie
REAL,    DIMENSION(:,:,:,:,:,:),  ALLOCATABLE :: ZW6      ! contains temporal serie to write
REAL,    DIMENSION(:,:),          ALLOCATABLE :: ZTRAJT   ! localization of the
REAL, DIMENSION(:,:,:,:),         ALLOCATABLE :: ZSV, ZN0, ZSIG, ZRG
REAL, DIMENSION(:,:,:),           ALLOCATABLE :: ZRHO
!
INTEGER, DIMENSION(:),            ALLOCATABLE :: IGRID    ! grid indicator
CHARACTER(LEN=  8)                            :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YUNIT    ! physical unit
!
INTEGER :: IPROC    ! number of variables records
INTEGER :: JPROC    ! loop counter
INTEGER :: JRR      ! loop counter
INTEGER :: JSV      ! loop counter
INTEGER :: IKU, IK  ! loop counter
CHARACTER(LEN=2)  :: INDICE
INTEGER           :: I
!
!----------------------------------------------------------------------------
!
IF (TPROFILER%X(II)==XUNDEF) RETURN
IF (TPROFILER%Y(II)==XUNDEF) RETURN
IKU = SIZE(TPROFILER%W,2)    !nbre de niveaux sur la verticale SIZE(TPROFILER%W,2)
!
IPROC = 24 + SIZE(TPROFILER%R,4) + SIZE(TPROFILER%SV,4)
IF (LDIAG_IN_RUN) IPROC = IPROC + 15
IF (LORILAM) IPROC = IPROC + JPMODE*3
IF (LDUST) IPROC = IPROC + NMODE_DST*3
IF (LDUST .OR. LORILAM .OR. LSALT) IPROC=IPROC+NAER
IF (SIZE(TPROFILER%TKE  )>0) IPROC = IPROC + 1
!
ALLOCATE (ZTRAJT(  SIZE(TPROFILER%TIME),1))
ALLOCATE (ZWORK6(1,1,IKU,SIZE(TPROFILER%TIME),1,IPROC))
ALLOCATE (YCOMMENT(IPROC))
ALLOCATE (YTITLE  (IPROC))
ALLOCATE (YUNIT   (IPROC))
ALLOCATE (IGRID   (IPROC))
!
ZTRAJT  (:,1) = TPROFILER%TIME(:)
!
IGRID  = 1
YGROUP = TPROFILER%NAME(II)
!
!----------------------------------------------------------------------------
DO IK=1, IKU
!
JPROC=0
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'Th'
YUNIT    (JPROC) = 'K'
YCOMMENT (JPROC) = 'Potential temperature' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%TH(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'Thv'
YUNIT    (JPROC) = 'K'
YCOMMENT (JPROC) = 'Virtual Potential temperature' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%THV(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'VISI'
YUNIT    (JPROC) = 'km'
YCOMMENT (JPROC) = 'Visibility'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%VISI(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'VISIKUN'
YUNIT    (JPROC) = 'km'
YCOMMENT (JPROC) = 'Visibility Kunkel'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%VISIKUN(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'RARE'
YUNIT    (JPROC) = 'dBZ'
YCOMMENT (JPROC) = 'Radar reflectivity'       
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%RARE(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'P'
YUNIT    (JPROC) = 'Pascal'
YCOMMENT (JPROC) = 'Pressure' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%P(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'ALT'
YUNIT    (JPROC) = 'm'
YCOMMENT (JPROC) = 'Altitude' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%ZZ(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'LON'
YUNIT    (JPROC) = 'degree'
YCOMMENT (JPROC) = 'Longitude'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%LON(II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'LAT'
YUNIT    (JPROC) = 'degree'
YCOMMENT (JPROC) = 'Latitude'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%LAT(II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'ZON_WIND'
YUNIT    (JPROC) = 'm s-1'
YCOMMENT (JPROC) = 'Zonal wind'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%ZON(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'MER_WIND'
YUNIT    (JPROC) = 'm s-1'
YCOMMENT (JPROC) = 'Meridional wind'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%MER(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'FF'           
YUNIT    (JPROC) = 'm s-1'
YCOMMENT (JPROC) = 'Wind intensity'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%FF(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'DD'        
YUNIT    (JPROC) = 'degree'
YCOMMENT (JPROC) = 'Wind direction'
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%DD(:,IK,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'W'
YUNIT    (JPROC) = 'm s-1'
YCOMMENT (JPROC) = 'Air vertical speed' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%W(:,IK,II)
!
IF (LDIAG_IN_RUN) THEN
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'T2m'
  YUNIT    (JPROC) = 'K'
  YCOMMENT (JPROC) = '2-m temperature' 
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%T2M(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'Q2m'
  YUNIT    (JPROC) = 'kg kg-1'
  YCOMMENT (JPROC) = '2-m humidity' 
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%Q2M(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'HU2m'
  YUNIT    (JPROC) = 'percent'
  YCOMMENT (JPROC) = '2-m relative humidity' 
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%HU2M(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'zon10m'
  YUNIT    (JPROC) = 'm s-1'
  YCOMMENT (JPROC) = '10-m zonal wind' 
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%ZON10M(:,II)
  !       
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'mer10m'
  YUNIT    (JPROC) = 'm s-1'
  YCOMMENT (JPROC) = '10-m meridian wind' 
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%MER10M(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'RN' 
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Net radiation'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%RN(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'H' 
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Sensible heat flux'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%H(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'LE' 
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Total Latent heat flux'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%LE(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'G' 
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Storage heat flux'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%GFLUX(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'SWD'
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Downward short-wave radiation'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SWD(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'SWU'
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Upward short-wave radiation'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SWU(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'LWD'
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Downward long-wave radiation'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%LWD(:,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'LWU'
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Upward long-wave radiation'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%LWU(:,II)
  !
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'TKE_DISS'
  YUNIT    (JPROC) = 'm2 s-2'
  YCOMMENT (JPROC) = 'TKE dissipation rate'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%TKE_DISS(:,IK,II)
  !
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'LEI' 
  YUNIT    (JPROC) = 'W m-2'
  YCOMMENT (JPROC) = 'Solid Latent heat flux'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%LEI(:,II)  
!
ENDIF
!
DO JRR=1,SIZE(TPROFILER%R,4)
  JPROC = JPROC+1
  YUNIT    (JPROC) = 'kg kg-1'
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%R(:,IK,II,JRR)
  IF (JRR==1) THEN
    YTITLE   (JPROC) = 'Rv'
    YCOMMENT (JPROC) = 'Water vapor mixing ratio' 
  ELSE IF (JRR==2) THEN
    YTITLE   (JPROC) = 'Rc'
    YCOMMENT (JPROC) = 'Liquid cloud water mixing ratio' 
  ELSE IF (JRR==3) THEN
    YTITLE   (JPROC) = 'Rr'
    YCOMMENT (JPROC) = 'Rain water mixing ratio' 
  ELSE IF (JRR==4) THEN
    YTITLE   (JPROC) = 'Ri'
    YCOMMENT (JPROC) = 'Ice cloud water mixing ratio' 
  ELSE IF (JRR==5) THEN
    YTITLE   (JPROC) = 'Rs'
    YCOMMENT (JPROC) = 'Snow mixing ratio' 
  ELSE IF (JRR==6) THEN
    YTITLE   (JPROC) = 'Rg'
    YCOMMENT (JPROC) = 'Graupel mixing ratio' 
  ELSE IF (JRR==7) THEN
    YTITLE   (JPROC) = 'Rh'
    YCOMMENT (JPROC) = 'Hail mixing ratio' 
  END IF
END DO
JPROC = JPROC + 1
YTITLE   (JPROC) = 'Rhod'
YUNIT    (JPROC) = 'kg m-3'
YCOMMENT (JPROC) = 'Density of dry air in moist' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%RHOD(:,IK,II)
!
IF (SIZE(TPROFILER%TKE,1)>0) THEN
  JPROC = JPROC+1
  YTITLE   (JPROC) = 'Tke'
  YUNIT    (JPROC) = 'm2 s-2'
  YCOMMENT (JPROC) = 'Turbulent kinetic energy' 
  ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%TKE(:,IK,II)
END IF
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'IWV'
YUNIT    (JPROC) = 'kg m-2'
YCOMMENT (JPROC) = 'Integrated Water Vapour' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%IWV(:,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'ZTD'
YUNIT    (JPROC) = 'm'
YCOMMENT (JPROC) = 'Zenith Tropospheric Delay' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%ZTD(:,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'ZWD'
YUNIT    (JPROC) = 'm'
YCOMMENT (JPROC) = 'Zenith Wet Delay' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%ZWD(:,II)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'ZHD'
YUNIT    (JPROC) = 'm'
YCOMMENT (JPROC) = 'Zenith Hydrostatic Delay' 
ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%ZHD(:,II)
!
IF (SIZE(TPROFILER%SV,4)>=1) THEN
  ! User scalar variables
  DO JSV = 1,NSV_USER
    JPROC = JPROC+1
    WRITE (YTITLE(JPROC),FMT='(A2,I3.3)')   'Sv',JSV
    YUNIT    (JPROC) = 'kg kg-1'
    YCOMMENT (JPROC) = ' ' 
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV)
  END DO
 ! Passive pollutant  scalar variables
  DO JSV = NSV_PPBEG,NSV_PPEND
    JPROC = JPROC+1
    WRITE (YTITLE(JPROC),FMT='(A2,I3.3)')   'Sv',JSV
    YUNIT    (JPROC) = ''
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV)
  END DO
 ! microphysical C2R2 scheme scalar variables
  DO JSV = NSV_C2R2BEG,NSV_C2R2END
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))
    YUNIT    (JPROC) = 'm-3'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV)
  END DO
  ! microphysical C3R5 scheme additional scalar variables
  DO JSV = NSV_C1R3BEG,NSV_C1R3END
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))
    YUNIT    (JPROC) = 'm-3'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV)
  END DO
  ! LIMA variables
  DO JSV=NSV_LIMA_BEG,NSV_LIMA_END
    JPROC = JPROC+1
    YUNIT    (JPROC) = 'kg-1'
    YCOMMENT (JPROC) = ' '
    IF (JSV==NSV_LIMA_NC) YTITLE(JPROC)=TRIM(CLIMA_WARM_NAMES(1))//'T' 
    IF (JSV==NSV_LIMA_NR) YTITLE(JPROC)=TRIM(CLIMA_WARM_NAMES(2))//'T' 
    IF (JSV .GE. NSV_LIMA_CCN_FREE .AND. JSV .LT. NSV_LIMA_CCN_ACTI) THEN
        WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_CCN_FREE + 1)
        YTITLE(JPROC)=TRIM(CLIMA_WARM_NAMES(3))//INDICE//'T'
    ENDIF
    IF (JSV .GE. NSV_LIMA_CCN_ACTI .AND. JSV .LT. NSV_LIMA_CCN_ACTI + NMOD_CCN) THEN
        WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_CCN_ACTI + 1)
        YTITLE(JPROC)=TRIM(CLIMA_WARM_NAMES(4))//INDICE//'T'
    ENDIF
    IF (JSV .EQ. NSV_LIMA_SCAVMASS) YTITLE(JPROC)=TRIM(CAERO_MASS(1))//'T'
    IF (JSV==NSV_LIMA_NI) YTITLE(JPROC)=TRIM(CLIMA_COLD_NAMES(1))//'T' 
    IF (JSV .GE. NSV_LIMA_IFN_FREE .AND. JSV .LT. NSV_LIMA_IFN_NUCL) THEN
        WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_IFN_FREE + 1)
        YTITLE(JPROC)=TRIM(CLIMA_COLD_NAMES(2))//INDICE//'T'
    ENDIF
    IF (JSV .GE. NSV_LIMA_IFN_NUCL .AND. JSV .LT. NSV_LIMA_IFN_NUCL + NMOD_IFN) THEN
        WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_IFN_NUCL + 1)
        YTITLE(JPROC)=TRIM(CLIMA_COLD_NAMES(3))//INDICE//'T'
    ENDIF
    I = 0
    IF (JSV .GE. NSV_LIMA_IMM_NUCL .AND. JSV .LT. NSV_LIMA_IMM_NUCL + NMOD_IMM) THEN
        I = I + 1
        WRITE(INDICE,'(I2.2)')(NINDICE_CCN_IMM(I))
        YTITLE(JPROC)=TRIM(CLIMA_COLD_NAMES(4))//INDICE//'T'
    ENDIF
    IF (JSV .EQ. NSV_LIMA_HOM_HAZE) YTITLE(JPROC)=TRIM(CLIMA_COLD_NAMES(5))//'T'

    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV)
  END DO 
  ! electrical scalar variables
  DO JSV = NSV_ELECBEG,NSV_ELECEND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))
    YUNIT    (JPROC) = 'C'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV)
  END DO
  ! chemical scalar variables
  DO JSV=NSV_CHEMBEG,NSV_CHEMEND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CNAMES(JSV-NSV_CHEMBEG+1))
    YUNIT    (JPROC) = 'ppb'
    WRITE(YCOMMENT (JPROC),'(A5,A3,I3.3)') 'T(s) ','SVT',JSV
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV) * 1.E9
  END DO
  IF ((LORILAM).AND. .NOT.(ANY(TPROFILER%P(:,IK,II) == 0.))) THEN
    ALLOCATE (ZSV(1,1,SIZE(TPROFILER%TIME),NSV_AER)) 
    ALLOCATE (ZRHO(1,1,SIZE(TPROFILER%TIME))) 
    ALLOCATE (ZN0(1,1,SIZE(TPROFILER%TIME),JPMODE)) 
    ALLOCATE (ZRG(1,1,SIZE(TPROFILER%TIME),JPMODE)) 
    ALLOCATE (ZSIG(1,1,SIZE(TPROFILER%TIME),JPMODE)) 
    ZSV(1,1,:,1:NSV_AER) = TPROFILER%SV(:,IK,II,NSV_AERBEG:NSV_AEREND)
    IF (SIZE(TPROFILER%R,4) >0) THEN
      ZRHO(1,1,:) = 0.
      DO JRR=1,SIZE(TPROFILER%R,4)
        ZRHO(1,1,:) = ZRHO(1,1,:) + TPROFILER%R(:,IK,II,JRR)
      ENDDO
      ZRHO(1,1,:) = TPROFILER%TH(:,IK,II) * ( 1. + XRV/XRD*TPROFILER%R(:,IK,II,1) )  &
                                          / ( 1. + ZRHO(1,1,:)                ) 
    ELSE
      ZRHO(1,1,:) = TPROFILER%TH(:,IK,II)
    ENDIF
    ZRHO(1,1,:) =  TPROFILER%P(:,IK,II) / &
                  (XRD *ZRHO(1,1,:) *((TPROFILER%P(:,IK,II)/XP00)**(XRD/XCPD)) )
    CALL PPP2AERO(ZSV,ZRHO, PSIG3D=ZSIG, PRG3D=ZRG, PN3D=ZN0)
    DO JSV=1,JPMODE
      ! mean radius
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'AERRGA',JSV
      YUNIT    (JPROC) = 'um'
      WRITE(YCOMMENT(JPROC),'(A18,I1)')'RG (nb) AERO MODE ',JSV
      ZWORK6 (1,1,IK,:,1,JPROC) = ZRG(1,1,:,JSV)
      ! standard deviation
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A7,I1)')'AERSIGA',JSV
      YUNIT    (JPROC) = '  '
      WRITE(YCOMMENT(JPROC),'(A16,I1)')'SIGMA AERO MODE ',JSV
      ZWORK6 (1,1,IK,:,1,JPROC) = ZSIG(1,1,:,JSV)
      ! particles number
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'AERN0A',JSV
      YUNIT    (JPROC) = 'm-3'
      WRITE(YCOMMENT(JPROC),'(A13,I1)')'N0 AERO MODE ',JSV
      ZWORK6 (1,1,IK,:,1,JPROC) = ZN0(1,1,:,JSV)
    ENDDO
    DEALLOCATE (ZSV,ZRHO) 
    DEALLOCATE (ZN0,ZRG,ZSIG) 
  END IF
  ! dust scalar variables
  DO JSV = NSV_DSTBEG,NSV_DSTEND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CDUSTNAMES(JSV-NSV_DSTBEG+1))
    YUNIT    (JPROC) = 'ppb'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV) * 1.E9
  END DO
  IF ((LDUST).AND. .NOT.(ANY(TPROFILER%P(:,IK,II) == 0.))) THEN
    ALLOCATE (ZSV(1,1,SIZE(TPROFILER%TIME),NSV_DST)) 
    ALLOCATE (ZRHO(1,1,SIZE(TPROFILER%TIME))) 
    ALLOCATE (ZN0(1,1,SIZE(TPROFILER%TIME),NMODE_DST)) 
    ALLOCATE (ZRG(1,1,SIZE(TPROFILER%TIME),NMODE_DST)) 
    ALLOCATE (ZSIG(1,1,SIZE(TPROFILER%TIME),NMODE_DST)) 
    ZSV(1,1,:,1:NSV_DST) = TPROFILER%SV(:,IK,II,NSV_DSTBEG:NSV_DSTEND)
    IF (SIZE(TPROFILER%R,4) >0) THEN
      ZRHO(1,1,:) = 0.
      DO JRR=1,SIZE(TPROFILER%R,4)
        ZRHO(1,1,:) = ZRHO(1,1,:) + TPROFILER%R(:,IK,II,JRR)
      ENDDO
      ZRHO(1,1,:) = TPROFILER%TH(:,IK,II) * ( 1. + XRV/XRD*TPROFILER%R(:,IK,II,1) )  &
                                          / ( 1. + ZRHO(1,1,:)                ) 
    ELSE
      ZRHO(1,1,:) = TPROFILER%TH(:,IK,II)
    ENDIF
    ZRHO(1,1,:) =  TPROFILER%P(:,IK,II) / &
                  (XRD *ZRHO(1,1,:) *((TPROFILER%P(:,IK,II)/XP00)**(XRD/XCPD)) )
    CALL PPP2DUST(ZSV,ZRHO, PSIG3D=ZSIG, PRG3D=ZRG, PN3D=ZN0)
    DO JSV=1,NMODE_DST
      ! mean radius
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'DSTRGA',JSV
      YUNIT    (JPROC) = 'um'
      WRITE(YCOMMENT(JPROC),'(A18,I1)')'RG (nb) DUST MODE ',JSV
      ZWORK6 (1,1,IK,:,1,JPROC) = ZRG(1,1,:,JSV)
      ! standard deviation
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A7,I1)')'DSTSIGA',JSV
      YUNIT    (JPROC) = '  '
      WRITE(YCOMMENT(JPROC),'(A16,I1)')'SIGMA DUST MODE ',JSV
      ZWORK6 (1,1,IK,:,1,JPROC) = ZSIG(1,1,:,JSV)
      ! particles number
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'DSTN0A',JSV
      YUNIT    (JPROC) = 'm-3'
      WRITE(YCOMMENT(JPROC),'(A13,I1)')'N0 DUST MODE ',JSV
      ZWORK6 (1,1,IK,:,1,JPROC) = ZN0(1,1,:,JSV)
    ENDDO
    DEALLOCATE (ZSV,ZRHO) 
    DEALLOCATE (ZN0,ZRG,ZSIG) 
  END IF
  ! sea salt scalar variables
  DO JSV = NSV_SLTBEG,NSV_SLTEND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CSALTNAMES(JSV-NSV_SLTBEG+1))
    YUNIT    (JPROC) = 'ppb'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%SV(:,IK,II,JSV) * 1.E9
  END DO
  IF (LDUST .OR. LORILAM .OR. LSALT) THEN
  DO JSV = 1,NAER
    JPROC = JPROC+1
    WRITE(YTITLE(JPROC),'(A6,I1)')'AEREXT',JSV
    YUNIT    (JPROC) = ' '
    YCOMMENT (JPROC) = 'Aerosol Extinction'
    ZWORK6 (1,1,IK,:,1,JPROC) = TPROFILER%AER(:,IK,II,JSV) 
  END DO
  ENDIF
ENDIF
!
END DO
!
!----------------------------------------------------------------------------
!

ALLOCATE (ZW6(1,1,IKU,SIZE(TPROFILER%TIME),1,JPROC))
ZW6 = ZWORK6(:,:,:,:,:,:JPROC)
DEALLOCATE(ZWORK6)

CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"CART",IGRID(:JPROC), TPROFILER%DATIME,&
                   ZW6,ZTRAJT,YTITLE(:JPROC),YUNIT(:JPROC),YCOMMENT(:JPROC),     &
                   .TRUE.,.TRUE.,.FALSE.,                                        &
                   KIL=1,KIH=1,KJL=1,KJH=1,KKL=1,KKH=IKU                         )
!
DEALLOCATE (ZTRAJT  )
DEALLOCATE (ZW6     )
DEALLOCATE (YCOMMENT)
DEALLOCATE (YTITLE  )
DEALLOCATE (YUNIT   )
DEALLOCATE (IGRID   )
!----------------------------------------------------------------------------
END SUBROUTINE PROFILER_DIACHRO_n
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
END SUBROUTINE WRITE_PROFILER_n
