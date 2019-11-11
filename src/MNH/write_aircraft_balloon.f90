!MNH_LIC Copyright 2000-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###########################
MODULE MODI_WRITE_AIRCRAFT_BALLOON
!      ###########################
!
INTERFACE
!
      SUBROUTINE WRITE_AIRCRAFT_BALLOON(TPDIAFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA), INTENT(IN) :: TPDIAFILE ! file to write
!
END SUBROUTINE WRITE_AIRCRAFT_BALLOON
!
END INTERFACE
!
END MODULE MODI_WRITE_AIRCRAFT_BALLOON
!
!     ##########################################
      SUBROUTINE WRITE_AIRCRAFT_BALLOON(TPDIAFILE)
!     ##########################################
!
!
!!****  *WRITE_AIRCRAFT_BALLOON* - write the balloon and aircraft trajectories and records
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
!!      Valery Masson             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original 15/05/2000
!!     10/01/2011 adding IMI, the model number
!!     March, 2013 :  C.Lac : add vertical profiles
!!     July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!     Oct 2016 : G.Delautier LIMA
!!     August 2016 (M.Leriche) Add mass concentration of aerosol species
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!     P. Wautelet 29/01/2019: bug: moved an instruction later (to prevent access to a not allocated array)
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
USE MODD_AIRCRAFT_BALLOON
USE MODD_CH_M9_n,         ONLY: CNAMES
USE MODD_CH_AEROSOL,      ONLY: CAERONAMES, LORILAM, NSP, NCARB, NSOA,    & 
                                JPMODE, JP_AER_BC, JP_AER_OC, JP_AER_DST, &
                                JP_AER_H2O, JP_AER_SO4, JP_AER_NO3,       &
                                JP_AER_NH3, JP_AER_SOA1, JP_AER_SOA2,     &
                                JP_AER_SOA3, JP_AER_SOA4, JP_AER_SOA5,    &
                                JP_AER_SOA6, JP_AER_SOA7, JP_AER_SOA8,    &
                                JP_AER_SOA9, JP_AER_SOA10
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_ICE_C1R3_DESCR,  ONLY: C1R3NAMES
USE MODD_ELEC_DESCR,      ONLY: CELECNAMES
USE MODD_LG,              ONLY: CLGNAMES
USE MODD_DUST,            ONLY: CDUSTNAMES, LDUST, NMODE_DST
USE MODD_SALT,            ONLY: CSALTNAMES
USE MODD_NSV
USE MODD_DIAG_IN_RUN
USE MODD_PARAM_LIMA_WARM, ONLY: CLIMA_WARM_NAMES, CAERO_MASS
USE MODD_PARAM_LIMA_COLD, ONLY: CLIMA_COLD_NAMES
USE MODD_PARAM_LIMA     , ONLY: NINDICE_CCN_IMM,NMOD_CCN,NMOD_IFN,NMOD_IMM
!
USE MODE_MODELN_HANDLER
USE MODE_DUST_PSD
USE MODE_AERO_PSD
!
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA), INTENT(IN) :: TPDIAFILE ! file to write
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER :: IMI    ! current model index
!
!----------------------------------------------------------------------------
!
IMI=GET_CURRENT_MODEL_INDEX()
!
CALL FLYER_DIACHRO(TBALLOON1)
CALL FLYER_DIACHRO(TBALLOON2)
CALL FLYER_DIACHRO(TBALLOON3)
CALL FLYER_DIACHRO(TBALLOON4)
CALL FLYER_DIACHRO(TBALLOON5)
CALL FLYER_DIACHRO(TBALLOON6)
CALL FLYER_DIACHRO(TBALLOON7)
CALL FLYER_DIACHRO(TBALLOON8)
CALL FLYER_DIACHRO(TBALLOON9)
!
CALL FLYER_DIACHRO(TAIRCRAFT1)
CALL FLYER_DIACHRO(TAIRCRAFT2)
CALL FLYER_DIACHRO(TAIRCRAFT3)
CALL FLYER_DIACHRO(TAIRCRAFT4)
CALL FLYER_DIACHRO(TAIRCRAFT5)
CALL FLYER_DIACHRO(TAIRCRAFT6)
CALL FLYER_DIACHRO(TAIRCRAFT7)
CALL FLYER_DIACHRO(TAIRCRAFT8)
CALL FLYER_DIACHRO(TAIRCRAFT9)
CALL FLYER_DIACHRO(TAIRCRAFT10)
CALL FLYER_DIACHRO(TAIRCRAFT11)
CALL FLYER_DIACHRO(TAIRCRAFT12)
CALL FLYER_DIACHRO(TAIRCRAFT13)
CALL FLYER_DIACHRO(TAIRCRAFT14)
CALL FLYER_DIACHRO(TAIRCRAFT15)
CALL FLYER_DIACHRO(TAIRCRAFT16)
CALL FLYER_DIACHRO(TAIRCRAFT17)
CALL FLYER_DIACHRO(TAIRCRAFT18)
CALL FLYER_DIACHRO(TAIRCRAFT19)
CALL FLYER_DIACHRO(TAIRCRAFT20)
CALL FLYER_DIACHRO(TAIRCRAFT21)
CALL FLYER_DIACHRO(TAIRCRAFT22)
CALL FLYER_DIACHRO(TAIRCRAFT23)
CALL FLYER_DIACHRO(TAIRCRAFT24)
CALL FLYER_DIACHRO(TAIRCRAFT25)
CALL FLYER_DIACHRO(TAIRCRAFT26)
CALL FLYER_DIACHRO(TAIRCRAFT27)
CALL FLYER_DIACHRO(TAIRCRAFT28)
CALL FLYER_DIACHRO(TAIRCRAFT29)
CALL FLYER_DIACHRO(TAIRCRAFT30)
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
SUBROUTINE FLYER_DIACHRO(TPFLYER)
!
TYPE(FLYER),        INTENT(IN)       :: TPFLYER
!
!*      0.2  declaration of local variables for diachro
!
REAL, DIMENSION(:,:),         ALLOCATABLE :: ZTRAJT ! localization of the
REAL, DIMENSION(:,:,:),       ALLOCATABLE :: ZTRAJX ! temporal series
REAL, DIMENSION(:,:,:),       ALLOCATABLE :: ZTRAJY ! in t,x,y and z.
REAL, DIMENSION(:,:,:),       ALLOCATABLE :: ZTRAJZ !
REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ZWORK6 ! contains temporal serie
REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ZW6    ! contains temporal serie to write
REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ZWORKZ6! contains temporal serie
REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ZWZ6   ! contains temporal serie
REAL, DIMENSION(:,:,:,:),     ALLOCATABLE :: ZSV, ZN0, ZSIG, ZRG
REAL, DIMENSION(:,:,:,:,:),   ALLOCATABLE :: ZPTOTA
REAL, DIMENSION(:,:,:),       ALLOCATABLE :: ZRHO
!
INTEGER, DIMENSION(:),            ALLOCATABLE :: IGRID    ! grid indicator
CHARACTER(LEN=  8)                            :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YUNIT    ! physical unit
!
INTEGER :: IPROC    ! number of variables records
INTEGER :: JPROC    ! loop counter
INTEGER, DIMENSION(:),            ALLOCATABLE :: IGRIDZ   ! grid indicator
CHARACTER(LEN=  8)                            :: YGROUPZ  ! group title
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YCOMMENTZ! comment string
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YTITLEZ  ! title
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YUNITZ   ! physical unit
INTEGER :: IPROCZ   ! number of variables records
INTEGER :: JPROCZ   ! loop counter
INTEGER :: JRR      ! loop counter
INTEGER :: JSV      ! loop counter
INTEGER :: JPT      ! loop counter
INTEGER :: IKU, IK
CHARACTER(LEN=2)  :: INDICE
INTEGER           :: I
INTEGER :: JLOOP
!
!----------------------------------------------------------------------------
!
IF (TPFLYER%NMODEL==0) RETURN
IF (ALL(TPFLYER%X==XUNDEF)) RETURN
IF (COUNT(TPFLYER%X/=XUNDEF)<=1) RETURN
IF ( IMI /= TPFLYER%NMODEL ) RETURN
!
IKU = SIZE(TPFLYER%RTZ,2) !number of vertical levels
!
IPROC = 20 + SIZE(TPFLYER%R,2) + SIZE(TPFLYER%SV,2) &
       + 2 + SIZE(TPFLYER%SVW_FLUX,2)
IPROCZ = SIZE(TPFLYER%RTZ,2)+ SIZE(TPFLYER%RZ,2)+ SIZE(TPFLYER%RZ,3)+  SIZE(TPFLYER%CRARE,2)+ &
         SIZE(TPFLYER%CRARE_ATT,2)+ SIZE(TPFLYER%WZ,2) + SIZE(TPFLYER%FFZ,2)+ &
         SIZE(TPFLYER%IWCZ,2)+ SIZE(TPFLYER%LWCZ,2) + SIZE(TPFLYER%CIZ,2) + &
         SIZE(TPFLYER%ZZ,2)

IF (NSV_LIMA_BEG/=NSV_LIMA_END) IPROCZ= IPROCZ+ SIZE(TPFLYER%CCZ,2) + SIZE(TPFLYER%CRZ,2)
IF (SIZE(TPFLYER%TKE  )>0) IPROC = IPROC + 1
IF (LDIAG_IN_RUN) IPROC = IPROC + 1
IF (LORILAM) IPROC = IPROC + JPMODE*3
IF (LDUST) IPROC = IPROC + NMODE_DST*3
IF (SIZE(TPFLYER%TSRAD)>0) IPROC = IPROC + 1
!
ALLOCATE (ZTRAJT(  SIZE(TPFLYER%TIME),1))
ALLOCATE (ZTRAJX(1,SIZE(TPFLYER%TIME),1))
ALLOCATE (ZTRAJY(1,SIZE(TPFLYER%TIME),1))
ALLOCATE (ZTRAJZ(1,SIZE(TPFLYER%TIME),1))
ALLOCATE (ZWORK6(1,1,1,SIZE(TPFLYER%TIME),1,IPROC))
ALLOCATE (YCOMMENT(IPROC))
ALLOCATE (YTITLE  (IPROC))
ALLOCATE (YUNIT   (IPROC))
ALLOCATE (IGRID   (IPROC))
ALLOCATE (ZWORKZ6(1,1,IKU,SIZE(TPFLYER%TIME),1,IPROCZ))
ALLOCATE (YCOMMENTZ(IPROCZ))
ALLOCATE (YTITLEZ (IPROCZ))
ALLOCATE (YUNITZ  (IPROCZ))
ALLOCATE (IGRIDZ  (IPROCZ))

!
ZTRAJT  (:,1) = TPFLYER%TIME
ZTRAJX(1,:,1) = TPFLYER%X
ZTRAJY(1,:,1) = TPFLYER%Y
ZTRAJZ(1,:,1) = TPFLYER%Z
!
!
IGRID  = 1
YGROUP = TPFLYER%TITLE
IGRIDZ = 1
YGROUPZ = TRIM(TPFLYER%TITLE)//"Z"
!
!----------------------------------------------------------------------------
JPROC = 0
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'ZS'
YUNIT    (JPROC) = 'm'
YCOMMENT (JPROC) = 'orography'
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%ZS(:)
!
IF (TPFLYER%ALTDEF) THEN
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'P'
  YUNIT    (JPROC) = 'Pascal'
  YCOMMENT (JPROC) = 'pressure' 
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%P(:)
ELSE
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'Z'
  YUNIT    (JPROC) = 'm'
  YCOMMENT (JPROC) = 'altitude' 
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%Z(:)
ENDIF
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'LON'
YUNIT    (JPROC) = 'degree'
YCOMMENT (JPROC) = 'longitude'
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%XLON(:)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'LAT'
YUNIT    (JPROC) = 'degree'
YCOMMENT (JPROC) = 'latitude'
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%YLAT(:)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'ZON_WIND'
YUNIT    (JPROC) = 'm s-1'
YCOMMENT (JPROC) = 'zonal wind'
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%ZON(:)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'MER_WIND'
YUNIT    (JPROC) = 'm s-1'
YCOMMENT (JPROC) = 'meridian wind'
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%MER(:)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'W'
YUNIT    (JPROC) = 'm s-1'
YCOMMENT (JPROC) = 'air vertical speed' 
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%W(:)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'Th'
YUNIT    (JPROC) = 'K'
YCOMMENT (JPROC) = 'potential temperature' 
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%TH(:)
!
DO JRR=1,SIZE(TPFLYER%R,2)
  JPROC = JPROC+1
  YUNIT    (JPROC) = 'kg kg-1'
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%R(:,JRR)
  IF (JRR==1) THEN
    YTITLE   (JPROC) = 'Rv'
    YCOMMENT (JPROC) = 'water vapor mixing ratio' 
  ELSE IF (JRR==2) THEN
    YTITLE   (JPROC) = 'Rc'
    YCOMMENT (JPROC) = 'liquid cloud water mixing ratio' 
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
!
!add cloud liquid water content in g/m3 to compare to measurements from FSSP
!IF (.NOT.(ANY(TPFLYER%P(:) == 0.))) THEN
ALLOCATE (ZRHO(1,1,SIZE(TPFLYER%TIME)))
IF (SIZE(TPFLYER%R,2) >1) THEN !cloud water is present
  ZRHO(1,1,:) = 0.
  DO JRR=1,SIZE(TPFLYER%R,2)
    ZRHO(1,1,:) = ZRHO(1,1,:) + TPFLYER%R(:,JRR)
  ENDDO
  ZRHO(1,1,:) = TPFLYER%TH(:) * ( 1. + XRV/XRD*TPFLYER%R(:,1) )  &
                                / ( 1. + ZRHO(1,1,:)              )
  DO JPT=1,SIZE(TPFLYER%TIME)
    IF (TPFLYER%P(JPT) == 0.) THEN
      ZRHO(1,1,JPT) = 0.
    ELSE
      ZRHO(1,1,JPT) =  TPFLYER%P(JPT) / &
               (XRD *ZRHO(1,1,JPT) *((TPFLYER%P(JPT)/XP00)**(XRD/XCPD))  )
    ENDIF
  ENDDO
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'LWC'
  YUNIT    (JPROC) = 'g m-3'
  YCOMMENT (JPROC) = 'cloud liquid water content'
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%R(:,2)*ZRHO(1,1,:)*1.E3
  DEALLOCATE (ZRHO)
ENDIF
!ENDIF
!
IF (SIZE(TPFLYER%TKE)>0) THEN
  JPROC = JPROC+1
  YTITLE   (JPROC) = 'Tke'
  YUNIT    (JPROC) = 'm2 s-2'
  YCOMMENT (JPROC) = 'Turbulent kinetic energy' 
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%TKE(:)
END IF
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'H_FLUX'
YUNIT    (JPROC) = 'W m-2'
YCOMMENT (JPROC) = 'sensible flux' 
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%THW_FLUX(:)
!
JPROC = JPROC + 1
YTITLE   (JPROC) = 'LE_FLUX'
YUNIT    (JPROC) = 'W m-2'
YCOMMENT (JPROC) = 'latent flux' 
ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%RCW_FLUX(:)
!
DO JSV=1,SIZE(TPFLYER%SVW_FLUX,2)
  JPROC = JPROC + 1
  YTITLE   (JPROC) = 'SV_FLUX'
  YUNIT    (JPROC) = 'SVUNIT m s-1'
  YCOMMENT (JPROC) = 'scalar flux' 
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SVW_FLUX(:,JSV)
END DO
IF (LDIAG_IN_RUN) THEN
  JPROC = JPROC+1
  YTITLE   (JPROC) = 'Tke_Diss'
  YUNIT    (JPROC) = 'm2 s-2'
  YCOMMENT (JPROC) = 'TKE dissipation rate' 
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%TKE_DISS(:)
ENDIF
!
IF (SIZE(TPFLYER%SV,2)>=1) THEN
  ! User scalar variables
  DO JSV = 1,NSV_USER
    JPROC = JPROC+1
    WRITE (YTITLE(JPROC),FMT='(A2,I3.3)')   'Sv',JSV
    YUNIT    (JPROC) = 'kg kg-1'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV)
  END DO
  ! microphysical C2R2 scheme scalar variables
  DO JSV = NSV_C2R2BEG,NSV_C2R2END
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))
    YUNIT    (JPROC) = 'm-3'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV)
  END DO
  ! microphysical C3R5 scheme additional scalar variables
  DO JSV = NSV_C1R3BEG,NSV_C1R3END
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))
    YUNIT    (JPROC) = 'm-3'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV)
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
    IF (JSV .EQ. NSV_LIMA_SCAVMASS) THEN
      YTITLE(JPROC)=TRIM(CAERO_MASS(1))//'T'
      YUNIT    (JPROC) = 'kg kg-1'
    ENDIF
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
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV)
  END DO 
  ! electrical scalar variables
  DO JSV = NSV_ELECBEG,NSV_ELECEND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))
    YUNIT    (JPROC) = 'C'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV)
  END DO
  ! chemical scalar variables
  DO JSV = NSV_CHEMBEG,NSV_CHEMEND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CNAMES(JSV-NSV_CHEMBEG+1))
    YUNIT    (JPROC) = 'ppb'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV) * 1.E9
  END DO
  ! LiNOX passive tracer
  DO JSV = NSV_LNOXBEG,NSV_LNOXEND
    JPROC = JPROC+1
    WRITE (YTITLE(JPROC),FMT='(A5)') 'LiNOx'
    YUNIT    (JPROC) = 'ppb'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV) * 1.E9
  END DO
  ! aerosol scalar variables
  DO JSV = NSV_AERBEG,NSV_AEREND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CAERONAMES(JSV-NSV_AERBEG+1))
    YUNIT    (JPROC) = 'ppb'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV) * 1.E9
  END DO
  IF ((LORILAM).AND. .NOT.(ANY(TPFLYER%P(:) == 0.))) THEN

    ALLOCATE (ZSV(1,1,SIZE(TPFLYER%TIME),NSV_AER)) 
    ALLOCATE (ZRHO(1,1,SIZE(TPFLYER%TIME))) 
    ALLOCATE (ZN0(1,1,SIZE(TPFLYER%TIME),JPMODE)) 
    ALLOCATE (ZRG(1,1,SIZE(TPFLYER%TIME),JPMODE)) 
    ALLOCATE (ZSIG(1,1,SIZE(TPFLYER%TIME),JPMODE)) 
    ALLOCATE (ZPTOTA(1,1,SIZE(TPFLYER%TIME),NSP+NCARB+NSOA,JPMODE))    
    ZSV(1,1,:,1:NSV_AER) = TPFLYER%SV(:,NSV_AERBEG:NSV_AEREND)
    IF (SIZE(TPFLYER%R,2) >0) THEN
      ZRHO(1,1,:) = 0.
      DO JRR=1,SIZE(TPFLYER%R,2)
        ZRHO(1,1,:) = ZRHO(1,1,:) + TPFLYER%R(:,JRR)
      ENDDO
      ZRHO(1,1,:) = TPFLYER%TH(:) * ( 1. + XRV/XRD*TPFLYER%R(:,1) )  &
                                  / ( 1. + ZRHO(1,1,:)                ) 
    ELSE
      ZRHO(1,1,:) = TPFLYER%TH(:)
    ENDIF
    ZRHO(1,1,:) =  TPFLYER%P(:) / &
                  (XRD *ZRHO(1,1,:) *((TPFLYER%P(:)/XP00)**(XRD/XCPD))  )
    ZSIG = 0.
    ZRG = 0.
    ZN0 = 0.
    ZPTOTA = 0.
    DO JPT=1,SIZE(TPFLYER%TIME) ! prevent division by zero if ZSV = 0.
      IF (ALL(ZSV(1,1,JPT,:)/=0.)) THEN
        CALL PPP2AERO(ZSV,ZRHO, PSIG3D=ZSIG, PRG3D=ZRG, PN3D=ZN0, PCTOTA=ZPTOTA)
      ENDIF
    ENDDO
    DO JSV=1,JPMODE
      ! mean radius
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'AERRGA',JSV
      YUNIT    (JPROC) = 'um'
      WRITE(YCOMMENT(JPROC),'(A18,I1)')'RG (nb) AERO MODE ',JSV
      ZWORK6 (1,1,1,:,1,JPROC) = ZRG(1,1,:,JSV)
      ! standard deviation
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A7,I1)')'AERSIGA',JSV
      YUNIT    (JPROC) = '  '
      WRITE(YCOMMENT(JPROC),'(A16,I1)')'SIGMA AERO MODE ',JSV
      ZWORK6 (1,1,1,:,1,JPROC) = ZSIG(1,1,:,JSV)
      ! particles number
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'AERN0A',JSV
      YUNIT    (JPROC) = 'm-3'
      WRITE(YCOMMENT(JPROC),'(A13,I1)')'N0 AERO MODE ',JSV
      ZWORK6 (1,1,1,:,1,JPROC) = ZN0(1,1,:,JSV)
      ! mass concentration in microg/m3
      ! sulfate
      JPROC = JPROC + 1
      WRITE(YTITLE(JPROC),'(A4,I1)')'MSO4',JSV
      YUNIT    (JPROC) = 'ug m-3'
      WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SO4 AEROSOL MODE ',JSV
      ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SO4,JSV)
      ! nitrate
      JPROC = JPROC + 1
      WRITE(YTITLE(JPROC),'(A4,I1)')'MNO3',JSV
      YUNIT    (JPROC) = 'ug m-3'
      WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS NO3 AEROSOL MODE ',JSV
      ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_NO3,JSV)
      ! amoniac
      JPROC = JPROC + 1
      WRITE(YTITLE(JPROC),'(A4,I1)')'MNH3',JSV
      YUNIT    (JPROC) = 'ug m-3'
      WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS NH3 AEROSOL MODE ',JSV
      ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_NH3,JSV)
      ! water
      JPROC = JPROC + 1
      WRITE(YTITLE(JPROC),'(A4,I1)')'MH2O',JSV
      YUNIT    (JPROC) = 'ug m-3'
      WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS H2O AEROSOL MODE ',JSV
      ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_H2O,JSV)
      IF (NSOA .EQ. 10) THEN
        ! SOA1
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA1',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA1 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA1,JSV)
        ! SOA2
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA2',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA2 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA2,JSV)
        ! SOA3
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA3',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA3 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA3,JSV)
        ! SOA4
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA4',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA4 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA4,JSV)
        ! SOA5
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA5',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA5 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA5,JSV)
        ! SOA6
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA6',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA6 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA6,JSV)
        ! SOA7
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA7',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA7 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA7,JSV)
        ! SOA8
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA8',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA8 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA8,JSV)
        ! SOA9
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA9',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA9 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA9,JSV)
        ! SOA10
        JPROC = JPROC + 1
        WRITE(YTITLE(JPROC),'(A4,I1)')'MSOA10',JSV
        YUNIT    (JPROC) = 'ug m-3'
        WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS SOA10 AEROSOL MODE ',JSV
        ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_SOA10,JSV)
      ENDIF
      ! OC
      JPROC = JPROC + 1
      WRITE(YTITLE(JPROC),'(A4,I1)')'MOC',JSV
      YUNIT    (JPROC) = 'ug m-3'
      WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS OC AEROSOL MODE ',JSV
      ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_OC,JSV)
      ! BC
      JPROC = JPROC + 1
      WRITE(YTITLE(JPROC),'(A4,I1)')'MBC',JSV
      YUNIT    (JPROC) = 'ug m-3'
      WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS BC AEROSOL MODE ',JSV
      ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_BC,JSV)
      ! dust
      JPROC = JPROC + 1
      WRITE(YTITLE(JPROC),'(A4,I1)')'MDUST',JSV
      YUNIT    (JPROC) = 'ug m-3'
      WRITE(YCOMMENT(JPROC),'(A22,I1)')'MASS DUST AEROSOL MODE ',JSV
      ZWORK6(1,1,1,:,1,JPROC) = ZPTOTA(1,1,:,JP_AER_DST,JSV)      
    ENDDO
    DEALLOCATE (ZSV,ZRHO) 
    DEALLOCATE (ZN0,ZRG,ZSIG,ZPTOTA) 
  END IF
! dust scalar variables
  DO JSV = NSV_DSTBEG,NSV_DSTEND
    JPROC = JPROC+1
    YTITLE(JPROC)= TRIM(CDUSTNAMES(JSV-NSV_DSTBEG+1))
    YUNIT    (JPROC) = 'ppb'
    YCOMMENT (JPROC) = ' '
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV) * 1.E9
  END DO
  IF ((LDUST).AND. .NOT.(ANY(TPFLYER%P(:) == 0.))) THEN
    ALLOCATE (ZSV(1,1,SIZE(TPFLYER%TIME),NSV_DST)) 
    ALLOCATE (ZRHO(1,1,SIZE(TPFLYER%TIME))) 
    ALLOCATE (ZN0(1,1,SIZE(TPFLYER%TIME),NMODE_DST)) 
    ALLOCATE (ZRG(1,1,SIZE(TPFLYER%TIME),NMODE_DST)) 
    ALLOCATE (ZSIG(1,1,SIZE(TPFLYER%TIME),NMODE_DST)) 
    ZSV(1,1,:,1:NSV_DST) = TPFLYER%SV(:,NSV_DSTBEG:NSV_DSTEND)
    IF (SIZE(TPFLYER%R,2) >0) THEN
      ZRHO(1,1,:) = 0.
      DO JRR=1,SIZE(TPFLYER%R,2)
        ZRHO(1,1,:) = ZRHO(1,1,:) + TPFLYER%R(:,JRR)
      ENDDO
      ZRHO(1,1,:) = TPFLYER%TH(:) * ( 1. + XRV/XRD*TPFLYER%R(:,1) )  &
                                          / ( 1. + ZRHO(1,1,:)                ) 
    ELSE
      ZRHO(1,1,:) = TPFLYER%TH(:)
    ENDIF
    ZRHO(1,1,:) =  TPFLYER%P(:) / &
                  (XRD *ZRHO(1,1,:) *((TPFLYER%P(:)/XP00)**(XRD/XCPD)) )
    CALL PPP2DUST(ZSV,ZRHO, PSIG3D=ZSIG, PRG3D=ZRG, PN3D=ZN0)
    DO JSV=1,NMODE_DST
      ! mean radius
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'DSTRGA',JSV
      YUNIT    (JPROC) = 'um'
      WRITE(YCOMMENT(JPROC),'(A18,I1)')'RG (nb) DUST MODE ',JSV
      ZWORK6 (1,1,1,:,1,JPROC) = ZRG(1,1,:,JSV)
      ! standard deviation
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A7,I1)')'DSTSIGA',JSV
      YUNIT    (JPROC) = '  '
      WRITE(YCOMMENT(JPROC),'(A16,I1)')'SIGMA DUST MODE ',JSV
      ZWORK6 (1,1,1,:,1,JPROC) = ZSIG(1,1,:,JSV)
      ! particles number
      JPROC = JPROC+1
      WRITE(YTITLE(JPROC),'(A6,I1)')'DSTN0A',JSV
      YUNIT    (JPROC) = 'm-3'
      WRITE(YCOMMENT(JPROC),'(A13,I1)')'N0 DUST MODE ',JSV
      ZWORK6 (1,1,1,:,1,JPROC) = ZN0(1,1,:,JSV)
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
    ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%SV(:,JSV) * 1.E9
  END DO
ENDIF
!
IF (SIZE(TPFLYER%TSRAD)>0) THEN
  JPROC = JPROC+1
  YTITLE   (JPROC) = 'Tsrad'
  YUNIT    (JPROC) = 'K'
  YCOMMENT (JPROC) = 'Radiative Surface Temperature'
  ZWORK6 (1,1,1,:,1,JPROC) = TPFLYER%TSRAD(:)
END IF
!
DO IK=1, IKU
!
  JPROCZ=0
!
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'Rt'
  YUNITZ   (JPROCZ) = 'kg kg-1'
  YCOMMENTZ(JPROCZ) = '1D Total hydrometeor mixing ratio'
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%RTZ(:,IK)
!
  DO JRR=1,SIZE(TPFLYER%RZ,3)
    JPROCZ = JPROCZ+1
    YUNITZ    (JPROCZ) = 'kg kg-1'
    ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%RZ(:,IK,JRR)
    IF (JRR==1) THEN
      YTITLEZ   (JPROCZ) = 'Rv'
      YCOMMENTZ (JPROCZ) = '1D water vapor mixing ratio' 
    ELSE IF (JRR==2) THEN
      YTITLEZ   (JPROCZ) = 'Rc'
      YCOMMENTZ (JPROCZ) = '1D liquid cloud water mixing ratio' 
    ELSE IF (JRR==3) THEN
      YTITLEZ   (JPROCZ) = 'Rr'
      YCOMMENTZ (JPROCZ) = '1D Rain water mixing ratio' 
    ELSE IF (JRR==4) THEN
      YTITLEZ   (JPROCZ) = 'Ri'
      YCOMMENTZ (JPROCZ) = '1D Ice cloud water mixing ratio' 
    ELSE IF (JRR==5) THEN
      YTITLEZ   (JPROCZ) = 'Rs'
      YCOMMENTZ (JPROCZ) = '1D Snow mixing ratio' 
    ELSE IF (JRR==6) THEN
      YTITLEZ   (JPROCZ) = 'Rg'
      YCOMMENTZ (JPROCZ) = '1D Graupel mixing ratio' 
    ELSE IF (JRR==7) THEN
      YTITLEZ   (JPROCZ) = 'Rh'
      YCOMMENTZ (JPROCZ) = '1D Hail mixing ratio' 
    END IF
  END DO
!
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'FF'
  YUNITZ   (JPROCZ) = 'm s-1'
  YCOMMENTZ(JPROCZ) = 'Horizontal wind'
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%FFZ(:,IK)
!
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'IWC'
  YUNITZ   (JPROCZ) = 'kg m-3'
  YCOMMENTZ(JPROCZ) = 'Ice water content'
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%IWCZ(:,IK)
!
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'LWC'
  YUNITZ   (JPROCZ) = 'kg m-3'
  YCOMMENTZ(JPROCZ) = 'Liquid water content'
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%LWCZ(:,IK)
!
  IF (NSV_LIMA_BEG/=NSV_LIMA_END) THEN
    JPROCZ = JPROCZ + 1
    YTITLEZ  (JPROCZ) = 'CIT'
    YUNITZ   (JPROCZ) = 'm-3'
    YCOMMENTZ(JPROCZ) = 'Ice concentration'
    ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%CIZ(:,IK)
  ELSE
    JPROCZ = JPROCZ + 1
    YTITLEZ  (JPROCZ) = 'CCLOUDT'
    YUNITZ   (JPROCZ) = 'kg-1'
    YCOMMENTZ(JPROCZ) = 'liquid cloud concentration'
    ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%CCZ(:,IK)
!
    JPROCZ = JPROCZ + 1
    YTITLEZ  (JPROCZ) = 'CRAINT'
    YUNITZ   (JPROCZ) = 'kg-1'
    YCOMMENTZ(JPROCZ) = 'Rain concentration'
    ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%CRZ(:,IK)
!
    JPROCZ = JPROCZ + 1
    YTITLEZ  (JPROCZ) = 'CICET'
    YUNITZ   (JPROCZ) = 'kg-1'
    YCOMMENTZ(JPROCZ) = 'Ice concentration'
    ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%CIZ(:,IK)
 ENDIF
!
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'RARE'
  YUNITZ   (JPROCZ) = 'dBZ'
  YCOMMENTZ(JPROCZ) = '1D cloud radar reflectivity'
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%CRARE(:,IK)
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'RAREatt'
  YUNITZ   (JPROCZ) = 'dBZ'
  YCOMMENTZ(JPROCZ) = '1D cloud radar attenuated reflectivity' 
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%CRARE_ATT(:,IK)
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'W'
  YUNITZ   (JPROCZ) = 'm s-1'
  YCOMMENTZ(JPROCZ) = '1D vertical velocity' 
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%WZ(:,IK)
  JPROCZ = JPROCZ + 1
  YTITLEZ  (JPROCZ) = 'Z'
  YUNITZ   (JPROCZ) = 'm'
  YCOMMENTZ(JPROCZ) = '1D altitude above sea'
  ZWORKZ6 (1,1,IK,:,1,JPROCZ) = TPFLYER%ZZ(:,IK)
END DO
!----------------------------------------------------------------------------
!
ALLOCATE (ZW6(1,1,1,SIZE(TPFLYER%TIME),1,JPROC))
ZW6  = ZWORK6(:,:,:,:,:,:JPROC)
DEALLOCATE(ZWORK6)
ALLOCATE (ZWZ6(1,1,IKU,SIZE(TPFLYER%TIME),1,JPROCZ))
ZWZ6 = ZWORKZ6(:,:,:,:,:,:JPROCZ)
DEALLOCATE(ZWORKZ6)
!
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"RSPL",IGRID, TPFLYER%DATIME, ZW6, &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,                               &
                   PTRAJX=ZTRAJX, PTRAJY=ZTRAJY, PTRAJZ=ZTRAJZ                 )
!
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUPZ,"CART",IGRIDZ, TPFLYER%DATIME,  &
                   ZWZ6,ZTRAJT,YTITLEZ,YUNITZ,YCOMMENTZ,                     &
                   .TRUE.,.TRUE.,.FALSE.,                                    &
                   KIL=1,KIH=1,KJL=1,KJH=1,KKL=1,KKH=IKU                     )

DEALLOCATE (ZTRAJT)
DEALLOCATE (ZTRAJX)
DEALLOCATE (ZTRAJY)
DEALLOCATE (ZTRAJZ)
DEALLOCATE (ZW6)    
DEALLOCATE (YCOMMENT)
DEALLOCATE (YTITLE  )
DEALLOCATE (YUNIT   )
DEALLOCATE (IGRID   )
DEALLOCATE (ZWZ6)    
DEALLOCATE (YCOMMENTZ)
DEALLOCATE (YTITLEZ )
DEALLOCATE (YUNITZ  )
DEALLOCATE (IGRIDZ  )
!----------------------------------------------------------------------------
END SUBROUTINE FLYER_DIACHRO
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
END SUBROUTINE WRITE_AIRCRAFT_BALLOON
