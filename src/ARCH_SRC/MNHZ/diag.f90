!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG3 2007/11/20 11:57:53
!-----------------------------------------------------------------
!     ############
      PROGRAM DIAG
!     ############
!
!!****
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!          MODD_DYN
!!          MODD_CONF
!!          MODD_PARAMETERS
!!          MODD_CONF_n
!!          MODD_DYN_n
!!          MODD_DIM_n
!!          MODD_ADV_n
!!          MODD_FIELD_n
!!          MODD_GRID_n
!!          MODD_LBC_n
!!          MODD_PARAM_n
!!          MODD_REF_n
!!          MODD_LUNIT_n
!!          MODD_OUT_n
!!          MODD_TIME_n
!!          MODD_TURB_n
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      E. Richard                  * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!  15/03/99    (V. Masson)    call to PHYS_PARAM1 and new PGD fields
!!  08/06/01    (F. Meleux)    call to diagnostics for chemistry
!!  13/07/01    (J Stein)      trajectory computation
!!  15/10/01    (I.Mallet)     allow namelists in different orders
!!  05/02/02    (G.Jaubert)    aircraft and balloon computation
!!  15/03/02    (F Solmon)     replace NMODE by NRAD_DIAG
!!  29/11/02    (J.-P. Pinty)  add C3R5, ICE2, ICE4, ELEC
!!  15/04/03    (J.-P. Chaboureau) add LRAD_SUBG_COND
!!  01/2004     (Masson)       surface externalization
!!  19/03/2008  (J.Escobar)    rename INIT to INIT_MNH --> grib problem
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_VAR_ll
!
USE MODD_DYN
USE MODD_CONF
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_LUNIT
USE MODD_PARAMETERS
!
USE MODD_DIAG_FLAG
USE MODD_STO_FILE
!
USE MODD_CONF_n
USE MODD_DIM_n
USE MODD_DYN_n
USE MODD_ADV_n
USE MODD_DEEP_CONVECTION_n
USE MODD_FIELD_n
USE MODD_GRID_n
USE MODD_LBC_n
USE MODD_PARAM_n
USE MODD_PARAM_CONVECT_n
USE MODD_PARAM_RAD_n
USE MODD_RADIATIONS_n
USE MODD_TURB_n
USE MODD_REF_n
USE MODD_LUNIT_n
USE MODD_OUT_n
USE MODD_TIME_n
USE MODD_NESTING, ONLY: CDAD_NAME
USE MODD_NSV
USE MODD_AIRCRAFT_BALLOON
USE MODD_TIME
USE MODD_GR_FIELD_n
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_BUDGET
!
USE MODN_DIAG_BLANK
!
USE MODI_PHYS_PARAM_n
USE MODI_WRITE_LFIFM1_FOR_DIAG
USE MODI_WRITE_LFIFM1_FOR_DIAG_SUPP
USE MODI_COMPUTE_R00
USE MODI_WRITE_LFIFMN_FORDIACHRO_n
USE MODI_AIRCRAFT_BALLOON
USE MODI_WRITE_AIRCRAFT_BALLOON
!
USE MODE_POS
USE MODE_TIME
USE MODE_FM
USE MODE_IO_ll
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODD_AIRCRAFT_BALLOON
USE MODD_PROFILER_n
USE MODD_STATION_n
!
!JUANZ
USE MODE_MNH_TIMING
!JUANZ
!     
IMPLICIT NONE
!
!*       0.1   declarations of local variables
!
TYPE(DATE_TIME)   :: TXDTBAL   ! current time and date for BALLOON and AIRCRAFT trajectories
CHARACTER (LEN=28), DIMENSION(1) :: YINIFILE ! names of the INPUT FM-file
CHARACTER (LEN=28) :: YFMFILE   ! name of the OUTPUT FM-file
CHARACTER (LEN=28) :: YFMDIAC   ! name of the OUTPUT Diachronic file
CHARACTER (LEN=5)  :: YSUFFIX   ! character string for the OUTPUT FM-file number
!CHARACTER (LEN=32) :: YDESFM    ! name of the desfm part of this FM-file
CHARACTER (LEN=4)  :: YRAD      ! initial flag to call to radiation schemes
CHARACTER (LEN=4)  :: YDCONV    ! initial flag to call to convection schemes
CHARACTER (LEN=4)  :: YTURB     ! initial flag to call to turbulence schemes
CHARACTER (LEN=40) :: YFMT,YFMT2! format for cpu analysis printing
INTEGER  :: IRESP               ! return code in FM routines
INTEGER  :: ILUOUT0             ! Logical unit number for the output listing
REAL,DIMENSION(2)     :: ZTIME0,ZTIME1,ZTIME2,ZRAD,ZDCONV,ZGROUND,ZTURB,ZCHEM,ZTIME_BU ! CPU time 
REAL,DIMENSION(2)     :: ZSTART,ZINIT,ZWRIT,ZBALL,ZPHYS,ZSURF,ZWRITS,ZTRAJ ! storing variables
INTEGER :: INPRAR               ! number of articles predicted  in
                                !  the LFIFM file
INTEGER :: ITYPE=1              ! type of file (conv2dia and transfer)
LOGICAL :: GCLOSE_OUT = .FALSE. ! conditional closure of the OUTPUT FM-file
INTEGER :: ISTEPBAL   ! loop indice for balloons and aircraft
INTEGER :: ILUNAM      ! Logical unit numbers for the namelist file
                       ! and for output_listing file
CHARACTER (LEN=9) :: YNAM ! name of the namelist file
INTEGER        :: JF =0   !  loop index
LOGICAL :: GFOUND         ! Return code when searching namelist
LOGICAL, DIMENSION(:,:),ALLOCATABLE     :: GMASKkids ! kids domains mask
!
INTEGER :: IIU, IJU, IKU
INTEGER :: IINFO_ll               ! return code for _ll routines 
!
NAMELIST/NAM_DIAG/ CISO, LVAR_RS, LVAR_LS,   &
                   NCONV_KF, NRAD_3D, CRAD_SAT, NRTTOVINFO, LRAD_SUBG_COND,  &
                   LVAR_TURB,LTURBFLX,LTURBDIAG,XDTSTEP,  &
                   LVAR_MRW, LVAR_MRSV, LVAR_FRC, &
                   LTPZH, LMOIST_V, LMOIST_E, LCOREF, &
                   LVORT, LDIV, LMEAN_POVO, XMEAN_POVO, &
                   LGEO, LAGEO, LWIND_ZM, LMSLP, LTHW, &
                   LCLD_COV, LVAR_PR, LTOTAL_PR, LMEAN_PR, XMEAN_PR, &
                   NCAPE, LBV_FR, LRADAR, LBLTOP, LTRAJ, &
                   LDIAG,XDIAG,LCHEMDIAG,XCHEMLAT,XCHEMLON, &
                   LAIRCRAFT_BALLOON,NTIME_AIRCRAFT_BALLOON,&
                   XSTEP_AIRCRAFT_BALLOON,&
                   XLAT_BALLOON,XLON_BALLOON,XALT_BALLOON,&
                   LC2R2, LC3R5, LELECDIAG
!
NAMELIST/NAM_DIAG_FILE/ YINIFILE, YSUFFIX
NAMELIST/NAM_STO_FILE/ CFILES, NSTART_SUPP
!
!-------------------------------------------------------------------------------
!
!*       0.0   Initializations
!              ---------------
!
CALL SECOND_MNH2(ZTIME1)
ZTIME0=ZTIME1
!
CALL GOTO_MODEL(1)
!
CALL VERSION
CPROGRAM='DIAG  '
!
CALL INITIO_ll()
!
! initialization of logical for the diagnostics
!
CISO='PREVTK'
LVAR_RS=.TRUE.
LVAR_LS=.FALSE.
NCONV_KF=-1
NRAD_3D=-1
CRAD_SAT='                                            '
LRAD_SUBG_COND=.TRUE.
NRTTOVINFO(:,:)=NUNDEF
LVAR_TURB=.FALSE.
LTURBFLX=.FALSE.
LTURBDIAG=.FALSE.
XDTSTEP=XUNDEF
LVAR_MRW=.FALSE.
LVAR_MRSV=.FALSE.
LTPZH=.FALSE.
LMOIST_V=.FALSE.
LMOIST_E=.FALSE.
LCOREF=.FALSE.
LVORT=.FALSE.
LDIV=.FALSE.
LMEAN_POVO=.FALSE.
XMEAN_POVO(1)=15000
XMEAN_POVO(2)=50000
LGEO=.FALSE.
LAGEO=.FALSE.
LWIND_ZM=.FALSE.
LMSLP=.FALSE.
LTHW=.FALSE.
LCLD_COV=.FALSE.
LVAR_PR=.FALSE.
LTOTAL_PR=.FALSE.
LMEAN_PR=.FALSE.
XMEAN_PR(1:2)=1.
NCAPE=-1
LBV_FR=.FALSE.
LRADAR=.FALSE.
LBLTOP=.FALSE.
LVAR_FRC=.FALSE.
LCHEMDIAG=.FALSE.
XCHEMLAT(:)=XUNDEF
XCHEMLON(:)=XUNDEF
LTRAJ=.FALSE.
!
LAIRCRAFT_BALLOON=.FALSE.
NTIME_AIRCRAFT_BALLOON=NUNDEF
XSTEP_AIRCRAFT_BALLOON=XUNDEF
XLAT_BALLOON(:)=XUNDEF
XLON_BALLOON(:)=XUNDEF
XALT_BALLOON(:)=XUNDEF
!
LDIAG(:)=.FALSE.
XDIAG(:)=XUNDEF
!
YINIFILE(:) = '                         '
YSUFFIX='_DIAG'
!
CFILES(:) = '                         '
NSTART_SUPP(:) = NUNDEF
CFILES_STA(:) = 'INIT_SV'
!
!-------------------------------------------------------------------------------
!
!*       1.0   Namelist reading
!              ----------------
!
PRINT*, ' '
PRINT*, '*********************************************************************'
PRINT*, '*********************************************************************'
PRINT*, ' '
!
YNAM  = 'DIAG1.nam'
CALL OPEN_ll (UNIT=ILUNAM,FILE=YNAM,IOSTAT=IRESP,STATUS="OLD",ACTION='READ', &
     FORM="FORMATTED",POSITION="REWIND",MODE=GLOBAL)
!
PRINT*, 'READ THE DIAG.NAM FILE'
!
CALL POSNAM(ILUNAM,'NAM_DIAG',GFOUND)
IF (GFOUND) THEN
  READ(UNIT=ILUNAM,NML=NAM_DIAG)
  PRINT*, '  namelist NAM_DIAG read'
END IF
!
CALL POSNAM(ILUNAM,'NAM_DIAG_BLANK',GFOUND)
IF (GFOUND) THEN
  READ(UNIT=ILUNAM,NML=NAM_DIAG_BLANK)
  PRINT*, '  namelist NAM_DIAG_BLANK read'
END IF
!
CALL POSNAM(ILUNAM,'NAM_DIAG_FILE',GFOUND)
IF (GFOUND) THEN
  READ(UNIT=ILUNAM,NML=NAM_DIAG_FILE)
  PRINT*, '  namelist NAM_DIAG_FILE read'
END IF
!
CALL POSNAM(ILUNAM,'NAM_STO_FILE',GFOUND)
IF (GFOUND) THEN
  READ(UNIT=ILUNAM,NML=NAM_STO_FILE)
  PRINT*, '  namelist NAM_STO_FILE read'
END IF
!
CALL CLOSE_ll(YNAM)
!
CINIFILE = YINIFILE(1)
!
IF (LTRAJ) THEN
  JF=1
  DO WHILE (LEN_TRIM(CFILES(JF))/=0)
    JF=JF+1
  END DO
  !
  IF (JF/=1) THEN
    IF (CINIFILE==CFILES(JF-1)) THEN
      PRINT*,'*****************************'
      PRINT*,'* INITIAL FILE NON TREATED  *'
      PRINT*,'*****************************'
      STOP
    END IF
  END IF
!
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       2.0   MESONH file
!              -----------
!
IF ( LEN_TRIM(CINIFILE)==0 ) STOP
!
PRINT*, ' '
PRINT*, '****************************************'
PRINT*, 'Treatment of file: ',CINIFILE
PRINT*, '****************************************'
!
INPRAR = 24 +2*(4+NRR+NSV)
YFMFILE=ADJUSTL(ADJUSTR(CINIFILE)//YSUFFIX)
CALL FMOPEN_ll(YFMFILE,'WRITE',CLUOUT,INPRAR,ITYPE,NVERB,INPRAR,IRESP)
COUTFMFILE=YFMFILE
! YDESFM=ADJUSTL(ADJUSTR(YFMFILE)//'.des')
! CALL WRITE_DESFM1(1,YDESFM,CLUOUT)
!
CALL SECOND_MNH2(ZTIME2)
ZSTART=ZTIME2-ZTIME1
ZTIME1=ZTIME2
!-------------------------------------------------------------------------------
!
!*       3.0   Fields initialization
!
CALL INIT_MNH
PRINT*, ' '
PRINT*, 'DIAG AFTER INIT'
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=NKMAX+2*JPVEXT
!
!* allocation of variables used 
!
ALLOCATE(GMASKkids  (IIU,IJU))
GMASKkids(:,:)=.FALSE.
! 
CALL INI_DIAG_IN_RUN(IIU,IJU,IKU,LFLYER,LSTATION,LPROFILER)
!
CALL SECOND_MNH2(ZTIME2)
ZINIT =ZTIME2-ZTIME1
ZTIME1=ZTIME2
!-------------------------------------------------------------------------------
!
!*       4.0    Stores the fields in MESONH files if necessary
!
CALL WRITE_LFIFM1_FOR_DIAG(YFMFILE,CDAD_NAME(1))
PRINT*, ' '
PRINT*, 'DIAG AFTER WRITE_LFIFM1_FOR_DIAG'
PRINT*, ' '
!
CALL SECOND_MNH2(ZTIME2)
ZWRIT =ZTIME2-ZTIME1
ZTIME1=ZTIME2
!-------------------------------------------------------------------------------
!
!*       4.1    BALLOON and AIRCRAFT
!
IF ( LAIRCRAFT_BALLOON ) THEN
  YFMDIAC=ADJUSTL(ADJUSTR(CINIFILE)//'BAL')
  CALL FMOPEN_ll(YFMDIAC,'WRITE',CLUOUT,INPRAR,ITYPE,NVERB,INPRAR,IRESP)
  PRINT*, ' '
  PRINT*, 'DIAG AFTER OPEN DIACHRONIC FILE'
  PRINT*, ' '
  CALL WRITE_LFIFMN_FORDIACHRO_n(YFMDIAC)
!
  TXDTBAL%TDATE%YEAR = TDTCUR%TDATE%YEAR
  TXDTBAL%TDATE%MONTH = TDTCUR%TDATE%MONTH
  TXDTBAL%TDATE%DAY = TDTCUR%TDATE%DAY
  TXDTBAL%TIME = TDTCUR%TIME - NTIME_AIRCRAFT_BALLOON/2.
  CALL SUBTRACT_TO_DATE(TXDTBAL%TDATE%YEAR, &
                          TXDTBAL%TDATE%MONTH,&
                          TXDTBAL%TDATE%DAY,  &
                          TXDTBAL%TIME        )
!
  DO ISTEPBAL=1,NTIME_AIRCRAFT_BALLOON,INT(XSTEP_AIRCRAFT_BALLOON)
!
    CALL AIRCRAFT_BALLOON(CLUOUT,XSTEP_AIRCRAFT_BALLOON, &
                      TDTEXP, TDTMOD, TDTCUR, TXDTBAL, &
                      XXHAT, XYHAT, XZZ, XMAP, XLONORI, XLATORI, &
                      XUT, XVT, XWT, XPABST, XTHT, XRT, XSVT, XTKET, XTSRAD)
!
!-----------------------------
!
    TXDTBAL%TIME=TXDTBAL%TIME + XSTEP_AIRCRAFT_BALLOON
    CALL ADD_FORECAST_TO_DATE(TXDTBAL%TDATE%YEAR, &
                          TXDTBAL%TDATE%MONTH,&
                          TXDTBAL%TDATE%DAY,  &
                          TXDTBAL%TIME        )
  ENDDO
  CALL WRITE_AIRCRAFT_BALLOON(YFMDIAC)
  CALL MENU_DIACHRO(YFMDIAC,CLUOUT,'END')
  CALL FMCLOS_ll(YFMDIAC,'KEEP',CLUOUT,IRESP)
  PRINT*, ' '
  PRINT*, 'DIAG AFTER CLOSE DIACHRONIC FILE'
  PRINT*, ' '
END IF
!
CALL SECOND_MNH2(ZTIME2)
ZBALL =ZTIME2-ZTIME1
ZTIME1=ZTIME2
!-------------------------------------------------------------------------------
!
!*       5.0   Call to physics
!
!* initialise the source terms
!
XRUS (:,:,:) = 0.
XRVS (:,:,:) = 0.
XRWS (:,:,:) = 0.
XRTHS(:,:,:) = 0.
IF (NRR>=1)             XRRS  (:,:,:,:) = 0.
IF (NSV>=1)             XRSVS (:,:,:,:) = 0.
IF (CTURB /= 'NONE')    XRTKES(:,:,:) = 0.
!
!* stores the initial flags

!
YTURB  = CTURB
YDCONV = CDCONV
YRAD   = CRAD
!
!* turbulence scheme
!
LTURB_DIAG=LTURBDIAG
LTURB_FLX =LTURBFLX
! no need to recompute the turbulent tendencies.
IF ( .NOT. LTURB_FLX .AND. .NOT. LTURB_DIAG ) THEN
  CTURB  = 'NONE'
END IF
! no way to compute the turbulent tendencies.
IF ( ( LTURB_FLX .OR. LTURB_DIAG ) .AND. CSTORAGE_TYPE/='MT' ) THEN
  CTURB  = 'NONE'
  PRINT*, '******************* WARNING in DIAG ***********************'
  PRINT*, ' '
  PRINT*, 'You wanted to compute turbulence fluxes or diagnostics,'
  PRINT*, 'But the initial file comes from PREP_REAL_CASE.'
  PRINT*, 'Therefore, the boundary layer turbulence is meaningless.'
  PRINT*, 'Turbulence fluxes and diagnostics will NOT be computed'
  PRINT*, 'Please make your turbulence diagnostics from a meso-NH file'
  PRINT*, 'coming from a MESO-NH simulation.'
END IF
!
!* convective scheme
!
IF (NCONV_KF == -1) CDCONV = 'NONE'
!
IF (NCONV_KF >= 0 ) THEN
  CALL  SM_PRINT_TIME(TDTCUR,CLUOUT0,'CURRENT TIME ')
  CALL  SM_PRINT_TIME(TDTDCONV,CLUOUT0,'LAST CONVECTION CALL')
  CDCONV='KAFR'
  LDIAGCONV= .TRUE.
  TDTDCONV=TDTCUR
END IF
!
!* radiation scheme
!
IF (NRAD_3D == -1) CRAD = 'NONE'
!
IF (NRAD_3D >= 0) THEN
  IF (YRAD=='ECMW') THEN
    ! radiative fields are already initialized by INIT
    CRAD = 'NONE'
  ELSE
    CRAD = 'ECMW'
  ENDIF
  IF (NRAD_3D >= 1) THEN
    NRAD_DIAG = NRAD_3D
    CRAD = 'ECMW'    ! radiation scheme is called to compute extra diags
  END IF
END IF
!
IF (LEN_TRIM(CRAD_SAT) /= 0 .AND. YRAD/='ECMW') THEN
  CRAD = 'ECMW'
END IF
!
!
IF ( CTURB /= 'NONE' .OR. CDCONV /= 'NONE' .OR. CRAD /= 'NONE' ) THEN
  IF (CSTORAGE_TYPE/='MT') THEN
    IF (XDTSTEP==XUNDEF) THEN
      PRINT*, ' '
      PRINT*, '******************* WARNING in DIAG ***********************'
      PRINT*, ' '
      PRINT*, 'You asked for diagnostics that need to call the physics monitor:'
      PRINT*, ' be aware of the time step used'
      PRINT*, 'you can modify it with XDTSTEP in namelist NAM_DIAG'
      PRINT*, ' '
    ELSE
      XTSTEP=XDTSTEP
    END IF
  END IF
  PRINT*,' XTSTEP= ', XTSTEP
  PRINT*, ' '
  PRINT*, 'DIAG BEFORE PHYS_PARAM1: CTURB=',CTURB,' CDCONV=',CDCONV, &
          ' CRAD=',CRAD
END IF
!
!* call to physics monitor
!
GCLOSE_OUT=.TRUE.
ZRAD=0.
ZDCONV=0.
ZGROUND=0.
ZTURB=0.
ZCHEM=0.
XTIME_LES=0.
XTIME_LES_BU_PROCESS=0.
XTIME_BU_PROCESS=0.
!
CALL PHYS_PARAM_n(1,XTSTEP,YFMFILE,GCLOSE_OUT,ZRAD,ZDCONV,ZGROUND,ZTURB,ZCHEM, &
                 ZTIME_BU,GMASKkids)
PRINT*, 'DIAG AFTER PHYS_PARAM1'
!
!* restores the initial flags
!
CTURB  = YTURB
CDCONV = YDCONV
CRAD   = YRAD
!
CALL SECOND_MNH2(ZTIME2)
ZPHYS =ZTIME2-ZTIME1
ZTIME1=ZTIME2
!-------------------------------------------------------------------------------
!
!*       6.0    Surface diagnostics
!
IF (CSURF=='EXTE') THEN
  CALL GOTO_SURFEX(1)
  CALL WRITE_SURF_ATM_n('MESONH','ALL')
  CALL DIAG_SURF_ATM_n('MESONH')
  CALL WRITE_DIAG_SURF_ATM_n('MESONH','ALL')
  PRINT*, ' '
  PRINT*, 'DIAG AFTER WRITE_DIAG_SURF_ATM_n'
ENDIF
!
CALL SECOND_MNH2(ZTIME2)
ZSURF =ZTIME2-ZTIME1
ZTIME1=ZTIME2
!
!-------------------------------------------------------------------------------
!
!*       7.0    Stores other fields in MESONH files if necessary
!
CALL WRITE_LFIFM1_FOR_DIAG_SUPP(YFMFILE)
PRINT*, ' '
PRINT*, 'DIAG AFTER WRITE_LFIFM1_FOR_DIAG_SUPP'
!
CALL SECOND_MNH2(ZTIME2)
ZWRITS=ZTIME2-ZTIME1
ZTIME1=ZTIME2
!-------------------------------------------------------------------------------
!
!*       8.0    Initial positions computation (back into simulation segments)
!
IF (LTRAJ .AND. JF/=1) CALL COMPUTE_R00(YFMFILE)
!
CALL SECOND_MNH2(ZTIME2)
ZTRAJ =ZTIME2-ZTIME1
ZTIME1=ZTIME2
!-------------------------------------------------------------------------------
!
!*       9.0    Closes the FM files
!
IF (GCLOSE_OUT) THEN
  GCLOSE_OUT=.FALSE.
  CALL FMCLOS_ll(YFMFILE,'KEEP',CLUOUT,IRESP)
END IF
!
CALL CLOSE_ll (CLUOUT,IOSTAT=IRESP)
!
CALL SECOND_MNH2(ZTIME2)
ZTIME2=ZTIME2-ZTIME0
!-------------------------------------------------------------------------------
!
CALL FMLOOK_ll(CLUOUT0,CLUOUT0,ILUOUT0,IRESP)
!
WRITE(ILUOUT0,*) '+--------------------------------------------------------------+'
WRITE(ILUOUT0,*) '|                                                              |'
WRITE(ILUOUT0,*) '|            COMPUTING TIME ANALYSIS in DIAG                   |'
WRITE(ILUOUT0,*) '|                                                              |'
WRITE(ILUOUT0,*) '|--------------------------------------------------------------|'
WRITE(ILUOUT0,*) '|                     |                    |                   |'
WRITE(ILUOUT0,*) '|    ROUTINE NAME     |      CPU-TIME      |   PERCENTAGE %    |'
WRITE(ILUOUT0,*) '|                     |                    |                   |'
WRITE(ILUOUT0,*) '|---------------------| -------------------|-------------------|'
WRITE(ILUOUT0,*) '|                     |                    |                   |'
YFMT='(A,F9.3,A,F9.3,A)'
YFMT2='(A,A4,A,F9.3)'
WRITE(ILUOUT0,YFMT) '|        START        |     ',ZSTART,'      |     ',100.*ZSTART/ZTIME2,'     |'
WRITE(ILUOUT0,YFMT) '|        INIT         |     ',ZINIT,'      |     ',100.*ZINIT/ZTIME2,'     |'
WRITE(ILUOUT0,YFMT) '|        WRIT         |     ',ZWRIT,'      |     ',100.*ZWRIT/ZTIME2,'     |'
WRITE(ILUOUT0,YFMT) '|        BALL         |     ',ZBALL,'      |     ',100.*ZBALL/ZTIME2,'     |'
WRITE(ILUOUT0,YFMT) '|        PHYS         |     ',ZPHYS,'      |     ',100.*ZPHYS/ZTIME2,'     |'
IF (ZRAD(1)>0.) &
  WRITE(ILUOUT0,YFMT2) '|          ',CRAD,'       |     ',ZRAD
IF (ZDCONV(1)>0.) &
  WRITE(ILUOUT0,YFMT2) '|          ',CDCONV,'       |     ',ZDCONV
IF (ZGROUND(1)>0.) &
  WRITE(ILUOUT0,YFMT2) '|          ',CSURF,'       |     ',ZGROUND
IF (ZTURB(1)>0.) &
  WRITE(ILUOUT0,YFMT2) '|          ',CTURB,'       |     ',ZTURB
IF (ZCHEM(1)>0.) &
  WRITE(ILUOUT0,'(A,F9.3)') '|       CHEM          |     ',ZCHEM
WRITE(ILUOUT0,YFMT) '|        SURF         |     ',ZSURF,'      |     ',100.*ZSURF/ZTIME2,'     |'
WRITE(ILUOUT0,YFMT) '|        WRITS        |     ',ZWRITS,'      |     ',100.*ZWRITS/ZTIME2,'     |'
WRITE(ILUOUT0,YFMT) '|        TRAJ         |     ',ZTRAJ,'      |     ',100.*ZTRAJ/ZTIME2,'     |'
WRITE(ILUOUT0,*) '|                     |                    |                   |'
WRITE(ILUOUT0,*) '|---------------------| -------------------|-------------------|'
!
!
CALL CLOSE_ll (CLUOUT0,IOSTAT=IRESP)
PRINT*, ' '
PRINT*, '***************************** **************'
PRINT*, '*            EXIT  DIAG CORRECTLY          *'
PRINT*, '**************************** ***************'
PRINT*, '  (see time analysis in ',TRIM(CLUOUT0),' )'
PRINT*, ' '
!-------------------------------------------------------------------------------
!
!*      10.    FINALIZE THE PARALLEL SESSION
!              -----------------------------
!
CALL END_PARA_ll(IINFO_ll)
!
!-------------------------------------------------------------------------------
END PROGRAM DIAG
