!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!!    ########################
      MODULE MODI_AER_MONITOR_n
!!    ########################
!!
!
INTERFACE
!!
SUBROUTINE AER_MONITOR_n(KTCOUNT,PTSTEP, KLUOUT, KVERB, KCLD)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KTCOUNT    ! iteration count
REAL,  INTENT(IN)   :: PTSTEP     ! Double timestep 
INTEGER, INTENT(IN) :: KLUOUT     ! unit for output listing count
INTEGER, INTENT(IN) :: KVERB      ! verbosity level
LOGICAL, INTENT(IN) :: KCLD       ! conditionnal call for dust wet deposition
END SUBROUTINE AER_MONITOR_n
!!
END INTERFACE
!!
END MODULE MODI_AER_MONITOR_n
!!
!!    ####################################################### 
      SUBROUTINE AER_MONITOR_n(KTCOUNT,PTSTEP, KLUOUT, KVERB, KCLD)
!!    #######################################################
!!
!!*** *AER_MONITOR_n*  monitor of the dust sea salt module
!!
!!    PURPOSE
!!    -------
!!       The purpose of this subroutine is to control the aerosol module
!!    i.e. to pass the meteorological parameters from MesoNH to its chemical
!!    part and to call the different subroutines (calculation of rate constants,
!!    photolysis rates, stiff solver,..)
!!
!!    METHOD
!!    ------
!!       The calculation  of the aerosols terms is performed using a loop
!!    over all spatial dimensions. 
!!
!!       For each single grid point, all necessary meteorological parameters are
!!    passed into the chemical core system (variable TZM). This variable is
!!    then passed on to the subroutines that calculate the reaction and
!!    photolysis rates. Then the chemical solver is called. As the chemistry
!!    part works with different units than MesoNH (MesoNH uses mixing ratio,
!!    the chemisty part uses molec/cm3) some unit conversion is also performed.
!!
!!       Temporal integration is performed over a double timestep 2*TSTEP
!!    (except in the case of a cold start). If the timestep of MesoNH
!!    is too large for the chemical solver, several smaller steps can
!!    be taken using the NCH_SUBSTEPS parameter.
!!    One option of temporal discretization is implemented:
!!    "SPLIT"  : from XRSVS the scalar variable at t+dt is calculated and
!!               given as input to the solver; the result is rewritten 
!!               into XRSVS; this corresponds to applying first only dynamics
!!               and then only chemistry; this option assures positivity, but
!!               degrades the order of the temporal integration.
!!               In fact, an overhead of a factor two is produced here.
!!               A future solution will be to calculate the dynamics
!!               of the scalar variables not using leapfrog, but forward
!!               temporal integration.
!!
!!    REFERENCE
!!    ---------
!!    Book 1, 2, 3 of MesoNH-chemistry
!!
!!    AUTHOR
!!    ------
!!    P. Tulet from ch_monitorn.f90
!!
!!    MODIFICATIONS
!!    -------------
!
!!    Bielli S. (02/2019) : Sea salt : significant sea wave height influences salt emission; 5 salt modes
!!    EXTERNAL
!!    --------
!
USE MODD_LUNIT_n
USE MODD_NSV
USE MODD_CH_MNHC_n, ONLY : NCH_VEC_LENGTH
USE MODE_ll
USE MODE_DUST_PSD
USE MODE_SALT_PSD
USE MODE_MODELN_HANDLER
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_FIELD_n,   ONLY: XSVT,      &! scalar variable at t
                          XPABST,    &! pressure
                          XRSVS       ! source of scalar variable
!!
USE MODD_REF_n,     ONLY: XRHODREF,  &! dry density for ref. state
                          XRHODJ      ! ( rhod J ) = dry density
!!
USE MODD_PARAMETERS,ONLY: JPHEXT,    &! number of horizontal External points
                          JPVEXT      ! number of vertical External points
!!
USE MODD_CST,       ONLY: XMD,       &! Molar mass of dry air
                          XPI
                          
! parameters of the namelist to come
!
USE MODD_VAR_ll
USE MODD_DUST
USE MODD_SALT
USE MODD_FIELD_n,   ONLY: XTHT, XPABST, XRRS, XRT
USE MODD_GRID_n,    ONLY: XZZ
USE MODD_LBC_n, ONLY: CLBCX, &!X-direction LBC type at left(1)
                              ! and right(2) boundaries
                      CLBCY   ! Y-direction LBC type at left(1)
                              ! and right(2) boundaries
USE MODD_CLOUDPAR_n, ONLY: NSPLITR  ! Nb of required small time step integration
                                    ! for rain sedimentation computation
USE MODD_CONF,      ONLY: L1D, L2D, NVERB
USE MODD_CONF_n,    ONLY: LUSERC,&    ! Logical to use clouds
                          LUSERV,&    ! Logical to use wapor water
                          LUSERR,&    ! Logical to use rain water
                          NRR         ! Total number of water variables
USE MODD_PARAM_n,    ONLY: CCLOUD
USE MODD_PRECIP_n, ONLY: XEVAP3D
USE MODI_SUM_ll
USE MODI_SEDIM_DUST
USE MODI_SEDIM_SALT
USE MODI_DUST_FILTER
USE MODI_SALT_FILTER
USE MODI_AER_WET_DEP_KMT_WARM
USE MODI_MNHGET_SURF_PARAM_n

!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN) :: KTCOUNT    ! iteration count
REAL,  INTENT(IN)   ::  PTSTEP    ! Double timestep except 
                                  ! for the first time step (single one)
INTEGER, INTENT(IN) :: KLUOUT     ! unit for output listing count
INTEGER, INTENT(IN) :: KVERB      ! verbosity level
LOGICAL, INTENT(IN) :: KCLD       ! conditionnal call for dust wet deposition
!
!*      0.2    declarations of local variables
!
INTEGER :: JI,JJ,JK,JL,JM,JN,JKAQ   ! loop counters
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSVT, ZRGDST,ZSIGDST,ZSVDST,ZNDST, ZVMASSMIN
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZRGSLT,ZSIGSLT,ZSVSLT,ZNSLT
REAL, DIMENSION(:,:),     ALLOCATABLE  :: ZSEA, ZTOWN
REAL, DIMENSION(:,:,:),   ALLOCATABLE  :: ZRCS, ZRRS
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZDENSITY
REAL, DIMENSION(:),       ALLOCATABLE  :: ZMASSMIN, ZINIRADIUS
REAL    :: ZSIGMIN 
INTEGER :: IMODEIDX
!

INTEGER             :: IIU  ! Upper dimension in x direction
INTEGER             :: IJU  ! Upper dimension in y direction
INTEGER             :: IKU  ! Upper dimension in z direction
INTEGER             :: IIB  ! indice I Beginning in x direction
INTEGER             :: IJB  ! indice J Beginning in y direction
INTEGER             :: IKB  ! indice K Beginning in z direction
INTEGER             :: IIE  ! indice I End       in x direction
INTEGER             :: IJE  ! indice J End       in y direction
INTEGER             :: IKE  ! indice K End       in z direction
INTEGER             :: JSV  ! loop index for SV
INTEGER             :: IMI  ! model index
!-------------------------------------------------------------------------------
!
!
!*       1.    Prepare monitor
!              ---------------
!
!*       1.1   compute dimensions of arrays
!
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IMI = GET_CURRENT_MODEL_INDEX()
!
IKU = SIZE(XRSVS,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
!*       1.2   calculate timestep variables
!
! ++ JORIS DEBUG ++
IF (NVERB == 10) WRITE(*,*) 'dans aer_monitorn.f90 1.'
! -- JORIS DEBUG --
!
! ++ PIERRE / MARINE SSA DUST - MODIF ++
!  XRSVS(:,:,:,NSV_DSTBEG:NSV_DSTEND) = &
!                      MAX(XRSVS(:,:,:,NSV_DSTBEG:NSV_DSTEND), 0.)  
! -- PIERRE / MARINE SSA DUST - MODIF --
!
!*       2.   Sedimentation of aerosols 
!              ------------------------
!*       2.1   Sedimentation of  dusts
IF (LDUST.AND.LSEDIMDUST) THEN
!
  ALLOCATE(ZSVT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NSV_DSTEND-NSV_DSTBEG+1))
  DO JSV = NSV_DSTBEG, NSV_DSTEND
    ZSVT(:,:,:,JSV-NSV_DSTBEG+1) = XRSVS(:,:,:,JSV) * PTSTEP / XRHODJ(:,:,:)
  ENDDO
! ++ PIERRE / MARINE SSA DUST - MODIF ++
 CALL DUST_FILTER(ZSVT,XRHODREF)
!  CALL DUST_FILTER(ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,:),&
!                   XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE)) 
! -- PIERRE / MARINE SSA DUST - MODIF --
  CALL SEDIM_DUST(XTHT(IIB:IIE,IJB:IJE,IKB:IKE), PTSTEP,&
                  XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE), &
                  XPABST(IIB:IIE,IJB:IJE,IKB:IKE), &
                  XZZ(IIB:IIE,IJB:IJE,IKB:IKE+1),    &
                  ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,:)) !ppp (concentration)
!
DO JSV = NSV_DSTBEG, NSV_DSTEND
    XRSVS(IIB:IIE,IJB:IJE,IKB:IKE,JSV) = ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,JSV-NSV_DSTBEG+1)  *&
                                         XRHODJ(IIB:IIE,IJB:IJE,IKB:IKE) /  PTSTEP
END DO
!
  DEALLOCATE(ZSVT)
END IF 
!
!*       2.1   Sedimentation of  Sea salt
IF ((LSALT).AND.(LSEDIMSALT)) THEN
!
  ALLOCATE(ZSVT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NSV_SLTEND-NSV_SLTBEG+1))
  DO JSV = NSV_SLTBEG, NSV_SLTEND
    ZSVT(:,:,:,JSV-NSV_SLTBEG+1) = XRSVS(:,:,:,JSV) * PTSTEP / XRHODJ(:,:,:)
  ENDDO

! ++ JORIS DEBUG ++
  CALL SALT_FILTER(ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                   XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE))
! 
  CALL SEDIM_SALT(XTHT(IIB:IIE,IJB:IJE,IKB:IKE),PTSTEP,&
                  XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE), &
                  XPABST(IIB:IIE,IJB:IJE,IKB:IKE), &
                  XZZ(IIB:IIE,IJB:IJE,IKB:IKE+1),    &
                  ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,:))  !ppp (concentration)
! -- JORIS DEBUG --

DO JSV = NSV_SLTBEG, NSV_SLTEND
 XRSVS(IIB:IIE,IJB:IJE,IKB:IKE,JSV) = &
             ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,JSV-NSV_SLTBEG+1)  *&
             XRHODJ(IIB:IIE,IJB:IJE,IKB:IKE) /  PTSTEP
END DO
!
  DEALLOCATE(ZSVT)
END IF 

IF (LDUST .AND. LDEPOS_DST(IMI) .AND. KCLD) THEN
!-------------------------------------------------------------------------------
!*       3.    Dust / Cloud / Rain interactions 
!              ------------------------------------
!        


ALLOCATE(ZSIGDST(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_DST))
ALLOCATE(ZRGDST(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_DST))  
ALLOCATE(ZSVDST(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_DST*3))
ALLOCATE(ZNDST(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_DST))
ALLOCATE(ZSEA(SIZE(XSVT,1),SIZE(XSVT,2)))
ALLOCATE(ZTOWN(SIZE(XSVT,1),SIZE(XSVT,2)))
ALLOCATE(ZMASSMIN(NMODE_DST))
ALLOCATE(ZINIRADIUS(NMODE_DST))
ALLOCATE(ZVMASSMIN(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_DST))
ALLOCATE(ZSVT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),SIZE(XSVT,4)))
ALLOCATE(ZRCS(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3)))
ALLOCATE(ZRRS(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3)))
ALLOCATE(ZDENSITY(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_DST))
!

ZSVDST(:,:,:,:) = 0.
ZSVT(:,:,:,:) = 0.
ZSEA(:,:) = 0.
ZTOWN(:,:) = 0.
ZRCS(:,:,:) = XRRS(:,:,:,2) * PTSTEP / XRHODJ(:,:,:) 
ZRRS(:,:,:) = XRRS(:,:,:,3) * PTSTEP / XRHODJ(:,:,:)
ZDENSITY(:,:,:,:) = XDENSITY_DUST

!
!     3.1 Minimum mass to transfer between dry mass or in-cloud droplets

DO JN=1,NMODE_DST
  IMODEIDX = JPDUSTORDER(JN)
   IF (CRGUNITD=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX) * EXP(-3.*(LOG(XINISIG(IMODEIDX)))**2)
   ELSE
    ZINIRADIUS(JN) = XINIRADIUS(IMODEIDX)
   END IF
   IF (LVARSIG) THEN
    ZSIGMIN = XSIGMIN
   ELSE
    ZSIGMIN = XINISIG(IMODEIDX)
   ENDIF
   ZMASSMIN(JN) = XN0MIN(IMODEIDX) * (ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(ZSIGMIN)**2)
! volume/um3 =>  #/molec_{air}
   ZVMASSMIN(:,:,:,JN)=  ZMASSMIN(JN) * XMD * XPI * 4./3. * XDENSITY_DUST  / &
                 (XMOLARWEIGHT_DUST*XM3TOUM3*XRHODREF(:,:,:))
ENDDO
  
!

!     3.2 Derive moment from aerosol moments sources
!         from moments

  DO JSV=1,SIZE(XRSVS,4)
    ZSVT(:,:,:,JSV) = XRSVS(:,:,:,JSV) * PTSTEP / XRHODJ(:,:,:) 
  ENDDO    

!     3.3 Compute and store Standard deviation and mean radius 
!         from moments
  CALL PPP2DUST(ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_DSTBEG:NSV_DSTEND), &
                XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),                   &
                PSIG3D=ZSIGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),           &
                PRG3D=ZRGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),             &
                PN3D=ZNDST(IIB:IIE,IJB:IJE,IKB:IKE,:))

!     3.4 Compute acquous aerosol mass vector from moment scalar vector
!
  DO JSV= 1, NMODE_DST
   IF (LVARSIG) THEN
    ZSVDST(:,:,:,JSV) = ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*3)
   ELSE IF (LRGFIX_DST) THEN
    ZSVDST(:,:,:,JSV) = ZSVT(:,:,:,NSV_DSTBEG+JSV-1)
   ELSE
    ZSVDST(:,:,:,JSV) = ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*2)
   ENDIF
  ENDDO    
  DO JSV=1,NSV_DSTDEP
    ZSVDST(:,:,:,NMODE_DST+JSV) = ZSVT(:,:,:,NSV_DSTDEPBEG-1+JSV)
  ENDDO
  ZSVDST(:,:,:,:) = MAX(ZSVDST(:,:,:,:), 0.)
  ZSVT(:,:,:,:) = MAX(ZSVT(:,:,:,:), 0.)

!     3.5 Mass transfer between dry mass and in-cloud mass aerosols
SELECT CASE (CCLOUD)
  
  CASE ('KESS','REVE','ICE3','ICE4')
! One moment cloud scheme
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP,                     &
                              XZZ(IIB:IIE,IJB:IJE,IKB:IKE),        &
                              XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),   &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,2),      &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,3),      &
                              ZRCS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZRRS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZSVDST(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XTHT(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              XPABST(IIB:IIE,IJB:IJE,IKB:IKE),     &
                              ZRGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XEVAP3D(IIB:IIE,IJB:IJE,IKB:IKE),    &
                              NMODE_DST,                           &
                              ZDENSITY(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                              ZVMASSMIN(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                              PSEA=ZSEA(IIB:IIE,IJB:IJE),          &
                              PTOWN=ZTOWN(IIB:IIE,IJB:IJE))
  CASE ('KHKO','C2R2','C3R5')
! Two moment cloud scheme
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP,                     &
                              XZZ(IIB:IIE,IJB:IJE,IKB:IKE),        &
                              XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),   &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,2),      &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,3),      &
                              ZRCS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZRRS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZSVDST(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XTHT(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              XPABST(IIB:IIE,IJB:IJE,IKB:IKE),     &
                              ZRGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XEVAP3D(IIB:IIE,IJB:IJE,IKB:IKE),    &
                              NMODE_DST,                           &
                              ZDENSITY(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                              ZVMASSMIN(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                              PCCT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_C2R2BEG+1),&
                              PCRT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_C2R2BEG+2) )
!++th++ 05/05/17 ajout LIMA
CASE ('LIMA')
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP,                     &
                              XZZ(IIB:IIE,IJB:IJE,IKB:IKE),        &
                              XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),   &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,2),      &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,3),      &
                              ZRCS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZRRS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZSVDST(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XTHT(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              XPABST(IIB:IIE,IJB:IJE,IKB:IKE),     &
                              ZRGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XEVAP3D(IIB:IIE,IJB:IJE,IKB:IKE),    &
                              NMODE_DST,                           &
                              ZDENSITY(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                              ZVMASSMIN(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                              PCCT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_LIMA_NC),&
                              PCRT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_LIMA_NR) )
!--th--
END SELECT

!     3.5 Compute return to moment vector
  DO JSV=1,NMODE_DST
   IF (LVARSIG) THEN
    ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*3) = ZSVDST(:,:,:,JSV)
   ELSE IF (LRGFIX_DST) THEN
    ZSVT(:,:,:,NSV_DSTBEG+JSV-1) = ZSVDST(:,:,:,JSV)
   ELSE
    ZSVT(:,:,:,NSV_DSTBEG+1+(JSV-1)*2) = ZSVDST(:,:,:,JSV)
   ENDIF
  ENDDO    
                             
!     3.5 Return to lognormal distribution (compute M0 and M6 using RG, SIG and
!     new mass from M3)
  CALL DUST2PPP(ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_DSTBEG:NSV_DSTEND), &
                XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),                   &
                ZSIGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),                  &
                ZRGDST(IIB:IIE,IJB:IJE,IKB:IKE,:))
!  
!     3.6 Return to source term
  DO JSV=NSV_DSTBEG,NSV_DSTEND
     XRSVS(IIB:IIE,IJB:IJE,IKB:IKE,JSV) =  ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,JSV) * &
                                           XRHODJ(IIB:IIE,IJB:IJE,IKB:IKE) / PTSTEP
  ENDDO
  DO JSV=1,NSV_DSTDEP
    XRSVS(IIB:IIE,IJB:IJE,IKB:IKE,NSV_DSTDEPBEG-1+JSV)=ZSVDST(IIB:IIE,IJB:IJE,IKB:IKE,NMODE_DST+JSV) *&
                                              XRHODJ(IIB:IIE,IJB:IJE,IKB:IKE) / PTSTEP
  ENDDO

  DEALLOCATE(ZSIGDST)
  DEALLOCATE(ZRGDST)  
  DEALLOCATE(ZSVDST)
  DEALLOCATE(ZNDST)
  DEALLOCATE(ZSEA)
  DEALLOCATE(ZTOWN)
  DEALLOCATE(ZMASSMIN)
  DEALLOCATE(ZINIRADIUS)
  DEALLOCATE(ZVMASSMIN)
  DEALLOCATE(ZSVT)
  DEALLOCATE(ZRCS)
  DEALLOCATE(ZRRS)
  DEALLOCATE(ZDENSITY)
END IF
IF (LSALT .AND. LDEPOS_SLT(IMI) .AND. KCLD) THEN

!-------------------------------------------------------------------------------

!*       4.    Sea Salt / Cloud / Rain interactions 
!              ------------------------------------
ALLOCATE(ZSIGSLT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_SLT))
ALLOCATE(ZRGSLT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_SLT))  
ALLOCATE(ZSVSLT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_SLT*3))
ALLOCATE(ZNSLT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_SLT))
ALLOCATE(ZSEA(SIZE(XSVT,1),SIZE(XSVT,2)))
ALLOCATE(ZTOWN(SIZE(XSVT,1),SIZE(XSVT,2)))
ALLOCATE(ZMASSMIN(NMODE_SLT))
ALLOCATE(ZINIRADIUS(NMODE_SLT))
ALLOCATE(ZVMASSMIN(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_SLT))
ALLOCATE(ZSVT(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),SIZE(XSVT,4)))
ALLOCATE(ZRCS(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3)))
ALLOCATE(ZRRS(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3)))
ALLOCATE(ZDENSITY(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),NMODE_SLT))
!

ZSVSLT(:,:,:,:) = 0.
ZSVT(:,:,:,:) = 0.
ZSEA(:,:) = 0.
ZTOWN(:,:) = 0.
ZRCS(:,:,:) = XRRS(:,:,:,2) * PTSTEP / XRHODJ(:,:,:) 
ZRRS(:,:,:) = XRRS(:,:,:,3) * PTSTEP / XRHODJ(:,:,:)
ZDENSITY(:,:,:,:) = XDENSITY_SALT

!
!     4.1 Minimum mass to transfer between dry mass or in-cloud droplets
! ++ PIERRE / MARINE SSA DUST - MODIF ++
DO JN=1,NMODE_SLT
  IMODEIDX = JPSALTORDER(JN)
   IF (CRGUNITD=="MASS") THEN
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX) * EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
   ELSE
    ZINIRADIUS(JN) = XINIRADIUS_SLT(IMODEIDX)
   END IF
   IF (LVARSIG) THEN
    ZSIGMIN = XSIGMIN_SLT
   ELSE
    ZSIGMIN = XINISIG_SLT(IMODEIDX)
   ENDIF
   ZMASSMIN(JN) = XN0MIN_SLT(IMODEIDX) * (ZINIRADIUS(JN)**3)*EXP(4.5 * LOG(ZSIGMIN)**2)
! volume/um3 =>  #/molec_{air}
   ZVMASSMIN(:,:,:,JN)=  ZMASSMIN(JN) * XMD * XPI * 4./3. * XDENSITY_SALT  / &
                 (XMOLARWEIGHT_SALT*XM3TOUM3_SALT*XRHODREF(:,:,:))
ENDDO
! -- PIERRE / MARINE SSA DUST - MODIF --
!

!     4.2 Derive moment from aerosol moments sources
!         from moments
  XRSVS(:,:,:,NSV_SLTDEPBEG:NSV_SLTDEPEND) = &
                      MAX(XRSVS(:,:,:,NSV_SLTDEPBEG:NSV_SLTDEPEND), 0.)

  DO JSV=1,SIZE(XRSVS,4)
    ZSVT(:,:,:,JSV) = XRSVS(:,:,:,JSV) * PTSTEP / XRHODJ(:,:,:) 
  ENDDO    

!     4.3 Compute and store Standard deviation and mean radius 
!         from moments
  CALL PPP2SALT(ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_SLTBEG:NSV_SLTEND), &
                XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),                   &
                PSIG3D=ZSIGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),           &
                PRG3D=ZRGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),             &
                PN3D=ZNSLT(IIB:IIE,IJB:IJE,IKB:IKE,:))

!     4.4 Compute acquous aerosol mass vector from moment scalar vector
!
  DO JSV= 1, NMODE_SLT
   IF (LVARSIG) THEN
    ZSVSLT(:,:,:,JSV) = ZSVT(:,:,:,NSV_SLTBEG+1+(JSV-1)*3)
   ELSE IF (LRGFIX_SLT) THEN
    ZSVSLT(:,:,:,JSV) = ZSVT(:,:,:,NSV_SLTBEG+JSV-1)
   ELSE
    ZSVSLT(:,:,:,JSV) = ZSVT(:,:,:,NSV_SLTBEG+1+(JSV-1)*2)
   ENDIF
  ENDDO    
  DO JSV=1,NSV_SLTDEP
    ZSVSLT(:,:,:,NMODE_SLT+JSV) = XRSVS(:,:,:,NSV_SLTDEPBEG-1+JSV) *&
                                            PTSTEP / XRHODJ(:,:,:)
  ENDDO

!     4.5 Mass transfer between dry mass and in-cloud mass aerosols
SELECT CASE (CCLOUD)
  
  CASE ('KESS','REVE','ICE3','ICE4')
! One moment cloud scheme
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP,                    &
                              XZZ(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),  &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,2),     &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,3),     &
                              ZRCS(IIB:IIE,IJB:IJE,IKB:IKE),      &
                              ZRRS(IIB:IIE,IJB:IJE,IKB:IKE),      &
                              ZSVSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),  &
                              XTHT(IIB:IIE,IJB:IJE,IKB:IKE),      &
                              XPABST(IIB:IIE,IJB:IJE,IKB:IKE),    &
                              ZRGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),  &
                              XEVAP3D(IIB:IIE,IJB:IJE,IKB:IKE),   &
                              NMODE_SLT,                           &
                              ZDENSITY(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                              ZVMASSMIN(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                              PSEA=ZSEA(IIB:IIE,IJB:IJE),          &
                              PTOWN=ZTOWN(IIB:IIE,IJB:IJE))
  CASE ('KHKO','C2R2','C3R5')
! Two moment cloud scheme
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP,                     &
                              XZZ(IIB:IIE,IJB:IJE,IKB:IKE),        &
                              XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),   &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,2),      &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,3),      &
                              ZRCS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZRRS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZSVSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XTHT(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              XPABST(IIB:IIE,IJB:IJE,IKB:IKE),     &
                              ZRGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XEVAP3D(IIB:IIE,IJB:IJE,IKB:IKE),    &
                              NMODE_SLT,                           &
                              ZDENSITY(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                              ZVMASSMIN(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                              PCCT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_C2R2BEG+1),&
                              PCRT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_C2R2BEG+2) )
!++th++05/05/17 ajout LIMA
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP,                     &
                              XZZ(IIB:IIE,IJB:IJE,IKB:IKE),        &
                              XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),   &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,2),      &
                              XRT(IIB:IIE,IJB:IJE,IKB:IKE,3),      &
                              ZRCS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZRRS(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              ZSVSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XTHT(IIB:IIE,IJB:IJE,IKB:IKE),       &
                              XPABST(IIB:IIE,IJB:IJE,IKB:IKE),     &
                              ZRGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
                              XEVAP3D(IIB:IIE,IJB:IJE,IKB:IKE),    &
                              NMODE_SLT,                           &
                              ZDENSITY(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                              ZVMASSMIN(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                              PCCT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_LIMA_NC),&
                              PCRT=ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_LIMA_NR) )
!--th--
END SELECT

!     4.5 Compute return to moment vector
  DO JSV=1,NMODE_SLT
   IF (LVARSIG) THEN
    ZSVT(:,:,:,NSV_SLTBEG+1+(JSV-1)*3) = ZSVSLT(:,:,:,JSV)
   ELSE IF (LRGFIX_SLT) THEN
    ZSVT(:,:,:,NSV_SLTBEG+JSV-1) = ZSVSLT(:,:,:,JSV)
   ELSE
    ZSVT(:,:,:,NSV_SLTBEG+1+(JSV-1)*2) = ZSVSLT(:,:,:,JSV)
   ENDIF
  ENDDO    
                             
!     4.5 Return to lognormal distribution (compute M0 and M6 using RG, SIG and
!     new mass from M3)
  CALL SALT2PPP(ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_SLTBEG:NSV_SLTEND), &
                XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),                   &
                ZSIGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),                  &
                ZRGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:))
!  
!     4.6 Return to source term
  DO JSV=NSV_SLTBEG,NSV_SLTEND
     XRSVS(:,:,:,JSV) =  ZSVT(:,:,:,JSV) * XRHODJ(:,:,:) / PTSTEP
  ENDDO
  DO JSV=1,NSV_SLTDEP
     XRSVS(:,:,:,NSV_SLTDEPBEG-1+JSV) =  ZSVSLT(:,:,:,NMODE_SLT+JSV) *&
                                              XRHODJ(:,:,:) / PTSTEP
  ENDDO

  DEALLOCATE(ZSIGSLT)
  DEALLOCATE(ZRGSLT)  
  DEALLOCATE(ZSVSLT)
  DEALLOCATE(ZNSLT)
  DEALLOCATE(ZSEA)
  DEALLOCATE(ZTOWN)
  DEALLOCATE(ZMASSMIN)
  DEALLOCATE(ZINIRADIUS)
  DEALLOCATE(ZVMASSMIN)
  DEALLOCATE(ZSVT)
  DEALLOCATE(ZRCS)
  DEALLOCATE(ZRRS)
  DEALLOCATE(ZDENSITY)
END IF

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AER_MONITOR_n
