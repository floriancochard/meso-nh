!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_ONE_WAY_n
!     ###################
!
INTERFACE 
!
      SUBROUTINE ONE_WAY_n(KDAD,PTSTEP,KMI,KTCOUNT,                      &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,     &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,     &
                    KDXRATIO,KDYRATIO,KDTRATIO,                          &
                    HLBCX,HLBCY,KRIMX,KRIMY,                             &
                    KKLIN_LBXU,PCOEFLIN_LBXU,KKLIN_LBYU,PCOEFLIN_LBYU,   &
                    KKLIN_LBXV,PCOEFLIN_LBXV,KKLIN_LBYV,PCOEFLIN_LBYV,   &
                    KKLIN_LBXW,PCOEFLIN_LBXW,KKLIN_LBYW,PCOEFLIN_LBYW,   &
                    KKLIN_LBXM,PCOEFLIN_LBXM,KKLIN_LBYM,PCOEFLIN_LBYM,   &
                    OSTEADY_DMASS,HCLOUD,OUSECHAQ,OUSECHIC,              &
                    PLBXUM,PLBYUM,PLBXVM,PLBYVM,PLBXWM,PLBYWM,           &
                    PLBXTHM,PLBYTHM,                                     &
                    PLBXTKEM,PLBYTKEM,                                   &
                    PLBXRM,PLBYRM,PLBXSVM,PLBYSVM,                       &
                    PDRYMASST,PDRYMASSS,                                 &
                    PLBXUS,PLBYUS,PLBXVS,PLBYVS,PLBXWS,PLBYWS,           &
                    PLBXTHS,PLBYTHS,                                     &
                    PLBXTKES,PLBYTKES,                                   &
                    PLBXRS,PLBYRS,PLBXSVS,PLBYSVS                        )
!
!
INTEGER,          INTENT(IN)    :: KDAD     !  Number of the DAD model
REAL,             INTENT(IN)    :: PTSTEP   !  Time step
INTEGER,          INTENT(IN)    :: KMI      ! model number

INTEGER,          INTENT(IN)    :: KTCOUNT  !  Temporal loop COUNTer
                                            ! (=1 at the segment beginning)
!
                                    ! interpolation coefficients 
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
INTEGER,   INTENT(IN)  :: KDTRATIO   !  Time step resolution RATIO
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
INTEGER,          INTENT(IN)    :: KRIMX,KRIMY ! size of the RIM area
!  coefficients for the vertical interpolation of the LB fields
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXU,KKLIN_LBYU 
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXU,PCOEFLIN_LBYU
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXV,KKLIN_LBYV 
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXV,PCOEFLIN_LBYV
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXW,KKLIN_LBYW 
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXW,PCOEFLIN_LBYW
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXM,KKLIN_LBYM 
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXM,PCOEFLIN_LBYM
!
LOGICAL,           INTENT(IN)  :: OSTEADY_DMASS ! Md evolution logical switch
CHARACTER (LEN=4), INTENT(IN)  :: HCLOUD        ! Indicator of the cloud scheme
LOGICAL,           INTENT(IN)  :: OUSECHAQ      ! logical for aqueous phase chemistry
LOGICAL,           INTENT(IN)  :: OUSECHIC      ! logical for ice phase chemistry
!  
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLBXUM,PLBXVM,PLBXWM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLBYUM,PLBYVM,PLBYWM 
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLBXTHM ,PLBYTHM  ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLBXTKEM,PLBYTKEM ! Theta, TKE
REAL, DIMENSION(:,:,:,:),INTENT(IN)  :: PLBXRM  ,PLBYRM   ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(IN)  :: PLBXSVM ,PLBYSVM  ! in x and y-dir.
!
REAL,             INTENT(INOUT) :: PDRYMASST     ! Mass of dry air Md 
REAL,             INTENT(INOUT) :: PDRYMASSS     !  Md source 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PLBXUS,PLBXVS,PLBXWS ! Large Scale source terms 
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PLBYUS,PLBYVS,PLBYWS
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLBXTHS ,PLBYTHS  ! Large Scale fields sources
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLBXTKES,PLBYTKES ! Theta, TKE
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PLBXRS  ,PLBYRS   ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PLBXSVS ,PLBYSVS  ! in x and y-dir.
!
END SUBROUTINE ONE_WAY_n
!
END INTERFACE
!
END MODULE MODI_ONE_WAY_n
!
!     ####################################################################
SUBROUTINE ONE_WAY_n(KDAD,PTSTEP,KMI,KTCOUNT,                            &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,     &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,     &
                    KDXRATIO,KDYRATIO,KDTRATIO,                          &
                    HLBCX,HLBCY,KRIMX,KRIMY,                             &
                    KKLIN_LBXU,PCOEFLIN_LBXU,KKLIN_LBYU,PCOEFLIN_LBYU,   &
                    KKLIN_LBXV,PCOEFLIN_LBXV,KKLIN_LBYV,PCOEFLIN_LBYV,   &
                    KKLIN_LBXW,PCOEFLIN_LBXW,KKLIN_LBYW,PCOEFLIN_LBYW,   &
                    KKLIN_LBXM,PCOEFLIN_LBXM,KKLIN_LBYM,PCOEFLIN_LBYM,   &
                    OSTEADY_DMASS,HCLOUD,OUSECHAQ,OUSECHIC,              &
                    PLBXUM,PLBYUM,PLBXVM,PLBYVM,PLBXWM,PLBYWM,           &
                    PLBXTHM,PLBYTHM,                                     &
                    PLBXTKEM,PLBYTKEM,                                   &
                    PLBXRM,PLBYRM,PLBXSVM,PLBYSVM,                       &
                    PDRYMASST,PDRYMASSS,                                 &
                    PLBXUS,PLBYUS,PLBXVS,PLBYVS,PLBXWS,PLBYWS,           &
                    PLBXTHS,PLBYTHS,                                     &
                    PLBXTKES,PLBYTKES,                                   &
                    PLBXRS,PLBYRS,PLBXSVS,PLBYSVS                        )
!     ####################################################################
!
!!****  *ONE_WAY_n* - Refreshing of a nested model Large Scale sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of ONE_WAY$n is to 'refresh' Large Scale sources
!!    of all the prognostic variables of the current nested model when the
!!    current time step is in phase with its outer (DAD) model $n.
!!      It also computes the dry mass at time t and the corresponding source,
!!    by integration of the dry density of outer model $n over the inner domain.
!
!
!!**  METHOD
!!    ------
!!      The basic task consists in interpolating fields from outer model $n
!!    to present inner model, using horizontal Bikhardt interpolation.
!!      The dry density reads:
!!                                          P
!!                       Rhod = --------------------------
!!                               Rd*Theta*PI*(1+rv*Rv/Rd)
!!      The total dry mass is deduced from Rhod and the Jacobian (the integration
!!    is performed on the inner points):
!!
!!                         Md =  SUM   rhod J
!!                              i,j,k
!!      Caution:  J is deduced from RHODJ and RHODREF depending on the system of
!!      -------  equations (CEQNSYS) with:
!!               RHODJ =  RHODeff * J
!!                    and RHODeff = RHODREF                  if CEQNSYS = MAE or LHE
!!                                           THVREF*(1+RVREF)
!!                     or RHODeff = RHODREF* --------------- if CEQNSYS = DUR
!!                                               TH00
!!
!!    EXTERNAL
!!    --------
!!
!!        Function  VER_INTERP_LIN : performs the vertical interpolation
!!
!!        Subroutine  BIKHARDT : performs the horizontal interpolation
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: JPHEXT,JPVEXT
!!
!!      Module MODD_CST: XRD,XRV,XCPD,XP00,XTH00
!!
!!      Module MODD_CONF: CEQNSYS
!!
!!      Module MODD_FIELD$n : XUT,XVT,XWT,XRT,XTHT,XPABST
!!
!!      Module MODD_REF$n   : XRHODJ, XRVREF,XTHVREF, XRHODREF
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    J. P. Lafore  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     22/10/96
!!    J. P. Lafore   20/01/97 nesting procedure for DRYMASS
!!    J. P. Lafore   21/10/97 DRYMASS correction to account for the
!!                            equations system (CEQNSYS)
!!    J. Stein       22/12/97 use the LB fields for lbc
!!    J. Stein       08/04/99 merge uvw_ls_nesting and scalar_ls_nesting
!!    P. Jabouille   19/04/00 parallelisation (without use of LB comlib routines)
!!    J.-P. Pinty    02/11/00 modify the LB*SVS for the C2R2 scheme
!!    J.-P. Pinty    29/11/02 modify the LB*SVS for the C3R5 scheme
!!                            and add ICE2, ICE4, CELEC
!!    O.Geoffroy     03/2006  Add KHKO scheme
!!                   05/2006  Remove EPS
!!    M. Leriche     11/2009  modify the LB*SVS for the aqueous phase chemistry
!!                   07/2010  idem for ice phase chemical species
!!    Bosseur & Filippi 07/2013 Adds Forefire
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      J.Escobar : 18/12/2015 : Correction of bug in bound in // for NHALO <>1 
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODD_PARAMETERS
USE MODD_NESTING
USE MODD_CST
USE MODD_REF_n        ! modules relative to the outer model $n
USE MODD_FIELD_n
USE MODD_CONF
USE MODD_PARAM_n
USE MODD_CH_MNHC_n, ONLY: LUSECHAQ, LUSECHIC
USE MODD_NSV
!
USE MODI_BIKHARDT
USE MODI_VER_INTERP_LIN
USE MODI_SET_CONC_RAIN_C2R2
USE MODI_SET_CONC_ICE_C1R3
USE MODI_SET_CHEMAQ_1WAY
USE MODI_SET_CONC_LIMA
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
INTEGER,          INTENT(IN)    :: KDAD     !  Number of the DAD model
REAL,             INTENT(IN)    :: PTSTEP   !  Time step
INTEGER,          INTENT(IN)    :: KMI      ! model number
INTEGER,          INTENT(IN)    :: KTCOUNT  !  Temporal loop COUNTer
                                            ! (=1 at the segment beginning)
                                    ! interpolation coefficients
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
INTEGER,   INTENT(IN)  :: KDTRATIO   !  Time step resolution RATIO
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
INTEGER,          INTENT(IN)    :: KRIMX,KRIMY ! size of the RIM area
!  coefficients for the vertical interpolation of the LB fields
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXU,KKLIN_LBYU
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXU,PCOEFLIN_LBYU
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXV,KKLIN_LBYV
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXV,PCOEFLIN_LBYV
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXW,KKLIN_LBYW
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXW,PCOEFLIN_LBYW
INTEGER, DIMENSION(:,:,:), INTENT(  IN ) :: KKLIN_LBXM,KKLIN_LBYM
REAL,    DIMENSION(:,:,:), INTENT(  IN ) :: PCOEFLIN_LBXM,PCOEFLIN_LBYM
!
LOGICAL,           INTENT(IN)  :: OSTEADY_DMASS ! Md evolution logical switch
CHARACTER (LEN=4), INTENT(IN)  :: HCLOUD        ! Indicator of the cloud scheme
LOGICAL,           INTENT(IN)  :: OUSECHAQ      ! logical for aqueous phase chemistry
LOGICAL,           INTENT(IN)  :: OUSECHIC      ! logical for ice phase chemistry
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLBXUM,PLBXVM,PLBXWM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLBYUM,PLBYVM,PLBYWM
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLBXTHM ,PLBYTHM  ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLBXTKEM,PLBYTKEM ! Theta, TKE
REAL, DIMENSION(:,:,:,:),INTENT(IN)  :: PLBXRM  ,PLBYRM   ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(IN)  :: PLBXSVM ,PLBYSVM  ! in x and y-dir.
!
REAL,             INTENT(INOUT) :: PDRYMASST     ! Mass of dry air Md
REAL,             INTENT(INOUT) :: PDRYMASSS     !  Md source
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PLBXUS,PLBXVS,PLBXWS ! Large Scale source terms
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PLBYUS,PLBYVS,PLBYWS
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLBXTHS ,PLBYTHS  ! Large Scale fields sources
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLBXTKES,PLBYTKES ! Theta, TKE
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PLBXRS  ,PLBYRS   ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(OUT) :: PLBXSVS ,PLBYSVS  ! in x and y-dir.
!
!
!*       0.2   declarations of local variables
!
REAL                   :: ZTIME                   ! Interpolation length
INTEGER                :: IIB,IIE,IJB,IJE,IIU,IJU
INTEGER                :: ILBX,ILBY,ILBX2,ILBY2
REAL     ::   ZBIGTSTEP    ! time step of the dad model ($n)
REAL     ::   ZRV_O_RD     ! = Rv /  Rd
REAL     ::   ZRD_O_CPD    ! = Rd /  Cpd
REAL     ::  ZDRYMASST,ZDRYMASSM
!REAL,   DIMENSION(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3)) :: ZJ,ZRHOD
REAL,   DIMENSION(:,:,:), ALLOCATABLE :: ZJ,ZRHOD
REAL,   DIMENSION(:,:,:), ALLOCATABLE  :: ZWORK
LOGICAL  :: GVERT_INTERP
!
INTEGER           :: IRR,ISV_USER          !  Number of moist and scalar variables
INTEGER           :: JRR,JSV          !  Loop index
!
! reduced array for the interpolation coefficients
REAL, DIMENSION(:,:,:), ALLOCATABLE ::  ZCOEFLIN_LBXM_RED,ZCOEFLIN_LBYM_RED
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IKLIN_LBXM_RED,IKLIN_LBYM_RED
!
INTEGER :: IINFO_ll, IDIMX, IDIMY
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTUT, ZTVT, ZTWT, ZTTHT, ZTTKET
REAL, DIMENSION(:,:,:,:), ALLOCATABLE ::ZTRT,ZTSVT
!
CHARACTER(LEN=4)                    :: ZINIT_TYPE
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCONCT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCHEMT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCHEMTI
!
INTEGER :: N_VAR_LIMA
!
!-------------------------------------------------------------------------------
!
!*      0.   INITIALISATION
!            --------------
CALL GOTO_MODEL(KDAD)
ALLOCATE(ZJ(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3)))
ALLOCATE(ZRHOD(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3)))
!
!!$CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
IIB=1
IIE=IIU
IJB=1
IJE=IJU

ALLOCATE(ZWORK(IIB:IIE,IJB:IJE,SIZE(PLBXTHM,3)))  ! can be smaller than child extended subdomain
! LS_FORCING routine can not correctly manage extra halo zone
! LB will be filled only with one layer halo zone for the moment
!
!
ZRV_O_RD  = XRV / XRD
ZRD_O_CPD = XRD / XCPD
!
ZTIME = PTSTEP * (1+KDTRATIO)
GVERT_INTERP=.TRUE.
ZJ(:,:,:)   =0.
ZRHOD(:,:,:)=0.
!
IRR=MIN(SIZE(XRT,4),SIZE(PLBXRM,4))
ISV_USER=MIN(NSV_USER_A(KDAD),NSV_USER_A(KMI))
!
IF (LWEST_ll() .AND. LEAST_ll()) THEN
  ALLOCATE (ZCOEFLIN_LBXM_RED(2,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
  ALLOCATE (  IKLIN_LBXM_RED(2,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
ELSE
  ALLOCATE (ZCOEFLIN_LBXM_RED(1,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
  ALLOCATE (  IKLIN_LBXM_RED(1,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
ENDIF
!
IF (LSOUTH_ll() .AND. LNORTH_ll()) THEN
  ALLOCATE (ZCOEFLIN_LBYM_RED(SIZE(PLBYTHM,1),2,SIZE(PLBYTHM,3)))
  ALLOCATE (  IKLIN_LBYM_RED(SIZE(PLBYTHM,1),2,SIZE(PLBYTHM,3)))
ELSE
  ALLOCATE (ZCOEFLIN_LBYM_RED(SIZE(PLBYTHM,1),1,SIZE(PLBYTHM,3)))
  ALLOCATE (  IKLIN_LBYM_RED(SIZE(PLBYTHM,1),1,SIZE(PLBYTHM,3)))
ENDIF
!

IF(LWEST_ll()) THEN
  ZCOEFLIN_LBXM_RED(1,:,:)=PCOEFLIN_LBXM(1,:,:)
  IKLIN_LBXM_RED(1,:,:)=KKLIN_LBXM(1,:,:)
ENDIF
IF(LEAST_ll()) THEN
  ZCOEFLIN_LBXM_RED(SIZE(ZCOEFLIN_LBXM_RED,1),:,:) = &
             PCOEFLIN_LBXM(SIZE(PCOEFLIN_LBXM,1),:,:)

  IKLIN_LBXM_RED(SIZE(IKLIN_LBXM_RED,1),:,:) = &
             KKLIN_LBXM(SIZE(IKLIN_LBXM_RED,1),:,:)
ENDIF
IF ( SIZE(PLBYTHM,2) /= 0 ) THEN
  IF(LSOUTH_ll()) THEN
    ZCOEFLIN_LBYM_RED(:,1,:)=PCOEFLIN_LBYM(:,1,:)
    IKLIN_LBYM_RED(:,1,:)=KKLIN_LBYM(:,1,:)
  ENDIF
  IF(LNORTH_ll()) THEN
    ZCOEFLIN_LBYM_RED(:,SIZE(ZCOEFLIN_LBYM_RED,2),:) = &
               PCOEFLIN_LBYM(:,SIZE(PCOEFLIN_LBYM,2),:)
    IKLIN_LBYM_RED(:,SIZE(IKLIN_LBYM_RED,2),:) = &
               KKLIN_LBYM(:,SIZE(IKLIN_LBYM_RED,2),:)
  ENDIF
END IF
!
!
!
!*      1 GATHER LB FIELD FOR THE CHILD MODEL KMI
!
!       1.1  Must be on the father model to call get_child_dim
!
CALL GO_TOMODEL_ll(KDAD, IINFO_ll)
CALL GET_CHILD_DIM_ll(KMI, IDIMX, IDIMY, IINFO_ll)
!
!
!-------------------------------------------------------------------------------
!
ALLOCATE(ZTUT(IDIMX,IDIMY,SIZE(XUT,3)))
ZTUT(:,:,:)=0

ALLOCATE(ZTVT(IDIMX,IDIMY,SIZE(XVT,3)))
ZTVT(:,:,:)=0

ALLOCATE(ZTWT(IDIMX,IDIMY,SIZE(XWT,3)))
ZTWT(:,:,:)=0

ALLOCATE(ZTTHT(IDIMX,IDIMY,SIZE(XTHT,3)))
ZTTHT(:,:,:)=0
IF (SIZE(XTKET) /= 0) ALLOCATE(ZTTKET(IDIMX,IDIMY,SIZE(XTKET,3)))
IF (IRR /= 0) ALLOCATE(ZTRT(IDIMX,IDIMY,SIZE(XRT,3),IRR))
IF (NSV_A(KMI)/= 0) ALLOCATE(ZTSVT(IDIMX,IDIMY,SIZE(XSVT,3),NSV_A(KMI)))
!
!         1.3  Specify the ls "source" fields and receiver fields
!
CALL SET_LSFIELD_1WAY_ll(XUT,ZTUT,KMI)
CALL SET_LSFIELD_1WAY_ll(XVT,ZTVT,KMI)
CALL SET_LSFIELD_1WAY_ll(XWT,ZTWT,KMI)
CALL SET_LSFIELD_1WAY_ll(XTHT,ZTTHT,KMI)
IF (ALLOCATED(ZTTKET)) CALL SET_LSFIELD_1WAY_ll(XTKET,ZTTKET,KMI)
!
DO JRR=1,IRR
  CALL SET_LSFIELD_1WAY_ll(XRT(:,:,:,JRR),ZTRT(:,:,:,JRR),KMI)
ENDDO
!
IF (ALLOCATED(ZTSVT)) ZTSVT=0.  ! to treat ISV_USER+1:NSV_USER_A(KMI) section
! USERs scalar variables
DO JSV=1,ISV_USER
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV),ZTSVT(:,:,:,JSV),KMI)
ENDDO
!
!  Checking if it is necessary to compute the Nc and Nr
!  concentrations to use the C2R2(or KHKO) microphysical scheme
!  (FATHER does not use C2R2(or KHKO) and CHILD uses C2R2(or KHKO))
!
IF (HCLOUD=="C2R2" .OR. HCLOUD=="KHKO") THEN
  IF  (CCLOUD/="NONE" .AND. CCLOUD/="C2R2" .AND. CCLOUD/="KHKO") THEN
    ZINIT_TYPE="NONE"
    ALLOCATE(ZCONCT(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),3))
    IF (CCLOUD == "REVE") THEN
      ZINIT_TYPE = "INI1"
    ELSE IF (CCLOUD == "KESS" ) THEN
      ZINIT_TYPE = "INI2"
    END IF
    CALL SET_CONC_RAIN_C2R2(ZINIT_TYPE,XRHODREF,XRT,ZCONCT)
    DO JSV=1,3
      CALL SET_LSFIELD_1WAY_ll(ZCONCT(:,:,:,JSV),&
           &ZTSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
    ENDDO
  ELSE
    DO JSV=1,NSV_C2R2_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KDAD)),&
           &ZTSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
    END DO
  ENDIF
ENDIF
!
!  Checking also if it is necessary to compute the Ni
!  concentrations to use the C3R5 microphysical scheme
!  (FATHER does not use C3R5 and CHILD uses C3R5)
!
IF (HCLOUD=="C3R5") THEN
  IF ( CCLOUD(1:3)=="ICE" ) THEN
    ZINIT_TYPE="NONE"
    ALLOCATE(ZCONCT(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),5))
    IF (CCLOUD == "REVE") THEN
      ZINIT_TYPE = "INI1"
    ELSE IF (CCLOUD == "KESS" ) THEN
      ZINIT_TYPE = "INI2"
    END IF
    CALL SET_CONC_RAIN_C2R2(ZINIT_TYPE,XRHODREF,XRT,ZCONCT)
    DO JSV=1,3
      CALL SET_LSFIELD_1WAY_ll(ZCONCT(:,:,:,JSV),&
           &ZTSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
    ENDDO
    ZINIT_TYPE="INI3"
    CALL SET_CONC_ICE_C1R3 (XRHODREF,XRT,ZCONCT)
    DO JSV=4,5
      CALL SET_LSFIELD_1WAY_ll(ZCONCT(:,:,:,JSV),&
           &ZTSVT(:,:,:,JSV-4+NSV_C1R3BEG_A(KMI)),KMI) 
    ENDDO
  ELSE
    DO JSV=1,NSV_C2R2_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KDAD)),&
           &ZTSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
    END DO
    DO JSV=1,NSV_C1R3_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C1R3BEG_A(KDAD)),&
           &ZTSVT(:,:,:,JSV-1+NSV_C1R3BEG_A(KMI)),KMI)
    END DO
  ENDIF
ENDIF
!
! LIMA Scheme
!
IF (HCLOUD=="LIMA"  ) THEN
   IF (CCLOUD/="LIMA") THEN
      ALLOCATE(ZCONCT(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),NSV_LIMA_A(KMI)))
      IF (CCLOUD == "REVE") THEN
         ZINIT_TYPE = "INI1"
      ELSE
         ZINIT_TYPE = "NONE"
      END IF
      CALL SET_CONC_LIMA (ZINIT_TYPE,XRHODREF,XRT,ZCONCT)
      DO JSV=1,NSV_LIMA_A(KMI)
         CALL SET_LSFIELD_1WAY_ll(ZCONCT(:,:,:,JSV),&
              &ZTSVT(:,:,:,JSV-1+NSV_LIMA_BEG_A(KMI)),KMI)
      ENDDO   
   ELSE
      IF (NSV_LIMA_A(KMI)/=NSV_LIMA_A(KDAD)) CALL ABORT
      DO JSV=1,NSV_LIMA_A(KMI)
         CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LIMA_BEG_A(KDAD)),&
              &ZTSVT(:,:,:,JSV-1+NSV_LIMA_BEG_A(KMI)),KMI)
      END DO
   END IF
ENDIF
!
! electrical variables
!
DO JSV=1,MIN(NSV_ELEC_A(KMI),NSV_ELEC_A(KDAD))
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_ELECBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_ELECBEG_A(KMI)),KMI)
END DO
!
! chemical Scalar variables
!  Checking if it is necessary to compute the Caq
!  concentrations to use the aqueous phase chemistry
!  (FATHER does not use aqueous phase chemistry and CHILD uses it)
!
IF (OUSECHAQ) THEN
  IF (.NOT.(LUSECHAQ)) THEN
    ALLOCATE(ZCHEMT(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),&
                   NSV_CHEM_A(KMI)))
    CALL SET_CHEMAQ_1WAY(XRHODREF,&
         XSVT(:,:,:,NSV_CHEMBEG_A(KDAD):NSV_CHEMEND_A(KDAD)),ZCHEMT)
    DO JSV=1,NSV_CHEM_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(ZCHEMT(:,:,:,JSV),&
         &ZTSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
    ENDDO
  ELSE
    DO JSV=1,NSV_CHEM_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KDAD)),&
           &ZTSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
    END DO
  ENDIF
ELSE
!
  DO JSV=1,NSV_CHEM_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
  END DO
ENDIF
!
!  Checking if it is necessary to compute the Cic
!  concentrations to use the ice phase chemistry
!  (FATHER does not use ice phase chemistry and CHILD uses it)
!
IF (OUSECHIC) THEN
  IF (.NOT.(LUSECHIC)) THEN
    ALLOCATE(ZCHEMTI(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),&
                   NSV_CHIC_A(KMI)))
    ZCHEMTI(:,:,:,:) = 0.
    DO JSV=1,NSV_CHIC_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(ZCHEMTI(:,:,:,JSV),&
       &ZTSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
    ENDDO
  ELSE
    DO JSV=1,NSV_CHIC_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KDAD)),&
           &ZTSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
    END DO
  ENDIF
ELSE
  DO JSV=1,NSV_CHIC_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KDAD)),&
         &ZTSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
  END DO
ENDIF
!
!
! dust Scalar variables
!
DO JSV=1,NSV_DST_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_DSTBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_DSTBEG_A(KMI)),KMI)
END DO
!
!
! sea salt Scalar variables
!
DO JSV=1,NSV_SLT_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_SLTBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_SLTBEG_A(KMI)),KMI)
END DO
!
! Orilam aerosol Scalar variables
!
DO JSV=1,NSV_AER_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_AERBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_AERBEG_A(KMI)),KMI)
END DO
!
! lagrangian variables
!
DO JSV=1,NSV_LG_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LGBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_LGBEG_A(KMI)),KMI)
END DO
!
! Passive pollutants    
!
DO JSV=1,NSV_PP_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_PPBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_PPBEG_A(KMI)),KMI)
END DO
#ifdef MNH_FOREFIRE

!
! ForeFire variables     
!
DO JSV=1,NSV_FF_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_FFBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_FFBEG_A(KMI)),KMI)
END DO
#endif
!
! Conditional sampling  
!
DO JSV=1,NSV_CS_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CSBEG_A(KDAD)),&
       &ZTSVT(:,:,:,JSV-1+NSV_CSBEG_A(KMI)),KMI)
END DO
!
! Communication
!
CALL LS_FORCING_ll(KMI,IINFO_ll)
CALL GO_TOMODEL_ll(KMI,IINFO_ll)
CALL UNSET_LSFIELD_1WAY_ll()
IF (ALLOCATED(ZCONCT)) DEALLOCATE(ZCONCT)
IF (ALLOCATED(ZCHEMT)) DEALLOCATE(ZCHEMT)
IF (ALLOCATED(ZCHEMTI)) DEALLOCATE(ZCHEMTI)
!
!*      1.   U FIELD TREATMENT
!            -----------------
PLBXUS=0.
PLBYUS=0.
!
!*      1.1  Horizontal Bikhardt interpolation
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,2,       &
                                           HLBCX,HLBCY,ZTUT,ZWORK)
DEALLOCATE(ZTUT)
!
ILBX2=SIZE(PLBXUM,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF (LWEST_ll() .AND. ILBX/=0) THEN
  PLBXUS(1:ILBX,IJB:IJE,:)=ZWORK(IIB+1:IIB+ILBX,IJB:IJE,:)  !  C grid
ENDIF
!
IF (LEAST_ll()  .AND. ILBX/=0) THEN
  PLBXUS(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
ENDIF
!
ILBY2=SIZE(PLBYUM,2)
IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
  ILBY=ILBY2/2
ELSE
  ILBY=ILBY2
ENDIF
!
IF(LSOUTH_ll( ) .AND. ILBY/=0) THEN
  PLBYUS(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB:IJB-1+ILBY,:)
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBYUS(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!
!
!*      1.2  Vertical interpolation and tendency computation
!
IF ( SIZE(PLBXUS,1) /= 0 ) THEN
  IF( GVERT_INTERP )   PLBXUS(:,:,:) =               &
     VER_INTERP_LIN(PLBXUS(:,:,:), KKLIN_LBXU(:,:,:),PCOEFLIN_LBXU(:,:,:))
  PLBXUS(:,:,:) = (PLBXUS(:,:,:) - PLBXUM(:,:,:)) / ZTIME
END IF
!
IF ( SIZE(PLBYUS,1) /= 0 ) THEN
  IF( GVERT_INTERP ) PLBYUS(:,:,:) =                  &
     VER_INTERP_LIN(PLBYUS(:,:,:), KKLIN_LBYU(:,:,:),PCOEFLIN_LBYU(:,:,:))
  PLBYUS(:,:,:) = (PLBYUS(:,:,:) - PLBYUM(:,:,:)) / ZTIME
END IF
!
!-------------------------------------------------------------------------------
!
!*      2.   V FIELD TREATMENT
!            -----------------
PLBXVS=0.
PLBYVS=0.
!
!*      2.1  Horizontal Bikhardt interpolation
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,3,       &
                                           HLBCX,HLBCY,ZTVT,ZWORK)
DEALLOCATE(ZTVT)
!
ILBX2=SIZE(PLBXVM,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  PLBXVS(1:ILBX,IJB:IJE,:)=ZWORK(IIB:IIB-1+ILBX,IJB:IJE,:)
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  PLBXVS(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
ENDIF
!
ILBY2=SIZE(PLBYVM,2)
IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
  ILBY=ILBY2/2
ELSE
  ILBY=ILBY2
ENDIF
!
IF(LSOUTH_ll( ) .AND. ILBY/=0) THEN
  PLBYVS(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB+1:IJB+ILBY,:)  !  C grid
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBYVS(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!

!*      2.2  Vertical interpolation and tendency computation
!
IF ( SIZE(PLBXVS,1) /= 0 ) THEN
  IF( GVERT_INTERP ) PLBXVS(:,:,:) =                                  &
     VER_INTERP_LIN(PLBXVS(:,:,:),  KKLIN_LBXV(:,:,:),PCOEFLIN_LBXV(:,:,:))
  PLBXVS(:,:,:) = (PLBXVS(:,:,:) - PLBXVM(:,:,:)) / ZTIME
END IF
!
IF ( SIZE(PLBYVS,1) /= 0 ) THEN
  IF( GVERT_INTERP ) PLBYVS(:,:,:) =                        &
     VER_INTERP_LIN(PLBYVS(:,:,:),  KKLIN_LBYV(:,:,:),PCOEFLIN_LBYV(:,:,:))
    PLBYVS(:,:,:) = (PLBYVS(:,:,:) - PLBYVM(:,:,:)) / ZTIME
END IF
!
!-------------------------------------------------------------------------------
!
!*      3.   W FIELD TREATMENT
!            -----------------
PLBXWS=0.
PLBYWS=0.
!
!*      3.1  Horizontal Bikhardt interpolation
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,4,       &
                                           HLBCX,HLBCY,ZTWT,ZWORK)
DEALLOCATE(ZTWT)
!
ILBX2=SIZE(PLBXWM,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  PLBXWS(1:ILBX,IJB:IJE,:)=ZWORK(IIB:IIB-1+ILBX,IJB:IJE,:)
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  PLBXWS(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
ENDIF
!
ILBY2=SIZE(PLBYWM,2)
IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
  ILBY=ILBY2/2
ELSE
  ILBY=ILBY2
ENDIF
!
IF(LSOUTH_ll( ) .AND. ILBY/=0) THEN
  PLBYWS(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB:IJB-1+ILBY,:)
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBYWS(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!
!
!
!*      3.2  Vertical interpolation and tendency computation
!
IF ( SIZE(PLBXWS,1) /= 0 ) THEN
  IF( GVERT_INTERP ) PLBXWS(:,:,:) =                         &
     VER_INTERP_LIN(PLBXWS(:,:,:), KKLIN_LBXW(:,:,:),PCOEFLIN_LBXW(:,:,:))
  PLBXWS(:,:,:) = (PLBXWS(:,:,:) - PLBXWM(:,:,:)) / ZTIME
END IF
!
IF ( SIZE(PLBYWS,1) /= 0 ) THEN
  IF( GVERT_INTERP ) PLBYWS(:,:,:) =                         &
     VER_INTERP_LIN(PLBYWS(:,:,:), KKLIN_LBYW(:,:,:),PCOEFLIN_LBYW(:,:,:))
  PLBYWS(:,:,:) = (PLBYWS(:,:,:) - PLBYWM(:,:,:)) / ZTIME
END IF
!
!
!-------------------------------------------------------------------------------
!
!*      4.   COMPUTE LARGE SCALE DRY MASS Md SOURCE
!            --------------------------------------
CALL GO_TOMODEL_ll(KDAD, IINFO_ll)
!
IF(OSTEADY_DMASS)  PDRYMASSS   =  0.
!
IF(.NOT. OSTEADY_DMASS) THEN
!
!*       4.0 inner model mask preparation relative to the outer model
!
!
!*       4.1 compute the jacobian J
!
  IF ( CEQNSYS == 'DUR' ) THEN
    IF ( SIZE(XRVREF,1) == 0 ) THEN
      ZJ(:,:,:) = XRHODJ(:,:,:)*XTH00/(XRHODREF(:,:,:)*XTHVREF(:,:,:))
    ELSE
      ZJ(:,:,:) = XRHODJ(:,:,:)*XTH00/(XRHODREF(:,:,:)*XTHVREF(:,:,:) &
                                                 *(1.+XRVREF(:,:,:)))
    END IF
  ELSEIF ( CEQNSYS == 'MAE' .OR. CEQNSYS == 'LHE' ) THEN
    ZJ(:,:,:) = XRHODJ(:,:,:)/XRHODREF(:,:,:)
  END IF
!
!*       4.2 computing of dry density at t
!
  IF(SIZE(XRT,4) == 0) THEN
                        ! dry air case
!                         ------------
      ZRHOD(:,:,:) = XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**ZRD_O_CPD       &
                                  /(XRD*XTHT(:,:,:))
  ELSE                  ! moist air case
!                         --------------
      ZRHOD(:,:,:) = XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**ZRD_O_CPD       &
                                  /(XRD*XTHT(:,:,:)*(1.+ZRV_O_RD*XRT(:,:,:,1)))
  ENDIF
!
!*       4.3 computing of the dry mass at t
!
!
  ZDRYMASST = SUM3D_ll (ZJ(:,:,:)*ZRHOD(:,:,:),IINFO_ll,NXOR_ALL(KMI)+JPHEXT,NYOR_ALL(KMI)+JPHEXT, &
            1+JPVEXT,NXEND_ALL(KMI)-JPHEXT,NYEND_ALL(KMI)-JPHEXT,SIZE(XRHODJ,3)-JPVEXT)
!
!
!*       4.4 normal processing (not at the segment beginning)
!
  IF( KTCOUNT /= 1 ) THEN
    ZBIGTSTEP =  PTSTEP*KDTRATIO                 ! time step of the dad model
    ZDRYMASSM =  PDRYMASST - PDRYMASSS*ZBIGTSTEP ! backward integration over this time step
    PDRYMASST =  ZDRYMASST
    PDRYMASSS = (PDRYMASST - ZDRYMASSM) / ZBIGTSTEP
  ELSE
!
!*       4.5 segment beginning (we have first to recover the dry mass at T-DT)
!
    IF(SIZE(XRT,4) == 0) THEN
                          ! dry air case
!                           ------------
      ZRHOD(:,:,:) = XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**ZRD_O_CPD/(XRD*XTHT(:,:,:))
    ELSE                  ! moist air case
!                           --------------
      ZRHOD(:,:,:) = XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**ZRD_O_CPD/(XRD*XTHT(:,:,:) &
                                                       *(1.+ZRV_O_RD*XRT(:,:,:,1)))
    ENDIF
!
!
    ZDRYMASSM = SUM3D_ll (ZJ(:,:,:)*ZRHOD(:,:,:),IINFO_ll,NXOR_ALL(KMI)+JPHEXT,NYOR_ALL(KMI)+JPHEXT, &
            1+JPVEXT,NXEND_ALL(KMI)-JPHEXT,NYEND_ALL(KMI)-JPHEXT,SIZE(XRHODJ,3)-JPVEXT)
!
    PDRYMASST =  ZDRYMASST
    PDRYMASSS = (PDRYMASST - ZDRYMASSM) / (PTSTEP*KDTRATIO)
  ENDIF
!
END IF
DEALLOCATE(ZJ,ZRHOD)
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
!-------------------------------------------------------------------------------
!
!*      5.   COMPUTE LARGE SCALE SOURCES FOR POTENTIAL TEMPERATURE
!            -----------------------------------------------------
!
!
CALL COMPUTE_LB_M(PLBXTHM,PLBYTHM,PLBXTHS,PLBYTHS,ZTTHT,XTH00)
DEALLOCATE(ZTTHT)
!
!
!-------------------------------------------------------------------------------
!
!*      6.   COMPUTE LARGE SCALE SOURCES FOR TURBULENT KINETIC ENERGY
!            --------------------------------------------------------
!
IF (SIZE(XTKET,3) == 0 .OR. SIZE(PLBXTKEM,3) == 0) THEN
  PLBXTKES(:,:,:) = 0.                      ! turbulence not activated
  PLBYTKES(:,:,:) = 0.
ELSE
  CALL COMPUTE_LB_M(PLBXTKEM,PLBYTKEM,PLBXTKES,PLBYTKES,ZTTKET)
  DEALLOCATE(ZTTKET)
END IF
!
!-------------------------------------------------------------------------------
!
!*      7.   COMPUTE LARGE SCALE SOURCES FOR MOIST VARIABLES
!            -----------------------------------------------
!

IF (IRR == 0) THEN
  PLBXRS(:,:,:,:) = 0.                      ! water cycle not activated
  PLBYRS(:,:,:,:) = 0.
ELSE
  DO JRR = 1,IRR
    CALL COMPUTE_LB_M(PLBXRM(:,:,:,JRR),PLBYRM(:,:,:,JRR),   &
                   PLBXRS(:,:,:,JRR),PLBYRS(:,:,:,JRR),ZTRT(:,:,:,JRR))
  END DO
  DEALLOCATE(ZTRT)
!
  IF ( SIZE(PLBXRS,1) /= 0 ) PLBXRS(:,:,:,IRR+1:SIZE(PLBXRS,4)) = 0.
  IF ( SIZE(PLBYRS,1) /= 0 ) PLBYRS(:,:,:,IRR+1:SIZE(PLBYRS,4)) = 0.
!
END IF
!
!-------------------------------------------------------------------------------
!
!*      8.   COMPUTE LARGE SCALE SOURCES FOR SCALAR VARIABLES
!            ------------------------------------------------
!
IF (NSV_A(KMI) > 0) THEN
  ! Users scalar variables
  DO JSV = 1,NSV_A(KMI)
    CALL COMPUTE_LB_M(PLBXSVM(:,:,:,JSV),PLBYSVM(:,:,:,JSV),  &
         PLBXSVS(:,:,:,JSV),PLBYSVS(:,:,:,JSV),ZTSVT(:,:,:,JSV))
  END DO
  DEALLOCATE(ZTSVT)
ELSE
  PLBXSVS(:,:,:,:) = 0.
  PLBYSVS(:,:,:,:) = 0.
END IF
!
DEALLOCATE(ZWORK)
DEALLOCATE(ZCOEFLIN_LBXM_RED,ZCOEFLIN_LBYM_RED,IKLIN_LBXM_RED,IKLIN_LBYM_RED)
!
!------------------------------------------------------------------------------
CALL GOTO_MODEL(KMI)
!
CONTAINS
!
!
! #############################################################
  SUBROUTINE COMPUTE_LB_M(PLBXM,PLBYM,PLBXS,PLBYS,PTFIELD,PTH00)
! #############################################################
!
REAL, DIMENSION(:,:,:), INTENT(IN) ::  PLBXM,PLBYM !LB fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PLBXS,PLBYS ! LB source terms
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTFIELD     ! ls forcing array
REAL, OPTIONAL, INTENT(IN) :: PTH00 ! reference temperature
!
!
IF(PRESENT(PTH00)) THEN
  PLBXS=PTH00 ! to avoid undefined computation
  PLBYS=PTH00
ELSE
  PLBXS=0.
  PLBYS=0.
ENDIF
!
!*    Horizontal Bikhardt interpolation
!
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,       &
                                           HLBCX,HLBCY,PTFIELD,ZWORK)
!
ILBX2=SIZE(PLBXM,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  PLBXS(1:ILBX,IJB:IJE,:)=ZWORK(IIB:IIB-1+ILBX,IJB:IJE,:)
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  PLBXS(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
ENDIF
!
ILBY2=SIZE(PLBYM,2)
IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
  ILBY=ILBY2/2
ELSE
  ILBY=ILBY2
ENDIF
!
IF(LSOUTH_ll( ) .AND. ILBY/=0) THEN
  PLBYS(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB:IJB-1+ILBY,:)
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBYS(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!
!
!*  Vertical interpolation and tendency
!
!
IF ( SIZE(PLBXS,1) /= 0 ) THEN
  IF( GVERT_INTERP ) THEN
    IF ( ILBX == KRIMX+JPHEXT ) THEN
      PLBXS(:,:,:) = VER_INTERP_LIN(PLBXS(:,:,:),  &
                          KKLIN_LBXM(:,:,:),PCOEFLIN_LBXM(:,:,:))
    ELSE
      PLBXS(:,:,:) = VER_INTERP_LIN(PLBXS(:,:,:),  &
                               IKLIN_LBXM_RED(:,:,:),ZCOEFLIN_LBXM_RED(:,:,:))
    ENDIF
  ENDIF
  PLBXS(:,:,:) = (PLBXS(:,:,:) - PLBXM(:,:,:)) / ZTIME
END IF
!
!
IF ( SIZE(PLBYS,1) /= 0 ) THEN
  IF( GVERT_INTERP ) THEN
    IF ( ILBY == KRIMY+JPHEXT ) THEN
      PLBYS(:,:,:) = VER_INTERP_LIN(PLBYS(:,:,:),  &
                          KKLIN_LBYM(:,:,:),PCOEFLIN_LBYM(:,:,:))
    ELSE
      PLBYS(:,:,:) = VER_INTERP_LIN(PLBYS(:,:,:),  &
                          IKLIN_LBYM_RED(:,:,:),ZCOEFLIN_LBYM_RED(:,:,:))
    ENDIF
  ENDIF
  PLBYS(:,:,:) = (PLBYS(:,:,:) - PLBYM(:,:,:)) / ZTIME
END IF
!
END SUBROUTINE  COMPUTE_LB_M
!
!
END SUBROUTINE ONE_WAY_n
