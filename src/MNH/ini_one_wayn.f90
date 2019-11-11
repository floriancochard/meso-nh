!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_INI_ONE_WAY_n
!     #######################
!
INTERFACE 
!
      SUBROUTINE INI_ONE_WAY_n( KDAD,HLUOUT,PTSTEP,KMI,KTCOUNT,          &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,     &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,     &
                    KDXRATIO,KDYRATIO,KDTRATIO,                          &
                    HLBCX,HLBCY,KRIMX,KRIMY,                             &
                    KKLIN_LBXU,PCOEFLIN_LBXU,KKLIN_LBYU,PCOEFLIN_LBYU,   &
                    KKLIN_LBXV,PCOEFLIN_LBXV,KKLIN_LBYV,PCOEFLIN_LBYV,   &
                    KKLIN_LBXW,PCOEFLIN_LBXW,KKLIN_LBYW,PCOEFLIN_LBYW,   &
                    KKLIN_LBXM,PCOEFLIN_LBXM,KKLIN_LBYM,PCOEFLIN_LBYM,   &
                    HCLOUD, OUSECHAQ, OUSECHIC,                          &
                    PLBXUM,PLBYUM,PLBXVM,PLBYVM,PLBXWM,PLBYWM,           &
                    PLBXTHM,PLBYTHM,                                     &
                    PLBXTKEM,PLBYTKEM,                                   &
                    PLBXRM,PLBYRM,PLBXSVM,PLBYSVM                        )
!
!
INTEGER,          INTENT(IN)    :: KDAD     !  Number of the DAD model
CHARACTER (LEN=*),INTENT(IN)    :: HLUOUT   ! name of the output-listing
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
CHARACTER (LEN=4), INTENT(IN)           :: HCLOUD        ! Indicator of the cloud scheme
LOGICAL,           INTENT(IN)           :: OUSECHAQ      ! logical for aqueous phase chemistry
LOGICAL,           INTENT(IN)           :: OUSECHIC      ! logical for ice phase chemistry
!  
REAL, DIMENSION(:,:,:), INTENT(OUT)    :: PLBXUM,PLBXVM,PLBXWM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(OUT)    :: PLBYUM,PLBYVM,PLBYWM 
REAL, DIMENSION(:,:,:),  INTENT(OUT)  :: PLBXTHM ,PLBYTHM  ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:),  INTENT(OUT)  :: PLBXTKEM,PLBYTKEM ! Theta, TKE
REAL, DIMENSION(:,:,:,:),INTENT(OUT)  :: PLBXRM  ,PLBYRM   ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(OUT)  :: PLBXSVM ,PLBYSVM  ! in x and y-dir.
!
END SUBROUTINE INI_ONE_WAY_n
!
END INTERFACE
!
END MODULE MODI_INI_ONE_WAY_n
!

!     ####################################################################
SUBROUTINE INI_ONE_WAY_n(KDAD,HLUOUT,PTSTEP,KMI,KTCOUNT,                 &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,     &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,     &
                    KDXRATIO,KDYRATIO,KDTRATIO,                          &
                    HLBCX,HLBCY,KRIMX,KRIMY,                             &
                    KKLIN_LBXU,PCOEFLIN_LBXU,KKLIN_LBYU,PCOEFLIN_LBYU,   &
                    KKLIN_LBXV,PCOEFLIN_LBXV,KKLIN_LBYV,PCOEFLIN_LBYV,   &
                    KKLIN_LBXW,PCOEFLIN_LBXW,KKLIN_LBYW,PCOEFLIN_LBYW,   &
                    KKLIN_LBXM,PCOEFLIN_LBXM,KKLIN_LBYM,PCOEFLIN_LBYM,   &
                    HCLOUD,OUSECHAQ,OUSECHIC,                            &
                    PLBXUM,PLBYUM,PLBXVM,PLBYVM,PLBXWM,PLBYWM,           &
                    PLBXTHM,PLBYTHM,                                     &
                    PLBXTKEM,PLBYTKEM,                                   &
                    PLBXRM,PLBYRM,PLBXSVM,PLBYSVM                        )
!     ####################################################################
!
!!****  *INI_ONE_WAY$n* - INItializing a nested model Large Scale sources
!!
!!    PURPOSE
!!    -------
!!      The purpose of INI_ONE_WAY$n is to initialize Large Scale sources
!!    of all the prognostic variables of the current nested model when the
!!    current time step is in phase with its outer (DAD) model $n.
!
!
!!**  METHOD
!!    ------
!!      The basic task consists in interpolating fields from outer model $n
!!    to present inner model, using horizontal Bikhardt interpolation.
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
!!      Module MODD_FIELD$n : XUM,XVM,XWM,XRM,XTHM
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    J. P. Lafore and J. Stein  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     22/09/99
!!    J.-P. Pinty    29/11/02 modify the LB*SVS for the C3R5 scheme
!!                            and add ICE2, ICE4, CELEC
!!    Modification   03/2006   (O.Geoffroy) add KHKO schem
!!    Modification   05/2006   Remove KEPS
!!    M. Leriche     11/2009  modify the LB*SVS for the aqueous phase chemistry
!!                   07/2010  idem for ice phase chemical species
!!    Bosseur & Filippi 07/2013 Adds Forefire
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 18/12/2015 : Correction of bug in bound in // for NHALO <>1 
!!      B.VIE   2016 : LIMA
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_CST
USE MODD_FIELD_n      ! modules relative to the outer model $n
USE MODD_PARAM_n
USE MODD_CH_MNHC_n, ONLY: LUSECHAQ, LUSECHIC
USE MODD_REF_n
USE MODD_NSV
!
USE MODI_BIKHARDT
USE MODI_VER_INTERP_LIN
USE MODI_SET_CONC_RAIN_C2R2
USE MODI_SET_CONC_ICE_C1R3
USE MODI_SET_CHEMAQ_1WAY
!
USE MODI_SET_CONC_LIMA
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
INTEGER,          INTENT(IN)    :: KDAD     !  Number of the DAD model
CHARACTER (LEN=*),INTENT(IN)    :: HLUOUT   ! name for output-listing
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
CHARACTER (LEN=4), INTENT(IN)  :: HCLOUD        ! Indicator of the cloud scheme
LOGICAL,           INTENT(IN)  :: OUSECHAQ      ! logical for aqueous phase
LOGICAL,           INTENT(IN)  :: OUSECHIC      ! logical for ice phase chemistry
!
REAL, DIMENSION(:,:,:), INTENT(OUT)    :: PLBXUM,PLBXVM,PLBXWM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(OUT)    :: PLBYUM,PLBYVM,PLBYWM
REAL, DIMENSION(:,:,:),  INTENT(OUT)  :: PLBXTHM ,PLBYTHM  ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:),  INTENT(OUT)  :: PLBXTKEM,PLBYTKEM ! Theta, TKE
REAL, DIMENSION(:,:,:,:),INTENT(OUT)  :: PLBXRM  ,PLBYRM   ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(OUT)  :: PLBXSVM ,PLBYSVM  ! in x and y-dir.
!
!
!*       0.2   declarations of local variables
!
INTEGER                :: IIB,IIE,IJB,IJE,IIU,IJU
INTEGER                :: ILBX,ILBY,ILBX2,ILBY2
REAL,   DIMENSION(:,:,:), ALLOCATABLE  :: ZWORK
LOGICAL  :: GVERT_INTERP
!
INTEGER           :: IRR,ISV_USER          !  Number of moist and user scalar variables
INTEGER           :: JRR,JSV          !  Loop index
!
! reduced array for the interpolation coefficients
REAL, DIMENSION(:,:,:), ALLOCATABLE ::  ZCOEFLIN_LBXM_RED,ZCOEFLIN_LBYM_RED
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IKLIN_LBXM_RED,IKLIN_LBYM_RED
!
! Variables used for LS communications
INTEGER :: IINFO_ll, IDIMX, IDIMY
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTUM, ZTVM, ZTWM, ZTTHM, ZTTKEM
REAL, DIMENSION(:,:,:,:), ALLOCATABLE ::ZTRM,ZTSVM
!
CHARACTER(LEN=4)                    :: ZINIT_TYPE ! type of C2R2 initialisation
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCONCM  ! C2R2 concentrations
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCHEMM  ! chemical concentrations
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCHEMMI  ! chemical ice phase concentrations
!-------------------------------------------------------------------------------
!
!*      0.   INITIALISATION
! 
CALL GOTO_MODEL(KDAD)
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
!
GVERT_INTERP=.TRUE.
!
IRR=MIN(SIZE(XRT,4),SIZE(PLBXRM,4))
ISV_USER=MIN(NSV_USER_A(KDAD),NSV_USER_A(KMI))
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
!*      1 GATHER LS FIELD FOR THE CHILD MODEL KMI
!
!       1.1  Must be on the father model to call get_child_dim
!
CALL GO_TOMODEL_ll(KDAD, IINFO_ll)
CALL GET_CHILD_DIM_ll(KMI, IDIMX, IDIMY, IINFO_ll)
!
!       1.2  Allocate array which will receive coarse grid points
!
ALLOCATE(ZTUM(IDIMX,IDIMY,SIZE(XUT,3)))
ALLOCATE(ZTVM(IDIMX,IDIMY,SIZE(XVT,3)))
ALLOCATE(ZTWM(IDIMX,IDIMY,SIZE(XWT,3)))
ALLOCATE(ZTTHM(IDIMX,IDIMY,SIZE(XTHT,3)))
IF (SIZE(XTKET) /= 0) ALLOCATE(ZTTKEM(IDIMX,IDIMY,SIZE(XTKET,3)))
IF (IRR /= 0) ALLOCATE(ZTRM(IDIMX,IDIMY,SIZE(XRT,3),IRR))
IF (NSV_A(KMI)/= 0) ALLOCATE(ZTSVM(IDIMX,IDIMY,SIZE(XRT,3),NSV_A(KMI)))
!
!         1.3  Specify the ls "source" fields and receiver fields
!
CALL SET_LSFIELD_1WAY_ll(XUT,ZTUM,KMI)
CALL SET_LSFIELD_1WAY_ll(XVT,ZTVM,KMI)
CALL SET_LSFIELD_1WAY_ll(XWT,ZTWM,KMI)
CALL SET_LSFIELD_1WAY_ll(XTHT,ZTTHM,KMI)
IF (ALLOCATED(ZTTKEM)) CALL SET_LSFIELD_1WAY_ll(XTKET,ZTTKEM,KMI)
!
DO JRR=1,IRR
  CALL SET_LSFIELD_1WAY_ll(XRT(:,:,:,JRR),ZTRM(:,:,:,JRR),KMI)
ENDDO
!
! USERs scalar variables
!
IF (ALLOCATED(ZTSVM)) ZTSVM=0.
DO JSV=1,ISV_USER
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV),ZTSVM(:,:,:,JSV),KMI)
ENDDO
!  Checking if it is necessary to compute the Nc and Nr
!  concentrations to use the C2R2 microphysical scheme
!  (FATHER does not use C2R2(or KHKO) and CHILD uses C2R2(or KHKO))
IF ( HCLOUD=="C2R2" .OR. HCLOUD=="KHKO" ) THEN    
 IF ( CCLOUD/="NONE" .AND. CCLOUD/="C2R2" .AND. CCLOUD/="KHKO" ) THEN
  ZINIT_TYPE="NONE"
  ALLOCATE(ZCONCM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),3))
  IF (CCLOUD == "REVE") THEN
    ZINIT_TYPE = "INI1"
  ELSE IF (CCLOUD == "KESS" ) THEN
    ZINIT_TYPE = "INI2"
  END IF
  CALL SET_CONC_RAIN_C2R2(ZINIT_TYPE,XRHODREF,XRT,ZCONCM)
  DO JSV=1,3
    CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  ENDDO
 ELSE
  DO JSV=1,NSV_C2R2_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  END DO
 ENDIF
ENDIF
!
!  Checking also if it is necessary to compute the Ni
!  concentrations to use the C3R5 microphysical scheme
!  (FATHER does not use C3R5 and CHILD uses C3R5)
!
IF (HCLOUD=="C3R5"  ) THEN
 IF (CCLOUD(1:3)=="ICE") THEN
  ZINIT_TYPE="NONE"
  ALLOCATE(ZCONCM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),5))
  IF (CCLOUD == "REVE") THEN
      ZINIT_TYPE = "INI1"
  ELSE IF (CCLOUD == "KESS" ) THEN
    ZINIT_TYPE = "INI2"
  END IF
  CALL SET_CONC_RAIN_C2R2(ZINIT_TYPE,XRHODREF,XRT,ZCONCM)
  DO JSV=1,3
    CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  ENDDO
  ZINIT_TYPE="INI3"
  CALL SET_CONC_ICE_C1R3 (XRHODREF,XRT,ZCONCM)
  DO JSV=4,5
    CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV), &
         ZTSVM(:,:,:,JSV-4+NSV_C1R3BEG_A(KMI)),KMI)
  ENDDO
 ELSE
  DO JSV=1,NSV_C2R2_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  END DO
  DO JSV=1,NSV_C1R3_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C1R3BEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_C1R3BEG_A(KMI)),KMI)
  END DO
 ENDIF
ENDIF
!
!  Checking if it is necessary to compute the Nc, Nr, Ni
!  concentrations to use the LIMA microphysical scheme
!  (FATHER does not use LIMA and CHILD uses LIMA)
!
IF (HCLOUD=="LIMA"  ) THEN
   IF (CCLOUD/="LIMA") THEN
      ALLOCATE(ZCONCM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),NSV_LIMA_A(KMI)))
      IF (CCLOUD == "REVE") THEN
         ZINIT_TYPE = "INI1"
      ELSE
         ZINIT_TYPE = "NONE"
      END IF
      CALL SET_CONC_LIMA (ZINIT_TYPE,XRHODREF,XRT,ZCONCM)
      DO JSV=1,NSV_LIMA_A(KMI)
         CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV),&
              &ZTSVM(:,:,:,JSV-1+NSV_LIMA_BEG_A(KMI)),KMI)
      ENDDO   
   ELSE
      IF (NSV_LIMA_A(KMI)/=NSV_LIMA_A(KDAD)) CALL ABORT
      DO JSV=1,NSV_LIMA_A(KMI)
         CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LIMA_BEG_A(KDAD)),&
              &ZTSVM(:,:,:,JSV-1+NSV_LIMA_BEG_A(KMI)),KMI)
      END DO
   END IF
ENDIF
!
! electrical variables
!
DO JSV=1,MIN(NSV_ELEC_A(KMI),NSV_ELEC_A(KDAD))
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_ELECBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_ELECBEG_A(KMI)),KMI)
END DO
!
! chemical Scalar variables
!  Checking if it is necessary to compute the Caq
!  concentrations to use the aqueous phase chemistry
!  (FATHER does not use aqueous phase chemistry and CHILD uses it)
!
IF (OUSECHAQ) THEN
  IF (.NOT.(LUSECHAQ)) THEN
    ALLOCATE(ZCHEMM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),&
                   NSV_CHEM_A(KMI)))
    CALL SET_CHEMAQ_1WAY(XRHODREF,&
         XSVT(:,:,:,NSV_CHEMBEG_A(KDAD):NSV_CHEMEND_A(KDAD)),ZCHEMM)
    DO JSV=1,NSV_CHEM_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(ZCHEMM(:,:,:,JSV),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
    ENDDO
  ELSE
    DO JSV=1,NSV_CHEM_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KDAD)),&
            &ZTSVM(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
    END DO
  ENDIF
ELSE
  DO JSV=1,NSV_CHEM_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
  END DO
ENDIF
!  Checking if it is necessary to compute the Cic
!  concentrations to use the ice phase chemistry
!  (FATHER does not use ice phase chemistry and CHILD uses it)
!
IF (OUSECHIC) THEN
  IF (.NOT.(LUSECHIC)) THEN
    ALLOCATE(ZCHEMMI(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),&
                   NSV_CHIC_A(KMI)))
    ZCHEMMI(:,:,:,:) = 0.
    DO JSV=1,NSV_CHIC_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(ZCHEMMI(:,:,:,JSV),&
       &ZTSVM(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
    ENDDO
  ELSE
    DO JSV=1,NSV_CHIC_A(KMI)
       CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
    END DO
  ENDIF
ELSE
  DO JSV=1,NSV_CHIC_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
  END DO
ENDIF
!
!
! lagrangian variables
DO JSV=1,NSV_LG_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LGBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_LGBEG_A(KMI)),KMI)
END DO
!
! NOX                     
DO JSV=1,NSV_LNOX_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LNOXBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_LNOXBEG_A(KMI)),KMI)
END DO
!
! Dust Scalar variables
DO JSV=1,NSV_DST_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_DSTBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_DSTBEG_A(KMI)),KMI)
END DO
!
! Moist Dust Scalar variables
DO JSV=1,NSV_DSTDEP_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_DSTDEPBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_DSTDEPBEG_A(KMI)),KMI)
END DO

! Sea Salt Scalar variables
DO JSV=1,NSV_SLT_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_SLTBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_SLTBEG_A(KMI)),KMI)
END DO
!
! Moist Sea Salt Scalar variables
DO JSV=1,NSV_SLTDEP_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_SLTDEPBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_SLTDEPBEG_A(KMI)),KMI)
END DO
!
!
! Passive pollutant      
DO JSV=1,NSV_PP_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_PPBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_PPBEG_A(KMI)),KMI)
END DO
#ifdef MNH_FOREFIRE
! ForeFire variables      
DO JSV=1,NSV_FF_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_FFBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_FFBEG_A(KMI)),KMI)
END DO
#endif
!        1.4  Communication
!
CALL LS_FORCING_ll(KMI,IINFO_ll)
!
!        1.5  Back to the (current) child model
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
CALL UNSET_LSFIELD_1WAY_ll()
IF (ALLOCATED(ZCONCM)) DEALLOCATE(ZCONCM)
IF (ALLOCATED(ZCHEMM)) DEALLOCATE(ZCHEMM)
IF (ALLOCATED(ZCHEMMI)) DEALLOCATE(ZCHEMMI)
!
!
!-------------------------------------------------------------------------------
!
!*      1.   U FIELD TREATMENT
!            -----------------
!
!*      1.1  Horizontal Bikhardt interpolation
!
PLBXUM=0.
PLBYUM=0.
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,2,       &
                                           HLBCX,HLBCY,ZTUM,ZWORK)
DEALLOCATE(ZTUM)
!
ILBX2=SIZE(PLBXUM,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  PLBXUM(1:ILBX,IJB:IJE,:)=ZWORK(IIB+1:IIB+ILBX,IJB:IJE,:)  !  C grid
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  PLBXUM(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
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
  PLBYUM(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB:IJB-1+ILBY,:)
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBYUM(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!
!*      1.2  Vertical interpolation
!
IF ( SIZE(PLBXUM,1) /= 0 .AND. GVERT_INTERP) THEN
  PLBXUM(:,:,:) = VER_INTERP_LIN(PLBXUM(:,:,:),   &
                                   KKLIN_LBXU(:,:,:),PCOEFLIN_LBXU(:,:,:))
END IF
!
IF ( SIZE(PLBYUM,1) /= 0 .AND. GVERT_INTERP) THEN
  PLBYUM(:,:,:) = VER_INTERP_LIN(PLBYUM(:,:,:),   &
                                   KKLIN_LBYU(:,:,:),PCOEFLIN_LBYU(:,:,:))
END IF
!
!-------------------------------------------------------------------------------
!
!*      2.   V FIELD TREATMENT
!            -----------------
!
!*      2.1  Horizontal Bikhardt interpolation
!
PLBXVM=0.
PLBYVM=0.
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,3,       &
                                           HLBCX,HLBCY,ZTVM,ZWORK)
DEALLOCATE(ZTVM)
!
ILBX2=SIZE(PLBXVM,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  PLBXVM(1:ILBX,IJB:IJE,:)=ZWORK(IIB:IIB-1+ILBX,IJB:IJE,:)
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  PLBXVM(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
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
  PLBYVM(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB+1:IJB+ILBY,:)  !  C grid
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBYVM(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!
!*      1.2  Vertical interpolation
!
IF ( SIZE(PLBXVM,1) /= 0 .AND. GVERT_INTERP) THEN
  PLBXVM(:,:,:) = VER_INTERP_LIN(PLBXVM(:,:,:),   &
                                   KKLIN_LBXV(:,:,:),PCOEFLIN_LBXV(:,:,:))
END IF
!
IF ( SIZE(PLBYVM,1) /= 0 .AND. GVERT_INTERP) THEN
  PLBYVM(:,:,:) = VER_INTERP_LIN(PLBYVM(:,:,:),   &
                                   KKLIN_LBYV(:,:,:),PCOEFLIN_LBYV(:,:,:))
END IF

!-------------------------------------------------------------------------------
!
!*      3.   W FIELD TREATMENT
!            -----------------
!
!*      3.1  Horizontal Bikhardt interpolation
!
PLBXWM=0.
PLBYWM=0.
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,4,       &
                                           HLBCX,HLBCY,ZTWM,ZWORK)
DEALLOCATE(ZTWM)
!
ILBX2=SIZE(PLBXWM,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  PLBXWM(1:ILBX,IJB:IJE,:)=ZWORK(IIB:IIB-1+ILBX,IJB:IJE,:)
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  PLBXWM(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
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
  PLBYWM(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB:IJB-1+ILBY,:)
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBYWM(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!
!*      1.2  Vertical interpolation
!
IF ( SIZE(PLBXWM,1) /= 0 .AND. GVERT_INTERP) THEN
  PLBXWM(:,:,:) = VER_INTERP_LIN(PLBXWM(:,:,:),   &
                                   KKLIN_LBXW(:,:,:),PCOEFLIN_LBXW(:,:,:))
END IF
!
IF ( SIZE(PLBYWM,1) /= 0 .AND. GVERT_INTERP) THEN
  PLBYWM(:,:,:) = VER_INTERP_LIN(PLBYWM(:,:,:),   &
                                   KKLIN_LBYW(:,:,:),PCOEFLIN_LBYW(:,:,:))
END IF
!
!
!
!-------------------------------------------------------------------------------
!
!*      5.   COMPUTE LARGE SCALE SOURCES FOR POTENTIAL TEMPERATURE
!            -----------------------------------------------------
!
CALL COMPUTE_LB_M(PLBXTHM,PLBYTHM,ZTTHM,XTH00)
!
DEALLOCATE(ZTTHM)
!
!
!-------------------------------------------------------------------------------
!
!*      6.   COMPUTE LARGE SCALE SOURCES FOR TURBULENT KINETIC ENERGY
!            --------------------------------------------------------
!
!
IF (SIZE(XTKET,3) == 0 .OR. SIZE(PLBXTKEM,3) == 0) THEN
  PLBXTKEM(:,:,:) = 0.                      ! turbulence not activated
  PLBYTKEM(:,:,:) = 0.
ELSE
  CALL COMPUTE_LB_M(PLBXTKEM,PLBYTKEM,ZTTKEM)
  DEALLOCATE(ZTTKEM)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      7.   COMPUTE LARGE SCALE SOURCES FOR MOIST VARIABLES
!            -----------------------------------------------
!
!
IF (IRR == 0 ) THEN
  PLBXRM(:,:,:,:) = 0.                      ! water cycle not activated
  PLBYRM(:,:,:,:) = 0.
ELSE
  DO JRR = 1,IRR
    CALL COMPUTE_LB_M(PLBXRM(:,:,:,JRR),PLBYRM(:,:,:,JRR),ZTRM(:,:,:,JRR))
  END DO
  DEALLOCATE(ZTRM)
!
  IF ( SIZE(PLBXRM,1) /= 0 ) PLBXRM(:,:,:,IRR+1:SIZE(PLBXRM,4)) = 0.
  IF ( SIZE(PLBYRM,1) /= 0 ) PLBYRM(:,:,:,IRR+1:SIZE(PLBYRM,4)) = 0.
!
END IF
!
!-------------------------------------------------------------------------------
!
!*      8.   COMPUTE LARGE SCALE SOURCES FOR SCALAR VARIABLES
!            ------------------------------------------------
!
!
IF (NSV_A(KMI) > 0) THEN
  DO JSV = 1,NSV_A(KMI)
    CALL COMPUTE_LB_M(PLBXSVM(:,:,:,JSV),PLBYSVM(:,:,:,JSV),ZTSVM(:,:,:,JSV))
  END DO
  DEALLOCATE(ZTSVM)
ELSE
  PLBXSVM(:,:,:,:) = 0.
  PLBYSVM(:,:,:,:) = 0.
END IF
!
DEALLOCATE(ZWORK)
DEALLOCATE(ZCOEFLIN_LBXM_RED,ZCOEFLIN_LBYM_RED,IKLIN_LBXM_RED,IKLIN_LBYM_RED)
!
CALL GOTO_MODEL(KMI)
!------------------------------------------------------------------------------
!
CONTAINS
!
!
!     ################################################
      SUBROUTINE COMPUTE_LB_M(PLBX,PLBY,PTFIELD,PTH00)
!     ################################################
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PLBX,PLBY ! source term
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PTFIELD   ! ls forcing array
REAL, OPTIONAL, INTENT(IN)          :: PTH00 ! reference temperature
!
IF(PRESENT(PTH00)) THEN
  PLBX=PTH00 ! to avoid undefined computation
  PLBY=PTH00
ELSE
  PLBX=0.
  PLBY=0.
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
ILBX2=SIZE(PLBX,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  PLBX(1:ILBX,IJB:IJE,:)=ZWORK(IIB:IIB-1+ILBX,IJB:IJE,:)
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  PLBX(ILBX2-ILBX+1:ILBX2,IJB:IJE,:)=ZWORK(IIE+1-ILBX:IIE,IJB:IJE,:)
ENDIF
!
ILBY2=SIZE(PLBY,2)
IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
  ILBY=ILBY2/2
ELSE
  ILBY=ILBY2
ENDIF
!
IF(LSOUTH_ll( ) .AND. ILBY/=0) THEN
  PLBY(IIB:IIE,1:ILBY,:)=ZWORK(IIB:IIE,IJB:IJB-1+ILBY,:)
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  PLBY(IIB:IIE,ILBY2-ILBY+1:ILBY2,:)=ZWORK(IIB:IIE,IJE+1-ILBY:IJE,:)
ENDIF
!
!*    Vertical interpolation
!
IF ( SIZE(PLBX,1) /= 0 .AND. GVERT_INTERP) THEN
  IF ( ILBX == KRIMX+JPHEXT ) THEN
    PLBX(:,:,:) = VER_INTERP_LIN(PLBX(:,:,:),   &
                             KKLIN_LBXM(:,:,:),PCOEFLIN_LBXM(:,:,:))
  ELSE
    PLBX(:,:,:) = VER_INTERP_LIN(PLBX(:,:,:),  &
                          IKLIN_LBXM_RED(:,:,:),ZCOEFLIN_LBXM_RED(:,:,:))
  END IF
END IF
!
IF ( SIZE(PLBY,1) /= 0 .AND. GVERT_INTERP) THEN
  IF ( ILBY == KRIMY+JPHEXT ) THEN
    PLBY(:,:,:) = VER_INTERP_LIN(PLBY(:,:,:),   &
                                   KKLIN_LBYM(:,:,:),PCOEFLIN_LBYM(:,:,:))
  ELSE
    PLBY(:,:,:) = VER_INTERP_LIN(PLBY(:,:,:),   &
                          IKLIN_LBYM_RED(:,:,:),ZCOEFLIN_LBYM_RED(:,:,:))
  END IF
END IF
!
!
END SUBROUTINE  COMPUTE_LB_M
!
END SUBROUTINE INI_ONE_WAY_n
