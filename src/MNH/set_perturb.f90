!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_SET_PERTURB
!     #######################
!
INTERFACE
!
SUBROUTINE SET_PERTURB(TPEXPREFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA), INTENT(IN)  :: TPEXPREFILE ! input data file
!
END SUBROUTINE SET_PERTURB
!
END INTERFACE
!
END MODULE MODI_SET_PERTURB
!
!     ##############################
      SUBROUTINE SET_PERTURB(TPEXPREFILE)
!     ##############################
!
!!****  *SET_PERTURB* - add a perturbation to a larger scale state field
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to add a perturbation to a balanced
!     mass and wind fields. This perturbation is added either to the theta 
!     field to initiate either a convective development, or to the wind field
!     to initiate a normal mode selection or ...
!
!!**  METHOD
!!    ------
!!      
!!      An analytical function is used to prescribe the theta perturbation 
!!    variation and a stream function is used to set the wind perturbation in
!!    not to create a divergent field.
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains parameters
!!        JPHEXT : number of horizontal external points
!!        JPVEXT : number of vertical external points
!!      Module MODD_GRID1 : contains grid variables
!!        XXHAT :  Position x in the conformal
!!                 plane or on the cartesian plane
!!        XYHAT :  Position y in the conformal
!!                 plane or on the cartesian plane
!!        XZHAT :  Position z in the grid without orography
!!      Module MODD_CST : contains physical constants
!!        XPI   :  pi 
!!      Module MODD_FIELD1 : contains prognostic variables
!!        XTHT  : potential temperature 
!!        XUT   : x-component of the wind
!!        XVT   : y-component of the wind
!!      Module MODD_CONF   : model configuration
!!        L2D   : switch for 2D configuration
!!      Module MODD_DIM1   : model dimensions
!!        NIMAX,NJMAX
!!      Module MODD_LBC1   : model lateral boudary conditions
!!        CLBCX(2),CLBCY(2)
!!
!!    REFERENCE
!!    ---------
!!
!!     Book2 subroutine SET_PERTURB
!!
!!    AUTHOR
!!    ------
!!  	J.Stein         * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              6/12/94
!!      I.Mallet and J.Stein  07/01/95 add the perturbation for the wind field
!!      J. Stein              16/03/95 remove R from the historical variables
!!      J. Stein              16/03/95 remove MODD_REF1 which is not used
!!      P. Mascart            01/09/95 change the bubble center computation
!!      J. Stein+F. Lohou     06/02/96 add the white noise perturbation
!!      J.-P. Pinty           06/12/96 add the argument list to control the 
!!                                     perturbation type, amplitude and location
!!      J-P Pinty P Jabouille 30/08/00 parallelization of the code
!!      I. Mallet                06/06 cleaning (namelist inside the routine)
!!      J.Escobar             25/03/2012 optim. parallelization of White noise
!!      J.Escobar             27/03/2012 force identical random seed & correct XOR/YOR global shift 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      C.Lac, V.Masson       1/2018 : White noise in the LBC
!!      Q.Rodier              10/2018 : move allocate(ZWHITE) for NKWH>2
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONF
USE MODD_DIM_n
USE MODD_FIELD_n
USE MODD_GRID_n
USE MODD_IO_ll,   ONLY: TFILEDATA
USE MODD_LBC_n
USE MODD_LUNIT_n, ONLY: CLUOUT, TLUOUT
USE MODD_LSFIELD_n
USE MODD_PARAMETERS
USE MODD_REF_n
!
USE MODE_FM
USE MODE_GATHER_ll
USE MODE_ll
USE MODE_MPPDB
USE MODE_POS
USE MODE_REPRO_SUM
!
USE MODI_GET_HALO
USE MODI_SCATTER
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA), INTENT(IN)  :: TPEXPREFILE ! input data file
!
!*       0.2   Declarations of local variables
!
INTEGER :: ILUPRE,IRESP ! logical unit number of the  EXPRE and FM return code
INTEGER :: ILUOUT       ! Logical unit number for output-listing
LOGICAL :: GFOUND                  ! Return code when searching namelist
!
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZDIST ! distance
                                   ! from the center of the 
                                   ! disturbance               
REAL :: ZCENTERX,ZCENTERY          ! coordinates of the disturbance center
INTEGER          :: IIB,IJB,IKB    ! Begining useful area  in x,y,z directions
INTEGER          :: IIE,IJE,IKE    ! End useful area in x,y,z directions
INTEGER          :: IIU,IJU,IKU    ! Array sizes in i,j,k directions
INTEGER          :: JI,JJ,JK       ! Loop indexes along the x,y,z directions
INTEGER          :: IJ0            ! J value of the center of the perturb.
INTEGER          :: IJMIN,IJMAX    ! J values for the y-width of the perturb.
REAL,DIMENSION(:,:),ALLOCATABLE :: ZPHI_ll        ! stream function
REAL,DIMENSION(:,:),ALLOCATABLE :: ZPU_ll,ZPV_ll  ! hor. wind components 
REAL,DIMENSION(:,:),ALLOCATABLE :: ZPU,ZPV
REAL             :: ZX, ZY, ZRANDOM, ZT1, ZT2, ZAL, ZBL, ZCL, ZDL, ZVAR
INTEGER          :: IKX, IKY, JX, JY              ! loop index in Fourier space
REAL,DIMENSION(:,:),ALLOCATABLE :: ZWHITE         ! white noise perturbation
REAL,DIMENSION(:,:),ALLOCATABLE :: ZWHITE_ll      ! global ZWHITE
REAL,DIMENSION(:,:),ALLOCATABLE :: ZCX_ll,ZSX_ll  ! Fourier coefficient for x
REAL,DIMENSION(:,:),ALLOCATABLE :: ZCY_ll,ZSY_ll  ! and y directions
!
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZT    ! temperature
REAL, DIMENSION(SIZE(XXHAT),SIZE(XYHAT),SIZE(XZHAT)) :: ZHU   ! rel. humidity
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZXHAT_ll    !  Position x in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:), ALLOCATABLE   :: ZYHAT_ll    !   Position y in the conformal
                                                 ! plane (array on the complete domain)
INTEGER :: IIU_ll,IJU_ll ! horizontal,vertical size of the extended global domain 
INTEGER :: IIB_ll,IJB_ll ! global coordinate of the physical global domain 
INTEGER :: IIE_ll,IJE_ll ! 
!
!*       0.3  Namelist declarations
!
CHARACTER(LEN=2) :: CPERT_KIND='TH'! Perturbation type
REAL             :: XAMPLITH=1.5   ! Perturbation amplitude maximum for THETA
REAL             :: XAMPLIRV=0.0   ! Perturbation amplitude maximum for R_V
REAL             :: XAMPLIUV=1.083367E6
                                   ! Perturbation amplitude maximum for U and V
                                   ! initially XAMPLIUV=1.7E5/(2*SIN(XPI/40.))
REAL             :: XAMPLIWH=0.1   ! Perturbation amplitude maximum for the
                                   !                                WHite noise
LOGICAL          :: LWH_LBXU=.FALSE.! White noise in inflow and outflow LBC of U
LOGICAL          :: LWH_LBYV=.FALSE.! White noise in inflow and outflow LBC of V                                   
INTEGER          :: NKWH=2         ! Upper level of the layer
                                   ! where white noise is applied
LOGICAL          :: LSET_RHU=.TRUE.! Conservation of the Relative HUmidity when
                                   !                                     TRUE
REAL             :: XCENTERZ=2000. ! Height of the maximum of TH perturbation
REAL             :: XRADX=  10000. ! X-Radius of the perturbation
REAL             :: XRADY=  10000. ! Y-Radius of the perturbation
REAL             :: XRADZ=   2200. ! Z-Radius of the perturbation
!
REAL             :: ZOMEGA
INTEGER, PARAMETER                 :: i_seed_param = 26032012
INTEGER, DIMENSION(:), ALLOCATABLE :: i_seed
INTEGER                            :: ni_seed
!
INTEGER                            :: IXOR,IYOR,JI_ll,JJ_ll
!
NAMELIST/NAM_PERT_PRE/CPERT_KIND,XAMPLITH,       &! Perturbation parameters
                      XAMPLIRV,XCENTERZ,XRADX,   &!
                      XRADY,XRADZ,LSET_RHU,      &
                      XAMPLIUV,XAMPLIWH,NKWH,LWH_LBXU,LWH_LBYV
!-------------------------------------------------------------------------------
!
!*	 1.     PROLOGUE : 
!               ----------
!
!*       1.1   Retrieve unit numbers and read namelist
!
ILUPRE = TPEXPREFILE%NLU
ILUOUT = TLUOUT%NLU
!
CALL POSNAM(ILUPRE,'NAM_PERT_PRE',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUPRE,NML=NAM_PERT_PRE)
!
!  Distribue/set identical parameter  seed to all the proc
!
CALL RANDOM_SEED(SIZE=ni_seed) ! get size of seed
ALLOCATE(i_seed(ni_seed))
i_seed(:) = i_seed_param !  
CALL RANDOM_SEED(PUT=i_seed)
!
!
!*       1.2
!
IIU=SIZE(XXHAT)
IJU=SIZE(XYHAT)
IKU=SIZE(XZHAT)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
!
IIB=JPHEXT+1
IIE=IIU-JPHEXT
IJB=JPHEXT+1
IJE=IJU-JPHEXT
!
IIU_ll=NIMAX_ll+2*JPHEXT
IJU_ll=NJMAX_ll+2*JPHEXT
ALLOCATE(ZXHAT_ll(IIU_ll),ZYHAT_ll(IJU_ll))
CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP)
CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP)
!  
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IIB_ll=1+JPHEXT
IJB_ll=1+JPHEXT
IIE_ll=IIU_ll-JPHEXT
IJE_ll=IJU_ll-JPHEXT
!
CALL GET_OR_ll('B',IXOR,IYOR)
!
!-------------------------------------------------------------------------------
!
!*	 2.     COMPUTE THE PERTURBATION ON THETA : 
!               -----------------------------------
!
SELECT CASE(CPERT_KIND)
  CASE('TH')
    ZDIST(:,:,:) = 2.
! C grid shift
    ZCENTERX=(ZXHAT_ll(2)+ZXHAT_ll(IIU_ll))*0.5
    ZCENTERY=(ZYHAT_ll(2)+ZYHAT_ll(IJU_ll))*0.5
!
    DO JK =IKB,IKE
      DO JJ = IJB,IJE
        DO JI = IIB,IIE
          ZDIST(JI,JJ,JK) = SQRT(                         &
          (( (XXHAT(JI)+XXHAT(JI+1))*0.5 - ZCENTERX ) / XRADX)**2 + &
          (( (XYHAT(JJ)+XYHAT(JJ+1))*0.5 - ZCENTERY ) / XRADY)**2 + &
          (( (XZHAT(JK)+XZHAT(JK+1))*0.5 - XCENTERZ ) / XRADZ)**2   &
                                )
        END DO
      END DO
    END DO
!
    CALL MPPDB_CHECK3D(ZDIST,"SET_PERTURB::ZDIST",PRECISION)
    !CALL MPPDB_CHECK3D(XTHM,"SET_PERTURB::XTHM",PRECISION)
!
    IF ( LSET_RHU) THEN
       ZT(:,:,:)  = 0.0
       ZHU(:,:,:) = 0.0
      WHERE ( ZDIST(:,:,:) <= 1.)
!
! save the actual relative humidity
!
        ZT(:,:,:) =XTHT(:,:,:)*(XPABST(:,:,:)/XP00) **(XRD/XCPD)
        ZHU(:,:,:)  = (XRT(:,:,:,1)*XPABST(:,:,:)) / (((XMV/XMD)+XRT(:,:,:,1))* &
                        EXP( XALPW - XBETAW/ZT(:,:,:) - XGAMW*ALOG(ZT(:,:,:)) ))
!
! set the perturbation for Theta
!
        XTHT(:,:,:) = XTHT(:,:,:) + XAMPLITH * COS (XPI*0.5*ZDIST(:,:,:)) **2
!
! compute the new vapor pressure (stored in ZHU!)
!
        ZT(:,:,:) =XTHT(:,:,:)*(XPABST(:,:,:)/XP00) **(XRD/XCPD)
        ZHU(:,:,:) = ZHU(:,:,:)*                                             &
                     EXP( XALPW - XBETAW/ZT(:,:,:) - XGAMW*ALOG(ZT(:,:,:)) )
!
! set the perturbation for r_v such that the relative humidity is conserved
!
        XRT(:,:,:,1) = (XMV/XMD) * ZHU(:,:,:) / ( XPABST(:,:,:) - ZHU(:,:,:) )
      END WHERE
      CALL MPPDB_CHECK3D(ZT,"SET_PERTURB::ZT",PRECISION)
      CALL MPPDB_CHECK3D(ZHU,"SET_PERTURB::ZHU",PRECISION)
   ELSE
      WHERE ( ZDIST(:,:,:) <= 1.) 
        XTHT(:,:,:)  = XTHT(:,:,:)  + XAMPLITH * COS (XPI*0.5*ZDIST(:,:,:)) **2
        XRT(:,:,:,1) = XRT(:,:,:,1) + XAMPLIRV * COS (XPI*0.5*ZDIST(:,:,:)) **2 
      END WHERE
    END IF
    CALL MPPDB_CHECK3D(XRT(:,:,:,1),"SET_PERTURB::XRT",PRECISION)
    !CALL MPPDB_CHECK3D(XTHM,"SET_PERTURB::XTHM",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*	 3.     COMPUTE THE PERTURBATION ON THE WIND : 
!               --------------------------------------
!
  CASE('UV')
    IJ0   = IJU_ll/2
    IJMIN = IJ0 - NINT(IJU_ll*3./8.)
    IJMAX = IJ0 + NINT(IJU_ll*3./8.)
    ALLOCATE(ZPHI_ll(IIU_ll,IJU_ll))
    ALLOCATE(ZPU_ll(IIU_ll,IJU_ll))
    ALLOCATE(ZPV_ll(IIU_ll,IJU_ll))
!
!*       3.1   stream function at vorticity point :
!
    DO JI = 1,IIU_ll
      DO JJ = 1,IJU_ll
        ZPHI_ll(JI,JJ) = XAMPLIUV*EXP(-((ZYHAT_ll(JJ)-ZYHAT_ll(IJ0))/XRADY)**2)       &
                             * COS(2.*XPI/XRADX*ZXHAT_ll(JI))
      END DO
    END DO
!
!*       3.2   Wind perturbation :
!
    DO JI = 1,IIU_ll
      DO JJ = IJMIN,IJMAX
        ZPU_ll(JI,JJ) = (ZPHI_ll(JI,JJ+1)-ZPHI_ll(JI,JJ)) / (-ZYHAT_ll(JJ+1)+ZYHAT_ll(JJ) )
        ZPV_ll(JI,JJ) = (ZPHI_ll(JI+1,JJ)-ZPHI_ll(JI,JJ)) / ( ZXHAT_ll(JI+1)-ZXHAT_ll(JI) )
      END DO
    END DO
!
    DO JJ = 1,IJMIN-1
      ZPU_ll(:,JJ) = ZPU_ll(:,IJMIN)
      ZPV_ll(:,JJ) = ZPV_ll(:,IJMIN)
    END DO
    DO JJ = IJMAX+1,IJU_ll
      ZPU_ll(:,JJ) = ZPU_ll(:,IJMAX)
      ZPV_ll(:,JJ) = ZPV_ll(:,IJMAX)
    END DO
!
!*       3.3   Add to the U and V fields :
!
    ALLOCATE(ZPU(IIU,IJU),ZPV(IIU,IJU))
    CALL SCATTER(ZPU_ll,ZPU)
    CALL SCATTER(ZPV_ll,ZPV)
    DEALLOCATE(ZPU_ll,ZPV_ll)
    DO JI = 1,IIU
      DO JJ = 1,IJU
        XUT(JI,JJ,:) = XUT(JI,JJ,:) + ZPU(JI,JJ)
      END DO
    END DO
    DO JJ = 1,IJU
      DO JI = 1,IIU
        XVT(JI,JJ,:) = XVT(JI,JJ,:) + ZPV(JI,JJ)
      END DO
    END DO
    DEALLOCATE(ZPU,ZPV)
!
!
  CASE('WH','WW')  !    white noise is computed on global domain
                   ! J.Escobar optim => need only identical random on all domain 
!
 ALLOCATE(ZWHITE(IIU,IJU))
!
 DO JK = IKB, NKWH
     IKX = (NIMAX_ll+1)/2
     ZX  = 2*XPI/NIMAX_ll
     ALLOCATE(ZCX_ll(IIU_ll,IKX))
     ALLOCATE(ZSX_ll(IIU_ll,IKX))
     IKY = (NJMAX_ll+1)/2
     ALLOCATE(ZCY_ll(IJU_ll,IKY))
     ALLOCATE(ZSY_ll(IJU_ll,IKY))

!
    DO JX = 1, IKX
      DO JI = IIB_ll,IIE_ll
        ZCX_ll(JI,JX) = COS(ZX*JX*(JI-IIB_ll))
        ZSX_ll(JI,JX) = SIN(ZX*JX*(JI-IIB_ll))
      END DO
    END DO
      IKY = (NJMAX_ll+1)/2
      ZY  = 2*XPI/NJMAX_ll
      DO JY = 1, IKY
       DO JJ = IJB_ll,IJE_ll
          ZCY_ll(JJ,JY) = COS(ZY*JY*(JJ-IJB_ll))
          ZSY_ll(JJ,JY) = SIN(ZY*JY*(JJ-IJB_ll))
        END DO
      END DO
!
    ZWHITE(:,:) = 0.0
! 
    DO JY = 1, IKY
      DO JX = 1, IKX
        CALL RANDOM_NUMBER(ZRANDOM)
        ZT1 = 2*XPI*ZRANDOM
        ZAL = SIN(ZT1)
        ZBL = COS(ZT1)
        IF ( .NOT. L2D ) THEN
          CALL RANDOM_NUMBER(ZRANDOM)
          ZT2 = 2*XPI*ZRANDOM
          ZCL = SIN(ZT2)
          ZDL = COS(ZT2)
        ELSE
          ZCL=0.
          ZDL=0.
        END IF
!
        DO JJ = IJB,IJE
           JJ_ll = JJ  + IYOR-1
          DO JI = IIB,IIE
             JI_ll = JI + IXOR-1
            ZWHITE(JI,JJ) = ZWHITE(JI,JJ)+             &
               (ZBL+ZDL)*ZSX_ll(JI_ll,JX)*ZCY_ll(JJ_ll,JY) +       &
               (ZAL+ZCL)*ZCX_ll(JI_ll,JX)*ZCY_ll(JJ_ll,JY) +       &
               (ZCL-ZAL)*ZSX_ll(JI_ll,JX)*ZSY_ll(JJ_ll,JY) +       &
               (ZBL-ZDL)*ZCX_ll(JI_ll,JX)*ZSY_ll(JJ_ll,JY)
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(ZCX_ll,ZSX_ll,ZCY_ll,ZSY_ll)
!
    ZVAR= SUM_DD_R2_ll( (ZWHITE(IIB:IIE,IJB:IJE))**2 )/FLOAT(NIMAX_ll*NJMAX_ll)
    CALL MPPDB_CHECK2D(ZWHITE,"SET_PERTURB::ZWHITE",PRECISION)
    ZWHITE(:,:) = ZWHITE(:,:)/SQRT(ZVAR)
!
!!$    CALL GET_HALO(XTHM)
    CALL MPPDB_CHECK2D(XTHM(:,:,JK),"SET_PERTURB::XTHM",PRECISION)
    CALL MPPDB_CHECK2D(ZWHITE,"SET_PERTURB::ZWHITE",PRECISION)

    IF (CPERT_KIND=='WH') THEN ! white noise on theta
      XTHT(:,:,JK) = XTHT(:,:,JK) +  XAMPLIWH * ZWHITE(:,:)
      CALL MPPDB_CHECK2D(XTHT(:,:,JK),"SET_PERTURB::XTHT",PRECISION)      
    ELSE  ! white noise on wind
      XWT(:,:,JK) = XWT(:,:,JK) +  XAMPLIWH * ZWHITE(:,:)
      XUT(:,:,JK) = XUT(:,:,JK) +  XAMPLIWH * ZWHITE(:,:)
      XVT(:,:,JK) = XVT(:,:,JK) +  XAMPLIWH * ZWHITE(:,:)
    ENDIF
!
 END DO
!
! white noise for inflow/outflow U field in X direction                   
!
 IF (LWH_LBXU) THEN
   ALLOCATE(ZWHITE_ll(IIU_ll,IJU_ll))
   CALL GATHERALL_FIELD_ll('XY',ZWHITE,ZWHITE_ll,IRESP)
   DO JK=1,MIN(IKU,IIU_ll)
     DO JI=1,SIZE(XLBXUM,1)
      XLBXUM(JI,:,JK) = XLBXUM(JI,:,JK) + XAMPLIWH * ZWHITE_ll(JK,:)
     END DO
   END DO
   DEALLOCATE(ZWHITE_ll)
 END IF
!
! white noise for inflow/outflow V field in Y direction                   
!
 IF (LWH_LBYV) THEN
   ALLOCATE(ZWHITE_ll(IIU_ll,IJU_ll))
   CALL GATHERALL_FIELD_ll('XY',ZWHITE,ZWHITE_ll,IRESP)
   DO JK=1,MIN(IKU,IJU_ll)
     DO JJ=1,SIZE(XLBXVM,2)
      XLBXVM(:,JJ,JK) = XLBXVM(:,JJ,JK) + XAMPLIWH * ZWHITE_ll(JK,:)
     END DO
   END DO
   DEALLOCATE(ZWHITE_ll)
 END IF

 DEALLOCATE(ZWHITE)

 CALL GET_HALO(XTHT)
 CALL GET_HALO(XUT)
 CALL GET_HALO(XVT)
 CALL GET_HALO(XWT)

 CALL MPPDB_CHECK3D(XTHT,"SET_PERTURB::XTHT",PRECISION)
 CALL MPPDB_CHECK3D(XUT,"SET_PERTURB::XUT",PRECISION)
 CALL MPPDB_CHECK3D(XVT,"SET_PERTURB::XVT",PRECISION)
 CALL MPPDB_CHECK3D(XWT,"SET_PERTURB::XWT",PRECISION)
 CALL MPPDB_CHECK3D(XRT(:,:,:,1),"SET_PERTURB::XRT",PRECISION)

!-------------------------------------------------------------------------------
!
!*	 4.     LATERAL BOUNDARY CONDITIONS :
!               -----------------------------
!   	
!
  IF (CPERT_KIND=='WH') THEN ! white noise on theta
    IF (LWEST_ll() .AND. CLBCX(1)/='CYCL')  &
     XTHT(IIB-1,:,IKB) = 2. * XTHT(IIB,:,IKB) - XTHT(IIB+1,:,IKB)
    IF (LEAST_ll() .AND. CLBCX(1)/='CYCL')   &
     XTHT(IIE+1,:,IKB) = 2. * XTHT(IIE,:,IKB) - XTHT(IIE-1,:,IKB)
    IF (LSOUTH_ll() .AND. CLBCY(1)/='CYCL') &
     XTHT(:,IJB-1,IKB) = 2. * XTHT(:,IJB,IKB) - XTHT(:,IJB+1,IKB)
    IF (LNORTH_ll() .AND. CLBCY(1)/='CYCL') &
     XTHT(:,IJE+1,IKB) = 2. * XTHT(:,IJE,IKB) - XTHT(:,IJE-1,IKB)
  ELSE ! white noise on wind
    IF (LWEST_ll() .AND. CLBCX(1)/='CYCL')  &
     XWT(IIB-1,:,IKB) = 2. * XWT(IIB,:,IKB) - XWT(IIB+1,:,IKB)
    IF (LEAST_ll() .AND. CLBCX(1)/='CYCL')   &
     XWT(IIE+1,:,IKB) = 2. * XWT(IIE,:,IKB) - XWT(IIE-1,:,IKB)    
    IF (LSOUTH_ll() .AND. CLBCY(1)/='CYCL') &
     XWT(:,IJB-1,IKB) = 2. * XWT(:,IJB,IKB) - XWT(:,IJB+1,IKB)
    IF (LNORTH_ll() .AND. CLBCY(1)/='CYCL') &
     XWT(:,IJE+1,IKB) = 2. * XWT(:,IJE,IKB) - XWT(:,IJE-1,IKB)
!     
    IF (LWEST_ll() .AND. CLBCX(1)/='CYCL')  &
     XUT(IIB-1,:,IKB) = 2. * XUT(IIB,:,IKB) - XUT(IIB+1,:,IKB)
    IF (LEAST_ll() .AND. CLBCX(1)/='CYCL')   &
     XUT(IIE+1,:,IKB) = 2. * XUT(IIE,:,IKB) - XUT(IIE-1,:,IKB)
    IF (LSOUTH_ll() .AND. CLBCY(1)/='CYCL') &
     XUT(:,IJB-1,IKB) = 2. * XUT(:,IJB,IKB) - XUT(:,IJB+1,IKB)
    IF (LNORTH_ll() .AND. CLBCY(1)/='CYCL') &
     XUT(:,IJE+1,IKB) = 2. * XUT(:,IJE,IKB) - XUT(:,IJE-1,IKB)
!
    IF (LWEST_ll() .AND. CLBCX(1)/='CYCL')  &
     XVT(IIB-1,:,IKB) = 2. * XVT(IIB,:,IKB) - XVT(IIB+1,:,IKB)
    IF (LEAST_ll() .AND. CLBCX(1)/='CYCL')   &
     XVT(IIE+1,:,IKB) = 2. * XVT(IIE,:,IKB) - XVT(IIE-1,:,IKB)
    IF (LSOUTH_ll() .AND. CLBCY(1)/='CYCL') &
     XVT(:,IJB-1,IKB) = 2. * XVT(:,IJB,IKB) - XVT(:,IJB+1,IKB)
    IF (LNORTH_ll() .AND. CLBCY(1)/='CYCL') &
     XVT(:,IJE+1,IKB) = 2. * XVT(:,IJE,IKB) - XVT(:,IJE-1,IKB)
  ENDIF
!
!
  CASE('SH')  !    Shock (Burger's Equation)
!
    ZOMEGA = 2.0*XPI/FLOAT(IIE-IIB)
    DO JI = IIB, IIE
      XUT(JI,:,:) = XUT(JI,:,:) + XAMPLIUV*SIN( ZOMEGA*FLOAT(JI-IIB) )
    END DO
    XVT(:,:,:) = 0.0
    XWT(:,:,:) = 0.0
!
!
END SELECT
!
DEALLOCATE(ZXHAT_ll,ZYHAT_ll)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_PERTURB 
