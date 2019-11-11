!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_EDDYUV_FLUX_n
!     #########################
!
INTERFACE
!
      SUBROUTINE EDDYUV_FLUX_n (KMI,KTCOUNT,PVM,PTHM,PRHODJ,PRHODREF,PPABSM,PRVS,PVU_FLUX_M)
!
!
INTEGER,               INTENT(IN)    :: KMI   ! Model index
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! iteration count
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PVM
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ            ! dry density *J
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF          ! dry density 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PVU_FLUX_M

!
END SUBROUTINE EDDYUV_FLUX_n

!
END INTERFACE
!
END MODULE MODI_EDDYUV_FLUX_n

!     #############################################################################################
      SUBROUTINE EDDYUV_FLUX_n (KMI,KTCOUNT,PVM,PTHM,PRHODJ,PRHODREF,PPABSM,PRVS,PVU_FLUX_M)
!     #############################################################################################
!!!
!!
!!    PURPOSE
!!    -------
!!    Barotropic fluxes parameterization (v'u') for 2D transect (latitude, altitude) model
!!
!!**  METHOD
!!    ------ 
!!    The points in the domain where barotropic instability is reached are first searched for using the
!!    criterion that meridional gradient of absolute vorticity be negative. 
!!    The test is then made to check if these consist of large areas or isolated points.
!!    If large area are found, a diffusion of angular momentum is applied using a K gradient formulation, 
!!    v'M'=Kyy * (dM/dt). The main purpose of the routine is the determination of the
!!    Kyy coefficient from DM/DY. 
!!    A linear relationship is assumed between the integral of the instability and the
!!    diffusion coefficient. Ultimately, the tendency is derived for the zonal
!!    momentum
!! 
!!
!!    REFERENCE
!!     ---------
!!    Luz and Hourdin (2001) : 
!!    Peyrillé et al. (2007) : 
!!   
!!    AUTHOR
!!    ------
!!	  P.Peyrille          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  18/02/04
!!      /05/12 M.Tomasini Grid-nesting
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O

USE MODE_ll
!
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_BUDGET
!
USE MODI_SHUMAN
USE MODI_BUDGET
USE MODD_CST
!
USE MODD_DIM_n
USE MODD_CONF
USE MODD_CONF_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_METRICS_n
!USE MODD_TIME
!USE MODD_TIME_n
!USE MODD_DYN_n
USE MODD_CURVCOR_n
USE MODI_GRADIENT_M
USE MODI_GRADIENT_W
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_SHUMAN
USE MODE_GRIDPROJ
USE MODI_MEAN_Z
USE MODD_LUNIT_n
USE MODD_LATZ_EDFLX
!
USE MODE_MSG
!
IMPLICIT NONE
!       0.01 Arguments declarations
!
INTEGER,                  INTENT(IN)    :: KMI   ! Model index
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! iteration count
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ            ! dry density of
!                                 ! anelastic reference state * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF          ! dry density 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABSM
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PVM
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHM
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS

REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PVU_FLUX_M

!*       0.02   Declarations of local variables :
!
INTEGER:: IIB,IJB,IKB        ! Begining useful area  in x,y directions
INTEGER:: IIE,IJE,IKE        ! End useful area in x,y directions
INTEGER:: JI,JK,JJ
INTEGER:: IIU,IKU,IJU
!
! 3D array x,y,z
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZCORIOZ,ZVOZ,ZABVOR
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZABVOR_DX,ZKVU
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZDIV_UV,ZP,ZAPV,ZLAT3D
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZUV_FLUX
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZMT,ZMT_DX,ZDAPV_DX
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZINT_SB
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZMT_FLUX ! Angular momentum flux
!
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZDX3D
!
! 2D arrays
INTEGER, DIMENSION(:,:), ALLOCATABLE :: J_NEG ! array containing index of negative pv gradient
REAL, DIMENSION(:,:), ALLOCATABLE    :: ZK,ZMODUL
!
! 1D array in x 
INTEGER , DIMENSION(:), ALLOCATABLE  :: JI_MIN,JI_MAX
REAL, DIMENSION(:), ALLOCATABLE      :: ZSB
REAL, DIMENSION(:), ALLOCATABLE      :: ZALPHA
!
REAL                                 :: ZB,ZA
REAL                                 :: ZFHMW,ZPHI0
REAL                                 :: JNEG
REAL                                 :: JNEG2
INTEGER                              :: JIDIFF
INTEGER                              :: JP
INTEGER                              :: JNI ! Nb of point where neg Pv gradiebnt  depiected 
INTEGER                              :: JNB ! nb of areas of pv gradioent neg
INTEGER                              :: JPI  ! loop index
INTEGER                              :: IRESP
CHARACTER(LEN=100)                   :: YMSG
!
!-------------------------------------------------------------------------
!
!*       1.     COMPUTES THE DOMAIN DIMENSIONS
!               ------------------------------
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)

IKU=SIZE(XZHAT) ! nb points sur Z
IKB=1+JPVEXT
IKE=IKU-JPVEXT

IIU = SIZE(XXHAT)
IJU = SIZE(XYHAT)

! allocation of temporary array
! ---------------------------------------------------------
! 1D arrays X direction
ALLOCATE(JI_MIN(IIU))     ; JI_MIN=0
ALLOCATE(JI_MAX(IIU))     ; JI_MAX=0
ALLOCATE(ZSB(IIU))        ; ZSB=0.0
ALLOCATE(ZALPHA(IIU)) ; ZALPHA= 0.0
! 2D Arrays X,Z
ALLOCATE(J_NEG(IIU,IKU))  ; J_NEG=0
! 2D Arrays (I, Nb instab zones at max) 
ALLOCATE(ZMODUL(IIU,IIU)) ; ZMODUL=0.0
ALLOCATE(ZK(IIU,IIU))     ; ZK=0.0
!
! 3D arrays
ALLOCATE(ZUV_FLUX(IIU,IJU,IKU )) ; ZUV_FLUX=0.0
ALLOCATE(ZCORIOZ(IIU,IJU,IKU)) ; ZCORIOZ=0.0
ALLOCATE(ZVOZ(IIU,IJU,IKU))    ; ZVOZ=0.0
ALLOCATE(ZABVOR(IIU,IJU,IKU))  ; ZABVOR=0.0
!
ALLOCATE(ZKVU(IIU,IJU,IKU )) ; ZKVU=0.0
ALLOCATE(ZABVOR_DX(IIU,IJU,IKU )) ; ZABVOR_DX=0.0
ALLOCATE(ZDIV_UV(IIU,IJU,IKU ))  ; ZDIV_UV=0.0
ALLOCATE(ZP(IIU,IJU,IKU ))       ; ZP=0.0
ALLOCATE(ZAPV(IIU,IJU,IKU ))     ; ZAPV=0.0
ALLOCATE(ZLAT3D(IIU,IJU,IKU ))   ; ZLAT3D=0.0
ALLOCATE(ZMT(IIU,IJU,IKU ))      ; ZMT=0.0
ALLOCATE(ZMT_DX(IIU,IJU,IKU ))   ; ZMT_DX=0.0
ALLOCATE(ZDAPV_DX(IIU,IJU,IKU )) ; ZDAPV_DX=0.0
ALLOCATE(ZINT_SB(IIU,IJU,IKU ))     ; ZINT_SB=0.0
ALLOCATE(ZMT_FLUX(IIU,IJU,IKU ))     ; ZMT_FLUX=0.0
!
ALLOCATE(ZDX3D(IIU,IJU,IKU ))     ; ZDX3D=0.0

! -------------------------------------------------------------------------
!           2.0  FIRST COMPUTATIONS
! -------------------------------------------------------------------------

!Coriolis parmeter and relative vorticity
ZCORIOZ(:,:,:)= SPREAD(XCORIOZ(:,:),3,IKU)
ZLAT3D(:,:,:) = SPREAD(XLAT(:,:),3,IKU)
! 
! relative vorticity
ZVOZ(:,:,:)=GX_V_UV(1,IKU,1,PVM,XDXX,XDZZ,XDZX)
ZVOZ(:,:,2)=ZVOZ(:,:,3)
ZVOZ(:,:,1)=ZVOZ(:,:,3)
!
! abs. vort -> mass point
ZABVOR(:,:,:)=MXF(ZVOZ(:,:,:)) + ZCORIOZ(:,:,:)
! deabvor_dx --> flux pt
ZABVOR_DX(:,:,:)= GX_M_U(1,IKU,1,ZABVOR,XDXX,XDZZ,XDZX)

! Pressure divided by fluid density
ZP(:,:,:) = PPABSM(:,:,:) / PRHODREF(:,:,:)                         
! 
! ABS potential vorticity
ZAPV(:,:,:) = ZABVOR(:,:,:)/ ZP(:,:,:)
ZDAPV_DX(:,:,:)= GX_M_U(1,IKU,1,ZAPV,XDXX,XDZZ,XDZX)

! APV --> flux pt
ZAPV(:,:,:) = MXM(ZAPV(:,:,:))
!
! Zonal abs momentum
ZMT(:,:,:) = ( XRADIUS * COS( (XPI/180.0)*ZLAT3D(:,:,:) )* (- PVM(:,:,:)) ) +  & 
   ( XOMEGA *((XRADIUS)**(2))*((COS( (XPI/180.0)*ZLAT3D(:,:,:) ))**(2)) )
! Dmt/DX
ZMT_DX(:,:,:) = GX_M_U(1,IKU,1,ZMT,XDXX,XDZZ,XDZX)

! ! --------------------------------------------------------------------------
! 
!
!
!      3.0  Beginning of loop for each level
!       -------------------------------
 J_NEG(:,:)=0  ! array containing index of negative pv gradient

! Loop over vertical levels to depict neg. pv graident
DO JK=IKB,IKE
!                         
 JNI=0
 ZALPHA(:)=0.0
 ZK(:,:)=0.0
 !
 !      3.1 Find  points where there is barotropic instability
 !            i.e. zdapv_dx < 0.
 !                 -----------------------------------
 !
  DO JI=IIB,IIE
   IF  (ZDAPV_DX(JI,2,JK) < 0.0) THEN
    JNI=JNI+1       
    J_NEG(JNI,JK) = JI
   END IF
  END DO
 !
 ! JNEG=0 if stable
 !     = JI if unstable
 ! 
 !       3.2 For each area, find first and last points
 !       ------------------------------------------------------
 !
 IF (JNI > 1) THEN
! ==========================================================
  JNB=0 ! nb de zone a gradients negatifs
  JI_MIN(:)=IIE
  JI_MAX(:)=IIB
  ZSB(:) = 0.0
  JNEG=0
  JNEG2=0

  ! ---------------------------------------------------------
   DO JPI=1,JNI-2
   ! If there is only 1 pt, we do not do the calculations
    JNEG = J_NEG(JPI+1,JK)-J_NEG(JPI,JK)
    JNEG2 = J_NEG(JPI+2,JK)-J_NEG(JPI+1,JK)
    !
    ! 1) Case of an encoutered area : JNEG=1 (two consecutives points)
    ! a. 1st area
     IF ( (JNEG==1).AND.(JNB == 0) )THEN
          JNB=JNB+1
          JI_MIN(JNB) = MIN(JI_MIN(JNB),J_NEG(JPI,JK))
          JI_MAX(JNB) = MAX(JI_MAX(JNB),J_NEG((JPI+1),JK))
  !   b) nth area 
     ELSE IF ( (JNEG==1).AND.(JNB /= 0) )THEN
          ! Same area JNB is not incremented
          JI_MIN(JNB) = MIN(JI_MIN(JNB),J_NEG(JPI,JK))
          JI_MAX(JNB) =  J_NEG((JPI+1),JK)
     ELSE IF ((JNEG>1) .AND. (JNEG2==1))  THEN 
          ! New area, JNB is incremented
          JNB=JNB+1
          JI_MIN(JNB) = J_NEG(JPI+1,JK)
          JI_MAX(JNB) = J_NEG((JPI+2),JK)
          !
     END IF
    END DO  

! 
!       3.3 Compute diffusion coefficient for each points
!           ---------------------------------------------------
!
  ZDX3D=SPREAD(SPREAD(XDXHAT(:),2,IJU),3,IKU)
  ZINT_SB(:,:,:) = ZDAPV_DX(:,:,:) * ZDX3D(:,:,:)


  ! ---------------------------------------------
  ! PhD formulation
  ! ---------------------------------------------

  DO JI=IIB,IIE
    DO JP=1,JNB
  !     Test  of JI_MIN_IMAX value
     IF ((JI_MAX(JP) > IIE) .OR. (JI_MAX(JP)<JI_MIN(JP)) & 
                           .OR. (JI_MAX(JP)<IIB)) THEN
       CALL PRINT_MSG(NVERB_FATAL,'GEN','EDDYUV_FLUX_n','wrong value of JI_MAX')
     END IF
!
     IF ((JI_MIN(JP) > IIE) .OR. (JI_MIN(JP)<IIB)) THEN
       CALL PRINT_MSG(NVERB_FATAL,'GEN','EDDYUV_FLUX_n','wrong value of JI_MIN')
     END IF

     JIDIFF=JI_MAX(JP) - JI_MIN(JP)
     IF ( MOD(JIDIFF,2)== 0 ) THEN 
!
       ZFHMW  = 0.5 * ( ( 0.5 * (XLAT(JI_MAX(JP),2) + XLAT((JI_MAX(JP)-1),2)) )  - &
                 ( 0.5 * (XLAT(JI_MIN(JP),2) + XLAT((JI_MIN(JP)-1),2)) ) ) 

       ZPHI0  = 0.5 * ( ( 0.5 * (XLAT(JI_MAX(JP),2) + XLAT((JI_MAX(JP)-1),2)) )  + &
                 ( 0.5 * (XLAT(JI_MIN(JP),2) + XLAT((JI_MIN(JP)-1),2)) ) )
!
     ELSE
       ZFHMW  = 0.5 * (XLAT(JI_MAX(JP),2) - XLAT(JI_MIN(JP),2) )
       ZPHI0  = 0.5 * (XLAT(JI_MAX(JP),2) + XLAT(JI_MIN(JP),2) ) 
     END IF 
!
  IF ( ZFHMW /= 0.0 ) THEN
  ! Exponential fction
  ZMODUL(JI,1:JNB)= & 
  EXP(- ( (( (0.5*(XLAT(JI,IJB)+XLAT(JI-1,IJB))) - ZPHI0 )/(2.*ZFHMW) )**2) )
  END IF
!
  ZB = XUV_FLX1
  ZA = XUV_FLX2
!
! calcul SB pour les points correpondants à des zones d'instazbilite
! Si les points sont ok 
! Integration sur Z de ZINT_SB

  ZSB(JP)= SUM(ZINT_SB(JI_MIN(JP):JI_MAX(JP),IJB,JK) )
!
  ZALPHA(JP) =  ( (-ZB) * ZSB(JP)) + ZA
  ZK(JI,JP)  = ZALPHA(JP) * ZMODUL(JI,JP)
  END DO

  ! Diffusion coeff= horizontal integral of Zk
  ZKVU(JI,:,JK)= SUM(ZK(JI,1:JNB))

!
 ! Control of ZKVU value
 IF (ZKVU(JI,IJB,JK) < 0.0 ) THEN
   WRITE(YMSG,*) 'ZKVU(',JI,',',IJB,',',JK,') < 0.0'
   CALL PRINT_MSG(NVERB_FATAL,'GEN','EDDYUV_FLUX_n',YMSG)
 ENDIF           

 END DO ! end of loot JI
 
 ! IF NO ZONE OF INSTABLITY IS DETECTED
 ELSE
  ZKVU(:,:,JK) =0.0
 END IF  

!
!       ENd of vertical loop
!       --------------------
END DO
!       --------------------
!
!
!     3.4 Compute  eddy momentum u'v' flux from abs momentum flux M'u'
!  ---------------------------------------------
!
ZMT_FLUX(:,:,:) =   - ZKVU(:,:,:)*ZMT_DX(:,:,:)
!
ZMT_FLUX(IIB,:,:)  = ZMT_FLUX(IIB+1,:,:)
ZMT_FLUX(1,:,:)    = ZMT_FLUX(IIB,:,:)
ZMT_FLUX(IIE,:,:)  = ZMT_FLUX(IIE-1,:,:) 
ZMT_FLUX(IIU,:,:)  = ZMT_FLUX(IIE,:,:) 

ZUV_FLUX(:,:,:) = ( - ZMT_FLUX(:,:,:)) / & 
       (XRADIUS *COS( (XPI/180.0)*MXM(ZLAT3D(:,:,:)) )) 
!
!
!      4.0 BOUNDARY CONDITIONS        
!     --------------------------

ZUV_FLUX(IIB,:,:)  = ZUV_FLUX(IIB+1,:,:)
ZUV_FLUX(1,:,:)    = ZUV_FLUX(IIB,:,:)
ZUV_FLUX(IIE,:,:)  = ZUV_FLUX(IIE-1,:,:) 
ZUV_FLUX(IIU,:,:)  = ZUV_FLUX(IIE,:,:) 

!       2.1  Lissage temporel et stockage flux(t-dt)
!            --------------------------------------

     IF (KTCOUNT > 1) THEN
             
ZUV_FLUX(:,:,:) =  0.5*(PVU_FLUX_M(:,:,:) + ZUV_FLUX(:,:,:))

PVU_FLUX_M(:,:,:) = ZUV_FLUX(:,:,:)
     ELSE
                        
PVU_FLUX_M(:,:,:) = ZUV_FLUX(:,:,:)
      ENDIF
! 
! ---------------------------------------------------------------------        
!       3.   INTEGRATION IN THE SOURCE OF V
!             -----------------------------
!
! Take the divergence of the momentum flux 
ZDIV_UV(:,:,:) = GX_U_M(1,IKU,1,ZUV_FLUX,XDXX,XDZZ,XDZX)

! Lateral boundary conditions
ZDIV_UV(IIB,:,:)=0.0
ZDIV_UV(1,:,:)=0.0
ZDIV_UV(IIE,:,:)=0.0
ZDIV_UV(IIU,:,:)=0.0
!
PRVS(:,:,:) =  PRVS(:,:,:) - PRHODJ(:,:,:) * ZDIV_UV(:,:,:)
!
! DEALLOCATE
DEALLOCATE(ZCORIOZ)
DEALLOCATE(ZVOZ)
DEALLOCATE(ZABVOR)
DEALLOCATE(J_NEG)
DEALLOCATE(JI_MIN)
DEALLOCATE(JI_MAX)
DEALLOCATE(ZSB)
DEALLOCATE(ZUV_FLUX)
DEALLOCATE(ZKVU)
DEALLOCATE(ZABVOR_DX)
DEALLOCATE(ZDIV_UV)
DEALLOCATE(ZP)
DEALLOCATE(ZAPV)
DEALLOCATE(ZLAT3D)
DEALLOCATE(ZMT)
DEALLOCATE(ZMT_DX)
DEALLOCATE(ZDAPV_DX)
DEALLOCATE(ZINT_SB)
DEALLOCATE(ZMT_FLUX)
DEALLOCATE(ZALPHA)
!
DEALLOCATE(ZK)
DEALLOCATE(ZMODUL)
DEALLOCATE(ZDX3D)


END SUBROUTINE EDDYUV_FLUX_n
