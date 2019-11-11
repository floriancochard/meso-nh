!MNH_LIC Copyright 1997-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SPAWN_LS_n
!     ####################
!
INTERFACE 
!
      SUBROUTINE SPAWN_LS_n(KDAD, PTSTEP, KMI,                       &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                    KDXRATIO,KDYRATIO,         &
                    HLBCX,HLBCY,PZZ,PZHAT,OSLEVE,PLEN1,PLEN2,        &
                    PCOEFLIN_LBXM,                                   &
                                  PLSTHM,PLSRVM,                     &
                                  PLSUM,PLSVM,PLSWM,PLSZWSM,         &
                                  PLSTHS,PLSRVS,                     &
                                  PLSUS,PLSVS,PLSWS,PLSZWSS          )
!
INTEGER,   INTENT(IN)  :: KDAD      ! number of the DAD model
REAL,             INTENT(IN)    :: PTSTEP   !  Time step
INTEGER,          INTENT(IN)    :: KMI      ! model number
                                    ! interpolation coefficients 
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ   ! small scale vertical grid
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! height of the Gal-Chen levels
LOGICAL,                INTENT(IN) :: OSLEVE! flag for Sleve coordinate
REAL,                   INTENT(IN) :: PLEN1 ! Decay scale for smooth topography
REAL,                   INTENT(IN) :: PLEN2 ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:,:), INTENT(IN) :: PCOEFLIN_LBXM ! coefficient used for
                                                    ! vertical interpolation
!  
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLSTHM,PLSRVM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLSUM,PLSVM,PLSWM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)    :: PLSZWSM          ! Large Scale fields at t-dt

REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PLSTHS,PLSRVS ! Large Scale source terms 
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PLSUS,PLSVS,PLSWS ! Large Scale source terms 
REAL, DIMENSION(:,:),   INTENT(OUT)   :: PLSZWSS
!
END SUBROUTINE SPAWN_LS_n
!
END INTERFACE
!
END MODULE MODI_SPAWN_LS_n
!
!
!     ################################################################
      SUBROUTINE SPAWN_LS_n (KDAD,PTSTEP,KMI,                        &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                    KDXRATIO,KDYRATIO,                               &
                    HLBCX,HLBCY,PZZ,PZHAT,OSLEVE,PLEN1,PLEN2,        &
                    PCOEFLIN_LBXM,                                   &
                                  PLSTHM,PLSRVM,                     &
                                  PLSUM,PLSVM,PLSWM,PLSZWSM,         &
                                  PLSTHS,PLSRVS,                     &
                                  PLSUS,PLSVS,PLSWS,PLSZWSS          )
!     ################################################################
!
!!****  *SPAWN_LS_n* - Refresh of the Large Scale sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of SPAWN_CPL$   is to 'refresh' Large Scale sources
!!    of all the LS fields of the current nested model when the current time
!!    step corresponds to a coupling instant common to all the nested models.
!!    The informations to be interpolated on the fine mesh grid of model _n is
!!    taken from is taken from the corresponding LS fields of its outer (DAD) 
!!    model DAD(_n). 
!
!
!!**  METHOD
!!    ------
!!      The basic task consists in interpolating fields from outer model _n
!!    to present inner model, using horizontal Bikhardt interpolation. The 
!!    DAD(_n) fields are passed by the module MODD_FIELD_n and the _n fields 
!!    are going out from the subroutine by argument.  
!!      We first interpolate the future insatnt and then, we compute the LS
!!    tendency for model _n.
!!
!!    EXTERNAL
!!    --------
!!
!!      Subroutine BIKHARDT : performs the horizontal interpolation
!!
!!      Subroutine COEF_VER_INTERP_LIN  : computes the coefficients for vertical
!!                                        interpolation
!!      Function   VER_INTERP_LIN  : performs the vertical interpolation
!!                                        interpolation
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_DYN: XTSTEP_MODEL1, NCPL_TIMES
!!
!!      Module MODD_FIELD_n: XLSUM,XLSVM,XLSWM,XLSUS,XLSVS,XLSWS,
!!                           XLSTHM,XLSRVM,XLSTHS,XLSRVS
!!
!!      Module MODD_GRID_n : XZS,XZZ
!!
!!      Module MODD_PARAMETERS : JPHEXT,JPVEXT
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    J. Stein and J. P. Lafore  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     22/12/97
!!    P. Jabouille   19/04/00 parallelisation
!!      J.Escobar : 18/12/2015 : Correction of bug in bound in // for NHALO <>1 
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 20/03/2019: fixes: wrong order of the dummy arguments + double deallocate
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_PARAMETERS
USE MODD_DYN 
USE MODD_LSFIELD_n
USE MODD_FIELD_n      ! modules relative to the outer model _n
USE MODD_GRID_n   
!
USE MODI_BIKHARDT
USE MODI_COEF_VER_INTERP_LIN
USE MODI_VER_INTERP_LIN
USE MODI_SHUMAN
USE MODI_VERT_COORD
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
INTEGER,          INTENT(IN)    :: KDAD     ! number of the DAD model
REAL,             INTENT(IN)    :: PTSTEP   !  Time step
INTEGER,          INTENT(IN)    :: KMI      ! model number
                                    ! interpolation coefficients
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ   ! small scale vertical grid
REAL, DIMENSION(:),     INTENT(IN) :: PZHAT ! height of the Gal-Chen levels
LOGICAL,                INTENT(IN) :: OSLEVE! flag for Sleve coordinate
REAL,                   INTENT(IN) :: PLEN1 ! Decay scale for smooth topography
REAL,                   INTENT(IN) :: PLEN2 ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:,:), INTENT(IN) :: PCOEFLIN_LBXM ! coefficient used for
                                                    ! vertical interpolation
!
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLSTHM,PLSRVM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLSUM,PLSVM,PLSWM ! Large Scale fields at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)    :: PLSZWSM          ! Large Scale fields at t-dt
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PLSTHS,PLSRVS ! Large Scale source terms
REAL, DIMENSION(:,:,:), INTENT(OUT)   :: PLSUS,PLSVS,PLSWS ! Large Scale source terms
REAL, DIMENSION(:,:),   INTENT(OUT)   :: PLSZWSS
!
!
!*       0.2   declarations of local variables
!
REAL                   :: ZTIME                   ! Interpolation length
INTEGER                :: IIB,IJB,IIE,IJE,IIU,IJU,IKU
LOGICAL                :: GVERT_INTERP
INTEGER, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: IKLIN
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZCOEFLIN
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2))                         :: ZZSLS
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2))                         :: ZZSMTLS
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZZSS
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZZLS1
                                    ! Gal Chen grid
                                    ! with the smooth LS orography on the SS grid
                                    ! at the W location then at the U location
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZZLS2
                                    ! Gal Chen grid
                                    ! with the smooth LS orography on the SS grid
                                    ! at the mass location then at the V location
!
TYPE(LIST_ll), POINTER :: TZLSFIELD_ll   ! list of metric coefficient fields
! Variables used for LS communications
INTEGER :: IINFO_ll, IDIMX, IDIMY
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTLSUM, ZTLSVM, ZTLSWM, ZTLSTHM, ZTLSRVM
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTZS,ZZS
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTZWS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTZSMT,ZZSMT
REAL, DIMENSION(:,:,:), ALLOCATABLE :: Z1,Z2,Z3,Z4,Z5
REAL, DIMENSION(:,:),   ALLOCATABLE :: Z6
!
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
GVERT_INTERP = .TRUE.
ZTIME = (NCPL_TIMES(NCPL_CUR,1) - NCPL_TIMES(NCPL_CUR-1,1)) * XTSTEP_MODEL1
!
ALLOCATE(ZTLSUM(IDIMX,IDIMY,SIZE(PLSUM,3)))
ALLOCATE(ZTLSVM(IDIMX,IDIMY,SIZE(PLSVM,3)))
ALLOCATE(ZTLSWM(IDIMX,IDIMY,SIZE(PLSWM,3)))
ALLOCATE(ZTLSTHM(IDIMX,IDIMY,SIZE(PLSTHM,3)))
ALLOCATE(ZTLSRVM(IDIMX,IDIMY,SIZE(PLSRVM,3)))
ALLOCATE(ZTZWS(IDIMX,IDIMY))
!
IF(GVERT_INTERP) THEN
  ALLOCATE(ZTZS  (IDIMX,IDIMY,1))
  ALLOCATE(ZTZSMT(IDIMX,IDIMY,1))
  ALLOCATE(ZZS(SIZE(XZS,1),SIZE(XZS,2),1))
  ZZS(:,:,1)=XZS(:,:)
  ALLOCATE(ZZSMT(SIZE(XZS,1),SIZE(XZS,2),1))
  ZZSMT(:,:,1)=XZSMT(:,:)
END IF
!
ALLOCATE(Z1(SIZE(XLSUM,1),SIZE(XLSUM,2),SIZE(XLSUM,3)))
ALLOCATE(Z2(SIZE(XLSUM,1),SIZE(XLSUM,2),SIZE(XLSUM,3)))
ALLOCATE(Z3(SIZE(XLSUM,1),SIZE(XLSUM,2),SIZE(XLSUM,3)))
ALLOCATE(Z4(SIZE(XLSUM,1),SIZE(XLSUM,2),SIZE(XLSUM,3)))
ALLOCATE(Z5(SIZE(XLSUM,1),SIZE(XLSUM,2),SIZE(XLSUM,3)))
ALLOCATE(Z6(SIZE(XLSUM,1),SIZE(XLSUM,2)))
!
Z1=XLSUM+XLSUS*ZTIME
CALL SET_LSFIELD_1WAY_ll(Z1, ZTLSUM, KMI)
Z2=XLSVM+XLSVS*ZTIME
CALL SET_LSFIELD_1WAY_ll(Z2, ZTLSVM, KMI)
Z3=XLSWM+XLSVS*ZTIME
CALL SET_LSFIELD_1WAY_ll(Z3, ZTLSWM, KMI)
Z4=XLSTHM+XLSTHS*ZTIME
CALL SET_LSFIELD_1WAY_ll(Z4, ZTLSTHM, KMI)
IF ( SIZE(PLSRVM,1) /= 0 ) THEN
  Z5=XLSRVM+XLSRVS*ZTIME
  CALL SET_LSFIELD_1WAY_ll(Z5, ZTLSRVM, KMI)
ENDIF
IF ( SIZE(PLSZWSM,1) /= 0 ) THEN
  Z6=XLSZWSM+XLSZWSS*ZTIME
  CALL SET_LSFIELD_1WAY_ll(Z6, ZTZWS, KMI)
ENDIF
!
IF ( GVERT_INTERP ) THEN
   CALL SET_LSFIELD_1WAY_ll(ZZS,   ZTZS,   KMI)
   CALL SET_LSFIELD_1WAY_ll(ZZSMT, ZTZSMT, KMI)
ENDIF
!
!        1.4  Communication
!
CALL LS_FORCING_ll(KMI, IINFO_ll)
!
DEALLOCATE(Z1,Z2,Z3,Z4,Z5,Z6)
!
!        1.5  Back to the (current) child model
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
CALL UNSET_LSFIELD_1WAY_ll()
!
IF ( GVERT_INTERP ) THEN
  !
  ! the orography of the DAD model is interpolated on the fine mesh model
  !
  ZZSLS=0.
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTZS(:,:,1),ZZSLS(IIB:IIE,IJB:IJE))
  !
  ZZSMTLS=0.
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTZSMT(:,:,1),ZZSMTLS(IIB:IIE,IJB:IJE))
  !
  CALL VERT_COORD(OSLEVE,ZZSLS,ZZSMTLS,PLEN1,PLEN2,PZHAT,ZZLS1)
  !
END IF
!
! duration between the two coupling instants
!
ZTIME = (NCPL_TIMES(NCPL_CUR,1) - NCPL_TIMES(NCPL_CUR-1,1)) * XTSTEP_MODEL1
!
!
!-------------------------------------------------------------------------------
!
!*      2.  COMPUTE LARGE SCALE FIELDS ON THE W GRID
!           ----------------------------------------
!
!
!*      2.1  Horizontal Bikhardt interpolation on the W grid
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,4,         &
               HLBCX,HLBCY,ZTLSWM,PLSWS(IIB:IIE,IJB:IJE,:))
!
!*      2.2  Vertical linear interpolation on the W grid
!
!
IF ( GVERT_INTERP ) THEN
  !
  CALL COEF_VER_INTERP_LIN(ZZLS1,PZZ,IKLIN,ZCOEFLIN)
  !
  PLSWS(:,:,:)=VER_INTERP_LIN(PLSWS,IKLIN,ZCOEFLIN)
  !
END IF
!
! compute the LS tendencies
!
PLSWS(:,:,:)   = (PLSWS(:,:,:) - PLSWM(:,:,:)) / ZTIME
!
!------------------------------------------------------------------------------
!
!
!*      3.   COMPUTE LARGE SCALE FIELDS ON THE MASS GRID
!            -------------------------------------------
!
!*      3.1  Horizontal Bikhardt interpolation on the mass grid
!
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
               HLBCX,HLBCY,ZTLSTHM,PLSTHS(IIB:IIE,IJB:IJE,:))

!
IF ( SIZE(PLSRVM,1) /= 0 ) THEN
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTLSRVM,PLSRVS(IIB:IIE,IJB:IJE,:))
END IF
IF ( SIZE(PLSZWSM,1) /= 0 ) THEN
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTZWS,PLSZWSS(IIB:IIE,IJB:IJE))
END IF
!
!*      3.2  Vertical linear interpolation on the mass grid
!
!
IF ( GVERT_INTERP ) THEN
  IKU = SIZE(PZZ,3)
  !
  ZZLS2=MZF(1,IKU,1,ZZLS1)
  ZZLS2(:,:,IKU)=2.*ZZLS2(:,:,IKU-1)-ZZLS2(:,:,IKU-2)
  ZZSS=MZF(1,IKU,1,PZZ)
  ZZSS(:,:,IKU)=2.*ZZSS(:,:,IKU-1)-ZZSS(:,:,IKU-2)
  !
  CALL COEF_VER_INTERP_LIN(ZZLS2,ZZSS,IKLIN,ZCOEFLIN)
  !
  PLSTHS(:,:,:)=VER_INTERP_LIN(PLSTHS,IKLIN,ZCOEFLIN)
  !
  IF ( SIZE(PLSRVM,1) /= 0 ) THEN
    PLSRVS(:,:,:)  = VER_INTERP_LIN(PLSRVS,IKLIN,ZCOEFLIN)
  END IF
  !
END IF
!
! compute the LS tendencies
!
PLSTHS(:,:,:)   = (PLSTHS(:,:,:) - PLSTHM(:,:,:)) / ZTIME
!
IF ( SIZE(PLSRVM,1) /= 0 ) THEN
  PLSRVS(:,:,:)   = (PLSRVS(:,:,:) - PLSRVM(:,:,:)) / ZTIME
END IF
IF ( SIZE(PLSZWSM,1) /= 0 ) THEN
  PLSZWSS(:,:)   = (PLSZWSS(:,:) - PLSZWSM(:,:)) / ZTIME
END IF
!
!------------------------------------------------------------------------------
!
!*      4.  COMPUTE LARGE SCALE FIELDS ON THE U GRID
!           ----------------------------------------
!
!
!*      4.1  Horizontal Bikhardt interpolation on the U grid
!
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,2,         &
               HLBCX,HLBCY,ZTLSUM,PLSUS(IIB:IIE,IJB:IJE,:))
!
IF ( GVERT_INTERP ) THEN
  !
  ZZLS1=MXM(ZZLS2)
  ZZLS1(1,:,:)=2.*ZZLS1(2,:,:)-ZZLS1(3,:,:)
  ZZSS=MXM(ZZSS)
  ZZSS(1,:,:)=2.*ZZSS(2,:,:)-ZZSS(3,:,:)
!
  CALL COEF_VER_INTERP_LIN(ZZLS1,ZZSS,IKLIN,ZCOEFLIN)
  !
  PLSUS(:,:,:)=VER_INTERP_LIN(PLSUS,IKLIN,ZCOEFLIN)
  !
END IF
!
! compute the LS tendencies
!
PLSUS(:,:,:)   = (PLSUS(:,:,:) - PLSUM(:,:,:)) / ZTIME
!
!
!------------------------------------------------------------------------------
!
!
!*      5.  COMPUTE LARGE SCALE FIELDS ON THE V GRID
!           ----------------------------------------
!
!
!*      5.1  Horizontal Bikhardt interpolation on the V grid
!
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,3,         &
               HLBCX,HLBCY,ZTLSVM,PLSVS(IIB:IIE,IJB:IJE,:))
IF ( GVERT_INTERP ) THEN
!
  !
  ZZLS1=MYM(ZZLS2)
  ZZLS1(:,1,:)=2.*ZZLS1(:,2,:)-ZZLS1(:,3,:)
  ZZSS=MZF(1,IKU,1,PZZ)
  ZZSS(:,:,IKU)=2.*ZZSS(:,:,IKU-1)-ZZSS(:,:,IKU-2)
  ZZSS=MYM(ZZSS)
  ZZSS(:,1,:)=2.*ZZSS(:,2,:)-ZZSS(:,3,:)
!
!
  CALL COEF_VER_INTERP_LIN(ZZLS1,ZZSS,IKLIN,ZCOEFLIN)
  !
  PLSVS(:,:,:)=VER_INTERP_LIN(PLSVS,IKLIN,ZCOEFLIN)
  !
END IF
!
! compute the LS tendencies
!
PLSVS(:,:,:)   = (PLSVS(:,:,:) - PLSVM(:,:,:)) / ZTIME
!
DEALLOCATE(ZTLSUM,ZTLSVM,ZTLSWM,ZTLSTHM,ZTLSRVM,ZTZWS)
IF(GVERT_INTERP) DEALLOCATE(ZTZS,ZZS)
IF(GVERT_INTERP) DEALLOCATE(ZTZSMT,ZZSMT)
!
NULLIFY(TZLSFIELD_ll)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSUS)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSVS)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSWS)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSTHS)
IF(SIZE(PLSRVM) /= 0) CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSRVS)
IF(SIZE(PLSZWSM) /= 0) CALL ADD2DFIELD_ll(TZLSFIELD_ll, PLSZWSS)
CALL UPDATE_HALO_ll(TZLSFIELD_ll,IINFO_ll)
CALL CLEANLIST_ll(TZLSFIELD_ll)
!------------------------------------------------------------------------------
!
CALL GOTO_MODEL(KMI)
!
END SUBROUTINE SPAWN_LS_n   
