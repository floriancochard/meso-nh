!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_INI_SPAWN_LS_n
!     #########################
!
INTERFACE
!
      SUBROUTINE INI_SPAWN_LS_n(KDAD,PTSTEP,KMI,                          &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,      &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,      &
                    KDXRATIO,KDYRATIO,              &
                    HLBCX,HLBCY,PZZ,PZHAT,                                &
                    OSLEVE,PLEN1,PLEN2,                                   &
                      PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,            &
                      PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS,PLSZWSS,            &
                      KKLIN_LBXU,PCOEFLIN_LBXU,KKLIN_LBYU,PCOEFLIN_LBYU,  &
                      KKLIN_LBXV,PCOEFLIN_LBXV,KKLIN_LBYV,PCOEFLIN_LBYV,  &
                      KKLIN_LBXW,PCOEFLIN_LBXW,KKLIN_LBYW,PCOEFLIN_LBYW,  &
                      KKLIN_LBXM,PCOEFLIN_LBXM,KKLIN_LBYM,PCOEFLIN_LBYM   )
!
INTEGER,          INTENT(IN)    :: KDAD     ! model number of the DAD model
REAL,             INTENT(IN)    :: PTSTEP   !  Time step
INTEGER,          INTENT(IN)    :: KMI      ! model number
!
                                    ! interpolation coefficients
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ  ! small scale vertical grid
REAL, DIMENSION(:), INTENT(IN) :: PZHAT ! height of the Gal-Chen levels
LOGICAL,            INTENT(IN) :: OSLEVE! flag for Sleve coordinate
REAL,               INTENT(IN) :: PLEN1 ! Decay scale for smooth topography
REAL,               INTENT(IN) :: PLEN2 ! Decay scale for small-scale topography deviation
!
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
!
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSUM,PLSVM,PLSWM ! Large Scale fields
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSTHM, PLSRVM    ! at t-dt
REAL, DIMENSION(:,:),   INTENT(  OUT) :: PLSZWSM    ! at t-dt
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSUS,PLSVS,PLSWS ! Large Scale source
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSTHS, PLSRVS    ! terms
REAL, DIMENSION(:,:),   INTENT(  OUT) :: PLSZWSS    ! source terms
!  coefficients for the vertical interpolation of the LB fields
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXU,KKLIN_LBYU
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXU,PCOEFLIN_LBYU
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXV,KKLIN_LBYV
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXV,PCOEFLIN_LBYV
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXW,KKLIN_LBYW
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXW,PCOEFLIN_LBYW
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXM,KKLIN_LBYM
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXM,PCOEFLIN_LBYM
!
END SUBROUTINE INI_SPAWN_LS_n
!
END INTERFACE
!
END MODULE MODI_INI_SPAWN_LS_n
!
!     ################################################################
      SUBROUTINE INI_SPAWN_LS_n(KDAD,PTSTEP,KMI,                          &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,      &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,      &
                    KDXRATIO,KDYRATIO,              &
                    HLBCX,HLBCY,PZZ,PZHAT,                                &
                    OSLEVE,PLEN1,PLEN2,                                   &
                      PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,            &
                      PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS,PLSZWSS,            &
                      KKLIN_LBXU,PCOEFLIN_LBXU,KKLIN_LBYU,PCOEFLIN_LBYU,  &
                      KKLIN_LBXV,PCOEFLIN_LBXV,KKLIN_LBYV,PCOEFLIN_LBYV,  &
                      KKLIN_LBXW,PCOEFLIN_LBXW,KKLIN_LBYW,PCOEFLIN_LBYW,  &
                      KKLIN_LBXM,PCOEFLIN_LBXM,KKLIN_LBYM,PCOEFLIN_LBYM   )
!     ################################################################
!
!!****  *INI_SPAWN_LS_n* - Compute by interpolation the Large Scale fields
!!****                      for a nested model   
!!
!!    PURPOSE
!!    -------

!!      The purpose of INI_SPAWN_LS_n is to interpolate the Large Scale fields
!!    used by the outermost model ( model 1) on the fine mesh grid of the $n 
!!    model. This interpolation is performed first on the "horizontal" model
!!    level, then a vertical interpolation is performed if the grid is refined
!!    compared to the model 1. The coefficient are stored only for the RIM area
!!    in order to save memory and will be used for the LB fields.
!
!
!!**  METHOD
!!    ------
!!      The Large scale field are taken from the DAD model of modeln $n in order
!!    to always have enough interpolation points. These fields arrived in the
!!    subroutine by argument. At the contrary, the fields of the model number
!!    $n, are passed by the module MODD_FIELD$n.
!!      The interpolation on the model level K is performed by the Bikhart
!!    interpolation and the vertical one is a simple linear interpolation used
!!    if the orography is more detailled than the original model 1 one.
!!
!!    EXTERNAL
!!    --------
!!       Subroutine BIKHARDT      : performs the horizontal interpolation
!!
!!       Function VER_INTERP_LIN  : performs the vertical interpolation
!!
!!       Subroutine COEF_VER_INTERP_LIN: computes the vertical interpolation coef.
!!
!!       Function   MZF,MXM,MYM   : averages operators
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: JPHEXT,JPVEXT,JPCPLFILEMAX
!!
!!      Module MODD_DYN: XTSTEP_MODEL1,NCPL_TIMES
!!
!!      Module MODD_CONF: L2D
!!
!!      Module MODD_LSFIELD$n : XLSTHM,XLSTHS,XLSRVM,XLSRVS,XLSUM,XLSUS,
!!                              XLSVM,XLSVS,XLSWM,XLSWS
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
!!                   22/01/98  ( J. Stein ) add the vertical interpolation
!!                   09/07/98  ( J. Stein ) bug in the storage of the interp
!!                                          coeff for U
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 18/12/2015 : Correction of bug in bound in // for NHALO <>1 
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
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
USE MODD_CONF
USE MODD_FIELD_n     ! modules relative to the outer model $n
USE MODD_LSFIELD_n
USE MODD_GRID_n   
!
USE MODI_BIKHARDT
USE MODI_SHUMAN
USE MODI_COEF_VER_INTERP_LIN
USE MODI_VER_INTERP_LIN
USE MODI_VERT_COORD
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
INTEGER,          INTENT(IN)    :: KDAD     ! model number of the DAD model
REAL,             INTENT(IN)    :: PTSTEP   !  Time step
INTEGER,          INTENT(IN)    :: KMI      ! model number
!
                                    ! interpolation coefficients
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
REAL, DIMENSION(:,:,:), INTENT(IN) :: PZZ  ! small scale vertical grid
REAL, DIMENSION(:), INTENT(IN) :: PZHAT ! height of the Gal-Chen levels
LOGICAL,            INTENT(IN) :: OSLEVE! flag for Sleve coordinate
REAL,               INTENT(IN) :: PLEN1 ! Decay scale for smooth topography
REAL,               INTENT(IN) :: PLEN2 ! Decay scale for small-scale topography deviation
!
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
!
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSUM,PLSVM,PLSWM ! Large Scale fields
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSTHM, PLSRVM    ! at t-dt
REAL, DIMENSION(:,:),   INTENT(  OUT) :: PLSZWSM     ! LS at t-dt
REAL, DIMENSION(:,:),   INTENT(  OUT) :: PLSZWSS    ! source terms
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSUS,PLSVS,PLSWS ! Large Scale source
REAL, DIMENSION(:,:,:), INTENT(  OUT) :: PLSTHS, PLSRVS    ! terms
!  coefficients for the vertical interpolation of the LB fields
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXU,KKLIN_LBYU
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXU,PCOEFLIN_LBYU
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXV,KKLIN_LBYV
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXV,PCOEFLIN_LBYV
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXW,KKLIN_LBYW
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXW,PCOEFLIN_LBYW
INTEGER, DIMENSION(:,:,:), INTENT(  OUT) :: KKLIN_LBXM,KKLIN_LBYM
REAL,    DIMENSION(:,:,:), INTENT(  OUT) :: PCOEFLIN_LBXM,PCOEFLIN_LBYM
!
!*       0.2   declarations of local variables
!
!
REAL                   :: ZTIME                   ! Interpolation length
INTEGER                :: JI                      ! Loop index
INTEGER, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: IKLIN
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZCOEFLIN
INTEGER                :: ILBX,ILBY,ILBX2,ILBY2
INTEGER                :: IIB,IJB,IIE,IJE,IIU,IJU,IKU
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),1)                       :: ZZSLS
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),1)                       :: ZZSMTLS
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZZSS
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZZLS1
                                    ! Gal Chen grid
                                    ! with the smooth LS orography on the SS grid
                                    ! at the W location then at the U location
REAL   , DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PCOEFLIN_LBXM,3))   :: ZZLS2
                                    ! Gal Chen grid
                                    ! with the smooth LS orography on the SS grid
                                    ! at the mass location then at the V location
LOGICAL                :: GVERT_INTERP
!
TYPE(LIST_ll), POINTER :: TZLSFIELD_ll   ! list of LS fields
! Variables used for LS communications
INTEGER :: IINFO_ll, IDIMX, IDIMY
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTLSUM, ZTLSVM, ZTLSWM, ZTLSTHM, ZTLSRVM
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTLSUS, ZTLSVS, ZTLSWS, ZTLSTHS, ZTLSRVS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTZS,ZZS 
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTZSMT,ZZSMT 
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTZWS,ZTZWSS
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
!
ALLOCATE(ZTLSUM(IDIMX,IDIMY,SIZE(PLSUM,3)))
ALLOCATE(ZTLSVM(IDIMX,IDIMY,SIZE(PLSVM,3)))
ALLOCATE(ZTLSWM(IDIMX,IDIMY,SIZE(PLSWM,3)))
ALLOCATE(ZTLSTHM(IDIMX,IDIMY,SIZE(PLSTHM,3)))
IF(SIZE(PLSZWSM) /= 0) ALLOCATE(ZTZWS(IDIMX,IDIMY))
IF(SIZE(PLSRVM) /= 0) ALLOCATE(ZTLSRVM(IDIMX,IDIMY,SIZE(PLSRVM,3)))
!
IF(GVERT_INTERP) THEN
  ALLOCATE(ZTZS(IDIMX,IDIMY,1))
  ALLOCATE(ZTZSMT(IDIMX,IDIMY,1))
  ALLOCATE(ZZS(SIZE(XZS,1),SIZE(XZS,2),1))
  ZZS(:,:,1)=XZS(:,:)
  ALLOCATE(ZZSMT(SIZE(XZS,1),SIZE(XZS,2),1))
  ZZSMT(:,:,1)=XZSMT(:,:)
END IF
!
IF ( SIZE(PLSTHS,1) /= 0 ) THEN
  ALLOCATE(ZTLSUS(IDIMX,IDIMY,SIZE(PLSUS,3)))
  ALLOCATE(ZTLSVS(IDIMX,IDIMY,SIZE(PLSVS,3)))
  ALLOCATE(ZTLSWS(IDIMX,IDIMY,SIZE(PLSWS,3)))
  ALLOCATE(ZTLSTHS(IDIMX,IDIMY,SIZE(PLSTHS,3)))
ENDIF
IF ( SIZE(PLSZWSS) /= 0 ) ALLOCATE(ZTZWSS(IDIMX,IDIMY))
IF ( SIZE(PLSRVS) /= 0 ) ALLOCATE(ZTLSRVS(IDIMX,IDIMY,SIZE(PLSRVS,3)))
!
!         1.3  Specify the ls "source" fields and receiver fields
!
CALL SET_LSFIELD_1WAY_ll(XLSUM, ZTLSUM, KMI)
IF ( SIZE(PLSZWSM,1) /= 0 ) &
  CALL SET_LSFIELD_1WAY_ll(XLSZWSM, ZTZWS, KMI)
CALL SET_LSFIELD_1WAY_ll(XLSVM, ZTLSVM, KMI)
CALL SET_LSFIELD_1WAY_ll(XLSWM, ZTLSWM, KMI)
CALL SET_LSFIELD_1WAY_ll(XLSTHM, ZTLSTHM, KMI)
IF ( SIZE(PLSRVM,1) /= 0 ) &
  CALL SET_LSFIELD_1WAY_ll(XLSRVM, ZTLSRVM, KMI)
!
IF ( SIZE(PLSTHS,1) /= 0 ) THEN
  CALL SET_LSFIELD_1WAY_ll(XLSUS, ZTLSUS, KMI)
  CALL SET_LSFIELD_1WAY_ll(XLSVS, ZTLSVS, KMI)
  CALL SET_LSFIELD_1WAY_ll(XLSWS, ZTLSWS, KMI)
  CALL SET_LSFIELD_1WAY_ll(XLSTHS, ZTLSTHS, KMI)
  IF ( SIZE(PLSRVM,1) /= 0 ) &
    CALL SET_LSFIELD_1WAY_ll(XLSRVS, ZTLSRVS, KMI)
  IF ( SIZE(PLSZWSM,1) /= 0 ) &
    CALL SET_LSFIELD_1WAY_ll(XLSZWSS, ZTZWSS, KMI)
ENDIF
!
IF ( GVERT_INTERP ) THEN
  CALL SET_LSFIELD_1WAY_ll(ZZS, ZTZS, KMI)
  CALL SET_LSFIELD_1WAY_ll(ZZSMT, ZTZSMT, KMI)
END IF
!
!        1.4  Communication
!
CALL LS_FORCING_ll(KMI, IINFO_ll)
!
!        1.5  Back to the (current) child model
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
CALL UNSET_LSFIELD_1WAY_ll()
!
!*      2.  COMPUTE THE GRIDS AT THE W LOCALIZATION AND THE COUPLING INSTANTS
!           -----------------------------------------------------------------
!
IF ( GVERT_INTERP ) THEN
  !
  ! the orography of the DAD model is interpolated on the fine mesh model
  !
  ZZSLS=0.
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTZS,ZZSLS(IIB:IIE,IJB:IJE,:)        )
  CALL MPPDB_CHECK3D(ZZSLS,"INI_SPAWN_LS::ZZSLS",PRECISION)
  !
  ZZSMTLS=0.
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTZS,ZZSMTLS(IIB:IIE,IJB:IJE,:)      )
  !
  !
  CALL VERT_COORD(OSLEVE,ZZSLS(:,:,1),ZZSMTLS(:,:,1),PLEN1,PLEN2,PZHAT,ZZLS1)
  !
END IF
!
IF ( SIZE(PLSTHS,1) /= 0 ) THEN
  DO JI=1,JPCPLFILEMAX
    IF ( NCPL_TIMES(JI,1) /= NUNDEF ) THEN
      NCPL_TIMES(JI,KMI) =  NINT( ((NCPL_TIMES(JI,1)-2)*XTSTEP_MODEL1)   &
                 / PTSTEP ) + 2
    END IF
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*      2.  COMPUTE LARGE SCALE FIELDS ON THE W GRID
!           ----------------------------------------
!
!
!*      2.1  Horizontal Bikhardt interpolation on the W grid
!
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,4,         &
               HLBCX,HLBCY,ZTLSWM,PLSWM(IIB:IIE,IJB:IJE,:))
!
IF ( SIZE(PLSWS,1) /= 0 ) THEN
  ! temporal distance between the two first coupling instants
  ZTIME = (NCPL_TIMES(NCPL_CUR,1)-2) * XTSTEP_MODEL1
  !
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,4,         &
                 HLBCX,HLBCY,ZTLSWM+ZTIME*ZTLSWS,PLSWS(IIB:IIE,IJB:IJE,:))
  !
END IF
!
!
!*      2.2  Vertical linear interpolation on the W grid
!
!
IF ( GVERT_INTERP ) THEN
!
  CALL COEF_VER_INTERP_LIN(ZZLS1,PZZ,IKLIN,ZCOEFLIN)
  !
  PLSWM(:,:,:)=VER_INTERP_LIN(PLSWM,IKLIN,ZCOEFLIN)
  !
  IF ( SIZE(PLSTHS,1) /= 0 ) THEN
    PLSWS(:,:,:)=VER_INTERP_LIN(PLSWS,IKLIN,ZCOEFLIN)
  END IF
!
!*       2.3 Store the interpolation coefficients in the RIM area
!
  ILBX2=SIZE(KKLIN_LBXW,1)
  IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
    ILBX=ILBX2/2
  ELSE
    ILBX=ILBX2
  ENDIF
!
  IF(LWEST_ll( )) THEN
    KKLIN_LBXW(1:ILBX,:,:)   =   IKLIN(IIB:IIB-1+ILBX,:,:)
    PCOEFLIN_LBXW(1:ILBX,:,:)=ZCOEFLIN(IIB:IIB-1+ILBX,:,:)
  ENDIF
!
  IF(LEAST_ll( )) THEN
    KKLIN_LBXW(ILBX2-ILBX+1:ILBX2,:,:)   =    IKLIN(IIE+1-ILBX:IIE,:,:)
    PCOEFLIN_LBXW(ILBX2-ILBX+1:ILBX2,:,:)=ZCOEFLIN (IIE+1-ILBX:IIE,:,:)
  ENDIF
!
  IF (.NOT. L2D ) THEN
    ILBY2=SIZE(KKLIN_LBYW,2)
    IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
      ILBY=ILBY2/2
    ELSE
      ILBY=ILBY2
    ENDIF
!
    IF(LSOUTH_ll( )) THEN
      KKLIN_LBYW(:,1:ILBY,:)   =   IKLIN(:,IJB:IJB-1+ILBY,:)
      PCOEFLIN_LBYW(:,1:ILBY,:)=ZCOEFLIN(:,IJB:IJB-1+ILBY,:)
    ENDIF
!
    IF(LNORTH_ll( )) THEN
      KKLIN_LBYW(:,ILBY2-ILBY+1:ILBY2,:)   =   IKLIN(:,IJE+1-ILBY:IJE,:)
      PCOEFLIN_LBYW(:,ILBY2-ILBY+1:ILBY2,:)=ZCOEFLIN(:,IJE+1-ILBY:IJE,:)
    ENDIF
  END IF
!
END IF
!
!
!*      2.4  LS tendencies
!
!
IF ( SIZE(PLSVS,1) /= 0 ) THEN
  PLSWS(:,:,:)   = (PLSWS(:,:,:) - PLSWM(:,:,:)) / ZTIME
END IF
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
               HLBCX,HLBCY,ZTLSTHM,PLSTHM(IIB:IIE,IJB:IJE,:))
!
IF ( SIZE(PLSRVM,1) /= 0 ) THEN
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTLSRVM,PLSRVM(IIB:IIE,IJB:IJE,:))
END IF
!
IF ( SIZE(PLSZWSM,1) /= 0 ) THEN
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTZWS,PLSZWSM(IIB:IIE,IJB:IJE))
END IF

IF ( SIZE(PLSTHS,1) /= 0 ) THEN
  !
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                 HLBCX,HLBCY,ZTLSTHM+ZTIME*ZTLSTHS,PLSTHS(IIB:IIE,IJB:IJE,:))
  !
  IF ( SIZE(PLSRVM,1) /= 0 ) THEN
    CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                   PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                   2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                   HLBCX,HLBCY,ZTLSRVM+ZTIME*ZTLSRVS,PLSRVS(IIB:IIE,IJB:IJE,:))
  !
  END IF
  !
  IF ( SIZE(PLSZWSM,1) /= 0 ) THEN
    CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                   PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                   2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
                   HLBCX,HLBCY,ZTZWS+ZTIME*ZTZWSS,PLSZWSS(IIB:IIE,IJB:IJE))
  !
  END IF
END IF
!
!*      3.2  Vertical linear interpolation on the mass grid
!
!
IF ( GVERT_INTERP ) THEN
  !
  IKU = SIZE(PZZ,3)
  !
  ZZLS2=MZF(1,IKU,1,ZZLS1)
  ZZLS2(:,:,IKU)=2.*ZZLS2(:,:,IKU-1)-ZZLS2(:,:,IKU-2)
  ZZSS=MZF(1,IKU,1,PZZ)
  ZZSS(:,:,IKU)=2.*ZZSS(:,:,IKU-1)-ZZSS(:,:,IKU-2)
  !
  CALL COEF_VER_INTERP_LIN(ZZLS2,ZZSS,IKLIN,ZCOEFLIN)
  !
  PLSTHM(:,:,:)=VER_INTERP_LIN(PLSTHM,IKLIN,ZCOEFLIN)
  !
  IF ( SIZE(PLSRVM,1) /= 0 ) THEN
    PLSRVM(:,:,:)  = VER_INTERP_LIN(PLSRVM,IKLIN,ZCOEFLIN)
  END IF
  !
  IF ( SIZE(PLSTHS,1) /= 0 ) THEN
    PLSTHS(:,:,:)=VER_INTERP_LIN(PLSTHS,IKLIN,ZCOEFLIN)
  END IF
  !
  IF ( SIZE(PLSRVS,1) /= 0  ) THEN
    PLSRVS(:,:,:)  = VER_INTERP_LIN(PLSRVS,IKLIN,ZCOEFLIN)
  END IF
!
!*       3.3 Store the interpolation coef. in the RIM area on the mass grid
!
  ILBX2=SIZE(KKLIN_LBXM,1)
  IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
    ILBX=ILBX2/2
  ELSE
    ILBX=ILBX2
  ENDIF
!
  IF(LWEST_ll( )) THEN
    KKLIN_LBXM(1:ILBX,:,:)   =   IKLIN(IIB:IIB-1+ILBX,:,:)
    PCOEFLIN_LBXM(1:ILBX,:,:)=ZCOEFLIN(IIB:IIB-1+ILBX,:,:)
  ENDIF
!
  IF(LEAST_ll( )) THEN
    KKLIN_LBXM(ILBX2-ILBX+1:ILBX2,:,:)   =    IKLIN(IIE+1-ILBX:IIE,:,:)
    PCOEFLIN_LBXM(ILBX2-ILBX+1:ILBX2,:,:)=ZCOEFLIN (IIE+1-ILBX:IIE,:,:)
  ENDIF
!
  IF (.NOT. L2D ) THEN
    ILBY2=SIZE(KKLIN_LBYM,2)
    IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
      ILBY=ILBY2/2
    ELSE
      ILBY=ILBY2
    ENDIF
!
    IF(LSOUTH_ll( )) THEN
      KKLIN_LBYM(:,1:ILBY,:)   =   IKLIN(:,IJB:IJB-1+ILBY,:)
      PCOEFLIN_LBYM(:,1:ILBY,:)=ZCOEFLIN(:,IJB:IJB-1+ILBY,:)
    ENDIF
!
    IF(LNORTH_ll( )) THEN
      KKLIN_LBYM(:,ILBY2-ILBY+1:ILBY2,:)   =   IKLIN(:,IJE+1-ILBY:IJE,:)
      PCOEFLIN_LBYM(:,ILBY2-ILBY+1:ILBY2,:)=ZCOEFLIN(:,IJE+1-ILBY:IJE,:)
    ENDIF
  END IF
  !
END IF
!
!*      3.4  LS tendencies on the mass grid
!
!
IF ( SIZE(PLSTHS,1) /= 0 ) THEN
  !
  PLSTHS(:,:,:)  = (PLSTHS(:,:,:) - PLSTHM(:,:,:)) / ZTIME
  !
  IF ( SIZE(PLSRVM,1) /= 0 ) THEN
    PLSRVS(:,:,:)  = (PLSRVS(:,:,:) - PLSRVM(:,:,:)) / ZTIME
  END IF
  !
  IF ( SIZE(PLSZWSM,1) /= 0 ) THEN
    PLSZWSS(:,:)  = (PLSZWSS(:,:) - PLSZWSM(:,:)) / ZTIME
  END IF
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
               HLBCX,HLBCY,ZTLSUM,PLSUM(IIB:IIE,IJB:IJE,:))
CALL MPPDB_CHECK3D(PLSUM,"INI_SPAWN_LS::PLSUM",PRECISION)

!
IF ( SIZE(PLSUS,1) /= 0 ) THEN
  !
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,2,         &
                 HLBCX,HLBCY,ZTLSUM+ZTIME*ZTLSUS,PLSUS(IIB:IIE,IJB:IJE,:))
!
END IF
!
!*      4.2  Vertical linear interpolation on the U grid
!
!
IF ( GVERT_INTERP ) THEN
!
  ZZLS1=MXM(ZZLS2)
  ZZLS1(1,:,:)=2.*ZZLS1(2,:,:)-ZZLS1(3,:,:)
  ZZSS=MXM(ZZSS)
  ZZSS(1,:,:)=2.*ZZSS(2,:,:)-ZZSS(3,:,:)
!
!
  CALL COEF_VER_INTERP_LIN(ZZLS1,ZZSS,IKLIN,ZCOEFLIN)
  !
  PLSUM(:,:,:)=VER_INTERP_LIN(PLSUM,IKLIN,ZCOEFLIN)
  !
  IF ( SIZE(PLSUS,1) /= 0 ) THEN
    PLSUS(:,:,:)=VER_INTERP_LIN(PLSUS,IKLIN,ZCOEFLIN)
  END IF
  !
!
!
!*       4.3 Store the interpolation coefficients in the RIM area on the U grid
!
  ILBX2=SIZE(KKLIN_LBXU,1)
  IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
    ILBX=ILBX2/2
  ELSE
    ILBX=ILBX2
  ENDIF
!
!
  IF(LWEST_ll( )) THEN
    KKLIN_LBXU(1:ILBX,:,:)   =   IKLIN(IIB+1:IIB+ILBX,:,:)  !  C grid
    PCOEFLIN_LBXU(1:ILBX,:,:)=ZCOEFLIN(IIB+1:IIB+ILBX,:,:)
  ENDIF
!
  IF(LEAST_ll( )) THEN
    KKLIN_LBXU(ILBX2-ILBX+1:ILBX2,:,:)   =    IKLIN(IIE+1-ILBX:IIE,:,:)
    PCOEFLIN_LBXU(ILBX2-ILBX+1:ILBX2,:,:)=ZCOEFLIN (IIE+1-ILBX:IIE,:,:)
  ENDIF
!
  IF (.NOT. L2D ) THEN
    ILBY2=SIZE(KKLIN_LBYU,2)
    IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
      ILBY=ILBY2/2
    ELSE
      ILBY=ILBY2
    ENDIF
!
    IF(LSOUTH_ll( )) THEN
      KKLIN_LBYU(:,1:ILBY,:)   =   IKLIN(:,IJB:IJB-1+ILBY,:)
      PCOEFLIN_LBYU(:,1:ILBY,:)=ZCOEFLIN(:,IJB:IJB-1+ILBY,:)
    ENDIF
!
    IF(LNORTH_ll( )) THEN
      KKLIN_LBYU(:,ILBY2-ILBY+1:ILBY2,:)   =   IKLIN(:,IJE+1-ILBY:IJE,:)
      PCOEFLIN_LBYU(:,ILBY2-ILBY+1:ILBY2,:)=ZCOEFLIN(:,IJE+1-ILBY:IJE,:)
    ENDIF
  END IF
!
END IF
!
!*      4.4  LS tendencies on the U grid
!
!
IF ( SIZE(PLSUS,1) /= 0 ) THEN
  !
  PLSUS(:,:,:)  = (PLSUS(:,:,:) - PLSUM(:,:,:)) / ZTIME
  !
END IF
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
               HLBCX,HLBCY,ZTLSVM,PLSVM(IIB:IIE,IJB:IJE,:))
!
IF ( SIZE(PLSTHS,1) /= 0 ) THEN
  !
  CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
                 PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
                 2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,3,         &
                 HLBCX,HLBCY,ZTLSVM+ZTIME*ZTLSVS,PLSVS(IIB:IIE,IJB:IJE,:))
!
END IF
!
!
!*      5.2  Vertical linear interpolation on the V grid
!
!
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
  PLSVM(:,:,:)=VER_INTERP_LIN(PLSVM,IKLIN,ZCOEFLIN)
  !
  IF ( SIZE(PLSVS,1) /= 0 ) THEN
    PLSVS(:,:,:)=VER_INTERP_LIN(PLSVS,IKLIN,ZCOEFLIN)
  END IF
!
!
!*       5.3 Store the interpolation coefficients in the RIM area on the V grid
!
  ILBX2=SIZE(KKLIN_LBXV,1)
  IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
    ILBX=ILBX2/2
  ELSE
    ILBX=ILBX2
  ENDIF
!
  IF(LWEST_ll( )) THEN
    KKLIN_LBXV(1:ILBX,:,:)   =   IKLIN(IIB:IIB-1+ILBX,:,:)
    PCOEFLIN_LBXV(1:ILBX,:,:)=ZCOEFLIN(IIB:IIB-1+ILBX,:,:)
  ENDIF
!
  IF(LEAST_ll( )) THEN
    KKLIN_LBXV(ILBX2-ILBX+1:ILBX2,:,:)   =    IKLIN(IIE+1-ILBX:IIE,:,:)
    PCOEFLIN_LBXV(ILBX2-ILBX+1:ILBX2,:,:)=ZCOEFLIN (IIE+1-ILBX:IIE,:,:)
  ENDIF
!
  IF (.NOT. L2D ) THEN
    ILBY2=SIZE(KKLIN_LBYV,2)
    IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
      ILBY=ILBY2/2
    ELSE
      ILBY=ILBY2
    ENDIF
!
    IF(LSOUTH_ll( )) THEN
      KKLIN_LBYV(:,1:ILBY,:)   =   IKLIN(:,IJB+1:IJB+ILBY,:)  !  C grid
      PCOEFLIN_LBYV(:,1:ILBY,:)=ZCOEFLIN(:,IJB+1:IJB+ILBY,:)
    ENDIF
!
    IF(LNORTH_ll( )) THEN
      KKLIN_LBYV(:,ILBY2-ILBY+1:ILBY2,:)   =   IKLIN(:,IJE+1-ILBY:IJE,:)
      PCOEFLIN_LBYV(:,ILBY2-ILBY+1:ILBY2,:)=ZCOEFLIN(:,IJE+1-ILBY:IJE,:)
    ENDIF
  END IF
!
END IF
!
!
!*      5.4  LS tendencies on the V grid
!
!
IF ( SIZE(PLSVS,1) /= 0 ) THEN
  !
  PLSVS(:,:,:)  = (PLSVS(:,:,:) - PLSVM(:,:,:)) / ZTIME
  !
END IF
!
!
DEALLOCATE(ZTLSUM,ZTLSVM,ZTLSWM,ZTLSTHM)
IF(SIZE(PLSRVM) /= 0) DEALLOCATE(ZTLSRVM)
IF(SIZE(PLSZWSM) /= 0) DEALLOCATE(ZTZWS)
!
IF(GVERT_INTERP) DEALLOCATE(ZTZS,ZZS)
IF(GVERT_INTERP) DEALLOCATE(ZTZSMT,ZZSMT)
!
IF ( SIZE(PLSTHS,1) /= 0 ) DEALLOCATE(ZTLSUS,ZTLSVS,ZTLSWS,ZTLSTHS)
IF ( SIZE(PLSRVS,1) /= 0 ) DEALLOCATE(ZTLSRVS)
IF ( SIZE(PLSZWSS,1) /= 0 ) DEALLOCATE(ZTZWSS)
!------------------------------------------------------------------------------
NULLIFY(TZLSFIELD_ll)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSUM)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSVM)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSWM)
CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSTHM)
IF(SIZE(PLSRVM) /= 0) CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSRVM)
IF(SIZE(PLSZWSM) /= 0) CALL ADD2DFIELD_ll(TZLSFIELD_ll, PLSZWSM)
IF ( SIZE(PLSTHS,1) /= 0 ) THEN
  CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSUS)
  CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSVS)
  CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSWS)
  CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSTHS)
ENDIF
IF ( SIZE(PLSRVS,1) /= 0 ) CALL ADD3DFIELD_ll(TZLSFIELD_ll, PLSRVS)
IF ( SIZE(PLSZWSS,1) /= 0 ) CALL ADD2DFIELD_ll(TZLSFIELD_ll, PLSZWSS)
CALL UPDATE_HALO_ll(TZLSFIELD_ll,IINFO_ll)
CALL CLEANLIST_ll(TZLSFIELD_ll)
!
CALL GOTO_MODEL(KMI)
!
END SUBROUTINE INI_SPAWN_LS_n 
