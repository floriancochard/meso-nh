!MNH_LIC Copyright 1997-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#######################
MODULE MODI_SPAWN_PRESSURE2
!#######################
!
INTERFACE
!
      SUBROUTINE SPAWN_PRESSURE2(KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,   &
                                PZZ_LS,PZZ,PTHVT, PPABST                    )
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ_LS     ! purely interpolated alt.
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ        ! model 2 altitudes
!                                                  !   model 2
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHVT      ! virt. pot. temp. at t
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PPABST    ! model 2 pressure a t
!
END SUBROUTINE SPAWN_PRESSURE2
!
END INTERFACE
!
END MODULE MODI_SPAWN_PRESSURE2
!
!     #######################################################################
      SUBROUTINE SPAWN_PRESSURE2(KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,   &
                                PZZ_LS,PZZ,PTHVT, PPABST                    )
!     #######################################################################
!
!!****  *SPAWN_PRESSURE2 * - subroutine generating the model 2 pressure
!!                      fields, consistently with the spawning model 1,
!!                      and the model 2 thermodynamic variables.
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!      The model 2 variables are transmitted by argument (P or K prefixes),
!!    while the ones of model 1 are declared through calls to MODD_...
!!    (X or N prefixes)
!!
!!      For the case where the resolution ratio between models is 1,
!!    the horizontal interpolation becomes a simple equality.
!!      For the general case where resolution ratio is not egal to one,
!!    fields are interpolated using 2 types of interpolations:
!!                 1. Clark and Farley (JAS 1984) on 9 points
!!                 2. Bikhardt on 16 points
!!
!!       In case of resolution change:
!!
!!    1. Model 1 top pressure is computed
!!    2. Hydrostatic pressure is computed from top in model 1
!!    3. Difference between absolute pressure and hyd. pressure is kept
!!    4. Model top pressure is interpolated
!!    5. Difference between the pressures is interpolated
!!    6. Hydrostatic pressure is computed in model 2
!!    7. Absolute pressure is recovered as sum of hydrostatic pressure
!!         and difference between absolute pressure and hydrostatic pressure
!!
!!    EXTERNAL
!!    --------
!!
!!      Routine BIKHARDT      : to perform horizontal interpolations
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains parameters
!!      Module MODD_CONF       : contains NVERB
!!      Module MODD_CONF1      : contains CONF_MODEL(1)%NRR (total Number of moist variables)
!!      Module MODD_FIELD1     : contains pronostic variables of model 1
!!      Module MODD_GRID1      : contains grid variables
!!
!!    REFERENCE
!!    ---------
!!
!!       Book1 of the documentation
!!       SUBROUTINE SPAWN_PRESSURE2 (Book2 of the documentation)
!!
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    10/07/97
!!                  14/09/97 (V. Masson) use of thv as dummy argument
!!      Modification 20/05/06 Remove Clark and Farley interpolation
!!                  2014     (M.Faivre) parallelization
!!                  10/02/15 (M.Moge) correction of M.Faivre's parallelization attempt
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                  05/03/2018 (J.Escobar) bypass gridnesting special case KD(X/Y)RATIO == 1 not parallelized
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_PARAMETERS       ! Declarative modules
USE MODD_CONF
USE MODD_CST
!
USE MODD_CONF_n, ONLY: CONF_MODEL
USE MODD_LBC_n,  ONLY: LBC_MODEL
USE MODD_LUNIT_n,ONLY: LUNIT_MODEL
USE MODD_REF_n,  ONLY: REF_MODEL
!
USE MODD_BIKHARDT_n
USE MODD_VER_INTERP_LIN
USE MODD_SPAWN
!
USE MODI_SHUMAN
USE MODI_BIKHARDT
USE MODI_COMPUTE_EXNER_FROM_TOP
USE MODI_COEF_VER_INTERP_LIN
USE MODI_VER_INTERP_LIN
!
USE MODE_MODELN_HANDLER
USE MODE_ll
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ_LS     ! purely interpolated alt.
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ        ! model 2 altitudes
!                                                  !   model 2
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHVT      ! virt. pot. temp. at t
!
REAL, DIMENSION(:,:,:),   INTENT(OUT) :: PPABST    ! model 2 pressure a t
!
!*       0.2    Declarations of local variables
!
INTEGER             :: JRR       ! Loop index for moist variables
INTEGER             :: IIU       ! Upper dimension in x direction (inner model)
INTEGER             :: IJU       ! Upper dimension in y direction (inner model)
INTEGER             :: IIU1      ! Upper dimension in x direction (outer model)
INTEGER             :: IJU1      ! Upper dimension in y direction (outer model)
INTEGER             :: IKE,IKU   ! vertical indexes for models 1 and 2
INTEGER             :: JCOUNT      ! iterative loop counter
INTEGER             :: JINSTANT    ! 1: time t-dt, 2: time t
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: & ! MODEL 1 VARIABLES
 ZEXN1,       & ! Exner functions at mass points               at t or t-dt
 ZTHV1,       & ! virtual potential temperature at mass points at t or t-dt
 ZHYDEXN1,    & ! hydrostatic Exner functions at mass points   at t or t-dt
 ZSUMR          ! sum of water mixing ratios (at t-dt or t)  
REAL, DIMENSION(SIZE(XTHT1,1),SIZE(XTHT1,2)) ::              & ! MODEL 1 VARIABLES
 ZHYDEXNTOP1    ! model top Exner functions                    at t or t-dt
!$20140709
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZHYDEXNTOP1_C
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZHYDEXN1_C
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEXN1_C
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: & ! MODEL 2 VARIABLES
 ZGRIDA,      & ! mass point altitudes with purely interpoled orography
 ZGRIDB,      & ! mass point altitudes with new orography
 ZTHV2,       & ! virtual potential temperature at mass points at t or t-dt
 ZHYDEXN2,    & ! hydrostatic Exner functions at mass points   at t or t-dt
 ZEXNMHEXN2     ! Exner function minus hydrostatic Ex. f.      at t or t-dt
REAL, DIMENSION(SIZE(PTHVT,1),SIZE(PTHVT,2)) ::            & ! MODEL 2 VARIABLES
 ZHYDEXNTOP2    ! model top Exner functions                    at t or t-dt
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORK
INTEGER  :: IMI
INTEGER  :: JI, IDIMX_C,IDIMY_C
INTEGER  :: IINFO_ll
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE:
!              ---------
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL GOTO_MODEL(2)
!
IIU = SIZE(PTHVT,1)
IJU = SIZE(PTHVT,2)
IIU1= SIZE(XTHT1,1)
IJU1= SIZE(XTHT1,2)
IKU=SIZE(PZZ,3)
IKE=IKU-JPVEXT
!
!-------------------------------------------------------------------------------
!
!*       2.    NO CHANGE OF RESOLUTION:
!              -----------------------
!
!
!!$IF (KDXRATIO == 1 .AND. KDYRATIO == 1 ) THEN
!!$!
!!$  PPABST  (:,:,:)   =  FIELD_MODEL(1)%XPABST  (KXOR:KXEND,KYOR:KYEND,:)
!!$!
!!$  CALL GOTO_MODEL(IMI) 
!!$  RETURN
!!$!
!!$END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    GENERAL CASE: CHANGE OF RESOLUTION
!              ----------------------------------
!
!*       2.1      Model 1 Pi and thetav
!                 ---------------------
!
  ALLOCATE(ZEXN1(IIU1,IJU1,IKU))
  ALLOCATE(ZTHV1(IIU1,IJU1,IKU))
  ALLOCATE(ZSUMR(IIU1,IJU1,IKU))
  ZSUMR(:,:,:) = 0.
  DO JRR=1,CONF_MODEL(1)%NRR
    ZSUMR(:,:,:) = ZSUMR(:,:,:) + XRT1(:,:,:,JRR)
  END DO
    !
  ZEXN1(:,:,:)=(XPABST1(:,:,:)/XP00)**(XRD/XCPD)
  IF (CONF_MODEL(1)%LUSERV) THEN
    ZTHV1(:,:,:)=XTHT1(:,:,:)*(1.+XRV/XRD*XRT1(:,:,:,1))/(1.+ZSUMR)
  ELSE
    ZTHV1(:,:,:)=XTHT1(:,:,:)
  END IF
  DEALLOCATE(ZSUMR)
!
!*       2.2      Model 1 top Exner function (guess)
!                 --------------------------
!
!* hydrostatism is supposed verified at the highest level.
!  since hydrostatic computation at mass points from flux points is non linear,
!  an iterative process is necessary to retrieve the correct hydrostatic Exner
!  pressure at model top.
!
!* guess
  ZHYDEXNTOP1(:,:)=0.5*(2.*ZEXN1(:,:,IKE)-REF_MODEL(1)%XEXNREF(:,:,IKE)+REF_MODEL(1)%XEXNREF(:,:,IKE+1))
!
  ALLOCATE(ZHYDEXN1(IIU1,IJU1,IKU))
  ALLOCATE(ZWORK   (IIU1,IJU1,IKU))
!
!* iterative loop
  DO JCOUNT=1,10
    CALL COMPUTE_EXNER_FROM_TOP(ZTHV1(:,:,IKE-1:IKE+1),XZZ1(:,:,IKE-1:IKE+1),   &
                                ZHYDEXNTOP1(:,:),                              &
                                ZWORK(:,:,IKE-1:IKE+1),ZHYDEXN1(:,:,IKE-1:IKE+1))
    ZHYDEXNTOP1(:,:)=ZHYDEXNTOP1(:,:)+ZEXN1(:,:,IKE)-ZHYDEXN1(:,:,IKE)
    IF (ALL(ABS(ZEXN1(:,:,IKE)-ZHYDEXN1(:,:,IKE))<1.E-10)) EXIT
  END DO
!
!*       2.3      Model 1 hydrostatic pressure
!                 ----------------------------
!
  CALL COMPUTE_EXNER_FROM_TOP(ZTHV1,XZZ1,ZHYDEXNTOP1,ZWORK,ZHYDEXN1)
!
  DEALLOCATE(ZTHV1)
  DEALLOCATE(ZWORK)
  !
  CALL GOTO_MODEL(1)
  CALL GO_TOMODEL_ll(1, IINFO_ll)
  CALL GET_CHILD_DIM_ll(2, IDIMX_C, IDIMY_C, IINFO_ll)
  ! 2D
  ALLOCATE(ZHYDEXNTOP1_C(IDIMX_C,IDIMY_C))
  ZHYDEXNTOP1_C=0.
  ! 3D
  ALLOCATE(ZHYDEXN1_C(IDIMX_C,IDIMY_C,SIZE(ZHYDEXN1,3)))
  ALLOCATE(ZEXN1_C(IDIMX_C,IDIMY_C,SIZE(ZEXN1,3)))
  ZHYDEXN1_C =0.
  ZEXN1_C    =0.
  !
  CALL SET_LSFIELD_1WAY_ll(ZHYDEXNTOP1,ZHYDEXNTOP1_C,2)
  !
  CALL LS_FORCING_ll(2, IINFO_ll, .TRUE.)
  CALL GO_TOMODEL_ll(2, IINFO_ll)
  CALL GOTO_MODEL(2)
  CALL UNSET_LSFIELD_1WAY_ll()
  ! 3D
  DO JI=1,SIZE(ZEXN1,3)
    CALL GOTO_MODEL(1)
    CALL GO_TOMODEL_ll(1, IINFO_ll)
    !
    CALL SET_LSFIELD_1WAY_ll(ZHYDEXN1(:,:,JI),ZHYDEXN1_C(:,:,JI),2)
    CALL SET_LSFIELD_1WAY_ll(ZEXN1(:,:,JI),ZEXN1_C(:,:,JI),2)
    !
    CALL LS_FORCING_ll(2, IINFO_ll, .TRUE.)
    CALL GO_TOMODEL_ll(2, IINFO_ll)
    CALL GOTO_MODEL(2)
    CALL UNSET_LSFIELD_1WAY_ll()
  ENDDO
!
!if the child grid is the whole father grid, we first need to extrapolate
!the data on a "pseudo halo" before doing BIKHARDT interpolation
!  CALL EXTRAPOL_ON_PSEUDO_HALO(ZHYDEXNTOP1_C)
!  CALL EXTRAPOL_ON_PSEUDO_HALO(ZHYDEXN1_C)
!  CALL EXTRAPOL_ON_PSEUDO_HALO(ZEXN1_C)
!
!*       2.4      model top Exner function interpolation
!                 --------------------------------------
!
     CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                    XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                    2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,4,     &
                    LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,ZHYDEXNTOP1_C,ZHYDEXNTOP2)
     CALL MPPDB_CHECK2D(ZHYDEXNTOP2,"SPAWN_PRESS2:ZHYDEXNTOP2",PRECISION)
!
!
!*       2.5      interpolation of pi-hyd pi
!                 --------------------------
!
  ALLOCATE(ZEXNMHEXN2(IIU,IJU,IKU))
!
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO,KDYRATIO,1,     &
                   LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,           &
                   ZEXN1_C-ZHYDEXN1_C,ZEXNMHEXN2)
    CALL MPPDB_CHECK3D(ZEXNMHEXN2,"SPAWN_PRESS2:ZEXNMHEXN2",PRECISION)
!
  DEALLOCATE(ZEXN1)
  DEALLOCATE(ZHYDEXN1)
  DEALLOCATE(ZEXN1_C)
  DEALLOCATE(ZHYDEXN1_C)
  DEALLOCATE(ZHYDEXNTOP1_C)
!
!* vertical interpolation
!
  ALLOCATE(ZGRIDA(IIU,IJU,IKU))
  ALLOCATE(ZGRIDB(IIU,IJU,IKU))
  ZGRIDA(:,:,:)=MZF(1,IKU,1,PZZ_LS(:,:,:))
  ZGRIDA(:,:,IKU)=2.*ZGRIDA(:,:,IKU-1)-ZGRIDA(:,:,IKU-2)
  ZGRIDB(:,:,:)=MZF(1,IKU,1,PZZ(:,:,:))
  ZGRIDB(:,:,IKU)=2.*ZGRIDB(:,:,IKU-1)-ZGRIDB(:,:,IKU-2)
  CALL COEF_VER_INTERP_LIN(ZGRIDA(:,:,:),ZGRIDB(:,:,:))
!
  ZEXNMHEXN2(:,:,:)=VER_INTERP_LIN(ZEXNMHEXN2(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
!
  DEALLOCATE(ZGRIDA)
  DEALLOCATE(ZGRIDB)
  DEALLOCATE(XCOEFLIN)
  DEALLOCATE(NKLIN)
!
!*       2.6      Model 2 thetav
!                 --------------
!
  ALLOCATE(ZTHV2(IIU,IJU,IKU))
  ZTHV2(:,:,:)=PTHVT(:,:,:)
!
!*       2.7      Model 2 hydrostatic pressure
!                 ----------------------------
!
  ALLOCATE(ZHYDEXN2(IIU,IJU,IKU))
  ALLOCATE(ZWORK   (IIU,IJU,IKU))
!
  CALL COMPUTE_EXNER_FROM_TOP(ZTHV2,PZZ,ZHYDEXNTOP2,ZWORK,ZHYDEXN2)
!
  DEALLOCATE(ZTHV2)
  DEALLOCATE(ZWORK)
!
!*       2.8      Model 2 pressure
!                 ----------------
!
  PPABST(:,:,:)=XP00*(ZEXNMHEXN2(:,:,:)+ZHYDEXN2(:,:,:))**(XCPD/XRD)
!
  DEALLOCATE(ZEXNMHEXN2)
  DEALLOCATE(ZHYDEXN2)
!
!-------------------------------------------------------------------------------
!
CALL GOTO_MODEL(IMI)
!
END SUBROUTINE SPAWN_PRESSURE2
!
