!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ############################
      MODULE MODI_FREE_ATM_PROFILE
!     ############################
INTERFACE
      SUBROUTINE FREE_ATM_PROFILE(TPFILE,PVAR_MX,PZMASS_MX,PZS_LS,PZSMT_LS,PCLIMGR,&
                           PF_FREE,PZ_FREE)
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
TYPE(TFILEDATA),          INTENT(IN)  :: TPFILE    ! File characteristics
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR_MX   ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! mass points altitude
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZSMT_LS  ! large scale smooth orography
REAL,                     INTENT(IN)  :: PCLIMGR   ! climatological gradient
!                                                  ! near the ground
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PF_FREE   ! mean profile of the
!                                                  ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PZ_FREE   ! discretization in x,y,z
!                                                  ! of the profile on the
!                                                  ! flat grid where zs is the
!                                                  ! minimum of both orographies
END SUBROUTINE FREE_ATM_PROFILE
END INTERFACE
END MODULE MODI_FREE_ATM_PROFILE
!     ##############################################################
      SUBROUTINE FREE_ATM_PROFILE(TPFILE,PVAR_MX,PZMASS_MX,PZS_LS,PZSMT_LS,PCLIMGR,&
                            PF_FREE,PZ_FREE)
!     ##############################################################
!
!!****  *FREE_ATM_PROFILE* - Computation of the profile of the free atmospheres
!!                           i.e. without the Boundary layer structures
!!
!!    PURPOSE
!!    -------
!!    This routine computes the profile used for the shift of a variable
!!    and the altitude of the discretization points of this profile.
!
!!    CAUTION:
!!    The shift profile is only defined on the inner vertical points of the grid.
!!
!!**  METHOD
!!    ------
!!    The profile is discretized on the vertical GS grid defined by
!!    the MESO-NH level array XZHAT and by a constant orography,
!!    corresponding to the minimum of the Arpege and MESO-NH orographies.
!!    If necessary, the profile is extrapolated under the minimum
!!    altitude of the Arpege orography with a climatological vertical
!!    gradient PCLIMGR (uniform on the whole domain).
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_GRID1     : contains grid variables for model1
!!         XZS   : orography of MESO-NH
!!         XZHAT : GS levels
!!      Module MODD_PARAMETERS
!!         JPVEXT
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/08/97
!!      C.Lac  04/2016  Modification of the free atm gradient when the top of
!!                      the model is too low
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_GRID_n
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_VER_INTERP_LIN
!
USE MODI_COEF_VER_INTERP_LIN
USE MODI_PGDFILTER
USE MODI_VER_INTERP_LIN
USE MODI_VERT_COORD
!
USE MODE_FIELD, ONLY: TFIELDDATA, TYPEINT, TYPEREAL
USE MODE_FMWRIT
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
TYPE(TFILEDATA),          INTENT(IN)  :: TPFILE    ! File characteristics
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR_MX   ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! mass points altitude
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZSMT_LS  ! large scale smooth orography
REAL,                     INTENT(IN)  :: PCLIMGR   ! climatological gradient
!                                                  ! near the ground
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PF_FREE   ! mean profile of the
!                                                  ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PZ_FREE   ! discretization in x,y,z
!                                                  ! of the profile on the
!                                                  ! flat grid where zs is the
!                                                  ! minimum of both orographies
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
REAL, PARAMETER :: XT_CLIM_GRAD = -0.0065! climatological temperature gradient
INTEGER :: ILUOUT0                       ! output listing logical unit
INTEGER :: IIU, IJU, IKB, IKE, IKU       ! array dimensions and phys. bound.
INTEGER :: JI, JJ, JK                    ! loop counters

INTEGER :: IK350
INTEGER :: IK2000,IK3000,IK4000,IK5000   ! levels just under 2000m, 3000m,
!                                        ! 4000m and 5000m above ground.
REAL    :: ZFREEGR                       ! guess of free atmosphere gradient
REAL    :: ZMIN, ZMAX                    ! lower and upper values of the
                                         ! gradients to verify the test
REAL, DIMENSION(SIZE(PZMASS_MX,3)) :: ZGR_MX      ! gradients along a vertical
LOGICAL, DIMENSION(SIZE(PZMASS_MX,3)) :: GTEST, & ! tests on the gradients
                                         GFREE    ! .T. : belongs to free atm.
INTEGER :: IK_BLTOP                               ! level just above the BL
REAL, DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2),SIZE(PZMASS_MX,3)) &
                                      :: ZF_FREE_MX ! profile of free atmosphere
                                                    ! on mixed grid
REAL, DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2)) &
                                      :: ZFREE_GR   ! gradient of low free atm.
INTEGER, DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2)) &
                                      :: IK_BL_TOP  ! level just below the top
                                                    ! of boundary layer
INTEGER, DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2)) &
                                      :: IWK_BL_TOP ! work array
REAL, DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2)) &
                                      :: ZK_BL_TOP  ! as K_BL_TOP but real
INTEGER                               :: IIMIN, IIMAX, IJMIN, IJMAX
REAL, DIMENSION(SIZE(PZMASS_MX,1),SIZE(PZMASS_MX,2)) &
                                      :: Z2D ! field to be recorded
REAL, DIMENSION(SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3)) &
                                      :: Z3D ! field to be recorded
REAL, DIMENSION(SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3)) &
                                      :: ZZMASS ! MESO-NH output mass grid
TYPE(TFIELDDATA)  :: TZFIELD
!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
IIU=SIZE(PZMASS_MX,1)
IJU=SIZE(PZMASS_MX,2)
IKB=JPVEXT+1
IKE=SIZE(PZMASS_MX,3)-JPVEXT
IKU=SIZE(PZMASS_MX,3)
!
!*       1.   Computation of the altitude of the grid for the shift profile
!             ----------------------------------------------------------------
!
CALL VERT_COORD(LSLEVE,PZS_LS,PZSMT_LS,XLEN1,XLEN2,XZHAT,PZ_FREE)
!
!-------------------------------------------------------------------------------
!
!* Computations in the following loop are performed on the mixed grid
!
DO JI=1,IIU
  DO JJ=1,IJU
!
!-------------------------------------------------------------------------------
!
!*       2.   Gradient at all levels
!             ----------------------
!
    ZGR_MX(1:IKE)=( PVAR_MX  (JI,JJ,2:IKE+1) &
                   -PVAR_MX  (JI,JJ,1:IKE  ))&
                 /( PZMASS_MX(JI,JJ,2:IKE+1) &
                   -PZMASS_MX(JI,JJ,1:IKE  ))
!
    ZGR_MX(IKE+1:IKU) = ZGR_MX(IKE)
!
!-------------------------------------------------------------------------------
!
!*       3.   Gradient 5000m above ground
!             ---------------------------
!
!*       3.1  index of level just under 2000m, 3000m and 5000m
!             ------------------------------------------------
!
!* the limits are set in case of high orography
!
      IK5000=MAX(IKB+2,MIN(IKE-2,                                            &
             COUNT(PZMASS_MX(JI,JJ,:)<MIN(5000.+PZS_LS(JI,JJ),               &
                                          MAX(6000.,2000.+PZS_LS(JI,JJ))) )  ))
      IK4000=MAX(IKB+2,MIN(IKE-2,                                            &
             COUNT(PZMASS_MX(JI,JJ,:)<MIN(4000.+PZS_LS(JI,JJ),               &
                                          MAX(5000.,1500.+PZS_LS(JI,JJ))) )  ))
      IK3000=MAX(IKB+2,MIN(IKE-2,                                            &
             COUNT(PZMASS_MX(JI,JJ,:)<MIN(3000.+PZS_LS(JI,JJ),               &
                                          MAX(4000.,1000.+PZS_LS(JI,JJ))) )  ))
      IK2000=MAX(IKB+2,MIN(IKE-2,                                            &
             COUNT(PZMASS_MX(JI,JJ,:)<    2000.+PZS_LS(JI,JJ)             )  ))
      IK2000=MIN(IK2000,IK3000-2)
      IK350 =MAX(IKB+2,MIN(IKE-2,                                            &
             COUNT(PZMASS_MX(JI,JJ,:)<    350. +PZS_LS(JI,JJ)             )  ))

!
!*       3.2  guess of free atm. gradient
!             ---------------------------
!
  IF (IK350/=IK5000) THEN
      ZFREEGR=( PVAR_MX  (JI,JJ,IK5000) &
               -PVAR_MX  (JI,JJ,IK350 ))&
             /( PZMASS_MX(JI,JJ,IK5000) &
               -PZMASS_MX(JI,JJ,IK350 ))
  ELSE
      ZFREEGR=PCLIMGR
  END IF
!
!-------------------------------------------------------------------------------
!
!*       4.   Tests on the gradients
!             ----------------------
!
!* We test the gradients separately.
!  The use of the climatological gradient is used for guesses in
!  very low stable profiles
!
      ZMIN=MAX(0.25*ZFREEGR,0.25*PCLIMGR)
      ZMAX=MAX(1.75*ZFREEGR,1.75*PCLIMGR)
!
      GTEST(:)=( ZMIN<ZGR_MX(:) .AND. ZGR_MX(:)<ZMAX )
!
!-------------------------------------------------------------------------------
!
!*       5.   Determination of the status of the point
!             ----------------------------------------
!
!* We use the following criteria:
!
!* if 3 or less following points do not verify the test, and if they are
!  surrounded by series of 3 points verifiing it, they are supposed to be
!  included in the free atmosphere (small structure in altitude).
!
!* on the contrary, series of 3 points who do not verify the test are in
!  the boundary layer
!
!* exception: above 2000m from ground, three in four following points verifiing
!  the test are always supposed in the free atmosphere, whatever the number of
!  points not verifiing it above them.
!
!* The points 5000m above ground are supposed to verify the test and
!  to be in the free atm.
!
      GTEST(IK5000:)=.TRUE.
      GFREE(IK5000:)=.TRUE.
      GFREE(:IK5000-1)=.FALSE.
!
!
!
!* Exploration begins 5000m above ground.
!
      DO JK=IK5000-1,IK350,-1
!
!
!* Point is clearly in free atmosphere
!
        IF (GTEST(JK) .AND. GFREE(JK+1)) THEN
          GFREE(JK)=.TRUE.
          CYCLE
        END IF
!
!* 2 following points verify the test
!  Does one find 1 other points in the two above verifying the test?
!
        IF (GTEST(JK) .AND. GTEST(JK+1)) THEN
          IF (GTEST(JK+2) .OR. GTEST(JK+3)) THEN
            GFREE(JK:)=.TRUE.
            CYCLE
          END IF
        END IF
!
!* Point does not verify the test
!  Does one find 2 other points just above it not verifying the test?
!  If yes and under 3000m from ground ---> Beginning of boundary layer.
!  End of the searching.
!
        IF (.NOT. GTEST(JK) .AND. .NOT. GTEST(JK+1) .AND. .NOT. GTEST(JK+2) &
                            .AND. JK <= IK2000    ) EXIT
!
!* Other cases: treated in further loop iteration
!
      END DO
!
!-------------------------------------------------------------------------------
!
!*       6.   Top of boundary layer
!             ---------------------
!
!* one level is added to remove the beginings of the boundary layer
!
      IK_BLTOP = COUNT(.NOT. GFREE) + 2
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!*       7.   End of loop on the guesses of free gradient
!             -------------------------------------------
!
      IK_BL_TOP(JI,JJ) = IK_BLTOP
!
!
  END DO
END DO
!-------------------------------------------------------------------------------
!
!*       8.   top of the boundary layer
!             -------------------------
!
!*       8.1  remove spurious values of boundary layer top level
!             --------------------------------------------------
!
!* small areas of high BL, or limits between two very different areas are
!  modified
!
IWK_BL_TOP(:,:)=IK_BL_TOP(:,:)
ZK_BL_TOP(:,:)=FLOAT(IK_BL_TOP(:,:))
CALL MPPDB_CHECK2D(ZK_BL_TOP,"FREE_ATM_PROFILE:8.1:ZK_BL_TOP",PRECISION)
!
!!$DO JI=1,IIU
!!$  DO JJ=1,IJU
!!$    IIMIN=MAX(JI-2,1)
!!$    IIMAX=MIN(JI+2,IIU)
!!$    IJMIN=MAX(JJ-2,1)
!!$    IJMAX=MIN(JJ+2,IJU)
!!$    IF (IWK_BL_TOP(JI,JJ) >=  SUM(IK_BL_TOP(IIMIN:IIMAX,IJMIN:IJMAX))  &
!!$                           / ((IIMAX-IIMIN+1)*(IJMAX-IJMIN+1)) +3 )    &
!!$    IWK_BL_TOP(JI,JJ)      =  SUM(IK_BL_TOP(IIMIN:IIMAX,IJMIN:IJMAX))  &
!!$                           / ((IIMAX-IIMIN+1)*(IJMAX-IJMIN+1))
!!$  END DO
!!$END DO
!!$IK_BL_TOP(:,:)=IWK_BL_TOP(:,:)

ZK_BL_TOP(:,:)=FLOAT(IK_BL_TOP(:,:))
CALL MPPDB_CHECK2D(ZK_BL_TOP,"FREE_ATM_PROFILE:8.2:ZK_BL_TOP",PRECISION)
!
!*       8.2  spatial filtering is applied (4 times) for boundary layer top
!             -------------------------------------------------------------
!
ZK_BL_TOP(:,:)=FLOAT(IK_BL_TOP(:,:))
CALL PGDFILTER(ZK_BL_TOP(:,:),4)
CALL MPPDB_CHECK2D(ZK_BL_TOP,"FREE_ATM_PROFILE:ZK_BL_TOP",PRECISION)
IK_BL_TOP(:,:)=NINT(ZK_BL_TOP(:,:))
!
!-------------------------------------------------------------------------------
!
!*       9.   Extrapolation of profile from top of the boundary layer
!             -------------------------------------------------------
!
!*       9.1  Gradient above boundary layer
!             -----------------------------
!
!* We use the points in free atmosphere up to 5000m
!
DO JI=1,IIU
  DO JJ=1,IJU
    IK5000=MAX(IKB+2,MIN(IKE-2,                                            &
           COUNT(PZMASS_MX(JI,JJ,:)<MIN(5000.+PZS_LS(JI,JJ),               &
                                        MAX(6000.,2000.+PZS_LS(JI,JJ))) )  ))
    IK3000=MAX(IKB+2,MIN(IKE-2,                                            &
           COUNT(PZMASS_MX(JI,JJ,:)<MIN(3000.+PZS_LS(JI,JJ),               &
                                        MAX(4000.,1000.+PZS_LS(JI,JJ))) )  ))
    IF (IK3000==IK5000) THEN
      IK3000=IK5000-1
    END IF
    IF (IK3000>IKB+1) THEN
      ZFREE_GR(JI,JJ) = (PVAR_MX  (JI,JJ,IK5000) - PVAR_MX  (JI,JJ,MIN(IK_BL_TOP(JI,JJ),IK3000)))&
                      / (PZMASS_MX(JI,JJ,IK5000) - PZMASS_MX(JI,JJ,MIN(IK_BL_TOP(JI,JJ),IK3000)))
    ELSE
      ZFREE_GR(JI,JJ) = XT_CLIM_GRAD
    END IF
  END DO
END DO
!
!*       9.2  spatial filtering is applied (8 times) for the gradient
!             -------------------------------------------------------
!
CALL PGDFILTER(ZFREE_GR(:,:),4)
!
!*       9.3  free atmosphere profile is computed
!             -----------------------------------
!
DO JI=1,IIU
  DO JJ=1,IJU
!
    ZF_FREE_MX(JI,JJ,IK_BL_TOP(JI,JJ):) = PVAR_MX(JI,JJ,IK_BL_TOP(JI,JJ):)
!
    ZF_FREE_MX(JI,JJ,:IK_BL_TOP(JI,JJ)-1) = ZF_FREE_MX(JI,JJ,IK_BL_TOP(JI,JJ)) &
                    + ZFREE_GR(JI,JJ) * (  PZMASS_MX(JI,JJ,:IK_BL_TOP(JI,JJ)-1)&
                                         - PZMASS_MX(JI,JJ, IK_BL_TOP(JI,JJ)  ) )
!
  END DO
END DO
!
!*       9.5  profile is modified in case of change of zs upwards
!             ---------------------------------------------------
!
!* We need to have the constant free atmosphere gradient also above the boundary
!  layer, in order to produce a correct shift. The added height where the
!  gradient apply is equal to the difference betweeen the two orographies.
!
DO JI=1,IIU
  DO JJ=1,IJU
    IF (PZS_LS(JI,JJ)<XZS(JI,JJ))                                           &
    IK_BL_TOP(JI,JJ) = COUNT(                                               &
                       PZMASS_MX(JI,JJ,:)<PZMASS_MX(JI,JJ,IK_BL_TOP(JI,JJ)) &
                                         +XZS(JI,JJ)-PZS_LS(JI,JJ)          )
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*      10.   Interpolation to the grid for the free profiles
!             -----------------------------------------------
!
!* With the following interpolation routine, the points under the surface
!  of the mixed grid are linearly extrapolated. Therefore, the gradient
!  of the free atmosphere is also kept for such points.
!
IF (CPROGRAM /= 'DIAG  ') THEN
  CALL COEF_VER_INTERP_LIN(PZMASS_MX(:,:,:),PZ_FREE(:,:,:),OLEUG=.TRUE.)
  PF_FREE(:,:,:)=VER_INTERP_LIN(ZF_FREE_MX(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
ELSE
  PF_FREE(:,:,:)=0.
END IF
!
!-------------------------------------------------------------------------------
!
!*      11.   Prints
!             ------
!
IF (CPROGRAM == 'DIAG  ' ) THEN
!
!*      11.1  Writing of height of boundary layer
!             -----------------------------------
  DO JI=1,IIU
    DO JJ=1,IJU
      Z2D(JI,JJ) = PZMASS_MX(JI,JJ,IK_BL_TOP(JI,JJ)) - PZS_LS(JI,JJ)
    END DO
  END DO
  TZFIELD%CMNHNAME   = 'HBLTOP'
  TZFIELD%CSTDNAME   = 'atmosphere_boundary_layer_thickness'
  TZFIELD%CLONGNAME  = 'HBLTOP'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Height of Boundary Layer TOP'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,Z2D)
!
!*      11.2  Writing of level of boundary layer top
!             --------------------------------------
!
  Z2D(:,:) = IK_BL_TOP(:,:)
  TZFIELD%CMNHNAME   = 'KBLTOP'
  TZFIELD%CSTDNAME   = 'model_level_number_at_top_of_atmosphere_boundary_layer'
  TZFIELD%CLONGNAME  = 'KBLTOP'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Index of Boundary Layer TOP'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,IK_BL_TOP)
END IF
!
IF (CPROGRAM /= 'DIAG  ' .AND. CPROGRAM /= 'IDEAL ' ) THEN
!
!*      11.3  Writing of free atmosphere gradient
!             -----------------------------------
!
  Z2D(:,:)=ZFREE_GR(:,:)
!
  TZFIELD%CMNHNAME   = 'FREE_ATM_GR'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'FREE_ATM_GR'
  TZFIELD%CUNITS     = 'K m-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'Free atmosphere gradient'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,Z2D)
!
!*      11.4  Writing of free atmosphere 3D profiles
!             --------------------------------------
!
  ZZMASS(:,:,1:IKE) = 0.5* XZZ(:,:,1:IKE) + 0.5* XZZ(:,:,2:IKE+1)
  ZZMASS(:,:,IKE+1) = 1.5* XZZ(:,:,IKE+1) - 0.5* XZZ(:,:,IKE)
!
  CALL COEF_VER_INTERP_LIN(PZ_FREE(:,:,:),ZZMASS(:,:,:),OLEUG=.TRUE.)
  Z3D(:,:,:)=VER_INTERP_LIN(PF_FREE(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
!
  TZFIELD%CMNHNAME   = 'THV_FREE'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'THV_FREE'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_THV_FREE'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,Z3D)
!
END IF
!-------------------------------------------------------------------------------
!
WRITE(ILUOUT0,*) 'Routine FREE_ATM_PROFILE completed'
!
END SUBROUTINE FREE_ATM_PROFILE
