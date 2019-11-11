!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################
MODULE MODI_SPAWN_GRID2
!######################
!
INTERFACE
!
     SUBROUTINE SPAWN_GRID2 (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,         &
                             PLONOR,PLATOR,PXHAT,PYHAT,PZHAT,PZTOP,           &
                             OSLEVE,PLEN1,PLEN2,                              &
                             PZS,PZSMT,PZS_LS,PZSMT_LS,                       &
                             TPDTMOD,TPDTCUR                                  )
!
USE MODD_TIME
!
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL,                 INTENT(INOUT) :: PLATOR            ! Latitude of the origine point
REAL,                 INTENT(INOUT) :: PLONOR            ! Longitude of the origine point
REAL, DIMENSION(:),   INTENT(INOUT) :: PXHAT,PYHAT,PZHAT ! positions x,y,z in the
                                     ! conformal plane or on the cartesian plane
REAL,                 INTENT(OUT)   :: PZTOP             ! model top (m)
LOGICAL,              INTENT(OUT)   :: OSLEVE            ! flag for SLEVE coordinate
REAL,                 INTENT(OUT)   :: PLEN1             ! Decay scale for smooth topography
REAL,                 INTENT(OUT)   :: PLEN2             ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:), INTENT(INOUT) :: PZS               ! orography
REAL, DIMENSION(:,:), INTENT(INOUT) :: PZSMT             ! smooth orography
REAL, DIMENSION(:,:), INTENT(OUT)   :: PZS_LS            ! interpolated orography
REAL, DIMENSION(:,:), INTENT(OUT)   :: PZSMT_LS          ! interpolated smooth orography
!
!
TYPE (DATE_TIME),     INTENT(INOUT) :: TPDTMOD  ! Date and Time of MODel beginning
TYPE (DATE_TIME),     INTENT(INOUT) :: TPDTCUR  ! CURent date and time
!
END SUBROUTINE SPAWN_GRID2
!
END INTERFACE
!
END MODULE MODI_SPAWN_GRID2
!
!
!     #########################################################################
     SUBROUTINE SPAWN_GRID2 (KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,         &
                             PLONOR,PLATOR,PXHAT,PYHAT,PZHAT,PZTOP,           &
                             OSLEVE,PLEN1,PLEN2,                              &
                             PZS,PZSMT,PZS_LS,PZSMT_LS,                       &
                             TPDTMOD,TPDTCUR                                  )
!     #########################################################################
!
!!****  *SPAWN_GRID2 * - subroutine to define spatial and temporal grid.
!!
!!    PURPOSE
!!    -------
!!
!!      This routine defines the information necessary to generate the model 2
!!    grid, consistently with the spawning model 1.
!!      The longitude and latitude of the model 2 origine are computed from
!!    the model 1. Then the grid in the conformal projection and terrain
!!    following coordinates (XHAT,YHAT and ZHAT) and orography, are interpolated
!!    from the model 1 grid and orography knowledge.
!!
!!      Date and time are set as for model 1.
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
!!    grid and orography are interpolated as follows:
!!         - linear interpolation for XHAT and YHAT
!!         - identity for ZHAT (no vertical spawning)
!!            2 types of interpolations can be used:
!!                 1. Clark and Farley (JAS 1984) on 9 points
!!                 2. Bikhardt on 16 points
!!
!!    EXTERNAL
!!    --------
!!
!!      Module MODE_TIME : contains SM_PRINT_TIME routine
!!      Routine BIKHARDT2     : to perform horizontal interpolations
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS : contains parameters
!!      Module MODD_CONF       : contains models configuration
!!      Module MODD_GRID1      : contains grid variables
!!      Module MODD_TIME1      : contains date and time of model 1
!!                              and uses MODD_TIME
!!      Module MODD_GR_FIELD1  : contains surface variables
!!
!!      Module MODD_LUNIT2     : contains unit numbers of model 2 files
!!
!!    REFERENCE
!!    ---------
!!
!!       Book1 of the documentation
!!       PROGRAM SPAWN_GRID2 (Book2 of the documentation)
!!
!!
!!    AUTHOR
!!    ------
!!
!!       J.P. Lafore     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     12/01/95
!!      Modification 05/07/95  (Lafore) Different resolution ratio case introduction
!!      Modification 31/01/96  (Lafore) Update for MASDEV2_2 version and corrections
!!                                      for the different resolution ratio case
!!      Modification 19/02/96  (Lafore) introduction of the Bikhardt interpolation
!!      Modification 19/03/96  (Lafore) interpolation of surface variables
!!      Modification 10/06/96 (V.Masson) remove the loops in case of no resolution change
!!      Modification 10/06/96 (V.Masson) interpolation computations performed in
!!                                       independant routines
!!      Modification 19/06/96 (V.Masson) case of integer input land sea mask
!!      Modification 02/10/96 (V.Masson) iterative method for zs computation
!!      Modification 21/11/96 (Lafore)   move from BIKHARDT2 to BIKHARDT routine
!!      Modification 16/07/97 (V.Masson) bug in test of positivity for zs
!!      Modification 17/07/97 (V.Masson) purely interpolated zs (PZS_LS)
!!      Modification 10/10/97 (V.Masson) bug on boundaries for zs procedure
!!      Modification 20/04/99 (J. Stein) bug on the last point if the whole
!!                             domain is used (2D case along y for instance
!!      Modification 15/03/99 (V.Masson) cover types
!!      Modification 04/07/01 (J.Stein)  convergence test set to 1 millimeter for XZS1
!!      Modification 05/09/05 (J. Escobar) change INTENT(OUT) --> INTENT(INOUT)
!!                             to avoid problem when Input parameter and GRID1 parameter
!!                             are exactly the same !!!
!!      Modification 20/05/06 Remove Clark and Farley interpolation
!!      Modification 24/02/15 (M.Moge) parallelization
!!      Modification 10/06/15 (M.Moge) bug fix for reproductibility
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      J.Escobar 05/03/2018 : bypass gridnesting special case KD(X/Y)RATIO == 1 not parallelized
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_PARAMETERS       ! Declarative modules
USE MODD_CONF
USE MODD_SPAWN
!
USE MODD_GRID, ONLY: XLONORI,XLATORI 
USE MODD_LBC_n,     ONLY: LBC_MODEL
!
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_BIKHARDT_n
USE MODD_VAR_ll
USE MODE_ll
USE MODE_FM
USE MODE_IO_ll
USE MODE_TIME
USE MODE_GRIDPROJ
!
USE MODI_BIKHARDT
USE MODI_SPAWN_ZS
!
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
!
REAL,                 INTENT(INOUT) :: PLATOR            ! Latitude of the origine point
REAL,                 INTENT(INOUT) :: PLONOR            ! Longitude of the origine point
REAL, DIMENSION(:),   INTENT(INOUT) :: PXHAT,PYHAT,PZHAT ! positions x,y,z in the
                                     ! conformal plane or on the cartesian plane
REAL,                 INTENT(OUT)   :: PZTOP             ! model top (m)
LOGICAL,              INTENT(OUT)   :: OSLEVE            ! flag for SLEVE coordinate
REAL,                 INTENT(OUT)   :: PLEN1             ! Decay scale for smooth topography
REAL,                 INTENT(OUT)   :: PLEN2             ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:), INTENT(OUT)   :: PZS               ! orography
REAL, DIMENSION(:,:), INTENT(OUT)   :: PZSMT             ! smooth orography
REAL, DIMENSION(:,:), INTENT(OUT)   :: PZS_LS            ! interpolated orography
REAL, DIMENSION(:,:), INTENT(OUT)   :: PZSMT_LS          ! interpolated smooth orography
!
!
TYPE (DATE_TIME),     INTENT(INOUT) :: TPDTMOD  ! Date and Time of MODel beginning
TYPE (DATE_TIME),     INTENT(INOUT) :: TPDTCUR  ! CURent date and time
!
!*       0.2    Declarations of local variables for print on FM file
!
INTEGER :: ILUOUT   ! Logical unit number for the output listing
INTEGER :: IRESP    ! Return codes in FM routines
!
REAL :: ZPOND1,ZPOND2               ! interpolation coefficients
INTEGER             :: IIU_C       ! Upper dimension in x direction
INTEGER             :: IJU_C       ! Upper dimension in y direction
INTEGER             :: IIB_C       ! indice I Beginning in x direction
INTEGER             :: IJB_C       ! indice J Beginning in y direction
!
INTEGER             :: IIU       ! Upper dimension in x direction
INTEGER             :: IJU       ! Upper dimension in y direction
INTEGER             :: IIB,IIE   ! indice I Beginning/End in x direction
INTEGER             :: IJB,IJE   ! indice J Beginning/End in y direction
INTEGER             :: IIS,IJS   ! indices I and J in x and y dir. for scalars
INTEGER             :: JI,JEPSX  ! Loop index in x direction
INTEGER             :: JJ,JEPSY  ! Loop index in y direction
REAL, DIMENSION(:), ALLOCATABLE :: ZXHAT_EXTENDED, ZYHAT_EXTENDED
INTEGER             :: IXSIZE1_F,IYSIZE1_F    ! sizes of the XHAT and YHAT arrays
!
CHARACTER (LEN=40)  :: YTITLE    ! Title for time print
INTEGER             :: IMI
INTEGER             :: IINFO_ll
INTEGER             :: IXOR_F, IYOR_F, IXEND_F, IYEND_F
INTEGER             :: IXOR_ll, IYOR_ll
INTEGER             :: IXDIM, IYDIM
REAL, DIMENSION(1)  :: PXMAX, PYMAX, PXMIN, PYMIN
INTEGER             :: DELTA_JI,JI_MIN,JI_MAX,  DELTA_JJ,JJ_MIN,JJ_MAX
REAL                :: ZMIN
INTEGER             :: IDIMX_C, IDIMY_C
REAL, DIMENSION(:,:), ALLOCATABLE :: ZXHAT_2D_EXTENDED_F, ZYHAT_2D_EXTENDED_F
REAL, DIMENSION(:), ALLOCATABLE :: ZXHAT_EXTENDED_C, ZYHAT_EXTENDED_C
REAL, DIMENSION(:,:), ALLOCATABLE :: ZXHAT_2D_C, ZYHAT_2D_C
REAL, DIMENSION(:,:), ALLOCATABLE :: ZXHAT_2D_F, ZYHAT_2D_F
LOGICAL             :: GCYCLIC_EXTRAPOL
!-------------------------------------------------------------------------------
!
!
!*       1.    PROLOGUE:
!              ---------
! 
IMI = GET_CURRENT_MODEL_INDEX()
CALL GOTO_MODEL(2)
!
!*       1.1   Interpolation method
!
!
!*       1.1   computes dimensions of arrays and other indices
!
IIU_C = SIZE(PXHAT)
IJU_C = SIZE(PYHAT)
IIB_C = 1+JPHEXT
IJB_C = 1+JPHEXT
!
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
CALL GET_FEEDBACK_COORD_ll(IXOR_F,IYOR_F,IXEND_F,IYEND_F,IINFO_ll)
!
CALL GO_TOMODEL_ll(1,IINFO_ll)
CALL GET_OR_ll('B',IXOR_ll,IYOR_ll)
CALL GET_DIM_EXT_ll('B',IXDIM,IYDIM)
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
!
IF (IXOR_F>0 .and. IYOR_F>0 .and. &
    IXEND_F>0 .and. IYEND_F>0) THEN
   IXOR_F = IXOR_F-JPHEXT
   IYOR_F = IYOR_F-JPHEXT
   IXEND_F= IXEND_F+JPHEXT
   IYEND_F= IYEND_F+JPHEXT
ELSE
   IXOR_F = 1!4!2
   IXEND_F= 1!4!10
   IYOR_F = -10!4!2
   IYEND_F= -10!4!10
ENDIF
!$
!
!*       1.2  recovers logical unit number of output listing
!
ILUOUT = TLUOUT%NLU
!
!*       1.3  checks that model 2 domain is included in the one of model 1
IF ( (IXEND_F) > SIZE(XXHAT1) )  THEN   
  WRITE(ILUOUT,FMT=*) 'SPAWN_MODEL2:  MODEL 2 DOMAIN OUTSIDE THE MODEL1 DOMAIN  ',  &
                  ' IXOR_F = ', IXOR_F,' IXEND_F = ', IXEND_F,                      &
                  ' IIU of model1 = ',SIZE(XXHAT1)
  !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SPAWN_GRID2','MODEL 2 DOMAIN OUTSIDE THE MODEL1 DOMAIN')
END IF 
IF ( (IYEND_F) > SIZE(XYHAT1) )  THEN  
  WRITE(ILUOUT,FMT=*) 'SPAWN_MODEL2:  MODEL 2 DOMAIN OUTSIDE THE MODEL1 DOMAIN  ',  &
                  ' IYOR_F = ', IYOR_F,' IYEND_F = ', IYEND_F,                  &
                  ' IJU of model1 = ',SIZE(XYHAT1)
  !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SPAWN_GRID2','MODEL 2 DOMAIN OUTSIDE THE MODEL1 DOMAIN')
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZATION OF THE GRID OF MODEL 2:
!              --------------------------------------
!
PZTOP    = XZTOP1
PZHAT(:) = XZHAT1(:) 
OSLEVE   = LSLEVE1
PLEN1    = XLEN11
PLEN2    = XLEN21
!
!!$IF (KDXRATIO == 1 .AND. KDYRATIO == 1 ) THEN
!!$!
!!$!*       2.1   special case of spawning - no change of resolution :
!!$!$ in our case we don't get them here !
!!$  PXHAT(:) = GRID_MODEL(1)%XXHAT(KXOR:KXEND)
!!$  PYHAT(:) = GRID_MODEL(1)%XYHAT(KYOR:KYEND)
!!$  PZS  (:,:) = GRID_MODEL(1)%XZS  (KXOR:KXEND,KYOR:KYEND)
!!$  PZS_LS(:,:)= GRID_MODEL(1)%XZS  (KXOR:KXEND,KYOR:KYEND)
!!$  PZSMT   (:,:) = GRID_MODEL(1)%XZSMT(KXOR:KXEND,KYOR:KYEND)
!!$  PZSMT_LS(:,:) = GRID_MODEL(1)%XZSMT(KXOR:KXEND,KYOR:KYEND)
!!$!
!!$ELSE
!
!*       2.2  general case - change of resolution :
!
!*       2.2.1 linear interpolation for XHAT and YHAT
  GCYCLIC_EXTRAPOL = .FALSE.
!
!     XHAT
!
!JUAN A REVOIR TODO_JPHEXT
! <<<<<<< spawn_grid2.f90
  IXSIZE1_F=SIZE(XXHAT1)
  IYSIZE1_F=SIZE(XYHAT1)
! before the interpolation of XXHAT into PXHAT, we need to use LS_FORCING_ll
! to communicate the values on the subdomains of the son grid to the appropriate processes
! LS_FORCING_ll does not work on 1D arrays, so we have to construct a temporary pseudo-2D array
  ALLOCATE(ZXHAT_2D_F(IXSIZE1_F,IYSIZE1_F))
  ZXHAT_2D_F(:,:) = SPREAD(XXHAT1(:),DIM=2,NCOPIES=IYSIZE1_F)
  CALL GOTO_MODEL(1)
  CALL GO_TOMODEL_ll(1, IINFO_ll)
  CALL GET_CHILD_DIM_ll(IMI, IDIMX_C, IDIMY_C, IINFO_ll)
  !allocation of the 1D and pseudo-2D arrays on child grid
  ALLOCATE(ZXHAT_EXTENDED_C(IDIMX_C+1))
  ALLOCATE(ZXHAT_2D_C(IDIMX_C,IDIMY_C))
  CALL SET_LSFIELD_1WAY_ll(ZXHAT_2D_F, ZXHAT_2D_C, IMI)
  CALL LS_FORCING_ll(IMI, IINFO_ll,.TRUE.,GCYCLIC_EXTRAPOL)
  CALL GO_TOMODEL_ll(IMI, IINFO_ll)
  CALL GOTO_MODEL(IMI)
  CALL UNSET_LSFIELD_1WAY_ll()
! initialization of ZXHAT_EXTENDED_C
! Remark : we take the 2nd row of ZXHAT_2D_C because the first one is the "pseudo halo" added for spawning
!          and may be uninitialized or an extrapolation of the second row
  ZXHAT_EXTENDED_C(1:IDIMX_C)=ZXHAT_2D_C(:,2)
! extrapolation on the extra point
  ZXHAT_EXTENDED_C(IDIMX_C+1)= 2.*ZXHAT_EXTENDED_C(IDIMX_C)-ZXHAT_EXTENDED_C(IDIMX_C-1)     !TODO : faire un update_nhalo1D
! interpolation on the child grid
  PXHAT(:)=0.
  !on the west halo of the son model
  DO JI = 1,JPHEXT
    DO JEPSX=1,KDXRATIO
      ZPOND2 = FLOAT(KDXRATIO-JEPSX)/FLOAT(KDXRATIO)
      ZPOND1 = 1.-ZPOND2
      IF( JPHEXT+1-(JI-1)*KDXRATIO-JEPSX > 0 ) THEN
        PXHAT(JPHEXT+1-(JI-1)*KDXRATIO-JEPSX) = ZPOND1*ZXHAT_EXTENDED_C(JPHEXT+1-JI+1) &
                      + ZPOND2*ZXHAT_EXTENDED_C(JPHEXT+1-JI+2)
      ENDIF
    ENDDO
  ENDDO
  !on the physical domain of the son model
  DO JI = 1,IDIMX_C-2*(JPHEXT+1)  !the physical size of the son model in the father grid
    DO JEPSX = 1,KDXRATIO
      ZPOND2 = FLOAT(JEPSX-1)/FLOAT(KDXRATIO)
      ZPOND1 = 1.-ZPOND2
      PXHAT(JPHEXT+JEPSX+(JI-1)*KDXRATIO) = ZPOND1*ZXHAT_EXTENDED_C(JI+IIB_C) & 
            + ZPOND2*ZXHAT_EXTENDED_C(JI+IIB_C+1)
    ENDDO
  ENDDO
  !on the east halo of the son model
  DO JI = 1,JPHEXT
    DO JEPSX=1,KDXRATIO
      ZPOND1 = FLOAT(KDXRATIO-JEPSX+1)/FLOAT(KDXRATIO)
      ZPOND2 = 1.-ZPOND1
      IF( SIZE(PXHAT)-JPHEXT+(JI-1)*KDXRATIO+JEPSX <= SIZE(PXHAT) ) THEN
        PXHAT(SIZE(PXHAT)-JPHEXT+(JI-1)*KDXRATIO+JEPSX) = ZPOND1*ZXHAT_EXTENDED_C(IDIMX_C-JPHEXT+JI-1) &
            + ZPOND2*ZXHAT_EXTENDED_C(IDIMX_C-JPHEXT+JI)
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(ZXHAT_2D_F)
  DEALLOCATE(ZXHAT_EXTENDED_C)
  DEALLOCATE(ZXHAT_2D_C)
!
!     YHAT
!
! before the interpolation of XXHAT into PXHAT, we need to use LS_FORCING_ll
! to communicate the values on the subdomains of the son grid to the appropriate processes
! LS_FORCING_ll does not work on 1D arrays, so we have to construct a temporary pseudo-2D array
  ALLOCATE(ZYHAT_2D_F(IXSIZE1_F,IYSIZE1_F))
  ZYHAT_2D_F(:,:) = SPREAD(XYHAT1(:),DIM=1,NCOPIES=IXSIZE1_F)
  CALL GOTO_MODEL(1)
  CALL GO_TOMODEL_ll(1, IINFO_ll)
  CALL GET_CHILD_DIM_ll(IMI, IDIMX_C, IDIMY_C, IINFO_ll)
  !allocation of the 1D and pseudo-2D arrays on child grid
  ALLOCATE(ZYHAT_EXTENDED_C(IDIMY_C+1))
  ALLOCATE(ZYHAT_2D_C(IDIMX_C,IDIMY_C))
  CALL SET_LSFIELD_1WAY_ll(ZYHAT_2D_F, ZYHAT_2D_C, IMI)
  CALL LS_FORCING_ll(IMI, IINFO_ll,.TRUE.,GCYCLIC_EXTRAPOL)
  CALL GO_TOMODEL_ll(IMI, IINFO_ll)
  CALL GOTO_MODEL(IMI)
  CALL UNSET_LSFIELD_1WAY_ll()
! initialization of ZXHAT_EXTENDED_C
  ZYHAT_EXTENDED_C(1:IDIMY_C)=ZYHAT_2D_C(1,:)
! extrapolation on the extra point
  ZYHAT_EXTENDED_C(IDIMY_C+1)= 2.*ZYHAT_EXTENDED_C(IDIMY_C)-ZYHAT_EXTENDED_C(IDIMY_C-1)
  PYHAT(:)=0.
  !on the south halo of the son model
  DO JJ = 1,JPHEXT
    DO JEPSY=1,KDYRATIO
      ZPOND2 = FLOAT(KDXRATIO-JEPSY)/FLOAT(KDYRATIO)
      ZPOND1 = 1.-ZPOND2
      IF( JPHEXT+1-(JJ-1)*KDYRATIO-JEPSY > 0 ) THEN
        PYHAT(JPHEXT+1-(JJ-1)*KDYRATIO-JEPSY) = ZPOND1*ZYHAT_EXTENDED_C(JPHEXT+1-JJ+1) &
             + ZPOND2*ZYHAT_EXTENDED_C(JPHEXT+1-JJ+2)
      ENDIF
    ENDDO
  ENDDO
  !on the physical domain of the son model
  DO JJ = 1,IDIMY_C-2*(JPHEXT+1)  !the physical size of the son model in the father grid
    DO JEPSY = 1,KDYRATIO
      ZPOND2 = FLOAT(JEPSY-1)/FLOAT(KDYRATIO)
      ZPOND1 = 1.-ZPOND2
      PYHAT(JPHEXT+JEPSY+(JJ-1)*KDYRATIO) = ZPOND1*ZYHAT_EXTENDED_C(JJ+JPHEXT+1) &
            + ZPOND2*ZYHAT_EXTENDED_C(JJ+JPHEXT+1+1)
    ENDDO
  ENDDO
  !on the north halo of the son model
  DO JJ = 1,JPHEXT
    DO JEPSY=1,KDYRATIO
      ZPOND1 = FLOAT(KDYRATIO-JEPSY+1)/FLOAT(KDYRATIO)
      ZPOND2 = 1.-ZPOND1
      IF( SIZE(PYHAT)-JPHEXT+(JJ-1)*KDYRATIO+JEPSY <= SIZE(PYHAT) ) THEN
        PYHAT(SIZE(PYHAT)-JPHEXT+(JJ-1)*KDYRATIO+JEPSY) = ZPOND1*ZYHAT_EXTENDED_C(IDIMY_C-JPHEXT+JJ-1) &
             + ZPOND2*ZYHAT_EXTENDED_C(IDIMY_C-JPHEXT+JJ)
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(ZYHAT_2D_F)
  DEALLOCATE(ZYHAT_EXTENDED_C)
  DEALLOCATE(ZYHAT_2D_C)
!!$=======
!!$  IXSIZE1=SIZE(XXHAT1)
!!$  ALLOCATE(ZXHAT_EXTENDED(IXSIZE1+1))
!!$  ZXHAT_EXTENDED(1:IXSIZE1)=XXHAT1(:)
!!$  ZXHAT_EXTENDED(IXSIZE1+1)=2.*XXHAT1(IXSIZE1)-XXHAT1(IXSIZE1-1)
!!$  DO JEPSX = 1,KDXRATIO
!!$    ZPOND2 = FLOAT(JEPSX-1)/FLOAT(KDXRATIO)
!!$    ZPOND1 = 1.-ZPOND2
!!$    DO JI = KXOR,KXEND
!!$      IIS = IIB+JEPSX-1+(JI-KXOR-JPHEXT)*KDXRATIO
!!$!
!!$      IF (1 <= IIS .AND. IIS <= IIU)                   &
!!$      PXHAT(IIS) = ZPOND1*ZXHAT_EXTENDED(JI) +ZPOND2*ZXHAT_EXTENDED(JI+1)
!!$    END DO
!!$  END DO
!!$  DEALLOCATE(ZXHAT_EXTENDED)
!!$!
!!$  IYSIZE1=SIZE(XYHAT1)
!!$  ALLOCATE(ZYHAT_EXTENDED(IYSIZE1+1))
!!$  ZYHAT_EXTENDED(1:IYSIZE1)=XYHAT1(:)
!!$  ZYHAT_EXTENDED(IYSIZE1+1)=2.*XYHAT1(IYSIZE1)-XYHAT1(IYSIZE1-1)
!!$  DO JEPSY = 1,KDYRATIO
!!$    ZPOND2 = FLOAT(JEPSY-1)/FLOAT(KDYRATIO)
!!$    ZPOND1 = 1.-ZPOND2
!!$    DO JJ = KYOR,KYEND
!!$      IJS = IJB+JEPSY-1+(JJ-KYOR-JPHEXT)*KDYRATIO
!!$!
!!$      IF (1 <= IJS .AND. IJS <= IJU)                   &
!!$      PYHAT(IJS) = ZPOND1*ZYHAT_EXTENDED(JJ) +ZPOND2*ZYHAT_EXTENDED(JJ+1)
!!$    END DO
!!$  END DO
!!$  DEALLOCATE(ZYHAT_EXTENDED)
!!$>>>>>>> 1.3.4.2.18.2.2.1
!
!
!*       2.2.2  interpolation of ZS performed later
!
!!$END IF
!
PLONOR = XLONORI
PLATOR = XLATORI
!
!-------------------------------------------------------------------------------
!
!*       3.    INITIALIZATION OF ZS and ZSMT:
!              ------------------------------
CALL SPAWN_ZS(IXOR_F,IXEND_F,IYOR_F,IYEND_F,KDXRATIO,KDYRATIO,IDIMX_C,IDIMY_C,  &
              LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,XZS1,  PZS,  'ZS    ',PZS_LS)
CALL SPAWN_ZS(IXOR_F,IXEND_F,IYOR_F,IYEND_F,KDXRATIO,KDYRATIO,IDIMX_C,IDIMY_C,    &
              LBC_MODEL(1)%CLBCX,LBC_MODEL(1)%CLBCY,XZSMT1,PZSMT,'ZSMT  ',PZSMT_LS)
!
CALL MPPDB_CHECK2D(PZS,"SPAWN_GRID2:PZS",PRECISION)
CALL MPPDB_CHECK2D(PZSMT,"SPAWN_GRID2:PZSMT",PRECISION)
!$
!-------------------------------------------------------------------------------
!
!*       4.    INITIALIZATION OF MODEL 2 DATE AND TIME:
!              ----------------------------------------
!
TPDTMOD = TDTCUR1
TPDTCUR = TDTCUR1
!
YTITLE='OUTER MODEL : CURRENT DATE AND TIME '
CALL SM_PRINT_TIME(TDTCUR1, TLUOUT, YTITLE)
YTITLE='SPAWNED MODEL : DATE AND TIME BEGINNING'
CALL SM_PRINT_TIME(TPDTMOD, TLUOUT, YTITLE)
YTITLE='SPAWNED MODEL : CURRENT DATE AND TIME '
CALL SM_PRINT_TIME(TPDTCUR, TLUOUT, YTITLE)
!
!-------------------------------------------------------------------------------
CALL GOTO_MODEL(IMI)
!
END SUBROUTINE SPAWN_GRID2
