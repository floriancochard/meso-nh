!MNH_LIC Copyright 2005-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!###################
MODULE MODI_SPAWN_ZS
!###################
!
INTERFACE
!
     SUBROUTINE SPAWN_ZS (KXOR,KXEND,KYOR,KYEND,KDXRATIO,KDYRATIO,KDIMX_C,KDIMY_C,&
                          HLBCX,HLBCY,PZS1_F,PZS2_C,HFIELD,PZS2_LS                )
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
INTEGER,   INTENT(IN)  :: KDIMX_C    ! dimension (X dir) of local son subdomain in father grid
INTEGER,   INTENT(IN)  :: KDIMY_C    ! dimension (Y dir) of local son subdomain in father grid
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
REAL, DIMENSION(:,:), INTENT(IN)  :: PZS1_F    ! model 1 orography
REAL, DIMENSION(:,:), INTENT(OUT) :: PZS2_C    ! interpolated orography with iterative correction
CHARACTER(LEN=6),     INTENT(IN)  :: HFIELD ! name of the field to nest
REAL, DIMENSION(:,:), INTENT(OUT),OPTIONAL  :: PZS2_LS ! interpolated orography
!
END SUBROUTINE SPAWN_ZS
!
END INTERFACE
!
END MODULE MODI_SPAWN_ZS
!
!
!     #############################################################################
     SUBROUTINE SPAWN_ZS (KXOR,KXEND,KYOR,KYEND,KDXRATIO,KDYRATIO,KDIMX_C,KDIMY_C,&
                          HLBCX,HLBCY,PZS1_F,PZS2_C,HFIELD,PZS2_LS                )
!     #############################################################################
!
!!****  *SPAWN_ZS * - subroutine to spawn zs field
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
!!      Routine BIKHARDT2     : to perform horizontal interpolations
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS : contains parameters
!!      Module MODD_CONF       : contains models configuration
!!      Module MODD_LUNIT2     : contains unit numbers of model 2 files
!!
!!    REFERENCE
!!    ---------
!!
!!       Book1 of the documentation
!!       PROGRAM SPAWN_ZS (Book2 of the documentation)
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson    * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     12/01/05
!      Modification    20/05/06 Remove Clark and Farley interpolation
!      Modification    2014 M.Faivre : parallelizattion attempt
!      Modification    10/02/15 M. Moge : paralellization
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_PARAMETERS, ONLY : JPHEXT       ! Declarative modules
USE MODD_CONF,       ONLY : NVERB
USE MODD_LUNIT_n,    ONLY: TLUOUT
!
USE MODD_BIKHARDT_n
!
USE MODI_BIKHARDT
USE MODI_ZS_BOUNDARY
!
USE MODE_MODELN_HANDLER
!
USE MODE_FM
USE MODE_MPPDB
USE MODD_VAR_ll
USE MODE_ll
USE MODD_LBC_n
USE MODD_NESTING
USE MODE_EXCHANGE_ll
USE MODE_EXTRAPOL
!
IMPLICIT NONE
!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$20140624 renaming for VARS :
!  frame=Father -> _F when DIMS = IDIMX,Y
!  projection from |grid1 to |grid2 :
!  obtained with SET_LSFIELD_1WAYn + LS_FORCING 
!  frame=Son    -> _C when DIMS = IOR,END
!  projection from |grid2 to |grid1 :
!  obtained with SET_LSFIELD_2WAYn + LS_FEEDBACK
!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,   INTENT(IN)  :: KXOR,KXEND !  horizontal position (i,j) of the ORigin and END
INTEGER,   INTENT(IN)  :: KYOR,KYEND ! of the model 2 domain, relative to model 1
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction Resolution ratio
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between model 2 and model 1
INTEGER,   INTENT(IN)  :: KDIMX_C    ! dimension (X dir) of local son subdomain in father grid
INTEGER,   INTENT(IN)  :: KDIMY_C    ! dimension (Y dir) of local son subdomain in father grid
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
REAL, DIMENSION(:,:), INTENT(IN)  :: PZS1_F    ! model 1 orography
REAL, DIMENSION(:,:), INTENT(OUT) :: PZS2_C    ! interpolated orography with iterative correction
CHARACTER(LEN=6),     INTENT(IN)  :: HFIELD ! name of the field to nest
REAL, DIMENSION(:,:), INTENT(OUT),OPTIONAL  :: PZS2_LS ! interpolated orography
!
!*       0.2    Declarations of local variables for print on FM file
!
INTEGER :: ILUOUT   ! Logical unit number for the output listing
INTEGER :: IRESP    ! Return codes in FM routines
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS2_LS ! interpolated orography
REAL, DIMENSION(:,:), ALLOCATABLE :: PZS1_C ! model 1 orography resticted to the grid of model 2
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS1_C ! zs of model 1 at iteration n or n+1 in GRID2
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS2_C ! averaged zs of model 2 at iteration n
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDZS_C ! difference between PZS1 and ZZS2
!$20140617 ZTZS1 result of SET_LSFIELD_1WAY_ll(PZS1)
!JUAN REAL, DIMENSION(:,:), ALLOCATABLE :: ZTZS1_C
!$20140704 ZDZS_3D to use MAX_ll(array3D arg)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDZS_3D
!$
INTEGER                :: IXMIN, IXMAX ! indices to interpolate the
                                    ! modified orography on model 2
                                    ! domain to model 1 grid
!JUAN A REVOIR TODO_JPHEXT /!\ /!\
! <<<<<<< spawn_zs.f90
INTEGER                :: JI,JEPSX  ! Loop index in x direction
INTEGER                :: JJ,JEPSY  ! Loop index in y direction
INTEGER                :: JCOUNTER  ! counter for iterative method
REAL                   :: ZRELAX             ! relaxation factor
INTEGER                :: JMAXITER = 2000    ! maximum number of iterations
!
INTEGER, DIMENSION(2)  :: IZSMAX
INTEGER                :: IMI     ! current model index
!$20140604
INTEGER                :: KMI,IDIMX_C,IDIMY_C
!$20140602
INTEGER                :: PZS1_FSIZE1
INTEGER                :: PZS1_FSIZE2
!$20140603
INTEGER                :: IINFO_ll
!$20140619
TYPE(LIST_ll), POINTER :: TZFIELDS_ll => NULL()  ! list of fields to exchange
!$20140623
INTEGER                :: IXOR_F,IXEND_F
INTEGER                :: IYOR_F,IYEND_F
INTEGER                :: KDXRATIO_C, KDYRATIO_C
!$20140704
!$20140711 not INT, REAL !!
REAL                   :: ZMAXVAL
REAL                   :: LOCMAXVAL
!$20140801
INTEGER                :: IORX, IORY, IIBINT,IJBINT,IIEINT,IJEINT
INTEGER                :: IXOR_C_ll, IXEND_C_ll  ! origin and end of the local subdomain of the child model 2
INTEGER                :: IYOR_C_ll, IYEND_C_ll  ! relative to the father model 1
INTEGER                :: IINFO  ! return code of // routines
!
TYPE(LIST_ll), POINTER :: TZZSFIELD_ll   ! list of fields to exchange
TYPE(HALO2LIST_ll), POINTER :: TZZSHALO2_ll   ! needed for update_halo2_ll

REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS1CHILDGRID_C  ! copy of ZZS1_C extended to the whole child domain
INTEGER                :: JI2INF,JI2SUP
INTEGER                :: JJ2INF,JJ2SUP
INTEGER               :: IXSIZE,IYSIZE
!
KDXRATIO_C=KDXRATIO
KDYRATIO_C=KDYRATIO
ZMAXVAL=1000.
!$
!!$=======
!!$INTEGER             :: JI,JEPSX  ! Loop index in x direction
!!$INTEGER             :: JJ,JEPSY  ! Loop index in y direction
!!$INTEGER             :: JCOUNTER  ! counter for iterative method
!!$REAL                :: ZRELAX             ! relaxation factor
!!$INTEGER             :: JMAXITER = 2000    ! maximum number of iterations
!!$!
!!$INTEGER, DIMENSION(2) :: IZSMAX
!!$INTEGER               :: IMI     ! current model index
!!$!
!!$INTEGER               :: IXSIZE,IYSIZE
!!$INTEGER                          :: INFO_ll                ! error return code
!!$>>>>>>> 1.1.4.1.18.2.2.1
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE:
!              ---------
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL GOTO_MODEL(IMI)
!
!
!*       1.2  recovers logical unit number of output listing
!
ILUOUT = TLUOUT%NLU
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTATION:
!              ---------
!
!*       2.1   Purely interpolated zs:
!              -----------------------
!
ALLOCATE(ZZS2_LS(SIZE(PZS2_C,1),SIZE(PZS2_C,2)))
PZS1_FSIZE1=SIZE(PZS1_F,1)
PZS1_FSIZE2=SIZE(PZS1_F,2)
!
!
!
! This is one way to do it, but then the compute load are not well balanced
! Each process computes the interpolation on the intersection of the global
! child model with its part of the father model
!CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
!               XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
!               KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,1,       &
!               HLBCX,HLBCY,PZS1_F,ZZS2_LS)
!
! We want instead that each process computes the interpolation for its part
! of the child model
! Then we have to communicate the values of PZS1_F on each subdomain of the
! child model to the corresponding process
! before calling BIKHARDT on the local subdomain of the child model
!
!*      1 GATHER LS FIELD FOR THE CHILD MODEL KMI
!       1.1  Must be on the father model to call get_child_dim
CALL GOTO_MODEL(NDAD(IMI))
!$20140623 KMI is DAD, IMI=son !!
!$20140623 use IMI not KMI
CALL GO_TOMODEL_ll(NDAD(IMI), IINFO_ll)
IDIMX_C = KDIMX_C! + 2*(JPHEXT+1) !KXEND-KXOR-1
IDIMY_C = KDIMY_C! + 2*(JPHEXT+1) !KYEND-KYOR-1
!CALL GET_CHILD_DIM_ll(IMI, IDIMX_C, IDIMY_C, IINFO_ll)
!
!         1.3  Specify the ls "source" fields and receiver fields
!
ALLOCATE(ZZS1_C(IDIMX_C,IDIMY_C))
ZZS1_C(:,:)=0.
CALL SET_LSFIELD_1WAY_ll(PZS1_F, ZZS1_C, IMI)
CALL MPPDB_CHECK2D(PZS1_F,"SPAWN_ZS:PZS1_F",PRECISION)
!        1.4  Communication
CALL LS_FORCING_ll(IMI, IINFO_ll, .TRUE.)
!        1.5  Back to the (current) child model
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
CALL GOTO_MODEL(IMI)
CALL UNSET_LSFIELD_1WAY_ll()
!
!if the child grid is the whole father grid, we first need to extrapolate
!the data on a "pseudo halo" before doing BIKHARDT interpolation
!CALL EXTRAPOL_ON_PSEUDO_HALO(ZZS1_C)
! <<<<<<< spawn_zs.f90
CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
               XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
               2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO_C,KDYRATIO_C,1, &
               HLBCX,HLBCY,ZZS1_C,ZZS2_LS)
!!$=======
!
!!$  CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
!!$                 XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
!!$                 KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,1,       &
!!$                 HLBCX,HLBCY,PZS1,ZZS2_LS)
!!$>>>>>>> 1.1.4.1.18.2.2.1
  CALL MPPDB_CHECK2D(ZZS2_LS,"SPAWN_ZS::ZZS2_LS",PRECISION)
!
!*       4.2   New zs:
!              -------
!
!* we use an iterative method to insure the equality between large-scale
!  orography and the average of fine orography, and to be sure not to have
!  spurious cliffs near the coast (multiplication by xland during the
!  iterative process).
!
IF (KDXRATIO/=1 .OR. KDYRATIO/=1) THEN
!
!* allocations
! <<<<<<< spawn_zs.f90
   IXSIZE = IDIMX_C-2*(JPHEXT+1)
   IYSIZE = IDIMY_C-2*(JPHEXT+1)
   ALLOCATE(ZZS2_C(IXSIZE,IYSIZE))
   ALLOCATE(ZDZS_C(IXSIZE,IYSIZE))
   ALLOCATE(PZS1_C(IXSIZE,IYSIZE))
!!$=======
!!$  IXSIZE = KXEND-KXOR - 2*JPHEXT + 1
!!$  IYSIZE = KYEND-KYOR - 2*JPHEXT + 1
!!$>>>>>>> 1.1.4.1.18.2.2.1
!
!* constants
!
  ZRELAX=16./13.     ! best relaxation for infinite aspect ratio.
                     ! for dx=2, one should take 32./27. !!!
!
!* initializations of initial state
!
  JCOUNTER=0
  PZS1_C(:,:) = ZZS1_C(JPHEXT+2:IDIMX_C-JPHEXT-1,JPHEXT+2:IDIMY_C-JPHEXT-1)
  PZS2_C=0.
  CALL MPPDB_CHECK2D(PZS2_C,"SPAWN_ZSbefBKAT:PZS2",PRECISION)
!
!* iterative loop
!
  DO
!
!* interpolation
!
! <<<<<<< spawn_zs.f90
    CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4,  &
                   XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4,  &
                   2,2,IDIMX_C-1,IDIMY_C-1,KDXRATIO_C,KDYRATIO_C,1,  &
                   HLBCX,HLBCY,ZZS1_C,PZS2_C)
!!$=======
!!$      CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
!!$                     XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
!!$                     KXOR,KYOR,KXEND,KYEND,KDXRATIO,KDYRATIO,1,       &
!!$                     HLBCX,HLBCY,ZZS1,PZS2)
      CALL MPPDB_CHECK2D(PZS2_C,"SPAWN_ZS::PZS2_C/LOOP",PRECISION)
!!$>>>>>>> 1.1.4.1.18.2.2.1
    JCOUNTER=JCOUNTER+1
!
!* if orography is positive, it stays positive
!
! <<<<<<< spawn_zs.f90
     DO JI=1,IXSIZE
       DO JJ=1,IYSIZE
         IF (PZS1_C(JI,JJ)>-1.E-15) THEN
            JI2INF = (JI-1)*KDXRATIO_C+1+JPHEXT
            JI2SUP = JI*KDXRATIO_C+JPHEXT
            JJ2INF = (JJ-1)*KDYRATIO_C+1+JPHEXT
            JJ2SUP = JJ*KDYRATIO_C+JPHEXT
            PZS2_C(JI2INF:JI2SUP,JJ2INF:JJ2SUP)= MAX(PZS2_C(JI2INF:JI2SUP,JJ2INF:JJ2SUP),0.)
         ENDIF
       END DO
     END DO
!!$=======
!!$    DO JI=1,IXSIZE
!!$      DO JJ=1,IYSIZE
!!$        IF (PZS1(JI+KXOR,JJ+KYOR)>-1.E-15)                            &
!!$          PZS2((JI-1)*KDXRATIO+1+JPHEXT:JI*KDXRATIO+JPHEXT,             &
!!$               (JJ-1)*KDYRATIO+1+JPHEXT:JJ*KDYRATIO+JPHEXT)  =          &
!!$            MAX( PZS2((JI-1)*KDXRATIO+1+JPHEXT:JI*KDXRATIO+JPHEXT,      &
!!$                      (JJ-1)*KDYRATIO+1+JPHEXT:JJ*KDYRATIO+JPHEXT), 0.)
!!$      END DO
!!$    END DO
    CALL MPPDB_CHECK2D(PZS2_C,"SPAWN_ZS::PZS2_C/POS",PRECISION)
!!$>>>>>>> 1.1.4.1.18.2.2.1
!
!* computation of new averaged orography
! <<<<<<< spawn_zs.f90
   ZZS2_C(:,:) = 0.
   DO JI=1,IXSIZE
     DO JJ=1,IYSIZE
       DO JEPSX = (JI-1)*KDXRATIO_C+1+JPHEXT, JI*KDXRATIO_C+JPHEXT
         DO JEPSY = (JJ-1)*KDYRATIO_C+1+JPHEXT, JJ*KDYRATIO_C+JPHEXT
           ZZS2_C(JI,JJ) = ZZS2_C(JI,JJ) + PZS2_C(JEPSX,JEPSY)
         END DO
       END DO
     END DO
!!$=======
!!$>>>>>>> 1.1.4.1.18.2.2.1
    END DO
! <<<<<<< spawn_zs.f90
    ZZS2_C(:,:) = ZZS2_C(:,:) / (KDXRATIO_C*KDYRATIO_C)
    !
    ZDZS_C(:,:)=PZS1_C(:,:)-ZZS2_C(:,:)
!!$=======
!!$    ZZS2(:,:) = ZZS2(:,:) / (KDXRATIO*KDYRATIO)
!!$    !
!!$    ZDZS(:,:)=PZS1(KXOR+JPHEXT:KXEND-JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT)-ZZS2(:,:)
!!$>>>>>>> 1.1.4.1.18.2.2.1
!
!* test to end the iterative process
!
    ALLOCATE(ZDZS_3D(SIZE(ZDZS_C,1),SIZE(ZDZS_C,2),1))  ! WARNING : this is highly inefficient, this copy is unecessary
    ZDZS_3D(:,:,1)=ZDZS_C(:,:)                          ! We could write a function MAX2D_ll or use a POINTER for ZDZS_3D
    LOCMAXVAL=MAXVAL(ABS(ZDZS_C))
    CALL MPI_ALLREDUCE(LOCMAXVAL,ZMAXVAL,1,MPI_PRECISION,MPI_MAX,NMNH_COMM_WORLD,IINFO_ll)
    IF (ZMAXVAL<1.E-3) THEN
      EXIT
    ENDIF
!
    IF (JCOUNTER>=JMAXITER) THEN
      WRITE(ILUOUT,FMT=*) 'SPAWN_ZS: convergence of ',TRIM(HFIELD), &
                          ' NOT obtained after',JCOUNTER,' iterations'
      WRITE(ILUOUT,FMT=*) TRIM(HFIELD),                             &
         ' is modified to insure egality of large scale and averaged fine field'
! <<<<<<< spawn_zs.f90
      DO JI=1,IXSIZE
        DO JJ=1,IYSIZE
          DO JEPSX = (JI-1)*KDXRATIO_C+1+JPHEXT, JI*KDXRATIO_C+JPHEXT
            DO JEPSY = (JJ-1)*KDYRATIO_C+1+JPHEXT, JJ*KDYRATIO_C+JPHEXT
              PZS2_C(JEPSX,JEPSY) = PZS2_C(JEPSX,JEPSY) + ZDZS_C(JI,JJ)
!!$>>>>>>> 1.1.4.1.18.2.2.1
            END DO
          END DO
        END DO
      END DO
      !
      EXIT
    END IF
!
!* prints
!
    IF (NVERB >=7) THEN
      IZSMAX=MAXLOC(ABS(ZDZS_C(:,:)))
      IF (MOD(JCOUNTER,500)==1) THEN
        WRITE(ILUOUT,FMT='(A4,1X,A4,1X,A2,1X,A2,1X,A12,1X,A12,1X,A12)')      &
                          'n IT','nDIV','I1','J1','   ZS1','   ZS2','   DZS'
        WRITE(ILUOUT,FMT='(I4,1X,I4,1X,I2,1X,I2,1X,F12.7,1X,F12.7,1X,F12.7)')&
                            JCOUNTER,COUNT(ABS(ZDZS_C(:,:))>=1.E-3),          &
                            IZSMAX(1)+2,IZSMAX(2)+2,                          &
!                            PZS1(KXOR+IZSMAX(1),KYOR+IZSMAX(2)),              &
!                            ZTZS1_C(2+IZSMAX(1),2+IZSMAX(2)),                 &
                            ZZS2_C(IZSMAX(1),IZSMAX(2)),                      &
                            ZDZS_C(IZSMAX(1),IZSMAX(2))
      ENDIF
    END IF
!
!* correction of coarse orography
!
! <<<<<<< spawn_zs.f90
    ZZS1_C(JPHEXT+2:IDIMX_C-JPHEXT-1,JPHEXT+2:IDIMY_C-JPHEXT-1) = &
    ZZS1_C(JPHEXT+2:IDIMX_C-JPHEXT-1,JPHEXT+2:IDIMY_C-JPHEXT-1) + ZRELAX * ZDZS_C(:,:)

    ! update the Halo
    ! UPDATE_HALO_ll routines only work with fields of the size of the subdomain
    ! so we have to copy the values we want to update in a temporary field ZZS1CHILDGRID_C
    ALLOCATE(ZZS1CHILDGRID_C(SIZE(PZS2_C,1)+2,SIZE(PZS2_C,2)+2))
    ! TODO : renommer ZZS1CHILDGRID_C avec un nom plus explicite
    ZZS1CHILDGRID_C = 0.
    ! west boundary of ZZS1_C
    DO JI=1,JPHEXT+1
      DO JJ=1,IDIMY_C
        ZZS1CHILDGRID_C(JI,JJ) = ZZS1_C(JI,JJ)  ! distant value, not on local physical domain
        ZZS1CHILDGRID_C(JI+JPHEXT+1,JJ) = ZZS1_C(JI+JPHEXT+1,JJ) ! local value, on local physical domain
      END DO
    END DO
    ! east boundary of ZZS1_C
    DO JI=1,JPHEXT+1
      DO JJ=1,IDIMY_C
        ZZS1CHILDGRID_C(SIZE(PZS2_C,1)+2-JI+1,JJ) = ZZS1_C(IDIMX_C-JI+1,JJ)  ! distant value, not on local physical domain
        ZZS1CHILDGRID_C(SIZE(PZS2_C,1)+2-JI+1-JPHEXT-1,JJ) = ZZS1_C(IDIMX_C-JI+1-JPHEXT-1,JJ) ! local value, on local physical domain
      END DO
    END DO
    ! south boundary of ZZS1_C
    DO JI=1,IDIMX_C
      DO JJ=1,JPHEXT+1
        ZZS1CHILDGRID_C(JI,JJ) = ZZS1_C(JI,JJ)  ! distant value, not on local physical domain
        ZZS1CHILDGRID_C(JI,JJ+JPHEXT+1) = ZZS1_C(JI,JJ+JPHEXT+1) ! local value, on local physical domain
      END DO
    END DO
    ! north boundary of ZZS1_C
    DO JI=1,IDIMX_C
      DO JJ=1,JPHEXT+1
        ZZS1CHILDGRID_C(JI,SIZE(PZS2_C,2)+2-JJ+1) = ZZS1_C(JI,IDIMY_C-JJ+1)  ! distant value, not on local physical domain
        ZZS1CHILDGRID_C(JI,SIZE(PZS2_C,2)+2-JJ+1-JPHEXT-1) = ZZS1_C(JI,IDIMY_C-JJ+1-JPHEXT-1) ! local value, on local physical domain
      END DO
    END DO
    ! If we leave the north-east corner with zero values, UPDATE_HALO_EXTENDED_ll will
    ! cause errors on the south-east and north-west internal border of the neigbouring processes
    DO JI=1,JPHEXT+1
      DO JJ=1,JPHEXT+1
        ZZS1CHILDGRID_C(SIZE(PZS2_C,1)+2-JI+1-JPHEXT-1,SIZE(PZS2_C,2)+2-JJ+1-JPHEXT-1) &
         = ZZS1_C(IDIMX_C-JI+1-JPHEXT-1,IDIMY_C-JJ+1-JPHEXT-1) ! local value, on local physical domain
      END DO
    END DO
    !
    NULLIFY(TZZSFIELD_ll)
    CALL ADD2DFIELD_ll(TZZSFIELD_ll, ZZS1CHILDGRID_C)
    CALL UPDATE_HALO_EXTENDED_ll(TZZSFIELD_ll,IINFO)
    CALL CLEANLIST_ll(TZZSFIELD_ll)
    ! west and east boundaries - distant points
    DO JI=1,JPHEXT+1
      DO JJ=JPHEXT+1,IDIMY_C-JPHEXT+1
        ZZS1_C(JI,JJ) = ZZS1CHILDGRID_C(JI,JJ)
        ZZS1_C(IDIMX_C-JI+1,JJ) = ZZS1CHILDGRID_C(SIZE(PZS2_C,1)+2-JI+1,JJ)
      END DO
    END DO
    ! north and south boundaries - distant points
    DO JI=JPHEXT+1,IDIMX_C-JPHEXT+1
      DO JJ=1,JPHEXT+1
        ZZS1_C(JI,JJ) = ZZS1CHILDGRID_C(JI,JJ)
        ZZS1_C(JI,IDIMY_C-JJ+1) = ZZS1CHILDGRID_C(JI,SIZE(PZS2_C,2)+2-JJ+1)
      END DO
    END DO
    ! "corner" halo
    DO JI=1,JPHEXT+1
      DO JJ=1,JPHEXT+1
        ZZS1_C(JI,JJ) = ZZS1CHILDGRID_C(JI,JJ)
        ZZS1_C(IDIMX_C-JI+1,JJ) = ZZS1CHILDGRID_C(SIZE(PZS2_C,1)+2-JI+1,JJ)
        ZZS1_C(JI,IDIMY_C-JJ+1) = ZZS1CHILDGRID_C(JI,SIZE(PZS2_C,2)+2-JJ+1)
        ZZS1_C(IDIMX_C-JI+1,IDIMY_C-JJ+1) = ZZS1CHILDGRID_C(SIZE(PZS2_C,1)+2-JI+1,SIZE(PZS2_C,2)+2-JJ+1)
      END DO
    END DO
    ! corner points - distant points
    ! we have to treat the halo points in the corner separately to have correct values
    ! in the intersection of the halos (points (1,1), (1,2), (2,1), (2,2), (IDIMX_C,IDIMY_C), etc.)
!!$=======
!!$    ZZS1(KXOR+JPHEXT:KXEND-JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT) =                             &
!!$         ZZS1(KXOR+JPHEXT:KXEND-JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT) + ZRELAX * ZDZS(:,:)
!!$>>>>>>> 1.1.4.1.18.2.2.1
    !
    ! extrapolations (X direction)
    !
! <<<<<<< spawn_zs.f90
    ! TODO: utiliser JPHEXT dans une boucle pour generaliser au cas ou le halo est plus grand que 1
    IF(KXOR==1 .AND. KXEND==SIZE(PZS1_F,1) .AND. HLBCX(1)=='CYCL' ) THEN
      !c'est pris en compte et deja fait par UPDATE_HALO_ll et UPDATE_HALO2_ll ? --------> NON
!!$=======
!!$    IF(KXOR==1 .AND. KXEND==SIZE(PZS1,1) .AND. HLBCX(1)=='CYCL' ) THEN
!!$      ZZS1(KXOR,KYOR+JPHEXT:KYEND-JPHEXT)  = ZZS1(KXEND-JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT)
!!$      ZZS1(KXEND,KYOR+JPHEXT:KYEND-JPHEXT) = ZZS1(KXOR+JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT)
!!$>>>>>>> 1.1.4.1.18.2.2.1
    ELSE
! <<<<<<< spawn_zs.f90
      IF ( LWEST_ll() ) THEN
        ZZS1_C(1+JPHEXT,1:IDIMY_C) = 2. * ZZS1_C(2+JPHEXT,1:IDIMY_C)  - ZZS1_C(3+JPHEXT,1:IDIMY_C)
        ZZS1_C(JPHEXT,1:IDIMY_C)   = 2. * ZZS1_C(1+JPHEXT,1:IDIMY_C)  - ZZS1_C(2+JPHEXT,1:IDIMY_C)
      ENDIF
      IF ( LEAST_ll() ) THEN
          ZZS1_C(IDIMX_C-JPHEXT,1:IDIMY_C)   = 2. * ZZS1_C(IDIMX_C-JPHEXT-1,1:IDIMY_C) - ZZS1_C(IDIMX_C-JPHEXT-2,1:IDIMY_C)
          ZZS1_C(IDIMX_C-JPHEXT+1,1:IDIMY_C) = 2. * ZZS1_C(IDIMX_C-JPHEXT,1:IDIMY_C)   - ZZS1_C(IDIMX_C-JPHEXT-1,1:IDIMY_C)
      ENDIF
!!$=======
!!$      ZZS1(KXOR+JPHEXT-1,KYOR+JPHEXT:KYEND-JPHEXT) =                                       &
!!$        2. * ZZS1(KXOR+JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT)  - ZZS1(KXOR+JPHEXT+1,KYOR+JPHEXT:KYEND-JPHEXT)
!!$      IF(KXOR>1)                                                        &
!!$      ZZS1(KXOR+JPHEXT-2,KYOR+JPHEXT:KYEND-JPHEXT) =                                     &
!!$        2. * ZZS1(KXOR+JPHEXT-1,KYOR+JPHEXT:KYEND-JPHEXT)  - ZZS1(KXOR+JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT)
!!$      ZZS1(KXEND-JPHEXT+1,KYOR+JPHEXT:KYEND-JPHEXT) =                                      &
!!$        2. * ZZS1(KXEND-JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT) - ZZS1(KXEND-JPHEXT-1,KYOR+JPHEXT:KYEND-JPHEXT)
!!$      IF(KXEND<SIZE(PZS1,1))                                             &
!!$      ZZS1(KXEND-JPHEXT+2,KYOR+JPHEXT:KYEND-JPHEXT) =                                    &
!!$        2. * ZZS1(KXEND-JPHEXT+1,KYOR+JPHEXT:KYEND-JPHEXT) - ZZS1(KXEND-JPHEXT,KYOR+JPHEXT:KYEND-JPHEXT)
!!$>>>>>>> 1.1.4.1.18.2.2.1
    END IF
    !
    ! extrapolations (Y direction)
    !
! <<<<<<< spawn_zs.f90
!    IXMIN=MAX(KXOR-1,1)
!    IXMAX=MIN(KXEND+1,SIZE(PZS1_F,1))
    IF(KYOR==1 .AND. KYEND==SIZE(PZS1_F,2) .AND. HLBCY(1)=='CYCL' ) THEN
      !c'est pris en compte et deja fait par UPDATE_HALO_ll et UPDATE_HALO2_ll ? --------> NON
!!$=======
!!$    IXMIN=MAX(KXOR-1,1)
!!$    IXMAX=MIN(KXEND+1,SIZE(PZS1,1))
!!$    IF(KYOR==1 .AND. KYEND==SIZE(PZS1,2) .AND. HLBCY(1)=='CYCL' ) THEN
!!$      ZZS1(IXMIN:IXMAX,KYOR)  = ZZS1(IXMIN:IXMAX,KYEND-JPHEXT)
!!$      ZZS1(IXMIN:IXMAX,KYEND) = ZZS1(IXMIN:IXMAX,KYOR+JPHEXT)
!!$>>>>>>> 1.1.4.1.18.2.2.1
    ELSE
! <<<<<<< spawn_zs.f90
      IF ( LSOUTH_ll() ) THEN
        ZZS1_C(1:IDIMX_C,1+JPHEXT) = 2. * ZZS1_C(1:IDIMX_C,2+JPHEXT)  - ZZS1_C(1:IDIMX_C,3+JPHEXT)
        ZZS1_C(1:IDIMX_C,JPHEXT)   = 2. * ZZS1_C(1:IDIMX_C,1+JPHEXT)  - ZZS1_C(1:IDIMX_C,2+JPHEXT)
      ENDIF
      IF ( LNORTH_ll() ) THEN
        ZZS1_C(1:IDIMX_C,IDIMY_C-JPHEXT)   = 2. * ZZS1_C(1:IDIMX_C,IDIMY_C-JPHEXT-1) - ZZS1_C(1:IDIMX_C,IDIMY_C-JPHEXT-2)
        ZZS1_C(1:IDIMX_C,IDIMY_C-JPHEXT+1) = 2. * ZZS1_C(1:IDIMX_C,IDIMY_C-JPHEXT)   - ZZS1_C(1:IDIMX_C,IDIMY_C-JPHEXT-1)
      ENDIF
!!$=======
!!$      ZZS1(IXMIN:IXMAX,KYOR+JPHEXT-1) =                                       &
!!$        2. * ZZS1(IXMIN:IXMAX,KYOR+JPHEXT)  - ZZS1(IXMIN:IXMAX,KYOR+JPHEXT+1)
!!$      IF(KYOR>1)                                                     &
!!$      ZZS1(IXMIN:IXMAX,KYOR+JPHEXT-2) =                                     &
!!$        2. * ZZS1(IXMIN:IXMAX,KYOR+JPHEXT-1)    - ZZS1(IXMIN:IXMAX,KYOR+JPHEXT)
!!$      ZZS1(IXMIN:IXMAX,KYEND-JPHEXT+1) =                                      &
!!$        2. * ZZS1(IXMIN:IXMAX,KYEND-JPHEXT) - ZZS1(IXMIN:IXMAX,KYEND-JPHEXT-1)
!!$      IF(KYEND<SIZE(PZS1,2))                                          &
!!$      ZZS1(IXMIN:IXMAX,KYEND-JPHEXT+2) =                                    &
!!$        2. * ZZS1(IXMIN:IXMAX,KYEND-JPHEXT+1)   - ZZS1(IXMIN:IXMAX,KYEND-JPHEXT)
!!$>>>>>>> 1.1.4.1.18.2.2.1
    END IF
!
    DEALLOCATE(ZZS1CHILDGRID_C)
    DEALLOCATE(ZDZS_3D)
!
  END DO
!
  CALL ZS_BOUNDARY(PZS2_C,ZZS2_LS)
  JI=0
  CALL MPPDB_CHECK2D(PZS2_C,"SPAWN_ZSend:PZS2",PRECISION)
!
  WRITE(ILUOUT,FMT=*) 'convergence of ',TRIM(HFIELD),' obtained after ', &
                      JCOUNTER,' iterations'
!
  DEALLOCATE(ZZS2_C)
  DEALLOCATE(ZDZS_C)
  DEALLOCATE(ZZS1_C)
END IF
!
IF (PRESENT(PZS2_LS)) PZS2_LS(:,:)=ZZS2_LS(:,:)
DEALLOCATE(ZZS2_LS)
!
CALL GOTO_MODEL(IMI)
CALL GO_TOMODEL_ll(IMI,IINFO_ll)
!-------------------------------------------------------------------------------
END SUBROUTINE SPAWN_ZS
!
