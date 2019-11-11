!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ############################## 
      MODULE MODI_LS_COUPLING
!     ##############################
!
INTERFACE
!
      SUBROUTINE LS_COUPLING (PTSTEP, OSTEADY_DMASS, HCONF,                  &
            HGETTKEM,                                                        &
            HGETRVM,HGETRCM,HGETRRM,HGETRIM,                                 &
            HGETRSM,HGETRGM,HGETRHM,HGETSVM,OCH_INIT_FIELD,KSV,              &
            KIMAX_ll,KJMAX_ll,                                               &
            KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,               &
            KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
            KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &
            PLSUM,PLSVM,PLSWM, PLSTHM,PLSRVM,PLSZWSM,PDRYMASST,              &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,            &
            PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS,PLSZWSS,PDRYMASSS,               &
            PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,            &
            PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS             )
!
REAL,             INTENT(IN)       :: PTSTEP      ! Time step
LOGICAL,          INTENT(IN)       :: OSTEADY_DMASS ! Md evolution logical switch
!
CHARACTER(LEN=*), INTENT(IN) :: HCONF  ! configuration var. linked to FMfile
                                                  ! Get indicators at t-dt for
CHARACTER(LEN=*), INTENT(IN) :: HGETTKEM          !  tke
CHARACTER(LEN=*), INTENT(IN) :: HGETRVM,HGETRCM,HGETRRM   !  vapor, cloud, rain
CHARACTER(LEN=*), INTENT(IN) :: HGETRIM,HGETRSM   !  ice, snow
CHARACTER(LEN=*), INTENT(IN) :: HGETRGM,HGETRHM   !  graupel, hail
CHARACTER(LEN=*),                                 &
    DIMENSION(:), INTENT(IN) :: HGETSVM           !  scalar variables
LOGICAL,          INTENT(IN) :: OCH_INIT_FIELD
INTEGER,                INTENT(IN)        :: KSV  ! number of scalar var.
INTEGER,               INTENT(IN)   :: KIMAX_ll  !  Dimensions  in x direction 
                                                 ! of the physical domain,
INTEGER,               INTENT(IN)   :: KJMAX_ll  !  Dimensions  in y direction 
                                                 ! of the physical domain,
! sizes of the West-east total LB area
INTEGER, INTENT(IN) :: KSIZELBX_ll,KSIZELBXU_ll      ! for T,V,W and u 
INTEGER, INTENT(IN) :: KSIZELBXTKE_ll                ! for TKE
INTEGER, INTENT(IN) :: KSIZELBXR_ll,KSIZELBXSV_ll    ! for Rx and SV    
! sizes of the North-south total LB area
INTEGER, INTENT(IN) :: KSIZELBY_ll,KSIZELBYV_ll      ! for T,U,W  and v
INTEGER, INTENT(IN):: KSIZELBYTKE_ll                ! for TKE
INTEGER, INTENT(IN) :: KSIZELBYR_ll,KSIZELBYSV_ll    ! for Rx and SV 
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLSUM,PLSVM,PLSWM ! Large Scale 
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLSTHM, PLSRVM    ! fields at t-dt
REAL, DIMENSION(:,:),    INTENT(IN)  :: PLSZWSM              ! fields at t-dt
REAL,                    INTENT(IN)  :: PDRYMASST         ! Mass of dry air Md
! larger scale fields for Lateral Boundary condition
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTKEM
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
!
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLSUS,PLSVS,PLSWS  ! Large Scale 
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLSTHS,PLSRVS      ! source terms
REAL, DIMENSION(:,:),    INTENT(OUT) :: PLSZWSS              ! source terms
REAL,                    INTENT(OUT) :: PDRYMASSS          !  Md source 
! larger scale fields sources for Lateral Boundary condition
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXUS,PLBXVS,PLBXWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYUS,PLBYVS,PLBYWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTKES          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTKES
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBXRS  ,PLBXSVS  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBYRS  ,PLBYSVS  ! in x and y-dir.
!
!
END SUBROUTINE LS_COUPLING
!
END INTERFACE
!
END MODULE MODI_LS_COUPLING
!
!
!
!
!     ######################################################################
      SUBROUTINE LS_COUPLING (PTSTEP, OSTEADY_DMASS, HCONF,                  &
            HGETTKEM,                                                        &
            HGETRVM,HGETRCM,HGETRRM,HGETRIM,                                 &
            HGETRSM,HGETRGM,HGETRHM,HGETSVM,OCH_INIT_FIELD,                  &
            KSV,KIMAX_ll,KJMAX_ll,                                           &
            KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,               &
            KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
            KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &
            PLSUM,PLSVM,PLSWM, PLSTHM,PLSRVM,PLSZWSM,PDRYMASST,              &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,            &
            PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS,PLSZWSS,PDRYMASSS,               &
            PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,            &
            PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS             )
!     ######################################################################
!
!!****  *LS_COUPLING* - Refreshing of model 1 Large Scale sources 
!!****                         of scalar fields.
!!
!!    PURPOSE
!!    -------
!!      The purpose of LS_COUPLING is to 'refresh' Large Scale sources 
!!    of scalar fields when the current time step reachs a coupling one.
!!    It only concerns model 1.
!
!
!!**  METHOD
!!    ------
!!      If the next time step reachs a coupling one, the Large Scale sources
!!    are 'refreshed' :
!!      - grid dimensions of initial file and coupling file are checked;
!!      - large scale scalar fields are read in the coupling file, which is closed;
!!      - large scale sources are computed from the current large scale fields
!!    and those just read (same computation as in subroutine INI_CPL) .
!!
!!    EXTERNAL
!!    --------
!!      IO_READ_FIELD : to read data in file
!!      IO_FILE_CLOSE_ll : to close a file
!!      INI_LS      : to initialize larger scale fields
!!      INI_LB      : to initialize "2D" surfacic LB fields 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS 
!!         JPHEXT : Horizontal EXTernal points number
!!         JPVEXT : Vertical EXTernal points number
!!
!!      Module MODD_CTURB :
!!         XTKEMIN : mimimum value for the TKE
!!
!!      Module MODD_DYN   
!!         NCPL_NBR    : NumBeR of CouPLing files
!!         NCPL_CUR    : Number of the CURrent CouPLing file              
!!         NCPL_TIMES  : List of the values of the time index in the temporal
!!                      model loop where large scale sources are computed
!!
!!      Module MODD_LUNIT1 
!!         TCPLFILE    :  CouPLing FILEs
!!
!!    REFERENCE
!!    ---------
!!    
!!
!!    AUTHOR
!!    ------
!!    J. P. Lafore  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     03/09/96   The previous routine SET_COUPLING have been splitted
!!                             in 2 routines (UVW_LS_COUPLING and LS_COUPLING),
!!                             and the temporal advance have been removed.
!!                              Correction of the LS sources names (removing of R).
!!                   25/09/96   (V. Masson) test for coupling made in MODEL1
!!                   11/03/97   (Lafore)   "surfacic" LS fieds  introduction
!!                   10/04/97   (Lafore)   proper treatment of minima for LS-fields
!!                   30/10/97   (Stein)    use the Lhorelax information for LB
!!                   20/04/98   (Josse)    temporal evolution of SST
!!                   22/09/98   (Ducrocq) //,  and introduce INI_LS and INI_LB
!!                   01/2004    (Masson) removes SST forcing (surface externalization)
!!                   05/2006    Remove KEPS
!!                   2/2014     (escobar) add paspol for Forefire
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!! 
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
USE MODD_PARAMETERS
USE MODD_CTURB
USE MODD_DYN
USE MODD_LUNIT_n
USE MODD_NSV
#ifdef MNH_FOREFIRE
USE MODD_PASPOL
#endif
USE MODD_CH_MNHC_n
!
USE MODE_FM
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_MSG
!
USE MODI_INI_LS
USE MODI_INI_LB
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
REAL,             INTENT(IN)       :: PTSTEP      ! Time step
LOGICAL,          INTENT(IN)       :: OSTEADY_DMASS ! Md evolution logical switch
CHARACTER(LEN=*), INTENT(IN) :: HCONF  ! configuration var. linked to FMfile
!
                                                  ! Get indicators at t-dt for
CHARACTER(LEN=*), INTENT(IN) :: HGETTKEM          !  tke, eps
CHARACTER(LEN=*), INTENT(IN) :: HGETRVM,HGETRCM,HGETRRM   !  vapor, cloud, rain
CHARACTER(LEN=*), INTENT(IN) :: HGETRIM,HGETRSM   !  ice, snow
CHARACTER(LEN=*), INTENT(IN) :: HGETRGM,HGETRHM   !  graupel, hail
CHARACTER(LEN=*),                                 &
    DIMENSION(:), INTENT(IN) :: HGETSVM           !  scalar variables
LOGICAL,          INTENT(IN) :: OCH_INIT_FIELD
INTEGER,                INTENT(IN)        :: KSV  ! number of scalar var.
INTEGER,               INTENT(IN)   :: KIMAX_ll  !  Dimensions  in x direction 
                                                 ! of the physical domain,
INTEGER,               INTENT(IN)   :: KJMAX_ll  !  Dimensions  in y direction 
                                                 ! of the physical domain,
! sizes of the West-east total LB area
INTEGER, INTENT(IN) :: KSIZELBX_ll,KSIZELBXU_ll      ! for T,V,W and u 
INTEGER, INTENT(IN) :: KSIZELBXTKE_ll                ! for TKE 
INTEGER, INTENT(IN) :: KSIZELBXR_ll,KSIZELBXSV_ll    ! for Rx and SV    
! sizes of the North-south total LB area
INTEGER, INTENT(IN) :: KSIZELBY_ll,KSIZELBYV_ll      ! for T,U,W  and v
INTEGER, INTENT(IN):: KSIZELBYTKE_ll                ! for TKE 
INTEGER, INTENT(IN) :: KSIZELBYR_ll,KSIZELBYSV_ll    ! for Rx and SV 
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLSUM,PLSVM,PLSWM ! Large Scale 
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLSTHM, PLSRVM    ! fields at t-dt
REAL, DIMENSION(:,:),    INTENT(IN)  :: PLSZWSM
REAL,                    INTENT(IN)  :: PDRYMASST         ! Mass of dry air Md
! larger scale fields for Lateral Boundary condition
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PLBYTKEM                  
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(IN) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
!
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLSUS,PLSVS,PLSWS  ! Large Scale 
REAL, DIMENSION(:,:,:),  INTENT(OUT) :: PLSTHS,PLSRVS      ! source terms
REAL,                    INTENT(OUT) :: PDRYMASSS          !  Md source 
! larger scale fields sources for Lateral Boundary condition
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXUS,PLBXVS,PLBXWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYUS,PLBYVS,PLBYWS ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTHS              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTKES          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTKES
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBXRS  ,PLBXSVS  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBYRS  ,PLBYSVS  ! in x and y-dir.
REAL, DIMENSION(:,:),    INTENT(OUT) :: PLSZWSS              ! source terms
!
!
!*       0.2   declarations of local variables
!
INTEGER                :: ILUOUT                     !  Logical unit number
                                                     ! associated with TLUOUT 
INTEGER                :: IRESP
REAL                   :: ZLENG                      ! Interpolation length
INTEGER                :: IIMAX,IJMAX,IKMAX       !  Dimensions  of the physical 
                                                  ! part of the arrays stored in
                                                  ! coupling file
INTEGER                :: IKU             !  Dimensions of arrays in 
                                                  ! initial file
LOGICAL               :: GLSOURCE ! switch for the source term (for ini_ls and ini_lb)
!CHARACTER(LEN=4), DIMENSION(KSV)    :: YGETSVM
!
!-------------------------------------------------------------------------------
!
!*      1.   COMPUTE LARGE SCALE SOURCES
!            ---------------------------
!
!*      1.1  Check dimensions
!
CALL IO_READ_FIELD(TCPLFILE(NCPL_CUR)%TZFILE,'IMAX',IIMAX)
CALL IO_READ_FIELD(TCPLFILE(NCPL_CUR)%TZFILE,'JMAX',IJMAX)
CALL IO_READ_FIELD(TCPLFILE(NCPL_CUR)%TZFILE,'KMAX',IKMAX)
!
IKU=SIZE(PLSTHM,3)
!
IF ( (IIMAX/=KIMAX_ll) .OR. (IJMAX/=KJMAX_ll)                    &
                             .OR. (IKMAX/=(IKU-2*JPVEXT)) ) THEN
  ILUOUT = TLUOUT%NLU
  WRITE(ILUOUT,FMT=9003)
  WRITE(ILUOUT,*) 'THE GRIDS ARE DIFFERENT IN THE INITIAL FILE :'
  WRITE(ILUOUT,*) KIMAX_ll,'*',KJMAX_ll,'*',IKU-2*JPVEXT
  WRITE(ILUOUT,*) 'AND IN THE COUPLING FILE :',TCPLFILE(NCPL_CUR)%TZFILE%CNAME
  WRITE(ILUOUT,*) IIMAX,'*',IJMAX,'*',IKMAX
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','LS_COUPLING','')
END IF
!
!*      1.2  Read  new scalar Large Scale fields and update sources
!
GLSOURCE=.TRUE.
ZLENG = (NCPL_TIMES(NCPL_CUR,1) - NCPL_TIMES(NCPL_CUR-1,1)) * PTSTEP 
!
CALL INI_LS(TCPLFILE(NCPL_CUR)%TZFILE,HGETRVM,GLSOURCE,PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS, &
            PLSZWSS, PDRYMASSS,PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,PDRYMASST,ZLENG,OSTEADY_DMASS)

!
!
!*      1.3   Read  new  LB sources and update sources
!
GLSOURCE=.TRUE.
!
!YGETSVM(1:KSV) = HGETSVM(1:KSV)
!IF ( LUSECHEM .AND. (.NOT. OCH_INIT_FIELD) ) &
!    YGETSVM(NSV_CHEMBEG: NSV_CHEMEND) = 'INIT'
!IF (HCONF == 'RESTA')  THEN
#ifdef MNH_FOREFIRE
!   IF (LPASPOL) YGETSVM(NSV_PPBEG: NSV_PPEND) = 'INIT'
#endif
!   IF (NSV_USER /= 0) YGETSVM(1/NSV_USER) = 'INIT'
!   IF (NSV_C2R2 /= 0) YGETSVM(NSV_C2R2BEG: NSV_C2R2END) = 'INIT'
!   IF (NSV_C1R3 /= 0) YGETSVM(NSV_C1R3BEG: NSV_C1R3END) = 'INIT'
!   IF (NSV_ELEC /= 0) YGETSVM(NSV_ELECBEG: NSV_ELECEND) = 'INIT'
!   IF (NSV_LG   /= 0) YGETSVM(NSV_LGBEG: NSV_LGEND) = 'INIT'
!   IF (NSV_LNOX /= 0) YGETSVM(NSV_LNOXBEG: NSV_LNOXEND) = 'INIT'
!   IF (NSV_DST  /= 0) YGETSVM(NSV_DSTBEG: NSV_DSTEND) = 'INIT'
!   IF (NSV_SLT  /= 0) YGETSVM(NSV_SLTBEG: NSV_SLTEND) = 'INIT'
!   IF (NSV_DSTDEP /= 0) YGETSVM(NSV_DSTDEPBEG: NSV_DSTDEPEND) = 'INIT'
!   IF (NSV_SLTDEP /= 0) YGETSVM(NSV_SLTDEPBEG: NSV_SLTDEPEND) = 'INIT'
!   IF (NSV_PP   /= 0) YGETSVM(NSV_PPBEG: NSV_PPEND) = 'INIT'
!   IF (NSV_CS   /= 0) YGETSVM(NSV_CSBEG: NSV_CSEND) = 'INIT'
#ifdef MNH_FOREFIRE
!   IF (NSV_FF   /= 0) YGETSVM(NSV_FFBEG: NSV_FFEND) = 'INIT'
#endif
!END IF
!
CALL INI_LB(TCPLFILE(NCPL_CUR)%TZFILE,GLSOURCE,KSV,                   &
     KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,               &
     KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
     KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &
     HGETTKEM,HGETRVM,HGETRCM,HGETRRM,HGETRIM,HGETRSM,                &
     HGETRGM,HGETRHM,HGETSVM,                                         &
     PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,            &
     PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS,            &
     PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
     PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,            &
     ZLENG)
!
!
!*      1.4  Close the coupling file
!
CALL IO_FILE_CLOSE_ll(TCPLFILE(NCPL_CUR)%TZFILE)
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!
!*      2.    FORMATS
!             -------
!
9003  FORMAT(/,'FATAL ERROR IN LS_COUPLING FOR MODEL ', I2,' : ',/, &
             '--------------------------------------' )
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LS_COUPLING
