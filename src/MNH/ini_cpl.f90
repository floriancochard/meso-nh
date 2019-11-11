!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_INI_CPL
!     ###################
!
INTERFACE
!
      SUBROUTINE INI_CPL(KSTOP,PTSTEP,OSTEADY_DMASS,HCONF,                   &
            HGETTKEM,                                                        &
            HGETRVM,HGETRCM,HGETRRM,HGETRIM,                                 &
            HGETRSM,HGETRGM,HGETRHM,HGETSVM,OCH_INIT_FIELD,                  &
            KSV,KIMAX_ll,KJMAX_ll,                                           &
            KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,               &
            KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
            KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &
            PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,PDRYMASST,               &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,            &
            PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS,PLSZWSS,PDRYMASSS,               &
            PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,            &
            PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS             )
!
INTEGER,          INTENT(IN)       :: KSTOP       ! Number of time steps for
                                                  ! current segment
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
REAL, DIMENSION(:,:),    INTENT(IN)  :: PLSZWSM           ! Significant wave height
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PLSTHM, PLSRVM    ! fields at t-dt
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
REAL, DIMENSION(:,:),    INTENT(OUT) :: PLSZWSS            ! Significant wave height
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
END SUBROUTINE INI_CPL
!
END INTERFACE
!
END MODULE MODI_INI_CPL
!
!
!     #####################################################################
      SUBROUTINE INI_CPL(KSTOP,PTSTEP,OSTEADY_DMASS,HCONF,                   &
            HGETTKEM,                                                        &
            HGETRVM,HGETRCM,HGETRRM,HGETRIM,                                 &
            HGETRSM,HGETRGM,HGETRHM,HGETSVM,OCH_INIT_FIELD,                  &
            KSV,KIMAX_ll,KJMAX_ll,                                           &
             KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,              &
            KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                   &
            KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,           &   
            PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,PDRYMASST,               &
            PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,            &
            PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,            &
            PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS,PLSZWSS,PDRYMASSS,               &
            PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,            &
            PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS             )
!     #####################################################################
!
!!****  *INI_CPL * - routine to initialize variables for the coupling
!!                   of the model 1 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize specific variables for the
!     coupling of the model 1, from current time and large scale fields read in
!     the coupling files.
!      
!!
!!**  METHOD
!!    ------
!!      First, all the coupling files are opened and current times are read to
!!    initialize the coupling time array. 
!!      Then, several tests are performed :
!!      - since the calendar is not taken into account yet, dates of the coupling
!!    times and the current date (which corresponds of this segment beginning 
!!    date) must be the same, otherwise the initialization is stopped;
!!      - if the coupling files are not specified in a chronological order, the
!!    initialization is stopped;
!!      - the coupling files, for which current time occurs before the current
!!    time (which corresponds of this segment beginning time), are released;
!!      - the number of the current coupling file is initialized as the first
!!    one following in time, the initial file;
!!      - the actual number of coupling files are reinitialized and the coupling
!!    times are converted in number of time step;
!!      - if the segment end is greater than the last coupling time, the 
!!    initialization is stopped;
!!      - finally, names of the useful coupling files are written in the 
!!    output-listing.
!!      Secondly, after checking grid dimensions, large scale fields are read in
!!    the first coupling file and source terms are computed. Note that this 
!!    computation occurs also in subroutine SET_COUPLING.
!!      The same treatment is also performed to couple the DRYMASS.
!!         
!!
!!    EXTERNAL
!!    --------
!!      IO_READ_FIELD: to read data in LFI_FM file
!!      IO_FILE_CLOSE_ll : to close a FM-file
!!      INI_LS      : to initialize larger scale fields
!!      INI_LB      : to initialize "2D" surfacic LB fields 
!!      DATETIME_DISTANCE : compute the temporal distance in seconds between 2 dates
!!
!!      Module MODE_TIME : contains SM_PRINT_TIME routine
!!                         and uses module MODD_TIME (for definition
!!                         of types DATE_TIME and DATE)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS 
!!         JPHEXT : Horizontal external points number
!!         JPVEXT : Vertical external points number
!!         JPCPLFILEMAX : Maximum allowed number of coupling files
!!
!!      Module MODD_CONF   
!!         NVERB      : Level of informations on output-listing
!!
!!      Module MODD_CTURB :
!!         XTKEMIN : mimimum value for the TKE
!!
!!      Module MODD_DYN   
!!         XSEGLEN     : Duration of segment (in seconds)
!!
!!         NCPL_NBR    : NumBeR of CouPLing files
!!         NCPL_CUR    : Number of the CURrent CouPLing file              
!!         NCPL_TIMES  : List of the values of the time index in the temporal
!!                      model loop where large scale sources are computed
!!
!!      Module MODD_NESTING : NDTRATIO  
!!
!!      Module MODD_TIME  
!!         TDTCPL       : Time and Date of the CouPLing files
!!
!!      Module MODD_TIME1 
!!         TDTCUR      : CURrent Time and Date (in the initial file)
!!
!!      Module MODD_LUNIT1 
!!         CCPLFILE    :  Names of the CouPLing FILEs 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_CPL)
!!      
!!
!!    AUTHOR
!!    ------
!!	I.Mallet       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                10/03/95 
!!      Modifications  (Mallet)  14/09/95  to take into account different days  
!!      Modifications  (Lafore)  20/09/95  Introduction of the DRYMASS coupling
!!                     (Stein)   02/01/96  add the calendar 
!!                     (Lafore)  11/03/97  "surfacic" LS fieds  introduction
!!      Modification   (Lafore)  10/04/97   proper treatment of minima for LS-fields
!!      Modification   (Stein)   22/12/97  new LS fields  
!!      Modification   (Josse)   20/04/98  temporal evolution of SST
!!      Modification   (Ducrocq)  22/09/98  //,  and introduce INI_LS and INI_LB
!!      Modification   (Masson)  01/2004 surface externalization, removes 
!!                                       SST forcing
!!      Modification             05/2006  Remove KEPS
!!                     (Escobar) 2/2014 add Forefire
!!                     (J.Escobar) 26/03/2014 bug in init of NSV_USER on RESTA case
!!                     (P.Wautelet)28/03/2018 replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CH_MNHC_n
USE MODD_CONF
USE MODD_CTURB
USE MODD_DYN
USE MODD_LUNIT_n, ONLY: CCPLFILE, TCPLFILE, TLUOUT
USE MODD_NSV
USE MODD_PARAMETERS
USE MODD_TIME_n
! #ifdef MNH_FOREFIRE
! USE MODD_PASPOL
! USE MODD_FOREFIRE
! USE MODD_FOREFIRE_n
! #endif
!
USE MODE_DATETIME
USE MODE_FM, ONLY: IO_FILE_OPEN_ll, IO_FILE_CLOSE_ll
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY : IO_FILE_ADD2LIST
USE MODE_MSG
USE MODD_NESTING
USE MODE_TIME
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
INTEGER,          INTENT(IN)       :: KSTOP       ! Number of time steps for
                                                  ! current segment
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
REAL, DIMENSION(:,:),    INTENT(IN)  :: PLSZWSM           ! Significant wave height at t-dt
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
REAL, DIMENSION(:,:),    INTENT(OUT) :: PLSZWSS           ! Significant wave height
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
!*       0.2   declarations of local variables
!
INTEGER                :: ILUOUT                     !  Logical unit number
INTEGER                :: IRESP
CHARACTER (LEN=40)     :: YTITLE                     !  Title for date print 
INTEGER                :: JCI                        !  Loop index on number of
                                                     ! coupling files
CHARACTER (LEN=2)      :: YCI                        !  String for coupling files
                                                     ! index
REAL                   :: ZLENG                      ! Interpolation length
LOGICAL                :: GEND                    !  Logical to see if coupling
                                                  ! times respect the chronolo.
                                                  ! order or the segment length
INTEGER                :: IIMAX,IJMAX,IKMAX       !  Dimensions  of the physical 
                                                  ! part of the arrays stored in
                                                  ! coupling file
INTEGER                :: IKU             !  Dimensions of arrays in 
                                                  ! initial file

INTEGER                :: ICPLEND                 ! number of the last cpl file
LOGICAL, DIMENSION(JPCPLFILEMAX)    :: GSKIP      ! array to skip or not after
                                                  ! a cpl file
REAL                   :: ZDIST                   ! temporal distance in secunds
                                                  ! between 2 dates
LOGICAL :: GLSOURCE ! switch for the source term (for ini_ls and ini_lb)
!CHARACTER(LEN=4), DIMENSION(KSV)    :: YGETSVM
!
!-------------------------------------------------------------------------------
!
!*       1.    SOME INITIALIZATIONS
!              --------------------
!
GSKIP(:)=.FALSE.
GEND=.FALSE.
NCPL_TIMES(:,:) = NUNDEF
ILUOUT = TLUOUT%NLU
!
!-------------------------------------------------------------------------------
!
!*       2.    CHECK COUPLING FILES DATES
!              --------------------------
!
DO JCI=1,NCPL_NBR
  WRITE(YCI,'(I2.0)') JCI
  CALL IO_FILE_ADD2LIST(TCPLFILE(JCI)%TZFILE,CCPLFILE(JCI),'UNKNOWN','READ',KLFITYPE=2,KLFIVERB=NVERB)
  CALL IO_FILE_OPEN_ll(TCPLFILE(JCI)%TZFILE,KRESP=IRESP)
  IF (IRESP /= 0) THEN
    CALL PRINT_MSG(NVERB_FATAL,'IO','INI_CPL','problem when opening coupling file '//TRIM(YCI))
  END IF
!
!*       2.1   Read current time in coupling files
!
  CALL IO_READ_FIELD(TCPLFILE(JCI)%TZFILE,'DTCUR',TDTCPL(JCI))
!
!*       2.2   Check chronological order
!
  CALL DATETIME_DISTANCE(TDTCUR,TDTCPL(JCI),ZDIST)
  !
  IF ( ZDIST <=0. ) THEN
    WRITE(ILUOUT,FMT=9002) 1
    WRITE(ILUOUT,*) 'YOUR COUPLING FILE ',JCI,' IS PREVIOUS TO THE DATE &
               & CORRESPONDING TO THE BEGINNING OF THE SEGMENT. IT WILL &
               & NOT BE TAKEN INTO ACCOUNT.'
    YTITLE='CURRENT DATE AND TIME IN THE INITIAL FILE'
    CALL SM_PRINT_TIME(TDTCUR,TLUOUT,YTITLE)
    YTITLE='CURRENT DATE AND TIME OF THE FILE'//YCI
    CALL SM_PRINT_TIME(TDTCPL(JCI),TLUOUT,YTITLE)
    GSKIP(JCI)=.TRUE.      ! flag to skip after this coupling file
  ELSE
    NCPL_TIMES(JCI,1) = NINT( ZDIST / PTSTEP ) + 2
  END IF
  !
  IF (JCI > 1) THEN
    CALL DATETIME_DISTANCE(TDTCPL(JCI-1),TDTCPL(JCI),ZDIST)
    !
    IF ( ZDIST < 0. ) THEN
      WRITE(ILUOUT,FMT=9003) 1
      WRITE(ILUOUT,*) 'YOU MUST SPECIFY THE COUPLING FILES IN A CHRONOLOGICAL &
                       & ORDER'
      YTITLE='CURRENT DATE AND TIME OF THE FILE'//YCI
      CALL SM_PRINT_TIME(TDTCPL(JCI),TLUOUT,YTITLE)
      WRITE(YCI,'(I2.0)') JCI-1        
      YTITLE='CURRENT DATE AND TIME OF THE FILE'//YCI
      CALL SM_PRINT_TIME(TDTCPL(JCI-1),TLUOUT,YTITLE)
      GEND=.TRUE.           ! error flag set to true   
    END IF
  !   
  END IF
END DO
!  exit when a fatal error has been encountered
IF ( GEND ) THEN
  RETURN
END IF
!
!*       2.3   Find the current coupling file and the last one
!
NCPL_CUR=1
DO JCI=1,NCPL_NBR
  IF( GSKIP(JCI) ) THEN
    NCPL_CUR = NCPL_CUR +1
  END IF 
  ! 
END DO
!
ICPLEND=NCPL_NBR
DO JCI=NCPL_CUR,NCPL_NBR
  IF( NCPL_TIMES(JCI,1) > KSTOP ) THEN
    ICPLEND = JCI 
    EXIT
  END IF
END DO
!
IF (ICPLEND==NCPL_NBR .AND. NCPL_TIMES(NCPL_NBR,1) < KSTOP) THEN
  WRITE(ILUOUT,FMT=9003) 1
  WRITE(ILUOUT,*) 'THE SEGMENT END IS GREATER THAN THE LAST COUPLING TIME'
  WRITE(ILUOUT,*) 'SPECIFY ANOTHER COUPLING FILE OR DECREASE THE SEGMENT LENGTH.'
  WRITE(ILUOUT,*) 'PLEASE, REFER TO THE USER GUIDE TO OBTAIN MORE INFORMATIONS'
  WRITE(ILUOUT,*) 'ON THE TEMPORAL GRID.'
  YTITLE='CURRENT DATE AND TIME OF THE LAST FILE'
  CALL SM_PRINT_TIME(TDTCPL(NCPL_NBR),TLUOUT,YTITLE)
  YTITLE='DATE AND TIME OF THE BEGINNING OF THE SEGMENT YOU WANT TO BE PERFORMED'
  CALL SM_PRINT_TIME(TDTCUR,TLUOUT,YTITLE)
  WRITE(ILUOUT,*) 'XSEGLEN = ', XSEGLEN 
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_CPL','')
END IF
! save the right number of useful coupling files
NCPL_NBR = ICPLEND - NCPL_CUR + 1
!         
!
!*      2.4   Write in output-listing the useful coupling files
!
!
WRITE(ILUOUT,*) 'FOR THE EVOLUTION OF LARGE SCALE FIELDS THE SEGMENT REQUIRES ',&
                & NCPL_NBR,' COUPLING FILES :'
DO JCI=NCPL_CUR,NCPL_CUR+NCPL_NBR-1
  WRITE(ILUOUT,*) TCPLFILE(JCI)%TZFILE%CNAME,' UNTIL TIME STEP : ',NCPL_TIMES(JCI,1)
END DO
!-------------------------------------------------------------------------------
!
!*      3.    FIRST COMPUTATION OF LARGE SCALE SOURCES
!             ----------------------------------------
!
!*      3.1   Read dimensions in coupling file and checks with initial file
!
CALL IO_READ_FIELD(TCPLFILE(NCPL_CUR)%TZFILE,'IMAX',IIMAX)
CALL IO_READ_FIELD(TCPLFILE(NCPL_CUR)%TZFILE,'JMAX',IJMAX)
CALL IO_READ_FIELD(TCPLFILE(NCPL_CUR)%TZFILE,'KMAX',IKMAX)
!
IKU=SIZE(PLSUM,3)
!
IF ( (IIMAX/=KIMAX_ll) .OR. (IJMAX/=KJMAX_ll)                      &
                             .OR. (IKMAX/=(IKU-2*JPVEXT)) ) THEN
  WRITE(ILUOUT,FMT=9003) 1
  WRITE(ILUOUT,*) 'THE GRIDS ARE DIFFERENT IN THE INITIAL FILE :'
  WRITE(ILUOUT,*) KIMAX_ll,'*',KJMAX_ll,'*',IKU-2*JPVEXT
  WRITE(ILUOUT,*) 'AND IN THE COUPLING FILE :',TCPLFILE(NCPL_CUR)%TZFILE%CNAME
  WRITE(ILUOUT,*) IIMAX,'*',IJMAX,'*',IKMAX
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_CPL','')
END IF
!
!*      3.2   Initialize  the large scale sources
!
GLSOURCE=.TRUE.
ZLENG = (NCPL_TIMES(NCPL_CUR,1)-2) * PTSTEP
!
CALL INI_LS(TCPLFILE(NCPL_CUR)%TZFILE,HGETRVM,GLSOURCE,PLSUS,PLSVS,PLSWS,PLSTHS,PLSRVS, &
            PLSZWSS,PDRYMASSS,PLSUM,PLSVM,PLSWM,PLSTHM,PLSRVM,PLSZWSM,PDRYMASST,ZLENG,  &
            OSTEADY_DMASS)
!
!
!*      3.2   Initialize  the LB sources
!
!YGETSVM(1:KSV) = HGETSVM(1:KSV)
!IF ( LUSECHEM .AND. (.NOT. OCH_INIT_FIELD) )  &
!   YGETSVM(NSV_CHEMBEG: NSV_CHEMEND) = 'INIT'
!IF (HCONF == 'RESTA')  THEN
#ifdef MNH_FOREFIRE
!   IF (LPASPOL) YGETSVM(NSV_PPBEG: NSV_PPEND) = 'INIT'
!   IF (LFOREFIRE) YGETSVM(NSV_FFBEG: NSV_FFEND) = 'INIT'
#endif
!   IF (NSV_USER /= 0) YGETSVM(1:NSV_USER) = 'INIT'
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
!END IF
 
GLSOURCE=.TRUE.
CALL INI_LB(TCPLFILE(NCPL_CUR)%TZFILE,GLSOURCE,KSV,                  &
     KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,              &
     KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                  &
     KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,          &
     HGETTKEM,HGETRVM,HGETRCM,HGETRRM,HGETRIM,HGETRSM,               &
     HGETRGM,HGETRHM,HGETSVM,                                        &
     PLBXUS,PLBXVS,PLBXWS,PLBXTHS,PLBXTKES,PLBXRS,PLBXSVS,           &
     PLBYUS,PLBYVS,PLBYWS,PLBYTHS,PLBYTKES,PLBYRS,PLBYSVS,           &
     PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,           &
     PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,           &
     ZLENG)
!
!*      3.5   Close the coupling file
!
CALL IO_FILE_CLOSE_ll(TCPLFILE(NCPL_CUR)%TZFILE)
!!-------------------------------------------------------------------------------
!
!*      6.    FORMATS
!             -------
!
9000  FORMAT(/,'NOTE  IN INI_CPL FOR MODEL ', I2, ' : ',/, &
             '--------------------------------')
9001  FORMAT(/,'CAUTION ERROR IN INI_CPL FOR MODEL ', I2,' : ',/, &
             '----------------------------------------' )
9002  FORMAT(/,'WARNING IN INI_CPL FOR MODEL ', I2,' : ',/, &
             '----------------------------------' )
9003  FORMAT(/,'FATAL ERROR IN INI_CPL FOR MODEL ', I2,' : ',/, &
             '--------------------------------------' )
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_CPL
