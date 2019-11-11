!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############################
      MODULE MODI_INI_DEEP_CONVECTION
!     ###############################
!
INTERFACE
!
      SUBROUTINE INI_DEEP_CONVECTION(TPINIFILE,HLUOUT,OINIDCONV,TPDTCUR,                &
                                     KCOUNTCONV,PDTHCONV,PDRVCONV,PDRCCONV,             &
                                     PDRICONV,PPRCONV,PPRSCONV,PPACCONV,                &
                                     PUMFCONV,PDMFCONV,PMFCONV,PPRLFLXCONV,PPRSFLXCONV, &
                                     PCAPE,KCLTOPCONV,KCLBASCONV,                       &
                                     TPDTDCONV, HGETSVCONV, PDSVCONV,                   &
                                     OCH_CONV_LINOX, PIC_RATE, PCG_RATE,                &
                                     PIC_TOTAL_NUMBER, PCG_TOTAL_NUMBER                 )
!
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_TIME
!
TYPE(TFILEDATA),        INTENT(IN) :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN) :: HLUOUT    ! name for output-listing
                                                !  of nested models
LOGICAL,                INTENT(IN) :: OINIDCONV ! switch to initialize or read
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR   ! Current date and time
CHARACTER (LEN=*),      INTENT(IN) :: HGETSVCONV ! GET indicator for SVCONV
!
TYPE (DATE_TIME),       INTENT(OUT):: TPDTDCONV ! date and time of the 
                                                ! last deep convection call
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCOUNTCONV! convective counter(recompute
                                                ! tendency or keep it
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDTHCONV  ! convective theta tendency (K/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDRVCONV  ! convective r_v tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDRCCONV  ! convective r_c tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDRICONV  ! convective r_i tendency (1/s)
REAL, DIMENSION(:,:),   INTENT(OUT):: PPRCONV   ! total (liquid+solid) surf.
                                                ! precipitation tendency (m/s)
REAL, DIMENSION(:,:),   INTENT(OUT):: PPRSCONV  ! solid surface
                                                ! precipitation tendency (m/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PUMFCONV  ! updraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDMFCONV  ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMFCONV   ! convective mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PPRLFLXCONV!liquid precip flux (m/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PPRSFLXCONV!solid  precip flux (m/s)
REAL, DIMENSION(:,:),   INTENT(OUT):: PCAPE     ! CAPE (J)
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCLTOPCONV! convective cloud top level 
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCLBASCONV! convective cloud base level
REAL, DIMENSION(:,:),   INTENT(OUT):: PPACCONV  ! accumulated convective
                                                ! precipitation (m)
REAL, DIMENSION(:,:,:,:),INTENT(OUT):: PDSVCONV ! conv. tracer tendencies (1/s)
LOGICAL,                INTENT(IN)    :: OCH_CONV_LINOX ! Flag to compute LiNOx
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PIC_TOTAL_NUMBER ! Total number of IC 
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCG_TOTAL_NUMBER ! Total number of CG
!
END SUBROUTINE INI_DEEP_CONVECTION
!
END INTERFACE
!
END MODULE MODI_INI_DEEP_CONVECTION
!     ###################################################################################
      SUBROUTINE INI_DEEP_CONVECTION(TPINIFILE,HLUOUT,OINIDCONV,TPDTCUR,                &
                                     KCOUNTCONV,PDTHCONV,PDRVCONV,PDRCCONV,             &
                                     PDRICONV,PPRCONV,PPRSCONV,PPACCONV,                &
                                     PUMFCONV,PDMFCONV,PMFCONV,PPRLFLXCONV,PPRSFLXCONV, &
                                     PCAPE,KCLTOPCONV,KCLBASCONV,                       &
                                     TPDTDCONV, HGETSVCONV, PDSVCONV,                   &
                                     OCH_CONV_LINOX, PIC_RATE, PCG_RATE,                &
                                     PIC_TOTAL_NUMBER, PCG_TOTAL_NUMBER                 )
!     ###################################################################################
!
!!**** Routine to initialize the convective tendencies and the
!!     convective counter
!!
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INI_DEEP_CONVECTION )
!!
!!    AUTHOR
!!    ------
!!	  P. Bechtold      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!      Pinty J.-P. 29/01/97  Unit conversion for PPACCONV and PPRCONV
!!      Bechtold P. 24/01/98  Initialisation of tracer tendencies
!!      Asencio N.  13/08/98   parallel code: ILENG no longer used
!!      Bechtold P. 11/12/98  Drop PWSUBCONV and add surf solid precip. +
!!                            diagnostics
!!      D.Gazen       22/01/01 use MODD_NSV and add names to scalar variables
!!      P.Jabouille   04/04/02 add PMFCONV used for subgrid condensation
!!                    for a correct restart this variable has to be writen in FM file
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_TIME
USE MODD_CONVPAR
USE MODD_CH_M9_n,         ONLY: CNAMES
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_ICE_C1R3_DESCR,  ONLY : C1R3NAMES
USE MODD_ELEC_DESCR,      ONLY : CELECNAMES
USE MODD_LG,              ONLY: CLGNAMES
USE MODD_NSV, ONLY : NSV,NSV_USER,NSV_CHEMBEG,NSV_CHEMEND,NSV_C2R2BEG,NSV_C2R2END, &
                     NSV_LGBEG,NSV_LGEND,NSV_LNOXBEG,NSV_LNOXEND, &
                     NSV_DSTBEG,NSV_DSTEND, NSV_AERBEG,NSV_AEREND, &
                     NSV_SLTBEG,NSV_SLTEND, NSV_PPBEG,NSV_PPEND, &
                     NSV_C1R3BEG,NSV_C1R3END, NSV_ELECBEG,NSV_ELECEND
USE MODD_CH_AEROSOL, ONLY : CAERONAMES
USE MODD_DUST, ONLY : CDUSTNAMES
USE MODD_SALT, ONLY : CSALTNAMES
!
USE MODE_FIELD
USE MODE_FM
USE MODE_FMREAD
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(TFILEDATA),        INTENT(IN) :: TPINIFILE ! Initial file
CHARACTER (LEN=*),      INTENT(IN) :: HLUOUT    ! name for output-listing
                                                !  of nested models
LOGICAL,                INTENT(IN) :: OINIDCONV ! switch to initialize or read
TYPE (DATE_TIME),       INTENT(IN) :: TPDTCUR   ! Current date and time
CHARACTER (LEN=*),      INTENT(IN) :: HGETSVCONV ! GET indicator for SVCONV
!
TYPE (DATE_TIME),       INTENT(OUT):: TPDTDCONV ! date and time of the 
                                                ! last deep convection call
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCOUNTCONV! convective counter(recompute
                                                ! tendency or keep it
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDTHCONV  ! convective theta tendency (K/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDRVCONV  ! convective r_v tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDRCCONV  ! convective r_c tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDRICONV  ! convective r_i tendency (1/s)
REAL, DIMENSION(:,:),   INTENT(OUT):: PPRCONV   ! total (liquid+solid) surf.
                                                ! precipitation tendency (m/s)
REAL, DIMENSION(:,:),   INTENT(OUT):: PPACCONV  ! accumulated convective
                                                ! precipitation (m)
REAL, DIMENSION(:,:),   INTENT(OUT):: PPRSCONV  ! solid surface
                                                ! precipitation tendency (m/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PUMFCONV  ! updraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PDMFCONV  ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PMFCONV   ! convective mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PPRLFLXCONV!liquid precip flux (m/s)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PPRSFLXCONV!solid  precip flux (m/s)
REAL, DIMENSION(:,:),   INTENT(OUT):: PCAPE     ! CAPE (J)
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCLTOPCONV! convective cloud top level 
INTEGER, DIMENSION(:,:),INTENT(OUT):: KCLBASCONV! convective cloud base level
REAL, DIMENSION(:,:,:,:),INTENT(OUT):: PDSVCONV ! conv. tracer tendencies (1/s)
LOGICAL,                INTENT(IN)    :: OCH_CONV_LINOX ! Flag to compute LiNOx
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PIC_TOTAL_NUMBER ! Total number of IC 
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCG_TOTAL_NUMBER ! Total number of CG
!
!
!
!*       0.2   declarations of local variables
!
!
INTEGER          :: IID
INTEGER          :: IRESP
INTEGER          :: JSV     ! number of tracers
TYPE(TFIELDDATA) :: TZFIELD
!
!-------------------------------------------------------------------------------
!
!*       1. INITIALIZE CONSTANTS USED IN DEEP CONVECTION PARAMETERIZATION
!	        -------------------------------------------------------------
!
! call of INI_CONVPAR is now in routine CONVECTION
!
!
!*       2. INITIALIZE CONVECTIVE TENDENCIES
!	        --------------------------------
!
PUMFCONV(:,:,:)  = 0.0
PDMFCONV(:,:,:)  = 0.0
PMFCONV(:,:,:)   = 0.0  ! warning, restart may be incorrect
PPRLFLXCONV(:,:,:)=0.0
PPRSFLXCONV(:,:,:)=0.0
PCAPE(:,:)       = 0.0
KCLTOPCONV(:,:)  = 0
KCLBASCONV(:,:)  = 0
!
IF ( OINIDCONV ) THEN
  TPDTDCONV        = TPDTCUR
  KCOUNTCONV(:,:)  = 1
  PDTHCONV(:,:,:)  = 0.0
  PDRVCONV(:,:,:)  = 0.0
  PDRCCONV(:,:,:)  = 0.0
  PDRICONV(:,:,:)  = 0.0
  PPRCONV(:,:)     = 0.0
  PPRSCONV(:,:)    = 0.0
  PPACCONV(:,:)    = 0.0
  PDSVCONV(:,:,:,:) = 0.0
  IF ( OCH_CONV_LINOX ) THEN
    PIC_RATE(:,:) = 0.
    PCG_RATE(:,:) = 0.
    PIC_TOTAL_NUMBER(:,:) = 0.
    PCG_TOTAL_NUMBER(:,:) = 0.
  END IF
!
ELSE
!
  CALL IO_READ_FIELD(TPINIFILE,'DTDCONV',  TPDTDCONV)
  CALL IO_READ_FIELD(TPINIFILE,'COUNTCONV',KCOUNTCONV)
  CALL IO_READ_FIELD(TPINIFILE,'DTHCONV',  PDTHCONV)
  CALL IO_READ_FIELD(TPINIFILE,'DRVCONV',  PDRVCONV)
  CALL IO_READ_FIELD(TPINIFILE,'DRCCONV',  PDRCCONV)
  CALL IO_READ_FIELD(TPINIFILE,'DRICONV',  PDRICONV)
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PRCONV',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm hour-1'
  CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PPRCONV)
  PPRCONV=PPRCONV/(1000.*3600.) ! conversion into m/s units
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PRSCONV',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm hour-1'
  CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PPRSCONV)
  PPRSCONV=PPRSCONV/(1000.*3600.) ! conversion into m/s units
!
  CALL FIND_FIELD_ID_FROM_MNHNAME('PACCONV',IID,IRESP)
  TZFIELD = TFIELDLIST(IID)
  TZFIELD%CUNITS = 'mm'
  CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PPACCONV)
  PPACCONV=PPACCONV/1000.       ! conversion into m unit
!
  IF ( OCH_CONV_LINOX ) THEN
    CALL IO_READ_FIELD(TPINIFILE,'IC_RATE',    PIC_RATE)
    CALL IO_READ_FIELD(TPINIFILE,'CG_RATE',    PCG_RATE)
    CALL IO_READ_FIELD(TPINIFILE,'IC_TOTAL_NB',PIC_TOTAL_NUMBER)
    CALL IO_READ_FIELD(TPINIFILE,'CG_TOTAL_NB',PCG_TOTAL_NUMBER)
  END IF
!
!
 SELECT CASE(HGETSVCONV)
  CASE('READ')
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 's-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = 1, NSV_USER
      WRITE(TZFIELD%CMNHNAME,'(A7,I3.3)')'DSVCONV',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_C2R2BEG, NSV_C2R2END
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_C1R3BEG, NSV_C1R3END
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_ELECBEG, NSV_ELECEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_CHEMBEG, NSV_CHEMEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(UPCASE(CNAMES(JSV-NSV_CHEMBEG+1)))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_AERBEG, NSV_AEREND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(UPCASE(CAERONAMES(JSV-NSV_AERBEG+1)))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_LNOXBEG,NSV_LNOXEND
      TZFIELD%CMNHNAME   = 'DSVCONV_LINOX'
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_DSTBEG, NSV_DSTEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(UPCASE(CDUSTNAMES(JSV-NSV_DSTBEG+1)))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_SLTBEG, NSV_SLTEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(UPCASE(CSALTNAMES(JSV-NSV_SLTBEG+1)))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_LGBEG, NSV_LGEND
      TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CLGNAMES(JSV-NSV_LGBEG+1))
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
    DO JSV = NSV_PPBEG, NSV_PPEND
      WRITE(TZFIELD%CMNHNAME,'(A7,I3.3)')'DSVCONV',JSV
      TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
      WRITE(TZFIELD%CCOMMENT,'(A6,A7,I3.3)')'X_Y_Z_','DSVCONV',JSV
      CALL IO_READ_FIELD(TPINIFILE,TZFIELD,PDSVCONV(:,:,:,JSV))
    END DO
 END SELECT
!
!
END IF
!  
CONTAINS
FUNCTION UPCASE(HSTRING)

CHARACTER(LEN=*)            :: HSTRING
CHARACTER(LEN=LEN(HSTRING)) :: UPCASE

INTEGER :: JC
INTEGER, PARAMETER :: IAMIN = IACHAR("a")
INTEGER, PARAMETER :: IAMAJ = IACHAR("A")

DO JC=1,LEN(HSTRING)
  IF (HSTRING(JC:JC) >= "a" .AND. HSTRING(JC:JC) <= "z") THEN
      UPCASE(JC:JC) = ACHAR(IACHAR(HSTRING(JC:JC)) - IAMIN + IAMAJ)
  ELSE
      UPCASE(JC:JC) = HSTRING(JC:JC)
  END IF
END DO

END FUNCTION UPCASE
!
END SUBROUTINE INI_DEEP_CONVECTION
