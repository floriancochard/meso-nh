!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications
! G. TANGUY 19/05/2014 : correctoin DATIME in case of time average
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------
!#######################
MODULE MODE_LES_DIACHRO
!#######################
!
USE MODD_LUNIT
!
CONTAINS
!
!---------------------------------------------------------------------
!
!########################################################
SUBROUTINE MAKE_NORM(HUNIT, KC, ODIV, PA_NORM, PLES_NORM)
!########################################################
!
USE MODD_PARAMETERS
USE MODD_LES
IMPLICIT NONE
!
CHARACTER(LEN=50),       INTENT(IN)    :: HUNIT     ! physical unit of field
INTEGER,                 INTENT(IN)    :: KC        ! character counter
LOGICAL,                 INTENT(INOUT) :: ODIV      ! flag to make a division
REAL, DIMENSION(:,:,:,:),INTENT(INOUT) :: PA_NORM   ! normalized field
REAL, DIMENSION(:),      INTENT(IN)    :: PLES_NORM ! normalization coefficient

INTEGER :: JK ! z counter
INTEGER :: JT ! time counter
INTEGER :: JP ! process counter
INTEGER :: JN ! variable number counter (larger than 1 only for scalar var.)
CHARACTER(LEN=50) :: YUNIT
!
!------------------------------------------------------------------------
!
YUNIT=HUNIT
!
IF ( ANY(PLES_NORM(:)==0.) ) THEN
  PA_NORM(:,:,:,:)=XUNDEF
  ODIV=.FALSE.
  RETURN
END IF

DO JN=1,SIZE(PA_NORM,4)
  IF (YUNIT(KC+1:KC+1)=='3') THEN
    IF (ODIV) THEN
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JN)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JN) = PA_NORM(JK,JT,JP,JN) * PLES_NORM(JT)**3
          END DO
        END DO
      END DO
    ELSE
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JN)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JN) = PA_NORM(JK,JT,JP,JN) / PLES_NORM(JT)**3
          END DO
        END DO
      END DO
    END IF
  ELSE IF  (YUNIT(KC+1:KC+1)=='2') THEN
    IF (ODIV) THEN
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JN)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JN) = PA_NORM(JK,JT,JP,JN) * PLES_NORM(JT)**2
          END DO
        END DO
      END DO
    ELSE
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JN)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JN) = PA_NORM(JK,JT,JP,JN) / PLES_NORM(JT)**2
          END DO
        END DO
      END DO
    END IF
  ELSE
    IF (ODIV) THEN
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JN)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JN) = PA_NORM(JK,JT,JP,JN) * PLES_NORM(JT)
          END DO
        END DO
      END DO
    ELSE
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JN)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JN) = PA_NORM(JK,JT,JP,JN) / PLES_NORM(JT)
          END DO
        END DO
      END DO
    END IF
  END IF
END DO
!
ODIV=.FALSE.
!
END SUBROUTINE MAKE_NORM
!
!###########################################################
SUBROUTINE MAKE_NORM_SV(HUNIT, KC, ODIV, PA_NORM, PLES_NORM)
!###########################################################
!
USE MODD_PARAMETERS
USE MODD_LES
IMPLICIT NONE
!
CHARACTER(LEN=50),       INTENT(IN)    :: HUNIT     ! physical unit of field
INTEGER,                 INTENT(IN)    :: KC        ! character counter
LOGICAL,                 INTENT(INOUT) :: ODIV      ! flag to make a division
REAL, DIMENSION(:,:,:,:),INTENT(INOUT) :: PA_NORM   ! normalized field
REAL, DIMENSION(:,:),    INTENT(IN)    :: PLES_NORM ! normalization coefficient

INTEGER :: JK ! z counter
INTEGER :: JT ! time counter
INTEGER :: JP ! process counter
INTEGER :: JSV! scalar variables counter
CHARACTER(LEN=50) :: YUNIT
!
!-------------------------------------------------------------------
!
YUNIT=HUNIT
!
DO JSV=1,SIZE(PLES_NORM,2)
  IF (ANY(PLES_NORM(:,JSV)==0.)) THEN
    PA_NORM(:,:,:,JSV)=XUNDEF
    CYCLE
  END IF
  IF (YUNIT(KC+1:KC+1)=='3') THEN
    IF (ODIV) THEN
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JSV)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JSV) = PA_NORM(JK,JT,JP,JSV) * PLES_NORM(JT,JSV)**3
          END DO
        END DO
      END DO
    ELSE
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JSV)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JSV) = PA_NORM(JK,JT,JP,JSV) / PLES_NORM(JT,JSV)**3
          END DO
        END DO
      END DO
    END IF
  ELSE IF  (YUNIT(KC+1:KC+1)=='2') THEN
    IF (ODIV) THEN
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JSV)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JSV) = PA_NORM(JK,JT,JP,JSV) * PLES_NORM(JT,JSV)**2
            END DO
        END DO
      END DO
    ELSE
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JSV)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JSV) = PA_NORM(JK,JT,JP,JSV) / PLES_NORM(JT,JSV)**2
            END DO
        END DO
      END DO
    END IF
  ELSE
    IF (ODIV) THEN
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JSV)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JSV) = PA_NORM(JK,JT,JP,JSV) * PLES_NORM(JT,JSV)
            END DO
        END DO
      END DO
    ELSE
      DO JP=1,SIZE(PA_NORM,3)
        DO JT=1,NLES_CURRENT_TIMES
          DO JK=1,NLES_K
            IF (PA_NORM(JK,JT,JP,JSV)/=XUNDEF) &
            PA_NORM(JK,JT,JP,JSV) = PA_NORM(JK,JT,JP,JSV) / PLES_NORM(JT,JSV)
          END DO
        END DO
      END DO
    END IF
  END IF
END DO
!
ODIV=.FALSE.
!
END SUBROUTINE MAKE_NORM_SV
!
!     ###################################################
      SUBROUTINE LES_NORM_4D(HUNIT, PA_LES, PA_NORM, OSV)
!     ###################################################
!
!
!!****  *LES_NORM* normalizes a field according to the chosen normalization
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         30/10/00
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_LES
!
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
CHARACTER(LEN=*),            INTENT(IN)  :: HUNIT   ! physical unit of field
REAL,    DIMENSION(:,:,:,:), INTENT(IN)  :: PA_LES  ! field
!
REAL,    DIMENSION(:,:,:,:), INTENT(OUT) :: PA_NORM ! normalized field
LOGICAL, OPTIONAL,           INTENT(IN)  :: OSV     ! flag for scalar variables
!
!       0.2  declaration of local variables
!
INTEGER :: JC    ! character counter
LOGICAL :: GDIV  ! flag to make a division
INTEGER :: IKG   ! number of 'kg' in the field unit
!
CHARACTER(LEN=50) :: YUNIT
!
!------------------------------------------------------------------------------
YUNIT=HUNIT//'                                              '
!
PA_NORM = PA_LES
!
IKG=0
!
DO JC=1,50
  IF (YUNIT(JC:JC)=='g') THEN
    IKG=IKG+1
  ELSE IF (YUNIT(JC:JC)==' ') THEN
    EXIT
  END IF
END DO
!
GDIV=.FALSE.
!
DO JC=1,49
  !
  SELECT CASE (YUNIT(JC:JC))
    CASE('m')
      CALL MAKE_NORM(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_M)
    CASE('g')
      IF (IKG==1) THEN
        CALL MAKE_NORM(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_RHO)
      ELSE
        IF (PRESENT(OSV)) THEN
          IF (OSV) THEN
            IF (.NOT. GDIV) &
            CALL MAKE_NORM_SV(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_SV)
            GDIV=.FALSE.
          ELSE
            IF (.NOT. GDIV) &
            CALL MAKE_NORM(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_RV)
            GDIV=.FALSE.
          END IF
        ELSE 
          IF (.NOT. GDIV) &
          CALL MAKE_NORM(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_RV)
          GDIV=.FALSE.
        END IF
      END IF
    CASE('K')
      CALL MAKE_NORM(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_K)
    CASE('s')
      CALL MAKE_NORM(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_S)
    CASE('a')
      CALL MAKE_NORM(YUNIT, JC, GDIV, PA_NORM, XLES_NORM_P)
    CASE('/')
      GDIV=.TRUE.
    CASE(' ')
      EXIT
  END SELECT
END DO
!
WHERE(PA_NORM==XUNDEF) PA_NORM = PA_LES
!
END SUBROUTINE LES_NORM_4D
!
!------------------------------------------------------------------------------
!
!     ###################################################
      SUBROUTINE LES_NORM_3D(HUNIT, PA_LES, PA_NORM, OSV)
!     ###################################################
!
!
!!****  *LES_NORM* normalizes a field according to the chosen normalization
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         30/10/00
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
CHARACTER(LEN=*),          INTENT(IN)  :: HUNIT   ! physical unit of field
REAL,    DIMENSION(:,:,:), INTENT(IN)  :: PA_LES  ! field
!
REAL,    DIMENSION(:,:,:), INTENT(OUT) :: PA_NORM ! normalized field
LOGICAL, OPTIONAL,         INTENT(IN)  :: OSV     ! flag for scalar variables
!
!       0.2  declaration of local variables
!
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZA_LES
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZA_NORM
!
INTEGER :: JSV
!------------------------------------------------------------------------------
!
IF (PRESENT(OSV)) THEN
  IF (OSV) THEN
    ALLOCATE(ZA_LES (SIZE(PA_LES,1),SIZE(PA_LES,2),1,SIZE(PA_LES,3)))
    ALLOCATE(ZA_NORM(SIZE(PA_LES,1),SIZE(PA_LES,2),1,SIZE(PA_LES,3)))
    DO JSV=1,SIZE(PA_LES,3)
      ZA_LES (:,:,1,JSV) = PA_LES(:,:,JSV)
    END DO
    CALL LES_NORM_4D(HUNIT, ZA_LES, ZA_NORM, OSV)
    DO JSV=1,SIZE(PA_LES,3)
      PA_NORM(:,:,JSV) = ZA_NORM(:,:,1,JSV)
    END DO
  ELSE
    ALLOCATE(ZA_LES (SIZE(PA_LES,1),SIZE(PA_LES,2),SIZE(PA_LES,3),1))
    ALLOCATE(ZA_NORM(SIZE(PA_LES,1),SIZE(PA_LES,2),SIZE(PA_LES,3),1))
    ZA_LES (:,:,:,1) = PA_LES(:,:,:)
    CALL LES_NORM_4D(HUNIT, ZA_LES, ZA_NORM, OSV)
    PA_NORM(:,:,:) = ZA_NORM(:,:,:,1)
  END IF
ELSE
  ALLOCATE(ZA_LES (SIZE(PA_LES,1),SIZE(PA_LES,2),SIZE(PA_LES,3),1))
  ALLOCATE(ZA_NORM(SIZE(PA_LES,1),SIZE(PA_LES,2),SIZE(PA_LES,3),1))
  ZA_LES (:,:,:,1) = PA_LES(:,:,:)
  CALL LES_NORM_4D(HUNIT, ZA_LES, ZA_NORM)
  PA_NORM(:,:,:) = ZA_NORM(:,:,:,1)
END IF
!
DEALLOCATE(ZA_LES)
DEALLOCATE(ZA_NORM)
!
END SUBROUTINE LES_NORM_3D
!
!     ##############################################
      SUBROUTINE LES_NORM_2D(HUNIT, PA_LES, PA_NORM)
!     ##############################################
!
!
!!****  *LES_NORM* normalizes a field according to the chosen normalization
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         30/10/00
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN)  :: HUNIT   ! physical unit of field
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA_LES  ! field
!
REAL,    DIMENSION(:,:), INTENT(OUT) :: PA_NORM ! normalized field
!
!       0.2  declaration of local variables
!
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZA_LES
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZA_NORM
!
!------------------------------------------------------------------------------
!
ALLOCATE(ZA_LES  (SIZE(PA_LES,1),SIZE(PA_LES,2),1,1))
ALLOCATE(ZA_NORM (SIZE(PA_LES,1),SIZE(PA_LES,2),1,1))
!
ZA_LES (:,:,1,1) = PA_LES(:,:)
CALL LES_NORM_4D(HUNIT, ZA_LES, ZA_NORM)
PA_NORM(:,:) = ZA_NORM(:,:,1,1)
!
DEALLOCATE(ZA_LES)
DEALLOCATE(ZA_NORM)
!
END SUBROUTINE LES_NORM_2D
!
!------------------------------------------------------------------------------
!
!###################################
SUBROUTINE LES_Z_NORM(OAVG,PTRAJZ,PWORK6)
!###################################
!
!* this subroutine interpolates the normalized field PWORK6 to the
!  vertical normalized coordinate.
!
USE MODD_PARAMETERS, ONLY : XUNDEF, JPVEXT
USE MODD_LES
!
USE MODI_COEF_VER_INTERP_LIN
USE MODI_VER_INTERP_LIN
!
IMPLICIT NONE
!
LOGICAL,                      INTENT(IN)    :: OAVG   ! time average ?
REAL, DIMENSION(:,:,:),       INTENT(INOUT) :: PTRAJZ ! vertical grid
REAL, DIMENSION(:,:,:,:,:,:), INTENT(INOUT) :: PWORK6 ! field data array
!
REAL,    DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3))  :: ZTRAJZ
! initial grid
REAL,    DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3), &
                   SIZE(PWORK6,1),SIZE(PWORK6,2),SIZE(PWORK6,6) ) :: ZWORK6
! initial data
REAL,    DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3))  :: ZNORMZ
! grid in normalized coordinates
!
REAL,    DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3), &
                   SIZE(PWORK6,1),SIZE(PWORK6,2),SIZE(PWORK6,6) ) :: ZNORM6
! data interpolated on normalized vertical grid
!
REAL,    DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3)+2*JPVEXT)  :: ZZ
REAL,    DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3)+2*JPVEXT)  :: ZW
INTEGER, DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3))           :: IKLIN
REAL,    DIMENSION(SIZE(PWORK6,4),SIZE(PWORK6,5),SIZE(PWORK6,3))           :: ZCOEFLIN
!
INTEGER :: JK, JT, JN, JP
INTEGER :: ITEMP_MEAN_START, ITEMP_MEAN_END
!
REAL    :: ZMAX_NORM_M
!------------------------------------------------------------------------------
!
!* normalization height (usually maximum BL height)
!
IF (OAVG) THEN
  ITEMP_MEAN_START = COUNT( XLES_CURRENT_TRAJT(:,1)<=XLES_TEMP_MEAN_START ) + 1
  ITEMP_MEAN_END   = COUNT( XLES_CURRENT_TRAJT(:,1)<=XLES_TEMP_MEAN_END   )
  IF (ITEMP_MEAN_START > ITEMP_MEAN_END) THEN
    ITEMP_MEAN_START = 1
    ITEMP_MEAN_END   = NLES_CURRENT_TIMES
  END IF
ELSE
  ITEMP_MEAN_START = 1
  ITEMP_MEAN_END   = NLES_CURRENT_TIMES
END IF
!
ZMAX_NORM_M = MAXVAL(ABS(XLES_NORM_M(ITEMP_MEAN_START:ITEMP_MEAN_END)))
!
IF (ZMAX_NORM_M<=0.) THEN
  PWORK6(:,:,:,:,:,:) = XUNDEF
  RETURN
END IF
!
!* moves K index in third position
!
DO JK=1,NLES_K
  DO JN=1,SIZE(PWORK6,5)
    DO JT=1,SIZE(PWORK6,4)
      ZTRAJZ(JT,JN,JK) = PTRAJZ(JK,1,JN)
      ZWORK6(JT,JN,JK,1,1,:) = PWORK6(1,1,JK,JT,JN,:)
    END DO
  END DO
END DO
!
!* computes the stretching due to the use of a normalized grid
!
DO JK=1,NLES_K
  DO JN=1,SIZE(PWORK6,5)
    DO JT=1,SIZE(PWORK6,4)
      ZNORMZ(JT,JN,JK) = (ZTRAJZ(JT,JN,JK)-XLES_CURRENT_ZS)/ZMAX_NORM_M * XLES_NORM_M(JT)  &
                         + XLES_CURRENT_ZS
    END DO
  END DO
END DO
!
!
IF (NLES_K>1) THEN
!
!* computes the interpolation coefficients
!
  ZZ(:,:,JPVEXT+1:JPVEXT+NLES_K) = ZTRAJZ(:,:,:)
  DO JK=1,JPVEXT
    ZZ(:,:,              JK) = ZTRAJZ(:,:,1)      - (ZTRAJZ(:,:,2)     -ZTRAJZ(:,:,1)       )*(JPVEXT+1-JK)
    ZZ(:,:,NLES_K+JPVEXT+JK) = ZTRAJZ(:,:,NLES_K) + (ZTRAJZ(:,:,NLES_K)-ZTRAJZ(:,:,NLES_K-1))*          JK
  END DO

  CALL COEF_VER_INTERP_LIN(ZZ,ZNORMZ,IKLIN,ZCOEFLIN)
!
!* performs the interpolation
!
  DO JP=1,SIZE(PWORK6,6)
    ZW = XUNDEF
    ZW(:,:,JPVEXT+1:JPVEXT+NLES_K) = ZWORK6(:,:,:,1,1,JP)
    ZNORM6(:,:,:,1,1,JP) = VER_INTERP_LIN(ZW,IKLIN,ZCOEFLIN)
  END DO
!
ELSE
  ZNORM6 = ZWORK6
END IF
!
!* puts the normalized grid and data in diachro arrays
!
PTRAJZ(:,:,:) = (PTRAJZ(:,:,:)-XLES_CURRENT_ZS)/ZMAX_NORM_M
!
DO JN=1,SIZE(PWORK6,5)
  DO JT=1,SIZE(PWORK6,4)
    DO JK=1,NLES_K
      PWORK6(1,1,JK,JT,JN,:) = ZNORM6(JT,JN,JK,1,1,:)
    END DO
  END DO
END DO
!
END SUBROUTINE LES_Z_NORM
!------------------------------------------------------------------------------
!
!########################################################
SUBROUTINE LES_TIME_AVG(PTRAJT,PWORK6,KRESP,PDATIME_AVG)
!########################################################
!
! this routine computes time averaging
!
! Modifications:
!  03/2018     (P.Wautelet)   replace ADD_FORECAST_TO_DATE by DATETIME_CORRECTDATE
!
USE MODD_LES
USE MODD_TYPE_DATE
!
USE MODE_DATETIME
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:),         POINTER     :: PTRAJT ! time
REAL, DIMENSION(:,:),         POINTER     :: PDATIME_AVG ! date
REAL, DIMENSION(:,:,:,:,:,:), POINTER     :: PWORK6 ! contains physical field
INTEGER,                      INTENT(OUT) :: KRESP  ! return code (0 is OK)
!------------------------------------------------------------------------------
INTEGER                                :: JT       ! time counter
INTEGER                                :: ITIME    ! nb of avg. points
INTEGER                                :: IAVG     ! nb of avg. periods
INTEGER                                :: JAVG     ! loop counter on avg. periods
REAL                                   :: ZLES_TEMP_MEAN_START ! initial and end times
REAL                                   :: ZLES_TEMP_MEAN_END   ! of one avergaing preiod
REAL, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: ZWORK6     ! contains averaged physical field
INTEGER :: JK                            ! vertical loop counter
INTEGER :: JP                            ! process loop counter
INTEGER :: JSV                           ! scalar loop counter
INTEGER :: JX                            ! first  spatial or spectral coordinate loop counter
INTEGER :: JY                            ! second spatial or spectral coordinate loop counter
REAL, DIMENSION(16)  :: ZDATIME_SAVE ! date
TYPE(DATE_TIME) :: TZDATE
!------------------------------------------------------------------------------
!
IF (     XLES_TEMP_MEAN_END==XUNDEF   &
    .OR. XLES_TEMP_MEAN_START==XUNDEF &
    .OR. XLES_TEMP_MEAN_STEP==XUNDEF  ) THEN
  KRESP=-1
  RETURN
END IF
!
IAVG=INT(XLES_TEMP_MEAN_END-1.E-10-XLES_TEMP_MEAN_START)/XLES_TEMP_MEAN_STEP + 1
IF (IAVG<=0) THEN
  KRESP=-1
  RETURN
END IF
!
ZDATIME_SAVE(:)=PDATIME_AVG(:,1)
DEALLOCATE(PTRAJT)
DEALLOCATE(PDATIME_AVG)
!
ALLOCATE (PTRAJT(IAVG,1))
ALLOCATE (PDATIME_AVG(16,IAVG))
ALLOCATE (ZWORK6(SIZE(PWORK6,1),SIZE(PWORK6,2),NLES_K,IAVG,SIZE(PWORK6,5),SIZE(PWORK6,6)))
!
ZWORK6(:,:,:,:,:,:) = 0.
!
PDATIME_AVG(1,:)=ZDATIME_SAVE(1)
PDATIME_AVG(2,:)=ZDATIME_SAVE(2)
PDATIME_AVG(3,:)=ZDATIME_SAVE(3)
PDATIME_AVG(4,:)=ZDATIME_SAVE(4)
PDATIME_AVG(5,:)=ZDATIME_SAVE(5)
PDATIME_AVG(6,:)=ZDATIME_SAVE(6)
PDATIME_AVG(7,:)=ZDATIME_SAVE(7)
PDATIME_AVG(8,:)=ZDATIME_SAVE(8)
PDATIME_AVG(9,:)=ZDATIME_SAVE(9)
PDATIME_AVG(10,:)=ZDATIME_SAVE(10)
PDATIME_AVG(11,:)=ZDATIME_SAVE(11)
PDATIME_AVG(12,:)=ZDATIME_SAVE(12)
!
DO JAVG=1,IAVG
  ZLES_TEMP_MEAN_START=XLES_TEMP_MEAN_START + (JAVG-1) * XLES_TEMP_MEAN_STEP
  ZLES_TEMP_MEAN_END  =MIN(XLES_TEMP_MEAN_END, ZLES_TEMP_MEAN_START + XLES_TEMP_MEAN_STEP)
  !
  DO JP=1,SIZE(PWORK6,6)
    DO JSV=1,SIZE(PWORK6,5)
      DO JK=1,SIZE(PWORK6,3)
        DO JY=1,SIZE(PWORK6,2)
          DO JX=1,SIZE(PWORK6,1)
            ITIME=0
            DO JT=1,NLES_CURRENT_TIMES
              IF ( XLES_CURRENT_TRAJT(JT,1) >= ZLES_TEMP_MEAN_START .AND. &
                 XLES_CURRENT_TRAJT(JT,1) <= ZLES_TEMP_MEAN_END) THEN
                IF (PWORK6(JX,JY,JK,JT,JSV,JP) /= XUNDEF) THEN
                 ZWORK6(JX,JY,JK,JAVG,JSV,JP) =  ZWORK6(JX,JY,JK,JAVG,JSV,JP) &
                                              + PWORK6(JX,JY,JK,JT,JSV,JP)
                 ITIME=ITIME+1
                END IF
              END IF
            END DO
            IF (ITIME >= 1) THEN
                    ZWORK6(JX,JY,JK,JAVG,JSV,JP)= &
                            ZWORK6(JX,JY,JK,JAVG,JSV,JP) / ITIME
            END IF
            IF (ITIME == 0) THEN
                   ZWORK6(JX,JY,JK,JAVG,JSV,JP)= XUNDEF 
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
  PTRAJT(JAVG,1)=(ZLES_TEMP_MEAN_START+ZLES_TEMP_MEAN_END)/2.
  TZDATE%TDATE%YEAR  = PDATIME_AVG(5,JAVG)
  TZDATE%TDATE%MONTH = PDATIME_AVG(6,JAVG)
  TZDATE%TDATE%DAY   = PDATIME_AVG(7,JAVG)
  TZDATE%TIME        = PDATIME_AVG(8,JAVG)+PTRAJT(JAVG,1)
  CALL DATETIME_CORRECTDATE(TZDATE)
  PDATIME_AVG(13,JAVG) = TZDATE%TDATE%YEAR
  PDATIME_AVG(14,JAVG) = TZDATE%TDATE%MONTH
  PDATIME_AVG(15,JAVG) = TZDATE%TDATE%DAY
  PDATIME_AVG(16,JAVG) = TZDATE%TIME
END DO
!
DEALLOCATE(PWORK6)
ALLOCATE(PWORK6(SIZE(ZWORK6,1),SIZE(ZWORK6,2),NLES_K,IAVG,SIZE(ZWORK6,5),SIZE(ZWORK6,6)))
PWORK6 = ZWORK6
DEALLOCATE(ZWORK6)

KRESP = 0

END SUBROUTINE LES_TIME_AVG
!------------------------------------------------------------------------------
!
!########################################################
SUBROUTINE LES_DIACHRO(TPDIAFILE,HGROUP,HCOMMENT,HUNIT,PFIELD,HAVG)
!########################################################
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),      INTENT(IN)       :: TPDIAFILE! file to write
CHARACTER(LEN=*),     INTENT(IN)       :: HGROUP   ! group title
CHARACTER(LEN=*),     INTENT(IN)       :: HCOMMENT ! comment string
CHARACTER(LEN=*),     INTENT(IN)       :: HUNIT    ! physical unit
REAL, DIMENSION(:,:), INTENT(IN)       :: PFIELD
CHARACTER(LEN=1),     INTENT(IN)       :: HAVG     ! flag to compute avg.
!
!*      0.2  declaration of local variables for diachro
!
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJX ! localization of the temporal
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJY ! series in x,y and z. remark:
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJZ ! x and y are not used for LES
REAL,    DIMENSION(:,:),   POINTER     :: ZTRAJT ! time
REAL,    DIMENSION(:,:),   POINTER     :: ZDATIME ! date
!
INTEGER, DIMENSION(1)                  :: IGRID    ! grid indicator
CHARACTER(LEN= 10)                     :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(1)       :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(1)       :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(1)       :: YUNIT    ! physical unit
!
REAL,    DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2))    &
                                       :: ZFIELD     ! normalized field
INTEGER                                :: IRESP    ! return code
!
REAL, DIMENSION(:,:,:,:,:,:), POINTER  :: ZWORK6 ! contains physical field
!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
INTEGER :: JK                            ! vertical loop counter
!
LOGICAL :: GAVG                          ! flag to compute time averagings
LOGICAL :: GNORM                         ! flag to compute normalizations
!
!-------------------------------------------------------------------------------
!
GAVG =(HAVG=='A' .OR. HAVG=='H')
GNORM=(HAVG=='E' .OR. HAVG=='H')
!
IF (GAVG .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
!
ZFIELD=PFIELD
IF (GNORM) CALL LES_NORM_2D(HUNIT, PFIELD, ZFIELD)
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
ALLOCATE (ZTRAJX(1,1,1))
ALLOCATE (ZTRAJY(1,1,1))
ALLOCATE (ZTRAJZ(NLES_K,1,1))
!
ALLOCATE(ZWORK6(1,1,NLES_K,NLES_CURRENT_TIMES,1,1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))
!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
ZTRAJX(:,:,:) = (IIL+IIH)/2
ZTRAJY(:,:,:) = (IJL+IJH)/2
IKL=NLES_LEVELS(1)
IKH=NLES_LEVELS(NLES_K)
DO JK=1,NLES_K
  ZTRAJZ(JK,1,1) = XLES_CURRENT_Z(JK)
END DO
IGRID(1)=1
YCOMMENT(1) = HCOMMENT
!
YUNIT (1) = HUNIT
YGROUP    = HGROUP
!
ZWORK6(1,1,:,:,1,1) = ZFIELD (:,:)
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)=XLES_CURRENT_DATIME(:,:)
!
!* normalization of vertical dimension
!
IF (GNORM) THEN
  IF (HUNIT(1:1)/=' ') YUNIT='-'
  CALL LES_Z_NORM(GAVG,ZTRAJZ,ZWORK6)
END IF
!
!* time average
!
IRESP = 0
IF (GAVG) CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
!
IF (HAVG/=' ')  YGROUP=HAVG//'_'//YGROUP
YTITLE(1) = YGROUP
!
!*      2.0  Writing of the profile
!            ----------------------
!
IF (IRESP==0 .AND. ANY(ZWORK6/=XUNDEF)) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SSOL",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH,                                    &
                   PTRAJX=ZTRAJX,PTRAJY=ZTRAJY,PTRAJZ=ZTRAJZ                   )
!
!
!*      3.0  Deallocations
!            -------------
!
DEALLOCATE (ZTRAJX)
DEALLOCATE (ZTRAJY)
DEALLOCATE (ZTRAJZ)
DEALLOCATE (ZTRAJT)
DEALLOCATE (ZWORK6)
DEALLOCATE (ZDATIME)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO
!-------------------------------------------------------------------------------
!###########################################################
SUBROUTINE LES_DIACHRO_SV(TPDIAFILE,HGROUP,HCOMMENT,HUNIT,PFIELD,HAVG)
!###########################################################
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),        INTENT(IN)       :: TPDIAFILE! file to write
CHARACTER(LEN=*),       INTENT(IN)       :: HGROUP   ! group title
CHARACTER(LEN=*),       INTENT(IN)       :: HCOMMENT ! comment string
CHARACTER(LEN=*),       INTENT(IN)       :: HUNIT    ! physical unit
REAL, DIMENSION(:,:,:), INTENT(IN)       :: PFIELD
CHARACTER(LEN=1),       INTENT(IN)       :: HAVG     ! flag to compute avg.
!
!*      0.2  declaration of local variables for diachro
!
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJX ! localization of the temporal
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJY ! series in x,y and z. remark:
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJZ ! x and y are not used for LES
REAL,    DIMENSION(:,:),   POINTER     :: ZTRAJT ! time
REAL,    DIMENSION(:,:),   POINTER     :: ZDATIME ! date
!
INTEGER, DIMENSION(1)                  :: IGRID    ! grid indicator
CHARACTER(LEN= 10)                     :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(1)       :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(1)       :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(1)       :: YUNIT    ! physical unit
REAL,    DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3))    &
                                       :: ZFIELD     ! normalized field
INTEGER                                :: IRESP    ! return code
!
REAL, DIMENSION(:,:,:,:,:,:), POINTER  :: ZWORK6 ! contains physical field
!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
INTEGER :: JK                            ! vertical loop counter
INTEGER :: JSV                           ! scalar loop counter
!
LOGICAL :: GAVG                          ! flag to compute time averagings
LOGICAL :: GNORM                         ! flag to compute normalizations
!
!-------------------------------------------------------------------------------
!
GAVG =(HAVG=='A' .OR. HAVG=='H')
GNORM=(HAVG=='E' .OR. HAVG=='H')
!
IF (GAVG .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
!
ZFIELD=PFIELD
IF (GNORM) CALL LES_NORM_3D(HUNIT, PFIELD, ZFIELD, .TRUE.)
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
ALLOCATE (ZTRAJX(1,1,SIZE(PFIELD,3)))
ALLOCATE (ZTRAJY(1,1,SIZE(PFIELD,3)))
ALLOCATE (ZTRAJZ(NLES_K,1,SIZE(PFIELD,3)))
ALLOCATE(ZWORK6(1,1,NLES_K,NLES_CURRENT_TIMES,SIZE(PFIELD,3),1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))
!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
ZTRAJX(:,:,:) = (IIL+IIH)/2
ZTRAJY(:,:,:) = (IJL+IJH)/2
IKL=NLES_LEVELS(1)
IKH=NLES_LEVELS(NLES_K)
DO JK=1,NLES_K
  ZTRAJZ(JK,:,:) = XLES_CURRENT_Z(JK)
END DO
IGRID(1)=1
YCOMMENT(1) = HCOMMENT
!
YUNIT (1) = HUNIT
YGROUP    = HGROUP
!
ZWORK6(1,1,:,:,:,1) = ZFIELD (:,:,:)
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)=XLES_CURRENT_DATIME(:,:)
!
IF (GNORM) THEN
  IF (HUNIT(1:1)/=' ') YUNIT='-'
  CALL LES_Z_NORM(GAVG,ZTRAJZ,ZWORK6)
END IF
!
!* time average
!
IRESP = 0
IF (GAVG) CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
!
IF (HAVG/=' ')  YGROUP=HAVG//'_'//YGROUP
YTITLE(1) = YGROUP
!
!*      2.0  Writing of the profile
!            ----------------------
!
!
IF (IRESP==0 .AND. ANY(ZWORK6/=XUNDEF)) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SSOL",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH,                                    &
                   PTRAJX=ZTRAJX,PTRAJY=ZTRAJY,PTRAJZ=ZTRAJZ                   )
!
!
!*      3.0  Deallocations
!            -------------
!
DEALLOCATE (ZTRAJX)
DEALLOCATE (ZTRAJY)
DEALLOCATE (ZTRAJZ)
DEALLOCATE (ZTRAJT)
DEALLOCATE(ZWORK6)
DEALLOCATE(ZDATIME)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO_SV
!-------------------------------------------------------------------------------
!#####################################################################
SUBROUTINE LES_DIACHRO_MASKS(TPDIAFILE,HGROUP,HTITLE,HCOMMENT,HUNIT,PFIELD,HAVG)
!#####################################################################
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),                    INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=*),                   INTENT(IN) :: HGROUP   ! group title
CHARACTER(LEN=*), DIMENSION(:),     INTENT(IN) :: HTITLE   ! title
CHARACTER(LEN=*), DIMENSION(:),     INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=*),                   INTENT(IN) :: HUNIT    ! physical unit
REAL,             DIMENSION(:,:,:), INTENT(IN) :: PFIELD
CHARACTER(LEN=1),                   INTENT(IN) :: HAVG     ! flag to compute avg.
!
!*      0.2  declaration of local variables for diachro
!
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJX ! localization of the temporal
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJY ! series in x,y and z. remark:
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJZ ! x and y are not used for LES
REAL,    DIMENSION(:,:),   POINTER     :: ZTRAJT ! time
REAL,    DIMENSION(:,:),   POINTER     :: ZDATIME ! date
!
INTEGER,            DIMENSION(SIZE(PFIELD,3)) :: IGRID    ! grid indicator
CHARACTER(LEN= 10)                            :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(SIZE(PFIELD,3)) :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(SIZE(PFIELD,3)) :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(SIZE(PFIELD,3)) :: YUNIT    ! physical unit
REAL,    DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3))    &
                                              :: ZFIELD     ! normalized field
INTEGER                                       :: IRESP    ! return code
!
REAL, DIMENSION(:,:,:,:,:,:), POINTER         :: ZWORK6 ! contains physical field
!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
INTEGER :: JK                            ! vertical loop counter
INTEGER :: JMASK                         ! Mask loop counter
!
LOGICAL :: GAVG                          ! flag to compute time averagings
LOGICAL :: GNORM                         ! flag to compute normalizations
!
!-------------------------------------------------------------------------------
!
GAVG =(HAVG=='A' .OR. HAVG=='H')
GNORM=(HAVG=='E' .OR. HAVG=='H')
!
IF (GAVG .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
!
ZFIELD=PFIELD
IF (GNORM) CALL LES_NORM_3D(HUNIT, PFIELD, ZFIELD)
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
ALLOCATE (ZTRAJX(1,1,1))
ALLOCATE (ZTRAJY(1,1,1))
ALLOCATE (ZTRAJZ(NLES_K,1,1))
ALLOCATE(ZWORK6(1,1,NLES_K,NLES_CURRENT_TIMES,1,SIZE(PFIELD,3)))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))

!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
ZTRAJX(:,:,:) = (IIL+IIH)/2
ZTRAJY(:,:,:) = (IJL+IJH)/2
IKL=NLES_LEVELS(1)
IKH=NLES_LEVELS(NLES_K)
DO JK=1,NLES_K
  ZTRAJZ(JK,:,1) = XLES_CURRENT_Z(JK)
END DO
!
IGRID(:)=1
!
YCOMMENT(:) = HCOMMENT(:)
YUNIT   (:) = HUNIT
YGROUP      = HGROUP
!
ZWORK6(1,1,:,:,1,:) = ZFIELD (:,:,:)
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)= XLES_CURRENT_DATIME(:,:)

!
IF (GNORM) THEN
  IF (HUNIT(1:1)/=' ') YUNIT='-'
  CALL LES_Z_NORM(GAVG,ZTRAJZ,ZWORK6)
END IF
!
!
!* time average
!
IRESP = 0
IF (GAVG) CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
!
IF (HAVG/=' ')  YGROUP=HAVG//'_'//YGROUP
YTITLE  (:) = YGROUP//HTITLE(:)
!
!
!*      2.0  Writing of the profile
!            ----------------------
!
IF (IRESP==0 .AND. ANY(ZWORK6/=XUNDEF)) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SSOL",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH,                                    &
                   PTRAJX=ZTRAJX,PTRAJY=ZTRAJY,PTRAJZ=ZTRAJZ                   )
!
!
!*      3.0  Deallocations
!            -------------
!
DEALLOCATE (ZTRAJX)
DEALLOCATE (ZTRAJY)
DEALLOCATE (ZTRAJZ)
DEALLOCATE (ZTRAJT)
DEALLOCATE(ZWORK6)
DEALLOCATE(ZDATIME)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO_MASKS
!-------------------------------------------------------------------------------
!########################################################################
SUBROUTINE LES_DIACHRO_SV_MASKS(TPDIAFILE,HGROUP,HTITLE,HCOMMENT,HUNIT,PFIELD,HAVG)
!########################################################################
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),                      INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=*),                     INTENT(IN) :: HGROUP   ! group title
CHARACTER(LEN=*), DIMENSION(:),       INTENT(IN) :: HTITLE   ! title
CHARACTER(LEN=*), DIMENSION(:),       INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=*),                     INTENT(IN) :: HUNIT    ! physical unit
REAL,             DIMENSION(:,:,:,:), INTENT(IN) :: PFIELD
CHARACTER(LEN=1),                     INTENT(IN) :: HAVG     ! flag to compute avg.
!
!*      0.2  declaration of local variables for diachro
!
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJX ! localization of the temporal
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJY ! series in x,y and z. remark:
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJZ ! x and y are not used for LES
REAL,    DIMENSION(:,:),   POINTER     :: ZTRAJT ! time
REAL,    DIMENSION(:,:),   POINTER     :: ZDATIME! date
!
INTEGER,            DIMENSION(SIZE(PFIELD,3)) :: IGRID    ! grid indicator
CHARACTER(LEN= 10)                            :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(SIZE(PFIELD,3)) :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(SIZE(PFIELD,3)) :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(SIZE(PFIELD,3)) :: YUNIT    ! physical unit
REAL, DIMENSION(SIZE(PFIELD,1),SIZE(PFIELD,2),SIZE(PFIELD,3),SIZE(PFIELD,4))&
                                              :: ZFIELD   ! normalized field
INTEGER                                       :: IRESP    ! return code
!
REAL, DIMENSION(:,:,:,:,:,:), POINTER         :: ZWORK6 ! contains physical field
!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
INTEGER :: JK                            ! vertical loop counter
INTEGER :: JP                            ! process loop counter
INTEGER :: JSV                           ! scalar loop counter
INTEGER :: JMASK                         ! mask loop counter
!
LOGICAL :: GAVG                          ! flag to compute time averagings
LOGICAL :: GNORM                         ! flag to compute normalizations
!
!-------------------------------------------------------------------------------
!
GAVG =(HAVG=='A' .OR. HAVG=='H')
GNORM=(HAVG=='E' .OR. HAVG=='H')
!
IF (GAVG .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
!
ZFIELD=PFIELD
IF (GNORM) CALL LES_NORM_4D(HUNIT, PFIELD, ZFIELD, .TRUE.)
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
ALLOCATE (ZTRAJX(1,1,SIZE(PFIELD,4)))
ALLOCATE (ZTRAJY(1,1,SIZE(PFIELD,4)))
ALLOCATE (ZTRAJZ(NLES_K,1,SIZE(PFIELD,4)))
ALLOCATE(ZWORK6(1,1,NLES_K,NLES_CURRENT_TIMES,SIZE(PFIELD,4),SIZE(PFIELD,3)))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))
!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
ZTRAJX(:,:,:) = (IIL+IIH)/2
ZTRAJY(:,:,:) = (IJL+IJH)/2
IKL=NLES_LEVELS(1)
IKH=NLES_LEVELS(NLES_K)
DO JK=1,NLES_K
  ZTRAJZ(JK,:,:) = XLES_CURRENT_Z(JK)
END DO
IGRID(:)=1
!
YCOMMENT(:) = HCOMMENT(:)
YUNIT   (:) = HUNIT
YGROUP      = HGROUP
!
!
DO JSV=1,SIZE(PFIELD,4)
  DO JP=1,SIZE(PFIELD,3)
    ZWORK6(1,1,:,:,JSV,JP) = ZFIELD (:,:,JP,JSV)
  END DO
END DO
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)= XLES_CURRENT_DATIME(:,:)
!
IF (GNORM) THEN
  IF (HUNIT(1:1)/=' ') YUNIT='-'
  CALL LES_Z_NORM(GAVG,ZTRAJZ,ZWORK6)
END IF
!n
!
!* time average
!
IRESP = 0
IF (GAVG) CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
!
IF (HAVG/=' ')  YGROUP=HAVG//'_'//YGROUP
YTITLE  (:) = YGROUP//HTITLE(:)
!
!*      2.0  Writing of the profile
!            ----------------------
!
!
IF (IRESP==0 .AND. ANY(ZWORK6/=XUNDEF)) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SSOL",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH,                                    &
                   PTRAJX=ZTRAJX,PTRAJY=ZTRAJY,PTRAJZ=ZTRAJZ                   )
!
!
!*      3.0  Deallocations
!            -------------
!
DEALLOCATE (ZTRAJX)
DEALLOCATE (ZTRAJY)
DEALLOCATE (ZTRAJZ)
DEALLOCATE (ZTRAJT)
DEALLOCATE(ZWORK6)
DEALLOCATE(ZDATIME)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO_SV_MASKS
!-------------------------------------------------------------------------------

!#############################################################
SUBROUTINE LES_DIACHRO_SURF(TPDIAFILE,HGROUP,HCOMMENT,HUNIT,PFIELD,HAVG)
!#############################################################
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),      INTENT(IN)       :: TPDIAFILE! file to write
CHARACTER(LEN=*),     INTENT(IN)       :: HGROUP   ! group title
CHARACTER(LEN=*),     INTENT(IN)       :: HCOMMENT ! comment string
CHARACTER(LEN=*),     INTENT(IN)       :: HUNIT    ! physical unit
REAL, DIMENSION(:),   INTENT(IN)       :: PFIELD
CHARACTER(LEN=1),     INTENT(IN)       :: HAVG     ! flag to compute avg.
!
!
!*      0.2  declaration of local variables for diachro
!
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJX ! localization of the temporal
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJY ! series in x,y and z. remark:
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJZ ! x and y are not used for LES
REAL,    DIMENSION(:,:),   POINTER     :: ZTRAJT ! time
REAL,    DIMENSION(:,:),   POINTER     :: ZDATIME ! DATE
!
INTEGER, DIMENSION(1)                  :: IGRID    ! grid indicator
CHARACTER(LEN= 10)                     :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(1)       :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(1)       :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(1)       :: YUNIT    ! physical unit
INTEGER                                :: IRESP    ! return code
!
REAL, DIMENSION(:,:,:,:,:,:), POINTER  :: ZWORK6 ! contains physical field
!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
!
LOGICAL :: GAVG                          ! flag to compute time averagings
LOGICAL :: GNORM                         ! flag to compute normalizations
!-------------------------------------------------------------------------------
!
GAVG =(HAVG=='A' .OR. HAVG=='H')
GNORM=(HAVG=='E' .OR. HAVG=='H')
!

IF (GAVG .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
!
IF (GNORM) RETURN
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
ALLOCATE (ZTRAJX(1,1,1))
ALLOCATE (ZTRAJY(1,1,1))
ALLOCATE (ZTRAJZ(1,1,1))
ALLOCATE(ZWORK6(1,1,1,NLES_CURRENT_TIMES,1,1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))
!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
ZTRAJX(:,:,:) = (IIL+IIH)/2
ZTRAJY(:,:,:) = (IJL+IJH)/2
IKL=NLES_LEVELS(1)
IKH=NLES_LEVELS(1)
ZTRAJZ(1,1,1) = XLES_CURRENT_Z(1)
IGRID(1)=1
YCOMMENT(1) = HCOMMENT
!
YUNIT (1) = HUNIT
YGROUP    = HGROUP
!
ZWORK6(1,1,1,:,1,1) = PFIELD (:)
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)=XLES_CURRENT_DATIME(:,:)
!
!* time average
!
IRESP = 0
IF (GAVG) CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
!
IF (HAVG/=' ')  YGROUP=HAVG//'_'//YGROUP
YTITLE(1) = HGROUP
!
!*      2.0  Writing of the profile
!            ----------------------
!
IF (IRESP==0) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SSOL",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH,                                    &
                   PTRAJX=ZTRAJX,PTRAJY=ZTRAJY,PTRAJZ=ZTRAJZ                   )
!
!
!*      3.0  Deallocations
!            -------------
!
DEALLOCATE (ZTRAJX)
DEALLOCATE (ZTRAJY)
DEALLOCATE (ZTRAJZ)
DEALLOCATE (ZTRAJT)
DEALLOCATE(ZWORK6)
DEALLOCATE(ZDATIME)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO_SURF
!-------------------------------------------------------------------------------
!################################################################
SUBROUTINE LES_DIACHRO_SURF_SV(TPDIAFILE,HGROUP,HCOMMENT,HUNIT,PFIELD,HAVG)
!################################################################
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),        INTENT(IN)       :: TPDIAFILE! file to write
CHARACTER(LEN=*),       INTENT(IN)       :: HGROUP   ! group title
CHARACTER(LEN=*),       INTENT(IN)       :: HCOMMENT ! comment string
CHARACTER(LEN=*),       INTENT(IN)       :: HUNIT    ! physical unit
REAL, DIMENSION(:,:),   INTENT(IN)       :: PFIELD
CHARACTER(LEN=1),       INTENT(IN)       :: HAVG     ! flag to compute avg.
!
!*      0.2  declaration of local variables for diachro
!
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJX ! localization of the temporal
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJY ! series in x,y and z. remark:
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: ZTRAJZ ! x and y are not used for LES
REAL,    DIMENSION(:,:),   POINTER     :: ZTRAJT ! time
REAL,    DIMENSION(:,:),   POINTER     :: ZDATIME ! date
INTEGER, DIMENSION(1)                  :: IGRID    ! grid indicator
CHARACTER(LEN= 10)                     :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(1)       :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(1)       :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(1)       :: YUNIT    ! physical unit
INTEGER                                :: IRESP    ! return code
!
REAL, DIMENSION(:,:,:,:,:,:), POINTER  :: ZWORK6 ! contains physical field
!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
!
LOGICAL :: GAVG                          ! flag to compute time averagings
LOGICAL :: GNORM                         ! flag to compute normalizations
!-------------------------------------------------------------------------------
!
GAVG =(HAVG=='A' .OR. HAVG=='H')
GNORM=(HAVG=='E' .OR. HAVG=='H')
!
IF (GAVG .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
!
IF (GNORM) RETURN
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
ALLOCATE (ZTRAJX(1,1,SIZE(PFIELD,2)))
ALLOCATE (ZTRAJY(1,1,SIZE(PFIELD,2)))
ALLOCATE (ZTRAJZ(1,1,SIZE(PFIELD,2)))
ALLOCATE(ZWORK6(1,1,1,NLES_CURRENT_TIMES,SIZE(PFIELD,2),1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))

!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
ZTRAJX(:,:,:) = (IIL+IIH)/2
ZTRAJY(:,:,:) = (IJL+IJH)/2
IKL=NLES_LEVELS(1)
IKH=NLES_LEVELS(1)
ZTRAJZ(1,1,:) = XLES_CURRENT_Z(1)
IGRID(1)=1
YCOMMENT(1) = HCOMMENT
!
YUNIT (1) = HUNIT
YGROUP    = HGROUP
!
IRESP = 0
ZWORK6(1,1,1,:,:,1) = PFIELD (:,:)
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)=XLES_CURRENT_DATIME(:,:)
!

!
!* time average
!
IF (GAVG) CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
!
!
IF (HAVG/=' ')  YGROUP=HAVG//'_'//YGROUP
YTITLE(1) = HGROUP
!
!*      2.0  Writing of the profile
!            ----------------------
!
IF (IRESP==0) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SSOL",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH,                                    &
                   PTRAJX=ZTRAJX,PTRAJY=ZTRAJY,PTRAJZ=ZTRAJZ                   )
!
!
!*      3.0  Deallocations
!            -------------
!
DEALLOCATE (ZTRAJX)
DEALLOCATE (ZTRAJY)
DEALLOCATE (ZTRAJZ)
DEALLOCATE (ZTRAJT)
DEALLOCATE(ZWORK6)
DEALLOCATE(ZDATIME)
!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO_SURF_SV
!-------------------------------------------------------------------------------
!#####################################################################
SUBROUTINE LES_DIACHRO_2PT(TPDIAFILE,HGROUP,HCOMMENT,HUNIT,PFIELDX,PFIELDY,HAVG)
!#####################################################################
!
!* Modification 01/04/03 (V. Masson) safer use of ZWORK6 with loops
!
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODD_CONF
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),                    INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=*),                   INTENT(IN) :: HGROUP   ! group title
CHARACTER(LEN=*),                   INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=*),                   INTENT(IN) :: HUNIT    ! physical unit
REAL,             DIMENSION(:,:,:), INTENT(IN) :: PFIELDX
REAL,             DIMENSION(:,:,:), INTENT(IN) :: PFIELDY
CHARACTER(LEN=1),                   INTENT(IN) :: HAVG     ! flag to compute avg.
!
!*      0.2  declaration of local variables for diachro
!
!
INTEGER,            DIMENSION(1) :: IGRID    ! grid indicator
CHARACTER(LEN= 10)               :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(1) :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(1) :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(1) :: YUNIT    ! physical unit
REAL, DIMENSION(SIZE(PFIELDX,1),SIZE(PFIELDX,2)) :: ZAVG_FIELDX
REAL, DIMENSION(SIZE(PFIELDY,1),SIZE(PFIELDY,2)) :: ZAVG_FIELDY
INTEGER                          :: JT       ! time counter
INTEGER                          :: JK       ! level counter
INTEGER                          :: IRESP    ! return code
REAL, DIMENSION(:,:),POINTER     :: ZTRAJT   ! time
REAL, DIMENSION(:,:),POINTER     :: ZDATIME   ! date

!
REAL, DIMENSION(:,:,:,:,:,:), POINTER :: ZWORK6 ! contains physical field
!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
!
CHARACTER(len=6) :: YSTRING
!
LOGICAL :: GAVG                          ! flag to compute time averagings
!-------------------------------------------------------------------------------
!
IF (HAVG/=' '.AND. HAVG/='A') RETURN
!
GAVG=(HAVG=='A')
!
IF (GAVG .AND. (XLES_TEMP_MEAN_START==XUNDEF .OR. XLES_TEMP_MEAN_END==XUNDEF)) RETURN
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
ALLOCATE(ZWORK6(SIZE(PFIELDX,1),1,NSPECTRA_K,NLES_CURRENT_TIMES,2,1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))
!
IGRID(:)=1
!
YUNIT (:) = HUNIT
!
IKL=1
IKH=NSPECTRA_K
!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = 1
IJH = 1
!
YGROUP    = 'CI_'//HGROUP
YTITLE(:) = YGROUP
WRITE(YSTRING,FMT="(I6.6)") NINT( XLES_CURRENT_DOMEGAX )
YCOMMENT(:) = " DOMEGAX="//YSTRING//' '//HCOMMENT
!
IRESP = 0
DO JT=1,SIZE(PFIELDX,3)
    DO JK=1,SIZE(PFIELDX,2)
      ZWORK6(:,1,JK,JT,1,1) = PFIELDX (:,JK,JT)
      ZWORK6(:,1,JK,JT,2,1) = 0.
    END DO
END DO
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)= XLES_CURRENT_DATIME(:,:)
!* time average
!
IF (GAVG) THEN
  CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
  YGROUP    = 'T_'//YGROUP
END IF
!
!
!*      2.0  Writing of the profile
!            ----------------------
!
IF (IRESP==0) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SPXY",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH                                     )
!
!
DEALLOCATE (ZTRAJT)
DEALLOCATE(ZDATIME)
DEALLOCATE(ZWORK6)
!
IF (L2D) RETURN
!
ALLOCATE(ZWORK6(1,SIZE(PFIELDY,1),NSPECTRA_K,NLES_CURRENT_TIMES,2,1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))
!
IIL = 1
IIH = 1
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
!
DO JT=1,SIZE(PFIELDY,3)
    DO JK=1,SIZE(PFIELDY,2)
      ZWORK6(1,:,JK,JT,1,1) = PFIELDY (:,JK,JT)
      ZWORK6(1,:,JK,JT,2,1) = 0.
    END DO
END DO
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)= XLES_CURRENT_DATIME(:,:)
!
YGROUP    = 'CJ_'//HGROUP
YTITLE(:) = YGROUP
WRITE(YSTRING,FMT="(I6.6)") NINT( XLES_CURRENT_DOMEGAY )
YCOMMENT(:) = " DOMEGAY="//YSTRING//' '//HCOMMENT
!
!
!* time average
!
IF (GAVG) THEN
  CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
  YGROUP    = 'T_'//YGROUP
END IF
!
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SPXY",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH                                     )
!
DEALLOCATE (ZTRAJT)
DEALLOCATE(ZWORK6)
DEALLOCATE(ZDATIME)

!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO_2PT
!------------------------------------------------------------------------------
!
!#####################################################################
SUBROUTINE LES_DIACHRO_SPEC(TPDIAFILE,HGROUP,HCOMMENT,HUNIT,PSPECTRAX,PSPECTRAY)
!#####################################################################
!
!* Modification 01/04/03 (V. Masson) safer use of ZWORK6 with loops
!
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LES
USE MODD_GRID
USE MODD_CONF
USE MODI_WRITE_DIACHRO
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(TFILEDATA),                      INTENT(IN) :: TPDIAFILE! file to write
CHARACTER(LEN=*),                     INTENT(IN) :: HGROUP   ! group title
CHARACTER(LEN=*),                     INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=*),                     INTENT(IN) :: HUNIT    ! physical unit
REAL,             DIMENSION(:,:,:,:), INTENT(IN) :: PSPECTRAX! spectra in x
REAL,             DIMENSION(:,:,:,:), INTENT(IN) :: PSPECTRAY! and y directions
!
!*      0.2  declaration of local variables for diachro
!
!
INTEGER,            DIMENSION(1) :: IGRID    ! grid indicator
CHARACTER(LEN= 10)               :: YGROUP   ! group title
CHARACTER(LEN=100), DIMENSION(1) :: YCOMMENT ! comment string
CHARACTER(LEN=100), DIMENSION(1) :: YTITLE   ! title
CHARACTER(LEN=100), DIMENSION(1) :: YUNIT    ! physical unit
INTEGER                          :: IRESP    ! return code
!
REAL, DIMENSION(:,:,:,:,:,:), POINTER     :: ZWORK6 ! contains physical field
REAL, DIMENSION(:,:),         POINTER     :: ZTRAJT ! time
REAL, DIMENSION(:,:),         POINTER     :: ZDATIME ! date

!
INTEGER :: IIL, IIH, IJL, IJH, IKL, IKH  ! cartesian area relatively to the
!                                        ! entire domain
!
CHARACTER(len=6) :: YSTRING
INTEGER          :: JT       ! time counter
INTEGER          :: JK       ! level counter
!
!-------------------------------------------------------------------------------
!
!*      1.0  Initialization of diachro variables for LES (z,t) profiles
!            ----------------------------------------------------------
!
IGRID(:)=1
!
YUNIT (:) = HUNIT
!
IKL=1
IKH=NSPECTRA_K
!
!*      2.0  Writing of the profile
!            ----------------------
!* spectra in X direction
!
ALLOCATE(ZWORK6(SIZE(PSPECTRAX,1),1,NSPECTRA_K,NLES_CURRENT_TIMES,2,1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))

!
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)= XLES_CURRENT_DATIME(:,:)
!
IIL = NLES_CURRENT_IINF
IIH = NLES_CURRENT_ISUP
IJL = 1
IJH = 1
!
DO JT=1,SIZE(PSPECTRAX,4)
  DO JK=1,SIZE(PSPECTRAX,3)
    ZWORK6(:,1,JK,JT,1,1) = PSPECTRAX (:,1,JK,JT)
    ZWORK6(:,1,JK,JT,2,1) = PSPECTRAX (:,2,JK,JT)
  END DO
END DO
!
YGROUP    = 'SI_'//HGROUP
YTITLE(:) = YGROUP
WRITE(YSTRING,FMT="(I6.6)") NINT( XLES_CURRENT_DOMEGAX )
YCOMMENT(:) = " DOMEGAX="//YSTRING//' '//HCOMMENT
!
!
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SPXY",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH                                     )
!
!
!* time average
!
IRESP=0
CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
YGROUP    = 'T_'//YGROUP
!
IF (IRESP==0) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SPXY",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH                                     )
DEALLOCATE(ZWORK6)
DEALLOCATE(ZTRAJT)
DEALLOCATE(ZDATIME)
!
!* spectra in Y direction
!

IF (L2D) RETURN
!
ALLOCATE(ZWORK6(1,SIZE(PSPECTRAY,1),NSPECTRA_K,NLES_CURRENT_TIMES,2,1))
ALLOCATE(ZTRAJT(NLES_CURRENT_TIMES,1))
ALLOCATE(ZDATIME(16,NLES_CURRENT_TIMES))
!
ZTRAJT(:,:) = XLES_CURRENT_TRAJT(:,:)
ZDATIME(:,:)= XLES_CURRENT_DATIME(:,:)
!
IIL = 1
IIH = 1
IJL = NLES_CURRENT_JINF
IJH = NLES_CURRENT_JSUP
!
DO JT=1,SIZE(PSPECTRAY,4)
  DO JK=1,SIZE(PSPECTRAY,3)
    ZWORK6(1,:,JK,JT,1,1) = PSPECTRAY (:,1,JK,JT)
    ZWORK6(1,:,JK,JT,2,1) = PSPECTRAY (:,2,JK,JT)
  END DO
END DO
!
YGROUP    = 'SJ_'//HGROUP
YTITLE(:) = YGROUP
WRITE(YSTRING,FMT="(I6.6)") NINT( XLES_CURRENT_DOMEGAY )
YCOMMENT(:) = " DOMEGAY="//YSTRING//' '//HCOMMENT
!
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SPXY",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH                                     )
!
!
!* time average
!
CALL LES_TIME_AVG(ZTRAJT,ZWORK6,IRESP,ZDATIME)
YGROUP    = 'T_'//YGROUP
!
IF (IRESP==0) &
CALL WRITE_DIACHRO(TPDIAFILE,TLUOUT0,YGROUP,"SPXY",IGRID,ZDATIME, ZWORK6,     &
                   ZTRAJT,YTITLE,YUNIT,YCOMMENT,.FALSE.,.FALSE.,.FALSE.,   &
                   IIL,IIH,IJL,IJH,IKL,IKH                                     )                   
!
DEALLOCATE(ZWORK6)
DEALLOCATE(ZTRAJT)
DEALLOCATE(ZDATIME)

!
!-------------------------------------------------------------------------------
END SUBROUTINE LES_DIACHRO_SPEC

!-------------------------------------------------------------------------------
END MODULE MODE_LES_DIACHRO
