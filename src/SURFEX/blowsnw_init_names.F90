!!    ###########################################
SUBROUTINE BLOWSNW_INIT_NAMES (KLUOUT,HSV,KSNWEQ,KSV_SNWBEG,KSV_SNWEND,    &
                             K2DSNWEQ,K2D_SNWBEG,K2D_SNWEND   )
!!    ###########################################
!!
!!*** *SNW_INIT_NAMES*
!!
!!    PURPOSE
!!    -------
!!      Read and filter all chemical species into the CSV array
!!     initialize NSV_SNWBEG and  NSV_SNWEND index for the begin and the ending chemical index
!!     
!!
!!    REFERENCE
!!    ---------
!!    Modified dst_init_names (february 2005)    
!!
!!    AUTHOR
!!    ------
!!    Vincent VIONNET <vincent.vionnet@meteo.fr>
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!

INTEGER,                         INTENT(IN)  :: KLUOUT   ! output listing channel
CHARACTER(LEN=*), DIMENSION(:),  INTENT(IN)  :: HSV      ! name of chemical species
                                                         ! with character # for chemistry
INTEGER,                         INTENT(OUT) :: KSNWEQ     ! number of blowing snow related variables
INTEGER,                         INTENT(OUT) :: KSV_SNWBEG  ! first blowing snow related scalar
INTEGER,                         INTENT(OUT) :: KSV_SNWEND  ! last  blowing snow related scalar
INTEGER,                         INTENT(OUT) :: K2DSNWEQ     ! number of 2D blowing snow related variables
INTEGER,                         INTENT(OUT) :: K2D_SNWBEG  ! first 2D blowing snow related scalar
INTEGER,                         INTENT(OUT) :: K2D_SNWEND  ! last  2D blowing snow related scalar
!
!*      0.2    declarations of local variables
INTEGER :: JSV  !! loop on scalar variables
CHARACTER(LEN=4) :: YRC1
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

!Initialize output variables
IF (LHOOK) CALL DR_HOOK('BLOWSNW_INIT_NAMES',0,ZHOOK_HANDLE)

KSNWEQ = 0
KSV_SNWBEG=0
KSV_SNWEND=0

K2DSNWEQ = 0
K2D_SNWBEG=0
K2D_SNWEND=0

DO JSV=1, SIZE(HSV)
  YRC1= HSV(JSV)(1:4)
  IF (YRC1 == 'SNWM') THEN
     KSNWEQ = KSNWEQ + 1
     IF (KSNWEQ == 1) THEN
        KSV_SNWBEG=JSV
     ENDIF !Check on first time
  ELSE IF (YRC1 == 'SNWC') THEN
     K2DSNWEQ = K2DSNWEQ + 1
     IF (K2DSNWEQ == 1) THEN
        K2D_SNWBEG=JSV
     ENDIF !Check on first time
  ELSE !Not snow variables
     !DO NOTHING
  ENDIF
ENDDO
!
! Set the output list of scalar to the input list of scalars

! Get the index of the last blowing snow relevant tracer
KSV_SNWEND = KSV_SNWBEG + KSNWEQ -1
! Get the index of the last 2D blowing snow relevant tracer
K2D_SNWEND = K2D_SNWBEG + K2DSNWEQ -1

!
IF (LHOOK) CALL DR_HOOK('BLOWSNW_INIT_NAMES',1,ZHOOK_HANDLE)
!
END SUBROUTINE BLOWSNW_INIT_NAMES
