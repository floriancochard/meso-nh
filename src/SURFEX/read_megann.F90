!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_MEGAN_n(MSF, U, HPROGRAM)
!     #################################
!
!!****  *READ_MEGAN_n* - routine to READ dummy surface fields
!!
!!    PURPOSE
!!    -------
!!
!!    AUTHOR
!!    ------
!!      P. Tulet & M. Leriche  *LACy & LA*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     06/2017
!!      Modification 07/2017    J. Pianezze adaptation to SurfEx v8
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_MEGAN_SURF_FIELDS_n, ONLY : MEGAN_SURF_FIELDS_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODI_READ_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(MEGAN_SURF_FIELDS_t), INTENT(INOUT) :: MSF
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM     ! 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: JMEGAN         ! loop counter
!
 CHARACTER(LEN=20 ):: YSTRING20      ! string
 CHARACTER(LEN=3  ):: YSTRING03      ! string
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       2.     Number of dummy fields :
!               ----------------------
!
IF (LHOOK) CALL DR_HOOK('READ_MEGAN_N',0,ZHOOK_HANDLE)
!
YRECFM='MEGAN_GR_NBR'
YCOMMENT=' '
!
 CALL READ_SURF(HPROGRAM,YRECFM,MSF%NMEGAN_NBR,IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!*       3.     Dummy fields :
!               ------------
!
ALLOCATE(MSF%CMEGAN_NAME(MSF%NMEGAN_NBR))
ALLOCATE(MSF%CMEGAN_AREA(MSF%NMEGAN_NBR))
ALLOCATE(MSF%XMEGAN_FIELDS(U%NSIZE_FULL,MSF%NMEGAN_NBR))
MSF%CMEGAN_NAME(:) = ' '
MSF%CMEGAN_AREA(:) = ' '
!
!
DO JMEGAN=1,MSF%NMEGAN_NBR
  !
  WRITE(YRECFM,FMT='(A8,I3.3,A1)') 'MEGAN_GR',JMEGAN,' '
  CALL READ_SURF(HPROGRAM,YRECFM,MSF%XMEGAN_FIELDS(:,JMEGAN),IRESP,HCOMMENT=YCOMMENT)
  !
  !
  YSTRING20=YCOMMENT(21:40)
  YSTRING03=YCOMMENT(41:43)
  !
  MSF%CMEGAN_NAME(JMEGAN) = YSTRING20
  MSF%CMEGAN_AREA(JMEGAN) = YSTRING03
  WRITE(YRECFM,FMT='(A10,I2.2)') 'MEGAN_NAME',JMEGAN
  CALL READ_SURF(HPROGRAM,YRECFM,MSF%CMEGAN_NAME(JMEGAN),IRESP,HCOMMENT=YCOMMENT)
  !
END DO
!
IF (LHOOK) CALL DR_HOOK('READ_MEGAN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_MEGAN_n
