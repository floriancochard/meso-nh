!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITESURF_MEGAN_n(HSELECT, MSF, HPROGRAM)
!     ##########################################
!
!!****  *WRITESURF_MEGAN_n* - routine to write dummy surface fields
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
!!      Original      06/2017
!!      J. Pianezze   07/2017  adapatation tu SurfEx v8.0
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_MEGAN_SURF_FIELDS_n, ONLY : MEGAN_SURF_FIELDS_t
!
USE MODI_WRITE_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
TYPE(MEGAN_SURF_FIELDS_t), INTENT(INOUT) :: MSF
 CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM     ! 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: JMEGAN         ! loop counter
!
CHARACTER(LEN=20) :: YSTRING20      ! string
CHARACTER(LEN=3 ) :: YSTRING03      ! string
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=LEN_HREC) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     Number of megan fields :
!               ----------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_MEGAN_N',0,ZHOOK_HANDLE)
!
YRECFM='MEGAN_GR_NBR'
YCOMMENT=' '
!
 CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,MSF%NMEGAN_NBR,IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!*       2.     MEGAN fields :
!               ------------
!
DO JMEGAN=1,MSF%NMEGAN_NBR
  !
  WRITE(YRECFM,'(A8,I3.3,A1)') 'MEGAN_GR',JMEGAN,' '
  YSTRING20=MSF%CMEGAN_NAME(JMEGAN)
  YSTRING03=MSF%CMEGAN_AREA(JMEGAN)
  YCOMMENT='X_Y_'//YRECFM//YSTRING20//YSTRING03//  &
             '                                                             '  
  CALL WRITE_SURF(HSELECT,HPROGRAM,YRECFM,MSF%XMEGAN_FIELDS(:,JMEGAN),IRESP,HCOMMENT=YCOMMENT)

  WRITE(YRECFM,'(A10,I2.2)') 'MEGAN_NAME',JMEGAN
  CALL WRITE_SURF(HSELECT, HPROGRAM,YRECFM,MSF%CMEGAN_NAME(JMEGAN),IRESP,HCOMMENT=YCOMMENT)
END DO
IF (LHOOK) CALL DR_HOOK('WRITESURF_MEGAN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_MEGAN_n
