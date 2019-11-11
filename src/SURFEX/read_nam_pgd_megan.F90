!SURFEX_LIC Copyright 1994-2014 Meteo-France 
!SURFEX_LIC This is part of the SURFEX software governed by the CeCILL-C  licence
!SURFEX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!SURFEX_LIC for details. version 1.
!     #########
      SUBROUTINE READ_NAM_PGD_MEGAN(HPROGRAM, KMEGAN_NBR, HMEGAN_NAME, HMEGAN_AREA, &
                                    HMEGAN_ATYPE, HMEGAN_FILE, HMEGAN_FILETYPE      )  
!     ##############################################################
!
!!**** *READ_NAM_PGD_MEGAN* reads namelist NAM_MEGAN_PGD
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
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
!!
!!      P. Tulet & M. Leriche  *LACy & LA*
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    06/2017
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!
USE MODE_POS_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),                   INTENT(IN)   :: HPROGRAM     ! Type of program
INTEGER,                            INTENT(OUT)  :: KMEGAN_NBR
!                          ! number of megan pgd fields chosen by user
 CHARACTER(LEN=20), DIMENSION(1000), INTENT(OUT)  :: HMEGAN_NAME
!                          ! name of the megan pgd fields (for information)
 CHARACTER(LEN=3),  DIMENSION(1000), INTENT(OUT)  :: HMEGAN_AREA
!                          ! areas where megan pgd fields are defined
!                          ! 'ALL' : everywhere
!                          ! 'SEA' : where sea exists
!                          ! 'LAN' : where land exists
!                          ! 'WAT' : where inland water exists
!                          ! 'NAT' : where natural or agricultural areas exist
!                          ! 'TWN' : where town areas exist
!                          ! 'STR' : where streets are present
!                          ! 'BLD' : where buildings are present
 CHARACTER(LEN=3),  DIMENSION(1000), INTENT(OUT)  :: HMEGAN_ATYPE    ! avg type for megan pgd fields
!                                                                   ! 'ARI' , 'INV'
 CHARACTER(LEN=28), DIMENSION(1000), INTENT(OUT)  :: HMEGAN_FILE     ! data files
 CHARACTER(LEN=6),  DIMENSION(1000), INTENT(OUT)  :: HMEGAN_FILETYPE ! type of these files
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                           :: ILUOUT    ! output listing logical unit
INTEGER                           :: ILUNAM    ! namelist file logical unit
LOGICAL                           :: GFOUND    ! flag when namelist is present
!
!*    0.3    Declaration of namelists
!            ------------------------
!
INTEGER                             :: NMEGAN_NBR
!                          ! number of megan pgd fields chosen by user
 CHARACTER(LEN=20), DIMENSION(1000)  :: CMEGAN_NAME
!                          ! name of the megan pgd fields (for information)
 CHARACTER(LEN=3),  DIMENSION(1000)  :: CMEGAN_AREA
!                          ! areas where megan pgd fields are defined
!                          ! 'ALL' : everywhere
!                          ! 'SEA' : where sea exists
!                          ! 'LAN' : where land exists
!                          ! 'WAT' : where inland water exists
!                          ! 'NAT' : where natural or agricultural areas exist
!                          ! 'TWN' : where town areas exist
!                          ! 'STR' : where streets are present
!                          ! 'BLD' : where buildings are present
 CHARACTER(LEN=3),  DIMENSION(1000)  :: CMEGAN_ATYPE    ! avg type for megan pgd fields
!                                                      ! 'ARI' , 'INV'
 CHARACTER(LEN=28), DIMENSION(1000)  :: CMEGAN_FILE     ! data files
 CHARACTER(LEN=6),  DIMENSION(1000)  :: CMEGAN_FILETYPE ! type of these files
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_MEGAN_PGD/ NMEGAN_NBR, CMEGAN_NAME, CMEGAN_AREA,       &
                          CMEGAN_ATYPE, CMEGAN_FILE, CMEGAN_FILETYPE  
!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_MEGAN',0,ZHOOK_HANDLE)
NMEGAN_NBR = 0
!
CMEGAN_NAME     = "                    "
CMEGAN_FILE     = "                            "
CMEGAN_FILETYPE = "      "
CMEGAN_AREA     = "ALL"
CMEGAN_ATYPE    = "ARI"
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
 CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
!
 CALL POSNAM(ILUNAM,'NAM_MEGAN_PGD',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_MEGAN_PGD)
!
 CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
!-------------------------------------------------------------------------------
!
!*    3.      Fills output arguments
!             ----------------------
!
KMEGAN_NBR         = NMEGAN_NBR
HMEGAN_NAME(:)     = CMEGAN_NAME(:)
HMEGAN_AREA(:)     = CMEGAN_AREA(:)
HMEGAN_ATYPE(:)    = CMEGAN_ATYPE(:)
HMEGAN_FILE(:)     = CMEGAN_FILE(:)
HMEGAN_FILETYPE(:) = CMEGAN_FILETYPE(:)
IF (LHOOK) CALL DR_HOOK('READ_NAM_PGD_MEGAN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_NAM_PGD_MEGAN
