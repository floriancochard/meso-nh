!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE PGD_MEGAN(DTCO, UG, U, USS, MSF, HPROGRAM, OCH_BIOEMIS)
!     ##############################################################
!
!!**** *PGD_MEGAN* monitor for averaging and interpolations of physiographic fields
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
!!    Original    15/06/2017
!!    Modified    06/07/2017  J. Pianezze : adapatation to SurfEx v8.0
!!                27/04/2018  J.Escobar : missing USE MODI_GET_SURF_MASK_n
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_MEGAN_SURF_FIELDS_n,ONLY : MEGAN_SURF_FIELDS_t
!
USE MODD_PGD_GRID,           ONLY : NL
USE MODD_PGDWORK,            ONLY : CATYPE
USE MODD_SURF_PAR,           ONLY : XUNDEF
!
USE MODI_GET_LUOUT
USE MODI_PGD_FIELD
USE MODI_READ_NAM_PGD_MEGAN
USE MODI_UNPACK_SAME_RANK
USE MODI_GET_SURF_SIZE_n
USE MODI_GET_SURF_MASK_n
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
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
TYPE(MEGAN_SURF_FIELDS_t), INTENT(INOUT) :: MSF
!
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM     ! Type of program
LOGICAL,             INTENT(OUT)   :: OCH_BIOEMIS     ! emission flag

!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                           :: ILUOUT    ! output listing logical unit
INTEGER                           :: JNBR      ! loop counter on dummy fields
INTEGER                           :: ILU, IL_SEA, IL_LAND, IL
!
!*    0.3    Declaration of namelists
!            ------------------------
!
INTEGER                             :: IMEGAN_NBR
CHARACTER(LEN=20), DIMENSION(1000)  :: YMEGAN_NAME
CHARACTER(LEN=3),  DIMENSION(1000)  :: YMEGAN_AREA
CHARACTER(LEN=3),  DIMENSION(1000)  :: CMEGAN_ATYPE    ! avg type for dummy pgd fields
!                                                      ! 'ARI' , 'INV'
CHARACTER(LEN=28), DIMENSION(1000)  :: CMEGAN_FILE     ! data files
CHARACTER(LEN=6),  DIMENSION(1000)  :: CMEGAN_FILETYPE ! type of these files
REAL, DIMENSION(:), ALLOCATABLE     :: ZMEGAN_FIELD, ZMEGAN_FIELDS
INTEGER, DIMENSION(:), ALLOCATABLE  :: IMASK
CHARACTER(LEN=6)                    :: YMASK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_MEGAN',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
 CALL READ_NAM_PGD_MEGAN(HPROGRAM, IMEGAN_NBR, YMEGAN_NAME, YMEGAN_AREA, &
                        CMEGAN_ATYPE, CMEGAN_FILE, CMEGAN_FILETYPE      )  
!
MSF%NMEGAN_NBR     = IMEGAN_NBR
!
ALLOCATE(MSF%CMEGAN_NAME(MSF%NMEGAN_NBR))
ALLOCATE(MSF%CMEGAN_AREA(MSF%NMEGAN_NBR))
MSF%CMEGAN_NAME(:) = YMEGAN_NAME(1:MSF%NMEGAN_NBR)
MSF%CMEGAN_AREA(:) = YMEGAN_AREA(1:MSF%NMEGAN_NBR)
!
!-------------------------------------------------------------------------------
!
!*    3.      Allocation
!             ----------
!
ALLOCATE(MSF%XMEGAN_FIELDS(NL,MSF%NMEGAN_NBR))
 CALL GET_SURF_SIZE_n(DTCO, U,'LAND', IL_LAND)
 CALL GET_SURF_SIZE_n(DTCO, U,'SEA   ',IL_SEA)
!
ALLOCATE(ZMEGAN_FIELDS (NL))
!
!-------------------------------------------------------------------------------
OCH_BIOEMIS = MSF%NMEGAN_NBR > 0
!-------------------------------------------------------------------------------
!
!
!*    4.      Computations
!             ------------
!
DO JNBR=1,MSF%NMEGAN_NBR

  CATYPE = CMEGAN_ATYPE(JNBR)
  SELECT CASE (MSF%CMEGAN_AREA(JNBR))
    CASE ('LAN')
      IL = IL_LAND
      YMASK='LAND  '
    CASE ('SEA')
      IL = IL_SEA
      YMASK='SEA   '
    CASE ('ALL')
      IL = NL
      YMASK='FULL  '
    CASE DEFAULT
      CALL ABOR1_SFX('PGD_MEGAN (1): MEGAN AREA NOT SUPPORTED')
  END SELECT
  ALLOCATE(ZMEGAN_FIELD (IL))
  ALLOCATE(IMASK(IL))
!
  CALL PGD_FIELD(DTCO, UG, U, USS, &
                 HPROGRAM,MSF%CMEGAN_NAME(JNBR),MSF%CMEGAN_AREA(JNBR),CMEGAN_FILE(JNBR), &
                 CMEGAN_FILETYPE(JNBR),XUNDEF,ZMEGAN_FIELD(:)              )  
  CATYPE = 'ARI'
!
!*    4.2     Expends field on all surface points
  ILU=0
  CALL GET_SURF_MASK_n(DTCO, U, &
                       YMASK,IL,IMASK,ILU,ILUOUT)
  CALL UNPACK_SAME_RANK(IMASK,ZMEGAN_FIELD(:),ZMEGAN_FIELDS(:))
  DEALLOCATE(ZMEGAN_FIELD)
  DEALLOCATE(IMASK)
!
!*    4.3      Weights field on all surface points 
!              (zero weight where field is not defined)
  SELECT CASE (MSF%CMEGAN_AREA(JNBR))
    CASE ('LAN')
      MSF%XMEGAN_FIELDS(:,JNBR) = (U%XNATURE(:)+U%XTOWN(:))*ZMEGAN_FIELDS(:)
    CASE ('SEA')
      MSF%XMEGAN_FIELDS(:,JNBR) = U%XSEA*ZMEGAN_FIELDS(:)
    CASE ('ALL')
      MSF%XMEGAN_FIELDS(:,JNBR) = ZMEGAN_FIELDS(:)
    CASE DEFAULT
      CALL ABOR1_SFX('PGD_MEGAN (2): MEGAN AREA NOT SUPPORTED')
  END SELECT

END DO

DEALLOCATE(ZMEGAN_FIELDS)

IF (LHOOK) CALL DR_HOOK('PGD_MEGAN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_MEGAN
