!     ######spl
      MODULE MODI_READMNMX_FT_PVKT
!     ############################
!
INTERFACE
!
SUBROUTINE READMNMX_FT_PVKT(HCARIN,PMN,PMX)
CHARACTER(LEN=*) :: HCARIN
REAL             :: PMN, PMX
END SUBROUTINE READMNMX_FT_PVKT
!
END INTERFACE
END MODULE MODI_READMNMX_FT_PVKT
!     ######spl
      SUBROUTINE READMNMX_FT_PVKT(HCARIN,PMN,PMX)
!     ###########################################
!
!!****  *READMNMX_FT_PVKT* - 
!!
!!    PURPOSE
!!    -------
!      
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
!!      Module
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       2/09/96
!!      Updated   PM   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
REAL             :: PMN, PMX
!
!*       0.1   Local variables
!              ---------------

INTEGER           :: IMASK
INTEGER           :: J,JM
LOGICAL           :: GOKMN, GOKMX
!REAL,DIMENSION(:),ALLOCATABLE  :: ZFTMN, ZFTMX
!CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE  :: YFTMN, YFTMX
CHARACTER(LEN=LEN(HCARIN)) :: YCARIN, YCARIN2

!
!------------------------------------------------------------------------------
GOKMN=.FALSE.
GOKMX=.FALSE.
YCARIN(1:LEN(YCARIN))=' '
HCARIN=ADJUSTL(HCARIN)
YCARIN=HCARIN
if(nverbia >0)then
  print *,' **READMNMX_FT_PVKT YCARIN ',YCARIN(1:LEN_TRIM(YCARIN))
endif
IMASK=INDEX(YCARIN,'MASK')
IF(IMASK /=0)THEN
DO J=1,LEN(YCARIN)
 IF(YCARIN(J:J) == ' ')THEN
   JM=J-1
   EXIT
 ENDIF
ENDDO
YCARIN(1:LEN(YCARIN))=' '
YCARIN=HCARIN(JM+2:LEN_TRIM(HCARIN))
YCARIN=ADJUSTL(YCARIN)
ENDIF
JM=0
DO J=1,LEN(YCARIN)
 IF(YCARIN(J:J) == ' ')THEN
   JM=J-1
   EXIT
 ENDIF
ENDDO
IF(JM /= 0)THEN
  YCARIN2(1:LEN(YCARIN2))=' '
  YCARIN2=YCARIN(1:JM)
  YCARIN(1:LEN(YCARIN))=' '
  YCARIN=ADJUSTL(YCARIN2)
ENDIF
!

if(nverbia >0)then
  print *,' **READMNMX_FT_PVKT JM,NBFTMN,NBFTMX ',JM,NBFTMN,NBFTMX
endif
IF(NBFTMN == 0)THEN
  GOKMN=.FALSE.
  print *,' AUCUN MIN USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
ELSE
  DO J=1,NBFTMN
!   IF(YCARIN(1:LEN_TRIM(YCARIN)) == CFTMN(J)(1:LEN_TRIM(YCARIN)))THEN
    IF(YCARIN(1:LEN_TRIM(YCARIN)) == CFTMN(J))THEN
      PMN=XFTMN(J)
      print *,' MIN ENREGISTRE SOUS LA FORME XPVMIN_',YCARIN(1:LEN_TRIM(YCARIN)),' UTILISE: ',PMN
      GOKMN=.TRUE.
      EXIT
    ENDIF
  ENDDO
  IF(.NOT.GOKMN)THEN
    print *,' AUCUN MIN USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
  ENDIF
ENDIF
!
IF(NBFTMX == 0)THEN
  GOKMX=.FALSE.
  print *,' AUCUN MAX USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
ELSE
  DO J=1,NBFTMX
!   IF(YCARIN(1:LEN_TRIM(YCARIN)) == CFTMX(J)(1:LEN_TRIM(YCARIN)))THEN
    IF(YCARIN(1:LEN_TRIM(YCARIN)) == CFTMX(J))THEN
      PMX=XFTMX(J)
      print *,' MAX ENREGISTRE SOUS LA FORME XPVMAX_',YCARIN(1:LEN_TRIM(YCARIN)),' UTILISE: ',PMX
      GOKMX=.TRUE.
      EXIT
    ENDIF
  ENDDO
  IF(.NOT.GOKMX)THEN
    print *,' AUCUN MAX USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
  ENDIF
ENDIF
IF(.NOT.GOKMN .OR. .NOT.GOKMX)THEN
  LOK=.FALSE.
  print *,' CALCUL AUTOMATIQUE DES BORNES POUR : ',YCARIN(1:LEN_TRIM(YCARIN))
ELSE
  LOK=.TRUE.
ENDIF
if(nverbia >0)then
  print *,' **READMNMX_FT_PVKT LOK ',LOK
endif
RETURN
END SUBROUTINE READMNMX_FT_PVKT
