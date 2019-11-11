!     ######spl
      MODULE  MODI_READXISOLEVP
!     #########################
!
INTERFACE
!
SUBROUTINE READXISOLEVP(HCARIN,K,PISOLEVP)
INTEGER          :: K
CHARACTER(LEN=*) :: HCARIN
REAL,DIMENSION(:):: PISOLEVP
END SUBROUTINE READXISOLEVP
!
END INTERFACE
!
END MODULE MODI_READXISOLEVP
!     ######spl
      SUBROUTINE READXISOLEVP(HCARIN,K,PISOLEVP)
!     ##########################################
!
!!****  *READXISOLEVP* - 
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

INTEGER           :: K
CHARACTER(LEN=*)  :: HCARIN
REAL,DIMENSION(:) :: PISOLEVP
!
!*       0.1   Local variables
!              ---------------

INTEGER           :: IMASK
INTEGER           :: J,JM
CHARACTER(LEN=LEN(HCARIN)) :: YCARIN, YCARIN2

!
!------------------------------------------------------------------------------
YCARIN(1:LEN(YCARIN))=' '
HCARIN=ADJUSTL(HCARIN)
YCARIN=HCARIN
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
LISOLEVP=.FALSE.
IF(NBISOLEVP == 0)THEN
  LISOLEVP=.FALSE.
  print *,' AUCUNE VALEUR USER ENREGISTREE POUR : ',YCARIN(1:LEN_TRIM(YCARIN))&
  ,' sous la forme XISOLEV_PROC= '
ELSE
  DO J=1,NBISOLEVP
    IF(YCARIN(1:LEN_TRIM(YCARIN)) == CISOLEVP(J)(1:LEN_TRIM(CISOLEVP(J))))THEN
      K=NLENP(J)
      PISOLEVP(1:NLENP(J))=XISOLEVP(1:NLENP(J),J)
      LISOLEVP=.TRUE.
      IF(NVERBIA >= 5)THEN
        print *,' READXISOLEVP NLENP PISOLEVP ',K,PISOLEVP(1:NLENP(J))
      ENDIF
      EXIT
    ENDIF
  ENDDO
  IF(.NOT.LISOLEVP)THEN
    print *,' AUCUNE VALEUR USER ENREGISTREE POUR : ',YCARIN(1:LEN_TRIM(YCARIN))&
    ,' sous la forme XISOLEV_PROC= '
  ELSE
     print *,' UTILISATION DES VALEURS ENREGISTREES sous la forme XISOLEV_PROC= '
     print *,' POUR : ',YCARIN(1:LEN_TRIM(YCARIN))
     print *,PISOLEVP(1:K-1)
  ENDIF
ENDIF
!
IF(.NOT.LISOLEVP)THEN
  print *,' UTILISATION DES VALEURS DE XISOLEV= (si elles existent) POUR : ',YCARIN(1:LEN_TRIM(YCARIN))
  DO J=1,SIZE(XISOLEV,1)
    IF(XISOLEV(J) == 9999.)THEN
      print *,XISOLEV(1:J-1)
      EXIT
    ENDIF
  ENDDO
ENDIF
RETURN
END SUBROUTINE READXISOLEVP
