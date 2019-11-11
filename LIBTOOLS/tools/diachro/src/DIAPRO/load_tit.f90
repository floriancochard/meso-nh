!     ######spl
      MODULE MODI_LOAD_TIT
!     ####################
!
INTERFACE
!
SUBROUTINE LOAD_TIT(HCARIN,KIND)
CHARACTER(LEN=*)  :: HCARIN
INTEGER           :: KIND
END SUBROUTINE LOAD_TIT
!
END INTERFACE
END MODULE MODI_LOAD_TIT
!     ######spl
      SUBROUTINE LOAD_TIT(HCARIN,KIND)
!     ################################
!
!!****  *LOAD_TIT* - 
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
!!      Original       24/11/95
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TIT
USE MODI_RESOLV_TIT

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KIND
!
!*       0.1   Local variables
!              ---------------

INTEGER           :: INDGUIL1, INDGUIL2 ,INDM
INTEGER           :: ILEN
INTEGER           :: J,JM
CHARACTER(LEN=8)  :: YTEM

!
!------------------------------------------------------------------------------
INDM=KIND
IF(HCARIN(KIND:KIND+6) == 'LTITDEF')THEN
  DO J=KIND+8,LEN(HCARIN)
    IF(HCARIN(J:J) /= '=' .AND. HCARIN(J:J) /= '.' &
       .AND. HCARIN(J:J) /= ' ')THEN
    JM=J
    EXIT
    ENDIF
  ENDDO
  IF(HCARIN(JM:JM) == 'T')THEN
    LTITDEF=.TRUE.
    CALL RESOLV_TIT('CTITALL',YTEM)
  ENDIF
  IF(HCARIN(JM:JM) == 'F')LTITDEF=.FALSE.
  RETURN
ENDIF
INDGUIL1=INDEX(HCARIN,'"')
IF(INDGUIL1 == 0)THEN
  INDGUIL1=INDEX(HCARIN,"'")
ENDIF
ILEN=LEN_TRIM(HCARIN)
INDGUIL2=INDEX(HCARIN(INDGUIL1+1:ILEN),'"')
IF(INDGUIL2 == 0)THEN
  INDGUIL2=INDEX(HCARIN(INDGUIL1+1:ILEN),"'")
ENDIF
INDGUIL2=INDGUIL1+INDGUIL2
!print *,' **load_tit INDGUIL1,INDGUIL2 ',INDGUIL1,INDGUIL2

SELECT CASE(HCARIN(INDM:INDM+5))
  CASE('CTITT1')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITT1=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITT2')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITT2=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITT3')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
!   print *,' **load_tit HCARIN et LEN CTITT3 ',HCARIN,LEN(HCARIN),&
!  LEN(CTITT3),CTITT3
    CTITT3=HCARIN(INDGUIL1+1:INDGUIL2-1)
!   print *,' **load_tit HCARIN et LEN CTITT3 ',HCARIN,LEN(HCARIN),&
!  LEN(CTITT3),CTITT3
  CASE('CTITB1')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITB1=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITB2')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITB2=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITB3')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITB3=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITYT')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITYT=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITYM')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITYM=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITYB')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITYB=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITXL')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITXL=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITXM')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITXM=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITXR')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITXR=HCARIN(INDGUIL1+1:INDGUIL2-1)
END SELECT
SELECT CASE(HCARIN(INDM:INDM+7))
  CASE('CTITVAR1')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR1=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITVAR2')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR2=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITVAR3')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR3=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITVAR4')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR4=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITVAR5')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR5=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITVAR6')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR6=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITVAR7')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR7=HCARIN(INDGUIL1+1:INDGUIL2-1)
  CASE('CTITVAR8')
    KIND=999
    IF(INDGUIL1 == 0 .OR. INDGUIL2 == 0)THEN
    print *,' Le TITRE doit etre entre guillemets ou quotes. Corrigez le et '
    print *,' Rentrez le a nouveau '
    ENDIF
    CTITVAR8=HCARIN(INDGUIL1+1:INDGUIL2-1)
END SELECT
RETURN
END SUBROUTINE LOAD_TIT
