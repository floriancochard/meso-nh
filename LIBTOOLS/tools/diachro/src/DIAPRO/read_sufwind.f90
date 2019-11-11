!     ######spl
      SUBROUTINE READ_SUFWIND(HGROUP)
!     ###############################
!
!!****  *READ_SUFWIND* - 
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
!!      Original       29/01/98
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

CHARACTER(LEN=*) :: HGROUP
!
!*       0.1   Local variables
!              ---------------

!
INTEGER                    ::   J, IND, ILENGP, I
CHARACTER(LEN=LEN(HGROUP)) :: YGROUP
!------------------------------------------------------------------------------
YGROUP=HGROUP
ILENGP=LEN_TRIM(YGROUP)
CSUFWIND='  '
NSUFWIND=0
DO J=1,1
  I=7
  IND=INDEX(YGROUP,'DIRUMVM')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'DIRUTVT')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  I=6
  IND=INDEX(YGROUP,'DDUMVM')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'DDUTVT')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  I=5
  IND=INDEX(YGROUP,'MUMVM')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'MUTVT')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'ULMWM')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'ULTWT')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  I=4
  IND=INDEX(YGROUP,'UMVM')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'UTVT')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  I=3
  IND=INDEX(YGROUP,'ULM')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'ULT')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'VTM')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
  IND=INDEX(YGROUP,'VTT')
  IF(IND /= 0)THEN
    IF(ILENGP == I)THEN
    ELSE IF((ILENGP-I) == 1)THEN
      CSUFWIND(1:1)=YGROUP(IND+I:IND+I)
      NSUFWIND=1
    ELSE IF((ILENGP-I) == 2)THEN
      CSUFWIND(1:2)=YGROUP(IND+I:IND+I+1)
      NSUFWIND=2
    ENDIF
    EXIT
  ENDIF
ENDDO
!print *,' YGROUP CSUFWIND NSUFWIND ',YGROUP,CSUFWIND,NSUFWIND
!

!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE READ_SUFWIND
