!     ######spl
      SUBROUTINE AXELOGPRES(PHMIN,PHMAX)
!     ##################################
!
!!****  *AXELOGPRES* - 
!!****    
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
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
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/10/2000
!!      Updated   PM   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODD_RESOLVCAR
USE MODD_PVT
!
IMPLICIT NONE
!
!*       0.1  Dummy arguments and results
!
REAL :: PHMIN,PHMAX
!
!*       0.2  Local variables
!
INTEGER             :: J, JA, ID
!
REAL :: ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
CHARACTER(LEN=5) :: YCAR
!
!-------------------------------------------------------------------------------
!
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
IF(LPRESY)THEN
IF(XPMAX /= 0. .AND. XPMIN /= 0. .AND. XPINT /= 0.)THEN
  IF(XPMIN < 1300.)THEN
    XPMAX=XPMAX*100.
    XPMIN=XPMIN*100.
    XPINT=XPINT*100.
  ENDIF
DO J=INT(XPMIN),INT(XPMAX),-INT((ABS(XPINT)))
  IF(FLOAT(J) >= ANINT(ZWT) .AND. FLOAT(J) <= ANINT(ZWB))THEN
  YCAR=' '
  IF(XPINT > 1000.)THEN
    WRITE(YCAR,'(F5.0)')FLOAT(J)/100.
  ELSE
    WRITE(YCAR,'(F5.0)')FLOAT(J)
  ENDIF
  YCAR=ADJUSTR(YCAR)
  CALL PLCHHQ(ZWL-ZWL/110.,FLOAT(J),YCAR,13.,0.,1.)
  CALL FRSTPT(ZWL,FLOAT(J))
  CALL VECTOR(ZWL+(ZWR-ZWL)/(ZVR-ZVL)*.015,FLOAT(J))
  ENDIF
ENDDO
ELSE
  IF(PHMIN < 1300)THEN
    PHMIN=PHMIN*100
    PHMAX=PHMAX*100
  ENDIF
DO J=INT(PHMIN),INT(PHMAX),-10000
  IF(FLOAT(J) >= ANINT(ZWT) .AND. FLOAT(J) <= ANINT(ZWB))THEN
  YCAR=' '
  IF(PHMAX > 1300.)THEN
    WRITE(YCAR,'(F5.0)')FLOAT(J)/100.
  ELSE
    WRITE(YCAR,'(F5.0)')FLOAT(J)
  ENDIF
  YCAR=ADJUSTR(YCAR)
  print *,' **axelogpres PHMIN,PHMAX ',PHMIN,PHMAX
  print *,' **axelogpres ZWL-ZWL/20.,FLOAT(J),YCAR ',ZWL-ZWL/20.,FLOAT(J),YCAR 
  CALL PLCHHQ(ZWL-ZWL/100.,FLOAT(J),YCAR,13.,0.,1.)
  CALL FRSTPT(ZWL,FLOAT(J))
  CALL VECTOR(ZWL+(ZWR-ZWL)/(ZVR-ZVL)*.015,FLOAT(J))
  ENDIF
ENDDO
ENDIF
ELSE
ENDIF
!*        2.     EXIT
!                ----
!
RETURN
END SUBROUTINE AXELOGPRES
