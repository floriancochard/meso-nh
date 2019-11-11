!     ######spl
      SUBROUTINE FACTIMP
!     #################
!
!!****  *FACTIMP* -  Impression du facteur a * ou + ou -
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
!!      Module MODD_MEMCV : CDIRCUR
!!
!!      Module MODN_RESOLVCAR
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!     NCAR Graphics Technical documentation, UNIX version 3.2,
!!     Scientific computing division, NCAR/UCAR, Boulder, USA.
!!      Volume 1: Fundamentals, Vers. 1, May 1993
!!      Volume 2: Contouring and mapping tutorial, Vers. 2, May 1993
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/10/99
!!      Updated   PM  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_MEMCV
USE MODD_TYPE_AND_LH
USE MODD_RESOLVCAR
USE MODD_EXPR

IMPLICIT NONE
!
!*       0.1   Local variables
!              ---------------

INTEGER :: J, ILEN, IL
INTEGER :: ID
CHARACTER(LEN=500) :: YCAR200
REAL :: ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
!
!------------------------------------------------------------------------------
!
!*       1.   
!             -----------------------------------------
IL = 0
!
IF (NVERBIA >=5) THEN
   print*, 'FACTIMP ',NSUPERDIA,CFACT(1:NSUPERDIA),LEN_TRIM(CFACT(1))
ENDIF
!
YCAR200(1:LEN(YCAR200))=' '
IF(NSUPERDIA == 1)THEN
  IF(NOPE(1) /= 0)THEN
    CALL PLCHHQ(.99,.032,CDIRCUR(NPARG:NPARD),.010,0.,+1.)
    IF(NMULTDIV(1) /= 0)THEN
      CALL PLCHHQ(.002,.93,CMULTDIV(1),.007,0.,-1.)
    ENDIF
  ELSE
    IF(NMULTDIV(1) /= 0)THEN
      CALL PLCHHQ(.99,.032,CMULTDIV(1),.010,0.,+1.)
    ENDIF
  ENDIF
ELSE
!JD Juillet 2009
  ILEN=0
!JD Juillet 2009        
  DO J = 1,NSUPERDIA
    IF(NOPE(J) /= 0) THEN
      NOPEL=NOPEL+1
      IF(NOPEL == 1)THEN
      IL=LEN_TRIM(CFACT(J))

!JD Juillet 2009
        IF(IL > 0)THEN
!JD Juillet 2009
          YCAR200(1:IL)=CFACT(J)(1:IL)
	  ILEN=LEN_TRIM(YCAR200)
	  ILEN=ILEN+3
!JD Juillet 2009
        ENDIF
!JD Juillet 2009        
      ELSE
	IL=LEN_TRIM(CFACT(J))
!JD Juillet 2009
        IF(IL > 0)THEN
!JD Juillet 2009        
        IF (NVERBIA >=5) THEN
          print*, 'FACTIMP ',J,IL,ILEN,CFACT(J)
        END IF
        YCAR200(ILEN:ILEN-1+IL)=CFACT(J)(1:IL)
	ILEN=LEN_TRIM(YCAR200)
	ILEN=ILEN+3
!JD Juillet 2009
        ENDIF
!JD Juillet 2009        
      ENDIF
      IF(NMULTDIV(J) /= 0)THEN
        ILEN=ILEN-2
        IL=LEN_TRIM(CMULTDIV(J))
        YCAR200(ILEN:ILEN-1+IL)=CMULTDIV(J)(1:IL)
        ILEN=LEN_TRIM(YCAR200)
        ILEN=ILEN+3
      ENDIF
    ENDIF
  ENDDO
!JD Juillet 2009
  IF(ILEN > 3)THEN
!JD Juillet 2009
    ILEN=ILEN-3
!JD Juillet 2009
   ENDIF
!JD Juillet 2009  
  CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
  if(nverbia >0)then
    print *,' ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT ',ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
  endif
!JD Juillet 2009
  IF(ILEN > 0)THEN
!JD Juillet 2009
  IF(CTYPE == 'MASK')THEN
    IF(ILEN > 100)THEN
      YCAR200(97:100)='....'
      CALL PLCHHQ(.002,.93,YCAR200(1:100),.007,0.,-1.)
    ELSE
      CALL PLCHHQ(.002,.93,YCAR200(1:ILEN),.007,0.,-1.)
    ENDIF
  ELSE
    IF(LVARNPVUSER)THEN
        CALL PLCHHQ(.02,.935,YCAR200(1:ILEN),.007,0.,-1.)
    ELSE      
      IF(ILEN > 100)THEN
        J=INDEX(YCAR200(100-IL:100),')')
        YCAR200(100-IL+J+1:100-IL+J+3)='...'
        CALL PLCHHQ(.002,.93,YCAR200(1:100-IL+J+3),.007,0.,-1.)
      ELSE
        CALL PLCHHQ(.002,.93,YCAR200(1:ILEN),.007,0.,-1.)
      ENDIF
    ENDIF
  ENDIF
!JD Juillet 2009
  ENDIF
!JD Juillet 2009  
ENDIF
!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE FACTIMP
