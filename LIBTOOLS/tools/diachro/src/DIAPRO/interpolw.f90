!     ######spl
      SUBROUTINE INTERPOLW(PZZU, PZZW, PSTRU, PSTRW)
!     ####################
!
!!****  *INTERPOLW* - Defines the display window for a cartesian model
!!
!!    PURPOSE
!!    -------
!       Interpolation des composantes du vent pour les streamlines en CV
!
!
!!**  METHOD
!!    ------
!!
!
!!
!!    EXTERNAL
!!    --------
!!      SET      : defines NCAR window and viewport in normalized and user
!!                 coordinates
!!      LABMOD   : defines axis label format
!!      GRIDAL   : draws axis divisions and ticks
!!      PERIM    : draws a perimeter box for the current plot
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
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
!!      Original       10/04/02
!!      Updated   PM
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
USE MODD_DIM1
USE MODD_GRID1
USE MODN_NCAR
!
IMPLICIT NONE
!
!
!*       0.1   Commons
!
COMMON/LOGI/LVERT,LHOR,LPT,LXABS
COMMON/TEMV/XZWORKZ,XZZDS,NINX,NINY
#include "big.h"
REAL,DIMENSION(N2DVERTX,2500):: XZWORKZ
!REAL,DIMENSION(1000,400):: XZWORKZ
REAL,DIMENSION(N2DVERTX):: XZZDS
!REAL,DIMENSION(1000):: XZZDS
INTEGER :: NINX, NINY
LOGICAL :: LVERT, LHOR, LPT, LXABS
!
!*       0.2   Dummy arguments and results
!
REAL,DIMENSION(:,:) :: PZZU, PZZW, PSTRU, PSTRW
!
!*       0.3   Local variables
!

REAL :: ZZ, ZPASZ, ZR, ZT, ZMX
!REAL,DIMENSION(:),ALLOCATABLE,SAVE  :: ZW
INTEGER :: I, J, K, IPASZ
INTEGER :: ISZ, ITER,ID,IE
!
!-------------------------------------------------------------------------------
IPASZ=NZSTR
ZMX=0.
DO K=1,NINY
DO I=1,NINX
  IF(XZWORKZ(I,K) /= XSPVAL)ZMX=MAX(ZMX,XZWORKZ(I,K))
ENDDO
ENDDO
ZPASZ=ZMX /(IPASZ-1)
if (nverbia >0)then
print *,' IPASZ ZPASZ MAXVAL(XZWORKZ) ',IPASZ,ZPASZ,ZMX
endif
IF(ALLOCATED(XZSTR))DEALLOCATE(XZSTR)
ALLOCATE(XZSTR(IPASZ))
XZSTR(1)=0.
DO J=2,IPASZ
XZSTR(J)=XZSTR(J-1)+ZPASZ
ENDDO
if (nverbia >0)then
print *,' **interpolw IPASZ XZSTR ',IPASZ
print *,XZSTR
endif
!!!!!PROVI
!IF(IPASZ == 100)THEN
! J=NINY-1
! ZT=XZWORKZ(20,J)-XZWORKZ(20,2)
! print *,' I=20        XZWORKZ(20,I)             DIFF              rap '
! print 102
! DO I=J,2,-1
! print 103,I,XZWORKZ(20,I),(XZWORKZ(20,I)-XZWORKZ(20,I-1)),(XZWORKZ(20,I)-XZWORKZ(20,I-1))/ZT,(XZWORKZ(1,I)-XZWORKZ(1,I-1))
! 103 format(1X,I3,4(E15.8,5X))
! ENDDO
!ENDIF
!!!!!PROVI

PSTRU=XSPVAL
PSTRW=XSPVAL
DO J=1,IPASZ
  ZZ=XZSTR(J)
DO I=1,NINX
DO K=1,NINY-1
  IF(ZZ < XZWORKZ(I,2) .OR. ZZ > XZWORKZ(I,NINY-1))THEN
    EXIT
  ELSEIF(ZZ == XZWORKZ(I,K))THEN
    IF(PZZU(I,K) == XSPVAL .OR. (PZZW(I,K) == XSPVAL))THEN
    EXIT
    ELSE
    IF(PZZU(I,K) /= XSPVAL)THEN
    PSTRU(I,J)=PZZU(I,K) 
    ENDIF
    IF(PZZW(I,K) /= XSPVAL)THEN
    PSTRW(I,J)=PZZW(I,K) 
    ENDIF
    EXIT
    ENDIF
  ELSEIF(ZZ > XZWORKZ(I,K) .AND. ZZ < XZWORKZ(I,K+1))THEN
    IF(XZWORKZ(I,K+1)-XZWORKZ(I,K) /= 0.)THEN
      IF(PZZU(I,K) == XSPVAL .OR. PZZW(I,K) == XSPVAL .OR. &
      PZZU(I,K+1) == XSPVAL .OR. PZZW(I,K+1) == XSPVAL )THEN
if (nverbia >0)then
        print *,'**interpolw I K PZZU(I,K),PZZU(I,K+1),PZZW(I,K), PZZW(I,K+1) ',&
        I,K,PZZU(I,K),PZZU(I,K+1),PZZW(I,K), PZZW(I,K+1)
endif
        EXIT
      ELSE
 
        ZR=(ZZ-XZWORKZ(I,K))/(XZWORKZ(I,K+1)-XZWORKZ(I,K))
        PSTRU(I,J)=PZZU(I,K) + ZR*(PZZU(I,K+1)-PZZU(I,K))
        PSTRW(I,J)=PZZW(I,K) + ZR*(PZZW(I,K+1)-PZZW(I,K))
        EXIT
      ENDIF
    ELSE
      IF(PZZU(I,K) == XSPVAL .OR. (PZZW(I,K) == XSPVAL))THEN
        EXIT
      ELSE
        IF(PZZU(I,K) /= XSPVAL)THEN
          PSTRU(I,J)=PZZU(I,K) 
        ENDIF
        IF(PZZW(I,K) /= XSPVAL)THEN
          PSTRW(I,J)=PZZW(I,K) 
        ENDIF
        EXIT
      ENDIF
    ENDIF
  ENDIF
ENDDO
ENDDO
ENDDO
if (nverbia >0)then
print *,' **interpolw sortie PSTRU,PSTRW '
ISZ=SIZE(PSTRU,1)
ITER=ISZ/5
IF(ITER*5 > ISZ)ITER=ITER+1
DO I=1,ITER
ID=(I-1)*5 +1
IE=ID+4
print 101,ID,ID+1,ID+2,ID+3,ID+4
print 102
DO J=IPASZ,1,-1
print 100,J,PSTRW(ID:IE,J),XZSTR(J)
!print 100,J,PSTRU(ID:IE,J),XZSTR(J)
ENDDO
print 102
ENDDO
endif
100 FORMAT(I3,5E13.6,E12.5)
101 FORMAT(8X,I3,4(10X,I3),10X,'XZSTR')
102 FORMAT(78('*'))
!-----------------------------------------------------------------------------
!
!*      2.   EXIT
!            ----
!
RETURN
END SUBROUTINE  INTERPOLW
