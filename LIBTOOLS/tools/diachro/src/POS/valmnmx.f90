!     ######spl
      SUBROUTINE VALMNMX(PMIN,PMAX)
!     #############################
!
!!****  *VALMNMX* - Dans le cadre des profils, determination automatique
!                   des bornes min et max. 
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      None
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       14/03/95
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1  Dummy arguments
!          

REAL        :: PMIN,PMAX
!
!*       0.2  local variables
!          
REAL        :: ZMN,ZMX,Z
REAL        :: ZVAL, ZABSVAL,ZJJB, ZJJT
REAL,DIMENSION(38)  :: ZDIXPM
INTEGER     :: J, JJ, ISIGNVAL, IJJM
!
!-------------------------------------------------------------------------------
ZDIXPM(1)=1.E-1;ZDIXPM(2)=1.E-2;ZDIXPM(3)=1.E-3;ZDIXPM(4)=1.E-4;ZDIXPM(5)=1.E-5
ZDIXPM(6)=1.E-6;ZDIXPM(7)=1.E-7;ZDIXPM(8)=1.E-8;ZDIXPM(9)=1.E-9
ZDIXPM(10)=1.E-10;ZDIXPM(11)=1.E-11;ZDIXPM(12)=1.E-12;ZDIXPM(13)=1.E-13
ZDIXPM(14)=1.E-14;ZDIXPM(15)=1.E-15;ZDIXPM(16)=1.E-16
ZDIXPM(17)=1.E-17;ZDIXPM(18)=1.E-18;ZDIXPM(19)=1.E-19
ZDIXPM(20)=1.E-20;ZDIXPM(21)=1.E-21;ZDIXPM(22)=1.E-22
ZDIXPM(23)=1.E-23;ZDIXPM(24)=1.E-24;ZDIXPM(25)=1.E-25
ZDIXPM(26)=1.E-26;ZDIXPM(27)=1.E-27;ZDIXPM(28)=1.E-28
ZDIXPM(29)=1.E-29;ZDIXPM(30)=1.E-30;ZDIXPM(31)=1.E-31
ZDIXPM(32)=1.E-32;ZDIXPM(33)=1.E-33;ZDIXPM(34)=1.E-34
ZDIXPM(35)=1.E-35;ZDIXPM(36)=1.E-36;ZDIXPM(37)=1.E-37
ZDIXPM(38)=1.E-38

! Juillet 99 pour correction sur station du resultat de la fonction ANINT
! pour les valeurs > a 2**31-1
Z=HUGE(1)

DO J=1,2
  IF(J == 1)ZVAL=PMIN
  IF(J == 2)ZVAL=PMAX
  ISIGNVAL=SIGN(1.,ZVAL)
  ZABSVAL=ABS(ZVAL)
! Rectification en Juin 99 pour tenir compte de la capacite des entiers
! sur station
! Juillet 99 pour correction sur station du resultat de la fonction ANINT
! pour les valeurs > a 2**31-1
    IF(ZABSVAL >= Z )THEN
      SELECT CASE(ISIGNVAL)
        CASE(1)
          IF(J == 1)ZMN=AINT(ZABSVAL-1.)
          IF(J == 2)ZMX=AINT(ZABSVAL+1.)
        CASE(-1)
          IF(J == 1)ZMN=AINT(ZABSVAL+1.)
          IF(J == 2)ZMX=AINT(ZABSVAL-1.)
      END SELECT
    ELSE IF(ZABSVAL >= 1. .AND. ZABSVAL < Z)THEN
!   IF(ZABSVAL >= 1.)THEN
      SELECT CASE(ISIGNVAL)
        CASE(1)
          IF(J == 1)ZMN=ANINT(ZABSVAL-1.)
          IF(J == 2)ZMX=ANINT(ZABSVAL+1.)
        CASE(-1)
          IF(J == 1)ZMN=ANINT(ZABSVAL+1.)
          IF(J == 2)ZMX=ANINT(ZABSVAL-1.)
      END SELECT
    ELSE IF(ZABSVAL >=1.E-37 .AND. ZABSVAL <1.)THEN
      SELECT CASE(ISIGNVAL)
        CASE(1)
        IF(ZABSVAL >= ZDIXPM(1) .AND. ZABSVAL < 1.)THEN
          DO JJ=1,9
          ZJJT=(JJ+1)*.1
          ZJJB=JJ*.1
          IF(ZABSVAL >= ZJJB .AND. ZABSVAL < ZJJT)EXIT
          ENDDO
          IF(J == 1)ZMN=ZJJB
          IF(J == 2)ZMX=ZJJT
        ELSE
          DO JJ=1,37
            IF(ZABSVAL >= ZDIXPM(JJ+1) .AND. ZABSVAL < ZDIXPM(JJ))EXIT
          ENDDO
          IJJM=JJ
          IF(J == 1)ZMN=ZDIXPM(IJJM+1)
          IF(J == 2)ZMX=ZDIXPM(IJJM)
        ENDIF
        CASE(-1)
        IF(ZABSVAL >= ZDIXPM(1) .AND. ZABSVAL < 1.)THEN
          DO JJ=1,9
          ZJJT=(JJ+1)*.1
          ZJJB=JJ*.1
          IF(ZABSVAL >= ZJJB .AND. ZABSVAL < ZJJT)EXIT
          ENDDO
          IF(J == 1)ZMN=ZJJT
          IF(J == 2)ZMX=ZJJB
        ELSE
          DO JJ=1,37
            IF(ZABSVAL >= ZDIXPM(JJ+1) .AND. ZABSVAL < ZDIXPM(JJ))EXIT
          ENDDO
          IJJM=JJ
          IF(J == 1)ZMN=ZDIXPM(IJJM)
          IF(J == 2)ZMX=ZDIXPM(IJJM+1)
        ENDIF
      END SELECT
    ELSE
      IF(J == 1)ZMN=0.
      IF(J == 2)ZMX=0.
    END IF
IF(J == 1)ZMN=ZMN*ISIGNVAL
IF(J == 1)PMIN=ZMN
IF(J == 2)ZMX=ZMX*ISIGNVAL
IF(J == 2)PMAX=ZMX
ENDDO
RETURN
END SUBROUTINE VALMNMX

