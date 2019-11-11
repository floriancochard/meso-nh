!     ######spl
      MODULE MODI_RESOLVI
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVI(HCARIN,KI,KOUT)
CHARACTER(LEN=*)  :: HCARIN
INTEGER           :: KI, KOUT
END SUBROUTINE RESOLVI
!
END INTERFACE
END MODULE MODI_RESOLVI
!     ##################################
      SUBROUTINE RESOLVI(HCARIN,KI,KOUT)
!     ##################################
!
!!****  *RESOLVI* - 
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
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KI, KOUT
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=8) :: YC8
INTEGER          :: ILENC
INTEGER          :: J,JM, I

!
!------------------------------------------------------------------------------
ILENC=LEN_TRIM(HCARIN)

DO J=KI,ILENC
  IF(HCARIN(J:J) == '=')EXIT
ENDDO

JM=J+1
YC8='        '
I=0

DO J=JM,ILENC
  IF(HCARIN(J:J) == '0'.OR.HCARIN(J:J) == '1'.OR.HCARIN(J:J) == '2'  &
    .OR.HCARIN(J:J) == '3'.OR.HCARIN(J:J) == '4'.OR.HCARIN(J:J) == '5' &
    .OR.HCARIN(J:J) == '6'.OR.HCARIN(J:J) == '7'.OR.HCARIN(J:J) == '8' &
    .OR.HCARIN(J:J) == '9')THEN
    YC8(1:1)=HCARIN(J:J)
    I=1
    IF(J+I > ILENC)EXIT
    IF(HCARIN(J+1:J+1) /= '0' .AND. HCARIN(J+1:J+1) /= '1' .AND.  &
       HCARIN(J+1:J+1) /= '2' .AND. HCARIN(J+1:J+1) /= '3' .AND.  &
       HCARIN(J+1:J+1) /= '4' .AND. HCARIN(J+1:J+1) /= '5' .AND.  &
       HCARIN(J+1:J+1) /= '6' .AND. HCARIN(J+1:J+1) /= '7' .AND.  &
       HCARIN(J+1:J+1) /= '8' .AND. HCARIN(J+1:J+1) /= '9')THEN
       EXIT
    ELSE
      YC8(2:2)=HCARIN(J+1:J+1)
      I=2
      IF(J+I > ILENC)EXIT
      IF(HCARIN(J+2:J+2) /= '0' .AND. HCARIN(J+2:J+2) /= '1' .AND.  &
	 HCARIN(J+2:J+2) /= '2' .AND. HCARIN(J+2:J+2) /= '3' .AND.  &
	 HCARIN(J+2:J+2) /= '4' .AND. HCARIN(J+2:J+2) /= '5' .AND.  &
	 HCARIN(J+2:J+2) /= '6' .AND. HCARIN(J+2:J+2) /= '7' .AND.  &
	 HCARIN(J+2:J+2) /= '8' .AND. HCARIN(J+2:J+2) /= '9')THEN
	 EXIT
      ELSE
	YC8(3:3)=HCARIN(J+2:J+2)
	I=3
        IF(J+I > ILENC)EXIT
	IF(HCARIN(J+3:J+3) /= '0' .AND. HCARIN(J+3:J+3) /= '1' .AND.  &
	   HCARIN(J+3:J+3) /= '2' .AND. HCARIN(J+3:J+3) /= '3' .AND.  &
	   HCARIN(J+3:J+3) /= '4' .AND. HCARIN(J+3:J+3) /= '5' .AND.  &
	   HCARIN(J+3:J+3) /= '6' .AND. HCARIN(J+3:J+3) /= '7' .AND.  &
	   HCARIN(J+3:J+3) /= '8' .AND. HCARIN(J+3:J+3) /= '9')THEN
	   EXIT
        ELSE
  	  YC8(4:4)=HCARIN(J+3:J+3)
  	  I=4
          IF(J+I > ILENC)EXIT
  	  IF(HCARIN(J+4:J+4) /= '0' .AND. HCARIN(J+4:J+4) /= '1' .AND.  &
  	     HCARIN(J+4:J+4) /= '2' .AND. HCARIN(J+4:J+4) /= '3' .AND.  &
  	     HCARIN(J+4:J+4) /= '4' .AND. HCARIN(J+4:J+4) /= '5' .AND.  &
  	     HCARIN(J+4:J+4) /= '6' .AND. HCARIN(J+4:J+4) /= '7' .AND.  &
  	     HCARIN(J+4:J+4) /= '8' .AND. HCARIN(J+4:J+4) /= '9')THEN
  	     EXIT
          ELSE
  	    YC8(5:5)=HCARIN(J+4:J+4)
  	    I=5
            IF(J+I > ILENC)EXIT
  	    IF(HCARIN(J+5:J+5) /= '0' .AND. HCARIN(J+5:J+5) /= '1' .AND.  &
  	       HCARIN(J+5:J+5) /= '2' .AND. HCARIN(J+5:J+5) /= '3' .AND.  &
  	       HCARIN(J+5:J+5) /= '4' .AND. HCARIN(J+5:J+5) /= '5' .AND.  &
  	       HCARIN(J+5:J+5) /= '6' .AND. HCARIN(J+5:J+5) /= '7' .AND.  &
  	       HCARIN(J+5:J+5) /= '8' .AND. HCARIN(J+5:J+5) /= '9')THEN
  	       EXIT
            ELSE
  	      YC8(6:6)=HCARIN(J+5:J+5)
  	      I=6
              IF(J+I > ILENC)EXIT
  	      IF(HCARIN(J+6:J+6) /= '0' .AND. HCARIN(J+6:J+6) /= '1' .AND.  &
  	         HCARIN(J+6:J+6) /= '2' .AND. HCARIN(J+6:J+6) /= '3' .AND.  &
  	         HCARIN(J+6:J+6) /= '4' .AND. HCARIN(J+6:J+6) /= '5' .AND.  &
  	         HCARIN(J+6:J+6) /= '6' .AND. HCARIN(J+6:J+6) /= '7' .AND.  &
  	         HCARIN(J+6:J+6) /= '8' .AND. HCARIN(J+6:J+6) /= '9')THEN
  	         EXIT
              ELSE
  	        YC8(7:7)=HCARIN(J+6:J+6)
  	        I=7
                IF(J+I > ILENC)EXIT
  	        IF(HCARIN(J+7:J+7) /= '0' .AND. HCARIN(J+7:J+7) /= '1' .AND.  &
  	           HCARIN(J+7:J+7) /= '2' .AND. HCARIN(J+7:J+7) /= '3' .AND.  &
  	           HCARIN(J+7:J+7) /= '4' .AND. HCARIN(J+7:J+7) /= '5' .AND.  &
  	           HCARIN(J+7:J+7) /= '6' .AND. HCARIN(J+7:J+7) /= '7' .AND.  &
  	           HCARIN(J+7:J+7) /= '8' .AND. HCARIN(J+7:J+7) /= '9')THEN
  	           EXIT
                ELSE
  	          YC8(8:8)=HCARIN(J+7:J+7)
  	          I=8
                  IF(J+I > ILENC)EXIT
  	          IF(HCARIN(J+8:J+8) /= '0' .AND. HCARIN(J+8:J+8) /= '1' .AND. &
  	             HCARIN(J+8:J+8) /= '2' .AND. HCARIN(J+8:J+8) /= '3' .AND. &
  	             HCARIN(J+8:J+8) /= '4' .AND. HCARIN(J+8:J+8) /= '5' .AND. &
  	             HCARIN(J+8:J+8) /= '6' .AND. HCARIN(J+8:J+8) /= '8' .AND. &
  	             HCARIN(J+8:J+8) /= '8' .AND. HCARIN(J+8:J+8) /= '9')THEN
  	             EXIT
                  ELSE
	            print *,' PB AVEC LA VALEUR FOURNIE  ', &
	            HCARIN(J-1:J+9),' VERIFIEZ LA ET RENTREZ LA A NOUVEAU ', &
		    '(8 chiffres MAXIMUM)'
                    KOUT=999999999
		    RETURN
	          ENDIF
	        ENDIF
	      ENDIF
	    ENDIF
	  ENDIF
	ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDDO

IF(I == 0)THEN
print *,' ABSENCE DE VALEUR. VERIFIEZ ET RENTREZ LA A NOUVEAU '
KOUT=999999999
RETURN
ENDIF
READ(YC8(1:I),*)KOUT
IF(HCARIN(J-1:J-1) == '-')KOUT=KOUT*(-1)
    
RETURN
END SUBROUTINE RESOLVI  
!     ######spl
      MODULE MODI_RESOLVIARRAY
!     ########################
!
INTERFACE
!
SUBROUTINE RESOLVIARRAY(HCARIN,KIND,KOUT,KIARRAY)
CHARACTER(LEN=*)  :: HCARIN
INTEGER          :: KIND, KIARRAY
INTEGER,DIMENSION(:)             :: KOUT
END SUBROUTINE RESOLVIARRAY
!
END INTERFACE
END MODULE MODI_RESOLVIARRAY
!     #################################################
      SUBROUTINE RESOLVIARRAY(HCARIN,KIND,KOUT,KIARRAY)
!     #################################################
!
!!****  *RESOLVIARRAY* - 
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
USE MODD_RESOLVCAR
USE MODN_PARA

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KIND, KIARRAY
INTEGER,DIMENSION(:)             :: KOUT
!
!*       0.1   Local variables
!              ---------------

INTEGER           :: ILENC
INTEGER           :: J,JM, JMF
INTEGER           :: INBV, IND9999

!
!------------------------------------------------------------------------------
ILENC=LEN_TRIM(HCARIN)
KOUT=9999

DO J=KIND,ILENC
  IF(HCARIN(J:J) == '=')EXIT
ENDDO

JM=J+1
DO J=1,10
  IF(HCARIN(JM:JM) == ' ')THEN
    JM=JM+1
  ELSE
    EXIT
  ENDIF
ENDDO

IND9999=INDEX(HCARIN(JM:ILENC),'9999.')
IF(IND9999 == 0)THEN
  IND9999=INDEX(HCARIN(JM:ILENC),'9999')
ENDIF
IF(IND9999 == 0)THEN
  JMF=ILENC
ELSE
  JMF=IND9999+JM-1+3
ENDIF
INBV=0
DO J=JM,JMF
  IF(HCARIN(J:J) == ',')THEN
    INBV=INBV+1
  ENDIF
ENDDO

IF(IND9999 == 0)THEN
  INBV=INBV+1
ENDIF
READ(HCARIN(JM:JMF),*)(KOUT(J),J=1,INBV)
KIARRAY=INBV
IF(NVERBIA >= 5)THEN
  print *,' RESOLVIARRAY ',INBV,(KOUT(J),J=1,INBV)
ENDIF
RETURN
END SUBROUTINE RESOLVIARRAY
!     ######spl
      MODULE MODI_RESOLVK
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVK(HCARIN,KINDK,KJ)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDK, KJ
END SUBROUTINE RESOLVK
!
END INTERFACE
!
END MODULE MODI_RESOLVK
!     ###################################
      SUBROUTINE RESOLVK(HCARIN,KINDK,KJ)
!     ###################################
!
!!****  *RESOLVK* - 
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
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDK, KJ
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=80) :: YCART
CHARACTER(LEN=20) :: YCAR
INTEGER          :: ILENC, ILENCART
INTEGER          :: INDKF, INDTO, INDBY, INDV, INDVM
INTEGER          :: ICAS, J

!
!------------------------------------------------------------------------------
INDKF = 0
INDTO = 0
INDBY = 0
INDV  = 0
ICAS = 0

NBLVLKDIA(KJ,:)=0
NLVLKDIA(:,KJ,:)=0
LVLKDIALL(KJ,:)=.FALSE.

IF(KINDK == 0)THEN
  LVLKDIALL(KJ,:) = .TRUE.
  RETURN
END IF

ILENC = LEN(HCARIN)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INDTO = INDEX(HCARIN(KINDK+3:ILENC),'_TO_')
  INDBY = INDEX(HCARIN(KINDK+3:ILENC),'_BY_')
  INDKF = INDEX(HCARIN(KINDK+3:ILENC),'_')
  IF(INDTO /= 0)THEN
  IF(INDKF < INDTO)THEN
!
! ICAS = 1  Niveau K unique ou plusieurs separes par des virgules
!
    INDTO=0;INDBY=0
    ICAS = 1
  ELSE IF(INDKF == INDTO)THEN
!
! ICAS = 3  Niv1 _TO_ Nivn _BY_ Nivx
!
    IF(INDBY /= 0)THEN
      DO J=INDTO+4+KINDK+3,INDBY+KINDK+3
        IF(HCARIN(J:J) == '_')THEN
          IF(HCARIN(J:J+3) == '_BY_')THEN
            EXIT
          ELSE
            INDBY=0
            EXIT
          END IF
        END IF
      ENDDO
    END IF
    IF(INDBY /= 0)THEN
      INDKF=INDEX(HCARIN(KINDK+3+INDBY+4:ILENC),'_')
      IF(INDKF /= 0)INDKF=INDKF+INDBY+4
      ICAS = 3
      LKINCRDIA(KJ,:) = .TRUE.
    ELSE
!
! ICAS = 2  Niv1 _TO_ Nivn
!
      INDKF=INDEX(HCARIN(KINDK+3+INDTO+4:ILENC),'_')
      IF(INDKF /= 0)INDKF=INDKF+INDTO+4
      ICAS = 2
      LKINCRDIA(KJ,:) = .TRUE.
    END IF
  END IF
  ELSE
    ICAS = 1
  END IF
IF(INDKF == 0)THEN
  INDKF = ILENC
ELSE
  INDKF = INDKF+KINDK+3-1-1
END IF


YCART(1:LEN(YCART))=' '
YCAR(1:LEN(YCAR))=' '
!
! Extraction de la partie Niveaux K dans YCART(1:ILENCART)
!
!print *,' KINDK INDKF ',KINDK,INDKF
YCART = ADJUSTL(HCARIN(KINDK+3:INDKF))
ILENCART = LEN_TRIM(YCART)
!print *,' YCART ',ILENCART,' ',YCART

! Recherche a nouveau des chaines de car. _TO_ , _BY_ et d'une virgule
! par rapport au debut de YCART

INDTO = INDEX(YCART,'_TO_')
INDBY = INDEX(YCART,'_BY_')
INDV = INDEX(YCART(1:ILENCART),',')
IF(ICAS == 1 .AND. INDV == 0)ICAS=0
!
! Expression des Niveaux K par mots-cles (LVLKALL ou LVLK1....)
!
IF(YCART(1:7) == 'LVLKALL')THEN
  LVLKDIALL(KJ,:) = .TRUE.
  if(nverbia >0)then
  print *,' RESOLVK LVLKALL '
  print *,' LVLKDIALL ',LVLKDIALL(KJ,1)
  print *,' NBLVLKDIA ',NBLVLKDIA(KJ,1)
  print *,' NLVLKDIA ',(NLVLKDIA(J,KJ,1),J=1,NBLVLKDIA(KJ,1))
  endif
  RETURN

ELSE IF(YCART(1:4) == 'LVLK')THEN
!print *,' YCART(1:4) ',YCART(1:4),' ICAS ',ICAS

  NBLVLKDIA(KJ,:)=NBLVLKDIA(KJ,:)+1
  SELECT CASE(ICAS)
    CASE(1)
!print *,' INDV YCART(5:5) ',INDV,YCART(5:5)
      IF(INDV-4-1 == 1)READ(YCART(5:5),'(I1)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      IF(INDV-4-1 == 2)READ(YCART(5:6),'(I2)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      DO J = 1,100
        INDVM=INDV
        INDV=0
        INDV=INDEX(YCART(INDVM+1:ILENCART),',')
        IF(INDV == 0)THEN
          NBLVLKDIA(KJ,:)=NBLVLKDIA(KJ,:)+1
          IF(ILENCART-(INDVM+4) == 1)READ(YCART(INDVM+4+1:INDVM+4+1),'(I1)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
          IF(ILENCART-(INDVM+4) == 2)READ(YCART(INDVM+4+1:ILENCART),'(I2)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
          NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
          EXIT
        ELSE
          INDV=INDV+INDVM
          NBLVLKDIA(KJ,:)=NBLVLKDIA(KJ,:)+1
          IF(INDV-(INDVM+4)-1 == 1)READ(YCART(INDVM+4+1:INDVM+4+1),'(I1)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
          IF(INDV-(INDVM+4)-1 == 2)READ(YCART(INDVM+4+1:INDV-1),'(I2)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
          NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
        END IF
      ENDDO   
      
    CASE(2)
      IF(INDTO-4-1 == 1)READ(YCART(5:5),'(I1)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      IF(INDTO-4-1 == 2)READ(YCART(5:6),'(I2)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      NBLVLKDIA(KJ,:)=NBLVLKDIA(KJ,:)+1
      IF(ILENCART-(INDTO+3+4) == 1)READ(YCART(INDTO+3+4+1:INDTO+3+4+1),'(I1)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      IF(ILENCART-(INDTO+3+4) == 2)READ(YCART(INDTO+3+4+1:ILENCART),'(I2)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
! 1 seul temps
    CASE DEFAULT
      IF(ILENCART-4 == 1)READ(YCART(5:5),'(I1)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      IF(ILENCART-4 == 2)READ(YCART(5:6),'(I2)')NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)

  END SELECT
  if(nverbia >0)then
  print *,' RESOLVK ICAS '
  print *,' LVLKDIALL ',LVLKDIALL(KJ,1)
  print *,' NBLVLKDIA ',NBLVLKDIA(KJ,1)
  print *,' NLVLKDIA ',(NLVLKDIA(J,KJ,1),J=1,NBLVLKDIA(KJ,1))
  endif
  RETURN
ELSE

!
! Expression des Niveaux K en numerique
!
  IF(INDV == 0)THEN

! Cas  _TO_  _BY_

    IF(INDTO /= 0)THEN
      YCAR = ADJUSTL(YCART(1:INDTO-1))
      NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
      CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
      NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      IF(INDBY /= 0)THEN
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:INDBY-1))
        NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
        NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDBY+4:ILENCART))
        NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
        NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      ELSE
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:ILENCART))
        NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
        NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      END IF
    ELSE

! Cas un seul niveau en fin de chaine de car. HCARIN ou au milieu

      IF(ILENCART > 9)THEN
	print *,' PB ecriture temps '
	STOP
      ELSE
	YCAR = ADJUSTL(YCART(1:ILENCART))
	NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
        NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      END IF

    END IF

  ELSE

! Presence de virgules

    YCAR = ADJUSTL(YCART(1:INDV-1))
    NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
    CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
    NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
    DO J = 1,100
      INDVM=INDV
      INDV = 0
      YCAR(1:LEN(YCAR))=' '
      INDV = INDEX(YCART(INDVM+1:ILENCART),',')
!     print *,' INDV ',INDV
      IF(INDV == 0)THEN
	YCAR = ADJUSTL(YCART(INDVM+1:ILENCART))
	NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
        NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
	EXIT
      ELSE
        INDV=INDV+INDVM
	YCAR = ADJUSTL(YCART(INDVM+1:INDV-1))
	NBLVLKDIA(KJ,:) = NBLVLKDIA(KJ,:)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1))
        NLVLKDIA(NBLVLKDIA(KJ,:),KJ,:)=NLVLKDIA(NBLVLKDIA(KJ,1),KJ,1)
      END IF
    ENDDO


  END IF
!
END IF
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
if(nverbia >0)then
print *,' RESOLVK '
print *,' LVLKDIALL ',LVLKDIALL(KJ,1)
print *,' NBLVLKDIA ',NBLVLKDIA(KJ,1)
print *,' NLVLKDIA ',(NLVLKDIA(J,KJ,1),J=1,NBLVLKDIA(KJ,1))
endif
RETURN
END SUBROUTINE RESOLVK  
!     ######spl
      MODULE MODI_RESOLVN
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVN(HCARIN,KINDN,KJ)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDN, KJ
END SUBROUTINE RESOLVN
!
END INTERFACE
!
END MODULE MODI_RESOLVN
!     ###################################
      SUBROUTINE RESOLVN(HCARIN,KINDN,KJ)
!     ###################################
!
!!****  *RESOLVN* - 
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
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDN, KJ
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=80) :: YCART
CHARACTER(LEN=20) :: YCAR
INTEGER          :: ILENC, ILENCART
INTEGER          :: INDPF, INDTO, INDBY, INDV, INDVM
INTEGER          :: ICAS, J

!
!------------------------------------------------------------------------------
INDPF = 0
INDTO = 0
INDBY = 0
INDV  = 0
ICAS = 0

NBNDIA(KJ)=0
NNDIA(:,KJ)=0
LNDIALL(KJ)=.FALSE.

IF(KINDN == 0)THEN
  LNDIALL(KJ) = .TRUE.
  RETURN
END IF

ILENC = LEN(HCARIN)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INDTO = INDEX(HCARIN(KINDN+3:ILENC),'_TO_')
  INDBY = INDEX(HCARIN(KINDN+3:ILENC),'_BY_')
  INDPF = INDEX(HCARIN(KINDN+3:ILENC),'_')
  IF(INDTO /= 0)THEN
  IF(INDPF < INDTO)THEN
!
! ICAS = 1  Num. unique ou separes par des virgules
!
    INDTO=0;INDBY=0
    ICAS = 1
  ELSE IF(INDPF == INDTO)THEN
!
! ICAS = 3  Proc1 _TO_ Procn _BY_ Procx
!
    IF(INDBY /= 0)THEN
      DO J=INDTO+4+KINDN+3,INDBY+KINDN+3
        IF(HCARIN(J:J) == '_')THEN
          IF(HCARIN(J:J+3) == '_BY_')THEN
            EXIT
          ELSE
            INDBY=0
            EXIT
          END IF
        END IF
      ENDDO
    END IF
    IF(INDBY /= 0)THEN
      INDPF=INDEX(HCARIN(KINDN+3+INDBY+4:ILENC),'_')
      IF(INDPF /= 0)INDPF=INDPF+INDBY+4
      ICAS = 3
      LPINCRDIA(KJ) = .TRUE.
    ELSE
!
! ICAS = 2  Num1 _TO_ Numn
!
      INDPF=INDEX(HCARIN(KINDN+3+INDTO+4:ILENC),'_')
      IF(INDPF /= 0)INDPF=INDPF+INDTO+4
      ICAS = 2
      LPINCRDIA(KJ) = .TRUE.
    END IF
  END IF
  ELSE
    ICAS = 1
  END IF
IF(INDPF == 0)THEN
  INDPF = ILENC
ELSE
  INDPF = INDPF+KINDN+3-1-1
END IF


YCART(1:LEN(YCART))=' '
YCAR(1:LEN(YCAR))=' '
!
! Extraction de la partie Numeros (masques ou traj.) dans YCART(1:ILENCART)
!
!print *,' KINDN INDPF ',KINDN,INDPF
YCART = ADJUSTL(HCARIN(KINDN+3:INDPF))
ILENCART = LEN_TRIM(YCART)
!print *,' YCART ',ILENCART,' ',YCART

! Recherche a nouveau des chaines de car. _TO_ , _BY_ et d'une virgule
! par rapport au debut de YCART

INDTO = INDEX(YCART,'_TO_')
INDBY = INDEX(YCART,'_BY_')
INDV = INDEX(YCART(1:ILENCART),',')
IF(ICAS == 1 .AND. INDV == 0)ICAS=0
!
! Expression des Numeros par mots-cles (NALL ou N1....)
!
IF(YCART(1:4) == 'NALL')THEN
  LNDIALL(KJ) = .TRUE.
  if (nverbia>0) then
  print *,' RESOLVN NALL '
  print *,' LNDIALL ',LNDIALL(KJ)
  print *,' NBNDIA ',NBNDIA(KJ)
  print *,' NNDIA ',(NNDIA(J,KJ),J=1,NBNDIA(KJ))
  endif
  RETURN

ELSE IF(YCART(1:1) == 'N')THEN
!print *,' YCART(1:1) ',YCART(1:1),' ICAS ',ICAS

  NBNDIA(KJ)=NBNDIA(KJ)+1
  SELECT CASE(ICAS)
    CASE(1)
!print *,' INDV YCART(2:2) ',INDV,YCART(2:2)
      IF(INDV-1-1 == 1)READ(YCART(2:2),'(I1)')NNDIA(NBNDIA(KJ),KJ)
      IF(INDV-1-1 == 2)READ(YCART(2:3),'(I2)')NNDIA(NBNDIA(KJ),KJ)
      DO J = 1,100
        INDVM=INDV
        INDV=0
        INDV=INDEX(YCART(INDVM+1:ILENCART),',')
        IF(INDV == 0)THEN
          NBNDIA(KJ)=NBNDIA(KJ)+1
          IF(ILENCART-(INDVM+1) == 1)READ(YCART(INDVM+1+1:INDVM+1+1),'(I1)')NNDIA(NBNDIA(KJ),KJ)
          IF(ILENCART-(INDVM+1) == 2)READ(YCART(INDVM+1+1:ILENCART),'(I2)')NNDIA(NBNDIA(KJ),KJ)
          EXIT
        ELSE
          INDV=INDV+INDVM
          NBNDIA(KJ)=NBNDIA(KJ)+1
          IF(INDV-(INDVM+1)-1 == 1)READ(YCART(INDVM+1+1:INDVM+1+1),'(I1)')NNDIA(NBNDIA(KJ),KJ)
          IF(INDV-(INDVM+1)-1 == 2)READ(YCART(INDVM+1+1:INDV-1),'(I2)')NNDIA(NBNDIA(KJ),KJ)
        END IF
      ENDDO   
      
    CASE(2)
      IF(INDTO-1-1 == 1)READ(YCART(2:2),'(I1)')NNDIA(NBNDIA(KJ),KJ)
      IF(INDTO-1-1 == 2)READ(YCART(2:3),'(I2)')NNDIA(NBNDIA(KJ),KJ)
      NBNDIA(KJ)=NBNDIA(KJ)+1
      IF(ILENCART-(INDTO+3+1) == 1)READ(YCART(INDTO+3+1+1:INDTO+3+1+1),'(I1)')NNDIA(NBNDIA(KJ),KJ)
      IF(ILENCART-(INDTO+3+1) == 2)READ(YCART(INDTO+3+1+1:ILENCART),'(I2)')NNDIA(NBNDIA(KJ),KJ)
! 1 seul temps
    CASE DEFAULT
      IF(ILENCART-1 == 1)READ(YCART(2:2),'(I1)')NNDIA(NBNDIA(KJ),KJ)
      IF(ILENCART-1 == 2)READ(YCART(2:3),'(I2)')NNDIA(NBNDIA(KJ),KJ)

  END SELECT
  if (nverbia>0) then
  print *,' RESOLVN ICAS '
  print *,' LNDIALL ',LNDIALL(KJ)
  print *,' NBNDIA ',NBNDIA(KJ)
  print *,' NNDIA ',(NNDIA(J,KJ),J=1,NBNDIA(KJ))
  endif
  RETURN
ELSE

!
! Expression des numeros en numerique
!
  IF(INDV == 0)THEN

! Cas  _TO_  _BY_

    IF(INDTO /= 0)THEN
      YCAR = ADJUSTL(YCART(1:INDTO-1))
      NBNDIA(KJ) = NBNDIA(KJ)+1
      CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
      IF(INDBY /= 0)THEN
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:INDBY-1))
        NBNDIA(KJ) = NBNDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDBY+4:ILENCART))
        NBNDIA(KJ) = NBNDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
      ELSE
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:ILENCART))
        NBNDIA(KJ) = NBNDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
      END IF
    ELSE

! Cas un seul processus en fin de chaine de car. HCARIN ou au milieu

      IF(ILENCART > 9)THEN
	print *,' PB ecriture temps '
	STOP
      ELSE
	YCAR = ADJUSTL(YCART(1:ILENCART))
	NBNDIA(KJ) = NBNDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
      END IF

    END IF

  ELSE

! Presence de virgules

    YCAR = ADJUSTL(YCART(1:INDV-1))
    NBNDIA(KJ) = NBNDIA(KJ)+1
    CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
    DO J = 1,100
      INDVM=INDV
      INDV = 0
      YCAR(1:LEN(YCAR))=' '
      INDV = INDEX(YCART(INDVM+1:ILENCART),',')
!     print *,' INDV ',INDV
      IF(INDV == 0)THEN
	YCAR = ADJUSTL(YCART(INDVM+1:ILENCART))
	NBNDIA(KJ) = NBNDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
	EXIT
      ELSE
        INDV=INDV+INDVM
	YCAR = ADJUSTL(YCART(INDVM+1:INDV-1))
	NBNDIA(KJ) = NBNDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NNDIA(NBNDIA(KJ),KJ))
      END IF
    ENDDO


  END IF
!
END IF
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
if (nverbia>0) then
print *,' end of RESOLVN '
print *,' LNDIALL ',LNDIALL(KJ)
print *,' NBNDIA ',NBNDIA(KJ)
print *,' NNDIA ',(NNDIA(J,KJ),J=1,NBNDIA(KJ))
endif
RETURN
END SUBROUTINE RESOLVN  
!     ######spl
      MODULE MODI_RESOLVON
!     ####################
!
INTERFACE
!
SUBROUTINE RESOLVON(HCARIN,KINDON)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDON
END SUBROUTINE RESOLVON
!
END INTERFACE
!
END MODULE MODI_RESOLVON
!     ##################################
      SUBROUTINE RESOLVON(HCARIN,KINDON)
!     ##################################
!
!!****  *RESOLVON* - 
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
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDON
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=LEN_TRIM(HCARIN)) :: YCARIN
INTEGER          :: ILENC, INDON, INDONM, ILONMS, INDMINUS, INDPLUS
INTEGER          :: J
LOGICAL          :: OMINUS, OPLUS

!
!------------------------------------------------------------------------------
OMINUS=LMINUS
OPLUS=LPLUS
LSUPERDIA=.TRUE.
NSUPERDIA=NSUPERDIA+1
ILENC=LEN_TRIM(HCARIN)
CARSUP(NSUPERDIA)(1:KINDON-1)=HCARIN(1:KINDON-1)
INDONM=KINDON
IF(LMINUS)THEN
  ILONMS=7
  NBPM=NBPM+1
  NUMPM(NBPM)=2
ELSE IF(LPLUS)THEN
  ILONMS=6
  NBPM=NBPM+1
  NUMPM(NBPM)=1
ELSE
  ILONMS=4
  NBPM=NBPM+1
  NUMPM(NBPM)=3
ENDIF
DO J=1,100
YCARIN(1:LEN(YCARIN))=' '
YCARIN(1:ILENC-INDONM-ILONMS+1)=ADJUSTL(HCARIN(INDONM+ILONMS:ILENC))
INDONM=INDONM+(ILONMS-1)
INDON=INDEX(YCARIN,'_ON_')
INDMINUS=INDEX(YCARIN,'_MINUS_')
INDPLUS=INDEX(YCARIN,'_PLUS_')
IF(INDON == 0)THEN
  IF(INDMINUS == 0)THEN
    IF(INDPLUS == 0)THEN
    ELSE
      INDON=INDPLUS
      NBPM=NBPM+1
      NUMPM(NBPM)=1
      ILONMS=6
    ENDIF
  ELSE
    IF(INDPLUS == 0)THEN
      INDON=INDMINUS
      NBPM=NBPM+1
      NUMPM(NBPM)=2
      ILONMS=7
    ELSE
      IF(INDMINUS < INDPLUS)THEN
	INDON=INDMINUS
        NBPM=NBPM+1
        NUMPM(NBPM)=2
	ILONMS=7
      ELSE
        INDON=INDPLUS
        NBPM=NBPM+1
        NUMPM(NBPM)=1
        ILONMS=6
      ENDIF
    ENDIF
  ENDIF

ELSE

! INDON =/= 0

  IF(INDMINUS == 0 .AND. INDPLUS == 0)THEN
    NBPM=NBPM+1
    NUMPM(NBPM)=3
    ILONMS=4
  ELSE
    IF(INDMINUS == 0)THEN
      IF(INDON < INDPLUS)THEN
        NBPM=NBPM+1
        NUMPM(NBPM)=3
        ILONMS=4
      ELSE
        INDON=INDPLUS
        NBPM=NBPM+1
        NUMPM(NBPM)=1
        ILONMS=6
      ENDIF
    ELSE
      IF(INDPLUS == 0)THEN
        IF(INDON < INDMINUS)THEN
          NBPM=NBPM+1
          NUMPM(NBPM)=3
          ILONMS=4
        ELSE
          INDON=INDMINUS
          NBPM=NBPM+1
          NUMPM(NBPM)=2
          ILONMS=7
        ENDIF
      ELSE
! ON + et -
        IF(INDON < INDMINUS .AND. INDON < INDPLUS)THEN
          NBPM=NBPM+1
          NUMPM(NBPM)=3
          ILONMS=4
        ELSE IF(INDMINUS < INDON .AND. INDMINUS < INDPLUS)THEN
          INDON=INDMINUS
          NBPM=NBPM+1
          NUMPM(NBPM)=2
          ILONMS=7
        ELSE IF(INDPLUS < INDON .AND. INDPLUS < INDMINUS)THEN
          INDON=INDPLUS
          NBPM=NBPM+1
          NUMPM(NBPM)=1
          ILONMS=6
        ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDIF
IF(INDON == 0)THEN
  NSUPERDIA=NSUPERDIA+1
  CARSUP(NSUPERDIA)(1:LEN_TRIM(YCARIN))=ADJUSTL(YCARIN(1:LEN_TRIM(YCARIN)))
EXIT
ELSE
  NSUPERDIA=NSUPERDIA+1
  CARSUP(NSUPERDIA)(1:INDON-1)=ADJUSTL(YCARIN(1:INDON-1))
  INDONM=INDONM+INDON
ENDIF
ENDDO
NBPMT=0
DO J=1,NBPM
  IF(NUMPM(J) == 1 .OR. NUMPM(J) == 2)THEN
    NBPMT=NBPMT+1
  ENDIF
ENDDO
LMINUS=OMINUS
LPLUS=OPLUS
!print *,' resolvon NBPM NUMPM ',NBPM,NUMPM(1:NBPM)
if(nverbia >0)then
print *,'resolvon NBPM NUMPM ',NBPM,NUMPM(1:NBPM)
endif
RETURN
END SUBROUTINE RESOLVON  
!     ######spl
      MODULE MODI_RESOLVP
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVP(HCARIN,KINDP,KJ)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDP, KJ
END SUBROUTINE RESOLVP
!
END INTERFACE
!
END MODULE MODI_RESOLVP
!     ###################################
      SUBROUTINE RESOLVP(HCARIN,KINDP,KJ)
!     ###################################
!
!!****  *RESOLVP* - 
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
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDP, KJ
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=80) :: YCART
CHARACTER(LEN=20) :: YCAR
INTEGER          :: ILENC, ILENCART
INTEGER          :: INDPF, INDTO, INDBY, INDV, INDVM
INTEGER          :: ICAS, J

!
!------------------------------------------------------------------------------
INDPF = 0
INDTO = 0
INDBY = 0
INDV  = 0
ICAS = 0

NBPROCDIA(KJ)=0
NPROCDIA(:,KJ)=0
LPROCDIALL(KJ)=.FALSE.

IF(KINDP == 0)THEN
  LPROCDIALL(KJ) = .TRUE.
  RETURN
END IF

ILENC = LEN(HCARIN)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INDTO = INDEX(HCARIN(KINDP+3:ILENC),'_TO_')
  INDBY = INDEX(HCARIN(KINDP+3:ILENC),'_BY_')
  INDPF = INDEX(HCARIN(KINDP+3:ILENC),'_')
  IF(INDTO /= 0)THEN
  IF(INDPF < INDTO)THEN
!
! ICAS = 1  Proc. unique ou separes par des virgules
!
    INDTO=0;INDBY=0
    ICAS = 1
  ELSE IF(INDPF == INDTO)THEN
!
! ICAS = 3  Proc1 _TO_ Procn _BY_ Procx
!
    IF(INDBY /= 0)THEN
      DO J=INDTO+4+KINDP+3,INDBY+KINDP+3
        IF(HCARIN(J:J) == '_')THEN
          IF(HCARIN(J:J+3) == '_BY_')THEN
            EXIT
          ELSE
            INDBY=0
            EXIT
          END IF
        END IF
      ENDDO
    END IF
    IF(INDBY /= 0)THEN
      INDPF=INDEX(HCARIN(KINDP+3+INDBY+4:ILENC),'_')
      IF(INDPF /= 0)INDPF=INDPF+INDBY+4
      ICAS = 3
      LPINCRDIA(KJ) = .TRUE.
    ELSE
!
! ICAS = 2  Proc1 _TO_ Procn
!
      INDPF=INDEX(HCARIN(KINDP+3+INDTO+4:ILENC),'_')
      IF(INDPF /= 0)INDPF=INDPF+INDTO+4
      ICAS = 2
      LPINCRDIA(KJ) = .TRUE.
    END IF
  END IF
  ELSE
    ICAS = 1
  END IF
IF(INDPF == 0)THEN
  INDPF = ILENC
ELSE
  INDPF = INDPF+KINDP+3-1-1
END IF


YCART(1:LEN(YCART))=' '
YCAR(1:LEN(YCAR))=' '
!
! Extraction de la partie Processus dans YCART(1:ILENCART)
!
!print *,' KINDP INDPF ',KINDP,INDPF
YCART = ADJUSTL(HCARIN(KINDP+3:INDPF))
ILENCART = LEN_TRIM(YCART)
!print *,' YCART ',ILENCART,' ',YCART

! Recherche a nouveau des chaines de car. _TO_ , _BY_ et d'une virgule
! par rapport au debut de YCART

INDTO = INDEX(YCART,'_TO_')
INDBY = INDEX(YCART,'_BY_')
INDV = INDEX(YCART(1:ILENCART),',')
IF(ICAS == 1 .AND. INDV == 0)ICAS=0
!
! Expression des Processus par mots-cles (PROCALL ou PROC1....)
!
IF(YCART(1:7) == 'PROCALL')THEN
  LPROCDIALL(KJ) = .TRUE.
  print *,' RESOLVP PROCALL '
  print *,' LPROCDIALL ',LPROCDIALL(KJ)
  print *,' NBPROCDIA ',NBPROCDIA(KJ)
  print *,' NPROCDIA ',(NPROCDIA(J,KJ),J=1,NBPROCDIA(KJ))
  RETURN

ELSE IF(YCART(1:4) == 'PROC')THEN
!print *,' YCART(1:4) ',YCART(1:4),' ICAS ',ICAS

  NBPROCDIA(KJ)=NBPROCDIA(KJ)+1
  SELECT CASE(ICAS)
    CASE(1)
!print *,' INDV YCART(5:5) ',INDV,YCART(5:5)
      IF(INDV-4-1 == 1)READ(YCART(5:5),'(I1)')NPROCDIA(NBPROCDIA(KJ),KJ)
      IF(INDV-4-1 == 2)READ(YCART(5:6),'(I2)')NPROCDIA(NBPROCDIA(KJ),KJ)
      DO J = 1,100
        INDVM=INDV
        INDV=0
        INDV=INDEX(YCART(INDVM+1:ILENCART),',')
        IF(INDV == 0)THEN
          NBPROCDIA(KJ)=NBPROCDIA(KJ)+1
          IF(ILENCART-(INDVM+4) == 1)READ(YCART(INDVM+4+1:INDVM+4+1),'(I1)')NPROCDIA(NBPROCDIA(KJ),KJ)
          IF(ILENCART-(INDVM+4) == 2)READ(YCART(INDVM+4+1:ILENCART),'(I2)')NPROCDIA(NBPROCDIA(KJ),KJ)
          EXIT
        ELSE
          INDV=INDV+INDVM
          NBPROCDIA(KJ)=NBPROCDIA(KJ)+1
          IF(INDV-(INDVM+4)-1 == 1)READ(YCART(INDVM+4+1:INDVM+4+1),'(I1)')NPROCDIA(NBPROCDIA(KJ),KJ)
          IF(INDV-(INDVM+4)-1 == 2)READ(YCART(INDVM+4+1:INDV-1),'(I2)')NPROCDIA(NBPROCDIA(KJ),KJ)
        END IF
      ENDDO   
      
    CASE(2)
      IF(INDTO-4-1 == 1)READ(YCART(5:5),'(I1)')NPROCDIA(NBPROCDIA(KJ),KJ)
      IF(INDTO-4-1 == 2)READ(YCART(5:6),'(I2)')NPROCDIA(NBPROCDIA(KJ),KJ)
      NBPROCDIA(KJ)=NBPROCDIA(KJ)+1
      IF(ILENCART-(INDTO+3+4) == 1)READ(YCART(INDTO+3+4+1:INDTO+3+4+1),'(I1)')NPROCDIA(NBPROCDIA(KJ),KJ)
      IF(ILENCART-(INDTO+3+4) == 2)READ(YCART(INDTO+3+4+1:ILENCART),'(I2)')NPROCDIA(NBPROCDIA(KJ),KJ)
! 1 seul temps
    CASE DEFAULT
      IF(ILENCART-4 == 1)READ(YCART(5:5),'(I1)')NPROCDIA(NBPROCDIA(KJ),KJ)
      IF(ILENCART-4 == 2)READ(YCART(5:6),'(I2)')NPROCDIA(NBPROCDIA(KJ),KJ)

  END SELECT
  print *,' RESOLVP ICAS '
  print *,' LPROCDIALL ',LPROCDIALL(KJ)
  print *,' NBPROCDIA ',NBPROCDIA(KJ)
  print *,' NPROCDIA ',(NPROCDIA(J,KJ),J=1,NBPROCDIA(KJ))
  RETURN
ELSE

!
! Expression des processus en numerique
!
  IF(INDV == 0)THEN

! Cas  _TO_  _BY_

    IF(INDTO /= 0)THEN
      YCAR = ADJUSTL(YCART(1:INDTO-1))
      NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
      CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
      IF(INDBY /= 0)THEN
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:INDBY-1))
        NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDBY+4:ILENCART))
        NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
      ELSE
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:ILENCART))
        NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
      END IF
    ELSE

! Cas un seul processus en fin de chaine de car. HCARIN ou au milieu

      IF(ILENCART > 9)THEN
	print *,' PB ecriture temps '
	STOP
      ELSE
	YCAR = ADJUSTL(YCART(1:ILENCART))
	NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
      END IF

    END IF

  ELSE

! Presence de virgules

    YCAR = ADJUSTL(YCART(1:INDV-1))
    NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
    CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
    DO J = 1,100
      INDVM=INDV
      INDV = 0
      YCAR(1:LEN(YCAR))=' '
      INDV = INDEX(YCART(INDVM+1:ILENCART),',')
!     print *,' INDV ',INDV
      IF(INDV == 0)THEN
	YCAR = ADJUSTL(YCART(INDVM+1:ILENCART))
	NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
	EXIT
      ELSE
        INDV=INDV+INDVM
	YCAR = ADJUSTL(YCART(INDVM+1:INDV-1))
	NBPROCDIA(KJ) = NBPROCDIA(KJ)+1
        CALL CARINT(YCAR(1:LEN_TRIM(YCAR)),NPROCDIA(NBPROCDIA(KJ),KJ))
      END IF
    ENDDO


  END IF
!
END IF
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
print *,' RESOLVP '
print *,' LPROCDIALL ',LPROCDIALL(KJ)
print *,' NBPROCDIA ',NBPROCDIA(KJ)
print *,' NPROCDIA ',(NPROCDIA(J,KJ),J=1,NBPROCDIA(KJ))
RETURN
END SUBROUTINE RESOLVP  
!     ######spl
      MODULE MODI_RESOLVX
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVX(HCARIN,KIND,POUT)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KIND
REAL             :: POUT
END SUBROUTINE RESOLVX
!
END INTERFACE
!
END MODULE MODI_RESOLVX
!     ########################################
      SUBROUTINE RESOLVX(HCARIN,KIND,POUT)
!     ########################################
!
!!****  *RESOLVX* - 
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
USE MODD_RESOLVCAR
USE MODN_PARA

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KIND
REAL             :: POUT
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=15) :: YC15
INTEGER           :: ILENC, ILENC15
INTEGER           :: J,JM

!
!------------------------------------------------------------------------------
ILENC=LEN_TRIM(HCARIN)

DO J=KIND,ILENC
  IF(HCARIN(J:J) == '=')EXIT
ENDDO

JM=J+1
DO J=1,10
  IF(HCARIN(JM:JM) == ' ')THEN
    JM=JM+1
  ELSE
    EXIT
  ENDIF
ENDDO
YC15='               '

DO J=JM,ILENC
  IF(HCARIN(J:J) == '0'.OR.HCARIN(J:J) == '1'.OR.HCARIN(J:J) == '2'  &
    .OR.HCARIN(J:J) == '3'.OR.HCARIN(J:J) == '4'.OR.HCARIN(J:J) == '5' &
    .OR.HCARIN(J:J) == '6'.OR.HCARIN(J:J) == '7'.OR.HCARIN(J:J) == '8' &
    .OR.HCARIN(J:J) == '9'.OR.HCARIN(J:J) == '.' .OR. &
    HCARIN(J:J) == '+' .OR.HCARIN(J:J) == '-' .OR.HCARIN(J:J) == 'E' &
    .OR.HCARIN(J:J) == 'e')THEN
    YC15(1:1)=HCARIN(J:J)
     IF(J+1 > ILENC)EXIT
    IF(HCARIN(J+1:J+1) /= '0' .AND. HCARIN(J+1:J+1) /= '1' .AND.  &
       HCARIN(J+1:J+1) /= '2' .AND. HCARIN(J+1:J+1) /= '3' .AND.  &
       HCARIN(J+1:J+1) /= '4' .AND. HCARIN(J+1:J+1) /= '5' .AND.  &
       HCARIN(J+1:J+1) /= '6' .AND. HCARIN(J+1:J+1) /= '7' .AND.  &
       HCARIN(J+1:J+1) /= '8' .AND. HCARIN(J+1:J+1) /= '9' .AND.  &
       HCARIN(J+1:J+1) /= '+' .AND. HCARIN(J+1:J+1) /= '-' .AND.  &
       HCARIN(J+1:J+1) /= 'E' .AND. HCARIN(J+1:J+1) /= 'e' .AND.  &
       HCARIN(J+1:J+1) /= '.')THEN
       EXIT
    ELSE
      YC15(2:2)=HCARIN(J+1:J+1)
      IF(J+2 > ILENC)EXIT
      IF(HCARIN(J+2:J+2) /= '0' .AND. HCARIN(J+2:J+2) /= '1' .AND.  &
	 HCARIN(J+2:J+2) /= '2' .AND. HCARIN(J+2:J+2) /= '3' .AND.  &
	 HCARIN(J+2:J+2) /= '4' .AND. HCARIN(J+2:J+2) /= '5' .AND.  &
	 HCARIN(J+2:J+2) /= '6' .AND. HCARIN(J+2:J+2) /= '7' .AND.  &
	 HCARIN(J+2:J+2) /= '8' .AND. HCARIN(J+2:J+2) /= '9' .AND.  &
	 HCARIN(J+2:J+2) /= '+' .AND. HCARIN(J+2:J+2) /= '-' .AND.  &
	 HCARIN(J+2:J+2) /= 'E' .AND. HCARIN(J+2:J+2) /= 'e' .AND.  &
	 HCARIN(J+2:J+2) /= '.')THEN
	 EXIT
      ELSE
	YC15(3:3)=HCARIN(J+2:J+2)
        IF(J+3 > ILENC)EXIT
	IF(HCARIN(J+3:J+3) /= '0' .AND. HCARIN(J+3:J+3) /= '1' .AND.  &
	   HCARIN(J+3:J+3) /= '2' .AND. HCARIN(J+3:J+3) /= '3' .AND.  &
	   HCARIN(J+3:J+3) /= '4' .AND. HCARIN(J+3:J+3) /= '5' .AND.  &
	   HCARIN(J+3:J+3) /= '6' .AND. HCARIN(J+3:J+3) /= '7' .AND.  &
	   HCARIN(J+3:J+3) /= '8' .AND. HCARIN(J+3:J+3) /= '9' .AND.  &
	   HCARIN(J+3:J+3) /= '+' .AND. HCARIN(J+3:J+3) /= '-' .AND.  &
	   HCARIN(J+3:J+3) /= 'E' .AND. HCARIN(J+3:J+3) /= 'e' .AND.  &
	   HCARIN(J+3:J+3) /= '.')THEN
	   EXIT
        ELSE
  	  YC15(4:4)=HCARIN(J+3:J+3)
          IF(J+4 > ILENC)EXIT
  	  IF(HCARIN(J+4:J+4) /= '0' .AND. HCARIN(J+4:J+4) /= '1' .AND. &
  	     HCARIN(J+4:J+4) /= '2' .AND. HCARIN(J+4:J+4) /= '3' .AND. &
  	     HCARIN(J+4:J+4) /= '4' .AND. HCARIN(J+4:J+4) /= '5' .AND. &
  	     HCARIN(J+4:J+4) /= '6' .AND. HCARIN(J+4:J+4) /= '7' .AND. &
  	     HCARIN(J+4:J+4) /= '8' .AND. HCARIN(J+4:J+4) /= '9' .AND. &
  	     HCARIN(J+4:J+4) /= '+' .AND. HCARIN(J+4:J+4) /= '-' .AND. &
  	     HCARIN(J+4:J+4) /= 'E' .AND. HCARIN(J+4:J+4) /= 'e' .AND. &
	     HCARIN(J+4:J+4) /= '.')THEN
  	     EXIT
           ELSE
  	     YC15(5:5)=HCARIN(J+4:J+4)
             IF(J+5 > ILENC)EXIT
  	     IF(HCARIN(J+5:J+5) /= '0' .AND. HCARIN(J+5:J+5) /= '1' .AND. &
  	        HCARIN(J+5:J+5) /= '2' .AND. HCARIN(J+5:J+5) /= '3' .AND. &
  	        HCARIN(J+5:J+5) /= '4' .AND. HCARIN(J+5:J+5) /= '5' .AND. &
  	        HCARIN(J+5:J+5) /= '6' .AND. HCARIN(J+5:J+5) /= '7' .AND. &
  	        HCARIN(J+5:J+5) /= '8' .AND. HCARIN(J+5:J+5) /= '9' .AND. &
  	        HCARIN(J+5:J+5) /= '+' .AND. HCARIN(J+5:J+5) /= '-' .AND. &
  	        HCARIN(J+5:J+5) /= 'E' .AND. HCARIN(J+5:J+5) /= 'e' .AND. &
		HCARIN(J+5:J+5) /= '.')THEN
  	        EXIT
              ELSE
  	        YC15(6:6)=HCARIN(J+5:J+5)
                IF(J+6 > ILENC)EXIT
  	        IF(HCARIN(J+6:J+6) /= '0' .AND. HCARIN(J+6:J+6) /= '1' .AND. &
  	           HCARIN(J+6:J+6) /= '2' .AND. HCARIN(J+6:J+6) /= '3' .AND. &
  	           HCARIN(J+6:J+6) /= '4' .AND. HCARIN(J+6:J+6) /= '5' .AND. &
  	           HCARIN(J+6:J+6) /= '6' .AND. HCARIN(J+6:J+6) /= '7' .AND. &
  	           HCARIN(J+6:J+6) /= '8' .AND. HCARIN(J+6:J+6) /= '9' .AND. &
  	           HCARIN(J+6:J+6) /= '+' .AND. HCARIN(J+6:J+6) /= '-' .AND. &
  	           HCARIN(J+6:J+6) /= 'E' .AND. HCARIN(J+6:J+6) /= 'e' .AND. &
		   HCARIN(J+6:J+6) /= '.')THEN
  	           EXIT
                 ELSE
  	           YC15(7:7)=HCARIN(J+6:J+6)
                   IF(J+7 > ILENC)EXIT
  	           IF(HCARIN(J+7:J+7) /= '0' .AND. HCARIN(J+7:J+7) /= '1' .AND.&
  	              HCARIN(J+7:J+7) /= '2' .AND. HCARIN(J+7:J+7) /= '3' .AND.&
  	              HCARIN(J+7:J+7) /= '4' .AND. HCARIN(J+7:J+7) /= '5' .AND.&
  	              HCARIN(J+7:J+7) /= '6' .AND. HCARIN(J+7:J+7) /= '7' .AND.&
  	              HCARIN(J+7:J+7) /= '8' .AND. HCARIN(J+7:J+7) /= '9' .AND.&
  	              HCARIN(J+7:J+7) /= '+' .AND. HCARIN(J+7:J+7) /= '-' .AND.&
  	              HCARIN(J+7:J+7) /= 'E' .AND. HCARIN(J+7:J+7) /= 'e' .AND.&
		      HCARIN(J+7:J+7) /= '.')THEN
  	              EXIT
                    ELSE
		      YC15(8:8)=HCARIN(J+7:J+7)
                      IF(J+8 > ILENC)EXIT
		      IF(HCARIN(J+8:J+8) /= '0' .AND.  &
			 HCARIN(J+8:J+8) /= '1' .AND.  &
                         HCARIN(J+8:J+8) /= '2' .AND.  &
			 HCARIN(J+8:J+8) /= '3' .AND.  &
		         HCARIN(J+8:J+8) /= '4' .AND.  &
			 HCARIN(J+8:J+8) /= '5' .AND.  &
                         HCARIN(J+8:J+8) /= '6' .AND.  &
			 HCARIN(J+8:J+8) /= '7' .AND.  &
                         HCARIN(J+8:J+8) /= '8' .AND.  &
			 HCARIN(J+8:J+8) /= '9' .AND.  &
			 HCARIN(J+8:J+8) /= '+' .AND.  &
			 HCARIN(J+8:J+8) /= '-' .AND.  &
			 HCARIN(J+8:J+8) /= 'E' .AND.  &
			 HCARIN(J+8:J+8) /= 'e' .AND.  &
                         HCARIN(J+8:J+8) /= '.')THEN
			 EXIT
                       ELSE
		         YC15(9:9)=HCARIN(J+8:J+8)
                         IF(J+9 > ILENC)EXIT
		         IF(HCARIN(J+9:J+9) /= '0' .AND.  &
			    HCARIN(J+9:J+9) /= '1' .AND.  &
                            HCARIN(J+9:J+9) /= '2' .AND.  &
			    HCARIN(J+9:J+9) /= '3' .AND.  &
		            HCARIN(J+9:J+9) /= '4' .AND.  &
			    HCARIN(J+9:J+9) /= '5' .AND.  &
                            HCARIN(J+9:J+9) /= '6' .AND.  &
			    HCARIN(J+9:J+9) /= '7' .AND.  &
                            HCARIN(J+9:J+9) /= '8' .AND.  &
			    HCARIN(J+9:J+9) /= '9' .AND.  &
			    HCARIN(J+9:J+9) /= '+' .AND.  &
			    HCARIN(J+9:J+9) /= '-' .AND.  &
			    HCARIN(J+9:J+9) /= 'E' .AND.  &
			    HCARIN(J+9:J+9) /= 'e' .AND.  &
                            HCARIN(J+9:J+9) /= '.')THEN
			    EXIT
                          ELSE
		            YC15(10:10)=HCARIN(J+9:J+9)
                            IF(J+10 > ILENC)EXIT
		            IF(HCARIN(J+10:J+10) /= '0' .AND.  &
			      HCARIN(J+10:J+10) /= '1' .AND.  &
                              HCARIN(J+10:J+10) /= '2' .AND.  &
			      HCARIN(J+10:J+10) /= '3' .AND.  &
		              HCARIN(J+10:J+10) /= '4' .AND.  &
			      HCARIN(J+10:J+10) /= '5' .AND.  &
                              HCARIN(J+10:J+10) /= '6' .AND.  &
			      HCARIN(J+10:J+10) /= '7' .AND.  &
                              HCARIN(J+10:J+10) /= '8' .AND.  &
			      HCARIN(J+10:J+10) /= '9' .AND.  &
			      HCARIN(J+10:J+10) /= '+' .AND.  &
			      HCARIN(J+10:J+10) /= '-' .AND.  &
			      HCARIN(J+10:J+10) /= 'E' .AND.  &
			      HCARIN(J+10:J+10) /= 'e' .AND.  &
                              HCARIN(J+10:J+10) /= '.')THEN
			      EXIT
                            ELSE
		              YC15(11:11)=HCARIN(J+10:J+10)
                               IF(J+11 > ILENC)EXIT
		              IF(HCARIN(J+11:J+11) /= '0' .AND.  &
                                 HCARIN(J+11:J+11) /= '1' .AND.  &
                                 HCARIN(J+11:J+11) /= '2' .AND.  &
			         HCARIN(J+11:J+11) /= '3' .AND.  &
		                 HCARIN(J+11:J+11) /= '4' .AND.  &
			         HCARIN(J+11:J+11) /= '5' .AND.  &
                                 HCARIN(J+11:J+11) /= '6' .AND.  &
			         HCARIN(J+11:J+11) /= '7' .AND.  &
                                 HCARIN(J+11:J+11) /= '8' .AND.  &
			         HCARIN(J+11:J+11) /= '9' .AND.  &
			         HCARIN(J+11:J+11) /= '+' .AND.  &
			         HCARIN(J+11:J+11) /= '-' .AND.  &
			         HCARIN(J+11:J+11) /= 'E' .AND.  &
			         HCARIN(J+11:J+11) /= 'e' .AND.  &
                                 HCARIN(J+11:J+11) /= '.')THEN
			         EXIT
                               ELSE
   		                 YC15(12:12)=HCARIN(J+11:J+11)
                                 IF(J+12 > ILENC)EXIT
   		                 IF(HCARIN(J+12:J+12) /= '0' .AND.  &
                                   HCARIN(J+12:J+12) /= '1' .AND.  &
                                   HCARIN(J+12:J+12) /= '2' .AND.  &
				   HCARIN(J+12:J+12) /= '3' .AND.  &
   		                   HCARIN(J+12:J+12) /= '4' .AND.  &
   			           HCARIN(J+12:J+12) /= '5' .AND.  &
                                   HCARIN(J+12:J+12) /= '6' .AND.  &
   			           HCARIN(J+12:J+12) /= '7' .AND.  &
                                   HCARIN(J+12:J+12) /= '8' .AND.  &
   			           HCARIN(J+12:J+12) /= '9' .AND.  &
   			           HCARIN(J+12:J+12) /= '+' .AND.  &
   			           HCARIN(J+12:J+12) /= '-' .AND.  &
   			           HCARIN(J+12:J+12) /= 'E' .AND.  &
   			           HCARIN(J+12:J+12) /= 'e' .AND.  &
                                   HCARIN(J+12:J+12) /= '.')THEN
   			           EXIT
                                 ELSE
   		                 YC15(13:13)=HCARIN(J+12:J+12)
                                 IF(J+13 > ILENC)EXIT
   		                 IF(HCARIN(J+13:J+13) /= '0' .AND.  &
                                   HCARIN(J+13:J+13) /= '1' .AND.  &
                                   HCARIN(J+13:J+13) /= '2' .AND.  &
				   HCARIN(J+13:J+13) /= '3' .AND.  &
   		                   HCARIN(J+13:J+13) /= '4' .AND.  &
   			           HCARIN(J+13:J+13) /= '5' .AND.  &
                                   HCARIN(J+13:J+13) /= '6' .AND.  &
   			           HCARIN(J+13:J+13) /= '7' .AND.  &
                                   HCARIN(J+13:J+13) /= '8' .AND.  &
   			           HCARIN(J+13:J+13) /= '9' .AND.  &
   			           HCARIN(J+13:J+13) /= '+' .AND.  &
   			           HCARIN(J+13:J+13) /= '-' .AND.  &
   			           HCARIN(J+13:J+13) /= 'E' .AND.  &
   			           HCARIN(J+13:J+13) /= 'e' .AND.  &
                                   HCARIN(J+13:J+13) /= '.')THEN
   			           EXIT
                                 ELSE
   		                 YC15(14:14)=HCARIN(J+13:J+13)
                                 IF(J+14 > ILENC)EXIT
   		                 IF(HCARIN(J+14:J+14) /= '0' .AND.  &
                                   HCARIN(J+14:J+14) /= '1' .AND.  &
                                   HCARIN(J+14:J+14) /= '2' .AND.  &
				   HCARIN(J+14:J+14) /= '3' .AND.  &
   		                   HCARIN(J+14:J+14) /= '4' .AND.  &
   			           HCARIN(J+14:J+14) /= '5' .AND.  &
                                   HCARIN(J+14:J+14) /= '6' .AND.  &
   			           HCARIN(J+14:J+14) /= '7' .AND.  &
                                   HCARIN(J+14:J+14) /= '8' .AND.  &
   			           HCARIN(J+14:J+14) /= '9' .AND.  &
   			           HCARIN(J+14:J+14) /= '+' .AND.  &
   			           HCARIN(J+14:J+14) /= '-' .AND.  &
   			           HCARIN(J+14:J+14) /= 'E' .AND.  &
   			           HCARIN(J+14:J+14) /= 'e' .AND.  &
                                   HCARIN(J+14:J+14) /= '.')THEN
   			           EXIT
                                 ELSE
   		                 YC15(15:15)=HCARIN(J+11:J+11)
                                 IF(J+15 > ILENC)EXIT
   		                 IF(HCARIN(J+15:J+15) /= '0' .AND.  &
                                   HCARIN(J+15:J+15) /= '1' .AND.  &
                                   HCARIN(J+15:J+15) /= '2' .AND.  &
				   HCARIN(J+15:J+15) /= '3' .AND.  &
   		                   HCARIN(J+15:J+15) /= '4' .AND.  &
   			           HCARIN(J+15:J+15) /= '5' .AND.  &
                                   HCARIN(J+15:J+15) /= '6' .AND.  &
   			           HCARIN(J+15:J+15) /= '7' .AND.  &
                                   HCARIN(J+15:J+15) /= '8' .AND.  &
   			           HCARIN(J+15:J+15) /= '9' .AND.  &
   			           HCARIN(J+15:J+15) /= '+' .AND.  &
   			           HCARIN(J+15:J+15) /= '-' .AND.  &
   			           HCARIN(J+15:J+15) /= 'E' .AND.  &
   			           HCARIN(J+15:J+15) /= 'e' .AND.  &
                                   HCARIN(J+15:J+15) /= '.')THEN
   			           EXIT
                                 ELSE
	                           print *,' PB AVEC LA VALEUR FOURNIE ', &
	                           HCARIN(J:J+15),' ARRET PG. VERIFIEZ SA VALEUR '
	                         ENDIF
	                         ENDIF
	                         ENDIF
	                         ENDIF
	                       ENDIF
                            ENDIF
		       ENDIF
	             ENDIF
                   ENDIF
	         ENDIF
	      ENDIF
	   ENDIF
	ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDDO


YC15=ADJUSTL(YC15)
ILENC15 = LEN_TRIM(YC15)
!print *, ' ILENC15 ',ILENC15,YC15
READ(YC15,*)POUT

RETURN
END SUBROUTINE RESOLVX
!     ######spl
      MODULE MODI_RESOLVXISOLEV
!     #########################
!
INTERFACE
!
SUBROUTINE RESOLVXISOLEV(HCARIN,KIND,POUT)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KIND
REAL,DIMENSION(:) :: POUT
END SUBROUTINE RESOLVXISOLEV
!
END INTERFACE
!
END MODULE MODI_RESOLVXISOLEV
!     ########################################
      SUBROUTINE RESOLVXISOLEV(HCARIN,KIND,POUT)
!     ########################################
!
!!****  *RESOLVXISOLEV* - 
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
USE MODD_RESOLVCAR
USE MODN_PARA

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KIND
REAL,DIMENSION(:)             :: POUT
!
!*       0.1   Local variables
!              ---------------

INTEGER           :: ILENC
INTEGER           :: J,JM, JMF
INTEGER           :: INBV, IND9999

!
!------------------------------------------------------------------------------
ILENC=LEN_TRIM(HCARIN)
POUT=9999.

DO J=KIND,ILENC
  IF(HCARIN(J:J) == '=')EXIT
ENDDO

JM=J+1
DO J=1,10
  IF(HCARIN(JM:JM) == ' ')THEN
    JM=JM+1
  ELSE
    EXIT
  ENDIF
ENDDO

IND9999=INDEX(HCARIN(JM:ILENC),'9999.')
JMF=IND9999+JM-1+3
INBV=0
IF(NVERBIA >= 5)THEN
  print *,' RESOLVXISOLEV carin: ',ind9999,jm,jmf,HCARIN(JM:JMF)
ENDIF
DO J=JM,JMF
  IF(HCARIN(J:J) == ',')THEN
    INBV=INBV+1
  ENDIF
ENDDO

READ(HCARIN(JM:JMF),*)(POUT(J),J=1,INBV+1)
IF(NVERBIA >= 5)THEN
  print *,' RESOLVXISOLEV ',INBV+1,(POUT(J),J=1,INBV+1)
ENDIF
RETURN
END SUBROUTINE RESOLVXISOLEV
!     ######spl
      MODULE MODI_RESOLVT
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVT(HCARIN,KINDT,KJ)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDT, KJ
END SUBROUTINE RESOLVT
!
END INTERFACE
!
END MODULE MODI_RESOLVT
!     ##################################
      SUBROUTINE RESOLVT(HCARIN,KINDT,KJ)
!     ###################################
!
!!****  *RESOLVT* - 
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
!!      Original       24/11/95
!!      Updated   PM   02/12/94
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
INTEGER          :: KINDT, KJ
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=LEN(HCARIN)) :: YCART
CHARACTER(LEN=20) :: YCAR
INTEGER          :: ILENC, ILENCART
INTEGER          :: INDTF, INDTO, INDBY, INDV, INDVM
INTEGER          :: ICAS, J

!
!------------------------------------------------------------------------------
INDTF = 0
INDTO = 0
INDBY = 0
INDV  = 0
ICAS = 0

NBTIMEDIA(KJ,:)=0
NTIMEDIA(:,KJ,:)=0
XTIMEDIA(:,KJ,:)=0.
LTIMEDIALL(KJ,:)=.FALSE.
LTINCRDIA(KJ,:)=.FALSE.

IF(KINDT == 0)THEN
  LTIMEDIALL(KJ,:) = .TRUE.
  RETURN
END IF

ILENC = LEN(HCARIN)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INDTO = INDEX(HCARIN(KINDT+3:ILENC),'_TO_')
INDBY = INDEX(HCARIN(KINDT+3:ILENC),'_BY_')
INDTF = INDEX(HCARIN(KINDT+3:ILENC),'_')
IF(INDTO /= 0)THEN
  IF(INDTF < INDTO)THEN
!
! ICAS = 1  Temps unique ou separes par des virgules
!
    INDTO=0;INDBY=0
    ICAS = 1
  ELSE IF(INDTF == INDTO)THEN
!
! ICAS = 3  Temps1 _TO_ Tempsn _BY_ Tempsx
!
    IF(INDBY /= 0)THEN
      DO J=INDTO+4+KINDT+3,INDBY+KINDT+3
        IF(HCARIN(J:J) == '_')THEN
          IF(HCARIN(J:J+3) == '_BY_')THEN
            EXIT
          ELSE
            INDBY=0
            EXIT
          END IF
        END IF
      ENDDO
    END IF
    IF(INDBY /= 0)THEN
      INDTF=INDEX(HCARIN(KINDT+3+INDBY+4:ILENC),'_')
      IF(INDTF /= 0)INDTF=INDTF+INDBY+4
      ICAS = 3
      LTINCRDIA(KJ,:) = .TRUE.
    ELSE
!
! ICAS = 2  Temps1 _TO_ Tempsn
!
      INDTF=INDEX(HCARIN(KINDT+3+INDTO+4:ILENC),'_')
      IF(INDTF /= 0)INDTF=INDTF+INDTO+4
      ICAS = 2
      LTINCRDIA(KJ,:) = .TRUE.
    END IF
  END IF
ELSE
  ICAS = 1
END IF

IF(INDTF == 0)THEN
  INDTF = ILENC
ELSE
  INDTF = INDTF+KINDT+3-1-1
END IF


YCART(1:LEN(YCART))=' '
YCAR(1:LEN(YCAR))=' '
!
! Extraction de la partie Temps dans YCART(1:ILENCART)
!
YCART = ADJUSTL(HCARIN(KINDT+3:INDTF))
ILENCART = LEN_TRIM(YCART)
if (nverbia >0) then
  print *,' ICAS KINDT INDTF ',ICAS,KINDT,INDTF
  print *,' YCART ',ILENCART,' ',YCART
endif

! Recherche a nouveau des chaines de car. _TO_ , _BY_ et d'une virgule
! par rapport au debut de YCART

INDTO = INDEX(YCART,'_TO_')
INDBY = INDEX(YCART,'_BY_')
INDV = INDEX(YCART(1:ILENCART),',')
IF(ICAS == 1 .AND. INDV == 0)ICAS=0
!
! Expression du temps par mots-cles (TIMEALL ou TIME1....)
!
IF(YCART(1:7) == 'TIMEALL')THEN
  LTIMEDIALL(KJ,:) = .TRUE.
if (nverbia >0) then
  print *,' RESOLVT TIMEALL '
  print *,' LTIMEDIALL(KJ,1) ',LTIMEDIALL(KJ,1)
  print *,' NBTIMEDIA(KJ,1) ',NBTIMEDIA(KJ,1)
  print *,' NTIMEDIA ',(NTIMEDIA(J,KJ,1),J=1,NBTIMEDIA(KJ,1))
  print *,' XTIMEDIA ',(XTIMEDIA(J,KJ,1),J=1,NBTIMEDIA(KJ,1))
endif
  RETURN

ELSE IF(YCART(1:4) == 'TIME')THEN
!print *,' YCART(1:4) ',YCART(1:4),' ICAS ',ICAS

  NBTIMEDIA(KJ,:)=NBTIMEDIA(KJ,:)+1
  SELECT CASE(ICAS)
    CASE(1)
!print *,' INDV YCART(5:5) ',INDV,YCART(5:5)
      READ(YCART(5:INDV-1),*)NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      NTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      DO J = 1,100
        INDVM=INDV
        INDV=0
        INDV=INDEX(YCART(INDVM+1:ILENCART),',')
        IF(INDV == 0)THEN
          NBTIMEDIA(KJ,:)=NBTIMEDIA(KJ,:)+1
          READ(YCART(INDVM+4+1:ILENCART),*)NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
	  NTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
          EXIT
        ELSE
          INDV=INDV+INDVM
          NBTIMEDIA(KJ,:)=NBTIMEDIA(KJ,:)+1
          READ(YCART(INDVM+4+1:INDV-1),*)NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
	  NTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
        END IF
      ENDDO   
      
    CASE(2)
      READ(YCART(5:INDTO-1),*)NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      NTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      NBTIMEDIA(KJ,:)=NBTIMEDIA(KJ,:)+1
      READ(YCART(INDTO+3+4+1:ILENCART),*)NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      NTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
! 1 seul temps
    CASE DEFAULT
      READ(YCART(5:ILENCART),*)NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      NTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=NTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)

  END SELECT
if (nverbia >0) then
  print *,' RESOLVT ICAS '
  print *,' LTIMEDIALL(KJ,1) ',LTIMEDIALL(KJ,1)
  print *,' NBTIMEDIA(KJ,1) ',NBTIMEDIA(KJ,1)
  print *,' NTIMEDIA ',(NTIMEDIA(J,KJ,1),J=1,NBTIMEDIA(KJ,1))
  print *,' XTIMEDIA ',(XTIMEDIA(J,KJ,1),J=1,NBTIMEDIA(KJ,1))
endif
  RETURN
ELSE

!
! Expression du temps en numerique
!
  IF(INDV == 0)THEN

! Cas  _TO_  _BY_

    IF(INDTO /= 0)THEN
      YCAR = ADJUSTL(YCART(1:INDTO-1))
      NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
      CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
      XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      IF(INDBY /= 0)THEN
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:INDBY-1))
        NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
	XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDBY+4:ILENCART))
        NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
	XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      ELSE
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:ILENCART))
        NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
	XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      END IF
    ELSE

! Cas un seul temps en fin de chaine de car. HCARIN ou au milieu

      IF(ILENCART > 9)THEN
	print *,' PB ecriture temps '
	STOP
      ELSE
	YCAR = ADJUSTL(YCART(1:ILENCART))
	NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
	XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      END IF

    END IF

  ELSE

! Presence de virgules

    YCAR = ADJUSTL(YCART(1:INDV-1))
    NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
    CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
    XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
    DO J = 1,100
      INDVM=INDV
      INDV = 0
      YCAR(1:LEN(YCAR))=' '
      INDV = INDEX(YCART(INDVM+1:ILENCART),',')
!     print *,' INDV ',INDV
      IF(INDV == 0)THEN
	YCAR = ADJUSTL(YCART(INDVM+1:ILENCART))
	NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
	XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
	EXIT
      ELSE
        INDV=INDV+INDVM
	YCAR = ADJUSTL(YCART(INDVM+1:INDV-1))
	NBTIMEDIA(KJ,:) = NBTIMEDIA(KJ,:)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1))
	XTIMEDIA(NBTIMEDIA(KJ,:),KJ,:)=XTIMEDIA(NBTIMEDIA(KJ,1),KJ,1)
      END IF
    ENDDO


  END IF
!
END IF
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
if (nverbia >0) then
print *,' end of RESOLVT '
print *,' LTIMEDIALL(KJ,1) ',LTIMEDIALL(KJ,1)
print *,' NBTIMEDIA(KJ,1) ',NBTIMEDIA(KJ,1)
print *,' NTIMEDIA ',(NTIMEDIA(J,KJ,1),J=1,NBTIMEDIA(KJ,1))
print *,' XTIMEDIA ',(XTIMEDIA(J,KJ,1),J=1,NBTIMEDIA(KJ,1))
endif
RETURN
END SUBROUTINE RESOLVT  
!     ######spl
      MODULE MODI_RESOLVL
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVL(HCARIN,K,OLOGIC)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: K
LOGICAL         :: OLOGIC
END SUBROUTINE RESOLVL
!
END INTERFACE
!
END MODULE MODI_RESOLVL
!     ###################################
      SUBROUTINE RESOLVL(HCARIN,K,OLOGIC)
!     ###################################
!
!!****  *RESOLVL* - 
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
USE MODD_RESOLVCAR
USE MODN_NCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: K
LOGICAL         :: OLOGIC
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=3) :: YC3
INTEGER          :: ILENC
INTEGER          :: J,JM, I

!
!------------------------------------------------------------------------------
ILENC=LEN_TRIM(HCARIN)

if(nverbia >0)then
  print *,' HCARIN K RESOLVL ',HCARIN,K
endif
DO J=K,ILENC
  IF(HCARIN(J:J) == '=')EXIT
ENDDO

JM=J+1
YC3='   '
I=0

if(nverbia >0)then
print *,' RESOLVL JM,ILENC ',JM,ILENC
endif
DO J=JM,ILENC
  IF(HCARIN(J:J) == 'T'.OR.HCARIN(J:J) == 'F')THEN
    YC3(1:1)=HCARIN(J:J)
    I=1
    EXIT
  ENDIF
ENDDO

IF(I == 0)THEN
print *,' PB AVEC LA VALEUR FOURNIE DE ',HCARIN(1:JM-2),'  ', &
 HCARIN(1:LEN_TRIM(HCARIN)),' VERIFIEZ SA VALEUR '
  RETURN
ENDIF
if(nverbia >0)then
  print *,' RESOLVL YC3 ',YC3
endif
IF(I == 1)READ(YC3(1:1),'(L1)')OLOGIC
print *,HCARIN(1:JM-2),' FOURNI ',OLOGIC
RETURN
END SUBROUTINE RESOLVL  
!     ######spl
      MODULE MODI_RESOLVZ
!     ###################
!
INTERFACE
!
SUBROUTINE RESOLVZ(HCARIN,KINDZ,KJ)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDZ, KJ
END SUBROUTINE RESOLVZ
!
END INTERFACE
!
END MODULE MODI_RESOLVZ
!     ###################################
      SUBROUTINE RESOLVZ(HCARIN,KINDZ,KJ)
!     ###################################
!
!!****  *RESOLVZ* - 
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
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
INTEGER          :: KINDZ, KJ
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=80) :: YCART
CHARACTER(LEN=20) :: YCAR
INTEGER          :: ILENC, ILENCART
INTEGER          :: INDTF, INDTO, INDBY, INDV, INDVM
INTEGER          :: ICAS, J

!
!------------------------------------------------------------------------------
INDTF = 0
INDTO = 0
INDBY = 0
INDV  = 0
ICAS = 0

NBLVLZDIA(KJ)=0
NLVLZDIA(:,KJ)=0
XLVLZDIA(:,KJ)=0.
LZINCRDIA(KJ)=.FALSE.

ILENC = LEN(HCARIN)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INDTO = INDEX(HCARIN(KINDZ+3:ILENC),'_TO_')
  INDBY = INDEX(HCARIN(KINDZ+3:ILENC),'_BY_')
  INDTF = INDEX(HCARIN(KINDZ+3:ILENC),'_')
  IF(INDTO /= 0)THEN
  IF(INDTF < INDTO)THEN
!
! ICAS = 1  Niveau Z unique ou plusieurs separes par des virgules
!
    INDTO=0;INDBY=0
    ICAS = 1
  ELSE IF(INDTF == INDTO)THEN
!
! ICAS = 3  NivZ1 _TO_ NivZn _BY_ NivZx
!
    IF(INDBY /= 0)THEN
      DO J=INDTO+4+KINDZ+3,INDBY+KINDZ+3
        IF(HCARIN(J:J) == '_')THEN
          IF(HCARIN(J:J+3) == '_BY_')THEN
            EXIT
          ELSE
            INDBY=0
            EXIT
          END IF
        END IF
      ENDDO
    END IF
    IF(INDBY /= 0)THEN
      INDTF=INDEX(HCARIN(KINDZ+3+INDBY+4:ILENC),'_')
      IF(INDTF /= 0)INDTF=INDTF+INDBY+4
      ICAS = 3
      LZINCRDIA(KJ) = .TRUE.
    ELSE
!
! ICAS = 2  NivZ1 _TO_ NivZn
!
      INDTF=INDEX(HCARIN(KINDZ+3+INDTO+4:ILENC),'_')
      IF(INDTF /= 0)INDTF=INDTF+INDTO+4
      ICAS = 2
      LZINCRDIA(KJ) = .TRUE.
    END IF
  END IF
  ELSE
    ICAS = 1
  END IF
IF(INDTF == 0)THEN
  INDTF = ILENC
ELSE
  INDTF = INDTF+KINDZ+3-1-1
END IF


YCART(1:LEN(YCART))=' '
YCAR(1:LEN(YCAR))=' '
!
! Extraction de la partie Niveaux Z dans YCART(1:ILENCART)
!
!print *,' KINDZ INDTF ',KINDZ,INDTF
YCART = ADJUSTL(HCARIN(KINDZ+3:INDTF))
ILENCART = LEN_TRIM(YCART)
!print *,' YCART ',ILENCART,' ',YCART

! Recherche a nouveau des chaines de car. _TO_ , _BY_ et d'une virgule
! par rapport au debut de YCART

INDTO = INDEX(YCART,'_TO_')
INDBY = INDEX(YCART,'_BY_')
INDV = INDEX(YCART(1:ILENCART),',')
IF(ICAS == 1 .AND. INDV == 0)ICAS=0
!
! Expression des niveaux Z par mots-cles (LVLZ1....)
!
IF(YCART(1:4) == 'LVLZ')THEN
!print *,' YCART(1:4) ',YCART(1:4),' ICAS ',ICAS

  NBLVLZDIA(KJ)=NBLVLZDIA(KJ)+1
  SELECT CASE(ICAS)
    CASE(1)
!print *,' INDV YCART(5:5) ',INDV,YCART(5:5)
      IF(INDV-4-1 == 1)READ(YCART(5:5),'(I1)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
      IF(INDV-4-1 == 2)READ(YCART(5:6),'(I2)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
      DO J = 1,100
        INDVM=INDV
        INDV=0
        INDV=INDEX(YCART(INDVM+1:ILENCART),',')
        IF(INDV == 0)THEN
          NBLVLZDIA(KJ)=NBLVLZDIA(KJ)+1
          IF(ILENCART-(INDVM+4) == 1)READ(YCART(INDVM+4+1:INDVM+4+1),'(I1)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
          IF(ILENCART-(INDVM+4) == 2)READ(YCART(INDVM+4+1:ILENCART),'(I2)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
          EXIT
        ELSE
          INDV=INDV+INDVM
          NBLVLZDIA(KJ)=NBLVLZDIA(KJ)+1
          IF(INDV-(INDVM+4)-1 == 1)READ(YCART(INDVM+4+1:INDVM+4+1),'(I1)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
          IF(INDV-(INDVM+4)-1 == 2)READ(YCART(INDVM+4+1:INDV-1),'(I2)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
        END IF
      ENDDO   
      
    CASE(2)
      IF(INDTO-4-1 == 1)READ(YCART(5:5),'(I1)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
      IF(INDTO-4-1 == 2)READ(YCART(5:6),'(I2)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
      NBLVLZDIA(KJ)=NBLVLZDIA(KJ)+1
      IF(ILENCART-(INDTO+3+4) == 1)READ(YCART(INDTO+3+4+1:INDTO+3+4+1),'(I1)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
      IF(ILENCART-(INDTO+3+4) == 2)READ(YCART(INDTO+3+4+1:ILENCART),'(I2)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
! 1 seul temps
    CASE DEFAULT
      IF(ILENCART-4 == 1)READ(YCART(5:5),'(I1)')NLVLZDIA(NBLVLZDIA(KJ),KJ)
      IF(ILENCART-4 == 2)READ(YCART(5:6),'(I2)')NLVLZDIA(NBLVLZDIA(KJ),KJ)

  END SELECT
  print *,' RESOLVZ ICAS '
  print *,' NBLVLZDIA ',NBLVLZDIA(KJ)
  print *,' NLVLZDIA ',(NLVLZDIA(J,KJ),J=1,NBLVLZDIA(KJ))
  print *,' XLVLZDIA ',(XLVLZDIA(J,KJ),J=1,NBLVLZDIA(KJ))
  RETURN
ELSE

!
! Expression des niveaux Z en numerique
!
  IF(INDV == 0)THEN

! Cas  _TO_  _BY_

    IF(INDTO /= 0)THEN
      YCAR = ADJUSTL(YCART(1:INDTO-1))
      NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
      CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
      IF(INDBY /= 0)THEN
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:INDBY-1))
        NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDBY+4:ILENCART))
        NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
      ELSE
        YCAR(1:LEN(YCAR))=' '
        YCAR = ADJUSTL(YCART(INDTO+4:ILENCART))
        NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
      END IF
    ELSE

! Cas un seul niveau Z en fin de chaine de car. HCARIN ou au milieu

      IF(ILENCART > 9)THEN
	print *,' PB ecriture temps '
	STOP
      ELSE
	YCAR = ADJUSTL(YCART(1:ILENCART))
	NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
      END IF

    END IF

  ELSE

! Presence de virgules

    YCAR = ADJUSTL(YCART(1:INDV-1))
    NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
    CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
    DO J = 1,100
      INDVM=INDV
      INDV = 0
      YCAR(1:LEN(YCAR))=' '
      INDV = INDEX(YCART(INDVM+1:ILENCART),',')
!     print *,' INDV ',INDV
      IF(INDV == 0)THEN
	YCAR = ADJUSTL(YCART(INDVM+1:ILENCART))
	NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
	EXIT
      ELSE
        INDV=INDV+INDVM
	YCAR = ADJUSTL(YCART(INDVM+1:INDV-1))
	NBLVLZDIA(KJ) = NBLVLZDIA(KJ)+1
        CALL CAREAL(YCAR(1:LEN_TRIM(YCAR)),XLVLZDIA(NBLVLZDIA(KJ),KJ))
      END IF
    ENDDO


  END IF
!
END IF
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
print *,' RESOLVZ '
print *,' NBLVLZDIA ',NBLVLZDIA(KJ)
print *,' NLVLZDIA ',(NLVLZDIA(J,KJ),J=1,NBLVLZDIA(KJ))
print *,' XLVLZDIA ',(XLVLZDIA(J,KJ),J=1,NBLVLZDIA(KJ))
RETURN
END SUBROUTINE RESOLVZ  
