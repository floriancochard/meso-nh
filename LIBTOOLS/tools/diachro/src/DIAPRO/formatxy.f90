!     ######spl
      SUBROUTINE FORMATXY(PWL,PWR,PWW,PWT)
!     ####################################
!
!!****  *FORMATXY* - 
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
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!  MODIFIED (I. Mallet 06/02) PWB to PWW
!!                            otherwise PWB is changed to 1 by cpp on HP
!!
!
USE MODD_RESOLVCAR
!
IMPLICIT NONE
!
!*      0.1    Dummy arguments 
!
REAL               :: PWL,PWR,PWW,PWT
!
!*      0.2    local variables 
!
!
INTEGER            :: J
CHARACTER(LEN=10)  :: FORMAX, FORMAY
CHARACTER(LEN=10)  :: YFORMAX
!
!-------------------------------------------------------------------------------
YFORMAX(1:LEN(YFORMAX))=' '
IF(NHISTORY(NLOOPSUPER) == 3)THEN
  DO J=1,MAX(1,NLOOPSUPER-1)
    IF(NHISTORY(J) == 1)THEN
      CALL GAGETC('XLF',YFORMAX)
      print *,' ** formatxy FORMAX ',YFORMAX
      EXIT
    ENDIF
  ENDDO
ENDIF
!
! PWR /= 0.
! ******************************************************************
  IF(PWR /= 0.)THEN
! ******************************************************************
  IF(LOG10(ABS(PWR)) >= 6. .OR. LOG10(ABS(PWR)) <= -1.)THEN

    FORMAX='          '
    IF(LFMTAXEX)THEN
      FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
    ELSE
      FORMAX='(E8.2)'
    ENDIF
  
! ------------------------------------------------------------------
! PWT /= 0.
    IF(PWT /= 0.)THEN
      FORMAY='          '
      IF(LOG10(ABS(PWT)) >= 6. .OR. LOG10(ABS(PWT)) <= -1.)THEN
        IF(LFMTAXEY)THEN
          FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
        ELSE
          FORMAY='(E8.2)'
        ENDIF
!       CALL LABMOD('(E8.2)','(E8.2)',0,0,10,10,0,0,0)
      ELSE
        IF(ABS(PWT-PWW) < 1.)THEN
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(F8.2)'
          ENDIF
!         CALL LABMOD('(E8.2)','(F8.2)',0,0,10,10,0,0,0)
        ELSE
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(F8.1)'
          ENDIF
!         CALL LABMOD('(E8.2)','(F8.1)',0,0,10,10,0,0,0)
        ENDIF
      ENDIF
      CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!     CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
! ------------------------------------------------------------------
    ELSE
! PWT == 0.
      FORMAY='          '
      IF(LOG10(ABS(PWW)) >= 6. .OR. LOG10(ABS(PWW)) <= -1.)THEN
        IF(LFMTAXEY)THEN
          FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
        ELSE
          FORMAY='(E8.2)'
        ENDIF
!       CALL LABMOD('(E8.2)','(E8.2)',0,0,10,10,0,0,0)
      ELSE
        IF(ABS(PWT-PWW) < 1.)THEN
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(F8.2)'
          ENDIF
!         CALL LABMOD('(E8.2)','(F8.2)',0,0,10,10,0,0,0)
        ELSE
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(F8.1)'
          ENDIF
!         CALL LABMOD('(E8.2)','(F8.1)',0,0,10,10,0,0,0)
        ENDIF
      ENDIF
      CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!     CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
    ENDIF
! ------------------------------------------------------------------
  
  ELSE
  
    IF(ABS(PWR-PWL) < 1.)THEN
      FORMAX='          '
      IF(LFMTAXEX)THEN
        FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
      ELSE
        FORMAX='(F8.2)'
      ENDIF
! ------------------------------------------------------------------
! PWT /= 0.
      IF(PWT /= 0.)THEN
        FORMAY='          '
        IF(LOG10(ABS(PWT)) >= 6. .OR. LOG10(ABS(PWT)) <= -1.)THEN
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(E8.2)'
          ENDIF
!         CALL LABMOD('(F8.2)','(E8.2)',0,0,10,10,0,0,0)
        ELSE
          IF(ABS(PWT-PWW) < 1.)THEN
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.2)'
            ENDIF
!           CALL LABMOD('(F8.2)','(F8.2)',0,0,10,10,0,0,0)
          ELSE
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.1)'
            ENDIF
!           CALL LABMOD('(F8.2)','(F8.1)',0,0,10,10,0,0,0)
          ENDIF
        ENDIF
        CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!       CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
! ------------------------------------------------------------------
      ELSE
! PWT == 0.
        FORMAY='          '
        IF(LOG10(ABS(PWW)) >= 6. .OR. LOG10(ABS(PWW)) <= -1.)THEN
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(E8.2)'
          ENDIF
!         CALL LABMOD('(F8.2)','(E8.2)',0,0,10,10,0,0,0)
        ELSE
          IF(ABS(PWT-PWW) < 1.)THEN
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.2)'
            ENDIF
!           CALL LABMOD('(F8.2)','(F8.2)',0,0,10,10,0,0,0)
          ELSE
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.1)'
            ENDIF
!           CALL LABMOD('(F8.2)','(F8.1)',0,0,10,10,0,0,0)
          ENDIF
        ENDIF
        CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!       CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
      ENDIF
! ------------------------------------------------------------------
  
    ELSE

      FORMAX='          '
      IF(LFMTAXEX)THEN
        FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
      ELSE
        FORMAX='(F8.1)'
      ENDIF
  
! ------------------------------------------------------------------
! PWT /= 0.
      IF(PWT /= 0.)THEN
        FORMAY='          '
        IF(LOG10(ABS(PWT)) >= 6. .OR. LOG10(ABS(PWT)) <= -1.)THEN
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(E8.2)'
          ENDIF
!         CALL LABMOD('(F8.1)','(E8.2)',0,0,10,10,0,0,0)
        ELSE
          IF(ABS(PWT-PWW) < 1.)THEN
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.2)'
            ENDIF
!           CALL LABMOD('(F8.1)','(F8.2)',0,0,10,10,0,0,0)
          ELSE
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.1)'
            ENDIF
!           CALL LABMOD('(F8.1)','(F8.1)',0,0,10,10,0,0,0)
          ENDIF
        ENDIF
        CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!       CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
  
! ------------------------------------------------------------------
      ELSE
! PWT == 0.
        FORMAY='          '
        IF(LOG10(ABS(PWW)) >= 6. .OR. LOG10(ABS(PWW)) <= -1.)THEN
          IF(LFMTAXEY)THEN
            FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
          ELSE
            FORMAY='(E8.2)'
          ENDIF
!         CALL LABMOD('(F8.1)','(E8.2)',0,0,10,10,0,0,0)
        ELSE
          IF(ABS(PWT-PWW) < 1.)THEN
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.2)'
            ENDIF
!           CALL LABMOD('(F8.1)','(F8.2)',0,0,10,10,0,0,0)
          ELSE
            IF(LFMTAXEY)THEN
              FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
            ELSE
              FORMAY='(F8.1)'
            ENDIF
!           CALL LABMOD('(F8.1)','(F8.1)',0,0,10,10,0,0,0)
          ENDIF
        ENDIF
        CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!       CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)
      ENDIF
! ------------------------------------------------------------------
    ENDIF
  
  ENDIF

! ******************************************************************
  ELSE
! ******************************************************************
! PWR = 0
          IF(LOG10(ABS(PWR-PWL)) >= 6. .OR. LOG10(ABS(PWR-PWL)) <= -1.)THEN

            FORMAX='          '
            IF(LFMTAXEX)THEN
              FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
            ELSE
              FORMAX='(E8.2)'
            ENDIF
            FORMAY='          '
	    IF(LOG10(ABS(PWT-PWW)) >= 6. .OR. LOG10(ABS(PWT-PWW)) <= -1.)THEN
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(E8.2)'
              ENDIF
!             CALL LABMOD('(E8.2)','(E8.2)',0,0,10,10,0,0,0)
	    ELSE IF(ABS(PWT-PWW) <1.)THEN
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(F8.2)'
              ENDIF
!             CALL LABMOD('(E8.2)','(F8.2)',0,0,10,10,0,0,0)
	    ELSE
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(F8.1)'
              ENDIF
!             CALL LABMOD('(E8.2)','(F8.1)',0,0,10,10,0,0,0)
	    ENDIF
            CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!           CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)

	  ELSE IF(ABS(PWR-PWL) < 1.)THEN

            FORMAX='          '
            IF(LFMTAXEX)THEN
              FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
            ELSE
              FORMAX='(F8.2)'
            ENDIF
            FORMAY='          '
	    IF(LOG10(ABS(PWT-PWW)) >= 6. .OR. LOG10(ABS(PWT-PWW)) <= -1.)THEN
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(E8.2)'
              ENDIF
!             CALL LABMOD('(F8.2)','(E8.2)',0,0,10,10,0,0,0)
	    ELSE IF(ABS(PWT-PWW) <1.)THEN
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(F8.2)'
              ENDIF
!             CALL LABMOD('(F8.2)','(F8.2)',0,0,10,10,0,0,0)
	    ELSE
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(F8.1)'
              ENDIF
!             CALL LABMOD('(F8.2)','(F8.1)',0,0,10,10,0,0,0)
	    ENDIF
            CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!           CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)

	  ELSE

            FORMAX='          '
            IF(LFMTAXEX)THEN
              FORMAX="("//CFMTAXEX(1:LEN_TRIM(CFMTAXEX))//")"
            ELSE
              FORMAX='(F8.1)'
            ENDIF
            FORMAY='          '
	    IF(LOG10(ABS(PWT-PWW)) >= 6. .OR. LOG10(ABS(PWT-PWW)) <= -1.)THEN
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(E8.2)'
              ENDIF
!             CALL LABMOD('(F8.1)','(E8.2)',0,0,10,10,0,0,0)
	    ELSE IF(ABS(PWT-PWW) <1.)THEN
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(F8.2)'
              ENDIF
!             CALL LABMOD('(F8.1)','(F8.2)',0,0,10,10,0,0,0)
	    ELSE
              IF(LFMTAXEY)THEN
                FORMAY="("//CFMTAXEY(1:LEN_TRIM(CFMTAXEY))//")"
              ELSE
                FORMAY='(F8.1)'
              ENDIF
!             CALL LABMOD('(F8.1)','(F8.1)',0,0,10,10,0,0,0)
	    ENDIF
            CALL LABMOD(FORMAX,FORMAY,0,0,NSZLBX,NSZLBY,0,0,0)
!           CALL LABMOD(FORMAX,FORMAY,0,0,10,10,0,0,0)

	  ENDIF

! ******************************************************************
  ENDIF
! ******************************************************************
! Prise en compte d'une superposition d'un PH=CV+K sur une CV pour des
! labels interieurs
IF(NHISTORY(NLOOPSUPER) == 3)THEN
  DO J=1,MAX(1,NLOOPSUPER-1)
    IF(NHISTORY(J) == 1)THEN
      CALL LABMOD(YFORMAX,FORMAY,0,0,NSZLBX,NSZLBY,-25,0,0)
!     CALL LABMOD(YFORMAX,FORMAY,0,0,10,10,-25,0,0)
      EXIT
    ENDIF
  ENDDO
ENDIF
!!----------------------------------------------------------------------------
RETURN
!
!*       4.     EXIT
!               ----
!
END SUBROUTINE  FORMATXY
