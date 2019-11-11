!     ######spl
      SUBROUTINE COLOR_FORDIACHRO(KN,KTYPE)
!     ###############################
!
!!****  *COLOR_FORDIACHRO* - Definition d'une table de couleurs en RGB
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     
!!    
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       16/01/95
!!      Updated   PM 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
IMPLICIT NONE
!
!*       0.1  dummy arguments
!          
INTEGER          :: KN
INTEGER          :: KTYPE   ! (1=color; 2=grey)
!
!*       0.2  local variables
!          
REAL :: ZHUE, ZHUES, ZL, ZLS
REAL :: ZRED, ZGREEN, ZBLUE


INTEGER          :: J, JJ
INTEGER          :: ICNT, ISTA, IER, INB, IWK
INTEGER          :: INBB
!
!-------------------------------------------------------------------------------
CALL GQOPS(ISTA)
CALL GQACWK(1,IER,INB,IWK)
CALL GQOPWK(1,IER,INB,IWK)
IF(KN <= 0)RETURN
!
IF(LINVWB)THEN
  CALL GSCR(IWK,0,1.,1.,1.) ! BACKGROUND COLOR (black)
  CALL GSCR(IWK,1,0.,0.,0.) ! First foreground color is white
ELSE
  CALL GSCR(IWK,0,0.,0.,0.) ! BACKGROUND COLOR (black)
  CALL GSCR(IWK,1,1.,1.,1.) ! First foreground color is white
ENDIF
!
DO JJ=1,INB
  CALL GQOPWK(JJ,IER,INBB,IWK)
  IF(IWK == 9)THEN
    CYCLE
  ELSE
    CALL GSCR(IWK,2,.75,.75,.75) ! Second foreground color is gray
  ENDIF
ENDDO
!
! Choose other foreground colors spaced equally around the spectrum
ICNT=0
ZHUES=360./KN
ZLS=100./KN
DO J=1,KN
  ZHUE=J*ZHUES
  ZL=J*ZLS
  IF(KTYPE==1) THEN
    !full colors
    CALL HLSRGB(ZHUE,50.,100.,ZRED,ZGREEN,ZBLUE)
!   CALL HLSRGB(ZHUE,55.,95.,ZRED,ZGREEN,ZBLUE)
!   CALL HLSRGB(ZHUE,60.,75.,ZRED,ZGREEN,ZBLUE)
    IF(ZHUE.LE.36.)THEN
      DO JJ=1,INB
        CALL GQOPWK(JJ,IER,INBB,IWK)
        IF(IWK == 9)THEN
          CYCLE
        ELSE
          CALL GSCR(IWK,KN+3-J,ZRED,ZGREEN,ZBLUE)
        ENDIF
      ENDDO
      ICNT=ICNT+1
    ELSE
      DO JJ=1,INB
        CALL GQOPWK(JJ,IER,INBB,IWK)
        IF(IWK == 9)THEN
          CYCLE
        ELSE
          CALL GSCR(IWK,J-ICNT+2,ZRED,ZGREEN,ZBLUE)
        ENDIF
      ENDDO
    END IF
  ELSE IF(KTYPE==2) THEN
    !greys (S=0.)
    CALL HLSRGB(ZHUE,ZL,0.,ZRED,ZGREEN,ZBLUE)
    DO JJ=1,INB
      CALL GQOPWK(JJ,IER,INBB,IWK)
      IF(IWK == 9)THEN
        CYCLE
      ELSE
        CALL GSCR(IWK,J-ICNT+2,ZRED,ZGREEN,ZBLUE)
      ENDIF
    ENDDO
END IF
ENDDO
!
RETURN
END SUBROUTINE COLOR_FORDIACHRO
