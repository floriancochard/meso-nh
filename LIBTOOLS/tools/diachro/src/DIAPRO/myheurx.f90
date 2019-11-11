!     ######spl
      SUBROUTINE MYHEURX(KITVXJ,KITVXN,KITVYJ,KITVYN,I1,I2,I3,Z1,Z2)
!     ####################
!
!!****  *MYHEURX* - 
!!
!!    PURPOSE
!!    -------
!       
!     
!
!!**  METHOD
!!    ------
!!      NCAR routines are called to select a display window 
!!    corresponding to the post-processed section of the model 
!!    arrays (NIINFxNISUP).(NJINFxNJSUP)
!!     
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
!!      Module MODD_COORD  : declares gridpoint coordinates (TRACE use)
!!       XXX,XXY  : coordinate values for all the MESO-NH grids
!!
!!      Module MODD_NMGRID  : declares global variable  NMGRID
!!         NMGRID      : Current MESO-NH grid indicator
!!
!!      Module MODD_DIM1 : contains dimensions of data arrays
!!         NIINF, NISUP : lower and upper bounds of arrays
!!                        to be plotted in x direction
!!         NJINF, NJSUP : lower and upper bounds of arrays
!!                        to be plotted in y direction
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
!!      Original       25/04/02
!!      Updated   PM   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
USE MODD_CTL_AXES_AND_STYL
USE MODD_DIM1
USE MODN_NCAR
!
IMPLICIT NONE
!
INTEGER :: KITVXJ,KITVXN,KITVYJ,KITVYN,I1,I2,I3
REAL    :: Z1,Z2
!

REAL :: ZWL, ZWR, ZWB, ZWT
REAL :: ZWLL, ZWRR, ZWBB, ZWTT
REAL :: ZVL, ZVR, ZVB, ZVT
REAL :: ZH, ZJ, ZJJ,ZINT, ZINTT, ZWBBB
INTEGER :: ID, IDD ,J
CHARACTER(LEN=2)  :: YC2
CHARACTER(LEN=3)  :: YC3
CHARACTER(LEN=4)  :: YC4
CHARACTER(LEN=10)  :: FORMAX, FORMAY
!
!-------------------------------------------------------------------------------
!
!*       1.    DISPLAY WINDOW SETTING AND DRAWING
!              ----------------------------------
!
!-----------------------------------------------------------------------------
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
!CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL/3600.,ZWR/3600.,ZWB,ZWT,ID)

!!!!!!!Avril 2002
  IF(LMYHEURX)THEN
    ZH=NHEURXGRAD*3600.
  ELSE
!!!!!!!Avril 2002

  IF((ZWR-ZWL)/3600. > 24.)THEN
    ZH=10800.
  ELSE
    ZH=3600.
  ENDIF
!!!!!!!Avril 2002
  ENDIF
!!!!!!!Avril 2002

  DO J=INT(ZWL),INT(ZWR)
    ZJ=J
!     print *,' ZJ, ',ZJ
    IF(MOD(ZJ,ZH) == 0.)THEN
!     print *,' ZJ,ZH,ZWB,ZWT ',ZJ,ZH,ZWB,ZWT
      IF(I1 /= -1 .AND. I1 /= 0)THEN
      
      CALL FRSTPT(ZJ,ZWB)
      CALL VECTOR(ZJ,ZWB+(ZWT-ZWB)/90.)
      CALL FRSTPT(ZJ,ZWT)
      CALL VECTOR(ZJ,ZWT-(ZWT-ZWB)/90.)
      
      ENDIF
!!!!!!!Avril 2002
  IF(LMYHEURX)THEN
    ZJJ=ZJ/ZH*NHEURXGRAD
    ZINTT=NHEURXLBL
  ELSE
!!!!!!!Avril 2002


      IF(ZH == 10800.)THEN
        ZJJ=ZJ/ZH*3.
        ZINTT=6.
      ELSE
        ZJJ=ZJ/ZH
        ZINTT=3.
      ENDIF
   !!!!!!!Avril 2002
  ENDIF
!!!!!!!Avril 2002

      CALL GSCLIP(0)
      ZWBBB=ZWB-((ZWT-ZWB)/40)
!     print *,' ZWB ZWT ZWBBB ',ZWB,ZWT,ZWBBB
      

      IF(I1 == 1 .AND. .NOT.LNOLABELX)THEN
      IF(MOD(ZJJ,ZINTT) == 0.)THEN
        IF(LFACTAXEX)THEN
          ZJJ=ZJJ*XFACTAXEX
        ENDIF
        IF(ZJJ < 1.)THEN
          YC4='    '
          WRITE(YC4,'(F4.2)')ZJJ
          CALL PLCHHQ(ZJ,ZWBBB,YC4,.010,0.,0.)

        ELSEIF(ZJJ < 10.)THEN
          YC2='  '
          WRITE(YC2,'(F2.0)')ZJJ
          CALL PLCHHQ(ZJ,ZWBBB,YC2,.010,0.,0.)
        ELSEIF(ZJJ < 100.)THEN
          YC3='   '
          WRITE(YC3,'(F3.0)')ZJJ
          CALL PLCHHQ(ZJ,ZWBBB,YC3,.010,0.,0.)
        ELSE
          YC4='    '
          WRITE(YC4,'(F4.0)')ZJJ
          CALL PLCHHQ(ZJ,ZWBBB,YC4,.010,0.,0.)
        ENDIF
      ENDIF
      ENDIF

    ENDIF
ENDDO
!!! Inutile IMPLEMENTE SEULEMENT EN CV
 CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
 print *,'**myheurx ZWLL,ZWRR,ZWBB,ZWTT ',ZWLL,ZWRR,ZWBB,ZWTT
   IF(LFACTAXEX)THEN
     IF(LFACTAXEY)THEN
       CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL*XFACTAXEX,ZWRR*XFACTAXEX,&
                ZWBB*XFACTAXEY,ZWTT*XFACTAXEY,IDD)
     ELSE
       CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL*XFACTAXEX,ZWRR*XFACTAXEX,&
                ZWBB,ZWTT,IDD)
     ENDIF
   ELSEIF(LFACTAXEY)THEN
       CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
                ZWBB*XFACTAXEY,ZWTT*XFACTAXEY,IDD)
   ELSEIF(LAXEXUSER)THEN
     IF(LAXEYUSER)THEN
       CALL SET(ZVL,ZVR,ZVB,ZVT,XAXEXUSERD,XAXEXUSERF,&
                XAXEYUSERD,XAXEYUSERF,IDD)
     ELSE
       CALL SET(ZVL,ZVR,ZVB,ZVT,XAXEXUSERD,XAXEXUSERF,&
                ZWBB,ZWTT,IDD)
     ENDIF
   ELSEIF(LAXEYUSER)THEN
       CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,&
                XAXEYUSERD,XAXEYUSERF,IDD)
   ENDIF
!!! Inutile IMPLEMENTE SEULEMENT EN CV
! Mars 2001


! Mars 2001
 print *,'**myheurx ZWLL,ZWRR,ZWBB,ZWTT ',ZWLL,ZWRR,ZWBB,ZWTT
!CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
 CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWLL,ZWRR,ZWBB,ZWTT,IDD)
 print *,'**myheurx ZWLL,ZWRR,ZWBB,ZWTT ',ZWLL,ZWRR,ZWBB,ZWTT
!IF(LFACTAXEX)THEN
!CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL/3600.*XFACTAXEX,ZWRR/3600.*XFACTAXEX,ZWBB,ZWTT,IDD)
!ELSE

 CALL SET(ZVL,ZVR,ZVB,ZVT,ZWLL/3600.,ZWRR/3600.,ZWBB,ZWTT,IDD)
 print *,'**myheurx ZWLL/3600,ZWRR/3600,ZWBB,ZWTT ',ZWLL/3600,ZWRR/3600,ZWBB,ZWTT
!ENDIF
!Avril 2002
    IF(LNOLABELX .AND. LNOLABELY)THEN
      IF(I1 /= -1)THEN
      CALL GRIDAL(0,0,KITVYJ,KITVYN,0,0,I3,Z1,Z2)
      ELSE
      CALL GRIDAL(0,0,KITVYJ,KITVYN,I1,0,I3,Z1,Z2)
      ENDIF
!     CALL GRIDAL(0,0,KITVYJ,KITVYN,I1,I2,I3,Z1,Z2)
    ELSEIF(LNOLABELX .AND. .NOT.LNOLABELY)THEN
      IF(I1 /= -1)THEN
      CALL GRIDAL(0,0,KITVYJ,KITVYN,0,I2,I3,Z1,Z2)
      ELSE
      CALL GRIDAL(0,0,KITVYJ,KITVYN,I1,I2,I3,Z1,Z2)
      ENDIF
!     CALL GRIDAL(0,0,KITVYJ,KITVYN,I1,I2,I3,Z1,Z2)
    ELSEIF(.NOT.LNOLABELX .AND. LNOLABELY)THEN
      IF(I1 /= -1)THEN
        CALL GRIDAL(0,0,KITVYJ,KITVYN,0,0,I3,Z1,Z2)
      ELSE
        CALL GRIDAL(0,0,KITVYJ,KITVYN,I1,0,I3,Z1,Z2)
      ENDIF
!     CALL GRIDAL(0,0,KITVYJ,KITVYN,I1,I2,I3,Z1,Z2)
    ELSE
      IF(I1 == 1)THEN
      CALL GRIDAL(0,0,KITVYJ,KITVYN,0,I2,I3,Z1,Z2)
      ELSE
      CALL GRIDAL(0,0,KITVYJ,KITVYN,I1,I2,I3,Z1,Z2)
      ENDIF
    ENDIF
!Avril 2002
    CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
! ENDIF
      CALL GSCLIP(1)

!
!*      2.   EXIT
!            ----
!
RETURN
END SUBROUTINE  MYHEURX
