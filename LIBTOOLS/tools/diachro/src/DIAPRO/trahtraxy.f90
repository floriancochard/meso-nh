!     #################
      SUBROUTINE TRAHTRAXY(KLOOP,PTEMCV,HTEXTE)
!     #################
!
!!****  *TRAHTRAXY* - 
!!                                                            
!!
!!    PURPOSE
!!    -------
!        Trace PH (tableaux 1D scalaires  y compris MUMVM et DIRUMVM)
!        dans traceh_fordiachro
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
!!      Original       30/11/01
!!      Updated   PM  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!

USE MODD_NMGRID
USE MODD_COORD
USE MODD_DEFCV
USE MODD_TIT  
USE MODD_TYPE_AND_LH
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO
USE MODN_PARA
USE MODN_NCAR
USE MODI_RESOLV_TIT
USE MODI_RESOLV_TITY

IMPLICIT NONE
!
INTERFACE
      SUBROUTINE TRAXY(PTEMX,PTEMY,KLOOP,HTITX,HTITY,PTIMED,PTIMEF)
      INTEGER    :: KLOOP
      REAL,DIMENSION(:)  :: PTEMX, PTEMY
      REAL               :: PTIMED, PTIMEF
      CHARACTER(LEN=*) :: HTITX, HTITY
      END SUBROUTINE TRAXY
END INTERFACE
!
!
!*       0.1   Dummy arguments
!
INTEGER           :: KLOOP
REAL,DIMENSION(:,:)         :: PTEMCV
CHARACTER(LEN=40) :: HTEXTE
!
!*       0.1   Local variables
!
!
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZTEMCV
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZTEMX, ZTEMY
REAL              :: ZTIMED, ZTIMEF
REAL              :: ZXPOSTITT1, ZXYPOSTITT1
REAL              :: ZXPOSTITT2, ZXYPOSTITT2
REAL              :: ZXPOSTITT3, ZXYPOSTITT3
REAL              :: ZXPOSTITT4, ZXYPOSTITT4
REAL              :: ZXPOSTITB1, ZXYPOSTITB1
REAL              :: ZXPOSTITB2, ZXYPOSTITB2
REAL              :: ZXPOSTITB3, ZXYPOSTITB3
REAL              :: ZXPOSTITB4, ZXYPOSTITB4
!
CHARACTER(LEN=16) :: YTITX,YTITY
CHARACTER(LEN=40) :: YTEXTE,YTEM
CHARACTER(LEN=80) :: YCARCOU
!
!-------------------------------------------------------------------------------
!
!*      1. 
!              ----------------------------
!
YTEXTE=HTEXTE
!!!!!!!!!!!!! Supprime le 30/11/01
! Appel a TRAXY pour le trace du PH
       IF(ALLOCATED(ZTEMX))THEN
	 DEALLOCATE(ZTEMX)
       ENDIF
       IF(ALLOCATED(ZTEMY))THEN
	 DEALLOCATE(ZTEMY)
       ENDIF
       IF(ALLOCATED(ZTEMCV))THEN
	 DEALLOCATE(ZTEMCV)
       ENDIF
       ALLOCATE(ZTEMCV(SIZE(PTEMCV,1),SIZE(PTEMCV,2)))
       ZTEMCV(:,:)=PTEMCV(:,:)
       ALLOCATE(ZTEMX(SIZE(ZTEMCV,1)))
       ALLOCATE(ZTEMY(SIZE(ZTEMCV,1)))
       IF(SIZE(ZTEMCV,2) == 1)THEN
         ZTEMY(:)=ZTEMCV(:,1)
       ELSE
         ZTEMY(:)=ZTEMCV(:,MAX(2,NKL))
       ENDIF
       ZTEMX(:)=XDS(1:NLMAX,NMGRID)
        WHERE(ZTEMY == XSPVAL)
	  ZTEMY=1.E36
        END WHERE
       YTITX(1:LEN(YTITX))=' '
       YTITY(1:LEN(YTITX))=' '
       YTITX='X(M)'
       YTITY=CUNITGAL(1:LEN(CUNITGAL))
       ZTIMED=XTRAJT(NLOOPT,1)
       ZTIMEF=ZTIMED
       IF(NVERBIA > 0)THEN
	 print *,' TRACEH AV TRAXY KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF',&
	 KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF
       ENDIF
       CALL TRAXY(ZTEMX,ZTEMY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)

        IF(KLOOP == 1)THEN

	  IF(LDATFILE)CALL DATFILE_FORDIACHRO
	  CALL RESOLV_TIMES(NLOOPT)
	  YTEM(1:LEN(YTEM))=' '
! CTITVAR1
	  CALL RESOLV_TIT('CTITVAR1',YTEM)
	  IF(CTITVAR1 == 'DEFAULT')THEN
	    CALL PLCHHQ(.99,.007,YTEXTE(1:LEN_TRIM(YTEXTE)),.011,0.,+1.)
          ELSE IF(YTEM /= ' ')THEN
	    CALL PLCHHQ(.99,.007,YTEM(1:LEN_TRIM(YTEM)),.011,0.,+1.)
	  ENDIF
! CTITT1
	  YCARCOU(1:LEN(YCARCOU))=' '
	  YTEM(1:LEN(YTEM))=' '
	  CALL RESOLV_TIT('CTITT1',YTEM)
	  ZXPOSTITT1=.002
          ZXYPOSTITT1=.98
          IF(XPOSTITT1 /= 0.)THEN
            ZXPOSTITT1=XPOSTITT1
	  ENDIF
	  IF(XYPOSTITT1 /= 0.)THEN
	    ZXYPOSTITT1=XYPOSTITT1
	  ENDIF

          IF(XIDEBCOU.NE.-999.)THEN
            IF(LDEFCV2CC)THEN           !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      IF(LDEFCV2IND)THEN
	        WRITE(YCARCOU,1018)NIDEBCV,NJDEBCV,NIFINCV,NJFINCV
	      ELSE IF(LDEFCV2LL)THEN
	        WRITE(YCARCOU,1019)XIDEBCVLL,XJDEBCVLL,XIFINCVLL,XJFINCVLL
	      ELSE
	        WRITE(YCARCOU,1020)XIDEBCV,XJDEBCV,XIFINCV,XJFINCV
	      ENDIF
            ELSE                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              IF(XIDEBCOU < 99999.)THEN
                IF(XJDEBCOU < 99999.)THEN
                  WRITE(YCARCOU,1011)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
                ELSE
                  WRITE(YCARCOU,1013)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
                END IF
              ELSE
                IF(XJDEBCOU < 99999.)THEN
                  WRITE(YCARCOU,1014)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
                ELSE
                  WRITE(YCARCOU,1015)XIDEBCOU,XJDEBCOU,NLANGLE,NLMAX
                END IF
              END IF
            ENDIF                       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          ELSE
            WRITE(YCARCOU,1010)NIDEBCOU,NJDEBCOU,NLANGLE,NLMAX
          ENDIF
	  IF(CTITT1 == 'DEFAULT')THEN
	    IF(XSZTITT1 /= 0.)THEN
	      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU(1:LEN_TRIM(YCARCOU)),XSZTITT1,0.,-1.)
!             CALL PLCHHQ(.002,.98,YCARCOU(1:LEN_TRIM(YCARCOU)),XSZTITT1,0.,-1.)
	    ELSE
	      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YCARCOU(1:LEN_TRIM(YCARCOU)),.012,0.,-1.)
!             CALL PLCHHQ(.002,.98,YCARCOU(1:LEN_TRIM(YCARCOU)),.012,0.,-1.)
	    ENDIF
          ELSE IF(YTEM /= ' ')THEN
	    IF(XSZTITT1 /= 0.)THEN
	      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM(1:LEN_TRIM(YTEM)),XSZTITT1,0.,-1.)
!             CALL PLCHHQ(.002,.98,YTEM(1:LEN_TRIM(YTEM)),XSZTITT1,0.,-1.)
	    ELSE
	      CALL PLCHHQ(ZXPOSTITT1,ZXYPOSTITT1,YTEM(1:LEN_TRIM(YTEM)),.012,0.,-1.)
!             CALL PLCHHQ(.002,.98,YTEM(1:LEN_TRIM(YTEM)),.012,0.,-1.)
	    ENDIF
	  ENDIF
! CTITT2
	  YTEM(1:LEN(YTEM))=' '
	  CALL RESOLV_TIT('CTITT2',YTEM)
	  ZXPOSTITT2=.002
          ZXYPOSTITT2=.95
          IF(XPOSTITT2 /= 0.)THEN
            ZXPOSTITT2=XPOSTITT2
	  ENDIF
	  IF(XYPOSTITT2 /= 0.)THEN
	    ZXYPOSTITT2=XYPOSTITT2
	  ENDIF
          IF(YTEM /= ' ')THEN
	    IF(XSZTITT2 /= 0.)THEN
	      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM(1:LEN_TRIM(YTEM)),XSZTITT2,0.,-1.)
!             CALL PLCHHQ(.002,.95,YTEM(1:LEN_TRIM(YTEM)),XSZTITT2,0.,-1.)
	    ELSE
	      CALL PLCHHQ(ZXPOSTITT2,ZXYPOSTITT2,YTEM(1:LEN_TRIM(YTEM)),.008,0.,-1.)
!             CALL PLCHHQ(.002,.95,YTEM(1:LEN_TRIM(YTEM)),.008,0.,-1.)
	    ENDIF
	  ENDIF
! CTITT3
	  YTEM(1:LEN(YTEM))=' '
	  CALL RESOLV_TIT('CTITT3',YTEM)
	  ZXPOSTITT3=.002
          ZXYPOSTITT3=.93
          IF(XPOSTITT3 /= 0.)THEN
            ZXPOSTITT3=XPOSTITT3
	  ENDIF
	  IF(XYPOSTITT3 /= 0.)THEN
	    ZXYPOSTITT3=XYPOSTITT3
	  ENDIF
          IF(YTEM /= ' ')THEN
	    IF(XSZTITT3 /= 0.)THEN
	      CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM(1:LEN_TRIM(YTEM)),XSZTITT3,0.,-1.)
!             CALL PLCHHQ(.002,.93,YTEM(1:LEN_TRIM(YTEM)),XSZTITT3,0.,-1.)
	    ELSE
	      CALL PLCHHQ(ZXPOSTITT3,ZXYPOSTITT3,YTEM(1:LEN_TRIM(YTEM)),.008,0.,-1.)
!             CALL PLCHHQ(.002,.93,YTEM(1:LEN_TRIM(YTEM)),.008,0.,-1.)
	    ENDIF
	  ENDIF

        ENDIF
!!!!!!!!!!!!! Supprime le 30/11/01
1010 FORMAT('Horiz. profile IDEB=',I3,' JDEB=',I3,' ANG.=',I3,' NBPTS=',I3)
1011 FORMAT('Horiz. profile XDEB=',F6.0,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I3)
1013 FORMAT('Horiz. profile XDEB=',F6.0,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I3)
1014 FORMAT('Horiz. profile XDEB=',E7.2,' YDEB=',F6.0,' ANG.=',I3,' NBPTS=',I3)
1015 FORMAT('Horiz. profile XDEB=',E6.2,' YDEB=',E7.2,' ANG.=',I3,' NBPTS=',I3)
1018 FORMAT('Horiz. profile IND I,J (BEGIN)-(END)=(',I3,',',I3,')-(',I3,',',I3,')')
1019 FORMAT('Horiz. profile LAT,LON (BEGIN)-(END)=(',F4.1,',',F5.1,')-(',F4.1,',',F5.1,')')
1020 FORMAT('Horiz. profile CONF. COORD.(BEGIN)-(END)=(',F8.0,',',F8.0,')-(',F8.0,',',F8.0,')')
!
!
!------------------------------------------------------------------------------
!
!*      2.    EXIT
!             ----
!
RETURN
END SUBROUTINE TRAHTRAXY
