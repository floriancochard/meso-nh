!     ######spl
      SUBROUTINE DATFILE_FORDIACHRO
!     #############################
!
!!****  *DATFILE_FORDIACHRO* - Recupere la date du run du graphique et l'inscrit sur
!                   le dessin ainsi que le nom du fichier traite
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
!!      Original       19/09/95
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_OUT
USE MODD_FILES_DIACHRO
USE MODD_RESOLVCAR
USE MODD_TYPE_AND_LH
USE MODD_ALLOC_FORDIACHRO
!
IMPLICIT NONE
!
!*       0.1  dummy argument
!
!          
!
!
!*       0.1  local variables
!          
!
CHARACTER(LEN=8) :: YTIM8, YTEM8
CHARACTER(LEN=9) :: YTEM9
#if defined(HPPA)
CHARACTER(LEN=9) :: YDAT8
#else
#if defined(LINUX) || defined (O2000) 
CHARACTER(LEN=9) :: YDAT8
CHARACTER(LEN=10) :: YTIM10
#else
#if defined(VPP)
CHARACTER(LEN=8) :: YDAT8
#endif
#endif
#endif
INTEGER          :: J, JM, ID
INTEGER,DIMENSION(3) :: ITIM
REAL             :: ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT
!-------------------------------------------------------------------------------
#if defined(HPPA)
CALL DATE(YDAT8)
CALL TIME(YTIM8)
#else
#if defined(LINUX) || defined (O2000) 
CALL DATE_AND_TIME(YDAT8,YTIM10)
#else
#if defined(VPP)
CALL ITIME(ITIM)
YTIM8='        '
WRITE(YTIM8,'(I2,I2,I2)')ITIM
CALL DATE_AND_TIME(YDAT8,YTIM8)
YTEM8='        '
#endif
#endif
#endif

!!!!!!!!!!! Date
YTEM9='        '
#if defined(HPPA)
YTEM9(1:2)=YDAT8(1:2)
#else
#if defined(LINUX) || defined (O2000) 
YTEM9(1:2)=YDAT8(7:8)
#else
#if defined(VPP)
YTEM8(1:2)=YDAT8(7:8)
YTEM8(4:5)=YDAT8(4:5)
#endif
#endif
#endif
YTEM9(3:3)='/'
#if defined(HPPA)
YTEM9(4:6)=YDAT8(4:6)
#else
#if defined(LINUX) || defined (O2000) 
YTEM9(4:5)=YDAT8(5:6)
#else
#if defined(VPP)
YTEM8(3:3)='/'
YTEM8(6:6)='/'
#endif
#endif
#endif

#if defined(HPPA)
YTEM9(7:7)='/'
#else
#if defined(LINUX) || defined (O2000) 
YTEM9(6:6)='/'
#else
#if defined(VPP)
YTEM8(7:8)=YDAT8(1:2)
#endif
#endif
#endif
#if defined(HPPA)
YTEM9(8:9)=YDAT8(8:9)
#else
#if defined(LINUX) 
YTEM9(7:8)=YDAT8(3:4)
YTEM9(9:9)='/'
#if defined (O2000) 
YTEM9(7:8)=YDAT8(1:2)
YTEM9(9:9)='/'
#endif
#endif
#endif
#if defined(VPP)
YDAT8=YTEM8
#else
#if defined(HPPA)
YDAT8=YTEM9(1:9)
#else
YDAT8=YTEM9(1:8)
#endif
#endif

!!!!!!!!!!! Time
YTEM8='        '
#if defined(HPPA)
YTEM8(1:2)=YTIM8(1:2)
#else
#if defined(LINUX) || defined (O2000) 
YTEM8(1:2)=YTIM10(1:2)
#else
#if defined(VPP)
YTEM8(4:5)=YTIM8(3:4)
#endif
#endif
#endif
YTEM8(3:3)='H'

#if defined(HPPA)
YTEM8(4:5)=YTIM8(4:5)
#else
#if defined(LINUX) || defined (O2000) 
YTEM8(4:5)=YTIM10(3:4)
#else
#if defined(VPP)
YTEM8(7:8)=YTIM8(5:6)
#endif
#endif
#endif
YTEM8(6:6)='M'

#if defined(HPPA)
YTEM8(7:8)=YTIM8(7:8)
#else
#if defined(LINUX) || defined (O2000) 
YTEM8(7:8)=YTIM10(5:6)
#endif
#endif

YTIM8=YTEM8
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
CALL SET(0.,1.,0.,1.,0.,1.,0.,1.,1)
#if defined(HPPA)
CALL PLCHHQ(0.80,0.99,YDAT8,.008,0.,-1.)
#else
#if defined(LINUX) || defined (O2000) 
CALL PLCHHQ(0.80,0.99,YDAT8(1:LEN_TRIM(YDAT8)),.008,0.,-1.)
#else
#if defined(VPP)
CALL PLCHHQ(0.78,0.99,YDAT8,.008,0.,-1.)
#endif
#endif
#endif
#if defined(HPPA)
CALL PLCHHQ(0.99,0.99,YTIM8,.008,0.,+1.)
#else
#if defined(LINUX) || defined (O2000) 
CALL PLCHHQ(0.99,0.99,YTIM8(1:LEN_TRIM(YTIM8)),.008,0.,+1.)
#else
#if defined(VPP)
CALL PLCHHQ(0.90,0.99,YTIM8,.008,0.,-1.)
#endif
#endif
#endif
!
! Modifs for diachro
!
DO J=1,NBFILES
  IF(NUMFILES(J) == NUMFILECUR)THEN
    JM=J
    EXIT
  ENDIF
ENDDO
#if defined(HPPA)
CALL PLCHHQ(0.80,.97,CFILEDIAS(JM),.008,0.,-1.)
#else
#if defined(LINUX) || defined (O2000) 
CALL PLCHHQ(0.80,.97,CFILEDIAS(JM)(1:LEN_TRIM(CFILEDIAS(JM))),.008,0.,-1.)
#else
#if defined(VPP)
CALL PLCHHQ(0.78,.97,CFILEDIAS(JM),.008,0.,-1.)
#endif
#endif
#endif
IF(ALLOCATED(XVAR))THEN
IF(SIZE(XVAR,6) > 1 )THEN
  CALL PLCHHQ(0.99,.95,CGROUP(1:LEN_TRIM(CGROUP)),.008,0.,+1.)
ENDIF
ENDIF
IF(CTYPE == 'MASK')THEN
  CALL PLCHHQ(0.99,.95,CGROUP(1:LEN_TRIM(CGROUP)),.008,0.,+1.)
ENDIF
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
RETURN
END SUBROUTINE DATFILE_FORDIACHRO
