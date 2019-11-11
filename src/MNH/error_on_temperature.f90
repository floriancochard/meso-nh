!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
MODULE MODI_ERROR_ON_TEMPERATURE
!###############################
INTERFACE
      SUBROUTINE ERROR_ON_TEMPERATURE(PT_LS,PP_LS,PPABS,PPS_LS,PPS)
!
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PT_LS   ! temperature at begining
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PP_LS   ! pressure at begining
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PPABS   ! pressure at end
REAL,DIMENSION(:,:),   INTENT(IN)    :: PPS_LS  ! large-scale surface pressure
REAL,DIMENSION(:,:),   INTENT(IN)    :: PPS     ! Meso-NH surface pressure
!
END SUBROUTINE ERROR_ON_TEMPERATURE
END INTERFACE
END MODULE MODI_ERROR_ON_TEMPERATURE
!     ######spl
      SUBROUTINE ERROR_ON_TEMPERATURE(PT_LS,PP_LS,PPABS,PPS_LS,PPS)
!     ############################################################
!
!!****  *ERROR_ON_TEMPERATURE* - computes RMS on temperature between begining
!!                               and end of prep_real_case
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/07/96
!!                  25/10/96 (V. Masson) add deallocations
!!                  28/10/96 (V. Masson) bug in sums in bias and rms computations
!!                                       and does not use external points
!!                  11/07/97 (V. Masson) absolute pressure
!!                  26/08/97 (V. Masson) call to new linear vertical
!!                                       interpolation routine
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_FIELD_n, ONLY: XTHT
USE MODD_LUNIT,   ONLY: TLUOUT0
USE MODD_REF_n
USE MODD_VER_INTERP_LIN
!
USE MODI_COEF_VER_INTERP_LIN
USE MODI_VER_INTERP_LIN
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PT_LS   ! temperature at begining
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PP_LS   ! pressure at begining
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PPABS   ! pressure at end
REAL,DIMENSION(:,:),   INTENT(IN)    :: PPS_LS  ! large-scale surface pressure
REAL,DIMENSION(:,:),   INTENT(IN)    :: PPS     ! Meso-NH surface pressure
!
!*       0.2   Declaration of local variables
!              ------------------------------
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZT_LS_P  ! large scale T at pressure levels
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZT_P     ! meso-NH T at pressure levels
REAL,DIMENSION(:),     ALLOCATABLE:: ZBIAS    ! bias T between begining and end (at pressure levels)
REAL,DIMENSION(:),     ALLOCATABLE:: ZSIG2    ! sig2 on T between begining and end (at pressure levels)
REAL,DIMENSION(:),     ALLOCATABLE:: ZTRMS    ! rms on T between begining and end (at pressure levels)
REAL,DIMENSION(:),     ALLOCATABLE:: ZPLEVELS ! pressure levels where RMS are computed
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZP1,ZT1  ! work arrays
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZP2,ZT2  ! work arrays
LOGICAL,DIMENSION(:,:,:), ALLOCATABLE :: GMASK! .T. where pressure level is significant
!
INTEGER                           :: ILUOUT0
INTEGER                           :: IIU,IJU,ILU,IKU
INTEGER                           :: JP
!
!-------------------------------------------------------------------------------
!
!*     1.     Initializations
!             ---------------
!
IIU=SIZE(PT_LS,1)
IJU=SIZE(PT_LS,2)
ILU=SIZE(PT_LS,3)
IKU=SIZE(XTHT,3)
!
ALLOCATE(ZPLEVELS(20))
ZPLEVELS(:) = (/ (5000.*JP,JP=20,1,-1) /)
!
ALLOCATE(ZT_LS_P(IIU,IJU,20))
ALLOCATE(ZT_P(IIU,IJU,20))
ALLOCATE(ZTRMS(20))
ALLOCATE(ZBIAS(20))
ALLOCATE(ZSIG2(20))
!
!-------------------------------------------------------------------------------
!
!*     2.     Interpolations on pressure levels for begining field
!             ----------------------------------------------------
!
ALLOCATE(ZP1(IIU,IJU,ILU))
ALLOCATE(ZT1(IIU,IJU,ILU))
ALLOCATE(ZP2(IIU,IJU,20))
ALLOCATE(ZT2(IIU,IJU,20))
!
ZT1(:,:,:)=PT_LS(:,:,ILU:1:-1)
ZP1(:,:,:)=PP_LS(:,:,ILU:1:-1)
ZP2(:,:,:)=SPREAD(SPREAD(ZPLEVELS(20:1:-1),1,IIU),2,IJU)
CALL COEF_VER_INTERP_LIN(ZP1(:,:,:),ZP2(:,:,:))
ZT2(:,:,:)=VER_INTERP_LIN(ZT1(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
ZT_LS_P(:,:,:)=ZT2(:,:,20:1:-1)
!
DEALLOCATE(ZP1)
DEALLOCATE(ZT1)
DEALLOCATE(ZP2)
DEALLOCATE(ZT2)
!
!-------------------------------------------------------------------------------
!
!*     3.     Interpolations on pressure levels for ending field
!             --------------------------------------------------
!
ALLOCATE(ZP1(IIU,IJU,IKU))
ALLOCATE(ZT1(IIU,IJU,IKU))
ALLOCATE(ZP2(IIU,IJU,20))
ALLOCATE(ZT2(IIU,IJU,20))
!
ZT1(:,:,:)=XTHT(:,:,IKU:1:-1)*(PPABS(:,:,IKU:1:-1)/XP00)**(XRD/XCPD)
ZP1(:,:,:)=PPABS(:,:,IKU:1:-1)
ZP2(:,:,:)=SPREAD(SPREAD(ZPLEVELS(20:1:-1),1,IIU),2,IJU)
CALL COEF_VER_INTERP_LIN(ZP1(:,:,:),ZP2(:,:,:))
ZT2(:,:,:)=VER_INTERP_LIN(ZT1(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
ZT_P(:,:,:)=ZT2(:,:,20:1:-1)
!
DEALLOCATE(ZP1)
DEALLOCATE(ZT1)
DEALLOCATE(ZP2)
DEALLOCATE(ZT2)
!-------------------------------------------------------------------------------
!
!*     3.     Temperature RMS computation
!             ---------------------------
!
ALLOCATE(GMASK(IIU,IJU,20))
GMASK=   SPREAD(SPREAD(ZPLEVELS(1:20),1,IIU),2,IJU) < SPREAD(PPS_LS(:,:),3,20) &
   .AND. SPREAD(SPREAD(ZPLEVELS(1:20),1,IIU),2,IJU) < SPREAD(PPS(:,:),3,20)
GMASK( 1 , : ,:)=.FALSE.
GMASK(IIU, : ,:)=.FALSE.
GMASK( : , 1 ,:)=.FALSE.
GMASK( : ,IJU,:)=.FALSE.
!
DO JP =1, 20
  IF ( COUNT(GMASK(:,:,JP)) > 1 ) THEN
    ZBIAS(JP)= SUM( (ZT_LS_P(:,:,JP)-ZT_P(:,:,JP)          )   ,MASK=GMASK(:,:,JP)) / COUNT (GMASK(:,:,JP))
    ZSIG2(JP)= SUM( (ZT_LS_P(:,:,JP)-ZT_P(:,:,JP)-ZBIAS(JP))**2,MASK=GMASK(:,:,JP)) / COUNT (GMASK(:,:,JP))
    ZTRMS(JP)= SQRT(ZBIAS(JP)*ZBIAS(JP)+ZSIG2(JP))
  ELSE
    ZTRMS(JP)=0.
  END IF
END DO
!
!-------------------------------------------------------------------------------
!
!*     4.     Prints
!             ------
!
ILUOUT0 = TLUOUT0%NLU
!
WRITE(ILUOUT0,*) ''
WRITE(ILUOUT0,*) 'Temperature RMS between begin and end of PREP_REAL_CASE :'
WRITE(ILUOUT0,*) ''
DO JP=20,1,-1
  WRITE(ILUOUT0,'(6Hlevel ,F5.0,7H hPa : ,F5.3,2H K)') ZPLEVELS(JP)/100.,ZTRMS(JP)
END DO
WRITE(ILUOUT0,*) ''
!
!-------------------------------------------------------------------------------
!
!*     5.     Deallocations
!             -------------
!
DEALLOCATE(ZPLEVELS)
DEALLOCATE(ZT_LS_P)
DEALLOCATE(ZT_P)
DEALLOCATE(ZTRMS)
DEALLOCATE(ZBIAS)
DEALLOCATE(ZSIG2)
DEALLOCATE(GMASK)
!-------------------------------------------------------------------------------
!
WRITE(ILUOUT0,*) 'Routine ERROR_ON_TEMPERATURE completed'
!
END SUBROUTINE ERROR_ON_TEMPERATURE
