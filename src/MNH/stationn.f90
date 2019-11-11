!MNH_LIC Copyright 2002-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##########################
MODULE MODI_STATION_n
!      ##########################
!
INTERFACE
!
      SUBROUTINE STATION_n(PTSTEP,                               &
                           TPDTEXP, TPDTMOD, TPDTSEG, TPDTCUR,   &
                           PXHAT, PYHAT, PZ,                     &
                           PU, PV, PW, PTH, PR, PSV, PTKE,       &
                           PTS,PP ) 
!
USE MODD_TYPE_DATE
!
REAL,                     INTENT(IN)     :: PTSTEP ! time step
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTEXP! experiment date and time
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTMOD! model start date and time
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTSEG! segment date and time
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTCUR! current date and time
REAL, DIMENSION(:),       INTENT(IN)     :: PXHAT  ! x coordinate
REAL, DIMENSION(:),       INTENT(IN)     :: PYHAT  ! y coordinate
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PZ     ! z array
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PU     ! horizontal wind X component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PV     ! horizontal wind Y component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PW     ! vertical wind
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTH    ! potential temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PR     ! water mixing ratios
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PSV    ! Scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTKE   ! turbulent kinetic energy
REAL, DIMENSION(:,:),     INTENT(IN)     :: PTS    ! surface temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PP     ! pressure
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE STATION_n
!
END INTERFACE
!
END MODULE MODI_STATION_n
!
!     ########################################################
      SUBROUTINE STATION_n(PTSTEP,                           &
                       TPDTEXP, TPDTMOD, TPDTSEG, TPDTCUR,   &
                       PXHAT, PYHAT, PZ,                     &
                       PU, PV, PW, PTH, PR, PSV, PTKE,   &
                       PTS, PP )
!     ########################################################
!
!
!!****  *STATION_n* - (advects and) stores 
!!                                stations/s in the model
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
!!
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
!!      Pierre TULET / Valery Masson             * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original 15/02/2002
!!     A. Lemonsu 19/11/2002
!!     P.Aumond 01/07/2011 : Add model levels
!!     C.Lac       04/2013 : Correction on the vertical levels
!!     C.Lac       04/2013 : Add I/J positioning                   
!!     P.Wautelet 28/03/2018 : Replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_TYPE_DATE
USE MODD_STATION_n
USE MODD_SUB_STATION_n
USE MODD_DIAG_IN_RUN
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_GRID
USE MODD_TIME
USE MODD_CONF
!
USE MODE_DATETIME
USE MODE_ll
!
USE MODI_WATER_SUM
!
!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
REAL,                     INTENT(IN)     :: PTSTEP ! time step
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTEXP! experiment date and time
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTMOD! model start date and time
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTSEG! segment date and time
TYPE(DATE_TIME),          INTENT(IN)     :: TPDTCUR! current date and time
REAL, DIMENSION(:),       INTENT(IN)     :: PXHAT  ! x coordinate
REAL, DIMENSION(:),       INTENT(IN)     :: PYHAT  ! y coordinate
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PZ     ! z array
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PU     ! horizontal wind X component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PV     ! horizontal wind Y component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PW     ! vertical wind
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTH    ! potential temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PR     ! water mixing ratios
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PSV    ! Scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTKE   ! turbulent kinetic energy
REAL, DIMENSION(:,:),     INTENT(IN)     :: PTS    ! surface temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PP     ! pressure
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER :: IIB        ! current processor domain sizes
INTEGER :: IJB        ! 
INTEGER :: IIE        !    
INTEGER :: IJE        !   
INTEGER :: IIU        ! 
INTEGER :: IJU        ! 
REAL    :: ZTIMEEXP   ! 
!
REAL, DIMENSION(SIZE(PXHAT))        :: ZXHATM ! mass point coordinates
REAL, DIMENSION(SIZE(PYHAT))        :: ZYHATM ! mass point coordinates
!
REAL, DIMENSION(SIZE(PSV,1),SIZE(PSV,2),SIZE(PSV,3),SIZE(PSV,4))  :: ZWORK   ! 
!
LOGICAL       :: GSTORE                       ! storage occurs at this time step
!
!
INTEGER :: IN       ! time index
INTEGER :: JSV      ! loop counter
!
REAL    :: ZU_STAT     ! horizontal wind speed at station location (along x)
REAL    :: ZV_STAT     ! horizontal wind speed at station location (along y)
REAL    :: ZGAM        ! rotation between meso-nh base and spherical lat-lon base.
!
INTEGER :: IINFO_ll   ! return code
INTEGER :: IRESP      ! return code
INTEGER :: I          ! loop for stations
INTEGER :: J          ! loop for levels

!
!----------------------------------------------------------------------------
!
!*      2.   PRELIMINARIES
!            -------------
!
!*      2.1  Indices
!            -------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
!
!*      2.2  Interpolations of model variables to mass points
!            ------------------------------------------------
!
IIU=SIZE(PXHAT)
IJU=SIZE(PYHAT)
!
ZXHATM(1:IIU-1)=0.5*PXHAT(1:IIU-1)+0.5*PXHAT(2:IIU  )
ZXHATM(  IIU  )=1.5*PXHAT(  IIU  )-0.5*PXHAT(  IIU-1)
!
ZYHATM(1:IJU-1)=0.5*PYHAT(1:IJU-1)+0.5*PYHAT(2:IJU  )
ZYHATM(  IJU  )=1.5*PYHAT(  IJU  )-0.5*PYHAT(  IJU-1)
!
!----------------------------------------------------------------------------
!
!*      3.4  instant of storage
!            ------------------
!
IF ( TSTATION%T_CUR == XUNDEF ) TSTATION%T_CUR = TSTATION%STEP - PTSTEP
!
TSTATION%T_CUR = TSTATION%T_CUR + PTSTEP
!
CALL DATETIME_DISTANCE(TDTEXP,TDTSEG,ZTIMEEXP)
IF ( TSTATION%T_CUR >= TSTATION%STEP - 1.E-10 ) THEN
     GSTORE = .TRUE.
     TSTATION%T_CUR = TSTATION%T_CUR - TSTATION%STEP
     TSTATION%N_CUR = TSTATION%N_CUR + 1
     IN = TSTATION%N_CUR
ELSE
     GSTORE = .FALSE.
END IF
!
IF (GSTORE) THEN
  !
     TSTATION%TIME(IN)      = (IN-1) * TSTATION%STEP + ZTIMEEXP
     TSTATION%DATIME( 1,IN) = TPDTEXP%TDATE%YEAR
     TSTATION%DATIME( 2,IN) = TPDTEXP%TDATE%MONTH
     TSTATION%DATIME( 3,IN) = TPDTEXP%TDATE%DAY
     TSTATION%DATIME( 4,IN) = TPDTEXP%TIME
     TSTATION%DATIME( 5,IN) = TPDTSEG%TDATE%YEAR
     TSTATION%DATIME( 6,IN) = TPDTSEG%TDATE%MONTH
     TSTATION%DATIME( 7,IN) = TPDTSEG%TDATE%DAY
     TSTATION%DATIME( 8,IN) = TPDTSEG%TIME
     TSTATION%DATIME( 9,IN) = TPDTMOD%TDATE%YEAR
     TSTATION%DATIME(10,IN) = TPDTMOD%TDATE%MONTH
     TSTATION%DATIME(11,IN) = TPDTMOD%TDATE%DAY
     TSTATION%DATIME(12,IN) = TPDTMOD%TIME
     TSTATION%DATIME(13,IN) = TPDTCUR%TDATE%YEAR
     TSTATION%DATIME(14,IN) = TPDTCUR%TDATE%MONTH
     TSTATION%DATIME(15,IN) = TPDTCUR%TDATE%DAY
     TSTATION%DATIME(16,IN) = TPDTCUR%TIME
END IF
!
!
!----------------------------------------------------------------------------
!
!*      4.   STATION POSITION
!            --------------
!
!*      4.0  initialization of processor test
!            --------------------------------
IF (GSTATFIRSTCALL) THEN
 GSTATFIRSTCALL=.FALSE.
!
 IF (.NOT.(ASSOCIATED(ZTHIS_PROCS))) ALLOCATE(ZTHIS_PROCS(NUMBSTAT))
!
 IF (.NOT.(ASSOCIATED(II)))     ALLOCATE(II(NUMBSTAT))
 IF (.NOT.(ASSOCIATED(IJ)))     ALLOCATE(IJ(NUMBSTAT))
 IF (.NOT.(ASSOCIATED(IV)))     ALLOCATE(IV(NUMBSTAT))
 IF (.NOT.(ASSOCIATED(IU)))     ALLOCATE(IU(NUMBSTAT))
 IF (.NOT.(ASSOCIATED(ZXCOEF))) ALLOCATE(ZXCOEF(NUMBSTAT))
 IF (.NOT.(ASSOCIATED(ZUCOEF))) ALLOCATE(ZUCOEF(NUMBSTAT))
 IF (.NOT.(ASSOCIATED(ZYCOEF))) ALLOCATE(ZYCOEF(NUMBSTAT))
 IF (.NOT.(ASSOCIATED(ZVCOEF))) ALLOCATE(ZVCOEF(NUMBSTAT))

 ZXCOEF(:)  =XUNDEF
 ZUCOEF(:)  =XUNDEF
 ZYCOEF(:)  =XUNDEF
 ZVCOEF(:)  =XUNDEF

!
 DO I=1,NUMBSTAT
!
  ZTHIS_PROCS(I)=0.
!
!*      4.1  X position
!            ----------
!
  IU(I)=COUNT( PXHAT (:)<=TSTATION%X(I) )
  II(I)=COUNT( ZXHATM(:)<=TSTATION%X(I) )
!
  IF (II(I)<=IIB-1   .AND. LWEST_ll() .AND. .NOT. L1D) TSTATION%ERROR(I)=.TRUE.
  IF (II(I)>=IIE     .AND. LEAST_ll() .AND. .NOT. L1D) TSTATION%ERROR(I)=.TRUE.
!
!
!*      4.2  Y position
!            ----------
!
  IV(I)=COUNT( PYHAT (:)<=TSTATION%Y(I) )
  IJ(I)=COUNT( ZYHATM(:)<=TSTATION%Y(I) )
!
  IF (IJ(I)<=IJB-1   .AND. LSOUTH_ll() .AND. .NOT. L1D) TSTATION%ERROR(I)=.TRUE.
  IF (IJ(I)>=IJE     .AND. LNORTH_ll() .AND. .NOT. L1D) TSTATION%ERROR(I)=.TRUE.
!
!
!*      4.3  Position of station according to processors
!            -------------------------------------------
!
  IF (IU(I)>=IIB .AND. IU(I)<=IIE .AND. IV(I)>=IJB .AND. IV(I)<=IJE) ZTHIS_PROCS(I)=1.
  IF (L1D) ZTHIS_PROCS(I)=1.
!
!
!*      4.4  Computations only on correct processor
!            --------------------------------------
  IF ( LSTATLAT ) THEN
    ZXCOEF(I) = 0.
    ZYCOEF(I) = 0.
    ZUCOEF(I) = 0.         
    ZVCOEF(I) = 0.
    IF (ZTHIS_PROCS(I) >0. .AND. .NOT. L1D) THEN
!----------------------------------------------------------------------------
!
!*      6.1  Interpolation coefficient for X
!            -------------------------------
!
       ZXCOEF(I) = (TSTATION%X(I) - ZXHATM(II(I))) / (ZXHATM(II(I)+1) - ZXHATM(II(I)))
!
!
!
!*      6.2  Interpolation coefficient for y
!            -------------------------------
!
       ZYCOEF(I) = (TSTATION%Y(I) - ZYHATM(IJ(I))) / (ZYHATM(IJ(I)+1) - ZYHATM(IJ(I)))
!
!-------------------------------------------------------------------
!
!*      7.   INITIALIZATIONS FOR INTERPOLATIONS OF U AND V
!            ---------------------------------------------
!
!*      7.1  Interpolation coefficient for X (for U)
!            -------------------------------
!
       ZUCOEF(I) = (TSTATION%X(I) - PXHAT(IU(I))) / (PXHAT(IU(I)+1) - PXHAT(IU(I)))
!
!
!*      7.2  Interpolation coefficient for y (for V)
!            -------------------------------
!
       ZVCOEF(I) = (TSTATION%Y(I) - PYHAT(IV(I))) / (PYHAT(IV(I)+1) - PYHAT(IV(I)))
!
!

    END IF
  END IF
 ENDDO
END IF
!----------------------------------------------------------------------------
!
!*      8.   DATA RECORDING
!            --------------
!
IF (GSTORE) THEN

  IF (TSTATION%TIME(IN) /= XUNDEF) THEN     

 DO I=1,NUMBSTAT                  
     !
     IF ((ZTHIS_PROCS(I)==1.).AND.(.NOT. TSTATION%ERROR(I))) THEN
       IF (TSTATION%K(I)/= XUNDEF) THEN
         J = TSTATION%K(I)
       ELSE  ! suppose TSTATION%Z(I) /= XUNDEF
        J=1
        DO WHILE ((STATION_INTERP_2D(PZ(:,:,J))-STATION_INTERP_2D(PZ(:,:,2))) &
        < TSTATION%Z(I))
         J = J + 1
        END DO
        IF (((STATION_INTERP_2D(PZ(:,:,J))-STATION_INTERP_2D(PZ(:,:,2)))-TSTATION%Z(I))>&
       (TSTATION%Z(I)-(STATION_INTERP_2D(PZ(:,:,J-1))-STATION_INTERP_2D(PZ(:,:,2))))) THEN
         J=J-1
        ENDIF
       END IF
      !
      ZGAM                  = (XRPK * (TSTATION%LON(I) - XLON0) - XBETA)*(XPI/180.)
      IF ( LSTATLAT ) THEN
       ZU_STAT               = STATION_INTERP_2D_U(PU(:,:,J))
       ZV_STAT               = STATION_INTERP_2D_V(PV(:,:,J))
      ELSE
       ZU_STAT               = PU(TSTATION%I(I),TSTATION%J(I),J)
       ZV_STAT               = PV(TSTATION%I(I),TSTATION%J(I),J)
      END IF
      !
      TSTATION%ZON (IN,I)   =   ZU_STAT     * COS(ZGAM) + ZV_STAT     * SIN(ZGAM)
      TSTATION%MER (IN,I)   = - ZU_STAT     * SIN(ZGAM) + ZV_STAT     * COS(ZGAM)
      IF ( LSTATLAT ) THEN
        TSTATION%W   (IN,I)   = STATION_INTERP_2D(PW(:,:,J))
        TSTATION%TH  (IN,I)   = STATION_INTERP_2D(PTH(:,:,J))
        TSTATION%P   (IN,I)   = STATION_INTERP_2D(PP(:,:,J))
      !
        DO JSV=1,SIZE(PR,4)
         TSTATION%R   (IN,I,JSV) = STATION_INTERP_2D(PR(:,:,J,JSV))
        END DO
      !
        DO JSV=1,SIZE(PSV,4)
         TSTATION%SV  (IN,I,JSV) = STATION_INTERP_2D(PSV(:,:,J,JSV))
        END DO
      !
        IF (SIZE(PTKE)>0) TSTATION%TKE  (IN,I) = STATION_INTERP_2D(PTKE(:,:,J))
        IF (SIZE(PTS) >0) TSTATION%TSRAD(IN,I) = STATION_INTERP_2D(PTS)
        TSTATION%ZS(I)      = STATION_INTERP_2D(PZ(:,:,1+JPVEXT))
      !
        IF (LDIAG_IN_RUN) THEN
          TSTATION%ZON10M(IN,I) = STATION_INTERP_2D(XCURRENT_ZON10M)
          TSTATION%MER10M(IN,I) = STATION_INTERP_2D(XCURRENT_MER10M)
          TSTATION%T2M   (IN,I) = STATION_INTERP_2D(XCURRENT_T2M   )
          TSTATION%Q2M   (IN,I) = STATION_INTERP_2D(XCURRENT_Q2M   )
          TSTATION%HU2M  (IN,I) = STATION_INTERP_2D(XCURRENT_HU2M  )
          TSTATION%RN    (IN,I) = STATION_INTERP_2D(XCURRENT_RN    ) 
          TSTATION%H     (IN,I) = STATION_INTERP_2D(XCURRENT_H     ) 
          TSTATION%LE    (IN,I) = STATION_INTERP_2D(XCURRENT_LE    ) 
          TSTATION%LEI   (IN,I) = STATION_INTERP_2D(XCURRENT_LEI   ) 
          TSTATION%GFLUX (IN,I) = STATION_INTERP_2D(XCURRENT_GFLUX ) 
          TSTATION%SWD   (IN,I) = STATION_INTERP_2D(XCURRENT_SWD   ) 
          TSTATION%SWU   (IN,I) = STATION_INTERP_2D(XCURRENT_SWU   ) 
          TSTATION%LWD   (IN,I) = STATION_INTERP_2D(XCURRENT_LWD   ) 
          TSTATION%LWU   (IN,I) = STATION_INTERP_2D(XCURRENT_LWU   ) 
          TSTATION%SWDIR (IN,I) = STATION_INTERP_2D(XCURRENT_SWDIR ) 
          TSTATION%SWDIFF(IN,I) = STATION_INTERP_2D(XCURRENT_SWDIFF)
          TSTATION%DSTAOD(IN,I) = STATION_INTERP_2D(XCURRENT_DSTAOD)
          TSTATION%SFCO2 (IN,I) = STATION_INTERP_2D(XCURRENT_SFCO2 ) 
        ENDIF
       ELSE
        TSTATION%W   (IN,I)   = PW(TSTATION%I(I),TSTATION%J(I),J)
        TSTATION%TH  (IN,I)   = PTH(TSTATION%I(I),TSTATION%J(I),J)
        TSTATION%P   (IN,I)   = PP(TSTATION%I(I),TSTATION%J(I),J)
      !
        DO JSV=1,SIZE(PR,4)
         TSTATION%R   (IN,I,JSV) = PR(TSTATION%I(I),TSTATION%J(I),J,JSV)
        END DO
      !
        DO JSV=1,SIZE(PSV,4)
         TSTATION%SV  (IN,I,JSV) = PSV(TSTATION%I(I),TSTATION%J(I),J,JSV)
        END DO
      !
        IF (SIZE(PTKE)>0) TSTATION%TKE  (IN,I) = PTKE(TSTATION%I(I),TSTATION%J(I),J)
        IF (SIZE(PTS) >0) TSTATION%TSRAD(IN,I) = PTS(TSTATION%I(I),TSTATION%J(I))
        TSTATION%ZS(I)      = PZ(TSTATION%I(I),TSTATION%J(I),1+JPVEXT)
      !
        IF (LDIAG_IN_RUN) THEN
          TSTATION%ZON10M(IN,I) = XCURRENT_ZON10M(TSTATION%I(I),TSTATION%J(I))
          TSTATION%MER10M(IN,I) = XCURRENT_MER10M(TSTATION%I(I),TSTATION%J(I))
          TSTATION%T2M   (IN,I) = XCURRENT_T2M(TSTATION%I(I),TSTATION%J(I))
          TSTATION%Q2M   (IN,I) = XCURRENT_Q2M(TSTATION%I(I),TSTATION%J(I))
          TSTATION%HU2M  (IN,I) = XCURRENT_HU2M(TSTATION%I(I),TSTATION%J(I))
          TSTATION%RN    (IN,I) = XCURRENT_RN(TSTATION%I(I),TSTATION%J(I))
          TSTATION%H     (IN,I) = XCURRENT_H(TSTATION%I(I),TSTATION%J(I))
          TSTATION%LE    (IN,I) = XCURRENT_LE(TSTATION%I(I),TSTATION%J(I))
          TSTATION%LEI   (IN,I) = XCURRENT_LEI(TSTATION%I(I),TSTATION%J(I))
          TSTATION%GFLUX (IN,I) = XCURRENT_GFLUX(TSTATION%I(I),TSTATION%J(I))
          TSTATION%SWD   (IN,I) = XCURRENT_SWD(TSTATION%I(I),TSTATION%J(I))
          TSTATION%SWU   (IN,I) = XCURRENT_SWU(TSTATION%I(I),TSTATION%J(I))
          TSTATION%LWD   (IN,I) = XCURRENT_LWD(TSTATION%I(I),TSTATION%J(I))
          TSTATION%LWU   (IN,I) = XCURRENT_LWU(TSTATION%I(I),TSTATION%J(I))
          TSTATION%SWDIR (IN,I) = XCURRENT_SWDIR(TSTATION%I(I),TSTATION%J(I))
          TSTATION%SWDIFF(IN,I) = XCURRENT_SWDIFF(TSTATION%I(I),TSTATION%J(I))         
          TSTATION%DSTAOD(IN,I) = XCURRENT_DSTAOD(TSTATION%I(I),TSTATION%J(I))
          TSTATION%SFCO2 (IN,I) = XCURRENT_SFCO2(TSTATION%I(I),TSTATION%J(I))
        ENDIF
       ENDIF
      !
    END IF
!
!----------------------------------------------------------------------------
!
!*     11.   EXCHANGE OF INFORMATION BETWEEN PROCESSORS
!            ------------------------------------------
!
!*     11.2  data stored
!            -----------
!
  CALL DISTRIBUTE_STATION(TSTATION%X   (I))
  CALL DISTRIBUTE_STATION(TSTATION%Y   (I))
  CALL DISTRIBUTE_STATION(TSTATION%Z   (I))
  CALL DISTRIBUTE_STATION(TSTATION%LON (I))
  CALL DISTRIBUTE_STATION(TSTATION%LAT (I))    
  CALL DISTRIBUTE_STATION(TSTATION%ZON (IN,I))
  CALL DISTRIBUTE_STATION(TSTATION%MER (IN,I))
  CALL DISTRIBUTE_STATION(TSTATION%W   (IN,I))
  CALL DISTRIBUTE_STATION(TSTATION%P   (IN,I))
  CALL DISTRIBUTE_STATION(TSTATION%TH  (IN,I))
  DO JSV=1,SIZE(PR,4)
    CALL DISTRIBUTE_STATION(TSTATION%R   (IN,I,JSV))
  END DO
  DO JSV=1,SIZE(PSV,4)
    CALL DISTRIBUTE_STATION(TSTATION%SV  (IN,I,JSV))
  END DO
  IF (SIZE(PTKE)>0) CALL DISTRIBUTE_STATION(TSTATION%TKE  (IN,I))
  IF (SIZE(PTS) >0) CALL DISTRIBUTE_STATION(TSTATION%TSRAD(IN,I))
  CALL DISTRIBUTE_STATION(TSTATION%ZS  (I))
  IF (LDIAG_IN_RUN) THEN
    CALL DISTRIBUTE_STATION(TSTATION%T2M    (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%Q2M    (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%HU2M   (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%ZON10M (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%MER10M (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%RN     (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%H      (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%LE     (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%LEI    (IN,I))    
    CALL DISTRIBUTE_STATION(TSTATION%GFLUX  (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%SWD    (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%SWU    (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%LWD    (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%LWU    (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%SWDIR  (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%SWDIFF (IN,I))    
    CALL DISTRIBUTE_STATION(TSTATION%DSTAOD (IN,I))
    CALL DISTRIBUTE_STATION(TSTATION%SFCO2  (IN,I))
  ENDIF
  !
 ENDDO
  !
 END IF
  !
END IF
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
CONTAINS
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!
FUNCTION STATION_INTERP_2D(PA) RESULT(PB)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PA
REAL                             :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
     JI=1
     JJ=1 
ELSEIF (L1D) THEN
     JI=2
     JJ=2
ELSE
     JI=II(I)
     JJ=IJ(I)
END IF
!
!
IF ((JI .GE. 1).AND. (JI .LE. SIZE(PA,1)) .AND. &
    (JJ .GE. 1).AND. (JJ .LE. SIZE(PA,2))) &
PB = (1.-ZYCOEF(I)) * (1.-ZXCOEF(I)) *  PA(JI,JJ)    + &
     (1.-ZYCOEF(I)) *    (ZXCOEF(I)) *  PA(JI+1,JJ)  + &
     (   ZYCOEF(I)) * (1.-ZXCOEF(I)) *  PA(JI,JJ+1)  + &
     (   ZYCOEF(I)) *    (ZXCOEF(I)) *  PA(JI+1,JJ+1)
!
END FUNCTION STATION_INTERP_2D
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! MODIFS
FUNCTION STATION_INTERP_2D_U(PA) RESULT(PB)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PA
REAL                             :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
     JI=1
     JJ=1
ELSEIF (L1D) THEN
     JI=2
     JJ=2
ELSE
     JI=II(I)
     JJ=IJ(I)
END IF
!
IF ((JI .GE. 1).AND. (JI .LE. SIZE(PA,1)) .AND. &
    (JJ .GE. 1).AND. (JJ .LE. SIZE(PA,2))) &
PB = (1.- ZYCOEF(I)) * (1.-ZUCOEF(I)) * PA(JI  ,JJ  ) &
   + (1.- ZYCOEF(I)) * (   ZUCOEF(I)) * PA(JI+1,JJ  ) &
   + (    ZYCOEF(I)) * (1.-ZUCOEF(I)) * PA(JI  ,JJ+1) &
   + (    ZYCOEF(I)) * (   ZUCOEF(I)) * PA(JI+1,JJ+1)
!
END FUNCTION STATION_INTERP_2D_U
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
! MODIFS
FUNCTION STATION_INTERP_2D_V(PA) RESULT(PB)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PA
REAL                             :: PB
!
INTEGER :: JI, JJ
!
IF (SIZE(PA,1)==2) THEN
  JI=1
  JJ=1
ELSEIF (L1D) THEN
     JI=2
     JJ=2  
ELSE
  JI=II(I)
  JJ=IJ(I)
END IF
!
IF ((JI .GT. 0).AND. (JI .LT. SIZE(PA,1)) .AND. &
    (JJ .GT. 0).AND. (JJ .LT. SIZE(PA,2))) &
PB = (1.- ZVCOEF(I)) * (1.-ZXCOEF(I)) * PA(JI  ,JJ  ) &
   + (1.- ZVCOEF(I)) * (   ZXCOEF(I)) * PA(JI+1,JJ  ) &
   + (    ZVCOEF(I)) * (1.-ZXCOEF(I)) * PA(JI  ,JJ+1) &
   + (    ZVCOEF(I)) * (   ZXCOEF(I)) * PA(JI+1,JJ+1)
!
END FUNCTION STATION_INTERP_2D_V
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE DISTRIBUTE_STATION(PAS)
!
REAL, INTENT(INOUT) :: PAS
!
PAS = PAS * ZTHIS_PROCS(I)
CALL REDUCESUM_ll(PAS,IINFO_ll)
!
END SUBROUTINE DISTRIBUTE_STATION
!----------------------------------------------------------------------------
!
END SUBROUTINE STATION_n
