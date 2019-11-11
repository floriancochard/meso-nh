!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ###########################
      MODULE MODI_HOR_INTERP_4PTS
!     ###########################
INTERFACE HOR_INTERP_4PTS
      SUBROUTINE HOR_INTERP_4PTS_2D(PX1,PY1,PFIELD1,PX2,PY2,PFIELD2)
      
!
REAL,   DIMENSION(:),INTENT(IN)       :: PX1       ! x of each grid mesh.
REAL,   DIMENSION(:),INTENT(IN)       :: PY1       ! y of each grid mesh.
REAL,   DIMENSION(:,:),INTENT(IN)     :: PFIELD1   ! field on grid mesh
!
REAL,   DIMENSION(:,:),INTENT(IN)     :: PX2       ! x of each new grid mesh.
REAL,   DIMENSION(:,:),INTENT(IN)     :: PY2       ! y of each new grid mesh.
REAL,   DIMENSION(:,:),INTENT(OUT)    :: PFIELD2   ! field on new grid mesh
!
END SUBROUTINE HOR_INTERP_4PTS_2D
!
      SUBROUTINE HOR_INTERP_4PTS_3D(PX1,PY1,PFIELD1,PX2,PY2,PFIELD2)
      
!
REAL,   DIMENSION(:),INTENT(IN)       :: PX1       ! x of each grid mesh.
REAL,   DIMENSION(:),INTENT(IN)       :: PY1       ! y of each grid mesh.
REAL,   DIMENSION(:,:,:),INTENT(IN)   :: PFIELD1   ! field on grid mesh
!
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PX2       ! x of each new grid mesh.
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PY2       ! y of each new grid mesh.
REAL,   DIMENSION(:,:,:),INTENT(OUT)  :: PFIELD2   ! field on new grid mesh
!
END SUBROUTINE HOR_INTERP_4PTS_3D
END INTERFACE
END MODULE MODI_HOR_INTERP_4PTS
!
!
!     ##############################
      MODULE MODI_HOR_INTERP_4PTS_3D
!     ##############################
INTERFACE HOR_INTERP_4PTS_3D
      SUBROUTINE HOR_INTERP_4PTS_3D(PX1,PY1,PFIELD1,PX2,PY2,PFIELD2)
      
!
REAL,   DIMENSION(:),INTENT(IN)       :: PX1       ! x of each grid mesh.
REAL,   DIMENSION(:),INTENT(IN)       :: PY1       ! y of each grid mesh.
REAL,   DIMENSION(:,:,:),INTENT(IN)   :: PFIELD1   ! field on grid mesh
!
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PX2       ! x of each new grid mesh.
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PY2       ! y of each new grid mesh.
REAL,   DIMENSION(:,:,:),INTENT(OUT)  :: PFIELD2   ! field on new grid mesh
!
END SUBROUTINE HOR_INTERP_4PTS_3D
END INTERFACE
END MODULE MODI_HOR_INTERP_4PTS_3D
!
!     ##############################################################
      SUBROUTINE HOR_INTERP_4PTS_3D(PX1,PY1,PFIELD1,PX2,PY2,PFIELD2)
!     ##############################################################
!
!!**** *HOR_INTERP_4PTS* interpolates horizontally a 3D field from a
!!                       REGULAR horizontal grid to any other grid
!!
!!    PURPOSE
!!    -------
!!
!!
!!    METHOD
!!    ------
!!   
!!    Bogus value of input field is XUNDEF
!!
!!    The routine uses only the points with physical values for interpolation:
!!       4pts available: interpolations linear in the 2 directions
!!       3pts available: plane interpolation
!!       2pts available: linear interpolation
!!       1pt  available: copy
!!
!!    Bogus value returned where field could not be interpolated is XUNDEF
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    19/03/95
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!USE MODE_FM
!USE MODD_LUNIT
USE MODD_PARAMETERS, ONLY: XUNDEF
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL,   DIMENSION(:),INTENT(IN)       :: PX1       ! x of each grid mesh.
REAL,   DIMENSION(:),INTENT(IN)       :: PY1       ! y of each grid mesh.
REAL,   DIMENSION(:,:,:),INTENT(IN)   :: PFIELD1   ! field on grid mesh
!
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PX2       ! x of each new grid mesh.
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PY2       ! y of each new grid mesh.
REAL,   DIMENSION(:,:,:),INTENT(OUT)  :: PFIELD2   ! field on new grid mesh
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                 :: ILUOUT0              ! logical unit
INTEGER                 :: IRESP                ! return codes
INTEGER                 :: JK
INTEGER                 :: IIU,IJU,IIOUT,IJOUT,II,IJ
INTEGER                 :: JI,JJ,JIOUT,JJOUT
REAL                    :: ZEPS
REAL :: ZXA,ZXB,ZXC,ZXD,ZYA,ZYB,ZYC,ZYD,ZA,ZB,ZC,ZD
REAL, DIMENSION(3) :: ZX,ZY,ZF
INTEGER :: JLOOP
REAL :: ZDET,ZALPHA,ZBETA,ZGAMMA
!
!-------------------------------------------------------------------------------
!
!*    1.     Initializations
!            ---------------
!
print *,'HOR_INTERP_4PTS: old grid',SIZE(PX1),SIZE(PY1), &
                         'new grid ',SIZE(PX2,1),SIZE(PY2,2)
!CALL FMLOOK_ll(CLUOUT0,CLUOUT0,ILUOUT0,IRESP)
IIOUT=SIZE(PX2,1)
IJOUT=SIZE(PY2,2)
IIU=SIZE(PX1)
IJU=SIZE(PY1)
ZEPS=1.E-10
!
!-------------------------------------------------------------------------------
!
!print * ,' avant boucle JK= k i j iold jold',SIZE(PFIELD1,3) ,IIOUT,IJOUT,SIZE(PX1),SIZE(PY1)
!print *, 'PX1(fin), PY1(fin)', PX1(SIZE(PX1)), PY1(SIZE(PY1))
!print *, 'PX2(IIOUT,IJOUT), PY2(IIOUT,IJOUT)', PX2(IIOUT,IJOUT), PY2(IIOUT,IJOUT)
DO JK=1,SIZE(PFIELD1,3)
  DO JIOUT=1,IIOUT
    DO JJOUT=1,IJOUT
      II=COUNT(PX1(:)<PX2(JIOUT,JJOUT))
      IJ=COUNT(PY1(:)<PY2(JIOUT,JJOUT))
      IF ( II<1 .OR. II>=IIU .OR. IJ<1 .OR. IJ>=IJU) THEN
        PFIELD2(JIOUT,JJOUT,:)=XUNDEF
        !print *,'pt nouvelle grille hors ancienne grille i j nbi nbj:',JIOUT,JJOUT,II,IJ
        !print *,'PX2(JIOUT,JJOUT),PY2(JIOUT,JJOUT)' ,PX2(JIOUT,JJOUT),PY2(JIOUT,JJOUT) 
        CYCLE
      END IF
!
      !print *,' valeur non indef i j nbi nbj:',JIOUT,JJOUT, II,IJ
      ZXA=PX1(II)
      ZXB=PX1(II)
      ZXC=PX1(II+1)
      ZXD=PX1(II+1)
!
      ZYA=PY1(IJ)
      ZYB=PY1(IJ+1)
      ZYC=PY1(IJ)
      ZYD=PY1(IJ+1)
!
      ZA=PFIELD1(II,IJ,JK)
      ZB=PFIELD1(II,IJ+1,JK)
      ZC=PFIELD1(II+1,IJ,JK)
      ZD=PFIELD1(II+1,IJ+1,JK)
!
      IF (ALL(ABS(PFIELD1(II:II+1,IJ:IJ+1,JK)-XUNDEF)<ZEPS) ) THEN
        !print * ,' 4 points a indef  :', PFIELD1(II:II+1,IJ:IJ+1,JK)
        PFIELD2(JIOUT,JJOUT,JK)=XUNDEF
        CYCLE
      ELSE IF (ALL(ABS(PFIELD1(II:II+1,IJ:IJ+1,JK)-XUNDEF)>=ZEPS) ) THEN
        ZALPHA=ZA+(ZB-ZA)*(PY2(JIOUT,JJOUT)-ZYA)/(ZYB-ZYA)
        ZBETA =ZC+(ZD-ZC)*(PY2(JIOUT,JJOUT)-ZYC)/(ZYD-ZYC)
        PFIELD2(JIOUT,JJOUT,JK)=ZALPHA+(ZBETA-ZALPHA)*(PX2(JIOUT,JJOUT)-ZXA)/(ZXC-ZXA)
      ELSE
        JLOOP=0
        DO JI=II,II+1
          DO JJ=IJ,IJ+1
            IF (ABS(PFIELD1(JI,JJ,JK)-XUNDEF)>ZEPS) THEN
              JLOOP=JLOOP+1
              ZX(JLOOP)=PX1(JI)
              ZY(JLOOP)=PY1(JJ)
              ZF(JLOOP)=PFIELD1(JI,JJ,JK)
            END IF
          END DO
        END DO
        IF (JLOOP==1) THEN
          PFIELD2(JIOUT,JJOUT,JK)=ZF(1)
        ELSE IF (JLOOP==2) THEN
          IF (ABS(ZX(1)-ZX(2))>ZEPS) THEN
            PFIELD2(JIOUT,JJOUT,JK)=ZF(1)+(ZF(2)-ZF(1))*(PX2(JIOUT,JJOUT)-ZX(1))/(ZX(2)-ZX(1))
          ELSE
            PFIELD2(JIOUT,JJOUT,JK)=ZF(1)+(ZF(2)-ZF(1))*(PY2(JIOUT,JJOUT)-ZY(1))/(ZY(2)-ZY(1))
          END IF
        ELSE IF (JLOOP==3) THEN
          ZDET=(ZX(1)-ZX(3))*(ZY(2)-ZY(3))-(ZX(2)-ZX(3))*(ZY(1)-ZY(3))
          ZALPHA=( (ZF(1)-ZF(3))*(ZY(2)-ZY(3))-(ZF(2)-ZF(3))*(ZY(1)-ZY(3)) )/ZDET
          ZBETA=-( (ZF(1)-ZF(3))*(ZX(2)-ZX(3))-(ZF(2)-ZF(3))*(ZX(1)-ZX(3)) )/ZDET
          ZGAMMA=ZF(1)-ZALPHA*ZX(1)-ZBETA*ZY(1)
          PFIELD2(JIOUT,JJOUT,JK)=ZALPHA*PX2(JIOUT,JJOUT) &
                                   +ZBETA *PY2(JIOUT,JJOUT) &
                                   +ZGAMMA
        END IF
      END IF
    END DO
  END DO
END DO
print *, 'fin routine HOR_INTERP_4PTS_3D'
!-------------------------------------------------------------------------------
!
!WRITE(ILUOUT0,*) ' Routine HOR_INTERP_4PTS completed'
!
END SUBROUTINE HOR_INTERP_4PTS_3D
!
!     ##############################################################
      SUBROUTINE HOR_INTERP_4PTS_2D(PX1,PY1,PFIELD1,PX2,PY2,PFIELD2)
!     ##############################################################
!
!!**** *HOR_INTERP_4PTS* interpolates horizontally a 2D field from a
!!                       REGULAR horizontal grid to any other grid
!!
!!    PURPOSE
!!    -------
!!
!!
!!    METHOD
!!    ------
!!   
!!    Bogus value of input field is XUNDEF
!!
!!    The routine uses only the points with physical values for interpolation:
!!       4pts available: interpolations linear in the 2 directions
!!       3pts available: plane interpolation
!!       2pts available: linear interpolation
!!       1pt  available: copy
!!
!!    Bogus value returned where field could not be interpolated is XUNDEF
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    04/07/96
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODI_HOR_INTERP_4PTS_3D
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL,   DIMENSION(:),INTENT(IN)       :: PX1       ! x of each grid mesh.
REAL,   DIMENSION(:),INTENT(IN)       :: PY1       ! y of each grid mesh.
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PFIELD1   ! field on grid mesh
!
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PX2       ! x of each new grid mesh.
REAL,   DIMENSION(:,:),  INTENT(IN)   :: PY2       ! y of each new grid mesh.
REAL,   DIMENSION(:,:),  INTENT(OUT)  :: PFIELD2   ! field on new grid mesh
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL, DIMENSION(SIZE(PFIELD1,1),SIZE(PFIELD1,2),1) :: ZFIELD1
REAL, DIMENSION(SIZE(PFIELD2,1),SIZE(PFIELD2,2),1) :: ZFIELD2
!
!-------------------------------------------------------------------------------
!
ZFIELD1(:,:,1)=PFIELD1(:,:)
CALL HOR_INTERP_4PTS_3D(PX1(:),PY1(:),ZFIELD1(:,:,:), &
                        PX2(:,:),PY2(:,:),ZFIELD2(:,:,:))
PFIELD2(:,:)=ZFIELD2(:,:,1)
!-------------------------------------------------------------------------------
!
END SUBROUTINE HOR_INTERP_4PTS_2D
