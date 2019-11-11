!     ######spl
MODULE MODI_INTERPXYZ
INTERFACE
!     #####################################################################
      SUBROUTINE INTERPXYZ(PAX,PAY,PAZ,PCHAMP,        &
                           PX,PY,PZ,                  &
                           PXOR,PYOR,PDX,PDY,         &
                           PZL,OTRAJ_GROUP,           &
                           PRESX,PRESY,PRESZ,PRESCHAMP)
!     #####################################################################
!
!
! entrees
!
REAL, DIMENSION(:,:,:),    INTENT(IN)     :: PAX,PAY,PAZ,PCHAMP
                                                                 !
                                                                 !
                                                                 !
REAL,                      INTENT(INOUT)     :: PX,PY,PZ            !
REAL,                      INTENT(IN)     :: PXOR,PYOR,PDX,PDY   !
REAL, DIMENSION(:,:,:),    INTENT(IN)     :: PZL                 !
LOGICAL,                   INTENT(IN)     :: OTRAJ_GROUP
!
! sorties
!
REAL,                      INTENT(OUT)    :: PRESX,PRESY,PRESZ,PRESCHAMP
!
!
END SUBROUTINE INTERPXYZ
!
END INTERFACE
!
END MODULE MODI_INTERPXYZ
!     ######spl
      SUBROUTINE INTERPXYZ(PAX,PAY,PAZ,PCHAMP,        &
                           PX,PY,PZ,                  &
                           PXOR,PYOR,PDX,PDY,         &
                           PZL,OTRAJ_GROUP,           &
                           PRESX,PRESY,PRESZ,PRESCHAMP)
!     #####################################################################
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BUT DE LA ROUTINE : interpoler les trois champs (3D) LG?M 
! (ou trois champs 3D quelconques ecrits sur les points de masse) 
! en un point M, de coordonnees cartesiennes (x,y,z) 
! a priori non-situe sur un point de grille. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           
!
!
!
!!!!!!!!!!!!!!!!!!!!!!
! Declarations
!!!!!!!!!!!!!!!!!!!!!!
!
IMPLICIT NONE
!
! entrees
!
REAL, DIMENSION(:,:,:),    INTENT(IN)     :: PAX,PAY,PAZ,PCHAMP
                                                                 !
                                                                 !
                                                                 !
REAL,                      INTENT(INOUT)     :: PX,PY,PZ            !
REAL,                      INTENT(IN)     :: PXOR,PYOR,PDX,PDY   !
REAL, DIMENSION(:,:,:),    INTENT(IN)     :: PZL                 !
LOGICAL,                   INTENT(IN)     :: OTRAJ_GROUP
!
! sorties
!
REAL,                      INTENT(OUT)    :: PRESX,PRESY,PRESZ,PRESCHAMP
!
! locales
!
INTEGER                              :: II,IJ,IK,JK             !
INTEGER                              :: IKU                     !
REAL                                 :: ZEPS1,ZEPS2,ZEPS3       !
REAL                                 :: ZXREL,ZYREL             !
REAL, DIMENSION(SIZE(PZL,3))         :: ZZLXY                   !
!
!
! initialisations des variables locales
!
IKU=SIZE(PZL,3)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Recherche de la maille contenant le point M(PX,PY,PZ) -> II,IJ,IK
!    Position de M au sein de la maille                    -> ZEPS1,ZEPS2,ZEPS3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! 1.a partie horizontale
!
ZXREL=(PX-PXOR)/PDX+2
ZYREL=(PY-PYOR)/PDY+2
!
II=FLOOR(ZXREL)
IJ=FLOOR(ZYREL)
!
ZEPS1=ZXREL-REAL(II)
ZEPS2=ZYREL-REAL(IJ)
!
!
! 1.b partie verticale
!
! 1.b.1 altitude des niveaux du modele sur la verticale (PX,PY)
!
DO JK=1,IKU
  ZZLXY(JK)=ZEPS2*(ZEPS1*(PZL(II+1,IJ+1,JK))+(1-ZEPS1)*(PZL(II,IJ+1,JK)))     &
             + (1-ZEPS2)*(ZEPS1*(PZL(II+1,IJ,JK))+(1-ZEPS1)*(PZL(II,IJ,JK)))
ENDDO
!
IK=999
DO JK=2,IKU
  IF (ZZLXY(JK).GE.PZ) THEN
    IK=JK-1
    EXIT 
  ENDIF
ENDDO
!
IF (IK==1) THEN
  print *,'la particule est sous le sol'
  print *,' on la remonte a zs + dz/2 = ', ZZLXY(2)
  PZ=ZZLXY(2)
ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!Emergency exit!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (IK==999) THEN
   PRINT*,'PROBLEME AU POINT',II,IJ
   PRINT*,'XREL, YREL, Z =',ZXREL,ZYREL,PZ
   PRINT*,'ZZLXY(IKU)',ZZLXY(IKU)
   STOP
END IF   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ZEPS3=(PZ-ZZLXY(IK))/(ZZLXY(IK+1)-ZZLXY(IK))
!
!------------------------------------------------------------------------------
!
!*    2. INTERPOLATION DES CHAMPS
!
PRESX=  ZEPS3 *                                                             &
      (  ZEPS2*(ZEPS1*(PAX(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PAX(II,IJ+1,IK+1)))  &
       + (1-ZEPS2)*(ZEPS1*(PAX(II+1,IJ,IK+1))+(1-ZEPS1)*(PAX(II,IJ,IK+1)))  &
      )                                                                     &    
      + (1-ZEPS3) *                                                         &
      (  ZEPS2*(ZEPS1*(PAX(II+1,IJ+1,IK))+(1-ZEPS1)*(PAX(II,IJ+1,IK)))      &
       + (1-ZEPS2)*(ZEPS1*(PAX(II+1,IJ,IK))+(1-ZEPS1)*(PAX(II,IJ,IK)))      &
      )
!
PRESY=  ZEPS3 *                                                             &
      (  ZEPS2*(ZEPS1*(PAY(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PAY(II,IJ+1,IK+1)))  &
       + (1-ZEPS2)*(ZEPS1*(PAY(II+1,IJ,IK+1))+(1-ZEPS1)*(PAY(II,IJ,IK+1)))  &
      )                                                                     &    
      + (1-ZEPS3) *                                                         &
      (  ZEPS2*(ZEPS1*(PAY(II+1,IJ+1,IK))+(1-ZEPS1)*(PAY(II,IJ+1,IK)))      &
       + (1-ZEPS2)*(ZEPS1*(PAY(II+1,IJ,IK))+(1-ZEPS1)*(PAY(II,IJ,IK)))      &
      )
!
PRESZ=  ZEPS3 *                                                             &
      (  ZEPS2*(ZEPS1*(PAZ(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PAZ(II,IJ+1,IK+1)))  &
       + (1-ZEPS2)*(ZEPS1*(PAZ(II+1,IJ,IK+1))+(1-ZEPS1)*(PAZ(II,IJ,IK+1)))  &
      )                                                                     &    
      + (1-ZEPS3) *                                                         &
      (  ZEPS2*(ZEPS1*(PAZ(II+1,IJ+1,IK))+(1-ZEPS1)*(PAZ(II,IJ+1,IK)))      &
       + (1-ZEPS2)*(ZEPS1*(PAZ(II+1,IJ,IK))+(1-ZEPS1)*(PAZ(II,IJ,IK)))      &
      )
IF (OTRAJ_GROUP) THEN
  PRESCHAMP=  ZEPS3 *                                                         &
        (  ZEPS2*(ZEPS1*(PCHAMP(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PCHAMP(II,IJ+1,IK+1)))  &
         + (1-ZEPS2)*(ZEPS1*(PCHAMP(II+1,IJ,IK+1))+(1-ZEPS1)*(PCHAMP(II,IJ,IK+1)))  &
        )                                                                     &    
        + (1-ZEPS3) *                                                         &
        (  ZEPS2*(ZEPS1*(PCHAMP(II+1,IJ+1,IK))+(1-ZEPS1)*(PCHAMP(II,IJ+1,IK)))      &
         + (1-ZEPS2)*(ZEPS1*(PCHAMP(II+1,IJ,IK))+(1-ZEPS1)*(PCHAMP(II,IJ,IK)))      &
        )
ENDIF
!
!------------------------------------------------------------------------------
!
!
END SUBROUTINE INTERPXYZ
