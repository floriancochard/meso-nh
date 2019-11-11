!     ######spl
      SUBROUTINE TRACIRCLE(PU,PV,PP,PLW)
!     ###################################
!
!!****  *TRACIRCLE* - 
!!
!!    PURPOSE
!!    -------
!        Trace de cercles concentriques (pour materialiser par ex la
!      portee de radar(s))     
!
!!**  METHOD
!!    ------
!!      L utilisateur fournit :
!!    Le centre du cercle  en latitude / longitude et
!!    son(ses) rayon(s) en metres 
!!    Conversion en coordonnees normalisees et trace des segments successifs
!!    du(des) cercle(s) 
!!
!!    EXTERNAL
!!    --------
!!      SET      : defines NCAR window and viewport in normalized and user
!!                 coordinates
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_RADAR 
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
!!      Original       23/04/03
!!      Updated   PM   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RADAR
!
IMPLICIT NONE
!
REAL :: ZWL, ZWR, ZWB, ZWT
REAL :: ZVL, ZVR, ZVB, ZVT
REAL :: PU, PV      ! Coord. conformes centre du cercle
REAL :: PP, PLW     ! Rayon et epaisseur du trai du cercle
REAL :: ZXCN, ZYCN  ! coord normalisees du centre du cercle
REAL :: ZRN         ! Rayon en coord normalisees <-> PP
REAL :: ZXA, ZYA, ZDTR, ZANG, ZSINA, ZCOSA, ZXB, ZYB, ZWIDTH, ZPPKM
REAL :: ZX30, ZY30, ZX60,ZY60, ZX90,ZY90, ZX120,ZY120, ZX150,ZY150,&
ZX180,ZY180, ZX210,ZY210, ZX240,ZY240, ZX270,ZY270,ZX300,ZY300, ZX330,ZY330,&
ZX360,ZY360
INTEGER :: ID, IER, J 
CHARACTER(LEN=4) :: YC
!
!-------------------------------------------------------------------------------
!
!*       1.    SAUVEGARDE FENETRE COURANTE
!              ---------------------------
!
CALL GETSET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
CALL GQLWSC(IER,ZWIDTH)
!
!*       2.    CALCUL DE COORDONNEES NORMALISEES et DEF. NOUVELLE FENETRE
!              ----------------------------------------------------------
! Calcul des coordonnees normalisees du centre du cercle et de la dim du rayon
ZXCN=ZVL+((PU-ZWL)*(ZVR-ZVL)/(ZWR-ZWL))
ZYCN=ZVB+((PV-ZWB)*(ZVT-ZVB)/(ZWT-ZWB))
ZRN=PP*(ZVR-ZVL)/(ZWR-ZWL)
CALL SET(ZVL,ZVR,ZVB,ZVT,ZVL,ZVR,ZVB,ZVT,1)
!
!*       3.    TRACE DU CERCLE
!              ---------------
!
CALL SFLUSH
IF(PLW == 0. .OR. PLW == 9999.)PLW=2.
CALL GSLWSC(PLW)
ZXA=ZXCN
ZYA=ZYCN+ZRN
CALL FRSTPT(ZXA,ZYA)
ZDTR=3.141592654/180.
CALL GSCLIP(1)
DO J=1,360
  ZANG=J*ZDTR
  ZSINA=SIN(ZANG)
  ZCOSA=COS(ZANG)
  IF(J == 90)ZSINA=1.
  IF(J == 90)ZCOSA=0.
  IF(J == 360)ZSINA=0.
  IF(J == 360)ZCOSA=1.
  ZXB=ZRN*ZSINA
  ZYB=ZRN*ZCOSA
  ZXB=ZXCN+ZXB
  ZYB=ZYCN+ZYB
  CALL VECTOR(ZXB,ZYB)
  ZXA=ZXB
  ZYA=ZYB
  IF(LRADIST)THEN
  ZPPKM=PP/1000.
  WRITE(YC,'(I4)')NINT(ZPPKM)
  YC=ADJUSTL(YC)
  IF(J == 90 .OR. J == 270)THEN
  CALL GSCLIP(0)
    IF(J == 90 .AND. ZXA > ZVR)THEN
    ELSE
      CALL PLCHHQ(ZXA,ZYA,YC(1:LEN_TRIM(YC)),.008,0.,0.)
    ENDIF
  CALL GSCLIP(1)
  ELSEIF(J == 180)THEN
    CALL PLCHHQ(ZXA,ZYA-.005,YC(1:LEN_TRIM(YC)),.008,0.,0.)
  ELSEIF(J == 360)THEN
    CALL PLCHHQ(ZXA,ZYA+.005,YC(1:LEN_TRIM(YC)),.008,0.,0.)
  ENDIF
  ENDIF
  IF(J == 30)THEN
    ZX30=ZXA; ZY30=ZYA
  ELSEIF(J == 60)THEN
    ZX60=ZXA; ZY60=ZYA
  ELSEIF(J == 90)THEN
    ZX90=ZXA; ZY90=ZYA
  ELSEIF(J == 120)THEN
    ZX120=ZXA; ZY120=ZYA
  ELSEIF(J == 150)THEN
    ZX150=ZXA; ZY150=ZYA
  ELSEIF(J == 180)THEN
    ZX180=ZXA; ZY180=ZYA
  ELSEIF(J == 210)THEN
    ZX210=ZXA; ZY210=ZYA
  ELSEIF(J == 240)THEN
    ZX240=ZXA; ZY240=ZYA
  ELSEIF(J == 270)THEN
    ZX270=ZXA; ZY270=ZYA
  ELSEIF(J == 300)THEN
    ZX300=ZXA; ZY300=ZYA
  ELSEIF(J == 330)THEN
    ZX330=ZXA; ZY330=ZYA
  ELSEIF(J == 360)THEN
    ZX360=ZXA; ZY360=ZYA
  ENDIF
ENDDO
CALL SFLUSH
CALL GSCLIP(1)
!
!*       4.    TRACE DES RAYONS
!              ----------------
!
IF(LRADRAY)THEN
  CALL SFLUSH
  CALL GSLN(2)
  CALL GSLWSC(2.)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX30,ZY30)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX60,ZY60)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX90,ZY90)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX120,ZY120)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX150,ZY150)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX180,ZY180)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX210,ZY210)
  CALL SFLUSH
  CALL GSLN(2)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX240,ZY240)
  CALL SFLUSH
  CALL GSLN(2)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX270,ZY270)
  CALL SFLUSH
  CALL GSLN(2)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX300,ZY300)
  CALL SFLUSH
  CALL GSLN(2)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX330,ZY330)
  CALL SFLUSH
  CALL GSLN(2)
  CALL FRSTPT(ZXCN,ZYCN)
  CALL VECTOR(ZX360,ZY360)
  CALL SFLUSH
ENDIF
!
CALL GSCLIP(0)
!
!*       5.    RESTORATION FENETRE COURANTE
!              ----------------------------
CALL SET(ZVL,ZVR,ZVB,ZVT,ZWL,ZWR,ZWB,ZWT,ID)
CALL GSLWSC(ZWIDTH)
CALL GSLN(1)

!
!*       6.   EXIT
!            ----
!
RETURN
END SUBROUTINE  TRACIRCLE
