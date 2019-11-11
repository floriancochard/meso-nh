!----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:./s.int2lalo.f90, Version:1.8, Date:03/06/05, Last modified:01/10/15
!-----------------------------------------------------------------
!     ######spl
MODULE MODI_INT2LALO
!###################
!
INTERFACE
      SUBROUTINE INT2LALO(HHORTYPE,P3D,PLATLON,PSVAL,PLALO)
!
CHARACTER(LEN=4),     INTENT(IN) :: HHORTYPE ! type of horizontal interpolation
REAL,DIMENSION(:,:,:),INTENT(IN) :: P3D  ! input 3d array s->n, w->e
REAL,DIMENSION(4),    INTENT(IN) :: PLATLON  ! NSWE target domain bounds (milliDEGS)
REAL,                 INTENT(IN) :: PSVAL    ! value for missing data
REAL,DIMENSION(:,:,:),INTENT(OUT):: PLALO    ! output interpolated LAT/LON field
!
END SUBROUTINE INT2LALO
END INTERFACE
END MODULE MODI_INT2LALO
!
!####################################################
SUBROUTINE INT2LALO(HHORTYPE,P3D,PLATLON,PSVAL,PLALO)
!####################################################
!
!!    PURPOSE
!!    -------
!       Interpolates data from a conformal grid to a lat/lon grid
!
!!**  METHOD
!!    ------
!!       Input is the conformal data (S->N scanning) and lat/lon domain 
!!      definition.
!!       Output is the lat/lon data in N->S scanning (required).
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,        ONLY : XRADIUS,XPI
USE MODD_GRID,       ONLY : XLONORI,XLATORI
USE MODD_PARAMETERS, ONLY : XUNDEF,JPHEXT
USE MODD_DIM1,       ONLY : NIMAX,NJMAX,NKMAX
USE MODD_GRID1,      ONLY : XXHAT,XYHAT
!
USE MODE_GRIDPROJ
!
IMPLICIT NONE
!
!*       0.1   Arguments
!
CHARACTER(LEN=4),     INTENT(IN) :: HHORTYPE ! type of horizontal interpolation
REAL,DIMENSION(:,:,:),INTENT(IN) :: P3D      ! input 3d array s->n, w->e
REAL,DIMENSION(4),    INTENT(IN) :: PLATLON  ! NSWE target domain bounds (milliDEGS)
REAL,                 INTENT(IN) :: PSVAL    ! value for missing data
REAL,DIMENSION(:,:,:),INTENT(OUT):: PLALO    ! output interpolated LAT/LON field
                                             !with N->S scanning
!
!*       0.2   Local variables
!
REAL :: ZLONW,ZLONE,ZLATN,ZLATS   !(degres)
REAL :: ZXHAT,ZYHAT
REAL :: ZDX,ZDY                  ! TARGET INCREMENTS IN LAT/LON
REAL :: ZLA,ZLO,ZAB,ZCD,ZI,ZJ,ZXR,ZYR
REAL :: ZEPS
INTEGER :: JX,JY,JK,INX,INY,II,IJ,IK,IM,IN
!
!------------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
!
ZEPS=1.E-10
PLALO(:,:,:)= PSVAL
INX=SIZE(PLALO,1) ; INY=SIZE(PLALO,2) ; IK=SIZE(PLALO,3)
ZLONW=PLATLON(3)/1000.  ; ZLONE=PLATLON(4)/1000. ; ZLATN=PLATLON(1)/1000. ; ZLATS=PLATLON(2)/1000.
!
ZDX=(ZLONE-ZLONW)/(INX-1)
IF (ZDX<0) ZDX=(ZLONE-ZLONW+360.)/(INX-1)
ZDY=(ZLATN-ZLATS)/(INY-1)
print*, 'INT2LALO: target increments:',ZDX,ZDY
PRINT*, 'INT2LALO: target increments:',ZDX,ZDY
!
!------------------------------------------------------------------------------
!
!*       2.    INTERPOLATION
!              -------------
!
!print*,'av interp.: ',minval(P3D),minloc(p3d),maxval(p3d),maxloc(p3d)
DO JK=1,IK
  DO JY=1,INY ; DO JX=1,INX
    ZLO=MOD(ZLONW+ZDX*(JX-1),360.)
    ZLA=ZLATN-ZDY*(JY-1)          ! output has N->S scanning
 !   print*,ZLO,ZLA,JX,JY
    CALL SM_XYHAT(XLATORI,XLONORI,ZLA,ZLO,ZXHAT,ZYHAT)
    II=MAX(MIN(COUNT(XXHAT(:)<ZXHAT),NIMAX+JPHEXT),1+JPHEXT)
    IJ=MAX(MIN(COUNT(XYHAT(:)<ZYHAT),NJMAX+JPHEXT),1+JPHEXT)
    ZI=(ZXHAT-XXHAT(II))/(XXHAT(II+1)-XXHAT(II))+FLOAT(II)
    ZJ=(ZYHAT-XYHAT(IJ))/(XYHAT(IJ+1)-XYHAT(IJ))+FLOAT(IJ)
    !
!!!!!!!!!!!!!!! PLUSE DE DECALAGE D INDICES ENTRE P3D ET LE TABLEAUX MNH!!!!!!!!! 
    IM=II
    IN=IJ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    IF (HHORTYPE=='NEAR') THEN 
    ! NEARest neighbour method on conformal plane
      IF (      (II>=1+JPHEXT) .AND. (II<=NIMAX+JPHEXT) &
         .AND.  (IJ>=1+JPHEXT) .AND. (IJ<=NJMAX+JPHEXT) ) THEN
        PLALO(JX,JY,JK)=P3D(IM,IN,JK) ! take nearest-neighbour value
      ENDIF
    !
    ELSEIF (HHORTYPE == 'BILI') THEN
    ! LInear interpolation method on conformal plane
      IF (      (ZI>=1+JPHEXT) .AND. (ZI<=NIMAX+JPHEXT) &
         .AND.  (ZJ>=1+JPHEXT) .AND. (ZJ<=NJMAX+JPHEXT) ) THEN
        IF (ALL(ABS(P3D(IM:IM+1,IN:IN+1,JK)-XUNDEF)>=ZEPS) ) THEN
        ! take the 4 surrounding values and apply bilinear interpolation
          ZXR=ZI-REAL(II) ; ZYR=ZJ-REAL(IJ)  ! coordinates inside rectangle
          ZAB= (1.-ZXR)*P3D(IM,  IN,  JK)  &
              +    ZXR *P3D(IM+1,IN,  JK)
          ZCD= (1.-ZXR)*P3D(IM,  IN+1,JK)  &
              +    ZXR *P3D(IM+1,IN+1,JK)
          PLALO(JX,JY,JK)= (1.-ZYR)*ZAB + ZYR*ZCD
        ENDIF
      ENDIF
    ELSE
      print*, 'Horizontal type interpolation unknown ',HHORTYPE
    ENDIF
  ENDDO ; ENDDO 
ENDDO
!
!------------------------------------------------------------------------------
!
!*       3.    EXTENSION
!              ---------
!
!IF (IK/=1) CALL EXTENDLAM
!print*,'ap interp.: ',minval(Plalo),minloc(plalo),maxval(plalo),maxloc(plalo)
!
RETURN
CONTAINS 
!-----------------------------
SUBROUTINE EXTENDLAM
!  PURPOSE: EXTEND AN INTERPOLATED LAT/LON FIELD OUTSIDE THE kAl MODEL
!           DOMAIN BY REMOVING ALL ITS UNDEFINED VALUES.
!  METHOD: REPLACED ALL UNDEFINED VALUES BY AVERAGE OF DEFINED VALUES.
!
REAL ZS
INTEGER IPOP,JI,JJ,JK
! 
! COMPUTE AVERAGE OF DEFINED VALUES
ZS=0. ; IPOP=0
DO JK=1,IK
  DO JJ=1,INY ; DO JI=1,INX
    IF(PLALO(JI,JJ,JK)/=PSVAL)THEN
      ZS=ZS+PLALO(JI,JJ,JK)
      IPOP=IPOP+1
    ENDIF
  ENDDO ; ENDDO
ENDDO
ZS=ZS/(FLOAT(IPOP)+TINY(ZS))
!
! Replace ALL UNDEFINED VALUES BY THE AVERAGE
DO JK=1,IK
  DO JJ=1,INY ; DO JI=1,INX
    IF(PLALO(JI,JJ,JK)==PSVAL) PLALO(JI,JJ,JK)=ZS
  ENDDO ; ENDDO
ENDDO
!  
RETURN
END SUBROUTINE EXTENDLAM

END SUBROUTINE INT2LALO
