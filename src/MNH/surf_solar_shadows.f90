!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 param 2006/05/18 13:07:25
!-----------------------------------------------------------------
!    ##############################
     MODULE MODI_SURF_SOLAR_SHADOWS 
!    ##############################
!
INTERFACE 
!
      SUBROUTINE SURF_SOLAR_SHADOWS ( PMAP, PXHAT, PYHAT,                     &
                  PCOSZEN, PSINZEN, PAZIMSOL,                                 &
                  PZS, PZS_XY, PDIRSWDT, PDIRSRFSWD                           )
!
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PMAP         ! map factor
REAL, DIMENSION(:),       INTENT(IN) :: PXHAT        ! X coordinate
REAL, DIMENSION(:),       INTENT(IN) :: PYHAT        ! Y coordinate
REAL, DIMENSION(:,:),     INTENT(IN) :: PCOSZEN ! COS(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PSINZEN ! SIN(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PAZIMSOL! azimuthal solar angle
                                                ! (radian, from North, clockwise)
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS        ! (resolved) model orography
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS_XY     ! orography at vort. points
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PDIRSWDT   ! SW Flux received by
                                                     ! each subgrid triangle
REAL, DIMENSION(:,:,:),   INTENT(OUT)  :: PDIRSRFSWD ! SuRF. DIRect SW Flux
!
END SUBROUTINE SURF_SOLAR_SHADOWS  
!
END INTERFACE
!
END MODULE MODI_SURF_SOLAR_SHADOWS
!     #########################################################################
      SUBROUTINE SURF_SOLAR_SHADOWS ( PMAP, PXHAT, PYHAT,                     &
                  PCOSZEN, PSINZEN, PAZIMSOL,                                 &
                  PZS, PZS_XY, PDIRSWDT, PDIRSRFSWD                           )
!     #########################################################################
!
!!****  * SURF_SOLAR_SHADOWS * - computes the modifications to the downwards
!!                           direct solar flux at the surface, due to
!!                           orientation and shape of this surface.
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
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
!!
!!    AUTHOR
!!    ------
!!	V. Masson      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/01/02
!!      V. Masson   01/03/03 add multiple wavelenghts
!!      V. Masson   04/01/12 standard definition of azimuthal angle
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS, ONLY : XUNDEF, JPHEXT
USE MODD_CST,        ONLY : XPI, XRADIUS
USE MODD_CONF,       ONLY : LCARTESIAN
USE MODD_SHADOWS_n,  ONLY : XZS_ll, XZS_XY_ll, XXHAT_ll, &
                            XYHAT_ll, XZS_MAX_ll, XZS_MAX_ll
!
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
!
REAL, DIMENSION(:,:),     INTENT(IN) :: PMAP         ! map factor
REAL, DIMENSION(:),       INTENT(IN) :: PXHAT        ! X coordinate
REAL, DIMENSION(:),       INTENT(IN) :: PYHAT        ! Y coordinate
REAL, DIMENSION(:,:),     INTENT(IN) :: PCOSZEN ! COS(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PSINZEN ! SIN(zenithal solar angle)
REAL, DIMENSION(:,:),     INTENT(IN) :: PAZIMSOL! azimuthal solar angle
                                                ! (radian, from North, clockwise)
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS        ! (resolved) model orography
REAL, DIMENSION(:,:),     INTENT(IN) :: PZS_XY     ! orography at vort. points
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PDIRSWDT   ! SW Flux received by
                                                     ! each subgrid triangle
REAL, DIMENSION(:,:,:),   INTENT(OUT)  :: PDIRSRFSWD ! SuRF. DIRect SW Flux
!
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: JI, JJ
INTEGER :: IRESP
REAL    :: ZAZIM    ! azimuthal solar angle
REAL    :: ZCOSAZIM ! azimuthal solar angle cosine
REAL    :: ZSINAZIM ! azimuthal solar angle sine
REAL    :: ZCOSZEN  ! zenithal solar angle cosine
REAL    :: ZSINZEN  ! zenithal solar angle sine
INTEGER :: JLOOP
INTEGER :: JT
!
REAL                   :: ZX       ! current X, Y and Z coordinates of
REAL                   :: ZY       ! the sun beam
REAL                   :: ZZ       !
REAL                   :: ZA       ! X, Y coordinates of a point
REAL                   :: ZB       ! projected according to sun dir.
REAL                   :: ZDX      ! current distance of the sun beam 
REAL                   :: ZDY      ! to the next X or Y side box limit
REAL                   :: ZZI      ! initial position of the sun beam
REAL                   :: ZZCURV   ! offset due to earth curve surface
REAL, DIMENSION(3)     :: ZXT  ! X, Y and Z coordinates of corners
REAL, DIMENSION(3)     :: ZYT  ! of a triangle
REAL, DIMENSION(3)     :: ZZT  !
REAL, DIMENSION(3)     :: ZAT  ! X, Y coordinates of projected corners
REAL, DIMENSION(3)     :: ZBT  ! of a triangle
LOGICAL                :: GF   ! interception flag
!
!*       0.3   DECLARATIONS OF LOCAL VARIABLES FOR ENTIRE PARALLELIZED DOMAIN
!
INTEGER :: IIU_ll, IJU_ll
INTEGER :: IIB_ll, IIE_ll, IJB_ll, IJE_ll
INTEGER :: IIMAX_ll, IJMAX_ll
INTEGER :: JI_ll, JJ_ll               ! current grid mesh
INTEGER :: II_ll, IJ_ll               ! current intercepting grid mesh
INTEGER :: IIOR_ll, IJOR_ll           ! origin of current processor
!                                     ! in entire domain
!
!
!-------------------------------------------------------------------------------
!
!*       1.    DIMENSION INITIALIZATIONS
!              -------------------------
!
CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
IIU_ll = IIMAX_ll+2*JPHEXT
IJU_ll = IJMAX_ll+2*JPHEXT
IIB_ll = 1       +  JPHEXT
IIE_ll = IIU_ll  -  JPHEXT
IJB_ll = 1       +  JPHEXT
IJE_ll = IJU_ll  -  JPHEXT
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
CALL GET_OR_ll('B',IIOR_ll,IJOR_ll)
!
ZZCURV=0.
!
PDIRSRFSWD=XUNDEF
!
!-------------------------------------------------------------------------------
!
!*       2.    LOOP ON THIS PROCESSOR GRID MESH
!              --------------------------------
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
!
!* If zenithal angle greater than Pi/2, sun is down.
!
    IF (PCOSZEN(JI,JJ)<0.) CYCLE
!
!
!* If zenithal angle is vertical, there is no shadow
!
    IF (PCOSZEN(JI,JJ)==1.) CYCLE
!
!-------------------------------------------------------------------------------
!
!*       3.    DEFINITION OF SOLAR ANGLES
!              --------------------------
!
!
!* zenithal angle between 0 and Pi/2. Equals zero for sun at zenith.
!
    ZCOSZEN =PCOSZEN(JI,JJ)
!
!* Azimuthal angle is set between 0 and 2*Pi, being equal to zero
!  in the East direction, Pi/2 in the North direction, etc...
!
    ZSINZEN =PSINZEN(JI,JJ)
    ZAZIM=PAZIMSOL(JI,JJ)
    !
!
    ZAZIM = - ZAZIM + XPI/2.
!
!* cosine and sine of azimuthal solar angle
!
    ZCOSAZIM=COS(ZAZIM)
    ZSINAZIM=SIN(ZAZIM)
!
!* indices of grid point relative to the entire domain
!
    JI_ll = JI + IIOR_ll - 1
    JJ_ll = JJ + IJOR_ll - 1
!
!-------------------------------------------------------------------------------
!
!*       4.    EXPLORATION OF INTERCEPTING OROGRAPHY FOLLOWING SUN BEAM
!              --------------------------------------------------------
!
    DO JT=1,4
!
!* slope already in its own shadow
!
      IF (ALL(PDIRSWDT(JI,JJ,JT,:)==0.)) CYCLE
!
!* coordinates of the point in the center of the triangle
!
      SELECT CASE (JT)
       CASE (1)
         ZX=(5.*PXHAT(JI)+PXHAT(JI+1))/6.
         ZY=0.5*(PYHAT(JJ)+PYHAT(JJ+1))
         ZZ=(PZS(JI,JJ)+PZS_XY(JI,JJ)+PZS_XY(JI,JJ+1))/3.
       CASE (2)
         ZX=0.5*(PXHAT(JI)+PXHAT(JI+1))
         ZY=(5.*PYHAT(JJ+1)+PYHAT(JJ))/6.
         ZZ=(PZS(JI,JJ)+PZS_XY(JI,JJ+1)+PZS_XY(JI+1,JJ+1))/3.
       CASE (3)
         ZX=(5.*PXHAT(JI+1)+PXHAT(JI))/6.
         ZY=0.5*(PYHAT(JJ)+PYHAT(JJ+1))
         ZZ=(PZS(JI,JJ)+PZS_XY(JI+1,JJ)+PZS_XY(JI+1,JJ+1))/3.
       CASE (4)
         ZX=0.5*(PXHAT(JI)+PXHAT(JI+1))
         ZY=(5.*PYHAT(JJ)+PYHAT(JJ+1))/6.
         ZZ=(PZS(JI,JJ)+PZS_XY(JI,JJ)+PZS_XY(JI+1,JJ))/3.
      END SELECT
!
!* projection of this point according to sun direction
!
      CALL PROJ_SOLAR(ZX,ZY,ZZ,ZA,ZB)
!
!* starts exploration at this point
!
      II_ll = JI_ll
      IJ_ll = JJ_ll
!
      ZZI = ZZ
!
!
!* loop following sun beam until interception
!
      DO JLOOP=1,2*(IIU_ll+IJU_ll)
!
!* go to the next grid point encountered by the sun beam
!
        IF (ZSINAZIM>0. .AND. ZCOSAZIM>0.) THEN
          ZDY=(XYHAT_ll(IJ_ll+1)-ZY)/ZSINAZIM
          ZDX=(XXHAT_ll(II_ll+1)-ZX)/ZCOSAZIM
        ELSE IF (ZSINAZIM<0. .AND. ZCOSAZIM>0.) THEN
          ZDY=(XYHAT_ll(IJ_ll)  -ZY)/ZSINAZIM
          ZDX=(XXHAT_ll(II_ll+1)-ZX)/ZCOSAZIM
        ELSE IF (ZSINAZIM>0. .AND. ZCOSAZIM<0.) THEN
          ZDY=(XYHAT_ll(IJ_ll+1)-ZY)/ZSINAZIM
          ZDX=(XXHAT_ll(II_ll)  -ZX)/ZCOSAZIM
        ELSE IF (ZSINAZIM<0. .AND. ZCOSAZIM<0.) THEN
          ZDY=(XYHAT_ll(IJ_ll)  -ZY)/ZSINAZIM
          ZDX=(XXHAT_ll(II_ll)  -ZX)/ZCOSAZIM
        ELSE IF (ZSINAZIM==0. .AND. ZCOSAZIM<0.) THEN
          ZDX=-(XXHAT_ll(II_ll)  -ZX)
          ZDY=2.*ZDX
        ELSE IF (ZSINAZIM==0. .AND. ZCOSAZIM>0.) THEN
          ZDX=(XXHAT_ll(II_ll+1)-ZX)
          ZDY=2.*ZDX
        ELSE IF (ZSINAZIM<0. .AND. ZCOSAZIM==0.) THEN
          ZDY=-(XYHAT_ll(IJ_ll)  -ZY)
          ZDX=2.*ZDY
        ELSE IF (ZSINAZIM>0. .AND. ZCOSAZIM==0.) THEN
          ZDY=(XYHAT_ll(IJ_ll+1)-ZY)
          ZDX=2.*ZDY
        END IF
!
        IF (ZDY<ZDX) THEN
          ZX = ZX + ZDY * ZCOSAZIM
          ZY = ZY + ZDY * ZSINAZIM
          ZZ = ZZ + ZDY * ZCOSZEN / PMAP(JI,JJ) / ZSINZEN
          IF (ZSINAZIM>0.) THEN
            IJ_ll = IJ_ll + 1
          ELSE
            IJ_ll = IJ_ll - 1
          END IF
        ELSE
          ZX = ZX + ZDX * ZCOSAZIM
          ZY = ZY + ZDX * ZSINAZIM
          ZZ = ZZ + ZDX * ZCOSZEN / PMAP(JI,JJ) / ZSINZEN
          IF (ZCOSAZIM>0.) THEN
            II_ll = II_ll + 1
          ELSE
            II_ll = II_ll - 1
          END IF
        END IF
!
!* effect of curvature of the earth surface
!
        IF (.NOT. LCARTESIAN) ZZCURV =  ((ZZ-ZZI)*ZSINZEN/ZCOSZEN)**2 &
                                      / (2.*XRADIUS)
!
!* sun beam goes to high
!
        IF (ZZ > XZS_MAX_ll-ZZCURV) EXIT
!
!* sun beam goes outside of the domain
!
        IF (     II_ll<IIB_ll .OR. II_ll>IIE_ll &
            .OR. IJ_ll<IJB_ll .OR. IJ_ll>IJE_ll ) EXIT
!
!
!* treatment where the sun beam is currently located
!
        IF (.NOT. ( XZS_ll   (II_ll  ,IJ_ll  ) -ZZCURV < ZZ .AND. &
                    XZS_XY_ll(II_ll  ,IJ_ll  ) -ZZCURV < ZZ .AND. &
                    XZS_XY_ll(II_ll+1,IJ_ll  ) -ZZCURV < ZZ .AND. &
                    XZS_XY_ll(II_ll  ,IJ_ll+1) -ZZCURV < ZZ .AND. &
                    XZS_XY_ll(II_ll+1,IJ_ll+1) -ZZCURV < ZZ  )    ) THEN
!
!* exploration of the 4 triangles in the grid mesh encoutered by the sun beam
!
!* left triangle
!
          ZXT(1) =0.5*(XXHAT_ll(II_ll)+XXHAT_ll(II_ll+1))
          ZYT(1) =0.5*(XYHAT_ll(IJ_ll)+XYHAT_ll(IJ_ll+1))
          ZZT(1) =XZS_ll   (II_ll          ,IJ_ll          ) - ZZCURV
          ZXT(2) =XXHAT_ll(II_ll)
          ZYT(2) =XYHAT_ll(IJ_ll)
          ZZT(2) =XZS_XY_ll(II_ll  ,IJ_ll  ) - ZZCURV
          ZXT(3) =XXHAT_ll(II_ll)
          ZYT(3) =XYHAT_ll(IJ_ll+1)
          ZZT(3) =XZS_XY_ll(II_ll  ,IJ_ll+1) - ZZCURV
          CALL PROJ_SOLAR(ZXT(1),ZYT(1),ZZT(1),ZAT(1),ZBT(1))
          CALL PROJ_SOLAR(ZXT(2),ZYT(2),ZZT(2),ZAT(2),ZBT(2))
          CALL PROJ_SOLAR(ZXT(3),ZYT(3),ZZT(3),ZAT(3),ZBT(3))
          CALL SOLAR_INTERC(ZAT(:),ZBT(:),ZA,ZB,GF)
          IF (GF) THEN
            PDIRSWDT(JI,JJ,JT,:)=0.
            EXIT
          END IF
!
!* top triangle
!
          ZAT(2) =ZAT(3)
          ZBT(2) =ZBT(3)
          ZXT(3) =XXHAT_ll(II_ll+1)
          ZYT(3) =XYHAT_ll(IJ_ll+1)
          ZZT(3) =XZS_XY_ll(II_ll+1,IJ_ll+1) - ZZCURV
          CALL PROJ_SOLAR(ZXT(3),ZYT(3),ZZT(3),ZAT(3),ZBT(3))
          CALL SOLAR_INTERC(ZAT(:),ZBT(:),ZA,ZB,GF)
          IF (GF) THEN
            PDIRSWDT(JI,JJ,JT,:)=0.
            EXIT
          END IF
!
!* right triangle
!
          ZAT(2) =ZAT(3)
          ZBT(2) =ZBT(3)
          ZXT(3) =XXHAT_ll(II_ll+1)
          ZYT(3) =XYHAT_ll(IJ_ll)
          ZZT(3) =XZS_XY_ll(II_ll+1,IJ_ll  ) - ZZCURV
          CALL PROJ_SOLAR(ZXT(3),ZYT(3),ZZT(3),ZAT(3),ZBT(3))
          CALL SOLAR_INTERC(ZAT(:),ZBT(:),ZA,ZB,GF)
          IF (GF) THEN
            PDIRSWDT(JI,JJ,JT,:)=0.
            EXIT
          END IF
!
!* bottom triangle
!
          ZAT(2) =ZAT(3)
          ZBT(2) =ZBT(3)
          ZXT(3) =XXHAT_ll(II_ll)
          ZYT(3) =XYHAT_ll(IJ_ll)
          ZZT(3) =XZS_XY_ll(II_ll  ,IJ_ll  ) - ZZCURV
          CALL PROJ_SOLAR(ZXT(3),ZYT(3),ZZT(3),ZAT(3),ZBT(3))
          CALL SOLAR_INTERC(ZAT(:),ZBT(:),ZA,ZB,GF)
          IF (GF) THEN
            PDIRSWDT(JI,JJ,JT,:)=0.
            EXIT
          END IF
        END IF
      END DO
!
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       6.    FLUX FOR THE GRID MESH RECEIVED OVER AN HORIZONTAL SURFACE
!              ----------------------------------------------------------
!
!* there is no problem of summation, because only fluxes projected
!  over horizontal surfaces are considered here.
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    PDIRSRFSWD(JI,JJ,:) = (  PDIRSWDT(JI,JJ,1,:) &
                           + PDIRSWDT(JI,JJ,2,:) &
                           + PDIRSWDT(JI,JJ,3,:) &
                           + PDIRSWDT(JI,JJ,4,:))&
                          * 0.25
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
SUBROUTINE PROJ_SOLAR(PX,PY,PZ,PA,PB)
REAL, INTENT(IN)  :: PX  ! X coordinate (conformal plane)
REAL, INTENT(IN)  :: PY  ! Y coordinate (conformal plane)
REAL, INTENT(IN)  :: PZ  ! Z coordinate (conformal plane)
REAL, INTENT(OUT) :: PA  ! X coordinate (plane perp. sun)
REAL, INTENT(OUT) :: PB  ! Y coordinate (plane perp. sun)

PA =   ZCOSZEN*ZCOSAZIM*PX + ZCOSZEN*ZSINAZIM*PY - ZSINZEN*PZ
PB = -         ZSINAZIM*PX +         ZCOSAZIM*PY
!
END SUBROUTINE PROJ_SOLAR
!
SUBROUTINE SOLAR_INTERC(PA,PB,PAD,PBD,OF)
REAL, DIMENSION(3), INTENT(IN)  :: PA  ! X coordinate 
REAL, DIMENSION(3), INTENT(IN)  :: PB  ! Y coordinate 
REAL,               INTENT(IN)  :: PAD ! Z coordinate
REAL,               INTENT(IN)  :: PBD ! X coordinate
LOGICAL,            INTENT(OUT) :: OF  ! true if interception occurs
!
REAL :: ZC1, ZC2, ZC3 ! coefficients of the c1 A + c2 B + c3 = 0 line
!                     ! defined by two corners of the triangle.
!
OF=.FALSE.
!
!* first side
!
ZC1=(PB(1)-PB(2))
ZC2=-(PA(1)-PA(2))
ZC3=-PB(1)*PA(2)+PB(2)*PA(1)
!
IF ( (ZC1*PAD+ZC2*PBD+ZC3)*(ZC1*PA(3)+ZC2*PB(3)+ZC3) <0.) RETURN
!
!* second side
!
ZC1=(PB(3)-PB(2))
ZC2=-(PA(3)-PA(2))
ZC3=-PB(3)*PA(2)+PB(2)*PA(3)
!
IF ( (ZC1*PAD+ZC2*PBD+ZC3)*(ZC1*PA(1)+ZC2*PB(1)+ZC3) <0.) RETURN
!
!* third side
!
ZC1=(PB(3)-PB(1))
ZC2=-(PA(3)-PA(1))
ZC3=-PB(3)*PA(1)+PB(1)*PA(3)
!
IF ( (ZC1*PAD+ZC2*PBD+ZC3)*(ZC1*PA(2)+ZC2*PB(2)+ZC3) <0.) RETURN
!
OF=.TRUE.
!
END SUBROUTINE SOLAR_INTERC
!
END SUBROUTINE SURF_SOLAR_SHADOWS

