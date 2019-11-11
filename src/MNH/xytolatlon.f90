!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ###################################################
       MODULE MODI_XYTOLATLON
!      ###################################################
!
INTERFACE
!
SUBROUTINE SM_XYTOLATLON_A(PLAT0,PLON0,PRPK,PLATOR,PLONOR,  &
                           PXA,PYA,PLAT,PLON,KNX,KNY)
USE MODD_LUNIT
USE MODD_CST
!
REAL,                 INTENT(IN) :: PLATOR ! Latitude of the first point 
REAL,                 INTENT(IN) :: PLONOR ! Longitude of the first point
REAL, DIMENSION(KNX,KNY), INTENT(IN) :: PXA,PYA
REAL, DIMENSION(KNX,KNY), INTENT(OUT):: PLAT,PLON
REAL,                 INTENT(IN) :: PLAT0 ! Latitude of the origine point
REAL,                 INTENT(IN) :: PLON0 ! Longitude of the origine point
REAL,                 INTENT(IN) :: PRPK ! Longitude of the origine point
INTEGER,              INTENT(IN) :: KNX
INTEGER,              INTENT(IN) :: KNY
END SUBROUTINE SM_XYTOLATLON_A
!
END INTERFACE
!
END MODULE  MODI_XYTOLATLON
!
!      ###################################################
       SUBROUTINE SM_XYTOLATLON_A(PLAT0,PLON0,PRPK,PLATOR,PLONOR,  &
                                  PXA,PYA,PLAT,PLON,KNX,KNY)
!      ###################################################
!
!!****  *SM_LATLON_A * - Routine to compute geographical coordinates
!!
!!     PURPOSE
!!     -------
!        This routine computes the latitude and longitude of
!      an array given in cartesian conformal coordinates
!        Five map projections are available: 
!      - polar-stereographic from south pole  (XRPK=1),
!      - lambert conformal from south pole  (0<XRPK<1),
!      - mercator                             (XRPK=0),
!      - lambert conformal from north pole (-1<XRPK<0),
!      - polar-stereographic from north pole  (XRPK=-1).
!
!
!!**   METHOD
!!     ------
!!       Spherical earth approximation is used. Longitude origin is 
!!     set in Greenwich, and is positive eastwards. An anticlockwise 
!!     rotation of XBETA degrees is applied to the conformal frame 
!!     with respect to the geographical directions.
!!
!!       WARNING: ALL INPUT AND OUTPUT ANGLES ARE IN DEGREES...
!!
!!     EXTERNAL
!!     --------
!!       None
!!
!!     EXPLICIT ARGUMENTS
!!     ------------------
!!       PXHAT,PYHAT(:)  : 1D arrays of the "velocity" gridpoints
!!                         cartesian conformal coordinates (meters,input).
!!       PLATOR   : Latitude of the (1,1) point of the "mass" grid
!!                      (degrees,input);
!!       PLONOR   : Longitude of the (1,1) point of the "mass" grid
!!                      (degrees,input);
!!       PXHATM   : conformal coordinate x  (meters, mass-grid, input)
!!       PYHATM   : conformal coordinate y  (meters, mass-grid, input)
!!       PLAT    : latitude                (degrees, mass-grid, output)
!!       PLON    : longitude               (degrees, mass-grid, output)
!!
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST      : contains Physical constants
!!          XPI           : Pi;    
!!          XRADIUS       : Earth radius (meters);
!!
!!
!!     REFERENCE
!!     ---------
!!      Asencio N. et al., 1994, "Le projet de modele non-hydrostatique
!!            commun CNRM-LA, specifications techniques", 
!!            Note CNRM/GMME, 26, 139p, (Chapter 2).
!!      Ducrocq V., 1994, "Generation de la grille dans le modele",
!!            Note interne MNH, 5 mai, 3p.
!!      Joly A., 1992, "Geographic parameters for ARPEGE/ALADIN",
!!            Internal note ARPEGE/ALADIN, february 27,28p.
!!      Levallois J., 1970, "Geodesie generale", Tome 2, Collection
!!             de l'IGN, Eyrolles, Paris, 408p.
!!       
!!     AUTHOR
!!     ------
!!      P.M.       *LA*
!!
!!     MODIFICATION
!!     ------------
!!       Original  PM  24/05/94
!!       Updated   PM  27/07/94
!!       Updated   VD  23/08/94
!!       Updated   VM  24/10/95 projection from north pole (XRPK<0) and 
!!                              longitudes set between XLON0-180. and XLON0+180.
!!       Adapted  V.Bousquet 07/01/99 to xytolatlon for grib case
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*     0.     DECLARATIONS
!             ------------
!
USE MODD_CST
USE MODD_GRID
!
IMPLICIT NONE
!
!*     0.1    Declarations of arguments and results
!
INTEGER,               INTENT(IN) :: KNX
INTEGER,               INTENT(IN) :: KNY
!
REAL,                 INTENT(IN) :: PLATOR ! Latitude of the origine point
REAL,                 INTENT(IN) :: PLONOR ! Longitude of the origine point
REAL, DIMENSION(KNX,KNY), INTENT(IN) :: PXA,PYA   
				! given conformal coordinates of the 
				! processed points (meters);
REAL, DIMENSION(KNX,KNY), INTENT(OUT):: PLAT,PLON    
				! returned geographic latitudes and 
				! longitudes of the processed points 
				! (degrees).
REAL,                 INTENT(IN) :: PLAT0 ! Latitude of the origine point
REAL,                 INTENT(IN) :: PLON0 ! Longitude of the origine point
REAL,                 INTENT(IN) :: PRPK ! Longitude of the origine point
!
!*     0.2    Declarations of local variables
! 
REAL, DIMENSION(SIZE(PYA,1),SIZE(PYA,2)) :: ZYA
REAL :: ZRPK,ZLAT0,ZLON0,ZLATOR,ZLONOR
REAL :: ZRDSDG,ZCLAT0,ZSLAT0,ZCLATOR,ZSLATOR
REAL :: ZRO0,ZGA0 ,ZCGAM,ZSGAM
REAL :: ZXP,ZYP,ZEPSI,ZT1,ZRACLAT0
REAL :: ZBETA
!
REAL, DIMENSION(SIZE(PXA,1),SIZE(PXA,2)) :: ZATA,ZRO2,ZT2,ZXMI0,ZYMI0
!
INTEGER :: IRESP
!
!
!--------------------------------------------------------------------------------
!
!*     1.     Preliminary calculations for all projections
!             --------------------------------------------
!
ZRDSDG = XPI/180.         ! Degree to radian conversion factor
ZEPSI  = 10.*EPSILON(1.)      ! A small number
!
! By definition, (PLONOR,PLATOR) are the geographical 
! coordinates, and (ZXBM0,ZYBM0) the conformal cartesian 
! coordinates of the (1,1) point of the "mass" grid.
!
!-------------------------------------------------------------------------------
!
!*     2.     POLAR STEREOGRAPHIC AND LAMBERT CONFORMAL CASES
!             -----------------------------------------------
!                   (XRPK=1 P-stereo, 0<XRPK<1 Lambert)
!
IF(PRPK /= 0.) THEN
!
  IF (PRPK<0.) THEN     ! projection from north pole
    ZRPK=-PRPK
    ZBETA=-XBETA
    ZLAT0=-PLAT0
    ZLON0=PLON0+180.
    ZLATOR=-PLATOR
    ZLONOR=PLONOR+180.
    ZYA(:,:)=-PYA(:,:)
  ELSE                  ! projection from south pole
    ZRPK=PRPK
    ZBETA=XBETA
    ZLAT0=PLAT0
    ZLON0=PLON0
    ZLATOR=PLATOR
    ZLONOR=PLONOR
    ZYA(:,:)=PYA(:,:)
  ENDIF    
!
!*     2.1    Preliminary calculations
!
  ZCLAT0  = COS(ZRDSDG*ZLAT0)
  ZSLAT0  = SIN(ZRDSDG*ZLAT0)
  ZCLATOR = COS(ZRDSDG*ZLATOR)
  ZSLATOR = SIN(ZRDSDG*ZLATOR)
  ZRO0    = (XRADIUS/ZRPK)*(ABS(ZCLAT0))**(1.-ZRPK)     &
          * ((1.+ZSLAT0)*ABS(ZCLATOR)/(1.+ZSLATOR))**ZRPK
  ZGA0    = (ZRPK*(ZLONOR-ZLON0)-ZBETA)*ZRDSDG
  ZXP     = -ZRO0*SIN(ZGA0)
  ZYP     = ZRO0*COS(ZGA0)
!
!*    2.2    Longitude
!
  WHERE  (ABS(ZYA(:,:)-ZYP) < ZEPSI    &
     .AND.ABS(PXA(:,:)-ZXP) < ZEPSI)
    ZATA(:,:) = 0.
  ELSEWHERE
    ZATA(:,:) = ATAN2(-(ZXP-PXA(:,:)),(ZYP-ZYA(:,:)))/ZRDSDG
  END WHERE
  !
  PLON(:,:) = (ZBETA+ZATA(:,:))/ZRPK+ZLON0
!
!*   2.3     Latitude
!
  ZRO2(:,:) = (PXA(:,:)-ZXP)**2+(ZYA(:,:)-ZYP)**2
  ZT1       = (XRADIUS*(ABS(ZCLAT0))**(1.-ZRPK))**(2./ZRPK)   &
            * (1+ZSLAT0)**2
  ZT2(:,:)  = (ZRPK**2*ZRO2(:,:))**(1./ZRPK)
  !
  PLAT(:,:) = (XPI/2.-ACOS((ZT1-ZT2(:,:))/(ZT1+ZT2(:,:))))/ZRDSDG
!
  IF (PRPK<0.) THEN     ! projection from north pole
    PLAT(:,:)=-PLAT(:,:)
    PLON(:,:)=PLON(:,:)+180.
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*  3.        MERCATOR PROJECTION WITH ROTATION
!             ---------------------------------
!                       (XRPK=0)
!
ELSE
!
!*  3.1       Preliminary calculations
!
  ZCGAM    = COS(-ZRDSDG*XBETA)
  ZSGAM    = SIN(-ZRDSDG*XBETA)
  ZRACLAT0 = XRADIUS*COS(ZRDSDG*PLAT0)
!
!*  3.2       Longitude
!
  ZXMI0(:,:) = PXA(:,:)
  ZYMI0(:,:) = PYA(:,:)
  !
  PLON(:,:) = (ZXMI0(:,:)*ZCGAM+ZYMI0(:,:)*ZSGAM)  / (ZRACLAT0*ZRDSDG)+PLONOR
!
!*  3.3       Latitude
!
  ZT1       = ALOG(TAN(XPI/4.+PLATOR*ZRDSDG/2.))
  ZT2(:,:)  = (-ZXMI0*ZSGAM+ZYMI0*ZCGAM) / ZRACLAT0
  !
  PLAT(:,:) = (-XPI/2.+2.*ATAN(EXP(ZT1+ZT2(:,:))))/ZRDSDG
!
!---------------------------------------------------------------------------------
!
!*  4.        EXIT
!             ----
!
END IF
PLON(:,:)=PLON(:,:)+NINT((PLON0-PLON(:,:))/360.)*360.
RETURN
!---------------------------------------------------------------------------------
END SUBROUTINE SM_XYTOLATLON_A 
!---------------------------------------------------------------------------------
