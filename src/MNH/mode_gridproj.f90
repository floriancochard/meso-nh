!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODE_GRIDPROJ
!     ####################
!
!!****  *MODE_GRIDPROJ*  -   module routine  SM_GRIDPROJ
!!
!!      PURPOSE
!!      -------
!         This executable module packages a set of cartographic
!       module-procedures: 
!
!       SM_GRIDPROJ : to  compute the Jacobian in the case of 
!                    conformal projection;
!       SM_LATLON   : to compute geographic  from conformal
!                    cartesian coordinates;
!       SM_XYHAT    : to compute conformal cartesian from
!                     geographic coordinates;
!       LATREF2     : to compute the second reference latitude
!                    in the case of Lambert conformal projection
!
!!
!!**    IMPLICIT ARGUMENTS
!!      ------------------
!!           NONE
!!
!!      AUTHOR
!!      ------
!!          P.M.         *LA*
!!
!!      MODIFICATION
!!      ------------
!!          Original  24/05/94
!!                    05/02/15   M.Moge (LA-CNRS)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!    
!------------------------------------------------------------------------------
!
!*                0.  DECLARATIONS
USE MODE_MPPDB
USE MODD_CONF
!------------
!------------------------------------------------------------------------------
!
INTERFACE SM_LATLON
   MODULE PROCEDURE SM_LATLON_A,SM_LATLON_S
END INTERFACE
INTERFACE SM_XYHAT
   MODULE PROCEDURE SM_XYHAT_A,SM_XYHAT_S
END INTERFACE
!
CONTAINS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*                1.  ROUTINE  SM_GRIDPROJ   
!                     --------------------
!-------------------------------------------------------------------------------
!      ####################################################################
       SUBROUTINE SM_GRIDPROJ(PXHAT,PYHAT,PZHAT,PZS,                  &
                              OSLEVE,PLEN1,PLEN2,PZSMT,PLATOR,PLONOR, &
                              PMAP,PLAT,PLON,PDXHAT,PDYHAT,PZZ,PJ     )
!      ####################################################################
!
!!*****  *SM_GRIDPROJ * - Computes Jacobian J, map factor M,
!!    horizontal grid-meshes, latitude and longitude  at the 
!!    "mass" point locations.
!!
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute the Jacobian (J) at the
!     "mass" point location in the case of a conformal projection.
!     The map factor of the projection, the horizontal mesh-sizezs, and the 
!     the geograpical locations are also computed in the course of 
!     this calculation.
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
!!       The height, and the correction for spherical earth are first computed.
!!     Next, the conformal horizontal locations, the geographical coordinates 
!!     and the map factor  are derived at the "mass" grid-points.  
!!     The same formula can (hopefully) be used for all the projections cases
!!     (see Joly, 1992).
!!
!!       WARNING: ALL INPUT AND OUTPUT ANGLES ARE IN DEGREES...
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    EXPLICIT ARGUMENTS (not required, but given for convenience)
!!    ------------------
!!       PXHAT   : conformal coordinate x  (meters, u-grid, input)
!!       PYHAT   : conformal coordinate y  (meters, v-grid, input)
!!       PZHAT   : Gal-chen altitude  zhat (meters, w-grid, input)
!!       PZS     : topography              (meters, masss-grid, input)
!!       PLATOR  : Latitude of the origine point (degrees, mass grid, input)
!!       PLONOR  : Longitude of the origine point (degrees, mass grid, input)
!!       PMAP    : map scale               (no-unit, mass-grid, output)
!!       PLAT    : latitude                (degrees, mass-grid, output)
!!       PLON    : longitude               (degrees, mass-grid, output)
!!       PDXHAT  : local x mesh size       (meters, u-grid, output)
!!       PDYHAT  : local y mesh size       (meters, v-grid, output)
!!       PZZ     : true altitude z         (meters, w-grid, output)
!!       PJ      : jacobian                (no-unit, mass-grid, output)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_CONF      : contains declaration of configuration variables 
!!          LTHINSHELL   : Logical for thinshell approximation
!!          NVERB        : Listing verbosity control
!!
!!       Module MODD_CST       : contains Physical constants
!!          XPI          : Pi;    
!!          XRADIUS      : Earth radius (meters);
!!
!!       Module MODD_PARAMETERS: contains declaration of parameter variables 
!!          JPHEXT       : horizontal depth of arrays border 
!!          JPVEXT       : vertical   depth of arrays border
!!
!!       Module MODD_GRID      : contains spatial grid variables
!!          XLAT0   : map reference latitude  (degrees)
!!          XRPK    : projection parameter    (no-unit)
!!                            
!!
!!    REFERENCE
!!    ---------
!!      Asencio N. et al., 1994, "Le projet de modele non-hydrostatique
!!            commun CNRM-LA, specifications techniques", 
!!            Note CNRM/GMME, 26, 139p, (Chapter 2).
!!      Ducrocq V., 1994, "Generation de la grille dans le modele",
!!            Note interne MNH, 5 mai, 3p.
!!      Joly A., 1992, "Geographic parameters for ARPEGE/ALADIN",
!!            Internal note ARPEGE/ALADIN, february 27,28p.
!!      Levallois J., 1970, "Geodesie generale", Tome 2, Collection
!!             de l'IGN, Eyrolles, Paris, 408p.
!!      (chapters 2 and 3)
!!
!!
!!    AUTHOR
!!    ------
!!      P. Mascart        * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    PM  20/06/94  (from SM_GRIDCART by V. Ducrocq)
!!      Updated     PM  26/07/94
!!      Updated     VD  23/08/94
!!                      14/04/95  (Masson) bug in the ZYHTAM computation 
!!                      24/10/95  (Masson) controls during PMAP computation and
!!                                         projection from north pole (XPRK<0)
!!                      14/03/96  (Masson) enforce  -180<LONOR<+180     
!!                      01/11/96  (Mallet) bug for the MAP FACTOR computation
!!      Sleve coordinate        G. Zangler  *LA*             nov 2005
!!      MPPDB_CHECK     05/02/15   M.Moge (LA-CNRS)
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll, LIST1D_ll
USE MODD_LUNIT_n,     ONLY : TLUOUT
!
USE MODD_CONF          
USE MODD_CST          
USE MODD_PARAMETERS 
USE MODD_GRID      
!
USE MODI_VERT_COORD
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:),      INTENT(IN) :: PXHAT,PYHAT,PZHAT ! Positions x,y,z in 
                                                         ! the cartesian plane
REAL, DIMENSION(:,:),    INTENT(IN) :: PZS               ! Orography
LOGICAL,                INTENT(IN)  :: OSLEVE            ! flag for SLEVE coordinate
REAL,                   INTENT(IN)  :: PLEN1             ! Decay scale for smooth topography
REAL,                   INTENT(IN)  :: PLEN2             ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:),   INTENT(IN)  :: PZSMT             ! smooth orography
REAL,                    INTENT(IN) :: PLATOR            ! Latitude of the 
                                                         ! origine point 
REAL,                    INTENT(IN) :: PLONOR            ! Longitude of the 
                                                         ! origine point 
REAL, DIMENSION(:),     INTENT(OUT) :: PDXHAT            ! Local meshlength in 
                                                         ! x direction
REAL, DIMENSION(:),     INTENT(OUT) :: PDYHAT            ! Local meshlength in 
                                                         ! y direction 
REAL, DIMENSION(:,:),   INTENT(OUT) :: PMAP              ! Local map scale
                                                         ! of mass gridpoints
REAL, DIMENSION(:,:),   INTENT(OUT) :: PLAT,PLON         ! Latitude and longitude
                                                         ! of mass gridpoints
REAL, DIMENSION(:,:,:), INTENT(OUT):: PZZ                ! True altitude of the
                                                         ! w grid-point
REAL, DIMENSION(:,:,:), INTENT(OUT):: PJ                 ! Jacobian of the
                                                         ! GCS transformation
                                                         ! of mass gridpoints
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PXHAT,1),SIZE(PYHAT,1),SIZE(PZHAT,1)):: ZDZ ! Local z
                                                                 ! meshsize
REAL, DIMENSION(SIZE(PXHAT,1),SIZE(PYHAT,1),SIZE(PZHAT,1)):: ZAPZOA2 ! Spherical
                                                                     ! earth factor
                                                                     ! for J  
REAL, DIMENSION(SIZE(PXHAT,1),SIZE(PYHAT,1)):: ZXHATM   ! X and Y mass point
REAL, DIMENSION(SIZE(PXHAT,1),SIZE(PYHAT,1)):: ZYHATM   ! conformal coordinates
REAL ZRDSDG                              ! Radian to Degree conversion factor
REAL ZCLAT0,ZSLAT0                       ! Cos and Sin of XLAT0
REAL,DIMENSION(SIZE(PLAT,1),SIZE(PLAT,2)) :: ZLAT
REAL                                      :: ZRPK,ZLAT0
!
TYPE(LIST1D_ll), POINTER :: TZHALO1_ll ! pointer for the list of 1D and
TYPE(LIST_ll),   POINTER :: TZHALO_ll  ! 2D fields that have to be updated immediatly
INTEGER                  :: IINFO_ll   ! return code of the // routine
!
INTEGER      :: IIU,IJU,IKU      ! Uupper bounds of PXHAT,PYHAT,PZHAT  
INTEGER      :: IIE,IJE,IKE      ! End of usefull area of PXHAT,PYHAT,PZHAT  
INTEGER      :: IIB,IJB,IKB      ! Begining of usefull area of PXHAT,PYHAT,PZHAT  
INTEGER      :: IDELTA1          ! Switch=0 if thin shell approximation
INTEGER      :: ILUOUT,IRESP     ! Unit number for prints, FM error code 
INTEGER      :: JKLOOP           ! Index for control prints
!
INTEGER      :: JI,JJ,JK         ! loop index
!-------------------------------------------------------------------------------
!
!*       1.    RETRIEVE LOGICAL UNIT NUMBER FOR OUTPUT-LISTING AND DIMENSIONS 
!              --------------------------------------------------------------
!
ILUOUT = TLUOUT%NLU
!
IIU = UBOUND(PXHAT,1)     
IJU = UBOUND(PYHAT,1)    
IKU = UBOUND(PZHAT,1)       
IIE = IIU-JPHEXT
IJE = IJU-JPHEXT
IKE = IKU-JPVEXT
IIB = 1+JPHEXT
IJB = 1+JPHEXT
IKB = 1+JPVEXT
NULLIFY(TZHALO_ll)
NULLIFY(TZHALO1_ll)
!
IF(NVERB >= 10) THEN                               !Value control
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: IIU,IJU,IKU=',IIU,IJU,IKU
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: IIE,IJE,IKE=',IIE,IJE,IKE
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: IIB,IJB,IKB=',IIB,IJB,IKB
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTES Z   (W LEVEL)
!              ----------------------
!
CALL VERT_COORD(OSLEVE,PZS,PZSMT,PLEN1,PLEN2,PZHAT,PZZ)
!
IF(NVERB >= 10) THEN                               !Value control
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: Some PZS values:'
  WRITE(ILUOUT,*) PZS(1,1),PZS(IIU/2,IJU/2),PZS(IIU,IJU)  
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: Some PZZ values:'
  DO JKLOOP=1,IKU
    WRITE(ILUOUT,*) PZZ(1,1,JKLOOP),PZZ(IIU/2,IJU/2,JKLOOP), &
                    PZZ(IIU,IJU,JKLOOP)  
  END DO
  ! cancel MPPDB_CHECK if cprog=='SPAWN  '
  IF(CPROGRAM/='SPAWN ')&
  CALL MPPDB_CHECK3D(PZZ,"GRIDPROJ:PZZ",PRECISION)
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.   COMPUTE SPHERICAL EARTH FACTOR (MASS LEVEL)
!             --------------------------------------------
!
!     NOTE: In this routine LCARTESIAN is ALWAYS .F.
!           Hence, IDELTA2 is always set to 1
!
IF     (LTHINSHELL) IDELTA1=0                ! THIN SHELL APPROX.
IF(.NOT.LTHINSHELL) IDELTA1=1                ! NO THIN SHELL APPROX.
!
IF(NVERB >= 10) THEN                         !Value control
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: LTHINSHELL, IDELTA1=',LTHINSHELL,IDELTA1
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: XRADIUS=',XRADIUS
ENDIF
!
! For the time being, an inline implementation of  MZF
! is provided here.
!
ZAPZOA2(:,:,1:IKU-1) = (.5*((XRADIUS+IDELTA1*PZZ(:,:,1:IKU-1))    &
                     + (XRADIUS+IDELTA1*PZZ(:,:,2:IKU)))          &
                     /XRADIUS)**2
ZAPZOA2(:,:,IKU)     = 2.*ZAPZOA2(:,:,IKU-1)-ZAPZOA2(:,:,IKU-2)
!
IF(NVERB >= 10) THEN                               !Value control
  WRITE(ILUOUT,*) 'SM_GRIDPROJ: Some ZAPZOA2 values:'
  DO JKLOOP=1,IKU
    WRITE(ILUOUT,*) ZAPZOA2(1,1,JKLOOP),ZAPZOA2(IIU/2,IJU/2,JKLOOP), &
                    ZAPZOA2(IIU,IJU,JKLOOP)  
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.   COMPUTE ZXHAT AND ZYHAT AT MASS POINTS
!              -------------------------------------
!
ZXHATM(:,:) = 0.
ZYHATM(:,:) = 0.
ZXHATM(1:IIU-1,1) = .5*(PXHAT(1:IIU-1)+PXHAT(2:IIU))
ZXHATM(IIU,1)     = 2.*PXHAT(IIU)-ZXHATM(IIU-1,1)
ZXHATM(:,2:IJU)   = SPREAD(ZXHATM(:,1),2,IJU-1)
! cancel MPPDB_CHECK if cprog=='SPAWN  '
IF(CPROGRAM/='SPAWN ')&
CALL MPPDB_CHECK2D(ZXHATM,"GRIDPROJ:ZXHATM",PRECISION)
!
ZYHATM(1,1:IJU-1) = .5*(PYHAT(1:IJU-1)+PYHAT(2:IJU))
ZYHATM(1,IJU)     = 2.*PYHAT(IJU)-ZYHATM(1,IJU-1)
ZYHATM(2:IIU,:)   = SPREAD(ZYHATM(1,:),1,IIU-1)
! cancel MPPDB_CHECK if cprog=='SPAWN  '
IF(CPROGRAM/='SPAWN ')&
CALL MPPDB_CHECK2D(ZYHATM,"GRIDPROJ:ZYHATM",PRECISION)
! ZXHATM and ZXHATM have to be updated
CALL ADD2DFIELD_ll(TZHALO_ll,ZXHATM)
CALL ADD2DFIELD_ll(TZHALO_ll,ZYHATM)
CALL UPDATE_HALO_ll(TZHALO_ll,IINFO_ll)
CALL CLEANLIST_ll(TZHALO_ll)
!
!-------------------------------------------------------------------------------
!
!*       5.   COMPUTE LOCAL MESH-SIZES AT MASS POINTS
!              --------------------------------------
!
PDXHAT(1:IIU-1)  = PXHAT(2:IIU) - PXHAT(1:IIU-1)
PDXHAT(IIU)      = PDXHAT(IIU-1)
!
PDYHAT(1:IJU-1)  = PYHAT(2:IJU) - PYHAT(1:IJU-1)
PDYHAT(IJU)      = PDYHAT(IJU-1)
!
ZDZ(:,:,1:IKU-1) = PZZ(:,:,2:IKU) - PZZ(:,:,1:IKU-1)
ZDZ(:,:,IKU)     = ZDZ(:,:,IKU-1)
!  PDXHAT and PDYHAT have to be updated
CALL ADD1DFIELD_ll("X",TZHALO1_ll,PDXHAT)
CALL ADD1DFIELD_ll("Y",TZHALO1_ll,PDYHAT)
CALL UPDATE_1DHALO_ll(TZHALO1_ll,IINFO_ll)
DEALLOCATE(TZHALO1_ll)
CALL MPPDB_CHECK3D(ZDZ,"GRIDPROJ:ZDZ",PRECISION)
!
!-----------------------------------------------------------------------------
!
!*       5.   COMPUTE LATITUDES AND LONGITUDES AT MASS POINTS
!              -------------------------------------------------
!
CALL SM_LATLON(PLATOR,PLONOR,ZXHATM,ZYHATM,PLAT,PLON)
!
!-----------------------------------------------------------------------------
!
!*       6.  COMPUTE  MAP FACTOR AT MASS POINTS
!             -----------------------------------
!
  IF (XRPK<0.) THEN     ! projection from north pole 
    ZRPK=-XRPK
    ZLAT0=-XLAT0
    ZLAT(:,:)=-PLAT(:,:)
  ELSE                  ! projection from south pole
    ZRPK=XRPK
    ZLAT0=XLAT0
    ZLAT(:,:)=PLAT(:,:)
  ENDIF    
!
ZRDSDG = XPI/180.
ZCLAT0 = COS(ZRDSDG*ZLAT0)
ZSLAT0 = SIN(ZRDSDG*ZLAT0)
!
IF ((ABS(ZRPK-1.)>1.E-10).AND. (ANY(ABS(COS(ZRDSDG*ZLAT))<1.E-10))) THEN
  WRITE(ILUOUT,*) 'Error in SM_GRIDPROJ : '
  WRITE(ILUOUT,*) 'pole in the domain, but not with stereopolar projection'
!callabortstop
CALL ABORT
  STOP
ENDIF
!
IF (ABS(ZCLAT0)<1.E-10 .AND. (ABS(ZRPK-1.)<1.E-10)) THEN
  PMAP(:,:) = (1.+ZSLAT0)/(1.+SIN(ZRDSDG*ZLAT(:,:)))
ELSE
  WHERE (ABS(COS(ZRDSDG*ZLAT(:,:)))>1.E-10)
    PMAP(:,:) = ((ZCLAT0/COS(ZRDSDG*ZLAT(:,:)))**(1.-ZRPK))      &
              * ((1.+ZSLAT0)/(1.+SIN(ZRDSDG*ZLAT(:,:))))**ZRPK
  ELSEWHERE
    PMAP(:,:) = (1.+ZSLAT0)/(1.+SIN(ZRDSDG*ZLAT(:,:)))
  ENDWHERE
END IF
!
CALL MPPDB_CHECK2D(PMAP,"GRIDPROJ:PMAP",PRECISION)
!
IF(NVERB >= 10) THEN                               !Value control
  WRITE(ILUOUT,*) 'Some PMAP values:'
  WRITE(ILUOUT,*) PMAP(1,1),PMAP(IIU/2,IJU/2),PMAP(IIU,IJU)  
END IF
!
!-------------------------------------------------------------------------------

!
!*       8.    COMPUTE J AT MASS POINTS
!              -------------------------
!
DO JK=1,IKU ; DO JJ=1,IJU ; DO JI=1,IIU
  PJ(JI,JJ,JK)  = ZAPZOA2(JI,JJ,JK) * (1.0/PMAP(JI,JJ)**2)  &
                * PDXHAT(JI) * PDYHAT(JJ) * ZDZ(JI,JJ,JK) 
ENDDO ; ENDDO ; ENDDO
!
  CALL MPPDB_CHECK3D(PJ,"GRIDPROJ:PJ",PRECISION)
!
RETURN
!-----------------------------------------------------------------------------
END SUBROUTINE SM_GRIDPROJ
!-----------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------
!
!
!
!*              2.   ROUTINE SM_LATLON_S  (Scalar Version)
!                    -------------------
!----------------------------------------------------------------------------
!      #################################################
       SUBROUTINE SM_LATLON_S(PLATOR,PLONOR,PXHATM,PYHATM,PLAT,PLON)
!      #################################################
!
!!****  *SM_LATLON_S * - Routine to compute geographical coordinates
!!
!!     PURPOSE
!!     -------
!        This routine computes the latitude and longitude of
!      a single point from  the cartesian conformal coordinates
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
!!       PLATOR   : Latitude of the x=0, y=0 point (degrees,input);
!!       PLONOR   : Longitude of the x=0, y=0 point (degrees,input);
!!       PXHATM   : conformal coordinate x  (meters, mass-grid, input)
!!       PYHATM   : conformal coordinate y  (meters, mass-grid, input)
!!       PLAT     : latitude                (degrees, mass-grid, output)
!!       PLON     : longitude               (degrees, mass-grid, output)
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST        : contains Physical constants  
!!          XPI        : Pi;    
!!          XRADIUS    : Earth radius (meters);
!!
!!       Module MODD_GRID       : contains spatial grid variables
!!          XLON0,XLAT0  : Reference latitude and longitude for 
!!                            the conformal projection (degrees);
!!          XBETA        : Rotation angle of the conformal frame
!!                            with respect to the geographical  
!!                            north (degrees);
!!          XRPK         : Projection constant (0 Mercator,
!!                            0<XRPK<1 Lambert, 1 Polar-stereographic)
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
!!       Updated   VM     01/04 LONOR,LATOR refer to the x=0,y=0 point
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
REAL,               INTENT(IN) :: PLATOR ! Latitude of the origine point
REAL,               INTENT(IN) :: PLONOR ! Longitude of the origine point
REAL,               INTENT(IN) :: PXHATM,PYHATM ! given conformal coordinates of the 
                                                ! proccessed point (meters);
REAL,               INTENT(OUT):: PLAT,PLON ! returned geographic latitude and 
                                            ! longitude of the processed point 
                                            ! (degrees).
!
!*     0.2    Declarations of local variables
! 
REAL :: ZRPK,ZBETA,ZLAT0,ZLON0,ZLATOR,ZLONOR,ZYHATM
REAL :: ZRDSDG,ZCLAT0,ZSLAT0,ZCLATOR,ZSLATOR
REAL :: ZXBM0,ZYBM0,ZRO0,ZGA0 
REAL :: ZXP,ZYP,ZEPSI,ZT1,ZCGAM,ZSGAM,ZRACLAT0
!
REAL :: ZATA,ZRO2,ZT2,ZXMI0,ZYMI0
!
!--------------------------------------------------------------------------------
!
!*     1.     PRELIMINARY CALCULATIONS FOR ALL PROJECTIONS
!             --------------------------------------------
!
ZRDSDG = XPI/180.             ! Degree to radian conversion factor
ZEPSI  = 10.*EPSILON(1.)      ! A small number
!
! By definition, (PLONOR,PLATOR) are the geographical 
! coordinates, and (ZXBM0,ZYBM0) the conformal cartesian 
! coordinates x=0, y=0 of the grid.
!
!ZXBM0 = 0.5*(PXHAT(1)+PXHAT(2))
!ZYBM0 = 0.5*(PYHAT(1)+PYHAT(2))
ZXBM0 = 0.
ZYBM0 = 0.
!
!--------------------------------------------------------------------------------
!
!*     2.     POLAR STEREOGRAPHIC AND LAMBERT CONFORMAL CASES
!             -----------------------------------------------
!                   (XRPK=1 P-stereo, 0<XRPK<1 Lambert)
!
IF (XRPK /= 0.) THEN
!
  IF (XRPK<0.) THEN     ! projection from north pole
    ZRPK=-XRPK
    ZBETA=-XBETA
    ZLAT0=-XLAT0
    ZLON0=XLON0+180.
    ZLATOR=-PLATOR
    ZLONOR=PLONOR+180.
    ZYHATM=-PYHATM
    ZYBM0=-ZYBM0
  ELSE                  ! projection from south pole
    ZRPK=XRPK
    ZBETA=XBETA
    ZLAT0=XLAT0
    ZLON0=XLON0
    ZLATOR=PLATOR
    ZLONOR=PLONOR
    ZYHATM=PYHATM
  ENDIF    
!
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
  ZXP     = ZXBM0-ZRO0*SIN(ZGA0)
  ZYP     = ZYBM0+ZRO0*COS(ZGA0)
!
!*    2.2    Longitude
!
  IF(ABS(ZYHATM-ZYP) < ZEPSI.AND.ABS(PXHATM-ZXP) < ZEPSI)THEN
    ZATA = 0.
  ELSE
    ZATA = ATAN2(-(ZXP-PXHATM),(ZYP-ZYHATM))/ZRDSDG
  ENDIF
  !
  PLON = (ZBETA+ZATA)/ZRPK+ZLON0
!
!*   2.3     Latitude
!
  ZRO2 = (PXHATM-ZXP)**2+(ZYHATM-ZYP)**2
  ZT1  = (XRADIUS*(ABS(ZCLAT0))**(1.-ZRPK))**(2./ZRPK)   &
       * (1+ZSLAT0)**2
  ZT2  = (ZRPK**2*ZRO2)**(1./ZRPK)
  !
  PLAT = (XPI/2.-ACOS((ZT1-ZT2)/(ZT1+ZT2)))/ZRDSDG
!
  IF (XRPK<0.) THEN     ! projection from north pole 
    PLAT=-PLAT
    PLON=PLON+180.
  ENDIF    
!
!---------------------------------------------------------------------------------
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
  ZRACLAT0 = XRADIUS*COS(ZRDSDG*XLAT0)
!
!*  3.2       Longitude
!
  ZXMI0 = PXHATM-ZXBM0
  ZYMI0 = PYHATM-ZYBM0
  !
  PLON  = (ZXMI0*ZCGAM+ZYMI0*ZSGAM)/(ZRACLAT0*ZRDSDG)+PLONOR
!
!*  3.3       Latitude
!
  ZT1  = LOG(TAN(XPI/4.+PLATOR*ZRDSDG/2.))
  ZT2  = (-ZXMI0*ZSGAM+ZYMI0*ZCGAM)/ZRACLAT0
  !
  PLAT = (-XPI/2.+2.*ATAN(EXP(ZT1+ZT2)))/ZRDSDG
!
!---------------------------------------------------------------------------------
!
!*  4.        EXIT
!             ----
!
END IF
PLON=PLON+NINT((XLON0-PLON)/360.)*360.
RETURN
!--------------------------------------------------------------------------------
END SUBROUTINE SM_LATLON_S
!-------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------
!
!*              3.   ROUTINE SM_LATLON_A  (Array Version )
!                    -------------------
!--------------------------------------------------------------------------------
!      ###################################################
       SUBROUTINE SM_LATLON_A(PLATOR,PLONOR,  &
                              PXHATM,PYHATM,PLAT,PLON)
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
!!       PLATOR   : Latitude of the x=0, y=0 point (degrees,input);
!!       PLONOR   : Longitude of the x=0, y=0 point (degrees,input);
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
!!       Module MODD_GRID     : contains spatial grid variables
!!          XLON0,XLAT0   : Reference latitude and longitude for 
!!                          the conformal projection (degrees);
!!          XBETA         : Rotation angle of the conformal frame
!!                          with respect to the geographical  
!!                          north (degrees);
!!          XRPK          : Projection constant (0 Mercator,
!!                          0<XRPK<1 Lambert, 1 Polar-stereographic);
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
!!       Updated   VM     01/04 LONOR,LATOR refer to the x=0,y=0 point
!!       MPPDB_CHECK      05/02/15   M.Moge (LA-CNRS)
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
REAL,                 INTENT(IN) :: PLATOR ! Latitude of the origine point
REAL,                 INTENT(IN) :: PLONOR ! Longitude of the origine point
REAL, DIMENSION(:,:), INTENT(IN) :: PXHATM,PYHATM   
                                           ! given conformal coordinates of the 
                                           ! processed points (meters);
REAL, DIMENSION(:,:), INTENT(OUT):: PLAT,PLON    
                                           ! returned geographic latitudes and 
                                           ! longitudes of the processed points 
                                           ! (degrees).
!
!*     0.2    Declarations of local variables
! 
REAL, DIMENSION(SIZE(PYHATM,1),SIZE(PYHATM,2)) :: ZYHATM
REAL :: ZRPK,ZBETA,ZLAT0,ZLON0,ZLATOR,ZLONOR
REAL :: ZRDSDG,ZCLAT0,ZSLAT0,ZCLATOR,ZSLATOR
REAL :: ZXBM0,ZYBM0,ZRO0,ZGA0 
REAL :: ZXP,ZYP,ZEPSI,ZT1,ZCGAM,ZSGAM,ZRACLAT0
!
REAL, DIMENSION(SIZE(PXHATM,1),SIZE(PXHATM,2)) :: ZATA,ZRO2,ZT2,ZXMI0,ZYMI0
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
! coordinates x=0, y=0 of the grid.
!
ZXBM0 = 0.
ZYBM0 = 0.
!
!-------------------------------------------------------------------------------
!
!*     2.     POLAR STEREOGRAPHIC AND LAMBERT CONFORMAL CASES
!             -----------------------------------------------
!                   (XRPK=1 P-stereo, 0<XRPK<1 Lambert)
!
IF(XRPK /= 0.) THEN
!
  IF (XRPK<0.) THEN     ! projection from north pole
    ZRPK=-XRPK
    ZBETA=-XBETA
    ZLAT0=-XLAT0
    ZLON0=XLON0+180.
    ZLATOR=-PLATOR
    ZLONOR=PLONOR+180.
    ZYHATM(:,:)=-PYHATM(:,:)
    ZYBM0=-ZYBM0
  ELSE                  ! projection from south pole
    ZRPK=XRPK
    ZBETA=XBETA
    ZLAT0=XLAT0
    ZLON0=XLON0
    ZLATOR=PLATOR
    ZLONOR=PLONOR
    ZYHATM(:,:)=PYHATM(:,:)
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
  ZXP     = ZXBM0-ZRO0*SIN(ZGA0)
  ZYP     = ZYBM0+ZRO0*COS(ZGA0)
!
!*    2.2    Longitude
!
  WHERE  (ABS(ZYHATM(:,:)-ZYP) < ZEPSI    &
     .AND.ABS(PXHATM(:,:)-ZXP) < ZEPSI)
    ZATA(:,:) = 0.
  ELSEWHERE
    ZATA(:,:) = ATAN2(-(ZXP-PXHATM(:,:)),(ZYP-ZYHATM(:,:)))/ZRDSDG
  END WHERE
  !
! cancel MPPDB_CHECK if cprog=='SPAWN  '
  IF(CPROGRAM/='SPAWN ')&
  CALL MPPDB_CHECK2D(ZATA,"GRIDPROJ:ZATA",PRECISION)
!
  PLON(:,:) = (ZBETA+ZATA(:,:))/ZRPK+ZLON0
!
!*   2.3     Latitude
!
  ZRO2(:,:) = (PXHATM(:,:)-ZXP)**2+(ZYHATM(:,:)-ZYP)**2
! cancel MPPDB_CHECK if cprog=='SPAWN  '
  IF(CPROGRAM/='SPAWN ')&
  CALL MPPDB_CHECK2D(ZRO2,"GRIDPROJ:ZRO2",PRECISION)
!
  ZT1       = (XRADIUS*(ABS(ZCLAT0))**(1.-ZRPK))**(2./ZRPK)   &
            * (1+ZSLAT0)**2
  ZT2(:,:)  = (ZRPK**2*ZRO2(:,:))**(1./ZRPK)
  !
  PLAT(:,:) = (XPI/2.-ACOS((ZT1-ZT2(:,:))/(ZT1+ZT2(:,:))))/ZRDSDG
!
  IF (XRPK<0.) THEN     ! projection from north pole
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
  ZRACLAT0 = XRADIUS*COS(ZRDSDG*XLAT0)
!
!*  3.2       Longitude
!
  ZXMI0(:,:) = PXHATM(:,:)-ZXBM0
! cancel MPPDB_CHECK if cprog=='SPAWN  '
  IF(CPROGRAM/='SPAWN ')&
  CALL MPPDB_CHECK2D(ZXMI0,"GRIDPROJ:ZXMI0",PRECISION)
  !
  ZYMI0(:,:) = PYHATM(:,:)-ZYBM0
! cancel MPPDB_CHECK if cprog=='SPAWN  '
  IF(CPROGRAM/='SPAWN ')&
  CALL MPPDB_CHECK2D(ZYMI0,"GRIDPROJ:ZYMI0",PRECISION)
  !
  PLON(:,:) = (ZXMI0(:,:)*ZCGAM+ZYMI0(:,:)*ZSGAM)     &
            / (ZRACLAT0*ZRDSDG)+PLONOR
!
!*  3.3       Latitude
!
  ZT1       = ALOG(TAN(XPI/4.+PLATOR*ZRDSDG/2.))
  ZT2(:,:)  = (-ZXMI0(:,:)*ZSGAM+ZYMI0(:,:)*ZCGAM)/ZRACLAT0
  !
  PLAT(:,:) = (-XPI/2.+2.*ATAN(EXP(ZT1+ZT2(:,:))))/ZRDSDG
!
!---------------------------------------------------------------------------------
!
!*  4.        EXIT
!             ----
!
END IF
PLON(:,:)=PLON(:,:)+NINT((XLON0-PLON(:,:))/360.)*360.
RETURN
!---------------------------------------------------------------------------------
END SUBROUTINE SM_LATLON_A
!---------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------
!
!*              4.   ROUTINE SM_XYHAT_S (Scalar Version )
!                    ------------------
!--------------------------------------------------------------------------------
!      ##################################################
       SUBROUTINE SM_XYHAT_S(PLATOR,PLONOR,         &
                             PLAT,PLON,PXHATM,PYHATM)
!      ##################################################
!
!!****  *SM_XYHAT_S * - Routine to compute conformal coordinates
!!
!!     PURPOSE
!!     -------
!        This routine computes the cartesian conformal coordinates 
!      of a single point from  its latitude and longitude
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
!!       PLATOR   : Latitude of the x=0, y=0 point (degrees,input);
!!       PLONOR   : Longitude of the x=0, y=0 point (degrees,input);
!!       PXHATM   : conformal coordinate x  (meters, mass-grid, input)
!!       PYHATM   : conformal coordinate y  (meters, mass-grid, input)
!!       PLAT    : latitude                (degrees, mass-grid, output)
!!       PLON    : longitude               (degrees, mass-grid, output)
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST     : contains Physical constants
!!          XPI          : Pi;    
!!          XRADIUS      : Earth radius (meters);
!!
!!       Module MODD_GRID    : contains spatial grid variables
!!          XLON0,XLAT0  : Reference latitude and longitude for 
!!                         the conformal projection (degrees);
!!          XBETA        : Rotation angle of the conformal frame
!!                         with respect to the geographical  
!!                         north (degrees);
!!          XRPK         : Projection constant (0 Mercator,
!!                         0<XRPK<1 Lambert, 1 Polar-stereographic);
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
!!       Updated   VM     01/04 LONOR,LATOR refer to the x=0,y=0 point
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
REAL,               INTENT(IN) :: PLATOR ! Latitude of the origine point
REAL,               INTENT(IN) :: PLONOR ! Longitude of the origine point
REAL,               INTENT(IN) :: PLAT,PLON 
                                         ! given geographic latitude and 
                                         ! longitude of the processed point 
                                         ! (degrees).
REAL,               INTENT(OUT):: PXHATM,PYHATM 
                                         ! returned conformal coordinates of 
                                         ! the processed point (meters);
!
!*     0.2    Declarations of local variables
! 
REAL :: ZRPK,ZBETA,ZLAT0,ZLON0,ZLATOR,ZLONOR
REAL :: ZLAT,ZLON
REAL :: ZRDSDG,ZCLAT0,ZSLAT0,ZCLATOR,ZSLATOR
REAL :: ZXBM0,ZYBM0,ZRO0,ZGA0 
REAL :: ZXP,ZYP,ZCGAM,ZSGAM,ZRACLAT0,ZXE,ZYE
!
REAL :: ZCLAT,ZSLAT,ZRO,ZGA,ZXPR,ZYPR
!
!--------------------------------------------------------------------------------
!
!*     1.     PRELIMINARY CALCULATION FOR ALL PROJECTIONS
!             -------------------------------------------
!
ZRDSDG = XPI/180.         ! Degree to radian conversion factor
!
! By definition, (PLONOR,PLATOR) are the geographical 
! coordinates of the x=0, y=0 point.
!
ZXBM0 = 0.
ZYBM0 = 0.
!
ZLON=PLON
ZLON=ZLON+NINT((XLON0-ZLON)/360.)*360.
!
ZLONOR=PLONOR
ZLONOR=ZLONOR+NINT((XLON0-ZLONOR)/360.)*360.
!---------------------------------------------------------------------------------
!
!*     2.     POLAR STEREOGRAPHIC AND LAMBERT CONFORMAL CASES
!             -----------------------------------------------
!                   (XRPK=1 P-stereo, 0<XRPK<1 Lambert)
!
IF(XRPK /= 0.) THEN
!
  IF (XRPK<0.) THEN     ! projection from north pole
    ZRPK=-XRPK
    ZBETA=-XBETA
    ZLAT0=-XLAT0
    ZLON0=XLON0+180.
    ZLATOR=-PLATOR
    ZLONOR=ZLONOR+180.
    ZLAT=-PLAT
    ZLON=ZLON+180.
    ZYBM0=-ZYBM0
  ELSE                  ! projection from south pole
    ZRPK=XRPK
    ZBETA=XBETA
    ZLAT0=XLAT0
    ZLON0=XLON0
    ZLATOR=PLATOR
    ZLONOR=ZLONOR
    ZLAT=PLAT
    ZLON=ZLON
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
  ZXP     = ZXBM0-ZRO0*SIN(ZGA0)
  ZYP     = ZYBM0+ZRO0*COS(ZGA0)
!
!*    2.2    Conformal coordinates in meters
!
  ZCLAT  = COS(ZRDSDG*ZLAT)
  ZSLAT  = SIN(ZRDSDG*ZLAT)
  ZRO    = (XRADIUS/ZRPK)*(ABS(ZCLAT0))**(1.-ZRPK)    &
         * ((1.+ZSLAT0)*ABS(ZCLAT)/(1.+ZSLAT))**ZRPK
  ZGA    = (ZRPK*(ZLON-ZLON0)-ZBETA)*ZRDSDG
!
  PXHATM = ZXP+ZRO*SIN(ZGA)
  PYHATM = ZYP-ZRO*COS(ZGA)
!
  IF (XRPK<0.) THEN     ! projection from north pole
    PYHATM=-PYHATM
  ENDIF
!
!
!------------------------------------------------------------------------------
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
  ZRACLAT0 = XRADIUS*COS(ZRDSDG*XLAT0)
  ZXE      = ZXBM0*ZCGAM+ZYBM0*ZSGAM            &
             - ZRACLAT0*(PLONOR-XLON0)*ZRDSDG  
  ZYE      =-ZXBM0*ZSGAM+ZYBM0*ZCGAM            &
            - ZRACLAT0*LOG(TAN(XPI/4.+PLATOR*ZRDSDG/2.))
!
!*  3.2       Conformal coordinates
!
  ZXPR   = ZRACLAT0*(ZLON-XLON0)*ZRDSDG+ZXE
  ZYPR   = ZRACLAT0*LOG(TAN(XPI/4.+PLAT*ZRDSDG/2.))+ZYE
  !
  PXHATM = ZXPR*ZCGAM-ZYPR*ZSGAM
  PYHATM = ZXPR*ZSGAM+ZYPR*ZCGAM
!
!-------------------------------------------------------------------------------
!
!*  4.        EXIT
!             ----
!
END IF
RETURN
!-------------------------------------------------------------------------------
END SUBROUTINE SM_XYHAT_S
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!*              5.   ROUTINE SM_XYHAT_A (Array Version )
!                    ------------------
!-------------------------------------------------------------------------------
!      ################################################
       SUBROUTINE SM_XYHAT_A(PLATOR,PLONOR,         &
                             PLAT,PLON,PXHATM,PYHATM)
!      ################################################
!
!!****  *SM_XYHAT_A * - Routine to compute conformal coordinates
!!
!!
!!     PURPOSE
!!     -------
!        This routine computes the cartesian conformal coordinates 
!      of an array given in latitude-longitude coordinates
!        Three map projections are available: 
!      - polar-stereographic (XRPK=1),
!      - lambert conformal  (0<XRPK<1),
!      - mercator (XRPK=0).
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
!!       PLATOR   : Latitude of the x=0, y=0 point (degrees,input);
!!       PLONOR   : Longitude of the x=0, y=0 point (degrees,input);
!!       PXHATM   : conformal coordinate x  (meters, mass-grid, input)
!!       PYHATM   : conformal coordinate y  (meters, mass-grid, input)
!!       PLAT    : latitude                (degrees, mass-grid, output)
!!       PLON    : longitude               (degrees, mass-grid, output)
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST         : contains Physical constants
!!          XPI         : Pi;    
!!          XRADIUS     : Earth radius (meters);
!!
!!       Module MODD_GRID        : contains spatial grid variables
!!          XLON0,XLAT0 : Reference latitude and longitude for 
!!                        the conformal projection (degrees);
!!          XBETA       : Rotation angle of the conformal frame
!!                        with respect to the geographical  
!!                        north (degrees);
!!          XRPK        : Projection constant (0 Mercator,
!!                        0<XRPK<1 Lambert, 1 Polar-stereographic);
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
!!       Original PM  24/05/94
!!       Updated  PM  27/07/94
!!       Updated  VD  23/08/94
!!       Updated  VM  24/10/95 projection from north pole (XRPK<0) and 
!!                             longitudes set between XLON0-180. and XLON0+180.
!!       Updated  VM     01/04 LONOR,LATOR refer to the x=0,y=0 point
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
REAL,                INTENT(IN) :: PLATOR ! Latitude of the origine point
REAL,                INTENT(IN) :: PLONOR ! Longitude of the origine point
REAL, DIMENSION(:,:), INTENT(IN):: PLAT,PLON     
          ! given geographic latitude and 
          ! longitude of the processed  
          ! array (degrees).
REAL, DIMENSION(:,:), INTENT(OUT):: PXHATM,PYHATM  
          ! returned conformal coordinates of 
          ! the processed array (meters);
!
!*     0.2    Declarations of local variables
! 
REAL,DIMENSION(SIZE(PLAT,1),SIZE(PLAT,2)) :: ZLAT,ZLON
REAL :: ZRPK,ZBETA,ZLAT0,ZLON0,ZLATOR,ZLONOR
REAL :: ZRDSDG,ZCLAT0,ZSLAT0,ZCLATOR,ZSLATOR
REAL :: ZXBM0,ZYBM0,ZRO0,ZGA0 
REAL :: ZXP,ZYP,ZCGAM,ZSGAM,ZRACLAT0,ZXE,ZYE
!
REAL,DIMENSION(SIZE(PLAT,1),SIZE(PLAT,2)) :: ZCLAT,ZSLAT,ZRO,ZGA,ZXPR,ZYPR
!
!
!-------------------------------------------------------------------------------
!
!*     1.     PRELIMINARY CALCULATION FOR ALL PROJECTIONS
!             -------------------------------------------
!
ZRDSDG = XPI/180.         ! Degree to radian conversion factor
!
! By definition, (PLONOR,PLATOR) are the geographical 
! coordinates of the x=0, y=0 point.
!
ZXBM0 = 0.
ZYBM0 = 0.
!
ZLON(:,:)=PLON(:,:)
ZLON(:,:)=ZLON(:,:)+NINT((XLON0-ZLON(:,:))/360.)*360.
!
ZLONOR=PLONOR
ZLONOR=ZLONOR+NINT((XLON0-ZLONOR)/360.)*360.
!------------------------------------------------------------------------------
!
!*     2.     POLAR SEREOGRAPHIC AND LAMBERT CONFORMAL CASES
!             ----------------------------------------------
!                   (XRPK=1 P-stereo, 0<XRPK<1 Lambert)
!
IF(XRPK  /=  0.) THEN
!
  IF (XRPK<0.) THEN     ! projection from north pole
    ZRPK=-XRPK
    ZBETA=-XBETA
    ZLAT0=-XLAT0
    ZLON0=XLON0+180.
    ZLATOR=-PLATOR
    ZLONOR=ZLONOR+180.
    ZLAT(:,:)=-PLAT(:,:)
    ZLON(:,:)=ZLON(:,:)+180.
    ZYBM0=-ZYBM0
  ELSE                  ! projection from south pole
    ZRPK=XRPK
    ZBETA=XBETA
    ZLAT0=XLAT0
    ZLON0=XLON0
    ZLATOR=PLATOR
    ZLONOR=ZLONOR
    ZLAT(:,:)=PLAT(:,:)
    ZLON(:,:)=ZLON(:,:)
  ENDIF    
!
!*     2.1    Preliminary calculations
!
  ZCLAT0  = COS(ZRDSDG*ZLAT0)
  ZSLAT0  = SIN(ZRDSDG*ZLAT0)
  ZCLATOR = COS(ZRDSDG*ZLATOR)
  ZSLATOR = SIN(ZRDSDG*ZLATOR)
  ZRO0    = (XRADIUS/ZRPK)*(ABS(ZCLAT0))**(1.-ZRPK)    &
          * ((1.+ZSLAT0)*ABS(ZCLATOR)/(1.+ZSLATOR))**ZRPK
  ZGA0    = (ZRPK*(ZLONOR-ZLON0)-ZBETA)*ZRDSDG
  ZXP     = ZXBM0-ZRO0*SIN(ZGA0)
  ZYP     = ZYBM0+ZRO0*COS(ZGA0)
!
!*    2.2    Conformal coordinates in meters
!
  ZCLAT(:,:)  = COS(ZRDSDG*ZLAT(:,:))
  ZSLAT(:,:)  = SIN(ZRDSDG*ZLAT(:,:))
  ZRO(:,:)    = (XRADIUS/ZRPK)*(ABS(ZCLAT0))**(1.-ZRPK)    &
               * ((1.+ZSLAT0)*ABS(ZCLAT(:,:))/(1.+ZSLAT(:,:)))**ZRPK
  ZGA(:,:)    = (ZRPK*(ZLON(:,:)-ZLON0)-ZBETA)*ZRDSDG
!
  PXHATM(:,:) = ZXP+ZRO(:,:)*SIN(ZGA(:,:))
  PYHATM(:,:) = ZYP-ZRO(:,:)*COS(ZGA(:,:))
!
  IF (XRPK<0.) THEN     ! projection from north pole
    PYHATM(:,:)=-PYHATM(:,:)
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
  ZRACLAT0 = XRADIUS*COS(ZRDSDG*XLAT0)
  ZXE      = ZXBM0*ZCGAM+ZYBM0*ZSGAM            &
            - ZRACLAT0*(PLONOR-XLON0)*ZRDSDG  
  ZYE      =-ZXBM0*ZSGAM+ZYBM0*ZCGAM            &
            - ZRACLAT0*LOG(TAN(XPI/4.+PLATOR*ZRDSDG/2.))
!
!*  3.2       Conformal coordinates
!
  ZXPR(:,:)   = ZRACLAT0*(ZLON(:,:)-XLON0)*ZRDSDG+ZXE
  ZYPR(:,:)   = ZRACLAT0*LOG(TAN(XPI/4.+PLAT(:,:)*ZRDSDG/2.))+ZYE
  !
  PXHATM(:,:) = ZXPR(:,:)*ZCGAM-ZYPR(:,:)*ZSGAM
  PYHATM(:,:) = ZXPR(:,:)*ZSGAM+ZYPR(:,:)*ZCGAM
!
!-------------------------------------------------------------------------------
!
!*  4.        EXIT
!             ----
!
END IF
RETURN
!-------------------------------------------------------------------------------
END SUBROUTINE SM_XYHAT_A
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!
!*              6.   FUNCTION LATREF2
!                    -----------------
!-------------------------------------------------------------------------------
!      #############################
       FUNCTION LATREF2(PLAT0,PRPK)
!      #############################
!
!!****  *LATREF2 * - returns the Lambert second reference latitude
!!
!!     PURPOSE
!!     -------
!        This routine computes the second reference latitude 
!      of a Lambert conformal projection for given projection
!      parameter PRPK and primary reference latitude PLAT0.
!        This second latitude is used in US and UK to define
!      the secant Lambert projection (as a substitute for the
!      cone constant PRPK used in France by IGN).
!        This latitude is required to call the NCAR map projection
!      package with the Lambert option.
!
!!**   METHOD
!!     ------
!!       The so-called "constant of the cone" equation is solved 
!!     using a simple Newton-Raphson iteration. The spherical earth 
!!     approximation is used.   
!!
!!       WARNING: ALL INPUT AND OUTPUT ANGLES ARE IN DEGREES...
!!
!!     EXTERNAL
!!     --------
!!       None
!!
!!     EXPLICIT ARGUMENTS
!!     -------------------
!!       PRPK    : projection factor       (no-unit, input)  
!!       PLAT0   : map reference latitude  (degrees, input)             
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!       Module MODD_CST      : contains Physical constants
!!          XPI        : Pi;    
!!
!!       Module MODD_LUNIT    : contains logical unit names
!!          CLUOUT0    : Output listing file name
!!
!!     REFERENCE
!!     ---------
!!      Joly A., 1992, "Geographic parameters for ARPEGE/ALADIN",
!!            Internal note ARPEGE/ALADIN, february 27,28p.
!!      Levallois J., 1970, "Geodesie generale", Tome 2, Collection
!!             de l'IGN, Eyrolles, Paris, 408p.
!!      Pearson F. II, 1990,"Map projections: theory and applications",
!!             CRC Press, Boca Raton, Florida, 372p. (Chapter 5).
!!       
!!     AUTHOR
!!     ------
!!      P.M.       *LA*
!!
!!     MODIFICATION
!!     ------------
!!       Original PM  24/05/94
!!       Updated  PM  27/07/94
!!       Updated  VD  25/08/94
!!       Updated  VM  24/10/95 projection from north pole (XRPK<0)
!!       Updated  VM  08/10/96 output-listing choice
!!       Updated  IM  27/11/03 special case if projection plane is tangent
!!
!-------------------------------------------------------------------------------
!
!*     0.     DECLARATIONS
!             ------------
!
USE MODD_CST
USE MODD_LUNIT_n
!
IMPLICIT NONE
!
!*     0.1    Declarations of arguments and results
!
REAL,INTENT(IN):: PLAT0,PRPK    ! Given first standard latitude (degrees)
                                ! and projection parameter (cone 
                                ! constant) for the Lambert conformal
                                ! projection used.
REAL :: LATREF2                 ! Returned latitude of the second 
                                ! reference (or standard) parallel
                                ! of the projection.
!
!*     0.2    Declarations of local variables
!
REAL    :: ZRPK
REAL    :: ZRDSDG,ZEPSI,ZLAT0,ZLAT,ZDLAT,ZGLAT,ZGPRSG
INTEGER :: ITER,ITERMAX
INTEGER :: ILUOUT,IRESP  
!
!-------------------------------------------------------------------------------
!
!*     1.     PRELIMINARY CALCULATIONS
!             ------------------------
!
ZRDSDG  = XPI/180.         ! Degree to radian conversion factor
ZEPSI   = 10.*EPSILON(1.)   ! a small number
ITERMAX = 10              ! number of iteration allowed
!
IF (PRPK ==SIN(PLAT0*ZRDSDG)) THEN     ! projection plane tangent to the sphere
  LATREF2 = PLAT0
ELSE                                   !       "          intersect the sphere
  ZLAT0   = PLAT0*ZRDSDG      ! Switch to radians
!
  ZLAT    = XPI-4.*ATAN(SQRT((1.-PRPK)/(1.+PRPK)))-ZLAT0    
  ITER    = 0                  ! Choose the side of the nice root
  ZDLAT   = 0.                ! and sets up for the loop
!
!
  IF (PRPK<0.) THEN    ! projection from north pole
    ZRPK=-PRPK
    ZLAT0=-ZLAT0
    ZLAT=-ZLAT
  ELSE                 ! projection from south pole
    ZRPK=PRPK
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*     2.     NEWTON-RAPHSON LOOP
!             -------------------
  DO
    ITER   = ITER+1
    ZLAT   = ZLAT+ZDLAT
    ZGLAT  =(COS(ZLAT)/COS(ZLAT0))*                         & 
            (((1.+SIN(ZLAT))/(1.+SIN(ZLAT0)))**(ZRPK/(1.-ZRPK)))
    ZGPRSG = ((ZRPK/(1.-ZRPK))*(COS(ZLAT)/(1.+SIN(ZLAT)))    &
          - (SIN(ZLAT)/COS(ZLAT)))*ZGLAT
    ZDLAT  = (1.-ZGLAT)/ZGPRSG
    !
    IF((ABS(ZGLAT-1.) <= ZEPSI).OR.(ITER >= ITERMAX))   EXIT
  END DO
!
  IF (PRPK<0.) ZLAT=-ZLAT
  LATREF2  = ZLAT/ZRDSDG     ! Degrees restored
!
ENDIF
!-------------------------------------------------------------------------------
!
!*  3.        EXIT
!             ----
!
IF(ITER <= ITERMAX)  RETURN
!
ILUOUT = TLUOUT%NLU
!
WRITE(ILUOUT,*) ' Error in function LATREF2 (module MODE_GRIDPROJ)'
WRITE(ILUOUT,*) ' Function fails to converge after ',ITER,' iterations.'
WRITE(ILUOUT,*) ' LATREF2=',LATREF2,' Residual=',ZGLAT-1.,           &
                ' ZEPSI=',ZEPSI,' Last increment=',ZDLAT/ZRDSDG
WRITE(ILUOUT,*) ' JOB ABORTS...'
!callabortstop
CALL ABORT
STOP
!-------------------------------------------------------------------------------
END FUNCTION LATREF2
!-------------------------------------------------------------------------------
!
END MODULE MODE_GRIDPROJ
