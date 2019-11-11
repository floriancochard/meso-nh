!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #############################
      MODULE MODI_CHANGE_GRIBEX_VAR
!     #############################
INTERFACE
      SUBROUTINE CHANGE_GRIBEX_VAR(PA_LS,PB_LS,PP00_LS,PPS_LS,PZS_LS, &
                          PT_LS,PQ_LS,PPMASS_LS,PZFLUX_LS,PZMASS_LS,  &
                          PTHV_LS,PR_LS,PHU_LS,PU_LS,PV_LS,PW_LS      )
!
REAL,DIMENSION(:),     INTENT(IN) :: PA_LS    ! function A in definition of eta
REAL,DIMENSION(:),     INTENT(IN) :: PB_LS    ! function B in definition of eta
REAL,                  INTENT(IN) :: PP00_LS  ! reference pressure in eta
REAL,DIMENSION(:,:),   INTENT(IN) :: PPS_LS   ! surface pressure 
REAL,DIMENSION(:,:),   INTENT(IN) :: PZS_LS   ! large scale orography 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PT_LS    ! temperature 
REAL,DIMENSION(:,:,:,:),INTENT(IN):: PQ_LS    ! specific ratio of hydrometeors
!
REAL,DIMENSION(:,:,:), INTENT(OUT):: PPMASS_LS! pressure
REAL,DIMENSION(:,:,:), INTENT(OUT):: PZFLUX_LS  ! altitude of the pressure levels
REAL,DIMENSION(:,:,:), INTENT(OUT):: PZMASS_LS  ! altitude of the mass levels
REAL,DIMENSION(:,:,:), INTENT(OUT):: PTHV_LS    ! thetav
REAL,DIMENSION(:,:,:,:), INTENT(OUT):: PR_LS      ! water mixing ratios
REAL,DIMENSION(:,:,:), INTENT(OUT):: PHU_LS     ! relative humidity
!
REAL,DIMENSION(:,:,:), INTENT(IN),OPTIONAL :: PU_LS      ! pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(IN),OPTIONAL :: PV_LS      ! pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(OUT),OPTIONAL:: PW_LS      ! vertical wind component
!
END SUBROUTINE CHANGE_GRIBEX_VAR
END INTERFACE
END MODULE MODI_CHANGE_GRIBEX_VAR
!     #########################################################################
      SUBROUTINE CHANGE_GRIBEX_VAR(PA_LS,PB_LS,PP00_LS,PPS_LS,PZS_LS, &
                          PT_LS,PQ_LS,PPMASS_LS,PZFLUX_LS,PZMASS_LS,  &
                          PTHV_LS,PR_LS,PHU_LS,PU_LS,PV_LS,PW_LS      )
!     #########################################################################
!
!!****  *CHANGE_GRIBEX_VAR* - changes the variables and computes the
!!                                       grib grid altitude
!!
!!    PURPOSE
!!    -------
!!    This routine computes:
!!  1 the pressure and Exner function on the grib grid
!!  2 the Exner function on the mass points of the grib grid
!!  3 the vapor mixing ratio and the virtual potential temperature at the mass
!!    of the grib grid
!!
!!  4 the altitude of the pressure points of the grib grid
!!  5 the altitude of the mass points of the grib grid
!!  6 the projection from grib contravariant wind components to covariant
!!     wind components
!!
!!**  METHOD
!!    ------
!!
!!  1 p= A*p00 + B*ps
!!    PI=(p/p00)**(Rd/Cpd)
!!
!!  2 PI is computed on mass points
!!
!!      ~        PI(l)-PI(l+1)
!!     PI(l) = -----------------
!!             lnPI(l)-lnPI(l+1)
!!
!!      ~                              ~
!!     Pi at l=1 is computed supposing p(l=1) = p(l=1)/2
!!     (It is intrinsiquely supposed that the grib file top eta level
!!      (l=0) is at p=0)
!!
!!  3 rv=1/(1/q-1)
!!              ~
!!    thetav=T/PI*(1+Rv/Rd*rv)/(1+rw)
!!
!!  4 the altitude of the pressure points PZFLUX_LS is computed by integration of the
!!    hydrostatic relation from bottom (l=ILU) to top (l=1)
!!
!!                   Cpd thetav(l)
!!    z  - z    = - -------------- (PI   - PI )
!!     l    l-1           g           l+1    l
!!
!!  5 the altitude of the mass points PZMASS_LS is computed as:
!!
!!    ~          Cpd (0.75thetav(l)+0.25thetav(l-1)     ~
!!    z  = z  - ----------------------------------- * ( PI(l) - PI(l) )
!!     l    l                  g
!!
!!
!!  6 the true vertical wind speed is integrated according to formula:
!!    (with some necessary approximations)
!!
!!                            /l                    ->
!!    rho(l) g w(l) =        | ( DIV | (A'p00+B'ps) U ) d eta
!!                          /top     |eta
!!
!!                            /surf                 ->
!!                   -B(l) * | ( DIV | (A'p00+B'ps) U ) d eta
!!                          /top     |eta
!!                                -> ->
!!                  + rho(surf) g U.grad(zs)
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    Module MODI_SHUMAN     : interface for Shuman operators
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models.
!!         NVERB      : verbosity level for output-listing
!!      Module MODD_CONF1
!!         NRR
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_CST       : contains physical constants
!!         XRD : gas constant for dry air
!!         XRV : gas constant for vapor
!!         XP00: reference pressure
!!         XCPD: specific heat for dry air
!!         XG  : gravity constant
!!         XRADIUS : earth radius
!!      Module MODD_REF       : contains anelastic reference state variables.
!!         XEXNTOP
!!      Module MODD_PARAMETERS
!!         JPHEXT
!!         JPVEXT
!!      Module MODD_GRID1
!!         XXHAT
!!         XYHAT
!!         XZHAT
!!         XMAP
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
!!      Original    31/01/96
!!      Modification
!!         Masson   25/05/96 set the correct high for the upper grib level
!!         Masson   20/08/96 bug in thetav definition
!!         Masson   20/10/96 move deallocations in calling routine
!!         Masson   06/12/96 add air temperature at ground
!!         Masson   12/12/96 add vertical wind component
!!         Masson   12/06/97 add relative humidity
!!         J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!         Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!         Pergaud  : 2018 add GFS
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST
USE MODD_GRID_n
USE MODD_LUNIT, ONLY: CLUOUT0, TLUOUT0
USE MODD_PARAMETERS
USE MODD_REF
!
USE MODE_THERMO
!
USE MODI_SHUMAN
USE MODI_WATER_SUM
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL,DIMENSION(:),     INTENT(IN) :: PA_LS    ! function A in definition of eta
REAL,DIMENSION(:),     INTENT(IN) :: PB_LS    ! function B in definition of eta
REAL,                  INTENT(IN) :: PP00_LS  ! reference pressure in eta
REAL,DIMENSION(:,:),   INTENT(IN) :: PPS_LS   ! surface pressure 
REAL,DIMENSION(:,:),   INTENT(IN) :: PZS_LS   ! large scale orography 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PT_LS    ! temperature 
REAL,DIMENSION(:,:,:,:),INTENT(IN):: PQ_LS    ! specific ratio of hydrometeors
!
REAL,DIMENSION(:,:,:), INTENT(OUT):: PPMASS_LS! pressure
REAL,DIMENSION(:,:,:), INTENT(OUT):: PZFLUX_LS  ! altitude of the pressure levels
REAL,DIMENSION(:,:,:), INTENT(OUT):: PZMASS_LS  ! altitude of the mass levels
REAL,DIMENSION(:,:,:), INTENT(OUT):: PTHV_LS    ! thetav
REAL,DIMENSION(:,:,:,:), INTENT(OUT):: PR_LS      ! water mixing ratios
REAL,DIMENSION(:,:,:), INTENT(OUT):: PHU_LS     ! relative humidity
!
REAL,DIMENSION(:,:,:), INTENT(IN),OPTIONAL :: PU_LS      ! pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(IN),OPTIONAL :: PV_LS      ! pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(OUT),OPTIONAL:: PW_LS      ! vertical wind component
!
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER                                                  :: IIU,IJU,ILU
INTEGER                                                  :: JL,JRR
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZPFLUX_LS
!                        ! pressure on the grib grid.
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZEXNFLUX_LS,     &
                                                            ZEXNMASS_LS
!                        ! local Exner function at pressure and mass points
!                        ! on the grib grid
REAL,DIMENSION(:), ALLOCATABLE                           :: ZA_PRIME_LS,     &
                                                            ZB_PRIME_LS
!                        ! derivative of A,B function
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZDPDETA_LS
!                        ! term dp/deta = A'p00+B'ps
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZSURFCOR_LS
!                        ! surface correction for w computation
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZHDIV_LS
!                        ! horizontal divergence for w computation
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZINTDIV_LS
!                        ! integration of horizontal divergence
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZRHO_LS
!                        ! mass density at pressure points
REAL,DIMENSION(:,:,:), ALLOCATABLE                       :: ZRHOMASS_LS
!                        ! mass density at mass points
!
INTEGER             :: IIB,IIE,IJB,IJE ! interior domaine bound
INTEGER             :: JI
!-------------------------------------------------------------------------------
!
IIU=SIZE(PT_LS,1)
IJU=SIZE(PT_LS,2)
ILU=SIZE(PT_LS,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTATION OF PRESSURE AND EXNER FUNCTION
!              ------------------------------------------
!
IF (SIZE(PB_LS)/=0) THEN   ! hybrid level

  ALLOCATE(ZPFLUX_LS(IIU,IJU,ILU),ZEXNFLUX_LS(IIU,IJU,ILU))
!
  ZPFLUX_LS(:,:,:)=SPREAD(SPREAD(PA_LS,1,IIU),2,IJU)*PP00_LS &
                +SPREAD(SPREAD(PB_LS,1,IIU),2,IJU)*SPREAD(PPS_LS,3,ILU)
!
  ZEXNFLUX_LS(:,:,:)=(ZPFLUX_LS(:,:,:)/XP00)**(XRD/XCPD)
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTATION OF EXNER FUNCTION AT MASS POINT
!              -------------------------------------------
!
  ALLOCATE(ZEXNMASS_LS(IIU,IJU,ILU))
!
  ZEXNMASS_LS(:,:,1:ILU-1)=(ZEXNFLUX_LS(:,:,1:ILU-1)-ZEXNFLUX_LS(:,:,2:ILU))            &
                      /(LOG(ZEXNFLUX_LS(:,:,1:ILU-1))-LOG(ZEXNFLUX_LS(:,:,2:ILU)))
!
  ZEXNMASS_LS(:,:,ILU)   =(ZPFLUX_LS(:,:,ILU)/2./XP00)**(XRD/XCPD)
!
  PPMASS_LS(:,:,:)=XP00*(ZEXNMASS_LS(:,:,:))**(XCPD/XRD)

ELSE
  PPMASS_LS(:,:,:)=SPREAD(SPREAD(PA_LS,1,IIU),2,IJU)
  ALLOCATE(ZEXNMASS_LS(IIU,IJU,ILU))
  ZEXNMASS_LS(:,:,:)=(PPMASS_LS(:,:,:)/XP00)**(XRD/XCPD)

  ALLOCATE(ZEXNFLUX_LS(IIU,IJU,ILU))
  ZEXNFLUX_LS(:,:,1:ILU-1)=(ZEXNMASS_LS(:,:,1:ILU-1)-ZEXNMASS_LS (:,:,2:ILU))            &
                        /(LOG(ZEXNMASS_LS(:,:,1:ILU-1))-LOG(ZEXNMASS_LS (:,:,2:ILU)))
  ZEXNFLUX_LS(:,:,ILU)   =(PPMASS_LS(:,:,ILU)/2./XP00)**(XRD/XCPD)
END IF!-------------------------------------------------------------------------------
!
!*       3.    COMPUTATION OF RELATIVE HUMIDITY
!              --------------------------------
!
!*       3.1   ECMWF or ARPEGE/ALADIN case
!              ---------------------------
!
IF (ANY( ABS(PA_LS(:))>1.E-12)) THEN
  PHU_LS(:,:,:)= 100.*PPMASS_LS(:,:,:)                         &
                /(XRD/XRV*(1./MAX(PQ_LS(:,:,:,1),1.E-12)-1.)+1.) &
                /SM_FOES(PT_LS(:,:,:))
!
!*       3.2   PERIDOT case
!              ------------
ELSE
  PHU_LS(:,:,:)= 100.*PPMASS_LS(:,:,:)                         &
                /(XRD/XRV*(1./MAX(PQ_LS(:,:,:,1),1.E-12)-1.)+1.) &
                /(610.78*EXP(17.269*(PT_LS(:,:,:)-273.16)/(PT_LS(:,:,:)-35.86)))
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.    COMPUTATION OF RV AND THETAV
!              ----------------------------
!
DO JRR=2,SIZE(PR_LS,4)
  PR_LS(:,:,:,JRR) =  1. / (1./MAX(PQ_LS(:,:,:,JRR),1.E-12) - 1.)
END DO
!
PR_LS(:,:,:,1)=SM_PMR_HU(PPMASS_LS(:,:,:),                              &
                         PT_LS(:,:,:)*(1.+(XRV/XRD-1.)*PQ_LS(:,:,:,1)), &
                         PHU_LS(:,:,:),PR_LS(:,:,:,:),KITERMAX=100)
!
PTHV_LS(:,:,:)=PT_LS(:,:,:)/ZEXNMASS_LS(:,:,:) &
              *(1.+XRV/XRD*PR_LS(:,:,:,1))/(1.+WATER_SUM(PR_LS(:,:,:,:)))
!
!-------------------------------------------------------------------------------
!
!*       5.    COMPUTATION OF THE ALTITUDE OF THE PRESSURE POINTS
!              --------------------------------------------------
!
PZFLUX_LS(:,:,1)=PZS_LS(:,:)
DO JL=1,ILU-1
  PZFLUX_LS(:,:,JL+1)=PZFLUX_LS(:,:,JL)-XCPD*PTHV_LS(:,:,JL)/XG *           &
                               (ZEXNFLUX_LS(:,:,JL+1)-ZEXNFLUX_LS(:,:,JL))
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       6.    COMPUTATION OF THE ALTITUDE OF THE MASS POINTS
!              ----------------------------------------------
!
PZMASS_LS(:,:,2:ILU)  =PZFLUX_LS(:,:,2:ILU)    -   XCPD/XG                    &
                     *(0.75*PTHV_LS(:,:,2:ILU)+0.25*PTHV_LS(:,:,1:ILU-1)) &
                     *(ZEXNMASS_LS(:,:,2:ILU)-ZEXNFLUX_LS(:,:,2:ILU))
PZMASS_LS(:,:,1)     =PZFLUX_LS(:,:,2)         -   XCPD/XG                    &
                     *(0.75*PTHV_LS(:,:,1)+0.25*PTHV_LS(:,:,2))           &
                     *(ZEXNMASS_LS(:,:,1)-ZEXNFLUX_LS(:,:,2))
!
!-------------------------------------------------------------------------------
!
!*       8.    VERTICAL WIND
!              -------------
!
IF (PRESENT(PW_LS) .AND. SIZE(PB_LS)==0) THEN
  PW_LS(:,:,:)=0.  ! NCEP case not treated

ELSEIF (PRESENT(PW_LS)) THEN
!
!*       8.0   allocations
!              -----------
!
  ALLOCATE(ZSURFCOR_LS(IIU,IJU,ILU))
  ALLOCATE(ZRHOMASS_LS(IIU,IJU,ILU),ZRHO_LS(IIU,IJU,ILU))
  ALLOCATE(ZA_PRIME_LS(ILU),ZB_PRIME_LS(ILU))
  ALLOCATE(ZDPDETA_LS(IIU,IJU,ILU))
  ALLOCATE(ZHDIV_LS(IIU,IJU,ILU),ZINTDIV_LS(IIU,IJU,ILU))
!
!*       8.1   correction because of orography
!              -------------------------------
!
  ZSURFCOR_LS(:,:,:)=0.
  ZSURFCOR_LS(IIB:IIE,IJB:IJE,:)= PU_LS(IIB:IIE,IJB:IJE,:)                     &
                           *SPREAD(                                            &
                            (PZS_LS(IIB+1:IIE+1,IJB:IJE)-PZS_LS(IIB-1:IIE-1,IJB:IJE))    &
                            /SPREAD(XXHAT(IIB+1:IIE+1)-XXHAT(IIB-1:IIE-1),2,IJE-IJB+1)       &
                            *XMAP(IIB:IIE,IJB:IJE)                             &
                            ,3,ILU)                                            &
                          + PV_LS(IIB:IIE,IJB:IJE,:)                           &
                           *SPREAD(                                            &
                            (PZS_LS(IIB:IIE,IJB+1:IJE+1)-PZS_LS(IIB:IIE,IJB-1:IJE-1))    &
                           /SPREAD(XYHAT(IJB+1:IJE+1)-XYHAT(IJB-1:IJE-1),1,IIE-IIB+1)        &
                           *XMAP(IIB:IIE,IJB:IJE)                              &
                            ,3,ILU)
!
  DO JI=1,JPHEXT
     ZSURFCOR_LS(IIB-JI,:,1)=2.*ZSURFCOR_LS(IIB-JI+1,:,1)-ZSURFCOR_LS(IIB-JI+2,:,1)
     ZSURFCOR_LS(IIE+JI,:,1)=2.*ZSURFCOR_LS(IIE+JI-1,:,1)-ZSURFCOR_LS(IIE+JI-2,:,1)
     ZSURFCOR_LS(:,IJB-JI,1)=2.*ZSURFCOR_LS(:,IJB-JI+1,1)-ZSURFCOR_LS(:,IJB-JI+2,1)
     ZSURFCOR_LS(:,IJE+JI,1)=2.*ZSURFCOR_LS(:,IJE+JI-1,1)-ZSURFCOR_LS(:,IJE+JI-2,1)
  END DO
!!$  ZSURFCOR_LS( 1 , : ,1)=2.*ZSURFCOR_LS(  2  ,  :  ,1)-ZSURFCOR_LS(  3  ,  :  ,1)
!!$  ZSURFCOR_LS(IIU, : ,1)=2.*ZSURFCOR_LS(IIU-1,  :  ,1)-ZSURFCOR_LS(IIU-2,  :  ,1)
!!$  ZSURFCOR_LS( : , 1 ,1)=2.*ZSURFCOR_LS(  :  ,  2  ,1)-ZSURFCOR_LS(  :  ,  3  ,1)
!!$  ZSURFCOR_LS( : ,IJU,1)=2.*ZSURFCOR_LS(  :  ,IJU-1,1)-ZSURFCOR_LS(  :  ,IJU-2,1)
!
!*       8.2   mass density
!              ------------
!
  ZRHOMASS_LS(:,:,:)=PPMASS_LS(:,:,:)                                &
    /(XRD*PT_LS(:,:,:)*(1.+XRV/XRD*PR_LS(:,:,:,1))/(1.+WATER_SUM(PR_LS(:,:,:,:))) )
!
  ZRHO_LS(:,:,2:ILU)=0.5*(ZRHOMASS_LS(:,:,2:ILU)+ZRHOMASS_LS(:,:,1:ILU-1))
  ZRHO_LS(:,:,1)    =1.5*ZRHOMASS_LS(:,:,1)-0.5*ZRHOMASS_LS(:,:,2)
!
!*       8.3   derivatives of eta function (eta going from 1 to ILU)
!              ---------------------------
!
  ZA_PRIME_LS(1:ILU-1)= PA_LS(2:ILU)-PA_LS(1:ILU-1)
  ZA_PRIME_LS(ILU)    =             -PA_LS(  ILU)
  ZB_PRIME_LS(1:ILU-1)= PB_LS(2:ILU)-PB_LS(1:ILU-1)
  ZB_PRIME_LS(ILU)    =             -PB_LS(  ILU)
!
!*       8.4   term  (A'p00 + B'ps)
!              --------------------
!
  ZDPDETA_LS(:,:,:)= SPREAD(SPREAD(ZA_PRIME_LS(:),1,IIU),2,IJU)*PP00_LS  &
                    +SPREAD(SPREAD(ZB_PRIME_LS(:),1,IIU),2,IJU)          &
                     *SPREAD(PPS_LS(:,:),3,ILU)
!
!*       8.5   horizontal divergence  d(A'p00 + B'ps)u/dx + d(A'p00 + B'ps)v/dy
!              ---------------------
!
  ZHDIV_LS(:,:,:)=0.
  ZHDIV_LS(IIB:IIE,IJB:IJE,:)=(ZDPDETA_LS(IIB+1:IIE+1,IJB:IJE,:)*PU_LS(IIB+1:IIE+1,IJB:IJE,:)       &
                              -ZDPDETA_LS(IIB-1:IIE-1,IJB:IJE,:)*PU_LS(IIB-1:IIE-1,IJB:IJE,:))  &
                             /SPREAD(SPREAD(XXHAT(IIB+1:IIE+1)-XXHAT(IIB-1:IIE-1),2,IJE-IJB+1),3,ILU)  &
                             *SPREAD(XMAP(IIB:IIE,IJB:IJE),3,ILU)                        &
                            +(ZDPDETA_LS(IIB:IIE,IJB+1:IJE+1,:)*PV_LS(IIB:IIE,IJB+1:IJE+1,:)        &
                             -ZDPDETA_LS(IIB:IIE,IJB-1:IJE-1,:)*PV_LS(IIB:IIE,IJB-1:IJE-1,:))   &
                            /SPREAD(SPREAD(XYHAT(IJB+1:IJE+1)-XYHAT(IJB-1:IJE-1),1,IIE-IIB+1),3,ILU)   &
                            *SPREAD(XMAP(IIB:IIE,IJB:IJE),3,ILU)
!
  DO JI=1,JPHEXT
     ZHDIV_LS(IIB-JI,:,:)=2.*ZHDIV_LS(IIB-JI+1,:,:)-ZHDIV_LS(IIB-JI+2,:,:)
     ZHDIV_LS(IIE+JI,:,:)=2.*ZHDIV_LS(IIE+JI-1,:,:)-ZHDIV_LS(IIE+JI-2,:,:)
     ZHDIV_LS(:,IJB-JI,:)=2.*ZHDIV_LS(:,IJB-JI+1,:)-ZHDIV_LS(:,IJB-JI+2,:)
     ZHDIV_LS(:,IJE+JI,:)=2.*ZHDIV_LS(:,IJE+JI-1,:)-ZHDIV_LS(:,IJE+JI-2,:)
  END DO
!!$  ZHDIV_LS( 1 , : ,:)=2.*ZHDIV_LS(  2  ,  :  ,:)-ZHDIV_LS(  3  ,  :  ,:)
!!$  ZHDIV_LS(IIU, : ,:)=2.*ZHDIV_LS(IIU-1,  :  ,:)-ZHDIV_LS(IIU-2,  :  ,:)
!!$  ZHDIV_LS( : , 1 ,:)=2.*ZHDIV_LS(  :  ,  2  ,:)-ZHDIV_LS(  :  ,  3  ,:)
!!$  ZHDIV_LS( : ,IJU,:)=2.*ZHDIV_LS(  :  ,IJU-1,:)-ZHDIV_LS(  :  ,IJU-2,:)
!
!*       8.6   Integration of horizontal divergence
!              ------------------------------------
!
  ZINTDIV_LS(:,:,ILU)=-ZHDIV_LS(:,:,ILU)
!
  DO JL=ILU-1,1,-1
    ZINTDIV_LS(:,:,JL)= ZINTDIV_LS(:,:,JL+1) - ZHDIV_LS(:,:,JL)
  END DO
!
!*       8.7   Vertical velocity at all levels
!              -------------------------------
!
  DO JL=1,ILU
    PW_LS(:,:,JL)=1./ZRHO_LS(:,:,JL)/XG                            &
              * ( ZINTDIV_LS(:,:,JL) - PB_LS(JL)*ZINTDIV_LS(:,:,1) &
                + ZRHO_LS(:,:,1)*PB_LS(JL)*ZSURFCOR_LS(:,:,JL)*XG)
  END DO
END IF
!-------------------------------------------------------------------------------
!
WRITE(TLUOUT0%NLU,*) 'Routine CHANGE_GRIBEX_VAR completed'
!
END SUBROUTINE CHANGE_GRIBEX_VAR
