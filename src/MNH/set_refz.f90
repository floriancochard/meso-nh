!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_SET_REFZ
!     ####################
INTERFACE
      SUBROUTINE SET_REFZ(PTHV,PRV)
!
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PTHV ! thetav on the MESO-NH grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRV  ! rv on the MESO-NH grid
!
END SUBROUTINE SET_REFZ
END INTERFACE
END MODULE MODI_SET_REFZ
!     #############################
      SUBROUTINE SET_REFZ(PTHV,PRV)
!     #############################
!
!!****  *SET_REFZ* - computation of the anelastic reference state profiles
!!
!!    PURPOSE
!!    -------
!!    This routine computes the profiles of virtual potential temperature
!!    and vapor mixing ratio for the anelastic reference state from the
!!    virtual potential temperature, the mixing ratio and the reference
!!    state Exner function at model top.
!!
!!**  METHOD
!!    ------
!!                    
!!  1 The profiles of thetav_ref and rv_ref are computed on the mass
!!    points of the GS grid defined by XZHAT and no orography. 
!!    For each level of this vertical grid, the horizontal mean are computed 
!!    in the function ZSECT.
!!    If necessary (min(zs)>0) the profiles are extrapolated to z=0:
!!     for thetav with a uniform gradient of 3.5E-3 K/m
!!     for rv with a exponential function a*exp(-bz) computed by using the
!!     two lowest points availables.
!!
!!    CAUTION: all the profiles in this routine are computed on the whole
!!    vertical, i.e. also for the non-physical points (with the same 
!!    extrapolation as above for levels under the ground, and with linear
!!    extrapolation for the upper level).
!!    
!!  2 PI_ref is computed by integration of the hydrostatic relation for the
!!    reference state variables from top (XEXNTOP) to bottom. 
!!     ~
!!  3 PI_ref, the exner function at mass levels is computed as follows and
!!    linearly extrapolated with zhat for the uppest non-physical level.
!!
!!      ~           PIref(k+1)-PIref(k)
!!     PIref(k) = -----------------------
!!                lnPIref(k+1)-lnPIref(k)
!!
!!   
!!  3 rhod_ref is deduced from the relation
!!
!!                      ~     (Cpd/Rd -1)
!!                p00 (PI_ref)
!!    rhod_ref=  _________________________
!!
!!                Rd thetav_ref (1+rv_ref)
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    subroutine COMPUTE_EXNER_FROM_TOP : to compute hydrostatic Exner function
!!    subroutine ZSECT :to compute the mean of a 3D field at a constant altitude
!!    MZF              : Shuman operator
!!
!!    Module MODI_SHUMAN    : interface for Shuman operators
!!    module MODI_ZSECT     :contains interface for subroutine ZSECT
!!    Module MODI_COMPUTE_EXNER_FROM_TOP
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models. 
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_REF       : contains anelastic reference state variables
!!         XEXNTOP : reference state Exner function at model top
!!         XRHODREFZ: reference state profile of rhod
!!         XTHVREFZ: reference state profile of thetav
!!      Module MODD_GRID1     : contains grid variables for model1
!!         XZZ   : altitude of the w points 
!!         XZHAT : GS levels 
!!      Module MODD_CST       : contains physical constants
!!         XG  : gravity constant
!!         XCPD: specific heat for dry air at constant pressure
!!         XP00: reference pressure
!!         XRD : gas constant for dry air
!!      Module MODD_PARAMETERS
!!         JPVEXT
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
!!      Original    12/12/94
!!                  06/03/96 (V. Masson) add call to COMPUTE_EXNER_FROM_TOP
!!                  29/10/96 (V. Masson) add positivity control on rv
!!                  26/10/10 (G. Tanguy) add control on rv : if equal to 0
!!                                       at MINLEVEL, keep 0 for lower levels
!!                                       (for ideal case)
!!                  2014 (M.Faivre)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_GRID_n
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_REF
!
USE MODE_MPPDB
!
USE MODI_COMPUTE_EXNER_FROM_TOP
USE MODI_SHUMAN
USE MODI_ZSECT
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PTHV ! thetav on the MESO-NH grid
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PRV  ! rv on the MESO-NH grid
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER                       :: ILUOUT0
INTEGER                       :: IRESP
INTEGER                       :: IKB,IKE,IKU         ! vertical limit indices
INTEGER                       :: JK                  ! vertical loop 
INTEGER                       :: IMINLEVEL           ! first level with physical
!                                                    ! field mean after ZSECT
REAL,DIMENSION(SIZE(XZHAT))   :: ZRREFZ,ZEXNMASSREFZ&! rv_ref, PI_ref at mass 
!                                                    ! point
                                ,ZEXNREFZ            ! PI_ref at flux point
REAL,DIMENSION(SIZE(XZZ,1),SIZE(XZZ,2),SIZE(XZZ,3))  &
                              :: ZZMASS              ! altitude of mass points
REAL                          :: ZTHVCLIMGR          ! gradient near the ground
!                                                    ! for extrapolation of thv
REAL                          :: ZCOEFA,ZCOEFB       ! coefficients of the 
!                                                    ! extrapolating function
!                                                    ! for rv
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATIONS
!              ---------------
!
IKU=SIZE(XZHAT)
IKB=JPVEXT+1
IKE=IKU-JPVEXT
!-------------------------------------------------------------------------------
!
!*       2.    ALTITUDE OF THE MASS POINTS
!              ---------------------------
!
ZZMASS(:,:,:)=MZF(1,IKU,1,XZZ(:,:,:))
ZZMASS(:,:,IKU)=1.5*XZZ(:,:,IKU)-0.5*XZZ(:,:,IKU-1)
!
!20131024 check zzmass and pthv
CALL MPPDB_CHECK3D(PTHV,"SET_REFZ:PTHV",PRECISION)
CALL MPPDB_CHECK3D(ZZMASS,"SET_REFZ:ZZMASS",PRECISION)
!                                                         
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTATION OF THE REFERENCE STATE PROFILE thetav_ref
!              -----------------------------------------------------
!
ZTHVCLIMGR=3.5E-3
ALLOCATE(XTHVREFZ(IKU))
XTHVREFZ(:)=-999.
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
!ocl scalar
!!!!!!!!!!!!!!!!  FUJI  compiler directive !!!!!!!!!!
DO JK=IKB,IKU-1
    XTHVREFZ(JK)= ZSECT(0.5*(XZHAT(JK)+XZHAT(JK+1)),             &
                        ZZMASS(:,:,IKB:IKE+1),PTHV(:,:,IKB:IKE+1))
END DO
    XTHVREFZ(IKU)=XTHVREFZ(IKU-1)                                 &
                 +(XTHVREFZ(IKU-1)-XTHVREFZ(IKU-2))               &
                 *(XZHAT(IKU)-XZHAT(IKU-1))                       &
                 /(0.5*(XZHAT(IKU)-XZHAT(IKU-2)) )                
!
IMINLEVEL=COUNT(XTHVREFZ(:)==-999.)+1
WHERE (XTHVREFZ==-999.)
    XTHVREFZ(:)=XTHVREFZ(IMINLEVEL)                                       &
                +ZTHVCLIMGR*0.5*( XZHAT(:)+EOSHIFT(XZHAT(:),1)            &
                                 -XZHAT(IMINLEVEL)-XZHAT(IMINLEVEL+1)) 
END WHERE
XTHVREFZ(1)=XTHVREFZ(2)
!
!-------------------------------------------------------------------------------
!                                                         
!*       4.    COMPUTATION OF THE REFERENCE STATE PROFILE rv_ref
!              -------------------------------------------------
!
ZRREFZ(:)=-999.
!ocl scalar
DO JK=IKB,IKU-1
    ZRREFZ(JK)= ZSECT(0.5*(XZHAT(JK)+XZHAT(JK+1)),              &
                      ZZMASS(:,:,IKB:IKE+1),PRV(:,:,IKB:IKE+1))
END DO
    ZRREFZ(IKU)=ZRREFZ(IKU-1)                                   &
               +(ZRREFZ(IKU-1)-ZRREFZ(IKU-2))                   &
               *(XZHAT(IKU)-XZHAT(IKU-1))                       &
               /(0.5*(XZHAT(IKU)-XZHAT(IKU-2)) )  
!
IMINLEVEL=COUNT(ZRREFZ(:)==-999.)+1
IF (ZRREFZ(IMINLEVEL)==0) THEN
  WHERE (ZRREFZ==-999.)
    ZRREFZ(:)=0
  END WHERE
ELSE
  ZCOEFB=-(LOG(ZRREFZ(IMINLEVEL+1))-LOG(ZRREFZ(IMINLEVEL))) &
       /(0.5*(XZHAT(IMINLEVEL+2)-XZHAT(IMINLEVEL)))
  ZCOEFA=ZRREFZ(IMINLEVEL)*EXP(ZCOEFB*0.5*(XZHAT(IMINLEVEL+1)+XZHAT(IMINLEVEL)))
  WHERE (ZRREFZ==-999.)
    ZRREFZ(:)=ZCOEFA*EXP(-ZCOEFB*0.5*(XZHAT(:)+EOSHIFT(XZHAT(:),1)))
  END WHERE
ENDIF
!
ZRREFZ(1)=ZRREFZ(2)
ZRREFZ(:)=MAX(ZRREFZ(:),0.)
!-------------------------------------------------------------------------------
!                             
!*       5.    COMPUTATION OF PI_ref
!              ---------------------
!
CALL COMPUTE_EXNER_FROM_TOP(XTHVREFZ,XZHAT,XEXNTOP,ZEXNREFZ,ZEXNMASSREFZ)
!
!-------------------------------------------------------------------------------
!                             
!*       6.    COMPUTATION OF rhod_ref
!              -----------------------
!
ALLOCATE(XRHODREFZ(IKU))
XRHODREFZ(:)=(XP00*(ZEXNMASSREFZ(:))**(XCPD/XRD-1.))                           &
            /(XRD*XTHVREFZ(:)*(1.+ZRREFZ(:)))
!
XRHODREFZ(1)=XRHODREFZ(2)
!-------------------------------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
WRITE(ILUOUT0,*)
WRITE(ILUOUT0,*) 'REFERENCE STATE VARIABLES WITHOUT OROGRAPHY'
WRITE(ILUOUT0,*) 'level   altitude     theta          rhod      pressure'
DO JK=1,IKU
  WRITE(ILUOUT0,9000) ' ',JK,'  ',XZHAT(JK),' ',XTHVREFZ(JK),' ',XRHODREFZ(JK),' ', XP00*(ZEXNREFZ(JK))**(XCPD/XRD)
ENDDO
9000 FORMAT(A1,I3.3,A2,F12.5,A1,F12.7,A1,F12.8,A1,F12.5)
!
!-------------------------------------------------------------------------------
!
WRITE(ILUOUT0,*) 'Routine SET_REFZ completed'
!
END SUBROUTINE SET_REFZ
