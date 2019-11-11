!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############################
      MODULE MODI_EDDYUV_FLUX_ONE_WAY_n
!     ###############################
!
INTERFACE
!
      SUBROUTINE EDDYUV_FLUX_ONE_WAY_n (KMI,KTCOUNT,KDXRATIO,KDYRATIO,HLBCX,HLBCY)
!
!
INTEGER, INTENT(IN) :: KMI     ! Model index
INTEGER, INTENT(IN) :: KTCOUNT ! iteration count
!
INTEGER, INTENT(IN) :: KDXRATIO   ! x and y-direction resolution RATIO
INTEGER, INTENT(IN) :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
!
END SUBROUTINE EDDYUV_FLUX_ONE_WAY_n
!
END INTERFACE
!
END MODULE MODI_EDDYUV_FLUX_ONE_WAY_n
!
!     ##################################################################################
      SUBROUTINE EDDYUV_FLUX_ONE_WAY_n (KMI,KTCOUNT,KDXRATIO,KDYRATIO,HLBCX,HLBCY)
!     ##################################################################################
!
!!    PURPOSE
!!    -------
!!      In case of 2D transect (latitude, altitude) grid-nesting models
!!      Barotropic fluxes (v'u') from the model 1 interpolated for the son models
!!
!!**  METHOD
!!    ------
!!
!!    IMPLICIT ARGUMENT
!!    -----------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	  M.Tomasini          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  07/07/11
!!      J.Escobar 2/05/2016 : bug in use of global/local bounds for call of BIKHARDT
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!     ##################################################################################
!
USE MODD_DEF_EDDYUV_FLUX_n
USE MODD_FIELD_n,               ONLY:XRVS
USE MODD_GRID_n
USE MODD_REF_n,                 ONLY:XRHODJ

USE MODD_METRICS_n
USE MODI_GRADIENT_U
!
! For the horizontal interpolation
USE MODI_BIKHARDT
USE MODD_BIKHARDT_n
USE MODD_NESTING
!
USE MODE_FIELD, ONLY : TFIELDLIST, FIND_FIELD_ID_FROM_MNHNAME
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KMI     ! Model index
INTEGER, INTENT(IN) :: KTCOUNT ! iteration count
!
INTEGER, INTENT(IN) :: KDXRATIO   ! x and y-direction resolution RATIO
INTEGER, INTENT(IN) :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions

!
!*       0.2   Declarations of local variables :
!
INTEGER:: IIB,IJB        ! Begining useful area  in x,y directions
INTEGER:: IIE,IJE        ! End useful area in x,y directions

INTEGER:: ISYNCHRO       ! model synchronic index relative to the model 1
                         ! = 1 for the first time step in phase with the model 1
                         ! = 0 for the last  time step (out of phase)
INTEGER:: JMI            ! Models loop
INTEGER:: IDTRATIO_KMI_1 ! Ratio between the time step of the son and the model 1

REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFLUX2 ! Work array=Dad interpolated flux field
                                              ! on the son grid
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDIV_UV! Work array=DIV of ZFLUX2
INTEGER :: IKU
!
INTEGER :: IDIMX,IDIMY
INTEGER :: IID, IRESP
!-------------------------------------------------------------------------------
!
!
! test of temporal synchronisation between the model 1 and the son KMI
!
IKU=SIZE(XZHAT)
IDTRATIO_KMI_1=1
DO JMI=2,KMI
   IDTRATIO_KMI_1=IDTRATIO_KMI_1*NDTRATIO(JMI)
END DO
ISYNCHRO = MODULO (KTCOUNT, IDTRATIO_KMI_1)
!
IF (ISYNCHRO==1 .OR. IDTRATIO_KMI_1 == 1) THEN

   CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
   ALLOCATE(ZFLUX2(SIZE(XRVS,1),SIZE(XRVS,2),SIZE(XRVS,3)))
   ALLOCATE(ZDIV_UV(SIZE(XRVS,1),SIZE(XRVS,2),SIZE(XRVS,3)))
   ZDIV_UV = 0.

   ! v'u' (EDDY_FLUX_MODEL(1)%XVU_FLUX_M) of model1 interpolation on the son grid put into ZFLUX2
   ZFLUX2 = 0.
   !IDIMX = SIZE(EDDYUV_FLUX_MODEL(1)%XVU_FLUX_M,1)
   !IDIMY = SIZE(EDDYUV_FLUX_MODEL(1)%XVU_FLUX_M,2)
   CALL FIND_FIELD_ID_FROM_MNHNAME('VU_FLX',IID,IRESP)
   IDIMX = SIZE(TFIELDLIST(IID)%TFIELD_X3D(1)%DATA,1)
   IDIMY = SIZE(TFIELDLIST(IID)%TFIELD_X3D(1)%DATA,2)
   CALL BIKHARDT (XBMX1,XBMX2,XBMX3,XBMX4,XBMY1,XBMY2,XBMY3,XBMY4, &
                  XBFX1,XBFX2,XBFX3,XBFX4,XBFY1,XBFY2,XBFY3,XBFY4, &
                  2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,&
!                  HLBCX,HLBCY,EDDYUV_FLUX_MODEL(1)%XVU_FLUX_M,ZFLUX2)
                  HLBCX,HLBCY,TFIELDLIST(IID)%TFIELD_X3D(1)%DATA,ZFLUX2)

   ! Lateral boundary conditions 
   ZFLUX2(IIB,:,:)    = ZFLUX2(IIB+1,:,:)
   ZFLUX2(1,:,:)      = ZFLUX2(IIB,:,:)
   ZFLUX2(IIE,:,:)    = ZFLUX2(IIE-1,:,:) 
   ZFLUX2(IIE+1,:,:)  = ZFLUX2(IIE,:,:)

   ZDIV_UV(:,:,:) = GX_U_M(1,IKU,1,ZFLUX2,XDXX,XDZZ,XDZX)

   ! Lateral boundary conditions
   ZDIV_UV(IIB,:,:)  =0.0
   ZDIV_UV(1,:,:)    =0.0
   ZDIV_UV(IIE,:,:)  =0.0
   ZDIV_UV(IIE+1,:,:)=0.0

   XRVS_EDDY_FLUX(:,:,:) = - XRHODJ(:,:,:)* ZDIV_UV(:,:,:)
        
   DEALLOCATE(ZFLUX2)

ENDIF

! COMPUTE NEW zonal wind at each son time step
! ---------------------------------------------
XRVS(:,:,:) = XRVS(:,:,:) + XRVS_EDDY_FLUX(:,:,:)

END SUBROUTINE EDDYUV_FLUX_ONE_WAY_n
