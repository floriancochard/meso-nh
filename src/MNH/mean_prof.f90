!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      MODULE MODI_MEAN_PROF
!     #####################
INTERFACE
      SUBROUTINE MEAN_PROF(PVAR_MX,PZMASS_MX,PZS_LS,PCLIMGR,&
                           PF_FREE,PZ_FREE)
!
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR_MX   ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! mass points altitude
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,                     INTENT(IN)  :: PCLIMGR   ! climatological gradient
!                                                  ! near the ground
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PF_FREE   ! mean profile of the
!                                                  ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PZ_FREE   ! discretization in x,y,z
!                                                  ! of the profile on the
!                                                  ! flat grid where zs is the
!                                                  ! minimum of both orographies
END SUBROUTINE MEAN_PROF
END INTERFACE
END MODULE MODI_MEAN_PROF
!     ##############################################################
      SUBROUTINE MEAN_PROF(PVAR_MX,PZMASS_MX,PZS_LS,PCLIMGR,&
                            PF_FREE,PZ_FREE)
!     ##############################################################
!
!!****  *MEAN_PROF* - Computation of the profile of the free atmospheres
!!                           i.e. without the Boundary layer structures
!!
!!    PURPOSE
!!    -------
!!    This routine computes the profile used for the shift of a variable
!!    and the altitude of the discretization points of this profile.
!
!!    CAUTION:
!!    The shift profile is only defined on the inner vertical points of the grid.
!!
!!**  METHOD
!!    ------
!!    The profile is discretized on the vertical GS grid defined by
!!    the MESO-NH level array XZHAT and by a constant orography,
!!    corresponding to the minimum of the Arpege and MESO-NH orographies.
!!    If necessary, the profile is extrapolated under the minimum
!!    altitude of the Arpege orography with a climatological vertical 
!!    gradient PCLIMGR (uniform on the whole domain).
!!
!!    EXTERNAL
!!    --------
!!
!!    function ZSECT   : to compute the mean of a 3D field at a constant 
!!                       altitude
!!    Module MODI_ZSECT: contains interface for function ZSECT
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CONF      : contains configuration variables for all models. 
!!         NVERB : verbosity level for output-listing
!!      Module MODD_LUNIT     :  contains logical unit names for all models
!!         CLUOUT0 : name of output-listing
!!      Module MODD_GRID1     : contains grid variables for model1
!!         XZS   : orography of MESO-NH
!!         XZHAT : GS levels
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
!!      Original    26/08/97
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_GRID_n
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
!
USE MODI_ZSECT
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR_MX   ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZMASS_MX ! mass points altitude
REAL,   DIMENSION(:,:),   INTENT(IN)  :: PZS_LS    ! large scale orography
REAL,                     INTENT(IN)  :: PCLIMGR   ! climatological gradient
!                                                  ! near the ground
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PF_FREE   ! mean profile of the
!                                                  ! thermodynamical field
REAL,   DIMENSION(:,:,:), INTENT(OUT) :: PZ_FREE   ! discretization in x,y,z
!                                                  ! of the profile on the
!                                                  ! flat grid where zs is the
!                                                  ! minimum of both orographies
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: ILEVEL,IKB,IKE,ILB,ILE,JK
REAL    :: ZMIN
REAL, DIMENSION(SIZE(PF_FREE,3)) :: ZF_FREE
REAL, DIMENSION(SIZE(PZ_FREE,3)) :: ZZ_FREE
!-------------------------------------------------------------------------------
!
!
IKB=JPVEXT+1
IKE=SIZE(PZ_FREE,3)-JPVEXT
!
!*       1.   Computation of the altitude of the GS grid for the shift profile
!             ----------------------------------------------------------------
!
ZMIN=MIN(MINVAL(PZS_LS),MINVAL(XZS))
ZZ_FREE(1:IKE)=ZMIN+0.5*(XZHAT(1:IKE)+XZHAT(2:IKE+1))*(1.-ZMIN/XZHAT(IKE+1))
ZZ_FREE(IKE+1)=2.*XZHAT(IKE+1)-ZZ_FREE(IKE)
!
!-------------------------------------------------------------------------------
!
!*       2.   Computation of the shift profile
!             --------------------------------
!
!*       2.1  Defined values
!             --------------
!
ZF_FREE(:)=-999.
DO JK=IKB,IKE
    ZF_FREE(JK)=ZSECT(ZZ_FREE(JK),PZMASS_MX(:,:,IKB:IKE+1),&
                         PVAR_MX(:,:,IKB:IKE+1))
END DO
!
!*       2.2  Low levels values (extrapolation with constant gradient)
!             --------------------------------------------------------
!
ILEVEL=0
DO JK=1,IKE
  IF (ABS(ZF_FREE(JK)+999.)<1.E-10) THEN
    ILEVEL=JK+1
  ELSE
    EXIT
  ENDIF
ENDDO
!
DO JK=1,ILEVEL-1
    ZF_FREE(JK)=ZF_FREE(ILEVEL)&
                +PCLIMGR*(ZZ_FREE(JK)-ZZ_FREE(ILEVEL))
ENDDO
!
!*       2.3  Upper levels values (linear extrapolation)
!             ------------------------------------------
!
ILEVEL=IKE+1
DO JK=IKE+1,1,-1
  IF (ABS(ZF_FREE(JK)+999.)<1.E-10) THEN
    ILEVEL=JK-1
  ELSE
    EXIT
  ENDIF
ENDDO
!
DO JK=IKE+1,ILEVEL+1
    ZF_FREE(JK)=ZF_FREE(ILEVEL)                          &
                  +(ZF_FREE(ILEVEL)-ZF_FREE(ILEVEL-1))   &
                  /(ZZ_FREE(ILEVEL)-ZZ_FREE(ILEVEL-1))   &
                  *(ZZ_FREE(JK)-ZZ_FREE(ILEVEL))
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       3.   3D output profiles arrays
!             -------------------------
!
PZ_FREE(:,:,:)=SPREAD(SPREAD(ZZ_FREE(:),1,SIZE(PZ_FREE,1)),2,SIZE(PZ_FREE,2))
PF_FREE(:,:,:)=SPREAD(SPREAD(ZF_FREE(:),1,SIZE(PF_FREE,1)),2,SIZE(PF_FREE,2))
!
!-------------------------------------------------------------------------------
!
WRITE(TLUOUT0%NLU,*) 'Routine MEAN_PROF completed'
!
END SUBROUTINE MEAN_PROF
