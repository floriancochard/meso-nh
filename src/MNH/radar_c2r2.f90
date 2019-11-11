!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ##########################
       MODULE MODI_RADAR_C2R2 
!      ##########################
!
INTERFACE
      SUBROUTINE RADAR_C2R2(PRT,PSVT,PRHODREF,PTEMP,PRARE,  &
                               PVDOP,PHHRE,PVVRE)
!
REAL,  DIMENSION(:,:,:,:), INTENT(IN)  :: PRT  ! microphysical  mix. ratios at t
REAL,  DIMENSION(:,:,:,:), INTENT(IN)  :: PSVT  ! microphysical concent. at t
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PRHODREF ! density of the ref. state
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PTEMP    ! air temperature
!
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRARE! radar reflectivity in dBZ
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PVDOP! radar Doppler fall speed
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PHHRE! hori. pol. reflectivity (mm6/m3)
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PVVRE! vert. pol. reflectivity (mm6/m3)
!
END SUBROUTINE RADAR_C2R2
!
END INTERFACE
!
END MODULE MODI_RADAR_C2R2
!     #########################################################################
      SUBROUTINE RADAR_C2R2(PRT,PSVT,PRHODREF,PTEMP,PRARE,  &
                               PVDOP,PHHRE,PVVRE)
!     #########################################################################
!
!!****  *RADAR_C2R2 * - computes some pertinent radar parameters
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the equivalent reflectivity,
!!    the Doppler reflectivity and the H and V polarized reflectivities of a 
!!    mixed phase cloud.
!!
!!**  METHOD
!!    ------
!!      The reflectivities are computed using the n(D) * D**6 formula. The 
!!    equivalent reflectiviy is the sum of the reflectivity produced by the
!!    the raindrops and the equivalent reflectivities of the ice crystals.
!!    The latter are computed using the melted diameter. The Doppler 
!!    reflectivity is the 'fall speed'-moment of individual particle
!!    reflectivity. Ice crystal are assumed to have no preferred orientation.
!!    the Z_VV formula is taken from Brandes et al. (MWR, 1995).
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XRHOLW               ! Liquid water density
!!      Module MODD_RAIN_DESCR
!!      Module MODD_RAIN_PARAM
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine RADAR_C2R2 )
!!      Smith P.L., 1984: Equivalent Radar Reflectivity Factors for Snow and
!!                        Ice Particles, JCAM, 23, 1258-1260.
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/05/96
!!                  03/12/96 change the arg. list
!!                  13/10/98 remove function statement
!!                  2014 G.Delautier : remplace MODD_RAIN_C2R2_PARAM par MODD_RAIN_C2R2_KHKO_PARAM
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_PARAM_C2R2
USE MODD_RAIN_C2R2_KHKO_PARAM
USE MODD_RAIN_C2R2_DESCR
USE MODD_REF
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,  DIMENSION(:,:,:,:), INTENT(IN)  :: PRT  ! microphysical  mix. ratios at t
REAL,  DIMENSION(:,:,:,:), INTENT(IN)  :: PSVT  ! microphysical concent. at t
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PRHODREF ! density of the ref. state
REAL,  DIMENSION(:,:,:),   INTENT(IN)  :: PTEMP    ! air temperature
!
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PRARE! radar reflectivity in dBZ
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PVDOP! radar Doppler fall speed
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PHHRE! hori. pol. reflectivity (mm6/m3)
REAL,  DIMENSION(:,:,:), INTENT(OUT) :: PVVRE! vert. pol. reflectivity (mm6/m3)
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB           ! Coordinates of the first physical points along z
INTEGER :: IND           ! Number of interval to integrate the kernels
REAL :: ZALPHA, ZNU, ZP  ! Parameters to compute the value of the p_moment
       			 ! of the generalized Gamma function
REAL :: ZDINFTY          ! Factor used to define the "infinite" diameter
!
REAL :: ZCXR=-1.0                     ! for rain N ~ 1/N_0 
                                      ! (in Kessler parameterization)
REAL :: ZSLOPE, ZINTERCEPT, ZEXPONENT ! parameters defining the mean axis ratio
                				      ! functionnal
REAL :: ZDMELT_FACT                   ! factor used to compute the equivalent
			                	      ! melted diameter
				                      ! water reflectivity (from Smith, JCAM 84)
REAL :: ZEXP                          ! anciliary parameter
REAL :: ZRHO00                        ! Surface reference air density
!
LOGICAL, DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: GRAIN
REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZLBDA 
                                      ! slope distribution parameter
REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZW
REAL,    DIMENSION(SIZE(PTEMP,1),SIZE(PTEMP,2),SIZE(PTEMP,3)) :: ZREFL_MELT_CONV
REAL, DIMENSION(:), ALLOCATABLE                               :: ZZLBDA
REAL, DIMENSION(:), ALLOCATABLE                               :: ZZVVRE
INTEGER                                                       :: JLBDA
INTEGER                                                       :: NLBDA
REAL                                                          :: ZFRAC_WATER
!
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZWLBDR
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZWLBDC
!
!-------------------------------------------------------------------------------
!
!
!*       1.     FUNCTION STATEMENTS
!   	        -------------------
!
!
!*       1.1    p_moment of the Generalized GAMMA function
!
! Recall that MOMG(ZALPHA,ZNU,ZP)=GAMMA(ZNU+ZP/ZALPHA)/GAMMA(ZNU)
!
!
!-------------------------------------------------------------------------------
!
!
!        2.     INTIALIZE OUTPUT LISTING AND OTHER ARRAYS
!               -----------------------------------------
!
!
PRARE(:,:,:) = 0.0    ! radar reflectivity
PVDOP(:,:,:) = 0.0    ! radar Doppler fall speed
PHHRE(:,:,:) = 0.0    ! Horizontally polarized radar reflectivity
PVVRE(:,:,:) = 0.0    ! Vertically polarized radar reflectivity
!
!-------------------------------------------------------------------------------
!
!
!*       3.     RAINDROPS
!               ---------
!
IF (SIZE(PRT,4) >= 3) THEN
!
!*       1.     COMPUTE THE SLOPE PARAMETER ZLBDR**6
!               ----------------------------------------
!
IKB = 1+JPVEXT
!
ALLOCATE(ZWLBDR(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)))
ALLOCATE(ZWLBDC(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)))
!
ZRHO00 = XP00/(XRD*XTHVREFZ(IKB))
ZWLBDR(:,:,:)  = XMNH_HUGE
!
WHERE (PRT(:,:,:,3)>XRTMIN(3) .AND. PSVT(:,:,:,3)>XCTMIN(3))
  ZWLBDR(:,:,:) = XLBR * PSVT(:,:,:,3) / (PRHODREF(:,:,:) * PRT(:,:,:,3))
  ZWLBDR(:,:,:)  = ZWLBDR(:,:,:)**2
END WHERE
!
ZWLBDC(:,:,:)  = XMNH_HUGE
WHERE (PRT(:,:,:,2)>XRTMIN(2) .AND. PSVT(:,:,:,2)>XCTMIN(2))
  ZWLBDC(:,:,:) = XLBC * PSVT(:,:,:,2) / (PRHODREF(:,:,:) * PRT(:,:,:,2))
  ZWLBDC(:,:,:)  = ZWLBDC(:,:,:)**2
END WHERE
!
    PRARE(:,:,:) = 1.E18*PSVT(:,:,:,3)/ZWLBDR(:,:,:)*       &
                         MOMG(XALPHAR,XNUR,6.0)
WHERE (PRARE(:,:,:) >= 1.0)
    PVDOP(:,:,:) = XCR*(ZRHO00/PRHODREF(:,:,:))**XCEXVT *             &
                   (ZWLBDR(:,:,:)**(-XDR/6.))                            &
                 *MOMG(XALPHAR,XNUR,6.0+XDR)/MOMG(XALPHAR,XNUR,6.0)
    PRARE(:,:,:) = 10.* LOG10(PRARE(:,:,:))
END WHERE
    PHHRE(:,:,:) = PRARE(:,:,:)
    PVVRE(:,:,:) = PRARE(:,:,:)
!
DEALLOCATE(ZWLBDR)
DEALLOCATE(ZWLBDC)
!
END IF
!
!
!-------------------------------------------------------------------------------
!
CONTAINS
!------------------------------------------------------------------------------
!
  FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA
!
  IMPLICIT NONE
!
  REAL     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL     :: PP     ! order of the moment
  REAL     :: PMOMG  ! result: moment of order ZP
!
!------------------------------------------------------------------------------
!
!
  PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
  END FUNCTION MOMG
!
!------------------------------------------------------------------------------
!
END SUBROUTINE RADAR_C2R2
