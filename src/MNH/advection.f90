!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #####################
      MODULE MODI_ADVECTION
!     #####################
!
INTERFACE
      SUBROUTINE ADVECTION (HUVW_ADV_SCHEME,HMET_ADV_SCHEME,HSV_ADV_SCHEME,    &
                           KLITER, HLBCX, HLBCY,KRR, KSV, KTCOUNT,             &
                           PTSTEP_MET, PTSTEP_SV,                              & 
                           PUM, PVM, PWM, PTHM, PRM, PTKEM, PSVM,              &
                           PUT, PVT, PWT, PTHT, PRT, PTKET, PSVT,              &
                           PRHODJ, PDXX, PDYY, PDZZ, PDZX, PDZY,               &
                           PRUS,PRVS, PRWS, PRTHS, PRRS, PRTKES, PRSVS         )
!
!
CHARACTER(LEN=6),         INTENT(IN)    :: HMET_ADV_SCHEME, & ! Control of the 
                                           HSV_ADV_SCHEME,  & ! scheme applied 
                                           HUVW_ADV_SCHEME     ! to the selected
                                                              ! variables 
!
INTEGER,                  INTENT(IN)    :: KLITER        ! Iteration number for
                                                         ! the MPDATA scheme
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! iteration count
REAL,                     INTENT(IN)    :: PTSTEP_MET !  Effective time step for
                                                ! meteorological scalar variables 
                                                ! (depending on advection scheme)
REAL,                     INTENT(IN)    :: PTSTEP_SV !  Effective time step for
                                                ! tracer scalar variables 
                                                ! (depending on advection scheme)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRM , PSVM
                                                  ! Variables at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUM, PVM, PWM
                                                  ! Variables at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT , PVT  , PWT, PRHODJ
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT , PSVT
                                                  ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  !  metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS , PRVS  , PRWS
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS , PRSVS
                                                  ! Sources terms 
!
!
END SUBROUTINE ADVECTION
!
END INTERFACE
!
END MODULE MODI_ADVECTION 
!     ##########################################################################
      SUBROUTINE ADVECTION (HUVW_ADV_SCHEME,HMET_ADV_SCHEME,HSV_ADV_SCHEME,    &
                           KLITER, HLBCX, HLBCY,KRR, KSV, KTCOUNT,             &
                           PTSTEP_MET, PTSTEP_SV,                              & 
                           PUM, PVM, PWM, PTHM, PRM, PTKEM, PSVM,              &
                           PUT, PVT, PWT, PTHT, PRT, PTKET, PSVT,              &
                           PRHODJ, PDXX, PDYY, PDZZ, PDZX, PDZY,               &
                           PRUS,PRVS, PRWS, PRTHS, PRRS, PRTKES, PRSVS         )
!     ##########################################################################
!
!!****  *ADVECTION * - routine to call the specialized advection routines
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to control the advection routines.
!!    For that, it is first necessary to compute the metric coefficients
!!    and the contravariant components of the momentum.
!!
!!**  METHOD
!!    ------
!!      The advection of momenta is calculated using a centred (second order) 
!!    scheme. Three schemes can be used to calculated the advection of a 
!!    scalar: centred (2nd) (ADVECSCALAR), Flux-Corrected Transport Scalar
!!    (FCT_SCALAR) and a Multidimensional Positive Definite Advection Transport
!!    Algorithm (MPDATA).
!!      Once the scheme is selected, it is applied to the following group of
!!    variables: METeorologicals (temperature, water substances, TKE,
!!    dissipation TKE) and Scalar Variables. It is possible to select different
!!    advection schemes for each group of variables.
!!
!!    EXTERNAL
!!    --------
!!      Functions MXM,MYM,MZM  : computes the averages along the 3 directions
!!      CONTRAV              : computes the contravariant components.
!!      ADVECUVW             : computes the advection terms for momentum.
!!      ADVECSCALAR          : computes the advection terms for scalar fields.
!!      ADD3DFIELD_ll        : add a field to 3D-list
!!      ADVEC_4TH_ORDER      : 4th order advection scheme
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book1 and book2 ( routine ADVECTION )
!!
!!    AUTHOR
!!    ------
!!	J.-P. Pinty      * Laboratoire d'Aerologie*
!!	J.-P. Lafore     * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/07/94 
!!                  01/04/95 (Ph. Hereil J. Nicolau) add the model number
!!                  23/10/95 (J. Vila and JP Lafore) advection schemes scalar
!!                  16/01/97 (JP Pinty)              change presentation 
!!                  30/04/98 (J. Stein P Jabouille)  extrapolation for the cyclic
!!                                                   case and parallelisation
!!                  24/06/99 (P Jabouille)           case of NHALO>1
!!                  25/10/05 (JP Pinty)              4th order scheme
!!                  24/04/06 (C.Lac)                 Split scalar and passive
!!                                                   tracer routines
!!                  08/06    (T.Maric)               PPM scheme
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!-------------------------------------------------------------------------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER(LEN=6),         INTENT(IN)    :: HMET_ADV_SCHEME, & ! Control of the 
                                           HSV_ADV_SCHEME,  & ! scheme applied 
                                           HUVW_ADV_SCHEME     ! to the selected
                                                              ! variables 
!
INTEGER,                  INTENT(IN)    :: KLITER        ! Iteration number for
                                                         ! the MPDATA scheme
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! iteration count
REAL,                     INTENT(IN)    :: PTSTEP_MET !  Effective time step for
                                                ! meteorological scalar variables 
                                                ! (depending on advection scheme)
REAL,                     INTENT(IN)    :: PTSTEP_SV !  Effective time step for
                                                ! tracer scalar variables 
                                                ! (depending on advection scheme)
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PRM , PSVM
                                                  ! Variables at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUM, PVM, PWM
                                                  ! Variables at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT , PVT  , PWT
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET, PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT , PSVT
                                                  ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  !  metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS , PRVS, PRWS
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS , PRSVS
                                                  ! Sources terms 
!
!
! ROUTINE TO REMOVE
!-------------------------------------------------------------------------------
!
END SUBROUTINE ADVECTION
