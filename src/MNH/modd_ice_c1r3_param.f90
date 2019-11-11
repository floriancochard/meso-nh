!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##########################
      MODULE MODD_ICE_C1R3_PARAM
!     ##########################
!
!!****  *MODD_ICE_C1R3_PARAM* - declaration of some microphysical factors
!!                              extensively used in the C1R3 cold scheme.
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24/11/00
!!      Jean-Pierre PINTY    29/ 6/01  Bug in RCHONI and RVHNCI
!!      Jean-Pierre PINTY    29/ 6/01  Add RHHONI process (freezing haze part.)
!!      Jean-Pierre PINTY    05/04/02  Change effective radius parameters
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
!
REAL,SAVE :: XFSEDRI,XFSEDCI,                  & ! Constants for sedimentation
    	     XFSEDS, XEXSEDS,                  & ! fluxes of R and C for Ice,
    	     XFSEDG, XEXSEDG                     ! of R for Snow and Graupel
!
REAL,SAVE :: XNUC_DEP,XEXSI_DEP,XEX_DEP,       & ! Constants for heterogeneous
    	     XNUC_CON,XEXTT_CON,XEX_CON,       & ! ice nucleation : DEP et CON
    	     XMNU0                               ! mass of nucleated ice crystal
!
REAL,SAVE :: XRHOI_HONH,XCEXP_DIFVAP_HONH,     & ! Constants for homogeneous
             XCOEF_DIFVAP_HONH,XRCOEF_HONH,    & ! haze freezing : HHONI
             XCRITSAT1_HONH,XCRITSAT2_HONH,    &
             XTMIN_HONH,XTMAX_HONH,            &
             XDLNJODT1_HONH,XDLNJODT2_HONH,    &
             XC1_HONH,XC2_HONH,XC3_HONH
!
REAL,SAVE :: XC_HONC,XR_HONC,                  & ! Constants for homogeneous
    	     XTEXP1_HONC,XTEXP2_HONC,          & ! droplet freezing : CHONI
             XTEXP3_HONC,XTEXP4_HONC,          &
             XTEXP5_HONC
!
REAL,SAVE :: XDICNVS_LIM, XLBDAICNVS_LIM,      & ! Constants for pristine ice
             XC0DEPIS,XC1DEPIS,                & ! deposition conversion to
             XR0DEPIS,XR1DEPIS                   ! snow : ICNVS
!
REAL,SAVE :: XCSCNVI_MAX, XLBDASCNVI_MAX,      &
             XRHORSMIN,                        &
             XDSCNVI_LIM, XLBDASCNVI_LIM,      & ! Constants for snow
             XC0DEPSI,XC1DEPSI,                & ! sublimation conversion to
             XR0DEPSI,XR1DEPSI                   ! pristine ice : SCNVI
!
REAL,SAVE :: XSCFAC,                           & ! Constants for the Bergeron
    	     X0DEPI,X2DEPI,                    & ! Findeisen process and
             X0DEPS,X1DEPS,XEX0DEPS,XEX1DEPS,  & ! deposition : DEP on S and
    	     X0DEPG,X1DEPG,XEX0DEPG,XEX1DEPG     !                  on G
!
REAL,SAVE :: XKER_ZRNIC_A1,XKER_ZRNIC_A2         ! Long-Zrnic Kernels
!
REAL,SAVE :: XSELFI,XCOLEXII                     ! Constants for pristine ice
                                                 ! self-collection
!
REAL,SAVE :: XAUTO3, XAUTO4,                   & ! Constants for pristine ice
             XLAUTS,   XLAUTS_THRESHOLD,       & ! autoconversion : AUT
             XITAUTS, XITAUTS_THRESHOLD,       &
             XTEXAUTI
!
REAL,SAVE :: XCOLEXIS,                         & ! Constants for snow 
    	     XAGGS_CLARGE1,XAGGS_CLARGE2,      & ! aggregation : AGG
             XAGGS_RLARGE1,XAGGS_RLARGE2
!
REAL,SAVE :: XHMTMIN,XHMTMAX,XHM1,XHM2,        & ! Constants for the
             XHM_YIELD,XHM_COLLCS,XHM_FACTS,   & !      revised
                       XHM_COLLCG,XHM_FACTG,   & ! Hallett-Mossop process
    	     XGAMINC_HMC_BOUND_MIN,            & ! Min val. of Lbda_c for HMC
    	     XGAMINC_HMC_BOUND_MAX,            & ! Max val. of Lbda_c for HMC
             XHMSINTP1,XHMSINTP2,              & ! (this is no more used !)
             XHMLINTP1,XHMLINTP2
!
REAL,SAVE :: XDCSLIM,XCOLCS,                   & ! Constants for the riming of
    	     XEXCRIMSS,XCRIMSS,                & ! the aggregates : RIM
    	     XEXCRIMSG,XCRIMSG,                & !
    	     XEXSRIMCG,XSRIMCG,                & !
    	     XGAMINC_BOUND_MIN,                & ! Min val. of Lbda_s for RIM
    	     XGAMINC_BOUND_MAX,                & ! Max val. of Lbda_s for RIM
    	     XRIMINTP1,XRIMINTP2                 ! Csts for lin. interpol. of 
                                                 ! the tab. incomplete Gamma law
INTEGER,SAVE :: NGAMINC                          ! Number of tab. Lbda_s
REAL, DIMENSION(:), SAVE, ALLOCATABLE          &
                       :: XGAMINC_RIM1,        & ! Tab. incomplete Gamma funct.
                          XGAMINC_RIM2,        & ! for XDS+2 and for XBS
                          XGAMINC_HMC            ! and for the HM process
!
REAL,SAVE :: XFRACCSS,                         & ! Constants for the accretion 
             XLBRACCS1,XLBRACCS2,XLBRACCS3,    & ! raindrops onto the aggregates
             XFSACCRG,                         & ! ACC (processes RACCSS and
             XLBSACCR1,XLBSACCR2,XLBSACCR3,    & !                SACCRG)
    	     XACCLBDAS_MIN,                    & ! Min val. of Lbda_s for ACC
    	     XACCLBDAS_MAX,                    & ! Max val. of Lbda_s for ACC
    	     XACCLBDAR_MIN,                    & ! Min val. of Lbda_r for ACC
    	     XACCLBDAR_MAX,                    & ! Max val. of Lbda_r for ACC
    	     XACCINTP1S,XACCINTP2S,            & ! Csts for bilin. interpol. of 
    	     XACCINTP1R,XACCINTP2R               !   Lbda_s and Lbda_r in the 
                        			 ! XKER_RACCSS and XKER_SACCRG
                          			 !            tables
INTEGER,SAVE :: NACCLBDAS,                     & ! Number of Lbda_s values and
    	        NACCLBDAR                        !   of Lbda_r values in the
                        			 ! XKER_RACCSS and XKER_SACCRG
                        			 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
       			 :: XKER_RACCSS,       & ! Normalized kernel for RACCSS
       			    XKER_RACCS,        & ! Normalized kernel for RACCS
       			    XKER_SACCRG          ! Normalized kernel for SACCRG 
REAL,SAVE :: XFSCVMG                             ! Melting-conversion factor of
                                                 ! the aggregates
!
REAL,SAVE :: XCOLIR,                           & ! Constants for rain contact
    	     XEXRCFRI,XRCFRI,                  & ! freezing : CFR
    	     XEXICFRR,XICFRR                     !
!
REAL,SAVE :: XFCDRYG,                          & ! Constants for the dry growth
             XCOLCG,                           & ! of the graupeln :
             XCOLIG,XCOLEXIG,XFIDRYG,          & ! 
             XCOLSG,XCOLEXSG,XFSDRYG,          & !             RCDRYG
             XLBSDRYG1,XLBSDRYG2,XLBSDRYG3,    & !             RIDRYG
             XFRDRYG,                          & !             RSDRYG
             XLBRDRYG1,XLBRDRYG2,XLBRDRYG3,    & !             RRDRYG
             XDRYLBDAR_MIN,                    & ! Min val. of Lbda_r for DRY
    	     XDRYLBDAR_MAX,                    & ! Max val. of Lbda_r for DRY
             XDRYLBDAS_MIN,                    & ! Min val. of Lbda_s for DRY
    	     XDRYLBDAS_MAX,                    & ! Max val. of Lbda_s for DRY
    	     XDRYLBDAG_MIN,                    & ! Min val. of Lbda_g for DRY
    	     XDRYLBDAG_MAX,                    & ! Max val. of Lbda_g for DRY
    	     XDRYINTP1R,XDRYINTP2R,            & ! Csts for bilin. interpol. of 
    	     XDRYINTP1S,XDRYINTP2S,            & ! Lbda_r, Lbda_s and Lbda_g in
    	     XDRYINTP1G,XDRYINTP2G               ! the XKER_SDRYG and XKER_RDRYG
                                                 !            tables
INTEGER,SAVE :: NDRYLBDAR,                     & ! Number of Lbda_r,
    	        NDRYLBDAS,                     & !        of Lbda_s and
    	        NDRYLBDAG                        !        of Lbda_g values in
    	                                         ! the XKER_SDRYG and XKER_RDRYG
                        			 !            tables
REAL,DIMENSION(:,:), SAVE, ALLOCATABLE         &
                         :: XKER_SDRYG,        & ! Normalized kernel for SDRYG
       			    XKER_RDRYG           ! Normalized kernel for RDRYG 
REAL,SAVE :: XCONCI_MAX                          ! Limitation of the pristine 
                                   ! ice concentration (init and grid-nesting) 
REAL,SAVE :: XFREFFI  ! Factor to compute the cloud ice effective radius
!
END MODULE MODD_ICE_C1R3_PARAM 
