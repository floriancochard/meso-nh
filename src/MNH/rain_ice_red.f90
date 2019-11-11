!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
       MODULE MODI_RAIN_ICE_RED
!      ########################
!
INTERFACE
      SUBROUTINE RAIN_ICE_RED ( OSEDIC,HSEDIM, HSUBG_AUCV_RC, OWARM, KKA, KKU, KKL,   &
                            PTSTEP, KRR, LDMICRO, PEXN,             &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC,PINPRR, PINPRR3D, PEVAP3D,           &
                            PINPRS, PINPRG, PSIGS, PINDEP, PRAINFR, PSEA, PTOWN,  &
                            PRHT, PRHS, PINPRH, PFPR                              )
!
!
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Switch for rc->rr Subgrid autoconversion
                                        ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
!
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index  
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step
                                                   ! (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL, DIMENSION(:,:,:), INTENT(IN)   :: LDMICRO ! mask to limit computation
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source

!
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(:,:,:),INTENT(OUT)      :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:,:), INTENT(OUT)     :: PRAINFR! Rain fraction            
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA ! Sea Mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN! Fraction that is town 
REAL, DIMENSION(:,:,:), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)      :: PINPRH! Hail instant precip
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)  :: PFPR ! upper-air precipitation fluxes
!
END SUBROUTINE RAIN_ICE_RED
END INTERFACE
END MODULE MODI_RAIN_ICE_RED
!     ######spl
      SUBROUTINE RAIN_ICE_RED ( OSEDIC,HSEDIM, HSUBG_AUCV_RC, OWARM, KKA, KKU, KKL,   &
                            PTSTEP, KRR, LDMICRO, PEXN,             &
                            PDZZ, PRHODJ, PRHODREF, PEXNREF, PPABST, PCIT, PCLDFR,&
                            PTHT, PRVT, PRCT, PRRT, PRIT, PRST,                   &
                            PRGT, PTHS, PRVS, PRCS, PRRS, PRIS, PRSS, PRGS,       &
                            PINPRC,PINPRR, PINPRR3D, PEVAP3D,           &
                            PINPRS, PINPRG, PSIGS, PINDEP, PRAINFR, PSEA, PTOWN,  &
                            PRHT, PRHS, PINPRH, PFPR                              )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the slow microphysical sources
!!    which can be computed explicitly
!!
!!
!!**  METHOD
!!    ------
!!      The autoconversion computation follows Kessler (1969).
!!      The sedimentation rate is computed with a time spliting technique and
!!    an upstream scheme, written as a difference of non-advective fluxes. This
!!    source term is added to the future instant ( split-implicit process ).
!!      The others microphysical processes are evaluated at the central instant
!!    (split-explicit process ): autoconversion, accretion and rain evaporation.
!!      These last 3 terms are bounded in order not to create negative values
!!    for the water species at the future instant.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!      Module MODD_CONF :
!!          CCONF configuration of the model for the first time step
!!      Module MODD_CST
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XMD,XMV            ! Molecular weight for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!          XCL                ! Cl (liquid)
!!          XCI                ! Ci (solid)
!!          XTT                ! Triple point temperature
!!          XLVTT              ! Vaporization heat constant
!!          XALPW,XBETAW,XGAMW ! Constants for saturation vapor pressure
!!                               function over liquid water
!!          XALPI,XBETAI,XGAMI ! Constants for saturation vapor pressure
!!                               function over solid ice
!!      Module MODD_BUDGET:
!!         NBUMOD       : model in which budget is calculated
!!         CBUTYPE      : type of desired budget
!!                          'CART' for cartesian box configuration
!!                          'MASK' for budget zone defined by a mask
!!                          'NONE'  ' for no budget
!!         NBUPROCCTR   : process counter used for each budget variable
!!         LBU_RTH      : logical for budget of RTH (potential temperature)
!!                        .TRUE. = budget of RTH
!!                        .FALSE. = no budget of RTH
!!         LBU_RRV      : logical for budget of RRV (water vapor)
!!                        .TRUE. = budget of RRV
!!                        .FALSE. = no budget of RRV
!!         LBU_RRC      : logical for budget of RRC (cloud water)
!!                        .TRUE. = budget of RRC
!!                        .FALSE. = no budget of RRC
!!         LBU_RRI      : logical for budget of RRI (cloud ice)
!!                        .TRUE. = budget of RRI
!!                        .FALSE. = no budget of RRI
!!         LBU_RRR      : logical for budget of RRR (rain water)
!!                        .TRUE. = budget of RRR
!!                        .FALSE. = no budget of RRR
!!         LBU_RRS      : logical for budget of RRS (aggregates)
!!                        .TRUE. = budget of RRS
!!                        .FALSE. = no budget of RRS
!!         LBU_RRG      : logical for budget of RRG (graupeln)
!!                        .TRUE. = budget of RRG
!!                        .FALSE. = no budget of RRG
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1 and Book2 of documentation ( routine RAIN_ICE )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/11/95
!!      (J.Viviand) 04/02/97  debug accumulated prcipitation & convert
!!                            precipitation rate in m/s
!!      (J.-P. Pinty) 17/02/97  add budget calls
!!      (J.-P. Pinty) 17/11/97  set ice sedim. for cirrus ice, reset RCHONI
!!                              and RRHONG, reverse order for DEALLOCATE
!!      (J.-P. Pinty) 11/02/98  correction of the air dynamical viscosity and
!!                              add advance of the budget calls
!!      (J.-P. Pinty) 18/05/98  correction of the air density in the RIAUTS
!!                              process
!!      (J.-P. Pinty) 18/11/98  split the main routine
!!      (V. Masson)   18/11/98  bug in IVEC1 and IVEC2 upper limits
!!      (J. Escobar & J.-P. Pinty)
!!                    11/12/98  contains and rewrite count+pack
!!      (J. Stein & J.-P. Pinty)
!!                    14/10/99  correction for very small RIT
!!      (J. Escobar & J.-P. Pinty)
!!                    24/07/00  correction for very samll m.r. in
!!                              the sedimentation subroutine
!!      (M. Tomasini) 11/05/01  Autoconversion of rc into rr modification to take
!!                              into account the subgrid variance
!!                              (cf Redelsperger & Sommeria JAS 86)
!!      (G. Molinie)  21/05/99  bug in RRCFRIG process, RHODREF**(-1) missing
!!                              in RSRIMCG
!!      (G. Molinie & J.-P. Pinty)
!!                    21/06/99  bug in RACCS process
!!      (P. Jabouille) 27/05/04 safety test for case where esw/i(T)> pabs (~Z>40km)
!!      (J-.P. Chaboureau) 12/02/05  temperature depending ice-to-snow autocon-
!                              version threshold (Chaboureau and Pinty GRL 2006)
!!      (J.-P. Pinty) 01/01/O1  add the hail category and correction of the
!!                              wet growth rate of the graupeln
!!      (S.Remy & C.Lac) 06/06 Add the cloud sedimentation
!!      (S.Remy & C.Lac) 06/06 Sedimentation becoming the last process
!!      to settle the precipitating species created during the current time step
!!      (S.Remy & C.Lac) 06/06 Modification of the algorithm of sedimentation
!!      to settle n times the precipitating species created during Dt/n instead
!!      of Dt
!!      (C.Lac) 11/06 Optimization of the sedimentation loop for NEC
!!      (J.Escobar) 18/01/2008 Parallel Bug in Budget when IMICRO >= 1
!!                  --> Path inhibit this test by IMICRO >= 0 allway true
!!      (Y.Seity) 03/2008 Add Statistic sedimentation
!!      (Y.Seity) 10/2009 Added condition for the raindrop accretion of the aggregates
!!         into graupeln process (5.2.6) to avoid negative graupel mixing ratio
!!      (V.Masson, C.Lac) 09/2010 Correction in split sedimentation for
!!                                reproducibility
!!      (S. Riette) Oct 2010 Better vectorisation of RAIN_ICE_SEDIMENTATION_STAT
!!      (Y. Seity), 02-2012  add possibility to run with reversed vertical levels
!!      (L. Bengtsson), 02-2013 Passing in land/sea mask and town fraction in
!!                      order to use different cloud droplet number conc. over
!!                      land, sea and urban areas in the cloud sedimentation.
!!      (D. Degrauwe), 2013-11: Export upper-air precipitation fluxes PFPR.
!!      (S. Riette) Nov 2013 Protection against null sigma
!!      Juan 24/09/2012: for BUG Pgi rewrite PACK function on mode_pack_pgi
!!      (C. Lac) FIT temporal scheme : instant M removed
!!      (JP Pinty), 01-2014 : ICE4 : partial reconversion of hail to graupel
!!              July, 2015 (O.Nuissier/F.Duffourg) Add microphysics diagnostic for
!!                                      aircraft, ballon and profiler
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      C.Lac : 10/2016 : add droplet deposition
!!      C.Lac : 01/2017 : correction on droplet deposition
!!      J.Escobar : 10/2017 : for real*4 , limit exp() in RAIN_ICE_SLOW with XMNH_HUGE_12_LOG
!!      (C. Abiven, Y. Léauté, V. Seigner, S. Riette) Phasing of Turner rain subgrid param
!!      (S. Riette) Source code split into several files
!!                  02/2019 C.Lac add rain fraction as an output field
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODD_PARAM_ICE
USE MODD_BUDGET
USE MODD_LES
USE MODI_BUDGET
USE MODI_ICE4_RAINFR_VERT
USE MODI_ICE4_SEDIMENTATION_STAT
USE MODI_ICE4_SEDIMENTATION_SPLIT
USE MODI_ICE4_NUCLEATION_WRAPPER
USE MODI_ICE4_TENDENCIES
USE MODE_FMWRIT
USE MODE_ll
USE MODE_MSG
!
#ifdef MNH_PGI
USE MODE_PACK_PGI
#endif
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
LOGICAL,                  INTENT(IN)    :: OSEDIC ! Switch for droplet sedim.
CHARACTER(LEN=4),         INTENT(IN)    :: HSEDIM ! Sedimentation scheme
CHARACTER(LEN=4),         INTENT(IN)    :: HSUBG_AUCV_RC ! Kind of Subgrid autoconversion method
LOGICAL,                  INTENT(IN)    :: OWARM   ! .TRUE. allows raindrops to
                                                   !   form by warm processes
                                                   !      (Kessler scheme)
INTEGER,                  INTENT(IN)    :: KKA   !near ground array index
INTEGER,                  INTENT(IN)    :: KKU   !uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL   !vert. levels type 1=MNH -1=ARO
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL, DIMENSION(:,:,:), INTENT(IN)   :: LDMICRO ! mask to limit computation
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ    ! Layer thikness (m)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS    ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRIS    ! Pristine ice m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRSS    ! Snow/aggregate m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRGS    ! Graupel m.r. source
!
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRC! Cloud instant precip
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINDEP  ! Cloud instant deposition
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(:,:,:),INTENT(OUT)      :: PINPRR3D! Rain inst precip 3D
REAL, DIMENSION(:,:,:), INTENT(OUT)     :: PEVAP3D! Rain evap profile
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(:,:), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(:,:,:), INTENT(OUT)     :: PRAINFR! Rain fraction            
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA ! Sea Mask
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN! Fraction that is town 
REAL, DIMENSION(:,:,:), OPTIONAL,  INTENT(IN)    :: PRHT    ! Hail m.r. at t
REAL, DIMENSION(:,:,:), OPTIONAL,  INTENT(INOUT) :: PRHS    ! Hail m.r. source
REAL, DIMENSION(:,:), OPTIONAL, INTENT(OUT)      :: PINPRH! Hail instant precip
REAL, DIMENSION(:,:,:,:), OPTIONAL, INTENT(OUT)  :: PFPR ! upper-air precipitation fluxes
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB           !  Define the domain where is
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IIT           !
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IJT           !
INTEGER :: IKB, IKTB, IKT!
INTEGER :: IKE, IKTE     !
!
!For packing
INTEGER :: IMICRO ! Case r_x>0 locations
INTEGER, DIMENSION(COUNT(LDMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                             :: JL       ! and PACK intrinsics
!
!Arrays for nucleation call outisde of LDMICRO points
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) :: ZW ! work array
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3)) :: ZT ! Temperature
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: &
                                  & ZZ_RVHENI_MR, & ! heterogeneous nucleation mixing ratio change
                                  & ZZ_RVHENI       ! heterogeneous nucleation
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: ZZ_LVFACT, ZZ_LSFACT
!
!Diagnostics
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) ::            &
                                                         & ZHLC_HCF3D,& ! HLCLOUDS cloud fraction in high water content part
                                                         & ZHLC_LCF3D,& ! HLCLOUDS cloud fraction in low water content part
                                                         & ZHLC_HRC3D,& ! HLCLOUDS cloud water content in high water content
                                                         & ZHLC_LRC3D   ! HLCLOUDS cloud water content in low water content
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2)) :: ZINPRI ! Pristine ice instant precip
!
!Packed variables
REAL, DIMENSION(COUNT(LDMICRO)) :: ZRVT,     & ! Water vapor m.r. at t
                                 & ZRCT,     & ! Cloud water m.r. at t
                                 & ZRRT,     & ! Rain water m.r. at t
                                 & ZRIT,     & ! Pristine ice m.r. at t
                                 & ZRST,     & ! Snow/aggregate m.r. at t
                                 & ZRGT,     & ! Graupel m.r. at t
                                 & ZRHT,     & ! Hail m.r. at t
                                 & ZCIT,     & ! Pristine ice conc. at t
                                 & ZTHT,     & ! Potential temperature
                                 & ZRHODREF, & ! RHO Dry REFerence
                                 & ZZT,      & ! Temperature
                                 & ZPRES,    & ! Pressure
                                 & ZEXN,     & ! EXNer Pressure
                                 & ZLSFACT,  & ! L_s/(Pi*C_ph)
                                 & ZLVFACT,  & ! L_v/(Pi*C_ph)
                                 & ZSIGMA_RC,& ! Standard deviation of rc at time t
                                 & ZCF,      & ! Cloud fraction
                                 & ZHLC_HCF, & ! HLCLOUDS : fraction of High Cloud Fraction in grid
                                 & ZHLC_LCF, & ! HLCLOUDS : fraction of Low  Cloud Fraction in grid
                                               !    note that ZCF = ZHLC_HCF + ZHLC_LCF
                                 & ZHLC_HRC, & ! HLCLOUDS : LWC that is High LWC in grid
                                 & ZHLC_LRC    ! HLCLOUDS : LWC that is Low  LWC in grid
                                               !    note that ZRC = ZHLC_HRC + ZHLC_LRC
!
!Output packed tendencies (for budgets only)
REAL, DIMENSION(COUNT(LDMICRO)) :: ZRVHENI_MR, & ! heterogeneous nucleation mixing ratio change
                                 & ZRCHONI, & ! Homogeneous nucleation
                                 & ZRRHONG_MR, & ! Spontaneous freezing mixing ratio change
                                 & ZRVDEPS, & ! Deposition on r_s,
                                 & ZRIAGGS, & ! Aggregation on r_s
                                 & ZRIAUTS, & ! Autoconversion of r_i for r_s production
                                 & ZRVDEPG, & ! Deposition on r_g
                                 & ZRCAUTR,  & ! Autoconversion of r_c for r_r production
                                 & ZRCACCR, & ! Accretion of r_c for r_r production
                                 & ZRREVAV, & ! Evaporation of r_r
                                 & ZRIMLTC_MR, & ! Cloud ice melting mixing ratio change
                                 & ZRCBERI, & ! Bergeron-Findeisen effect
                                 & ZRHMLTR, & ! Melting of the hailstones
                                 & ZRSMLTG, & ! Conversion-Melting of the aggregates
                                 & ZRCMLTSR, & ! Cloud droplet collection onto aggregates by positive temperature
                                 & ZRRACCSS, ZRRACCSG, ZRSACCRG, & ! Rain accretion onto the aggregates
                                 & ZRCRIMSS, ZRCRIMSG, ZRSRIMCG, ZRSRIMCG_MR, & ! Cloud droplet riming of the aggregates
                                 & ZRICFRRG, ZRRCFRIG, ZRICFRR, & ! Rain contact freezing
                                 & ZRCWETG, ZRIWETG, ZRRWETG, ZRSWETG, &  ! Graupel wet growth
                                 & ZRCDRYG, ZRIDRYG, ZRRDRYG, ZRSDRYG, &  ! Graupel dry growth
                                 & ZRWETGH, & ! Conversion of graupel into hail
                                 & ZRWETGH_MR, & ! Conversion of graupel into hail, mr change
                                 & ZRGMLTR, & ! Melting of the graupel
                                 & ZRCWETH, ZRIWETH, ZRSWETH, ZRGWETH, ZRRWETH, & ! Dry growth of hailstone
                                 & ZRCDRYH, ZRIDRYH, ZRSDRYH, ZRRDRYH, ZRGDRYH, & ! Wet growth of hailstone
                                 & ZRDRYHG    ! Conversion of hailstone into graupel
!
!Output packed total mixing ratio change (for budgets only)
REAL, DIMENSION(COUNT(LDMICRO)) :: ZTOT_RVHENI, & ! heterogeneous nucleation mixing ratio change
                                 & ZTOT_RCHONI, & ! Homogeneous nucleation
                                 & ZTOT_RRHONG, & ! Spontaneous freezing mixing ratio change
                                 & ZTOT_RVDEPS, & ! Deposition on r_s,
                                 & ZTOT_RIAGGS, & ! Aggregation on r_s
                                 & ZTOT_RIAUTS, & ! Autoconversion of r_i for r_s production
                                 & ZTOT_RVDEPG, & ! Deposition on r_g
                                 & ZTOT_RCAUTR,  & ! Autoconversion of r_c for r_r production
                                 & ZTOT_RCACCR, & ! Accretion of r_c for r_r production
                                 & ZTOT_RREVAV, & ! Evaporation of r_r
                                 & ZTOT_RCRIMSS, ZTOT_RCRIMSG, ZTOT_RSRIMCG, & ! Cloud droplet riming of the aggregates
                                 & ZTOT_RIMLTC, & ! Cloud ice melting mixing ratio change
                                 & ZTOT_RCBERI, & ! Bergeron-Findeisen effect
                                 & ZTOT_RHMLTR, & ! Melting of the hailstones
                                 & ZTOT_RSMLTG, & ! Conversion-Melting of the aggregates
                                 & ZTOT_RCMLTSR, & ! Cloud droplet collection onto aggregates by positive temperature
                                 & ZTOT_RRACCSS, ZTOT_RRACCSG, ZTOT_RSACCRG, & ! Rain accretion onto the aggregates
                                 & ZTOT_RICFRRG, ZTOT_RRCFRIG, ZTOT_RICFRR, & ! Rain contact freezing
                                 & ZTOT_RCWETG, ZTOT_RIWETG, ZTOT_RRWETG, ZTOT_RSWETG, &  ! Graupel wet growth
                                 & ZTOT_RCDRYG, ZTOT_RIDRYG, ZTOT_RRDRYG, ZTOT_RSDRYG, &  ! Graupel dry growth
                                 & ZTOT_RWETGH, & ! Conversion of graupel into hail
                                 & ZTOT_RGMLTR, & ! Melting of the graupel
                                 & ZTOT_RCWETH, ZTOT_RIWETH, ZTOT_RSWETH, ZTOT_RGWETH, ZTOT_RRWETH, & ! Dry growth of hailstone
                                 & ZTOT_RCDRYH, ZTOT_RIDRYH, ZTOT_RSDRYH, ZTOT_RRDRYH, ZTOT_RGDRYH, & ! Wet growth of hailstone
                                 & ZTOT_RDRYHG    ! Conversion of hailstone into graupel
!
!For time- or mixing-ratio- splitting
REAL, DIMENSION(COUNT(LDMICRO)) :: Z0RVT,     &   ! Water vapor m.r. at the beginig of the current loop
                                 & Z0RCT,     &   ! Cloud water m.r. at the beginig of the current loop
                                 & Z0RRT,     &   ! Rain water m.r. at the beginig of the current loop
                                 & Z0RIT,     &   ! Pristine ice m.r. at the beginig of the current loop
                                 & Z0RST,     &   ! Snow/aggregate m.r. at the beginig of the current loop
                                 & Z0RGT,     &   ! Graupel m.r. at the beginig of the current loop
                                 & Z0RHT,     &   ! Hail m.r. at the beginig of the current loop
                                 & ZA_TH, ZA_RV, ZA_RC, ZA_RR, ZA_RI, ZA_RS, ZA_RG, ZA_RH, &
                                 & ZB_TH, ZB_RV, ZB_RC, ZB_RR, ZB_RI, ZB_RS, ZB_RG, ZB_RH
!
!To take into acount external tendencies inside the splitting
REAL, DIMENSION(COUNT(LDMICRO)) :: ZEXT_RV,   &   ! External tendencie for rv
                                   ZEXT_RC,   &   ! External tendencie for rc
                                   ZEXT_RR,   &   ! External tendencie for rr
                                   ZEXT_RI,   &   ! External tendencie for ri
                                   ZEXT_RS,   &   ! External tendencie for rs
                                   ZEXT_RG,   &   ! External tendencie for rg
                                   ZEXT_RH,   &   ! External tendencie for rh
                                   ZEXT_TH,   &   ! External tendencie for th
                                   ZEXT_WW        ! Working array
LOGICAL :: LEXT_TEND
!
INTEGER, DIMENSION(COUNT(LDMICRO)) :: IITER ! Number of iterations done (with real tendencies computation)
INTEGER :: INB_ITER_MAX ! Maximum number of iterations (with real tendencies computation)
REAL, DIMENSION(COUNT(LDMICRO)) :: ZTIME,    & ! Current integration time (starts with 0 and ends with PTSTEP)
                                 & ZMAXTIME, & ! Time on which we can apply the current tendencies
                                 & ZTIME_THRESHOLD, & ! Time to reach threshold
                                 & ZTIME_LASTCALL     ! Integration time when last tendecies call has been done
LOGICAL, DIMENSION(COUNT(LDMICRO)) :: LLCOMPUTE ! Points where we must compute tendenceis
LOGICAL :: LSOFT ! Must we really compute tendencies or only adjust them to new T variables
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2)):: GDEP
REAL :: ZTSTEP ! length of sub-timestep in case of time splitting
REAL :: ZINV_TSTEP ! Inverse ov PTSTEP
REAL, DIMENSION(COUNT(LDMICRO), 6) :: ZRS_TEND
REAL, DIMENSION(COUNT(LDMICRO), 6) :: ZRG_TEND
REAL, DIMENSION(COUNT(LDMICRO), 8) :: ZRH_TEND
!
!For total tendencies computation
REAL, DIMENSION(SIZE(PTHT,1),SIZE(PTHT,2),SIZE(PTHT,3)) :: &
        &ZW_RVS, ZW_RCS, ZW_RRS, ZW_RIS, ZW_RSS, ZW_RGS, ZW_RHS, ZW_THS
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IIT=SIZE(PDZZ,1)
IJT=SIZE(PDZZ,2)
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL
IKT=SIZE(PDZZ,3)
IKTB=1+JPVEXT
IKTE=IKT-JPVEXT
!
ZINV_TSTEP=1./PTSTEP
LEXT_TEND=.TRUE.
!
ZT(:,:,:) = PTHT(:,:,:) * PEXN(:,:,:)
! LSFACT and LVFACT without exner
IF(KRR==7) THEN
  ZZ_LSFACT(:,:,:)=(XLSTT+(XCPV-XCI)*(ZT(:,:,:)-XTT))   &
                   /( XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
                   + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)+PRHT(:,:,:)))
  ZZ_LVFACT(:,:,:)=(XLVTT+(XCPV-XCL)*(ZT(:,:,:)-XTT))   &
                   /( XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
                   + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)+PRHT(:,:,:)))
ELSE
  ZZ_LSFACT(:,:,:)=(XLSTT+(XCPV-XCI)*(ZT(:,:,:)-XTT))   &
                   /( XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
                   + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)))
  ZZ_LVFACT(:,:,:)=(XLVTT+(XCPV-XCL)*(ZT(:,:,:)-XTT))   &
                   /( XCPD + XCPV*PRVT(:,:,:) + XCL*(PRCT(:,:,:)+PRRT(:,:,:))   &
                   + XCI*(PRIT(:,:,:)+PRST(:,:,:)+PRGT(:,:,:)))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(.NOT. LSEDIM_AFTER) THEN
  !
  !*       2.1     sedimentation
  !
  IF(HSEDIM=='STAT') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA, PTOWN,  &
                                  &PINPRH=PINPRH, PRHT=PRHS*PTSTEP, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC,  PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA, PTOWN,  &
                                  &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !No negativity correction here as we apply sedimentation on PR.S*PTSTEP variables
  ELSEIF(HSEDIM=='SPLI') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC,XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN,  &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC,XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN, &
                                   &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !We correct negativities with conservation
    !SPLI algorith uses a time-splitting. Inside the loop a temporary m.r. is used.
    !   It is initialized with the m.r. at T and is modified by two tendencies:
    !   sedimentation tendency and an external tendency which represents all other
    !   processes (mainly advection and microphysical processes). If both tendencies
    !   are negative, sedimentation can remove a specie at a given sub-timestep. From
    !   this point sedimentation stops for the remaining sub-timesteps but the other tendency
    !   will be still active and will lead to negative values.
    !   We could prevent the algorithm to not consume too much a specie, instead we apply
    !   a correction here.
    CALL CORRECT_NEGATIVITIES(KRR, PRVS, PRCS, PRRS, &
                             &PRIS, PRSS, PRGS, &
                             &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)
  ELSE
    WRITE(*,*) ' STOP'
    WRITE(*,*) ' NO SEDIMENTATION SCHEME FOR HSEDIM=', HSEDIM
    CALL PRINT_MSG(NVERB_FATAL,'GEN','RAIN_ICE_RED','')
  END IF
  !
  !*       2.2     budget storage
  !
  IF (LBUDGET_RC .AND. OSEDIC) &
                  CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:), 7 , 'SEDI_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:), 8 , 'SEDI_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET (PRIS(:,:,:)*PRHODJ(:,:,:), 9 , 'SEDI_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET (PRSS(:,:,:)*PRHODJ(:,:,:), 10, 'SEDI_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET (PRGS(:,:,:)*PRHODJ(:,:,:), 11, 'SEDI_BU_RRG')
  IF ( KRR == 7 .AND. LBUDGET_RH) &
                  CALL BUDGET (PRHS(:,:,:)*PRHODJ(:,:,:), 12, 'SEDI_BU_RRH')
  IF ( LBUDGET_RC .AND. LDEPOSC ) &
   CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:),7 ,'DEPO_BU_RRC')
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     PACKING
!               --------
!  optimization by looking for locations where
!  the microphysical fields are larger than a minimal value only !!!
!
IMICRO=0
IF(COUNT(LDMICRO)/=0) IMICRO=RAIN_ICE_COUNTJV(LDMICRO(:,:,:), IIT, IJT, IKT, SIZE(I1), I1(:), I2(:), I3(:))
!Packing
IF(IMICRO>0) THEN
  DO JL=1, IMICRO
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
    ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
    ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
    ZCF(JL) = PCLDFR(I1(JL),I2(JL),I3(JL))
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZTHT(JL) = PTHT(I1(JL),I2(JL),I3(JL))
    ZPRES(JL) = PPABST(I1(JL),I2(JL),I3(JL))
    ZEXN(JL) = PEXN(I1(JL),I2(JL),I3(JL))
  ENDDO
  IF(LEXT_TEND) THEN
    DO JL=1, IMICRO
      ZEXT_RV(JL) = PRVS(I1(JL),I2(JL),I3(JL)) - ZRVT(JL)*ZINV_TSTEP
      ZEXT_RC(JL) = PRCS(I1(JL),I2(JL),I3(JL)) - ZRCT(JL)*ZINV_TSTEP
      ZEXT_RR(JL) = PRRS(I1(JL),I2(JL),I3(JL)) - ZRRT(JL)*ZINV_TSTEP
      ZEXT_RI(JL) = PRIS(I1(JL),I2(JL),I3(JL)) - ZRIT(JL)*ZINV_TSTEP
      ZEXT_RS(JL) = PRSS(I1(JL),I2(JL),I3(JL)) - ZRST(JL)*ZINV_TSTEP
      ZEXT_RG(JL) = PRGS(I1(JL),I2(JL),I3(JL)) - ZRGT(JL)*ZINV_TSTEP
      ZEXT_TH(JL) = PTHS(I1(JL),I2(JL),I3(JL)) - ZTHT(JL)*ZINV_TSTEP
      !The th tendency is not related to a mixing ratio change, there is no exn/exnref issue here
    ENDDO
  ENDIF
  IF(HSUBG_AUCV_RC=='PDF ' .AND. CSUBG_PR_PDF=='SIGM') THEN
    DO JL=1, IMICRO
      ZSIGMA_RC(JL) = PSIGS(I1(JL),I2(JL),I3(JL))*2.
!     ZSIGMA_RC(JL) = MAX(PSIGS(I1(JL),I2(JL),I3(JL)) * 2., 1.E-12)
    ENDDO
  ENDIF
  IF(KRR==7) THEN
    DO JL=1, IMICRO
      ZRHT(JL) = PRHT(I1(JL),I2(JL),I3(JL))
    ENDDO
    IF(LEXT_TEND) THEN
      DO JL=1, IMICRO
        ZEXT_RH(JL) = PRHS(I1(JL),I2(JL),I3(JL)) - ZRHT(JL)*ZINV_TSTEP
      ENDDO
    ENDIF
  ELSE
    ZRHT(:)=0.
    IF(LEXT_TEND) ZEXT_RH(:)=0.
  ENDIF
  IF(LBU_ENABLE) THEN
    ZTOT_RVHENI(:)=0.
    ZTOT_RCHONI(:)=0.
    ZTOT_RRHONG(:)=0.
    ZTOT_RVDEPS(:)=0.
    ZTOT_RIAGGS(:)=0.
    ZTOT_RIAUTS(:)=0.
    ZTOT_RVDEPG(:)=0.
    ZTOT_RCAUTR(:)=0.
    ZTOT_RCACCR(:)=0.
    ZTOT_RREVAV(:)=0.
    ZTOT_RCRIMSS(:)=0.
    ZTOT_RCRIMSG(:)=0.
    ZTOT_RSRIMCG(:)=0.
    ZTOT_RIMLTC(:)=0.
    ZTOT_RCBERI(:)=0.
    ZTOT_RHMLTR(:)=0.
    ZTOT_RSMLTG(:)=0.
    ZTOT_RCMLTSR(:)=0.
    ZTOT_RRACCSS(:)=0.
    ZTOT_RRACCSG(:)=0.
    ZTOT_RSACCRG(:)=0.
    ZTOT_RICFRRG(:)=0.
    ZTOT_RRCFRIG(:)=0.
    ZTOT_RICFRR(:)=0.
    ZTOT_RCWETG(:)=0.
    ZTOT_RIWETG(:)=0.
    ZTOT_RRWETG(:)=0.
    ZTOT_RSWETG(:)=0.
    ZTOT_RCDRYG(:)=0.
    ZTOT_RIDRYG(:)=0.
    ZTOT_RRDRYG(:)=0.
    ZTOT_RSDRYG(:)=0.
    ZTOT_RWETGH(:)=0.
    ZTOT_RGMLTR(:)=0.
    ZTOT_RCWETH(:)=0.
    ZTOT_RIWETH(:)=0.
    ZTOT_RSWETH(:)=0.
    ZTOT_RGWETH(:)=0.
    ZTOT_RRWETH(:)=0.
    ZTOT_RCDRYH(:)=0.
    ZTOT_RIDRYH(:)=0.
    ZTOT_RSDRYH(:)=0.
    ZTOT_RRDRYH(:)=0.
    ZTOT_RGDRYH(:)=0.
    ZTOT_RDRYHG(:)=0.
  ENDIF
ENDIF
!-------------------------------------------------------------------------------
!
!*       4.     LOOP
!               ----
!
!Maximum number of iterations
!We only count real iterations (those for which we *compute* tendencies)
INB_ITER_MAX=NMAXITER
IF(XTSTEP_TS/=0.)THEN
  INB_ITER_MAX=MAX(1, INT(PTSTEP/XTSTEP_TS)) !At least the number of iterations needed for the time-splitting
  ZTSTEP=PTSTEP/INB_ITER_MAX
  INB_ITER_MAX=MAX(NMAXITER, INB_ITER_MAX) !Fot the case XMRSTEP/=0. at the same time
ENDIF
IITER(:)=0
ZTIME(:)=0. ! Current integration time (all points may have a different integration time)
DO WHILE(ANY(ZTIME(:)<PTSTEP)) ! Loop to *really* compute tendencies
  IF(XMRSTEP/=0.) THEN
    ! In this case we need to remember the mixing ratios used to compute the tendencies
    ! because when mixing ratio has evolved more than a threshold, we must re-compute tendecies
    Z0RVT(:)=ZRVT(:)
    Z0RCT(:)=ZRCT(:)
    Z0RRT(:)=ZRRT(:)
    Z0RIT(:)=ZRIT(:)
    Z0RST(:)=ZRST(:)
    Z0RGT(:)=ZRGT(:)
    Z0RHT(:)=ZRHT(:)
  ENDIF
  IF(XTSTEP_TS/=0.) THEN
    ! In this case we need to remember the time when tendencies were computed
    ! because when time has evolved more than a limit, we must re-compute tendecies
    ZTIME_LASTCALL(:)=ZTIME(:)
  ENDIF
  LLCOMPUTE(:)=ZTIME(:)<PTSTEP ! Compuation only for points for which integration time has not reached the timestep
  LSOFT=.FALSE. ! We *really* compute the tendencies
  WHERE(LLCOMPUTE(:))
    IITER(:)=IITER(:)+1
  END WHERE
  DO WHILE(ANY(LLCOMPUTE(:))) ! Loop to adjust tendencies when we cross the 0°C or when a specie disappears
    ZZT(:) = ZTHT(:) * ZEXN(:)
    IF(KRR==7) THEN
      ZLSFACT(:)=(XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))   &
                 /( (XCPD + XCPV*ZRVT(:) + XCL*(ZRCT(:)+ZRRT(:))   &
                 + XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)+ZRHT(:)))*ZEXN(:) )
      ZLVFACT(:)=(XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))   &
                 /( (XCPD + XCPV*ZRVT(:) + XCL*(ZRCT(:)+ZRRT(:))   &
                 + XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)+ZRHT(:)))*ZEXN(:) )
    ELSE
      ZLSFACT(:)=(XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))   &
                 /( (XCPD + XCPV*ZRVT(:) + XCL*(ZRCT(:)+ZRRT(:))   &
                 + XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)))*ZEXN(:) )
      ZLVFACT(:)=(XLVTT+(XCPV-XCL)*(ZZT(:)-XTT))   &
                 /( (XCPD + XCPV*ZRVT(:) + XCL*(ZRCT(:)+ZRRT(:))   &
                 + XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)))*ZEXN(:) )
    ENDIF
    !
    !***       4.1 Tendecies computation
    !
    ! Tendencies are *really* computed when LSOFT==.FALSE. and only adjusted otherwise
    CALL ICE4_TENDENCIES(IMICRO, IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKT, KKL, &
                        &KRR, LSOFT, LLCOMPUTE, &
                        &OWARM, CSUBG_RC_RR_ACCR, CSUBG_RR_EVAP, HSUBG_AUCV_RC, CSUBG_PR_PDF, &
                        &ZEXN, ZRHODREF, ZLVFACT, ZLSFACT, LDMICRO, I1, I2, I3, &
                        &ZPRES, ZCF, ZSIGMA_RC, &
                        &ZCIT, &
                        &ZZT, ZTHT, &
                        &ZRVT, ZRCT, ZRRT, ZRIT, ZRST, ZRGT, ZRHT, PRRT, &
                        &ZRVHENI_MR, ZRRHONG_MR, ZRIMLTC_MR, ZRSRIMCG_MR, &
                        &ZRCHONI, ZRVDEPS, ZRIAGGS, ZRIAUTS, ZRVDEPG, &
                        &ZRCAUTR, ZRCACCR, ZRREVAV, &
                        &ZRCRIMSS, ZRCRIMSG, ZRSRIMCG, ZRRACCSS, ZRRACCSG, ZRSACCRG, ZRSMLTG, ZRCMLTSR, &
                        &ZRICFRRG, ZRRCFRIG, ZRICFRR, ZRCWETG, ZRIWETG, ZRRWETG, ZRSWETG, &
                        &ZRCDRYG, ZRIDRYG, ZRRDRYG, ZRSDRYG, ZRWETGH, ZRWETGH_MR, ZRGMLTR, &
                        &ZRCWETH, ZRIWETH, ZRSWETH, ZRGWETH, ZRRWETH, &
                        &ZRCDRYH, ZRIDRYH, ZRSDRYH, ZRRDRYH, ZRGDRYH, ZRDRYHG, ZRHMLTR, &
                        &ZRCBERI, &
                        &ZRS_TEND, ZRG_TEND, ZRH_TEND, &
                        &ZA_TH, ZA_RV, ZA_RC, ZA_RR, ZA_RI, ZA_RS, ZA_RG, ZA_RH, &
                        &ZB_TH, ZB_RV, ZB_RC, ZB_RR, ZB_RI, ZB_RS, ZB_RG, ZB_RH, &
                        &ZHLC_HCF, ZHLC_LCF, ZHLC_HRC, ZHLC_LRC, PRAINFR)
    ! External tendencies
    IF(LEXT_TEND) THEN
      ZA_TH(:) = ZA_TH(:) + ZEXT_TH(:)
      ZA_RV(:) = ZA_RV(:) + ZEXT_RV(:)
      ZA_RC(:) = ZA_RC(:) + ZEXT_RC(:)
      ZA_RR(:) = ZA_RR(:) + ZEXT_RR(:)
      ZA_RI(:) = ZA_RI(:) + ZEXT_RI(:)
      ZA_RS(:) = ZA_RS(:) + ZEXT_RS(:)
      ZA_RG(:) = ZA_RG(:) + ZEXT_RG(:)
      ZA_RH(:) = ZA_RH(:) + ZEXT_RH(:)
    ENDIF
    !
    !***       4.2 Integration time
    !
    ! If we can, we will use these tendecies until the end of the timestep
    ZMAXTIME(:)=0.
    WHERE(LLCOMPUTE(:))
      ZMAXTIME(:)=PTSTEP-ZTIME(:) ! Remaining time until the end of the timestep
    ENDWHERE

    !We need to adjust tendencies when temperature reaches 0
    IF(LFEEDBACKT) THEN
      !Is ZB_TH enough to change temperature sign?
      WHERE( (ZTHT(:) - XTT/ZEXN(:)) * (ZTHT(:) + ZB_TH(:) - XTT/ZEXN(:)) < 0. )
        ZMAXTIME(:)=0.
      ENDWHERE
      !Can ZA_TH make temperature change of sign?
      ZTIME_THRESHOLD(:)=-1.
      WHERE(ABS(ZA_TH(:))>1.E-20)
        ZTIME_THRESHOLD(:)=(XTT/ZEXN(:) - ZB_TH(:) - ZTHT(:))/ZA_TH(:)
      ENDWHERE
      WHERE(ZTIME_THRESHOLD(:)>0.)
        ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
      ENDWHERE
    ENDIF

    !We need to adjust tendencies when a specy disappears
    !When a specy is missing, only the external tendencies can be negative (and we must keep track of it)
    WHERE(ZA_RV(:)<-1.E-20 .AND. ZRVT(:)>XRTMIN(1))
      ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RV(:)+ZRVT(:))/ZA_RV(:))
    END WHERE
    WHERE(ZA_RC(:)<-1.E-20 .AND. ZRCT(:)>XRTMIN(2))
      ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RC(:)+ZRCT(:))/ZA_RC(:))
    END WHERE
    WHERE(ZA_RR(:)<-1.E-20 .AND. ZRRT(:)>XRTMIN(3))
      ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RR(:)+ZRRT(:))/ZA_RR(:))
    END WHERE
    WHERE(ZA_RI(:)<-1.E-20 .AND. ZRIT(:)>XRTMIN(4))
      ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RI(:)+ZRIT(:))/ZA_RI(:))
    END WHERE
    WHERE(ZA_RS(:)<-1.E-20 .AND. ZRST(:)>XRTMIN(5))
      ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RS(:)+ZRST(:))/ZA_RS(:))
    END WHERE
    WHERE(ZA_RG(:)<-1.E-20 .AND. ZRGT(:)>XRTMIN(6))
      ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RG(:)+ZRGT(:))/ZA_RG(:))
    END WHERE
    IF(KRR==7) THEN
      WHERE(ZA_RH(:)<-1.E-20 .AND. ZRHT(:)>XRTMIN(7))
        ZMAXTIME(:)=MIN(ZMAXTIME(:), -(ZB_RH(:)+ZRHT(:))/ZA_RH(:))
      END WHERE
    ENDIF

    !We stop when the end of the timestep is reached
    WHERE(PTSTEP-ZTIME(:)-ZMAXTIME(:)<=0.)
      LLCOMPUTE(:)=.FALSE.
    ENDWHERE

    !We must recompute tendencies when the end of the sub-timestep is reached
    IF(XTSTEP_TS/=0.) THEN
      WHERE(IITER(:)<INB_ITER_MAX .AND. ZTIME(:)+ZMAXTIME(:)>ZTIME_LASTCALL(:)+ZTSTEP)
        ZMAXTIME(:)=ZTIME_LASTCALL(:)-ZTIME(:)+ZTSTEP
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE
    ENDIF

    !We must recompute tendencies when the maximum allowed change is reached
    !When a specy is missing, only the external tendencies can be active and we do not want to recompute
    !the microphysical tendencies when external tendencies are negative (results won't change because specy was already missing)
    IF(XMRSTEP/=0.) THEN
      ZTIME_THRESHOLD(:)=-1.
      WHERE(IITER(:)<INB_ITER_MAX .AND. ABS(ZA_RV(:))>1.E-20)
        ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RV(:))*XMRSTEP+Z0RVT(:)-ZRVT(:)-ZB_RV(:))/ZA_RV(:)
      ENDWHERE
      WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
           &(ZRVT(:)>XRTMIN(1) .OR. ZA_RV(:)>0.))
        ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE

      ZTIME_THRESHOLD(:)=-1.
      WHERE(IITER(:)<INB_ITER_MAX .AND. ABS(ZA_RC(:))>1.E-20)
        ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RC(:))*XMRSTEP+Z0RCT(:)-ZRCT(:)-ZB_RC(:))/ZA_RC(:)
      ENDWHERE
      WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
           &(ZRCT(:)>XRTMIN(2) .OR. ZA_RC(:)>0.))
        ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE

      ZTIME_THRESHOLD(:)=-1.
      WHERE(IITER(:)<INB_ITER_MAX .AND. ABS(ZA_RR(:))>1.E-20)
        ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RR(:))*XMRSTEP+Z0RRT(:)-ZRRT(:)-ZB_RR(:))/ZA_RR(:)
      ENDWHERE
      WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
           &(ZRRT(:)>XRTMIN(3) .OR. ZA_RR(:)>0.))
        ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE

      ZTIME_THRESHOLD(:)=-1.
      WHERE(IITER(:)<INB_ITER_MAX .AND. ABS(ZA_RI(:))>1.E-20)
        ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RI(:))*XMRSTEP+Z0RIT(:)-ZRIT(:)-ZB_RI(:))/ZA_RI(:)
      ENDWHERE
      WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
           &(ZRIT(:)>XRTMIN(4) .OR. ZA_RI(:)>0.))
        ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE

      ZTIME_THRESHOLD(:)=-1.
      WHERE(IITER(:)<INB_ITER_MAX .AND. ABS(ZA_RS(:))>1.E-20)
        ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RS(:))*XMRSTEP+Z0RST(:)-ZRST(:)-ZB_RS(:))/ZA_RS(:)
      ENDWHERE
      WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
           &(ZRST(:)>XRTMIN(5) .OR. ZA_RS(:)>0.))
        ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE

      ZTIME_THRESHOLD(:)=-1.
      WHERE(IITER(:)<INB_ITER_MAX .AND. ABS(ZA_RG(:))>1.E-20)
        ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RG(:))*XMRSTEP+Z0RGT(:)-ZRGT(:)-ZB_RG(:))/ZA_RG(:)
      ENDWHERE
      WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
           &(ZRGT(:)>XRTMIN(6) .OR. ZA_RG(:)>0.))
        ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE

      IF(KRR==7) THEN
        ZTIME_THRESHOLD(:)=-1.
        WHERE(IITER(:)<INB_ITER_MAX .AND. ABS(ZA_RH(:))>1.E-20)
          ZTIME_THRESHOLD(:)=(SIGN(1., ZA_RH(:))*XMRSTEP+Z0RHT(:)-ZRHT(:)-ZB_RH(:))/ZA_RH(:)
        ENDWHERE
        WHERE(ZTIME_THRESHOLD(:)>=0. .AND. ZTIME_THRESHOLD(:)<ZMAXTIME(:) .AND. &
             &(ZRHT(:)>XRTMIN(7) .OR. ZA_RH(:)>0.))
          ZMAXTIME(:)=MIN(ZMAXTIME(:), ZTIME_THRESHOLD(:))
          LLCOMPUTE(:)=.FALSE.
        ENDWHERE
      ENDIF

      WHERE(IITER(:)<INB_ITER_MAX .AND. MAX(ABS(ZB_RV(:)), ABS(ZB_RC(:)), ABS(ZB_RR(:)), ABS(ZB_RI(:)), &
                                            ABS(ZB_RS(:)), ABS(ZB_RG(:)), ABS(ZB_RH(:)))>XMRSTEP)
        ZMAXTIME(:)=0.
        LLCOMPUTE(:)=.FALSE.
      ENDWHERE
    ENDIF
    !
    !***       4.3 New values of variables for next iteration
    !
    ZTHT=ZTHT+ZA_TH(:)*ZMAXTIME(:)+ZB_TH(:)
    ZRVT=ZRVT+ZA_RV(:)*ZMAXTIME(:)+ZB_RV(:)
    ZRCT=ZRCT+ZA_RC(:)*ZMAXTIME(:)+ZB_RC(:)
    ZRRT=ZRRT+ZA_RR(:)*ZMAXTIME(:)+ZB_RR(:)
    ZRIT=ZRIT+ZA_RI(:)*ZMAXTIME(:)+ZB_RI(:)
    ZRST=ZRST+ZA_RS(:)*ZMAXTIME(:)+ZB_RS(:)
    ZRGT=ZRGT+ZA_RG(:)*ZMAXTIME(:)+ZB_RG(:)
    IF(KRR==7) ZRHT=ZRHT+ZA_RH(:)*ZMAXTIME(:)+ZB_RH(:)
    WHERE(ZRIT(:)==0.)
      ZCIT(:) = 0.
    END WHERE
    !
    !***       4.4 Mixing ratio change due to each process
    !
    IF(LBU_ENABLE) THEN
      ZTOT_RVHENI(:)= ZTOT_RVHENI(:) +ZRVHENI_MR(:)
      ZTOT_RCHONI(:)= ZTOT_RCHONI(:) +ZRCHONI(:) *ZMAXTIME(:)
      ZTOT_RRHONG(:)= ZTOT_RRHONG(:) +ZRRHONG_MR(:)
      ZTOT_RVDEPS(:)= ZTOT_RVDEPS(:) +ZRVDEPS(:) *ZMAXTIME(:)
      ZTOT_RIAGGS(:)= ZTOT_RIAGGS(:) +ZRIAGGS(:) *ZMAXTIME(:)
      ZTOT_RIAUTS(:)= ZTOT_RIAUTS(:) +ZRIAUTS(:) *ZMAXTIME(:)
      ZTOT_RVDEPG(:)= ZTOT_RVDEPG(:) +ZRVDEPG(:) *ZMAXTIME(:)
      ZTOT_RCAUTR(:)= ZTOT_RCAUTR(:) +ZRCAUTR(:) *ZMAXTIME(:)
      ZTOT_RCACCR(:)= ZTOT_RCACCR(:) +ZRCACCR(:) *ZMAXTIME(:)
      ZTOT_RREVAV(:)= ZTOT_RREVAV(:) +ZRREVAV(:) *ZMAXTIME(:)
      ZTOT_RCRIMSS(:)=ZTOT_RCRIMSS(:)+ZRCRIMSS(:)*ZMAXTIME(:)
      ZTOT_RCRIMSG(:)=ZTOT_RCRIMSG(:)+ZRCRIMSG(:)*ZMAXTIME(:)
      ZTOT_RSRIMCG(:)=ZTOT_RSRIMCG(:)+ZRSRIMCG(:)*ZMAXTIME(:)+ZRSRIMCG_MR(:)
      ZTOT_RRACCSS(:)=ZTOT_RRACCSS(:)+ZRRACCSS(:)*ZMAXTIME(:)
      ZTOT_RRACCSG(:)=ZTOT_RRACCSG(:)+ZRRACCSG(:)*ZMAXTIME(:)
      ZTOT_RSACCRG(:)=ZTOT_RSACCRG(:)+ZRSACCRG(:)*ZMAXTIME(:)
      ZTOT_RSMLTG(:)= ZTOT_RSMLTG(:) +ZRSMLTG(:) *ZMAXTIME(:)
      ZTOT_RCMLTSR(:)=ZTOT_RCMLTSR(:)+ZRCMLTSR(:) *ZMAXTIME(:)
      ZTOT_RICFRRG(:)=ZTOT_RICFRRG(:)+ZRICFRRG(:)*ZMAXTIME(:)
      ZTOT_RRCFRIG(:)=ZTOT_RRCFRIG(:)+ZRRCFRIG(:)*ZMAXTIME(:)
      ZTOT_RICFRR(:)= ZTOT_RICFRR(:) +ZRICFRR(:) *ZMAXTIME(:)
      ZTOT_RCWETG(:)= ZTOT_RCWETG(:) +ZRCWETG(:) *ZMAXTIME(:)
      ZTOT_RIWETG(:)= ZTOT_RIWETG(:) +ZRIWETG(:) *ZMAXTIME(:)
      ZTOT_RRWETG(:)= ZTOT_RRWETG(:) +ZRRWETG(:) *ZMAXTIME(:)
      ZTOT_RSWETG(:)= ZTOT_RSWETG(:) +ZRSWETG(:) *ZMAXTIME(:)
      ZTOT_RWETGH(:)= ZTOT_RWETGH(:) +ZRWETGH(:) *ZMAXTIME(:)+ZRWETGH_MR(:)
      ZTOT_RCDRYG(:)= ZTOT_RCDRYG(:) +ZRCDRYG(:) *ZMAXTIME(:)
      ZTOT_RIDRYG(:)= ZTOT_RIDRYG(:) +ZRIDRYG(:) *ZMAXTIME(:)
      ZTOT_RRDRYG(:)= ZTOT_RRDRYG(:) +ZRRDRYG(:) *ZMAXTIME(:)
      ZTOT_RSDRYG(:)= ZTOT_RSDRYG(:) +ZRSDRYG(:) *ZMAXTIME(:)
      ZTOT_RGMLTR(:)= ZTOT_RGMLTR(:) +ZRGMLTR(:) *ZMAXTIME(:)
      ZTOT_RCWETH(:)= ZTOT_RCWETH(:) +ZRCWETH(:) *ZMAXTIME(:)
      ZTOT_RIWETH(:)= ZTOT_RIWETH(:) +ZRIWETH(:) *ZMAXTIME(:)
      ZTOT_RSWETH(:)= ZTOT_RSWETH(:) +ZRSWETH(:) *ZMAXTIME(:)
      ZTOT_RGWETH(:)= ZTOT_RGWETH(:) +ZRGWETH(:) *ZMAXTIME(:)
      ZTOT_RRWETH(:)= ZTOT_RRWETH(:) +ZRRWETH(:) *ZMAXTIME(:)
      ZTOT_RCDRYH(:)= ZTOT_RCDRYH(:) +ZRCDRYH(:) *ZMAXTIME(:)
      ZTOT_RIDRYH(:)= ZTOT_RIDRYH(:) +ZRIDRYH(:) *ZMAXTIME(:)
      ZTOT_RSDRYH(:)= ZTOT_RSDRYH(:) +ZRSDRYH(:) *ZMAXTIME(:)
      ZTOT_RRDRYH(:)= ZTOT_RRDRYH(:) +ZRRDRYH(:) *ZMAXTIME(:)
      ZTOT_RGDRYH(:)= ZTOT_RGDRYH(:) +ZRGDRYH(:) *ZMAXTIME(:)
      ZTOT_RDRYHG(:)= ZTOT_RDRYHG(:) +ZRDRYHG(:) *ZMAXTIME(:)
      ZTOT_RHMLTR(:)= ZTOT_RHMLTR(:) +ZRHMLTR(:) *ZMAXTIME(:)
      ZTOT_RIMLTC(:)= ZTOT_RIMLTC(:) +ZRIMLTC_MR(:)
      ZTOT_RCBERI(:)= ZTOT_RCBERI(:) +ZRCBERI(:) *ZMAXTIME(:)
    ENDIF
    !
    !***       4.5 Next loop
    !
    LSOFT=.TRUE. ! We try to adjust tendencies (inner while loop)
    ZTIME(:)=ZTIME(:)+ZMAXTIME(:)
  ENDDO
ENDDO









!-------------------------------------------------------------------------------
!
!*       5.     UNPACKING DIAGNOSTICS
!               ---------------------
!
IF(IMICRO>0) THEN
  ZW(:,:,:) = 0.
  ZHLC_HCF3D(:,:,:) = UNPACK(ZHLC_HCF(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))
  ZW(:,:,:) = 0.
  ZHLC_LCF3D(:,:,:) = UNPACK(ZHLC_LCF(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))
  ZW(:,:,:) = 0.
  ZHLC_HRC3D(:,:,:) = UNPACK(ZHLC_HRC(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))
  ZW(:,:,:) = 0.
  ZHLC_LRC3D(:,:,:) = UNPACK(ZHLC_LRC(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))
  PCIT(:,:,:) = UNPACK(ZCIT(:), MASK=LDMICRO(:,:,:), FIELD=PCIT(:,:,:))
ELSE
  PRAINFR(:,:,:)=0.
  ZHLC_HCF3D(:,:,:)=0.
  ZHLC_LCF3D(:,:,:)=0.
  ZHLC_HRC3D(:,:,:)=0.
  ZHLC_LRC3D(:,:,:)=0.
  PCIT(:,:,:) = 0.
ENDIF
IF(OWARM) THEN
  ZW(:,:,:)=0.
  PEVAP3D(:,:,:)=UNPACK(ZRREVAV(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))
ENDIF
!
!
!*       6.     COMPUTES THE SLOW COLD PROCESS SOURCES OUTSIDE OF LDMICRO POINTS
!               ----------------------------------------------------------------
!
CALL ICE4_NUCLEATION_WRAPPER(IIT, IJT, IKT, .NOT. LDMICRO, &
                             PTHT, PPABST, PRHODREF, PEXN, ZZ_LSFACT/PEXN, ZT, &
                             PRVT, &
                             PCIT, ZZ_RVHENI_MR)
ZZ_LSFACT(:,:,:)=ZZ_LSFACT(:,:,:)/PEXNREF(:,:,:)
ZZ_LVFACT(:,:,:)=ZZ_LVFACT(:,:,:)/PEXNREF(:,:,:)
ZZ_RVHENI(:,:,:) = MIN(PRVS(:,:,:), ZZ_RVHENI_MR(:,:,:)/PTSTEP)
PRIS(:,:,:)=PRIS(:,:,:)+ZZ_RVHENI(:,:,:)
PRVS(:,:,:)=PRVS(:,:,:)-ZZ_RVHENI(:,:,:)
PTHS(:,:,:)=PTHS(:,:,:) + ZZ_RVHENI(:,:,:)*ZZ_LSFACT(:,:,:)
!-------------------------------------------------------------------------------
!
!*       7.     UNPACKING AND TOTAL TENDENCIES
!               ------------------------------
!
!
!***     7.1    total tendencies limited by available species
!
! ZW_??S variables will contain the new S variables values
!
IF(LEXT_TEND) THEN
  !Z..T variables contain the exeternal tendency, we substract it
  ZRVT(:) = ZRVT(:) - ZEXT_RV(:) * PTSTEP
  ZRCT(:) = ZRCT(:) - ZEXT_RC(:) * PTSTEP
  ZRRT(:) = ZRRT(:) - ZEXT_RR(:) * PTSTEP
  ZRIT(:) = ZRIT(:) - ZEXT_RI(:) * PTSTEP
  ZRST(:) = ZRST(:) - ZEXT_RS(:) * PTSTEP
  ZRGT(:) = ZRGT(:) - ZEXT_RG(:) * PTSTEP
  IF (KRR==7) ZRHT(:) = ZRHT(:) - ZEXT_RH(:) * PTSTEP
  ZTHT(:) = ZTHT(:) - ZEXT_TH(:) * PTSTEP
ENDIF
!Tendencies computed from difference between old state and new state (can be negative)
ZW_RVS(:,:,:) = (UNPACK(ZRVT(:), MASK=LDMICRO(:,:,:), FIELD=PRVT(:,:,:)) - PRVT(:,:,:))*ZINV_TSTEP
ZW_RCS(:,:,:) = (UNPACK(ZRCT(:), MASK=LDMICRO(:,:,:), FIELD=PRCT(:,:,:)) - PRCT(:,:,:))*ZINV_TSTEP
ZW_RRS(:,:,:) = (UNPACK(ZRRT(:), MASK=LDMICRO(:,:,:), FIELD=PRRT(:,:,:)) - PRRT(:,:,:))*ZINV_TSTEP
ZW_RIS(:,:,:) = (UNPACK(ZRIT(:), MASK=LDMICRO(:,:,:), FIELD=PRIT(:,:,:)) - PRIT(:,:,:))*ZINV_TSTEP
ZW_RSS(:,:,:) = (UNPACK(ZRST(:), MASK=LDMICRO(:,:,:), FIELD=PRST(:,:,:)) - PRST(:,:,:))*ZINV_TSTEP
ZW_RGS(:,:,:) = (UNPACK(ZRGT(:), MASK=LDMICRO(:,:,:), FIELD=PRGT(:,:,:)) - PRGT(:,:,:))*ZINV_TSTEP
IF(KRR==7) THEN
  ZW_RHS(:,:,:) = (UNPACK(ZRHT(:), MASK=LDMICRO(:,:,:), FIELD=PRHT(:,:,:)) - PRHT(:,:,:))*ZINV_TSTEP
ELSE
  ZW_RHS(:,:,:) = 0.
ENDIF
ZW_THS(:,:,:) = (ZW_RCS(:,:,:)+ZW_RRS(:,:,:)                            )*ZZ_LVFACT(:,:,:) + &
              & (ZW_RIS(:,:,:)+ZW_RSS(:,:,:)+ZW_RGS(:,:,:)+ZW_RHS(:,:,:))*ZZ_LSFACT(:,:,:)
!We apply these tendencies to the S variables
ZW_RVS(:,:,:) = PRVS(:,:,:) + ZW_RVS(:,:,:)
ZW_RCS(:,:,:) = PRCS(:,:,:) + ZW_RCS(:,:,:)
ZW_RRS(:,:,:) = PRRS(:,:,:) + ZW_RRS(:,:,:)
ZW_RIS(:,:,:) = PRIS(:,:,:) + ZW_RIS(:,:,:)
ZW_RSS(:,:,:) = PRSS(:,:,:) + ZW_RSS(:,:,:)
ZW_RGS(:,:,:) = PRGS(:,:,:) + ZW_RGS(:,:,:)
IF(KRR==7) ZW_RHS(:,:,:) = PRHS(:,:,:) + ZW_RHS(:,:,:)
ZW_THS(:,:,:) = PTHS(:,:,:) + ZW_THS(:,:,:)
!We correct negativities with conservation
CALL CORRECT_NEGATIVITIES(KRR, ZW_RVS, ZW_RCS, ZW_RRS, &
                         &ZW_RIS, ZW_RSS, ZW_RGS, &
                         &ZW_THS, ZZ_LVFACT, ZZ_LSFACT, ZW_RHS)
!
!***     7.2    LBU_ENABLE case
!
IF(LBU_ENABLE) THEN
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RVHENI(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*ZZ_LSFACT(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'HENU_BU_RTH')
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'HENU_BU_RRV')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'HENU_BU_RRI')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCHONI(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'HON_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'HON_BU_RRC')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'HON_BU_RRI')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRHONG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'SFR_BU_RTH')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'SFR_BU_RRR')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'SFR_BU_RRG')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RVDEPS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*ZZ_LSFACT(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DEPS_BU_RTH')
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'DEPS_BU_RRV')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'DEPS_BU_RRS')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIAGGS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'AGGS_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'AGGS_BU_RRS')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIAUTS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'AUTS_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'AUTS_BU_RRS')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RVDEPG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PRVS(:,:,:) = PRVS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*ZZ_LSFACT(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DEPG_BU_RTH')
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'DEPG_BU_RRV')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'DEPG_BU_RRG')

  IF(OWARM) THEN
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCAUTR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'AUTO_BU_RRC')
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'AUTO_BU_RRR')

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCACCR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'ACCR_BU_RRC')
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'ACCR_BU_RRR')

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RREVAV(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
    PRVS(:,:,:) = PRVS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*ZZ_LVFACT(:,:,:)
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'REVA_BU_RTH')
    IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'REVA_BU_RRV')
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'REVA_BU_RRR')
  ENDIF

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCRIMSS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCRIMSG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSRIMCG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'RIM_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'RIM_BU_RRC')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'RIM_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'RIM_BU_RRG')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRACCSS(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRSS(:,:,:) = PRSS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRACCSG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSACCRG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'ACC_BU_RTH')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'ACC_BU_RRR')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'ACC_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'ACC_BU_RRG')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSMLTG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCMLTSR, MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'CMEL_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'CMEL_BU_RRG')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'CMEL_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'CMEL_BU_RRR')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RICFRRG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRCFRIG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RICFRR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'CFRZ_BU_RTH')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'CFRZ_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'CFRZ_BU_RRI')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'CFRZ_BU_RRG')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSWETG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'WETG_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'WETG_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'WETG_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'WETG_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'WETG_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'WETG_BU_RRG')

  IF(KRR==7) THEN
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RWETGH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'GHCV_BU_RRG')
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'GHCV_BU_RRH')
  END IF

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RRDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RSDRYG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DRYG_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'DRYG_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'DRYG_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'DRYG_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'DRYG_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'DRYG_BU_RRG')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RGMLTR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
  PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'GMLT_BU_RTH')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'GMLT_BU_RRR')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'GMLT_BU_RRG')

  IF(KRR==7) THEN
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RRWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RIWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RSWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RGWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'WETH_BU_RTH')
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'WETH_BU_RRC')
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'WETH_BU_RRR')
    IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'WETH_BU_RRI')
    IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'WETH_BU_RRS')
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'WETH_BU_RRH')

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RGWETH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'HGCV_BU_RRG')
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'HGCV_BU_RRH')

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RCDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RRDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RIDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RSDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRSS(:,:,:) = PRSS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RGDRYH(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRGS(:,:,:) = PRGS(:,:,:) - ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) + ZW(:,:,:)
    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RDRYHG(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRHS(:,:,:) = PRHS(:,:,:) - ZW(:,:,:)
    PRGS(:,:,:) = PRGS(:,:,:) + ZW(:,:,:)
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'DRYH_BU_RTH')
    IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'DRYH_BU_RRC')
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'DRYH_BU_RRR')
    IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'DRYH_BU_RRI')
    IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'DRYH_BU_RRS')
    IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'DRYH_BU_RRG')
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'DRYH_BU_RRH')

    ZW(:,:,:) = 0.
    ZW(:,:,:)=UNPACK(ZTOT_RHMLTR(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
    PRRS(:,:,:) = PRRS(:,:,:) + ZW(:,:,:)
    PRHS(:,:,:) = PRHS(:,:,:) - ZW(:,:,:)
    PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
    IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'HMLT_BU_RTH')
    IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'HMLT_BU_RRR')
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'HMLT_BU_RRH')
  ENDIF

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RIMLTC(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRIS(:,:,:) = PRIS(:,:,:) - ZW(:,:,:)
  PRCS(:,:,:) = PRCS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) - ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'IMLT_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'IMLT_BU_RRC')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'IMLT_BU_RRI')

  ZW(:,:,:) = 0.
  ZW(:,:,:)=UNPACK(ZTOT_RCBERI(:), MASK=LDMICRO(:,:,:), FIELD=ZW(:,:,:))*ZINV_TSTEP
  PRCS(:,:,:) = PRCS(:,:,:) - ZW(:,:,:)
  PRIS(:,:,:) = PRIS(:,:,:) + ZW(:,:,:)
  PTHS(:,:,:) = PTHS(:,:,:) + ZW(:,:,:)*(ZZ_LSFACT(:,:,:)-ZZ_LVFACT(:,:,:))
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'BERFI_BU_RTH')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'BERFI_BU_RRC')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'BERFI_BU_RRI')
ENDIF
!
!***     7.3    Final tendencies
!
PRVS(:,:,:) = ZW_RVS(:,:,:)
PRCS(:,:,:) = ZW_RCS(:,:,:)
PRRS(:,:,:) = ZW_RRS(:,:,:)
PRIS(:,:,:) = ZW_RIS(:,:,:)
PRSS(:,:,:) = ZW_RSS(:,:,:)
PRGS(:,:,:) = ZW_RGS(:,:,:)
IF (KRR==7) THEN
  PRHS(:,:,:) = ZW_RHS(:,:,:)
ENDIF
PTHS(:,:,:) = ZW_THS(:,:,:)
IF(LBU_ENABLE) THEN
  IF (LBUDGET_TH) CALL BUDGET(PTHS(:,:,:)*PRHODJ(:,:,:), 4, 'CORR_BU_RTH')
  IF (LBUDGET_RV) CALL BUDGET(PRVS(:,:,:)*PRHODJ(:,:,:), 6, 'CORR_BU_RRV')
  IF (LBUDGET_RC) CALL BUDGET(PRCS(:,:,:)*PRHODJ(:,:,:), 7, 'CORR_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET(PRRS(:,:,:)*PRHODJ(:,:,:), 8, 'CORR_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET(PRIS(:,:,:)*PRHODJ(:,:,:), 9, 'CORR_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET(PRSS(:,:,:)*PRHODJ(:,:,:), 10,'CORR_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET(PRGS(:,:,:)*PRHODJ(:,:,:), 11,'CORR_BU_RRG')
  IF (KRR==7) THEN
    IF (LBUDGET_RH) CALL BUDGET(PRHS(:,:,:)*PRHODJ(:,:,:), 12,'CORR_BU_RRH')
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       8.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!               -------------------------------------
!
IF(LSEDIM_AFTER) THEN
  !
  !*       8.1     sedimentation
  !
  IF(HSEDIM=='STAT') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA, PTOWN,  &
                                  &PINPRH=PINPRH, PRHT=PRHS*PTSTEP, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_STAT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                  &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ,&
                                  &PRHODREF, PPABST, PTHT, PRHODJ, &
                                  &PRCS, PRCS*PTSTEP, PRRS, PRRS*PTSTEP, PRIS, PRIS*PTSTEP,&
                                  &PRSS, PRSS*PTSTEP, PRGS, PRGS*PTSTEP,&
                                  &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                  &PSEA, PTOWN,  &
                                  &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !No negativity correction here as we apply sedimentation on PR.S*PTSTEP variables
  ELSEIF(HSEDIM=='SPLI') THEN
    !SR: It *seems* that we must have two separate calls for ifort
    IF(KRR==7) THEN
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN,  &
                                   &PINPRH=PINPRH, PRHT=PRHT, PRHS=PRHS, PFPR=PFPR)
    ELSE
      CALL ICE4_SEDIMENTATION_SPLIT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKTB, IKTE, IKT, KKL, &
                                   &PTSTEP, KRR, OSEDIC, LDEPOSC, XVDEPOSC, PDZZ, &
                                   &PRHODREF, PPABST, PTHT, PRHODJ, &
                                   &PRCS, PRCT, PRRS, PRRT, PRIS, PRIT, PRSS, PRST, PRGS, PRGT,&
                                   &PINPRC, PINDEP, PINPRR, ZINPRI, PINPRS, PINPRG, &
                                   &PSEA, PTOWN,  &
                                   &PFPR=PFPR)
    ENDIF
    PINPRS(:,:) = PINPRS(:,:) + ZINPRI(:,:)
    !We correct negativities with conservation
    !SPLI algorith uses a time-splitting. Inside the loop a temporary m.r. is used.
    !   It is initialized with the m.r. at T and is modified by two tendencies:
    !   sedimentation tendency and an external tendency which represents all other
    !   processes (mainly advection and microphysical processes). If both tendencies
    !   are negative, sedimentation can remove a specie at a given sub-timestep. From
    !   this point sedimentation stops for the remaining sub-timesteps but the other tendency
    !   will be still active and will lead to negative values.
    !   We could prevent the algorithm to not consume too much a specie, instead we apply
    !   a correction here.
    CALL CORRECT_NEGATIVITIES(KRR, PRVS, PRCS, PRRS, &
                             &PRIS, PRSS, PRGS, &
                             &PTHS, ZZ_LVFACT, ZZ_LSFACT, PRHS)
  ELSE
    WRITE(*,*) ' STOP'
    WRITE(*,*) ' NO SEDIMENTATION SCHEME FOR HSEDIM=', HSEDIM
    CALL PRINT_MSG(NVERB_FATAL,'GEN','RAIN_ICE_RED','')
  END IF
  !
  !*       8.2     budget storage
  !
  IF (LBUDGET_RC .AND. OSEDIC) &
                  CALL BUDGET (PRCS(:,:,:)*PRHODJ(:,:,:), 7 , 'SEDI_BU_RRC')
  IF (LBUDGET_RR) CALL BUDGET (PRRS(:,:,:)*PRHODJ(:,:,:), 8 , 'SEDI_BU_RRR')
  IF (LBUDGET_RI) CALL BUDGET (PRIS(:,:,:)*PRHODJ(:,:,:), 9 , 'SEDI_BU_RRI')
  IF (LBUDGET_RS) CALL BUDGET (PRSS(:,:,:)*PRHODJ(:,:,:), 10, 'SEDI_BU_RRS')
  IF (LBUDGET_RG) CALL BUDGET (PRGS(:,:,:)*PRHODJ(:,:,:), 11, 'SEDI_BU_RRG')
  IF ( KRR == 7 .AND. LBUDGET_RH) &
                  CALL BUDGET (PRHS(:,:,:)*PRHODJ(:,:,:), 12, 'SEDI_BU_RRH')
  !
  !sedimentation of rain fraction
  CALL ICE4_RAINFR_VERT(IIB, IIE, IIT, IJB, IJE, IJT, IKB, IKE, IKT, KKL, PRAINFR, PRRS(:,:,:)*PTSTEP)

ENDIF
!
!
CONTAINS
  FUNCTION RAIN_ICE_COUNTJV(LTAB, KIT, KJT, KKT, KSIZE, I1,I2,I3) RESULT(IC)
  !
  !*      0. DECLARATIONS
  !          ------------
  !
  IMPLICIT NONE
  !
  !*       0.2  declaration of local variables
  !
  !
  INTEGER, INTENT(IN) :: KIT, KJT, KKT, KSIZE
  LOGICAL, DIMENSION(KIT,KJT,KKT), INTENT(IN) :: LTAB ! Mask
  INTEGER, DIMENSION(KSIZE), INTENT(OUT) :: I1,I2,I3 ! Used to replace the COUNT and PACK
  INTEGER :: JI,JJ,JK,IC
  !
  !-------------------------------------------------------------------------------
  !
  IC = 0
  DO JK = 1, SIZE(LTAB,3)
    DO JJ = 1, SIZE(LTAB,2)
      DO JI = 1, SIZE(LTAB,1)
        IF(LTAB(JI,JJ,JK)) THEN
          IC = IC +1
          I1(IC) = JI
          I2(IC) = JJ
          I3(IC) = JK
        END IF
      END DO
    END DO
  END DO
  !
  !
  END FUNCTION RAIN_ICE_COUNTJV
  !
  !
  SUBROUTINE CORRECT_NEGATIVITIES(KRR, PRV, PRC, PRR, &
                                 &PRI, PRS, PRG, &
                                 &PTH, PLVFACT, PLSFACT, PRH)
  !
  IMPLICIT NONE
  !
  INTEGER,                INTENT(IN)    :: KRR
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRV, PRC, PRR, PRI, PRS, PRG, PTH
  REAL, DIMENSION(:,:,:), INTENT(IN)    :: PLVFACT, PLSFACT
  REAL, DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PRH
  !
  REAL, DIMENSION(SIZE(PRV,1), SIZE(PRV,2), SIZE(PRV,3)) :: ZW
  !
  !We correct negativities with conservation
  ! 1) deal with negative values for mixing ratio, except for vapor
  WHERE(PRC(:,:,:)<0.)
    PRV(:,:,:)=PRV(:,:,:)+PRC(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-PRC(:,:,:)*PLVFACT(:,:,:)
    PRC(:,:,:)=0.
  ENDWHERE
  WHERE(PRR(:,:,:)<0.)
    PRV(:,:,:)=PRV(:,:,:)+PRR(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-PRR(:,:,:)*PLVFACT(:,:,:)
    PRR(:,:,:)=0.
  ENDWHERE
  WHERE(PRI(:,:,:)<0.)
    PRV(:,:,:)=PRV(:,:,:)+PRI(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-PRI(:,:,:)*PLSFACT(:,:,:)
    PRI(:,:,:)=0.
  ENDWHERE
  WHERE(PRS(:,:,:)<0.)
    PRV(:,:,:)=PRV(:,:,:)+PRS(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-PRS(:,:,:)*PLSFACT(:,:,:)
    PRS(:,:,:)=0.
  ENDWHERE
  WHERE(PRG(:,:,:)<0.)
    PRV(:,:,:)=PRV(:,:,:)+PRG(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-PRG(:,:,:)*PLSFACT(:,:,:)
    PRG(:,:,:)=0.
  ENDWHERE
  IF(KRR==7) THEN
    WHERE(PRH(:,:,:)<0.)
      PRV(:,:,:)=PRV(:,:,:)+PRH(:,:,:)
      PTH(:,:,:)=PTH(:,:,:)-PRH(:,:,:)*PLSFACT(:,:,:)
      PRH(:,:,:)=0.
    ENDWHERE
  ENDIF
  ! 2) deal with negative vapor mixing ratio
  WHERE(PRV(:,:,:)<0. .AND. PRC(:,:,:)+PRI(:,:,:)>0.)
    ! for rc and ri, we keep ice fraction constant
    ZW(:,:,:)=MIN(1., -PRV(:,:,:)/(PRC(:,:,:)+PRI(:,:,:))) ! Proportion of rc+ri to convert into rv
    PTH(:,:,:)=PTH(:,:,:)-ZW(:,:,:)*(PRC(:,:,:)*PLVFACT(:,:,:)+PRI(:,:,:)*PLSFACT(:,:,:))
    PRV(:,:,:)=PRV(:,:,:)+ZW(:,:,:)*(PRC(:,:,:)+PRI(:,:,:))
    PRC(:,:,:)=(1.-ZW(:,:,:))*PRC(:,:,:)
    PRI(:,:,:)=(1.-ZW(:,:,:))*PRI(:,:,:)
  ENDWHERE
  WHERE(PRV(:,:,:)<0. .AND. PRR(:,:,:)>0.)
    ZW(:,:,:)=MIN(PRR(:,:,:), -PRV(:,:,:)) ! Quantity of rr to convert into rv
    PRV(:,:,:)=PRV(:,:,:)+ZW(:,:,:)
    PRR(:,:,:)=PRR(:,:,:)-ZW(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-ZW(:,:,:)*PLVFACT(:,:,:)
  ENDWHERE
  WHERE(PRV(:,:,:)<0. .AND. PRS(:,:,:)>0.)
    ZW(:,:,:)=MIN(PRS(:,:,:), -PRV(:,:,:)) ! Quantity of rs to convert into rv
    PRV(:,:,:)=PRV(:,:,:)+ZW(:,:,:)
    PRS(:,:,:)=PRS(:,:,:)-ZW(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-ZW(:,:,:)*PLSFACT(:,:,:)
  ENDWHERE
  WHERE(PRV(:,:,:)<0. .AND. PRG(:,:,:)>0.)
    ZW(:,:,:)=MIN(PRG(:,:,:), -PRV(:,:,:)) ! Quantity of rg to convert into rv
    PRV(:,:,:)=PRV(:,:,:)+ZW(:,:,:)
    PRG(:,:,:)=PRG(:,:,:)-ZW(:,:,:)
    PTH(:,:,:)=PTH(:,:,:)-ZW(:,:,:)*PLSFACT(:,:,:)
  ENDWHERE
  IF(KRR==7) THEN
    WHERE(PRV(:,:,:)<0. .AND. PRH(:,:,:)>0.)
      ZW(:,:,:)=MIN(PRH(:,:,:), -PRV(:,:,:)) ! Quantity of rh to convert into rv
      PRV(:,:,:)=PRV(:,:,:)+ZW(:,:,:)
      PRH(:,:,:)=PRH(:,:,:)-ZW(:,:,:)
      PTH(:,:,:)=PTH(:,:,:)-ZW(:,:,:)*PLSFACT(:,:,:)
    ENDWHERE
  ENDIF
  !
  !
  END SUBROUTINE CORRECT_NEGATIVITIES
!
  recursive subroutine quicksort(a, first, last)
  implicit none
  integer  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
  end subroutine quicksort

END SUBROUTINE RAIN_ICE_RED
