!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/07/12 18:36:21
!-----------------------------------------------------------------
!     ############################
      MODULE MODI_CH_CONVECT_LINOX
!     ############################
!
INTERFACE
!
      SUBROUTINE CH_CONVECT_LINOX( KLON, KLEV, PCH1, PCH1C,            &
                                   KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                   PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                   PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                   KFTSTEPS, PUTT, PRHODREF,           &
                                   OUSECHEM, PZZ, PIC_RATE, PCG_RATE   )
!
INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
!
REAL,DIMENSION(KLON,KLEV),INTENT(IN)   :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(KLON,KLEV),INTENT(OUT)  :: PCH1C! conv adjusted tracer concntr.
!
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)
!
REAL, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PUTT      ! updraft temperature (K)
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PRHODREF  
!
LOGICAL,                      INTENT(IN) :: OUSECHEM ! to indicate if chemistry is used
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL, DIMENSION(KLON),     INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(KLON),     INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
!
END SUBROUTINE CH_CONVECT_LINOX
!
END INTERFACE
!
END MODULE MODI_CH_CONVECT_LINOX
!
!    ###################################################################
      SUBROUTINE CH_CONVECT_LINOX( KLON, KLEV, PCH1, PCH1C,            &
                                   KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                   PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                   PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                   KFTSTEPS, PUTT, PRHODREF,           &
                                   OUSECHEM, PZZ, PIC_RATE, PCG_RATE   )
!     ##################################################################
!
!!**** Compute the production of NOx by lightning flashes inside deep convective
!!     clouds and its transport
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!      environmental values of NOx=NO+NO2 produced by lightning flashes.
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PCH1C-PCH1)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      main deep convection code
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    AUTHOR
!!    ------
!!      T. Fehr       * Laboratoire d'Aerologie *
!!      J.-P. Pinty   * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/07/03
!! 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_CONVPAREXT
USE MODD_NSV,     ONLY : NSV_CHEMBEG, NSV_CHEMEND
!
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
!
REAL,DIMENSION(KLON,KLEV),INTENT(IN)   :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(KLON,KLEV),INTENT(OUT)  :: PCH1C! conv adjusted tracer concntr.
!
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)
!
REAL, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PUTT      ! updraft temperature (K)
REAL, DIMENSION(KLON,KLEV),   INTENT(IN) :: PRHODREF  
!
LOGICAL,                      INTENT(IN) :: OUSECHEM ! to indicate if chemistry is used
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL, DIMENSION(KLON),     INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(KLON),     INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIB,IIE,IJB,IJE,IKB,IKE  
INTEGER :: IKS            ! vertical dimension
INTEGER :: JI,JJ          ! horizontal loop index
INTEGER :: JK, JKP, JKM   ! vertical loop index
INTEGER :: JSTEP          ! fractional time loop index
INTEGER :: JKLC, JKLD, JKLT, JKLB, JKLF, JKLFF, JKLP, JKMAX ! loop index 
                                                            ! of levels
INTEGER, DIMENSION(KLON) :: ICFL  ! index freezing level
INTEGER, DIMENSION(KLON) :: ICFFL ! index of -10 degree C level
!
REAL, DIMENSION(KLON,KLEV)     :: ZOMG ! compensat. subsidence (Pa/s)
REAL, DIMENSION(KLON,KLEV,2)     :: ZUCH1, ZDCH1 ! updraft/downdraft values
REAL, DIMENSION(KLON,KLEV,2)     :: ZCH1C
REAL, DIMENSION(KLON)          :: ZTIMEC  ! fractional convective time step
REAL, DIMENSION(KLON,KLEV)     :: ZTIMC! 2D work array for ZTIMEC
REAL, DIMENSION(KLON,KLEV)     :: ZCH1MFIN, ZCH1MFOUT
                                   ! work arrays for environm. compensat. mass
REAL, DIMENSION(KLON)          :: ZWORK1, ZWORK2, ZWORK3
!
REAL, DIMENSION(KLON)          :: ZCLOUD_DEPTH,       &!
                                  ZTFLASH_RATE,       &!
                                  ZCLOUD_ABOVE_FREEZ, &!
                                  ZBETA,              &!
                                  ZSUM_LMASS_IC,      &!
                                  ZSUM_LMASS_CG
REAL, DIMENSION(KLON,KLEV)     :: ZNOX_PROD_IC, ZNOX_PROD_CG, ZCH1 
REAL                           :: XFRAC_UIC, XFRAC_UCG
!
! -----------------------------------------------------------------------------
!
!*       0.1   Compute loop bounds
!              -------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB    = 1 + JCVEXB 
IKS    = KLEV
IKE    = KLEV - JCVEXT 
JKMAX  = MAXVAL( KCTL(:) )
!
!*      1.      Lightning characteristics
!               -------------------------
!
XFRAC_UIC = 1.0 ! Fractional part of IC in updrafts
XFRAC_UCG = 1.0 ! Fractional part of CG in updrafts
!
!*      1.1     Cloud depth
!               -----------
!
DO JI = 1, KLON
  JKLT = KCTL(JI) ! cloud top
  ZCLOUD_DEPTH(JI) = PZZ(JI,JKLT) - PZZ(JI,IKB)
END DO
!
! Lightning frequency for continental clouds: Price and Rind (1992)
!
ZTFLASH_RATE(:) = 0.97241*EXP( (2.76183/XRADIUS**2)*PDXDY(:) ) & ! scale factor
                       * 3.44E-5 * (ZCLOUD_DEPTH(:)*1.0E-3)**4.9
ZTFLASH_RATE(:) = ZTFLASH_RATE(:)/60. 
!
!*      1.2     IC/CG lightning ratio
!               ---------------------
!
ICFL(:) = KLCL(:)
DO JI = 1, KLON
  DO JK = MINVAL( KLCL(:) ), JKMAX
    IF( PUTT(JI,JK)<XTT ) CYCLE       ! T < 273 K  (bottom level for ICs)
    ICFL(JI) = JK 
  END DO
END DO
!
ICFFL(:) = KCTL(:)
DO JI = 1, KLON
  DO JK = JKMAX,  MINVAL( KLCL(:) ), -1
    IF( PUTT(JI,JK)>XTT-10.0 ) CYCLE  ! T > 263 K (top level for CGs)
    ICFFL(JI) = JK 
  END DO
END DO
!
DO JI = 1, KLON
  JKLT = KCTL(JI) ! cloud top
  JKLF = ICFL(JI) ! cloud freezing level
  ZCLOUD_ABOVE_FREEZ(JI) = (PZZ(JI,JKLT) - PZZ(JI,JKLF))*1.E-3 ! in km
END DO
!
ZBETA(:) = (((0.021*ZCLOUD_ABOVE_FREEZ(:) -  0.648) &
                   *ZCLOUD_ABOVE_FREEZ(:) +  7.493) &
                   *ZCLOUD_ABOVE_FREEZ(:) - 36.540) &
                   *ZCLOUD_ABOVE_FREEZ(:) + 63.090
ZBETA(:) = MIN( 48.7, MAX( 0.19,ZBETA(:) ) ) ! 5.5km < ZCLOUD_ABOVE_FREEZ < 14km
!
!*      1.3     Profiles of NOx production rates
!               --------------------------------
!
ZNOX_PROD_IC(:,:) = 0.0
ZNOX_PROD_CG(:,:) = 0.0
ZSUM_LMASS_IC(:) = 0.0
ZSUM_LMASS_CG(:) = 0.0
ZCH1(:,:) = PCH1(:,:) ! initialize with the original value
DO JI = 1, KLON
  JKLT = KCTL(JI) ! cloud top
  JKLF  = ICFL(JI)  ! cloud freezing level
  JKLFF = ICFFL(JI) ! cloud -10 C level
  IF( JKLT<=JKLF ) THEN
    PIC_RATE(JI) = 0.0
    PCG_RATE(JI) = 0.0
    CYCLE
  ENDIF
!
! IC_NOx production rate
!
  PIC_RATE(JI) = (ZBETA(JI) / (1.0+ZBETA(JI)))*ZTFLASH_RATE(JI)
  DO JK = JKLF, JKLT
    ZSUM_LMASS_IC(JI)   = ZSUM_LMASS_IC(JI) + PLMASS(JI,JK)
    ZNOX_PROD_IC(JI,JK) = PIC_RATE(JI) * 6.7E25 * PLMASS(JI,JK)
  END DO
  ZNOX_PROD_IC(JI,:) = ZNOX_PROD_IC(JI,:) / ZSUM_LMASS_IC(JI)
!
! CG_NOx production rate
!
  PCG_RATE(JI) = (1.0 / (1.0+ZBETA(JI)))*ZTFLASH_RATE(JI)
  DO JK = IKB, JKLFF
    ZSUM_LMASS_CG(JI)   = ZSUM_LMASS_CG(JI) + PLMASS(JI,JK)
    ZNOX_PROD_CG(JI,JK) = PCG_RATE(JI) * 6.7E26 * PLMASS(JI,JK)
  END DO
  ZNOX_PROD_CG(JI,:) = ZNOX_PROD_CG(JI,:) / ZSUM_LMASS_CG(JI)
END DO
!
!*      1.4     Update units (molecules/s => pp/s) and integrate
!               ------------------------------------------------
!
ZNOX_PROD_IC(:,IKB:IKE) = ZNOX_PROD_IC(:,IKB:IKE) &
                  *XMD/(XAVOGADRO*PLMASS(:,IKB:IKE))
ZNOX_PROD_CG(:,IKB:IKE) = ZNOX_PROD_CG(:,IKB:IKE) &
                  *XMD/(XAVOGADRO*PLMASS(:,IKB:IKE))
!
!*      2.      Updraft computations
!               --------------------
!
ZUCH1(:,:,:) = 0.
!
!*      2.1     Initialization  at LCL
!               ----------------------
!
DO JI = 1, KLON
  JKLD = KDPL(JI)
  JKLP = KPBL(JI)
  ZWORK1(JI) = .5 * ( PCH1(JI,JKLD) + PCH1(JI,JKLP) )
END DO
!
!*      2.2     Final updraft loop
!               ------------------
!
DO JK = MINVAL( KDPL(:) ), JKMAX
  JKP = MIN(JK + 1,JKMAX)
  DO JI = 1, KLON
    IF ( KDPL(JI) <= JK     .AND. KLCL(JI) > JK ) THEN
      ZUCH1(JI,JK,1)  = ZWORK1(JI) 
      ZUCH1(JI,JK,2)  = ZWORK1(JI) + PTIMEC(JI)*(XFRAC_UIC*ZNOX_PROD_IC(JI,JK)&
                                                +XFRAC_UCG*ZNOX_PROD_CG(JI,JK))
    END IF
    IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
      ZUCH1(JI,JKP,1) = ( PUMF(JI,JK)  * ZUCH1(JI,JK,1) +                     &
                          PUER(JI,JKP) * PCH1(JI,JK)                          &
                        )  / ( PUMF(JI,JKP) + PUDR(JI,JKP) )
      IF ( KLFS(JI) -1 >= JK ) THEN
        ZUCH1(JI,JKP,2) = ( PUMF(JI,JK)  * ZUCH1(JI,JK,2) +                   &
                            PUER(JI,JKP) * PCH1(JI,JK)  +                     &
                        PLMASS(JI,JKP) * (XFRAC_UIC*ZNOX_PROD_IC(JI,JKP)      &
                                         +XFRAC_UCG*ZNOX_PROD_CG(JI,JKP))     &
                        )  / ( PUMF(JI,JKP) + PUDR(JI,JKP) )
      ELSE                                                                    
        ZUCH1(JI,JKP,2) = ( PUMF(JI,JK)  * ZUCH1(JI,JK,2) +                   &
                            PUER(JI,JKP) * PCH1(JI,JK)  +                     &
                        PLMASS(JI,JKP) * (ZNOX_PROD_IC(JI,JKP)                &
                                         +ZNOX_PROD_CG(JI,JKP))               &
                        )  / ( PUMF(JI,JKP) + PUDR(JI,JKP) )

      END IF
    END IF
  END DO
END DO
!
!*      3.      Downdraft computations
!               ----------------------
!
ZDCH1(:,:,:) = 0.
!
!*      3.1     Initialization at the LFS and at the DBL
!               ----------------------------------------
!
DO JI = 1, KLON
  JK = KLFS(JI)
  ZDCH1(JI,JK,1) = PMIXF(JI)*PCH1(JI,JK) + ( 1. - PMIXF(JI) )*ZUCH1(JI,JK,1)
  ZDCH1(JI,JK,2) = PMIXF(JI)*PCH1(JI,JK) + ( 1. - PMIXF(JI) )*ZUCH1(JI,JK,2)  &
                   + PTIMEC(JI) * ((1-XFRAC_UIC)*ZNOX_PROD_IC(JI,JK) +        &
                                   (1-XFRAC_UCG)*ZNOX_PROD_CG(JI,JK))
END DO
!
!*      3.2     Final downdraft loop
!               --------------------
!
DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
  JKM = JK - 1
  DO JI = 1, KLON
    IF ( JK <= KLFS(JI) .AND. JKM >= KDBL(JI) ) THEN
      ZDCH1(JI,JKM,1) = ( PDMF(JI,JK)  * ZDCH1(JI,JK,1) -                     &
                          PDER(JI,JKM) * PCH1(JI,JK)                          &
                        )  / ( PDMF(JI,JKM) - PDDR(JI,JKM) ) 
      ZDCH1(JI,JKM,2) = ( PDMF(JI,JK) *ZDCH1(JI,JK,2) -                       &
                          PDER(JI,JKM)*PCH1(JI,JK)  -                         &
                        PLMASS(JI,JKM)*((1.0-XFRAC_UIC)*ZNOX_PROD_IC(JI,JKM)  &
                                        +(1.0-XFRAC_UCG)*ZNOX_PROD_CG(JI,JKM))&
                        )  / ( PDMF(JI,JKM) - PDDR(JI,JKM) ) 
    END IF
  END DO
END DO
!							   
!*      4.      Final closure (environmental) computations
!               ------------------------------------------
!
ZCH1C(:,:,1) = ZCH1(:,:) ! initialize adjusted envir. values followed by an 
ZCH1C(:,:,2) = ZCH1(:,:) ! initialize adjusted envir. values followed by an 
!
DO JK = IKB, IKE
  ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / XG ! environmental subsidence
END DO
!
ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
                                         ! to be an integer multiple of PTIMEC
WHERE ( PTIMEC(:) < 1. ) ZTIMEC(:) = 0.
ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
!
DO jj=1,2
  ZCH1MFIN(:,:)   = 0.
  ZCH1MFOUT(:,:)  = 0.
!
  DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop
    DO JK = IKB + 1, JKMAX
      JKP = MAX( IKB + 1, JK - 1 )
      ZWORK3(:) = ZOMG(:,JK)
      ZWORK1(:) = SIGN( 1., ZWORK3(:) )
      ZWORK2(:) = 0.5 * ( 1. + ZWORK1(:) )
      ZWORK1(:) = 0.5 * ( 1. - ZWORK1(:) )
      ZCH1MFIN(:,JK)  = - ZWORK3(:) * ZCH1C(:,JKP,jj) * ZWORK1(:)
      ZCH1MFOUT(:,JK) =   ZWORK3(:) * ZCH1C(:,JK,jj)  * ZWORK2(:)
      ZCH1MFIN(:,JKP) = ZCH1MFIN(:,JKP) + ZCH1MFOUT(:,JK) * ZWORK2(:)
      ZCH1MFOUT(:,JKP)= ZCH1MFOUT(:,JKP) + ZCH1MFIN(:,JK) * ZWORK1(:)
    END DO
!
    DO JK = IKB + 1, JKMAX
      ZCH1C(:,JK,jj) = ZCH1C(:,JK,jj) + ZTIMC(:,JK) / PLMASS(:,JK) *  (    &
                   ZCH1MFIN(:,JK) + PUDR(:,JK) * ZUCH1(:,JK,jj) +       &
                   PDDR(:,JK) * ZDCH1(:,JK,jj) - ZCH1MFOUT(:,JK) -      &
                    ( PUER(:,JK) + PDER(:,JK) ) * PCH1(:,JK)    )
      ZCH1C(:,JK,jj) = MAX( 0., ZCH1C(:,JK,jj) )
    END DO
  END DO ! Exit the fractional time step loop
END DO
!
!
!----------------------------------------------------------------------------
!
!*           8.8    Apply conservation correction
!                   -----------------------------
!
ZCH1MFIN(:,:)= ZCH1C(:,:,2)- ZCH1C(:,:,1)
!
!
! Compute vertical integrals
!
ZWORK1(:) = 0.
ZWORK2(:) = 0.
DO JI = 1, KLON
  JKP = KCTL(JI)
  IF(JKP < IKB+1)CYCLE
  DO JK = IKB+1, JKP
    ZWORK1(JI) = ZWORK1(JI) + (ZCH1C(JI,JK,1)-PCH1(JI,JK)) *        &
                              .5 * (PLMASS(JI,JK-1) + PLMASS(JI,JK+1))
    ZWORK2(JI) = ZWORK2(JI) +                                       &
                          .5 * (PLMASS(JI,JK-1) + PLMASS(JI,JK+1))
  END DO
  ZWORK1(JI) = ZWORK1(JI) / ZWORK2(JI)
END DO
!
! Mass error (integral must be zero)
!
! Apply uniform correction but assure positive mass at each level
!
DO JI = 1, KLON
  JKP = KCTL(JI)
  IF(JKP < IKB+1)CYCLE
  DO JK = IKB+1, JKP
    ZCH1C(JI,JK,1) = ZCH1C(JI,JK,1) - ZWORK1(JI)
  END DO
END DO
!
! Compute vertical integrals
!
ZWORK1(:) = 0.
DO JI = 1, KLON
  JKP = KCTL(JI)
  IF(JKP < IKB+1)CYCLE
  DO JK = IKB+1, JKP
    ZWORK1(JI) = ZWORK1(JI) + (ZCH1C(JI,JK,1)-PCH1(JI,JK)) *        &
                            .5 * (PLMASS(JI,JK-1) + PLMASS(JI,JK+1))
  END DO
END DO
!
! Mass error (integral must be zero)
!
!
! Add final
!
PCH1C(:,:) = ZCH1C(:,:,1) + ZCH1MFIN(:,:)
PCH1C(:,:) = MAX(0.,PCH1C(:,:))
!
!----------------------------------------------------------------------------
!
!
END SUBROUTINE CH_CONVECT_LINOX
