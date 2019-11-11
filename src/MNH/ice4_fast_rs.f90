!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_FAST_RS
INTERFACE
SUBROUTINE ICE4_FAST_RS(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAR, PLBDAS, &
                       &PT,  PRVT, PRCT, PRRT, PRST, &
                       &PRIAGGS, &
                       &PRCRIMSS, PRCRIMSG, PRSRIMCG, &
                       &PRRACCSS, PRRACCSG, PRSACCRG, PRSMLTG, &
                       &PRCMLTSR, &
                       &PRS_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RS, PA_RG)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIAGGS  ! r_i aggregation on r_s
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCRIMSS ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCRIMSG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSRIMCG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRACCSS ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRACCSG ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSACCRG ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSMLTG  ! Conversion-Melting of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCMLTSR ! Cloud droplet collection onto aggregates by positive temperature
REAL, DIMENSION(KSIZE, 6),    INTENT(INOUT) :: PRS_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
END SUBROUTINE ICE4_FAST_RS
END INTERFACE
END MODULE MODI_ICE4_FAST_RS
SUBROUTINE ICE4_FAST_RS(KSIZE, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAR, PLBDAS, &
                       &PT,  PRVT, PRCT, PRRT, PRST, &
                       &PRIAGGS, &
                       &PRCRIMSS, PRCRIMSG, PRSRIMCG, &
                       &PRRACCSS, PRRACCSG, PRSACCRG, PRSMLTG, &
                       &PRCMLTSR, &
                       &PRS_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RS, PA_RG)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast rs processes
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
USE MODD_PARAM_ICE, ONLY : LEVLIMIT, CSNOWRIMING
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIAGGS  ! r_i aggregation on r_s
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCRIMSS ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCRIMSG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSRIMCG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRACCSS ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRACCSG ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSACCRG ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSMLTG  ! Conversion-Melting of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCMLTSR ! Cloud droplet collection onto aggregates by positive temperature
REAL, DIMENSION(KSIZE, 6),    INTENT(INOUT) :: PRS_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCRIMS=1, IRCRIMSS=2, IRSRIMCG=3, IRRACCS=4, IRRACCSS=5, IRSACCRG=6
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GRIM, GACC, GMASK
INTEGER :: IGRIM, IGACC
REAL, DIMENSION(SIZE(PRHODREF)) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(SIZE(PRHODREF)) :: IVEC1, IVEC2
REAL, DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZW6, ZFREEZ_RATE
INTEGER :: JJ
!-------------------------------------------------------------------------------
!
!
!*       5.0    maximum freezing rate
!
ZFREEZ_RATE(:)=0.
GMASK(:)=PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:)
WHERE(GMASK(:))
  ZFREEZ_RATE(:)=PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
END WHERE
IF(LEVLIMIT) THEN
  WHERE(GMASK(:))
    ZFREEZ_RATE(:)=MIN(ZFREEZ_RATE(:), EXP(XALPI-XBETAI/PT(:)-XGAMI*ALOG(PT(:)))) ! min(ev, es_i(T))
  END WHERE
ENDIF
WHERE(GMASK(:))
  ZFREEZ_RATE(:)=PKA(:)*(XTT-PT(:)) +                              &
           (PDV(:)*(XLVTT+(XCPV-XCL)*(PT(:)-XTT)) &
                         *(XESTT-ZFREEZ_RATE(:))/(XRV*PT(:))           )
  ZFREEZ_RATE(:)=MAX(0.,                                               &
             (ZFREEZ_RATE(:) * ( X0DEPS*       PLBDAS(:)**XEX0DEPS +     &
                         X1DEPS*PCJ(:)*PLBDAS(:)**XEX1DEPS ) +   &
             PRIAGGS(:) *                            &
             (PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PT(:)))   ) ) / &
                        ( PRHODREF(:)*(XLMTT-XCL*(XTT-PT(:))) )   )
  !We must agregate, at least, the cold species
  !And we are only interested by the freezing rate of liquid species
  ZFREEZ_RATE(:)=MAX(ZFREEZ_RATE(:)-PRIAGGS(:), 0.)
END WHERE
!
!*       5.1    cloud droplet riming of the aggregates
!
GRIM(:) = PRCT(:)>XRTMIN(2) .AND. PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:)
!
! Collection of cloud droplets by snow: this rate is used for riming (T<0) and for conversion/melting (T>0)
IF(LDSOFT) THEN
  WHERE(.NOT. GRIM(:))
    PRS_TEND(:, IRCRIMS)=0.
    PRS_TEND(:, IRCRIMSS)=0.
    PRS_TEND(:, IRSRIMCG)=0.
  END WHERE
ELSE
  PRS_TEND(:, IRCRIMS)=0.
  PRS_TEND(:, IRCRIMSS)=0.
  PRS_TEND(:, IRSRIMCG)=0.
  IGRIM = COUNT(GRIM(:))
  !
  IF(IGRIM>0) THEN
    !
    !        5.1.1  select the PLBDAS
    !
    ZVEC1(1:IGRIM) = PACK( PLBDAS(:),MASK=GRIM(:) )
    !
    !        5.1.2  find the next lower indice for the PLBDAS in the geometrical
    !               set of Lbda_s used to tabulate some moments of the incomplete
    !               gamma function
    !
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( FLOAT(NGAMINC)-0.00001,           &
                          XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - FLOAT( IVEC2(1:IGRIM) )
    !
    !        5.1.3  perform the linear interpolation of the normalized
    !               "2+XDS"-moment of the incomplete gamma function
    !
    ZVEC1(1:IGRIM) =   XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(1:IGRIM),MASK=GRIM,FIELD=0.0 )
    !
    !        5.1.4  riming of the small sized aggregates
    !
    WHERE (GRIM(:))
      PRS_TEND(:, IRCRIMSS) = XCRIMSS * ZZW(:) * PRCT(:)                & ! RCRIMSS
                                      *   PLBDAS(:)**XEXCRIMSS &
                                      * PRHODREF(:)**(-XCEXVT)
    END WHERE
    !
    !        5.1.5  perform the linear interpolation of the normalized
    !               "XBS"-moment of the incomplete gamma function (XGAMINC_RIM2) and
    !               "XBG"-moment of the incomplete gamma function (XGAMINC_RIM4)
    !
    ZVEC1(1:IGRIM) =  XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = UNPACK( VECTOR=ZVEC1(1:IGRIM),MASK=GRIM,FIELD=0.0 )

    ZVEC1(1:IGRIM) =  XGAMINC_RIM4( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM4( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW2(:) = UNPACK( VECTOR=ZVEC1(1:IGRIM),MASK=GRIM,FIELD=0.0)
    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
    WHERE(GRIM(:))
      PRS_TEND(:, IRCRIMS)=XCRIMSG * PRCT(:)               & ! RCRIMS
                                   * PLBDAS(:)**XEXCRIMSG  &
                                   * PRHODREF(:)**(-XCEXVT)
      ZZW6(:) = PRS_TEND(:, IRCRIMS) - PRS_TEND(:, IRCRIMSS) ! RCRIMSG
    END WHERE

    IF(CSNOWRIMING=='M90 ')THEN
      !Murakami 1990
      WHERE(GRIM(:))
        PRS_TEND(:, IRSRIMCG)=XSRIMCG * PLBDAS(:)**XEXSRIMCG*(1.0-ZZW(:))
        PRS_TEND(:, IRSRIMCG)=ZZW6(:)*PRS_TEND(:, IRSRIMCG)/ &
                       MAX(1.E-20, &
                           XSRIMCG3*XSRIMCG2*PLBDAS(:)**XEXSRIMCG2*(1.-ZZW2(:)) - &
                           XSRIMCG3*PRS_TEND(:, IRSRIMCG))
      END WHERE
    ELSE
      PRS_TEND(:, IRSRIMCG)=0.
    END IF
  ENDIF
ENDIF
!
GRIM(:) = GRIM(:) .AND. PT(:)<XTT ! More restrictive GRIM mask to be used for riming by negative temperature only
PRCRIMSS(:)=0.
PRCRIMSG(:)=0.
PRSRIMCG(:)=0.
WHERE(GRIM(:))
  PRCRIMSS(:) = MIN(ZFREEZ_RATE(:), PRS_TEND(:, IRCRIMSS))
  ZFREEZ_RATE(:) = MAX(0., ZFREEZ_RATE(:)-PRCRIMSS(:))
  ZZW(:) = MIN(1., ZFREEZ_RATE(:) / MAX(1.E-20, PRS_TEND(:, IRCRIMS) - PRCRIMSS(:))) ! proportion we are able to freeze
  PRCRIMSG(:) = ZZW(:) * MAX(0., PRS_TEND(:, IRCRIMS) - PRCRIMSS(:)) ! RCRIMSG
  ZFREEZ_RATE(:) = MAX(0., ZFREEZ_RATE(:)-PRCRIMSG(:))
  PRSRIMCG(:) = ZZW(:) * PRS_TEND(:, IRSRIMCG)
END WHERE
WHERE(PRCRIMSG(:)<=0.)
  PRCRIMSG(:)=0.
  PRSRIMCG(:)=0.
END WHERE
PA_RC(:) = PA_RC(:) - PRCRIMSS(:)
PA_RS(:) = PA_RS(:) + PRCRIMSS(:)
PA_TH(:) = PA_TH(:) + PRCRIMSS(:)*(PLSFACT(:)-PLVFACT(:))
PA_RC(:) = PA_RC(:) - PRCRIMSG(:)
PA_RS(:) = PA_RS(:) - PRSRIMCG(:)
PA_RG(:) = PA_RG(:) + PRCRIMSG(:)+PRSRIMCG(:)
PA_TH(:) = PA_TH(:) + PRCRIMSG(:)*(PLSFACT(:)-PLVFACT(:))
!
!*       5.2    rain accretion onto the aggregates
!
GACC(:) = PRRT(:)>XRTMIN(3) .AND. PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GACC(:))
    PRS_TEND(:, IRRACCS)=0.
    PRS_TEND(:, IRRACCSS)=0.
    PRS_TEND(:, IRSACCRG)=0.
  END WHERE
ELSE
  PRS_TEND(:, IRRACCS)=0.
  PRS_TEND(:, IRRACCSS)=0.
  PRS_TEND(:, IRSACCRG)=0.
  IGACC = COUNT(GACC(:))
  IF(IGACC>0)THEN
    !
    !
    !        5.2.1  select the (PLBDAS,PLBDAR) couplet
    !
    ZVEC1(1:IGACC) = PACK( PLBDAS(:),MASK=GACC(:) )
    ZVEC2(1:IGACC) = PACK( PLBDAR(:),MASK=GACC(:) )
    !
    !        5.2.2  find the next lower indice for the PLBDAS and for the PLBDAR
    !               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
    !               tabulate the RACCSS-kernel
    !
    ZVEC1(1:IGACC) = MAX( 1.00001, MIN( FLOAT(NACCLBDAS)-0.00001,           &
                          XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + XACCINTP2S ) )
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - FLOAT( IVEC1(1:IGACC) )
    !
    ZVEC2(1:IGACC) = MAX( 1.00001, MIN( FLOAT(NACCLBDAR)-0.00001,           &
                          XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + XACCINTP2R ) )
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - FLOAT( IVEC2(1:IGACC) )
    !
    !        5.2.3  perform the bilinear interpolation of the normalized
    !               RACCSS-kernel
    !
    DO JJ = 1, IGACC
      ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                 - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGACC),MASK=GACC,FIELD=0.0 )
    !
    !        5.2.4  raindrop accretion on the small sized aggregates
    !
    WHERE(GACC(:))
      ZZW6(:) =                                                        & !! coef of RRACCS
            XFRACCSS*( PLBDAS(:)**XCXS )*( PRHODREF(:)**(-XCEXVT-1.) ) &
       *( XLBRACCS1/((PLBDAS(:)**2)               ) +                  &
          XLBRACCS2/( PLBDAS(:)    * PLBDAR(:)    ) +                  &
          XLBRACCS3/(               (PLBDAR(:)**2)) )/PLBDAR(:)**4
      PRS_TEND(:, IRRACCSS) =ZZW(:)*ZZW6(:)
    END WHERE
    !
    !        5.2.4b perform the bilinear interpolation of the normalized
    !               RACCS-kernel
    !
    DO JJ = 1, IGACC
      ZVEC3(JJ) =  (   XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                 - (   XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGACC),MASK=GACC(:),FIELD=0.0 )
    WHERE(GACC(:))
      PRS_TEND(:, IRRACCS) = ZZW(:)*ZZW6(:)
    END WHERE
    !        5.2.5  perform the bilinear interpolation of the normalized
    !               SACCRG-kernel
    !
    DO JJ = 1, IGACC
        ZVEC3(JJ) =  (  XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                            * ZVEC2(JJ) &
                   - (  XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                            * (ZVEC2(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGACC),MASK=GACC,FIELD=0.0 )
    !
    !        5.2.6  raindrop accretion-conversion of the large sized aggregates
    !               into graupeln
    !
    WHERE(GACC(:))
      PRS_TEND(:, IRSACCRG) = XFSACCRG*ZZW(:)*                    & ! RSACCRG
          ( PLBDAS(:)**(XCXS-XBS) )*( PRHODREF(:)**(-XCEXVT-1.) ) &
         *( XLBSACCR1/((PLBDAR(:)**2)               ) +           &
            XLBSACCR2/( PLBDAR(:)    * PLBDAS(:)    ) +           &
            XLBSACCR3/(               (PLBDAS(:)**2)) )/PLBDAR(:)
    END WHERE
  ENDIF
ENDIF
!
GACC(:) = GACC(:) .AND. PT(:)<XTT ! More restrictive GACC mask to be used for accretion by negative temperature only
PRRACCSS(:)=0.
PRRACCSG(:)=0.
PRSACCRG(:)=0.
WHERE(GACC(:))
  PRRACCSS(:) = MIN(ZFREEZ_RATE(:), PRS_TEND(:, IRRACCSS))
  ZFREEZ_RATE(:) = MAX(0., ZFREEZ_RATE(:)-PRRACCSS(:))
  ZZW(:) = MIN(1., ZFREEZ_RATE(:) / MAX(1.E-20, PRS_TEND(:, IRRACCS)-PRRACCSS(:))) ! proportion we are able to freeze
  PRRACCSG(:)=ZZW(:) * MAX(0., PRS_TEND(:, IRRACCS)-PRRACCSS(:))
  ZFREEZ_RATE(:) = MAX(0., ZFREEZ_RATE(:)-PRRACCSG(:))
  PRSACCRG(:)=ZZW(:) * PRS_TEND(:, IRSACCRG)
END WHERE
WHERE(PRRACCSG(:)<=0.)
  PRRACCSG(:)=0.
  PRSACCRG(:)=0.
END WHERE
PA_RR(:) = PA_RR(:) - PRRACCSS(:)
PA_RS(:) = PA_RS(:) + PRRACCSS(:)
PA_TH(:) = PA_TH(:) + PRRACCSS(:)*(PLSFACT(:)-PLVFACT(:))
PA_RR(:) = PA_RR(:) - PRRACCSG(:)
PA_RS(:) = PA_RS(:) - PRSACCRG(:)
PA_RG(:) = PA_RG(:) + PRRACCSG(:)+PRSACCRG(:)
PA_TH(:) = PA_TH(:) + PRRACCSG(:)*(PLSFACT(:)-PLVFACT(:))
!
!
!*       5.3    Conversion-Melting of the aggregates
!
GMASK(:)=PRST(:)>XRTMIN(5) .AND. PT(:)>XTT .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRSMLTG(:) = 0.
    PRCMLTSR(:) = 0.
  END WHERE
ELSE
  PRSMLTG(:) = 0.
  PRCMLTSR(:) = 0.
  WHERE(GMASK(:))
    PRSMLTG(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
  END WHERE
  IF(LEVLIMIT) THEN
    WHERE(GMASK(:))
      PRSMLTG(:)=MIN(PRSMLTG(:), EXP(XALPW-XBETAW/PT(:)-XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  WHERE(GMASK(:))
    PRSMLTG(:) =  PKA(:)*(XTT-PT(:)) +                                 &
               ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                           *(XESTT-PRSMLTG(:))/(XRV*PT(:))             )
    !
    ! compute RSMLT
    !
    PRSMLTG(:)  = XFSCVMG*MAX( 0.0,( -PRSMLTG(:) *             &
                         ( X0DEPS*       PLBDAS(:)**XEX0DEPS +     &
                           X1DEPS*PCJ(:)*PLBDAS(:)**XEX1DEPS ) -   &
                                   ( PRS_TEND(:, IRCRIMS) + PRS_TEND(:, IRRACCS) ) *       &
                            ( PRHODREF(:)*XCL*(XTT-PT(:))) ) /    &
                                           ( PRHODREF(:)*XLMTT ) )
    !
    ! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
    ! because the graupeln produced by this process are still icy!!!
    !
    ! When T < XTT, rc is collected by snow (riming) to produce snow and graupel
    ! When T > XTT, if riming was still enabled, rc would produce snow and graupel with snow becomming graupel (conversion/melting) and graupel becomming rain (melting)
    ! To insure consistency when crossint T=XTT, rc collected with T>XTT must be transformed in rain.
    ! rc cannot produce iced species with a positive temperature but is still collected with a good efficiency by snow
    PRCMLTSR(:) = PRS_TEND(:, IRCRIMS) ! both species are liquid, no heat is exchanged
  END WHERE
ENDIF
PA_RS(:) = PA_RS(:) - PRSMLTG(:)
PA_RG(:) = PA_RG(:) + PRSMLTG(:)
PA_RC(:) = PA_RC(:) - PRCMLTSR(:)
PA_RR(:) = PA_RR(:) + PRCMLTSR(:)

!
END SUBROUTINE ICE4_FAST_RS
