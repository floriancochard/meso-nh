!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_FAST_RG
INTERFACE
SUBROUTINE ICE4_FAST_RG(KSIZE, LDSOFT, LDCOMPUTE, KRR, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, PCIT, &
                       &PLBDAR, PLBDAS, PLBDAG, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                       &PRGSI, PRGSI_MR, &
                       &LDWETG, &
                       &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                       &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                       &PRG_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH, PB_RG, PB_RH)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
INTEGER,                      INTENT(IN)    :: KRR      ! Number of moist variable
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI    ! Graupel tendency by other processes
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI_MR ! Graupel mr change by other processes
LOGICAL, DIMENSION(KSIZE),    INTENT(OUT)   :: LDWETG   ! True where graupel grows in wet mode
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRRG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRCFRIG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRR  ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRWETGH  ! Conversion of graupel into hail
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRWETGH_MR ! Conversion of graupel into hail, mr change
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGMLTR  ! Melting of the graupel
REAL, DIMENSION(KSIZE, 6),    INTENT(INOUT) :: PRG_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RH
END SUBROUTINE ICE4_FAST_RG
END INTERFACE
END MODULE MODI_ICE4_FAST_RG
SUBROUTINE ICE4_FAST_RG(KSIZE, LDSOFT, LDCOMPUTE, KRR, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, PCIT, &
                       &PLBDAR, PLBDAS, PLBDAG, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                       &PRGSI, PRGSI_MR, &
                       &LDWETG, &
                       &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                       &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                       &PRG_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH, PB_RG, PB_RH)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast rg processes
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
USE MODD_PARAM_ICE, ONLY : LEVLIMIT, LNULLWETG, LWETGPOST, LCRFLIMIT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
INTEGER,                      INTENT(IN)    :: KRR      ! Number of moist variable
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI    ! Graupel tendency by other processes
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI_MR ! Graupel mr change by other processes
LOGICAL, DIMENSION(KSIZE),    INTENT(OUT)   :: LDWETG   ! True where graupel grows in wet mode
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRRG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRCFRIG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRR  ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRWETGH  ! Conversion of graupel into hail
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRWETGH_MR ! Conversion of graupel into hail, mr change
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGMLTR  ! Melting of the graupel
REAL, DIMENSION(KSIZE, 6),    INTENT(INOUT) :: PRG_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RH
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCDRYG=1, IRIDRYG=2, IRIWETG=3, IRSDRYG=4, IRSWETG=5, IRRDRYG=6
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GDRY, LLDRYG, GMASK
INTEGER :: IGDRY
REAL, DIMENSION(SIZE(PRHODREF)) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(SIZE(PRHODREF)) :: IVEC1, IVEC2
REAL, DIMENSION(SIZE(PRHODREF)) :: ZZW, &
                                   ZRDRYG_INIT, & !Initial dry growth rate of the graupeln
                                   ZRWETG_INIT !Initial wet growth rate of the graupeln
INTEGER :: JJ
!
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing
!
GMASK(:)=PRIT(:)>XRTMIN(4) .AND. PRRT(:)>XRTMIN(3) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRICFRRG(:)=0.
    PRRCFRIG(:)=0.
    PRICFRR(:)=0.
  ENDWHERE
ELSE
  PRICFRRG(:)=0.
  PRRCFRIG(:)=0.
  PRICFRR(:)=0.
  WHERE(GMASK(:))
    PRICFRRG(:) = XICFRR*PRIT(:)                & ! RICFRRG
                                 *PLBDAR(:)**XEXICFRR    &
                                 *PRHODREF(:)**(-XCEXVT)
    PRRCFRIG(:) = XRCFRI*PCIT(:)                & ! RRCFRIG
                                 * PLBDAR(:)**XEXRCFRI    &
                                 * PRHODREF(:)**(-XCEXVT-1.)
  END WHERE
  ZZW(:)=1.
  IF(LCRFLIMIT) THEN
    WHERE(GMASK(:))
      !Comparison between heat to be released (to freeze rain) and heat sink (rain and ice temperature change)
      !ZZW is the proportion of process that can take place
      ZZW(:) = MAX(0., MIN(1., (PRICFRRG(:)*XCI+PRRCFRIG(:)*XCL)*(XTT-PT(:)) / &
                               MAX(1.E-20, XLVTT*PRRCFRIG(:))))
    ENDWHERE
  ENDIF
  PRRCFRIG(:) = ZZW(:) * PRRCFRIG(:) !Part of rain that can be freezed
  PRICFRR(:) = (1-ZZW(:)) * PRICFRRG(:) !Part of collected pristine ice converted to rain
  PRICFRRG(:) = ZZW(:) * PRICFRRG(:) !Part of collected pristine ice that lead to graupel
ENDIF
PA_RI(:) = PA_RI(:) - PRICFRRG(:) - PRICFRR(:)
PA_RR(:) = PA_RR(:) - PRRCFRIG(:) + PRICFRR(:)
PA_RG(:) = PA_RG(:) + PRICFRRG(:) + PRRCFRIG(:)
PA_TH(:) = PA_TH(:) + (PRRCFRIG(:) - PRICFRR(:))*(PLSFACT(:)-PLVFACT(:))
!
!
!*       6.3    compute the graupel growth
!
! Wet and dry collection of rc and ri on graupel
GMASK(:)=PRGT(:)>XRTMIN(6) .AND. PRCT(:)>XRTMIN(2) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRG_TEND(:, IRCDRYG)=0.
  END WHERE
ELSE
  PRG_TEND(:, IRCDRYG)=0.
  WHERE(GMASK(:))
    ZZW(:)=PLBDAG(:)**(XCXG-XDG-2.) * PRHODREF(:)**(-XCEXVT)
    PRG_TEND(:, IRCDRYG)=XFCDRYG * PRCT(:) * ZZW(:)
  END WHERE
ENDIF
GMASK(:)=PRGT(:)>XRTMIN(6) .AND. PRIT(:)>XRTMIN(4) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRG_TEND(:, IRIDRYG)=0.
    PRG_TEND(:, IRIWETG)=0.
  END WHERE
ELSE
  PRG_TEND(:, IRIDRYG)=0.
  PRG_TEND(:, IRIWETG)=0.
  WHERE(GMASK(:))
    ZZW(:)=PLBDAG(:)**(XCXG-XDG-2.) * PRHODREF(:)**(-XCEXVT)
    PRG_TEND(:, IRIDRYG)=XFIDRYG*EXP(XCOLEXIG*(PT(:)-XTT))*PRIT(:)*ZZW(:)
    PRG_TEND(:, IRIWETG)=PRG_TEND(:, IRIDRYG) / (XCOLIG*EXP(XCOLEXIG*(PT(:)-XTT)))
  END WHERE
ENDIF

! Wet and dry collection of rs on graupel (6.2.1)
GDRY(:)=PRST(:)>XRTMIN(5) .AND. PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GDRY(:))
    PRG_TEND(:, IRSDRYG)=0.
    PRG_TEND(:, IRSWETG)=0.
  END WHERE
ELSE
  PRG_TEND(:, IRSDRYG)=0.
  PRG_TEND(:, IRSWETG)=0.
  IGDRY=COUNT(GDRY(:))
  IF(IGDRY>0)THEN
    !
    !*       6.2.3  select the (PLBDAG,PLBDAS) couplet
    !
    ZVEC1(1:IGDRY)=PACK(PLBDAG(:), MASK=GDRY(:))
    ZVEC2(1:IGDRY)=PACK(PLBDAS(:), MASK=GDRY(:))
    !
    !*       6.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
    !               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
    !               tabulate the SDRYG-kernel
    !
    ZVEC1(1:IGDRY)=MAX(1.00001, MIN(FLOAT(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G*LOG(ZVEC1(1:IGDRY))+XDRYINTP2G))
    IVEC1(1:IGDRY)=INT(ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY)=ZVEC1(1:IGDRY)-FLOAT(IVEC1(1:IGDRY))
    !
    ZVEC2(1:IGDRY)=MAX(1.00001, MIN( FLOAT(NDRYLBDAS)-0.00001,           &
                          XDRYINTP1S*LOG(ZVEC2(1:IGDRY))+XDRYINTP2S))
    IVEC2(1:IGDRY)=INT(ZVEC2(1:IGDRY))
    ZVEC2(1:IGDRY)=ZVEC2(1:IGDRY)-FLOAT(IVEC2(1:IGDRY))
    !
    !*       6.2.5  perform the bilinear interpolation of the normalized
    !               SDRYG-kernel
    !
    DO JJ=1, IGDRY
      ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * ZVEC1(JJ) &
                 - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         *(ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:)=UNPACK(VECTOR=ZVEC3(1:IGDRY), MASK=GDRY(:), FIELD=0.0)
    !
    WHERE(GDRY(:))
      PRG_TEND(:, IRSWETG)=XFSDRYG*ZZW(:)                         & ! RSDRYG
                                    / XCOLSG &
                  *(PLBDAS(:)**(XCXS-XBS))*( PLBDAG(:)**XCXG )    &
                  *(PRHODREF(:)**(-XCEXVT-1.))                    &
                       *( XLBSDRYG1/( PLBDAG(:)**2              ) + &
                          XLBSDRYG2/( PLBDAG(:)   * PLBDAS(:)   ) + &
                          XLBSDRYG3/(               PLBDAS(:)**2))
      PRG_TEND(:, IRSDRYG)=PRG_TEND(:, IRSWETG)*XCOLSG*EXP(XCOLEXSG*(PT(:)-XTT))
    END WHERE
  ENDIF
ENDIF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
GDRY(:)=PRRT(:)>XRTMIN(3) .AND. PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GDRY(:))
    PRG_TEND(:, IRRDRYG)=0.
  END WHERE
ELSE
  PRG_TEND(:, IRRDRYG)=0.
  IGDRY=COUNT(GDRY(:))
  !
  IF(IGDRY>0) THEN
    !
    !*       6.2.8  select the (PLBDAG,PLBDAR) couplet
    !
    ZVEC1(1:IGDRY)=PACK(PLBDAG(:), MASK=GDRY(:))
    ZVEC2(1:IGDRY)=PACK(PLBDAR(:), MASK=GDRY(:))
    !
    !*       6.2.9  find the next lower indice for the PLBDAG and for the PLBDAR
    !               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
    !               tabulate the RDRYG-kernel
    !
    ZVEC1(1:IGDRY)=MAX(1.00001, MIN( FLOAT(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G*LOG(ZVEC1(1:IGDRY))+XDRYINTP2G))
    IVEC1(1:IGDRY)=INT(ZVEC1(1:IGDRY))
    ZVEC1(1:IGDRY)=ZVEC1(1:IGDRY)-FLOAT(IVEC1(1:IGDRY))
    !
    ZVEC2(1:IGDRY)=MAX(1.00001, MIN( FLOAT(NDRYLBDAR)-0.00001,           &
                          XDRYINTP1R*LOG(ZVEC2(1:IGDRY))+XDRYINTP2R))
    IVEC2(1:IGDRY)=INT(ZVEC2(1:IGDRY))
    ZVEC2(1:IGDRY)=ZVEC2(1:IGDRY)-FLOAT(IVEC2(1:IGDRY))
    !
    !*       6.2.10 perform the bilinear interpolation of the normalized
    !               RDRYG-kernel
    !
    DO JJ=1, IGDRY
      ZVEC3(JJ)= (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         *(ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:)=UNPACK(VECTOR=ZVEC3(1:IGDRY), MASK=GDRY, FIELD=0.)
    !
    WHERE(GDRY(:))
      PRG_TEND(:, IRRDRYG) = XFRDRYG*ZZW(:)                    & ! RRDRYG
                        *( PLBDAR(:)**(-4) )*( PLBDAG(:)**XCXG ) &
                               *( PRHODREF(:)**(-XCEXVT-1.) )   &
                    *( XLBRDRYG1/( PLBDAG(:)**2              ) + &
                       XLBRDRYG2/( PLBDAG(:)   * PLBDAR(:)   ) + &
                       XLBRDRYG3/(               PLBDAR(:)**2) )
    END WHERE
  ENDIF
ENDIF

ZRDRYG_INIT(:)=PRG_TEND(:, IRCDRYG)+PRG_TEND(:, IRIDRYG)+PRG_TEND(:, IRSDRYG)+PRG_TEND(:, IRRDRYG)

!Freezing rate
ZRWETG_INIT(:)=0.
GMASK(:)=PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:)
WHERE(GMASK(:))
  ZRWETG_INIT(:)=PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
END WHERE
IF(LEVLIMIT) THEN
  WHERE(GMASK(:))
    ZRWETG_INIT(:)=MIN(ZRWETG_INIT(:), EXP(XALPI-XBETAI/PT(:)-XGAMI*ALOG(PT(:)))) ! min(ev, es_i(T))
  END WHERE
ENDIF
WHERE(GMASK(:))
  ZRWETG_INIT(:)=PKA(:)*(XTT-PT(:)) +                              &
           (PDV(:)*(XLVTT+(XCPV-XCL)*(PT(:)-XTT)) &
                         *(XESTT-ZRWETG_INIT(:))/(XRV*PT(:))           )
  ZRWETG_INIT(:)=MAX(0.,                                               &
             (ZRWETG_INIT(:) * ( X0DEPG*       PLBDAG(:)**XEX0DEPG +     &
                         X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) +   &
             (PRG_TEND(:, IRIWETG)+PRG_TEND(:, IRSWETG) ) *                            &
             (PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PT(:)))   ) ) / &
                        ( PRHODREF(:)*(XLMTT-XCL*(XTT-PT(:))) )   )
  !We must agregate, at least, the cold species
  ZRWETG_INIT(:)=MAX(ZRWETG_INIT(:), PRG_TEND(:, IRIWETG)+PRG_TEND(:, IRSWETG))
END WHERE

!Growth mode
LDWETG(:)=GMASK(:) .AND. &
         &MAX(0., ZRDRYG_INIT(:)-PRG_TEND(:, IRIDRYG)-PRG_TEND(:, IRSDRYG))>= &
         &MAX(0., ZRWETG_INIT(:)-PRG_TEND(:, IRIWETG)-PRG_TEND(:, IRSWETG))
IF(LNULLWETG) THEN
  LDWETG(:)=LDWETG(:) .AND. ZRDRYG_INIT(:)>0.
ELSE
  LDWETG(:)=LDWETG(:) .AND. ZRWETG_INIT(:)>0.
ENDIF
IF(.NOT. LWETGPOST) LDWETG(:)=LDWETG(:) .AND. PT(:)<XTT

LLDRYG(:)=GMASK(:) .AND. PT(:)<XTT .AND. ZRDRYG_INIT(:)>0. .AND. &
    &MAX(0., ZRDRYG_INIT(:)-PRG_TEND(:, IRIDRYG)-PRG_TEND(:, IRSDRYG))<&
    &MAX(0., ZRWETG_INIT(:)-PRG_TEND(:, IRIWETG)-PRG_TEND(:, IRSWETG))

! Part of ZRWETG to be converted into hail
! Graupel can be produced by other processes instantaneously (inducing a mixing ratio change, PRGSI_MR) or
! as a tendency (PRWETGH)
PRWETGH(:)=0.
PRWETGH_MR(:)=0.
IF(KRR==7) THEN
  WHERE(LDWETG(:))
    !assume a linear percent of conversion of produced graupel into hail
    PRWETGH(:)=(MAX(0., PRGSI(:)+PRICFRRG(:)+PRRCFRIG(:))+ZRWETG_INIT(:))*ZRDRYG_INIT(:)/(ZRWETG_INIT(:)+ZRDRYG_INIT(:))
    PRWETGH_MR(:)=MAX(0., PRGSI_MR(:))*ZRDRYG_INIT(:)/(ZRWETG_INIT(:)+ZRDRYG_INIT(:))
  END WHERE
ENDIF

PRCWETG(:)=0.
PRIWETG(:)=0.
PRSWETG(:)=0.
PRRWETG(:)=0.
WHERE(LDWETG(:))
  !Aggregated minus collected
  PRRWETG(:)=-(PRG_TEND(:, IRIWETG)+PRG_TEND(:, IRSWETG)+PRG_TEND(:, IRCDRYG)-ZRWETG_INIT(:))
  PRCWETG(:)=PRG_TEND(:, IRCDRYG)
  PRIWETG(:)=PRG_TEND(:, IRIWETG)
  PRSWETG(:)=PRG_TEND(:, IRSWETG)
END WHERE
PRCDRYG(:)=0.
PRIDRYG(:)=0.
PRRDRYG(:)=0.
PRSDRYG(:)=0.
WHERE(LLDRYG(:))
  PRCDRYG(:)=PRG_TEND(:, IRCDRYG)
  PRRDRYG(:)=PRG_TEND(:, IRRDRYG)
  PRIDRYG(:)=PRG_TEND(:, IRIDRYG)
  PRSDRYG(:)=PRG_TEND(:, IRSDRYG)
END WHERE
PA_RC(:) = PA_RC(:) - PRCWETG(:)
PA_RI(:) = PA_RI(:) - PRIWETG(:)
PA_RS(:) = PA_RS(:) - PRSWETG(:)
PA_RG(:) = PA_RG(:) + PRCWETG(:) + PRIWETG(:) + PRSWETG(:) + PRRWETG(:)
PA_RR(:) = PA_RR(:) - PRRWETG(:)
PA_TH(:) = PA_TH(:) + (PRCWETG(:) + PRRWETG(:))*(PLSFACT(:)-PLVFACT(:))
PA_RG(:) = PA_RG(:) - PRWETGH(:)
PA_RH(:) = PA_RH(:) + PRWETGH(:)
PB_RG(:) = PB_RG(:) - PRWETGH_MR(:)
PB_RH(:) = PB_RH(:) + PRWETGH_MR(:)
PA_RC(:) = PA_RC(:) - PRCDRYG(:)
PA_RI(:) = PA_RI(:) - PRIDRYG(:)
PA_RS(:) = PA_RS(:) - PRSDRYG(:)
PA_RR(:) = PA_RR(:) - PRRDRYG(:)
PA_RG(:) = PA_RG(:) + PRCDRYG(:) + PRIDRYG(:) + PRSDRYG(:) + PRRDRYG(:)
PA_TH(:) = PA_TH(:) + (PRCDRYG(:)+PRRDRYG(:))*(PLSFACT(:)-PLVFACT(:))

!
!*       6.5    Melting of the graupeln
!
GMASK(:)=PRGT(:)>XRTMIN(6) .AND. PT(:)>XTT .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRGMLTR(:) = 0.
  END WHERE
ELSE
  PRGMLTR(:) = 0.
  WHERE(GMASK(:))
    PRGMLTR(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
  END WHERE
  IF(LEVLIMIT) THEN
    WHERE(GMASK(:))
      PRGMLTR(:)=MIN(PRGMLTR(:), EXP(XALPW-XBETAW/PT(:)-XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  WHERE(GMASK(:))
    PRGMLTR(:) =  PKA(:)*(XTT-PT(:)) +                                 &
               ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                           *(XESTT-PRGMLTR(:))/(XRV*PT(:))             )
  END WHERE
  WHERE(GMASK(:))
    !
    ! compute RGMLTR
    !
    PRGMLTR(:)  = MAX( 0.0,( -PRGMLTR(:) *                     &
                           ( X0DEPG*       PLBDAG(:)**XEX0DEPG +     &
                             X1DEPG*PCJ(:)*PLBDAG(:)**XEX1DEPG ) -   &
                         ( PRG_TEND(:, IRCDRYG)+PRG_TEND(:, IRRDRYG) ) *       &
                               ( PRHODREF(:)*XCL*(XTT-PT(:))) ) /    &
                                             ( PRHODREF(:)*XLMTT ) )
  END WHERE
ENDIF
PA_RR(:) = PA_RR(:) + PRGMLTR(:)
PA_RG(:) = PA_RG(:) - PRGMLTR(:)
PA_TH(:) = PA_TH(:) - PRGMLTR(:)*(PLSFACT(:)-PLVFACT(:))

!
END SUBROUTINE ICE4_FAST_RG
