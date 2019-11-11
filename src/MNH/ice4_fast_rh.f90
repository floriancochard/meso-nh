!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODI_ICE4_FAST_RH
INTERFACE
SUBROUTINE ICE4_FAST_RH(KSIZE, LDSOFT, LDCOMPUTE, LDWETG, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAS, PLBDAG, PLBDAR, PLBDAH, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                       &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                       &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                       &PRH_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH)
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDWETG   ! True where graupel grows in wet mode
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the rain      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAH   ! Slope parameter of the hail      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHT     ! Hail m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRGWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRGDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRDRYHG  ! Conversion of hailstone into graupel
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRHMLTR  ! Melting of the hailstones
REAL, DIMENSION(KSIZE, 8),    INTENT(INOUT) :: PRH_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
END SUBROUTINE ICE4_FAST_RH
END INTERFACE
END MODULE MODI_ICE4_FAST_RH
SUBROUTINE ICE4_FAST_RH(KSIZE, LDSOFT, LDCOMPUTE, LDWETG, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAS, PLBDAG, PLBDAR, PLBDAH, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                       &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                       &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                       &PRH_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast rh process
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
USE MODD_PARAM_ICE, ONLY : LEVLIMIT, LNULLWETH, LWETHPOST, LCONVHG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDWETG   ! True where graupel grows in wet mode
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the rain      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAH   ! Slope parameter of the hail      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHT     ! Hail m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRGWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRGDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRDRYHG  ! Conversion of hailstone into graupel
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRHMLTR  ! Melting of the hailstones
REAL, DIMENSION(KSIZE, 8),    INTENT(INOUT) :: PRH_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCWETH=1, IRRWETH=2, IRIDRYH=3, IRIWETH=4, IRSDRYH=5, IRSWETH=6, IRGDRYH=7, IRGWETH=8
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GHAIL, GWET, GMASK, LLWETH, LLDRYH
INTEGER :: IHAIL, IGWET
REAL, DIMENSION(SIZE(PRHODREF)) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(SIZE(PRHODREF)) :: IVEC1, IVEC2
REAL, DIMENSION(SIZE(PRHODREF)) :: ZZW, &
                                   ZRDRYH_INIT, ZRWETH_INIT, &
                                   ZRDRYHG
INTEGER :: JJ
!
!-------------------------------------------------------------------------------
!
!
!*       7.2    compute the Wet and Dry growth of hail
!
GMASK(:)=PRHT(:)>XRTMIN(7) .AND. PRCT(:)>XRTMIN(2) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRH_TEND(:, IRCWETH)=0.
  END WHERE
ELSE
  PRH_TEND(:, IRCWETH)=0.
  WHERE(GMASK(:))
    ZZW(:) = PLBDAH(:)**(XCXH-XDH-2.0) * PRHODREF(:)**(-XCEXVT)
    PRH_TEND(:, IRCWETH)=XFWETH * PRCT(:) * ZZW(:)    ! RCWETH
  END WHERE
ENDIF
GMASK(:)=PRHT(:)>XRTMIN(7) .AND. PRIT(:)>XRTMIN(4) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRH_TEND(:, IRIWETH)=0.
    PRH_TEND(:, IRIDRYH)=0.
  END WHERE
ELSE
  PRH_TEND(:, IRIWETH)=0.
  PRH_TEND(:, IRIDRYH)=0.
  WHERE(GMASK(:))
    ZZW(:) = PLBDAH(:)**(XCXH-XDH-2.0) * PRHODREF(:)**(-XCEXVT)
    PRH_TEND(:, IRIWETH)=XFWETH * PRIT(:) * ZZW(:)   ! RIWETH
    PRH_TEND(:, IRIDRYH)=PRH_TEND(:, IRIWETH)*(XCOLIH*EXP(XCOLEXIH*(PT(:)-XTT)))   ! RIDRYH
  END WHERE
ENDIF

!
!*       7.2.1  accretion of aggregates on the hailstones
!
GWET(:) = PRHT(:)>XRTMIN(7) .AND. PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GWET(:))
    PRH_TEND(:, IRSWETH)=0.
    PRH_TEND(:, IRSDRYH)=0.
  END WHERE
ELSE
  PRH_TEND(:, IRSWETH)=0.
  PRH_TEND(:, IRSDRYH)=0.
  IGWET=COUNT(GWET(:))
  IF(IGWET>0)THEN
    !
    !*       7.2.3  select the (PLBDAH,PLBDAS) couplet
    !
    ZVEC1(1:IGWET) = PACK( PLBDAH(:),MASK=GWET(:) )
    ZVEC2(1:IGWET) = PACK( PLBDAS(:),MASK=GWET(:) )
    !
    !*       7.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
    !               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
    !               tabulate the SWETH-kernel
    !
    ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAH)-0.00001,           &
                          XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
    !
    ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAS)-0.00001,           &
                          XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + XWETINTP2S ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
    !
    !*       7.2.5  perform the bilinear interpolation of the normalized
    !               SWETH-kernel
    !
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                  - ( XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                  - XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGWET),MASK=GWET,FIELD=0.0 )
    !
    WHERE(GWET(:))
      PRH_TEND(:, IRSWETH)=XFSWETH*ZZW(:)                       & ! RSWETH
                    *( PLBDAS(:)**(XCXS-XBS) )*( PLBDAH(:)**XCXH )  &
                       *( PRHODREF(:)**(-XCEXVT-1.) )               &
                       *( XLBSWETH1/( PLBDAH(:)**2              ) + &
                          XLBSWETH2/( PLBDAH(:)   * PLBDAS(:)   ) + &
                          XLBSWETH3/(               PLBDAS(:)**2) )
      PRH_TEND(:, IRSDRYH)=PRH_TEND(:, IRSWETH)*(XCOLSH*EXP(XCOLEXSH*(PT(:)-XTT)))
    END WHERE
  ENDIF
ENDIF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
GWET(:) = PRHT(:)>XRTMIN(7) .AND. PRGT(:)>XRTMIN(6) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GWET(:))
    PRH_TEND(:, IRGWETH)=0.
    PRH_TEND(:, IRGDRYH)=0.
  END WHERE
ELSE
  PRH_TEND(:, IRGWETH)=0.
  PRH_TEND(:, IRGDRYH)=0.
  IGWET=COUNT(GWET(:))
  IF(IGWET>0)THEN
    !
    !*       7.2.8  select the (PLBDAH,PLBDAG) couplet
    !
    ZVEC1(1:IGWET) = PACK( PLBDAH(:),MASK=GWET(:) )
    ZVEC2(1:IGWET) = PACK( PLBDAG(:),MASK=GWET(:) )
    !
    !*       7.2.9  find the next lower indice for the PLBDAH and for the PLBDAG
    !               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
    !               tabulate the GWETH-kernel
    !
    ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAG)-0.00001,           &
                          XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
    !
    ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAG)-0.00001,           &
                          XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + XWETINTP2G ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
    !
    !*       7.2.10 perform the bilinear interpolation of the normalized
    !               GWETH-kernel
    !
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                - (  XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                        * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGWET),MASK=GWET,FIELD=0.0 )
    !
    WHERE(GWET(:))
      PRH_TEND(:, IRGWETH)=XFGWETH*ZZW(:)                       & ! RGWETH
                    *( PLBDAG(:)**(XCXG-XBG) )*( PLBDAH(:)**XCXH )  &
                       *( PRHODREF(:)**(-XCEXVT-1.) )               &
                       *( XLBGWETH1/( PLBDAH(:)**2              ) + &
                          XLBGWETH2/( PLBDAH(:)   * PLBDAG(:)   ) + &
                          XLBGWETH3/(               PLBDAG(:)**2) )
      PRH_TEND(:, IRGDRYH)=PRH_TEND(:, IRGWETH)
    END WHERE
    !When graupel grows in wet mode, graupel is wet (!) and collection efficiency must remain the same
    WHERE(GWET(:) .AND. .NOT. LDWETG(:))
      PRH_TEND(:, IRGDRYH)=PRH_TEND(:, IRGDRYH)*(XCOLGH*EXP(XCOLEXGH*(PT(:)-XTT)))
    END WHERE
  END IF
ENDIF
!
!*       7.2.11  accretion of raindrops on the hailstones
!
GWET(:) = PRHT(:)>XRTMIN(7) .AND. PRRT(:)>XRTMIN(3) .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GWET(:))
    PRH_TEND(:, IRRWETH)=0.
  END WHERE
ELSE
  PRH_TEND(:, IRRWETH)=0.
  IGWET=COUNT(GWET(:))
  IF(IGWET>0)THEN
    !
    !*       7.2.12  select the (PLBDAH,PLBDAR) couplet
    !
    ZVEC1(1:IGWET)=PACK(PLBDAH(:), MASK=GWET(:))
    ZVEC2(1:IGWET)=PACK(PLBDAR(:), MASK=GWET(:))
    !
    !*       7.2.13 find the next lower indice for the PLBDAH and for the PLBDAR
    !               in the geometrical set of (Lbda_h,Lbda_r) couplet use to
    !               tabulate the RWETH-kernel
    !
    ZVEC1(1:IGWET)=MAX(1.00001, MIN( FLOAT(NWETLBDAH)-0.00001,           &
                          XWETINTP1H*LOG(ZVEC1(1:IGWET))+XWETINTP2H))
    IVEC1(1:IGWET)=INT(ZVEC1(1:IGWET))
    ZVEC1(1:IGWET)=ZVEC1(1:IGWET)-FLOAT(IVEC1(1:IGWET))
    !
    ZVEC2(1:IGWET)=MAX(1.00001, MIN( FLOAT(NWETLBDAR)-0.00001,           &
                          XWETINTP1R*LOG(ZVEC2(1:IGWET))+XWETINTP2R))
    IVEC2(1:IGWET)=INT(ZVEC2(1:IGWET))
    ZVEC2(1:IGWET)=ZVEC2(1:IGWET)-FLOAT(IVEC2(1:IGWET))
    !
    !*       7.2.14 perform the bilinear interpolation of the normalized
    !               RWETH-kernel
    !
    DO JJ=1, IGWET
      ZVEC3(JJ)= (  XKER_RWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  XKER_RWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         *(ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:)=UNPACK(VECTOR=ZVEC3(1:IGWET), MASK=GWET, FIELD=0.)
    !
    WHERE(GWET(:))
      PRH_TEND(:, IRRWETH) = XFRWETH*ZZW(:)                    & ! RRWETH
                        *( PLBDAR(:)**(-4) )*( PLBDAH(:)**XCXH ) &
                               *( PRHODREF(:)**(-XCEXVT-1.) )   &
                    *( XLBRWETH1/( PLBDAH(:)**2              ) + &
                       XLBRWETH2/( PLBDAH(:)   * PLBDAR(:)   ) + &
                       XLBRWETH3/(               PLBDAR(:)**2) )
    END WHERE
  ENDIF
ENDIF
!
ZRDRYH_INIT(:)=PRH_TEND(:, IRCWETH)+PRH_TEND(:, IRIDRYH)+PRH_TEND(:, IRSDRYH)+PRH_TEND(:, IRRWETH)+PRH_TEND(:, IRGDRYH)
!
!*       7.3    compute the Wet growth of hail
!
GHAIL(:) = PRHT(:)>XRTMIN(7) .AND. LDCOMPUTE(:)
ZRWETH_INIT(:)=0.
WHERE(GHAIL(:))
  ZRWETH_INIT(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
END WHERE
IF(LEVLIMIT) THEN
  WHERE(GHAIL(:))
    ZRWETH_INIT(:) = MIN(ZRWETH_INIT(:), EXP(XALPI-XBETAI/PT(:)-XGAMI*ALOG(PT(:)))) ! min(ev, es_i(T))
  END WHERE
ENDIF
WHERE(GHAIL(:))
  ZRWETH_INIT(:) = PKA(:)*(XTT-PT(:)) +                                 &
            ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                        *(XESTT-ZRWETH_INIT(:))/(XRV*PT(:))             )
  !
  ! compute RWETH
  !
  ZRWETH_INIT(:)  =  MAX(0.,  ( ZRWETH_INIT(:) * ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &
                            X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH ) +   &
               ( PRH_TEND(:, IRIWETH)+PRH_TEND(:, IRSWETH)+PRH_TEND(:, IRGWETH) ) *                  &
               ( PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PT(:)))   ) ) / &
                     ( PRHODREF(:)*(XLMTT-XCL*(XTT-PT(:))) ) )
  ZRWETH_INIT(:)=MAX(ZRWETH_INIT(:), PRH_TEND(:, IRIWETH)+PRH_TEND(:, IRSWETH)+PRH_TEND(:, IRGWETH))
END WHERE
!
!*       7.4    Select Wet or Dry case
!
!Wet case
LLWETH(:)=GHAIL(:) .AND. MAX(0., ZRDRYH_INIT(:)-PRH_TEND(:, IRIDRYH)-PRH_TEND(:, IRSDRYH)-PRH_TEND(:, IRGDRYH))>= &
                       & MAX(0., ZRWETH_INIT(:)-PRH_TEND(:, IRIWETH)-PRH_TEND(:, IRSWETH)-PRH_TEND(:, IRGWETH))
IF(LNULLWETH) THEN
  LLWETH(:)=LLWETH(:) .AND. ZRDRYH_INIT(:)>0.
ELSE
  LLWETH(:)=LLWETH(:) .AND. ZRWETH_INIT(:)>0.
ENDIF
IF(.NOT. LWETHPOST) LLWETH(:)=LLWETH(:) .AND. PT(:)<XTT
LLDRYH(:)=GHAIL(:) .AND. PT(:)<XTT .AND. ZRDRYH_INIT(:)>0. .AND. &
                      & MAX(0., ZRDRYH_INIT(:)-PRH_TEND(:, IRIDRYH)-PRH_TEND(:, IRSDRYH))< &
                      & MAX(0., ZRWETH_INIT(:)-PRH_TEND(:, IRIWETH)-PRH_TEND(:, IRSWETH))
!
PRCWETH(:)=0.
PRIWETH(:)=0.
PRSWETH(:)=0.
PRGWETH(:)=0.
PRRWETH(:)=0.
WHERE (LLWETH(:))
  PRCWETH(:) = PRH_TEND(:, IRCWETH)
  PRIWETH(:) = PRH_TEND(:, IRIWETH)
  PRSWETH(:) = PRH_TEND(:, IRSWETH)
  PRGWETH(:) = PRH_TEND(:, IRGWETH)
  !Collected minus aggregated
  PRRWETH(:) = ZRWETH_INIT(:) - PRH_TEND(:, IRIWETH) - PRH_TEND(:, IRSWETH) - PRH_TEND(:, IRGWETH) - PRH_TEND(:, IRCWETH)
END WHERE

PRCDRYH(:) = 0.
PRIDRYH(:) = 0.
PRSDRYH(:) = 0.
PRRDRYH(:) = 0.
PRGDRYH(:) = 0.
PRDRYHG(:) = 0.
ZRDRYHG(:)=0.
IF(LCONVHG)THEN
  WHERE(LLDRYH(:))
    ZRDRYHG(:)=ZRDRYH_INIT(:)*ZRWETH_INIT(:)/(ZRDRYH_INIT(:)+ZRWETH_INIT(:))
  END WHERE
ENDIF
WHERE(LLDRYH(:)) ! Dry
  PRCDRYH(:) = PRH_TEND(:, IRCWETH)
  PRIDRYH(:) = PRH_TEND(:, IRIDRYH)
  PRSDRYH(:) = PRH_TEND(:, IRSDRYH)
  PRRDRYH(:) = PRH_TEND(:, IRRWETH)
  PRGDRYH(:) = PRH_TEND(:, IRGDRYH)
  PRDRYHG(:) = ZRDRYHG(:)
END WHERE
PA_RC(:) = PA_RC(:) - PRCWETH(:)
PA_RI(:) = PA_RI(:) - PRIWETH(:)
PA_RS(:) = PA_RS(:) - PRSWETH(:)
PA_RG(:) = PA_RG(:) - PRGWETH(:)
PA_RH(:) = PA_RH(:) + PRCWETH(:)+PRIWETH(:)+PRSWETH(:)+PRGWETH(:)+PRRWETH
PA_RR(:) = PA_RR(:) - PRRWETH(:)
PA_TH(:) = PA_TH(:) + (PRRWETH(:)+PRCWETH(:))*(PLSFACT(:)-PLVFACT(:))
PA_RC(:) = PA_RC(:) - PRCDRYH(:)
PA_RI(:) = PA_RI(:) - PRIDRYH(:)
PA_RS(:) = PA_RS(:) - PRSDRYH(:)
PA_RR(:) = PA_RR(:) - PRRDRYH(:)
PA_RG(:) = PA_RG(:) - PRGDRYH(:) + PRDRYHG(:)
PA_RH(:) = PA_RH(:) + PRCDRYH(:)+PRIDRYH(:)+PRSDRYH(:)+PRRDRYH(:)+PRGDRYH(:) - PRDRYHG(:)
PA_TH(:) = PA_TH(:) + (PRCDRYH(:)+PRRDRYH(:))*(PLSFACT(:)-PLVFACT(:))
!
!*       7.5    Melting of the hailstones
!
GMASK(:)=PRHT(:)>XRTMIN(7) .AND. PT(:)>XTT .AND. LDCOMPUTE(:)
IF(LDSOFT) THEN
  WHERE(.NOT. GMASK(:))
    PRHMLTR(:) = 0.
  END WHERE
ELSE
  PRHMLTR(:) = 0.0
  WHERE(GMASK(:))
    PRHMLTR(:) = PRVT(:)*PPRES(:)/((XMV/XMD)+PRVT(:)) ! Vapor pressure
  END WHERE
  IF(LEVLIMIT) THEN
    WHERE(GMASK(:))
      PRHMLTR(:)=MIN(PRHMLTR(:), EXP(XALPW-XBETAW/PT(:)-XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  WHERE(GMASK(:))
    PRHMLTR(:) = PKA(:)*(XTT-PT(:)) +                              &
           ( PDV(:)*(XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT )) &
                           *(XESTT-PRHMLTR(:))/(XRV*PT(:))         )
  END WHERE
  WHERE(GMASK(:))
    !
    ! compute RHMLTR
    !
    PRHMLTR(:)  = MAX( 0.0,( -PRHMLTR(:) *                     &
                           ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &
                             X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH ) -   &
                         ( PRH_TEND(:, IRCWETH)+PRH_TEND(:, IRRWETH) )*        &
                               ( PRHODREF(:)*XCL*(XTT-PT(:))) ) /    &
                                             ( PRHODREF(:)*XLMTT ) )
  END WHERE
END IF
PA_RR(:) = PA_RR(:) + PRHMLTR(:)
PA_RH(:) = PA_RH(:) - PRHMLTR(:)
PA_TH(:) = PA_TH(:) - PRHMLTR(:)*(PLSFACT(:)-PLVFACT(:))
!
!
END SUBROUTINE ICE4_FAST_RH
