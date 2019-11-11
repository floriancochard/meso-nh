!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ####################
       MODULE MODI_SUBL_BLOWSNOW
!      ####################
!
INTERFACE
      SUBROUTINE SUBL_BLOWSNOW(PZZ, PRHODJ , PRHODREF, PEXNREF , PPABST,   &
                         PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PSVT,    &
                         PTHS, PRVS, PSVS,PSNWSUBL3D,PVGK)

REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t

REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:,:,:),   INTENT(IN)  :: PSVT    ! Blowing snow concentration
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PVGK    ! Mass averaged settling velocity



REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:,:),   INTENT(INOUT) :: PSVS  ! Blowing snow source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSNWSUBL3D ! Blowing snow sublimation flux (kg/m3/s)



END SUBROUTINE SUBL_BLOWSNOW
END INTERFACE
END MODULE MODI_SUBL_BLOWSNOW

      SUBROUTINE SUBL_BLOWSNOW(PZZ, PRHODJ , PRHODREF, PEXNREF , PPABST,    &
                         PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PSVT,     &
                         PTHS, PRVS, PSVS,PSNWSUBL3D,PVGK)
!!   MODIFICATIONS
!!    -------------
!!  Philippe Wautelet 28/05/2018: corrected truncated integer division (1*10**(-6) -> 1E-6)

USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CSTS_BLOWSNOW
USE MODD_BLOWSNOW

USE MODI_GAMMA
USE MODI_GAMMA_INC
USE MODI_GAMMA_INC_LOW

USE MODE_BLOWSNOW_PSD

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t

REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel/hail m.r. at t
REAL, DIMENSION(:,:,:,:),   INTENT(INOUT)  :: PSVT    ! Drifting snow concentration at t
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PVGK    ! ! Mass averaged settling velocity


REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS    ! Water vapor m.r. source
REAL, DIMENSION(:,:,:,:),   INTENT(INOUT) :: PSVS  ! Drifting snow source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSNWSUBL3D ! Drifting snow sublimation flux (kg/m3/s)
!
!*       0.2   Declarations of local variables :
!
!
INTEGER :: JN            ! Loop index for numerical integration 
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB           !
INTEGER :: IKE           !
!
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   ::  ZBETA
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3))   ::  ZT
REAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: ZW ! work array
LOGICAL,    DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))   &
                                  :: GSUBL ! Test where to compute sublimation

REAL, DIMENSION(:), ALLOCATABLE :: ZRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:), ALLOCATABLE :: ZRIT    ! Pristine ice m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRST    ! Snow/aggregate m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRGT    ! Graupel m.r. at t

REAL, DIMENSION(:), ALLOCATABLE :: ZRVS    ! Water vapor m.r. source

REAL, DIMENSION(:,:), ALLOCATABLE :: ZSVT    ! Drifting snow m.r. at t
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSVS    ! Drifting snow m.r. source


REAL, DIMENSION(:), ALLOCATABLE :: ZTHS    ! Theta source

INTEGER, DIMENSION(:), ALLOCATABLE :: NMAX    ! Maximum index for numerical integration


REAL, DIMENSION(:), ALLOCATABLE &
         ::  ZZT,      & ! Temperature
             ZRHODREF, & ! RHO Dry REFerence
             ZPRES,    & ! Pressure
             ZKA,      & ! Thermal conductivity of the air
             ZDV,      & ! Diffusivity of water vapor in the air
             ZUSI,     & ! Undersaturation over ice
             ZZW,      & ! Work array  
             ZAI,      & ! Denominator in Thorpe and Masson (66) formulation 
             ZAA,      & ! Constant in Carrier's equation for settling velocity
             ZBB,      & ! Constant in Carrier's equation for settling velocity
             ZR1,      & ! 1st limit radius in Mitchell's formulation 
             ZR2,      & ! 2nd limit radius in Mitchell's formulation 
             ZAM1,     & ! Constant in Mitchell's fall speed : v = a * R^b 
             ZAM2,     & ! Constant in Mitchell's fall speed : v = a * R^b 
             ZAM3,     & ! Constant in Mitchell's fall speed : v = a * R^b 
             ZZBETA,   & ! Scale parameter
             ZEXNREF,  & ! EXNer Pressure REFerence
             ZMU,      & ! Air kinematic viscosity
             ZZZ,      & ! Height
             ZLSFACT,  & ! L_s/(Pi_ref*C_ph)     
             ZSNWSUBL, & ! Snow sublimation rate kg.m{-3}.s{-1}
             ZVGK        ! Mass averaged settling velocity


REAL                                                       :: ZGAM,ZVEL_CARRIER,ZR,ZVEL_VENT
REAL                                                       :: ZW_M0, ZNU , ZMASS
REAL                                                       :: ZSUM_SUBL,ZNUM,ZMOB,ZTEMP
REAL                                                       :: ZDELTAR
REAL                                                       :: ZGAM_BM1,ZGAM_BM2,ZGAM_BM3
REAL                                                       :: ZR_EFF

INTEGER                           :: IMICRO
INTEGER , DIMENSION(SIZE(GSUBL)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL, JSV       ! and PACK intrinsics
LOGICAL :: LNONEFFICIENT
LOGICAL :: LSUBL_PIEKTUK
LOGICAL :: LSUBL_ALPINE3D
!
!              Initialize variables
ZDELTAR = 1e-6                   ! Bin size (m)
ZGAM     = GAMMA(XALPHA_SNOW)
ZGAM_BM1 = GAMMA(1.5*XBM1+XALPHA_SNOW+1)
ZGAM_BM2 = GAMMA(1.5*XBM2+XALPHA_SNOW+1)
ZGAM_BM3 = GAMMA(1.5*XBM3+XALPHA_SNOW+1)
!
LSUBL_PIEKTUK = .TRUE.  ! Compute sublimation according to PIEKTUK (Dery and Yau, 1999)
                        ! Use mass-averaged settling velocity as ventilation
                        ! velocity 
                        ! Save computational time compared to numerical
                        ! integration of Carrier's or Mitchell's formulation
!
LSUBL_ALPINE3D    = .FALSE.  ! Compute sublimation using the method of reprsentative
                        ! radius implemented in Alpine 3D (Groot and al, 2011)



!              Air Temperature
ZT(:,:,:) = PTHT(:,:,:) * ( PPABST(:,:,:) / XP00 ) ** (XRD/XCPD)

!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PZZ,2) - JPHEXT
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!

!
!-------------------------------------------------------------------------------
!
!*       2. COMPUTE THE BLOWINGG SNOW SUBLIMATION
!
! Optimization by looking for locations where
! the blowing snow fields are larger than a minimal value only !!!
!
!     compute parameters of the snow particle distribution
!
CALL PPP2SNOW(PSVT, PRHODREF, PBET3D=ZBETA)
!

GSUBL(:,:,:) = .FALSE.

GSUBL(IIB:IIE,IJB:IJE,IKB:IKE) =                               &
                PSVT(IIB:IIE,IJB:IJE,IKB:IKE,1)>10 .AND. &
                PSVT(IIB:IIE,IJB:IJE,IKB:IKE,2)>1e-20

!GSUBL(IIB:IIE,IJB:IJE,IKB:IKE) =                               &
!                PSVT(IIB:IIE,IJB:IJE,IKB:IKE,1)>0. .OR. &
!                PSVT(IIB:IIE,IJB:IJE,IKB:IKE,2)>0.

IMICRO = COUNTJV( GSUBL(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 0 ) THEN
  ALLOCATE(ZRVT(IMICRO))
  ALLOCATE(ZRCT(IMICRO))
  ALLOCATE(ZRRT(IMICRO))
  ALLOCATE(ZRIT(IMICRO))
  ALLOCATE(ZRST(IMICRO))
  ALLOCATE(ZRGT(IMICRO)) 
  ALLOCATE(ZRVS(IMICRO))
  ALLOCATE(ZSVT(IMICRO, NBLOWSNOW3D )) 
  ALLOCATE(ZSVS(IMICRO, NBLOWSNOW3D ))
  ALLOCATE(ZTHS(IMICRO))
  ALLOCATE(ZZT(IMICRO))
  ALLOCATE(ZRHODREF(IMICRO))
  ALLOCATE(ZPRES(IMICRO))
  ALLOCATE(ZZBETA(IMICRO))
  ALLOCATE(ZEXNREF(IMICRO))
  ALLOCATE(ZZZ(IMICRO))
  ALLOCATE(ZVGK(IMICRO))
  ALLOCATE(ZSNWSUBL(IMICRO))


  DO JL=1,IMICRO
        ZRVT(JL)     = PRVT(I1(JL),I2(JL),I3(JL))
        ZRCT(JL)     = PRCT(I1(JL),I2(JL),I3(JL))
        ZRRT(JL)     = PRRT(I1(JL),I2(JL),I3(JL))
        ZRIT(JL)     = PRIT(I1(JL),I2(JL),I3(JL))
        ZRST(JL)     = PRST(I1(JL),I2(JL),I3(JL))
        ZRGT(JL)     = PRGT(I1(JL),I2(JL),I3(JL))
        ZRVS(JL)     = PRVS(I1(JL),I2(JL),I3(JL))
        ZSVT(JL,:)   = PSVT(I1(JL),I2(JL),I3(JL),:)
        ZSVS(JL,:)   = PSVS(I1(JL),I2(JL),I3(JL),:)
        ZTHS(JL)     = PTHS(I1(JL),I2(JL),I3(JL))
        ZZT(JL)      = ZT(I1(JL),I2(JL),I3(JL))
        ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
        ZPRES(JL)    = PPABST(I1(JL),I2(JL),I3(JL))
        ZZBETA(JL)   = ZBETA(I1(JL),I2(JL),I3(JL))
        ZEXNREF(JL)  = PEXNREF(I1(JL),I2(JL),I3(JL))
        ZZZ(JL)      = PZZ(I1(JL),I2(JL),I3(JL))
        ZVGK(JL)   = PVGK(I1(JL),I2(JL),I3(JL))
        ZSNWSUBL(JL)   = PSNWSUBL3D(I1(JL),I2(JL),I3(JL))
   END DO
  ALLOCATE(ZZW(IMICRO))
  ALLOCATE(ZUSI(IMICRO))
  ZZW(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) )
  ZUSI(:) = ZRVT(:)*( ZPRES(:)-ZZW(:) ) / ( (XMV/XMD) * ZZW(:) ) - 1.0
                                                      ! Undersaturation over ice
 
  ALLOCATE(ZLSFACT(IMICRO))
    ZZW(:)  = ZEXNREF(:)*( XCPD+XCPV*ZRVT(:)+XCL*(ZRCT(:)+ZRRT(:)) &
                                    +XCI*(ZRIT(:)+ZRST(:)+ZRGT(:)) )
    ZLSFACT(:) = (XLSTT+(XCPV-XCI)*(ZZT(:)-XTT))/ZZW(:) ! L_s/(Pi_ref*C_ph)                                                     
  ALLOCATE(ZKA(IMICRO))
  ALLOCATE(ZDV(IMICRO))
  ALLOCATE(ZMU(IMICRO))
  ALLOCATE(ZAI(IMICRO))
  ALLOCATE(ZAA(IMICRO))
  ALLOCATE(ZBB(IMICRO))
  ALLOCATE(ZR1(IMICRO))
  ALLOCATE(ZR2(IMICRO))
  ALLOCATE(ZAM1(IMICRO))
  ALLOCATE(ZAM2(IMICRO))
  ALLOCATE(ZAM3(IMICRO))
  ALLOCATE(NMAX(IMICRO))


CALL SNOW_SUBL

  ZW(:,:,:) = PRVS(:,:,:)
  PRVS(:,:,:) = UNPACK( ZRVS(:),MASK=GSUBL(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PTHS(:,:,:)
  PTHS(:,:,:) = UNPACK( ZTHS(:),MASK=GSUBL(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PSVS(:,:,:,1)
  PSVS(:,:,:,1) = UNPACK( ZSVS(:,1),MASK=GSUBL(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PSVS(:,:,:,2)
  PSVS(:,:,:,2) = UNPACK( ZSVS(:,2),MASK=GSUBL(:,:,:),FIELD=ZW(:,:,:) )
!  ZW(:,:,:) = PSVS(:,:,:,3)
!  PSVS(:,:,:,3) = UNPACK( ZSVS(:,3),MASK=GSUBL(:,:,:),FIELD=ZW(:,:,:) )
  ZW(:,:,:) = PSNWSUBL3D(:,:,:)
  PSNWSUBL3D(:,:,:) = UNPACK( ZSNWSUBL(:),MASK=GSUBL(:,:,:),FIELD=ZW(:,:,:) )


  DEALLOCATE(ZRVT)
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZRIT)
  DEALLOCATE(ZRST)
  DEALLOCATE(ZRGT) 
  DEALLOCATE(ZRVS)
  DEALLOCATE(ZSVT) 
  DEALLOCATE(ZSVS)
  DEALLOCATE(ZTHS)
  DEALLOCATE(ZZT)
  DEALLOCATE(ZRHODREF)
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZKA)
  DEALLOCATE(ZDV)
  DEALLOCATE(ZUSI)
  DEALLOCATE(ZZW)
  DEALLOCATE(ZAI)
  DEALLOCATE(ZAA)
  DEALLOCATE(ZBB)
  DEALLOCATE(ZR1)
  DEALLOCATE(ZR2)
  DEALLOCATE(ZAM1)
  DEALLOCATE(ZAM2)
  DEALLOCATE(ZAM3)
  DEALLOCATE(ZZBETA)
  DEALLOCATE(NMAX)
  DEALLOCATE(ZEXNREF)
  DEALLOCATE(ZLSFACT)
  DEALLOCATE(ZZZ)
  DEALLOCATE(ZMU)
  DEALLOCATE(ZSNWSUBL)
  DEALLOCATE(ZVGK)


END IF
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!
CONTAINS

!
!-------------------------------------------------------------------------------
!

SUBROUTINE SNOW_SUBL

IMPLICIT NONE

! Sutherland's equation for kinematic viscosity
ZMU(:)=1.8325d-5*416.16/(ZZT(:)+120)*(ZZT(:)/296.16)*SQRT(ZZT(:)/296.16)/ZRHODREF(:)
! Thermal conductivity of the air
  ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
! Diffusivity of water vapor in the air. 
  ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
!
!*       Compute the denominator in the Thorpe and Masson (66) equation
!  
ZAI(:) = EXP( XALPI - XBETAI/ZZT(:) - XGAMI*ALOG(ZZT(:) ) ) ! es_i
ZAI(:) = ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) ) / (ZKA(:)*ZZT(:))                  &
             *( ( XLSTT + (XCPV-XCI)*(ZZT(:)-XTT) ) / (XRV*ZZT(:)) - 1.)  &
             + (XRV*ZZT(:)) / (ZDV(:)*ZAI(:))

IF(LSUBL_ALPINE3D) THEN
         ZR_EFF = 73.5e-6  ! Effective radisus computed following the Swiss
                          ! method. This effective radius give the same total 
                          ! sublimation for a equal concentration an ensemble of
                          ! gamma distributed particles with rm = 35e-6 m and
                          ! alpha=3 
! Compute coefficient for settling velocity following Carrier (1953)
        ZAA(:) = 6.203*ZMU(:)/2.
        ZBB(:) = 5.516*XRHOLI/(4.*ZRHODREF(:))*XG
       DO JL=1,IMICRO
          ZSUM_SUBL = 0.
          ZUSI(JL) = MIN(ZUSI(JL), 0.)   !Only the undersaturation over ice is considered.
          ! Ventilation velocity taken as settling velocity of particle of mean
          ! radius
          ZVEL_VENT = - ZAA(JL)/ZR_EFF+((ZAA(JL)/ZR_EFF)**2+ZBB(JL)*ZR_EFF)**0.5
!       Nusselt Number                                                 
            ZNU    =    NUSSELT(ZR_EFF,ZMU(JL),ZVEL_VENT)
!       Rate of change of mass for a subliming ice sphere of radius ZR_EFF            
            ZMASS  =    2*XPI*ZR_EFF*ZNU*ZUSI(JL)/ZAI(JL)
!       Integration over the radius spectrum 
        ZSUM_SUBL = ZMASS*ZSVT(JL,2)/(4./3.*XPI*XRHOLI*ZR_EFF**2)

    ZSUM_SUBL =  MIN( ZRVS(JL),ZSUM_SUBL)*(0.5+SIGN(0.5,ZSUM_SUBL)) &
                   - MIN(ZSVS(JL,2),ABS(ZSUM_SUBL))*(0.5-SIGN(0.5,ZSUM_SUBL))
    ZSUM_SUBL=MIN(0.,ZSUM_SUBL) ! Sink of snow               
!   Change in concentration rate  Sn = Sb*N/qb (Dery and Yau,2000)
    ZNUM      = ZSUM_SUBL*ZSVT(JL,1)/ZSVT(JL,2)
!   Change in mobility index : value mob*M3 is reduced according to reduction of
!   M3 due to sublimation so that mob is constant due to sublimation
!    ZMOB      = ZSUM_SUBL*ZSVT(JL,3)/ZSVT(JL,2)
!   Update tendencies for snow particles, water vapour and potential temperature 
    ZSVS(JL,2) = ZSVS(JL,2) + ZSUM_SUBL             ! Particle mixing ratio
    ZSVS(JL,1) = ZSVS(JL,1) + ZNUM                  ! Particles number
!    ZSVS(JL,3) = ZSVS(JL,3) + ZMOB
    ZRVS(JL) = ZRVS(JL) - ZSUM_SUBL                 ! Water vapour
    ZTHS(JL) = ZTHS(JL) + ZSUM_SUBL*ZLSFACT(JL)     ! Potential temperature
    ZSNWSUBL(JL) = ZSNWSUBL(JL)+ZSUM_SUBL*ZRHODREF(JL)  ! Sublimation rate kg/m3/s

END DO
ELSE IF(LSUBL_PIEKTUK) THEN 
DO JL=1,IMICRO
   ZSUM_SUBL=0.
   ZUSI(JL) = MIN(ZUSI(JL), 0.)   !Only the undersaturation over ice is considered.
!       Ventilation velocity as mass averaged settling velocity            
   ZVEL_VENT = ZVGK(JL)      
!       Nusselt Number using mean radius of particle size distribution  
   ZNU    =    NUSSELT(XALPHA_SNOW*ZZBETA(JL),ZMU(JL),ZVEL_VENT)
   ! mass averaged sublimation rate follows Dery and Yan (1999) and avoids
   ! numerical integration over the particle spectrum
   ZSUM_SUBL = ZSVT(JL,2)*ZNU*ZUSI(JL)/   &
                (ZAI(JL)*2*XRHOLI*(XALPHA_SNOW*ZZBETA(JL))**2)
!   Restriction of ZSUM_SUBL   
    ZTEMP=ZSUM_SUBL
    ZSUM_SUBL =  MIN( ZRVS(JL),ZSUM_SUBL)*(0.5+SIGN(0.5,ZSUM_SUBL)) &
                   - MIN(ZSVS(JL,2),ABS(ZSUM_SUBL))*(0.5-SIGN(0.5,ZSUM_SUBL))
    ZSUM_SUBL=MIN(0.,ZSUM_SUBL) ! Sink of snow               
    IF(ZSUM_SUBL>0) THEN  
       write(*,*) 'Warning Subl',JL,'Subl',ZSUM_SUBL,'TEMP',ZTEMP
       write(*,*) 'Warning Subl ZSVT',ZSVT(JL,1),ZSVT(JL,2)
       write(*,*) 'Warning vap',ZRVS(JL),'ZSVs',ZSVS(JL,2)
    END IF
!   Change in concentration rate  Sn = Sb*N/qb (Dery and Yau,2000)
    ZNUM      = ZSUM_SUBL*ZSVT(JL,1)/ZSVT(JL,2)
!   Change in mobility index : value mob*M3 is reduced according to reduction of
!   M3 due to sublimation so that mob is constant due to sublimation
!    ZMOB      = ZSUM_SUBL*ZSVT(JL,3)/ZSVT(JL,2)
!   Update tendencies for snow particles, water vapour and potential temperature 
    ZSVS(JL,2) = ZSVS(JL,2) + ZSUM_SUBL             ! Particle mixing ratio
    ZSVS(JL,1) = ZSVS(JL,1) + ZNUM                  ! Particles number
!    ZSVS(JL,3) = ZSVS(JL,3) + ZMOB
    ZRVS(JL) = ZRVS(JL) - ZSUM_SUBL                 ! Water vapour
    ZTHS(JL) = ZTHS(JL) + ZSUM_SUBL*ZLSFACT(JL)     ! Potential temperature
    ZSNWSUBL(JL) = ZSNWSUBL(JL)+ZSUM_SUBL*ZRHODREF(JL)  ! Sublimation rate kg/m3/s

END DO

ELSE        
 !
!*       Compute the constants in Carrier equation
!
IF(CSNOWSEDIM=='CARR') THEN
ZAA(:) = 6.203*ZMU(:)/2.
ZBB(:) = 5.516*XRHOLI/(4.*ZRHODREF(:))*XG
NMAX=GET_INDEX(ZZBETA(:),ZDELTAR)

DO JL=1,IMICRO
   ZSUM_SUBL=0.
   ZUSI(JL) = MIN(ZUSI(JL), 0.)   !Only the undersaturation over ice is considered.  
   DO JN=1,NMAX(JL)
            ZR = 1E-6+(JN-0.5)*ZDELTAR
!       Carrier settling velocity
            ZVEL_CARRIER = - ZAA(JL)/ZR+((ZAA(JL)/ZR)**2+ZBB(JL)*ZR)**0.5
!       Weight of the corresponding bin following the gamma distribution            
            ZW_M0=ZSVT(JL,1)*ZR**(XALPHA_SNOW-1)*exp(-ZR/ZZBETA(JL))/(ZZBETA(JL))**XALPHA_SNOW*ZGAM
!       Ventilation velocity as a sum of settling velocity and relative
!       turbulent velocity fluctuations            
            ZVEL_VENT = ZVEL_CARRIER!+TURB_FLUC(ZR,ZMU(JL),ZVEL_CARRIER,ZRHODREF(JL),    & 
!                                                 ZZZ(JL),ZVMOD(JL))
!       Nusselt Number                                                 
            ZNU    =    NUSSELT(ZR,ZMU(JL),ZVEL_VENT)
!       Rate of change of mass for a subliming ice sphere            
            ZMASS  =    2*XPI*ZR*ZNU*ZUSI(JL)/ZAI(JL)
!       Integration over the radius spectrum            
            ZSUM_SUBL = ZSUM_SUBL+ZW_M0*ZMASS*ZDELTAR 
    END DO
!   Restriction of ZSUM_SUBL    
    ZSUM_SUBL =  MIN( ZRVS(JL),ZSUM_SUBL)*(0.5+SIGN(0.5,ZSUM_SUBL)) &
                   - MIN( ZSVS(JL,2),ABS(ZSUM_SUBL) )*(0.5-SIGN(0.5,ZSUM_SUBL))
!   Change in concentration rate  Sn = Sb*N/qb (Dery and Yau,2000)
    ZNUM      = ZSUM_SUBL*ZSVT(JL,1)/ZSVT(JL,2)
!   Change in mobility index : value mob*M3 is reduced according to reduction of
!   M3 due to sublimation so that mob is constant due to sublimation
!    ZMOB      = ZSUM_SUBL*ZSVT(JL,3)/ZSVT(JL,2)
!   Update tendencies for snow particles, water vapour and potential temperature 
    ZSVS(JL,2) = ZSVS(JL,2) + ZSUM_SUBL             ! Particle mixing ratio
    ZSVS(JL,1) = ZSVS(JL,1) + ZNUM                  ! Particles number
!    ZSVS(JL,3) = ZSVS(JL,3) + ZMOB
    ZRVS(JL) = ZRVS(JL) - ZSUM_SUBL                 ! Water vapour
    ZTHS(JL) = ZTHS(JL) + ZSUM_SUBL*ZLSFACT(JL)     ! Potential temperature
    ZSNWSUBL(JL) = ZSUM_SUBL*ZRHODREF(JL)           ! Sublimation rate kg/m3/s
END DO
END IF

IF(CSNOWSEDIM=='MITC') THEN
LNONEFFICIENT = .FALSE.
!        write(*,*) 'MITC'
  ! Compute limit radius for integration of Mitchell's formulation
ZR1(:)=RLIM(ZMU,ZRHODREF,XBESTL_1)
ZR2(:)=RLIM(ZMU,ZRHODREF,XBESTL_2)
  ! Compute parameter avr for integration of Mitchell's formulation
ZAM1(:)=AVR(XAM1,XBM1,ZRHODREF,ZMU)
ZAM2(:)=AVR(XAM2,XBM2,ZRHODREF,ZMU)
ZAM3(:)=AVR(XAM3,XBM3,ZRHODREF,ZMU)

DO JL=1,IMICRO
    ZUSI(JL) = MIN(ZUSI(JL), 0.)  !Only the undersaturation over ice is considered.
!                               no water deposition on blown snow particles
IF(LNONEFFICIENT) THEN
    ZSUM_SUBL = 2*XPI*ZUSI(JL)*ZSVT(JL,1)/ZAI(JL)*(XANU*ZZBETA(JL)*XALPHA_SNOW  +      &
                 XBNU/ZGAM*(2/ZMU(JL))**0.5*(                                          & 
                 ZZBETA(JL)**(1.5*XBM1+1)*ZAM1(JL)**0.5*ZGAM_BM1*                        &
                       GAMMA_INC(1.5*XBM1+XALPHA_SNOW+1,ZR1(JL)/ZZBETA(JL)) +            &
                 ZZBETA(JL)**(1.5*XBM2+1)*ZAM2(JL)**0.5*ZGAM_BM2*                        &
                       (GAMMA_INC(1.5*XBM2+XALPHA_SNOW+1,ZR2(JL)/ZZBETA(JL))-            &
                       GAMMA_INC(1.5*XBM2+XALPHA_SNOW+1,ZR1(JL)/ZZBETA(JL)))+            & 
                 ZZBETA(JL)**(1.5*XBM3+1)*ZAM3(JL)**0.5*ZGAM_BM3*                        & 
                   (1.-GAMMA_INC(1.5*XBM3+XALPHA_SNOW+1,ZR2(JL)/ZZBETA(JL)))))
ELSE
    ZSUM_SUBL = 2*XPI*ZUSI(JL)*ZSVT(JL,1)/ZAI(JL)*(XANU*ZZBETA(JL)*XALPHA_SNOW  +      &
                 XBNU/ZGAM*(2/ZMU(JL))**0.5*(                                          & 
                 ZZBETA(JL)**(1.5*XBM1+1)*ZAM1(JL)**0.5*                        &
                       GAMMA_INC_LOW(1.5*XBM1+XALPHA_SNOW+1,ZR1(JL)/ZZBETA(JL)) +            &
                 ZZBETA(JL)**(1.5*XBM2+1)*ZAM2(JL)**0.5*                        &
                       (GAMMA_INC_LOW(1.5*XBM2+XALPHA_SNOW+1,ZR2(JL)/ZZBETA(JL))-            &
                       GAMMA_INC_LOW(1.5*XBM2+XALPHA_SNOW+1,ZR1(JL)/ZZBETA(JL)))+            & 
                 ZZBETA(JL)**(1.5*XBM3+1)*ZAM3(JL)**0.5*                        & 
                   (ZGAM_BM3-GAMMA_INC_LOW(1.5*XBM3+XALPHA_SNOW+1,ZR2(JL)/ZZBETA(JL)))))
END IF
!   Restriction of ZSUM_SUBL   
    ZTEMP=ZSUM_SUBL
    ZSUM_SUBL =  MIN( ZRVS(JL),ZSUM_SUBL)*(0.5+SIGN(0.5,ZSUM_SUBL)) &
                   - MIN(ZSVS(JL,2),ABS(ZSUM_SUBL))*(0.5-SIGN(0.5,ZSUM_SUBL))
    ZSUM_SUBL=MIN(0.,ZSUM_SUBL) ! Sink of snow               
    IF(ZSUM_SUBL>0) THEN  
       write(*,*) 'Warning Subl',JL,'Subl',ZSUM_SUBL,'TEMP',ZTEMP
       write(*,*) 'Warning Subl ZSVT',ZSVT(JL,1),ZSVT(JL,2)
       write(*,*) 'Warning vap',ZRVS(JL),'ZSVs',ZSVS(JL,2)
    END IF
!   Change in concentration rate  Sn = Sb*N/qb (Dery and Yau,2000)
    ZNUM      = ZSUM_SUBL*ZSVT(JL,1)/ZSVT(JL,2)
!   Change in mobility index : value mob*M3 is reduced according to reduction of
!   M3 due to sublimation so that mob is constant due to sublimation
!    ZMOB      = ZSUM_SUBL*ZSVT(JL,3)/ZSVT(JL,2)
!   Update tendencies for snow particles, water vapour and potential temperature 
    ZSVS(JL,2) = ZSVS(JL,2) + ZSUM_SUBL             ! Particle mixing ratio
    ZSVS(JL,1) = ZSVS(JL,1) + ZNUM                  ! Particles number
!    ZSVS(JL,3) = ZSVS(JL,3) + ZMOB
    ZRVS(JL) = ZRVS(JL) - ZSUM_SUBL                 ! Water vapour
    ZTHS(JL) = ZTHS(JL) + ZSUM_SUBL*ZLSFACT(JL)     ! Potential temperature
    ZSNWSUBL(JL) = ZSUM_SUBL*ZRHODREF(JL)           ! Sublimation rate kg/m3/s
END DO


END IF
END IF

END SUBROUTINE SNOW_SUBL
!
!-------------------------------------------------------------------------------
!
FUNCTION GET_INDEX(PBETA,PDELTAR) RESULT(KMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the upper index in numerical integration of Carrier's formulation 
!     Index equals to 5* mean radius
!
!
USE MODD_BLOWSNOW,     ONLY : XALPHA_SNOW


!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN)                                  :: PDELTAR      ! (-)
REAL, DIMENSION(:), INTENT(IN)                                  :: PBETA ! (kg/m3)

!
INTEGER, DIMENSION(SIZE(PBETA,1)) :: KMAX ! (-)
!

KMAX(:)=int(PBETA(:)*XALPHA_SNOW*5/PDELTAR)


END FUNCTION GET_INDEX
!
!-------------------------------------------------------------------------------
!
FUNCTION RLIM(PMU,PRHODREF,PBEST_LIM) RESULT(PRLIM)
!
!!    PURPOSE
!!    -------
!     Calculate the radius of a sperical particle for a given Best Number 
!
!
USE MODD_CSTS_BLOWSNOW,     ONLY : XRHOLI,XG
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)                                  :: PRHODREF ! (kg/m3)
REAL, DIMENSION(:), INTENT(IN)                                  :: PMU      ! (m2/s)
REAL,               INTENT(IN)                                  :: PBEST_LIM! (-)

!
REAL, DIMENSION(SIZE(PMU,1)) :: PRLIM ! (m)
!
PRLIM(:)=(3./32.*PRHODREF(:)/(XRHOLI*XG)*PMU(:)**2.*PBEST_LIM)**0.333333333

END FUNCTION RLIM

FUNCTION AVR(PARE,PBRE,PRHODREF,PMU) RESULT(PAVR)
!
!!    PURPOSE
!!    -------
!     Calculate the parameter av_r in KC02 formulation (Eq. 3.1)
!
!
USE MODD_CSTS_BLOWSNOW,     ONLY : XRHOLI,XG


!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN)                                  :: PARE      ! (-)
REAL,               INTENT(IN)                                  :: PBRE      ! (-)
REAL, DIMENSION(:), INTENT(IN)                                  :: PRHODREF ! (kg/m3)
REAL, DIMENSION(:), INTENT(IN)                                  :: PMU      ! (m2/s)

!
REAL, DIMENSION(SIZE(PMU,1)) :: PAVR ! (-)
!


PAVR(:)=2.**(3.*PBRE-1.)*PARE*PMU(:)**(1.-2.*PBRE)*(4./3.*XRHOLI/PRHODREF(:)*XG)**PBRE

END FUNCTION AVR
!
!-------------------------------------------------------------------------------
!
FUNCTION TURB_FLUC(PR,PMU,PCARRIER,PRHODREF,PZZ,PVMOD) RESULT(PSIG)
!
!!    PURPOSE
!!    -------
!     Calculate the relative turbulent velocity fluctuations for a given radius. 
!     Used to compute the ventilation velocity.
!     Formulation based on Dover (1993)
!
USE MODD_CSTS
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                                  :: PR ! (m)
REAL, INTENT(IN)                                  :: PMU ! (m2/s)
REAL, INTENT(IN)                                  :: PCARRIER ! (m/s)
REAL, INTENT(IN)                                  :: PRHODREF ! (kg/m3)
REAL, INTENT(IN)                                  :: PZZ ! (m)
REAL, INTENT(IN)                                  :: PVMOD ! (m/s)
!
REAL              :: PSIG ! (m/s)
!
!
!*      0.2    declaration of local variables
!
REAL                  :: ZFCRI1,ZFCRI2,ZFCRI
REAL                  :: ZS0,ZSIGU,ZSIGV,ZSIGW,ZUSTAR
!
!
!*      1    Calculate critical frequency
!
ZFCRI1 = 9*PRHODREF*PMU/(4*XPI*PR**2*XRHOLI)
ZFCRI2 = 0.363*PRHODREF*PCARRIER/(XPI*PR*XRHOLI)
ZFCRI = MAX(ZFCRI1,ZFCRI2) 
!
!*      2    Calculate variances of the horizontal and vertical velocity components
!
ZS0   = ZFCRI*PZZ/PVMOD
ZSIGU = 4.77 *ZUSTAR**2/ (1+33*ZS0)**0.66666
ZSIGV = 2.76 *ZUSTAR**2/ (1+9.5*ZS0)**0.66666
ZSIGW = 1.31 *ZUSTAR**2/ (1+3.12*ZS0)**0.66666

PSIG = (ZSIGU+ZSIGV+ZSIGW)**0.5
END FUNCTION TURB_FLUC
!
FUNCTION NUSSELT(PR,PMU,PVEL_VENT) RESULT(PNU)
!
!!    PURPOSE
!!    -------
!     Calculate the Nusselt number for a given particle radius 
!     Formulation based on Lee (1975)
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                                  :: PR ! (m)
REAL, INTENT(IN)                                  :: PMU ! (m2/s)
REAL, INTENT(IN)                                  :: PVEL_VENT ! (m/s)
!
REAL              :: PNU ! (m/s)
!
!
!*      0.2    declaration of local variables
!
REAL                  :: ZRE
!
!
!*      1    Calculate Reynolds number 
!
ZRE = 2*PR*PVEL_VENT/(PMU)
!
!*      2    Calculate Nusselt number
!
IF(ZRE<10) THEN
        PNU = 1.79+0.606*ZRE**0.5
ELSE
        PNU = 1.88+0.580*ZRE**0.5
END IF

END FUNCTION NUSSELT
!
!-------------------------------------------------------------------------------
!

!
!-------------------------------------------------------------------------------
!

 FUNCTION COUNTJV(LTAB,I1,I2,I3) RESULT(IC)
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
!
LOGICAL, DIMENSION(:,:,:) :: LTAB ! Mask
INTEGER, DIMENSION(:) :: I1,I2,I3 ! Used to replace the COUNT and PACK
INTEGER :: JI,JJ,JK,IC
!
!-------------------------------------------------------------------------------
!
IC = 0
DO JK = 1,SIZE(LTAB,3)
  DO JJ = 1,SIZE(LTAB,2)
    DO JI = 1,SIZE(LTAB,1)
      IF( LTAB(JI,JJ,JK) ) THEN
        IC = IC +1
        I1(IC) = JI
        I2(IC) = JJ
        I3(IC) = JK
      END IF
    END DO
  END DO
END DO
!
END FUNCTION COUNTJV
!
!-------------------------------------------------------------------------------
END SUBROUTINE SUBL_BLOWSNOW



