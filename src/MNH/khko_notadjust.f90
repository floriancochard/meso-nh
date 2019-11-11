!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_KHKO_NOTADJUST
!     ##########################
!
INTERFACE
!
      SUBROUTINE KHKO_NOTADJUST(KRR, KTCOUNT, TPFILE, HRAD, OCLOSE_OUT,             &
                                PTSTEP, PRHODJ, PPABSM,  PPABST, PRHODREF, PZZ,     &
                                PTHT,PRVT,PRCT,PRRT,                                &
                                PTHS, PRVS, PRCS, PRRS, PCCS, PCNUCS, PSAT,         &
                                PCLDFR, PSRCS, PNPRO,PSSPRO                          )
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KTCOUNT      ! Number of moist variables
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE   ! Output file
CHARACTER*4,              INTENT(IN)    :: HRAD     ! Radiation scheme name
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                    ! the OUTPUT FM-file
                                                    ! Condensation
REAL,                     INTENT(IN)    :: PTSTEP   ! Double Time step
                                                    ! (single if cold start)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PPABSM  ! Absolute Pressure at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PPABST  ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! Reference density
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT     ! Water vapor m.r. 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT     ! Cloud water m.r. 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT     ! Rain  water m.r. 
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS     ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS     ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS     ! Rain  water m.r. source
                                                    ! the nucleation param.
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS     ! Cloud water conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCNUCS   ! Nucl. aero. conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNPRO    ! Nucl. aero. conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSSPRO   ! Nucl. aero. conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSAT     ! Super saturation  source
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction
!
!
END SUBROUTINE KHKO_NOTADJUST
!
END INTERFACE
!
END MODULE MODI_KHKO_NOTADJUST
!
!     ################################################################################
      SUBROUTINE KHKO_NOTADJUST(KRR, KTCOUNT, TPFILE, HRAD, OCLOSE_OUT,             &
                                PTSTEP, PRHODJ, PPABSM,  PPABST, PRHODREF, PZZ,     &
                                PTHT,PRVT,PRCT,PRRT,                                &
                                PTHS, PRVS, PRCS, PRRS, PCCS, PCNUCS, PSAT,         &
                                PCLDFR, PSRCS, PNPRO,PSSPRO                          )
!     ################################################################################
!
!!****  * -  compute pseudo-prognostic of supersaturation according to Thouron
!                                                                     et al. 2012
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!
!!    REFERENCE
!!    ---------
!!
!!      Thouron, O., J.-L. Brenguier, and F. Burnet, Supersaturation calculation
!!      in large eddy simulation models for prediction of the droplet number
!!      concentration, Geosci. Model Dev., 5, 761-772, 2012.
!!
!!    AUTHOR
!!    ------
!!
!!      O. Thouron   * CNRM Meteo-France* : 
!!
!!    MODIFICATIONS
!!    -------------
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!   M.Mazoyer : 04/16 : New dummy arguments
!!   M.Mazoyer : 10/2016 New KHKO output fields
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_BUDGET
USE MODD_CONF
USE MODD_CST
USE MODD_IO_ll,           ONLY: TFILEDATA
USE MODD_LUNIT_n,         ONLY: TLUOUT
USE MODD_NSV,             ONLY: NSV_C2R2BEG
USE MODD_PARAMETERS
USE MODD_RAIN_C2R2_DESCR, ONLY: XRTMIN

!
USE MODE_FIELD,           ONLY: TFIELDDATA,TYPEREAL
USE MODE_FMWRIT
USE MODE_IO_ll
USE MODE_MSG
!
USE MODI_BUDGET
USE MODI_PROGNOS
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KTCOUNT      ! Number of moist variables
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE   ! Output file
CHARACTER*4,              INTENT(IN)    :: HRAD     ! Radiation scheme name
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT ! Conditional closure of 
                                                    ! the OUTPUT FM-file
                                                    ! Condensation
REAL,                     INTENT(IN)    :: PTSTEP   ! Double Time step
                                                    ! (single if cold start)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT     ! Water vapor m.r. 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT     ! Cloud water m.r. 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT     ! Rain  water m.r. 
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ   ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABSM   ! Absolute Pressure at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST   ! Absolute Pressure at t     
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! Reference density

REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PTHS     ! Theta source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRVS     ! Water vapor m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRCS     ! Cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRRS     ! Rain  water m.r. source
                                                    ! the nucleation param.
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCCS     ! Cloud water conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCNUCS   ! Nucl. aero. conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PNPRO    ! Nucl. aero. conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSSPRO   ! Nucl. aero. conc. source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSAT     ! pseudo-prognostic S
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PSRCS   ! Second-order flux
                                                   ! s'rc'/2Sigma_s2 at time t+1
                                                   ! multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(OUT)   :: PCLDFR  ! Cloud fraction
!
!
!*       0.2   Declarations of local variables :
!
!
!
INTEGER             :: IRESP      ! Return code of FM routines
INTEGER             :: ILUOUT     ! Logical unit of output listing 

! For Activation :                       
LOGICAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) &
            :: GNUCT  ! Test where to compute the HEN process
INTEGER , DIMENSION(SIZE(GNUCT))  :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
REAL, DIMENSION(:), ALLOCATABLE ::ZPRES,ZRHOD,ZRR,ZTT,ZRV,ZRC,ZS0,ZCCN,ZCCL, &
                                  ZZDZ, ZZLV, ZZCPH
!
INTEGER :: INUCT
INTEGER :: IIB           !  Define the domain where 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1             ! Work vectors for 
                                                        ! interpolations
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1             ! Vectors of indices for

REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3)) ::&
                       ZEXNT,ZEXNS,ZT,ZRVSAT,ZWORK,ZLV,ZCPH, ZW1,        &
                       ZACT, ZDZ
INTEGER           :: JK            ! For loop
TYPE(TFIELDDATA)  :: TZFIELD

!-------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
ILUOUT = TLUOUT%NLU
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE QUANTITIES WITH THE GUESS OF THE FUTURE INSTANT
!               -------------------------------------------------------
!
!*       2.1    remove negative non-precipitating negative water
!               ------------------------------------------------
!
IF (ANY(PRCS(:,:,:) < 0. .OR. PCCS(:,:,:) < 0.)) THEN
  WRITE(ILUOUT,*) 'C2R2_ADJUST beginning:  negative values of PRCS or PCCS'
  WRITE(ILUOUT,*) '  location of minimum of PRCS:', MINLOC(PRCS(:,:,:))
  WRITE(ILUOUT,*) ' value of minimum   :', MINVAL(PRCS(:,:,:))
  WRITE(ILUOUT,*) '  location of minimum of PCCS:', MINLOC(PCCS(:,:,:))
  WRITE(ILUOUT,*) ' value of minimum   :', MINVAL(PCCS(:,:,:))
END IF
!
IF (ANY(PRCS(:,:,:)+PRVS(:,:,:) < 0.) .AND. NVERB>5) THEN
  WRITE(ILUOUT,*) 'C2R2_ADJUST:  negative values of total water (reset to zero)'
  WRITE(ILUOUT,*) '  location of minimum:', MINLOC(PRCS(:,:,:)+PRVS(:,:,:))
  WRITE(ILUOUT,*) '  value of minimum   :', MINVAL(PRCS(:,:,:)+PRVS(:,:,:))
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','KHKO_NOTADJUST','')
END IF
!
ZACT(:,:,:)   = PCCS(:,:,:)
!
!*       2.2    estimate the Exner function at t+1 and t respectively
!
ZEXNS(:,:,:)=((2.* PPABST(:,:,:)-PPABSM(:,:,:))/XP00 )**(XRD/XCPD)
ZEXNT(:,:,:)=(PPABST(:,:,:)/XP00 )**(XRD/XCPD)
!sources terms *dt
PRRS(:,:,:)   = PRRS(:,:,:) * PTSTEP
PRVS(:,:,:)   = PRVS(:,:,:) * PTSTEP
PRCS(:,:,:)   = PRCS(:,:,:) * PTSTEP
PSAT(:,:,:)   = (PSAT(:,:,:)* PTSTEP)-1.0
PCNUCS(:,:,:) = PCNUCS(:,:,:)* PTSTEP
PCCS(:,:,:)   = PCCS(:,:,:) * PTSTEP
!state temperature at t+dt
PTHS(:,:,:)   = PTHS(:,:,:) * PTSTEP * ZEXNS(:,:,:)

!state temperature at t
ZT(:,:,:)=PTHT(:,:,:)*ZEXNT(:,:,:)
!Lv and Cph at t
ZLV(:,:,:) = XLVTT+(XCPV-XCL)*(ZT(:,:,:)-XTT)                
ZCPH(:,:,:)= XCPD+XCPV*PRVT(:,:,:)+XCL*(PRCT(:,:,:)+PRRT(:,:,:))
!dz
DO JK=1,IKE                 
  ZDZ(:,:,JK)=PZZ(:,:,JK+1)-PZZ(:,:,JK)
END DO
!
!*       2.3    compute the latent heat of vaporization Lv(T*) at t+1
!
!Removed negligible values
WHERE ( ((PRCS(:,:,:).LT.XRTMIN(2)) .AND. (PSAT(:,:,:).LT.0.0)) .OR. &
        ((PRCS(:,:,:).GT.0.0)       .AND. (PCCS(:,:,:).LE.0.0)))
 PTHS(:,:,:)  = PTHS(:,:,:)-(ZLV(:,:,:)/ZCPH(:,:,:))*PRCS(:,:,:)
 PRVS(:,:,:)  = PRVS(:,:,:)+PRCS(:,:,:)
 PRCS(:,:,:)  = 0.0
 PCNUCS(:,:,:)= 0.0
 PCCS(:,:,:)  = 0.0
ENDWHERE
!
!selection of mesh where condensation/evaportion/activation is performed
GNUCT(:,:,:) = .FALSE.
GNUCT(IIB:IIE,IJB:IJE,IKB:IKE) = PSAT(IIB:IIE,IJB:IJE,IKB:IKE)>0.0 .OR.   &
                                 PRCS(IIB:IIE,IJB:IJE,IKB:IKE)>0.0
INUCT = COUNTJV( GNUCT(:,:,:),I1(:),I2(:),I3(:))
!3D array to 1D array
!
IF( INUCT >= 1 ) THEN
  ALLOCATE(ZPRES(INUCT))
  ALLOCATE(ZRHOD(INUCT))
  ALLOCATE(ZRR(INUCT))
  ALLOCATE(ZTT(INUCT))
  ALLOCATE(ZRV(INUCT))
  ALLOCATE(ZRC(INUCT))
  ALLOCATE(ZS0(INUCT))
  ALLOCATE(ZCCN(INUCT))
  ALLOCATE(ZCCL(INUCT))
  ALLOCATE(ZZDZ(INUCT))
  ALLOCATE(ZZLV(INUCT))
  ALLOCATE(ZZCPH(INUCT))
  DO JL=1,INUCT
   ZPRES(JL) = 2. * PPABST(I1(JL),I2(JL),I3(JL)) - PPABSM(I1(JL),I2(JL),I3(JL))
   ZRHOD(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
   ZRR(JL)   = PRRS(I1(JL),I2(JL),I3(JL))
   ZTT(JL)   = PTHS(I1(JL),I2(JL),I3(JL))
   ZRV(JL)   = PRVS(I1(JL),I2(JL),I3(JL))
   ZRC(JL)   = PRCS(I1(JL),I2(JL),I3(JL))
   ZS0(JL)   = PSAT(I1(JL),I2(JL),I3(JL))
   ZCCN(JL)  = PCNUCS(I1(JL),I2(JL),I3(JL))
   ZCCL(JL)  = PCCS(I1(JL),I2(JL),I3(JL))
   ZZDZ(JL)=ZDZ(I1(JL),I2(JL),I3(JL))
   ZZLV(JL)=ZLV(I1(JL),I2(JL),I3(JL))
   ZZCPH(JL)=ZCPH(I1(JL),I2(JL),I3(JL))
  ENDDO
  !
  !Evaporation/Condensation/activation
  CALL PROGNOS(PTSTEP,ZZDZ,ZZLV,ZZCPH,ZPRES,ZRHOD,  &
                ZRR,ZTT,ZRV,ZRC,ZS0,ZCCN,ZCCL)
  !
!1D array to 3D array
  ZWORK(:,:,:)  = PCNUCS(:,:,:)
  PCNUCS(:,:,:) = UNPACK( ZCCN(:) ,MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:)  )
  ZWORK(:,:,:)  = PCCS(:,:,:)
  PCCS(:,:,:)   = UNPACK( ZCCL(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  ZWORK(:,:,:)  = PTHS(:,:,:)
  PTHS(:,:,:)   = UNPACK( ZTT(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  ZWORK(:,:,:)  = PRVS(:,:,:)
  PRVS(:,:,:)   = UNPACK( ZRV(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  ZWORK(:,:,:)  = PRCS(:,:,:)
  PRCS(:,:,:)   = UNPACK( ZRC(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  ZWORK(:,:,:)  = PSAT(:,:,:)
  PSAT(:,:,:)   = UNPACK( ZS0(:),MASK=GNUCT(:,:,:),FIELD=ZWORK(:,:,:) )
  DEALLOCATE(ZPRES)
  DEALLOCATE(ZRHOD)
  DEALLOCATE(ZRR)
  DEALLOCATE(ZTT)
  DEALLOCATE(ZRV)
  DEALLOCATE(ZRC)
  DEALLOCATE(ZS0)
  DEALLOCATE(ZCCN)
  DEALLOCATE(ZCCL)
  DEALLOCATE(ZZDZ)
!
ENDIF
!
!Computation of saturation in the meshes where there is no
!condensation/evaporation/activation
WHERE(.NOT.GNUCT(:,:,:) )
 ZRVSAT(:,:,:) = EXP(XALPW-XBETAW/PTHS(:,:,:)-XGAMW*ALOG(PTHS(:,:,:)))
 !rvsat
 ZRVSAT(:,:,:) = (XMV / XMD)*ZRVSAT(:,:,:)/((2.* PPABST(:,:,:)-PPABSM(:,:,:))-ZRVSAT(:,:,:))
 PSAT(:,:,:)   = (PRVS(:,:,:)/ZRVSAT(:,:,:))-1D0
ENDWHERE
!
!removed erroneous concentration
WHERE (PRCS(:,:,:) .EQ. 0.0)
 PCNUCS(:,:,:) = 0.0
 PCCS(:,:,:)   = 0.0
ENDWHERE
!
!source terms /dt
PRRS(:,:,:)   = PRRS(:,:,:)/PTSTEP
PTHS(:,:,:)   = PTHS(:,:,:)/PTSTEP/ZEXNS(:,:,:)
PRVS(:,:,:)   = PRVS(:,:,:)/PTSTEP
PRCS(:,:,:)   = PRCS(:,:,:)/PTSTEP
PSAT(:,:,:)   = (PSAT(:,:,:)+1.0)/PTSTEP
PCNUCS(:,:,:) = PCNUCS(:,:,:)/PTSTEP
PCCS(:,:,:)   = PCCS(:,:,:)/PTSTEP
!
IF (ANY(PRCS(:,:,:)+PRVS(:,:,:) < 0.) .AND. NVERB>5) THEN
  WRITE(*,*) 'KHKO_NOTADJUST:  negative values of total water (reset to zero)'
  WRITE(*,*) '  location of minimum:', MINLOC(PRCS(:,:,:)+PRVS(:,:,:))
  WRITE(*,*) '  value of minimum   :', MINVAL(PRCS(:,:,:)+PRVS(:,:,:))
  CALL PRINT_MSG(NVERB_FATAL,'GEN','KHKO_NOTADJUST','')
END IF
!
!number of activated aerosols
ZACT(:,:,:)=max(PCCS(:,:,:)*PTSTEP-ZACT(:,:,:)*PTSTEP,0.0)
!saturation after condensation/evaporation processes
ZWORK(:,:,:)=(PSAT(:,:,:)*PTSTEP)-1.0
!
!*              compute the cloud fraction PCLDFR
!
WHERE (PRCS(:,:,:) > 0. )
    ZW1(:,:,:)  = 1.
ELSEWHERE
    ZW1(:,:,:)  = 0. 
ENDWHERE 
IF ( SIZE(PSRCS,3) /= 0 ) THEN
    PSRCS(:,:,:) = ZW1(:,:,:) 
END IF
!
IF ( HRAD /= 'NONE' ) THEN
     PCLDFR(:,:,:) = ZW1(:,:,:)
END IF
!
  PSSPRO(:,:,:) = 0.0
  PSSPRO(:,:,:) = ZWORK(:,:,:) 
  PNPRO(:,:,:) = 0.0
  PNPRO(:,:,:) = ZACT(:,:,:) 
!
IF ( OCLOSE_OUT ) THEN
  TZFIELD%CMNHNAME   = 'SURSAT'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'SURSAT'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_NEB'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK)
  !
  TZFIELD%CMNHNAME   = 'ACT_OD'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ACT_OD'
  TZFIELD%CUNITS     = '1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_NEB'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZACT)
END IF
!
!*       7.  STORE THE BUDGET TERMS
!            ----------------------
!
!
IF (LBUDGET_RV) CALL BUDGET (PRVS(:,:,:) * PRHODJ(:,:,:),6,'COND_BU_RRV')
IF (LBUDGET_RC) CALL BUDGET (PRCS(:,:,:) * PRHODJ(:,:,:),7,'COND_BU_RRC')
IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:) * PRHODJ(:,:,:),4,'COND_BU_RTH')
IF (LBUDGET_SV) THEN
  CALL BUDGET (PCNUCS(:,:,:) * PRHODJ(:,:,:),13+(NSV_C2R2BEG-1),'CEVA_BU_RSV') ! RCN
  CALL BUDGET (PCCS(:,:,:) * PRHODJ(:,:,:),14+(NSV_C2R2BEG-1),'CEVA_BU_RSV') ! RCC
END IF
!
 CONTAINS
!!
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
END SUBROUTINE KHKO_NOTADJUST
