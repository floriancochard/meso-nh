!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:17
!-----------------------------------------------------------------
!     #####################
      MODULE MODI_CALCSOUND
!     #####################
INTERFACE
        SUBROUTINE CALCSOUND(PPRESS,PTEMPE,PRV,   &
                             PCAPEP,PCINP,PDCAPE,PCAPEPMAX,PCINPMAX )
!
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPRESS  ! Pressure in hPa
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTEMPE  ! Temperature in Celcius
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRV     ! vapour mixing ratio in g/kg
REAL, DIMENSION(:,:,:), INTENT(OUT) ::PCAPEP
REAL, DIMENSION(:,:,:), INTENT(OUT) ::PCINP
REAL, DIMENSION(:,:,:), INTENT(OUT) ::PDCAPE
REAL, DIMENSION(:,:),   INTENT(OUT) ::PCAPEPMAX
REAL, DIMENSION(:,:),   INTENT(OUT) ::PCINPMAX
!
END SUBROUTINE CALCSOUND
END INTERFACE
END MODULE MODI_CALCSOUND
!       #############################################################
        SUBROUTINE CALCSOUND(PPRESS,PTEMPE,PRV,   &
                             PCAPEP,PCINP,PDCAPE,PCAPEPMAX,PCINPMAX )
!!                             PZVKE,PZFCL,PZFCLMAX,PZVKEMAX,PALTMAX)
!       #############################################################
!
!!****
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to calculate various properties of
!       samples of air raised or lowered to different levels.
!!
!!**  METHOD
!!    ------
!!        The horizontal dimensions of model arrays are splitted in arrays of
!!      1000 columns. If there is at least 1000 elements, computation is
!!      made in a static way, otherwise in a dynamical way.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      C. Lac, V. Ducrocq  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!     Original  from K. Emanuel
!!     J. Stein  Jan. 2001  optimisation by splitting arrays in 1000 columns
!!     C.Lac     May 2015 correction in downdraft loop 
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PPRESS  ! Pressure in hPa
REAL, DIMENSION(:,:,:), INTENT(IN) :: PTEMPE  ! Temperature in Celcius
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRV     ! Vapor mixing ratio in g/kg
REAL, DIMENSION(:,:,:), INTENT(OUT) ::PCAPEP
REAL, DIMENSION(:,:,:), INTENT(OUT) ::PCINP
REAL, DIMENSION(:,:,:), INTENT(OUT) ::PDCAPE
REAL, DIMENSION(:,:),   INTENT(OUT) ::PCAPEPMAX
REAL, DIMENSION(:,:),   INTENT(OUT) ::PCINPMAX
!
INTEGER, PARAMETER :: NPARAM=1000
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PPRESS,1)*SIZE(PPRESS,2),SIZE(PPRESS,3)) :: ZPRESS2D, &
                    ZTEMPE2D,ZRV2D,ZCAPEP2D,ZCINP2D,ZDCAPE2D
REAL, DIMENSION(SIZE(PPRESS,1)*SIZE(PPRESS,2)) :: ZCAPEMAX2D,ZCINMAX2D
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZP,ZT,ZR,ZCAPEP,ZCINP,ZDCAPE
REAL, ALLOCATABLE, DIMENSION(:)     :: ZCAPEMAX,ZCINMAX
INTEGER :: IIU,IJU,IKU, IDIM,JLOOP
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATIONS
!              ---------------
!
!
IIU=SIZE(PPRESS,1)
IJU=SIZE(PPRESS,2)
IKU=SIZE(PPRESS,3)
!
CALL TRANSF_3D_2D(PPRESS,ZPRESS2D)
CALL TRANSF_3D_2D(PTEMPE,ZTEMPE2D)
CALL TRANSF_3D_2D(PRV,ZRV2D)
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTATION
!              -----------
!
! loops of NPARAM points
!
DO JLOOP=1,1000000
  IF (NPARAM*JLOOP< IIU*IJU) THEN
    !
    IF (.NOT. ALLOCATED(ZP)) THEN
      ALLOCATE(ZP(NPARAM,IKU))
      ALLOCATE(ZT(NPARAM,IKU))
      ALLOCATE(ZR(NPARAM,IKU))
      ALLOCATE(ZCAPEP(NPARAM,IKU))
      ALLOCATE(ZCINP(NPARAM,IKU))
      ALLOCATE(ZDCAPE(NPARAM,IKU))
      ALLOCATE(ZCAPEMAX(NPARAM))
      ALLOCATE(ZCINMAX(NPARAM))
    ENDIF
    !
    ZP(1:NPARAM,:)=ZPRESS2D((JLOOP-1)*NPARAM +1:JLOOP*NPARAM,:)
    ZT(1:NPARAM,:)=ZTEMPE2D((JLOOP-1)*NPARAM +1:JLOOP*NPARAM,:)
    ZR(1:NPARAM,:)=ZRV2D((JLOOP-1)*NPARAM +1:JLOOP*NPARAM,:)
    !
    CALL CALC(NPARAM,IKU,ZP,ZT,ZR,ZCAPEP,ZCINP,ZDCAPE,ZCAPEMAX,ZCINMAX)
    !
    ZCAPEP2D((JLOOP-1)*NPARAM +1:JLOOP*NPARAM,:) = ZCAPEP(1:NPARAM,:)
    ZCINP2D ((JLOOP-1)*NPARAM +1:JLOOP*NPARAM,:) = ZCINP(1:NPARAM,:)
    ZDCAPE2D((JLOOP-1)*NPARAM +1:JLOOP*NPARAM,:) = ZDCAPE(1:NPARAM,:)
    ZCAPEMAX2D((JLOOP-1)*NPARAM +1:JLOOP*NPARAM) = ZCAPEMAX(1:NPARAM)
    ZCINMAX2D ((JLOOP-1)*NPARAM +1:JLOOP*NPARAM) = ZCINMAX(1:NPARAM)
    !
  ELSE
    IF (ALLOCATED(ZP)) THEN
      DEALLOCATE(ZP)
      DEALLOCATE(ZT)
      DEALLOCATE(ZR)
      DEALLOCATE(ZCAPEP)
      DEALLOCATE(ZCINP)
      DEALLOCATE(ZDCAPE)
      DEALLOCATE(ZCAPEMAX)
      DEALLOCATE(ZCINMAX)
    ENDIF
    !
    IDIM=IIU*IJU-NPARAM*(JLOOP-1)
    ALLOCATE(ZP(1:IDIM,IKU))
    ALLOCATE(ZT(1:IDIM,IKU))
    ALLOCATE(ZR(1:IDIM,IKU))
    ALLOCATE(ZCAPEP(1:IDIM,IKU))
    ALLOCATE(ZCINP(1:IDIM,IKU))
    ALLOCATE(ZDCAPE(1:IDIM,IKU))
    ALLOCATE(ZCAPEMAX(1:IDIM))
    ALLOCATE(ZCINMAX(1:IDIM))
    !
    ZP(1:IDIM,:)=ZPRESS2D(NPARAM*(JLOOP-1)+1:IIU*IJU,:)
    ZT(1:IDIM,:)=ZTEMPE2D(NPARAM*(JLOOP-1)+1:IIU*IJU,:)
    ZR(1:IDIM,:)=ZRV2D(NPARAM*(JLOOP-1)+1:IIU*IJU,:)
    !
    CALL CALC(IDIM,IKU,ZP,ZT,ZR,ZCAPEP,ZCINP,ZDCAPE,ZCAPEMAX,ZCINMAX)
    !
    ZCAPEP2D(NPARAM*(JLOOP-1)+1:IIU*IJU,:) = ZCAPEP(1:IDIM,:)
    ZCINP2D (NPARAM*(JLOOP-1)+1:IIU*IJU,:) = ZCINP(1:IDIM,:)
    ZDCAPE2D(NPARAM*(JLOOP-1)+1:IIU*IJU,:) = ZDCAPE(1:IDIM,:)
    ZCAPEMAX2D(NPARAM*(JLOOP-1)+1:IIU*IJU) = ZCAPEMAX(1:IDIM)
    ZCINMAX2D (NPARAM*(JLOOP-1)+1:IIU*IJU) = ZCINMAX(1:IDIM)
    !
    DEALLOCATE(ZP)
    DEALLOCATE(ZT)
    DEALLOCATE(ZR)
    DEALLOCATE(ZCAPEP)
    DEALLOCATE(ZCINP)
    DEALLOCATE(ZDCAPE)
    DEALLOCATE(ZCAPEMAX)
    DEALLOCATE(ZCINMAX)
    !
    EXIT
    !
  ENDIF
ENDDO
!
!  back to 3D and 2D arrays
!
CALL TRANSF_2D_3D(ZCAPEP2D,PCAPEP)
CALL TRANSF_2D_3D(ZCINP2D,PCINP)
CALL TRANSF_2D_3D(ZDCAPE2D,PDCAPE)
CALL TRANSF_1D_2D(ZCAPEMAX2D,PCAPEPMAX)
CALL TRANSF_1D_2D(ZCINMAX2D,PCINPMAX)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
SUBROUTINE TRANSF_3D_2D(PTAB3D,PTAB2D)
REAL, DIMENSION (:,:)    :: PTAB2D
REAL, DIMENSION (:,:,:)  :: PTAB3D
!
INTEGER                  :: JIJ,JI,JJ,JK
!
DO JK=1,SIZE(PTAB3D,3)
  JIJ = 0
  DO JJ=1,SIZE(PTAB3D,2)
    DO JI=1,SIZE(PTAB3D,1)
     JIJ = JIJ + 1
     PTAB2D(JIJ,JK) = PTAB3D(JI,JJ,JK)
    ENDDO
  ENDDO
ENDDO
!
END SUBROUTINE TRANSF_3D_2D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE TRANSF_2D_3D(PTAB2D,PTAB3D)
REAL, DIMENSION (:,:)    :: PTAB2D
REAL, DIMENSION (:,:,:)  :: PTAB3D
!
INTEGER                  :: JIJ,JI,JJ,JK
!
DO JK=1,SIZE(PTAB3D,3)
  JIJ = 0
  DO JJ=1,SIZE(PTAB3D,2)
    DO JI=1,SIZE(PTAB3D,1)
      JIJ = JIJ + 1
      PTAB3D(JI,JJ,JK) = PTAB2D(JIJ,JK)
    ENDDO
  ENDDO
ENDDO
!
END SUBROUTINE TRANSF_2D_3D
!
!-------------------------------------------------------------------------------
!
SUBROUTINE TRANSF_1D_2D(PTAB1D,PTAB2D)
REAL, DIMENSION (:)      :: PTAB1D
REAL, DIMENSION (:,:)    :: PTAB2D
!
INTEGER                  :: JIJ, JI, JJ
!
JIJ = 0
DO JJ=1,SIZE(PTAB2D,2)
  DO JI=1,SIZE(PTAB2D,1)
    JIJ = JIJ + 1
    PTAB2D(JI,JJ) = PTAB1D(JIJ)
  ENDDO
ENDDO
!
END SUBROUTINE TRANSF_1D_2D
!
!-------------------------------------------------------------------------------
!
!       ###########################################################
SUBROUTINE  CALC(KIU,KKU,PP,PT,PR,   &
                            PCAPEP,PCINP,PDCAPE,PCAPEPMAX,PCINPMAX)
!       ###########################################################
!
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CST
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
INTEGER, INTENT(IN) :: KIU,KKU
REAL, DIMENSION(KIU,KKU), INTENT(IN)    :: PP    ! Pressure in hPa
REAL, DIMENSION(KIU,KKU), INTENT(INOUT) :: PT    ! Temperature in Celcius
REAL, DIMENSION(KIU,KKU), INTENT(INOUT) :: PR    ! vapor mixing ratio in g/kg
REAL, DIMENSION(KIU,KKU), INTENT(OUT)   :: PCAPEP
REAL, DIMENSION(KIU,KKU), INTENT(OUT)   :: PCINP
REAL, DIMENSION(KIU,KKU), INTENT(OUT)   :: PDCAPE
REAL, DIMENSION(KIU),     INTENT(OUT)   :: PCAPEPMAX
REAL, DIMENSION(KIU),     INTENT(OUT)   :: PCINPMAX
!
INTEGER, PARAMETER :: NBITER=4
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(KIU,KKU) ::  ZEV, ZES ! vapor pressure and saturation vapor pressure
REAL, DIMENSION(KIU,KKU,KKU) :: ZTVPDIF, ZTLP, ZTLVP
REAL, DIMENSION(KIU,KKU) :: ZTVD
REAL, DIMENSION(KIU,KKU) :: ZPAP
!
REAL :: ZCPVMCL ! CPV - CL
REAL :: ZEPS ! Rd/Rv
REAL, DIMENSION(KIU) :: ZRS ! saturation vapor mixing ratio
REAL, DIMENSION(KIU) :: ZALV ! Latent heat
REAL, DIMENSION(KIU) :: ZSP ! total entropy conserved under psueudo-adiabatic transformations
REAL, DIMENSION(KIU) :: ZAH ! enthalpie
REAL, DIMENSION(KIU) :: ZEM,ZSLOPE,ZTG,ZRG,ZALV1,ZAHG,ZTC,ZENEW
REAL, DIMENSION(KIU) :: ZEG,ZSPD,ZRGD0,ZTGD0,ZPM,ZTVM,ZSPG,ZTVDIFM
REAL, DIMENSION(KIU) :: ZRH ! relative humidity
REAL, DIMENSION(KIU) :: ZPLCL ! pressure of condensation level
REAL, DIMENSION(KIU) :: ZCHI,ZSUM,ZRG0,ZTG0,ZSLP,ZCPW,ZSUM2
INTEGER :: JKLOOP,J,K,JH
INTEGER, DIMENSION(KIU) ::  ICB,INBP, IMAX
!
!-------------------------------------------------------------------------------
!
!*       1.     INITIALIZATIONS
!               ---------------
!
!*       1.1   Assign values of thermodynamic constants
!
ZCPVMCL=2320.0
ZEPS=XRD/XRV
!
PR =PR*0.001
ZEV=PR*PP/(ZEPS+PR)
ZES=6.112*EXP(17.67*PT/(243.5+PT))  ! Eq (4.4.14)
PT =PT+XTT ! PT is now in Kelvin
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTATION OF PROPERTIES OF SAMPLES OF AIR
!               -------------------------------------------
!
! Begin outer loop, which cycles through parcel origin levels I
!
DO  JKLOOP=1,KKU     !  loop 2 on vertical levels
! Calculation limited to air parcels origin below 100 hPa
  IF (MINVAL(PP(:,JKLOOP))<100.) EXIT
!
!*       2.1   Various conserved parcel quantities
!
  ZRS(:)=ZEPS*ZES(:,JKLOOP)/(PP(:,JKLOOP)-ZES(:,JKLOOP))
  ZALV(:)=XLVTT-ZCPVMCL*(PT(:,JKLOOP)-XTT)                   ! Eq (4.4.4)
  WHERE(ZEV(:,JKLOOP)<1.0E-6)
    ZEM(:)=1.E-6
  ELSEWHERE
    ZEM(:)=ZEV(:,JKLOOP)
  ENDWHERE
  ! pseudo-adiabatic entropy
  ZSP(:)=XCPD*LOG(PT(:,JKLOOP))                         &
         - XRD*LOG(PP(:,JKLOOP)-ZEV(:,JKLOOP))  &
         + ZALV*PR(:,JKLOOP)/PT(:,JKLOOP)-PR(:,JKLOOP)*XRV*LOG(ZEM(:)/ZES(:,JKLOOP))
  ! enthalpy
  ZAH(:)=(XCPD+PR(:,JKLOOP)*XCL)*PT(:,JKLOOP)+ZALV*PR(:,JKLOOP)   ! Eq (4.5.23)
!
!*       2.2  Temperature and mixing ratio of the parcel at
!             level JKLOOP saturated by a wet bulb process
!
  ZSLOPE(:)=XCPD+ZALV(:)*ZALV(:)*ZRS(:)/(XRV*PT(:,JKLOOP)*PT(:,JKLOOP))
  ZTG(:)=PT(:,JKLOOP)
  ZRG(:)=ZRS(:)
  DO  J=1,NBITER
    ZALV1(:)=XLVTT-ZCPVMCL*(ZTG(:)-XTT)
    ZAHG(:)=(XCPD+XCL*ZRG(:))*ZTG(:)+ZALV1(:)*ZRG(:)
    ZTG(:)=ZTG(:)+(ZAH(:)-ZAHG(:))/ZSLOPE(:)
    ZTC(:)=ZTG(:)-XTT
    ZENEW(:)=6.112*EXP(17.67*ZTC(:)/(243.5+ZTC(:)))
    ZRG(:)=ZEPS*ZENEW(:)/(PP(:,JKLOOP)-ZENEW(:))
  ENDDO
!
!*       2.3  Calculate conserved variable at top of downdraft
!
  ZEG(:)=ZRG(:)*PP(:,JKLOOP)/(ZEPS+ZRG(:))
  ZSPD(:)=XCPD*LOG(ZTG(:))-XRD*LOG(PP(:,JKLOOP)-ZEG(:))+ ZALV1(:)*ZRG(:)/ZTG(:)
  ZTVD(:,JKLOOP)=ZTG(:)*(1.+ZRG(:)/ZEPS)/(1.+ZRG(:))                      &
                -PT(:,JKLOOP)*(1.+PR(:,JKLOOP)/ZEPS)/ (1.+PR(:,JKLOOP))
  WHERE(PP(:,JKLOOP).LT.100.0)
    ZTVD(:,JKLOOP)=0.0
  ENDWHERE
  ZRGD0(:)=ZRG(:)
  ZTGD0(:)=ZTG(:)
!
!*       2.4   Find lifted condensation pressure
!
  ZRH(:)=PR(:,JKLOOP)/ZRS(:)
  WHERE(ZRH(:)>1.)
    ZRH(:)=1.0
  ENDWHERE
  ZCHI(:)=PT(:,JKLOOP)/(1669.0-122.0*ZRH(:)-PT(:,JKLOOP))
  ZPLCL(:)=1.0
  WHERE(ZRH(:).GT.0.0)
    ZPLCL=PP(:,JKLOOP)*(ZRH**ZCHI)
  ENDWHERE
!
!*       2.5     Begin updraft loop
!
  ZSUM(:)=0.0
  ZRG0(:)=PR(:,JKLOOP)
  ZTG0(:)=PT(:,JKLOOP)
  DO  J=JKLOOP,KKU     ! inner loop  of the ascent for parcel JKLOOP
    !
    ! estimates of the rates of change of the entropies
    !  with temperature at constant pressure
    !
    ZRS(:)=ZEPS*ZES(:,J)/(PP(:,J)-ZES(:,J))
    ZALV=XLVTT-ZCPVMCL*(PT(:,J)-XTT)
    ZSLP=(XCPD+ZRS*XCL+ZALV*ZALV*ZRS/(XRV*PT(:,J)*PT(:,J)))/PT(:,J)
    !
    ! lifted parcel temperature below its LCL
    DO JH=1,KIU
    !
      IF(PP(JH,J).GE.ZPLCL(JH)) THEN
        ZTLP(JH,JKLOOP,J)=PT(JH,JKLOOP)*(PP(JH,J)/PP(JH,JKLOOP))**(XRD/XCPD) !  dry adiabat
        ZTLVP(JH,JKLOOP,J)=ZTLP(JH,JKLOOP,J)*(1.+PR(JH,JKLOOP)/ZEPS)/(1.+PR(JH,JKLOOP))
                                                                 ! vapor mixing of parcel JKLOOP
        ZTVPDIF(JH,JKLOOP,J)=ZTLVP(JH,JKLOOP,J)-PT(JH,J)*(1.+PR(JH,J)/ZEPS)/(1.+PR(JH,J)) ! Tvp -Tva
      ELSE
      !
      ! iteratively calculate lifted parcel temperature and mixing ratios
      ! for both reversible and pseudo-adiabatic ascent:
      !  do pseudo-adiabatic ascent
        ZTG(JH)=PT(JH,J)
        ZRG(JH)=ZRS(JH)
        DO  K=1,NBITER
          ZCPW(JH)=0.0
          IF(J.GT.1)THEN
            ZCPW(JH)=ZSUM(JH)+XCL*0.5*(ZRG0(JH)+ZRG(JH))*(LOG(ZTG(JH))-LOG(ZTG0(JH)))
          END IF
          ZEM(JH)=ZRG(JH)*PP(JH,J)/(ZEPS+ZRG(JH))
          ZALV(JH)=XLVTT-ZCPVMCL*(ZTG(JH)-XTT)
          ZSPG(JH)=XCPD*LOG(ZTG(JH))-XRD*LOG(PP(JH,J)-ZEM(JH))+ZCPW(JH)+ZALV(JH)*ZRG(JH)/ZTG(JH)
          ZTG(JH)=ZTG(JH)+(ZSP(JH)-ZSPG(JH))/ZSLP(JH)
          ZTC(JH)=ZTG(JH)-XTT
          ZENEW(JH)=6.112*EXP(17.67*ZTC(JH)/(243.5+ZTC(JH)))
          ZRG(JH)=ZEPS*ZENEW(JH)/(PP(JH,J)-ZENEW(JH))
        END DO
        ZTLVP(JH,JKLOOP,J)=ZTG(JH)*(1.+ZRG(JH)/ZEPS)/(1.+ZRG(JH))
        ZTVPDIF(JH,JKLOOP,J)=ZTLVP(JH,JKLOOP,J)-PT(JH,J)*(1.+PR(JH,J)/ZEPS)/(1.+PR(JH,J))
        ZRG0(JH)=ZRG(JH)
        ZTG0(JH)=ZTG(JH)
        ZSUM(JH)=ZCPW(JH)
      END IF
    END DO  ! end of the loop for the horiz. index
  END DO   ! end of inner loop of the ascent for parcel JKLOOP
  IF(JKLOOP.EQ.1) CYCLE
!
!*       2.5     Begin downdraft loop
!
  ZSUM2=0.0
  DO  J=JKLOOP-1,1,-1 ! loop 3 from top to bottom
  !
  ! estimate of the rate of change of entropy
  !  with temperature at constant pressure
  !
    ZRS=ZEPS*ZES(:,J)/(PP(:,J)-ZES(:,J))
    ZALV=XLVTT-ZCPVMCL*(PT(:,J)-XTT)
    ZSLP=(XCPD+ZRS*XCL+ZALV*ZALV*ZRS/(XRV*PT(:,J)*PT(:,J)))/PT(:,J)
    ZTG=PT(:,J)
    ZRG=ZRS
    !
    ! downdraft temperature
    !
    DO  K=1,NBITER
     DO JH=1,KIU
      ZEM(JH)=ZRG(JH)*PP(JH,J)/(ZEPS+ZRG(JH))
      IF (PP(JH,J) >= ZEM(JH)) THEN
       ZCPW(JH)=ZSUM2(JH)+XCL*0.5*(ZRGD0(JH)+ZRG(JH))*(LOG(ZTG(JH))-LOG(ZTGD0(JH)))
       ZALV(JH)=XLVTT-ZCPVMCL*(ZTG(JH)-XTT)
       ZSPG(JH)=XCPD*LOG(ZTG(JH))-XRD*LOG(PP(JH,J)-ZEM(JH))+ZCPW(JH)+ZALV(JH)*ZRG(JH)/ZTG(JH)
       ZTG(JH)=ZTG(JH)+(ZSPD(JH)-ZSPG(JH))/ZSLP(JH)
       ZTC(JH)=ZTG(JH)-XTT
       ZENEW(JH)=6.112*EXP(17.67*ZTC(JH)/(243.5+ZTC(JH)))
       ZRG(JH)=ZEPS*ZENEW(JH)/(PP(JH,J)-ZENEW(JH))
      END IF
     END DO
    END DO
    ZSUM2=ZCPW
    ZTGD0=ZTG
    ZRGD0=ZRG
    ZTLVP(:,JKLOOP,J)=ZTG*(1.+ZRG/ZEPS)/(1.+ZRG)
    ZTVPDIF(:,JKLOOP,J)=ZTLVP(:,JKLOOP,J)-PT(:,J)*(1.+PR(:,J)/ZEPS)/(1.+PR(:,J))
    WHERE(PP(:,JKLOOP).LT.100.0)
      ZTVPDIF(:,JKLOOP,J)=0.0
    ENDWHERE
    WHERE(ZTVPDIF(:,JKLOOP,J)>0)
      ZTVPDIF(:,JKLOOP,J)=0.0
    ENDWHERE
  END DO ! loop 3 from top to bottom
END DO  ! end loop 2 on vertical levels
!
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTATION OF DCAPE, PA, NA (from pseudo-adiabatic ascent), CAPE
!              ----------------------------------------------------------------
!
!
  PCAPEP(:,:)=0.0
  PDCAPE(:,:)=0.0
  ZPAP(:,:)=0.0
  PCINP(:,:)=0.0
DO  JKLOOP=1,KKU ! loop 4
! Calculation limited to air parcels origin below 100 hPa
  IF (MINVAL(PP(:,JKLOOP))<100.) EXIT
!
!*       3.1 lifted condensation pressure
!
  ZRS=ZEPS*ZES(:,JKLOOP)/(PP(:,JKLOOP)-ZES(:,JKLOOP))  !saturation vapor mixing ratio
  ZRH=PR(:,JKLOOP)/ZRS  ! relative humidity
  WHERE(ZRH>1)
    ZRH=1.
  ENDWHERE
  ZCHI=PT(:,JKLOOP)/(1669.0-122.0*ZRH-PT(:,JKLOOP))
  ZPLCL=1.0
  WHERE(ZRH.GT.0.0)
    ZPLCL=PP(:,JKLOOP)*(ZRH**ZCHI)
  ENDWHERE
!
!*       3.2 lifted condensation level and maximum level of positive buoyancy
!
  ICB=KKU ! condensation level
  INBP=1  ! level of neutral buoyancy for pseudo-adiabatic ascent
  DO  J=KKU,JKLOOP,-1
    DO JH=1,KIU
      IF(PP(JH,J).LT.ZPLCL(JH)) ICB(JH)=MIN(ICB(JH),J)
      IF(ZTVPDIF(JH,JKLOOP,J).GT.0.0) INBP(JH)=MAX(INBP(JH),J)
    END DO
  END DO
  DO JH=1,KIU
    IMAX(JH)=MAX(INBP(JH),JKLOOP)
  END DO
  DO JH=1,KIU
    ZTVPDIF(JH,JKLOOP,IMAX(JH))=0.0
  END DO
!
!*       3.3 updraft loops
!
  DO JH=1,KIU
    IF(INBP(JH).GT.JKLOOP)THEN
      DO J=JKLOOP+1,INBP(JH)
        ZTVM(JH)=0.5*(ZTVPDIF(JH,JKLOOP,J)+ZTVPDIF(JH,JKLOOP,J-1))
        ZPM(JH)=0.5*(PP(JH,J)+PP(JH,J-1))
        IF(ZTVM(JH).LE.0.0)THEN
          PCINP(JH,JKLOOP)=PCINP(JH,JKLOOP)-XRD*ZTVM(JH)*(PP(JH,J-1)-PP(JH,J))/ZPM(JH)
        ELSE
          ZPAP(JH,JKLOOP)=ZPAP(JH,JKLOOP)+XRD*ZTVM(JH)*(PP(JH,J-1)-PP(JH,J))/ZPM(JH)
        END IF
      END DO
      PCAPEP(JH,JKLOOP)=ZPAP(JH,JKLOOP)-PCINP(JH,JKLOOP)
    END IF
  ENDDO  ! loop on the horiz. index
!
!*       3.4 find DCAPE
!
  IF(JKLOOP.EQ.1) CYCLE
  DO  J=JKLOOP-1,1,-1
    ZTVDIFM=ZTVPDIF(:,JKLOOP,J+1)
    IF(JKLOOP.EQ.(J+1)) ZTVDIFM=ZTVD(:,JKLOOP)
    ZTVM=0.5*(ZTVPDIF(:,JKLOOP,J)+ZTVDIFM)
    ZPM=0.5*(PP(:,J)+PP(:,J+1))
    WHERE(ZTVM.LT.0.0)
      PDCAPE(:,JKLOOP)=PDCAPE(:,JKLOOP)-XRD*ZTVM*(PP(:,J)-PP(:,J+1))/ZPM
    ENDWHERE
  END DO
!
END DO   ! end loop 4
!
!*       3.5 find CAPEMAX, CINMAX
!
PCAPEPMAX = 0.0
PCINPMAX  = 0.0
!
DO JKLOOP=1,KKU
  WHERE (PCAPEP(:,JKLOOP) >  PCAPEPMAX(:) )
    PCAPEPMAX= PCAPEP(:,JKLOOP)
    PCINPMAX = PCINP(:,JKLOOP)
  ENDWHERE
END DO
!
END SUBROUTINE CALC
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE CALCSOUND
