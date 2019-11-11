!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/10/18 12:20:58
!-----------------------------------------------------------------
!!   ##############################
     MODULE MODI_CH_AER_DEPOS
!!   ##############################
!!
INTERFACE
!
SUBROUTINE CH_AER_DEPOS(NSPLITR, PTSTEP, PZZ, PRHODREF, PRT, PRS,&
                        PRHODJ, PSVT, PMI, PTHT, PPABST, PEVAP3D, PSEDA)

IMPLICIT NONE

INTEGER,  INTENT(IN) :: NSPLITR
REAL,     INTENT(IN) :: PTSTEP
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF, PZZ, PRHODJ 
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PRT, PRS, PMI
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST, PEVAP3D
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSEDA, PSVT


END SUBROUTINE CH_AER_DEPOS
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_DEPOS
!!
!!   #######################################
SUBROUTINE CH_AER_DEPOS(NSPLITR, PTSTEP, PZZ, PRHODREF, PRT, PRS,&
                        PRHODJ, PSVT, PMI, PTHT, PPABST, PEVAP3D, PSEDA)
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!    Pierre TULET (GMEI) 
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
! Entry variables:
!
!!
!!   IMPLICIT ARGUMENTS
!
USE MODD_CH_AEROSOL
USE MODE_AERO_PSD
USE MODD_CST,       ONLY: XMD,       &! Molar mass of dry air
                          XPI
USE MODD_NSV,  ONLY : NSV_C2R2BEG, NSV_AERBEG,NSV_AEREND,NSV_AER, &
                      NSV_AERDEP, NSV_AERDEPBEG
USE MODD_PARAM_n,    ONLY: CCLOUD
USE MODI_AER_WET_DEP_KMT_WARM

!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
INTEGER,  INTENT(IN) :: NSPLITR
REAL,     INTENT(IN) :: PTSTEP
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PTHT,PRHODREF, PZZ, PRHODJ 
REAL,  DIMENSION(:,:,:,:),  INTENT(IN)    :: PRT, PRS, PMI
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PPABST, PEVAP3D
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSEDA, PSVT
!
!*       0.2   Declarations of local variables
!
INTEGER :: JSV, JN, JJ
REAL    :: ZRGMIN, ZINIRADIUSI, ZINIRADIUSJ
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE)   :: ZRG, ZSIG, ZN, ZRHOP
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE*3) :: ZSVAER
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NSV_AER)  :: ZSVT
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE*3) :: ZPM, ZFLUXSED, ZFLUXMAX, ZPMOLD
REAL,  DIMENSION(JPMODE*3)     :: ZPMIN
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE)   :: ZMASSMIN, ZMI
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3)) :: ZSUM, ZRRS, ZRCS
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NSP+NCARB+NSOA,JPMODE) :: ZCTOTA, ZCCTOT
REAL,  DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2)) :: ZSEA, ZTOWN
!-------------------------------------------------------------------------------
!
!*       1.1   compute dimensions of arrays
!
!
ZFLUXSED(:,:,:,:) = 0.
ZRCS(:,:,:) = PRS(:,:,:,2) * PTSTEP / PRHODJ(:,:,:) 
ZRRS(:,:,:) = PRS(:,:,:,3) * PTSTEP / PRHODJ(:,:,:)
DO JN=NSV_AERBEG,NSV_AEREND
ZSVT(:,:,:,JN-NSV_AERBEG+1) = PSVT(:,:,:,JN)
ENDDO
!
IF (CRGUNIT=="MASS") THEN
    ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
    ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
    ZINIRADIUSI = XINIRADIUSI 
    ZINIRADIUSJ = XINIRADIUSJ
END IF
!
!Get minimum values possible
ZPMIN(1) = XN0IMIN
ZRGMIN = ZINIRADIUSI
ZPMIN(2) = ZPMIN(1) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGIMIN)**2) 
ZPMIN(3) = ZPMIN(1) * (ZRGMIN**6)*EXP(18. * LOG(XSIGIMIN)**2)
ZPMIN(4) = XN0JMIN
ZRGMIN = ZINIRADIUSJ
ZPMIN(5) = ZPMIN(4) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGJMIN)**2) 
ZPMIN(6) = ZPMIN(4) * (ZRGMIN**6)*EXP(18. * LOG(XSIGJMIN)**2)
!
!
CALL  PPP2AERO(ZSVT, PRHODREF, &
               PSIG3D=ZSIG, PRG3D=ZRG, PCTOTA=ZCTOTA, PM3D=ZPM)
!
ZPMOLD(:,:,:,:) = ZPM(:,:,:,:)
!
ZRHOP(:,:,:,:)   = 0.
ZCCTOT(:,:,:,:,:)= 0.
ZMI(:,:,:,:) = 0.
ZSEA(:,:) = 0. 
ZTOWN(:,:) = 0. 

DO JN=1,JPMODE
  ZSUM(:,:,:)=0.
  DO JJ=1,NSP+NCARB+NSOA
    ZSUM(:,:,:)=ZSUM(:,:,:)+ZCTOTA(:,:,:,JJ,JN)/XRHOI(JJ)
  END DO
  DO JJ=1,NSP+NCARB+NSOA
    ZCCTOT(:,:,:,JJ,JN) = ZCTOTA(:,:,:,JJ,JN)/XRHOI(JJ)/ZSUM(:,:,:)
    ZRHOP(:,:,:,JN)=ZRHOP(:,:,:,JN)+ZCCTOT(:,:,:,JJ,JN)*XRHOI(JJ)
    ZMI(:,:,:,JN)= ZMI(:,:,:,JN) + PMI(:,:,:,JJ)
  ENDDO
  ZMI(:,:,:,JN)= ZMI(:,:,:,JN) / (NSP+NCARB+NSOA)
  ZMASSMIN(:,:,:,JN)= ZPMIN(JN*3-1) * XMD * XPI * 4./3. * ZRHOP(:,:,:,JN)  / &
                     (ZMI(:,:,:,JN)*1.d15*PRHODREF(:,:,:))
ENDDO
!
!     Compute acquous aerosol mass vector from moment scalar vector
DO JN=1,JPMODE*2
   ZSVAER(:,:,:,JPMODE+JN) = PSVT(:,:,:,NSV_AERDEPBEG-1+JN)
ENDDO

DO JN=1,JPMODE
   ZSVAER(:,:,:,JN) =   ZPM(:,:,:,JN*3-1) * XMD * XPI * &
                        4./3. * ZRHOP(:,:,:,JN)  /      &
                       (ZMI(:,:,:,JN)*1.d15*PRHODREF(:,:,:))
ENDDO

!     3.5 Mass transfer between dry mass and in-cloud mass aerosols
SELECT CASE (CCLOUD)
  
  CASE ('KESS','REVE','ICE3','ICE4')
! One moment cloud scheme
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP, PZZ, PRHODREF,    &
                              PRT(:,:,:,2), PRT(:,:,:,3), ZRCS,  &
                              ZRRS, ZSVAER, PTHT, PPABST,  ZRG,  &
                              PEVAP3D, JPMODE,ZRHOP, ZMASSMIN,   &
                              PSEA=ZSEA, PTOWN=ZTOWN)
  CASE ('KHKO','C2R2','C3R5')
! Two moment cloud scheme
  CALL AER_WET_DEP_KMT_WARM  (NSPLITR, PTSTEP, PZZ, PRHODREF,    &
                              PRT(:,:,:,2), PRT(:,:,:,3), ZRCS,  &
                              ZRRS, ZSVAER, PTHT, PPABST,  ZRG,  &
                              PEVAP3D, JPMODE,ZRHOP, ZMASSMIN,   &
                              PSEA=ZSEA, PTOWN=ZTOWN,            &
                              PCCT=PSVT(:,:,:,NSV_C2R2BEG+1),    &
                              PCRT=PSVT(:,:,:,NSV_C2R2BEG+2) )

END SELECT

! 
DO JN=1,JPMODE*2
   PSVT(:,:,:,NSV_AERDEPBEG-1+JN) = ZSVAER(:,:,:,JPMODE+JN)
ENDDO

DO JN=1,JPMODE
   ! Return from dry mass on m3
   ZPM(:,:,:,JN*3-1) = ZSVAER(:,:,:,JN)  *                    &
                       ZMI(:,:,:,JN)*1.d15*PRHODREF(:,:,:) /  &
                      (XMD * XPI * 4./3. * ZRHOP(:,:,:,JN) )

   ! Compute  m0 and m6 using m3, Rg and sigma
   ZPM(:,:,:,JN*3-2) = ZPM(:,:,:,JN*3-1)/                     &
              ( (ZRG(:,:,:,JN)**3)*EXP(4.5 * LOG(ZSIG(:,:,:,JN))**2) ) 

   ZPM(:,:,:,JN*3)   =  ZPM(:,:,:,JN*3-2)*(ZRG(:,:,:,JN)**6) * &
                        EXP(18 *(LOG(ZSIG(:,:,:,JN)))**2)
ENDDO

! Filter minimum values
!DO JN=1,JPMODE
!  WHERE ((ZPM(:,:,:,NM0(JN)) .LE. ZPMIN(NM0(JN))).OR.&
!         (ZPM(:,:,:,NM3(JN)) .LE. ZPMIN(NM3(JN))).OR.& 
!         (ZPM(:,:,:,NM6(JN)) .LE. ZPMIN(NM6(JN)))) 
!    ZPM(:,:,:,NM0(JN)) = ZPMOLD(:,:,:,NM0(JN)) 
!    ZPM(:,:,:,NM3(JN)) = ZPMOLD(:,:,:,NM3(JN)) 
!    ZPM(:,:,:,NM6(JN)) = ZPMOLD(:,:,:,NM6(JN)) 
!  END WHERE
!ENDDO 

  ! Compute final sedimentation tendency with adding acqueous sedimentation
DO JN=1,JPIN
 PSEDA(:,:,:,JN) = PSEDA(:,:,:,JN) + (ZPM(:,:,:,JN) - ZPMOLD(:,:,:,JN)) / PTSTEP
END DO
!
END SUBROUTINE CH_AER_DEPOS
