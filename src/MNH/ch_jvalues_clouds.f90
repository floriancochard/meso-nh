!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/06/23 11:35:48
!-----------------------------------------------------------------
!      #############################
       MODULE MODI_CH_JVALUES_CLOUDS
!      #############################
        INTERFACE
         SUBROUTINE CH_JVALUES_CLOUDS(PZENITH, KPJVMAX, KJVLEVEL, PAZ, PLWC, PFCLD)
          IMPLICIT NONE
          INTEGER, INTENT(IN)                   :: KPJVMAX, KJVLEVEL
          REAL, DIMENSION(:,:),     INTENT(IN)  :: PZENITH
          REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PLWC, PAZ
          REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PFCLD
         END SUBROUTINE CH_JVALUES_CLOUDS
        END INTERFACE
       END MODULE MODI_CH_JVALUES_CLOUDS
!-----------------------------------------------------------------------------
!      #####################################################################
         SUBROUTINE CH_JVALUES_CLOUDS(PZENITH, KPJVMAX, KJVLEVEL, PAZ, PLWC, PFCLD)
!      #####################################################################
! 
!!
!!*** *CH_JVALUES_CLOUDS* 
!
!!    PURPOSE Application of a cloud correction factor to the clear-sky values of 
!!            photolysis rate constants 
!!    -------
!!
!!**  METHOD
!!    ------
!!    The method used to correct for cloud cover is taken from RADM 
!!    (Chang et al., 1987; Madronich, 1987)
!!    The correction of clear-sky values depends on whether the location 
!!    is below, above or within the cloud. The equation allows for enhancement 
!!    of photolysis rates above the cloud due to the reflected 
!!    radiation from the cloud. Below cloud photolysis rates will be lower 
!!    than the clear-sky values due to the reduced transmission of radiation 
!!    through the cloud.
!!    Within the cloud, the cloud correction factor is a simple linear 
!!    interpolation of the below cloud factor at cloud base to the above 
!!    cloud factor at cloud top. 
!!       
!!    REFERENCE
!!    ---------
!!   
!!
!!    AUTHOR
!!    ------
!!    Celine Mari (LA)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original 21/01/00
!!    05/03/05  P. Tulet (CNRM/GMEI) Update for Arome
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CST, ONLY        : XRHOLW
USE MODD_PARAMETERS
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
!!
IMPLICIT NONE
INTEGER, INTENT(IN)                   :: KPJVMAX, KJVLEVEL
REAL, DIMENSION(:,:),     INTENT(IN)  :: PZENITH
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PLWC, PAZ
REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: PFCLD
!!
!!    LOCAL VARIABLES
!!    ---------------
!!
INTEGER                                             :: JK, ILOOP, JI,JJ
INTEGER                                             :: IIU,IJU,IKU,IKB,IKE 
INTEGER,DIMENSION(SIZE(PLWC,1),SIZE(PLWC,2))          :: ICL_TOP_ST, ICL_BASE_ST
INTEGER                                             :: ICL_TOP, ICL_BASE
REAL,DIMENSION(SIZE(PLWC,1),SIZE(PLWC,2))             :: ZTAU, ZTR, ZFCLDB
REAL,DIMENSION(SIZE(PLWC,1),SIZE(PLWC,2),KPJVMAX)     :: ZFCLDA, X1, X2
REAL,DIMENSION(SIZE(PLWC,1),SIZE(PLWC,2))             :: ZLWCBAR, ZLWP
REAL                                                :: ZLWCMIN     
REAL,DIMENSION(SIZE(PLWC,1),SIZE(PLWC,2))             :: ZCOSZEN     ! cosine of zenithal angle
LOGICAL,DIMENSION(SIZE(PLWC,1),SIZE(PLWC,2))          :: GMASK, GCLOUD
REAL, DIMENSION(KPJVMAX)                            :: ALPHA       ! reaction dependant coefficient
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
IIU = SIZE(PLWC,1)
IJU = SIZE(PLWC,2)
IKU = SIZE(PLWC,3)
IKB = 1+JPVEXT
!
IKE = IKU-JPVEXT

!
!
!           1. Reaction dependent coefficient
!           ---------------------------------
!
!
ALPHA(1:KPJVMAX)=1.
!! NO2
ALPHA(2)=1.2
!! O3
ALPHA(3)=0.7
!! NO3
ALPHA(5)=1.3
ALPHA(6)=1.3
!! HCHO
ALPHA(14)=1.1
!! ALD
ALPHA(17)=0.9
!
!            2. Cloud Top and Cloud Base 
!            ---------------------------
!
  ICL_TOP_ST(:,:) = 1 
  ZLWCMIN  = 1.E-4  
!
  GCLOUD(:,:)=.TRUE.
  GMASK(:,:)=.TRUE.
  !
  DO JK=KJVLEVEL-1,1,-1
    !
    ! *  Cloud top 
    !
    WHERE ( (GMASK(:,:)) .AND. (PLWC(:,:,JK) > ZLWCMIN) ) 
      GMASK(:,:)=.FALSE.
      ICL_TOP_ST(:,:)=JK
    ENDWHERE
    !
  END DO
  !
  WHERE (ICL_TOP_ST(:,:).EQ.1) 
      ICL_TOP_ST(:,:)=0.
      GCLOUD(:,:)=.FALSE.
  ENDWHERE
    !
    ! * Cloud base
    ! 
   ICL_BASE_ST(:,:)=IKE-IKB+1
   GMASK(:,:)=.TRUE.
    !
  DO JK=1,IKE-IKB+1,1
     WHERE ( (GMASK(:,:)) .AND. (PLWC(:,:,JK) > ZLWCMIN) ) 
            GMASK(:,:)=.FALSE.
            ICL_BASE_ST(:,:)=JK
     END WHERE
  END DO
  WHERE (ICL_BASE_ST(:,:).EQ.IKE-IKB+1)
      ICL_BASE_ST(:,:)=0.
      GCLOUD(:,:)=.FALSE.
  END WHERE
  WHERE (ICL_BASE_ST(:,:).EQ.ICL_TOP_ST(:,:))
      GCLOUD(:,:)=.FALSE.
  END WHERE
    !
!
!
!            3. Cosinus of the zenith angle 
!            ------------------------------
!
!
ZCOSZEN(:,:) = COS(PZENITH(:,:))
ZCOSZEN(:,:)=MAX(ZCOSZEN(:,:),0.5)
!
PFCLD(:,:,:,:)=1.0
ZTR(:,:)=0.
ZTAU(:,:)=0.
ZLWCBAR(:,:)=0.
!
DO JI = 1, IIU                      
 DO JJ = 1, IJU                
  IF (GCLOUD(JI,JJ))  THEN
!
!            4.  Mean LWC of the cloud
!            -------------------------
!
     ICL_BASE=ICL_BASE_ST(JI,JJ)
     ICL_TOP=ICL_TOP_ST(JI,JJ)
     !
     DO JK=ICL_BASE,ICL_TOP,1
      ZLWCBAR(JI,JJ) = ZLWCBAR(JI,JJ) + PLWC(JI,JJ,JK) 
     END DO
     ZLWCBAR(JI,JJ) =  ZLWCBAR(JI,JJ) /(ICL_TOP-ICL_BASE)
!
!            Liquid water path in g/m2
!
     ZLWP(JI,JJ) = ZLWCBAR(JI,JJ) * 1.E3 * ( PAZ(JI,JJ,ICL_TOP) - PAZ(JI,JJ,ICL_BASE) )
!
!            5. Cloud optical depth (Stephens, JAS, 1978) 
!            -------------------------------------------------------------
!
   IF (ZLWP(JI,JJ).GE.10.) THEN
      ZTAU(JI,JJ) = 10.**(0.2633 + 1.7095 * LOG(LOG10(ZLWP(JI,JJ)))) 
   ELSE
      ZTAU(JI,JJ) = 0.
   END IF
!
!            6. Energy transmission coefficient 
!            ----------------------------------
!
   ZTR(JI,JJ) = (5.0 - EXP(-ZTAU(JI,JJ))) / (4.0 + 0.42 * ZTAU(JI,JJ))
!
  END IF
 END DO
END DO
!
DO JI = 1, IIU 
 DO JJ = 1, IJU 
  !
  IF (GCLOUD(JI,JJ).AND.(ZTAU(JI,JJ).GT.5.0)) THEN 
     !
     ICL_BASE=ICL_BASE_ST(JI,JJ)
     ICL_TOP=ICL_TOP_ST(JI,JJ)
!
!            7. Calculate cloud correction factor
!            ------------------------------------
!
    !
    ! Below cloud 
    !
      ZFCLDB(JI,JJ) = 1.0 + (1.6 * ZTR(JI,JJ) *  ZCOSZEN(JI,JJ) - 1.0)
    !
    ! Above cloud 
    !
     DO ILOOP = 1, KPJVMAX
      ZFCLDA(JI,JJ,ILOOP) = 1.0 + (1.0 + ALPHA(ILOOP)*(1.0 - ZTR(JI,JJ)) * &
                            ZCOSZEN(JI,JJ) - 1.0)
    !
      X1(JI,JJ,ILOOP) = (ZFCLDA(JI,JJ,ILOOP) - ZFCLDB(JI,JJ))/ &
                  (PAZ(JI,JJ,ICL_TOP)-PAZ(JI,JJ,ICL_BASE))
      X2(JI,JJ,ILOOP) = (PAZ(JI,JJ,ICL_TOP)*ZFCLDB(JI,JJ)-PAZ(JI,JJ,ICL_BASE)* &
                   ZFCLDA(JI,JJ,ILOOP))/ &
                  (PAZ(JI,JJ,ICL_TOP)-PAZ(JI,JJ,ICL_BASE))
    !
    !In cloud layer: linear interpolation
    !
    DO JK = 1, IKE-IKB+1
    !
      PFCLD(JI,JJ,JK,ILOOP) = X1(JI,JJ,ILOOP) *  PAZ(JI,JJ,JK) + X2(JI,JJ,ILOOP)
      IF ( JK .LE. ICL_BASE) PFCLD(JI,JJ,JK,ILOOP) = ZFCLDB(JI,JJ)
      IF ( JK .GE. ICL_TOP)  PFCLD(JI,JJ,JK,ILOOP) = ZFCLDA(JI,JJ,ILOOP)
    !
    ENDDO
      PFCLD(JI,JJ,KJVLEVEL,ILOOP) = PFCLD(JI,JJ,IKE-IKB+1,ILOOP)
    ENDDO ! ILOOP
   END IF
    !
 END DO
END DO
!
!
END SUBROUTINE CH_JVALUES_CLOUDS
!

