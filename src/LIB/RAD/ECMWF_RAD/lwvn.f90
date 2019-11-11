!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE LWVN &
 &( KIDIA, KFDIA, KLON  , KLEV , KUAER &
 &, PABCU, PDBSL, PGA   , PGB &
 &, PADJD, PADJU, PCNTRB, PDBDT, PDWFSU &
 &)

!**** *LWVN*   - L.W., VERTICAL INTEGRATION, NEARBY LAYERS

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
!           TO GIVE LONGWAVE FLUXES OR RADIANCES

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1)  ; ABSORBER AMOUNTS
! PDBSL  : (KLON,KLEV*2)       ; SUB-LAYER PLANCK FUNCTION GRADIENT
! PGA, PGB                     ; PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PADJ.. : (KLON,KLEV+1)       ; CONTRIBUTION OF ADJACENT LAYERS
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PDBDT  : (KLON,NUA,KLEV)     ; LAYER PLANCK FUNCTION GRADIENT
! PDWFSU : (KLON,NSIL)         ; SPECTRAL DOWNWARD FLUX AT SURFACE

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
!     CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE

!     EXTERNALS.
!     ----------

!          *LWTT*

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!        JJ Morcrette 97-04-18 Revised Continuum + Surf.Emissiv.
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOELW    , ONLY : NSIL     ,NIPD     ,NTRA     ,NUA      ,&
            &NG1      ,NG1P1    ,WG1
!
USE MODI_LWTT
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KUAER



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------


REAL_B :: PABCU(KLON,NUA,3*KLEV+1)&
  &,  PDBSL(KLON,NSIL,KLEV*2)&
  &,  PGA(KLON,NIPD,2,KLEV)   , PGB(KLON,NIPD,2,KLEV)

REAL_B :: PADJD(KLON,KLEV+1)      , PADJU(KLON,KLEV+1)&
  &,  PCNTRB(KLON,KLEV+1,KLEV+1)&
  &,  PDBDT(KLON,NSIL,KLEV)   , PDWFSU(KLON,NSIL)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZTT(KLON,NTRA), ZTT1(KLON,NTRA), ZTT2(KLON,NTRA),  ZUU(KLON,NUA)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IBS, IDD, IM12, IMU, IND, INU, IXD, IXU,&
             &JA, JG, JK, JK1, JK2, JL, JNU

!     LOCAL REAL SCALARS
REAL_B :: ZWTR, ZWTR1, ZWTR2, ZWTR3, ZWTR4, ZWTR5, ZWTR6


!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PADJD(JL,JK) = _ZERO_
    PADJU(JL,JK) = _ZERO_
  ENDDO
ENDDO

!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------

DO JA = 1 , NTRA
  DO JL = KIDIA,KFDIA
    ZTT (JL,JA) = _ONE_
    ZTT1(JL,JA) = _ONE_
    ZTT2(JL,JA) = _ONE_
  ENDDO
ENDDO

DO JA = 1 , NUA
  DO JL = KIDIA,KFDIA
    ZUU(JL,JA) = _ZERO_
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------


!*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
!                  ---------------------------------

DO JK = 1 , KLEV

!*         2.1.1   DOWNWARD LAYERS
!                  ---------------

  IM12 = 2 * (JK - 1)
  IND = (JK - 1) * NG1P1 + 1
  IXD = IND
  INU = JK * NG1P1 + 1
  IXU = IND

  DO JG = 1 , NG1
    IBS = IM12 + JG
    IDD = IXD + JG

    DO JA = 1 , KUAER
      DO JL = KIDIA,KFDIA
        ZUU(JL,JA) = PABCU(JL,JA,IND) - PABCU(JL,JA,IDD)
      ENDDO
    ENDDO


    CALL LWTT &
     &( KIDIA        , KFDIA        , KLON &
     &, PGA(1,1,1,JK), PGB(1,1,1,JK)&
     &, ZUU          , ZTT &
     &)

    DO JL = KIDIA,KFDIA
      ZWTR1=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)
      ZWTR2=PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
      ZWTR3=PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
      ZWTR4=PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
      ZWTR5=PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)
      ZWTR6=PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZWTR=ZWTR1+ZWTR2+ZWTR3+ZWTR4+ZWTR5+ZWTR6
      PADJD(JL,JK) = PADJD(JL,JK) + ZWTR * WG1(JG)
      IF (JK == 1) THEN
        PDWFSU(JL,1)=PDWFSU(JL,1)+WG1(JG)*ZWTR1
        PDWFSU(JL,2)=PDWFSU(JL,2)+WG1(JG)*ZWTR2
        PDWFSU(JL,3)=PDWFSU(JL,3)+WG1(JG)*ZWTR3
        PDWFSU(JL,4)=PDWFSU(JL,4)+WG1(JG)*ZWTR4
        PDWFSU(JL,5)=PDWFSU(JL,5)+WG1(JG)*ZWTR5
        PDWFSU(JL,6)=PDWFSU(JL,6)+WG1(JG)*ZWTR6
      ENDIF
    ENDDO

!*         2.1.2   UPWARD LAYERS
!                  -------------

    IMU = IXU + JG
    DO JA = 1 , KUAER
      DO JL = KIDIA,KFDIA
        ZUU(JL,JA) = PABCU(JL,JA,IMU) - PABCU(JL,JA,INU)
      ENDDO
    ENDDO


    CALL LWTT &
     &( KIDIA        , KFDIA        , KLON &
     &, PGA(1,1,1,JK), PGB(1,1,1,JK)&
     &, ZUU          , ZTT &
     &)

    DO JL = KIDIA,KFDIA
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)&
       &+PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)&
       &+PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)&
       &+PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)&
       &+PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)&
       &+PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      PADJU(JL,JK+1) = PADJU(JL,JK+1) + ZWTR * WG1(JG)
    ENDDO

  ENDDO

  DO JL = KIDIA,KFDIA
    PCNTRB(JL,JK,JK+1) = PADJD(JL,JK)
    PCNTRB(JL,JK+1,JK) = PADJU(JL,JK+1)
    PCNTRB(JL,JK  ,JK) = _ZERO_
  ENDDO

ENDDO

DO JK = 1 , KLEV
  JK2 = 2 * JK
  JK1 = JK2 - 1

  DO JNU = 1 , NSIL
    DO JL = KIDIA,KFDIA
      PDBDT(JL,JNU,JK) = PDBSL(JL,JNU,JK1) + PDBSL(JL,JNU,JK2)
    ENDDO
  ENDDO
ENDDO

!-----------------------------------------------------------------------

RETURN
END SUBROUTINE LWVN
