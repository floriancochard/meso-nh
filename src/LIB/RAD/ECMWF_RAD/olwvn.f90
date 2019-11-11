!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE OLWVN ( KIDIA,KFDIA,KLON,KLEV,KUAER &
     &  , PABCU,PDBSL,PGA,PGB &
     &  , PADJD,PADJU,PCNTRB,PDBDT )
!
!**** *LWVN*   - L.W., VERTICAL INTEGRATION, NEARBY LAYERS
!
!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
!           TO GIVE LONGWAVE FLUXES OR RADIANCES
!
!**   INTERFACE.
!     ----------
!
!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU : (KLON,NUA,3*KLEV+1); ABSORBER AMOUNTS
! PDBSL  : (KLON,KLEV*2)     ; SUB-LAYER PLANCK FUNCTION GRADIENT
! PGA, PGB                     ; PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PADJ.. : (KLON,KLEV+1)     ; CONTRIBUTION OF ADJACENT LAYERS
! PCNTRB : (KLON,KLEV+1,KLEV+1); CLEAR-SKY ENERGY EXCHANGE MATRIX
! PDBDT  : (KLON,NUA,KLEV)   ; LAYER PLANCK FUNCTION GRADIENT
!
!        IMPLICIT ARGUMENTS :   NONE
!        --------------------
!
!     METHOD.
!     -------
!
!          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
!     CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE
!
!     EXTERNALS.
!     ----------
!
!          *LWTT*
!
!     REFERENCE.
!     ----------
!
!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
!
!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!-----------------------------------------------------------------------
!      IMPLICIT LOGICAL (L)
!
!#include "yomcst.h"
!#include "yoerad.h"
!#include "yoeaer.h"
!#include "yoelw.h"
!#include "yoesw.h"
!#include "yoerdu.h"


#include "tsmbkind.h"

USE YOEOLW   , ONLY : NISP     ,NIPD     ,NTRA     ,NUA      ,&
            &NG1      ,NG1P1    ,WG1


IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KUAER


!-----------------------------------------------------------------------
!
!*       0.1   ARGUMENTS
!              ---------
!
!
REAL_B :: PABCU(KLON,NUA,3*KLEV+1) &
     &  ,  PDBSL(KLON,NISP,KLEV*2) &
     &  , PGA(KLON,8,2,KLEV) , PGB(KLON,8,2,KLEV)
!
REAL_B :: PADJD(KLON,KLEV+1), PADJU(KLON,KLEV+1) &
     &  ,  PCNTRB(KLON,KLEV+1,KLEV+1) &
     &  ,  PDBDT(KLON,NISP,KLEV)
!
!-----------------------------------------------------------------------
!
!*       0.2   LOCAL ARRAYS
!              ------------
!
INTEGER_M :: ITX(KLON)
!
REAL_B :: ZGLAYD(KLON),ZGLAYU(KLON) &
     &  ,  ZTT(KLON,NTRA), ZTT1(KLON,NTRA), ZTT2(KLON,NTRA) &
     &  ,  ZUU(KLON,NUA)
!

!     LOCAL INTEGER SCALARS
INTEGER_M :: IBS, IDD, IM12, IMU, IND, INU, IXD, IXU,&
             &JA, JG, JK, JK1, JK2, JL, JNU

!     LOCAL REAL SCALARS
REAL_B :: ZWTR, ZWTR1, ZWTR2, ZWTR3, ZWTR4, ZWTR5, ZWTR6

!-----------------------------------------------------------------------
!
!*         1.    INITIALIZATION
!                --------------
!
!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------
!
DO JK = 1 , KLEV+1
  DO JL = KIDIA,KFDIA
    PADJD(JL,JK) = 0.
    PADJU(JL,JK) = 0.
  END DO
END DO
!
!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------
!
DO JA = 1 , NTRA
  DO JL = KIDIA,KFDIA
    ZTT (JL,JA) = 1.0
    ZTT1(JL,JA) = 1.0
    ZTT2(JL,JA) = 1.0
  END DO
END DO
!
DO JA = 1 , NUA
  DO JL = KIDIA,KFDIA
    ZUU(JL,JA) = 0.
 END DO
END DO
!
!     ------------------------------------------------------------------
!
!*         2.      VERTICAL INTEGRATION
!                  --------------------
!
!
!*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
!                  ---------------------------------
!
DO JK = 1 , KLEV
!
!*         2.1.1   DOWNWARD LAYERS
!                  ---------------
!
  IM12 = 2 * (JK - 1)
  IND = (JK - 1) * NG1P1 + 1
  IXD = IND
  INU = JK * NG1P1 + 1
  IXU = IND
!
  DO JL = KIDIA,KFDIA
    ZGLAYD(JL) = 0.
    ZGLAYU(JL) = 0.
  END DO
!
  DO JG = 1 , NG1
    IBS = IM12 + JG
    IDD = IXD + JG
    DO JA = 1 , KUAER
      DO JL = KIDIA,KFDIA
        ZUU(JL,JA) = PABCU(JL,JA,IND) - PABCU(JL,JA,IDD)
      END DO
    END DO
!
!
    CALL LWTT ( KIDIA,KFDIA,KLON &
    &          , PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT )
!
    
    DO JL = KIDIA,KFDIA
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10) &
      &    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11) &
      &    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12) &
      &    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13) &
      &    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14) &
      &    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYD(JL)=ZGLAYD(JL)+ZWTR*WG1(JG)
    END DO
!
!*         2.1.2   DOWNWARD LAYERS
!                  ---------------
!
    IMU = IXU + JG
    DO JA = 1 , KUAER
      DO JL = KIDIA,KFDIA
        ZUU(JL,JA) = PABCU(JL,JA,IMU) - PABCU(JL,JA,INU)
      END DO
    END DO
!
!
    CALL LWTT ( KIDIA,KFDIA,KLON &
    &          , PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT )
!
    DO JL = KIDIA,KFDIA
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10) &
      &    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11) &
      &    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12) &
      &    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13) &
      &    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14) &
      &    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15) 
      ZGLAYU(JL)=ZGLAYU(JL)+ZWTR*WG1(JG)
    END DO
!
  END DO
!
  DO JL = KIDIA,KFDIA
    PADJD(JL,JK) = ZGLAYD(JL)
    PCNTRB(JL,JK,JK+1) = ZGLAYD(JL)
    PADJU(JL,JK+1) = ZGLAYU(JL)
    PCNTRB(JL,JK+1,JK) = ZGLAYU(JL)
    PCNTRB(JL,JK  ,JK) = 0.0
  END DO
!
END DO
!
DO JK = 1 , KLEV
  JK2 = 2 * JK
  JK1 = JK2 - 1
  DO JNU = 1 , NISP
    DO JL = KIDIA,KFDIA
      PDBDT(JL,JNU,JK) = PDBSL(JL,JNU,JK1) + PDBSL(JL,JNU,JK2)
    END DO
  END DO
END DO
!
RETURN
END SUBROUTINE OLWVN
