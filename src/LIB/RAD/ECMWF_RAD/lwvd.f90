!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2 2003/02/19 13:36:35
!-----------------------------------------------------------------
SUBROUTINE LWVD &
 &( KIDIA,  KFDIA, KLON , KLEV  , KTRAER &
 &, PABCU,  PDBDT &
 &, PGA  ,  PGB &
 &, PCNTRB, PDISD, PDISU, PDWFSU &
 &)

!**** *LWVD*   - L.W., VERTICAL INTEGRATION, DISTANT LAYERS

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU  : (KLON,NUA,3*KLEV+1) ; ABSORBER AMOUNTS
! PDBDT  : (KLON,KLEV)         ; LAYER PLANCK FUNCTION GRADIENT
! PGA, PGB                     ; PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PCNTRB : (KLON,KLEV+1,KLEV+1); ENERGY EXCHANGE MATRIX
! PDIS.. : (KLON,KLEV+1)       ; CONTRIBUTION BY DISTANT LAYERS
! PDWFSU : (KLON,NSIL)         ; SPECTRAL DOWNWARD FLUX AT SURFACE

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
!     CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE

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
!        JJ Morcrette 97-04-18 Revised continuum + Surf. Emissiv.
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOELW    , ONLY : NSIL     ,NIPD     ,NTRA     ,NUA      ,NG1P1
!
USE MODI_LWTTm
!
IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KTRAER



!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------


REAL_B :: PABCU(KLON,NUA,3*KLEV+1)&
  &,  PDBDT(KLON,NSIL,KLEV)&
  &,  PGA(KLON,NIPD,2,KLEV)   , PGB(KLON,NIPD,2,KLEV)

REAL_B :: PCNTRB(KLON,KLEV+1,KLEV+1)&
  &,  PDISD(KLON,KLEV+1)      , PDISU(KLON,KLEV+1)&
  &,  PDWFSU(KLON,NSIL)

!-----------------------------------------------------------------------

!*       0.2   LOCAL ARRAYS
!              ------------

REAL_B :: ZTT(KLON,NTRA), ZTT1(KLON,NTRA), ZTT2(KLON,NTRA)

!     LOCAL INTEGER SCALARS
INTEGER_M :: IJKL, IKD1, IKD2, IKJ, IKJP1, IKM1, IKN,&
             &IKP1, IKU1, IKU2, ITT, JA, JK, JKJ, JL, JLK

!     LOCAL REAL SCALARS
REAL_B :: ZWW, ZWW1, ZWW2, ZWW3, ZWW4, ZWW5, ZWW6


!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

DO JK = 1, KLEV+1
  DO JL = KIDIA,KFDIA
    PDISD(JL,JK) = _ZERO_
    PDISU(JL,JK) = _ZERO_
  ENDDO
ENDDO

!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------

DO JA = 1, NTRA
  DO JL = KIDIA,KFDIA
    ZTT (JL,JA) = _ONE_
    ZTT1(JL,JA) = _ONE_
    ZTT2(JL,JA) = _ONE_
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------


!*         2.2     CONTRIBUTION FROM DISTANT LAYERS
!                  ---------------------------------


!*         2.2.1   DISTANT AND ABOVE LAYERS
!                  ------------------------


!*         2.2.2   FIRST UPPER LEVEL
!                  -----------------

DO JK = 1 , KLEV-1
  IKP1=JK+1
  IKN=(JK-1)*NG1P1+1
  IKD1= JK  *NG1P1+1

  CALL LWTTM &
   &( KIDIA         , KFDIA          , KLON &
   &, PGA(1,1,1,JK) , PGB(1,1,1,JK)&
   &, PABCU(1,1,IKN), PABCU(1,1,IKD1), ZTT1 &
   &)



!*         2.2.3   HIGHER UP
!                  ---------

  ITT=1
  DO JKJ=IKP1,KLEV
    IF(ITT == 1) THEN
      ITT=2
    ELSE
      ITT=1
    ENDIF
    IKJP1=JKJ+1
    IKD2= JKJ  *NG1P1+1

    IF(ITT == 1) THEN
      CALL LWTTM &
       &( KIDIA         , KFDIA          , KLON &
       &, PGA(1,1,1,JKJ), PGB(1,1,1,JKJ)&
       &, PABCU(1,1,IKN), PABCU(1,1,IKD2), ZTT1 &
       &)


    ELSE
      CALL LWTTM &
       &( KIDIA         , KFDIA          , KLON &
       &, PGA(1,1,1,JKJ), PGB(1,1,1,JKJ)&
       &, PABCU(1,1,IKN), PABCU(1,1,IKD2), ZTT2 &
       &)


    ENDIF

    DO JA = 1, KTRAER
      DO JL = KIDIA,KFDIA
        ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*_HALF_
      ENDDO
    ENDDO

    DO JL = KIDIA,KFDIA
      ZWW1=PDBDT(JL,1,JKJ)*ZTT(JL,1)          *ZTT(JL,10)
      ZWW2=PDBDT(JL,2,JKJ)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
      ZWW3=PDBDT(JL,3,JKJ)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
      ZWW4=PDBDT(JL,4,JKJ)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
      ZWW5=PDBDT(JL,5,JKJ)*ZTT(JL,3)          *ZTT(JL,14)
      ZWW6=PDBDT(JL,6,JKJ)*ZTT(JL,6)          *ZTT(JL,15)
      ZWW=ZWW1+ZWW2+ZWW3+ZWW4+ZWW5+ZWW6
      PDISD(JL,JK)=PDISD(JL,JK)+ZWW
      PCNTRB(JL,JK,IKJP1)=ZWW
      IF (JK == 1) THEN
        PDWFSU(JL,1)=PDWFSU(JL,1)+ZWW1
        PDWFSU(JL,2)=PDWFSU(JL,2)+ZWW2
        PDWFSU(JL,3)=PDWFSU(JL,3)+ZWW3
        PDWFSU(JL,4)=PDWFSU(JL,4)+ZWW4
        PDWFSU(JL,5)=PDWFSU(JL,5)+ZWW5
        PDWFSU(JL,6)=PDWFSU(JL,6)+ZWW6
      ENDIF
    ENDDO


  ENDDO
ENDDO


!*         2.2.4   DISTANT AND BELOW LAYERS
!                  ------------------------


!*         2.2.5   FIRST LOWER LEVEL
!                  -----------------

DO JK=3,KLEV+1
  IKN=(JK-1)*NG1P1+1
  IKM1=JK-1
  IKJ=JK-2
  IKU1= IKJ  *NG1P1+1


  CALL LWTTM &
   &( KIDIA          , KFDIA         , KLON &
   &, PGA(1,1,1,IKJ) , PGB(1,1,1,IKJ)&
   &, PABCU(1,1,IKU1), PABCU(1,1,IKN), ZTT1 &
   &)




!*         2.2.6   DOWN BELOW
!                  ----------

  ITT=1
  DO JLK=1,IKJ
    IF(ITT == 1) THEN
      ITT=2
    ELSE
      ITT=1
    ENDIF
    IJKL=IKM1-JLK
    IKU2=(IJKL-1)*NG1P1+1


    IF(ITT == 1) THEN
      CALL LWTTM &
       &( KIDIA          , KFDIA          , KLON &
       &, PGA(1,1,1,IJKL), PGB(1,1,1,IJKL)&
       &, PABCU(1,1,IKU2), PABCU(1,1,IKN) , ZTT1 &
       &)


    ELSE
      CALL LWTTM &
       &( KIDIA          , KFDIA          , KLON &
       &, PGA(1,1,1,IJKL), PGB(1,1,1,IJKL)&
       &, PABCU(1,1,IKU2), PABCU(1,1,IKN) , ZTT2 &
       &)


    ENDIF

    DO JA = 1, KTRAER
      DO JL = KIDIA,KFDIA
        ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*_HALF_
      ENDDO
    ENDDO

    DO JL = KIDIA,KFDIA
      ZWW=PDBDT(JL,1,IJKL)*ZTT(JL,1)          *ZTT(JL,10)&
       &+PDBDT(JL,2,IJKL)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)&
       &+PDBDT(JL,3,IJKL)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)&
       &+PDBDT(JL,4,IJKL)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)&
       &+PDBDT(JL,5,IJKL)*ZTT(JL,3)          *ZTT(JL,14)&
       &+PDBDT(JL,6,IJKL)*ZTT(JL,6)          *ZTT(JL,15)
      PDISU(JL,JK)=PDISU(JL,JK)+ZWW
      PCNTRB(JL,JK,IJKL)=ZWW
    ENDDO


  ENDDO
ENDDO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE LWVD
