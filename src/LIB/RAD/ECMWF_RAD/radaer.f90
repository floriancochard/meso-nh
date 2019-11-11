!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! ECMWF_RAD2_OLD2NEW 2003/02/19 13:36:50
!-----------------------------------------------------------------
!OPTIONS XOPT(HSFUN)
SUBROUTINE RADAER ( KIDIA , KFDIA , KLON , KTDIA , KLEV &
                  &, PAPRS , PTH ,&
                  PCVDAES ,PCVDAEL ,PCVDAEU ,PCVDAED, &
                  PAESEA,  PAELAN, PAEURB,PAEDES &
                  &, PAER                           )

!
!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------

!     J.-F. MAHFOUF    E.C.M.W.F.    97/07/04

!     Aadapted from

!     J.-J. MORCRETTE  E.C.M.W.F.    91/03/15

!     MODIFICATIONS.
!     --------------
!     J.-J. MORCRETTE  E.C.M.W.F.    93/03/15   OPERATIONAL CLIMATOLOGY
!     JJMorcrette  99-05-25     Revised aerosols
!-----------------------------------------------------------------------

#include "tsmbkind.h"

USE YOEAERD  , ONLY : CVDAES   ,CVDAEL   ,CVDAEU   ,CVDAED   ,&
            &RCAEOPS  ,RCAEOPL  ,RCAEOPU  ,RCAEOPD  ,RCTRBGA  ,&
            &RCVOBGA  ,RCSTBGA  ,RCTRPT   ,RAESC    ,RAESS    ,&
            &RAELC    ,RAELS    ,RAEUC    ,RAEUS    ,RAEDC    ,&
            &RAEDS



IMPLICIT NONE


!     DUMMY INTEGER SCALARS
INTEGER_M :: KFDIA
INTEGER_M :: KIDIA
INTEGER_M :: KLEV
INTEGER_M :: KLON
INTEGER_M :: KTDIA




!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

REAL_B :: PAPRS(KLON,KLEV+1)&
  &,  PTH  (KLON,KLEV+1),PCVDAES(KLON,KLEV+1) ,PCVDAEL(KLON,KLEV+1) ,PCVDAEU(KLON,KLEV+1) ,&
       PCVDAED(KLON,KLEV+1)

REAL_B :: PAESEA(KLON),PAELAN(KLON),PAEURB(KLON),PAEDES(KLON)

REAL_B :: PAER(KLON,KLEV,6)
!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------



REAL_B :: ZDPN (KLON) ,ZDPO (KLON)

REAL_B :: ZAEQSN(KLON),ZAEQSO(KLON),ZAEQLN(KLON),ZAEQLO(KLON)
REAL_B :: ZAEQUN(KLON),ZAEQUO(KLON),ZAEQDN(KLON),ZAEQDO(KLON)
REAL_B :: ZAETRN(KLON),ZAETRO(KLON)


!     LOCAL INTEGER SCALARS
INTEGER_M ::   JK, JL

!     LOCAL REAL SCALARS
REAL_B :: ZAETR, ZDPNMO


!     ------------------------------------------------------------------

!*         1.     AEROSOL PARAMETERS COMPUTATIONS
!                 -------------------------------

! Aerosols come from files containing the climatologies.

!*         2.      VERTICAL DISTRIBUTION
!*                 ---------------------

DO JL=KIDIA,KFDIA
  ZDPO(JL)=PAPRS (JL,1)
  
  ZAEQSO(JL)=RCAEOPS*PAESEA(JL)*PCVDAES(JL,1)
  ZAEQLO(JL)=RCAEOPL*PAELAN(JL)*PCVDAEL(JL,1)
  ZAEQUO(JL)=RCAEOPU*PAEURB(JL)*PCVDAEU(JL,1)
  ZAEQDO(JL)=RCAEOPD*PAEDES(JL)*PCVDAED(JL,1)
  ZAETRO(JL)=_ONE_
  
ENDDO
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZDPN(JL)=PAPRS (JL,JK+1)
    
    ZAEQSN(JL)=RCAEOPS*PAESEA(JL)*PCVDAES(JL,JK+1)
    ZAEQLN(JL)=RCAEOPL*PAELAN(JL)*PCVDAEL(JL,JK+1)
    ZAEQUN(JL)=RCAEOPU*PAEURB(JL)*PCVDAEU(JL,JK+1)
    ZAEQDN(JL)=RCAEOPD*PAEDES(JL)*PCVDAED(JL,JK+1)

    IF (_HALF_*(PAPRS(JL,JK)+PAPRS(JL,JK+1)) < 999._JPRB) THEN
! for models with top above 10hPa
      ZAETRN(JL)=_ONE_
      ZAETRO(JL)=_ONE_
 
    ELSE
    
     ZAETRN(JL)=ZAETRO(JL)*(MIN(_ONE_,PTH(JL,JK)/PTH(JL,JK+1)))**RCTRPT
    ENDIF

    ZAETR=SQRT  (ZAETRN(JL)*ZAETRO(JL))
   
    ZDPNMO    =ZDPN(JL)-ZDPO(JL)
!
    PAER(JL,JK,1)=(_ONE_-ZAETR)*(RCTRBGA*ZDPNMO+ ZAEQLN(JL)-ZAEQLO(JL))
!
  
    PAER(JL,JK,2)=(_ONE_-ZAETR)*(ZAEQSN(JL)-ZAEQSO(JL))
    PAER(JL,JK,3)=(_ONE_-ZAETR)*(ZAEQDN(JL)-ZAEQDO(JL))
    PAER(JL,JK,4)=(_ONE_-ZAETR)*(ZAEQUN(JL)-ZAEQUO(JL))
    PAER(JL,JK,5)= ZAETR * RCVOBGA*ZDPNMO
    PAER(JL,JK,6)= ZAETR * RCSTBGA*ZDPNMO
!         AADS(JL,JK)=MAX(RCAEADM, (RCAEADK(1)*PAER(JL,JK,1)
!           + RCAEADK(2)*PAER(JL,JK,2)+RCAEADK(3)*PAER(JL,JK,3))/ZDPNMO)
    
!**** **************************************************
!**** **************************************************
  ENDDO
  DO JL=KIDIA,KFDIA
    ZDPO(JL)=ZDPN(JL)
   
    ZAEQSO(JL)=ZAEQSN(JL)
    ZAEQLO(JL)=ZAEQLN(JL)
    ZAEQUO(JL)=ZAEQUN(JL)
    ZAEQDO(JL)=ZAEQDN(JL)
    ZAETRO(JL)=ZAETRN(JL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------
!MODIF essai sortie AER
!!$open (unit=99, file='testradaer',status='replace')
!!$write(99,*) 'Aerosols 1=Continental   2=Maritime   3=Desert     4=Urban     5=Volcanic  6=Stratos.Bckgnd'
!!$    DO jk=1,KLEV
!!$      write(99,*) jk,(PAER(JL,JK,:))
!!$9106     FORMAT(1x,'Aeros. ',I3,6E13.5)
!!$    ENDDO   
!!$close(99)

RETURN

END SUBROUTINE RADAER






