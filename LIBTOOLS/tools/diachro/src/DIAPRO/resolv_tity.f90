!     ######spl
      MODULE MODI_RESOLV_TITY
!     ######################
!
INTERFACE
!
SUBROUTINE RESOLV_TITY(HTIT,PVL,PVR,PVB,PVT,HOUT)
REAL             :: PVL, PVR, PVB, PVT
CHARACTER(LEN=*)  :: HTIT, HOUT
END SUBROUTINE RESOLV_TITY
!
END INTERFACE
END MODULE MODI_RESOLV_TITY
!     ######spl
      SUBROUTINE RESOLV_TITY(HTIT,PVL,PVR,PVB,PVT,HOUT)
!     #################################################
!
!!****  *RESOLV_TITY* - 
!!
!!    PURPOSE
!!    -------
!      
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
!!      Module
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       24/11/95
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
USE MODD_TIT

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

REAL             :: PVL, PVR, PVB, PVT
CHARACTER(LEN=*) :: HTIT, HOUT
!
!*       0.1   Local variables
!              ---------------

CHARACTER(LEN=LEN(HOUT)) :: YTEM
INTEGER                  :: ILEN, INBV, J, IM
INTEGER,DIMENSION(10)    :: IJM
REAL                     :: ZSIZC, ZM
REAL              :: ZXPOSTITYT, ZXYPOSTITYT
REAL              :: ZXPOSTITYM, ZXYPOSTITYM
REAL              :: ZXPOSTITYB, ZXYPOSTITYB

!
!------------------------------------------------------------------------------
YTEM=HOUT
IF(.NOT.LTITDEF)THEN
SELECT CASE(HTIT)
  CASE('CTITYT')
    IF(CTITYT == 'WHITE' .OR. CTITYT == 'white' .OR.  &
      CTITYT == 'BLANC' .OR. CTITYT == 'blanc')THEN
      YTEM(1:LEN(YTEM))=' '
      CTITALL='NODEFAULT'
    ELSE  IF(CTITYT == 'DEFAULT' .OR.   CTITYT == 'default' .OR. &
             CTITYT == 'DEFAUT' .OR.  CTITYT == 'defaut')THEN
    ELSE
      YTEM=CTITYT
      CTITALL='NODEFAULT'
    ENDIF
  CASE('CTITYM')
    IF(CTITYM == 'WHITE' .OR. CTITYM == 'white' .OR.  &
      CTITYM == 'BLANC' .OR. CTITYM == 'blanc')THEN
      YTEM(1:LEN(YTEM))=' '
      CTITALL='NODEFAULT'
    ELSE  IF(CTITYM == 'DEFAULT' .OR.   CTITYM == 'default' .OR. &
             CTITYM == 'DEFAUT' .OR.  CTITYM == 'defaut')THEN
    ELSE
      YTEM=CTITYM
      CTITALL='NODEFAULT'
    ENDIF
  CASE('CTITYB')
    IF(CTITYB == 'WHITE' .OR. CTITYB == 'white' .OR.  &
      CTITYB == 'BLANC' .OR. CTITYB == 'blanc')THEN
      YTEM(1:LEN(YTEM))=' '
      CTITALL='NODEFAULT'
    ELSE  IF(CTITYB == 'DEFAULT' .OR.   CTITYB == 'default' .OR. &
             CTITYB == 'DEFAUT' .OR.  CTITYB == 'defaut')THEN
    ELSE
      YTEM=CTITYB
      CTITALL='NODEFAULT'
    ENDIF
END SELECT
ENDIF
YTEM=ADJUSTL(YTEM)
ILEN=LEN_TRIM(YTEM)
IJM=0
INBV=1; IJM(INBV)=0
DO J =1,ILEN
  IF(YTEM(J:J) == ';')THEN
    INBV=INBV+1
    IJM(INBV)=J
  ENDIF
ENDDO
INBV=INBV+1
IJM(INBV)=ILEN+1
ZSIZC=(.9-.1)/50.
!ZSIZC=(PVT-PVB)/50.
print*,PVL,PVT,PVR,PVB
DO J=2,INBV
SELECT CASE(HTIT)
  CASE('CTITYT')
     IF (L90TITYT) THEN
     ZXPOSTITYT=MAX(PVL-0.03,0.)
     ZXYPOSTITYT=PVT-J*ZSIZC
     IF(XPOSTITYT /= 0.)THEN
       ZXPOSTITYT=XPOSTITYT
     ENDIF
     IF(XYPOSTITYT /= 0.)THEN
       ZXYPOSTITYT=XYPOSTITYT
     ENDIF
      IF(XSZTITYT /= 0.)THEN       
        CALL PLCHHQ(ZXPOSTITYT,ZXYPOSTITYT,YTEM(IJM(J-1)+1:IJM(J)-1),XSZTITYT,90.,0.)
      ELSE
        CALL PLCHHQ(ZXPOSTITYT,ZXYPOSTITYT,YTEM(IJM(J-1)+1:IJM(J)-1),ZSIZC/2.,90.,0.)
      ENDIF
     ELSE
     ZXPOSTITYT=MAX(PVL-.12,0.)
     ZXYPOSTITYT=PVT-J*ZSIZC  
     IF(XPOSTITYT /= 0.)THEN
       ZXPOSTITYT=XPOSTITYT
     ENDIF
     IF(XYPOSTITYT /= 0.)THEN
       ZXYPOSTITYT=XYPOSTITYT
     ENDIF

      IF(XSZTITYT /= 0.)THEN       
        CALL PLCHHQ(ZXPOSTITYT,ZXYPOSTITYT,YTEM(IJM(J-1)+1:IJM(J)-1),XSZTITYT,0.,-1.)
      ELSE
         CALL PLCHHQ(ZXPOSTITYT,ZXYPOSTITYT,YTEM(IJM(J-1)+1:IJM(J)-1),ZSIZC/2.,0.,-1.)
      ENDIF
     ENDIF
  CASE('CTITYM')
    ZM=(PVB+PVT)/2.
    IM=(INBV-1)/2
    IF(IM /= 0)THEN
      IM=INBV-1-IM-J
    ENDIF

    IF (L90TITYM) THEN
     ZXPOSTITYM=MAX(PVL-0.03,0.)
     ZXYPOSTITYM=ZM+IM*ZSIZC
     IF(XPOSTITYM /= 0.)THEN
       ZXPOSTITYM=XPOSTITYM
     ENDIF
     IF(XYPOSTITYM /= 0.)THEN
       ZXYPOSTITYM=XYPOSTITYM
     ENDIF
        
      IF(XSZTITYM /= 0.)THEN      
         CALL PLCHHQ(ZXPOSTITYM,ZXYPOSTITYM,YTEM(IJM(J-1)+1:IJM(J)-1),XSZTITYM,90.,0.)
      ELSE
         CALL PLCHHQ(ZXPOSTITYM,ZXYPOSTITYM,YTEM(IJM(J-1)+1:IJM(J)-1),ZSIZC/2.,90.,0.)
      ENDIF      
    ELSE
     ZXPOSTITYM=MAX(PVL-.12,0.)
     ZXYPOSTITYM=ZM+IM*ZSIZC
     IF(XPOSTITYM /= 0.)THEN
       ZXPOSTITYM=XPOSTITYM
     ENDIF
     IF(XYPOSTITYM /= 0.)THEN
       ZXYPOSTITYM=XYPOSTITYM
     ENDIF     
      IF(XSZTITYM /= 0.)THEN      
          CALL PLCHHQ(ZXPOSTITYM,ZXYPOSTITYM,YTEM(IJM(J-1)+1:IJM(J)-1),XSZTITYM,0.,-1.)
      ELSE
          CALL PLCHHQ(ZXPOSTITYM,ZXYPOSTITYM,YTEM(IJM(J-1)+1:IJM(J)-1),ZSIZC/2.,0.,-1.)
      ENDIF       
    ENDIF
  CASE('CTITYB')
    IF (L90TITYB) THEN
     ZXPOSTITYB=MAX(PVL-0.03,0.)
     ZXYPOSTITYB=PVB+(INBV-J)*ZSIZC
     IF(XPOSTITYB /= 0.)THEN
       ZXPOSTITYB=XPOSTITYB
     ENDIF
     IF(XYPOSTITYB /= 0.)THEN
       ZXYPOSTITYB=XYPOSTITYB
     ENDIF       
      IF(XSZTITYB /= 0.)THEN      
        CALL PLCHHQ(ZXPOSTITYB,ZXYPOSTITYB,YTEM(IJM(J-1)+1:IJM(J)-1),XSZTITYB,90.,0.)
      ELSE
        CALL PLCHHQ(ZXPOSTITYB,ZXYPOSTITYB,YTEM(IJM(J-1)+1:IJM(J)-1),ZSIZC/2.,90.,0.)
      ENDIF       
    ELSE
     ZXPOSTITYB=MAX(PVL-.12,0.)
     ZXYPOSTITYB=PVB+(INBV-J)*ZSIZC
     IF(XPOSTITYB /= 0.)THEN
       ZXPOSTITYB=XPOSTITYB
     ENDIF
     IF(XYPOSTITYB /= 0.)THEN
       ZXYPOSTITYB=XYPOSTITYB
     ENDIF           
      IF(XSZTITYB /= 0.)THEN      
         CALL PLCHHQ(ZXPOSTITYB,ZXYPOSTITYB,YTEM(IJM(J-1)+1:IJM(J)-1),XSZTITYB,0.,-1.)
      ELSE
         CALL PLCHHQ(ZXPOSTITYB,ZXYPOSTITYB,YTEM(IJM(J-1)+1:IJM(J)-1),ZSIZC/2.,0.,-1.)
      ENDIF       
    ENDIF
END SELECT
ENDDO
RETURN
END SUBROUTINE RESOLV_TITY
