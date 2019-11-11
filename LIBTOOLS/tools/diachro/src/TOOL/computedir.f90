!     ######spl
      MODULE MODI_COMPUTEDIR
!     #########################
!
INTERFACE
!
SUBROUTINE COMPUTEDIR(KITER,KJTER,KIUB1,KIUB2,KISKIP,PDIRU,PDIRV,PLO)
!
INTEGER           :: KITER, KJTER, KIUB1, KIUB2, KISKIP
REAL,DIMENSION(:,:)         :: PDIRU, PDIRV
REAL,DIMENSION(:,:),OPTIONAL ::PLO
!
END SUBROUTINE COMPUTEDIR
!
END INTERFACE
!
END MODULE MODI_COMPUTEDIR
!
!     #################
      SUBROUTINE COMPUTEDIR(KITER,KJTER,KIUB1,KIUB2,KISKIP,PDIRU,PDIRV,PLO)
!     #################
!
!!****  *COMPUTEDIR* - 
!!                                                            
!!
!!    PURPOSE
!!    -------
!        Trace PH (tableaux 1D scalaires  y compris MUMVM et DIRUMVM)
!        dans traceh_fordiachro
!
!!**  METHOD
!!    ------
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
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       30/11/01
!!      Updated   PM  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_GRIDPROJ
USE MODD_RESOLVCAR, ONLY: LCV, NVERBIA
USE MODD_GRID, ONLY: XBETA, XRPK, XLON0
USE MODD_COORD, ONLY: XDSX, XDSY
USE MODD_GRID1, ONLY: XXHAT, XYHAT
USE MODD_GRID, ONLY: XLATORI, XLONORI
USE MODN_NCAR, ONLY: XSPVAL
!
IMPLICIT NONE
!
!
COMMON/TEMH/XZZX,XZZY,NIIMAX,NIJMAX
#include "big.h"
INTEGER              :: NIIMAX,NIJMAX
REAL,DIMENSION(N2DVERTX) :: XZZX
REAL,DIMENSION(400)  :: XZZY
!
!
!*       0.1   Dummy arguments
!
INTEGER           :: KITER, KJTER, KIUB1, KIUB2, KISKIP
REAL,DIMENSION(:,:)         :: PDIRU, PDIRV
REAL,DIMENSION(:,:),OPTIONAL ::PLO
!
!*       0.1   Local variables
!
!
INTEGER           :: JILOOP, JJLOOP
!
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZX, ZYY,ZLAT,ZLON,ZLO
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: ZZY
REAL              :: ZRPK, ZBETA, ZLON0
!
!
!-------------------------------------------------------------------------------
!
!*      1. 
!              ----------------------------
!
if(nverbia > 0)then
print *,' **entree computedir KIUB1,KIUB2,KITER,KJTER,KISKIP ',KIUB1,KIUB2,KITER,KJTER,KISKIP
print *,' **entree computedir size(PDIRU) ',size(PDIRU,1),size(PDIRU,2)
print *,' **entree computedir PDIRU PDIRV ',PDIRU,PDIRV
endif
!
ALLOCATE(ZLO(KITER,KJTER))
IF (PRESENT (PLO) ) THEN
   if(nverbia > 0)then
     print *,' **computedir : utilisation du lat lon passe en argument'
   endif
   ZLO=PLO
ELSE
 ! calcule ZLO en fonction de XXHAT et XYHAT
!! Supprime en nov 2001 Appel routine COMPUTEDIR
  ALLOCATE(ZX(KITER,1),ZZY(KJTER))
  IF(LCV)THEN
    ZX(:,1)=XDSX(1:KIUB1:KISKIP,1)
  ELSE
    ZX(:,1)=XZZX(1:KIUB1:KISKIP)
    ZZY=XZZY(1:KIUB2:KISKIP)
  ENDIF
  ALLOCATE(ZYY(KITER,1),ZLAT(KITER,1),ZLON(KITER,1))
  DO JJLOOP=1,KJTER
    DO JILOOP=1,KITER
      IF(LCV)THEN
        ZYY(JILOOP,1)=XDSY(JILOOP,1)
      ELSE
        ZYY(JILOOP,1)=ZZY(JJLOOP)
      ENDIF
    ENDDO
    CALL SM_LATLON_A(XLATORI,XLONORI,ZX,ZYY,ZLAT,ZLON)
    ZLO(:,JJLOOP)=ZLON(:,1)
  ENDDO
!if(nverbia > 0)then
!print *,' **computedir  ZX,ZZY,ZYY ',ZX,ZZY,ZYY
!endif
ENDIF
! fin de if (PRESENT (PLO) )
!
if(nverbia > 0)then
print*,'** computedir LO ',KITER,KJTER,ZLO
endif

where(PDIRU /= xspval .AND. PDIRV /= xspval)
    PDIRU=ATAN2(PDIRV,PDIRU)*180./ACOS(-1.)
endwhere
!if(nverbia > 0)then
!print *,' **computedir  PDIRU EN DEG. ',PDIRU
!endif
if(nverbia > 0)then
  print *,' PDIRU 1,1 KITER/2,1 1,KJTER/2 KITER/2,KJTER/2 KITER,KJTER 22,29 '
  print *,PDIRU(1,1),  PDIRU(KITER/2,1), PDIRU(1,KJTER/2), PDIRU(KITER/2,KJTER/2), &
  PDIRU(KITER,KJTER),PDIRU(22,29)
endif
!
ZRPK=XRPK
ZBETA=XBETA
ZLON0=XLON0
where(PDIRU /= xspval .AND. PDIRV /= xspval)
  PDIRU=PDIRU - (ZRPK*(ZLO-ZLON0)-ZBETA) + 90.
endwhere
!
!if(nverbia > 0)then
!print *,' **computedir  PDIRU suite ',PDIRU
!print *,' **computedir  ZRPK,ZBETA,ZLON0 ',ZRPK,ZBETA,ZLON0
!endif
! 
WHERE(PDIRU < 0.)PDIRU=PDIRU+360.
WHERE(PDIRU > 360. .AND. PDIRU /= XSPVAL)PDIRU=PDIRU-360.
if(nverbia > 0)then
   print *,' PDIRU 1,1 KITER/2,1 1,KJTER/2 KITER/2,KJTER/2 KITER,KJTER '
   print *,PDIRU(1,1),  PDIRU(KITER/2,1), PDIRU(1,KJTER/2), PDIRU(KITER/2,KJTER/2), &
   PDIRU(KITER,KJTER)
endif
!
where(PDIRU /= xspval .AND. PDIRV /= xspval)
  PDIRV=360.-PDIRU
elsewhere
  PDIRV=XSPVAL
endwhere
!if(nverbia > 0)then
!print *,' **computedir  PDIRV EN DEG. ',PDIRV
!endif
if(nverbia > 0)then
  print *,' PDIRV 1,1 KITER/2,1 1,KJTER/2 KITER/2,KJTER/2 KITER,KJTER '
  print *,PDIRV(1,1),  PDIRV(KITER/2,1), PDIRV(1,KJTER/2), PDIRV(KITER/2,KJTER/2), &
  PDIRV(KITER,KJTER)
endif
!! Supprime en nov 2001 Appel routine COMPUTEDIR
IF (PRESENT (PLO) ) THEN
  DEALLOCATE(ZLO)
ELSE
  DEALLOCATE(ZX,ZZY,ZYY,ZLAT,ZLON,ZLO)
ENDIF
!
!------------------------------------------------------------------------------
!
!*      2.    EXIT
!             ----
RETURN
!
END SUBROUTINE COMPUTEDIR
