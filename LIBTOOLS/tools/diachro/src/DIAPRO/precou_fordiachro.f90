!     ######spl
      MODULE MODI_PRECOU_FORDIACHRO
!     #############################
!
INTERFACE
!
SUBROUTINE PRECOU_FORDIACHRO(PWORK3D,PTEMCV)
REAL,DIMENSION(:,:,:)  :: PWORK3D
REAL,DIMENSION(:,:)    :: PTEMCV
END SUBROUTINE PRECOU_FORDIACHRO
!
END INTERFACE
!
END MODULE MODI_PRECOU_FORDIACHRO

      SUBROUTINE PRECOU_FORDIACHRO(PWORK3D,PTEMCV)
!     ############################################
!
!!****  *PRECOU_FORDIACHRO* - Preliminary calculation for vertical cross-sections of
!!****             basis set prognostic Meso-NH variables
!!
!!    PURPOSE
!!    -------
!!      
!       When a verical cross-section is requested, this routine allocates
!     2D work arrays to to store the interpolated fields produced by the
!     COUPE routine. 
!
!!**  METHOD
!!    ------
!!      Array allocation and call to the COUPE vertical plane interpolator 
!!
!!      WARNING: This program section is exceptionally boring, 
!!               I fell asleep twice updating it.
!!
!!    EXTERNAL
!!    --------
!!      COUPE    : interpolates the model data onto the vertical 
!!                 cross-section plane requested by the user.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module MODN_PARA: Defines NAM_DOMAIN_POS namelist (former PARA common)
!!          NLMAX            :  Number of points horizontally along
!!                              the vertical section
!!          Module MODD_DIM1 : contains dimensions of data arrays
!!              NKMAX      : z array dimension
!!
!!     Module MODD_CVERT:  Declares work arrays for vertical cross-sections
!!          XWORKZ   : working array for true altitude storage (all grids)
!!          XWZ      : working array for topography (all grids)
!!
!!      Module MODD_OUT    : Defines a log. unit for printing
!!          NIMAXT   :  Size of the displayed window within a
!!          NJMAXT   :                MESO-NH field arrays
!!
!!     Module MODD_PARAMETERS :  Contains array border depths
!!          JPVEXT   : Vertical external points number
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   15/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
! modules MesoNH
USE MODD_CONF, ONLY: L2D
USE MODD_DIM1, ONLY: NIMAX,NJMAX,NKMAX,NIINF,NISUP,NJINF,NJSUP
USE MODD_GRID1, ONLY: XZZ
! modules diaprog
USE MODN_PARA
USE MODD_TYPE_AND_LH
USE MODD_NMGRID
USE MODN_NCAR
USE MODD_CVERT
USE MODD_NMGRID
USE MODD_PARAMETERS
USE MODD_RESOLVCAR
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_ALLOC_FORDIACHRO
USE MODD_PVT
USE MODD_MEMGRIUV
USE MODI_COMPUTEDIR

IMPLICIT NONE
!
!*       0.1    Interface declarations
!
INTERFACE
      SUBROUTINE COUPE_FORDIACHRO(PTABI,PTABO,K)
      REAL,DIMENSION(:,:)      :: PTABI
      REAL,DIMENSION(:)        :: PTABO
      INTEGER :: K
      END SUBROUTINE COUPE_FORDIACHRO  
END INTERFACE
INTERFACE
      SUBROUTINE ROTA(PTEM1,PTEMV)
      REAL, DIMENSION(:,:),  INTENT(INOUT) :: PTEM1
      REAL, DIMENSION(:,:),  INTENT(INOUT) :: PTEMV
      END SUBROUTINE ROTA
END INTERFACE
INTERFACE
      SUBROUTINE COUPEUW_FORDIACHRO(PTABI,PTABO,K,KCOMP)
      REAL,DIMENSION(:,:)      :: PTABI
      REAL,DIMENSION(:)        :: PTABO
      INTEGER                  ::  K    
      INTEGER                  ::  KCOMP 
      END SUBROUTINE COUPEUW_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE ROTAUW(PTEM1,PTEMV)
      REAL, DIMENSION(:),  INTENT(INOUT) :: PTEM1
      REAL, DIMENSION(:),  INTENT(INOUT) :: PTEMV
      END SUBROUTINE ROTAUW
END INTERFACE
!
COMMON/TEMH/XZZX,XZZY,NIIMAX,NIJMAX
#include "big.h"
REAL,DIMENSION(N2DVERTX) :: XZZX
REAL,DIMENSION(N2DVERTX) :: XZZY
INTEGER :: NIIMAX, NIJMAX

!
!*      0.12    Dummy arguments
!
REAL,DIMENSION(:,:,:)  :: PWORK3D
REAL,DIMENSION(:,:)    :: PTEMCV
!
!*      0.2     Local variables
!
INTEGER :: IIU,IJU,IKU, JKLOOP, IKB, IKE, IWKU
INTEGER :: IUI, IUJ
INTEGER :: ITER, JTER, IUB1, IUB2, ISKIP
INTEGER,SAVE :: IPRESM, ITPRESY
!
!
REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: ZWORK3D, ZWORK3W
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: ZTEM1, ZTEMV, ZTEMW
REAL,DIMENSION(:),ALLOCATABLE,SAVE   :: ZTEM2, ZTEMVR, ZTEMWR
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE     :: ZX
REAL,DIMENSION(:),ALLOCATABLE,SAVE       :: ZZY
!
!-----------------------------------------------------------------------------
!
!*       1.      SETS ARRAY SIZES AND ALLOCATES ARRAYS
!                -------------------------------------
!
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
IKB=1+JPVEXT
IKE=IKU-JPVEXT
IWKU=SIZE(PWORK3D,3)
!
! Dedicated work arrays for vertical cross sections; last index is
! NMGRID grid selector.
! XWORZ contains true altitudes, for all grids
! XWZ   contains topography, for all grids 
!
if(nverbia > 0)then
  print *,' **precou IKU AV ALLOCATE(XWORKZ NLMAX ',IKU,NLMAX
endif
IF(ALLOCATED(XWORKZ))THEN
  IF (SIZE(XWORKZ,1) /= NLMAX)THEN
    DEALLOCATE(XWORKZ)
    ALLOCATE(XWORKZ(NLMAX,IKU,7))
  ENDIF
ELSEIF(.NOT.ALLOCATED(XWORKZ))THEN
!ELSE
  ALLOCATE(XWORKZ(NLMAX,IKU,7))
if(nverbia > 0)then
  print *,' **precou IKU AP ALLOCATE(XWORKZ NLMAX ',IKU,NLMAX
endif
ENDIF
if(nverbia > 0)then
  print *,' **precou IKU AV ALLOCATE(XWZ NLMAX ',IKU,NLMAX
! print *,' **precou  ALLOCATE(XWZ size(XWZ,1)et 2 ',size(XWZ,1),size(XWZ,2)
endif
IF(ALLOCATED(XWZ))THEN
  IF(SIZE(XWZ,1) /= NLMAX)THEN
  DEALLOCATE(XWZ)
  ALLOCATE(XWZ(NLMAX,7))
  ENDIF
ELSE IF(.NOT.ALLOCATED(XWZ))THEN
  ALLOCATE(XWZ(NLMAX,7))
ENDIF
! Oct 2000 prise en compte PH issus du 2D horiz. 
! Volontairement place apres ALLOCATE XWORKZ sinon pb
IF(IWKU == 1)THEN
  IKB=1; IKE=1; IKU=1
if(nverbia > 0)then
  print *,' **precou IKU AP ALLOCATE(XWORKZ NLMAX ',IKU,NLMAX
  print *,' **precou  sizePTEMCV ',size(PTEMCV,1),size(PTEMCV,2)
endif
ENDIF
!
! Local work arrays
!
ALLOCATE(ZTEM1(1:IIU,1:IJU))
!ALLOCATE(ZTEM1(1:NIH-NIL+1,1:NJH-NJL+1))
ALLOCATE(ZTEM2(NLMAX))
! Janvier 2001 + LDIRWIND et LUMVM et LUTVT et LSUMVM et LSUTVT
IF(LULM .OR. LULT .OR.LVTM .OR. LVTT .OR. LULMWM .OR. LULTWT .OR. &
   LMUMVM .OR. LMUTVT .OR. LMLSUMVM .OR. LMLSUTVT .OR. LDIRWIND .OR. &
   !LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT)THEN
   LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT .OR. &
   (LDIRWT .AND. .NOT.LDIRWIND).OR.(LDIRWM .AND. .NOT.LDIRWIND) )THEN
  ALLOCATE(ZWORK3D(SIZE(PWORK3D,1),SIZE(PWORK3D,2),SIZE(PWORK3D,3)))
  ALLOCATE(ZTEMV(1:IIU,1:IJU))
ENDIF
IF(LULMWM .OR. LULTWT)THEN
  ALLOCATE(ZWORK3W(SIZE(PWORK3D,1),SIZE(PWORK3D,2),SIZE(PWORK3D,3)))
  ALLOCATE(ZTEMW(1:IIU,1:IJU))
  ALLOCATE(ZTEMVR(NLMAX),ZTEMWR(NLMAX))
  IF(ALLOCATED(XWCV))DEALLOCATE(XWCV)
  ALLOCATE(XWCV(SIZE(PTEMCV,1),SIZE(PTEMCV,2)))
ENDIF
! Janvier 2001 + LDIRWIND et LUMVM et LUTVT et LSUMVM et LSUTVT
!IF(LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT .OR. LDIRWIND)THEN
IF(LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT .OR. LDIRWIND .OR. &
   (LDIRWT .AND. .NOT.LDIRWIND).OR.(LDIRWM .AND. .NOT.LDIRWIND) )THEN
  ALLOCATE(ZTEMVR(NLMAX))
  IF(ALLOCATED(XWCV))DEALLOCATE(XWCV)
  ALLOCATE(XWCV(SIZE(PTEMCV,1),SIZE(PTEMCV,2)))
ENDIF
!
!------------------------------------------------------------------------------

XWORKZ(:,:,:)=0.
XWZ(:,:)=0.
PTEMCV=XSPVAL
IF(ALLOCATED(XWCV))THEN
  XWCV=XSPVAL
ENDIF
!
!*     2.        GETS VERTICAL CROSS-SECTION DATA THROUGH INTERPOLTATION
!                -------------------------------------------------------
! Prise en compte du 2D horizontal NON je prefere allouer correctement XWORKZ
!IF(IKU /= 1)THEN
CALL COMPCOORD_FORDIACHRO(NMGRID)
!ENDIF
IF(NVERBIA > 0)THEN
  print *,' ** PRECOU AP COMPCOORD_FORDIACHRO NMGRID ',NMGRID
  print *,' ** PRECOU Entree NPROFILE ',NPROFILE
ENDIF
print*, LUMVM,LDIRWIND,LDIRWM,LDIRWT

IF(LPRESY)THEN
  IF(NMGRID /= 1 .AND. SIZE(XPRES,1) /= 1 .AND. SIZE(XPRES,2) /= 1 .AND. &
     SIZE(XPRES,3) /= 1)THEN
    LPRESYT=.TRUE.
    print *,' ** PRECOU Appel volontaire INTERP_GRIDS NMGRID courant ',NMGRID,' IGRID de PR = 1 '
    CALL INTERP_GRIDS(0)
    LPRESYT=.FALSE.
  ENDIF
  XZZ(:,:,:)=XPRES(:,:,:,NLOOPT,1,1)
  print *,' ** PRECOU Remplacement volontaire de XZZ par XPRES(:,:,:,NLOOPT,1,1)'
! XZZ(:,:,:)=ALOG10(XZZ(:,:,:))
  IF(LPVT)THEN
    IF(.NOT.LTINCRDIA(NLOOPSUPER,1))THEN
      IF(NLOOPT == NTIMEDIA(1,NLOOPSUPER,1))THEN
        IF(ALLOCATED(XPRESM))THEN
          DEALLOCATE(XPRESM)
        ENDIF
        ALLOCATE(XPRESM(NBTIMEDIA(NLOOPSUPER,1),IKU))
        ITPRESY=0
      ELSE IF(NLOOPT == NTIMEDIA(NBTIMEDIA(NLOOPSUPER,1),NLOOPSUPER,1))THEN
      ENDIF
    ELSE
      IF(NLOOPT == NTIMEDIA(1,NLOOPSUPER,1))THEN
        IF(ALLOCATED(XPRESM))THEN
          DEALLOCATE(XPRESM)
        ENDIF
        IPRESM=(NTIMEDIA(2,NLOOPSUPER,1)-NTIMEDIA(1,NLOOPSUPER,1))/ &
        NTIMEDIA(3,NLOOPSUPER,1)+1
        ALLOCATE(XPRESM(IPRESM,IKU))
        ITPRESY=0
      ELSEIF(NLOOPT == NTIMEDIA(2,NLOOPSUPER,1))THEN
      ENDIF
    ENDIF
  ENDIF
ENDIF

!!!essai nov 2001
IF((LULM .OR. LULT .OR.LVTM .OR. LVTT) .AND. .NOT.(LCH .AND.LCV))THEN
!IF(LULM .OR. LULT .OR.LVTM .OR. LVTT)THEN
!!!essai nov 2001

  ZWORK3D=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
             NJINF-NJL+1:NJSUP-NJL+1, &
            :,NLOOPT,1,1)
  DO JKLOOP=1,IKU
    ZTEM1(:,:)=0.
    ZTEMV(:,:)=0.

    IF(JKLOOP <MAX(IKB,NKL) .OR. JKLOOP> MIN(NKH,IKE))THEN
    ELSE
      ZTEM1(NIL:NIH,NJL:NJH)=PWORK3D(:,:,JKLOOP-NKL+1)
      ZTEMV(NIL:NIH,NJL:NJH)=ZWORK3D(:,:,JKLOOP-NKL+1)
      CALL ROTA(ZTEM1,ZTEMV)

      IF(LULM .OR. LULT)THEN
        CALL COUPE_FORDIACHRO(ZTEM1,ZTEM2,JKLOOP)
      ELSE
        CALL COUPE_FORDIACHRO(ZTEMV,ZTEM2,JKLOOP)
      ENDIF

      PTEMCV(:,JKLOOP)=ZTEM2(:)
!     IF(LULM)THEN
!      print *,'LULM ZTEM2 JKLOOP ',JKLOOP
!      print *,ZTEM2
!     ENDIF
    ENDIF

  ENDDO

ELSE IF(LULMWM .OR. LULTWT)THEN

  NMGRID=1
  CALL COMPCOORD_FORDIACHRO(NMGRID)
! CALL COMPCOORD_FORDIACHRO(1)

  ZWORK3D=XV(NIINF-NIL+1:NISUP-NIL+1, &
             NJINF-NJL+1:NJSUP-NJL+1, &
            :,NLOOPT,1,1)
  ZWORK3W=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
             NJINF-NJL+1:NJSUP-NJL+1, &
            :,NLOOPT,1,1)
! On place la composante W aux points de masse
  ZWORK3W(:,:,1:IWKU-1)=.5*(ZWORK3W(:,:,1:IWKU-1)+ZWORK3W(:,:,2:IWKU))
  ZWORK3W(:,:,IWKU)=2.*ZWORK3W(:,:,IWKU-1)-ZWORK3W(:,:,IWKU-2)

  DO JKLOOP=1,IKU
    ZTEM1(:,:)=0.
    ZTEMV(:,:)=0.
    ZTEMW(:,:)=0.

    IF(JKLOOP <MAX(IKB,NKL) .OR. JKLOOP> MIN(NKH,IKE))THEN
    ELSE

      ZTEM1(NIL:NIH,NJL:NJH)=PWORK3D(:,:,JKLOOP-NKL+1)
      ZTEMV(NIL:NIH,NJL:NJH)=ZWORK3D(:,:,JKLOOP-NKL+1)
      ZTEMW(NIL:NIH,NJL:NJH)=ZWORK3W(:,:,JKLOOP-NKL+1)

      CALL COUPEUW_FORDIACHRO(ZTEM1,ZTEM2,JKLOOP,1)

!  Janvier 2001 ..PROVISOIRE
!     L2D=.FALSE.
      IF(L2D)THEN
! 2D // axe X
	ZTEMVR=ZTEMV(NIDEBCOU:NIDEBCOU+NLMAX-1,NJDEBCOU)
      ELSE
	CALL COUPEUW_FORDIACHRO(ZTEMV,ZTEMVR,JKLOOP,2)
      ENDIF

      CALL ROTAUW(ZTEM2,ZTEMVR)
      PTEMCV(:,JKLOOP)=ZTEM2

      CALL COUPEUW_FORDIACHRO(ZTEMW,ZTEMWR,JKLOOP,3)
      XWCV(:,JKLOOP)=ZTEMWR

    ENDIF
  ENDDO

! Janvier 2001 + LDIRWIND et LUMVM et LUTVT et LSUMVM et LSUTVT
!ELSE IF(LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT .OR. LDIRWIND)THEN
!! essai nov 2001
ELSE IF(LUMVM .OR. LUTVT .OR. LSUMVM .OR. LSUTVT .OR. &
!(LDIRWIND .AND. .NOT.(LCV .AND.LCH)))THEN
(LDIRWIND .AND. .NOT.(LCV .AND.LCH)) .OR. &
(LDIRWM .AND. .NOT.LDIRWIND)         .OR. &
(LDIRWT .AND. .NOT.LDIRWIND)          )THEN

  ZWORK3D=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
             NJINF-NJL+1:NJSUP-NJL+1, &
            :,NLOOPT,1,1)
! On positionne les 2 composantes aux points de masse

  IUI=SIZE(PWORK3D,1)
  IUJ=SIZE(PWORK3D,2)
print*, NGRIU,NGRIV,IKU,IUI,IUJ
!! Nov 2001 sauf si ce n'est deja fait
  IF(NGRIU == 1 .AND. NGRIV == 1)THEN
    print *,' ** Precou NGRIU=',NGRIU,' NGRIV=',NGRIV,' pas de repositionnement sur la grille de masse (deja fait) GRP=',CGROUP
  ELSE
!! Nov 2001 sauf si ce n'est deja fait
  PWORK3D(1:IUI-1,:,:)=0.5*(PWORK3D(2:IUI,:,:)+PWORK3D(1:IUI-1,:,:))
  PWORK3D(IUI,:,:)=2*PWORK3D(IUI-1,:,:)-PWORK3D(IUI-2,:,:)
  ZWORK3D(:,1:IUJ-1,:)=0.5*(ZWORK3D(:,2:IUJ,:)+ZWORK3D(:,1:IUJ-1,:))
  ZWORK3D(:,IUJ,:)=2*ZWORK3D(:,IUJ-1,:)-ZWORK3D(:,IUJ-2,:)
!! Nov 2001 sauf si ce n'est deja fait
  ENDIF
!! Nov 2001 sauf si ce n'est deja fait
  DO JKLOOP=1,IKU
    ZTEM1(:,:)=0.
    ZTEMV(:,:)=0.

    IF(JKLOOP <MAX(IKB,NKL) .OR. JKLOOP> MIN(NKH,IKE))THEN
    ELSE

      ZTEM1(NIL:NIH,NJL:NJH)=PWORK3D(:,:,JKLOOP-NKL+1)
      ZTEMV(NIL:NIH,NJL:NJH)=ZWORK3D(:,:,JKLOOP-NKL+1)
	if(nverbia > 5)then
	  print*,'** PRECOU Composante U av coupe'
	endif

      CALL COUPE_FORDIACHRO(ZTEM1,ZTEM2,JKLOOP)
!     CALL COUPEUW_FORDIACHRO(ZTEM1,ZTEM2,JKLOOP,1)
      PTEMCV(:,JKLOOP)=ZTEM2
	if(nverbia > 0)then
	  print *,' ** PRECOU Composante U ap coupe, K= ',JKLOOP
	endif

!  Janvier 2001 ..PROVISOIRE
!     L2D=.FALSE.
      IF(L2D)THEN
! 2D // axe X
	ZTEMVR=ZTEMV(NIDEBCOU:NIDEBCOU+NLMAX-1,NJDEBCOU)
      ELSE
	if(nverbia > 5)then
	  print *,' ** PRECOU Composante V AV coupe'
	endif
	CALL COUPE_FORDIACHRO(ZTEMV,ZTEMVR,JKLOOP)
!       CALL COUPEUW_FORDIACHRO(ZTEMV,ZTEMVR,JKLOOP,2)
	if(nverbia > 0)then
	  print *,' ** PRECOU Composante V ap coupe, K= ',JKLOOP
	endif
      ENDIF

      XWCV(:,JKLOOP)=ZTEMVR
    ENDIF
  ENDDO
!! 30 nov 2001
!     IF(LDIRWIND)THEN
     IF(LDIRWIND .OR.  &
        (LDIRWM .AND. .NOT.LDIRWIND)     .OR. &
        (LDIRWT .AND. .NOT.LDIRWIND)          ) THEN
      IUB1=SIZE(XWCV,1)
      IUB2=SIZE(XWCV,2)
      ISKIP=1
      ITER=IUB1; JTER=IUB2
      IF(ALLOCATED(ZX))THEN
        DEALLOCATE(ZX)
      ENDIF
      IF(ALLOCATED(ZZY))THEN
        DEALLOCATE(ZZY)
      ENDIF
      ALLOCATE(ZX(ITER,1),ZZY(JTER))
      ZX(:,1)=XZZX(1:IUB1:ISKIP)
      ZZY=XZZY(1:IUB2:ISKIP)
!! DEc 2001
!!Fev 2002
      IF(LDIRWIND .AND. (LCH .OR. LFT .OR. LPVKT ))THEN
!     IF(LCH .OR. LFT .OR. LPVKT)THEN
!!Fev 2002
!! DEc 2001
      CALL COMPUTEDIR(ITER,JTER,IUB1,IUB2,ISKIP,PTEMCV,XWCV)
      PTEMCV(:,:)=XWCV(:,:)
!! DEc 2001
      ENDIF
!! DEc 2001
      IF ( (LDIRWM .AND. .NOT.LDIRWIND)     .OR. &
           (LDIRWT .AND. .NOT.LDIRWIND)          )THEN
     print*,'precou av dd ',MINVAL(PTEMCV),MAXVAL(PTEMCV),MINVAL(XWCV),MAXVAL(XWCV)
      CALL COMPUTEDIR(ITER,JTER,IUB1,IUB2,ISKIP,PTEMCV,XWCV)
      PTEMCV(:,:)=XWCV(:,:)
     print*,'precou ap dd ',MINVAL(PTEMCV),MAXVAL(PTEMCV)
      ENDIF
     ENDIF
!! 30 nov 2001

!!essai Nov 2001 -> PH traites ds traceh_fordiachro
ELSE IF((LMUMVM .OR. LMUTVT .OR. LMLSUMVM .OR. LMLSUTVT) .AND. &
        (.NOT.(LCH.AND.LCV)))THEN
!ELSE IF(LMUMVM .OR. LMUTVT .OR. LMLSUMVM .OR. LMLSUTVT)THEN
!!essai Nov 2001


  CALL COMPCOORD_FORDIACHRO(NMGRID)
  ZWORK3D=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
             NJINF-NJL+1:NJSUP-NJL+1, &
            :,NLOOPT,1,1)

! On positionne les 2 composantes aux points de masse

  if(nverbia > 0 .AND. size(PWORK3D,1) >= 12  .AND. &
  size(PWORK3D,2) >= 7 .AND. size(PWORK3D,3) >= 9)THEN
  print *,' ***PRECOU IK=9 I=8A12 J=3A7 U Grid 2 et V GRID 3 '
  print *,PWORK3D(8:12,3,9)
  print *,PWORK3D(8:12,4,9)
  print *,PWORK3D(8:12,5,9)
  print *,PWORK3D(8:12,6,9)
  print *,PWORK3D(8:12,7,9),' *******'
  print *,ZWORK3D(8:12,3,9)
  print *,ZWORK3D(8:12,4,9)
  print *,ZWORK3D(8:12,5,9)
  print *,ZWORK3D(8:12,6,9)
  print *,ZWORK3D(8:12,7,9),' *******'
  endif
  IUI=SIZE(PWORK3D,1)
  IUJ=SIZE(PWORK3D,2)
!! Nov 2001 sauf si ce n'est deja fait
  IF(NGRIU == 1 .AND. NGRIV == 1)THEN
    print *,' ** Precou NGRIU=',NGRIU,' NGRIV=',NGRIV,' pas de repositionnement sur la grille de masse (deja fait) GRP=',CGROUP
  ELSE
!! Nov 2001 sauf si ce n'est deja fait
  PWORK3D(1:IUI-1,:,:)=0.5*(PWORK3D(2:IUI,:,:)+PWORK3D(1:IUI-1,:,:))
  PWORK3D(IUI,:,:)=2*PWORK3D(IUI-1,:,:)-PWORK3D(IUI-2,:,:)
  ZWORK3D(:,1:IUJ-1,:)=0.5*(ZWORK3D(:,2:IUJ,:)+ZWORK3D(:,1:IUJ-1,:))
  ZWORK3D(:,IUJ,:)=2*ZWORK3D(:,IUJ-1,:)-ZWORK3D(:,IUJ-2,:)
!! Nov 2001 sauf si ce n'est deja fait
  ENDIF
!! Nov 2001 sauf si ce n'est deja fait
  if(nverbia > 0 .AND. size(PWORK3D,1) >= 12 .AND. &
   size(PWORK3D,2) >= 7 .AND. size(PWORK3D,3) >= 9)THEN
  print *,' ***PRECOU IK=9 I=8A12 J=3A7 U et V Grille 1 '
  print *,PWORK3D(8:12,3,9)
  print *,PWORK3D(8:12,4,9)
  print *,PWORK3D(8:12,5,9)
  print *,PWORK3D(8:12,6,9)
  print *,PWORK3D(8:12,7,9),' *******'
  print *,ZWORK3D(8:12,3,9)
  print *,ZWORK3D(8:12,4,9)
  print *,ZWORK3D(8:12,5,9)
  print *,ZWORK3D(8:12,6,9)
  print *,ZWORK3D(8:12,7,9),' *******'
  endif
  PWORK3D=PWORK3D*PWORK3D
  ZWORK3D=ZWORK3D*ZWORK3D
  PWORK3D=SQRT(PWORK3D+ZWORK3D)

  DO JKLOOP=1,IKU
    ZTEM1(:,:)=0.
    IF(JKLOOP <MAX(IKB,NKL) .OR. JKLOOP> MIN(NKH,IKE))THEN
    ELSE
      ZTEM1(NIL:NIH,NJL:NJH)=PWORK3D(:,:,JKLOOP-NKL+1)
  !   ZTEM1(:,:)=PWORK3D(:,:,JKLOOP-NKL+1)
    ENDIF
    CALL COUPE_FORDIACHRO(ZTEM1,ZTEM2,JKLOOP)
    PTEMCV(:,JKLOOP)=ZTEM2(:)
  
  !print *,' JKLOOP NKL NKH ',JKLOOP,NKL,NKH,'   ZTEM2'
  !print *,ZTEM2
  ENDDO

ELSE
IF(NVERBIA > 0)THEN
  print *,' ** PRECOU AV DO JKLOOP=1,IKU'
ENDIF

DO JKLOOP=1,IKU
    ZTEM1(:,:)=0.
! Ajout Avril 2001

!!Nov 2001
  IF(IKU == 1 )THEN
! IF(IKU == 1 .AND. LKCP)THEN
!!Nov 2001
    
    ZTEM1(NIL:NIH,NJL:NJH)=PWORK3D(:,:,1)
    IF(NVERBIA > 5)THEN
       print *,' ** PRECOU LKCP=',LKCP,' IKU=',IKU,' ZTEM1(NIL:NIH,NJL:NJH)'
       print *,ZTEM1(NIL:NIH,NJL:NJH)
    ENDIF

  ELSE

  IF(JKLOOP <MAX(IKB,NKL) .OR. JKLOOP> MIN(NKH,IKE))THEN
  ELSE
    ZTEM1(NIL:NIH,NJL:NJH)=PWORK3D(:,:,JKLOOP-NKL+1)
!   ZTEM1(:,:)=PWORK3D(:,:,JKLOOP-NKL+1)
  ENDIF
IF(NVERBIA > 5)THEN
  IF(JKLOOP == MAX(2,NKL) .OR. IKU == 1)THEN
  print *,' ** PRECOU DS DO JKLOOP=1,IKU  AV COUPE, JKLOOP',JKLOOP
  print *,' ** PRECOU  AV COUPE, ZTEM2 ',ZTEM2
  ENDIF
ENDIF

  ENDIF

  CALL COUPE_FORDIACHRO(ZTEM1,ZTEM2,JKLOOP)
  PTEMCV(:,JKLOOP)=ZTEM2(:)

IF(NVERBIA > 5)THEN
  IF(JKLOOP == MAX(2,NKL) .OR. IKU == 1)THEN
print *,' JKLOOP NKL NKH ',JKLOOP,NKL,NKH,'   ZTEM2'
print *,ZTEM2
  ENDIF
ENDIF
ENDDO
IF(NVERBIA > 0)THEN
  print *,' **Sortie PRECOU (XWORKZ) ',SIZE(XWORKZ,1),SIZE(XWORKZ,2),&
  SIZE(XWORKZ,3)
! print *,' **Sortie PRECOU  XWORKZ(NPROFILE,:,NMGRID) ',XWORKZ(NPROFILE,:,NMGRID)
ENDIF
IF(LPRESY .AND. LPVT)THEN
  ITPRESY=ITPRESY+1
  XPRESM(ITPRESY,:)=XWORKZ(NPROFILE,:,NMGRID)
ENDIF

ENDIF

!print *,' ** precou AV DEALLOCATE(ZTEM1,ZTEM2) '
DEALLOCATE(ZTEM1,ZTEM2)
!print *,' ** precou AP DEALLOCATE(ZTEM1,ZTEM2) '
IF(ALLOCATED(ZTEMWR))THEN
  DEALLOCATE(ZTEMWR)
ENDIF
IF(ALLOCATED(ZTEMVR))THEN
  DEALLOCATE(ZTEMVR)
ENDIF
IF(ALLOCATED(ZTEMW))THEN
  DEALLOCATE(ZTEMW)
ENDIF
IF(ALLOCATED(ZWORK3W))THEN
  DEALLOCATE(ZWORK3W)
ENDIF
IF(ALLOCATED(ZTEMV))THEN
  DEALLOCATE(ZTEMV)
ENDIF
IF(ALLOCATED(ZWORK3D))THEN
  DEALLOCATE(ZWORK3D)
ENDIF
if(nverbia > 0)then
 print *,' ** precou FIN'
endif
!
!----------------------------------------------------------------------------
!
!*       3.      EXIT
!                ----
!
END SUBROUTINE  PRECOU_FORDIACHRO
