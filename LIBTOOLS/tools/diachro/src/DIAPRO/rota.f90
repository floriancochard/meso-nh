!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.rota.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE ROTA(PTEM1,PTEMV)
!     ############################
!
!!****  *ROTA* - For the vertical oblique cross-sections, rotates the wind
!!****           components from the model frame to the section natural frame
!!
!!    PURPOSE
!!    -------
!       In the case of oblique vertical cross-sections, computes the
!       longitudinal and transverse components of the wind with respect
!       to the section plane.
!
!!**  METHOD
!!    ------
!!    To make a physically meanigfull rotation, the u and v components
!!    of the wind are interpolated back to be colocated at the mass gridpoint.
!!
!!    EXTERNAL
!!    --------
!!      COS  ! trigonometric functions
!!      SIN  !
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist (former PARA common)
!!         NLANGLE :  Angle between X Meso-NH axis and
!!                    cross-section direction in degrees
!!                    (Integer value anticlockwise)
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   13/01/95
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODN_PARA
USE MODD_DEFCV
USE MODD_MEMGRIUV
USE MODD_RESOLVCAR
!
IMPLICIT NONE
!
!*       0.1  Dummy arguments and results
!
                                              ! On entry, model x-y components
                                              ! of the wind. 1 stands for U,
REAL, DIMENSION(:,:),  INTENT(INOUT) :: PTEM1 ! V stands for V. On return, 
REAL, DIMENSION(:,:),  INTENT(INOUT) :: PTEMV ! longitudinal, transverse
                                              ! wind components with respect
                                              ! to the current olblique
                                              ! vertical section plane.
!
!*       0.2  Local variables
!
INTEGER             :: IWIU, IWJU
INTEGER             :: J, JA
!
REAL                :: ZU, ZV
REAL                :: ZRANGLE, ZCANGLE, ZSANGLE
!
!-------------------------------------------------------------------------------
!
!*        1.   COMPUTING THE LONGITUDINAL AND TRANSVERSE COMPONENTS
!              ----------------------------------------------------
!
!*        1.1  Array sizes calculations
!
IWIU=SIZE(PTEM1,1)
IWJU=SIZE(PTEM1,2)
!
!*        1.2  Wind component are interpolated back to the mass point:
!*             only colocated u-v components can be mixed in a rotation
!*             in a physically meaningfull way to obtain lonitudinal and
!*             transverse components
!
!! Nov 2001 sauf si ce n'est deja fait
IF(NGRIU == 1 .AND. NGRIV == 1)THEN
    print *,' ** Rota  NGRIU=',NGRIU,' NGRIV=',NGRIV,' pas de repositionnement sur la grille de masse (deja fait) GRP=',CGROUP
ELSE
!! Nov 2001 sauf si ce n'est deja fait
PTEM1(1:IWIU-1,:)=.5*(PTEM1(1:IWIU-1,:)+PTEM1(2:IWIU,:))
PTEM1(IWIU,:)=2.*PTEM1(IWIU-1,:)-PTEM1(IWIU-2,:)
PTEMV(:,1:IWJU-1)=.5*(PTEMV(:,1:IWJU-1)+PTEMV(:,2:IWJU))
PTEMV(:,IWJU)=2.*PTEMV(:,IWJU-1)-PTEMV(:,IWJU-2)
!! Nov 2001 sauf si ce n'est deja fait
ENDIF
!! Nov 2001 sauf si ce n'est deja fait
!
!*       1.3   Rotation to the natural frame of the oblique section
!
!!! Essai Nov 2001 pour prise en compte PH A suivre... 29/11/2001
IF(((LCH.AND.LULM).OR.(LCH.AND.LULT).OR.(LCH.AND.LVTM).OR. &
   (LCH.AND.LVTT)) .AND. .NOT.LCV)THEN
!IF((LCH.AND.LULM).OR.(LCH.AND.LULT).OR.(LCH.AND.LVTM).OR. &
!   (LCH.AND.LVTT))THEN
!!! Essai Nov 2001 pour prise en compte PH A suivre... 29/11/2001
  ZRANGLE=XANGULVT*ACOS(-1.)/180.
ELSE
IF(LDEFCV2CC)THEN
  ZRANGLE=XANGLECV
ELSE
  ZRANGLE=FLOAT(NLANGLE)*ACOS(-1.)/180.   ! NLANGLE is the section direction
ENDIF
ENDIF
ZCANGLE=COS(ZRANGLE)
ZSANGLE=SIN(ZRANGLE)
!!! Essai Nov 2001 pour prise en compte PH A suivre... 29/11/2001
IF(((LCH.AND.LULM).OR.(LCH.AND.LULT).OR.(LCH.AND.LVTM).OR. &
   (LCH.AND.LVTT)) .AND. .NOT.LCV)THEN
!IF((LCH.AND.LULM).OR.(LCH.AND.LULT).OR.(LCH.AND.LVTM).OR. &
!   (LCH.AND.LVTT))THEN
!!! Essai Nov 2001 pour prise en compte PH A suivre... 29/11/2001
  IF(XANGULVT == 0. .OR. XANGULVT == 180.)ZSANGLE=0.
  IF(XANGULVT == 90. .OR. XANGULVT == 270.)ZCANGLE=0.
ELSE
IF(.NOT.LDEFCV2CC)THEN
  IF(NLANGLE.EQ.0.OR.NLANGLE.EQ.180)ZSANGLE=0.
  IF(NLANGLE.EQ.90.OR.NLANGLE.EQ.270)ZCANGLE=0.
ELSE
  IF(XANGLECV == 0. .OR. XANGLECV/ACOS(-1.)*180. == 180.)ZSANGLE=0.
  IF(XANGLECV/ACOS(-1.)*180. == 90. .OR.XANGLECV/ACOS(-1.)*180. == 270.)ZCANGLE=0.
ENDIF
ENDIF
IF(nverbia > 0)THEN
  print *,' ** rota XANGULVT,ZSANGLE,ZCANGLE ',XANGULVT,ZSANGLE,ZCANGLE
endif
DO J=1,IWIU
DO JA=1,IWJU
ZU=PTEM1(J,JA)
ZV=PTEMV(J,JA)
PTEM1(J,JA)=ZU*ZCANGLE+ZV*ZSANGLE
PTEMV(J,JA)=-ZU*ZSANGLE+ZV*ZCANGLE
ENDDO
ENDDO
!
!*       1.4   Rotated components re-interpolated back to their nominal
!*             Meso-NH locations
!
! Suppression debut Avril 99 a la demande de Joel, Nicole et les autres..
!PTEM1(2:IWIU,:)=.5*(PTEM1(1:IWIU-1,:)+PTEM1(2:IWIU,:))
!PTEM1(1,:)=2.*PTEM1(2,:)-PTEM1(3,:)
!PTEMV(:,2:IWJU)=.5*(PTEMV(:,1:IWJU-1)+PTEMV(:,2:IWJU))
!PTEMV(:,1)=2.*PTEMV(:,2)-PTEMV(:,3)
!
!------------------------------------------------------------------------------
!
!*        2.     EXIT
!                ----
!
RETURN
END SUBROUTINE ROTA
