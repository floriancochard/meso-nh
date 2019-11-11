!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/post/s.rotauw.f90, Version:1.2, Date:98/10/02, Last modified:98/06/04
!-----------------------------------------------------------------
!     ######spl
      SUBROUTINE ROTAUW(PTEM1,PTEMV)
!     ##############################
!
!!****  *ROTAUW* - For the vertical oblique cross-sections, rotates the wind
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
!
IMPLICIT NONE
!
!*       0.1  Dummy arguments and results
!
                                              ! On entry, model x-y components
                                              ! of the wind. 1 stands for U,
REAL, DIMENSION(:),  INTENT(INOUT) :: PTEM1   ! V stands for V. On return, 
REAL, DIMENSION(:),  INTENT(INOUT) :: PTEMV   ! longitudinal, transverse
                                              ! wind components with respect
                                              ! to the current olblique
                                              ! vertical section plane.
!
!*       0.2  Local variables
!
INTEGER             :: IWIU
INTEGER             :: J
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
!
!*       1.2   Rotation to the natural frame of the oblique section
!
IF(LDEFCV2CC)THEN
  ZRANGLE=XANGLECV
ELSE
ZRANGLE=FLOAT(NLANGLE)*ACOS(-1.)/180.   ! NLANGLE is the section direction
ENDIF
ZCANGLE=COS(ZRANGLE)
ZSANGLE=SIN(ZRANGLE)
IF(.NOT.LDEFCV2CC)THEN
  IF(NLANGLE.EQ.0.OR.NLANGLE.EQ.180)ZSANGLE=0.
  IF(NLANGLE.EQ.90.OR.NLANGLE.EQ.270)ZCANGLE=0.
ELSE
  IF(XANGLECV == 0. .OR. XANGLECV/ACOS(-1.)*180. == 180.)ZSANGLE=0.
  IF(XANGLECV/ACOS(-1.)*180. == 90. .OR.XANGLECV/ACOS(-1.)*180. == 270.)ZCANGLE=0.
ENDIF
DO J=1,IWIU
ZU=PTEM1(J)
ZV=PTEMV(J)
PTEM1(J)=ZU*ZCANGLE+ZV*ZSANGLE
PTEMV(J)=-ZU*ZSANGLE+ZV*ZCANGLE
ENDDO
!
!------------------------------------------------------------------------------
!
!*        2.     EXIT
!                ----
!
RETURN
END SUBROUTINE ROTAUW
