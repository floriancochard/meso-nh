!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_real 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_WINDCALC
!     ####################
INTERFACE
!
      SUBROUTINE WINDCALC(PUIN,PVIN,KXCEN,KYCEN,PUOUT,PVOUT)
!
REAL, DIMENSION(:,:,:),  INTENT(IN)   :: PUIN  ! U or tangential wind component
REAL, DIMENSION(:,:,:),  INTENT(IN)   :: PVIN  ! V or radial wind component
INTEGER, DIMENSION(:), INTENT(IN)     :: KXCEN
INTEGER, DIMENSION(:), INTENT(IN)     :: KYCEN
REAL, DIMENSION(:,:,:),  INTENT(OUT)  :: PUOUT ! U or tangential component
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL :: PVOUT ! V or radial component
!
END SUBROUTINE WINDCALC
!
END INTERFACE
!
END MODULE MODI_WINDCALC
!
!
!
!     ######################################################
      SUBROUTINE WINDCALC(PUIN,PVIN,KXCEN,KYCEN,PUOUT,PVOUT)
!     ######################################################
!
!!****  *WINDCALC* - interpolation between cylindrical grid and cartesian grid
!!
!!    PURPOSE
!!    -------
!       If the input arrays are in cylindrical grid (tangential and radial wind
!     components), they are interpolated into cartesian grid
!     (model grid components).
!       If the input arrays are in cartesian grid (model frid components),
!     they are interpolated into cylindrical grid (tangential and radial wind
!     components).
!       The conversion formula are the same..
!       (xcen,ycen) is the vortex center at level pk
!
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      R. Rogers             * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY: XPI
USE MODD_GRID_n, ONLY: XXHAT,XYHAT
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:),  INTENT(IN)   :: PUIN  ! U or tangential wind component
REAL, DIMENSION(:,:,:),  INTENT(IN)   :: PVIN  ! V or radial wind component
INTEGER, DIMENSION(:), INTENT(IN)     :: KXCEN
INTEGER, DIMENSION(:), INTENT(IN)     :: KYCEN
REAL, DIMENSION(:,:,:),  INTENT(OUT)  :: PUOUT ! U or tangential component
REAL, DIMENSION(:,:,:),  INTENT(OUT), OPTIONAL :: PVOUT ! V or radial component
!
!
!*       0.2   Declarations of local variables
!
INTEGER          :: JI,JJ,JP     ! Loop indexes along the x,y,z directions
INTEGER          :: IX,IY,IP     ! dimensions along the x,y,z directions
REAL             :: ZDELTAX,ZDELTAY,ZDELTAR,ZXK,ZYK
REAL             :: ZR2,ZR,ZPHI
!
!-------------------------------------------------------------------------------
!
!*	 1.     COMPUTATION OF THE FIELD 
!               ------------------------
!
IX=SIZE(PUIN,1) 
IY=SIZE(PUIN,2) 
IP=SIZE(PUIN,3)
!
ZDELTAX = XXHAT(3) - XXHAT(2)
ZDELTAY = XYHAT(3) - XYHAT(2)
ZDELTAR = MAX(ZDELTAX,ZDELTAY)
!
DO JP = 1, IP
  DO JJ = 1, IY
    DO JI = 1, IX
      ZXK = XXHAT(JI) - XXHAT(KXCEN(JP))
      ZYK = XYHAT(JJ) - XYHAT(KYCEN(JP))
      ZR2 = ZXK * ZXK + ZYK * ZYK
      ZR = SQRT(ZR2)
      IF (ZR == 0. .OR. ZR < ZDELTAR) THEN
        PUOUT(JI,JJ,JP) = 0.
        IF (PRESENT(PVOUT)) PVOUT(JI,JJ,JP) =  0.
      ELSE
        ZPHI = ATAN2(ZYK,ZXK)
        IF (ZPHI < 0.) ZPHI = ZPHI + 2*XPI
        PUOUT(JI,JJ,JP) = PVIN(JI,JJ,JP) * COS(ZPHI) -   &
                          PUIN(JI,JJ,JP) * SIN(ZPHI)
        IF (PRESENT(PVOUT)) PVOUT(JI,JJ,JP) =  PUIN(JI,JJ,JP) * COS(ZPHI) + &
                                               PVIN(JI,JJ,JP) * SIN(ZPHI)
      END IF
    END DO
  END DO
END DO
! 
!-------------------------------------------------------------------------------
END SUBROUTINE WINDCALC
