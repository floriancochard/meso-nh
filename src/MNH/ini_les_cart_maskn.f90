!MNH_LIC Copyright 1994-2017 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!#############################
MODULE MODI_INI_LES_CART_MASKn
!#############################
!
!
!
INTERFACE
!
      SUBROUTINE  INI_LES_CART_MASK_n(KMI, PXHAT_ll,PYHAT_ll,                &
                                     KLES_IINF,KLES_JINF,KLES_ISUP,KLES_JSUP)
!

INTEGER,            INTENT(IN)  :: KMI      ! son model index
REAL, DIMENSION(:), INTENT(IN)  :: PXHAT_ll ! son model X coordinate
REAL, DIMENSION(:), INTENT(IN)  :: PYHAT_ll ! son model X coordinate
INTEGER,            INTENT(OUT) :: KLES_IINF ! limits of the cartesian
INTEGER,            INTENT(OUT) :: KLES_JINF ! mask in son model
INTEGER,            INTENT(OUT) :: KLES_ISUP ! domain
INTEGER,            INTENT(OUT) :: KLES_JSUP !
!
END SUBROUTINE INI_LES_CART_MASK_n
!
END INTERFACE
!
END MODULE MODI_INI_LES_CART_MASKn
!
!     ########################################################################
      SUBROUTINE  INI_LES_CART_MASK_n(KMI, PXHAT_ll,PYHAT_ll,                &
                                      KLES_IINF,KLES_JINF,KLES_ISUP,KLES_JSUP)
!     ########################################################################
!
!
!!****  *INI_LES_CART_MASK_n* initializes the LES cartesian mask
!!
!!    PURPOSE
!!    -------
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
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!      Modification     01/02/01 (D.Gazen) add module MODD_NSV for NSV variable
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      P.Wautelet: 19/10/2017 : IO: removed extern_userio.f90
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODE_ll
USE MODE_GATHER_ll 
USE MODE_MODELN_HANDLER
!
USE MODD_CONF
USE MODD_PARAMETERS
!
USE MODD_GRID_n
USE MODD_LES
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
INTEGER,            INTENT(IN)  :: KMI      ! son model index
REAL, DIMENSION(:), INTENT(IN)  :: PXHAT_ll ! son model X coordinate
REAL, DIMENSION(:), INTENT(IN)  :: PYHAT_ll ! son model X coordinate
INTEGER,            INTENT(OUT) :: KLES_IINF ! limits of the cartesian
INTEGER,            INTENT(OUT) :: KLES_JINF ! mask in son model
INTEGER,            INTENT(OUT) :: KLES_ISUP ! domain
INTEGER,            INTENT(OUT) :: KLES_JSUP !
!
!
!       0.2  declaration of local variables
!
!
INTEGER :: IIMAX_ll  ! total physical father domain I size
INTEGER :: IJMAX_ll  ! total physical father domain J size
!
INTEGER :: IIB_ll    ! son domain index
INTEGER :: IIE_ll    ! son domain index
INTEGER :: IJB_ll    ! son domain index
INTEGER :: IJE_ll    ! son domain index
!
INTEGER :: JI, JJ    ! loop counters
!
REAL    :: ZX, ZY    ! coordinates of msak boundaries
!
INTEGER :: IINFO_ll, IRESP
!
REAL, DIMENSION(:), ALLOCATABLE :: ZXHAT_ll ! father model coordinates
REAL, DIMENSION(:), ALLOCATABLE :: ZYHAT_ll !
INTEGER :: IMI
!
IMI = GET_CURRENT_MODEL_INDEX()
!-------------------------------------------------------------------------------
!
!*      1.   Coordinate of father model
!            --------------------------
!
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll) 
ALLOCATE(ZXHAT_ll(IIMAX_ll+ 2 * JPHEXT))
ALLOCATE(ZYHAT_ll(IJMAX_ll+ 2 * JPHEXT))
CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP)
CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP)
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
!-------------------------------------------------------------------------------
!
IIB_ll=1+JPHEXT
IIE_ll=SIZE(PXHAT_ll)-JPHEXT
IJB_ll=1+JPHEXT
IJE_ll=SIZE(PYHAT_ll)-JPHEXT
!
!-------------------------------------------------------------------------------
!
!*      2.   X limit of LES cartesian mask
!            -----------------------------
!
!* left limit
!
ZX = ZXHAT_ll(NLESn_IINF(IMI))
IF (PXHAT_ll(IIB_ll)>ZX) THEN
  KLES_IINF=IIB_ll ! father mask starts left of son domain
ELSE IF (PXHAT_ll(IIE_ll+1)<ZX) THEN
  CALL MASK_OVER_ALL_DOMAIN
  RETURN
ELSE
  DO JI=IIB_ll,IIE_ll
    IF (ABS(PXHAT_ll(JI)-ZX) <= (PXHAT_ll(JI+1)-PXHAT_ll(JI))/2. ) THEN
      KLES_IINF=JI
    END IF
  END DO
END IF
!
!* right limit
!
ZX = ZXHAT_ll(NLESn_ISUP(IMI)+1)
IF (PXHAT_ll(IIB_ll)>ZX) THEN
  CALL MASK_OVER_ALL_DOMAIN
  RETURN
ELSE IF (PXHAT_ll(IIE_ll+1)<ZX) THEN
  KLES_ISUP=IIE_ll ! father mask ends right of son domain
ELSE
  DO JI=IIB_ll,IIE_ll
    IF (ABS(PXHAT_ll(JI+1)-ZX) <= (PXHAT_ll(JI+1)-PXHAT_ll(JI))/2. ) THEN
      KLES_ISUP=JI
    END IF
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*      3.   Y limit of LES cartesian mask
!            -----------------------------
!
!* bottom limit
!
ZY = ZYHAT_ll(NLESn_JINF(IMI))
IF (PYHAT_ll(IJB_ll)>ZY) THEN
  KLES_JINF=IJB_ll ! father mask starts under the son domain
ELSE IF (PYHAT_ll(IJE_ll+1)<ZY) THEN
  CALL MASK_OVER_ALL_DOMAIN
  RETURN
ELSE
  DO JJ=IJB_ll,IJE_ll
    IF (ABS(PYHAT_ll(JJ)-ZY) <= (PYHAT_ll(JJ+1)-PYHAT_ll(JJ))/2. ) THEN
      KLES_JINF=JJ
    END IF
  END DO
END IF
!
!* top limit
!
ZY = ZYHAT_ll(NLESn_JSUP(IMI)+1)
IF (PYHAT_ll(IJB_ll)>ZY) THEN
  CALL MASK_OVER_ALL_DOMAIN
  RETURN
ELSE IF (PYHAT_ll(IJE_ll+1)<ZY) THEN
  KLES_JSUP=IJE_ll ! father mask ends over the son domain
ELSE
  DO JJ=IJB_ll,IJE_ll
    IF (ABS(PYHAT_ll(JJ+1)-ZY) <= (PYHAT_ll(JJ+1)-PYHAT_ll(JJ))/2. ) THEN
      KLES_JSUP=JJ
    END IF
  END DO
END IF
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZXHAT_ll)
DEALLOCATE(ZYHAT_ll)
!-------------------------------------------------------------------------------
!
  CONTAINS
!
  SUBROUTINE MASK_OVER_ALL_DOMAIN
    KLES_IINF=IIB_ll ! father mask not in son domain, so all domain is taken
    KLES_ISUP=IIE_ll
    KLES_JINF=IJB_ll
    KLES_JSUP=IJE_ll
    DEALLOCATE(ZXHAT_ll)
    DEALLOCATE(ZYHAT_ll)
  END SUBROUTINE MASK_OVER_ALL_DOMAIN
!
END SUBROUTINE INI_LES_CART_MASK_n   

