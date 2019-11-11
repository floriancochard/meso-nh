!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ###########################
      MODULE MODI_ADVECUVW_WENO_K
!     ###########################
!
INTERFACE
!
      SUBROUTINE ADVECUVW_WENO_K(HLBCX, HLBCY, KWENO_ORDER, PUT, PVT, PWT,     &
                             PRUCT, PRVCT, PRWCT, PRUS, PRVS, PRWS, TPHALO2LIST)
!
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
INTEGER,                  INTENT(IN)    :: KWENO_ORDER   ! Order of the WENO
                                                         ! scheme (3 or 5)
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVCT !  components
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRWCT ! of momentum
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUT, PVT, PWT        ! U,V,W at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS     ! Source terms
!
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list for diffusion
!
END SUBROUTINE ADVECUVW_WENO_K
!
END INTERFACE
!
END MODULE MODI_ADVECUVW_WENO_K
!
!     ##########################################################################
      SUBROUTINE ADVECUVW_WENO_K(HLBCX, HLBCY, KWENO_ORDER, PUT, PVT, PWT,     &
                             PRUCT, PRVCT, PRWCT, PRUS, PRVS, PRWS, TPHALO2LIST)
!     ##########################################################################
!
!!    AUTHOR
!!    ------
!!
!!
!!    MODIFICATIONS
!!    ------------- 
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 tests
!!		  T.Lunet	 	02/10/2014: add get_halo for WENO 5
!!					suppress comment of NHALO=1 tests
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
!
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
USE MODI_SHUMAN
USE MODI_ADVEC_WENO_K_1_AUX
USE MODI_ADVEC_WENO_K_2_AUX
USE MODI_ADVEC_WENO_K_3_AUX
!
USE MODD_CONF,        ONLY : NHALO
USE MODE_MPPDB
USE MODI_GET_HALO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X direction LBC type
CHARACTER (LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y direction LBC type
INTEGER,                  INTENT(IN)    :: KWENO_ORDER   ! Order of the WENO
                                                         ! scheme (3 or 5)
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRUCT  ! contravariant
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRVCT  !  components
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRWCT  ! of momentum
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUT, PVT, PWT     ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS     ! Source terms
!
TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST ! list for diffusion
!
!*       0.2   Declarations of local variables :
!
TYPE(HALO2LIST_ll), POINTER :: TZHALO2_UT,TZHALO2_VT,TZHALO2_WT

TYPE(LIST_ll), POINTER :: TZHALO2_ZMEAN
INTEGER                     :: IINFO_ll    ! return code of parallel routine
!
REAL, DIMENSION(SIZE(PUT,1), SIZE(PUT,2), SIZE(PUT,3)) :: ZMEAN, ZWORK
!
INTEGER :: K_SCHEME
INTEGER :: IKU
INTEGER :: IWORK
!
!------------------------- ADVECTION OF MOMENTUM ------------------------------
!
!
TZHALO2_UT => TPHALO2LIST                   ! 1rst add3dfield in model_n
TZHALO2_VT => TPHALO2LIST%NEXT              ! 2nd  add3dfield in model_n
TZHALO2_WT => TPHALO2LIST%NEXT%NEXT         ! 3rst add3dfield in model_n
!
IKU=SIZE(PUT,3)
!      -------------------------------------------------------
!
SELECT CASE(KWENO_ORDER)
!
CASE(1) ! WENO 1
!
!  U component
!
  PRUS = PRUS - DXM(UP_UX(PUT,MXF(PRUCT)))
!
  PRUS = PRUS - DYF(UP_MY(PUT,MXM(PRVCT)))
!
  PRUS = PRUS - DZF(1,IKU,1,UP_MZ(PUT,MXM(PRWCT)))
!
! V component
!
  PRVS = PRVS - DXF(UP_MX(PVT,MYM(PRUCT)))
!
  PRVS = PRVS - DYM(UP_VY(PVT,MYF(PRVCT)))
!
  PRVS = PRVS - DZF(1,IKU,1,UP_MZ(PVT,MYM(PRWCT)))
!
! W component
!
  PRWS = PRWS - DXF(UP_MX(PWT,MZM(1,IKU,1,PRUCT)))
!
  PRWS = PRWS - DYF(UP_MY(PWT,MZM(1,IKU,1,PRVCT)))
!
  PRWS = PRWS - DZM(1,IKU,1,UP_WZ(PWT,MZF(1,IKU,1,PRWCT)))
!
!
CASE(3) ! WENO 3
!
! U component
!
  ZWORK = MXF(PRUCT)
  CALL ADVEC_WENO_K_2_UX(HLBCX, PUT, ZWORK, ZMEAN, TZHALO2_UT%HALO2)
  PRUS = PRUS - DXM(ZMEAN)
  
!   
  IF (.NOT.L2D) THEN
    ZWORK = MXM(PRVCT)
    CALL ADVEC_WENO_K_2_MY(HLBCY, PUT, ZWORK, ZMEAN, TZHALO2_UT%HALO2)
    PRUS = PRUS - DYF(ZMEAN)
  END IF
!
  PRUS = PRUS - DZF(1,IKU,1,WENO_K_2_MZ(PUT, MXM(PRWCT)))
!
! V component
!
  IF (.NOT.L2D) THEN
    ZWORK = MYM(PRUCT)
    CALL ADVEC_WENO_K_2_MX(HLBCX, PVT, ZWORK, ZMEAN, TZHALO2_VT%HALO2)
    PRVS = PRVS - DXF(ZMEAN)
!   
    ZWORK = MYF(PRVCT)
    CALL ADVEC_WENO_K_2_VY(HLBCY, PVT, ZWORK, ZMEAN, TZHALO2_VT%HALO2)
    PRVS = PRVS - DYM(ZMEAN)
!
    PRVS = PRVS - DZF(1,IKU,1,WENO_K_2_MZ(PVT, MYM(PRWCT)))
  END IF
!
! W component
!
  ZWORK = MZM(1,IKU,1,PRUCT)
  CALL ADVEC_WENO_K_2_MX(HLBCX, PWT, ZWORK, ZMEAN, TZHALO2_WT%HALO2)
  PRWS = PRWS - DXF(ZMEAN)
!
  IF (.NOT.L2D) THEN
    ZWORK = MZM(1,IKU,1,PRVCT)
    CALL ADVEC_WENO_K_2_MY(HLBCY, PWT, ZWORK, ZMEAN, TZHALO2_WT%HALO2)
    PRWS = PRWS - DYF(ZMEAN)
  END IF
!
  PRWS = PRWS - DZM(1,IKU,1,WENO_K_2_WZ(PWT,MZF(1,IKU,1,PRWCT)))
!
!
CASE(5) ! WENO 5
!
! U component
!
  ZWORK = MXF(PRUCT)
  CALL ADVEC_WENO_K_3_UX(HLBCX, PUT, ZWORK, ZMEAN) 
  CALL GET_HALO(ZMEAN)! Update HALO 
  PRUS = PRUS - DXM(ZMEAN)
!
  IF (.NOT.L2D) THEN! 3D Case
   ZWORK = MXM(PRVCT)     
   CALL ADVEC_WENO_K_3_MY(HLBCY, PUT, ZWORK, ZMEAN)
   CALL GET_HALO(ZMEAN)! Update HALO 
   PRUS = PRUS - DYF(ZMEAN)
  END IF
!
  ZMEAN = WENO_K_3_MZ(PUT, MXM(PRWCT))
  CALL GET_HALO(ZMEAN)! Update HALO - maybe not necessary (T.Lunet)
  PRUS = PRUS - DZF(1,IKU,1,ZMEAN) 
!
! V component, only called in 3D case
!
  IF (.NOT.L2D) THEN
!
    ZWORK = MYM(PRUCT)
    CALL ADVEC_WENO_K_3_MX(HLBCX, PVT, ZWORK, ZMEAN) 
    CALL GET_HALO(ZMEAN)! Update HALO 
    PRVS = PRVS - DXF(ZMEAN)
!
    ZWORK = MYF(PRVCT)
    CALL ADVEC_WENO_K_3_VY(HLBCY, PVT, ZWORK, ZMEAN)
    CALL GET_HALO(ZMEAN)! Update HALO 
    PRVS = PRVS - DYM(ZMEAN)
!
    ZMEAN = WENO_K_3_MZ(PVT, MYM(PRWCT))
    CALL GET_HALO(ZMEAN)! Update HALO - maybe not necessary (T.Lunet)
    PRVS = PRVS - DZF(1,IKU,1,ZMEAN) 
!
  END IF
!
! W component
!
  ZWORK = MZM(1,IKU,1,PRUCT)
  CALL ADVEC_WENO_K_3_MX(HLBCX, PWT, ZWORK, ZMEAN)
  CALL GET_HALO(ZMEAN)! Update HALO
  PRWS = PRWS - DXF(ZMEAN)
!
  IF (.NOT.L2D) THEN! 3D Case
    ZWORK = MZM(1,IKU,1,PRVCT)
    CALL ADVEC_WENO_K_3_MY(HLBCY, PWT, ZWORK, ZMEAN)
    CALL GET_HALO(ZMEAN)! Update HALO
    PRWS = PRWS - DYF(ZMEAN)
  END IF
!
  ZMEAN = WENO_K_3_WZ(PWT,MZF(1,IKU,1,PRWCT))
  CALL GET_HALO(ZMEAN)! Update HALO - maybe not necessary (T.Lunet)
  PRWS = PRWS - DZM(1,IKU,1,ZMEAN)
!
!
END SELECT
!             ---------------------------------
!
END SUBROUTINE ADVECUVW_WENO_K

