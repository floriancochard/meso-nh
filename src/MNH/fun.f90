!MNH_LIC Copyright 1994-2017 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############
      MODULE MODI_FUN
!     ###############
!
INTERFACE
!
FUNCTION FUNUYZ(KJU,KKU,PZHAT)
INTEGER, INTENT(IN)      :: KJU,KKU  ! array size in y and z directions 
REAL,DIMENSION(KKU),INTENT(IN) :: PZHAT
REAL, DIMENSION(KJU,KKU) :: FUNUYZ   ! U(y,z)
END FUNCTION FUNUYZ
!
FUNCTION FUNUY(KJU)
INTEGER, INTENT(IN)    :: KJU     ! array size in the y direction
REAL, DIMENSION(KJU)   :: FUNUY   ! U(y)
END FUNCTION FUNUY
!
FUNCTION FUNVXZ(KIU,KKU,PZHAT)
INTEGER, INTENT(IN)      :: KIU,KKU  ! array size in x and z directions 
REAL,DIMENSION(KKU),INTENT(IN) :: PZHAT
REAL, DIMENSION(KIU,KKU) :: FUNVXZ   ! V(x,z)
END FUNCTION FUNVXZ
!
FUNCTION FUNVX(KIU)
INTEGER, INTENT(IN)    :: KIU     ! array size in the x direction
REAL, DIMENSION(KIU)   :: FUNVX   ! V(x)
END FUNCTION FUNVX
!
END INTERFACE
!
END MODULE MODI_FUN
!
!
!     ########################
      FUNCTION FUNUYZ(KJU,KKU,PZHAT)
!     ########################
!
!!****  *FUNUYZ* - function  of the coordinates Y, Z
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to return a 2D array  of size KJU*KKU
!     containing the x component of the wind U and thus its variation along the 
!     y and z directions.
!
!!**  METHOD
!!    ------
!!      
!!      An analytical function is used to prescribe the U  variation along the 
!!    y and z directions ( J and K index). 
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_GRID1 : contains grid variables
!!        XYHAT :  Position y in the conformal
!!                 plane or on the cartesian plane
!!        XZHAT :  Position z without orography
!!
!!    REFERENCE
!!    ---------
!!
!!     NONE
!!
!!    AUTHOR
!!    ------
!!	J.Stein         * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    2/12/94
!!      P.Jabouille 28/08/00  parallelization
!!      G.Tanguy    26/10/10  add  PZHAT to aplly funuyz 
!!                            to a field on a vertical grid different
!!                            from XZHAT (mixed grid)
!!      P.Wautelet  19/10/2017 removed extern_userio.f90
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_GATHER_ll
USE MODD_GRID_n
USE MODD_DIM_n
USE MODD_PARAMETERS
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
INTEGER, INTENT(IN)      :: KJU,KKU  ! array size in y and z directions 
REAL,DIMENSION(KKU),INTENT(IN) :: PZHAT
REAL, DIMENSION(KJU,KKU) :: FUNUYZ   ! U(y,z)
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IJ0,IK0 ! Jet center 
INTEGER :: JJ ,JK  ! Loop index
INTEGER :: IINFO_ll, IJU_ll  ! parallel variables
REAL    :: ZWIDTHY ! Width of the jet along the y direction
REAL    :: ZWIDTHZ ! Width of the jet along the z direction
REAL, DIMENSION(:), ALLOCATABLE :: ZYHAT_ll
!
!-------------------------------------------------------------------------------
!
!*	 1.     COMPUTE FUNUYZ
!	        -------------
! 
IJU_ll=NJMAX_ll+2*JPHEXT
IJ0=IJU_ll/2
IK0=KKU/2
ALLOCATE(ZYHAT_ll(IJU_ll))
CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IINFO_ll)
ZWIDTHY =ZYHAT_ll(IJ0+IJU_ll/5)-ZYHAT_ll(IJ0)
ZWIDTHZ =PZHAT(IK0+KKU/5)-PZHAT(IK0)
DO JJ = 1,KJU-1
  DO JK = 1,KKU
    FUNUYZ(JJ,JK) = 1./COSH(                                    &
     (( (XYHAT(JJ)+XYHAT(JJ+1))*0.5-ZYHAT_ll(IJ0))/ZWIDTHY) **2 &
    +((  PZHAT(JK)                 -   PZHAT(IK0))/ZWIDTHZ) **2 )
  END DO
END DO
DEALLOCATE(ZYHAT_ll)
FUNUYZ(KJU,:)=2.*FUNUYZ(KJU-1,:)-FUNUYZ(KJU-2,:) !simple extrapolation  
                                                 !for the last point
!
!-------------------------------------------------------------------------------
!
END FUNCTION FUNUYZ
!
!     ###################
      FUNCTION FUNUY(KJU)
!     ###################
!
!!****  *FUNUY* - function  of the coordinate Y
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to return a 1D array  of size KJU
!     containing the x component of the wind U and thus its variation along the 
!     y direction
!
!!**  METHOD
!!    ------
!!      
!!      An analytical function is used to prescribe the U  variation along the 
!!    y direction ( J index). The result is located at the W-point. 
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_GRID1 : contains grid variables
!!        XYHAT :  Position y in the conformal
!!                 plane or on the cartesian plane
!!
!!    REFERENCE
!!    ---------
!!
!!      None
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                                         30/08/94
!!      change the way to get the size of the function   2 /12/94   (J.Stein)
!!      parallelization                                  04/09/00   (P Jabouille)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_GRID_n
USE MODD_DIM_n
USE MODD_PARAMETERS
USE MODE_GATHER_ll
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
INTEGER, INTENT(IN)    :: KJU     ! array size in the y direction
REAL, DIMENSION(KJU)   :: FUNUY   ! U(y)
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IJ0      ! Jet center
INTEGER :: JJ       ! Loop index
INTEGER :: IINFO_ll, IJU_ll  ! parallel variables
REAL    :: ZWIDTH   ! Width of the jet
REAL, DIMENSION(:), ALLOCATABLE ::  ZYHAT_ll
!-------------------------------------------------------------------------------
!
!*	 1.     COMPUTE FUNUY
!	        -------------
IJU_ll=NJMAX_ll+2*JPHEXT
IJ0=IJU_ll/2
ALLOCATE(ZYHAT_ll(IJU_ll))
CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IINFO_ll)
ZWIDTH=ZYHAT_ll(IJ0+IJU_ll/5)-ZYHAT_ll(IJ0)
DO JJ = 1,KJU-1
 FUNUY(JJ) = 1./COSH(((XYHAT(JJ)+XYHAT(JJ+1))*0.5-ZYHAT_ll(IJ0))/ZWIDTH)
END DO
FUNUY(KJU)=2.*FUNUY(KJU-1)-FUNUY(KJU-2) !simple extrapolation for the last point
!
!-------------------------------------------------------------------------------
!
END FUNCTION FUNUY
!
!     ########################
      FUNCTION FUNVXZ(KIU,KKU,PZHAT)
!     ########################
!
!!****  *FUNVXZ* - function  of the coordinates X, Z
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to return a 2D array  of size KIU*KKU
!     containing the y component of the wind V and thus its variation along the 
!     x and z direction.
!
!!**  METHOD
!!    ------
!!      
!!      An analytical function is used to prescribe the V  variation along the 
!!    x and z directions ( I and K index). 
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_GRID1 : contains grid variables
!!        XXHAT :  Position x in the conformal
!!                 plane or on the cartesian plane
!!        XZHAT :  Position z
!!
!!    REFERENCE
!!    ---------
!!
!!     NONE
!!
!!    AUTHOR
!!    ------
!!	J.Stein         * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     2/12/94
!!      P.Jabouille 28/08/00  parallelization
!!      G.Tanguy    26/10/10  add  PZHAT to aplly funuyz 
!!                            to a field on a vertical grid different
!!                            from XZHAT (mixed grid)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_GATHER_ll
USE MODD_GRID_n
USE MODD_DIM_n
USE MODD_PARAMETERS
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
INTEGER, INTENT(IN)      :: KIU,KKU  ! array size in x and z directions 
REAL,DIMENSION(KKU),INTENT(IN) :: PZHAT
REAL, DIMENSION(KIU,KKU) :: FUNVXZ   ! V(x,z)
!
!*       0.2   Declarations of local variables :
!
INTEGER :: II0,IK0   ! Jet center 
INTEGER :: JI,JK     ! Loop index
INTEGER :: IINFO_ll, IIU_ll  ! parallel variables
REAL    :: ZWIDTHX ! Width of the jet along the x direction
REAL    :: ZWIDTHZ ! Width of the jet along the z direction
REAL, DIMENSION(:), ALLOCATABLE :: ZXHAT_ll
!
!-------------------------------------------------------------------------------
!
!*	 1.     COMPUTE FUNVXZ
!	        -------------
! 
IIU_ll=NIMAX_ll+2*JPHEXT
II0=IIU_ll/2
IK0=KKU/2
ALLOCATE(ZXHAT_ll(IIU_ll))
CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IINFO_ll)
ZWIDTHX=ZXHAT_ll(II0+IIU_ll/5)-ZXHAT_ll(II0)
ZWIDTHZ=PZHAT(IK0+KKU/5)-PZHAT(IK0)
DO JI = 1,KIU-1
  DO JK = 1,KKU
    FUNVXZ(JI,JK) = 1./COSH(                                   &
     (( (XXHAT(JI)+XXHAT(JI+1))*0.5-ZXHAT_ll(II0))/ZWIDTHX)**2 &
    +((  PZHAT (JK)                 - PZHAT (IK0))/ZWIDTHZ)**2 )
  END DO
END DO
FUNVXZ(KIU,:)=2.*FUNVXZ(KIU-1,:)-FUNVXZ(KIU-2,:) !simple extrapolation for the last point
!
!-------------------------------------------------------------------------------
!
END FUNCTION FUNVXZ
!
!
!     ###################
      FUNCTION FUNVX(KIU)
!     ###################
!
!!****  *FUNVX* - function  of the coordinate X
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to return a 1D array  of size KIU
!     containing the y component of the wind V and thus its variation along the 
!     x direction
!
!!**  METHOD
!!    ------
!!      
!!      An analytical function is used to prescribe the V  variation along the 
!!    x direction ( I index). 
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_GRID1 : contains grid variables
!!        XXHAT :  x position  in the conformal
!!                 plane or on the cartesian plane
!!
!!    REFERENCE
!!    ---------
!!
!!      None
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                                         30/08/94
!!      change the way to get the size and the function   2 /12/94   (J.Stein)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_GATHER_ll
USE MODD_GRID_n
USE MODD_DIM_n
USE MODD_PARAMETERS
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
INTEGER, INTENT(IN)    :: KIU     ! array size in the x direction
REAL, DIMENSION(KIU)   :: FUNVX   ! V(x)
!
!*       0.2   Declarations of local variables :
!
INTEGER :: II0   ! Jet center
INTEGER :: JI    ! Loop index
INTEGER :: IINFO_ll, IIU_ll  ! parallel variables
REAL    :: ZWIDTH ! Width of the jet
REAL, DIMENSION(:), ALLOCATABLE :: ZXHAT_ll
!-------------------------------------------------------------------------------
!
!*	 1.     COMPUTE FUNUY
!	        -------------
IIU_ll=NIMAX_ll+2*JPHEXT
II0=IIU_ll/2
ALLOCATE(ZXHAT_ll(IIU_ll))
CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IINFO_ll)
ZWIDTH=ZXHAT_ll(II0+IIU_ll/5)-ZXHAT_ll(II0)
DO JI = 1,KIU
 FUNVX(JI)=1./COSH(((XXHAT(JI)+XXHAT(JI))*0.5-ZXHAT_ll(II0))/ZWIDTH)
END DO
DEALLOCATE(ZXHAT_ll)
!-------------------------------------------------------------------------------
!
END FUNCTION FUNVX
