!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 operators 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_GRADIENT_UV
!     #######################
!
INTERFACE
!
!     
FUNCTION GX_UV_V(KKA,KKU,KL,PA,PDXX,PDZZ,PDZX)      RESULT(PGX_UV_V)
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UV point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_UV_V ! result V point
!
END FUNCTION GX_UV_V
!
!     
FUNCTION GY_UV_U(KKA,KKU,KL,PA,PDYY,PDZZ,PDZY)      RESULT(PGY_UV_U)
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UV point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_UV_U ! result U point
!
END FUNCTION GY_UV_U
!
END INTERFACE
!
END MODULE MODI_GRADIENT_UV
!
!
!
!
!     #########################################################
      FUNCTION GX_UV_V(KKA,KKU,KL,PA,PDXX,PDZZ,PDZX)      RESULT(PGX_UV_V)
!     #########################################################
!
!!****  *GX_UV_V* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian X
!!                          direction for a variable placed at the 
!!                          UV point and the result is placed at
!!                          the V point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the X cartesian direction for a field PA placed at the 
!     UV point. The result is placed at the V point.
!
!
!                        (          ________________z )
!                        (          (_____________x ) )
!                    1   (          (___y (dzm(PA)) ) ) 
!      PGX_UV_V =  ----  (dxf(PA) - (d*zx (-------) ) )
!                  ___x,y(          (     ( ___y  ) ) )
!                  d*xx  (          (     ( d*zz  ) ) )     
!
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MXF,MZF         : Shuman functions (mean operators)
!!      DXF,DZF         : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (GRAD_CAR operators)
!!      A Turbulence scheme for the Meso-NH model (Chapter 6)
!!
!!    AUTHOR
!!    ------
!!      Valery Masson      *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/07/02
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN
USE MODD_CONF
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UV point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_UV_V ! result V point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GX_UV_V
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGX_UV_V(:,:,:)= ( DXF(PA)        -                         &
                    MZF(KKA,KKU,KL, MXF( MYM(PDZX)*DZM(KKA,KKU,KL,PA)/MYM(PDZZ) ) ) &
                   ) / MXF(MYM(PDXX))
ELSE
  PGX_UV_V(:,:,:)= DXF(PA) /  MXF(MYM(PDXX))
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GX_UV_V
!
! 
!     #########################################################
      FUNCTION GY_UV_U(KKA,KKU,KL,PA,PDYY,PDZZ,PDZY)      RESULT(PGY_UV_U)
!     #########################################################
!
!!****  *GY_UV_U* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          UV point and the result is placed at
!!                          the U point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     UV point. The result is placed at the U point.
!
!
!
!                        (          ________________z )
!                        (          (_____________y ) )
!                    1   (          (___x (dzm(PA)) ) ) 
!      PGY_UV_U =  ----  (dyf(PA) - (d*zy (-------) ) )
!                  ___x,y(          (     ( ___x  ) ) )
!                  d*yy  (          (     ( d*zz  ) ) )     
!
!                        (          __________________z )
!                        (          (       ________y ) )
!                    1   (          (___x,y (dzm(PA)) ) ) 
!      PGY_UV_U =  ----  (dyf(PA) - (d*zy   (-------) ) )
!                  ___x,y(          (       ( ___x  ) ) )
!                  d*yy  (          (       ( d*zz  ) ) )     
!       
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MXM,MYM,MZF     : Shuman functions (mean operators)
!!      DYM,DZM         : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (GRAD_CAR operators)
!!      A Turbulence scheme for the Meso-NH model (Chapter 6)
!!
!!    AUTHOR
!!    ------
!!      Valery Masson      *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/07/02
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN
USE MODD_CONF
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UV point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_UV_U ! result U point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_UV_U
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGY_UV_U(:,:,:)= ( DYF(PA)        -                         &
                    MZF(KKA,KKU,KL, MYF( MXM(PDZY)*DZM(KKA,KKU,KL,PA)/MXM(PDZZ) ) ) &
                   ) / MYF(MXM(PDYY))
ELSE
  PGY_UV_U(:,:,:)= DYF(PA) / MYF(MXM(PDYY))
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GY_UV_U
