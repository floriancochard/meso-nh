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
      MODULE MODI_GRADIENT_UW
!     #######################
!
INTERFACE
!
!     
FUNCTION GX_UW_W(KKA,KKU,KL,PA,PDXX,PDZZ,PDZX)      RESULT(PGX_UW_W)
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_UW_W ! result W point
!
END FUNCTION GX_UW_W
!
!     
FUNCTION GZ_UW_U(KKA,KKU,KL,PA,PDZZ)      RESULT(PGZ_UW_U)
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_UW_U ! result U point
!
END FUNCTION GZ_UW_U
!
END INTERFACE
!
END MODULE MODI_GRADIENT_UW
!
!
!
!
!     #########################################################
      FUNCTION GX_UW_W(KKA,KKU,KL,PA,PDXX,PDZZ,PDZX)      RESULT(PGX_UW_W)
!     #########################################################
!
!!****  *GX_UW_W* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian X
!!                          direction for a variable placed at the 
!!                          UW point and the result is placed at
!!                          the W point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the X cartesian direction for a field PA placed at the 
!     UW point. The result is placed at the W point.
!
!
!                        (          __________________z )
!                        (          ( _____________x  ) )
!                        (          ( ___z            ) ) 
!                    1   (          ( d*zx  dzf(PA)   ) )
!      PGX_UW_W =  ----  (dxf(PA) - (-----------------) )
!                  ___x,z(          (       ___z      ) )
!                  d*xx  (          (       d*zz      ) )
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
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise*
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_UW_W ! result W point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GX_UW_W
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGX_UW_W(:,:,:)= ( DXF(PA)        -                          &
                    MZF(KKA,KKU,KL, MXF(MZF(KKA,KKU,KL, PDZX)*DZF(KKA,KKU,KL,PA)) / MZF(KKA,KKU,KL,PDZZ) )  &
                  ) / MXF(MZM(KKA,KKU,KL,PDXX))
ELSE
  PGX_UW_W(:,:,:)= DXF(PA) /  MXF(MZM(KKA,KKU,KL,PDXX))
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GX_UW_W
!
! 
!     ###############################################
      FUNCTION GZ_UW_U(KKA,KKU,KL,PA,PDZZ)      RESULT(PGZ_UW_U)
!     ###############################################
!
!!****  *GZ_UW_U* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the 
!!                          UW point and the result is placed at
!!                          the U point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Z cartesian direction for a field PA placed at the 
!     UW point. The result is placed at the U point.
!
!
!                    dzf(PA) 
!      PGZ_UW_U =   --------  
!                    ____x,z
!                    d*zz   
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
!!      MXM,MZF    : Shuman functions (mean operators)
!!      DZF        : Shuman functions (finite difference operators)
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
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the UW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_UW_U ! result U point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_UW_U
!              ---------------------
!
PGZ_UW_U(:,:,:)= DZF(KKA,KKU,KL,PA) / MXM(MZF(KKA,KKU,KL,PDZZ))
!
!----------------------------------------------------------------------------
!
END FUNCTION GZ_UW_U
