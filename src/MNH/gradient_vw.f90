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
      MODULE MODI_GRADIENT_VW
!     #######################
!
INTERFACE
!
!     
FUNCTION GY_VW_W(KKA,KKU,KL,PA,PDYY,PDZZ,PDZY)      RESULT(PGY_VW_W)
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the VW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_VW_W ! result W point
!
END FUNCTION GY_VW_W
!
!     
FUNCTION GZ_VW_V(KKA,KKU,KL,PA,PDZZ)      RESULT(PGZ_VW_V)
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the VW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_VW_V ! result V point
!
END FUNCTION GZ_VW_V
!
END INTERFACE
!
END MODULE MODI_GRADIENT_VW
!
!
!
!
!     #########################################################
      FUNCTION GY_VW_W(KKA,KKU,KL,PA,PDYY,PDZZ,PDZY)      RESULT(PGY_VW_W)
!     #########################################################
!
!!****  *GY_VW_W* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          VW point and the result is placed at
!!                          the W point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     VW point. The result is placed at the W point.
!
!
!                        (          __________________z )
!                        (          ( _____________y  ) )
!                        (          ( ___z            ) ) 
!                    1   (          ( d*zy  dzf(PA)   ) )
!      PGY_VW_W =  ----  (dyf(PA) - (-----------------) )
!                  ___y,z(          (       ___z      ) )
!                  d*yy  (          (       d*zz      ) )
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
!!      MYF,MZF         : Shuman functions (mean operators)
!!      DYF,DZF         : Shuman functions (finite difference operators)
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
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the VW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_VW_W ! result W point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_VW_W
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGY_VW_W(:,:,:)= ( DYF(PA)        -                          &
                    MZF(KKA,KKU,KL, MYF(MZF(KKA,KKU,KL,PDZY)*DZF(KKA,KKU,KL,PA)) / MZF(KKA,KKU,KL,PDZZ) )  &
                  ) / MYF(MZM(KKA,KKU,KL,PDYY))
ELSE
  PGY_VW_W(:,:,:)= DYF(PA) /  MYF(MZM(KKA,KKU,KL,PDYY))
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GY_VW_W
!
! 
!     ###############################################
      FUNCTION GZ_VW_V(KKA,KKU,KL,PA,PDZZ)      RESULT(PGZ_VW_V)
!     ###############################################
!
!!****  *GZ_VW_V* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the 
!!                          VW point and the result is placed at
!!                          the V point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Z cartesian direction for a field PA placed at the 
!     VW point. The result is placed at the V point.
!
!
!                    dzf(PA) 
!      PGZ_VW_V =   --------  
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
!!      MYM,MZF    : Shuman functions (mean operators)
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
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the VW point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_VW_V ! result V point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_VW_V
!              ---------------------
!
PGZ_VW_V(:,:,:)= DZF(KKA,KKU,KL,PA) / MYM(MZF(KKA,KKU,KL,PDZZ))
!
!----------------------------------------------------------------------------
!
END FUNCTION GZ_VW_V
