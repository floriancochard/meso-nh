!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 interpol 2006/05/18 13:07:25
!-----------------------------------------------------------------
!######################## 
MODULE MODI_INI_BIKHARDT_n
!########################
!
INTERFACE
!
SUBROUTINE INI_BIKHARDT_n (KDXRATIO,KDYRATIO,KMI)
!
INTEGER, INTENT(IN) :: KDXRATIO            !  x and y-direction resolution RATIO
INTEGER, INTENT(IN) :: KDYRATIO            !      between model 2 and model 1
INTEGER, INTENT(IN) :: KMI                 !  model index to work with 
!
END SUBROUTINE INI_BIKHARDT_n
!
END INTERFACE
!
END MODULE MODI_INI_BIKHARDT_n
!
!
!     #################################################
      SUBROUTINE INI_BIKHARDT_n (KDXRATIO,KDYRATIO,KMI)
!     #################################################
!
!!****  *INI_BIKHARDT_n * - subroutine to initialize coefficients for
!!                          Bikhardt interpolation
!!
!!    PURPOSE
!!    -------
!!
!!      Initializes coefficients for Bikhardt interpolation.
!!      In case of gridnesting or spawning, _n refers to the INNER model.
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      Module MODD_BIKHARDT_n : contains Bikhardt coefficients
!!
!!    REFERENCE
!!    ---------
!!
!!       Routine INI_BIKHARDT_n (Book2 of the documentation)
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       V. Masson     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     10/06/96 
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_BIKHARDT_n
USE MODE_MODELN_HANDLER
!
IMPLICIT NONE
!
!*       0.1    Declarations of dummy arguments :
!
!
INTEGER, INTENT(IN) :: KDXRATIO            !  x and y-direction resolution RATIO
INTEGER, INTENT(IN) :: KDYRATIO            !      between model 2 and model 1
INTEGER, INTENT(IN) :: KMI                 !  model index to work with 

!
!*       0.2    Declarations of local variables :
!
REAL                :: ZX, ZY
INTEGER             :: JI             ! Loop index in x direction
INTEGER             :: JJ             ! Loop index in y direction      
INTEGER             :: IMI            ! current model index
!
!-------------------------------------------------------------------------------
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL GOTO_MODEL(KMI)
!
!*       1.     Allocation of Bikhardt interpolation coefficients :
!
ALLOCATE (XBMX1(KDXRATIO))
ALLOCATE (XBMX2(KDXRATIO))
ALLOCATE (XBMX3(KDXRATIO))
ALLOCATE (XBMX4(KDXRATIO))
ALLOCATE (XBMY1(KDYRATIO))
ALLOCATE (XBMY2(KDYRATIO))
ALLOCATE (XBMY3(KDYRATIO))
ALLOCATE (XBMY4(KDYRATIO))
ALLOCATE (XBFX1(KDXRATIO))
ALLOCATE (XBFX2(KDXRATIO))
ALLOCATE (XBFX3(KDXRATIO))
ALLOCATE (XBFX4(KDXRATIO))
ALLOCATE (XBFY1(KDYRATIO))
ALLOCATE (XBFY2(KDYRATIO))
ALLOCATE (XBFY3(KDYRATIO))
ALLOCATE (XBFY4(KDYRATIO))
!
!*       2.     Bikhardt interpolation coefficients computation :
!
DO JI = 1,KDXRATIO
  ZX = FLOAT(JI-1)/FLOAT(KDXRATIO)
  XBFX1(JI) = -0.5*ZX*ZX*ZX  +    ZX*ZX  -0.5*ZX
  XBFX2(JI) =  1.5*ZX*ZX*ZX  -2.5*ZX*ZX           +1.
  XBFX3(JI) = -1.5*ZX*ZX*ZX  +2.0*ZX*ZX  +0.5*ZX
  XBFX4(JI) =  0.5*ZX*ZX*ZX  -0.5*ZX*ZX
!
  IF (MOD(KDXRATIO,2).EQ.0)  ZX = ZX + .5/FLOAT(KDXRATIO)
  XBMX1(JI) = -0.5*ZX*ZX*ZX  +    ZX*ZX  -0.5*ZX
  XBMX2(JI) =  1.5*ZX*ZX*ZX  -2.5*ZX*ZX           +1.
  XBMX3(JI) = -1.5*ZX*ZX*ZX  +2.0*ZX*ZX  +0.5*ZX
  XBMX4(JI) =  0.5*ZX*ZX*ZX  -0.5*ZX*ZX
!
END DO
!
DO JJ = 1,KDYRATIO
  ZY = FLOAT(JJ-1)/FLOAT(KDYRATIO)
  XBFY1(JJ) = -0.5*ZY*ZY*ZY  +    ZY*ZY  -0.5*ZY
  XBFY2(JJ) =  1.5*ZY*ZY*ZY  -2.5*ZY*ZY           +1.
  XBFY3(JJ) = -1.5*ZY*ZY*ZY  +2.0*ZY*ZY  +0.5*ZY
  XBFY4(JJ) =  0.5*ZY*ZY*ZY  -0.5*ZY*ZY
!
  IF (MOD(KDYRATIO,2).EQ.0)  ZY = ZY + .5/FLOAT(KDYRATIO)
  XBMY1(JJ) = -0.5*ZY*ZY*ZY  +    ZY*ZY  -0.5*ZY
  XBMY2(JJ) =  1.5*ZY*ZY*ZY  -2.5*ZY*ZY           +1.
  XBMY3(JJ) = -1.5*ZY*ZY*ZY  +2.0*ZY*ZY  +0.5*ZY
  XBMY4(JJ) =  0.5*ZY*ZY*ZY  -0.5*ZY*ZY
!
END DO
!
CALL GOTO_MODEL(IMI)
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_BIKHARDT_n
