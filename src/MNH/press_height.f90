!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!     ########################
      MODULE MODI_PRESS_HEIGHT
!     ########################
!
INTERFACE
!
FUNCTION PRESS_HEIGHT(PHEIGHT,PTHV,PPGROUND,PTHVGROUND,ZGROUND) RESULT(PP) 
REAL, DIMENSION(:), INTENT(IN) :: PHEIGHT ! height
REAL, DIMENSION(:), INTENT(IN) :: PTHV ! Thetav
REAL,               INTENT(IN) :: PPGROUND    ! pressure at ground level
REAL,               INTENT(IN) :: ZGROUND     ! orography
REAL,               INTENT(IN) :: PTHVGROUND  ! Thetav at at ground level
!
REAL, DIMENSION(SIZE(PHEIGHT))            :: PP   ! Pressure
!
END FUNCTION PRESS_HEIGHT
!
END INTERFACE
!
END MODULE MODI_PRESS_HEIGHT
!
!
!
!
!     ###########################################################################
      FUNCTION PRESS_HEIGHT(PHEIGHT,PTHV,PPGROUND,PTHVGROUND,PZGROUND) RESULT(PP) 
!     ###########################################################################
!
!!****  *PRESS_HEIGHT* - function to deduce  pressure  from height 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to deduce pressure from height
!     and potential virtual  temperature profiles by integration of the 
!     hydrostatic relation from the ground level. Pressure, potential virtual 
!     temperature  and height are known at ground level.  
!
!      
!
!!**  METHOD
!!    ------
!!
!!      The hydrostatic relation is integrated from ground level to
!!      upper levels : 
!!         
!!                   g        1
!!         d Exn = - ---  ----------  d z    
!!                   Cpd    -------z 
!!                          thetav
!!
!!       Then the pressure is deduced : 
!!        
!!          P =   Poo (Exn)**(Cpd /Rd)   (<----> Exn = (P/P00) ** (Rd/Cpd) )
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST        : contains physical constants
!!        XRD   : Gas constant for dry air
!!        XCPD  : Specific heat for dry air at constant pressure
!!        XG    : Gravity constant 
!!        XP00  : reference pressure
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (function PRESS_HEIGHT)
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/08/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST 
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
REAL, DIMENSION(:), INTENT(IN) :: PHEIGHT ! height
REAL, DIMENSION(:), INTENT(IN) :: PTHV    ! Thetav
REAL,               INTENT(IN) :: PPGROUND    ! pressure at ground level
REAL,               INTENT(IN) :: PZGROUND     ! orography
REAL,               INTENT(IN) :: PTHVGROUND  ! Thetav at at ground level
!
REAL, DIMENSION(SIZE(PHEIGHT)) :: PP   ! Pressure
!
!*       0.2   Declarations of local variables :
!
REAL                      :: ZGSCPD, ZCPDSRD ! g/Cpd, Cpd/Rd
REAL, DIMENSION(SIZE(PHEIGHT)) :: ZEXN            ! Exn = ( Poo/ P) ** Rd/Cpd
REAL                      :: ZEXNGROUND      ! Exner function  at ground level
INTEGER                   :: JK              ! Loop index
!-------------------------------------------------------------------------------
!
!*	 1.     COMPUTE PRESSURE
!	        ----------------
!
ZGSCPD     = XG / XCPD
ZCPDSRD   = XCPD / XRD
!
ZEXNGROUND = (PPGROUND/XP00) ** (1/ZCPDSRD) ! compute Exn at ground level
ZEXN(1) = ZEXNGROUND + 2. * ZGSCPD * (PZGROUND  - PHEIGHT(1))           &
                                    /(PTHVGROUND +PTHV(1))
!
DO JK = 2, SIZE(PP)
  ZEXN(JK) = ZEXN(JK-1) + 2. * ZGSCPD * ( PHEIGHT(JK-1)- PHEIGHT(JK) ) &
                                        /(PTHV(JK)+PTHV(JK-1))
END DO
!
PP(:) = XP00 *(ZEXN(:)**ZCPDSRD)
!
!-------------------------------------------------------------------------------
!
END FUNCTION PRESS_HEIGHT
