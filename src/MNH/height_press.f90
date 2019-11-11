!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_ideal 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ########################
      MODULE MODI_HEIGHT_PRESS
!     ########################
!
INTERFACE
!
FUNCTION HEIGHT_PRESS(PP,PTHV,PPGROUND,PTHVGROUND,ZGROUND) RESULT(PHEIGHT) 
REAL, DIMENSION(:), INTENT(IN) :: PP   ! Pressure
REAL, DIMENSION(:), INTENT(IN) :: PTHV ! Thetav
REAL,               INTENT(IN) :: PPGROUND    !  pressure at ground level
REAL,               INTENT(IN) :: ZGROUND     ! orography
REAL,               INTENT(IN) :: PTHVGROUND  ! Thetav at ground level
!
REAL, DIMENSION(SIZE(PP))            :: PHEIGHT ! height
!
END FUNCTION HEIGHT_PRESS
!
END INTERFACE
!
END MODULE MODI_HEIGHT_PRESS
!
!
!
!     ###########################################################################
      FUNCTION HEIGHT_PRESS(PP,PTHV,PPGROUND,PTHVGROUND,PZGROUND) RESULT(PHEIGHT) 
!     ###########################################################################
!
!!****  *HEIGHT_PRESS* - function to deduce height from pressure 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to deduce for an atmospheric column, the
!     height from pressure and potential virtual  temperature profiles by 
!     integration of the hydrostatic relation from the ground level.
!     Pressure, potential virtual temperature  and height are known at ground 
!     level.  
!
!!**  METHOD
!!    ------
!!      The reduced pressure Exn =  ( P/ Poo) ** Rd/Cpd is first computed.
!!    Then, the hydrostatic relation is integrated from ground level to
!!    upper levels : 
!!                     
!!                Cpd   -------z       
!!        dz  = - ---   thetav    d Exn  
!!                 g    
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
!!
!!    REFERENCE
!!    ---------
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
USE MODD_CST    ! declarative modules
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
REAL, DIMENSION(:), INTENT(IN) :: PP          ! Pressure
REAL, DIMENSION(:), INTENT(IN) :: PTHV        ! Thetav
REAL,               INTENT(IN) :: PPGROUND    !  pressure at ground level
REAL,               INTENT(IN) :: PZGROUND    ! orography
REAL,               INTENT(IN) :: PTHVGROUND  ! Thetav at ground level
!
REAL, DIMENSION(SIZE(PP))      :: PHEIGHT     ! height (result)
!
!*       0.2   Declarations of local variables :
!
REAL                      :: ZCPDOG, ZRDOCPD ! Cpd/g , Rd/Cpd
REAL, DIMENSION(SIZE(PP)) :: ZEXN            ! Exn = ( Poo/ P) ** Rd/Cpd
REAL                      :: ZEXNGROUND      ! Exn at ground level
INTEGER                   :: JK              ! Loop index
!-------------------------------------------------------------------------------
!
!*	 1.    COMPUTE HEIGHT 
!	        ------------
ZCPDOG     = XCPD / XG
ZRDOCPD    = XRD  / XCPD
!
ZEXN(:)    = (PP(:)   /XP00) ** ZRDOCPD ! compute Exn
ZEXNGROUND = (PPGROUND/XP00) ** ZRDOCPD ! compute Exn at ground level
!  
PHEIGHT(1) = PZGROUND + (PTHVGROUND +PTHV(1))*0.5 * ZCPDOG * (ZEXNGROUND -ZEXN(1))
!   
DO JK = 2, SIZE(PP)
  PHEIGHT(JK) = PHEIGHT(JK-1) + (PTHV(JK)+PTHV(JK-1))*0.5 * ZCPDOG  &
                  * ( ZEXN(JK-1)- ZEXN(JK) )
END DO
!-------------------------------------------------------------------------------
!
END FUNCTION HEIGHT_PRESS
