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
!     #############################
      MODULE MODI_THETAVPU_THETAVPM
!     #############################
!
INTERFACE
!
FUNCTION THETAVPU_THETAVPM(PPM,PPU,PTHVM) RESULT(PTHVU) 
!
REAL, DIMENSION(:), INTENT(IN) :: PPM   ! Pressure ot mass levels
REAL, DIMENSION(:), INTENT(IN) :: PPU   ! Pressure at wind levels
REAL, DIMENSION(:), INTENT(IN) :: PTHVM ! Virtual potential temperature at mass
                                        ! levels
!
REAL, DIMENSION(SIZE(PPU)) :: PTHVU     ! Virtual potential temperature at wind
                                        ! levels
END FUNCTION THETAVPU_THETAVPM
!
END INTERFACE
!
END MODULE MODI_THETAVPU_THETAVPM
!
!
!
!     #########################################################
      FUNCTION THETAVPU_THETAVPM(PPM,PPU,PTHVM) RESULT(PTHVU) 
!     #########################################################
!
!!****  *THETAVPU_THETAVPM* - function to deduce virtual potential  temperature 
!!                            at  given pressure levels
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to interpolate (or extrapolate) 
!     virtual potential temperature at pressure PPU (at wind level), 
!     from   virtual potential temperature at pressure PPM (mass level).
!      
!     NB: PPU(1) is always smaller or equal to PPM(1)
!
!!**  METHOD
!!    ------
!!      The virtual potential temperature is linearly interpolated in Logarithm
!!    of pressure. 
!!      When  wind levels are higher than mass levels, the  virtual potential 
!!   temperature is extrapolated at constant gradient :
!!        d thetav 
!!        -------- = constant
!!        d Log P
!!     
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (function THETAVPU_THETAVPM)
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
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
REAL, DIMENSION(:), INTENT(IN) :: PPM   ! Pressure ot mass levels
REAL, DIMENSION(:), INTENT(IN) :: PPU   ! Pressure at wind levels
REAL, DIMENSION(:), INTENT(IN) :: PTHVM ! Virtual potential temperature at mass
                                        ! levels
!
REAL, DIMENSION(SIZE(PPU)) :: PTHVU     ! Virtual potential temperature at wind
                                        ! levels
!
!*       0.2   Declarations of local variables :
!
INTEGER  :: ILEVELU,ILEVELM        ! number of wind levels and mass levels
REAL :: ZDPSDPM,ZDP1SDPM,ZDP2SDPM  ! working variables
INTEGER :: JK,JKM ! loop indexes 
!-------------------------------------------------------------------------------
!
!*	 1.     COMPUTE THETAV AT WIND LEVELS
!	        ------------------------------
!
ILEVELU=SIZE(PPU)  ! Retrieve number of wind levels
ILEVELM=SIZE(PPM)  ! Retrieve number of mass levels
!
DO JK = 1,ILEVELU
  IF (PPU(JK) <= PPM(ILEVELM)) THEN     ! Extrapolation for upper wind levels
    ZDPSDPM  = LOG(PPU(JK)/PPM(ILEVELM)) &
             / LOG(PPM(ILEVELM)/PPM(ILEVELM-1))  
    PTHVU(JK) = PTHVM(ILEVELM) + (PTHVM(ILEVELM)-PTHVM(ILEVELM-1)) * ZDPSDPM                          
  ELSE 
    DO JKM=1,ILEVELM-1                  ! Interpolation 
      IF ((PPU(JK) <= PPM(JKM)).AND.(PPU(JK) > PPM(JKM+1)))  THEN   
        ZDP1SDPM = LOG(PPU(JK)/PPM(JKM))/LOG(PPM(JKM+1)/PPM(JKM))
        ZDP2SDPM = 1. - ZDP1SDPM
        PTHVU(JK) = (PTHVM(JKM) * ZDP2SDPM) + (PTHVM(JKM+1) * ZDP1SDPM)
      END IF 
    END DO       
  END IF 
END DO
!-------------------------------------------------------------------------------
!
END FUNCTION THETAVPU_THETAVPM
