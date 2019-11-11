!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2007/02/22 09:35:07
!-----------------------------------------------------------------
!    #######################
     MODULE MODI_INI_RADCONF
!    #######################
!
INTERFACE
!
     SUBROUTINE INI_RADCONF (HLW,KSWB,OSUBG_COND)
!
CHARACTER (LEN=*), INTENT (IN) :: HLW  ! choice in LW scheme
INTEGER,           INTENT (IN) :: KSWB ! number of SW band 
LOGICAL,           INTENT(IN)  :: OSUBG_COND ! Switch for sub-grid condensation
!
END SUBROUTINE INI_RADCONF
!
END INTERFACE
!
END MODULE MODI_INI_RADCONF
!
!
!     ############################################
      SUBROUTINE INI_RADCONF (HLW,KSWB,OSUBG_COND)
!     ############################################
!
!     PURPOSE.
!     --------
!     Initialisation of grid-independant quantities for ECMWF radiation routines
!
!     METHOD.
!     -------
!     The values of variable contained in specific radiation module "yo===" are mainly      
!     initialised through the call of the external subroutines. However, some variables  
!     which are initially fixed in the radiation driver of ECMWF model, has been explicitely  
!     set in the present subroutine (for MesoNH interface). This was done only for the variable 
!     used later in  MesoNH adaptation of ECMWF code. Consequently,  some variable which are 
!     declared in "yo---" module are not defined , but they have no futher impact in MNH context.   
!
!     EXTERNALS.
!     ----------
!     SEE the corresponding call !
!
!     AUTHORS.
!     --------
!     F.Solmon, March 2002
! MODIF
!!                   02/2018 Q.Libois ECRAD
!
!----------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE OYOMCST   , ONLY : RG, RCPD, RDAY      
!
USE OYOERAD   , ONLY : NMODE, NOVLP, NLW, NSW, &
                      NTSW, LRRTM, LONEWSW,LRADIP, LRADLP, LINHOM    
!   
USE YOERDU   , ONLY : NUAER, NTRAER, RCDAY, R10E, REPLOG, &
                      REPSC, REPSCO, REPSCQ, REPSCT, REPSCW, DIFF
!
USE MODI_OSURDI
USE MODI_SULWN
USE MODI_OSURRTAB
USE MODI_OSURRTPK
USE MODI_OSURRTRF
USE MODI_OSURRTFTR
!USE MODI_ORRTM_KGB1
!USE MODI_ORRTM_KGB2
!USE MODI_ORRTM_KGB3
!USE MODI_ORRTM_KGB4
!USE MODI_ORRTM_KGB5
!USE MODI_ORRTM_KGB6
!USE MODI_ORRTM_KGB7
!USE MODI_ORRTM_KGB8
!USE MODI_ORRTM_KGB9
!USE MODI_ORRTM_KGB10
!USE MODI_ORRTM_KGB11
!USE MODI_ORRTM_KGB12
!USE MODI_ORRTM_KGB13
!USE MODI_ORRTM_KGB14
!USE MODI_ORRTM_KGB15
!USE MODI_ORRTM_KGB16
!USE MODI_ORRTM_INIT_140GP
USE MODI_SUCST
USE MODI_SUSWN
USE MODI_SUAERL
USE MODI_SUAERH
USE MODI_SUAERSN
USE MODI_SUCLOPN
USE MODI_SUCLD
!
IMPLICIT NONE
! 
!*      0.1  declarations of arguments
!
CHARACTER (LEN=*), INTENT (IN) :: HLW 

INTEGER,           INTENT (IN) :: KSWB
LOGICAL,           INTENT(IN)  :: OSUBG_COND ! Switch for sub-grid condensation
!
!-------------------------------------------------------------------
!
!Cloud properties
!
LRADIP =.TRUE. ! used in default calculation of effective radius
LRADLP =.TRUE. !
LINHOM =.FALSE. !
!
! cloud overlap assumption
!
IF (.NOT. OSUBG_COND) THEN
  NOVLP=5
ELSE IF (HLW =='MORC') THEN
  NOVLP=5
ELSE
  NOVLP=6   
END IF
!
! "MORCR LW" code  
!
NMODE  = 0  !   
NLW    = 6  !
!NUAER  = 31 !
!NTRAER = 19 !
!
! "RRTM LW"
!
IF(HLW =='RRTM') THEN
  LRRTM  =.TRUE.
ELSE
  LRRTM =.FALSE.
END IF
!
! SW code
!
LONEWSW=.TRUE.
NTSW   = 6  !
NSW    = KSWB
!
! default values used in the radiation code
!
REPSC  = 1.E-12
REPSCO = 1.E-12
REPSCQ = 1.E-12
REPSCT = 1.E-12
REPSCW = 1.E-12
REPLOG = 1.E-12
!
! CALL to the ECMW initialisation routines
!
CALL OSURDI
CALL SULWN
CALL OSURRTAB
CALL OSURRTPK 
CALL OSURRTRF
CALL OSURRTFTR
!
IF(LRRTM) THEN
  CALL ORRTM_KGB1
  CALL ORRTM_KGB2
  CALL ORRTM_KGB3
  CALL ORRTM_KGB4
  CALL ORRTM_KGB5
  CALL ORRTM_KGB6
  CALL ORRTM_KGB7
  CALL ORRTM_KGB8
  CALL ORRTM_KGB9
  CALL ORRTM_KGB10
  CALL ORRTM_KGB11
  CALL ORRTM_KGB12
  CALL ORRTM_KGB13
  CALL ORRTM_KGB14
  CALL ORRTM_KGB15
  CALL ORRTM_KGB16
!
!*** mji ***
!Initialization routine for RRTM
! Reduce absorption coefficient data from 256 to 140 g-points
!
  CALL ORRTM_INIT_140GP
END IF
!
!- basic constants      
!
CALL SUCST 
!
RCDAY= RDAY * RG/RCPD
R10E=  0.4342945000000000
DIFF=  1.66!! 
!
CALL SUSWN    ( NTSW, NSW )
CALL SUAERL
CALL SUAERH
CALL SUAERSN  ( NTSW, NSW )
CALL SUCLOPN  ( NTSW, NSW, 60 )
!
!  constant for clouds
!
CALL SUCLD 
!
END SUBROUTINE INI_RADCONF








