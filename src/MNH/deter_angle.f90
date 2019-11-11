!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 diag 2006/05/18 13:07:25
!-----------------------------------------------------------------
!    ########################
     MODULE MODI_DETER_ANGLE
!    ########################
!
INTERFACE
!
     SUBROUTINE DETER_ANGLE(KGEO, KDLON, PULAT, PULON, PVIEW)
!
INTEGER,                INTENT(IN) :: KGEO  
INTEGER,                INTENT(IN) :: KDLON
REAL, DIMENSION(:),     INTENT(IN) :: PULAT
REAL, DIMENSION(:),     INTENT(IN) :: PULON
REAL, DIMENSION(:),     INTENT(OUT):: PVIEW
!
END SUBROUTINE DETER_ANGLE
END INTERFACE
END MODULE MODI_DETER_ANGLE
!     ########################################################
      SUBROUTINE DETER_ANGLE(KGEO, KDLON, PULAT, PULON, PVIEW)
!     ########################################################
!
!**** *DETER_ANGLE* COMPUTES THE VIEWING ANGLE FOR GEOSTATIONARY SATELLITES
!     
!     PURPOSE. 
!     -------- 
!     THIS ROUTINE COMPUTES FOR GIVEN LATITUDES AND LONGITUDES 
!     THE COSECANT OF THE ANGLE UNDER WHICH A GEOSTATIONARY SATELLITE IS SEEN
!     
!     METHOD.
!     -------
!     
!     FORMULA FROM "ISCCP DOCUMENTATION"
!     
!     REFERENCE. 
!     ---------- 
!     
!     W.M.O. - TD 58 ,
!     AND RADIATION PART OF THE MODEL'S DOCUMENTATION
!     
!     EXTERNALS.
!     ---------- 
!     NONE
!     
!     AUTHOR.
!     -------
!     J.-J. Morcrette  *ECMWF*
!     
!     MODIFICATIONS.
!     --------------
!     Original : 88-12-15 (Routine RADGEO)
!     Modified : 97-01-21 (Remy Roca LMD)
!     Modified : 29/03/00 (J-P Chaboureau) rename in DETER_ANGLE for MESONH
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
!
IMPLICIT NONE
!     
!*       0.1  DECLARATIONS OF DUMMY ARGUMENTS
!     
INTEGER,                INTENT(IN) :: KGEO  
INTEGER,                INTENT(IN) :: KDLON
REAL, DIMENSION(:),     INTENT(IN) :: PULAT  ! GEOGRAPHICAL LATITUDES
REAL, DIMENSION(:),     INTENT(IN) :: PULON  ! GEOGRAPHICAL LONGITUDES
REAL, DIMENSION(:),     INTENT(OUT):: PVIEW  ! INVERSE OF COSINE OF VIEWING ANGLE    
!
!*       0.2   LOCAL VARIABLES     
!     
INTEGER :: ISATEL, JL
REAL    :: ZRDSDG, ZALTSA, ZLON, ZCSPSI, ZS2, ZS, ZMUV
REAL, DIMENSION(6) ::  ZGALT  ! ALTITUDE ABOVE EARTH'S SURFACE
                              ! OF GEOSTATIONARY SATELLITE
REAL, DIMENSION(6) ::  ZGNAD  ! LONGITUDE OF NADIR OF GEOSTATIONARY SATELLITE
!----------------------------------------------------------------------------
!     
!*       1.   INITIALIZATION
!              -------------
!     
ZRDSDG = XPI/180.             ! Degree to radian conversion factor
!
!----------------------------------------------------------------------------
!
!*       2.1   GOES EAST SATELLITE
!              -------------------
ISATEL = 1
ZGALT(ISATEL) = 35793000.
ZGNAD(ISATEL) = 2*XPI - 75. * ZRDSDG
!
!*       2.2   GOES WEST SATELLITE
!              -------------------
ISATEL = 2
ZGALT(ISATEL) = 0.
ZGNAD(ISATEL) = 0.
!
!*       2.3   G.M.S. SATELLITE
!              ----------------
ISATEL = 3
ZGALT(ISATEL) = 0.
ZGNAD(ISATEL) = 0.
!
!*       2.4   INDSAT SATELLITE
!              ----------------
ISATEL = 4
ZGALT(ISATEL) = 35793000.
ZGNAD(ISATEL) = 65. * ZRDSDG
!
!*       2.5   METEOSAT SATELLITE
!              ------------------
ISATEL = 5
ZGALT(ISATEL) = 35793000.
ZGNAD(ISATEL) = 0. 
!
!*       2.6   MSG-1 SATELLITE
!              ------------------
ISATEL = 6
ZGALT(ISATEL) = 35793000.
ZGNAD(ISATEL) = 2*XPI - 3.4 * ZRDSDG
!
!----------------------------------------------------------------------------
!
PVIEW(:) = XUNDEF
DO JL=1,KDLON
  ZALTSA=ZGALT(KGEO)
  ZLON  =ZGNAD(KGEO)-PULON(JL)*ZRDSDG
  ZCSPSI = COS (PULAT(JL)*ZRDSDG) * COS (ZLON)
  ZS2 = XRADIUS*XRADIUS + (XRADIUS+ZALTSA)*(XRADIUS+ZALTSA) &
                      - 2.*XRADIUS*(XRADIUS+ZALTSA)*ZCSPSI
  ZS = SQRT ( ZS2 )
  ZMUV = (XRADIUS+ZALTSA)/ZS * ZCSPSI - XRADIUS/ZS
  IF (ZMUV > 0.) PVIEW(JL)=1./ZMUV
END DO
!
END SUBROUTINE DETER_ANGLE
