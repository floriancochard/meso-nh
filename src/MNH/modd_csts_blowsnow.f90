!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ######################
        MODULE MODD_CSTS_BLOWSNOW
!      ######################
!!
!!     PURPOSE
!!     -------
!!
!!     Declaration of BLOWSNOW constants
!!
!!     METHOD
!!     ------
!!
!!
!!     REFERENCE
!!     ---------
!!     none
!!
!!
!!     AUTHOR
!!     ------
!!     V.Vionnet (GMME)               
!!
!!
!!     MODIFICATIONS
!!     -------------
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
!
USE MODD_CST, ONLY :    &
       XPI              & !Definition of pi
      ,XBOLTZ           & ! Boltzman constant 
      ,XRHOLI           & ![kg/m3] ice density
      ,XG               & ! Gravity constant
      ,XP00             & ! Reference pressure
      ,XMD              & ![kg/mol] molar weight of air
      ,XRD              & ! Gaz constant for dry air
      ,XCPD               !  Cpd (dry air)
!
IMPLICIT NONE
!
 !      Parameters used in Mitchell (96) parameterization for settling velocity  
      REAL,PARAMETER    :: XAM1     = 0.04394
      REAL,PARAMETER    :: XAM2     = 0.06049
      REAL,PARAMETER    :: XAM3     = 0.2072
      REAL,PARAMETER    :: XBM1     = 0.970
      REAL,PARAMETER    :: XBM2     = 0.831
      REAL,PARAMETER    :: XBM3     = 0.638
      REAL,PARAMETER    :: XBESTL_1 = 10.0
      REAL,PARAMETER    :: XBESTL_2 = 585.
 !      Parameters used in Nusselt Number computation 
      REAL,PARAMETER    :: XANU     = 1.88
      REAL,PARAMETER    :: XBNU     = 0.58

END MODULE MODD_CSTS_BLOWSNOW
