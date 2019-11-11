!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##########################
      MODULE MODD_ICE_C1R3_DESCR
!     ##########################
!
!!****  *MODD_ICE_C1R3_DESCR* - declaration of the microphysical descriptive
!!                              constants for use in the warm and cold schemes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the microphysical
!     constants. This includes the descriptive parameters for the raindrop and
!     the ice crystal habits and the parameters relevant of the dimensional
!     distributions.
!
!         m(D)    = XAx * D**XBx      : Mass-MaxDim relationship
!         v(D)    = XCx * D**XDx      : Fallspeed-MaxDim relationship
!         N(Lbda) = XCCx * Lbda**XCXx : NumberConc-Slopeparam relationship
!         XF0x, XF1x, XF2x            : Ventilation factors
!         XC1x                        : Shape parameter for deposition
!
!              and
!
!         XALPHAx, XNUx                        : Generalized GAMMA law 
!         Lbda = XLBx * (r_x*rho_dref)**XLBEXx : Slope parameter of the 
!                                                distribution law
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_ICE_C1R3_DESCR)
!!          
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24/11/00
!!       J.-P. Pinty   29/11/02 add C1R3NAMES
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
REAL,SAVE :: XCEXVT               ! air density fall speed correction
!
REAL,SAVE :: XAI,XBI,XC_I,XDI         ,XF0I,XF2I,XC1I ! Cloud ice      charact.
REAL,SAVE ::                           XF0IS,XF1IS    ! (large Di vent. coef.)
REAL,SAVE :: XAS,XBS,XCS,XDS,XCCS,XCXS,XF0S,XF1S,XC1S ! Snow/agg.      charact.
REAL,SAVE :: XAG,XBG,XCG,XDG,XCCG,XCXG,XF0G,XF1G,XC1G ! Graupel        charact.
!
REAL,DIMENSION(:),SAVE,ALLOCATABLE :: XRTMIN
                                      ! Min values of the mixing ratios
REAL,DIMENSION(:),SAVE,ALLOCATABLE :: XCTMIN
                                      ! Min values of the drop concentrations
!
REAL,SAVE :: XLBEXI,XLBI              ! Prist. ice     distribution parameters
REAL,SAVE :: XLBEXS,XLBS              ! Snow/agg.      distribution parameters
REAL,SAVE :: XLBEXG,XLBG              ! Graupel        distribution parameters 
!
REAL,SAVE :: XLBDAS_MAX,XLBDAG_MAX    ! Max values allowed for the shape
                                      ! parameters (snow,graupeln)
!
CHARACTER(LEN=10),DIMENSION(2),PARAMETER &
                                   :: C1R3NAMES=(/'CICE  ','CIN   '/)
                                       ! basenames of the SV articles stored
                                       ! in the binary files
!
END MODULE MODD_ICE_C1R3_DESCR 
