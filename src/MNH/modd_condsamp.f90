!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_CONDSAMP
!     ##################
!-------------------------------------------------------------------------------
!***	MODD_CONDSAMP  Declaration of conditional sampling tracers 
!
!!    AUTHOR
!!    ------
!	           : C.Lac                               
!	Creation   : 01/06/2011
!
!-------------------------------------------------------------------------------
!
!
!*    0. DECLARATIONS
!        ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
LOGICAL            :: LCONDSAMP = .FALSE.  ! Switch to activate conditional sampling
!
INTEGER, PARAMETER :: JPCSMAX = 3     
!
INTEGER                            :: NCONDSAMP            ! Number of conditional
                                                           ! sampling tracers
REAL,         DIMENSION(JPCSMAX)   :: XRADIO               ! Radioactive decay period
REAL,         DIMENSION(JPCSMAX)   :: XSCAL                ! Scaling factor
REAL                               :: XHEIGHT_BASE         ! Distance below the
                              !         cloud base where the 2nd tracer is emitted
REAL                               :: XDEPTH_BASE          ! Depth in which the
                              !         2nd tracer is emitted
REAL                               :: XHEIGHT_TOP          ! Distance above the
                              !         cloud top  where the 3rd tracer is emitted
REAL                               :: XDEPTH_TOP           ! Depth in which the
                              !         3rd tracer is emitted
!
END MODULE MODD_CONDSAMP
