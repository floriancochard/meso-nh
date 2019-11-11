!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_NESTING
!     ###################
!
!!****  *MODD_NESTING* - declaration of gridnesting configuration variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the variables
!     which concern the gridnesting configuration of all models.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS  :
!!         JPMODELMAX : Maximum allowed  number of nested models
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_NESTING)
!!
!!    AUTHOR
!!    ------
!!	J.P. Lafore   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/08/95
!!      updated     29/07/96  (J.P. Lafore) MY_NAME(m) introduction
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 08/02/2019: add missing NULL association for pointers
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
                            ! resolution RATIO between models m and its father NDAD(m)
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDXRATIO_ALL        ! in x-direction 
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDYRATIO_ALL        ! in y-direction 
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDTRATIO            ! in Time 
!
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NXOR_ALL, NYOR_ALL  ! horizontal position (i,j) of the
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NXEND_ALL,NYEND_ALL ! ORigin and END of model m 
                                                     ! relative to its father NDAD(m)    
!
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDAD ! model number of the father of each model "m"
REAL,SAVE,     DIMENSION(JPMODELMAX) :: XWAY ! model m interactive nesting level with its father NDAD(m)
!
                                                            !   MeSsaGes concerning 
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_IF  ! var. Interpolation at Flux
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_IS  ! and Scalar location
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_AVR ! AVeRage
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_END ! timestep END
                                                            !   MeSsaGes concerning
INTEGER,SAVE,  DIMENSION(JPMODELMAX,JPMODELMAX) :: NMSG_AVR_END ! AVeRage END
!
CHARACTER(LEN=NFILENAMELGTMAX),SAVE,   DIMENSION(JPMODELMAX) :: CMY_NAME,CDAD_NAME
                                                  ! names of the initial FM-Files
                                                  ! then generic names of output FM-Files
                                                  ! of each model "m"
                                                  ! and of its DAD model
                                                  ! (read and written on the LFI parts)
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NDT_2_WAY ! number of times the time step
              ! of model n used for the relaxation time of the 2_WAY grid-nesting
              ! interaction  i.e. Tau = NDT_2_WAY * XTSTEP


INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NIMAX_NEST, NJMAX_NEST  ! local sizes of model m
INTEGER,SAVE,  DIMENSION(JPMODELMAX) :: NIMAX_NEST_ll, NJMAX_NEST_ll  ! globcal sizes of model m
LOGICAL,SAVE,  DIMENSION(JPMODELMAX) :: L1D_NEST         ! Logical for 1D model version of model m
LOGICAL,SAVE,  DIMENSION(JPMODELMAX) :: L2D_NEST         ! Logical for 2D model version of model m
LOGICAL,SAVE,  DIMENSION(JPMODELMAX) :: LPACK_NEST       ! Logical to compress 1D or 2D FM files of model m
!
TYPE REAL_FIELD2D_ALL
    REAL, DIMENSION(:,:), POINTER :: XFIELD2D => NULL()
END TYPE REAL_FIELD2D_ALL

TYPE REAL_FIELD1D_ALL
    REAL, DIMENSION(:), POINTER :: XFIELD1D => NULL()
END TYPE REAL_FIELD1D_ALL
!
TYPE(REAL_FIELD2D_ALL), DIMENSION(JPMODELMAX), TARGET :: TXZS   ! orography of model m
TYPE(REAL_FIELD2D_ALL), DIMENSION(JPMODELMAX), TARGET :: TXZSMT   ! smooth orography for SLEVE coordinate of model m
TYPE(REAL_FIELD1D_ALL), DIMENSION(JPMODELMAX), TARGET :: TXXHAT   ! Position x in the
                                         ! conformal or cartesian plane of model m
TYPE(REAL_FIELD1D_ALL), DIMENSION(JPMODELMAX), TARGET :: TXYHAT   ! Position y in the
                                         ! conformal or cartesian plane of model m



END MODULE MODD_NESTING
