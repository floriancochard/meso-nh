!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_SPAWN
!     #################
!
!!****  *MODD_SPAWN* - declaration of configuration variables between model 2
!!                          and model 1 for spawning purpose
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the variables
!     which concern the configuration between model 2 and model 1 during spawning
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!       
!!    AUTHOR
!!    ------
!!	P. Jabouille   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/07/99
!!      Modification 08/04/04   (G.Jaubert) Spawning 1 option
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE
!
IMPLICIT NONE
!
INTEGER,SAVE   :: NDXRATIO,NDYRATIO ! x and y-direction resolution RATIO between models 2 and 1
!
INTEGER,SAVE   :: NXSIZE,NYSIZE     ! number of model 1 grid points in x and y-directions
                                    ! in the model 2 physical domain
INTEGER,SAVE   :: NXOR, NYOR        ! horizontal position (i,j) of the
INTEGER,SAVE   :: NXEND,NYEND       ! ORigin and END of model 2 relative to model 1 
!
CHARACTER (LEN=28) :: CDOMAIN       ! input fm-file for grid definition
!
LOGICAL            :: LBAL_ONLY     ! logical switch for spawning 1 with
                                    ! balance calculation only and
                                    ! no modification of the model definition
                                    ! v.s its DAD except the DAD name
!
CHARACTER (LEN=28) :: CDADINIFILE ! DAD fm-file for initial file
                                    ! if LBAL_ONLY=T
CHARACTER (LEN=28) :: CDADSPAFILE ! DAD fm-file for spawning file
                                    ! if LBAL_ONLY=T
!
! Pointers to the data of the mesh #1
REAL,DIMENSION(:),    SAVE,POINTER :: XXHAT1  => NULL()
REAL,DIMENSION(:),    SAVE,POINTER :: XYHAT1  => NULL()
REAL,DIMENSION(:),    SAVE,POINTER :: XZHAT1  => NULL()
REAL,                 SAVE,POINTER :: XZTOP1  => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XZS1    => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XZSMT1  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XZZ1    => NULL()
LOGICAL,              SAVE,POINTER :: LSLEVE1 => NULL()
REAL,                 SAVE,POINTER :: XLEN11  => NULL()
REAL,                 SAVE,POINTER :: XLEN21  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XTHT1   => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XTKET1  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XPABST1 => NULL()
REAL,DIMENSION(:,:,:,:),SAVE,POINTER :: XRT1  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XUT1    => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XVT1    => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XWT1    => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XZWS1   => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XSRCT1  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XSIGS1  => NULL()
TYPE(DATE_TIME),      SAVE,POINTER :: TDTCUR1 => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XLSUM1  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XLSVM1  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XLSWM1  => NULL()
REAL,DIMENSION(:,:)  ,SAVE,POINTER :: XLSZWSM1=> NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XLSTHM1 => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XLSRVM1 => NULL()
!
REAL,DIMENSION(:,:),  SAVE,POINTER :: XINPRC1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XACPRC1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XINDEP1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XACDEP1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XINPRR1   => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XINPRR3D1 => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XEVAP3D1  => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XACPRR1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XINPRS1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XACPRS1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XINPRG1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XACPRG1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XINPRH1   => NULL()
REAL,DIMENSION(:,:),  SAVE,POINTER :: XACPRH1   => NULL()
!
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XVU_FLUX_M1  => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XVTH_FLUX_M1 => NULL()
REAL,DIMENSION(:,:,:),SAVE,POINTER :: XWTH_FLUX_M1 => NULL()
!
END MODULE MODD_SPAWN
