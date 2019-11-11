!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:./s.modn_outfile.f90, Version:1.3, Date:03/06/05, Last modified:01/10/19
!-----------------------------------------------------------------
!     ####################
      MODULE  MODN_OUTFILE
!     ####################
!
!!****  *MODN_OUTFILE* - defines the three namelists controling conversion
!!
!!    PURPOSE
!!    -------
!      This declarative module defines the NAM_OUTFILE, NAM_OUTHOR, NAM_OUTVER
!     namelists, which contains the parameters controling the grib or Vis5D
!     coding of fields.
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Ducrocq   
!!
!!    MODIFICATIONS
!!    -------------
!!     original        20/03/97
!!     modifications   20/02/01 (I.Mallet) merge with JPChaboureau Vis5D files
!!     modifications   20/10/01 ( " ) split in 3 namelists, add horizontal interpolation
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!
USE MODD_CONF, ONLY : NVERB
IMPLICIT NONE
!
CHARACTER(LEN=28) :: CMNHFILE  ! Name of the input FM file
CHARACTER(LEN=3)  :: COUTFILETYPE  ! Type of the outfile (GRB or V5D)
! Common characteristics
LOGICAL::  LAGRID      !If  T., fields are interpolated on an arakawa A-grid
                       ! (mass grid) else they are on the mesonh grids
CHARACTER(LEN=4)  :: CHORTYPE  ! Type of horizontal grid
      ! NONE: MesoNH grid
      ! NEAR: nearest-neighbour interpolation
      ! BILI: bilinear interpolation
REAL,DIMENSION(4) :: XLATLON   ! NSWE target domain bounds (in degrees)
CHARACTER(LEN=1)  :: CLEVTYPE  ! Type of vertical levels in output file
      ! GRB: P=pressure levels, K=native coordinate of MESO-NH
      ! V5D: P=pressure levels, Z=z levels, K=native coordinate of lowest point
CHARACTER(LEN=6)  :: CLEVLIST  ! How vertical levels are specified
      ! 'MANUAL' list of levels in free format
      ! 'FUNCTN' list of levels in next 3 variables
REAL :: XVLMIN,XVLMAX,XVLINT ! minimum, maximum and increment values
                           ! for the vertical grid
                           ! (used only if CLEVTYPE='P' or 'Z')
! Grib characteristics
LOGICAL ::  LLMULTI   !If .T., a multigrib file is produced, else monogrib files
!
!*     0.1  Namelist NAM_OUTFILE
!
NAMELIST/NAM_OUTFILE/CMNHFILE,COUTFILETYPE, &
                     NVERB,LLMULTI
NAMELIST/NAM_OUTHOR/LAGRID,CHORTYPE,XLATLON
NAMELIST/NAM_OUTVER/CLEVTYPE,CLEVLIST, &
                    XVLMIN,XVLMAX,XVLINT
!
END MODULE MODN_OUTFILE
