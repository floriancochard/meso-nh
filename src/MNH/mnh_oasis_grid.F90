!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##########
MODULE MODI_MNH_OASIS_GRID
!     ##########
!
INTERFACE 
!
      SUBROUTINE MNH_OASIS_GRID(OD_MASTER,KD_LCOMM)
!
      LOGICAL, INTENT(IN) :: OD_MASTER                  ! MASTER process or not
      INTEGER, INTENT(IN) :: KD_LCOMM                   ! Model local communicator
!
      END SUBROUTINE MNH_OASIS_GRID
!
END INTERFACE
!
END MODULE MODI_MNH_OASIS_GRID
!
!     ####################################################################
SUBROUTINE MNH_OASIS_GRID(OD_MASTER,KD_LCOMM)
!     ####################################################################
!
!
!!****  *MNH_OASIS_GRID*
!!
!!    PURPOSE
!!    -------
!!    Define the grids for OASIS coupling
!!    The grids definition for the hydrological coupling part has to be coded
!!    (cf. sfx_oasis_prep.F90).
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	J. Pianezze   *LOPS*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2016
!-------------------------------------------------------------------------------
!
!*     0.     Declarations
!             ------------
!
!
#ifdef CPLOASIS
USE MOD_OASIS
#endif
!
USE MODD_PARAMETERS,  ONLY : XUNDEF
USE MODD_IO_SURF_MNH, ONLY : NHALO
USE MODD_CST,         ONLY : XPI, XRADIUS
USE MODD_DIM_n,       ONLY : NIMAX_ll, NJMAX_ll, NIMAX, NJMAX
USE MODD_PARAMETERS,  ONLY : JPHEXT
USE MODD_GRID_n,      ONLY : XLAT, XLON 
USE MODD_VAR_ll,      ONLY : NMNH_COMM_WORLD
USE MODD_SFX_OASIS
USE MODD_MNH_SURFEX_n
USE MODD_MPIF
!
USE MODI_GET_FRAC_N
!
USE MODE_GATHER_ll
USE MODE_ll
!
!*      0.1   Declarations of argument
! ------------------------------------------------
LOGICAL, INTENT(IN) :: OD_MASTER                  ! MASTER process or not
INTEGER, INTENT(IN) :: KD_LCOMM                   ! Model local communicator
!
!*      0.2  Declaration of local parameters
! ------------------------------------------------
!
INTEGER,           PARAMETER  :: INC = 4    ! Number of grid-cell corners
!
CHARACTER(LEN=4),  PARAMETER  :: YSFX_LAND = 'slan'
CHARACTER(LEN=4),  PARAMETER  :: YSFX_LAKE = 'slak'
CHARACTER(LEN=4),  PARAMETER  :: YSFX_SEA  = 'ssea'
!
!*      0.2  Declaration of local variables
! ------------------------------------------------
INTEGER :: JI, JJ ! loop index
INTEGER :: IIU, IJU, ILU
INTEGER :: IIU_ll, IJU_ll
INTEGER :: IERROR
!
REAL :: ZPHI, ZRADEG
REAL :: ZXBOX, ZYBOX
!
REAL, DIMENSION(NIMAX_ll,NJMAX_ll,INC) :: ZCLON, ZCLAT
REAL, DIMENSION(NIMAX_ll,NJMAX_ll)    :: ZAREA
REAL, DIMENSION(NIMAX_ll,NJMAX_ll)    :: ZDX, ZDY
REAL, ALLOCATABLE, DIMENSION(:,:)    :: ZSEA, ZWATER
REAL, ALLOCATABLE, DIMENSION(:,:)    :: ZTOWN, ZNATURE
!
REAL, ALLOCATABLE, DIMENSION(:)    :: ZWATER1D
REAL, ALLOCATABLE, DIMENSION(:)    :: ZNATURE1D
REAL, ALLOCATABLE, DIMENSION(:)    :: ZTOWN1D
REAL, ALLOCATABLE, DIMENSION(:)    :: ZSEA1D
!
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZLAT_GLOBAL, ZLON_GLOBAL
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZSEA_GLOBAL, ZWATER_GLOBAL
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZTOWN_GLOBAL, ZNATURE_GLOBAL
!
INTEGER, DIMENSION(NIMAX_ll,NJMAX_ll) :: ZMASK_LAND
INTEGER, DIMENSION(NIMAX_ll,NJMAX_ll) :: ZMASK_LAKE
INTEGER, DIMENSION(NIMAX_ll,NJMAX_ll) :: ZMASK_SEA
!
!-------------------------------------------------------------------------------
!
!*     1.     Initialize :
!             ------------
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
IIU_ll=NIMAX_ll+2*JPHEXT
IJU_ll=NJMAX_ll+2*JPHEXT
ILU = (IIE-IIB+1+2*NHALO)*(IJE-IJB+1+2*NHALO)
!
ALLOCATE(ZLAT_GLOBAL(IIU_ll,IJU_ll,1))
ALLOCATE(ZLON_GLOBAL(IIU_ll,IJU_ll,1))
!
ALLOCATE(ZSEA_GLOBAL   (IIU_ll,IJU_ll,1))
ALLOCATE(ZWATER_GLOBAL (IIU_ll,IJU_ll,1))
ALLOCATE(ZNATURE_GLOBAL(IIU_ll,IJU_ll,1))
ALLOCATE(ZTOWN_GLOBAL  (IIU_ll,IJU_ll,1))
!
ALLOCATE(ZSEA    (IIU,IJU))
ALLOCATE(ZWATER  (IIU,IJU))
ALLOCATE(ZNATURE (IIU,IJU))
ALLOCATE(ZTOWN   (IIU,IJU))
!
ALLOCATE(ZSEA1D    ( ILU ))
ALLOCATE(ZWATER1D  ( ILU ))
ALLOCATE(ZNATURE1D ( ILU ))
ALLOCATE(ZTOWN1D   ( ILU ))
!
!-------------------------------------------------------------------------------
!
!*     2.     Get grid definition :
!             ---------------------
!
!
!*     2.1    Get lat/lon :
!             -------------
CALL GATHER_XYFIELD(XLON,ZLON_GLOBAL(:,:,1),1,NMNH_COMM_WORLD)
CALL GATHER_XYFIELD(XLAT,ZLAT_GLOBAL(:,:,1),1,NMNH_COMM_WORLD)
!
!*     2.2    Get corners :
!             -------------
ZXBOX=ZLON_GLOBAL(2,2,1)-ZLON_GLOBAL(2,1,1)
ZCLON(:,:,1)=ZLON_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)+ZXBOX/2.
ZCLON(:,:,2)=ZLON_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)-ZXBOX/2.
ZCLON(:,:,3)=ZLON_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)-ZXBOX/2.
ZCLON(:,:,4)=ZLON_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)+ZXBOX/2.
!
ZYBOX=ZLAT_GLOBAL(2,2,1)-ZLAT_GLOBAL(1,2,1)
ZCLAT(:,:,1)=ZLAT_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)+ZYBOX/2.
ZCLAT(:,:,2)=ZLAT_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)+ZYBOX/2.
ZCLAT(:,:,3)=ZLAT_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)-ZYBOX/2.
ZCLAT(:,:,4)=ZLAT_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)-ZYBOX/2.
!
!
!*     2.3    Compute area :
!             -------------
! Radian to degree
ZRADEG = XPI/180.0
!
DO JJ=1,NJMAX_ll
  DO JI =1,NIMAX_ll
    ZPHI = ZLAT_GLOBAL(JI+1,JJ+1,1)
!
    ZDX(JI,JJ)=XRADIUS*COS(ZPHI*ZRADEG)*(ZLON_GLOBAL(JI+1,JJ+1,1)-ZLON_GLOBAL(JI,JJ+1,1))*ZRADEG
    ZDY(JI,JJ)=XRADIUS*(ZLAT_GLOBAL(JI+1,JJ+1,1)-ZLAT_GLOBAL(JI+1,JJ,1))*ZRADEG
!
  ENDDO
ENDDO
!
ZAREA(:,:)=ZDY(:,:)*ZDX(:,:)
!
!-------------------------------------------------------------------------------
!
!*     3.     Compute masks
!             -------------
!
CALL GET_FRAC_n(YSURF_CUR%U,'MESONH',ILU,ZSEA1D,ZWATER1D,ZNATURE1D,ZTOWN1D)
CALL REMOVE_HALO(ZSEA1D,ZSEA)
CALL REMOVE_HALO(ZWATER1D,ZWATER)
CALL REMOVE_HALO(ZNATURE1D,ZNATURE)
CALL REMOVE_HALO(ZTOWN1D,ZTOWN)
!
CALL GATHER_XYFIELD(ZSEA,ZSEA_GLOBAL(:,:,1),1,NMNH_COMM_WORLD)
CALL GATHER_XYFIELD(ZWATER,ZWATER_GLOBAL(:,:,1),1,NMNH_COMM_WORLD)
CALL GATHER_XYFIELD(ZNATURE,ZNATURE_GLOBAL(:,:,1),1,NMNH_COMM_WORLD)
CALL GATHER_XYFIELD(ZTOWN,ZTOWN_GLOBAL(:,:,1),1,NMNH_COMM_WORLD)
!
!*       3.1    Mask for Land surface :
!               ----------------------- 
ZMASK_LAND(:,:)=1
DO JJ=1,NJMAX_ll
  DO JI=1,NIMAX_ll
    IF ( (ZNATURE_GLOBAL(JI+1,JJ+1,1)+ZTOWN_GLOBAL(JI+1,JJ+1,1)) /= 0 ) ZMASK_LAND(JI,JJ)=0
  ENDDO
ENDDO
!
!*       3.2    Mask for Lake surface :
!               ----------------------- 
ZMASK_LAKE(:,:)=1
DO JJ=1,NJMAX_ll
  DO JI=1,NIMAX_ll
    IF ( ZWATER_GLOBAL(JI+1,JJ+1,1) /= 0 ) ZMASK_LAKE(JI,JJ)=0
  ENDDO
ENDDO
!
!*       3.3    Mask for sea/water/wave surface :
!               ---------------------------------
ZMASK_SEA(:,:)=1
DO JJ=1,NJMAX_ll
  DO JI=1,NIMAX_ll
    IF ( ZSEA_GLOBAL(JI+1,JJ+1,1) /= 0 ) ZMASK_SEA(JI,JJ)=0
  ENDDO ! XSEA=0: land, XSEA=1: sea, XSEA=2: sea-ice
ENDDO
!
CALL MPI_BARRIER(KD_LCOMM, IERROR)
!
!-------------------------------------------------------------------------------
!
!*    4.     Write grid definition :
!            -----------------------
!
!
IF (OD_MASTER) THEN
!
  CALL OASIS_START_GRIDS_WRITING(IERROR)
!
!*       4.1    Grid definition for Land surface :
!               ----------------------------------
!
  IF (LCPL_LAND) THEN
    CALL OASIS_WRITE_GRID (YSFX_LAND, NIMAX_ll, NJMAX_ll,             &
                           ZLON_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1), &
                           ZLAT_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)  )
    CALL OASIS_WRITE_CORNER(YSFX_LAND, NIMAX_ll, NJMAX_ll, INC, ZCLON, ZCLAT)   
    CALL OASIS_WRITE_AREA(YSFX_LAND, NIMAX_ll, NJMAX_ll, ZAREA)
    CALL OASIS_WRITE_MASK(YSFX_LAND, NIMAX_ll, NJMAX_ll, ZMASK_LAND)
  ENDIF
!
!*       4.2    Grid definition for lake surface :
!               ----------------------------------
!
  IF (LCPL_LAKE) THEN
    CALL OASIS_WRITE_GRID (YSFX_LAKE, NIMAX_ll, NJMAX_ll,             &
                           ZLON_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1), &
                           ZLAT_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)  )
    CALL OASIS_WRITE_CORNER(YSFX_LAKE, NIMAX_ll, NJMAX_ll, INC, ZCLON, ZCLAT)   
    CALL OASIS_WRITE_AREA(YSFX_LAKE, NIMAX_ll, NJMAX_ll, ZAREA)
    CALL OASIS_WRITE_MASK(YSFX_LAKE, NIMAX_ll, NJMAX_ll, ZMASK_LAKE)
  ENDIF
!
!
!*       4.3    Grid definition for sea/water :
!               -------------------------------
  IF (LCPL_SEA .OR. LCPL_WAVE) THEN
    CALL OASIS_WRITE_GRID (YSFX_SEA, NIMAX_ll, NJMAX_ll,             &
                           ZLON_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1), &
                           ZLAT_GLOBAL(2:NIMAX_ll+1,2:NJMAX_ll+1,1)  )
    CALL OASIS_WRITE_CORNER(YSFX_SEA, NIMAX_ll, NJMAX_ll, INC, ZCLON, ZCLAT)   
    CALL OASIS_WRITE_AREA(YSFX_SEA, NIMAX_ll, NJMAX_ll, ZAREA)
    CALL OASIS_WRITE_MASK(YSFX_SEA, NIMAX_ll, NJMAX_ll, ZMASK_SEA)
  ENDIF
!
!-------------------------------------------------------------------------------
!
!*    5.     Terminate grid definition :
!            -----------------------
  CALL OASIS_TERMINATE_GRIDS_WRITING()
!
ENDIF
!
DEALLOCATE(ZLAT_GLOBAL)
DEALLOCATE(ZLON_GLOBAL)
!
DEALLOCATE(ZSEA_GLOBAL)
DEALLOCATE(ZWATER_GLOBAL)
DEALLOCATE(ZNATURE_GLOBAL)
DEALLOCATE(ZTOWN_GLOBAL)
!
DEALLOCATE(ZSEA)
DEALLOCATE(ZWATER)
DEALLOCATE(ZNATURE)
DEALLOCATE(ZTOWN)
!
DEALLOCATE(ZSEA1D)
DEALLOCATE(ZWATER1D)
DEALLOCATE(ZNATURE1D)
DEALLOCATE(ZTOWN1D)
!
CALL MPI_BARRIER(KD_LCOMM, IERROR)
!
!==============================================================================
!
CONTAINS
!
SUBROUTINE REMOVE_HALO(PFIELD,POUT)
!
REAL, DIMENSION(:),   INTENT(IN)  :: PFIELD
REAL, DIMENSION(:,:), INTENT(OUT) :: POUT
!
INTEGER :: JI, JJ
!
POUT=XUNDEF
!
DO JJ=IJB,IJE
  DO JI=IIB,IIE
    POUT(JI,JJ) = PFIELD( JI-IIB+1 + NHALO + (JJ-IJB+NHALO)*(IIE-IIB+1+2*NHALO))
  END DO
END DO
!
END SUBROUTINE REMOVE_HALO
!
END SUBROUTINE MNH_OASIS_GRID
!!======================================================================
