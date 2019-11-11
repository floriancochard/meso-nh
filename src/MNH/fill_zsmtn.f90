!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_FILL_ZSMTn
!     ######################
!
INTERFACE 
!
      SUBROUTINE FILL_ZSMT_n(HFIELD,PFIELD,KSON)
!
CHARACTER(LEN=6),         INTENT(IN) :: HFIELD ! name of the field to nest
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PFIELD
INTEGER,                  INTENT(IN) :: KSON   ! son model index
!
END SUBROUTINE FILL_ZSMT_n
!
END INTERFACE
!
END MODULE MODI_FILL_ZSMTn
!
!
!
!     ##########################################
      SUBROUTINE FILL_ZSMT_n(HFIELD,PFIELD,KSON)
!     ##########################################
!
!!****  *FILL_ZSMT_n* - fill the working array for nesting of pgd files
!!                          with KSON model index
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        12/01/05
!!      Modification    20/05/06 Remove Clark and Farley interpolation
!!        M.Moge        01/2016 bug fix for parallel execution
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_GRID_n,  ONLY : XZSMT
USE MODD_LBC_n,   ONLY : CLBCX,CLBCY
USE MODD_NESTING
USE MODD_PARAMETERS
!
USE MODI_INI_BIKHARDT_n
USE MODI_SPAWN_ZS
USE MODE_MODELN_HANDLER
!
USE MODE_SPLITTING_ll, ONLY : SPLIT2, DEF_SPLITTING2
USE MODD_VAR_ll, ONLY : NPROC, IP, YSPLITTING, NMNH_COMM_WORLD
USE MODD_STRUCTURE_ll, ONLY : ZONE_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
CHARACTER(LEN=6),     INTENT(IN)     :: HFIELD ! name of the field to fill
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PFIELD
INTEGER,              INTENT(IN) :: KSON   ! son model index
!
!*       0.2   declarations of local variables
!-------------------------------------------------------------------------------
INTEGER :: IMI ! current model index (DAD index)
!
! Dummy pointers needed to correct an ifort Bug
CHARACTER(LEN=4), DIMENSION(:), POINTER :: DPTR_CLBCX,DPTR_CLBCY
REAL, DIMENSION(:,:),  POINTER          :: DPTR_XZSMT
INTEGER :: IINFO_ll
INTEGER :: IXSIZE, IYSIZE        ! sizes of global son domain in father grid
INTEGER :: IXSIZE_F, IYSIZE_F    ! sizes of global father domain
TYPE(ZONE_ll), DIMENSION(:), ALLOCATABLE :: TZSPLITTING
INTEGER :: IXOR,IXEND,IYOR,IYEND ! limits of extended  domain of KSON model in its father's grid
INTEGER :: IDIMX, IDIMY  ! dimensions of extended son subdomain in father's grid + one point in each direction
INTEGER :: IXDOMAINS, IYDOMAINS               ! number of subdomains in X and Y directions
LOGICAL :: GPREM                              ! needed for DEF_SPLITTING2, true if NPROC is a prime number
!
!*       1.    initializations
!              ---------------
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL GOTO_MODEL(KSON)
CALL GO_TOMODEL_ll(KSON, IINFO_ll)
!
! get sizes of global son domain in father grid
IXSIZE = NXEND_ALL(KSON) - NXOR_ALL (KSON) + 1 - 2*JPHEXT
IYSIZE = NYEND_ALL(KSON) - NYOR_ALL (KSON) + 1 - 2*JPHEXT
! get splitting of current model KMI in father grid
IXSIZE_F = NXEND_ALL(NDAD(KSON)) - NXOR_ALL (NDAD(KSON)) + 1 - 2*JPHEXT
IYSIZE_F = NYEND_ALL(NDAD(KSON)) - NYOR_ALL (NDAD(KSON)) + 1 - 2*JPHEXT
ALLOCATE(TZSPLITTING(NPROC))
! we want the same domain partitioning for the child domain and for the father domain
CALL DEF_SPLITTING2(IXDOMAINS,IYDOMAINS,IXSIZE_F,IYSIZE_F,NPROC,GPREM)
CALL SPLIT2 ( IXSIZE, IYSIZE, 1, NPROC, TZSPLITTING, YSPLITTING, IXDOMAINS, IYDOMAINS )
! get coords of extended domain of KSON in its father's grid
IXOR  = NXOR_ALL(KSON) + TZSPLITTING(IP)%NXOR  -1 - JPHEXT 
IXEND = NXOR_ALL(KSON) + TZSPLITTING(IP)%NXEND -1 + JPHEXT 
IYOR  = NYOR_ALL(KSON) + TZSPLITTING(IP)%NYOR  -1 - JPHEXT 
IYEND = NYOR_ALL(KSON) + TZSPLITTING(IP)%NYEND -1 + JPHEXT
!
!IDIMX = IXEND - IXOR - 1
!IDIMY = IYEND - IYOR - 1
IDIMX = IXEND - IXOR + 1 +2*1 ! + 2*JPHEXT
IDIMY = IYEND - IYOR + 1 +2*1 ! + 2*JPHEXT
!
CALL INI_BIKHARDT_n(NDXRATIO_ALL(KSON),NDYRATIO_ALL(KSON),KSON)
!
!-------------------------------------------------------------------------------
!
!*       2.    interpolation of dad field
!              --------------------------
!
DPTR_CLBCX=>CLBCX
DPTR_CLBCY=>CLBCY
DPTR_XZSMT=>XZSMT
CALL SPAWN_ZS(NXOR_ALL(KSON),NXEND_ALL(KSON),NYOR_ALL(KSON),NYEND_ALL(KSON), &
              NDXRATIO_ALL(KSON),NDYRATIO_ALL(KSON),IDIMX,IDIMY,DPTR_CLBCX,DPTR_CLBCY,         &
              PFIELD,DPTR_XZSMT,HFIELD                             )
!-------------------------------------------------------------------------------
!
CALL GOTO_MODEL(IMI)
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
!
END SUBROUTINE FILL_ZSMT_n
