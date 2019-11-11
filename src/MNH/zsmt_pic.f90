!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_ZSMT_PIC
!     ######################
!
INTERFACE 
!
      SUBROUTINE ZSMT_PIC(KSLEVE,PSMOOTH_ZS)
!
INTEGER,             INTENT(IN)  :: KSLEVE     ! number of iterations
REAL,                INTENT(IN)  :: PSMOOTH_ZS ! optional uniform smooth orography for SLEVE coordinate
!
END SUBROUTINE ZSMT_PIC
!
END INTERFACE
!
END MODULE MODI_ZSMT_PIC
!
!
!
!     #############################
      SUBROUTINE ZSMT_PIC(KSLEVE,PSMOOTH_ZS)
!     #############################
!
!!****  *ZSMT_PIC* computes smoothed orography for SLEVE coordinate
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
!!	G. Zangler      * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        nov 2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_LUNIT,      ONLY : CLUOUT0
USE MODD_PARAMETERS, ONLY : XUNDEF
USE MODD_GRID_n,     ONLY : XZS,XZSMT
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,             INTENT(IN)  :: KSLEVE     ! number of iterations
REAL,                INTENT(IN)  :: PSMOOTH_ZS ! optional uniform smooth orography for SLEVE coordinate
!
!
!*       0.2   declarations of local variables
!
!
INTEGER :: JN         ! loop counter on iterations
INTEGER :: JI         ! loop counter on X coordinate
INTEGER :: JJ         ! loop counter on Y coordinate
!
INTEGER :: IIU        ! number of points in X direction
INTEGER :: IJU        ! number of points in Y direction
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
INTEGER :: IINFO_ll                     ! error report of parallel routines
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZS! orography
!-------------------------------------------------------------------------------
IF (PSMOOTH_ZS /= XUNDEF) THEN
!-------------------------------------------------------------------------------
!
!*       1.    Case of uniform smooth orography prescribed
!              -------------------------------------------
!
  XZSMT = PSMOOTH_ZS
!
!-------------------------------------------------------------------------------
ELSE
!-------------------------------------------------------------------------------
!
!*       3.    Computes smoothed orography
!              ---------------------------
!
  NULLIFY(TZFIELDS_ll)
  CALL GET_DIM_EXT_ll ('B',IIU,IJU)
  ALLOCATE(ZZS(IIU,IJU))
!
  XZSMT = XZS
  ZZS   = XZS
!
! ZZS is here used as storage field

  DO JN = 1,KSLEVE
  !
    DO JI = 2,IIU-1
      ZZS(JI,:) = 0.5*XZSMT(JI,:)+0.25*(XZSMT(JI+1,:)+XZSMT(JI-1,:))
    ENDDO
    ZZS(1,:) = XZSMT(1,:)
    ZZS(IIU,:) = XZSMT(IIU,:)
    DO JJ = 2,IJU-1
      XZSMT(:,JJ) = 0.5*ZZS(:,JJ)+0.25*(ZZS(:,JJ+1)+ZZS(:,JJ-1))
    ENDDO
    XZSMT(:,1) = ZZS(:,1)
    XZSMT(:,IJU) = ZZS(:,IJU)
    CALL ADD2DFIELD_ll(TZFIELDS_ll,ZZS)
    CALL ADD2DFIELD_ll(TZFIELDS_ll,XZSMT)
    CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
    CALL CLEANLIST_ll(TZFIELDS_ll)
  ENDDO
!
!-------------------------------------------------------------------------------
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE ZSMT_PIC
