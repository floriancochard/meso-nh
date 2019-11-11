!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_UPDATE_METRICS
!     ###################
INTERFACE
!
SUBROUTINE UPDATE_METRICS(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ)
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X boundary type
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y boundary type
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDXX  ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDYY  ! metric coefficient dyy 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZX  ! metric coefficient dzx 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZY  ! metric coefficient dzy 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZZ  ! metric coefficient dzz  
!
END SUBROUTINE UPDATE_METRICS
!
END INTERFACE
!
END MODULE MODI_UPDATE_METRICS
!
!
!
!     #################################################################
      SUBROUTINE UPDATE_METRICS(HLBCX,HLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ)
!     #################################################################
!
!!****  *UPDATE_METRICS* - routine to set external points for metric coefficients
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
!!      Book2 of documentation (routine UPDATE_METRICS)
!!
!!    AUTHOR
!!    ------
!!	V. Masson        * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    april 2006
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!         
USE MODD_CONF
USE MODD_PARAMETERS
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! X boundary type
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! Y boundary type
!
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDXX  ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDYY  ! metric coefficient dyy 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZX  ! metric coefficient dzx 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZY  ! metric coefficient dzy 
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PDZZ  ! metric coefficient dzz  
!
!*       0.2   declarations of local variables
!
INTEGER             :: IIB      ! First physical index in x direction
INTEGER             :: IJB      ! First physical index in y direction
INTEGER             :: JI       ! loop index
!
TYPE(LIST_ll), POINTER :: TZMETRICS_ll   ! list of fields to exchange
INTEGER                :: IINFO_ll       ! return code of parallel routine
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS :
!              ----------------------------
IIB = 1 + JPHEXT
IJB = 1 + JPHEXT
!
NULLIFY(TZMETRICS_ll)
!
!-------------------------------------------------------------------------------
!
!*       3.  UPDATE EXTERNAL POINTS OF GLOBAL DOMAIN:
!            ---------------------------------------
!
IF ( HLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
  PDXX(IIB-1,:,:) = PDXX(IIB,:,:)
  PDZX(IIB-1,:,:) = PDZX(IIB,:,:)
END IF
IF ( HLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
  DO JI=1,SIZE(PDYY,1)
    PDYY(JI,IJB-1,:) = PDYY(JI,IJB,:)
    PDZY(JI,IJB-1,:) = PDZY(JI,IJB,:)
  END DO
END IF

!-------------------------------------------------------------------------------
!
!*       2.  UPDATE HALOs :
!            -------------
!
!
!!$IF(NHALO == 1) THEN
  CALL ADD3DFIELD_ll(TZMETRICS_ll,PDXX)
  CALL ADD3DFIELD_ll(TZMETRICS_ll,PDYY)
  CALL ADD3DFIELD_ll(TZMETRICS_ll,PDZX)
  CALL ADD3DFIELD_ll(TZMETRICS_ll,PDZY)
  CALL ADD3DFIELD_ll(TZMETRICS_ll,PDZZ)
  CALL UPDATE_HALO_ll(TZMETRICS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZMETRICS_ll)
!!$END IF

!-----------------------------------------------------------------------------
END SUBROUTINE UPDATE_METRICS
