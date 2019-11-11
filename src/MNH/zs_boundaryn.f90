!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #######################
      MODULE MODI_ZS_BOUNDARY_n
!     #######################
!
INTERFACE 
!
      SUBROUTINE ZS_BOUNDARY_n (KMI,                                      & 
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,      &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,      &
                    KDXRATIO,KDYRATIO,              &
                    HLBCX,HLBCY,                                          &
                    PZS                                                   )
!
INTEGER,              INTENT(IN)     :: KMI ! Model index
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
!
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PZS ! orography of the fine mesh model
!
END SUBROUTINE ZS_BOUNDARY_n
!
END INTERFACE
!
END MODULE MODI_ZS_BOUNDARY_n
!
!     #####################################################################
      SUBROUTINE ZS_BOUNDARY_n (KMI,                                      & 
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,      &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,      &
                    KDXRATIO,KDYRATIO,              &
                    HLBCX,HLBCY,                                          &
                    PZS                                                   )
!     #####################################################################
!
!!****  *ZS_BOUNDARY_n* - interpolate the orography of the DAD model to the
!!                        CHILD model
!!
!!    PURPOSE
!!    -------
!!      The purpose of ZS_BOUNDARY$n is to perform the Bikhardt interpolation
!!    from the DAD model orography toward the fine scale model KMI. 
!
!
!!**  METHOD
!!    ------
!!      
!!       We use the Bikhardt interpolation scheme to compute the orography of
!!    the fine-mesh model from the value of its DAD orography. This intermediate
!!    smooth orography is only used to modify the orography of the fine-mesh
!!    model for the marginal points ( 1, IIU, 1, IJU). The fine-scale orography
!!    is imported in this subroutine by argument and the DAD orography by module
!!    MODD_GRID$n 
!!
!!    EXTERNAL
!!    --------
!!       BIKHARDT : horizontal interpolation scheme using 4 points in x and y
!!                  directions
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODULE MODD_GRID$n:   XZS
!!
!!    REFERENCE
!!    ---------
!!    
!!
!!    AUTHOR
!!    ------
!!    J. Stein  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     1/2/99 
!!                 
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_GRID_n       ! contains the DAD model informations
!
USE MODI_BIKHARDT
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments 
! 
!
!
INTEGER,              INTENT(IN)     :: KMI ! Model index
REAL, DIMENSION(:), INTENT(IN) :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL, DIMENSION(:), INTENT(IN) :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
INTEGER,   INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,   INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN) :: HLBCY   ! boundary conditions
!
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PZS ! orography of the fine mesh model
!
!*       0.2   declarations of local variables
!
!REAL   , DIMENSION(SIZE(PZS,1),SIZE(PZS,2))  :: ZZSLS ! interpolated orography 
                                                      ! at high resolution 
REAL   , DIMENSION(SIZE(PZS,1),SIZE(PZS,2),1)                         :: ZZSLS
INTEGER    :: IIB,IJB,IIE,IJE
INTEGER    :: IIU,IJU  ! array sizes in x and y directions of the CHILD model 
!
! Variables used for LS communications
TYPE(LIST_ll), POINTER :: TZLSFIELD_ll   ! list of fields to exchange
INTEGER :: IINFO_ll, IDIMX, IDIMY
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTZS,ZZS   !!! provisoire
!REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTZS
INTEGER :: IMI !Current model index
!-------------------------------------------------------------------------------
!
!*      0.   INITIALISATION
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IIB=IIB-1
IIE=IIE+1
IJB=IJB-1
IJE=IJE+1
!
!*      1 GATHER LS FIELD FOR THE CHILD MODEL KMI
!
!       1.1  Must be on the father model to call get_child_dim
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
CALL GET_CHILD_DIM_ll(KMI, IDIMX, IDIMY, IINFO_ll)
!
!       1.2  Allocate array which will receive coarse grid points
!
ALLOCATE(ZTZS(IDIMX,IDIMY,1))
ALLOCATE(ZZS(SIZE(XZS,1),SIZE(XZS,2),1))
ZZS(:,:,1)=XZS(:,:)
!
!         1.3  Specify the ls "source" fields and receiver fields
!
CALL SET_LSFIELD_1WAY_ll(ZZS, ZTZS, KMI)
!
!
!        1.4  Communication
!
CALL LS_FORCING_ll(KMI, IINFO_ll)
!
!        1.5  Back to the (current) child model
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
CALL UNSET_LSFIELD_1WAY_ll()
!
!
!*       1.    BIKARDT INTERPOLATION
!              ---------------------
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,IDIMX-1,IDIMY-1,KDXRATIO,KDYRATIO,1,         &
               HLBCX,HLBCY,ZTZS,ZZSLS(IIB:IIE,IJB:IJE,:))
!
DEALLOCATE(ZTZS,ZZS)
!-------------------------------------------------------------------------------
!
!*       2.    SET THE OROGRAPHY AT THE MARGINAL POINTS
!              ----------------------------------------
!
!
IIU=SIZE(PZS,1)
IJU=SIZE(PZS,2)
!
IF(HLBCX(1)/='CYCL' .AND. LWEST_ll()) PZS(1,IJB:IJE)   = ZZSLS(1,IJB:IJE,1)
IF(HLBCX(2)/='CYCL' .AND. LEAST_ll()) PZS(IIU,IJB:IJE) = ZZSLS(IIU,IJB:IJE,1)
IF(HLBCY(1)/='CYCL' .AND. LSOUTH_ll()) PZS(IIB:IIE,1)   = ZZSLS(IIB:IIE,1,1)
IF(HLBCY(2)/='CYCL' .AND. LNORTH_ll()) PZS(IIB:IIE,IJU) = ZZSLS(IIB:IIE,IJU,1)
!
NULLIFY(TZLSFIELD_ll)
CALL ADD2DFIELD_ll(TZLSFIELD_ll, PZS)
CALL UPDATE_HALO_ll(TZLSFIELD_ll,IINFO_ll)
CALL CLEANLIST_ll(TZLSFIELD_ll)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE ZS_BOUNDARY_n
