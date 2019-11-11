!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 newsrc 2006/12/12 15:06:50
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_MPDATA_SCALAR
!     #########################
INTERFACE
      SUBROUTINE MPDATA_SCALAR   ( KLITER, HLBCX, HLBCY, KSV,             &
                                   PTSTEP, PRHODJ, PSVM, PSVT,            &
                                   PRUCT, PRVCT, PRWCT, PRSVS             )                
!
INTEGER,                  INTENT(IN)    :: KLITER  ! Number of iterations MPDATA 
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVM
                                                ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT ! Contravariants
                                                ! components of the momentum
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT
                                                ! Variables at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS
                                                ! Sources terms
END SUBROUTINE MPDATA_SCALAR 
!
END INTERFACE
!
END MODULE MODI_MPDATA_SCALAR 
!
!
!     ########################################################################
      SUBROUTINE MPDATA_SCALAR   ( KLITER, HLBCX, HLBCY, KSV,             &
                                   PTSTEP, PRHODJ, PSVM, PSVT,            &
                                   PRUCT, PRVCT, PRWCT, PRSVS             )                
!     ########################################################################
!
!!****  *MPDATA_SCALAR* - routine to compute the advection tendancies of the scalar 
!!                 fields using an upstream scheme. The excesive numerical
!!                 correction of the scheme is corrected by means of an
!!                 antidiffusive velocity.
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the total advection 
!!    tendencies of all the scalar fields using the MPDATA scheme.
!!
!!**  METHOD
!!    ------
!!    MPDATA solves the advection of a quantity in the following way
!!              1.- 1st iteration. Upstream scheme.
!!                  The quantity is advected by the contravariant
!!                  velocities using an upstream scheme.
!!              2.- 2nd and next iterations. The excessive diffusion
!!                  of the upstream scheme is corrected by defining
!!                  the antidiffusive velocities (ANTI_DIFF routine)
!!                  and using for each iteration the upstream scheme.
!!    EXTERNAL
!!    --------
!!      ADD3DFIELD_ll  : add a field to 3D-list
!!      CLEANLIST_ll   : deallocate a list
!!      UPDATE_HALO_ll : update internal halos
!!      DXF,DYF,DZF    : Shuman operators
!!      FXM,FYM,FZM    : Flux operators 
!!      ANTI_DIFF      : antidiffusion
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1 of documentation ( MPDATA scheme )
!!
!!    AUTHOR
!!    ------
!!      J. Vila-Guerau      * Meteo France*
!!      J.-P. Lafore        * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/10/95 
!!      J.-P. Pinty & J.-M. Cohard *LA*   Add the budget calls
!!      J.    Stein                       include the cyclic case
!!      P. Jabouille                      parallelization
!!      V. Masson   06/11/02              updates the budget calls
!!      C.Lac                             Split meteorological scalar and tracer
!!                                        variables routines
!!      P.Tulet                           Upstream condition for aerosols
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
USE MODD_BUDGET
USE MODD_PARAMETERS
USE MODD_NSV, ONLY : NSV_DSTBEG, NSV_DSTEND, NSV_AERBEG, NSV_AEREND,&
                     NSV_SLTBEG, NSV_SLTEND
!
USE MODI_SHUMAN
USE MODI_FLUX
USE MODI_ANTI_DIFF
USE MODI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KLITER        ! Number of iterations MPDATA 
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KSV  ! Number of Scalar Variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVM
                                                  ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT
                                                  ! Contravariants components momentum
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS
                                                  ! Sources terms
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JLITER            ! Loop index for MPDATA iterations 
INTEGER :: JSV               ! Loop index for Scalar Variables
!
INTEGER:: IIB,IJB            ! Begining useful area  in x,y,z directions
INTEGER:: IIE,IJE            ! End useful area in x,y,z directions
INTEGER:: IKU
!  
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3))   :: ZGUESS  ! Guess 
                             ! variable (to be removed in the future !)
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZRAUCT
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZRAVCT
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZRAWCT
                             ! Antidiffusive contravariant component of the
                             ! momentum
!  
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZFADV   !   used
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZFADVU  !   for
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZFADVV  !  budget
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZFADVW  ! purpose
REAL, DIMENSION(SIZE(PRUCT,1),SIZE(PRUCT,2),SIZE(PRUCT,3)):: ZRVARS  !   only
!
INTEGER                :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
NULLIFY(TZFIELDS_ll)
!
!*       0. PROLOGUE
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
IKU=SIZE(PSVM,3)
!
!
!-------------------------------------------------------------------------------
!
!*        1.- case with KSV Scalar Variables
!             --------------------------------
  DO JSV=1,KSV
!
    CALL ADD3DFIELD_ll(TZFIELDS_ll, PRSVS(:,:,:,JSV))
    ZRVARS(:,:,:) = PRSVS(:,:,:,JSV)
    ZFADVU(:,:,:) = -DXF(FXM( PSVM(:,:,:,JSV),PRUCT(:,:,:) )  )
    ZFADVV(:,:,:) = -DYF(FYM( PSVM(:,:,:,JSV),PRVCT(:,:,:) )  )
    ZFADVW(:,:,:) = -DZF(1,IKU,1,FZM( PSVM(:,:,:,JSV),PRWCT(:,:,:) )  )
!
    PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) + ZFADVU(:,:,:) + ZFADVV(:,:,:) +  &
                                                          ZFADVW(:,:,:)
!
    ZRAUCT(:,:,:)=PRUCT(:,:,:)
    ZRAVCT(:,:,:)=PRVCT(:,:,:)
    ZRAWCT(:,:,:)=PRWCT(:,:,:) 
!                                                                

! ANTI_DIFF of MPDATA not suported by aerosols variables
  IF ((.NOT.((JSV .GE. NSV_AERBEG).AND.(JSV .LE. NSV_AEREND))).AND.&
      (.NOT.((JSV .GE. NSV_DSTBEG).AND.(JSV .LE. NSV_DSTEND))).AND.&
      (.NOT.((JSV .GE. NSV_SLTBEG).AND.(JSV .LE. NSV_SLTEND)))) THEN
  
    DO JLITER=2,KLITER
      CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
!
      CALL ANTI_DIFF(HLBCX,HLBCY,PTSTEP,PRHODJ,PRSVS(:,:,:,JSV),ZRAUCT,ZRAVCT,ZRAWCT)
!
      ZGUESS(:,:,:)=PTSTEP*PRSVS(:,:,:,JSV)/PRHODJ(:,:,:)
!
!
      ZFADV(:,:,:) = -DXF(FXM( ZGUESS(:,:,:),ZRAUCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVU(:,:,:) = ZFADVU(:,:,:) + ZFADV(:,:,:) 
      PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) + ZFADV(:,:,:)
!
      ZFADV(:,:,:) = -DYF(FYM( ZGUESS(:,:,:),ZRAVCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVV(:,:,:) = ZFADVV(:,:,:) + ZFADV(:,:,:) 
      PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) + ZFADV(:,:,:)
!
      ZFADV(:,:,:) = -DZF(1,IKU,1,FZM( ZGUESS(:,:,:),ZRAWCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVW(:,:,:) = ZFADVW(:,:,:) + ZFADV(:,:,:) 
      PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) + ZFADV(:,:,:)
    END DO
  END IF
!
    CALL CLEANLIST_ll(TZFIELDS_ll)
!
    IF (LBUDGET_SV)  THEN
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVU(:,:,:)
      CALL BUDGET (ZRVARS,JSV+12,'ADVX_BU_RSV')
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVV(:,:,:)
      CALL BUDGET (ZRVARS,JSV+12,'ADVY_BU_RSV')
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVW(:,:,:)
      CALL BUDGET (ZRVARS,JSV+12,'ADVZ_BU_RSV')
    END IF
!
  END DO

!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MPDATA_SCALAR  
