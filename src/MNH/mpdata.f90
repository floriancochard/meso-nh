!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/12/12 15:06:37
!-----------------------------------------------------------------
!     ##################
      MODULE MODI_MPDATA
!     ##################
INTERFACE
      SUBROUTINE MPDATA   (KLITER, HLBCX, HLBCY, KRR,                   &
                           PTSTEP, PRHODJ, PTHM, PRM, PTKEM,            &
                           PTHT, PRT, PTKET,                            &
                           PRUCT, PRVCT, PRWCT,                         &
                           PRTHS, PRRS, PRTKES                          )
!
INTEGER,                  INTENT(IN)    :: KLITER  ! Number of iterations MPDATA 
!
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
!
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM 
                                                ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT ! Contravariants
                                                ! components of the momentum
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET, PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
                                                ! Sources terms
END SUBROUTINE MPDATA 
!
END INTERFACE
!
END MODULE MODI_MPDATA 
!
!
!     ########################################################################
      SUBROUTINE MPDATA   (KLITER, HLBCX, HLBCY, KRR,                   &
                           PTSTEP, PRHODJ, PTHM, PRM, PTKEM,            &
                           PTHT, PRT, PTKET,                            &
                           PRUCT, PRVCT, PRWCT,                         &
                           PRTHS, PRRS, PRTKES                          )
!     ########################################################################
!
!!****  *MPDATA* - routine to compute the advection tendancies of the scalar 
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
!!      05/2006                           Remove EPS
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
INTEGER,                  INTENT(IN)    :: KRR  ! Number of moist variables
!
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHM, PTKEM
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRM 
                                                  ! Variables at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRUCT,PRVCT,PRWCT
                                                  ! Contravariants components momentum
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PTKET, PRHODJ
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT 
                                                  ! Variables at t
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS, PRTKES
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS 
                                                  ! Sources terms
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JLITER            ! Loop index for MPDATA iterations 
INTEGER :: JRR               ! Loop index for moist variables
!
INTEGER:: IIB,IJB            ! Begining useful area  in x,y,z directions
INTEGER:: IIE,IJE            ! End useful area in x,y,z directions
INTEGER:: IKU
!  
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3))   :: ZGUESS  ! Guess 
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
CHARACTER (LEN=3) , DIMENSION(7)                          :: YRX
CHARACTER (LEN=20)                                        :: YBURX
LOGICAL           , DIMENSION(7)                          :: LBUDGET_R
!
INTEGER                :: IINFO_ll      ! return code of parallel routine
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!JUAN : init of TZFIELDS_ll
NULLIFY(TZFIELDS_ll)
!
!*       0.3 PROLOGUE
!
CALL GET_PHYSICAL_ll(IIB,IJB,IIE,IJE)
IKU=SIZE(PTHM,3)
!
YRX(1) = 'RRV'
YRX(2) = 'RRC'
YRX(3) = 'RRR'
YRX(4) = 'RRI'
YRX(5) = 'RRS'
YRX(6) = 'RRG'
YRX(7) = 'RRH'
!
LBUDGET_R(1) = LBUDGET_RV
LBUDGET_R(2) = LBUDGET_RC
LBUDGET_R(3) = LBUDGET_RR
LBUDGET_R(4) = LBUDGET_RI
LBUDGET_R(5) = LBUDGET_RS
LBUDGET_R(6) = LBUDGET_RG
LBUDGET_R(7) = LBUDGET_RH
!
!
!-------------------------------------------------------------------------------
!
!
!
!*       1.  Thermodynamical variable
!             -----------------------
!
  CALL ADD3DFIELD_ll(TZFIELDS_ll, PRTHS)
!
!* 1st iteration (upstream scheme)
!
  ZRVARS(:,:,:) = PRTHS(:,:,:)
  ZFADVU(:,:,:) = -DXF(FXM( PTHM(:,:,:),PRUCT(:,:,:) )  )
  ZFADVV(:,:,:) = -DYF(FYM( PTHM(:,:,:),PRVCT(:,:,:) )  )
  ZFADVW(:,:,:) = -DZF(1,IKU,1,FZM( PTHM(:,:,:),PRWCT(:,:,:) )  )
!
  PRTHS(:,:,:) = PRTHS(:,:,:) + ZFADVU(:,:,:) + ZFADVV(:,:,:) + ZFADVW(:,:,:)
!
!
!* Iterations greater than 1 
!
  ZRAUCT(:,:,:)=PRUCT(:,:,:)
  ZRAVCT(:,:,:)=PRVCT(:,:,:)
  ZRAWCT(:,:,:)=PRWCT(:,:,:) 
!                                                                
  DO JLITER=2,KLITER
!
!   update halo (and possibly periodize) the guess of the future time
!
    CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
!
    CALL ANTI_DIFF(HLBCX,HLBCY,PTSTEP,PRHODJ,PRTHS,ZRAUCT,ZRAVCT,ZRAWCT)          
!
    ZGUESS(:,:,:)=PTSTEP*PRTHS/PRHODJ(:,:,:)
!
    ZFADV(:,:,:) = -DXF(FXM( ZGUESS(:,:,:),ZRAUCT(:,:,:) )  )
    IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
    IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
    IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
    IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
    ZFADVU(:,:,:) = ZFADVU(:,:,:) + ZFADV(:,:,:) 
    PRTHS(:,:,:) = PRTHS(:,:,:) + ZFADV(:,:,:)
!
    ZFADV(:,:,:) = -DYF(FYM( ZGUESS(:,:,:),ZRAVCT(:,:,:) )  )
    IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
    IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
    IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
    IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
    ZFADVV(:,:,:) = ZFADVV(:,:,:) + ZFADV(:,:,:) 
    PRTHS(:,:,:) = PRTHS(:,:,:) + ZFADV(:,:,:)
!
    ZFADV(:,:,:) = -DZF(1,IKU,1,FZM( ZGUESS(:,:,:),ZRAWCT(:,:,:) )  )
    IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
    IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
    IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
    IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
    ZFADVW(:,:,:) = ZFADVW(:,:,:) + ZFADV(:,:,:) 
    PRTHS(:,:,:) = PRTHS(:,:,:) + ZFADV(:,:,:)
! 
  END DO
!
  CALL CLEANLIST_ll(TZFIELDS_ll)
!
  IF (LBUDGET_TH)  THEN
    ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVU(:,:,:)
    CALL BUDGET (ZRVARS,4,'ADVX_BU_RTH')
    ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVV(:,:,:)
    CALL BUDGET (ZRVARS,4,'ADVY_BU_RTH')
    ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVW(:,:,:)
    CALL BUDGET (ZRVARS,4,'ADVZ_BU_RTH')
  END IF
!
!-------------------------------------------------------------------------------
!
!*        2.  Case with KRR water variables 
!             -----------------------------
!
  DO JRR=1,KRR
    CALL ADD3DFIELD_ll(TZFIELDS_ll, PRRS(:,:,:,JRR))
    ZRVARS(:,:,:) = PRRS(:,:,:,JRR)
    ZFADVU(:,:,:) = -DXF(FXM( PRM(:,:,:,JRR),PRUCT(:,:,:) )  )
    ZFADVV(:,:,:) = -DYF(FYM( PRM(:,:,:,JRR),PRVCT(:,:,:) )  )
    ZFADVW(:,:,:) = -DZF(1,IKU,1,FZM( PRM(:,:,:,JRR),PRWCT(:,:,:) )  )
!
    PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR) + ZFADVU(:,:,:) + ZFADVV(:,:,:) +  &
                                     ZFADVW(:,:,:)
!
    ZRAUCT(:,:,:)=PRUCT(:,:,:)
    ZRAVCT(:,:,:)=PRVCT(:,:,:)
    ZRAWCT(:,:,:)=PRWCT(:,:,:) 
!                                                                
    DO JLITER=2,KLITER
      CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
!
      CALL ANTI_DIFF(HLBCX,HLBCY,PTSTEP,PRHODJ,PRRS(:,:,:,JRR),ZRAUCT,ZRAVCT,ZRAWCT)
!
      ZGUESS(:,:,:)=PTSTEP*PRRS(:,:,:,JRR)/PRHODJ(:,:,:)
!
!
      ZFADV(:,:,:) = -DXF(FXM( ZGUESS(:,:,:),ZRAUCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVU(:,:,:) = ZFADVU(:,:,:) + ZFADV(:,:,:) 
      PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR) + ZFADV(:,:,:)
!
      ZFADV(:,:,:) = -DYF(FYM( ZGUESS(:,:,:),ZRAVCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVV(:,:,:) = ZFADVV(:,:,:) + ZFADV(:,:,:) 
      PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR) + ZFADV(:,:,:)
!
      ZFADV(:,:,:) = -DZF(1,IKU,1,FZM( ZGUESS(:,:,:),ZRAWCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVW(:,:,:) = ZFADVW(:,:,:) + ZFADV(:,:,:) 
      PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR) + ZFADV(:,:,:)
    END DO
!
    CALL CLEANLIST_ll(TZFIELDS_ll)
!
    IF (LBUDGET_R(JRR)) THEN
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVU(:,:,:)
      YBURX = 'ADVX_BU_'//YRX(JRR)
      CALL BUDGET (ZRVARS(:,:,:),JRR+5 ,YBURX)
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVV(:,:,:)
      YBURX = 'ADVY_BU_'//YRX(JRR)
      CALL BUDGET (ZRVARS(:,:,:),JRR+5 ,YBURX)
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVW(:,:,:)
      YBURX = 'ADVZ_BU_'//YRX(JRR)
      CALL BUDGET (ZRVARS(:,:,:),JRR+5 ,YBURX)
    END IF
  END DO
!
!
!-------------------------------------------------------------------------------
!
!*        3.  TKE variable
!             -------------
  IF (SIZE(PTKET,1) /= 0) THEN
!
    CALL ADD3DFIELD_ll(TZFIELDS_ll, PRTKES)
    ZRVARS(:,:,:) = PRTKES(:,:,:)
    ZFADVU(:,:,:) = -DXF(FXM( PTKEM(:,:,:),PRUCT(:,:,:) )  )
    ZFADVV(:,:,:) = -DYF(FYM( PTKEM(:,:,:),PRVCT(:,:,:) )  )
    ZFADVW(:,:,:) = -DZF(1,IKU,1,FZM( PTKEM(:,:,:),PRWCT(:,:,:) )  )
!
    PRTKES(:,:,:) = PRTKES(:,:,:) + ZFADVU(:,:,:) + ZFADVV(:,:,:) + ZFADVW(:,:,:)
!
    ZRAUCT(:,:,:)=PRUCT(:,:,:)
    ZRAVCT(:,:,:)=PRVCT(:,:,:)
    ZRAWCT(:,:,:)=PRWCT(:,:,:) 
!                                                                
    DO JLITER=2,KLITER
      CALL UPDATE_HALO_ll(TZFIELDS_ll, IINFO_ll)
!
      CALL ANTI_DIFF(HLBCX,HLBCY,PTSTEP,PRHODJ,PRTKES,ZRAUCT,ZRAVCT,ZRAWCT)
!
      ZGUESS(:,:,:)=PTSTEP*PRTKES/PRHODJ(:,:,:)
!
!
      ZFADV(:,:,:) = -DXF(FXM( ZGUESS(:,:,:),ZRAUCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVU(:,:,:) = ZFADVU(:,:,:) + ZFADV(:,:,:) 
      PRTKES(:,:,:) = PRTKES(:,:,:) + ZFADV(:,:,:)
!
      ZFADV(:,:,:) = -DYF(FYM( ZGUESS(:,:,:),ZRAVCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVV(:,:,:) = ZFADVV(:,:,:) + ZFADV(:,:,:) 
      PRTKES(:,:,:) = PRTKES(:,:,:) + ZFADV(:,:,:)
!
      ZFADV(:,:,:) = -DZF(1,IKU,1,FZM( ZGUESS(:,:,:),ZRAWCT(:,:,:) )  )
      IF(LWEST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIB,:,:)=0.
      IF(LEAST_ll()  .AND. HLBCX(1) /= 'CYCL') ZFADV(IIE,:,:)=0.
      IF(LSOUTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJB,:)=0.
      IF(LNORTH_ll() .AND. HLBCY(1) /= 'CYCL') ZFADV(:,IJE,:)=0.
      ZFADVW(:,:,:) = ZFADVW(:,:,:) + ZFADV(:,:,:) 
      PRTKES(:,:,:) = PRTKES(:,:,:) + ZFADV(:,:,:)
    END DO
!
    CALL CLEANLIST_ll(TZFIELDS_ll)
!  
    IF (LBUDGET_TKE)  THEN
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVU(:,:,:)
      CALL BUDGET (ZRVARS,5,'ADVX_BU_RTKE')
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVV(:,:,:)
      CALL BUDGET (ZRVARS,5,'ADVY_BU_RTKE')
      ZRVARS(:,:,:) = ZRVARS(:,:,:) + ZFADVW(:,:,:)
      CALL BUDGET (ZRVARS,5,'ADVZ_BU_RTKE')
    END IF
  END IF
!
!-------------------------------------------------------------------------------
!  
!
END SUBROUTINE MPDATA  
