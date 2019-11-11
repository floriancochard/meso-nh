!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 solver 2006/05/18 14:41:09
!-----------------------------------------------------------------
!####################
MODULE MODI_MASS_LEAK
!####################
!
INTERFACE
!
          SUBROUTINE MASS_LEAK (PDXX,PDYY,HLBCX,HLBCY,PLINMASS,PRHODJ,PRUS,PRVS)
!
REAL, DIMENSION(:,:,:),         INTENT(IN) :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),         INTENT(IN) :: PDYY    ! metric coefficient dyy
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX   ! type of lateral boundary
!                                                     ! condition (i=IB, i=IE+1)
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY   ! type of lateral boundary
!                                                     ! condition (j=JB, j=JE+1)
REAL,                          INTENT(IN) :: PLINMASS ! lineic mass through open
!                                                     ! boundaries
REAL, DIMENSION(:,:,:),        INTENT(IN)  :: PRHODJ  ! rhodref*J
REAL, DIMENSION(:,:,:),        INTENT(INOUT) :: PRUS  ! Horizontal
REAL, DIMENSION(:,:,:),        INTENT(INOUT) :: PRVS  ! momentum  tendencies
!
END SUBROUTINE MASS_LEAK
!
END INTERFACE
!
END MODULE MODI_MASS_LEAK
!
!     #####################################################################
      SUBROUTINE MASS_LEAK(PDXX,PDYY,HLBCX,HLBCY,PLINMASS,PRHODJ,PRUS,PRVS)
!     #####################################################################
!
!!***   *MASS_LEAK* - assures global non-divergence condition
!!
!!    PURPOSE
!!    -------
!!
!!      This routine changes the horizontal reference dry mass fluxes through
!!    the open boundary condition faces to set the global divergence in the
!!    model domain to zero.
!!
!!**  METHOD
!!    ------
!!
!!  1) The leak term is computed as:
!!
!!          --                     -- IE+1   --                     -- JE+1
!!          | JE KE                 |        | IE KE                 |
!!          | _  _                  |        | _  _                  |
!!          | \  \    1        _ _  |        | \  \    1        _ _  |
!!   ZLEAK= | /  /   --- PRUS dydz  |   +    | /  /   --- PRVS dxdz  |
!!          | -  -   dxx            |        | -  -   dyy            |
!!          | JB KB                 |        | IB KB                 |
!!          --                     -- i=IB   --                     -- j=JB
!!
!!  2) Then the correction wind value ZRUSTOP is found as
!!               ZLEAK
!!    ZRUSTOP= ----------
!!              PLINMASS
!!
!!     where PLINMASS is the lineic mass through the open lateral boundaries.
!!
!!  3) This horizontal wind correction is applied on the normal wind of
!!     open lateral boundaries only.
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS : contains declaration of parameter variables
!!
!!         JPHEXT   : Horizontal external points number
!!         JPVEXT   : Vertical external points number
!!
!!    REFERENCE
!!    ---------
!!
!!      book2
!!
!!    AUTHOR
!!    ------
!!
!!      V. Masson   Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     09/02/95
!!      Modification 20/10/97  (J.P.Lafore) introduction of 'DAVI' type of lbc
!!                  15/06/98  (D.Lugato, R.Guivarch) Parallelisation
!!                   05/06    Suppression of Davies type of lbc
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
USE MODE_ll
!JUAN
USE MODE_REPRO_SUM
!JUAN
!
IMPLICIT NONE
!
!*       0.1   declarations of dummy arguments
!
REAL, DIMENSION(:,:,:),        INTENT(IN) :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),        INTENT(IN) :: PDYY    ! metric coefficient dyy
CHARACTER(LEN=4), DIMENSION(2),INTENT(IN) :: HLBCX   ! type of lateral boundary
!                                                    ! condition (i=IB, i=IE+1)
CHARACTER(LEN=4), DIMENSION(2),INTENT(IN) :: HLBCY   ! type of lateral boundary
!                                                    ! condition (j=JB, j=JE+1)
REAL,                          INTENT(IN) :: PLINMASS! lineic mass through open
!                                                    ! boundaries
REAL, DIMENSION(:,:,:),        INTENT(IN) :: PRHODJ  ! rhodref*J
REAL, DIMENSION(:,:,:),        INTENT(INOUT) :: PRUS ! Horizontal
REAL, DIMENSION(:,:,:),        INTENT(INOUT) :: PRVS ! momentum  tendencies
!
!*       0.2   declarations of local variables
!
!JUAN16
REAL                               :: ZLEAK     ! total leak of mass
REAL, ALLOCATABLE, DIMENSION (:,:) :: ZLEAK_W_2D , ZLEAK_E_2D , ZLEAK_S_2D , ZLEAK_N_2D
!JUAN16

REAL                :: ZUSTOP     ! wind correction!
INTEGER             :: IIB        ! indice I Beginning in x direction
INTEGER             :: IJB        ! indice J Beginning in y direction
INTEGER             :: IKB        ! indice K Beginning in z direction
INTEGER             :: IIE        ! indice I End       in x direction
INTEGER             :: IJE        ! indice J End       in y direction
INTEGER             :: IKE        ! indice K End       in z direction
INTEGER             :: JI         ! Loop index in x direction
INTEGER             :: JJ         ! Loop index in y direction
INTEGER             :: JK         ! Loop index in z direction
!
INTEGER             :: IINFO_ll   ! return code of parallel routine

!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES:
!              ----------------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
IKB = 1 + JPVEXT
IKE = SIZE(PRUS,3) - JPVEXT
!
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTATION OF LEAK
!              -------------------
!
ZLEAK=0.
!
IF( HLBCX(1) /= 'CYCL' ) THEN
   ALLOCATE( ZLEAK_W_2D(IIB:IIB,IJB:IJE))
   ZLEAK_W_2D = 0.0
   IF (LWEST_ll()) THEN
      DO JK=IKB,IKE
         DO JJ=IJB,IJE
            ZLEAK_W_2D(IIB,JJ) = ZLEAK_W_2D(IIB,JJ) - 1./PDXX(IIB,JJ,JK) *PRUS(IIB,JJ,JK)
         END DO
      END DO      
   END IF
   ZLEAK =         SUM_DD_R2_ll(ZLEAK_W_2D)
!
  ALLOCATE( ZLEAK_E_2D(IIE+1:IIE+1,IJB:IJE))
  ZLEAK_E_2D = 0.0
  IF (LEAST_ll()) THEN
     DO JK=IKB,IKE
        DO JJ=IJB,IJE
           ZLEAK_E_2D(IIE+1,JJ) = ZLEAK_E_2D(IIE+1,JJ) + 1./PDXX(IIE+1,JJ,JK)*PRUS(IIE+1,JJ,JK)
        END DO
     END DO
  END IF
  ZLEAK = ZLEAK + SUM_DD_R2_ll(ZLEAK_E_2D)  
!
END IF
!
IF( HLBCY(1) /= 'CYCL' ) THEN
   ALLOCATE( ZLEAK_S_2D(IIB:IIE,IJB:IJB))
   ZLEAK_S_2D = 0.0
   IF (LSOUTH_ll()) THEN
      DO JI=IIB,IIE
         DO JK=IKB,IKE
            ZLEAK_S_2D(JI,IJB) = ZLEAK_S_2D(JI,IJB) - 1./PDYY(JI,IJB,JK) *PRVS(JI,IJB,JK)
         END DO
      END DO
   END IF
   ZLEAK = ZLEAK + SUM_DD_R2_ll(ZLEAK_S_2D)  
   !
   ALLOCATE( ZLEAK_N_2D(IIB:IIE,IJE+1:IJE+1))  
   ZLEAK_N_2D = 0.0
   IF (LNORTH_ll()) THEN
      DO JI=IIB,IIE
         DO JK=IKB,IKE
            ZLEAK_N_2D(JI,IJE+1) = ZLEAK_N_2D(JI,IJE+1)  + 1./PDYY(JI,IJE+1,JK)*PRVS(JI,IJE+1,JK)
         END DO
      END DO
   END IF
   ZLEAK = ZLEAK + SUM_DD_R2_ll(ZLEAK_N_2D)  
   !
END IF
!
!CALL REDUCESUM_ll(ZLEAK,IINFO_ll)	! we do the reducesum_ll in SUM_DD_R2_ll so we do not do it here
!
!-------------------------------------------------------------------------------
!
!*       3.    CORRECTION OF WIND ON OPEN BOUNDARIES
!              -------------------------------------
!
ZUSTOP=ZLEAK
ZUSTOP=ZUSTOP/PLINMASS
!
IF (HLBCX(1)=='OPEN' .AND. LWEST_ll() ) &
 PRUS(IIB,:,:)=PRUS(IIB,:,:)+ZUSTOP*0.5*(PRHODJ(IIB,:,:)+PRHODJ(IIB-1,:,:))
!
IF (HLBCX(2)=='OPEN' .AND. LEAST_ll() ) &
 PRUS(IIE+1,:,:)=PRUS(IIE+1,:,:)-ZUSTOP*0.5*(PRHODJ(IIE+1,:,:)+PRHODJ(IIE,:,:))
!
IF (HLBCY(1)=='OPEN' .AND. LSOUTH_ll() ) &
 PRVS(:,IJB,:)=PRVS(:,IJB,:)+ZUSTOP*0.5*(PRHODJ(:,IJB,:)+PRHODJ(:,IJB-1,:))
!
IF (HLBCY(2)=='OPEN' .AND. LNORTH_ll() ) &
 PRVS(:,IJE+1,:)=PRVS(:,IJE+1,:)-ZUSTOP*0.5*(PRHODJ(:,IJE+1,:)+PRHODJ(:,IJE,:))
!
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MASS_LEAK
