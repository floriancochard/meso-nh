
!###########################################################
     SUBROUTINE BLOWSNW_DEP (K2D_SNWBEG, K2D_SNWEND,PTA,PPA,              &
             PRHODREF,PSVT,PSFSNW)
!###########################################################
  !-------------------------------------------------------------------------------
  !
  !*       0.    DECLARATIONS
  !              ------------
  !
  USE MODE_BLOWSNW_SURF
  USE MODI_BLOWSNW_VELGRAV1D
  USE MODD_CSTS, ONLY : XG,XRHOLI, XPI 
  USE MODD_BLOWSNW_SURF

  !
  IMPLICIT NONE
  !
  !*       0.1   Declarations of dummy arguments :
  !
   INTEGER,      INTENT(IN)              :: K2D_SNWBEG, K2D_SNWEND ! index of first and last blowing snow 2D variable sent to MNH
   REAL, DIMENSION(:),     INTENT(IN)    :: PTA        ! air temperature
   REAL, DIMENSION(:),     INTENT(IN)    :: PPA          ! air pressure
   REAL, DIMENSION(:),     INTENT(IN)    :: PRHODREF   ! air density
   REAL, DIMENSION(:,:),   INTENT(INOUT)    :: PSVT       ! blowing snow concentration
                                                              ! in surface units (__/m3)
  REAL, DIMENSION(:,:),   INTENT(INOUT)    :: PSFSNW      ! blowing snow concentration fluxes
                                                             



 !
  !
  !*       0.2   Declarations of local variables :
  !
      REAL, DIMENSION(SIZE(PSVT,1),1, SIZE(PSVT,2)) :: ZSVT
      REAL, DIMENSION(SIZE(PSVT,1),1) :: ZBET, ZRG,ZTA,ZPA
      REAL, DIMENSION(SIZE(PSVT,1),1,SIZE(PSVT,2)) :: ZVGK   ! Terminal fallspeed (m/s)

      INTEGER                 :: JN


! Save scalars in local array
DO JN=1,SIZE(PSVT,2)
    ZSVT(:,1,JN)=MAX(PSVT(:,JN),1E-60)
END DO
ZTA(:,1)=PTA(:)
ZPA(:,1)=PPA(:)

! Get gamma parameter distribution : scale factor and mean radius at the first
!                                    level in the atmosphere
CALL SNOWMOMENT2SIZE(ZSVT, PBETA1D=ZBET, PRG1D=ZRG ) 
! Compute mass-weighted terminal fall speed based on particles distribution and
! atmospheric conditions.
CALL BLOWSNW_VELGRAV1D(ZBET, ZRG, ZTA, PRHODREF,ZPA, ZVGK)

! Compute sedimentation fluxes of blowing snow variables

PSVT(:,K2D_SNWBEG) = PSVT(:,1)  * ZVGK(:,1,1)
PSVT(:,K2D_SNWBEG+1) = PSVT(:,2)  * ZVGK(:,1,2)


PSFSNW(:,1:2) = PSVT(:,K2D_SNWBEG:(K2D_SNWBEG+1))


END SUBROUTINE BLOWSNW_DEP
