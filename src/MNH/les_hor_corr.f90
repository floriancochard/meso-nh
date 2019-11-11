!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 les 2006/05/18 13:07:25
!-----------------------------------------------------------------
!#######################
MODULE MODI_LES_HOR_CORR
!#######################
!
INTERFACE
!
      SUBROUTINE  LES_HOR_CORR(PA,PB,HLBCX,HLBCY,PRX,PRY)

REAL,             DIMENSION(:,:), INTENT(IN)   :: PA       ! first variable
REAL,             DIMENSION(:,:), INTENT(IN)   :: PB       ! second variable
CHARACTER(LEN=4), DIMENSION(2),   INTENT(IN)   :: HLBCX    ! boundary condition
CHARACTER(LEN=4), DIMENSION(2),   INTENT(IN)   :: HLBCY    ! in x and y
REAL,             DIMENSION(:),   INTENT(OUT)  :: PRX      ! 2 points corr.
REAL,             DIMENSION(:),   INTENT(OUT)  :: PRY      ! in x and y dir.

END SUBROUTINE LES_HOR_CORR
!
END INTERFACE
!
END MODULE MODI_LES_HOR_CORR
!
!
!     ###################################################
      SUBROUTINE  LES_HOR_CORR(PA,PB,HLBCX,HLBCY,PRX,PRY)
!     ###################################################
!
!
!!****  *LES_HOR_CORR* computes the current time-step 2 points HORizontal
!!                     CORRelations diagnostics, for a two horizontal fields
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!    The fields are already gathered between the processors: fields
!!    PA and PB are on the entire physical simulation domain (or cartesian
!!    subdomain if any).
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!                       01/04/03 (V. Masson) bug in non-cyclic case
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_PARAMETERS
USE MODD_CONF

!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!

REAL,             DIMENSION(:,:), INTENT(IN)   :: PA       ! first variable
REAL,             DIMENSION(:,:), INTENT(IN)   :: PB       ! second variable
CHARACTER(LEN=4), DIMENSION(2),   INTENT(IN)   :: HLBCX    ! boundary condition
CHARACTER(LEN=4), DIMENSION(2),   INTENT(IN)   :: HLBCY    ! in x and y
REAL,             DIMENSION(:),   INTENT(OUT)  :: PRX      ! 2 points corr.
REAL,             DIMENSION(:),   INTENT(OUT)  :: PRY      ! in x and y dir.
!
!       0.2  declaration of local variables
!
INTEGER :: IIMAX_ll ! fields dimensions
INTEGER :: IJMAX_ll ! in x and y dir.
!
INTEGER :: JI       ! loop counter along x direction
INTEGER :: JJ       ! loop counter along y direction
INTEGER :: JL       ! loop counter for wavenumbers increments
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZA  ! extended A field
REAL, DIMENSION(:,:), ALLOCATABLE :: ZB  ! extended B field
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZW1 ! working array
REAL, DIMENSION(:),   ALLOCATABLE :: ZW2 ! working array
!
!-------------------------------------------------------------------------------
!
!* intersection of the level with orography (when used of diagnostics
!  in real cases)
!
IF (ANY(PA(:,:)==XUNDEF .OR. PB(:,:)==XUNDEF)) THEN
  PRX(:)=XUNDEF
  PRY(:)=XUNDEF
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
IIMAX_ll = SIZE(PA,1)
IJMAX_ll = SIZE(PA,2)
!
!-------------------------------------------------------------------------------
!
!       1.   computation of correlations (x direction)
!            -----------------------------------------
!
IF (HLBCX(1)=='CYCL') THEN
!
!* cyclic case
!  -----------
!
  !* extension of the data arrays
  !
  ALLOCATE(ZA(IIMAX_ll+(IIMAX_ll+1)/2,IJMAX_ll))
  ALLOCATE(ZB(IIMAX_ll+(IIMAX_ll+1)/2,IJMAX_ll))
  !
  ZA(   1      :  IIMAX_ll      ,:) = PA(:              ,:)
  ZB(   1      :  IIMAX_ll      ,:) = PB(:              ,:)
  ZA(IIMAX_ll+1:IIMAX_ll+(IIMAX_ll+1)/2,:) = PA(:(IIMAX_ll+1)/2,:)
  ZB(IIMAX_ll+1:IIMAX_ll+(IIMAX_ll+1)/2,:) = PB(:(IIMAX_ll+1)/2,:)
  !
  !* computation of correlations in both directions along axis i
  !
  ALLOCATE(ZW1(IIMAX_ll/2+1))
  ALLOCATE(ZW2(IIMAX_ll/2+1))
  !
  ZW1(:) = 0.
  ZW2(:) = 0.
  !
  DO JL=1,IIMAX_ll/2+1
    !
    !* computation for wave length JL
    !
    DO JI=1,IIMAX_ll
      DO JJ=1,IJMAX_ll
        ZW1(JL) = ZW1(JL) + ZA(JI     ,JJ)*ZB(JI+JL-1,JJ)
        ZW2(JL) = ZW2(JL) + ZA(JI+JL-1,JJ)*ZB(JI     ,JJ)
      END DO
    END DO
    !
    !* merging of correlations in both directions along axis i
    !
    PRX(JL) = ( ZW1(JL) + ZW2(JL) ) / (2.*IIMAX_ll) / IJMAX_ll
    !
    !* complete the large wave length by symetry
    !
    IF (JL>1) PRX(IIMAX_ll+2-JL) = PRX(JL)
    !
  END DO
  !
  DEALLOCATE(ZA)
  DEALLOCATE(ZB)
  DEALLOCATE(ZW1)
  DEALLOCATE(ZW2)
  !
ELSE
!
!* no cyclic case
!  --------------
!
  !* computation of correlations in both directions along axis i
  !
  ALLOCATE(ZW1(IIMAX_ll))
  ALLOCATE(ZW2(IIMAX_ll))
  !
  ZW1(:) = 0.
  ZW2(:) = 0.
  !
  DO JL=1,IIMAX_ll
    !
    !* computation for wave length JL
    !
    DO JI=1,IIMAX_ll-(JL-1)
      DO JJ=1,IJMAX_ll
        ZW1(JL) = ZW1(JL) + PA(JI     ,JJ)*PB(JI+JL-1,JJ)
        ZW2(JL) = ZW2(JL) + PA(JI+JL-1,JJ)*PB(JI     ,JJ)
      END DO
    END DO
    !
    !* merging of correlations in both directions along axis i
    !
    PRX(JL) = ( ZW1(JL) + ZW2(JL) ) / (2.*(IIMAX_ll-(JL-1))) / IJMAX_ll
    !
  END DO
  !
  DEALLOCATE(ZW1)
  DEALLOCATE(ZW2)
  !
END IF
!
!-------------------------------------------------------------------------------
!
!       2.   computation of correlations (y direction)
!            -----------------------------------------
!
IF (L2D) THEN
  PRY(:) = XUNDEF
  RETURN
END IF
!
IF (HLBCY(1)=='CYCL') THEN
!
!* cyclic case
!  -----------
!
  !* extension of the data arrays
  !
  ALLOCATE(ZA(IIMAX_ll,IJMAX_ll+(IJMAX_ll+1)/2))
  ALLOCATE(ZB(IIMAX_ll,IJMAX_ll+(IJMAX_ll+1)/2))
  !
  ZA(:,   1      :  IJMAX_ll      ) = PA(:,:              )
  ZB(:,   1      :  IJMAX_ll      ) = PB(:,:              )
  ZA(:,IJMAX_ll+1:IJMAX_ll+(IJMAX_ll+1)/2) = PA(:,:(IJMAX_ll+1)/2)
  ZB(:,IJMAX_ll+1:IJMAX_ll+(IJMAX_ll+1)/2) = PB(:,:(IJMAX_ll+1)/2)
  !
  !* computation of correlations in both directions along axis i
  !
  ALLOCATE(ZW1(IJMAX_ll/2+1))
  ALLOCATE(ZW2(IJMAX_ll/2+1))
  !
  ZW1(:) = 0.
  ZW2(:) = 0.
  !
  DO JL=1,IJMAX_ll/2+1
    !
    !* computation for wave length JL
    !
    DO JJ=1,IJMAX_ll
      DO JI=1,IIMAX_ll
        ZW1(JL) = ZW1(JL) + ZA(JI,JJ     )*ZB(JI,JJ+JL-1)
        ZW2(JL) = ZW2(JL) + ZA(JI,JJ+JL-1)*ZB(JI,JJ     )
      END DO
    END DO
    !
    !* merging of correlations in both directions along axis i
    !
    PRY(JL) = ( ZW1(JL) + ZW2(JL) ) / (2.*IJMAX_ll) / IIMAX_ll
    !
    !* complete the large wave length by symetry
    !
    IF (JL>1) PRY(IJMAX_ll+2-JL) = PRY(JL)
    !
  END DO
  !
  DEALLOCATE(ZA)
  DEALLOCATE(ZB)
  DEALLOCATE(ZW1)
  DEALLOCATE(ZW2)
  !
ELSE
!
!* no cyclic case
!  --------------
!
  !* computation of correlations in both directions along axis i
  !
  ALLOCATE(ZW1(IJMAX_ll))
  ALLOCATE(ZW2(IJMAX_ll))
  !
  ZW1(:) = 0.
  ZW2(:) = 0.
  !
  DO JL=1,IJMAX_ll
    !
    !* computation for wave length JL
    !
    DO JJ=1,IJMAX_ll-(JL-1)
      DO JI=1,IIMAX_ll
        ZW1(JL) = ZW1(JL) + PA(JI,JJ     )*PB(JI,JJ+JL-1)
        ZW2(JL) = ZW2(JL) + PA(JI,JJ+JL-1)*PB(JI,JJ     )
      END DO
    END DO
    !
    !* merging of correlations in both directions along axis i
    !
    PRY(JL) = ( ZW1(JL) + ZW2(JL) ) / (2.*(IJMAX_ll-(JL-1))) / IIMAX_ll
    !
  END DO
  !
  DEALLOCATE(ZW1)
  DEALLOCATE(ZW2)
  !
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_HOR_CORR
