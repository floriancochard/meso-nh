!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########################################
      SUBROUTINE CANOPY_GRID_UPDATE(KI,PH,PZFORC,SB,OLOG_GRID)
!     #########################################
!
!!****  *CANOPY_GRID_UPDATE* - set the upper levels at and just below forcing level
!!                        
!!
!!    PURPOSE
!!    -------
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2006 
!!      S. Riette   Oct 2010 Vectorisation
!-------------------------------------------------------------------------------
!
USE MODD_CANOPY_n, ONLY : CANOPY_t
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CANOPY_GRID
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                  INTENT(IN)    :: KI        ! number of horizontal points
REAL, DIMENSION(KI),      INTENT(IN)    :: PH        ! maximum canopy height                 (m)
REAL, DIMENSION(KI),      INTENT(IN)    :: PZFORC    ! height of wind forcing                (m)
!
TYPE(CANOPY_t), INTENT(INOUT) :: SB
!
LOGICAL, OPTIONAL,        INTENT(IN)    :: OLOG_GRID ! true if logarithmic grid is used (blowing snow scheme)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER, DIMENSION(KI)      :: IL     ! latest level below forcing height
INTEGER, DIMENSION(KI,SB%NLVL) :: ILEVEL ! to test if level is high enough
!
INTEGER :: ICOUNT                 ! number of layers above forcing height, these must be changed
INTEGER :: JLAYER                 ! loop counter on layers
INTEGER :: JI                     ! loop counter on points
REAL    :: ZZTOP                  ! altitude of top of the grid of the initial level
!                                 ! just below forcing height
REAL    :: ZDZ                    ! difference of height between new levels
LOGICAL :: GLOG_GRID
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CANOPY_GRID_UPDATE',0,ZHOOK_HANDLE)
!
IF(ALL(SB%XZ(:,SB%NLVL)==PZFORC(:)) .AND. LHOOK) CALL DR_HOOK('CANOPY_GRID_UPDATE',1,ZHOOK_HANDLE)
IF(ALL(SB%XZ(:,SB%NLVL)==PZFORC(:))) RETURN
!
IF (PRESENT(OLOG_GRID)) THEN
  GLOG_GRID = OLOG_GRID
ELSE
  GLOG_GRID = .FALSE.
END IF
!
!-------------------------------------------------------------------------------
!
!*    1.  set upper level to forcing height
!         ---------------------------------
!
SB%XZ(:,SB%NLVL) = PZFORC(:)
!
!*    2.  all canopy levels remaining above forcing height are relocated below
!         --------------------------------------------------------------------
!
! determination of levels below forcing height, low enough
!
IF(GLOG_GRID) THEN
  ! Logarithmic grid when blowing snow scheme is used
  ! Lower height for flux is taken at 30 cm above the snow surface
  SB%XZF(:,2) = 0.3
  SB%XZF(:,SB%NLVL) = 0.95*PZFORC(:)
  DO JI=1,KI
    ZDZ = (LOG(SB%XZF(JI,SB%NLVL))-LOG(SB%XZF(JI,2)))/(SB%NLVL-2)
    DO JLAYER=2,SB%NLVL
      SB%XZF(JI,JLAYER) = EXP((JLAYER-2)*ZDZ+LOG(SB%XZF(JI,2)))
      SB%XZ(JI,JLAYER-1) = (SB%XZF(JI,JLAYER)+SB%XZF(JI,JLAYER-1))/2
    END DO
  END DO
ELSE
  ILEVEL=0
  DO JI=1,KI
    DO JLAYER=1,SB%NLVL-1
      IF( PZFORC(JI) > SB%XZF(JI,JLAYER+1) + 0.25 * SB%XDZ(JI,JLAYER) .AND. &
          SB%XZ(JI,JLAYER) < PZFORC(JI) ) ILEVEL(JI,JLAYER) = JLAYER
    ENDDO
    ! determination of latest level from the ones selected before
    IL(JI)=MAXVAL(ILEVEL(JI,1:SB%NLVL-1))
    !
    ICOUNT = SB%NLVL-IL(JI)-1
    !
    !* determination grid top of this level
    ZZTOP = SB%XZF(JI,IL(JI)+1) ! ZZTOP=0 for IL=0
    ZDZ   = 2. * ( SB%XZ(JI,SB%NLVL)-ZZTOP ) / ( 2*ICOUNT+1 )
    DO JLAYER=1,ICOUNT
      SB%XZ(JI,JLAYER+IL(JI)) = ZZTOP + (JLAYER-0.5) * ZDZ
    END DO
  END DO
ENDIF
!
!*    3.  New grid characteristics
!         ------------------------
!
 CALL CANOPY_GRID(KI,SB)
!
!
!*    5.  at least one canopy level in addition to forcing level must be above canopy top
!         -------------------------------------------------------------------------------
!
DO JI=1,KI
  !
  !* tests if the level below forcing height is high enough above canopy
  IF(SB%XZF(JI,SB%NLVL-1) < PH(JI) ) THEN
    !
    !* sets bottom of grid box that is below the forcing level one at canopy height
    !
    SB%XZF(JI,SB%NLVL-1) = PH(JI)
    !
    !* rebuilds vertical grid from the bottom of each grid
    !
    SB%XZ(JI,SB%NLVL-2) = 0.5 * ( SB%XZF(JI,SB%NLVL-2) + SB%XZF(JI,SB%NLVL-1) )
    SB%XZ(JI,SB%NLVL-1) = ( 2.* SB%XZF(JI,SB%NLVL-1) + SB%XZ (JI,SB%NLVL) ) /3.
  END IF
END DO
!
!*    6.  Final grid characteristics
!         --------------------------
!
 CALL CANOPY_GRID(KI,SB)
!
IF (LHOOK) CALL DR_HOOK('CANOPY_GRID_UPDATE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE CANOPY_GRID_UPDATE
