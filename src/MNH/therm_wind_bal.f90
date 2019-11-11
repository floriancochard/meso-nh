!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_THERM_WIND_BAL
!     ##########################
INTERFACE
!
      SUBROUTINE THERM_WIND_BAL(PZMAX,PCORIO,PR_M,PRADMAX,  &
                                PZ,PTH,                          &
                                PTH_DEV)
!
REAL, INTENT(IN)        :: PZMAX! Altitude for which the wind vanishes (m)
                                !for the call to HOLLAND_VT
REAL, INTENT(IN)        :: PCORIO ! Coriolis parameter at vortex center
REAL, INTENT(IN)        :: PR_M   ! Distance from the considered
                                  !point to the center of circulation
REAL, INTENT(IN)        :: PRADMAX ! Outermost boundary of the vortex
REAL, DIMENSION(:), INTENT(IN):: PZ     ! altitude of points of reference profil
REAL, DIMENSION(:), INTENT(IN):: PTH    ! THETA reference state profil
REAL, DIMENSION(:),INTENT(OUT):: PTH_DEV ! Potential temperature deviation 
                                         !from thermal wind balance equation
!
END SUBROUTINE THERM_WIND_BAL
END INTERFACE
!
END MODULE MODI_THERM_WIND_BAL
!
!
!
!     ###################################################################
      SUBROUTINE THERM_WIND_BAL(PZMAX,PCORIO,PR_M,PRADMAX,  &
                                PZ,PTH,                          &
                                PTH_DEV)
!     ###################################################################
!
!!****  *THERM_WIND_BAL* - This subroutine computes a potential temperature
!!                         perturbation which balances a symmetric
!!                          tropical like-vortex
!!                           
!!
!!    PURPOSE
!!    -------
!         This subroutine calculates............ (TO BE COMPLETED)
!
!
!
!!**  METHOD
!!    ------
!!           ..................................................
!!
!!          TO BE COMPLETED
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!     **??** subroutine THERM_WIND_BAL
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      F. Roux               * L.A  *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              10/12/03
!!      F. Roux               22/03/04 computation of Vt according zeta 
!!      Modification 01/02/08 (D.Barbary) Use NAM_HURR_PARAM for Holland's parameters
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY: XG
USE MODD_GRID_n, ONLY: XXHAT,XYHAT
USE MODI_HOLLAND_VT  ! module which compute Vt(r,z)
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, INTENT(IN)        :: PZMAX! Altitude for which the wind vanishes
                                !for the call to HOLLAND_VT
REAL, INTENT(IN)        :: PCORIO! Coriolis parameter at vortex center
REAL, INTENT(IN)        :: PR_M  ! Distance from the considered
                                 !point to the center of circulation
REAL, INTENT(IN)        :: PRADMAX ! Outermost boundary of the vortex
REAL, DIMENSION(:), INTENT(IN):: PZ     ! altitude of points of reference profil
REAL, DIMENSION(:), INTENT(IN):: PTH    ! THETA reference state profil
REAL, DIMENSION(:),INTENT(OUT):: PTH_DEV ! Potential temperature deviation 
                                         !from thermal wind balance equation
!
!*       0.2   Declarations of local variables
!
REAL                          :: ZDELTAX,ZDELTAY,ZDELTAR 
REAL                          :: ZRIR,ZRIR_M 
REAL,DIMENSION(SIZE(PTH_DEV)) :: ZVT_M,ZVT_TOP,ZVT_BOTOM
REAL,DIMENSION(SIZE(PTH_DEV)) :: ZZETA_TOP,ZZETA_BOTOM
REAL,DIMENSION(SIZE(PTH_DEV)) :: ZDVT_DZETA,ZDTH_DR
INTEGER                       :: INR,JR,ILU
!
!-------------------------------------------------------------------------------
!
!*    1. INITIALIZATIONS
!        ---------------
!
ZDELTAX = XXHAT(3) - XXHAT(2)
ZDELTAY = XYHAT(3) - XYHAT(2)
ZDELTAR = MIN(MAX(ZDELTAX,ZDELTAY),5000.)
!
INR = (PRADMAX - PR_M)/ZDELTAR
ILU = SIZE(PTH_DEV)
PTH_DEV(:) = 0.
!
!-------------------------------------------------------------------------------
!
!*    2. COMPUTE THE POTENTIAL TEMPERATURE DEVIATION AT A GIVEN R RADIUS 
!             AND A GIVEN Z ALTITUDE
!        Suppose une integration radiale dans la direction I afin d avoir
!        coriolis constant dans la boucle radiale
!
DO JR = 0, INR-1
  ! Starting integration from Rmax to the given Rm
  ZRIR = PRADMAX - JR * ZDELTAR
  ZRIR_M = ZRIR - ZDELTAR/2.
  !
  !
  ! Compute Vt_m at the point rir_m and zeta_m 
  CALL HOLLAND_VT(PZMAX,PZ(1:ILU),ZRIR_M,PCORIO,ZVT_M(1:ILU))
  !
  ! Compute Vt_m at the point rir_m and zeta_top
  ZVT_TOP(2:ILU-1) = ZVT_M(3:ILU)
  ZZETA_TOP(2:ILU-1) = PZ(3:ILU)
  !
  ! Compute Vt_m at the point rir_m and zeta_botom
  ZVT_BOTOM(2:ILU-1) = ZVT_M(1:ILU-2)
  ZZETA_BOTOM(2:ILU-1) = PZ(1:ILU-2)
  !
  ! Finally, this last execution computes the potential temperature
  ! perturbation at the point rir_m and z_m 
  ZDVT_DZETA(2:ILU-1) =  ( ZVT_TOP(2:ILU-1)   - ZVT_BOTOM(2:ILU-1)    ) &
                        /( ZZETA_TOP(2:ILU-1) - ZZETA_BOTOM(2:ILU-1)  )
  !
  ZDTH_DR(2:ILU-1) =  (PTH(2:ILU-1) / XG)                     &
                    * (2. * ZVT_M(2:ILU-1) / ZRIR_M + PCORIO) &
                    * ZDVT_DZETA(2:ILU-1)                     ! dVT/dzeta term 
  !
  PTH_DEV(2:ILU-1) = PTH_DEV(2:ILU-1) - ZDTH_DR(2:ILU-1) * ZDELTAR
  !
END DO  
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE THERM_WIND_BAL
