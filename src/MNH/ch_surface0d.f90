!MNH_LIC Copyright 1999-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################################################
      SUBROUTINE CH_SURFACE0D(PTSIMUL, PDTACT, PCONC, KEQ)
!!    ####################################################
!!
!!*** *CH_SURFACE0D*
!!
!!    PURPOSE
!!    -------
!!      calculate surface deposition /emission fluxes
!!
!!    METHOD
!!    ------
!!      The fluxes are calculates after the chemical solver and
!!    are applied to the resulting vector C^*, which means we are doing time
!!    splitting between chemistry and surface fluxes; the Euler implicit
!!    method is applied to achieve the best stability properties:
!!
!!    C(n+1) = ( C^* + DT*F/H ) / ( 1 + DT*v_dep/H )
!!
!!      The deposition velocities will be calculated using Pierre Tulet's
!!    Wesley scheme (iterfaced to CH_SURF_FLUXn)
!!      The different surface parameters will be updated from a namelist
!!    that is continously read from CHCONTROL1.nam (keyword: SURFDATA)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 03/03/99
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!  Philippe Wautelet: 10/01/2019: use newunit argument to open files
!!  Philippe Wautelet: 22/01/2019: use standard FLUSH statement instead of non standard intrinsics
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_OPEN_INPUT
USE MODI_CH_READ_VECTOR
!UPG
!USE MODI_CH_DEPOSITION_FLUX1
!UPG
USE MODI_CH_EMISSION_FLUX0D
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D, ONLY: XDTOUT, XTBEGIN, NVERB
USE MODD_CH_M9_n,    ONLY: CNAMES
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
REAL, DIMENSION(KEQ), INTENT(INOUT) :: PCONC ! concentration vector
!
!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
INTEGER, SAVE :: ILU = -1 ! unit number for IO
INTEGER :: JI
REAL, DIMENSION(KEQ) :: ZDEPO  ! deposition velocity
REAL, DIMENSION(KEQ) :: ZEMIS  ! emission flux after multiplication 
			       ! with emission factor
REAL :: ZEMISFACTOR            ! emission factor
!
REAL          :: ZHEIGHT              ! the height of the box
!
LOGICAL, SAVE :: LSFIRST = .TRUE. 
INTEGER, SAVE :: ZSTNEXTOUT
!
!==========================================================================
! INTERFACE TO CH_DEPOSITION_FLUXn (to be cleaned up later)
!
REAL, DIMENSION(1,1,1,KEQ)            :: ZSVT     ! scalar variables at t
REAL, DIMENSION(1,1)                  :: ZVMOD    ! wind module
REAL, DIMENSION(1,1,1)                :: ZRHODREF ! air density
REAL, DIMENSION(1,1)                  :: ZTS      ! surface temperature
REAL, DIMENSION(1,1)                  :: ZRS      ! stomatal resistance
REAL, DIMENSION(1,1)                  :: ZRESA    ! aerodynamic resistance
!
REAL, DIMENSION(1,1,1)                :: ZUSTAR_PATCH, ZHU_PATCH,ZPSN_PATCH
REAL, DIMENSION(1,1)                  :: ZUSTAR_SEA, &
                                         ZUSTAR_WATER, ZUSTAR_TOWN, &
                                         ZRESA_SEA, ZRESA_WATER, ZRESA_TOWN
!
REAL, DIMENSION(1,1)                  :: ZTHT2D   ! potential temperature at the surface  
REAL, DIMENSION(1,1)                  :: ZPABST2D ! pressure at the surface  
REAL, DIMENSION(1,1)                  :: ZSST, ZCLAY, ZSAND, ZLAND  
REAL, DIMENSION(1,1)                  :: ZZ0VEG, ZVEG, ZLAI   
REAL, DIMENSION(1,1)                  :: ZPSN, ZPSNV, ZPSNG, ZCD, ZHU   
!                       ZCD = drag coefficients for momentum
!                       ZCLAY, ZSAND = percentage of clay and
!                                      sand of the soil
!                       ZLAND = land-surface mask
!                       ZVEG = fraction of vegetation
!                       ZLAI = leaf area index
!                       ZSST = sea surface temperature
!                       ZTS = land surface temperature
!                       ZHU = soil humidity
!                       ZPSN = fraction of the grid covered by snow
!                       ZPSNG = fraction of the bare ground
!                               covered by snow
!                       ZPSNV = fraction of the vegetation
!                               covered by snow
!                       ZZ0VEG = roughness length for momentum
REAL, DIMENSION(1,1,KEQ) :: ZSFSV        ! flux of scalar variables
!==========================================================================
!
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
!*    1.   CALCULATE EMISSION FLUX
!          -----------------------
!
!*     1.1 read emission data and calculate the flux
!
CALL CH_EMISSION_FLUX0D(PTSIMUL, ZEMIS, "CHCONTROL1.nam", 6, NVERB)

!*     1.2 calculate some variables for later use
!
! convert m2 to cm2, days to seconds and Mole to molecules
! ZEMIS(1:KEQ) = (6.022136E23/1E4/86400.)*ZEMIS(1:KEQ)
!
! convert ppp*m/s to molecules/cm2/s
! based on 1 ppp*m/s = (Na/Vmol)*m/s = (6.022136E23/22.414E-3) molecules/m2/s
ZEMIS(1:KEQ) = (6.022136E23/224.14)*ZEMIS(1:KEQ)
!
!*    2.   CALCULATE DEPOSITION FLUXES WITH WESLEY 
!          ---------------------------------------
!
CALL CH_SURFACE0D_SETPARAM
!
!*    2.2  call CH_DEPOSITION_FLUXn
!
ZSVT(1,1,1,:) = PCONC(:) ! scalar variables at t
!
!UPG
!TO BE UPDATE WITH EXTERNALIZATION *********pierre tulet***********
!CALL CH_DEPOSITION_FLUX1(ZSVT, ZSFSV, 6, NVERB,                   &
!                      ZUSTAR_PATCH,ZUSTAR_SEA,ZUSTAR_WATER,ZUSTAR_TOWN,&
!                      ZRESA_SEA, ZRESA_WATER, ZRESA_TOWN,         &
!                      ZHU_PATCH,ZPSN_PATCH)
!UPG
!
! results in  ZSFSV(1,1,:) ! flux of scalar variables (ppp*m/s)
! but we do not use them here, we rather take the deposition velocity 
! directly from the module (variable XVDEPT)
!
! convert deposition velocity from m/s to cm/s
!UPG
!TO BE UPDATE WITH EXTERNALIZATION *********pierre tulet***********
!ZDEPO(:) = 100.*XVDEPT(1,1,:)
!UPG
!
! apply the emission factor
ZEMIS(:) = ZEMIS(1:KEQ) * ZEMISFACTOR
!
!
!*    3.   FLUX CALCULATION
!          ----------------
!
PCONC(:) = (PCONC(:) + PDTACT*ZEMIS(:)/(100.*ZHEIGHT)) &
	 / (1. + PDTACT*ZDEPO(:)/(100.*ZHEIGHT))
!
!* generate a makeplot-type output file
!
IF (LSFIRST) THEN
  ZSTNEXTOUT = XTBEGIN
  OPEN(NEWUNIT = ILU,            &
       FILE    = "SURFACE.out", &
       FORM    = "FORMATTED",   &
       STATUS  = "UNKNOWN"      )
  WRITE(UNIT=ILU, FMT='(A)') 'parameters from the Wesley scheme'
  WRITE(UNIT=ILU, FMT='(A)') 'quick & dirty coding :-)'
  WRITE(UNIT=ILU, FMT='(2I5,A)') 1+2*KEQ+21, KEQ, ' 0 0'
  WRITE(UNIT=ILU, FMT='(A)') 'XTSIMUL'
  DO JI = 1, KEQ
    WRITE (UNIT=ILU, FMT='(2A)') "VDEP-", CNAMES(JI)
  END DO
  DO JI = 1, KEQ
    WRITE (UNIT=ILU, FMT='(2A)') "EMIS-", CNAMES(JI)
  END DO
  WRITE(UNIT=ILU, FMT='(A)') 'VMOD'
  WRITE(UNIT=ILU, FMT='(A)') 'RHODREF'
  WRITE(UNIT=ILU, FMT='(A)') 'RS'
  WRITE(UNIT=ILU, FMT='(A)') 'RESA'
  WRITE(UNIT=ILU, FMT='(A)') 'THT2D'
  WRITE(UNIT=ILU, FMT='(A)') 'PABST2D'
  WRITE(UNIT=ILU, FMT='(A)') 'CD'
  WRITE(UNIT=ILU, FMT='(A)') 'CLAY'
  WRITE(UNIT=ILU, FMT='(A)') 'SAND'
  WRITE(UNIT=ILU, FMT='(A)') 'LAND'
  WRITE(UNIT=ILU, FMT='(A)') 'VEG'
  WRITE(UNIT=ILU, FMT='(A)') 'LAI'
  WRITE(UNIT=ILU, FMT='(A)') 'SST'
  WRITE(UNIT=ILU, FMT='(A)') 'TS'
  WRITE(UNIT=ILU, FMT='(A)') 'HU'
  WRITE(UNIT=ILU, FMT='(A)') 'PSN'
  WRITE(UNIT=ILU, FMT='(A)') 'PSNG'
  WRITE(UNIT=ILU, FMT='(A)') 'PSNV'
  WRITE(UNIT=ILU, FMT='(A)') 'Z0VEG'
  WRITE(UNIT=ILU, FMT='(A)') 'HEIGHT'
  WRITE(UNIT=ILU, FMT='(A)') 'EMISFACTOR'
  WRITE(UNIT=ILU, FMT='(A)') '(5E16.8)'
  LSFIRST = .FALSE.
END IF

IF (PTSIMUL >= ZSTNEXTOUT) THEN 
  ZSTNEXTOUT = ZSTNEXTOUT + XDTOUT
  WRITE(UNIT=ILU, FMT= '(5E16.8)') PTSIMUL, &
        (ZDEPO(JI),JI=1,KEQ), &
        (ZEMIS(JI),JI=1,KEQ), &
        ZVMOD(1,1), &   
        ZRHODREF(1,1,1), &   
        ZRS(1,1), &   
        ZRESA(1,1), &   
        ZTHT2D(1,1), &   
        ZPABST2D(1,1), &   
        ZCD(1,1), &   
        ZCLAY(1,1), &   
        ZSAND(1,1), &   
        ZLAND(1,1), &   
        ZVEG(1,1), &   
        ZLAI(1,1), &   
        ZSST(1,1), &   
        ZTS(1,1), &   
        ZHU(1,1), &   
        ZPSN(1,1), &   
        ZPSNG(1,1), &   
        ZPSNV(1,1), &   
        ZZ0VEG(1,1), &   
        ZHEIGHT, &
        ZEMISFACTOR
  FLUSH(UNIT=ILU)
END IF

RETURN
!
!----------------------------------------------------------------------------
CONTAINS
!
SUBROUTINE CH_SURFACE0D_SETPARAM
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_FM,    ONLY: IO_FILE_CLOSE_ll,IO_FILE_OPEN_ll
!
IMPLICIT NONE
!
! set the input parameters for CH_DEPOSITION_FLUX that are given
! in the 3-D model by MesoNH
REAL, SAVE :: VMOD    = 10.     ! wind module
REAL, SAVE :: RHODREF = 1.122   ! air density
REAL, SAVE :: RS      = 200.    ! stomatal resistance
REAL, SAVE :: RESA    = -999.   ! aerodynamic resistance
REAL, SAVE :: THT2D   = 298.    ! potential temperature at the surface  
REAL, SAVE :: PABST2D = 1000E2  ! pressure at the surface  
REAL, SAVE :: CD      = 0.001   ! drag coefficients for momentum
REAL, SAVE :: CLAY    = 0.5     ! percentage of clay of the soil
REAL, SAVE :: SAND    = 0.5     ! percentage of sand of the soil
REAL, SAVE :: LAND    = 1.0     ! land-surface mask
REAL, SAVE :: VEG     = 0.3     ! fraction of vegetation
REAL, SAVE :: LAI     = 2.5     ! leaf area index
REAL, SAVE :: SST     = 298.    ! sea surface temperature
REAL, SAVE :: TS      = 298.    ! land surface temperature
REAL, SAVE :: HU      = 0.75    ! soil humidity
REAL, SAVE :: PSN     = 0.0     ! fraction of the grid covered by snow
REAL, SAVE :: PSNG    = 0.0     ! fraction of the bare ground covered by snow
REAL, SAVE :: PSNV    = 0.0     ! fraction of the vegetation covered by snow
REAL, SAVE :: Z0VEG   = 0.10    ! roughness length for momentum
REAL, SAVE :: HEIGHT  = 1000.   ! box height in meter
REAL, SAVE :: EMISFACTOR  = 1.  ! emission factor
REAL, SAVE :: TVALID  = 0.      ! duration of validity for the current setting
REAL, SAVE :: TNEXTUPDATE = -1E20 ! time for next update from namelist
INTEGER, SAVE :: ISURFIO        ! I/O channel for surface namelist reading
!
LOGICAL, SAVE :: LSFIRSTCALL = .TRUE. 
REAL, SAVE :: USTAR_PATCH, USTAR_SEA, USTAR_WATER, USTAR_TOWN,&
              RESA_SEA, RESA_WATER, RESA_TOWN,                &
              HU_PATCH, PSN_PATCH
TYPE(TFILEDATA),POINTER,SAVE :: TZFILE => NULL()
!
NAMELIST /NAM_SURF/ VMOD, RHODREF, RS, THT2D, PABST2D, CD,              &
                    CLAY, SAND, LAND, VEG, LAI, SST, TS, HU,            &
                    PSN, PSNG, PSNV, Z0VEG, HEIGHT, EMISFACTOR, TVALID, &
                    USTAR_PATCH, USTAR_SEA, USTAR_WATER, USTAR_TOWN,    &
                    RESA_SEA, RESA_WATER, RESA_TOWN,                    &
                    HU_PATCH, PSN_PATCH
!
IF (LSFIRSTCALL) THEN
!
  LSFIRSTCALL = .FALSE.
  CALL CH_OPEN_INPUT("SURFACE.nam", "SURFDATA", TZFILE, 6, NVERB)
  ISURFIO = TZFILE%NLU
  CALL IO_FILE_CLOSE_ll(TZFILE)
!
END IF
!
IF (PTSIMUL .GE. TNEXTUPDATE) THEN
  PRINT *, "updating surface variables from file ",TRIM(TZFILE%CNAME)
  CALL IO_FILE_OPEN_ll(TZFILE)
  READ(ISURFIO,NAM_SURF)
  CALL IO_FILE_CLOSE_ll(TZFILE)
  ! PRINT *, "current setting is:"
  ! WRITE(*,NAM_SURF)
  TNEXTUPDATE = PTSIMUL + TVALID
  ! PRINT *, "next update of surface parameters at T=", TNEXTUPDATE
END IF
!
!     * calculate drag coefficient and aerodynamic resistance
!
!
IF (LAND .LE. 0.0) THEN
  CD = 4.4E-4 * (VMOD**0.55)  ! Stull, p266
END IF
!
IF (VMOD .LE. 0.0) THEN
  RESA = 1E10
ELSE
  RESA = 1. / (VMOD * CD)
END IF
!
ZVMOD    = VMOD    ! wind module
ZRHODREF = RHODREF ! air density
ZRS      = RS      ! stomatal resistance
ZRESA    = RESA    ! aerodynamic resistance
ZTHT2D   = THT2D   ! potential temperature at the surface  
ZPABST2D = PABST2D ! pressure at the surface  
ZCD      = CD      ! drag coefficients for momentum
ZCLAY    = CLAY    ! percentage of clay of the soil
ZSAND    = SAND    ! percentage of sand of the soil
ZLAND    = LAND    ! land-surface mask
ZVEG     = VEG     ! fraction of vegetation
ZLAI     = LAI     ! leaf area index
ZSST     = SST     ! sea surface temperature
ZTS      = TS      ! land surface temperature
ZHU      = HU      ! soil humidity
ZPSN     = PSN     ! fraction of the grid covered by snow
ZPSNG    = PSNG    ! fraction of the bare ground covered by snow
ZPSNV    = PSNV    ! fraction of the vegetation covered by snow
ZZ0VEG   = Z0VEG   ! roughness length for momentum
ZEMISFACTOR = EMISFACTOR  ! emission factor
ZHEIGHT  = HEIGHT  ! box height in meter
ZUSTAR_PATCH = USTAR_PATCH
ZUSTAR_SEA = USTAR_SEA
ZUSTAR_WATER = USTAR_WATER
ZUSTAR_TOWN = USTAR_TOWN
ZRESA_SEA = RESA_SEA
ZRESA_WATER = RESA_WATER
ZRESA_TOWN = RESA_TOWN
ZHU_PATCH = HU_PATCH
ZPSN_PATCH = PSN_PATCH
!
END SUBROUTINE CH_SURFACE0D_SETPARAM
!
END SUBROUTINE CH_SURFACE0D
