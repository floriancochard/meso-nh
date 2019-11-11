!ORILAM_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!ORILAM_LIC This is part of the ORILAM software governed by the CeCILL-C licence
!ORILAM_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!ORILAM_LIC for details.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!   ########################
     MODULE MODI_CH_AER_EQSAM
!!   ########################
!!
INTERFACE
!!
SUBROUTINE CH_AER_EQSAM(PAER,PRH, PPRESSURE, PTEMP)
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(INOUT) :: PAER
REAL, DIMENSION(:),   INTENT(IN)    :: PRH, PPRESSURE, PTEMP
END SUBROUTINE CH_AER_EQSAM
!!
END INTERFACE
!!
END MODULE MODI_CH_AER_EQSAM
!!
!!    ##########################################################
      SUBROUTINE CH_AER_EQSAM(PAER,PRH, PPRESSURE, PTEMP )
!!    ##########################################################
!!
!!    PURPOSE
!!    -------
!!
!!    Interface between ORILAM and EQSAM for calculate the aerosol chemical speciation and water content.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    P. Tulet ( Meteo France / GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!
!      interface
!      ---------
!      call  eqsam_v03d(yi,yo,nca,nco,iopt,loop,imax,ipunit,in)
!
!      yi = input  array (imax, nca)
!      yo = output array (imax, nco)
!      imax = max loop (e.g. time steps)
!      nca >= 11
!      nc0 >= 35
!      iopt = 1 metastable 
!      iopt = 2 solids 
!      iopt = 3 hysteresis (metastable/solids) for online calculations
!      iopt = 31 hysteresis lower branch 
!      iopt = 32 hysteresis upper branch 
!      ipunit = I/O unit (can be skipped)
!      in = array        (can be skipped)
!         
!      method
!      ------
!      equilibrium / internal mixture assumption / aw=rh
!      System: NH3,NH4+/H2SO4+,HSO4-,SO4--/HNO3,NO3-, HCl,Cl-/Na+, H2O 
!              (K+,Ca++,Mg++)
!      external
!      --------
!      program    eqmd.f90    (driver only needed for the box model version)
!      subroutine gribio.f90  (provides diagnostics output in grib/binary/ascii format)
!      
!      references
!      ---------
!      Swen Metzger Ph.D Thesis, University Utrecht, 2000.
!         http://www.library.uu.nl/digiarchief/dip/diss/1930853/inhoud.htm
!
!      Metzger, S. M., F. J. Dentener, J. Lelieveld, and S. N. Pandis, 
!         GAS/AEROSOL PARTITIONING I: A COMPUTATIONALLY EFFICIENT MODEL, 
!         J Geophys. Res., 107, D16, 10.1029/2001JD001102, 2002
!         http://www.agu.org/journals/jd/jd0216/2001JD001102/index.html
!      Metzger, S. M., F. J. Dentener, A. Jeuken, and M. Krol, J. Lelieveld, 
!         GAS/AEROSOL PARTITIONING II: GLOBAL MODELING RESULTS, 
!         J Geophys. Res., 107, D16, 10.1029/2001JD001103, 2002.
!         http://www.agu.org/journals/jd/jd0216/2001JD001103/index.html
!_________________________________________________________________________________________________
!
!***************************************************************
!!
!!   EXTERNAL
!!   -------
!!
USE MODI_eqsam_v03d_sub
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
IMPLICIT NONE
!!...........ARGUMENTS and their descriptions
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PAER
REAL, DIMENSION(:),   INTENT(IN)    :: PRH, PPRESSURE, PTEMP
!
! PAER(:,1) :: H2SO4    in micrograms / m**3
! PAER(:,2) :: NH3(g)   in micrograms / m**3
! PAER(:,3) :: HNO3(g)  in micrograms / m**3
! PAER(:,4) :: H2O(a)   in micrograms / m**3
! PAER(:,5) :: NO3(a)   in micrograms / m**3
! PAER(:,6) :: NH4(a)   in micrograms / m**3
!
!...........PARAMETERS and their descriptions:

REAL, PARAMETER ::  ZMWH2O = 18.0           ! molecular weight for water      
REAL, PARAMETER ::  ZMWNO3 = 62.0049        ! molecular weight for NO3
REAL, PARAMETER ::  ZMWHNO3 = 63.01287      ! molecular weight for HNO3
REAL, PARAMETER ::  ZMWSO4 = 96.0576        ! molecular weight for SO4
REAL, PARAMETER ::  ZMH2SO4 = 98.07354      ! molecular weight for H2SO4
REAL, PARAMETER ::  ZMWNH3 = 17.03061       ! molecular weight for NH3
REAL, PARAMETER ::  ZMWNH4 = 18.03858       ! molecular weight for NH4
REAL, PARAMETER ::  ZMWAIR = 28.964         ! molecular weight for AIR

INTEGER :: NCA,NCO,IOPT

REAL,   ALLOCATABLE, DIMENSION(:,:)  :: ZYI, ZYO
!
!-----------------------------------------------------------------------------
NCA=11  ! Number of INPUT
NCO=35  ! Number of OUTPUT
IOPT=1  ! Options for metastable, solids, hysteresis (metastable/solids) 
!
ALLOCATE(ZYI(SIZE(PAER,1), NCA))
ALLOCATE(ZYO(SIZE(PAER,1), NCO))
ZYI(:,:) =0.
ZYO(:,:) =0.
 
      ZYI(:,1)=PTEMP(:)                                   ! T                      [K]
      ZYI(:,2)=PRH(:)                                     ! RH                     [0-1]
      ZYI(:,3)=PAER(:,2) / ZMWNH3 + PAER(:,6) / ZMWNH4    ! NH3  (g) + NH4+  (p)   [umol/m^3 air]
      ZYI(:,4)=PAER(:,1) / ZMH2SO4                        ! H2SO4    + SO4-- (p)   [umol/m^3 air]
      ZYI(:,5)=PAER(:,5) / ZMWNO3 + PAER(:,3) / ZMWHNO3   ! HNO3 (g) + NO3-  (p)   [umol/m^3 air]
      ZYI(:,6)= 0. ! Na+ (ss  + xsod) (a)   [umol/m^3 air]
      ZYI(:,7)= 0. ! HCl  (g) + Cl-   (p)   [umol/m^3 air]
      ZYI(:,8)= 0. ! K+   (p) from Dust     [umol/m^3 air]
      ZYI(:,9)= 0. ! Ca++ (p) from Dust     [umol/m^3 air]
      ZYI(:,10)=0. ! Mg++ (p) from Dust     [umol/m^3 air]
      ZYI(:,11)= PPRESSURE(:)*1E-2 !  p     [hPa]

      CALL eqsam_v03d_sub(ZYI,ZYO, NCA, NCO, IOPT, SIZE(PAER,1), SIZE(PAER,1)) 
      !PAER(:,1) = ZYO(:,21) * ZMH2SO4
      PAER(:,2) = ZYO(:,10) * ZMWNH3
      PAER(:,3) = ZYO(:,9)  * ZMWHNO3
      PAER(:,4) = ZYO(:,12) ! * ZMWH2O note H2O is in micro-g/m3
      PAER(:,5) = ZYO(:,20) * ZMWNO3 
      PAER(:,6) = ZYO(:,19) * ZMWNH4

DEALLOCATE(ZYO)
DEALLOCATE(ZYI)
!
END SUBROUTINE CH_AER_EQSAM

