!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_CONVECTION
!     ######################
!
INTERFACE
!
    SUBROUTINE CONVECTION( PDTCONV, HDCONV, HSCONV, OREFRESH_ALL, ODOWN, KICE, &
                           OSETTADJ, PTADJD, PTADJS, ODIAGCONV,                &
                           KENSM,                                              &
                           PPABST, PZZ, PDXDY,                                 &
                           PTHT, PRVT, PRCT, PRIT, PUT, PVT, PWT, PTKECLS,     &
                           KCOUNT, PTHTEN, PRVTEN, PRCTEN, PRITEN,             &
                           PPRTEN, PPRSTEN,                                    &
                           PUMF, PDMF, PMF, PPRLFLX, PPRSFLX, PCAPE, KCLTOP, KCLBAS,&
                           OCHTRANS, PCH1, PCH1TEN,                            &
                           OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX,            &
                           ODUST, OSALT,                                       &
                           PRHODREF, PIC_RATE, PCG_RATE                        )
!
REAL,                     INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                 ! calls of the deep convection
                                                 ! scheme
CHARACTER(LEN=4),         INTENT(IN) :: HDCONV   ! kind of deep convection
CHARACTER(LEN=4),         INTENT(IN) :: HSCONV   ! kind of shallow convection
LOGICAL,                  INTENT(IN) :: OREFRESH_ALL ! refresh or not all 
                                                 ! tendencies  at every call
LOGICAL,                  INTENT(IN) :: ODOWN    ! take or not convective
                                                 ! downdrafts into account
INTEGER,                  INTENT(IN) :: KICE     ! flag for ice ( 1 = yes, 
                                                 !                0 = no ice )
LOGICAL,                  INTENT(IN) :: OSETTADJ ! logical to set convective
                                                 ! adjustment time by user
REAL,                     INTENT(IN) :: PTADJD   ! user defined deep adjustment time
REAL,                     INTENT(IN) :: PTADJS   ! user defined shal. adjustment time
LOGICAL,                  INTENT(IN) :: ODIAGCONV! activate or not convection diagnostics
INTEGER,                  INTENT(IN) :: KENSM    ! number of additional deep convection calls
                                                 ! for ensemble (presently limited to 3)
                                                 ! KENSM=0 corresponds to base run with
                                                 ! 1 deep and 1 shallow call
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT     ! grid scale theta at time t  
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRVT     ! grid scale water vapor at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRCT     ! grid scale r_c at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRIT     ! grid scale r_i at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUT      ! grid scale horiz. wind u at t 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PVT      ! grid scale horiz. wind v at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PWT      ! grid scale vertical velocity at t
REAL, DIMENSION(:,:),     INTENT(IN) :: PTKECLS  ! TKE CLS
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABST   ! grid scale pressure at t 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ      ! height of model layer (m) 
REAL, DIMENSION(:,:),     INTENT(IN) :: PDXDY    ! grid area (m2)
!   
INTEGER, DIMENSION(:,:), INTENT(INOUT):: KCOUNT  ! convective counter (recompute
                                                 ! tendency or keep it)
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PTHTEN   ! convective theta tendency (K/s)
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRVTEN   ! convective r_v tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRCTEN   ! convective r_c tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRITEN   ! convective r_i tendency (1/s)
REAL, DIMENSION(:,:),   INTENT(INOUT):: PPRTEN   ! total (liquid+solid) surf 
                                                 ! precipitation tendency (m/s)
REAL, DIMENSION(:,:),   INTENT(INOUT):: PPRSTEN  ! solid surf 
                                                 ! precipitation tendency (m/s)
!
! Tracers:
LOGICAL,                  INTENT(IN) :: OCHTRANS ! flag to compute convective
                                                 ! transport for chemical tracer
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PCH1     ! grid scale chemical species
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PCH1TEN! chemical convective tendency
                                                 ! (1/s)
! Chemical Tracers:
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODREF   ! grid scale density
LOGICAL,                INTENT(IN) :: OUSECHEM, OCH_CONV_SCAV
LOGICAL,                INTENT(IN) :: OCH_CONV_LINOX, ODUST, OSALT
REAL, DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
! Diagnostic variables:
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PUMF   ! updraft mass flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PDMF   ! downdraft mass flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PMF    ! convective mass flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PPRSFLX! solid precip flux
REAL, DIMENSION(:,:),   INTENT(INOUT)  :: PCAPE  ! CAPE (Joule)
INTEGER, DIMENSION(:,:),INTENT(INOUT)  :: KCLTOP ! cloud top level
INTEGER, DIMENSION(:,:),INTENT(INOUT)  :: KCLBAS ! cloud base level
                                                 ! they are given a value of
                                                 ! 0 if no convection
!
END SUBROUTINE CONVECTION
!
END INTERFACE
!
END MODULE MODI_CONVECTION
!
!   ############################################################################
    SUBROUTINE CONVECTION( PDTCONV, HDCONV, HSCONV, OREFRESH_ALL, ODOWN, KICE, &
                           OSETTADJ, PTADJD, PTADJS, ODIAGCONV,                &
                           KENSM,                                              &
                           PPABST, PZZ, PDXDY,                                 &
                           PTHT, PRVT, PRCT, PRIT, PUT, PVT, PWT, PTKECLS,     &
                           KCOUNT, PTHTEN, PRVTEN, PRCTEN, PRITEN,             &
                           PPRTEN, PPRSTEN,                                    &
                           PUMF, PDMF, PMF, PPRLFLX, PPRSFLX, PCAPE, KCLTOP, KCLBAS,&
                           OCHTRANS, PCH1, PCH1TEN,                            &
                           OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX,            &
                           ODUST, OSALT,                                       &
                           PRHODREF, PIC_RATE, PCG_RATE                        )
!   ############################################################################
!
!!**** Interface routine to the fast MNH convection code developed for ARPEGE
!!     having a structure typical for operational routines
!!
!!
!!    PURPOSE
!!    -------
!!      The routine interfaces the MNH convection code as developed for operational
!!      forecast models like ARPEGE or HIRLAM with the typical MNH array structure
!!      Calls the deep and/or shallow convection routine
!!
!!
!!**  METHOD
!!    ------
!!      We essentially just have to skip 3D tables in 2D tables and vice versa
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!    DEEP_CONVECTION
!!    SHALLOW_CONVECTION
!!    INI_CONVPAR, INI_CONVPAR_E1
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_PARAMETERS
!!          JPHEXT, JPVEXT       ! extra levels on the horizontal
!!                               ! and vertical boundaries
!!
!!      Module MODI_SHUMAN
!!          MZF, MXF, MYF        ! vertical and horizontal average operators
!!         
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/12/98 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CST
USE MODI_SHUMAN
USE MODI_DEEP_CONVECTION
USE MODI_SHALLOW_CONVECTION
USE MODI_INI_CONVPAR
USE MODI_INI_CONVPAR_E1
USE MODI_INI_CONVPAR_SHAL

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,                     INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                 ! calls of the deep convection
                                                 ! scheme
CHARACTER(LEN=4),         INTENT(IN) :: HDCONV   ! kind of deep convection
CHARACTER(LEN=4),         INTENT(IN) :: HSCONV   ! kind of shallow convection
LOGICAL,                  INTENT(IN) :: OREFRESH_ALL ! refresh or not all 
                                                 ! tendencies  at every call
LOGICAL,                  INTENT(IN) :: ODOWN    ! take or not convective
                                                 ! downdrafts into account
INTEGER,                  INTENT(IN) :: KICE     ! flag for ice ( 1 = yes, 
                                                 !                0 = no ice )
LOGICAL,                  INTENT(IN) :: OSETTADJ ! logical to set convective
                                                 ! adjustment time by user
REAL,                     INTENT(IN) :: PTADJD   ! user defined deep adjustment time
REAL,                     INTENT(IN) :: PTADJS   ! user defined shal. adjustment time
LOGICAL,                  INTENT(IN) :: ODIAGCONV! activate or not convection diagnostics
INTEGER,                  INTENT(IN) :: KENSM    ! number of additional deep convection calls
                                                 ! for ensemble (presently limited to 3)
                                                 ! KENSM=0 corresponds to base run 
                                                 ! 1 deep + 1 shallow call 
!
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PTHT     ! grid scale theta at time t  
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRVT     ! grid scale water vapor at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRCT     ! grid scale r_c at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRIT     ! grid scale r_i at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PUT      ! grid scale horiz. wind u at t 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PVT      ! grid scale horiz. wind v at t
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PWT      ! grid scale vertical velocity at t
REAL, DIMENSION(:,:),     INTENT(IN) :: PTKECLS  ! TKE CLS
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PPABST   ! grid scale pressure at t 
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PZZ      ! height of model layer (m) 
REAL, DIMENSION(:,:),     INTENT(IN) :: PDXDY    ! grid area (m2)
!   
INTEGER, DIMENSION(:,:), INTENT(INOUT):: KCOUNT  ! convective counter(recompute
                                                 ! tendency or keep it
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PTHTEN   ! convective theta tendency (K/s)
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRVTEN   ! convective r_v tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRCTEN   ! convective r_c tendency (1/s)
REAL, DIMENSION(:,:,:), INTENT(INOUT):: PRITEN   ! convective r_i tendency (1/s)
REAL, DIMENSION(:,:),   INTENT(INOUT):: PPRTEN   ! total surf precipitation tendency (m/s)
REAL, DIMENSION(:,:),   INTENT(INOUT):: PPRSTEN  ! solid surf 
                                                 ! precipitation tendency (m/s)
!
! Tracers:
LOGICAL,                  INTENT(IN) :: OCHTRANS ! flag to compute convective
                                                 ! transport for chemical tracer
REAL, DIMENSION(:,:,:,:), INTENT(IN) :: PCH1     ! grid scale chemical species
REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PCH1TEN! chemical convective tendency
                                                 ! (1/s)
! Chemical Tracers:
REAL, DIMENSION(:,:,:), INTENT(IN) :: PRHODREF   ! grid scale density
LOGICAL,                INTENT(IN) :: OUSECHEM, OCH_CONV_SCAV
LOGICAL,                INTENT(IN) :: OCH_CONV_LINOX, ODUST, OSALT
REAL, DIMENSION(:,:), OPTIONAL,  INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(:,:), OPTIONAL,  INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
!					 
! Diagnostic variables:
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PUMF   ! updraft mass flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PDMF   ! downdraft mass flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PMF    ! convective mass flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PPRSFLX! solid precip flux
REAL, DIMENSION(:,:),   INTENT(INOUT)  :: PCAPE  ! CAPE (J/kg)
INTEGER, DIMENSION(:,:),INTENT(INOUT)  :: KCLTOP ! cloud top level
INTEGER, DIMENSION(:,:),INTENT(INOUT)  :: KCLBAS ! cloud base level
                                                 ! they are given a value of
                                                 ! 0 if no convection
!						 
!
!*       0.2   Declarations of local variables :
!
INTEGER  :: IIB, IIE, IJB, IJE ! horizontal loop bounds
INTEGER  :: IKB, IKE           ! vertical loop bounds
INTEGER  :: IIU, IJU, IKU, ILON, ICH1    ! dimension
INTEGER  :: JI, JJ, JK, JN, IIND ! loop index
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))  :: ZZZ, ZUT, ZVT,&
                                                         ZWT, ZWORK
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPABS, ZZ, ZU, ZV, ZW,           &
                                      ZTT, ZPI, ZRV, ZRC, ZRI,         &
                                      ZTTEN, ZRVTEN, ZRCTEN, ZRITEN,   &
                                      ZPRLFLX, ZPRSFLX,  ZUMF, ZDMF 
REAL, DIMENSION(:), ALLOCATABLE    :: ZTKECLS                                     
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1, ZCH1TEN
REAL, DIMENSION(:), ALLOCATABLE    :: ZDXDY, ZPRLTEN, ZPRSTEN, ZCAPE, ZTIMEC
INTEGER, DIMENSION(:), ALLOCATABLE :: ICOUNT, ICLBAS, ICLTOP
INTEGER                            :: IIDIA, IFDIA, IBDIA, ITDIA
!
! special for shallow convection
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,  &
                                      ZUMFS
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1TENS
INTEGER, DIMENSION(:), ALLOCATABLE :: ICLBASS, ICLTOPS
REAL                               :: ZRDOCP
! Chemical tracer (unit conversion)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRHODREF
! lightning IC and CG
REAL, DIMENSION(:), ALLOCATABLE    :: ZIC_RATE, ZCG_RATE
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZIC_RATEE, ZCG_RATEE
!
!*       0.3   Declarations of additional Ensemble fields:
!
INTEGER                                :: IENS     ! number of allowed additional ensemble members
                                                   ! baserun: 2 (1 deep + 1 shallow)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZTTENE   ! convective temperat. tendency (K/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZRVTENE  ! convective r_v tendency (1/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZRCTENE  ! convective r_c tendency (1/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZRITENE  ! convective r_i tendency (1/s)
REAL, DIMENSION(:,:),   ALLOCATABLE    :: ZPRLTENE ! liquid surf precipitation tendency (m/s)
REAL, DIMENSION(:,:),   ALLOCATABLE    :: ZPRSTENE ! solid surf precipitation tendency (m/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZUMFE    ! updraft mass flux   (kg/s m2)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZDMFE    ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZPRLFLXE ! liquid precip flux  (m/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE    :: ZPRSFLXE ! solid precip flux   (m/s)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZCH1TENE ! chemical convective tendency
REAL, DIMENSION(:),     ALLOCATABLE    :: ZEDUMMY  ! field not to be recomputed by ensemble
INTEGER, DIMENSION(:),  ALLOCATABLE    :: IEDUMMY  ! field not to be recomputed by ensemble
REAL, DIMENSION(:),     ALLOCATABLE    :: ZWEIGHT  ! weighting factor for ensemble members
REAL                                   :: ZSUM     ! sum of weighting factors
!
!-------------------------------------------------------------------------------------------
!
IIU = SIZE(PTHT,1)
IJU = SIZE(PTHT,2)
IKU = SIZE(PTHT,3)
ICH1= SIZE(PCH1,4)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
ILON = (IIE-IIB+1) * (IJE-IJB+1)
!
IIDIA = 1        ! start index of horizontal convection computations
IFDIA = ILON     ! end   index of horizontal convection computations
IBDIA = IKB      ! start index of vertical convection computations
ITDIA = 1+JPVEXT ! limit vertical convection computations
                 ! to IKU + 1 - ITDIA
!
ZRDOCP = XRD / XCPD
!
!
!*       0.5  Allocate 2D (horizontal, vertical) arrays
!             -----------------------------------------
!
ALLOCATE( ZPABS(ILON,IKU) )
ALLOCATE( ZZ(ILON,IKU) ) 
ALLOCATE( ZTT(ILON,IKU) )
ALLOCATE( ZRV(ILON,IKU) ) 
ALLOCATE( ZPI(ILON,IKU) )
ALLOCATE( ZRC(ILON,IKU) )
ALLOCATE( ZRI(ILON,IKU) ) 
ALLOCATE( ZU(ILON,IKU) )
ALLOCATE( ZV(ILON,IKU) ) 
ALLOCATE( ZW(ILON,IKU) )
ALLOCATE( ZTKECLS(ILON) )
ALLOCATE( ZDXDY(ILON) )
ALLOCATE( ZTTEN(ILON,IKU) )
ALLOCATE( ZRVTEN(ILON,IKU) ) 
ALLOCATE( ZRCTEN(ILON,IKU) )
ALLOCATE( ZRITEN(ILON,IKU) ) 
ALLOCATE( ZPRLTEN(ILON) )
ALLOCATE( ZPRSTEN(ILON) ) 
ALLOCATE( ZCH1(ILON,IKU,ICH1) )
ALLOCATE( ZCH1TEN(ILON,IKU,ICH1) ) 
ZCH1TEN(:,:,:)=0.
ALLOCATE( ZPRLFLX(ILON,IKU))
ALLOCATE( ZPRSFLX(ILON,IKU) ) 
ALLOCATE( ZUMF(ILON,IKU) )
ALLOCATE( ZDMF(ILON,IKU) ) 
ALLOCATE( ICLBAS(ILON) )
ALLOCATE( ICLTOP(ILON) ) 
ALLOCATE( ZCAPE(ILON) )
ALLOCATE( ICOUNT(ILON) )
ALLOCATE( ZTIMEC(ILON) )
ALLOCATE( ZRHODREF(ILON,IKU) ) ;  ZRHODREF = 0.0
!
!
ALLOCATE( ZIC_RATE(ILON) ) ; ZIC_RATE  = 0.0
ALLOCATE( ZCG_RATE(ILON) ) ; ZCG_RATE  = 0.0
!
ALLOCATE( ZTTENS(ILON,IKU) )
ALLOCATE( ZRVTENS(ILON,IKU) ) 
ALLOCATE( ZRCTENS(ILON,IKU) )
ALLOCATE( ZRITENS(ILON,IKU) ) 
ALLOCATE( ZCH1TENS(ILON,IKU,ICH1) ) 
ZCH1TENS(:,:,:)=0.
ALLOCATE( ZUMFS(ILON,IKU) )
ALLOCATE( ICLBASS(ILON) )
ALLOCATE( ICLTOPS(ILON) )
!
!* Allocation of additional ensemble members
IENS = MIN( KENSM, 3 )
IF ( IENS > 0 ) THEN
     ALLOCATE( ZTTENE(ILON,IKU,IENS) )
     ALLOCATE( ZRVTENE(ILON,IKU,IENS) )
     ALLOCATE( ZRCTENE(ILON,IKU,IENS) )
     ALLOCATE( ZRITENE(ILON,IKU,IENS) )
     ALLOCATE( ZUMFE(ILON,IKU,IENS) )
     ALLOCATE( ZDMFE(ILON,IKU,IENS) )
     ALLOCATE( ZCH1TENE(ILON,IKU,ICH1,IENS) )
     ZCH1TENE(:,:,:,:)=0.
     ALLOCATE( ZPRLFLXE(ILON,IKU,IENS) )
     ALLOCATE( ZPRSFLXE(ILON,IKU,IENS) )
     ALLOCATE( ZPRLTENE(ILON,IENS) )
     ALLOCATE( ZPRSTENE(ILON,IENS) )
     ALLOCATE( ZEDUMMY(ILON) )
     ALLOCATE( IEDUMMY(ILON) )
     ALLOCATE( ZWEIGHT(IENS) )
     IEDUMMY(:) = 0
     ALLOCATE( ZIC_RATEE(ILON,IENS) ) ; ZIC_RATEE = 0.0
     ALLOCATE( ZCG_RATEE(ILON,IENS) ) ; ZCG_RATEE = 0.0
END IF
!
!----------------------------------------------------------------------------
!
!
!*       1.  Center all fields on thermo levels
!            ----------------------------------
!
ZWORK(:,:,:) = MZF(1,IKU,1, PZZ(:,:,:) ) 
ZZZ(:,:,:)   = ZWORK(:,:,:) 
ZWORK(:,:,:) = MZF(1,IKU,1, PWT(:,:,:) ) 
ZWT(:,:,:)   = ZWORK(:,:,:)
ZWORK(:,:,:) = MXF( PUT(:,:,:) ) 
ZUT(:,:,:)   = ZWORK(:,:,:)
ZWORK(:,:,:) = MYF( PVT(:,:,:) ) 
ZVT(:,:,:)   = ZWORK(:,:,:)
!
!----------------------------------------------------------------------------
!
!
!*       3.  Skip 3D arrays in 2D arrays as used in ARPEGE, ECMWF
!            ----------------------------------------------------
!
DO JK = 1,IKU
  IIND = 1
    DO JI = IIB,IIE
    DO JJ = IJB,IJE
       ZPABS(IIND,JK) = PPABST(JI,JJ,JK)
       ZZ(IIND,JK)    = ZZZ(JI,JJ,JK)
       ZU(IIND,JK)    = ZUT(JI,JJ,JK)
       ZV(IIND,JK)    = ZVT(JI,JJ,JK)
       ZW(IIND,JK)    = ZWT(JI,JJ,JK)
       ZPI(IIND,JK)   = ( XP00 / ZPABS(IIND,JK) ) ** ZRDOCP
       ZTT(IIND,JK)   = PTHT(JI,JJ,JK) / ZPI(IIND,JK)
       ZRV(IIND,JK)   = PRVT(JI,JJ,JK)
       ZRC(IIND,JK)   = PRCT(JI,JJ,JK)
       ZRI(IIND,JK)   = PRIT(JI,JJ,JK)
       IIND = IIND + 1
    END DO
  END DO
END DO
IIND = 1
DO JI = IIB,IIE
  DO JJ = IJB,IJE
     ZDXDY(IIND)    = PDXDY(JI,JJ)
     ZTKECLS(IIND)  = PTKECLS(JI,JJ)
     ICOUNT(IIND)   = KCOUNT(JI,JJ)
     IIND = IIND + 1
  END DO
END DO
ZCH1(:,:,:)=0.
IF ( OCHTRANS ) THEN
  DO JK = 1,IKU
    IIND = 1
    DO JI = IIB,IIE
      DO JJ = IJB,IJE
        ZCH1(IIND,JK,:) = PCH1(JI,JJ,JK,:)
        IIND = IIND + 1
      END DO
    END DO
  END DO
END IF
IF ( OCH_CONV_SCAV .OR. OCH_CONV_LINOX) THEN
  DO JK = 1,IKU
    IIND = 1
    DO JI = IIB,IIE
      DO JJ = IJB,IJE
        ZRHODREF(IIND,JK) = PRHODREF(JI,JJ,JK)
        IIND = IIND + 1
      END DO
    END DO
  END DO
END IF
!
!*       4.a  Call deep convection routine
!             ----------------------------
!
IF ( HDCONV == 'KAFR' ) THEN
!
! 1. Base version                                                                                  !  
  CALL INI_CONVPAR
!
  IF ( OSETTADJ ) ZTIMEC(:) = PTADJD
!
  CALL DEEP_CONVECTION( ILON, IKU, IIDIA, IFDIA, IBDIA, ITDIA,          &
                        PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,   & 
                        ZPABS, ZZ, ZDXDY, ZTIMEC,                       &
                        ZTT, ZRV, ZRC, ZRI, ZU, ZV, ZW,                 &
                        ICOUNT, ZTTEN, ZRVTEN, ZRCTEN, ZRITEN,          &
                        ZPRLTEN, ZPRSTEN,                               &
                        ICLTOP, ICLBAS, ZPRLFLX, ZPRSFLX,               & 
                        ZUMF, ZDMF, ZCAPE,                              &
                        OCHTRANS, ICH1, ZCH1, ZCH1TEN,                  &
                        OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX,        &
                        ODUST, OSALT,                                   &
                        ZRHODREF, ZIC_RATE, ZCG_RATE                    )

!
!  2. Additional Ensemble members
!
    IF ( IENS > 0 ) THEN
!
    CALL INI_CONVPAR_E1
!
!* first member - changes in MODD_CONVPAR (cloud radius of 500 m)
!
    CALL DEEP_CONVECTION( ILON, IKU, IIDIA, IFDIA, IBDIA, ITDIA,                                 &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,                          &
                          ZPABS,  ZZ, ZDXDY, ZTIMEC,                                             &
                          ZTT, ZRV, ZRC, ZRI, ZU, ZV, ZW,                                        &
                          IEDUMMY, ZTTENE(:,:,1), ZRVTENE(:,:,1), ZRCTENE(:,:,1), ZRITENE(:,:,1),&
                          ZPRLTENE(:,1), ZPRSTENE(:,1),                                          &
                          IEDUMMY, IEDUMMY, ZPRLFLXE(:,:,1), ZPRSFLXE(:,:,1),                    &
                          ZUMFE(:,:,1), ZDMFE(:,:,1), ZEDUMMY,                                   &
                          OCHTRANS, ICH1, ZCH1, ZCH1TENE(:,:,:,1),                               &
                          OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX, ODUST, OSALT,                 &
                          ZRHODREF, ZIC_RATEE(:,1), ZCG_RATEE(:,1)                               )
    END IF
!
    IF ( IENS > 1 ) THEN
!
      CALL INI_CONVPAR
!
!* second member (positive vertical velocity perturb for Trigger)
!
      CALL DEEP_CONVECTION( ILON, IKU, IIDIA, IFDIA, IBDIA, ITDIA,                                &
                           PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,                          &
                           ZPABS, ZZ, ZDXDY, ZTIMEC,                                              &
                           ZTT, ZRV, ZRC, ZRI, ZU, ZV, ZW*1.5+1.E-4,                              &
                           IEDUMMY, ZTTENE(:,:,2), ZRVTENE(:,:,2), ZRCTENE(:,:,2), ZRITENE(:,:,2),&
                           ZPRLTENE(:,2), ZPRSTENE(:,2),                                          &
                           IEDUMMY, IEDUMMY, ZPRLFLXE(:,:,2), ZPRSFLXE(:,:,2),                    &
                           ZUMFE(:,:,2), ZDMFE(:,:,2), ZEDUMMY,                                   &
                           OCHTRANS, ICH1, ZCH1, ZCH1TENE(:,:,:,2),                               &
                           OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX, ODUST, OSALT,                 &
                           ZRHODREF, ZIC_RATEE(:,2), ZCG_RATEE(:,2)                               )
    END IF
!
    IF ( IENS > 2 ) THEN
!
!* third member (positive vertical velocity perturb for Trigger)
!
      CALL DEEP_CONVECTION( ILON, IKU, IIDIA, IFDIA, IBDIA, ITDIA,                                &
                           PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,                          &
                           ZPABS, ZZ, ZDXDY, ZTIMEC,                                              &
                           ZTT, ZRV, ZRC, ZRI, ZU, ZV, ZW*.5-1.E-4,                               &
                           IEDUMMY, ZTTENE(:,:,3), ZRVTENE(:,:,3), ZRCTENE(:,:,3), ZRITENE(:,:,3),&
                           ZPRLTENE(:,3), ZPRSTENE(:,3),                                          &
                           IEDUMMY, IEDUMMY, ZPRLFLXE(:,:,3), ZPRSFLXE(:,:,3),                    &
                           ZUMFE(:,:,3), ZDMFE(:,:,3), ZEDUMMY,                                   &
                           OCHTRANS, ICH1, ZCH1, ZCH1TENE(:,:,:,3),                               &
                           OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX, ODUST, OSALT,                 &
                           ZRHODREF, ZIC_RATEE(:,3), ZCG_RATEE(:,3)                               )
    END IF
!
ELSE                        
  ICOUNT(:)=0
  ZTTEN(:,:)=0.
  ZRVTEN(:,:)=0.
  ZRCTEN(:,:)=0.
  ZRITEN(:,:)=0.
  ZUMF(:,:)=0
  ZDMF(:,:)=0
  ICLTOP(:)=0
  ICLBAS(:)=0
  ZCH1TEN(:,:,:)=0.
  ZPRLTEN(:)=0.
  ZPRSTEN(:)=0.
  ZPRLFLX(:,:)=0.
  ZPRSFLX(:,:)=0.
  ZCAPE(:)=0.
  ZIC_RATE(:) = 0.0
  ZCG_RATE(:) = 0.0
END IF
!
!*       4.b  Call shallow convection routine
!             -------------------------------
!
IF ( HSCONV == 'KAFR' ) THEN
!
  IF ( HDCONV == 'NONE' ) CALL INI_CONVPAR          
  CALL INI_CONVPAR_SHAL
!
  CALL SHALLOW_CONVECTION( ILON, IKU, IIDIA, IFDIA, IBDIA, ITDIA,      &
                           PDTCONV, KICE, OSETTADJ, PTADJS,            &
                           ZPABS, ZZ,ZTKECLS,                          &
                           ZTT, ZRV, ZRC, ZRI, ZW,                     &
                           ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,          &
                           ICLTOPS, ICLBASS, ZUMFS,                    &
                           OCHTRANS, ICH1, ZCH1, ZCH1TENS              )
ELSE                        
  ZTTENS(:,:)=0.
  ZRVTENS(:,:)=0.
  ZRCTENS(:,:)=0.
  ZRITENS(:,:)=0.
  ZUMFS(:,:)=0
  ICLTOPS(:)=0
  ICLBASS(:)=0
  ZCH1TENS(:,:,:)=0.
END IF
!
!----------------------------------------------------------------------------
!
!*       5.  Reskip 2D convective tendencies in 3D arrays
!            Add deep and shallow tendencies, and - if activated - ensemble memebers
!            -----------------------------------------------------------------------
!
! Ensemble members for deep
ZSUM = 1.
IF ( IENS > 0 ) THEN
    IF ( IENS == 1 ) ZWEIGHT(:) = .5
    IF ( IENS >  1 ) ZWEIGHT(:) = 1.
    DO JN = 1, IENS
       ZTTEN(:,:)  = ZTTEN(:,:)  + ZWEIGHT(JN) * ZTTENE(:,:,JN)
       ZRVTEN(:,:) = ZRVTEN(:,:) + ZWEIGHT(JN) * ZRVTENE(:,:,JN)
       ZRCTEN(:,:) = ZRCTEN(:,:) + ZWEIGHT(JN) * ZRCTENE(:,:,JN)
       ZRITEN(:,:) = ZRITEN(:,:) + ZWEIGHT(JN) * ZRITENE(:,:,JN)
       ZPRLFLX(:,:)= ZPRLFLX(:,:)+ ZWEIGHT(JN) * ZPRLFLXE(:,:,JN)
       ZPRSFLX(:,:)= ZPRSFLX(:,:)+ ZWEIGHT(JN) * ZPRSFLXE(:,:,JN)
       ZUMF(:,:)   = ZUMF(:,:)   + ZWEIGHT(JN) * ZUMFE(:,:,JN)
       ZDMF(:,:)   = ZDMF(:,:)   + ZWEIGHT(JN) * ZDMFE(:,:,JN)
       ZPRLTEN(:)  = ZPRLTEN(:)  + ZWEIGHT(JN) * ZPRLTENE(:,JN)
       ZPRSTEN(:)  = ZPRSTEN(:)  + ZWEIGHT(JN) * ZPRSTENE(:,JN)
       IF ( OCHTRANS )  &
         & ZCH1TEN(:,:,:) = ZCH1TEN(:,:,:) + ZWEIGHT(JN) * ZCH1TENE(:,:,:,JN)
       IF ( OCH_CONV_LINOX ) THEN
         ZIC_RATE(:)  = ZIC_RATE(:)  + ZWEIGHT(JN) * ZIC_RATEE(:,JN)
         ZCG_RATE(:)  = ZCG_RATE(:)  + ZWEIGHT(JN) * ZCG_RATEE(:,JN)
       END IF
    END DO
!
    ZSUM = 1. / ( 1. + SUM( ZWEIGHT(:) ) )
END IF
!
! Add deep and shallow and reskip
DO JK = 1,IKU
  IIND = 1
  DO JI = IIB,IIE
    DO JJ = IJB,IJE
      PTHTEN(JI,JJ,JK) = ( ZTTEN(IIND,JK) * ZSUM + ZTTENS(IIND,JK) ) * ZPI(IIND,JK)
      PRVTEN(JI,JJ,JK) = ZRVTEN(IIND,JK)  * ZSUM + ZRVTENS(IIND,JK)
      PRCTEN(JI,JJ,JK) = ZRCTEN(IIND,JK)  * ZSUM + ZRCTENS(IIND,JK)
      PRITEN(JI,JJ,JK) = ZRITEN(IIND,JK)  * ZSUM + ZRITENS(IIND,JK)
      IIND = IIND + 1
    END DO
  END DO
END DO
IIND = 1
DO JI = IIB,IIE
  DO JJ = IJB,IJE
    PPRTEN(JI,JJ)  = ( ZPRLTEN(IIND) + ZPRSTEN(IIND) ) * ZSUM 
    PPRSTEN(JI,JJ) = ZPRSTEN(IIND) * ZSUM
    KCOUNT(JI,JJ)  = ICOUNT(IIND)
    IIND = IIND + 1
  END DO
END DO
IF ( OCHTRANS ) THEN
  DO JK = 1,IKU
    IIND = 1
    DO JI = IIB,IIE
      DO JJ = IJB,IJE
        PCH1TEN(JI,JJ,JK,:) = ZCH1TEN(IIND,JK,:) * ZSUM + ZCH1TENS(IIND,JK,:)
        IIND = IIND + 1
      END DO
    END DO
  END DO
END IF
!
IF ( OCH_CONV_LINOX ) THEN
  IIND = 1
  DO JI = IIB,IIE
    DO JJ = IJB,IJE
      PIC_RATE(JI,JJ) = ZIC_RATE(IIND) * ZSUM 
      PCG_RATE(JI,JJ) = ZCG_RATE(IIND) * ZSUM 
      IIND = IIND + 1
    END DO
  END DO
END IF
!
!convective mass flux for subgrid condensation
IF(SIZE(PMF) /= 0) THEN
  DO JK = 1,IKU
    IIND = 1
    DO JI = IIB,IIE
      DO JJ = IJB,IJE
        PMF(JI,JJ,JK)   = ( ZUMF(IIND,JK) + ZDMF(IIND,JK) ) * ZSUM + ZUMFS(IIND,JK) 
        IIND = IIND + 1
      END DO
    END DO
  END DO
ENDIF
!
! Diagnostics:
IF ( ODIAGCONV ) THEN
  DO JK = 1,IKU
    IIND = 1
    DO JI = IIB,IIE
      DO JJ = IJB,IJE
        PPRLFLX(JI,JJ,JK)= ZPRLFLX(IIND,JK) * ZSUM
        PPRSFLX(JI,JJ,JK)= ZPRSFLX(IIND,JK) * ZSUM
        PUMF(JI,JJ,JK)   = ZUMF(IIND,JK) * ZSUM + ZUMFS(IIND,JK)
        PDMF(JI,JJ,JK)   = ZDMF(IIND,JK) * ZSUM
        IIND = IIND + 1
      END DO
    END DO
  END DO
  IIND = 1
  DO JI = IIB,IIE
    DO JJ = IJB,IJE
      PCAPE(JI,JJ)   = ZCAPE(IIND)
      KCLTOP(JI,JJ)  = MAX(ICLTOP(IIND), ICLTOPS(IIND))
      KCLBAS(JI,JJ)  = MAX(ICLBAS(IIND), ICLBASS(IIND))
      IIND = IIND + 1
    END DO
  END DO
END IF
!
DEALLOCATE( ICLBASS )
DEALLOCATE( ICLTOPS )
DEALLOCATE( ZUMFS )
DEALLOCATE( ZCH1TENS ) 
DEALLOCATE( ZRCTENS )
DEALLOCATE( ZRITENS ) 
DEALLOCATE( ZTTENS )
DEALLOCATE( ZRVTENS ) 
!
DEALLOCATE( ZTIMEC )
DEALLOCATE( ZCAPE )
DEALLOCATE( ICOUNT )
DEALLOCATE( ICLBAS )
DEALLOCATE( ICLTOP ) 
DEALLOCATE( ZUMF )
DEALLOCATE( ZDMF ) 
DEALLOCATE( ZPRLFLX)
DEALLOCATE( ZPRSFLX ) 
DEALLOCATE( ZCH1 )
DEALLOCATE( ZRHODREF )
DEALLOCATE( ZCH1TEN ) 
DEALLOCATE( ZRCTEN )
DEALLOCATE( ZRITEN ) 
DEALLOCATE( ZPRLTEN )
DEALLOCATE( ZPRSTEN ) 
DEALLOCATE( ZTTEN )
DEALLOCATE( ZRVTEN ) 
DEALLOCATE( ZW )
DEALLOCATE( ZDXDY )
DEALLOCATE( ZU )
DEALLOCATE( ZV ) 
DEALLOCATE( ZRC )
DEALLOCATE( ZRI ) 
DEALLOCATE( ZPI )
DEALLOCATE( ZTT )
DEALLOCATE( ZRV ) 
DEALLOCATE( ZPABS )
DEALLOCATE( ZZ )
DEALLOCATE( ZIC_RATE )
DEALLOCATE( ZCG_RATE )
!
!* Deallocation of additional ensemble members
IF ( IENS > 0 ) THEN
     DEALLOCATE( ZTTENE )
     DEALLOCATE( ZRVTENE )
     DEALLOCATE( ZRCTENE )
     DEALLOCATE( ZRITENE )
     DEALLOCATE( ZUMFE )
     DEALLOCATE( ZDMFE )
     DEALLOCATE( ZCH1TENE )
     DEALLOCATE( ZPRLFLXE )
     DEALLOCATE( ZPRSFLXE )
     DEALLOCATE( ZPRLTENE )
     DEALLOCATE( ZPRSTENE )
     DEALLOCATE( ZEDUMMY )
     DEALLOCATE( IEDUMMY )
     DEALLOCATE( ZWEIGHT )
     DEALLOCATE( ZIC_RATEE )
     DEALLOCATE( ZCG_RATEE )
END IF
!
END SUBROUTINE CONVECTION
!
!----------------------------------------------------------------------------


