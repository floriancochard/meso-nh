!MNH_LIC Copyright 1998-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
MODULE MODI_PRESSURE_IN_PREP
!     ######################
!
INTERFACE
!
      SUBROUTINE PRESSURE_IN_PREP(PDXX,PDYY,PDZX,PDZY,PDZZ)
!
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDXX     ! metric coefficient dxx
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDYY     ! metric coefficient dyy 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDZX     ! metric coefficient dzx 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDZY     ! metric coefficient dzy 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz  
!
END SUBROUTINE PRESSURE_IN_PREP
!
END INTERFACE
!
END MODULE MODI_PRESSURE_IN_PREP
!
!     #####################################################
      SUBROUTINE PRESSURE_IN_PREP(PDXX,PDYY,PDZX,PDZY,PDZZ)
!     #####################################################
!
!!****  *PRESSURE_IN_PREP* - calls the pressure solver in prep_real_case and
!!                           checks the result
!!
!!    PURPOSE
!!    -------
!!
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!      
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    22/12/98
!!      parallelization                                   18/06/00 (Jabouille)
!!                  2014 M.Faivre
!!               08/2015 M.Moge    removing UPDATE_HALO_ll on XUT, XVT, XRHODJ in part 4
!!   J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODE_MSG
!
USE MODI_ANEL_BALANCE_n
USE MODI_GDIV
USE MODI_SHUMAN
!
USE MODD_CONF            ! declaration modules
USE MODD_CONF_n
USE MODD_LUNIT
USE MODD_DIM_n
USE MODD_GRID_n
USE MODD_LBC_n
USE MODD_PARAMETERS
USE MODD_FIELD_n, ONLY: XUT,XVT,XWT
USE MODD_DYN_n
USE MODD_REF_n
USE MODD_CST
USE MODE_MPPDB
USE MODE_EXTRAPOL
!
IMPLICIT NONE
!
!*       0.1   Declaration of dummy arguments
!              ------------------------------
!
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDXX     ! metric coefficient dxx
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDYY     ! metric coefficient dyy 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDZX     ! metric coefficient dzx 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDZY     ! metric coefficient dzy 
REAL,DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz  
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
REAL,DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: ZU   ! U
REAL,DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: ZV   ! V
REAL,DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: ZW   ! W
REAL,DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: ZRU  ! U * rho * J
REAL,DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: ZRV  ! V * rho * J
REAL,DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: ZRW  ! W * rho * J
REAL,DIMENSION(SIZE(PDXX,1),SIZE(PDXX,2),SIZE(PDXX,3)):: ZDIV ! residual divergence
!
!* file management variables and counters
!
INTEGER      :: ILUOUT0  ! logical unit for listing file
INTEGER      :: IRESP    ! error code
INTEGER      :: IKB, IKE ! inner limits in Z direction
INTEGER      :: IKU
INTEGER      :: IINFO_ll
REAL         :: ZMAXRES
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
REAL                  :: ZMAXVAL,ZRESIDUAL
INTEGER, DIMENSION(3) :: IMAXLOC
INTEGER               :: I,J,K
!-------------------------------------------------------------------------------
!
!*       1.    Initialisations
!              ---------------
!
ILUOUT0 = TLUOUT0%NLU
!
IKB=1+JPVEXT
IKE=NKMAX+JPVEXT
IKU=IKE+JPVEXT
!
ZU(:,:,:) = XUT(:,:,:)
ZV(:,:,:) = XVT(:,:,:)
ZW(:,:,:) = XWT(:,:,:)
!
NULLIFY(TZFIELDS_ll)
!
!-------------------------------------------------------------------------------
!
!*       2.    Loop
!              ----
!
DO
  XUT(:,:,:) = ZU(:,:,:)
  XVT(:,:,:) = ZV(:,:,:)
  XWT(:,:,:) = ZW(:,:,:)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XUT)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XVT)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XWT)
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)

  CALL MPPDB_CHECK3D(XUT,"PREP::XUM",PRECISION)
  CALL MPPDB_CHECK3D(XVT,"PREP::XVM",PRECISION)
  CALL MPPDB_CHECK3D(XWT,"PREP::XWM",PRECISION)
  CALL MPPDB_CHECK3D(XRHODJ,"PREP::XRHODJ",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       3.    ANELASTIC CORRECTION
!              --------------------
!
  WRITE(ILUOUT0,*) ' '
  WRITE(ILUOUT0,*) 'CPRESOPT = ',CPRESOPT
  WRITE(ILUOUT0,*) 'NITR     = ',NITR
  IF (CPRESOPT=='RICHA') &
    WRITE(ILUOUT0,*) 'XRELAX   = ',XRELAX
!
  CALL ANEL_BALANCE_n(ZRESIDUAL)
!
!-------------------------------------------------------------------------------
!
!*       4.    compute the residual divergence
!              -------------------------------
!20140225 forgot this update_halo
!20131112 check 1st XUT
CALL MPPDB_CHECK3D(XUT,"PressInP4-beforeupdhalo::XUT",PRECISION)
CALL MPPDB_CHECK3D(XVT,"PressInP4-beforeupdhalo::XVT",PRECISION)
!CALL ADD3DFIELD_ll(TZFIELDS_ll, XUT)
!CALL ADD3DFIELD_ll(TZFIELDS_ll, XVT)
!CALL ADD3DFIELD_ll(TZFIELDS_ll, XRHODJ)
!  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
!    CALL CLEANLIST_ll(TZFIELDS_ll)
!
  ZRU(:,:,:) = XUT(:,:,:) * MXM(XRHODJ)
  ZRV(:,:,:) = XVT(:,:,:) * MYM(XRHODJ)
  ZRW(:,:,:) = XWT(:,:,:) * MZM(1,IKU,1,XRHODJ)
!
  CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRU)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRV)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRW)
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
  CALL GDIV(CLBCX,CLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,ZRU,ZRV,ZRW,ZDIV)
CALL MPPDB_CHECK3D(XUT,"PressInP4-afterupdhalo::XUT",PRECISION)
CALL MPPDB_CHECK3D(XVT,"PressInP4-afterupdhalo::XVT",PRECISION)
!
!20131125 add extrapol on ZRU
CALL EXTRAPOL('W',ZRU)
CALL MPPDB_CHECK3D(ZRU,"PressInP4-afterextrapol W::ZRU",PRECISION)
!
!20131126 add extrapol on ZRV
CALL EXTRAPOL('S',ZRV)
!
  IF ( CEQNSYS=='DUR' ) THEN
    IF ( SIZE(XRVREF,1) == 0 ) THEN
      ZDIV=ZDIV/XRHODJ/XTH00*XRHODREF*XTHVREF
    ELSE
      ZDIV=ZDIV/XRHODJ/XTH00*XRHODREF*XTHVREF*(1.+XRVREF)
    END IF
  ELSEIF( CEQNSYS=='MAE' .OR. CEQNSYS=='LHE' ) THEN
    ZDIV=ZDIV/XRHODJ*XRHODREF
  END IF
!
!-------------------------------------------------------------------------------
!
!*       5.    modifies the solver parameters if necessary
!              -------------------------------------------
!
  IF (LRES) THEN
     ZMAXRES = XRES
            ELSE
     ZMAXRES = XRES_PREP
  END IF
!
  ZMAXVAL = ZRESIDUAL
  IF (  ZMAXVAL > ZMAXRES) THEN
!!$  IF (LRES) THEN
!!$     ZMAXRES = XRES
!!$            ELSE
!!$     ZMAXRES = 1.E-08
!!$  END IF
!!$!
!!$  IF ( MAX_ll( ABS(ZDIV),IINFO_ll) > ZMAXRES ) THEN
    IF (CPRESOPT=='RICHA') THEN
      IF (XRELAX>0.75) THEN
        XRELAX  = XRELAX - 0.1
      ELSE
        XRELAX  = 1.
        NITR = NITR + 4
      END IF
    ELSE
      NITR = NITR + 4
    END IF
!
    IF (NITR>40) THEN
      WRITE(ILUOUT0,*) ' '
      WRITE(ILUOUT0,*) '******************************************************************************'
      WRITE(ILUOUT0,*) '*                                                                            *'
      IF (CPROGRAM=='REAL  ') WRITE(ILUOUT0,*) &
                       '*                        WARNING in PREP_REAL_CASE                           *'
      IF (CPROGRAM=='IDEAL  ') WRITE(ILUOUT0,*) &
                       '*                        WARNING in PREP_IDEAL_CASE                          *'
      WRITE(ILUOUT0,*) '*                                                                            *'
      WRITE(ILUOUT0,*) '******************************************************************************'
      WRITE(ILUOUT0,*) ' '
      WRITE(ILUOUT0,*) 'The pressure solver was unable to converge for your case.'
      WRITE(ILUOUT0,*) 'This can be due to : '
      WRITE(ILUOUT0,*) '                a locally very steep orography'
      WRITE(ILUOUT0,*) '                a too high vertical stretching'
      WRITE(ILUOUT0,*) '                a too high horizontal stretching on the sphere (large domains)'
      WRITE(ILUOUT0,*) '                a too high horizontal stretching in conformal plane'
      WRITE(ILUOUT0,*) '                an error in your input velocity and thermodynamical fields'
      WRITE(ILUOUT0,*) ' '
      WRITE(ILUOUT0,*) '******************************************************************************'
      WRITE(ILUOUT0,*) '******************************************************************************'
      WRITE(ILUOUT0,*) ' '
      WRITE(ILUOUT0,*) 'The model will probably not be able to run with this initial or coupling file.'
      WRITE(ILUOUT0,*) ' '
      WRITE(ILUOUT0,*) '******************************************************************************'
      WRITE(ILUOUT0,*) '******************************************************************************'
      WRITE(ILUOUT0,*) ' '
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','PRESSURE_IN_PREP','')
    END IF
  ELSE
!*       7.    Happy conclusion
!              ----------------
!
    WRITE(ILUOUT0,*) ' '
    IF (.NOT. LCARTESIAN) THEN
      WRITE(ILUOUT0,*) 'Horizontal stretching in conformal plane: '
      WRITE(ILUOUT0,*) ' map factor varies between ', MINVAL(XMAP),' and ', MAXVAL(XMAP)
      WRITE(ILUOUT0,*) ' '
    ENDIF
    WRITE(ILUOUT0,*) '***************************************************'
    WRITE(ILUOUT0,*) 'The variables CPRESOPT = ',CPRESOPT
    WRITE(ILUOUT0,*) '              NITR     = ',NITR
    IF (CPRESOPT=='RICHA') &
      WRITE(ILUOUT0,'(A27,F3.1)') '               XRELAX   =  ',XRELAX
    WRITE(ILUOUT0,*) ' '
    WRITE(ILUOUT0,*) 'LOOK correct for the pressure problem in your case.'
    WRITE(ILUOUT0,*) 'They are stored in the coupling file, and will be'
    WRITE(ILUOUT0,*) 'the default for the model run.'
    WRITE(ILUOUT0,*) '***************************************************'
    WRITE(ILUOUT0,*) ' '
    EXIT
  END IF
!
!
!*       6.    Next try
!              --------
!
END DO
!
!
END SUBROUTINE PRESSURE_IN_PREP
