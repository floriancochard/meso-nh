! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_FM
USE MODE_ll
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
USE MODD_FIELD_n, ONLY: XUM,XVM,XWM
USE MODD_DYN_n
USE MODD_REF_n
USE MODD_CST
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
INTEGER      :: IIB, IIE ! inner limits in X direction
INTEGER      :: IJB, IJE ! inner limits in Y direction
INTEGER      :: IKB, IKE ! inner limits in Z direction
INTEGER      :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
!JUAN
REAL                  :: ZMAXVAL,ZRESIDUAL
INTEGER, DIMENSION(3) :: IMAXLOC
INTEGER               :: I,J,K
!JUAN
!-------------------------------------------------------------------------------
!
!*       1.    Initialisations
!              ---------------
!
CALL FMLOOK_ll(CLUOUT0,CLUOUT0,ILUOUT0,IRESP)
!
IIB=1+JPHEXT
IIE=NIMAX+JPHEXT
IJB=1+JPHEXT
IJE=NJMAX+JPHEXT
IKB=1+JPVEXT
IKE=NKMAX+JPVEXT
!
ZU(:,:,:) = XUM(:,:,:)
ZV(:,:,:) = XVM(:,:,:)
ZW(:,:,:) = XWM(:,:,:)
!
NULLIFY(TZFIELDS_ll)
!
!-------------------------------------------------------------------------------
!
!*       2.    Loop
!              ----
!
DO
  XUM(:,:,:) = ZU(:,:,:)
  XVM(:,:,:) = ZV(:,:,:)
  XWM(:,:,:) = ZW(:,:,:)
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
  CALL ANEL_BALANCE_n('M',ZRESIDUAL)
!
!-------------------------------------------------------------------------------
!
!*       4.    compute the residual divergence
!              -------------------------------
!
  ZRU(:,:,:) = XUM(:,:,:) * MXM(XRHODJ)
  ZRV(:,:,:) = XVM(:,:,:) * MYM(XRHODJ)
  ZRW(:,:,:) = XWM(:,:,:) * MZM(XRHODJ)
!
  CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRU)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRV)
  CALL ADD3DFIELD_ll(TZFIELDS_ll, ZRW)
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
  CALL GDIV(CLBCX,CLBCY,PDXX,PDYY,PDZX,PDZY,PDZZ,ZRU,ZRV,ZRW,ZDIV)
!
  IF ( CEQNSYS=='DUR' ) THEN
    IF ( SIZE(XRVREF,1) == 0 ) THEN
      ZDIV=ZDIV/XRHODJ/XTH00*XRHODREF*XTHVREF
    ELSE
      ZDIV=ZDIV/XRHODJ/XTH00*XRHODREF*XTHVREF*(1.+XRVREF)
    END IF
    ZMAXVAL=MAX_ll( ABS(ZDIV),IINFO_ll)
    IMAXLOC=MAXLOC( ABS(ZDIV(IIB:IIE,IJB:IJE,IKB:IKE)))
    WRITE(ILUOUT0,*) 'JUAN1 residual divergence / 2 DT = ',   &
    ZMAXVAL,  ' located at ',    &
    IMAXLOC
    WRITE(ILUOUT0,*) 'JUAN1 residual divergence / 2 DT = ',   &
    ZRESIDUAL
    DO K=1,size(ZDIV,3)
    DO J=1,size(ZDIV,2)
    DO I=1,size(ZDIV,1)
     IF ( ABS(ZDIV(I,J,K)) .EQ. ZMAXVAL ) THEN
	PRINT*,"I=",I," J=",J," K=",K," ZMAXVAL=",ZDIV(I,J,K)
	PRINT*,"SI=",size(ZDIV,1)," SJ=",size(ZDIV,2)," SK=",size(ZDIV,3)
     ENDIF
    ENDDO
    ENDDO
    ENDDO
  ELSEIF( CEQNSYS=='MAE' .OR. CEQNSYS=='LHE' ) THEN
    ZDIV=ZDIV/XRHODJ*XRHODREF
    ZMAXVAL=MAX_ll( ABS(ZDIV),IINFO_ll)
    IMAXLOC=MAXLOC( ABS(ZDIV(IIB:IIE,IJB:IJE,IKB:IKE)))
    WRITE(ILUOUT0,*) 'JUAN2 residual divergence / 2 DT = ',   &
    ZMAXVAL,  ' located at ',    &
    IMAXLOC
  END IF
!
!-------------------------------------------------------------------------------
!
!*       5.    modifies the solver parameters if necessary
!              -------------------------------------------
!
  ZMAXVAL = ZRESIDUAL
  IF ( ZMAXVAL > 1.E-8 ) THEN
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
      STOP
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
