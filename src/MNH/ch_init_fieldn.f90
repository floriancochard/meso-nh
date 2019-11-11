!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/07/06 15:28:51
!-----------------------------------------------------------------
!!    ########################### 
      MODULE MODI_CH_INIT_FIELD_n
!!    ########################### 
!!
INTERFACE
!!
SUBROUTINE CH_INIT_FIELD_n(KMI, KLUOUT, KVERB)
!!
IMPLICIT NONE
!!
INTEGER, INTENT(IN)  :: KMI      ! model index
INTEGER, INTENT(IN)  :: KLUOUT   ! output listing channel
INTEGER, INTENT(IN)  :: KVERB    ! verbosity level
!!
!!
END SUBROUTINE CH_INIT_FIELD_n
!!
END INTERFACE
!!
END MODULE MODI_CH_INIT_FIELD_n
!!
!!    ##############################################
      SUBROUTINE CH_INIT_FIELD_n(KMI, KLUOUT, KVERB)
!!    ##############################################
!!
!!*** *CH_INIT_FIELD_n*
!!
!!    PURPOSE
!!    -------
!        initialize MesoNH scalar variables
!!
!!**  METHOD
!!    ------
!!       The subroutine CH_FIELD_VALUE_n returns for each grid-point 
!!    (LAT,LON,ZZ) and each species a corresponding initial value, either
!!    in part/part or in molec/cm3. If necessary, that initial value is 
!!    then converted to mixing ratio (part/part).
!!    The variables at time t and t-dt are given identic values.
!!    Presently, there is only a 1D initialization (homogeneous in x-y)
!!    available. For more sophisticated initializations, the subroutine
!!    CH_FIELD_VALUE_n may be modified by the user. The character parameter
!!    CCH_INIT_FIELD_OPT may be used in order to pass user specific information
!!    on to that subroutine. These subroutines have been duplicated in order
!!    to allow future inclusion of model dependant parameters (like an
!!    initialization that depends on variables stored in MODD_FIELD_n)
!!
!!    REFERENCE
!!    ---------
!!    book 2 of MesoNH
!!
!!    AUTHOR
!!    ------
!!    K. Suhre   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/11/95
!!    05/08/96 (K. Suhre) restructured
!!    11/08/98 (N. Asencio) add parallel code
!!    09/03/99 (V. Crassier) speed up 1-D initialization by reducing a 3-D
!!                           loop to a 1-D loop with 10m precision
!!    09/12/99 (K. Suhre) add missing update halo and a fix for MAXVAL pbs.
!!    09/01/01 (P. Tulet) initialize chemical constant (molar mass, henry 
!!                        specific constant and biological reactivity
!!    22/01/01 (D. Gazen) add NSV_CHEMBEG and NSV_CHEMEND indices to handle SV
!!    04/06/07 (M. Leriche & JP Pinty) add pH initialization
!!    20/04/10 (M. Leriche) remove pH initialization to ini_modeln
!!
!!    EXTERNAL
!!    --------
!!     GET_DIM_EXT_ll : get extended sub-domain sizes
!!     GET_INDICE_ll  : get physical sub-domain bounds
!!
!!------------------------------------------------------------------------------
!!
USE MODI_CH_FIELD_VALUE_n ! returns value of chemical species at each grid point
USE MODI_CH_INIT_CONST_n
USE MODI_CH_AER_EQM_INIT_n
USE MODE_ll
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_GRID_n,      ONLY : XZZ,       &! height z
                             XLAT,XLON   ! latitude and longitude
USE MODD_REF_n,       ONLY : XRHODREF,  &! dry density of ref. state
                             XRHODJ      ! ( rhod J ) = dry density
USE MODD_LBC_n
USE MODD_NSV,         ONLY : NSV_CHEM, NSV_CHEMBEG,NSV_CHEMEND, &
                             NSV_AER, NSV_AERBEG,NSV_AEREND
USE MODD_CST,         ONLY : XMD, XAVOGADRO

USE MODD_FIELD_n,     ONLY : XSVT       ! scalar variable at t
USE MODD_PARAMETERS,  ONLY : JPVEXT, JPHEXT  ! number of External points
USE MODD_ARGSLIST_ll, ONLY : LIST_ll  ! for update_halo
USE MODD_CH_CONST_n                   ! for Chemical constants
USE MODD_CONF,        ONLY : CPROGRAM, L1D, L2D 
USE MODD_CONF_n,      ONLY : NRRL
USE MODD_CH_MNHC_n
USE MODD_CH_M9_n,     ONLY : CNAMES, NEQ
USE MODD_CH_AEROSOL
USE MODD_CH_AERO_n
USE MODD_LSFIELD_n,   ONLY : XLBXSVM, XLBYSVM
USE MODD_DYN_n,       ONLY : NRIMX,NRIMY
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)  :: KMI      ! model index
INTEGER, INTENT(IN)  :: KLUOUT   ! output listing channel
INTEGER, INTENT(IN)  :: KVERB    ! verbosity level
!
!*      0.2    declarations local variables
!
INTEGER :: JI, JJ, JK, JN        ! loop control variables
CHARACTER(LEN=3) :: YUNIT              ! units of returned initial values
                                       ! "CON" = molec./cm3
                                       ! "MIX" = mixing ratio
REAL    :: ZDEN2MOL
        !  ZDEN2MOL = 6.0221367E+23 * 1E-6 / 28.9644E-3
        !  conversion factor density to mol/cm3
        !  n_molec (moelc./cm3):  M = 1E-6*RHO(kg/m3) * XAVOGADRO / XMD

REAL, ALLOCATABLE, DIMENSION(:)   :: ZHEIGHT   !Height lookup table
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZSVINIT   !Species concentration lookup table  
REAL, ALLOCATABLE, DIMENSION(:,:) :: ZSVINITA  !Aerosols species concentration lookup table  

INTEGER             :: ILEVMAX   !Maximum height level
INTEGER             :: JLEV      !Current height level

INTEGER             :: IIU  ! Upper dimension in x direction
INTEGER             :: IJU  ! Upper dimension in y direction
INTEGER             :: IKU  ! Upper dimension in z direction
INTEGER             :: IIB  ! indice I Beginning in x direction
INTEGER             :: IJB  ! indice J Beginning in y direction
INTEGER             :: IKB  ! indice K Beginning in z direction
INTEGER             :: IIE  ! indice I End       in x direction
INTEGER             :: IJE  ! indice J End       in y direction
INTEGER             :: IKE  ! indice K End       in z direction
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll ! pointer for the list of 3D fields
INTEGER :: IINFO_ll  ! Return code of //routines
INTEGER :: IOR, JOR, IEND, JEND, KINFO, NIU,NJU, ILBX, ILBY, IRIMX, IRIMY
!
!-------------------------------------------------------------------------------
!
!*       0.    PROLOGUE
!              --------
!
NULLIFY(TZFIELDS_ll)
!
!*       1.    PREPARE INITIALIZATION
!              ----------------------

IF (CORGANIC == TRIM("MPMPO") .OR. CORGANIC == TRIM("PUN") .OR. CORGANIC == TRIM("EQSAM2")) THEN
  IF ((CCH_SCHEME .EQ. TRIM("NONE")) .OR. (CCH_SCHEME .EQ. TRIM("RELACS"))&
   .OR. (CCH_SCHEME .EQ. TRIM("RACM"))) THEN
  WRITE(KLUOUT,FMT=*) '**********************************************'
  WRITE(KLUOUT,FMT=*) 'WARNING : NO SOA !!!!'
  WRITE(KLUOUT,FMT=*) 'YOU WANT TO USE SOA GAS PARTICLE BALANCE'
  WRITE(KLUOUT,FMT=*) 'BUT THE SCHEME NEED TO BE CACM or RELACS 2'
  WRITE(KLUOUT,FMT=*) 'CORGANIC HAS BEEN SET TO NONE'
  WRITE(KLUOUT,FMT=*) 'OTHERWISE COMPILE THE CORRECT SCHEME BEFORE'
  WRITE(KLUOUT,FMT=*) '**********************************************'
  CORGANIC = "NONE"
  END IF
END IF
!
!*       1.1   compute dimensions of arrays
!

CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_GLOBALDIMS_ll(IIE,IJE)
IIU = IIE + 2 * JPHEXT
IJU = IJE + 2 * JPHEXT
CALL GET_INTERSECTION_ll(1+JPHEXT, 1+JPHEXT, IIU-JPHEXT , IJU-JPHEXT, IOR, JOR, IEND, JEND, "EXTE", KINFO)
IKB = 1 + JPVEXT
IKU = SIZE(XSVT,3)
IKE = IKU - JPVEXT
CALL GET_DIM_EXT_ll('B',NIU,NJU)
!
!        1.1.1 find maximum height level
ILEVMAX=INT(MAXVAL(XZZ(:,:,IKE)/10.))+1
! the following print serves to break compiler optimization with MAXVAL
! (pb. on OS2000 with option -O3 for example, Peter Bechtold had
! similar surprises with MAXVAL on Fuji VPP700 in the convection scheme)
WRITE(KLUOUT,*) "CH_INIT_FIELD_n: ILEVMAX =",ILEVMAX 
ALLOCATE(ZHEIGHT(ILEVMAX))
ALLOCATE(ZSVINIT(ILEVMAX,NEQ))
ALLOCATE(ZSVINITA(ILEVMAX,NSV_AER))
!
!*       1.2   compute conversion factor kg/m3 --> molec/cm3
!
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD
!
!
!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZE T FIELDS AND CONVERT CONC. TO MIXING RATIO
!              -----------------------
!

YUNIT="MIX"

IF (LORILAM) THEN
  IF (.NOT.(ASSOCIATED(XN3D)))      ALLOCATE(XN3D(SIZE(XSVT,1),SIZE(XSVT,2),IKU,JPMODE))
  IF (.NOT.(ASSOCIATED(XRG3D)))     ALLOCATE(XRG3D(SIZE(XSVT,1),SIZE(XSVT,2),IKU,JPMODE))
  IF (.NOT.(ASSOCIATED(XSIG3D)))    ALLOCATE(XSIG3D(SIZE(XSVT,1),SIZE(XSVT,2),IKU,JPMODE))
  IF (.NOT.(ASSOCIATED(XRHOP3D)))   ALLOCATE(XRHOP3D(SIZE(XSVT,1),SIZE(XSVT,2),IKU,JPMODE))
  IF (.NOT.(ASSOCIATED(XM3D)))      ALLOCATE(XM3D(SIZE(XSVT,1),SIZE(XSVT,2),IKU,JPMODE*3))
  IF (.NOT.(ASSOCIATED(XSEDA)))     ALLOCATE(XSEDA(SIZE(XSVT,1),SIZE(XSVT,2),IKU,JPMODE*3))
  IF (.NOT.(ASSOCIATED(XCTOTA3D))) &
    ALLOCATE(XCTOTA3D(SIZE(XSVT,1),SIZE(XSVT,2),IKU,NSP+NCARB+NSOA,JPMODE))
  IF (.NOT.(ASSOCIATED(XVDEPAERO))) ALLOCATE(XVDEPAERO(SIZE(XSVT,1),SIZE(XSVT,2),JPIN))
  IF (.NOT.(ALLOCATED(XFAC)))       ALLOCATE(XFAC(NSP+NSOA+NCARB))
  IF (.NOT.(ALLOCATED(XRHOI)))      ALLOCATE(XRHOI(NSP+NSOA+NCARB))
  IF (.NOT.(ASSOCIATED(XFRAC))) THEN
    ALLOCATE(XFRAC(SIZE(XSVT,1),SIZE(XSVT,2),IKU,NEQ))
    XFRAC(:,:,:,:) = 0.
  END IF
  IF (.NOT.(ASSOCIATED(XMI))) THEN
    ALLOCATE(XMI(SIZE(XSVT,1),SIZE(XSVT,2),IKU,NSP+NCARB+NSOA))
  END IF
END IF
!
!*          print info for user
IF ((LCH_INIT_FIELD).AND.(CPROGRAM/='DIAG  ')) THEN
!
  WRITE(KLUOUT,*) "CH_INIT_FIELD_n will now initialize XSVT fields"
!
!
  jlev_loop : DO JLEV=1,ILEVMAX
    ZHEIGHT=REAL(JLEV-1)*10.
    jn_loop : DO JN = 1, NEQ
      ZSVINIT(JLEV,JN) = &
           CH_FIELD_VALUE_n(ZHEIGHT(JLEV), "LLZ", &
           CNAMES(JN), YUNIT, KLUOUT, KVERB)
      ! "LLZ" identifies the type of x-y-z values passed on to
      ! CH_FIELD_VALUE_n ("LLZ"=lon-lat-Z)
      ! in future developpements, "IJK" may be used in order
      ! to pass the grid indices rather than coordinates
    END DO jn_loop
  END DO jlev_loop

  XSVT(:,:,:,NSV_CHEMBEG:NSV_CHEMEND) = 0.
  jk_loop : DO JK = IKB, IKE
    jj_loop : DO JJ = JOR, JEND
      ji_loop : DO JI = IOR, IEND

        JLEV=INT(MAX(XZZ(JI,JJ,JK),0.)/10.)+1
        XSVT(JI,JJ,JK,NSV_CHEMBEG:NSV_CHEMEND) = ZSVINIT(JLEV,:)

      END DO ji_loop
    END DO jj_loop
  END DO jk_loop
  DO JN = NSV_CHEMBEG,NSV_CHEMEND
    DO JK=1,JPVEXT
      XSVT(:,:,IKB-JPVEXT,JN) = XSVT(:,:,IKB,JN)
      XSVT(:,:,IKE+JPVEXT,JN) = XSVT(:,:,IKE,JN)
    END DO
  END DO
  !
  IF (YUNIT .EQ. "CON") THEN
    WRITE(KLUOUT,*) "CH_INIT_FIELD_n: converting initial values to mixing ratio"
    DO JN = NSV_CHEMBEG,NSV_CHEMEND
      XSVT(:,:,:,JN) = XSVT(:,:,:,JN)/(XRHODREF(:,:,:)*ZDEN2MOL)
    ENDDO
  ELSE
    WRITE(KLUOUT,*)"CH_INIT_FIELD_n: initial values are used as is (mixing ratio)"
  ENDIF
!
!
  IF (LORILAM) THEN
  jlev_loop2 : DO JLEV=1,ILEVMAX
    ZHEIGHT=REAL(JLEV-1)*10.
    jn_loop2 : DO JN = 1, NSV_AER
      ZSVINITA(JLEV,JN) = &
           CH_FIELD_VALUE_n(ZHEIGHT(JLEV), "LLZ", &
           CAERONAMES(JN), YUNIT, KLUOUT, KVERB)
      ! "LLZ" identifies the type of x-y-z values passed on to
      ! CH_FIELD_VALUE_n ("LLZ"=lon-lat-Z)
      ! in future developpements, "IJK" may be used in order
      ! to pass the grid indices rather than coordinates
    END DO jn_loop2
  END DO jlev_loop2
  !
  XSVT(:,:,:,NSV_AERBEG:NSV_AEREND) = 0.
  jk_loop2 : DO JK = IKB, IKE
    jj_loop2 : DO JJ = JOR, JEND
      ji_loop2 : DO JI = IOR, IEND

        JLEV=INT(MAX(XZZ(JI,JJ,JK),0.)/10.)+1
        XSVT(JI,JJ,JK,NSV_AERBEG:NSV_AEREND) = ZSVINITA(JLEV,:)

      END DO ji_loop2
    END DO jj_loop2
  END DO jk_loop2
  DO JN = NSV_AERBEG,NSV_AEREND
    DO JK=1,JPVEXT
      XSVT(:,:,IKB-JPVEXT,JN) = XSVT(:,:,IKB,JN)
      XSVT(:,:,IKE+JPVEXT,JN) = XSVT(:,:,IKE,JN)
    END DO
  END DO
  !
  IF (YUNIT .EQ. "CON") THEN
    WRITE(KLUOUT,*) "CH_INIT_FIELD_n (ORILAM): converting initial values to mixing ratio"
    DO JN = NSV_AERBEG,NSV_AEREND
      XSVT(:,:,:,JN) = XSVT(:,:,:,JN)/(XRHODREF(:,:,:)*ZDEN2MOL)
    ENDDO
  ELSE
    WRITE(KLUOUT,*)"CH_INIT_FIELD_n (ORILAM): initial values are used as is (mixing ratio)"
  ENDIF
  !
  ENDIF !LORILAM
  !
ENDIF
!
!
DO JN = NSV_CHEMBEG,NSV_CHEMEND
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XSVT(:,:,:,JN))
END DO
DO JN = NSV_AERBEG,NSV_AEREND
  CALL ADD3DFIELD_ll(TZFIELDS_ll, XSVT(:,:,:,JN))
END DO
!
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
!-------------------------------------------------------------------------------
!
!*        3.    INITIALIZE CHEMICAL CONSTANTS
!
CALL CH_INIT_CONST_n(KLUOUT, KVERB)       
!
!-------------------------------------------------------------------------------
!
!*        4.   INITIALIZE AEROSOLS
!              -------------------
!
IF (LORILAM) THEN
  CALL CH_AER_EQM_INIT_n(XSVT(:,:,:,NSV_CHEMBEG:NSV_CHEMEND),&
                         XSVT(:,:,:,NSV_AERBEG:NSV_AEREND),&
                         XM3D,XRHOP3D,XSIG3D,& 
                         XRG3D,XN3D, XRHODREF, XCTOTA3D) 
  DO JN = 1,JPIN
      XM3D(:,:,IKB-JPVEXT,JN) =  XM3D(:,:,IKB,JN)
      XM3D(:,:,IKE+JPVEXT,JN) =  XM3D(:,:,IKE,JN)
  END DO
  !
  DO JN = NSV_CHEMBEG,NSV_CHEMEND
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XSVT(:,:,:,JN))
  END DO
  DO JN = NSV_AERBEG,NSV_AEREND
    CALL ADD3DFIELD_ll(TZFIELDS_ll, XSVT(:,:,:,JN))
  END DO
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
END IF
!
!
!-------------------------------------------------------------------------------
!
!*        5.    INITIALIZE LB IN CASE OF LCH_INIT_FIELD
!               ---------------------------------------
!
IF ((LCH_INIT_FIELD).AND.(CPROGRAM/='DIAG  ').AND.(KMI .EQ. 1)) THEN
  ILBX=SIZE(XLBXSVM,1)
  ILBY=SIZE(XLBYSVM,2)
  IRIMX = INT(ILBX/2)
  IRIMY = INT(ILBY/2)
  DO JN = NSV_CHEMBEG,NSV_CHEMEND
    IF(LWEST_ll() .AND. .NOT. L1D) &
      XLBXSVM(1:IRIMX+1,        :,:,JN) = XSVT(1:IRIMX+1,        :,:,JN)
    IF(LEAST_ll() .AND. .NOT. L1D) &
      XLBXSVM(ILBX-IRIMX:ILBX,:,:,JN)   = XSVT(NIU-IRIMX:NIU,    :,:,JN)
    IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
      XLBYSVM(:,1:IRIMY+1,        :,JN) = XSVT(:,1:IRIMY+1,      :,JN)
    IF(LNORTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
      XLBYSVM(:,ILBY-IRIMY:ILBY,:,JN)   = XSVT(:,NJU-IRIMY:NJU,  :,JN)
  END DO
  IF (LORILAM) THEN
  DO JN = NSV_AERBEG,NSV_AEREND
    IF(LWEST_ll() .AND. .NOT. L1D) &
      XLBXSVM(1:IRIMX+1,        :,:,JN) = XSVT(1:IRIMX+1,        :,:,JN)
    IF(LEAST_ll() .AND. .NOT. L1D) &
      XLBXSVM(ILBX-IRIMX:ILBX,:,:,JN)   = XSVT(NIU-IRIMX:NIU,    :,:,JN)
    IF(LSOUTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
      XLBYSVM(:,1:IRIMY+1,        :,JN) = XSVT(:,1:IRIMY+1,      :,JN)
    IF(LNORTH_ll() .AND. .NOT. L1D .AND. .NOT. L2D) &
      XLBYSVM(:,ILBY-IRIMY:ILBY,:,JN)   = XSVT(:,NJU-IRIMY:NJU,  :,JN)
  END DO
  ENDIF
!
ENDIF
!
!
END SUBROUTINE CH_INIT_FIELD_n
