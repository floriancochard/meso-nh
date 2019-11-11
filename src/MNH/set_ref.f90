!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!###################
MODULE MODI_SET_REF
!###################
!
INTERFACE
!
      SUBROUTINE SET_REF(KMI,TPINIFILE,                                    &
                         PZZ,PZHAT,PJ,PDXX,PDYY,HLBCX,HLBCY,               &
                         PREFMASS,PMASS_O_PHI0,PLINMASS,                   &
                         PRHODREF,PTHVREF,PRVREF,PEXNREF,PRHODJ            )
!
USE MODD_IO_ll, ONLY : TFILEDATA
!
INTEGER,                INTENT(IN)  :: KMI       ! Model index 
TYPE(TFILEDATA),        INTENT(IN)  :: TPINIFILE ! Initial file
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PZZ       ! Height  of the w levels
                                                 ! with orography                                          
REAL, DIMENSION(:),     INTENT(IN)  :: PZHAT     ! Height of the w levels
           ! in the transformed space (GCS transf.) or without orography                                         
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! Jacobian 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! metric coefficient dyy
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! type of lateral boundary 
!                                                   ! condition (i=IB, i=IE+1)
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! type of lateral boundary 
                                                    ! condition (j=JB, j=JE+1)  
REAL,                   INTENT(OUT) :: PREFMASS  ! Mass of the ref. atmosphere
                                          !  contained in the simulation domain 
REAL,                   INTENT(OUT) :: PMASS_O_PHI0 ! normalization constant 
                                          ! used in the PHI0 computation   
REAL,                   INTENT(OUT) :: PLINMASS  ! lineic mass through open 
                                                 ! boundaries                                     
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PRHODREF  ! rhod for reference state  
                                                 !  with orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PTHVREF   ! Thetav for reference state                                          
                                                 !  with orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PRVREF    ! mixing ratio for the reference
                                                 ! state with orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNREF   ! Exner function for reference                                     
                                                 !  state  with orography                                     
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PRHODJ    ! rhod J                                           
!                                            
END SUBROUTINE SET_REF
!
END INTERFACE
!
END MODULE MODI_SET_REF
!
!
!     #########################################################################
      SUBROUTINE SET_REF(KMI,TPINIFILE,                                &
                         PZZ,PZHAT,PJ,PDXX,PDYY,HLBCX,HLBCY,           &
                         PREFMASS,PMASS_O_PHI0,PLINMASS,               &
                         PRHODREF,PTHVREF,PRVREF,PEXNREF,PRHODJ        )
!     #########################################################################
!
!!****  *SET_REF* - routine to set reference state for anelastic approximation
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set reference state for anelastic
!     approximation
!
!!**  METHOD
!!    ------
!!      The vertical profiles of thetav, rho (dry) for anelastic reference
!!    state and the Exner function at model top are read in LFIFM file
!!    (TPINIFILE). These vertical profiles do not take
!!    into account the orography. Since these vertical profiles are the same for 
!!    all nested models, they are only read at the first call  by INI_MODEL1 
!!    (i.e. KMI=1). Variables in module MODD_REF are therefore initialized during
!!    the initialization of model 1.
!!      Then, the 3D reference state which takes into account the orography is
!!    deduced from these vertical profiles by a linear interpolation in height
!!    for virtual potential temperature and density.
!!    The Exner function is computed by integration of hydrostatic relation
!!    from model top.
!!      Then, rho J is computed and the total mass of reference atmosphere is 
!!    diagnozed.
!!      
!!      The lineic mass is computed on the faces with open lateral boundary.
!!
!!    EXTERNAL
!!    --------   
!!      FMREAD      : to read data in LFIFM file 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS : contains declaration of parameter variables
!!       
!!         JPHEXT   : Horizontal external points number
!!         JPVEXT   : Vertical external points number
!!
!!      Module MODD_REF  : contains declaration of reference state variables
!!                         without orography
!!
!!        XRHODREFZ: rhod for reference state  without orography
!!        XTHVREFZ : thetav for reference state without orography
!!        XEXNTOP  : Exner function at model top
!!
!!      Module MODD_CST  : contains physical constants
!!      
!!        XRD      : Gaz constant for dry air Rd
!!        XCPD     : Specific heat at constant pressure for dry air Cp
!!        XP00     : Reference pressure 
!!        XG       : gravity constant
!!        XCPD     : specific heat for dry air
!!
!!      Module MODD_CONF   : contains configuration variables
!!      
!!       NVERB  : verbosity level
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine SET_REF)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        30/06/94 
!!      Modification    02/11/94  (J.Stein)  change the extrapolation for thvref
!!      Modification    02/02/95  (J.P.Lafore) Total mass computation to
!!                                 diagnoze the absolute pressure function Phi0
!!      Modification    01/02/95  (J.Stein)  change the extrapolation for 
!!                                 rhodref
!!      Modification    09/02/95  (V.Masson) computation of lineic mass for
!!                                 open boundaries
!!      Modification    30/10/96  (V.Masson) add prints
!!      Modification    02/02/95  (J.P.Lafore) Introduction of 2 anelastic systems:
!!                                  Modified Anelastic Equation and one derived 
!!                                 from Durran (1989), ANE and DUR respectively
!!      Modification    20/10/97  (J.P.Lafore) introduction of 'DAVI' type of lbc
!!      Modification    14/08/98  (V. Ducrocq) //
!!      Modification    14/07/01  (V. MASSON)  LNEUTRAL case
!!                      13/09/01  (J. Stein)  change the option for the
!!                                point under the ground
!!      Modification    03/12/02  (P. Jabouille)  add no thinshell condition
!!      Modification    05/06     Remove the 'DAVI' type of lbc
!!      Modification    07/13     (J.Colin) Special case for LBOUSS=T 
!!      Modification    07/13     (M.Moge) calling UPDATE_HALO_ll on PRHODJ, PRVREF, 
!!                                PRHODREF, PEXNREF, PTHVREF after computation
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_CONF
USE MODD_CST
USE MODD_IO_ll,   ONLY : TFILEDATA
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAMETERS
USE MODD_REF
!
USE MODE_FMREAD
USE MODE_ll
USE MODE_MPPDB
USE MODE_REPRO_SUM
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!  
INTEGER,                INTENT(IN)  :: KMI       ! Model index 
TYPE(TFILEDATA),        INTENT(IN)  :: TPINIFILE ! Initial file
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PZZ       ! Height  of the w levels
                                                 ! with orography                                          
REAL, DIMENSION(:),     INTENT(IN)  :: PZHAT     ! Height of the w levels
           ! in the transformed space (GCS transf.) or without orography                                         
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PJ        ! Jacobian 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX      ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY      ! metric coefficient dyy
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCX ! type of lateral boundary 
!                                                   ! condition (i=IB, i=IE+1)
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN) :: HLBCY ! type of lateral boundary 
                                                    ! condition (j=JB, j=JE+1) 
                     
REAL,                   INTENT(OUT) :: PREFMASS  ! Mass of the ref. atmosphere
                                          !  contained in the simulation domain 
REAL,                   INTENT(OUT) :: PMASS_O_PHI0 ! normalization constant 
                                          ! used in the PHI0 computation   
REAL,                   INTENT(OUT) :: PLINMASS  ! lineic mass through open 
                                                 ! boundaries                                     
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PRHODREF  ! rhod for reference state  
                                                 !  with orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PTHVREF   ! Thetav for reference state                                          
                                                 !  with orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PRVREF    ! mixing ratio for the reference
                                                 ! state with orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PEXNREF   ! Exner function for reference                                     
                                                 !  state  with orography
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PRHODJ    ! rhod J                                           
!
!*       0.2   declarations of local variables
!
INTEGER             :: ILUOUT                    ! Unit number for prints
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZZM 
                                                 ! height of the mass levels
                                                 ! with orography
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) :: ZRHOREF
                                                 ! Reference density
REAL, DIMENSION(SIZE(PZZ,3))    :: ZZHATM        ! height of the mass levels
               ! in the transformed space (GCS transf.) or without orography 
!
INTEGER             :: IIU        ! Upper dimension in x direction
INTEGER             :: IJU        ! Upper dimension in y direction
INTEGER             :: IKU        ! Upper dimension in z direction
INTEGER             :: IIB        ! indice I Beginning in x direction
INTEGER             :: IJB        ! indice J Beginning in y direction
INTEGER             :: IKB        ! indice K Beginning in z direction
INTEGER             :: IIE        ! indice I End       in x direction 
INTEGER             :: IJE        ! indice J End       in y direction 
INTEGER             :: IKE        ! indice K End       in z direction 
INTEGER             :: JI         ! Loop index in x direction
INTEGER             :: JJ         ! Loop index in y direction      
INTEGER             :: JK         ! Loop index in z direction
INTEGER             :: JKS        ! Loop index
INTEGER             :: IKS        ! index of 1D level just above 3d level at I,J
INTEGER             :: JKLOOP     ! Loop index 
INTEGER             :: IINFO_ll   ! return status of the // routines
REAL       :: ZGSCPD              ! = g/Cpd   
REAL       :: ZCVD_O_RD           ! = Cvd /  Rd
REAL       :: ZCVD_O_RDCPD        ! = Cvd /  (Rd*Cpd)
REAL       :: ZDZ1SDZ,ZDZ2SDZ     ! working arrays
REAL       :: ZD1                 ! DELTA1 (switch 0/1) for thinshell approximation
!JUAN16
REAL, ALLOCATABLE, DIMENSION (:,:) :: ZREFMASS_2D , ZMASS_O_PHI0_2D   
REAL, ALLOCATABLE, DIMENSION (:,:) :: ZLINMASS_W_2D , ZLINMASS_E_2D ,  ZLINMASS_S_2D ,  ZLINMASS_N_2D
!REAL                              :: ZREFMASS , ZMASS_O_PHI0   , ZLINMASS     ! total leak of mass
!JUAN16
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()   ! list of fields to exchange
!
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES:
!              ----------------------------------------------
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IKU = SIZE(PEXNREF,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
ILUOUT = TLUOUT%NLU
!
!*       2.    READ REFERENCE STATE WITHOUT OROGRAPHY IN LFIFM FILE
!              ----------------------------------------------------
!
IF (KMI == 1) THEN
  CALL IO_READ_FIELD(TPINIFILE,'RHOREFZ',XRHODREFZ)
  CALL IO_READ_FIELD(TPINIFILE,'THVREFZ',XTHVREFZ)
  CALL IO_READ_FIELD(TPINIFILE,'EXNTOP', XEXNTOP)
!
  LNEUTRAL=.FALSE.
  IF (MAXVAL(XTHVREFZ(IKB:IKE))-MINVAL(XTHVREFZ(IKB:IKE)) < 1.E-10) LNEUTRAL=.TRUE.
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.    SET REFERENCE STATE WITH OROGRAPHY 
!              ----------------------------------
!
!
!*       3.1    Compute  level and  height of mass position
! 
DO JK = 1,IKU-1
  ZZM(:,:,JK) = 0.5*(PZZ(:,:,JK) + PZZ(:,:,JK+1))
  ZZHATM(JK)  = 0.5*(PZHAT(JK)+PZHAT(JK+1))
END DO
ZZHATM(IKU) = 2.* PZHAT(IKU) -ZZHATM(IKU-1)
ZZM(:,:,IKU)    = 2.* PZZ(:,:,IKU)   -ZZM(:,:,IKU-1)
! ZZM(:,:,IKU) is always smaller than or equal ZZHATM(IKU)
!
!
CALL MPPDB_CHECK3D(ZZM,"SET_REF::ZZM",PRECISION)
!
!*       3.2    Interpolation 
!  
DO JI = 1,SIZE(PZZ,1)
  DO JJ = 1,SIZE(PZZ,2)
!
    DO JK = 1,IKU
!
      IF (ZZM(JI,JJ,JK) >= ZZHATM(IKU)) THEN     ! copy out when  
        PTHVREF(JI,JJ,JK) =  XTHVREFZ(IKU)       ! ZZM(IKU)= ZZHATM(IKU)  
        PRHODREF(JI,JJ,JK) =  XRHODREFZ(IKU)     ! (in case zs=0.)
! 
      ELSE              ! search levels on the mass grid without orography
        IF (ZZM(JI,JJ,JK) < ZZHATM(2)) THEN
           IKS=3
        ELSE
          SEARCH : DO JKS = 3,IKU
            IF((ZZM(JI,JJ,JK) >= ZZHATM(JKS-1)).AND.(ZZM(JI,JJ,JK) < ZZHATM(JKS))) &
            THEN          ! interpolation with the values on the grid without
                          ! orography
              IKS=JKS
              EXIT SEARCH
            END IF
          END DO SEARCH
        END IF
        ZDZ1SDZ = (ZZM(JI,JJ,JK)-ZZHATM(IKS-1)) / (ZZHATM(IKS)-ZZHATM(IKS-1))
        ZDZ2SDZ = 1. - ZDZ1SDZ
        PTHVREF(JI,JJ,JK) = ( ZDZ1SDZ* XTHVREFZ(IKS) )  &
                          + (ZDZ2SDZ* XTHVREFZ(IKS-1) )      
        PRHODREF(JI,JJ,JK)= ( ZDZ1SDZ* XRHODREFZ(IKS) ) &
                          + (ZDZ2SDZ* XRHODREFZ(IKS-1) )
      END IF
    END DO
  END DO
END DO
!
!   change the extrapolation option for the thvref field to be consistent with
!   the extrapolation option for the flottability at the ground and for rhodref
!   to be consistent with the extrapolation to compute a divergence 
PTHVREF(:,:,IKB-1) = PTHVREF(:,:,IKB)
PRHODREF(:,:,IKB-1) = PRHODREF(:,:,IKB)

CALL MPPDB_CHECK3D(PTHVREF,"SET_REF::PTHVREF",PRECISION)
CALL MPPDB_CHECK3D(PRHODREF,"SET_REF::PRHODREF",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       4.    COMPUTE EXNER FUNCTION
!              ----------------------
!
IF (LCARTESIAN .OR. LTHINSHELL) THEN
  ZD1=0.
ELSE
  ZD1=1.
ENDIF
!
ZGSCPD = XG/XCPD
!
PEXNREF(:,:,IKE)=(XEXNTOP*(1.+ZD1*2./7.*(PZZ(:,:,IKE+1)-ZZM(:,:,IKE))/  &
                                (XRADIUS+(PZZ(:,:,IKE+1)+ZZM(:,:,IKE))/2.))  &
  + ZGSCPD/PTHVREF(:,:,IKE)*(PZZ(:,:,IKE+1)-ZZM(:,:,IKE)))/ &
(1.-ZD1*2./7.*(PZZ(:,:,IKE+1)-ZZM(:,:,IKE))/(XRADIUS+(PZZ(:,:,IKE+1)+ZZM(:,:,IKE))/2.))
!
DO JK = IKE-1, 1, -1
  PEXNREF(:,:,JK)=(PEXNREF(:,:,JK+1)*(1.+ZD1*2./7.*(ZZM(:,:,JK+1) -ZZM(:,:,JK))/ &
                                           (XRADIUS+PZZ(:,:,JK+1)))+ &
  2.*ZGSCPD/(PTHVREF(:,:,JK+1)+PTHVREF(:,:,JK))*(ZZM(:,:,JK+1) -ZZM(:,:,JK)))/&
  (1.-ZD1*2./7.*(ZZM(:,:,JK+1) -ZZM(:,:,JK))/(XRADIUS+PZZ(:,:,JK+1)))
END DO
!
DO JK = IKE+1, IKU
  PEXNREF(:,:,JK)=(PEXNREF(:,:,JK-1)*(1.+ZD1*2./7.*(ZZM(:,:,JK-1) -ZZM(:,:,JK))/  &
                                           (XRADIUS+PZZ(:,:,JK)))+ &
  2.*ZGSCPD/(PTHVREF(:,:,JK-1)+PTHVREF(:,:,JK))*(ZZM(:,:,JK-1) -ZZM(:,:,JK)))/&
  (1.-ZD1*2./7.*(ZZM(:,:,JK-1) -ZZM(:,:,JK))/ (XRADIUS+PZZ(:,:,JK)))
END DO
!
CALL MPPDB_CHECK3D(PEXNREF,"SET_REF::PEXNREF",PRECISION)
!-------------------------------------------------------------------------------
!
!*       5.    SET RHODJ  AND  REFERENCE DENSITY
!              ---------------------------------
!
!
ZCVD_O_RD = (XCPD / XRD) - 1.
IF (LBOUSS) THEN
  ZRHOREF(:,:,:) = PRHODREF(:,:,:)
ELSE
  ZRHOREF(:,:,:) = PEXNREF(:,:,:) ** ZCVD_O_RD * XP00 / ( XRD * PTHVREF(:,:,:) )
  ZRHOREF(:,:,1)=ZRHOREF(:,:,2)  ! this avoids to obtain erroneous values for
END IF
                               ! rv at this last point
!
IF ( CEQNSYS == 'DUR' ) THEN
  IF ( SIZE(PRVREF,1) == 0 ) THEN
    PRHODJ(:,:,:) =  PRHODREF(:,:,:)* PJ(:,:,:) * PTHVREF(:,:,:)  &
                   / XTH00
  ELSE
    PRVREF(:,:,:) = ( ZRHOREF(:,:,:)/PRHODREF(:,:,:) ) - 1.
    PRHODJ(:,:,:) =  PRHODREF(:,:,:)* PJ(:,:,:) * PTHVREF(:,:,:)   &
                   * (1. + PRVREF(:,:,:)) / XTH00
  END IF
ELSEIF ( CEQNSYS == 'MAE' .OR. CEQNSYS == 'LHE' ) THEN
  PRHODJ(:,:,:) =  PRHODREF(:,:,:)* PJ(:,:,:)
END IF
!
! update halo of PRHODJ and PRVREF for future use ( notably in anel_balance_n )
!
NULLIFY( TZFIELDS_ll )
CALL ADD3DFIELD_ll(TZFIELDS_ll,PRHODJ)
IF ( SIZE(PRVREF,1) /= 0 ) CALL ADD3DFIELD_ll(TZFIELDS_ll,PRVREF)
CALL ADD3DFIELD_ll(TZFIELDS_ll,PRHODREF)
CALL ADD3DFIELD_ll(TZFIELDS_ll,PEXNREF)
CALL ADD3DFIELD_ll(TZFIELDS_ll,PTHVREF)
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
CALL MPPDB_CHECK3D(ZRHOREF,"SET_REF::ZRHOREF",PRECISION)
IF ( SIZE(PRVREF,1) /= 0 ) CALL MPPDB_CHECK3D(PRVREF,"SET_REF::PRVREF",PRECISION)
CALL MPPDB_CHECK3D(PRHODJ,"SET_REF::PRHODJ",PRECISION)


!
!*       6.     COMPUTES THE TOTAL MASS OF REFERENCE ATMOSPHERE   
!	        -----------------------------------------------
!
IF (CEQNSYS == "LHE" ) THEN
  ZCVD_O_RDCPD = ZCVD_O_RD / XCPD
  !
  ALLOCATE(ZREFMASS_2D(IIB:IIE,IJB:IJE))
  ALLOCATE(ZMASS_O_PHI0_2D(IIB:IIE,IJB:IJE))
  ZREFMASS_2D     = 0.
  ZMASS_O_PHI0_2D = 0.
  DO JK = IKB,IKE
    DO JJ = IJB,IJE
      DO JI = IIB,IIE
        ZREFMASS_2D(JI,JJ)  = ZREFMASS_2D(JI,JJ) + ZRHOREF (JI,JJ,JK) * PJ(JI,JJ,JK) ! Reference density
        ZMASS_O_PHI0_2D(JI,JJ) = ZMASS_O_PHI0_2D(JI,JJ) + ZRHOREF(JI,JJ,JK) / PTHVREF(JI,JJ,JK) &
                          * ZCVD_O_RDCPD * PJ(JI,JJ,JK) / PEXNREF(JI,JJ,JK)     
      END DO
    END DO
  END DO
! 
!JUAN16 
!!$  CALL REDUCESUM_ll(ZREFMASS,IINFO_ll)
!!$  CALL REDUCESUM_ll(ZMASS_O_PHI0,IINFO_ll)
  PREFMASS     =  SUM_DD_R2_ll(ZREFMASS_2D)
  PMASS_O_PHI0 =  SUM_DD_R2_ll(ZMASS_O_PHI0_2D)
!JUAN16
!
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       6.    COMPUTATION OF LINEIC MASS
!              --------------------------
!
PLINMASS=0.
!
IF ( HLBCX(1)=='OPEN' ) THEN
   ALLOCATE(ZLINMASS_W_2D(IIB:IIB,IJB:IJE))
   ZLINMASS_W_2D = 0.0
   IF  (LWEST_ll(HSPLITTING='B')) THEN
      DO JJ=IJB,IJE
         DO JK=IKB,IKE
            ZLINMASS_W_2D(IIB,JJ)=ZLINMASS_W_2D(IIB,JJ)+1./PDXX(IIB,JJ,JK) &
                 *0.5*(PRHODJ(IIB,JJ,JK)+PRHODJ(IIB-1,JJ,JK))
         ENDDO
      ENDDO
   ENDIF
   PLINMASS =         SUM_DD_R2_ll(ZLINMASS_W_2D)   
!
   ALLOCATE(ZLINMASS_E_2D(IIE+1:IIE+1,IJB:IJE))
   ZLINMASS_E_2D = 0.0
   IF (LEAST_ll(HSPLITTING='B')) THEN
      DO JJ=IJB,IJE
         DO JK=IKB,IKE
            ZLINMASS_E_2D(IIE+1,JJ)=ZLINMASS_E_2D(IIE+1,JJ)+1./PDXX(IIE+1,JJ,JK) &
                 *0.5*(PRHODJ(IIE+1,JJ,JK)+PRHODJ(IIE,JJ,JK))
         ENDDO
      ENDDO
   ENDIF
   PLINMASS = PLINMASS +  SUM_DD_R2_ll(ZLINMASS_E_2D)
!
ENDIF
IF ( HLBCY(1)=='OPEN' ) THEN
   ALLOCATE( ZLINMASS_S_2D(IIB:IIE,IJB:IJB))
   ZLINMASS_S_2D = 0.0
   IF (LSOUTH_ll(HSPLITTING='B')) THEN
      DO JI=IIB,IIE
         DO JK=IKB,IKE
            ZLINMASS_S_2D(JI,IJB)=ZLINMASS_S_2D(JI,IJB)+1./PDYY(JI,IJB,JK) &
                 *0.5*(PRHODJ(JI,IJB,JK)+PRHODJ(JI,IJB-1,JK))
         ENDDO
      ENDDO
   ENDIF
   PLINMASS = PLINMASS +  SUM_DD_R2_ll(ZLINMASS_S_2D)
!
   ALLOCATE( ZLINMASS_N_2D(IIB:IIE,IJE+1:IJE+1))  
   ZLINMASS_N_2D = 0.0
   IF (LNORTH_ll(HSPLITTING='B')) THEN
      DO JI=IIB,IIE
         DO JK=IKB,IKE
            ZLINMASS_N_2D(JI,IJE+1)=ZLINMASS_N_2D(JI,IJE+1)+1./PDYY(JI,IJE+1,JK) &
                 *0.5*(PRHODJ(JI,IJE+1,JK)+PRHODJ(JI,IJE,JK))
         ENDDO
      ENDDO
   ENDIF
   PLINMASS = PLINMASS +  SUM_DD_R2_ll(ZLINMASS_N_2D)
!
END IF
!
CALL MPPDB_CHECK3D(PRHODREF,"SET_REF::PRHODREF",PRECISION)
CALL MPPDB_CHECK3D(PTHVREF,"SET_REF::PTHVREF",PRECISION)
CALL MPPDB_CHECK3D(PRVREF,"SET_REF::PRVREF",PRECISION)
CALL MPPDB_CHECK3D(PEXNREF,"SET_REF::PEXNREF",PRECISION)
CALL MPPDB_CHECK3D(PRHODJ,"SET_REF::PRHODJ",PRECISION)
!
!-------------------------------------------------------------------------------
!
!*       7.    PRINT ON OUTPUT-LISTING
!              -----------------------
!
IF(NVERB >= 5 ) THEN                               !Value control
  WRITE(ILUOUT,*) 'SET_REF : PLINMASS     = ',PLINMASS
END IF
!
IF(NVERB >= 10) THEN                               !Value control
!
  WRITE(ILUOUT,*) 'SET_REF: XTHVREFZ values:'
  WRITE(ILUOUT,*)  XTHVREFZ
!
  WRITE(ILUOUT,*) 'SET_REF: XRHODREFZ values:'
  WRITE(ILUOUT,*)  XRHODREFZ
!
  WRITE(ILUOUT,*) 'SET_REF: XEXNTOP'
  WRITE(ILUOUT,*)  XEXNTOP
!
  WRITE(ILUOUT,*) 'SET_REF: Some PTHVREF values:'
  DO JKLOOP=1,IKU,5
    WRITE(ILUOUT,*) PTHVREF(1,1,JKLOOP),PTHVREF(IIU/2,IJU/2,JKLOOP), &
     PTHVREF(IIU,IJU,JKLOOP)  
  END DO
!
  WRITE(ILUOUT,*) 'SET_REF: Some PRHODREF values:'
  DO JKLOOP=1,IKU,5
    WRITE(ILUOUT,*) PRHODREF(1,1,JKLOOP),PRHODREF(IIU/2,IJU/2,JKLOOP), &
     PRHODREF(IIU,IJU,JKLOOP)  
  END DO
  WRITE(ILUOUT,*) 'SET_REF: Some PEXNREF values:'
  DO JKLOOP=1,IKU,5
    WRITE(ILUOUT,*) PEXNREF(1,1,JKLOOP),PEXNREF(IIU/2,IJU/2,JKLOOP), &
     PEXNREF(IIU,IJU,JKLOOP)  
  END DO
  WRITE(ILUOUT,*) 'SET_REF: Some PRHODJ values:'
  DO JKLOOP=1,IKU,5
    WRITE(ILUOUT,*) PRHODJ(1,1,JKLOOP),PRHODJ(IIU/2,IJU/2,JKLOOP), &
     PRHODJ(IIU,IJU,JKLOOP)  
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_REF
