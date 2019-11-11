!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ####################################
       MODULE MODI_CH_AQUEOUS_TMICKESS
!      ####################################
!
INTERFACE
      SUBROUTINE CH_AQUEOUS_TMICKESS( PTSTEP, PRHODREF, PRHODJ, PRTMIN_AQ, &
                                      PRCT, PRRT, PRCS, PRRS,              &
                                      PCSVT, PCRSVS, PRSVT, PRRSVS         )
!
REAL,                     INTENT(IN)    :: PTSTEP    ! Time step          
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rainwater m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCS    ! cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rainwater m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PCSVT    ! cloud water aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PCRSVS   ! cloud water aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRSVT    ! Rainwater aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRSVS   ! Rainwater aq. species source
!
END SUBROUTINE CH_AQUEOUS_TMICKESS
END INTERFACE
END MODULE MODI_CH_AQUEOUS_TMICKESS
!
!     ######################################################################
      SUBROUTINE CH_AQUEOUS_TMICKESS( PTSTEP, PRHODREF, PRHODJ, PRTMIN_AQ, &
                                      PRCT, PRRT, PRCS, PRRS,              &
                                      PCSVT, PCRSVS, PRSVT, PRRSVS         )
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources
!!    corresponding to collision/coalescence processes (autoconversion + accretion)
!!    for the Kessler cloud microphysics parameterization (see slow_terms)
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!
!!    REFERENCE
!!    ---------
!!      Book1 of the documentation ( routine CH_AQUEOUS_TMICKESS )
!!
!!    AUTHOR
!!    ------
!!      C. Mari J.P. Pinty M. Leriche      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/08
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,      ONLY : JPHEXT,    &! number of horizontal External points
                                 JPVEXT      ! number of vertical External points
USE MODD_CLOUDPAR,        ONLY : XC1RC, XC2RC, &! autoconversion param.
                                 XCEXRA, XCRA, &! accrection param.
                                 XCEXVT         ! constant in the rain drop fall velocity
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,                     INTENT(IN)    :: PTSTEP    ! Time step          
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rainwater m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCS    ! cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rainwater m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PCSVT   ! cloud water aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PCRSVS  ! cloud water aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRSVT   ! Rainwater aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRSVS  ! Rainwater aq. species source
!ch_monitorn.f90
!*       0.2   Declarations of local variables :
!
INTEGER :: JLC, JLR      ! Loop index for cloud water and rainwater aq. species
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !
!
INTEGER :: IMICRO        ! case number of r_x>0 locations
LOGICAL, DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: GMICRO   ! where to compute mic. processes
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRCS   ! Cloud water m.r. source phys.tendency
!REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
!                                :: ZRRS   ! Rain water m.r. source phys. tendency
REAL,    DIMENSION(SIZE(PCRSVS,1),SIZE(PCRSVS,2),SIZE(PCRSVS,3),SIZE(PCRSVS,4)) &
                                :: ZZCRSVS   ! Cloud water aq. species source
REAL,    DIMENSION(SIZE(PRRSVS,1),SIZE(PRRSVS,2),SIZE(PRRSVS,3),SIZE(PRRSVS,4)) &
                                :: ZZRRSVS   ! Rain water aq.  species source
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZW       ! work array
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZCW      ! work array
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRW      ! work array
REAL, DIMENSION(:),   ALLOCATABLE :: ZRCT   ! Cloud water m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRRT   ! Rain water m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZZRCS   ! Cloud water m.r. source
!REAL, DIMENSION(:),   ALLOCATABLE :: ZZRRS   ! Rain water m.r. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCSVT  ! Cloud water aq. species at t
!REAL, DIMENSION(:,:), ALLOCATABLE :: ZRSVT  ! Rain water aq. species at t
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCRSVS ! Cloud water aq. species source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRRSVS ! Rain water aq. species source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZW1  ! Work array 
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZRHODREF, & ! RHO Dry REFerence
                                   ZZW           ! Work array
!
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
!  
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PRCT,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
!!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
ZRCS(:,:,:)  = PRCS(:,:,:) / PRHODJ(:,:,:)
!ZRRS(:,:,:)  = PRRS(:,:,:) / PRHODJ(:,:,:)
!
DO JLC= 1, SIZE(PCRSVS,4)
  ZZCRSVS(:,:,:,JLC) = PCRSVS(:,:,:,JLC) / PRHODJ(:,:,:)
ENDDO
DO JLR= 1, SIZE(PRRSVS,4)
  ZZRRSVS(:,:,:,JLR) = PRRSVS(:,:,:,JLR) / PRHODJ(:,:,:)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       3.     OPTIMIZATION: looking for locations where lwc lwr > min value
!   	        -------------------------------------------------------------
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                                                           &
      (PRCT(IIB:IIE,IJB:IJE,IKB:IKE)>PRTMIN_AQ*1.e3/PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE)) .OR. &
      (PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>PRTMIN_AQ*1.e3/PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE))
!
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 1 ) THEN
  ALLOCATE(ZRCT(IMICRO))
  ALLOCATE(ZRRT(IMICRO))  
  ALLOCATE(ZCSVT(IMICRO,SIZE(PCSVT,4)))
!  ALLOCATE(ZRSVT(IMICRO,SIZE(PRSVT,4)))
  ALLOCATE(ZZRCS(IMICRO)) 
!  ALLOCATE(ZZRRS(IMICRO)) 
  ALLOCATE(ZCRSVS(IMICRO,SIZE(PCRSVS,4)))
  ALLOCATE(ZRRSVS(IMICRO,SIZE(PRRSVS,4)))
  ALLOCATE(ZRHODREF(IMICRO)) 
  ALLOCATE(ZZW(IMICRO))
  ALLOCATE(ZZW1(IMICRO,SIZE(PCSVT,4)))
  DO JL=1,IMICRO   
    ZCSVT(JL,:) = PCSVT(I1(JL),I2(JL),I3(JL),:)
    ZCRSVS(JL,:) = ZZCRSVS(I1(JL),I2(JL),I3(JL),:)
!    ZRSVT(JL,:) = PRSVT(I1(JL),I2(JL),I3(JL),:)
    ZRRSVS(JL,:) = ZZRRSVS(I1(JL),I2(JL),I3(JL),:)
!
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
!
    ZZRCS(JL) = ZRCS(I1(JL),I2(JL),I3(JL))
!    ZZRRS(JL) = ZRRS(I1(JL),I2(JL),I3(JL))
!
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
  ENDDO
!
!
!-------------------------------------------------------------------------------
!
!*       4.     COMPUTES THE SLOW WARM PROCESS SOURCES
!   	        --------------------------------------
!
!*       4.1    compute the autoconversion of r_c for r_r production: RCAUTR
!
    ZZW(:) = 0.0
    ZZW1(:,:) = 0.0
!
  DO JL=1,IMICRO
    IF ( (ZRCT(JL)>0.0) .AND. (ZZRCS(JL)>0.0) ) THEN
      ZZW(JL) = MIN( ZZRCS(JL),XC1RC*MAX( ZRCT(JL)-XC2RC/ZRHODREF(JL),0.0))
!
      ZZW1(JL,:) = ZZW(JL) * ZCSVT(JL,:)/ZRCT(JL)
      ZZW1(JL,:) = MAX(MIN(ZZW1(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
      ZCRSVS(JL,:)  = ZCRSVS(JL,:) - ZZW1(JL,:)
      ZRRSVS(JL,:)  = ZRRSVS(JL,:) + ZZW1(JL,:)
    END IF
  END DO 
!
!*       4.2    compute the accretion of r_c for r_r production: RCACCR
!
    ZZW(:) = 0.0
    ZZW1(:,:) = 0.0
! 
  DO JL = 1,IMICRO
    IF( (ZRCT(JL)>0.0) .AND. (ZRRT(JL)>0.0) .AND. (ZZRCS(JL)>0.0) ) THEN
      ZZW(JL) = MIN( ZZRCS(JL),XCRA * ZRCT(JL) * ZRRT(JL)**XCEXRA       &
                                    * ZRHODREF(JL)**( XCEXRA - XCEXVT ) )
!
      ZZW1(JL,:) = ZZW(JL) * ZCSVT(JL,:)/ZRCT(JL)
      ZZW1(JL,:) = MAX(MIN(ZZW1(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
      ZCRSVS(JL,:)  = ZCRSVS(JL,:) - ZZW1(JL,:)
      ZRRSVS(JL,:)  = ZRRSVS(JL,:) + ZZW1(JL,:)
    END IF
  END DO
!
!
!*       4.3    compute the evaporation of r_r: RREVAV
!
! calculated by the kinetic mass transfer equation (BASIC.f90)
!
!
!-------------------------------------------------------------------------------
!
!*       5.     UNPACK RESULTS AND DEALLOCATE ARRAYS
!   	        ------------------------------------
!
 DO JLC= 1, SIZE(PCRSVS,4)
    ZCW(:,:,:) = ZZCRSVS(:,:,:,JLC)
    ZZCRSVS(:,:,:,JLC) = UNPACK(ZCRSVS(:,JLC), MASK=GMICRO(:,:,:), FIELD=ZCW(:,:,:))
    PCRSVS(:,:,:,JLC) = ZZCRSVS(:,:,:,JLC) * PRHODJ(:,:,:)
 END DO
 DO JLR= 1, SIZE(PRRSVS,4)
    ZRW(:,:,:) = ZZRRSVS(:,:,:,JLR)
    ZZRRSVS(:,:,:,JLR) = UNPACK(ZRRSVS(:,JLR), MASK=GMICRO(:,:,:), FIELD=ZRW(:,:,:))
    PRRSVS(:,:,:,JLR) = ZZRRSVS(:,:,:,JLR) * PRHODJ(:,:,:)
 END DO
!
  DEALLOCATE(ZRHODREF) 
  DEALLOCATE(ZZW) 
  DEALLOCATE(ZZW1)
!  DEALLOCATE(ZZRRS)
  DEALLOCATE(ZZRCS)
  DEALLOCATE(ZCRSVS)
  DEALLOCATE(ZRRSVS)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZCSVT)
!  DEALLOCATE(ZRSVT)
!
END IF
!
!
!-------------------------------------------------------------------------------
!
!
CONTAINS
!
!
  FUNCTION COUNTJV(LTAB,I1,I2,I3) RESULT(IC)
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
!
LOGICAL, DIMENSION(:,:,:) :: LTAB ! Mask
INTEGER, DIMENSION(:) :: I1,I2,I3 ! Used to replace the COUNT and PACK
INTEGER :: JI,JJ,JK,IC
!  
!-------------------------------------------------------------------------------
!
IC = 0
DO JK = 1,SIZE(LTAB,3)
  DO JJ = 1,SIZE(LTAB,2)
    DO JI = 1,SIZE(LTAB,1)
      IF( LTAB(JI,JJ,JK) ) THEN
        IC = IC +1
        I1(IC) = JI
        I2(IC) = JJ
        I3(IC) = JK
      END IF
    END DO
  END DO
END DO
!
END FUNCTION COUNTJV
!
END SUBROUTINE CH_AQUEOUS_TMICKESS
