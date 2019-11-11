!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ################################
       MODULE MODI_CH_AQUEOUS_SEDIM2MOM
!      ################################
!
INTERFACE
      SUBROUTINE CH_AQUEOUS_SEDIM2MOM (KSPLITR, HCLOUD, PTSTEP, PRTMIN_AQ,&
                                       PZZ, PRHODREF, PRHODJ, PRRT, PRRS, &
                                       PCRT, PCRS, PSVT, PRSVS, PINPRR    )
!
CHARACTER(LEN=*),         INTENT(IN)    :: HCLOUD ! kind if cloud
INTEGER,                  INTENT(IN)    :: KSPLITR ! splitting factor
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT    ! Rain water C at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS    ! Rain water C. source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Precip. aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS   ! Precip. aq. species source
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRR  ! instantaneaous precip.
!
END SUBROUTINE CH_AQUEOUS_SEDIM2MOM
END INTERFACE
END MODULE MODI_CH_AQUEOUS_SEDIM2MOM
!
!     #####################################################################
      SUBROUTINE CH_AQUEOUS_SEDIM2MOM (KSPLITR, HCLOUD, PTSTEP, PRTMIN_AQ,&
                                       PZZ, PRHODREF, PRHODJ, PRRT, PRRS, &
                                       PCRT, PCRS, PSVT, PRSVS, PINPRR    )
!     #####################################################################
!
!!****  * -  compute the explicit microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the sedimentation of chemical
!!   species in the raindrops for the C2R2 and C3R5 cloud microphysical schemes.
!!      The sedimentation rates are computed with a time spliting technique: 
!!    an upstream scheme, written as a difference of non-advective fluxes. 
!!    This source term is added to the next coming time step (split-implicit 
!!    process). see rain_c2r2.f90
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
!!      Module MODD_CONF :
!!          CCONF configuration of the model for the first time step
!!
!!    REFERENCE
!!    ---------
!!      Book1 of the documentation ( routine CH_AQUEOUS_SEDIMC2R2 )
!!
!!    AUTHOR
!!    ------
!!      M. Leriche & J.P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/10/08
!!   2014 G.Delautier : remplace MODD_RAIN_C2R2_PARAM par MODD_RAIN_C2R2_KHKO_PARAM
!!   12/15 M.Leriche : compute instantaneous rain at the surface 
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1      
!!   01/16 M. Leriche : Fusion C2R2 and KHKO
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_CST,             ONLY : XRHOLW, XPI
USE MODD_RAIN_C2R2_DESCR, ONLY : XRTMIN, XCTMIN
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
CHARACTER(LEN=*),         INTENT(IN)    :: HCLOUD ! kind if cloud
INTEGER,                  INTENT(IN)    :: KSPLITR ! splitting factor
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step          
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rain water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rain water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCRT    ! Rain water C at t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCRS    ! Rain water C. source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Precip. aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS   ! Precip. aq. species source
REAL, DIMENSION(:,:),     INTENT(OUT)   :: PINPRR  ! instantaneaous precip.
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK,JI,JJ            ! Vertical loop index for the rain sedimentation 
INTEGER :: JN            ! Temporal loop index for the rain sedimentation
INTEGER :: IIB           !  Define the domain where is 
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           ! 
INTEGER :: IJE           !
INTEGER :: IKB           ! 
INTEGER :: IKE           !
!
REAL    :: ZTSPLITR      ! Small time step for rain sedimentation
!
INTEGER :: ISEDIM ! Case number of sedimentation
LOGICAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3)) &
                                :: GSEDIM ! where to compute the SED processes
!REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
!                                :: ZRRS  ! Rain water m.r. source phys.tendency (*dt)
!REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
!                                :: ZCRS  ! Rain water C source phys.tendency
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZW     ! work array
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZMVRR, ZVRR ! sedimentation velocity for rain
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZWSED  ! sedimentation fluxes for chem.
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZWSEDR  ! sedimentation fluxes
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZRR_SEDIM       ! Drain/Dt sur ZTSPLIT
REAL,    DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2),SIZE(PZZ,3))   &
                                :: ZSV_SEDIM_FACT  ! Cumul des Dsv/DT
REAL, DIMENSION(:), ALLOCATABLE :: ZRRT  ! Rain water m.r. at t
REAL, DIMENSION(:), ALLOCATABLE :: ZRHODREF, & ! RHO Dry REFerence
                                   ZZVRR,    & ! sedimentation velocity for rain
                                   ZZW1,ZZW2  ! Work array
!
INTEGER , DIMENSION(SIZE(GSEDIM)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PZZ,3) - JPVEXT
PINPRR(:,:) = 0. ! initialize instantaneous precip.
!
!-------------------------------------------------------------------------------
!
!!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!ZRRS(:,:,:) = PRRS(:,:,:) / PRHODJ(:,:,:)
!ZCRS(:,:,:) = PCRS(:,:,:) / PRHODJ(:,:,:)
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTE THE SEDIMENTATION (RS) SOURCE
!	        -------------------------------------
!
!*       3.1    time splitting loop initialization
!
ZTSPLITR = PTSTEP / FLOAT(KSPLITR)       ! Small time step
!
!
!*       3.2    compute the sedimentation velocities for rain
!
ZMVRR(:,:,:) = 0.
ZVRR(:,:,:) = 0.
WHERE (PCRT(:,:,:) > XCTMIN(3) .and. PRRT(:,:,:)>XRTMIN(3) )
  ZMVRR(:,:,:) = ((3. * PRHODREF(:,:,:)*PRRT(:,:,:))/   &
                  (4. * XPI *XRHOLW*PCRT(:,:,:)))**0.333 ! in m
  ZVRR(:,:,:) = 0.012 * 1.0E6 * ZMVRR(:,:,:) - 0.2 ! velocity for mixing ratio
END WHERE
WHERE (ZVRR(:,:,:) .lt. 0.0)
  ZVRR(:,:,:) = 0.0
END WHERE
!
!*       3.4    compute the fluxes
! 
!  optimization by looking for locations where
!  the precipitating fields are larger than a minimal value only !!!
!
ZSV_SEDIM_FACT(:,:,:) = 1.0
DO JN = 1 , KSPLITR
! 
  GSEDIM(:,:,:) = .FALSE.
  GSEDIM(IIB:IIE,IJB:IJE,IKB:IKE) = PRRT(IIB:IIE,IJB:IJE,IKB:IKE) > XRTMIN(3) .AND. &
                                    PCRT(IIB:IIE,IJB:IJE,IKB:IKE) > XCTMIN(3)
  ISEDIM = COUNTJV( GSEDIM(:,:,:),I1(:),I2(:),I3(:))
  IF( ISEDIM >= 1 ) THEN
    IF( JN==1 ) THEN
      ZW(:,:,:) = 0.0
      DO JK = IKB , IKE
        ZW(:,:,JK) =ZTSPLITR/(PZZ(:,:,JK+1)-PZZ(:,:,JK))
      END DO
    END IF
    ALLOCATE(ZRHODREF(ISEDIM))
    ALLOCATE(ZZVRR(ISEDIM))
    DO JL=1,ISEDIM
      ZRHODREF(JL) =  PRHODREF(I1(JL),I2(JL),I3(JL))
      ZZVRR(JL) = ZVRR(I1(JL),I2(JL),I3(JL))
    ENDDO
    ALLOCATE(ZZW1(ISEDIM)) ; ZZW1(:) = 0.0
!
!*      for rain/drizzle
!
    ZWSED(:,:,:) = 0.0
    ZWSEDR(:,:,:) = 0.0
    ZZW1(:) = ZZVRR(:) * ZRHODREF(:)
    ZWSED(:,:,:) = UNPACK( ZZW1(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
    ZRR_SEDIM(:,:,:) = 0.0
    DO JK = IKB , IKE
      ZRR_SEDIM(:,:,JK) = ZW(:,:,JK)*(ZWSED(:,:,JK+1)-ZWSED(:,:,JK))&
                          /PRHODREF(:,:,JK)
    END DO
!*      compute instantaneous accumulated precipitations
    IF ( JN == 1 ) THEN
      ALLOCATE(ZZW2(ISEDIM)) ; ZZW2(:) = 0.0
      ALLOCATE(ZRRT(ISEDIM))
      DO JL=1,ISEDIM
        ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
      ENDDO
      WHERE(ZRRT(:) > XRTMIN(3) )
        ZZW2(:) = ZZVRR(:) * ZRRT(:) * ZRHODREF(:)
      ENDWHERE
      ZWSEDR(:,:,:) = UNPACK( ZZW2(:),MASK=GSEDIM(:,:,:),FIELD=0.0 )
      DEALLOCATE(ZRRT)
      DEALLOCATE(ZZW2)
      PINPRR(:,:) = ZWSEDR(:,:,IKB)/XRHOLW              ! in m/s
    ENDIF
! end compute accumulated precip.
    DEALLOCATE(ZRHODREF)
    DEALLOCATE(ZZVRR)
    DEALLOCATE(ZZW1)
    ZSV_SEDIM_FACT(:,:,:) =   ZSV_SEDIM_FACT(:,:,:) * (1.0 + ZRR_SEDIM(:,:,:))
!                       (1.0 + ZRR_SEDIM(:,:,:)/MAX(ZZRRS(:,:,:),XRTMIN_AQ))
  END IF      
END DO
!
! Apply the rain sedimentation rate to the WR_xxx aqueous species
!
DO JL= 1, SIZE(PRSVS,4)
  PRSVS(:,:,:,JL) = MAX( 0.0,ZSV_SEDIM_FACT(:,:,:)*PRSVS(:,:,:,JL) )
END DO
!
CONTAINS
!
!-------------------------------------------------------------------------------
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
INTEGER :: IC
INTEGER :: JI,JJ,JK
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
!-------------------------------------------------------------------------------
!
END SUBROUTINE CH_AQUEOUS_SEDIM2MOM
