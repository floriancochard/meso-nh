!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      ############################
       MODULE MODI_CH_AQUEOUS_CHECK
!      ############################
!
INTERFACE
      SUBROUTINE CH_AQUEOUS_CHECK (PTSTEP, PRHODREF, PRHODJ,PRRS, PRSVS, KRRL,  &
                                   KRR, KEQ, KEQAQ, HNAMES, PRTMIN_AQ, OUSECHIC )
!
REAL,                     INTENT(IN)    :: PTSTEP    ! Timestep  
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS    ! water m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS   ! S.V. source
!
INTEGER,                  INTENT(IN)    :: KRRL    ! Number of liq. variables
INTEGER,                  INTENT(IN)    :: KRR     ! Number of water variables
INTEGER,                  INTENT(IN)    :: KEQ     ! Number of chem. spec.
INTEGER,                  INTENT(IN)    :: KEQAQ   ! Number of liq. chem. spec.
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
!
END SUBROUTINE CH_AQUEOUS_CHECK
END INTERFACE
END MODULE MODI_CH_AQUEOUS_CHECK 
!
!     ###########################################################################
      SUBROUTINE CH_AQUEOUS_CHECK (PTSTEP, PRHODREF, PRHODJ,PRRS, PRSVS, KRRL,  &
                                   KRR, KEQ, KEQAQ, HNAMES, PRTMIN_AQ, OUSECHIC )
!     ###########################################################################
!
!!****  * -  Check the coherence between the mixing ratio of water and the
!!           concentrations of aqueous species
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to nullify the concentration of aqueous
!!    species in place where the mixing ratio of the corresponding water
!!    contents are very low. The residual aqueous concentrations are lost.
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
!!      Book1 of the documentation ( routine CH_AQUEOUS_CHECK )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/11/07
!!      21/11/07 (M. Leriche) correct threshold for aqueous phase chemistry
!!      20/09/10 (M. Leriche) add ice phase chemical species
!!      04/11/13 (M. Leriche) add transfer back to the gas phase if evaporation
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,ONLY: JPHEXT,    &! number of horizontal External points
                          JPVEXT      ! number of vertical External points
USE MODD_NSV,       ONLY : NSV_CHACBEG, NSV_CHACEND, NSV_CHICBEG, NSV_CHICEND, &
                           NSV_CHGSBEG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,                     INTENT(IN)    :: PTSTEP    ! Timestep  
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS    ! water m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS   ! S.V. source
!
INTEGER,                  INTENT(IN)    :: KRRL    ! Number of liq. variables
INTEGER,                  INTENT(IN)    :: KRR     ! Number of water variables
INTEGER,                  INTENT(IN)    :: KEQ     ! Number of chem. spec.
INTEGER,                  INTENT(IN)    :: KEQAQ   ! Number of liq. chem. spec.
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR           ! Loop index for the moist variables
INTEGER :: JSV, JSV2     ! Loop index for the aqueous/ice concentrations
!
INTEGER :: INOCLOUD      ! Case number no cloud water
INTEGER :: INORAIN       ! Case number no rainwater
INTEGER :: IWATER        ! Case number aqueous species
INTEGER :: IICE          ! Case number ice phase species
LOGICAL, DIMENSION(SIZE(PRRS,1),SIZE(PRRS,2),SIZE(PRRS,3)) &
                                   :: GNOCLOUD ! where to compute
LOGICAL, DIMENSION(SIZE(PRRS,1),SIZE(PRRS,2),SIZE(PRRS,3)) &
                                   :: GNORAIN ! where to compute
LOGICAL, DIMENSION(SIZE(PRRS,1),SIZE(PRRS,2),SIZE(PRRS,3)) &
                                   :: GWATER ! where to compute
LOGICAL, DIMENSION(SIZE(PRRS,1),SIZE(PRRS,2),SIZE(PRRS,3)) &
                                   :: GICE   ! where to compute
REAL,    DIMENSION(SIZE(PRRS,1),SIZE(PRRS,2),SIZE(PRRS,3),SIZE(PRRS,4)) &
                                   :: ZRRS
REAL,    DIMENSION(:), ALLOCATABLE :: ZWORK, ZWORK2  ! work array
INTEGER, DIMENSION(3)              :: ISV_BEG, ISV_END
!
REAL                               :: ZRTMIN_AQ
!
INTEGER , DIMENSION(SIZE(GNOCLOUD)) :: I1NC,I2NC,I3NC ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GNORAIN)) :: I1NR,I2NR,I3NR ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GWATER)) :: I1W,I2W,I3W ! Used to replace the COUNT
INTEGER , DIMENSION(SIZE(GICE))   :: I1I,I2I,I3I
INTEGER                           :: JL       ! and PACK intrinsics
!
!-------------------------------------------------------------------------------
!
!*       1.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
DO JRR = 2, KRRL+1
  ZRRS(:,:,:,JRR)  = PRRS(:,:,:,JRR) / PRHODJ(:,:,:)
END DO
IF (OUSECHIC) THEN
  DO JRR = KRRL+1, KRR
    ZRRS(:,:,:,JRR)  = PRRS(:,:,:,JRR) / PRHODJ(:,:,:)
  END DO
ENDIF  
!
!-------------------------------------------------------------------------------
!
!*       2.     COMPUTE THE CHECK (RS) SOURCE
!	        -----------------------------
!
!*       2.1    threshold for the aqueous phase species
!
ZRTMIN_AQ = PRTMIN_AQ / PTSTEP
!
!*       2.2    bounds of the aqueous phase species
!
IF( KRRL==1 ) THEN
  ISV_BEG(2) = NSV_CHACBEG
  ISV_END(2) = NSV_CHACEND
ELSE
  ISV_BEG(2) = NSV_CHACBEG
  ISV_BEG(3) = NSV_CHACBEG+KEQAQ/2
  ISV_END(2) = ISV_BEG(3)-1
  ISV_END(3) = NSV_CHACEND
END IF
!
!*      3.    TRANSFER BACK TO THE GAS PHASE IF EVAPORATION
!             ---------------------------------------------
!
GNOCLOUD(:,:,:)=.FALSE.
WHERE(ZRRS(:,:,:,2)<=(ZRTMIN_AQ*1.e3/PRHODREF(:,:,:)))  !cloud
  GNOCLOUD(:,:,:)=.TRUE.
ENd WHERE
INOCLOUD = COUNTJV( GNOCLOUD(:,:,:),I1NC(:),I2NC(:),I3NC(:))
IF (INOCLOUD >=1 ) THEN
  ALLOCATE(ZWORK(INOCLOUD))
  ZWORK(:) = 0.
  ALLOCATE(ZWORK2(INOCLOUD))
  ZWORK2(:) = 0.
  DO JSV = 1, KEQ-KEQAQ  ! gas phase species
    DO JL = 1, INOCLOUD
      ZWORK(JL) = PRSVS(I1NC(JL),I2NC(JL),I3NC(JL),NSV_CHGSBEG-1+JSV)
    ENDDO
    DO JSV2 = KEQ-KEQAQ + 1, KEQ - KEQAQ/2 !cloud
      DO JL = 1, INOCLOUD
        ZWORK2(JL) = MAX(PRSVS(I1NC(JL),I2NC(JL),I3NC(JL),NSV_CHGSBEG-1+JSV2),0.)
      ENDDO
      IF ((TRIM(HNAMES(JSV))) == (TRIM(HNAMES(JSV2)(4:32))).AND.(ANY(ZWORK2(:)>0))) THEN
!        print*,'evaporation of cloud for chemistry'
        ZWORK(:) = ZWORK(:) + ZWORK2(:)   
      ENDIF
    END DO
    PRSVS(:,:,:,NSV_CHGSBEG-1+JSV) = UNPACK( ZWORK(:),MASK=GNOCLOUD(:,:,:), &
                                            FIELD=PRSVS(:,:,:,NSV_CHGSBEG-1+JSV) )
  END DO
  DEALLOCATE(ZWORK)
  DEALLOCATE(ZWORK2)
END IF
IF( KRRL==2 ) THEN
GNORAIN(:,:,:)=.FALSE.
WHERE(ZRRS(:,:,:,3)<=(ZRTMIN_AQ*1.e3/PRHODREF(:,:,:)))  !rain
  GNORAIN(:,:,:)=.TRUE.
ENd WHERE
INORAIN = COUNTJV( GNORAIN(:,:,:),I1NR(:),I2NR(:),I3NR(:))
IF (INORAIN >=1 ) THEN
  ALLOCATE(ZWORK(INORAIN))
  ZWORK(:) = 0.
  ALLOCATE(ZWORK2(INORAIN))
  ZWORK2(:) = 0.
  DO JSV = 1, KEQ-KEQAQ  ! gas phase species
    DO JL = 1, INORAIN
      ZWORK(JL) = PRSVS(I1NR(JL),I2NR(JL),I3NR(JL),NSV_CHGSBEG-1+JSV)
    ENDDO
    DO JSV2 = KEQ-KEQAQ/2 + 1, KEQ !rain
      DO JL = 1, INORAIN
        ZWORK2(JL) = MAX(PRSVS(I1NR(JL),I2NR(JL),I3NR(JL),NSV_CHGSBEG-1+JSV2),0.)
      ENDDO
      IF ((TRIM(HNAMES(JSV))) == (TRIM(HNAMES(JSV2)(4:32))).AND.(ANY(ZWORK2(:)>0.))) THEN
!        print*,'evaporation of rain for chemistry'
        ZWORK(:) = ZWORK(:) + ZWORK2(:)   
      ENDIF
    END DO
    PRSVS(:,:,:,NSV_CHGSBEG-1+JSV) = UNPACK( ZWORK(:),MASK=GNORAIN(:,:,:), &
                                            FIELD=PRSVS(:,:,:,NSV_CHGSBEG-1+JSV) )
  END DO
  DEALLOCATE(ZWORK)
  DEALLOCATE(ZWORK2)
END IF
END IF
!
!*       4.     FILTER OUT THE AQUEOUS SPECIES WHEN MICROPHYSICS<ZRTMIN_AQ
!	        --------------------------------------------------------
!
DO JRR = 2, KRRL+1
  GWATER(:,:,:) = .FALSE.
  WHERE (ZRRS(:,:,:,JRR)>(ZRTMIN_AQ*1.e3/PRHODREF(:,:,:)))
    GWATER(:,:,:)=.TRUE.
  END WHERE
!
  IWATER = COUNTJV( GWATER(:,:,:),I1W(:),I2W(:),I3W(:))
  IF( IWATER >= 1 ) THEN
    ALLOCATE(ZWORK(IWATER))
    DO JSV = ISV_BEG(JRR), ISV_END(JRR)
      DO JL = 1, IWATER
        ZWORK(JL) = PRSVS(I1W(JL),I2W(JL),I3W(JL),JSV)
      END DO
      PRSVS(:,:,:,JSV) = 0.0
      PRSVS(:,:,:,JSV) = UNPACK( ZWORK(:),MASK=GWATER(:,:,:),FIELD=0.0 )
    END DO
    DEALLOCATE(ZWORK)
  ELSE
    DO JSV = ISV_BEG(JRR), ISV_END(JRR)
      PRSVS(:,:,:,JSV) = 0.0
    ENDDO
  END IF
END DO
!
!
!*       5.     FILTER OUT THE ICE PHASE SPECIES WHEN MICROPHYSICS<ZRTMIN_AQ
!	        ------------------------------------------------------------
!
IF (OUSECHIC) THEN
  DO JRR = KRRL+1, KRR
    GICE(:,:,:) = .FALSE.
    WHERE (ZRRS(:,:,:,JRR)>(ZRTMIN_AQ*1.e3/PRHODREF(:,:,:)))
      GICE(:,:,:)=.TRUE.
    END WHERE
  ENDDO
!
  IICE = COUNTJV( GICE(:,:,:),I1I(:),I2I(:),I3I(:))
  IF( IICE >= 1 ) THEN
    ALLOCATE(ZWORK(IICE))
    DO JSV = NSV_CHICBEG, NSV_CHICEND
      DO JL = 1, IICE
        ZWORK(JL) = PRSVS(I1I(JL),I2I(JL),I3I(JL),JSV)
      END DO
      PRSVS(:,:,:,JSV) = 0.0
      PRSVS(:,:,:,JSV) = UNPACK( ZWORK(:),MASK=GICE(:,:,:),FIELD=0.0 )
    END DO
    DEALLOCATE(ZWORK)
  ELSE
    DO JSV = NSV_CHICBEG, NSV_CHICEND
      PRSVS(:,:,:,JSV) = 0.0
    ENDDO
  ENDIF
ENDIF
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
  FUNCTION COUNTJV(GTAB,I1,I2,I3) RESULT(IC)
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!*       0.2  declaration of local variables
!
!
LOGICAL, DIMENSION(:,:,:) :: GTAB ! Mask
INTEGER, DIMENSION(:) :: I1,I2,I3 ! Used to replace the COUNT and PACK
INTEGER :: JI,JJ,JK,IC
!  
!-------------------------------------------------------------------------------
!
IC = 0
DO JK = 1,SIZE(GTAB,3)
  DO JJ = 1,SIZE(GTAB,2)
    DO JI = 1,SIZE(GTAB,1)
      IF( GTAB(JI,JJ,JK) ) THEN
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
END SUBROUTINE CH_AQUEOUS_CHECK 
