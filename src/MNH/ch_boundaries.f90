!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!########################
MODULE MODI_CH_BOUNDARIES
!########################
!
INTERFACE
!
      SUBROUTINE CH_BOUNDARIES (HLBCX,HLBCY,                 &
                                PUT,PVT,PSVBT,PSVMIN         )
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)  :: HLBCX,HLBCY  
REAL, DIMENSION(:,:,:),       INTENT(INOUT) :: PSVBT
REAL, DIMENSION(:,:,:),         INTENT(IN) :: PUT,PVT
REAL                                        :: PSVMIN
!
END SUBROUTINE CH_BOUNDARIES
!
END INTERFACE
!
END MODULE MODI_CH_BOUNDARIES
!
!
!     ####################################################################
SUBROUTINE CH_BOUNDARIES (HLBCX,HLBCY,                           &
                          PUT,PVT,PSVBT,PSVMIN                   )
!     ####################################################################
!
!!****  *CH_BOUNDARIES* - routine to prepare the Lateral Boundary Conditions for
!!                 scalar variables at a scalar localization relative to the 
!!                 considered boundary.
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!      Only for 'OPEN' case  Boundary Condition type
!!   
!!    EXTERNAL 
!!    --------  
!!    GET_INDICE_ll  : get physical sub-domain bounds
!!    LWEAST_ll,LEAST_ll,LNORTH_ll,LSOUTH_ll : position functions
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------  
!!      Module MODD_PARAMETERS : 
!!        JPHEXT ,JPVEXT 
!!
!!      Module MODD_CONF :
!!        CCONF
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	P. Tulet    * LA *
!!
!!    MODIFICATIONS
!!    -------------
!!     Original        06/06/00
!!     06/06/00 (C. Mari) embedded into mesonh routines
!!     15/02/01 (P. Tulet) update for MOCAGE lateral boundary conditions
!!     10/02/17 (M. Leriche) prevent negative values
!!	
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!         
USE MODD_PARAMETERS
USE MODD_CONF
USE MODD_GRID_n,    ONLY : XZZ
!
!
USE MODE_ll
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments
!
!
CHARACTER(LEN=4), DIMENSION(2), INTENT(IN)  :: HLBCX,HLBCY  
REAL, DIMENSION(:,:,:),     INTENT(INOUT)   :: PSVBT
REAL, DIMENSION(:,:,:),          INTENT(IN) :: PUT,PVT
REAL                                        :: PSVMIN
!
!
!
!*       0.2   declarations of local variables
!
REAL,    SAVE, DIMENSION(:,:,:), ALLOCATABLE :: ZZZS, ZZZN, ZZZW, ZZZE
INTEGER, SAVE, DIMENSION(:,:,:), ALLOCATABLE :: IZZW, IZZE, IZZS, IZZN
INTEGER :: II,IJ,IK,IKK   ! loop index for LS Scalar Variable
INTEGER :: IIU,IJU,IKU    ! array size in first, second and third dimensions
INTEGER :: IIB,IIE,IJB,IJE! index values for the Beginning or the End
! of the physical domain in x and y directions
INTEGER :: IKB, IKE
LOGICAL, SAVE     :: GFIRSTCALLS = .TRUE.
LOGICAL, SAVE     :: GFIRSTCALLN = .TRUE.
LOGICAL, SAVE     :: GFIRSTCALLW = .TRUE.
LOGICAL, SAVE     :: GFIRSTCALLE = .TRUE.
REAL              :: ZZZ
!
!*       1.    COMPUTE DIMENSIONS OF ARRAYS AND OTHER INDICES:
!              ----------------------------------------------
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PUT,3) - JPVEXT
IIU = SIZE(PUT,1)
IJU = SIZE(PUT,2)
IKU = SIZE(PUT,3)
!
IF (LWEST_ll()) THEN
  IF (.NOT.(ALLOCATED(ZZZW))) ALLOCATE(ZZZW(1,IJU,IKU))
  IF (.NOT.(ALLOCATED(IZZW))) ALLOCATE(IZZW(1,IJU,IKU))
END IF

IF (LEAST_ll()) THEN
  IF (.NOT.(ALLOCATED(ZZZE))) ALLOCATE(ZZZE(1,IJU,IKU))
  IF (.NOT.(ALLOCATED(IZZE))) ALLOCATE(IZZE(1,IJU,IKU))
END IF

IF (LNORTH_ll()) THEN
  IF (.NOT.(ALLOCATED(ZZZN))) ALLOCATE(ZZZN(IIU,1,IKU))
  IF (.NOT.(ALLOCATED(IZZN))) ALLOCATE(IZZN(IIU,1,IKU))
END IF

IF (LSOUTH_ll()) THEN
  IF (.NOT.(ALLOCATED(ZZZS))) ALLOCATE(ZZZS(IIU,1,IKU))
  IF (.NOT.(ALLOCATED(IZZS))) ALLOCATE(IZZS(IIU,1,IKU))
END IF
!
!-------------------------------------------------------------------------------
!
!*             LBC FILLING IN THE X DIRECTION (LEFT WEST SIDE):   
!              ------------------------------------------------
!
IF (LWEST_ll( ) .AND. HLBCX(1)=='OPEN') THEN
  !
  IF (GFIRSTCALLW) THEN
    !
    GFIRSTCALLW = .FALSE.
    ZZZW(1,:,:)= -999.
    IZZW(1,:,:)= 0
    !
    DO IJ= IJB,IJE
      DO IK= IKB,IKE
        ZZZ=XZZ(IIB-1,IJ,IK)
        DO IKK=IKB,IKE
          IF ((ZZZ.GT.XZZ(IIB+1,IJ,IKK)).AND.  &
               (ZZZ.LT.XZZ(IIB+1,IJ,IKK+1))) THEN
            ZZZW(1,IJ,IK)= &
                 (ZZZ - XZZ(IIB+1,IJ,IKK)) / &
                 (XZZ(IIB+1,IJ,IKK+1)-XZZ(IIB+1,IJ,IKK))
            IZZW(1,IJ,IK)= IKK
          END IF
        END DO
      END DO
    END DO
  END IF
  !
  DO IJ=1, IJU
    DO IK= 1,IKU
      IF ( PUT(IIB,IJ,IK) <= 0. ) THEN         !  OUTFLOW condition
        PSVBT(IIB-1,IJ,IK) = &
             MAX(PSVMIN,2.*PSVBT (IIB,IJ,IK) -PSVBT (IIB+1,IJ,IK))
      ELSE                                     !  INFLOW  condition
        IF (ZZZW(1,IJ,IK) > -999.) THEN
          PSVBT(IIB-1,IJ,IK) = MAX(PSVMIN,PSVBT(IIB+1,IJ,IZZW(1,IJ,IK))+&
               (PSVBT(IIB+1,IJ,IZZW(1,IJ,IK)+1)-&
               PSVBT(IIB+1,IJ,IZZW(1,IJ,IK))) * ZZZW(1,IJ,IK))
          !
        ELSE
          PSVBT(IIB-1,IJ,IK) =  PSVBT(IIB+1,IJ,IK)
        END IF
      END IF
    END DO
  END DO
END IF
!-------------------------------------------------------------------------------
!
!*             LBC FILLING IN THE X DIRECTION (RIGHT EAST SIDE): 
!              -------------------------------------------------
!
IF (LEAST_ll( ) .AND. HLBCX(1)=='OPEN') THEN
  ! 
  IF (GFIRSTCALLE) THEN
    !
    GFIRSTCALLE = .FALSE.
    ZZZE(1,:,:)= -999.
    IZZE(1,:,:)= 0
    !
    DO IJ= IJB,IJE
      DO IK= IKB,IKE
        ZZZ=XZZ(IIE+1,IJ,IK)
        DO IKK=IKB,IKE
          IF ((ZZZ.GT.XZZ(IIE-1,IJ,IKK)).AND. &
               (ZZZ.LT.XZZ(IIE-1,IJ,IKK+1))) THEN
            ZZZE(1,IJ,IK)= &
                 (ZZZ - XZZ(IIE-1,IJ,IKK)) &
                 / ( XZZ(IIE-1,IJ,IKK+1)-XZZ(IIE-1,IJ,IKK))
            IZZE(1,IJ,IK)= IKK
          END IF
        END DO
      END DO
    END DO
    !
  END IF
  !
  DO IJ=1,IJU
    DO IK=1,IKU
      IF ( PUT(IIE+1,IJ,IK) >= 0. ) THEN         !  OUTFLOW condition
        PSVBT(IIE+1,IJ,IK) = &
             MAX(PSVMIN,2.*PSVBT (IIE,IJ,IK) -PSVBT (IIE-1,IJ,IK))
      ELSE                                       !  INFLOW  condition
        IF (ZZZE(1,IJ,IK) > -999.) THEN
          PSVBT(IIE+1,IJ,IK) = &
               MAX(PSVMIN,PSVBT(IIE-1,IJ,IZZE(1,IJ,IK))+&
               (PSVBT(IIE-1,IJ,IZZE(1,IJ,IK)+1)-&
               PSVBT(IIE-1,IJ,IZZE(1,IJ,IK))) * ZZZE(1,IJ,IK))
          !
        ELSE
          PSVBT(IIE+1,IJ,IK) =  PSVBT(IIE-1,IJ,IK)
        END IF
      END IF
    END DO
  END DO
END IF
!-------------------------------------------------------------------------------
!
!*             LBC FILLING IN THE Y DIRECTION (BOTTOM SOUTH SIDE): 
!              ---------------------------------------------------
!
IF (LSOUTH_ll( ) .AND. HLBCY(1)=='OPEN') THEN
  !
  !
  IF (GFIRSTCALLS) THEN
    !
    GFIRSTCALLS = .FALSE.
    ZZZS(:,1,:)= -999.
    IZZS(:,1,:)= 0
    !
    DO II= IIB,IIE
      DO IK= IKB,IKE
        ZZZ=XZZ(II,IJB-1,IK)
        DO IKK=IKB,IKE
          IF ((ZZZ.GT.XZZ(II,IJB+1,IKK)).AND.&
               (ZZZ.LT.XZZ(II,IJB+1,IKK+1))) THEN
            ZZZS(II,1,IK)= &
                 (ZZZ - XZZ(II,IJB+1,IKK)) / &
                 (XZZ(II,IJB+1,IKK+1)-XZZ(II,IJB+1,IKK))
            IZZS(II,1,IK)= IKK
          END IF
        END DO
      END DO
    END DO
    !
  END IF
  !
  DO II=1,IIU
    DO IK=1,IKU
      IF ( PVT(II,IJB,IK) <= 0. ) THEN         !  OUTFLOW condition
        PSVBT(II,IJB-1,IK) = &
             MAX(PSVMIN,2.*PSVBT (II,IJB,IK) -PSVBT (II,IJB+1,IK))
      ELSE                                     !  INFLOW  condition
        IF (ZZZS(II,1,IK) > -999.) THEN
          PSVBT(II,IJB-1,IK) = &
               MAX(PSVMIN,PSVBT(II,IJB+1,IZZS(II,1,IK))+&
               (PSVBT(II,IJB+1,IZZS(II,1,IK)+1)-&
               PSVBT(II,IJB+1,IZZS(II,1,IK))) * ZZZS(II,1,IK))
          !
        ELSE
          PSVBT(II,IJB-1,IK) =  PSVBT(II,IJB+1,IK)
        END IF
      END IF
    END DO
  END DO
  !
END IF
!-------------------------------------------------------------------------------
!
!*             LBC FILLING IN THE Y DIRECTION (TOP NORTH SIDE): 
!              ------------------------------------------------
!
IF (LNORTH_ll( ) .AND. HLBCY(2)=='OPEN') THEN
  !
  !
  IF (GFIRSTCALLN) THEN
    !
    GFIRSTCALLN = .FALSE.
    ZZZN(:,1,:)= -999.
    IZZN(:,1,:)= 1
    !
    DO II= IIB,IIE
      DO IK= IKB,IKE
        ZZZ=XZZ(II,IJE+1,IK)
        DO IKK=IKB,IKE
          IF ((ZZZ.GT.XZZ(II,IJE-1,IKK)).AND. &
               (ZZZ.LT.XZZ(II,IJE-1,IKK+1))) THEN
            ZZZN(II,1,IK)= &
                 (ZZZ - XZZ(II,IJE-1,IKK)) / &
                 ( XZZ(II,IJE-1,IKK+1)-XZZ(II,IJE-1,IKK))
            IZZN(II,1,IK)= IKK
          END IF
        END DO
      END DO
    END DO
    !
  END IF
  !
  DO II=1,IIU
    DO IK=1,IKU
      IF ( PVT(II,IJE+1,IK) >= 0. ) THEN         !  OUTFLOW condition
        PSVBT(II,IJE+1,IK) = &
             MAX(PSVMIN,2.*PSVBT (II,IJE,IK) -PSVBT (II,IJE-1,IK))
      ELSE                                       !  INFLOW  condition
        IF (ZZZN(II,1,IK) > -999.) THEN
          PSVBT(II,IJE+1,IK) = &
               MAX(PSVMIN,PSVBT(II,IJE-1,IZZN(II,1,IK))+&
               (PSVBT(II,IJE-1,IZZN(II,1,IK)+1)-&
               PSVBT(II,IJE-1,IZZN(II,1,IK))) * ZZZN(II,1,IK))
          !
        ELSE
          PSVBT(II,IJE+1,IK) =  PSVBT(II,IJE-1,IK)
        END IF
      END IF
    END DO
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CH_BOUNDARIES
