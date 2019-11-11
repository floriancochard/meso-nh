!MNH_LIC Copyright 2000-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################################
      MODULE MODI_UV_TO_ZONAL_AND_MERID
!     #################################
INTERFACE UV_TO_ZONAL_AND_MERID
      SUBROUTINE UV_TO_ZONAL_AND_MERID3D(PU,PV,KGRID,PZC,PMC,TPFILE,TZFIELDS)
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_FIELD, ONLY: TFIELDDATA
!
REAL, DIMENSION(:,:,:),                  INTENT(IN)  :: PU       ! Input U component
REAL, DIMENSION(:,:,:),                  INTENT(IN)  :: PV       ! Input V component
INTEGER,                                 INTENT(IN)  :: KGRID    ! Grid positions of components
REAL, DIMENSION(:,:,:),        OPTIONAL, INTENT(OUT) :: PZC      ! Output U component
REAL, DIMENSION(:,:,:),        OPTIONAL, INTENT(OUT) :: PMC      ! Output V component
TYPE(TFILEDATA),               OPTIONAL, INTENT(IN)  :: TPFILE   ! Output file
TYPE(TFIELDDATA),DIMENSION(2), OPTIONAL, INTENT(IN)  :: TZFIELDS ! Fields characteristics
!
END SUBROUTINE UV_TO_ZONAL_AND_MERID3D
!
      SUBROUTINE UV_TO_ZONAL_AND_MERID2D(PU,PV,KGRID,PZC,PMC,TPFILE,TZFIELDS)
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_FIELD, ONLY: TFIELDDATA
!
REAL, DIMENSION(:,:),                    INTENT(IN)  :: PU       ! Input U component
REAL, DIMENSION(:,:),                    INTENT(IN)  :: PV       ! Input V component
INTEGER,                                 INTENT(IN)  :: KGRID    ! Grid positions of components
REAL, DIMENSION(:,:),          OPTIONAL, INTENT(OUT) :: PZC      ! Output U component
REAL, DIMENSION(:,:),          OPTIONAL, INTENT(OUT) :: PMC      ! Output V component
TYPE(TFILEDATA),               OPTIONAL, INTENT(IN)  :: TPFILE   ! Output file
TYPE(TFIELDDATA),DIMENSION(2), OPTIONAL, INTENT(IN)  :: TZFIELDS ! Fields characteristics
!
END SUBROUTINE UV_TO_ZONAL_AND_MERID2D
!
END INTERFACE
END MODULE MODI_UV_TO_ZONAL_AND_MERID
!
!     ###################################
      MODULE MODI_UV_TO_ZONAL_AND_MERID3D
!     ###################################
INTERFACE
!
      SUBROUTINE UV_TO_ZONAL_AND_MERID3D(PU,PV,KGRID,PZC,PMC,TPFILE,TZFIELDS)
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_FIELD, ONLY: TFIELDDATA
!
REAL, DIMENSION(:,:,:),                  INTENT(IN)  :: PU       ! Input U component
REAL, DIMENSION(:,:,:),                  INTENT(IN)  :: PV       ! Input V component
INTEGER,                                 INTENT(IN)  :: KGRID    ! Grid positions of components
REAL, DIMENSION(:,:,:),        OPTIONAL, INTENT(OUT) :: PZC      ! Output U component
REAL, DIMENSION(:,:,:),        OPTIONAL, INTENT(OUT) :: PMC      ! Output V component
TYPE(TFILEDATA),               OPTIONAL, INTENT(IN)  :: TPFILE   ! Output file
TYPE(TFIELDDATA),DIMENSION(2), OPTIONAL, INTENT(IN)  :: TZFIELDS ! Fields characteristics
!
END SUBROUTINE UV_TO_ZONAL_AND_MERID3D
END INTERFACE
END MODULE MODI_UV_TO_ZONAL_AND_MERID3D
!
!     ##########################################
      SUBROUTINE UV_TO_ZONAL_AND_MERID3D(PU,PV,KGRID,PZC,PMC,TPFILE,TZFIELDS)
!     ##########################################
!
!!****  *UV_TO_ZONAL_AND_MERID* - compute the zonal and meridien components
!!                                of input wind, and return or write them
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      I. Mallet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/11/00
!!      N.Asencio   10/09/03 no pointer for PZC,PMC (no pointer in SHUMAN)
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_DIM_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IO_ll,   ONLY: TFILEDATA,NVERB_ERROR,NVERB_INFO,NVERB_WARNING
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAMETERS
!
USE MODE_FIELD,   ONLY: TFIELDDATA
USE MODE_FMWRIT
USE MODE_MSG
!
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:),                  INTENT(IN)  :: PU       ! Input U component
REAL, DIMENSION(:,:,:),                  INTENT(IN)  :: PV       ! Input V component
INTEGER,                                 INTENT(IN)  :: KGRID    ! Grid positions of components
REAL, DIMENSION(:,:,:),        OPTIONAL, INTENT(OUT) :: PZC      ! Output U component
REAL, DIMENSION(:,:,:),        OPTIONAL, INTENT(OUT) :: PMC      ! Output V component
TYPE(TFILEDATA),               OPTIONAL, INTENT(IN)  :: TPFILE   ! Output file
TYPE(TFIELDDATA),DIMENSION(2), OPTIONAL, INTENT(IN)  :: TZFIELDS ! Fields characteristics
!
!*      0.2    declarations of local variables
!
INTEGER                            :: IKU
REAL                               :: ZRAD_O_DG
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZWORK2
REAL, DIMENSION(:,:,:),ALLOCATABLE :: ZWORK3
REAL, DIMENSION(:,:,:),ALLOCATABLE :: ZZC,ZMC
!
INTEGER           :: IGRID          ! grid indicator
INTEGER           :: ILUOUT         ! logical unit for output listing
!-----------------------------------------------------------------
!
ILUOUT = TLUOUT%NLU
!
IKU=SIZE(PU,3)
ALLOCATE(ZWORK2(SIZE(XLON,1),SIZE(XLON,2)))
ALLOCATE(ZWORK3(SIZE(XLON,1),SIZE(XLON,2),IKU))
!
ALLOCATE(ZZC(SIZE(XLON,1),SIZE(XLON,2),IKU))
ALLOCATE(ZMC(SIZE(XLON,1),SIZE(XLON,2),IKU))
!
ZRAD_O_DG = XPI/180.
IF (LCARTESIAN) THEN                 ! cartesian geometry
  ZWORK2(:,:) = -XBETA *ZRAD_O_DG
ELSE                                 ! conformal projection
  ZWORK2(:,:) = XRPK * (XLON(:,:) -XLON0) * ZRAD_O_DG -(XBETA *ZRAD_O_DG)
END IF
ZWORK3(:,:,:) = SPREAD( ZWORK2(:,:),DIM=3,NCOPIES=IKU )
DEALLOCATE(ZWORK2)
!
! Zonal and Meridien components of wind
!
IF (KGRID==23) THEN
  WRITE(ILUOUT,*) '- zonal and meridian components of winds are computed'
  CALL PRINT_MSG(NVERB_INFO,'GEN','UV_TO_ZONAL_AND_MERID3D','zonal and meridian components of winds are computed')
  WHERE ( (PU(:,:,:) /= XUNDEF) .AND. (PV(:,:,:) /= XUNDEF) ) 
    ZZC(:,:,:) =  PU(:,:,:) *MXM(COS(ZWORK3(:,:,:))) &
              + MYF(MXM(PV(:,:,:))) *MXM(SIN(ZWORK3(:,:,:)))
    ZMC(:,:,:) = - MXF(MYM(PU(:,:,:))) *MYM(SIN(ZWORK3(:,:,:))) &
               + PV(:,:,:) *MYM(COS(ZWORK3(:,:,:)))
  ELSEWHERE
    ZZC(:,:,:) = XUNDEF
    ZMC(:,:,:) = XUNDEF
  END WHERE
ELSE IF (KGRID==11) THEN
  WRITE(ILUOUT,*) '- zonal and meridian components of winds are computed'
  CALL PRINT_MSG(NVERB_INFO,'GEN','UV_TO_ZONAL_AND_MERID3D','zonal and meridian components of winds are computed')
  WHERE ( (PU(:,:,:) /= XUNDEF) .AND. (PV(:,:,:) /= XUNDEF) ) 
    ZZC(:,:,:) = PU(:,:,:) *COS(ZWORK3(:,:,:)) +PV(:,:,:) *SIN(ZWORK3(:,:,:))
    ZMC(:,:,:) = - PU(:,:,:) *SIN(ZWORK3(:,:,:)) +PV(:,:,:) *COS(ZWORK3(:,:,:))
  ELSEWHERE
    ZZC(:,:,:) = XUNDEF
    ZMC(:,:,:) = XUNDEF
  END WHERE
ELSE IF (KGRID==0) THEN
!
! in this case, input winds are ZONal and MERidien
!          and, output ones are in MesoNH grid (mass points)
  WRITE(ILUOUT,*) '- components of winds are replaced in MesoNH grid'
  CALL PRINT_MSG(NVERB_INFO,'GEN','UV_TO_ZONAL_AND_MERID3D','components of winds are replaced in MesoNH grid')
  WHERE ( (PU(:,:,:) /= XUNDEF) .AND. (PV(:,:,:) /= XUNDEF) ) 
    ZZC(:,:,:) = COS(ZWORK3(:,:,:))* PU(:,:,:) - SIN(ZWORK3(:,:,:))* PV(:,:,:)
    ZMC(:,:,:) = SIN(ZWORK3(:,:,:))* PU(:,:,:) + COS(ZWORK3(:,:,:))* PV(:,:,:)
  ELSEWHERE
    ZZC(:,:,:) = XUNDEF
    ZMC(:,:,:) = XUNDEF
  END WHERE
ELSE
  WRITE(ILUOUT,*) '- warning in uv_to_zonal_and_merid: no computation for GRID= ',KGRID
  CALL PRINT_MSG(NVERB_ERROR,'GEN','UV_TO_ZONAL_AND_MERID3D','invalid KGRID value')
  RETURN
END IF
!
IF(PRESENT(TPFILE)) THEN
  IF(.NOT.PRESENT(TZFIELDS)) THEN
    CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID3D','TPFILE provided without TZFIELDS')
    RETURN
  END IF
  !
  IF ( KGRID==23 ) THEN
    IF ( TZFIELDS(1)%NGRID/=2 .OR. TZFIELDS(2)%NGRID/=3 ) THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID3D','inconsistent values for TZFIELDS(x)%NGRID')
    END IF
  ELSE IF ( KGRID==0 .OR. KGRID==11 ) THEN
    IF ( TZFIELDS(1)%NGRID/=1 .OR. TZFIELDS(2)%NGRID/=1 ) THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID3D','inconsistent values for TZFIELDS(x)%NGRID')
    END IF
  END IF
  !
  IF ( TZFIELDS(1)%CDIR/='XY' .OR. TZFIELDS(2)%CDIR/='XY' ) THEN
    CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID3D','inconsistent values for TZFIELDS(x)%HDIR')
  END IF
  !
  CALL IO_WRITE_FIELD(TPFILE,TZFIELDS(1),ZZC(:,:,:))
  CALL IO_WRITE_FIELD(TPFILE,TZFIELDS(2),ZMC(:,:,:))
ELSE IF (PRESENT(PZC).AND.PRESENT(PMC)) THEN
  PZC(:,:,:)=ZZC(:,:,:)
  PMC(:,:,:)=ZMC(:,:,:)
ELSE
  WRITE(ILUOUT,*) '- warning in uv_to_zonal_and_merid3d: bad optional arguments'
  CALL PRINT_MSG(NVERB_WARNING,'GEN','UV_TO_ZONAL_AND_MERID3D','bad optional arguments')
  RETURN
END IF
!
!-------------------------------------------------------------------------------
DEALLOCATE(ZWORK3)
DEALLOCATE(ZZC,ZMC)
!
END SUBROUTINE UV_TO_ZONAL_AND_MERID3D
!
!
!     ##########################################
      SUBROUTINE UV_TO_ZONAL_AND_MERID2D(PU,PV,KGRID,PZC,PMC,TPFILE,TZFIELDS)
!     ##########################################
!
!!****  *UV_TO_ZONAL_AND_MERID* - compute the zonal and meridien components
!!                                of input wind, and return or write them
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      I. Mallet   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/11/00
!!      I. Mallet   11/09/03 call to UV_ZONAL_AND_MERID3D
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_IO_ll,   ONLY: TFILEDATA,NVERB_WARNING
USE MODD_LUNIT_n, ONLY: TLUOUT
!
USE MODE_FIELD,   ONLY: TFIELDDATA
USE MODE_FMWRIT
USE MODE_MSG
!
USE MODI_UV_TO_ZONAL_AND_MERID3D
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:),                    INTENT(IN)  :: PU       ! Input U component
REAL, DIMENSION(:,:),                    INTENT(IN)  :: PV       ! Input V component
INTEGER,                                 INTENT(IN)  :: KGRID    ! Grid positions of components
REAL, DIMENSION(:,:),          OPTIONAL, INTENT(OUT) :: PZC      ! Output U component
REAL, DIMENSION(:,:),          OPTIONAL, INTENT(OUT) :: PMC      ! Output V component
TYPE(TFILEDATA),               OPTIONAL, INTENT(IN)  :: TPFILE   ! Output file
TYPE(TFIELDDATA),DIMENSION(2), OPTIONAL, INTENT(IN)  :: TZFIELDS ! Fields characteristics
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2),1) :: ZU3D,ZV3D
REAL, DIMENSION(SIZE(PU,1),SIZE(PU,2),1) :: ZZC3D,ZMC3D
INTEGER           :: ILUOUT         ! logical unit for output listing
!
!-----------------------------------------------------------------
!
ILUOUT = TLUOUT%NLU
!
ZU3D(:,:,1)=PU(:,:)
ZV3D(:,:,1)=PV(:,:)
!
CALL UV_TO_ZONAL_AND_MERID3D(ZU3D,ZV3D,KGRID,PZC=ZZC3D,PMC=ZMC3D)
!
IF(PRESENT(TPFILE)) THEN
  IF(.NOT.PRESENT(TZFIELDS)) THEN
    CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID2D','TPFILE provided without TZFIELDS')
    RETURN
  END IF
  !
  IF ( KGRID==23 ) THEN
    IF ( TZFIELDS(1)%NGRID/=2 .OR. TZFIELDS(2)%NGRID/=3 ) THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID2D','inconsistent values for TZFIELDS(x)%NGRID')
    END IF
  ELSE IF ( KGRID==0 .OR. KGRID==11 ) THEN
    IF ( TZFIELDS(1)%NGRID/=1 .OR. TZFIELDS(2)%NGRID/=1 ) THEN
      CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID2D','inconsistent values for TZFIELDS(x)%NGRID')
    END IF
  END IF
  !
  IF ( TZFIELDS(1)%CDIR/='XY' .OR. TZFIELDS(2)%CDIR/='XY' ) THEN
    CALL PRINT_MSG(NVERB_WARNING,'IO','UV_TO_ZONAL_AND_MERID2D','inconsistent values for TZFIELDS(x)%HDIR')
  END IF
  !
  CALL IO_WRITE_FIELD(TPFILE,TZFIELDS(1),ZZC3D(:,:,1))
  CALL IO_WRITE_FIELD(TPFILE,TZFIELDS(2),ZMC3D(:,:,1))
ELSE IF (PRESENT(PZC).AND.PRESENT(PMC)) THEN
  PZC(:,:)=ZZC3D(:,:,1)
  PMC(:,:)=ZMC3D(:,:,1)
ELSE
  WRITE(ILUOUT,*) '- warning in uv_to_zonal_and_merid2d: bad optional arguments'
  CALL PRINT_MSG(NVERB_WARNING,'GEN','UV_TO_ZONAL_AND_MERID2D','bad optional arguments')
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE UV_TO_ZONAL_AND_MERID2D
