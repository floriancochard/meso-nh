!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 surfex 2006/05/18 13:07:25
!-----------------------------------------------------------------
!      #####################
MODULE MODI_DETECT_FIELD_MNH
!      #####################
!
INTERFACE
!
      SUBROUTINE DETECT_FIELD_MNH(HPROGRAM,KI,KJ,PFIELD,OITSHERE)
!
CHARACTER (LEN=6),  INTENT(IN) :: HPROGRAM   ! program
!
INTEGER,            INTENT(IN) :: KI         ! 1st size of array
INTEGER,            INTENT(IN) :: KJ         ! 2nd size of array
REAL, DIMENSION(KI,KJ), INTENT(IN)::PFIELD ! array containing the data field
!
LOGICAL    , INTENT(OUT)         :: OITSHERE  ! T --> PFIELD is non zero somewhere
!

END SUBROUTINE DETECT_FIELD_MNH
!
END INTERFACE
!
END MODULE MODI_DETECT_FIELD_MNH
!
!     ################################################
SUBROUTINE DETECT_FIELD_MNH(HPROGRAM,KI,KJ,PFIELD,OITSHERE)
!     ################################################
!
!!****  *DETECT_FIELD_MNH* - routine to check if a field is non-zero
!!        (MESONH files)
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!    
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!
!!    AUTHOR
!!    ------
!!	S.Malardel       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original      02/2003
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
!
USE MODE_FM
USE MODE_ll
USE MODE_IO_ll

USE MODD_PARAMETERS,     ONLY : XUNDEF, JPHEXT
!
USE MODD_IO_SURF_MNH, ONLY : NMASK, NIU, NJU, NIB, NJB, NIE, NJE

USE MODI_UNPACK_1D_2D
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
CHARACTER (LEN=6),  INTENT(IN) :: HPROGRAM   ! program
!
INTEGER,            INTENT(IN) :: KI         ! 1st size of array
INTEGER,            INTENT(IN) :: KJ         ! 2nd size of array
REAL, DIMENSION(KI,KJ), INTENT(IN)::PFIELD ! array containing the data field
!
LOGICAL    , INTENT(OUT)         :: OITSHERE  ! T --> PFIELD is non zero somewhere
!
!
!
!*       0.2   declarations of local variables
!
!
REAL, DIMENSION(:,:,:),ALLOCATABLE    :: ZWORK      ! 3D array

INTEGER                               :: IPATCH     ! number of tiles
INTEGER                               :: JPATCH     ! loop index
!
INTEGER            :: IIMAX_ll,IJMAX_ll ! Number of points of physical domain
!
INTEGER            :: IRESP         !   error code
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
CALL GET_GLOBALDIMS_ll ( IIMAX_ll,IJMAX_ll)

IPATCH=SIZE(PFIELD,2)

ALLOCATE(ZWORK(NIU,NJU,IPATCH))

      DO JPATCH=1,IPATCH
         ZWORK=XUNDEF
         CALL UNPACK_1D_2D(NMASK,PFIELD(:,JPATCH),ZWORK(NIB:NIE,NJB:NJE,JPATCH))
      END DO

       OITSHERE = MAX_ll( ZWORK(:,:,:) * (XUNDEF-ZWORK(:,:,:))            &
           ,IRESP,1+JPHEXT,1+JPHEXT,1,IIMAX_ll+JPHEXT,IJMAX_ll+JPHEXT,&
           IPATCH              )>0.

DEALLOCATE(ZWORK)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DETECT_FIELD_MNH
