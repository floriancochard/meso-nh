!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##########
MODULE MODI_MNH_OASIS_DEFINE
!     ##########
!
INTERFACE 
!
      SUBROUTINE MNH_OASIS_DEFINE(HPROGRAM,IP)
!
      CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM   ! program calling surf. schemes
      INTEGER,            INTENT(IN)  :: IP         ! nb current proc.
!
      END SUBROUTINE MNH_OASIS_DEFINE
!
END INTERFACE
!
END MODULE MODI_MNH_OASIS_DEFINE
!
!     ####################################################################
SUBROUTINE MNH_OASIS_DEFINE(HPROGRAM,IP)
!     ####################################################################
!
!
!!****  *MNH_OASIS_DEFINE*
!!
!!    PURPOSE
!!    -------
!!    Define the mpi partition for OASIS coupling
!!
!!**  METHOD
!!    ------
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
!!    AUTHOR
!!    ------
!!	J. Pianezze   *LPO*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2014
!!                  11/2016  Correction WENO5
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
#ifdef CPLOASIS
USE MODI_SFX_OASIS_DEFINE
USE MODD_DIM_n, ONLY : NIMAX, NJMAX, NIMAX_ll, NJMAX_ll
USE MOD_OASIS
USE MODD_MNH_SURFEX_n
#endif
!
USE MODD_CONF, ONLY : NVERB,NHALO
USE MODD_PARAMETERS, ONLY : JPHEXT
!
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM   ! program calling surf. schemes
INTEGER,            INTENT(IN)  :: IP         ! nb current proc.
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
#ifdef CPLOASIS
INTEGER, DIMENSION(:), ALLOCATABLE :: IPARAL
INTEGER, DIMENSION(:), ALLOCATABLE :: ISEG_SIZE
INTEGER, DIMENSION(:), ALLOCATABLE :: ISEG_OFFSET
!
INTEGER :: JI, JSEG
INTEGER :: ISEGMENT, INPAR, INPTS
INTEGER :: IIOR, IJOR
#endif
!
!-------------------------------------------------------------------------------
#ifdef CPLOASIS
!-------------------------------------------------------------------------------
!
CALL GET_OR_ll('B',IIOR,IJOR)
!
!*       1.     Define ORANGE parallel partitions:
!               ----------------------------------
!
! Number of segments for this proc
!
ISEGMENT=NJMAX
!  
INPAR=2+2*ISEGMENT
!
! Local offset and extent for this proc
!
ALLOCATE(ISEG_SIZE  (ISEGMENT))
ALLOCATE(ISEG_OFFSET(ISEGMENT))
ALLOCATE(IPARAL(INPAR))
!
! OASIS orange partition
!
IPARAL(CLIM_STRATEGY) = CLIM_ORANGE
!
! Number of proc segments for OASIS
!
IPARAL(2) = ISEGMENT
!
! Local offset and extent for OASIS
!
JI=2
INPTS=0
DO JSEG=1,ISEGMENT
   JI=JI+1
   IF (LWEST_ll() .AND. LSOUTH_ll()) THEN
     IPARAL(JI) = (IIOR - JPHEXT) + NIMAX_ll*(IJOR - JPHEXT) + NIMAX_ll*(JSEG - 1)
   ELSE IF (LWEST_ll() .AND. .NOT. LSOUTH_ll()) THEN
     IPARAL(JI) = (IIOR - JPHEXT) + NIMAX_ll*(IJOR + NHALO - 2*JPHEXT) + NIMAX_ll*(JSEG - 1)
   ELSE IF (LSOUTH_ll() .AND. .NOT. LWEST_ll()) THEN
     IPARAL(JI) = (IIOR + NHALO - 2*JPHEXT) + NIMAX_ll*(IJOR - JPHEXT) + NIMAX_ll*(JSEG - 1)
   ELSE IF (.NOT. LSOUTH_ll() .AND. .NOT. LWEST_ll()) THEN
     IPARAL(JI) = (IIOR + NHALO - 2*JPHEXT) + NIMAX_ll*(IJOR + NHALO - 2*JPHEXT) + NIMAX_ll*(JSEG - 1)
   END IF

   JI=JI+1
   IPARAL(JI) = NIMAX
ENDDO
!  
INPTS=NIMAX*NJMAX
!
DEALLOCATE(ISEG_SIZE  )
DEALLOCATE(ISEG_OFFSET)
!
!
!*       2.     Put definitions for exchange of coupling fields :
!               -------------------------------------------------
!
CALL SFX_OASIS_DEFINE(YSURF_CUR%IM%O, YSURF_CUR%U, &
                      HPROGRAM,INPTS, IPARAL        )
!
DEALLOCATE(IPARAL)
!
!-------------------------------------------------------------------------------
#endif
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNH_OASIS_DEFINE
