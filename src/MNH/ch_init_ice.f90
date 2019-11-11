!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
! $Source: /home/cvsroot/MNH-V5-1-4/src/SRC_CHIMAQ/ch_init_icen.f90
!-----------------------------------------------------------------
!     ############################
      MODULE MODI_CH_INIT_ICE
!     ############################
!
INTERFACE
!
      SUBROUTINE CH_INIT_ICE(OUSECHIC, OCH_RET_ICE,      &
                             HNAMES, HICNAMES, KEQ, KEQAQ)
!
INTEGER,                  INTENT(IN)    :: KEQ   ! Number of chem. spec.
INTEGER,                  INTENT(IN)    :: KEQAQ   ! Number of liq. chem. spec.
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
LOGICAL,                  INTENT(IN)    :: OCH_RET_ICE ! flag for retention in ice
!
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES ! name of chem. species
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HICNAMES ! name of ice chem. species
!
END SUBROUTINE CH_INIT_ICE
!
END INTERFACE
!
END MODULE MODI_CH_INIT_ICE
!
!
!      ######################################################
       SUBROUTINE CH_INIT_ICE(OUSECHIC, OCH_RET_ICE,      &
                              HNAMES, HICNAMES, KEQ, KEQAQ)
!      ######################################################
!!
!!***  *CH_INIT_ICE*
!!
!!    PURPOSE
!!    -------
!      Initialize module aqueous ice chemistry (index, )
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCES
!!    ----------
!!    MesoNH-chemistry book 3
!!
!!    AUTHOR
!!    ------
!!    M. Leriche, P. Tulet
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 11/12/15
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    None
!----------------------------------------------------------------------------
!
USE MODD_CH_ICE_n, ONLY : NINDEXGI, NINDEXWI, NINDEXWG
USE MODD_NSV, ONLY : NSV_CHGSBEG, NSV_CHGSEND, &
                     NSV_CHACBEG, NSV_CHACEND, &
                     NSV_CHICBEG, NSV_CHICEND
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN)    :: KEQ   ! Number of chem. spec.
INTEGER,                  INTENT(IN)    :: KEQAQ   ! Number of liq. chem. spec.
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
LOGICAL,                  INTENT(IN)    :: OCH_RET_ICE ! flag for retention in ice
!
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES ! name of chem. species
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HICNAMES ! name of ice chem. species
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JI,JJ, JLI, JLW, JLG   ! loop counters
!----------------------------------------------------------------------------
!
!  PREPARE INDEX ARRAY FOR ICE PHASE CHEMISTRY
!  -------------------------------------------
!

ALLOCATE (NINDEXGI(NSV_CHICEND-NSV_CHICBEG+1))
ALLOCATE (NINDEXWI(NSV_CHICEND-NSV_CHICBEG+1))
ALLOCATE (NINDEXWG(NSV_CHACEND-NSV_CHACBEG+KEQAQ/2+1))
NINDEXGI(:) = 0
NINDEXWI(:) = 0
NINDEXWG(:) = 0
IF (OUSECHIC) THEN
  DO JLI = 1, NSV_CHICEND-NSV_CHICBEG+1
    DO JLG = 1, NSV_CHGSEND-NSV_CHGSBEG+1
      IF ( TRIM(HICNAMES(JLI)(4:32)) == TRIM(HNAMES(JLG)) ) THEN
         NINDEXGI(JLI) = JLG
         EXIT
      ENDIF
    ENDDO
    DO JLW = KEQ-KEQAQ+1, KEQ-KEQAQ/2  ! loop over cloud chem. species
      IF ( TRIM(HICNAMES(JLI)(4:32)) == TRIM(HNAMES(JLW)(4:32))) THEN
        NINDEXWI(JLI) = JLW - (KEQ-KEQAQ)
        EXIT
      ENDIF
    ENDDO
  ENDDO
ELSE IF (.NOT.(OCH_RET_ICE)) THEN
    DO JLW = KEQ-KEQAQ+1, KEQ-KEQAQ/2  ! loop over cloud chem. species
      DO JLG = 1, NSV_CHGSEND-NSV_CHGSBEG+1
        IF ( TRIM(HNAMES(JLW)(4:32)) == TRIM(HNAMES(JLG)) ) THEN
          NINDEXWG(JLW-(KEQ-KEQAQ)) = JLG
          EXIT
        ENDIF
      ENDDO
    ENDDO
END IF
!
 
END SUBROUTINE CH_INIT_ICE

