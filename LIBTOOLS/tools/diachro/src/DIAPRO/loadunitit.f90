!     ######spl
      MODULE  MODI_LOADUNITIT
!     ##############################
!
INTERFACE
!
SUBROUTINE LOADUNITIT(KJ,K)
INTEGER :: KJ,K
END SUBROUTINE LOADUNITIT
!
END INTERFACE
!
END MODULE MODI_LOADUNITIT
!     ######spl
      SUBROUTINE LOADUNITIT(KJ,K)
!     ################################
!
!!****  *LOADUNITIT* - 
!!
!!    PURPOSE
!!    -------
!      
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
!!      Module
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       24/11/95
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_NMGRID
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER          :: KJ,K
!
!*       0.1   Local variables
!              ---------------
! !------------------------------------------------------------------------------
IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
   LSUMVM .OR. LSUTVT .OR.  &
   LDIRWM .OR. LDIRWT .OR. &
   LMLSUMVM .OR. LMLSUTVT .OR. LVTM .OR. LVTT .OR. &
   LULM .OR. LULT)THEN
   CTITGAL=ADJUSTL(CGROUP)
   CTITGAL=ADJUSTL(CTITGAL)
   CUNITGAL(1:LEN(CUNITGAL))=' '
   NMGRID=1
ELSE

  CTITGAL=ADJUSTL(CTITRE(NPROCDIA(KJ,K)))
  CUNITGAL=ADJUSTL(CUNITE(NPROCDIA(KJ,K)))
  CTITGAL=ADJUSTL(CTITGAL)
  IF(CTITGAL(1:LEN_TRIM(CTITGAL)) == 'ZSBIS')THEN
    CTITGAL(1:LEN_TRIM(CTITGAL))=' '
    CTITGAL='ZS'
  ENDIF
  IF(CTITGAL(1:LEN_TRIM(CTITGAL)) == 'ZSMTBIS')THEN
    CTITGAL(1:LEN_TRIM(CTITGAL))=' '
    CTITGAL='ZSMT'
  ENDIF
  CTITGAL=ADJUSTL(CTITGAL)
  CUNITGAL=ADJUSTL(CUNITGAL)
  CUNITGAL(INDEX(CUNITGAL,' '):LEN(CUNITGAL))=' '
  NMGRID=NGRIDIA(NPROCDIA(KJ,K))

ENDIF

IF(NMGRID <1 .OR. NMGRID >7)THEN
  PRINT *,' VALEUR NMGRID ABERRANTE: ',NMGRID, &
        '        FORCEE A        :  1'
  NMGRID=1
ENDIF
RETURN
END SUBROUTINE LOADUNITIT
