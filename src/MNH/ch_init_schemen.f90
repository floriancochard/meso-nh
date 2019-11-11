!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ############################
      MODULE MODI_CH_INIT_SCHEME_n
!     ############################
!
INTERFACE
!
      SUBROUTINE CH_INIT_SCHEME_n(KMI,OUSECHAQ,OUSECHIC, &
                                 OCH_PH,KLUOUT,KVERB     )
!
INTEGER,           INTENT(IN) :: KMI
LOGICAL,           INTENT(IN) :: OUSECHAQ
LOGICAL,           INTENT(IN) :: OUSECHIC
LOGICAL,           INTENT(IN) :: OCH_PH
INTEGER,           INTENT(IN) :: KLUOUT
INTEGER,           INTENT(IN) :: KVERB
!
!
END SUBROUTINE CH_INIT_SCHEME_n
!
END INTERFACE
!
END MODULE MODI_CH_INIT_SCHEME_n
!
!
!      ####################################################
       SUBROUTINE CH_INIT_SCHEME_n(KMI,OUSECHAQ,OUSECHIC, &
                                  OCH_PH,KLUOUT,KVERB     )
!      ####################################################
!!
!!***  *CH_INIT_SCHEME*
!!
!!    PURPOSE
!!    -------
!      Initialize module MODD_CH_M9_n variables (now MesoNH variables) from 
!     internal scheme constants defined in module MODD_CH_M9_SCHEME by 
!!    BASIC.f90 (scheme dependant).
!!
!!**  METHOD
!!    ------
!!    Names of variables are identical in both MODD_CH_M9_n and MODD_CH_M9_SCHEME
!!    modules to minimize sources modifications and keep users habits.
!!    So we rename variables/constant from MODD_CH_M9_SCHEME by prefixing them 
!!    with I_ (I means scheme Internal).
!!
!!    REFERENCES
!!    ----------
!!    MesoNH-chemistry book 3
!!
!!    AUTHOR
!!    ------
!!    D. Gazen
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 19/10/03
!!    12/04/07 (M. Leriche) add aqueous chemistry
!!    15/07/10 (M. Leriche) add ice phase chemistry
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    None
!----------------------------------------------------------------------------
!
USE MODD_CH_M9_n !! the mesonh interface for chemical variables
USE MODD_CH_M9_SCHEME, ONLY :  TACCS, &
     & I_CNAMES=>CNAMES,               &          
     & I_CREACS=>CREACS,               &
     & I_CFULLREACS=>CFULLREACS
!
USE MODI_CH_INIT_CCS
!
IMPLICIT NONE
!
!*       0.   DECLARATIONS
!        -----------------
!
INTEGER,           INTENT(IN) :: KMI
LOGICAL,           INTENT(IN) :: OUSECHAQ
LOGICAL,           INTENT(IN) :: OUSECHIC
LOGICAL,           INTENT(IN) :: OCH_PH
INTEGER,           INTENT(IN) :: KLUOUT
INTEGER,           INTENT(IN) :: KVERB
!
INTEGER :: JI,JJ   ! loop counters
LOGICAL :: GCO2
!----------------------------------------------------------------------------
!
!*       1.   INITIALISATION
!        -------------------

! Initialize BASIC for model KMI
!
CALL CH_INIT_CCS(KMI,OUSECHAQ,OCH_PH,KLUOUT,KVERB)
!
!Left member belongs to MODD_CH_M9_n module variables
!
NEQ           = TACCS(KMI)%NEQ
NEQAQ         = TACCS(KMI)%NEQAQ
NREAC         = TACCS(KMI)%NREAC
NMETEOVARS    = TACCS(KMI)%NMETEOVARS
NNONZEROTERMS = TACCS(KMI)%NNONZEROTERMS
! 
CNAMES=>I_CNAMES(1:NEQ)
CREACS=>I_CREACS(1:NREAC)
CFULLREACS=>I_CFULLREACS(1:NREAC)
IF (OUSECHIC) THEN
!add ice phase chemistry 
  GCO2 = .false.
  DO JJ = NEQ-NEQAQ+1, NEQ-NEQAQ/2
    IF (TRIM(CNAMES(JJ)(4:32)) == 'CO2') THEN
      GCO2 = .true.
      EXIT
    ENDIF
  ENDDO 
  IF (GCO2) THEN
    ALLOCATE(CICNAMES(NEQAQ/2-1))
    GCO2 = .false.
    DO JI = 1, NEQAQ/2 -1
      JJ = NEQ-NEQAQ+JI
      IF (TRIM(CNAMES(JJ)(4:32)) == 'CO2') GCO2 = .true.
      IF (GCO2) THEN
        CICNAMES(JI) = 'IC_'//TRIM(CNAMES(JJ+1)(4:32))
      ELSE
        CICNAMES(JI) = 'IC_'//TRIM(CNAMES(JJ)(4:32))
      ENDIF
    ENDDO
  ELSE
    ALLOCATE(CICNAMES(NEQAQ/2)) 
    DO JI = 1, NEQAQ/2
      JJ = NEQ-NEQAQ+JI
      CICNAMES(JI) = 'IC_'//TRIM(CNAMES(JJ)(4:32))
    ENDDO
  ENDIF
ENDIF
!
WRITE(KLUOUT,*) 'CH_INIT_SCHEME_n: KMI, NEQ, NEQAQ, NREAC, NMETEOVARS, NNONZEROTERMS = ',&
     & KMI, NEQ, NEQAQ, NREAC, NMETEOVARS, NNONZEROTERMS
 
END SUBROUTINE CH_INIT_SCHEME_n

