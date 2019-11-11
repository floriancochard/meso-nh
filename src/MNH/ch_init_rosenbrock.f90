!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##############################
      MODULE MODI_CH_INIT_ROSENBROCK
!     ##############################
!
INTERFACE
!
      SUBROUTINE CH_INIT_ROSENBROCK(KMI,KLUOUT)
!
INTEGER, INTENT(IN) :: KMI
INTEGER, INTENT(IN) :: KLUOUT
!
!
END SUBROUTINE CH_INIT_ROSENBROCK
!
END INTERFACE
!
END MODULE MODI_CH_INIT_ROSENBROCK
!
!
!      #####################################
       SUBROUTINE CH_INIT_ROSENBROCK(KMI,KLUOUT)
!      #####################################
!!
!!***  *CH_INIT_ROSENBROCK*
!!
!!    PURPOSE
!!    -------
!      initialize the sparse index vector of module MODD_CH_ROSENBROCK
!!
!!**  METHOD
!!    ------
!!    
!!   
!!    REFERENCES
!!    ----------
!!    MesoNH-chemistry book 3 and kpp documentation of D. Sandu
!!
!!    AUTHOR
!!    ------
!!    J.-P. Pinty    *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 05/06/07
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    None
!!
!-----------------------------------------------------------------------------
!
!
USE MODI_CH_SPARSE
!
USE MODD_CH_M9_n, ONLY: NEQ, NEQAQ, NNONZEROTERMS
USE MODD_CH_ROSENBROCK_n
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: KMI
INTEGER, INTENT(IN) :: KLUOUT
!
!
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISPARSE
INTEGER :: ISPARSEDIM
INTEGER :: JL, JLL ! Loop indexes
!
LOGICAL :: GCHECK_DIAG ! Used to check the diagonal element of each row
!
INTEGER :: INAQ    ! Running index
INTEGER, DIMENSION(:), ALLOCATABLE :: NSPARSE_IROW_NAQ_WORK, &
                                      NSPARSE_ICOL_NAQ_WORK
!
!----------------------------------------------------------------------------
!
!*       1.  SPARSE JACOBIAN 
!        -------------------
!
ALLOCATE(ISPARSE(2,NNONZEROTERMS))
ISPARSEDIM = NNONZEROTERMS
CALL CH_SPARSE(KMI,ISPARSE, ISPARSEDIM)
NSPARSEDIM = ISPARSEDIM
!
ALLOCATE(NSPARSE_IROW(NSPARSEDIM))
ALLOCATE(NSPARSE_ICOL(NSPARSEDIM))
NSPARSE_IROW(1:NSPARSEDIM) = ISPARSE(1,1:NSPARSEDIM)
NSPARSE_ICOL(1:NSPARSEDIM) = ISPARSE(2,1:NSPARSEDIM)
DEALLOCATE(ISPARSE)
!
ALLOCATE(NSPARSE_CROW(NEQ+1))
ALLOCATE(NSPARSE_DIAG(NEQ+1))
!
!
          DO JL = 1, NEQ
JLL_Loop1:  DO JLL = 1, NSPARSEDIM
              IF( NSPARSE_IROW(JLL)==JL ) THEN
                NSPARSE_CROW(JL) = JLL
                EXIT JLL_Loop1
              END IF
            END DO JLL_Loop1 
            GCHECK_DIAG = .FALSE.
JLL_Loop2:  DO JLL = 1, NSPARSEDIM
              IF( NSPARSE_IROW(JLL)==JL .AND. NSPARSE_ICOL(JLL)==JL ) THEN
                NSPARSE_DIAG(JL) = JLL
                GCHECK_DIAG = .TRUE.
                EXIT JLL_Loop2
              END IF
            END DO JLL_Loop2
            IF( .NOT.GCHECK_DIAG ) THEN
              WRITE(KLUOUT,*)"THE PROGRAM STOPPED IN CH_INIT_ROSENBROCK AS NO"
              WRITE(KLUOUT,*)"DIAGONAL ELEMENT IS FOUND FOR CHEMICAL COMPOUND"
              WRITE(KLUOUT,*)"NUMBER: ",JLL," IN THE JACOBIAN MATRIX !!!"
              WRITE(KLUOUT,*)"PLEASE MODIFY AND REPROCESS THE CHEMICAL SYSTEM"
              STOP
            ENDIF
          END DO
!
NSPARSE_CROW(NEQ+1) = NSPARSEDIM + 1
NSPARSE_DIAG(NEQ+1) = NSPARSEDIM + 1
!
!
WRITE(KLUOUT,*) 'CH_INIT_ROSENBROCK done : NEQ = ', NEQ
WRITE(KLUOUT,*) 'CH_INIT_ROSENBROCK done : NSPARSEDIM = ', NSPARSEDIM
!
! Now explore the subset of NonAQueous species
!
IF( NEQAQ>0 ) THEN
  NEQ_NAQ = NEQ - NEQAQ
!
  ALLOCATE(NSPARSE_IROW_NAQ_WORK(NSPARSEDIM))
  ALLOCATE(NSPARSE_ICOL_NAQ_WORK(NSPARSEDIM))
  INAQ = 0
  DO JLL = 1, NSPARSEDIM
    IF( NSPARSE_IROW(JLL)<=NEQ_NAQ .AND. NSPARSE_ICOL(JLL)<=NEQ_NAQ ) THEN 
      INAQ = INAQ + 1
      NSPARSE_IROW_NAQ_WORK(INAQ) = NSPARSE_IROW(JLL)
      NSPARSE_ICOL_NAQ_WORK(INAQ) = NSPARSE_ICOL(JLL)
    END IF
  END DO
  NSPARSEDIM_NAQ = INAQ
  ALLOCATE(NSPARSE_IROW_NAQ(NSPARSEDIM_NAQ))
  ALLOCATE(NSPARSE_ICOL_NAQ(NSPARSEDIM_NAQ))
  NSPARSE_IROW_NAQ(1:NSPARSEDIM_NAQ) = NSPARSE_IROW_NAQ_WORK(1:NSPARSEDIM_NAQ)
  NSPARSE_ICOL_NAQ(1:NSPARSEDIM_NAQ) = NSPARSE_ICOL_NAQ_WORK(1:NSPARSEDIM_NAQ)
  DEALLOCATE(NSPARSE_IROW_NAQ_WORK)
  DEALLOCATE(NSPARSE_ICOL_NAQ_WORK)
!
  ALLOCATE(NSPARSE_CROW_NAQ(NEQ_NAQ+1))
  ALLOCATE(NSPARSE_DIAG_NAQ(NEQ_NAQ+1))
!
!
          DO JL = 1, NEQ_NAQ
JLL_Loop3:  DO JLL = 1, NSPARSEDIM_NAQ
              IF( NSPARSE_IROW_NAQ(JLL)==JL ) THEN
                NSPARSE_CROW_NAQ(JL) = JLL
                EXIT JLL_Loop3
              END IF
            END DO JLL_Loop3 
JLL_Loop4:  DO JLL = 1, NSPARSEDIM_NAQ
              IF( NSPARSE_IROW_NAQ(JLL)==JL.AND.NSPARSE_ICOL_NAQ(JLL)==JL ) THEN
                NSPARSE_DIAG_NAQ(JL) = JLL
                EXIT JLL_Loop4
              END IF
            END DO JLL_Loop4
          END DO
!
  NSPARSE_CROW_NAQ(NEQ_NAQ+1) = NSPARSEDIM_NAQ + 1
  NSPARSE_DIAG_NAQ(NEQ_NAQ+1) = NSPARSEDIM_NAQ + 1
!
!
  WRITE(KLUOUT,*) 'CH_INIT_ROSENBROCK done : NEQ_NAQ = ', NEQ_NAQ
  WRITE(KLUOUT,*) 'CH_INIT_ROSENBROCK done : NSPARSEDIM_NAQ = ', NSPARSEDIM_NAQ
!
END IF
!
END SUBROUTINE CH_INIT_ROSENBROCK
