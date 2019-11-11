!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ##############################
      SUBROUTINE CH_INIT_DIAGNOSTICS
!!    ##############################
!!
!!****  *CH_INIT_DIAGNOSTICS*
!!
!!    PURPOSE
!!    -------
!!    prepare diagnostics
!!
!!    METHOD
!!    ------
!!    the diagnostics file CDIAGFILE (default: "BOX.DIAG") will be opened
!!    and the headder will be written, containing general information
!!    on the saved variables, actually the following variables will be saved:
!!    CONC, reaction rates, PROD, LOSS, TERMS;
!!    these variables are  calculated by CH_PRODLOSS, CH_TERMS, and CH_GET_RATES
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 14/05/95
!!    27/07/96 (K. Suhre) restructured
!!    01/12/03  (Gazen)   change Chemical scheme interface
!!    Philippe Wautelet: 10/01/2019: use newunit argument to open files
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D, ONLY: CDIAGFILE, CRUNID, CDIAGFORMAT, &
                           NDIAGIO, XTSIMUL
USE MODD_CH_M9_n,    ONLY: NEQ, NREAC, NNONZEROTERMS, CNAMES, CREACS

USE MODE_MODELN_HANDLER
!!
!!
!!
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!        -----------------------------
!
!*       0.2 Declaration of local variables
!        ----------------------------------
!
CHARACTER*8  :: YDATE  ! for retrieval of date and time
CHARACTER*10 :: YTIME  ! dito
!
INTEGER :: JI
!INTEGER, DIMENSION(2, NNONZEROTERMS) :: IINDEX ! indices of non-zero terms
INTEGER, DIMENSION(:,:), ALLOCATABLE  :: IINDEX ! indices of non-zero terms
INTEGER :: IMI   ! model index
!
IMI = GET_CURRENT_MODEL_INDEX()
! Allocate IINDEX
ALLOCATE(IINDEX(2,NNONZEROTERMS))
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
! open the file
OPEN(NEWUNIT =  NDIAGIO,    &
     FILE    =  CDIAGFILE,  &
     FORM    = "FORMATTED", &
     STATUS  = "UNKNOWN"    )
PRINT *, "CH_INIT_DIAGNOSTICS: opening unit ", NDIAGIO, " for file ", CDIAGFILE
!
! write the headder
CALL DATE_AND_TIME(YDATE, YTIME)
WRITE(NDIAGIO,'(4A)')'MODEL0D DIAGNOSTICS VERSION 1.1 AT ',YDATE,":",YTIME
WRITE(NDIAGIO,'(A)')  CRUNID
!
! number of variables and their names
WRITE(NDIAGIO,'(4I10)') 1+NEQ+NREAC+2*NEQ+NNONZEROTERMS, &
			NEQ, NREAC, NNONZEROTERMS
!
! diagnostics time
WRITE(NDIAGIO, '(A)') "XTSIMUL"
!
! the concentrations
DO JI = 1, NEQ
  WRITE(NDIAGIO, '(A)') CNAMES(JI)
ENDDO
!
! the reaction rates
DO JI = 1, NREAC
  WRITE(NDIAGIO, '(A)') CREACS(JI)
ENDDO
!
! the production terms
DO JI = 1, NEQ
  WRITE(NDIAGIO, '(A)') CNAMES(JI) // "-PROD"
ENDDO
!
! the destruction terms terms
DO JI = 1, NEQ
  WRITE(NDIAGIO, '(A)') CNAMES(JI) // "-LOSS"
ENDDO
!
! get the indices of the explicit contributions
! (only non-zero terms will be stored)
CALL CH_NONZEROTERMS(IMI,IINDEX, NNONZEROTERMS)
DO JI = 1, NNONZEROTERMS
  WRITE(NDIAGIO, '(A)') CNAMES(IINDEX(1,JI)) // "-" // CREACS(IINDEX(2,JI))
ENDDO
DO JI = 1, NNONZEROTERMS
  WRITE(NDIAGIO, '(2I5)') IINDEX(1,JI), IINDEX(2,JI)
ENDDO
!
! write format of following data
WRITE(NDIAGIO,'(A)') CDIAGFORMAT
!
DEALLOCATE(IINDEX)
!
END SUBROUTINE CH_INIT_DIAGNOSTICS
