!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ##########################
      MODULE MODI_CH_INIT_OUTPUT
!!    ##########################
!!
!
INTERFACE
SUBROUTINE CH_INIT_OUTPUT(TPM)
USE MODD_CH_M9_n, ONLY: METEOTRANSTYPE
IMPLICIT NONE
TYPE(METEOTRANSTYPE), INTENT(IN) :: TPM
END SUBROUTINE CH_INIT_OUTPUT
END INTERFACE
END MODULE MODI_CH_INIT_OUTPUT
!!
!!    ############################## 
      SUBROUTINE CH_INIT_OUTPUT(TPM)
!!    ##############################
!!
!!*** *CH_INIT_OUTPUT*
!!
!!    PURPOSE
!!    -------
!!    prepare regular output of results
!!
!!**  METHOD
!!    ------
!!    the result file CRESULTFILE (default: "BOX.RESULT") will be opened
!!    and the headder will be written, containing general information
!!    on the saved variables, presently the following variables will be saved:
!!    - concentration of all prognostic variables
!!    - the reaction constants and photolysis rates
!!    - the meteo variables
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 25/04/95
!!    27/07/96 (K. Suhre) restructured
!!    01/12/03 (D. Gazen) change Chemical scheme interface
!!    Philippe Wautelet: 10/01/2019: use newunit argument to open files
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D,  ONLY: CRESULTFILE, CRUNID, CRESULTFORMAT, &
                            NRESULTIO, XTSIMUL
USE MODD_CH_M9_n,     ONLY: NEQ, NREAC, NMETEOVARS, CNAMES, CREACS, METEOTRANSTYPE
USE MODD_CH_AEROSOL
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
TYPE(METEOTRANSTYPE), INTENT(IN) :: TPM  ! the meteo variables
!
!*       0.2  declaration of local variables
!
CHARACTER*8            :: YDATE  ! for retrieval of date and time
CHARACTER*10           :: YTIME  ! dito
INTEGER                :: JI     ! loop control
INTEGER                :: NAERO
!
!------------------------------------------------------------------------------
!
!*       1.   OPEN OUTPUT FILE
!        ----------------------
!
OPEN(NEWUNIT = NRESULTIO,   &
     FILE    = CRESULTFILE, &
     FORM    = "FORMATTED", &
     STATUS  = "UNKNOWN"    )
PRINT *, "CH_INIT_OUTPUT: opening unit ", NRESULTIO, " for file ", CRESULTFILE
!
!*       2.   WRITE HEADDER
!        ------------------
!
!*       2.1  write two lines of comment (date, time, runid)
!
IF (LORILAM) THEN
NAERO=SIZE(CAERONAMES)
ELSE
NAERO=0
END IF

CALL DATE_AND_TIME(YDATE, YTIME)
WRITE(NRESULTIO,'(4A)')'MODEL0D RESULTS VERSION 1.1 AT ',YDATE,":",YTIME
WRITE(NRESULTIO,'(A)')  CRUNID
!
!*       2.2  write number of variables and their names
!
WRITE(NRESULTIO,'(4I6)') 1+NEQ+NAERO+NREAC+NMETEOVARS, NEQ, NREAC, NMETEOVARS
!
!*       2.3  write simulation time
!
WRITE(NRESULTIO, '(A)') "XTSIMUL"
!
!*       2.4  write chemical concentrations
!
DO JI = 1, NEQ
  WRITE(NRESULTIO, '(A)') CNAMES(JI)
ENDDO
DO JI = 1, NAERO
  WRITE(NRESULTIO, '(A)') CAERONAMES(JI)
ENDDO
!
!*       2.5  write reaction and photolysis rates
!
DO JI = 1, NREAC
  WRITE(NRESULTIO, '(A)') CREACS(JI)
ENDDO
!
!*       2.6  write meteo variables
!
DO JI = 1, NMETEOVARS
  WRITE(NRESULTIO, '(A)') TPM%CMETEOVAR(JI)
ENDDO
!
!*       2.7  write data output format
!
WRITE(NRESULTIO,'(A)') CRESULTFORMAT
!
END SUBROUTINE CH_INIT_OUTPUT
