!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ######################### 
      MODULE MODI_CH_WRITE_CHEM
!!    ######################### 
!!
!
INTERFACE
SUBROUTINE CH_WRITE_CHEM(PCONC, PAERO, PRHODREF, HFILE)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN) :: PCONC ! PCONC: concentration vector 
REAL, DIMENSION(:), INTENT(IN) :: PAERO ! PCONC: concentration vector 
REAL, DIMENSION(:), INTENT(IN) :: PRHODREF ! air density
CHARACTER(LEN=*), INTENT(IN)     :: HFILE ! HFILE: name of the file for write
END SUBROUTINE CH_WRITE_CHEM
END INTERFACE
END MODULE MODI_CH_WRITE_CHEM
!!
!!    ############################# 
      SUBROUTINE CH_WRITE_CHEM(PCONC, PAERO, PRHODREF, HFILE)
!!    #############################
!!
!!*** *CH_WRITE_CHEM*
!!
!!    PURPOSE
!!    -------
!!    write a set of values to a file
!!
!!**  METHOD
!!    ------
!!    NEQ values will be written to a file,
!!    the format of the file will be as follows:
!!    line 1:     'name1'     value1
!!      ...
!!    line NEQ:   'nameNEQ'     valueNEQ
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
!!    Original 21/04/95
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!    M.Leriche 2015 : masse molaire Black carbon à 12 g/mol
!!    Philippe Wautelet: 10/01/2019: use newunit argument to open files
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9_n, ONLY:      NEQ, CNAMES
USE MODD_CH_AEROSOL
!!
!!    EXPLICIT ARGUMENTS
!!    ----------------
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(INOUT) :: PCONC ! PCONC: concentration vector 
REAL, DIMENSION(:), INTENT(INOUT) :: PAERO ! PCONC: concentration vector 
REAL, DIMENSION(:), INTENT(IN) :: PRHODREF ! air density
CHARACTER(LEN=*), INTENT(IN)     :: HFILE ! HFILE: name of the file for write
!
!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
INTEGER :: ILU ! unit number for IO
INTEGER :: JI, IAER
REAL :: ZMD
REAL, DIMENSION(NSP+NCARB+NSOA) :: ZMI ! aerosol molecular mass in g/mol

!
!*    EXECUTABLE STATEMENTS
!     ---------------------
IF (LORILAM) THEN
IAER=SIZE(PAERO)
ELSE
IAER=0
END IF
!
! open file 
PRINT *, 'CH_WRITE_CHEM: opening file ', HFILE
OPEN(NEWUNIT =  ILU,        &
     FILE    =  HFILE,      &
     FORM    = 'FORMATTED', &
     STATUS  = 'UNKNOWN'    )
!
DO JI = 1, NEQ
  WRITE(UNIT=ILU,FMT='(3A,E20.8)') "'", CNAMES(JI), "' ", PCONC(JI)*1.E9 ! convert ppp to ppb
ENDDO
IF (LORILAM) THEN
!Conversion  ppp to microgram/m3
ZMD    = 28.9644E-3
! Constants initialization
ZMI(:) = 250.
ZMI(JP_AER_SO4)  = 98.
ZMI(JP_AER_NO3)  = 63.
ZMI(JP_AER_NH3)  = 17.
ZMI(JP_AER_H2O)  = 18.
ZMI(JP_AER_BC)   = 12.
ZMI(JP_AER_DST)  = 100.
IF (NSOA == 10) THEN
ZMI(JP_AER_SOA1) = 88. 
ZMI(JP_AER_SOA2) = 180.
ZMI(JP_AER_SOA3) = 1.5374857E+02
ZMI(JP_AER_SOA4) = 1.9586780E+02
ZMI(JP_AER_SOA5) = 195.
ZMI(JP_AER_SOA6) = 195.
ZMI(JP_AER_SOA7) = 165.
ZMI(JP_AER_SOA8) = 195.
ZMI(JP_AER_SOA9) = 270.
ZMI(JP_AER_SOA10) = 210.
END IF

! mineral phase
  PAERO(JP_CH_SO4i) = PAERO(JP_CH_SO4i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SO4))
  PAERO(JP_CH_SO4j) = PAERO(JP_CH_SO4j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SO4))

  PAERO(JP_CH_NO3i) = PAERO(JP_CH_NO3i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_NO3))
  PAERO(JP_CH_NO3j) = PAERO(JP_CH_NO3j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_NO3))

  PAERO(JP_CH_NH3i) = PAERO(JP_CH_NH3i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_NH3))
  PAERO(JP_CH_NH3j) = PAERO(JP_CH_NH3j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_NH3))

! water
  PAERO(JP_CH_H2Oi) = PAERO(JP_CH_H2Oi)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_H2O))
  PAERO(JP_CH_H2Oj) = PAERO(JP_CH_H2Oj)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_H2O))

! primary organic carbon
  PAERO(JP_CH_OCi) = PAERO(JP_CH_OCi)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_OC))
  PAERO(JP_CH_OCj) = PAERO(JP_CH_OCj)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_OC))

! primary black carbon
  PAERO(JP_CH_BCi) = PAERO(JP_CH_BCi)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_BC))
  PAERO(JP_CH_BCj) = PAERO(JP_CH_BCj)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_BC))
!dust
  PAERO(JP_CH_DSTi) = PAERO(JP_CH_DSTi)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_DST))
  PAERO(JP_CH_DSTj) = PAERO(JP_CH_DSTj)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_DST))

IF (NSOA .EQ. 10) THEN
  PAERO(JP_CH_SOA1i) = PAERO(JP_CH_SOA1i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA1))
  PAERO(JP_CH_SOA1j) = PAERO(JP_CH_SOA1j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA1))
  PAERO(JP_CH_SOA2i) = PAERO(JP_CH_SOA2i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA2))
  PAERO(JP_CH_SOA2j) = PAERO(JP_CH_SOA2j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA2))
  PAERO(JP_CH_SOA3i) = PAERO(JP_CH_SOA3i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA3))
  PAERO(JP_CH_SOA3j) = PAERO(JP_CH_SOA3j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA3))
  PAERO(JP_CH_SOA4i) = PAERO(JP_CH_SOA4i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA4))
  PAERO(JP_CH_SOA4j) = PAERO(JP_CH_SOA4j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA4))
  PAERO(JP_CH_SOA5i) = PAERO(JP_CH_SOA5i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA5))
  PAERO(JP_CH_SOA5j) = PAERO(JP_CH_SOA5j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA5))
  PAERO(JP_CH_SOA6i) = PAERO(JP_CH_SOA6i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA6))
  PAERO(JP_CH_SOA6j) = PAERO(JP_CH_SOA6j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA6))
  PAERO(JP_CH_SOA7i) = PAERO(JP_CH_SOA7i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA7))
  PAERO(JP_CH_SOA7j) = PAERO(JP_CH_SOA7j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA7))
  PAERO(JP_CH_SOA8i) = PAERO(JP_CH_SOA8i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA8))
  PAERO(JP_CH_SOA8j) = PAERO(JP_CH_SOA8j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA8))
  PAERO(JP_CH_SOA9i) = PAERO(JP_CH_SOA9i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA9))
  PAERO(JP_CH_SOA9j) = PAERO(JP_CH_SOA9j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA9))
  PAERO(JP_CH_SOA10i) = PAERO(JP_CH_SOA10i)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA10))
  PAERO(JP_CH_SOA10j) = PAERO(JP_CH_SOA10j)/ZMD*1.E6 * &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA10))
END IF
END IF

DO JI = 1, IAER
  WRITE(UNIT=ILU,FMT='(3A,E20.8)') "'", CAERONAMES(JI), "' ", PAERO(JI)
ENDDO

!
! close file
CLOSE(UNIT=ILU)
!
END SUBROUTINE CH_WRITE_CHEM
