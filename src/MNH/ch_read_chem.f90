!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ######################## 
      MODULE MODI_CH_READ_CHEM
!!    ######################## 
!!
!
INTERFACE
SUBROUTINE CH_READ_CHEM(PCONC, PAERO, PRHODREF, HFILE)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(INOUT) :: PCONC ! gas concentration vector to be read
REAL, DIMENSION(:), INTENT(INOUT) :: PAERO ! aerosol concentration vector to be read
REAL, DIMENSION(:), INTENT(IN) :: PRHODREF ! air density
CHARACTER(LEN=*), INTENT(IN)      :: HFILE ! name of the file to be read from
END SUBROUTINE CH_READ_CHEM
END INTERFACE
!!
END MODULE MODI_CH_READ_CHEM
!!
!!    ##################################### 
      SUBROUTINE CH_READ_CHEM(PCONC, PAERO, PRHODREF, HFILE)
!!    #####################################
!!
!!*** *CH_READ_CHEM*
!!
!!    PURPOSE
!!    -------
!!    read a set of values from a file
!!
!!    METHOD
!!    ------
!!    NEQ values will be read from a file and their names will be checked,
!!    the format of the file has to be as follows:
!!    line 1:     'name1'     value1
!!      ...
!!    line NEQ:   'nameNEQ'     valueNEQ
!!
!!    if the filename is CHCONTROL1.nam, the file pointer will be placed
!!    after the keyword INITCHEM and an additional comment line will be read
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
!!    18/02/99 (K. Suhre) allow reading from namelist file
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!    M.Leriche 2015 : masse molaire Black carbon à 12 g/mol
!!    Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!    Philippe Wautelet: 10/01/2019: use newunit argument to open files
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_OPEN_INPUT
USE MODI_CH_READ_VECTOR
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_FM,    ONLY: IO_FILE_CLOSE_ll
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D, ONLY: NVERB
USE MODD_CH_M9_n, ONLY:      NEQ, CNAMES
USE MODD_CH_AEROSOL
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
!REAL, DIMENSION(NEQ), INTENT(OUT) :: PCONC ! concentration vector to be read
REAL, DIMENSION(:), INTENT(INOUT) :: PCONC ! concentration vector to be read
REAL, DIMENSION(:), INTENT(INOUT) :: PAERO ! aerosol concentration vector to be read
REAL, DIMENSION(:), INTENT(IN) :: PRHODREF ! air density
CHARACTER(LEN=*), INTENT(IN)      :: HFILE ! name of the file to be read from
!
!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
CHARACTER(LEN=32) :: YVARNAME
CHARACTER(LEN=80) :: YINPUT
INTEGER :: ILU ! unit number for IO
INTEGER :: JI, JJ, IIN
REAL :: ZMD
REAL, DIMENSION(NSP+NCARB+NSOA) :: ZMI ! aerosol molecular mass in g/mol
TYPE(TFILEDATA),POINTER :: TZFILE
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
TZFILE => NULL()
!
! if the namelist file is the input file, position file pointer after keyword
!
IF (HFILE(1:14) .EQ. "CHCONTROL1.nam") THEN
!
  PRINT *, "reading initial data from CHCONTROL1.nam"
!
! open the namelist file for input
!
  CALL CH_OPEN_INPUT("CHCONTROL1.nam", "INITCHEM", TZFILE, 6, NVERB)
  IIN = TZFILE%NLU
!
  CALL CH_READ_VECTOR(NEQ, CNAMES, PCONC, 0.0, IIN, 6, NVERB)
  IF (LORILAM) CALL CH_READ_VECTOR(SIZE(PAERO,1), CAERONAMES, PAERO, 0.0, IIN, 6, NVERB)
!
  CALL IO_FILE_CLOSE_ll(TZFILE)
!
ELSE
!
! open file 
!
  PRINT *, 'CH_READ_CHEM: opening file ', HFILE
  OPEN(NEWUNIT = ILU,         &
       FILE    = HFILE,       &
       FORM    = 'FORMATTED', &
       STATUS  = 'OLD')
!
! read line by line and check variable names
!
  outer_loop : DO JI = 1, NEQ
    READ(UNIT=ILU,FMT=*,END=999) YVARNAME, PCONC(JI)
    check_loop : DO JJ = 1, 32
      IF (YVARNAME(JJ:JJ).NE.CNAMES(JI)(JJ:JJ)) THEN
        PRINT *, 'CH_READ_CHEM: Error: variable names do not match:'
        PRINT *, 'CNAMES = >>>', CNAMES(JI), '<<<'
        PRINT *, 'read   = >>>', YVARNAME, '<<<'
!callabortstop
        CALL ABORT
        STOP 'Program stopped by CH_READ_CHEM'
      ENDIF
    ENDDO check_loop
  ENDDO outer_loop

!Conversion ppb to ppp
 PCONC(:) =  PCONC(:) * 1E-9
IF (LORILAM) THEN
  outer_loop2 : DO JI = 1, SIZE(PAERO,1)
    READ(UNIT=ILU,FMT=*,END=997) YVARNAME, PAERO(JI)
    check_loop2 : DO JJ = 1, 32
      IF (YVARNAME(JJ:JJ).NE.CAERONAMES(JI)(JJ:JJ)) THEN
        PRINT *, 'CH_READ_CHEM: Error: variable names do not match:'
        PRINT *, 'CAERONAMES = >>>', CAERONAMES(JI), '<<<'
        PRINT *, 'read   = >>>', YVARNAME, '<<<'
!callabortstop
        CALL ABORT
        STOP 'Program stopped by CH_READ_CHEM'
      ENDIF
    ENDDO check_loop2
  ENDDO outer_loop2
!Conversion  microgram/m3 to ppp
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
  PAERO(JP_CH_SO4i) = PAERO(JP_CH_SO4i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SO4))
  PAERO(JP_CH_SO4j) = PAERO(JP_CH_SO4j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SO4))

  PAERO(JP_CH_NO3i) = PAERO(JP_CH_NO3i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_NO3))
  PAERO(JP_CH_NO3j) = PAERO(JP_CH_NO3j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_NO3))

  PAERO(JP_CH_NH3i) = PAERO(JP_CH_NH3i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_NH3))
  PAERO(JP_CH_NH3j) = PAERO(JP_CH_NH3j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_NH3))

! water
  PAERO(JP_CH_H2Oi) = PAERO(JP_CH_H2Oi)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_H2O))
  PAERO(JP_CH_H2Oj) = PAERO(JP_CH_H2Oj)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_H2O))

! primary organic carbon
  PAERO(JP_CH_OCi) = PAERO(JP_CH_OCi)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_OC))
  PAERO(JP_CH_OCj) = PAERO(JP_CH_OCj)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_OC))

! primary black carbon
  PAERO(JP_CH_BCi) = PAERO(JP_CH_BCi)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_BC))
  PAERO(JP_CH_BCj) = PAERO(JP_CH_BCj)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_BC))

!dust
  PAERO(JP_CH_DSTi) = PAERO(JP_CH_DSTi)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_DST))
  PAERO(JP_CH_DSTj) = PAERO(JP_CH_DSTj)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_DST))

IF (NSOA .EQ. 10) THEN
  PAERO(JP_CH_SOA1i) = PAERO(JP_CH_SOA1i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA1))
  PAERO(JP_CH_SOA1j) = PAERO(JP_CH_SOA1j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA1))
  PAERO(JP_CH_SOA2i) = PAERO(JP_CH_SOA2i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA2))
  PAERO(JP_CH_SOA2j) = PAERO(JP_CH_SOA2j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA2))
  PAERO(JP_CH_SOA3i) = PAERO(JP_CH_SOA3i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA3))
  PAERO(JP_CH_SOA3j) = PAERO(JP_CH_SOA3j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA3))
  PAERO(JP_CH_SOA4i) = PAERO(JP_CH_SOA4i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA4))
  PAERO(JP_CH_SOA4j) = PAERO(JP_CH_SOA4j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA4))
  PAERO(JP_CH_SOA5i) = PAERO(JP_CH_SOA5i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA5))
  PAERO(JP_CH_SOA5j) = PAERO(JP_CH_SOA5j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA5))
  PAERO(JP_CH_SOA6i) = PAERO(JP_CH_SOA6i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA6))
  PAERO(JP_CH_SOA6j) = PAERO(JP_CH_SOA6j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA6))
  PAERO(JP_CH_SOA7i) = PAERO(JP_CH_SOA7i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA7))
  PAERO(JP_CH_SOA7j) = PAERO(JP_CH_SOA7j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA7))
  PAERO(JP_CH_SOA8i) = PAERO(JP_CH_SOA8i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA8))
  PAERO(JP_CH_SOA8j) = PAERO(JP_CH_SOA8j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA8))
  PAERO(JP_CH_SOA9i) = PAERO(JP_CH_SOA9i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA9))
  PAERO(JP_CH_SOA9j) = PAERO(JP_CH_SOA9j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA9))
  PAERO(JP_CH_SOA10i) = PAERO(JP_CH_SOA10i)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA10))
  PAERO(JP_CH_SOA10j) = PAERO(JP_CH_SOA10j)*ZMD*1.E-6 / &
                                  (PRHODREF(1) * ZMI(JP_AER_SOA10))
END IF

END IF
!
! close file
  CLOSE(UNIT=ILU)
!

END IF
!
RETURN
!
!-----------------------------------------------------------------------------
!
999  PRINT *, 'CH_READ_CHEM: Error: not enough variables defined in file', &
               HFILE
PRINT *, 'number of gas lines in that file should be ', NEQ, &
         ', but is ', JI-1
!callabortstop
CALL ABORT
STOP 'Program stopped by CH_READ_CHEM'
!
998 STOP "CH_READ_CHEM: ERROR - keyword INITCHEM not found in CHCONTROL1.nam"  
!
997  PRINT *, 'CH_READ_CHEM: Error: not enough variables defined in file', &
               HFILE
PRINT *, 'number of aerosols lines in that file should be ', SIZE(PAERO,1), &
         ', but is ', JI-1
!callabortstop
CALL ABORT
STOP 'Program stopped by CH_READ_CHEM'
!
END SUBROUTINE CH_READ_CHEM
