!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ######################### 
      MODULE MODI_CH_OPEN_INPUT
!!    ######################### 
!!
INTERFACE
!!
SUBROUTINE CH_OPEN_INPUT(HCHEM_INPUT_FILE,HKEYWORD,TPFILE,KLUOUT,KVERB)
!
USE MODD_IO_ll,            ONLY: TFILEDATA
!
IMPLICIT NONE
!
CHARACTER(LEN=*),        INTENT(IN)  :: HCHEM_INPUT_FILE ! general purpose input file
CHARACTER(LEN=*),        INTENT(IN)  :: HKEYWORD         ! keyword for positioning
TYPE(TFILEDATA),POINTER, INTENT(OUT) :: TPFILE           ! Newly opened file
INTEGER,                 INTENT(IN)  :: KLUOUT           ! output listing channel
INTEGER,                 INTENT(IN)  :: KVERB            ! verbosity level
END SUBROUTINE CH_OPEN_INPUT
!!
END INTERFACE
!!
END MODULE MODI_CH_OPEN_INPUT
!!
!!    ######################################################################### 
SUBROUTINE CH_OPEN_INPUT(HCHEM_INPUT_FILE,HKEYWORD,TPFILE,KLUOUT,KVERB)
!!    #########################################################################
!!
!!*** *CH_OPEN_INPUT*
!!
!!    PURPOSE
!!    -------
!       Open the general purpose chemical input file and position the
!     file pointer after the indicated keyword.
!!
!!**  METHOD
!!    ------
!!      An unused input channel is selected using OPEN_ll.
!!    The file HCHEM_INPUT_FILE will be rewinded
!!    at each call and data will be read in using (A8)-format until the
!!    given keyword is found. The following comment line will then
!!    be read and printed and the input channel number will be returned.
!!    After reading the needed data, the user must assure that the file
!!    will be closed and that the unit will be freed using CLOSE_ll.
!!
!!    REFERENCE
!!    ---------
!!    MesoNH book 2
!!
!!    AUTHOR
!!    ------
!!    K. Suhre   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 03/11/95
!!    05/08/96 (K. Suhre) restructured
!!    11/08/98 (N. Asencio) add parallel code
!!    Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_IO_ll,            ONLY: TFILEDATA
!
USE MODE_FM,               ONLY: IO_FILE_OPEN_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=*),        INTENT(IN)  :: HCHEM_INPUT_FILE ! general purpose input file
CHARACTER(LEN=*),        INTENT(IN)  :: HKEYWORD         ! keyword for positioning
TYPE(TFILEDATA),POINTER, INTENT(OUT) :: TPFILE           ! Newly opened file
INTEGER,                 INTENT(IN)  :: KLUOUT           ! output listing channel
INTEGER,                 INTENT(IN)  :: KVERB            ! verbosity level
!
!*      0.2    declarations of local variables
!
CHARACTER(LEN=79) :: YIN ! character string for line-by-line read
INTEGER :: ILU
INTEGER :: IRESP         ! return code from OPEN_ll
!
!-------------------------------------------------------------------------------
!
TPFILE => NULL()
!
!*       1.    FIND A FREE I/O CHANNEL
!              -----------------------
!
IF (KVERB >= 5) WRITE(KLUOUT,*) "CH_OPEN_INPUT: opening file ", HCHEM_INPUT_FILE
CALL IO_FILE_ADD2LIST(TPFILE,HCHEM_INPUT_FILE,'CHEMINPUT','READ',OOLD=.TRUE.)
CALL IO_FILE_OPEN_ll(TPFILE,KRESP=IRESP)
ILU = TPFILE%NLU
!
IF (IRESP /= 0) THEN
  WRITE(KLUOUT,*) "CH_OPEN_INPUT ERROR: unable to open file", HCHEM_INPUT_FILE
  WRITE(KLUOUT,*) "                     IO_FILE_OPEN_ll return code is: ", IRESP
  WRITE(KLUOUT,*) "                     the program will be stopped now"
  ! callabortstop
  CALL ABORT
  STOP
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    READ INPUT DATA UNTIL KEYWORD IS FOUND
!              --------------------------------------
!
! open the input file
REWIND(ILU)
WRITE(KLUOUT,*) "CH_OPEN_INPUT: opened channel ", ILU, &
                " for file ", HCHEM_INPUT_FILE
!
! read general comment line and print it
READ(ILU,"(A)") YIN
IF (KVERB >= 5) WRITE(KLUOUT, *) YIN
!
search_key : DO
  READ(ILU,"(A8)", END=100) YIN
  IF (HKEYWORD(1:8) .EQ. YIN(1:8)) EXIT search_key
ENDDO search_key
!
! read specific comment line and print it
IF (KVERB >= 5) WRITE(KLUOUT,*) &
   "Keyword ", HKEYWORD(1:8), " has been found, the specific comment is:"
READ(ILU,"(A)") YIN
IF (KVERB >= 5) WRITE(KLUOUT, *) YIN
!
RETURN
!
!-------------------------------------------------------------------------------
!
!*       3.    IF THE KEYWORD HAS NOT BEEN FOUND, STOP
!              ---------------------------------------
!
100 CONTINUE
WRITE(KLUOUT,*) "CH_OPEN_INPUT-Error: Keyword ", HKEYWORD(1:8), " not found."
! callabortstop
CALL ABORT
STOP "Program stopped"
!
END SUBROUTINE CH_OPEN_INPUT
