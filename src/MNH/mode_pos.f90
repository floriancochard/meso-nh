!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 mode 2006/05/18 13:07:25
!-----------------------------------------------------------------
!!    ###############
      MODULE MODE_POS
!!    ###############
!!
INTERFACE POS
!!
MODULE PROCEDURE POSNAM
MODULE PROCEDURE POSKEY
!!
END INTERFACE
!!
!!
CONTAINS
!!
!!    ##############################################
      SUBROUTINE POSNAM(KULNAM,HDNAML,OFOUND,KLUOUT)
!!    ##############################################
!!
!!*** *POSNAM*
!!
!!    PURPOSE
!!    -------
!     To position namelist file at correct place for reading
!     namelist CDNAML.
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENT
!!    -----------------
!!
!!    REFERENCE
!!    ----------
!!       ECMWF Research Department documentation of the IFS
!!
!!    AUTHOR
!!    -------
!!       Mats Hamrud *ECMWF*
!!
!!    MODIFICATIONS
!!    --------------
!!       Original : 22/06/93
!!       I. Mallet  15/10/01     adaptation to MesoNH (F90 norm)
!------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!*       0.1   Declarations of arguments
!
INTEGER,          INTENT(IN) :: KULNAM
CHARACTER(LEN=*), INTENT(IN) :: HDNAML
LOGICAL,          INTENT(OUT):: OFOUND
INTEGER, OPTIONAL,INTENT(IN) :: KLUOUT
!
!*       0.2   Declarations of local variables
!
CHARACTER(LEN=120) :: YLINE
CHARACTER(LEN=1)   :: YLTEST
INTEGER            :: ILEN,ILEY,IND1,IRET
INTEGER            :: J,JA
!
CHARACTER(LEN=1),DIMENSION(26) :: YLO=(/'a','b','c','d','e','f','g','h', &
     'i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'/)
CHARACTER(LEN=1),DIMENSION(26) :: YUP=(/'A','B','C','D','E','F','G','H', &
     'I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/)
!
!*       1.    POSITION FILE
!              -------------
!
REWIND(KULNAM)
ILEN=LEN(HDNAML)
IRET = 0
!
search_nam : DO
      YLINE=' '
      READ(UNIT=KULNAM,FMT='(A)',IOSTAT=IRET,END=100) YLINE
      IF (IRET /=0 ) THEN
        IF (PRESENT(KLUOUT)) &
         WRITE(KLUOUT,FMT=*) '-> error reading from unit',KULNAM,' line ',YLINE
      ELSE   
        ILEY=LEN(YLINE)
        DO J=1,ILEY
          DO JA=1,26
            IF (YLINE(J:J)==YLO(JA)) YLINE(J:J)=YUP(JA) 
          END DO
        END DO
        IND1=INDEX(YLINE,'&'//HDNAML)
        IF(IND1.NE.0) THEN
          YLTEST=YLINE(IND1+ILEN+1:IND1+ILEN+1)
          !IF((LLT(YLTEST,'0').OR.LGT(YLTEST,'9')).AND. &
          !   (LLT(YLTEST,'A').OR.LGT(YLTEST,'Z'))) EXIT search_nam
          IF(YLTEST == ' ') EXIT search_nam
        ENDIF
      ENDIF
ENDDO search_nam
!
BACKSPACE(KULNAM)
OFOUND=.TRUE.
IF (PRESENT(KLUOUT)) WRITE(KLUOUT,FMT=*) '-- namelist ',HDNAML,' read'
!
RETURN
!
! end of file: namelist name not found
100  CONTINUE
OFOUND=.FALSE.
IF (PRESENT(KLUOUT)) &
WRITE(KLUOUT,FMT=*)  &
'-- namelist ',HDNAML,' not found: default values used if required'
!------------------------------------------------------------------
END SUBROUTINE POSNAM
!!
!!
!!    ################################################
      SUBROUTINE POSKEY(KULNAM,KLUOUT,HKEYWD1,HKEYWD2)
!!    ################################################
!!
!!*** *POSKEY*
!!
!!    PURPOSE
!!    -------
!     To position namelist file at correct place after reading
!     keyword HKEYWD
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENT
!!    -----------------
!!
!!    REFERENCE
!!    ----------
!!
!!    AUTHOR
!!    -------
!!       I. Mallet *Meteo-France*
!!
!!    MODIFICATIONS
!!    --------------
!!       Original : 15/10/01
!------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!*       0.1   Declarations of arguments
!
INTEGER,                    INTENT(IN) :: KULNAM
INTEGER,                    INTENT(IN) :: KLUOUT
CHARACTER(LEN=*),           INTENT(IN) :: HKEYWD1
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: HKEYWD2
!
!*       0.2   Declarations of local variables
!
CHARACTER(LEN=120) :: YLINE
INTEGER            :: ILEN1,IRET
!
!
!*       1.    POSITION FILE
!              -------------
!
REWIND(KULNAM)
ILEN1=LEN(HKEYWD1)
IF (PRESENT(HKEYWD2)) ILEN2=LEN(HKEYWD2)
!
search_key : DO
      YLINE=' '
      READ(UNIT=KULNAM,FMT='(A)',IOSTAT=IRET,END=100) YLINE
      IF (IRET /=0 ) THEN
         WRITE(KLUOUT,FMT=*) '-> error when reading line from unit ',KULNAM
      ELSE
        YLINE=ADJUSTL(YLINE)
        IF (YLINE(1:ILEN1) .EQ. HKEYWD1(1:ILEN1)) EXIT search_key
      ENDIF
ENDDO search_key
!
WRITE(KLUOUT,FMT=*) '-- keyword ',HKEYWD1,' found'
!
RETURN
!
! end of file: keyword not found
100  CONTINUE
IF (.NOT.PRESENT(HKEYWD2)) THEN
  WRITE(KLUOUT,FMT=*)  '-- keyword ',HKEYWD1,' not found: program stop'
!callabortstop
CALL ABORT
  STOP
ELSE
!
!*       2.    SECOND KEYWORD: POSITION FILE
!              -----------------------------
!
  REWIND(KULNAM)
  search_key2 : DO
      YLINE=' '
      READ(UNIT=KULNAM,FMT='(A)',IOSTAT=IRET,END=101) YLINE
      IF (IRET /=0 ) THEN
        WRITE(KLUOUT,FMT=*) '-> error when reading line from unit ',KULNAM
      ELSE
        YLINE=ADJUSTL(YLINE)
        IF (YLINE(1:ILEN2) .EQ. HKEYWD2(1:ILEN2)) EXIT search_key2
      ENDIF
  ENDDO search_key2
  WRITE(KLUOUT,FMT=*) '-- keyword ',HKEYWD2,' found'
  RETURN
END IF
! end of file: scd keyword not found
101  CONTINUE
WRITE(KLUOUT,FMT=*)  '-- keyword ',HKEYWD2,' not found: program stop'
!callabortstop
CALL ABORT
STOP
!------------------------------------------------------------------
END SUBROUTINE POSKEY
!
END MODULE MODE_POS
