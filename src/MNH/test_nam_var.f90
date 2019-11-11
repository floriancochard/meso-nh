!MNH_LIC Copyright 1998-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!############################
MODULE MODI_TEST_NAM_VAR
!############################
!
INTERFACE TEST_NAM_VAR
!
      SUBROUTINE TEST_NAM_VARC0(KLUOUT,HNAME,HVAR,            &
                                     HVALUE1,HVALUE2,HVALUE3, &
                                     HVALUE4,HVALUE5,HVALUE6, &
                                     HVALUE7,HVALUE8,HVALUE9, &
                                     HVALUE10,HVALUE11,HVALUE12   )
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
CHARACTER(LEN=*) ,INTENT(IN)           ::HVAR     ! variable to test

CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE1  ! first possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE2  ! second possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE3  ! third possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE4  ! fourth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE5  ! fiveth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE6  ! sixth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE7  ! seventh possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE8  ! eightth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE9  ! nineth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE10 ! tenth possible value             
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE11 ! eleventh possible value             
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE12 ! twelveth possible value             
!
END SUBROUTINE TEST_NAM_VARC0
!
END INTERFACE
!
END MODULE MODI_TEST_NAM_VAR
!
!
!     #########################################################
      SUBROUTINE TEST_NAM_VARC0(KLUOUT,HNAME,HVAR,            &
                                     HVALUE1,HVALUE2,HVALUE3, &
                                     HVALUE4,HVALUE5,HVALUE6, &
                                     HVALUE7,HVALUE8,HVALUE9, &
                                     HVALUE10,HVALUE11,HVALUE12   )
!     #########################################################
!
!!****  *TEST_NAM_VARC0* - routine to test the value of a character var.
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     17/04/98
!!      10/2016 (C.Lac) Increase of the number of values
!!      P.Wautelet 22/01/2019: use standard FLUSH statement instead of non standard intrinsics
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
INTEGER,          INTENT(IN)           ::KLUOUT   ! output listing logical unit
CHARACTER(LEN=*) ,INTENT(IN)           ::HNAME    ! name of the variable to test
CHARACTER(LEN=*) ,INTENT(IN)           ::HVAR     ! variable to test

CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE1  ! first possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE2  ! second possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE3  ! third possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE4  ! fourth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE5  ! fiveth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE6  ! sixth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE7  ! seventh possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE8  ! eightth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE9  ! nineth possible value
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE10 ! tenth possible value           
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE11 ! eleventh possible value             
CHARACTER(LEN=*) ,INTENT(IN), OPTIONAL ::HVALUE12 ! twelveth possible value             
!
!*      0.2   Declarations of local variables
!
!
!-------------------------------------------------------------------------------
!
IF ( PRESENT (HVALUE1) ) THEN
  IF ( HVAR==HVALUE1 ) RETURN
END IF
!
IF ( PRESENT (HVALUE2) ) THEN
  IF ( HVAR==HVALUE2 ) RETURN
END IF
!
IF ( PRESENT (HVALUE3) ) THEN
  IF ( HVAR==HVALUE3 ) RETURN
END IF
!
IF ( PRESENT (HVALUE4) ) THEN
  IF ( HVAR==HVALUE4 ) RETURN
END IF
!
IF ( PRESENT (HVALUE5) ) THEN
  IF ( HVAR==HVALUE5 ) RETURN
END IF
!
IF ( PRESENT (HVALUE6) ) THEN
  IF ( HVAR==HVALUE6 ) RETURN
END IF
!
IF ( PRESENT (HVALUE7) ) THEN
  IF ( HVAR==HVALUE7 ) RETURN
END IF
!
IF ( PRESENT (HVALUE8) ) THEN
  IF ( HVAR==HVALUE8 ) RETURN
END IF
!
IF ( PRESENT (HVALUE9) ) THEN
  IF ( HVAR==HVALUE9 ) RETURN
END IF
!
IF ( PRESENT (HVALUE10) ) THEN  
  IF ( HVAR==HVALUE10 ) RETURN
END IF
!
IF ( PRESENT (HVALUE11) ) THEN  
  IF ( HVAR==HVALUE11 ) RETURN
END IF
!
IF ( PRESENT (HVALUE12) ) THEN  
  IF ( HVAR==HVALUE12 ) RETURN
END IF
!
!
!-------------------------------------------------------------------------------
!
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'FATAL ERROR:'
WRITE (KLUOUT,*) '-----------'
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Value "',HVAR,'" is not allowed for variable ',HNAME
WRITE (KLUOUT,*) ' '
WRITE (KLUOUT,*) 'Possible values are:'
IF ( PRESENT (HVALUE1) ) WRITE (KLUOUT,*) '"',HVALUE1,'"'
IF ( PRESENT (HVALUE2) ) WRITE (KLUOUT,*) '"',HVALUE2,'"'
IF ( PRESENT (HVALUE3) ) WRITE (KLUOUT,*) '"',HVALUE3,'"'
IF ( PRESENT (HVALUE4) ) WRITE (KLUOUT,*) '"',HVALUE4,'"'
IF ( PRESENT (HVALUE5) ) WRITE (KLUOUT,*) '"',HVALUE5,'"'
IF ( PRESENT (HVALUE6) ) WRITE (KLUOUT,*) '"',HVALUE6,'"'
IF ( PRESENT (HVALUE7) ) WRITE (KLUOUT,*) '"',HVALUE7,'"'
IF ( PRESENT (HVALUE8) ) WRITE (KLUOUT,*) '"',HVALUE8,'"'
IF ( PRESENT (HVALUE9) ) WRITE (KLUOUT,*) '"',HVALUE9,'"'
IF ( PRESENT (HVALUE10) ) WRITE (KLUOUT,*) '"',HVALUE10,'"'
IF ( PRESENT (HVALUE11) ) WRITE (KLUOUT,*) '"',HVALUE11,'"'
IF ( PRESENT (HVALUE12) ) WRITE (KLUOUT,*) '"',HVALUE12,'"'
FLUSH(unit=KLUOUT)
!
 !callabortstop
CALL ABORT
STOP
!-------------------------------------------------------------------------------
END SUBROUTINE TEST_NAM_VARC0
