!     ######spl
      MODULE MODI_RESOLV_TIT
!     ######################
!
INTERFACE
!
SUBROUTINE RESOLV_TIT(HTIT,HOUT)
CHARACTER(LEN=*)  :: HTIT, HOUT
END SUBROUTINE RESOLV_TIT
!
END INTERFACE
END MODULE MODI_RESOLV_TIT
!     ######spl
      SUBROUTINE RESOLV_TIT(HTIT,HOUT)
!     ################################
!
!!****  *RESOLV_TIT* - 
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
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO
USE MODD_TIT

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HTIT, HOUT
!
!*       0.1   Local variables
!              ---------------


!
!------------------------------------------------------------------------------
!print *,' RESOLV_TIT HTIT HOUT ',HTIT,HOUT
IF(.NOT.LTITDEF)THEN
  CTITALL='NODEFAULT'
ELSE
  CTITALL='DEFAULT'
ENDIF
!print *,' RESOLV_TIT CTITALL',CTITALL  
IF(CTITALL == 'DEFAULT' .OR. CTITALL == 'default' .OR. &
   CTITALL == 'DEFAUT'  .OR. CTITALL == 'defaut')THEN
   CTITT1(1:LEN(CTITT1))=' '
   CTITT1='DEFAULT'
   CTITT2(1:LEN(CTITT2))=' '
   CTITT2='DEFAULT'
   CTITT3(1:LEN(CTITT3))=' '
   CTITT3='DEFAULT'
   CTITB1(1:LEN(CTITB1))=' '
   CTITB1='DEFAULT'
   CTITB2(1:LEN(CTITB2))=' '
   CTITB2='DEFAULT'
   CTITB3(1:LEN(CTITB3))=' '
   CTITB3='DEFAULT'
   CTITYT(1:LEN(CTITYT))=' '
   CTITYT='DEFAULT'
   CTITYM(1:LEN(CTITYM))=' '
   CTITYM='DEFAULT'
   CTITYB(1:LEN(CTITYB))=' '
   CTITYB='DEFAULT'
   CTITXL(1:LEN(CTITXL))=' '
   CTITXL='DEFAULT'
   CTITXM(1:LEN(CTITXM))=' '
   CTITXM='DEFAULT'
   CTITXR(1:LEN(CTITXR))=' '
   CTITXR='DEFAULT'
   CTITVAR1(1:LEN(CTITVAR1))=' '
   CTITVAR1='DEFAULT'
   CTITVAR2(1:LEN(CTITVAR2))=' '
   CTITVAR2='DEFAULT'
   CTITVAR3(1:LEN(CTITVAR3))=' '
   CTITVAR3='DEFAULT'
   CTITVAR4(1:LEN(CTITVAR4))=' '
   CTITVAR4='DEFAULT'
   CTITVAR5(1:LEN(CTITVAR5))=' '
   CTITVAR5='DEFAULT'
   CTITVAR6(1:LEN(CTITVAR6))=' '
   CTITVAR6='DEFAULT'
   CTITVAR7(1:LEN(CTITVAR7))=' '
   CTITVAR7='DEFAULT'
   CTITVAR8(1:LEN(CTITVAR8))=' '
   CTITVAR8='DEFAULT'
ELSE
! print *,' HTIT '
  IF(.NOT.LTITDEF)THEN
  SELECT CASE(HTIT)
    CASE('CTITT1')
      IF(CTITT1 == 'WHITE' .OR. CTITT1 == 'white' .OR.  &
	 CTITT1 == 'BLANC' .OR. CTITT1 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITT1 == 'CCOMMENT' .OR.   CTITT1 == 'ccomment' .OR. &
               CTITT1 == 'COMMENT' .OR.  CTITT1 == 'comment')THEN
        CTITT1=ADJUSTL(ADJUSTR(CCOMMENT(NLOOPP)))
	HOUT=CTITT1
	CTITALL='NODEFAULT'
      ELSE  IF(CTITT1 == 'DEFAULT' .OR.   CTITT1 == 'default' .OR. &
               CTITT1 == 'DEFAUT' .OR.  CTITT1 == 'defaut')THEN
      ELSE
        HOUT=CTITT1
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITT2')
      IF(CTITT2 == 'WHITE' .OR. CTITT2 == 'white' .OR.  &
	 CTITT2 == 'BLANC' .OR. CTITT2 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITT2 == 'CCOMMENT' .OR.   CTITT2 == 'ccomment' .OR. &
               CTITT2 == 'COMMENT' .OR.  CTITT2 == 'comment')THEN
        CTITT2=ADJUSTL(ADJUSTR(CCOMMENT(NLOOPP)))
	HOUT=CTITT2
	CTITALL='NODEFAULT'
      ELSE  IF(CTITT2 == 'DEFAULT' .OR.   CTITT2 == 'default' .OR. &
               CTITT2 == 'DEFAUT' .OR.  CTITT2 == 'defaut')THEN
      ELSE
        HOUT=CTITT2
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITT3')
      IF(CTITT3 == 'WHITE' .OR. CTITT3 == 'white' .OR.  &
	 CTITT3 == 'BLANC' .OR. CTITT3 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITT3 == 'CCOMMENT' .OR.   CTITT3 == 'ccomment' .OR. &
               CTITT3 == 'COMMENT' .OR.  CTITT3 == 'comment')THEN
        CTITT3=ADJUSTL(ADJUSTR(CCOMMENT(NLOOPP)))
	HOUT=CTITT3
	CTITALL='NODEFAULT'
      ELSE  IF(CTITT3 == 'DEFAULT' .OR.   CTITT3 == 'default' .OR. &
               CTITT3 == 'DEFAUT' .OR.  CTITT3 == 'defaut')THEN
      ELSE
        HOUT=CTITT3
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITB1')
      IF(CTITB1 == 'WHITE' .OR. CTITB1 == 'white' .OR.  &
	 CTITB1 == 'BLANC' .OR. CTITB1 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITB1 == 'CCOMMENT' .OR.   CTITB1 == 'ccomment' .OR. &
               CTITB1 == 'COMMENT' .OR.  CTITB1 == 'comment')THEN
        CTITB1=ADJUSTL(ADJUSTR(CCOMMENT(NLOOPP)))
	HOUT=CTITB1
	CTITALL='NODEFAULT'
      ELSE  IF(CTITB1 == 'DEFAULT' .OR.   CTITB1 == 'default' .OR. &
               CTITB1 == 'DEFAUT' .OR.  CTITB1 == 'defaut')THEN
      ELSE
        HOUT=CTITB1
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
!     print *,' HOUT ',HOUT
      RETURN
    CASE('CTITB2')
      IF(CTITB2 == 'WHITE' .OR. CTITB2 == 'white' .OR.  &
	 CTITB2 == 'BLANC' .OR. CTITB2 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITB2 == 'CCOMMENT' .OR.   CTITB2 == 'ccomment' .OR. &
               CTITB2 == 'COMMENT' .OR.  CTITB2 == 'comment')THEN
        CTITB2=ADJUSTL(ADJUSTR(CCOMMENT(NLOOPP)))
	HOUT=CTITB2
	CTITALL='NODEFAULT'
      ELSE  IF(CTITB2 == 'DEFAULT' .OR.   CTITB2 == 'default' .OR. &
               CTITB2 == 'DEFAUT' .OR.  CTITB2 == 'defaut')THEN
      ELSE
        HOUT=CTITB2
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
!     print *,' HOUT ',HOUT
      RETURN
    CASE('CTITB3')
      IF(CTITB3 == 'WHITE' .OR. CTITB3 == 'white' .OR.  &
	 CTITB3 == 'BLANC' .OR. CTITB3 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITB3 == 'CCOMMENT' .OR.   CTITB3 == 'ccomment' .OR. &
               CTITB3 == 'COMMENT' .OR.  CTITB3 == 'comment')THEN
        CTITB3=ADJUSTL(ADJUSTR(CCOMMENT(NLOOPP)))
	HOUT=CTITB3
	CTITALL='NODEFAULT'
      ELSE  IF(CTITB3 == 'DEFAULT' .OR.   CTITB3 == 'default' .OR. &
               CTITB3 == 'DEFAUT' .OR.  CTITB3 == 'defaut')THEN
      ELSE
        HOUT=CTITB3
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
!     print *,' HOUT ',HOUT
      RETURN
    CASE('CTITYT')
      IF(CTITYT == 'WHITE' .OR. CTITYT == 'white' .OR.  &
	 CTITYT == 'BLANC' .OR. CTITYT == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITYT == 'DEFAULT' .OR.   CTITYT == 'default' .OR. &
               CTITYT == 'DEFAUT' .OR.  CTITYT == 'defaut')THEN
      ELSE
        HOUT=CTITYT
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITYM')
      IF(CTITYM == 'WHITE' .OR. CTITYM == 'white' .OR.  &
	 CTITYM == 'BLANC' .OR. CTITYM == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITYM == 'DEFAULT' .OR.   CTITYM == 'default' .OR. &
               CTITYM == 'DEFAUT' .OR.  CTITYM == 'defaut')THEN
      ELSE
        HOUT=CTITYM
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITYB')
      IF(CTITYB == 'WHITE' .OR. CTITYB == 'white' .OR.  &
	 CTITYB == 'BLANC' .OR. CTITYB == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITYB == 'DEFAULT' .OR.   CTITYB == 'default' .OR. &
               CTITYB == 'DEFAUT' .OR.  CTITYB == 'defaut')THEN
      ELSE
        HOUT=CTITYB
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITXL')
      IF(CTITXL == 'WHITE' .OR. CTITXL == 'white' .OR.  &
	 CTITXL == 'BLANC' .OR. CTITXL == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITXL == 'DEFAULT' .OR.   CTITXL == 'default' .OR. &
               CTITXL == 'DEFAUT' .OR.  CTITXL == 'defaut')THEN
      ELSE
        HOUT=CTITXL
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITXM')
      IF(CTITXM == 'WHITE' .OR. CTITXM == 'white' .OR.  &
	 CTITXM == 'BLANC' .OR. CTITXM == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITXM == 'DEFAULT' .OR.   CTITXM == 'default' .OR. &
               CTITXM == 'DEFAUT' .OR.  CTITXM == 'defaut')THEN
      ELSE
        HOUT=CTITXM
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITXR')
      IF(CTITXR == 'WHITE' .OR. CTITXR == 'white' .OR.  &
	 CTITXR == 'BLANC' .OR. CTITXR == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITXR == 'DEFAULT' .OR.   CTITXR == 'default' .OR. &
               CTITXR == 'DEFAUT' .OR.  CTITXR == 'defaut')THEN
      ELSE
        HOUT=CTITXR
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR1')
      IF(CTITVAR1 == 'WHITE' .OR. CTITVAR1 == 'white' .OR.  &
	 CTITVAR1 == 'BLANC' .OR. CTITVAR1 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR1 == 'DEFAULT' .OR.   CTITVAR1 == 'default' .OR. &
               CTITVAR1 == 'DEFAUT' .OR.  CTITVAR1 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR1
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR2')
      IF(CTITVAR2 == 'WHITE' .OR. CTITVAR2 == 'white' .OR.  &
	 CTITVAR2 == 'BLANC' .OR. CTITVAR2 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR2 == 'DEFAULT' .OR.   CTITVAR2 == 'default' .OR. &
               CTITVAR2 == 'DEFAUT' .OR.  CTITVAR2 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR2
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR3')
      IF(CTITVAR3 == 'WHITE' .OR. CTITVAR3 == 'white' .OR.  &
	 CTITVAR3 == 'BLANC' .OR. CTITVAR3 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR3 == 'DEFAULT' .OR.   CTITVAR3 == 'default' .OR. &
               CTITVAR3 == 'DEFAUT' .OR.  CTITVAR3 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR3
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR4')
      IF(CTITVAR4 == 'WHITE' .OR. CTITVAR4 == 'white' .OR.  &
	 CTITVAR4 == 'BLANC' .OR. CTITVAR4 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR4 == 'DEFAULT' .OR.   CTITVAR4 == 'default' .OR. &
               CTITVAR4 == 'DEFAUT' .OR.  CTITVAR4 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR4
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR5')
      IF(CTITVAR5 == 'WHITE' .OR. CTITVAR5 == 'white' .OR.  &
	 CTITVAR5 == 'BLANC' .OR. CTITVAR5 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR5 == 'DEFAULT' .OR.   CTITVAR5 == 'default' .OR. &
               CTITVAR5 == 'DEFAUT' .OR.  CTITVAR5 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR5
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR6')
      IF(CTITVAR6 == 'WHITE' .OR. CTITVAR6 == 'white' .OR.  &
	 CTITVAR6 == 'BLANC' .OR. CTITVAR6 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR6 == 'DEFAULT' .OR.   CTITVAR6 == 'default' .OR. &
               CTITVAR6 == 'DEFAUT' .OR.  CTITVAR6 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR6
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR7')
      IF(CTITVAR7 == 'WHITE' .OR. CTITVAR7 == 'white' .OR.  &
	 CTITVAR7 == 'BLANC' .OR. CTITVAR7 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR7 == 'DEFAULT' .OR.   CTITVAR7 == 'default' .OR. &
               CTITVAR7 == 'DEFAUT' .OR.  CTITVAR7 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR7
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
    CASE('CTITVAR8')
      IF(CTITVAR8 == 'WHITE' .OR. CTITVAR8 == 'white' .OR.  &
	 CTITVAR8 == 'BLANC' .OR. CTITVAR8 == 'blanc')THEN
        HOUT(1:LEN(HOUT))=' '
        CTITALL='NODEFAULT'
      ELSE  IF(CTITVAR8 == 'DEFAULT' .OR.   CTITVAR8 == 'default' .OR. &
               CTITVAR8 == 'DEFAUT' .OR.  CTITVAR8 == 'defaut')THEN
      ELSE
        HOUT=CTITVAR8
        CTITALL='NODEFAULT'
      ENDIF
!fuji      HOUT=ADJUSTL(HOUT)
      HOUT=TRIM(HOUT)
      RETURN
  END SELECT
ENDIF
ENDIF
RETURN
END  SUBROUTINE RESOLV_TIT
