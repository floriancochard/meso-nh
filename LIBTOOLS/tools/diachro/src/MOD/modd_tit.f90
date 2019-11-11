!     ######spl
      MODULE  MODD_TIT
!     ################
!
!!****  *MODD_TIT* - 
!!
!!    PURPOSE
!!    -------
!
!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     original  JD    08/01/96
!!     updated   PM   
!!
!-------------------------------------------------------------------------
!
!*     0.   Declarations
!           ------------
!

IMPLICIT NONE


! Nom fichier diachronique
!
CHARACTER(LEN=10),SAVE  :: CTITALL
CHARACTER(LEN=100),SAVE :: CTITT1, CTITT2, CTITT3, CTITB1, CTITB2, CTITB3
CHARACTER(LEN=100),SAVE :: CTITB3MEM
CHARACTER(LEN=40),SAVE  :: CTITYT, CTITYM, CTITYB, CTITXL, CTITXM, CTITXR
CHARACTER(LEN=40),SAVE  :: CTITVAR1, CTITVAR2, CTITVAR3, CTITVAR4, CTITVAR5
CHARACTER(LEN=40),SAVE  :: CTITVAR6, CTITVAR7, CTITVAR8
LOGICAL                 :: LTITDEF, LTITDEFM
REAL,SAVE               :: XSZTITT1=0., XSZTITT2=0., XSZTITT3=0.
REAL,SAVE               :: XPOSTITT1=0., XPOSTITT2=0., XPOSTITT3=0.
REAL,SAVE               :: XYPOSTITT1=0., XYPOSTITT2=0., XYPOSTITT3=0.
!
REAL,SAVE               :: XSZTITB1=0., XSZTITB2=0., XSZTITB3=0.
REAL,SAVE               :: XPOSTITB1=0., XPOSTITB2=0., XPOSTITB3=0.
REAL,SAVE               :: XYPOSTITB1=0., XYPOSTITB2=0., XYPOSTITB3=0.
!
REAL,SAVE               :: XSZTITYT=0., XSZTITYM=0., XSZTITYB=0.
REAL,SAVE               :: XPOSTITYT=0., XPOSTITYM=0., XPOSTITYB=0.
REAL,SAVE               :: XYPOSTITYT=0., XYPOSTITYM=0., XYPOSTITYB=0.
!
REAL,SAVE               :: XSZTITVAR1=0., XSZTITVAR2=0., XSZTITVAR3=0.
REAL,SAVE               :: XSZTITVAR4=0., XSZTITVAR5=0., XSZTITVAR6=0.
REAL,SAVE               :: XSZTITVAR7=0., XSZTITVAR8=0.
REAL,SAVE               :: XPOSTITVAR1=0., XPOSTITVAR2=0., XPOSTITVAR3=0.
REAL,SAVE               :: XPOSTITVAR4=0., XPOSTITVAR5=0., XPOSTITVAR6=0.
REAL,SAVE               :: XPOSTITVAR7=0., XPOSTITVAR8=0.
REAL,SAVE               :: XYPOSTITVAR1=0., XYPOSTITVAR2=0., XYPOSTITVAR3=0.
REAL,SAVE               :: XYPOSTITVAR4=0., XYPOSTITVAR5=0., XYPOSTITVAR6=0.
REAL,SAVE               :: XYPOSTITVAR7=0., XYPOSTITVAR8=0.

!
END MODULE MODD_TIT
