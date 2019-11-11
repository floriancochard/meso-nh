!     ######spl
      MODULE MODD_PT_FOR_CH_FORDIACHRO
!     ################################
!
!!****  *MODD_PT_FOR_CH_FORDIACHRO* - Global variables  NMT declaration
!!                                          XPRES  "
!!                                          XPHI   "
!!                                          XTH    "
!!
!!    PURPOSE
!!    -------
!!      This declarative module defines the NMT global variable, which
!!    takes the value 1 for variables at t-dt time  and 2 for variables
!!    at t time 
!!    XPRES contains the pressure value computed  either at t-dt or t
!!    times for constant pressure sections processing and RS.
!!    XTH contains either XTHM or XTHT
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/07/96                      
!!      Updated   PM 
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
INTEGER           :: NMT         
INTEGER           :: NLOOPT      

REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE,SAVE :: XPRES, XPHI, XTH
LOGICAL           :: LTHSTAB=.TRUE. ! flag to write possibly 'UNSTABLE THETA' message
! Ajout pour RS
! 'CART'
REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE,SAVE :: XU, XV, XRVJD
REAL,DIMENSION(:),ALLOCATABLE,SAVE :: XTIMRS
! 'CART' + 'RSPL'
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XTRS, XPRS, XURS, XVRS, XRVRS
! 'RSPL'
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XTIMRS2
INTEGER,DIMENSION(:),ALLOCATABLE,SAVE :: NST, NNST
!
! Ajout pour la composante W (Cas CV ULMWM et ULTWT)
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE :: XWCV
! 
! Septembre 2000 Pour =/= UMVM ou MUMVM ....
REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE,SAVE :: XUMEM, XVMEM
!
END MODULE MODD_PT_FOR_CH_FORDIACHRO
