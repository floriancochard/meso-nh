!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_title.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #################
      MODULE MODD_RADAR
!     #################
!
!!****  *MODD_RADAR* - Declare des variables concernant un (ou plusieurs)
!!                     radars
!!
!!    PURPOSE
!!    -------
!       Definit des variables pour localiser 1 ou +sieurs radars et
!     materialiser leur portee par un ou +esiurs cercles concentriques        
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
!!      Original    23/04/2003                    
!!      Updated  
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
LOGICAL           :: LRADAR=.FALSE., LRADIST=.FALSE., LRADRAY=.FALSE.
REAL              :: XLATRAD1=0., XLONRAD1=0.
REAL              :: XLATRAD2=0., XLONRAD2=0.
REAL              :: XLATRAD3=0., XLONRAD3=0.
REAL              :: XLATRAD4=0., XLONRAD4=0.
REAL,DIMENSION(6) :: XPORTRAD1,XPORTRAD2,XPORTRAD3,XPORTRAD4
REAL,DIMENSION(6) :: XLWRAD1=2.,XLWRAD2=2.,XLWRAD3=2.,XLWRAD4=2.
INTEGER           :: NPORTRAD1, NPORTRAD2, NPORTRAD3, NPORTRAD4
INTEGER           :: NLWRAD1, NLWRAD2, NLWRAD3, NLWRAD4
CHARACTER(LEN=1)  :: CSYMRAD1='+',CSYMRAD2='+',CSYMRAD3='+',CSYMRAD4='+'
END MODULE MODD_RADAR
