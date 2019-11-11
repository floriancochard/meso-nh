!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_rea_lfi.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_REA_LFI
!     ###################
!
!!****  *MODD_REA_LFI* - Defines a LFIFM file record
!!
!!    PURPOSE
!!    -------
!       This declarative module globally defines the set of variables 
!     controlling the recors of the LFIFM file.
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      - Fischer C., 1994, "File structure and content in the Meso-NH 
!!        model", Meso-nh internal note, CNRM/GMME,  July 5.
!!       
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/06/94                      
!!      Updated  PM 22/11/94
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
INTEGER           :: NRESP,NMELEV   ! NRESP  : return-code if a problem appears 
                                    !  opening the file
                                    ! NMELEV : level of message printing in
                                    !  LFI  routines 
INTEGER           :: NNPRAR,NFTYPE  ! NNPRAR : number of predicted articles
                                    ! NFTYPE : type of FM-file for FMCLOS 
INTEGER           :: NNINAR         !  number of articles initially present in
                                    !  the file
INTEGER           :: NGRID,NLENG    ! NGRID : grid indicator
                                    ! NLENG : length of the data field  
INTEGER           :: NLENCH         ! NLENCH : length of comment string 
!
CHARACTER(LEN=3)  :: CSTATU         ! Status of the file before the open
CHARACTER(LEN=16) :: CRECFM         ! Name of the article to be written
CHARACTER(LEN=100):: CCOMMENT       ! Comment string
!
LOGICAL           :: LFATER,LSTATS  ! LFATER : true if LFI-file manipulation 
                                    !  error is a fatal error 
                                    ! LSTATS : true if statistics of file
                                    !  manipulation sould be printed
END MODULE MODD_REA_LFI
