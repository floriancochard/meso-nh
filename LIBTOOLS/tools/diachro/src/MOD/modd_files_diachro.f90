!     ######spl
      MODULE  MODD_FILES_DIACHRO
!     ##########################
!
!!****  *MODD_FILES_DIACHRO* - 
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


! 
INTEGER,SAVE       :: NBGUIL
INTEGER,DIMENSION(180),SAVE     :: NMGUIL

!
! Nb fichiers ouverts
!
INTEGER,SAVE       :: NBFILES=0
! Numero x de _filex_ du fichier courant
INTEGER,SAVE       :: NUMFILECUR
! Memorisation du fichier courant dans le cas de diffrerence de 2 champs
INTEGER,SAVE       :: NUMFILECUR2
! cf JPNXFM (modd_fmdeclar)= limite Fortran (99=JPNXLU) -10
INTEGER,DIMENSION(90),SAVE     :: NUMFILES

!
! Indication traitement un seul fichier ou plusieurs fichiers simultanement
!
LOGICAL            :: LFIC1=.TRUE.
!
! Plusieurs fichiers simultanes
!
INTEGER,SAVE       :: NBSIMULT
INTEGER,DIMENSION(90),SAVE     :: NUMFILESIMULT
INTEGER,DIMENSION(90),SAVE     :: NINDFILESIMULT


! Nom fichiers diachroniques
!
CHARACTER(LEN=28),DIMENSION(90),SAVE :: CFILEDIAS                
!
! Listings associes au traitement diachronique
!
CHARACTER(LEN=16),DIMENSION(90),SAVE :: CLUOUTDIAS='OUT_DIA'               
!
! Numeros logiques des  listings et parametres d'ouverture des fichiers
!
INTEGER,DIMENSION(90),SAVE           :: NLUOUTDIAS, NNPRARDIAS, NFTYPEDIAS=2, NVERBDIAS,  &
                          NNINARDIAS, NRESPDIAS

!
END MODULE MODD_FILES_DIACHRO
