!     ######spl
      MODULE MODI_WRITE_DIMGRIDREF
!     ############################
!
INTERFACE
!
SUBROUTINE WRITE_DIMGRIDREF
END SUBROUTINE WRITE_DIMGRIDREF
!
END INTERFACE
!
END MODULE MODI_WRITE_DIMGRIDREF
!     ###########################
      SUBROUTINE WRITE_DIMGRIDREF
!     ###########################
!
!!****  *WRITE_DIMGRIDREF* - Ouverture du fichier diachronique et ecriture
!!          de l'entete
!! 
!!
!!    PURPOSE
!!    -------
! 
!
!!**  METHOD
!!    ------
!!      
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHORS
!!    -------
!!    J. Duron      *Lab. Aerologie* 
!!
!!    Copyright 1994,  Meteo-France and Laboratoire d'Aerologie
!!    All Rights Reserved
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/01/96 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIACHRO 
USE MODI_WRITE_LFIFM1_FORDIACHRO_CV

!
!*       0.1   Local variables
!
INTEGER           :: IRESP
!
!*       1.    Ouverture du fichier diachronique
!              ---------------------------------
!
CALL FMLOOK(CLUOUTDIA,CLUOUTDIA,NLUOUTDIA,IRESP)
IF (IRESP/=0) THEN
  ! ouverture du listing
  CALL FMATTR(CLUOUTDIA,CLUOUTDIA,NLUOUTDIA,NRESPDIA)
  OPEN(UNIT=NLUOUTDIA,FILE=CLUOUTDIA,FORM='FORMATTED')
END IF
!
WRITE(UNIT=NLUOUTDIA,FMT=1)CFILEDIA
1 FORMAT(' OPEN NEW DIACHRONIC FILE ',A28)

! Modif demandee par Nicole Asencio. 28/9/98
NFTYPEDIA=2
!NFTYPEDIA=0
NVERBDIA=5
CALL FMOPEN(CFILEDIA,'NEW',CLUOUTDIA,NNPRARDIA,NFTYPEDIA,NVERBDIA,NNINARDIA, &
	    NRESPDIA)

!
!*       2.    Fermeture du fichier descriptif correspondant
! et unite logique correspondante liberee
!              ----------------------------------------------
!
!  non, on ferme DES et LFI par FMCLOS a la fin du programme
!(On peut envisager d'y ecrire le DESFM des fichiers d'entree)
!
!*       3.    Ecriture des dimensions, parametres de grille, etat de ref...
!              ----------------------------------------------------------
!
CALL WRITE_LFIFM1_FORDIACHRO_CV(CFILEDIA)

!
!------------------------------------------------------------------------------
!
!*      4.    EPILOGUE
!             --------

RETURN

END SUBROUTINE WRITE_DIMGRIDREF
