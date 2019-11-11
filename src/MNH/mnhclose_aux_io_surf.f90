!MNH_LIC Copyright 2003-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_MNHCLOSE_AUX_IO_SURF
!     #########################
INTERFACE
      SUBROUTINE MNHCLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
!
CHARACTER(LEN=28), INTENT(IN), OPTIONAL :: HFILE    ! file to close
CHARACTER(LEN=6),  INTENT(IN), OPTIONAL :: HFILETYPE! type of file to close

!
END SUBROUTINE MNHCLOSE_AUX_IO_SURF
!
END INTERFACE
END MODULE MODI_MNHCLOSE_AUX_IO_SURF
!
!     #######################################################
      SUBROUTINE MNHCLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
!     #######################################################
!
!!****  *MNHCLOSE_AUX_IO_SURF* - routine to close IO files in MESONH universe 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	S.Malardel   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2003 
!!      J.Escobar : 19/04/2016 : Pb IOZ/NETCDF , missing OPARALLELIO=.FALSE. for PGD files
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 surfex 2006/05/23 15:47:28
!-----------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
USE MODE_FM
USE MODE_IO_ll

USE MODD_IO_SURF_MNH, ONLY : TPINFILE, CACTION, NMASK_ALL, NMASK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=28), INTENT(IN), OPTIONAL :: HFILE    ! file to close
CHARACTER(LEN=6),  INTENT(IN), OPTIONAL :: HFILETYPE! type of file to close

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! return-code if a problem appears
!
!-------------------------------------------------------------------------------
!
IF (CACTION=='OPEN  ') THEN
  CALL IO_FILE_CLOSE_ll(TPINFILE,OPARALLELIO=.FALSE.)
  CACTION='      '
END IF
!
DEALLOCATE(NMASK_ALL)
DEALLOCATE(NMASK)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MNHCLOSE_AUX_IO_SURF
