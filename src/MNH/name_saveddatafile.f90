!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_pgd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##############################
      MODULE MODI_NAME_SAVEDDATAFILE
!     ##############################
INTERFACE
      SUBROUTINE NAME_SAVEDDATAFILE(HFILENAME,HFIELD,HSAVEDDATAFILE,ODATASAVE)
!
CHARACTER(LEN=28), INTENT(IN)  :: HFILENAME      ! Name of the global file.
CHARACTER(LEN=10), INTENT(IN)  :: HFIELD         ! field being treated
CHARACTER(LEN=28), INTENT(OUT) :: HSAVEDDATAFILE ! name of the file containing 
!                                                ! the primary data, but on the
!                                                ! MESONH domain only
LOGICAL,           INTENT(INOUT) :: ODATASAVE    ! flag to save data
!                                                ! on the local file.
!
END SUBROUTINE NAME_SAVEDDATAFILE
END INTERFACE
END MODULE MODI_NAME_SAVEDDATAFILE
!
!
!     #############################################################|##########
      SUBROUTINE NAME_SAVEDDATAFILE(HFILENAME,HFIELD,HSAVEDDATAFILE,ODATASAVE)
!     ########################################################################
!
!!**** *NAME_SAVEDDATAFILE* 
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!
!!    EXTERNAL
!!    --------
!!                           
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    23/07/98
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PGDGRID
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=28), INTENT(IN)  :: HFILENAME      ! Name of the global file.
CHARACTER(LEN=10), INTENT(IN)  :: HFIELD         ! field being treated
CHARACTER(LEN=28), INTENT(OUT) :: HSAVEDDATAFILE ! name of the file containing 
!                                                ! the primary data, but on the
!                                                ! MESONH domain only
LOGICAL,           INTENT(INOUT) :: ODATASAVE    ! flag to save data
!                                                ! on the local file.
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
CHARACTER(LEN=6) :: YSTRING 
!----------------------------------------------------------------------------
WRITE(YSTRING,'(I3.3,F3.2)')  INT(XPGDXHAT(2)-XPGDXHAT(1)+1.E-10) ,&
                                  XPGDXHAT(2)-XPGDXHAT(1)+1.E-10   &
                            - INT(XPGDXHAT(2)-XPGDXHAT(1)+1.E-10)
HSAVEDDATAFILE='data.pgd'//YSTRING//'km.'//ADJUSTL(HFIELD)
!
IF (HSAVEDDATAFILE==HFILENAME) ODATASAVE=.FALSE.
!-------------------------------------------------------------------------------
!
END SUBROUTINE NAME_SAVEDDATAFILE
