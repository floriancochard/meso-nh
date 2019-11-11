!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ################################################################
      SUBROUTINE INI_ELEC(KMI,HINIFILE,HLUOUT,PTSTEP,PDZMIN,KSPLITR, &
                          PDXX,PDYY,PDZZ,PDZX,PDZY                   )
!     ################################################################
!
!!****  *INI_ELEC* - routine to initialize the electrical parameters 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize the variables
!     of the atmospheric electricity scheme
!
!!**  METHOD
!!    ------
!!      The initialization of the scheme is performed as follows :
!!   
!!    EXTERNAL
!!    --------
!!      CLEANLIST_ll : deaalocate a list
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_ELEC)
!!      
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     29/11/02
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODE_IO_ll
USE MODE_FM
USE MODE_FMREAD
!
USE MODD_LUNIT_n
USE MODD_NSV, ONLY : NSV,NSV_ELEC,NSV_ELECBEG,NSV_ELECEND
USE MODD_PARAMETERS
USE MODD_CST
USE MODD_CONF
USE MODD_DYN
USE MODD_REF
USE MODD_TIME
!
USE MODN_CONF_n
!
USE MODI_INI_CLOUD
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
INTEGER,            INTENT(IN)   :: KMI      ! Model Index 
CHARACTER (LEN=*) , INTENT(IN)   :: HINIFILE ! name of the initial file
CHARACTER (LEN=*),  INTENT(IN)   :: HLUOUT   ! name for output-listing
                                             !  of nested models
REAL,               INTENT(IN)   :: PTSTEP   ! Time STEP 
!
REAL,               INTENT(IN)   :: PDZMIN   ! minimun vertical mesh size
INTEGER,            INTENT(IN)   :: KSPLITR  ! Number of small time step
                                             ! integration for  rain
                                             ! sedimendation
!
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZX    ! metric coefficient dzx
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
!
!*       0.2   declarations of local variables
!
INTEGER             :: IRESP   ! Return code of FM routines 
INTEGER             :: ILUOUT  ! Logical unit number of output-listing
!
!
!-------------------------------------------------------------------------------
!
!*       0.    PROLOGUE
!              --------
!
!
PRINT *,' INI_ELEC IS NOT YET DEVELOPPED'
!
!callabortstop
CALL ABORT
STOP
!
!-------------------------------------------------------------------------------
!
! 
END SUBROUTINE INI_ELEC
