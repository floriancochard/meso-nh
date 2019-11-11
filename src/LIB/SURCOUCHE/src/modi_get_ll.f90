!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!     ##################
      MODULE MODI_GET_ll 
!     ##################
!
INTERFACE
!
!!     ##################################################
       SUBROUTINE  GET_DIM_EXT_ll( HSPLIT, KXDIM, KYDIM )
!!     ##################################################
!
  CHARACTER*1, INTENT(IN) :: HSPLIT
  INTEGER, INTENT(OUT) :: KXDIM, KYDIM
!
       END SUBROUTINE GET_DIM_EXT_ll
!
!!     ###################################################
       SUBROUTINE  GET_DIM_PHYS_ll( HSPLIT, KXDIM, KYDIM )
!!     ###################################################
!
  CHARACTER*1, INTENT(IN) :: HSPLIT
  INTEGER, INTENT(OUT) :: KXDIM, KYDIM
!
       END SUBROUTINE GET_DIM_PHYS_ll
!
!!     ##########################################
       SUBROUTINE GET_OR_ll( HSPLIT, KXOR, KYOR )
!!     ##########################################
!
  CHARACTER*1, INTENT(IN) :: HSPLIT
  INTEGER, INTENT(OUT) :: KXOR, KYOR
!
       END SUBROUTINE GET_OR_ll
!
!!     ####################################################
       SUBROUTINE GET_INDICE_ll( KXOR, KYOR, KXEND, KYEND )
!!     ####################################################
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
!
       END SUBROUTINE GET_INDICE_ll
!
!!     ##########################################
       SUBROUTINE GET_GLOBALDIMS_ll(KIMAX, KJMAX)
!!     ##########################################
!
  INTEGER, INTENT(OUT) :: KIMAX, KJMAX ! current model dimensions
!
       END SUBROUTINE GET_GLOBALDIMS_ll
!
!!     ######################################################
       SUBROUTINE GET_PHYSICAL_ll( KXOR, KYOR, KXEND, KYEND )
!!     ######################################################
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
!
       END SUBROUTINE GET_PHYSICAL_ll
!
!!     ###############################################################
       SUBROUTINE GET_INTERSECTION_ll( KXOR, KYOR, KXEND, KYEND, &
                                       KXORI, KYORI, KXENDI, KYENDI, &
                                       HDOM, KINFO, KIP )
!!     ###############################################################
!
  CHARACTER(LEN=4), INTENT(IN) :: HDOM    ! 'EXTE' for extended subdomain
                                          ! 'PHYS' for physical subdomain
  INTEGER, INTENT(IN)  :: KXOR, KYOR, &   ! Coordinates of the
                          KXEND, KYEND    ! region
!
  INTEGER, INTENT(OUT) :: KXORI, KYORI, & ! Global Coordinates
                          KXENDI, KYENDI  ! of the intersection
  INTEGER, INTENT(OUT) :: KINFO           ! Returned Info
  INTEGER, INTENT(IN), OPTIONAL:: KIP     ! Processor number
                                          ! (or subdomain number)
       END SUBROUTINE GET_INTERSECTION_ll
!
END INTERFACE
!
INTERFACE GET_GLOBALSLICE_ll
!
!!      ####################################################################
        SUBROUTINE GET_1DGLOBALSLICE_ll( PARRAY, HDIR, KLOC, PGLOBALSLICE, &
                                         KB, KE, KERR )
!!      ####################################################################
!
  REAL, DIMENSION(:,:), TARGET, INTENT(IN) :: PARRAY
  CHARACTER(LEN=1), INTENT(IN)      :: HDIR
  INTEGER, INTENT(IN)               :: KLOC
  REAL, DIMENSION(:), INTENT(OUT)   :: PGLOBALSLICE
  INTEGER, OPTIONAL    :: KB, KE, KERR
!
       END SUBROUTINE GET_1DGLOBALSLICE_ll
!
!!     ####################################################################
       SUBROUTINE GET_2DGLOBALSLICE_ll( PARRAY, HDIR, KLOC, PGLOBALSLICE, &
                                        KB, KE, KKB, KKE, KERR )
!!     ####################################################################
!
  REAL, DIMENSION(:,:,:) :: PARRAY
  CHARACTER(LEN=1)       :: HDIR
  INTEGER                :: KLOC
  REAL, DIMENSION(:,:)   :: PGLOBALSLICE
  INTEGER, OPTIONAL    :: KB, KE, KKB, KKE, KERR
!
       END SUBROUTINE GET_2DGLOBALSLICE_ll
!
END INTERFACE
!
INTERFACE GET_SLICE_ll
!
!!     #####################################################################
       SUBROUTINE GET_1DSLICE_ll( PARRAY, HDIR, KLOC, PSLICE, KB, KE, KERR )
!!     #####################################################################
!
  REAL, DIMENSION(:,:), TARGET, INTENT(IN) :: PARRAY
  CHARACTER(LEN=1), INTENT(IN)      :: HDIR
  INTEGER, INTENT(IN)               :: KLOC
  REAL, DIMENSION(:), INTENT(OUT)   :: PSLICE
  INTEGER, OPTIONAL    :: KB, KE, KERR
!
       END SUBROUTINE GET_1DSLICE_ll
!
!!     ########################################################
       SUBROUTINE GET_2DSLICE_ll( PARRAY, HDIR, KLOC, PSLICE, &
                                  KB, KE, KKB, KKE, KERR )
!!     ########################################################
!
  REAL, DIMENSION(:,:,:) :: PARRAY
  CHARACTER(LEN=1)       :: HDIR
  INTEGER                :: KLOC
  REAL, DIMENSION(:,:)   :: PSLICE
  INTEGER, OPTIONAL    :: KB, KE, KKB, KKE, KERR
!
       END SUBROUTINE GET_2DSLICE_ll
!
END INTERFACE
!
END MODULE MODI_GET_ll
