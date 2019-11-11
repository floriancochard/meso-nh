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

!     ###################
      MODULE MODI_INIT_ll
!     ###################
!!
INTERFACE 
!
!      #########################################
       SUBROUTINE SET_SPLITTING_ll( HSPLITTING )
!      #########################################
!
  CHARACTER(LEN=*) :: HSPLITTING
!
       END SUBROUTINE SET_SPLITTING_ll
!
!      ##################################
       SUBROUTINE SET_LBX_ll( KLBX, KMI )
!      ##################################
!
  CHARACTER(LEN=*) :: KLBX
  INTEGER :: KMI
!
       END SUBROUTINE SET_LBX_ll
!
!      ##################################
       SUBROUTINE SET_LBY_ll( KLBY, KMI )
!      ##################################
!
  CHARACTER(LEN=*) :: KLBY
  INTEGER :: KMI
!
       END SUBROUTINE SET_LBY_ll
!
!      ############################################
       SUBROUTINE SET_LBSIZEX_ll( KNBRIM, KRIMTAB )
!      ############################################
!
  INTEGER :: KNBRIM
  INTEGER, DIMENSION(:) :: KRIMTAB
!
       END SUBROUTINE SET_LBSIZEX_ll
!
!      ############################################
       SUBROUTINE SET_LBSIZEY_ll( KNBRIM, KRIMTAB )
!      ############################################
!
  INTEGER :: KNBRIM
  INTEGER, DIMENSION(:) :: KRIMTAB
!
       END SUBROUTINE SET_LBSIZEY_ll
!
!      ###################################
       SUBROUTINE SET_DIM_ll( KX, KY, KZ )
!      ###################################
!
 INTEGER :: KX,KY,KZ
!
        END SUBROUTINE SET_DIM_ll
!
!      #######################################################
       SUBROUTINE SET_JP_ll( KMODELMAX, KHEXT, KVEXT, KPHALO )
!      #######################################################
!
  INTEGER :: KMODELMAX, KHEXT, KVEXT, KPHALO
!
       END SUBROUTINE SET_JP_ll
!
!      ########################################
       SUBROUTINE SET_XRATIO_ll( KXRATIO, KMI )
!      ########################################
!
  INTEGER :: KXRATIO, KMI
!
       END SUBROUTINE SET_XRATIO_ll
!
!      ########################################
       SUBROUTINE SET_YRATIO_ll( KYRATIO, KMI )
!      ########################################
!
  INTEGER :: KYRATIO, KMI
!
       END SUBROUTINE SET_YRATIO_ll
!
!      ##################################
       SUBROUTINE SET_DAD_ll( KDAD, KMI )
!      ##################################
!
  INTEGER :: KDAD, KMI
!
       END SUBROUTINE SET_DAD_ll
!
!      ##################################
       SUBROUTINE SET_XOR_ll( KXOR, KMI )
!      ##################################
!
  INTEGER :: KXOR, KMI
!
       END SUBROUTINE SET_XOR_ll
!
!      ####################################
       SUBROUTINE SET_XEND_ll( KXEND, KMI )
!      ####################################
!
  INTEGER :: KXEND, KMI
!
       END SUBROUTINE SET_XEND_ll
!
!      ##################################
       SUBROUTINE SET_YOR_ll( KYOR, KMI )
!      ##################################
!
  INTEGER :: KYOR, KMI
!
       END SUBROUTINE SET_YOR_ll
!
!      ####################################
       SUBROUTINE SET_YEND_ll( KYEND, KMI )
!      ####################################
!
  INTEGER :: KYEND, KMI
!
       END SUBROUTINE SET_YEND_ll
!
!      ########################
       SUBROUTINE SET_DAD0_ll()
!      ########################
!
       END SUBROUTINE SET_DAD0_ll
!
!      #######################
       SUBROUTINE INIT_LB_ll()
!      #######################
!
       END SUBROUTINE INIT_LB_ll
!
!      #######################
       SUBROUTINE SET_LB_FIELD_ll(HLBTYPE, PFIELD, PLBXFIELD, PLBYFIELD, IIB, IJB, IIE, IJE, &
        SHIFTWEST, SHIFTEAST, SHIFTSOUTH, SHIFTNORTH )
  !
  CHARACTER(LEN=*),INTENT(IN) :: HLBTYPE ! LB type : 'LB','LBU'
  REAL, DIMENSION(:,:,:), INTENT(IN)  :: PFIELD      ! field on the whole domain (or subdomain)
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLBXFIELD    ! LB field - X direction
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLBYFIELD    ! LB field - Y direction
  !beginning and end of the local physical subdomain
  INTEGER, INTENT(IN)   :: IIB            ! indice I Beginning in x direction
  INTEGER, INTENT(IN)   :: IJB            ! indice J Beginning in y direction
  INTEGER, INTENT(IN)   :: IIE            ! indice I End       in x direction
  INTEGER, INTENT(IN)   :: IJE            ! indice J End       in y direction
  INTEGER, INTENT(IN)   :: SHIFTWEST, SHIFTEAST, SHIFTSOUTH, SHIFTNORTH ! shifting applied to the indices copied from PFIELD in each direction
                                                                        ! it is used for LBXUM et LBXVM
                                                                        ! I do not know why...
  !
!      #######################
!
       END SUBROUTINE SET_LB_FIELD_ll
!
!      ###################################
        SUBROUTINE INI_PARA_ll( KINFO_ll )
!      ###################################
!
  INTEGER, INTENT(OUT) :: KINFO_ll
!
       END SUBROUTINE INI_PARA_ll
!
!     ###################################
       SUBROUTINE END_PARA_ll( KINFO_ll )
!     ###################################
!
  INTEGER, INTENT(OUT) :: KINFO_ll
!
       END SUBROUTINE END_PARA_ll
!
END INTERFACE
!
END MODULE MODI_INIT_ll
