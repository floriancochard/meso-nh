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
!!    Authors
!!    -------
!
!     R. Guivarch, D. Lugato    * CERFACS *
!     Ph. Kloos                 * CERFACS - CNRM *
!
!!    Modifications
!!    -------------
!     Original 01/05/98
!     Juan 19/08/2005:  modification INTENT -> INTENT(INOUT)
!-----------------------------------------------------------------

!
!!     #################################################
       SUBROUTINE ADD_FIELD2_ll( TPLIST_ll, TPHALO2_ll )
!!     #################################################
!
  USE MODE_ARGSLIST2_ll, ONLY : E_ADD_FIELD2_ll => ADD_FIELD2_ll
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2_ll, HALO2LIST_ll
!
  TYPE(HALO2LIST_ll), POINTER :: TPLIST_ll
  TYPE(HALO2_ll), TARGET      :: TPHALO2_ll
!
  CALL E_ADD_FIELD2_ll( TPLIST_ll, TPHALO2_ll )
!
       END SUBROUTINE ADD_FIELD2_ll
!
!!     ########################################################
       SUBROUTINE DEL_FIELD2_ll( TPLIST_ll, TPHALO2_ll, KINFO )
!!     ########################################################
!
  USE MODE_ARGSLIST2_ll, ONLY : E_DEL_FIELD2_ll => DEL_FIELD2_ll
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2_ll, HALO2LIST_ll
!
  TYPE(HALO2LIST_ll), POINTER :: TPLIST_ll
  TYPE(HALO2_ll), TARGET      :: TPHALO2_ll
  INTEGER                     :: KINFO
!
  CALL E_DEL_FIELD2_ll( TPLIST_ll, TPHALO2_ll, KINFO )
!
       END SUBROUTINE DEL_FIELD2_ll

!
!!     ################################################
       SUBROUTINE ADD1DFIELD_ll( HDIR, TPLIST, PFIELD )
!!     ################################################
!
  USE MODE_ARGSLIST_ll, ONLY : E_ADD1DFIELD_ll => ADD1DFIELD_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST1D_ll
!
  CHARACTER(LEN=1), INTENT(IN) :: HDIR
  TYPE(LIST1D_ll), POINTER     :: TPLIST
  REAL, DIMENSION(:), TARGET   :: PFIELD
!
  CALL E_ADD1DFIELD_ll( HDIR, TPLIST, PFIELD )
!
       END SUBROUTINE ADD1DFIELD_ll
!
!!     ##########################################
       SUBROUTINE ADD2DFIELD_ll( TPLIST, PFIELD )
!!     ##########################################
!
  USE MODE_ARGSLIST_ll, ONLY : E_ADD2DFIELD_ll => ADD2DFIELD_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER       :: TPLIST
  REAL, DIMENSION(:,:), TARGET :: PFIELD
!
  CALL E_ADD2DFIELD_ll( TPLIST, PFIELD )
!
       END SUBROUTINE ADD2DFIELD_ll
!
!!     ##########################################
       SUBROUTINE ADD3DFIELD_ll( TPLIST, PFIELD )
!!     ##########################################
!
  USE MODE_ARGSLIST_ll, ONLY : E_ADD3DFIELD_ll => ADD3DFIELD_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER         :: TPLIST
  REAL, DIMENSION(:,:,:), TARGET :: PFIELD
!
  CALL E_ADD3DFIELD_ll( TPLIST, PFIELD )
!
       END SUBROUTINE ADD3DFIELD_ll
!
!!     #################################################
       SUBROUTINE DEL1DFIELD_ll( TPLIST, PFIELD, KINFO )
!!     #################################################
!
  USE MODE_ARGSLIST_ll, ONLY : E_DEL1DFIELD_ll => DEL1DFIELD_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST1D_ll
!
  TYPE(LIST1D_ll), POINTER     :: TPLIST
  REAL, DIMENSION(:), TARGET   :: PFIELD
  INTEGER, INTENT(OUT)         :: KINFO
!
  CALL E_DEL1DFIELD_ll( TPLIST, PFIELD, KINFO )
!
       END SUBROUTINE DEL1DFIELD_ll
!
!!     #################################################
       SUBROUTINE DEL2DFIELD_ll( TPLIST, PFIELD, KINFO )
!!     #################################################
!
  USE MODE_ARGSLIST_ll, ONLY : E_DEL2DFIELD_ll => DEL2DFIELD_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER       :: TPLIST
  REAL, DIMENSION(:,:), TARGET :: PFIELD
  INTEGER, INTENT(OUT)         :: KINFO
!
  CALL E_DEL2DFIELD_ll( TPLIST, PFIELD, KINFO )
!
       END SUBROUTINE DEL2DFIELD_ll
!
!!     #################################################
       SUBROUTINE DEL3DFIELD_ll( TPLIST, PFIELD, KINFO )
!!     #################################################
!
  USE MODE_ARGSLIST_ll, ONLY : E_DEL3DFIELD_ll => DEL3DFIELD_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER         :: TPLIST
  REAL, DIMENSION(:,:,:), TARGET :: PFIELD
  INTEGER, INTENT(OUT)           :: KINFO
!
  CALL E_DEL3DFIELD_ll( TPLIST, PFIELD, KINFO )
!
       END SUBROUTINE DEL3DFIELD_ll
!
!!     #################################
       SUBROUTINE CLEANLIST_ll( TPLIST )
!!     #################################
!
  USE MODE_ARGSLIST_ll, ONLY : E_CLEANLIST_ll => CLEANLIST_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER     :: TPLIST
!
  CALL E_CLEANLIST_ll( TPLIST )
!
       END SUBROUTINE CLEANLIST_ll
!
!!     ##################################
       SUBROUTINE END_PARA_ll( KINFO_ll )
!!     ##################################
!
  USE MODE_INIT_ll, ONLY : E_END_PARA_ll => END_PARA_ll
!
  INTEGER, INTENT(OUT) :: KINFO_ll
!
  CALL E_END_PARA_ll( KINFO_ll )
!
       END SUBROUTINE END_PARA_ll
!
!!     ##################################################
       SUBROUTINE  GET_DIM_EXT_ll( HSPLIT, KXDIM, KYDIM )
!!     ##################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_DIM_EXT_ll => GET_DIM_EXT_ll
!
  CHARACTER*1, INTENT(IN) :: HSPLIT
  INTEGER, INTENT(OUT)    :: KXDIM, KYDIM
!
  CALL E_GET_DIM_EXT_ll( HSPLIT, KXDIM, KYDIM )
!
       END SUBROUTINE GET_DIM_EXT_ll
!
!!     ###################################################
       SUBROUTINE  GET_DIM_PHYS_ll( HSPLIT, KXDIM, KYDIM )
!!     ###################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_DIM_PHYS_ll => GET_DIM_PHYS_ll
!
  CHARACTER*1, INTENT(IN) :: HSPLIT
  INTEGER, INTENT(OUT)    :: KXDIM, KYDIM
!
  CALL E_GET_DIM_PHYS_ll( HSPLIT, KXDIM, KYDIM )
!
       END SUBROUTINE GET_DIM_PHYS_ll
!
!!     ##########################################
       SUBROUTINE GET_OR_ll( HSPLIT, KXOR, KYOR )
!!     ##########################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_OR_ll => GET_OR_ll
!
  CHARACTER*1, INTENT(IN) :: HSPLIT
  INTEGER, INTENT(OUT)    :: KXOR, KYOR
!
  CALL E_GET_OR_ll( HSPLIT, KXOR, KYOR )
!
       END SUBROUTINE GET_OR_ll
!
!!     ####################################################
       SUBROUTINE GET_INDICE_ll( KXOR, KYOR, KXEND, KYEND )
!!     ####################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_INDICE_ll => GET_INDICE_ll
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
!
  CALL E_GET_INDICE_ll( KXOR, KYOR, KXEND, KYEND )
!
       END SUBROUTINE GET_INDICE_ll
!
!!     ############################################
       SUBROUTINE GET_GLOBALDIMS_ll( KIMAX, KJMAX )
!!     ############################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_GLOBALDIMS_ll => GET_GLOBALDIMS_ll
!
  INTEGER, INTENT(OUT) :: KIMAX, KJMAX
!
  CALL E_GET_GLOBALDIMS_ll( KIMAX, KJMAX )
!
       END SUBROUTINE GET_GLOBALDIMS_ll
!
!!     ######################################################
       SUBROUTINE GET_PHYSICAL_ll( KXOR, KYOR, KXEND, KYEND )
!!     ######################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_PHYSICAL_ll => GET_PHYSICAL_ll
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
!
  CALL E_GET_PHYSICAL_ll( KXOR, KYOR, KXEND, KYEND )
!
       END SUBROUTINE GET_PHYSICAL_ll
!
!!     ###############################################################
       SUBROUTINE GET_INTERSECTION_ll( KXOR, KYOR, KXEND, KYEND, &
                                       KXORI, KYORI, KXENDI, KYENDI, &
                                       HDOM, KINFO, KIP )
!!     ###############################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_INTERSECTION_ll => GET_INTERSECTION_ll
!
  CHARACTER(LEN=4), INTENT(IN)  :: HDOM
  INTEGER, INTENT(IN)           :: KXOR, KYOR, KXEND, KYEND
  INTEGER, INTENT(OUT)          :: KXORI, KYORI, KXENDI, KYENDI
  INTEGER, INTENT(OUT)          :: KINFO
  INTEGER, INTENT(IN), OPTIONAL :: KIP
!
  CALL E_GET_INTERSECTION_ll( KXOR, KYOR, KXEND, KYEND, &
                              KXORI, KYORI, KXENDI, KYENDI, &
                              HDOM, KINFO, KIP )
!
       END SUBROUTINE GET_INTERSECTION_ll
!
!!  ####################################################################
    SUBROUTINE GET_1DGLOBALSLICE_ll( PARRAY, HDIR, KLOC, PGLOBALSLICE, &
                                     KB, KE, KERR )
!!  ####################################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_1DGLOBALSLICE_ll => GET_1DGLOBALSLICE_ll
!
  REAL, DIMENSION(:,:), TARGET, INTENT(IN) :: PARRAY
  CHARACTER(LEN=1), INTENT(IN)             :: HDIR
  INTEGER, INTENT(IN)                      :: KLOC
  REAL, DIMENSION(:), INTENT(OUT)          :: PGLOBALSLICE
  INTEGER, OPTIONAL                        :: KB, KE, KERR
!
  CALL E_GET_1DGLOBALSLICE_ll( PARRAY, HDIR, KLOC, PGLOBALSLICE, KB, KE, KERR )
!
       END SUBROUTINE GET_1DGLOBALSLICE_ll
!
!!  ####################################################################
    SUBROUTINE GET_2DGLOBALSLICE_ll( PARRAY, HDIR, KLOC, PGLOBALSLICE, &
                                     KB, KE, KKB, KKE, KERR )
!!  ####################################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_2DGLOBALSLICE_ll => GET_2DGLOBALSLICE_ll
!
  REAL, DIMENSION(:,:,:) :: PARRAY
  CHARACTER(LEN=1)       :: HDIR
  INTEGER                :: KLOC
  REAL, DIMENSION(:,:)   :: PGLOBALSLICE
  INTEGER, OPTIONAL      :: KB, KE, KKB, KKE, KERR
!
  CALL E_GET_2DGLOBALSLICE_ll( PARRAY, HDIR, KLOC, PGLOBALSLICE, &
                               KB, KE, KKB, KKE, KERR )
!
       END SUBROUTINE GET_2DGLOBALSLICE_ll
!
!!  #####################################################################
    SUBROUTINE GET_1DSLICE_ll( PARRAY, HDIR, KLOC, PSLICE, KB, KE, KERR )
!!  #####################################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_1DSLICE_ll => GET_1DSLICE_ll
!
  REAL, DIMENSION(:,:), TARGET, INTENT(IN) :: PARRAY
  CHARACTER(LEN=1), INTENT(IN)             :: HDIR
  INTEGER, INTENT(IN)                      :: KLOC
  REAL, DIMENSION(:), INTENT(OUT)          :: PSLICE
  INTEGER, OPTIONAL                        :: KB, KE, KERR
!
  CALL E_GET_1DSLICE_ll( PARRAY, HDIR, KLOC, PSLICE, KB, KE, KERR )
!
       END SUBROUTINE GET_1DSLICE_ll
!
!!  ########################################################
    SUBROUTINE GET_2DSLICE_ll( PARRAY, HDIR, KLOC, PSLICE, &
                               KB, KE, KKB, KKE, KERR )
!!  ########################################################
!
  USE MODE_TOOLS_ll, ONLY : E_GET_2DSLICE_ll => GET_2DSLICE_ll
!
  REAL, DIMENSION(:,:,:) :: PARRAY
  CHARACTER(LEN=1)       :: HDIR
  INTEGER                :: KLOC
  REAL, DIMENSION(:,:)   :: PSLICE
  INTEGER, OPTIONAL      :: KB, KE, KKB, KKE, KERR
!
  CALL E_GET_2DSLICE_ll( PARRAY, HDIR, KLOC, PSLICE, &
                         KB, KE, KKB, KKE, KERR )
!
       END SUBROUTINE GET_2DSLICE_ll
!
!!     ##################################
       SUBROUTINE INI_PARA_ll( KINFO_ll )
!!     ##################################
!
  USE MODE_INIT_ll, ONLY : E_INI_PARA_ll => INI_PARA_ll
!
  INTEGER, INTENT(OUT) :: KINFO_ll
!
  CALL E_INI_PARA_ll( KINFO_ll )
!
       END SUBROUTINE INI_PARA_ll
!
!     #########################################
      SUBROUTINE SET_SPLITTING_ll( HSPLITTING )
!     #########################################
!
  USE MODE_INIT_ll, ONLY : E_SET_SPLITTING_ll=>SET_SPLITTING_ll
!
  CHARACTER(LEN=*) :: HSPLITTING
!
  CALL E_SET_SPLITTING_ll(HSPLITTING)
!
       END SUBROUTINE SET_SPLITTING_ll
!
!     ##################################
      SUBROUTINE SET_LBX_ll( KLBX, KMI )
!     ##################################
!
  USE MODE_INIT_ll, ONLY : E_SET_LBX_ll=>SET_LBX_ll
!
  CHARACTER(LEN=*) :: KLBX
  INTEGER          :: KMI
!
  CALL E_SET_LBX_ll(KLBX, KMI)
!
       END SUBROUTINE SET_LBX_ll
!
!     ##################################
      SUBROUTINE SET_LBY_ll( KLBY, KMI )
!     ###################################
!
  USE MODE_INIT_ll, ONLY : E_SET_LBY_ll=>SET_LBY_ll
!
  CHARACTER(LEN=*) :: KLBY
  INTEGER          :: KMI
!
  CALL E_SET_LBY_ll(KLBY, KMI)
!
       END SUBROUTINE SET_LBY_ll
!
!     ############################################
      SUBROUTINE SET_LBSIZEX_ll( KNBRIM, KRIMTAB )
!     ############################################
!
  USE MODE_LB_ll, ONLY : E_SET_LBSIZEX_ll => SET_LBSIZEX_ll
!
  INTEGER               :: KNBRIM
  INTEGER, DIMENSION(:) :: KRIMTAB
!
  CALL E_SET_LBSIZEX_ll(KNBRIM, KRIMTAB)
!
       END SUBROUTINE SET_LBSIZEX_ll
!
!     ############################################
      SUBROUTINE SET_LBSIZEY_ll( KNBRIM, KRIMTAB )
!     ############################################
!
  USE MODE_LB_ll, ONLY : E_SET_LBSIZEY_ll => SET_LBSIZEY_ll
!
  INTEGER               :: KNBRIM
  INTEGER, DIMENSION(:) :: KRIMTAB
!
  CALL E_SET_LBSIZEY_ll(KNBRIM, KRIMTAB)
!
       END SUBROUTINE SET_LBSIZEY_ll
!
!     ###################################
      SUBROUTINE SET_DIM_ll( KX, KY, KZ )
!     ###################################
!
  USE MODE_INIT_ll, ONLY : E_SET_DIM_ll=>SET_DIM_ll
!
  INTEGER :: KX,KY,KZ
!
  CALL E_SET_DIM_ll(KX, KY, KZ)
!
       END SUBROUTINE SET_DIM_ll
!
!     #######################################################
      SUBROUTINE SET_JP_ll( KMODELMAX, KHEXT, KVEXT, KPHALO )
!     #######################################################
!
  USE MODE_INIT_ll, ONLY : E_SET_JP_ll=>SET_JP_ll
  IMPLICIT NONE 
!
  INTEGER :: KMODELMAX, KHEXT, KVEXT, KPHALO
!
  CALL E_SET_JP_ll(KMODELMAX, KHEXT, KVEXT, KPHALO)
!
       END SUBROUTINE SET_JP_ll
!
!     ########################################
      SUBROUTINE SET_XRATIO_ll( KXRATIO, KMI )
!     ########################################
!
  USE MODE_INIT_ll, ONLY : E_SET_XRATIO_ll=>SET_XRATIO_ll
!
  INTEGER :: KXRATIO, KMI
!
  CALL E_SET_XRATIO_ll(KXRATIO, KMI)
!
       END SUBROUTINE SET_XRATIO_ll
!
!     ########################################
      SUBROUTINE SET_YRATIO_ll( KYRATIO, KMI )
!     ########################################
!
  USE MODE_INIT_ll, ONLY : E_SET_YRATIO_ll=>SET_YRATIO_ll
!
 INTEGER :: KYRATIO, KMI
!
  CALL E_SET_YRATIO_ll(KYRATIO, KMI)
!
       END SUBROUTINE SET_YRATIO_ll
!
!     ##################################
      SUBROUTINE SET_DAD_ll( KDAD, KMI )
!     ##################################
!
  USE MODE_INIT_ll, ONLY : E_SET_DAD_ll=>SET_DAD_ll
!
  INTEGER :: KDAD, KMI
!
  CALL E_SET_DAD_ll(KDAD, KMI)
!
       END SUBROUTINE SET_DAD_ll
!
!     ##################################
      SUBROUTINE SET_XOR_ll( KXOR, KMI )
!     ##################################
!
  USE MODE_INIT_ll, ONLY : E_SET_XOR_ll=>SET_XOR_ll
!
  INTEGER :: KXOR, KMI
!
  CALL E_SET_XOR_ll(KXOR, KMI)
!
       END SUBROUTINE SET_XOR_ll
!
!     ####################################
      SUBROUTINE SET_XEND_ll( KXEND, KMI )
!     ####################################
!
  USE MODE_INIT_ll, ONLY : E_SET_XEND_ll=>SET_XEND_ll
!
  INTEGER :: KXEND, KMI
!
  CALL E_SET_XEND_ll(KXEND, KMI)
! 
       END SUBROUTINE SET_XEND_ll
!
!     ##################################
      SUBROUTINE SET_YOR_ll( KYOR, KMI )
!     ##################################
!
  USE MODE_INIT_ll, ONLY : E_SET_YOR_ll=>SET_YOR_ll
!
  INTEGER :: KYOR, KMI
!
  CALL E_SET_YOR_ll(KYOR, KMI)
!
       END SUBROUTINE SET_YOR_ll
!
!     ####################################
      SUBROUTINE SET_YEND_ll( KYEND, KMI )
!     ####################################
!
  USE MODE_INIT_ll, ONLY : E_SET_YEND_ll=>SET_YEND_ll
!
  INTEGER :: KYEND, KMI
!
  CALL E_SET_YEND_ll(KYEND, KMI)
!
       END SUBROUTINE SET_YEND_ll
!
!     ########################
      SUBROUTINE SET_DAD0_ll()
!     ########################
!
  USE MODE_INIT_ll, ONLY : E_SET_DAD0_ll=>SET_DAD0_ll
!
  CALL E_SET_DAD0_ll()
!
       END SUBROUTINE SET_DAD0_ll
!
!     #######################
      SUBROUTINE INIT_LB_ll()
!     #######################
!
  USE MODE_LB_ll, ONLY : E_INIT_LB_ll => INIT_LB_ll
!
  CALL E_INIT_LB_ll()
!
       END SUBROUTINE INIT_LB_ll
!
!     #######################
      SUBROUTINE SET_LB_FIELD_ll(HLBTYPE, PFIELD, PLBXFIELD, PLBYFIELD, IIB, IJB, IIE, IJE, &
        SHIFTWEST, SHIFTEAST, SHIFTSOUTH, SHIFTNORTH )
!     #######################
!
  USE MODE_LB_ll, ONLY : E_SET_LB_FIELD_ll => SET_LB_FIELD_ll
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
  CALL E_SET_LB_FIELD_ll(HLBTYPE, PFIELD, PLBXFIELD, PLBYFIELD, IIB, IJB, IIE, IJE, &
        SHIFTWEST, SHIFTEAST, SHIFTSOUTH, SHIFTNORTH )
!
       END SUBROUTINE SET_LB_FIELD_ll
!
!!     ###################################
       FUNCTION LNORTH_ll( K, HSPLITTING )
!!     ###################################
!
  USE MODE_TOOLS_ll, ONLY : E_LNORTH_ll => LNORTH_ll
!
  LOGICAL                           :: LNORTH_ll
  INTEGER, INTENT(IN), OPTIONAL     :: K
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING
!
  LNORTH_ll=E_LNORTH_ll( K, HSPLITTING )
!
       END FUNCTION LNORTH_ll
!
!!     ##################################
       FUNCTION LWEST_ll( K, HSPLITTING )
!!     ##################################
!
  USE MODE_TOOLS_ll, ONLY : E_LWEST_ll => LWEST_ll
!
  LOGICAL                           :: LWEST_ll  
  INTEGER, INTENT(IN), OPTIONAL     :: K
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING
!
  LWEST_ll=E_LWEST_ll( K, HSPLITTING )
!
       END FUNCTION LWEST_ll
!
!!     ##################################
       FUNCTION LEAST_ll( K, HSPLITTING )
!!     ##################################
!
  USE MODE_TOOLS_ll, ONLY : E_LEAST_ll => LEAST_ll
!
  LOGICAL                           :: LEAST_ll
  INTEGER, INTENT(IN), OPTIONAL     :: K
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING
!
  LEAST_ll=E_LEAST_ll( K, HSPLITTING )
!
       END FUNCTION LEAST_ll
!
!!     ###########################################
       FUNCTION LSOUTH_ll( K, HSPLITTING )
!!     ###########################################
!
  USE MODE_TOOLS_ll, ONLY : E_LSOUTH_ll => LSOUTH_ll
!
  LOGICAL                           :: LSOUTH_ll
  INTEGER, INTENT(IN), OPTIONAL     :: K
  CHARACTER*1, INTENT(IN), OPTIONAL :: HSPLITTING
!
  LSOUTH_ll=E_LSOUTH_ll( K, HSPLITTING )
!
       END FUNCTION LSOUTH_ll
!
!     ###############################################
      SUBROUTINE GET_MODEL_NUMBER_ll( KMODEL_NUMBER )
!     ###############################################
!
  USE MODE_NEST_ll, ONLY : E_GET_MODEL_NUMBER_ll => GET_MODEL_NUMBER_ll
!
  INTEGER :: KMODEL_NUMBER
!
  CALL E_GET_MODEL_NUMBER_ll( KMODEL_NUMBER )
!
      END SUBROUTINE GET_MODEL_NUMBER_ll
!
!     ####################################################
      SUBROUTINE GET_CHILD_DIM_ll( KCHILD, KX, KY, KINFO )
!     ####################################################
!
  USE MODE_NEST_ll, ONLY : E_GET_CHILD_DIM_ll => GET_CHILD_DIM_ll
!
  INTEGER, INTENT(IN) :: KCHILD
!
  INTEGER, INTENT(OUT) :: KX, KY
!
  INTEGER, INTENT(OUT) :: KINFO
!
  CALL E_GET_CHILD_DIM_ll( KCHILD, KX, KY, KINFO )
!
       END SUBROUTINE GET_CHILD_DIM_ll
!
!     ###################################################################
      SUBROUTINE GET_FEEDBACK_COORD_ll( KXOR, KYOR, KXEND, KYEND, KINFO )
!     ###################################################################
!
  USE MODE_NEST_ll, ONLY : E_GET_FEEDBACK_COORD_ll => GET_FEEDBACK_COORD_ll
!
  INTEGER, INTENT(OUT) :: KXOR, KYOR, KXEND, KYEND
  INTEGER, INTENT(OUT) :: KINFO
!
  CALL E_GET_FEEDBACK_COORD_ll( KXOR, KYOR, KXEND, KYEND, KINFO )
!
       END SUBROUTINE GET_FEEDBACK_COORD_ll
!
!     #############################################################
      SUBROUTINE SET_LS2DFIELD_1WAY_ll( P2DFIELD, PTFIELD, KMODEL )
!     #############################################################
!
  USE MODE_LS_ll, ONLY : E_SET_LS2DFIELD_1WAY_ll => SET_LS2DFIELD_1WAY_ll
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
  INTEGER, INTENT(IN)                      :: KMODEL
!
  CALL E_SET_LS2DFIELD_1WAY_ll( P2DFIELD, PTFIELD, KMODEL )
!
       END SUBROUTINE SET_LS2DFIELD_1WAY_ll
!
!     #############################################################
      SUBROUTINE SET_LS3DFIELD_1WAY_ll( P3DFIELD, PTFIELD, KMODEL )
!     #############################################################
!
  USE MODE_LS_ll, ONLY : E_SET_LS3DFIELD_1WAY_ll => SET_LS3DFIELD_1WAY_ll
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
  INTEGER, INTENT(IN)                        :: KMODEL
!
  CALL E_SET_LS3DFIELD_1WAY_ll( P3DFIELD, PTFIELD, KMODEL )
!
       END SUBROUTINE SET_LS3DFIELD_1WAY_ll
!
!     ##################################
      SUBROUTINE UNSET_LSFIELD_1WAY_ll()
!     ##################################
!
  USE MODE_LS_ll, ONLY : E_UNSET_LSFIELD_1WAY_ll => UNSET_LSFIELD_1WAY_ll
!
  CALL E_UNSET_LSFIELD_1WAY_ll()
!
       END SUBROUTINE UNSET_LSFIELD_1WAY_ll
!
!     #####################################################
      SUBROUTINE SET_LS2DFIELD_2WAY_ll( P2DFIELD, PTFIELD )
!     #####################################################
!
  USE MODE_LS_ll, ONLY : E_SET_LS2DFIELD_2WAY_ll => SET_LS2DFIELD_2WAY_ll
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
!
  CALL E_SET_LS2DFIELD_2WAY_ll( P2DFIELD, PTFIELD )
!
       END SUBROUTINE SET_LS2DFIELD_2WAY_ll
!
!     #####################################################
      SUBROUTINE SET_LS3DFIELD_2WAY_ll( P3DFIELD, PTFIELD )
!     #####################################################
!
  USE MODE_LS_ll, ONLY : E_SET_LS3DFIELD_2WAY_ll => SET_LS3DFIELD_2WAY_ll
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
!
  CALL E_SET_LS3DFIELD_2WAY_ll( P3DFIELD, PTFIELD )
!
       END SUBROUTINE SET_LS3DFIELD_2WAY_ll
!
!     ##########################################
      SUBROUTINE UNSET_LSFIELD_2WAY_ll( KMODEL )
!     ##########################################
!
  USE MODE_LS_ll, ONLY : E_UNSET_LSFIELD_2WAY_ll => UNSET_LSFIELD_2WAY_ll
!
  INTEGER, INTENT(IN) :: KMODEL
!
  CALL E_UNSET_LSFIELD_2WAY_ll( KMODEL )
!
       END SUBROUTINE UNSET_LSFIELD_2WAY_ll
!
!     #########################################
      SUBROUTINE LS_FORCING_ll( KCHILD, KINFO, OEXTRAPOL, OCYCLIC_EXTRAPOL )
!     #########################################
!
  USE MODE_LS_ll, ONLY : E_LS_FORCING_ll => LS_FORCING_ll
!
  INTEGER, INTENT(IN)  :: KCHILD
  INTEGER, INTENT(OUT) :: KINFO
  LOGICAL, OPTIONAL, INTENT(IN) :: OEXTRAPOL
  LOGICAL, OPTIONAL, INTENT(IN) :: OCYCLIC_EXTRAPOL
!
  IF ( PRESENT(OEXTRAPOL) .AND. PRESENT(OCYCLIC_EXTRAPOL) ) THEN
    CALL E_LS_FORCING_ll(  KCHILD, KINFO, OEXTRAPOL, OCYCLIC_EXTRAPOL )
  ELSEIF ( PRESENT(OEXTRAPOL) ) THEN
    CALL E_LS_FORCING_ll(  KCHILD, KINFO, OEXTRAPOL )
  ELSE
    CALL E_LS_FORCING_ll(  KCHILD, KINFO )
  ENDIF
!
       END SUBROUTINE LS_FORCING_ll
!
!     ###################################
      SUBROUTINE LS_FEEDBACK_ll( KINFO  )
!     ###################################
!
  USE MODE_LS_ll, ONLY : E_LS_FEEDBACK_ll => LS_FEEDBACK_ll
!
  INTEGER, INTENT(OUT) :: KINFO
!
  CALL E_LS_FEEDBACK_ll( KINFO )
!
       END SUBROUTINE LS_FEEDBACK_ll
!
!     ##############################################################
      SUBROUTINE SET_LB2DFIELD_ll( P2DFIELD, PTFIELD, KFINELBSIZE, &
                                   HSIDE, KMODEL )
!     ##############################################################
!
  USE MODE_LB_ll, ONLY : E_SET_LB2DFIELD_ll => SET_LB2DFIELD_ll
!
  REAL, DIMENSION(:,:), INTENT(IN), TARGET :: P2DFIELD, PTFIELD
  INTEGER, INTENT(IN)                      :: KFINELBSIZE, KMODEL
  CHARACTER(LEN=*), INTENT(IN)             :: HSIDE
!
  CALL E_SET_LB2DFIELD_ll( P2DFIELD, PTFIELD, KFINELBSIZE, &
                           HSIDE, KMODEL ) 
!
       END SUBROUTINE SET_LB2DFIELD_ll
!
!     ##############################################################
      SUBROUTINE SET_LB3DFIELD_ll( P3DFIELD, PTFIELD, KFINELBSIZE, &
                                   HSIDE, KMODEL )
!     ##############################################################
!
  USE MODE_LB_ll, ONLY : E_SET_LB3DFIELD_ll => SET_LB3DFIELD_ll
!
  REAL, DIMENSION(:,:,:), INTENT(IN), TARGET :: P3DFIELD, PTFIELD
  INTEGER, INTENT(IN)                        :: KFINELBSIZE, KMODEL
  CHARACTER(LEN=*), INTENT(IN)               :: HSIDE
!
  CALL E_SET_LB3DFIELD_ll( P3DFIELD, PTFIELD, KFINELBSIZE, &
                           HSIDE, KMODEL )
!
       END SUBROUTINE SET_LB3DFIELD_ll
!
!     #############################
      SUBROUTINE UNSET_LBFIELD_ll()
!     #############################
!
  USE MODE_LB_ll, ONLY : E_UNSET_LBFIELD_ll => UNSET_LBFIELD_ll
!
  CALL E_UNSET_LBFIELD_ll()
!
       END SUBROUTINE UNSET_LBFIELD_ll
!
!     #########################################
      SUBROUTINE LB_FORCING_ll( KCHILD, KINFO )
!     #########################################
!
  USE MODE_LB_ll, ONLY : E_LB_FORCING_ll => LB_FORCING_ll
!
  INTEGER, INTENT(IN)  :: KCHILD
  INTEGER, INTENT(OUT) :: KINFO
!
   CALL E_LB_FORCING_ll( KCHILD, KINFO )
!
       END SUBROUTINE LB_FORCING_ll
!
!!     ###########################################################
       FUNCTION LBFINE2COARSE( KRATIO, KLBSIZE ) RESULT( KCOARSE )
!!     ###########################################################
!
  USE MODE_NEST_ll, ONLY : E_LBFINE2COARSE => LBFINE2COARSE
  IMPLICIT NONE
!
  INTEGER :: KCOARSE
  INTEGER :: KRATIO, KLBSIZE
!
  KCOARSE = E_LBFINE2COARSE( KRATIO, KLBSIZE )
!
       END FUNCTION LBFINE2COARSE
!
!     #########################################
      SUBROUTINE GO_TOMODEL_ll( KMODEL, KINFO )
!     #########################################
!
  USE MODE_NEST_ll, ONLY : E_GO_TOMODEL_ll => GO_TOMODEL_ll
  IMPLICIT NONE
!
  INTEGER :: KMODEL, KINFO
!
  CALL E_GO_TOMODEL_ll( KMODEL, KINFO )
!
       END SUBROUTINE GO_TOMODEL_ll
!
!!     ########################################################
       SUBROUTINE REMAP_2WAY_X_ll( PFIELDIN, PFIELDOUT, KINFO )
!!     ########################################################
!
  USE MODE_EXCHANGE_ll, ONLY : E_REMAP_2WAY_X_ll => REMAP_2WAY_X_ll
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PFIELDIN
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT
  INTEGER                             :: KINFO
!
  CALL E_REMAP_2WAY_X_ll( PFIELDIN, PFIELDOUT, KINFO )
!
       END SUBROUTINE REMAP_2WAY_X_ll
!
!     ########################################################
      SUBROUTINE REMAP_X_2WAY_ll( PFIELDIN, PFIELDOUT, KINFO )
!     ########################################################
!
  USE MODE_EXCHANGE_ll, ONLY : E_REMAP_X_2WAY_ll => REMAP_X_2WAY_ll
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PFIELDIN
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT
  INTEGER                             :: KINFO
!
  CALL E_REMAP_X_2WAY_ll( PFIELDIN, PFIELDOUT, KINFO )
!
       END SUBROUTINE REMAP_X_2WAY_ll
!
!     #####################################################
      SUBROUTINE REMAP_X_Y_ll( PFIELDIN, PFIELDOUT, KINFO )
!     #####################################################
!
  USE MODE_EXCHANGE_ll, ONLY : E_REMAP_X_Y_ll => REMAP_X_Y_ll
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PFIELDIN
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT
  INTEGER                             :: KINFO
!
  CALL E_REMAP_X_Y_ll( PFIELDIN, PFIELDOUT, KINFO )
!
       END SUBROUTINE REMAP_X_Y_ll
!
!     #####################################################
      SUBROUTINE REMAP_Y_X_ll( PFIELDIN, PFIELDOUT, KINFO )
!     #####################################################
!
  USE MODE_EXCHANGE_ll, ONLY : E_REMAP_Y_X_ll => REMAP_Y_X_ll
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PFIELDIN
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFIELDOUT
  INTEGER                             :: KINFO
!
  CALL E_REMAP_Y_X_ll( PFIELDIN, PFIELDOUT, KINFO )
!
       END SUBROUTINE REMAP_Y_X_ll
!
!!     #######################################################
       FUNCTION EXTRACT_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                            KXEND, KYEND, KZEND )
!!     #######################################################
!
  USE MODE_SUM_ll, ONLY : E_EXTRACT_ll => EXTRACT_ll
!
  REAL, DIMENSION(:,:,:), POINTER    :: EXTRACT_ll
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  INTEGER, INTENT(OUT)               :: KINFO
  INTEGER, OPTIONAL, INTENT(IN)      :: KXOR, KYOR, KZOR, &
                                        KXEND, KYEND, KZEND
!
  EXTRACT_ll => E_EXTRACT_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                              KXEND, KYEND, KZEND )
!
       END FUNCTION EXTRACT_ll
!
!!     ###########################################################
       FUNCTION SUM1D_ll( PFIELD, KDIR, KINFO, KXOR, KYOR, KZOR, &
                           KXEND, KYEND, KZEND )
!!     ###########################################################
!
  USE MODE_SUM_ll, ONLY : E_SUM1D_ll => SUM1D_ll
!
  REAL, DIMENSION(:,:), POINTER      :: SUM1D_ll
  INTEGER, INTENT(IN)                :: KDIR
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  INTEGER, INTENT(OUT)               :: KINFO
  INTEGER, OPTIONAL, INTENT(IN)      :: KXOR, KYOR, KZOR, &
                                        KXEND, KYEND, KZEND
!
  SUM1D_ll => E_SUM1D_ll( PFIELD, KDIR, KINFO, KXOR, KYOR, KZOR, &
                        KXEND, KYEND, KZEND )
!
       END FUNCTION SUM1D_ll
!
!!     ###################################################################
       FUNCTION SUM2D_ll( PFIELD, KDIR1, KDIR2, KINFO, KXOR, KYOR, KZOR, &
                          KXEND, KYEND, KZEND )
!!     ###################################################################
!
  USE MODE_SUM_ll, ONLY : E_SUM2D_ll => SUM2D_ll
!
  REAL, DIMENSION(:), POINTER        :: SUM2D_ll
  INTEGER, INTENT(IN)                :: KDIR1, KDIR2
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  INTEGER, INTENT(OUT)               :: KINFO
  INTEGER, OPTIONAL, INTENT(IN)      :: KXOR, KYOR, KZOR, &
                                        KXEND, KYEND, KZEND
!
  SUM2D_ll => E_SUM2D_ll( PFIELD, KDIR1, KDIR2, KINFO, KXOR, KYOR, KZOR, &
                          KXEND, KYEND, KZEND )
!
       END FUNCTION SUM2D_ll
!
!!     #####################################################
       FUNCTION SUM3D_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                          KXEND, KYEND, KZEND )
!!     #####################################################
!
  USE MODE_SUM_ll, ONLY : E_SUM3D_ll => SUM3D_ll
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  INTEGER, INTENT(OUT)               :: KINFO
  INTEGER, OPTIONAL, INTENT(IN)      :: KXOR, KYOR, KZOR, &
                                        KXEND, KYEND, KZEND
!
  SUM3D_ll = E_SUM3D_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                          KXEND, KYEND, KZEND )
!
       END FUNCTION SUM3D_ll
!
!!     #######################################################################
       FUNCTION SUM_1DFIELD_ll( PFIELD, HDIR, KOR, KEND, KERR ) RESULT( ZSUM )
!!     #######################################################################
!
  USE MODE_SUM_ll, ONLY : E_SUM_1DFIELD_ll => SUM_1DFIELD_ll
!
  REAL, DIMENSION(:), INTENT(IN) :: PFIELD
  CHARACTER(LEN=1)               :: HDIR
  INTEGER, OPTIONAL, INTENT(OUT) :: KERR
  INTEGER, OPTIONAL, INTENT(IN)  :: KOR, KEND
  REAL                           :: ZSUM
!
  ZSUM = E_SUM_1DFIELD_ll( PFIELD, HDIR, KOR, KEND, KERR )
!
       END FUNCTION SUM_1DFIELD_ll
!
!!     ###################################################
       FUNCTION MAX_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                             KXEND, KYEND, KZEND )
!!     ###################################################
!
  USE MODE_SUM_ll, ONLY : E_MAX_ll => MAX_ll
!
  REAL                               :: MAX_ll
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  INTEGER, INTENT(OUT)               :: KINFO
  INTEGER, OPTIONAL, INTENT(IN)      :: KXOR, KYOR, KZOR, &
                                        KXEND, KYEND, KZEND
!
  MAX_ll = E_MAX_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                     KXEND, KYEND, KZEND )
!
       END FUNCTION MAX_ll
!
!!     ###################################################
       FUNCTION MIN_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                             KXEND, KYEND, KZEND )
!!     ###################################################
!
  USE MODE_SUM_ll, ONLY : E_MIN_ll => MIN_ll
!
  REAL                               :: MIN_ll
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  INTEGER, INTENT(OUT)               :: KINFO
  INTEGER, OPTIONAL, INTENT(IN)      :: KXOR, KYOR, KZOR, &
                                         KXEND, KYEND, KZEND
!
  MIN_ll = E_MIN_ll( PFIELD, KINFO, KXOR, KYOR, KZOR, &
                             KXEND, KYEND, KZEND )
!
       END FUNCTION MIN_ll
!
!!     ###########################################
       FUNCTION SUMMASK_ll( PFIELD, OMASK, KINFO )
!!     ###########################################
!
  USE MODE_SUM_ll, ONLY : E_SUMMASK_ll => SUMMASK_ll
!
  REAL, DIMENSION(:), POINTER         :: SUMMASK_ll
  REAL, DIMENSION(:,:,:), INTENT(IN)  :: PFIELD
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASK
  INTEGER, INTENT(OUT)                :: KINFO
!
  SUMMASK_ll => E_SUMMASK_ll( PFIELD, OMASK, KINFO )
!
       END FUNCTION SUMMASK_ll
!
!!     ###############################################
       FUNCTION SUMMASKCOMP_ll( PFIELD, OMASK, KINFO )
!!     ###############################################
!
  USE MODE_SUM_ll, ONLY : E_SUMMASKCOMP_ll => SUMMASKCOMP_ll
!
  REAL                                :: SUMMASKCOMP_ll
  REAL, DIMENSION(:,:,:), INTENT(IN)  :: PFIELD
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASK
  INTEGER, INTENT(OUT)                :: KINFO
!
  SUMMASKCOMP_ll = E_SUMMASKCOMP_ll( PFIELD, OMASK, KINFO )
!
       END FUNCTION SUMMASKCOMP_ll
!
!!     #############################################
       SUBROUTINE SUM_DIM1_ll( PFIELD, PRES, KINFO )
!!     #############################################
!
  USE MODE_SUM_ll, ONLY : E_SUM_DIM1_ll => SUM_DIM1_ll
!
  REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD
  REAL, DIMENSION(:,:), INTENT(OUT)  :: PRES
  INTEGER, INTENT(OUT)               :: KINFO
!
  CALL E_SUM_DIM1_ll( PFIELD, PRES, KINFO )
!
       END SUBROUTINE SUM_DIM1_ll
!
!!     #############################################
       SUBROUTINE SUM_DIM2_ll( PFIELD, PRES, KINFO )
!!     #############################################
!
  USE MODE_SUM_ll, ONLY : E_SUM_DIM2_ll => SUM_DIM2_ll
!
  REAL, DIMENSION(:), INTENT(IN)  :: PFIELD
  REAL, DIMENSION(:), INTENT(OUT) :: PRES
  INTEGER, INTENT(OUT)            :: KINFO
!
  CALL E_SUM_DIM2_ll( PFIELD, PRES, KINFO )
!
       END SUBROUTINE SUM_DIM2_ll
!
!!     #######################################################
       FUNCTION GMAXLOC3D_ll( PARRAY, MASK ) RESULT( KMAXLOC )
!!     #######################################################
!
  USE MODE_SUM2_ll, ONLY : E_GMAXLOC3D_ll => GMAXLOC3D_ll
!
  REAL, DIMENSION(:,:,:), INTENT(IN)              :: PARRAY
  LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: MASK
  INTEGER, DIMENSION(3)                           :: KMAXLOC
!
  KMAXLOC = E_GMAXLOC3D_ll( PARRAY, MASK )
!
       END FUNCTION GMAXLOC3D_ll
!
!!     ##############################################################
       FUNCTION GMAXLOC2D_ll( PARRAY, KDIMS, MASK ) RESULT( KMAXLOC )
!!     ##############################################################
!
  USE MODE_SUM2_ll, ONLY : E_GMAXLOC2D_ll => GMAXLOC2D_ll
!
  REAL, DIMENSION(:,:), INTENT(IN)              :: PARRAY
  LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK
  INTEGER, DIMENSION(2)                         :: KMAXLOC
  INTEGER, DIMENSION(2), OPTIONAL               :: KDIMS
!
  KMAXLOC = E_GMAXLOC2D_ll( PARRAY, KDIMS, MASK )
!
       END FUNCTION GMAXLOC2D_ll
!
!!     ##############################################################
       FUNCTION GMAXLOC1D_ll( PARRAY, KDIMS, MASK ) RESULT( KMAXLOC )
!!     ##############################################################
!
  USE MODE_SUM2_ll, ONLY : E_GMAXLOC1D_ll => GMAXLOC1D_ll
!
  REAL, DIMENSION(:), INTENT(IN)              :: PARRAY
  LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: MASK
  INTEGER                                     :: KMAXLOC
  INTEGER,  OPTIONAL                          :: KDIMS
!
  KMAXLOC = E_GMAXLOC1D_ll( PARRAY, KDIMS, MASK )
!
       END FUNCTION GMAXLOC1D_ll
!
!!     #######################################################
       FUNCTION GMINLOC3D_ll( PARRAY, MASK ) RESULT( KMINLOC )
!!     #######################################################
!
  USE MODE_SUM2_ll, ONLY : E_GMINLOC3D_ll => GMINLOC3D_ll
!
  REAL, DIMENSION(:,:,:), INTENT(IN)              :: PARRAY
  LOGICAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL :: MASK
  INTEGER, DIMENSION(3)                           :: KMINLOC
!
  KMINLOC = E_GMINLOC3D_ll( PARRAY, MASK )
!
       END FUNCTION GMINLOC3D_ll
!
!!     ##############################################################
       FUNCTION GMINLOC2D_ll( PARRAY, KDIMS, MASK ) RESULT( KMINLOC )
!!     ##############################################################
!
  USE MODE_SUM2_ll, ONLY : E_GMINLOC2D_ll => GMINLOC2D_ll
!
  REAL, DIMENSION(:,:), INTENT(IN)              :: PARRAY
  LOGICAL, DIMENSION(:,:), INTENT(IN), OPTIONAL :: MASK
  INTEGER, DIMENSION(2)                         :: KMINLOC
  INTEGER, DIMENSION(2), OPTIONAL               :: KDIMS
!
  KMINLOC = E_GMINLOC2D_ll( PARRAY, KDIMS, MASK )
!
       END FUNCTION GMINLOC2D_ll
!
!!     ##############################################################
       FUNCTION GMINLOC1D_ll( PARRAY, KDIMS, MASK ) RESULT( KMINLOC )
!!     ##############################################################
!
  USE MODE_SUM2_ll, ONLY : E_GMINLOC1D_ll => GMINLOC1D_ll
!
  REAL, DIMENSION(:), INTENT(IN)              :: PARRAY
  LOGICAL, DIMENSION(:), INTENT(IN), OPTIONAL :: MASK
  INTEGER                                     :: KMINLOC
  INTEGER,  OPTIONAL                          :: KDIMS
!
  KMINLOC = E_GMINLOC1D_ll( PARRAY, KDIMS, MASK )
!
       END FUNCTION GMINLOC1D_ll
!

!!     ##########################################
       SUBROUTINE REDUCE_SUM_0DD_ll( PRES, KINFO )
!!     ##########################################
!
         USE MODD_REPRO_SUM , ONLY : DOUBLE_DOUBLE
         USE MODE_SUM_ll    , ONLY : E_REDUCE_SUM_0DD_ll => REDUCE_SUM_0DD_ll
         !
         TYPE(DOUBLE_DOUBLE) , INTENT(INOUT)  :: PRES
         INTEGER             , INTENT(OUT)    :: KINFO 
         !
         CALL E_REDUCE_SUM_0DD_ll(PRES, KINFO)
!
       END SUBROUTINE REDUCE_SUM_0DD_ll

!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_0D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_0D_ll => REDUCE_SUM_0D_ll
!
  REAL, INTENT(INOUT)  :: PRES
  INTEGER, INTENT(OUT) :: KINFO
!
  CALL E_REDUCE_SUM_0D_ll(PRES, KINFO)
!
       END SUBROUTINE REDUCE_SUM_0D_ll
!

!!     ##########################################
       SUBROUTINE REDUCE_SUM_1DD_ll( PRES, KINFO )
!!     ##########################################
!
         USE MODD_REPRO_SUM  , ONLY :  DOUBLE_DOUBLE
         USE MODE_SUM_ll     , ONLY : E_REDUCE_SUM_1DD_ll => REDUCE_SUM_1DD_ll
         !
         TYPE(DOUBLE_DOUBLE), DIMENSION(:), INTENT(INOUT) :: PRES
         INTEGER                          , INTENT(OUT)   :: KINFO
         !
         CALL E_REDUCE_SUM_1DD_ll( PRES, KINFO )
!

       END SUBROUTINE REDUCE_SUM_1DD_ll

!!     ##########################################
       SUBROUTINE REDUCE_SUM_1D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_1D_ll => REDUCE_SUM_1D_ll
!
  REAL, DIMENSION(:), INTENT(INOUT) :: PRES
  INTEGER, INTENT(OUT)              :: KINFO
!
  CALL E_REDUCE_SUM_1D_ll( PRES, KINFO )
!
       END SUBROUTINE REDUCE_SUM_1D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_2D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_2D_ll => REDUCE_SUM_2D_ll
!
  REAL, DIMENSION(:,:), INTENT(INOUT) :: PRES
  INTEGER, INTENT(OUT)                :: KINFO
!
  CALL E_REDUCE_SUM_2D_ll( PRES, KINFO )
!
       END SUBROUTINE REDUCE_SUM_2D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_3D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_3D_ll => REDUCE_SUM_3D_ll
!
  REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PRES
  INTEGER, INTENT(OUT)                  :: KINFO
!
  CALL E_REDUCE_SUM_3D_ll(PRES, KINFO)
!
       END SUBROUTINE REDUCE_SUM_3D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I0D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_I0D_ll => REDUCE_SUM_I0D_ll
!
  INTEGER, INTENT(INOUT)  :: PRES
  INTEGER, INTENT(OUT) :: KINFO
!
  CALL E_REDUCE_SUM_I0D_ll(PRES, KINFO)
!
       END SUBROUTINE REDUCE_SUM_I0D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I1D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_I1D_ll => REDUCE_SUM_I1D_ll
!
  INTEGER, DIMENSION(:), INTENT(INOUT) :: PRES
  INTEGER, INTENT(OUT)              :: KINFO
!
  CALL E_REDUCE_SUM_I1D_ll( PRES, KINFO )
!
       END SUBROUTINE REDUCE_SUM_I1D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I2D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_I2D_ll => REDUCE_SUM_I2D_ll
!
  INTEGER, DIMENSION(:,:), INTENT(INOUT) :: PRES
  INTEGER, INTENT(OUT)                :: KINFO
!
  CALL E_REDUCE_SUM_I2D_ll( PRES, KINFO )
!
       END SUBROUTINE REDUCE_SUM_I2D_ll
!
!!     ##########################################
       SUBROUTINE REDUCE_SUM_I3D_ll( PRES, KINFO )
!!     ##########################################
!
  USE MODE_SUM_ll, ONLY : E_REDUCE_SUM_I3D_ll => REDUCE_SUM_I3D_ll
!
  INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: PRES
  INTEGER, INTENT(OUT)                  :: KINFO
!
  CALL E_REDUCE_SUM_I3D_ll(PRES, KINFO)
!
       END SUBROUTINE REDUCE_SUM_I3D_ll
!

!!     ##########################################
       SUBROUTINE UPDATE_HALO_ll( TPLIST, KINFO )
!!     ##########################################
!
  USE MODE_EXCHANGE_ll, ONLY : E_UPDATE_HALO_ll => UPDATE_HALO_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  TYPE(LIST_ll), POINTER :: TPLIST
  INTEGER                :: KINFO
!
  CALL E_UPDATE_HALO_ll( TPLIST, KINFO )
!
       END SUBROUTINE UPDATE_HALO_ll
!
!!     ############################################
       SUBROUTINE UPDATE_1DHALO_ll( TPLIST, KINFO )
!!     ############################################
!
  USE MODE_EXCHANGE_ll, ONLY : E_UPDATE_1DHALO_ll => UPDATE_1DHALO_ll
!
  USE MODD_ARGSLIST_ll, ONLY : LIST1D_ll
!
  TYPE(LIST1D_ll), POINTER :: TPLIST
  INTEGER, INTENT(OUT)     :: KINFO
!
  CALL E_UPDATE_1DHALO_ll( TPLIST, KINFO )
!
       END SUBROUTINE UPDATE_1DHALO_ll
!
!!     ############################################################
       SUBROUTINE UPDATE_BOUNDARIES_ll( HDIRECTION, TPLIST, KINFO )
!!     ############################################################
!
  USE MODE_BOUNDARIES_ll, ONLY : E_UPDATE_BOUNDARIES_ll => UPDATE_BOUNDARIES_ll
  USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
  CHARACTER*2, INTENT(IN) :: HDIRECTION
  TYPE(LIST_ll), POINTER  :: TPLIST
  INTEGER                 :: KINFO
!
  CALL E_UPDATE_BOUNDARIES_ll( HDIRECTION, TPLIST, KINFO )
!
       END SUBROUTINE UPDATE_BOUNDARIES_ll
!
!!     ####################################################################
       SUBROUTINE INIT_HALO2_ll( TPHALO2LIST, KNBVAR, KDIMX, KDIMY, KDIMZ )
!!     ####################################################################
!
  USE MODE_EXCHANGE2_ll, ONLY : E_INIT_HALO2_ll => INIT_HALO2_ll
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll
!
  TYPE(HALO2LIST_ll), POINTER :: TPHALO2LIST
  INTEGER                     :: KNBVAR
  INTEGER                     :: KDIMX, KDIMY, KDIMZ
!
  CALL E_INIT_HALO2_ll(TPHALO2LIST, KNBVAR, KDIMX, KDIMY, KDIMZ)
!
       END SUBROUTINE INIT_HALO2_ll
!
!!     ########################################################
       SUBROUTINE UPDATE_HALO2_ll( TPLIST, TPLISTHALO2, KINFO )
!!     ########################################################
!
  USE MODE_EXCHANGE2_ll, ONLY : E_UPDATE_HALO2_ll => UPDATE_HALO2_ll
!
  USE MODD_ARGSLIST_ll, ONLY : HALO2LIST_ll,LIST_ll
!
  TYPE(LIST_ll), POINTER      :: TPLIST
  TYPE(HALO2LIST_ll), POINTER :: TPLISTHALO2
  INTEGER                     :: KINFO
!
  CALL E_UPDATE_HALO2_ll( TPLIST, TPLISTHALO2, KINFO )
!
       END SUBROUTINE UPDATE_HALO2_ll
!
!     #########################
      SUBROUTINE SCATTER(P1,P2)
!     #########################
!
USE MODE_SCATTER_ll
USE MODD_MPIF
USE MODD_VAR_ll , ONLY : NMNH_COMM_WORLD
!
IMPLICIT NONE
!
!  INCLUDE 'mpif.h'
!
REAL, DIMENSION(:,:),  INTENT(IN)    :: P1
REAL, DIMENSION(:,:),  INTENT(OUT)   :: P2
!
CALL SCATTER_XYFIELD(P1,P2,1,NMNH_COMM_WORLD)
!
END SUBROUTINE SCATTER
