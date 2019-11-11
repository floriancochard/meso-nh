!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_nest_pgd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ############################
      SUBROUTINE INIT_HORGRID_ll_n
!     ############################
!
!!****  *INIT_HORGRID_ll_n* - to initialize //
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        11/2004
!!    10/10/2011  J.Escobar call INI_PARAZ_ll
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT, JPMODELMAX
USE MODD_DIM_n,      ONLY : NIMAX, NJMAX
USE MODD_CONF
USE MODD_NESTING,    ONLY : NDAD
!
!JUANZ
USE MODE_SPLITTINGZ_ll
!JUANZ

!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!*       0.2   declarations of local variables
!
INTEGER :: IINFO_ll ! return code of // routines
INTEGER :: IMI ! Current model index
!
!-------------------------------------------------------------------------------
!
IMI = GET_CURRENT_MODEL_INDEX()
CALL SET_DAD0_ll()
CALL SET_DAD_ll(NDAD(IMI),IMI)
CALL SET_DIM_ll(NIMAX, NJMAX, IMI)
CALL SET_LBX_ll('OPEN',IMI)
CALL SET_LBY_ll('OPEN', IMI)
!JUANZ CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_HORGRID_ll_n
