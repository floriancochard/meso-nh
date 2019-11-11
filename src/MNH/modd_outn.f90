!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ################
      MODULE MODD_OUT_n
!     ################
!
!!****  *MODD_OUT$n* - declaration of informations on the instants for the 
!!      outputs
!!
!!    PURPOSE
!!    -------
!       The purpose of this  module is to declare the variables
!     describing the instants for the outputs realized by one nested model.         
!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of Meso-NH documentation (module MODD_OUTn)
!!          
!!    AUTHOR
!!    ------
!!	J.Stein      *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/10/94                      
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX, JPOUTMAX
USE MODD_IO_ll, ONLY:TOUTBAK
IMPLICIT NONE

TYPE OUT_t
!
  INTEGER             :: NBAK_NUMB=0, NOUT_NUMB=0 ! Number of outputs and backups performed by model n
  TYPE(TOUTBAK),DIMENSION(:),POINTER :: TBACKUPN=>NULL(),TOUTPUTN=>NULL()
! Lists of the outputs and backups
!
!
END TYPE OUT_t

TYPE(OUT_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: OUT_MODEL

INTEGER, POINTER :: NBAK_NUMB=>NULL(), NOUT_NUMB=>NULL()
TYPE(TOUTBAK),DIMENSION(:),POINTER :: TBACKUPN=>NULL(), TOUTPUTN=>NULL()

CONTAINS

SUBROUTINE OUT_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Current model is set to model KTO
NBAK_NUMB=>OUT_MODEL(KTO)%NBAK_NUMB
NOUT_NUMB=>OUT_MODEL(KTO)%NOUT_NUMB
TBACKUPN=>OUT_MODEL(KTO)%TBACKUPN
TOUTPUTN=>OUT_MODEL(KTO)%TOUTPUTN

END SUBROUTINE OUT_GOTO_MODEL

END MODULE MODD_OUT_n
