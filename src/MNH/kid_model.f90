!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 adiab 2006/05/18 13:07:25
!-----------------------------------------------------------------
!####################
MODULE MODI_KID_MODEL
!####################
!
INTERFACE
!
RECURSIVE SUBROUTINE KID_MODEL(KMODEL,KTEMP_MODEL,OEXIT)
!
INTEGER, INTENT(IN) :: KMODEL           ! Model number
INTEGER, INTENT(IN) :: KTEMP_MODEL      ! temporal loop index for model KMODEL
LOGICAL, INTENT(INOUT) :: OEXIT         ! switch for the last time step
!
END SUBROUTINE KID_MODEL
!
END INTERFACE
!
END MODULE MODI_KID_MODEL
!
!#######################################################
RECURSIVE SUBROUTINE KID_MODEL(KMODEL,KTEMP_MODEL,OEXIT)
!#######################################################
!
!!****  *KID_MODEL * -recursive call of the kid models of model KMODEL
!!
!!    PURPOSE
!!    -------
!! 
!!      This subroutine call for every kid model of the model KMODEL, the  
!!    time step computation, performed by MODEL$n and the call to kids of
!!    the kid model (KID_MODEL).
!!    
!!    
!!    
!!**  METHOD
!!    ------
!!      The recursivity, allowed by Fortran 90, is used to explore all the 
!!    possible configurations of the different models, regarding the nesting.
!!    The time step computation is called for the kid models JKID and the call
!!    to their own kid models is performed by the KID_MODEL call.
!!
!!    EXTERNAL
!!    --------
!!
!!       subroutine MODEL: choose the right MODELn to be called
!! 
!!       subroutine KID_MODEL: recursive function which calls the kid models
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    MODD_CONF       :    NMODEL
!!    MODD_NESTING    :    NDAD, NDTRATIO
!!                   
!!
!!    REFERENCE
!!    ---------
!!
!!    NONE
!!
!!    AUTHOR
!!    ------
!!
!!       J. STEIN  * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    09/04/99 
!!
!!
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CONF
USE MODD_NESTING
!
USE MODI_MODEL_n
USE MODE_MODELN_HANDLER
!
USE MODE_ll
!
!
!*       0.1    declarations of arguments
!
INTEGER, INTENT(IN) :: KMODEL           ! Model number
INTEGER, INTENT(IN) :: KTEMP_MODEL      ! temporal loop index for model KMODEL
LOGICAL, INTENT(INOUT) :: OEXIT         ! switch for the last time step
!
!*       0.2    declarations of local variables
!
!
INTEGER :: JKID          ! loop index on the kid models
INTEGER :: JTEMP_KID     ! nested temporal loop for the kid model
INTEGER :: ITEMP_LOOP    ! number of the temporal iteration for the kid model
LOGICAL :: GEXIT         ! return value of the EXIT signal from MODEL
INTEGER :: IINFO_ll      ! return code of // routines
!
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION
!              --------------
!
DO JKID=KMODEL+1,NMODEL   
  !
  IF ( NDAD(JKID)==KMODEL ) THEN
    !
    DO JTEMP_KID=1,NDTRATIO(JKID)
      !
      ! compute the time step number for this model
      ITEMP_LOOP=(KTEMP_MODEL-2)*NDTRATIO(JKID) + JTEMP_KID +1
      !
      ! switch for the last time step
      GEXIT= OEXIT .AND. (JTEMP_KID==NDTRATIO(JKID))
      !
      ! call the model$n corresponding to JKID
      CALL GO_TOMODEL_ll(JKID,IINFO_ll)
      CALL GOTO_MODEL(JKID)
      CALL MODEL_n(ITEMP_LOOP,GEXIT)
      !
      ! call to the kid models of model JKID
      CALL KID_MODEL(JKID,ITEMP_LOOP,GEXIT)
      !
    END DO
    !
  END IF
  !
END DO
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE KID_MODEL  
  
