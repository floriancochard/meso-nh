!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/set_mask.f90,v $ $Revision: 1.2.2.1.2.1.18.2 $
! MASDEV4_7 budget 2006/09/08 10:35:15
!-----------------------------------------------------------------
!     ###################
      SUBROUTINE SET_MASK
!     ###################
!
!!****SET_MASK** -routine to define the mask 
!!                           
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to test the occurence or not of the
!     different criteria, used to compute the budgets. It also updates the 
!     number of occurence of the different criteria.
!
!!**  METHOD
!!    ------
!!      According to each criterion associated to one zone, the mask is
!!    set to TRUE at each point where the criterion is confirmed, at each 
!!    time step of the model. Finally, The number of occurence of this criteria is 
!!    increased by 1 and stored in the array XBUSURF. 
!!    Caution : The mask is defined on the inner domain.
!!      
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!         LBU_MASK   : logical array mask defining the zones
!!         NBUTIME    : number of the budget step
!!         XBUSURF    : mask tracer array (surface array) 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (routine BUDGET)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Nicolau       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/02/95
!!      Modification 10/11/97  (P.Jabouille) : computation made only in the inner domain
!!      Modification 18/06/99  (N.Asencio) : //  , computation are performed on the extended
!!                                           domain but logical array mask is initialized
!!                                           to FALSE outside the physical domain
!!                   02/02/2017 (J.Escobar & JPP ) bug for 1 model only <-> remove unneeded FIELD_MODEL%
!---------------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODE_ll
USE MODD_FIELD_n , ONLY : XWT , XRT
!
!
IMPLICIT NONE
!  
!  
!*       0.2   Declarations of local variables :
!
INTEGER                    :: IIB,IJB       ! Lower bounds of the physical
                                            ! sub-domain in x and y directions
INTEGER                    :: IIE,IJE       ! Upper bounds of the physical
                                            ! sub-domain in x and y directions
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
!              ---------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!*	 2.     DEFINITION OF THE MASK
!               ----------------------
!  initialization to FALSE on the extended subdomain
LBU_MASK(:,:,:)=.FALSE.
!
!  computing on the physical subdomain 
!==============================================================================
! Change the following lines to set the criterion for each of the NBUMASK masks
! 
! 1st mask on vertical velocity at level k=10
LBU_MASK(IIB:IIE,IJB:IJE,1)=XWT(IIB:IIE,IJB:IJE,10)>0.
!
!2rd mask on rain mixing ratio at level k=2
IF (NBUMASK>=2) &
  LBU_MASK(IIB:IIE,IJB:IJE,2)=XRT(IIB:IIE,IJB:IJE,2,3)>1.E-8
!
!==============================================================================
!
!*	 3.     INCREASE IN SURFACE ARRAY
!               -------------------------
!
WHERE (LBU_MASK(:,:,:))
  XBUSURF(:,:,:,NBUTIME)=XBUSURF(:,:,:,NBUTIME)+1.
END WHERE
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_MASK
