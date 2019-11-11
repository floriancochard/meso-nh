!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 budget 2006/05/18 13:07:25
!-----------------------------------------------------------------
!#########################
 MODULE MODI_MASK_COMPRESS
!#########################
!
INTERFACE
!
FUNCTION MASK_COMPRESS(PVARS) RESULT(PCOMPRESS)
!
USE MODD_BUDGET
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS     ! Source 
REAL, DIMENSION(NBUKMAX,NBUMASK)   :: PCOMPRESS ! result
!
END FUNCTION MASK_COMPRESS
!
END INTERFACE
!
END MODULE MODI_MASK_COMPRESS
!     ###############################################
      FUNCTION MASK_COMPRESS(PVARS) RESULT(PCOMPRESS) 
!     ###############################################
!
!!****  *MASK_COMPRESS* - function to compress the Source in MASK case. 
!!                           
!!
!!    PURPOSE
!!    -------
!       This function compresses the Source PVARS of the VARiable
!     VAR whose budget is analysed. Compressions in I and J directions
!     are automatically achieved according an horizontal logical mask 
!     (LBU_MASK). The compression in K direction is controlled by a logical
!     switch (LBU_KCP).
!    
!
!!**  METHOD
!!    ------
!!       The source PVARS is first transfered in a local array whose 
!!    vertical dimension corresponds to the budget one. Then compressions
!!    in I and J are achieved. At last, the compression in K direction is
!!    achieved or not depending on the logical switch LBU_KCP.
!!     
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!           
!!           LBU_KCP    : switch for compression in K direction
!!           NBUKL      : lowest K indice value of the budget box
!!           NBUKH      : highest K indice value of the budget box
!!           NBUKMAX    : dimension along K of the budget tabular
!!           NBUMASK    : number of mask zones where the budget is performed
!!           LBU_MASK   : logical array mask definig each zone where the
!!                        budget is to be computed
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (function CART_COMPRESS)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Nicolau       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/02/95
!!      Modification 10/11/97  (P.Jabouille) : computation made only in the inner domain
!!      Modification 18/06/99  (N. Asencio)  : // add the bounds of the physical
!!                                               sub-domain in x and y directions
!!                                            permutation of the two dimensions of the
!!                                            ZCOMP array : first is K, second is mask number
!!      J.Escobar       02/10/2015 modif for JPHEXT(JPVEXT) variable  
!----------------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_BUDGET
USE MODE_ll
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS     ! Source 
REAL, DIMENSION(NBUKMAX,NBUMASK)   :: PCOMPRESS ! result
!
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION (SIZE(LBU_MASK,1),SIZE(LBU_MASK,2),NBUKH-NBUKL+1) :: ZVARS   ! Work array
REAL, DIMENSION (NBUKH-NBUKL+1,NBUMASK)                           :: ZCOMP   ! Work array
!
INTEGER                    :: IK,IM         ! loop indexes
INTEGER                    :: IIB,IJB       ! Lower bounds of the physical
                                            ! sub-domain in x and y directions
INTEGER                    :: IIE,IJE       ! Upper bounds of the physical
                                            ! sub-domain in x and y directions
! 
!
!-------------------------------------------------------------------------------------
!
!*       1.    COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
!              ---------------------------------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
!*	 2.     SOURCE TRANSFERT IN A LOCAL ARRAY 
!	        ---------------------------------
!
ZVARS=0.
DO IK=1, NBUKH-NBUKL+1
  ZVARS(IIB:IIE,IJB:IJE,IK)=PVARS(IIB:IIE, IJB:IJE , IK+NBUKL+JPVEXT-1)
END DO
!
!-------------------------------------------------------------------------------
!
!*	 3.     COMPRESSIONS IN I AND J DIRECTIONS
!	        ----------------------------------
!
DO IM=1, NBUMASK
  DO IK=1, NBUKH-NBUKL+1
    ZCOMP(IK,IM)=SUM(SUM(ZVARS(IIB:IIE, IJB:IJE, IK),1,LBU_MASK(IIB:IIE, IJB:IJE, IM)),1)
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*	 4.     COMPRESSION IN K DIRECTION
!	        --------------------------
IF (LBU_KCP) THEN
  PCOMPRESS(1,:)=SUM(ZCOMP,1)
!
ELSE
  PCOMPRESS=ZCOMP
!
END IF
!
!--------------------------------------------------------------------------------
!
END FUNCTION MASK_COMPRESS
