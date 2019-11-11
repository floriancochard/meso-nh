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
!#############################
 MODULE MODI_END_MASK_COMPRESS
!#############################
!
INTERFACE
!
FUNCTION END_MASK_COMPRESS(PVARS) RESULT(PCOMPRESS)
!
USE MODD_BUDGET
!
REAL, DIMENSION(:,:,:), INTENT(IN)                   :: PVARS     ! Source
REAL, DIMENSION(NBUKMAX,NBUWRNB,NBUMASK)             :: PCOMPRESS ! result
!
END FUNCTION END_MASK_COMPRESS
!
END INTERFACE
!
END MODULE MODI_END_MASK_COMPRESS
!
!
!
!     ###################################################
      FUNCTION END_MASK_COMPRESS(PVARS) RESULT(PCOMPRESS) 
!     ###################################################
!
!!****  *END_MASK_COMPRESS* - function to end the  compress of  the Source in MASK cases. 
!!                           
!!
!!    PURPOSE
!!    -------
!       This function performs the global sum with the local sums PVARS
!     This reducing  is controlled by one logical switch for the budget in Z direction 
!     ( LBU_KCP).
!
!!**  METHOD
!!    ------
!!      The local sums PVARS are transfered in a local variable. 
!!      Then  the local sum are reduced . 
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!           LBU_KCP   : switch for compression in K direction
!!           NBUKMAX   : first dimension of the budget array ( number of grid points along
!!                       K direction)
!!           NBUWRNB   : second dimension of the budget array ( number of buffered times)
!!           NBUMASK   : third dimension of the budget array ( number of mask zones)
!!          
!!
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHOR
!!    ------
!!	N. Asencio       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     22/06/99   The budget are reset after the write, so the
!!                              local arrays ZVAR2D and ZVAR3D may not be used
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_PARAMETERS
USE MODE_ll
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments and result :
!
REAL, DIMENSION(:,:,:), INTENT(IN)                   :: PVARS     ! Source 
REAL, DIMENSION(NBUKMAX,NBUWRNB,NBUMASK)             :: PCOMPRESS ! result
!
!*       0.2   Declarations of local variables :
! 
! 
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZVAR2D        ! Work array
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZVAR3D        ! Work array
INTEGER :: IINFO_ll  ! return status code of the interface routines
! 
!
!-------------------------------------------------------------------------------
!
!*       1.     GLOBAL SUMS
!	        ------------------------------------
!                                                                                            !                                 
!
IF (LBU_KCP) THEN
  IF (CBUTYPE=='MASK' ) THEN
    ALLOCATE(ZVAR2D(NBUWRNB,NBUMASK))
    ZVAR2D(:,:)=PVARS(1,:,:)
  ELSE  ! the processor has a empty intersection with the MASK region
    ZVAR2D(:,:)=0.
  ENDIF
  CALL REDUCESUM_ll(ZVAR2D,IINFO_ll)
  PCOMPRESS(1,:,:)=ZVAR2D(:,:)
  DEALLOCATE(ZVAR2D)
!
ELSE IF (.NOT.LBU_KCP) THEN
  ALLOCATE(ZVAR3D(NBUKMAX,NBUWRNB,NBUMASK))
  IF (CBUTYPE=='MASK' ) THEN
    ZVAR3D(:,:,:)=PVARS(:,:,:)
  ELSE ! the processor has a empty intersection with the MASK region
    ZVAR3D(:,:,:)=0.
  ENDIF 
  CALL REDUCESUM_ll(ZVAR3D,IINFO_ll)
  PCOMPRESS(:,:,:)=ZVAR3D(:,:,:)
  DEALLOCATE(ZVAR3D)
!
!
END IF
!
END FUNCTION END_MASK_COMPRESS                           
