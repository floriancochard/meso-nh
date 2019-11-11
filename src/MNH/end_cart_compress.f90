!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!#############################
 MODULE MODI_END_CART_COMPRESS
!#############################
!
INTERFACE
!
FUNCTION END_CART_COMPRESS(PVARS) RESULT(PCOMPRESS)
!
USE MODD_BUDGET
!
REAL, DIMENSION(:,:,:), INTENT(IN)       :: PVARS     ! Source
REAL, DIMENSION(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX) :: PCOMPRESS ! result
!
END FUNCTION END_CART_COMPRESS
!
END INTERFACE
!
END MODULE MODI_END_CART_COMPRESS
!
!
!
!     ###################################################
      FUNCTION END_CART_COMPRESS(PVARS) RESULT(PCOMPRESS) 
!     ###################################################
!
!!****  *CART_COMPRESS* - function to end the  compress of  the Source in CART case. 
!!                           
!!
!!    PURPOSE
!!    -------
!       This function performs the global sum with the local sums PVARS
!     This reducing  is controlled by 3 
!     logical switches for the budget in I,J and K directions (LBU_ICP,
!     LBU_JCP, LBU_KCP), in the budget box described by the lowest and
!     highest values of the I,J and K indices.  
!
!!**  METHOD
!!    ------
!!      The local sums PVARS are : 
!!   i)  for 0D and 1D(z) arrays 
!!                 transfered in a local variable. This variable  is zero for CBUTYP=SKIP 
!!   ii) for  1D (x or y) and 2D(x or y, z) arrays 
!!       transfered in an array of dimensions the extended subdomain size in this direction
!!                 Zero is put outside the cart region
!!  
!!       Then  the local sum are reduced in an array covering the extended domain. 
!!   The section corresponding to the CART region is kept
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!           LBU_ICP   : switch for compression in I direction
!!           LBU_JCP   : switch for compression in J direction
!!           LBU_KCP   : switch for compression in K direction
!!           NBUIL     : lowest I indice value of the global budget box
!!           NBUJL     : lowest J indice value of the global budget box
!!           NBUSIL     : lowest I indice value of the local budget box
!!           NBUSJL     : lowest J indice value of the local budget box
!!           NBUKL     : lowest K indice value of the budget box
!!           NBUIH     : highest I indice value of the global budget box
!!           NBUJH     : highest J indice value of the global budget box
!!           NBUSIH     : highest I indice value of the local budget box
!!           NBUSJH     : highest J indice value of the local budget box
!!           NBUKH     : highest K indice value of the budget box
!!           NBUIMAX_ll   : dimension along I of the global budget array
!!           NBUJMAX_ll   : dimension along J of the global budget array
!!           NBUKMAX   : dimension along K of the budget array
!!          
!!
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     4/06/99 
!!      J.Escobar       02/10/2015 modif for JPHEXT(JPVEXT) variable 
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
REAL, DIMENSION(:,:,:), INTENT(IN)       :: PVARS     ! Source 
REAL, DIMENSION(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX) :: PCOMPRESS ! result
!
!*       0.2   Declarations of local variables :
! 
! 
REAL, ALLOCATABLE, DIMENSION(:) ::     ZVAR1D 
REAL, ALLOCATABLE, DIMENSION(:,:) ::   ZWORK2D ! Work array
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZVAR2D
INTEGER :: IIU,IJU  ! extended dimension of the subdomain
INTEGER :: IIU_ll,IJU_ll,IIMAX_ll, IJMAX_ll ! extended dimensions and size of the physical global domain
INTEGER :: IINFO_ll  ! return status code of the interface routines
! 
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPRESSIONS IN I,J AND K DIRECTIONS
!	        ------------------------------------
!                                                                                            !                                 
!
!
IF (LBU_ICP.AND.LBU_JCP) THEN
  ALLOCATE(ZVAR1D(SIZE(PVARS,3)))
!*JPC*20081209*BUGPARALLEL
! IF (CBUTYPE=='CART') THEN
  IF (CBUTYPE=='CART'.AND.SIZE(PVARS)/=0) THEN
!*JPC*20081209*BUGPARALLEL
    ZVAR1D(:)=PVARS(1,1,:)
  ELSE ! the processor has a empty intersection with the CART region
    ZVAR1D(:)=0.
  ENDIF 
  CALL REDUCESUM_ll(ZVAR1D,IINFO_ll)
  PCOMPRESS(1,1,:)=ZVAR1D(:)
  DEALLOCATE(ZVAR1D)
!
ELSE IF (LBU_ICP.AND..NOT.LBU_JCP) THEN ! compress along x direction,
                                        ! the result is an array along y,z directions
  CALL GET_DIM_EXT_ll('B',IIU,IJU)
  ALLOCATE(ZVAR2D(1,IJU,SIZE(PVARS,3)))
  CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
  IJU_ll= IJMAX_ll + 2 * JPHEXT
  ALLOCATE (ZWORK2D(IJU_ll,SIZE(PVARS,3)))
  ZVAR2D = 0.
!*JPC*20081209*BUGPARALLEL
! IF (CBUTYPE=='CART') THEN
  IF (CBUTYPE=='CART'.AND.SIZE(PVARS)/=0) THEN
!*JPC*20081209*BUGPARALLEL
    ZVAR2D(1,NBUSJL:NBUSJH,:)=PVARS(1,:,:)
  ENDIF 
  CALL SUM_DIM1_ll(ZVAR2D,ZWORK2D,IINFO_ll)
  PCOMPRESS(1,:,:)=ZWORK2D(NBUJL+JPHEXT:NBUJH+JPHEXT,:)
  DEALLOCATE(ZVAR2D,ZWORK2D)
!
ELSE IF (.NOT.LBU_ICP.AND.LBU_JCP) THEN  ! compress along y direction,
                                         ! the result is an array along x,z directions
  CALL GET_DIM_EXT_ll('B',IIU,IJU)
  ALLOCATE(ZVAR2D(IIU,1,SIZE(PVARS,3))) 
  CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
  IIU_ll= IIMAX_ll + 2 * JPHEXT
  ALLOCATE (ZWORK2D(IIU_ll,SIZE(PVARS,3)))
  ZVAR2D = 0.
!*JPC*20081209*BUGPARALLEL
! IF (CBUTYPE=='CART') THEN
  IF (CBUTYPE=='CART'.AND.SIZE(PVARS)/=0) THEN
!*JPC*20081209*BUGPARALLEL
    ZVAR2D(NBUSIL:NBUSIH,1,:)=PVARS(:,1,:)
  ENDIF
  CALL SUM_DIM1_ll(ZVAR2D,ZWORK2D,IINFO_ll)
  PCOMPRESS(:,1,:)=ZWORK2D(NBUIL+JPHEXT:NBUIH+JPHEXT,:)
  DEALLOCATE(ZVAR2D,ZWORK2D)
!
END IF
!
END FUNCTION END_CART_COMPRESS                           
