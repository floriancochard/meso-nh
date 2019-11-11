!MNH_LIC Copyright 1996-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!###################
MODULE MODI_RMS_AT_Z
!###################
INTERFACE
      SUBROUTINE RMS_AT_Z(PF1,PZS1,PZ1,PF2,PZS2,PZ2,HTITLE)
!
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PF1     ! field 1
REAL,DIMENSION(:,:),   INTENT(IN)    :: PZS1    ! zs for field 1
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PZ1     ! z for field 1
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PF2     ! field 2
REAL,DIMENSION(:,:),   INTENT(IN)    :: PZS2    ! zs for field 2
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PZ2     ! z for field 2
CHARACTER(LEN=80),     INTENT(IN)    :: HTITLE  ! title for print message
!
END SUBROUTINE RMS_AT_Z
END INTERFACE
END MODULE MODI_RMS_AT_Z
!
!     ##########################################################
      SUBROUTINE RMS_AT_Z(PF1,PZS1,PZ1,PF2,PZS2,PZ2,HTITLE)
!     ##########################################################
!
!!****  *RMS_AT_Z* - computes RMS on theta between begining
!!                               and end of subroutines of prep_real_case
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/07/96
!!                  25/10/96 (V. Masson) add deallocations
!!                  28/10/96 (V. Masson) bug in sums in bias and rms computations
!!                                       and does not use external points
!!                  26/08/97 (V. Masson) call to new linear vertical
!!                                       interpolation routine
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_FIELD_n
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_REF_n
USE MODD_VER_INTERP_LIN
!
USE MODI_COEF_VER_INTERP_LIN
USE MODI_VER_INTERP_LIN
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PF1     ! field 1
REAL,DIMENSION(:,:),   INTENT(IN)    :: PZS1    ! zs for field 1
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PZ1     ! z for field 1
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PF2     ! field 2
REAL,DIMENSION(:,:),   INTENT(IN)    :: PZS2    ! zs for field 2
REAL,DIMENSION(:,:,:), INTENT(IN)    :: PZ2     ! z for field 2
CHARACTER(LEN=80),     INTENT(IN)    :: HTITLE  ! title for print message
!
!*       0.2   Declaration of local variables
!              ------------------------------
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZF1_z    ! field 1 at z levels
REAL,DIMENSION(:,:,:), ALLOCATABLE:: ZF2_z    ! field 2 at z levels
REAL,DIMENSION(:),     ALLOCATABLE:: ZBIAS    ! bias between begining and end (at pressure levels)
REAL,DIMENSION(:),     ALLOCATABLE:: ZSIG2    ! sig2 between begining and end (at pressure levels)
REAL,DIMENSION(:),     ALLOCATABLE:: ZRMS     ! rms between begining and end (at pressure levels)
REAL,DIMENSION(:),     ALLOCATABLE:: ZZLEVELS ! z levels where RMS are computed
LOGICAL,DIMENSION(:,:,:), ALLOCATABLE :: GMASK! .T. where pressure level is significant
!
INTEGER                           :: IRESP, ILUOUT0
INTEGER                           :: IIU,IJU,IKU
INTEGER                           :: JZ
!
!-------------------------------------------------------------------------------
!
!*     1.     Initializations
!             ---------------
!
IIU=SIZE(PF1,1)
IJU=SIZE(PF1,2)
IKU=SIZE(PF2,3)
!
ALLOCATE(ZZLEVELS(40))
ZZLEVELS(:) = (/ (500.*JZ,JZ=1,40) /)
!
ALLOCATE(ZF1_Z(IIU,IJU,40))
ALLOCATE(ZF2_Z(IIU,IJU,40))
ALLOCATE(ZRMS(40))
ALLOCATE(ZBIAS(40))
ALLOCATE(ZSIG2(40))
!
!-------------------------------------------------------------------------------
!
!*     2.     Interpolations on pressure levels for begining field
!             ----------------------------------------------------
!
CALL COEF_VER_INTERP_LIN(PZ1(:,:,:),SPREAD(SPREAD(ZZLEVELS(1:40),1,IIU),2,IJU))
ZF1_Z(:,:,:)=VER_INTERP_LIN(PF1(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
!
!-------------------------------------------------------------------------------
!
!*     3.     Interpolations on pressure levels for ending field
!             --------------------------------------------------
!
CALL COEF_VER_INTERP_LIN(PZ2(:,:,:),SPREAD(SPREAD(ZZLEVELS(1:40),1,IIU),2,IJU))
ZF2_Z(:,:,:)=VER_INTERP_LIN(PF2(:,:,:),NKLIN(:,:,:),XCOEFLIN(:,:,:))
!
!-------------------------------------------------------------------------------
!
!*     3.     Temperature RMS computation
!             ---------------------------
!
ALLOCATE(GMASK(IIU,IJU,40))
GMASK=   SPREAD(SPREAD(ZZLEVELS(1:40),1,IIU),2,IJU) > SPREAD(PZS1(:,:),3,40) &
   .AND. SPREAD(SPREAD(ZZLEVELS(1:40),1,IIU),2,IJU) > SPREAD(PZS2(:,:),3,40)
GMASK( 1 , : ,:)=.FALSE.
GMASK(IIU, : ,:)=.FALSE.
GMASK( : , 1 ,:)=.FALSE.
GMASK( : ,IJU,:)=.FALSE.
!
DO JZ =1, 40
  IF ( COUNT(GMASK(:,:,JZ)) > 1 ) THEN
    ZBIAS(JZ)= SUM( (ZF1_Z(:,:,JZ)-ZF2_Z(:,:,JZ)          )   ,MASK=GMASK(:,:,JZ)) / COUNT (GMASK(:,:,JZ))
    ZSIG2(JZ)= SUM( (ZF1_Z(:,:,JZ)-ZF2_Z(:,:,JZ)-ZBIAS(JZ))**2,MASK=GMASK(:,:,JZ)) / COUNT (GMASK(:,:,JZ))
    ZRMS(JZ)= SQRT(ZBIAS(JZ)*ZBIAS(JZ)+ZSIG2(JZ))
  ELSE
    ZRMS(JZ)=0.
  END IF
END DO
!
!-------------------------------------------------------------------------------
!
!*     4.     Prints
!             ------
!
ILUOUT0 = TLUOUT0%NLU
!
WRITE(ILUOUT0,*) ''
WRITE(ILUOUT0,*) HTITLE
WRITE(ILUOUT0,*) ''
DO JZ=1,40
  WRITE(ILUOUT0,'(6Hlevel ,F6.0,5H m : ,F9.3)') ZZLEVELS(JZ),ZRMS(JZ)
END DO
WRITE(ILUOUT0,*) ''
!
!-------------------------------------------------------------------------------
!
!*     5.     Deallocations
!             -------------
!
DEALLOCATE(ZZLEVELS)
DEALLOCATE(ZF1_Z)
DEALLOCATE(ZF2_Z)
DEALLOCATE(ZRMS)
DEALLOCATE(ZBIAS)
DEALLOCATE(ZSIG2)
DEALLOCATE(GMASK)
!
!-------------------------------------------------------------------------------
!
!WRITE(ILUOUT0,*) 'Routine RMS_AT_Z completed'
!
END SUBROUTINE RMS_AT_Z
