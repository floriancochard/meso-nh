!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!######################
MODULE MODE_NEIGHBORAVG
!######################
!
CONTAINS
!
SUBROUTINE BLOCKAVG(PMATIN,KDX,KDY,PMATOUT)
!    #####################################################################
!
!!    PURPOSE
!!    -------
!!     AVERAGE PMATIN FIELDS BY BLOCK OF DX TIMES DY
!!
!!    AUTHOR
!!    ------
!!      J. Escobar       *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_PARAMETERS, ONLY : JPHEXT
USE MODD_LBC_N     , ONLY : CLBCX, CLBCY
USE MODE_ll        , ONLY : GET_INDICE_ll, GET_OR_ll 
USE MODE_ll        , ONLY : LNORTH_ll, LSOUTH_ll, LEAST_ll, LWEST_ll  
USE MODI_GET_HALO  , ONLY : GET_HALO

IMPLICIT NONE
!
REAL,    INTENT(IN),DIMENSION(:,:,:)      :: PMATIN
INTEGER, INTENT(IN)                       :: KDX,KDY
REAL,    INTENT(INOUT),DIMENSION(:,:,:)   :: PMATOUT
!
!*       0.1    declarations of local variables
!
REAL, DIMENSION(SIZE(PMATIN,1),SIZE(PMATIN,2),SIZE(PMATIN,3)) :: ZSUM, ZTMP 

INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: IIBDX, IJBDY, IDX, IDY
INTEGER :: IDIM1, IDIM2, IDIM3
INTEGER :: II, IJ, IXOR_ll, IYOR_ll
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION
!              --------------
!
! Init the local grid info 

CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_OR_ll('B', IXOR_ll, IYOR_ll)

! shift of the first local pts to DX/DY grid pts

IDX = MOD(IIB+IXOR_ll-1-JPHEXT-1,KDX) 
IDY = MOD(IJB+IYOR_ll-1-JPHEXT-1,KDY)

IF ( IDX .NE. 0 ) IDX = KDX - IDX
IF ( IDY .NE. 0 ) IDY = KDY - IDY

IIBDX = IIB + IDX
IJBDY = IJB + IDY

! Init Sum

CALL GET_HALO(PMATIN)
ZSUM = 0.0
ZTMP = PMATIN

ZSUM(IIBDX:IIE:KDX,:,:) = PMATIN(IIBDX:IIE:KDX,:,:)

! Do the sum on mod(X,KDX) column

DO II = 1, KDX-1
   IF (LEAST_ll() .AND. CLBCX(2)/='CYCL') ZTMP(IIE+1,:,:) = 0.0
   CALL GET_HALO(ZTMP)
   ZTMP(IIB:IIE,:,:) = ZTMP(IIB+1:IIE+1,:,:)
   ZSUM(IIBDX:IIE:KDX,:,:) = ZSUM(IIBDX:IIE:KDX,:,:) + ZTMP(IIBDX:IIE:KDX,:,:)
END DO

ZTMP =  ZSUM

! DO the sum on mod(Y,KDY) raw

DO IJ = 1, KDY-1
   IF (LNORTH_ll() .AND. CLBCY(2)/='CYCL') ZTMP(:,IJE+1,:) = 0.0
   CALL GET_HALO(ZTMP)
   ZTMP(IIB:IIE,IJB:IJE,:) = ZTMP(IIB:IIE,IJB+1:IJE+1,:)
   ZSUM(IIBDX:IIE:KDX,IJBDY:IJE:KDY,:) = ZSUM(IIBDX:IIE:KDX,IJBDY:IJE:KDY,:) &
                                       + ZTMP(IIBDX:IIE:KDX,IJBDY:IJE:KDY,:)
END DO

ZTMP = 0.0

! Dispatch sum on mod(Y,KDY) raw

ZTMP(IIBDX:IIE:KDX,IJBDY:IJE:KDY,:) = ZSUM(IIBDX:IIE:KDX,IJBDY:IJE:KDY,:)

DO IJ = 1, KDY-1
   IF (LSOUTH_ll() .AND. CLBCY(1)/='CYCL') ZTMP(:,IJB-1,:) = 0.0
   CALL GET_HALO(ZTMP)
   ZTMP(IIB:IIE,IJB:IJE,:) = ZTMP(IIB:IIE,IJB-1:IJE-1,:)
   ZTMP(IIBDX:IIE:KDX,IJBDY:IJE:KDY,:) = ZSUM(IIBDX:IIE:KDX,IJBDY:IJE:KDY,:)
END DO


! Dispatch sum on mod(X,KDX) column

ZSUM(IIBDX:IIE:KDX,:,:) =  ZTMP(IIBDX:IIE:KDX,:,:)

DO II = 1, KDX-1
   IF (LWEST_ll() .AND. CLBCX(1)/='CYCL') ZTMP(IIB-1,:,:) = 0.0
   CALL GET_HALO(ZTMP)
   ZTMP(IIB:IIE,:,:) = ZTMP(IIB-1:IIE-1,:,:)
   ZTMP(IIBDX:IIE:KDX,:,:) = ZSUM(IIBDX:IIE:KDX,:,:)
END DO

CALL GET_HALO(ZTMP)
PMATOUT(:,:,:) = ZTMP(:,:,:) / float(KDX*KDY)

END SUBROUTINE BLOCKAVG

SUBROUTINE MOVINGAVG(PMATIN,KDX,KDY,PMATOUT)
!    #####################################################################
!
!!    PURPOSE
!!    -------
!!     MOVING AVERAGE PMATIN FIELDS OVER OF DX TIMES DY
!!
!!    AUTHOR
!!    ------
!!      J. Escobar       *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_PARAMETERS, ONLY : JPHEXT
USE MODD_LBC_N     , ONLY : CLBCX,CLBCY
USE MODE_ll        , ONLY : GET_INDICE_ll ,  GET_OR_ll 
USE MODE_ll        , ONLY : LNORTH_ll , LSOUTH_ll, LEAST_ll , LWEST_ll  
USE MODI_GET_HALO  , ONLY :  GET_HALO
    
IMPLICIT NONE
!
REAL,    INTENT(IN),DIMENSION(:,:,:)      :: PMATIN
INTEGER, INTENT(IN)                       :: KDX,KDY
REAL,    INTENT(INOUT),DIMENSION(:,:,:)   :: PMATOUT   
!local var
    
REAL, DIMENSION(SIZE(PMATIN,1),SIZE(PMATIN,2),SIZE(PMATIN,3)) :: ZSUMP1 , ZSUMM1 , ZTMP , ZTMP2 
    
INTEGER :: IIB,IIE,IJB,IJE
    
INTEGER ::  II, IJ, IXOR_ll, IYOR_ll , ISX
    
    ! Init the local grid info 
    
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
        
    ! Init Sum
    
ZTMP = PMATIN

    !  shift input tab to -KDX , -KDY

CALL  SHIFTXY(ZTMP,-KDX,-KDY)

ZSUMP1 = 0.0

ISX = +1
DO IJ = 1 , 2*KDY +1 
! Do the sum on X+/-1 
  DO II = 1 , 2*KDX +1 
    ZSUMP1 = ZSUMP1 + ZTMP 
    IF ( II .NE.  2*KDY +1 ) THEN
      CALL  SHIFTXY(ZTMP,ISX,0) 
    END IF
  END DO
! Do the sum on Y+1 
  CALL  SHIFTXY(ZTMP,0,1)
  ISX = - ISX
END DO
    
PMATOUT(:,:,:) = ZSUMP1(:,:,:) / FLOAT((1+2*KDX)*(1+2*KDY))
  
END SUBROUTINE MOVINGAVG

SUBROUTINE SHIFTXY(PTAB,KDX,KDY)

USE MODE_ll        , ONLY : GET_INDICE_ll ,  GET_OR_ll 
USE MODE_ll        , ONLY : LNORTH_ll , LSOUTH_ll, LEAST_ll , LWEST_ll  
USE MODI_GET_HALO  , ONLY : GET_HALO

USE MODD_LBC_N     , ONLY : CLBCX,CLBCY
    
IMPLICIT NONE

!----- Argument 

REAL ,DIMENSION(:,:,:)  :: PTAB
INTEGER, INTENT(IN)         :: KDX,KDY

!----- Local var 

INTEGER :: IIB,IIE,IJB,IJE
INTEGER :: ISDX,ISDY,IADX,IADY
INTEGER :: II, IJ

!------ init 

CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
   
ISDX = SIGN(1,KDX)
ISDY = SIGN(1,KDY)
IADX = ABS(KDX)
IADY = ABS(KDY)
   
!---------- shift in X
DO II = 1 , IADX
  IF (LEAST_ll() .AND. CLBCX(2)/='CYCL'  ) PTAB(IIE+1,:,:) = 0.0
  IF (LWEST_ll() .AND. CLBCX(1)/='CYCL'  ) PTAB(IIB-1,:,:) = 0.0
  CALL GET_HALO(PTAB)
  PTAB(IIB:IIE,:,:) = PTAB(IIB+ISDX:IIE+ISDX,:,:)
END DO

!---------- shitf in Y

DO IJ = 1 , IADY
  IF (LNORTH_ll() .AND. CLBCY(2)/='CYCL') PTAB(:,IJE+1,:) = 0.0
  IF (LSOUTH_ll() .AND. CLBCY(1)/='CYCL') PTAB(:,IJE+1,:) = 0.0
  CALL GET_HALO(PTAB)
  PTAB(IIB:IIE,IJB:IJE,:) = PTAB(IIB:IIE,IJB+ISDY:IJE+ISDY,:)
END DO

CALL GET_HALO(PTAB)

END SUBROUTINE SHIFTXY

END MODULE MODE_NEIGHBORAVG

