!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 convert 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ############################################################
      MODULE MODI_ZINTER
!     ############################################################
!
INTERFACE
      SUBROUTINE ZINTER(PVMNH,PZGMNH,PVZL,PLZL,KLT,KLN,KKU,KKB,KNP,PUNDEF)
!
REAL,DIMENSION(KLT,KLN,KKU),INTENT(IN) :: PVMNH
REAL,DIMENSION(KLT,KLN,KKU),INTENT(IN) :: PZGMNH 
REAL,DIMENSION(KLT,KLN,KNP),INTENT(OUT):: PVZL 
REAL,DIMENSION(KNP),INTENT(IN)         :: PLZL
REAL,INTENT(IN)         :: PUNDEF
!
INTEGER,INTENT(IN)           :: KLT,KLN,KKU,KKB,KNP
!
END SUBROUTINE ZINTER
END INTERFACE
END MODULE MODI_ZINTER
!
!------------------------------------------------------------------------------
!
!     #############################################################
      SUBROUTINE ZINTER(PVMNH,PZGMNH,PVZL,PLZL,KLN,KLT,KKU,KKB,KNP,PUNDEF)
!     #############################################################
!
!
!!****  *ZINTER * - routine to linearly interpolate
!!
!!     PURPOSE
!!     -------
!    This routine interpolates an input field on Gal-Chen grid, linearly in 
!    another Z-grid (regular or not).
!
!!**   METHOD
!!     ------
!!
!!
!!     EXTERNAL
!!     --------
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!      None
!!
!!     REFERENCE
!!     ---------
!!      Research manual 2 ECMWF forecast model, 1988, Ref M1.6/3
!!      "adiabatic part", Appendix 6 postprocessing
!!      Section 3.  Vertical interpolation, p. A6.5-6
!!      Section 3.4 Extrapolation, pp. A6.6-7
!!
!!     AUTHOR
!!     ------
!!       P. Mascart     * LA *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original       22/04/96
!!       Modification   11/02/99 Chaboureau - some simplifications
!! ----------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!*       0.1  Declaration of arguments 
!
INTEGER,INTENT(IN)           :: KLT  ! number of points in the 1st dimension
INTEGER,INTENT(IN)           :: KLN  ! number of points in the 2nd dimension
INTEGER,INTENT(IN)           :: KKU  ! number of vertical levels 
INTEGER,INTENT(IN)           :: KKB  ! 1st level above ground    
                                     !  KKB: ground ; KKU: top 
REAL,DIMENSION(KLT,KLN,KKU),INTENT(IN) :: PVMNH
REAL,DIMENSION(KLT,KLN,KKU),INTENT(IN) :: PZGMNH 
INTEGER,INTENT(IN)           :: KNP  ! number of new vertical levels
                                     !  1: base ; KNP: top
REAL,DIMENSION(KNP),INTENT(IN)   :: PLZL ! list of the new vertical levels
REAL,DIMENSION(KLT,KLN,KNP),INTENT(OUT):: PVZL ! interpolated output field
REAL,INTENT(IN)         :: PUNDEF  ! undefined value
!
!!  PVMNH  = tableau du champ donne au points masse Meso-NH
!!  PZGMNH = altitude geopotentiel au point masse Meso-NH
!!
!
!*       0.2  Declaration of local variables
!
REAL      :: ZSLOPE
INTEGER   :: JI,JJ,JKZL,JK
INTEGER   :: IKD
!
!------------------------------------------------------------------------------
!
!*       1.   INTERPOLATION
!             -------------
!
OX: DO  JI =1,KLT
  OY:  DO   JJ =1,KLN
    PLEV:  DO   JKZL=1,KNP
      !
      !   i) Zones flagging
      !
      IKD=0
      IF(PLZL(JKZL).GE.PZGMNH(JI,JJ,KKU))             IKD=10*KKU
      DO  JK  =KKU-1,KKB,-1
         IF((PZGMNH(JI,JJ,JK+1).GT.PLZL(JKZL)).AND.   &
           (PLZL(JKZL).GE.PZGMNH(JI,JJ,JK)))          IKD=JK
      END DO
      IF(PLZL(JKZL).LT.PZGMNH(JI,JJ,KKB))             IKD=-10*KKU
      !
      !   ii) Regular points interpolation
      !
      IF(ABS(IKD).NE.(10*KKU)) THEN
         ZSLOPE=(PLZL(JKZL)-PZGMNH(JI,JJ,IKD))       &
              /(PZGMNH(JI,JJ,IKD+1)-PZGMNH(JI,JJ,IKD))
         PVZL(JI,JJ,JKZL)=PVMNH(JI,JJ,IKD)                &
              +ZSLOPE*(PVMNH(JI,JJ,IKD+1)-PVMNH(JI,JJ,IKD))
      ELSE
      !
      !   iii) No extrapolation below the ground and above the top
      !
         PVZL(JI,JJ,JKZL)=PUNDEF
      ENDIF
    END DO PLEV
  END DO OY 
END DO OX
!
!
END SUBROUTINE ZINTER
