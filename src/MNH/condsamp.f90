!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     ######spl
     MODULE MODI_CONDSAMP
!    ################## 
!
INTERFACE
!
      SUBROUTINE CONDSAMP (PTSTEP, PSFSV, KLUOUT, KVERB)
IMPLICIT NONE
REAL,  INTENT(IN)            ::  PTSTEP  ! Time step
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSFSV ! surface flux of scalars
INTEGER, INTENT(IN)          :: KLUOUT     ! unit for output listing count
INTEGER, INTENT(IN)          :: KVERB      ! verbosity level
!
END SUBROUTINE CONDSAMP
!
END INTERFACE
!
END MODULE MODI_CONDSAMP
!     ######spl
      SUBROUTINE CONDSAMP (PTSTEP, PSFSV, KLUOUT, KVERB)
!     ############################################################
!
!
!
!!****  *PASPOL* -
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to release tracers for conditional
!!       samplings according to Couvreux et al. (2010)
!
!!**  METHOD
!!    ------
!!    
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      F.Couvreux, C.Lac         * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!! --------------------------------------------------------------------------
!       
!!    EXTERNAL
!!    --------
!!
USE MODD_PARAMETERS , ONLY : JPVEXT
USE MODD_NSV        , ONLY : NSV_CSBEG, NSV_CSEND, NSV_CS
USE MODD_CONF_n     , ONLY : LUSERC
USE MODD_FIELD_n    , ONLY : XSVT, XRT,XRSVS
USE MODD_GRID_n     , ONLY : XZHAT
USE MODD_REF_n      , ONLY : XRHODJ
USE MODD_DYN        , ONLY : XTSTEP_MODEL1
USE MODD_CONDSAMP
USE MODE_ll
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL,  INTENT(IN)            ::  PTSTEP  ! Time step
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PSFSV ! surface flux of scalars
INTEGER, INTENT(IN)          :: KLUOUT     ! unit for output listing count
INTEGER, INTENT(IN)          :: KVERB      ! verbosity level
!
!*      0.2    declarations of local variables
!

INTEGER :: IIB,IIE,IJB,IJE, IKB, IKE
INTEGER :: IIU, IJU, IKU                      ! dimensional indexes
INTEGER :: JK,JSV,IBOT,ITOP ! Loop indice
INTEGER :: IINFO_ll       ! return code of parallel routine
REAL, DIMENSION(SIZE(XRT,1),SIZE(XRT,2),SIZE(XRT,3)) :: ZRT
REAL, DIMENSION(SIZE(XSVT,1),SIZE(XSVT,2),SIZE(XSVT,3),SIZE(XSVT,4)) :: ZSVT
!
!--------------------------------------------------------------------------------------
!
!
!*	0. Initialisation
!
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU = SIZE(XRT,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
ZSVT(:,:,:,:) = XSVT(:,:,:,:)
!
!
!
!*	1.  INITIALIZATION OF CONDITIONAL SAMPLING TRACERS
!	    ----------------------------------------------
! on veut initialiser le 1er traceur a la surface tout le temps
! le 2E si cloud 100m en dessous de cloud base
! le 3eme si cloud 100m au dessus de cloud top
!
IBOT=0
ITOP=0

IF (( NSV_CS >= 2) .AND.LUSERC .AND.  MAX_ll(XRT(:,:,:,2),IINFO_ll) > 1.E-6 )  THEN
   ! calcul de la base et du sommet des nuages
   ! on ne considere que l'eau liquide car que pour nuages de couche limite
   DO JK=1,IKE
    ZRT(:,:,:) = SPREAD(XRT(:,:,JK,2),3,IKU)
    IF ((MAX_ll(ZRT(:,:,:),IINFO_ll) > 1.E-6).AND.(IBOT == 0)) IBOT=JK
    IF ( MAX_ll(ZRT(:,:,:),IINFO_ll) > 1.E-6) ITOP=JK
   END DO
   IF (KVERB >= 10) THEN
    WRITE(KLUOUT,'(A)') ' '
    WRITE(KLUOUT,'(A,F7.1)') 'Base nuage  : ',XZHAT(IBOT)
    WRITE(KLUOUT,'(A,F7.1)') 'Sommet nuage: ',XZHAT(ITOP)
    WRITE(KLUOUT,'(A,I3.1)') 'JK Base   : ',IBOT
    WRITE(KLUOUT,'(A,I3.1)') 'JK Sommet : ',ITOP
   END IF
   !
END IF

DO JSV=NSV_CSBEG, NSV_CSEND
 !
 IF (JSV== NSV_CSBEG ) THEN 
  ! emission en surface
  PSFSV(IIB:IIE,IJB:IJE,JSV) = 1.
 ENDIF

 IF ((JSV == NSV_CSBEG + 1 ).AND.(IBOT > 2)) THEN
    ! emission XHEIGHT_BASE(m) below the base on XDEPTH_BASE(m)
    !
    DO JK=1,IKE
     IF ((XZHAT(JK) > XZHAT(IBOT) - XHEIGHT_BASE - XDEPTH_BASE/2. ).AND. &
         (XZHAT(JK) < XZHAT(IBOT) - XHEIGHT_BASE + XDEPTH_BASE/2. )) THEN
         ZSVT(IIB:IIE,IJB:IJE,JK,JSV) =  &
           XSVT(IIB:IIE,IJB:IJE,JK,JSV)+1.  
     END IF
    END DO
 END IF    

 IF ((JSV == NSV_CSBEG + 2 ).AND.(ITOP > 2)) THEN
    ! emission XHEIGHT_TOP(m) above the top on XDEPTH_TOP(m)
    !
    DO JK=1,IKE
     IF ((XZHAT(JK) > XZHAT(ITOP) + XHEIGHT_TOP - XDEPTH_TOP/2. ).AND. &
         (XZHAT(JK) < XZHAT(ITOP) + XHEIGHT_TOP + XDEPTH_TOP/2. )) THEN
         ZSVT(IIB:IIE,IJB:IJE,JK,JSV) = &
           XSVT(IIB:IIE,IJB:IJE,JK,JSV)+1. 
     END IF
    END DO
 END IF    !
!
END DO            
         !
!
! correction d'eventuelle concentration négative
WHERE (ZSVT(:,:,:,NSV_CSBEG:NSV_CSEND) <0.0) &
       ZSVT(:,:,:,NSV_CSBEG:NSV_CSEND)=0.0
!
!
!  2: Radioactive decrease            
!
DO JSV=NSV_CSBEG, NSV_CSEND
   ZSVT(:,:,:,JSV) = ZSVT(:,:,:,JSV) *         &
           EXP(-1.*XTSTEP_MODEL1/XRADIO(JSV-NSV_CSBEG+1))
   XRSVS(:,:,:,JSV) = XRSVS(:,:,:,JSV) + &
        XRHODJ(:,:,:)*(ZSVT(:,:,:,JSV)-XSVT(:,:,:,JSV))/PTSTEP
END DO
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CONDSAMP
