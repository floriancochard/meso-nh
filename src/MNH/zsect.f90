!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!     #################
      MODULE MODI_ZSECT
!     #################
INTERFACE
      FUNCTION ZSECT(PZSECT,PZ,PVAR) RESULT(PHORSECT)
!
REAL,                     INTENT(IN)  :: PZSECT  ! altitude of the cross section
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR    ! physical fields
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZ      ! discretization in z of the field
REAL                                  :: PHORSECT! mean on the cross section
!                                                ! on the INNER points
END FUNCTION ZSECT
END INTERFACE
END MODULE MODI_ZSECT
!
!     ###############################################
      FUNCTION ZSECT(PZSECT,PZ,PVAR) RESULT(PHORSECT)
!     ###############################################
!
!!****  *ZSECT* - Mean on horizontal planes
!!
!!    PURPOSE
!!    -------
!!    This routine computes the mean of a 3D field on the cross section
!!    at constant altitude PZSECT for the inner horizontal points.
!!
!!    CAUTION:
!!    The level numbers must increase from bottom to top.
!!
!!**  METHOD
!!    ------
!     The vertical interpolations on PVAR are linear with z.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!         JPHEXT
!!
!!    REFERENCE
!!    ---------
!!
!!    Book 2
!!
!!    AUTHOR
!!    ------
!!	
!     V.Masson  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/12/94
!!      J. Escobar  24/03/2012 modif for reprod sum
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!JUAN REALZ
USE MODE_ll
USE MODE_REPRO_SUM
!JUAN REALZ
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL,                     INTENT(IN)  :: PZSECT  ! altitude of the cross section
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PVAR    ! physical fields
REAL,   DIMENSION(:,:,:), INTENT(IN)  :: PZ      ! discretization in z of the field
REAL                                  :: PHORSECT! mean on the cross section
!                                                ! on the INNER points
!
!*       0.2   Declaration of local variables
!              ------------------------------
INTEGER                                                  :: JI,JJ
INTEGER                                                  :: IIB,IIE,IJB,IJE
REAL,   DIMENSION(SIZE(PVAR,1),SIZE(PVAR,2))             :: ZVARZSECT
INTEGER,DIMENSION(SIZE(PVAR,1),SIZE(PVAR,2))             :: ILEVEL
LOGICAL,DIMENSION(SIZE(PVAR,1),SIZE(PVAR,2),SIZE(PVAR,3)):: GLEVEL
LOGICAL,DIMENSION(SIZE(PVAR,1),SIZE(PVAR,2))             :: GMASK2D,& !inner points.
                           GMASK ! defines the inner points where the cross section
!                                ! is over the relief (points of interest).
!JUAN REALZ
REAL                                                     :: ZCOUNT
INTEGER                                                  :: IINFO_ll
!JUAN REALZ
!
!-------------------------------------------------------------------------------
!
!*       1.    Determination of the inner points of the horizontal domain
!              ----------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
GMASK2D=RESHAPE((/ ((JI>=IIB.AND.JI<=IIE.AND.JJ>=IJB.AND.JJ<=IJE             &
                    ,JI=1,SIZE(PVAR,1)),JJ=1,SIZE(PVAR,2)) /),               &
                    (/ SIZE(PVAR,1),SIZE(PVAR,2) /) )
!
!-------------------------------------------------------------------------------
!
!*       2.    Determination of the level just under ZSECT
!              --------------------------------------------
!
GLEVEL(:,:,:)=PZ(:,:,:)<=PZSECT
ILEVEL(:,:)=COUNT(GLEVEL,3)
!
!-------------------------------------------------------------------------------
!
!*       3.    Vertical interpolation
!              ----------------------
!
GMASK(:,:)=(ILEVEL(:,:)>0).AND.(ILEVEL(:,:)<SIZE(PVAR,3)).AND.(GMASK2D(:,:))
DO JI=1,SIZE(PVAR,1)
  DO JJ=1,SIZE(PVAR,2)
    IF (GMASK(JI,JJ)) THEN
      ZVARZSECT(JI,JJ)=PVAR(JI,JJ,ILEVEL(JI,JJ)) +              &  
       (PZSECT-PZ(JI,JJ,ILEVEL(JI,JJ)))                         &
      *( (PVAR(JI,JJ,ILEVEL(JI,JJ)+1)-PVAR(JI,JJ,ILEVEL(JI,JJ))) &
      /(PZ(JI,JJ,ILEVEL(JI,JJ)+1)-PZ(JI,JJ,ILEVEL(JI,JJ))) )    

!!$      /(PZ(JI,JJ,ILEVEL(JI,JJ)+1)-PZ(JI,JJ,ILEVEL(JI,JJ)))    &
!!$      *(PVAR(JI,JJ,ILEVEL(JI,JJ)+1)-PVAR(JI,JJ,ILEVEL(JI,JJ))) 
      ELSE
         ZVARZSECT(JI,JJ)= 0.0
    ENDIF
  ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4.    Horizontal mean
!              ---------------
!
   PHORSECT = SUM_DD_R2_ll(ZVARZSECT) ! mask included with 0.0 value

   ZCOUNT   = FLOAT(COUNT(GMASK))
   CALL REDUCESUM_ll(ZCOUNT,IINFO_ll)

IF (ZCOUNT > 0.0 ) THEN
   PHORSECT = PHORSECT / ZCOUNT
ELSE
   PHORSECT=-999.
END IF
!
!-------------------------------------------------------------------------------
!
END FUNCTION ZSECT
