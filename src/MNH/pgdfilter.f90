!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_PGDFILTER
!     ####################
INTERFACE
      SUBROUTINE PGDFILTER(PPGDARRAY,KPGDFILTER)
!
REAL,             DIMENSION(:,:), INTENT(INOUT) :: PPGDARRAY  ! pgd field
INTEGER,                          INTENT(IN)    :: KPGDFILTER ! iteration number
!
END SUBROUTINE PGDFILTER
END INTERFACE
END MODULE MODI_PGDFILTER
!
!
!     ##########################################
      SUBROUTINE PGDFILTER(PPGDARRAY,KPGDFILTER)
!     ##########################################
!
!!**** *PGDFILTER* add a Laplacian to filter orographic signal
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!    An iterative method is used, adding each time the discretized Laplacian
!!    to the point value.
!!    Note that only points where the mask is not 0. are
!!    modified, taking into account only such points in the filtering.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson          Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    21/03/96
!!    Modification
!!        25/05/96 (V Masson) remove useless ZMASKIJ
!!        28/11/96 (V Masson) test on point localisation itself
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PARAMETERS
!JUAN REALZ
USE MODE_ll
USE MODE_MPPDB
!JUAN REALZ
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL,             DIMENSION(:,:), INTENT(INOUT) :: PPGDARRAY  ! pgd field
INTEGER,                          INTENT(IN)    :: KPGDFILTER ! iteration number
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL,DIMENSION(0:SIZE(PPGDARRAY,1)+1,0:SIZE(PPGDARRAY,2)+1) :: ZARRAY
REAL,DIMENSION(  SIZE(PPGDARRAY,1)  ,  SIZE(PPGDARRAY,2)  ) :: ZARRAY_ll
!                                                  ! modified pgd field
INTEGER :: JI,JJ,JITER,IIU,IJU
REAL    :: ZK                    ! filter efficiency coefficient (0.=< ZK =<1.)
INTEGER, DIMENSION(2,2) :: ICOEF ! 0 or 1 according to field definition area;
!                                ! 1: used for filtering
!JUAN REALZ
INTEGER                :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll => NULL()  ! list of fields to exchange
INTEGER                :: IIB,IJB,IIE,IJE
!JUAN REALZ
!----------------------------------------------------------------------------
!
!*       1.     Initialisations
!               ---------------
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)

CALL ADD2DFIELD_ll(TZFIELDS_ll,PPGDARRAY )
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL MPPDB_CHECK2D(PPGDARRAY,"PGDFILTER:PPGDARRAY",PRECISION)

ZK=1.
ZARRAY(:,:)=0.
ZARRAY(1:IIU,1:IJU)=PPGDARRAY(:,:)
ZARRAY(0,:)    =2.*ZARRAY(1,:)  -ZARRAY(2,:)
ZARRAY(IIU+1,:)=2.*ZARRAY(IIU,:)-ZARRAY(IIU-1,:)
ZARRAY(:,0)    =2.*ZARRAY(:,1)  -ZARRAY(:,2)
ZARRAY(:,IJU+1)=2.*ZARRAY(:,IJU)-ZARRAY(:,IJU-1)
!
!*       2.     Iterative loop
!               --------------
!
ZARRAY_ll(1:IIU,1:IJU) = ZARRAY(1:IIU,1:IJU)
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZARRAY_ll )
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
ZARRAY(1:IIU,1:IJU) = ZARRAY_ll(1:IIU,1:IJU)
DO JITER=1,KPGDFILTER
  DO JI= IIB-1,IIE+1
    DO JJ= IJB-1,IJE+1
      IF ( ZARRAY(JI,JJ)==XUNDEF ) CYCLE
      ICOEF(:,:)=0
      IF (JI>IIB-1  ) THEN
        IF(PPGDARRAY(JI-1,JJ)/=XUNDEF)   ICOEF(1,1)=1
      END IF
      IF (JI<IIE+1) THEN
        IF(PPGDARRAY(JI+1,JJ)/=XUNDEF)   ICOEF(2,1)=1
      END IF
      IF (JJ>IJB-1  ) THEN
        IF(PPGDARRAY(JI,JJ-1)/=XUNDEF)   ICOEF(1,2)=1
      END IF
      IF (JJ<IJE+1) THEN
        IF(PPGDARRAY(JI,JJ+1)/=XUNDEF)   ICOEF(2,2)=1
      END IF
      IF (ANY( ICOEF == 1 ) )  THEN
        PPGDARRAY(JI,JJ)= ZARRAY(JI,JJ)                      &
           + ZK*0.125                                        &
            * ( ICOEF(1,1)*ZARRAY(JI-1,JJ)                   &
               +ICOEF(2,1)*ZARRAY(JI+1,JJ)                   &
               +ICOEF(1,2)*ZARRAY(JI,JJ-1)                   &
               +ICOEF(2,2)*ZARRAY(JI,JJ+1)                   &
               -SUM(ICOEF)*ZARRAY(JI,JJ)  )
      END IF
    ENDDO
  ENDDO
  ZARRAY_ll(1:IIU,1:IJU)=PPGDARRAY(:,:)
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)  
  ZARRAY(1:IIU,1:IJU) = ZARRAY_ll(1:IIU,1:IJU)
  CALL MPPDB_CHECK2D(PPGDARRAY,"PGDFILTER:LOOP:PPGDARRAY",PRECISION)
ENDDO
!
!-------------------------------------------------------------------------------
!
CALL CLEANLIST_ll(TZFIELDS_ll)
END SUBROUTINE PGDFILTER
