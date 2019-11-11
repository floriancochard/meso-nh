!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_pgd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_LAKEFIT
!     ###################
INTERFACE
      SUBROUTINE LAKEFIT(PPGDZS,PLAKE)
!
REAL,             DIMENSION(:,:), INTENT(INOUT) :: PPGDZS ! orography
REAL,             DIMENSION(:,:), INTENT(IN)    :: PLAKE  ! lake fraction
!
END SUBROUTINE LAKEFIT
END INTERFACE
END MODULE MODI_LAKEFIT
!
!
!     ################################
      SUBROUTINE LAKEFIT(PPGDZS,PLAKE)
!     ################################
!
!!**** *LAKEFIT* builds horizontal lake surfaces
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    An iterative method is used, looking each time at the nearest points
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
!!    Original    27/07/95
!!                19/12/97 (V. Masson) uses lake fraction directly
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL,             DIMENSION(:,:), INTENT(INOUT) :: PPGDZS ! orography
REAL,             DIMENSION(:,:), INTENT(IN)    :: PLAKE  ! lake fraction
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL,DIMENSION(0:SIZE(PPGDZS,1)+1,0:SIZE(PPGDZS,2)+1) :: ZZS,ZZSn! modified 
                                                                 ! orographies
REAL,DIMENSION(0:SIZE(PLAKE,1)+1,0:SIZE(PLAKE,2)+1)   :: ZLAKE
INTEGER :: JI,JJ
REAL    :: ZEPS=1.E-6
!----------------------------------------------------------------------------
!
!*       1.     Initialisations
!               ---------------
!
ZLAKE(:,:)=0.
ZLAKE(1:SIZE(PLAKE,1),1:SIZE(PLAKE,2))=PLAKE(:,:)
ZZS(:,:)=100000.
ZZS(1:SIZE(PPGDZS,1),1:SIZE(PPGDZS,2))=PPGDZS(:,:)
!
ZZSn(:,:)=ZZS(:,:)
!
!*       2.     Iterative loop
!               --------------
!
iter: DO
        DO JI=1,SIZE(PPGDZS,1)
          DO JJ=1,SIZE(PPGDZS,2)
            IF (ZLAKE(JI,JJ)>1.-ZEPS) THEN
              IF (ZLAKE(JI-1,JJ)>1.-ZEPS) ZZS(JI,JJ)=MIN(ZZS(JI,JJ),ZZS(JI-1,JJ))
              IF (ZLAKE(JI+1,JJ)>1.-ZEPS) ZZS(JI,JJ)=MIN(ZZS(JI,JJ),ZZS(JI+1,JJ))
              IF (ZLAKE(JI,JJ-1)>1.-ZEPS) ZZS(JI,JJ)=MIN(ZZS(JI,JJ),ZZS(JI,JJ-1))
              IF (ZLAKE(JI,JJ+1)>1.-ZEPS) ZZS(JI,JJ)=MIN(ZZS(JI,JJ),ZZS(JI,JJ+1))
            ENDIF
          ENDDO
        ENDDO
        IF (ALL(ZZS(:,:)==ZZSn(:,:))) EXIT iter
        ZZSn(:,:)=ZZS(:,:)
      ENDDO iter
!
PPGDZS(:,:)=ZZS(1:SIZE(PPGDZS,1),1:SIZE(PPGDZS,2))
!-------------------------------------------------------------------------------
!
END SUBROUTINE LAKEFIT
