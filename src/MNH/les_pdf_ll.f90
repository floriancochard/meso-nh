!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 les 2006/10/16 16:28:01
!-----------------------------------------------------------------
!      ###################
MODULE MODI_LES_PDF_ll
!      ###################
!
INTERFACE LES_PDF_ll
!
      SUBROUTINE LES_PDF_ll_2D(PA, OMASK,PMIN,PMAX,PA_PDF_ll)
!
REAL,    DIMENSION(:,:),   INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),   INTENT(IN)  :: OMASK
REAL,                      INTENT(IN)  :: PMIN
REAL,                      INTENT(IN)  :: PMAX
!
REAL,    DIMENSION(:), INTENT(OUT) :: PA_PDF_ll
!
END SUBROUTINE LES_PDF_ll_2D
!
      SUBROUTINE LES_PDF_ll_3D(PA, OMASK,PMIN,PMAX,PA_PDF_ll)

REAL,    DIMENSION(:,:,:),   INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:,:),   INTENT(IN)  :: OMASK
REAL,                        INTENT(IN)  :: PMIN
REAL,                        INTENT(IN)  :: PMAX
!
REAL,    DIMENSION(:,:), INTENT(OUT) :: PA_PDF_ll
!
END SUBROUTINE LES_PDF_ll_3D
!
      SUBROUTINE LES_PDF_ll_3DM(PA, OMASK,PMIN,PMAX,PA_PDF_ll)

REAL,    DIMENSION(:,:,:),   INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),     INTENT(IN)  :: OMASK
REAL,                        INTENT(IN)  :: PMIN
REAL,                        INTENT(IN)  :: PMAX
!
REAL,    DIMENSION(:,:), INTENT(OUT) :: PA_PDF_ll
!
END SUBROUTINE LES_PDF_ll_3DM
END INTERFACE
!
END MODULE MODI_LES_PDF_ll
!
!     #################################################################
      SUBROUTINE LES_PDF_ll_2D(PA, OMASK,PMIN,PMAX,PA_PDF_ll )
!     #################################################################
!
!
!!****  *LES_PDF_ll* computes the average of one field on all processors
!!
!!    PURPOSE
!!    -------
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
!!    F. Couvreux
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         10/02/03
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_LES
USE MODD_PARAMETERS
USE MODE_LL
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:),    INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),    INTENT(IN)  :: OMASK
REAL,                       INTENT(IN)  :: PMIN
REAL,                       INTENT(IN)  :: PMAX
!
REAL,    DIMENSION(:),      INTENT(OUT) :: PA_PDF_ll
!
!       0.2  declaration of local variables
INTEGER                              :: JPDF !pdf counter
REAL,    DIMENSION(NPDF+1)             :: ZLIMITPDF !limits of pdf intervals
REAL,    DIMENSION(NPDF)             :: ZA_SUM 
INTEGER                              :: JI, JJ             ! loop counters
INTEGER                              :: INFO_ll
INTEGER                              :: II, IJ             !number of hor. points
!
!-------------------------------------------------------------------------------
!
!introduction de la distribution equivalente telle que la moyenne de la
!distribution calculee soit egale a la moyenne simulee
IF (NPDF > 0) THEN
DO JPDF=1,NPDF+1
   ZLIMITPDF(JPDF)=PMIN+(PMAX-PMIN)/REAL(NPDF)*(JPDF-1)
ENDDO
II=SIZE(PA,1)
IJ=SIZE(PA,2)
DO JPDF=1,NPDF
  ZA_SUM(JPDF)=0.
 DO JJ=1,IJ
  DO JI=1,II
      IF ( OMASK(JI,JJ).AND.(PA(JI,JJ)/=XUNDEF).AND.&
          (PA(JI,JJ)>=ZLIMITPDF(JPDF)).AND. &
           (PA(JI,JJ)<ZLIMITPDF(JPDF+1)) ) THEN
         ZA_SUM(JPDF)  = ZA_SUM(JPDF) + 1
      END IF
  END DO
 END DO
END DO   
CALL REDUCESUM_ll(ZA_SUM(:),INFO_ll)
PA_PDF_ll(:)=ZA_SUM(:)  

END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_PDF_ll_2D
!
!     #################################################################
      SUBROUTINE LES_PDF_ll_3D(PA, OMASK,PMIN,PMAX,PA_PDF_ll)
!     #################################################################
!
!
!!****  *LES_PDF_ll* computes the average of one field on all processors
!!
!!    PURPOSE
!!    -------
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
!!      F. Couvreux
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         10/02/03
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_LES
USE MODD_PARAMETERS
USE MODE_LL

IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:),   INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:,:),   INTENT(IN)  :: OMASK
REAL,                        INTENT(IN)  :: PMIN
REAL,                        INTENT(IN)  :: PMAX
!
REAL,    DIMENSION(:,:),     INTENT(OUT) :: PA_PDF_ll
!
!
!       0.2  declaration of local variables
!
INTEGER                              :: JPDF !pdf counter
REAL,    DIMENSION(NPDF+1)             :: ZLIMITPDF !limits of pdf intervals
REAL,    DIMENSION(SIZE(PA,3),NPDF)  :: ZA_SUM 
INTEGER                              :: JI, JJ,JK             ! loop counters
INTEGER                              :: INFO_ll
INTEGER                              :: II,IJ,IK!number of hor. and vert. levels
!
!-------------------------------------------------------------------------------
!
!
II=SIZE(PA,1)
IJ=SIZE(PA,2)
IK=SIZE(PA,3)
IF (NPDF > 0) THEN
 DO JPDF=1,NPDF+1
   ZLIMITPDF(JPDF)=PMIN+(PMAX-PMIN)/REAL(NPDF)*(JPDF-1)
 ENDDO
!
DO JPDF=1,NPDF
 DO JK=1,IK
  ZA_SUM(JK,JPDF)=0.
  DO JJ=1,IJ
    DO JI=1,II
        IF ( OMASK(JI,JJ,JK).AND.(PA(JI,JJ,JK)/=XUNDEF).AND.&
	     (PA(JI,JJ,JK)>=ZLIMITPDF(JPDF)).AND.&
              (PA(JI,JJ,JK)<ZLIMITPDF(JPDF+1))) THEN
             ZA_SUM(JK,JPDF)  = ZA_SUM(JK,JPDF) + 1
        END IF
    END DO
  END DO
 END DO
END DO
END IF
CALL REDUCESUM_ll(ZA_SUM(:,:),INFO_ll)
PA_PDF_ll(:,:)=ZA_SUM(:,:)  
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_PDF_ll_3D
!
!
!     #################################################################
      SUBROUTINE LES_PDF_ll_3DM(PA, OMASK,PMIN,PMAX,PA_PDF_ll)
!     #################################################################
!
!
!!****  *LES_PDF_ll* computes the average of one field on all processors
!!
!!    PURPOSE
!!    -------
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
!!      F. Couvreux
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         10/02/03
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_LES
USE MODD_PARAMETERS
USE MODE_LL

IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
REAL,    DIMENSION(:,:,:),   INTENT(IN)  :: PA
LOGICAL, DIMENSION(:,:),     INTENT(IN)  :: OMASK
REAL,                        INTENT(IN)  :: PMIN
REAL,                        INTENT(IN)  :: PMAX
!
REAL,    DIMENSION(:,:),     INTENT(OUT) :: PA_PDF_ll
!
!
!       0.2  declaration of local variables
!
INTEGER                              :: JPDF !pdf counter
REAL,    DIMENSION(NPDF+1)             :: ZLIMITPDF !limits of pdf intervals
REAL,    DIMENSION(SIZE(PA,3),NPDF)  :: ZA_SUM 
INTEGER                              :: JI, JJ,JK             ! loop counters
INTEGER                              :: INFO_ll
INTEGER                              :: II,IJ,IK!number of hor. and vert. levels
!
!-------------------------------------------------------------------------------
!
!
II=SIZE(PA,1)
IJ=SIZE(PA,2)
IK=SIZE(PA,3)
IF (NPDF > 0) THEN
 DO JPDF=1,NPDF+1
   ZLIMITPDF(JPDF)=PMIN+(PMAX-PMIN)/REAL(NPDF)*(JPDF-1)
 ENDDO
!
DO JPDF=1,NPDF
 DO JK=1,IK
  ZA_SUM(JK,JPDF)=0.
  DO JJ=1,IJ
    DO JI=1,II
        IF ( OMASK(JI,JJ).AND.(PA(JI,JJ,JK)/=XUNDEF).AND.&
	     (PA(JI,JJ,JK)>=ZLIMITPDF(JPDF)).AND.&
              (PA(JI,JJ,JK)<ZLIMITPDF(JPDF+1))) THEN
             ZA_SUM(JK,JPDF)  = ZA_SUM(JK,JPDF) + 1
        END IF
    END DO
  END DO
 END DO
END DO
END IF
CALL REDUCESUM_ll(ZA_SUM(:,:),INFO_ll)
PA_PDF_ll(:,:)=ZA_SUM(:,:)  
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_PDF_ll_3DM
