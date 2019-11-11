!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 prep_real 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_BARNES_FILTER
!     #########################
INTERFACE
!
      SUBROUTINE BARNES_FILTER(PDATIN,KX,PLAMBDA,PGRDATA)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDATIN ! 3D array of the field to filter 
INTEGER,                INTENT(IN) :: KX     !space coefficient
REAL,                   INTENT(IN) :: PLAMBDA!slope coefficient
REAL, DIMENSION(:,:,:), INTENT(OUT):: PGRDATA! 3D array of the filtered field
!
END SUBROUTINE BARNES_FILTER
!
END INTERFACE
!
END MODULE MODI_BARNES_FILTER
!
!
!
!     ###################################################
      SUBROUTINE BARNES_FILTER(PDATIN,KX,PLAMBDA,PGRDATA) 
!     ###################################################
!
!!****  *BARNES_FILTER* - separate a perturbed fiels from large-scale analyses
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to remove characterictics of a vortex 
!       from large-scale analyse (ECMWF,Arpege,...) using a 
!       filtering technique. This method can be used on any of pronostic
!       variables (wind,temperature,pressure,moisture,...)
!
!!**  METHOD
!!    ------
!!      
!!        This technique uses an exponential weighting function and requires
!!    only one iteration. Development of this scheme can found in the
!!    following:
!!              BARNES (64), Mon., Wea., Rev.
!!              BARNES (73), NOAA/ERL TECH MEMO.
!!        This first filtering (pass low) discriminates the total disturbed 
!!    field from the total field. The scheme takes into account an
!!    exponential weighting function.
!!     
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	O. Nuissier           * L.A. *
!!      R. Rogers             * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODD_GRID1, ONLY: XXHAT,XYHAT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDATIN ! 3D array of the field to filter 
                                             ! in the weighting function:
INTEGER,                INTENT(IN) :: KX     !space coefficient
REAL,                   INTENT(IN) :: PLAMBDA!slope coefficient
REAL, DIMENSION(:,:,:), INTENT(OUT):: PGRDATA! 3D array of the filtered field
!
!*       0.2   Declarations of local variables
!
REAL, DIMENSION(SIZE(PDATIN,1),SIZE(PDATIN,2))  :: ZRSQ,ZTEST,ZWGT
REAL                                       :: ZKMAP,ZEXP
                                           ! intermediate sums of wgted data
REAL, DIMENSION(SIZE(PDATIN,1),SIZE(PDATIN,2),SIZE(PDATIN,3)):: ZSUMONE,ZSUMTWO
INTEGER                                    :: IMAXI,IMINI
INTEGER                                    :: IMAXJ,IMINJ
INTEGER                                    :: IIX,IIY,IIZ !sizes of input arrays
INTEGER                                    :: JI,JJ,JK,JI2,JJ2  ! Loop indexes 
!
!-------------------------------------------------------------------------------
!
!*	 1.     SET CONSTANTS 
!               -------------
!
IIX=SIZE(PDATIN,1)
IIY=SIZE(PDATIN,2)
IIZ=SIZE(PDATIN,3)
!
ZKMAP  =KX*KX        ! en nb de points
ZEXP = -1. / (ZKMAP * PLAMBDA)
!
!-------------------------------------------------------------------------------
!
!*       2.     BEGIN THE LOOP THROUGH THE GRID POINTS, THE WEIGHTING
!           SCHEME VARIES DEPENDING ON WHETHER THIS IS THE INITIAL PASS-1
!           OR THE CORRECTION PASS-2 (only the 2nd is used here)
!           -------------------------------------------------------------
!
ZSUMONE(:,:,:) = 0.
ZSUMTWO(:,:,:) = 0.
PGRDATA(:,:,:) = 0.
!
!DO JK = 1, IIZ
  DO JJ = 1, IIY
    DO JI = 1, IIX
!
!*	 2.1    Define max and min values to limit the number of the grid points
!            to be considered
!
      IMAXI = MIN(JI+KX,IIX)
      IMINI = MAX(JI-KX,1)
      IMAXJ = MIN(JJ+KX,IIY)
      IMINJ = MAX(JJ-KX,1)
!
!*       2.2    Compute the grid points values with the weighting function
!
      DO JJ2 = IMINJ, IMAXJ
        DO JI2 = IMINI, IMAXI
          ZRSQ(JI2,JJ2) = (JI-JI2)**2 + (JJ-JJ2)**2
          ZTEST(JI2,JJ2) = 0.5 + SIGN(0.5,ZKMAP-ZRSQ(JI2,JJ2)) 
          ! better for vectorization than: IF (ZRSQ>ZKMAP) CYCLE
!          ZWGT(JI2,JJ2) = EXP(- ZRSQ(JI2,JJ2) / (ZKMAP * PLAMBDA))*ZTEST(JI2,JJ2)
          ZWGT(JI2,JJ2) = EXP(ZRSQ(JI2,JJ2) * ZEXP) * ZTEST(JI2,JJ2)
        END DO
      END DO

      DO JK = 1, IIZ
        DO JJ2 = IMINJ, IMAXJ
          DO JI2 = IMINI, IMAXI
            ZSUMTWO(JI2,JJ2,JK) = ZSUMTWO(JI2,JJ2,JK) + ZWGT(JI2,JJ2)
            ZSUMONE(JI2,JJ2,JK) = ZSUMONE(JI2,JJ2,JK) + ZWGT(JI2,JJ2) * PDATIN(JI,JJ,JK)
          END DO
        END DO
      END DO
!      
    END DO
  END DO
!END DO
!
!*       2.3    Fill the output array
!
DO JK = 1, IIZ
  DO JJ = 1, IIY 
    DO JI = 1, IIX
      IF (ZSUMTWO(JI,JJ,JK) /= 0.) THEN
        PGRDATA(JI,JJ,JK) = ZSUMONE(JI,JJ,JK) / ZSUMTWO(JI,JJ,JK)
      ENDIF
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE BARNES_FILTER
