!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/05/18 13:07:25
!-----------------------------------------------------------------
!      #############################################
       MODULE MODI_CH_JVALUES_n
!      #############################################
!
       INTERFACE
       SUBROUTINE CH_JVALUES_n(KVECNPT, KIINDEX, KJINDEX, KMODELLEVEL, KRATES, PJVALUES, PRATES)
       IMPLICIT NONE
       INTEGER, DIMENSION(:), INTENT(IN)  :: KMODELLEVEL
       INTEGER,               INTENT(IN)  :: KVECNPT 
       INTEGER, DIMENSION(:), INTENT(IN)  :: KIINDEX,KJINDEX 
                                                     ! current model grid point
       INTEGER,               INTENT(IN)  :: KRATES          
                                                     ! dimension of PRATES
       REAL, DIMENSION(KVECNPT,KRATES), INTENT(OUT) :: PRATES          
                                                     ! photolysis rates 
       REAL,DIMENSION(:,:,:,:), INTENT(IN) :: PJVALUES    ! Tuv coefficients

       END SUBROUTINE CH_JVALUES_n
       END INTERFACE
       END MODULE MODI_CH_JVALUES_n
!-----------------------------------------------------------------------------
!      ###############################################################
       SUBROUTINE CH_JVALUES_n(KVECNPT, KIINDEX, KJINDEX, KMODELLEVEL, KRATES, PJVALUES, PRATES)
!      ###############################################################
! 
!!
!!*** *CH_JVALUES* extracts photolysis rates
!!
!!    PURPOSE
!!    -------
!!      Extract photolysis rates from MODD_CH_JVALUES at a given model level
!!
!!**  METHOD
!!    ------
!!      After the photolysis rates have been set with CH_UPDATE_JVALUES_n,
!!    the different photolysis rates may be extracted on the model levels.
!       
!!    REFERENCE
!!    ---------
!!    MesoNH documentation
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!    
!!    MODIFICATIONS
!!    -------------
!!    Original 05/03/97
!!    C. Mari  20/03/01 3D version of J values + vectorisation
!!    Modification   01/12/03  (Gazen)   change Chemical scheme interface
!!
!!------------------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_INIT_JVALUES, ONLY: JPJVMAX
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
INTEGER, DIMENSION(:),   INTENT(IN)  :: KMODELLEVEL
INTEGER,                 INTENT(IN)  :: KVECNPT 
INTEGER, DIMENSION(:),   INTENT(IN)  :: KIINDEX,KJINDEX 
                                      ! current model grid point
INTEGER,                 INTENT(IN)  :: KRATES          
                                      ! dimension of PRATES
REAL, DIMENSION(KVECNPT,KRATES), INTENT(OUT) :: PRATES          
                                      ! photolysis rates 
REAL,DIMENSION(:,:,:,:), INTENT(IN) :: PJVALUES    ! Tuv coefficients
!!
!!    LOCAL VARIABLES
!!    ---------------
INTEGER :: JI,JN
!!
!------------------------------------------------------------------------------
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
DO JI = 1, KRATES
 DO JN = 1, KVECNPT
   PRATES(JN,JI) = PJVALUES(KIINDEX(JN),KJINDEX(JN),KMODELLEVEL(JN),JI)
 ENDDO
ENDDO
!
END SUBROUTINE CH_JVALUES_n
