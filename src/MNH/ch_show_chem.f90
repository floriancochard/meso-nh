!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 chimie 2006/06/30 11:42:50
!-----------------------------------------------------------------
!!    ######################## 
      MODULE MODI_CH_SHOW_CHEM
!!    ######################## 
!!
!
INTERFACE
!!
SUBROUTINE CH_SHOW_CHEM(PCONC, HNAMES)
IMPLICIT NONE
!REAL, DIMENSION(NEQ), INTENT(IN) :: PCONC ! PCONC: concentration vector
REAL, DIMENSION(:), INTENT(IN) :: PCONC ! PCONC: concentration vector
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES ! name of each species
END SUBROUTINE CH_SHOW_CHEM
!!
END INTERFACE
!!
END MODULE MODI_CH_SHOW_CHEM
!!
!!    ############################# 
      SUBROUTINE CH_SHOW_CHEM(PCONC, HNAMES)
!!    #############################
!!
!!****  *CH_SHOW_CHEM*
!!
!!    PURPOSE
!!    -------
!!    print a set of values to the screen
!!
!!    METHOD
!!    ------
!!    NEQ values will be written to the screen,
!!    the format will be as follows:
!!    line 1:     'name1'     value1
!!      ...
!!    line NEQ:   'nameNEQ'     valueNEQ
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 21/04/95
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
!!    PCONC: concentration vector to be read
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN) :: PCONC ! PCONC: concentration vector
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES ! name of each species
!
!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
INTEGER :: JI
!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
DO JI = 1, SIZE(PCONC)
  PRINT '(3A,E20.8)', "'", HNAMES(JI), "' ", PCONC(JI)
ENDDO
!
END SUBROUTINE CH_SHOW_CHEM
