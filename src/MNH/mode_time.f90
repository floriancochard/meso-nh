!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODE_TIME
!     ####################
!
!!****  *MODE_TIME* -  module for time routines 
!!
!!    PURPOSE
!!    -------
!       The purpose of this executive module  is to package 
!     the routines SM_PRINT_TIME 
!    
!      
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_TIME : contains definition of types for time variables          
!!                          and time variable for all model
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/07/94 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_TIME
!
IMPLICIT NONE
!-------------------------------------------------------------------------------
!
CONTAINS
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
!*       1.   ROUTINE SM_PRINT_TIME
!             ---------------------
!-------------------------------------------------------------------------------
!     #####################################################
      SUBROUTINE SM_PRINT_TIME(TPDATETIME,TPOUTFILE,HTITLE)
!     #####################################################
!
!!****  *SM_PRINT_TIME * - routine to print a variable of type DATE_TIME
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to print a variable of type DATE_TIME
!      
!
!!**  METHOD
!!    ------
!!       The date and time are printed with or without a title.
!!   If it is an idealized case, no date is printed (only time).
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_TYPE_TIME : contains definition of types for time variables
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                        07/07/94 
!!      updated    V. Ducrocq           23/08/94                   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DATE_TIME),   INTENT(IN)           :: TPDATETIME  ! Date and time variable
TYPE(TFILEDATA),   INTENT(IN)           :: TPOUTFILE   ! Output listing file
CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: HTITLE      ! Title for Date and time
                                                       ! variable 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IHOUR,IMINUTE
REAL    :: ZSECOND,ZREMAIN
INTEGER :: ILUOUT
!-------------------------------------------------------------------------------
!
!*       1.    CONVERT TIME IN HOURS,MINUTES AND SECONDS :
!              ------------------------------------------
!
IHOUR   = INT(TPDATETIME%TIME/3600.)
ZREMAIN = MOD(TPDATETIME%TIME,3600.)
IMINUTE = INT(ZREMAIN/60.)
ZSECOND = MOD(ZREMAIN,60.)
!
!-------------------------------------------------------------------------------
!
!*       2.    PRINT ON OUTPUT-LISTING
!              -----------------------
!
ILUOUT = TPOUTFILE%NLU
!
IF (PRESENT(HTITLE)) THEN
  IF ((TPDATETIME%TDATE%YEAR < 0).OR.(TPDATETIME%TDATE%MONTH < 0).OR.    &
     (TPDATETIME%TDATE%DAY < 0) ) THEN 
    WRITE(UNIT=ILUOUT,FMT='(1X,A," :",2X,I2.2,"H",I2.2,"M", &
         & F5.2,"S")') HTITLE, IHOUR,IMINUTE,ZSECOND
  ELSE
    WRITE(UNIT=ILUOUT,FMT='(1X,A," :",I4.4,I2.2,I2.2,2X,I2.2,"H",I2.2,"M", &
         & F5.2,"S")') HTITLE, TPDATETIME%TDATE, IHOUR,IMINUTE,ZSECOND
  END IF
ELSE
  IF ((TPDATETIME%TDATE%YEAR < 0).OR.(TPDATETIME%TDATE%MONTH < 0).OR.    &
     (TPDATETIME%TDATE%DAY < 0) ) THEN 
    WRITE(UNIT=ILUOUT,FMT='(1X,2X,I2.2,"H",I2.2,"M", &
         & F5.2,"S")') IHOUR,IMINUTE,ZSECOND  
  ELSE 
    WRITE(UNIT=ILUOUT,FMT='(1X,I4.4,I2.2,I2.2,2X,I2.2,"H",I2.2,"M", &
         & F5.2,"S")') TPDATETIME%TDATE, IHOUR,IMINUTE,ZSECOND  
  END IF
END IF
!-------------------------------------------------------------------------------
!
END SUBROUTINE SM_PRINT_TIME
!-------------------------------------------------------------------------------
!
END MODULE MODE_TIME
