!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!--------------- C. Fischer 30/09/94
!      @(#) Lib:/mesonh/sources/init/s.set_dim.f90, Version:1.9, Date:98/06/23, Last modified:98/06/04
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_SET_DIM
!     ###################
!
INTERFACE
!
SUBROUTINE SET_DIM(HINIFILE,HLUOUT,KIINF,KISUP,KJINF,KJSUP,         &
                   KIMAX,KJMAX,KKMAX)
CHARACTER (LEN=*), INTENT(IN)  :: HINIFILE ! Name of the initial file 
CHARACTER (LEN=*), INTENT(IN)  :: HLUOUT   ! name for output-listing
                                           !  of nested models
INTEGER,         INTENT(INOUT) :: KIINF    !  Lower bound  in x direction of the 
                                           ! arrays in Initialization or in 
                                           ! Post-processing subroutines
INTEGER,         INTENT(INOUT) :: KISUP    !  Upper bound  in x direction of the 
                                           ! arraysin Initialization or in 
                                           ! Post-processing subroutines
INTEGER,        INTENT(INOUT)  :: KJINF    !  Lower bound  in y direction of the 
                                           ! arrays in Initialization or in 
                                           ! Post-processing subroutines
INTEGER,        INTENT(INOUT)  :: KJSUP    !  Upper bound  in y direction of the 
                                           ! arraysin Initialization or in 
                                           ! Post-processing subroutines
INTEGER,           INTENT(OUT) :: KIMAX    !  Dimension in x direction of the 
                                           ! arrays  stored in LFIFM file
INTEGER,           INTENT(OUT) :: KJMAX    !  Dimension in y direction of the 
                                           ! arrays  stored in LFIFM file
INTEGER,           INTENT(OUT) :: KKMAX    !  Dimension in z direction of the  
                                           ! arrays stored in LFIFM file
END  SUBROUTINE SET_DIM
!
END INTERFACE
!
END MODULE MODI_SET_DIM
!
!
!
!     ##############################################################
      SUBROUTINE SET_DIM(HINIFILE,HLUOUT,KIINF,KISUP,KJINF,KJSUP,         &
                         KIMAX,KJMAX,KKMAX)
!     ##############################################################
!
!!****  *SET_DIM* - routine to set model dimensions
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set dimensions of the model
!
!
!!**  METHOD
!!    ------
!!      The dimensions KIMAX,KJMAX,KKMAX are read in initial file. 
!!      Then, the horizontal dimensions of arrays are deduced :
!!        - If it is a segment achievement  configuration (CCONF='START' or 
!!     'RESTA'), the horizontal dimensions of the arrays are :
!!            KIINF=1, KISUP=KIMAX+2*JPHEXT
!!            KJINF=1, KJSUP=KJMAX+2*JPHEXT
!!        - If it is a postprocessing configuration (CCONF='POSTP'), 
!!     an horizontal window is possible ; KIINF, KISUP, 
!!     KJINF,KJSUP are the values read in EXSEG file, except when :
!!             * KIINF is  greater than KIMAX + 2*JPHEXT . Then it is set
!!     equal to KIMAX + 2*JPHEXT
!!             * KISUP is  greater than KIMAX + 2*JPHEXT . Then it is set
!!     equal to KIMAX + 2*JPHEXT
!!             * KJINF is  greater than KJMAX + 2*JPHEXT . Then it is set
!!     equal to KJMAX + 2*JPHEXT
!!             * KJSUP is  greater than KJMAX + 2*JPHEXT . Then it is set
!!     equal to KJMAX + 2*JPHEXT
!!             * KIINF or KISUP is less or equal to zero. It means that there
!!     is no window in x direction. Then, KIINF is set equal to 1 and  KISUP 
!!      is set equal to KIMAX + 2*JPHEXT.
!!             * KJINF or KJSUP is less or equal to zero. It means that there
!!     is no window in x direction. Then, KJINF is set equal to 1 and  KJSUP 
!!      is set equal to KJMAX + 2*JPHEXT.
!!
!!             
!!      
!!    EXTERNAL
!!    --------   
!!      FMREAD      : to read data in LFIFM file 
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_PARAMETERS : contains declaration of parameter variables
!!
!!        JPHEXT : Horizontal external points number
!!        JPVEXT : Vertical  external points number
!!
!!      Module MODD_CONF  : contains declaration of configuration variables
!!
!!         CCONF      : configuration of models
!!                          'START' for start configuration
!!                          'RESTA' for restart configuration
!!                          'POSTP' for post-processing configuration
!!         NVERB      : Level of informations on output-listing
!!                          0 for minimum  prints
!!                          5 for intermediate level of prints
!!                         10 for maximum  prints 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine SET_DIM)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    14/06/94 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
USE MODD_PARAMETERS
USE MODD_CONF
!
USE MODI_FMREAD
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
CHARACTER (LEN=*), INTENT(IN)  :: HINIFILE ! Name of the initial file 
CHARACTER (LEN=*), INTENT(IN)  :: HLUOUT   ! name for output-listing
                                           !  of nested models
INTEGER,         INTENT(INOUT) :: KIINF    !  Lower bound  in x direction of the 
                                           ! arrays in Initialization or in 
                                           ! Post-processing subroutines
INTEGER,         INTENT(INOUT) :: KISUP    !  Upper bound  in x direction of the 
                                           ! arraysin Initialization or in 
                                           ! Post-processing subroutines
INTEGER,        INTENT(INOUT)  :: KJINF    !  Lower bound  in y direction of the 
                                           ! arrays in Initialization or in 
                                           ! Post-processing subroutines
INTEGER,        INTENT(INOUT)  :: KJSUP    !  Upper bound  in y direction of the 
                                           ! arraysin Initialization or in 
                                           ! Post-processing subroutines
INTEGER,           INTENT(OUT) :: KIMAX    !  Dimension in x direction of 
                                           ! the physical part of the  
                                           ! arrays  stored in LFIFM file
INTEGER,           INTENT(OUT) :: KJMAX    !  Dimension in y direction of the 
                                           !  physical part of the  
                                           ! arrays  stored in LFIFM file
INTEGER,           INTENT(OUT) :: KKMAX    !  Dimension in z direction of the  
                                           !  physical part of the  
                                           ! arrays stored in LFIFM file
!
!*       0.2   declarations of local variables
!
INTEGER             :: ILENG,IGRID,ILENCH,IRESP  !   File 
CHARACTER (LEN=16)  :: YRECFM                    ! management
CHARACTER (LEN=100) :: YCOMMENT                  ! variables  
INTEGER             :: ILUOUT                    ! Logical unit number for
                                                 ! output-listing
!
!-------------------------------------------------------------------------------
!
!*       1.    READ DIMENSIONS OF ARRAYS IN LFIFM FILE
!              ---------------------------------------
!
YRECFM='IMAX'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,KIMAX,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='JMAX'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,KJMAX,IGRID,ILENCH,YCOMMENT,IRESP)
!
YRECFM='KMAX'
ILENG=1
CALL FMREAD(HINIFILE,YRECFM,HLUOUT,ILENG,KKMAX,IGRID,ILENCH,YCOMMENT,IRESP)
!
!
!-------------------------------------------------------------------------------
!
!*       2.    SET DIMENSIONS FOR ARRAY IN INITIALIZATION OR POST-PROCESSING
!              -------------------------------------------------------------
!
IF (CCONF == 'POSTP') THEN 
!              
  IF ((KIINF <= 0).OR.(KISUP <= 0)) THEN ! this condition corresponds to a 
    KIINF = 1                            ! post-processing case where the whole 
    KISUP = KIMAX + 2*JPHEXT             ! simulation domain must be considered
                                         ! along the x direction
  ELSE                                                           
    KIINF = MIN(KIINF,KIMAX+2*JPHEXT)    ! post-processing case with an 
    KISUP = MIN(KISUP,KIMAX+2*JPHEXT)    ! explicit window
  END IF 
!                                                        
  IF ((KJINF <= 0 ).OR.(KJSUP <= 0 )) THEN                       
    KJINF = 1
    KJSUP = KJMAX + 2* JPHEXT 
  ELSE                                                           
    KJINF = MIN(KJINF,KJMAX+2*JPHEXT) 
    KJSUP = MIN(KJSUP,KJMAX+2*JPHEXT)
  END IF 
!                                                        
ELSE 
!                                                            
  KIINF = 1                             ! case corresponding to a simulation
  KISUP = KIMAX + 2* JPHEXT
  KJINF = 1
  KJSUP = KJMAX+ 2* JPHEXT
!
END IF                                                           
!
!-------------------------------------------------------------------------------
!
!*       3.    PRINT DIMENSIONS ON OUTPUT_LISTING
!              ----------------------------------
!
CALL FMLOOK(HLUOUT,HLUOUT,ILUOUT,IRESP)
IF(KIINF > KISUP) THEN
  WRITE(UNIT=ILUOUT,FMT="(' THE PROGRAM STOPS IN THE SET_DIM SUBROUTINE ',/,&
              & 'BECAUSE THE WINDOW BOUNDS ARE NOT CONSISTENT ',/,          &
              & 'KIINF =',I5,' KISUP =',I5,' KJINF =',I5,' KJSUP =',I5)")   &
                 KIINF,KISUP,KJINF,KJSUP 
  STOP
END IF
! 
IF (NVERB >= 5) THEN
  WRITE(UNIT=ILUOUT,FMT="(' DIMENSIONS INITIALIZED BY SET_GRID :',/,       &
              & 'KIMAX =',I5,' KJMAX =',I5,' KKMAX =',I5,/,                 &
              & 'KIINF =',I5,' KISUP =',I5,' KJINF =',I5,' KJSUP =',I5)")   &
                 KIMAX,KJMAX,KKMAX,KIINF,KISUP,KJINF,KJSUP
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_DIM  
