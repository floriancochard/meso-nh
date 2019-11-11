!MNH_LIC Copyright 2001-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ################################ 
      MODULE MODI_CH_INIT_CONST_n
!!    ################################ 
!!
INTERFACE
!!
SUBROUTINE CH_INIT_CONST_n (KLUOUT, KVERB)

IMPLICIT NONE
!!

INTEGER,                  INTENT(IN)  :: KLUOUT   ! output listing channel
INTEGER,                  INTENT(IN)  :: KVERB    ! verbosity level
!!
END SUBROUTINE CH_INIT_CONST_n
!!
END INTERFACE
!!
END MODULE MODI_CH_INIT_CONST_n
!!
!!    ###########################################
       SUBROUTINE CH_INIT_CONST_n(KLUOUT, KVERB)
!!    ###########################################
!!
!!*** *CH_INIT_CONST_n*
!!
!!    PURPOSE
!!    -------
!      Read Henry Specific constant,  Molecular Mass and Biological reactivity
!!     factor.
!!     
!!
!!**  METHOD
!!    ------
!
!!    Chemical constant will be read from
!!    the general purpose input file (variable CCHEM_INPUT_FILE). 
!!
!!       chemical molecular diffusivity MASS_MOL
!!       molecular reactivity factor REA_FACT     
!!       molecular effective Henry constant HENRY_SP
!!
!!
!!    REFERENCE
!!    ---------
!!    
!!    AUTHOR
!!    ------
!!    P. Tulet    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 16/02/01
!  P. Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet: 25/02/2019: bug correction for the file unit numbers
!!    EXTERNAL
!!    --------
USE MODI_CH_OPEN_INPUT  ! open the general purpose ASCII input file
USE MODD_IO_ll,         ONLY: TFILEDATA
USE MODE_FM,            ONLY: IO_FILE_CLOSE_ll

!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9_n,   ONLY: NEQ,               &! number of chemical variables
                          CNAMES              ! names of the chemical species
USE MODD_CH_MNHC_n, ONLY: CCHEM_INPUT_FILE    ! general purpose ASCII input file

USE MODD_CH_CONST_n
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!

INTEGER,                  INTENT(IN)  :: KLUOUT  
                          ! output listing channel
INTEGER,                  INTENT(IN)  :: KVERB   
                          ! verbosity level
!
!*      0.2    declarations of local variables
!
INTEGER           :: ICHANNEL    
                          ! I/O channel for file input
CHARACTER(LEN=40) :: YFORMAT    
                          ! format for input
CHARACTER(LEN=40) :: YOUTFORMAT = '(A32,2E15.5)'
                          ! format for output
!

INTEGER :: IMASS          ! number of molecular diffusivity to be read
CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE :: YMASSMOLNAME !species names
REAL             , DIMENSION(:), ALLOCATABLE :: ZMASSMOLVAL
                          ! molecular diffusivity value
!
INTEGER :: IREACT         ! number of chemical reactivity factor to be read
CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE :: YREACTNAME !species names
REAL             , DIMENSION(:), ALLOCATABLE :: ZREACTVAL 
                          ! chemical reactivity factor value
!
INTEGER :: IHENRY         ! number of chemical Henry constant to be read
CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE :: YHENRYNAME !species names
REAL             , DIMENSION(:,:), ALLOCATABLE :: ZHENRYVAL
                          !chemical Henry constant value
!
INTEGER :: IFAIL          ! return code from CLOSE_ll
INTEGER :: JI, JN, JNREAL ! loop control variables
INTEGER :: INACT          ! array pointer
TYPE(TFILEDATA),POINTER :: TZFILE
!-------------------------------------------------------------------------------
TZFILE => NULL()
!
!*       1.    ALLOCATE FIELD
!              --------------
!
IF (.NOT. ASSOCIATED(XSREALMASSMOLVAL)) ALLOCATE( XSREALMASSMOLVAL(NEQ) )
IF (.NOT. ASSOCIATED(XSREALREACTVAL))   ALLOCATE( XSREALREACTVAL(NEQ) )
IF (.NOT. ASSOCIATED(XSREALHENRYVAL))   ALLOCATE( XSREALHENRYVAL(NEQ,2) )
!
!*       2.  read chemical molecular diffusivity MASS_MOL
!
! open input file
IF (KVERB >= 5) WRITE(KLUOUT,*) &
   "CH_INIT_CONST: reading  molar mass"
CALL CH_OPEN_INPUT(CCHEM_INPUT_FILE, "MASS_MOL", TZFILE, KLUOUT, KVERB)
ICHANNEL = TZFILE%NLU
!
! read number of molecular diffusivity IMASS
READ(ICHANNEL, *) IMASS
IF (KVERB >= 5) WRITE(KLUOUT,*) "number of molecular diffusivity: ", IMASS
!
! read data input format
READ(ICHANNEL,"(A)") YFORMAT
IF (KVERB >= 5) WRITE(KLUOUT,*) "input format is: ", YFORMAT
!
! allocate fields
ALLOCATE(YMASSMOLNAME(IMASS))
ALLOCATE(ZMASSMOLVAL(IMASS))
!
! read molecular diffusivity
DO JI = 1, IMASS
  READ(ICHANNEL,YFORMAT) YMASSMOLNAME(JI), ZMASSMOLVAL(JI)
  WRITE(KLUOUT,YFORMAT) YMASSMOLNAME(JI), ZMASSMOLVAL(JI)
END DO
!
! close file
CALL IO_FILE_CLOSE_ll(TZFILE)
TZFILE => NULL()
!
!
IF (KVERB >= 10) THEN
  WRITE(KLUOUT,'(A)') '----------------------------------------------------'
  WRITE(KLUOUT,'(A)') 'MASS_MOL'
  WRITE(KLUOUT,'(A)') 'molecular mass (in g/mol) for molecular diffusion'
  WRITE(KLUOUT,'(I4)') NEQ
  WRITE(KLUOUT,'(A)') YOUTFORMAT
END IF
XSREALMASSMOLVAL(:) = 50. ! default molecular mass 
DO JNREAL = 1, NEQ
  INACT = 0
  search_loop3 : DO JN = 1, IMASS
    IF (CNAMES(JNREAL) .EQ. YMASSMOLNAME(JN)) THEN
      INACT = JN
      EXIT search_loop3
    END IF
  END DO search_loop3
  IF (INACT .NE. 0) XSREALMASSMOLVAL(JNREAL) = ZMASSMOLVAL(INACT)
  IF (KVERB >= 10) THEN
    WRITE(KLUOUT,YOUTFORMAT) CNAMES(JNREAL), XSREALMASSMOLVAL(JNREAL)
  END IF
END DO
!
!-----------------------------------------------------------------------------
!
!*       3.  read molecular reactivity factor REA_FACT
!
! open input file
IF (KVERB >= 5) WRITE(KLUOUT,*) &
   "CH_INIT_CONST: reading  reactivity factor "
CALL CH_OPEN_INPUT(CCHEM_INPUT_FILE, "REA_FACT", TZFILE, KLUOUT, KVERB)
ICHANNEL = TZFILE%NLU
!
! read number of molecular diffusivity IREACT
READ(ICHANNEL, *) IREACT
IF (KVERB >= 5) WRITE(KLUOUT,*) "number of reactivity factor : ", IREACT
!
! read data input format
READ(ICHANNEL,"(A)") YFORMAT
IF (KVERB >= 5) WRITE(KLUOUT,*) "input format is: ", YFORMAT
!
! allocate fields
ALLOCATE(YREACTNAME(IREACT))
ALLOCATE(ZREACTVAL(IREACT))
! read reactivity factor 
DO JI = 1, IREACT
  READ(ICHANNEL,YFORMAT) YREACTNAME(JI), ZREACTVAL(JI)
  WRITE(KLUOUT,YFORMAT) YREACTNAME(JI), ZREACTVAL(JI)
END DO
!
! close file
CALL IO_FILE_CLOSE_ll(TZFILE)
TZFILE => NULL()
!
!
IF (KVERB >= 10) THEN
  WRITE(KLUOUT,'(A)') '----------------------------------------------------'
  WRITE(KLUOUT,'(A)') 'REA_FACT'
  WRITE(KLUOUT,'(A)') 'reactivity factor'
  WRITE(KLUOUT,'(I4)') NEQ
  WRITE(KLUOUT,'(A)') YOUTFORMAT
END IF
XSREALREACTVAL(:) = 0.0 ! default (high surface resistance)
DO JNREAL = 1, NEQ
  INACT = 0
  search_loop4 : DO JN = 1, IREACT
    IF (CNAMES(JNREAL) .EQ. YREACTNAME(JN)) THEN
      INACT = JN
      EXIT search_loop4
    END IF
  END DO search_loop4
  IF (INACT .NE. 0) XSREALREACTVAL(JNREAL) = ZREACTVAL(INACT)
  IF (KVERB >= 10) THEN
    WRITE(KLUOUT,YOUTFORMAT) CNAMES(JNREAL), XSREALREACTVAL(JNREAL)
  END IF
END DO
!
!-----------------------------------------------------------------------------
!
!*       4.  read molecular effective  Henry constant HENRY_SP
!
! open input file
IF (KVERB >= 5) WRITE(KLUOUT,*) &
   "CH_INIT_CONST: reading effective Henry constant", &
   " and its temperature correction "
CALL CH_OPEN_INPUT(CCHEM_INPUT_FILE, "HENRY_SP", TZFILE, KLUOUT, KVERB)
ICHANNEL = TZFILE%NLU
!
! read number of molecular diffusivity IHENRY
READ(ICHANNEL, *) IHENRY
IF (KVERB >= 5) WRITE(KLUOUT,*) "number of reactivity factor : ", IHENRY
!
! read data input format
READ(ICHANNEL,"(A)") YFORMAT
IF (KVERB >= 5) WRITE(KLUOUT,*) "input format is: ", YFORMAT
!
! allocate fields
ALLOCATE(YHENRYNAME(IHENRY))
ALLOCATE(ZHENRYVAL(IHENRY,2))
!
! read reactivity factor 
DO JNREAL = 1, IHENRY
  READ(ICHANNEL,YFORMAT) YHENRYNAME(JNREAL), ZHENRYVAL(JNREAL,1),&
                         ZHENRYVAL(JNREAL,2)
  WRITE(KLUOUT,YFORMAT) YHENRYNAME(JNREAL), ZHENRYVAL(JNREAL,1),&
                         ZHENRYVAL(JNREAL,2)
END DO
!
! close file
CALL IO_FILE_CLOSE_ll(TZFILE)
TZFILE => NULL()
!
IF (KVERB >= 10) THEN
  WRITE(KLUOUT,'(A)') '----------------------------------------------------'
  WRITE(KLUOUT,'(A)') 'HENRY_SP'
  WRITE(KLUOUT,'(A)') 'Henrys law constants factor / exponent'
  WRITE(KLUOUT,'(I4)') NEQ
  WRITE(KLUOUT,'(A)') YOUTFORMAT
END IF
XSREALHENRYVAL(:,1) = 1E-8 ! no deposition; low Henry constant
XSREALHENRYVAL(:,2) = 0. ! 
DO JNREAL = 1, NEQ
  INACT = 0
  search_loop5 : DO JN = 1, IHENRY
    IF (CNAMES(JNREAL) .EQ. YHENRYNAME(JN)) THEN
      INACT = JN
      EXIT search_loop5
    END IF
  END DO search_loop5
  IF (INACT .NE. 0) XSREALHENRYVAL(JNREAL,1) = ZHENRYVAL(INACT,1)
  IF (INACT .NE. 0) XSREALHENRYVAL(JNREAL,2) = ZHENRYVAL(INACT,2)
  IF (KVERB >= 10) THEN
    WRITE(KLUOUT,YOUTFORMAT) CNAMES(JNREAL), &
      XSREALHENRYVAL(JNREAL,1),&
      XSREALHENRYVAL(JNREAL,2)
  END IF
END DO
!
!
END SUBROUTINE CH_INIT_CONST_n
