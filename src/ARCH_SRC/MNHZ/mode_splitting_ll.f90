!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

!     ########################
MODULE MODE_SPLITTING_ll
  !     ########################
  !
  !!    Purpose
  !!    -------
  !
  !     The purpose of this module is to provide subroutines 
  !     for the splitting of the domain.
  !
  !!    Routines Of The User Interface
  !!    ------------------------------
  !
  !     None
  !
  !!    Reference
  !!    ---------
  !
  !     User Interface for Meso-NH parallel package
  !     Ph. Kloos, L. Giraud, R. Guivarch, D. Lugato
  !
  !!    Authors
  !!    -------
  !
  !     R. Guivarch, D. Lugato    * CERFACS *
  !     Ph. Kloos                 * CERFACS - CNRM *
  !
  !
  !!    Implicit Arguments
  !!    ------------------
  !
  !!    Implicit Arguments
  !!    ------------------
  !
  !     Module MODD_STRUCTURE_ll
  !       type ZONE_LL
  !
  !     Module MODD_VAR_ll
  !       JPHEXT
  !
  !------------------------------------------------------------------------------
  !
  !***********************************************************************
  !
  ! TYPE PRIME_DECOMPOS
  !
  !***********************************************************************
  !
  TYPE PRIME_DECOMPOS
     !
     INTEGER :: NUMBER            ! prime number
     INTEGER :: POWER	     ! power of the prime number
     !
  END TYPE PRIME_DECOMPOS
  !
  !
  INTEGER, PARAMETER :: PRIME_MAX=30    ! maximum number of prime numbers
  !
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  !
  !     ########################
  SUBROUTINE INIT_TAB(TAB)
    !     ########################
    !
    !!****  *INIT_TAB* - routine which initializes an array to decompose 
    !                    an integer with prime numbers.
    !!    Purpose
    !!    -------
    !     The purpose of this routine is to fill an array for the prime
    !     numbers decomposition. We use the first 30 prime numbers and it is
    !     sufficiant because of the maximum processors used.
    !
    !!**  Method
    !!    ------
    !     The array TAB is filled with the first 30 prime numbers with the
    !     power 0.
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    MODIFICATIONS
    !!    -------------
    !     Original 01/05/98
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PRIME_DECOMPOS), DIMENSION(PRIME_MAX) :: TAB
    !
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: J ! loop control variable
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1. ARRAY INITIALIZATION
    !
    TAB(1)%NUMBER = 2
    TAB(2)%NUMBER = 3
    TAB(3)%NUMBER = 5
    TAB(4)%NUMBER = 7
    TAB(5)%NUMBER = 11
    TAB(6)%NUMBER = 13
    TAB(7)%NUMBER = 17
    TAB(8)%NUMBER = 19
    TAB(9)%NUMBER = 23
    TAB(10)%NUMBER = 29
    TAB(11)%NUMBER = 31
    TAB(12)%NUMBER = 37
    TAB(13)%NUMBER = 41
    TAB(14)%NUMBER = 43
    TAB(15)%NUMBER = 53
    TAB(16)%NUMBER = 59
    TAB(17)%NUMBER = 61
    TAB(18)%NUMBER = 67
    TAB(19)%NUMBER = 71
    TAB(20)%NUMBER = 73
    TAB(21)%NUMBER = 79
    TAB(22)%NUMBER = 83
    TAB(23)%NUMBER = 89
    TAB(24)%NUMBER = 97
    TAB(25)%NUMBER = 101
    TAB(26)%NUMBER = 103
    TAB(27)%NUMBER = 107
    TAB(28)%NUMBER = 109
    TAB(29)%NUMBER = 113
    TAB(30)%NUMBER = 127
    !
    DO J=1,PRIME_MAX
       TAB(J)%POWER = 0
    END DO
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE INIT_TAB
  !
  !     ################################
  SUBROUTINE DECOMPOSE(TAB,N,PREM)
    !     ################################
    !
    !!****  *DECOMPOSE* - routine which decompose an integer with prime numbers.
    !
    !!    Purpose
    !!    -------
    !     this routine fill the array TAB with the prime number decomposition of N.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PRIME_DECOMPOS), DIMENSION(PRIME_MAX) :: TAB
    ! decomposition of N
    INTEGER :: N                       ! integer to decompose
    LOGICAL :: PREM                    ! True if N is prime
    !
    !*       0.2   declarations of local variables
    !
    INTRINSIC MOD
    INTEGER :: INTERM,I
    LOGICAL :: AUX
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    fILL-IN OF TAB
    !
    INTERM = N
    PREM = .TRUE.
    AUX = .FALSE.
    I = 1
    DO WHILE ((TAB(I)%NUMBER).LE.N)
       DO WHILE (MOD(INTERM,TAB(I)%NUMBER) .EQ. 0) 
          INTERM = INTERM/TAB(I)%NUMBER
          TAB(I)%POWER = TAB(I)%POWER + 1
          IF (AUX) PREM = .FALSE.
          AUX = .TRUE.
       END DO
       I = I + 1
    END DO
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE DECOMPOSE
  !
  !     #############################################################
  SUBROUTINE DEF_SPLITTING(TAB,X_DOMAINS,Y_DOMAINS,X_DIM,Y_DIM)
    !     #############################################################
    !
    !!****  *DEF_SPLITTING* - routine which define the splitting of a domain
    !                         this routine is slower than the following one.
    !!
    !!    Purpose
    !!    -------
    !     This routine define the splitting of a domain.
    !     The maximum dimension will be cut more than the other one.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !     DISPLAY_ZONE - routine which displays a structured type ZONE_ll
    !!
    !!    Implicit Arguments
    !!    ------------------
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    TYPE(PRIME_DECOMPOS), DIMENSION(PRIME_MAX) :: TAB ! number of procs
    INTEGER :: X_DOMAINS,Y_DOMAINS ! number of domains in each direction
    INTEGER :: X_DIM,Y_DIM         ! dimensions of the domain
    !
    !
    !*       0.2   declarations of local variables
    !
    TYPE(PRIME_DECOMPOS), DIMENSION(PRIME_MAX) :: INTER
    INTEGER :: I
    INTEGER :: X_LOCDIM, Y_LOCDIM  ! local dimensions of the subdomains
    !
    !-------------------------------------------------------------------------------
    !
    !*       1. COMPUTE X_DOMAINS AND Y_DOMAINS
    !
    INTER = TAB
    X_DOMAINS = 1
    Y_DOMAINS = 1
    X_LOCDIM = X_DIM / X_DOMAINS
    Y_LOCDIM = Y_DIM / Y_DOMAINS
    DO I=SIZE(TAB),1,-1
       DO WHILE (INTER(I)%POWER .NE. 0) 
          IF (X_LOCDIM.GE.Y_LOCDIM) THEN
             X_DOMAINS = X_DOMAINS*INTER(I)%NUMBER
             X_LOCDIM = X_DIM / X_DOMAINS
          ELSE
             Y_DOMAINS = Y_DOMAINS*INTER(I)%NUMBER
             Y_LOCDIM = Y_DIM / Y_DOMAINS
          END IF
          INTER(I)%POWER = INTER(I)%POWER - 1
       END DO
    END DO
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE DEF_SPLITTING
  !
  !     #######################################################################
  SUBROUTINE DEF_SPLITTING2(X_DOMAINS,Y_DOMAINS,X_DIM,Y_DIM,NB_PROC,PREM)
    !     #######################################################################
    !
    !!****  *DEF_SPLITTING2* - routine which define the splitting of a domain
    !                          with a faster method.
    !!
    !!    Purpose
    !!    -------
    !     This routine define the splitting of a domain.
    !     The greater dimension will be cut more than the other one.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    MODIFICATIONS
    !!    -------------
    !     Original 01/05/98
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    INTEGER :: NB_PROC             ! number of processors
    INTEGER :: X_DOMAINS,Y_DOMAINS ! number of domains in each direction
    INTEGER :: X_DIM,Y_DIM         ! dimensions of the domain
    LOGICAL :: PREM                ! true if nbproc is a prime number
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: I
    INTRINSIC MOD,SQRT,INT,REAL
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    FIND THE GREATER DIVISER OF NB_PROC
    !
    I= INT(SQRT(REAL(NB_PROC)))
    DO WHILE (MOD(NB_PROC,I).NE.0)
       I = I - 1
    END DO
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    AFFECT THE GREATER DIMENSION WITH THE GREATER DIVISER
    !
    IF (X_DIM.GE.Y_DIM) THEN
       IF (I.GE.(NB_PROC/I)) THEN
          X_DOMAINS = I
          Y_DOMAINS = NB_PROC/I
       ELSE
          X_DOMAINS = NB_PROC/I
          Y_DOMAINS = I
       END IF
    ELSE
       IF (I.GE.(NB_PROC/I)) THEN
          X_DOMAINS = NB_PROC/I
          Y_DOMAINS = I
       ELSE
          X_DOMAINS = I
          Y_DOMAINS = NB_PROC/I
       END IF
    END IF
    !
    !-------------------------------------------------------------------------------
    !
    !*	 3.    IS NB_PROC A PRIME NUMBER?
    !
    IF ((X_DOMAINS.EQ.1).OR.(Y_DOMAINS.EQ.1)) THEN
       PREM = .TRUE.
    ELSE
       PREM = .FALSE.
    END IF
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE DEF_SPLITTING2
  !
  !     ###################################################################
  SUBROUTINE CARTESIAN(TPROC,NB_PROC,X_DIM,Y_DIM,X_DOMAINS,Y_DOMAINS)
    !     ###################################################################
    !
    !!****  *CARTESIAN* - routine which splits a domain if NB_PROC 
    !  		      is not a prime number
    !
    !!    Purpose
    !!    -------
    !     this routine fills the elements of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    TYPE(ZONE_LL), DIMENSION(:)  :: TPROC
    INTEGER :: NB_PROC,X_DIM,Y_DIM,X_DOMAINS,Y_DOMAINS
    !
    !*       0.2   declarations of local variables
    !
    INTRINSIC SQRT,MOD
    INTEGER :: N_X,N_Y,RESTX,RESTY,AUX
    INTEGER :: SOMX,SOMY,INTX,INTY,IDOM,I,J
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    COMPUTE THE AVERAGE DIMENSION
    !
    N_X = X_DIM/X_DOMAINS
    N_Y = Y_DIM/Y_DOMAINS
    RESTX = MOD(X_DIM,X_DOMAINS)
    RESTY = MOD(Y_DIM,Y_DOMAINS)
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    FILL THE FIELDS OF TPROC
    ! 
    SOMY = 1
    DO I=1,Y_DOMAINS
       SOMX = 1
       INTY = 0
       RESTX = MOD(X_DIM,X_DOMAINS)
       IF(I .GE. (Y_DOMAINS-RESTY+1)) INTY = 1
       ! 
       DO J=1,X_DOMAINS
          IDOM = (I-1)*X_DOMAINS+J
          TPROC(IDOM)%NUMBER = IDOM
          INTX = 0
          IF(RESTX .GT. 0) INTX = 1
          TPROC(IDOM)%NXOR = SOMX
          TPROC(IDOM)%NXEND = SOMX + N_X + INTX - 1
          SOMX = SOMX + N_X + INTX
          RESTX = RESTX - 1
          TPROC(IDOM)%NYOR = SOMY
          TPROC(IDOM)%NYEND = SOMY + N_Y + INTY -1
       END DO
       !
       SOMY = SOMY + N_Y + INTY
    END DO
    ! 
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE CARTESIAN
  !
  !     ###########################################
  SUBROUTINE SPLIT(X_DIM,Y_DIM,NB_PROC,TPROC)
    !     ###########################################
    !
    !!****  *SPLIT* - routine which splits a domain in NB_PROC sub-domains
    !                     by using DEFINE_SPLITTING.
    !
    !!    Purpose
    !!    -------
    !     this routine fills the fields of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !     Module MODD_VAR_ll
    !       JPHEXT
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    INTEGER :: NB_PROC,X_DIM,Y_DIM
    TYPE(ZONE_LL), DIMENSION(NB_PROC)  :: TPROC
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: X_DOMAINS,Y_DOMAINS
    TYPE(PRIME_DECOMPOS), DIMENSION(PRIME_MAX) :: TAB
    LOGICAL :: PREM
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    PRIME NUMBER DECOMPOSITION OF NB_PROC
    !
    CALL INIT_TAB(TAB)
    CALL DECOMPOSE(TAB,NB_PROC,PREM)
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    FILL THE FIELDS OF TPROC
    ! 
    IF ((PREM).AND.(NB_PROC.GT.2)) THEN
       ! 
       CALL INIT_TAB(TAB)
       CALL DECOMPOSE(TAB,NB_PROC-1,PREM)
       !        CALL DISPLAY_TAB(TAB)
       CALL DEF_SPLITTING(TAB,X_DOMAINS,Y_DOMAINS,X_DIM,Y_DIM)
       !
       ! the last proc is a vertical sub-domain
       !
       TPROC(NB_PROC)%NUMBER = NB_PROC
       TPROC(NB_PROC)%NXOR = X_DIM-X_DIM/NB_PROC 
       TPROC(NB_PROC)%NXEND = X_DIM
       TPROC(NB_PROC)%NYOR = 1
       TPROC(NB_PROC)%NYEND = Y_DIM 
       !
       ! cartesian splitting with NB_PROC-1
       !
       CALL CARTESIAN(TPROC,NB_PROC-1,X_DIM-X_DIM/NB_PROC-1,   &
            Y_DIM,X_DOMAINS,Y_DOMAINS)
       ! 
    ELSE
       !
       CALL DEF_SPLITTING(TAB,X_DOMAINS,Y_DOMAINS,X_DIM,Y_DIM)
       CALL CARTESIAN(TPROC,NB_PROC,X_DIM,Y_DIM,X_DOMAINS,Y_DOMAINS)
    END IF
    !
    !-------------------------------------------------------------------------------
    !
    !*	 3.     SHIFT FROM PHYSICAL TO EXTENDED DOMAIN
    !
    TPROC(:)%NXOR = TPROC(:)%NXOR + JPHEXT
    TPROC(:)%NYOR = TPROC(:)%NYOR + JPHEXT
    TPROC(:)%NXEND = TPROC(:)%NXEND + JPHEXT
    TPROC(:)%NYEND = TPROC(:)%NYEND + JPHEXT
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE SPLIT
  !
  !
  !     #################################################################
  SUBROUTINE SPLIT2(X_DIM,Y_DIM,Z_DIM,NB_PROC,TPROC,HSPLITTING)
    !     #################################################################
    !
    !!****  *SPLIT2* - routine which splits a domain in NB_PROC sub-domains
    !                     by using DEFINE_SPLITTING2.
    !
    !!    Purpose
    !!    -------
    !     this routine fills the fields of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !     Module MODD_VAR_ll
    !       JPHEXT
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !     R.Guivarch 29/11/99 : x and y splitting : HSPLITTING
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT, JPVEXT
    !
    IMPLICIT NONE
    !
    !
    !*       0.1   declarations of arguments
    !
    INTEGER, INTENT(IN) :: NB_PROC,X_DIM,Y_DIM,Z_DIM
    CHARACTER*10, INTENT(IN) :: HSPLITTING ! kind of splitting
    TYPE(ZONE_LL), INTENT(OUT), DIMENSION(NB_PROC)  :: TPROC
    !
    !*       0.2   declarations of local variables
    !
    INTEGER :: X_DOMAINS,Y_DOMAINS
    LOGICAL :: PREM
    !
    !-------------------------------------------------------------------------------
    !
    !*	 1.    FIND THE SPLITTING
    !
    CALL DEF_SPLITTING2(X_DOMAINS,Y_DOMAINS,X_DIM,Y_DIM,NB_PROC,PREM)
    !
    !-------------------------------------------------------------------------------
    !
    !*	 2.    FILL THE FIELDS OF TPROC
    !
    IF(HSPLITTING.EQ."BSPLITTING") THEN
       IF ((PREM).AND.(NB_PROC.GT.2)) THEN
          CALL DEF_SPLITTING2(X_DOMAINS,Y_DOMAINS,X_DIM-X_DIM/NB_PROC-1,  &
               Y_DIM,NB_PROC-1,PREM)
          !
          ! the last proc is a vertical sub-domain
          !
          TPROC(NB_PROC)%NUMBER = NB_PROC
          TPROC(NB_PROC)%NXOR = X_DIM-X_DIM/NB_PROC
          TPROC(NB_PROC)%NXEND = X_DIM
          TPROC(NB_PROC)%NYOR = 1
          TPROC(NB_PROC)%NYEND = Y_DIM
          !
          ! cartesian splitting with NB_PROC-1
          !
          CALL CARTESIAN(TPROC,NB_PROC-1,X_DIM-X_DIM/NB_PROC-1,   &
               Y_DIM,X_DOMAINS,Y_DOMAINS)
          !
       ELSE
          !
          CALL CARTESIAN(TPROC,NB_PROC,X_DIM,Y_DIM,X_DOMAINS,Y_DOMAINS)
       END IF
    ELSEIF(HSPLITTING.EQ."XSPLITTING") THEN
       !
       CALL CARTESIAN(TPROC,NB_PROC,X_DIM,Y_DIM,NB_PROC,1)
       !
    ELSE
       !
       CALL CARTESIAN(TPROC,NB_PROC,X_DIM,Y_DIM,1,NB_PROC)
       !
    END IF
    !JUAN
    TPROC(:)%NZOR = 1
    TPROC(:)%NZEND = Z_DIM 
    !JUAN
    !
    !*       3.     shift from physical to extended domain
    !
    TPROC(:)%NXOR = TPROC(:)%NXOR + JPHEXT
    TPROC(:)%NYOR = TPROC(:)%NYOR + JPHEXT
    TPROC(:)%NZOR  = TPROC(:)%NZOR + JPVEXT
    TPROC(:)%NXEND = TPROC(:)%NXEND + JPHEXT
    TPROC(:)%NYEND = TPROC(:)%NYEND + JPHEXT
    TPROC(:)%NZEND = TPROC(:)%NZEND + JPVEXT
    !
    !-------------------------------------------------------------------------------
    !
  END SUBROUTINE SPLIT2
  !
  !     #################################################################
  SUBROUTINE SPLIT3(XDIMINT, YDIMINT, NPROC, TPROC)
    !     #################################################################
    !
    !!****  *SPLIT3* - routine which splits a domain in NB_PROC sub-domains
    !!                 when NPROC equals 2, 4 or 8
    !
    !!    Purpose
    !!    -------
    !     this routine fills the fields of TPROC.
    !
    !!**  Method
    !!    ------
    !
    !!    External
    !!    --------
    !
    !!    Implicit Arguments
    !!    ------------------
    !
    !     Module MODD_STRUCTURE_ll
    !       type ZONE_LL
    !
    !     Module MODD_VAR_ll
    !       JPHEXT
    !
    !!    Reference
    !!    ---------
    !
    !!    Author
    !!    ------
    !     D. Lugato    * CERFACS *
    !
    !!    Modifications
    !!    -------------
    !     Original 01/05/98
    !!    -------
    !     this routine fills the fields of TPROC.
    !
    !-------------------------------------------------------------------------------
    !
    !*       0.    DECLARATIONS
    !
    USE MODD_STRUCTURE_ll, ONLY : ZONE_LL
    USE MODD_PARAMETERS_ll, ONLY : JPHEXT
    !
    IMPLICIT NONE
    !
    !*       0.1   declarations of arguments
    !
    INTEGER :: XDIMINT, YDIMINT, NPROC
    TYPE(ZONE_ll), DIMENSION(NPROC) :: TPROC
    !
    !*       0.2   declarations of local variables
    !
    INTEGER ::XDIM, YDIM
    REAL :: RPROC
    !
    !-------------------------------------------------------------------------------
    !
    XDIM = XDIMINT + 2
    YDIM = YDIMINT + 2
    !
    TPROC(1)%NXOR  = 2
    RPROC = NPROC
    TPROC(1)%NXEND = (XDIM-2)/(1+(ALOG(RPROC)/ALOG(2.)))+1
    TPROC(1)%NYOR  = 2
    TPROC(1)%NYEND = YDIM - 1
    !
    TPROC(2)%NXOR  = TPROC(1)%NXEND+1
    TPROC(2)%NYOR  = 2
    !
    IF (NPROC == 2) TPROC(2)%NYEND = YDIM - 1
    !
    IF (NPROC == 2 .OR. NPROC == 4) THEN
       TPROC(2)%NXEND = XDIM - 1
    ELSEIF (NPROC == 8) THEN
       TPROC(2)%NXEND = (XDIM - 2)/2 + 1
    ENDIF
    !
    IF (NPROC == 8) THEN
       !
       TPROC(3)%NXOR = TPROC(2)%NXEND+1
       TPROC(3)%NXEND = (XDIM - 2)*3/4 
       TPROC(5)%NXOR = TPROC(3)%NXOR
       TPROC(5)%NXEND = TPROC(3)%NXEND
       TPROC(7)%NXOR = TPROC(3)%NXOR
       TPROC(7)%NXEND = TPROC(3)%NXEND
       TPROC(4)%NXOR = TPROC(3)%NXEND + 1
       TPROC(4)%NXEND = XDIM - 1
       TPROC(6)%NXOR = TPROC(4)%NXOR
       TPROC(6)%NXEND = TPROC(4)%NXEND
       TPROC(8)%NXOR = TPROC(4)%NXOR
       TPROC(8)%NXEND = TPROC(4)%NXEND
       !
       TPROC(2)%NYEND = YDIM - 1 
       TPROC(3)%NYOR = 2
       TPROC(3)%NYEND = (YDIM - 2)/3 
       TPROC(5)%NYOR = TPROC(3)%NYEND + 1
       TPROC(5)%NYEND = (YDIM - 2)*2/3 
       TPROC(7)%NYOR = TPROC(5)%NYEND + 1
       TPROC(7)%NYEND = YDIM -1
       TPROC(4)%NYOR = TPROC(3)%NYOR
       TPROC(4)%NYEND = TPROC(3)%NYEND
       TPROC(6)%NYOR = TPROC(5)%NYOR
       TPROC(6)%NYEND = TPROC(5)%NYEND
       TPROC(8)%NYOR = TPROC(7)%NYOR
       TPROC(8)%NYEND = TPROC(7)%NYEND
       !
    ENDIF
    !
    IF (NPROC == 4) THEN
       !
       TPROC(2)%NYEND = (YDIM-2)/3
       TPROC(3)%NYOR = TPROC(2)%NYEND+1
       TPROC(3)%NYEND = (YDIM-2)*2/3
       TPROC(4)%NYOR = TPROC(3)%NYEND+1
       TPROC(4)%NYEND = YDIM - 1
       !
       TPROC(3)%NXOR = TPROC(2)%NXOR
       TPROC(3)%NXEND = TPROC(2)%NXEND
       TPROC(4)%NXOR = TPROC(2)%NXOR
       TPROC(4)%NXEND = TPROC(2)%NXEND
       !
    ENDIF
    !
  END SUBROUTINE SPLIT3
  !
END MODULE MODE_SPLITTING_ll

