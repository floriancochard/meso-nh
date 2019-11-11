!MNH_LIC Copyright 2010-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!---------------------------------------------------------------------
MODULE MODE_GATHERALL_BORDER_ll

CONTAINS

SUBROUTINE GATHERALL_BORDER_ll(HCODE,PA,PA_ll,KINFO_ll)
!!
!!***   *GATHERALL_BORDER_ll* - gathers one of the 4 (North, East, south or West)
!!                              border of the total domain
!!
!!    PURPOSE
!!    -------
!!
!!       gathers one of the 4 (North, East, south or West)
!!       border of the total domain.
!!
!!**  METHOD
!!    ------
!!
!!      The method MUST BE OPTIMIZED by someone knowing exchange parallel
!!      functions. For the time being a complete 2D field is exchange, while
!!      1 vector (along X or Y) exchange is really needed !!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. Masson   Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     08/2010
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!

USE MODD_PARAMETERS
USE MODE_ll
USE MODE_GATHER_ll
!
!---------------------------------------------------------------------
!
!*  dummy arguments
CHARACTER(LEN=1),   INTENT(IN)  :: HCODE    ! border treated 'N', 'E', 'S', 'W'
REAL, DIMENSION(:), INTENT(IN)  :: PA       ! border on each processor
REAL, DIMENSION(:), INTENT(OUT) :: PA_ll    ! true border of all domain
INTEGER,            INTENT(OUT) :: KINFO_ll ! return code

!---------------------------------------------------------------------

!
!* local variables
INTEGER :: IIU
INTEGER :: IJU
INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: IIU_ll
INTEGER :: IJU_ll
INTEGER :: IIB_ll, IIE_ll, IJB_ll, IJE_ll
INTEGER :: IIMAX_ll, IJMAX_ll
INTEGER :: JI_ll, JJ_ll
REAL, DIMENSION(:,:), ALLOCATABLE :: ZA_ll
!---------------------------------------------------------------------
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
IIU = IIE+JPHEXT
IJU = IJE+JPHEXT

CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
IIU_ll = IIMAX_ll+2*JPHEXT
IJU_ll = IJMAX_ll+2*JPHEXT
IIB_ll = 1       +  JPHEXT
IIE_ll = IIU_ll  -  JPHEXT
IJB_ll = 1       +  JPHEXT
IJE_ll = IJU_ll  -  JPHEXT
!
! Cette partie doit etre optimisee une fois que l'on disposera des fonctions de
! parallelisation adequates. Actuellement, il n'est pas possible d'utiliser les
! GATHERALL_ll('XX' ou 'YY'   car on ne dispose des valeurs de la variable (PA)
! que sur le processeur (par exemple si HCODE='N') touchant le bord nord, 
! et pas sur l'ensemble des processeurs au sud de celui-ci.
!
! Donc pour le moment, on rappatrie l'ensemble des points brutalement, même si
! au final on ne traite que la premiere(derniere) ligne ou premiere(derniere) 
! colonne du domaine global.
!
ALLOCATE(ZA_ll(IIU_ll,IJU_ll))

IF (HCODE=='E') THEN
  CALL GATHERALL_FIELD_ll('XY',SPREAD(PA,1,IIU),ZA_ll,KINFO_ll)
  PA_ll(:) =  ZA_ll(IIE_ll,:)
END IF

IF (HCODE=='W') THEN
  CALL GATHERALL_FIELD_ll('XY',SPREAD(PA,1,IIU),ZA_ll,KINFO_ll)
  PA_ll(:) =  ZA_ll(IIB_ll,:)
END IF

IF (HCODE=='N') THEN
  CALL GATHERALL_FIELD_ll('XY',SPREAD(PA,2,IJU),ZA_ll,KINFO_ll)
  PA_ll(:) = ZA_ll(:,IJE_ll)
END IF

IF (HCODE=='S') THEN
  CALL GATHERALL_FIELD_ll('XY',SPREAD(PA,2,IJU),ZA_ll,KINFO_ll)
  PA_ll(:) = ZA_ll(:,IJB_ll)
END IF
!
DEALLOCATE(ZA_ll)

!---------------------------------------------------------------------
END SUBROUTINE GATHERALL_BORDER_ll


END MODULE MODE_GATHERALL_BORDER_ll
