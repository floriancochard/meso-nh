!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #############################################################
#ifdef MNH_PARALLEL
      SUBROUTINE SPLIT_GRID_PARAMETERX1_MNH(HGRID,HREC,KDIM,KSIZE,KIMAX_ll,KJMAX_ll,KHALO,PFIELD,PFIELD_SPLIT)
#else
      SUBROUTINE SPLIT_GRID_PARAMETERX1_MNH(HGRID,HREC,KDIM,KSIZE,PFIELD,PFIELD_SPLIT)
#endif
!     #############################################################
!
!!****  * - routine to split a real array on the splitted grid 
!
!	Modifications
!  M.Moge  10/02/2015  Using local subdomain for parallel execution
!  M.Moge  01/03/2015  Using KIMAX_ll,KJMAX_ll,KHALO for the call to SPLIT_GRID in subroutine PGD_GRID
!
!
#ifdef MNH_PARALLEL
#else
USE MODD_IO_SURF_MNH, ONLY : NHALO
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=10),      INTENT(IN) :: HGRID       ! grid type
CHARACTER(LEN=6),       INTENT(IN) :: HREC        ! name of the parameter
INTEGER,                INTENT(IN) :: KDIM        ! size of PFIELD
INTEGER,                INTENT(IN) :: KSIZE       ! size of PFIELD_SPLIT
#ifdef MNH_PARALLEL
INTEGER,                INTENT(IN) :: KIMAX_ll    !(global) dimension of the domain - X direction
INTEGER,                INTENT(IN) :: KJMAX_ll    !(global) dimension of the domain - Y direction
INTEGER,                INTENT(IN) :: KHALO ! size of the Halo
#endif
REAL, DIMENSION(KDIM ), INTENT(IN) :: PFIELD      ! real field for complete grid
REAL, DIMENSION(KSIZE), INTENT(OUT):: PFIELD_SPLIT! real field for splitted grid
!
!*      0.2   Declarations of local variables
!
INTEGER :: JI, JJ
INTEGER :: IIB, IIE, IJB, IJE
INTEGER :: NIMAX, NJMAX  !(local) dimensions of the subdomain
INTEGER :: IXOR, IYOR    !origin of local subdomain
!
REAL, DIMENSION(:), ALLOCATABLE :: ZX, ZY
REAL                            :: ZDX, ZDY
!
INTEGER :: IINDEX
#ifdef MNH_PARALLEL
#else
INTEGER :: KIMAX_ll
INTEGER :: KJMAX_ll
INTEGER :: KHALO
!
KHALO = NHALO
CALL GET_GLOBALDIMS_ll (KIMAX_ll,KJMAX_ll)
#endif
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_OR_ll ('B',IXOR,IYOR)
NIMAX=IIE-IIB+1
NJMAX=IJE-IJB+1
!

!
IF (HREC=='XX' .OR. HREC=='DX') THEN
  ALLOCATE(ZX(NIMAX+2*KHALO))
  ZX(1+KHALO:NIMAX+KHALO) = PFIELD(IXOR:IXOR+NIMAX-1)
  IF (HREC=='DX') THEN
    ZDX = PFIELD(1)
    DO JI=KHALO,1,-1
      ZX(JI) = ZDX
      ZX(NIMAX+2*KHALO-JI+1) = ZDX
    END DO
  ELSE IF (HREC=='XX') THEN
    ZDX = 0.
    IF (NIMAX>1) ZDX = PFIELD(2) - PFIELD(1)
    IF (NIMAX==1) ZDX = PFIELD(1) ! in 1D conf, one assumes that grid
                                     ! is located between X=DX/2 and X=3DX/2
    DO JI=KHALO,1,-1
      ZX(JI) = ZX(JI+1) - ZDX
      ZX(NIMAX+2*KHALO-JI+1) = ZX(NIMAX+2*KHALO-JI) + ZDX
    END DO
  END IF

!
  DO JJ=1,NJMAX+2*KHALO
    DO JI=1,NIMAX+2*KHALO
      IINDEX = JI+(JJ-1)*(NIMAX+2*KHALO)
      PFIELD_SPLIT(IINDEX) = ZX(JI)
    END DO
  END DO
  DEALLOCATE(ZX)
  
ELSEIF (HREC=='YY' .OR. HREC=='DY') THEN
  ALLOCATE(ZY(NJMAX+2*KHALO))
  DO JJ=1+KHALO,NJMAX+KHALO
    ZY(JJ) = PFIELD(1+(IYOR+JJ-1-1-KHALO)*KIMAX_ll)
  END DO
  IF (HREC=='DY') THEN
    ZDY = PFIELD(1)
    DO JJ=KHALO,1,-1
      ZY(JJ) = ZDY
      ZY(NJMAX+2*KHALO-JJ+1) = ZDY
    END DO
  ELSE IF (HREC=='YY') THEN
    ZDY = 0.
    IF (NJMAX>1) ZDY = PFIELD(1+KIMAX_ll) - PFIELD(1)
    IF (NJMAX==1) ZDY = PFIELD(1) ! in 1D or 2D conf, one assumes that grid
                                     ! is located between Y=DY/2 and Y=3DY/2
    DO JJ=KHALO,1,-1
      ZY(JJ) = ZY(JJ+1) - ZDY
      ZY(NJMAX+2*KHALO-JJ+1) = ZY(NJMAX+2*KHALO-JJ) + ZDY
    END DO
  END IF

  DO JJ=1,NJMAX+2*KHALO
    DO JI=1,NIMAX+2*KHALO
      IINDEX = JI+(JJ-1)*(NIMAX+2*KHALO)
      PFIELD_SPLIT(IINDEX) = ZY(JJ)
    END DO
  END DO
  DEALLOCATE(ZY)

END IF  
!
!-------------------------------------------------------------------------------
END SUBROUTINE SPLIT_GRID_PARAMETERX1_MNH
!
!
!     #############################################################
#ifdef MNH_PARALLEL
      SUBROUTINE SPLIT_GRID_PARAMETERN0_MNH(HGRID,HREC,KHALO,KFIELD,KFIELD_SPLIT)
#else
      SUBROUTINE SPLIT_GRID_PARAMETERN0_MNH(HGRID,HREC,KFIELD,KFIELD_SPLIT)
#endif
!     #############################################################
!
!!****  * - routine to define an integer related to splitted grid
!
!
!
USE MODE_ll
!
#ifdef MNH_PARALLEL
#else
USE MODD_IO_SURF_MNH, ONLY : NHALO
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=10), INTENT(IN) :: HGRID        ! grid type
CHARACTER(LEN=6),  INTENT(IN) :: HREC         ! name of the parameter
#ifdef MNH_PARALLEL
INTEGER,           INTENT(IN) :: KHALO        ! size of the Halo
#endif
INTEGER,           INTENT(IN) :: KFIELD       ! integer scalar for complete grid
INTEGER,           INTENT(OUT):: KFIELD_SPLIT ! integer scalar for splitted grid
!*      0.2   Declarations of local variables
!
INTEGER :: IIB, IIE, IJB, IJE
#ifdef MNH_PARALLEL
#else
INTEGER :: KHALO
!
KHALO = NHALO
#endif
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
IF (HREC=='IMAX') KFIELD_SPLIT = IIE-IIB+1 + 2*KHALO
IF (HREC=='JMAX') KFIELD_SPLIT = IJE-IJB+1 + 2*KHALO
!
!-------------------------------------------------------------------------------
END SUBROUTINE SPLIT_GRID_PARAMETERN0_MNH
