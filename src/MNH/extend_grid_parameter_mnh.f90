!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #############################################################
      SUBROUTINE EXTEND_GRID_PARAMETERX1_MNH(HGRID,HREC,KDIM,KSIZE,KIMAX,KJMAX,PFIELD,PFIELD_EXTEND)
!     #############################################################
!
!!****  * - routine to extend a real splitted array on SURFEX halo
!
!    Author
!  M.Moge  01/03/2015 
! 
!    Modifications
!  06/2015   (M.Moge) using MPI_BCAST to have the same space step on all processes + calling UPDATE_NHALO1D
!  07/2015   (M.Moge) initializing ZY and ZY to zero
!  08/2015   (M.Moge) bug fix in the call to UPDATE_NHALO1D : IIMAX_ll instead of IJMAX_ll
!
USE MODD_IO_SURF_MNH, ONLY : NHALO
USE MODD_VAR_ll, ONLY : NPROC, IP, MPI_PRECISION, NMNH_COMM_WORLD
USE MODD_MPIF
USE MODE_TOOLS_ll, ONLY : INTERSECTION, LWEST_ll, LEAST_ll, LNORTH_ll, LSOUTH_ll
USE MODI_UPDATE_NHALO1D
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=10),      INTENT(IN) :: HGRID       ! grid type
CHARACTER(LEN=6),       INTENT(IN) :: HREC        ! name of the parameter
INTEGER,                INTENT(IN) :: KDIM        ! size of PFIELD
INTEGER,                INTENT(IN) :: KSIZE       ! size of PFIELD_EXTEND
INTEGER,                INTENT(IN) :: KIMAX    !(local) dimension of the domain - X direction
INTEGER,                INTENT(IN) :: KJMAX    !(local) dimension of the domain - Y direction
REAL, DIMENSION(KDIM ), INTENT(IN) :: PFIELD      ! real field for complete grid
REAL, DIMENSION(KSIZE), INTENT(OUT):: PFIELD_EXTEND! real field for splitted grid
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
INTEGER                         :: IINFO_ll
INTEGER :: IIMAX_ll, IJMAX_ll
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_OR_ll ('B',IXOR,IYOR)
NIMAX=IIE-IIB+1
NJMAX=IJE-IJB+1
!
CALL GET_GLOBALDIMS_ll(IIMAX_ll, IJMAX_ll)
!
IF (HREC=='XX' .OR. HREC=='DX') THEN
  ALLOCATE(ZX(NIMAX+2*NHALO))
  ZX(:) = 0.
  ZX(1+NHALO:NIMAX+NHALO) = PFIELD(1:NIMAX)
  IF (HREC=='DX') THEN
    ZDX = PFIELD(1)
    DO JI=NHALO,1,-1
      ZX(JI) = ZDX
      ZX(NIMAX+2*NHALO-JI+1) = ZDX
    END DO
  ELSE IF (HREC=='XX') THEN
    ! UPDATE_NHALO1D should be called after the extrapolation on the outer boundaries
    ! to avoid any problem when the halo is larger than the neighbouring subdomain
    ! but this case is not handled and will lead to an ABORT in UPDATE_NHALO1D
    CALL UPDATE_NHALO1D( NHALO, ZX, IIMAX_ll, IJMAX_ll, IIB+IXOR-1,IIE+IXOR-1,IJB+IYOR-1,IJE+IYOR-1, HREC)
    ZDX = 0.
    IF (NIMAX>1) ZDX = PFIELD(2) - PFIELD(1)
    IF (NIMAX==1) ZDX = PFIELD(1) ! in 1D conf, one assumes that grid
                                     ! is located between X=DX/2 and X=3DX/2
    CALL MPI_BCAST(ZDX, 1, MPI_PRECISION, 0, NMNH_COMM_WORLD, IINFO_ll)
    IF( LWEST_ll() ) THEN
      DO JI=NHALO,1,-1
	ZX(JI) = ZX(JI+1) - ZDX
      END DO
    ENDIF
    IF( LEAST_ll() ) THEN
      DO JI=NHALO,1,-1
	ZX(NIMAX+2*NHALO-JI+1) = ZX(NIMAX+2*NHALO-JI) + ZDX
      END DO
    ENDIF
  END IF

!
  DO JJ=1,NJMAX+2*NHALO
    DO JI=1,NIMAX+2*NHALO
      IINDEX = JI+(JJ-1)*(NIMAX+2*NHALO)
      PFIELD_EXTEND(IINDEX) = ZX(JI)
    END DO
  END DO
  DEALLOCATE(ZX)
  
ELSEIF (HREC=='YY' .OR. HREC=='DY') THEN
  ALLOCATE(ZY(NJMAX+2*NHALO))
  ZY(:) = 0.
  DO JJ=1+NHALO,NJMAX+NHALO
    ZY(JJ) = PFIELD(1+(JJ-NHALO-1)*KIMAX)
  END DO
  IF (HREC=='DY') THEN
    ZDY = PFIELD(1)
    DO JJ=NHALO,1,-1
      ZY(JJ) = ZDY
      ZY(NJMAX+2*NHALO-JJ+1) = ZDY
    END DO
  ELSE IF (HREC=='YY') THEN
    ! UPDATE_NHALO1D should be called after the extrapolation on the outer boundaries
    ! to avoid any problem when the halo is larger than the neighbouring subdomain
    ! but this case is not handled and will lead to an ABORT  in UPDATE_NHALO1D
    CALL UPDATE_NHALO1D( NHALO, ZY, IIMAX_ll, IJMAX_ll, IIB+IXOR-1,IIE+IXOR-1,IJB+IYOR-1,IJE+IYOR-1, HREC)
    ZDY = 0.
    IF (NJMAX>1) ZDY = PFIELD(1+KIMAX) - PFIELD(1)
    IF (NJMAX==1) ZDY = PFIELD(1) ! in 1D or 2D conf, one assumes that grid
                                     ! is located between Y=DY/2 and Y=3DY/2
    CALL MPI_BCAST(ZDY, 1, MPI_PRECISION, 0, NMNH_COMM_WORLD, IINFO_ll)
    IF( LSOUTH_ll() ) THEN
      DO JJ=NHALO,1,-1
	ZY(JJ) = ZY(JJ+1) - ZDY
      END DO
    ENDIF
    IF( LNORTH_ll() ) THEN
      DO JJ=NHALO,1,-1
	ZY(NJMAX+2*NHALO-JJ+1) = ZY(NJMAX+2*NHALO-JJ) + ZDY
      END DO
    ENDIF
  END IF

  DO JJ=1,NJMAX+2*NHALO
    DO JI=1,NIMAX+2*NHALO
      IINDEX = JI+(JJ-1)*(NIMAX+2*NHALO)
      PFIELD_EXTEND(IINDEX) = ZY(JJ)
    END DO
  END DO
  DEALLOCATE(ZY)

END IF  
!
!-------------------------------------------------------------------------------
END SUBROUTINE EXTEND_GRID_PARAMETERX1_MNH
!
!
!     #############################################################
      SUBROUTINE EXTEND_GRID_PARAMETERN0_MNH(HGRID,HREC,KFIELD,KFIELD_EXTEND)
!     #############################################################
!
!!****  * - routine to "extend" an integer related to splitted grid on SURFEX halo
!
!
!
USE MODE_ll
!
USE MODD_IO_SURF_MNH, ONLY : NHALO
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=10), INTENT(IN) :: HGRID        ! grid type
CHARACTER(LEN=6),  INTENT(IN) :: HREC         ! name of the parameter
INTEGER,           INTENT(IN) :: KFIELD       ! integer scalar for complete grid
INTEGER,           INTENT(OUT):: KFIELD_EXTEND ! integer scalar for splitted grid
!*      0.2   Declarations of local variables
!
INTEGER :: IIB, IIE, IJB, IJE
!-------------------------------------------------------------------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
IF (HREC=='IMAX') KFIELD_EXTEND = IIE-IIB+1 + 2*NHALO
IF (HREC=='JMAX') KFIELD_EXTEND = IJE-IJB+1 + 2*NHALO
!
!-------------------------------------------------------------------------------
END SUBROUTINE EXTEND_GRID_PARAMETERN0_MNH
