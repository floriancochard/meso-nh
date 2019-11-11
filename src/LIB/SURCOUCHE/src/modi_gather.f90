!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for CVS information
!-----------------------------------------------------------------
! $Source$
! $Name$ 
! $Revision$ 
! $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------

MODULE MODI_GATHER_ll
!
INTERFACE GATHERALL_FIELD_ll
  SUBROUTINE GATHERALL_X1(HDIR,PSEND,PRECV,KRESP)
  CHARACTER(LEN=*),  INTENT(IN) :: HDIR
  REAL,DIMENSION(:), INTENT(IN) :: PSEND
  REAL,DIMENSION(:), INTENT(INOUT):: PRECV
  INTEGER,           INTENT(INOUT):: KRESP
  END SUBROUTINE GATHERALL_X1
  
  SUBROUTINE GATHERALL_X2(HDIR,PSEND,PRECV,KRESP)
  CHARACTER(LEN=*),    INTENT(IN) :: HDIR
  REAL,DIMENSION(:,:), INTENT(IN) :: PSEND
  REAL,DIMENSION(:,:), INTENT(INOUT):: PRECV
  INTEGER,             INTENT(INOUT):: KRESP
  END SUBROUTINE GATHERALL_X2
  
  SUBROUTINE GATHERALL_X3(HDIR,PSEND,PRECV,KRESP)
  CHARACTER(LEN=*),      INTENT(IN) :: HDIR
  REAL,DIMENSION(:,:,:), INTENT(IN) :: PSEND
  REAL,DIMENSION(:,:,:), INTENT(INOUT):: PRECV
  INTEGER,               INTENT(INOUT):: KRESP
  END SUBROUTINE GATHERALL_X3
  
  SUBROUTINE GATHERALL_N1(HDIR,KSEND,KRECV,KRESP)
  CHARACTER(LEN=*),     INTENT(IN) :: HDIR
  INTEGER,DIMENSION(:), INTENT(IN) :: KSEND
  INTEGER,DIMENSION(:), INTENT(INOUT):: KRECV
  INTEGER,              INTENT(INOUT):: KRESP
  END SUBROUTINE GATHERALL_N1
  
  SUBROUTINE GATHERALL_N2(HDIR,KSEND,KRECV,KRESP)
  CHARACTER(LEN=*),       INTENT(IN) :: HDIR
  INTEGER,DIMENSION(:,:), INTENT(IN) :: KSEND
  INTEGER,DIMENSION(:,:), INTENT(INOUT):: KRECV
  INTEGER,                INTENT(INOUT):: KRESP
  END SUBROUTINE GATHERALL_N2
END INTERFACE
!
END MODULE MODI_GATHER_ll
