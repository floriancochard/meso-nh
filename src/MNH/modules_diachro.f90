!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_MENU_DIACHRO
!     #########################
!
INTERFACE
!
SUBROUTINE MENU_DIACHRO(TPDIAFILE,HGROUP)
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),  INTENT(IN) :: TPDIAFILE    ! file to write
CHARACTER(LEN=*), INTENT(IN) :: HGROUP
!
END SUBROUTINE MENU_DIACHRO
!
END INTERFACE
!
END MODULE MODI_MENU_DIACHRO
!-----------------------------------------------------------------
!     #####################################
      MODULE MODI_WRITE_LFIFMN_FORDIACHRO_n
!     #####################################
!
INTERFACE
!
SUBROUTINE WRITE_LFIFMN_FORDIACHRO_n(TPFILE)
USE MODD_IO_ll, ONLY: TFILEDATA
TYPE(TFILEDATA),INTENT(IN) :: TPFILE
END SUBROUTINE WRITE_LFIFMN_FORDIACHRO_n
!
END INTERFACE
!
END MODULE MODI_WRITE_LFIFMN_FORDIACHRO_n
!     #########################
      MODULE MODI_WRITE_DIACHRO
!     #########################
!
INTERFACE
!
SUBROUTINE WRITE_DIACHRO(TPDIAFILE,TPLUOUTDIA,HGROUP,HTYPE,          &
      KGRID,PDATIME,PVAR,PTRAJT,                                     &
      HTITRE,HUNITE,HCOMMENT,OICP,OJCP,OKCP,KIL,KIH,KJL,KJH,KKL,KKH, &
      PTRAJX,PTRAJY,PTRAJZ,PMASK)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),              INTENT(IN)          :: TPDIAFILE    ! file to write
TYPE(TFILEDATA),              INTENT(IN)          :: TPLUOUTDIA
CHARACTER(LEN=*),             INTENT(IN)          :: HGROUP, HTYPE
INTEGER,DIMENSION(:),         INTENT(IN)          :: KGRID
REAL,DIMENSION(:,:),          INTENT(IN)          :: PDATIME
REAL,DIMENSION(:,:,:,:,:,:),  INTENT(IN)          :: PVAR
REAL,DIMENSION(:,:),          INTENT(IN)          :: PTRAJT
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN)          :: HTITRE, HUNITE, HCOMMENT
LOGICAL,                      INTENT(IN),OPTIONAL :: OICP, OJCP, OKCP
INTEGER,                      INTENT(IN),OPTIONAL :: KIL, KIH
INTEGER,                      INTENT(IN),OPTIONAL :: KJL, KJH
INTEGER,                      INTENT(IN),OPTIONAL :: KKL, KKH
REAL,DIMENSION(:,:,:),        INTENT(IN),OPTIONAL :: PTRAJX
REAL,DIMENSION(:,:,:),        INTENT(IN),OPTIONAL :: PTRAJY
REAL,DIMENSION(:,:,:),        INTENT(IN),OPTIONAL :: PTRAJZ
REAL,DIMENSION(:,:,:,:,:,:),  INTENT(IN),OPTIONAL :: PMASK
!
END SUBROUTINE WRITE_DIACHRO
!
END INTERFACE
!
END MODULE MODI_WRITE_DIACHRO
