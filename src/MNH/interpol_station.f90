!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!###########################
MODULE MODI_INTERPOL_STATION
!###########################
INTERFACE INTERPOL_STATION
      SUBROUTINE INTERPOL_STATION_3D(PGPS,PXHAT,PYHAT,KI,KJ,&
                               PXHAT_GPS,PYHAT_GPS,PSTAT_GPS)
REAL, DIMENSION(:,:,:),  INTENT(IN) :: PGPS ! GPS parameter array defined on the whole domain
REAL, DIMENSION(:),    INTENT(IN) :: PXHAT !  x positions of the MESO-NH grid 
REAL, DIMENSION(:),    INTENT(IN) :: PYHAT !  y positions of the MESO-NH grid 
INTEGER,               INTENT(IN) :: KI ! the closest south-west grid point 
INTEGER,               INTENT(IN) :: KJ ! the closest south-west grid point  
REAL,                  INTENT(IN) :: PXHAT_GPS ! x positions of the GPS station
REAL,                  INTENT(IN) :: PYHAT_GPS ! y positions of the GPS station
REAL,   DIMENSION(:), INTENT(OUT) :: PSTAT_GPS !  value of the GPS parameter at the station
END SUBROUTINE INTERPOL_STATION_3D
! 
      SUBROUTINE INTERPOL_STATION_2D(PGPS,PXHAT,PYHAT,KI,KJ,&
                               PXHAT_GPS,PYHAT_GPS,PSTAT_GPS)
REAL, DIMENSION(:,:),  INTENT(IN) :: PGPS ! GPS parameter array defined on the whole domain
REAL, DIMENSION(:),    INTENT(IN) :: PXHAT !  x positions of the MESO-NH grid 
REAL, DIMENSION(:),    INTENT(IN) :: PYHAT !  y positions of the MESO-NH grid 
INTEGER,               INTENT(IN) :: KI ! the closest south-west grid point 
INTEGER,               INTENT(IN) :: KJ ! the closest south-west grid point  
REAL,                  INTENT(IN) :: PXHAT_GPS ! x positions of the GPS station
REAL,                  INTENT(IN) :: PYHAT_GPS ! y positions of the GPS station
REAL,                 INTENT(OUT) :: PSTAT_GPS !  value of the GPS parameter at the station
END SUBROUTINE INTERPOL_STATION_2D
!
END INTERFACE
END MODULE MODI_INTERPOL_STATION
!##############################
MODULE MODI_INTERPOL_STATION_3D
!##############################
INTERFACE
!
      SUBROUTINE  INTERPOL_STATION_3D(PGPS,PXHAT,PYHAT,KI,KJ,&
                                PXHAT_GPS,PYHAT_GPS,PSTAT_GPS)
REAL, DIMENSION(:,:,:),  INTENT(IN) :: PGPS ! GPS parameter array defined on the whole domain
REAL, DIMENSION(:),    INTENT(IN) :: PXHAT !  x positions of the MESO-NH grid 
REAL, DIMENSION(:),    INTENT(IN) :: PYHAT !  y positions of the MESO-NH grid 
INTEGER,               INTENT(IN) :: KI ! the closest south-west grid point 
INTEGER,               INTENT(IN) :: KJ ! the closest south-west grid point  
REAL,                  INTENT(IN) :: PXHAT_GPS ! x positions of the GPS station
REAL,                  INTENT(IN) :: PYHAT_GPS ! y positions of the GPS station
REAL,   DIMENSION(:), INTENT(OUT) :: PSTAT_GPS !  value of the GPS parameter at the station
END SUBROUTINE INTERPOL_STATION_3D
END INTERFACE
!
END MODULE MODI_INTERPOL_STATION_3D
!
!     #######################################################
      SUBROUTINE INTERPOL_STATION_3D(PGPS,PXHAT,PYHAT,KI,KJ,&
                               PXHAT_GPS,PYHAT_GPS,PSTAT_GPS)
!     #######################################################
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL, DIMENSION(:,:,:),  INTENT(IN) :: PGPS ! GPS parameter array defined on the whole domain
REAL, DIMENSION(:),    INTENT(IN) :: PXHAT !  x positions of the MESO-NH grid 
REAL, DIMENSION(:),    INTENT(IN) :: PYHAT !  y positions of the MESO-NH grid 
INTEGER,               INTENT(IN) :: KI ! the closest south-west grid point 
INTEGER,               INTENT(IN) :: KJ ! the closest south-west grid point  
REAL,                  INTENT(IN) :: PXHAT_GPS ! x positions of the GPS station
REAL,                  INTENT(IN) :: PYHAT_GPS ! y positions of the GPS station
REAL,   DIMENSION(:), INTENT(OUT) :: PSTAT_GPS !  value of the GPS parameter at the station
INTEGER                           :: IKL ,JK     !
!
!*       0.2   Declarations of local variables :
!
REAL :: ZALPHA,ZBETA
!-------------------------------------------------------------------------------
!
IKL=SIZE (PGPS,3)
!
DO JK=1,IKL  
  ZALPHA=PGPS(KI,KJ,JK)+ (PGPS(KI,KJ+1,JK)-PGPS(KI,KJ,JK))* &
       ( (PYHAT_GPS-PYHAT(KJ))/(PYHAT(KJ+1)-PYHAT(KJ)))
  ZBETA=PGPS(KI+1,KJ,JK)+ (PGPS(KI+1,KJ+1,JK)-PGPS(KI+1,KJ,JK))* &
       ( (PYHAT_GPS-PYHAT(KJ))/(PYHAT(KJ+1)-PYHAT(KJ)))
  PSTAT_GPS(JK)=ZALPHA + (ZBETA-ZALPHA) * ((PXHAT_GPS-PXHAT(KI))/(PXHAT(KI+1)-PXHAT(KI)))
END DO
! 
END SUBROUTINE INTERPOL_STATION_3D
!
!     #######################################################
      SUBROUTINE INTERPOL_STATION_2D(PGPS,PXHAT,PYHAT,KI,KJ,&
                               PXHAT_GPS,PYHAT_GPS,PSTAT_GPS)
!     #######################################################
!
!-------------------------------------------------------------------------------!
!*       0.1   Declarations of dummy arguments :
!
USE MODI_INTERPOL_STATION_3D
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:),  INTENT(IN) :: PGPS ! GPS parameter array defined on the whole domain
REAL, DIMENSION(:),    INTENT(IN) :: PXHAT !  x positions of the MESO-NH grid 
REAL, DIMENSION(:),    INTENT(IN) :: PYHAT !  y positions of the MESO-NH grid 
INTEGER,               INTENT(IN) :: KI ! the closest south-west grid point 
INTEGER,               INTENT(IN) :: KJ ! the closest south-west grid point  
REAL,                  INTENT(IN) :: PXHAT_GPS ! x positions of the GPS station
REAL,                  INTENT(IN) :: PYHAT_GPS ! y positions of the GPS station
REAL,    INTENT(OUT) :: PSTAT_GPS !  value of the GPS parameter at the station
REAL, DIMENSION(SIZE(PGPS,1), SIZE(PGPS,2),1) :: ZFIELD1
REAL, DIMENSION(1) :: ZFIELD2
!
!-------------------------------------------------------------------------------!
ZFIELD1(:,:,1)=PGPS(:,:)
CALL INTERPOL_STATION_3D(ZFIELD1,PXHAT,PYHAT,KI,KJ,PXHAT_GPS,PYHAT_GPS,ZFIELD2)
PSTAT_GPS=ZFIELD2(1)
!
END SUBROUTINE INTERPOL_STATION_2D
