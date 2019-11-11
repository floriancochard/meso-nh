!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 diag 2006/11/17 14:04:01
!-----------------------------------------------------------------
!    ########################
     MODULE MODI_MAKE_RADSAT
!    ########################
!
INTERFACE MAKE_RADSAT
!
     SUBROUTINE MAKE_RADSAT(KYEARF, KMONTHF, KDAYF, PSECF, &
                            KGEO, KLON, PRADB, PRADF)
!
INTEGER, INTENT(IN) :: KYEARF  ! year of Final date
INTEGER, INTENT(IN) :: KMONTHF ! month of Final date
INTEGER, INTENT(IN) :: KDAYF   ! day of Final date
REAL,    INTENT(IN) :: PSECF   ! number of seconds since date at 00 UTC 
INTEGER, INTENT(IN) :: KGEO  
INTEGER, INTENT(IN) :: KLON
REAL, DIMENSION(:,:),  INTENT(IN)  :: PRADB
REAL, DIMENSION(:,:),  INTENT(OUT) :: PRADF
!
END SUBROUTINE MAKE_RADSAT
END INTERFACE
END MODULE MODI_MAKE_RADSAT
!     #########################################################################
!     #########################################################################
      SUBROUTINE MAKE_RADSAT(KYEARF, KMONTHF, KDAYF, PSECF, &
                             KGEO, KLON, PRADB, PRADF)
!     #########################################################################
!
!!****  *MAKE_RADSAT* -
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to convert radiances to brightness temperatures
!     taking into account the filter function of each channel (IR and WV). 
!     
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      J.-P. Chaboureau       *L.A.*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    29/03/00
!       Modified    15/04/03 (J-P Chaboureau) use INDSAT keyword for METEOSAT 63deg service
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_RAD_TRANSF
!
IMPLICIT NONE
!
!*       0.1   DECLARATIONS OF DUMMY ARGUMENTS :
!
INTEGER, INTENT(IN) :: KYEARF  ! year of Final date
INTEGER, INTENT(IN) :: KMONTHF ! month of Final date
INTEGER, INTENT(IN) :: KDAYF   ! day of Final date
REAL,    INTENT(IN) :: PSECF   ! number of seconds since date at 00 UTC 
INTEGER, INTENT(IN) :: KGEO  
INTEGER, INTENT(IN) :: KLON
REAL, DIMENSION(:,:), INTENT(IN)  :: PRADB
REAL, DIMENSION(:,:), INTENT(OUT) :: PRADF
!
!*       0.2   DECLARATIONS OF LOCAL VARIABLES
!
INTEGER :: JI,JJ,JK,JKRAD ! loop indexes
INTEGER, DIMENSION(JPGEOST,JPSAT)         :: IBEGDATETOT
INTEGER, DIMENSION(JPGEOST,JPSAT)         :: IENDDATETOT
INTEGER                                   :: IDATE, IFLAG
REAL, DIMENSION(2,JPCAN,JPGEOST,JPSAT)    :: ZCOEFTOT
REAL, DIMENSION(JPV10,JPCAN,JPGEOST,JPSAT) :: ZFILTTOT
REAL, DIMENSION(KLON,JPWVINT)         :: ZFILT
REAL, DIMENSION(KLON,2,JPCAN)         :: ZCOEF
REAL                                  :: ZRAD
!
!----------------------------------------------------------------------------
!
!*       1   GOES EAST SATELLITE
!
! Source : Remy Roca LMD (roca@lmd.polytechnique.fr)
IBEGDATETOT(1,:)  = (/ 0, 0, 0, 0, 0, 0, 0/)
IENDDATETOT(1,:)  = (/ 2099999999, 0, 0, 0, 0, 0, 0/)
ZCOEFTOT(:,:,1,:) = 0.
ZFILTTOT(:,:,1,:) = 0.
!GOES-8 (Mar 95 - Apr 98) 
ZCOEFTOT(:,1,1,1) = (/ 6.80985, -1359.2   /)
ZCOEFTOT(:,2,1,1) = (/ 8.38525, -2140.27  /)
ZFILTTOT(:,1,1,1) = (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &
                       0.308532, 0.977172, 0.985759, 0.845067, 0.00813178/)
ZFILTTOT(:,2,1,1) = (/ 0.00000, 0.00000, 0.801756, 0.811982, 0.714856, &
                       0.00000, 0.00000, 0.00000, 0.00000, 0.00000 /)     

!----------------------------------------------------------------------------
!
!*       2   GOES WEST SATELLITE
!
IBEGDATETOT(2,:)  = 0
IENDDATETOT(2,:)  = 0
ZCOEFTOT(:,:,2,:) = 0.
ZFILTTOT(:,:,2,:) = 0.
!
!----------------------------------------------------------------------------
!
!*       3   G.M.S. SATELLITE
!
IBEGDATETOT(3,:)  = 0
IENDDATETOT(3,:)  = 0
ZCOEFTOT(:,:,3,:) = 0.
ZFILTTOT(:,:,3,:) = 0.
!
!----------------------------------------------------------------------------
!
!*       4   INDSAT SATELLITE
!
IBEGDATETOT(4,:)  = 0
IENDDATETOT(4,:)  = 0
ZCOEFTOT(:,:,4,:) = 0.
ZFILTTOT(:,:,4,:) = 0.
!*       4   METEOSAT 63deg SERVICE INSTEAD OF INDSAT SATELLITE
IBEGDATETOT(4,:)  = (/ 0, 0, 0, 0, 0, 0, 0/)
IENDDATETOT(4,:)  = (/ 2099999999, 0, 0, 0, 0, 0, 0/)
ZCOEFTOT(:,:,4,:) = 0.
ZFILTTOT(:,:,4,:) = 0.
!METEOSAT-5 (1 Apr 99 - present)  same filter as METEOSAT-5 (4 Feb 94 - 13 Feb 97)
ZCOEFTOT(:,1,4,1) = (/ 6.76691, -1279.36 /)
ZCOEFTOT(:,2,4,1) = (/ 9.303  , -2272.82 /)
ZFILTTOT(:,1,4,1) = (/ 0.02279, 0.21925, 0.39500, 0.49198, 0.67889, &
                       0.87578, 0.98972, 0.37100, 0.02501, 0.00000 /)
ZFILTTOT(:,2,4,1) = (/ 0.00000, 0.02495, 0.40000, 0.55000, 0.64000, &
                       0.85000, 0.91741, 1.00000, 0.33000, 0.01000 /)
!
!----------------------------------------------------------------------------
!
!*       5   METEOSAT SATELLITE - IR and WV channels
!
!source EUMETSAT http://www.eumetsat.de/en/mtp/satellites.html#operation
IBEGDATETOT(5,:)  = (/ 0, 0,  1988081100, 1990041900,  &
                       1994020400, 1997021309, 1998060300/)
IENDDATETOT(5,:)  = (/ 0, 1988081100, 1990041900, 1994020400, &
                       1997021309, 1998060300, 2099999999/)
!source EUMETSAT http://www.eumetsat.de/en/area3/mpef/calib.html
!METEOSAT-1 irrelevant
ZCOEFTOT(:,1,5,1) = (/ 0., 0. /)
ZCOEFTOT(:,2,5,1) = (/ 0., 0. /)
ZFILTTOT(:,1:2,5,1) = 0.
!METEOSAT-2 (Jul 83 - Aug 88)
ZCOEFTOT(:,1,5,2) = (/ 6.1401, -1267.0000 /)
ZCOEFTOT(:,2,5,2) = (/ 8.7698, -2180.5000 /)
ZFILTTOT(:,1:2,5,2) = 0.
!METEOSAT-3 (Aug 88 - Jul 94)
ZCOEFTOT(:,1,5,3) = (/ 6.1694, -1262.7000 /)
ZCOEFTOT(:,2,5,3) = (/ 8.8812, -2167.9000 /)
ZFILTTOT(:,1:2,5,3) = 0.
! Source : Remy Roca LMD (roca@lmd.polytechnique.fr)
! Corrige le 7 juillet 2000. Les coefficients de conversion
! radiance-temperature sont dorenavant corrects. Pour plus
! d'informations: Roca, These de l'universite de Paris 7, 2000.
!METEOSAT-4 (19 Apr 90 - 4 Feb 94)
ZCOEFTOT(:,1,5,4) = (/ 7.078, -1390.02 /)
ZCOEFTOT(:,2,5,4) = (/ 9.133, -2262.45 /)
ZFILTTOT(:,1,5,4) = (/ 0.01399, 0.15769, 0.33500, 0.51432, 0.83281, &
                       0.99841, 0.891201, 0.30700, 0.02181, 0.00000 /)
ZFILTTOT(:,2,5,4) = (/ 0.00000, 0.01897, 0.38600, 0.54900, 0.43800, &
                       0.94830, 0.82135, 0.66000, 0.21600, 0.01000 /)
!METEOSAT-5 (4 Feb 94 - 13 Feb 97)
ZCOEFTOT(:,1,5,5) = (/ 6.76691, -1279.36 /)
ZCOEFTOT(:,2,5,5) = (/ 9.303  , -2272.82 /)
ZFILTTOT(:,1,5,5) = (/ 0.02279, 0.21925, 0.39500, 0.49198, 0.67889, &
                       0.87578, 0.98972, 0.37100, 0.02501, 0.00000 /)
ZFILTTOT(:,2,5,5) = (/ 0.00000, 0.02495, 0.40000, 0.55000, 0.64000, &
                       0.85000, 0.91741, 1.00000, 0.33000, 0.01000 /)
!METEOSAT-6 (13 Feb 97 - 3 Jun 98)
ZCOEFTOT(:,1,5,6) = (/ 6.79713, -1274.26 /)
ZCOEFTOT(:,2,5,6) = (/ 9.14985, -2271.18 /)
ZFILTTOT(:,1,5,6) = (/ 0.01919, 0.19007, 0.40200, 0.59991, 0.91245, &
                       0.99741, 0.85842, 0.27600, 0.01881, 0.00000 /)
ZFILTTOT(:,2,5,6) = (/ 0.00000, 0.02022, 0.36000, 0.46900, 0.43900, &
                       0.94880, 0.86402, 0.72500, 0.22400, 0.00000 /)
!METEOSAT-7 (3 Jun 98 - present)
ZCOEFTOT(:,1,5,7) = (/ 7.00093, -1263.14 /)
ZCOEFTOT(:,2,5,7) = (/ 9.22152, -2233.02 /)
ZFILTTOT(:,1,5,7) = (/ 0.10831, 0.48301, 0.75100, 0.85672, 0.98435, &
                       0.99461, 0.91660, 0.41100, 0.01721, 0.00000 /)
ZFILTTOT(:,2,5,7) = (/ 0.00009, 0.08544, 0.61000, 0.68800, 0.80100, &
                       0.74559, 0.84384, 0.73400, 0.14600, 0.00000 /)
!
!----------------------------------------------------------------------------
!
IDATE = KYEARF*1000000 + KMONTHF*10000 + KDAYF*100 + NINT(PSECF/3600.)
IFLAG = 0
DO JI=1,JPSAT
  IF ( IDATE >= IBEGDATETOT(KGEO,JI) .AND. IDATE < IENDDATETOT(KGEO,JI) )THEN
    IFLAG=JI
    PRINT *,'SATELLITE #',KGEO,' FILTER SELECTED FOR #',IFLAG
  END IF
END DO
!
IF (IFLAG==0) THEN
  PRINT *,'NO SATELLITE FILTER FOUND FOR #',KGEO
  PRINT *,' BRIGHTNESS TEMPERATURES will be UNDEFINED'
  PRADF(:,:) = XUNDEF
!
ELSE
  IF( ALL(ZCOEFTOT(:,:, KGEO, IFLAG)==0.) .OR. &
      ALL(ZFILTTOT(:,1, KGEO, IFLAG)==0.) .OR. &
      ALL(ZFILTTOT(:,2, KGEO, IFLAG)==0.)      ) THEN
    PRINT *,'ALL COEFFICIENTS ARE NULL for ',IDATE
    IF(KGEO==5) THEN
      IFLAG=4
      PRINT *,' COEFFICIENTS of METEOSAT-4 are nevertheless used'
      PRINT *,'=> but DO NOT COMPARE to SATELLITE DATA  !!'
    ELSE  
      PRINT *,' BRIGHTNESS TEMPERATURES will be UNDEFINED'
      PRADF(:,:) = XUNDEF
      IFLAG=0
    END IF
  END IF
END IF
!
IF (IFLAG/=0) THEN
  DO JI=1,KLON
    ZCOEF(JI,:,:)   = ZCOEFTOT(:,:, KGEO, IFLAG)
    ZFILT(JI,1:10)  = ZFILTTOT(:,1, KGEO, IFLAG)
    ZFILT(JI,11:20) = ZFILTTOT(:,2, KGEO, IFLAG)
  END DO
  !
  DO JK=1,JPCAN
    DO JI=1,KLON
      ZRAD = 0.
      !Spectral integration
      DO JKRAD=1,JPV10
        JJ   = JKRAD + (JK-1)*JPV10
        ZRAD = ZRAD  + PRADB(JI,JJ)*ZFILT(JI,JJ)
      END DO
      !Conversion from radiance to brightness temperature
      IF (ZRAD > 0.) THEN
        PRADF(JI,JK) = ZCOEF(JI,2,JK) / ( LOG(ZRAD) - ZCOEF(JI,1,JK) )
      ELSE
        PRADF(JI,JK) = XUNDEF
      END IF
    END DO
  END DO
  !
END IF
!
END SUBROUTINE MAKE_RADSAT
