!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.

!     ########################
      MODULE MODI_READ_ASC_LATPRESS
!     ########################
INTERFACE
      SUBROUTINE READ_ASC_LATPRESS(HFILENAME,KLEV,PLAT,PLEV,PTHFRC,          &  
                                   PRVFRC)
      
!
CHARACTER(LEN=28), INTENT(IN) :: HFILENAME     ! Name of the field file.
INTEGER , INTENT(IN)      :: KLEV
REAL , DIMENSION(:)   , INTENT(OUT)      :: PLAT
REAL , DIMENSION(:)   , INTENT(OUT)      :: PLEV
REAL , DIMENSION(:,:) , INTENT(OUT)      :: PTHFRC
REAL , DIMENSION(:,:) , INTENT(OUT)      :: PRVFRC
!
!
END SUBROUTINE READ_ASC_LATPRESS
END INTERFACE
END MODULE MODI_READ_ASC_LATPRESS
!
!
!     ##############################################################
      SUBROUTINE READ_ASC_LATPRESS(HFILENAME, KLEV,PLAT,PLEV,PTHFRC,        &  
                                   PRVFRC)
!     ##############################################################
!
!!**** *READ_ASCLLV* reads a binary latlonvalue file and call treatment 
!!                   subroutine
!!
!!    PURPOSE
!!    -------
!!    Reads ascii files to set advective and relaxation forcing values.
!!    NB : Files must be lat,lev, th_frc, rv_frc
!!       
!!    AUTHOR
!!    ------
!!    P. Peyrille 
!!
!!    MODIFICATION
!!    ------------
!!
!!  Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!


!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=28), INTENT(IN) :: HFILENAME     ! Name of the field file.
INTEGER , INTENT(IN)      :: KLEV
REAL , DIMENSION(:)   , INTENT(OUT)      :: PLAT
REAL , DIMENSION(:)   , INTENT(OUT)      :: PLEV
REAL , DIMENSION(:,:) , INTENT(OUT)      :: PTHFRC
REAL , DIMENSION(:,:) , INTENT(OUT)      :: PRVFRC
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER      :: KUNIT                       ! logical unit
!
INTEGER      :: ILUOUT                     ! output listing
INTEGER ::  JI,JK
INTEGER :: IIB,IIE,IJB,IJE
!----------------------------------------------------------------------------
!
!*    1.      Open the file
!             -------------
!
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)

!
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
OPEN(NEWUNIT=KUNIT,FILE=HFILENAME)
!
!*    3.     Reading of a data point
!            -----------------------
!
DO JI=IIB,IIE
  DO JK=1,KLEV
    READ(UNIT=KUNIT,FMT=*,END=99) PLAT(JI),PLEV(JK),PTHFRC(JI,JK),PRVFRC(JI,JK)
  END DO
END DO
!
!----------------------------------------------------------------------------

!
!*    8.    Closing of the data file
!           ------------------------
!
99 CLOSE(UNIT=KUNIT)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_ASC_LATPRESS
