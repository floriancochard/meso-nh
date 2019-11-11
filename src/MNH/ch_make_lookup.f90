!MNH_LIC Copyright 1999-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ######################
       PROGRAM CH_MAKE_LOOKUP
!      ######################
!!
!!***  *CH_MAKE_LOOKUP*
!!
!!    PURPOSE
!!    -------
!      create a lookup table for J-Values used by MesoNH-chemistry
!!
!!**  METHOD
!!    ------
!!      The radiative transfer model TUV is called for a dicrete number
!!    of times. The lookup table is fixed for one latitude, longitude, date, 
!!    surface albedo and ozone column dobson. TUV itself requires
!!    a number of input files that will reside in directories
!!    DATAX, DATA0 and DATA4 (to be created at runtime by tar xvf TUVDATA.tar).
!!    The program may be piloted via a namelistfile LOOKUP1.nam. 
!!    The format of the lookup table is as follows:
!!
!!==============================================================================
!! COMMENT LINE
!! NPHOTO  (NUMBER OF REACTIONS)
!! NLEVEL, DZ  (NUMBER OF VERTICAL POINTS, HEIGHT INCREMENT IN M)
!! NTIME,  DT  (NUMBER OF TEMPORAL POINTS, TIME INCREMENT IN H)
!! REACTION 1
!!  ...
!! REACTION NPHOTO
!! FORTRAN FORMAT FOR THE FOLLOWING DATA
!! HEIGHT1  TIME1  J1 ... JN 
!! HEIGHT1  TIME2  J1 ... JN 
!! ...
!! HEIGHT1  TIMEN  J1 ... JN 
!! HEIGHT2  TIME1  J1 ... JN 
!! ...
!! HEIGHTM  TIMEN  J1 ... JN 
!!==============================================================================
!!
!!    REFERENCES
!!    ----------
!!    MesoNH-chemistry book 3
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 01/03/99
!!    Philippe Wautelet: 10/01/2019: use newunit argument to open files
!!
!!    EXTERNAL
!!    --------
!!    TUV39.f (Fortran 77 code from S. Madronich)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    None
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
!
IMPLICIT NONE
!
REAL    :: ALAT, ALONG    ! LATITUDE AND LONGITUDE
REAL    :: ALBNEW, DOBNEW ! SURFACE ALBEDO AND O3 COLUMN DOBSON
INTEGER :: IDATE          ! DATE IN FORMAT YYMMDD
!
INTEGER, PARAMETER :: NTIME  = 96 + 1 ! TEMPORAL DISCRETIZATION
INTEGER, PARAMETER :: NLEVEL = 30 + 1 ! VERTICAL DISCRETIZATION
REAL, DIMENSION(NLEVEL) :: AZ, LWC
REAL, DIMENSION(NTIME)  :: ATIME
REAL :: DZ, DT
REAL :: ZMAX = 30E3       ! MAXIMUM HEIGHT FOR WHICH J-VALUES WILL BE COMPUTED
!
!      J VALUE STORAGE
INTEGER, PARAMETER                    :: NJOUT = 21
REAL, DIMENSION(NLEVEL,NJOUT)         :: JOUT
REAL, DIMENSION(NLEVEL, NJOUT, NTIME) :: JDATA
CHARACTER*40, DIMENSION(NJOUT)        :: JLABELOUT
!
CHARACTER*120 :: HEADDER
REAL          :: UT
INTEGER       :: ILU ! unit number for IO
INTEGER       :: I, J, K, NJIO
CHARACTER*40  :: YFMT = '(2F11.2,5E11.4/99(7E11.4/))'
!
! NAMELIST for options
NAMELIST /NAM_TUV/ ALAT, ALONG, IDATE, ALBNEW, DOBNEW
!
!
!*       1.   INITIALISATION
!        -------------------
!
!      initialize az and atime
DZ = ZMAX / FLOAT(NLEVEL - 1)
DO J = 1, NLEVEL
  AZ(J) = FLOAT(J-1) * DZ
  LWC(J)= 0.0
ENDDO
DT = 24.00 / FLOAT(NTIME - 1)
DO I = 1, NTIME
  ATIME(I) = FLOAT(I-1) * DT
ENDDO
!
!      initialize default values
ALAT = 45.
ALONG = 0.
IDATE = 970621
ALBNEW = -1.
DOBNEW = -1.
!
!      read the namelist
WRITE(*,*) "TRYING TO OPEN FILE LOOKUP1.nam ..."
OPEN(NEWUNIT=ILU,FILE="LOOKUP1.nam",STATUS="OLD",FORM="FORMATTED", ERR=100)
READ(UNIT=ILU,NML=NAM_TUV)
WRITE(*,*) "NAMELIST NAM_TUV INITIALIZED TO:"
WRITE(*,NAM_TUV)
CLOSE(UNIT=ILU)
GOTO 200
!
100    CONTINUE
!      read from input if no namelist file exists
WRITE(*,*) "NO FILE LOOKUP1.NAM FOUND, TRYING INTERACTIVE ..."
PRINT *, "ENTER LATITUDE (+N, -S)"
READ(*,*) ALAT
PRINT *, "ENTER LONGITUDE (+E, -W)"
READ(*,*) ALONG
PRINT *, "ENTER DATE (YYMMDD)"
READ(*,*) IDATE
PRINT *, "ENTER SURFACE ALBEDO (IF <0 -->DATAX/ALBEDO.DAT)"
READ(*,*) ALBNEW
PRINT *, "ENTER TOTAL COLUMN DOBSON (IF <0 NO SCALING)"
READ(*,*) DOBNEW
!
200    CONTINUE
!
WRITE(HEADDER,'(A,F6.1,A,F6.1,A,I6,A,F6.3,A,F8.2)') &
      "TUV39,LAT=", ALAT, ",LON=", ALONG,           &
      ",IDATE(YYMMDD)=", IDATE,                     &
      ",ALB=", ALBNEW,                              &
      ",DOB=", DOBNEW
PRINT *, HEADDER
!
!*       2.   CALL TUV 3.5
!        -----------------
!
DO I = 1, NTIME
  UT = ATIME(I)
  PRINT *, I, "  FROM  ", NTIME, "  DONE ... UT = ", UT
  NJIO = NJOUT
  CALL TUVMAIN(ALAT, ALONG, IDATE, UT, ALBNEW, DOBNEW, NLEVEL, AZ, LWC, &
              NJIO, JOUT, JLABELOUT)
  DO J = 1, NLEVEL
    DO K = 1, NJIO
      JDATA(J, K, I) = JOUT(J,K) 
    ENDDO
  ENDDO
ENDDO
!
!*       3.   WRITE LOOKUP TABLE FILE
!        -----------------
!
OPEN(NEWUNIT=ILU,FILE="PHOTO.TUV39",STATUS="UNKNOWN",FORM="FORMATTED")
WRITE(ILU,'(A)') HEADDER
WRITE(UNIT=ILU, FMT='(I6,A)') NJIO,   "  NUMBER OF PHOTOLYSIS REACTIONS"
WRITE(UNIT=ILU, FMT='(I6, F10.0, A)') NLEVEL, DZ, &
      "  NUMBER OF LEVELS, VERTICAL INCREMENT (M)"
WRITE(UNIT=ILU, FMT='(I6, F10.4, A)') NTIME,  DT, &
      "  NUMBER OF TEMPORAL RECS, TEMPORAL INCREMENT (H)"
WRITE(UNIT=ILU, FMT='(A)') (JLABELOUT(K), K=1, NJIO)
WRITE(UNIT=ILU, FMT='(A)') YFMT
DO J = 1, NLEVEL
  DO I = 1, NTIME
    WRITE(UNIT=ILU, FMT=YFMT) AZ(J), ATIME(I), (JDATA(J,K,I), K=1, NJIO)
  ENDDO
ENDDO
CLOSE(UNIT=ILU)
!
PRINT *, 'Lookup table file PHOTO.TUV39 has been generated.'
PRINT *, 'CH_MAKE_LOOKUP ended correctly.'
!
END PROGRAM CH_MAKE_LOOKUP
