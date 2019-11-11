!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_LOCHEAD
!     ###################
INTERFACE
      SUBROUTINE LOCHEAD(PLATMIN,PLATMAX,PLONMIN,PLONMAX,                    &
                         PGLBLATMIN,PGLBLATMAX,PGLBLONMIN,PGLBLONMAX,        &
                         KGLBNBLAT,KGLBNBLON,PCUTVAL,KSHIFT,KMAX,            &
                         HSAVEDDATAFILE,ODATASAVE                            )
!      
REAL,              INTENT(IN) :: PLATMIN       ! min latitude  of the local field.
REAL,              INTENT(IN) :: PLATMAX       ! max latitude  of the local field.
REAL,              INTENT(IN) :: PLONMIN       ! min longitude of the local field.
REAL,              INTENT(IN) :: PLONMAX       ! min longitude of the local field.
REAL,              INTENT(IN) :: PGLBLATMIN    ! min latitude  of the global file
REAL,              INTENT(IN) :: PGLBLATMAX    ! max latitude  of the global file
REAL,              INTENT(IN) :: PGLBLONMIN    ! min longitude of the global file
REAL,              INTENT(IN) :: PGLBLONMAX    ! max longitude of the global file
INTEGER,           INTENT(IN) :: KGLBNBLAT     ! number of latitude  rows in global file
INTEGER,           INTENT(IN) :: KGLBNBLON     ! number of longitude rows in global file
REAL,              INTENT(IN) :: PCUTVAL       ! special value in data file
INTEGER,           INTENT(OUT):: KSHIFT        ! shift applied to longitude array
INTEGER,           INTENT(OUT):: KMAX          ! maximum index of new longitude 
!                                              ! array in local area
CHARACTER(LEN=28), INTENT(IN) :: HSAVEDDATAFILE! Name of the local field file
LOGICAL,           INTENT(IN) :: ODATASAVE     ! flag to save data on the local
!                                              ! field file
!
END SUBROUTINE LOCHEAD
END INTERFACE
END MODULE MODI_LOCHEAD
!
!
!     ################################################################
      SUBROUTINE LOCHEAD(PLATMIN,PLATMAX,PLONMIN,PLONMAX,                    &
                         PGLBLATMIN,PGLBLATMAX,PGLBLONMIN,PGLBLONMAX,        &
                         KGLBNBLAT,KGLBNBLON,PCUTVAL,KSHIFT,KMAX,            &
                         HSAVEDDATAFILE,ODATASAVE                            )
!     ################################################################
!
!!**** *LOCHEAD* writes the head of a the local 'latlon' file.
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson              Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    29/08/95
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_IO_ll,            ONLY : TFILEDATA
!
USE MODE_IO_MANAGE_STRUCT, ONLY : IO_FILE_FIND_BYNAME
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL,              INTENT(IN) :: PLATMIN       ! min latitude  of the local area
REAL,              INTENT(IN) :: PLATMAX       ! max latitude  of the local area
REAL,           INTENT(INOUT) :: PLONMIN       ! min longitude of the local area
REAL,           INTENT(INOUT) :: PLONMAX       ! min longitude of the local area
REAL,              INTENT(IN) :: PGLBLATMIN    ! min latitude  of the global file
REAL,              INTENT(IN) :: PGLBLATMAX    ! max latitude  of the global file
REAL,              INTENT(IN) :: PGLBLONMIN    ! min longitude of the global file
REAL,              INTENT(IN) :: PGLBLONMAX    ! max longitude of the global file
INTEGER,           INTENT(IN) :: KGLBNBLAT     ! number of latitude  rows in global file
INTEGER,           INTENT(IN) :: KGLBNBLON     ! number of longitude rows in global file
REAL,              INTENT(IN) :: PCUTVAL       ! special value in data file
INTEGER,           INTENT(OUT):: KSHIFT        ! shift applied to longitude array
INTEGER,           INTENT(OUT):: KMAX          ! maximum index of new longitude 
!                                              ! array in local area
CHARACTER(LEN=28), INTENT(IN) :: HSAVEDDATAFILE! Name of the local field file
LOGICAL,           INTENT(IN) :: ODATASAVE     ! flag to save data on the local
!                                              ! field file
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER                    :: ISAVE       ! logical unit
INTEGER                    :: IRESP      ! return code
REAL, DIMENSION(KGLBNBLAT) :: Z1         ! latitudes  of global field
REAL, DIMENSION(KGLBNBLON) :: Z2         ! longitudes of global field
REAL                       :: ZDLAT      ! latitude  mesh in the data file
REAL                       :: ZDLON      ! longitude mesh in the data file
REAL                       :: Z1MIN      ! min latitude  of the local file
REAL                       :: Z1MAX      ! max latitude  of the local file
REAL                       :: Z2MIN      ! min longitude of the local file
REAL                       :: Z2MAX      ! max longitude of the local file
INTEGER, DIMENSION(1)      :: INB1       ! number of lines   in local file
INTEGER, DIMENSION(1)      :: INB2       ! number of columns in local file
INTEGER                    :: JLAT       ! loop control
INTEGER                    :: JLON       ! loop control
TYPE(TFILEDATA),POINTER    :: TZFILE
!-------------------------------------------------------------------------------
!
IF (ODATASAVE) THEN
  CALL IO_FILE_FIND_BYNAME(HSAVEDDATAFILE,TZFILE,IRESP)
  ISAVE = TZFILE%NLU
END IF
!
ZDLAT=(PGLBLATMAX-PGLBLATMIN)/KGLBNBLAT
ZDLON=(PGLBLONMAX-PGLBLONMIN)/KGLBNBLON
!
Z1(:)=(/ (PGLBLATMAX-(JLAT-0.5)*ZDLAT, JLAT=1,KGLBNBLAT) /)
Z2(:)=(/ (PGLBLONMIN+(JLON-0.5)*ZDLON, JLON=1,KGLBNBLON) /)
!
IF (MINVAL(Z1)>PLATMAX .OR. MAXVAL(Z1)<PLATMIN) THEN
  KSHIFT=0
  RETURN
END IF
Z1MIN=MINVAL(Z1,MASK=(Z1>PLATMIN))-0.5*ZDLAT
Z1MAX=MAXVAL(Z1,MASK=(Z1<PLATMAX))+0.5*ZDLAT
INB1(:)=MINLOC(Z1,MASK=(Z1>PLATMIN)) &
       -MAXLOC(Z1,MASK=(Z1<PLATMAX)) +1
!
!*     Computations on longitudes, shift of longitudes below PLONMIN
!
IF (  (PLONMAX+NINT((PGLBLONMIN-180.-PLONMIN)/360.)*360.<PGLBLONMIN)   &
  .AND.(PLONMIN+NINT((PGLBLONMIN+180.-PLONMIN)/360.)*360.>PGLBLONMAX)  ) THEN
  KMAX=0
  KSHIFT=0
ELSE 
  IF (PLONMAX+NINT((PGLBLONMIN-180.-PLONMIN)/360.)*360.>PGLBLONMIN) THEN
    Z2(:)=Z2(:)+NINT((PGLBLONMIN-180.-PLONMIN)/360.)*360.
  ELSE
    Z2(:)=Z2(:)+NINT((PGLBLONMIN+180.-PLONMIN)/360.)*360.
  END IF
  KSHIFT=COUNT(Z2(:)<PLONMIN)
  WHERE(Z2<PLONMIN)
    Z2=Z2+360.
  ENDWHERE
  Z2=CSHIFT(Z2,SHIFT=KSHIFT)
  INB2(:)=MAXLOC(Z2,MASK=(Z2<PLONMAX))
  KMAX=INB2(1)
  Z2MIN=MINVAL(Z2,MASK=(Z2>PLONMIN))-0.5*ZDLON
  Z2MAX=MAXVAL(Z2,MASK=(Z2<PLONMAX))+0.5*ZDLON
END IF
!
!-------------------------------------------------------------------------------
!
IF ( (KMAX>0) .AND. (INB1(1)>0) .AND. (ODATASAVE) ) THEN
  WRITE(ISAVE,*) 'local file ',HSAVEDDATAFILE
  WRITE(ISAVE,'(A8,F13.8)') 'nodata: ',PCUTVAL
  WRITE(ISAVE,'(A7,F13.8)') 'north: ',Z1MAX
  WRITE(ISAVE,'(A7,F13.8)') 'south: ',Z1MIN
  WRITE(ISAVE,'(A7,F13.8)') 'east: ', Z2MAX
  WRITE(ISAVE,'(A7,F13.8)') 'west: ', Z2MIN
  WRITE(ISAVE,'(A6,I7)') 'rows: ', INB1(1)
  WRITE(ISAVE,'(A6,I7)') 'cols: ', INB2(1)
END IF
!
!-------------------------------------------------------------------------------
END SUBROUTINE LOCHEAD
