!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    #############################
      SUBROUTINE WRITE_TS1D
!!    #############################

!!    PURPOSE
!!    -------
!!    write quick and dirty 1D time series

!!    METHOD
!!    ------
!!    write at each time step variables in PICASSO.EXE format

!!    REFERENCE
!!    ---------
!!    none

!!    AUTHOR
!!    ------

!!    PURPOSE
!!    -------
!!    write quick and dirty 1D time series

!!    METHOD
!!    ------
!!    write at each time step variables in PICASSO.EXE format

!!    REFERENCE
!!    ---------
!!    none

!!    AUTHOR
!!    ------
!!    K. Suhre + C. Mari

!!    MODIFICATIONS
!!    -------------
!!    Original 12/10/95
!!    Add U, V, W, all water variables and all tracers to timeseries (KS)
!!     27/10/95 KS: add air density
!!     09/02/96 KS: change time to TDTSEG%TIME + ISCOUNT * XTSTEP
!!                 & grid identifier all set to 1 (mass levels)
!!     01/03/96 KS: remove air density
!!                  write variables every 900s to disk, regardless of XTSTEP
!!     02/03/96 KS: write time in hours rather than seconds
!!     24/04/96 KS: add PTIME to parameterlist
!!     29/07/96 KS: use XDUMMY1 variable from MODD_BLANK to define write-tstep
!!     31/07/96 KS: clip values that are smaller than 1E-40 to zero, since
!!                  otherwise strange behavior of PV-WAVE is possible
!!                  add printing of date and time and personalized comment 
!!     07/08/96 KS: add printing of chemical species names and convert to ppt
!!                  print also XRHODREF so that conversiopn ppt to mol/cm3
!!                  could later be made 
!!     20/02/97 KS: add cloud fraction XCLDFR
!!     05/03/97 KS: add writing of JNO2 and JO1D
!!     12/11/99 KS: add parallel open and close
!!     30/11/99 KS: add namelist parameters XCH_TS1D_TSTEP, CCH_TS1D_COMMENT,
!!                  and CCH_TS1D_FILENAME
!!     30/05/04 Tulet : update to be use with diag : now variables XCHEMLAT and
!!                      XCHEMLON could gives some profiles if specified in
!!                      namelist DIAG1.nam. Otherwise treatment is active only in
!!                      1D case
!!     28/03/2018 P. Wautelet: replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!!    EXTERNAL
!!    --------
!!    OPEN_ll and CLOSE_ll  ! attribute a free I/O unit and close it again


!!    IMPLICIT ARGUMENTS
!!    ------------------

USE MODE_DATETIME
USE MODE_FM,              ONLY: IO_FILE_CLOSE_ll,IO_FILE_OPEN_ll
USE MODE_IO_MANAGE_STRUCT,ONLY: IO_FILE_ADD2LIST
USE MODE_IO_ll
USE MODE_GRIDPROJ
USE MODE_ll
!
USE MODD_LUNIT_n,         ONLY: CINIFILE
USE MODD_NSV,             ONLY: NSV,NSV_CHEMBEG,NSV_CHEMEND,  &
                                NSV_AERBEG,NSV_AEREND
USE MODD_CH_AEROSOL,      ONLY: CAERONAMES, LORILAM
USE MODD_DYN_n,           ONLY: XTSTEP       ! time-step of the model
USE MODD_DIM_n,           ONLY: NKMAX        ! # of points in Z of the physical grid
USE MODD_IO_ll,           ONLY: TFILEDATA
USE MODD_PARAMETERS,      ONLY: JPVEXT   ! vertical external points number
USE MODD_GRID,      ONLY: XLATORI,XLONORI
USE MODD_GRID_n,    ONLY: XXHAT,XYHAT,XZZ
USE MODD_TIME_n,    ONLY: TDTCUR     ! Current Time and Date
USE MODD_TIME,      ONLY: TDTEXP     ! Begining Simulation Time and Date

USE MODD_CONF,      ONLY: CPROGRAM,L1D   ! Configuration of models

!!

USE MODD_FIELD_n,         ONLY: XUT,XVT,XWT,&! wind speed at time t 
                                XTHT,      &! theta
                                XRT,       &! water variables
                                XSVT,      &! scalar variables
                                XTKET,     &! Kinetic energy (rho e)
                                XCLDFR      ! cloud fraction

USE MODD_CONF_n,          ONLY: NRR          ! Total number of water variables
USE MODD_REF_n,           ONLY: XRHODREF     ! dry air density of the reference state
USE MODD_CH_M9_n,         ONLY: CNAMES      ! names of chemical species
USE MODD_CH_MNHC_n,       ONLY: LUSECHEM,   &! flag that indicates if chemistry is used
                                XCH_TS1D_TSTEP,   &! time between two call to write_ts1d
                                CCH_TS1D_COMMENT, &! comment for write_ts1d
                                CCH_TS1D_FILENAME  ! filename for write_ts1d files
USE MODD_PARAM_n,         ONLY: CRAD         ! Kind of radiation parameterization
                                             ! 'NONE' if no parameterization
USE MODD_RADIATIONS_n,    ONLY: XDTHRAD ! radiative heating/cooling rate
USE MODD_CH_JVALUES_n,    ONLY: XJVALUES   ! Jvalues and 
USE MODD_CH_INIT_JVALUES, ONLY:JPJVMAX ! their number
USE MODD_PARAMETERS,      ONLY: XUNDEF
USE MODD_DIAG_FLAG,       ONLY: LCHEMDIAG, XCHEMLAT, XCHEMLON
USE MODI_TRANSFER_FILE

IMPLICIT NONE
!!    EXPLICIT ARGUMENTS
!!    ------------------
!     none

!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
REAL  :: ZTIME
LOGICAL, SAVE     :: GSFIRSTCALL = .TRUE.
INTEGER, SAVE     :: ISIO1D                 ! IO-channel 
CHARACTER(LEN=80), SAVE :: YSIO1DDEF        ! name of def-file
CHARACTER(LEN=80), SAVE :: YSIO1DDAT        ! name of dat-file
CHARACTER(LEN=40) :: YFORM = "(E15.8)"      ! data output format
INTEGER, SAVE     :: ISCOUNT = 0            ! timestep counter
INTEGER           :: JK, JL, JN             ! loop counter
INTEGER           :: IIU,IJU                ! domain size
INTEGER           :: NBPROF                 ! number of profile
INTEGER           :: IUTIME                 ! hour in UT
INTEGER, SAVE     :: ISSKIP                 ! # of timesteps to skip
INTEGER           :: IINDEX, JINDEX         ! index of each profile
REAL              :: ZXINDEX, ZYINDEX       ! distance from origin of each profile
REAL              :: ZX, ZY                 ! poisition of each profile
REAL, DIMENSION(SIZE(XCHEMLAT)) ::  ZLAT, ZLON
TYPE(TFILEDATA),POINTER,SAVE :: TZFILE
!
CHARACTER*8  :: YDATE  ! for retrieval of date and time
CHARACTER*10 :: YTIME  ! dito
CHARACTER*13 :: YLATLON  ! dito
CHARACTER*13 :: YCLATLON  ! dito
CHARACTER*4  :: YCYEAR  ! current year
CHARACTER*2  :: YCMONTH ! current month
CHARACTER*2  :: YCDAY   ! current day
CHARACTER*5  :: YCTIME  ! current time
!!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
TZFILE => NULL()
NBPROF = 0

CALL DATETIME_DISTANCE(TDTEXP,TDTCUR,ZTIME)
ZTIME = ZTIME + TDTEXP%TIME
CALL GET_DIM_EXT_ll ('B',IIU,IJU)

IF ((CPROGRAM =='DIAG  ').AND.(LCHEMDIAG)) THEN
  WRITE(YCYEAR,'(I4.4)')  TDTCUR%TDATE%YEAR
  WRITE(YCMONTH,'(I2.2)') TDTCUR%TDATE%MONTH
  WRITE(YCDAY,'(I2.2)')   TDTCUR%TDATE%DAY
  IUTIME = INT(3600*MOD( 24.0+MOD((TDTCUR%TIME)/3600.,24.0),24.0 ))  ! TU 
  WRITE(YCTIME,'(I5.5)')  IUTIME


  DO JN=1,SIZE(XCHEMLAT)
    IF ((XCHEMLAT(JN) /= XUNDEF).AND.(XCHEMLON(JN) /= XUNDEF)) THEN 
      NBPROF = NBPROF + 1
      ZLAT(NBPROF) = XCHEMLAT(JN)
      ZLON(NBPROF) = XCHEMLON(JN)
    END IF
  END DO
END IF

IF (L1D) NBPROF = 1

DO JN=1,NBPROF
  IF ((CPROGRAM =='DIAG  ').AND.(LCHEMDIAG)) THEN
    CALL SM_XYHAT(XLATORI,XLONORI, &
                  ZLAT(JN), ZLON(JN), ZX, ZY )

    ZXINDEX=(ZX-XXHAT(1)) * IIU/(XXHAT(IIU)-XXHAT(1))
    ZYINDEX=(ZY-XYHAT(1)) * IJU/(XYHAT(IJU)-XYHAT(1))
    IINDEX=INT(ZXINDEX)
    JINDEX=INT(ZYINDEX)


    IF ((ABS(ZLAT(JN)) .GE. 10.).AND.(ABS(ZLON(JN)) .GE. 10.)) THEN 
       WRITE(YLATLON,'(F6.3,A1,F6.3)') ABS(ZLAT(JN)),'-',ABS(ZLON(JN))
       WRITE(YCLATLON,'(A13)') TRIM(YLATLON)
    END IF
    IF ((ABS(ZLAT(JN)) .LT. 10.).AND.(ABS(ZLON(JN)) .GE. 10.)) THEN
       WRITE(YLATLON,'(F5.3,A1,F6.3)') ABS(ZLAT(JN)),'-',ABS(ZLON(JN))
       WRITE(YCLATLON,'(A12)') TRIM(YLATLON)
    END IF
    IF ((ABS(ZLAT(JN)) .LT. 10.).AND.(ABS(ZLON(JN)) .LT. 10.)) THEN
       WRITE(YLATLON,'(F5.3,A1,F5.3)') ABS(ZLAT(JN)),'-',ABS(ZLON(JN))
       WRITE(YCLATLON,'(A11)') TRIM(YLATLON)
    END IF
    IF ((ABS(ZLAT(JN)) .GE. 10.).AND.(ABS(ZLON(JN)) .LT. 10.)) THEN
       WRITE(YLATLON,'(F6.3,A1,F5.3)') ABS(ZLAT(JN)),'-',ABS(ZLON(JN))
       WRITE(YCLATLON,'(A12)') TRIM(YLATLON)
    END IF

    CCH_TS1D_FILENAME = ADJUSTR(YCLATLON) // ":" // YCYEAR // YCMONTH // YCDAY // "-" // YCTIME
    CCH_TS1D_COMMENT = "Fichier issu de " // CINIFILE
    YSIO1DDEF = "def." // ADJUSTL(CCH_TS1D_FILENAME)
    YSIO1DDAT = "dat." // ADJUSTL(CCH_TS1D_FILENAME)

  END IF

  IF (L1D) THEN
    YSIO1DDEF = "def." // CCH_TS1D_FILENAME
    YSIO1DDAT = "dat." // CCH_TS1D_FILENAME
    IINDEX= 2
    JINDEX= 2
  END IF

  IF ((IINDEX >= 1).AND.(IINDEX <= IIU).AND.&
      (JINDEX >= 1).AND.(JINDEX <= IJU)) THEN  
  ! write picasso def-file
    IF (GSFIRSTCALL) THEN
      CALL IO_FILE_ADD2LIST(TZFILE,YSIO1DDEF,'TXT','WRITE')
      CALL IO_FILE_OPEN_ll(TZFILE,HPOSITION='REWIND',HSTATUS='NEW')
      ISIO1D = TZFILE%NLU

    ! write comment
      CALL DATE_AND_TIME(YDATE, YTIME)
      WRITE(ISIO1D,*) YDATE, "@", YTIME, ": ", CCH_TS1D_COMMENT

    ! write # of points in Z
    WRITE(ISIO1D,'(I5)') NKMAX 

    ! write # of variables
    WRITE(ISIO1D,'(I5)') 10 + NRR + NSV

    ! write grid identifier (all variables are assumed to be on mass points)
    WRITE(ISIO1D,'(999A1)') ("1", JL = 1, 10 + NRR + NSV)

    ! write variable names
    WRITE(ISIO1D,'(A)') "time (h)" 
    WRITE(ISIO1D,'(A)') "height (m)" 
    WRITE(ISIO1D,'(A)') "U (m/s)"
    WRITE(ISIO1D,'(A)') "V (m/s)"
    WRITE(ISIO1D,'(A)') "W (m/s)"
    WRITE(ISIO1D,'(A)') "THT (K)"
    WRITE(ISIO1D,'(A)') "TKE"
    WRITE(ISIO1D,'(A)') "RHODREF (kg/m3)"
    WRITE(ISIO1D,'(A)') "XDTHRAD radiat. heating/cooling rate"
    WRITE(ISIO1D,'(A)') "XCLDFR (cloud fraction)"
    WRITE(ISIO1D,'(A)') "JNO2 (1/s)"
    WRITE(ISIO1D,'(A)') "JO1D (1/s)"
    DO JL = 1, NRR 
      WRITE(ISIO1D,'(A, I1, A)') "Q", JL, " water variable (kg/kg)"
    ENDDO

    DO JL=1,NSV
      IF (LUSECHEM .AND. JL >= NSV_CHEMBEG .AND. JL <= NSV_CHEMEND) THEN
        WRITE(ISIO1D,'(A)') TRIM(CNAMES(JL-NSV_CHEMBEG+1)) // " (ppt)"
      ELSE IF(LORILAM .AND. JL >= NSV_AERBEG .AND. JL <= NSV_AEREND) THEN
        WRITE(ISIO1D,'(A)') TRIM(CAERONAMES(JL-NSV_AERBEG+1)) // " (ppt)"
      ELSE
        WRITE(ISIO1D,'(A, I3.3, A)') "CHEM", JL, " scalar variable (kg/kg)" 
      END IF
    END DO
  
    CALL IO_FILE_CLOSE_ll(TZFILE)
    TZFILE => NULL()
    CALL TRANSFER_FILE('fujitransfer.x','NIL',YSIO1DDEF)

    ! open picasso dat-file
    CALL IO_FILE_ADD2LIST(TZFILE,YSIO1DDAT,'TXT','WRITE')
    CALL IO_FILE_OPEN_ll(TZFILE,HPOSITION='REWIND',HSTATUS='NEW')
    ISIO1D = TZFILE%NLU

    ! calculate ISSKIP
    IF (L1D) THEN
      ISSKIP = MAX(1,IFIX(XCH_TS1D_TSTEP/XTSTEP))
    END IF

    IF ((CPROGRAM =='DIAG  ').AND.(LCHEMDIAG)) ISSKIP = 1
      PRINT *, "WRITE_TS1D: every ", ISSKIP, " time steps will be stored (!)"
    END IF
  END IF

  IF (MOD(ISCOUNT,ISSKIP).EQ.0) THEN
  ! write a single data set at time t (in hours)
    CALL WRITECLIP ( ZTIME / 3600. )

  ! print height variable
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      CALL WRITECLIP ( XZZ(IINDEX,JINDEX,JK) )
    ENDDO

    ! print U
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      CALL WRITECLIP ( XUT(IINDEX,JINDEX,JK) )
    ENDDO

    ! print V
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      CALL WRITECLIP ( XVT(IINDEX,JINDEX,JK) )
    ENDDO

    ! print W
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      CALL WRITECLIP ( XWT(IINDEX,JINDEX,JK) )
    ENDDO

    ! print THETA
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      CALL WRITECLIP ( XTHT(IINDEX,JINDEX,JK) )
    ENDDO

    ! print TKE
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      CALL WRITECLIP ( XTKET(IINDEX,JINDEX,JK) )
    ENDDO

    ! print RHODREF
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      CALL WRITECLIP ( XRHODREF(IINDEX,JINDEX,JK) )
    ENDDO

    ! print XDTHRAD
    IF (CRAD .NE. "NONE") THEN
      DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
        CALL WRITECLIP ( XDTHRAD(IINDEX,JINDEX,JK) )
      END DO
    ELSE
      DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
        CALL WRITECLIP ( 0.0 )
      END DO
    END IF

    ! print XCLDFR (the field is only allocated when radiation is used)
    IF (CRAD .NE. "NONE") THEN
      DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
        CALL WRITECLIP ( XCLDFR(IINDEX,JINDEX,JK) )
      END DO
    ELSE
      DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
        CALL WRITECLIP ( 0.0 )
      END DO
    END IF

    ! print JNO2
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      ! the grid index of XJVALUES is starting with 1 for the physical grid
      CALL WRITECLIP ( XJVALUES(IINDEX,JINDEX,JK-JPVEXT,2) )
    ENDDO

    ! print JO1D
    DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
      ! the grid index of XJVALUES is starting with 1 for the physical grid
      CALL WRITECLIP ( XJVALUES(IINDEX,JINDEX,JK-JPVEXT,3) )
    ENDDO

    ! print water variables
    DO JL = 1, NRR
      DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
        CALL WRITECLIP ( XRT(IINDEX,JINDEX,JK,JL) )
      ENDDO
    ENDDO

    ! print scalar variables
    DO JL = 1, NSV
      IF (JL>=NSV_CHEMBEG .AND. JL<=NSV_CHEMEND) THEN
        DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
        ! convert ppp to ppt
          CALL WRITECLIP ( XSVT(IINDEX,JINDEX,JK,JL) * 1E12 )
        ENDDO
      ELSE
        DO JK = NKMAX + JPVEXT , JPVEXT + 1, -1
          CALL WRITECLIP ( XSVT(IINDEX,JINDEX,JK,JL) )
        ENDDO
      END IF
    ENDDO

    IF ((CPROGRAM =='DIAG  ').AND.(LCHEMDIAG)) THEN
      CALL IO_FILE_CLOSE_ll(TZFILE)
      TZFILE => NULL()
      CALL TRANSFER_FILE('fujitransfer.x','NIL',YSIO1DDAT)
    END IF
 
    IF (L1D) THEN
      GSFIRSTCALL = .FALSE.
      CALL TRANSFER_FILE('fujitransfer.x','NIL',YSIO1DDAT)
    END IF

  END IF
ENDDO

! increase timestep counter
ISCOUNT = ISCOUNT + 1


CONTAINS

SUBROUTINE WRITECLIP(PVAR)
USE MODD_CST, ONLY : XMNH_TINY
IMPLICIT NONE
REAL, INTENT(IN) :: PVAR

! clip very small numbers to zero, due to problems with PV-WAVE read
IF (ABS(PVAR) <=  XMNH_TINY ) THEN
  WRITE(ISIO1D,FMT=YFORM) 0.0
ELSE
  WRITE(ISIO1D,FMT=YFORM) PVAR
END IF

RETURN
END SUBROUTINE WRITECLIP


END SUBROUTINE WRITE_TS1D
