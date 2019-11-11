!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############################
      MODULE MODI_COMPUTE_R00
!     ###############################
!
INTERFACE
SUBROUTINE COMPUTE_R00(TPFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! Output file
!
END SUBROUTINE COMPUTE_R00
END INTERFACE
END MODULE MODI_COMPUTE_R00
!
!     ###############################
      SUBROUTINE COMPUTE_R00(TPFILE)
!     ###############################
!
!!**** 
!!
!!    PURPOSE
!!    -------
!     
!!**  METHOD 
!!    ------
!!    
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!     MODD_STO_FILE : CFILES
!!     MODD_GRID1 : XZZ, XXHAT,XYHAT
!!     MODD_LUNIT1: CINIFILE
!!     MODD_FIELD1: XSVM
!!     MODD_CONF : NVERB
!!     MODD_PARAMETERS : NUNDEF
!!     
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!	 F. Gheusi and J. Stein  * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!     J. Stein  Jan. 2001  add supplementary starts and spare some memory
!!     October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                              change of YCOMMENT
!!     Mai 2016 (G.Delautier) replace LG?M by LG?T
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!    
USE MODD_FIELD_n
USE MODD_GRID_n
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_LUNIT_n
USE MODD_GRID_n
USE MODD_STO_FILE
USE MODD_CONF
USE MODD_PARAMETERS
USE MODD_NSV,            ONLY : NSV_LGBEG,NSV_LGEND
!
USE MODI_SHUMAN
!
USE MODD_VAR_ll
!
USE MODE_FIELD
USE MODE_FM
USE MODE_FMWRIT
USE MODE_FMREAD
USE MODE_IO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY : IO_FILE_ADD2LIST
USE MODE_ll
USE MODE_MSG
USE MODD_TYPE_DATE
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! Output file
!
!*       0.2   declarations of local variables
!
INTEGER                            :: IFILECUR,JFILECUR,NIU,NJU,NKU
INTEGER                            :: NFILES,JLOOP
REAL                               :: ZXOR,ZYOR,ZDX,ZDY
REAL                               :: ZSPVAL
REAL, ALLOCATABLE, DIMENSION(:,:,:):: ZX0, ZY0, ZZ0        ! origin of the 
       ! particules colocated with the mesh-grid points read in the file
REAL, ALLOCATABLE, DIMENSION(:,:,:):: ZX00, ZY00, ZZ00, ZZL ! cumulative
       ! origin for more than one restart of the tracers 
REAL, ALLOCATABLE, DIMENSION(:,:,:):: ZTH0          ! same fields 
       ! for Theta as for the coordinates of the origin
REAL, ALLOCATABLE, DIMENSION(:,:,:):: ZRV0          ! same fields 
       ! for Rv as for the coordinates of the origin
REAL, ALLOCATABLE, DIMENSION(:,:,:):: ZWORK1,ZWORK2,ZWORK3       
TYPE(DATE_TIME)                    :: TDTCUR_START
CHARACTER(LEN=24)                  :: YDATE 
INTEGER                            :: IHOUR, IMINUTE
REAL                               :: ZSECOND, ZREMAIN
LOGICAL                            :: GSTART
INTEGER                            :: INBR_START
REAL                               :: ZXMAX,ZYMAX,ZZMAX  ! domain extrema
INTEGER, DIMENSION(100)            :: NBRFILES 
INTEGER                            :: IKU
TYPE(TFIELDDATA)                   :: TZFIELD
TYPE(TFILEDATA),POINTER            :: TZTRACFILE
!
!-------------------------------------------------------------------------------
!
!*       1.0    INITIALIZATION
!               --------------
!
ZSPVAL=-1.E+11
IKU=SIZE(XZHAT)
!
!-------------------------------------------------------------------------------
!
!*       2.0    FIND THE FILE TO BE TREATED AND THE INIT-SV FILES
!               -------------------------------------------------
!
! Search the number of the file to be treated
IFILECUR=0
DO JFILECUR=1,100
  IF (CINIFILE==CFILES(JFILECUR)) THEN
    IFILECUR=JFILECUR
    EXIT
  END IF
END DO
!
IF (IFILECUR==0) THEN
  PRINT*,'PROBLEM WITH THE FOLLOWING FILE: ',CINIFILE
  PRINT*,CFILES
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','COMPUTE_R00','')
ENDIF
!
! Search the number of the files(NFILES), where the Lagrangian tracers 
!have been reinitialized 
NFILES=0
DO JFILECUR=IFILECUR+1,100
  IF (LEN_TRIM(CFILES(JFILECUR)) /= 0) THEN
    NFILES= NFILES +1 
    NBRFILES(NFILES)=JFILECUR       ! contains the number of the files where
                                    ! the Lag. tracers have been restarted
  ENDIF
END DO
!
! compute the number of supplementary cumulative starts
INBR_START=1
DO JLOOP=1,NFILES-1
  IF (NSTART_SUPP(JLOOP)/=NUNDEF .AND. NSTART_SUPP(JLOOP)> IFILECUR ) THEN
    INBR_START=INBR_START+1
  END IF
END DO
!
!-------------------------------------------------------------------------------
!
!*       3.0    ALLOCATIONS OF THE ARRAYS AND CONVERSIONS
!               -----------------------------------------
!
!
NIU=SIZE(XZZ,1)
NJU=SIZE(XZZ,2)
NKU=SIZE(XZZ,3)
!
ALLOCATE(ZX0(NIU,NJU,NKU))
ALLOCATE(ZY0(NIU,NJU,NKU))
ALLOCATE(ZZ0(NIU,NJU,NKU))
ALLOCATE(ZWORK1(NIU,NJU,NKU))
ALLOCATE(ZWORK2(NIU,NJU,NKU))
ALLOCATE(ZWORK3(NIU,NJU,NKU))
ALLOCATE(ZX00(NIU,NJU,NKU))
ALLOCATE(ZY00(NIU,NJU,NKU))
ALLOCATE(ZZ00(NIU,NJU,NKU))
ALLOCATE(ZZL(NIU,NJU,NKU))
ALLOCATE(ZTH0(NIU,NJU,NKU))
ALLOCATE(ZRV0(NIU,NJU,NKU))
! initial values
ZXOR=0.5 * (XXHAT(2)+XXHAT(3)) 
ZYOR=0.5 * (XYHAT(2)+XYHAT(3))
ZDX= XXHAT(3)-XXHAT(2)
ZDY= XYHAT(3)-XYHAT(2)
ZZL=MZF(1,IKU,1,XZZ)
ZZL(:,:,NKU)=2*XZZ(:,:,NKU)-ZZL(:,:,NKU-1)
ZXMAX=ZXOR+(NIU-3)*ZDX
ZYMAX=ZYOR+(NJU-3)*ZDY
ZZMAX=ZZL(2,2,NKU-1)
!  conversion from km to meters
ZXOR=ZXOR*1.E-3
ZYOR=ZYOR*1.E-3
ZDX=ZDX*1.E-3
ZDY=ZDY*1.E-3
ZZL(:,:,:)=ZZL(:,:,:)*1.E-3
ZXMAX=ZXMAX*1.E-3
ZYMAX=ZYMAX*1.E-3
ZZMAX=ZZMAX*1.E-3
!
ZX00(:,:,:)=XSVT(:,:,:,NSV_LGBEG)*1.E-3   ! ZX0 in km
ZY00(:,:,:)=XSVT(:,:,:,NSV_LGBEG+1)*1.E-3 ! ZY0 in km
ZZ00(:,:,:)=XSVT(:,:,:,NSV_LGEND)*1.E-3   ! ZZ0 in km
!
IF (L2D) THEN
  WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
          ZZ00>ZZMAX)
    ZX00=ZSPVAL
    ZZ00=ZSPVAL
  END WHERE
ELSE
  WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
          ZY00<ZYOR .OR. ZY00>ZYMAX .OR. &
	      ZZ00>ZZMAX)
    ZX00=ZSPVAL
    ZY00=ZSPVAL
    ZZ00=ZSPVAL
  END WHERE
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.0    COMPUTE THE ORIGIN STEP BY STEP
!               -------------------------------
!
!
! General loop for the files where a reinitialisation of the tracers 
! is performed
DO JFILECUR=1,NFILES
  !
  CALL IO_FILE_ADD2LIST(TZTRACFILE,CFILES(NBRFILES(JFILECUR)),'UNKNOWN','READ',KLFITYPE=2,KLFIVERB=NVERB)
  CALL IO_FILE_OPEN_ll(TZTRACFILE)
!
!*       4.1  check if this file is a start instant
!
  GSTART=.FALSE.
  DO JLOOP=1,NFILES
    IF (NBRFILES(JFILECUR)==NSTART_SUPP(JLOOP) .OR. JFILECUR==NFILES) THEN
      INBR_START=INBR_START-1
      GSTART=.TRUE.
      EXIT
    END IF
  ENDDO
!
!*       4.2 read the potential temp or the water vapor at the start instant      
!
  IF (GSTART) THEN
    !
    CALL IO_READ_FIELD(TZTRACFILE,'DTCUR',TDTCUR_START)
    IHOUR   = INT(TDTCUR_START%TIME/3600.)
    ZREMAIN = MOD(TDTCUR_START%TIME,3600.)
    IMINUTE = INT(ZREMAIN/60.)
    ZSECOND = MOD(ZREMAIN,60.)
    WRITE(YDATE,FMT='(1X,I4.4,I2.2,I2.2,2X,I2.2,"H",I2.2,"M", &
         & F5.2,"S")') TDTCUR_START%TDATE, IHOUR,IMINUTE,ZSECOND  
    !
    CALL IO_READ_FIELD(TZTRACFILE,'THT',ZTH0(:,:,:))
    !
    CALL IO_READ_FIELD(TZTRACFILE,'RVT',ZRV0(:,:,:))
    ZRV0(:,:,:)=ZRV0(:,:,:)*1.E+3  ! ZRV0 in g/kg
    !
  END IF
!
!*       4.3  store the X0,Y0,Z0 field for the current start before 
!             computing the new origin
!
  IF (GSTART) THEN
    PRINT *,'INBR_START',INBR_START,' NBRFILES(JFILECUR)',NBRFILES(JFILECUR)
    WRITE(TZFIELD%CMNHNAME,'(A2,I2.2)')'X0',INBR_START
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'km'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)//YDATE
    PRINT *,'YCOMMENT = ',TRIM(TZFIELD%CCOMMENT)
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZX00(:,:,:))
    !
    WRITE(TZFIELD%CMNHNAME,'(A2,I2.2)')'Y0',INBR_START
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'km'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)//YDATE
    PRINT *,'YCOMMENT = ',TRIM(TZFIELD%CCOMMENT)
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZY00(:,:,:))
    !
    WRITE(TZFIELD%CMNHNAME,'(A2,I2.2)')'Z0',INBR_START
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'km'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)//YDATE
    PRINT *,'YCOMMENT = ',TRIM(TZFIELD%CCOMMENT)
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZZ00(:,:,:))
  END IF
!
!
!*       4.6   compute and store potential temp and water vapor at the origin
!
  IF (GSTART) THEN
    !
    CALL INTERPXYZ(ZX00,ZY00,ZZ00,     &
                   ZTH0,ZWORK1         )
    !
    CALL INTERPXYZ(ZX00,ZY00,ZZ00,     &
                   ZRV0,ZWORK2         )
    !
  END IF
!
  IF (GSTART) THEN
    !
    WRITE(TZFIELD%CMNHNAME,'(A3,I2.2)')'TH0',INBR_START
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)//YDATE
    PRINT *,'YCOMMENT = ',TRIM(TZFIELD%CCOMMENT)
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK1(:,:,:))
    !
    WRITE(TZFIELD%CMNHNAME,'(A3,I2.2)')'RV0',INBR_START
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'g kg-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(TZFIELD%CMNHNAME)//YDATE
    PRINT *,'YCOMMENT = ',TRIM(TZFIELD%CCOMMENT)
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZWORK2(:,:,:))
  ENDIF
!*       4.4   compute the origin of the particules using one more segment
!
  IF (JFILECUR /= NFILES) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = '' !Unknown comment
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    TZFIELD%CMNHNAME   = 'LGXT'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    CALL IO_READ_FIELD(TZTRACFILE,TZFIELD,ZX0)
    ZX0(:,:,:)=ZX0(:,:,:)*1.E-3   ! ZX0 in km
    !
    TZFIELD%CMNHNAME   = 'LGYT'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    CALL IO_READ_FIELD(TZTRACFILE,TZFIELD,ZY0)
    ZY0(:,:,:)=ZY0(:,:,:)*1.E-3   ! ZY0 in km
    !
    TZFIELD%CMNHNAME   = 'LGZT'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    CALL IO_READ_FIELD(TZTRACFILE,TZFIELD,ZZ0)
    ZZ0(:,:,:)=ZZ0(:,:,:)*1.E-3   ! ZZ0 in km
    !
    ! old position of the set of particles
    ZWORK1=ZX00
    ZWORK2=ZY00
    ZWORK3=ZZ00
    !
    IF (L2D) THEN
      CALL INTERPXYZ(ZWORK1,ZWORK2,ZWORK3,         &
                     ZX0,ZX00,ZZ0,ZZ00             )
    ELSE
      CALL INTERPXYZ(ZWORK1,ZWORK2,ZWORK3,         &
                     ZX0,ZX00,ZY0,ZY00,ZZ0,ZZ00    )
    END IF
    !
    IF (L2D) THEN
      WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
              ZZ00>ZZMAX)
        ZX00=ZSPVAL
        ZZ00=ZSPVAL
      END WHERE
    ELSE
      WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
              ZY00<ZYOR .OR. ZY00>ZYMAX .OR. &
              ZZ00>ZZMAX)
        ZX00=ZSPVAL
        ZY00=ZSPVAL
        ZZ00=ZSPVAL
      END WHERE
    END IF
    !
  END IF
!
!*       4.5   close the input file
!
  CALL IO_FILE_CLOSE_ll(TZTRACFILE)
!
END DO
!
PRINT*, ' '
PRINT*, 'DIAG AFTER ORIGIN COMPUTATIONS AND STORAGE'
!
!-------------------------------------------------------------------------------
!
!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
SUBROUTINE INTERPXYZ(PX,PY,PZ,PIN1,POUT1,PIN2,POUT2,PIN3,POUT3)
!
!
!*      0. DECLARATIONS
!          ------------
!
!*       0.1  declaration of arguments
!
REAL, INTENT(IN),  DIMENSION(:,:,:)           :: PX,PY,PZ
REAL, INTENT(IN),  DIMENSION(:,:,:)           :: PIN1
REAL, INTENT(OUT), DIMENSION(:,:,:)           :: POUT1
REAL, INTENT(IN),  DIMENSION(:,:,:), OPTIONAL :: PIN2,PIN3
REAL, INTENT(OUT), DIMENSION(:,:,:), OPTIONAL :: POUT2,POUT3   
!
!*       0.2  declaration of local variables
!
INTEGER  :: JI,JJ,JK,JKK    ! loop index
INTEGER  :: II,IJ,IK        ! grid index for the interpolation
REAL     :: ZXREL,ZYREL     ! fractional grid index for the interpolation
REAL, DIMENSION(SIZE(PIN1,3)) :: ZZLXY ! vertical grid at the interpolated point
REAL     :: ZEPS1,ZEPS2,ZEPS3          ! coeff. for the interpolation
REAL     :: ZX,ZY,ZZ
LOGICAL  :: GEXT
!
!-------------------------------------------------------------------------------
!
DO JK=1,NKU
  DO JJ=1,NJU
    DO JI=1,NIU
      !
      ZX=PX(JI,JJ,JK) 
      ZY=PY(JI,JJ,JK)
      ZZ=PZ(JI,JJ,JK)
      !
      ! remove external points
      IF (L2D) THEN
        GEXT=(ZX==ZSPVAL).OR.(ZZ==ZSPVAL)
      ELSE
        GEXT=(ZX==ZSPVAL).OR.(ZY==ZSPVAL).OR.(ZZ==ZSPVAL)
      END IF
      IF (GEXT) THEN
        POUT1(JI,JJ,JK) = ZSPVAL
        IF (PRESENT(PIN2)) THEN
          POUT2(JI,JJ,JK) = ZSPVAL
        END IF
        IF (PRESENT(PIN3)) THEN
          POUT3(JI,JJ,JK) = ZSPVAL
        ENDIF
        !
        CYCLE
        !
      END IF
      !
      ZXREL=(ZX-ZXOR)/ZDX+2
      ZYREL=(ZY-ZYOR)/ZDY+2
      !
      II=FLOOR(ZXREL)
      IJ=FLOOR(ZYREL)
      !
      ZEPS1=ZXREL-REAL(II)
      ZEPS2=ZYREL-REAL(IJ)
      IF (L2D) ZEPS2=0.
      !
      DO JKK=1,NKU
        ZZLXY(JKK)=ZEPS2*(ZEPS1*(ZZL(II+1,IJ+1,JKK))+(1-ZEPS1)*(ZZL(II,IJ+1,JKK)))     &
             + (1-ZEPS2)*(ZEPS1*(ZZL(II+1,IJ,JKK))+(1-ZEPS1)*(ZZL(II,IJ,JKK)))
      ENDDO
      !
      IK=999
      DO JKK=2,NKU
        IF (ZZLXY(JKK).GE.ZZ) THEN
          IK=JKK-1
          EXIT 
        ENDIF
      ENDDO
      !
      IF (IK==999) THEN
        PRINT*,'PROBLEM AT POINT',II,IJ
        PRINT*,'XREL, YREL, Z =',ZXREL,ZYREL,ZZ
        PRINT*,'ZZLXY(NKU)',ZZLXY(NKU)
!callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','COMPUTE_R00','')
      END IF 
      !
      ZEPS3=(ZZ-ZZLXY(IK))/(ZZLXY(IK+1)-ZZLXY(IK))
      !
      POUT1(JI,JJ,JK) =                                                       & 
        ZEPS3 *                                                               &
      (  ZEPS2*(ZEPS1*(PIN1(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PIN1(II,IJ+1,IK+1)))  &
       + (1-ZEPS2)*(ZEPS1*(PIN1(II+1,IJ,IK+1))+(1-ZEPS1)*(PIN1(II,IJ,IK+1)))  &
      )                                                                       & 
      + (1-ZEPS3) *                                                           &
      (  ZEPS2*(ZEPS1*(PIN1(II+1,IJ+1,IK))+(1-ZEPS1)*(PIN1(II,IJ+1,IK)))      &
       + (1-ZEPS2)*(ZEPS1*(PIN1(II+1,IJ,IK))+(1-ZEPS1)*(PIN1(II,IJ,IK)))      &
      )
      IF (PRESENT(POUT2)) THEN
        POUT2(JI,JJ,JK) =                                                     & 
          ZEPS3 *                                                             &
        (  ZEPS2*(ZEPS1*(PIN2(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PIN2(II,IJ+1,IK+1)))&
         + (1-ZEPS2)*(ZEPS1*(PIN2(II+1,IJ,IK+1))+(1-ZEPS1)*(PIN2(II,IJ,IK+1)))&
        )                                                                     & 
        + (1-ZEPS3) *                                                         &
        (  ZEPS2*(ZEPS1*(PIN2(II+1,IJ+1,IK))+(1-ZEPS1)*(PIN2(II,IJ+1,IK)))    &
         + (1-ZEPS2)*(ZEPS1*(PIN2(II+1,IJ,IK))+(1-ZEPS1)*(PIN2(II,IJ,IK)))    &
        )
      ENDIF
        !
      IF (PRESENT(POUT3)) THEN
        POUT3(JI,JJ,JK) =                                                     & 
          ZEPS3 *                                                             &
        (  ZEPS2*(ZEPS1*(PIN3(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PIN3(II,IJ+1,IK+1)))&
         + (1-ZEPS2)*(ZEPS1*(PIN3(II+1,IJ,IK+1))+(1-ZEPS1)*(PIN3(II,IJ,IK+1)))&
        )                                                                     &
        + (1-ZEPS3) *                                                         &
        (  ZEPS2*(ZEPS1*(PIN3(II+1,IJ+1,IK))+(1-ZEPS1)*(PIN3(II,IJ+1,IK)))    &
         + (1-ZEPS2)*(ZEPS1*(PIN3(II+1,IJ,IK))+(1-ZEPS1)*(PIN3(II,IJ,IK)))    &
        )
      ENDIF
      !
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INTERPXYZ
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE COMPUTE_R00 
