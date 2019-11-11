!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
      MODULE MODI_SET_RELFRC
!     ###################
!
INTERFACE
!
SUBROUTINE SET_RELFRC(TPEXPREFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA), INTENT(IN)  :: TPEXPREFILE ! input data file
!
END SUBROUTINE SET_RELFRC
!
END INTERFACE
!
END MODULE MODI_SET_RELFRC
!
!     ##########################
      SUBROUTINE SET_RELFRC(TPEXPREFILE)
!     ##########################
!
!!*** *SET_RELFRC * -  to initialize Advective or relaxation forcing fields for different months
!!    
!!
!!
!!    Author : P. Peyrille
!!    PURPOSE
!!    -------
!!      Same purpose as set_frc.f90 but for 2D forcing. 
!!      The routine reads the external ascii files and performs vertical
!!     interpolations to get the forcing fields on the model levels.
!! 
!!     Caution : File must be written under the form (lat, lev, th_frc, rv_frc) 
!!                See read_ascp, read_lat_press for more details
!!
!!**  METHOD
!!    ------
!!      Pressure or height level data are considered.  The forcing data must be
!       given on the same horizontal grid as MNH domain but can have different
!!      vertical resolution. 
!!
!!     Depending on YREL value (P or Z REL), an interpolation is
!!     made by assuming the ground is at z=0 and the first pressure level is at the
!!     surface.
!!
!!     Once ZREL_FRC is found in prep_ideal_case, the namlist is read as follows:
!!
!!
!!     NRELFRC            ! Number of time-dependent forcing files to be read 
!      YREL               ! type of relaxation file (PREL, or ZREL)
!      NPRESSLEV_REL      ! nb of levels  in relxaxation file 
!      DATE NB 1
!     File 1 for mean atm profile for vertical interpolation of relax forcing from  P=> to Z levels
!     File 1 for relaxation forcing
!      DATE NB 2
!     File 2 for mean relaxation profile atm
!     File 2 for relaxation forcing
!
!
!     EXTERNALS:
!!      MODI_READ_ASCP         : Reads a 1D field in an ascii file
!!      MODI_READ_ASC_LATPRESS : Reads a 2D field (lat,lev) ------
!! 
!!    MODIFICATIONS
!!    -------------
!!      03/02/10 (Tomasini) USE MODD_RELFRC_n for grid-nesting
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      P.Wautelet  28/03/2018 : use overloaded comparison operator for date_time
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST 
USE MODD_FRC
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IO_ll, ONLY : TFILEDATA
USE MODD_LUNIT_n
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_REF
USE MODD_RELFRC_n
! 
USE MODE_DATETIME
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
USE MODE_THERMO 
!
USE MODD_DIM_n 
USE MODI_HEIGHT_PRESS
USE MODI_PRESS_HEIGHT
USE MODI_READ_ASC_LATPRESS
USE MODI_READ_ASCP
USE MODI_THETAVPU_THETAVPM 
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments :
!
TYPE(TFILEDATA), INTENT(IN)  :: TPEXPREFILE ! input data file
!
!*       0.2   Declarations of local variables :
!
INTEGER :: ILUPRE,IRESP ! logical unit number of the  EXPRE and FM return code
INTEGER :: ILUOUT       ! Logical unit number for output-listing   
INTEGER :: JKT,JL,JK,JI! Loop control   
INTEGER :: IIU,IIB,IIE !dimensions du modele
INTEGER :: IJU,IJB,IJE !dimensions du modele
INTEGER :: IKU,IKB,IKE !dimensions du modele
!
REAL    :: ZRVSRD    ! Rv/Rd
REAL    :: ZDZ1SDH, ZDZ2SDH, ZDZSDH
INTEGER :: JKLEV
!
REAL    :: ZZGROUNDF    ! height and Pressure at ground
CHARACTER(LEN=48)   :: CFNAM_MEANVAR_REL
CHARACTER(LEN=28) :: CFNAM_REL
!
REAL, DIMENSION(:), ALLOCATABLE::     ZHEIGHTMFR,ZHEIGHTFR,ZTHVUFR,ZTHVUF
REAL, DIMENSION(:), ALLOCATABLE::     ZZHATM
REAL, DIMENSION(:), ALLOCATABLE::     ZPRESS_REL
REAL, DIMENSION(:), ALLOCATABLE::     ZTHDFR,ZTHVFR,ZRVFR
!
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZRVREL1D,ZTHREL1D,ZVREL1D
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZTHREL,ZRVREL,ZVREL
!
REAL, DIMENSION(:), ALLOCATABLE::     ZLAT_FRC
CHARACTER(LEN=6)                :: YREL         ! choice of zfrc or pfrc
!-------------------------------------------------------------------------------
!
print*,"!*	 1.     PROLOGUE : RETRIEVE LOGICAL UNIT NUMBERS" 
!	        ----------------------------------------
!                           
ILUPRE = TPEXPREFILE%NLU
ILUOUT = TLUOUT%NLU
!
!-------------------------------------------------------------------------------
!
print*,"!*	 2.     COMPUTE FORCING FIELDS PROFILES"
!	        -------------------------------
ZRVSRD  = XRV/XRD
!
!
print*,"!	     2.1	Compute array size and allocate memory"
!
READ(ILUPRE,*) NRELFRC            ! Number of time-dependent forcing soundings
READ(ILUPRE,*) YREL               ! type of relaxation file (PREL, or ZREL)
READ(ILUPRE,*) NPRESSLEV_REL      ! nb of levels for low leves forcing=nb lev in 2D model
 
! CAUTION: number of forcing times is limited by the WRITE format 99(8E10.3)
!          and also by the name of forcing variables (format I3.3)
!          You have to modify those if you need more forcing times :-(
!
IF (NRELFRC > 99*8) THEN
  CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_RELFRC','maximum forcing times NRELFRC is 99*8')
END IF
!

IIU=SIZE(XXHAT)
IJU=SIZE(XYHAT)
IKU=SIZE(XZHAT)

CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
! allocation !
!
print*,"! temporary forcing alloation"
!
! For relaxation forcing
ALLOCATE(ZTHREL(IIU,NPRESSLEV_REL,NRELFRC))
ALLOCATE(ZRVREL(IIU,NPRESSLEV_REL,NRELFRC))
ALLOCATE(ZVREL(IIU,NPRESSLEV_REL,NRELFRC))
ALLOCATE(ZLAT_FRC(IIE))
ALLOCATE(ZRVFR(NPRESSLEV_REL))
!
ALLOCATE(ZTHDFR(NPRESSLEV_REL))
ALLOCATE(ZTHVUFR(NPRESSLEV_REL))
ALLOCATE(ZTHVFR(NPRESSLEV_REL))
ALLOCATE(ZPRESS_REL(NPRESSLEV_REL))
ALLOCATE(ZHEIGHTMFR(NPRESSLEV_REL))
ALLOCATE(ZHEIGHTFR(NPRESSLEV_REL))

!
! allocations pour le module moddb_advfrc
!  Adv forcing
! 
! For reading in PRE_IDEA1.nam
ALLOCATE(ZZHATM(IKU))
!
! relaxation profile
ALLOCATE(ZRVREL1D(IIU,IKU,NRELFRC))
ALLOCATE(ZVREL1D(IIU,IKU,NRELFRC))
ALLOCATE(ZTHREL1D(IIU,IKU,NRELFRC))
!
ALLOCATE(XTHREL(IIU,IJU,IKU,NRELFRC))
ALLOCATE(XRVREL(IIU,IJU,IKU,NRELFRC))
ALLOCATE(TDTRELFRC(NRELFRC))
!
! initilisation
XTHREL(:,:,:,:)=0.
XRVREL(:,:,:,:) = 0.
!
!
 !
!  3.   READ ASCII FILES FOR RELAXTAION FORCING 
!	        -------------------------------

DO JKT = 1,NRELFRC


!   Reading the date and the filename  of mean state and frc
    READ(ILUPRE,*) TDTRELFRC(JKT)%TDATE%YEAR,  &
                   TDTRELFRC(JKT)%TDATE%MONTH, &
                   TDTRELFRC(JKT)%TDATE%DAY,   &
                   TDTRELFRC(JKT)%TIME
!                   ; Read filenames
    READ(ILUPRE,*)  CFNAM_MEANVAR_REL  
    READ(ILUPRE,*)  CFNAM_REL
    !
!

! 3.3 READ AND INTERPOLATE CLIMATOLOGICAL THETA AND RV PROFILE (2D)
  ! ----------------------------------------------------------------
  ! Call read_asc_latpress for relaxation file 
  ! Given on pressure levels




  IF (YREL=='PREL2D') THEN
   CALL READ_ASCP(CFNAM_MEANVAR_REL,NPRESSLEV_REL,ZTHDFR,ZRVFR)
   CALL READ_ASC_LATPRESS(CFNAM_REL,NPRESSLEV_REL,ZLAT_FRC,ZPRESS_REL,ZTHREL(:,:,JKT), & 
                         ZRVREL(:,:,JKT))
 !
   ZTHVFR(:) = ZTHDFR(:) * (1.+ZRVSRD*ZRVFR(:))/(1.+ZRVFR(:))
 !

   ! Compute the heights at the mass levels of the RS
  !
   ZZGROUNDF=0.

    ZHEIGHTMFR(:) = HEIGHT_PRESS(ZPRESS_REL,ZTHVFR,ZPRESS_REL(1),ZTHVFR(1),ZZGROUNDF)
   ! Compute thetav and heights of the wind levels
    ZTHVUFR(:) = THETAVPU_THETAVPM(ZPRESS_REL,ZPRESS_REL,ZTHVFR)
    ZHEIGHTFR(:) = HEIGHT_PRESS(ZPRESS_REL,ZTHVUFR,ZPRESS_REL(1),ZTHVFR(1),ZZGROUNDF)
   !
  !
  ELSE IF (YREL=='ZREL2D') THEN
   CALL READ_ASC_LATPRESS(CFNAM_REL,NPRESSLEV_REL,ZLAT_FRC,ZHEIGHTFR(:),ZTHREL(:,:,JKT), & 
                         ZRVREL(:,:,JKT))

  END IF
  !
  ZZHATM(1:IKU-1) = 0.5*(XZHAT(2:IKU)+XZHAT(1:IKU-1))
  ZZHATM(IKU) = 2.*XZHAT(IKU)-ZZHATM(IKU-1)
!

  DO JK = 1,IKU
    IF (ZZHATM(JK) <= ZHEIGHTFR(1)) THEN
!
 ! extrapolation below the first level
      ZDZSDH  = (ZZHATM(JK)-ZHEIGHTFR(1)) / (ZHEIGHTFR(2)-ZHEIGHTFR(1))
      !
      ZRVREL1D(IIB:IIE,JK,JKT) = ZRVREL(IIB:IIE,1,JKT) + & 
      (ZRVREL(IIB:IIE,2,JKT) - ZRVREL(IIB:IIE,1,JKT)) * ZDZSDH
      ZVREL1D(IIB:IIE,JK,JKT) = ZVREL(IIB:IIE,1,JKT) + & 
      (ZVREL(IIB:IIE,2,JKT) - ZVREL(IIB:IIE,1,JKT)) * ZDZSDH
      ZTHREL1D(IIB:IIE,JK,JKT) = ZtHREL(IIB:IIE,1,JKT) + & 
      (ZTHREL(IIB:IIE,2,JKT) - ZTHREL(IIB:IIE,1,JKT)) * ZDZSDH

    ELSE IF (ZZHATM(JK) > ZHEIGHTFR(NPRESSLEV_REL) ) THEN
!
 ! extrapolation above the last level
      ZDZSDH  = (ZZHATM(JK) - ZHEIGHTFR(NPRESSLEV_REL)) /      &
                (ZHEIGHTFR(NPRESSLEV_REL) - ZHEIGHTFR(NPRESSLEV_REL-1))
      !
      ZRVREL1D(IIB:IIE,JK,JKT)  = ZRVREL(IIB:IIE,NPRESSLEV_REL,JKT) +  & 
      (ZRVREL(IIB:IIE,NPRESSLEV_REL,JKT)-ZRVREL(IIB:IIE,NPRESSLEV_REL-1,JKT)) * ZDZSDH
      ZVREL1D(IIB:IIE,JK,JKT)  = ZVREL(IIB:IIE,NPRESSLEV_REL,JKT) +  & 
      (ZVREL(IIB:IIE,NPRESSLEV_REL,JKT)-ZVREL(IIB:IIE,NPRESSLEV_REL-1,JKT)) * ZDZSDH
      ZTHREL1D(IIB:IIE,JK,JKT)  = ZTHREL(IIB:IIE,NPRESSLEV_REL,JKT) + & 
      (ZTHREL(IIB:IIE,NPRESSLEV_REL,JKT)-ZTHREL(IIB:IIE,NPRESSLEV_REL-1,JKT)) * ZDZSDH

    ELSE
!
 ! interpolation between first and last levels
!
      DO JKLEV = 1,NPRESSLEV_REL-1
        IF ( (ZZHATM(JK) > ZHEIGHTFR(JKLEV)).AND. &
             (ZZHATM(JK) <= ZHEIGHTFR(JKLEV+1))   ) THEN
          ZDZ1SDH = (ZZHATM(JK) - ZHEIGHTFR(JKLEV)) /  &
                    (ZHEIGHTFR(JKLEV+1)-ZHEIGHTFR(JKLEV))
          ZDZ2SDH = 1.- ZDZ1SDH
          ZRVREL1D(IIB:IIE,JK,JKT)  = ZRVREL(IIB:IIE,JKLEV,JKT)*ZDZ2SDH & 
          + ZRVREL(IIB:IIE,JKLEV+1,JKT)*ZDZ1SDH
!
          ZTHREL1D(IIB:IIE,JK,JKT)  = ZTHREL(IIB:IIE,JKLEV,JKT)*ZDZ2SDH & 
          + ZTHREL(IIB:IIE,JKLEV+1,JKT)*ZDZ1SDH
          !

        END IF
      END DO
    END IF
  END DO 

!  . 2.2  END CASES , AND TREAT BOUNDARY POINTS
! -----------------------------------------------------------------
!
! Arrays that will be stored in the moddb_advfrc module
!
  ! Expand arrays to 3D
  !
  XTHREL=SPREAD(ZTHREL1D,2,IJU)
  XRVREL=SPREAD(ZRVREL1D,2,IJU)
  ! Fill the first and last point of arrays
  !
  XTHREL(1,:,:,:) = XTHREL(IIB,:,:,:)
  XTHREL(IIU,:,:,:) = XTHREL(IIE,:,:,:)
  XRVREL(1,:,:,:) =XRVREL (IIB,:,:,:)
  XRVREL(IIU,:,:,:) = XRVREL(IIE,:,:,:)



END DO ! End of loop in time


!-----------------------------------------------------------------------------
!
!*        3.    PRINT FORCING FIELDS
!         --------------------------
!

WRITE(UNIT=ILUOUT,FMT='(" THERE ARE ",I2," REL FORCING AT:")') NRELFRC
DO JL = 1 , NRELFRC
  WRITE(UNIT=ILUOUT,FMT='(F9.0, "s, date:", I3, "/", I3, "/", I5)') &
    TDTRELFRC(JL)%TIME, &
    TDTRELFRC(JL)%TDATE%DAY,   &
    TDTRELFRC(JL)%TDATE%MONTH, &
    TDTRELFRC(JL)%TDATE%YEAR
END DO
!
DO JKT = 2,NRELFRC-1
  IF ( TDTRELFRC(JKT-1) >= TDTRELFRC(JKT) ) THEN
    WRITE(ILUOUT,*) &
      "SET_FRC ERROR: sounding", JKT-1, " is given for a later time than", JKT
    WRITE(ILUOUT,*) &
      "               soundings have to be entered in increasing temporal order"
    WRITE(ILUOUT,*) "SOUNDING TIME ", JKT-1, " IS: "
    WRITE(ILUOUT,*) TDTRELFRC(JKT-1)%TDATE%YEAR,  &
                    TDTRELFRC(JKT-1)%TDATE%MONTH, &
                    TDTRELFRC(JKT-1)%TDATE%DAY,   &
                    TDTRELFRC(JKT-1)%TIME
    WRITE(ILUOUT,*) "SOUNDING TIME ", JKT, " IS: "
    WRITE(ILUOUT,*) TDTRELFRC(JKT)%TDATE%YEAR,  &
                    TDTRELFRC(JKT)%TDATE%MONTH, &
                    TDTRELFRC(JKT)%TDATE%DAY,   &
                    TDTRELFRC(JKT)%TIME
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_RELFRC','')
  END IF
END DO

!DEALLOCATE ARRAYS
! 
DEALLOCATE(ZTHVFR)
DEALLOCATE(ZTHDFR)
DEALLOCATE(ZTHVUFR)
DEALLOCATE(ZRVFR)

DEALLOCATE(ZHEIGHTFR)
! pour lecture dans PREIDEA
DEALLOCATE(ZZHATM)
DEALLOCATE(ZPRESS_REL)

DEALLOCATE(ZRVREL1D)
DEALLOCATE(ZVREL1D)
DEALLOCATE(ZTHREL1D)

!-------------------------------------------------------------------------------
END SUBROUTINE SET_RELFRC
