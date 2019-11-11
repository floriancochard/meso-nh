!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
      MODULE MODI_SETADVFRC
!     ###################
!
INTERFACE
!
SUBROUTINE SET_ADVFRC(TPEXPREFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
TYPE(TFILEDATA), INTENT(IN)  :: TPEXPREFILE ! input data file
!
END SUBROUTINE SET_ADVFRC
!
END INTERFACE
!
END MODULE MODI_SETADVFRC
!
!     ##################################
      SUBROUTINE SET_ADVFRC(TPEXPREFILE)
!     ##################################
!
!!*** *SET_ADVFRC * -  to initialize Advective forcing fields for different months
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
!!     Depending on YADV  value (P or Z ADV), an interpolation is
!!     made by assuming the ground is at z=0 and the first pressure level is at the
!!     surface.
!!
!!     Once ZADV_FRC is found in prep_ideal_case, the namlist is read as follows:
!!
!!
!!     NADVFRC            ! Number of time-dependent forcing files to be read 
!      YADV               ! type of advection file (PADV for pressure levels, or ZADV for height)
!      NPRESSLEV_ADV      ! nb of levels in the  advective forcing file
!      DATE NB 1
!     File 1 for mean atm profile for vertical interpolation of adv forcing from  P=> to Z levels
!     File 1 for Adv forcing
!      DATE NB 2
!     File 2 for mean adv profile atm
!     File 2 for Adv forcing
!
!
!     EXTERNALS:
!!      MODI_READ_ASCP         : Reads a 1D field in an ascii file
!!      MODI_READ_ASC_LATPRESS : Reads a 2D field (lat,lev) ------
!! 
!!    MODIFICATIONS
!!    -------------
!!      03/02/10 (Tomasini) USE MODDB_ADVFRC_n for grid-nesting
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      P.Wautelet  28/03/2018 : use overloaded comparison operator for date_time
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ADVFRC_n
USE MODD_CONF
USE MODD_CST 
USE MODD_DIM_n 
USE MODD_FRC
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IO_ll,      ONLY: TFILEDATA
USE MODD_LUNIT_n,    ONLY: TLUOUT
USE MODD_PARAMETERS, ONLY: JPHEXT, JPVEXT
USE MODD_REF
! 
USE MODE_DATETIME
USE MODE_IO_ll
USE MODE_MSG
USE MODE_THERMO
!
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
CHARACTER(LEN=48)   :: CFNAM_MEANVAR_ADV
CHARACTER(LEN=48) :: CFNAM_ADV
!
REAL, DIMENSION(:), ALLOCATABLE::     ZHEIGHTMF,ZHEIGHTF,ZTHVUF
REAL, DIMENSION(:), ALLOCATABLE::     ZZHATM
REAL, DIMENSION(:), ALLOCATABLE::     ZTHDF,ZRVF,ZPRESS_ADV,ZTHVF
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZTHFRC,ZRVFRC
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZDRVFRC1D,ZDTHFRC1D,ZDVFRC1D
!
!
REAL, DIMENSION(:), ALLOCATABLE::     ZLAT_FRC
CHARACTER(LEN=6)                :: YADV         ! choice of zfrc or pfrc
!-------------------------------------------------------------------------------
!
print*,"!*	 1.     PROLOGUE : RETRIEVE LOGICAL UNIT NUMBERS "
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
!	     2.1	Compute array size and allocate memory
!
READ(ILUPRE,*) NADVFRC            ! Number of time-dependent forcing soundings
READ(ILUPRE,*) YADV               ! type of advection file (PADV, or ZADV)
READ(ILUPRE,*) NPRESSLEV_ADV      ! nb of levels for low leves forcing=nb lev in 2D model
 
! CAUTION: number of forcing times is limited by the WRITE format 99(8E10.3)
!          and also by the name of forcing variables (format I3.3)
!          You have to modify those if you need more forcing times :-(
!
IF (NADVFRC > 99*8) CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_ADVFRC','maximum forcing times NADVFRC is 99*8')
!
IIU=SIZE(XXHAT)
IJU=SIZE(XYHAT)
IKU=SIZE(XZHAT)
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB= 1+ JPVEXT 
IKE= IKU-JPVEXT 
!
! allocation !
!
! temporary forcing alloation
ALLOCATE(ZPRESS_ADV(NPRESSLEV_ADV))
ALLOCATE(ZHEIGHTF(NPRESSLEV_ADV))
ALLOCATE(ZTHDF(NPRESSLEV_ADV))
ALLOCATE(ZRVF(NPRESSLEV_ADV))
ALLOCATE(ZTHVF(NPRESSLEV_ADV))
!
ALLOCATE(ZTHFRC(IIU,NPRESSLEV_ADV,NADVFRC))
ALLOCATE(ZRVFRC(IIU,NPRESSLEV_ADV,NADVFRC))
!
ALLOCATE(ZHEIGHTMF(NPRESSLEV_ADV))
ALLOCATE(ZTHVUF(NPRESSLEV_ADV))
! allocations pour le module moddb_advfrc
!  Adv forcing
ALLOCATE(TDTADVFRC(NADVFRC))
ALLOCATE(XDTHFRC(IIU,IJU,IKU,NADVFRC))
ALLOCATE(XDRVFRC(IIU,IJU,IKU,NADVFRC))
! 
! For reading in PRE_IDEA1.nam
ALLOCATE(ZZHATM(IKU))
ALLOCATE(ZDRVFRC1D(IIU,IKU,NADVFRC))
ALLOCATE(ZDTHFRC1D(IIU,IKU,NADVFRC))
ALLOCATE(ZDVFRC1D(IIU,IKU,NADVFRC))
!
ALLOCATE(ZLAT_FRC(IIE))
!
! initilisation
XDRVFRC(:,:,:,:) = 0.
XDTHFRC(:,:,:,:) = 0.
!
!
print*,"!  3.   READ ASCII FILES FOR ADVECTIVE FORCING "
!	        -------------------------------

DO JKT = 1,NADVFRC
!
!   Reading the date and the filename  of mean state and frc
    READ(ILUPRE,*) TDTADVFRC(JKT)%TDATE%YEAR,  &
                   TDTADVFRC(JKT)%TDATE%MONTH, &
                   TDTADVFRC(JKT)%TDATE%DAY,   &
                   TDTADVFRC(JKT)%TIME
!                   ; Read filenames
    READ(ILUPRE,*)  CFNAM_MEANVAR_ADV  
    READ(ILUPRE,*)  CFNAM_ADV
    !
!
print*,"  ! 3.2  READ AND INTERPOLATE FORCING"
  !   ---------------------------------------------------------------------------------
  !! 3.2.1  Case ADV forcing on height levels
  !! ZHEIGHTF is read in the external file
  print*,"YADV=",YADV
  IF (YADV=='ZADV2D') THEN
          print*,"call READ_ASC_LATPRESS"
    CALL READ_ASC_LATPRESS(CFNAM_ADV,NPRESSLEV_ADV,ZLAT_FRC,ZHEIGHTF(:),ZTHFRC(:,:,JKT), & 
                         ZRVFRC(:,:,JKT))
                         print*,"apres READ_ASC_LATPRESS"
  ELSE IF (YADV=='PADV2D') THEN
    CALL READ_ASCP(CFNAM_MEANVAR_ADV,NPRESSLEV_ADV,ZTHDF,ZRVF)
  !! 3.2.1  Case ADV forcing on pressure levels
    CALL READ_ASC_LATPRESS(CFNAM_ADV,NPRESSLEV_ADV,ZLAT_FRC,ZPRESS_ADV,ZTHFRC(:,:,JKT), & 
                         ZRVFRC(:,:,JKT))
 !
   ZTHVF(:) = ZTHDF(:) * (1.+ZRVSRD*ZRVF(:))/(1.+ZRVF(:))
 !

   ! Compute the heights at the mass levels of the RS assume surface at z=0m
   ZZGROUNDF=0.
   ZHEIGHTMF(:) = HEIGHT_PRESS(ZPRESS_ADV,ZTHVF,ZPRESS_ADV(1),ZTHVF(1),ZZGROUNDF)
   ! Compute thetav and heights of the wind levels
    ZTHVUF(:) = THETAVPU_THETAVPM(ZPRESS_ADV,ZPRESS_ADV,ZTHVF)
    ZHEIGHTF(:) = HEIGHT_PRESS(ZPRESS_ADV,ZTHVUF,ZPRESS_ADV(1),ZTHVF(1),ZZGROUNDF)
   !
  END IF
 !
  ZZHATM(1:IKU-1) = 0.5*(XZHAT(2:IKU)+XZHAT(1:IKU-1))
  ZZHATM(IKU) = 2.*XZHAT(IKU)-ZZHATM(IKU-1)

print*,"  !! 3.2.2  Vertical interpolation"
  DO JK = 1,IKU
    IF (ZZHATM(JK) <= ZHEIGHTF(1)) THEN
!
print*,"! extrapolation below the first level"
!   
   
      ZDZSDH  = (ZZHATM(JK)-ZHEIGHTF(1)) / (ZHEIGHTF(2)-ZHEIGHTF(1))
      ZDRVFRC1D(IIB:IIE,JK,JKT) = ZRVFRC(IIB:IIE,1,JKT) + & 
      (ZRVFRC(IIB:IIE,2,JKT) - ZRVFRC(IIB:IIE,1,JKT)) * ZDZSDH
      ZDTHFRC1D(IIB:IIE,JK,JKT) = ZTHFRC(IIB:IIE,1,JKT) + & 
      (ZTHFRC(IIB:IIE,2,JKT) - ZTHFRC(IIB:IIE,1,JKT)) * ZDZSDH
    ELSE IF (ZZHATM(JK) > ZHEIGHTF(NPRESSLEV_ADV) ) THEN
!
print*,"! extrapolation above the last level"
!
      ZDZSDH  = (ZZHATM(JK) - ZHEIGHTF(NPRESSLEV_ADV)) /      &
                (ZHEIGHTF(NPRESSLEV_ADV) - ZHEIGHTF(NPRESSLEV_ADV-1))
      ZDRVFRC1D(IIB:IIE,JK,JKT)  = ZRVFRC(IIB:IIE,NPRESSLEV_ADV,JKT) +  & 
      (ZRVFRC(IIB:IIE,NPRESSLEV_ADV,JKT)-ZRVFRC(IIB:IIE,NPRESSLEV_ADV-1,JKT)) * ZDZSDH
      ZDTHFRC1D(IIB:IIE,JK,JKT)  = ZTHFRC(IIB:IIE,NPRESSLEV_ADV,JKT) + & 
      (ZTHFRC(IIB:IIE,NPRESSLEV_ADV,JKT)-ZTHFRC(IIB:IIE,NPRESSLEV_ADV-1,JKT)) * ZDZSDH
    ELSE
!
print*,"! interpolation between first and last levels"
!
      DO JKLEV = 1,NPRESSLEV_ADV-1
        IF ( (ZZHATM(JK) > ZHEIGHTF(JKLEV)).AND. &
             (ZZHATM(JK) <= ZHEIGHTF(JKLEV+1))   ) THEN
          ZDZ1SDH = (ZZHATM(JK) - ZHEIGHTF(JKLEV)) /  &
                    (ZHEIGHTF(JKLEV+1)-ZHEIGHTF(JKLEV))
          ZDZ2SDH = 1.- ZDZ1SDH
          ZDRVFRC1D(IIB:IIE,JK,JKT)  = ZRVFRC(IIB:IIE,JKLEV,JKT)*ZDZ2SDH & 
          + ZRVFRC(IIB:IIE,JKLEV+1,JKT)*ZDZ1SDH
          ZDTHFRC1D(IIB:IIE,JK,JKT)  = ZTHFRC(IIB:IIE,JKLEV,JKT)*ZDZ2SDH & 
          + ZTHFRC(IIB:IIE,JKLEV+1,JKT)*ZDZ1SDH
        END IF
      END DO
    END IF
  END DO 



print*,"!  . 2.2  END CASES , AND TREAT BOUNDARY POINTS"
! -----------------------------------------------------------------
!
! Arrays that will be stored in the moddb_advfrc module
!
  ! Expand arrays to 3D
  XDTHFRC=SPREAD(ZDTHFRC1D,2,IJU)
  XDRVFRC=SPREAD(ZDRVFRC1D,2,IJU)
  !
  ! Fill the first and last point of arrays
  XDTHFRC(1,:,:,:) = XDTHFRC(IIB,:,:,:)
  XDTHFRC(IIU,:,:,:) = XDTHFRC(IIE,:,:,:)
  XDRVFRC(1,:,:,:) =XDRVFRC (IIB,:,:,:)
  XDRVFRC(IIU,:,:,:) = XDRVFRC(IIE,:,:,:)
  !
END DO ! End of loop in time


!-----------------------------------------------------------------------------
!
print*,"!*        3.    PRINT FORCING FIELDS"
!         --------------------------
!

WRITE(UNIT=ILUOUT,FMT='(" THERE ARE ",I2," ADV FORCING AT:")') NADVFRC
DO JL = 1 , NADVFRC
  WRITE(UNIT=ILUOUT,FMT='(F9.0, "s, date:", I3, "/", I3, "/", I5)') &
    TDTADVFRC(JL)%TIME, &
    TDTADVFRC(JL)%TDATE%DAY,   &
    TDTADVFRC(JL)%TDATE%MONTH, &
    TDTADVFRC(JL)%TDATE%YEAR
END DO
!
DO JKT = 2,NADVFRC-1
  IF ( TDTADVFRC(JKT-1) >= TDTADVFRC(JKT) ) THEN
    WRITE(ILUOUT,*) &
      "SET_FRC ERROR: sounding", JKT-1, " is given for a later time than", JKT
    WRITE(ILUOUT,*) &
      "               soundings have to be entered in increasing temporal order"
    WRITE(ILUOUT,*) "SOUNDING TIME ", JKT-1, " IS: "
    WRITE(ILUOUT,*) TDTADVFRC(JKT-1)%TDATE%YEAR,  &
                    TDTADVFRC(JKT-1)%TDATE%MONTH, &
                    TDTADVFRC(JKT-1)%TDATE%DAY,   &
                    TDTADVFRC(JKT-1)%TIME
    WRITE(ILUOUT,*) "SOUNDING TIME ", JKT, " IS: "
    WRITE(ILUOUT,*) TDTADVFRC(JKT)%TDATE%YEAR,  &
                    TDTADVFRC(JKT)%TDATE%MONTH, &
                    TDTADVFRC(JKT)%TDATE%DAY,   &
                    TDTADVFRC(JKT)%TIME
 !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SET_ADVFRC','')
  END IF
END DO

!DEALLOCATE ARRAYS
! 
DEALLOCATE(ZTHVF)
DEALLOCATE(ZTHDF)
DEALLOCATE(ZRVF)

DEALLOCATE(ZPRESS_ADV)
!
DEALLOCATE(ZTHFRC)
DEALLOCATE(ZRVFRC)
DEALLOCATE(ZLAT_FRC)
DEALLOCATE(ZHEIGHTMF)
DEALLOCATE(ZHEIGHTF)
DEALLOCATE(ZTHVUF)
! pour lecture dans PREIDEA
DEALLOCATE(ZZHATM)
DEALLOCATE(ZDRVFRC1D)
DEALLOCATE(ZDTHFRC1D)
DEALLOCATE(ZDVFRC1D)


!-------------------------------------------------------------------------------
END SUBROUTINE SET_ADVFRC
