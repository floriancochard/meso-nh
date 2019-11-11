!MNH_LIC Copyright 2001-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##########################
      MODULE MODI_REMOVAL_VORTEX
!     ##########################
INTERFACE
!
SUBROUTINE REMOVAL_VORTEX(PZS_LS,PU_LS,PV_LS,PT_LS,PQ_LS,PPS_LS)
!
REAL,DIMENSION(:,:), INTENT(IN)      :: PZS_LS ! Large-scale orography
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PU_LS  ! Pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PV_LS  ! Pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PT_LS  ! Temperature
REAL,DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PQ_LS  ! Specific humidity
REAL,DIMENSION(:,:), OPTIONAL, INTENT(INOUT)   :: PPS_LS ! Pressure (mass points)
!
!
END SUBROUTINE REMOVAL_VORTEX
!
END INTERFACE
!
END MODULE MODI_REMOVAL_VORTEX
!
!
!     ################################################################
      SUBROUTINE REMOVAL_VORTEX(PZS_LS,PU_LS,PV_LS,PT_LS,PQ_LS,PPS_LS)
!     ################################################################
!
!!**  *REMOVAL_VORTEX* - remove an analyzed vortex from large-scale fields
!!
!!    PURPOSE
!!    -------
!      The purpose of this routine is to remove characterictics of a vortex
!       from large-scale analyse fields (ECMWF,Arpege,...) using a
!       filtering technique. This method can be used on any of pronostic
!       variables (wind,temperature,pressure,moisture,...). 
!
!!**  METHOD
!!    ------
!!
!!     Two successive filters are applied: the low-pass filter of Barnes
!!     removes the basic fields from the large-scale analyse, next
!!     the smoothing function of Kurihara (1993) is used to extract the 
!!     analyzed vortex from the total disturbed fields (total fields -
!!     basic fields). Finally, the non-hurricane disturbed fields (total
!!     disturbed fields - analyzed vortex) is added to the basic fields to
!!     obtain the environmental fields. 
!!     
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
!!	O. Nuissier         * L.A *
!!      R. Rogers           * NOAA/AOML/HRD *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original              01/12/01
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_HURR_CONF
USE MODD_CONF, ONLY: NVERB
USE MODD_CST, ONLY: XPI
USE MODD_PARAMETERS, ONLY: JPHEXT,XUNDEF
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_HURR_FIELD_n
USE MODD_DIM_n, ONLY: NIMAX,NJMAX
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_GRID_n, ONLY: XXHAT,XYHAT
!
USE MODE_GRIDPROJ
USE MODE_MSG
!
USE MODI_SHUMAN
USE MODI_BARNES_FILTER
USE MODI_WINDCALC
USE MODI_POLAR_CALC
USE MODI_POLAR_MEAN
USE MODI_EXTRACT_VORTEX
!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments
!
REAL,DIMENSION(:,:),   INTENT(IN)    :: PZS_LS ! Large-scale orography
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PU_LS  ! Pseudo zonal wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PV_LS  ! Pseudo meridian wind component
REAL,DIMENSION(:,:,:), INTENT(INOUT) :: PT_LS  ! Temperature
REAL,DIMENSION(:,:,:), OPTIONAL, INTENT(INOUT) :: PQ_LS  ! Specific humidity
REAL,DIMENSION(:,:), OPTIONAL, INTENT(INOUT)   :: PPS_LS ! Pressure (mass points)
!
!*       0.2   Declarations of local variables
!
INTEGER                                                   :: ILUOUT0
!
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),SIZE(PT_LS,3)) :: ZUMASS,ZVMASS
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),SIZE(PT_LS,3)) :: ZUDIS,ZVDIS
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),SIZE(PT_LS,3)) :: ZTDIS
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),1)             :: ZPDIS
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),SIZE(PT_LS,3)) :: ZQDIS
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),SIZE(PT_LS,3)) :: ZUVORCART,ZVVORCART
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),SIZE(PT_LS,3)) :: ZTVORCART
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),1)             :: ZPVORCART
REAL,DIMENSION(SIZE(PT_LS,1),SIZE(PT_LS,2),SIZE(PT_LS,3)) :: ZQVORCART
INTEGER,DIMENSION(SIZE(PT_LS,3))                          :: IICEN,IJCEN
INTEGER                                                   :: IICEN2,IJCEN2
REAL,DIMENSION(SIZE(PT_LS,3))                             :: ZUR0MOY,ZVR0MOY
REAL,DIMENSION(SIZE(PT_LS,3))                             :: ZTR0MOY
REAL,DIMENSION(1)                                         :: ZPR0MOY
REAL,DIMENSION(SIZE(PT_LS,3))                             :: ZQR0MOY
REAL,DIMENSION(NPHIL)                                     :: ZR01D,ZR01D_EXT
REAL,DIMENSION(NPHIL,SIZE(PT_LS,3))                       :: ZR0WIND,ZR0TEMP
REAL,DIMENSION(NPHIL,1)                                   :: ZR0PRES
REAL,DIMENSION(NPHIL,SIZE(PT_LS,3))                       :: ZR0HUM
REAL,DIMENSION(:,:,:), ALLOCATABLE             :: ZUTCYL
!REAL,DIMENSION(:,:),   ALLOCATABLE             :: ZALTDISCYL
INTEGER          :: IIU,IJU,IP   ! dimensions of arrays
INTEGER          :: JI,JJ,JK     ! Loop indexes along the x,y,z directions
INTEGER          :: II,IJ        ! Indexes of the point for the 1st guess
INTEGER          :: JR,JPHI      ! Loop indexes along radial and azimuth
INTEGER          :: IIMIN,IJMIN
INTEGER          :: IAVGWINDI,IAVGWINDJ
INTEGER          :: IRADMAX0,IRNEWMIN,IRNEWMAX
INTEGER,DIMENSION(NPHIL) :: IRADMAX
REAL             :: ZI,ZJ   ! Fractional indexes of the point for the 1st guess
REAL             :: ZXHAT,ZYHAT
REAL             :: ZSPEED,ZSPEEDMIN
REAL             :: ZGRADVT,ZVTCYLMAX
REAL             :: ZDELTAX,ZDELTAY  ! grid meshes
REAL             :: ZDELTAR
INTEGER          :: IRESP   ! Return code of FM-routines
!
!-------------------------------------------------------------------------------
!
!*	 1. INITIALIZATIONS
!           ---------------
!             
ILUOUT0 = TLUOUT0%NLU
IIU= SIZE(PT_LS,1)
IJU= SIZE(PT_LS,2)
IP= SIZE(PT_LS,3)
ZDELTAX= XXHAT(3) - XXHAT(2)
ZDELTAY= XYHAT(3) - XYHAT(2)
ZDELTAR= MAX(ZDELTAX,ZDELTAY)
!
IF (XLATGUESS == XUNDEF .AND. XLONGUESS == XUNDEF) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','REMOVAL_VORTEX',&
                 'You do not specify a first guess position of the cyclone (XLATGUESS, XLONGUESS)')
END IF
IF (XBOXWIND == XUNDEF) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','REMOVAL_VORTEX',&
                 'ou do not specify a radius to determine the cyclone center inside (XBOXWIND,km)')
END IF
IF (XRADGUESS == XUNDEF) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','REMOVAL_VORTEX',&
                 'You do not specify a guess for the filtered radius (XRADGUESS,km)')
END IF
!
! Allocations
ALLOCATE(XUBASIC(IIU,IJU,IP),XVBASIC(IIU,IJU,IP),XTBASIC(IIU,IJU,IP))
IF (PRESENT(PPS_LS)) ALLOCATE(XPBASIC(IIU,IJU,1))
IF (PRESENT(PQ_LS)) ALLOCATE(XQBASIC(IIU,IJU,IP))
!
! Interpolate U and V onto mass grid
ZUMASS(:,:,:) = MXF(PU_LS(:,:,:))
ZVMASS(:,:,:) = MYF(PV_LS(:,:,:))
!
!-------------------------------------------------------------------------------
!
!*       2. APPLY THE LOW-PASS FILTER OF BARNES
!           -----------------------------------
!
WRITE(ILUOUT0,'(A)')' Begin of REMOVAL_VORTEX routine'
!
CALL BARNES_FILTER(ZUMASS(:,:,:),NK,XLAMBDA,XUBASIC(:,:,:))
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'low-pass filter of Barner applied on U'
!
CALL BARNES_FILTER(ZVMASS(:,:,:),NK,XLAMBDA,XVBASIC(:,:,:))
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'low-pass filter of Barner applied on V'
!
CALL BARNES_FILTER(PT_LS(:,:,:),NK,XLAMBDA,XTBASIC(:,:,:))
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'low-pass filter of Barner applied on T'
!
! disturbance fields
ZUDIS(:,:,:) = ZUMASS(:,:,:) - XUBASIC(:,:,:)
ZVDIS(:,:,:) = ZVMASS(:,:,:) - XVBASIC(:,:,:)
ZTDIS(:,:,:) = PT_LS(:,:,:)  - XTBASIC(:,:,:)
!
IF (PRESENT(PPS_LS)) THEN
  ZPDIS(:,:,1)=PPS_LS(:,:)
  CALL BARNES_FILTER(ZPDIS(:,:,:),NK,XLAMBDA,XPBASIC(:,:,:))
  ZPDIS(:,:,1) = PPS_LS(:,:) - XPBASIC(:,:,1) 
  IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'low-pass filter of Barner applied on P'
END IF
IF (PRESENT(PQ_LS)) THEN
  CALL BARNES_FILTER(PQ_LS(:,:,:),NK,XLAMBDA,XQBASIC(:,:,:))
  ZQDIS(:,:,:) = PQ_LS(:,:,:) - XQBASIC(:,:,:)
  IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'low-pass filter of Barner applied on Q'
END IF
!-------------------------------------------------------------------------------
!
!*       3. LOCATE ANALYSED VORTEX CENTER (dynamical method)
!           -----------------------------
!
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'localizing the position of the fix given'
!
CALL SM_XYHAT(XLATORI,XLONORI,XLATGUESS,XLONGUESS,ZXHAT,ZYHAT)
II=MAX(MIN(COUNT(XXHAT(:)<ZXHAT),IIU-1),1)
IJ=MAX(MIN(COUNT(XYHAT(:)<ZYHAT),IJU-1),1)
ZI=(ZXHAT-XXHAT(II))/(XXHAT(II+1)-XXHAT(II))+FLOAT(II)
ZJ=(ZYHAT-XYHAT(IJ))/(XYHAT(IJ+1)-XYHAT(IJ))+FLOAT(IJ)
IIMIN = INT(ZI)
IJMIN = INT(ZJ)
IF (NVERB>=5) WRITE(ILUOUT0,'(A,I3,A,I3)')' equivalent indexes in the Meso-NH grid: I= ',IIMIN,' J= ',IJMIN
!
IF ( (IIMIN<1) .OR. (IIMIN>IIU+1) .OR. &
     (IJMIN<1) .OR. (IJMIN>IJU+1)      ) THEN
 !callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','REMOVAL_VORTEX','The first guess position of the fix is not in the Meso-NH domain')
END IF
!
IAVGWINDI = NINT(XBOXWIND*1000./ZDELTAX)
IAVGWINDJ = NINT(XBOXWIND*1000./ZDELTAY)
!
DO JK = 1, IP 
  ZSPEEDMIN = 100.
  DO JJ = MAX(IJMIN-IAVGWINDJ,1+JPHEXT),MIN(IJMIN+IAVGWINDJ,IJU-JPHEXT)
    DO JI = MAX(IIMIN-IAVGWINDI,1+JPHEXT),MIN(IIMIN+IAVGWINDI,IIU-JPHEXT)
! center based on dynamical wind method
      ZSPEED = (ZUMASS(JI,JJ,JK)**2. + ZVMASS(JI,JJ,JK)**2.) **0.5
      IF (ZSPEED <= ZSPEEDMIN) THEN
        ZSPEEDMIN = ZSPEED
        IICEN2 = JI
        IJCEN2= JJ
      END IF
    END DO
  END DO
IICEN(JK)=IICEN2
IJCEN(JK)=IJCEN2
!
IF (NVERB>=5) WRITE(ILUOUT0,'(A,I3,A,I3,A,I3)') &
  'dynamical center computed on sigma level ',JK,': ICENDYN= ',IICEN2, ' JCENDYN= ',IJCEN2
END DO
IF (NVERB>=5) WRITE(ILUOUT0,'(A,I3)') 'use dynamical center from sigma level ',NLEVELR0
IICEN(:)=IICEN(NLEVELR0)
IJCEN(:)=IJCEN(NLEVELR0)
!
!-------------------------------------------------------------------------------
!
!*       4. COMPUTE THE TANGENTIEL WIND TO APPLY THE KURIHARA et al. CRITERIA
!           -----------------------------------------------------------------
!
IRADMAX0 = (XRADGUESS*1000. / ZDELTAR) + 1
IRADMAX(:)=IRADMAX0
!
! calculate tangential wind component in model grid
IF (NDIAG_FILT>=2) THEN
  ALLOCATE(XVTDIS(IIU,IJU,IP))
  CALL WINDCALC(ZUDIS,ZVDIS,IICEN,IJCEN,XVTDIS)
END IF
!
! calculate tangential wind component in polar grid
ALLOCATE(ZUTCYL(IRADMAX0,NPHIL,IP))
CALL POLAR_CALC(ZUDIS,IICEN,IJCEN,ZUTCYL,PVAR2IN=ZVDIS,KRADMAX=IRADMAX)
!
! Reverse the signe of tangentiel wind when it is a southern hemisphere
IF (XLATGUESS<0.) THEN
  IF (ASSOCIATED (XVTDIS)) XVTDIS(:,:,:) = -XVTDIS(:,:,:)
  ZUTCYL(:,:,:) = -ZUTCYL(:,:,:)
END IF
!
!ALLOCATE(ZALTDISCYL(IRADMAX0,NPHIL))
!CALL POLAR_CALC(PZS_LS,IICEN(2),IJCEN(2),ZALTDISCYL)
!
!-------------------------------------------------------------------------------
!
!*       5.  COMPUTE THE RADIUS OF THE FILTER DOMAIN
!            ---------------------------------------
!
! using the same criteria than Kurihara et al. (1995)
! The criteria are applied for the grib level sigma=2 (at about 850 hPa)
!
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'Compute the radius of the filter domain'
!
DO JPHI = 1, NPHIL 
  IRNEWMAX = IRADMAX0 
  IRNEWMIN = 2
  ! test on orography
  !DO JR = 1, IRADMAX0
  !   IF (ZALTDISCYL(JR,JPHI) > 1.)   THEN
  !     IRNEWMAX = JR
  !     EXIT
  !   END IF
  !END DO
  ZVTCYLMAX = 0.
  DO JR = 1, IRNEWMAX
    IF (ZUTCYL(JR,JPHI,NLEVELR0) > ZVTCYLMAX) THEN
      ZVTCYLMAX = ZUTCYL(JR,JPHI,NLEVELR0)
      IRNEWMIN = MAX(2,JR)
    END IF
  END DO
  DO JR = IRNEWMIN, IRNEWMAX
    ZGRADVT = (ZUTCYL(JR-1,JPHI,NLEVELR0) - ZUTCYL(JR,JPHI,NLEVELR0))  /  ZDELTAR
    IF (     (ZUTCYL(JR,JPHI,NLEVELR0)<6. .AND. ZGRADVT < 4.E-6 )    &
        .OR. (ZUTCYL(JR,JPHI,NLEVELR0)<3.                       )    ) THEN
      IRNEWMAX = JR - 1
      EXIT
    END IF
  END DO
  ZR01D(JPHI) = IRNEWMAX * ZDELTAR
  ZR01D_EXT(JPHI) = MIN(IRNEWMAX*1.25,IRADMAX(JPHI)*1.) * ZDELTAR
  IF (NVERB>=10) WRITE(ILUOUT0,'(A,I3,A,I3,A,F6.1,A)')' for sector ',JPHI,&
               ' radius is: ',IRNEWMAX,' grid points (',ZR01D(JPHI)/1000.,'km)'
END DO
!
DEALLOCATE(ZUTCYL)
!
!*     5.1 Extend the radius of the filter domain to 1.25 of R0 (except for P)
!
ZR0WIND(:,:) = SPREAD(ZR01D_EXT(:),DIM=2,NCOPIES=IP)
ZR0TEMP(:,:) = ZR0WIND(:,:)
!
ZR0PRES(:,1) = ZR01D(:)
ZR0HUM(:,:) = ZR0WIND(:,:)
!
!-------------------------------------------------------------------------------
!
!*     6. COMPUTING THE MEAN VALUES
!         ------------------------
!
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'Mean value along the periphery of the filter radius'
!
CALL POLAR_MEAN(ZUDIS,ZR0WIND,IICEN,IJCEN,ZUR0MOY)
IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' mean of U computed'
!
CALL POLAR_MEAN(ZVDIS,ZR0WIND,IICEN,IJCEN,ZVR0MOY)
IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' mean of V computed'
!
CALL POLAR_MEAN(ZTDIS,ZR0TEMP,IICEN,IJCEN,ZTR0MOY)
IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' mean of T computed'
!
IF (PRESENT(PPS_LS)) THEN
  CALL POLAR_MEAN(ZPDIS,ZR0PRES,IICEN,IJCEN,ZPR0MOY)
  IF (NVERB>=10) WRITE(ILUOUT0,'(A,F10.2)')' mean of P computed: ',ZPR0MOY(1)
ENDIF
IF (PRESENT(PQ_LS)) THEN
  CALL POLAR_MEAN(ZQDIS,ZR0HUM,IICEN,IJCEN,ZQR0MOY)
  IF (NVERB>=10) WRITE(ILUOUT0,'(A,E10.5)')' mean of Q anomaly computed: ',ZQR0MOY(1)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      7. EXTRACT THE TROPICAL CYCLONE FROM THE TOTAL DISTURBANCE
!          -------------------------------------------------------
!
IF (NVERB>=5) WRITE(ILUOUT0,'(A)')'Separate global vortex from non-hurricane disturbance'
!
CALL EXTRACT_VORTEX(ZUDIS,ZR0WIND,ZUR0MOY,IICEN,IJCEN,ZUVORCART)
IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' U extracted'
!
CALL EXTRACT_VORTEX(ZVDIS,ZR0WIND,ZVR0MOY,IICEN,IJCEN,ZVVORCART)
IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' V extracted'
!
CALL EXTRACT_VORTEX(ZTDIS,ZR0TEMP,ZTR0MOY,IICEN,IJCEN,ZTVORCART)
IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' T extracted'
!
IF (PRESENT(PPS_LS)) THEN
  CALL EXTRACT_VORTEX(ZPDIS,ZR0PRES,ZPR0MOY,IICEN,IJCEN,ZPVORCART)
  IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' P extracted'
ENDIF
IF (PRESENT(PQ_LS)) THEN
  CALL EXTRACT_VORTEX(ZQDIS,ZR0HUM,ZQR0MOY,IICEN,IJCEN,ZQVORCART)
  IF (NVERB>=10) WRITE(ILUOUT0,'(A)')' Q extracted'
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      8. COMPUTE THE ENVIRONMENTAL FIELDS
!          --------------------------------
!
ALLOCATE(XUENV(IIU,IJU,IP),XVENV(IIU,IJU,IP),XTENV(IIU,IJU,IP))
IF (PRESENT(PPS_LS)) ALLOCATE(XPENV(IIU,IJU))
IF (PRESENT(PQ_LS)) ALLOCATE(XQENV(IIU,IJU,IP))
!
WHERE (ZUVORCART(:,:,:) /= XUNDEF)
  XUENV(:,:,:) = XUBASIC(:,:,:) + ZUDIS(:,:,:) -ZUVORCART(:,:,:)
  XVENV(:,:,:) = XVBASIC(:,:,:) + ZVDIS(:,:,:) -ZVVORCART(:,:,:)
  XTENV(:,:,:) = XTBASIC(:,:,:) + ZTDIS(:,:,:) -ZTVORCART(:,:,:)
ELSEWHERE
  XUENV(:,:,:) = XUBASIC(:,:,:) + ZUDIS(:,:,:)
  XVENV(:,:,:) = XVBASIC(:,:,:) + ZVDIS(:,:,:)
  XTENV(:,:,:) = XTBASIC(:,:,:) + ZTDIS(:,:,:)
END WHERE  
IF (PRESENT(PPS_LS)) THEN
  WHERE (ZPVORCART(:,:,1) /= XUNDEF)
    XPENV(:,:)   = XPBASIC(:,:,1) + ZPDIS(:,:,1) -ZPVORCART(:,:,1)
  ELSEWHERE
    XPENV(:,:)   = XPBASIC(:,:,1) + ZPDIS(:,:,1)
  END WHERE
ENDIF
IF (PRESENT(PQ_LS)) THEN
  WHERE (ZQVORCART(:,:,:) /= XUNDEF)
    XQENV(:,:,:) = XQBASIC(:,:,:) + ZQDIS(:,:,:) - ZQVORCART(:,:,:)
  ELSEWHERE
    XQENV(:,:,:) = XQBASIC(:,:,:) + ZQDIS(:,:,:)
  END WHERE
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       9. EPILOG
!            ------
!
!*       9.1 Interpolate back to U and V grid
!
XUBASIC(:,:,:) = MXM(XUBASIC(:,:,:))
XUENV(:,:,:) = MXM(XUENV(:,:,:))
!
XVBASIC(:,:,:) = MYM(XVBASIC(:,:,:))
XVENV(:,:,:) = MYM(XVENV(:,:,:))
!
!*       9.2 
!
IF (NVERB>=10) WRITE(ILUOUT0,'(A)')'Allocate arrays for storing every component computed'
!
ALLOCATE(XUTOT(IIU,IJU,IP),XVTOT(IIU,IJU,IP),XTTOT(IIU,IJU,IP))
IF (PRESENT(PPS_LS)) ALLOCATE(XPTOT(IIU,IJU))
IF (PRESENT(PQ_LS)) ALLOCATE(XQTOT(IIU,IJU,IP))
!
XUTOT(:,:,:) = PU_LS(:,:,:)
XVTOT(:,:,:) = PV_LS(:,:,:)
XTTOT(:,:,:) = PT_LS(:,:,:)
IF (PRESENT(PPS_LS)) XPTOT(:,:)   = PPS_LS(:,:)
IF (PRESENT(PQ_LS)) XQTOT(:,:,:) = PQ_LS(:,:,:) 
!
!*       9.3 Remplace analysed fields by environmental ones
!
! only where cyclonic disturbance has been removed 
WHERE (ZUVORCART(:,:,:) /= XUNDEF)
  PU_LS(:,:,:) = XUENV(:,:,:)
  PV_LS(:,:,:) = XVENV(:,:,:)
  PT_LS(:,:,:) = XTENV(:,:,:)
END WHERE
! elsewhere environmental fields = MXM(MXF(input fields))
!we keep rather input fields themselves
IF (PRESENT(PPS_LS)) THEN
  WHERE (ZPVORCART(:,:,1) /= XUNDEF)
    PPS_LS(:,:)  = XPENV(:,:)
  END WHERE
ENDIF
IF (PRESENT(PQ_LS)) THEN
  WHERE (ZQVORCART(:,:,:) /= XUNDEF)
    PQ_LS(:,:,:)  = XQENV(:,:,:)
  END WHERE
ENDIF
!
WRITE(ILUOUT0,'(A)')' End of REMOVAL_VORTEX routine'
!
!-------------------------------------------------------------------------------!
END SUBROUTINE REMOVAL_VORTEX 
