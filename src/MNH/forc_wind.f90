!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      %Z% Lib:%F%, Version:%I%, Date:%D%, Last modified:%E%
!-----------------------------------------------------------------
!     #########################
      MODULE MODI_FORC_WIND
!     #########################
!
INTERFACE

      SUBROUTINE FORC_WIND(PUT,PWT,PZZ,PTIME)
!
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PUT   ! U
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PWT   ! W
REAL, DIMENSION(:,:,:),   INTENT(IN)       :: PZZ   ! height
REAL,                     INTENT(IN)       :: PTIME
!
END SUBROUTINE FORC_WIND    
!
END INTERFACE
!
END MODULE MODI_FORC_WIND    
!
!
!
!     ############################################
      SUBROUTINE FORC_WIND(PUT,PWT,PZZ,PTIME)
!     ############################################
!
!!
!!    PURPOSE
!!    -------
!
!  this routine provides flow input to the model
!       ON INPUT:              
!                      mx     - x dimension of scalar arrays
!                      mz     - z dimension of scalar arrays
!                      dx,dz  - grid intervals in x and z direction
!                      PTIME   - current PTIME in seconds
!       ON OUTPUT:
!                     ux,uz - rho PTIMEs velocity components normalized
!                             by dx/dt or dz/dt
!
!     AUTHOR:  MARCIN SZUMOWSKI (marcin@ncar.ucar.edu) NOVEMBER 1995
!
!!**  METHOD
!!    ------
!     THIS PROGRZAM SIMULATES AN IDEALIZED FLOW THROUGH A SHALLOW
!     CONVECTIVE CELL. THE FLOW FIELD EVOLVING IN TIME IS DEFINED
!     USING ANALYTICAL FORMULATION FOR THE STREZAMFUNCTION.
!     IT PRODUCES AN UPDRAFT WITH TEMPORALLY EVOLVING MAXIMUM SPEED
!     AND THE TILT OF THE STREZAMLINES THROUGH THE UPDRAFT CORE.
!     THE MAGNITUDE OF THE UPDRAFT MAXIMUM VARIES FROM 1 M/S TO 8 M/S.
!     DEEPER INFLOW AND SHALLOWER OUTFLOW (BELOW THE TRADE WIND
!     INVERSION LEVEL - ZSCALE2) ARE REQUIRED TO PRODUCE AN UPDRAFT
!     MAXIMUM AT THE LEVEL APPROXIMATELY CONSISTENT WITH THE DUAL-
!     DOPPLER DATA. THE TILT IS OBTAINED THROUGH PRESCRIBING A LINEAR 
!     SHEAR IN THE DOMAIN. THE STREZAMLINE TILT VARIES FROM 0 TO 20 DEG.
!     THE RUN LASTS 50 MINUTES WITH A SPIN-UP TIME OF 15 MINUTES WHEN
!     THE UPDRAFT IS STEADY AT 1.0 M/S AND THE SHEAR STEADILY INCREASES
!     TO ITS PEAK VALUE (WHICH PRODUCES STREALINE TILT OF ABOUT 20 DEG).
!     THE LAST 20 MINUTES HAVE STEADY UPDRAFT (2 M/S) AND SHEAR. THIS
!     CORRESPONDS TO THE END OF THE "RAIN-OUT" PROCESS.
!
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!     AUTHOR:  MARCIN SZUMOWSKI (marcin@ncar.ucar.edu) NOVEMBER 1995
!
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_GRID
USE MODD_GRID_n
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PUT   ! U
REAL, DIMENSION(:,:,:),   INTENT(INOUT)    :: PWT   ! W
REAL, DIMENSION(:,:,:),   INTENT(IN)       :: PZZ   ! height
REAL,                     INTENT(IN)       :: PTIME
!
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JI,JK         ! Loop index for the moist variables
INTEGER :: IIE           !
INTEGER :: IIB           !
INTEGER :: IKB           !
INTEGER :: IKE           !
INTEGER :: NXMID    !
!
REAL    :: ZTIME1
REAL    :: ZSCALE1,ZSCALE2
REAL    :: ZAMPL,ZAMPL2
REAL    :: ZAMPL0,ZAMPL20,ZXSCALE,ZXSCALE0,ZAMPA,ZAMPB
REAL    :: ZAMP2A,ZAMP2B,ZTSCALE1,ZTSCALE2
REAL    :: ZTOP,ZXCEN,ZX0,ZZ1,ZZ2,ZXL,ZXX,ZZL,ZSH
!
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZPHI
!
!
!------------------------------------------------------------------------------
!
!
!
IIB = 1
IKB = 1
IIE = SIZE(PUT,1)
IKE = SIZE(PUT,3)-1
ALLOCATE (ZPHI(IIE+1,IKE+2))
!
!
! DETERMINE THE DEPTH OF THE INFLOW (ZSCALE1) AND OUTFLOW (ZSCALE2)
! (in METERS)
!
      ZSCALE1=1.7*1.e3
      ZSCALE2=2.7*1.e3
!
!
!C INITIAL DATA FOR THE STREZAMFUNCTION. ZAMPL IS WMAX IN M/S,
!C ZAMPL2 IS THE VALUE OF LINEAR SHEAR OVER THE ZSCALE2 DEPTH (M/S)
!C ZXSCALE IS WIDTH OF THE UPDRAFT IN M
!C ZAMPA AND ZAMPB CONTROL THE TEMPORAL FLUCTUATIONS OF WMAX 
!C ZAMP2A AND ZAMP2B CONTROL THE TEMPORAL FLUCTUATIONS OF SHEAR (TILT)
!C TSCALE1 AND TSCALE2 ARE PERIODS OF COSINE FLUCTUATIONS OF THE ABOVE
!
      ZAMPL0=1.0
      ZAMPL20=0.
      ZXSCALE0=1.8*1.e3
      ZAMPA=3.5
      ZAMPB=3.0
      ZAMP2A=0.6*1.e3
      ZAMP2B=0.5*1.e3
      ZTSCALE1=600.
      ZTSCALE2=900.
!
!  ZAMPL AND ZAMPL2 VARYING IN TIME
!
        if (PTIME.lt.300.) then
          ZAMPL=ZAMPL0
          ZAMPL2=ZAMPL20
          ZXSCALE=ZXSCALE0
        elseif (PTIME.ge.300..and.PTIME.lt.900.) then
          ZAMPL=ZAMPL0
          ZAMPL2=ZAMP2A*(cos(XPI*((PTIME-300.)/ZTSCALE1 - 1.)) + 1.)
          ZXSCALE=ZXSCALE0
        elseif (PTIME.ge.900..and.PTIME.le.1500.) then
          ZAMPL=ZAMPA*(cos(XPI*((PTIME-900.)/ZTSCALE1 + 1.)) +1.) + ZAMPL0
          ZAMPL2=ZAMP2A*(cos(XPI*((PTIME-300.)/ZTSCALE1 - 1.)) + 1.)
          ZXSCALE=ZXSCALE0
        elseif (PTIME.ge.1500..and.PTIME.lt.2100.) then
          ZAMPL=ZAMPB*(cos(XPI*(PTIME-1500.)/ZTSCALE2) +1.) + 2.*ZAMPL0
          ZAMPL2=ZAMP2B*(cos(XPI*((PTIME-1500.)/ZTSCALE2 - 1.)) + 1.)
          ZXSCALE=ZXSCALE0
        elseif (PTIME.ge.2100..and.PTIME.lt.2400.) then
         ZAMPL=ZAMPB*(cos(XPI*(PTIME-1500.)/ZTSCALE2) +1.) + 2.*ZAMPL0
         ZAMPL2=ZAMP2B*(cos(XPI*((PTIME-1500.)/ZTSCALE2 - 1.)) + 1.)
         ZXSCALE=ZXSCALE0
        else
	 ZTIME1 = 2400.
         ZAMPL=ZAMPB*(cos(XPI*(ZTIME1-1500.)/ZTSCALE2) +1.) + 2.*ZAMPL0
         ZAMPL2=ZAMP2B*(cos(XPI*((ZTIME1-1500.)/ZTSCALE2 - 1.)) + 1.)
         ZXSCALE=ZXSCALE0
        endif
!
! ADOPT THE ZAMPL FOR STREZAMFUNCTION CALCULATION TO BE W_MAX/K_x
      ZAMPL=ZAMPL/XPI*ZXSCALE
!
! DEFINE STREZAMFUNCTION AS A FUNCTION OF HEIGHT
! ALSO, CENTRALIZE THE UPDRAFT TO OCCUPY ONLY THE INNER ZXSCALE
!
      ZTOP=XZHAT(IKE+1)/ZSCALE1
      ZXCEN = 0.5 * XXHAT(IIE)
      NXMID=IIE/2+1
      ZX0=(XXHAT(IIE)-ZXSCALE)/2.
      print *,ZTOP,NXMID,ZX0
  DO JI=IIB,NXMID
    DO JK=IKB,IKE+1
      ZZ1=XZHAT(JK)/ZSCALE1
      ZZ2=XZHAT(JK)/ZSCALE2
      ZXL=2.*SQRT((XXHAT(JI)-ZXCEN)**2)
      ZXX=ZXL/ZXSCALE
      ZZL=2.*ZSCALE1
      ZPHI(JI,JK)=0.
      IF (ZXX.LE.1.) THEN
        IF(ZZ1.LT.1.) THEN
        ZPHI(JI,JK)=-cos(XPI*(XXHAT(JI)-ZX0)/ZXSCALE)*sin(XPI*XZHAT(JK)/ZZL)
        ZPHI(JI,JK)=ZPHI(JI,JK)*ZAMPL
        ELSEIF (ZZ1.GE.1..AND.ZZ2.LE.1.) THEN
        ZZL=2.*(ZSCALE2-ZSCALE1)
        ZPHI(JI,JK)=-cos(XPI*(XXHAT(JI)-ZX0)/ZXSCALE)                      &
                   *sin(XPI*(.5+(XZHAT(JK)-ZSCALE1)/ZZL))
        ZPHI(JI,JK)=ZPHI(JI,JK)*ZAMPL
        ENDIF
      ELSE
        IF(ZZ1.LT.1.) THEN
        ZPHI(JI,JK)=-sin(XPI*XZHAT(JK)/ZZL)
        ZPHI(JI,JK)=ZPHI(JI,JK)*ZAMPL
        ELSEIF (ZZ1.GE.1..AND.ZZ2.LE.1.) THEN
        ZZL=2.*(ZSCALE2-ZSCALE1)
        ZPHI(JI,JK)=-sin(XPI*(.5+(XZHAT(JK)-ZSCALE1)/ZZL))
        ZPHI(JI,JK)=ZPHI(JI,JK)*ZAMPL 
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!
! APPLY THE SYMMETRY CONDITION
      DO JI=NXMID,IIE+1
        DO JK=IKB,IKE+1
          ZPHI(JI,JK)=-ZPHI((IIE+1)-JI+1,JK)
        ENDDO
      ENDDO
!
!     ADD LINEAR SHEAR TO PRODUCE A WEAK TILT OF THE UPDRAFT
      DO JK=IKB,IKE+1
        ZSH=XZHAT(JK)/ZSCALE2
       DO JI=IIB,IIE+1
        ZPHI(JI,JK)=ZPHI(JI,JK) - ZAMPL2*.5*ZSH**2.
       ENDDO
      ENDDO

!c calculate rho*vel by derivation of streamfunction and normalize
!c rho*ux velocity:
      DO JI=IIB,IIE
        DO JK=IKB,IKE
         PUT(JI,:,JK)=-(ZPHI(JI,JK+1)-ZPHI(JI,JK))&
                           /(PZZ(JI,2,JK+1)-PZZ(JI,2,JK))
        ENDDO
      ENDDO
      PUT(:,:,IKE+1) = PUT(:,:,IKE)
!c rho*uz velocity
      DO JK=IKB,IKE+1
        DO JI=IIB,IIE
          PWT(JI,:,JK)=(ZPHI(JI+1,JK)-ZPHI(JI,JK))/XDXHAT(JI) 
        ENDDO
      ENDDO
      DEALLOCATE(ZPHI)
write(*,'(A7,F5.0,A10,F5.2,A10,F5.2)')' ptime=', ptime,' max(PUT)=', maxval(put), ' max(PWT)=', maxval(pwt)
!cc
      RETURN
      END SUBROUTINE FORC_WIND
