!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home/cvsroot/MNH-VX-Y-Z/src/MNH/exchange.f90,v $ $Revision: 1.2.2.2.2.2.16.1.2.5.2.1 $ $Date: 2015/12/01 15:26:23 $
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ####################
      MODULE MODI_EXCHANGE
!     ####################
!
INTERFACE
!
!     ##############################################################################
      SUBROUTINE EXCHANGE (PTSTEP,KRR,KSV,PRHODJ,TPFIELDS_ll,                      &
                           PRUS,PRVS,PRWS,PRTHS,PRRS,PRTKES,PRSVS                  )
!     ##############################################################################
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
!
REAL,                     INTENT(IN) :: PTSTEP            !  Time step
INTEGER,                  INTENT(IN) :: KRR               !  Number of water var.
INTEGER,                  INTENT(IN) :: KSV               !  Number of scal. var.
                                                          ! (=1 at the segment beginning)
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ            ! (Rho) dry * Jacobian
TYPE(LIST_ll), POINTER               :: TPFIELDS_ll       ! list of fields to exchange
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS,PRVS,PRWS,   &!
                                           PRTHS,PRTKES       ! Source terms
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS,PRSVS         !
!
END SUBROUTINE EXCHANGE
!
END INTERFACE
!
END MODULE MODI_EXCHANGE
!
!
!
!     #######################################################################
      SUBROUTINE EXCHANGE (PTSTEP,KRR,KSV,PRHODJ,TPFIELDS_ll,               &
                           PRUS,PRVS,PRWS,PRTHS,PRRS,PRTKES,PRSVS           )
!     #######################################################################
!
!!****  * EXCHANGE* - update the halo of each subdomains for the variables at time step t+dt
!!
!!    PURPOSE
!!    -------
!!
!!    The purpose of EXCHANGE is to transform the source terms in the variables at time step t+dt
!!    and update the halo of each subdomains. This routine also takes into account the
!!    cyclic conditions
!
!!**  METHOD
!!    ------
!!    The source term is multipied by twice the time step (except for the first time step)
!!    and divided by ( rhod J ) to obtain the value of the variables at
!!    time step t+dt. The halos of these fields are updated with the values computed by the
!!    neighbor subdomains. Cyclic conditions are treated during this exchange.
!!
!!    EXTERNAL
!!    --------
!!       UPDATE_HALO_ll  :   routine to update the halo
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!          MODD_CONF
!!
!!    REFERENCE
!!    ---------
!!    Book2 of documentation
!!
!!    AUTHOR
!!    ------
!!    P. Jabouille  Meteo France
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    original     18/09/98
!!      05/2006   Remove KEPS
!!      10/2009 (C.Lac) FIT for variables advected by PPM
!!      05/2014 (C.Lac) Correction of negative values of chemical
!!                   tracers moved from ch_monitor to the end of the time step
!!      11/2014 (G.Delautier) Call correction of negative values only if LUSECHEM 
!!      16/02/16 (M. Leriche) conserve total mass for gas phase chem.
!!                   species only, remove negative values without mass 
!!                   corrections for aq. phase and ice phase (lost mass neglig.)
!!      25/08/16 (M.Leriche) remove negative values for aerosols and conserve
!!                   total mass for chemical species in aerosols
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!            
!-----------------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODE_ll
!
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_GRID_n
USE MODD_NSV
USE MODD_BUDGET,      ONLY : LBUDGET_SV
USE MODD_CST,         ONLY : XMNH_TINY
USE MODD_LUNIT_n,     ONLY : TLUOUT
USE MODI_SHUMAN
USE MODI_SUM_ll
USE MODI_BUDGET
USE MODD_CH_MNHC_n,   ONLY : LUSECHEM, LUSECHAQ, LUSECHIC
USE MODD_CH_AEROSOL,  ONLY : LORILAM, NM6_AER
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
!
IMPLICIT NONE
!
!*      0.1  DECLARATIONS OF ARGUMENTS
!
REAL,                     INTENT(IN) :: PTSTEP            !  Time step
INTEGER,                  INTENT(IN) :: KRR               !  Number of water var.
INTEGER,                  INTENT(IN) :: KSV               !  Number of scal. var.
REAL, DIMENSION(:,:,:),   INTENT(IN) :: PRHODJ            ! (Rho) dry * Jacobian
TYPE(LIST_ll), POINTER               :: TPFIELDS_ll       ! list of fields to exchange
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS,PRVS,PRWS,   &!
                                           PRTHS,PRTKES       ! Source terms
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS,PRSVS         !
!
!!
!*      0.2  DECLARATIONS OF LOCAL VARIABLES
!
INTEGER   :: IINFO_ll              ! return code of parallel routine
INTEGER   :: JRR,JSV              ! loop counters
!
INTEGER   :: IKU
INTEGER   :: ILUOUT         ! logical unit numbers of output-listing
INTEGER   :: IRESP          ! IRESP  : return-code if a problem appears
                                    !in LFI subroutines at the open of the file
REAL      :: ZRATIO, ZMASSTOT, ZMASSPOS
!------------------------------------------------------------------------------
!
IKU=SIZE(XZHAT)
ILUOUT = TLUOUT%NLU
!
!*       1.     TRANSFORMS THE SOURCE TERMS INTO PROGNOSTIC VARIABLES
!               -----------------------------------------------------
!
!        1.a Momentum variables
!
PRUS(:,:,:) = PRUS(:,:,:)*PTSTEP / MXM(PRHODJ)
PRVS(:,:,:) = PRVS(:,:,:)*PTSTEP / MYM(PRHODJ)
PRWS(:,:,:) = PRWS(:,:,:)*PTSTEP / MZM(1,IKU,1,PRHODJ)
!
!        1.b Meteorological scalar variables
!
PRTHS(:,:,:) = PRTHS(:,:,:)*PTSTEP/PRHODJ
DO JRR=1,KRR
  PRRS(:,:,:,JRR) = PRRS(:,:,:,JRR)*PTSTEP/PRHODJ
END DO
IF (SIZE(PRTKES,1) /= 0) PRTKES(:,:,:) = PRTKES(:,:,:)*PTSTEP/PRHODJ
!
!        1.c Tracer scalar variables
!
!      REMOVE NEGATIVE VALUES OF CHEM SCALAR
!
IF (LUSECHEM) THEN
  DO JSV=NSV_CHGSBEG,NSV_CHGSEND
    IF ( MIN_ll( PRSVS(:,:,:,JSV), IINFO_ll) < 0.0 ) THEN
!
! compute the total mass 
!
      ZMASSTOT = MAX( 0. , SUM3D_ll( PRSVS(:,:,:,JSV), IINFO_ll ) )
!
! remove the negative values
!
      PRSVS(:,:,:,JSV) = MAX(0., PRSVS(:,:,:,JSV) )
!
! compute the new total mass
!
      ZMASSPOS = MAX(XMNH_TINY,SUM3D_ll( PRSVS(:,:,:,JSV), IINFO_ll ) )
!
! correct again in such a way to conserve the total mass 
!
      ZRATIO = ZMASSTOT / ZMASSPOS
      PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) * ZRATIO
!
      WRITE(ILUOUT,*)'CHEMISTRY',JSV,'HAS NEGATIVE VALUES'
      WRITE(ILUOUT,*)'GAS SOURCES IS CORRECTED BY RATIO',ZRATIO
    END IF
  END DO
  IF (LUSECHAQ) THEN
    DO JSV =  NSV_CHACBEG, NSV_CHACEND
      IF ( MIN_ll( PRSVS(:,:,:,JSV), IINFO_ll) < 0.0 ) THEN
! remove the negative values
        PRSVS(:,:,:,JSV) = MAX(0., PRSVS(:,:,:,JSV) )
        WRITE(ILUOUT,*)'CHEMISTRY',JSV,'HAS NEGATIVE VALUES'
        WRITE(ILUOUT,*)'CLOUD OR RAIN SOURCE IS SET TO ZERO'
      END IF
    END DO
  ENDIF
!
  IF (LUSECHIC) THEN
    DO JSV =  NSV_CHICBEG, NSV_CHICEND
      IF ( MIN_ll( PRSVS(:,:,:,JSV), IINFO_ll) < 0.0 ) THEN
! remove the negative values
        PRSVS(:,:,:,JSV) = MAX(0., PRSVS(:,:,:,JSV) )
        WRITE(ILUOUT,*)'CHEMISTRY',JSV,'HAS NEGATIVE VALUES'
        WRITE(ILUOUT,*)'ICE PHASE SOURCE IS SET TO ZERO'
      END IF
    END DO
  ENDIF
!
  IF (LBUDGET_SV) THEN
    DO JSV=NSV_CHEMBEG,NSV_CHEMEND
      CALL BUDGET(PRSVS(:,:,:,JSV),JSV+12,'NEGA_BU_RSV')
    ENDDO
  ENDIF
!
! aerosol sv
  IF (LORILAM) THEN
    DO JSV=NSV_AERBEG,NSV_AEREND-2-NM6_AER ! keep chem. species only
      IF ( MIN_ll( PRSVS(:,:,:,JSV), IINFO_ll) < 0.0 ) THEN
!
! compute the total mass
!
        ZMASSTOT = MAX( 0. , SUM3D_ll( PRSVS(:,:,:,JSV), IINFO_ll ) )
!
! remove the negative values
!
        PRSVS(:,:,:,JSV) = MAX(0., PRSVS(:,:,:,JSV) )
!
! compute the new total mass
!
        ZMASSPOS = MAX(XMNH_TINY,SUM3D_ll( PRSVS(:,:,:,JSV), IINFO_ll ) )
!
! correct again in such a way to conserve the total mass 
!
        ZRATIO = ZMASSTOT / ZMASSPOS
        PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV) * ZRATIO
!
        WRITE(ILUOUT,*)'CHEMISTRY',JSV,'HAS NEGATIVE VALUES'
        WRITE(ILUOUT,*)'AP SOURCES IS CORRECTED BY RATIO',ZRATIO
      END IF
    END DO
!
    DO JSV=NSV_AEREND-2-NM6_AER+1,NSV_AEREND
      IF ( MIN_ll( PRSVS(:,:,:,JSV), IINFO_ll) < 0.0 ) THEN
! remove the negative values for M0 and M6
         PRSVS(:,:,:,JSV) = MAX(0., PRSVS(:,:,:,JSV) )
         WRITE(ILUOUT,*)'CHEMISTRY',JSV,'HAS NEGATIVE VALUES'
         WRITE(ILUOUT,*)'AP MOMENT SOURCES IS SET TO ZERO'
      END IF
    END DO
    IF (LBUDGET_SV) THEN
      DO JSV=NSV_AERBEG,NSV_AEREND
        CALL BUDGET(PRSVS(:,:,:,JSV),JSV+12,'NEGA_BU_RSV')
      ENDDO
    ENDIF
  ENDIF
ENDIF
!
DO JSV=1,KSV
  PRSVS(:,:,:,JSV) = PRSVS(:,:,:,JSV)*PTSTEP/PRHODJ
END DO
!
IF(LBLOWSNOW) THEN
   DO JSV=1,(NBLOWSNOW_2D)
       XRSNWCANOS(:,:,JSV) = XRSNWCANOS(:,:,JSV)*PTSTEP/PRHODJ(:,:,1)
   END DO
END IF
!
!
!------------------------------------------------------------------------------
!
!*      2      UPDATE THE FIRST LAYER OF THE HALO
!              ----------------------------------
!
CALL UPDATE_HALO_ll(TPFIELDS_ll, IINFO_ll)
!
END SUBROUTINE EXCHANGE
