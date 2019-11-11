!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source: /home//MESONH/MNH-V4-6-5/src/SRC_CHIMAQ/ch_set_ph.f90
!-----------------------------------------------------------------
!!    #######################
      MODULE MODI_CH_SET_PH
!!    #######################
!!
INTERFACE
!!
SUBROUTINE CH_SET_PH(KLUOUT,KMI,KVECNPT,PCONC,PPH,KVERB,KRR)
!!
IMPLICIT NONE
!!
INTEGER, INTENT(IN)                      :: KLUOUT
INTEGER, INTENT(IN)                      :: KMI
INTEGER, INTENT(IN)                      :: KVECNPT
REAL,    INTENT(IN),    DIMENSION(:,:)   :: PCONC
REAL,    INTENT(INOUT), DIMENSION(:)     :: PPH
INTEGER, INTENT(IN)                      :: KVERB
INTEGER, INTENT(IN)                      :: KRR
!!
!!
END SUBROUTINE CH_SET_PH
!!
END INTERFACE
!!
END MODULE MODI_CH_SET_PH
!!
!!    ############################################################
      SUBROUTINE CH_SET_PH(KLUOUT,KMI,KVECNPT,PCONC,PPH,KVERB,KRR)
!!    ############################################################
!!
!!*** *CH_SET_PH*
!!
!!    PURPOSE
!!    -------
!        drives the pH calculation
!!
!!**  METHOD
!!    ------
!!       call cloud mask and, then call subroutine resolving pH, which
!!    only solve pH on points where the liquid water content is superior to
!!    a threshold value as selected by the cloud mask subroutine. Return
!!    vector to ch_set_rates_n and, then to ch_monitor_n.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    M. Leriche   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 05/06/07
!!    M. Leriche 11/07/07 extension to the rain drops
!!    Juan 24/09/2012: for BUG Pgi rewrite PACK function on mode_pack_pgi
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODI_CH_CLOUD_MASK
USE MODI_CH_RAIN_MASK
USE MODI_CH_SOLVE_PH
!
USE MODD_CH_M9_n,      ONLY: NEQ, NEQAQ
USE MODD_CH_M9_SCHEME, ONLY: CCSTYPE, TACCS
!
#ifdef MNH_PGI
USE MODE_PACK_PGI
#endif
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER, INTENT(IN)                      :: KLUOUT
INTEGER, INTENT(IN)                      :: KMI
INTEGER, INTENT(IN)                      :: KVECNPT
REAL,    INTENT(IN),    DIMENSION(:,:)   :: PCONC
REAL,    INTENT(INOUT), DIMENSION(:)     :: PPH
INTEGER, INTENT(IN)                      :: KVERB
INTEGER, INTENT(IN)                      :: KRR
!
!*      0.2    declarations local variables
!
TYPE(CCSTYPE), POINTER :: TZK
!
LOGICAL, DIMENSION(KVECNPT)  :: GLW
        ! cloud mask for pH and aqueous phase chemistry
INTEGER :: ILW !number of points cloud mask
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCHEM_CONC
REAL, DIMENSION(:),   ALLOCATABLE :: ZCHEM_T
REAL, DIMENSION(:),   ALLOCATABLE :: ZCHEM_LW
REAL, DIMENSION(:),   ALLOCATABLE :: ZCHEM_PH
!
INTEGER :: JL
!
!-------------------------------------------------------------------------------
!
TZK=>TACCS(KMI) ! to get the temperature and the lwc
!
SELECT CASE (KRR)
  CASE(2) ! cloud water 
    CALL CH_CLOUD_MASK(KMI,KVECNPT,TZK%LWC,GLW,ILW)
  CASE(3) ! rain water 
    CALL CH_RAIN_MASK(KMI,KVECNPT,TZK%LWR,GLW,ILW)
  CASE DEFAULT
    WRITE(UNIT=KLUOUT,FMT='("CH_SET_PH: bad setting for KRR")')
    RETURN
END SELECT
!
IF (ILW>=1) THEN
  ALLOCATE(ZCHEM_CONC(ILW,NEQAQ/2))
  ALLOCATE(ZCHEM_T(ILW))
  ALLOCATE(ZCHEM_LW(ILW))
  ALLOCATE(ZCHEM_PH(ILW))
  ZCHEM_PH(:) = PACK( PPH(:), MASK=GLW(:) )
  ZCHEM_T(:)  = PACK( TZK%T, MASK=GLW(:) )
  SELECT CASE (KRR)
    CASE(2) ! cloud water 
      DO JL = 1, NEQAQ/2
        ZCHEM_CONC(:,JL) = PACK( PCONC(:,NEQ-NEQAQ+JL), MASK=GLW(:) )
      END DO
      ZCHEM_LW(:)   = PACK( TZK%LWC, MASK=GLW(:) )
    CASE(3) ! rain water 
      DO JL = 1, NEQAQ/2
        ZCHEM_CONC(:,JL) = PACK( PCONC(:,NEQ-NEQAQ/2+JL), MASK=GLW(:) )
      END DO
      ZCHEM_LW(:)   = PACK( TZK%LWR, MASK=GLW(:) )
    CASE DEFAULT
      WRITE(UNIT=KLUOUT,FMT='("CH_SET_PH: bad setting for KRR")')
      RETURN
  END SELECT
!
  CALL CH_SOLVE_PH(KLUOUT,ZCHEM_CONC,ZCHEM_T,ZCHEM_LW,ILW,ZCHEM_PH,KRR)
!
  PPH(:) = UNPACK( ZCHEM_PH(:), MASK=GLW(:), FIELD=0.0 )
!
  DEALLOCATE(ZCHEM_CONC)
  DEALLOCATE(ZCHEM_T)
  DEALLOCATE(ZCHEM_LW)
  DEALLOCATE(ZCHEM_PH)
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CH_SET_PH
