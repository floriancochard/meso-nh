!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!##################
 MODULE MODI_BUDGET
!##################
!
INTERFACE
!
SUBROUTINE BUDGET(PVARS,KBUDN,HBUVAR)
!
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS    ! Source 
INTEGER               , INTENT(IN) :: KBUDN    ! variable number
CHARACTER (LEN=*)    , INTENT(IN) :: HBUVAR   ! Identifier of the Budget of the
                                               ! variable that is considered 
!
END SUBROUTINE BUDGET
!
END INTERFACE
!
END MODULE MODI_BUDGET
!     #####################################
      SUBROUTINE BUDGET(PVARS,KBUDN,HBUVAR)
!     #####################################
!
!!****  *BUDGET* - routine to call the BUDGET routine. 
!!                           
!!
!!    PURPOSE
!!    -------
!        This routine selects the variable RVAR, the budget of which is 
!     processed in the inner routine BUDGET_CASE.  !
!!**  METHOD
!!    ------
!!       
!!     
!!
!!    EXTERNAL
!!    --------
!!      CART_COMPRESS 
!!      MASK_COMPRESS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!         XBURU       : budget array of the variable RU
!!         XBURV       : budget array of the variable RV
!!         XBURW       : budget array of the variable RW
!!         XBURTH      : budget array of the variable RTH
!!         XBURTKE     : budget array of the variable RTKE
!!         XBURRV      : budget array of the variable RRV
!!         XBURRC      : budget array of the variable RRC
!!         XBURRR      : budget array of the variable RRR
!!         XBURRI      : budget array of the variable RRI
!!         XBURRS      : budget array of the variable RRS
!!         XBURRG      : budget array of the variable RRG
!!         XBURRH      : budget array of the variable RRH
!!         XBURTKE     : budget array of the variable RTKE
!!         XBURSV(x)   : budget array of the variable RSVx
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!  	J. Nicolau       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/08/94
!!      J. Stein    26/06/96  add the 'OF','NO' option  
!!      J.-P. Pinty 12/12/96  simplifies the coding
!!      V. Masson   06/10/02  add LES budgets
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BUDGET
USE MODD_LUNIT
!USE MODD_CONF_n
USE MODD_CONF, ONLY : LCHECK
USE MODD_NSV,  ONLY : NSV
USE MODD_LES
!
USE MODE_FM
USE MODE_IO_ll
USE MODE_MSG
!
USE MODI_LES_BUDGET
USE MODI_CART_COMPRESS
USE MODI_MASK_COMPRESS
!
USE MODE_MPPDB
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments :
!
INTEGER               , INTENT(IN) :: KBUDN    ! variable number
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVARS    ! source of the variable 
CHARACTER (LEN=*)     , INTENT(IN) :: HBUVAR   ! Identifier of the Budget of the
                                               ! variable that is considered 
INTEGER  :: IBUSV   ! Index of the SV 
!
INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
INTEGER  :: IRESP   ! Return code of FM-routines
!
REAL     :: ZTIME1  ! CPU time counter
REAL     :: ZTIME2  ! CPU time counter
!
REAL     :: XPRECISION ! for reproductibility checks

!-------------------------------------------------------------------------------
!
!* Reproductivity checks
!  Warning: requires an adaptation of the makefile in order to run two runs in
!  parallel for comparison
!
XPRECISION = 1E-10
IF (LCHECK) THEN
  print*,'BUDGET :',HBUVAR
  CALL MPPDB_CHECK3D(PVARS,HBUVAR,XPRECISION)
END IF
!
!
!* call to LES budgets
!
IF (LLES_CALL) CALL LES_BUDGET(PVARS,KBUDN,HBUVAR)
!
!* call to prognostic variables budgets
!
IF (.NOT. LBU_ENABLE) RETURN
!
SELECT CASE (KBUDN)
  CASE (1) 
    IF (.NOT. LBU_RU) RETURN 
  CASE (2) 
    IF (.NOT. LBU_RV) RETURN 
  CASE (3) 
    IF (.NOT. LBU_RW) RETURN 
  CASE (4) 
    IF (.NOT. LBU_RTH) RETURN 
  CASE (5) 
    IF (.NOT. LBU_RTKE) RETURN 
  CASE (6) 
    IF (.NOT. LBU_RRV) RETURN 
  CASE (7) 
    IF (.NOT. LBU_RRC) RETURN 
  CASE (8) 
    IF (.NOT. LBU_RRR) RETURN 
  CASE (9) 
    IF (.NOT. LBU_RRI) RETURN 
  CASE (10) 
    IF (.NOT. LBU_RRS) RETURN 
  CASE (11) 
    IF (.NOT. LBU_RRG) RETURN 
  CASE (12) 
    IF (.NOT. LBU_RRH) RETURN 
  CASE (13:) 
    IF (.NOT. LBU_RSV) RETURN 
END SELECT
!
!-------------------------------------------------------------------------------
!
CALL SECOND_MNH(ZTIME1)
!
SELECT CASE (KBUDN)
!
  CASE (1)  !            ==>  RU BUDGET
    CALL BUDGET_CASE(XBURU)
!
  CASE (2)  !            ==>  RV BUDGET
    CALL BUDGET_CASE(XBURV)
!
  CASE (3)  !            ==>  RW BUDGET
    CALL BUDGET_CASE(XBURW)
!
  CASE (4)  !            ==>  RTH BUDGET
    CALL BUDGET_CASE(XBURTH)
!
  CASE (5)  !            ==>  RTKE BUDGET
    CALL BUDGET_CASE(XBURTKE)
!
  CASE (6)  !            ==>  RRV BUDGET
    CALL BUDGET_CASE(XBURRV)
!
  CASE (7)  !            ==>  RRC BUDGET
    CALL BUDGET_CASE(XBURRC)
!
  CASE (8)  !            ==>  RRR BUDGET
    CALL BUDGET_CASE(XBURRR)
!
  CASE (9)  !            ==>  RRI BUDGET
    CALL BUDGET_CASE(XBURRI)
!
  CASE (10) !            ==>  RRS BUDGET
    CALL BUDGET_CASE(XBURRS)
!
  CASE (11) !            ==>  RRG BUDGET
    CALL BUDGET_CASE(XBURRG)
!
  CASE (12) !            ==>  RRH BUDGET
    CALL BUDGET_CASE(XBURRH)
!
  CASE (13:)!            ==>  RSVx BUDGET
    IBUSV = KBUDN - 12
    IF( IBUSV <= NSV ) THEN 
      CALL BUDGET_CASE(XBURSV(:,:,:,:,IBUSV))
    ELSE
      ILUOUT0 = TLUOUT0%NLU
      WRITE(UNIT=ILUOUT0,FMT='("BUDGET: SCALAR VARIABLE",I2," IS ABSENT !!")') &
                                IBUSV
      WRITE(UNIT=ILUOUT0,FMT='("CHECK FOR THE CALL BUDGET OF THAT VARIABLE")')
!callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','BUDGET','')
    END IF
END SELECT
!
CALL SECOND_MNH(ZTIME2)
!
XTIME_BU_PROCESS = XTIME_BU_PROCESS + ZTIME2 - ZTIME1
XTIME_BU = XTIME_BU + ZTIME2 - ZTIME1
!
!----------------------------------------------------------------------
CONTAINS
!----------------------------------------------------------------------
!     ###############################
      SUBROUTINE BUDGET_CASE(PBURVAR)
!     ###############################
!
!!****  *BUDGET_CASE* - routine to call the BUDGET_CASE routine. 
!!                           
!!
!!    PURPOSE
!!    -------
!        This routine chooses the right call to the functions CART_COMPRESS
!     or MASK_COMPRESS (which realize the compression of the source PVARS
!     in the different directions) and achieves in function of HACTION (which
!     determines the operations to be executed) the budget for the variable 
!     corresponding to the number KBUDN. The budget process counter is
!     incremented by NBUINC depending on the number of active processes in the 
!     model.
!
!!**  METHOD
!!    ------
!!       
!!     
!!
!!    EXTERNAL
!!    --------
!!      CART_COMPRESS 
!!      MASK_COMPRESS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!         CBUACTION   : type of operation
!!         CBUTYPE     : budget type (CART,MASK or NONE)
!!         NBUTIME     : number of the budget step
!!         NBUPROCCTR  : process counter for each budget variable
!!         PBURVAR     : budget array of the variable RVAR
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/12/96
!!      Modification 24/06/99 N. Asencio  : budget // , the dimensions of the
!!                                          budget arrays are implicit
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
  USE MODI_CART_COMPRESS
  USE MODI_MASK_COMPRESS
!
  IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments :
!
 REAL, DIMENSION(:,:,:,:), INTENT(INOUT):: PBURVAR  ! budget of variable RVAR
!
!*       0.2   Declarations of local variables :
  CHARACTER (LEN=99) ::   YBUVAR_ADJUSTED           ! Adjusted string
  CHARACTER (LEN=99) ::   YBUCOMMENT_ADJUSTED       ! Adjusted string
  CHARACTER (LEN=99) ::   YBUVAR                    ! local string
  CHARACTER (LEN=99) ::   YBUCOMMENT                ! local string

  INTEGER            ::   ILEN                      ! Number of non-blank char.
!
!
!*       1.     SECURITY TEST
!               -------------
!
  YBUVAR      =   HBUVAR
  YBUCOMMENT  =   CBUCOMMENT(KBUDN,NBUPROCCTR(KBUDN))
  YBUVAR_ADJUSTED     = ADJUSTR(YBUVAR)
  YBUCOMMENT_ADJUSTED = ADJUSTR(YBUCOMMENT)
  ILEN =  LEN_TRIM( ADJUSTL(YBUVAR))
!
  IF( CBUACTION(KBUDN,NBUCTR_ACTV(KBUDN))/='NO'.AND. &
      CBUACTION(KBUDN,NBUCTR_ACTV(KBUDN))/='OF'.AND. &
      CBUACTION(KBUDN,NBUCTR_ACTV(KBUDN))/='CC'     ) THEN
    IF( YBUVAR_ADJUSTED(100-ILEN:99) /= YBUCOMMENT_ADJUSTED(100-ILEN:99) &
                                             .OR. ILEN==0 ) THEN
      ILUOUT0 = TLUOUT0%NLU
      WRITE(UNIT=ILUOUT0,FMT='("BUDGET: WRONG BUDGET IDENTIFICATION !!")')
      WRITE(UNIT=ILUOUT0,FMT='("BUDGET: PRESENT  VARIABLE: ",I2)') KBUDN
      WRITE(UNIT=ILUOUT0,FMT='("BUDGET: PRESENT  IDENTIFIER: ",A99)') &
                                    YBUVAR_ADJUSTED
      WRITE(UNIT=ILUOUT0,FMT='("BUDGET: EXPECTED IDENTIFIER: ",A99)') &
                            YBUCOMMENT_ADJUSTED
      WRITE(UNIT=ILUOUT0,FMT='("PLEASE CHECK THE CALL BUDGET OF THE VARIABLE")')
      WRITE(UNIT=ILUOUT0,FMT='("AND THE BUDGET PROCESS ORDER IN INI_BUDGET !")')
!callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','BUDGET','')
    END IF
  END IF
!
! Budget integration in case of successful test
!
  SELECT CASE (CBUTYPE)
!
!*	     2.     "CART" CASE
!               -----------
!
    CASE ('CART')
!
      SELECT CASE (CBUACTION(KBUDN,NBUCTR_ACTV(KBUDN)))
!
!*	     2.1    Budget beginning : initial fields
!               filled in budget tabulars (NBUPROCCTR=1)
!
        CASE('IG')            
          PBURVAR(:,:,:,1)=CART_COMPRESS(PVARS)
!
!*	     2.2    average tendancy filled every time
!               step in budget tabulars (NBUPROCCTR=3)
!            
        CASE('ES')          
          PBURVAR(:,:,:,3)=PBURVAR(:,:,:,3)+CART_COMPRESS(PVARS)/NBUSTEP
!
!*    	 2.3    Cumul of the sources 
!
        CASE('CC')
          PBURVAR(:,:,:,2)=CART_COMPRESS(PVARS)
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
!
!*	     2.4    Difference in order to compute the budget
!                   for the process NBUPROCCTR                 
!
        CASE('DD')
          PBURVAR(:,:,:,NBUPROCCTR(KBUDN))= PBURVAR(:,:,:,NBUPROCCTR(KBUDN)) &
                                          + CART_COMPRESS(PVARS)             &
                                          - PBURVAR(:,:,:,2)          
          NBUPROCCTR(KBUDN)=NBUPROCCTR(KBUDN)+1
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
!
!*	     2.5    Difference in order to compute the budget for the
!               process NBUPROCCTR and Cumul of the sources (NBUPROCCTR=2)
!
        CASE('DC')
          PBURVAR(:,:,:,NBUPROCCTR(KBUDN)) = PBURVAR(:,:,:,NBUPROCCTR(KBUDN))&
                                           + CART_COMPRESS(PVARS)            &
                                           - PBURVAR(:,:,:,2)          
          PBURVAR(:,:,:,2)=CART_COMPRESS(PVARS)
          NBUPROCCTR(KBUDN)=NBUPROCCTR(KBUDN)+1
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
        CASE('NO')
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
        CASE('OF')
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
          RETURN
      END SELECT
!        
!*	     3.    "MASK" CASE
!               -----------
!
    CASE ('MASK')
!
      SELECT CASE (CBUACTION(KBUDN,NBUCTR_ACTV(KBUDN)))            
!
!*	     3.1    Budget beginning : initial fields
!               filled in budget tabulars (NBUPROC=1)
!
        CASE('IG')
          PBURVAR(:,NBUTIME,:,1) = MASK_COMPRESS(PVARS)
!
!*	     3.2    average tendancy filled every time
!                 step in budget tabulars (NBUPROCCTR=3)
!    
        CASE('ES')      
          PBURVAR(:,NBUTIME,:,3) = PBURVAR(:,NBUTIME,:,3)   &
                                 + MASK_COMPRESS(PVARS)/NBUSTEP
!
!*	     3.3    Cumul of the sources 
!
        CASE('CC')
          PBURVAR(:,NBUTIME,:,2)=MASK_COMPRESS(PVARS)
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
!
!*	     3.4    Difference in order to compute the budget
!               for the process NBUPROCCTR                 
!
        CASE('DD')
          PBURVAR(:,NBUTIME,:,NBUPROCCTR(KBUDN))                      &
                             = PBURVAR(:,NBUTIME,:,NBUPROCCTR(KBUDN)) &
                             + MASK_COMPRESS(PVARS)                   &
                             - PBURVAR(:,NBUTIME,:,2)          
          NBUPROCCTR(KBUDN)=NBUPROCCTR(KBUDN)+1
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)              &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
!
!*       3.5    Difference in order to compute the budget for the
!               process NBUPROCCTR and Cumul of the sources (NBUPROCCTR=2)
!
        CASE('DC')
          PBURVAR(:,NBUTIME,:,NBUPROCCTR(KBUDN))                      &
                             = PBURVAR(:,NBUTIME,:,NBUPROCCTR(KBUDN)) &
                                               +MASK_COMPRESS(PVARS)  &
                                               -PBURVAR(:,NBUTIME,:,2)
          PBURVAR(:,NBUTIME,:,2)=MASK_COMPRESS(PVARS)
          NBUPROCCTR(KBUDN)=NBUPROCCTR(KBUDN)+1
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
        CASE('NO')
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
        CASE('OF')
!
! advance the process counter
!
          NBUCTR_ACTV(KBUDN) = NBUCTR_ACTV(KBUDN)             &
                             + NBUINC(KBUDN,NBUCTR_ACTV(KBUDN))
          RETURN
      END SELECT          
  END SELECT
!
  END SUBROUTINE BUDGET_CASE
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE BUDGET
