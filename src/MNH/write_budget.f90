!MNH_LIC Copyright 1995-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!######################## 
 MODULE MODI_WRITE_BUDGET
!########################
!
INTERFACE
!
      SUBROUTINE WRITE_BUDGET(TPDIAFILE,TPDTCUR, &
                              TPDTMOD,PTSTEP, KSV)
!
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODD_TYPE_DATE
!
TYPE(TFILEDATA),    INTENT(IN) :: TPDIAFILE    ! file to write
TYPE (DATE_TIME),   INTENT(IN) :: TPDTCUR      ! Current date and time
TYPE (DATE_TIME),   INTENT(IN) :: TPDTMOD      ! Creation date and time
REAL,               INTENT(IN) :: PTSTEP       ! time step
INTEGER,            INTENT(IN) :: KSV          ! Number of Scalar Variables
!
END SUBROUTINE WRITE_BUDGET  
!
END INTERFACE
!
END MODULE MODI_WRITE_BUDGET
!
!
!
!     ############################################
      SUBROUTINE WRITE_BUDGET(TPDIAFILE,TPDTCUR, &
                              TPDTMOD,PTSTEP, KSV)
!     ############################################
!
!!****  *WRITE_BUDGET* - routine to write a LFIFM file for the budget. 
!!                           
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to write an initial LFIFM File 
!     of name YFILEDIA//'.lfi' with the FM routines. This routine is 
!     temporary because the budget terms had to be stored in the diachronic
!     MesoNH-files, not yet developped. 
!
!!**  METHOD
!!    ------
!!      The data are written in the LFIFM file :
!!        - dimensions
!!        - budget arrays
!!        - tracer array in mask case
!!
!!      The localization on the model grid is also indicated :
!!
!!        IGRID = 1 for mass grid point
!!        IGRID = 2 for U grid point
!!        IGRID = 3 for V grid point
!!        IGRID = 4 for w grid point
!!        IGRID = 0 for meaningless case
!!            
!!      
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!       Module MODD_BUDGET
!!         
!!         CBUTYPE     : Budget type (CART,MASK,SKIP or NONE)
!!         CBURECORD   : name of output recording files for the budgets
!!         CBUCOMMENT  : name of a process for a budget
!!         NBUPROCNBR  : number of processes for each variable
!!         NBUTIME     : number of the budget time intervals ('MASK' case)
!!         NBUWRNB      : number of budget steps when the budget is written
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
!!         XBURSV      : budget array of the variable RSVx
!!         
!!
!!    REFERENCE
!!    ---------
!!      Book2 of MESO-NH documentation (routine WRITE_BUDGET)
!!
!!
!!    AUTHOR
!!    ------
!!	J. Nicolau       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     27/02/95
!!      J. Stein     9/9/96     add the writings in the diachronic file 
!!      J.-P. Pinty  18/12/96   clarify the coding
!!      J.-P. Pinty  18/03/97   correction for the SVx
!!      V. Gouget M. Chong J.-P. Lafore  add the BURHODJ, TSTEP and BULEN
!!                   10/02/98              and writes in physical units
!!      V. Ducrocq   07/06/99   //
!!      N. Asencio   18/06/99  // budget with MASK case 
!!                             delete ZTORE arrays no longer used, so delete
!!                             KIU,KJU,KKU arguments
!!                             the mask is written once with a FMWRIT call outside
!!                             write_diachro: its name is MASK_(value of NBUTSHIFT).MASK
!!                             MENU_DIACHRO must be called after FMWRIT to be read in
!!                             read_diachro.
!!                             NBUTSHIFT is incremented at the beginning of the routine
!!                             The dimensions of the XBUR.. arrays are : first one
!!                             is the dimension along K, second one is the time, the
!!                             third one is the number of the masks.
!!      October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!      24/03/2014  (J.Escobar ) miss placed deallocate in RSV budget
!!      11/09/2015  (C.Lac)    Correction due to FIT temporal scheme
!!      28/03/2018  (P.Wautelet) Replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    
!              ------------
USE MODD_BUDGET
USE MODD_IO_ll,   ONLY: TFILEDATA
USE MODD_LUNIT_n, ONLY: TLUOUT
!
USE MODE_DATETIME
USE MODE_FIELD,   ONLY: TFIELDDATA,TYPEREAL
USE MODE_FMWRIT
USE MODE_TIME
!
USE MODI_END_CART_COMPRESS
USE MODI_END_MASK_COMPRESS
USE MODI_MENU_DIACHRO
USE MODI_WRITE_DIACHRO
!
!
IMPLICIT NONE
!  
!  
!*       0.1   Declarations of arguments :
!
TYPE(TFILEDATA),    INTENT(IN) :: TPDIAFILE    ! file to write
TYPE (DATE_TIME),   INTENT(IN) :: TPDTCUR      ! Current date and time
TYPE (DATE_TIME),   INTENT(IN) :: TPDTMOD      ! Creation date and time
REAL,               INTENT(IN) :: PTSTEP       ! time step
INTEGER,            INTENT(IN) :: KSV          ! Number of Scalar Variables
!  
!*       0.2   Declarations of local variables :
!
CHARACTER(LEN=NMNHNAMELGTMAX) :: YRECFM        ! Name of the article to be written
INTEGER           :: JT,JPROC,JMASK
!
!
REAL, ALLOCATABLE  , DIMENSION(:,:,:,:,:,:) :: ZWORK, ZWORKT,  ZWORKMASK  ! local array 
            ! conformal to what is asked by the diachro format for the fields
            ! and for the masks
LOGICAL :: GNOCOMPRESS !  If  TRUE : no compress along x and y direction in the CART option
REAL,    ALLOCATABLE              , DIMENSION(:)  :: ZCONVERT      ! unit conversion coefficient
REAL,    ALLOCATABLE              , DIMENSION(:,:):: ZWORKTEMP     ! time
INTEGER, ALLOCATABLE              , DIMENSION(:)  :: IWORKGRID     ! grid label
CHARACTER (LEN=99),  ALLOCATABLE  , DIMENSION(:)  :: YBUCOMMENT    ! comment   
CHARACTER (LEN=100), ALLOCATABLE  , DIMENSION(:)  :: YWORKCOMMENT  ! comment   
CHARACTER (LEN=100), ALLOCATABLE  , DIMENSION(:)  :: YWORKUNIT     ! comment
CHARACTER (LEN=9)                                 :: YGROUP_NAME   ! group name                                    
CHARACTER(LEN=28)                                 :: YFILEDIA
REAL,    ALLOCATABLE              , DIMENSION(:,:):: ZWORKDATIME   ! global time
                                                                   !     info
INTEGER                                           :: JSV           ! loop index
                                                                   ! over the 
                                                                   ! KSV  SVx
INTEGER :: IP
TYPE(TFIELDDATA) :: TZFIELD
!
!-------------------------------------------------------------------------------
!
YFILEDIA = TPDIAFILE%CNAME
!
!*	 1.     write TSTEP and BULEN
!	        ---------------------
!
TZFIELD%CMNHNAME   = 'TSTEP'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = 'TSTEP'
TZFIELD%CUNITS     = 's'
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = 'Time step'
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,PTSTEP)
!
TZFIELD%CMNHNAME   = 'BULEN'
TZFIELD%CSTDNAME   = ''
TZFIELD%CLONGNAME  = 'BULEN'
TZFIELD%CUNITS     = 's'
TZFIELD%CDIR       = '--'
TZFIELD%CCOMMENT   = 'Time step'
TZFIELD%NGRID      = 0
TZFIELD%NTYPE      = TYPEREAL
TZFIELD%NDIMS      = 0
TZFIELD%LTIMEDEP   = .FALSE.
CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,XBULEN)
!
!*   1.1   initialize NBUTSHIFT
!           ---------------------
!
NBUTSHIFT = NBUTSHIFT+1
!
!
SELECT CASE (CBUTYPE)
!
!-------------------------------------------------------------------------------
!
!*	 2.     'CART' CASE
!	        -----------
!
  CASE('CART','SKIP')
    GNOCOMPRESS=(.NOT.LBU_ICP .AND. .NOT.LBU_JCP)
!
!*	 2.1    Initialization
!
    ALLOCATE(ZWORKTEMP(1,1))
    ALLOCATE(ZWORKDATIME(16,1))
!
    ZWORKDATIME(1,1)=TDTEXP%TDATE%YEAR
    ZWORKDATIME(2,1)=TDTEXP%TDATE%MONTH
    ZWORKDATIME(3,1)=TDTEXP%TDATE%DAY
    ZWORKDATIME(4,1)=TDTEXP%TIME
    ZWORKDATIME(5,1)=TDTSEG%TDATE%YEAR
    ZWORKDATIME(6,1)=TDTSEG%TDATE%MONTH
    ZWORKDATIME(7,1)=TDTSEG%TDATE%DAY
    ZWORKDATIME(8,1)=TDTSEG%TIME
    ZWORKDATIME(9,1)=TPDTMOD%TDATE%YEAR
    ZWORKDATIME(10,1)=TPDTMOD%TDATE%MONTH
    ZWORKDATIME(11,1)=TPDTMOD%TDATE%DAY
    ZWORKDATIME(12,1)=TPDTMOD%TIME
!
    CALL DATETIME_DISTANCE(TDTEXP,TPDTCUR,ZWORKTEMP(1,1))
!
    ZWORKTEMP(1,1)=ZWORKTEMP(1,1)+(1.-NBUSTEP*0.5)*PTSTEP
!
    ZWORKDATIME(13,1)=TDTEXP%TDATE%YEAR
    ZWORKDATIME(14,1)=TDTEXP%TDATE%MONTH
    ZWORKDATIME(15,1)=TDTEXP%TDATE%DAY
    ZWORKDATIME(16,1)=TDTEXP%TIME+ZWORKTEMP(1,1)
!
!*	 2.2    storage of the budgets array
!
!*	 2.2.1  RU budget
!
    IF (LBU_RU) THEN   
!                         XBURHODJU and RU  budgets 
!
      IP=1
! unit conversion for  RU budgets 
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RU
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURU(:,:,:,JPROC) * ZCONVERT(JPROC)    &
                                                    / XBURHODJU(:,:,:)
        END DO
      ELSE
        ALLOCATE(ZWORK(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,1)) ! global budget of RhodjU
        ZWORK(:,:,:,1,1,1)=END_CART_COMPRESS(XBURHODJU(:,:,:))
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RU
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURU(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)   &
                                                    / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
!
!  RU budgets storage
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 'm s-2'; YWORKUNIT(1:3) = 'm s-1'
      YWORKCOMMENT(:)    = 'Budget of momentum along X axis'
      IWORKGRID(:)       = 2

      WRITE(YGROUP_NAME,FMT="('UU___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID,  &
                         ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(1, :),  &
                         YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
                         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
! XBURHODJU storage
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORK(NBUIMAX,NBUJMAX,NBUKMAX,1,1,1)) ! local budget of RHODJU
        ZWORK(:,:,:,1,1,1) = XBURHODJU(:,:,:)
      END IF
      ALLOCATE(YBUCOMMENT(1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
      YBUCOMMENT(1)      = 'RhodJX'
      YWORKUNIT(1)       = 'kg'
      YWORKCOMMENT(1)    = 'RhodJ for momentum along X axis'
      IWORKGRID(1)       = 2
      WRITE(YGROUP_NAME,FMT="('RJX__',I4.4)") NBUTSHIFT
! 
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME,'CART', IWORKGRID,   &
                         ZWORKDATIME, ZWORK, ZWORKTEMP, YBUCOMMENT,         &
                         YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
                         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORK, YBUCOMMENT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 2.2.2  RV budget
!
    IF (LBU_RV) THEN 
                       ! XBURHODJV and RV budgets
!
      IP=2
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RV
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURV(:,:,:,JPROC) * ZCONVERT(JPROC)    &
                                                    / XBURHODJV(:,:,:)
        END DO
      ELSE
        ALLOCATE(ZWORK(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,1)) ! global budget of RhodjV
        ZWORK(:,:,:,1,1,1)=END_CART_COMPRESS(XBURHODJV(:,:,:))
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RV
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURV(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
                                                    / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
!
!  RV budgets storage
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
                                !
      YWORKUNIT(:)       = 'm s-2'; YWORKUNIT(1:3) = 'm s-1'
      YWORKCOMMENT(:)    = 'Budget of momentum along Y axis'
      IWORKGRID(:)       = 3
      WRITE(YGROUP_NAME,FMT="('VV___',I4.4)") NBUTSHIFT
                                !
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!                     XBURHODJV storage
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORK(NBUIMAX,NBUJMAX,NBUKMAX,1,1,1)) ! local budget of RHODJV
        ZWORK(:,:,:,1,1,1) = XBURHODJV(:,:,:)
      END IF
      ALLOCATE(YBUCOMMENT(1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
      YBUCOMMENT(1)      = 'RhodJY'
      YWORKUNIT(1)       = 'kg'
      YWORKCOMMENT(1)    = 'RhodJ for momentum along Y axis'
      IWORKGRID(1)       = 3
      WRITE(YGROUP_NAME,FMT="('RJY__',I4.4)") NBUTSHIFT
! 
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME,'CART', IWORKGRID, &
           ZWORKDATIME, ZWORK, ZWORKTEMP, YBUCOMMENT,         &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(YBUCOMMENT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
      DEALLOCATE(ZWORK)
    END IF
!
!
!*	 2.2.3  RW budget
!
    IF (LBU_RW) THEN
! 
      IP=3
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RW
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURW(:,:,:,JPROC) * ZCONVERT(JPROC)    &
               / XBURHODJW(:,:,:)
        END DO
      ELSE
        ALLOCATE(ZWORK(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,1)) ! global budget of RhodjW
        ZWORK(:,:,:,1,1,1)=END_CART_COMPRESS(XBURHODJW(:,:,:))
                                !
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RW
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURW(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
               / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
!
!  RW budgets storage
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 'm s-2'; YWORKUNIT(1:3) = 'm s-1'
      YWORKCOMMENT(:)    = 'Budget of momentum along Z axis'
      IWORKGRID(:)       = 4
      WRITE(YGROUP_NAME,FMT="('WW___',I4.4)") NBUTSHIFT
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!                     XBURHODJW storage
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORK(NBUIMAX,NBUJMAX,NBUKMAX,1,1,1)) ! local budget of RHODJW
        ZWORK(:,:,:,1,1,1) = XBURHODJW(:,:,:)
      END IF
      ALLOCATE(YBUCOMMENT(1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
      YBUCOMMENT(1)      = 'RhodJZ'
      YWORKUNIT(1)       = 'kg'
      YWORKCOMMENT(1)    = 'RhodJ for momentum along Z axis'
      IWORKGRID(1)       = 4
      WRITE(YGROUP_NAME,FMT="('RJZ__',I4.4)") NBUTSHIFT
! 
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME,'CART', IWORKGRID, &
           ZWORKDATIME, ZWORK, ZWORKTEMP, YBUCOMMENT,         &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(YBUCOMMENT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
      DEALLOCATE(ZWORK)
    END IF
  
!
!*	 2.2.3'  XBURHODJ storage for Scalars
!
    IF (LBU_RTH .OR. LBU_RTKE .OR. LBU_RRV .OR. LBU_RRC .OR. LBU_RRR .OR. &
        LBU_RRI .OR. LBU_RRS  .OR. LBU_RRG .OR. LBU_RRH .OR. LBU_RSV      ) THEN
!      
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORK(NBUIMAX,NBUJMAX,NBUKMAX,1,1,1)) ! local budget of RHODJ
        ZWORK(:,:,:,1,1,1) = XBURHODJ(:,:,:)
      ELSE
        ALLOCATE(ZWORK(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,1)) ! global budget of RodhjW
        ZWORK(:,:,:,1,1,1)=END_CART_COMPRESS(XBURHODJ(:,:,:))
      END IF
      ALLOCATE(YBUCOMMENT(1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
      YBUCOMMENT(1)      = 'RhodJS'
      YWORKUNIT(1)       = 'kg'
      YWORKCOMMENT(1)    = 'RhodJ for Scalars variables'
      IWORKGRID(1)       = 1
      WRITE(YGROUP_NAME,FMT="('RJS__',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORK, ZWORKTEMP, YBUCOMMENT,         &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      IF (GNOCOMPRESS) THEN
        DEALLOCATE(ZWORK, YBUCOMMENT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
      ELSE
        DEALLOCATE(YBUCOMMENT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
      ENDIF
!
    ENDIF
!
!*	 2.2.4  RTH budget
!
    IF (LBU_RTH) THEN
!  RTH budgets storage
      IP=4
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RTH
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURTH(:,:,:,JPROC) * ZCONVERT(JPROC)    &
               / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RTH
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURTH(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)   &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 'K s-1' ; YWORKUNIT(1:3) = 'K'
      YWORKCOMMENT(:)    = 'Budget of potential temperature'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('TH___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.5  RTKE budget
!
    IF (LBU_RTKE) THEN
!  RTKE budgets storage
      IP=5
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RTKE
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURTKE(:,:,:,JPROC) * ZCONVERT(JPROC)    &
             / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RTKE
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURTKE(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 'm2 s-3' ; YWORKUNIT(1:3) = 'm2 s-1'
      YWORKCOMMENT(:)    = 'Budget of turbulent kinetic energy'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('TK___',I4.4)") NBUTSHIFT 
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF   
!
!*	 2.2.6  RRV budget
!
    IF (LBU_RRV) THEN
!  RRV budgets storage
      IP=6
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RTKE
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURRV(:,:,:,JPROC) * ZCONVERT(JPROC)    &
             / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RTKE
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURRV(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of water vapor mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RV___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.7  RRC budget
!
    IF (LBU_RRC) THEN
!  RRV budgets storage
      IP=7
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RRC
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURRC(:,:,:,JPROC) * ZCONVERT(JPROC)    &
             / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RRC
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURRC(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of cloud water mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RC___',I4.4)") NBUTSHIFT 
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.8  RRR budget
!
    IF (LBU_RRR) THEN
      IP=8
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RRR
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURRR(:,:,:,JPROC) * ZCONVERT(JPROC)    &
             / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RRR
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURRR(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of rain water mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RR___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.9  RRI budget
!
    IF (LBU_RRI) THEN
      IP=9
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RRI
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURRI(:,:,:,JPROC) * ZCONVERT(JPROC)    &
               / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RRI
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURRI(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of cloud ice mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RI___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.10  RRS budget
!
    IF (LBU_RRS) THEN
      IP=10
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RRS
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURRS(:,:,:,JPROC) * ZCONVERT(JPROC)    &
             / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RRS
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURRS(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of snow/aggregate mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RS___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.11  RRG budget
!
    IF (LBU_RRG) THEN
      IP=11
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RRG
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURRG(:,:,:,JPROC) * ZCONVERT(JPROC)    &
             / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RRG
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURRG(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of graupel mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RG___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.12  RRH budget
!
    IF (LBU_RRH) THEN
      IP=12
      ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
      IF (GNOCOMPRESS) THEN
        ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RRH
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = XBURRH(:,:,:,JPROC) * ZCONVERT(JPROC)    &
             / XBURHODJ(:,:,:)
        END DO
      ELSE
!
        ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RRH
!      
        DO JPROC=1,NBUPROCNBR(IP)
          ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURRH(:,:,:,JPROC))
          ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)  &
             / ZWORK(:,:,:,1,1,1)
        END DO
      ENDIF
      DEALLOCATE(ZCONVERT)
! 
      ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
      YWORKUNIT(:)       = 's-1' ; YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of hail mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RH___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
           ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
           YWORKUNIT, YWORKCOMMENT,                           &
           LBU_ICP, LBU_JCP, LBU_KCP,                         & 
           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
    END IF
!
!*	 2.2.13  RSV budget
!
    IF (LBU_RSV) THEN
      DO JSV = 1,KSV
        IP=12+JSV
        ALLOCATE(ZCONVERT(NBUPROCNBR(IP)))
        ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
        ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
        ZCONVERT(4:NBUPROCNBR(IP)) = 1.
! 
        IF (GNOCOMPRESS) THEN
          ALLOCATE(ZWORKT(NBUIMAX,NBUJMAX,NBUKMAX,1,1,NBUPROCNBR(IP))) ! local budget of RRH
          DO JPROC=1,NBUPROCNBR(IP)
            ZWORKT(:,:,:,1,1,JPROC) = XBURSV(:,:,:,JPROC,JSV) * ZCONVERT(JPROC)  &
               / XBURHODJ(:,:,:)
          END DO
        ELSE
!
          ALLOCATE(ZWORKT(NBUIMAX_ll,NBUJMAX_ll,NBUKMAX,1,1,NBUPROCNBR(IP))) ! global budget of RRH
!      
          DO JPROC=1,NBUPROCNBR(IP)
            ZWORKT(:,:,:,1,1,JPROC) = END_CART_COMPRESS(XBURSV(:,:,:,JPROC,JSV))
            ZWORKT(:,:,:,1,1,JPROC) = ZWORKT(:,:,:,1,1,JPROC)* ZCONVERT(JPROC)    &
               / ZWORK(:,:,:,1,1,1)
          END DO
         DEALLOCATE(ZWORK)
        ENDIF
        DEALLOCATE(ZCONVERT)
!   
        ALLOCATE(YWORKUNIT(NBUPROCNBR(IP)))
        ALLOCATE(YWORKCOMMENT(NBUPROCNBR(IP)))
        ALLOCATE(IWORKGRID(NBUPROCNBR(IP)))
!
        YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = '  '
        DO JT = 1,NBUPROCNBR(IP)
          WRITE(YWORKCOMMENT(JT),FMT="('Budget of SVx=',I3.3)") JSV
        END DO
        IWORKGRID(:)       = 1
        WRITE(YGROUP_NAME,FMT="('SV',I3.3,I4.4)") JSV,NBUTSHIFT
!
        CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, 'CART', IWORKGRID, &
             ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(IP, :),   &
             YWORKUNIT, YWORKCOMMENT,                           &
             LBU_ICP, LBU_JCP, LBU_KCP,                         & 
             NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH           )
        DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
      END DO
    END IF
!
    IF (ALLOCATED(ZWORK)) DEALLOCATE(ZWORK)
    DEALLOCATE (ZWORKTEMP, ZWORKDATIME)
!-------------------------------------------------------------------------------
!
!*	 3.     'MASK' CASE
!	        -----------
!
  CASE('MASK')
    ALLOCATE(ZWORKTEMP(NBUWRNB,1))
    ALLOCATE(ZWORKDATIME(16,NBUWRNB))
    ALLOCATE(ZWORKMASK(SIZE(XBUSURF,1),SIZE(XBUSURF,2),1,NBUWRNB,NBUMASK,1))
!
! local array
    DO JMASK=1,NBUMASK
      DO JT=1,NBUWRNB
        ZWORKMASK(:,:,1,JT,JMASK,1) = XBUSURF(:,:,JMASK,JT)
      END DO
    END DO
!
    ZWORKDATIME(1,:)=TDTEXP%TDATE%YEAR
    ZWORKDATIME(2,:)=TDTEXP%TDATE%MONTH
    ZWORKDATIME(3,:)=TDTEXP%TDATE%DAY
    ZWORKDATIME(4,:)=TDTEXP%TIME
    ZWORKDATIME(5,:)=TDTSEG%TDATE%YEAR
    ZWORKDATIME(6,:)=TDTSEG%TDATE%MONTH
    ZWORKDATIME(7,:)=TDTSEG%TDATE%DAY
    ZWORKDATIME(8,:)=TDTSEG%TIME
    ZWORKDATIME(9,:)=TPDTMOD%TDATE%YEAR
    ZWORKDATIME(10,:)=TPDTMOD%TDATE%MONTH
    ZWORKDATIME(11,:)=TPDTMOD%TDATE%DAY
    ZWORKDATIME(12,:)=TPDTMOD%TIME
!
    CALL DATETIME_DISTANCE(TDTEXP,TPDTCUR,ZWORKTEMP(NBUWRNB,1))
!
    ZWORKTEMP(NBUWRNB,1)=ZWORKTEMP(NBUWRNB,1)+(1.-NBUSTEP*0.5)*PTSTEP
!
    ZWORKDATIME(13,NBUWRNB)=TDTEXP%TDATE%YEAR
    ZWORKDATIME(14,NBUWRNB)=TDTEXP%TDATE%MONTH
    ZWORKDATIME(15,NBUWRNB)=TDTEXP%TDATE%DAY
    ZWORKDATIME(16,NBUWRNB)=TDTEXP%TIME+ZWORKTEMP(NBUWRNB,1)
    DO JT=1,NBUWRNB-1
      ZWORKTEMP(JT,1) = ZWORKTEMP(NBUWRNB,1)-NBUSTEP*PTSTEP*(NBUWRNB-JT)
      ZWORKDATIME(13,JT)=TDTEXP%TDATE%YEAR
      ZWORKDATIME(14,JT)=TDTEXP%TDATE%MONTH
      ZWORKDATIME(15,JT)=TDTEXP%TDATE%DAY
      ZWORKDATIME(16,JT)=TDTEXP%TIME  + ZWORKTEMP(JT,1)
    END DO
!
!*     3.1    storage of the masks  array
!
        WRITE(TZFIELD%CMNHNAME,FMT="('MASK_',I4.4,'.MASK')") NBUTSHIFT
        TZFIELD%CSTDNAME   = ''
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CUNITS     = ''
        TZFIELD%CDIR       = 'XY'
        WRITE(TZFIELD%CCOMMENT,FMT="('X_Y_MASK',I4.4)") NBUTSHIFT
        TZFIELD%NGRID      = 1
        TZFIELD%NTYPE      = TYPEREAL
        TZFIELD%NDIMS      = 6
        TZFIELD%LTIMEDEP   = .FALSE.
        CALL IO_WRITE_FIELD(TPDIAFILE,TZFIELD,ZWORKMASK(:,:,:,:,:,:))
        WRITE(YRECFM,FMT="('MASK_',I4.4)") NBUTSHIFT
        CALL MENU_DIACHRO(TPDIAFILE,YRECFM)
        DEALLOCATE(ZWORKMASK)
!
!*	 3.2    storage of the budgets array
!
!*	 3.2.1  RU budget
!
    IF (LBU_RU) THEN
                       ! XBURHODJU storage
!
      ALLOCATE(ZWORK(1,1,NBUKMAX,NBUWRNB,NBUMASK,1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
      ZWORK(1,1,:,:,:,1) = END_MASK_COMPRESS(XBURHODJU(:,:,:))
      WHERE  (ZWORK(1,1,:,:,:,1) <= 0.)
          ZWORK(1,1,:,:,:,1)=-999.
      END WHERE
      YWORKUNIT(:)       = 'kg'
      YWORKCOMMENT(:)    = 'RhodJ for momentum along X axis'
      IWORKGRID(:)       = 2
      WRITE(YGROUP_NAME,FMT="('RJX__',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                         ZWORKDATIME, ZWORK, ZWORKTEMP, CBUCOMMENT(1, :),   &
                         YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
                         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE( YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
!                 unit conversion of RU budgets and storage
!                 -----------------------------------------
!
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(1)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(1)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(1)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(1)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(1)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       =  PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(1)) = 1.
      DO JPROC=1,NBUPROCNBR(1)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURU(:,:,:,JPROC)) &
                               * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT, ZWORK)
      
!
      YWORKUNIT(:)       = 'm s-2'; YWORKUNIT(1:3) = 'm s-1'
      YWORKCOMMENT(:)    = 'Budget of momentum along X axis'
      IWORKGRID(:)       = 2
      WRITE(YGROUP_NAME,FMT="('UU___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                         ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(1, :),   &
                         YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.2  RV budget
!
    IF (LBU_RV) THEN
                       ! XBURHODJV storage
!
      ALLOCATE(ZWORK(1,1,NBUKMAX,NBUWRNB,NBUMASK,1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
      ZWORK(1,1,:,:,:,1)= END_MASK_COMPRESS( XBURHODJV(:,:,:))
      WHERE ( ZWORK(1,1,:,:,:,1) <= 0.)
        ZWORK(1,1,:,:,:,1)=-999.
      END WHERE
      YWORKUNIT(:)       = 'kg'
      YWORKCOMMENT(:)    = 'RhodJ for momentum along Y axis'
      IWORKGRID(:)       = 3
      WRITE(YGROUP_NAME,FMT="('RJY__',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORK, ZWORKTEMP, CBUCOMMENT(1, :),   &
	       		 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE( YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
!                 unit conversion of RU budgets and storage
!
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(2)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(2)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(2)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(2)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(2)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(2)) = 1.
      DO JPROC=1,NBUPROCNBR(2)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURV  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1) 
      END DO
      DEALLOCATE(ZCONVERT, ZWORK)
!
      YWORKUNIT(:)       = 'm s-2'; YWORKUNIT(1:3) = 'm s-1'
      YWORKCOMMENT(:)    = 'Budget of momentum along Y axis'
      IWORKGRID(:)       = 3
      WRITE(YGROUP_NAME,FMT="('VV___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                         ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(2, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.3  RW budget
!
    IF (LBU_RW) THEN
                       ! XBURHODJW storage
!
      ALLOCATE(ZWORK(1,1,NBUKMAX,NBUWRNB,NBUMASK,1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
      ZWORK(1,1,:,:,:,1)=END_MASK_COMPRESS(XBURHODJW(:,:,:))
      WHERE (ZWORK(1,1,:,:,:,1) <= 0.)
        ZWORK(1,1,:,:,:,1)=-999.
      END WHERE
      YWORKUNIT(:)       = 'kg'
      YWORKCOMMENT(:)    = 'RhodJ for momentum along Z axis'
      IWORKGRID(:)       = 4
      WRITE(YGROUP_NAME,FMT="('RJZ__',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORK, ZWORKTEMP, CBUCOMMENT(1, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE( YWORKUNIT, YWORKCOMMENT, IWORKGRID)
!
!                 unit conversion of RU budgets and storage
!
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(3)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(3)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(3)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(3)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(3)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(3)) = 1.
      DO JPROC=1,NBUPROCNBR(3)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURW (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT, ZWORK)
!
      YWORKUNIT(:)       = 'm s-2'; YWORKUNIT(1:3) = 'm s-1'
      YWORKCOMMENT(:)    = 'Budget of momentum along Z axis'
      IWORKGRID(:)       = 4
      WRITE(YGROUP_NAME,FMT="('WW___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(3, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.3'  XBURHODJ storage for Scalars
!
    IF (LBU_RTH .OR. LBU_RTKE .OR. LBU_RRV .OR. LBU_RRC .OR. LBU_RRR .OR. &
        LBU_RRI .OR. LBU_RRS  .OR. LBU_RRG .OR. LBU_RRH .OR. LBU_RSV      ) THEN
!
      ALLOCATE(ZWORK(1,1,NBUKMAX,NBUWRNB,NBUMASK,1))
      ALLOCATE(YWORKUNIT(1))
      ALLOCATE(YWORKCOMMENT(1))
      ALLOCATE(IWORKGRID(1))
!
        ZWORK(1,1,:,:,:,1) = END_MASK_COMPRESS(XBURHODJ(:,:,:))
        WHERE (ZWORK(1,1,:,:,:,1) <= 0.)
         ZWORK(1,1,:,:,:,1)=-999.
        END WHERE
      YWORKUNIT(:)       = 'kg'
      YWORKCOMMENT(:)    = 'RhodJ for Scalars'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RJS__',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORK, ZWORKTEMP, CBUCOMMENT(1, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE( YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.4  RTH budget
!
    IF (LBU_RTH) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(4)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(4)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(4)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(4)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(4)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(4)) = 1.
      DO JPROC=1,NBUPROCNBR(4)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURTH  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 'K s-1' ; YWORKUNIT(1:3) = 'K'
      YWORKCOMMENT(:)    = 'Budget of potential temperature'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('TH___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(4, :),   &
	         	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.5  RTKE budget
!
    IF (LBU_RTKE) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(5)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(5)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(5)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(5)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(5)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(5)) = 1.
      DO JPROC=1,NBUPROCNBR(5)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURTKE (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 'm2 s-3' ; YWORKUNIT(1:3) = 'm2 s-2'
      YWORKCOMMENT(:)    = 'Budget of turbulent kinetic energy'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('TK___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(5, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF   
!
!*	 3.2.6  RRV budget
!
    IF (LBU_RRV) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(6)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(6)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(6)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(6)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(6)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(6)) = 1.
      DO JPROC=1,NBUPROCNBR(6)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURRV  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of water vapor mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RV___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(6, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.7  RRC budget
!
    IF (LBU_RRC) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(7)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(7)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(7)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(7)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(7)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(7)) = 1.
      DO JPROC=1,NBUPROCNBR(7)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURRC  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of cloud water mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RC___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(7, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.8  RRR budget
!
    IF (LBU_RRR) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(8)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(8)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(8)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(8)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(8)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(8)) = 1.
      DO JPROC=1,NBUPROCNBR(8)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURRR  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of rain water mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RR___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(8, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.9  RRI budget
!
    IF (LBU_RRI) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(9)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(9)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(9)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(9)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(9)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(9)) = 1.
      DO JPROC=1,NBUPROCNBR(9)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURRI  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of cloud ice mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RI___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(9, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.10  RRS budget
!
    IF (LBU_RRS) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(10)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(10)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(10)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(10)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(10)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(10)) = 1.
      DO JPROC=1,NBUPROCNBR(10)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURRS  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of snow/aggregate mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RS___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(10, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.11  RRG budget
!
    IF (LBU_RRG) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(11)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(11)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(11)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(11)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(11)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(11)) = 1.
      DO JPROC=1,NBUPROCNBR(11)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURRG  (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT )
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of graupel mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RG___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(11, :),   &
	         	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.12  RRH budget
!
    IF (LBU_RRH) THEN
      ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(12)))
      ALLOCATE(YWORKUNIT(NBUPROCNBR(12)))
      ALLOCATE(YWORKCOMMENT(NBUPROCNBR(12)))
      ALLOCATE(IWORKGRID(NBUPROCNBR(12)))
!
      ALLOCATE(ZCONVERT(NBUPROCNBR(12)))
      ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
      ZCONVERT(4:NBUPROCNBR(12)) = 1.
      DO JPROC=1,NBUPROCNBR(12)
        ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURRH (:,:,:,JPROC)) &
                         * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
      END DO
      DEALLOCATE(ZCONVERT)
!
      YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = 'kg kg-1'
      YWORKCOMMENT(:)    = 'Budget of hail mixing ratio'
      IWORKGRID(:)       = 1
      WRITE(YGROUP_NAME,FMT="('RH___',I4.4)") NBUTSHIFT
!
      CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID, &
                	 ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(12, :),   &
	        	 YWORKUNIT, YWORKCOMMENT,                           &
                         LBU_ICP, LBU_JCP, LBU_KCP,                         & 
		         NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
      DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
    END IF
!
!*	 3.2.13  RSV budget
!
    IF (LBU_RSV) THEN
      DO JSV = 1,KSV
        ALLOCATE(ZWORKT(1,1,NBUKMAX,NBUWRNB,NBUMASK,NBUPROCNBR(12+JSV)))
        ALLOCATE(YWORKUNIT(NBUPROCNBR(12+JSV)))
        ALLOCATE(YWORKCOMMENT(NBUPROCNBR(12+JSV)))
        ALLOCATE(IWORKGRID(NBUPROCNBR(12+JSV)))
!
        ALLOCATE(ZCONVERT(NBUPROCNBR(12+JSV)))
        ZCONVERT(1:2)     = PTSTEP * REAL(NBUSTEP)
        ZCONVERT(3)       = PTSTEP * REAL(NBUSTEP)
        ZCONVERT(4:NBUPROCNBR(12+JSV)) = 1.
        DO JPROC=1,NBUPROCNBR(12+JSV)
          ZWORKT(1,1,:,:,:,JPROC) = END_MASK_COMPRESS( XBURSV  (:,:,:,JPROC,JSV)) &
                           * ZCONVERT(JPROC) / ZWORK(1,1,:,:,:,1)
        END DO
        DEALLOCATE(ZCONVERT)
!
        YWORKUNIT(:)       = 's-1' ;  YWORKUNIT(1:3) = '  '
        DO JT = 1,NBUPROCNBR(12+JSV)
          WRITE(YWORKCOMMENT(JT),FMT="('Budget of SVx=',I3.3)") JSV
        END DO
        IWORKGRID(:)       = 1
        WRITE(YGROUP_NAME,FMT="('SV',I3.3,I4.4)") JSV,NBUTSHIFT
!
        CALL WRITE_DIACHRO(TPDIAFILE, TLUOUT, YGROUP_NAME, CBUTYPE, IWORKGRID,  &
                  	   ZWORKDATIME, ZWORKT, ZWORKTEMP, CBUCOMMENT(12+JSV,:),&
		           YWORKUNIT, YWORKCOMMENT,                            &
                           LBU_ICP, LBU_JCP, LBU_KCP,                          &
		           NBUIL, NBUIH, NBUJL, NBUJH, NBUKL, NBUKH )
        DEALLOCATE(ZWORKT, YWORKUNIT, YWORKCOMMENT, IWORKGRID)
      END DO
    END IF
!
    IF (LBU_RTH .OR. LBU_RTKE .OR. LBU_RRV .OR. LBU_RRC .OR. LBU_RRR .OR. &
        LBU_RRI .OR. LBU_RRS  .OR. LBU_RRG .OR. LBU_RRH .OR. LBU_RSV      ) THEN
      DEALLOCATE(ZWORK)
    END IF
!
  DEALLOCATE (ZWORKTEMP, ZWORKDATIME)
!
END SELECT   
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_BUDGET
