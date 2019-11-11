!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ########################
      SUBROUTINE INI_SERIES_n
!     ########################
!
!!****  *INI_SERIES_n* - routine to Initialize and Allocate arrays of diagnostics
!     for diachro files (temporal series)  for model _n
!!                           
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!       NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq from contributions of M. Tomasini, S. Donier,... * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    4/03/2002
!!                Oct. 2011 : (P.Le Moigne) Surface series
!!      June 2016: P. Wautelet: corrected writes
!!      Nov. 2017: J.-P. Chaboureau: fix a bug in dimension check
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*    0. Declaration
!     --------------
!
USE MODE_ll
USE MODE_IO_ll
USE MODE_MSG
USE MODE_MODELN_HANDLER
!
USE MODD_TIME        ! Experiment and segment times (TDTEXP and TDTSEG)
USE MODD_CONF
USE MODD_DYN, ONLY: XSEGLEN
USE MODD_SERIES
USE MODD_SERIES_n
USE MODD_PARAMETERS
USE MODD_CONF_n, ONLY: LUSERV,LUSERC,LUSERR,LUSERI,LUSERS,LUSERG,LUSERH
USE MODD_DIM_n, ONLY: NKMAX
USE MODD_DYN_n, ONLY: XTSTEP,NSTOP
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAM_n, ONLY: CSURF
USE MODD_PRECIP_n, ONLY: XINPRR,XINPRS,XINPRG
USE MODD_TIME_n
! 
USE MODI_MNHGET_SURF_PARAM_n
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!*       0.2     Local variables
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSEA !sea/ocean fraction
LOGICAL :: GMASKLANDSEA   ! local for LMASKLANDSEA
INTEGER :: IIMAX_ll  ! total physical domain I size
INTEGER :: IJMAX_ll  ! total physical domain J size
INTEGER :: IIE_ll,IJE_ll,IIB_ll,IJB_ll  ! begin and end of total physical domain
INTEGER :: IIB,IJB         ! Begin of physical dimensions of the sub-domain
INTEGER :: IIE,IJE,IKE     ! End   of physical dimensions of the sub-domain
INTEGER :: IIU,IJU ! begin and end of the sub-domain
INTEGER :: IINFO_ll ! return status of the interface routine
INTEGER :: IIDIM1 ! I size of the slice
INTEGER :: JJ,JI ! loop indices
INTEGER :: ISB1,ISB2,ISB3
INTEGER :: ISER
CHARACTER (LEN=5), DIMENSION(3) :: YSUF
INTEGER  :: ILUOUT ! Logical unit number for output-listing
INTEGER  :: IRESP   ! Return code of FM-routines
INTEGER :: IMI ! Current model index
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATIONS
!              ---------------
!
IMI = GET_CURRENT_MODEL_INDEX()
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE) 
IKE = NKMAX + JPVEXT
CALL GET_DIM_EXT_ll ('B', IIU,IJU)
CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
IIB_ll = 1 ! + JPHEXT
IJB_ll = 1 ! + JPHEXT
IIE_ll = IIMAX_ll ! + JPHEXT
IJE_ll = IJMAX_ll ! + JPHEXT
!
ILUOUT = TLUOUT%NLU
!
!
!*       1.1   Dimensions
!
IF ( ( NFREQSERIES*XTSTEP < XSEGLEN )       .AND.                        &
     ( ( NIBOXL      < IIB_ll         )      .OR.                        &
       ( NIBOXL      > NIBOXH         )      .OR.                        &
       ( NIBOXH      > IIE_ll         )      .OR.                        &
       ( NJBOXL      < IJB_ll         )      .OR.                        &
       ( NJBOXL      > NJBOXH         )      .OR.                        &
       ( NJBOXH      > IJE_ll         )      .OR.                        &
       ( NBJSLICE   >  IJE_ll         )      .OR.                        &
       ( NBJSLICE   > SIZE(NJSLICEL(:),1)  ) .OR.                        &
       ( NJSLICEL(1) < IJB_ll         )      .OR.                        &
       ( NJSLICEL(1) > IJE_ll         )      .OR.                        &
       ( NJSLICEL(1) > NJSLICEH(1)    )      .OR.                        &
       ( NJSLICEH(1) < IJB_ll         )      .OR.                        &
       ( NJSLICEH(1) > IJE_ll         )      .OR.                        &
       ( NKCLS       > IKE            )      .OR.                        &
       ( NKCLA       > IKE            )      .OR.                        &
       ( NKLOW       > IKE            )      .OR.                        &
       ( NKMID       > IKE            )      .OR.                        &
       ( NKUP        > IKE            ) )                           ) THEN

   WRITE(ILUOUT,FMT=*) '**********************************************'
   WRITE(ILUOUT,FMT=*) 'STOP in INI_SERIESn due to incompatible dimension:'
   WRITE(ILUOUT,FMT=*) ' IIB,IJB,IIE,IJE,IKE= ',IIB_ll,IJB_ll,IIE_ll,IJE_ll,IKE
   WRITE(ILUOUT,FMT=*) ' NIBOXL,NIBOXH,NJBOXL,NJBOXH= ',                 &
                        NIBOXL,NIBOXH,NJBOXL,NJBOXH
   WRITE(ILUOUT,FMT=*) ' NBJSLICE,NJSLICEL(1),NJSLICEH(1)= ',            &
                        NBJSLICE,NJSLICEL(1),NJSLICEH(1)
   WRITE(ILUOUT,FMT=*) ' NKCLS,NKCLA,NKLOW,NKMID,NKUP= ',                &
                        NKCLS,NKCLA,NKLOW,NKMID,NKUP
   WRITE(ILUOUT,FMT=*) '**********************************************'
!callabortstop
   CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SERIES_n','')
END IF
!
ALLOCATE(LINBOX(IIU,IJU))
LINBOX(:,:)=.FALSE.
CALL GET_INTERSECTION_ll(NIBOXL+JPHEXT,NJBOXL+JPHEXT,NIBOXH+JPHEXT,NJBOXH+JPHEXT, &
     NIBOXSL,NJBOXSL,NIBOXSH,NJBOXSH,"PHYS",IINFO_ll)
IF  ( IINFO_ll /= 1 ) THEN !
  DO JI=IIB,IIE
    DO JJ=IJB,IJE
     LINBOX(JI,JJ)=((JI>=NIBOXSL).AND.(JI<=NIBOXSH).AND. &
           (JJ>=NJBOXSL).AND.(JJ<=NJBOXSH))
    END DO
  END DO
  LSERIES1=.TRUE.
  LSERIES2=.TRUE.
  WRITE(UNIT=ILUOUT,FMT=*) 'INI_SERIESn: intersection with horizontal box'
ELSE ! the intersection is void
  LSERIES1=.FALSE.
  LSERIES2=.FALSE.
ENDIF
!
IIDIM1=0
IF (NBJSLICE > 0 ) THEN
  ALLOCATE(LSERIES3(NBJSLICE))
  DO JJ=1,NBJSLICE
    CALL GET_INTERSECTION_ll(IIB_ll+JPHEXT,NJSLICEL(JJ)+JPHEXT,IIE_ll+JPHEXT,NJSLICEH(JJ)+JPHEXT, &
    NISL(JJ),NJSLICESL(JJ),NISH(JJ),NJSLICESH(JJ),"PHYS",IINFO_ll)
    IF  ( IINFO_ll /= 1 ) THEN !
      LSERIES3(JJ)=.TRUE.
      IIDIM1=MAX(IIDIM1,NISH(JJ) -  NISL(JJ) + 1)
      IF (IIDIM1.EQ.0) THEN
        WRITE(UNIT=ILUOUT,FMT=*) 'STOP in INI_SERIESn: VOID INTERSECTION for slice ',JJ
        WRITE(ILUOUT,*) ' NJSLICEL=', NJSLICEL(JJ),'NJSLICEH=',NJSLICEH(JJ)
        WRITE(ILUOUT,*) ' NISL=',NISL(JJ),'NJSLICESL=',NJSLICESL(JJ),'NISH=',NISH(JJ),'NJSLICESH=',NJSLICESH(JJ)
!callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SERIES_n','')
      END IF
      WRITE(UNIT=ILUOUT,FMT=*) 'INI_SERIESn: intersection with slice ',JJ
    ELSE ! the intersection is void
      LSERIES3(JJ)=.FALSE.
      IIDIM1=MAX(IIDIM1,0)
    ENDIF
  END DO
END IF
!
LDOSERIES=(LSERIES1).OR.(LSERIES2)
DO JJ=1,NBJSLICE
  IF (LSERIES3(JJ)) LDOSERIES=.TRUE.
END DO
!
GMASKLANDSEA = LMASKLANDSEA
LMASKLANDSEA = GMASKLANDSEA .AND.(CSURF=='EXTE')
IF (LMASKLANDSEA .NEQV. GMASKLANDSEA) WRITE(UNIT=ILUOUT,FMT=*)  &
  'INI_SERIESn: CSURF=',CSURF,', LMASKLANDSEA set to ',LMASKLANDSEA
IF (LMASKLANDSEA) THEN
  ALLOCATE(LINBOXL(IIU,IJU))
  ALLOCATE(LINBOXS(IIU,IJU))
  LINBOXL(:,:)=.FALSE.
  LINBOXS(:,:)=.FALSE.
  ALLOCATE(ZSEA(IIU,IJU))
  CALL MNHGET_SURF_PARAM_n(PSEA=ZSEA)
  WHERE((ZSEA(:,:).LE.0.5).AND.LINBOX) 
    LINBOXL(:,:)=.TRUE.  ! LINBOX on Land zone
  END WHERE
  WHERE((ZSEA(:,:).GT.0.5).AND.LINBOX)
    LINBOXS(:,:)=.TRUE.  ! LINBOX on Sea zone
  END WHERE
  DEALLOCATE(ZSEA)
END IF
!
!*       1.2   Allocate arrays of diagnostics
!              ------------------------------
!
! Temporal series (t)
NSTEMP_SERIE1 = 0
! SURFACE FIELDS
IF (LSURF)           NSTEMP_SERIE1 = NSTEMP_SERIE1 +5
! end SURFACE FIELDS
IF (SIZE(XINPRR)/=0) NSTEMP_SERIE1 = NSTEMP_SERIE1 +2
IF (LUSERV)          NSTEMP_SERIE1 = NSTEMP_SERIE1 +1
IF (LUSERC)          NSTEMP_SERIE1 = NSTEMP_SERIE1 +1
IF (LUSERR)          NSTEMP_SERIE1 = NSTEMP_SERIE1 +1
IF (LUSERI)          NSTEMP_SERIE1 = NSTEMP_SERIE1 +1
IF (LUSERS)          NSTEMP_SERIE1 = NSTEMP_SERIE1 +1
IF (LUSERG)          NSTEMP_SERIE1 = NSTEMP_SERIE1 +1
IF (LUSERH)          NSTEMP_SERIE1 = NSTEMP_SERIE1 +1
IF (LMASKLANDSEA)    NSTEMP_SERIE1 = NSTEMP_SERIE1 *3
NAVER1 = NSTEMP_SERIE1
IF (LWMINMAX) THEN
  NSTEMP_SERIE1 = NSTEMP_SERIE1 +2
  IF (LMASKLANDSEA) NSTEMP_SERIE1 = NSTEMP_SERIE1 +4
END IF
!
! Temporal series (z,t)
NSTEMP_SERIE2 = 3             ! WT,THT,PABST
IF (LUSERV)          NSTEMP_SERIE2 = NSTEMP_SERIE2 +1
IF (LUSERC)          NSTEMP_SERIE2 = NSTEMP_SERIE2 +1
IF (LUSERR)          NSTEMP_SERIE2 = NSTEMP_SERIE2 +1
IF (LUSERI)          NSTEMP_SERIE2 = NSTEMP_SERIE2 +1
IF (LUSERS)          NSTEMP_SERIE2 = NSTEMP_SERIE2 +1
IF (LUSERG)          NSTEMP_SERIE2 = NSTEMP_SERIE2 +1
IF (LMASKLANDSEA)    NSTEMP_SERIE2 = NSTEMP_SERIE2 *3
!
! Temporal series (x,t)
NSTEMP_SERIE3 = 3            ! Ucls,Wcla,Wlow-up
IF (LUSERV)          NSTEMP_SERIE3 = NSTEMP_SERIE3 +2      ! RVcls,RVmid
IF (LUSERC)          NSTEMP_SERIE3 = NSTEMP_SERIE3 +1      ! RC0_up
IF (LUSERR)          NSTEMP_SERIE3 = NSTEMP_SERIE3 +1      ! RRcls
!
ALLOCATE( CSTITLE1   (NSTEMP_SERIE1) )
ALLOCATE( CSTITLE2   (NSTEMP_SERIE2) )
ALLOCATE( CSTITLE3   (NSTEMP_SERIE3) )
ALLOCATE( CSUNIT1    (NSTEMP_SERIE1) )
ALLOCATE( CSUNIT2    (NSTEMP_SERIE2) )
ALLOCATE( CSUNIT3    (NSTEMP_SERIE3) )
ALLOCATE( NSGRIDD1   (NSTEMP_SERIE1) )
ALLOCATE( NSGRIDD2   (NSTEMP_SERIE2) )
ALLOCATE( NSGRIDD3   (NSTEMP_SERIE3) )
!
NSTEMP3=NSTEMP_SERIE3 * NBJSLICE
!
NSNBSTEPT = NSTOP / NFREQSERIES
IF ( MOD(NSTOP,NFREQSERIES) /= 0 ) NSNBSTEPT = NSNBSTEPT + 1
!
ALLOCATE( XSSERIES1 ( 1 ,1,  1  ,NSNBSTEPT,1,NSTEMP_SERIE1) )
ALLOCATE( XSSERIES2 ( 1 ,1,NKMAX,NSNBSTEPT,1,NSTEMP_SERIE2) )
ALLOCATE( XSSERIES3 (IIDIM1,1,  1  ,NSNBSTEPT,1,NSTEMP3      ) )
!
ALLOCATE( CSCOMMENT1 (NSTEMP_SERIE1) )
ALLOCATE( CSCOMMENT2 (NSTEMP_SERIE2) )
ALLOCATE( CSCOMMENT3 (NSTEMP_SERIE3) )
!
ALLOCATE( XSTRAJT  (NSNBSTEPT,1)  )
ALLOCATE( XSDATIME (16,NSNBSTEPT) )
!
XSSERIES1(:,:,:,:,:,:)=0.
XSSERIES2(:,:,:,:,:,:)=0.
XSSERIES3(:,:,:,:,:,:)=0.
!
WRITE(ILUOUT,*) 'NSNBSTEPT=',NSNBSTEPT,' NSTOP=', NSTOP, &
                ' NFREQSERIES=',NFREQSERIES
WRITE(ILUOUT,*) 'NSTEMP_SERIE1=',NSTEMP_SERIE1,' NSTEMP_SERIE2=',NSTEMP_SERIE2
WRITE(ILUOUT,*) 'NSTEMP_SERIE3=',NSTEMP_SERIE3, ' NSTEMP3=',NSTEMP3
WRITE(ILUOUT,*) 'NIBOXL=',NIBOXL,' NJBOXL=',NJBOXL, &
                ' NIBOXH=',NIBOXH,' NJBOXH=',NJBOXH
WRITE(ILUOUT,*) 'NIBOXSL=',NIBOXSL,' NJBOXSL=',NJBOXSL, &
                ' NIBOXSH=',NIBOXSH,' NJBOXSH=',NJBOXSH
WRITE(ILUOUT,*) 'NJSLICEL=', NJSLICEL(1:NBJSLICE),' NJSLICEH=',NJSLICEH(1:NBJSLICE)
WRITE(ILUOUT,*) 'NISL=',NISL(1:NBJSLICE),' NJSLICESL=',NJSLICESL(1:NBJSLICE), &
                ' NISH=',NISH(1:NBJSLICE),' NJSLICESH=',NJSLICESH(1:NBJSLICE)
WRITE(ILUOUT,*) 'IIDIM1=',IIDIM1,' NBJSLICE=',NBJSLICE
!
!-------------------------------------------------------------------------------
!
!*    2. Variable units of the diagnostics arrays
!     -------------------------------------------
!
ISER=1
IF (LMASKLANDSEA) ISER=3
YSUF(1)='-GLOB'
YSUF(2)='-LAND'
YSUF(3)='-SEA '
!
!*       2.1   Temporal series t
!              -----------------
!
DO JI=1,NSTEMP_SERIE1
  WRITE(CSCOMMENT1(JI),'("TEMPORAL SERIE (t)   : ",A,".",I0,".",A)') CEXP,IMI,CSEG
  IF (JI==1) WRITE(ILUOUT,FMT=*) CSCOMMENT1(JI)
END DO
!
NSGRIDD1(:)=1
!
ISB1=0
DO JI=1,ISER
  ! total surface explicit precipitations
  IF (SIZE(XINPRR)/=0) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='INPRT'//YSUF(JI) ; CSUNIT1(ISB1)='MM/DAY'
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='ACPRT'//YSUF(JI) ; CSUNIT1(ISB1)='MM'
  END IF
  ! Mixing ratios
  IF (LUSERV) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='RVT'//YSUF(JI) ; CSUNIT1(ISB1)='KG/M2'
  END IF
  IF (LUSERC) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='RCT'//YSUF(JI) ; CSUNIT1(ISB1)='KG/M2'
  END IF
  IF (LUSERR) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='RRT'//YSUF(JI) ; CSUNIT1(ISB1)='KG/M2'
  ENDIF
  IF (LUSERI) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='RIT'//YSUF(JI) ; CSUNIT1(ISB1)='KG/M2'
  END IF
  IF (LUSERS) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='RST'//YSUF(JI) ; CSUNIT1(ISB1)='KG/M2'
  END IF
  IF (LUSERG) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='RGT'//YSUF(JI) ; CSUNIT1(ISB1)='KG/M2'
  END IF
  IF (LUSERH) THEN
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='RHT'//YSUF(JI) ; CSUNIT1(ISB1)='KG/M2'
  END IF
! SURFACE FIELDS
  IF (LSURF) THEN
     ISB1=ISB1+1 ; CSTITLE1(ISB1)='TS_WATER'//YSUF(JI) ; CSUNIT1(ISB1)='K'
     ISB1=ISB1+1 ; CSTITLE1(ISB1)='T_MNW_WATER'//YSUF(JI) ; CSUNIT1(ISB1)='K'
     ISB1=ISB1+1 ; CSTITLE1(ISB1)='T_BOT_WATER'//YSUF(JI) ; CSUNIT1(ISB1)='K'
     ISB1=ISB1+1 ; CSTITLE1(ISB1)='CT_WATER'//YSUF(JI) ; CSUNIT1(ISB1)='-'
     ISB1=ISB1+1 ; CSTITLE1(ISB1)='HML_WATER'//YSUF(JI) ; CSUNIT1(ISB1)='m'
  ENDIF
  ! end SURFACE FIELDS
END DO
!
IF (LWMINMAX) THEN
    DO JI=1,ISER
    ! Max of vertical speed
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='WMAX'//YSUF(JI) ; CSUNIT1(ISB1)='M/S' ; NSGRIDD1(ISB1)=4
    ! Min of vertical speed
    ISB1=ISB1+1 ; CSTITLE1(ISB1)='WMIN'//YSUF(JI) ; CSUNIT1(ISB1)='M/S' ; NSGRIDD1(ISB1)=4
  END DO
END IF
!
IF (ISB1.NE.NSTEMP_SERIE1) THEN
  WRITE(ILUOUT,FMT=*) 'STOP in INI_SERIESn:'
  WRITE(UNIT=ILUOUT,FMT=*) ' NUMBER OF SERIES1 DIFFERS FROM ALLOC, ISB1=', &
                            ISB1,' NSTEMP_SERIE1=',NSTEMP_SERIE1
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SERIES_n','')
END IF
!
!*       2.2   Temporal series (z,t)
!              ---------------------
!
DO JI=1,NSTEMP_SERIE2
  WRITE(CSCOMMENT2(JI),'("TEMPORAL SERIE (z,t) : ",A,".",I0,".",A)') CEXP,IMI,CSEG
  IF (JI==1) WRITE(ILUOUT,FMT=*) CSCOMMENT2(JI)
END DO
!
NSGRIDD2(:)=1
!
ISB2=0
DO JI=1,ISER
  ! Vertical velocity
  ISB2=ISB2+1 ; CSTITLE2(ISB2)='WT'//YSUF(JI) ; CSUNIT2(ISB2)='M/S' ; NSGRIDD2(ISB2)=4
  ! Potential temperature
  ISB2=ISB2+1 ; CSTITLE2(ISB2)='THT'//YSUF(JI) ; CSUNIT2(ISB2)='K'
  ! Pressure
  ISB2=ISB2+1 ; CSTITLE2(ISB2)='PABST'//YSUF(JI) ; CSUNIT2(ISB2)='Pa'
  ! Mixing ratios
  IF (LUSERV) THEN
    ISB2=ISB2+1 ; CSTITLE2(ISB2)='RVT'//YSUF(JI) ; CSUNIT2(ISB2)='KG/KG'
  END IF
  IF (LUSERC) THEN
    ISB2=ISB2+1 ; CSTITLE2(ISB2)='RCT'//YSUF(JI) ; CSUNIT2(ISB2)='KG/KG'
  END IF
  IF (LUSERR) THEN
    ISB2=ISB2+1 ; CSTITLE2(ISB2)='RRT'//YSUF(JI) ; CSUNIT2(ISB2)='KG/KG'
  END IF
  IF (LUSERI) THEN
    ISB2=ISB2+1 ; CSTITLE2(ISB2)='RIT'//YSUF(JI) ; CSUNIT2(ISB2)='KG/KG'
  END IF
  IF (LUSERS) THEN
    ISB2=ISB2+1 ; CSTITLE2(ISB2)='RST'//YSUF(JI) ; CSUNIT2(ISB2)='KG/KG'
  END IF
  IF (LUSERG) THEN
    ISB2=ISB2+1 ; CSTITLE2(ISB2)='RGT'//YSUF(JI) ; CSUNIT2(ISB2)='KG/KG'
  END IF
END DO
!
IF (ISB2.NE.NSTEMP_SERIE2) THEN
  WRITE(ILUOUT,FMT=*) 'STOP in INI_SERIESn:'
  WRITE(ILUOUT,FMT=*) ' NUMBER OF SERIES2 DIFFERS FROM ALLOC, ISB2=',ISB2, &
                      ' NSTEMP_SERIE2=',NSTEMP_SERIE2
!callabortstop
   CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SERIES_n','')
END IF
!
!*    2.3 Temporal series (x,t)
!
DO JI=1,NSTEMP_SERIE3
  WRITE(CSCOMMENT3(JI),'("TEMPORAL SERIE (x,t) : ",A,".",I0,".",A)') CEXP,IMI,CSEG 
  IF (JI==1) WRITE(ILUOUT,FMT=*) CSCOMMENT3(JI)
END DO

!
WRITE(CSKCLS,'(I3.3)') NKCLS
WRITE(CSKCLA,'(I3.3)') NKCLA
WRITE(CSKUP, '(I3.3)') NKUP
WRITE(CSKLOW,'(I3.3)') NKLOW
WRITE(CSKMID,'(I3.3)') NKMID
!
NSGRIDD3(:)=1
!
ISB3=0
!
! U in CLS
ISB3=ISB3+1; CSTITLE3(ISB3)='UCLS'//CSKCLS ; CSUNIT3(ISB3)='M/S' ; NSGRIDD3(ISB3)=2
! W in CLA
ISB3=ISB3+1; CSTITLE3(ISB3)='WCLA'//CSKCLA ; CSUNIT3(ISB3)='M/S' ; NSGRIDD3(ISB3)=4
! W averaged in mid troposphere (between KLOW and KUP)
ISB3=ISB3+1; CSTITLE3(ISB3)='W'   //CSKLOW//'-'//CSKUP ; CSUNIT3(ISB3) ='M/S' ; NSGRIDD3(ISB3)=4
! mixing ratios
IF (LUSERV) THEN
  ! RV in CLS
  ISB3=ISB3+1; CSTITLE3(ISB3)='RVCLS'//CSKCLS ; CSUNIT3(ISB3)='KG/KG'
  ! RV in MID troposphere
  ISB3=ISB3+1; CSTITLE3(ISB3)='RVMID'//CSKMID ; CSUNIT3(ISB3)='KG/KG'
END IF
IF (LUSERC) THEN
  ! RC averaged between ground and KUP
  ISB3=ISB3+1    ; CSTITLE3(ISB3)='RC'//'0-'//CSKUP ; CSUNIT3(ISB3)='KG/KG'
END IF
IF (LUSERR) THEN
  ! RR in CLS
  ISB3=ISB3+1    ; CSTITLE3(ISB3)='RR'//CSKCLS      ; CSUNIT3(ISB3)='KG/KG'
END IF
!
IF (ISB3.NE.NSTEMP_SERIE3) THEN
  WRITE(ILUOUT,FMT=*) 'STOP in INI_SERIESn:'
  WRITE(ILUOUT,FMT=*) ' NUMBER OF SERIES3 DIFFERS FROM ALLOC, ISB3=',ISB3, &
                      ' NTEMP_SERIE3=',NSTEMP_SERIE3
!callabortstop
  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_SERIES_n','')
END IF
!
!-------------------------------------------------------------------------------
!
!*      3.   Initialization of time
!            ----------------------
!
NSCOUNTD=0                         ! Counting the nb of temporal series outputs
!
XSDATIME( 1,:)= TDTEXP%TDATE%YEAR
XSDATIME( 2,:)= TDTEXP%TDATE%MONTH
XSDATIME( 3,:)= TDTEXP%TDATE%DAY
XSDATIME( 4,:)= TDTEXP%TIME
XSDATIME( 5,:)= TDTSEG%TDATE%YEAR
XSDATIME( 6,:)= TDTSEG%TDATE%MONTH
XSDATIME( 7,:)= TDTSEG%TDATE%DAY
XSDATIME( 8,:)= TDTSEG%TIME
XSDATIME( 9,:)= TDTMOD%TDATE%YEAR
XSDATIME(10,:)= TDTMOD%TDATE%MONTH
XSDATIME(11,:)= TDTMOD%TDATE%DAY
XSDATIME(12,:)= TDTMOD%TIME
!
END SUBROUTINE INI_SERIES_n 
