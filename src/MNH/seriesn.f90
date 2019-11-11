!MNH_LIC Copyright 2002-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      SUBROUTINE SERIES_n
!     ###################
!
!!****  *SERIES* - routine to compute diagnostics for diachro files (temporal series)
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
!      V. Ducrocq from contributions of M. Tomasini, S. Donier,... * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    4/03/2002
!!         Oct. 2011 : (P.Le Moigne) Surface series
!!  06/2016     (G.Delautier) phasage surfex 8
!!  01/2018      (G.Delautier) SURFEX 8.1
!!  03/2018     (P.Wautelet)   replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODD_SERIES
USE MODD_SERIES_n
USE MODD_PARAMETERS
USE MODD_CONF, ONLY: NVERB
USE MODD_REF, ONLY: XRHODREFZ
USE MODD_TIME, ONLY: TDTEXP
USE MODD_TYPE_DATE
USE MODD_CONF_n,   ONLY: LUSERV,LUSERC,LUSERR,LUSERI,LUSERS,LUSERG,LUSERH
USE MODD_FIELD_n,  ONLY: XTHT,XWT,XUT,XPABST,XRT
USE MODD_GRID_n,   ONLY: XZZ
USE MODD_LUNIT_n,  ONLY: TLUOUT
USE MODD_PRECIP_n, ONLY: XINPRC,XINPRR,XINPRS,XINPRG,XINPRH, &
                         XACPRC,XACPRR,XACPRS,XACPRG,XACPRH
USE MODD_TIME_n, ONLY: TDTCUR
! SURFACE FIELDS
USE MODI_GET_SURF_VAR_n
!
USE MODE_DATETIME
USE MODE_IO_ll
USE MODE_ll
USE MODE_MSG
!
USE MODD_MNH_SURFEX_n
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!*       0.2     Local variables
!
INTEGER                              :: IIB,IJB,IKB     ! Begin of physical dimensions
INTEGER                              :: IIE,IJE,IKE     ! End   of physical dimensions
INTEGER                              :: IIU,IJU,IKU
INTEGER                              :: IKMAX     !
INTEGER                              :: IISL, IISH      ! Lower and upper  limits along x axis
                                                        !   of y-slice for the sub-domain
INTEGER                              :: IJSL, IJSH      ! Lower and upper  limits along y axis
                                                        ! of y-slice for the sub-domain
INTEGER                              :: IIDIM1,IKM
INTEGER                              :: JI,JK,JS        ! Loop indices
INTEGER                              :: ISB1,ISB2,ISB3  ! the current process for
                                                        ! 1st,2nd and 3rd group respectively
REAL,    DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2))  :: ZDIA
INTEGER  :: ILUOUT ! Logical unit number for output-listing
INTEGER  :: IRESP  ! Return code of FM-routines
INTEGER  :: ISER
CHARACTER (LEN=5), DIMENSION(3) :: YSUF
LOGICAL, DIMENSION(SIZE(LINBOX,1),SIZE(LINBOX,2),3) :: GINBOX
TYPE (DATE_TIME)  :: TZDTCUR ! current date and time
!SURFACE FIELDS
REAL,    DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2))  :: ZTS, ZTMNW, ZTBOT, ZCT,ZHML
INTEGER                                        :: ILOOP, JLOOP, KI
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSERIES
!

!-------------------------------------------------------------------------------
!
!*      1.   INITIALIZATIONS
!            ---------------
!
ILUOUT = TLUOUT%NLU
!
IKU=SIZE(XTHT,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_DIM_EXT_ll ('B', IIU,IJU)
IKB=1+JPVEXT
IKE = IKU - JPVEXT
IKMAX= IKU - 2*JPVEXT
!
ISER=1
WRITE(ILUOUT,*)'LMASKLANDSEA=',LMASKLANDSEA
IF (LMASKLANDSEA) ISER=3
YSUF(1)='-GLOB'
YSUF(2)='-LAND'
YSUF(3)='-SEA '
WRITE(ILUOUT,*)'YSUF=',YSUF
GINBOX(:,:,1)=LINBOX(:,:)
IF (LMASKLANDSEA) THEN
  GINBOX(:,:,2)=LINBOXL(:,:)
  GINBOX(:,:,3)=LINBOXS(:,:)
END IF
!SURFACE FIELDS
IF (LSURF) THEN
   KI=(IIE-IIB+1)*(IJE-IJB+1)
   ALLOCATE(ZSERIES(KI,5))
   CALL GET_SURF_VAR_n(YSURF_CUR%FM,YSURF_CUR%IM,YSURF_CUR%SM,YSURF_CUR%TM, &
                       YSURF_CUR%WM,YSURF_CUR%DUO,YSURF_CUR%DU,YSURF_CUR%UG,& 
                       YSURF_CUR%U, YSURF_CUR%USS,&
                       'MESONH',KI,5,PSERIES=ZSERIES)
   ZTS(:,:)=0.
   ZTMNW(:,:)=0.
   ZTBOT(:,:)=0.
   ZCT(:,:)=0.
   ZHML(:,:)=0.
   DO JLOOP=IJB,IJE
      DO ILOOP=IIB,IIE
         ZTS(ILOOP,JLOOP)=ZSERIES(ILOOP-1+(IIE-1)*(JLOOP-2),1)
         ZTMNW(ILOOP,JLOOP)=ZSERIES(ILOOP-1+(IIE-1)*(JLOOP-2),2)
         ZTBOT(ILOOP,JLOOP)=ZSERIES(ILOOP-1+(IIE-1)*(JLOOP-2),3)
         ZCT(ILOOP,JLOOP)=ZSERIES(ILOOP-1+(IIE-1)*(JLOOP-2),4)
         ZHML(ILOOP,JLOOP)=ZSERIES(ILOOP-1+(IIE-1)*(JLOOP-2),5)
      ENDDO
   ENDDO
   IF(NVERB==10) THEN
      DO JLOOP=IJB-1,IJE+1
         DO ILOOP=IIB-1,IIE+1
            WRITE(ILUOUT,*)'TS   : ',ILOOP,JLOOP,ZTS(ILOOP,JLOOP)
            WRITE(ILUOUT,*)'TMNW : ',ILOOP,JLOOP,ZTMNW(ILOOP,JLOOP)
            WRITE(ILUOUT,*)'TBOT : ',ILOOP,JLOOP,ZTBOT(ILOOP,JLOOP)
            WRITE(ILUOUT,*)'CT   : ',ILOOP,JLOOP,ZCT(ILOOP,JLOOP)
            WRITE(ILUOUT,*)'HML  : ',ILOOP,JLOOP,ZHML(ILOOP,JLOOP)
         ENDDO
      ENDDO
   ENDIF
ENDIF
!end SURFACE FIELDS
!
IF(NVERB>=5) WRITE(ILUOUT,*) &
             'SERIESn: sizes ', IIB,IJB,IIE,IJE,IKB,IKE,IIU,IJU,IKU,IKMAX
!
!*      1.1  initializes current time and time counter
!
NSCOUNTD=NSCOUNTD+1
!
TZDTCUR=TDTCUR
!
CALL DATETIME_DISTANCE(TDTEXP,TZDTCUR,XSTRAJT(NSCOUNTD,1))
XSDATIME(13,NSCOUNTD)= TZDTCUR%TDATE%YEAR
XSDATIME(14,NSCOUNTD)= TZDTCUR%TDATE%MONTH
XSDATIME(15,NSCOUNTD)= TZDTCUR%TDATE%DAY
XSDATIME(16,NSCOUNTD)= TZDTCUR%TIME
!
!-------------------------------------------------------------------------------
!
!*      2.   Temporal series t
!            -----------------
!
IF (LSERIES1) THEN
   ISB1=0
   DO JI=1,ISER
!
!*      2.1  Instant total explicit precipitation (conversion: m/s -> mm/day)
!
     IF (SIZE(XINPRR)/= 0) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='INPRT'//YSUF(JI)) THEN
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( XINPRR(IIB:IIE,IJB:IJE),      &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )  &
                                           * 3.6E6*24.
         IF (SIZE(XINPRS)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
              (SUM( XINPRS(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
              ) * 3.6E6*24.
         IF (SIZE(XINPRG)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
              (SUM( XINPRG(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
              ) * 3.6E6*24.
         IF (SIZE(XINPRH)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
              (SUM( XINPRH(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
              ) * 3.6E6*24.
         IF (SIZE(XINPRC)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
              (SUM( XINPRC(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
              ) * 3.6E6*24.
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'INPRT','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
!
!*      2.2  Accumulated total explicit precipitation (conversion: m -> mm)
!
     IF (SIZE(XINPRR)/= 0) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='ACPRT'//YSUF(JI)) THEN
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( XACPRR(IIB:IIE,IJB:IJE),      &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )  &
                                           * 1000.
         IF (SIZE(XINPRS)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
               SUM( XACPRS(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
               * 1000.
         IF (SIZE(XINPRG)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
               SUM( XACPRG(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
               * 1000.
         IF (SIZE(XINPRH)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
               SUM( XACPRH(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
               * 1000.
         IF (SIZE(XINPRC)/= 0) XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) =    &
                               XSSERIES1(1,1,1,NSCOUNTD,1,ISB1) +    &
               SUM( XACPRC(IIB:IIE,IJB:IJE), MASK=GINBOX(IIB:IIE,IJB:IJE,JI) ) &
               * 1000.
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'ACPRT','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
!
!*      2.3   Mixing ratios
!
     IF (LUSERV) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='RVT'//YSUF(JI)) THEN
         ZDIA   (:,:)  =0.
         ! Average RV
         DO JK=IKB,IKE
           ZDIA(:,:)= XRT(:,:,JK,1)*(XZZ(:,:,JK+1)-XZZ(:,:,JK))*XRHODREFZ(JK)&
                       +ZDIA(:,:)
         END DO
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( ZDIA(IIB:IIE,IJB:IJE),       &
                                         MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'RVT  ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
     !
     IF (LUSERC) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='RCT'//YSUF(JI)) THEN
         ZDIA   (:,:)  =0.
         ! Average RC
         DO JK=IKB,IKE
           ZDIA(:,:)= XRT(:,:,JK,2)*(XZZ(:,:,JK+1)-XZZ(:,:,JK))*XRHODREFZ(JK)&
                       +ZDIA(:,:)
         END DO
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( ZDIA(IIB:IIE,IJB:IJE),     &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'RCT  ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
     !
     IF (LUSERR) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='RRT'//YSUF(JI)) THEN
         ZDIA   (:,:)  =0.
         ! Average RR
         DO JK=IKB,IKE
           ZDIA(:,:)= XRT(:,:,JK,3)*(XZZ(:,:,JK+1)-XZZ(:,:,JK))*XRHODREFZ(JK)&
                      +ZDIA(:,:)
         END DO
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( ZDIA(IIB:IIE,IJB:IJE),     &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'RRT  ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
     !
     IF (LUSERI) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='RIT'//YSUF(JI)) THEN
         ZDIA   (:,:)  =0.
         ! Average RI
         DO JK=IKB,IKE
           ZDIA(:,:)= XRT(:,:,JK,4)*(XZZ(:,:,JK+1)-XZZ(:,:,JK))*XRHODREFZ(JK)&
                      +ZDIA(:,:)
         END DO
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( ZDIA(IIB:IIE,IJB:IJE),     &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'RIT  ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
     !
     IF (LUSERS) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='RST'//YSUF(JI)) THEN
         ZDIA   (:,:)  =0.
         ! Average RS
         DO JK=IKB,IKE
           ZDIA(:,:)= XRT(:,:,JK,5)*(XZZ(:,:,JK+1)-XZZ(:,:,JK))*XRHODREFZ(JK)&
                      +ZDIA(:,:)
         END DO
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( ZDIA(IIB:IIE,IJB:IJE),     &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'RST  ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
     !
     IF (LUSERG) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='RGT'//YSUF(JI)) THEN
         ZDIA   (:,:)  =0.
         ! Average RG
         DO JK=IKB,IKE
           ZDIA(:,:)= XRT(:,:,JK,6)*(XZZ(:,:,JK+1)-XZZ(:,:,JK))*XRHODREFZ(JK)&
                      +ZDIA(:,:)
         END DO
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( ZDIA(IIB:IIE,IJB:IJE),     &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'RGT  ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
     !
     IF (LUSERH) THEN
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='RHT'//YSUF(JI)) THEN
         ZDIA   (:,:)  =0.
         ! Average RG
         DO JK=IKB,IKE
           ZDIA(:,:)= XRT(:,:,JK,7)*(XZZ(:,:,JK+1)-XZZ(:,:,JK))*XRHODREFZ(JK)&
                      +ZDIA(:,:)
         END DO
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM( ZDIA(IIB:IIE,IJB:IJE),     &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'RHT  ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     END IF
!SURFACE FIELDS
     IF (LSURF) THEN
        ISB1=ISB1+1
        IF (TRIM(CSTITLE1(ISB1))=='TS_WATER'//YSUF(JI)) THEN
          XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM(ZTS(IIB:IIE,IJB:IJE),           &
                                                MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
        END IF
        ISB1=ISB1+1
        IF (TRIM(CSTITLE1(ISB1))=='T_MNW_WATER'//YSUF(JI)) THEN
          XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM(ZTMNW(IIB:IIE,IJB:IJE),           &
                                                MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
        END IF
        ISB1=ISB1+1
        IF (TRIM(CSTITLE1(ISB1))=='T_BOT_WATER'//YSUF(JI)) THEN
          XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM(ZTBOT(IIB:IIE,IJB:IJE),           &
                                                MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
        END IF
        ISB1=ISB1+1
        IF (TRIM(CSTITLE1(ISB1))=='CT_WATER'//YSUF(JI)) THEN
          XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM(ZCT(IIB:IIE,IJB:IJE),           &
                                                MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
        END IF
        ISB1=ISB1+1
        IF (TRIM(CSTITLE1(ISB1))=='HML_WATER'//YSUF(JI)) THEN
          XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)= SUM(ZHML(IIB:IIE,IJB:IJE),           &
                                                MASK=GINBOX(IIB:IIE,IJB:IJE,JI) )
        END IF
     ENDIF
!end SURFACE FIELDS
   END DO
!
!*      2.4   Vertical velocity
!
   DO JI=1,ISER
     IF (LWMINMAX) THEN
       ! Max vertical velocity
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='WMAX'//YSUF(JI)) THEN
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)=MAXVAL(XWT(IIB:IIE,IJB:IJE,IKB:IKE), &
               MASK=SPREAD(GINBOX(IIB:IIE,IJB:IJE,JI),DIM=3,NCOPIES=IKMAX))
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'WMAX ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
       ! Min vertical velocity
       ISB1=ISB1+1
       IF (TRIM(CSTITLE1(ISB1))=='WMIN'//YSUF(JI)) THEN
         XSSERIES1(1,1,1,NSCOUNTD,1,ISB1)=MINVAL(XWT(IIB:IIE,IJB:IJE,IKB:IKE), &
              MASK=SPREAD(GINBOX(IIB:IIE,IJB:IJE,JI),DIM=3,NCOPIES=IKMAX))
       ELSE
         WRITE(UNIT=ILUOUT,FMT=1) 'WMIN ','1',CSTITLE1(ISB1)
 !callabortstop
         CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
       END IF
     ENDIF
   END DO
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      3.   Temporal series (z,t)
!            ---------------------
!
IF (LSERIES2) THEN
  ISB2=0
  DO JI=1,ISER
!
!*      3.1  Vertical velocity
!
    ISB2=ISB2+1
    IF (TRIM(CSTITLE2(ISB2))=='WT'//YSUF(JI)) THEN
      DO JK=IKB,IKE
        IKM=JK-IKB +1
        XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM( XWT   (IIB:IIE,IJB:IJE,JK  ), &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
      END DO
    ELSE
      WRITE(UNIT=ILUOUT,FMT=1) 'WT   ','2',CSTITLE2(ISB2)
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
    END IF
!
!*      3.2  Potential Temperature
!
    ISB2=ISB2+1
    IF (TRIM(CSTITLE2(ISB2))=='THT'//YSUF(JI)) THEN
      DO JK=IKB,IKE
        IKM=JK-IKB +1
        XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM( XTHT  (IIB:IIE,IJB:IJE,JK  ), &
                                           MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
      END DO
    ELSE
      WRITE(UNIT=ILUOUT,FMT=1) 'THT  ','2',CSTITLE2(ISB2)
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
    END IF
!
!*      3.3  Pressure
!
    ISB2=ISB2+1
    IF (TRIM(CSTITLE2(ISB2))=='PABST'//YSUF(JI)) THEN
      DO JK=IKB,IKE
        IKM=JK-IKB +1
         XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM(XPABST (IIB:IIE,IJB:IJE,JK  ), &
                                            MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
      END DO
    ELSE
      WRITE(UNIT=ILUOUT,FMT=1) 'PASBT','2',CSTITLE2(ISB2)
 !callabortstop
      CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
    END IF
!
!*      3.4  Vapor
!
    IF (LUSERV) THEN
      ISB2=ISB2+1
      IF (TRIM(CSTITLE2(ISB2))=='RVT'//YSUF(JI)) THEN
        DO JK=IKB,IKE
          IKM=JK-IKB +1
          XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM(  XRT  (IIB:IIE,IJB:IJE,JK,1),&
                                             MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
        END DO
      ELSE
        WRITE(UNIT=ILUOUT,FMT=1) 'RVT  ','2',CSTITLE2(ISB2)
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
      END IF
    END IF
!
!*      3.5  Cloud water
!
    IF (LUSERC) THEN
      ISB2=ISB2+1
      IF (TRIM(CSTITLE2(ISB2))=='RCT'//YSUF(JI)) THEN
        DO JK=IKB,IKE
          IKM=JK-IKB +1
          XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM(  XRT  (IIB:IIE,IJB:IJE,JK,2),&
                                             MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
        END DO
      ELSE
        WRITE(UNIT=ILUOUT,FMT=1) 'RCT  ','2',CSTITLE2(ISB2)
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
      END IF
    END IF
!
!*      3.6  Rainwater
!
    IF (LUSERR) THEN
      ISB2=ISB2+1
      IF (TRIM(CSTITLE2(ISB2))=='RRT'//YSUF(JI)) THEN
        DO JK=IKB,IKE
          IKM=JK-IKB +1
           XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM(  XRT (IIB:IIE,IJB:IJE,JK,3),&
                                              MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
        END DO
      ELSE
        WRITE(UNIT=ILUOUT,FMT=1) 'RRT  ','2',CSTITLE2(ISB2)
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
      END IF
    END IF
!
!*      3.7  Primary ice
!
    IF (LUSERI) THEN
      ISB2=ISB2+1
      IF (TRIM(CSTITLE2(ISB2))=='RIT'//YSUF(JI)) THEN
        DO JK=IKB,IKE
          IKM=JK-IKB +1
          XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM(  XRT (IIB:IIE,IJB:IJE,JK,4), &
                                             MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
        END DO
      ELSE
        WRITE(UNIT=ILUOUT,FMT=1) 'RIT  ','2',CSTITLE2(ISB2)
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
      END IF
    END IF
!
!*      3.7  Snow
!
    IF (LUSERS) THEN
      ISB2=ISB2+1
      IF (TRIM(CSTITLE2(ISB2))=='RST'//YSUF(JI)) THEN
        DO JK=IKB,IKE
          IKM=JK-IKB +1
          XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM(  XRT (IIB:IIE,IJB:IJE,JK,5), &
                                             MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
        END DO
      ELSE
        WRITE(UNIT=ILUOUT,FMT=1) 'RST  ','2',CSTITLE2(ISB2)
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
      END IF
    END IF
!
!*      3.8  Graupel
!
    IF (LUSERG) THEN
      ISB2=ISB2+1
      IF (TRIM(CSTITLE2(ISB2))=='RGT'//YSUF(JI)) THEN
        DO JK=IKB,IKE
          IKM=JK-IKB +1
          XSSERIES2(1,1,IKM,NSCOUNTD,1,ISB2)=SUM(  XRT (IIB:IIE,IJB:IJE,JK,6), &
                                             MASK=GINBOX(IIB:IIE,IJB:IJE,JI))
        END DO
      ELSE
        WRITE(UNIT=ILUOUT,FMT=1) 'RGT  ','2',CSTITLE2(ISB2)
 !callabortstop
        CALL PRINT_MSG(NVERB_FATAL,'GEN','SERIES_n','')
      END IF
    END IF
!
  END DO
END IF
!
1 FORMAT("SERIES: ORDER IS NOT CORRECT FOR ",A5,".CSTITLE",A1,'=',A10)
!-------------------------------------------------------------------------------
!
!*      4.   Temporal series (x,t)
!            ---------------------
!
IIDIM1=0
IF (NBJSLICE > 0 ) THEN
  !
  DO JS=1,NBJSLICE
    IF (LSERIES3(JS)) THEN
      IIDIM1=NISH(JS) - NISL(JS) + 1
      IISL=NISL(JS)
      IISH=NISH(JS)
      IJSL=NJSLICESL(JS)
      IJSH=NJSLICESH(JS)
      !
      ISB3= (JS-1) * NSTEMP_SERIE3
!
!*      4.1  U at level KCLS
!
      ISB3=ISB3+1
      XSSERIES3(1:IIDIM1,1,1,NSCOUNTD,1,ISB3)=                         &
                          SUM( XUT (IISL:IISH,IJSL:IJSH,NKCLS+JPVEXT),  DIM=2 )
!
!*      4.2  W at level KCLA
!
      ISB3=ISB3+1
      XSSERIES3(1:IIDIM1,1,1,NSCOUNTD,1,ISB3)=                         &
                          SUM( XWT (IISL:IISH,IJSL:IJSH,NKCLA+JPVEXT),  DIM=2 )
!
!*      4.3  Mean W between KLOW and KUP
!
      ISB3=ISB3+1
      ZDIA (:,:)=0.
      DO JK=NKLOW+JPVEXT,NKUP+JPVEXT
        ZDIA(:,:)= XWT(:,:,JK)*(XZZ(:,:,JK+1)-XZZ(:,:,JK)) + ZDIA(:,:)
      END DO
      ZDIA(:,:)= ZDIA(:,:)/(XZZ(:,:,NKUP+JPVEXT+1)-XZZ(:,:,NKLOW+JPVEXT))
      XSSERIES3(1:IIDIM1,1,1,NSCOUNTD,1,ISB3)=                         &
                              SUM( ZDIA (IISL:IISH,IJSL:IJSH) , DIM=2 )
!
      IF (LUSERV) THEN
!
!*      4.4  RV  at level KCLS
!
        ISB3=ISB3+1
        XSSERIES3(1:IIDIM1,1,1,NSCOUNTD,1,ISB3)=                       &
                          SUM( XRT (IISL:IISH,IJSL:IJSH,NKCLS+JPVEXT,1),DIM=2 )
!
!*      4.5  Rv at level KMID
!
        ISB3=ISB3+1
        XSSERIES3(1:IIDIM1,1,1,NSCOUNTD,1,ISB3)=                       &
                          SUM( XRT (IISL:IISH,IJSL:IJSH,NKMID+JPVEXT,1),DIM=2 )
!
      END IF
!
!*      4.6  Mean rc between 0 and KUP
!
      IF (LUSERC) THEN
        ISB3=ISB3+1
        ZDIA (:,:)=0.
        DO JK=IKB,NKUP+JPVEXT
          ZDIA(:,:)= XRT(:,:,JK,2)*(XZZ(:,:,JK+1)-XZZ(:,:,JK)) + ZDIA(:,:)
        END DO
        ZDIA(:,:)= ZDIA(:,:)/(XZZ(:,:,NKUP+JPVEXT+1)-XZZ(:,:,IKB))
        XSSERIES3(1:IIDIM1,1,1,NSCOUNTD,1,ISB3)=                     &
                             SUM( ZDIA (IISL:IISH,IJSL:IJSH) , DIM=2 )
      END IF
!
!*      4.7  Rr at level KCLS
!
      IF (LUSERR) THEN
        ISB3=ISB3+1
        XSSERIES3(1:IIDIM1,1,1,NSCOUNTD,1,ISB3)=                       &
                         SUM( XRT(IISL:IISH,IJSL:IJSH,NKCLS+JPVEXT,3) , DIM=2 )
      END IF
    END IF
  END DO
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SERIES_n 
