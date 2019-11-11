!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODI_PASPOL
!    ################## 
!
INTERFACE
!
      SUBROUTINE PASPOL (PTSTEP, PSFSV, KLUOUT, KVERB, OCLOSE_OUT, TPFILE)
!
USE MODD_IO_ll, ONLY: TFILEDATA
!
IMPLICIT NONE
!
REAL,                   INTENT(IN)    :: PTSTEP     ! Double timestep except 
                                                    ! for the first time step (single one)
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSFSV      ! surface flux of scalars
INTEGER,                INTENT(IN)    :: KLUOUT     ! unit for output listing count
INTEGER,                INTENT(IN)    :: KVERB      ! verbosity level
LOGICAL,                INTENT(IN)    :: OCLOSE_OUT ! conditional closure of the OUTPUT FM-file
TYPE(TFILEDATA),        INTENT(IN)    :: TPFILE     ! Output file
!
END SUBROUTINE PASPOL
!
END INTERFACE
!
END MODULE MODI_PASPOL
!     ######spl
      SUBROUTINE PASPOL (PTSTEP, PSFSV, KLUOUT, KVERB, OCLOSE_OUT, TPFILE)
!     ############################################################
!
!
!
!!****  *PASPOL* -
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to release a passive tracer 
!
!!**  METHOD
!!    ------
!!    
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      M.Bouzom, C.Lac         * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    C.Lac 30/09/2010 Bugs for reproducibility : position of release +
!!                            GET_INDICE_ll replaced by GET_PHYSICAL_ll +
!!                            remove the diffusion at the release 
!!    C.Lac 11/11 Remove instant M
!!    P.Wautelet 28/03/2018 Replace TEMPORAL_DIST by DATETIME_DISTANCE
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!! --------------------------------------------------------------------------
!       
!!    EXTERNAL
!!    --------
!!
USE MODD_PARAMETERS
USE MODD_NSV
USE MODD_CST
USE MODD_IO_ll, ONLY: TFILEDATA
USE MODE_GRIDPROJ
USE MODD_PASPOL
USE MODD_CTURB
USE MODI_SHUMAN
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_M
USE MODE_FMWRIT
USE MODE_ll
USE MODE_FM
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_DYN_n
USE MODD_CONF
USE MODD_GRID
USE MODD_CONF_n
USE MODD_REF_n
USE MODD_FIELD_n
!
USE MODD_PASPOL_n
USE MODD_GRID_n
USE MODD_TIME_n
USE MODD_SUB_PASPOL_n
USE MODD_TYPE_DATE
!
USE MODE_DATETIME
USE MODE_FIELD, ONLY: TFIELDDATA,TYPEREAL
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL,                   INTENT(IN)    :: PTSTEP     ! Double timestep except 
                                                    ! for the first time step (single one)
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSFSV      ! surface flux of scalars
INTEGER,                INTENT(IN)    :: KLUOUT     ! unit for output listing count
INTEGER,                INTENT(IN)    :: KVERB      ! verbosity level
LOGICAL,                INTENT(IN)    :: OCLOSE_OUT ! conditional closure of the OUTPUT FM-file
TYPE(TFILEDATA),        INTENT(IN)    :: TPFILE     ! Output file
!
!*      0.2    declarations of local variables
!
INTEGER :: IIB,IIE,IJB,IJE             !intersection of the physical global domain
! with the local extended subdomain (in local indices)
INTEGER :: IKB, IKE
INTEGER :: IIU, IJU, IKU               ! dimensional indexes
INTEGER :: JK ! Loop indice
INTEGER :: JSV,IBOT,ITOP,II,IJ,IP
REAL    :: ZSRCX,ZSRCY,ZSRCI,ZSRCJ,ZTOP,ZGROUND,ZX,ZY,ZSOMME
INTEGER :: I1YY,I1MM,I1DD,I1HH,I1MN,I1SS,I2YY,I2MM,I2DD,I2HH,I2MN,I2SS
INTEGER :: I3YY,I3MM,I3DD,I3HH,I3MN,I3SS,I4YY,I4MM,I4DD,I4HH,I4MN,I4SS
REAL    :: Z1SEC,Z2SEC,Z3SEC,Z4SEC,ZT1,ZT2,ZT3,ZT4
REAL    :: ZCBOT,ZCTOP,ZEPAIS,ZENBAS,ZENHAUT,ZDEPUIS,ZRATE
!
REAL    :: ZSURF  ! Surface of a grid box            
!
!
REAL    :: ZP, ZTH, ZT, ZRHO, ZMASAIR
!REAL,  DIMENSION(-1:2,-1:2)  :: ZP      ! Pressure                    (Pa)
!REAL,  DIMENSION(-1:2,-1:2)  :: ZTH     ! Potential temperature        (K)
!REAL,  DIMENSION(-1:2,-1:2)  :: ZT      ! Temperature                  (K)
!REAL,  DIMENSION(-1:2,-1:2)  :: ZRHO    ! Rho                       (g/m3)
!REAL,  DIMENSION(-1:2,-1:2)  :: ZMASAIR ! Masse of air in the grid box (g)
!
!REAL, DIMENSION( 0:1, 0:1) :: Z4PT
!REAL, DIMENSION(-1:1,-1:1) :: Z9PT
!REAL                       :: ZDELTAX,ZDELTAY
!INTEGER                    :: J4PTI,J4PTJ,J9PTI,J9PTJ
!
REAL,  DIMENSION(:,:,:), ALLOCATABLE :: ZRHOM  ! 
REAL,  DIMENSION(:,:,:), ALLOCATABLE :: ZTEMPO, ZSVT ! Work arrays
!
TYPE(DATE_TIME)   :: TZDATE1,TZDATE2,TZDATE3,TZDATE4,TZDATE
TYPE(TFIELDDATA)  :: TZFIELD
!
!
!--------------------------------------------------------------------------------------
!
!
!*	0. Initialisation
!
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_PHYSICAL_ll (IIB,IJB,IIE,IJE)
!
IKU=SIZE(XTHT,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
ALLOCATE( ZRHOM(IIU,IJU,IKU) )
ALLOCATE( ZSVT(IIU,IJU,IKU) )
!
ZSURF = (XXHAT(2)-XXHAT(1))*(XYHAT(2)-XYHAT(1))   ! Surface d'une maille.
!
!
!*	1.  INITIALIZATION OF PASSIVE POLLUTANT  
!	    -------------------------------------
!

IF (GPPFIRSTCALL) THEN
  GPPFIRSTCALL = .FALSE.
   !
  ALLOCATE( GBEGEMIS(NSV_PP) )
  ALLOCATE( IDEBYY(NSV_PP) )
  ALLOCATE( IDEBMM(NSV_PP) )
  ALLOCATE( IDEBDD(NSV_PP) )
  ALLOCATE( ZDEBSS(NSV_PP) )
  ALLOCATE( ZCHRO2(NSV_PP) )
  ALLOCATE( ZCHRO3(NSV_PP) )
  ALLOCATE( ZCHRO4(NSV_PP) )
  ALLOCATE( IPIGI(NSV_PP) )
  ALLOCATE( IPIGJ(NSV_PP) )
  ALLOCATE( ZQTOT(NSV_PP) )
  ALLOCATE( ZRATEINC(NSV_PP) )
  ALLOCATE( ZRATECST(NSV_PP) )
  ALLOCATE( ZRATEDEC(NSV_PP) )
  ALLOCATE( CNIVO(NSV_PP) )
  ALLOCATE( ZDIFHORZ(-1:2,-1:2,NSV_PP) )
  ALLOCATE( ZREPVERT(IKU,NSV_PP) )
   !
   GBEGEMIS(:) = .FALSE.  ! L'emission n'a pas commence.
   !
 !
   DO JSV=1,NSV_PP
      !
    IF (KVERB >= 10) THEN
      WRITE(KLUOUT,'(A)')   ' '
      WRITE(KLUOUT,'(A)') '******************************'
      WRITE(KLUOUT,'(A)') '	    EMIS_PASPOL'
      WRITE(KLUOUT,'(A)') '  Initialization of sources	 '
      WRITE(KLUOUT,'(A)') '******************************'
      WRITE(KLUOUT,'(A,I3.1)')  'IJU       : ',IJU
      WRITE(KLUOUT,'(A,I3.1)')  'Rejet Nu  : ',JSV
      WRITE(KLUOUT,'(A,F11.6)') 'Longitude : ',XPPLON(JSV)
      WRITE(KLUOUT,'(A,F11.6)') 'Latitude  : ',XPPLAT(JSV)
    END IF
      !
      !
      !*	1.1 Position du rejet.
      !
      ! On calcule les coordonnees cartesiennes (ZSRCX,ZSRCY) en metres,
      ! puis les indices fractionnaires (ZSRCI,ZSRCJ) et entiers
      ! (IPIGI,IPIGJ) du point de rejet dans le domaine de travail global.
      !
      CALL SM_XYHAT(XLATORI,XLONORI,XPPLAT(JSV),XPPLON(JSV),ZSRCX,ZSRCY)
      II=MAX(MIN(COUNT(XXHAT(:)<ZSRCX),IIU-1),1)
      IJ=MAX(MIN(COUNT(XYHAT(:)<ZSRCY),IJU-1),1)
      ZSRCI=(ZSRCX-XXHAT(II))/(XXHAT(II+1)-XXHAT(II))+FLOAT(II)
      ZSRCJ=(ZSRCY-XYHAT(IJ))/(XYHAT(IJ+1)-XYHAT(IJ))+FLOAT(IJ)
      !
      IPIGI(JSV)=INT(ZSRCI)
      IPIGJ(JSV)=INT(ZSRCJ)
      !
    IF (KVERB >= 10) THEN
      WRITE(KLUOUT,'(A,I3.1)')  'IJU       : ',IJU
      WRITE(KLUOUT,'(A,F12.4)') 'Zsrc X (m) : ',ZSRCX
      WRITE(KLUOUT,'(A,F12.4)') 'Zsrc Y (m) : ',ZSRCY
      WRITE(KLUOUT,'(A,F9.3)')  'Ind Rel X  : ',ZSRCI
      WRITE(KLUOUT,'(A,F9.3)')  'Ind Rel Y  : ',ZSRCJ
      WRITE(KLUOUT,'(A,I5.1)')  'Ind X	: ',IPIGI(JSV)
      WRITE(KLUOUT,'(A,I5.1)')  'Ind Y	: ',IPIGJ(JSV)
      !
      WRITE(KLUOUT,'(A,F9.3)')  &
        'Ind X : ',(ZSRCX-XXHAT(1))/(XXHAT(2)-XXHAT(1))+1.0
      WRITE(KLUOUT,'(A,F9.3)')  &
        'Ind Y : ',(ZSRCY-XYHAT(1))/(XYHAT(2)-XYHAT(1))+1.0
     END IF
      !
      !
      IF (IPIGI(JSV).GE.IIB.AND.IPIGI(JSV).LE.IIE.AND.  &
          IPIGJ(JSV).GE.IJB.AND.IPIGJ(JSV).LE.IJE) THEN
         !
        IF (KVERB >= 10) &        
         WRITE(KLUOUT,'(A,I3.1)')  'IJU       : ',IJU
         WRITE(KLUOUT,'(A)') 'La source est dans le domaine de travail courant.'
         !
         !
         !*	1.2 Dispersion autour de la source.
         !
         ! On commence par dispatcher le polluant sur les 4 points
         ! entourant la source (tableau Z4PT).
!
!        ZDELTAX   =  DMOD(ZSRCI,1.0)
!        ZDELTAY   =  DMOD(ZSRCJ,1.0)
!        Z4PT(0,0) = (1.0-ZDELTAX)*(1.0-ZDELTAY)
!        Z4PT(1,0) = (    ZDELTAX)*(1.0-ZDELTAY)
!        Z4PT(0,1) = (1.0-ZDELTAX)*(    ZDELTAY)
!        Z4PT(1,1) = (    ZDELTAX)*(    ZDELTAY)
         !
         ! Pour chaque point precedemment trouve, on opere une 
         ! diffusion sur les 9 points autour avec une formulation
         ! en cos(DX)*cos(DY) (tableau Z9PT). Les resultats sont
         ! integres dans une matrice 4*4 (tableau ZDIFHOR).
         !
!        ZDIFHORZ(:,:,JSV) = 0.
         !
!        DO J4PTI=0,1
!           DO J4PTJ=0,1
               !
!              Z9PT(:,:)=0.0
               !
!              IF (CPPINIT(JSV)=='9PT') THEN
                  !
!                 DO J9PTI= -1,1
!                    ZX=ABS(FLOAT(J9PTI)*XPI/4.) 
!                    DO J9PTJ= -1,1
!                       ZY=ABS(FLOAT(J9PTJ)*XPI/4.)
!                       Z9PT(J9PTI,J9PTJ)=COS(ZX)*COS(ZY)
!                    END DO
!                 END DO
!                 ZSOMME=SUM(Z9PT)
!                 Z9PT(:,:)=Z9PT(:,:)/ZSOMME
                  !
!              ELSE
                  !
                  ! On emet seulement au point courant.
                  !
!                 Z9PT(:,:)=0.0
!                 Z9PT(0,0)=1.0
                  !
!              ENDIF
               !
!        DO J9PTI=-1,1
!                 DO J9PTJ=-1,1
!                    ZDIFHORZ(J4PTI+J9PTI,J4PTJ+J9PTJ,JSV)       &
!                       = ZDIFHORZ(J4PTI+J9PTI,J4PTJ+J9PTJ,JSV)  &
!                       + Z9PT(J9PTI,J9PTJ)*Z4PT(J4PTI,J4PTJ)
!                 END DO
!              END DO
               !
!           END DO
!        END DO
         !
         !
         !
         !*	1.3 Chronologie du rejet.
         !
         ! Eclatement des dates caracteristiques et calcul du
         ! nombre de secondes depuis 0UTC (equiv. TDTCUR%TIME).
         !
         READ(CPPT1(JSV),'(I4,5I2)') I1YY,I1MM,I1DD,I1HH,I1MN,I1SS
         READ(CPPT2(JSV),'(I4,5I2)') I2YY,I2MM,I2DD,I2HH,I2MN,I2SS
         READ(CPPT3(JSV),'(I4,5I2)') I3YY,I3MM,I3DD,I3HH,I3MN,I3SS
         READ(CPPT4(JSV),'(I4,5I2)') I4YY,I4MM,I4DD,I4HH,I4MN,I4SS
         Z1SEC=FLOAT(I1SS+I1MN*60+I1HH*3600)
         Z2SEC=FLOAT(I2SS+I2MN*60+I2HH*3600)
         Z3SEC=FLOAT(I3SS+I3MN*60+I3HH*3600)
         Z4SEC=FLOAT(I4SS+I4MN*60+I4HH*3600)
         !
         ! Chrono relative au debut du rejet en secondes.
         !
         ZT1=0.
         TZDATE1%TDATE%YEAR=I1YY;TZDATE1%TDATE%MONTH=I1MM;TZDATE1%TDATE%DAY=I1DD;TZDATE1%TIME=Z1SEC
         TZDATE2%TDATE%YEAR=I2YY;TZDATE2%TDATE%MONTH=I2MM;TZDATE2%TDATE%DAY=I2DD;TZDATE2%TIME=Z2SEC
         TZDATE3%TDATE%YEAR=I3YY;TZDATE3%TDATE%MONTH=I3MM;TZDATE3%TDATE%DAY=I3DD;TZDATE3%TIME=Z3SEC
         TZDATE4%TDATE%YEAR=I4YY;TZDATE4%TDATE%MONTH=I4MM;TZDATE4%TDATE%DAY=I4DD;TZDATE4%TIME=Z4SEC
         CALL DATETIME_DISTANCE(TZDATE1,TZDATE2,ZT2)
         CALL DATETIME_DISTANCE(TZDATE1,TZDATE3,ZT3)
         CALL DATETIME_DISTANCE(TZDATE1,TZDATE4,ZT4)
         !
         ! On met de cote le debut du rejet sous forme pratique ainsi
         ! que la chronologie relative.
!
         IDEBYY(JSV)    = I1YY
         IDEBMM(JSV)    = I1MM
         IDEBDD(JSV)    = I1DD
         ZDEBSS(JSV)    = Z1SEC
         ZCHRO2(JSV)    = ZT2
         ZCHRO3(JSV)    = ZT3
         ZCHRO4(JSV)    = ZT4
         !
         IF (KVERB >= 10) THEN
           WRITE(KLUOUT,'(A,A14)')   'Debut emission     : ',CPPT1(JSV)
           WRITE(KLUOUT,'(A,A14)')   'Debut pallier      : ',CPPT2(JSV)
           WRITE(KLUOUT,'(A,A14)')   'Fin   pallier      : ',CPPT3(JSV)
           WRITE(KLUOUT,'(A,A14)')   'Fin   emission     : ',CPPT4(JSV)
           WRITE(KLUOUT,'(A,F10.0)') 'Debut          (s) : ',ZT1
           WRITE(KLUOUT,'(A,F10.0)') 'Debut pallier  (s) : ',ZT2
           WRITE(KLUOUT,'(A,F10.0)') 'Fin   pallier  (s) : ',ZT3
           WRITE(KLUOUT,'(A,F10.0)') 'Fin   emission (s) : ',ZT4
         END IF
         !
         !
         !*	1.4 Debit (g/s) et gradients de debit (g/s/s).
         !
         ! Les gradients de debit ne sont calcules que si ils sont significatifs.
         !
         ZRATECST(JSV) =   XPPMASS(JSV)*2.0/(ZT4-ZT1+ZT3-ZT2)
         ZRATEINC(JSV) =   0.0
         ZRATEDEC(JSV) =   0.0
         IF (ZT2.NE.ZT1) ZRATEINC(JSV) =  ZRATECST(JSV)/(ZT2-ZT1)
         IF (ZT4.NE.ZT3) ZRATEDEC(JSV) = -ZRATECST(JSV)/(ZT4-ZT3)
         !
         IF (KVERB >= 10) THEN
           WRITE(KLUOUT,'(A)') ' '
           WRITE(KLUOUT,'(A,E12.4)') 'Quantite totale   (g)     : ',XPPMASS(JSV)
           WRITE(KLUOUT,'(A,E12.4)') 'Debit grad montee (g.s-1) : ',ZRATEINC(JSV)
           WRITE(KLUOUT,'(A,E12.4)') 'Debit palier      (g.s-2) : ',ZRATECST(JSV)
           WRITE(KLUOUT,'(A,E12.4)') 'Debit grad desc   (g.s-2) : ',ZRATEDEC(JSV)
         END IF
         !
         !
         !*	1.5 Colonne contaminee.
         !
         !
         ! Calcul des couches auxquelles appartiennent la base et le sommet.
         !
         ZTOP    = XZHAT(IKU)
         ZGROUND = XZS(IPIGI(JSV),IPIGJ(JSV))           ! Altitude sol.
         ZCBOT  = XPPBOT(JSV)*ZTOP/(ZTOP-ZGROUND)
         ZCTOP   = XPPTOP(JSV)*ZTOP/(ZTOP-ZGROUND)
         IBOT   = COUNT( XZHAT(:)<=ZCBOT )
         ITOP    = COUNT( XZHAT(:)<=ZCTOP  )
         IF (KVERB >= 10) THEN
           WRITE(KLUOUT,'(A)') ' '
           WRITE(KLUOUT,'(A,F7.1)') 'Base      : ',XPPBOT(JSV)
           WRITE(KLUOUT,'(A,F7.1)') 'Sommet    : ',XPPTOP(JSV)
           WRITE(KLUOUT,'(A,I3.1)') 'JK Base   : ',IBOT
           WRITE(KLUOUT,'(A,I3.1)') 'JK Sommet : ',ITOP
         END IF
         !
         ZREPVERT(:,JSV)=0.0
         CNIVO(JSV) = 'ALT'
         IF ( IBOT.EQ.ITOP ) THEN
            ZREPVERT(IBOT,JSV)=1.0
            IF (IBOT.EQ.IKB) CNIVO(JSV)='SRF'
         ELSE
            ZEPAIS = ZCTOP-ZCBOT
            DO JK=IBOT,ITOP
               ZENBAS             = MAX(ZCBOT,XZHAT(JK))
               ZENHAUT            = MIN(ZCTOP, XZHAT(JK+1))
               ZREPVERT(JK,JSV) = (ZENHAUT-ZENBAS)/ZEPAIS
               IF (KVERB >= 10) &
                WRITE(KLUOUT,'(A,I3.1,A,F7.4)') 'Niveau:',JK,' Poids :',ZREPVERT(JK,JSV)
            END DO
         ENDIF
         IF (KVERB >= 10) &
          WRITE(KLUOUT,'(A,A3)') 'Nature de la source : ',CNIVO(JSV)
         !
         !
         !
      ENDIF  ! Appartenance au domaine.
   !
   !
   END DO
   !
!
   !
   !
   !
   !* 1.7 Quantite totale emise par SV.
   !
   !
   ZQTOT(:)=0.0
   DO JSV=1,NSV_PP
      ZQTOT(JSV)=ZQTOT(JSV)+XPPMASS(JSV)
   END DO
   !
   !
   !
   IF (KVERB >= 10) THEN
     WRITE(KLUOUT,'(A)') ' '
     WRITE(KLUOUT,'(A)') '**************************************'
     WRITE(KLUOUT,'(A)') ' '
     WRITE(KLUOUT,'(A)') '             EMIS_PASPOL'
     WRITE(KLUOUT,'(A)') '  FIN de l initialisation des sources  '
     WRITE(KLUOUT,'(A)') ' '
     WRITE(KLUOUT,'(A)') '***************************************'
     WRITE(KLUOUT,'(A)') ' '
   END IF
   !
ENDIF 
!
!
!
!
!
!*	2.  EMISSIONS.
!	    ----------
!
!*	2.1 Date-heure courante sous forme plus pratique.
!
WHERE (XSVT(:,:,:,NSV_PPBEG:NSV_PPEND) <0.0) &
        XSVT(:,:,:,NSV_PPBEG:NSV_PPEND)=0.0
!
DO JSV=1,NSV_PP
   !
   II=IPIGI(JSV)    
   IJ=IPIGJ(JSV)
   IP= NSV_PPBEG + JSV - 1
   !
   ZSVT(:,:,:) = XSVT(:,:,:,IP)
   IF ( (II.GE.IIB).AND.(II.LE.IIE).AND.(IJ.GE.IJB).AND.(IJ.LE.IJE) ) THEN
      !
      !
      !*	2.2 Distance temporelle DEPUIS le debut de rejet.
      !
      TZDATE%TDATE%YEAR=IDEBYY(JSV);TZDATE%TDATE%MONTH=IDEBMM(JSV);TZDATE%TDATE%DAY=IDEBDD(JSV);TZDATE%TIME=ZDEBSS(JSV)
      CALL DATETIME_DISTANCE(TZDATE,TDTCUR,ZDEPUIS)
      !
      !
      !* 2.3 Si la source emet.
      !
      IF ( (ZDEPUIS.NE.XUNDEF).AND.(ZDEPUIS.LE.ZCHRO4(JSV)) ) THEN
         !
         !
         !* 2.4 Debit en g/s.
         !
         ZRATE=ZRATECST(JSV)
         IF (ZDEPUIS.LT.ZCHRO2(JSV)) ZRATE=ZRATEINC(JSV)*ZDEPUIS
         IF (ZDEPUIS.GT.ZCHRO3(JSV)) ZRATE=ZRATECST(JSV)    &
                                          +ZRATEDEC(JSV)*(ZDEPUIS-ZCHRO3(JSV))
         !
         !	  
         IF (CNIVO(JSV).EQ.'SRF') THEN
            !
            !
            !*	2.5 Emission a la surface.
            !	    On passe par les SFSV en (g/g).m.s-1
            !
            ZP = XPABST(II,IJ,IKB)
            ZTH = XTHT(II,IJ,IKB)
            !
            ZT = ZTH *( (ZP/XP00)**(XRD/XCPD) )
            ZRHO = ZP /(XRD*ZT)*1000.0              ! en g.m-3
            ZMASAIR = ZRHO *ZSURF*(XZHAT(IKB+1)-XZHAT(IKB))   ! en g
!
            PSFSV(II,IJ,IP) = ZRATE/ZMASAIR  &
                 *(XZHAT(IKB+1)-XZHAT(IKB))
            !
            IF (.NOT.GBEGEMIS(JSV)) GBEGEMIS(JSV) = .TRUE.
            !
         ELSE
            !
            !
            !*	2.6 Emission en altitude.
            !	    On modifie directement les XVT en g/g.
            !
            DO JK=1,IKE
               !
               ZP = XPABST(II,IJ,JK)
               ZTH = XTHT(II,IJ,JK)
               !
               ZT = ZTH *((ZP/XP00)**(XRD/XCPD))
               ZRHO = ZP /(XRD*ZT)*1000.0           ! en g.m-3
               ZMASAIR = ZRHO *ZSURF*(XZHAT(JK+1)-XZHAT(JK))  ! en g
               !
               XSVT(II,IJ,JK,IP)      &
                  = XSVT(II,IJ,JK,IP) &
                  + ZRATE/ZMASAIR*XTSTEP           &
                  * ZREPVERT(JK,JSV)
               !
            END DO

            !
            !
            IF (.NOT.GBEGEMIS(JSV)) THEN
               XRSVS(:,:,:,IP) =  XRSVS(:,:,:,IP) &
                                 +XRHODJ(:,:,:)*XSVT(:,:,:,IP)/PTSTEP
               GBEGEMIS(JSV)= .TRUE.
            ELSE
               XRSVS(:,:,:,IP) = XRSVS(:,:,:,IP) &
                +XRHODJ(:,:,:)*(XSVT(:,:,:,IP)-ZSVT(:,:,:))/PTSTEP
            ENDIF
            !
            !
         ENDIF    ! TEST Nature de la source ('SRF' ou 'ALT')
      ENDIF       ! TEST Rejet actif.
   ELSE
         WRITE(KLUOUT,'(A,I3.1)')  'IJU       : ',IJU
         WRITE(KLUOUT,'(A)') 'La source n est pas dans le domaine de travail courant.'
   ENDIF          ! TEST Rejet dans le domaine de travail courant.
!
END DO            ! BOUCLE sur les rejets.
!
!
!
!
!*	3.   CALCUL DES CTA.
!       ---------------
!
!*	3.1 Calcul de la masse volumique de l'air en Kg/m3.
!
ZRHOM(:,:,:)=XPABST(:,:,:)/(XRD*XTHT(:,:,:)*((XPABST(:,:,:)/XP00)**(XRD/XCPD)))
!
!
!*	3.2 Passage en g/m3.
!
ZRHOM(:,:,:)=ZRHOM(:,:,:)*1000.0
!
!
!* 3.3 Calcul du CTA.
!
DO JSV=1,NSV_PP
   IP= NSV_PPBEG + JSV - 1
   XATC(:,:,:,JSV)=XATC(:,:,:,JSV) &
       + XSVT(:,:,:,IP)*ZRHOM(:,:,:)*XTSTEP/ZQTOT(JSV)
END DO
!
!
!*	3.4 Ecriture conditionnelle.
!
IF (OCLOSE_OUT) THEN
  ALLOCATE( ZTEMPO(IIU,IJU,IKU) )
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'm-3'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV=1,NSV_PP
    ZTEMPO(:,:,:)=XATC(:,:,:,JSV)
    !
    WRITE(TZFIELD%CMNHNAME,'(A3,I3.3)')'ATC',JSV+NSV_PPBEG-1
    TZFIELD%CLONGNAME = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A3,I3.3)')'X_Y_Z_','ATC',JSV+NSV_PPBEG-1
    !
    CALL IO_WRITE_FIELD(TPFILE,TZFIELD,ZTEMPO)
  END DO
  !
  DEALLOCATE(ZTEMPO)
ENDIF
!
DEALLOCATE(ZRHOM, ZSVT)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PASPOL
