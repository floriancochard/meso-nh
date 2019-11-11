PROGRAM COMPUTE_R00       
!     ###############################
!
! ce programme est la version PC du programme compute_r00.f90 de mesoNH utilisee
! dans DIAG pouvant tourner sur PC seul afin de pouvoir se passer du
! super-calculateur pour reconstituer a loisir des lachers de particules
! arbitraires.
!
! on garde la structure Fortran 90 et les routines d'interpolation mais on 
! saisit les noms des fichiers a travers le fichier compute_r00.nam
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!  modules commun pour la lecture et l'ecriture
!
!                    NIMAX,NJMAX,NKMAX, NIINF, NISUP
USE MODD_DIM1
!                    grille : XXDXHAT(:,1:7) et XXX(:,1:7), XXZS(:,:,1:7)
USE MODD_COORD
!                    ref grille: XLON0,XLAT0,XBETA,XRPK
USE MODD_GRID
!                    descriptif grille: XXHAT(:) ,XLAT(:,:),XDXHAT(:),XMAP(:,:)
!                    ,XZS(:,:),XZZ(:,:,:) ,XCOSSLOPE(:,:),XDIRCOSXW(:,:)
USE MODD_GRID1
!                    XVAR(i,j,k,,,), XMASK,XTRAJ ,XDATIME(16,t)
USE MODD_ALLOC_FORDIACHRO
!                    NBFILES + nom des fichiers CFILEDIAS, CLUOUTDIAS
USE MODD_DIACHRO, ONLY:CFILEDIA,CLUOUTDIA, &
                   NLUOUTDIA,NRESPDIA,NNPRARDIA,NFTYPEDIA,NVERBDIA, NNINARDIA
!
USE MODI_WRITEVAR
!
IMPLICIT NONE
!
TYPE DATE
INTEGER :: YEAR
INTEGER :: MONTH
INTEGER :: DAY
END TYPE DATE
!
TYPE DATE_TIME
TYPE (DATE) :: TDATE
REAL :: TIME
END TYPE DATE_TIME 
!
!
CHARACTER (LEN=28) :: HFMFILE   ! name of the OUTPUT FM-file
CHARACTER (LEN=31) :: HFMFILE_sto 
!
!*       0.2   declarations of local variables
!
INTEGER  :: IRESP                ! return code in FM routines
INTEGER  :: INPRAR               ! number of articles predicted  in
                                 !  the LFIFM file
INTEGER  :: ININAR               ! number of articles  present in
                                 !  the LFIFM file
INTEGER  :: ITYPE                ! type of file (conv2dia and transfer)
!
! **** la longueur du nom ne doit pas depasser 13 car. si le fichier
! contient des groupes a un seul PROCessus, ou 9 si plusieurs PROCessus ****
CHARACTER (LEN=13)                 :: YRECFM
CHARACTER (LEN=100)                :: YCOMMENT
!
INTEGER                            :: IFILECUR,JFILECUR,NIU,NJU,NKU,IGRID,ILENCH
INTEGER                            :: NFILES,JLOOP
REAL                               :: ZXOR,ZYOR,ZDX,ZDY
REAL                               :: ZSPVAL
REAL, ALLOCATABLE, DIMENSION(:,:,:):: ZX0, ZY0, ZZ0        ! origin of the 
       ! particules colocated with the mesh-grid points read in the file
REAL, ALLOCATABLE, DIMENSION(:,:,:):: ZX00, ZY00, ZZ00, ZZL ! cumulative
       ! origin for more than one restart of the tracers 
REAL, ALLOCATABLE, DIMENSION(:,:,:,:):: ZWORK 
TYPE(DATE_TIME)                    :: TDTCUR_START
CHARACTER(LEN=24)                  :: YDATE 
INTEGER                            :: IHOUR, IMINUTE
REAL                               :: ZSECOND, ZREMAIN
LOGICAL                            :: GSTART
INTEGER                            :: INBR_START
REAL                               :: ZXMAX,ZYMAX,ZZMAX  ! domain extrema
INTEGER, DIMENSION(100)            :: NBRFILES
!  declarations supplementaires
INTEGER                            :: iret     ! code de retour de lecture  
CHARACTER (LEN=3), SAVE :: CNAME_SUP
!-----------------------------------------------------------------------
!  definitions des noms de fichiers venant de modd_sto_file de Meso-NH
!  et definition de la namelist prise dans diag.f90
CHARACTER (LEN=28), SAVE :: CFILES(100)       ! names of the files to be treated
CHARACTER (LEN=28), SAVE :: CFILES_STA(100)   ! status of these files 'INIT_SV'
                                              ! if a restart of the lagrangian
                                              ! tracers has been performed
INTEGER           , SAVE :: NSTART_SUPP(100)  ! supplementary starts 
                                              ! for the lagrangian trajectories
!-----------------------------------------------------------------------
!                                              
!  article supplementaire
CHARACTER (LEN=28), SAVE :: CFIELD_LAG(100)   ! tableau de noms de record devant
CHARACTER (LEN=100),DIMENSION(:),ALLOCATABLE  :: YUNITE
! etre etudies lagrangiennement THM RVM RRM...
INTEGER                  :: NUNDEF, inbr_field, NFILES_tot, k, ifield
LOGICAL                  :: L2D
CHARACTER (LEN=3), SAVE :: CFLAGFILE
!
!-----------------------------------------------------------------------
!
NAMELIST/NAM_STO_FILE/ CFILES, NSTART_SUPP
!-----------------------------------------------------------------------
!
NAMELIST/NAM_FIELD /CFIELD_LAG
!
!-------------------------------------------------------------------------------
!
!*       1.0    Lecture des noms des fichiers et initialisation
!               -----------------------------------------------
!
! ouverture du fichier contenant les noms des fichiers diachroniques a
! traiter
!
! ecrire les fichiers dans le meme ordre que pour DIAG1.nam (cf doc Gheusi +
!   Stein) dans NAM_STO_FILES i.e. ordre inverse chrono
!
!
open (unit=104,FILE='compute_r00.nam',FORM='FORMATTED')
!
!
nverbdia=1
ITYPE=2
ZSPVAL=-1.E+11
NUNDEF=-9999
CFILES(:) = '                         '
NSTART_SUPP(:) = NUNDEF
CFILES_STA(:) = 'INIT_SV'
CFIELD_LAG(:) = '                         '
CNAME_SUP='SAM'
!
READ(104,NML=NAM_STO_FILE)
!
READ(104,NML=NAM_FIELD)
!
close(104)
!
!-------------------------------------------------------------------------------
!
!*       2.0    FIND THE FILE TO BE TREATED AND THE INIT-SV FILES
!               -------------------------------------------------
!
!
! determination du nombre de champs de var a traiter lagrangiennement
inbr_field=0
DO JLOOP=1,100
  IF (LEN_TRIM(CFIELD_LAG(JLOOP))/= 0) THEN
    inbr_field=inbr_field+1
  END IF
END DO
!
!
! recherche du nombre total de fichier a traiter
NFILES_tot=0
DO JFILECUR=1,100
  IF (LEN_TRIM(CFILES(JFILECUR)) /= 0 ) THEN
    NFILES_tot= NFILES_tot +1 
  ENDIF
END DO
!
! ouverture des fichiers
do jfilecur=1,NFILES_tot
  CFLAGFILE='OPE'
  CALL READVAR('ZSBIS',CFILES(jfilecur),CFLAGFILE,nverbdia,iret)
end do
!
if (nverbdia>0) then
  print *,'nbre de fichiers a traiter',NFILES_tot
  print *,'nombre de champs de var a traiter lagrangiennement',inbr_field
end if
!
!**************************************************
!**************************************************
! pour coller a la version MESONH je prends les memes noms
! cette boucle correspond au traitement de diag pour chacun des fichiers 
! a traiter
do ifilecur=1,NFILES_tot
HFMFILE=CFILES(ifilecur)
print *,'fichier traite HFMFILE = ',HFMFILE
!**************************************************
!**************************************************
!   rem on n'indente pas la boucle ifilecur
!pour garder le code commun avec compute_r00.f90 sur VPP
!
!
! Search the number of the files(NFILES), where the Lagrangian tracers 
!have been reinitialized 
NFILES=0
DO JFILECUR=IFILECUR+1,100
  IF (LEN_TRIM(CFILES(JFILECUR)) /= 0 .AND.        &
      CFILES_STA(JFILECUR) == 'INIT_SV') THEN
    NFILES= NFILES +1 
    NBRFILES(NFILES)=JFILECUR       ! contains the number of the files where
                                    ! the Lag. tracers have been restarted
  ENDIF
END DO
!
! compute the number of supplementary cumulative starts
INBR_START=1
DO JLOOP=1,NFILES-1
  IF (NSTART_SUPP(JLOOP)/=NUNDEF .AND. NSTART_SUPP(JLOOP)> IFILECUR ) THEN
    INBR_START=INBR_START+1
  END IF
END DO
!
if (nverbdia >0) then
  print *,'INBR_START = ',INBR_START,' pour le fichier ',IFILECUR
end if
!-------------------------------------------------------------------------------
!
!*       3.0    ALLOCATIONS OF THE ARRAYS AND CONVERSIONS
!               -----------------------------------------
!
NIU=SIZE(XZZ,1)
NJU=SIZE(XZZ,2)
NKU=SIZE(XZZ,3)
if (nju==3) then
  L2D=.TRUE.
else
  L2D=.FALSE.
end if
if (nverbdia >0) print *,'L2D = ',L2D
!
if (.NOT. allocated(ZX0)) then  ! pas d'indentation pour garder la possibilite 
                               ! de faire un diff des compute_r00
ALLOCATE(ZX0(NIU,NJU,NKU))
ALLOCATE(ZY0(NIU,NJU,NKU))
ALLOCATE(ZZ0(NIU,NJU,NKU))
ALLOCATE(ZWORK(NIU,NJU,NKU,inbr_field+3))
ALLOCATE(YUNITE(inbr_field))
ALLOCATE(ZX00(NIU,NJU,NKU))
ALLOCATE(ZY00(NIU,NJU,NKU))
ALLOCATE(ZZ00(NIU,NJU,NKU))
ALLOCATE(ZZL(NIU,NJU,NKU))
!
end if
! initial values
ZXOR=0.5 * (XXHAT(2)+XXHAT(3)) 
ZYOR=0.5 * (XYHAT(2)+XYHAT(3))
ZDX= XXHAT(3)-XXHAT(2)
ZDY= XYHAT(3)-XYHAT(2)
!ZZL=MZF(XZZ)
do k=1,nku-1
  zzl(:,:,k)=(XZZ(:,:,k)+XZZ(:,:,k+1))*0.5
end do
ZZL(:,:,NKU)=2*XZZ(:,:,NKU)-ZZL(:,:,NKU-1)
ZXMAX=ZXOR+(NIU-3)*ZDX
ZYMAX=ZYOR+(NJU-3)*ZDY
ZZMAX=ZZL(2,2,NKU-1)
!  conversion from m to km
ZXOR=ZXOR*1.E-3
ZYOR=ZYOR*1.E-3
ZDX=ZDX*1.E-3
ZDY=ZDY*1.E-3
ZZL(:,:,:)=ZZL(:,:,:)*1.E-3
ZXMAX=ZXMAX*1.E-3
ZYMAX=ZYMAX*1.E-3
ZZMAX=ZZMAX*1.E-3
!
CALL READVAR('LGXM',CFILES(ifilecur),CFLAGFILE,nverbdia,iret)
ZX00(:,:,:)=XVAR(:,:,:,1,1,1)   
CALL READVAR('LGYM',CFILES(ifilecur),CFLAGFILE,nverbdia,iret)
ZY00(:,:,:)=XVAR(:,:,:,1,1,1)  
CALL READVAR('LGZM',CFILES(ifilecur),CFLAGFILE,nverbdia,iret)
ZZ00(:,:,:)=XVAR(:,:,:,1,1,1) 
! what is the unit of Lag. var. (km after DIAG, m after MODEL) ?
IF (INDEX(CCOMMENT(1),'KM')/=0 .OR. &
    MAXVAL(ZZ00(:,:,:))<100.         ) THEN
  print*,'unit of Lagrangian variables in ',TRIM(CFILES(ifilecur)),' is KM' 
ELSE
  print*,'unit of Lagrangian variables in ',TRIM(CFILES(ifilecur)),' is M' 
  ZX00(:,:,:)=ZX00(:,:,:)*1.E-3  !  conversion from m to km
  ZY00(:,:,:)=ZY00(:,:,:)*1.E-3
  ZZ00(:,:,:)=ZZ00(:,:,:)*1.E-3
ENDIF
!
!
IF (L2D) THEN
  WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
          ZZ00>ZZMAX)
    ZX00=ZSPVAL
    ZZ00=ZSPVAL
  END WHERE
ELSE
  WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
          ZY00<ZYOR .OR. ZY00>ZYMAX .OR. &
	      ZZ00>ZZMAX)
    ZX00=ZSPVAL
    ZY00=ZSPVAL
    ZZ00=ZSPVAL
  END WHERE
END IF
!
!-------------------------------------------------------------------------------
!
!*       4.0    COMPUTE THE ORIGIN STEP BY STEP
!               -------------------------------
!
!
! General loop for the files where a reinitialisation of the tracers 
! is performed
DO JFILECUR=1,NFILES
  !
  !CALL FMOPEN_ll(CFILES(NBRFILES(JFILECUR)),'READ',CLUOUT,   &
  !               INPRAR,ITYPE,NVERB,ININAR,IRESP)
!
!*       4.1  check if this file is a start instant
!
  GSTART=.FALSE.
  DO JLOOP=1,NFILES
    IF (NBRFILES(JFILECUR)==NSTART_SUPP(JLOOP) .OR. JFILECUR==NFILES) THEN
      INBR_START=INBR_START-1
      GSTART=.TRUE.
      EXIT
    END IF
  ENDDO
  !
  if (nverbdia>0) then
   print *, 'fichier pour la reconstitution ',JFILECUR,' GSTART =',GSTART
  end if
!
!*       4.2 read the potential temp or the water vapor at the start instant      
!
  IF (GSTART) THEN
    !
    if(inbr_field>0) then
     do ifield=1,inbr_field
      YRECFM=CFIELD_LAG(ifield)
      CALL READVAR(YRECFM,CFILES(NBRFILES(JFILECUR)),CFLAGFILE &
                    ,nverbdia,iret)
      ZWORK(:,:,:,ifield)=XVAR(:,:,:,1,1,1)
      YUNITE(ifield)=CUNITE(1)
     end do
    else
      CALL READVAR('PABSM',CFILES(NBRFILES(JFILECUR)),CFLAGFILE &
                    ,nverbdia,iret)
    endif
    TDTCUR_START%TDATE%YEAR=XDATIME(5,1)
    TDTCUR_START%TDATE%MONTH=XDATIME(6,1)
    TDTCUR_START%TDATE%DAY=XDATIME(7,1)
    TDTCUR_START%TIME=XDATIME(8,1)
    IHOUR   = INT(TDTCUR_START%TIME/3600.)
    ZREMAIN = MOD(TDTCUR_START%TIME,3600.)
    IMINUTE = INT(ZREMAIN/60.)
    ZSECOND = MOD(ZREMAIN,60.)
    WRITE(YDATE,FMT='(1X,I4.4,I2.2,I2.2,2X,I2.2,"H",I2.2,"M", &
         & F5.2,"S")') TDTCUR_START%TDATE, IHOUR,IMINUTE,ZSECOND 
  END IF
!
!*       4.3  store the X0,Y0,Z0 field for the current start before 
!             computing the new origin
!
  IF (GSTART) THEN
    IGRID=1
    PRINT *,'INBR_START',INBR_START,' NBRFILES(JFILECUR)',NBRFILES(JFILECUR)
    WRITE(YRECFM,'(A2,I2.2)')'X0',INBR_START
    WRITE(YCOMMENT,'(A8,I2.2)')'X_Y_Z_X0',INBR_START
    CTITRE(1)=YRECFM
    CUNITE(1)='(KM)' 
    CCOMMENT(1)=YCOMMENT(1:10)//YDATE//' (KM)'
    PRINT *,'COMMENT = ',CCOMMENT(1)
    XVAR(:,:,:,1,1,1)=ZX00(:,:,:)
    CALL WRITEVAR(1,NIU,1,NJU,1,NKU,1,1,1,1,1,1, &
         YRECFM,HFMFILE,'OLD',CNAME_SUP,nverbdia,iret)
    !
    WRITE(YRECFM,'(A2,I2.2)')'Y0',INBR_START
    WRITE(YCOMMENT,'(A8,I2.2)')'X_Y_Z_Y0',INBR_START
    CTITRE(1)=YRECFM
    CCOMMENT(1)=YCOMMENT(1:10)//YDATE//' (KM)'
    CUNITE(1)='(KM)' 
    PRINT *,'COMMENT = ',CCOMMENT(1)
    XVAR(:,:,:,1,1,1)=ZY00(:,:,:)
    CALL WRITEVAR(1,NIU,1,NJU,1,NKU,1,1,1,1,1,1, &
         YRECFM,HFMFILE,'OLD',CNAME_SUP,nverbdia,iret)
    !
    WRITE(YRECFM,'(A2,I2.2)')'Z0',INBR_START
    WRITE(YCOMMENT,'(A8,I2.2)')'X_Y_Z_Z0',INBR_START
    CTITRE(1)=YRECFM
    CCOMMENT(1)=YCOMMENT(1:10)//YDATE//' (KM)'
    CUNITE(1)='(KM)' 
    PRINT *,'COMMENT = ',CCOMMENT(1)
    XVAR(:,:,:,1,1,1)=ZZ00(:,:,:)
    CALL WRITEVAR(1,NIU,1,NJU,1,NKU,1,1,1,1,1,1, &
         YRECFM,HFMFILE,'OLD',CNAME_SUP,nverbdia,iret)
  END IF
!
!*       4.4   compute the origin of the particules using one more segment
!
  IF (JFILECUR /= NFILES) THEN
    CALL READVAR('LGXM',CFILES(NBRFILES(JFILECUR)),  &
                  CFLAGFILE,nverbdia,iret)
    ZX0(:,:,:)=XVAR(:,:,:,1,1,1) 
    CALL READVAR('LGYM',CFILES(NBRFILES(JFILECUR)),  &
                  CFLAGFILE,nverbdia,iret)
    ZY0(:,:,:)=XVAR(:,:,:,1,1,1) 
    CALL READVAR('LGZM',CFILES(NBRFILES(JFILECUR)),  &
                  CFLAGFILE,nverbdia,iret)
    ZZ0(:,:,:)=XVAR(:,:,:,1,1,1)  
    ! what is the unit of Lag. var. (km after DIAG, m after MODEL) ?
    IF (INDEX(CCOMMENT(1),'KM')/=0 .OR. &
       MAXVAL(ZZ00(:,:,:))<100.         ) THEN
      print*,'unit of Lagrangian variables in ', &
             TRIM(CFILES(NBRFILES(jfilecur))),' is KM' 
    ELSE
      print*,'unit of Lagrangian variables in ', &
             TRIM(CFILES(NBRFILES(jfilecur))),' is M' 
      ZX00(:,:,:)=ZX00(:,:,:)*1.E-3  !  conversion from m to km
      ZY00(:,:,:)=ZY00(:,:,:)*1.E-3
      ZZ00(:,:,:)=ZZ00(:,:,:)*1.E-3
    ENDIF
    !
    ! old position of the set of particles
    ZWORK(:,:,:,inbr_field+1)=ZX00
    ZWORK(:,:,:,inbr_field+2)=ZY00
    ZWORK(:,:,:,inbr_field+3)=ZZ00
    !
    IF (L2D) THEN
      CALL INTERPXYZ(ZWORK(:,:,:,inbr_field+1),ZWORK(:,:,:,inbr_field+2),&
                     ZWORK(:,:,:,inbr_field+3),ZX0,ZX00,ZZ0,ZZ00             )
    ELSE
      CALL INTERPXYZ(ZWORK(:,:,:,inbr_field+1),ZWORK(:,:,:,inbr_field+2),&
                     ZWORK(:,:,:,inbr_field+3),ZX0,ZX00,ZY0,ZY00,ZZ0,ZZ00    )
    END IF
    !
    IF (L2D) THEN
      WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
              ZZ00>ZZMAX)
        ZX00=ZSPVAL
        ZZ00=ZSPVAL
      END WHERE
    ELSE
      WHERE ( ZX00<ZXOR .OR. ZX00>ZXMAX .OR. &
              ZY00<ZYOR .OR. ZY00>ZYMAX .OR. &
              ZZ00>ZZMAX)
        ZX00=ZSPVAL
        ZY00=ZSPVAL
        ZZ00=ZSPVAL
      END WHERE
    END IF
    !
    !
  END IF
!
!*       4.5   close the input file
!
  !!CALL FMCLOS_ll(CFILES(NBRFILES(JFILECUR)),'KEEP',CLUOUT,IRESP)
!
!
!*       4.6   compute and store potential temp and water vapor at the origin
!
  IF (GSTART) THEN
    !
    do ifield=1,inbr_field
    !
      CALL INTERPXYZ(ZX00,ZY00,ZZ00,     &
      ZWORK(:,:,:,ifield),ZWORK(:,:,:,inbr_field+1)         )
    !
    WRITE(YRECFM,'(A3,I2.2)')CFIELD_LAG(ifield),INBR_START
    CTITRE(1)=YRECFM
    print*,'CFIELD_LAG ',ifield,' TITRE= ',TRIM(CTITRE(1))
    WRITE(YCOMMENT,'(A6,A3,I2.2)')'X_Y_Z_',CFIELD_LAG(ifield),INBR_START
    CCOMMENT(1)=YCOMMENT(1:10)//YDATE//' (USI)'
    PRINT *,'COMMENT = ',TRIM(CCOMMENT(1))
    CUNITE(1)=YUNITE(ifield)
    PRINT *,'CUNIT = ',TRIM(CUNITE(1))
    XVAR(:,:,:,1,1,1)=ZWORK(:,:,:,ifield)
    CALL WRITEVAR(1,NIU,1,NJU,1,NKU,1,1,1,1,1,1, &
         YRECFM,HFMFILE,'OLD',CNAME_SUP,nverbdia,iret)
    !
    !
    end do
    !
  END IF
!
!
END DO
!
! fermeture du fichier diachronique
IF (GSTART) call WRITEVAR(1,NIU,1,NJU,1,NKU,1,1,1,1,1,1, &
                 YRECFM,HFMFILE,'CLO',CNAME_SUP,nverbdia,iret)
end do
!***********************************************
!***********************************************
!
PRINT*, ' '
PRINT*, 'COMPUTE_R00 AFTER ORIGIN COMPUTATIONS AND STORAGE'
!
!-------------------------------------------------------------------------------
!!
CONTAINS
!
!
!-------------------------------------------------------------------------------
!
!
SUBROUTINE INTERPXYZ(PX,PY,PZ,PIN1,POUT1,PIN2,POUT2,PIN3,POUT3)
!
!
!*      0. DECLARATIONS
!          ------------
!
!*       0.1  declaration of arguments
!
REAL, INTENT(IN),  DIMENSION(:,:,:)           :: PX,PY,PZ
REAL, INTENT(IN),  DIMENSION(:,:,:)           :: PIN1
REAL, INTENT(OUT), DIMENSION(:,:,:)           :: POUT1
REAL, INTENT(IN),  DIMENSION(:,:,:), OPTIONAL :: PIN2,PIN3
REAL, INTENT(OUT), DIMENSION(:,:,:), OPTIONAL :: POUT2,POUT3   
!
!*       0.2  declaration of local variables
!
INTEGER  :: JI,JJ,JK,JKK    ! loop index
INTEGER  :: II,IJ,IK        ! grid index for the interpolation
REAL     :: ZXREL,ZYREL     ! fractional grid index for the interpolation
REAL, DIMENSION(SIZE(PIN1,3)) :: ZZLXY ! vertical grid at the interpolated point
REAL     :: ZEPS1,ZEPS2,ZEPS3          ! coeff. for the interpolation
REAL     :: ZX,ZY,ZZ
LOGICAL  :: GEXT
!
!-------------------------------------------------------------------------------
!
DO JK=1,NKU
  DO JJ=1,NJU
    DO JI=1,NIU
      !
      ZX=PX(JI,JJ,JK) 
      ZY=PY(JI,JJ,JK)
      ZZ=PZ(JI,JJ,JK)
      !
      ! remove external points
      IF (L2D) THEN
        GEXT=(ZX==ZSPVAL).OR.(ZZ==ZSPVAL)
      ELSE
        GEXT=(ZX==ZSPVAL).OR.(ZY==ZSPVAL).OR.(ZZ==ZSPVAL)
      END IF
      IF (GEXT) THEN
        POUT1(JI,JJ,JK) = ZSPVAL
        IF (PRESENT(PIN2)) THEN
          POUT2(JI,JJ,JK) = ZSPVAL
        END IF
        IF (PRESENT(PIN3)) THEN
          POUT3(JI,JJ,JK) = ZSPVAL
        ENDIF
        !
        CYCLE
        !
      END IF
      !
      ZXREL=(ZX-ZXOR)/ZDX+2
      ZYREL=(ZY-ZYOR)/ZDY+2
      !
      II=FLOOR(ZXREL)
      IJ=FLOOR(ZYREL)
      !
      ZEPS1=ZXREL-REAL(II)
      ZEPS2=ZYREL-REAL(IJ)
      IF (L2D) ZEPS2=0.
      !
      DO JKK=1,NKU
        ZZLXY(JKK)=ZEPS2*(ZEPS1*(ZZL(II+1,IJ+1,JKK))+(1-ZEPS1)*(ZZL(II,IJ+1,JKK)))     &
             + (1-ZEPS2)*(ZEPS1*(ZZL(II+1,IJ,JKK))+(1-ZEPS1)*(ZZL(II,IJ,JKK)))
      ENDDO
      !
      IK=999
      DO JKK=2,NKU
        IF (ZZLXY(JKK).GE.ZZ) THEN
          IK=JKK-1
          EXIT 
        ENDIF
      ENDDO
      !
      IF (IK==999) THEN
        PRINT*,'PROBLEM AT POINT',II,IJ
        PRINT*,'XREL, YREL, Z =',ZXREL,ZYREL,ZZ
        PRINT*,'ZZLXY(NKU)',ZZLXY(NKU)
        STOP
      END IF 
      !
      ZEPS3=(ZZ-ZZLXY(IK))/(ZZLXY(IK+1)-ZZLXY(IK))
      !
      POUT1(JI,JJ,JK) =                                                       & 
        ZEPS3 *                                                               &
      (  ZEPS2*(ZEPS1*(PIN1(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PIN1(II,IJ+1,IK+1)))  &
       + (1-ZEPS2)*(ZEPS1*(PIN1(II+1,IJ,IK+1))+(1-ZEPS1)*(PIN1(II,IJ,IK+1)))  &
      )                                                                       & 
      + (1-ZEPS3) *                                                           &
      (  ZEPS2*(ZEPS1*(PIN1(II+1,IJ+1,IK))+(1-ZEPS1)*(PIN1(II,IJ+1,IK)))      &
       + (1-ZEPS2)*(ZEPS1*(PIN1(II+1,IJ,IK))+(1-ZEPS1)*(PIN1(II,IJ,IK)))      &
      )
      IF (PRESENT(POUT2)) THEN
        POUT2(JI,JJ,JK) =                                                     & 
          ZEPS3 *                                                             &
        (  ZEPS2*(ZEPS1*(PIN2(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PIN2(II,IJ+1,IK+1)))&
         + (1-ZEPS2)*(ZEPS1*(PIN2(II+1,IJ,IK+1))+(1-ZEPS1)*(PIN2(II,IJ,IK+1)))&
        )                                                                     & 
        + (1-ZEPS3) *                                                         &
        (  ZEPS2*(ZEPS1*(PIN2(II+1,IJ+1,IK))+(1-ZEPS1)*(PIN2(II,IJ+1,IK)))    &
         + (1-ZEPS2)*(ZEPS1*(PIN2(II+1,IJ,IK))+(1-ZEPS1)*(PIN2(II,IJ,IK)))    &
        )
      ENDIF
        !
      IF (PRESENT(POUT3)) THEN
        POUT3(JI,JJ,JK) =                                                     & 
          ZEPS3 *                                                             &
        (  ZEPS2*(ZEPS1*(PIN3(II+1,IJ+1,IK+1))+(1-ZEPS1)*(PIN3(II,IJ+1,IK+1)))&
         + (1-ZEPS2)*(ZEPS1*(PIN3(II+1,IJ,IK+1))+(1-ZEPS1)*(PIN3(II,IJ,IK+1)))&
        )                                                                     &
        + (1-ZEPS3) *                                                         &
        (  ZEPS2*(ZEPS1*(PIN3(II+1,IJ+1,IK))+(1-ZEPS1)*(PIN3(II,IJ+1,IK)))    &
         + (1-ZEPS2)*(ZEPS1*(PIN3(II+1,IJ,IK))+(1-ZEPS1)*(PIN3(II,IJ,IK)))    &
        )
      ENDIF
      !
    END DO
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INTERPXYZ
!
!-------------------------------------------------------------------------------
!
END program
