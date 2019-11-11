!     ######spl
      MODULE MODI_EXTRACT_AND_OPEN_FILES
!     ##################################
!
INTERFACE
!
SUBROUTINE EXTRACT_AND_OPEN_FILES(HCARIN,HCAROUT)
CHARACTER(LEN=*)    :: HCARIN
CHARACTER(LEN=*)    :: HCAROUT
END SUBROUTINE EXTRACT_AND_OPEN_FILES
!
END INTERFACE
!
END MODULE MODI_EXTRACT_AND_OPEN_FILES
!     ######spl
      SUBROUTINE EXTRACT_AND_OPEN_FILES(HCARIN,HCAROUT)
!     #################################################
!
!!****  *EXTRACT_AND_OPEN_FILES* - 
!!
!!    PURPOSE
!!    -------
!      
!
!!**  METHOD
!!    ------
!!     
!!     N.A.
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       06/06/94
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_FILES_DIACHRO ! NBGUIL
USE MODD_ALLOC_FORDIACHRO
USE MODD_RESOLVCAR
USE MODD_PARAMETERS,ONLY:JPHEXT
!USE MODD_DIM1
!USE MODN_PARA
!USE MODN_NCAR
USE MODI_CREATLINK
USE MODI_FMREAD
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------
!
CHARACTER(LEN=*)    :: HCARIN
CHARACTER(LEN=*)    :: HCAROUT
!
!*       0.1   Local variables
!              ---------------

!
CHARACTER(LEN=LEN_TRIM(HCARIN)) :: YCARIN
CHARACTER(LEN=28) :: YNAMFILE,YDUMMYFILE
CHARACTER(LEN=32) :: YDESFM   
CHARACTER(LEN=1)  :: YC1
CHARACTER(LEN=2)  :: YC2
INTEGER   ::   ILENC
INTEGER   ::   INCR, INDFI, INDQUI, IDIF, INDFIS, INDON
INTEGER   ::   ILUDES, IRESP, INUMFILECUR
INTEGER   ::   J, JJ, JM, JMM, JA, JME
INTEGER,DIMENSION(13),SAVE             :: IASF

INTEGER   ::   ISTA, IER, INB, IWK
INTEGER   ::   ILU, INUM, IRESP2
LOGICAL   ::   GPLUS
!INTEGER           :: IIINF, IJINF, IISUP, IJSUP
!REAL              :: ZIDEBCOU, ZJDEBCOU
CHARACTER(LEN=20) :: YCOMMENT
INTEGER           ::  ILENCH,ILENG,IGRID
!------------------------------------------------------------------------------
!
YCARIN = HCARIN
if(nverbia >0)then
  print *,' ENTREE EXTRACT LEN et YCARIN ',LEN(YCARIN),YCARIN
! print *,' ENTREE EXTRACT HCAROUT ',HCAROUT
endif
ILENC = LEN(YCARIN)
! En cas de superpositions ou presence _MINUS_ , on ne traite pas immediatement
INDON=INDEX(YCARIN,'_ON_')
IF(INDON == 0)THEN
  INDON=INDEX(YCARIN,'_MINUS_')
ENDIF
IF(INDON == 0)THEN
  INDON=INDEX(YCARIN,'_PLUS_')
ENDIF
IF(INDON /= 0)THEN
  HCAROUT(1:LEN(HCAROUT))=' '
  HCAROUT=YCARIN
  HCAROUT=ADJUSTL(HCAROUT)
!print *,' PRESENCE _ON_ HCAROUT ',HCAROUT
!print *,' YCARIN ',YCARIN(1:LEN_TRIM(YCARIN))
  RETURN
ENDIF
!
HCAROUT(1:LEN(HCAROUT))=' '
!print *,' HCARIN ',LEN(HCARIN)
!print *,' YCARIN ILENC ',ILENC,YCARIN
!
! Extraction des noms de fichiers
!
! Absence nom de fichier mais presence chaine _FILEx_ ou _FILExx_
!
if(nverbia >0)then
 print *,' ** EXTRACT NBGUILlemets= ',NBGUIL
endif
IF(NBGUIL == 0)THEN
  INDQUI=0
  INDQUI=INDEX(YCARIN,'_QUIT')
  IF(INDQUI == 0)THEN
    INDQUI=INDEX(YCARIN,'QUIT')
  ENDIF
  IF(INDQUI /= 0)THEN
! Fermeture des fichiers et arret du programme    
! Inutile pour les fichiers FM ouverts en lecture
    !DO J=1,NBFILES
      !CALL FMCLOS(CFILEDIAS(J),'KEEP',CLUOUTDIAS(J),NRESPDIAS(J))
      ! plante car le .des est deja ferme
    !ENDDO
    YDUMMYFILE=''
    CALL CREATLINK(' ',YDUMMYFILE,'CLEAN',NVERBIA)
    CALL FMLOOK('FICVAL','FICVAL',ILU,IRESP2)
    IF(IRESP2 == 0)THEN
      CLOSE(ILU)
    ENDIF
    CALL SFLUSH
    CALL GQOPS(ISTA)
    ! INB donne le nombre de stations ouvertes
    ! Eventuellement on ferme la WISS N9
    CALL GQOPWK(1,IER,INB,IWK)
if(nverbia >0)then
 print *,' ** EXTRACT nb de stations ouvertes INB= ',INB
endif
    IF(INB >1)THEN
      DO JJ=1,INB
        CALL GQOPWK(JJ,IER,INB,IWK)
        IF(IWK == 9)THEN
          CALL GCLWK(9)
          EXIT
        ENDIF
      ENDDO
    ENDIF
    ! INB donne le nombre de stations actives
    CALL GQACWK(1,IER,INB,IWK)
if(nverbia >0)then
 print *,' ** EXTRACT nb de stations actives INB= ',INB
endif
    IF(ISTA >1 .AND. INB > 1)THEN
      CALL GDAWK(2)
      CALL GCLWK(2)
    ENDIF
! CALL FRAME
    CALL NGPICT(1,1)
    CALL CLSGKS
if(nverbia >0)then
 print *,' ** EXTRACT AV RETURN'
endif
    RETURN
  ENDIF     ! fin de 'QUIT'
  !
  INDFI=0
  INDFI=INDEX(YCARIN,'_FILE')
  INUMFILECUR=NUMFILECUR
  IF(INDFI /= 0)THEN
    INDFIS=0
! On reutilise un fichier deja ouvert; on renvoit l'instruction sans la chaine
! _FILEx_ ou _FILExx_; on positionne le numero du fichier courant
! Cas numero suivant _FILE a 1 chiffre
    IF(YCARIN(INDFI+6:INDFI+6) == '_')THEN
      READ(YCARIN(INDFI+5:INDFI+5),'(I1)')NUMFILECUR
! Modif le 3/1/96. Pour conserver la chaine _FILEx_
!     HCAROUT(1:INDFI-1)=YCARIN(1:INDFI-1)
!     HCAROUT(INDFI:ILENC-7)=YCARIN(INDFI+7:ILENC)
      HCAROUT(1:ILENC)=YCARIN(1:ILENC)
      INDFIS=MIN(INDFI+6+1,ILENC)
! Cas numero suivant _FILE a 2 chiffres
    ELSE IF(YCARIN(INDFI+7:INDFI+7) == '_')THEN
      READ(YCARIN(INDFI+5:INDFI+6),'(I2)')NUMFILECUR
! Modif le 3/1/96. Pour conserver la chaine _FILEx_
!     HCAROUT(1:INDFI-1)=YCARIN(1:INDFI-1)
!     HCAROUT(INDFI:ILENC-8)=YCARIN(INDFI+8:ILENC)
      HCAROUT(1:ILENC)=YCARIN(1:ILENC)
      INDFIS=MIN(INDFI+7+1,ILENC)
    ENDIF
    
    JME=0
    DO JA=1,NBFILES
      IF(NUMFILES(JA) == NUMFILECUR)THEN
        JME=JA
      ENDIF
    ENDDO
    IF(JME==0) THEN
      PRINT*,'*PB avec la directive:'
      PRINT*,'  _file',NUMFILECUR,'_ n est pas associe a un nom de fichier'
      LPBREAD=.TRUE.
      RETURN
    ENDIF

!   IIINF=NIINF; IJINF=NJINF; IISUP=NISUP; IJSUP=NJSUP
!   ZIDEBCOU=XIDEBCOU; ZJDEBCOU=XJDEBCOU
!   CALL INI_CST
!   CALL READ_DIMGRIDREF(JME,CFILEDIAS(JME),CLUOUTDIAS(JME))
!   CALL INIDEF
!   NIMNMX=-1
!   LMINMAX=.TRUE.
!   CALL COMPCOORD_FORDIACHRO(0)
!   NIINF=IIINF; NJINF=IJINF; NISUP=IISUP; NJSUP=IJSUP
!   XIDEBCOU=ZIDEBCOU; XJDEBCOU=ZJDEBCOU
    IF (INUMFILECUR /= NUMFILECUR) THEN
      ! lecture de l en-tete si le fichier traite n est pas l ancien fichier
      ! courant      
      IF(NVERBIA>0) THEN
        print *,' ** EXTRACT avant lecture de l entete de ',TRIM(CFILEDIAS(JME))
      ENDIF
      CALL READ_FILEHEAD(JME,CFILEDIAS(JME),CLUOUTDIAS(JME))
    ENDIF

    INDFI=INDEX(YCARIN(INDFIS:ILENC),'_FILE')
    IF(INDFI == 0)THEN

      LFIC1=.TRUE.

    ELSE

      DO J=1,90  ! cf nb max de fic dans modd_files_diachro
        INDFI=INDEX(YCARIN(INDFIS:ILENC),'_FILE')

        IF(INDFI == 0)THEN
          EXIT

        ELSE

          LFIC1=.FALSE.
          INDFI=INDFIS+INDFI-1
          IF(J == 1)THEN
            NBSIMULT=1
            NUMFILESIMULT(:)=0
            NINDFILESIMULT(:)=0
            NUMFILESIMULT(NBSIMULT)=NUMFILECUR
          ENDIF
          NBSIMULT=NBSIMULT+1
          IF(YCARIN(INDFI+6:INDFI+6) == '_')THEN
            READ(YCARIN(INDFI+5:INDFI+5),'(I1)')NUMFILESIMULT(NBSIMULT)
            INDFIS=MIN(INDFI+6+1,ILENC)
          ELSE IF(YCARIN(INDFI+7:INDFI+7) == '_')THEN
            READ(YCARIN(INDFI+5:INDFI+6),'(I2)')NUMFILESIMULT(NBSIMULT)
            INDFIS=MIN(INDFI+7+1,ILENC)
          ENDIF

        ENDIF

      ENDDO

    ENDIF

    IF(.NOT.LFIC1)THEN
      DO J=1,NBSIMULT
        DO JA=1,NBFILES
          IF(NUMFILESIMULT(J) == NUMFILES(JA))THEN
            NINDFILESIMULT(J)=JA
            EXIT
          ENDIF
        ENDDO
        IF(NINDFILESIMULT(J)==0) THEN
          PRINT*,'*PB avec la directive:'
          PRINT*,'  _file',NUMFILECUR,'_ n est pas associe a un nom de fichier'
          LPBREAD=.TRUE.
          RETURN
        ENDIF
      ENDDO
    ENDIF

  ELSE
! Cas absence nom de fichier; on renvoit l'instruction telle quelle
    HCAROUT=ADJUSTL(YCARIN)
  ENDIF
  RETURN
ENDIF
!
! Presence d'au moins un nom de fichier
!
DO J=1,NBGUIL,2 !***********************************************************
!
  IF(YCARIN(NMGUIL(J)-1:NMGUIL(J)-1) /= '_')THEN
    print *,'*PB. UN GUILLEMET DOIT ETRE PRECEDE D UN _', &
    ' (Cas instruction _FILEx_)'
    print *,'ou ERREUR DANS LE NOM SYMBOLIQUE UTILISE. ', &
    ' VERIFIEZ LA SYNTAXE OU L''ORTHOGRAPHE DE VOS INSTRUCTIONS'
    LPBREAD=.TRUE.
    RETURN
  ENDIF
! Cas nom d'un processus
  IF(YCARIN(NMGUIL(J)-3:NMGUIL(J)-3) == '_' .AND. &
     YCARIN(NMGUIL(J)-2:NMGUIL(J)-2) == 'P')THEN
     CYCLE
  ELSE
! Cas nom d'un fichier
    INCR=1
    DO JJ=1,10
      INCR=INCR+1
      IF(YCARIN(NMGUIL(J)-INCR:NMGUIL(J)-INCR) == '_')EXIT
    ENDDO
!
! JM = indice debut chaine  _FILEx_  ou  _FILExx_
!
    JM=NMGUIL(J)-INCR;!print *,' JM ',JM
    IF(YCARIN(JM+1:JM+4) /= 'FILE')THEN
      print *,' CHAINE DE CARACTERES _FILEx_ ATTENDUE DEVANT LES GUILLEMETS', &
      '  ABSENTE. VERIFIEZ LA SYNTAXE DE VOS INSTRUCTIONS'
      STOP
    ELSE

      YNAMFILE(1:LEN(YNAMFILE))=' '
      YNAMFILE=ADJUSTL(YCARIN(NMGUIL(J)+1:NMGUIL(J+1)-1))
      IF(NVERBIA>0) THEN
        print *,' ** EXTRACT YNAMFILE ',YNAMFILE
      ENDIF

      IF(NBFILES == 0)THEN
!
! INIT GKS et ouverture du premier fichier
!
        IASF(:)=1
	CALL GQOPS(ISTA)
	IF(ISTA == 0)THEN
          CALL OPNGKS
	  CALL TABCOL_FORDIACHRO
	ENDIF
        CALL GSTXFP(-13,2)
        CALL GSASF(IASF)

	NBFILES=NBFILES+1
	CFILEDIAS(NBFILES)=ADJUSTL(YNAMFILE)
        IF (ABS(JM-NMGUIL(J))-1-1 == 4)THEN
  	  NUMFILES(NBFILES)=NBFILES
        ELSE IF (ABS(JM-NMGUIL(J))-1-1 == 5)THEN
  	  READ(YCARIN(NMGUIL(J)-2:NMGUIL(J)-2),'(I1)')NUMFILES(NBFILES)
        ELSE IF (ABS(JM-NMGUIL(J))-1-1 == 6)THEN
  	  READ(YCARIN(NMGUIL(J)-3:NMGUIL(J)-2),'(I2)')NUMFILES(NBFILES)
        ENDIF
	NUMFILECUR=NUMFILES(NBFILES)

! ouverture du listing
        CALL FMATTR(CLUOUTDIAS(NBFILES),CLUOUTDIAS(NBFILES), &
                    NLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
        OPEN(UNIT=NLUOUTDIAS(NBFILES),FILE=CLUOUTDIAS(NBFILES),FORM='FORMATTED')
        WRITE(UNIT=NLUOUTDIAS(NBFILES),FMT=1)NBFILES,' ',CFILEDIAS(NBFILES)
	1 FORMAT(' OPEN DIACHRONIC FILE ',I2.2,A,A28)

! Ouverture du fichier .lfi et fermeture du fichier .des correspondant
      IF(NVERBIA>0) THEN
        print *,' ** EXTRACT avant link et open premier fichier ',   &
                CFILEDIAS(NBFILES)
      ENDIF
        CALL CREATLINK('DIRLFI',CFILEDIAS(NBFILES),'CREAT',NVERBIA)
        CALL FMOPEN(CFILEDIAS(NBFILES),'OLD',CLUOUTDIAS(NBFILES), &
                    NNPRARDIAS(NBFILES),NFTYPEDIAS(NBFILES),NVERBDIAS(NBFILES),&
                    NNINARDIAS(NBFILES),NRESPDIAS(NBFILES))
        IF (NRESPDIAS(NBFILES) .NE. 0) THEN
          PRINT*,'*PB a l ouverture de ',CFILEDIAS(NBFILES)
          LPBREAD=.TRUE.
          RETURN
        ENDIF
	YDESFM(1:LEN(YDESFM))=' '
	YDESFM=ADJUSTL(ADJUSTR(CFILEDIAS(NBFILES))//'.des')
        CALL FMLOOK(YDESFM,YDESFM,ILUDES,IRESP)
        CLOSE(ILUDES)
        CALL FMFREE(YDESFM,CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))

! Modif le 3/1/96. Pour conserver la chaine _FILEx_
	IF(JM>=1)THEN
	  HCAROUT(1:NMGUIL(J)-1)=YCARIN(1:NMGUIL(J)-1)
        ENDIF
! READ JPHEXT
        CALL FMREAD(CFILEDIAS(NBFILES),'JPHEXT',CLUOUTDIAS(NBFILES),ILENG,JPHEXT,IGRID,ILENCH,YCOMMENT,NRESPDIAS(NBFILES))
      ELSE    ! NBFILES/=0
!
! Fichiers autres que le premier
!
	INUMFILECUR=NUMFILECUR
        NUMFILECUR=0
        DO JJ=1,NBFILES
          IF(YNAMFILE == CFILEDIAS(JJ))THEN
            PRINT*,'*PB avec la directive:'
            IF (NUMFILES(JJ)<10)THEN
              WRITE(YC1,'(I1)')NUMFILES(JJ)
              PRINT*,'  ce nom de fichier ',TRIM(YNAMFILE), &
                     ' est deja ouvert avec _FILE'//YC1,'_'
            ELSE
              WRITE(YC2,'(I2)')NUMFILES(JJ)
              PRINT*,'  ce nom de fichier ',TRIM(YNAMFILE), &
                     ' est deja ouvert avec _FILE'//YC2,'_'
            ENDIF
            LPBREAD=.TRUE.
            NUMFILECUR=INUMFILECUR
            RETURN
          END IF
        ENDDO

!       IF(INUMFILECUR /= NUMFILECUR)THEN
	IF(NUMFILECUR == 0)THEN
	  IF (ABS(JM-NMGUIL(J))-1-1 == 4)THEN       ! _file_
            ! pas d incrementation de NBFILES
	    NUMFILES(NBFILES)=NBFILES
	  ELSE IF (ABS(JM-NMGUIL(J))-1-1 == 5)THEN  ! _filex_
	    NBFILES=NBFILES+1
	    READ(YCARIN(NMGUIL(J)-2:NMGUIL(J)-2),'(I1)')NUMFILES(NBFILES)
	  ELSE IF (ABS(JM-NMGUIL(J))-1-1 == 6)THEN  ! _filexx_
	    NBFILES=NBFILES+1
	    READ(YCARIN(NMGUIL(J)-3:NMGUIL(J)-2),'(I2)')NUMFILES(NBFILES)
          ENDIF
          ! on ne passe pas dans la boucle pour _file_ car NBFILES=1
          !(sauf si _file_ et _filex_ melanges ...)
          DO JJ=1,NBFILES-1
            IF(NUMFILES(NBFILES)==NUMFILES(JJ))THEN
              PRINT*,'*PB avec la directive:'
              IF (NUMFILES(NBFILES)<10)THEN
                WRITE(YC1,'(I1)')NUMFILES(JJ)
                PRINT*,' _FILE'//YC1,'_ deja associe au ', &
                     'nom de fichier ',TRIM(CFILEDIAS(JJ))
              ELSE
                WRITE(YC2,'(I2)')NUMFILES(JJ)
                PRINT*,' _FILE'//YC2,'_ deja associe au ', &
                     'nom de fichier ',TRIM(CFILEDIAS(JJ))
              ENDIF
              NBFILES=NBFILES-1
              LPBREAD=.TRUE.
              NUMFILECUR=INUMFILECUR
              RETURN
            ENDIF
          ENDDO
          !
          NUMFILECUR=NUMFILES(NBFILES)
          CFILEDIAS(NBFILES)=ADJUSTL(YNAMFILE)

! Ouverture du fichier lfi et fermeture du fichier des correspondant
          CALL FMLOOK(CLUOUTDIAS(NBFILES),CLUOUTDIAS(NBFILES), &
                      NLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
          IF (NRESPDIAS(NBFILES) .NE. 0) THEN
            PRINT*,'*PB pour l ecriture dans ',CLUOUTDIAS(NBFILES)
            LPBREAD=.TRUE.
            RETURN
          ENDIF
          WRITE(UNIT=NLUOUTDIAS(NBFILES),FMT=1)NBFILES,' ',CFILEDIAS(NBFILES)

      IF(NVERBIA>0) THEN
        print *,' ** EXTRACT avant link et open fichier suivant'
      ENDIF
          CALL CREATLINK('DIRLFI',CFILEDIAS(NBFILES),'CREAT',NVERBIA)
          CALL FMOPEN(CFILEDIAS(NBFILES),'OLD',CLUOUTDIAS(NBFILES), &
                      NNPRARDIAS(NBFILES),NFTYPEDIAS(NBFILES),      &
                      NVERBDIAS(NBFILES),NNINARDIAS(NBFILES),NRESPDIAS(NBFILES))
          IF (NRESPDIAS(NBFILES) .NE. 0) THEN
            PRINT*,'*PB a l ouverture de ',CFILEDIAS(NBFILES)
            LPBREAD=.TRUE.
            RETURN
          ENDIF
          YDESFM(1:LEN(YDESFM))=' '
          YDESFM=ADJUSTL(ADJUSTR(CFILEDIAS(NBFILES))//'.des')
          CALL FMLOOK(YDESFM,YDESFM,ILUDES,IRESP)
          CLOSE(ILUDES)
          CALL FMFREE(YDESFM,CLUOUTDIAS(NBFILES),NRESPDIAS(NBFILES))
        ENDIF
        IF(NVERBIA>0) THEN
          print *,' ** EXTRACT fichier suivant numero: ',NUMFILECUR
        ENDIF

	IF(MAX(1,J-1) == 1)THEN
! Modif le 3/1/96. Pour conserver la chaine _FILEx_
	  IDIF=NMGUIL(J)-1-1
	  IF(IDIF >0)THEN
	    JMM=LEN_TRIM(HCAROUT)+1
! Modif le 3/1/96. Pour conserver la chaine _FILEx_
	    HCAROUT(JMM:JMM+IDIF)=YCARIN(1:NMGUIL(J)-1)
          ENDIF
	ELSE
! Modif le 3/1/96. Pour conserver la chaine _FILEx_
	  IDIF=NMGUIL(J)-1-(NMGUIL(MAX(1,J-1))+1)
          IF(IDIF >0)THEN
	    JMM=LEN_TRIM(HCAROUT)+1
! Modif le 3/1/96. Pour conserver la chaine _FILEx_
	    HCAROUT(JMM:JMM+IDIF)=YCARIN(NMGUIL(MAX(1,J-1))+1:NMGUIL(J)-1)
          ENDIF
	ENDIF

      ENDIF

      DO JA=1,NBFILES
        IF(NUMFILES(JA) == NUMFILECUR)THEN
          JME=JA
        ENDIF
      ENDDO
      IF(NVERBIA>0) THEN
        print *,' ** EXTRACT avant lecture de l entete de ',TRIM(CFILEDIAS(JME))
      ENDIF
    CALL READ_FILEHEAD(JME,CFILEDIAS(JME),CLUOUTDIAS(JME))
    LFIC1=.TRUE.

    ENDIF
      
  ENDIF
ENDDO     !***********************************************************


IDIF=ILENC-(NMGUIL(NBGUIL)+1)
!print *,' IDIF ILENC ',IDIF,ILENC,NMGUIL(NBGUIL)
IF(IDIF >0)THEN
  JMM=LEN_TRIM(HCAROUT)+1
  HCAROUT(JMM:JMM+IDIF)=YCARIN(NMGUIL(NBGUIL)+1:ILENC)
ENDIF
!
IF(nverbia >0)then
  print *,' END of EXTRACT_AND_OPEN_FILES HCAROUT ',TRIM(HCAROUT)
ENDIF
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE EXTRACT_AND_OPEN_FILES
