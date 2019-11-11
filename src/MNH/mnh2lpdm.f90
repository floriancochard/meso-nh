!MNH_LIC Copyright 2002-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
        PROGRAM MNH2LPDM
!	##############
!-----------------------------------------------------------------------------
!****	MNH2DIF COUPLAGE MESO-NH / SPRAY.
!
!	Auteur   :   Michel Bouzom, DP/SERV/ENV
!	Creation :   16.07.2002
!       Modification  : 07.01.2006 (T.LAUVAUX, adaptation LPDM)
!       Modification  : 04.01.2009 (F. BONNARDOT, DP/SER/ENV )
!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!
!-----------------------------------------------------------------------------
!
!
!
!*	0.  DECLARATIONS.
!	    -------------
!
!*	0.1 Modules.
!
USE MODD_CONF,             ONLY : CPROGRAM
USE MODD_IO_ll,            ONLY : TFILEDATA,TPTR2FILE
USE MODD_MNH2LPDM
!
USE MODE_FM,               ONLY: IO_FILE_OPEN_ll,IO_FILE_CLOSE_ll
USE MODE_IO_ll,            ONLY: INITIO_ll,SET_CONFIO_ll
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_FILE_ADD2LIST
USE MODE_MODELN_HANDLER
USE MODE_POS
!
USE MODI_MNH2LPDM_ECH
USE MODI_MNH2LPDM_INI
!
USE MODN_CONFIO
!
!
!*	0.2 Variables locales.
!
IMPLICIT NONE
!
CHARACTER(LEN=*),PARAMETER :: YFLOG = 'METEO.log'     ! Log filename
CHARACTER(LEN=*),PARAMETER :: YFNML = 'MNH2LPDM1.nam' ! Namelist filename
INTEGER,         PARAMETER :: IVERB = 5
!
INTEGER :: IFNML  ! Unit of namelist
INTEGER :: JFIC
LOGICAL :: GFOUND ! Return code when searching namelist
TYPE(TPTR2FILE),DIMENSION(JPMNHMAX) :: TZFMNH  ! MesoNH files
TYPE(TFILEDATA),POINTER :: TZDATEFILE  => NULL() ! Date file
TYPE(TFILEDATA),POINTER :: TZGRIDFILE  => NULL() ! Grid file
TYPE(TFILEDATA),POINTER :: TZMETEOFILE => NULL() ! Meteo file
TYPE(TFILEDATA),POINTER :: TZLOGFILE   => NULL() ! Log file
TYPE(TFILEDATA),POINTER :: TZNMLFILE   => NULL() ! Namelist file
!
!
!
!
!*	1.  INITIALISATION.
!	    ---------------
!
CPROGRAM='M2LPDM'
CALL GOTO_MODEL(1)
!
!*	1.1 Variables generales.
!
 CFMNH(:) = ''
!
!
!*	1.2 Initialisation routines LL.
!
CALL INITIO_ll()
!
!
!*	1.3 Ouverture du fichier log.
!
CALL IO_FILE_ADD2LIST(TZLOGFILE,YFLOG,'TXT','WRITE')
CALL IO_FILE_OPEN_ll(TZLOGFILE)
!
!
!*	1.4 Lecture des namelists.
!
CALL IO_FILE_ADD2LIST(TZNMLFILE,YFNML,'NML','READ')
CALL IO_FILE_OPEN_ll(TZNMLFILE)
IFNML = TZNMLFILE%NLU

READ(UNIT=IFNML,NML=NAM_TURB)
READ(UNIT=IFNML,NML=NAM_FIC)
print *,'Lecture de NAM_FIC OK.'

CALL POSNAM(IFNML,'NAM_CONFIO',GFOUND)
IF (GFOUND) THEN
  READ(UNIT=IFNML,NML=NAM_CONFIO)
END IF
LCDF4 = .FALSE.
LLFIOUT  = .FALSE.
LLFIREAD = .FALSE.
CALL SET_CONFIO_ll()
CALL IO_FILE_CLOSE_ll(TZNMLFILE)
!
!
!*	1.5 Comptage des FM a traiter.
!
IF (LEN_TRIM(CFMNH(1))>0) THEN
   NBMNH=1
   CALL IO_FILE_ADD2LIST(TZFMNH(1)%TZFILE,TRIM(CFMNH(1)),'UNKNOWN','READ',KLFITYPE=2,KLFIVERB=IVERB)
   DO WHILE (CFMNH(NBMNH+1).NE.'VIDE')
      NBMNH=NBMNH+1
      CALL IO_FILE_ADD2LIST(TZFMNH(NBMNH)%TZFILE,TRIM(CFMNH(NBMNH)),'UNKNOWN','READ',KLFITYPE=2,KLFIVERB=IVERB)
   END DO
   print *,NBMNH,' fichiers a traiter.'
ELSE
   STOP
ENDIF   
!
!
!
!
!*	2.  TRAITEMENTS.
!	    ------------
!
!*	2.1 Ouverture des fichiers METEO et GRILLE et DATE.
!
CALL IO_FILE_ADD2LIST(TZGRIDFILE,CFGRI,'TXT','WRITE')
CALL IO_FILE_OPEN_ll(TZGRIDFILE)
CALL IO_FILE_ADD2LIST(TZDATEFILE,CFDAT,'TXT','WRITE')
CALL IO_FILE_OPEN_ll(TZDATEFILE)
!
!
!*	2.2 Preparation du couplage.
!
CALL MNH2LPDM_INI(TZFMNH(1)%TZFILE,TZFMNH(NBMNH)%TZFILE,TZLOGFILE,TZGRIDFILE,TZDATEFILE)
!
!
!*	2.3 Traitement des echeances.
!
DO JFIC=1,NBMNH
   print*,"CFMTO(JFIC)=",CFMTO(JFIC)
   CALL IO_FILE_ADD2LIST(TZMETEOFILE,CFMTO(JFIC),'METEO','WRITE')
   CALL IO_FILE_OPEN_ll(TZMETEOFILE)
   CALL MNH2LPDM_ECH(TZFMNH(JFIC)%TZFILE,TZMETEOFILE)
   print*,"CLOSE_LL(CFMTO(JFIC)"
   CALL IO_FILE_CLOSE_ll(TZMETEOFILE)
   TZMETEOFILE => NULL()
END DO
!
!
!*	2.4 Fermeture des fichiers, METEO, GRILLE et LOG.
!
CALL IO_FILE_CLOSE_ll(TZGRIDFILE)
CALL IO_FILE_CLOSE_ll(TZDATEFILE)
CALL IO_FILE_CLOSE_ll(TZLOGFILE)
!
!
!
END PROGRAM MNH2LPDM
