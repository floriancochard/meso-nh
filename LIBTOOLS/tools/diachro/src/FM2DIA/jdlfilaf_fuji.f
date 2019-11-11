      SUBROUTINE JDLFILAF ( KREP, KNUMER, LDTOUT )
C****
C            !--------------------------------------------------!
C            !        Sous-programme du logiciel LFI            !
C            ! (Logiciel de Fichiers Indexes par nom d'article) !
C            !--------------------------------------------------!
C
C       - Version originale de LFI: Octobre 1989, auteur:
C                                   Jean CLOCHARD, METEO FRANCE.
C
C       - Aout 1991: Ajout de la notion de "facteur multiplicatif"
C         (on sait traiter un fichier dont la longueur d'article
C          "physique" est multiple de la longueur elementaire JPLARD),
C         et (sur option) toute la messagerie peut etre en anglais.
C
C       - Janvier 1996 : ajout ecriture dans 1 fichier de nom FICJD
C         du numero des enregistrements, de leur nom et de leur longueur
C         totale   (CCCCCCCCCCCCCCCCCC JDJD CCCCCCCCCCCCCCCCCCCCCCC)
C
C
C****
C        Sous-programme donnant, pour une unite logique ouverte au sens
C     du logiciel de fichiers indexes *LFI*, la Liste des Articles logi-
C     ques de donnees presents dans le Fichier, liste donnee toutefois
C     dans l'ordre PHYSIQUE ou ceux-ci figurent dans le fichier.
C        Sur option on donne aussi des renseignements sur les articles
C     (physiques) de gestion propres au logiciel, ainsi que sur les
C     trous repertories dans l'index.
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                LDTOUT (Entree) ==> Vrai si on doit donner les rensei-
C                                    gnements optionnels (qui ne concer-
C                                    nent pas directement les articles
C                                    logiques de donnees).
C
C
#include "lficom0.h"
C
C
C----- DESCRIPTION DES "PARAMETER" DU LOGICIEL DE FICHIERS INDEXES -----
C
C     JPDBLE= PRECISION UTILISE POUR LES ENTIERS
C                 * SI JPDBLE=8 COMPILER EN INTEGER 32 BITS
C                 * SI JPDBLE=4 COMPILER EN INTEGER 64 BITS
C
      INTEGER JPDBLE
C
      PARAMETER (JPDBLE=8)
C
C--- DESCRIPTIF DES TABLES CONCERNANT LES (PAIRES DE) PAGES D'INDEX ----
C                       ( ALIAS "P.P.I." )
C
C     CNOMAR = TABLE DES PAGES D'INDEX DE TYPE "NOMS D'ARTICLES"
C     MLGPOS = TABLE DES PAGES D'INDEX DE TYPE "LONGUEUR/POSITION"
C     MRGPIF = TABLE DES RANGS DES P.P.I. DANS LEUR FICHIER RESPECTIF
C     MCOPIF = TABLE DE CORRESPONDANCE PAGES D'INDEX/UNITES LOGIQUES
C     MRGPIM = TABLE DES RANGS EN MEMOIRE DES P.P.I. AFFECTEES
C              ( DANS *MCOPIF,MRGPIF,CNOMAR,MLGPOS,LECRPI,LPHASP* )
C     LECRPI = VRAI SI LA PAGE D'INDEX CORRESP. DOIT ETRE (RE)ECRITE
C              (.,1) ==> PAGE "NOM", (.,2) ==> PAGE "LONGUEUR/POSITION"
C     LPHASP = VRAI SI LA PAGE D'INDEX "LONG/POS" EST PHASEE EN MEMOIRE
C              AVEC LA PAGE D'INDEX "NOM" CORRESPONDANTE
C
C---------------- VARIABLES "SIMPLES" GLOBALES -------------------------
C
C     NBFIOU = Nombre d'Unites Logiques ouvertes
C     NFACTM = Somme des Facteurs Multiplicatifs utilises
C     NIMESG = NIVEAU *GLOBAL* DE LA MESSAGERIE
C     NERFAG = NIVEAU DE FILTRAGE GLOBAL DES ERREURS FATALES
C     NISTAG = NIVEAU D'IMPRESSION GLOBAL DES STATISTIQUES
C     NPISAF = NBRE DE PAIRES DE PAGES D'INDEX SUPPLEMENTAIRES AFFECTEES
C     LMULTI = VRAI SI ON DOIT TRAVAILLER EN MODE MULTI-TACHES
C     LTAMLG = OPTION PAR DEFAUT D'UTILISATION DE LA MEMOIRE TAMPON EN
C              LECTURE; VRAIE ==> UTILISATION MAXIMUM
C     LTAMEG = CF. CI-DESSUS, EN ECRITURE
C     VERGLA = VERROU GLOBAL (EN MULTI-TASKING)
C     NULOFM = Nombre d'Unites LOgiques a Facteur Multiplicat. predefini
C     CHINCO = Nom par defaut d'une variable qui devrait etre CHaracter
C     NUIMEX = Nombre d'Unites LOgiques en cours d'IMport/EXport
C
C--------- DESCRIPTIF DES ELEMENTS CONCERNANT UNE UNITE LOGIQUE --------
C
C     NUMIND = TABLE D'ADRESSAGE INDIRECT DANS LES TABLEAUX CI-DESSOUS
C     NUMERO = NUMERO DE L'UNITE LOGIQUE
C     MFACTM = FACteur Multiplicatif de la longueur physique elementaire
C     CNOMFI = NOM eventuel du FIchier associe a l'unite logique
C     CNOMSY = Idem pour le systeme, ou a defaut pour l'utilisateur.
C     NLNOMF = LONGUEUR (CARACTERES) DU NOM EVENTUEL
C     NLNOMS = Longueur (en caracteres) du Nom SYSTEME eventuel
C     NDEROP = CODE DE LA DERNIERE ACTION EFFECTUEE
C     CSTAOP = 'STATUS' DE L'OUVERTURE
C     LNOUFI = VRAI SI LE FICHIER EST NOUVEAU (AU SENS DU LOGICIEL)
C     LMODIF =  "   "   "    "    A ETE MODIFIE DEPUIS L'OUVERTURE
C     NDERCO = DERNIER CODE-REPONSE (CORRESPONDANT A LA DERNIERE ACTION)
C     MTAMPD = PAGES DE DONNEES "TAMPON"
C     NUMAPD = NUMERO D'ARTICLE PHYSIQUE CORRESPONDANT A CES PAGES
C     LECRPD = VRAI SI LA PAGE DE DONNEES CORRESP. DOIT ETRE ECRITE
C     NLONPD = LONGUEUR DE PAGE DE DONNEES REELLEMENT REMPLIE
C     NDERPD = NUMERO DE LA DERNIERE PAGE DE DONNEES UTILISEE
C     NPODPI = RANG DE LA DERNIERE PAGE D'INDEX DANS LA TABLE *MRGPIM*
C     NALDPI = NOMBRE D'ARTICLES LOGIQUES DANS LA DERNIERE PAGE D'INDEX
C     NBLECT =    "   DE LECTURES          EFFECTUEES DEPUIS L'OUVERTURE
C     NBNECR =    "   "  NOUVELLES ECRITURES    "        "       "
C     NREESP =    "   "  "VRAIES" REECRITURES SUR PLACE  "       "
C     NREECO =    "   "  REECRITURES PLUS COURTES        "       "
C     NREELO =    "   "       "      PLUS LONGUES        "       "
C     NBRENO =    "   "  FOIS OU ON A RENOMME UN ARTICLE "       "
C     NBSUPP =    "   "   "  " "  " " SUPPRIME "    "    "       "
C     NBTROU =    "   "  TROUS D'INDEX CREES             "       "
C     NIVMES = NIVEAU DE LA MESSAGERIE
C     LERFAT = VRAI SI TOUTE ERREUR DOIT ETRE FATALE
C     LISTAT = OPTION D'IMPRESSION DES STATISTIQUES ( A LA FERMETURE )
C     VERRUE = VERROU DE L'UNITE LOGIQUE (EN MODE MULTI-TASKING)
C     NPPIMM = NBRE DE PAIRES DE PAGES D'INDEX EN MEMOIRE
C     MDES1D = TABLE CONTENANT LE 1ER ARTICLE ("DESCRIPTIF")
C     NTRULZ = NOMBRE DE TROUS D'INDEX DE LONGUEUR NULLE
C     NRFPTZ = RANG PREMIERE ARTICLE AYANT LA CARACTERISTIQUE CI-DESSUS
C     NRFDTZ =   "  DERNIER     "    "    "         "         "
C     NBREAD = NOMBRE DE "READ" FORTRAN REELLEMENT EXECUTES  (DEPUIS L'
C     NBWRIT =    "      "WRITE"   "        "         "       OUVERTURE)
C     NBMOLU = NOMBRE DE MOTS UTILISATEUR LUS   CORRECTEMENT (DEPUIS L'
C     NBMOEC =    "    "   "       "      ECRITS     "        OUVERTURE)
C     LTAMPL = OPTION D'UTILISATION MAXI DE LA MEMOIRE TAMPON EN LECTURE
C     LTAMPE =    "   "      "       "   "   "    "      "    " ECRITURE
C     NDERGF = RANG DANS LE FICHIER DU DERNIER ARTICLE LOGIQUE LU
C              ou dont on a demande les caracteristiques (LFICAS/LFICAP)
C     CNDERA = NOM de ce dernier article logique de donnees
C     NSUIVF = RANG DANS LE FICHIER DU PROCHAIN ARTICLE LOGIQUE A LIRE
C              "SEQUENTIELLEMENT"
C     NPRECF = RANG DANS LE FICHIER DU PROCHAIN ARTICLE LOGIQUE
C              "PRECEDENT" A LIRE
C     LMIMAL = VRAI SI ON DOIT RECALCULER LES LONGUEURS MINI. ET MAXI.
C              DES ARTICLES LOGIQUES DE DONNEES
C     NUMAPH = NUMero d'Article PHysique (pour messages d'erreur E/S).
C     NEXPOR = Rang eventuel (d'EXPORt) dans les tables MNUIEX,NDIMPL,
C     NIMPOR =  "      "     (d'IMPORt) NDEXPL,NREXPL,CNEXPL,NIMPEX...
C
C------------------------ VARIABLES DIVERSES ---------------------------
C
C     MULOFM = Table des Unites LOgiques avec Facteur Multip. predefini
C     MFACTU =   "    "  FActeurs mUltiplicatifs associes a ces Unites
C     MNUIEX =   "    "  Numeros d'Unites logiques en Import/EXport
C     NINIEX =   "   d'adressage INdirect dans MNUIEX
C     NDIMPL = Descripteurs IMPLicites d'import/export en memoire
C     NDEXPL =      "       EXPLicites "   "   /  "    "     "
C     CNIMPL = Profil des articles a description IMPLicite
C     NAEXPL = Nombre d'articles decrits EXPLicitement
C     CNEXPL = Noms des articles decrits dans NDEXPL
C     NREXPL = Rang  "      "       "      "  NDEXPL
C     NIMPEX = Numero d'unite logique associee a l'IMPort ou l'EXport.
C     NUTRAV =    "   "   "      "    de TRAVail pour import ou export.
C     NLAPFD = Longueur d'Article Physique du fichier d'export/import.
C     NXCNLD = Nb.maX. Caracteres/Nom d'article du logiciel LFI Distant.
C     NRCFMX = Rang de la config. Imp/eXport dans CFGMXD, NBMOSD, NBCASD
C     CFGMXD = ConFiGuration pour iMport/eXport des systemes Distants.
C     NBMOSD = Nombre de Bits par MOt       des systemes Distants.
C     NBCASD =    "   "    "   "  CAractere  "     "        "    .
C     CTYPMX = Liste des types de variables valides pour Import/eXport.
C
      CHARACTER*(JPNCPN) CNOMAR (JPNXNA*JPNXPI), CNDERA (JPNXFI), CHINCO
      CHARACTER*(JPLFTX) CNOMFI (JPNXFI), CNOMSY (JPNXFI), CLACTI
      CHARACTER CSTAOP (JPNXFI)*(JPLSTX), CLNSPR*(JPLSPX), CLMESS*132
      CHARACTER CNEXPL (JPXDAM,JPIMEX)*(JPNCPN), CTYPMX*(JPTYMX)
      CHARACTER CNIMPL (JPIMEX)*(JPXMET), CFGMXD (0:JPCFMX)*(JPXCCF)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CHARACTER*16 CFICJD,CFICJDOUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      COMMON /LFICHA/ CNOMAR, CNDERA, CNOMFI, CNOMSY, CSTAOP, CHINCO
     S              , CNEXPL, CNIMPL, CFGMXD, CTYPMX
C
      INTEGER NBFIOU, NFACTM, NIMESG, NERFAG, NISTAG, NPISAF, NULOFM
      INTEGER (KIND=JPDBLE) MLGPOS (JPLARD*JPNXPI)
      INTEGER (KIND=JPDBLE) MTAMPD (JPLARD*JPNPDF*JPNXFI)
      INTEGER (KIND=JPDBLE) MDES1D (JPLARD*JPNXFI)
      INTEGER MRGPIM (JPNPIA+JPNPIS,JPNXFI), NDERPD (JPNXFI)
      INTEGER MCOPIF (JPNXPI), MRGPIF (JPNXPI), NLNOMS (JPNXFI)
      INTEGER NUMERO (JPNXFI), NLNOMF (JPNXFI), NDERCO (JPNXFI)
      INTEGER NPODPI (JPNXFI), NUMAPH (0:JPNXFI)
      INTEGER NALDPI (JPNXFI), NBLECT (JPNXFI), NBNECR (JPNXFI)
      INTEGER NREESP (JPNXFI), NREECO (JPNXFI), NREELO (JPNXFI)
      INTEGER NIVMES (0:JPNXFI), NDEROP (JPNXFI), NPPIMM (JPNXFI)
      INTEGER NUMAPD (0:JPNPDF-1,JPNXFI), NLONPD (0:JPNPDF-1,JPNXFI)
      INTEGER NTRULZ (JPNXFI), NRFPTZ (JPNXFI), NRFDTZ (JPNXFI)
      INTEGER NBTROU (JPNXFI), NUMIND (JPNXFI), NBREAD (JPNXFI)
      INTEGER NBWRIT (JPNXFI), NBMOLU (JPNXFI), NBMOEC (JPNXFI)
      INTEGER NDERGF (JPNXFI), NSUIVF (JPNXFI), NPRECF (JPNXFI)
      INTEGER NBRENO (JPNXFI), NBSUPP (JPNXFI), MFACTM (0:JPNXFI)
      INTEGER MULOFM (JPXUFM), MFACTU (0:JPXUFM)
      INTEGER NIMPEX (JPIMEX), NUTRAV (JPIMEX), NBMOSD (0:JPCFMX)
      INTEGER NBCASD (0:JPCFMX), NLAPFD (JPIMEX)
      INTEGER MNUIEX (JPIMEX), NINIEX (JPIMEX), NDEXPL (JPDEXP,JPIMEX)
      INTEGER NDIMPL (JPDIMP,JPIMEX), NXCNLD (JPIMEX), NAEXPL (JPIMEX)
      INTEGER NEXPOR (JPNXFI), NIMPOR (JPNXFI), NUIMEX, NRCFMX (JPIMEX)
      INTEGER NREXPL (0:JPXDAM,JPIMEX)
C
      REAL VERRUE (JPNXFI), VERGLA
C
      LOGICAL LLFATA, LMULTI, LTAMLG, LTAMEG, LECRPI (JPNXPI,2)
      LOGICAL LTAMPL (JPNXFI), LTAMPE (JPNXFI), LMODIF (JPNXFI)
      LOGICAL LNOUFI (JPNXFI), LERFAT (0:JPNXFI), LISTAT (JPNXFI)
      LOGICAL LPHASP (JPNXPI), LECRPD (0:JPNPDF-1,JPNXFI)
      LOGICAL LMIMAL (JPNXFI)
C
      COMMON /LFIDIV/ NBFIOU, NIMESG, NERFAG, NISTAG, NPISAF, LMULTI
     S              , VERGLA, LTAMLG, LTAMEG, MRGPIM, MRGPIF, NUMIND
     S              , VERRUE, MLGPOS, MDES1D, MCOPIF, LECRPI, LPHASP
     S              , NUMERO, NLNOMF, LNOUFI, NDERCO, MTAMPD, NUMAPD
     S              , NPODPI, NALDPI, NBLECT, NBNECR, NREESP, NREECO
     S              , NREELO, NIVMES, LERFAT, LISTAT, NDEROP, LMODIF
     S              , NPPIMM, NRFPTZ, NRFDTZ, NTRULZ, NBREAD, NBWRIT
     S              , LECRPD, NLONPD, NDERPD, NBTROU, NBMOLU, NBMOEC
     S              , LTAMPL, LTAMPE, NDERGF, NSUIVF, NBRENO, NBSUPP
     S              , LMIMAL, NPRECF, MFACTM, NULOFM, MULOFM, MFACTU
     S              , NLNOMS, NFACTM, NUMAPH, NEXPOR, NIMPOR, NIMPEX
     S              , NUTRAV, NBMOSD, NBCASD, NLAPFD, NXCNLD, NUIMEX
     S              , MNUIEX, NINIEX, NDEXPL, NREXPL, NDIMPL, NAEXPL
     S              , NRCFMX
C
C
      INTEGER KREP, KNUMER, IMDESC, IREP, IRANG, INTROU, INBPIR, INBALO
      INTEGER INALDO, IFACTM, ILARPH, INALPP, INTPPI, INPPIM, INIMES, J
      INTEGER INAGES, IRESER, INUTIL, IPERTE, IPOSFI, IPOSDE, INEXCE
      INTEGER INABAL, INALDI, INTROI, INPIMD, INPIMF, INPILE, JRGPIF
      INTEGER IRGPFS, IRGPIM, IRANGM, IRPIMS, INALPI, ILONGA, IRECPI
      INTEGER IDERPU, IREC, IRETIN
C
      LOGICAL LDTOUT
C
C
C       FONCTION SERVANT A RENDRE FATALE OU NON UNE ERREUR DETECTEE,
C       A L'AIDE DU CODE-REPONSE COURANT, DU NIVEAU DE FILTRAGE GLOBAL,
C       ET DE L'OPTION D'ERREUR FATALE PROPRE AU FICHIER.
C       S'IL N'Y A PAS DE FICHIER (I5678=0, D'OU DIMENSIONNEMENT DE
C          *LERFAT*), LE NIVEAU DE FILTRAGE JOUE LE ROLE PRINCIPAL.
C
      INTEGER IXNIMS, I1234, I5678, I3456, IXC, IXM, IXT, IABCDE, IFGHIJ
      INTEGER IKLMNO, IPQRST, IUVWXY, IZABCD, IEFGHI
C
      LOGICAL LLMOER
C
      LLMOER (I1234,I5678)=I1234.EQ.-16.OR.
     S (I1234.NE.0.AND.(NERFAG.EQ.0.OR.(NERFAG.EQ.1.AND.LERFAT(I5678))))
C
C       FONCTION DONNANT LE PLUS HAUT NIVEAU DE MESSAGERIE ACCEPTABLE
C       POUR L'UNITE LOGIQUE DE RANG "I3456" .
C       (UTILISATION DES NIVEAUX DE MESSAGERIE GLOBAL ET PROPRE AU
C        FICHIER - MEME REMARQUE QUE CI-DESSUS SI I3456=0, POUR NIVMES)
C
      IXNIMS (I3456)=MIN0 (2,2*NIMESG,MAX0 (2*NIMESG-2,NIVMES(I3456)))
C
C       Fonctions servant a l'adressage 1D dans les tableaux CNOMAR,
C     MLGPOS et MDES1D, MTAMPD.
C
      IXC (IABCDE,IFGHIJ) = IABCDE + JPNXNA * ( IFGHIJ - 1 )
      IXM (IKLMNO,IPQRST) = IKLMNO + JPLARD * ( IPQRST - 1 )
      IXT (IUVWXY,IZABCD,IEFGHI) = IUVWXY + JPLARD *
     S ( MFACTM(IEFGHI) * IZABCD + JPNPDF * ( IEFGHI - 1 ) )
C
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, PUIS INITIALISATIONS.
C-----------------------------------------------------------------------
C
      IREP=0
      IRANG=0
      CLNSPR='LFILAF'
      print *,' jdlfilaf BALISE 1 KNUMER,IRANG',KNUMER,IRANG
      CALL LFINUM (KNUMER,IRANG)
      print *,' jdlfilaf BALISE 1 Bis LMULTI',LMULTI,KNUMER,IRANG
C
      IF (IRANG.EQ.0) THEN
        IREP=-1
        GOTO 1001
      ENDIF
C
      IF (LMULTI) CALL LFIVER (VERRUE(IRANG),'ON')
      INTROU=MDES1D(IXM(JPNTRU,IRANG))+NBTROU(IRANG)
      INBPIR=MDES1D(IXM(JPNPIR,IRANG))
      INBALO=MDES1D(IXM(JPNALO,IRANG))
      INALDO=INBALO-INTROU
      print *,' MFACTM(0), MFACTM(1) ',MFACTM(0),MFACTM(1)
      IFACTM=MFACTM(IRANG)
      print *,' jdlfilaf BALISE 1 IRANG, IFACTM ',IRANG,IFACTM
      ILARPH=JPLARD*IFACTM
      INALPP=JPNAPP*IFACTM
C     INALPP=512
      print *,' jdlfilaf BALISE 1 INALPP',INALPP
C     INALPP=1
      INTPPI=(INBALO-1+INALPP)/INALPP
      INPPIM=NPPIMM(IRANG)
C
C         Envoi d'une banniere.
C
      WRITE (UNIT=*,FMT='(///)')
C
      IF (LFRANC) THEN
        WRITE (UNIT=CLMESS,FMT='(''Catalogue de l''''Unite Logique LFI''
     S ,I3,'' dans l''''ordre *PHYSIQUE* (sequentiel) des articles'')')
     S     KNUMER
      ELSE
        WRITE (UNIT=CLMESS,FMT='(''Catalog of LFI Logical Unit'',I3,
     S         '' in *PHYSICAL* (sequential) record order'')') KNUMER
      ENDIF
C
      INIMES=2
      LLFATA=.FALSE.
      CALL LFIEMS (KNUMER,INIMES,IREP,LLFATA,CLMESS,CLNSPR,CLACTI)
C**
C     2.  -  SUR OPTION, RENSEIGNEMENTS SUR LES ARTICLES "DE GESTION".
C            (ARTICLE DOCUMENTAIRE, PAIRES D'ARTICLES D'INDEX)
C-----------------------------------------------------------------------
C
      print *,' jdlfilaf BALISE 2'
      IF (LDTOUT) THEN
        INAGES=1+2*INBPIR
        IRESER=ILARPH*INAGES
C
        IF (LFRANC) THEN
          WRITE (UNIT=*,FMT='(//,TR1,I6,
     S           '' article(s) "physique(s)" de gestion,'',I6,
     S           '' mots chacun, occupant donc'',I7,'' mots; detail:'',
     S /,TR10,''Article documentaire de la position 1 a'',I6,/,TR10,I6,
     S'' paire(s) d''''articles d''''index prereserves, de la position''
     S           ,I6,'' a'',I7)')
     S         INAGES,ILARPH,IRESER,ILARPH,INBPIR,ILARPH+1,IRESER
        ELSE
          WRITE (UNIT=*,FMT='(//,TR1,I6,
     S           '' "physical" records for file handling,'',I6,
     S           '' words each, occupying then'',I7,'' words; detail:'',
     S /,TR10,''Documentary record from position 1 to'',I6,/,TR10,I6,
     S'' pair(s) of pre-reserved index records, from position''
     S           ,I6,'' to'',I7)')
     S         INAGES,ILARPH,IRESER,ILARPH,INBPIR,ILARPH+1,IRESER
        ENDIF
C
        IF (INTPPI.LT.INBPIR) THEN
          INUTIL=INBPIR-INTPPI
          IPERTE=ILARPH*INUTIL*2
C
          IF (LFRANC) THEN
            WRITE (UNIT=*,FMT='(/,TR10,5(''=''),''> Il y a'',I3,
     S '' paire(s) d''''articles d''''index inutilises, representant'',
     S             I8,'' mots'')') INUTIL,IPERTE
          ELSE
            WRITE (UNIT=*,FMT='(/,TR10,5(''=''),''> There is (are)'',I3,
     S '' pair(s) of unused index records, leading to a loss of'',
     S             I8,'' words'')') INUTIL,IPERTE
          ENDIF
C
        ELSEIF (INTPPI.EQ.INBPIR) THEN
C
          IF (LFRANC) THEN
            WRITE (UNIT=*,FMT='(TR15,5(''-''),TR3,''pas de paire '',
     S        ''d''''articles d''''index inutilises ni excedentaires'',
     S          TR3,5(''-''))')
          ELSE
            WRITE (UNIT=*,FMT='(TR15,5(''-''),TR3,''no pair of '',
     S        ''unused or overflow pages'',
     S          TR3,5(''-''))')
          ENDIF
C
        ELSEIF (INTPPI.EQ.(INBPIR+1)) THEN
          IPOSFI=ILARPH*(MDES1D(IXM(ILARPH,IRANG))+1)
          IPOSDE=IPOSFI-2*ILARPH+1
C
          IF (LFRANC) THEN
            WRITE (UNIT=*,FMT='(TR10,''une paire d''''articles '',
     S             ''d''''index excedentaires, de la position'',
     S             I9,'' a'',I9)')
     S      IPOSDE,IPOSFI
          ELSE
            WRITE (UNIT=*,FMT='(TR10,''one pair of overflow index '',
     S             ''pages ,from position'',
     S             I9,'' to'',I9)')
     S      IPOSDE,IPOSFI
          ENDIF
C
      print *,' jdlfilaf BALISE 3'
        ELSE
          INEXCE=INTPPI-INBPIR
C
          IF (LFRANC) THEN
            WRITE (UNIT=*,FMT='(TR10,I6,'' paires d''''articles '',
     S           ''d''''index excedentaires, des positions:'')') INEXCE
C
            DO 201 J=1,INEXCE
            IPOSFI=ILARPH*(MDES1D(IXM(ILARPH+1-J,IRANG))+1)
            IPOSDE=IPOSFI-2*ILARPH+1
            WRITE (UNIT=*,FMT='(TR20,I9,'' a'',I9)') IPOSDE,IPOSFI
  201       CONTINUE
C
          ELSE
            WRITE (UNIT=*,FMT='(TR10,I6,'' pairs of overflow index '',
     S           ''pages, from positions:'')') INEXCE
C
            DO 202 J=1,INEXCE
            IPOSFI=ILARPH*(MDES1D(IXM(ILARPH+1-J,IRANG))+1)
            IPOSDE=IPOSFI-2*ILARPH+1
            WRITE (UNIT=*,FMT='(TR20,I9,'' to'',I9)') IPOSDE,IPOSFI
  202       CONTINUE
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
      WRITE (UNIT=*,FMT='(//)')
C**
C     3.  -  RENSEIGNEMENTS INDIVIDUALISES SUR LES ARTICLES LOGIQUES.
C            (DONNEES, ET SUR OPTION TROUS REPERTORIES DANS L'INDEX)
C-----------------------------------------------------------------------
      print *,' jdlfilaf BALISE 4'
C
      IF (LFRANC) THEN
C
        IF (INBALO.EQ.0) THEN
          WRITE (UNIT=*,FMT='(/,TR10,5(''=''),''> L''''unite logique'',
     S I3,'' ne contient AUCUN ARTICLE LOGIQUE (ni donnees, ni trous)'',
     S           //)') KNUMER
          GOTO 1001
        ELSEIF (INBALO.EQ.INTROU) THEN
          WRITE (UNIT=*,FMT='(/,TR10,5(''=''),''> L''''unite logique'',
     S I3,'' ne contient QUE DES TROUS, pas de donnees)'',//)') KNUMER
          IF (.NOT.LDTOUT) GOTO 1001
        ENDIF
C
      ELSE
C
        IF (INBALO.EQ.0) THEN
          WRITE (UNIT=*,FMT='(/,TR10,5(''=''),''> The logical unit'',I3,
     S '' contains NO LOGICAL RECORD AT ALL (neither data, nor holes)'',
     S           //)') KNUMER
          GOTO 1001
        ELSEIF (INBALO.EQ.INTROU) THEN
          WRITE (UNIT=*,FMT='(/,TR10,5(''=''),''> The logical unit'',I3,
     S '' contains ONLY HOLES, no dat)'',//)') KNUMER
          IF (.NOT.LDTOUT) GOTO 1001
        ENDIF
C
      ENDIF
C*
C     3.1 -  BALAYAGE DES PAIRES D'ARTICLES D'INDEX, PAR ORDRE CROISSANT
C-----------------------------------------------------------------------
C
      INABAL=0
      INALDI=0
      INTROI=0
      INPIMD=2
      INPIMF=INPPIM
      IF (NPODPI(IRANG).EQ.2) INPIMD=3
      IF (NPODPI(IRANG).EQ.INPPIM) INPIMF=INPPIM-1
      INPILE=2
C
      DO 319 JRGPIF=1,INTPPI
      IRGPFS=JRGPIF+1
C
C        On fait en sorte que la P.A.I. concernee, ainsi que sa suivante
C     eventuelle, soient toutes les deux en memoire.
C
      IF (JRGPIF.EQ.INTPPI) THEN
        IRGPIM=MRGPIM(NPODPI(IRANG),IRANG)
        GOTO 314
C
      ELSEIF (JRGPIF.NE.1) THEN
C
C       Recherche de la P.A.I. dans les Paires de Pages d'Index memoire.
C
        DO 311 J=INPIMD,INPIMF
        IRGPIM=MRGPIM(J,IRANG)
C
        IF (MRGPIF(IRGPIM).EQ.JRGPIF) THEN
C
          IF (.NOT.LPHASP(IRGPIM)) THEN
C
            CALL LFIPHA (IREP,IRANG,IRGPIM,IRETIN)
C
            IF (IRETIN.EQ.1) THEN
              GOTO 903
            ELSEIF (IRETIN.EQ.2) THEN
              GOTO 904
            ELSEIF (IRETIN.NE.0) THEN
              GOTO 1001
            ENDIF
C
          ENDIF
C
          GOTO 312
C
        ENDIF
C
      print *,' jdlfilaf BALISE 5'
  311   CONTINUE
C
C          Mise en memoire de la Paire d'Articles d'Index cherchee.
C
        CALL LFIPIM (IREP,IRANG,IRANGM,IRGPIM,JRGPIF,IRGPFS,INPILE,
     S               IRETIN)
C
        IF (IRETIN.EQ.1) THEN
          GOTO 903
        ELSEIF (IRETIN.EQ.2) THEN
          GOTO 904
        ELSEIF (IRETIN.NE.0) THEN
          GOTO 1001
        ELSEIF (IRANGM.GT.INPPIM) THEN
          INPPIM=IRANGM
          INPIMF=INPPIM
        ENDIF
C
      ELSE
        IRGPIM=MRGPIM(1,IRANG)
C
      ENDIF
C
  312 CONTINUE
C
      IF (IRGPFS.EQ.INTPPI) THEN
        IRPIMS=MRGPIM(NPODPI(IRANG),IRANG)
C
      ELSE
C
C       Recherche de la P.A.I. dans les Paires de Pages d'Index memoire.
C
        DO 313 J=INPIMD,INPIMF
        IRPIMS=MRGPIM(J,IRANG)
C
        IF (MRGPIF(IRPIMS).EQ.IRGPFS) THEN
C
          IF (.NOT.LPHASP(IRPIMS)) THEN
C
            CALL LFIPHA (IREP,IRANG,IRPIMS,IRETIN)
C
            IF (IRETIN.EQ.1) THEN
              GOTO 903
            ELSEIF (IRETIN.EQ.2) THEN
              GOTO 904
            ELSEIF (IRETIN.NE.0) THEN
              GOTO 1001
            ENDIF
C
          ENDIF
C
          GOTO 314
C
        ENDIF
C
  313   CONTINUE
C
C          Mise en memoire de la Paire d'Articles d'Index cherchee.
C
      print *,' jdlfilaf BALISE 6'
        CALL LFIPIM (IREP,IRANG,IRANGM,IRPIMS,IRGPFS,JRGPIF,INPILE,
     S               IRETIN)
C
        IF (IRETIN.EQ.1) THEN
          GOTO 903
        ELSEIF (IRETIN.EQ.2) THEN
          GOTO 904
        ELSEIF (IRETIN.NE.0) THEN
          GOTO 1001
        ELSEIF (IRANGM.GT.INPPIM) THEN
          INPPIM=IRANGM
          INPIMF=INPPIM
        ENDIF
C
      ENDIF
C
  314 CONTINUE
      INALPI=MIN0 (INALPP,INBALO-INABAL)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(JRGPIF .EQ. 1)THEN
        CFICJD='FICJD'
        CFICJDOUT='FICJDOUT'
        CALL FMATTR(CFICJD,CFICJDOUT,IFICJD,IREP)
        OPEN(UNIT=IFICJD,FILE=CFICJD,FORM='FORMATTED')
      ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C        Balayage de la Paire d'Article d'Index concernee.
C
      DO 318 J=1,INALPI
C
      IF (CNOMAR(IXC(J,IRGPIM)).NE.' ') THEN
C
C              Il s'agit d'un article logique de donnees; en plus de ses
C         caracteristiques tabulees, on verifie s'il n'y a pas de la
C         place "perdue" juste derriere les donnees, place recuperable
C         eventuellement en cas de reecriture plus longue de l'article
C         logique.
C
        INALDI=INALDI+1
        ILONGA=MLGPOS(IXM(2*J-1,IRGPIM))
        IPOSDE=MLGPOS(IXM(2*J  ,IRGPIM))
        IPOSFI=IPOSDE+ILONGA-1
C
        IF (J.EQ.1.AND.JRGPIF.GT.INBPIR) THEN
C
C          Cas du premier article logique d'une P.A.I. excedentaire;
C     dans ce cas, la P.A.I. est situee derriere l'article logique,
C     en occupant deux articles physiques.
C
          IRECPI=MDES1D(IXM(ILARPH+1-(JRGPIF-INBPIR),IRANG))
          IDERPU=ILARPH*(IRECPI-1)
C
        ELSEIF (J.EQ.INALPI.AND.JRGPIF.EQ.INTPPI) THEN
C
C          Cas du dernier article logique du fichier, sans P.A.I. situee
C     derriere: la derniere position utilisable sans modifier le nombre
C     d'articles physiques du fichier correspond a la fin du dernier
C     article physique contenant des donnees, ou a la fin du dernier
C     article physique ecrit sur le fichier.
C
          IMDESC=MDES1D(IXM(JPNAPH,IRANG))
          IREC=MAX0 (1+(IPOSFI-1)/ILARPH,IMDESC)
          IDERPU=ILARPH*IREC
C
C          Si on arrive au test ci-dessous, on est sur que l'article lo-
C     gique n'est pas le dernier du fichier.
C
        ELSEIF (J.NE.INALPP) THEN
C
C          Cas general, ou l'article logique n'est pas le dernier de sa
C     (Paire de) Page(s) d'Index.
C
          IDERPU=MLGPOS(IXM(2*J+2,IRGPIM))-1
C
        ELSE
C
C          Cas particulier ou l'article logique est le dernier de sa
C     (Paire de) Page(s) d'Index.
C
          IDERPU=MLGPOS(IXM(2,IRPIMS))-1
        ENDIF
C
        IF (IDERPU.EQ.IPOSFI) THEN
C
          IF (LFRANC) THEN
            WRITE (UNIT=*,FMT='(I7,''-eme article de donnees: "'',A,
     S             ''",'',I7,'' mots, position'',I9,'' a'',I9)')
     S       INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,IPOSFI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            WRITE (UNIT=IFICJD,FMT='(I7,''  '',A,''  '',I8)')
     S       INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ELSE
            WRITE (UNIT=*,FMT='(I7,''-th data record: "'',A,''",'',I7,
     S             '' words, position'',I9,'' to'',I9)')
     S       INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,IPOSFI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            WRITE (UNIT=IFICJD,FMT='(I7,''  '',A,''  '',I8)')
     S       INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ENDIF
C
        ELSE
C
C           On visualise en plus la place "perdue" derriere l'article.
C
          IF (LFRANC) THEN
            WRITE (UNIT=*,FMT='(I7,''-eme article de donnees: "'',A,
     S             ''",'',I7,'' mots, position'',I9,'' a'',I9,'' <'',SP,
     S             I8,'' >'')')
     S   INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,IPOSFI,IDERPU-IPOSFI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            WRITE (UNIT=IFICJD,FMT='(I7,''  '',A,''  '',I8)')
     S       INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          ELSE
            WRITE (UNIT=*,FMT='(I7,''-th data record: "'',A,''",'',I7,
     S             '' words, position'',I9,'' to'',I9,'' <'',SP,
     S             I8,'' >'')')
     S   INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,IPOSFI,IDERPU-IPOSFI
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            WRITE (UNIT=IFICJD,FMT='(I7,''  '',A,''  '',I8)')
     S       INALDI,CNOMAR(IXC(J,IRGPIM)),ILONGA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          ENDIF
C
        ENDIF
C
      ELSEIF (LDTOUT) THEN
        INTROI=INTROI+1
        ILONGA=MLGPOS(IXM(2*J-1,IRGPIM))
        IPOSDE=MLGPOS(IXM(2*J  ,IRGPIM))
        IPOSFI=IPOSDE+ILONGA-1
C
        IF (LFRANC) THEN
          WRITE (UNIT=*,FMT='(TR1,5(''=''),''>'',T10,I6,
     S ''-eme TROU repertorie dans l''''index, longueur reutilisable:'',
     S         I7,'' mots, position'',I9,'' a'',I9)')
     S   INTROI,ILONGA,IPOSDE,IPOSFI
        ELSE
          WRITE (UNIT=*,FMT='(TR1,5(''=''),''>'',T10,I6,
     S ''-th HOLE cataloged within index, re-usable length:'',
     S         I7,'' words, position'',I9,'' to'',I9)')
     S   INTROI,ILONGA,IPOSDE,IPOSFI
        ENDIF
C
      ENDIF
C
  318 CONTINUE
C
      INABAL=INABAL+INALPI
  319 CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CLOSE(IFICJD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC JDJDJD CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCC didier
      call FMFREE(CFICJD,CFICJDOUT,IREP)
CCCCCCCCCCCCCCC didier
C*
C     3.2 -  ENVOI DE MESSAGES RECAPITULATIFS.
C-----------------------------------------------------------------------
C
      IF (LFRANC) THEN
C
        IF (LDTOUT) THEN
          WRITE (UNIT=*,FMT='(//,T5,8(''-''),TR3,I7,
     S           '' articles logiques de donnees et'',I6,
     S           '' trous repertories listes'',TR3,8(''-''),//)')
     S    INALDI,INTROI
        ELSE
          WRITE (UNIT=*,FMT='(//,T5,8(''-''),TR3,I7,
     S       '' articles logiques de donnees listes'',TR3,8(''-''),//)')
     S    INALDI
        ENDIF
C
      ELSE
C
        IF (LDTOUT) THEN
          WRITE (UNIT=*,FMT='(//,T5,8(''-''),TR3,I7,
     S           '' logical records of data and'',I6,
     S           '' holes within index listed'',TR3,8(''-''),//)')
     S    INALDI,INTROI
        ELSE
          WRITE (UNIT=*,FMT='(//,T5,8(''-''),TR3,I7,
     S       '' logical records of data listed'',TR3,8(''-''),//)')
     S    INALDI
        ENDIF
C
      ENDIF
C
      IF (INALDI.EQ.INALDO.AND.(.NOT.LDTOUT.OR.INTROI.EQ.INTROU)) THEN
C
        IF (LFRANC) THEN
          WRITE (UNIT=CLMESS,FMT=
     S     '(''Fin du catalogue de l''''Unite Logique'',I3,'' ---'',I7,
     S       '' Articles logiques en tout'')') KNUMER,INBALO
        ELSE
          WRITE (UNIT=CLMESS,FMT=
     S     '(''End of catalog of Logical Unit'',I3,'' ---'',I7,
     S       '' logical Records for whole file'')') KNUMER,INBALO
        ENDIF
C
        CALL LFIEMS (KNUMER,INIMES,IREP,LLFATA,CLMESS,CLNSPR,CLACTI)
        WRITE (UNIT=*,FMT='(///)')
      ELSE
        IREP=-16
      ENDIF
C
      GOTO 1001
C**
C     9.  - CI-DESSOUS, ETIQUETTES DE BRANCHEMENT EN CAS D'ERREUR E/S.
C-----------------------------------------------------------------------
C
  903 CONTINUE
      CLACTI='WRITE'
      GOTO 909
C
  904 CONTINUE
      CLACTI='READ'
C
  909 CONTINUE
C
C      AU CAS OU, ON FORCE LE CODE-REPONSE ENTREE/SORTIE A ETRE POSITIF.
C
      IREP=IABS (IREP)
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "LFIEMS" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
      IF (IRANG.NE.0) THEN
        NDEROP(IRANG)=18
        NDERCO(IRANG)=IREP
        IF (LMULTI) CALL LFIVER (VERRUE(IRANG),'OFF')
      ENDIF
      print *,' jdlfilaf BALISE 7'
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        RETURN
      ENDIF
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S    '', LDTOUT= '',L1)') KREP,KNUMER,LDTOUT
      CALL LFIEMS (KNUMER,INIMES,IREP,LLFATA,CLMESS,CLNSPR,CLACTI)
      print *,' jdlfilaf BALISE 8'
C
      RETURN
      END
