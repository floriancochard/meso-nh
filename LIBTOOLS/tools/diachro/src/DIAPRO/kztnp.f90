!     ######spl
      SUBROUTINE KZTNP(K)
!     ###################
!
!!****  *KZTNP* - 
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
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist 
!!                         (former NCAR common)
!!
!!       NIOFFD     : Label normalisation (=0 none, =/=0 active)
!!       NULBLL     : Nb of contours between 2 labelled contours
!!       NIOFFM     : =0    --> message at picture bottom
!!                    =/= 0 --> no message
!!       NIOFFP     : Special point value detection
!!                    (=0 none, =/=0 active)
!!       NHI        : Extrema detection
!!                    (=0 --> H+L, <0 nothing)
!!       NINITA     : For streamlimes
!!       NINITB     : Not yet implemented
!!       NIGRNC     : Not yet implemented
!!       NDOT       : Line style
!!                    (=0|1|1023|65535 --> solid lines;
!!                    <0 --> solid lines for positive values and
!!                    dotted lines(ABS(NDOT))for negative values;
!!                    >0 --> dotted lines(ABS(NDOT)) )
!!       NIFDC      : Coastline data style (0 none, 1 NCAR, 2 IGN)
!!       NLPCAR     : Number of land-mark points to be plotted
!!       NIMNMX     : Contour selection option
!!                    (=-1 Min, max and inc. automatically set;
!!                    =0 Min, max automatically set; inc. given;
!!                    >0 Min, max, inc. given by user)
!!       NISKIP     : Rate for drawing velocity vectors
!!       CTYPHOR    : Horizontal cross-section type
!!                    (='K' --> model level section;
!!                     ='Z' --> constant-altitude section;
!!                     ='P' --> isobar section (planned)
!!                     ='T' --> isentrope section (planned)
!!       XSPVAL     : Special value
!!       XSIZEL     : Label size
!!       XLATCAR, XLONCAR :  Lat. and Long. of land-mark points
!!       LXY        : If =.TRUE., plots  a grid-mesh stencil background
!!       LXZ        : If =.TRUE., plots  a model-level stencil background 
!!
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist 
!!                          (former PARA common)
!!
!!       XIDEBCOU, XJDEBCOU : Origin of a vertical cross-section
!!                            in cartesian (or conformal) real values
!!       XHMIN      : Altitude of the vert. cross-section
!!                    bottom (in meters above sea-level)
!!       XHMAX      : Altitude of the vert. cross-section
!!                    top (in meters above sea-level)
!!
!!
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
USE MODD_RESOLVCAR
USE MODD_MASK3D
USE MODD_ALLOC_FORDIACHRO
USE MODD_TYPE_AND_LH
USE MODN_NCAR    
USE MODN_PARA    

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------
INTEGER  :: K
!
!*       0.1   Local variables
!              ---------------

!
INTEGER   ::   J, JJ, JE
INTEGER   ::   IP1, IP2, IP3, IT
INTEGER   ::   JLOOPN, INDN, JF, JLOOPNF
INTEGER   ::   ILEN, INBGRA

REAL      ::   ZDIF
CHARACTER(LEN=8) :: YREP
!------------------------------------------------------------------------------
!
! Traitement des processus
!
IF(LPROCDIALL(K))THEN

  NBPROCDIA(K)=SIZE(XVAR,6)
  DO J=1,NBPROCDIA(K)
    NPROCDIA(J,K)=J
  ENDDO

ELSE

  IF(LPINCRDIA(K))THEN
    
    NPROCDIA(2,K)=MIN(NPROCDIA(2,K),SIZE(XVAR,6))

    IF(NBPROCDIA(K) == 2)THEN

      IP1=NPROCDIA(1,K)
      IP2=NPROCDIA(2,K)
      NBPROCDIA(K)=IP2-IP1+1
      JJ=0
      DO J=IP1,IP2
	JJ=JJ+1
	NPROCDIA(JJ,K)=J
      ENDDO

    ELSE IF(NBPROCDIA(K) == 3)THEN

      IP1=NPROCDIA(1,K)
      IP2=NPROCDIA(2,K)
      IP3=NPROCDIA(3,K)
      NBPROCDIA(K)=1
      DO J=2,100
	IP1=IP1+IP3
	IF(IP1 > IP2)EXIT
	NBPROCDIA(K)=NBPROCDIA(K)+1
	NPROCDIA(J,K)=IP1
      ENDDO

    ENDIF

  ENDIF

ENDIF

LPINCRDIA(K)=.FALSE.

IF(NBPROCDIA(K) == 0)THEN
  NPROCDIA(:,K)=0
ENDIF

!
! Traitement des numeros de masques et trajectoires 
!
IF(LNDIALL(K))THEN
  
  NBNDIA(K)=SIZE(XVAR,5)
  DO J=1,NBNDIA(K)
    NNDIA(J,K)=J
  ENDDO

ELSE

  IF(LNINCRDIA(K))THEN
   
    NNDIA(2,K)=MIN(NNDIA(2,K),SIZE(XVAR,5))

    IF(NBNDIA(K) == 2)THEN

      IP1=NNDIA(1,K)
      IP2=NNDIA(2,K)
      NBNDIA(K)=IP2-IP1+1
      JJ=0
      DO J=IP1,IP2
	JJ=JJ+1
	NNDIA(JJ,K)=J
      ENDDO

    ELSE IF(NBNDIA(K) == 3)THEN

      IP1=NNDIA(1,K)
      IP2=NNDIA(2,K)
      IP3=NNDIA(3,K)
      NBNDIA(K)=1
      DO J=2,100
	IP1=IP1+IP3
	IF(IP1 > IP2)EXIT
	NBNDIA(K)=NBNDIA(K)+1
	NNDIA(J,K)=IP1
      ENDDO

    ENDIF

  ENDIF

ENDIF

LNINCRDIA(K)=.FALSE.

IF(NBNDIA(K) == 0)THEN
  NNDIA(:,K)=0
ENDIF
!
! Traitement des temps
!
SELECT CASE(CTYPE)
  CASE('MASK','SSOL','SPXY')
    JLOOPNF=1
  CASE DEFAULT
    JLOOPNF=NBNDIA(K)
END SELECT

DO JLOOPN=1,JLOOPNF  ! Boucle sur les Num traj ou stations

SELECT CASE(CTYPE)
  CASE('MASK','SSOL','SPXY')
    INDN=1
  CASE DEFAULT
    INDN=NNDIA(JLOOPN,K)
END SELECT

SELECT CASE(CTYPE)
  CASE('CART','MASK','SPXY','SSOL')
    JF=SIZE(XVAR,4)
  CASE DEFAULT
    DO JE=SIZE(XTRAJT,1),1,-1
      IF(XTRAJT(JE,INDN) /= -1.E-15)THEN
	JF=JE
	EXIT
      ENDIF
    ENDDO
END SELECT

IF(LTIMEDIALL(K,INDN))THEN

  LTINCRDIA(K,INDN)=.TRUE.
  NBTIMEDIA(K,INDN)=3
  NTIMEDIA(1,K,INDN)=1
  NTIMEDIA(2,K,INDN)=JF
  NTIMEDIA(3,K,INDN)=1

  XTIMEDIA(1,K,INDN)=XTRAJT(NTIMEDIA(1,K,INDN),INDN)
  XTIMEDIA(2,K,INDN)=XTRAJT(NTIMEDIA(2,K,INDN),INDN)

ELSE

  IF(LTINCRDIA(K,INDN))THEN
! Incremental

    IF(NTIMEDIA(2,K,INDN) /=  0)THEN
      NTIMEDIA(2,K,INDN)=MIN(NTIMEDIA(2,K,INDN),JF)
    ENDIF

    IF(NBTIMEDIA(K,INDN) == 2)THEN

      IP1=NTIMEDIA(1,K,INDN)
      IP2=NTIMEDIA(2,K,INDN)
      IF(IP1 /=0 .AND. IP2 /=0)THEN
        NBTIMEDIA(K,INDN)=3
        NTIMEDIA(3,K,INDN)=1
        XTIMEDIA(1,K,INDN)=XTRAJT(NTIMEDIA(1,K,INDN),INDN)
        XTIMEDIA(2,K,INDN)=XTRAJT(NTIMEDIA(2,K,INDN),INDN)
! CONTROLER LA VALIDITE DES VALEURS

      ELSE

	DO J=1,JF
	  IF(XTIMEDIA(1,K,INDN) <= XTRAJT(J,INDN))EXIT
	ENDDO
	NTIMEDIA(1,K,INDN)=J
	DO J=1,JF
	  IF(XTIMEDIA(2,K,INDN) <= XTRAJT(J,INDN))EXIT
        ENDDO
	NTIMEDIA(2,K,INDN)=J
        NTIMEDIA(2,K,INDN)=MIN(NTIMEDIA(2,K,INDN),JF)
        NBTIMEDIA(K,INDN)=3
	NTIMEDIA(3,K,INDN)=1
      ENDIF

    ELSE IF(NBTIMEDIA(K,INDN) == 3)THEN

      IP1=NTIMEDIA(1,K,INDN)
      IP2=NTIMEDIA(2,K,INDN)
      IP3=NTIMEDIA(3,K,INDN)
      IF(IP1 /=0 .AND. IP2 /=0 .AND. IP3 /=0)THEN
        XTIMEDIA(1,K,INDN)=XTRAJT(NTIMEDIA(1,K,INDN),INDN)
        XTIMEDIA(2,K,INDN)=XTRAJT(NTIMEDIA(2,K,INDN),INDN)

      ELSE

	
	DO J=1,JF
	  IF(XTIMEDIA(1,K,INDN) <= XTRAJT(J,INDN))EXIT
	ENDDO
	NTIMEDIA(1,K,INDN)=J
	DO J=1,JF
	  IF(XTIMEDIA(2,K,INDN) <= XTRAJT(J,INDN))EXIT
        ENDDO
	NTIMEDIA(2,K,INDN)=J
        NTIMEDIA(2,K,INDN)=MIN(NTIMEDIA(2,K,INDN),JF)
	ZDIF=ABS(XTRAJT(2,INDN)-XTRAJT(3,INDN))
	IT=ANINT(XTIMEDIA(3,K,INDN)/ZDIF)
	NTIMEDIA(3,K,INDN)=IT
      ENDIF

    ENDIF

! Non incremental
  ELSE
    DO J=1,NBTIMEDIA(K,INDN)
      IF(NTIMEDIA(J,K,INDN) /= 0)THEN
	NTIMEDIA(J,K,INDN)=MIN(NTIMEDIA(J,K,INDN),JF)
	XTIMEDIA(J,K,INDN)=XTRAJT(NTIMEDIA(J,K,INDN),INDN)

      ELSE

	DO JJ=1,JF
	  IF(XTIMEDIA(J,K,INDN) <= XTRAJT(JJ,INDN))EXIT
        ENDDO
	NTIMEDIA(J,K,INDN)=JJ
	NTIMEDIA(J,K,INDN)=MIN(NTIMEDIA(J,K,INDN),JF)

      ENDIF
    ENDDO

  ENDIF

ENDIF
ENDDO      ! Fin boucle Num traj ou stations
!
! Traitement des niveaux de modele K
!
SELECT CASE(CTYPE)
  CASE('MASK')
! CASE('MASK','SSOL')
    JLOOPNF=1
  CASE DEFAULT
    JLOOPNF=NBNDIA(K)
END SELECT

DO JLOOPN=1,JLOOPNF  ! Boucle sur les Num traj ou stations

SELECT CASE(CTYPE)
  CASE('MASK')
! CASE('MASK','SSOL')
    INDN=1
  CASE DEFAULT
    INDN=NNDIA(JLOOPN,K)
END SELECT

SELECT CASE(CTYPE)
  CASE('CART','MASK','SPXY')
    JF=SIZE(XVAR,3)
  CASE('SSOL','DRST','RSPL','RAPL')
    DO JE=SIZE(XTRAJZ,1),1,-1
! Le 2eme indice (temps) est mis arbitrairement a 1 parce que la
! dimension en K pour le temps indice 1 est la meme que pour le
! temps indice n.
      IF(XTRAJZ(JE,1,INDN) /= -1.E-15)THEN
	JF=JE
	NKL=1
	NKH=JF
	EXIT
      ENDIF
    ENDDO
END SELECT

IF(LVLKDIALL(K,INDN))THEN

  NBLVLKDIA(K,INDN)=JF
  DO J=1,NBLVLKDIA(K,INDN)
    NLVLKDIA(J,K,INDN)=J+NKL-1
  ENDDO

ELSE

  IF(LKINCRDIA(K,INDN))THEN

    IF(NBLVLKDIA(K,INDN) == 2)THEN

      IP1=MAX(NLVLKDIA(1,K,INDN),NKL)
      IP2=MIN(NLVLKDIA(2,K,INDN),NKH)
      NBLVLKDIA(K,INDN)=IP2-IP1+1
      JJ=0
      DO J=IP1,IP2
	JJ=JJ+1
	NLVLKDIA(JJ,K,INDN)=J
      ENDDO

    ELSE IF(NBLVLKDIA(K,INDN) == 3)THEN

      IP1=MAX(NLVLKDIA(1,K,INDN),NKL)
      IP2=MIN(NLVLKDIA(2,K,INDN),NKH)
      IP3=NLVLKDIA(3,K,INDN)
      NLVLKDIA(1,K,INDN)=IP1
      NLVLKDIA(2,K,INDN)=IP2
      NBLVLKDIA(K,INDN)=1
      DO J=2,1000
	IP1=IP1+IP3
	IF(IP1 > IP2)EXIT
	NBLVLKDIA(K,INDN)=NBLVLKDIA(K,INDN)+1
	NLVLKDIA(J,K,INDN)=IP1
      ENDDO

    ENDIF

  ENDIF

ENDIF

LKINCRDIA(K,INDN)=.FALSE.

IF(NBLVLKDIA(K,INDN) == 0)THEN
  NLVLKDIA(:,K,INDN)=0
ENDIF
ENDDO      ! Fin boucle Num traj ou stations
!
! Traitement des altitudes Z
!
! On a directement  les altitudes en numerique en incremental ou non.
! Si (LZINCRDIA(K))  -->   NBLVLZDIA(K)=3
!                          XLVLZDIA(1:3,K)= extremes + increment
! Si (.NOT.LZINCRDIA(K))  -->    NBLVLZDIA(K)=N
!                                XLVLZDIA(1:N,K)=altitudes
!
!
! Positionnement de CTYPHOR
!
SELECT CASE(CTYPE)
  CASE('CART','MASK','SPXY')
    CTYPHOR(1:LEN(CTYPHOR))=' '
    IF(NBLVLKDIA(K,1) == 0 .AND. NBLVLZDIA(K) /=0 )THEN
      IF(LPR)THEN
        CTYPHOR='P'
      ELSE IF(LTK)THEN
        CTYPHOR='T'
      ELSE IF(LEV)THEN
        CTYPHOR='E'
      ELSE IF(LSV3)THEN
        CTYPHOR='V'
      ELSE
        CTYPHOR='Z'
      ENDIF
      LHORIZ=.TRUE.; LVERTI=.FALSE.
    ELSE IF(NBLVLKDIA(K,1) /= 0 .AND. NBLVLZDIA(K) ==0 )THEN
      CTYPHOR='K'
      LHORIZ=.TRUE.; LVERTI=.FALSE.

      IF(LTINCRDIA(K,1))THEN
	ILEN=(NTIMEDIA(2,K,1)-NTIMEDIA(1,K,1))/NTIMEDIA(3,K,1)+1
      ELSE
	ILEN=NBTIMEDIA(K,1)
      ENDIF

      INBGRA=NBPROCDIA(K)*NBLVLKDIA(K,1)*ILEN

      IF(INBGRA > 35 .AND. LCH .AND. CTYPE /= 'SPXY')THEN
	print *,'VOUS AVEZ DEMANDE: ',NBLVLKDIA(K,1),' NIVEAUX * ',  &
&       ILEN,' TEMPS * ',NBPROCDIA(K),' PROCESSUS = '
	print *,INBGRA,' GRAPHIQUES '
	print *,' EN ETES VOUS SUR ???? (y/n) '
	YREP(1:LEN(YREP))=' '
	READ(5,*)YREP
	SELECT CASE(YREP)
	  CASE('y','Y','o','O','yes','YES','oui','OUI')
	  CASE DEFAULT
	    LPBREAD=.TRUE.
	    print *,' VERIFIEZ LA SYNTAXE DE VOTRE DIRECTIVE ET RENTREZ LA A ',&
&           'NOUVEAU'
	END SELECT
      ENDIF

    ENDIF
  CASE DEFAULT
END SELECT



!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE KZTNP
