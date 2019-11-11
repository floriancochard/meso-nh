!     ######spl
      MODULE MODI_WRITE_OTHERSFIELDS
!     ##############################
!
INTERFACE
!
SUBROUTINE WRITE_OTHERSFIELDS(K,HFILEDIA,HLUOUTDIA,KX,KY,KZ)
INTEGER :: K
CHARACTER(LEN=*) :: HFILEDIA,HLUOUTDIA
INTEGER, INTENT(IN), OPTIONAL :: KX,KY,KZ
END SUBROUTINE WRITE_OTHERSFIELDS
!
END INTERFACE
!
END MODULE MODI_WRITE_OTHERSFIELDS
!     #############################################################
      SUBROUTINE WRITE_OTHERSFIELDS(K,HFILEDIA,HLUOUTDIA,KX,KY,KZ)
!     #############################################################
!
!!****  *WRITE_OTHERSFIELDS* - 
!! 
!!
!!    PURPOSE
!!    -------
! 
!
!!**  METHOD
!!    ------
!!      
!!
!!    REFERENCE
!!    ---------
!!     
!!
!!    AUTHORS
!!    -------
!!    J. Duron      *Lab. Aerologie* 
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    30/01/96 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE  MODD_DIMGRID_FORDIACHRO        
USE  MODD_OUT_DIA
USE  MODD_DIACHRO  
USE  MODD_ALLOC_FORDIACHRO  
USE  MODI_ALLOC_FORDIACHRO
USE  MODD_PARAMETERS
USE  MODD_DIM1
USE  MODD_TYPE_AND_LH
USE  MODD_RESOLVCAR, ONLY : CGROUP
USE  MODD_GRID
USE  MODD_CONF
USE  MODD_GRID1
USE  MODD_TIME1
USE MODD_TYPE_DATE
USE  MODI_WRITE_DIACHRO
USE  MODI_READ_DIACHRO
USE  MODI_RESOLV_UNITS
USE  MODI_TEMPORAL_DIST
USE  MODD_TIME
USE  MODI_FMREAD
USE  MODI_FMWRIT
!
IMPLICIT NONE
!
!*       0.1   Dummy arguments
!
INTEGER           :: K   ! Input file number
CHARACTER(LEN=*)  :: HFILEDIA, HLUOUTDIA
INTEGER, INTENT(IN), OPTIONAL :: KX,KY,KZ
!
!*       0.2   Local variables declarations
!
INTEGER           :: JJ, J, JA, I
INTEGER           :: ixyz, J1, J2, J3, I1, I2, I3
INTEGER           :: IIU, IJU, IKU
INTEGER           :: IGRID, ILENCH, IRESP
INTEGER           :: IPCENT
INTEGER           :: IMULT, ILYCOMM
INTEGER           :: ILUOUTDIA
!
CHARACTER(LEN=100):: YCOMMENT, YCAROUT
CHARACTER(LEN=20):: YCOMM
CHARACTER(LEN=16) :: YRECFM
!
REAL,DIMENSION(:),ALLOCATABLE  :: ZTAB
REAL,DIMENSION(:,:,:),ALLOCATABLE  :: ZTAB3, ZTABM3, Z3D
INTEGER,DIMENSION(3):: ITAB3  ! sizes of array ZTAB3
!
TYPE (DATE_TIME), SAVE :: TZDTEXP  ! to store exp. time when TT files
LOGICAL :: GPACK  ! to store LPACK
!----------------------------------------------------------------------------
!
!*       1.    INITIALISATION
!              --------------
!
GPACK=LPACK
! Duplication du profil au niveau des points de garde en 1D ou 2D
IF(NIMAX==1 .OR. NJMAX==1) LPACK=.FALSE.
!
ILENCH=LEN(YCOMMENT)
ILYCOMM=LEN(YCOMM)
YCOMM(1:ILYCOMM)='NOTHING'
!
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
print*,'IIU,IJU,IKU= ',IIU,IJU,IKU
!JDJDJDJD 291196
WRITE(NLUOUTD,*)' ******** WRITE_OTHERSFIELDS ENTREE CSTORAGE_TYPE ',CSTORAGE_TYPE
IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE == 'SU')THEN
  IKU=1
ENDIF
!JDJDJDJD 291196

CALL FMLOOK(HLUOUTDIA,HLUOUTDIA,ILUOUTDIA,IRESP)
!
! resolution degradee
ixyz=0
IF (PRESENT(KX)) THEN
IF (KX>1.AND.NIMAX/=1) ixyz=1 
ENDIF
IF (PRESENT(KY)) THEN 
IF (KY>1.AND.NJMAX/=1) ixyz=ixyz+10
ENDIF
IF (PRESENT(KZ)) THEN 
IF (KZ>1)              ixyz=ixyz+100
ENDIF
!
! NNB= NB d'articles a lire dans le fichier en cours de traitement en entree
! Mais en fait on prend comme ref. les articles du premier fichier
! (CF instruction IF(NNUMT(JJ,1....) en supposant que tous les fichiers
! traites ont la meme organisation (ce qui doit etre le cas sachant que
! ces fichiers sont differentes echeances d'un meme run)
!
DO JJ=1,NNB
!
  IF(NNUMT(JJ,1) /= 0)THEN
!
!----------------------------------------------------------------------------
!
!*       2.    TREATMENT ACCORDING THE VARIABLE SHAPE
!              --------------------------------------
!
! 130198 Introduction de IMULT pour prise en compte du 2D Vertical dont
! seul le plan central est enregistre
    IMULT=1
!
!*       2.0  
!
    IF(NSIZT(JJ,K) == IIU*IJU)THEN
! 051296 Modif pour tenir compte du 2D surfacique horizontal
      IKU=1
    ELSE IF(NSIZT(JJ,K)*(2*JPHEXT+NJMAX) == IIU*IJU)THEN
! 130198 Modif pour tenir compte du 2D Vertical filaire et surfacique; cas
! enregistrement du seul plan central
      IKU=1
      IMULT=2*JPHEXT+NJMAX
      WRITE(NLUOUTD,*)'***************************************************************'
      WRITE(NLUOUTD,*)' Variable 1D rencontree // X et enregistree dans le fichier',&
      &' diachronique ',CRECFM2T(JJ,K),' size origine et size enr.  ',NSIZT(JJ,K),NSIZT(JJ,K)*IMULT
      WRITE(NLUOUTD,*)' (Duplication du profil (<--> 2D filaire) au niveau des points de garde)'
      WRITE(NLUOUTD,*)'***************************************************************'
    ELSE IF(NSIZT(JJ,K)*(2*JPHEXT+NIMAX) == IIU*IJU)THEN
      WRITE(NLUOUTD,*)'***************************************************************'
      WRITE(NLUOUTD,*)' Variable 1D // Y non enregistree dans le fichier',&
     &' diachronique ',CRECFM2T(JJ,K),' size et IIU,IJU,IKU ',NSIZT(JJ,K),IIU,IJU,IKU
      WRITE(NLUOUTD,*)'***************************************************************'
      CYCLE
    ELSE IF(NSIZT(JJ,K)*(2*JPHEXT+NIMAX) == IIU*IJU/(2*JPHEXT+NJMAX))THEN
      IF(NIMAX==1 .AND. NJMAX==1) THEN
! 110906 Cas 0D Vertical ou seul le profil central est enregistre
!        Duplication du profil sur les points de garde
!        (rigoureusement, il faut dupliquer car type CART)
        IKU=1
        IMULT = (2*JPHEXT+NIMAX)*(2*JPHEXT+NJMAX)
        WRITE(NLUOUTD,*)'***************************************************************'
        WRITE(NLUOUTD,*)' Variable 0D enregistree dans le fichier',&
      &' diachronique ',CRECFM2T(JJ,K),' size origine et size enr.  ',NSIZT(JJ,K),NSIZT(JJ,K)*IMULT
        WRITE(NLUOUTD,*)' (Duplication du profil au niveau des points de garde...)'
        WRITE(NLUOUTD,*)'***************************************************************'
      ENDIF
    ELSE
      IKU=NKMAX+2*JPVEXT 
      IF(NSIZT(JJ,K)*(2*JPHEXT+NJMAX) == IIU*IJU*IKU)THEN
        IMULT=2*JPHEXT+NJMAX
        WRITE(NLUOUTD,*)'***************************************************************'
        WRITE(NLUOUTD,*)' Variable 2D Vertical // X et enregistree dans le fichier',&
      &' diachronique ',CRECFM2T(JJ,K),' size origine et size enr.  ',NSIZT(JJ,K),NSIZT(JJ,K)*IMULT
        WRITE(NLUOUTD,*)' (Duplication du plan au niveau des points de garde)'
        WRITE(NLUOUTD,*)'***************************************************************'
      ELSE IF(NSIZT(JJ,K)*(2*JPHEXT+NIMAX) == IIU*IJU*IKU)THEN
        WRITE(NLUOUTD,*)'***************************************************************'
        WRITE(NLUOUTD,*)' Variable 2D Vertical // Y non enregistree dans le fichier',&
     &' diachronique ',CRECFM2T(JJ,K),' size et IIU,IJU,IKU ',NSIZT(JJ,K),IIU,IJU,IKU
        WRITE(NLUOUTD,*)'***************************************************************'
        CYCLE
      !ELSE IF(NSIZT(JJ,K)*(2*JPHEXT+NIMAX)*(2*JPHEXT+NJMAX) == IIU*IJU*IKU)THEN
      !remplace par la ligne suivante car le membre de gauche peut etre tres grand
      ELSE IF(NSIZT(JJ,K)*(2*JPHEXT+NIMAX)==IIU*IJU*IKU/(2*JPHEXT+NJMAX) )THEN
        WRITE(NLUOUTD,*)'***************************************************************'
        IF(NIMAX==1 .AND. NJMAX==1) THEN
! 180703 Cas 1D Vertical ou seul le profil central est enregistre
!        Duplication du profil sur les points de garde
!        (rigoureusement, il faut dupliquer car type CART)
          IMULT = (2*JPHEXT+NIMAX)*(2*JPHEXT+NJMAX)
        WRITE(NLUOUTD,*)' Variable 1D Vertical enregistree dans le fichier',&
      &' diachronique ',CRECFM2T(JJ,K),' size origine et size enr.  ',NSIZT(JJ,K),NSIZT(JJ,K)*IMULT
          WRITE(NLUOUTD,*)' (Duplication du profil au niveau des points de garde...)'
        ELSE
          WRITE(NLUOUTD,*)' Variable 1D Vertical enregistree dans le fichier',&
      &' diachronique ',CRECFM2T(JJ,K),' size origine et size enr.  ',NSIZT(JJ,K),NSIZT(JJ,K)*IMULT
        ENDIF
        WRITE(NLUOUTD,*)'***************************************************************'
      ELSE
        IF(NSIZT(JJ,K) == IIU*IJU*IKU)THEN
! Variable 3D normale IKU= NKMAX+2*JPVEXT IMULT=1 On ne fait rien
        ELSE
          IF(NJMAX==1 .AND. GPACK) THEN
            IF(MOD(NSIZT(JJ,K) , IIU) == 0)THEN
! Variable 3D avec la 3eme dim  <= a IKU habituel et sans signification spatiale
              IKU=NSIZT(JJ,K)/IIU
              WRITE(NLUOUTD,*)'*********** 3D mais 3e dimension =/= de IKU *******************'
              WRITE(NLUOUTD,*)' Variable 3D enregistree dans le fichier diachronique ',&
       &CRECFM2T(JJ,K),' size et IIU,3e DIMENSION,IKU ',NSIZT(JJ,K),IIU,IKU,NKMAX+2*JPVEXT
              IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE == 'SU')THEN
                WRITE(NLUOUTD,*)' cas d un fichier physiographique: niveaux supplementaires de 1 a ',IKU
              ELSE
                WRITE(NLUOUTD,*)' consideree comme une matrice partielle en K dont seuls les niveaux 1 a ',IKU,' sont enregistres'
              END IF
!        Duplication du profil sur les points de garde
!        (rigoureusement, il faut dupliquer car type CART)
              IMULT = 2*JPHEXT+NJMAX
              WRITE(NLUOUTD,*)' (Duplication au niveau des points de garde)'
            ENDIF
          ELSE IF(MOD(NSIZT(JJ,K) , IIU*IJU) == 0)THEN
! Variable 3D avec la 3eme dim  <= a IKU habituel et sans signification spatiale
            IKU=NSIZT(JJ,K)/(IIU*IJU)
            WRITE(NLUOUTD,*)'*********** 3D mais 3e dimension =/= de IKU *******************'
            WRITE(NLUOUTD,*)' Variable 3D enregistree dans le fichier diachronique ',&
     &CRECFM2T(JJ,K),' size et IIU,IJU,3e DIMENSION,IKU ',NSIZT(JJ,K),IIU,IJU,IKU,NKMAX+2*JPVEXT
            IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE == 'SU')THEN
              WRITE(NLUOUTD,*)' cas d un fichier physiographique: niveaux supplementaires de 1 a ',IKU
            ELSE
              WRITE(NLUOUTD,*)' consideree comme une matrice partielle en K dont seuls les niveaux 1 a ',IKU,' sont enregistres'
            END IF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    !
!
! Allocation de la zone tampon de lecture
    ALLOCATE(ZTAB(NSIZT(JJ,K)))
    ! LPACK n intervient pas dans cette maniere de lire (ZTAB est 1D)
!
! Lecture de l'article concerne (CRECFM2T(JJ,K))
    CALL FMREAD(CNAMFILED(K),CRECFM2T(JJ,K),CLUOUTD,NSIZT(JJ,K), &
      ZTAB,IGRID,ILENCH,YCOMMENT,IRESP)
    YCOMMENT=ADJUSTL(ADJUSTR(YCOMMENT))
    CGROUP(1:LEN(CGROUP))=' '
    CGROUP=CRECFM2T(JJ,K)
    CGROUP=ADJUSTL(CGROUP)
!
! 051296 Modifs pour enregistrer le relief ZS egalement sous le nom ZSBIS
    IF(CGROUP(1:LEN_TRIM(CGROUP)) == 'ZS')THEN
      CRECFM2T(JJ,K)='ZSBIS'
      CGROUP='ZSBIS'
    ENDIF
! 120106 idem pour le smooth relief
    IF(CGROUP(1:LEN_TRIM(CGROUP)) == 'ZSMT')THEN
      CRECFM2T(JJ,K)='ZSMTBIS'
      CGROUP='ZSMTBIS'
    ENDIF
!
! Extraction des unites du champ commentaire
  YCAROUT(1:LEN(YCAROUT))=' '
  IF (LEN_TRIM(YCOMMENT)/=0) &
    CALL RESOLV_UNITS(YCOMMENT(1:LEN_TRIM(YCOMMENT)),YCAROUT)
!
!
!*       2.1  ++++3D + 2D H + 2D V et 1D // X+++++
!
! Traitement informations 3D et 2D Horiz. Sont considerees de type CART
! dans le fichier diachronique
! (En realite si on W en 2D, on recupere le 2D plan et filaire 
! (3D + 2D avec les points de garde) et si on W
! en 1D on recupere 1 profil vertical (3D avec les points de garde) et
! peut-etre 1 scalaire avec des points de garde horiz. (2D)) A VERIFIER

! 130198 Ajout 2D Vertical surfacique + filaire // X
!   IF(NSIZT(JJ,K) == IIU*IJU*IKU)THEN
    IF(NSIZT(JJ,K)*IMULT == IIU*IJU*IKU)THEN
      IF(IMULT /= 1)THEN
        IF(IMULT == (2*JPHEXT+NIMAX)*(2*JPHEXT+NJMAX))THEN
! 180703 Cas 1D Vertical ou seul le profil central est enregistre
!         si pas de duplication du profil sur les points de garde:
!         ITAB3(1)=1; ITAB3(2)=1; ITAB3(3)=IKU
!         ALLOCATE(ZTAB3(ITAB3(1),ITAB3(2),ITAB3(3)))
!         ZTAB3=RESHAPE(ZTAB,ITAB3)
!        il faut dupliquer car type CART:
          ITAB3(1)=1 ; ITAB3(2)=1 ; ITAB3(3)=IKU
          ALLOCATE(ZTABM3(ITAB3(1),ITAB3(2),ITAB3(3)))
          ZTABM3=RESHAPE(ZTAB,ITAB3)
          ITAB3(1)=2*JPHEXT+NIMAX ; ITAB3(2)=2*JPHEXT+NJMAX ; ITAB3(3)=IKU
          IF (ALLOCATED(ZTAB3)) DEALLOCATE(ZTAB3)
          ALLOCATE(ZTAB3(ITAB3(1),ITAB3(2),ITAB3(3)))
          DO J=1,ITAB3(2)
          DO I=1,ITAB3(1)
            ZTAB3(I,J,:)=ZTABM3(1,1,:)
          ENDDO
          ENDDO
          DEALLOCATE(ZTABM3)
        ELSE
! 130198 Cas 2D Vertical // X ou seul le plan central est enregistre
!        Duplication du plan sur les points de garde
          ITAB3(1)=IIU; ITAB3(2)=1; ITAB3(3)=IKU
          ALLOCATE(ZTABM3(ITAB3(1),ITAB3(2),ITAB3(3)))
          ZTABM3=RESHAPE(ZTAB,ITAB3)
          IF (ALLOCATED(ZTAB3)) DEALLOCATE(ZTAB3)
          ITAB3(1)=IIU; ITAB3(2)=IJU; ITAB3(3)=IKU
          ALLOCATE(ZTAB3(ITAB3(1),ITAB3(2),ITAB3(3)))
          DO J=1,ITAB3(2)
            ZTAB3(:,J,:)=ZTABM3(:,1,:)
          ENDDO
          DEALLOCATE(ZTABM3)
        END IF
      ELSE ! Variable 3D normale IKU= NKMAX+2*JPVEXT IMULT=1 
        ITAB3(1)=IIU; ITAB3(2)=IJU; ITAB3(3)=IKU
        IF (ALLOCATED(ZTAB3)) DEALLOCATE(ZTAB3)
        ALLOCATE(ZTAB3(ITAB3(1),ITAB3(2),ITAB3(3)))
        ZTAB3=RESHAPE(ZTAB,ITAB3)
      ENDIF
!
! Dans ce pg de conversion, on considere que chaque variable (prognostique,
! diagnostique, generique represente a elle seule un groupe a 1 processus
! (--> indice de processus = 1)
! On affecte (arbitrairement) le meme nom au groupe et au processus
      IF(K == 1)THEN      
        CTYPE='CART'
! resolution degradee
        IF (PRESENT(KX)) THEN
        IF (KX>1.AND.NIMAX/=1) ITAB3(1)=(IIU-1)/KX +1 
        ENDIF
        IF (PRESENT(KY)) THEN
        IF (KY>1.AND.NJMAX/=1) ITAB3(2)=(IJU-1)/KY +1 
        ENDIF
        IF (PRESENT(KZ)) THEN
        IF (KZ>1)              ITAB3(3)=(IKU-1)/KZ +1
        ENDIF
! Allocation des matrices utilisees dans le fichier diachronique (dernier
! argument = 1 pour ecriture; = 2 pour lecture; si =3, desallocation)
        CALL ALLOC_FORDIACHRO(ITAB3(1),ITAB3(2),ITAB3(3),NNBF,1,1,1)
! Initialisation de variables et matrices
        LICP=.FALSE. ; LJCP=.FALSE. ; LKCP=.FALSE.
        NIL=1 ; NJL=1 ; NKL=1
        NIH=ITAB3(1) ; NJH=ITAB3(2) ; NKH=ITAB3(3)
        XVAR(:,:,:,:,:,:)=0.
        XTRAJT(:,:)=0.
        CTITRE(:)(1:LEN(CTITRE))=' '
        CUNITE(:)(1:LEN(CUNITE))=' '
        CCOMMENT(:)(1:LEN(CCOMMENT))=' '
        XDATIME(:,:)=0.
      ENDIF
!
! Distinction 1er fichier et les suivants. Dans le premier cas on ecrit di-
! -rectement dans le fic. diachronique et apres les avoir reorganisees les
! informations lues. Dans les cas suivants, on relit d'abord les infos du
! fic. diachron. pour les augmenter des nouvelles  fraichement lues avant
! de les reecrire.
! NOTA on a pris la precaution de prevoir des le depart une taille d'article
! = a la dimension de la matrice traitee * par le nb de fichiers lus (NNBF)
!
      IF (K == 1)THEN                   !************************************
! resolution degradee
        ! in:  ZTAB3, taille:IIU(ou 1),IJU(ou 1),IKU
        ! out: XVAR,  taille:ITAB3
        SELECT CASE(ixyz)
        CASE (0)
          XVAR(:,:,:,K,1,1)=ZTAB3
        CASE (1)   !X
          DO J3=1,SIZE(ZTAB3,3)
          DO J2=1,SIZE(ZTAB3,2)
            XVAR(:,J2,J3,K,1,1)=ZTAB3(1:IIU:KX,J2,J3)
          END DO
          END DO
        CASE (10)  !Y
          DO J3=1,SIZE(ZTAB3,3)
          DO J1=1,SIZE(ZTAB3,1)
            XVAR(J1,:,J3,K,1,1)=ZTAB3(J1,1:IJU:KY,J3)
          END DO
          END DO
        CASE (11)  !X et Y
          DO J3=1,SIZE(ZTAB3,3)
            I2=0
            DO J2=1,SIZE(ZTAB3,2),KY
              I2=I2+1
              XVAR(:,I2,J3,K,1,1)=ZTAB3(1:IIU:KX,J2,J3)
            END DO
            IF (I2>SIZE(XVAR,2)) THEN
              print*,'cas xy: niveau ',J3,' debordement de tableau: ', &
                     I2,SIZE(XVAR,2)
              STOP
            ENDIF
          END DO
        CASE (100) !Z
          DO J2=1,SIZE(ZTAB3,2)
          DO J1=1,SIZE(ZTAB3,1)
            XVAR(J1,J2,:,K,1,1)=ZTAB3(J1,J2,1:IKU:KZ)
          END DO
          END DO
        CASE (101) !X et Z
          DO J2=1,SIZE(ZTAB3,2)
            I1=0
            DO J1=1,SIZE(ZTAB3,1),KX
              I1=I1+1
              XVAR(I1,J2,:,K,1,1)=ZTAB3(J1,J2,1:IKU:KZ)
            END DO
            IF (I1>SIZE(XVAR,1)) THEN
              print*,'cas xz: colonne ',J2,' debordement de tableau: ', &
                     I1,SIZE(XVAR,1)
              STOP
            ENDIF
          END DO
        CASE (110)  !Y et Z
          DO J1=1,SIZE(ZTAB3,1)
            I2=0
            DO J2=1,SIZE(ZTAB3,2),KY
              I2=I2+1
              XVAR(J1,I2,:,K,1,1)=ZTAB3(J1,I2,1:IKU:KZ)
              IF (I2>SIZE(XVAR,2)) THEN
                print*,'cas xy: ligne ',J1,' debordement de tableau: ', &
                       I2,SIZE(XVAR,2)
                STOP
              ENDIF
            END DO
          END DO
        CASE (111)  !X, Y et Z
          ALLOCATE(Z3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(ZTAB3,3)))
          !first X et Y
          DO J3=1,SIZE(ZTAB3,3)
            I2=0
            DO J2=1,SIZE(ZTAB3,2),KY
              I2=I2+1
              Z3D(:,I2,J3)=ZTAB3(1:IIU:KX,J2,J3)
            END DO
            IF (I2>SIZE(XVAR,2)) THEN
              print*,'cas xyz: niveau ',J3,' debordement de tableau: ', &
                     I2,SIZE(XVAR,2)
              STOP
            ENDIF
          END DO
          !then Z
          DO J2=1,SIZE(XVAR,2)
          DO J1=1,SIZE(XVAR,1)
            XVAR(J1,J2,:,K,1,1)=Z3D(J1,J2,1:IKU:KZ)
          END DO
          END DO
          DEALLOCATE(Z3D)
        END SELECT
!
! Le tps courant est transforme en temps relatif par / au debut de l'experience
        CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH, &
          TDTCUR%TDATE%DAY,TDTCUR%TIME,TDTEXP%TDATE%YEAR,    &
          TDTEXP%TDATE%MONTH,TDTEXP%TDATE%DAY,TDTEXP%TIME,           &
          XTRAJT(K,1))
        TZDTEXP=TDTEXP
        CTITRE(1)=CGROUP
        CUNITE(1)=ADJUSTL(YCAROUT)
        CCOMMENT(1)=YCOMMENT
        XDATIME(1,K)=TDTEXP%TDATE%YEAR; XDATIME(2,K)=TDTEXP%TDATE%MONTH
        XDATIME(3,K)=TDTEXP%TDATE%DAY;  XDATIME(4,K)=TDTEXP%TIME
        XDATIME(5,K)=TDTSEG%TDATE%YEAR; XDATIME(6,K)=TDTSEG%TDATE%MONTH
        XDATIME(7,K)=TDTSEG%TDATE%DAY;  XDATIME(8,K)=TDTSEG%TIME
        XDATIME(9,K)=TDTMOD%TDATE%YEAR; XDATIME(10,K)=TDTMOD%TDATE%MONTH
        XDATIME(11,K)=TDTMOD%TDATE%DAY; XDATIME(12,K)=TDTMOD%TIME
        XDATIME(13,K)=TDTCUR%TDATE%YEAR;XDATIME(14,K)=TDTCUR%TDATE%MONTH
        XDATIME(15,K)=TDTCUR%TDATE%DAY; XDATIME(16,K)=TDTCUR%TIME
!
! Ecriture dans le fichier diachronique
        NGRIDIA(1)=IGRID
        CALL WRITE_DIACHRO(HFILEDIA,HLUOUTDIA,CGROUP,CTYPE,NGRIDIA,XDATIME,XVAR, &
          XTRAJT,CTITRE,CUNITE,CCOMMENT, &
        LICP,LJCP,LKCP,NIL,NIH,NJL,NJH,NKL,NKH)
!
! Desallocation des matrices
        DEALLOCATE(ZTAB3)
        CALL ALLOC_FORDIACHRO(IIU,IJU,IKU,NNBF,1,1,3)
!
      ELSE                              !************************************
!
! On relit les infos deja enregistrees du fichier diachronique en connaissant
! le nom du groupe CGROUP=CRECFM2T(JJ,K)
        CALL READ_DIACHRO(CFILEDIA,CLUOUTDIA,CGROUP)
        SELECT CASE(ixyz)
        CASE (0)
          XVAR(:,:,:,K,1,1)=ZTAB3
        CASE (1)   !X
          DO J3=1,SIZE(ZTAB3,3)
          DO J2=1,SIZE(ZTAB3,2)
            XVAR(:,J2,J3,K,1,1)=ZTAB3(1:IIU:KX,J2,J3)
          END DO
          END DO
        CASE (10)  !Y
          DO J3=1,SIZE(ZTAB3,3)
          DO J1=1,SIZE(ZTAB3,1)
            XVAR(J1,:,J3,K,1,1)=ZTAB3(J1,1:IJU:KY,J3)
          END DO
          END DO
        CASE (11)  !X et Y
          DO J3=1,SIZE(ZTAB3,3)
            I2=0
            DO J2=1,SIZE(ZTAB3,2),KY
              I2=I2+1
              XVAR(:,I2,J3,K,1,1)=ZTAB3(1:IIU:KX,J2,J3)
            END DO
          END DO
        CASE (100) !Z
          DO J2=1,SIZE(ZTAB3,2)
          DO J1=1,SIZE(ZTAB3,1)
            XVAR(J1,J2,:,K,1,1)=ZTAB3(J1,J2,1:IKU:KZ)
          END DO
          END DO
        CASE (101) !X et Z
          DO J2=1,SIZE(ZTAB3,2)
            I1=0
            DO J1=1,SIZE(ZTAB3,1),KX
              I1=I1+1
              XVAR(I1,J2,:,K,1,1)=ZTAB3(J1,J2,1:IKU:KZ)
            END DO
          END DO
        CASE (110)  !Y et Z
          DO J1=1,SIZE(ZTAB3,1)
            I2=0
            DO J2=1,SIZE(ZTAB3,2),KY
              I2=I2+1
              XVAR(J1,I2,:,K,1,1)=ZTAB3(J1,I2,1:IKU:KZ)
            END DO
          END DO
        CASE (111)  !X, Y et Z
          ALLOCATE(Z3D(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(ZTAB3,3)))
          !first X et Y
          DO J3=1,SIZE(ZTAB3,3)
            I2=0
            DO J2=1,SIZE(ZTAB3,2),KY
              I2=I2+1
              Z3D(:,I2,J3)=ZTAB3(1:IIU:KX,J2,J3)
            END DO
          END DO
          !then Z
          DO J2=1,SIZE(XVAR,2)
          DO J1=1,SIZE(XVAR,1)
            XVAR(J1,J2,:,K,1,1)=Z3D(J1,J2,1:IKU:KZ)
          END DO
          END DO
          DEALLOCATE(Z3D)
        END SELECT
        CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH, &
          TDTCUR%TDATE%DAY,TDTCUR%TIME,TDTEXP%TDATE%YEAR,    &
          TDTEXP%TDATE%MONTH,TDTEXP%TDATE%DAY,TDTEXP%TIME,           &
        XTRAJT(K,1))
        IF (CSTORAGE_TYPE=='TT') THEN
          CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH, &
            TDTCUR%TDATE%DAY,TDTCUR%TIME,TZDTEXP%TDATE%YEAR,    &
            TZDTEXP%TDATE%MONTH,TZDTEXP%TDATE%DAY,TZDTEXP%TIME,        &
          XTRAJT(K,1))
          WRITE(NLUOUTD,*) &
          ' WRITE_OTHERSFIELDS calcul de XTRAJT par rapport au 1er fichier ',XTRAJT(K,1)
        END IF
        XDATIME(1,K)=TDTEXP%TDATE%YEAR; XDATIME(2,K)=TDTEXP%TDATE%MONTH
        XDATIME(3,K)=TDTEXP%TDATE%DAY;  XDATIME(4,K)=TDTEXP%TIME
        XDATIME(5,K)=TDTSEG%TDATE%YEAR; XDATIME(6,K)=TDTSEG%TDATE%MONTH
        XDATIME(7,K)=TDTSEG%TDATE%DAY;  XDATIME(8,K)=TDTSEG%TIME
        XDATIME(9,K)=TDTMOD%TDATE%YEAR; XDATIME(10,K)=TDTMOD%TDATE%MONTH
        XDATIME(11,K)=TDTMOD%TDATE%DAY; XDATIME(12,K)=TDTMOD%TIME
        XDATIME(13,K)=TDTCUR%TDATE%YEAR;XDATIME(14,K)=TDTCUR%TDATE%MONTH
        XDATIME(15,K)=TDTCUR%TDATE%DAY; XDATIME(16,K)=TDTCUR%TIME
!
        WRITE(ILUOUTDIA,*)' OTHERSFIELDS IGRID XVAR,XTRAJT,CTITRE,CUNITE,CCOMMENT'
        WRITE(ILUOUTDIA,*)IGRID,SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
        SIZE(XVAR,5),SIZE(XVAR,6),'  ',SIZE(XTRAJT,1),SIZE(XTRAJT,2),XTRAJT
        WRITE(ILUOUTDIA,*)(CTITRE(J)(1:LEN(CTITRE)),J=1,SIZE(CTITRE))
        WRITE(ILUOUTDIA,*)(CUNITE(J)(1:LEN(CUNITE)),J=1,SIZE(CUNITE))
        WRITE(ILUOUTDIA,*)(CCOMMENT(J)(1:LEN(CCOMMENT)),J=1,SIZE(CCOMMENT))

! Ecriture dans le fichier diachronique
        NGRIDIA(1)=IGRID
        CALL WRITE_DIACHRO(HFILEDIA,HLUOUTDIA,CGROUP,CTYPE,NGRIDIA,XDATIME,XVAR, &
        XTRAJT,CTITRE,CUNITE,CCOMMENT, &
        LICP,LJCP,LKCP,NIL,NIH,NJL,NJH,NKL,NKH)

! Desallocation des matrices
        DEALLOCATE(ZTAB3)
        CALL ALLOC_FORDIACHRO(IIU,IJU,IKU,NNBF,1,1,3)

      ENDIF                             !************************************
!
!
!*       2.2  ++++2D+++++
!
! Traitement des infos 2D (Traite avec le 3D)
!   ELSE IF(NSIZT(JJ,K) == IIU*IJU)THEN
!
!
!*       2.3  ++++1D // Z+++++
!
! Traitement des infos 1D
    ELSE IF(NSIZT(JJ,K) == IKU)THEN
      WRITE(NLUOUTD,*)'***************************************************************'
      WRITE(NLUOUTD,*)' Variable 1D rencontree et enregistree dans le fichier',&
      &' diachronique ',CGROUP,' size et IKU ',NSIZT(JJ,K),IKU
      WRITE(NLUOUTD,*)'***************************************************************'
      ITAB3(1)=1; ITAB3(2)=1; ITAB3(3)=IKU
      IF(ALLOCATED(ZTAB3))THEN
        DEALLOCATE(ZTAB3)
      ENDIF
      ALLOCATE(ZTAB3(ITAB3(1),ITAB3(2),ITAB3(3)))
      ZTAB3=RESHAPE(ZTAB,ITAB3)
!
      IF(K == 1)THEN      
        CTYPE='CART'
! resolution degradee
        IF (PRESENT(KZ)) THEN
        IF (KZ>1)              ITAB3(3)=(IKU-1)/KZ +1
        ENDIF
! Allocation des matrices utilisees dans le fichier diachronique (dernier
! argument = 1 pour ecriture; = 2 pour lecture; si =3, desallocation)
        CALL ALLOC_FORDIACHRO(ITAB3(1),ITAB3(2),ITAB3(3),NNBF,1,1,1)
! Initialisation de variables et matrices
        LICP=.FALSE. ; LJCP=.FALSE. ; LKCP=.FALSE.
        NIL=JPHEXT ; NJL=JPHEXT ; NKL=1
        NIH=JPHEXT ; NJH=JPHEXT ; NKH=ITAB3(3)
        XVAR(:,:,:,:,:,:)=0.
        XTRAJT(:,:)=0.
        CTITRE(:)(1:LEN(CTITRE))=' '
        CUNITE(:)(1:LEN(CUNITE))=' '
        CCOMMENT(:)(1:LEN(CCOMMENT))=' '
        XDATIME(:,:)=0
      ENDIF
!
! Distinction 1er fichier et les suivants. Dans le premier cas on ecrit di-
! -rectement dans le fic. diachronique et apres les avoir reorganisees les
! informations lues. Dans les cas suivants, on relit d'abord les infos du
! fic. diachron. pour les augmenter des nouvelles  fraichement lues avant
! de les reecrire.
! NOTA on a pris la precaution de prevoir des le depart une taille d'article
! = a la dimension de la matrice traitee * par le nb de fichiers lus (NNBF)
!
      IF (K == 1)THEN                   !************************************
        IF (PRESENT(KZ)) THEN
        DO J2=1,SIZE(ZTAB3,2)
        DO J1=1,SIZE(ZTAB3,1)
          XVAR(J1,J2,:,K,1,1)=ZTAB3(J1,J2,1:IKU:KZ)
        END DO
        END DO
        ELSE
        XVAR(:,:,:,K,1,1)=ZTAB3
        ENDIF
!
! Le tps courant est transforme en temps relatif par / au debut de l'experience
        CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH, &
          TDTCUR%TDATE%DAY,TDTCUR%TIME,TDTEXP%TDATE%YEAR,    &
          TDTEXP%TDATE%MONTH,TDTEXP%TDATE%DAY,TDTEXP%TIME,           &
          XTRAJT(K,1))
        CTITRE(1)=CGROUP
        CUNITE(1)=ADJUSTL(YCAROUT)
        CCOMMENT(1)=YCOMMENT
        XDATIME(1,K)=TDTEXP%TDATE%YEAR; XDATIME(2,K)=TDTEXP%TDATE%MONTH
        XDATIME(3,K)=TDTEXP%TDATE%DAY;  XDATIME(4,K)=TDTEXP%TIME
        XDATIME(5,K)=TDTSEG%TDATE%YEAR; XDATIME(6,K)=TDTSEG%TDATE%MONTH
        XDATIME(7,K)=TDTSEG%TDATE%DAY;  XDATIME(8,K)=TDTSEG%TIME
        XDATIME(9,K)=TDTMOD%TDATE%YEAR; XDATIME(10,K)=TDTMOD%TDATE%MONTH
        XDATIME(11,K)=TDTMOD%TDATE%DAY; XDATIME(12,K)=TDTMOD%TIME
        XDATIME(13,K)=TDTCUR%TDATE%YEAR;XDATIME(14,K)=TDTCUR%TDATE%MONTH
        XDATIME(15,K)=TDTCUR%TDATE%DAY; XDATIME(16,K)=TDTCUR%TIME

! Ecriture dans le fichier diachronique
        NGRIDIA(1)=IGRID
        CALL WRITE_DIACHRO(HFILEDIA,HLUOUTDIA,CGROUP,CTYPE,NGRIDIA,XDATIME,XVAR, &
        XTRAJT,CTITRE,CUNITE,CCOMMENT, &
        LICP,LJCP,LKCP,NIL,NIH,NJL,NJH,NKL,NKH)

! Desallocation des matrices
        DEALLOCATE(ZTAB3)
        CALL ALLOC_FORDIACHRO(1,1,IKU,NNBF,1,1,3)
!     
      ELSE                              !************************************
!
! On relit les infos deja enregistrees du fichier diachronique en connaissant
! le nom du groupe CGROUP=CRECFM2T(JJ,K)
        CALL READ_DIACHRO(CFILEDIA,CLUOUTDIA,CGROUP)
        IF (PRESENT(KZ)) THEN
        DO J2=1,SIZE(ZTAB3,2)
        DO J1=1,SIZE(ZTAB3,1)
          XVAR(J1,J2,:,K,1,1)=ZTAB3(J1,J2,1:IKU:KZ)
        END DO
        END DO
        ELSE
        XVAR(:,:,:,K,1,1)=ZTAB3
        ENDIF
        CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH, &
          TDTCUR%TDATE%DAY,TDTCUR%TIME,TDTEXP%TDATE%YEAR,    &
          TDTEXP%TDATE%MONTH,TDTEXP%TDATE%DAY,TDTEXP%TIME,           &
          XTRAJT(K,1))
        XDATIME(1,K)=TDTEXP%TDATE%YEAR; XDATIME(2,K)=TDTEXP%TDATE%MONTH
        XDATIME(3,K)=TDTEXP%TDATE%DAY;  XDATIME(4,K)=TDTEXP%TIME
        XDATIME(5,K)=TDTSEG%TDATE%YEAR; XDATIME(6,K)=TDTSEG%TDATE%MONTH
        XDATIME(7,K)=TDTSEG%TDATE%DAY;  XDATIME(8,K)=TDTSEG%TIME
        XDATIME(9,K)=TDTMOD%TDATE%YEAR; XDATIME(10,K)=TDTMOD%TDATE%MONTH
        XDATIME(11,K)=TDTMOD%TDATE%DAY; XDATIME(12,K)=TDTMOD%TIME
        XDATIME(13,K)=TDTCUR%TDATE%YEAR;XDATIME(14,K)=TDTCUR%TDATE%MONTH
        XDATIME(15,K)=TDTCUR%TDATE%DAY; XDATIME(16,K)=TDTCUR%TIME
  
        WRITE(ILUOUTDIA,*)' OTHERSFIELDS IGRID XVAR,XTRAJT,CTITRE,CUNITE,CCOMMENT'
        WRITE(ILUOUTDIA,*)IGRID,SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
        SIZE(XVAR,5),SIZE(XVAR,6),'  ',SIZE(XTRAJT,1),SIZE(XTRAJT,2)
        WRITE(ILUOUTDIA,*)(CTITRE(J)(1:LEN(CTITRE)),J=1,SIZE(CTITRE))
        WRITE(ILUOUTDIA,*)(CUNITE(J)(1:LEN(CUNITE)),J=1,SIZE(CUNITE))
        WRITE(ILUOUTDIA,*)(CCOMMENT(J)(1:LEN(CCOMMENT)),J=1,SIZE(CCOMMENT))
!
! Ecriture dans le fichier diachronique
        NGRIDIA(1)=IGRID
        CALL WRITE_DIACHRO(HFILEDIA,HLUOUTDIA,CGROUP,CTYPE,NGRIDIA,XDATIME,XVAR, &
        XTRAJT,CTITRE,CUNITE,CCOMMENT, &
        LICP,LJCP,LKCP,NIL,NIH,NJL,NJH,NKL,NKH)
!
! Desallocation des matrices
        DEALLOCATE(ZTAB3)
        CALL ALLOC_FORDIACHRO(1,1,IKU,NNBF,1,1,3)
!
      ENDIF                             !************************************
!
!
!*       2.4  ++++0D+++++
!
! Traitement des scalaires 'individuels'
    ELSE IF(NSIZT(JJ,K) == 1)THEN
!     WRITE(NLUOUTD,*)'***************************************************************'
!     WRITE(NLUOUTD,*)' Scalaire rencontre et non enregistre dans le fichier',&
!     WRITE(NLUOUTD,*)' Scalaire rencontre et enregistre dans le fichier',&
!     &' diachronique ',CGROUP,' size ',NSIZT(JJ,K)
!     WRITE(NLUOUTD,*)' Prevenir J.DURON . Mail: durj@aero.obs-mip.fr '
!     Prise en compte de certains temps 
!     WRITE(NLUOUTD,*)'***************************************************************'
      IPCENT=0
      IPCENT=INDEX(CRECFM2T(JJ,K),'%TIM')
      IF(IPCENT /= 0)THEN                      !===================
        CALL FMWRIT(HFILEDIA,CGROUP,HLUOUTDIA,NSIZT(JJ,K),ZTAB,IGRID,&
        ILYCOMM,YCOMM,IRESP)
!       ILENCH,YCOMMENT,IRESP)
        CALL ELIM(CRECFM2T(JJ,K))
        print *,' Impression pour controle ',CGROUP,ZTAB,' size ', &
        NSIZT(JJ,K)
      ELSE                                 !===================
!
        ITAB3(1)=1; ITAB3(2)=1; ITAB3(3)=1
        IF(ALLOCATED(ZTAB3))THEN
          DEALLOCATE(ZTAB3)
        ENDIF
        ALLOCATE(ZTAB3(ITAB3(1),ITAB3(2),ITAB3(3)))
        ZTAB3=RESHAPE(ZTAB,ITAB3)
!
        IF(K == 1)THEN      
          CTYPE='CART'
!
! Allocation des matrices utilisees dans le fichier diachronique (dernier
! argument = 1 pour ecriture; = 2 pour lecture; si =3, desallocation)
!
          CALL ALLOC_FORDIACHRO(ITAB3(1),ITAB3(2),ITAB3(3),NNBF,1,1,1)

! Initialisation de variables et matrices
          LICP=.FALSE. ; LJCP=.FALSE. ; LKCP=.FALSE.
          NIL=1 ; NJL=1 ; NKL=1
          NIH=1 ; NJH=1 ; NKH=1
          XVAR(:,:,:,:,:,:)=0.
          XTRAJT(:,:)=0.
          CTITRE(:)(1:LEN(CTITRE))=' '
          CUNITE(:)(1:LEN(CUNITE))=' '
          CCOMMENT(:)(1:LEN(CCOMMENT))=' '
          XDATIME(:,:)=0
        ENDIF
!
! Distinction 1er fichier et les suivants. Dans le premier cas on ecrit di-
! -rectement dans le fic. diachronique et apres les avoir reorganisees les
! informations lues. Dans les cas suivants, on relit d'abord les infos du
! fic. diachron. pour les augmenter des nouvelles  fraichement lues avant
! de les reecrire.
! NOTA on a pris la precaution de prevoir des le depart une taille d'article
! = a la dimension de la matrice traitee * par le nb de fichiers lus (NNBF)
!
        IF (K == 1)THEN                   !************************************
          XVAR(:,:,:,K,1,1)=ZTAB3
!
! Le tps courant est transforme en temps relatif par / au debut de l'experience
          CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH, &
            TDTCUR%TDATE%DAY,TDTCUR%TIME,TDTEXP%TDATE%YEAR,    &
            TDTEXP%TDATE%MONTH,TDTEXP%TDATE%DAY,TDTEXP%TIME,           &
            XTRAJT(K,1))
          CTITRE(1)=CGROUP
          CUNITE(1)=ADJUSTL(YCAROUT)
          CCOMMENT(1)=YCOMMENT
          XDATIME(1,K)=TDTEXP%TDATE%YEAR; XDATIME(2,K)=TDTEXP%TDATE%MONTH
          XDATIME(3,K)=TDTEXP%TDATE%DAY;  XDATIME(4,K)=TDTEXP%TIME
          XDATIME(5,K)=TDTSEG%TDATE%YEAR; XDATIME(6,K)=TDTSEG%TDATE%MONTH
          XDATIME(7,K)=TDTSEG%TDATE%DAY;  XDATIME(8,K)=TDTSEG%TIME
          XDATIME(9,K)=TDTMOD%TDATE%YEAR; XDATIME(10,K)=TDTMOD%TDATE%MONTH
          XDATIME(11,K)=TDTMOD%TDATE%DAY; XDATIME(12,K)=TDTMOD%TIME
          XDATIME(13,K)=TDTCUR%TDATE%YEAR;XDATIME(14,K)=TDTCUR%TDATE%MONTH
          XDATIME(15,K)=TDTCUR%TDATE%DAY; XDATIME(16,K)=TDTCUR%TIME
!
! Ecriture dans le fichier diachronique
          NGRIDIA(1)=IGRID
          CALL WRITE_DIACHRO(HFILEDIA,HLUOUTDIA,CGROUP,CTYPE,NGRIDIA,XDATIME,XVAR, &
            XTRAJT,CTITRE,CUNITE,CCOMMENT, &
            LICP,LJCP,LKCP,NIL,NIH,NJL,NJH,NKL,NKH)
!
! Desallocation des matrices
          DEALLOCATE(ZTAB3)
          CALL ALLOC_FORDIACHRO(1,1,1,NNBF,1,1,3)
!
        ELSE                              !************************************
!
! On relit les infos deja enregistrees du fichier diachronique en connaissant
! le nom du groupe CGROUP=CRECFM2T(JJ,K)
          CALL READ_DIACHRO(CFILEDIA,CLUOUTDIA,CGROUP)
          XVAR(:,:,:,K,1,1)=ZTAB3
          CALL TEMPORAL_DIST(TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH, &
            TDTCUR%TDATE%DAY,TDTCUR%TIME,TDTEXP%TDATE%YEAR,    &
            TDTEXP%TDATE%MONTH,TDTEXP%TDATE%DAY,TDTEXP%TIME,           &
            XTRAJT(K,1))
          XDATIME(1,K)=TDTEXP%TDATE%YEAR; XDATIME(2,K)=TDTEXP%TDATE%MONTH
          XDATIME(3,K)=TDTEXP%TDATE%DAY;  XDATIME(4,K)=TDTEXP%TIME
          XDATIME(5,K)=TDTSEG%TDATE%YEAR; XDATIME(6,K)=TDTSEG%TDATE%MONTH
          XDATIME(7,K)=TDTSEG%TDATE%DAY;  XDATIME(8,K)=TDTSEG%TIME
          XDATIME(9,K)=TDTMOD%TDATE%YEAR; XDATIME(10,K)=TDTMOD%TDATE%MONTH
          XDATIME(11,K)=TDTMOD%TDATE%DAY; XDATIME(12,K)=TDTMOD%TIME
          XDATIME(13,K)=TDTCUR%TDATE%YEAR;XDATIME(14,K)=TDTCUR%TDATE%MONTH
          XDATIME(15,K)=TDTCUR%TDATE%DAY; XDATIME(16,K)=TDTCUR%TIME

          WRITE(ILUOUTDIA,*)' OTHERSFIELDS IGRID XVAR,XTRAJT,CTITRE,CUNITE,CCOMMENT'
          WRITE(ILUOUTDIA,*)IGRID,SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
          SIZE(XVAR,5),SIZE(XVAR,6),'  ',SIZE(XTRAJT,1),SIZE(XTRAJT,2)
          WRITE(ILUOUTDIA,*)(CTITRE(J)(1:LEN(CTITRE)),J=1,SIZE(CTITRE))
          WRITE(ILUOUTDIA,*)(CUNITE(J)(1:LEN(CUNITE)),J=1,SIZE(CUNITE))
          WRITE(ILUOUTDIA,*)(CCOMMENT(J)(1:LEN(CCOMMENT)),J=1,SIZE(CCOMMENT))
!
! Ecriture dans le fichier diachronique
          NGRIDIA(1)=IGRID
          CALL WRITE_DIACHRO(HFILEDIA,HLUOUTDIA,CGROUP,CTYPE,NGRIDIA,XDATIME,XVAR, &
          XTRAJT,CTITRE,CUNITE,CCOMMENT, &
          LICP,LJCP,LKCP,NIL,NIH,NJL,NJH,NKL,NKH)
!
! Desallocation des matrices
          DEALLOCATE(ZTAB3)
          CALL ALLOC_FORDIACHRO(1,1,1,NNBF,1,1,3)

        ENDIF                             !************************************
      ENDIF                             !===============
!
!
!*       2.5  ++++  +++++
!
! Impression des infos non recensees
!
    ELSE
!     Some dates are taken into account
      IPCENT=0
      IPCENT=INDEX(CRECFM2T(JJ,K),'%TDA')
      IF(IPCENT /= 0)THEN                      !===================
        CALL FMWRIT(HFILEDIA,CGROUP,HLUOUTDIA,NSIZT(JJ,K),ZTAB,IGRID,&
        ILYCOMM,YCOMM,IRESP)
!       ILENCH,YCOMMENT,IRESP)
        CALL ELIM(CRECFM2T(JJ,K))
        print *,' Impression pour controle ',CGROUP,ZTAB,' size ',NSIZT(JJ,K)
      ELSE
        WRITE(NLUOUTD,*)'***************************************************************'
        WRITE(NLUOUTD,*)' Variable non prise en compte dans le fichier diachronique ',&
        CGROUP,' size ',NSIZT(JJ,K),' IIU IJU IKU ',IIU,IJU,IKU
        IF (LEN_TRIM(YCOMMENT) /=0) THEN
          WRITE(NLUOUTD,*)' YCOMMENT=',YCOMMENT(1:LEN_TRIM(YCOMMENT))
        ELSE
          WRITE(NLUOUTD,*)' YCOMMENT '
        ENDIF
        WRITE(NLUOUTD,*)'***************************************************************'
      ENDIF
    ENDIF
!
!
!*       2.6  ++++END+++++
!
!
    DEALLOCATE(ZTAB)
    IF(K == NNBF)THEN
      WRITE(ILUOUTDIA,*)CRECFM2T(JJ,K),' TREATED with size ', NSIZT(JJ,K)*K*IMULT
    ENDIF
!
!
!----------------------------------------------------------------------------
!
!*       3.    TREATMENT OF ELIMINATED VARIABLE
!              --------------------------------
!
  ELSE
    IPCENT=0
    IPCENT=INDEX(CRECFM2T(JJ,K),'%TIM')
    IF(IPCENT /= 0 .AND. K >1)THEN   
      IF(INDEX(CRECFM2T(JJ,K),'TDTEXP%TDA') /= 0 .OR.      &
         INDEX(CRECFM2T(JJ,K),'TDTEXP%TIM') /= 0 .OR.      &
         INDEX(CRECFM2T(JJ,K),'TDTSEG%TDA') /= 0 .OR.      &
         INDEX(CRECFM2T(JJ,K),'TDTSEG%TIM') /= 0 .OR.      &
         INDEX(CRECFM2T(JJ,K),'TDTMOD%TDA') /= 0 .OR.      &
         INDEX(CRECFM2T(JJ,K),'TDTMOD%TIM') /= 0 .OR.      &
         INDEX(CRECFM2T(JJ,K),'TDTCUR%TDA') /= 0 .OR.      &
         INDEX(CRECFM2T(JJ,K),'TDTCUR%TIM') /= 0)THEN
      ELSE
        ALLOCATE(ZTAB(NSIZT(JJ,K)))
        CALL FMREAD(CNAMFILED(K),CRECFM2T(JJ,K),CLUOUTD,NSIZT(JJ,K), &
        ZTAB,IGRID,ILENCH,YCOMMENT,IRESP)
        print *,' CRECFM2T(JJ,K)  K= ',CRECFM2T(JJ,K),K,'  non enr. volontairement .'
        DEALLOCATE(ZTAB)
      ENDIF
    ENDIF
  ENDIF
!
ENDDO
!
LPACK=GPACK
!----------------------------------------------------------------------------
RETURN
!
END SUBROUTINE WRITE_OTHERSFIELDS
