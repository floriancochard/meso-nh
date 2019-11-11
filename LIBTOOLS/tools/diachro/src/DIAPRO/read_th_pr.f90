!     ##################################################
      SUBROUTINE READ_TH_PR(HFILEDIA,HLUOUTDIA,KMT,KIND)
!     ##################################################
!
!!****  *READ_TH_PR* - 
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
!!      Original       08/01/97
!!      Updated   PM 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ALLOC_FORDIACHRO
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_SEVERAL_RECORDS
USE MODD_RESOLVCAR
USE MODD_FILES_DIACHRO
USE MODD_MASK3D

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER :: KMT, KIND
CHARACTER(LEN=*) :: HFILEDIA, HLUOUTDIA
!
!*       0.1   Local variables
!              ---------------

!
INTEGER   ::   J 
CHARACTER(LEN=12) :: YGP, YGPM
!------------------------------------------------------------------------------
!
! KIND=1 --> LTK=.TRUE. or LEV=.TRUE.
!
YGP='      '
YGPM='      '
IF(KIND == 1)THEN
  IF(KMT == 1)THEN
    IF(LTK .OR. LRS .OR. LRS1)THEN
      YGP='THM'
    ELSE IF(LEV)THEN
      YGP='POVOM'
    ELSE IF(LSV3)THEN
      IF(LXYZ00)THEN
        YGP=CGROUPSV3(1:LEN_TRIM(CGROUPSV3))
      ELSE
        YGP='LGZM'
      ENDIF
      YGPM=YGP
      CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,YGP)
      IF(LPBREAD .AND. .NOT.LXYZ00)THEN
        LPBREAD=.FALSE.
        YGP='SVM003'
        CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,YGP)
        IF(LPBREAD)THEN
          LPBREAD=.FALSE.
          YGP='SVM3'
        ENDIF
      ENDIF
    ENDIF
  ELSE IF(KMT == 2)THEN
    IF(LTK .OR. LRS .OR. LRS1)THEN
      YGP='THT'
    ELSE IF(LEV)THEN
      YGP='POVOT'
    ELSE IF(LSV3)THEN
      IF(LXYZ00)THEN
        YGP=CGROUPSV3(1:LEN_TRIM(CGROUPSV3))
      ELSE
        YGP='LGZT'
      ENDIF
      YGPM=YGP
      CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,YGP)
      IF(LPBREAD .AND. .NOT.LXYZ00)THEN
        LPBREAD=.FALSE.
        YGP='SVT003'
        CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,YGP)
        IF(LPBREAD)THEN
          LPBREAD=.FALSE.
          YGP='SVT3'
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  SELECT CASE(KMT)
    CASE(1)
      CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,YGP)
      IF(LPBREAD)THEN
!       LPBREAD=.FALSE.
        IF(LSV3)THEN
!         IF(.NOT.LXY00)THEN
	    IF(YGP /= YGPM)THEN
	      IF(INDEX(YGP,'00') == 0)THEN
                print *,' **READ_TH_PR requete peut-etre impossible.', YGPM, &
                       ', ',YGP(1:3)//'00'//YGP(4:4),' et ',YGP,' n''existent pas'
              ELSE
                print *,' **READ_TH_PR requete peut-etre impossible.',YGPM, &
                        ' et ',YGP,' n''existent pas'
              ENDIF
            ENDIF
!         ENDIF
        ELSE
          print *,' REQUETE IMPOSSIBLE .',YGP,' N''EXISTE PAS'
        ENDIF
        IF(.NOT.LSV3)THEN
          YGP(LEN_TRIM(YGP):LEN_TRIM(YGP))='T'
          print *,' **READ_TH_PR  Recherche de  ** ',YGP,' ** pour resoudre le pb'
        ENDIF
        RETURN
      ELSE
        print *,' **READ_TH_PR Utilisation de   ** ',YGP,' **'
      ENDIF
      IF(LGROUP)THEN
        CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,YGP)
      ENDIF
      IF(.NOT.LFIC1)THEN
        CALL REALLOC_AND_LOAD(YGP)
        IF(LPBREAD)THEN
!         LPBREAD=.FALSE.
          print *,' REQUETE IMPOSSIBLE .',YGP,' N''EXISTE PAS DANS', &
          ' L''UN DES FICHIERS '
	  IF(ALLOCATED(XVAR))THEN
	    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	  ENDIF
          RETURN
        ENDIF
      ENDIF
    CASE(2)
      CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,YGP)
      IF(LPBREAD)THEN
!       LPBREAD=.FALSE.
        IF(LSV3)THEN
!         IF(.NOT.LXY00)THEN
            IF(YGP /= YGPM)THEN
	      IF(INDEX(YGP,'00') == 0)THEN
                print *,' **READ_TH_PR requete peut-etre impossible. ',YGPM, &
                       ', ',YGP(1:3)//'00'//YGP(4:4),' et ',YGP,' n''existent pas'
              ELSE
                print *,' **READ_TH_PR requete peut-etre impossible. ',YGPM, &
                        ' et ',YGP,' n''existent pas'
              ENDIF
	    ENDIF
!         ENDIF
	ELSE
          print *,' REQUETE IMPOSSIBLE .',YGP,' N''EXISTE PAS'
        ENDIF
        IF(.NOT.LSV3)THEN
          YGP(LEN_TRIM(YGP):LEN_TRIM(YGP))='M'
	  print *,' **READ_TH_PR  Recherche de   ** ',YGP,' ** pour resoudre le pb'
        ENDIF
        RETURN
      ELSE
        print *,' **READ_TH_PR  Utilisation de   ** ',YGP,' **'
      ENDIF
      IF(LGROUP)THEN
        CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,YGP)
      ENDIF
      IF(.NOT.LFIC1)THEN
        CALL REALLOC_AND_LOAD(YGP)
        IF(LPBREAD)THEN
!         LPBREAD=.FALSE.
          print *,' REQUETE IMPOSSIBLE .',YGP,' N''EXISTE PAS DANS', &
          ' L''UN DES FICHIERS '
          IF(ALLOCATED(XVAR))THEN
            CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
          ENDIF
          RETURN
        ENDIF
      ENDIF
  END SELECT
  IF(ALLOCATED(XTH)) DEALLOCATE(XTH)
  ALLOCATE(XTH(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
  SIZE(XVAR,5),SIZE(XVAR,6)))
  XTH(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)
  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)

! KIND=2 --> LPR=.TRUE.
ELSE IF(KIND == 2)THEN

  SELECT CASE(KMT)
    CASE(1)
      CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,'PABSM')
!     CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,'PHIM')
      IF(LPBREAD)THEN
!       LPBREAD=.FALSE.
        print *,' REQUETE a priori IMPOSSIBLE . PABSM N''EXISTE PAS . '
	print *,' **READ_TH_PR  Recherche de  **  PABST  ** pour resoudre le pb'
        RETURN
      ELSE
        print *,' **READ_TH_PR  Utilisation de  ** PABSM **'
      ENDIF
      IF(LGROUP)THEN
        CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,'PABSM')
!       CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,'PHIM')
      ENDIF
      IF(.NOT.LFIC1)THEN
        CALL REALLOC_AND_LOAD('PABSM')
!       CALL REALLOC_AND_LOAD('PHIM')
        IF(LPBREAD)THEN
!         LPBREAD=.FALSE.
          print *,' REQUETE IMPOSSIBLE . PABSM N''EXISTE PAS DANS', &
          ' L''UN DES FICHIERS '
	  IF(ALLOCATED(XVAR))THEN
	    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	  ENDIF
          RETURN
        ENDIF
      ENDIF
    CASE(2)
      CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,'PABST')
!     CALL VERIF_GROUP(HFILEDIA,HLUOUTDIA,'PHIT')
      IF(LPBREAD)THEN
!       LPBREAD=.FALSE.
        print *,' REQUETE a priori IMPOSSIBLE . PABST N''EXISTE PAS . '
	print *,' **READ_TH_PR  Recherche de  **  PABSM  ** pour resoudre le pb'
         RETURN
       ELSE
         print *,' **READ_TH_PR  Utilisation de   ** PABST **'
       ENDIF
       IF(LGROUP)THEN
         CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,'PABST')
!        CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,'PHIT')
       ENDIF
       IF(.NOT.LFIC1)THEN
         CALL REALLOC_AND_LOAD('PABST')
!        CALL REALLOC_AND_LOAD('PHIT')
         IF(LPBREAD)THEN
!          LPBREAD=.FALSE.
           print *,' REQUETE IMPOSSIBLE . PABST N''EXISTE PAS DANS', &
           ' L''UN DES FICHIERS '
           IF(ALLOCATED(XVAR))THEN
             CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
           ENDIF
           RETURN
         ENDIF
       ENDIF
    END SELECT
    ALLOCATE(XPHI(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
        SIZE(XVAR,5),SIZE(XVAR,6)))
    XPHI(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)
    IF(.NOT.LRS .AND. .NOT.LRS1)THEN
      CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
    ENDIF
    IF(ALLOCATED(XPRES))THEN
      DEALLOCATE(XPRES)
    ENDIF
    ALLOCATE(XPRES(SIZE(XPHI,1),SIZE(XPHI,2),SIZE(XPHI,3),SIZE(XPHI,4), &
      SIZE(XPHI,5),SIZE(XPHI,6)))
    IF(SIZE(XPHI,5) /= 1 .OR. SIZE(XPHI,6) /= 1)THEN
      print *,' SIZE(XPHI,5) SIZE(XPHI,6) /= 1 ',SIZE(XPHI,5),SIZE(XPHI,6)
      print *,' CALCUL DE LA PRESSION IMPOSSIBLE. REQUETE NON TRAITEE '
      DEALLOCATE(XPHI,XPRES)
      LPBREAD=.TRUE.
      RETURN
    ENDIF
!!  Calcul de la pres/sion
!   Chargement de la pression
    DO J=1,SIZE(XPHI,4)
!     XPRES(:,:,:,J,1,1)=XP00*(XEXNREF(:,:,:)+XPHI(:,:,:,J,1,1) &
!                        /(XCPD*XTHVREF(:,:,:)))**(XCPD/XRD)
      XPRES(:,:,:,J,1,1)=XPHI(:,:,:,J,1,1)
    ENDDO
    DEALLOCATE(XPHI)
ENDIF
!
!-----------------------------------------------------------------------------
!
!*       2.       RETURNS
!                 -----
! 
RETURN
END SUBROUTINE READ_TH_PR
