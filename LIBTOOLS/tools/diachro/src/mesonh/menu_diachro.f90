!     ######spl
      MODULE MODI_MENU_DIACHRO
!     #########################
!
INTERFACE
!
SUBROUTINE MENU_DIACHRO(HFILEDIA,HLUOUTDIA,HGROUP)
CHARACTER(LEN=*) :: HGROUP
CHARACTER(LEN=*) :: HFILEDIA,HLUOUTDIA
END SUBROUTINE MENU_DIACHRO
!
END INTERFACE
!
END MODULE MODI_MENU_DIACHRO
!     ##################################################
      SUBROUTINE MENU_DIACHRO(HFILEDIA,HLUOUTDIA,HGROUP)
!     ##################################################
!
!!****  *MENU_DIACHRO* - Creation, ecriture (eventuellement lecture) de
!          l'enregistrement MENU_BUDGET  dans un fichier diachronique
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!     
!      A chaque ecriture d'un enregistrement dans un fichier diachronique,
!      cette routine est appelee pour memoriser le nom du groupe correspon-
!      -dant (passe en argument dans HGROUP)
!     Au terme des ecritures, elle est appelee avec HGROUP='END' qui
!     a pour effet d'ecrire dans le fichier diachronique le tableau contenant
!     le nom des groupes avec l'identificateur de record : MENU_BUDGET
!     Quand HGROUP='READ', l'enregistrement MENU_BUDGET est lu et la
!     liste des groupes enregistres est imprimee
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
!!      Original       08/01/96
!!      Updated   PM 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!

USE MODD_OUT_DIA
USE MODI_FMREAD 
USE MODI_FMWRIT 

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HGROUP
CHARACTER(LEN=*) :: HFILEDIA, HLUOUTDIA
!
!*       0.1   Local variables
!              ---------------

!
CHARACTER(LEN=16) :: YRECFM
CHARACTER(LEN=20) :: YCOMMENT
CHARACTER(LEN=16),DIMENSION(2000),SAVE    :: YGROUP 
!CHARACTER(LEN=16),DIMENSION(5000),SAVE    :: YGROUP  ! dans le conv2dia.select
INTEGER   ::   ILENG, ILENCH, IGRID, J, JJ, ILENDIM, IALREADY
INTEGER   ::   IRESPDIA
INTEGER,SAVE   ::   IGROUP=0
INTEGER,DIMENSION(:),ALLOCATABLE :: ITABCHAR
!------------------------------------------------------------------------------
!
IF(HGROUP == 'END')THEN

  IF(IGROUP == 0)THEN
    print *,' No record for the diachronic file'
    RETURN
  ENDIF
  IGRID=0
  ILENDIM=1
  ILENG=16*IGROUP
  ILENCH=LEN(YCOMMENT)
  YRECFM='MENU_BUDGET.DIM'
  CALL FMWRIT(HFILEDIA,YRECFM,HLUOUTDIA,ILENDIM,ILENG,&
  IGRID,ILENCH,YCOMMENT,IRESPDIA)

  YRECFM='MENU_BUDGET'
  ALLOCATE(ITABCHAR(ILENG))
  DO JJ=1,IGROUP
    DO J = 1,16
      ITABCHAR(16*(JJ-1)+J) = ICHAR(YGROUP(JJ)(J:J))
    ENDDO
  ENDDO
  CALL FMWRIT(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
  IGRID,ILENCH,YCOMMENT,IRESPDIA)
  DEALLOCATE(ITABCHAR)

ELSE IF(HGROUP == 'READ')THEN

  ILENDIM=1
  YRECFM='MENU_BUDGET.DIM'
  CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENDIM,ILENG,&
  IGRID,ILENCH,YCOMMENT,IRESPDIA)
  IF(IRESPDIA == -47)THEN
    print *,' No record MENU_BUDGET '
    RETURN
  ENDIF

  ALLOCATE(ITABCHAR(ILENG))
  YRECFM='MENU_BUDGET'
  CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
  IGRID,ILENCH,YCOMMENT,IRESPDIA)
  IGROUP=ILENG/16
  DO JJ=1,IGROUP
    DO J = 1,16
      YGROUP(JJ)(J:J)=CHAR(ITABCHAR(16*(JJ-1)+J))
    ENDDO
  ENDDO
  DO JJ=1,IGROUP
    WRITE(NLUOUTD,*)' ******** YGROUP :  ',YGROUP(JJ)
    !print *,' ******** YGROUP :  ',YGROUP(JJ)
  ENDDO
  print *,'****************************** GROUPS *****************************'
  print 100,(YGROUP(JJ),JJ=1,IGROUP)
100 FORMAT(1X,5A15)
  DEALLOCATE(ITABCHAR)

ELSE

  IALREADY=0
  IF(IGROUP > 1)THEN
    DO JJ=1,IGROUP
      IF(ADJUSTL(HGROUP) == YGROUP(JJ))IALREADY=1
    ENDDO
  ENDIF
  IF(IALREADY == 0)THEN
    IGROUP=IGROUP+1
    YGROUP(IGROUP)=ADJUSTL(HGROUP)
  ENDIF
ENDIF
!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE MENU_DIACHRO
