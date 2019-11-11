!     ######spl
      MODULE MODI_VERIF_GROUP
!     #######################
!
INTERFACE
!
SUBROUTINE VERIF_GROUP(HFILEDIA,HLUOUTDIA,HGROUP)
CHARACTER(LEN=*) :: HFILEDIA, HLUOUTDIA, HGROUP
END SUBROUTINE VERIF_GROUP
!
END INTERFACE
END MODULE MODI_VERIF_GROUP
!     ######spl
      SUBROUTINE VERIF_GROUP(HFILEDIA,HLUOUTDIA,HGROUP)
!     #################################################
!
!!****  *VERIF_GROUP* - 
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
!!      Original       08/01/96
!!      Updated   PM 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIACHRO
USE MODD_TYPE_AND_LH
USE MODD_RESOLVCAR
USE MODD_SEVERAL_RECORDS
USE MODN_NCAR
USE MODD_ALLOC_FORDIACHRO
USE MODI_REALLOC_AND_LOAD_RECORDS
USE MODI_FMREAD

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HFILEDIA, HLUOUTDIA,HGROUP
!
!*       0.1   Local variables
!              ---------------

!
CHARACTER(LEN=16) :: YRECFM
CHARACTER(LEN=8)  :: YNAM1, YNAM2, YNAM1M, YNAM2M
! Aout 99 Longueur YCOMMENT passee de 20 a 100
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER*1       :: Y1
CHARACTER*2       :: Y2
CHARACTER*3       :: Y3
CHARACTER*4       :: Y4
CHARACTER(LEN=16),DIMENSION(:),ALLOCATABLE:: YGROUP 
INTEGER   ::   ILENG, ILENCH, IGRID, J, JJ, JM, ILENDIM
INTEGER   ::   JM1, JM2, INCR1, INCR2
INTEGER   ::   IRESPDIA
INTEGER   ::   IMINUS, ILENGP, INBC2, INBC1
INTEGER,SAVE   ::   IGROUP=0
INTEGER,DIMENSION(:),ALLOCATABLE :: ITABCHAR
LOGICAL   ::   GPART
!------------------------------------------------------------------------------
!

GPART=.FALSE.
NBCNUM=0
NINCRNAM=1
CGPNAM(1:LEN(CGPNAM))=' '
CGPNAM1(1:LEN(CGPNAM1))=' '
CGPNAM2(1:LEN(CGPNAM2))=' '
YNAM1(1:LEN(YNAM1))=' '
YNAM2(1:LEN(YNAM2))=' '
YNAM1M(1:LEN(YNAM1M))=' '
YNAM2M(1:LEN(YNAM2M))=' '
print *,' VERIF_GROUP HGROUP ',HGROUP

ILENDIM=1
YRECFM='MENU_BUDGET.DIM'
CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENDIM,ILENG,&
IGRID,ILENCH,YCOMMENT,IRESPDIA)

IF(ALLOCATED(ITABCHAR))DEALLOCATE(ITABCHAR)
ALLOCATE(ITABCHAR(ILENG))
YRECFM='MENU_BUDGET'
CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
IGRID,ILENCH,YCOMMENT,IRESPDIA)
IGROUP=ILENG/16
IF(ALLOCATED(YGROUP))DEALLOCATE(YGROUP)
ALLOCATE(YGROUP(IGROUP))
print *,' ILENG ILENCH IGROUP ',ILENG,ILENCH,IGROUP

DO JJ=1,IGROUP
  DO J = 1,16
    YGROUP(JJ)(J:J)=CHAR(ITABCHAR(16*(JJ-1)+J))
  ENDDO
ENDDO
DEALLOCATE(ITABCHAR)
YRECFM=ADJUSTL(ADJUSTR(HGROUP)//'.TYPE')
ILENG=LEN(CTYPE)
ALLOCATE(ITABCHAR(ILENG))
CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
IGRID,ILENCH,YCOMMENT,IRESPDIA)
!******************************************************************************

IF(IRESPDIA == 0)THEN
!*************  A DEFINIR **********************
  LGROUP=.TRUE.
  RETURN
!******************************************************************************

ELSE IF(IRESPDIA == -47)THEN

  LGROUP=.FALSE.

! On decortique HGROUP

  ILENGP=LEN_TRIM(HGROUP)

!---------------------------------------------------
  IF(HGROUP(1:ILENGP) == 'PABSM' .OR. HGROUP(1:ILENGP) == 'PABST' .OR. &
     HGROUP(1:ILENGP) == 'THM'   .OR. HGROUP(1:ILENGP) == 'THT'   .OR. &
     HGROUP(1:ILENGP) == 'POVOM' .OR. HGROUP(1:ILENGP) == 'POVOT' .OR. &
     HGROUP(1:ILENGP) == 'SVM3' .OR. HGROUP(1:ILENGP) == 'SVM003' .OR. &
     HGROUP(1:ILENGP) == 'SVT3' .OR. HGROUP(1:ILENGP) == 'SVT003' .OR. &
     HGROUP(1:ILENGP) == 'LGZM' .OR. HGROUP(1:ILENGP) == 'LGZT'   )THEN
!   print *,' VERIF_GROUP PAS OK 1',HGROUP
     LPBREAD=.TRUE.
     RETURN
  ENDIF
!---------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(ILENGP > 1)THEN
    IF(ILENGP <= 4 ) THEN
      IF( HGROUP(ILENGP:ILENGP) == '0' .OR. HGROUP(ILENGP:ILENGP) == '1' .OR. &
          HGROUP(ILENGP:ILENGP) == '2' .OR. HGROUP(ILENGP:ILENGP) == '4' .OR. &
          HGROUP(ILENGP:ILENGP) == '5' .OR. HGROUP(ILENGP:ILENGP) == '6' .OR. &
          HGROUP(ILENGP:ILENGP) == '7' .OR. HGROUP(ILENGP:ILENGP) == '8' .OR. &
          HGROUP(ILENGP:ILENGP) == '9') THEN
          IF (HGROUP(1:2) == 'UM' .OR. HGROUP(1:2) == 'VM' .OR.&
              HGROUP(1:2) == 'WM' .OR. HGROUP(1:2) == 'UT' .OR.&
              HGROUP(1:2) == 'VT' .OR. HGROUP(1:2) == 'WT') THEN
                LPBREAD=.TRUE.
   !             print *,' VERIF_GROUP PAS OK 2',HGROUP
               RETURN
          ENDIF
      ENDIF
    ENDIF
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recherche d'un signe - a partir de la fin

  DO J=0,4
    IMINUS=INDEX(HGROUP(MAX(ILENGP-J,1):ILENGP),'-')
    IF(IMINUS /= 0)THEN
      JM=J
      EXIT
    ENDIF
  ENDDO

! Presence d'un signe moins

  IF(IMINUS /= 0)THEN

! Cas expression groupe sous la forme AA__0001-0099 (Donc LFIC1=.TRUE.)  ou
! sous la forme AA_b-c-

    IMINUS=ILENGP-JM+IMINUS-1

    IF(IMINUS == ILENGP)THEN     !00000000000000000000000000000000000000
! Pas d'intervalle mais presence d'un ou plusieurs signes -

      GPART=.TRUE.

    ELSE              !0000000000000000000000000000000000000

! Intervalle poossible

      JM1=0; JM2=0; INCR1=0; INCR2=0
      J=IMINUS-1 ;  JJ=IMINUS+1
      IF((HGROUP(J:J) == '0' .OR. HGROUP(J:J) == '1' .OR. HGROUP(J:J) == '2' &
      .OR. HGROUP(J:J) == '3' .OR. HGROUP(J:J) == '4' .OR. HGROUP(J:J) == '5' &
      .OR. HGROUP(J:J) == '6' .OR. HGROUP(J:J) == '7' .OR. HGROUP(J:J) == '8' &
      .OR. HGROUP(J:J) == '9') .AND.                     &
       (HGROUP(JJ:JJ) == '0' .OR. HGROUP(JJ:JJ) =='1' .OR. HGROUP(JJ:JJ) == '2' &
       .OR. HGROUP(JJ:JJ)=='3' .OR. HGROUP(JJ:JJ)=='4' .OR. HGROUP(JJ:JJ) == '5' &
       .OR. HGROUP(JJ:JJ)=='6' .OR. HGROUP(JJ:JJ)=='7' .OR. HGROUP(JJ:JJ) == '8' &
       .OR. HGROUP(JJ:JJ) == '9'))THEN
      
      INBC2=ILENGP-IMINUS
      READ(HGROUP(IMINUS+1:ILENGP),*)NAM2
      JM=0
      DO J=2,IMINUS-1
        IF(HGROUP(J:J) == '0' .OR. HGROUP(J:J) == '1' .OR. HGROUP(J:J) == '2'  &
        .OR. HGROUP(J:J) == '3' .OR. HGROUP(J:J) == '4' .OR. HGROUP(J:J) == '5' &
        .OR. HGROUP(J:J) == '6' .OR. HGROUP(J:J) == '7' .OR. HGROUP(J:J) == '8' &
        .OR. HGROUP(J:J) == '9')THEN
        JM=J
        EXIT
	ENDIF
      ENDDO

	INBC1=IMINUS-JM
! On memorise les infos pour realloc_several_records
	READ(HGROUP(JM:IMINUS-1),*)NAM1
	IF(INBC1-INBC2 == 0)NBCNUM=INBC1
	CGPNAM=HGROUP(1:JM-1)
	CGPNAM=ADJUSTL(CGPNAM)
	CGPNAM1=HGROUP(1:IMINUS-1)
	CGPNAM1=ADJUSTL(CGPNAM1)
	CGPNAM2=ADJUSTL(ADJUSTR(CGPNAM)//HGROUP(IMINUS+1:ILENGP))
	IF(LTYPE)RETURN
	CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,CGPNAM1)
        IF(LPBREAD)THEN
          print *,' VRAISEMBLABLEMENT PB AVEC LE NOM DU GROUPE : ', &
          HGROUP(1:ILENGP)
          RETURN
        ENDIF

	DO J=NAM1,NAM2

        SELECT CASE(NBCNUM)
          CASE(:1)
            IF(J < 10)THEN
      	      WRITE(Y1,'(I1)')J
      	      YNAM1=ADJUSTL(ADJUSTR(CGPNAM)//Y1)
            ELSE IF(J < 100)THEN
      	      WRITE(Y2,'(I2)')J
      	      YNAM1=ADJUSTL(ADJUSTR(CGPNAM)//Y2)
            ELSE IF(J < 1000)THEN
      	      WRITE(Y3,'(I3)')J
      	      YNAM1=ADJUSTL(ADJUSTR(CGPNAM)//Y3)
            ELSE
      	      WRITE(Y4,'(I4)')J
      	      YNAM1=ADJUSTL(ADJUSTR(CGPNAM)//Y4)
            ENDIF
          CASE(2)
            WRITE(Y2,'(I2.2)')J
            YNAM1=ADJUSTL(ADJUSTR(CGPNAM)//Y2)
          CASE(3)
            WRITE(Y3,'(I3.3)')J
            YNAM1=ADJUSTL(ADJUSTR(CGPNAM)//Y3)
          CASE(4)
            WRITE(Y4,'(I4.4)')J
            YNAM1=ADJUSTL(ADJUSTR(CGPNAM)//Y4)
        END SELECT

          YRECFM=ADJUSTL(ADJUSTR(YNAM1)//'.TYPE')
          YNAM1=ADJUSTL(YNAM1)
          ILENG=LEN(CTYPE)
          DEALLOCATE(ITABCHAR)
          ALLOCATE(ITABCHAR(ILENG))
          CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
          IGRID,ILENCH,YCOMMENT,IRESPDIA)
          IF(IRESPDIA == 0)THEN
    	    IF(JM1 == 0)THEN
              JM1=J
    	    ELSE
    	      INCR1=J-JM1
              EXIT
    	    ENDIF
          ENDIF

	ENDDO

	DO J=NAM2,NAM1,-1

        SELECT CASE(NBCNUM)
          CASE(:1)
            IF(J < 10)THEN
      	      WRITE(Y1,'(I1)')J
      	      YNAM2=ADJUSTL(ADJUSTR(CGPNAM)//Y1)
            ELSE IF(J < 100)THEN
      	      WRITE(Y2,'(I2)')J
      	      YNAM2=ADJUSTL(ADJUSTR(CGPNAM)//Y2)
            ELSE IF(J < 1000)THEN
      	      WRITE(Y3,'(I3)')J
      	      YNAM2=ADJUSTL(ADJUSTR(CGPNAM)//Y3)
            ELSE
      	      WRITE(Y4,'(I4)')J
      	      YNAM2=ADJUSTL(ADJUSTR(CGPNAM)//Y4)
            ENDIF
          CASE(2)
            WRITE(Y2,'(I2.2)')J
            YNAM2=ADJUSTL(ADJUSTR(CGPNAM)//Y2)
          CASE(3)
            WRITE(Y3,'(I3.3)')J
            YNAM2=ADJUSTL(ADJUSTR(CGPNAM)//Y3)
          CASE(4)
            WRITE(Y4,'(I4.4)')J
            YNAM2=ADJUSTL(ADJUSTR(CGPNAM)//Y4)
        END SELECT

          YRECFM=ADJUSTL(ADJUSTR(YNAM2)//'.TYPE')
          YNAM2=ADJUSTL(YNAM2)
          ILENG=LEN(CTYPE)
          DEALLOCATE(ITABCHAR)
          ALLOCATE(ITABCHAR(ILENG))
          CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
          IGRID,ILENCH,YCOMMENT,IRESPDIA)
          IF(IRESPDIA == 0)THEN
    	    IF(JM2 == 0)THEN
              JM2=J
    	    ELSE
    	      INCR2=JM2-J
              EXIT
    	    ENDIF
          ENDIF

	ENDDO

        IF(INCR1 /= 0 .AND. INCR1 == INCR2)THEN
          NINCRNAM=INCR1
        ELSE IF(INCR1 /= 0 .AND. INCR1 /= INCR2)THEN
          LPBREAD=.TRUE.
          print *,' Increment Numero Nom Groupe non constant : CAS NON PREVU '
	IF(ALLOCATED(XVAR))THEN
	  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	ENDIF
          RETURN
        ENDIF
    
    	CALL REALLOC_AND_LOAD_RECORDS(HFILEDIA,HLUOUTDIA)
	RETURN

      ELSE

	GPART=.TRUE.

      ENDIF

    ENDIF        !0000000000000000000000000000000000000

  ELSE
! Cas expression groupe sous la forme AA__  (Donc LFIC1=.TRUE. ou .FALSE.)

    GPART=.TRUE.
  ENDIF

  IF(GPART)THEN
! On essaie de rajouter 1, puis 2 puis 3 chiffres
    JM1=0; JM2=0; INCR1=0; INCR2=0
    DO J=1,9999
      IF(J <10)THEN
        WRITE(Y1,'(I1)')J
        YNAM1=ADJUSTL(ADJUSTR(HGROUP)//Y1)
      ELSE IF(J <=99)THEN
        WRITE(Y2,'(I2)')J
        YNAM1=ADJUSTL(ADJUSTR(HGROUP)//Y2)
      ELSE IF(J <= 999)THEN
        WRITE(Y3,'(I3)')J
        YNAM1=ADJUSTL(ADJUSTR(HGROUP)//Y3)
      ELSE
        WRITE(Y4,'(I4)')J
        YNAM1=ADJUSTL(ADJUSTR(HGROUP)//Y4)
      ENDIF
      YRECFM=ADJUSTL(ADJUSTR(YNAM1)//'.TYPE')
      YNAM1=ADJUSTL(YNAM1)
      ILENG=LEN(CTYPE)
      DEALLOCATE(ITABCHAR)
      ALLOCATE(ITABCHAR(ILENG))
      CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
      IGRID,ILENCH,YCOMMENT,IRESPDIA)
      IF(IRESPDIA == 0)THEN
	IF(JM1 == 0)THEN
          JM1=J
	  YNAM1M=YNAM1
	ELSE
	  INCR1=J-JM1
	  YNAM1=YNAM1M
          EXIT
	ENDIF
      ENDIF
    ENDDO
    IF(JM1 /= 0)THEN    !+++++++++++++++++++++++++++++++++++++
    DO J=9999,1,-1
      IF(J <10)THEN
        WRITE(Y1,'(I1)')J
        YNAM2=ADJUSTL(ADJUSTR(HGROUP)//Y1)
      ELSE IF(J <=99)THEN
        WRITE(Y2,'(I2)')J
        YNAM2=ADJUSTL(ADJUSTR(HGROUP)//Y2)
      ELSE IF(J <= 999)THEN
	WRITE(Y3,'(I3)')J
        YNAM2=ADJUSTL(ADJUSTR(HGROUP)//Y3)
      ELSE
        WRITE(Y4,'(I4)')J
        YNAM2=ADJUSTL(ADJUSTR(HGROUP)//Y4)
      ENDIF
      YRECFM=ADJUSTL(ADJUSTR(YNAM2)//'.TYPE')
      YNAM2=ADJUSTL(YNAM2)
      ILENG=LEN(CTYPE)
      DEALLOCATE(ITABCHAR)
      ALLOCATE(ITABCHAR(ILENG))
      CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
      IGRID,ILENCH,YCOMMENT,IRESPDIA)
      IF(IRESPDIA == 0)THEN
	IF(JM2 == 0)THEN
	  JM2=J
	  YNAM2M=YNAM2
	ELSE
	  INCR2=JM2-J
	  YNAM2=YNAM2M
	  EXIT
	ENDIF
      ENDIF
    ENDDO
    ENDIF        !+++++++++++++++++++++++++++++++++++++

    IF(INCR1 /= 0 .AND. INCR1 == INCR2)THEN
      NINCRNAM=INCR1
    ELSE IF(INCR1 /= 0 .AND. INCR1 /= INCR2)THEN
      LPBREAD=.TRUE.
      print *,' Increment Numero Nom Groupe non constant : CAS NON PREVU '
	IF(ALLOCATED(XVAR))THEN
	  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	ENDIF
      RETURN
    ENDIF

    IF(JM1 /= 0 .AND. JM2 /=0)THEN
! On memorise les infos pour realloc_several_records
      CGPNAM=HGROUP(1:LEN_TRIM(HGROUP))
      CGPNAM=ADJUSTL(CGPNAM)
      CGPNAM1=YNAM1
      CGPNAM1=ADJUSTL(CGPNAM1)
      CGPNAM2=YNAM2
      NAM1=JM1; NAM2=JM2
      IF(LTYPE)RETURN
      CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,CGPNAM1)
      CALL REALLOC_AND_LOAD_RECORDS(HFILEDIA,HLUOUTDIA)
      RETURN

    ELSE

! On essaie de rajouter une zone numerique sur 4 positions
      JM1=0; JM2=0; INCR1=0; INCR2=0
      DO J=1,9999
        WRITE(Y4,'(I4.4)')J
        YNAM1=ADJUSTL(ADJUSTR(HGROUP)//Y4)
        YRECFM=ADJUSTL(ADJUSTR(YNAM1)//'.TYPE')
        YNAM1=ADJUSTL(YNAM1)
        ILENG=LEN(CTYPE)
        DEALLOCATE(ITABCHAR)
        ALLOCATE(ITABCHAR(ILENG))
        CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
        IGRID,ILENCH,YCOMMENT,IRESPDIA)
        IF(IRESPDIA == 0)THEN
	  IF(JM1 == 0)THEN
            JM1=J
	    YNAM1M=YNAM1
	  ELSE
	    INCR1=J-JM1
	    YNAM1=YNAM1M
	    EXIT
          ENDIF
        ENDIF
      ENDDO
      IF(JM1 /= 0)THEN    !+++++++++++++++++++++++++++++++++++++
      DO J=9999,1,-1
        WRITE(Y4,'(I4.4)')J
        YNAM2=ADJUSTL(ADJUSTR(HGROUP)//Y4)
        YRECFM=ADJUSTL(ADJUSTR(YNAM2)//'.TYPE')
        YNAM2=ADJUSTL(YNAM2)
        ILENG=LEN(CTYPE)
        DEALLOCATE(ITABCHAR)
        ALLOCATE(ITABCHAR(ILENG))
        CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
        IGRID,ILENCH,YCOMMENT,IRESPDIA)
        IF(IRESPDIA == 0)THEN
	  IF(JM2 == 0)THEN
            JM2=J
	    YNAM2M=YNAM2
	  ELSE
	    INCR2=JM2-J
	    YNAM2=YNAM2M
	    EXIT
	  ENDIF
        ENDIF
      ENDDO
      ENDIF

      IF(INCR1 /= 0 .AND. INCR1 == INCR2)THEN
        NINCRNAM=INCR1
      ELSE IF(INCR1 /= 0 .AND. INCR1 /= INCR2)THEN
        LPBREAD=.TRUE.
        print *,' Increment Numero Nom Groupe non constant : CAS NON PREVU '
	IF(ALLOCATED(XVAR))THEN
	  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	ENDIF
        RETURN
      ENDIF
  
      IF(JM1 /= 0 .AND. JM2 /=0)THEN
! On memorise les infos pour realloc_several_records
        CGPNAM=HGROUP(1:LEN_TRIM(HGROUP))
        CGPNAM=ADJUSTL(CGPNAM)
        CGPNAM1=YNAM1
        CGPNAM1=ADJUSTL(CGPNAM1)
        CGPNAM2=YNAM2
!       print *,' 4 positions CGPNAM,CGPNAM1,CGPNAM2 ',CGPNAM,CGPNAM1,CGPNAM2
        NAM1=JM1; NAM2=JM2
        NBCNUM=4
        IF(LTYPE)RETURN
        CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,CGPNAM1)
        CALL REALLOC_AND_LOAD_RECORDS(HFILEDIA,HLUOUTDIA)
        RETURN

      ELSE

! On essaie de rajouter une zone numerique sur 3 positions
        JM1=0; JM2=0; INCR1=0; INCR2=0
        DO J=1,999
          WRITE(Y3,'(I3.3)')J
          YNAM1=ADJUSTL(ADJUSTR(HGROUP)//Y3)
          YRECFM=ADJUSTL(ADJUSTR(YNAM1)//'.TYPE')
          YNAM1=ADJUSTL(YNAM1)
          ILENG=LEN(CTYPE)
          DEALLOCATE(ITABCHAR)
          ALLOCATE(ITABCHAR(ILENG))
          CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
          IGRID,ILENCH,YCOMMENT,IRESPDIA)
          IF(IRESPDIA == 0)THEN
	    IF(JM1 == 0)THEN
              JM1=J
	      YNAM1M=YNAM1
	    ELSE
	      INCR1=J-JM1
	      YNAM1=YNAM1M
  	      EXIT
            ENDIF
          ENDIF
        ENDDO
        IF(JM1 /= 0)THEN    !+++++++++++++++++++++++++++++++++++++
        DO J=999,1,-1
          WRITE(Y3,'(I3.3)')J
          YNAM2=ADJUSTL(ADJUSTR(HGROUP)//Y3)
          YRECFM=ADJUSTL(ADJUSTR(YNAM2)//'.TYPE')
          YNAM2=ADJUSTL(YNAM2)
          ILENG=LEN(CTYPE)
          DEALLOCATE(ITABCHAR)
          ALLOCATE(ITABCHAR(ILENG))
          CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
          IGRID,ILENCH,YCOMMENT,IRESPDIA)
          IF(IRESPDIA == 0)THEN
	    IF(JM2 == 0)THEN
              JM2=J
	      YNAM2M=YNAM2
	    ELSE
              INCR2=JM2-J
	      YNAM2=YNAM2M
  	      EXIT
            ENDIF
          ENDIF
        ENDDO
        ENDIF

        IF(INCR1 /= 0 .AND. INCR1 == INCR2)THEN
          NINCRNAM=INCR1
        ELSE IF(INCR1 /= 0 .AND. INCR1 /= INCR2)THEN
          LPBREAD=.TRUE.
          print *,' Increment Numero Nom Groupe non constant : CAS NON PREVU '
	IF(ALLOCATED(XVAR))THEN
	  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	ENDIF
          RETURN
        ENDIF
    
        IF(JM1 /= 0 .AND. JM2 /=0)THEN
! On memorise les infos pour realloc_several_records
          CGPNAM=HGROUP(1:LEN_TRIM(HGROUP))
          CGPNAM=ADJUSTL(CGPNAM)
          CGPNAM1=YNAM1
          CGPNAM1=ADJUSTL(CGPNAM1)
          CGPNAM2=YNAM2
          NAM1=JM1; NAM2=JM2
          NBCNUM=3
!         print *,' 3 positions CGPNAM,CGPNAM1,CGPNAM2 ',CGPNAM,CGPNAM1,CGPNAM2
          IF(LTYPE)RETURN
          CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,CGPNAM1)
          CALL REALLOC_AND_LOAD_RECORDS(HFILEDIA,HLUOUTDIA)
          RETURN

        ELSE

! On essaie de rajouter une zone numerique sur 2 positions
          JM1=0; JM2=0; INCR1=0; INCR2=0
          DO J=1,99
            WRITE(Y2,'(I2.2)')J
            YNAM1=ADJUSTL(ADJUSTR(HGROUP)//Y2)
            YRECFM=ADJUSTL(ADJUSTR(YNAM1)//'.TYPE')
            YNAM1=ADJUSTL(YNAM1)
            ILENG=LEN(CTYPE)
            DEALLOCATE(ITABCHAR)
            ALLOCATE(ITABCHAR(ILENG))
            CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
            IGRID,ILENCH,YCOMMENT,IRESPDIA)
            IF(IRESPDIA == 0)THEN
	      IF(JM1 == 0)THEN
                JM1=J
		YNAM1M=YNAM1
	      ELSE
		INCR1=J-JM1
		YNAM1=YNAM1M
    	        EXIT
	      ENDIF
            ENDIF
          ENDDO
          IF(JM1 /= 0)THEN    !+++++++++++++++++++++++++++++++++++++
          DO J=99,1,-1
            WRITE(Y2,'(I2.2)')J
            YNAM2=ADJUSTL(ADJUSTR(HGROUP)//Y2)
            YRECFM=ADJUSTL(ADJUSTR(YNAM2)//'.TYPE')
            YNAM2=ADJUSTL(YNAM2)
            ILENG=LEN(CTYPE)
            DEALLOCATE(ITABCHAR)
            ALLOCATE(ITABCHAR(ILENG))
            CALL FMREAD(HFILEDIA,YRECFM,HLUOUTDIA,ILENG,ITABCHAR, &
            IGRID,ILENCH,YCOMMENT,IRESPDIA)
            IF(IRESPDIA == 0)THEN
              JM2=J
    	  EXIT
            ENDIF
          ENDDO
          ENDIF
      
          IF(INCR1 /= 0 .AND. INCR1 == INCR2)THEN
            NINCRNAM=INCR1
          ELSE IF(INCR1 /= 0 .AND. INCR1 /= INCR2)THEN
            LPBREAD=.TRUE.
            print *,' Increment Numero Nom Groupe non constant : CAS NON PREVU '
	IF(ALLOCATED(XVAR))THEN
	  CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
	ENDIF
            RETURN
          ENDIF
      
          IF(JM1 /= 0 .AND. JM2 /=0)THEN
! On memorise les infos pour realloc_several_records
            CGPNAM=HGROUP(1:LEN_TRIM(HGROUP))
            CGPNAM=ADJUSTL(CGPNAM)
            CGPNAM1=YNAM1
            CGPNAM1=ADJUSTL(CGPNAM1)
            CGPNAM2=YNAM2
            NAM1=JM1; NAM2=JM2
            NBCNUM=2
!           print *,' 2 positions CGPNAM,CGPNAM1,CGPNAM2 ',CGPNAM,CGPNAM1,CGPNAM2
            IF(LTYPE)RETURN
            CALL READ_DIACHRO(HFILEDIA,HLUOUTDIA,CGPNAM1)
            CALL REALLOC_AND_LOAD_RECORDS(HFILEDIA,HLUOUTDIA)
            RETURN

	  ELSE

	  ENDIF

	ENDIF

      ENDIF

    ENDIF

! ELSE

  ENDIF

    LPBREAD=.TRUE.
!************   Le tester dans le pg appelant **************
    IF(INDEX(HGROUP(1:ILENGP),'NPROFILE') /= 0)THEN
    RETURN
    ELSE
    print *,' PB AVEC LE NOM DU GROUPE ou DU PARAMETRE : ',HGROUP(1:ILENGP)
    print *,' VERIFIEZ ET RENTREZ A NOUVEAU VOTRE DIRECTIVE '
    RETURN
    ENDIF


ENDIF

!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE VERIF_GROUP
