!     ######spl
      SUBROUTINE CALUV_FORDIACHRO(KLOOP)
!     ##################################
!
!!****  *CALUV_FORDIACHRO* - Computes a wind,  and moisture
!!                sounding for the emagram mode
!!
!!    PURPOSE
!!    -------
!       For the emagram plots case only, reads U, V, and mix. ratio
!     from the Diachro file, and
!     relocates the results on the mass gridpoint, to obtain a colocated
!     emagram sounding data set.
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      MXM, MYM, MXF, MYF : Shuman averaging operators
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module MODI_SHUMAN      : Contains Shuman operator interfaces
!!
!!         MXM  : mean operator in x direction for a mass variable
!!         MYM  : mean operator in y direction for a mass variable
!!         MXF  : mean operator in x direction for a velocity variable
!!         MYF  : mean operator in y direction for a velocity variable
!!
!!      Module MODD_DIM1       : Contains dimensions
!!
!!         NIMAX,NJMAX,NKMAX :  x, y, and z array dimensions
!!
!!      Module MODD_PARAMETERS : Declares array border depths
!!
!!         JPHEXT   : Horizontal external points number
!!         JPVEXT   : Vertical external points number
!!
!!      Module MODD_LUNIT1     : Declares names and log. unit of files
!!
!!         CLUOUT   : Name of output_listing file
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     06/06/94
!!      Updated  PM  01/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODI_SHUMAN
USE MODI_VERIF_GROUP
USE MODI_REALLOC_AND_LOAD
USE MODD_DIM1
USE MODD_COORD
USE MODD_PARAMETERS
USE MODD_RESOLVCAR
USE MODD_FILES_DIACHRO
USE MODD_ALLOC_FORDIACHRO
USE MODD_SEVERAL_RECORDS
USE MODD_TYPE_AND_LH
USE MODD_PT_FOR_CH_FORDIACHRO

IMPLICIT NONE
!
!             Dummy arguments
!

INTEGER :: KLOOP
!
!             Local variables
!

INTEGER :: IIU, IJU, IKU
INTEGER :: IT, IN, IP
INTEGER :: J, JM, I, IXXX, IXXY
INTEGER :: IRS1, IRSP1, IRS2, IRSP2, IRS3, IRSP3
INTEGER :: JRS1, JRSP1, JRS2, JRSP2, JRS3, JRSP3

CHARACTER(LEN=16) :: YGROUP

REAL :: ZCIINF, ZCISUP, ZCJINF, ZCJSUP         
REAL,DIMENSION(:,:,:,:,:,:),ALLOCATABLE,SAVE :: ZMEANR, ZVAL
REAL,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE :: ZV

!-------------------------------------------------------------------------------
!
!*       1.     COMPUTES SIZES AND RE-ALLOCATES ARRAYS
!               --------------------------------------
IIU=NIMAX+2*JPHEXT
IJU=NJMAX+2*JPHEXT
IKU=NKMAX+2*JPVEXT
!
!
!
!-------------------------------------------------------------------------------
!
!*       2.     READS DATA FROM DIACHRO FILE
!	        ----------------------------
!
! 
NUMFILECUR=NFILESCUR(KLOOP)
DO J=1,NBFILES
  IF(NUMFILES(J) == NUMFILECUR)THEN
    JM=J
  ENDIF
ENDDO
DO J = 1,3
YGROUP(1:LEN(YGROUP))=' '
  IF(NMT == 1)THEN
    IF(J == 1)YGROUP = 'UM'
    IF(J == 2)YGROUP = 'VM'
    IF(J == 3)YGROUP = 'RVM'
  ELSE
    IF(J == 1)YGROUP = 'UT'
    IF(J == 2)YGROUP = 'VT'
    IF(J == 3)YGROUP = 'RVT'
  ENDIF
  YGROUP=ADJUSTL(YGROUP)

  CALL VERIF_GROUP(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
  IF(LPBREAD)THEN
    print *,YGROUP(1:LEN_TRIM(YGROUP)),' N''EXISTE PAS'
    EXIT
  ENDIF
  IF(LGROUP)THEN
    CALL READ_DIACHRO(CFILEDIAS(JM),CLUOUTDIAS(JM),YGROUP)
  ENDIF
  IF(.NOT.LFIC1)THEN
    CALL REALLOC_AND_LOAD(YGROUP)
    IF(LPBREAD)THEN
!     LPBREAD=.FALSE.
      print *,YGROUP(1:LEN_TRIM(YGROUP)),' N''EXISTE PAS DANS', &
      ' L''UN DES FICHIERS '
      EXIT
    ENDIF
  ENDIF

  IF(J == 1)THEN
    ALLOCATE(XU(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
		SIZE(XVAR,5),SIZE(XVAR,6)))
    XU(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
  ELSE IF(J == 2)THEN
    ALLOCATE(XV(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
		SIZE(XVAR,5),SIZE(XVAR,6)))
    XV(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)
    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
  ELSE
    ALLOCATE(XRVJD(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
		 SIZE(XVAR,5),SIZE(XVAR,6)))
    XRVJD(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)
! VOLONTAIREMENT Je ne desalloue pas parce besoin de XDATIME et desallocation
! dans le pg pal comme pour les autres cas.
!    CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)
  ENDIF
ENDDO
!
!
!------------------------------------------------------------------------------
!
!*        3.   RELOCATES THE EMAGRAM POINTS (when profile is defined with
!              ----------------------------  NIRS et NJRS)
!
!
IF(XIRS /= -999.)THEN

  IXXX=SIZE(XXX,1)
  IXXY=SIZE(XXY,1)

  DO I=1,IXXX-1
    IF(XIRSCC >= XXX(I,1) .AND. XIRSCC < XXX(I+1,1))THEN
      IRS1=I
      IRSP1=MIN(I+1,NIH)
      if(nverbia > 0)then
      print *,' XIRSCC,XXX(I,1),XXX(IRSP1,1) ',XIRSCC,XXX(I,1),XXX(IRSP1,1)
      endif
      EXIT
    ENDIF
  ENDDO

  DO J=1,IXXY-1
    IF(XJRSCC >= XXY(J,1) .AND. XJRSCC < XXY(J+1,1))THEN
      JRS1=J
      JRSP1=MIN(J+1,NJH)
      if(nverbia > 0)then
      print *,' XJRSCC,XXY(J,1),XXY(JRSP1,1) ',XJRSCC,XXY(J,1),XXY(JRSP1,1)
      endif
      EXIT
    ENDIF
  ENDDO

  DO I=1,IXXX-1
    IF(XIRSCC >= XXX(I,2) .AND. XIRSCC < XXX(I+1,2))THEN
      IRS2=I
      IRSP2=MIN(I+1,NIH)
      EXIT
    ENDIF
  ENDDO

  DO J=1,IXXY-1
    IF(XJRSCC >= XXY(J,2) .AND. XJRSCC < XXY(J+1,2))THEN
      JRS2=J
      JRSP2=MIN(J+1,NJH)
      EXIT
    ENDIF
  ENDDO

  DO I=1,IXXX-1
    IF(XIRSCC >= XXX(I,3) .AND. XIRSCC < XXX(I+1,3))THEN
      IRS3=I
      IRSP3=MIN(I+1,NIH)
      EXIT
    ENDIF
  ENDDO

  DO J=1,IXXY-1
    IF(XJRSCC >= XXY(J,3) .AND. XJRSCC < XXY(J+1,3))THEN
      JRS3=J
      JRSP3=MIN(J+1,NJH)
      EXIT
    ENDIF
  ENDDO

! Je mets toutes les informations du RS arbitrairement au point NIRS=2,NJRS=2
! qd le profil est defini avec XIRS et XJRS. Cela m'evite d'avoir a modifier
! la partie dans oper (ou je sauvegarde et restitue ap. le RS NIRS et NJRS)
  NIRS=2; NJRS=2

! Grille 1
  IF(IRS1 == IRSP1)THEN
    ZCIINF=0.
    ZCISUP=0.
  ELSE
    ZCIINF=(XXX(IRSP1,1)-XIRSCC)/MAX(1.E-10,(XXX(IRSP1,1)-XXX(IRS1,1)))
    ZCISUP=(XIRSCC-XXX(IRS1,1))/MAX(1.E-10,(XXX(IRSP1,1)-XXX(IRS1,1)))
  ENDIF
  IF(JRS1 == JRSP1)THEN
    ZCJINF=0.
    ZCJSUP=0.
  ELSE
    ZCJINF=(XXY(JRSP1,1)-XJRSCC)/MAX(1.E-10,(XXY(JRSP1,1)-XXY(JRS1,1)))
    ZCJSUP=(XJRSCC-XXY(JRS1,1))/MAX(1.E-10,(XXY(JRSP1,1)-XXY(JRS1,1)))
  ENDIF
  IF(NVERBIA == 10)THEN
    print *,' ZCIINF...',ZCIINF,ZCISUP,ZCJINF,ZCJSUP
    print *,' IRS1,JRS1,IRSP1,JRSP1 ',IRS1,JRS1,IRSP1,JRSP1
    print *,' TH 1 2 3 4 ',XTH(IRS1,JRS1,:,:,:,:)
    print *,' TH 1 2 3 4 ',XTH(IRSP1,JRS1,:,:,:,:)
    print *,' TH 1 2 3 4 ',XTH(IRS1,JRSP1,:,:,:,:)
    print *,' TH 1 2 3 4 ',XTH(IRSP1,JRSP1,:,:,:,:)
    print *,' PRES 1 2 3 4 ',XPRES(IRS1,JRS1,:,:,:,:)
    print *,' PRES 1 2 3 4 ',XPRES(IRSP1,JRS1,:,:,:,:)
    print *,' PRES 1 2 3 4 ',XPRES(IRS1,JRSP1,:,:,:,:)
    print *,' PRES 1 2 3 4 ',XPRES(IRSP1,JRSP1,:,:,:,:)
    print *,' RVJD 1 2 3 4 ',XRVJD(IRS1,JRS1,:,:,:,:)
    print *,' RVJD 1 2 3 4 ',XRVJD(IRSP1,JRS1,:,:,:,:)
    print *,' RVJD 1 2 3 4 ',XRVJD(IRS1,JRSP1,:,:,:,:)
    print *,' RVJD 1 2 3 4 ',XRVJD(IRSP1,JRSP1,:,:,:,:)
  ENDIF
  IF(NVERBIA == 10)THEN
    print *,' U 1 2 3 4 ',XU(IRS2,JRS2,:,:,:,:), &
    XU(IRSP2,JRS2,:,:,:,:),XU(IRS2,JRSP2,:,:,:,:),&
    XU(IRSP2,JRSP2,:,:,:,:)
    print *,' V 1 2 3 4 ',XV(IRS3,JRS3,:,:,:,:), &
    XV(IRSP3,JRS3,:,:,:,:),XV(IRS3,JRSP3,:,:,:,:),&
    XV(IRSP3,JRSP3,:,:,:,:)
  ENDIF

ALLOCATE(ZVAL(SIZE(XTH,1),SIZE(XTH,2),SIZE(XTH,3),SIZE(XTH,4),SIZE(XTH,5),SIZE(XTH,6)))
ALLOCATE(ZV(SIZE(XTH,3),SIZE(XTH,4),SIZE(XTH,5),SIZE(XTH,6)))
! XTH
! ZVAL(IRS1,JRS1,:,:,:,:)=ZCIINF*ZCJINF*XTH(IRS1,JRS1,:,:,:,:)+ &
  DO IP=1,SIZE(XTH,6)
  DO IN=1,SIZE(XTH,5)
  DO IT=1,SIZE(XTH,4)
  ZV(:,IT,IN,IP)=ZCIINF*ZCJINF*XTH(IRS1,JRS1,:,IT,IN,IP)+ &
       ZCIINF*ZCJSUP*XTH(IRS1,JRSP1,:,IT,IN,IP)+ &
       ZCISUP*ZCJINF*XTH(IRSP1,JRS1,:,IT,IN,IP)+ &
       ZCISUP*ZCJSUP*XTH(IRSP1,JRSP1,:,IT,IN,IP)
! ZV(:,IT,IN,IP)=ZVAL(IRS1,JRS1,:,IT,IN,IP)
  XTH(NIRS,NJRS,:,IT,IN,IP)=ZV(:,IT,IN,IP)
  print *,' XTH(NIRS,NJRS,:,IT,IN,IP) ',XTH(NIRS,NJRS,:,IT,IN,IP)
  ENDDO
  ENDDO
  ENDDO
! XPRES
  ZVAL(IRS1,JRS1,:,:,:,:)=ZCIINF*ZCJINF*XPRES(IRS1,JRS1,:,:,:,:)+ &
       ZCIINF*ZCJSUP*XPRES(IRS1,JRSP1,:,:,:,:)+ &
       ZCISUP*ZCJINF*XPRES(IRSP1,JRS1,:,:,:,:)+ &
       ZCISUP*ZCJSUP*XPRES(IRSP1,JRSP1,:,:,:,:)
  ZV(:,:,:,:)=ZVAL(IRS1,JRS1,:,:,:,:)
  XPRES(NIRS,NJRS,:,:,:,:)=ZV(:,:,:,:)
! XRVJD
  ZVAL(IRS1,JRS1,:,:,:,:)=ZCIINF*ZCJINF*XRVJD(IRS1,JRS1,:,:,:,:)+ &
       ZCIINF*ZCJSUP*XRVJD(IRS1,JRSP1,:,:,:,:)+ &
       ZCISUP*ZCJINF*XRVJD(IRSP1,JRS1,:,:,:,:)+ &
       ZCISUP*ZCJSUP*XRVJD(IRSP1,JRSP1,:,:,:,:)
  ZV(:,:,:,:)=ZVAL(IRS1,JRS1,:,:,:,:)
  XRVJD(NIRS,NJRS,:,:,:,:)=ZV(:,:,:,:)
! Grille 2
  IF(IRS2 == IRSP2)THEN
    ZCIINF=0.
    ZCISUP=0.
  ELSE
    ZCIINF=(XXX(IRSP2,2)-XIRSCC)/MAX(1.E-10,(XXX(IRSP2,2)-XXX(IRS2,2)))
    ZCISUP=(XIRSCC-XXX(IRS2,2))/MAX(1.E-10,(XXX(IRSP2,2)-XXX(IRS2,2)))
  ENDIF
  IF(JRS2 == JRSP2)THEN
    ZCJINF=0.
    ZCJSUP=0.
  ELSE
    ZCJINF=(XXY(JRSP2,2)-XJRSCC)/MAX(1.E-10,(XXY(JRSP2,2)-XXY(JRS2,2)))
    ZCJSUP=(XJRSCC-XXY(JRS2,2))/MAX(1.E-10,(XXY(JRSP2,2)-XXY(JRS2,2)))
  ENDIF
! XU
  ZVAL(IRS2,JRS2,:,:,:,:)=ZCIINF*ZCJINF*XU(IRS2,JRS2,:,:,:,:)+ &
       ZCIINF*ZCJSUP*XU(IRS2,JRSP2,:,:,:,:)+ &
       ZCISUP*ZCJINF*XU(IRSP2,JRS2,:,:,:,:)+ &
       ZCISUP*ZCJSUP*XU(IRSP2,JRSP2,:,:,:,:)
  ZV(:,:,:,:)=ZVAL(IRS2,JRS2,:,:,:,:)
  XU(NIRS,NJRS,:,:,:,:)=ZV(:,:,:,:)
! Grille 3
  IF(IRS3 == IRSP3)THEN
    ZCIINF=0.
    ZCISUP=0.
  ELSE
    ZCIINF=(XXX(IRSP3,3)-XIRSCC)/MAX(1.E-10,(XXX(IRSP3,3)-XXX(IRS3,3)))
    ZCISUP=(XIRSCC-XXX(IRS3,3))/MAX(1.E-10,(XXX(IRSP3,3)-XXX(IRS3,3)))
  ENDIF
  IF(JRS3 == JRSP3)THEN
    ZCJINF=0.
    ZCJSUP=0.
  ELSE
    ZCJINF=(XXY(JRSP3,3)-XJRSCC)/MAX(1.E-10,(XXY(JRSP3,3)-XXY(JRS3,3)))
    ZCJSUP=(XJRSCC-XXY(JRS3,3))/MAX(1.E-10,(XXY(JRSP3,3)-XXY(JRS3,3)))
  ENDIF

! XV
  ZVAL(IRS3,JRS3,:,:,:,:)=ZCIINF*ZCJINF*XV(IRS3,JRS3,:,:,:,:)+ &
       ZCIINF*ZCJSUP*XV(IRS3,JRSP3,:,:,:,:)+ &
       ZCISUP*ZCJINF*XV(IRSP3,JRS3,:,:,:,:)+ &
       ZCISUP*ZCJSUP*XV(IRSP3,JRSP3,:,:,:,:)
  ZV(:,:,:,:)=ZVAL(IRS3,JRS3,:,:,:,:)
  XV(NIRS,NJRS,:,:,:,:)=ZV(:,:,:,:)

  DEALLOCATE(ZVAL,ZV)

  IF(NVERBIA == 10)THEN
    print *,' TH,PRES,RVJD,U,V interpoles ',XTH(NIRS,NJRS,:,:,:,:),' ', &
    XPRES(NIRS,NJRS,:,:,:,:),' ',XRVJD(NIRS,NJRS,:,:,:,:),' ',&
    XU(NIRS,NJRS,:,:,:,:),' ',XV(NIRS,NJRS,:,:,:,:)
  ENDIF

ELSE

IF(.NOT.ALLOCATED(ZMEANR))THEN
  ALLOCATE(ZMEANR(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),SIZE(XVAR,4), &
		  SIZE(XVAR,5),SIZE(XVAR,6)))
END IF
! A CORRIGER (Fait le 6/1/97)
ZMEANR(1:IIU-1,:,:,:,:,:)=.5*(XU(1:IIU-1,:,:,:,:,:)+XU(2:IIU,:,:,:,:,:))
!ZMEANR(:,:,:,:,:,:)=MXF(XU)
ZMEANR(IIU,:,:,:,:,:)=2.*ZMEANR(IIU-1,:,:,:,:,:)-ZMEANR(IIU-2,:,:,:,:,:)
XU(:,:,:,:,:,:)=ZMEANR(:,:,:,:,:,:)
!
!ZMEANR(:,:,:,:,:,:)=MYF(XV)
! A CORRIGER (Fait le 6/1/97)
ZMEANR(:,1:IJU-1,:,:,:,:)=.5*(XV(:,1:IJU-1,:,:,:,:)+XV(:,2:IJU,:,:,:,:))
ZMEANR(:,IJU,:,:,:,:)=2.*ZMEANR(:,IJU-1,:,:,:,:)-ZMEANR(:,IJU-2,:,:,:,:)
XV(:,:,:,:,:,:)=ZMEANR(:,:,:,:,:,:)
!
!
!-----------------------------------------------------------------------------
!
!*      4.    EXIT
!             ----
!
DEALLOCATE(ZMEANR)

ENDIF
!
RETURN
END SUBROUTINE  CALUV_FORDIACHRO
