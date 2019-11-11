!     ######spl
      MODULE MODI_OPER_PROCESS
!     #########################
!
INTERFACE
!
SUBROUTINE OPER_PROCESS(KLOOP,HTYPE)
CHARACTER(LEN=*) :: HTYPE
INTEGER          :: KLOOP
END SUBROUTINE OPER_PROCESS
!
END INTERFACE
!
END MODULE MODI_OPER_PROCESS
!     ######spl
      SUBROUTINE OPER_PROCESS(KLOOP,HTYPE)
!     ####################################
!

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
USE MODD_FILES_DIACHRO
USE MODD_ALLOC_FORDIACHRO
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODI_TRACEH_FORDIACHRO
USE MODD_TYPE_AND_LH
USE MODD_DIM1
USE MODD_GRID1
USE MODD_GRID, ONLY:XLONORI,XLATORI
USE MODD_NMGRID
USE MODD_CVERT
USE MODD_MASK3D
USE MODD_TITLE
USE MODD_PARAMETERS
USE MODD_EXPERIM
USE MODN_NCAR    
USE MODN_PARA    
USE MODI_PRECOU_FORDIACHRO
USE MODI_TRACEV_FORDIACHRO
USE MODI_VARFCT
USE MODI_PVFCT
USE MODI_CLOSF
USE MODI_LOADUNITIT
USE MODI_TRAMASK
USE MODI_CONV2XY
USE MODI_TRAPRO_FORDIACHRO
USE MODD_COORD
USE MODD_CONF
USE MODD_SUPER
USE MODD_CST
USE MODD_PVT
USE MODD_DEFCV
USE MODD_MEMCV
USE MODE_GRIDPROJ

IMPLICIT NONE

INTERFACE
  SUBROUTINE COLVECT(KKU,PTEM2D)
  REAL, DIMENSION(:,:),  INTENT(IN) :: PTEM2D
  INTEGER   :: KKU
  END SUBROUTINE COLVECT
END INTERFACE
INTERFACE
	      SUBROUTINE IMCOU_FORDIACHRO(PTABV,PINT,HLEGEND,HTEXT)
	      REAL,DIMENSION(:,:) :: PTABV
	      REAL                :: PINT
	      CHARACTER(LEN=*)    :: HTEXT, HLEGEND
	      END SUBROUTINE IMCOU_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE INTERP_FORDIACHRO(KLREF,KD,KF,PTAB,PTABREF)
      REAL,DIMENSION(:,:,:), INTENT(IN)         :: PTAB 
      REAL,DIMENSION(SIZE(PTAB,1),SIZE(PTAB,2)) :: PTABREF
      INTEGER :: KLREF
      END SUBROUTINE INTERP_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE IMAGE_FORDIACHRO(PTAB,KLREF,PTABINT,KNHI,KNDOT,HTEXTE)
      CHARACTER(LEN=*)   :: HTEXTE
      REAL                :: PTABINT
      REAL,DIMENSION(:,:) :: PTAB
      INTEGER :: KNHI, KNDOT, KLREF
      END SUBROUTINE IMAGE_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE TSOUND_FORDIACHRO(PPRES,PPTEMP,PPQV,PPU,PPV,KNN,HEADER,HTEXTE, &
                        OMXRAT,  &
			OMIXRAT,ODOFRAME,OSAMPLEUV)
      REAL,DIMENSION(:)   ::  PPRES, PPTEMP, PPQV, PPU, PPV
      CHARACTER(LEN=*)               :: HEADER
      CHARACTER(LEN=*)               :: HTEXTE
      LOGICAL                        :: OMXRAT, OMIXRAT, ODOFRAME
      LOGICAL                        :: OSAMPLEUV
      END SUBROUTINE TSOUND_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE TRAXY(PTEMX,PTEMY,KLOOP,HTITX,HTITY,PTIMED,PTIMEF)
      INTEGER    :: KLOOP
      REAL,DIMENSION(:)  :: PTEMX, PTEMY
      REAL               :: PTIMED, PTIMEF
      CHARACTER(LEN=*) :: HTITX, HTITY
      END SUBROUTINE TRAXY
END INTERFACE
INTERFACE
      SUBROUTINE ROTA(PTEM1,PTEMV)
      REAL, DIMENSION(:,:),  INTENT(INOUT) :: PTEM1
      REAL, DIMENSION(:,:),  INTENT(INOUT) :: PTEMV
      END SUBROUTINE ROTA
END INTERFACE
INTERFACE
     SUBROUTINE CALUV_FORDIACHRO(KLOOP)
     INTEGER    :: KLOOP
     END SUBROUTINE CALUV_FORDIACHRO
END INTERFACE
COMMON/TEMV/XZWORKZ,XZZDS,NINX,NINY
COMMON/LOGI/LVERT,LHOR,LPT,LXABS
#include "big.h"
REAL,DIMENSION(N2DVERTX,2500) :: XZWORKZ
!REAL,DIMENSION(1000,400) :: XZWORKZ
!REAL,DIMENSION(200,200) :: XZWORKZ
REAL,DIMENSION(N2DVERTX)     :: XZZDS
!REAL,DIMENSION(1000)     :: XZZDS
!REAL,DIMENSION(200)     :: XZZDS
INTEGER                 :: NINX, NINY
LOGICAL                 :: LVERT, LHOR, LPT, LXABS
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*)  :: HTYPE
INTEGER           :: KLOOP

!
!*       0.1   Local variables
!              ---------------
!
INTEGER   ::   J, JJ
INTEGER   ::   II, IJ, IK, IKU, IKB, IKE, IIU, IJU
INTEGER   ::   JU, ILT
INTEGER   ::   JLOOPP, JLOOPN, JLOOPT, JLOOPK, JLOOPZ, JLOOPK1, JLOOPPF
INTEGER   ::   IZ, IN, ILOOPP
INTEGER   ::   JKLOOP
INTEGER   ::   ILENW, IJLT, ILENT, ILENU, ITIMEND
INTEGER   ::   ISUP, IJSUP, IINF, IJINF
INTEGER   ::   IIB, IIE, IJB, IJE
INTEGER   ::   INBK, INUMK, INUMK1
INTEGER   ::   INDN
INTEGER,SAVE   ::   ISEGM=0, ISEGD=0, ISEGMCOL, ICOLSEGM
INTEGER   ::   IJDEBCOU, IIDEBCOU
INTEGER   ::   IER, INB, IWK, IX, IY, ICOLI
INTEGER   ::   IDEFCV
INTEGER   ::   IINFCV, IISUPCV, IJINFCV, IJSUPCV
INTEGER,SAVE   ::   IIRS, IJRS
INTEGER   ::   IGRID

REAL      ::   ZLAT, ZLON
REAL      ::   ZX, ZY
REAL      ::   ZWL, ZWR, ZWB, ZWT
REAL      ::   ZTIMED, ZTIMEF
REAL      ::   ZZZXD, ZZZXF, ZZZYD, ZZZYF
REAL      ::   ZLW


REAL,DIMENSION(:,:,:),ALLOCATABLE,SAVE :: ZWORK3D, ZPROVI, ZWORK3V
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: ZTEMCV,ZTEM2D, ZWORKRS,ZPROVI2
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: ZTEM1, ZTEMV
REAL,DIMENSION(:),ALLOCATABLE,SAVE     :: ZWORK1D, ZWORKT, ZTEM1D, ZWORKZ, ZWORKZ2
REAL,DIMENSION(:),ALLOCATABLE,SAVE     :: ZTE, ZWO, ZWORKY
REAL,DIMENSION(:,:),ALLOCATABLE,SAVE   :: ZTE2, ZSTAB

CHARACTER(LEN=40)  :: YTEXTE
CHARACTER(LEN=LEN(CTITGAL))  :: YTITGAL
CHARACTER(LEN=2)   :: YC2
CHARACTER(LEN=1)   :: YC1
CHARACTER(LEN=16)  :: YTITX, YTITY, YTEM
CHARACTER(LEN=16)  :: YBID
INTEGER            :: IBID,IRESP

LOGICAL            :: GINVZ
LOGICAL            :: GMXRAT
LOGICAL            :: GII1, GIJ1, GCH
!------------------------------------------------------------------------------
!
YTEXTE(1:LEN(YTEXTE)) = ' '
YTEXTE=ADJUSTL(CGROUP)
CLEGEND(1:LEN(CLEGEND))=' '
!CLEGEND2(1:LEN(CLEGEND2))=' '
!CLEGEND2(1:7)='TIME = '
CTITGAL(1:LEN(CTITGAL))=' '
CUNITGAL(1:LEN(CUNITGAL))=' '
CTIMEC(1:LEN(CTIMEC))=' '
CTIMECS(1:LEN(CTIMECS))=' '
CTIMEC(1:7)='TIME = '
CTIMECS(1:7)='TIME = '
NLOOPT=0
LXABS=LXABSC
if(nverbia > 0)then
  print *,' **oper entree LPRESY,XHMIN,XHMAX ',LPRESY,XHMIN,XHMAX
endif

SELECT CASE(HTYPE)

!*****************************************************************************
!*****************************************************************************
    CASE('CART')

        IF(ALLOCATED(XVAR))THEN
	II=SIZE(XVAR,1)
	IJ=SIZE(XVAR,2)
	IK=SIZE(XVAR,3)  
  
        ELSE
          IF(LRS .OR. LRS1)THEN
            IF(ALLOCATED(XTH))THEN
              II=SIZE(XTH,1)
              IJ=SIZE(XTH,2)
              IK=SIZE(XTH,3)
            ENDIF
          ENDIF
        ENDIF
	if(nverbia > 0)then
	  print *,' **oper Entree II,IJ,IK,KLOOP ',II,IJ,IK,KLOOP
	endif

	IIB=1+JPHEXT; IIE=NIMAX+JPHEXT
	IJB=1+JPHEXT; IJE=NJMAX+JPHEXT
	IIU=NIMAX+2*JPHEXT
	IJU=NJMAX+2*JPHEXT
	IKU=NKMAX+2*JPVEXT
        IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE == 'SU')THEN
          IKU=1
        ENDIF
	IKB=1+JPVEXT; IKE=IKU-JPVEXT

! Traitement des RS
! *****************
	IF(LRS .OR. LRS1)THEN
!
! Cas LRS ou LRS1 et KLOOP = 1 --> Allocation de tableaux pour memoriser
! les infos utiles
!
    IF(KLOOP == 1)THEN

      IF(.NOT.LTINCRDIA(KLOOP,1))THEN
        IF(LRS)THEN
          ILENW=NBTIMEDIA(KLOOP,1)
        ELSE
          ILENW=NSUPERDIA
        ENDIF
      ELSE
        ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
        NTIMEDIA(3,KLOOP,1)+1
        if(nverbia >0)then
        print *,' **oper ilenw ',ILENW
        endif
      ENDIF
      ALLOCATE(XTRS(SIZE(XTH,3),ILENW))
      ALLOCATE(XPRS(SIZE(XTH,3),ILENW))
      ALLOCATE(XURS(SIZE(XTH,3),ILENW))
      ALLOCATE(XVRS(SIZE(XTH,3),ILENW))
      ALLOCATE(XRVRS(SIZE(XTH,3),ILENW))
      ALLOCATE(XTIMRS(ILENW))
      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
        IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
        ALLOCATE(XPRDAT(16,ILENW))
      ENDIF

    ENDIF
!
! Lecture de U V et RV; chargement dans les tableaux de
! travail puis desallocation des tableaux inutiles.
!
	IF(XIRS /= -999.)THEN
	  IIRS=NIRS
	  IJRS=NJRS
	ENDIF
        CALL CALUV_FORDIACHRO(KLOOP)
        if(nverbia >0)then
              print *,' **oper NIRS,NJRS ',NIRS,NJRS
        endif


	IF(.NOT.LTINCRDIA(KLOOP,1))THEN

	  DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
	    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	      NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
	      CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
	    ENDIF

            IF(LRS)THEN
	      XTRS(:,JLOOPT)=XTH(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correction suggeree par Joel et Isa en Decembre 98
!             XTRS(:,JLOOPT)=XTRS(:,JLOOPT)*XEXNREF(NIRS,NJRS,:)
	      XPRS(:,JLOOPT)=(XPRES(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)/ &
			     XP00)**(XRD/XCPD)
	      XTRS(:,JLOOPT)=XTRS(:,JLOOPT)*XPRS(:,JLOOPT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      XPRS(:,JLOOPT)=XPRES(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
            ELSE IF(LRS1)THEN
	      XTRS(:,KLOOP)=XTH(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correction suggeree par Joel et Isa en Decembre 98
!             XTRS(:,KLOOP)=XTRS(:,KLOOP)*XEXNREF(NIRS,NJRS,:)
	      XPRS(:,KLOOP)=(XPRES(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)/ &
			     XP00)**(XRD/XCPD)
	      XTRS(:,KLOOP)=XTRS(:,KLOOP)*XPRS(:,KLOOP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      XPRS(:,KLOOP)=XPRES(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
            ENDIF
          ENDDO

	ELSE

 	  II=0
 	  DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
	    II=II+1
        if(nverbia >0)then
            print *,' **oper JLOOPT II ',JLOOPT,II
        endif
	    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	      CALL LOAD_XPRDAT(II,JLOOPT)
	    ENDIF
	    XTRS(:,II)=XTH(NIRS,NJRS,:,JLOOPT,1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Correction suggeree par Joel et Isa en Decembre 98
!           XTRS(:,II)=XTRS(:,II)*XEXNREF(NIRS,NJRS,:)
            XPRS(:,II)=(XPRES(NIRS,NJRS,:,JLOOPT,1,1)/ &
			     XP00)**(XRD/XCPD)
	    XTRS(:,II)=XTRS(:,II)*XPRS(:,II)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	    XPRS(:,II)=XPRES(NIRS,NJRS,:,JLOOPT,1,1)
 	  ENDDO

	ENDIF

        IF(ALLOCATED(XTH))THEN
 	  DEALLOCATE(XTH)
        ENDIF
        IF(ALLOCATED(XPRES))THEN
 	  DEALLOCATE(XPRES)
        ENDIF

	IF(.NOT.LTINCRDIA(KLOOP,1))THEN

	  DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
            IF(LRS)THEN
	      XURS(:,JLOOPT)=XU(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
	      XVRS(:,JLOOPT)=XV(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
	      XRVRS(:,JLOOPT)=XRVJD(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
	      XTIMRS(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
            ELSE IF(LRS1)THEN
	      XURS(:,KLOOP)=XU(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
	      XVRS(:,KLOOP)=XV(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
	      XRVRS(:,KLOOP)=XRVJD(NIRS,NJRS,:,NTIMEDIA(JLOOPT,KLOOP,1),1,1)
	      XTIMRS(KLOOP)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
            ENDIF
	  ENDDO

	ELSE

 	  II=0
 	  DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
 	    II=II+1
	    XTIMRS(II)=XTRAJT(JLOOPT,1)
 	    XURS(:,II)=XU(NIRS,NJRS,:,JLOOPT,1,1)
 	    XVRS(:,II)=XV(NIRS,NJRS,:,JLOOPT,1,1)
 	    XRVRS(:,II)=XRVJD(NIRS,NJRS,:,JLOOPT,1,1)
 	  ENDDO
	ENDIF

	IF(ALLOCATED(XU))THEN
	  DEALLOCATE(XU)
	ENDIF
	IF(ALLOCATED(XV))THEN
	  DEALLOCATE(XV)
	ENDIF
	IF(ALLOCATED(XRVJD))THEN
	  DEALLOCATE(XRVJD)
	ENDIF
     

      GMXRAT=.TRUE.
      IF(XIRS == -999.)THEN
        IF(NIRS>99) THEN
          IF(NJRS>99) THEN
            WRITE(YTEXTE,'(''I='',I4,'' J='',I4)')NIRS,NJRS
          ELSE
            WRITE(YTEXTE,'(''I='',I4,'' J='',I2)')NIRS,NJRS
          ENDIF
        ELSE
          IF(NJRS>99) THEN
            WRITE(YTEXTE,'(''I='',I2,'' J='',I4)')NIRS,NJRS
          ELSE
            WRITE(YTEXTE,'(''I='',I2,'' J='',I2)')NIRS,NJRS
          ENDIF
        ENDIF
      ELSE
      WRITE(YTEXTE,'(''LAT='',F6.2,'' LON='',F6.2)')XIRS,XJRS
      ENDIF
      YTEXTE=ADJUSTL(YTEXTE)
      IF(NMT == 1)THEN
!       WRITE(CLEGEND(104:110),'(''UM-VM'')')
!       YTEXTE(1:5)='UM-VM'
        CLEGEND(104:108)='UM-VM'
      ELSE
!       WRITE(CLEGEND(104:110),'(''UT-VT'')')
!       YTEXTE(1:5)='UT-VT'
        CLEGEND(104:108)='UT-VT'
      ENDIF
	CALL TABCOL_FORDIACHRO
	CALL GSTXFP(-13,2)

      IF(KLOOP == 1 .AND. LRS)THEN

      DO JLOOPT=1,ILENW
        IF(LPRDAT .AND. ILENW > 1)THEN ! Juin 2001 Ajout des dates ds FICVAL 
! Pour distiller les dates une par une
! Si ILENW = 1 on ne fait rien . OK
	  IF(JLOOPT == 1)THEN
!!!dec 2001
           IF(ALLOCATED(XPRDAT))THEN
!!!dec 2001
	  IF(ALLOCATED(ZPROVI2))DEALLOCATE(ZPROVI2)
	    ALLOCATE(ZPROVI2(16,SIZE(XPRDAT,2)))
            ZPROVI2(:,:)=XPRDAT(:,:)
            DEALLOCATE(XPRDAT)
	    ALLOCATE(XPRDAT(16,1))
	    XPRDAT(:,1)=ZPROVI2(:,JLOOPT)
!!!dec 2001
	  ELSE
	    XPRDAT(:,1)=ZPROVI2(:,JLOOPT)
	  ENDIF
           ELSE
            print *,' *operA XPRDAT NON ALLOUE'
	   ENDIF
!!!dec 2001
        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

	IF(.NOT.LTINCRDIA(KLOOP,1))THEN
	  IF(NVERBIA > 0)THEN
          print *,' KLOOP,LRS,JLOOPT,NTIMEDIA(JLOOPT,KLOOP,1) ', &
          KLOOP,LRS,JLOOPT,NTIMEDIA(JLOOPT,KLOOP,1)
	  ENDIF
	  CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
	ELSE
	  II=NTIMEDIA(1,KLOOP,1)+(JLOOPT-1)*NTIMEDIA(3,KLOOP,1)
        if(nverbia >0)then
          print *,' **oper II de  RESOLV_TIMES(II) ',II
        endif
	  CALL RESOLV_TIMES(II)
	ENDIF
        CTIMEC(1:LEN(CTIMEC))=' '
        CTIMEC(1:3)='  ('
        WRITE(CTIMEC(4:11),'(F8.0)')XTIMRS(JLOOPT)
        CTIMEC(LEN_TRIM(CTIMEC)+1:LEN_TRIM(CTIMEC)+2)='s)'

        GMXRAT=.TRUE.

	DO J=IKB,IKE
	  IF(XRVRS(J,JLOOPT) <=0.)print *,' No dew point line drawn as nil or' &
				,' negative water values were found'
	ENDDO
	CALL GSCLIP(0)
	CALL TSOUND_FORDIACHRO(XPRS(IKB:IKE,JLOOPT),XTRS(IKB:IKE,JLOOPT),  &
		    XRVRS(IKB:IKE,JLOOPT),XURS(IKB:IKE,JLOOPT), &
		    XVRS(IKB:IKE,JLOOPT),IKE-IKB+1,CLEGEND,&
                     YTEXTE,GMXRAT,.TRUE.&
		    ,.FALSE.,.FALSE.)
	CALL GSCLIP(1)
!      CALL NGPICT(1,1)
!      CALL GQACWK(1,IER,INB,IWK)
!      IF(INB > 1)CALL NGPICT(2,3)
        CALL FRAME
      ENDDO
      IF(.NOT.ALLOCATED(XTRS))print *,' XTRS NON ALLOUE'
      IF(.NOT.ALLOCATED(XPRS))print *,' XPRS NON ALLOUE'
      IF(.NOT.ALLOCATED(XURS))print *,' XURS NON ALLOUE'
      IF(.NOT.ALLOCATED(XVRS))print *,' XVRS NON ALLOUE'
      IF(.NOT.ALLOCATED(XRVRS))print *,' XRVRS NON ALLOUE'
      IF(.NOT.ALLOCATED(XTIMRS))print *,' XTIMRS NON ALLOUE'
      if(nverbia > 0)then
        print *,' *operA AV DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS) '
      endif
      DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS)
      if(nverbia > 0)then
         print *,' *operA AP DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS) '
      endif
    ELSE IF(LRS1 .AND. KLOOP == NSUPERDIA)THEN

        GMXRAT=.TRUE.
! On met la date courante du 1er temps demande de la 1ere superposition
        CALL RESOLV_TIMES(NTIMEDIA(1,1,1))
	CALL GSCLIP(0)
	CALL TSOUND_FORDIACHRO(XPRS(IKB:IKE,1),XTRS(IKB:IKE,1),  &
		    XRVRS(IKB:IKE,1),XURS(IKB:IKE,1), &
		    XVRS(IKB:IKE,1),IKE-IKB+1,CLEGEND,YTEXTE,GMXRAT,.TRUE.&
		    ,.FALSE.,.FALSE.)
	CALL GSCLIP(1)
!       CALL NGPICT(1,1)
!       CALL GQACWK(1,IER,INB,IWK)
!       IF(INB > 1)CALL NGPICT(2,3)
        CALL FRAME
            print *,' *operB AV DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS) '
        DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS)
            print *,' *operB AP DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS) '
    ENDIF
    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

	IF(XIRS /= -999.)THEN
	  NIRS=IIRS
	  NJRS=IJRS
	ENDIF

	ELSE
!
! Infos autres que RS
! *******************

	  IF(II == 1 .AND. IJ == 1 .AND. IK == 1)THEN

! Cas compression bilan sur tous les axes ou scalaire unique  f(t)
! ****************************************************************


              IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
                IF(.NOT.LTINCRDIA(KLOOP,1))THEN
		  ILENW=NBTIMEDIA(KLOOP,1)
	        ELSE
		  ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1)+1
	        ENDIF
		ALLOCATE(XPRDAT(16,ILENW))
	      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 

              IF(.NOT.LTINCRDIA(KLOOP,1))THEN
		ALLOCATE(ZWORKT(NBTIMEDIA(KLOOP,1)))
		ALLOCATE(ZWORK1D(NBTIMEDIA(KLOOP,1)))
		DO JLOOPP=1,NBPROCDIA(KLOOP)
		  NLOOPP=NPROCDIA(JLOOPP,KLOOP)

		  CALL LOADUNITIT(JLOOPP,KLOOP)

		  DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		    NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
		    IF(JLOOPT == 1)CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

		    ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
		    ZWORK1D(JLOOPT)=XVAR(1,1,1,NTIMEDIA(JLOOPT,KLOOP,1),1,NPROCDIA(JLOOPP,KLOOP))
		  ENDDO
                  CALL VARFCT(ZWORKT,ZWORK1D,1)
		  IF(KLOOP == NSUPERDIA)CALL FRAME
		ENDDO
		DEALLOCATE(ZWORKT,ZWORK1D)
	      ELSE
		ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1)+1
		ALLOCATE(ZWORKT(ILENW))
		ALLOCATE(ZWORK1D(ILENW))
		DO JLOOPP=1,NBPROCDIA(KLOOP)
		  NLOOPP=NPROCDIA(JLOOPP,KLOOP)

		  CALL LOADUNITIT(JLOOPP,KLOOP)

		  IJLT=0
		  DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		    NLOOPT=JLOOPT
		    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))CALL RESOLV_TIMES(JLOOPT)
		    IJLT=IJLT+1
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      CALL LOAD_XPRDAT(IJLT,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ZWORKT(IJLT)=XTRAJT(JLOOPT,1)
		    ZWORK1D(IJLT)=XVAR(1,1,1,JLOOPT,1,NPROCDIA(JLOOPP,KLOOP))
		  ENDDO
                  CALL VARFCT(ZWORKT,ZWORK1D,1)
		  IF(KLOOP == NSUPERDIA)CALL FRAME
		ENDDO
		DEALLOCATE(ZWORKT,ZWORK1D)
              ENDIF
              IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		DEALLOCATE(XPRDAT)
              ENDIF !  Juin 2001 Ajout des dates ds FICVAL 

              IF(.NOT.LICP .AND. .NOT.LJCP .AND. .NOT.LKCP)THEN

!         Cas scalaire (Impression dim mat. modele et matrice(1,1,1)
!         ------------

	      ELSE IF(LICP .AND. LJCP .AND. LKCP)THEN

!         Cas bilan compresse (Impression dim mat. modele et matrice
!         -------------------  NIL:NIH,NJL:NJH,NKL:NKH)et
!                              matrice (1,1,1)

	      ENDIF

	  ELSE IF(II == 1 .AND. IJ == 1 .AND. IK /= 1)THEN

! Cas compression bilan sur axes X et Y ou PV -->  Profil vertical
! ****************************************************************
!
	      IDEFCV=0                      !%%%%%%%%%%%%%%%%%%%%%%%%%%
	      IF(LDEFCV2CC)THEN
	        LDEFCV2CC=.FALSE.
	        IDEFCV=1
	      ENDIF                         !%%%%%%%%%%%%%%%%%%%%%%%%%%
	      L1DT=.TRUE.
	      ALLOCATE(ZTEM1D(IKU),ZWORKZ(IKU))

              DO JLOOPP=1,NBPROCDIA(KLOOP)
	         NLOOPP=NPROCDIA(JLOOPP,KLOOP)

!!! Octobre 2001
                IF(JLOOPP > 1 .AND. LUMVMPV .AND. LPV)EXIT
!!! Octobre 2001
		IF(LPVKT .AND. NSUPERDIA>1)THEN
		  IF(NBPROCDIA(KLOOP)>1 .OR. NBLVLKDIA(KLOOP,1)>1)THEN
		    print *,' _PVKT_  SUPERPOSITIONS : '
!fuji    print *,'         On ne peut definir de part de d''autre '&
!fuji    &'de _ON_ qu''1 seul processus et 1 seul niveau'
		    print *,'         On ne peut definir de part de d''autre '
		    print *,'de _ON_ qu''1 seul processus et 1 seul niveau'
		    print *,' Nb de niveaux demandes   : ',NBLVLKDIA(KLOOP,1)
		    print *,' Nb de processus demandes : ',NBPROCDIA(KLOOP)
		    print *,' *** MODIFIEZ VOTRE DIRECTIVE *** '
		    EXIT
		  ENDIF
		ENDIF

! Modif AOUT 97
	        ZTEM1D(:)=XSPVAL; ZWORKZ(:)=0.
!               ZTEM1D(:)=0.; ZWORKZ(:)=0.

		  CALL LOADUNITIT(JLOOPP,KLOOP)
!!!!!Mars 2000
                  IF(LUMVM)THEN
		    NMGRID=1
		  ENDIF
                  IF(LUMVMPV)THEN
		    NMGRID=1
		  ENDIF
!!!!!Mars 2000

		  CALL COMPCOORD_FORDIACHRO(NMGRID)
! Expression temps non incrementale
		IF(.NOT.LTINCRDIA(KLOOP,1))THEN

                DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)

	          CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
! Chargement cas PV

	          ZTEM1D(NKL:NKH)=XVAR(1,1,: &
		  ,NTIMEDIA(JLOOPT,KLOOP,1),1,NPROCDIA(JLOOPP,KLOOP))

		  ZWORKZ(:)=XXZ(:,NMGRID)
!                 print * ,'**operoper NMGRID XXZ ',NMGRID
!                 print * ,XXZ(:,NMGRID)
		  IF(NIL /= 1 .OR. NJL /= 1)THEN
		    IF(LICP .OR. LJCP)THEN
!      	              print *,'**operoper LICP LJCP ',LICP,LJCP
		    ELSE
		    ZWORKZ(:)=XZZ(NIL,NJL,:)
		    ENDIF
		    IF(NKL == 1 .AND. NKH == IKU)THEN
		      ZTEM1D(1)=XSPVAL
		      ZTEM1D(IKU)=XSPVAL
		    ENDIF
		  ENDIF
                 

    	          IF(LPV)THEN
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	              IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		      ALLOCATE(XPRDAT(16,1))
		      CALL LOAD_XPRDAT(1,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

		    IF(LUMVMPV)THEN
		      LPV=.FALSE. ; LPVT=.TRUE.
		      IF(JLOOPP == 1)THEN
!!!! Octobre 2001
                        ILENW=1
                        ALLOCATE(ZTEM2D(1:IKU,ILENW))
                        ALLOCATE(ZWORKT(ILENW))
                        ZWORKT=NLOOPT
			IF(ALLOCATED(XTEM2D))THEN
			  DEALLOCATE(XTEM2D)
			ENDIF
		        ALLOCATE(XTEM2D(1:IKU,ILENW))
			XTEM2D=XSPVAL
			IF(ALLOCATED(XTEM2D2))THEN
			  DEALLOCATE(XTEM2D2)
			ENDIF
		        ALLOCATE(XTEM2D2(1:IKU,ILENW))
			XTEM2D2=XSPVAL
                        XTEM2D(:,1)=ZTEM1D
                        XTEM2D2(NKL:NKH,1)=XVAR(1,1,: &
                        ,NTIMEDIA(JLOOPT,KLOOP,1),1,NPROCDIA(JLOOPP+1,KLOOP))
			IF(NBPROCDIA(KLOOP) == 3)THEN
			  ZTEM2D=XSPVAL
                          ZTEM2D(NKL:NKH,1)=XVAR(1,1,: &
                          ,NTIMEDIA(JLOOPT,KLOOP,1),1,NPROCDIA(JLOOPP+2,KLOOP))
                          
                          CALL COLVECT(IKU,ZTEM2D)
                         ENDIF
                         CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                         IF(LUMVMPV)THEN
		           LPV=.TRUE. ; LPVT=.FALSE.
                         ENDIF
		         DEALLOCATE(ZTEM2D,ZWORKT)
                         IF(ALLOCATED(XTEM2D))THEN
		           DEALLOCATE(XTEM2D)
		         ENDIF
		         IF(ALLOCATED(XTEM2D2))THEN
		           DEALLOCATE(XTEM2D2)
		         ENDIF
                         LCOLPVT=.FALSE.
                       ENDIF

                    ELSE
!!!! Octobre 2001

		      CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)

		    ENDIF
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      DEALLOCATE(XPRDAT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    IF(KLOOP == NSUPERDIA)CALL FRAME
    	          ELSE IF(LPVT .OR. LPVKT)THEN
		    IF(JLOOPT == 1)THEN
		      ILENW=NBTIMEDIA(KLOOP,1)
		      ALLOCATE(ZTEM2D(1:IKU,ILENW))
		      ZTEM2D=XSPVAL
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	              IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		      ALLOCATE(XPRDAT(16,ILENW))
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
!!!!!Mars 2000
                      IF(LUMVM)THEN
			IF(ALLOCATED(XTEM2D))THEN
			  DEALLOCATE(XTEM2D)
			ENDIF
		        ALLOCATE(XTEM2D(1:IKU,ILENW))
			XTEM2D=XSPVAL
		      ENDIF

                      IF(LUMVMPV .AND. JLOOPP == 1)THEN
			IF(ALLOCATED(XTEM2D))THEN
			  DEALLOCATE(XTEM2D)
			ENDIF
		        ALLOCATE(XTEM2D(1:IKU,ILENW))
			XTEM2D=XSPVAL
			IF(ALLOCATED(XTEM2D2))THEN
			  DEALLOCATE(XTEM2D2)
			ENDIF
		        ALLOCATE(XTEM2D2(1:IKU,ILENW))
			XTEM2D2=XSPVAL
		      ENDIF
!!!!!Mars 2000
		      ALLOCATE(ZWORKT(ILENW))
		    ENDIF
		    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
                    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
		    ZTEM2D(NKL:NKH,JLOOPT)=  XVAR(1,1,:,  &
		      NTIMEDIA(JLOOPT,KLOOP,1),1,NPROCDIA(JLOOPP,KLOOP))
!!!!!Mars 2000
                    IF(LUMVM)THEN
		      XTEM2D(NKL:NKH,JLOOPT)= XU(1,1,:,  &
		        NTIMEDIA(JLOOPT,KLOOP,1),1,NPROCDIA(JLOOPP,KLOOP))
		    ENDIF
!!!!!Mars 2000
		    IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
		      XPVMIN=MINVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      XPVMAX=MAXVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      CALL VALMNMX(XPVMIN,XPVMAX)
                      IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
			XPVMIN=XPVMIN-1.
			XPVMAX=XPVMAX+1.
                      ENDIF
		      IF(NKL == 1 .AND. NKH == IKU)THEN
	                ZTEM2D(1,:)=XSPVAL
	                ZTEM2D(IKU,:)=XSPVAL
		      ENDIF

		      IF(LUMVMPV)THEN
			IF(JLOOPP == 1)THEN
! Memorisation de U
			  XTEM2D=ZTEM2D
			  CYCLE
			ELSEIF(JLOOPP == 2)THEN
			  IF(JLOOPP == NBPROCDIA(KLOOP))THEN
			    XTEM2D2=ZTEM2D
			  ELSE
			    XTEM2D2=ZTEM2D
			    CYCLE
			  ENDIF
			ELSEIF(JLOOPP == 3)THEN 
                          CALL COLVECT(IKU,ZTEM2D)
			ENDIF
		      ENDIF

		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                      IF(LPRDAT) DEALLOCATE(XPRDAT) ! Juin 2001 Ajout des dates ds FICVAL 
		      DEALLOCATE(ZTEM2D,ZWORKT)
                      IF(ALLOCATED(XTEM2D))THEN
		        DEALLOCATE(XTEM2D)
		      ENDIF
		      IF(ALLOCATED(XTEM2D2))THEN
		        DEALLOCATE(XTEM2D2)
		      ENDIF
                      LCOLPVT=.FALSE.
		      IF(.NOT.LPBREAD)THEN
		        IF(KLOOP == NSUPERDIA)CALL FRAME
		      ENDIF
		    ENDIF
    	          ENDIF
	        ENDDO
		ELSE
! Expression temps incrementale !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		  NLOOPT=JLOOPT
	          CALL RESOLV_TIMES(JLOOPT)
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)

	          ZTEM1D(NKL:NKH)=XVAR(1,1,: &
		  ,JLOOPT,1,NPROCDIA(JLOOPP,KLOOP))

		  ZWORKZ(:)=XXZ(:,NMGRID)
!                 print * ,'**operoper NMGRID XXZ ',NMGRID
!                 print * ,XXZ(:,NMGRID)
		  IF(NIL /= 1 .OR. NJL /= 1)THEN
		    IF(LICP .OR. LJCP)THEN
!                     print * ,'**operoper LICP, LJCP ',LICP, LJCP
		    ELSE
		    ZWORKZ(:)=XZZ(NIL,NJL,:)
		    ENDIF
		    IF(NKL == 1 .AND. NKH == IKU)THEN
		      ZTEM1D(1)=XSPVAL
		      ZTEM1D(IKU)=XSPVAL
		    ENDIF
		  ENDIF

    	          IF(LPV)THEN
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	              IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		      ALLOCATE(XPRDAT(16,1))
		      CALL LOAD_XPRDAT(1,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

		    IF(LUMVMPV)THEN
		      LPV=.FALSE. ; LPVT=.TRUE.
!!!Octobre 2001
		      IF(JLOOPP == 1)THEN
                        ILENW=1
                        ALLOCATE(ZTEM2D(1:IKU,ILENW))
                        ALLOCATE(ZWORKT(ILENW))
                        ZWORKT=NLOOPT
                        IF(ALLOCATED(XTEM2D))THEN
                          DEALLOCATE(XTEM2D)
                        ENDIF
                        ALLOCATE(XTEM2D(1:IKU,ILENW))
                        XTEM2D=XSPVAL
                        IF(ALLOCATED(XTEM2D2))THEN
                          DEALLOCATE(XTEM2D2)
                        ENDIF
                        ALLOCATE(XTEM2D2(1:IKU,ILENW))
                        XTEM2D2=XSPVAL
                        XTEM2D(:,1)=ZTEM1D
                        XTEM2D2(NKL:NKH,1)=XVAR(1,1,: &
                        ,JLOOPT,1,NPROCDIA(JLOOPP+1,KLOOP))
                        IF(NBPROCDIA(KLOOP) == 3)THEN
			  ZTEM2D=XSPVAL
                          ZTEM2D(NKL:NKH,1)=XVAR(1,1,: &
                          ,JLOOPT,1,NPROCDIA(JLOOPP+2,KLOOP))

                          CALL COLVECT(IKU,ZTEM2D)
                        ENDIF
                        CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                         IF(LUMVMPV)THEN
		           LPV=.TRUE. ; LPVT=.FALSE.
                         ENDIF
                        DEALLOCATE(ZTEM2D,ZWORKT)
                        IF(ALLOCATED(XTEM2D))THEN
                          DEALLOCATE(XTEM2D)
                        ENDIF
                        IF(ALLOCATED(XTEM2D2))THEN
                          DEALLOCATE(XTEM2D2)
                        ENDIF
                        LCOLPVT=.FALSE.
		      ENDIF

                    ELSE
!!!Octobre 2001
		      CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
!!!Octobre 2001
		    ENDIF
!!!Octobre 2001
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      DEALLOCATE(XPRDAT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    IF(KLOOP == NSUPERDIA)CALL FRAME

    	          ELSE IF(LPVT .OR. LPVKT)THEN

		    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		      ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1) +1
	              IF(NVERBIA > 0)THEN
                      print *,' OPER  NTIMEDIA(2,KLOOP,1) NTIMEDIA(1,KLOOP,1) NTIMEDIA(3,KLOOP,1) ILENW ', &
                      NTIMEDIA(2,KLOOP,1),NTIMEDIA(1,KLOOP,1),NTIMEDIA(3,KLOOP,1), &
		      ILENW, &
                      XTIMEDIA(2,KLOOP,1),XTIMEDIA(1,KLOOP,1),XTIMEDIA(3,KLOOP,1)
		      ENDIF

		      ITIMEND=NTIMEDIA(1,KLOOP,1) + &
		      (((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/  &
		      NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))

	              IF(NVERBIA > 0)THEN
		      print *,' ITIMEND  A',ITIMEND
		      ENDIF

		      IF(ALLOCATED(ZTEM2D))THEN
			DEALLOCATE(ZTEM2D)
		      ENDIF
		      IF(ALLOCATED(ZWORKT))THEN
			DEALLOCATE(ZWORKT)
		      ENDIF
		      ALLOCATE(ZTEM2D(1:IKU,ILENW))
		      ZTEM2D=XSPVAL
!!!!!Mars 2000
                      IF(LUMVM)THEN
			IF(ALLOCATED(XTEM2D))THEN
			  DEALLOCATE(XTEM2D)
			ENDIF
		        ALLOCATE(XTEM2D(1:IKU,ILENW))
			XTEM2D=XSPVAL
		      ENDIF

                      IF(LUMVMPV .AND. JLOOPP == 1)THEN
			IF(ALLOCATED(XTEM2D))THEN
			  DEALLOCATE(XTEM2D)
			ENDIF
		        ALLOCATE(XTEM2D(1:IKU,ILENW))
			XTEM2D=XSPVAL
			IF(ALLOCATED(XTEM2D2))THEN
			  DEALLOCATE(XTEM2D2)
			ENDIF
		        ALLOCATE(XTEM2D2(1:IKU,ILENW))
			XTEM2D2=XSPVAL
		      ENDIF
!!!!!Mars 2000
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,ILENW))
		      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                      ALLOCATE(ZWORKT(ILENW))
		      IJLT=0
		    ENDIF

		    IJLT=IJLT+1
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			CALL LOAD_XPRDAT(IJLT,NLOOPT)
		      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ZWORKT(IJLT)=XTRAJT(JLOOPT,1)
		      if(nverbia >0)then
!                      print *,' **oper AV ZTEM2D(NKL:NKH,IJLT)= '
		    endif
		    ZTEM2D(NKL:NKH,IJLT)= &
		    XVAR(1,1,:,  &
		    JLOOPT,1,NPROCDIA(JLOOPP,KLOOP))
		      if(nverbia >0)then
!                       print *,' **oper AP ZTEM2D(NKL:NKH,IJLT)= '
		      endif
!!!!!Mars 2000
                      IF(LUMVM)THEN
		        XTEM2D(NKL:NKH,IJLT)= &
		        XU(1,1,:,  &
			JLOOPT,1,NPROCDIA(JLOOPP,KLOOP))
		      ENDIF
!!!!!Mars 2000

!                   IF(JLOOPT == NTIMEDIA(2,KLOOP,1))THEN
		    IF(JLOOPT == ITIMEND)THEN
		      XPVMIN=MINVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      XPVMAX=MAXVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      CALL VALMNMX(XPVMIN,XPVMAX)
		      if(nverbia >0)then
		        print *,' **oper AP CALL VALMNMX(XPVMIN,XPVMAX)'
		      endif
                      IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
			XPVMIN=XPVMIN-1.
			XPVMAX=XPVMAX+1.
                      ENDIF
		      IF(NKL == 1 .AND. NKH == IKU)THEN
	                ZTEM2D(1,:)=XSPVAL
	                ZTEM2D(IKU,:)=XSPVAL
		      ENDIF

		      IF(LUMVMPV)THEN        !llllllllllllllllllll

			IF(JLOOPP == 1)THEN  !kkkkkkkkkkkkkkkkkkkkkkk
! Memorisation de U
			  XTEM2D=ZTEM2D
			  CYCLE
			ELSEIF(JLOOPP == 2)THEN !kkkkkkkkkkkkkkkkkkkkk
			  IF(JLOOPP == NBPROCDIA(KLOOP))THEN
			    XTEM2D2=ZTEM2D
			  ELSE
			    XTEM2D2=ZTEM2D
			    CYCLE
			  ENDIF
			ELSEIF(JLOOPP == 3)THEN !kkkkkkkkkkkkkkkkkkkkk
                          CALL COLVECT(IKU,ZTEM2D)
			ENDIF         !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
		      ENDIF           !llllllllllllllllllllllllllllllllll

		      if(nverbia >0)then
			print *,' ** oper AV CALL PVFCT xx'
		      endif
		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      DEALLOCATE(XPRDAT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		      DEALLOCATE(ZWORKT,ZTEM2D)
		      if(nverbia >0)then
			print *,' ** oper AP CALL PVFCT xx'
		      endif
                      IF(ALLOCATED(XTEM2D))THEN
		        DEALLOCATE(XTEM2D)
		      ENDIF
		      IF(ALLOCATED(XTEM2D2))THEN
		        DEALLOCATE(XTEM2D2)
		      ENDIF
			LCOLPVT=.FALSE.
		      IF(.NOT.LPBREAD)THEN
		        IF(KLOOP == NSUPERDIA)CALL FRAME
		      if(nverbia >0)then
			print *,' ** oper AP CALL FRAME xx'
		      endif
		      ENDIF

		    ENDIF         ! Fin if=ITIMEND
    	          ENDIF
	        ENDDO    ! fin boucle temporelle
		ENDIF    ! Tps increm ou non

	      ENDDO     !    Processus
	      DEALLOCATE(ZTEM1D,ZWORKZ)
            IF(.NOT.LICP .AND. .NOT.LJCP .AND. .NOT.LKCP)THEN
!
!  Cas PV enregistre comme tel 
!
	    ELSE IF(LICP .AND. LJCP .AND. .NOT.LKCP)THEN
! (Impression dim mat. modele et matrice(NIL:NIH,NJL:NJH,
!  NKL:NKH) et matrice(1,1,NKL:NKH)
	    ENDIF

	    IF(IDEFCV==1)THEN                !%%%%%%%%%%%%%%%%%%%%%%%%%%
	      LDEFCV2CC=.TRUE.
	      IDEFCV=0
	    ENDIF                           !%%%%%%%%%%%%%%%%%%%%%%%%%%


	  ELSE IF(II == 1 .AND. IJ /= 1 .AND. IK /= 1 .AND. LICP)THEN

! Cas compression bilan sur axe X -->  Plan vertical // Y
! *******************************************************
! (Impression dim mat. modele et matrice(NIL:NIH,NJL:NJH,
!  NKL:NKH) et matrice(1,NJL:NJH,NKL:NKH)
            LCVYZ=.TRUE.
	    IDEFCV=0                         !%%%%%%%%%%%%%%%%%%%%%%%%%%
	    IF(LDEFCV2CC)THEN
	      LDEFCV2CC=.FALSE.
	      IDEFCV=1
	    ENDIF                            !%%%%%%%%%%%%%%%%%%%%%%%%%%
	    IF(.NOT.L2DBY)THEN
	      IJINF=MAX(IJB,NJL)
	      IJSUP=MIN(IJE,NJH)
	      print *,' 2D Vertical // Y '
	      print *,' Limites J par defaut (L2DBY=.FALSE.)(par / au domaine integral de simulation, points de garde compris) :',&
&             ' MAX(IJB,NJL) - MIN(IJE,NJH) ',IJINF,' - ',IJSUP
              print *,' Si vous voulez selectionner les limites en J, mettez : ',&
&             'L2DBY=.TRUE.' 
              print *,' et definissez : NJDEBCOU=    NLMAX= '
	    ELSE
	      IJINF=NJDEBCOU     
	      IJSUP=NJDEBCOU+NLMAX-1
	      IJSUP=MIN(IJSUP,NJH)
	    ENDIF
	    ALLOCATE(ZTEM2D(1:IJSUP-IJINF+1,1:IKU))
	    NINX=IJSUP-IJINF+1
	    NINY=IKU
            NLMAX=NINX
            NLANGLE=90
            NJDEBCOU=IJINF
	    IIDEBCOU=-999
	    IF(NIDEBCOU /= NIL)THEN
	      IIDEBCOU=NIDEBCOU
              NIDEBCOU=NIL
!	      print *,' NIDEBCOU force a la valeur de NIL ',NIL,' pour ', &
!&            'obtention altitudes correctes '
!             print *,' AP utilisation, sera remis a la valeur precedente : ', &
!             IIDEBCOU
            ENDIF
	    LVERT=.TRUE.
	    LHOR=.FALSE.
	    LPT=LPXT
	    IF(NSUPERDIA > 1)THEN
		    IF(LMINUS .OR. LPLUS)THEN
		      IF(NBPM > 1)THEN
			DO JU=1,NBPM
			  IF(NUMPM(JU) == 3)THEN
		            LSUPER=.TRUE.
			    EXIT
			  ELSE
		            LSUPER=.FALSE.
			  ENDIF
			ENDDO
		      ELSE
		        LSUPER=.FALSE.
		      ENDIF
		    ELSE
		      LSUPER=.TRUE.
		    ENDIF
            ELSE
	      LSUPER=.FALSE.
	    ENDIF
	    IF(KLOOP == 1)NSUPER=0
	    DO JLOOPP=1,NBPROCDIA(KLOOP)      !--- LCVYZ-------------
	       NLOOPP=NPROCDIA(JLOOPP,KLOOP)
               NMGRID=NGRIDIA(NLOOPP)
	      IF(JLOOPP == 1)NSUPER=0

		  CALL LOADUNITIT(JLOOPP,KLOOP)

              ILENT=LEN_TRIM(CTITGAL)
	      ILENU=LEN_TRIM(CUNITGAL)
	      YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
	      YTEXTE(ILENT+1:ILENT+1)=' '
	      YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
              IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		ALLOCATE(XPRDAT(16,1))
              ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

	      IF(.NOT.LTINCRDIA(KLOOP,1))THEN

		DO JLOOPT=1,NBTIMEDIA(KLOOP,1)

		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	            CALL LOAD_XPRDAT(1,NLOOPT)
		  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		  CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
    	          IF(.NOT. LSUPER .OR. (LSUPER .AND. NSUPER == 0))THEN
    	            DO J=1,NINX
    		      XZWORKZ(J,1:IKU)=XXZ(:,NMGRID)
                    ENDDO
    	            XZZDS(1:NINX)=XXY(IJINF:IJSUP,NMGRID)
    	            ZWL=XZZDS(1); ZWR=XZZDS(NINX)
    	            IF((XHMIN == 0. .AND. XHMAX == 0.) .OR. (XHMAX<=XHMIN))THEN
    		      XHMIN=0.
    		      XHMAX=XZWORKZ(1,IKE)
                    ENDIF
    	            ZWB=XHMIN; ZWT=XHMAX
    	            CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,1)
    	            CALL GSCLIP(1)
    	            CALL CPSETI('SET',0)
    	            CALL CPSETI('MAP',4)
    	          ENDIF
		  ZTEM2D=XSPVAL
		  ZTEM2D(1:IJSUP-IJINF+1,NKL:NKH)=XVAR(1, &
		  IJINF-NJL+1:IJSUP-NJL+1,:,NTIMEDIA(JLOOPT,KLOOP,1),&
		  1,NPROCDIA(JLOOPP,KLOOP))
    	          IF(NKL < IKB)THEN
    		    ZTEM2D(:,1:IKB-1)=XSPVAL
                  ENDIF
    	          IF(NKH > IKE)THEN
    		    ZTEM2D(:,IKE+1:IKU)=XSPVAL
                  ENDIF
		  if(nverbia >0)THEN
		    print *,' ** oper appel imcou  Ytexte ',YTEXTE(1:LEN_TRIM(YTEXTE))
		  endif
		  IF(KLOOP == 1)NSUPER=0
                  CALL IMCOU_FORDIACHRO(ZTEM2D,XDIAINT,CLEGEND,YTEXTE(1:LEN_TRIM&
                  (YTEXTE)))
!                 IF(KLOOP == NSUPERDIA)CALL FRAME
		  IF(KLOOP == NSUPERDIA)THEN
		    CALL NGPICT(1,1)
		    CALL GQACWK(1,IER,INB,IWK)
		    IF(INB > 1)CALL NGPICT(2,3)
		  ENDIF
		ENDDO
	      ELSE
		DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),  &
			  NTIMEDIA(3,KLOOP,1)
		  NLOOPT=JLOOPT
                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	            CALL LOAD_XPRDAT(1,NLOOPT)
	          ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		  CALL RESOLV_TIMES(JLOOPT)
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)
    	          IF(.NOT. LSUPER .OR. (LSUPER .AND. NSUPER == 0))THEN
    	            DO J=1,NINX
    		      XZWORKZ(J,1:IKU)=XXZ(:,NMGRID)
                    ENDDO
    	            XZZDS(1:NINX)=XXY(IJINF:IJSUP,NMGRID)
    	            ZWL=XZZDS(1); ZWR=XZZDS(NINX)
    	            IF((XHMIN == 0. .AND. XHMAX == 0.) .OR. (XHMAX<=XHMIN))THEN
    		      XHMIN=0.
    		      XHMAX=XZWORKZ(1,IKE)
                    ENDIF
    	            ZWB=XHMIN; ZWT=XHMAX
    	            CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,1)
    	            CALL GSCLIP(1)
    	            CALL CPSETI('SET',0)
    	            CALL CPSETI('MAP',4)
    	          ENDIF
		  ZTEM2D=XSPVAL
		  ZTEM2D(1:IJSUP-IJINF+1,NKL:NKH)=XVAR(1, &
		  IJINF-NJL+1:IJSUP-NJL+1,:,JLOOPT,1, &
		  NPROCDIA(JLOOPP,KLOOP))
    	          IF(NKL < IKB)THEN
    		    ZTEM2D(:,1:IKB-1)=XSPVAL
                  ENDIF
    	          IF(NKH > IKE)THEN
    		    ZTEM2D(:,IKE+1:IKU)=XSPVAL
                  ENDIF
		  if(nverbia >0)THEN
		    print *,' ** oper appel imcou  Ytexte ',YTEXTE(1:LEN_TRIM(YTEXTE))
		  endif
                  IF(KLOOP ==1)NSUPER=0
		  CALL IMCOU_FORDIACHRO(ZTEM2D,XDIAINT,CLEGEND,YTEXTE(1:LEN_TRIM&
                  (YTEXTE)))
!                 IF(KLOOP == NSUPERDIA)CALL FRAME
		  IF(KLOOP == NSUPERDIA)THEN
		    CALL NGPICT(1,1)
		    CALL GQACWK(1,IER,INB,IWK)
		    IF(INB > 1)CALL NGPICT(2,3)
		  ENDIF
                ENDDO
              ENDIF
            ENDDO                             !--- LCVYZ-------------
            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	      DEALLOCATE(XPRDAT)
	    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	    DEALLOCATE(ZTEM2D)
	    IF(IIDEBCOU /= -999)THEN
	      NIDEBCOU=IIDEBCOU
	    ENDIF

	    IF(IDEFCV==1)THEN                 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      LDEFCV2CC=.TRUE.
	      IDEFCV=0
	    ENDIF                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  ELSE IF((II == 1 .OR. IIE-IIB == 0) .AND. IJ /= 1 .AND. IK == 1)THEN

! Cas compression bilan sur axes X et Z -->  Profil horizontal // Y
! mais a representer comme f(t)
! ********************************************************************
! (Impression dim mat. modele et matrice(NIL:NIH,NJL:NJH,
!  NKL:NKH) et matrice(1,NJL:NJH,1)
            print *,' Profil horizontal // Y'
            IINF=NIINF;ISUP=NISUP;IJINF=NJINF;IJSUP=NJSUP
	    if(nverbia > 0)then
	    print *,'IINF,ISUP,IJINF,IJSUP ',IINF,ISUP,IJINF,IJSUP
	    endif
	    IF(II == 1)THEN
	      GII1=.TRUE.
	    ELSE
	      GII1=.FALSE.
              LCH=.FALSE.
	    ENDIF

	    IF(GII1)THEN
              IF(.NOT.L2DBY)THEN
	        NIINF=1; NISUP=1
	        NJINF=MAX(IJB,NJL); NJSUP=MIN(IJE,NJH)
	        print *,' Profil horizontal // Y '
	        print *,' Limites J par defaut (L2DBY=.FALSE.) :',&
&               ' MAX(IJB,NJL) - MIN(IJE,NJH) ',NJINF,' - ',NJSUP
                print *,' Si vous voulez selectionner les limites en J, mettez : ',&
&               'L2DBY=.TRUE.' 
                print *,' et definissez : NJDEBCOU=    NLMAX= '
              ELSE
	        NIINF=1; NISUP=1
                NJINF=NJDEBCOU; NJSUP=NJDEBCOU+NLMAX-1
                NJINF=MAX(NJINF,NJL);NJSUP=MIN(NJSUP,NJH)
              ENDIF
	    ELSE
              IF(.NOT.L2DBY)THEN
	        NIINF=IIB; NISUP=NIINF
	        NJINF=MAX(IJB,NJL); NJSUP=MIN(IJE,NJH)
	        print *,' Profil horizontal // Y '
	        print *,' Limites J par defaut (L2DBY=.FALSE.) :',&
&               ' MAX(IJB,NJL) - MIN(IJE,NJH) ',NJINF,' - ',NJSUP
                print *,' Si vous voulez selectionner les limites en J, mettez : ',&
&               'L2DBY=.TRUE.' 
                print *,' et definissez : NJDEBCOU=    NLMAX= '
              ELSE
                NIINF=IIB; NISUP=NIINF
                NJINF=NJDEBCOU; NJSUP=NJDEBCOU+NLMAX-1
                NJINF=MAX(NJINF,NJL);NJSUP=MIN(NJSUP,NJH)
              ENDIF
	    ENDIF
	    ILENW=NJSUP-NJINF+1

	    ALLOCATE(ZWORK1D(ILENW),ZWORKY(ILENW))

	    DO JLOOPP=1,NBPROCDIA(KLOOP)
	      NLOOPP=NPROCDIA(JLOOPP,KLOOP)

	      YTITX(1:LEN(YTITX))=' '
	      YTITY(1:LEN(YTITY))=' '

		  CALL LOADUNITIT(JLOOPP,KLOOP)

	      YTITX='Y(M)'
	      YTITY=CUNITGAL(1:LEN_TRIM(CUNITGAL))

	      ZWORK1D(:)=0.; ZWORKY(:)=0.
	      IF(.NOT.LTINCRDIA(KLOOP,1))THEN
		
		DO JLOOPT=1,NBTIMEDIA(KLOOP,1)

		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)

		  IF(LPYT)THEN
		    IF(JLOOPT == 1)THEN
		      ILENW=NBTIMEDIA(KLOOP,1)
		      IX=NJSUP-NJINF+1
		      ALLOCATE(ZTEM2D(IX,ILENW))
		      ALLOCATE(ZWORKT(ILENW))
		      ZTEM2D=XSPVAL
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,ILENW))
		      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ENDIF

                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ZWORKT(JLOOPT)=XTRAJT(NLOOPT,1)
		    ZTEM2D(:,JLOOPT)=XVAR(NIINF,NJINF-NJL+1:NJSUP-NJL+1,1, &
		    NLOOPT,1,NLOOPP)
		    IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
		    IF(.NOT.LPBREAD)THEN
		      IF(KLOOP == NSUPERDIA)THEN
		        CALL NGPICT(1,1)
		        CALL GQACWK(1,IER,INB,IWK)
		        IF(INB > 1)CALL NGPICT(2,3)
		      ENDIF
		    ENDIF
		    DEALLOCATE(ZTEM2D,ZWORKT)
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			DEALLOCATE(XPRDAT)
		      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ENDIF

		  ELSE

		    ZWORK1D=XXY(NJINF:NJSUP,NMGRID)
		    ZWORKY=XVAR(NIINF,NJINF-NJL+1:NJSUP-NJL+1,1,NTIMEDIA(JLOOPT,KLOOP,1),1,NLOOPP)
		    ZTIMED=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
		    ZTIMEF=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
		    IF(JLOOPT == 1)THEN
		      IF(LDATFILE)CALL DATFILE_FORDIACHRO
		      CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                    ENDIF
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,1))
		        CALL LOAD_XPRDAT(1,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    CALL TRAXY(ZWORK1D,ZWORKY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			DEALLOCATE(XPRDAT)
		      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
!                   IF(KLOOP == NSUPERDIA)CALL FRAME
		    IF(KLOOP == NSUPERDIA)THEN
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
		    ENDIF

		  ENDIF
	        ENDDO

	      ELSE

		DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		  NLOOPT=JLOOPT

		  IF(LPYT)THEN
		    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		      ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
		      NTIMEDIA(3,KLOOP,1)+1
!                     print *,'oper verif ilenw ',ILENW
		      ITIMEND=NTIMEDIA(1,KLOOP,1)+(((NTIMEDIA(2,KLOOP,1)- &
		      NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
		      IX=NJSUP-NJINF+1
		      ALLOCATE(ZTEM2D(IX,ILENW))
		      ALLOCATE(ZWORKT(ILENW))
		      ZTEM2D=XSPVAL
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,ILENW))
		      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		      IJLT=0
		    ENDIF
                    IJLT=IJLT+1
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      CALL LOAD_XPRDAT(IJLT,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ZWORKT(IJLT)=XTRAJT(NLOOPT,1)
		    ZTEM2D(:,IJLT)=XVAR(NIINF,NJINF-NJL+1:NJSUP-NJL+1,1, &
		    NLOOPT,1,NLOOPP)
		    IF(JLOOPT == ITIMEND)THEN
		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
		    IF(.NOT.LPBREAD)THEN
		      IF(KLOOP == NSUPERDIA)THEN
		        CALL NGPICT(1,1)
		        CALL GQACWK(1,IER,INB,IWK)
		        IF(INB > 1)CALL NGPICT(2,3)
		      ENDIF
		    ENDIF
		    DEALLOCATE(ZTEM2D,ZWORKT)
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			DEALLOCATE(XPRDAT)
		      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    ENDIF

		  ELSE

		    ZWORK1D=XXY(NJINF:NJSUP,NMGRID)
		    ZWORKY=XVAR(NIINF,NJINF-NJL+1:NJSUP-NJL+1,1,JLOOPT,1,NLOOPP)
		    ZTIMED=XTRAJT(JLOOPT,1)
		    ZTIMEF=XTRAJT(JLOOPT,1)
                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,1))
		        CALL LOAD_XPRDAT(1,NLOOPT)
		    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    IF(JLOOPT == 1)THEN
		      IF(LDATFILE)CALL DATFILE_FORDIACHRO
		      CALL RESOLV_TIMES(JLOOPT)
                    ENDIF
		    CALL TRAXY(ZWORK1D,ZWORKY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)
!                   IF(KLOOP == NSUPERDIA)CALL FRAME
		    IF(KLOOP == NSUPERDIA)THEN
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
		    ENDIF
                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			DEALLOCATE(XPRDAT)
		      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

		  ENDIF
		ENDDO
	      ENDIF
	    ENDDO

	    DEALLOCATE(ZWORK1D,ZWORKY)

            NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP

	  ELSE IF((II /= 1 .AND. IIE /= IIB) .AND. (IJ /= 1 .AND. IJB /= IJE) .AND. IK == 1)THEN

! Cas compression bilan sur axe Z ou 2D hor.  -->  Plan horizontal
! ****************************************************************
! (Impression dim mat. modele et matrice(NIL:NIH,NJL:NJH,
!  NKL:NKH) et matrice(NIL:NIH,NJL:NJH,1)

	    LCHXY=.TRUE.
	    CALL RESOLV_NIJINF_NIJSUP

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CH  Allocation matrice 2D de reception des valeurs
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            ALLOCATE (ZTEM2D(1:NISUP-NIINF+1,1:NJSUP-NJINF+1))

! Ajout PH Oct 2000 + 1pt FT ou PVKT_k_1
	    IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
!! Nov 2001
               LDIRWM .OR. LDIRWT .OR. LDIRWIND .OR. &
!! Nov 2001
	       (LCH .AND. LCV) .OR. LFT .OR. LPVKT)THEN
	      ALLOCATE (ZWORK3D(1:NISUP-NIINF+1,1:NJSUP-NJINF+1,1))
	      IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
!! Nov 2001
               LDIRWM .OR. LDIRWT .OR. LDIRWIND )THEN 
!! Nov 2001
		NMGRID=1
	      ENDIF
	    ENDIF

	    DO JLOOPP=1,NBPROCDIA(KLOOP)      !--- LCHXY-------------
	      NLOOPP=NPROCDIA(JLOOPP,KLOOP)

		  CALL LOADUNITIT(JLOOPP,KLOOP)
              YTEXTE(1:LEN(YTEXTE)) = ' '
              ILENT=LEN_TRIM(CTITGAL)
	      ILENU=LEN_TRIM(CUNITGAL)
	      YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
	      YTEXTE(ILENT+1:ILENT+1)=' '
	      YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
           if(nverbia >0)then
             print *,' OPER TIT=',CTITGAL(1:ILENT),' UNIT=',CUNITGAL(1:ILENU),&
                     ' TEXTE=',TRIM(YTEXTE)
           endif
	      IF(.NOT.LTINCRDIA(KLOOP,1))THEN

		DO JLOOPT=1,NBTIMEDIA(KLOOP,1)

		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)

		  IF(LANIMT .AND. NISUP-NIINF /= 0 .AND. NJSUP-NJINF /= 0)THEN
		    IF(JLOOPT == 1)THEN
		      CALL FMFREE(YBID,YBID,IRESP)
                      if(nverbia >0)then
		      print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
                      endif

		      CALL FMATTR(YBID,YBID,IBID,IRESP)
		      CALL GOPWK(9,IBID,3)
!                     CALL GOPWK(9,20,3)
		      ISEGM=ISEGM+1
		      ISEGD=ISEGM
		      CALL GFLAS1(ISEGM)
                    ELSE
		      ISEGM=ISEGM+1
		      CALL GFLAS1(ISEGM)
                    ENDIF
                  ENDIF
		  IF((.NOT.LFT .AND. .NOT.LPVKT) .OR. (LFT .OR. LPVKT .OR. JLOOPT == 1))THEN
		  CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
		  ENDIF
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
! Ajout PH Oct 2000
!! Nov 2001
	          IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. LDIRWM &
                    .OR. LDIRWT .OR. LDIRWIND )THEN
!! Nov 2001
!                 IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT )THEN
		  ZWORK3D(:,:,1)=XU(NIINF-NIL+1:NISUP-NIL+1, &
				   NJINF-NJL+1:NJSUP-NJL+1, &
				   1,NTIMEDIA(JLOOPT,KLOOP,1),1,  &
				   NPROCDIA(JLOOPP,KLOOP))
                  ELSE IF((LCH .AND. LCV) .OR. LFT .OR. LPVKT)THEN
		  ZWORK3D(:,:,1)=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
				   NJINF-NJL+1:NJSUP-NJL+1, &
				   1,NTIMEDIA(JLOOPT,KLOOP,1),1,  &
				   NPROCDIA(JLOOPP,KLOOP))
		  ELSE
		  ZTEM2D(:,:)=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
				   NJINF-NJL+1:NJSUP-NJL+1, &
				   1,NTIMEDIA(JLOOPT,KLOOP,1),1,  &
				   NPROCDIA(JLOOPP,KLOOP))
		  ENDIF
		  IF(NSUPERDIA > 1)THEN
		    IF(LMINUS .OR. LPLUS)THEN
		      IF(NBPM > 1)THEN
			DO JU=1,NBPM
			  IF(NUMPM(JU) == 3)THEN
		            LSUPER=.TRUE.
			    EXIT
			  ELSE
		            LSUPER=.FALSE.
			  ENDIF
			ENDDO
		      ELSE
		        LSUPER=.FALSE.
		      ENDIF
		    ELSE
		      LSUPER=.TRUE.
		    ENDIF
		    IF(KLOOP == 1)NSUPER=0
		  ELSE
		    LSUPER=.FALSE.
		  ENDIF
		  CTYPHOR='K'

		  IF(NISUP-NIINF == 0 .OR. NJSUP-NJINF == 0)THEN

		    IF(LPXT .OR. LPYT)THEN
		      IF(JLOOPT == 1)THEN
			ILENW=NBTIMEDIA(KLOOP,1)
			IF(LPXT)THEN
			  IX=NISUP-NIINF+1
			ELSE IF(LPYT)THEN
			  IX=NJSUP-NJINF+1
			ENDIF
			ALLOCATE(ZPROVI2(IX,ILENW))
			ALLOCATE(ZWORKT(ILENW))
			ZPROVI2=XSPVAL
                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			  ALLOCATE(XPRDAT(16,ILENW))
                        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

		      ENDIF
                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			  CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
                        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		      ZWORKT(JLOOPT)=XTRAJT(NLOOPT,1)
		      IF(LPXT)THEN
			ZPROVI2(:,JLOOPT)=ZTEM2D(:,1)
		      ELSE IF(LPYT)THEN
			ZPROVI2(:,JLOOPT)=ZTEM2D(1,:)
		      ENDIF
		      IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
			CALL PVFCT(ZWORKT,ZPROVI2,KLOOP)
			IF(.NOT.LPBREAD)THEN
			  IF(KLOOP == NSUPERDIA)THEN
			    CALL NGPICT(1,1)
			    CALL GQACWK(1,IER,INB,IWK)
			    IF(INB > 1)CALL NGPICT(2,3)
		          ENDIF
                        ENDIF
		        DEALLOCATE(ZPROVI2,ZWORKT)
                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			  DEALLOCATE(XPRDAT)
                        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                      ENDIF
                      
		    ELSE
		    ALLOCATE(ZPROVI(SIZE(ZTEM2D,1),SIZE(ZTEM2D,2),1))
		    ZPROVI(:,:,1)=ZTEM2D(:,:)
                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			  ALLOCATE(XPRDAT(16,1))
			  CALL LOAD_XPRDAT(1,NLOOPT)
                        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    CALL TRACEH_FORDIACHRO(1,ZPROVI,KLOOP)
                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			  DEALLOCATE(XPRDAT)
                        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		    DEALLOCATE(ZPROVI)
		    ENDIF

		  ELSE

! Ajout PH Oct 2000
	            IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
!! Nov 2001
               LDIRWM .OR. LDIRWT .OR. LDIRWIND .OR. &
!! Nov 2001
		       (LCH .AND. LCV) .OR. LFT .OR. LPVKT)THEN

		      IF(LFT .OR. LPVKT)THEN
			ILENW=NBTIMEDIA(KLOOP,1)

			IF(JLOOPT == 1)THEN
			  ALLOCATE(ZWORKT(ILENW))
			  ALLOCATE(ZWORK1D(ILENW))
                          CALL VERIFLEN_FORDIACHRO
			  CALL MEMCV
			  IF(ALLOCATED(ZTEMCV))THEN
			    DEALLOCATE(ZTEMCV)
			  ENDIF
			  ALLOCATE(ZTEMCV(NLMAX,1))
                          IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			    IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			    ALLOCATE(XPRDAT(16,ILENW))
                          ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
                        ENDIF

			CALL PRECOU_FORDIACHRO(ZWORK3D,ZTEMCV)
			ZWORK1D(JLOOPT)=ZTEMCV(NPROFILE,1)
			ZWORKT(JLOOPT)=XTRAJT(NLOOPT,1)
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 

			IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
			  IF(LFT)THEN
			    CALL VARFCT(ZWORKT,ZWORK1D,1)
			  ELSEIF(LPVKT)THEN
			    ALLOCATE(ZPROVI2(1,SIZE(ZWORKT,1)))
			    ZPROVI2(1,:)=ZWORK1D
			    CALL PVFCT(ZWORKT,ZPROVI2,KLOOP)
			    DEALLOCATE(ZPROVI2)
			  ENDIF
			  DEALLOCATE(ZWORKT,ZWORK1D)
			  IF(ALLOCATED(ZTEMCV))THEN
			    DEALLOCATE(ZTEMCV)
			  ENDIF
                          IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			    DEALLOCATE(XPRDAT)
                          ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
			  IF(KLOOP == NSUPERDIA)THEN
			    CALL NGPICT(1,1)
			    CALL GQACWK(1,IER,INB,IWK)
			    IF(INB > 1)CALL NGPICT(2,3)
		          ENDIF
			ENDIF

		      ELSE
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			  ALLOCATE(XPRDAT(16,1))
			  CALL LOAD_XPRDAT(1,NLOOPT)
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		        CALL TRACEH_FORDIACHRO(1,ZWORK3D,KLOOP)
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  DEALLOCATE(XPRDAT)
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		      ENDIF

		    ELSE

                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,1))
			CALL LOAD_XPRDAT(1,NLOOPT)
                      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
                      CALL IMAGE_FORDIACHRO(ZTEM2D,1,XDIAINT,NHI,NDOT,YTEXTE(1:&
                                                              LEN_TRIM(YTEXTE)))
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
	                DEALLOCATE(XPRDAT)
                      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
	if(nverbia > 0)then
	  print *,' **oper AP IMAGE1 II,IJ,IK,KLOOP ',II,IJ,IK,KLOOP
	endif
                    ENDIF
		  ENDIF
                  IF(LANIMT .AND. NISUP-NIINF /= 0 .AND. NJSUP-NJINF /= 0)THEN
		    CALL GFLAS2
		    IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
		      DO JJ=ISEGD,ISEGM
			CALL GFLAS3(JJ)
		      ENDDO
		      CALL GCLWK(9)
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
                    ENDIF
		  ELSE IF(LPXT.OR.LPYT .OR. LFT .OR. LPVKT)THEN
		  ELSE
!                 IF(KLOOP == NSUPERDIA)CALL FRAME
		  IF(KLOOP == NSUPERDIA)THEN

                    ! Trace du domaine fils eventuellement
		    IF(LDOMAIN .AND. .NOT.LCV)THEN
                      ZZZXD=XXX(NDOMAINL,NMGRID)
                      ZZZXF=XXX(NDOMAINR,NMGRID)
                      ZZZYD=XXY(NDOMAINB,NMGRID)
                      ZZZYF=XXY(NDOMAINT,NMGRID)
                      CALL GSLWSC(XLWDOMAIN)
                      CALL FRSTPT(ZZZXD,ZZZYD)
                      CALL VECTOR(ZZZXF,ZZZYD)
                      CALL VECTOR(ZZZXF,ZZZYF)
                      CALL VECTOR(ZZZXD,ZZZYF)
                      CALL VECTOR(ZZZXD,ZZZYD)
		    ENDIF
                    ! Trace de segments eventuellement
		    IF(LSEGM .AND. .NOT.LCV)THEN
		      CALL GQPLCI(IER,ICOLI)
		      DO J=1,NCOLSEGM
      !IF(.NOT.LCOLAREA .AND. .NOT.LCOLINE .AND. NCOLSEGMS(J) > 1)THEN
      IF(NCOLSEGMS(J) > 1)THEN
	CALL TABCOL_FORDIACHRO
	print *,' appel a TABCOL_FORDIACHRO pour le trace de polynes'
      ENDIF
		      EXIT
		      ENDDO
		      CALL GSLWSC(XLWSEGM)
		      ISEGMCOL=0
                      if(nverbia > 0)then
                        print *,' **oper size((NSEGMS) ',size(NSEGMS)
                      endif
                      IGRID=NGRIDIA(NPROCDIA(JLOOPP,KLOOP))
		      DO J=1,SIZE(NSEGMS,1)
                      ! Conversion en coordonnees conformes
                        ZLAT=XSEGMS(J,1)
                        ZLON=XSEGMS(J,2)
                        IF (NSEGMS(J)==1) THEN           ! XSEGMS
                          IF (XCONFSEGMS(J,1)==0. .AND. XCONFSEGMS(J,2)==0.) &
                            CALL SM_XYHAT_S(XLATORI,XLONORI, &
                                            ZLAT,ZLON,                 &
                                            XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                        ELSE IF (NSEGMS(J)==-1) THEN     ! ISEGMS
                          NSEGMS(J)=1
                          II=MAX(MIN(INT(ZLAT),NIMAX+2*JPHEXT-1),1)
                          IJ=MAX(MIN(INT(ZLON),NJMAX+2*JPHEXT-1),1)
                          XCONFSEGMS(J,1)=XXX(II,IGRID) +  &
                             (ZLAT-FLOAT(II))*(XXX(II+1,IGRID) - XXX(II,IGRID) )
                          XCONFSEGMS(J,2)=XXY(IJ,IGRID) + &
                             (ZLON-FLOAT(IJ))*(XXY(IJ+1,IGRID) - XXY(IJ,IGRID) )
                        END IF
			IF(J == 1 .AND. NSEGMS(J) == 1) THEN
		          ISEGMCOL=ISEGMCOL+1
			  ICOLSEGM=NCOLSEGMS(ISEGMCOL)
		      IF((LCOLAREA .OR. LCOLINE) .AND. ICOLSEGM > 1)THEN
	print *,' Avec LCOLAREA=T ou LCOLINE=T , attention a la superposition des couleurs'
	!print *,' valeur trouvee: ',NCOLSEGMS,'FORCEE a 1 '
        print *,' pour les segments preferez NCOLSEGMS= 0 ou 1 '
		       !ICOLSEGM=1
		      ENDIF
		          CALL GSPLCI(ICOLSEGM)
		          CALL GSTXCI(ICOLSEGM)
                          CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
!!!!!
			ELSE IF(J > 1 .AND. NSEGMS(J) == 1 )THEN
                          IF( NSEGMS(J-1) == 0)THEN
                            ISEGMCOL=ISEGMCOL+1
                            ICOLSEGM=NCOLSEGMS(ISEGMCOL)
                            IF(J > 1)CALL SFLUSH
                      IF((LCOLAREA .OR. LCOLINE) .AND. ICOLSEGM > 1)THEN
	print *,' Avec LCOLAREA=T ou LCOLINE=T , attention a la superposition des couleurs'
        !print *,' valeur trouvee: ',NCOLSEGMS,'FORCEE a 1 '
        print *,' pour les segments preferez NCOLSEGMS= 0 ou 1 '
                       !ICOLSEGM=1
                      ENDIF
                            CALL GSPLCI(ICOLSEGM)
                            CALL GSTXCI(ICOLSEGM)
                            CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          ELSEIF(NSEGMS(J-1)== 1)THEN
                            CALL VECTOR(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          ENDIF
!!!!!
			ENDIF
		      ENDDO
		      CALL SFLUSH
		      CALL GSPLCI(ICOLI)
		      CALL GSTXCI(1)
		    ENDIF
                    ! Trace de la CV dans CH suivante(s) eventuellement
		    IF(LTRACECV .AND. .NOT.LCV)THEN
		      CALL GQLWSC(IER,ZLW)
		      CALL GSLWSC(XLWTRACECV)
		      CALL GSMKSC(2.)
                      if(nverbia > 0)then
                        print *,' **oper size((NSEGMS) for tracecv',size(NSEGMS)
                      endif
                      DO J=1,SIZE(NSEGMS,1)
                        ICOLSEGM=1
			IF(J == 1 .AND. NSEGMS(J) == 2) THEN
		          CALL GSPLCI(ICOLSEGM)
		          CALL GSTXCI(ICOLSEGM)
                          CALL GSMK(4)
                          CALL GPM(1,XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                        ELSE IF(J > 1 .AND. NSEGMS(J) == 2 )THEN
                          IF( NSEGMS(J-1) == 0)THEN
                            CALL SFLUSH
                            CALL GSPLCI(ICOLSEGM)
                            CALL GSTXCI(ICOLSEGM)
                            CALL GSMK(4)
                            CALL GPM(1,XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                            CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          ELSEIF(NSEGMS(J-1)== 2)THEN
                            CALL GSMK(5)
                            CALL GPM(1,XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                            CALL VECTOR(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          ENDIF
			ENDIF
		      ENDDO
		      CALL SFLUSH
		      CALL GSLWSC(ZLW)
		      CALL GSTXCI(1)
		    ENDIF
                    !
		    CALL NGPICT(1,1)
		    CALL GQACWK(1,IER,INB,IWK)
		    IF(INB > 1)CALL NGPICT(2,3)
		  ENDIF
		  ENDIF
		ENDDO
	      ELSE
		DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),  &
			  NTIMEDIA(3,KLOOP,1)
		  NLOOPT=JLOOPT
		  IF(LANIMT .AND. NJSUP-NJINF /= 0 .AND. NISUP-NIINF /=0)THEN
		    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		      CALL FMFREE(YBID,YBID,IRESP)
                      if(nverbia >0)then
		      print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
                      endif
		      CALL FMATTR(YBID,YBID,IBID,IRESP)
		      CALL GOPWK(9,IBID,3)
		      ISEGM=ISEGM+1
		      ISEGD=ISEGM
		      CALL GFLAS1(ISEGM)
		      ITIMEND=NTIMEDIA(1,KLOOP,1) + &
		      (((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/  &
		      NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
                    ELSE
		      ISEGM=ISEGM+1
		      CALL GFLAS1(ISEGM)
                    ENDIF
		  ENDIF
		  IF((.NOT.LFT .AND. .NOT.LPVKT) .OR. (LFT .OR. LPVKT .OR. JLOOPT == NTIMEDIA(1,KLOOP,1)))THEN
		    CALL RESOLV_TIMES(JLOOPT)
                  ENDIF
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)

! Ajout PH Oct 2000
!! Nov 2001
!                 IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT)THEN
                  IF(LDIRWM .OR. LDIRWT .OR. LDIRWIND .OR. &
	          LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT)THEN
!! Nov 2001
		    ZWORK3D(:,:,1)=XU(NIINF-NIL+1:NISUP-NIL+1, &
				   NJINF-NJL+1:NJSUP-NJL+1, &
				   1,JLOOPT,1,NPROCDIA(JLOOPP,KLOOP))
                  ELSEIF((LCH .AND. LCV) .OR. LFT .OR.LPVKT)THEN
		    ZWORK3D(:,:,1)=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
				   NJINF-NJL+1:NJSUP-NJL+1, &
				   1,JLOOPT,1,NPROCDIA(JLOOPP,KLOOP))

                  ELSE
		    ZTEM2D(:,:)=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
				   NJINF-NJL+1:NJSUP-NJL+1, &
				   1,JLOOPT,1,NPROCDIA(JLOOPP,KLOOP))
                  ENDIF
		  IF(NSUPERDIA > 1)THEN
!	    LSUPER=.TRUE.
		    IF(LMINUS .OR. LPLUS)THEN
		      IF(NBPM > 1)THEN
			DO JU=1,NBPM
			  IF(NUMPM(JU) == 3)THEN
		            LSUPER=.TRUE.
			    EXIT
			  ELSE
		            LSUPER=.FALSE.
			  ENDIF
			ENDDO
		      ELSE
		        LSUPER=.FALSE.
		      ENDIF
		    ELSE
		      LSUPER=.TRUE.
		    ENDIF
		    IF(KLOOP == 1)NSUPER=0
		  ELSE
		    LSUPER=.FALSE.
		  ENDIF
		  CTYPHOR='K'
		  IF(NISUP-NIINF == 0 .OR. NJSUP-NJINF == 0)THEN
		    IF(LPXT .OR. LPYT)THEN
		      IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
			ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))&
			       /NTIMEDIA(3,KLOOP,1)+1
	                IF(NVERBIA > 0)THEN
                        print *,'oper verif ilenw ',ILENW
                        ENDIF
			ITIMEND=NTIMEDIA(1,KLOOP,1)+(((NTIMEDIA(2,KLOOP,1)- &
			NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1))* &
			NTIMEDIA(3,KLOOP,1))
			IF(LPXT)THEN
			  IX=NISUP-NIINF+1
			ELSE IF(LPYT)THEN
			  IX=NJSUP-NJINF+1
			ENDIF
			ALLOCATE(ZPROVI2(IX,ILENW))
			ALLOCATE(ZWORKT(ILENW))
			ZPROVI2=XSPVAL
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			  ALLOCATE(XPRDAT(16,ILENW))
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
			IJLT=0
		      ENDIF
		      IJLT=IJLT+1
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			CALL LOAD_XPRDAT(IJLT,NLOOPT)
                      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		      ZWORKT(IJLT)=XTRAJT(NLOOPT,1)
		      IF(LPXT)THEN
			ZPROVI2(:,IJLT)=ZTEM2D(:,1)
		      ELSE IF(LPYT)THEN
			ZPROVI2(:,IJLT)=ZTEM2D(1,:)
		      ENDIF
		      IF(JLOOPT == ITIMEND)THEN
			CALL PVFCT(ZWORKT,ZPROVI2,KLOOP)
			IF(.NOT.LPBREAD)THEN
			  IF(KLOOP == NSUPERDIA)THEN
			    CALL NGPICT(1,1)
			    CALL GQACWK(1,IER,INB,IWK)
			    IF(INB > 1)CALL NGPICT(2,3)
		          ENDIF
                        ENDIF
		        DEALLOCATE(ZPROVI2,ZWORKT)
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  DEALLOCATE(XPRDAT)
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
                      ENDIF
		    ELSE
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		        IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		        ALLOCATE(XPRDAT(16,1))
		        CALL LOAD_XPRDAT(1,NLOOPT)
                      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		      ALLOCATE(ZPROVI(SIZE(ZTEM2D,1),SIZE(ZTEM2D,2),1))
		      ZPROVI(:,:,1)=ZTEM2D(:,:)
		      CALL TRACEH_FORDIACHRO(1,ZPROVI,KLOOP)
		      DEALLOCATE(ZPROVI)
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		        DEALLOCATE(XPRDAT)
                      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    ENDIF
		  ELSE
! Ajout PH Oct 2000 + Nov FT ou PVKT
	            IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
!! Nov 2001
               LDIRWM .OR. LDIRWT .OR. LDIRWIND .OR. &
!! Nov 2001
		       (LCH .AND. LCV ) .OR. LFT .OR. LPVKT)THEN

                      IF(LFT .OR. LPVKT)THEN
			ILENW=(NTIMEDIA(2,KLOOP,1)- &
			NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1)+1

			IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
			  ALLOCATE(ZWORKT(ILENW))
			  ALLOCATE(ZWORK1D(ILENW))
                          CALL VERIFLEN_FORDIACHRO
			  CALL MEMCV
			  IF(ALLOCATED(ZTEMCV))THEN
			    DEALLOCATE(ZTEMCV)
			  ENDIF
			  ALLOCATE(ZTEMCV(NLMAX,1))
                          IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			    IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			    ALLOCATE(XPRDAT(16,ILENW))
                          ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
			  ILT=0
                        ENDIF

			CALL PRECOU_FORDIACHRO(ZWORK3D,ZTEMCV)
			ILT=ILT+1
			ZWORK1D(ILT)=ZTEMCV(NPROFILE,1)
			ZWORKT(ILT)=XTRAJT(NLOOPT,1)
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  CALL LOAD_XPRDAT(ILT,NLOOPT)
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 

			IF(JLOOPT == NTIMEDIA(2,KLOOP,1))THEN
			  IF(LFT)THEN
			  CALL VARFCT(ZWORKT,ZWORK1D,1)
			  ELSEIF(LPVKT)THEN
			    ALLOCATE(ZPROVI2(1,SIZE(ZWORKT,1)))
			    ZPROVI2(1,:)=ZWORK1D
			    CALL PVFCT(ZWORKT,ZPROVI2,KLOOP)
			    DEALLOCATE(ZPROVI2)
			  ENDIF
			  DEALLOCATE(ZWORKT,ZWORK1D)
			  IF(ALLOCATED(ZTEMCV))THEN
			    DEALLOCATE(ZTEMCV)
			  ENDIF
                          IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			    DEALLOCATE(XPRDAT)
                          ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
			  IF(KLOOP == NSUPERDIA)THEN
			    CALL NGPICT(1,1)
			    CALL GQACWK(1,IER,INB,IWK)
			    IF(INB > 1)CALL NGPICT(2,3)
		          ENDIF
			ENDIF

		      ELSE
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			  ALLOCATE(XPRDAT(16,1))
			  CALL LOAD_XPRDAT(1,NLOOPT)
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		        CALL TRACEH_FORDIACHRO(1,ZWORK3D,KLOOP)
                        IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			  DEALLOCATE(XPRDAT)
                        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		      ENDIF

		    ELSE
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
	                IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		        ALLOCATE(XPRDAT(16,1))
		        CALL LOAD_XPRDAT(1,NLOOPT)
                      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		  if(nverbia >0)THEN
		    print *,' ** oper appel image  Ytexte ',YTEXTE(1:LEN_TRIM(YTEXTE))
		  endif
                      CALL IMAGE_FORDIACHRO(ZTEM2D,1,XDIAINT,NHI,NDOT,YTEXTE(1: &
                                                              LEN_TRIM(YTEXTE)))
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		        DEALLOCATE(XPRDAT)
                      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
	if(nverbia > 0)then
	  print *,' **oper AP IMAGE2 II,IJ,IK,KLOOP ',II,IJ,IK,KLOOP
	endif
		    ENDIF
		  ENDIF
		  IF(LANIMT .AND. NISUP-NIINF /= 0 .AND. NJSUP-NJINF /= 0)THEN
		    CALL GFLAS2
		    IF(JLOOPT == ITIMEND)THEN
		      DO JJ=ISEGD,ISEGM
                        CALL GFLAS3(JJ)
                      ENDDO 
		      CALL GCLWK(9)
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
		    ENDIF
		  ELSE IF(LPXT.OR.LPYT .OR. LFT .OR. LPVKT)THEN
		  ELSE
!                 IF(KLOOP == NSUPERDIA)CALL FRAME
		  IF(KLOOP == NSUPERDIA)THEN
                    ! Trace du domaine fils eventuellement
		    IF(LDOMAIN .AND. .NOT.LCV)THEN
                      ZZZXD=XXX(NDOMAINL,NMGRID)
                      ZZZXF=XXX(NDOMAINR,NMGRID)
                      ZZZYD=XXY(NDOMAINB,NMGRID)
                      ZZZYF=XXY(NDOMAINT,NMGRID)
                      CALL GSLWSC(XLWDOMAIN)
                      CALL FRSTPT(ZZZXD,ZZZYD)
                      CALL VECTOR(ZZZXF,ZZZYD)
                      CALL VECTOR(ZZZXF,ZZZYF)
                      CALL VECTOR(ZZZXD,ZZZYF)
                      CALL VECTOR(ZZZXD,ZZZYD)
		    ENDIF
                    ! Trace de segments eventuellement
		    IF(LSEGM .AND. .NOT.LCV)THEN
		      CALL GQPLCI(IER,ICOLI)
		      ICOLSEGM=NCOLSEGMS(1)
		      DO J=1,NCOLSEGM
      !IF(.NOT.LCOLAREA .AND. .NOT.LCOLINE .AND. NCOLSEGMS(J) > 1)THEN
      IF(NCOLSEGMS(J) > 1)THEN
	CALL TABCOL_FORDIACHRO
	print *,' appel a TABCOL_FORDIACHRO pour le trace de polynes'
      ENDIF
		      EXIT
		      ENDDO
		      CALL GSLWSC(XLWSEGM)
		      ISEGMCOL=0
                      if(nverbia > 0)then
                        print *,' **oper size2(NSEGMS) ',size(NSEGMS)
                      endif
                      IGRID=NGRIDIA(NPROCDIA(JLOOPP,KLOOP))
		      DO J=1,SIZE(NSEGMS,1)
                      ! Conversion en coordonnees conformes
                        ZLAT=XSEGMS(J,1)
                        ZLON=XSEGMS(J,2)
                        IF (NSEGMS(J)==1) THEN           ! XSEGMS
                          IF (XCONFSEGMS(J,1)==0. .AND. XCONFSEGMS(J,2)==0.) &
                            CALL SM_XYHAT_S(XLATORI,XLONORI, &
                                            ZLAT,ZLON,                 &
                                            XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                        ELSE IF (NSEGMS(J)==-1) THEN     ! ISEGMS
                          NSEGMS(J)=1
                          II=MAX(MIN(INT(ZLAT),NIMAX+2*JPHEXT-1),1)
                          IJ=MAX(MIN(INT(ZLON),NJMAX+2*JPHEXT-1),1)
                          XCONFSEGMS(J,1)=XXX(II,IGRID) +  &
                             (ZLAT-FLOAT(II))*(XXX(II+1,IGRID) - XXX(II,IGRID) )
                          XCONFSEGMS(J,2)=XXY(IJ,IGRID) + &
                             (ZLON-FLOAT(IJ))*(XXY(IJ+1,IGRID) - XXY(IJ,IGRID) )
                        END IF
			IF(J == 1 .AND. NSEGMS(J) == 1)THEN
			  ISEGMCOL=ISEGMCOL+1
			  ICOLSEGM=NCOLSEGMS(ISEGMCOL)
		      IF((LCOLAREA .OR. LCOLINE) .AND. ICOLSEGM > 1)THEN
	print *,' Avec LCOLAREA=T ou LCOLINE=T ,  attention a la superposition des couleurs'
	!print *,' valeur trouvee: ',NCOLSEGMS,'FORCEE a 1 '
        print *,' pour les segments preferez NCOLSEGMS= 0 ou 1 '
			!ICOLSEGM=1
		      ENDIF
		          CALL GSPLCI(ICOLSEGM)
		          CALL GSTXCI(ICOLSEGM)
                          CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
			ELSE IF(J > 1 .AND. NSEGMS(J) == 1 )THEN
                          IF(NSEGMS(J-1) == 0)THEN
                            ISEGMCOL=ISEGMCOL+1
                            ICOLSEGM=NCOLSEGMS(ISEGMCOL)
                            IF(J > 1)CALL SFLUSH
                      IF((LCOLAREA .OR. LCOLINE) .AND. ICOLSEGM > 1)THEN
	print *,' Avec LCOLAREA=T ou LCOLINE=T ,  attention a la superposition des couleurs'
        !print *,' valeur trouvee: ',NCOLSEGMS,'FORCEE a 1 '
        print *,' pour les segments preferez NCOLSEGMS= 0 ou 1 '
                        !ICOLSEGM=1
                      ENDIF
                            CALL GSPLCI(ICOLSEGM)
                            CALL GSTXCI(ICOLSEGM)
                            CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))

			  ELSEIF(NSEGMS(J-1)== 1)THEN
                            CALL VECTOR(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          ENDIF
			ENDIF
		      ENDDO
		      CALL SFLUSH
		      CALL GSPLCI(ICOLI)
		      CALL GSTXCI(1)
		    ENDIF
		     ! Trace de la CV dans CH suivante(s) eventuellement
		    IF(LTRACECV .AND. .NOT.LCV)THEN
		      CALL GQLWSC(IER,ZLW)
		      CALL GSLWSC(XLWTRACECV)
		      CALL GSMKSC(2.)
                      if(nverbia > 0)then
                        print *,' **oper size((NSEGMS) for tracecv2',size(NSEGMS)
                      endif
                      DO J=1,SIZE(NSEGMS,1)
                        ICOLSEGM=1
			IF(J == 1 .AND. NSEGMS(J) == 2) THEN
		          CALL GSPLCI(ICOLSEGM)
		          CALL GSTXCI(ICOLSEGM)
                          CALL GSMK(4)
                          CALL GPM(1,XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                        ELSE IF(J > 1 .AND. NSEGMS(J) == 2 )THEN
                          IF( NSEGMS(J-1) == 0)THEN
                            CALL SFLUSH
                            CALL GSPLCI(ICOLSEGM)
                            CALL GSTXCI(ICOLSEGM)
                            CALL GSMK(4)
                            CALL GPM(1,XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                            CALL FRSTPT(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          ELSEIF(NSEGMS(J-1)== 2)THEN
                            CALL GSMK(5)
                            CALL GPM(1,XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                            CALL VECTOR(XCONFSEGMS(J,1),XCONFSEGMS(J,2))
                          ENDIF
			ENDIF
		      ENDDO
                      CALL SFLUSH
		      CALL GSLWSC(ZLW)
		      CALL GSTXCI(1)
		    ENDIF
                    CALL NGPICT(1,1)
		    CALL GQACWK(1,IER,INB,IWK)
		    IF(INB > 1)CALL NGPICT(2,3)
		  ENDIF
		  ENDIF
                ENDDO
	      ENDIF
	    ENDDO                             !--- LCHXY-------------
	    DEALLOCATE(ZTEM2D)
	    IF(ALLOCATED(ZWORK3D))THEN
	      DEALLOCATE(ZWORK3D)
	    ENDIF

	  ELSE IF(II /= 1 .AND. (IJ == 1 .OR. IJE-IJB == 0) .AND. IK == 1)THEN

! Cas compression bilan sur axes Y et Z -->  Profil horizontal // X
! *****************************************************************
! (Impression dim mat. modele et matrice(NIL:NIH,NJL:NJH,
!  NKL:NKH) et matrice(NIL:NIH,1,1)

            print *,'  Profil horizontal // X'
            IINF=NIINF;ISUP=NISUP;IJINF=NJINF;IJSUP=NJSUP
	    print *,'IINF,ISUP,IJINF,IJSUP ',IINF,ISUP,IJINF,IJSUP
	    IF(IJ == 1)THEN
	      GIJ1=.TRUE.
	    ELSE
	      GIJ1=.FALSE.
              LCH=.FALSE.
	    ENDIF

	    IF(GIJ1)THEN
              IF(.NOT.L2DBX)THEN
	        NJINF=1; NJSUP=1
	        NIINF=MAX(IIB,NIL); NISUP=MIN(IIE,NIH)
                print *,' Limites I par defaut (L2DBX=.FALSE.) :',&
&             ' MAX(IIB,NIL) - MIN(IIE,NIH) ',NIINF,' - ',NISUP
                print *,' Si vous voulez selectionner les limites en I, mettez :',&
&             ' L2DBX=.TRUE.'
                print *,' et definissez : NIDEBCOU=    NLMAX= '
              ELSE
                NJINF=1;NJSUP=1
                NIINF=NIDEBCOU; NISUP=NIDEBCOU+NLMAX-1
                NIINF=MAX(NIINF,NIL);NISUP=MIN(NISUP,NIH)
              ENDIF
	    ELSE
              IF(.NOT.L2DBX)THEN
	        NJINF=IJB; NJSUP=IJE
	        NIINF=MAX(IIB,NIL); NISUP=MIN(IIE,NIH)
                print *,' Limites I par defaut (L2DBX=.FALSE.) :',&
&             ' MAX(IIB,NIL) - MIN(IIE,NIH) ',NIINF,' - ',NISUP
                print *,' Si vous voulez selectionner les limites en I, mettez :',&
&             ' L2DBX=.TRUE.'
                print *,' et definissez : NIDEBCOU=    NLMAX= '
              ELSE
	        NJINF=IJB; NJSUP=IJE
                NIINF=NIDEBCOU; NISUP=NIDEBCOU+NLMAX-1
                NIINF=MAX(NIINF,NIL);NISUP=MIN(NISUP,NIH)
              ENDIF
	    ENDIF
	    ILENW=NISUP-NIINF+1

	    ALLOCATE(ZWORK1D(ILENW),ZWORKY(ILENW))

	    DO JLOOPP=1,NBPROCDIA(KLOOP)
	      NLOOPP=NPROCDIA(JLOOPP,KLOOP)

	      YTITX(1:LEN(YTITX))=' '
	      YTITY(1:LEN(YTITY))=' '

		  CALL LOADUNITIT(JLOOPP,KLOOP)

	      YTITX='X(M)'
	      YTITY=CUNITGAL(1:LEN_TRIM(CUNITGAL))

	      ZWORK1D(:)=0.; ZWORKY(:)=0.
	      IF(.NOT.LTINCRDIA(KLOOP,1))THEN
		
		DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)

		  IF(LPXT)THEN
		    IF(JLOOPT == 1)THEN
		      ILENW=NBTIMEDIA(KLOOP,1)
		      IX=NISUP-NIINF+1
		      ALLOCATE(ZTEM2D(IX,ILENW))
		      ALLOCATE(ZWORKT(ILENW))
		      ZTEM2D=XSPVAL
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,ILENW))
		      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    ENDIF
		    ZWORKT(JLOOPT)=XTRAJT(NLOOPT,1)
                    IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		      CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
		    ENDIF !  Juin 2001 Ajout des dates ds FICVAL 

		    ZTEM2D(:,JLOOPT)=XVAR(NIINF-NIL+1:NISUP-NIL+1,NJINF,1, &
		    NLOOPT,1,NLOOPP)
		    IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
		      IF(.NOT.LPBREAD)THEN
		        IF(KLOOP == NSUPERDIA)THEN
		          CALL NGPICT(1,1)
		          CALL GQACWK(1,IER,INB,IWK)
		          IF(INB > 1)CALL NGPICT(2,3)
		        ENDIF
		      ENDIF
		      DEALLOCATE(ZTEM2D,ZWORKT)
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		        DEALLOCATE(XPRDAT)
		      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    ENDIF

		  ELSE

		    ZWORK1D=XXX(NIINF:NISUP,NMGRID)
		    ZWORKY=XVAR(NIINF-NIL+1:NISUP-NIL+1,NJINF,1,NTIMEDIA(JLOOPT,KLOOP,1),1,NLOOPP)
		    ZTIMED=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
		    ZTIMEF=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,1))
			CALL LOAD_XPRDAT(1,NLOOPT)
		      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    IF(JLOOPT == 1)THEN
		      IF(LDATFILE)CALL DATFILE_FORDIACHRO
		      CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                    ENDIF
		    CALL TRAXY(ZWORK1D,ZWORKY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			DEALLOCATE(XPRDAT)
		      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    IF(KLOOP == NSUPERDIA)THEN
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
		    ENDIF
		  ENDIF
	        ENDDO

	      ELSE

		DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		  NLOOPT=JLOOPT
		  IF(LPXT)THEN

		    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		      ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
		      NTIMEDIA(3,KLOOP,1)+1
	              IF(NVERBIA > 0)THEN
                      print *,'oper verif ilenw ',ILENW
		      ENDIF
		      ITIMEND=NTIMEDIA(1,KLOOP,1)+(((NTIMEDIA(2,KLOOP,1)- &
		      NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
		      IX=NISUP-NIINF+1
		      ALLOCATE(ZTEM2D(IX,ILENW))
		      ALLOCATE(ZWORKT(ILENW))
		      ZTEM2D=XSPVAL
                      IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
			IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			ALLOCATE(XPRDAT(16,ILENW))
		      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		      IJLT=0
		    ENDIF
                    IJLT=IJLT+1
                    IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		      CALL LOAD_XPRDAT(IJLT,NLOOPT)
		    ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    ZWORKT(IJLT)=XTRAJT(NLOOPT,1)
		    ZTEM2D(:,IJLT)=XVAR(NIINF-NIL+1:NISUP-NIL+1,NJINF,1, &
		    NLOOPT,1,NLOOPP)
		    IF(JLOOPT == ITIMEND)THEN
		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
		    IF(.NOT.LPBREAD)THEN
		      IF(KLOOP == NSUPERDIA)THEN
		        CALL NGPICT(1,1)
		        CALL GQACWK(1,IER,INB,IWK)
		        IF(INB > 1)CALL NGPICT(2,3)
		      ENDIF
		    ENDIF
                    IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		      DEALLOCATE(XPRDAT)
		    ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    DEALLOCATE(ZTEM2D,ZWORKT)
		    ENDIF

		  ELSE

		    ZWORK1D=XXX(NIINF:NISUP,NMGRID)
		    ZWORKY=XVAR(NIINF-NIL+1:NISUP-NIL+1,NJINF,1,JLOOPT,1,NLOOPP)
		    ZTIMED=XTRAJT(JLOOPT,1)
		    ZTIMEF=XTRAJT(JLOOPT,1)
		    IF(JLOOPT == 1)THEN
		      IF(LDATFILE)CALL DATFILE_FORDIACHRO
		      CALL RESOLV_TIMES(JLOOPT)
                    ENDIF
                    IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		      ALLOCATE(XPRDAT(16,1))
		      CALL LOAD_XPRDAT(1,NLOOPT)
		    ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    CALL TRAXY(ZWORK1D,ZWORKY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)
                    IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		      DEALLOCATE(XPRDAT)
		    ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		    IF(KLOOP == NSUPERDIA)THEN
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
		    ENDIF

		  ENDIF

		ENDDO
	      ENDIF
	    ENDDO

	    DEALLOCATE(ZWORK1D,ZWORKY)

            NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP

	  ELSE IF(II /= 1 .AND. IJ == 1 .AND. IK /= 1 .AND. LJCP)THEN

! Cas compression bilan sur axe Y -->  Plan vertical // X
! *******************************************************
! (Impression dim mat. modele et matrice(NIL:NIH,NJL:NJH,
!  NKL:NKH) et matrice(NIL:NIH,1,NKL:NKH)
	      IDEFCV=0                      !%%%%%%%%%%%%%%%%%%%%%%%%%%
	      IF(LDEFCV2CC)THEN
	        LDEFCV2CC=.FALSE.
	        IDEFCV=1
	      ENDIF                         !%%%%%%%%%%%%%%%%%%%%%%%%%%
            LCVXZ=.TRUE.
	    IF(.NOT.L2DBX)THEN
	      IINF=MAX(IIB,NIL)
	      ISUP=MIN(IIE,NIH)
	      print *,' 2D Vertical // X '
	      print *,' Limites I par defaut (L2DBX=.FALSE.)(par / au domaine integral de simulation,points de garde compris) :',&
&             ' MAX(IIB,NIL) - MIN(IIE,NIH) ',IINF,' - ',ISUP
              print *,' Si vous voulez selectionner les limites en I, mettez : ',&
&             'L2DBX=.TRUE.' 
              print *,' et definissez : NIDEBCOU=    NLMAX= '
	    ELSE
	      IINF=NIDEBCOU     
	      ISUP=NIDEBCOU+NLMAX-1
	      ISUP=MIN(ISUP,NIH)
	    ENDIF
	    ALLOCATE(ZTEM2D(1:ISUP-IINF+1,1:IKU))
	    NINX=ISUP-IINF+1
	    NINY=IKU
            NLMAX=NINX
            NLANGLE=0
            NIDEBCOU=IINF
	    IJDEBCOU=-999
	    IF(NJDEBCOU /= NJL)THEN
	      IJDEBCOU=NJDEBCOU
              NJDEBCOU=NJL
	      print *,' NJDEBCOU force a la valeur de NJL ',NJL,' pour ', &
&            'obtention altitudes correctes '
	      print *,' AP utilisation, sera remis a la valeur precedente : ', &
	      IJDEBCOU
            ENDIF
	    LVERT=.TRUE.
	    LHOR=.FALSE.
	    LPT=LPXT
	    IF(NSUPERDIA > 1)THEN
!      LSUPER=.TRUE.
		    IF(LMINUS .OR. LPLUS)THEN
		      IF(NBPM > 1)THEN
			DO JU=1,NBPM
			  IF(NUMPM(JU) == 3)THEN
		            LSUPER=.TRUE.
			    EXIT
			  ELSE
		            LSUPER=.FALSE.
			  ENDIF
			ENDDO
		      ELSE
		        LSUPER=.FALSE.
		      ENDIF
		    ELSE
		      LSUPER=.TRUE.
		    ENDIF
            ELSE
	      LSUPER=.FALSE.
	    ENDIF
	    IF(KLOOP == 1)NSUPER=0
	    DO JLOOPP=1,NBPROCDIA(KLOOP)      !--- LCVXZ-------------
	      NLOOPP=NPROCDIA(JLOOPP,KLOOP)

		  CALL LOADUNITIT(JLOOPP,KLOOP)

              ILENT=LEN_TRIM(CTITGAL)
	      ILENU=LEN_TRIM(CUNITGAL)
	      YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
	      YTEXTE(ILENT+1:ILENT+1)=' '
	      YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
	      IF(.NOT.LTINCRDIA(KLOOP,1))THEN
		DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
		  CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
      	          IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 0))THEN
!                   print *,' OPER LJCP .AND. SIZE(XZS,2) ',LJCP,SIZE(XZS,2)
		    IF(.NOT.LJCP .AND. SIZE(XZS,2) == 3)THEN
		      CALL COMPCOORD_FORDIACHRO(NMGRID)
		      IF(ALLOCATED(XWORKZ))THEN
			DEALLOCATE(XWORKZ)
		      ENDIF
		      IF(ALLOCATED(XDS))THEN
			DEALLOCATE(XDS)
		      ENDIF
		      IF(ALLOCATED(XWZ))THEN
			DEALLOCATE(XWZ)
		      ENDIF
		      ALLOCATE(XWORKZ(NLMAX,IKU,7))
		      ALLOCATE(XWZ(NLMAX,7))
		      ALLOCATE(XDS(NLMAX+100,7))
		      XDS(1:NLMAX,NMGRID)=XXX(IINF:ISUP,NMGRID)
		      XWORKZ(1:NLMAX,1:IKU,NMGRID)=XZZ(IINF:ISUP,NJDEBCOU,1:IKU)
		      XWZ(1:NLMAX,NMGRID)=XXZS(IINF:ISUP,2,NMGRID)
		    ENDIF
		    IF(.NOT.LJCP .AND. SIZE(XZS,2) == 3)THEN
                      DO J=1,NLMAX
		      XZWORKZ(J,1:IKU)=XWORKZ(J,1:IKU,NMGRID)
                      ENDDO
		    ELSE
      	              DO J=1,NINX
      	                XZWORKZ(J,1:IKU)=XXZ(:,NMGRID)
                      ENDDO
		    ENDIF
      	            XZZDS(1:NINX)=XXX(IINF:ISUP,NMGRID)
      	            ZWL=XZZDS(1); ZWR=XZZDS(NINX)
      	            IF((XHMIN == 0. .AND. XHMAX == 0.) .OR. (XHMAX<=XHMIN))THEN
      	              XHMIN=0.
      	              XHMAX=XZWORKZ(1,IKE)
                    ENDIF
!                   print *,' OPER XHMIN XHMAX ',XHMIN,XHMAX
      	            ZWB=XHMIN; ZWT=XHMAX
      	            CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,1)
      	            CALL GSCLIP(1)
      	            CALL CPSETI('SET',0)
      	            CALL CPSETI('MAP',4)
      	          ENDIF
		  ZTEM2D=XSPVAL
		  ZTEM2D(1:ISUP-IINF+1,NKL:NKH)=XVAR( &
		  IINF-NIL+1:ISUP-NIL+1,1,:,NTIMEDIA(JLOOPT,KLOOP,1),&
		  1,NPROCDIA(JLOOPP,KLOOP))
    	          IF(NKL < IKB)THEN
    		    ZTEM2D(:,1:IKB-1)=XSPVAL
                  ENDIF
    	          IF(NKH > IKE)THEN
    		    ZTEM2D(:,IKE+1:IKU)=XSPVAL
                  ENDIF
		  if(nverbia >0)THEN
		    print *,' ** oper appel imcou  Ytexte ',YTEXTE(1:LEN_TRIM(YTEXTE))
		  endif
                  IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		    IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		    ALLOCATE(XPRDAT(16,1))
		    CALL LOAD_XPRDAT(1,NLOOPT)
		  ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
                  IF(KLOOP == 1)NSUPER=0
		  CALL IMCOU_FORDIACHRO(ZTEM2D,XDIAINT,CLEGEND,YTEXTE(1: &
                                                        LEN_TRIM(YTEXTE)))
                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		    DEALLOCATE(XPRDAT)
		  ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
!                 IF(KLOOP == NSUPERDIA)CALL FRAME
		  IF(KLOOP == NSUPERDIA)THEN
		    CALL NGPICT(1,1)
		    CALL GQACWK(1,IER,INB,IWK)
		    IF(INB > 1)CALL NGPICT(2,3)
		  ENDIF
		ENDDO
	      ELSE
		DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),  &
			  NTIMEDIA(3,KLOOP,1)
		  NLOOPT=JLOOPT
		  CALL RESOLV_TIMES(JLOOPT)
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)
      	          IF(.NOT.LSUPER .OR. (LSUPER .AND. NSUPER == 0))THEN
		    IF(.NOT.LJCP .AND. SIZE(XZS,2) == 3)THEN
		      CALL COMPCOORD_FORDIACHRO(NMGRID)
		      IF(ALLOCATED(XWORKZ))THEN
			DEALLOCATE(XWORKZ)
		      ENDIF
		      IF(ALLOCATED(XDS))THEN
			DEALLOCATE(XDS)
		      ENDIF
		      IF(ALLOCATED(XWZ))THEN
			DEALLOCATE(XWZ)
		      ENDIF
		      ALLOCATE(XWORKZ(NLMAX,IKU,7))
		      ALLOCATE(XWZ(NLMAX,7))
		      ALLOCATE(XDS(NLMAX+100,7))
		      XDS(1:NLMAX,NMGRID)=XXX(IINF:ISUP,NMGRID)
		      XWORKZ(1:NLMAX,1:IKU,NMGRID)=XZZ(IINF:ISUP,NJDEBCOU,1:IKU)
		      XWZ(1:NLMAX,NMGRID)=XXZS(IINF:ISUP,2,NMGRID)
		    ENDIF
		    IF(.NOT.LJCP .AND. SIZE(XZS,2) == 3)THEN
		      XZWORKZ(1:NLMAX,1:IKU)=XWORKZ(1:NLMAX,1:IKU,NMGRID)
		    ELSE
      	              DO J=1,NINX
      	                XZWORKZ(J,1:IKU)=XXZ(:,NMGRID)
                      ENDDO
		    ENDIF
      	            XZZDS(1:NINX)=XXX(IINF:ISUP,NMGRID)
      	            ZWL=XZZDS(1); ZWR=XZZDS(NINX)
      	            IF((XHMIN == 0. .AND. XHMAX == 0.) .OR. (XHMAX<=XHMIN))THEN
      	              XHMIN=0.
      	              XHMAX=XZWORKZ(1,IKE)
                    ENDIF
!                   print *,' OPER 2 XHMIN XHMAX ',XHMIN,XHMAX
      	            ZWB=XHMIN; ZWT=XHMAX
      	            CALL SET(.1,.9,.1,.9,ZWL,ZWR,ZWB,ZWT,1)
      	            CALL GSCLIP(1)
      	            CALL CPSETI('SET',0)
      	            CALL CPSETI('MAP',4)
      	          ENDIF
		  ZTEM2D=XSPVAL
		  ZTEM2D(1:ISUP-IINF+1,NKL:NKH)=XVAR( &
		  IINF-NIL+1:ISUP-NIL+1,1,:,JLOOPT,1,NPROCDIA(  &
		  JLOOPP,KLOOP))
    	          IF(NKL < IKB)THEN
    		    ZTEM2D(:,1:IKB-1)=XSPVAL
                  ENDIF
    	          IF(NKH > IKE)THEN
    		    ZTEM2D(:,IKE+1:IKU)=XSPVAL
                  ENDIF
		  if(nverbia >0)THEN
		    print *,' ** oper appel imcou  Ytexte ',YTEXTE(1:LEN_TRIM(YTEXTE))
		  endif
                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		    IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		    ALLOCATE(XPRDAT(16,1))
		    CALL LOAD_XPRDAT(1,NLOOPT)
		  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                  IF(KLOOP == 1)NSUPER=0
		  CALL IMCOU_FORDIACHRO(ZTEM2D,XDIAINT,CLEGEND,YTEXTE(1: &
                                                        LEN_TRIM(YTEXTE)))
                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		    DEALLOCATE(XPRDAT)
		  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
!                 IF(KLOOP == NSUPERDIA)CALL FRAME
		  IF(KLOOP == NSUPERDIA)THEN
		    CALL NGPICT(1,1)
		    CALL GQACWK(1,IER,INB,IWK)
		    IF(INB > 1)CALL NGPICT(2,3)
		  ENDIF
                ENDDO
              ENDIF
            ENDDO                             !--- LCVXZ-------------
	    DEALLOCATE(ZTEM2D)

	    IF(IJDEBCOU /= -999)THEN
	      NJDEBCOU=IJDEBCOU
	    ENDIF

	    IF(IDEFCV==1)THEN                 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	      LDEFCV2CC=.TRUE.
	      IDEFCV=0
	    ENDIF                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  ELSE

! PAS DE COMPRESSION
! ******************

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CH Positionnement NIINF, NJINF, NISUP, NJSUP
! Defaut : NIINF=MAX(IIB,NIL), NJINF=MAX(IJB,NJL), NISUP=MIN(IIE,NIH), 
!          NJSUP=MIN(IJE,NJH)
! Sinon valeurs fournies par l'utilisateur dans les limites (NIL,NJL NIH,
! NJH)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CV Positionnement NIINF, NJINF, NISUP, NJSUP
! CV Positionnement LHORIZ et LVERTI
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	        IF(LCV)THEN
                  IINF=NIINF;ISUP=NISUP;IJINF=NJINF;IJSUP=NJSUP
                  IF(IINF == 0)THEN
                    GCH=LCH
                    LCH=.TRUE.
		    LCV=.FALSE.
                    CALL RESOLV_NIJINF_NIJSUP
                    LCH=GCH
		    LCV=.TRUE.
                    IINF=NIINF;ISUP=NISUP;IJINF=NJINF;IJSUP=NJSUP
        	  ENDIF
	          if(NVERBIA > 0)THEN
	            print *,' oper CV IINF,ISUP,IJINF,IJSUP ',IINF,ISUP,IJINF,IJSUP
	          endif
                  ! fichier 1D (points de garde dupliques dans conv2dia)
                  !pour eviter de definir la localisation du profil
                  IF (NIMAX==1 .AND. NJMAX==1) THEN
                    IF(NIDEBCOU==0 .OR. NIDEBCOU==999999999) NIDEBCOU=1+JPHEXT
                    IF(NJDEBCOU==0 .OR. NJDEBCOU==999999999) NJDEBCOU=1+JPHEXT
                    IF(NLMAX==0 .OR. NLMAX==999999999) NLMAX=2
                    IF(NLANGLE==0 .OR. NLANGLE==999999999) NLANGLE=0
                    IF(NPROFILE==0 .OR. NPROFILE==999999999) NPROFILE=1
                    LPOINTG=.TRUE.
                  ENDIF
	          if(NVERBIA > 0)THEN
	            print *,' oper CV NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,PROFILE '&
	                              ,NIDEBCOU,NJDEBCOU,NLMAX,NLANGLE,NPROFILE
	          endif
	        ENDIF
		CALL RESOLV_NIJINF_NIJSUP
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CH + CV Allocation matrice 3D de reception des valeurs
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            	ALLOCATE (ZWORK3D(1:NISUP-NIINF+1,1:NJSUP-NJINF+1, &
    		                  1:NKH-NKL+1))

                if(nverbia >0)then
        	print *,' NBPROCDIA(KLOOP) ',NBPROCDIA(KLOOP)
                endif

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle externe sur les processus
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        	DO JLOOPP=1,NBPROCDIA(KLOOP)
                  NLOOPP=NPROCDIA(JLOOPP,KLOOP)

    		  IF((LPVKT .OR. LPVKT1) .AND. NSUPERDIA>1)THEN
    		    IF(NBPROCDIA(KLOOP)>1 .OR. NBLVLKDIA(KLOOP,1)>1)THEN
    		      print *,' _PVKT_ (_PVKT1_)  SUPERPOSITIONS : '
!fuji  		      print *,'         On ne peut definir de part de d''autre '&
!fuji  		      &'de _ON_ qu''1 seul processus et 1 seul niveau'
    		      print *,'         On ne peut definir de part de d''autre '
    		      print *,'de _ON_ qu''1 seul processus et 1 seul niveau'
    		      print *,' Nb de niveaux demandes   : ',NBLVLKDIA(KLOOP,1)
    		      print *,' Nb de processus demandes : ',NBPROCDIA(KLOOP)
    		      print *,' *** MODIFIEZ VOTRE DIRECTIVE *** '
    		      EXIT
    		    ENDIF
    		  ENDIF

		  IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
		     LULMWM .OR. LULTWT .OR. LSUMVM .OR. LSUTVT .OR.  &
		     LDIRWM .OR. LDIRWT .OR. &
		     LMLSUMVM .OR. LMLSUTVT)THEN
		    NMGRID=1
                  ELSE IF(LULM .OR. LULT)THEN
! Avril 99 a la demande de Joel, Nicole et les autres
		    NMGRID=1
!                   NMGRID=2
                  ELSE IF(LVTM .OR. LVTT)THEN
! Avril 99 a la demande de Joel, Nicole et les autres
		    NMGRID=1
!                   NMGRID=3
                  ELSE
		    NMGRID=NGRIDIA(NPROCDIA(JLOOPP,KLOOP))
		    IF(NGRIDIAM /= 0 .AND. (NGRIDIAM /= NMGRID))THEN
		      print *,' ****oper NMGRID Av modif ',NMGRID
		      NMGRID=NGRIDIAM
		      print *,' ****oper NMGRID mis volontairement a la valeur de NGRIDIAM ',NGRIDIAM
		    ENDIF
		  ENDIF
		  IF(NMGRID <1 .OR. NMGRID >7)THEN
		    PRINT *,' VALEUR NMGRID ABERRANTE: ',NMGRID, &
                            '        FORCEE A        :  1'
                    NMGRID=1
                  ENDIF
		  IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
		     LULMWM .OR. LULTWT .OR. LULM .OR. LULT .OR.   &
		     LVTM .OR. LVTT .OR. LSUMVM .OR. LSUTVT .OR.   &
		     LDIRWM .OR. LDIRWT .OR. &
		     LMLSUMVM .OR. LMLSUTVT)THEN
                    CTITGAL=ADJUSTL(CGROUP)
		    CUNITGAL(1:LEN(CUNITGAL))=' '
                  ELSE
		    CTITGAL=ADJUSTL(CTITRE(NPROCDIA(JLOOPP,KLOOP)))
		    CUNITGAL=ADJUSTL(CUNITE(NPROCDIA(JLOOPP,KLOOP)))
                  ENDIF
                  if(nverbia >0)then
                    print *,' ++OPER++ CTITGAL,CUNITGAL ',CTITGAL,CUNITGAL
                  endif
		  CTITGAL=ADJUSTL(CTITGAL)
		  CUNITGAL=ADJUSTL(ADJUSTR(CUNITGAL))
                  IF(INDEX(CUNITGAL,' ') /= 0)THEN
		  CUNITGAL(INDEX(CUNITGAL,' '):LEN(CUNITGAL))=' '
                  ELSE
                  IF(LEN(CUNITGAL) > 8)Then
                  print *,' **oper DES caracteres bizarres ds le champ UNITE ',&
                  &' tronque a 8 caractères '
                  CUNITGAL(9:LEN(CUNITGAL))=' '
                  ELSE
                  ENDIF
                  ENDIF
                  if(nverbia >0)then
                    print *,' ++OPER++ CTITGAL,CUNITGAL ',CTITGAL,CUNITGAL
                  endif

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les numeros de masques ou trajectoires
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!       	print *,' NBNDIA(KLOOP) ',NBNDIA(KLOOP)

        	  DO JLOOPN=1,NBNDIA(KLOOP)
                if(nverbia >0)then
        	print *,' **oper JLOOPN ',JLOOPN
                endif

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  temps (Formulation sequentielle)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        	    IF(.NOT.LTINCRDIA(KLOOP,1))THEN

!       	      print *,' NBTIMEDIA(KLOOP,1) ',NBTIMEDIA(KLOOP,1)

        	      DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		        NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
                if(nverbia >0)then
        	print *,' **oper**A JLOOPT ',JLOOPT
                endif

			IF(LANIMT)THEN
			  IF(LPVT .OR. LPVKT .OR. LPVKT1)THEN
			    print *,' ANIMATION IMPOSSIBLE avec _PVT_ ou _PVKT_ ou _PVKT1_'
			    print *,' LANIMT remis a .FALSE. '
			    LANIMT=.FALSE.
			  ELSE IF(LPV .AND. NSUPERDIA>1)THEN
			    print *,' ANIMATION IMPOSSIBLE ', &
			    &'avec _PV_ et superpositions'
			    print *,' LANIMT remis a .FALSE. '
			    print *,' mais POSSIBLE sous la forme : ',&
&                           'GPE_PV__P_1 ou GPE_PV__P_1_T_300_TO_3600 '
			    print *,' PENSER a fournir les bornes dans ',&
&       		    'XPVMIN_proc= et XPVMAX_proc= et a les activer ',& 
&                	    'avec LMNMXUSER=T '
			    print *,' Rappel : proc=nom du processus tel ',&
&                           'qu''il est enregistre '
			    LANIMT=.FALSE.
			  ELSE
  			    IF(JLOOPT == 1)THEN
		              CALL FMFREE(YBID,YBID,IRESP)
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
			      CALL FMATTR(YBID,YBID,IBID,IRESP)
  			      CALL GOPWK(9,IBID,3)
  			      ISEGM=ISEGM+1
  			      ISEGD=ISEGM
  			      CALL GFLAS1(ISEGM)
  			    ELSE
  			      ISEGM=ISEGM+1
  			      CALL GFLAS1(ISEGM)
  			    ENDIF
			  ENDIF
			ENDIF

		        CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))

                        if(nverbia > 0)then
			  print *,' **oper LULM LCH LMUMVM,LDIRWM,LDIRWIND lig 2406 ',LULM,LCH,LMUMVM,LDIRWM,LDIRWIND
			endif
		        IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
		           LULMWM .OR. LULTWT .OR. LULM .OR. LULT .OR.   &
		           LVTM .OR. LVTT .OR. LSUMVM .OR. LSUTVT .OR.   &
			   LDIRWM .OR. LDIRWT .OR. &
			   LMLSUMVM .OR. LMLSUTVT)THEN
        	           ZWORK3D=XU(NIINF-NIL+1:NISUP-NIL+1, &
		                       NJINF-NJL+1:NJSUP-NJL+1, &
        	                     :,NTIMEDIA(JLOOPT,KLOOP,1),JLOOPN, &
				     NPROCDIA(JLOOPP,KLOOP))
!!!!! Avril 99 Ajout ULM et VTM en CH
                           IF((LCH.AND.LULM).OR.(LCH.AND.LVTM).OR. &
			   (LCH.AND.LULT).OR.(LCH.AND.LVTT))THEN
			     ALLOCATE(ZWORK3V(SIZE(ZWORK3D,1), &
			     SIZE(ZWORK3D,2),SIZE(ZWORK3D,3)))
			     ALLOCATE(ZTEM1(IIU,IJU),ZTEMV(IIU,IJU))
			     ZTEM1=0.
			     ZTEMV=0.
        	             ZWORK3V=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
		                       NJINF-NJL+1:NJSUP-NJL+1, &
        	                     :,NTIMEDIA(JLOOPT,KLOOP,1),JLOOPN, &
				     NPROCDIA(JLOOPP,KLOOP))
                             DO JKLOOP=1,IKU
			       IF(JKLOOP < MAX(IKB,NKL) .OR. &
				  JKLOOP > MIN(IKE,NKH))THEN
			       ELSE
			         ZTEM1(NIINF:NISUP,NJINF:NJSUP)= &
				 ZWORK3D(:,:,JKLOOP-NKL+1)
			         ZTEMV(NIINF:NISUP,NJINF:NJSUP)= &
				 ZWORK3V(:,:,JKLOOP-NKL+1)
!!!!essai Nov 2001 pour prise en compte PH 29/11/2001 .. A suivre
                                 CALL VERIFLEN_FORDIACHRO
!!!!essai Nov 2001
				 CALL ROTA(ZTEM1,ZTEMV)
				 ZWORK3D(:,:,JKLOOP-NKL+1)=ZTEM1(NIINF:NISUP,NJINF:NJSUP)
				 ZWORK3V(:,:,JKLOOP-NKL+1)=ZTEMV(NIINF:NISUP,NJINF:NJSUP)
			       ENDIF
			     ENDDO
			     IF(LVTM .OR. LVTT)THEN
			       ZWORK3D=ZWORK3V
			     ENDIF
                             DEALLOCATE(ZWORK3V,ZTEM1,ZTEMV)
			   ENDIF
!!!!! Avril 99 Ajout ULM et VTM en CH
                        ELSE
        	          ZWORK3D=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
		                       NJINF-NJL+1:NJSUP-NJL+1, &
        	                     :,NTIMEDIA(JLOOPT,KLOOP,1),JLOOPN, &
				     NPROCDIA(JLOOPP,KLOOP))
                        ENDIF
!print *,' OPER NIINF-NIL+1:NISUP-NIL+1  ',NIINF-NIL+1,NISUP-NIL+1
!print *,' OPER NJINF-NJL+1:NJSUP-NJL+1 ',NJINF-NJL+1,NJSUP-NJL+1
!print *,' OPER XVAR ',XVAR(NIINF-NIL+1,NJINF-NJL+1,1,JLOOPT,JLOOPN,JLOOPP)
!print *,' OPER XVAR ',XVAR(NISUP-NIL+1,NJSUP-NJL+1,SIZE(ZWORK3D,3),JLOOPT,JLOOPN,JLOOPP)
!                      WRITE(CLEGEND2(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
                       WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
!!!!!!!!!!!!!!!!!!!!!!!!!    CH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        	        IF(LCH)THEN
                if(nverbia >0)then
        	print *,' **oper** AP LCH ',LCH
                endif

          		  IF(NBLVLKDIA(KLOOP,1) == 0)THEN

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  altitudes Z (Formulation sequentielle)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

			    IF(.NOT.LZINCRDIA(KLOOP))THEN
          		    DO JLOOPZ=1,NBLVLZDIA(KLOOP)

			      IZ=XLVLZDIA(JLOOPZ,KLOOP)
! Pour LANIMK
			      XLOOPZ=XLVLZDIA(JLOOPZ,KLOOP)
	                      if(nverbia > 0)then
			      print *,' ***oper XLOOPZ ',XLOOPZ
			      endif
!                             XLOOPZ=IZ
			      IF(LANIMK)THEN
            		        IF(JLOOPZ == 1)THEN
		              CALL FMFREE(YBID,YBID,IRESP)
                              if(nverbia >0)then
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
                              endif
            		          CALL FMATTR(YBID,YBID,IBID,IRESP)
            		          CALL GOPWK(9,IBID,3)
            		          ISEGM=ISEGM+1
            		          ISEGD=ISEGM
            		          CALL GFLAS1(ISEGM)
                                ELSE
            		          ISEGM=ISEGM+1
            		          CALL GFLAS1(ISEGM)
                                ENDIF
			      ENDIF
! Pour LANIMK
			      IF(LPXT .OR. LPYT)THEN
				IF(JLOOPT == 1)THEN
				  IF(ALLOCATED(ZSTAB))THEN
				    DEALLOCATE(ZSTAB)
				  ENDIF
				  IX=NISUP-NIINF+1
				  IY=NJSUP-NJINF+1
				  ILENW=NBTIMEDIA(KLOOP,1)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			     IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			     ALLOCATE(XPRDAT(16,ILENW))
			   ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
				  IF(IX /= 1 .AND. IY /= 1)THEN
				    IF(LPXT)THEN
				    print *,' _PXT_ --> Profil horizontal // X f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP
				    ELSE IF(LPYT)THEN
				    print *,' _PYT_ --> Profil horizontal // Y f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP

                                    ENDIF

                                    LPBREAD=.TRUE.
				    RETURN
				  ELSE IF(IY == 1 .AND. IX /= 1)THEN
				    ALLOCATE(ZTEM2D(IX,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				  ELSE IF(IX == 1 .AND. IY /= 1)THEN
		                    ALLOCATE(ZTEM2D(IY,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
                                  ENDIF
				  ALLOCATE(ZSTAB(IX,IY))
				ENDIF
				CALL INTERP_FORDIACHRO(IZ,NKL,NKH,ZWORK3D,ZSTAB)
				ZWORKT(JLOOPT)=XTRAJT(NLOOPT,1)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			     CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
			   ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				IF(LPXT)THEN
				  ZTEM2D(:,JLOOPT)=ZSTAB(:,1)
				ELSE IF(LPYT)THEN
				  ZTEM2D(:,JLOOPT)=ZSTAB(1,:)
				ENDIF
				IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
	                          ILENU=LEN_TRIM(CUNITGAL)
				  ILENT=LEN(CUNITGAL)
				  IF(ILENT-ILENU-2+1 < 8)THEN
	                            IF(NVERBIA > 0)THEN
				    print *,' CUNITGAL ILENT-ILENU-2+1 < 8 ',CUNITGAL 
				    ENDIF
				  ELSE
				  IF(LEV)THEN
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A2,''='',I5)')'PV',IZ
				  ELSE IF(LSV3)THEN
				    IF(LXYZ00)THEN
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')CGROUPSV3(1:3),IZ
!			              WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'Z00',IZ
				    ELSE
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'SV3',IZ
				    ENDIF
				  ELSE
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A1,''='',I5)')CTYPHOR,IZ
				  ENDIF
				  ENDIF
				  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
		             DEALLOCATE(XPRDAT)
			   ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				  IF(.NOT.LPBREAD)THEN
				    IF(KLOOP == NSUPERDIA)THEN
              		              CALL NGPICT(1,1)
              		              CALL GQACWK(1,IER,INB,IWK)
              		              IF(INB > 1)CALL NGPICT(2,3)
				    ENDIF
				    DEALLOCATE(ZWORKT,ZTEM2D,ZSTAB)
				  ENDIF
				ENDIF
			      ELSE
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			     IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			     ALLOCATE(XPRDAT(16,1))
			     CALL LOAD_XPRDAT(1,NLOOPT)
			   ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
          		        CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
		             DEALLOCATE(XPRDAT)
			   ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	if(nverbia > 0)then
	  print *,' **oper AP TRACEH1 IZ II,IJ,IK,KLOOP ',IZ,II,IJ,IK,KLOOP
	endif
			      ENDIF
	        IF(LCV .AND. JLOOPZ == NBLVLZDIA(KLOOP))THEN    
		IINFCV=NIINF; IISUPCV=NISUP; IJINFCV=NJINF; IJSUPCV=NJSUP
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
	        IF(NVERBIA > 0)THEN
	        print *,'oper 1 NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
	        ENDIF
	        ENDIF

                              CALL CLOSF(JLOOPT,NBTIMEDIA(KLOOP,1), &
			      ISEGD,ISEGM,KLOOP)
	        IF(LCV .AND. JLOOPZ == NBLVLZDIA(KLOOP))THEN    
		  NIINF=IINFCV; NISUP=IISUPCV; NJINF=IJINFCV; NJSUP=IJSUPCV
	        ENDIF

          		    ENDDO

			    ELSE

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  altitudes Z (Formulation incrementale)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Mars 2000
                              XLOOPZ=XLVLZDIA(1,KLOOP)-XLVLZDIA(3,KLOOP)
!Mars 2000

          		    DO JLOOPZ=INT(XLVLZDIA(1,KLOOP)),INT(XLVLZDIA(2,KLOOP)), &
				      INT(XLVLZDIA(3,KLOOP))
			      IZ=JLOOPZ
! Pour LANIMK
!Mars 2000
                              XLOOPZ=XLOOPZ+XLVLZDIA(3,KLOOP)
	                      if(nverbia > 0)then
			      print *,' ***oper XLOOPZ ',XLOOPZ
			      endif
!                             XLOOPZ=IZ
!Mars 2000
			      IF(LANIMK)THEN
            		        IF(JLOOPZ == INT(XLVLZDIA(1,KLOOP)))THEN
		              CALL FMFREE(YBID,YBID,IRESP)
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
            		          CALL FMATTR(YBID,YBID,IBID,IRESP)
            		          CALL GOPWK(9,IBID,3)
            		          ISEGM=ISEGM+1
            		          ISEGD=ISEGM
            		          CALL GFLAS1(ISEGM)
                                ELSE
            		          ISEGM=ISEGM+1
            		          CALL GFLAS1(ISEGM)
                                ENDIF
			      ENDIF
! Pour LANIMK
  			      IF(LPXT .OR. LPYT)THEN
				IF(JLOOPT == 1)THEN
				  IF(ALLOCATED(ZSTAB))THEN
				    DEALLOCATE(ZSTAB)
				  ENDIF
				  IX=NISUP-NIINF+1
				  IY=NJSUP-NJINF+1
				  ILENW=NBTIMEDIA(KLOOP,1)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			     IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			     ALLOCATE(XPRDAT(16,ILENW))
			   ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
				  IF(IX /= 1 .AND. IY /= 1)THEN
				    IF(LPXT)THEN
				    print *,' _PXT_ --> Profil horizontal // X f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP
				    ELSE IF(LPYT)THEN
				    print *,' _PYT_ --> Profil horizontal // Y f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP

                                    ENDIF

                                    LPBREAD=.TRUE.
				    RETURN
				  ELSE IF(IY == 1 .AND. IX /= 1)THEN
				    ALLOCATE(ZTEM2D(IX,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				  ELSE IF(IX == 1 .AND. IY /= 1)THEN
		                    ALLOCATE(ZTEM2D(IY,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
                                  ENDIF
				  ALLOCATE(ZSTAB(IX,IY))
				ENDIF
				CALL INTERP_FORDIACHRO(IZ,NKL,NKH,ZWORK3D,ZSTAB)
				ZWORKT(JLOOPT)=XTRAJT(NLOOPT,1)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			     CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
			   ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
				IF(LPXT)THEN
				  ZTEM2D(:,JLOOPT)=ZSTAB(:,1)
				ELSE IF(LPYT)THEN
				  ZTEM2D(:,JLOOPT)=ZSTAB(1,:)
				ENDIF
				IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
	                          ILENU=LEN_TRIM(CUNITGAL)
				  ILENT=LEN(CUNITGAL)
				  IF(ILENT-ILENU-2+1 < 8)THEN
	                            IF(NVERBIA > 0)THEN
				    print *,' CUNITGAL ILENT-ILENU-2+1 < 8 ',CUNITGAL 
				    ENDIF
				  ELSE
				  IF(LEV)THEN
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A2,''='',I5)')'PV',IZ
				  ELSE IF(LSV3)THEN
				    IF(LXYZ00)THEN
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')CGROUPSV3(1:3),IZ
!	        		      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'Z00',IZ
				    ELSE
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'SV3',IZ
				    ENDIF
				  ELSE
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A1,''='',I5)')CTYPHOR,IZ
				  ENDIF
				  ENDIF
				  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
		             DEALLOCATE(XPRDAT)
			   ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
				  IF(.NOT.LPBREAD)THEN
				    IF(KLOOP == NSUPERDIA)THEN
              		              CALL NGPICT(1,1)
              		              CALL GQACWK(1,IER,INB,IWK)
              		              IF(INB > 1)CALL NGPICT(2,3)
				    ENDIF
				    DEALLOCATE(ZWORKT,ZTEM2D,ZSTAB)
				  ENDIF
				ENDIF
			      ELSE
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			     IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			     ALLOCATE(XPRDAT(16,1))
			     CALL LOAD_XPRDAT(1,NLOOPT)
			   ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
          		        CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		             DEALLOCATE(XPRDAT)
			   ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
	if(nverbia > 0)then
	  print *,' **oper AP TRACEH2 IZ II,IJ,IK,KLOOP ',IZ,II,IJ,IK,KLOOP
	endif
			      ENDIF
	        IF(LCV .AND. JLOOPZ == NINT(XLVLZDIA(2,KLOOP)))THEN
		IINFCV=NIINF; IISUPCV=NISUP; IJINFCV=NJINF; IJSUPCV=NJSUP
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
	        IF(NVERBIA > 0)THEN
	        print *,'oper 2 NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
		print *,' oper 2 JLOOPZ NINT(XLVLZDIA(2,KLOOP)) ',JLOOPZ,NINT(XLVLZDIA(2,KLOOP))
	        ENDIF
	        ENDIF

                              CALL CLOSF(JLOOPT,NBTIMEDIA(KLOOP,1), &
			      ISEGD,ISEGM,KLOOP)
	        IF(LCV .AND. JLOOPZ == NINT(XLVLZDIA(2,KLOOP)))THEN    
		  NIINF=IINFCV; NISUP=IISUPCV; NJINF=IJINFCV; NJSUP=IJSUPCV
	        ENDIF

          		    ENDDO

			    ENDIF
        
          		  ELSE
        
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les niveaux de modele (Formulation sequentielle)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          		    DO JLOOPK=1,NBLVLKDIA(KLOOP,1)
! Pour LANIMK
			      NLOOPK=JLOOPK
			      IF(LANIMK)THEN
            		        IF(JLOOPK == 1)THEN
		              CALL FMFREE(YBID,YBID,IRESP)
                              if(nverbia >0)then
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
                              endif
            		          CALL FMATTR(YBID,YBID,IBID,IRESP)
            		          CALL GOPWK(9,IBID,3)
            		          ISEGM=ISEGM+1
            		          ISEGD=ISEGM
            		          CALL GFLAS1(ISEGM)
                                ELSE
            		          ISEGM=ISEGM+1
            		          CALL GFLAS1(ISEGM)
                                ENDIF
			      ENDIF
! Pour LANIMK
                              IZ=NLVLKDIA(JLOOPK,KLOOP,1)
			      if(nverbia > 0)then
				print *,' **oper Niveau K transmis a INTERP ',IZ
                                print *,' **oper LPR,LTK,LEV,LSV3 ',LPR,LTK,LEV,LSV3
			      endif
  			      IF(LPXT .OR. LPYT)THEN
				IF(JLOOPT == 1)THEN
				  IF(ALLOCATED(ZSTAB))THEN
				    DEALLOCATE(ZSTAB)
				  ENDIF
				  IX=NISUP-NIINF+1
				  IY=NJSUP-NJINF+1
				  ILENW=NBTIMEDIA(KLOOP,1)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			     IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			     ALLOCATE(XPRDAT(16,ILENW))
			   ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
				  IF(IX /= 1 .AND. IY /= 1)THEN
				    IF(LPXT)THEN
				    print *,' _PXT_ --> Profil horizontal // X f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP
				    ELSE IF(LPYT)THEN
				    print *,' _PYT_ --> Profil horizontal // Y f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP

                                    ENDIF

                                    LPBREAD=.TRUE.
				    RETURN
				  ELSE IF(IY == 1 .AND. IX /= 1)THEN
				    ALLOCATE(ZTEM2D(IX,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				  ELSE IF(IX == 1 .AND. IY /= 1)THEN
		                    ALLOCATE(ZTEM2D(IY,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
                                  ENDIF
				  ALLOCATE(ZSTAB(IX,IY))
				ENDIF
				CALL INTERP_FORDIACHRO(IZ,NKL,NKH,ZWORK3D,ZSTAB)
				ZWORKT(JLOOPT)=XTRAJT(NLOOPT,1)
                                IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			          CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
			        ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
				IF(LPXT)THEN
				  ZTEM2D(:,JLOOPT)=ZSTAB(:,1)
				ELSE IF(LPYT)THEN
				  ZTEM2D(:,JLOOPT)=ZSTAB(1,:)
				ENDIF
				IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
	                          ILENU=LEN_TRIM(CUNITGAL)
				  ILENT=LEN(CUNITGAL)
				  IF(ILENT-ILENU-2+1 < 8)THEN
	                            IF(NVERBIA > 0)THEN
				    print *,' CUNITGAL ILENT-ILENU-2+1 < 8 ',CUNITGAL 
				    ENDIF
				  ELSE
				  IF(LEV)THEN
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A2,''='',I5)')'PV',IZ
				  ELSE IF(LSV3)THEN
				    IF(LXYZ00)THEN
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')CGROUPSV3(1:3),IZ
!  			              WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'Z00',IZ
				    ELSE
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'SV3',IZ
				    ENDIF
				  ELSE
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A1,''='',I5)')CTYPHOR,IZ
				  ENDIF
				  ENDIF
				  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                    DEALLOCATE(XPRDAT)
			          ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
				  IF(.NOT.LPBREAD)THEN
				    IF(KLOOP == NSUPERDIA)THEN
              		              CALL NGPICT(1,1)
              		              CALL GQACWK(1,IER,INB,IWK)
              		              IF(INB > 1)CALL NGPICT(2,3)
				    ENDIF
				    DEALLOCATE(ZWORKT,ZTEM2D,ZSTAB)
				  ENDIF
				ENDIF
			      ELSE
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			     IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			     ALLOCATE(XPRDAT(16,1))
			     CALL LOAD_XPRDAT(1,NLOOPT)
			   ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
          		        CALL TRACEH_FORDIACHRO(NLVLKDIA(JLOOPK, &
						     KLOOP,1),ZWORK3D,KLOOP)
                           IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		             DEALLOCATE(XPRDAT)
			   ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
          if(nverbia > 0)then
	  print *,' **oper AP TRACEH3 IZ II,IJ,IK,KLOOP ',NLVLKDIA(JLOOPK,KLOOP,1),II,IJ,IK,KLOOP
	  endif
                              ENDIF
	        IF(LCV .AND. JLOOPK == NBLVLKDIA(KLOOP,1))THEN
		IINFCV=NIINF; IISUPCV=NISUP; IJINFCV=NJINF; IJSUPCV=NJSUP
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
	        IF(NVERBIA > 0)THEN
	        print *,' oper 3 NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
	        ENDIF
	        ENDIF


                              CALL CLOSF(JLOOPT,NBTIMEDIA(KLOOP,1), &
			      ISEGD,ISEGM,KLOOP)
	        IF(LCV .AND. JLOOPK == NBLVLKDIA(KLOOP,1))THEN    
		  NIINF=IINFCV; NISUP=IISUPCV; NJINF=IJINFCV; NJSUP=IJSUPCV
	        ENDIF

          		    ENDDO

          		  ENDIF
!                         CALL ALLOC_FORDIACHRO(1,1,1,1,1,1,3)


!!!!!!!!!!!!!!!!!!!!!!!!!    CV    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			ELSE IF(LCV)THEN

                          IF(.NOT.LDEFCV2CC)THEN  !%%%%%%%%%%%%%%%%%%%%%%%%

			  IF(NLMAX <= 1 .OR. (NLANGLE<0 .OR. NLANGLE>360) .OR. &
			  (NIDEBCOU <=0 .AND. XIDEBCOU == -999.) .OR. &
			  (NJDEBCOU <=0 .AND. XJDEBCOU == -999.))THEN
			    PRINT *,' DEFINISSEZ D''ABORD NIDEBCOU, NJDEBCOU,',&
&                           ' NLMAX, NLANGLE (Pour CV + PV), PROFILE (Pour PV)'
                            PRINT *,'                  ou XIDEBCOU, XJDEBCOU'
			    PRINT *,' PUIS RENTREZ A NOUVEAU VOTRE DIRECTIVE '
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
			    PRINT *,' VALEURS ACTUELLES: '
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                           I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLMAX,NLANGLE,NPROFILE
                            IF(II == 1 .AND. .NOT.LICP .AND. IJ>1 .AND. IK>1)THEN
                            PRINT *,'DANS LE CAS CONSIDERE (CV // Y), si vous voulez ',&
                            'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                          & '' NLMAX='',I6,'' NLANGLE= 90 '')',NIL,NJL,NJH-NJL+1
                            ENDIF
                            IF(IJ == 1 .AND. .NOT.LJCP .AND. II > 1 .AND. IK >1)THEN
                            PRINT *,' DANS LE CAS CONSIDERE (CV // X), si vous voulez ',&
                            &'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                            &'' NLMAX='',I6,'' NLANGLE= 0 '')',NIL,NJl,NIH-NIL+1
                            ENDIF
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
			  ELSE
			    IF((.NOT.LPVT .AND. .NOT.LPVKT .AND. .NOT.LPVKT1) .OR. &
			      (LPVT .AND. JLOOPT==1) .OR.  &
			      (LPVKT .AND. JLOOPT==1) .OR.  &
			      (LPVKT1 .AND. JLOOPT==1))THEN  !!!!
                            IF(II == 1 .AND. .NOT.LICP .AND. IJ>1 .AND. IK>1)THEN
                            PRINT *,'DANS LE CAS CONSIDERE (CV // Y), si vous voulez ',&
                            'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                          & '' NLMAX='',I6,'' NLANGLE= 90 '')',NIL,NJL,NJH-NJL+1
                            ENDIF
                            IF(IJ == 1 .AND. .NOT.LJCP .AND. II > 1 .AND. IK >1)THEN
                            PRINT *,' DANS LE CAS CONSIDERE (CV // X), si vous voulez ',&
                            &'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                            &'' NLMAX='',I6,'' NLANGLE= 0 '')',NIL,NJl,NIH-NIL+1
                            ENDIF
			    PRINT *,' VALEURS DES PARAMETRES DE DEFINITION DE LA COUPE',&
&                           ' ou DU PROFIL :'
			    IF(XIDEBCOU == -999. .AND. XJDEBCOU == -999.)THEN
			      PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                             NJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
                            ELSE
			      PRINT '('' XIDEBCOU:'',F7.1,'' XJDEBCOU:'',F7.1,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',XIDEBCOU, &
&                             XJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    ENDIF     !!!!
			    ENDIF
                          ENDIF
			  IF((LPV.OR.LPVT.OR.LPVKT.OR.LPVKT1) .AND. NPROFILE > NLMAX)THEN
			    PRINT *,' PROFILE DOIT ETRE <= NLMAX '
			    print *,' NLMAX:',NLMAX,' PROFILE: ',NPROFILE
			    print *,' Valeur des autres informations utiles :'
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5, &
&                           '' NLANGLE:'',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLANGLE
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
			  ENDIF
			  IF((LPV.OR.LPVT.OR.LPVKT.OR.LPVKT1) .AND. NPROFILE <= 0)THEN
			    PRINT *,' PROFILE DOIT ETRE DEFINI.',&
			    &'Sa valeur actuelle: ',NPROFILE
			    print *,' Valeur des autres informations utiles :'
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                           I6,'' NLANGLE:'',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLMAX,NLANGLE
                            print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
			  ENDIF

			  ENDIF                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			  CALL VERIFLEN_FORDIACHRO
			  CALL MEMCV
		          ALLOCATE (ZTEMCV(NLMAX,1:IKU))
			  CALL PRECOU_FORDIACHRO(ZWORK3D,ZTEMCV)
			  IF(LPV)THEN
			    L1DT=.FALSE.
! Janvier 2001
			    IF(LUMVM.OR.LUTVT.OR.LSUMVM.OR.LSUTVT.OR.&
			       LDIRWIND)THEN
			    ILENT=LEN_TRIM(CTITGAL)
			    ILENU=LEN_TRIM(CUNITGAL)
			    YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
			    YTEXTE(ILENT+1:ILENT+1)=' '
			    YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
                            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			      ALLOCATE(XPRDAT(16,1))
			      CALL LOAD_XPRDAT(1,NLOOPT)
			    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                            CALL TRACEV_FORDIACHRO(ZTEMCV,KLOOP,YTEXTE(1: &
                                                         LEN_TRIM(YTEXTE)))
                            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
		              DEALLOCATE(XPRDAT)
			    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
			    ELSE
! Janvier 2001
			    ALLOCATE(ZTEM1D(IKU),ZWORKZ(IKU))
! Modif AOUT 97
			    ZTEM1D(:)=XSPVAL; ZWORKZ(:)=0.
!                           ZTEM1D(:)=0.; ZWORKZ(:)=0.
			    ZTEM1D(MAX(IKB,NKL):MIN(IKE,NKH))= &
			    ZTEMCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
			    ZWORKZ(:)=XWORKZ(NPROFILE,:,NMGRID)
                            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			      ALLOCATE(XPRDAT(16,1))
			      CALL LOAD_XPRDAT(1,NLOOPT)
			    ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
			    CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
                            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
		               DEALLOCATE(XPRDAT)
			    ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
			    ENDIF
                          ELSE IF(LPVT .OR. LPVKT.OR. LPVKT1)THEN
			    L1DT=.FALSE.
      		            IF(JLOOPT == 1)THEN
      		              ILENW=NBTIMEDIA(KLOOP,1)
                              IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			        IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			        ALLOCATE(XPRDAT(16,ILENW))
			      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
      		              ALLOCATE(ZTEM2D(1:IKU,ILENW))
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
!Fev 2002
                              IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT )THEN
!                             IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT &
!		      .OR.LDIRWIND)THEN
!Fev 2002
			      IF(ALLOCATED(XTEM2D))DEALLOCATE(XTEM2D)
			      IF(ALLOCATED(XTEM2D2))DEALLOCATE(XTEM2D2)
      		              ALLOCATE(XTEM2D(1:IKU,ILENW))
      		              ALLOCATE(XTEM2D2(1:IKU,ILENW))
			      XTEM2D=XSPVAL
			      XTEM2D2=XSPVAL
			      ENDIF
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
      		              ALLOCATE(ZWORKT(ILENW))
			      ALLOCATE(ZWORKZ2(IKU))
			      ZWORKZ2(:)=0.; ZWORKT(:)=0.; ZTEM2D(:,:)=0.
		              ZTEM2D=XSPVAL
			      ZWORKZ2(:)=XWORKZ(NPROFILE,:,NMGRID)
      		            ENDIF
                            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
			      CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
			    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
      		            ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
      		            ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),JLOOPT)= &
      		            ZTEMCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
!Fev 2002
                            IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT )THEN
!                           IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT &
!		      .OR.LDIRWIND)THEN
!Fev 2002
      		              XTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),JLOOPT)= &
      		              ZTEMCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
      		              XTEM2D2(MAX(IKB,NKL):MIN(IKE,NKH),JLOOPT)= &
      		              XWCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
                            ENDIF
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
      		            IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
      		              XPVMIN=MINVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
      		              XPVMAX=MAXVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
      		              CALL VALMNMX(XPVMIN,XPVMAX)
                              IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
      			        XPVMIN=XPVMIN-1.
      			        XPVMAX=XPVMAX+1.
                              ENDIF
      		              CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                              IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
		                DEALLOCATE(XPRDAT)
			      ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
      		              DEALLOCATE(ZTEM2D,ZWORKT,ZWORKZ2)
			      IF(ALLOCATED(XTEM2D))DEALLOCATE(XTEM2D)
			      IF(ALLOCATED(XTEM2D2))DEALLOCATE(XTEM2D2)
      		            ENDIF
			  ELSE
			    ILENT=LEN_TRIM(CTITGAL)
			    ILENU=LEN_TRIM(CUNITGAL)
			    YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
			    YTEXTE(ILENT+1:ILENT+1)=' '
			    YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
                            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
			      ALLOCATE(XPRDAT(16,1))
			      CALL LOAD_XPRDAT(1,NLOOPT)
			    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                            CALL TRACEV_FORDIACHRO(ZTEMCV,KLOOP,YTEXTE(1: &
                            LEN_TRIM(YTEXTE)))
                            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		              DEALLOCATE(XPRDAT)
			    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
			  ENDIF
			  IF((LCV .OR. LPV) .AND. .NOT. LPVT .AND. .NOT. LPVKT .AND. .NOT.LPVKT1)THEN
!!Fev 2002
                IF(JLOOPT ==  NBTIMEDIA(KLOOP,1))THEN
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
                ENDIF
!!Fev 2002
                              CALL CLOSF(JLOOPT,NBTIMEDIA(KLOOP,1), &
			      ISEGD,ISEGM,KLOOP)

			  ENDIF

			  DEALLOCATE(ZTEMCV)
			  DEALLOCATE(XWORKZ,XWZ)
			  IF(ALLOCATED(ZTEM1D))THEN
			    DEALLOCATE(ZTEM1D)
                          ENDIF
			  IF(ALLOCATED(ZWORKZ))THEN
			    DEALLOCATE(ZWORKZ)
                          ENDIF

        	        ENDIF
        	      ENDDO
			  IF((LPVT.AND..NOT.LPBREAD) .OR. LPVKT .OR. LPVKT1)THEN
!                           IF(KLOOP == NSUPERDIA)CALL FRAME
              		    IF(KLOOP == NSUPERDIA)THEN
              		      CALL NGPICT(1,1)
              		      CALL GQACWK(1,IER,INB,IWK)
              		      IF(INB > 1)CALL NGPICT(2,3)
              		    ENDIF
			  ENDIF
        
        	    ELSE
        
!       	      print *,' NBTIMEDIA(KLOOP,1) ',NBTIMEDIA(KLOOP,1)
!       	      print *,' NTIMEDIA(1 et 2,KLOOP,1) ',NTIMEDIA(1,KLOOP,1), &
!                     NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Boucle sur les  temps (Formulation incrementale)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        	      DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1), &
						  NTIMEDIA(3,KLOOP,1)
                        NLOOPT=JLOOPT
                if(nverbia >0)then
        	print *,' **oper**B JLOOPT ',JLOOPT
                endif

			IF(LANIMT)THEN
			  IF(LPVT .OR. LPVKT .OR. LPVKT1)THEN
			    print *,' ANIMATION IMPOSSIBLE avec _PVT_ ou _PVKT_ ou _PVKT1_'
			    print *,' LANIMT remis a .FALSE. '
			    LANIMT=.FALSE.
			  ELSE IF(LPV .AND. NSUPERDIA>1)THEN
			    print *,' ANIMATION IMPOSSIBLE ', &
			    &'avec _PV_ et superpositions'
			    print *,' LANIMT remis a .FALSE. '
			    print *,' mais POSSIBLE sous la forme : ',&
&                           'GPE_PV__P_1 ou GPE_PV__P_1_T_300_TO_3600 '
			    print *,' PENSER a fournir les bornes dans ',&
&        		    'XPVMIN_proc= et XPVMAX_proc= et a les activer ',& 
&                	    'avec LMNMXUSER=T '
			    print *,' Rappel : proc=nom du processus tel ',&
&                           'qu''il est enregistre '
			    LANIMT=.FALSE.
			  ELSE
			    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		              CALL FMFREE(YBID,YBID,IRESP)
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
			      CALL FMATTR(YBID,YBID,IBID,IRESP)
			      CALL GOPWK(9,IBID,3)
			      ISEGM=ISEGM+1
			      ISEGD=ISEGM
			      CALL GFLAS1(ISEGM)
		              ITIMEND=NTIMEDIA(1,KLOOP,1) + &
		              (((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/  &
		              NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
	                      if(nverbia > 0)then
			      print *,' ITIMEND ',ITIMEND
			      endif
			    ELSE
			      ISEGM=ISEGM+1
			      CALL GFLAS1(ISEGM)
			    ENDIF
			  ENDIF
			ENDIF

		        CALL RESOLV_TIMES(JLOOPT)
                        if(nverbia > 0)then
			  print *,' **oper LULM LCH LMUMVM LDIRWM LDIRWIND lig 3088 ',LULM,LCH,LMUMVM,LDIRWM,LDIRWIND
			endif
		        IF(LUMVM .OR. LUTVT .OR. LMUMVM .OR. LMUTVT .OR. &
		           LULMWM .OR. LULTWT .OR. LULM .OR. LULT .OR.   &
		           LVTM .OR. LVTT .OR. LSUMVM .OR. LSUTVT .OR.   &
			   LDIRWM .OR. LDIRWT .OR. &
			   LMLSUMVM .OR. LMLSUTVT)THEN
                        if(nverbia > 0)then
			  print *,' **oper ds test LULM LCH LMUMVM LDIRWM LDIRWIND lig 3096 ',LULM,LCH,LMUMVM,LDIRWM,LDIRWIND
			endif
        	          ZWORK3D=XU(NIINF-NIL+1:NISUP-NIL+1, &
		                       NJINF-NJL+1:NJSUP-NJL+1, &
        	                     :,JLOOPT,JLOOPN, &
				     NPROCDIA(JLOOPP,KLOOP))
!!!!! Avril 99 Ajout ULM et VTM en CH
                           IF((LCH.AND.LULM).OR.(LCH.AND.LVTM).OR. &
			   (LCH.AND.LULT).OR.(LCH.AND.LVTT))THEN
			     ALLOCATE(ZWORK3V(SIZE(ZWORK3D,1), &
			     SIZE(ZWORK3D,2),SIZE(ZWORK3D,3)))
			     ALLOCATE(ZTEM1(IIU,IJU),ZTEMV(IIU,IJU))
			     ZTEM1=0.
			     ZTEMV=0.
        	             ZWORK3V=XVAR(NIINF-NIL+1:NISUP-NIL+1, &
		                       NJINF-NJL+1:NJSUP-NJL+1, &
        	                     :,JLOOPT,JLOOPN, &
				     NPROCDIA(JLOOPP,KLOOP))
                             DO JKLOOP=1,IKU
			       IF(JKLOOP < MAX(IKB,NKL) .OR. &
				  JKLOOP > MIN(IKE,NKH))THEN
			       ELSE
			         ZTEM1(NIINF:NISUP,NJINF:NJSUP)= &
				 ZWORK3D(:,:,JKLOOP-NKL+1)
			         ZTEMV(NIINF:NISUP,NJINF:NJSUP)= &
				 ZWORK3V(:,:,JKLOOP-NKL+1)
				 if(nverbia > 0)then
				   print *,'** oper ZTEM1(NIINF,NJINF),&
                                   &ZTEM1(NISUP,NJSUP) av rota',&
				   ZTEM1(NIINF,NJINF),ZTEM1(NISUP,NJSUP)
				   print *,'** oper JKLOOP NKL ',&
				   JKLOOP,NKL
				 endif
!!!!essai Nov 2001 pour prise en compte PH 29/11/2001 .. A suivre
                                 CALL VERIFLEN_FORDIACHRO
!!!!essai Nov 2001
				 CALL ROTA(ZTEM1,ZTEMV)
				 if(nverbia > 0)then
				   print *,'** oper ZTEM1(NIINF,NJINF),&
                                 &  ZTEM1(NISUP,NJSUP) ap rota',&
				   ZTEM1(NIINF,NJINF),ZTEM1(NISUP,NJSUP)
				   print *,'** oper JKLOOP NKL ',&
				   JKLOOP,NKL
				 endif
				 ZWORK3D(:,:,JKLOOP-NKL+1)=ZTEM1(NIINF:NISUP,NJINF:NJSUP)
				 ZWORK3V(:,:,JKLOOP-NKL+1)=ZTEMV(NIINF:NISUP,NJINF:NJSUP)
			       ENDIF
			     ENDDO
			     IF(LVTM .OR. LVTT)THEN
			       ZWORK3D=ZWORK3V
			     ENDIF
                             DEALLOCATE(ZWORK3V,ZTEM1,ZTEMV)
			   ENDIF
!!!!! Avril 99 Ajout ULM et VTM en CH
                        ELSE
			  if(nverbia > 0)then
			  print *,' **oper lig 3149'
			  endif
        	          ZWORK3D=XVAR(NIINF-NIL+1:NISUP-NIL+1,    &
				     NJINF-NJL+1:NJSUP-NJL+1,    &
        	                     :,JLOOPT,JLOOPN,NPROCDIA(JLOOPP,KLOOP))
                        ENDIF
!                       WRITE(CLEGEND2(8:15),'(F8.0)')XTRAJT(JLOOPT,1)
                        WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)

!!!!!!!!!!!!!!!!!!!!!!!!!    CH    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        if(nverbia > 0)then
			  print *,' **oper  AV LCH lig 3166'
			endif

        	        IF(LCH)THEN

          		  IF(NBLVLKDIA(KLOOP,1) == 0)THEN

			    IF(.NOT.LZINCRDIA(KLOOP))THEN
          		    DO JLOOPZ=1,NBLVLZDIA(KLOOP)

                              IZ=XLVLZDIA(JLOOPZ,KLOOP)
! Pour LANIMK
!Mars 2000
			      XLOOPZ=XLVLZDIA(JLOOPZ,KLOOP)
	                      if(nverbia > 0)then
			      print *,' ***oper XLOOPZ ',XLOOPZ
			      endif
!                             XLOOPZ=IZ
!Mars 2000
			      IF(LANIMK)THEN
            		        IF(JLOOPZ == 1)THEN
		              CALL FMFREE(YBID,YBID,IRESP)
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
            		          CALL FMATTR(YBID,YBID,IBID,IRESP)
            		          CALL GOPWK(9,IBID,3)
            		          ISEGM=ISEGM+1
            		          ISEGD=ISEGM
            		          CALL GFLAS1(ISEGM)
                                ELSE
            		          ISEGM=ISEGM+1
            		          CALL GFLAS1(ISEGM)
                                ENDIF
			      ENDIF
! Pour LANIMK
  			      IF(LPXT .OR. LPYT)THEN
				IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
				  IF(ALLOCATED(ZSTAB))THEN
				    DEALLOCATE(ZSTAB)
				  ENDIF
				  IX=NISUP-NIINF+1
				  IY=NJSUP-NJINF+1
				  ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
				  NTIMEDIA(3,KLOOP,1)+1
	                          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                    IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		                    ALLOCATE(XPRDAT(16,ILENW))
                                  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				  ITIMEND=NTIMEDIA(1,KLOOP,1)+(((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
				  NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
				  IF(IX /= 1 .AND. IY /= 1)THEN
				    IF(LPXT)THEN
				    print *,' _PXT_ --> Profil horizontal // X f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP
				    ELSE IF(LPYT)THEN
				    print *,' _PYT_ --> Profil horizontal // Y f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP

                                    ENDIF

                                    LPBREAD=.TRUE.
				    RETURN
				  ELSE IF(IY == 1 .AND. IX /= 1)THEN
				    ALLOCATE(ZTEM2D(IX,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				    IJLT=0
				  ELSE IF(IX == 1 .AND. IY /= 1)THEN
		                    ALLOCATE(ZTEM2D(IY,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				    IJLT=0
                                  ENDIF
				  ALLOCATE(ZSTAB(IX,IY))
				ENDIF
				CALL INTERP_FORDIACHRO(IZ,NKL,NKH,ZWORK3D,ZSTAB)
				IJLT=IJLT+1
	                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			          CALL LOAD_XPRDAT(IJLT,NLOOPT)
                                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				ZWORKT(IJLT)=XTRAJT(NLOOPT,1)
				IF(LPXT)THEN
				  ZTEM2D(:,IJLT)=ZSTAB(:,1)
				ELSE IF(LPYT)THEN
				  ZTEM2D(:,IJLT)=ZSTAB(1,:)
				ENDIF
				IF(JLOOPT == ITIMEND)THEN
	                          ILENU=LEN_TRIM(CUNITGAL)
				  ILENT=LEN(CUNITGAL)
				  IF(ILENT-ILENU-2+1 < 8)THEN
	                            IF(NVERBIA > 0)THEN
				    print *,' CUNITGAL ILENT-ILENU-2+1 < 8 ',CUNITGAL 
				    ENDIF
				  ELSE
				  IF(LEV)THEN
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A2,''='',I5)')'PV',IZ
				  ELSE IF(LSV3)THEN
				    IF(LXYZ00)THEN
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')CGROUPSV3(1:3),IZ
!			              WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'Z00',IZ
				    ELSE
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'SV3',IZ
				    ENDIF
				  ELSE
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A1,''='',I5)')CTYPHOR,IZ
				  ENDIF
				  ENDIF
				  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
	                          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                    DEALLOCATE(XPRDAT)
                                  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				  IF(.NOT.LPBREAD)THEN
				    IF(KLOOP == NSUPERDIA)THEN
              		              CALL NGPICT(1,1)
              		              CALL GQACWK(1,IER,INB,IWK)
              		              IF(INB > 1)CALL NGPICT(2,3)
				    ENDIF
				    DEALLOCATE(ZWORKT,ZTEM2D,ZSTAB)
				  ENDIF
				ENDIF
			      ELSE
	if(nverbia > 0)then
	  print *,' **oper AP TRACEH4 IZ II,IJ,IK ',IZ,II,IJ,IK
	endif
	                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		                  ALLOCATE(XPRDAT(16,1))
		                  CALL LOAD_XPRDAT(1,NLOOPT)
                                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
          		        CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)
	                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                  DEALLOCATE(XPRDAT)
                                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	if(nverbia > 0)then
	  print *,' **oper AP TRACEH4 IZ II,IJ,IK ',IZ,II,IJ,IK
	endif
			      ENDIF
	        IF(LCV .AND. JLOOPZ == NBLVLZDIA(KLOOP))THEN
		IINFCV=NIINF; IISUPCV=NISUP; IJINFCV=NJINF; IJSUPCV=NJSUP
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
	        IF(NVERBIA > 0)THEN
	        print *,' oper 4 NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
	        ENDIF
	        ENDIF

                              CALL CLOSF(JLOOPT,ITIMEND, &
			      ISEGD,ISEGM,KLOOP)
	        IF(LCV .AND. JLOOPZ == NBLVLZDIA(KLOOP))THEN    
		  NIINF=IINFCV; NISUP=IISUPCV; NJINF=IJINFCV; NJSUP=IJSUPCV
	        ENDIF

          		    ENDDO

			    ELSE

!Mars 2000
                            XLOOPZ=XLVLZDIA(1,KLOOP)-XLVLZDIA(3,KLOOP)
!Mars 2000
			    DO JLOOPZ=INT(XLVLZDIA(1,KLOOP)),INT(XLVLZDIA(2,KLOOP)), &
                                      INT(XLVLZDIA(3,KLOOP))
                              IZ=JLOOPZ
! Pour LANIMK
!Mars 2000
			      XLOOPZ=XLOOPZ+XLVLZDIA(3,KLOOP)
	                      if(nverbia > 0)then
			      print *,' ***oper XLOOPZ ',XLOOPZ
			      endif
!                             XLOOPZ=IZ
!Mars 2000
			      IF(LANIMK)THEN
            		        IF(JLOOPZ == XLVLZDIA(1,KLOOP))THEN
		              CALL FMFREE(YBID,YBID,IRESP)
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
            		          CALL FMATTR(YBID,YBID,IBID,IRESP)
            		          CALL GOPWK(9,IBID,3)
            		          ISEGM=ISEGM+1
            		          ISEGD=ISEGM
            		          CALL GFLAS1(ISEGM)
                                ELSE
            		          ISEGM=ISEGM+1
            		          CALL GFLAS1(ISEGM)
                                ENDIF
			      ENDIF
! Pour LANIMK
  			      IF(LPXT .OR. LPYT)THEN
				IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
				  IF(ALLOCATED(ZSTAB))THEN
				    DEALLOCATE(ZSTAB)
				  ENDIF
				  IX=NISUP-NIINF+1
				  IY=NJSUP-NJINF+1
				  ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
				  NTIMEDIA(3,KLOOP,1)+1
	                          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                    IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		                    ALLOCATE(XPRDAT(16,ILENW))
                                  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				  ITIMEND=NTIMEDIA(1,KLOOP,1)+(((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
				  NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
				  IF(IX /= 1 .AND. IY /= 1)THEN
				    IF(LPXT)THEN
				    print *,' _PXT_ --> Profil horizontal // X f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP
				    ELSE IF(LPYT)THEN
				    print *,' _PYT_ --> Profil horizontal // Y f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP

                                    ENDIF

                                    LPBREAD=.TRUE.
				    RETURN
				  ELSE IF(IY == 1 .AND. IX /= 1)THEN
				    ALLOCATE(ZTEM2D(IX,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				    IJLT=0
				  ELSE IF(IX == 1 .AND. IY /= 1)THEN
		                    ALLOCATE(ZTEM2D(IY,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				    IJLT=0
                                  ENDIF
				  ALLOCATE(ZSTAB(IX,IY))
				ENDIF
				CALL INTERP_FORDIACHRO(IZ,NKL,NKH,ZWORK3D,ZSTAB)
				IJLT=IJLT+1
	                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			          CALL LOAD_XPRDAT(IJLT,NLOOPT)
                                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				ZWORKT(IJLT)=XTRAJT(NLOOPT,1)
				IF(LPXT)THEN
				  ZTEM2D(:,IJLT)=ZSTAB(:,1)
				ELSE IF(LPYT)THEN
				  ZTEM2D(:,IJLT)=ZSTAB(1,:)
				ENDIF
				IF(JLOOPT == ITIMEND)THEN
	                          ILENU=LEN_TRIM(CUNITGAL)
				  ILENT=LEN(CUNITGAL)
				  IF(ILENT-ILENU-2+1 < 8)THEN
	                            IF(NVERBIA > 0)THEN
				    print *,' CUNITGAL ILENT-ILENU-2+1 < 8 ',CUNITGAL 
				    ENDIF
				  ELSE
				  IF(LEV)THEN
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A2,''='',I5)')'PV',IZ
				  ELSE IF(LSV3)THEN
				    IF(LXYZ00)THEN
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')CGROUPSV3(1:3),IZ
! 			              WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'Z00',IZ
				    ELSE
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'SV3',IZ
				    ENDIF
				  ELSE
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A1,''='',I5)')CTYPHOR,IZ
				  ENDIF
				  ENDIF
				  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
	                          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                    DEALLOCATE(XPRDAT)
                                  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				  IF(.NOT.LPBREAD)THEN
				    IF(KLOOP == NSUPERDIA)THEN
              		              CALL NGPICT(1,1)
              		              CALL GQACWK(1,IER,INB,IWK)
              		              IF(INB > 1)CALL NGPICT(2,3)
				    ENDIF
				    DEALLOCATE(ZWORKT,ZTEM2D,ZSTAB)
				  ENDIF
				ENDIF

			      ELSE
	                        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		                  ALLOCATE(XPRDAT(16,1))
		                  CALL LOAD_XPRDAT(1,NLOOPT)
                                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

                                CALL TRACEH_FORDIACHRO(IZ,ZWORK3D,KLOOP)
	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  DEALLOCATE(XPRDAT)
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	if(nverbia > 0)then
	  print *,' **oper AP TRACEH5 IZ II,IJ,IK,KLOOP ',IZ,II,IJ,IK,KLOOP
	endif
			      ENDIF
	        IF(LCV .AND. JLOOPZ == NINT(XLVLZDIA(2,KLOOP)))THEN
		IINFCV=NIINF; IISUPCV=NISUP; IJINFCV=NJINF; IJSUPCV=NJSUP
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
	        IF(NVERBIA > 0)THEN
	        print *,' oper 5 NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
		print *,' oper 5 JLOOPZ NINT(XLVLZDIA(2,KLOOP)) ',JLOOPZ,NINT(XLVLZDIA(2,KLOOP))
	        ENDIF
	        ENDIF

                              CALL CLOSF(JLOOPT,ITIMEND, &
			      ISEGD,ISEGM,KLOOP)
	        IF(LCV .AND. JLOOPZ == NINT(XLVLZDIA(2,KLOOP)))THEN    
		  NIINF=IINFCV; NISUP=IISUPCV; NJINF=IJINFCV; NJSUP=IJSUPCV
	        ENDIF

                            ENDDO
			    ENDIF
        
          		  ELSE
        
          		    DO JLOOPK=1,NBLVLKDIA(KLOOP,1)
! Pour LANIMK
			      NLOOPK=JLOOPK
			      IF(LANIMK)THEN
            		        IF(JLOOPK == 1)THEN
		              CALL FMFREE(YBID,YBID,IRESP)
		              print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
            		          CALL FMATTR(YBID,YBID,IBID,IRESP)
            		          CALL GOPWK(9,IBID,3)
            		          ISEGM=ISEGM+1
            		          ISEGD=ISEGM
            		          CALL GFLAS1(ISEGM)
                                ELSE
            		          ISEGM=ISEGM+1
            		          CALL GFLAS1(ISEGM)
                                ENDIF
			      ENDIF
! Pour LANIMK

                              IZ=NLVLKDIA(JLOOPK,KLOOP,1)
			      if(nverbia > 0)then
				print *,' **oper Niveau K transmis a INTERP ',IZ
                                print *,' **oper LPR,LTK,LEV,LSV3 ',LPR,LTK,LEV,LSV3
			      endif

  			      IF(LPXT .OR. LPYT)THEN
				IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
				  IF(ALLOCATED(ZSTAB))THEN
				    DEALLOCATE(ZSTAB)
				  ENDIF
				  IX=NISUP-NIINF+1
				  IY=NJSUP-NJINF+1
				  ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
				  NTIMEDIA(3,KLOOP,1)+1
				  ITIMEND=NTIMEDIA(1,KLOOP,1)+(((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/ &
				  NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		  ALLOCATE(XPRDAT(16,ILENW))
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				  IF(IX /= 1 .AND. IY /= 1)THEN
				    IF(LPXT)THEN
				    print *,' _PXT_ --> Profil horizontal // X f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP
				    ELSE IF(LPYT)THEN
				    print *,' _PYT_ --> Profil horizontal // Y f(t) demande'
				    print *,' Definissez correctement NIINF,NISUP,NJINF,NJSUP. Valeurs actuelles :'
				    print *,' NIINF=',NIINF,' NISUP=',NISUP,' NJINF=',NJINF,' NJSUP=',NJSUP

                                    ENDIF

                                    LPBREAD=.TRUE.
				    RETURN
				  ELSE IF(IY == 1 .AND. IX /= 1)THEN
				    ALLOCATE(ZTEM2D(IX,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				    IJLT=0
				  ELSE IF(IX == 1 .AND. IY /= 1)THEN
		                    ALLOCATE(ZTEM2D(IY,ILENW))
		                    ALLOCATE(ZWORKT(ILENW))
		                    ZTEM2D=XSPVAL
				    IJLT=0
                                  ENDIF
				  ALLOCATE(ZSTAB(IX,IY))
				ENDIF
				CALL INTERP_FORDIACHRO(IZ,NKL,NKH,ZWORK3D,ZSTAB)
				IJLT=IJLT+1
	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			     CALL LOAD_XPRDAT(IJLT,NLOOPT)
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				ZWORKT(IJLT)=XTRAJT(NLOOPT,1)
				IF(LPXT)THEN
				  ZTEM2D(:,IJLT)=ZSTAB(:,1)
				ELSE IF(LPYT)THEN
				  ZTEM2D(:,IJLT)=ZSTAB(1,:)
				ENDIF
				IF(JLOOPT == ITIMEND)THEN
	                          ILENU=LEN_TRIM(CUNITGAL)
				  ILENT=LEN(CUNITGAL)
				  IF(ILENT-ILENU-2+1 < 8)THEN
	                            IF(NVERBIA > 0)THEN
				    print *,' CUNITGAL ILENT-ILENU-2+1 < 8 ',CUNITGAL 
				    ENDIF
				  ELSE
				  IF(LEV)THEN
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A2,''='',I5)')'PV',IZ
				  ELSE IF(LSV3)THEN
				    IF(LXYZ00)THEN
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')CGROUPSV3(1:3),IZ
!			              WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'Z00',IZ
				    ELSE
				      WRITE(CUNITGAL(ILENU+2:ILENT),'(A3,''='',I5)')'SV3',IZ
				    ENDIF
				  ELSE
				    WRITE(CUNITGAL(ILENU+2:ILENT),'(A1,''='',I5)')CTYPHOR,IZ
				  ENDIF
				  ENDIF
				  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  DEALLOCATE(XPRDAT)
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
				  IF(.NOT.LPBREAD)THEN
				    IF(KLOOP == NSUPERDIA)THEN
              		              CALL NGPICT(1,1)
              		              CALL GQACWK(1,IER,INB,IWK)
              		              IF(INB > 1)CALL NGPICT(2,3)
				    ENDIF
				    DEALLOCATE(ZWORKT,ZTEM2D,ZSTAB)
				  ENDIF
				ENDIF

			      ELSE

	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		  ALLOCATE(XPRDAT(16,1))
		CALL LOAD_XPRDAT(1,NLOOPT)
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                CALL TRACEH_FORDIACHRO(NLVLKDIA(JLOOPK,KLOOP,1), &
						     ZWORK3D,KLOOP)
	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  DEALLOCATE(XPRDAT)
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
          if(nverbia > 0)then
	  print *,' **oper AP TRACEH6 IZ II,IJ,IK,KLOOP ',NLVLKDIA(JLOOPK,KLOOP,1),II,IJ,IK,KLOOP
	  endif
                              ENDIF
	        IF(LCV .AND. JLOOPK == NBLVLKDIA(KLOOP,1))THEN
		IINFCV=NIINF; IISUPCV=NISUP; IJINFCV=NJINF; IJSUPCV=NJSUP
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
	        IF(NVERBIA > 0)THEN
	        print *,'oper 6 NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
	        ENDIF
	        ENDIF

                              CALL CLOSF(JLOOPT,ITIMEND, &
			      ISEGD,ISEGM,KLOOP)
	        IF(LCV .AND. JLOOPK == NBLVLKDIA(KLOOP,1))THEN    
		  NIINF=IINFCV; NISUP=IISUPCV; NJINF=IJINFCV; NJSUP=IJSUPCV
	        ENDIF

          		    ENDDO

          		  ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!    CV    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			ELSE IF(LCV)THEN
                        if(nverbia > 0)then
			  print *,' **oper  AP LCV lig 3570'
			endif

			  IF(.NOT.LDEFCV2CC)THEN     !%%%%%%%%%%%%%%%%%%%%%%%%

			  IF(NLMAX <= 1 .OR. (NLANGLE<0 .OR. NLANGLE>360) .OR. &
			  (NIDEBCOU <=0 .AND. XIDEBCOU == -999.) .OR. &
			  (NJDEBCOU <=0 .AND. XJDEBCOU == -999.))THEN
			    PRINT *,' DEFINISSEZ D''ABORD NIDEBCOU, NJDEBCOU,',&
&                           ' NLMAX, NLANGLE (Pour CV + PV), PROFILE (Pour PV)'
                            PRINT *,'                  ou XIDEBCOU, XJDEBCOU'
			    PRINT *,' PUIS RENTREZ A NOUVEAU VOTRE DIRECTIVE '
			    print *,' (Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T)'
			    PRINT *,' VALEURS ACTUELLES: '
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                           I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLMAX,NLANGLE,NPROFILE
                            IF(II == 1 .AND. .NOT.LICP .AND. IJ > 1 .AND. IK >1)THEN
                            PRINT *,' DANS LE CAS CONSIDERE (CV // Y), si vous voulez ',&
                            &'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                            &'' NLMAX='',I6,'' NLANGLE= 90 '')',NIL,NJl,NJH-NJL+1
                            ENDIF
                            IF(IJ == 1 .AND. .NOT.LJCP .AND. II > 1 .AND. IK >1)THEN
                            PRINT *,' DANS LE CAS CONSIDERE (CV // X), si vous voulez ',&
                            &'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                            &'' NLMAX='',I6,'' NLANGLE= 0 '')',NIL,NJl,NIH-NIL+1
                            ENDIF
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
                          ELSE
			    IF((.NOT.LPVT .AND. .NOT.LPVKT .AND. .NOT.LPVKT1) .OR. &
			      (LPVT .AND. JLOOPT==NTIMEDIA(1,KLOOP,1)) .OR.  &
			      (LPVKT .AND. JLOOPT==NTIMEDIA(1,KLOOP,1)) .OR.  &
			      (LPVKT1 .AND. JLOOPT==NTIMEDIA(1,KLOOP,1)))THEN  !!!!
			    IF(II == 1 .AND. .NOT.LICP .AND. IJ > 1 .AND. IK >1)THEN
                            PRINT *,' DANS LE CAS CONSIDERE (CV // Y), si vous voulez ',&
                            &'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                            &'' NLMAX='',I6,'' NLANGLE= 90 '')',NIL,NJl,NJH-NJL+1
                            ENDIF
                            IF(IJ == 1 .AND. .NOT.LJCP .AND. II > 1 .AND. IK >1)THEN
                            PRINT *,' DANS LE CAS CONSIDERE (CV // X), si vous voulez ',&
                            &'la totalite de la coupe, METTEZ: '
                            PRINT '('' NIDEBCOU='',I5,'' NJDEBCOU='',I5,&
                            &'' NLMAX='',I6,'' NLANGLE= 0 '')',NIL,NJl,NIH-NIL+1
                            ENDIF
                            PRINT *,' VALEURS DES PARAMETRES DE DEFINITION DE LA COUPE',&
&                           ' ou DU PROFIL :'
			    IF(XIDEBCOU == -999. .AND. XJDEBCOU == -999.)THEN
			      PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',NIDEBCOU, &
&                             NJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
                            ELSE
			      PRINT '('' XIDEBCOU:'',F7.1,'' XJDEBCOU:'',F7.1,'' NLMAX: '',&
&                             I6,'' NLANGLE:'',I5,'' PROFILE: '',I5)',XIDEBCOU, &
&                             XJDEBCOU,NLMAX,NLANGLE,NPROFILE
			    ENDIF   !!!!
			    ENDIF
                          ENDIF
                        if(nverbia > 0)then
			  print *,' **oper   lig 3613'
			endif
			  IF((LPV.OR.LPVT.OR.LPVKT .OR.LPVKT1) .AND. NPROFILE > NLMAX)THEN
			    PRINT *,' PROFILE DOIT ETRE <= NLMAX '
			    print *,' NLMAX:',NLMAX,' PROFILE: ',NPROFILE
			    print *,' Valeur des autres informations utiles :'
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5, &
&                           '' NLANGLE:'',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLANGLE
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
			  ENDIF
			  IF((LPV.OR.LPVT.OR.LPVKT.OR.LPVKT1) .AND. NPROFILE <= 0)THEN
			    PRINT *,' PROFILE DOIT ETRE DEFINI.',&
			    &'Sa valeur actuelle: ',NPROFILE
			    print *,' Valeur des autres informations utiles :'
			    PRINT '('' NIDEBCOU:'',I5,'' NJDEBCOU:'',I5,'' NLMAX: '',&
&                           I6,'' NLANGLE:'',I5)',NIDEBCOU, &
&                           NJDEBCOU,NLMAX,NLANGLE
			    print *,' ( Pour le 1D, mettre Obligatoirement ',&
&                           'NLMAX=2 et LPOINTG=T )'
                            IF(ALLOCATED(ZWORK3D))THEN
			      DEALLOCATE(ZWORK3D)
			      LPBREAD=.TRUE.
                            ENDIF
                            RETURN
			  ENDIF

			  ENDIF             !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if(nverbia > 0)then
			  print *,' **oper   lig 3649'
			endif

			  CALL VERIFLEN_FORDIACHRO
			  CALL MEMCV
		          ALLOCATE (ZTEMCV(NLMAX,1:IKU))
			  CALL PRECOU_FORDIACHRO(ZWORK3D,ZTEMCV)
		  if(nverbia >0)THEN
		    print *,' ** oper appel imcou  Ytexte ',YTEXTE(1:LEN_TRIM(YTEXTE))
		  endif
!                         CALL IMCOU_FORDIACHRO(ZTEMCV,XDIAINT,CLEGEND,YTEXTE( &
!                         1:LEN_TRIM(YTEXTE)))

			  IF(LPV)THEN
		            L1DT=.FALSE.
! Janvier 2001
			    IF(LUMVM.OR.LUTVT.OR.LSUMVM.OR.LSUTVT.OR.&
			       LDIRWIND)THEN
			    ILENT=LEN_TRIM(CTITGAL)
			    ILENU=LEN_TRIM(CUNITGAL)
			    YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
			    YTEXTE(ILENT+1:ILENT+1)=' '
			    YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
	                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		              IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		              ALLOCATE(XPRDAT(16,1))
		              CALL LOAD_XPRDAT(1,NLOOPT)
                            ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                            CALL TRACEV_FORDIACHRO(ZTEMCV,KLOOP,YTEXTE(1: &
                                                          LEN_TRIM(YTEXTE)))
	                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL
		                DEALLOCATE(XPRDAT)
                              ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

			    ELSE
! Janvier 2001
			    ALLOCATE(ZTEM1D(IKU),ZWORKZ(IKU))
! Modif AOUT 97
			    ZTEM1D(:)=XSPVAL; ZWORKZ(:)=0.
!                           ZTEM1D(:)=0.; ZWORKZ(:)=0.
			    ZTEM1D(MAX(IKB,NKL):MIN(IKE,NKH))= &
			    ZTEMCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
			    ZWORKZ(:)=XWORKZ(NPROFILE,:,NMGRID)
	                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		              IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		              ALLOCATE(XPRDAT(16,1))
	                      CALL LOAD_XPRDAT(1,NLOOPT)
                            ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
			    CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
	                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	                        DEALLOCATE(XPRDAT)
                              ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
			    ENDIF
			  ELSE IF(LPVT .OR. LPVKT.OR. LPVKT1)THEN
			    L1DT=.FALSE.
		            IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		              ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1) +1
		              ITIMEND=NTIMEDIA(1,KLOOP,1) + &
		              (((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/  &
		              NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
!                             print *,' ITIMEND ',ITIMEND
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
!Fev 2002
                              IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT )THEN
!                             IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT &
!		      .OR.LDIRWIND)THEN
!Fev 2002
			      IF(ALLOCATED(XTEM2D))DEALLOCATE(XTEM2D)
			      IF(ALLOCATED(XTEM2D2))DEALLOCATE(XTEM2D2)
      		              ALLOCATE(XTEM2D(1:IKU,ILENW))
      		              ALLOCATE(XTEM2D2(1:IKU,ILENW))
			      XTEM2D=XSPVAL
			      XTEM2D2=XSPVAL
			      ENDIF
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
		              ALLOCATE(ZTEM2D(1:IKU,ILENW))
                              ALLOCATE(ZWORKT(ILENW))
			      ALLOCATE(ZWORKZ2(IKU))
			      ZWORKZ2(:)=0.; ZWORKT(:)=0.; ZTEM2D(:,:)=0.
			      ZWORKZ2(:)=XWORKZ(NPROFILE,:,NMGRID)
		              ZTEM2D=XSPVAL
	                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		                ALLOCATE(XPRDAT(16,ILENW))
                              ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		              IJLT=0
		            ENDIF
		            IJLT=IJLT+1
	                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
			      CALL LOAD_XPRDAT(IJLT,NLOOPT)
                            ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		            ZWORKT(IJLT)=XTRAJT(JLOOPT,1)
		            ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),IJLT)= &
      		            ZTEMCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
!Fev 2002
                              IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT )THEN
!                             IF(LUMVM .OR.LUTVT .OR.LSUMVM .OR.LSUTVT &
!		      .OR.LDIRWIND)THEN
!Fev 2002
		               XTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),IJLT)= &
      		               ZTEMCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
		               XTEM2D2(MAX(IKB,NKL):MIN(IKE,NKH),IJLT)= &
      		               XWCV(NPROFILE,MAX(IKB,NKL):MIN(IKE,NKH))
			      ENDIF
! Janvier 2001 LUMVM + LDIRWIND + LMUMVM
!                           IF(JLOOPT == NTIMEDIA(2,KLOOP,1))THEN
  			    IF(JLOOPT == ITIMEND)THEN
		              XPVMIN=MINVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		              XPVMAX=MAXVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		              CALL VALMNMX(XPVMIN,XPVMAX)
                              IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
			        XPVMIN=XPVMIN-1.
			        XPVMAX=XPVMAX+1.
                              ENDIF
		              CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
	                      IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		                DEALLOCATE(XPRDAT)
                              ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
			      print *,' **oper lig 3735 AP PVFCT '
		              DEALLOCATE(ZTEM2D,ZWORKT,ZWORKZ2)

			      IF(ALLOCATED(XTEM2D))DEALLOCATE(XTEM2D)
			      IF(ALLOCATED(XTEM2D2))DEALLOCATE(XTEM2D2)
		            ENDIF
			  ELSE
			    ILENT=LEN_TRIM(CTITGAL)
			    ILENU=LEN_TRIM(CUNITGAL)
			    YTEXTE(1:ILENT)=CTITGAL(1:ILENT)
			    YTEXTE(ILENT+1:ILENT+1)=' '
			    YTEXTE(ILENT+2:ILENT+2+ILENU-1)=CUNITGAL(1:ILENU)
	                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		              IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		              ALLOCATE(XPRDAT(16,1))
		              CALL LOAD_XPRDAT(1,NLOOPT)
                            ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                            CALL TRACEV_FORDIACHRO(ZTEMCV,KLOOP,YTEXTE(1: &
                            LEN_TRIM(YTEXTE)))
	                    IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		              DEALLOCATE(XPRDAT)
                            ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
			  ENDIF
			  IF((LCV .OR. LPV) .AND. .NOT. LPVT .AND. .NOT. LPVKT .AND. .NOT. LPVKT1)THEN
!!Fev 2002
                          IF(JLOOPT == NTIMEDIA(2,KLOOP,1))THEN
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
			  ENDIF
!!Fev 2002
                              CALL CLOSF(JLOOPT,ITIMEND, &
			      ISEGD,ISEGM,KLOOP)

			  ENDIF
			  DEALLOCATE(ZTEMCV)
			  DEALLOCATE(XWORKZ,XWZ)
			  IF(ALLOCATED(ZTEM1D))THEN
			    DEALLOCATE(ZTEM1D)
                          ENDIF
			  IF(ALLOCATED(ZWORKZ))THEN
			    DEALLOCATE(ZWORKZ)
                          ENDIF


        	        ENDIF
        	      ENDDO
			  IF((LPVT.AND..NOT.LPBREAD) .OR. LPVKT .OR. LPVKT1)THEN
!                           IF(KLOOP == NSUPERDIA)CALL FRAME
              		    IF(KLOOP == NSUPERDIA)THEN
              		      CALL NGPICT(1,1)
              		      CALL GQACWK(1,IER,INB,IWK)
              		      IF(INB > 1)CALL NGPICT(2,3)
              		    ENDIF
			  ENDIF
        	    ENDIF
        	  ENDDO
        	ENDDO
	                     IF(NVERBIA > 0)THEN
			      print *,' **oper lig 3779  bien AP PVFCT '
			      endif
	        IF(LCV)THEN
                NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
	        IF(NVERBIA > 0)THEN
	        print *,'NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
	        ENDIF
	        ENDIF


	  ENDIF

        ENDIF
!*****************************************************************************
!*****************************************************************************
    CASE('MASK')

      II=SIZE(XVAR,1)
      IJ=SIZE(XVAR,2)
      IK=SIZE(XVAR,3)
      IKU=NKMAX+2*JPVEXT
      IF(CSTORAGE_TYPE == 'PG' .OR. CSTORAGE_TYPE == 'SU')THEN
        IKU=1
      ENDIF
      IKB=1+JPVEXT; IKE=IKU-JPVEXT
      IINF=NIINF;ISUP=NISUP;IJINF=NJINF;IJSUP=NJSUP
      IF(NVERBIA > 0)THEN
      print *,'IINF,ISUP,IJINF,IJSUP ',IINF,ISUP,IJINF,IJSUP
      ENDIF
!     print *,' MASK SIZ XVAR XMASK ',II,IJ,IK,SIZE(XVAR,4),SIZE(XVAR,5), &
!     SIZE(XVAR,6),SIZE(XMASK,1),SIZE(XMASK,2), &
!     SIZE(XMASK,3),SIZE(XMASK,4),SIZE(XMASK,5),SIZE(XMASK,6)

      IF(LCN .OR. LCNCUM .OR. LCNSUM)THEN
!
! Traitement des masques proprement dits (Mot-cle _MASK_ dans la directive)
!
!
! Determination des limites du masque
!
	IF(SIZE(XMASK,1) == NIMAX)THEN
          IF(NIMAX == 1)THEN
            NIINF=1; NISUP=1
            NIL=1; NIH=1
          ELSE
	    NIINF=1+JPHEXT
	    NISUP=NIINF+NIMAX-1
	    NIL=1+JPHEXT; NIH=NIMAX+NIL-1
          ENDIF
	ELSE IF(SIZE(XMASK,1) == NIMAX + 2*JPHEXT)THEN
	  NIINF=1+JPHEXT
	  NISUP=NIINF+NIMAX-1
	  NIL=1; NIH=NIMAX+2*JPHEXT
	ELSE
	  print *,' Taille des masques en X differente de IIU OU IMAX ', &
	  SIZE(XMASK,1)
	  print *,' PAS DE TRACE '
	  RETURN
	ENDIF
	IF(SIZE(XMASK,2) == NJMAX)THEN
          IF(NJMAX == 1)THEN
            NJINF=1; NJSUP=1
            NJL=1; NJH=1
          ELSE
	    NJINF=1+JPHEXT
	    NJSUP=NJINF+NJMAX-1
	    NJL=1+JPHEXT; NJH=NJMAX+NJL-1
          ENDIF
	ELSE IF(SIZE(XMASK,2) == NJMAX + 2*JPHEXT)THEN
	  NJINF=1+JPHEXT
	  NJSUP=NJINF+NJMAX-1
	  NJL=1; NJH=NJMAX+2*JPHEXT
	ELSE
	  print *,' Taille des masques en Y differente de IJU OU JMAX ', &
	  SIZE(XMASK,2)
	  print *,' PAS DE TRACE '
	ENDIF
	ALLOCATE(ZWORK3D(NISUP-NIINF+1,NJSUP-NJINF+1,1))
	ZWORK3D=0.
	CTYPHOR(1:LEN(CTYPHOR))=' '
	CTYPHOR='K'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  IF(LCN)THEN                  !!!!!!!!!!!!!!!!!!!!
! *****************************************
! Boucle externe sur les numeros de masques
! *****************************************
!
            DO JLOOPN=1,NBNDIA(KLOOP)         !......................1
              NLOOPN=NNDIA(JLOOPN,KLOOP)
	    NMGRID=1
	    YC1=' '; YC2='  '
	    IF(NNDIA(JLOOPN,KLOOP) < 10)THEN
	      WRITE(YC1,'(I1)')NNDIA(JLOOPN,KLOOP)
	      CTITGAL='MASK'//YC1
            ELSE
	      WRITE(YC2,'(I2)')NNDIA(JLOOPN,KLOOP)
	      CTITGAL='MASK'//YC2
	    ENDIF
	    CTITGAL=ADJUSTL(ADJUSTR(CTITGAL))
              IF(.NOT.LTINCRDIA(KLOOP,1))THEN
!
! ***********************************************
! Boucle sur les temps (Formulation sequentielle)
! ***********************************************
!
                DO JLOOPT=1,NBTIMEDIA(KLOOP,1)        !................2
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
! Juillet 2001
		  IF(LANIMT .AND. NISUP-NIINF /= 0 .AND. NJSUP-NJINF /= 0)THEN
		    IF(JLOOPT == 1)THEN
        	      CALL FMFREE(YBID,YBID,IRESP)
		      print *,' OPER FMFREE YBID IRESP ',YBID,IRESP

		      CALL FMATTR(YBID,YBID,IBID,IRESP)
		      CALL GOPWK(9,IBID,3)
!                     CALL GOPWK(9,20,3)
		      ISEGM=ISEGM+1
		      ISEGD=ISEGM
		      CALL GFLAS1(ISEGM)
                    ELSE
		      ISEGM=ISEGM+1
		      CALL GFLAS1(ISEGM)
                    ENDIF
                  ENDIF
! Juillet 2001
		  CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
	          ZWORK3D(:,:,1)=XMASK(NIINF-NIL+1:NISUP-NIL+1,  &
		  NJINF-NJL+1:NJSUP-NJL+1,1, &
		  NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP),1)
! Traitement cas 2D (--> masque filaire)
                  IF(NIINF == 1 .AND. NISUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(1:1,:,1:1),KLOOP)
                  ELSE IF(NJINF == 1 .AND. NJSUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(:,1:1,1:1),KLOOP)
                  ELSE
! Traitement cas 3D (--> masque surfacique)
	            CALL TRACEH_FORDIACHRO(1,ZWORK3D,KLOOP)
                  ENDIF
! Juillet 2001
                  IF(LANIMT .AND. NISUP-NIINF /= 0 .AND. NJSUP-NJINF /= 0)THEN
		    CALL GFLAS2
		    IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
		      DO JJ=ISEGD,ISEGM
			CALL GFLAS3(JJ)
		      ENDDO
		      CALL GCLWK(9)
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
                    ENDIF
                    ELSE
! Juillet 2001
		  IF(KLOOP == NSUPERDIA)CALL FRAME
! Juillet 2001
                    ENDIF
! Juillet 2001
	        ENDDO                               !................2
	      ELSE
!
! ***********************************************
! Boucle sur les temps (Formulation incrementale)
! ***********************************************
!
	        DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)!.3
		  NLOOPT=JLOOPT
! Juillet 2001
		  IF(LANIMT .AND. NJSUP-NJINF /= 0 .AND. NISUP-NIINF /=0)THEN
		    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
        	      CALL FMFREE(YBID,YBID,IRESP)
                      if(nverbia >0)then
		      print *,' OPER FMFREE YBID IRESP ',YBID,IRESP
                      endif
		      CALL FMATTR(YBID,YBID,IBID,IRESP)
                      if(nverbia >0)then
		      print *,' OPER FMATTR YBID IBID IRESP ',YBID,IBID,IRESP
                      endif
		      CALL GOPWK(9,IBID,3)
		      ISEGM=ISEGM+1
		      ISEGD=ISEGM
		      CALL GFLAS1(ISEGM)
		      ITIMEND=NTIMEDIA(1,KLOOP,1) + &
		      (((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/  &
		      NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
                    ELSE
		      ISEGM=ISEGM+1
		      print *,' OPER ISEGM ',ISEGM
		      CALL GFLAS1(ISEGM)
                    ENDIF
		  ENDIF
! Juillet 2001
		  CALL RESOLV_TIMES(JLOOPT)
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)
	          ZWORK3D(:,:,1)=XMASK(NIINF-NIL+1:NISUP-NIL+1,  &
		  NJINF-NJL+1:NJSUP-NJL+1,1,JLOOPT, &
		  NNDIA(JLOOPN,KLOOP),1)
! Traitement cas 2D (--> masque filaire)
                  IF(NIINF == 1 .AND. NISUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(1:1,:,1:1),KLOOP)
                  ELSE IF(NJINF == 1 .AND. NJSUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(:,1:1,1:1),KLOOP)
                  ELSE
! Traitement cas 3D (--> masque surfacique)
	            CALL TRACEH_FORDIACHRO(1,ZWORK3D,KLOOP)
                  ENDIF
! Juillet 2001
		  IF(LANIMT .AND. NISUP-NIINF /= 0 .AND. NJSUP-NJINF /= 0)THEN
		    CALL GFLAS2
		    IF(JLOOPT == ITIMEND)THEN
		      DO JJ=ISEGD,ISEGM
                        CALL GFLAS3(JJ)
                      ENDDO 
		      CALL GCLWK(9)
		      CALL NGPICT(1,1)
		      CALL GQACWK(1,IER,INB,IWK)
		      IF(INB > 1)CALL NGPICT(2,3)
		    ENDIF
                   ELSE
! Juillet 2001
		  IF(KLOOP == NSUPERDIA)CALL FRAME
! Juillet 2001
                    ENDIF
! Juillet 2001
	        ENDDO                                                          !.3
	      ENDIF
            ENDDO                             !......................1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ELSE IF(LCNCUM)THEN          !!!!!!!!!!!!!!!!!!!!

! *****************************************
! Boucle externe sur les numeros de masques
! *****************************************
!
            DO JLOOPN=1,NBNDIA(KLOOP)         !......................1
              NLOOPN=NNDIA(JLOOPN,KLOOP)
	    NMGRID=1
	    ZWORK3D=0.
	    YC1=' '; YC2='  '
	    IF(NNDIA(JLOOPN,KLOOP) < 10)THEN
	      WRITE(YC1,'(I1)')NNDIA(JLOOPN,KLOOP)
	      CTITGAL='MASK'//YC1
            ELSE
	      WRITE(YC2,'(I2)')NNDIA(JLOOPN,KLOOP)
	      CTITGAL='MASK'//YC2
	    ENDIF
	    CTITGAL=ADJUSTL(ADJUSTR(CTITGAL))
              IF(.NOT.LTINCRDIA(KLOOP,1))THEN
!
! ***********************************************
! Boucle sur les temps (Formulation sequentielle)
! ***********************************************
!
                IJLT=0
                DO JLOOPT=1,NBTIMEDIA(KLOOP,1)        !................2
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
		  IJLT=IJLT+1
		  CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                  IF(IJLT < 9)THEN
		    WRITE(CTIMECS(8*IJLT:8*IJLT+7),'(F8.0)')XTRAJT(  &
		    NTIMEDIA(JLOOPT,KLOOP,1),1)
                  ELSE IF(IJLT == 9)THEN
                    CTIMECS(8*IJLT:8*IJLT+4)='.....'
                  ENDIF
	          ZWORK3D(:,:,1)=ZWORK3D(:,:,1) + XMASK(NIINF-NIL+1:  &
		  NISUP-NIL+1,NJINF-NJL+1:NJSUP-NJL+1,1,  &
		  NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP),1)
!                 print *,' JLOOPT JLOOPN ZWORK3D ',JLOOPT,JLOOPN
!                 print *,ZWORK3D(:,:,1)
	        ENDDO                               !................2
! Traitement cas 2D (--> masque filaire)
                  IF(NIINF == 1 .AND. NISUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(1:1,:,1:1),KLOOP)
                  ELSE IF(NJINF == 1 .AND. NJSUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(:,1:1,1:1),KLOOP)
                    CALL EZXY(XXX(NIINF:NISUP,NMGRID),ZWORK3D(:,1,1), &
                    NISUP-NIINF+1,0)
                  ELSE
! Traitement cas 3D (--> masque surfacique)
	            CALL TRACEH_FORDIACHRO(1,ZWORK3D,KLOOP)
                  ENDIF
		  IF(KLOOP == NSUPERDIA)CALL FRAME
	      ELSE
!
! ***********************************************
! Boucle sur les temps (Formulation incrementale)
! ***********************************************
!
                IJLT=0
	        DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)!.3
		  NLOOPT=JLOOPT
		  IJLT=IJLT+1
		  CALL RESOLV_TIMES(JLOOPT)
                  IF(IJLT < 9)THEN
		    WRITE(CTIMECS(8*IJLT:8*IJLT+7),'(F8.0)')XTRAJT(  &
		    JLOOPT,1)
                  ELSE IF(IJLT == 9)THEN
                    CTIMECS(8*IJLT:8*IJLT+4)='.....'
                  ENDIF
	          ZWORK3D(:,:,1)=ZWORK3D(:,:,1) + XMASK(NIINF-NIL+1:  &
		  NISUP-NIL+1,NJINF-NJL+1:NJSUP-NJL+1,1,JLOOPT,  &
		  NNDIA(JLOOPN,KLOOP),1)
!                 print *,' JLOOPT JLOOPN ZWORK3D ',JLOOPT,JLOOPN
!                 print *,ZWORK3D(:,:,1)
	        ENDDO                                                          !.3
! Traitement cas 2D (--> masque filaire)
                  IF(NIINF == 1 .AND. NISUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(1:1,:,1:1),KLOOP)
                  ELSE IF(NJINF == 1 .AND. NJSUP == 1)THEN
                    CALL TRAMASK(ZWORK3D(:,1:1,1:1),KLOOP)
                  ELSE
! Traitement cas 3D (--> masque surfacique)
	            CALL TRACEH_FORDIACHRO(1,ZWORK3D,KLOOP)
                  ENDIF
		  IF(KLOOP == NSUPERDIA)CALL FRAME
	      ENDIF
            ENDDO                             !......................1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ELSE IF(LCNSUM)THEN          !!!!!!!!!!!!!!!!!!!!

! *****************************************
! Boucle externe sur les numeros de masques
! *****************************************
!
            DO JLOOPN=1,NBNDIA(KLOOP)         !......................1
              NLOOPN=NNDIA(JLOOPN,KLOOP)
	    NMGRID=1
	    YC1=' '; YC2='  '
	    IF(NNDIA(JLOOPN,KLOOP) < 10)THEN
	      WRITE(YC1,'(I1)')NNDIA(JLOOPN,KLOOP)
	      CTITGAL='MASK'//YC1
            ELSE
	      WRITE(YC2,'(I2)')NNDIA(JLOOPN,KLOOP)
	      CTITGAL='MASK'//YC2
	    ENDIF
	    CTITGAL=ADJUSTL(ADJUSTR(CTITGAL))
              IF(.NOT.LTINCRDIA(KLOOP,1))THEN
!
! ***********************************************
! Boucle sur les temps (Formulation sequentielle)
! ***********************************************
!
		ALLOCATE(ZWORK1D(NBTIMEDIA(KLOOP,1)))
		ALLOCATE(ZWORKT(NBTIMEDIA(KLOOP,1)))
                IJLT=0
                DO JLOOPT=1,NBTIMEDIA(KLOOP,1)        !................2
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
                  IJLT=IJLT+1
		  CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                  IF(IJLT < 9)THEN
		    WRITE(CTIMECS(8*IJLT:8*IJLT+7),'(F8.0)')XTRAJT(  &
		    NTIMEDIA(JLOOPT,KLOOP,1),1)
                  ELSE IF(IJLT == 9)THEN
                    CTIMECS(8*IJLT:8*IJLT+4)='.....'
                  ENDIF
	          ZWORK3D(:,:,1)=XMASK(NIINF-NIL+1:NISUP-NIL+1,  &
		  NJINF-NJL+1:NJSUP-NJL+1,1,NTIMEDIA(JLOOPT,KLOOP,1),  &
		  NNDIA(JLOOPN,KLOOP),1)
		  ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
		  ZWORK1D(JLOOPT)=SUM(ZWORK3D)
	        ENDDO                               !................2
		LFT1=.TRUE.
		CALL VARFCT(ZWORKT,ZWORK1D,1)
		  IF(KLOOP == NSUPERDIA)CALL FRAME
	      ELSE
!
! ***********************************************
! Boucle sur les temps (Formulation incrementale)
! ***********************************************
!
		ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1)+1
		ALLOCATE(ZWORKT(ILENW))
		ALLOCATE(ZWORK1D(ILENW))
                IJLT=0
	        DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)!.3
		  NLOOPT=JLOOPT
                  IJLT=IJLT+1
		  CALL RESOLV_TIMES(JLOOPT)
                  IF(IJLT < 9)THEN
		    WRITE(CTIMECS(8*IJLT:8*IJLT+7),'(F8.0)')XTRAJT(  &
		    JLOOPT,1)
                  ELSE IF(IJLT == 9)THEN
                    CTIMECS(8*IJLT:8*IJLT+4)='.....'
                  ENDIF
	          ZWORK3D(:,:,1)=XMASK(NIINF-NIL+1:NISUP-NIL+1,  &
		  NJINF-NJL+1:NJSUP-NJL+1,1,  &
		  JLOOPT,NNDIA(JLOOPN,KLOOP),1)
!                 print *,' OPER OPER JLOOPT ZWORK3D ',JLOOPT
!                 print *,ZWORK3D
! Correction AOUT 2001
		  ZWORKT(IJLT)=XTRAJT(JLOOPT,1)
		  ZWORK1D(IJLT)=SUM(ZWORK3D)
!                 ZWORKT(JLOOPT)=XTRAJT(JLOOPT,1)
!                 ZWORK1D(JLOOPT)=SUM(ZWORK3D)
	        ENDDO                                                          !.3
		LFT1=.TRUE.
		CALL VARFCT(ZWORKT,ZWORK1D,1)
		  IF(KLOOP == NSUPERDIA)CALL FRAME
	      ENDIF
	      DEALLOCATE(ZWORKT,ZWORK1D)
            ENDDO                             !......................1
	  ENDIF
	DEALLOCATE(ZWORK3D)
      ELSE
!
! Traitement des infos gerees par un masque:  PV
! Cas compression sur l'axe Z (Compressions en X et Y implicites)
! ***************************************************************
        DO JLOOPN=1,NBNDIA(KLOOP)         !......................1
              NLOOPN=NNDIA(JLOOPN,KLOOP)

   	  IF(LPVKT .AND. NSUPERDIA>1)THEN
            IF(NBPROCDIA(KLOOP)>1 .OR. NBLVLKDIA(KLOOP,1)>1 &
	    .OR. NBNDIA(KLOOP)>1)THEN
              print *,' _PVKT_  SUPERPOSITIONS : '
              print *,'         On ne peut definir de part de d''autre '&
             &'de _ON_ qu''1 seul processus, 1 seul niveau, 1 seule station '
              print *,' Nb de niveaux demandes   : ',NBLVLKDIA(KLOOP,1)
              print *,' Nb de processus demandes : ',NBPROCDIA(KLOOP)
              print *,' Nb de stations demandees : ',NBNDIA(KLOOP)
              print *,' *** MODIFIEZ VOTRE DIRECTIVE *** '
              EXIT
            ENDIF
          ENDIF

        YTITGAL(1:LEN(YTITGAL))=' '
	YC1=' '; YC2='  '
	IF(NLOOPN < 10)THEN
	  WRITE(YC1,'(I1)')NNDIA(JLOOPN,KLOOP)
	  YTITGAL='MASK'//YC1
        ELSE
	  WRITE(YC2,'(I2)')NNDIA(JLOOPN,KLOOP)
	  YTITGAL='MASK'//YC2
	ENDIF
	YTITGAL=ADJUSTL(ADJUSTR(YTITGAL))
        IF(II == 1 .AND. IJ == 1 .AND. IK  == 1)THEN
          IF(.NOT.LTINCRDIA(KLOOP,1))THEN
            ALLOCATE(ZWORKT(NBTIMEDIA(KLOOP,1)))
            ALLOCATE(ZWORK1D(NBTIMEDIA(KLOOP,1)))
            DO JLOOPP=1,NBPROCDIA(KLOOP)
              NLOOPP=NPROCDIA(JLOOPP,KLOOP)
	      NMGRID=NGRIDIA(NPROCDIA(JLOOPP,KLOOP))
		    IF(NGRIDIAM /= 0 .AND. (NMGRID /= NGRIDIAM))THEN
		      print *,' ****oper NMGRID Av modif ',NMGRID
		      NMGRID=NGRIDIAM
		      print *,' ****oper NMGRID mis volontairement a la valeur de NGRIDIAM ',NGRIDIAM
		    ENDIF
	        IF(NMGRID <1 .OR. NMGRID >7)THEN
	          PRINT *,' VALEUR NMGRID ABERRANTE: ',NMGRID, &
                    '        FORCEE A        :  1'
                  NMGRID=1
                ENDIF
              CTITGAL(1:LEN(CTITGAL))=' '
              CTITGAL=ADJUSTL(ADJUSTR(YTITGAL)//' '//ADJUSTL(CTITRE(NPROCDIA(JLOOPP,KLOOP))))
	      CUNITGAL=ADJUSTL(CUNITE(NPROCDIA(JLOOPP,KLOOP)))
	      CTITGAL=ADJUSTL(CTITGAL)
!             print *,' MASK PV JLOOPP, NPROCDIA CTITGAL ',JLOOPP,NPROCDIA(JLOOPP,KLOOP),' ',CTITGAL
	      CUNITGAL=ADJUSTL(CUNITGAL)
		  CUNITGAL(INDEX(CUNITGAL,' '):LEN(CUNITGAL))=' '
	      DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
	        CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
	        ZWORK1D(JLOOPT)=XVAR(1,1,1,NTIMEDIA(JLOOPT,KLOOP,1), &
		NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP,KLOOP))
	      ENDDO
              CALL VARFCT(ZWORKT,ZWORK1D,1)
              IF(KLOOP == NSUPERDIA)CALL FRAME
              DEALLOCATE(ZWORKT,ZWORK1D)
            ENDDO
          ELSE
            ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1)+1
            ALLOCATE(ZWORKT(ILENW))
            ALLOCATE(ZWORK1D(ILENW))
            DO JLOOPP=1,NBPROCDIA(KLOOP)
              NLOOPP=NPROCDIA(JLOOPP,KLOOP)
              NMGRID=NGRIDIA(NPROCDIA(JLOOPP,KLOOP))
		    IF(NGRIDIAM /= 0 .AND. (NMGRID /= NGRIDIAM))THEN
		      print *,' ****oper NMGRID Av modif ',NMGRID
		      NMGRID=NGRIDIAM
		      print *,' ****oper NMGRID mis volontairement a la valeur de NGRIDIAM ',NGRIDIAM
		    ENDIF
      	      IF(NMGRID <1 .OR. NMGRID >7)THEN
		PRINT *,' VALEUR NMGRID ABERRANTE: ',NMGRID, &
                    '        FORCEE A        :  1'
                NMGRID=1
              ENDIF
              CTITGAL(1:LEN(CTITGAL))=' '
              CTITGAL=ADJUSTL(ADJUSTR(YTITGAL)//' '//ADJUSTL(CTITRE(NPROCDIA(JLOOPP,KLOOP))))
              CUNITGAL=ADJUSTL(CUNITE(NPROCDIA(JLOOPP,KLOOP)))
              CTITGAL=ADJUSTL(CTITGAL)
              CUNITGAL=ADJUSTL(CUNITGAL)
		  CUNITGAL(INDEX(CUNITGAL,' '):LEN(CUNITGAL))=' '
!             print *,' MASK PV JLOOPP, NPROCDIA CTITGAL ',JLOOPP,NPROCDIA(JLOOPP,KLOOP),' ',CTITGAL
              IJLT=0
              DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		NLOOPT=JLOOPT
                CALL RESOLV_TIMES(JLOOPT)
                IJLT=IJLT+1
                ZWORKT(IJLT)=XTRAJT(JLOOPT,1)
                ZWORK1D(IJLT)=XVAR(1,1,1,JLOOPT,NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP,KLOOP))
	      ENDDO
              CALL VARFCT(ZWORKT,ZWORK1D,1)
              IF(KLOOP == NSUPERDIA)CALL FRAME
              DEALLOCATE(ZWORKT,ZWORK1D)
            ENDDO
          ENDIF
        ELSE IF(II == 1 .AND. IJ == 1 .AND. IK /= 1)THEN
! Pas de compression en Z
! ***********************
	      L1DT=.TRUE.
	      ALLOCATE(ZTEM1D(IKU),ZWORKZ(IKU))
              DO JLOOPP=1,NBPROCDIA(KLOOP)
!!! Octobre 2001
                IF(JLOOPP > 1 .AND. LUMVMPV .AND. LPV)EXIT
!!! Octobre 2001
                NLOOPP=NPROCDIA(JLOOPP,KLOOP)
! Modif AOUT 97
	        ZTEM1D(:)=XSPVAL; ZWORKZ(:)=0.
!               ZTEM1D(:)=0.; ZWORKZ(:)=0.
                NMGRID=NGRIDIA(NPROCDIA(JLOOPP,KLOOP))
		    IF(NGRIDIAM /= 0 .AND. (NMGRID /= NGRIDIAM))THEN
		      print *,' ****oper NMGRID Av modif ',NMGRID
		      NMGRID=NGRIDIAM
		      print *,' ****oper NMGRID mis volontairement a la valeur de NGRIDIAM ',NGRIDIAM
		    ENDIF
                IF(NMGRID <1 .OR. NMGRID >7)THEN
	          PRINT *,' VALEUR NMGRID ABERRANTE: ',NMGRID, &
                            '        FORCEE A        :  1'
                  NMGRID=1
                ENDIF
!!!!!!!!!!Octobre 2001
             IF(LUMVMPV)THEN
               NMGRID=1
             ENDIF
!!!!!!!!!!Octobre 2001
	        CALL COMPCOORD_FORDIACHRO(NMGRID)
                CTITGAL(1:LEN(CTITGAL))=' '
                CTITGAL=ADJUSTL(ADJUSTR(YTITGAL)//' '//ADJUSTL(CTITRE(NPROCDIA(JLOOPP,KLOOP))))
!               CTITGAL=ADJUSTL(CTITRE(NPROCDIA(JLOOPP,KLOOP)))
                CUNITGAL=ADJUSTL(CUNITE(NPROCDIA(JLOOPP,KLOOP)))
                CTITGAL=ADJUSTL(CTITGAL)
!             print *,' MASK PV JLOOPP, NPROCDIA CTITGAL ',JLOOPP,NPROCDIA(JLOOPP,KLOOP),' ',CTITGAL
!              print *,' MASK CTITRE ',CTITRE(NPROCDIA(JLOOPP,KLOOP))
                CUNITGAL=ADJUSTL(CUNITGAL)
		  CUNITGAL(INDEX(CUNITGAL,' '):LEN(CUNITGAL))=' '
! Expression temps non incrementale
		IF(.NOT.LTINCRDIA(KLOOP,1))THEN
                DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
	          CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
                  CTIMEC(16:16)='s'
	          ZTEM1D(NKL:NKH)=XVAR(1,1,: &
		  ,NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP),&
		  NPROCDIA(JLOOPP,KLOOP))
		  ZWORKZ(:)=XXZ(:,NMGRID)
!                 print *,' ZTEM1D '
!                 print *,ZTEM1D
!!!!!!!!!!Octobre 2001
                  
!!!!!!!!!!Octobre 2001
    	          IF(LPV)THEN
!!!!!!!!!!Octobre 2001
                    IF(LUMVMPV)THEN
		      LPV=.FALSE. ; LPVT=.TRUE.
                      IF(JLOOPP == 1)THEN
                        ILENW=1
                        ALLOCATE(ZTEM2D(1:IKU,ILENW))
                        ALLOCATE(ZWORKT(ILENW))
                        ZWORKT=NLOOPT
                        IF(ALLOCATED(XTEM2D))THEN
                           DEALLOCATE(XTEM2D)
                        ENDIF
                       ALLOCATE(XTEM2D(1:IKU,ILENW))
                       XTEM2D=XSPVAL
                       IF(ALLOCATED(XTEM2D2))THEN
                         DEALLOCATE(XTEM2D2)
                       ENDIF
                       ALLOCATE(XTEM2D2(1:IKU,ILENW))
                       XTEM2D2=XSPVAL
                       XTEM2D(:,1)=ZTEM1D
                       XTEM2D2(NKL:NKH,1)=XVAR(1,1,: &
                       ,NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP+1,KLOOP))
                        IF(NBPROCDIA(KLOOP) == 3)THEN
			  ZTEM2D=XSPVAL
                          ZTEM2D(NKL:NKH,1)=XVAR(1,1,: &
                        ,NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP+2,KLOOP))
                        CALL COLVECT(IKU,ZTEM2D)
                       ENDIF
    
                        CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
			IF(LUMVMPV)THEN
		          LPV=.TRUE. ; LPVT=.FALSE.
			ENDIF
                        DEALLOCATE(ZTEM2D,ZWORKT)
                        IF(ALLOCATED(XTEM2D))THEN
                          DEALLOCATE(XTEM2D)
                        ENDIF
                        IF(ALLOCATED(XTEM2D2))THEN
                          DEALLOCATE(XTEM2D2)
                        ENDIF
                        LCOLPVT=.FALSE.
                      ENDIF
                  ELSE
!!!!!!!!!!Octobre 2001
		    CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
!!!!!!!!!!Octobre 2001
                  ENDIF
!!!!!!!!!!Octobre 2001
		    IF(KLOOP == NSUPERDIA)CALL FRAME
    	          ELSE IF(LPVT .OR. LPVKT)THEN
		    IF(JLOOPT == 1)THEN
		      ILENW=NBTIMEDIA(KLOOP,1)
                      IF(ALLOCATED(ZTEM2D))THEN
                        DEALLOCATE(ZTEM2D)
                      ENDIF
                      IF(ALLOCATED(ZWORKT))THEN
                        DEALLOCATE(ZWORKT)
                      ENDIF
		      ALLOCATE(ZTEM2D(1:IKU,ILENW))
		      ZTEM2D=XSPVAL
		      ALLOCATE(ZWORKT(ILENW))
!!!!!!!!!!Octobre 2001
                     IF(LUMVM)THEN
                        IF(ALLOCATED(XTEM2D))THEN
                          DEALLOCATE(XTEM2D)
                        ENDIF
                        ALLOCATE(XTEM2D(1:IKU,ILENW))
                        XTEM2D=XSPVAL
                      ENDIF
                      IF(LUMVMPV .AND. JLOOPP == 1)THEN
                        IF(ALLOCATED(XTEM2D))THEN
                          DEALLOCATE(XTEM2D)
                        ENDIF
                        ALLOCATE(XTEM2D(1:IKU,ILENW))
                        XTEM2D=XSPVAL
                        IF(ALLOCATED(XTEM2D2))THEN
                          DEALLOCATE(XTEM2D2)
                        ENDIF
                        ALLOCATE(XTEM2D2(1:IKU,ILENW))
                        XTEM2D2=XSPVAL
                      ENDIF
!!!!!!!!!!Octobre 2001
		    ENDIF
		      ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
		      ZTEM2D(NKL:NKH,NTIMEDIA(JLOOPT,KLOOP,1))= &
		      XVAR(1,1,:,  &
		      NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP,KLOOP))
		    IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
		      XPVMIN=MINVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      XPVMAX=MAXVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      CALL VALMNMX(XPVMIN,XPVMAX)
                      IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
			XPVMIN=XPVMIN-1.
			XPVMAX=XPVMAX+1.
                      ENDIF
!!!!!!!!!!Octobre 2001
                      IF(LUMVMPV)THEN
                        IF(JLOOPP == 1)THEN
! Memorisation de U
                          XTEM2D=ZTEM2D
                          CYCLE
                        ELSEIF(JLOOPP == 2)THEN
                          IF(JLOOPP == NBPROCDIA(KLOOP))THEN
                            XTEM2D2=ZTEM2D
                          ELSE
                            XTEM2D2=ZTEM2D
                            CYCLE
                          ENDIF
                        ELSEIF(JLOOPP == 3)THEN
                          CALL COLVECT(IKU,ZTEM2D)
                        ENDIF
                      ENDIF

!!!!!!!!!!Octobre 2001
		      CALL COMPCOORD_FORDIACHRO(NMGRID)
		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
		      DEALLOCATE(ZTEM2D,ZWORKT)
!!!!!!!!!!Octobre 2001
                      IF(ALLOCATED(XTEM2D))THEN
                        DEALLOCATE(XTEM2D)
                      ENDIF
                      IF(ALLOCATED(XTEM2D2))THEN
                        DEALLOCATE(XTEM2D2)
                      ENDIF
                      LCOLPVT=.FALSE.

!!!!!!!!!!Octobre 2001
		      IF(.NOT.LPBREAD)THEN
		        IF(KLOOP == NSUPERDIA)CALL FRAME
		      ENDIF
		    ENDIF
    	          ENDIF
	        ENDDO
		ELSE
! Expression temps incrementale
!               print *,'NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1) ', &
!                        NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
!               print *,XTIMEDIA(1,KLOOP,1),XTIMEDIA(2,KLOOP,1),XTIMEDIA(3,KLOOP,1)
                DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
		  NLOOPT=JLOOPT
	          CALL RESOLV_TIMES(JLOOPT)
                  WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)
	          ZTEM1D(NKL:NKH)=XVAR(1,1,: &
		  ,JLOOPT,NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP,KLOOP))
		  ZWORKZ(:)=XXZ(:,NMGRID)
!                 print *,' ZTEM1D '
!                 print *,ZTEM1D

    	          IF(LPV)THEN
!!! Octobre 2001
                    IF(LUMVMPV)THEN
		      LPV=.FALSE. ; LPVT=.TRUE.
                      IF(JLOOPP == 1)THEN
                        ILENW=1
                        ALLOCATE(ZTEM2D(1:IKU,ILENW))
                        ALLOCATE(ZWORKT(ILENW))
                        ZWORKT=NLOOPT
                        IF(ALLOCATED(XTEM2D))THEN
                        DEALLOCATE(XTEM2D)
                        ENDIF
                        ALLOCATE(XTEM2D(1:IKU,ILENW))
                        XTEM2D=XSPVAL
                        IF(ALLOCATED(XTEM2D2))THEN
                        DEALLOCATE(XTEM2D2)
                        ENDIF
                        ALLOCATE(XTEM2D2(1:IKU,ILENW))
                        XTEM2D2=XSPVAL
                        XTEM2D(:,1)=ZTEM1D
                        XTEM2D2(NKL:NKH,1)=XVAR(1,1,: &
                        ,JLOOPT,NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP+1,KLOOP))
                        IF(NBPROCDIA(KLOOP) == 3)THEN
			  ZTEM2D=XSPVAL
                          ZTEM2D(NKL:NKH,1)=XVAR(1,1,: &
                          ,JLOOPT,NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP+2,KLOOP))
                        
                          CALL COLVECT(IKU,ZTEM2D)
                        ENDIF
                        CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
			IF(LUMVMPV)THEN
			  LPV=.TRUE. ; LPVT=.FALSE.
			ENDIF
                        DEALLOCATE(ZTEM2D,ZWORKT)
                        IF(ALLOCATED(XTEM2D))THEN
                          DEALLOCATE(XTEM2D)
                        ENDIF
                        IF(ALLOCATED(XTEM2D2))THEN
                          DEALLOCATE(XTEM2D2)
                        ENDIF
                        LCOLPVT=.FALSE.
                      ENDIF

                    ELSE
!!! Octobre 2001
		    CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
!!! Octobre 2001
                    ENDIF
!!! Octobre 2001
		    IF(KLOOP == NSUPERDIA)CALL FRAME

    	          ELSE IF(LPVT .OR. LPVKT)THEN

		    IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		      ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1) +1
		      ITIMEND=NTIMEDIA(1,KLOOP,1) + &
		      (((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/  &
		      NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
                      IF(NVERBIA > 0)THEN
		      print *,' ITIMEND ',ITIMEND
		      ENDIF
                      IF(ALLOCATED(ZTEM2D))THEN
                        DEALLOCATE(ZTEM2D)
                      ENDIF
                      IF(ALLOCATED(ZWORKT))THEN
                        DEALLOCATE(ZWORKT)
                      ENDIF

		      ALLOCATE(ZTEM2D(1:IKU,ILENW))
		      ZTEM2D=XSPVAL
                      ALLOCATE(ZWORKT(ILENW))
		      IJLT=0
!!!!!!!!!!Octobre 2001
                      IF(LUMVMPV .AND. JLOOPP == 1)THEN
                        IF(ALLOCATED(XTEM2D))THEN
                          DEALLOCATE(XTEM2D)
                        ENDIF
                        ALLOCATE(XTEM2D(1:IKU,ILENW))
                        XTEM2D=XSPVAL
                        IF(ALLOCATED(XTEM2D2))THEN
                          DEALLOCATE(XTEM2D2)
                        ENDIF
                        ALLOCATE(XTEM2D2(1:IKU,ILENW))
                        XTEM2D2=XSPVAL
                      ENDIF

!!!!!!!!!!Octobre 2001
		    ENDIF
		    IJLT=IJLT+1
		    ZWORKT(IJLT)=XTRAJT(JLOOPT,1)
		    ZTEM2D(NKL:NKH,IJLT)= &
		    XVAR(1,1,:,  &
		    JLOOPT,NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP,KLOOP))

!                   IF(JLOOPT == NTIMEDIA(2,KLOOP,1))THEN
		    IF(JLOOPT == ITIMEND)THEN
		      XPVMIN=MINVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      XPVMAX=MAXVAL(ZTEM2D(MAX(IKB,NKL):MIN(IKE,NKH),:))
		      CALL VALMNMX(XPVMIN,XPVMAX)
                      IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
			XPVMIN=XPVMIN-1.
			XPVMAX=XPVMAX+1.
                      ENDIF
!!!!!!!!!!Octobre 2001
                      IF(LUMVMPV)THEN        !llllllllllllllllllll

                        IF(JLOOPP == 1)THEN  !kkkkkkkkkkkkkkkkkkkkkkk
! Memorisation de U
                          XTEM2D=ZTEM2D
                          CYCLE
                        ELSEIF(JLOOPP == 2)THEN !kkkkkkkkkkkkkkkkkkkkk
                          IF(JLOOPP == NBPROCDIA(KLOOP))THEN
                            XTEM2D2=ZTEM2D
                          ELSE
                            XTEM2D2=ZTEM2D
                            CYCLE
                          ENDIF
                        ELSEIF(JLOOPP == 3)THEN !kkkkkkkkkkkkkkkkkkkkk
                          CALL COLVECT(IKU,ZTEM2D)
                        ENDIF         !kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
                      ENDIF           !llllllllllllllllllllllllllllllllll


  
!!!!!!!!!!Octobre 2001
		      CALL COMPCOORD_FORDIACHRO(NMGRID)
		      CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
		      DEALLOCATE(ZTEM2D,ZWORKT)
!!!!!!!!!!Octobre 2001
                      IF(ALLOCATED(XTEM2D))THEN
                        DEALLOCATE(XTEM2D)
                      ENDIF
                      IF(ALLOCATED(XTEM2D2))THEN
                        DEALLOCATE(XTEM2D2)
                      ENDIF
                      LCOLPVT=.FALSE.

!!!!!!!!!!Octobre 2001

		      IF(.NOT.LPBREAD)THEN
		        IF(KLOOP == NSUPERDIA)CALL FRAME
		      ENDIF
		    ENDIF
    	          ENDIF
	        ENDDO
		ENDIF
	      ENDDO
	      DEALLOCATE(ZTEM1D,ZWORKZ)
        ELSE
        ENDIF
	ENDDO
      ENDIF 
      NIINF=IINF;NISUP=ISUP;NJINF=IJINF;NJSUP=IJSUP
      IF(NVERBIA > 0)THEN
      print *,'NIINF,NISUP,NJINF,NJSUP ',NIINF,NISUP,NJINF,NJSUP
      ENDIF

!*****************************************************************************
!*****************************************************************************
    CASE('SSOL')
!
! ******************************************
! Boucle externe sur les numeros de stations
! ******************************************
!
      L1DT=.TRUE.
      DO JLOOPN=1,NBNDIA(KLOOP)
        NLOOPN=NNDIA(JLOOPN,KLOOP)

   	  IF(LPVKT .AND. NSUPERDIA>1)THEN
            IF(NBPROCDIA(KLOOP)>1 .OR. NBLVLKDIA(KLOOP,1)>1 &
	    .OR. NBNDIA(KLOOP)>1)THEN
            print *,' _PVKT_  SUPERPOSITIONS : '
            print *,'         On ne peut definir de part de d''autre '&
           &'de _ON_ qu''1 seul processus, 1 seul niveau, 1 seul masque'
            print *,' Nb de niveaux demandes   : ',NBLVLKDIA(KLOOP,1)
            print *,' Nb de processus demandes : ',NBPROCDIA(KLOOP)
            print *,' Nb de masques demandes   : ',NBNDIA(KLOOP)
            print *,' *** MODIFIEZ VOTRE DIRECTIVE *** '
            EXIT
          ENDIF
          ENDIF
        IK=SIZE(XVAR,3)
        ALLOCATE(ZTEM1D(IK),ZWORKZ(IK))
!
! Controle ordre des niveaux demandes. Eventuellement remise dans l'ordre 
! croissant
!
      INBK=NBLVLKDIA(KLOOP,NLOOPN)
      NKH=INBK
      IF(INBK > 1)THEN
      DO JLOOPK=1,INBK-1
	INUMK=NLVLKDIA(JLOOPK,KLOOP,NLOOPN)
        DO JLOOPK1=JLOOPK+1,INBK
	  INUMK1=NLVLKDIA(JLOOPK1,KLOOP,NLOOPN)
	  IF(INUMK < INUMK1)THEN
	    CYCLE
	  ELSE
	    NLVLKDIA(JLOOPK,KLOOP,NLOOPN)=INUMK1
	    NLVLKDIA(JLOOPK1,KLOOP,NLOOPN)=INUMK
	  ENDIF
	ENDDO
      ENDDO
      ENDIF
!
! Altitudes enregistees du niv 1 a n dans l'ordre croissant   --> GINVZ=.FALSE.
! Altitudes enregistees du niv n a 1 dans l'ordre decroissant --> GINVZ=.TRUE.
!
      IF(XTRAJZ(NLVLKDIA(1,KLOOP,NLOOPN),1,NNDIA(JLOOPN,KLOOP)) <  &
        XTRAJZ(NLVLKDIA(INBK,KLOOP,NLOOPN),1,NNDIA(JLOOPN,KLOOP)))THEN
        GINVZ=.FALSE.
      ELSE
        GINVZ=.TRUE.
! Remise des niveaux dans un ordre tel que les altitudes soient croissantes
! (/indices croissants)
	  NLVLKDIA(1:INBK,KLOOP,NLOOPN)=NLVLKDIA(INBK:1:-1,KLOOP,NLOOPN)
      ENDIF

!
! ************************
! Boucle sur les processus
! ************************
!
        DO JLOOPP=1,NBPROCDIA(KLOOP)
	  NLOOPP=NPROCDIA(JLOOPP,KLOOP)

		  CALL LOADUNITIT(JLOOPP,KLOOP)

          ZTEM1D(:)=0.; ZWORKZ(:)=0.

	  INDN=NNDIA(JLOOPN,KLOOP)

          IF(.NOT.LTINCRDIA(KLOOP,1))THEN      !----------------------- Tps
!
! Expression temps non incrementale
!
            DO JLOOPT=1,NBTIMEDIA(KLOOP,1)
              NLOOPT=NTIMEDIA(JLOOPT,KLOOP,1)
              CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,1))
              WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
	      ZTEM1D(1:IK)=XVAR(1,1,: &
	     ,NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP),NPROCDIA(JLOOPP,KLOOP))
              ZWORKZ(:)=XTRAJZ(:,1,NNDIA(JLOOPN,KLOOP))


    	      IF(LPV)THEN                                  !---LPV(KT)(1)-----
                IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
  	          IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		  ALLOCATE(XPRDAT(16,1))
		  CALL LOAD_XPRDAT(1,NLOOPT)
                ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
! PENSER A EXTRAIRE LES DIFFERENTS NIVEAUX

                ALLOCATE(ZTE(INBK),ZWO(INBK))

		  DO JLOOPK=1,NBLVLKDIA(KLOOP,NLOOPN)
		    ZTE(JLOOPK)=ZTEM1D(NLVLKDIA(JLOOPK,KLOOP,NLOOPN))
		    ZWO(JLOOPK)=ZWORKZ(NLVLKDIA(JLOOPK,KLOOP,NLOOPN))
		  ENDDO
	        DEALLOCATE(ZTEM1D,ZWORKZ)
	        ALLOCATE(ZTEM1D(SIZE(ZTE)))
	        ALLOCATE(ZWORKZ(SIZE(ZWO)))
	        ZTEM1D=ZTE; ZWORKZ=ZWO
	        DEALLOCATE(ZTE,ZWO)

                CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  DEALLOCATE(XPRDAT)
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

	        IF(KLOOP == NSUPERDIA)CALL FRAME

    	      ELSE IF(LPVT .OR. LPVKT .OR. LPVKT1)THEN     !---LPV(KT)(1)-----

                IF(JLOOPT == 1)THEN
		  ILENW=NBTIMEDIA(KLOOP,1)
                  ALLOCATE(ZTEM2D(1:IK,ILENW))
	          ZTEM2D=XSPVAL
	          ALLOCATE(ZWORKT(ILENW))
                  IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
  	            IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		    ALLOCATE(XPRDAT(16,ILENW))
                  ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		ENDIF

		IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
		ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	        ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,1),1)
	        ZTEM2D(1:IK,JLOOPT)= &
                XVAR(1,1,:,  &
                NTIMEDIA(JLOOPT,KLOOP,1),NNDIA(JLOOPN,KLOOP), &
		NPROCDIA(JLOOPP,KLOOP))

	        IF(JLOOPT == NBTIMEDIA(KLOOP,1))THEN
! PENSER A EXTRAIRE LES DIFFERENTS NIVEAUX

		  IF(ALLOCATED(XZSOL))THEN
		    DEALLOCATE(XZSOL)
                  ENDIF
                  ALLOCATE(ZTE2(INBK,ILENW),XZSOL(INBK))

		    DO JLOOPK=1,NBLVLKDIA(KLOOP,NLOOPN)
		      ZTE2(JLOOPK,:)=ZTEM2D(NLVLKDIA(JLOOPK,KLOOP,NLOOPN),:)
		      XZSOL(JLOOPK)=XTRAJZ(NLVLKDIA(JLOOPK,KLOOP,NLOOPN),1,  &
		      NNDIA(JLOOPN,KLOOP))
		    ENDDO

	          DEALLOCATE(ZTEM2D)
	          ALLOCATE(ZTEM2D(SIZE(ZTE2,1),SIZE(ZTE2,2)))
	          ZTEM2D=ZTE2
	          DEALLOCATE(ZTE2)
		  
		  XPVMIN=MINVAL(ZTEM2D)
	          XPVMAX=MAXVAL(ZTEM2D)
                  CALL VALMNMX(XPVMIN,XPVMAX)

                  IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
                    XPVMIN=XPVMIN-1.
                    XPVMAX=XPVMAX+1.
                  ENDIF

		  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
		  DEALLOCATE(ZTEM2D,ZWORKT)
	          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		    DEALLOCATE(XPRDAT)
                  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		  IF(.NOT.LPBREAD)THEN
	            IF(KLOOP == NSUPERDIA)CALL FRAME
	          ENDIF
	        ENDIF

    	      ENDIF                                        !---LPV(KT)(1)-----

            ENDDO  ! Fin Boucle Temps (Non incremental)

          ELSE                           !----------------------- Tps

! Expression temps incrementale

            DO JLOOPT=NTIMEDIA(1,KLOOP,1),NTIMEDIA(2,KLOOP,1),NTIMEDIA(3,KLOOP,1)
	      NLOOPT=JLOOPT
              CALL RESOLV_TIMES(JLOOPT)
              WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,1)
	      ZTEM1D(1:IK)=XVAR(1,1,: &
	      ,JLOOPT,INDN,NPROCDIA(JLOOPP,KLOOP))
	      ZWORKZ(:)=XTRAJZ(:,1,INDN)

    	      IF(LPV)THEN                                  !---LPV(KT)(1)-----
                IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
  	          IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
                  ALLOCATE(XPRDAT(16,1))
		  CALL LOAD_XPRDAT(1,NLOOPT)
                ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
! PENSER A EXTRAIRE LES DIFFERENTS NIVEAUX

                ALLOCATE(ZTE(INBK),ZWO(INBK))

		  DO JLOOPK=1,NBLVLKDIA(KLOOP,NLOOPN)
		    ZTE(JLOOPK)=ZTEM1D(NLVLKDIA(JLOOPK,KLOOP,NLOOPN))
		    ZWO(JLOOPK)=ZWORKZ(NLVLKDIA(JLOOPK,KLOOP,NLOOPN))
		  ENDDO

	        DEALLOCATE(ZTEM1D,ZWORKZ)
	        ALLOCATE(ZTEM1D(SIZE(ZTE)))
	        ALLOCATE(ZWORKZ(SIZE(ZWO)))
	        ZTEM1D=ZTE; ZWORKZ=ZWO
	        DEALLOCATE(ZTE,ZWO)

	        CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
	        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  DEALLOCATE(XPRDAT)
                ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	        IF(KLOOP == NSUPERDIA)CALL FRAME

    	      ELSE IF(LPVT .OR. LPVKT .OR.LPVKT1)THEN      !---LPV(KT)(1)-----

		IF(JLOOPT == NTIMEDIA(1,KLOOP,1))THEN
		  ILENW=(NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/NTIMEDIA(3,KLOOP,1) +1
                  IF(NVERBIA > 0)THEN
                  print *,' OPER  NTIMEDIA(2,KLOOP,1) NTIMEDIA(1,KLOOP,1) NTIMEDIA(3,KLOOP,1) ILENW ', &
                  NTIMEDIA(2,KLOOP,1),NTIMEDIA(1,KLOOP,1),NTIMEDIA(3,KLOOP,1), &
		  ILENW, &
                  XTIMEDIA(2,KLOOP,1),XTIMEDIA(1,KLOOP,1),XTIMEDIA(3,KLOOP,1)
		  ENDIF
		  ITIMEND=NTIMEDIA(1,KLOOP,1) + &
		  (((NTIMEDIA(2,KLOOP,1)-NTIMEDIA(1,KLOOP,1))/  &
		  NTIMEDIA(3,KLOOP,1))*NTIMEDIA(3,KLOOP,1))
		  if(nverbia > 0)then
		  print *,' ITIMEND B ',ITIMEND
		  endif
		  IF(ALLOCATED(ZTEM2D))THEN
		    DEALLOCATE(ZTEM2D)
		  ENDIF
	          ALLOCATE(ZTEM2D(1:IK,ILENW))
		  ZTEM2D=XSPVAL
		  IF(ALLOCATED(ZWORKT))THEN
		    DEALLOCATE(ZWORKT)
		  ENDIF
                  ALLOCATE(ZWORKT(ILENW))
		  IJLT=0
                  IF(LPRDAT)THEN !  Juin 2001 Ajout des dates ds FICVAL 
  	            IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
		    ALLOCATE(XPRDAT(16,ILENW))
                  ENDIF !  Juin 2001 Ajout des dates ds FICVAL 
		ENDIF

	        IJLT=IJLT+1
		IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  CALL LOAD_XPRDAT(IJLT,NLOOPT)
		ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                ZWORKT(IJLT)=XTRAJT(JLOOPT,1)
		ZTEM2D(1:IK,IJLT)= &
		XVAR(1,1,:,  &
                JLOOPT,INDN,NPROCDIA(JLOOPP,KLOOP))

!               IF(JLOOPT == NTIMEDIA(2,KLOOP,1))THEN
	        IF(JLOOPT == ITIMEND)THEN
! PENSER A EXTRAIRE LES DIFFERENTS NIVEAUX

		  IF(ALLOCATED(XZSOL))THEN
		    DEALLOCATE(XZSOL)
                  ENDIF
		  IF(ALLOCATED(ZTE2))THEN
		    DEALLOCATE(ZTE2)
                  ENDIF
                  ALLOCATE(ZTE2(INBK,ILENW),XZSOL(INBK))

		    DO JLOOPK=1,NBLVLKDIA(KLOOP,NLOOPN)
		      ZTE2(JLOOPK,:)=ZTEM2D(NLVLKDIA(JLOOPK,KLOOP,NLOOPN),:)
		      XZSOL(JLOOPK)=XTRAJZ(NLVLKDIA(JLOOPK,KLOOP,NLOOPN),1,  &
		      NNDIA(JLOOPN,KLOOP))
		    ENDDO
	                    
	          DEALLOCATE(ZTEM2D)
	          ALLOCATE(ZTEM2D(SIZE(ZTE2,1),SIZE(ZTE2,2)))
	          ZTEM2D=ZTE2
	          DEALLOCATE(ZTE2)
		  XPVMIN=MINVAL(ZTEM2D)
		  XPVMAX=MAXVAL(ZTEM2D)
	          CALL VALMNMX(XPVMIN,XPVMAX)

                  IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
                    XPVMIN=XPVMIN-1.
                    XPVMAX=XPVMAX+1.
                  ENDIF

                  CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
                  DEALLOCATE(ZTEM2D,ZWORKT)
	          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		    DEALLOCATE(XPRDAT)
                  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		  IF(.NOT.LPBREAD)THEN
	            IF(KLOOP == NSUPERDIA)CALL FRAME
		  ENDIF

                ENDIF

    	      ENDIF                                        !---LPV(KT)(1)-----

	    ENDDO  ! Fin Boucle Temps (Incremental)

          ENDIF                        !----------------------- Tps

	ENDDO  !  Fin Boucle Processus

	DEALLOCATE(ZTEM1D,ZWORKZ)

      ENDDO  !  Fin Boucle Num station



!*****************************************************************************
!*****************************************************************************
    CASE('SPXY')

      if(nverbia > 0)then
        print *,' **oper AV SUBSPXY '
      ENDIF
      CALL SUBSPXY(KLOOP)
      if(nverbia > 0)then
        print *,' **oper AP SUBSPXY '
      ENDIF
!*****************************************************************************
!*****************************************************************************
    CASE('DRST','RAPL')

      L1DT=.TRUE.
      DO JLOOPN=1,NBNDIA(KLOOP)
	
	NLOOPN=NNDIA(JLOOPN,KLOOP)

! Controle ordre des niveaux demandes. Eventuellement remise dans l'ordre 
! croissant pour verifier si les altitudes sont en ordre croissant ou
! decroissant (/aux indices croissant)
!
        INBK=NBLVLKDIA(KLOOP,NLOOPN)
	NKH=INBK
        IF(INBK > 1)THEN
        DO JLOOPK=1,INBK-1
  	INUMK=NLVLKDIA(JLOOPK,KLOOP,NLOOPN)
          DO JLOOPK1=JLOOPK+1,INBK
  	  INUMK1=NLVLKDIA(JLOOPK1,KLOOP,NLOOPN)
  	  IF(INUMK < INUMK1)THEN
  	    CYCLE
  	  ELSE
  	    NLVLKDIA(JLOOPK,KLOOP,NLOOPN)=INUMK1
  	    NLVLKDIA(JLOOPK1,KLOOP,NLOOPN)=INUMK
  	  ENDIF
  	  ENDDO
        ENDDO
        ENDIF
  !
  ! Altitudes enregistees du niv 1 a n dans l'ordre croissant   --> GINVZ=.FALSE.
  ! Altitudes enregistees du niv n a 1 dans l'ordre decroissant --> GINVZ=.TRUE.
  !
        IF(XTRAJZ(NLVLKDIA(1,KLOOP,NLOOPN),1,NNDIA(JLOOPN,KLOOP)) <  &
          XTRAJZ(NLVLKDIA(INBK,KLOOP,NLOOPN),1,NNDIA(JLOOPN,KLOOP)))THEN
          GINVZ=.FALSE.
        ELSE
          GINVZ=.TRUE.
! Remise des niveaux dans un ordre tel que les altitudes soient croissantes
	  NLVLKDIA(1:INBK,KLOOP,NLOOPN)=NLVLKDIA(INBK:1:-1,KLOOP,NLOOPN)
        ENDIF
  
!

	  IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN
	    ILENW=NBTIMEDIA(KLOOP,NLOOPN)
	  ELSE
	    ILENW=(NTIMEDIA(2,KLOOP,NLOOPN)-NTIMEDIA(1,KLOOP,NLOOPN))/ &
		   NTIMEDIA(3,KLOOP,NLOOPN)+1
	  ENDIF
          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	    IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
            ALLOCATE(XPRDAT(16,ILENW))
	  ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

        IF(LZT .OR. LPV .OR. LPVT .OR. LPVKT .OR. LPVKT1)THEN


	  IF(LZT)THEN
	    LPVKT1=.TRUE.
	    JLOOPPF=1
          ELSE
	    JLOOPPF=NBPROCDIA(KLOOP)
          ENDIF

! Boucle sur les processus

	  DO JLOOPP = 1,JLOOPPF
            NLOOPP=NPROCDIA(JLOOPP,KLOOP)
            if(nverbia >0)then
              print *, '***OPEROPER NLOOPP,JLOOPPF ', NLOOPP,JLOOPPF
            endif

	    CALL LATLONGRID

	    IK=NBLVLKDIA(KLOOP,NLOOPN)
	    ALLOCATE (ZTEM2D(1:IK,ILENW),ZWORKT(ILENW),ZWORKZ(IK))
	    IJLT=0

            IF(LZT)THEN
    	      CTITGAL='Altitude'
	      CUNITGAL='(M)'
    	    ELSE
              CTITGAL=ADJUSTL(CTITRE(NPROCDIA(JLOOPP,KLOOP)))
              CUNITGAL=ADJUSTL(CUNITE(NPROCDIA(JLOOPP,KLOOP)))
            ENDIF
            CTITGAL=ADJUSTL(CTITGAL)
            CUNITGAL=ADJUSTL(CUNITGAL)
		  CUNITGAL(INDEX(CUNITGAL,' '):LEN(CUNITGAL))=' '
    !

            IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN
  
    	      DO JLOOPT=1,NBTIMEDIA(KLOOP,NLOOPN)
	        NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
                IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	          NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
	          CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
	        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                CALL RESOLV_TIMES(NTIMEDIA(JLOOPT,KLOOP,NLOOPN))
    	        ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,NLOOPN), &
    				    NLOOPN)
		WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(NTIMEDIA(JLOOPT,KLOOP,NLOOPN), &
		NLOOPN)
    		DO JLOOPK=1,NBLVLKDIA(KLOOP,NLOOPN)

    		  IF(LZT)THEN

    		    ZTEM2D(JLOOPK,JLOOPT)=XTRAJZ(NLVLKDIA(JLOOPK,KLOOP, &
    		           NLOOPN),NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
                  ELSE
    		    ZTEM2D(JLOOPK,JLOOPT)=XVAR(1,1,NLVLKDIA(JLOOPK, &
		    KLOOP, &
    		    NLOOPN),NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN,NLOOPP)
                    if(nverbia > 0)then
                      print *,' **OPER modif JLOOPP en NLOOPPP '
                    endif
!ERRJD    NLOOPN),NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN,JLOOPP)
		    IF(LPV)THEN
		      ZWORKZ(JLOOPK)=XTRAJZ(NLVLKDIA(JLOOPK,KLOOP,NLOOPN), &
				     NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
		    ENDIF

    		  ENDIF

    		ENDDO

		IF(LPV)THEN
		  ALLOCATE(ZTEM1D(SIZE(ZTEM2D,1)))
		  ZTEM1D(:)=ZTEM2D(:,JLOOPT)
		  CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
		  IF(KLOOP == NSUPERDIA)CALL FRAME
		  DEALLOCATE(ZTEM1D)
	        ENDIF

    	      ENDDO
    
    	    ELSE
  
  	      DO JLOOPT=NTIMEDIA(1,KLOOP,NLOOPN),NTIMEDIA(2,KLOOP, &
			NLOOPN),NTIMEDIA(3,KLOOP,NLOOPN)
                NLOOPT=JLOOPT

                CALL RESOLV_TIMES(JLOOPT)
  	        IJLT=IJLT+1
                IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	          CALL LOAD_XPRDAT(IJLT,JLOOPT)
	        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
  	        ZWORKT(IJLT)=XTRAJT(JLOOPT,NLOOPN)
		WRITE(CTIMEC(8:15),'(F8.0)')XTRAJT(JLOOPT,NLOOPN)

  		DO JLOOPK=1,NBLVLKDIA(KLOOP,NLOOPN)

  		  IF(LZT)THEN
  		    ZTEM2D(JLOOPK,IJLT)=XTRAJZ(NLVLKDIA(JLOOPK,KLOOP, &
  		           NLOOPN),JLOOPT,NLOOPN)
                  ELSE
		    ZTEM2D(JLOOPK,IJLT)=XVAR(1,1,NLVLKDIA(JLOOPK,KLOOP,&
		    NLOOPN),JLOOPT,NLOOPN,NLOOPP)
!ERRJD              NLOOPN),JLOOPT,NLOOPN,JLOOPP)
		    IF(LPV)THEN
		      ZWORKZ(JLOOPK)=XTRAJZ(NLVLKDIA(JLOOPK,KLOOP,NLOOPN), &
				     JLOOPT,NLOOPN)
		    ENDIF
  		  ENDIF

  		ENDDO

		IF(LPV)THEN
		  ALLOCATE(ZTEM1D(SIZE(ZTEM2D,1)))
		  ZTEM1D(:)=ZTEM2D(:,IJLT)
		  CALL TRAPRO_FORDIACHRO(ZTEM1D,ZWORKZ,KLOOP)
		  IF(KLOOP == NSUPERDIA)CALL FRAME
		  DEALLOCATE(ZTEM1D)
	        ENDIF

  	      ENDDO
  
            ENDIF
  
  	    XPVMIN=MINVAL(ZTEM2D)
  	    XPVMAX=MAXVAL(ZTEM2D)
            CALL VALMNMX(XPVMIN,XPVMAX)
  
            IF(ABS(XPVMAX-XPVMIN) < 1.E-4)THEN
              XPVMIN=XPVMIN-1.
              XPVMAX=XPVMAX+1.
            ENDIF
  
	    IF(.NOT.LPV)THEN
              CALL PVFCT(ZWORKT,ZTEM2D,KLOOP)
              DEALLOCATE(ZTEM2D,ZWORKT)
	      IF(.NOT.LPBREAD)THEN
  	        IF(KLOOP == NSUPERDIA)CALL FRAME
	      ENDIF
	    ENDIF
	     
            IF(ALLOCATED(ZTEM2D))THEN
	      DEALLOCATE(ZTEM2D)
            ENDIF
            IF(ALLOCATED(ZWORKT))THEN
	      DEALLOCATE(ZWORKT)
            ENDIF
            IF(ALLOCATED(ZWORKZ))THEN
	      DEALLOCATE(ZWORKZ)
            ENDIF

	  ENDDO  ! Fin Boucle processus

        ELSE IF(LZTPVKT1)THEN
        ELSE IF(LXT .OR. LYT .OR. LXYDIA)THEN

          ALLOCATE(ZWORKT(ILENW),ZWORKY(ILENW))
          YTITX(1:LEN(YTITX))=' '
          YTITY(1:LEN(YTITY))=' '
          IJLT=0
          ILOOPP=NLOOPP
          NLOOPP=1
          CALL LATLONGRID
          NLOOPP=ILOOPP

          IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN

            DO JLOOPT=1,NBTIMEDIA(KLOOP,NLOOPN)
!! Octobre 2001
	      NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
!! Octobre 2001
              IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	        NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
	        CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
	      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

              IF(LXT .OR. LYT)THEN
                ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
                YTITX='TIME (sec)'
              ELSE IF(LXYDIA)THEN
                ZWORKT(JLOOPT)=XTRAJX(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)

		CALL CONV2XY(ZWORKT(JLOOPT), &
		XTRAJY(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN),ZX,ZY,11)

                YTITX='X'
              ENDIF
              YTITX=ADJUSTL(YTITX)
      
              IF(LXT)THEN
                ZWORKY(JLOOPT)=XTRAJX(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)

		CALL CONV2XY(ZWORKY(JLOOPT), &
		XTRAJY(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN),ZX,ZY,11)

                YTITY='X'
              ELSE IF(LXYDIA .OR. LYT)THEN
                ZWORKY(JLOOPT)=XTRAJY(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)

		CALL CONV2XY(XTRAJX(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN), &
		ZWORKY(JLOOPT),ZX,ZY,22)
!               IF(LCONV2XY .AND. NLATLON /= 0)THEN
!                 CALL SM_XYHAT_S(XLATORI,XLONORI, &
!                 XTRAJX(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN),ZWORKY(JLOOPT),ZX,ZY)
!                 ZWORKY(JLOOPT)=ZY
!               ENDIF
                YTITY='Y'
              ENDIF
              YTITY=ADJUSTL(YTITY)

            ENDDO
            
            ZTIMED=XTRAJT(NTIMEDIA(1,KLOOP,NLOOPN),NLOOPN)
            ZTIMEF=XTRAJT(NTIMEDIA(NBTIMEDIA(KLOOP,NLOOPN),KLOOP,NLOOPN),NLOOPN)

          ELSE

            DO JLOOPT=NTIMEDIA(1,KLOOP,NLOOPN),NTIMEDIA(2,KLOOP, &
                        NLOOPN),NTIMEDIA(3,KLOOP,NLOOPN)
!! Octobre 2001
               NLOOPT=JLOOPT
!! Octobre 2001

              IJLT=IJLT+1
            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	      CALL LOAD_XPRDAT(IJLT,JLOOPT)
	    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
              IF(LXT .OR. LYT)THEN
                ZWORKT(IJLT)=XTRAJT(JLOOPT,NLOOPN)
                YTITX='TIME (sec)'
              ELSE IF(LXYDIA)THEN
                ZWORKT(IJLT)=XTRAJX(1,JLOOPT,NLOOPN)
		CALL CONV2XY(ZWORKT(IJLT), &
		XTRAJY(1,JLOOPT,NLOOPN),ZX,ZY,11)
!               IF(LCONV2XY .AND. NLATLON /= 0)THEN
!                 CALL SM_XYHAT_S(XLATORI,XLONORI,ZWORKT(IJLT), &
!                 XTRAJY(1,JLOOPT,NLOOPN),ZX,ZY)
!                 ZWORKT(IJLT)=ZX
!               ENDIF
                YTITX='X'
              ENDIF
      
              IF(LXT)THEN
                ZWORKY(IJLT)=XTRAJX(1,JLOOPT,NLOOPN)
		CALL CONV2XY(ZWORKY(IJLT), &
		XTRAJY(1,JLOOPT,NLOOPN),ZX,ZY,11)
!               IF(LCONV2XY .AND. NLATLON /= 0)THEN
!                 CALL SM_XYHAT_S(XLATORI,XLONORI,ZWORKY(IJLT), &
!                 XTRAJY(1,JLOOPT,NLOOPN),ZX,ZY)
!                 ZWORKY(IJLT)=ZX
!               ENDIF
                YTITY='X'
              ELSE IF(LXYDIA .OR. LYT)THEN
                ZWORKY(IJLT)=XTRAJY(1,JLOOPT,NLOOPN)
		CALL CONV2XY(XTRAJX(1,JLOOPT,NLOOPN), &
		ZWORKY(IJLT),ZX,ZY,22)
!               IF(LCONV2XY .AND. NLATLON /= 0)THEN
!               CALL SM_XYHAT_S(XLATORI,XLONORI, &
!                 XTRAJX(1,JLOOPT,NLOOPN),ZWORKY(IJLT),ZX,ZY)
!                 ZWORKY(IJLT)=ZY
!               ENDIF
                YTITY='Y'
              ENDIF

            ENDDO

            ZTIMED=XTRAJT(NTIMEDIA(1,KLOOP,NLOOPN),NLOOPN)
            ZTIMEF=XTRAJT(NTIMEDIA(2,KLOOP,NLOOPN),NLOOPN)

          ENDIF

          CALL TRAXY(ZWORKT,ZWORKY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)

	  DEALLOCATE(ZWORKT,ZWORKY)
	  IF(KLOOP == NSUPERDIA)THEN
	    IF(LDATFILE)CALL DATFILE_FORDIACHRO
	    CALL FRAME
          ENDIF

        ENDIF
        IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
          DEALLOCATE(XPRDAT)
        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

      ENDDO  ! Fin Boucle Numeros DRST

!*****************************************************************************
!*****************************************************************************
    CASE('RSPL')

      DO JLOOPN=1,NBNDIA(KLOOP)

	NLOOPN=NNDIA(JLOOPN,KLOOP)

! Traitement des RS
! *****************
	IF(LRS .OR. LRS1)THEN
!
! Cas LRS ou LRS1 et KLOOP = 1 --> Allocation de tableaux pour memoriser
! les infos utiles
!
! LRS : pas de superpositions ; donc KLOOP=NSUPERDIA=1 . Boucle externe sur le
! Num. des RS (que l'on peut ou non preciser dans les directives) . Boucle
! interne sur les temps (que l'on peut ou non preciser) avant appel TSOUND.
!
! LRS1 : superpositions ; KLOOP varie . De part et d'autre de _ON_ on ne
! donne qu'1 station . Donc JLOOPN tjrs = 1
!
          IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN
            ILENW=NBTIMEDIA(KLOOP,NLOOPN)
          ELSE
            ILENW=(NTIMEDIA(2,KLOOP,NLOOPN)-NTIMEDIA(1,KLOOP,NLOOPN))/ &
            NTIMEDIA(3,KLOOP,NLOOPN)+1
          ENDIF

          NST(KLOOP)=ILENW
          IF(KLOOP == 1)THEN
!
! SIZE(XVAR,3) = normalement 1
!
            ALLOCATE(XTRS(SIZE(XVAR,3)*NSUPERDIA,ILENW))
            ALLOCATE(XPRS(SIZE(XVAR,3)*NSUPERDIA,ILENW))
            ALLOCATE(XURS(SIZE(XVAR,3)*NSUPERDIA,ILENW))
            ALLOCATE(XVRS(SIZE(XVAR,3)*NSUPERDIA,ILENW))
            ALLOCATE(XRVRS(SIZE(XVAR,3)*NSUPERDIA,ILENW))
            ALLOCATE(XTIMRS2(SIZE(XVAR,3)*NSUPERDIA,ILENW))
	    ALLOCATE(NST(SIZE(XVAR,3)*NSUPERDIA))
	    ALLOCATE(NNST(SIZE(XVAR,3)*NSUPERDIA))
            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
	      ALLOCATE(XPRDAT(16,ILENW))
	    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

          ENDIF

	  IF(KLOOP > 1 .AND. LRS1)THEN
	    IF(ILENW > SIZE(XTRS,2))THEN
	      ALLOCATE(ZWORKRS(SIZE(XTRS,1),SIZE(XTRS,2)))
	      ZWORKRS(:,:)=XTRS(:,:)
	      DEALLOCATE(XTRS)
	      ALLOCATE(XTRS(SIZE(ZWORKRS,1),ILENW))
	      XTRS(:,1:SIZE(ZWORKRS,2))=ZWORKRS(:,:)
	      DEALLOCATE(ZWORKRS)
	      ALLOCATE(ZWORKRS(SIZE(XPRS,1),SIZE(XPRS,2)))
	      ZWORKRS(:,:)=XPRS(:,:)
	      DEALLOCATE(XPRS)
	      ALLOCATE(XPRS(SIZE(ZWORKRS,1),ILENW))
	      XPRS(:,1:SIZE(ZWORKRS,2))=ZWORKRS(:,:)
	      DEALLOCATE(ZWORKRS)
	      ALLOCATE(ZWORKRS(SIZE(XURS,1),SIZE(XURS,2)))
	      ZWORKRS(:,:)=XURS(:,:)
	      DEALLOCATE(XURS)
	      ALLOCATE(XURS(SIZE(ZWORKRS,1),ILENW))
	      XURS(:,1:SIZE(ZWORKRS,2))=ZWORKRS(:,:)
	      DEALLOCATE(ZWORKRS)
	      ALLOCATE(ZWORKRS(SIZE(XVRS,1),SIZE(XVRS,2)))
	      ZWORKRS(:,:)=XVRS(:,:)
	      DEALLOCATE(XVRS)
	      ALLOCATE(XVRS(SIZE(ZWORKRS,1),ILENW))
	      XVRS(:,1:SIZE(ZWORKRS,2))=ZWORKRS(:,:)
	      DEALLOCATE(ZWORKRS)
	      ALLOCATE(ZWORKRS(SIZE(XRVRS,1),SIZE(XRVRS,2)))
	      ZWORKRS(:,:)=XRVRS(:,:)
	      DEALLOCATE(XRVRS)
	      ALLOCATE(XRVRS(SIZE(ZWORKRS,1),ILENW))
	      XRVRS(:,1:SIZE(ZWORKRS,2))=ZWORKRS(:,:)
	      DEALLOCATE(ZWORKRS)
	      ALLOCATE(ZWORKRS(SIZE(XTIMRS2,1),SIZE(XTIMRS2,2)))
	      ZWORKRS(:,:)=XTIMRS2(:,:)
	      DEALLOCATE(XTIMRS2)
	      ALLOCATE(XTIMRS2(SIZE(ZWORKRS,1),ILENW))
	      XTIMRS2(:,1:SIZE(ZWORKRS,2))=ZWORKRS(:,:)
	      DEALLOCATE(ZWORKRS)
	    ENDIF
	  ENDIF

	  NNST(KLOOP)=NLOOPN

! Dans XVAR PROC1=TCelsius  PROC2=PRES(Pls) PROC3=U PROC4=V PROC5=RCM

	  IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN

  	    DO JLOOPT=1,NBTIMEDIA(KLOOP,NLOOPN)

              NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
              IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
	      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

  	      XTRS(KLOOP,JLOOPT)=XVAR(1,1,1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN), &
				     NLOOPN,1)+XTT
  	      XPRS(KLOOP,JLOOPT)=XVAR(1,1,1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN), &
								     NLOOPN,2)
	      XURS(KLOOP,JLOOPT)=XVAR(1,1,1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN), &
								     NLOOPN,3)
	      XVRS(KLOOP,JLOOPT)=XVAR(1,1,1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN), &
							             NLOOPN,4)
	      XRVRS(KLOOP,JLOOPT)=XVAR(1,1,1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),&
								     NLOOPN,5)
	      XTIMRS2(KLOOP,JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
            ENDDO

	  ELSE

 	    II=0
 	    DO JLOOPT=NTIMEDIA(1,KLOOP,NLOOPN),NTIMEDIA(2,KLOOP,NLOOPN),NTIMEDIA(3,KLOOP,NLOOPN)
!! Octobre 2001
              NLOOPT=JLOOPT
!! Octobre 2001
	      II=II+1
              IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		CALL LOAD_XPRDAT(II,JLOOPT)
	      ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	      XTRS(KLOOP,II)=XVAR(1,1,1,JLOOPT,NLOOPN,1)+273.16
 	      XPRS(KLOOP,II)=XVAR(1,1,1,JLOOPT,NLOOPN,2)
	      XTIMRS2(KLOOP,II)=XTRAJT(JLOOPT,NLOOPN)
 	      XURS(KLOOP,II)=XVAR(1,1,1,JLOOPT,NLOOPN,3)
 	      XVRS(KLOOP,II)=XVAR(1,1,1,JLOOPT,NLOOPN,4)
 	      XRVRS(KLOOP,II)=XVAR(1,1,1,JLOOPT,NLOOPN,5)
 	    ENDDO

          ENDIF

          GMXRAT=.TRUE.

	  CLEGEND(104:106)='U-V'
!         YTEXTE(1:5)='U-V'
          WRITE(YTEXTE,'(''I='',I2,'' J='',I2)')NIRS,NJRS
	  CALL TABCOL_FORDIACHRO
	  CALL GSTXFP(-13,2)

          IF(LRS)THEN

            IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN
              IF(NVERBIA > 0)THEN
              print *,' KLOOP,LRS,JLOOPT,NTIMEDIA(1,KLOOP,NLOOPN) ', &
              KLOOP,LRS,JLOOPT,NTIMEDIA(1,KLOOP,NLOOPN)
	      ENDIF
              CALL RESOLV_TIMES(NTIMEDIA(1,KLOOP,NLOOPN))
            ELSE
              II=NTIMEDIA(1,KLOOP,NLOOPN)
              CALL RESOLV_TIMES(II)
            ENDIF
! CTIMEC(S) est determine ds OPER pour LRS et ds TSOUND pour LRS1
              CTIMECS(1:LEN(CTIMECS))=' '
              CTIMECS(1:3)='  ('
              WRITE(CTIMECS(4:11),'(F8.0)')XTIMRS2(1,1)
              CTIMECS(LEN_TRIM(CTIMECS)+1:LEN_TRIM(CTIMECS)+1)='-'
	      YTEM(1:LEN(YTEM))=' '
	      WRITE(YTEM(1:8),'(F8.0)')XTIMRS2(1,ILENW)
	      YTEM=ADJUSTL(YTEM)
	      IN=LEN_TRIM(CTIMECS)
	      II=LEN_TRIM(YTEM)
	      IN=IN+1
	      CTIMECS(IN:IN+II-1)=YTEM(1:II)
	      IN=IN+1
	      CTIMECS(IN:IN+1)='s)'

            GMXRAT=.TRUE.
   
	    DO J=1,SIZE(XRVRS,2)
	      IF(XRVRS(1,J) <=0.)print *,' No dew point line drawn as nil or' &
		,' negative water values were found'
	    ENDDO

	    CALL GSCLIP(0)
	    CALL TSOUND_FORDIACHRO(XPRS(1,:),XTRS(1,:),  &
            XRVRS(1,:),XURS(1,:), &
            XVRS(1,:),SIZE(XPRS,2),CLEGEND,YTEXTE,GMXRAT,.TRUE.&
	    ,.FALSE.,.FALSE.)
	    CALL GSCLIP(1)
	    CALL FRAME

            DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS2,NST,NNST)

          ELSE IF(LRS1 .AND. KLOOP == NSUPERDIA)THEN

            GMXRAT=.TRUE.
! On met la date courante du 1er temps demande de la 1ere superposition
            CALL RESOLV_TIMES(NTIMEDIA(1,1,NLOOPN))
	    CALL GSCLIP(0)
! Dans OPER on ne transmet que le 1er temps et les autres son charges dans
! TSOUND
	    CALL TSOUND_FORDIACHRO(XPRS(1,:),XTRS(1,:),  &
		    XRVRS(1,:),XURS(1,:), &
		    XVRS(1,:),NST(1),CLEGEND,YTEXTE,GMXRAT,.TRUE.&
		    ,.FALSE.,.FALSE.)
	    CALL GSCLIP(1)
	    CALL FRAME
            DEALLOCATE(XTRS,XPRS,XURS,XVRS,XRVRS,XTIMRS2,NST,NNST)
            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	      DEALLOCATE(XPRDAT)
	    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

          ENDIF

! Infos autres que RS
! *******************
 
        ELSE

          IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN
            ILENW=NBTIMEDIA(KLOOP,NLOOPN)
          ELSE
	    ILENW=(NTIMEDIA(2,KLOOP,NLOOPN)-NTIMEDIA(1,KLOOP,NLOOPN))/ &
				       NTIMEDIA(3,KLOOP,NLOOPN)+1
	  ENDIF
            IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	      IF(ALLOCATED(XPRDAT))DEALLOCATE(XPRDAT)
	      ALLOCATE(XPRDAT(16,ILENW))
	    ENDIF ! Juin 2001 Ajout des dates ds FICVAL 

	  IF(LFT .OR. LFT1)THEN

            ALLOCATE(ZWORKT(ILENW),ZWORK1D(ILENW))

	    DO JLOOPP=1,NBPROCDIA(KLOOP)
	      NLOOPP=NPROCDIA(JLOOPP,KLOOP)

	      CALL LATLONGRID

              CTITGAL=ADJUSTL(CTITRE(NPROCDIA(JLOOPP,KLOOP)))
	      CUNITGAL=ADJUSTL(CUNITE(NPROCDIA(JLOOPP,KLOOP)))
	      CTITGAL=ADJUSTL(CTITGAL)
	      CUNITGAL=ADJUSTL(CUNITGAL)
		  CUNITGAL(INDEX(CUNITGAL,' '):LEN(CUNITGAL))=' '

	      IF(.NOT. LTINCRDIA(KLOOP,NLOOPN))THEN

		DO JLOOPT=1,NBTIMEDIA(KLOOP,NLOOPN)
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		    CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
	          ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		  IF(JLOOPT == 1)CALL RESOLV_TIMES(NLOOPT)
		  ZWORKT(JLOOPT)=XTRAJT(NLOOPT,NLOOPN)
		  ZWORK1D(JLOOPT)=XVAR(1,1,1,NLOOPT,NLOOPN,NPROCDIA(JLOOPP, &
								    KLOOP))
		ENDDO
	      ELSE
		IJLT=0
		DO JLOOPT=NTIMEDIA(1,KLOOP,NLOOPN),NTIMEDIA(2,KLOOP,NLOOPN), &
			  NTIMEDIA(3,KLOOP,NLOOPN)
                  NLOOPT=JLOOPT
		  IJLT=IJLT+1
                  IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		    CALL LOAD_XPRDAT(IJLT,NLOOPT)
	          ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
		  IF(IJLT == 1)CALL RESOLV_TIMES(NLOOPT)
		  ZWORKT(IJLT)=XTRAJT(NLOOPT,NLOOPN)
		  ZWORK1D(IJLT)=XVAR(1,1,1,NLOOPT,NLOOPN,NPROCDIA(JLOOPP, &
								    KLOOP))
		ENDDO

	      ENDIF

	      CALL VARFCT(ZWORKT,ZWORK1D,1)
	      if(nverbia > 0)then
		print *,' ** oper RSPL AP VARFCT KLOOP NSUPERDIA ',KLOOP,&
		NSUPERDIA
	      endif
	      IF(KLOOP == NSUPERDIA)CALL FRAME

	    ENDDO

	    DEALLOCATE(ZWORKT,ZWORK1D)

	  ELSE IF(LZT .OR. LXT .OR. LYT .OR. LXYDIA)THEN

            ALLOCATE(ZWORKT(ILENW),ZWORKY(ILENW))
            YTITX(1:LEN(YTITX))=' '
            YTITY(1:LEN(YTITY))=' '
            IJLT=0

            ILOOPP=NLOOPP
            NLOOPP=1
            CALL LATLONGRID
            NLOOPP=ILOOPP
  
            IF(.NOT.LTINCRDIA(KLOOP,NLOOPN))THEN
  
              DO JLOOPT=1,NBTIMEDIA(KLOOP,NLOOPN)
!! Octobre 2001
		NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
!! Octobre 2001
                IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  NLOOPT=NTIMEDIA(JLOOPT,KLOOP,NLOOPN)
		  CALL LOAD_XPRDAT(JLOOPT,NLOOPT)
	        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
  
                IF(LZT .OR. LXT .OR. LYT)THEN
                  ZWORKT(JLOOPT)=XTRAJT(NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
                  YTITX='TIME (sec)'
                ELSE IF(LXYDIA)THEN
                  ZWORKT(JLOOPT)=XTRAJX(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
		IF(LCONV2XY .AND. NLATLON /= 0)THEN
		  CALL SM_XYHAT_S(XLATORI,XLONORI,ZWORKT(JLOOPT), &
		  XTRAJY(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN),&
		  ZX,ZY)
                  ZWORKT(JLOOPT)=ZX
		ENDIF
                  YTITX='X'
                ENDIF
                YTITX=ADJUSTL(YTITX)
        
                IF(LZT)THEN
                  ZWORKY(JLOOPT)=XTRAJZ(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
                  YTITY='Z'
                ELSE IF(LXT)THEN
                  ZWORKY(JLOOPT)=XTRAJX(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
		IF(LCONV2XY .AND. NLATLON /= 0)THEN
		  CALL SM_XYHAT_S(XLATORI,XLONORI,ZWORKY(JLOOPT), &
		  XTRAJY(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN),&
		  ZX,ZY)
                  ZWORKY(JLOOPT)=ZX
		ENDIF
                  YTITY='X'
                ELSE IF(LXYDIA .OR. LYT)THEN
                  ZWORKY(JLOOPT)=XTRAJY(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN)
		IF(LCONV2XY .AND. NLATLON /= 0)THEN
		  CALL SM_XYHAT_S(XLATORI,XLONORI, &
		 XTRAJX(1,NTIMEDIA(JLOOPT,KLOOP,NLOOPN),NLOOPN),ZWORKY(JLOOPT),&
		  ZX,ZY)
                  ZWORKY(JLOOPT)=ZY
		ENDIF
                  YTITY='Y'
                ENDIF
                YTITY=ADJUSTL(YTITY)
  
              ENDDO
              
              ZTIMED=XTRAJT(NTIMEDIA(1,KLOOP,NLOOPN),NLOOPN)
              ZTIMEF=XTRAJT(NTIMEDIA(NBTIMEDIA(KLOOP,NLOOPN),KLOOP,NLOOPN),NLOOPN)
  
            ELSE

              DO JLOOPT=NTIMEDIA(1,KLOOP,NLOOPN),NTIMEDIA(2,KLOOP, &
                        NLOOPN),NTIMEDIA(3,KLOOP,NLOOPN)

!! Octobre 2001
                NLOOPT=JLOOPT
!! Octobre 2001
                IJLT=IJLT+1
                IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
		  CALL LOAD_XPRDAT(IJLT,JLOOPT)
	        ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
                IF(LZT .OR. LXT .OR. LYT)THEN
                  ZWORKT(IJLT)=XTRAJT(JLOOPT,NLOOPN)
                  YTITX='TIME (sec)'
                ELSE IF(LXYDIA)THEN
                  ZWORKT(IJLT)=XTRAJX(1,JLOOPT,NLOOPN)
		IF(LCONV2XY .AND. NLATLON /= 0)THEN
		  CALL SM_XYHAT_S(XLATORI,XLONORI,ZWORKT(IJLT), &
		  XTRAJY(1,JLOOPT,NLOOPN),ZX,ZY)
                  ZWORKT(IJLT)=ZX
		ENDIF
                  YTITX='X'
                ENDIF
		YTITX=ADJUSTL(YTITX)
      
                IF(LZT)THEN
                  ZWORKY(IJLT)=XTRAJZ(1,JLOOPT,NLOOPN)
                  YTITY='Z'
                ELSE IF(LXT)THEN
                  ZWORKY(IJLT)=XTRAJX(1,JLOOPT,NLOOPN)
		IF(LCONV2XY .AND. NLATLON /= 0)THEN
		  CALL SM_XYHAT_S(XLATORI,XLONORI,ZWORKY(IJLT), &
		  XTRAJY(1,JLOOPT,NLOOPN),ZX,ZY)
                  ZWORKY(IJLT)=ZX
		ENDIF
                  YTITY='X'
                ELSE IF(LXYDIA .OR. LYT)THEN
                  ZWORKY(IJLT)=XTRAJY(1,JLOOPT,NLOOPN)
		IF(LCONV2XY .AND. NLATLON /= 0)THEN
		  CALL SM_XYHAT_S(XLATORI,XLONORI, &
		  XTRAJX(1,JLOOPT,NLOOPN),ZWORKY(IJLT),ZX,ZY)
                  ZWORKY(IJLT)=ZY
		ENDIF
                  YTITY='Y'
                ENDIF
                YTITY=ADJUSTL(YTITY)

              ENDDO

              ZTIMED=XTRAJT(NTIMEDIA(1,KLOOP,NLOOPN),NLOOPN)
              ZTIMEF=XTRAJT(NTIMEDIA(2,KLOOP,NLOOPN),NLOOPN)

            ENDIF

	    CALL TRAXY(ZWORKT,ZWORKY,KLOOP,YTITX,YTITY,ZTIMED,ZTIMEF)
	      if(nverbia > 0)then
		print *,' ** oper RSPL AP TRAXY KLOOP NSUPERDIA ',KLOOP,&
		NSUPERDIA
	      endif

	    DEALLOCATE(ZWORKT,ZWORKY)
	    IF(KLOOP == NSUPERDIA)THEN
	      IF(LDATFILE)CALL DATFILE_FORDIACHRO
	      CALL FRAME
	    ENDIF

	  ENDIF

          IF(LPRDAT)THEN ! Juin 2001 Ajout des dates ds FICVAL 
	    DEALLOCATE(XPRDAT)
          ENDIF ! Juin 2001 Ajout des dates ds FICVAL 
	ENDIF

      ENDDO   ! Fin boucle N. RS ou avions

!*****************************************************************************
!*****************************************************************************
!   CASE('RAPL')


END SELECT

IF(ALLOCATED(ZWORK3D))THEN
  DEALLOCATE(ZWORK3D)
ENDIF
IF(ALLOCATED(ZWORK1D))THEN
  DEALLOCATE(ZWORK1D)
ENDIF
IF(ALLOCATED(ZWORKT))THEN
  DEALLOCATE(ZWORKT)
ENDIF
IF(ALLOCATED(ZWORKZ))THEN
  DEALLOCATE(ZWORKZ)
ENDIF
IF(ALLOCATED(ZWORKZ2))THEN
  DEALLOCATE(ZWORKZ2)
ENDIF
IF(ALLOCATED(ZWORKRS))THEN
  DEALLOCATE(ZWORKRS)
ENDIF
IF(ALLOCATED(ZWORKY))THEN
  DEALLOCATE(ZWORKY)
ENDIF
IF(ALLOCATED(ZTEMCV))THEN
  DEALLOCATE(ZTEMCV)
ENDIF
IF(ALLOCATED(ZTEM2D))THEN
  DEALLOCATE(ZTEM2D)
ENDIF
IF(ALLOCATED(ZTEM1D))THEN
  DEALLOCATE(ZTEM1D)
ENDIF
IF(ALLOCATED(ZTE))THEN
  DEALLOCATE(ZTE)
ENDIF
IF(ALLOCATED(ZTE2))THEN
  DEALLOCATE(ZTE2)
ENDIF
IF(ALLOCATED(ZWO))THEN
  DEALLOCATE(ZWO)
ENDIF
!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
if(nverbia > 0)then
  print *,' **oper sortie LPRESY,XHMIN,XHMAX ',LPRESY,XHMIN,XHMAX
endif
RETURN
END SUBROUTINE OPER_PROCESS
