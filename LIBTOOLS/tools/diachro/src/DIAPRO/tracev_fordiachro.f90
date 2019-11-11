!     ######spl
      MODULE MODI_TRACEV_FORDIACHRO
!     #############################
!
INTERFACE
!
SUBROUTINE TRACEV_FORDIACHRO(PTEMCV,KLOOP,HTEXT)
INTEGER          :: KLOOP
CHARACTER(LEN=*) :: HTEXT
REAL,DIMENSION(:,:)   :: PTEMCV
END SUBROUTINE TRACEV_FORDIACHRO
!
END INTERFACE
END MODULE MODI_TRACEV_FORDIACHRO
!     ######spl
      SUBROUTINE TRACEV_FORDIACHRO(PTEMCV,KLOOP,HTEXT)
!     ################################################
!
!!****  *TRACEV_FORDIACHRO* - Manager for the horizontal cross-section plots
!!
!!    PURPOSE
!!    -------
!       In the case of horizontal cross-sections, call the interpolation and
!     display routines: - for contour plots
!                       - for vector arrow plots
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!    VALNGRID   : retrieves the NGRID grid number when given the variable name
!!    COMCOORD   : computes true altitudes corresponding to the NGRID value
!!    IMCOU      : contour plot manager for vertical cross-sections
!!    IMCOUV     : vector  plot manager for vertical cross-sections
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_TITLE  : Declares heading variables for the plots (TRACE)
!!         NCONT  :  Current plot number
!!
!!      Module MODD_NMGRID : declares global variable  NMGRID (TRACE)
!!         NMGRID  : Current MESO-NH grid indicator
!!
!!      Module MODN_NCAR : defines NAM_DIRTRA_POS namelist
!!                         (former NCAR common)
!!         NHI        : Extrema detection
!!                       (=0 --> H+L, <0 nothing)
!!         NDOT       : Line style
!!                       (=0|1|1023|65535 --> solid lines;
!!                       <0 --> solid lines for positive values and
!!                       dotted lines(ABS(NDOT))for negative values;
!!                       >0 --> dotted lines(ABS(NDOT)) )
!!         XPHINT     : Increment contour value for PHIM, PHIT
!!
!!      Module MODD_OUT  : defines various logical units and dimensions
!!         NIMAXT     : x-size of the displayed section of the MESO-NH arrays
!!         NJMAXT     : y-size of the displayed section of the MESO-NH arrays
!!         NKMAXT     : z-size of the displayed section of the MESO-NH arrays
!!
!!      Module MODN_PARA
!!         Module MODD_DIM1 : contains dimensions of data arrays
!!             NKMAX :  z array dimensions  
!!
!!      Module MODD_PARAMETERS : Contains array border depths
!!         JPVEXT   : Vertical external points number
!!
!!      Module MODD_SUPER   : defines plot overlay control variables 
!!         LSUPER   : =.TRUE. --> plot overlay is active
!!                    =.FALSE. --> plot overlay is not active
!!         NSUPER   : Rank of the current plot in the overlay
!!                    sequence. The initial plot is rank 1.
!!
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TITLE
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO
USE MODD_NMGRID
USE MODN_NCAR
USE MODD_OUT
USE MODD_DIM1
USE MODN_PARA
USE MODD_SUPER
USE MODD_PT_FOR_CH_FORDIACHRO
!USE MODI_IMCOU_FORDIACHRO

IMPLICIT NONE
!
!*      0.1   Interfaces declaration
!
INTERFACE
      SUBROUTINE IMCOU_FORDIACHRO(PTABV,PINT,HLEGEND,HTEXT)
      CHARACTER(LEN=*)   :: HTEXT, HLEGEND
      REAL                :: PINT
      REAL,DIMENSION(:,:) :: PTABV
      END SUBROUTINE IMCOU_FORDIACHRO
END INTERFACE
INTERFACE
      SUBROUTINE IMCOUV_FORDIACHRO(PU,PW,HLEGEND,HTEXT)
      REAL,DIMENSION(:,:) :: PU,PW
      CHARACTER(LEN=*)    :: HTEXT
      CHARACTER(LEN=*)    :: HLEGEND
      END SUBROUTINE IMCOUV_FORDIACHRO
END INTERFACE
!
!*       0.15 Dummy arguments
!          
INTEGER          :: KLOOP, JU
CHARACTER(LEN=*) :: HTEXT
REAL,DIMENSION(:,:)   :: PTEMCV
!
!*       0.2  local variables
!          

REAL,SAVE         :: ZHMIN, ZHMAX
INTEGER          :: IBEGITXT, IENDTXT
INTEGER          :: ISTA, IER, INB, IWK, J, II

CHARACTER(LEN=20) :: YNOM

!
!-------------------------------------------------------------------------------
!*       1.3    Interactive option selection and plot overlay management  
!
!
!!!!!!!!!!! 110797
IF(NLOOPSUPER == 1)THEN
ZHMIN=XHMIN; ZHMAX=XHMAX
if(nverbia > 0)then
  print *,' TRACEV ENTREE XHMIN XHMAX ZHMIN ZHMAX ',XHMIN,XHMAX,ZHMIN,ZHMAX
endif
ELSE
  IF(NBPMT > 0)THEN
  DO J=NLOOPSUPER,1,-1
    IF(NUMPM(J) /= 0 .AND. NUMPM(J) /= 1 .AND. NUMPM(J) /= 2)THEN
      II=1
      EXIT
    ELSE
!     print *,' J NUMPM(J) ',J,NUMPM(J)
      II=0
    ENDIF
  ENDDO
  IF(II == 0)THEN
    ZHMIN=XHMIN; ZHMAX=XHMAX
  if(nverbia > 0)then
    print *,' TRACEV ENTREE XHMIN XHMAX ZHMIN ZHMAX ',XHMIN,XHMAX,ZHMIN,ZHMAX
  endif
  ENDIF
  ENDIF
ENDIF
!!!!!!!!!!! 110797

YNOM=ADJUSTL(CGROUP)
IF(YNOM.EQ.'QUIT')THEN
  CALL GQOPS(ISTA)
  CALL GQACWK(1,IER,INB,IWK)
  IF(ISTA >1 .AND. INB >1)THEN
    CALL GDAWK(2)
    CALL GCLWK(2)
  ENDIF
! CALL FRAME
  CALL NGPICT(1,1)
  CALL CLSGKS
  STOP
END IF
IBEGITXT=1
IENDTXT=30

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
XLWIDTH=XLWDEF
!
!*       2.     PROCESSING OF ALL VARIABLES
!               ---------------------------
!
IF(LULMWM .OR. LULTWT)THEN

  CALL IMCOUV_FORDIACHRO(PTEMCV,XWCV,CLEGEND,HTEXT)

! Ajout Janvier 2001
ELSE IF(LUMVM .OR. LUTVT .OR. LDIRWIND .OR. LSUMVM .OR. LSUTVT)THEN
  if(nverbia > 0)then
  print *,' LUMVM LDIRWIND LSUMVM AV CALL IMCOUV_FORDIACHRO ds TRACEV ',LUMVM,LDIRWIND,LSUMVM
  endif
  CALL IMCOUV_FORDIACHRO(PTEMCV,XWCV,CLEGEND,HTEXT)
 
ELSE IF((LDIRWM .AND. .NOT.LDIRWIND) .OR. (LDIRWT .AND. .NOT.LDIRWIND))THEN
  if(nverbia > 0)then
  print *,' tracev LDIRWM LDIRWT LDIRWIND AV call imcou_fordiachro ',LDIRWM,LDIRWT,LDIRWIND
  print *,' tracev SIZE(PTEMCV) AV call imcou_fordiachro ',SIZE(PTEMCV,1),SIZE(PTEMCV,2)
  endif
  CALL IMCOU_FORDIACHRO(PTEMCV,XDIAINT,CLEGEND,HTEXT)
ELSE

  CALL COMPCOORD_FORDIACHRO(NMGRID)
!print *,' ZSTAB'
!print *,ZSTAB
!REGLER LE PB DE L'INTERVALLE
  CALL IMCOU_FORDIACHRO(PTEMCV,XDIAINT,CLEGEND,HTEXT)
!
ENDIF
IF(ALLOCATED(XWCV))THEN
  DEALLOCATE(XWCV)
ENDIF

IF(.NOT.LSUPER .OR. (LSUPER .AND. NLOOPSUPER == NSUPERDIA))THEN
  XHMIN=ZHMIN; XHMAX=ZHMAX
if(nverbia > 0)then
  print *,' TRACEV SORTIE XHMIN XHMAX ZHMIN ZHMAX ',XHMIN,XHMAX,ZHMIN,ZHMAX
endif
ENDIF
!
  RETURN
!------------------------------------------------------------------------------
!
!*     5.       EXIT
!               ----
!
!
END SUBROUTINE  TRACEV_FORDIACHRO
