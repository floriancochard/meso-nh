!     ######spl
      MODULE  MODI_LOAD_SEGMENTS
!     #########################
!
INTERFACE
!
SUBROUTINE LOAD_SEGMENTS(HCARIN,K)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: K
END SUBROUTINE LOAD_SEGMENTS
!
END INTERFACE
!
END MODULE MODI_LOAD_SEGMENTS
!     ######spl
      SUBROUTINE LOAD_SEGMENTS(HCARIN,K)
!     ################################
!
!!****  *LOAD_SEGMENTS* - 
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
!!      Original       02/08/00
!!      Updated   PM   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY : JPHEXT
USE MODD_DIM1, ONLY : NIMAX,NJMAX
USE MODD_RESOLVCAR
USE MODD_GRID1
USE MODD_ALLOC_FORDIACHRO, ONLY : NGRIDIA
USE MODD_COORD, ONLY : XXX,XXY
USE MODE_GRIDPROJ
USE MODI_RESOLVXISOLEV

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER          :: K
CHARACTER(LEN=*) :: HCARIN
!
!*       0.1   Local variables
!              ---------------
INTEGER          :: IK,JM,J,IL, IMIN,II,IJ
INTEGER, DIMENSION(1):: IMINA
REAL              :: ZX,ZY,ZLAT,ZLON
REAL,DIMENSION(200) :: ZTEM
!------------------------------------------------------------------------------
IL=LEN_TRIM(HCARIN)
! IK= indice du 1er XSEGMx trouve puis mis a jour ensuite si plusieurs/ligne
IK=K
! Exploration de toute la ligne au cas ou plusieurs definitions/ligne
ZTEM(:)=9999.
CALL RESOLVXISOLEV(HCARIN(1:IL),IK,ZTEM)
DO J=SIZE(ZTEM,1),1,-1
  IF(ZTEM(J) /= 9999.)THEN
    JM=J
    EXIT
  ENDIF
ENDDO
!
IMIN=1
DO
  ! NSEGMS=0 ou 1, ici on cherche deux 0 consecutifs
  IMINA(1:1)=MINLOC(NSEGMS(IMIN:))
  IMIN=IMINA(1)+(IMIN-1)+1
  IF (NSEGMS(IMIN)==0 .OR. IMIN==SIZE(NSEGMS,1)) EXIT
ENDDO
XSEGMS(IMIN:IMIN-1+JM/2,1)=ZTEM(1:JM-1:2)
XSEGMS(IMIN:IMIN-1+JM/2,2)=ZTEM(2:JM:2)
DO J=IMIN,IMIN-1+JM/2
  if(nverbia >0)then
  print *,' J XSEGMS(J,:) ','J= ',J,' ',XSEGMS(J,:)
  endif
  ZLAT=XSEGMS(J,1)
  ZLON=XSEGMS(J,2)
  IF(ZLAT /= 0. .OR. ZLON /= 0.)THEN
    IF(HCARIN(K:K)=='X') THEN   ! XSEGMS
      NSEGMS(J)=1
    ENDIF
    IF(HCARIN(K:K)=='I') THEN   ! ISEGMS
      NSEGMS(J)=-1
    ENDIF
  ENDIF
! Conversion en coordonnees conformes
!maintenant dans oper_process et closf (juste avant le trace)
! IF(HCARIN(K:K)=='X') THEN   ! XSEGMS
!   CALL SM_XYHAT_S(XLATOR,XLONOR,ZLAT,ZLON,ZX,ZY)
!   XCONFSEGMS(J,1)=ZX
!   XCONFSEGMS(J,2)=ZY
! ENDIF
! IF(HCARIN(K:K)=='I') THEN   ! ISEGMS
!   II=MAX(MIN(INT(ZLAT),NIMAX+2*JPHEXT-1),1)
!   IJ=MAX(MIN(INT(ZLON),NJMAX+2*JPHEXT-1),1)
!   ZX=XXX(II,NGRIDIA(1)) +  &
!      (ZLAT-FLOAT(II))*(XXX(II+1,NGRIDIA(1)) - XXX(II,NGRIDIA(1)) )
!   ZY=XXY(IJ,NGRIDIA(1)) + &
!      (ZLON-FLOAT(IJ))*(XXY(IJ+1,NGRIDIA(1)) - XXY(IJ,NGRIDIA(1)) )
!   XCONFSEGMS(J,1)=ZX
!   XCONFSEGMS(J,2)=ZY
! ENDIF
ENDDO
do j=1,size(xsegms,1)
if(nverbia >0)then
print *,' J XSEGM ', J ,XSEGMS(J,:) 
endif
enddo
RETURN
END SUBROUTINE LOAD_SEGMENTS
