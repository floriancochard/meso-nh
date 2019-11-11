!     ######spl
      MODULE  MODI_CONV2XY
!     ####################
!
INTERFACE
!
SUBROUTINE CONV2XY(PXX,PYY,PX,PY,K)
REAL  ::  PXX,PYY,PX,PY
INTEGER,INTENT(IN)          :: K
END SUBROUTINE CONV2XY
!
END INTERFACE
!
END MODULE MODI_CONV2XY
!     ######spl
      SUBROUTINE CONV2XY(PXX,PYY,PX,PY,K)
!     ###################################
!
!!****  *CONV2XY* - 
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
!!      Original       16/06/98
!!      Updated   PM  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_COORD
USE MODD_DIM1
USE MODD_CONF
USE MODD_GRID1
USE MODD_GRID, ONLY: XLONORI,XLATORI
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO
USE MODD_FILES_DIACHRO     
USE MODE_GRIDPROJ

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------
!
REAL  ::  PXX,PYY,PX,PY
INTEGER,INTENT(IN)  ::  K
!
!*       0.1   Local variables
!              ---------------
!
INTEGER  ::  J,JM,JMCUR
INTEGER  ::  IINF, IJINF, ISUP, IJSUP
LOGICAL  ::  GOK
! !------------------------------------------------------------------------------
GOK=.FALSE.
IINF=NIINF; ISUP=NISUP; IJINF=NJINF; IJSUP=NJSUP
IF(ALLOCATED(XXHAT))THEN
ELSE
  IF (NBFILES == 1)THEN
  ELSE
    DO J=1,NBFILES
      IF(NUMFILES(J)==NUMFILECUR)THEN
	JMCUR=J
	if(nverbia > 0)then
	print *,' CONV2XY J JMCUR ',J,JMCUR
	endif
	EXIT
      ENDIF
    ENDDO
    DO J=1,NBFILES
      IF(NUMFILES(J)==NUMFILECUR)THEN
	CYCLE
      ELSE
	JM=J
	if(nverbia > 0 )THEN
       	  print *,' CONV2XY JM,CFILEDIAS(JM) ',JM,CFILEDIAS(JM)
	ENDIF
	CALL READ_FILEHEAD(JM,CFILEDIAS(JM),CLUOUTDIAS(JM))
	IF(NIMAX /= 0)THEN
	  GOK=.TRUE.
	  EXIT
	ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF
IF(ALLOCATED(XXHAT))THEN
IF(LCONV2XY .AND. NLATLON /= 0)THEN
  CALL SM_XYHAT_S(XLATORI,XLONORI,PXX,PYY, &
  PX,PY)
  IF(K == 11)THEN
    PXX=PX
  ELSE IF(K == 12)THEN
    PXX=PY
  ELSE IF(K == 21)THEN
    PYY=PX
  ELSE IF(K == 22)THEN
    PYY=PY
  ENDIF
ENDIF
ELSE
  print *,' Absence d''entete dans les differents fichiers ouverts'
  print *,' Impossibilite de convertir les coordonnees geographiques en conformes '
  print *,' LCONV2XY remis a .FALSE. '
  LCONV2XY=.FALSE.
ENDIF
IF(GOK)THEN
  CALL READ_FILEHEAD(JMCUR,CFILEDIAS(JMCUR),CLUOUTDIAS(JMCUR))
ENDIF
NIINF=IINF; NISUP=ISUP; NJINF=IJINF; NJSUP=IJSUP
RETURN
END SUBROUTINE CONV2XY
