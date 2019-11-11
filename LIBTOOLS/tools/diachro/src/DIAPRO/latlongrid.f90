!     ######spl
      SUBROUTINE LATLONGRID
!     ################################
!
!!****  *LATLONGRID* - 
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
USE MODD_NMGRID
USE MODD_RESOLVCAR
USE MODD_ALLOC_FORDIACHRO

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

!
!*       0.1   Local variables
!              ---------------
! !------------------------------------------------------------------------------
NLATLON=0
NLATLON=INDEX(CCOMMENT(NLOOPP),'LATLON')
IF(NLATLON == 0)THEN
  NLATLON=INDEX(CCOMMENT(NLOOPP),'Latlon')
ENDIF
IF(NLATLON == 0)THEN
  NLATLON=INDEX(CCOMMENT(NLOOPP),'latlon')
ENDIF
IF(NLATLON == 0)THEN
  NLATLON=INDEX(CCOMMENT(NLOOPP),'LatLon')
ENDIF
IF(NVERBIA > 5)THEN
  print *,' NLATLON,CCOMMENT(NLOOPP) ',NLATLON,CCOMMENT(NLOOPP)
ENDIF
NMGRID=NGRIDIA(NLOOPP)
IF(NMGRID <1 .OR. NMGRID >7)THEN
  PRINT *,' VALEUR NMGRID ABERRANTE: ',NMGRID, &
        '        FORCEE A        :  1'
  NMGRID=1
ENDIF
RETURN
END SUBROUTINE LATLONGRID
