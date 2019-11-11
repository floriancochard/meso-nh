!     ######spl
      MODULE MODI_RESOLV_UNITS
!     #############################
!
INTERFACE
!
SUBROUTINE RESOLV_UNITS(HCARIN,HCAROUT)
CHARACTER(LEN=*) :: HCARIN
CHARACTER(LEN=*) :: HCAROUT
END SUBROUTINE  RESOLV_UNITS
!
END INTERFACE
END MODULE MODI_RESOLV_UNITS
!     #######################################
      SUBROUTINE RESOLV_UNITS(HCARIN,HCAROUT)
!     #######################################
!
!!****  *RESOLV_UNITS* - Extraction du champ unites

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
!!      Original       06/06/94
!!      Updated   PM   02/12/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_RESOLVCAR
USE MODD_CONF

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

CHARACTER(LEN=*) :: HCARIN
CHARACTER(LEN=*)         :: HCAROUT
!
!*       0.1   Local variables
!              ---------------

!
CHARACTER(LEN=1)         :: YC
CHARACTER(LEN=LEN(HCARIN)) :: YCARIN
INTEGER   ::   ILENC
               
INTEGER   ::   J, J1, J2, JJ
!------------------------------------------------------------------------------
!
YCARIN=HCARIN
ILENC = LEN(YCARIN)
!print *,' YCARIN ',LEN(YCARIN),YCARIN
J1=0; J2=0
J1=INDEX(YCARIN,'(')
DO J=ILENC,1,-1
  IF(YCARIN(J:J) == ')')THEN
  J2=J
  EXIT
  ENDIF
ENDDO
CGROUP=ADJUSTL(CGROUP)
!print *,'CGROUP ',CGROUP
IF(J2 < J1)THEN
  J2=LEN_TRIM(YCARIN)+1
ENDIF
IF(J1 == 0 .AND. J2 == 0)THEN
  IF(INDEX(YCARIN,CGROUP(1:LEN_TRIM(CGROUP))) /= 0 )THEN
    HCAROUT(1:LEN(HCAROUT))=' '
  ELSE
    HCAROUT=ADJUSTL(YCARIN)
  ENDIF
ELSE
  HCAROUT=ADJUSTL(YCARIN(J1+1:J2-1))
ENDIF
!print *,' HCAROUT ',HCAROUT
YCARIN(1:LEN(YCARIN))=' '
!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE RESOLV_UNITS
