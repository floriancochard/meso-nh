!     ######spl
      MODULE  MODI_LOAD_FMTAXES
!     #########################
!
INTERFACE
!
SUBROUTINE LOAD_FMTAXES(HCARIN,K)
CHARACTER(LEN=*) :: HCARIN
INTEGER          :: K
END SUBROUTINE LOAD_FMTAXES
!
END INTERFACE
!
END MODULE MODI_LOAD_FMTAXES
!     ######spl
      SUBROUTINE LOAD_FMTAXES(HCARIN,K)
!     ################################
!
!!****  *LOAD_FMTAXES* - 
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
USE MODD_RESOLVCAR

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER          :: K
CHARACTER(LEN=*) :: HCARIN
!
!*       0.1   Local variables
!              ---------------
INTEGER          :: IEGAL,IQ1,IQ2
! !------------------------------------------------------------------------------
!nverbia=6
IEGAL=INDEX(HCARIN,'=')
IQ2=LEN_TRIM(HCARIN)
IQ1=INDEX(HCARIN,'"')
IF(IQ1 == 0)THEN
  IQ1=INDEX(HCARIN,"'")
ENDIF
IF(IQ1 == 0 .OR. IQ1 == IQ2)THEN
  IQ1=IEGAL
ENDIF
IF(HCARIN(IQ2:IQ2) == "'" .OR. HCARIN(IQ2:IQ2) == '"')THEN
ELSE
  IQ2=IQ2+1
ENDIF
!print *,' HCARIN(K:IEGAL-1) ',HCARIN(K:IEGAL-1)
IF(HCARIN(K:IEGAL-1) == 'CFMTAXEX')THEN
  CFMTAXEX='          '
  CFMTAXEX=HCARIN(IQ1+1:IQ2-1)
  CFMTAXEX=ADJUSTL(CFMTAXEX)
! CFMTAXEX="'"//HCARIN(IQ1+1:IQ2-1)//"'"
  if(nverbia >0)then
    print *,' CFMTAXEX=',CFMTAXEX
  endif
ELSEIF(HCARIN(K:IEGAL-1) == 'CFMTAXEY')THEN
  CFMTAXEY='          '
  CFMTAXEY=HCARIN(IQ1+1:IQ2-1)
  CFMTAXEY=ADJUSTL(CFMTAXEY)
! CFMTAXEY="'"//HCARIN(IQ1+1:IQ2-1)//"'"
  if(nverbia >0)then
    print *,' CFMTAXEY=',CFMTAXEY
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 19/12/2008 : modification pour controler la taille et le format des labels !!
!! pour les retrotrajectoires                                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ELSEIF(HCARIN(K:IEGAL-1) == 'CFMTRTRAJ')THEN
  CFMTRTRAJ='          '
  CFMTRTRAJ=HCARIN(IQ1+1:IQ2-1)
  CFMTRTRAJ=ADJUSTL(CFMTRTRAJ)
  if(nverbia >0)then
    print *,' CFMTRTRAJ=',CFMTRTRAJ
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE
  print *, ' Erreur Passage ds LOAD_FMTAXES mais la variable n''est ni CFMTAXEX ni CFMTAXEY ni CFMTRTRAJ'
ENDIF
RETURN
END SUBROUTINE LOAD_FMTAXES
