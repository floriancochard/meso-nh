!     ######spl
      SUBROUTINE GENFORMAT_FORDIACHRO(PCLV,HLLBS)
!     ###########################################
!
!!****  *GENFORMAT* - Determination du format des valeurs d'isocontours en 
!                     legende
!!
!!    PURPOSE
!!    -------
!       Pour une valeur d'isocontour donnee, recherche le format le mieux
!       adapte pour cette valeur et l'ecrit dans une chaine de caracteres
!       suivant le dit format
!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      None
!!
!!    AUTHOR
!!    ------
!!	
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       25/01/95
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
COMMON/GENF/NBCU
!
!*       0.1  Dummy arguments
!          

REAL                :: PCLV
REAL                :: ZEPS
CHARACTER(LEN=*)    :: HLLBS
!
!*       0.2  local variables
!          
REAL :: ZALOG10
INTEGER :: I7,I8, NBCU
!
!-------------------------------------------------------------------------------
!print *,' ENTREE genformat PCLV HLLBS',PCLV,HLLBS
I7=0; I8=0
ZEPS=1.E-30
HLLBS(1:LEN(HLLBS))=' '
IF(PCLV == 0. .OR. (ABS(PCLV) >=0.01 .AND. ABS(PCLV) <= 1.))THEN
  WRITE(HLLBS,'(F6.3)')PCLV
  I7=6
ELSE
  ZALOG10=ALOG10(ABS(PCLV))
  IF(ZALOG10 < 0.)THEN
    IF(PCLV >= 0.)THEN
      WRITE(HLLBS,'(E7.2)')PCLV
      I7=7
    ELSE
      IF(ABS(ZALOG10) <= 10)THEN
        WRITE(HLLBS,'(E7.2E1)')PCLV
        I7=7
      ELSE
        WRITE(HLLBS,'(E8.2)')PCLV
        I8=8
      ENDIF
    ENDIF
  ELSE
    IF(ZALOG10 >= 5.)THEN
      IF(PCLV >= 0.)THEN
        WRITE(HLLBS,'(E7.2)')PCLV
        I7=7
      ELSE
        WRITE(HLLBS,'(E8.2)')PCLV
        I8=8
      ENDIF
    ENDIF
    IF(ZALOG10 >= 4. .AND. ZALOG10 < 5.)THEN
      WRITE(HLLBS,'(F7.0)')PCLV
      I7=7
    ENDIF
    IF(ZALOG10 < 4)THEN
    IF(ZALOG10 >= 3.-ZEPS .AND. ZALOG10 < 4.)WRITE(HLLBS,'(F6.0)')PCLV
    IF(ZALOG10 >= 2.-ZEPS .AND. ZALOG10 < 3.)WRITE(HLLBS,'(F6.1)')PCLV
    IF(ZALOG10 >= 1.-ZEPS .AND. ZALOG10 < 2.)WRITE(HLLBS,'(F6.2)')PCLV
    IF(ZALOG10 >= 0. .AND. ZALOG10 < 1.-ZEPS)WRITE(HLLBS,'(F6.3)')PCLV
      I7=6
    ENDIF
  END IF
END IF
HLLBS=ADJUSTL(HLLBS)
!print *,' SORTIE genformat PCLV HLLBS',PCLV,HLLBS
NBCU=MAX(I7,I8)
!print *,' NBCU ',NBCU
RETURN
END SUBROUTINE GENFORMAT_FORDIACHRO
