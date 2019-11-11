!     ######spl
      MODULE MODI_READREFINT_ISO 
!     ###########################
!
INTERFACE
!
SUBROUTINE READREFINT_ISO(HCARIN,PTABMN,PTABMX,PINT,PISOLEV)
CHARACTER(LEN=*)   :: HCARIN
REAL, INTENT(IN)   :: PTABMN,PTABMX
REAL               :: PINT
REAL, DIMENSION(:) :: PISOLEV
END SUBROUTINE READREFINT_ISO
!
END INTERFACE
END MODULE MODI_READREFINT_ISO
!     ######spl
      SUBROUTINE READREFINT_ISO(HCARIN,PTABMN,PTABMX,PINT,PISOLEV)
!     ###############################################
!
!!****  *READREFINT_ISO* - 
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
!!      Original       2/09/96
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

CHARACTER(LEN=*) :: HCARIN
REAL, INTENT(IN)   :: PTABMN,PTABMX
REAL               :: PINT
REAL, DIMENSION(:) :: PISOLEV
!
!*       0.1   Local variables
!              ---------------

INTEGER           :: IMASK,II,IIMIN,IIMAX,IIDEB,IIFIN,INBISO
INTEGER           :: J,JM
REAL              :: ZMEMINT,ZREF,ZVALMIN,ZVALMAX
LOGICAL           :: GOKREF, GOKINT
REAL, DIMENSION(SIZE(PISOLEV)) :: ZISOLEV
CHARACTER(LEN=LEN(HCARIN)) :: YCARIN, YCARIN2
!
!------------------------------------------------------------------------------
GOKREF=.FALSE.
GOKINT=.FALSE.
!
YCARIN(1:LEN(YCARIN))=' '
HCARIN=ADJUSTL(HCARIN)
YCARIN=HCARIN
IMASK=INDEX(YCARIN,'MASK')
IF(IMASK /=0)THEN
DO J=1,LEN(YCARIN)
 IF(YCARIN(J:J) == ' ')THEN
   JM=J-1
   EXIT
 ENDIF
ENDDO
YCARIN(1:LEN(YCARIN))=' '
YCARIN=HCARIN(JM+2:LEN_TRIM(HCARIN))
YCARIN=ADJUSTL(YCARIN)
ENDIF
JM=0
DO J=1,LEN(YCARIN)
 IF(YCARIN(J:J) == ' ')THEN
   JM=J-1
   EXIT
 ENDIF
ENDDO
IF(JM /= 0)THEN
  YCARIN2(1:LEN(YCARIN2))=' '
  YCARIN2=YCARIN(1:JM)
  YCARIN(1:LEN(YCARIN))=' '
  YCARIN=ADJUSTL(YCARIN2)
ENDIF
!
ZMEMINT=PINT
!
IF(NBISOREF == 0)THEN
  GOKREF=.FALSE.
  print *,' AUCUN REF USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
ELSE
  DO J=1,NBISOREF
    IF(YCARIN(1:LEN_TRIM(YCARIN)) == CISOREF(J)(1:LEN_TRIM(YCARIN)))THEN
      ZREF=XISOREFP(J)
      GOKREF=.TRUE.
      EXIT
    ENDIF
  ENDDO
  IF(.NOT.GOKREF)THEN
    print *,' AUCUN REF USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
  ENDIF
ENDIF
!
IF(NBISOINT == 0)THEN
  GOKINT=.FALSE.
  print *,' AUCUN INT USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
ELSE
  DO J=1,NBISOINT
    IF(YCARIN(1:LEN_TRIM(YCARIN)) == CISOINT(J)(1:LEN_TRIM(YCARIN)))THEN
      PINT=XISOINT(J)
      GOKINT=.TRUE.
      EXIT
    ENDIF
  ENDDO
  IF(.NOT.GOKINT)THEN
    print *,' AUCUN INT USER ENREGISTRE POUR :  ',YCARIN(1:LEN_TRIM(YCARIN))
  ENDIF
ENDIF
IF(.NOT.GOKREF .OR. .NOT.GOKINT)THEN
  LISOREF=.FALSE.
  print *,' UTILISATION DES VALEURS DE XISOREF,XDIAINT POUR : ',YCARIN(1:LEN_TRIM(YCARIN))
ELSE
  LISOREF=.TRUE.
ENDIF
!------------------------------------------------------------------------------

IF(.NOT. LISOREF)THEN
  PINT=XDIAINT
  IF(PINT == 0.)THEN
    PINT=ZMEMINT
  ENDIF
  ZREF=XISOREF
  IF (ZREF.LT.PTABMN .OR. ZREF.GT.PTABMX) THEN
if (nverbia>5) then
  print*,'TABmin-max= ',PTABMN,PTABMX
  print*,'ISO REF hors des valeurs extremes du champ = ',XISOREF
endif
    ZREF=0.5*(PTABMN+PTABMX)
if (nverbia>5) then
  print*,'ISO REF calcule = ',ZREF
endif
  ENDIF
ELSE
  LISOREF=.FALSE.
ENDIF
!------------------------------------------------------------------------------
ZISOLEV(:)=0.
ZVALMIN=ZREF ; ZVALMAX=ZREF
! ZISOLEV contient les valeurs des differentes isolignes a tracer
!rempli ainsi: ZREF -PINT +PINT -2.PINT +2.PINT ...
II=1 ; IIMIN=II ; IIMAX=II
ZISOLEV(1)=ZREF
DO J=1,SIZE(ZISOLEV)
  ZVALMIN=ZVALMIN-PINT
  IF (ZVALMIN.GT.PTABMN) THEN
    II=II+1
    ZISOLEV(II)=ZVALMIN
    IIMIN=II
  ENDIF
  ZVALMAX=ZVALMAX+PINT
  IF (ZVALMAX.LT.PTABMX) THEN
    II=II+1
    ZISOLEV(II)=ZVALMAX
    IIMAX=II
  ENDIF
ENDDO
if (nverbia>=5) then
  print*,'IIMIN,IIMAX,II= ',IIMIN,IIMAX,II
endif
if (nverbia>5) then
  print*,'ZISOLEV= ',ZISOLEV
endif
! 
! reordonne pour PISOLEV de la valeur min a la valeur max
INBISO=II
IF (INBISO.LE.2) THEN
  PISOLEV(1)=ZISOLEV(1)
  PISOLEV(2)=ZISOLEV(2)
ELSE
  II=1
  IF (IIMIN .GT. (IIMAX+1)) THEN   ! premiers min contigus
    DO J=IIMIN,IIMAX+1,-1
      PISOLEV(II)=ZISOLEV(J)
      II=II+1
    END DO
    IIDEB=IIMAX+1-2
  ELSE
    IIDEB=IIMIN
  ENDIF
  !
  IF (IIDEB.GT.0) THEN             ! traite les valeurs inf a ZREF
    ! une valeur sur 2 pour les min suivants
    DO J=IIDEB,2,-2
      PISOLEV(II)=ZISOLEV(J)
      II=II+1
    END DO
    IIFIN=MIN(IIMAX,IIMIN+1)
    ! une valeur sur 2 pour les premiers max
    DO J=1,IIFIN,2
      PISOLEV(II)=ZISOLEV(J)
      II=II+1
    END DO
  ELSE                             ! toutes les valeurs sont sup a ZREF
    IIFIN=0
  ENDIF
  !
  IF (IIMAX.GT.IIMIN+1) THEN       ! derniers max contigus
    DO J=IIFIN+1,IIMAX
      PISOLEV(II)=ZISOLEV(J)
      II=II+1
    ENDDO
  ENDIF
ENDIF
if (nverbia>5) then
  print*,'II= ',II
endif

END SUBROUTINE READREFINT_ISO
