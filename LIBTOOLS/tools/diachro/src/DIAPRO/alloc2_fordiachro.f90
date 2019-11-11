!     ######spl
      MODULE  MODI_ALLOC2_FORDIACHRO
!     ##############################
!
INTERFACE
!
SUBROUTINE ALLOC2_FORDIACHRO(KOP)
INTEGER :: KOP
END SUBROUTINE ALLOC2_FORDIACHRO
!
END INTERFACE
!
END MODULE MODI_ALLOC2_FORDIACHRO
!     ######spl
      SUBROUTINE ALLOC2_FORDIACHRO(KOP)
!     #################################
!
!!****  *ALLOC2_FORDIACHRO* - 
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
!!      Original       08/01/96
!!      Updated   PM 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ALLOC_FORDIACHRO
USE MODD_ALLOC2_FORDIACHRO
USE MODD_RESOLVCAR
USE MODD_PT_FOR_CH_FORDIACHRO
USE MODD_FILES_DIACHRO

IMPLICIT NONE
!
!*       0.1   Dummy arguments
!              ---------------

INTEGER  :: KOP
!
!*       0.1   Local variables
!              ---------------

!
!------------------------------------------------------------------------------
!
IF (KOP == 1)THEN

  ALLOCATE(XDATIME2(SIZE(XDATIME,1),SIZE(XDATIME,2)))
  XDATIME2(:,:)=XDATIME(:,:)
  ALLOCATE(XVAR2(SIZE(XVAR,1),SIZE(XVAR,2),SIZE(XVAR,3),  &
  SIZE(XVAR,4),SIZE(XVAR,5),SIZE(XVAR,6)))
  XVAR2(:,:,:,:,:,:)=XVAR(:,:,:,:,:,:)
! print *,' XVAR2 ',XVAR2
  IF(ALLOCATED(XU))THEN
  ALLOCATE(XUMEM(SIZE(XU,1),SIZE(XU,2),SIZE(XU,3),  &
  SIZE(XU,4),SIZE(XU,5),SIZE(XU,6)))
  XUMEM(:,:,:,:,:,:)=XU(:,:,:,:,:,:)
  if(nverbia > 0)THEN
    print *,' ** ALLOC2 XUMEM alloue'
  endif
  ENDIF
  IF(ALLOCATED(XV))THEN
  ALLOCATE(XVMEM(SIZE(XV,1),SIZE(XV,2),SIZE(XV,3),  &
  SIZE(XV,4),SIZE(XV,5),SIZE(XV,6)))
  XVMEM(:,:,:,:,:,:)=XV(:,:,:,:,:,:)
  ENDIF
  ALLOCATE(XTRAJT2(SIZE(XTRAJT,1),SIZE(XTRAJT,2)))
  XTRAJT2(:,:)=XTRAJT(:,:)
  ALLOCATE(NGRIDIA2(SIZE(NGRIDIA)))
  NGRIDIA2(:)=NGRIDIA(:)
  ALLOCATE(CTITRE2(SIZE(CTITRE)))
  CTITRE2(:)(1:LEN(CTITRE2))=' '
  CTITRE2(:)=CTITRE(:)
  ALLOCATE(CUNITE2(SIZE(CUNITE)))
  CUNITE2(:)(1:LEN(CUNITE2))=' '
  CUNITE2(:)=CUNITE(:)
  ALLOCATE(CCOMMENT2(SIZE(CCOMMENT)))
  CCOMMENT2(:)(1:LEN(CCOMMENT2))=' '
  CCOMMENT2(:)=CCOMMENT(:)

  IF(ALLOCATED(XTRAJX))THEN
    ALLOCATE(XTRAJX2(SIZE(XTRAJX,1),SIZE(XTRAJX,2),SIZE(XTRAJX,3)))
    XTRAJX2(:,:,:)=XTRAJX(:,:,:)
  ENDIF
  IF(ALLOCATED(XTRAJY))THEN
    ALLOCATE(XTRAJY2(SIZE(XTRAJY,1),SIZE(XTRAJY,2),SIZE(XTRAJY,3)))
    XTRAJY2(:,:,:)=XTRAJY(:,:,:)
  ENDIF
  IF(ALLOCATED(XTRAJZ))THEN
    ALLOCATE(XTRAJZ2(SIZE(XTRAJZ,1),SIZE(XTRAJZ,2),SIZE(XTRAJZ,3)))
    XTRAJZ2(:,:,:)=XTRAJZ(:,:,:)
  ENDIF

  IF (ALLOCATED(XMASK))THEN
    ALLOCATE(XMASK2(SIZE(XMASK,1),SIZE(XMASK,2),SIZE(XMASK,3), &
    SIZE(XMASK,4),SIZE(XMASK,5),SIZE(XMASK,6)))
    XMASK2(:,:,:,:,:,:)=XMASK(:,:,:,:,:,:)
  ENDIF
  NUMFILECUR2=NUMFILECUR

ELSE

  IF (ALLOCATED(XMASK2))THEN
    DEALLOCATE(XMASK2)
  ENDIF
  IF (ALLOCATED(XTRAJZ2))THEN
    DEALLOCATE(XTRAJZ2)
  ENDIF
  IF (ALLOCATED(XTRAJY2))THEN
    DEALLOCATE(XTRAJY2)
  ENDIF
  IF (ALLOCATED(XTRAJX2))THEN
    DEALLOCATE(XTRAJX2)
  ENDIF
  DEALLOCATE(CCOMMENT2,CUNITE2,CTITRE2)
  IF(ALLOCATED(NGRIDIA2))THEN
    DEALLOCATE(NGRIDIA2)
  ENDIF
  DEALLOCATE(XTRAJT2)
! DEALLOCATE(XVAR2,XTRAJT2,CTITRE2,CUNITE2,CCOMMENT2)
  IF(ALLOCATED(XVMEM))THEN
    DEALLOCATE(XVMEM)
  ENDIF
  IF(ALLOCATED(XUMEM))THEN
    DEALLOCATE(XUMEM)
  if(nverbia > 0)THEN
    print *,' ** ALLOC2 XUMEM desalloue'
  endif
  ENDIF
  DEALLOCATE(XVAR2)
  DEALLOCATE(XDATIME2)

ENDIF

!
!-----------------------------------------------------------------------------
!
!*       2.       EXITS
!                 -----
! 
RETURN
END SUBROUTINE ALLOC2_FORDIACHRO
