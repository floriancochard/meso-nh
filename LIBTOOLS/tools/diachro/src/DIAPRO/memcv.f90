!     #################
      SUBROUTINE MEMCV
!     #################
!
!!****  *MEMCV* - 
!!                                                            
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       15/03/99
!!      Updated   PM  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_NMGRID
USE MODD_COORD
USE MODD_MEMCV
USE MODN_PARA
USE MODD_RESOLVCAR, ONLY : XCONFSEGMS,NSEGMS,NVERBIA
!
IMPLICIT NONE
!
!*       0.1   Local variables
!
INTEGER           :: JILOOP,IMIN
INTEGER,DIMENSION(1):: IMINA
!
REAL,DIMENSION(:,:),ALLOCATABLE :: ZX, ZY
!
!-------------------------------------------------------------------------------
!
!*      1. 
!              ----------------------------
!
!
IF(NTRACECV > 0)THEN
  DO JILOOP=1,NTRACECV
    IF(XTRACECV(1,JILOOP)==XDSX(1,NMGRID) .AND. XTRACECV(2,JILOOP)==XDSX(NLMAX,NMGRID) .AND.&
       XYTRACECV(1,JILOOP)==XDSY(1,NMGRID) .AND. XYTRACECV(2,JILOOP)==XDSY(NLMAX,NMGRID))THEN
      RETURN
    ENDIF
  ENDDO
  ALLOCATE(ZX(2,SIZE(XTRACECV,2)))
  ALLOCATE(ZY(2,SIZE(XYTRACECV,2)))
  ZX=XTRACECV
  ZY=XYTRACECV
  NTRACECV=NTRACECV+1
  DEALLOCATE(XTRACECV)
  DEALLOCATE(XYTRACECV)
  ALLOCATE(XTRACECV(2,NTRACECV))
  ALLOCATE(XYTRACECV(2,NTRACECV))
  XTRACECV(:,1:NTRACECV-1)=ZX
  XYTRACECV(:,1:NTRACECV-1)=ZY
  XTRACECV(1,NTRACECV)=XDSX(1,NMGRID)
  XTRACECV(2,NTRACECV)=XDSX(NLMAX,NMGRID)
  XYTRACECV(1,NTRACECV)=XDSY(1,NMGRID)
  XYTRACECV(2,NTRACECV)=XDSY(NLMAX,NMGRID)
  DEALLOCATE(ZX)
  DEALLOCATE(ZY)
ELSE
  NTRACECV=NTRACECV+1
  IF(ALLOCATED(XTRACECV))THEN
    DEALLOCATE(XTRACECV)
  ENDIF
  IF(ALLOCATED(XYTRACECV))THEN
    DEALLOCATE(XYTRACECV)
  ENDIF
  ALLOCATE(XTRACECV(2,NTRACECV))
  ALLOCATE(XYTRACECV(2,NTRACECV))
  XTRACECV(1,NTRACECV)=XDSX(1,NMGRID)
  XTRACECV(2,NTRACECV)=XDSX(NLMAX,NMGRID)
  XYTRACECV(1,NTRACECV)=XDSY(1,NMGRID)
  XYTRACECV(2,NTRACECV)=XDSY(NLMAX,NMGRID)
ENDIF
! stockage dans segments pour trace de la CV dans CH suivante(s)
IF(LTRACECV) THEN
  IMIN=1
  DO
    ! NSEGMS=0 ou 1, ici on cherche deux 0 consecutifs
    IMINA(1:1)=MINLOC(NSEGMS(IMIN:)) 
    IMIN=IMINA(1)+(IMIN-1)+1
    IF (NSEGMS(IMIN)==0 .OR. IMIN==SIZE(NSEGMS,1)) EXIT
  ENDDO
  NSEGMS(IMIN)=2    
  XCONFSEGMS(IMIN,1)=XDSX(1,NMGRID)
  XCONFSEGMS(IMIN,2)=XDSY(1,NMGRID)
  NSEGMS(IMIN+1)=2    
  XCONFSEGMS(IMIN+1,1)=XDSX(NLMAX,NMGRID)
  XCONFSEGMS(IMIN+1,2)=XDSY(NLMAX,NMGRID)
END IF
!
!------------------------------------------------------------------------------
!
!*      2.    EXIT
!             ----
!
RETURN
END SUBROUTINE MEMCV
