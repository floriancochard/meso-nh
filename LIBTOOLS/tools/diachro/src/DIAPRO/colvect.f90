!     ##############################
      SUBROUTINE COLVECT(KKU,PTEM2D)
!     ##############################
!
!!****  *COLVECT* -  Couleur fleches par un autre parametre
!! Possible uniquement pour les profils verticaux de vecteurs vent horizontal
!! generes directement ds un fichier diachronique (CART + MASK)
!!****           
!!
!!    PURPOSE
!!    -------
!       
!      
!     
!
!!**  METHOD
!!    ------
!!    
!!   
!!
!!    EXTERNAL
!!    --------
!!      
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODN_PARA  : Defines NAM_DOMAIN_POS namelist (former PARA common)
!!         NLANGLE :  Angle between X Meso-NH axis and
!!                    cross-section direction in degrees
!!                    (Integer value anticlockwise)
!!
!!    REFERENCE
!!    ---------
!!
!!      MESO-NH User's Manual, TRACE Post Processing sections, Version 1.0:
!!       + Book1: Concepts and Fundamentals, to appear in 1994;
!!       + Book2: Technical Reference and Flowcharts, to appear in 1994;
!!       + Book3: Tutorial, November 1994.
!!
!!
!!    AUTHOR
!!    ------
!!      J. Duron    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       23/10/2001
!!      Updated   
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODN_NCAR
USE MODD_PVT
USE MODD_RESOLVCAR
!
IMPLICIT NONE
!
!*       0.1  Dummy arguments and results
!
                                              ! 
                                              !
REAL, DIMENSION(:,:),  INTENT(IN) :: PTEM2D !
                                              !
INTEGER   :: KKU                              !
                                              ! 
!
!*       0.2  Local variables
!
INTEGER             :: JILOOP, JJLOOP, JKLOOP
!
REAL                :: ZMXPARCOL, ZMNPARCOL, ZINTPARCOL
!
!-------------------------------------------------------------------------------
!
!*        1.   COMPUTING THE LONGITUDINAL AND TRANSVERSE COMPONENTS
!              ----------------------------------------------------
!
!*        1.1  
!
IF(ALLOCATED(NCOL2DUV))THEN
  DEALLOCATE(NCOL2DUV)
ENDIF
ALLOCATE(NCOL2DUV(SIZE(PTEM2D,2),KKU))
LCOLPVT=.TRUE.
NCOL2DUV=1
IF(LCOLUSERUV)THEN !:::::::::::::::::::::::::::
  DO JILOOP=1,SIZE(PTEM2D,1)
  DO JJLOOP=1,SIZE(PTEM2D,2)

  IF(PTEM2D(JILOOP,JJLOOP) /= XSPVAL)THEN
  IF(PTEM2D(JILOOP,JJLOOP) < XPARCOLUV(1))THEN
    NCOL2DUV(JJLOOP,JILOOP)=NINDCOLUV(1)
  ELSE IF(PTEM2D(JILOOP,JJLOOP) >= XPARCOLUV&
  (NBPARCOLUV))THEN
    NCOL2DUV(JJLOOP,JILOOP)=NINDCOLUV(NBCOLUV)
  ELSE
    DO JKLOOP=2,NBPARCOLUV
       IF(PTEM2D(JILOOP,JJLOOP) >= XPARCOLUV(&
       JKLOOP-1) .AND. PTEM2D(JILOOP,JJLOOP)<&
        XPARCOLUV(JKLOOP))then
         NCOL2DUV(JJLOOP,JILOOP)=NINDCOLUV(&
         JKLOOP)
         EXIT
       ENDIF
    ENDDO
  ENDIF
ENDIF

ENDDO
ENDDO
      ELSE               !:::::::::::::::::::::::::::
ZMXPARCOL=-1.e14
ZMNPARCOL=+1.e14
DO JILOOP=1,SIZE(PTEM2D,1)
DO JJLOOP=1,SIZE(PTEM2D,2)
                            
  IF(PTEM2D(JILOOP,JJLOOP) /= XSPVAL)THEN
    ZMXPARCOL=MAX(PTEM2D(JILOOP,JJLOOP),ZMXPARCOL)
    ZMNPARCOL=MIN(PTEM2D(JILOOP,JJLOOP),ZMNPARCOL)
  ENDIF
ENDDO
ENDDO
IF(ABS(ZMXPARCOL-ZMNPARCOL) >= 20)THEN
 ZMNPARCOL=ZMNPARCOL+1
 ZMXPARCOl=ZMXPARCOl-1
ENDIF
ZINTPARCOL=(ZMXPARCOL-ZMNPARCOL)/5.
XPARCOLUVSTD(1)=ZMNPARCOL
DO JILOOP=2,NBPARCOLUVSTD-1
  XPARCOLUVSTD(JILOOP)=XPARCOLUVSTD(JILOOP-1)+&
  ZINTPARCOL
ENDDO
XPARCOLUVSTD(NBPARCOLUVSTD)=ZMXPARCOL
if(nverbia > 0)then
print *,' **OPER_UV** XPARCOLUVSTD ',XPARCOLUVSTD
endif
DO JILOOP=1,SIZE(PTEM2D,1)
DO JJLOOP=1,SIZE(PTEM2D,2)

IF(PTEM2D(JILOOP,JJLOOP) /= XSPVAL)THEN
  IF(PTEM2D(JILOOP,JJLOOP) < XPARCOLUVSTD(1))THEN
    NCOL2DUV(JJLOOP,JILOOP)=NCOLUVSTD(1)
  ELSE IF(PTEM2D(JILOOP,JJLOOP) >= XPARCOLUVSTD&
  (NBPARCOLUVSTD))THEN
    NCOL2DUV(JJLOOP,JILOOP)=NCOLUVSTD(NBCOLUVSTD)
  ELSE
    DO JKLOOP=2,NBPARCOLUVSTD
       IF(PTEM2D(JILOOP,JJLOOP) >= XPARCOLUVSTD(&
       JKLOOP-1) .AND. PTEM2D(JILOOP,JJLOOP)<&
        XPARCOLUVSTD(JKLOOP))then
         NCOL2DUV(JJLOOP,JILOOP)=NCOLUVSTD(&
         JKLOOP)
         EXIT
       ENDIF
    ENDDO
  ENDIF
ENDIF


ENDDO
ENDDO

      ENDIF             !::::::::::::::::::::::::::::::::::
!
!*        1.2  
!*            
!*           
!*          
!
!
!*       1.3   
!
IF(nverbia > 0)THEN
 print *,' ** colvect '
endif
!
!------------------------------------------------------------------------------
!
!*        2.     EXIT
!                ----
!
RETURN
END SUBROUTINE COLVECT
