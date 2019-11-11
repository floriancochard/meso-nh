!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:./s.interp3d.f90, Version:1.1, Date:03/06/05, Last modified:01/10/10
!-----------------------------------------------------------------
!     ######spl
MODULE MODI_PINTER
!#################################
!
INTERFACE
      SUBROUTINE PINTER(PFIELD,KGRID,PSVAL,PPLEV,PFIELDAP,PPABSHO)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD      ! values of the field
INTEGER,                INTENT(IN) :: KGRID       ! Mesonh grid indicator
REAL,                   INTENT(IN) :: PSVAL       ! value for missing data
REAL, DIMENSION(:),     INTENT(IN) :: PPLEV       ! list of vertical levels(hPa)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PFIELDAP    ! values of the field on the pressure levels
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PPABSHO     ! abs. pressure when hor. interpolation
END SUBROUTINE  PINTER
END INTERFACE
END MODULE MODI_PINTER
!     ######spl
      SUBROUTINE PINTER(PFIELD,KGRID,PSVAL,PPLEV,PFIELDAP,PPABSHO)
!     #####################
!
!!****  *PINTER* - interpole 3D fields on pressure levels
!!                         
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      Functions MXF, MYF, MZF
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_FIELD1    : contains prognostics  variables 
!!         XPASBM
!!      Module MODD_GRID1
!!         XZZ
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	
!!      V.Ducrocq  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    21/03/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
! 
USE MODD_PARAMETERS
USE MODD_DIM1
!
USE MODI_SHUMAN ! interface modules 
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
REAL, DIMENSION(:,:,:), INTENT(IN) :: PFIELD      ! values of the field
INTEGER,                INTENT(IN) :: KGRID       ! Mesonh grid indicator
REAL,                   INTENT(IN) :: PSVAL       ! value for missing data
REAL, DIMENSION(:),     INTENT(IN) :: PPLEV       ! list of vertical levels(hPa)
REAL, DIMENSION(:,:,:), INTENT(OUT):: PFIELDAP    ! values of the field on the pressure levels
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PPABSHO    ! abs. pressure (when hor. interpolation IGRID=0)
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: JKP,JKLOOP,JJLOOP,JILOOP,IJ,II ! loop indices
INTEGER      :: IIE,IJE,IPU      ! End of usefull area 
INTEGER      :: IIB,IJB,IKB      ! Begining of usefull area 
REAL, DIMENSION(:,:,:),ALLOCATABLE :: ZPTH  ! pressure for grid points corresponding to KGRID type 
REAL  :: ZREF,ZXP,ZXM,ZDIXEPS ! pressure values and epsilon value
!-------------------------------------------------------------------------------
!
!*         1.    
!                 ------------
IPU=SIZE(PFIELDAP,3)
IKB=1 +JPVEXT
ZDIXEPS=10.*EPSILON(1.)
!
ALLOCATE(ZPTH(SIZE(PPABSHO,1),SIZE(PPABSHO,2),SIZE(PPABSHO,3)))
IIB=JPHEXT+1
IIE=NIMAX+JPHEXT
IJB=JPHEXT+1
IJE=NJMAX+JPHEXT
SELECT CASE (KGRID)
CASE(0)
  ZPTH(:,:,:)=PPABSHO(:,:,:)
  IIB=1
  IIE=SIZE(PPABSHO,1)
  IJB=1
  IJE=SIZE(PPABSHO,2)
CASE(1)
  ZPTH=PPABSHO
CASE(2)
  ZPTH(:,:,:)=MXM(PPABSHO(:,:,:))
  ZPTH(1,:,:)=2.*ZPTH(2,:,:) - ZPTH(3,:,:)
CASE(3)
  ZPTH(:,:,:)=MYM(PPABSHO(:,:,:))
    ZPTH(:,1,:)=2.*ZPTH(:,2,:) - ZPTH(:,3,:)
  CASE(4)
    ZPTH(:,:,:)=MZM(PPABSHO(:,:,:))
    ZPTH(:,:,1)=2.*ZPTH(:,:,2) - ZPTH(:,:,3)
  END SELECT
!
!
DO JKP= 1, IPU
   ZREF=ALOG10(PPLEV(JKP)*100.)
   DO JILOOP = IIB,IIE
      DO JJLOOP = IJB,IJE
         IJ=JJLOOP-IJB+1
         II=JILOOP-IIB+1
         PFIELDAP(II,IJ,JKP)=PSVAL
         DO JKLOOP = 1,NKMAX+2*JPVEXT
            IF (ZPTH(JILOOP,JJLOOP,JKLOOP)==PSVAL) CYCLE
            ZXM=ALOG10(ZPTH(JILOOP,JJLOOP,JKLOOP))
            ZXP=ALOG10(ZPTH(JILOOP,JJLOOP,MIN(NKMAX+2*JPVEXT,JKLOOP+1)))
            IF ((ZXP-ZREF)*(ZREF-ZXM) .GE.0.) THEN
               IF (JKLOOP+1 == IKB) THEN
                  CYCLE
               ELSE
                  GO TO 4
               ENDIF
            ELSE IF (ZXP.GE.ZXM-ZDIXEPS.AND.ZXP.LE.ZXM+ZDIXEPS.AND.  &
      ZREF.GE.ZXM-ZDIXEPS.AND.ZREF.LE.ZXM+ZDIXEPS) THEN
               IF(JKLOOP+1 == IKB)THEN
                  CYCLE
               ELSE
                  GO TO 4
               ENDIF
            END IF
         END DO
         GO TO 3
4     CONTINUE
!
!  We interpolate 
         PFIELDAP(II,IJ,JKP)= (PFIELD(II,IJ,JKLOOP)* (ZXP-ZREF)+ &
              PFIELD(II,IJ,MIN(NKMAX+2*JPVEXT,JKLOOP+1))* (ZREF-ZXM)) &
              / MIN(-1.E-08,(ZXP-ZXM))
         GO TO 3
3     CONTINUE
      END DO
   END DO
END DO
!
END SUBROUTINE PINTER
