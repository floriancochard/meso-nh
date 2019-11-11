!!   ########################
MODULE MODE_BLOWSNW_SURF
!!   ########################
!!

  USE MODD_BLOWSNW_SURF
  USE MODD_CSTS
  USE MODD_BLOWSNW_n


  IMPLICIT NONE
  PUBLIC

CONTAINS

SUBROUTINE SNOWMOMENT2SIZE(       &
     PSVT                         & !I [XX/m3] input scalar variables (moment of distribution)
     , PBETA1D                    & !O [m] scale factor of blowing snow particles distribution (gamma distribution)
     , PRG1D                      & !O [m] number mean radius of blowing snow particles distribution
     )
  !!   ############################################################
  !!
  !!
  !!    PURPOSE
  !!    -------
  !!    Translate the two moments M0 and M3 into
  !!    Values which can be understood more easily (R, beta)
  !!    At this point, M3 is in kg_{snow}/m3, M0 in #/m3
  !!
  !!   
  !!    REFERENCE
  !!    ---------
  !!    Based on mode_dst_surf.f90 (Tulet et al)
  !!
  !!    AUTHOR
  !!    ------
  !!    Vincent Vionnet
  !!
  !!    MODIFICATIONS
  !!    -------------
  !!
  !!    EXTERNAL
  !!    --------
  !!    None
  !!
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  USE PARKIND1  ,ONLY : JPRB

  IMPLICIT NONE
  !!
  !-------------------------------------------------------------------------------
  !
  !*       0.     DECLARATIONS
  !               ------------
  !
  !*      0.1    declarations of arguments
  !
  !INPUT
  REAL,       DIMENSION(:,:,:),  INTENT(IN)     :: PSVT      !I [ __ /m3] moments in surface units
  
  !OUTPUT
  REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PBETA1D   !O [m] scale factor deviation
  REAL,       DIMENSION(:,:),  OPTIONAL, INTENT(OUT)     :: PRG1D     !O [m] number median diameter
  !
  !
  !*      0.2    declarations local variables
  !
  REAL,DIMENSION(:,:,:), ALLOCATABLE  :: ZSV                 ! [snow particles moment concentration]
  REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZBETA               ! [-] standard deviation
  REAL,DIMENSION(:,:),   ALLOCATABLE  :: ZRG                 ! [um] number median diameter
  INTEGER                             :: JK   !Loop index 
  REAL,  DIMENSION(2)   :: ZPMIN
  REAL                  :: ZRGMIN
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_SURF:SNOWMOMENT2SIZE',0,ZHOOK_HANDLE)

  ALLOCATE (ZBETA(SIZE(PSVT,1),SIZE(PSVT,2)))
  ALLOCATE (ZRG(SIZE(PSVT,1),SIZE(PSVT,2)))
  ALLOCATE (ZSV(SIZE(PSVT,1), SIZE(PSVT,2),SIZE(PSVT,3)))

  !Save the moments in a local array
  ZSV(:,:,:) = MAX(PSVT(:,:,:), 1E-80)

!Get minimum values possible
  ZPMIN(1) = XN0MIN_SNW
  ZRGMIN   = XINIRADIUS_SNW
  ZPMIN(2) = 4*XPI*XRHOLI/3*(ZRGMIN/XEMIALPHA_SNW)**3.*(XEMIALPHA_SNW+2)*(XEMIALPHA_SNW+1)*XEMIALPHA_SNW*XN0MIN_SNW
  ZSV(:,:,1)=MAX(PSVT(:,:,1),ZPMIN(1))
  ZSV(:,:,2)=MAX(PSVT(:,:,2),ZPMIN(2))  
 
  DO JK=1,SIZE(PSVT,2)
    ZBETA(:,JK)=((3*ZSV(:,JK,2))/(4*XPI*XRHOLI*(XEMIALPHA_SNW+2)*(XEMIALPHA_SNW+1)*XEMIALPHA_SNW*ZSV(:,JK,1)))**(1./3.)

    ZRG(:,JK)  = ZBETA(:,JK)*XEMIALPHA_SNW

   END DO

  !Give the beta-values to the passed array
     IF(PRESENT(PBETA1D))THEN
        PBETA1D(:,:) = ZBETA(:,:)
     ENDIF

    !Get the number median radius
    IF(PRESENT(PRG1D))THEN
       PRG1D(:,:)= ZRG(:,:)
    ENDIF

    DEALLOCATE(ZSV)
    DEALLOCATE(ZRG)
    DEALLOCATE(ZBETA)

IF (LHOOK) CALL DR_HOOK('MODE_BLOWSNW_SURF:SNOWMOMENT2SIZE',1,ZHOOK_HANDLE)


END SUBROUTINE SNOWMOMENT2SIZE  

END MODULE MODE_BLOWSNW_SURF
