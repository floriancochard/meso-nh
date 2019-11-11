!-----------------------------------------------------------------
!--------------- special set of characters for SCCS information
!-----------------------------------------------------------------
!      @(#) Lib:/opt/local/MESONH/sources/modd/s.modd_field1_cv2d.f90, Version:1.9, Date:98/10/01, Last modified:98/06/04
!-----------------------------------------------------------------
!     #######################
      MODULE MODD_FIELD1_CV2D
!     #######################
!
!!****  *MODD_FIELD1_CV2D* - declaration of arrays for prognostic variables
!                            in case of vertical sections
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     arrays holding prognostic variables in vertical planes
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_FIELD1_CV2D), to appear in 1994 
!!
!!    AUTHOR
!!    ------
!!      JD    "LA"
!!
!!    MODIFICATIONS
!!    -------------
!!      Original          05/05/94                      
!!      Updated   PM      17/11/94  
!!                (Stein) 08/03/95 Change the historical variables
!!                (Stein) 25/07/97 AChange the pressure variables
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE

REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XUMCV,XVMCV,XWMCV ! U,V,W at time t-dt
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XUTCV,XVTCV,XWTCV ! U,V,W at time t   
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XRUSCV,XRVSCV,XRWSCV ! Source of 
                                                   ! (rho U), (rho V), (rho W) 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XTHMCV    ! theta at time t-dt
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XTHTCV    ! theta at time t
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XRTHSCV   ! Source of (rho theta)
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XTKEMCV   ! Kinetic energy at time t-dt 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XTKETCV   ! Kinetic energy at time t
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XRTKESCV  ! Source of kinetic energy
                                                   ! (rho e)
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XPABSMCV  ! Pressure variable 
                                                   ! at time t-dt
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XPABSTCV  ! Pressure variable 
                                                   ! at time t
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE :: XRMCV   ! Moist variables
                                                   ! at time t-dt
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE :: XRTCV   ! Moist variables 
                                                   ! at time t
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE :: XRRSCV  ! Source of Moist variables
                                                   ! (rho Rn) 
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE :: XSVMCV  ! Additionnal scalar
                                                   ! variables at time t-deltat
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE :: XSVTCV  ! Additionnal scalar
                                                   ! variables at time t
REAL,DIMENSION(:,:,:), ALLOCATABLE,SAVE :: XRSVSCV ! Source of Additionnal scal.
                                                   !  variables (rho Sn.) 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XULMCV 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XULTCV 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XVTMCV 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XVTTCV 

REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XWORKCV 

REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XLSUMCV   ! Larger scale fields at
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XLSVMCV   ! time t-deltat for
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XLSWMCV   ! U,V,W,TH and Rv
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XLSTHMCV   
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XLSRVMCV   

REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XULMWMUCV  ! U component for UW 
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XULTWTUCV  ! vectors plot
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XULMWMWCV  ! W component for UW
REAL,DIMENSION(:,:), ALLOCATABLE,SAVE :: XULTWTWCV  ! vectors plot
!
END MODULE MODD_FIELD1_CV2D
