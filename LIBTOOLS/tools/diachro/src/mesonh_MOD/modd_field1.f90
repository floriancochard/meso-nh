!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ###################
      MODULE MODD_FIELD1
!     ###################
!
!!****  *MODD_FIELD1* - declaration of prognostic variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     prognostic variables. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_FIELDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       05/05/94                      
!!      Modifications  03/01/95  (Lafore)  To add the dry mass variables Md  
!!                     09/03/95  (Stein)   eliminate R from the progn. var                    
!!                     15/03/95  (Stein)   add EPS variable
!!      Modifications  21/03/95  (Carriere) To add the subgrid condensation 
!!                                           related parameters
!!                     01/03/96  (J. Stein) add the cloud fraction
!!                     10/10/96  (J. Stein) add XSRCM and XSRCT
!!                     11/04/96  (J.-P. Pinty) add the ice concentration
!!                     25/07/97  (J. Stein) Change the variable pressure
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XUM,XVM,XWM ! U,V,W  at time t-dt
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XUT,XVT,XWT ! U,V,W  at time t
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XRUS,XRVS,XRWS ! Source of (rho U),
                                                     ! (rho V), (rho w) 
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XTHM     ! (rho theta) at time t-dt
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XTHT     ! (rho theta) at time t
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XRTHS    ! Source of (rho theta)
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XTKEM    ! Kinetic energy 
                                                     ! at time t-dt
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XTKET    ! Kinetic energy
                                                     ! at time t
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XRTKES   ! Source of kinetic energy
                                                     ! (rho e)
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XEPSM    ! Dissipation of TKE 
                                                     ! (eps) at time t-dt
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XEPST    ! Dissipation of TKE    
                                                     ! (eps) at time t
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XREPSS   ! Source of Dissipation  
                                                     ! of TKE (rho eps)
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XPABSM   ! absolute pressure at
                                                     ! time t-dt
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XPABST   ! absolute pressure at
                                                     ! time t
REAL,SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: XRM    ! Moist variables 
                                                     ! at time t-dt
REAL,SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: XRT    ! Moist variables (rho Rn) 
                                                     ! at time t
REAL,SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: XRRS   ! Source of Moist variables
                                                     ! (rho Rn) 
REAL,SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: XSVM   ! Additionnal scalar
                                                     ! variables at time t-dt
REAL,SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: XSVT   ! Additionnal scalar
                                                     ! variables at time t  
REAL,SAVE, DIMENSION(:,:,:,:), ALLOCATABLE :: XRSVS  ! Source of addi. scalar
                                                     !  variables (rho Sn.) 
REAL,SAVE                          ::   XDRYMASST    ! Mass of dry air Md
REAL,SAVE                          ::   XDRYMASSS    ! LS sources of Md
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XSRC     ! turbulent flux <s'Rc'>
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XSIGS    ! =sqrt(<s's'>) for the
                                                     ! Subgrid Condensation
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XCLDFR   ! cloud fraction
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XSRCM    ! turbulent flux <s'Rc'>
                                                     ! at t- delta t
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XSRCT    ! turbulent flux <s'Rc'>
                                                     ! at t
REAL,SAVE, DIMENSION(:,:,:), ALLOCATABLE :: XCIT     ! Pristine ice concentration
!
END MODULE MODD_FIELD1
