!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######################
      MODULE MODI_RELAXATION 
!     ######################
!
INTERFACE
!
      SUBROUTINE RELAXATION( OVE_RELAX,OVE_RELAX_GRD,                          &
                             OHORELAX_UVWTH,OHORELAX_RV,OHORELAX_RC,           &
                             OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,  &
                             OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,             &
                             OHORELAX_SVC2R2,OHORELAX_SVC1R3,                  &
                             OHORELAX_SVELEC,OHORELAX_SVLG,                    &
                             OHORELAX_SVCHEM,OHORELAX_SVCHIC, OHORELAX_SVAER,  &
                             OHORELAX_SVDST, OHORELAX_SVSLT, OHORELAX_SVPP,    &
                             OHORELAX_SVCS,OHORELAX_SVSNW,                     &
#ifdef MNH_FOREFIRE
                             OHORELAX_SVFF,                                    &
#endif
                             KTCOUNT,KRR,KSV,PTSTEP,PRHODJ,                    &
                             PUT, PVT, PWT, PTHT, PRT, PSVT, PTKET,            &
                             PLSUM, PLSVM, PLSWM, PLSTHM,                      &
                             PLBXUM, PLBXVM, PLBXWM, PLBXTHM, PLBXRM, PLBXSVM, &
                             PLBXTKEM,                                         &
                             PLBYUM, PLBYVM, PLBYWM, PLBYTHM, PLBYRM, PLBYSVM, &
                             PLBYTKEM,                                         &
                             KALBOT, PALK, PALKW,                              &
                             KALBAS, PALKBAS, PALKWBAS,                        &
                             OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX,         &
                             KRIMX,KRIMY,                                      &
                             PRUS,PRVS,PRWS,PRTHS,PRRS,PRSVS,PRTKES            )
!
LOGICAL,                  INTENT(IN)    :: OVE_RELAX  ! general switch for 
                                      ! vertical relaxation 
LOGICAL,                  INTENT(IN)    :: OVE_RELAX_GRD  ! general switch for 
                                      ! vertical relaxation (ground)                                      
LOGICAL,            INTENT(IN) :: OHORELAX_UVWTH  ! switch for the 
                       ! horizontal relaxation for U,V,W,TH
LOGICAL,            INTENT(IN) :: OHORELAX_RV     ! switch for the 
                       ! horizontal relaxation for Rv
LOGICAL,            INTENT(IN) :: OHORELAX_RC     ! switch for the 
                       ! horizontal relaxation for Rc
LOGICAL,            INTENT(IN) :: OHORELAX_RR     ! switch for the 
                       ! horizontal relaxation for Rr
LOGICAL,            INTENT(IN) :: OHORELAX_RI     ! switch for the 
                       ! horizontal relaxation for Ri
LOGICAL,            INTENT(IN) :: OHORELAX_RS     ! switch for the 
                       ! horizontal relaxation for Rs
LOGICAL,            INTENT(IN) :: OHORELAX_RG     ! switch for the 
                       ! horizontal relaxation for Rg
LOGICAL,            INTENT(IN) :: OHORELAX_RH     ! switch for the 
                       ! horizontal relaxation for Rh
LOGICAL,            INTENT(IN) :: OHORELAX_TKE    ! switch for the 
                       ! horizontal relaxation for tke
LOGICAL,DIMENSION(:),INTENT(IN):: OHORELAX_SV     ! switch for the 
                       ! horizontal relaxation for sv
LOGICAL,             INTENT(IN):: OHORELAX_SVC2R2 ! switch for the 
                       ! horizontal relaxation for c2r2 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVC1R3 ! switch for the 
                       ! horizontal relaxation for c1r3 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVELEC ! switch for the 
                       ! horizontal relaxation for elec variables
LOGICAL,             INTENT(IN):: OHORELAX_SVLG   ! switch for the 
                       ! horizontal relaxation for lg variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHEM ! switch for the 
                       ! horizontal relaxation for chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHIC ! switch for the 
                       ! horizontal relaxation for ice chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVAER  ! switch for the 
                       ! horizontal relaxation for aer variables
LOGICAL,             INTENT(IN):: OHORELAX_SVDST  ! switch for the 
                       ! horizontal relaxation for dst variables
LOGICAL,             INTENT(IN):: OHORELAX_SVSLT  ! switch for the 
                       ! horizontal relaxation for slt variables
LOGICAL,             INTENT(IN):: OHORELAX_SVPP   ! switch for the 
                       ! horizontal relaxation for passive scalar
#ifdef MNH_FOREFIRE
LOGICAL,             INTENT(IN):: OHORELAX_SVFF   ! switch for the 
                       ! horizontal relaxation for ForeFire variables 
#endif
LOGICAL,             INTENT(IN):: OHORELAX_SVCS   ! switch for the 
                       ! horizontal relaxation for conditional sampling
LOGICAL,             INTENT(IN):: OHORELAX_SVSNW  ! switch for the
                       ! horizontal relaxation for blowing snow variables
INTEGER,                  INTENT(IN)    :: KTCOUNT! Temporal loop counter       
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
INTEGER,                  INTENT(IN)    :: KRR    ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV    ! Number of scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! effective dry rho * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT, PWT    ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT             !  
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT              !   
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT             !         
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTKET
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSUM, PLSVM     !   Large Scale
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSWM, PLSTHM    !    variables
                                                            !     at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBXUM, PLBXVM, PLBXWM! Lateral
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBXTHM         ! boundarie fields
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBXRM          !    at  t-dt
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBXSVM         !    along the x         
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBXTKEM        !    direction
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBYUM, PLBYVM, PLBYWM! Lateral
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBYTHM         ! boundarie fields
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBYRM          !    at  t-dt
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBYSVM         !    along the y         
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBYTKEM        !    direction

!
INTEGER,                  INTENT(IN)    :: KALBOT !Vertical index corresponding
                                      ! to the absorbing layer base
REAL, DIMENSION(:),       INTENT(IN)    :: PALK   ! Function of the absorbing 
                                      ! layer damping coefficient defined for 
                                      ! u,v,and theta
REAL, DIMENSION(:),       INTENT(IN)    :: PALKW  ! Idem but defined for w    
!
INTEGER,                  INTENT(IN)    :: KALBAS !Vertical index corresponding
                                      ! to the absorbing layer base
REAL, DIMENSION(:),       INTENT(IN)    :: PALKBAS! Function of the absorbing 
                                      ! layer damping coefficient defined for 
                                      ! u,v,and theta
REAL, DIMENSION(:),       INTENT(IN)    :: PALKWBAS  ! Idem but defined for w    
!
LOGICAL, DIMENSION(:,:),  INTENT(IN)    :: OMASK_RELAX ! Mask for the locations
                                      ! where lateral relax. must be performed
!
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKURELAX ! Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKVRELAX ! coefficients for the
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKWRELAX ! u, v and mass locations
INTEGER,                  INTENT(IN)    :: KRIMX, KRIMY ! size of the hor. 
                                                        ! relaxation area
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS ! Source terms
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS            !
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS             !
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS            !
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTKES
!
END SUBROUTINE RELAXATION
!
END INTERFACE
!
END MODULE MODI_RELAXATION
!     ######spl
      SUBROUTINE RELAXATION( OVE_RELAX,OVE_RELAX_GRD,                          &
                             OHORELAX_UVWTH,OHORELAX_RV,OHORELAX_RC,           &
                             OHORELAX_RR,OHORELAX_RI,OHORELAX_RS,OHORELAX_RG,  &
                             OHORELAX_RH,OHORELAX_TKE,OHORELAX_SV,             &
                             OHORELAX_SVC2R2,OHORELAX_SVC1R3,                  &
                             OHORELAX_SVELEC,OHORELAX_SVLG,                    &
                             OHORELAX_SVCHEM,OHORELAX_SVCHIC, OHORELAX_SVAER,  &
                             OHORELAX_SVDST, OHORELAX_SVSLT, OHORELAX_SVPP,    &
                             OHORELAX_SVCS,OHORELAX_SVSNW,                     &
#ifdef MNH_FOREFIRE
                             OHORELAX_SVFF,                                    &
#endif
                             KTCOUNT,KRR,KSV,PTSTEP,PRHODJ,                    &
                             PUT, PVT, PWT, PTHT, PRT, PSVT, PTKET,            &
                             PLSUM, PLSVM, PLSWM, PLSTHM,                      &
                             PLBXUM, PLBXVM, PLBXWM, PLBXTHM, PLBXRM, PLBXSVM, &
                             PLBXTKEM,                                         &
                             PLBYUM, PLBYVM, PLBYWM, PLBYTHM, PLBYRM, PLBYSVM, &
                             PLBYTKEM,                                         &
                             KALBOT, PALK, PALKW,                              &
                             KALBAS, PALKBAS, PALKWBAS,                        &
                             OMASK_RELAX,PKURELAX, PKVRELAX, PKWRELAX,         &
                             KRIMX,KRIMY,                                      &
                             PRUS,PRVS,PRWS,PRTHS,PRRS,PRSVS,PRTKES            )
!     ##########################################################################
!
!!****  *RELAXATION * - routine to apply a Rayleigh damping at the top  
!!                       and in the outermost verticals of the domain
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to damp the u, v, w, and theta fields  
!!    near the top of the domain and in the outermost vertical planes
!!    (vapor field added) to prevent spurious reflections on the top and on
!!    the side boundaries.
!!
!!
!!**  METHOD
!!    ------
!!      Implicit Rayleigh damping for the absorbing layer 
!!    upper part:
!!                                     -Kal(zhat)
!!    -Kal(zhat)(phi(t+dt)-LSphi)= ------------------ (phi(t-dt)-LSphi)
!!                                  1 + 2dt Kal(zhat)
!!    lateral part:
!!                                  -Krelax
!!    -Krelax(phi(t+dt)-LSphi)= ------------------ (phi(t-dt)-LSphi)
!!                                1 + 2dt Krelax
!!
!!     phi being u, v, w, theta or qv
!!     LSphi are the corresponding Larger Scale values or to the lateral 
!!     boundaries fields.
!!     The different sources terms are stored for the budget computations.
!!
!!
!!    EXTERNAL
!!    --------
!!      MXM,MYM,MZM         : mean Shuman operators in the x,y,z directions
!!      BUDGET              : Stores the different budget components
!!      GET_INDICE_ll       : get physical sub-domain bounds
!!      GET_GLOBALDIMS_ll   : get physical global domain sizes
!!      GET_INTERSECTION_ll : get the indices of the intersection inside the
!!                            extended global domain
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!       Module MODD_PARAMETERS : JPVEXT
!!
!!       Module MODD_CONF       : CCONF 
!!
!!       Module MODD_BUDGET     : NBUMOD,LBU_ENABLE,
!!                                LBU_RU,LBU_RV,LBU_RW,LBU_RTH,LBU_RRV,
!!                                LBU_RRC,LBU_RRR,LBU_RRI,LBU_RRH,LBU_RRG,
!!                                LBU_RTKE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine RELAXATION )
!!
!!    AUTHOR
!!    ------
!!      E. Richard       * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         10/10/94 
!!      Modification
!!                 03/11/94 (J.Stein)     change the signs of the Rayleigh 
!!                 07/12/94 (J.-P. Pinty) merge ABS_LAYER and the former
!!                                              RELAXATION routines
!!                 07/12/94 (J.Stein)     add a 3D mask 
!!                 20/.3/95 (J.Stein)     remove R from the historical var.
!!                 01/04/95 (Ph. Hereil J. Nicolau) add the budget computation
!!                 16/10/95 (J. Stein)    change the budget calls 
!!                 19/12/96 (J.-P. Pinty) update the budget calls 
!!                 21/10/97 (J. Stein)    separate the vertical relaxation
!!                                        towards the LS fields from the lateral
!!                                        boundaries treatment
!!                 24/08/98 (Jabouille)   parallelization 
!!                 06/11/02 (V. Masson)  update the budget calls 
!!                 05/2006               Remove EPS
!!                 06/2011 (M.Chong)     Case of ELEC
!!                 11/2011 (C.Lac)       Adaptation to FIT temporal scheme
!!                 J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_CONF 
USE MODD_BUDGET 
USE MODD_NSV, ONLY : NSV_ELECBEG, NSV_ELECEND      
USE MODD_ELEC_DESCR, ONLY: LRELAX2FW_ION          
!
USE MODE_ll
!
USE MODI_SHUMAN     
USE MODI_BUDGET
USE MODE_EXTRAPOL
!
USE MODE_MPPDB
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
LOGICAL,                  INTENT(IN)    :: OVE_RELAX  ! general switch for 
                                      ! vertical relaxation 
LOGICAL,                  INTENT(IN)    :: OVE_RELAX_GRD  ! general switch for 
                                      ! vertical relaxation (ground)                                      
LOGICAL,            INTENT(IN) :: OHORELAX_UVWTH  ! switch for the 
                       ! horizontal relaxation for U,V,W,TH
LOGICAL,            INTENT(IN) :: OHORELAX_RV     ! switch for the 
                       ! horizontal relaxation for Rv
LOGICAL,            INTENT(IN) :: OHORELAX_RC     ! switch for the 
                       ! horizontal relaxation for Rc
LOGICAL,            INTENT(IN) :: OHORELAX_RR     ! switch for the 
                       ! horizontal relaxation for Rr
LOGICAL,            INTENT(IN) :: OHORELAX_RI     ! switch for the 
                       ! horizontal relaxation for Ri
LOGICAL,            INTENT(IN) :: OHORELAX_RS     ! switch for the 
                       ! horizontal relaxation for Rs
LOGICAL,            INTENT(IN) :: OHORELAX_RG     ! switch for the 
                       ! horizontal relaxation for Rg
LOGICAL,            INTENT(IN) :: OHORELAX_RH     ! switch for the 
                       ! horizontal relaxation for Rh
LOGICAL,            INTENT(IN) :: OHORELAX_TKE    ! switch for the 
                       ! horizontal relaxation for tke
LOGICAL,DIMENSION(:),INTENT(IN):: OHORELAX_SV     ! switch for the 
                       ! horizontal relaxation for sv
LOGICAL,             INTENT(IN):: OHORELAX_SVC2R2 ! switch for the 
                       ! horizontal relaxation for c2r2 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVC1R3 ! switch for the 
                       ! horizontal relaxation for c1r3 variables
LOGICAL,             INTENT(IN):: OHORELAX_SVELEC ! switch for the 
                       ! horizontal relaxation for elec variables
LOGICAL,             INTENT(IN):: OHORELAX_SVLG   ! switch for the 
                       ! horizontal relaxation for lg variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHEM ! switch for the 
                       ! horizontal relaxation for chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVCHIC ! switch for the
                       ! horizontal relaxation for ice chem variables
LOGICAL,             INTENT(IN):: OHORELAX_SVAER  ! switch for the 
                       ! horizontal relaxation for aer variables
LOGICAL,             INTENT(IN):: OHORELAX_SVDST  ! switch for the 
                       ! horizontal relaxation for dst variables
LOGICAL,             INTENT(IN):: OHORELAX_SVSLT  ! switch for the 
                       ! horizontal relaxation for slt variables
LOGICAL,             INTENT(IN):: OHORELAX_SVPP   ! switch for the 
                       ! horizontal relaxation for passive scalar
#ifdef MNH_FOREFIRE 
LOGICAL,             INTENT(IN):: OHORELAX_SVFF   ! switch for the 
                       ! horizontal relaxation for ForeFire variables 
#endif
LOGICAL,             INTENT(IN):: OHORELAX_SVCS   ! switch for the 
                       ! horizontal relaxation for conditional sampling
LOGICAL,             INTENT(IN):: OHORELAX_SVSNW  ! switch for the
                       ! horizontal relaxation for blowing snow variables                       
INTEGER,                  INTENT(IN)    :: KTCOUNT! Temporal loop counter       
REAL,                     INTENT(IN)    :: PTSTEP ! Time step
INTEGER,                  INTENT(IN)    :: KRR    ! Number of moist variables
INTEGER,                  INTENT(IN)    :: KSV    ! Number of scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ ! effective dry rho * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PUT, PVT, PWT    ! Variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT             !  
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT              !  
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT             !         
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTKET
!
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSUM, PLSVM     !   Large Scale
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLSWM, PLSTHM    !    variables
                                                            !     at t-dt
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBXUM, PLBXVM, PLBXWM! Lateral
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBXTHM         ! boundarie fields
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBXRM          !    at  t-dt
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBXSVM         !    along the x         
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBXTKEM        !    direction
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBYUM, PLBYVM, PLBYWM! Lateral
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBYTHM         ! boundarie fields
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBYRM          !    at  t-dt
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PLBYSVM         !    along the y         
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PLBYTKEM        !    direction

!
INTEGER,                  INTENT(IN)    :: KALBOT !Vertical index corresponding
                                      ! to the absorbing layer base
REAL, DIMENSION(:),       INTENT(IN)    :: PALK   ! Function of the absorbing 
                                      ! layer damping coefficient defined for 
                                      ! u,v,and theta
REAL, DIMENSION(:),       INTENT(IN)    :: PALKW  ! Idem but defined for w    
!
INTEGER,                  INTENT(IN)    :: KALBAS !Vertical index corresponding
                                      ! to the absorbing layer base
REAL, DIMENSION(:),       INTENT(IN)    :: PALKBAS! Function of the absorbing 
                                      ! layer damping coefficient defined for 
                                      ! u,v,and theta
REAL, DIMENSION(:),       INTENT(IN)    :: PALKWBAS  ! Idem but defined for w    
!
LOGICAL, DIMENSION(:,:),  INTENT(IN)    :: OMASK_RELAX ! Mask for the locations
                                      ! where lateral relax. must be performed
!
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKURELAX ! Horizontal relaxation
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKVRELAX ! coefficients for the
REAL, DIMENSION(:,:),     INTENT(IN)    :: PKWRELAX ! u, v and mass locations
INTEGER,                  INTENT(IN)    :: KRIMX, KRIMY ! size of the hor. 
                                                        ! relaxation area
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS ! Source terms
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS            !
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS             !
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS            !
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTKES
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK,JRR,JSV        ! Loop index 
INTEGER :: IIU_ll,IJU_ll     ! horizontal size of the extended global domain  
INTEGER :: IIB,IJB,IIE,IJE   ! Bounds of the physical sub-domain             
INTEGER :: IKU,IKE           ! size of the array in the z direction and index
                             ! of the last inner mass point
INTEGER    :: IIBINT,IJBINT  ! Corner of relaxation domain 
INTEGER    :: IIEINT,IJEINT  ! global indices translated to local indices
INTEGER    :: IINFO_ll        ! return code of parallel routine
INTEGER    :: IDIMLB         ! size of the LB field
!
REAL, DIMENSION(SIZE(PALK))          :: ZKV  ! Function of the upper absorbing 
                             ! layer damping coefficient for u,v,theta and qv
REAL, DIMENSION(SIZE(PALK))          :: ZKVW !Idem but for w 
!
REAL, DIMENSION(SIZE(PALKBAS))          :: ZKVBAS  ! Function of the upper absorbing 
                             ! layer damping coefficient for u,v,theta and qv
REAL, DIMENSION(SIZE(PALKBAS))          :: ZKVWBAS !Idem but for w
!
REAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: ZKHU,ZKHV,ZKHW,       &
                             ! Function of the lateral absorbing layer damping 
                             ! for u,v and mass points respectively
                                                     ZRHODJU,ZRHODJV,ZRHODJW, &
                             ! averages along x,y,z of the PRHODJ field
                                                     ZWORK
                             ! work array used to expand the LB fields
LOGICAL, DIMENSION(SIZE(PUT,1),SIZE(PUT,2),SIZE(PUT,3)) :: GMASK3D_RELAX ! 3D
                             ! mask for hor. relax.
LOGICAL, DIMENSION(7) :: GHORELAXR ! local array of logical
#ifdef MNH_FOREFIRE
LOGICAL, DIMENSION(13) :: GHORELAXSV! local array of logical
#else
LOGICAL, DIMENSION(12) :: GHORELAXSV! local array of logical
#endif
!  
!-------------------------------------------------------------------------------
!
!
!*       1.     PRELIMINARIES
!	        -------------
IKU=SIZE(PUT,3)
IKE=IKU-JPVEXT
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
CALL GET_GLOBALDIMS_ll(IIU_ll,IJU_ll)
IIU_ll=IIU_ll+2*JPHEXT
IJU_ll=IJU_ll+2*JPHEXT
!
ZRHODJU(:,:,:) = MXM(PRHODJ)
ZRHODJV(:,:,:) = MYM(PRHODJ)
ZRHODJW(:,:,:) = MZM(1,IKU,1,PRHODJ)
!
GHORELAXR(1) = OHORELAX_RV
GHORELAXR(2) = OHORELAX_RC
GHORELAXR(3) = OHORELAX_RR
GHORELAXR(4) = OHORELAX_RI
GHORELAXR(5) = OHORELAX_RS
GHORELAXR(6) = OHORELAX_RG
GHORELAXR(7) = OHORELAX_RH
!
GHORELAXSV(1) = OHORELAX_SVC2R2
GHORELAXSV(2) = OHORELAX_SVC1R3
GHORELAXSV(3) = OHORELAX_SVELEC
GHORELAXSV(4) = OHORELAX_SVLG
GHORELAXSV(5) = OHORELAX_SVCHEM
GHORELAXSV(6) = OHORELAX_SVAER
GHORELAXSV(7) = OHORELAX_SVDST
GHORELAXSV(8) = OHORELAX_SVSLT
GHORELAXSV(9) = OHORELAX_SVPP
GHORELAXSV(10) = OHORELAX_SVCS
GHORELAXSV(11) = OHORELAX_SVCHIC
GHORELAXSV(12) = OHORELAX_SVSNW
#ifdef MNH_FOREFIRE
GHORELAXSV(13) = OHORELAX_SVFF
#endif
!-------------------------------------------------------------------------------
!
!*       2.     RELAXATION IN THE UPPER LAYERS
!	        ------------------------------
!
IF(OVE_RELAX) THEN
!
!*       2.1     SET THE TOP-LEVEL DAMPING COEF. (UPSTREAM OR LEAPFROG)
!
! IF (KTCOUNT.EQ.1) THEN
    ZKV(:)  = PALK(:) /(1. - PTSTEP * PALK(:))
    ZKVW(:) = PALKW(:)/(1. - PTSTEP * PALKW(:))
! ELSE
!   ZKV(:)  = PALK(:)
!   ZKVW(:) = PALKW(:)
! ENDIF
!
!
!*       2.2     APPLIES THE DAMPING IN THE UPPERMOST LEVELS
!
  DO JK = KALBOT, IKE+1
!
    PRUS(:,:,JK)  = PRUS(:,:,JK)  - ZKV(JK)  *(PUT(:,:,JK)  -PLSUM(:,:,JK)  )&
                    * ZRHODJU(:,:,JK)
!
    PRVS(:,:,JK)  = PRVS(:,:,JK)  - ZKV(JK)  *(PVT(:,:,JK)  -PLSVM(:,:,JK)  )&
                    * ZRHODJV(:,:,JK)
!
    PRWS(:,:,JK)  = PRWS(:,:,JK)  - ZKVW(JK) *(PWT(:,:,JK)  -PLSWM(:,:,JK)  )&
                    * ZRHODJW(:,:,JK)
!
    PRTHS(:,:,JK) = PRTHS(:,:,JK) - ZKV(JK)  *(PTHT(:,:,JK) -PLSTHM(:,:,JK) )&
                    * PRHODJ(:,:,JK)
!
  END DO  
!
END IF   
!-------------------------------------------------------------------------------
!
!*       2.bis     RELAXATION IN THE GROUND LAYERS
!	        ------------------------------
!
IF(OVE_RELAX_GRD) THEN
!
!*       2.bis.1     SET THE GROUND-LEVEL DAMPING COEF. (UPSTREAM OR LEAPFROG)
!
  IF ((KTCOUNT.EQ.1) .AND. (CCONF.EQ.'START') ) THEN
    ZKVBAS(:)  = PALKBAS(:) /(1. - PTSTEP * PALKBAS(:))
    ZKVWBAS(:) = PALKWBAS(:)/(1. - PTSTEP * PALKWBAS(:))
  ELSE
    ZKVBAS(:)  = PALKBAS(:)
    ZKVWBAS(:) = PALKWBAS(:)
  ENDIF
!
!
!*       2.bis.2     APPLIES THE DAMPING IN THE GROUND LEVELS
!
  DO JK = 1,KALBAS
!
    PRUS(:,:,JK)  = PRUS(:,:,JK)  - ZKVBAS(JK)  *(PUT(:,:,JK)  -PLSUM(:,:,JK)  )&
                    * ZRHODJU(:,:,JK)
!
    PRVS(:,:,JK)  = PRVS(:,:,JK)  - ZKVBAS(JK)  *(PVT(:,:,JK)  -PLSVM(:,:,JK)  )&
                    * ZRHODJV(:,:,JK)
!
    PRWS(:,:,JK)  = PRWS(:,:,JK)  - ZKVWBAS(JK) *(PWT(:,:,JK)  -PLSWM(:,:,JK)  )&
                    * ZRHODJW(:,:,JK)
!
    PRTHS(:,:,JK) = PRTHS(:,:,JK) - ZKVBAS(JK)  *(PTHT(:,:,JK) -PLSTHM(:,:,JK) )&
                    * PRHODJ(:,:,JK)
!
  END DO  
!
END IF   

!  
!-------------------------------------------------------------------------------
!
!*       3.     RELAXATION IN THE OUTERMOST VERTICAL PLANES
!	        -------------------------------------------
!
!
!*       3.1    SET THE RIM ZONE DAMPING COEF. (UPSTREAM OR LEAPFROG)
!
IF ( ANY(GHORELAXR) .OR. ANY(GHORELAXSV) .OR. ANY(OHORELAX_SV) &
                    .OR. OHORELAX_UVWTH  .OR. OHORELAX_TKE     ) THEN 
! IF (KTCOUNT.EQ.1)  THEN
    DO JK=1,IKU
      ZKHU(:,:,JK) = PKURELAX(:,:) /(1. - PTSTEP * PKURELAX(:,:))
      ZKHV(:,:,JK) = PKVRELAX(:,:) /(1. - PTSTEP * PKVRELAX(:,:)) 
      ZKHW(:,:,JK) = PKWRELAX(:,:) /(1. - PTSTEP * PKWRELAX(:,:))
    END DO
! ELSE
!   DO JK=1,IKU
!     ZKHU(:,:,JK) = PKURELAX(:,:) 
!     ZKHV(:,:,JK) = PKVRELAX(:,:)
!     ZKHW(:,:,JK) = PKWRELAX(:,:)
!   END DO
!  END IF
ENDIF
!
!
!*       3.2    APPLIES THE DAMPING NEAR THE LATERAL BOUNDARIES
!
  ZWORK=0.
!
! special treatment is needed to expand LB for U and V because of the C grid 
IF ( OHORELAX_UVWTH ) THEN
!
  DO JK=1,IKU
    GMASK3D_RELAX(:,:,JK)=OMASK_RELAX(:,:)
  END DO
!
  IDIMLB = SIZE(PLBXUM,1)
  IF ( IDIMLB /= 0) THEN
    CALL GET_INTERSECTION_ll (1,1,KRIMX+JPHEXT+1,IJU_ll,                 &   ! +2
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(2:IIEINT,:,:) = PLBXUM(1:IIEINT-1,:,:)
    END IF
    CALL GET_INTERSECTION_ll (IIU_ll-KRIMX-JPHEXT+1,1,IIU_ll,IJU_ll,       &          ! -KRIMX 
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(IIBINT:IIE+JPHEXT,:,:) =  PLBXUM(IDIMLB-(IIE+JPHEXT-IIBINT):IDIMLB,:,:) ! +1
    END IF
  ENDIF
!
  IDIMLB = SIZE(PLBYUM,2)
  IF ( IDIMLB /= 0) THEN
    CALL GET_INTERSECTION_ll (1,1,IIU_ll,KRIMY+JPHEXT,                 &    ! +1
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(:,1:IJEINT,:) = PLBYUM(:,1:IJEINT,:)
    END IF
    CALL GET_INTERSECTION_ll (1,IJU_ll-KRIMY-JPHEXT+1,IIU_ll,IJU_ll,   &  ! -KRIMY
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(:,IJBINT:IJE+JPHEXT,:) =  PLBYUM(:,IDIMLB-(IJE+JPHEXT-IJBINT):IDIMLB,:) ! +1
    END IF
  END IF
  !
  CALL MPPDB_CHECK3DM("before PRUS relax:ZWORK",PRECISION,ZWORK)
  WHERE (GMASK3D_RELAX)
    PRUS(:,:,:)  = PRUS(:,:,:) - ZKHU(:,:,:)*(PUT(:,:,:)-ZWORK(:,:,:)) &
                       * ZRHODJU(:,:,:)
  END WHERE
 CALL MPPDB_CHECK3DM("after  PRUS relax:ZWORK",PRECISION,ZWORK)
!
  IDIMLB =  SIZE(PLBXVM,1)
  IF ( IDIMLB /= 0) THEN
    CALL GET_INTERSECTION_ll (1,1,KRIMX+JPHEXT,IJU_ll,                 &     ! +1
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(1:IIEINT,:,:) = PLBXVM(1:IIEINT,:,:)
    END IF
    CALL GET_INTERSECTION_ll (IIU_ll-KRIMX-JPHEXT+1,1,IIU_ll,IJU_ll,       & ! -KRIMX   
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(IIBINT:IIE+JPHEXT,:,:) = PLBXVM(IDIMLB-(IIE+JPHEXT-IIBINT):IDIMLB,:,:)     ! +1
    END IF
  ENDIF
!
  IDIMLB =  SIZE(PLBYVM,2)
  IF ( IDIMLB /= 0) THEN
    CALL GET_INTERSECTION_ll (1,1,IIU_ll,KRIMY+JPHEXT+1,                 & ! +2
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(:,2:IJEINT,:) = PLBYVM(:,1:IJEINT-1,:)
    END IF
    CALL GET_INTERSECTION_ll (1,IJU_ll-KRIMY-JPHEXT+1,IIU_ll,IJU_ll,     & ! -KRIMY
                               IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
    IF ( IINFO_ll == 0 ) THEN
      ZWORK(:,IJBINT:IJE+JPHEXT,:) = PLBYVM(:,IDIMLB-(IJE+JPHEXT-IJBINT):IDIMLB,:) ! +1
    END IF
  ENDIF
  !
  WHERE (GMASK3D_RELAX)
    PRVS(:,:,:)  = PRVS(:,:,:) - ZKHV(:,:,:)*(PVT(:,:,:)-ZWORK(:,:,:)) &
                       * ZRHODJV(:,:,:)
  END WHERE
  !
  IF (SIZE(PLBXWM,1) > 0) CALL EXPAND_LBX (PLBXWM,ZWORK)
  IF (SIZE(PLBYWM,2) > 0) CALL EXPAND_LBY (PLBYWM,ZWORK)
  !
  WHERE (GMASK3D_RELAX)
    PRWS(:,:,:)  = PRWS(:,:,:) - ZKHW(:,:,:)*(PWT(:,:,:)-ZWORK(:,:,:)) &
                       * ZRHODJW(:,:,:)
  END WHERE
  !
  !
  IF (SIZE(PLBXTHM,1) > 0) CALL EXPAND_LBX (PLBXTHM,ZWORK)
  IF (SIZE(PLBYTHM,2) > 0) CALL EXPAND_LBY (PLBYTHM,ZWORK)
  !
  WHERE (GMASK3D_RELAX)
    PRTHS(:,:,:)  = PRTHS(:,:,:) - ZKHW(:,:,:)*(PTHT(:,:,:)-ZWORK(:,:,:)) &
                       * PRHODJ(:,:,:)
  END WHERE
END IF
!
DO JRR = 1,KRR
  IF ( GHORELAXR(JRR) ) THEN
    IF (SIZE(PLBXRM,1) > 0) CALL EXPAND_LBX (PLBXRM(:,:,:,JRR),ZWORK)
    IF (SIZE(PLBYRM,2) > 0) CALL EXPAND_LBY (PLBYRM(:,:,:,JRR),ZWORK)
    WHERE (GMASK3D_RELAX)
      PRRS(:,:,:,JRR)  = PRRS(:,:,:,JRR) - ZKHW(:,:,:)*(PRT(:,:,:,JRR)-ZWORK(:,:,:)) &
                         * PRHODJ(:,:,:)
    END WHERE
  END IF
END DO
!
IF ( OHORELAX_TKE ) THEN
  IF (SIZE(PLBXTKEM,1) > 0) CALL EXPAND_LBX (PLBXTKEM,ZWORK)
  IF (SIZE(PLBYTKEM,2) > 0) CALL EXPAND_LBY (PLBYTKEM,ZWORK)
  WHERE (GMASK3D_RELAX)
    PRTKES(:,:,:)  = PRTKES(:,:,:) - ZKHW(:,:,:)*(PTKET(:,:,:)-ZWORK(:,:,:)) &
                       * PRHODJ(:,:,:)
  END WHERE
END IF  
!
DO JSV=1,KSV
  IF ( .NOT. LRELAX2FW_ION .OR. (JSV .NE. NSV_ELECBEG .AND. JSV .NE. NSV_ELECEND)) THEN 
    IF ( OHORELAX_SV(JSV) ) THEN
      IF (SIZE(PLBXSVM,1) > 0) CALL EXPAND_LBX (PLBXSVM(:,:,:,JSV),ZWORK)
      IF (SIZE(PLBYSVM,2) > 0) CALL EXPAND_LBY (PLBYSVM(:,:,:,JSV),ZWORK)
      WHERE (GMASK3D_RELAX)
        PRSVS(:,:,:,JSV)  = PRSVS(:,:,:,JSV) - ZKHW(:,:,:)*(PSVT(:,:,:,JSV)-ZWORK(:,:,:)) &
                         * PRHODJ(:,:,:)
      END WHERE
    END IF
  END IF
END DO
!
!
!-------------------------------------------------------------------------------
!
!*       3.     STORES FIELDS IN BUDGET ARRAYS
!	        ------------------------------
!
CALL EXTRAPOL('W ', PRUS)
IF (LBUDGET_U) CALL BUDGET  (PRUS,1,'REL_BU_RU')
IF (LBUDGET_V) CALL BUDGET  (PRVS,2,'REL_BU_RV')
IF (LBUDGET_W) CALL BUDGET  (PRWS,3,'REL_BU_RW')
IF (LBUDGET_TH) CALL BUDGET (PRTHS,4,'REL_BU_RTH')
IF (LBUDGET_TKE) CALL BUDGET (PRTKES,5,'REL_BU_RTKE')
IF (LBUDGET_RV)  CALL BUDGET (PRRS(:,:,:,1),6,'REL_BU_RRV')
IF (LBUDGET_RC)  CALL BUDGET (PRRS(:,:,:,2),7,'REL_BU_RRC')
IF (LBUDGET_RR)  CALL BUDGET (PRRS(:,:,:,3),8,'REL_BU_RRR')
IF (LBUDGET_RI)  CALL BUDGET (PRRS(:,:,:,4),9,'REL_BU_RRI')
IF (LBUDGET_RS)  CALL BUDGET (PRRS(:,:,:,5),10,'REL_BU_RRS')
IF (LBUDGET_RG)  CALL BUDGET (PRRS(:,:,:,6),11,'REL_BU_RRG')
IF (LBUDGET_RH)  CALL BUDGET (PRRS(:,:,:,7),12,'REL_BU_RRH')
IF (LBUDGET_SV) THEN
  DO JSV=1,KSV 
    CALL BUDGET (PRSVS(:,:,:,JSV),JSV+12,'REL_BU_RSV')
  END DO
END IF
!
CONTAINS
!     ######################################
      SUBROUTINE EXPAND_LB (PLBX,PLBY,PWORK)
!     ######################################
!
!!****  *EXPAND_LB * - internal routine to expand the lateral boundarie 
!!                     (non horizontal velocity) fields in a complete 3D array
!!    AUTHOR
!!    ------
!!      P. Jabouille     * * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    14/09/98 
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PLBX,PLBY !Lateral boundarie fields
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PWORK  !work array used to 
                                                !expand the LB fields
IDIMLB = SIZE(PLBX,1)
IF ( IDIMLB /= 0) THEN
  CALL GET_INTERSECTION_ll (1,1,KRIMX+1,IJU_ll,                 &
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(1:IIEINT,:,:) = PLBX(1:IIEINT,:,:)
  END IF
  CALL GET_INTERSECTION_ll (IIU_ll-KRIMX,1,IIU_ll,IJU_ll,       &
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(IIBINT:IIE+1,:,:) = PLBX(IDIMLB-(IIE+1-IIBINT):IDIMLB,:,:)
  END IF
ENDIF
!
IDIMLB = SIZE(PLBY,2)
IF ( IDIMLB /= 0) THEN
  CALL GET_INTERSECTION_ll (1,1,IIU_ll,KRIMY+1                 ,&
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(:,1:IJEINT,:) = PLBY(:,1:IJEINT,:)
  END IF
  CALL GET_INTERSECTION_ll (1,IJU_ll-KRIMY,IIU_ll,IJU_ll,       &
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(:,IJBINT:IJE+1,:) = PLBY(:,IDIMLB-(IJE+1-IJBINT):IDIMLB,:)
  END IF
END IF

END SUBROUTINE EXPAND_LB
!
!     ######################################
      SUBROUTINE EXPAND_LBX (PLBX,PWORK)
!     ######################################
!
!!****  *EXPAND_LB * - internal routine to expand the lateral boundarie 
!!                     (non horizontal velocity) fields in a complete 3D array
!!    AUTHOR
!!    ------
!!      P. Jabouille     * * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    14/09/98 
!!      P. Tulet    06/08/07 - split X and Y LB 
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PLBX !Lateral boundarie fields
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PWORK  !work array used to 
                                                !expand the LB fields
IDIMLB = SIZE(PLBX,1)
IF ( IDIMLB /= 0) THEN
  CALL GET_INTERSECTION_ll (1,1,KRIMX+JPHEXT,IJU_ll,                 & ! +1
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(1:IIEINT,:,:) = PLBX(1:IIEINT,:,:)
  END IF
  CALL GET_INTERSECTION_ll (IIU_ll-KRIMX-JPHEXT+1,1,IIU_ll,IJU_ll,       & ! -KRIMX
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(IIBINT:IIE+JPHEXT,:,:) = PLBX(IDIMLB-(IIE+JPHEXT-IIBINT):IDIMLB,:,:)  ! +1
  END IF
ENDIF
!
END SUBROUTINE EXPAND_LBX
!
!     ######################################
      SUBROUTINE EXPAND_LBY (PLBY,PWORK)
!     ######################################
!
!!****  *EXPAND_LB * - internal routine to expand the lateral boundarie 
!!                     (non horizontal velocity) fields in a complete 3D array
!!    AUTHOR
!!    ------
!!      P. Jabouille     * * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    14/09/98 
!!      P. Tulet    06/09/07 - split X and Y LB 
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments :
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PLBY !Lateral boundarie fields
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PWORK  !work array used to 
                                                !expand the LB fields
IDIMLB = SIZE(PLBY,2)
IF ( IDIMLB /= 0) THEN
  CALL GET_INTERSECTION_ll (1,1,IIU_ll,KRIMY+JPHEXT                 ,& ! +1
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(:,1:IJEINT,:) = PLBY(:,1:IJEINT,:)
  END IF
  CALL GET_INTERSECTION_ll (1,IJU_ll-KRIMY-JPHEXT+1,IIU_ll,IJU_ll,       & ! -KRIMY
                             IIBINT,IJBINT,IIEINT,IJEINT,"EXTE",IINFO_ll)
  IF ( IINFO_ll == 0 ) THEN
    PWORK(:,IJBINT:IJE+JPHEXT,:) = PLBY(:,IDIMLB-(IJE+JPHEXT-IJBINT):IDIMLB,:) ! +1
  END IF
END IF

END SUBROUTINE EXPAND_LBY
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RELAXATION
